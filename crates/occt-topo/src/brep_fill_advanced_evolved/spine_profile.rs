//! `GetSpineAndProfile` (`BRepFill_AdvancedEvolved.cxx:104-371`).
//!
//! Decides which side of the spine the profile must sweep on and flips
//! `myProfile` when needed (`cxx:240-252`, `:260-268`, `:352-366`). Three
//! control-flow details of the reference are reproduced verbatim:
//!
//! - the second corner tangent is evaluated on the **first** edge's adaptor
//!   (`anAC1.D1(aPar2, aP, aT2)`, `cxx:141`, `:229`, `:308`), `anAC2` only
//!   being used for `BRepAdaptor_Curve` construction;
//! - the solution scan tests `anExtr.Value()` — the shape-shape reference
//!   distance — rather than the per-solution distance (`cxx:186-195`);
//! - the second `SupportTypeShape2` read inside the `IsOnEdge/IsVertex` arm
//!   (`cxx:285`) is the raw one, not the `IsInFace` override applied to
//!   `anExtrType2` at `cxx:200-202`.

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_core::precision::{CONFUSION, INFINITE, SQUARE_CONFUSION};
use occt_geom::geom_api::project_point_on_surface;
use occt_geom::plane::GeomPlane;

use crate::abs::{Orientation, ShapeType};
use crate::brep_class::FaceClassifier;
use crate::brep_extrema_dist_shape_shape::DistShapeShape;
use crate::brep_extrema_solution_elem::SupportType;
use crate::brep_lib_make_face::make_face_from_wire;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::topopebrep::tool::adaptor::CurveAdaptor;

use super::maps::ordered_map_shapes_and_ancestors;
use super::AdvancedEvolved;

/// The two tangents `anAC1.D1(aPar1, ...)` / `anAC1.D1(aPar2, ...)` of
/// `cxx:136-141` (`:220-230`, `:299-308`): both taken from the adaptor of the
/// **first** ancestor edge, as the reference does.
fn corner_tangents(e1: &TopoShape, e2: &TopoShape, vertex: &TopoShape) -> Option<(GpVec, GpVec)> {
    let ac1 = CurveAdaptor::of_edge(&Edge(e1.clone()));
    let par1 = BRepTool::parameter_of_vertex(e1, vertex);
    let par2 = BRepTool::parameter_of_vertex(e2, vertex);
    let (_, t1) = ac1.d1(par1)?;
    let (_, t2) = ac1.d1(par2)?;
    Some((t1, t2))
}

/// `BRepLib_MakeFace(aPln, theProfile)` (`cxx:158`).
///
/// UNPORTED: `BRepBuilderAPI_MakeFace(const gp_Pln&, const TopoDS_Wire&)`
/// (`BRepBuilderAPI_MakeFace.cxx:236-245`) has no port. The face is built
/// directly on `pln` with the wire as its outer wire, which is what its
/// `Add(W)` arm ends up with once the wire is accepted; the reference's
/// `IsDone()` guard therefore has no failing case here.
fn make_face_on_plane(pln: &GpPln, wire: &Wire) -> Face {
    let builder = TopoBuilder::new();
    builder.make_face(Arc::new(GeomPlane::new(pln.clone())), &[wire.clone()])
}

impl AdvancedEvolved {
    /// `GetSpineAndProfile(theSpine, theProfile)` (`cxx:104-371`).
    pub(super) fn get_spine_and_profile(&mut self, spine: &Wire, profile: &Wire) {
        // `cxx:108-109`.
        self.spine = Some(spine.clone());
        self.profile = Some(profile.clone());

        // `cxx:111-112`: `TopExp::MapShapesAndAncestors(theProfile, VERTEX,
        // EDGE, aMVEP)`.
        let mvep =
            ordered_map_shapes_and_ancestors(&profile.0, ShapeType::Vertex, ShapeType::Edge);

        let mut n2 = GpVec::zero();
        let mut loc = GpPnt::new(0.0, 0.0, 0.0);

        // `cxx:117-149`: first profile corner with non-parallel tangents.
        for entry in &mvep {
            if entry.ancestors.len() < 2 {
                continue;
            }
            let e1 = entry.ancestors.first().expect("len >= 2");
            let e2 = entry.ancestors.last().expect("len >= 2");
            let Some((t1, t2)) = corner_tangents(e1, e2, &entry.key) else {
                continue;
            };
            n2 = t1.crossed(&t2);
            if n2.square_magnitude() > SQUARE_CONFUSION {
                loc = BRepTool::vertex_point(&Vertex(entry.key.clone()));
                break;
            }
        }

        // `cxx:151-168`.
        let mut extr = DistShapeShape::default();
        extr.load_s1(&spine.0);
        if n2.square_magnitude() > SQUARE_CONFUSION {
            let Ok(dir) = GpDir::from_vec(&n2) else {
                return;
            };
            let Ok(ax3) = GpAx3::from_point_dir(loc, dir) else {
                return;
            };
            let pln = GpPln::new(ax3);
            let face = make_face_on_plane(&pln, profile);
            extr.load_s2(&face.0);
        } else {
            extr.load_s2(&profile.0);
        }

        // `cxx:170-178`.
        if !extr.perform() {
            return;
        }
        let nb_sol = extr.nb_solution();
        if nb_sol < 1 {
            return;
        }

        // `cxx:180-195`: `anExtr.Value()` is the shape-shape reference
        // distance, so every solution compares equal and the last one wins.
        let mut dist_min = INFINITE;
        let mut idx_min = 0usize;
        for sol_id in 1..=nb_sol {
            let d = extr.value();
            if d > dist_min {
                continue;
            }
            dist_min = d;
            idx_min = sol_id;
        }

        // `cxx:197-202`.
        let mut extr_type2 = extr.support_type_shape2(idx_min);
        if dist_min < CONFUSION {
            extr_type2 = SupportType::IsInFace;
        }

        // `cxx:204-369`.
        match extr_type2 {
            SupportType::IsInFace => {
                // `cxx:206-254`.
                if extr.support_type_shape1(idx_min) == SupportType::IsVertex {
                    let Some(v) = extr.support_on_shape1(idx_min) else {
                        return;
                    };
                    let mves = ordered_map_shapes_and_ancestors(
                        &spine.0,
                        ShapeType::Vertex,
                        ShapeType::Edge,
                    );
                    let key = crate::topo_tools::shape_map_key(&v);
                    let Some(entry) = mves.iter().find(|e| {
                        crate::topo_tools::shape_map_key(&e.key) == key
                    }) else {
                        return;
                    };
                    let Some(e1) = entry.ancestors.first() else {
                        return;
                    };
                    let Some(e2) = entry.ancestors.last() else {
                        return;
                    };
                    let Some((t1, t2)) = corner_tangents(e1, e2, &v) else {
                        return;
                    };
                    // `cxx:232-234`: "find minimal sine".
                    let sq_t1 = t1.square_magnitude().max(1.0 / INFINITE);
                    let sq_t2 = t2.square_magnitude().max(1.0 / INFINITE);
                    let sq_sin1 = t1.cross_square_magnitude(&n2) / sq_t1;
                    let sq_sin2 = t2.cross_square_magnitude(&n2) / sq_t2;
                    // `cxx:238-253`.
                    if sq_sin1 < sq_sin2 {
                        if t1.dot(&n2) > 0.0 {
                            self.reverse_profile();
                        }
                    } else if t2.dot(&n2) > 0.0 {
                        self.reverse_profile();
                    }
                } else {
                    // `cxx:254-268`: the solution lies on a spine edge.
                    let Some(e) = extr.support_on_shape1(idx_min) else {
                        return;
                    };
                    let ac = CurveAdaptor::of_edge(&Edge(e.clone()));
                    let par = extr.par_on_edge_s1(idx_min);
                    let Some((_, t1)) = ac.d1(par) else {
                        return;
                    };
                    if t1.dot(&n2) > 0.0 {
                        self.reverse_profile();
                    }
                }
            }
            SupportType::IsOnEdge | SupportType::IsVertex => {
                // `cxx:271-278`: `BRepLib_MakeFace(theSpine, true)`.
                let Some(f_spine) = make_face_from_wire(spine, true) else {
                    return;
                };
                // `cxx:279-283`: `down_cast<Geom_Plane>(BRep_Tool::Surface)`;
                // OCCT dereferences the handle without a null check.
                let Some(surf) = BRepTool::face_surface(&f_spine) else {
                    return;
                };
                let Some(pln_spine) = surf.gp_pln() else {
                    return;
                };
                let n1 = GpVec::from_xyz(pln_spine.axis().direction().xyz());

                // `cxx:285-346`: the tangent of the profile at the solution.
                let tan_v: GpVec;
                if extr.support_type_shape2(idx_min) == SupportType::IsVertex {
                    // `cxx:287-336`.
                    let Some(v) = extr.support_on_shape2(idx_min) else {
                        return;
                    };
                    let mves = ordered_map_shapes_and_ancestors(
                        &profile.0,
                        ShapeType::Vertex,
                        ShapeType::Edge,
                    );
                    let key = crate::topo_tools::shape_map_key(&v);
                    let Some(entry) = mves.iter().find(|e| {
                        crate::topo_tools::shape_map_key(&e.key) == key
                    }) else {
                        return;
                    };
                    let Some(e1) = entry.ancestors.first() else {
                        return;
                    };
                    let Some(e2) = entry.ancestors.last() else {
                        return;
                    };
                    let Some((t1, t2)) = corner_tangents(e1, e2, &v) else {
                        return;
                    };
                    // `cxx:310-336`: "find maximal cosine".
                    let mut sq_t1 = t1.square_magnitude();
                    let mut sq_t2 = t2.square_magnitude();
                    if sq_t1 < SQUARE_CONFUSION {
                        sq_t1 = INFINITE;
                    }
                    if sq_t2 < SQUARE_CONFUSION {
                        sq_t2 = INFINITE;
                    }
                    let dp1 = t1.dot(&n1);
                    let dp2 = t2.dot(&n1);
                    tan_v = if dp1 * dp1 * sq_t2 > dp2 * dp2 * sq_t1 {
                        t1
                    } else {
                        t2
                    };
                } else {
                    // `cxx:337-346`.
                    let Some(e) = extr.support_on_shape2(idx_min) else {
                        return;
                    };
                    let ac = CurveAdaptor::of_edge(&Edge(e.clone()));
                    let par = extr.par_on_edge_s2(idx_min);
                    let Some((_, t)) = ac.d1(par) else {
                        return;
                    };
                    tan_v = t;
                }

                // `cxx:348-366`.
                let p_near = extr.point_on_shape2(idx_min);
                let state = classify_on_face(&f_spine, &p_near);
                if state != FaceState::Out {
                    if n1.dot(&tan_v) < 0.0 {
                        self.reverse_profile();
                    }
                } else if n1.dot(&tan_v) > 0.0 {
                    self.reverse_profile();
                }
            }
        }
    }

    /// `myProfile.Reverse()` (`cxx:242`, `:248`, `:263`, `:355`, `:363`).
    fn reverse_profile(&mut self) {
        if let Some(p) = &self.profile {
            self.profile = Some(Wire(p.0.oriented(Orientation::Reversed)));
        }
    }
}

/// `BRepClass_FaceClassifier(aFSpine, aPnear, Precision::Confusion())`
/// (`cxx:352-353`): the 3D point is projected to the face's parameters and the
/// 2D classifier runs on them (`BRepClass_FaceClassifier.cxx:47-56`).
fn classify_on_face(face: &Face, p: &GpPnt) -> FaceState {
    let mut cls = FaceClassifier::new();
    let Some(surf) = BRepTool::face_surface(face) else {
        return FaceState::Unknown;
    };
    let Some(pos) = project_point_on_surface(surf.as_ref(), p, CONFUSION) else {
        return FaceState::Unknown;
    };
    cls.perform(face, GpPnt2d::new(pos.u, pos.v), CONFUSION);
    cls.state()
}

