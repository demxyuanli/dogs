//! `Draft_Modification::Propagate` (`Draft_Modification_1.cxx:424-705`): fills
//! the edge and vertex records implied by the faces `InternalAdd` marked, then
//! invents a face for the free borders the neutral-plane cut leaves behind.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | edges / vertices / parameters of the modified faces | `cxx:434-473` | `propagate` |
//! | edges containing a modified vertex | `cxx:476-512` | `propagate` |
//! | faces containing a modified edge | `cxx:515-551` | `propagate` |
//! | faces for the free borders | `cxx:553-702` | `propagate` |
//!
//! # Translation notes
//!
//! * `BRep_Tool::Surface(F, L)` followed by `Transformed(L.Transformation())` is
//!   [`BRepTool::face_surface_world`]; the cylinder / cone domain widening of
//!   `cxx:530-542` tests the *untransformed* surface (`S->DynamicType()`) and
//!   trims the *transformed* one (`NewS`), exactly as the reference does.
//! * `GeomAdaptor_SurfaceOfLinearExtrusion` (`cxx:658-676`) is
//!   [`GeomAdaptorSurfaceOfLinearExtrusion`]; the `default` arm keeps the
//!   extrusion surface, not a plane.
//! * `Draft_VertexInfo`'s cursor loop (`cxx:607-633`) is
//!   [`DraftVertexInfo::init_edge_iterator`] / [`DraftVertexInfo::more_edge`] /
//!   [`DraftVertexInfo::edge`] / [`DraftVertexInfo::next_edge`].

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpVec};
use occt_core::precision::{Precision, CONFUSION};
use occt_geom::adaptor_surface_of_linear_extrusion::GeomAdaptorSurfaceOfLinearExtrusion;
use occt_geom::cylinder::GeomCylinder;
use occt_geom::plane::GeomPlane;
use occt_geom::rectangular_trimmed::GeomRectangularTrimmedSurface;
use occt_geom::surface_of_linear_extrusion::GeomSurfaceOfLinearExtrusion;
use occt_geom::Curve;
use occt_geom::Surface;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::brep_tools::uv_bounds;
use crate::draft::edge_info::DraftEdgeInfo;
use crate::draft::error_status::DraftErrorStatus;
use crate::draft::face_info::DraftFaceInfo;
use crate::draft::vertex_info::DraftVertexInfo;
use crate::shape::{Edge, Face, Vertex};
use crate::topo_tools::{explorer_shapes, is_same};

use super::{make_face, DraftModification};

impl DraftModification {
    /// `Propagate()` (`Draft_Modification_1.cxx:424-705`).
    ///
    /// Called by `Perform` (`_1.cxx:718`, see [`super::perform`]).
    pub(super) fn propagate(&mut self) -> bool {
        // `cxx:426-431`
        if self.bad_shape.is_some() {
            return false;
        }

        // `cxx:434-473`: every edge and vertex of the faces already recorded.
        for i in 1..=self.my_f_map.extent() {
            let fc = Face(self.my_f_map.find_key(i).clone());
            for e in explorer_shapes(&fc.0, ShapeType::Edge) {
                let e = Edge(e);
                // `cxx:449-453`
                if !self.my_e_map.contains(&e.0) {
                    self.my_e_map.add(&e.0, DraftEdgeInfo::with_new_geometry(true));
                }
                // `cxx:454`
                if let Some(ei) = self.my_e_map.change_from_key(&e.0) {
                    ei.add(&fc);
                }
                // `cxx:456-471`
                for v in explorer_shapes(&e.0, ShapeType::Vertex) {
                    let v = Vertex(v);
                    if !self.my_v_map.contains(&v.0) {
                        self.my_v_map.add(&v.0, DraftVertexInfo::new());
                    }
                    let param = BRepTool::parameter_of_vertex(&e.0, &v.0);
                    if let Some(vi) = self.my_v_map.change_from_key(&v.0) {
                        vi.add(&e);
                        *vi.change_parameter(&e) = param;
                    }
                }
            }
        }

        // `cxx:476-512`: every edge of the base shape that carries a modified
        // vertex joins the edge map with the `false` (no new geometry) seed.
        let base_edges: Vec<Edge> = explorer_shapes(&self.my_shape, ShapeType::Edge)
            .into_iter()
            .map(Edge)
            .collect();
        for i in 1..=self.my_v_map.extent() {
            let vt = Vertex(self.my_v_map.find_key(i).clone());
            for e in &base_edges {
                // `cxx:490-501`
                let mut found = false;
                for v in explorer_shapes(&e.0, ShapeType::Vertex) {
                    if is_same(&vt.0, &v) {
                        found = true;
                        break;
                    }
                }
                if !found {
                    continue;
                }
                // `cxx:502-510`
                if !self.my_e_map.contains(&e.0) {
                    self.my_e_map.add(&e.0, DraftEdgeInfo::with_new_geometry(false));
                }
                let param = BRepTool::parameter_of_vertex(&e.0, &vt.0);
                if let Some(vi) = self.my_v_map.change_from_key(&vt.0) {
                    vi.add(e);
                    *vi.change_parameter(e) = param;
                }
            }
        }

        // `cxx:515-551`: the faces of every modified edge join the face map
        // without a new geometry, keeping only the widened domain.
        for i in 1..=self.my_e_map.extent() {
            let ed = Edge(self.my_e_map.find_key(i).clone());
            let ancestors: Vec<Face> = self
                .my_ef_map
                .find_from_key(&ed.0)
                // `cxx:520`: `myEFMap.FindFromKey(Ed)` throws for an edge that
                // is not in the map (`Add` after `Perform` empties it, `cxx:71`).
                .expect("Standard_NoSuchObject: NCollection_IndexedDataMap::FindFromKey")
                .iter()
                .map(|s| Face(s.clone()))
                .collect();
            for f in ancestors {
                // `cxx:521-548`
                if !self.my_f_map.contains(&f.0) {
                    if let Some(new_s) = widened_basis_surface(&f) {
                        // `cxx:546-547`: `Draft_FaceInfo FInf(NewS, false)` -
                        // no root face is recorded here.
                        let fi = DraftFaceInfo::with_surface(Some(&new_s), false);
                        self.my_f_map.add(&f.0, fi);
                    }
                }
                // `cxx:549`
                if let Some(ei) = self.my_e_map.change_from_key(&ed.0) {
                    ei.add(&f);
                }
            }
        }

        // `cxx:553-702`: invent a face for the free borders.
        for i in 1..=self.my_e_map.extent() {
            let ek = Edge(self.my_e_map.find_key(i).clone());
            let (new_geom, geom_is_null, second_is_null, first_face) = {
                let einf = self.my_e_map.find_from_index(i);
                (
                    einf.new_geometry(),
                    einf.geometry().is_none(),
                    einf.second_face().is_none(),
                    einf.first_face().cloned(),
                )
            };
            // `cxx:557`
            if !(new_geom && geom_is_null && second_is_null) {
                continue;
            }
            let Some(first_face) = first_face else {
                continue;
            };
            // `cxx:560-562`
            let Some(s1) = BRepTool::face_surface_world(&first_face) else {
                continue;
            };
            // `cxx:566-571`
            let Some(mut c) = BRepTool::edge_curve_world(&ek) else {
                continue;
            };
            if c.is_geom_trimmed() {
                if let Some(b) = c.trim_basis() {
                    c = b;
                }
            }

            // `cxx:572`: `Handle(Geom_Surface) S2;`
            let s2: Option<Arc<dyn Surface>>;

            if s1.gp_pln().is_none() {
                // `cxx:573-601`
                if let Some(ax2) = conic_position(&c) {
                    // `cxx:578-580`: `gp_Ax3 thePl(down_cast<Geom_Conic>(C)->Position())`.
                    s2 = Some(Arc::new(GeomPlane::new(occt_core::gp::GpPln::new(
                        GpAx3::from_ax2(&ax2),
                    ))));
                } else if c.gp_line().is_some() {
                    // `cxx:580-594`
                    let axis = match s1.rectangular_trimmed_basis() {
                        Some(basis) => basis.gp_cylinder().map(|c| c.axis()),
                        None => s1.gp_cylinder().map(|c| c.axis()),
                    };
                    let Some(axis) = axis else {
                        self.bad_shape = Some(ek.0.clone());
                        self.err_stat = DraftErrorStatus::EdgeRecomputation;
                        break;
                    };
                    let they = GpVec::from_pnts(axis.location(), &c.d0(0.0));
                    let crossed = GpVec::from_xyz(axis.direction().xyz()).crossed(&they);
                    let Ok(axz) = GpDir::from_vec(&crossed) else {
                        self.bad_shape = Some(ek.0.clone());
                        self.err_stat = DraftErrorStatus::EdgeRecomputation;
                        break;
                    };
                    let Ok(ax3) = GpAx3::new(*axis.location(), axz, &axis.direction()) else {
                        self.bad_shape = Some(ek.0.clone());
                        self.err_stat = DraftErrorStatus::EdgeRecomputation;
                        break;
                    };
                    s2 = Some(Arc::new(GeomPlane::new(occt_core::gp::GpPln::new(ax3))));
                } else {
                    // `cxx:596-600`
                    self.bad_shape = Some(ek.0.clone());
                    self.err_stat = DraftErrorStatus::EdgeRecomputation;
                    break;
                }
            } else {
                // `cxx:603-687`: on a plane - extrude the neighbouring edge.
                let mut s2_built: Option<Arc<dyn Surface>> = None;
                for j in 1..=self.my_v_map.extent() {
                    let mut cursor_end = false;
                    let mut other_edge: Option<Edge> = None;
                    {
                        let vinf = self.my_v_map.change_from_index(j);
                        // `cxx:608-613`
                        vinf.init_edge_iterator();
                        while vinf.more_edge() {
                            if is_same(&vinf.edge().0, &ek.0) {
                                cursor_end = true;
                                break;
                            }
                            vinf.next_edge();
                        }
                    }
                    if !cursor_end {
                        continue;
                    }
                    {
                        let vinf = self.my_v_map.change_from_index(j);
                        // `cxx:617-629`: the first edge of the same vertex lying
                        // on another (or on none) of the two faces.
                        vinf.init_edge_iterator();
                        while vinf.more_edge() {
                            let edg = vinf.edge();
                            if !is_same(&edg.0, &ek.0) {
                                if let Some(ei) = self.my_e_map.find_from_key(&edg.0) {
                                    let not_first = ei
                                        .first_face()
                                        .map(|ff| !is_same(&ff.0, &first_face.0))
                                        .unwrap_or(true);
                                    let second_ok = match ei.second_face() {
                                        None => true,
                                        Some(sf) => !is_same(&sf.0, &first_face.0),
                                    };
                                    if not_first && second_ok {
                                        other_edge = Some(edg);
                                        break;
                                    }
                                }
                            }
                            vinf.next_edge();
                        }
                        cursor_end = vinf.more_edge();
                    }
                    if !cursor_end {
                        continue;
                    }
                    // `cxx:632-677`
                    let Some(other_edge) = other_edge else {
                        self.bad_shape = Some(ek.0.clone());
                        self.err_stat = DraftErrorStatus::EdgeRecomputation;
                        break;
                    };
                    let Some(mut c2) = BRepTool::edge_curve_world(&other_edge) else {
                        self.bad_shape = Some(ek.0.clone());
                        self.err_stat = DraftErrorStatus::EdgeRecomputation;
                        break;
                    };
                    if c2.is_geom_trimmed() {
                        if let Some(b) = c2.trim_basis() {
                            c2 = b;
                        }
                    }
                    let (hcur, direc) = if let Some(line) = c.gp_line() {
                        (c2.clone(), line.direction())
                    } else if let Some(line2) = c2.gp_line() {
                        (c.clone(), line2.direction())
                    } else {
                        self.bad_shape = Some(ek.0.clone());
                        self.err_stat = DraftErrorStatus::EdgeRecomputation;
                        break;
                    };
                    // `cxx:658-677`
                    let sle = GeomAdaptorSurfaceOfLinearExtrusion::new(hcur.clone(), direc);
                    match sle.get_type() {
                        occt_core::kernel::geomabs::SurfaceType::Plane => {
                            s2_built = sle
                                .plane()
                                .map(|pl| Arc::new(GeomPlane::new(pl)) as Arc<dyn Surface>);
                        }
                        occt_core::kernel::geomabs::SurfaceType::Cylinder => {
                            s2_built = sle
                                .cylinder()
                                .map(|cy| Arc::new(GeomCylinder::new(cy)) as Arc<dyn Surface>);
                        }
                        _ => {
                            s2_built = Some(Arc::new(GeomSurfaceOfLinearExtrusion::new(
                                hcur.clone(),
                                direc,
                            )) as Arc<dyn Surface>);
                        }
                    }
                    break;
                }
                if self.bad_shape.is_some() {
                    break;
                }
                s2 = s2_built;
            }

            // `cxx:688-699`
            if self.bad_shape.is_some() {
                break;
            }
            let Some(s2) = s2 else {
                self.bad_shape = Some(ek.0.clone());
                self.err_stat = DraftErrorStatus::EdgeRecomputation;
                break;
            };
            let new_face = make_face(s2.clone(), CONFUSION);
            let ei = self.my_e_map.change_from_index(i);
            ei.add(&new_face);
            let fi = DraftFaceInfo::with_surface(Some(&s2), false);
            self.my_f_map.add(&new_face.0, fi);
        }

        // `cxx:704`
        self.bad_shape.is_none()
    }
}

/// `down_cast<Geom_Conic>(C)->Position()` (`_1.cxx:578-579`): the `gp_Ax2` of a
/// circle, ellipse, parabola or hyperbola, `None` for every other curve.
fn conic_position(c: &Arc<dyn Curve>) -> Option<occt_core::gp::GpAx2> {
    if let Some(v) = c.gp_circ() {
        return Some(v.position());
    }
    if let Some(v) = c.gp_ellipse() {
        return Some(*v.position());
    }
    if let Some(v) = c.gp_parabola() {
        return Some(*v.position());
    }
    if let Some(v) = c.gp_hyperbola() {
        return Some(*v.position());
    }
    None
}

/// `cxx:526-547`: the surface of `F` with its cylinder / cone `v` domain
/// widened by ten spans, so the cut by the neutral plane stays inside the
/// domain. Answers `None` where `BRep_Tool::Surface` would have thrown.
fn widened_basis_surface(f: &Face) -> Option<Arc<dyn Surface>> {
    let s = BRepTool::face_surface(f)?;
    let new_s: Arc<dyn Surface> = Arc::from(
        s.transformed(&f.0.location().transformation()),
    );
    // `cxx:530-543`
    if s.gp_cylinder().is_some() || s.gp_cone().is_some() {
        let (_, _, vmin, vmax) = uv_bounds(f);
        if !Precision::is_negative_infinite(vmin) && !Precision::is_positive_infinite(vmax) {
            let deltav = 10.0 * (vmax - vmin);
            let trimmed = GeomRectangularTrimmedSurface::uv(
                new_s,
                0.0,
                2.0 * std::f64::consts::PI,
                vmin - deltav,
                vmax + deltav,
            );
            if let Ok(trimmed) = trimmed {
                return Some(Arc::new(trimmed));
            }
            return None;
        }
    }
    Some(new_s)
}
