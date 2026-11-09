//! `BRepOffset_Offset::Init(Path, Edge1, Edge2, [FirstEdge, LastEdge,]
//! Offset, Polynomial, Tol, Conti)` (`BRepOffset_Offset.cxx:1098-1553`): the
//! tube swept along `Path` between the two offset edges `Edge1` / `Edge2`.
//!
//! Control flow, all from the `.cxx`:
//!
//! * `:1122-1192` take the three curves, baking each edge's location, and
//!   detect a degenerated (`Radius() < Precision::Confusion()`) circular edge
//!   end (`C1Denerated` / `C2Denerated`). An edge without a 3D curve is
//!   described by its `CurveOnSurface` adaptor instead (`:1139-1150`,
//!   `:1169-1180`), which is the usual case here because the offset faces
//!   built by [`super::BRepOffsetOffset::init_face`] keep p-curves.
//! * `:1192-1200` `GeomFill_Pipe(HCP, HEdge1, HEdge2, |Offset|)` +
//!   `Perform(Tol, Polynomial, Conti)`.
//! * `:1201-1320` the pipe surface, its `ExchangeUV`, and the two straight
//!   p-curves `Edge1` / `Edge2` carry on it, with `SameRange` / `Range` /
//!   `SameParameter`.
//! * `:1322-1528` the two end edges `Edge3` / `Edge4`: re-use `FirstEdge` /
//!   `LastEdge` when the caller already built one, otherwise create an edge
//!   out of the two endpoint vertices, give it the iso p-curve of the pipe
//!   surface and, when it is not degenerated, a 3D curve through the
//!   file-local `ComputeCurve3d` (`compute_curve3d.rs`).
//! * `:1530-1552` the wire `[Edge1(reversed), Edge2, Edge4(reversed), Edge3]`
//!   (reversed under `ExchangeUV`), then the face orientation.
//!
//! UNPORTED, with the OCCT lines that keep it out:
//!
//! * the `Polynomial` approximation arm: the flag is forwarded to
//!   `GeomFill_Pipe::Perform` like `:1194`, but `GeomConvert_ApproxSurface`
//!   (used by the vertex overload, `:1623-1632`) has no counterpart, and every
//!   `BRepOffset_MakeOffset` call site uses the default `Polynomial = false`.
//! * `OCCT_DEBUG` blocks (`:1349-1358`, `:1431-1439`): debug output only.
//! * `BRepTools::Update(myFace)` (`:1548`): the port keeps no cached UV points
//!   nor triangulation.

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{APPROXIMATION, CONFUSION};
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::geom_fill::pipe::Pipe;
use occt_geom::{Curve, GeomTrimmedCurveBasis, Surface};
use occt_geom2d::{Curve2d, Geom2dLine};

use crate::abs::{Orientation, ShapeType};
use crate::brep_lib_build_curves3d::build_curve3d;
use crate::brep_lib_same_parameter::same_parameter_edge_inplace;
use crate::brep_lib_same_range::same_range;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, is_same};

use super::compute_curve3d::{
    compute_curve3d, update_edge_curve3d, update_edge_pcurve, update_edge_pcurves,
};
use super::{BRepOffsetOffset, BRepOffsetStatus};

/// `gp::DX2d()` (`gp.hxx`).
const DX2D: GpDir2d = GpDir2d { x: 1.0, y: 0.0 };
/// `gp::DY2d()` (`gp.hxx`).
const DY2D: GpDir2d = GpDir2d { x: 0.0, y: 1.0 };

/// `BRep_Tool::CurveOnSurface(E, C, S, L, First, Last)`: the first curve-on-
/// surface representation of `e`, with the surface returned for the caller to
/// bake the edge location into (`cxx:1142-1148`, `:1171-1177`).
fn curve_on_surface(e: &Edge) -> Option<(Arc<dyn Curve2d>, Arc<dyn Surface>, f64, f64)> {
    GeometryRegistry::global()
        .edge_pcurve_reps(&e.0)
        .into_iter()
        .next()
        .map(|(surf, c, f, l)| (c, surf, f, l))
}

/// `Geom_TrimmedCurve(C, f, l)` + `Transform(Loc)` on the edge's own curve:
/// the `[f, l]` view of the edge curve in world coordinates.
pub(super) fn trimmed_world(e: &Edge, c: &Arc<dyn Curve>, f: f64, l: f64) -> Arc<dyn Curve> {
    let trimmed: Arc<dyn Curve> = Arc::new(GeomTrimmedCurveBasis::new(c.clone(), f, l));
    if e.0.location().is_identity() {
        trimmed
    } else {
        Arc::from(trimmed.transformed(&e.0.location().transformation()))
    }
}

/// `Adaptor3d_CurveOnSurface(Geom2dAdaptor_Curve(C2d), GeomAdaptor_Surface(S))`
/// for a p-curve-only edge (`cxx:1142-1150`, `:1171-1180`), with the edge
/// location baked into the surface.
fn pcurve_adaptor(e: &Edge, c2d: Arc<dyn Curve2d>, surf: Arc<dyn Surface>) -> Arc<dyn Curve> {
    let s = if e.0.location().is_identity() {
        surf
    } else {
        Arc::from(surf.transformed(&e.0.location().transformation()))
    };
    Arc::new(CurveOnSurface::new(c2d, s))
}

/// The pipe p-curve of one end edge: a `Geom2d_Line` over the pipe surface
/// (`cxx:1219`, `:1228`, `:1274`, `:1283`, `:1385-1392`, `:1453-1457`,
/// `:1485-1489`).
fn end_pcurve(p: GpPnt2d, d: GpDir2d) -> Geom2dLine {
    Geom2dLine::from_pnt_dir(p, d)
}

impl BRepOffsetOffset {
    /// `BRepOffset_Offset::Init(Path, Edge1, Edge2, Offset, Polynomial, Tol,
    /// Conti)` (`cxx:1098-1107`): the `FirstEdge` / `LastEdge` overload with
    /// two null edges.
    pub fn init_pipe(
        &mut self,
        path: &Edge,
        edge1: &mut Edge,
        edge2: &mut Edge,
        offset: f64,
        polynomial: bool,
        tol: f64,
        conti: Shape,
    ) -> Result<(), String> {
        self.init_pipe_first_last(path, edge1, edge2, offset, None, None, polynomial, tol, conti)
    }

    /// `BRepOffset_Offset::Init(Path, Edge1, Edge2, Offset, FirstEdge,
    /// LastEdge, Polynomial, Tol, Conti)` (`cxx:1112-1553`).
    ///
    /// `path` is the source edge the tube runs along; `edge1` / `edge2` are its
    /// images on the two offset faces, which this call rewrites in place
    /// (`UpdateEdge` + `Range` + `SameRange` / `SameParameter`).
    #[allow(clippy::too_many_arguments)]
    pub fn init_pipe_first_last(
        &mut self,
        path: &Edge,
        edge1: &mut Edge,
        edge2: &mut Edge,
        offset: f64,
        first_edge: Option<&Edge>,
        last_edge: Option<&Edge>,
        polynomial: bool,
        tol: f64,
        conti: Shape,
    ) -> Result<(), String> {
        let mut c1_degenerated = false;
        let mut c2_degenerated = false;
        self.my_status = BRepOffsetStatus::Good;
        self.my_shape = Some(path.0.clone());

        // `cxx:1130-1135`: `CP = Geom_TrimmedCurve(Curve(Path), f, l)`, location
        // baked; `HCP = GeomAdaptor_Curve(CP)`.
        let (pf, pl) = BRepTool::edge_parameters(path);
        let cp = BRepTool::edge_curve(path)
            .ok_or("BRepOffset_Offset::Init(Path): the path has no 3D curve")?;
        let hcp = trimmed_world(path, &cp, pf, pl);

        // `cxx:1136-1161`: `Edge1`'s curve, p-curve adaptor when it has no 3D
        // curve, `C1Denerated` on a null-radius circle.
        let (e1f, e1l) = BRepTool::edge_parameters(edge1);
        let c1 = BRepTool::edge_curve(edge1);
        let c1_is3d = c1.is_some();
        let hedge1: Arc<dyn Curve> = match c1.as_ref() {
            None => {
                let (c2d, surf, ..) = curve_on_surface(edge1)
                    .ok_or("BRepOffset_Offset::Init(Path): Edge1 has no curve")?;
                pcurve_adaptor(edge1, c2d, surf)
            }
            Some(c1) => {
                if let Some(r) = c1.circle_radius() {
                    c1_degenerated = r < CONFUSION;
                }
                trimmed_world(edge1, c1, e1f, e1l)
            }
        };

        // `cxx:1164-1191`: the same for `Edge2`.
        let (e2f, e2l) = BRepTool::edge_parameters(edge2);
        let c2 = BRepTool::edge_curve(edge2);
        let c2_is3d = c2.is_some();
        let hedge2: Arc<dyn Curve> = match c2.as_ref() {
            None => {
                let (c2d, surf, ..) = curve_on_surface(edge2)
                    .ok_or("BRepOffset_Offset::Init(Path): Edge2 has no curve")?;
                pcurve_adaptor(edge2, c2d, surf)
            }
            Some(c2) => {
                if let Some(r) = c2.circle_radius() {
                    c2_degenerated = r < CONFUSION;
                }
                trimmed_world(edge2, c2, e2f, e2l)
            }
        };

        // `cxx:1193-1199`: `GeomFill_Pipe Pipe(HCP, HEdge1, HEdge2, |Offset|)`
        // then `Perform(Tol, Polynomial, Conti)` with the header defaults
        // `MaxDegree = 11`, `NbMaxSegment = 30` (`GeomFill_Pipe.hxx`).
        let mut pipe = Pipe::with_path_curve1_curve2_radius(hcp, hedge1, hedge2, offset.abs())
            .map_err(|e| format!("GeomFill_Pipe: {e}"))?;
        pipe.perform_with_tol(tol, polynomial, conti, 11, 30)
            .map_err(|e| format!("GeomFill_Pipe: {e}"))?;
        if !pipe.is_done() {
            return Err("GeomFill_Pipe : Cannot make a surface".into());
        }
        let error_pipe = pipe.error_on_surf();
        let s = pipe
            .surface()
            .ok_or("GeomFill_Pipe : Cannot make a surface")?;
        let exch_uv = pipe.exchange_uv();
        let (f1, l1) = s.u_range();
        let (f2, l2) = s.v_range();

        // `cxx:1207-1215`: `MakeFace` + `UpdateFace(myFace, S, Id, PathTol)`.
        let path_tol = BRepTool::edge_tolerance(path);
        let builder = TopoBuilder::new();
        let mut my_face = builder.make_face_with_tolerance(s.clone(), &[], path_tol);

        // `cxx:1217-1268`: `Edge1`'s p-curve on the pipe surface, the 3D curve
        // it gains when it had none, then its range.
        let (mut pc, mut u1, mut u2) = if exch_uv {
            (
                end_pcurve(GpPnt2d::new(0.0, f2), DX2D),
                f1,
                l1,
            )
        } else {
            (
                end_pcurve(GpPnt2d::new(f1, 0.0), DY2D),
                f2,
                l2,
            )
        };
        if !c1_is3d {
            let iso = if exch_uv {
                s.v_iso_curve(f2)
            } else {
                s.u_iso_curve(f1)
            }
            .ok_or("BRepOffset_Offset::Init(Path): Edge1 iso curve is null")?;
            update_edge_curve3d(edge1, &iso, BRepTool::edge_tolerance(edge1));
        } else if c1_degenerated {
            // `UpdateEdge(Edge1, Dummy, Id, Tol)` with a null curve, then
            // `myBuilder.Degenerated(Edge1, true)` (`cxx:1246-1248`).
            GeometryRegistry::global().remove_edge_curve3d(&edge1.0);
            GeometryRegistry::global().set_degenerated(&edge1.0, true);
        }

        let mut the_tol = path_tol.max(BRepTool::edge_tolerance(edge1) + error_pipe);
        let pc1_arc: Arc<dyn Curve2d> = Arc::new(pc.clone());
        update_edge_pcurve(edge1, &pc1_arc, &my_face, the_tol);

        // `cxx:1253-1260`.
        if !c1_is3d && !c1_degenerated {
            GeometryRegistry::global().set_same_range(&edge1.0, false);
            GeometryRegistry::global().set_edge_range(&edge1.0, u1, u2);
        }
        GeometryRegistry::global().set_pcurve_range(
            &edge1.0,
            GeometryRegistry::shape_key(&my_face.0),
            u1,
            u2,
        );
        same_range(&edge1.0, CONFUSION);

        // `cxx:1262-1268`: `SameParameter` only when the pipe is exact.
        if error_pipe == 0.0 {
            the_tol = the_tol.max(tol);
            GeometryRegistry::global().set_same_parameter(&edge1.0, false);
            same_parameter_edge_inplace(edge1, the_tol);
        }

        // `cxx:1270-1320`: the same for `Edge2`.
        let (pc2, u21, u22) = if exch_uv {
            (end_pcurve(GpPnt2d::new(0.0, l2), DX2D), f1, l1)
        } else {
            (end_pcurve(GpPnt2d::new(l1, 0.0), DY2D), f2, l2)
        };
        pc = pc2;
        u1 = u21;
        u2 = u22;
        if !c2_is3d {
            let iso = if exch_uv {
                s.v_iso_curve(l2)
            } else {
                s.u_iso_curve(l1)
            }
            .ok_or("BRepOffset_Offset::Init(Path): Edge2 iso curve is null")?;
            update_edge_curve3d(edge2, &iso, BRepTool::edge_tolerance(edge2));
        } else if c2_degenerated {
            GeometryRegistry::global().remove_edge_curve3d(&edge2.0);
            GeometryRegistry::global().set_degenerated(&edge2.0, true);
        }

        the_tol = path_tol.max(BRepTool::edge_tolerance(edge2) + error_pipe);
        let pc2_arc: Arc<dyn Curve2d> = Arc::new(pc);
        update_edge_pcurve(edge2, &pc2_arc, &my_face, the_tol);

        // `cxx:1305-1318`.
        GeometryRegistry::global().set_same_range(&edge2.0, false);
        if !c2_is3d && !c2_degenerated {
            GeometryRegistry::global().set_edge_range(&edge2.0, u1, u2);
        }
        GeometryRegistry::global().set_pcurve_range(
            &edge2.0,
            GeometryRegistry::shape_key(&my_face.0),
            u1,
            u2,
        );
        same_range(&edge2.0, CONFUSION);

        if error_pipe == 0.0 {
            the_tol = the_tol.max(tol);
            GeometryRegistry::global().set_same_parameter(&edge2.0, false);
            same_parameter_edge_inplace(edge2, the_tol);
        }

        // `cxx:1322-1333`: the end vertices.
        let (pv1, pv2) = edge_vertices(path);
        let is_closed = match (pv1.as_ref(), pv2.as_ref()) {
            (Some(a), Some(b)) => is_same(&a.0, &b.0),
            _ => false,
        };
        let (v1f, v1l) = edge_vertices(edge1);
        let (v2f, v2l) = edge_vertices(edge2);
        let start_degenerated = is_same_opt(v1f.as_ref(), v2f.as_ref());
        let end_degenerated = is_same_opt(v1l.as_ref(), v2l.as_ref());

        // `cxx:1336-1372`: `Edge3` (start edge).
        let mut e3_rev = false;
        let edge3 = match first_edge {
            Some(fe) => {
                let mut e = Edge(fe.0.oriented(Orientation::Forward));
                let (vvf, _) = edge_vertices(&e);
                let shares =
                    is_same_opt(vvf.as_ref(), v1f.as_ref()) || is_same_opt(vvf.as_ref(), v2f.as_ref());
                if !shares {
                    // `cxx:1360-1364`: the vertices are not shared, so build
                    // another edge and let the sewing fix it.
                    new_end_edge(&builder, v1f.as_ref(), v2f.as_ref())
                } else {
                    if !is_same_opt(vvf.as_ref(), v1f.as_ref()) {
                        e.0.reverse();
                        e3_rev = true;
                    }
                    e
                }
            }
            None => new_end_edge(&builder, v1f.as_ref(), v2f.as_ref()),
        };

        // `cxx:1372-1375`.
        let mut edge4 = if is_closed {
            Some(edge3.clone())
        } else {
            None
        };

        if is_closed {
            // `cxx:1379-1417`.
            let (mut l1_pc, mut l2_pc) = if exch_uv {
                (
                    end_pcurve(GpPnt2d::new(f1, 0.0), DY2D),
                    end_pcurve(GpPnt2d::new(l1, 0.0), DY2D),
                )
            } else {
                (
                    end_pcurve(GpPnt2d::new(0.0, f2), DX2D),
                    end_pcurve(GpPnt2d::new(0.0, l2), DX2D),
                )
            };
            if exch_uv {
                u1 = f2;
                u2 = l2;
            } else {
                u1 = f1;
                u2 = l1;
            }
            if e3_rev {
                l1_pc.reverse();
                l2_pc.reverse();
                let u = -u1;
                u1 = -u2;
                u2 = u;
            }
            let a1: Arc<dyn Curve2d> = Arc::new(l1_pc.clone());
            let a2: Arc<dyn Curve2d> = Arc::new(l2_pc);
            update_edge_pcurves(&edge3, &a1, &a2, &my_face, path_tol);
            GeometryRegistry::global().set_pcurve_range(
                &edge3.0,
                GeometryRegistry::shape_key(&my_face.0),
                u1,
                u2,
            );
            if start_degenerated {
                GeometryRegistry::global().set_degenerated(&edge3.0, true);
            } else if first_edge.is_none() {
                let l1_arc: Arc<dyn Curve2d> = Arc::new(l1_pc);
                compute_curve3d(&edge3, &l1_arc, &s, APPROXIMATION);
            }
        } else {
            // `cxx:1419-1451`: `Edge4` (end edge).
            let mut e4_rev = false;
            let e4 = match last_edge {
                Some(le) => {
                    let mut e = Edge(le.0.oriented(Orientation::Forward));
                    let (vvf, _) = edge_vertices(&e);
                    let shares = is_same_opt(vvf.as_ref(), v1l.as_ref())
                        || is_same_opt(vvf.as_ref(), v2l.as_ref());
                    if !shares {
                        new_end_edge(&builder, v1l.as_ref(), v2l.as_ref())
                    } else {
                        if !is_same_opt(vvf.as_ref(), v1l.as_ref()) {
                            e.0.reverse();
                            e4_rev = true;
                        }
                        e
                    }
                }
                None => new_end_edge(&builder, v1l.as_ref(), v2l.as_ref()),
            };

            // `cxx:1453-1482`: `Edge3`'s p-curve.
            let mut l1_pc = if exch_uv {
                u1 = f2;
                u2 = l2;
                end_pcurve(GpPnt2d::new(f1, 0.0), DY2D)
            } else {
                u1 = f1;
                u2 = l1;
                end_pcurve(GpPnt2d::new(0.0, f2), DX2D)
            };
            if e3_rev {
                l1_pc.reverse();
                let u = -u1;
                u1 = -u2;
                u2 = u;
            }
            let l1_arc: Arc<dyn Curve2d> = Arc::new(l1_pc.clone());
            update_edge_pcurve(&edge3, &l1_arc, &my_face, path_tol);
            GeometryRegistry::global().set_pcurve_range(
                &edge3.0,
                GeometryRegistry::shape_key(&my_face.0),
                u1,
                u2,
            );
            if start_degenerated {
                GeometryRegistry::global().set_degenerated(&edge3.0, true);
            } else if first_edge.is_none() {
                let l1_arc: Arc<dyn Curve2d> = Arc::new(l1_pc);
                compute_curve3d(&edge3, &l1_arc, &s, APPROXIMATION);
            }

            // `cxx:1484-1511`: `Edge4`'s p-curve.
            let mut l2_pc = if exch_uv {
                u1 = f2;
                u2 = l2;
                end_pcurve(GpPnt2d::new(l1, 0.0), DY2D)
            } else {
                u1 = f1;
                u2 = l1;
                end_pcurve(GpPnt2d::new(0.0, l2), DX2D)
            };
            if e4_rev {
                l2_pc.reverse();
                let u = -u1;
                u1 = -u2;
                u2 = u;
            }
            let l2_arc: Arc<dyn Curve2d> = Arc::new(l2_pc.clone());
            update_edge_pcurve(&e4, &l2_arc, &my_face, path_tol);
            GeometryRegistry::global().set_pcurve_range(
                &e4.0,
                GeometryRegistry::shape_key(&my_face.0),
                u1,
                u2,
            );
            if end_degenerated {
                GeometryRegistry::global().set_degenerated(&e4.0, true);
            } else if last_edge.is_none() {
                let l2_arc: Arc<dyn Curve2d> = Arc::new(l2_pc);
                compute_curve3d(&e4, &l2_arc, &s, APPROXIMATION);
            }
            edge4 = Some(e4);
        }

        // `cxx:1514-1528`: `BuildCurve3d` + `SameRange` + `SameParameter` when
        // the caller supplied the end edge.
        if first_edge.is_some() && !start_degenerated {
            let _ = build_curve3d(&edge3, path_tol, Shape::C1, 14, 0);
            GeometryRegistry::global().set_same_range(&edge3.0, false);
            GeometryRegistry::global().set_same_parameter(&edge3.0, false);
            same_parameter_edge_inplace(&edge3, tol);
        }
        if last_edge.is_some() && !end_degenerated {
            if let Some(e4) = edge4.as_ref() {
                let _ = build_curve3d(e4, path_tol, Shape::C1, 14, 0);
                GeometryRegistry::global().set_same_range(&e4.0, false);
                GeometryRegistry::global().set_same_parameter(&e4.0, false);
                same_parameter_edge_inplace(e4, tol);
            }
        }

        // `cxx:1530-1548`: the wire and the face.
        let mut w = Wire::new();
        builder.add_edge(&mut w, &Edge(edge1.0.oriented(Orientation::Reversed)));
        builder.add_edge(&mut w, &Edge(edge2.0.oriented(Orientation::Forward)));
        if let Some(e4) = edge4.as_ref() {
            let mut er = e4.clone();
            er.0.reverse();
            builder.add_edge(&mut w, &er);
        }
        builder.add_edge(&mut w, &edge3);
        if exch_uv {
            w.0.reverse();
        }

        builder.add_wire(&mut my_face, &w);
        if exch_uv {
            my_face.0.reverse();
        }

        // `cxx:1550-1552`: a reversed source edge reverses the face.
        if edge1.orientation() == Orientation::Reversed {
            my_face.0.reverse();
        }

        self.my_face = Some(my_face);
        Ok(())
    }
}

/// `v1.IsSame(v2)` tolerating a missing endpoint.
fn is_same_opt(a: Option<&crate::shape::Vertex>, b: Option<&crate::shape::Vertex>) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => is_same(&a.0, &b.0),
        _ => false,
    }
}

/// `myBuilder.MakeEdge(Edge)` + two `Add` calls (`cxx:1340-1342`,
/// `:1361-1363`, `:1420-1422`, `:1442-1444`).
fn new_end_edge(
    builder: &TopoBuilder,
    v1: Option<&crate::shape::Vertex>,
    v2: Option<&crate::shape::Vertex>,
) -> Edge {
    let mut e = Edge(TopoShape::new(ShapeType::Edge));
    if let (Some(a), Some(b)) = (v1, v2) {
        builder.add(&mut e.0, &a.0.oriented(Orientation::Forward));
        builder.add(&mut e.0, &b.0.oriented(Orientation::Reversed));
    }
    e
}

