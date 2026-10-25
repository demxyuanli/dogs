//! `TopOpeBRepDS_BuildTool` — the edge-update family.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepDS/TopOpeBRepDS_BuildTool.cxx`
//! (OCCT 8.0.0): `GetOrientedEdgeVertices` `:320-339`, `UpdateEdgeCurveTol`
//! `:343-468`, `UpdateEdge` `:1006-1028`, `Project` `:1031-1055`,
//! `Parameter(C, E, V)` `:1058-1065`, `PCurve(F, E, CDS, PC)` `:1243-1291`.

use std::sync::Arc;

use occt_core::elib::clib::adjust_periodic;
use occt_core::gp::{GpDir2d, GpTrsf2d, GpVec2d};
use occt_core::precision::{ANGULAR, PCONFUSION};
use occt_geom::curve::Curve;
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::explorer_shapes;

use super::super::curve::Curve as CurveData;
use super::BuildTool;

impl BuildTool {
    /// `GetOrientedEdgeVertices(E, Vmin, Vmax, Parmin, Parmax)` (`cxx:320-339`).
    ///
    /// `TopExp::Vertices(E, V1, V2)` walks the edge vertices in
    /// `TopoDS_Iterator` order, so the port takes the explorer result in
    /// order. The two parameters are written only when both vertices exist
    /// (`cxx:335-339`); OCCT leaves the caller's values untouched otherwise, so
    /// the port answers `None` for them.
    pub fn get_oriented_edge_vertices(
        &self,
        e: &Edge,
    ) -> Option<(TopoShape, TopoShape, Option<(f64, f64)>)> {
        let verts = explorer_shapes(&e.0, ShapeType::Vertex);
        let first = verts.first()?.clone();
        let last = verts.last()?.clone();
        let (vmin, vmax) = if e.0.orientation() == Orientation::Forward {
            (first, last)
        } else {
            (last, first)
        };
        let params = if !vmin.is_null() && !vmax.is_null() {
            Some((
                BRepTool::parameter_of_vertex(&e.0, &vmin),
                BRepTool::parameter_of_vertex(&e.0, &vmax),
            ))
        } else {
            None
        };
        Some((vmin, vmax, params))
    }

    /// `UpdateEdgeCurveTol(F1, F2, E, C3Dnew, tol3d, tol2d1, tol2d2, newtol,
    /// newparmin, newparmax)` (`cxx:343-468`).
    ///
    /// `F1`, `F2`, `tol3d`, `tol2d1` and `tol2d2` are unnamed in the definition
    /// (`cxx:344-356`) and never read; `newtol` is the constant `1.E-7`
    /// (`cxx:380`). The null-`C3Dnew` early return (`cxx:360-363`) is the
    /// caller's `None`. Returns `(newtol, newparmin, newparmax)`.
    pub fn update_edge_curve_tol(
        &self,
        e: &TopoShape,
        c3dnew: Option<Arc<dyn Curve>>,
    ) -> (f64, f64, f64) {
        let newtol = 1.0e-7;
        let Some(c3dnew) = c3dnew else {
            return (newtol, 0.0, 0.0);
        };
        let edge = Edge::wrap(e.clone()).expect("BuildTool::UpdateEdgeCurveTol: not an edge");

        // cxx:390-404
        let mut parmin = 0.0;
        let mut parmax = 0.0;
        let mut tolmin = newtol;
        let mut tolmax = newtol;
        let (vmin, vmax) = match self.get_oriented_edge_vertices(&edge) {
            Some((vmin, vmax, params)) => {
                if let Some((pmin, pmax)) = params {
                    parmin = pmin;
                    parmax = pmax;
                }
                (vmin, vmax)
            }
            None => (TopoShape::new(ShapeType::Vertex), TopoShape::new(ShapeType::Vertex)),
        };
        if let Some(v) = Vertex::wrap(vmin.clone()) {
            let t = BRepTool::vertex_tolerance(&v);
            if t > tolmin {
                tolmin = t;
            }
        }
        if let Some(v) = Vertex::wrap(vmax.clone()) {
            let t = BRepTool::vertex_tolerance(&v);
            if t > tolmax {
                tolmax = t;
            }
        }

        // cxx:409-428: a trimmed curve always takes its own range; so does a
        // B-spline; every other curve keeps the vertex parameters.
        let (newparmin, newparmax) =
            if c3dnew.is_geom_trimmed() || c3dnew.bspline_knots().is_some() {
                (c3dnew.first_parameter(), c3dnew.last_parameter())
            } else {
                (parmin, parmax)
            };

        // cxx:430-438
        let (p_min, p_max) = if vmin.orientation() == Orientation::Forward {
            (newparmin, newparmax)
        } else {
            (newparmax, newparmin)
        };
        if let (Some(vmin), Some(vmax)) = (Vertex::wrap(vmin.clone()), Vertex::wrap(vmax.clone())) {
            // `BRep_Builder::UpdateVertex(V, Par, E, Tol)` raises the vertex
            // tolerance as well (`BRep_TVertex::UpdateTolerance`,
            // `BRep_TVertex.cxx:29-35`).
            self.update_vertex(&vmin, p_min, &edge);
            GeometryRegistry::global().raise_vertex_tolerance(&vmin.0, tolmin);
            self.update_vertex(&vmax, p_max, &edge);
            GeometryRegistry::global().raise_vertex_tolerance(&vmax.0, tolmax);
        }

        // cxx:442
        self.curve3d(e, c3dnew.clone(), newtol);

        // cxx:446-466: every INTERNAL vertex is re-projected on the new curve.
        for vi in explorer_shapes(e, ShapeType::Vertex) {
            if vi.orientation() != Orientation::Internal {
                continue;
            }
            let Some(v) = Vertex::wrap(vi.clone()) else {
                continue;
            };
            let p = BRepTool::vertex_point(&v);
            let tolvi = BRepTool::vertex_tolerance(&v);
            let sols = occt_geom::extrema_pc::extrema_ext_pc_range(
                c3dnew.as_ref(),
                &p,
                newparmin,
                newparmax,
            );
            let mut best: Option<f64> = None;
            let mut best_d = f64::MAX;
            for s in &sols {
                if s.sq_dist < best_d {
                    best_d = s.sq_dist;
                    best = Some(s.u);
                }
            }
            if let Some(newpar) = best {
                self.update_vertex(&v, newpar, &edge);
                GeometryRegistry::global().raise_vertex_tolerance(&v.0, tolvi);
            }
        }

        (newtol, newparmin, newparmax)
    }

    /// `UpdateEdge(Ein, Eou)` (`cxx:1006-1028`).
    ///
    /// When `Eou`'s 3D curve is periodic and its own range is empty
    /// (`l2 <= f2`), the range is re-opened inside `Ein`'s range with
    /// `ElCLib::AdjustPeriodic(f1, l1, Precision::PConfusion(), f2n, l2n)`.
    /// A null curve on either side returns early (`cxx:1015-1018`).
    pub fn update_edge(&self, ein: &TopoShape, eou: &TopoShape) {
        let (Some(edge_in), Some(edge_out)) = (Edge::wrap(ein.clone()), Edge::wrap(eou.clone()))
        else {
            return;
        };
        if BRepTool::edge_curve(&edge_in).is_none() || BRepTool::edge_curve(&edge_out).is_none() {
            return;
        }
        let cout = BRepTool::edge_curve(&edge_out).unwrap();
        if cout.is_periodic() {
            let (f1, l1) = BRepTool::edge_parameters(&edge_in);
            let (f2, l2) = BRepTool::edge_parameters(&edge_out);
            let (mut f2n, mut l2n) = (f2, l2);
            if l2n <= f2n {
                adjust_periodic(f1, l1, PCONFUSION, &mut f2n, &mut l2n);
                self.range(eou, f2n, l2n);
            }
        }
    }

    /// `Project(C, V, p)` (`cxx:1031-1055`): the first minimum of
    /// `Extrema_ExtPC(P, GeomAdaptor_Curve(C))` whose point is within the
    /// vertex tolerance.
    pub fn project(c: &dyn Curve, v: &Vertex) -> Option<f64> {
        let p = BRepTool::vertex_point(v);
        let tol = BRepTool::vertex_tolerance(v);
        let (sols, done) = occt_geom::extrema_pc::extrema_ext_pc_range_done(
            c,
            &p,
            c.first_parameter(),
            c.last_parameter(),
        );
        if !done {
            return None;
        }
        for s in &sols {
            if s.is_min && p.distance(&s.point) <= tol {
                return Some(s.u);
            }
        }
        None
    }

    /// `Parameter(C, E, V)` (`cxx:1058-1065`).
    ///
    /// OCCT ignores `Project`'s result and hands `newparam` to `Parameter`
    /// either way — `newparam` is uninitialized when the projection fails, so
    /// the port skips the call in that case instead of writing garbage.
    pub fn parameter_with_curve(&self, c: &CurveData, e: &TopoShape, v: &TopoShape) {
        let Some(gc) = c.curve() else {
            return;
        };
        let Some(vert) = Vertex::wrap(v.clone()) else {
            return;
        };
        let Some(newparam) = Self::project(gc.as_ref(), &vert) else {
            return;
        };
        self.parameter(e, v, newparam);
    }

    /// `PCurve(F, E, CDS, PC)` (`cxx:1243-1291`).
    ///
    /// Before the pcurve is installed, a 2D-line pcurve whose parameter
    /// against the edge's 3D curve has drifted by more than
    /// `Precision::PConfusion()` is translated along the periodic direction of
    /// the face surface by `(CDSmin - Cf) * dir2d` (`cxx:1264-1284`).
    ///
    /// `CDSmin` is only read inside the `rangedef` guard in OCCT, but `deca`
    /// itself reads it unguarded (`cxx:1264`); the port evaluates `deca` only
    /// when the range is defined, which is the same branch `tran` selects.
    pub fn p_curve_with_ds(
        &self,
        f: &TopoShape,
        e: &TopoShape,
        cds: &CurveData,
        pc: Option<Arc<dyn Curve2d>>,
    ) {
        let Some(pc) = pc else {
            return;
        };
        let ff = Face::wrap(f.clone()).expect("BuildTool::PCurve: not a face");
        let ee = Edge::wrap(e.clone()).expect("BuildTool::PCurve: not an edge");

        let mut pct = pc;
        if let (Some((cdsmin, _)), Some(c)) = (cds.range(), BRepTool::edge_curve(&ee)) {
            let (cf, _) = BRepTool::edge_parameters(&ee);
            let deca = (cf - cdsmin).abs() > PCONFUSION;
            let line2d = pct.gp_lin2d();
            if deca && c.is_periodic() {
                if let (Some(line2d), Some(surf)) =
                    (line2d, GeometryRegistry::global().face_surface(&ff.0))
                {
                    let dir2d = *line2d.direction();
                    let dx2d = GpDir2d::new(1.0, 0.0).unwrap();
                    let dy2d = GpDir2d::new(0.0, 1.0).unwrap();
                    if surf.is_u_periodic() && dir2d.is_parallel(&dx2d, ANGULAR) {
                        pct = translate_2d(&pct, (cdsmin - cf) * dir2d.x(), 0.0);
                    } else if surf.is_v_periodic() && dir2d.is_parallel(&dy2d, ANGULAR) {
                        pct = translate_2d(&pct, 0.0, (cdsmin - cf) * dir2d.y());
                    }
                }
            }
        }

        // `TopOpeBRepDS_SetThePCurve(myBuilder, EE, FF, E.Orientation(), PCT)`
        // (`cxx:1288`, `:1147-1184`).
        self.set_the_pcurve(&ee, &ff, e.orientation(), &pct);
    }
}

/// `PCT->Translate(gp_Vec2d(du, dv))` (`cxx:1278`, `:1283`).
fn translate_2d(pc: &Arc<dyn Curve2d>, du: f64, dv: f64) -> Arc<dyn Curve2d> {
    let mut tr = GpTrsf2d::identity();
    tr.set_translation_vec(&GpVec2d::new(du, dv));
    Arc::from(pc.transformed(&tr))
}
