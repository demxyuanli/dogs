//! `BRepOffset_Tool::TryProject` (`BRepOffset_Tool.cxx:1945-2006`): for every
//! edge of `edges`, test whether it already lies on `f1`'s surface; on a hit
//! attach its pcurve to `f1`, orient it between the two faces and append it to
//! the two output lists.
//!
//! The returned flag is the outer `Ok`: it starts `true` (`cxx:1958`), any edge
//! that misses the surface sets it to `false` (`cxx:2001-2004`), and the loop
//! does not stop. `LInt1` / `LInt2` are cleared up front (`cxx:1955-1956`) and
//! hold only the hits.
//!
//! `BRep_Tool::Curve(E, L, f, l)` is [`BRepTool::edge_curve_world`], which
//! already applies the edge location, so the `C->Transform(L)` of
//! `cxx:1981-1984` has no counterpart here; the `Geom_TrimmedCurve(C, f, l)` of
//! `cxx:1979` is the [`trimmed_view`] handed to `IsOnSurface`.

use std::sync::Arc;

use occt_core::kernel::geomabs::Shape;
use occt_geom::{Curve, Surface};

use crate::abs::State;
use crate::brep_lib_build_curves3d::build_curves3d_tol;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;

use super::curve_checks::is_on_surface;
use super::pcurve_builder::{build_pcurves, trimmed_view};
use super::BRepOffsetTool;

impl BRepOffsetTool {
    /// `TryProject(F1, F2, Edges, LInt1, LInt2, Side, TolConf)`
    /// (`cxx:1945-2006`).
    ///
    /// `l_int1` / `l_int2` are the two out-lists (cleared first); the return
    /// value is OCCT's `Ok`.
    pub fn try_project(
        f1: &Face,
        f2: &Face,
        edges: &[Edge],
        l_int1: &mut Vec<Edge>,
        l_int2: &mut Vec<Edge>,
        side: State,
        tol_conf: f64,
    ) -> bool {
        l_int1.clear();
        l_int2.clear();

        // Port guard: `IsOnSurface(C, Bouchon, ...)` dereferences the surface
        // handle built at `cxx:1961`.
        let Some(bouchon) = BRepTool::face_surface(f1) else {
            return false;
        };
        let mut ok = true;

        for cur_e in edges {
            // `BRep_Tool::Curve(CurE, L, f, l)` (`cxx:1969`).
            let mut curve = BRepTool::edge_curve_world(cur_e);
            if curve.is_none() {
                // `cxx:1971-1978`: rebuild the 3D curve, then re-read it; an
                // edge that still has none is skipped (`continue`).
                let tol = BRepTool::edge_tolerance(cur_e);
                build_curves3d_tol(&cur_e.0, tol, Shape::C1, 14, 0);
                curve = BRepTool::edge_curve_world(cur_e);
            }
            let Some(curve) = curve else {
                continue;
            };
            if !Self::try_project_one(
                cur_e,
                curve,
                bouchon.clone(),
                f1,
                f2,
                l_int1,
                l_int2,
                side,
                tol_conf,
            ) {
                ok = false;
            }
        }
        ok
    }

    /// The `TryProject` loop body for one edge that has a 3D curve: the
    /// `IsOnSurface` test and the hit bookkeeping of `cxx:1984-2005`.
    #[allow(clippy::too_many_arguments)]
    fn try_project_one(
        cur_e: &Edge,
        c: Arc<dyn Curve>,
        bouchon: Arc<dyn Surface>,
        f1: &Face,
        f2: &Face,
        l_int1: &mut Vec<Edge>,
        l_int2: &mut Vec<Edge>,
        side: State,
        tol_conf: f64,
    ) -> bool {
        let (f, l) = BRepTool::edge_parameters(cur_e);
        // `C = new Geom_TrimmedCurve(C, f, l)` (`cxx:1979`).
        let trimmed = trimmed_view(c, f, l);

        // `cxx:1983-1985`: a miss is the outer `Ok = false`.
        let Some(tol_reached) = is_on_surface(trimmed.as_ref(), bouchon, tol_conf) else {
            return false;
        };

        // `B.UpdateEdge(CurE, TolReached)` (`cxx:1989`).
        GeometryRegistry::global().set_edge_tolerance(&cur_e.0, tol_reached);
        // `BuildPCurves(CurE, F1)` (`cxx:1990`): its only failure channel is
        // the `Standard_ConstructionError` of `cxx:480`, which OCCT lets
        // escape; the port drops the pcurve and keeps the hit.
        let _ = build_pcurves(cur_e, f1);
        // `OrientSection(CurE, F1, F2, O1, O2)` (`cxx:1991`).
        let Some((mut o1, mut o2)) = Self::orient_section(cur_e, f1, f2) else {
            return false;
        };
        // `cxx:1992-1996`.
        if side == State::Out {
            o1 = o1.reversed();
            o2 = o2.reversed();
        }
        // `cxx:1997-1998`.
        l_int1.push(Edge(cur_e.0.oriented(o1)));
        l_int2.push(Edge(cur_e.0.oriented(o2)));
        true
    }
}
