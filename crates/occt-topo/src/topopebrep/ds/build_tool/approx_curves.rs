//! `TopOpeBRepDS_BuildTool::ApproxCurves`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepDS/TopOpeBRepDS_BuildTool.cxx`
//! (OCCT 8.0.0): `ApproxCurves` `:471-653`.
//!
//! The member re-fits the degree-1 walking approximation of an intersection
//! curve on the two faces through `TopOpeBRepTool_CurveTool::MakeCurves`, and
//! leaves the old degree-1 curves in place when that fit fails
//! (`cxx:513-521`, "MSV Nov 12, 2001").
//!
//! UNPORTED: nothing.

use std::sync::Arc;

use occt_geom::curve::Curve;
use occt_geom2d::curve::Curve2d;

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::tool::curve_tool::CurveTool;

use super::super::curve::Curve as CurveData;
use super::BuildTool;

impl BuildTool {
    /// `ApproxCurves(C, E, inewC, HDS)` (`cxx:471-653`).
    ///
    /// `inewC` is `HDS->MakeCurve(C, newC1)`'s answer; the curves are stored
    /// into `HDS->ChangeCurve(inewC)`.
    pub fn approx_curves(
        &self,
        c: &CurveData,
        e: &TopoShape,
        inew_c: &mut i32,
        hds: &mut HDataStructure,
    ) {
        // cxx:476-478
        let mut new_c1 = CurveData::new();
        *inew_c = hds
            .make_curve(c, &mut new_c1)
            .expect("TopOpeBRepDS_HDataStructure::MakeCurve");

        // cxx:482-483: `TopoDS::Face(newC.Shape1())` / `(Shape2())`. A missing
        // shape leaves OCCT with a null face, whose `BRepAdaptor_Surface`
        // carries no surface and whose `BRepTools::UVBounds` is undefined; the
        // port needs both faces to build the adaptors, so it stops here.
        let (s1, s2) = {
            let new_c = hds.change_curve(*inew_c);
            (new_c.shape1().cloned(), new_c.shape2().cloned())
        };
        let (Some(f1), Some(f2)) = (
            s1.and_then(Face::wrap),
            s2.and_then(Face::wrap),
        ) else {
            return;
        };

        // cxx:485-487
        let c3d = c.curve().cloned();
        let pc1 = c.curve1();
        let pc2 = c.curve2();

        // cxx:492-495: Vmin, Vmax = bounding vertices of edge E and their
        // parameters. `GetOrientedEdgeVertices` leaves them at `0.` when a
        // vertex is null.
        let edge = Edge::wrap(e.clone());
        let (mut parmin, mut parmax) = (0.0, 0.0);
        if let Some(edge) = edge.as_ref() {
            if let Some((_, _, Some((pmin, pmax)))) = self.get_oriented_edge_vertices(edge) {
                parmin = pmin;
                parmax = pmax;
            }
        }

        // cxx:497-509: `myCurveTool.MakeCurves(...)` over `myCurveTool.myGeomTool`,
        // which `BuildTool` keeps as its own `geom_tool`.
        let curve_tool = CurveTool::with_geom_tool(self.geom_tool);
        let outcome =
            curve_tool.make_curves(parmin, parmax, c3d.clone(), pc1.clone(), pc2.clone(), &f1, &f2);

        // cxx:511-533
        let newtol;
        let newparmin;
        let newparmax;
        let c3dnew: Option<Arc<dyn Curve>>;
        let pc1new: Option<Arc<dyn Curve2d>>;
        let pc2new: Option<Arc<dyn Curve2d>>;
        if !outcome.done {
            // MSV Nov 12, 2001: if approx failed then leave old curves of
            // degree 1. OCCT re-reads `C3D` / `PC1` / `PC2` here, i.e. the
            // handles `MakeCurves` may have rewritten through
            // `CurveImprovement`'s `const_cast` (`cxx:516-530`).
            newtol = match edge.as_ref() {
                Some(edge) => BRepTool::edge_tolerance(edge),
                None => 0.0,
            };
            newparmin = parmin;
            newparmax = parmax;
            c3dnew = outcome.improved_c3d.clone().or(c3d);
            pc1new = outcome.improved_pc1.clone().or(pc1);
            pc2new = outcome.improved_pc2.clone().or(pc2);
        } else {
            // cxx:523-533: `UpdateEdgeCurveTol(F1, F2, E, C3Dnew,
            // tolreached3d, tolreached2d, tolreached2d, ...)`. `F1`, `F2` and
            // the three tolerances are unnamed in OCCT's definition and never
            // read (`cxx:343-356`), so the port passes only the edge and the
            // new curve.
            let (t, pmin, pmax) = self.update_edge_curve_tol(e, outcome.c3d_new.clone());
            newtol = t;
            newparmin = pmin;
            newparmax = pmax;
            c3dnew = outcome.c3d_new.clone();
            pc1new = outcome.pc1_new.clone();
            pc2new = outcome.pc2_new.clone();
        }

        // cxx:535-543
        let new_c = hds.change_curve(*inew_c);
        if let Some(c3dnew) = c3dnew {
            new_c.define_curve(Some(c3dnew), newtol, false);
            new_c.set_range(newparmin, newparmax);
        }
        if let Some(pc1new) = pc1new {
            new_c.set_curve1(Some(pc1new));
        }
        if let Some(pc2new) = pc2new {
            new_c.set_curve2(Some(pc2new));
        }
    }
}
