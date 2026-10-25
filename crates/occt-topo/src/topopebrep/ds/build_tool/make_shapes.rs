//! `TopOpeBRepDS_BuildTool` — the `Make*` family.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepDS/TopOpeBRepDS_BuildTool.cxx:150-290`
//! (`MakeVertex`, the four `MakeEdge` overloads, `MakeFace`, `MakeShell`,
//! `MakeSolid`) (OCCT 8.0.0). `MakeShell` / `MakeSolid` / `MakeVertex` stay in
//! `mod.rs`; the rest is here.
//!
//! `BRep_Builder::MakeEdge(E, C, Tol)` is the inline pair
//! `MakeEdge(E)` + `UpdateEdge(E, C, L, Tol)` (`BRep_Builder.lxx:33-38`):
//! a fresh `BRep_TEdge` whose empty curve list receives one
//! `BRep_Curve3D(C, C->FirstParameter(), C->LastParameter())`
//! (`BRep_Builder::UpdateCurves`, `BRep_Builder.cxx:58-102`), then
//! `BRep_TEdge::UpdateTolerance(Tol)` (`BRep_TEdge.lxx:33-37`, raise only).

use std::sync::Arc;

use occt_geom::curve::Curve;

use crate::abs::ShapeType;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;

use super::super::curve::Curve as CurveData;
use super::super::data_structure::DataStructure;
use super::super::surface::Surface as SurfaceData;
use super::BuildTool;

impl BuildTool {
    /// `MakeEdge(E, C)` (`cxx:157-210`).
    ///
    /// `addorigin` (`cxx:172`) and `setrange` (`cxx:173`) are local variables
    /// that are never assigned in OCCT, so both of their blocks
    /// (`cxx:175-209`) are unreachable; only the `MakeEdge` and the
    /// `Degenerated(E, true)` of the null-curve arm survive.
    pub fn make_edge(&self, e: &mut TopoShape, c: &CurveData) {
        let reg = GeometryRegistry::global();
        match c.curve().cloned() {
            None => {
                *e = TopoShape::new(ShapeType::Edge);
                // `BRep_Builder::Degenerated(E, true)` (`BRep_Builder.cxx:1073-1085`).
                reg.set_degenerated(e, true);
                reg.set_curve3d_null(e, true);
            }
            Some(gc) => self.make_edge_curve(e, gc, c.tolerance()),
        }
    }

    /// `MakeEdge(E, C, BDS)` (`cxx:214-246`).
    ///
    /// The null-curve arm builds a degenerated edge carrying the pcurve of
    /// `C.GetSCI1()` on the surface `BDS.Surface(iS)` (`cxx:230-239`). That
    /// pcurve is *surface*-keyed (`BRep_Builder::UpdateEdge(E, PC, Surface, L,
    /// Tol)`, `BRep_Builder.cxx:655-677`), and the port keys every pcurve by
    /// face, so this arm is UNPORTED past the edge creation.
    pub fn make_edge_with_ds(&self, e: &mut TopoShape, c: &CurveData, _bds: &DataStructure) {
        let reg = GeometryRegistry::global();
        match c.curve().cloned() {
            None => {
                *e = TopoShape::new(ShapeType::Edge);
                reg.set_degenerated(e, true);
                reg.set_curve3d_null(e, true);
            }
            Some(gc) => self.make_edge_curve(e, gc, c.tolerance()),
        }
    }

    /// `MakeEdge(E, C, Tol)` (`cxx:249-256`) over `BRep_Builder::MakeEdge(E, C,
    /// Tol)` (`BRep_Builder.lxx:33-38`).
    pub fn make_edge_curve(&self, e: &mut TopoShape, c: Arc<dyn Curve>, tol: f64) {
        let edge = Edge::new();
        *e = edge.0;
        let reg = GeometryRegistry::global();
        reg.set_edge_curve3d(e, c.clone(), tol);
        // `BRep_Curve3D(C, C->FirstParameter(), C->LastParameter())`.
        reg.set_ranges_all(e, c.first_parameter(), c.last_parameter());
    }

    /// `MakeEdge(E)` (`cxx:258-263`) over `BRep_Builder::MakeEdge(TopoDS_Edge&)`
    /// (`BRep_Builder.cxx:623-632`): a fresh `BRep_TEdge` with an empty curve
    /// list and no geometry.
    pub fn make_edge_empty(&self, e: &mut TopoShape) {
        *e = TopoShape::new(ShapeType::Edge);
    }

    /// `MakeFace(F, S)` (`cxx:272-277`) over `BRep_Builder::MakeFace(F,
    /// Surface, Tol)` (`BRep_Builder.cxx:508-536`).
    ///
    /// A null `S.Surface()` leaves a `BRep_TFace` whose surface is null; the
    /// port's `FaceGeom` always carries a surface, so the face is left without
    /// geometry (which is what `BRep_Tool::Surface` answers for both).
    pub fn make_face(&self, f: &mut TopoShape, s: &SurfaceData) {
        match s.surface().cloned() {
            Some(su) => {
                let face = TopoBuilder::new().make_face_with_tolerance(su, &[], s.tolerance());
                *f = face.0;
            }
            None => {
                *f = Face::new().0;
            }
        }
    }
}
