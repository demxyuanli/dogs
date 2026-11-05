//! `BRepTools_TrsfModification` (`BRepTools_TrsfModification.hxx:38-146`,
//! `BRepTools_TrsfModification.cxx:42-451`): the `gp_Trsf` driven modification,
//! and the reference implementation of [`Modification`].
//!
//! Geometry frames: the reference asks `BRep_Tool::Surface(F, L)` /
//! `BRep_Tool::Curve(E, L, f, l)` for the geometry **and** its location `L` and
//! transforms it by `L.Transformation().Inverted() * myTrsf * L.Transformation()`
//! (`cxx:85-88`, `cxx:286-289`). This port keeps the geometry of a shape in the
//! shape's own frame, so `L` here is the shape's own [`TopLocLocation`] and that
//! is exactly what is handed back: the modifier then finds
//! `myLoc.Predivided(S.Location())` to be the identity
//! (`BRepTools_Modifier.cxx:261`) and keeps the shape's frame, which is also
//! why the pcurves written by `NewCurve2d` stay valid in it.
//!
//! UNPORTED:
//! - `NewTriangulation` (`cxx:95-164`): the `myCopyMesh == true` path (`Copy`,
//!   `SetNode`, `TransformParameters`, mirrored triangles, `SetNormal`) is not
//!   ported; the trait default answers "not modified", which is the reference's
//!   answer while `myCopyMesh` is false (its constructor value, `cxx:42-46`).
//! - `NewPolygon` / `NewPolygonOnTriangulation` (`cxx:168-217`, `cxx:221-272`):
//!   the port keeps no edge polygons (see `modification.rs`).
//! - `NewCurve2d` (`cxx:313-409`) always answers "no new curve on surface". The
//!   reference re-parametrizes the pcurve through
//!   `Geom_Surface::ParametricTransformation` + `Geom2d_TrimmedCurve` +
//!   `GeomLib::GTransform` (`cxx:383-388`) and then hands it to
//!   `GeomLib::SameRange` (`cxx:407`). None of those three exist in this port
//!   (`Surface` has no `parametric_transformation`, and `GeomLib::GTransform` /
//!   `GeomLib::SameRange` are not ported), so no branch could produce the
//!   reference's result curve. The three exact "no new curve" early returns —
//!   null surface (`cxx:326-330`), planar surface (`cxx:331-335`), missing pcurve
//!   (`cxx:337-341`) — answer `false` in the reference as well, so the gap only
//!   shows on a non-planar surface under a transform that moves its parametric
//!   frame.

use std::sync::Arc;

use occt_core::gp::GpTrsf;
use occt_core::kernel::geomabs::Shape as GeomAbsShape;
use occt_core::toploc::TopLocLocation;
use occt_geom2d::curve::Curve2d;

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;

use super::modification::{
    geom_abs_shape_from_u8, Modification, NewCurveInfo, NewPointInfo, NewParameterInfo,
    NewSurfaceInfo,
};

/// `BRepTools_TrsfModification` (`BRepTools_TrsfModification.hxx:38-146`).
#[derive(Debug, Clone)]
pub struct TrsfModification {
    /// `myTrsf` (`BRepTools_TrsfModification.hxx:145`).
    my_trsf: GpTrsf,
    /// `myCopyMesh` (`BRepTools_TrsfModification.hxx:146`).
    my_copy_mesh: bool,
}

impl TrsfModification {
    /// `BRepTools_TrsfModification(const gp_Trsf& T)` (`cxx:42-46`).
    pub fn new(trsf: GpTrsf) -> Self {
        Self { my_trsf: trsf, my_copy_mesh: false }
    }

    /// `Trsf()` (`cxx:51-54`).
    pub fn trsf(&self) -> &GpTrsf {
        &self.my_trsf
    }

    /// `gp_Trsf& Trsf()` (`cxx:51-54`): the reference hands out a writable
    /// reference, so the port exposes the same handle.
    pub fn trsf_mut(&mut self) -> &mut GpTrsf {
        &mut self.my_trsf
    }

    /// `bool& IsCopyMesh()` (`cxx:59-62`) read back.
    pub fn is_copy_mesh(&self) -> bool {
        self.my_copy_mesh
    }

    /// The writable `bool& IsCopyMesh()` (`cxx:59-62`).
    pub fn set_copy_mesh(&mut self, v: bool) {
        self.my_copy_mesh = v;
    }

    /// `LT = L.Transformation().Inverted() * myTrsf * L.Transformation()`
    /// (`cxx:85-88`, `cxx:286-289`).
    fn local_trsf(&self, loc: &TopLocLocation) -> GpTrsf {
        if loc.is_identity() {
            return self.my_trsf.clone();
        }
        let l = loc.transformation();
        let inv = l.inverted().unwrap_or_else(|_| GpTrsf::identity());
        inv.multiplied(&self.my_trsf).multiplied(&l)
    }

    /// `Tol *= std::abs(myTrsf.ScaleFactor())` (`cxx:79-81`, `cxx:283-284`,
    /// `cxx:306-307`, `cxx:321-323`).
    fn scaled(&self, tol: f64) -> f64 {
        tol * self.my_trsf.scale_factor().abs()
    }
}

impl Modification for TrsfModification {
    /// `NewSurface` (`cxx:63-92`).
    fn new_surface(&self, f: &TopoShape) -> Option<NewSurfaceInfo> {
        // `S = BRep_Tool::Surface(F, L); if (S.IsNull()) return false;`
        // (`cxx:68-72`).
        let surface = GeometryRegistry::global().face_surface(f)?;
        let face = Face::wrap(f.clone())?;
        // `Tol = BRep_Tool::Tolerance(F) * |ScaleFactor|` (`cxx:79-81`).
        let tol = self.scaled(BRepTool::face_tolerance(&face));
        // `BRep_Tool::Surface(F, L)` reports the face's own location; the port
        // keeps no separate surface location, so this is it.
        let loc = f.location().clone();
        let lt = self.local_trsf(&loc);
        Some(NewSurfaceInfo {
            // `S = S->Transformed(LT)` (`cxx:90`).
            surface: Arc::from(surface.transformed(&lt)),
            location: loc,
            tolerance: tol,
            // `RevWires = false` (`cxx:82`).
            rev_wires: false,
            // `RevFace = myTrsf.IsNegative()` (`cxx:83`).
            rev_face: self.my_trsf.is_negative(),
        })
    }

    /// `NewCurve` (`cxx:276-297`): always `true`; a null curve stays null
    /// (`cxx:291-294` leaves `C` null and only the tolerance/L are written).
    fn new_curve(&self, e: &TopoShape) -> Option<NewCurveInfo> {
        let edge = Edge::wrap(e.clone())?;
        let reg = GeometryRegistry::global();
        // `Tol = BRep_Tool::Tolerance(E) * |ScaleFactor|` (`cxx:283-284`).
        let tol = self.scaled(BRepTool::edge_tolerance(&edge));
        let loc = e.location().clone();
        let lt = self.local_trsf(&loc);
        // `if (!C.IsNull()) C = C->Transformed(LT)` (`cxx:291-294`).
        let curve = reg.edge_curve(e).map(|c| Arc::from(c.transformed(&lt)));
        Some(NewCurveInfo { curve, location: loc, tolerance: tol })
    }

    /// `NewPoint` (`cxx:302-309`): the vertex point is already in world
    /// coordinates, so it takes the global `myTrsf`.
    fn new_point(&self, v: &TopoShape) -> Option<NewPointInfo> {
        let vertex = crate::shape::Vertex::wrap(v.clone())?;
        // `P = BRep_Tool::Pnt(V)` (`cxx:304`): `BRep_Tool.cxx:378-386` applies
        // the vertex location.
        let p = BRepTool::vertex_point_world(&vertex);
        // `Tol = BRep_Tool::Tolerance(V) * |ScaleFactor|` (`cxx:306-307`).
        let tol = self.scaled(BRepTool::vertex_tolerance(&vertex));
        // `P.Transform(myTrsf)` (`cxx:308`).
        Some(NewPointInfo { point: p.transformed(&self.my_trsf), tolerance: tol })
    }

    /// `NewCurve2d` (`cxx:313-409`).
    ///
    /// UNPORTED — see the module documentation:
    /// `ParametricTransformation`, `GeomLib::GTransform` and
    /// `GeomLib::SameRange` are not ported, so this answers the reference's
    /// `false` ("no new curve on surface").
    fn new_curve2d(
        &self,
        e: &TopoShape,
        f: &TopoShape,
        new_e: &TopoShape,
        new_f: &TopoShape,
    ) -> Option<Arc<dyn Curve2d>> {
        let _ = (e, f, new_e, new_f);
        None
    }

    /// `NewParameter` (`cxx:413-437`).
    fn new_parameter(&self, v: &TopoShape, e: &TopoShape) -> Option<NewParameterInfo> {
        // `if (V.IsNull()) return false;` (`cxx:419-422`): an infinite edge may
        // have a null vertex. This port has no null shape (`TopoShape::is_null`
        // is always false), so the test has no counterpart.
        let vertex = crate::shape::Vertex::wrap(v.clone())?;
        // `Tol = BRep_Tool::Tolerance(V) * |ScaleFactor|` (`cxx:425-426`).
        let tol = self.scaled(BRepTool::vertex_tolerance(&vertex));
        // `P = BRep_Tool::Parameter(V, E)` (`cxx:427`).
        let mut p = BRepTool::parameter_of_vertex(e, v);
        // `C = BRep_Tool::Curve(E, loc, f, l); if (!C.IsNull())
        //  P = C->TransformedParameter(P, myTrsf);` (`cxx:429-435`): the
        // **global** trsf, unlike the `LT` of the geometry writes.
        if let Some(c) = GeometryRegistry::global().edge_curve(e) {
            p = c.transformed_parameter(p, &self.my_trsf);
        }
        Some(NewParameterInfo { parameter: p, tolerance: tol })
    }

    /// `Continuity` (`cxx:441-451`): `BRep_Tool::Continuity(E, F1, F2)` of the
    /// **source** edge, i.e. the stored regularity representation
    /// (`BRep_Builder::Continuity`, `BRep_Builder.cxx:1012-1021`). The three
    /// `New*` arguments are ignored by the reference.
    fn continuity(
        &self,
        e: &TopoShape,
        f1: &TopoShape,
        f2: &TopoShape,
        new_e: &TopoShape,
        new_f1: &TopoShape,
        new_f2: &TopoShape,
    ) -> GeomAbsShape {
        let _ = (new_e, new_f1, new_f2);
        let (Some(edge), Some(face1), Some(face2)) =
            (Edge::wrap(e.clone()), Face::wrap(f1.clone()), Face::wrap(f2.clone()))
        else {
            return GeomAbsShape::C0;
        };
        geom_abs_shape_from_u8(crate::brep_lib_continuity::stored_continuity(&edge, &face1, &face2))
    }
}
