//! `ChFi3d_BuildPlane` — the tangent plane that replaces a support face whose
//! contact point has no face left to walk to.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:1113-1144`
//! (OCCT 8.0.0).
//!
//! When `PerformSetOfSurfOnElSpine` runs out of connected faces at a spine
//! extremity, it builds the plane tangent to the surface at the last known UV
//! point and continues the band on that plane (`_2.cxx:2514-2523`). The plane's
//! face carries the orientation of the face it replaces, and the UV point is
//! reset to `(0, 0)` because the plane is its own parameter space.
//!
//! Divergences that are not heuristics:
//! - `BRepLib_MakeFace(Pln, Precision::Confusion())` is
//!   [`TopoBuilder::make_face_plane`], the same whole-plane face the OCCT
//!   builder produces for a `gp_Pln`,
//! - `BRepLProp_SLProps(*HS, u, v, 1, 1.e-12)` is [`SlProps::new`] on the
//!   adaptor's surface, then `Value()` / `Normal()` with the same arguments,
//! - a null `DStr.Shape(SD->Index(ons))` or a shape that is not a face answers
//!   `Ok(())` without touching `hs`, exactly the `F.IsNull()` early return of
//!   `cxx:1122-1125`,
//! - the two dereferences OCCT takes without a guard — the null p-curve of
//!   `BRep_Tool::CurveOnSurface` (`cxx:1129-1131`) and the null normal of the
//!   props — answer the `Standard_Failure("ChFi3d_BuildPlane : echec .")` of
//!   `cxx:1143` here.

use occt_core::gp::{GpAx3, GpPln, GpPnt2d};
use occt_geom::geom_lprop::SlProps;

use crate::boptools_2d::curve_on_surface_range;
use crate::builder::TopoBuilder;
use crate::chfi_ds::SurfData;
use crate::shape::{Edge, Face};
use crate::topopebrep::ds::HDataStructure;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `ChFi3d_BuildPlane(DStr, HS, pons, SD, isfirst, ons)`
/// (`_2.cxx:1113-1144`).
///
/// `hs` is re-pointed at the tangent-plane face and `pons` reset to `(0, 0)`
/// on success; both are left untouched when the data structure holds no face
/// for the interference.
pub fn ch_fi3d_build_plane(
    ds: &HDataStructure,
    hs: &mut SurfaceAdaptor,
    pons: &mut GpPnt2d,
    sd: &SurfData,
    is_first: bool,
    ons: i32,
) -> Result<(), String> {
    // `const TopoDS_Face F = TopoDS::Face(DStr.Shape(SD->Index(ons)));
    // if (F.IsNull()) return;` (`cxx:1120-1125`).
    let Some(shape) = ds.shape(sd.index(ons), true) else {
        return Ok(());
    };
    let Some(f) = Face::wrap(shape) else {
        return Ok(());
    };

    // `if (SD->Vertex(isfirst, ons).IsOnArc())` (`cxx:1127`).
    let common_point = sd.vertex(is_first, ons);
    if common_point.is_on_arc() {
        let Some(arc) = common_point.arc().cloned() else {
            return Err("ChFi3d_BuildPlane : echec .".into());
        };
        let Some(arc) = Edge::wrap(arc) else {
            return Err("ChFi3d_BuildPlane : echec .".into());
        };
        let Some(param) = common_point.parameter_on_arc() else {
            return Err("ChFi3d_BuildPlane : echec .".into());
        };
        // `BRep_Tool::CurveOnSurface(...)` then `Hc->Value(Param).Coord(u, v)`
        // (`cxx:1128-1132`).
        let Some((pc, _, _)) = curve_on_surface_range(&arc, &f) else {
            return Err("ChFi3d_BuildPlane : echec .".into());
        };
        let uv = pc.d0(param);

        // `BRepLProp_SLProps theProp(*HS, u, v, 1, 1.e-12)` (`cxx:1133`).
        let Some(surface) = hs.surface().cloned() else {
            return Err("ChFi3d_BuildPlane : echec .".into());
        };
        let mut props = SlProps::new(surface.as_ref(), uv.x(), uv.y(), 1, 1.0e-12)
            .map_err(|_| "ChFi3d_BuildPlane : echec .".to_string())?;

        // `if (theProp.IsNormalDefined())` (`cxx:1134`).
        if props.is_normal_defined(surface.as_ref()) {
            // `new Geom_Plane(theProp.Value(), theProp.Normal())`
            // (`cxx:1136`): `gp_Pln(P, N)` is `gp_Ax3(P, N)`.
            let point = props.value();
            let normal = props
                .normal(surface.as_ref())
                .map_err(|_| "ChFi3d_BuildPlane : echec .".to_string())?;
            let Ok(ax3) = GpAx3::from_point_dir(point, normal) else {
                return Err("ChFi3d_BuildPlane : echec .".into());
            };
            let plane = GpPln::new(ax3);
            // `TopoDS_Face NewF = BRepLib_MakeFace(Pln,
            // Precision::Confusion()); NewF.Orientation(F.Orientation());`
            // (`cxx:1137-1138`).
            let builder = TopoBuilder::new();
            let mut new_face = builder.make_face_plane(&plane);
            new_face.0.set_orientation(f.0.orientation());
            // `pons.SetCoord(0., 0.); HS->Initialize(NewF); return;`
            // (`cxx:1139-1141`).
            pons.set_coord(0.0, 0.0);
            hs.initialize(&new_face);
            return Ok(());
        }
    }

    // `throw Standard_Failure("ChFi3d_BuildPlane : echec .")` (`cxx:1143`).
    Err("ChFi3d_BuildPlane : echec .".into())
}
