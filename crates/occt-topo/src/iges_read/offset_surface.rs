//! `IGESToBRep_TopoSurface::TransferOffsetSurface` (`cxx:1041-1215`).
//!
//! A finite basis that is already at least C1 becomes `Geom_OffsetSurface`
//! and `BRepLib_MakeFace`. A C0 basis runs `C0ShapeToC1Shape`
//! (`ShapeAlgo_AlgoContainer.cxx:390`). Knot removal that leaves one patch
//! is applied. A patch split runs `ShapeUpgrade_SplitSurface::Build`
//! (`ShapeUpgrade_SplitSurface.cxx:505-618`) and `ShapeFix_ComposeShell`.
//! An infinite elementary basis uses
//! `ShapeConstruct::ConvertSurfaceToBSpline` (`ShapeConstruct.cxx:181-185`).
//! An infinite extrusion converts its generatrix with `ConvertCurveToBSpline`
//! (`:188-235`). Any other infinite basis runs `Iso::MakeApprox` through
//! `UpdateInU` and `UpdateInV`. A finished square is stored, and the loop
//! returns the B-spline from `ConvertBS`. A thrown approximation decrements
//! `GeomAbs_Shape` and tries again (`ShapeConstruct.cxx:301`). A failed iso
//! cuts with `PrefAndRec::Value`. `CutSense` 0 keeps a stored square
//! (`AdvApp2Var_ApproxAFunc2Var.cxx:442`). `D3` stores `D3UUV` or
//! `D3UVV` (`GeomConvert_ApproxSurface.cxx:193`). `mmfmtb1_` transposes
//! the IsoU tables (`AdvApp2Var_ApproxF2var.cxx:5634`). `MakeApprox`
//! passes one `NDIMSE` of 3 (`AdvApp2Var_Iso.cxx:122`). The rational
//! Bezier `Segment` arm runs when weights are stored
//! (`Geom_BezierSurface.cxx:961`). `EvalD0`, `EvalD1`, `EvalD2`, and `EvalDN` use the rational B-spline
//! evaluator (`Geom_BezierSurface.cxx:1425`). `UIso` and `VIso` call
//! `BSplSLib::Iso` (`Geom_BezierSurface.cxx:1776`). A trimmed rational
//! Bezier becomes a rational B-spline (`GeomConvert_1.cxx:714`).
//! Curve `Segment` trims weights with the poles (`Geom_BezierCurve.cxx:395`).
//! `SetWeight` turns a curve or surface rational (`Geom_BezierCurve.cxx:455`,
//! `Geom_BezierSurface.cxx:1179`). `SetWeightCol` and `SetWeightRow`
//! rewrite one column or a row span (`Geom_BezierSurface.cxx:1219`,
//! `:1262`). `UReverse` and `VReverse` swap poles and weights
//! (`:1305`, `:1352`). `ExchangeUV` transposes poles and weights
//! (`Geom_BezierSurface.cxx:490`). `Increase` elevates degree through
//! `BSplSLib::IncreaseDegree` (`:532`). Insert and remove pole columns
//! and rows (`Geom_BezierSurface.cxx:652-937`). `SetPole` and column/row
//! setters (`:1047-1175`). Weighted construction (`:416`). Curve
//! `Increase`, insert/remove poles, and `Reverse`
//! (`Geom_BezierCurve.cxx:141`, `:210`, `:298`, `:352`).

use std::sync::Arc;

use occt_core::bspl::locate::{first_u_knot_index, last_u_knot_index};
use occt_core::precision::{CONFUSION, Precision};
use occt_core::toploc::TopLocLocation;
use occt_geom::{GeomOffsetSurface, Surface};

use crate::brep_tool::BRepTool;
use crate::iterator::cumulated_children;
use crate::shape::{Face, TopoShape};

use super::iges_data_iges_model::EntityData;
use super::iges_to_brep_basic_curve::{convert_location, entity_location};
use super::iges_to_brep_curve_and_surface::{is_topo_surface, Session};
use super::iges_to_brep_topo_surface::transfer_topo_surface;

pub(in crate::iges_read) fn transfer_offset_surface(
    session: &mut Session<'_>,
    index: usize,
) -> Result<TopoShape, String> {
    let (surface_index, distance) = match &session.model.entity(index)?.data {
        EntityData::OffsetSurface {
            surface, distance, ..
        } => (*surface, *distance),
        _ => return Err("IGES 140 has no offset-surface parameters".into()),
    };
    if !is_topo_surface(session.model, surface_index) {
        return Err("XSTEP_164: IGESToBRep_TopoSurface.cxx:1060-1066".into());
    }
    let iges_shape = match transfer_topo_surface(session, surface_index) {
        Ok(shape) => shape,
        Err(_) => {
            return Err("IGES_1156: IGESToBRep_TopoSurface.cxx:1070-1078".into());
        }
    };
    let face = basis_face(&iges_shape)?;
    // `BRep_Tool::Surface(face)` one-arg (`BRep_Tool.cxx`) copies the surface
    // through the face location. The two-arg form is commented out at cxx:1110.
    let mut geom = BRepTool::face_surface_world(&face)
        .ok_or_else(|| "BRep_Tool::Surface: IGESToBRep_TopoSurface.cxx:1115".to_string())?;
    let dist = distance * session.model.unit_factor;
    let basis: Arc<dyn Surface> = if geom.is_offset_surface() {
        // cxx:1121 calls `basisSrf->Offset()` while `basisSrf` is still null.
        // The live value is `geom140->Offset()`.
        let prev = geom.offset_distance().unwrap_or(0.0);
        let inner = geom.offset_basis_surface().ok_or_else(|| {
            "Geom_OffsetSurface::BasisSurface: IGESToBRep_TopoSurface.cxx:1120".to_string()
        })?;
        Arc::new(GeomOffsetSurface::new(inner, prev + dist))
    } else {
        if is_geom_abs_c0(geom.as_ref()) {
            // `C0ShapeToC1Shape` (`cxx:1126-1154`). A result that is not a
            // face is returned (`cxx:1136-1140`). A face that is still C0
            // is nullified (`cxx:1142-1149`).
            match crate::shape_upgrade::correct_c0_basis(&face, geom, dist.abs())? {
                crate::shape_upgrade::C0Basis::Shape(shape) => return Ok(shape),
                crate::shape_upgrade::C0Basis::Surface(surf) => {
                    geom = surf;
                    if is_geom_abs_c0(geom.as_ref()) {
                        return Err("IGES_1266: IGESToBRep_TopoSurface.cxx:1142-1149".into());
                    }
                }
            }
        }
        let (umin, umax) = geom.u_range();
        let (vmin, vmax) = geom.v_range();
        if Precision::is_infinite(umin)
            || Precision::is_infinite(umax)
            || Precision::is_infinite(vmin)
            || Precision::is_infinite(vmax)
        {
            // `BRepTools::UVBounds(face, umin, umax, vmin, vmax)`, then a
            // `Geom_RectangularTrimmedSurface` and `ConvertSurfaceToBSpline`
            // (`TopoSurface.cxx:1166-1181`):
            //
            //   BRepTools::UVBounds(face, umin, umax, vmin, vmax);
            //   Handle(Geom_RectangularTrimmedSurface) TS =
            //     new Geom_RectangularTrimmedSurface(geomSupport, umin, umax, vmin, vmax);
            //   Handle(Geom_BSplineSurface) BS =
            //     ShapeAlgo::AlgoContainer()->ConvertSurfaceToBSpline(TS, umin, umax, vmin, vmax);
            //
            // An unbounded basis reaches this point with infinite bounds because
            // `BRepTools::AddUVBounds` falls back to the surface's natural
            // `Bounds` for a face without edges (`BRepTools.cxx:139-153`), and
            // `Geom_Plane` answers `Bounds() == +-Precision::Infinite()`
            // (`Geom_Plane.cxx:178-185`).
            //
            // OCCT does not survive that call either. `ConvertSurfaceToBSpline`
            // peels the trim, finds a `Geom_ElementarySurface`, rebuilds
            // `Geom_RectangularTrimmedSurface(S, UF, UL, VF, VL)` and calls
            // `GeomConvert::SurfaceToBSplineSurface` (`ShapeConstruct.cxx:181-185`),
            // which opens with
            //   if (Precision::IsNegativeInfinite(UFirst) || Precision::IsPositiveInfinite(ULast)
            //       || Precision::IsNegativeInfinite(VFirst) || Precision::IsPositiveInfinite(VLast))
            //     throw Standard_DomainError(...);
            // (`GeomConvert_1.cxx:346-358`). The +-Precision::Infinite() bounds
            // taken from the natural range satisfy that test, so the entity
            // throws and its transfer produces nothing.
            //
            // The port reaches the same outcome one line earlier: the guard below
            // (`TopoSurface.cxx:1169`) returns the error the caller would have
            // seen as the thrown `Standard_DomainError`. It must stay early
            // rather than call `convert_surface_to_bspline`: the port's
            // `Surface::u_range`/`v_range` for the unbounded elementary surfaces
            // answer IEEE infinity (`occt-geom/src/plane.rs:42`,
            // `cylinder.rs:49`, `cone.rs:58`, `surface_of_linear_extrusion.rs:160`,
            // `line.rs:26`, `parabola.rs:26`, `hyperbola.rs:27`), and
            // `plane_patch` (`gtransform.rs:496-506`) would read that as a finite
            // 2.e+100-style range and return a gigantic patch instead of failing.
            //
            // Aligning those sentinels to `Precision::INFINITE` (what
            // `Geom_Surface::Bounds` returns) is a separate, wider migration:
            // applied to `GeomPlane` alone it takes `cargo test -p occt-topo
            // --lib` from 1281 passed to 1170/111 failed, because the port has
            // ~112 `f64::is_finite()` tests that use IEEE infinity to mean
            // "unbounded" (e.g. `bop_build_common::tests::box_shell_wraps_into_unit_solid`
            // reads a unit-box volume of 1.6e+201, `bean_face` range queries
            // return none). Tracked, not needed for this file.
            //
            // Verified on `_igesdump/oracle/OffsetPlaneHoleEdge.iges` (entity 140
            // over the unbounded plane entity 108): OCCT's 402 form-1 group has
            // 10 members (`402,10,3,53,81,97,113,129,145,161,177,193,...`) and
            // its result has 9 faces (`occt_probe`, faces=9, BRepCheck valid=1,
            // volume 134.4) -- exactly one member, this 140 over the unbounded
            // plane, produces nothing. The port now returns the same 9 faces.
            let (fu, lu, fv, lv) = crate::brep_tools::uv_bounds(&face);
            if Precision::is_infinite(fu)
                || Precision::is_infinite(lu)
                || Precision::is_infinite(fv)
                || Precision::is_infinite(lv)
            {
                return Err(
                    "Geom_RectangularTrimmedSurface: IGESToBRep_TopoSurface.cxx:1169".into(),
                );
            }
            let bs = super::gtransform::convert_surface_to_bspline(geom.as_ref(), fu, lu, fv, lv)?;
            if is_geom_abs_c0(&bs) {
                return Err("IGES_1265: IGESToBRep_TopoSurface.cxx:1174-1176".into());
            }
            geom = Arc::new(bs);
        }
        Arc::new(GeomOffsetSurface::new(geom, dist))
    };
    let built = crate::brep_lib_make_face::make_face_from_surface(basis, CONFUSION);
    // `HasTransf` uses `CompoundLocation` (`cxx:1196-1206`). This applies the
    // entity's own 124, the same gap as `finish_plane`.
    apply_location(session, index, built.0)
}

fn basis_face(shape: &TopoShape) -> Result<Face, String> {
    if shape.is_face() {
        return Face::wrap(shape.clone())
            .ok_or_else(|| "TopoDS::Face: IGESToBRep_TopoSurface.cxx:1085".to_string());
    }
    if shape.is_shell() {
        let kids = cumulated_children(shape);
        if let Some(first) = kids.into_iter().next() {
            if first.is_face() {
                return Face::wrap(first)
                    .ok_or_else(|| "TopoDS::Face: IGESToBRep_TopoSurface.cxx:1093".to_string());
            }
        }
    }
    Err("IGES_1156: IGESToBRep_TopoSurface.cxx:1099-1106".into())
}

/// `Geom_Surface::Continuity() == GeomAbs_C0`.
///
/// B-spline smoothness is `degree - MaxKnotMult` (`Geom_BSplineSurface.cxx:1167-1190`),
/// and `MaxKnotMult` is the max interior multiplicity (`BSplCLib.cxx:740-751`).
/// A swept surface uses the basis curve (`Geom_SweptSurface.cxx:26-28`).
/// The surface trait's `continuity` on a B-spline is the degree, so it is not used.
fn is_geom_abs_c0(surf: &dyn Surface) -> bool {
    if let Some(basis) = surf.rectangular_trimmed_basis() {
        return is_geom_abs_c0(basis.as_ref());
    }
    if surf.is_bspline_surface() {
        let Some(bs) = surf.osculating_bspline() else {
            return false;
        };
        let (_, um) = bs.distinct_knots_and_mults_u();
        let (_, vm) = bs.distinct_knots_and_mults_v();
        return direction_is_c0(bs.deg_u, bs.u_periodic, &um)
            || direction_is_c0(bs.deg_v, bs.v_periodic, &vm);
    }
    if surf.is_surface_of_revolution() || surf.is_surface_of_linear_extrusion() {
        return surf.continuity() == 0;
    }
    false
}

fn direction_is_c0(degree: usize, periodic: bool, mults: &[i32]) -> bool {
    let first = if periodic {
        1
    } else {
        first_u_knot_index(degree as i32, mults)
    };
    let last = if periodic {
        mults.len() as i32
    } else {
        last_u_knot_index(degree as i32, mults)
    };
    if last - first == 1 {
        return false;
    }
    let mut max_mult = 0i32;
    for i in (first + 1)..last {
        let idx = (i as usize).saturating_sub(1);
        if let Some(m) = mults.get(idx) {
            max_mult = max_mult.max(*m);
        }
    }
    degree as i32 - max_mult == 0
}

fn apply_location(
    session: &Session<'_>,
    index: usize,
    mut shape: TopoShape,
) -> Result<TopoShape, String> {
    if let Some(g) = entity_location(session.model, index)? {
        if let Some(tr) = convert_location(1.0e-4, &g, session.model.unit_factor) {
            if tr.form() != occt_core::gp::TrsfForm::Identity {
                shape.move_location(&TopLocLocation::composed(
                    &tr,
                    TopLocLocation::identity(),
                ));
            }
        }
    }
    Ok(shape)
}
