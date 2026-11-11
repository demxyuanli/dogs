//! The file-local statics of `BRepFeat_MakeCylindricalHole`.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeCylindricalHole.cxx`:
//!
//! | Reference | Rust |
//! |---|---|
//! | `Baryc` (`cxx:690-718`) | [`baryc`] |
//! | `BoxParameters` (`cxx:720-748`) | [`box_parameters`] |
//! | `GetOffset` (`cxx:750-777`) | [`get_offset`] |
//! | `CreateCyl` (`cxx:779-810`) | [`create_cyl`] |
//!
//! `CreateCyl` builds a `BRepPrim_Cylinder` (`gp_Ax2, R, H`). The port keeps
//! only the identity-frame cylinder (`primitives::BRepPrimCylinder`, the
//! `BRepPrim_OneAxis` branch at `BRepPrim_OneAxis.cxx:388-439`, `:488-503`,
//! `:660-684`), so a frame `a2` reuses that branch and places the result with
//! the `gp_Trsf` from the local system to the absolute one
//! (`gp_Trsf::SetTransformation(FromSystem, gp::XOY())`, `gp_Trsf.cxx:172-194`),
//! which is `TopoDS_Shape::Move` (`TopoDS_Shape.cxx:184-195`) and therefore
//! shares the TShapes. `BRepPrim_Cylinder` itself builds the same three faces
//! in that frame (`BRepPrim_Cylinder.hxx:34-84`).

use std::f64::consts::PI;

use occt_core::bnd::BndBox;
use occt_core::elib::clib::line_parameter_ax1;
use occt_core::gp::{GpAx1, GpAx2, GpAx3, GpPnt, GpTrsf, GpXyz};
use occt_core::precision::ANGULAR;
use occt_geom::cs_lib::{self, DerivativeStatus};

use crate::abs::ShapeType;
use crate::brep_bnd_lib;
use crate::brep_tool::BRepTool;
use crate::primitives::BRepPrimCylinder;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape};
use crate::topo_tools::{explorer_shapes, faces_of, shapes_of};

use crate::locope::pnt_face::PntFace;

/// `Baryc(S, B)` (`cxx:690-718`): the mean of 11 parameter samples of every
/// non-degenerated edge of `S`.
///
/// The reference dereferences the `BRep_Tool::Curve` handle unconditionally;
/// the port skips an edge without a 3d curve and answers the zero vector when
/// no edge contributed a sample (the reference divides by zero).
pub fn baryc(s: &TopoShape) -> GpPnt {
    let mut bar = GpXyz::zero();
    let mut nbp = 0usize;
    for shape in explorer_shapes(s, ShapeType::Edge) {
        let e = Edge(shape);
        if BRepTool::is_degenerated(&e) {
            continue;
        }
        let Some(c) = BRepTool::edge_curve_world(&e) else {
            continue;
        };
        let (first, last) = BRepTool::edge_parameters(&e);
        for i in 1..=11 {
            let prm = ((11.0 - i as f64) * first + (i as f64 - 1.0) * last) / 10.0;
            bar = bar.add(c.value(prm).xyz());
            nbp += 1;
        }
    }
    if nbp == 0 {
        return GpPnt::new(bar.x(), bar.y(), bar.z());
    }
    bar = bar.multiply_scalar(1.0 / nbp as f64);
    GpPnt::new(bar.x(), bar.y(), bar.z())
}

/// `BoxParameters(S, Axis, parmin, parmax)` (`cxx:720-748`): the extreme
/// parameters of the shape's bounding box along the axis.
///
/// `BRepBndLib::Add(S, B)` with the default `useTriangulation` (`cxx:725`).
/// The reference's `B.Get` raises on an empty box; the port answers `None`.
pub fn box_parameters(s: &TopoShape, axis: &GpAx1) -> Option<(f64, f64)> {
    let mut b = BndBox::new();
    brep_bnd_lib::add(s, &mut b, true);
    let (xmin, xmax, ymin, ymax, zmin, zmax) = b.get()?;
    let c = [xmin, xmax, ymin, ymax, zmin, zmax];
    let mut parmin = f64::MAX;
    let mut parmax = -f64::MAX;
    for i in 0..=1 {
        for j in 2..=3 {
            for k in 4..=5 {
                let p = GpPnt::new(c[i], c[j], c[k]);
                let param = line_parameter_ax1(axis, &p);
                parmin = parmin.min(param);
                parmax = parmax.max(param);
            }
        }
    }
    Some((parmin, parmax))
}

/// `GetOffset(PntInfo, Radius, Axis, outOff)` (`cxx:750-777`): the length the
/// cylinder must overshoot the hit point, `Radius * |tan(angle)|` between the
/// axis and the face normal. `None` is the reference's `return false`.
///
/// The reference constructs `BRepAdaptor_Surface` on `PntInfo.Face()`; the
/// port answers `None` when the face carries no surface (the reference
/// dereferences the adaptor's null handle).
pub fn get_offset(pnt_info: &PntFace, radius: f64, axis: &GpAx1) -> Option<f64> {
    let surface = BRepTool::face_surface_world(pnt_info.face())?;
    let (_, d1u, d1v) = surface.d1(pnt_info.u_parameter(), pnt_info.v_parameter());
    let (stat, norm) = cs_lib::normal(&d1u, &d1v, ANGULAR);
    if stat != DerivativeStatus::Done {
        return None;
    }
    let norm = norm?;
    let angle = axis.direction().angle(&norm);
    if (PI / 2.0 - angle).abs() < ANGULAR {
        return None;
    }
    Some(radius * angle.tan().abs())
}

/// `CreateCyl(PntInfoFirst, PntInfoLast, Radius, Axis, Cyl, CylTopF, CylBottF)`
/// (`cxx:779-810`).
pub fn create_cyl(
    first: &PntFace,
    last: &PntFace,
    radius: f64,
    axis: &GpAx1,
) -> Result<(Shell, Face, Face), String> {
    let first_par = first.parameter();
    let last_par = last.parameter();
    let height = last_par - first_par;
    let off_f = get_offset(first, radius, axis).unwrap_or(radius);
    let off_l = get_offset(last, radius, axis).unwrap_or(radius);

    // `cxx:802`: `gp_XYZ theOrig = PntInfoFirst.Pnt().XYZ() - offF * Axis.Direction().XYZ()`.
    let orig = first
        .pnt()
        .xyz()
        .subtract(&axis.direction().xyz().multiply_scalar(off_f));
    let pnt = GpPnt::new(orig.x(), orig.y(), orig.z());
    let a2 = GpAx2::from_axis(pnt, *axis.direction());
    let (solid, top, bottom) = cylinder(&a2, radius, height + off_f + off_l)?;
    let shell = shapes_of(&solid.0, ShapeType::Shell)
        .into_iter()
        .next()
        .ok_or_else(|| "Standard_DomainError: BRepPrim_Cylinder has no shell".to_string())?;
    Ok((Shell(shell), top, bottom))
}

/// `BRepPrim_Cylinder(a2, R, H)`: the cylinder of the frame `a2`, answering the
/// solid plus its top and bottom faces (`BRepPrim_Cylinder::TopFace` /
/// `BottomFace`).
///
/// `R <= 0` and `H <= 0` are `BRepPrim_Cylinder`'s `Standard_DomainError`.
///
/// The frame is **baked into the geometry**, not attached as a `TopLoc_Location`:
/// `BRepPrim_Cylinder(gp_Ax2, R, H)` builds its faces through `BRepPrim_OneAxis`
/// on the given `gp_Ax2` (`BRepPrim_Cylinder.cxx:34-64`), and every edge, vertex
/// and surface of that result is expressed in the absolute frame -- the shape it
/// answers has an identity location. Attaching a location instead (the earlier
/// implementation) hands `BOPAlgo_BOP` a tool whose faces live in a location the
/// boolean's shape set does not know, which fails with
/// `BOPAlgo_AlertUnknownShape` on a through hole and silently mis-builds the
/// `PartsOfTool` split on the `CreateCyl` path.
pub fn cylinder(a2: &GpAx2, radius: f64, height: f64) -> Result<(Solid, Face, Face), String> {
    if radius <= 0.0 || height <= 0.0 {
        return Err("Standard_DomainError: BRepPrim_Cylinder".to_string());
    }
    let prim = BRepPrimCylinder::make_cylinder(radius, height);
    let mut t = GpTrsf::identity();
    // `gp_Trsf::SetTransformation(FromSystem, gp::XOY())` (`gp_Trsf.cxx:172-194`).
    t.set_transformation_from_to(&a2.to_ax3(), &GpAx3::standard());
    let solid = Solid(crate::shape_ops::transformed_copy(&prim.solid.0, &t)?);
    // `primitives::BRepPrimCylinder::make_cylinder` stores bottom, top, lateral
    // in that order (`BRepPrim_OneAxis::BottomFace` / `TopFace` /
    // `LateralFace`), the same convention `brep_builder_api::make_cylinder_full`
    // relies on; `transformed_copy` keeps that order.
    let faces = faces_of(&solid.0);
    let bottom = faces
        .first()
        .cloned()
        .ok_or_else(|| "Standard_DomainError: BRepPrim_Cylinder has no bottom face".to_string())?;
    let top = faces
        .get(1)
        .cloned()
        .ok_or_else(|| "Standard_DomainError: BRepPrim_Cylinder has no top face".to_string())?;
    Ok((solid, top, bottom))
}
