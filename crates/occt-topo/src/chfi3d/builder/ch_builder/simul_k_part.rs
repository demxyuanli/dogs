//! `ChFi3d_ChBuilder::SimulKPart` — the simulation sections of a chamfer
//! `ChFiDS_SurfData`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_ChBuilder.cxx:655-716`
//! (OCCT 8.0.0). The chamfer body only fills the sections of a **plane** and a
//! **cone** support surface; every other type falls through the `default` arm
//! and stores a null section. This is the chamfer twin of
//! `ChFi3d_FilBuilder::SimulKPart` (`ChFi3d_FilBuilder.cxx:475-554`), which
//! covers cylinder / torus / sphere instead — the two `switch`es do not share a
//! single arm except the `default`.
//!
//! Translation notes:
//! - `DStr.Surface(SD->Surf()).Surface()` — a null surface leaves `sec` null in
//!   OCCT (the `switch` is skipped and `SetSimul(null)` is stored); the port
//!   answers `Err` because a null surface has no `GetType()` to switch on,
//! - `sec->ChangeValue(i)` with `1 <= i <= n` is bit-for-bit the port's
//!   `change_value`,
//! - `ElSLib::PlaneUIso` / `ElSLib::ConeUIso` are
//!   [`occt_core::elib::slib::plane_u_iso`] / `cone_u_iso`.

use std::sync::Arc;

use occt_core::elib::slib;
use occt_core::gp::GpPnt2d;

use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::chfi_ds::{CircSectionArray, SurfData};
use crate::topopebrep::ds::HDataStructure;

/// The `Value()` of one of the four `PCurveOnSurf`s of `SD` (`cxx:659-667`).
fn p_curve_point(sd: &SurfData, on_s: i32, first: bool) -> Result<GpPnt2d, String> {
    let interf = sd.interference(on_s);
    let pc = interf.p_curve_on_surf().ok_or_else(|| {
        format!("SimulKPart: PCurveOnSurf of interference {on_s} is null")
    })?;
    let par = if first {
        interf.first_parameter()
    } else {
        interf.last_parameter()
    };
    Ok(pc.value(par))
}

/// `(int)(36. * ang / M_PI + 1)`, clamped to at least two (`cxx:694-698`).
fn section_count(ang: f64) -> i32 {
    let n = (36.0 * ang / std::f64::consts::PI + 1.0) as i32;
    if n < 2 {
        2
    } else {
        n
    }
}

/// `ChFi3d_ChBuilder::SimulKPart(SD)` (`cxx:655-716`).
pub(super) fn simul_k_part_impl(
    sd: &mut SurfData,
    ds: &HDataStructure,
) -> Result<(), String> {
    // `occ::handle<Geom_Surface> S = DStr.Surface(SD->Surf()).Surface();`
    // (`cxx:658`).
    let surface = ds
        .ds()
        .surface(sd.surf())
        .surface()
        .ok_or_else(|| "SimulKPart: DStr.Surface(SD->Surf()).Surface() is null".to_string())?
        .clone();

    // `gp_Pnt2d p1f = SD->InterferenceOnS1().PCurveOnSurf()->Value(...)`
    // (`cxx:659-667`).
    let p1f = p_curve_point(sd, 1, true)?;
    let p1l = p_curve_point(sd, 1, false)?;
    let p2f = p_curve_point(sd, 2, true)?;
    let p2l = p_curve_point(sd, 2, false)?;

    let mut sec: Option<CircSectionArray> = None;
    match classify_surface(surface.as_ref()) {
        // `case GeomAbs_Plane:` (`cxx:672-686`).
        SurfaceKind::Plane => {
            let v1 = p1f.y();
            let v2 = p2f.y();
            let u1 = p1f.x().max(p2f.x());
            let u2 = p1l.x().min(p2l.x());
            let pl = surface
                .gp_pln()
                .ok_or_else(|| "SimulKPart: GeomAdaptor_Surface::Plane() on a plane".to_string())?;
            let mut arr = CircSectionArray::new(1, 2);
            arr.change_value(1)
                .ok_or_else(|| "SimulKPart: CircSection(1)".to_string())?
                .set_lin(&slib::plane_u_iso(&pl.position(), u1), v1, v2);
            arr.change_value(2)
                .ok_or_else(|| "SimulKPart: CircSection(2)".to_string())?
                .set_lin(&slib::plane_u_iso(&pl.position(), u2), v1, v2);
            sec = Some(arr);
        }
        // `case GeomAbs_Cone:` (`cxx:687-707`).
        SurfaceKind::Cone => {
            let v1 = p1f.y();
            let v2 = p2f.y();
            let u1 = p1f.x().max(p2f.x());
            let u2 = p1l.x().min(p2l.x());
            let ang = u2 - u1;
            let co = surface
                .gp_cone()
                .ok_or_else(|| "SimulKPart: GeomAdaptor_Surface::Cone() on a cone".to_string())?;
            let rad = co.radius;
            let sang = co.semi_angle;
            let n = section_count(ang);
            let mut arr = CircSectionArray::new(1, n);
            for i in 1..=n {
                let u = u1 + f64::from(i - 1) * (u2 - u1) / f64::from(n - 1);
                arr.change_value(i)
                    .ok_or_else(|| format!("SimulKPart: CircSection({i})"))?
                    .set_lin(&slib::cone_u_iso(&co.pos, rad, sang, u), v1, v2);
            }
            sec = Some(arr);
        }
        // `default: break;` (`cxx:708-709`): `sec` stays null.
        _ => {}
    }

    // `SD->SetSimul(sec);` (`cxx:711`): `sec` is null on the `default` arm, and
    // `SetSimul` stores the handle unconditionally.
    match sec {
        Some(arr) => sd.set_simul(Arc::new(arr)),
        None => sd.reset_simul(),
    }
    Ok(())
}
