//! `BRepFill_MultiLine::Curves` — the analytic 3D curve and the two p-curves
//! of a particular case.
//!
//! Source: `ModelingAlgorithms/TKBool/BRepFill/BRepFill_MultiLine.cxx:405-548`.
//!
//! `myKPart == 1` (`cxx:406-508`) means the bisector runs along an iso line of
//! both faces: the 3D curve is that iso trimmed between the two v (or u) values
//! of the bisector ends, and each p-curve is a line at the constant other
//! parameter. `myKPart == 2` (`cxx:509-547`) means the bisector is a straight
//! line and both faces are planes: the 3D curve is the bisector projected on
//! the first plane along `gp::DZ()` and each p-curve the projection of that
//! curve on its plane. Anything else is not a particular case and OCCT raises
//! (`hxx:72`).
//!
//! The trims of `cxx:433` / `cxx:439` / `cxx:457` / `cxx:463` use
//! `Geom_TrimmedCurve(C, U1, U2, Sense)` with the default
//! `theAdjustPeriodic = true`; [`GeomTrimmedCurveBasis::with_sense`] is that
//! constructor, and its `FirstParameter()` / `LastParameter()` are the ones the
//! p-curves of `cxx:442-449` / `cxx:466-473` are built from.
//!
//! UNPORTED: the `myKPart == 2` branch calls `GeomProjLib::ProjectOnPlane`
//! (`GeomProjLib.cxx:293-311` -> `ProjLib_ProjectOnPlane.cxx`), which has no
//! port; [`curves`] reports [`Error::UnportedProjectOnPlane`] instead.

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d};
use occt_geom::curve::Curve;
use occt_geom::trimmed::GeomTrimmedCurveBasis;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::brep_tool::BRepTool;

use super::{Error, MultiLine};

/// `gp::DX2d()` / `-gp::DX2d()`.
fn dx2d(sign: f64) -> GpDir2d {
    GpDir2d::new(sign, 0.0).expect("gp::DX2d")
}

/// `gp::DY2d()` / `-gp::DY2d()`.
fn dy2d(sign: f64) -> GpDir2d {
    GpDir2d::new(0.0, sign).expect("gp::DY2d")
}

/// `Curves(Curve, PCurve1, PCurve2)` (`cxx:405-548`).
pub(super) fn curves(
    m: &MultiLine,
) -> Result<(Arc<dyn Curve>, Arc<dyn Curve2d>, Arc<dyn Curve2d>), Error> {
    match m.k_part {
        1 => k_part_1(m),
        2 => Err(Error::UnportedProjectOnPlane),
        _ => Err(Error::NotParticularCase),
    }
}

/// `myKPart == 1` (`cxx:406-508`): the bisector is an iso of the two faces, so
/// the 3D curve and the p-curves are analytic.
fn k_part_1(
    m: &MultiLine,
) -> Result<(Arc<dyn Curve>, Arc<dyn Curve2d>, Arc<dyn Curve2d>), Error> {
    let (first, last) = (m.first_parameter(), m.last_parameter());

    // `P1 = ValueOnF1(myBis.FirstParameter()); P2 = ValueOnF1(...Last...);
    // PMil = ValueOnF1(0.5 * (First + Last));` (`cxx:410-418`).
    let mut p1 = m.value_on_f1(first)?;
    let mut p2 = m.value_on_f1(last)?;
    let mut p_mil = m.value_on_f1(0.5 * (first + last))?;

    // `S = BRep_Tool::Surface(myFace1, L)` transformed by `L`
    // (`cxx:420-425`).
    let s = BRepTool::face_surface_world(&m.face1).ok_or(Error::NullSurface)?;

    let sense: bool;
    let curve: Arc<dyn Curve>;
    let f: f64;
    let l: f64;
    let pcurve1: Arc<dyn Curve2d>;

    if !m.iso_u1 {
        // `cxx:428-450`: a `U = const` face-1 iso, trimmed along V.
        let iso = s.u_iso_curve(p_mil.x()).ok_or(Error::NullCurve)?;
        sense = p1.y() < p2.y();
        let trimmed = if sense {
            GeomTrimmedCurveBasis::with_sense(iso, p1.y(), p2.y(), sense)
        } else {
            GeomTrimmedCurveBasis::with_sense(iso, p2.y(), p1.y(), sense)
        }
        .map_err(|_| Error::TrimError)?;
        f = trimmed.first_parameter();
        l = trimmed.last_parameter();
        let line = if sense {
            Geom2dLine::from_pnt_dir(GpPnt2d::new(p_mil.x(), p1.y() - f), dy2d(1.0))
        } else {
            Geom2dLine::from_pnt_dir(GpPnt2d::new(p_mil.x(), p1.y() + f), dy2d(-1.0))
        };
        pcurve1 = Arc::new(Geom2dTrimmedCurve::new(Arc::new(line), f, l));
        curve = Arc::new(trimmed);
    } else {
        // `cxx:451-474`: a `V = const` face-1 iso, trimmed along U.
        let iso = s.v_iso_curve(p_mil.y()).ok_or(Error::NullCurve)?;
        sense = p1.x() < p2.x();
        let trimmed = if sense {
            GeomTrimmedCurveBasis::with_sense(iso, p1.x(), p2.x(), sense)
        } else {
            GeomTrimmedCurveBasis::with_sense(iso, p2.x(), p1.x(), sense)
        }
        .map_err(|_| Error::TrimError)?;
        f = trimmed.first_parameter();
        l = trimmed.last_parameter();
        let line = if sense {
            Geom2dLine::from_pnt_dir(GpPnt2d::new(p1.x() - f, p_mil.y()), dx2d(1.0))
        } else {
            Geom2dLine::from_pnt_dir(GpPnt2d::new(p1.x() + f, p_mil.y()), dx2d(-1.0))
        };
        pcurve1 = Arc::new(Geom2dTrimmedCurve::new(Arc::new(line), f, l));
        curve = Arc::new(trimmed);
    }

    // `P1 = ValueOnF2(...); P2 = ValueOnF2(...); PMil = ValueOnF2(...);`
    // (`cxx:477-480`).
    p1 = m.value_on_f2(first)?;
    p2 = m.value_on_f2(last)?;
    p_mil = m.value_on_f2(0.5 * (first + last))?;

    // `cxx:482-507`: the second p-curve, still trimmed with the *first* face's
    // `f` / `l` — the 3D curve is the same one.
    let pcurve2: Arc<dyn Curve2d> = if !m.iso_u2 {
        let sense = p1.y() < p2.y();
        let line = if sense {
            Geom2dLine::from_pnt_dir(GpPnt2d::new(p_mil.x(), p1.y() - f), dy2d(1.0))
        } else {
            Geom2dLine::from_pnt_dir(GpPnt2d::new(p_mil.x(), p1.y() + f), dy2d(-1.0))
        };
        Arc::new(Geom2dTrimmedCurve::new(Arc::new(line), f, l))
    } else {
        let sense = p1.x() < p2.x();
        let line = if sense {
            Geom2dLine::from_pnt_dir(GpPnt2d::new(p1.x() - f, p_mil.y()), dx2d(1.0))
        } else {
            Geom2dLine::from_pnt_dir(GpPnt2d::new(p1.x() + f, p_mil.y()), dx2d(-1.0))
        };
        Arc::new(Geom2dTrimmedCurve::new(Arc::new(line), f, l))
    };

    Ok((curve, pcurve1, pcurve2))
}
