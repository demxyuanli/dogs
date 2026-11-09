//! Entity 112. `IGESToBRep_BasicCurve::TransferSplineCurve`
//! (`IGESToBRep_BasicCurve.cxx:878-930`) and `Transfer2dSplineCurve`
//! (`cxx:934-971`). `IGESConvGeom::SplineCurveFromIGES`
//! (`IGESConvGeom.cxx:37-169`).
//!
//! The directory 124 block is commented out (`IGESConvGeom.cxx:153-164`), so
//! poles stay in file coordinates. `TransferBasicCurve` then scales the 3D
//! curve (`cxx:231`). The 2D path calls `TransferSplineCurve` with
//! `EpsGeom = PConfusion` and drops Z; it does not scale.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::curve::Curve;
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::curve::Curve2d;

use super::copious_data::increase_continuity_3d;
use super::iges_data_iges_model::{EntityData, IgesModel};
use super::iges_data_param_reader::ParamReader;
use super::iges_to_brep_basic_curve::{eps_geom, scale_curve, Curve2dSpan, CurveSpan};

/// `read.iges.bspline.continuity` default (`IGESData.cxx:66`).
const CONTINUITY: i32 = 1;

pub(in crate::iges_read) fn parse_spline_curve(pr: &mut ParamReader<'_>) -> Result<EntityData, String> {
    // `IGESGeom_ToolSplineCurve.cxx:61-180`.
    let spline_type = pr.read_integer()?;
    let continuity_degree = pr.read_integer()?;
    let nb_dimensions = pr.read_integer()?;
    let nb_segments = pr.read_integer()?;
    if nb_segments <= 0 {
        return Err("XSTEP_94: IGESGeom_ToolSplineCurve.cxx:86-90".into());
    }
    let n = nb_segments as usize;
    let mut breaks = Vec::with_capacity(n + 1);
    for _ in 0..=n {
        breaks.push(pr.read_real()?);
    }
    let mut x_poly = Vec::with_capacity(n);
    let mut y_poly = Vec::with_capacity(n);
    let mut z_poly = Vec::with_capacity(n);
    for _ in 0..n {
        x_poly.push(read4(pr)?);
        y_poly.push(read4(pr)?);
        z_poly.push(read4(pr)?);
    }
    let x_values = read4(pr)?;
    let y_values = read4(pr)?;
    let z_values = read4(pr)?;
    Ok(EntityData::SplineCurve {
        spline_type,
        continuity_degree,
        nb_dimensions,
        breaks,
        x_poly,
        y_poly,
        z_poly,
        x_values,
        y_values,
        z_values,
    })
}

pub(in crate::iges_read) fn transfer_spline_curve(
    model: &IgesModel,
    index: usize,
) -> Result<CurveSpan, String> {
    // `TransferSplineCurve` (`cxx:878-929`) then unit scale (`cxx:231`).
    let curve = spline_from_iges(model, index, eps_geom(model.resolution))?;
    let first = curve.first_parameter();
    let last = curve.last_parameter();
    let scaled = scale_curve(Arc::new(curve), model.unit_factor)?;
    Ok(CurveSpan {
        curve: scaled,
        first,
        last,
    })
}

pub(in crate::iges_read) fn transfer_2d_spline_curve(
    model: &IgesModel,
    index: usize,
) -> Result<Curve2dSpan, String> {
    // `Transfer2dSplineCurve` (`cxx:934-971`): 3D transfer at `PConfusion`, then XY.
    let curve = spline_from_iges(model, index, PCONFUSION)?;
    let degree = curve.degree();
    let (knots, mults) = curve.distinct_knots_and_mults();
    let flat = occt_core::bspl::banded_interp::knot_sequence(&knots, &mults, degree as i32);
    let xs: Vec<f64> = curve.poles.iter().map(|p| p.x()).collect();
    let ys: Vec<f64> = curve.poles.iter().map(|p| p.y()).collect();
    let bs = Geom2dBSplineCurve::new(xs, ys, flat, degree)
        .map_err(|e| format!("IGES 112 2d curve: {e} (IGESToBRep_BasicCurve.cxx:970)"))?;
    let first = bs.first_parameter();
    let last = bs.last_parameter();
    Ok(Curve2dSpan {
        curve: Arc::new(bs),
        first,
        last,
    })
}

fn spline_from_iges(
    model: &IgesModel,
    index: usize,
    epsgeom: f64,
) -> Result<GeomBSplineCurve, String> {
    let result = spline_curve_from_iges(model, index, epsgeom)?;
    let mut curve = result.curve;
    // `IncreaseCurveContinuity` (`BasicCurve.cxx:926-928`).
    let tol = CONFUSION.min(epsgeom);
    increase_continuity_3d(&mut curve, tol, CONTINUITY)?;
    Ok(curve)
}

struct Built {
    curve: GeomBSplineCurve,
}

fn spline_curve_from_iges(
    model: &IgesModel,
    index: usize,
    epsgeom: f64,
) -> Result<Built, String> {
    // `IGESConvGeom::SplineCurveFromIGES` (`cxx:37-169`).
    let (spline_type, nb_dimensions, breaks, x_poly, y_poly, z_poly) =
        match &model.entity(index)?.data {
            EntityData::SplineCurve {
                spline_type,
                nb_dimensions,
                breaks,
                x_poly,
                y_poly,
                z_poly,
                ..
            } => (
                *spline_type,
                *nb_dimensions,
                breaks.clone(),
                x_poly.clone(),
                y_poly.clone(),
                z_poly.clone(),
            ),
            _ => return Err("IGES 112 has no spline-curve parameters".into()),
        };
    let mut degree = spline_type;
    if degree > 3 {
        degree = 3;
    }
    if degree < 0 {
        return Err("Standard_RangeError: IGESConvGeom.cxx:99".into());
    }
    let nb_segs = breaks.len().saturating_sub(1);
    if nb_segs < 1 || x_poly.len() < nb_segs {
        return Err("XSTEP_246: IGESToBRep_BasicCurve.cxx:896-900".into());
    }
    let degree_us = degree as usize;
    let nb_knots = nb_segs + 1;
    let mut mults = vec![degree; nb_knots];
    if let Some(first) = mults.first_mut() {
        *first = degree + 1;
    }
    if let Some(last) = mults.last_mut() {
        *last = degree + 1;
    }
    let knots = breaks[..nb_knots].to_vec();
    let npoles = nb_segs * degree_us + 1;
    let mut bspoles = vec![GpPnt::new(0.0, 0.0, 0.0); npoles];
    // OCCT `ibspole` starts at `Lower() - 1` (`cxx:81`), 1-based.
    let mut ibspole = 0i32;
    for i in 1..=nb_segs {
        let di = breaks[i] - breaks[i - 1];
        let di2 = di * di;
        let di3 = di2 * di;
        let x = x_poly[i - 1];
        let y = y_poly[i - 1];
        let z = z_poly[i - 1];
        let (bz, cz, dz) = if nb_dimensions == 2 {
            (0.0, 0.0, 0.0)
        } else {
            (z[1], z[2], z[3])
        };
        let mut coeff = vec![GpPnt::new(0.0, 0.0, 0.0); degree_us + 1];
        if degree_us >= 1 {
            if degree_us >= 3 {
                coeff[3] = GpPnt::new(x[3] * di3, y[3] * di3, dz * di3);
            }
            if degree_us >= 2 {
                coeff[2] = GpPnt::new(x[2] * di2, y[2] * di2, cz * di2);
            }
            coeff[1] = GpPnt::new(x[1] * di, y[1] * di, bz * di);
            coeff[0] = GpPnt::new(x[0], y[0], z[0]);
        }
        let bzpoles = occt_core::bspl::plib::coefficients_poles(&coeff);
        if ibspole > 1 {
            let prev = bspoles[(ibspole as usize) - 1];
            let first = bzpoles[0];
            if prev.distance(&first) > epsgeom {
                bspoles[(ibspole as usize) - 1] = GpPnt::new(
                    (prev.x() + first.x()) * 0.5,
                    (prev.y() + first.y()) * 0.5,
                    (prev.z() + first.z()) * 0.5,
                );
            }
        }
        if i == 1 {
            ibspole += 1;
            if (ibspole as usize) > bspoles.len() || bzpoles.is_empty() {
                return Err("IGES_1175: IGESToBRep_BasicCurve.cxx:908-912".into());
            }
            bspoles[(ibspole as usize) - 1] = bzpoles[0];
        }
        for pole in bzpoles.iter().skip(1) {
            ibspole += 1;
            if (ibspole as usize) > bspoles.len() {
                return Err("IGES_1175: IGESToBRep_BasicCurve.cxx:908-912".into());
            }
            bspoles[(ibspole as usize) - 1] = *pole;
        }
    }
    if ibspole != npoles as i32 {
        return Err("IGES_1175: IGESToBRep_BasicCurve.cxx:908-912".into());
    }
    let curve = GeomBSplineCurve::from_poles_knots_mults(bspoles, knots, mults, degree_us)
        .map_err(|e| format!("IGES 112 curve: {e} (IGESToBRep_BasicCurve.cxx:216-223)"))?;
    Ok(Built { curve })
}

fn read4(pr: &mut ParamReader<'_>) -> Result<[f64; 4], String> {
    Ok([
        pr.read_real()?,
        pr.read_real()?,
        pr.read_real()?,
        pr.read_real()?,
    ])
}
