//! `IGESGeom_ConicArc` (`IGESGeom_ConicArc.cxx:121-364`) and
//! `IGESToBRep_BasicCurve::TransferConicArc` / `Transfer2dConicArc`
//! (`IGESToBRep_BasicCurve.cxx:307-713`).
//!
//! `TransferTopoBasicCurve` sets `ModeTransfer` false (`cxx:967-968`), so a
//! directory 124 is applied. `Transfer2dConicArc` skips that 124 when the
//! transformed axis is not parallel to Z within `1e-3` (`cxx:526-533`).

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::elib::{clib, clib2d};
use occt_core::gp::{
    GpAx2, GpAx22d, GpCirc, GpCirc2d, GpDir, GpDir2d, GpElips, GpElips2d, GpHypr, GpHypr2d, GpParab,
    GpParab2d, GpPnt, GpPnt2d, GpXyz,
};
use occt_core::precision::{CONFUSION, PCONFUSION, REAL_SMALL};
use occt_geom::curve::Curve;
use occt_geom::{GeomCircle, GeomEllipse, GeomHyperbola, GeomParabola};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_geom2d::{Geom2dCircle, Geom2dEllipse, Geom2dHyperbola, Geom2dParabola};

use super::iges_data_iges_model::{EntityData, IgesModel};
use super::iges_to_brep_basic_curve::{
    apply_gtrsf, entity_location, scale_curve, transformed_dir, Curve2dSpan, CurveSpan,
};

const FORM_EPS: f64 = 1.0e-8;
const FORM_EPS2: f64 = FORM_EPS * FORM_EPS;
const FORM_EPS4: f64 = FORM_EPS2 * FORM_EPS2;

struct Coeff {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    f: f64,
    zt: f64,
    start: [f64; 2],
    end: [f64; 2],
    form: i32,
}

pub(in crate::iges_read) fn transfer_conic_arc(
    model: &IgesModel,
    index: usize,
) -> Result<CurveSpan, String> {
    let coeff = load(model, index)?;
    if computed_form(&coeff) == 0 {
        return Err("IGES_1155: IGESToBRep_BasicCurve.cxx:318-323".into());
    }
    let placed = place_3d(model, index, &coeff)?;
    let frame = GpAx2::new(placed.center, placed.norm, placed.main)
        .map_err(|_| "IGES 104 frame is singular: gp_Ax2".to_string())?;
    let closed = coeff.start[0] == coeff.end[0] && coeff.start[1] == coeff.end[1];
    let (curve, first, last, scale_u) = if is_ellipse(&coeff)
        && (coeff.a - coeff.c).abs() <= PCONFUSION * PCONFUSION
        && coeff.b.abs() < PCONFUSION * PCONFUSION
    {
        let curve: Arc<dyn Curve> = Arc::new(GeomCircle::new(GpCirc::new(frame, placed.rmin)));
        let t1 = clib::circle_parameter(&frame, &placed.start);
        let mut t2 = clib::circle_parameter(&frame, &placed.end);
        if !closed {
            if t1 > t2 && (t1 - t2) > CONFUSION {
                t2 += 2.0 * PI;
            }
            if (t1 - t2).abs() <= CONFUSION {
                (curve, 0.0, 2.0 * PI, false)
            } else {
                (curve, t1, t2, false)
            }
        } else {
            (curve, 0.0, 2.0 * PI, false)
        }
    } else if is_parabola(&coeff) {
        let focal = placed.rmin / 2.0;
        let curve: Arc<dyn Curve> = Arc::new(GeomParabola::new(GpParab::new(frame, focal)));
        let parab = GpParab::new(frame, focal);
        let t1 = clib::parameter_parab(&parab, &placed.start);
        let t2 = clib::parameter_parab(&parab, &placed.end);
        if (t1 - t2).abs() <= CONFUSION {
            (curve, f64::NEG_INFINITY, f64::INFINITY, true)
        } else {
            (curve, t1, t2, true)
        }
    } else if is_ellipse(&coeff) {
        let curve: Arc<dyn Curve> =
            Arc::new(GeomEllipse::new(GpElips::new(frame, placed.rmax, placed.rmin)));
        let elips = GpElips::new(frame, placed.rmax, placed.rmin);
        if closed {
            (curve, 0.0, 2.0 * PI, false)
        } else {
            let t1 = clib::parameter_elips(&elips, &placed.start);
            let mut t2 = clib::parameter_elips(&elips, &placed.end);
            if t2 < t1 && (t1 - t2) > CONFUSION {
                t2 += 2.0 * PI;
            }
            if (t1 - t2).abs() <= CONFUSION {
                (curve, 0.0, 2.0 * PI, false)
            } else {
                (curve, t1, t2, false)
            }
        }
    } else {
        let curve: Arc<dyn Curve> =
            Arc::new(GeomHyperbola::new(GpHypr::new(frame, placed.rmax, placed.rmin)));
        let hypr = GpHypr::new(frame, placed.rmax, placed.rmin);
        let t1 = clib::parameter_hypr(&hypr, &placed.start);
        let t2 = clib::parameter_hypr(&hypr, &placed.end);
        if (t1 - t2).abs() <= PCONFUSION {
            (curve, f64::NEG_INFINITY, f64::INFINITY, false)
        } else if t1 > t2 {
            (curve, t2, t1, false)
        } else {
            (curve, t1, t2, false)
        }
    };
    let curve = scale_curve(curve, model.unit_factor)?;
    let factor = if scale_u { model.unit_factor.abs() } else { 1.0 };
    Ok(CurveSpan {
        curve,
        first: scale_param(first, factor),
        last: scale_param(last, factor),
    })
}

pub(in crate::iges_read) fn transfer_2d_conic_arc(
    model: &IgesModel,
    index: usize,
) -> Result<Curve2dSpan, String> {
    let coeff = load(model, index)?;
    if computed_form(&coeff) == 0 {
        return Err("IGES_1155: IGESToBRep_BasicCurve.cxx:510-514".into());
    }
    let axis_z = GpDir::from_axis(occt_core::gp::dir::DirAxis::Z);
    let loc = entity_location(model, index)?;
    let mut skip = false;
    if let Some(g) = loc.as_ref() {
        let transformed = transformed_dir(g, GpXyz::new(0.0, 0.0, 1.0))?;
        if !transformed.is_parallel_tol(&axis_z, 1.0e-3) {
            skip = true;
        }
    }
    let placed = if skip {
        place_local(&coeff)?
    } else {
        place_3d(model, index, &coeff)?
    };
    let mut frame = ax22d(placed.center, placed.main)?;
    if !skip {
        if let Some(g) = loc.as_ref() {
            let transformed = transformed_dir(g, GpXyz::new(0.0, 0.0, 1.0))?;
            if transformed.is_parallel_tol(&axis_z, 1.0e-3) && transformed.dot(&axis_z) < 0.0 {
                frame = reverse_ax22d(frame)?;
            }
        }
    }
    let start = GpPnt2d::new(placed.start.x(), placed.start.y());
    let end = GpPnt2d::new(placed.end.x(), placed.end.y());
    let closed = coeff.start[0] == coeff.end[0] && coeff.start[1] == coeff.end[1];
    if is_ellipse(&coeff)
        && (coeff.a - coeff.c).abs() <= PCONFUSION * PCONFUSION
        && coeff.b.abs() < PCONFUSION * PCONFUSION
    {
        let curve: Arc<dyn Curve2d> = Arc::new(Geom2dCircle::new(GpCirc2d::new(frame, placed.rmin)));
        if closed {
            return Ok(Curve2dSpan {
                curve,
                first: 0.0,
                last: 2.0 * PI,
            });
        }
        let circ = GpCirc2d::new(frame, placed.rmin);
        let t1 = clib2d::parameter_circ2d(&circ, &start);
        let mut t2 = clib2d::parameter_circ2d(&circ, &end);
        if t2 < t1 && (t1 - t2) > PCONFUSION {
            t2 += 2.0 * PI;
        }
        if (t1 - t2).abs() <= PCONFUSION {
            return Ok(Curve2dSpan {
                curve,
                first: 0.0,
                last: 2.0 * PI,
            });
        }
        return Ok(trimmed_2d(curve, t1, t2));
    }
    if is_parabola(&coeff) {
        let focal = placed.rmin / 2.0;
        let curve: Arc<dyn Curve2d> = Arc::new(Geom2dParabola::new(GpParab2d::new(frame, focal)));
        let parab = GpParab2d::new(frame, focal);
        let t1 = clib2d::parameter_parab2d(&parab, &start);
        let t2 = clib2d::parameter_parab2d(&parab, &end);
        if (t1 - t2).abs() <= PCONFUSION {
            return Ok(Curve2dSpan {
                curve,
                first: f64::NEG_INFINITY,
                last: f64::INFINITY,
            });
        }
        if t1 > t2 {
            return Ok(trimmed_2d(curve, t2, t1));
        }
        return Ok(trimmed_2d(curve, t1, t2));
    }
    if is_ellipse(&coeff) {
        let curve: Arc<dyn Curve2d> =
            Arc::new(Geom2dEllipse::new(GpElips2d::new(frame, placed.rmax, placed.rmin)));
        if closed {
            return Ok(Curve2dSpan {
                curve,
                first: 0.0,
                last: 2.0 * PI,
            });
        }
        let elips = GpElips2d::new(frame, placed.rmax, placed.rmin);
        let t1 = clib2d::parameter_elips2d(&elips, &start);
        let mut t2 = clib2d::parameter_elips2d(&elips, &end);
        if t2 < t1 && (t1 - t2) > PCONFUSION {
            t2 += 2.0 * PI;
        }
        if (t1 - t2).abs() <= PCONFUSION {
            return Ok(Curve2dSpan {
                curve,
                first: 0.0,
                last: 2.0 * PI,
            });
        }
        return Ok(trimmed_2d(curve, t1, t2));
    }
    let curve: Arc<dyn Curve2d> =
        Arc::new(Geom2dHyperbola::new(GpHypr2d::new(frame, placed.rmax, placed.rmin)));
    let hypr = GpHypr2d::new(frame, placed.rmax, placed.rmin);
    let t1 = clib2d::parameter_hypr2d(&hypr, &start);
    let t2 = clib2d::parameter_hypr2d(&hypr, &end);
    if (t1 - t2).abs() <= PCONFUSION {
        return Ok(Curve2dSpan {
            curve,
            first: f64::NEG_INFINITY,
            last: f64::INFINITY,
        });
    }
    if t1 > t2 {
        return Ok(trimmed_2d(curve, t2, t1));
    }
    Ok(trimmed_2d(curve, t1, t2))
}

fn trimmed_2d(curve: Arc<dyn Curve2d>, first: f64, last: f64) -> Curve2dSpan {
    Curve2dSpan {
        curve: Arc::new(Geom2dTrimmedCurve::new(curve, first, last)),
        first,
        last,
    }
}

fn scale_param(u: f64, factor: f64) -> f64 {
    if !u.is_finite() {
        u
    } else {
        u * factor
    }
}

struct Place {
    center: GpPnt,
    main: GpDir,
    norm: GpDir,
    rmin: f64,
    rmax: f64,
    start: GpPnt,
    end: GpPnt,
}

fn place_3d(model: &IgesModel, index: usize, coeff: &Coeff) -> Result<Place, String> {
    let def = computed_definition(coeff)?;
    let loc = entity_location(model, index)?;
    if let Some(g) = loc.as_ref() {
        let center = apply_gtrsf(g, GpPnt::new(def.xcen, def.ycen, coeff.zt));
        let main = transformed_dir(g, GpXyz::new(def.xax, def.yax, 0.0))?;
        let norm = transformed_dir(g, GpXyz::new(0.0, 0.0, 1.0))?;
        let start = apply_gtrsf(g, GpPnt::new(coeff.start[0], coeff.start[1], coeff.zt));
        let end = apply_gtrsf(g, GpPnt::new(coeff.end[0], coeff.end[1], coeff.zt));
        Ok(Place {
            center,
            main,
            norm,
            rmin: def.rmin,
            rmax: def.rmax,
            start,
            end,
        })
    } else {
        Ok(Place {
            center: GpPnt::new(def.xcen, def.ycen, coeff.zt),
            main: dir2(def.xax, def.yax)?,
            norm: GpDir::from_axis(occt_core::gp::dir::DirAxis::Z),
            rmin: def.rmin,
            rmax: def.rmax,
            start: GpPnt::new(coeff.start[0], coeff.start[1], coeff.zt),
            end: GpPnt::new(coeff.end[0], coeff.end[1], coeff.zt),
        })
    }
}

fn place_local(coeff: &Coeff) -> Result<Place, String> {
    let def = computed_definition(coeff)?;
    Ok(Place {
        center: GpPnt::new(def.xcen, def.ycen, coeff.zt),
        main: dir2(def.xax, def.yax)?,
        norm: GpDir::from_axis(occt_core::gp::dir::DirAxis::Z),
        rmin: def.rmin,
        rmax: def.rmax,
        start: GpPnt::new(coeff.start[0], coeff.start[1], coeff.zt),
        end: GpPnt::new(coeff.end[0], coeff.end[1], coeff.zt),
    })
}

fn ax22d(center: GpPnt, main: GpDir) -> Result<GpAx22d, String> {
    let dir = GpDir2d::new(main.x(), main.y()).map_err(|_| "IGES 104 2d axis is null".to_string())?;
    Ok(GpAx22d::from_xdir(GpPnt2d::new(center.x(), center.y()), dir))
}

fn reverse_ax22d(ax: GpAx22d) -> Result<GpAx22d, String> {
    let vy = ax.y_direction().reversed();
    GpAx22d::new(*ax.location(), *ax.x_direction(), vy)
        .map_err(|_| "IGES 104 reversed 2d frame is singular".to_string())
}

fn dir2(x: f64, y: f64) -> Result<GpDir, String> {
    GpDir::from_xyz(&GpXyz::new(x, y, 0.0)).map_err(|_| "IGES 104 main axis is null".to_string())
}

fn load(model: &IgesModel, index: usize) -> Result<Coeff, String> {
    let ent = model.entity(index)?;
    match &ent.data {
        EntityData::ConicArc {
            a,
            b,
            c,
            d,
            e,
            f,
            zt,
            start,
            end,
        } => Ok(Coeff {
            a: *a,
            b: *b,
            c: *c,
            d: *d,
            e: *e,
            f: *f,
            zt: *zt,
            start: *start,
            end: *end,
            form: ent.form_number,
        }),
        _ => Err("IGES 104 has no conic parameters".into()),
    }
}

/// `IGESGeom_ConicArc::ComputedFormNumber` (`cxx:121-149`).
fn computed_form(k: &Coeff) -> i32 {
    let q1 = k.a * (k.c * k.f - k.e * k.e / 4.0)
        + k.b / 2.0 * (k.e * k.d / 4.0 - k.b * k.f / 2.0)
        + k.d / 2.0 * (k.b * k.e / 4.0 - k.c * k.d / 2.0);
    let q2 = k.a * k.c - k.b * k.b / 4.0;
    let q3 = k.a + k.c;
    if q2 > FORM_EPS4 && q1 * q3 < 0.0 {
        return 1;
    }
    if q2 < -FORM_EPS4 && q1.abs() > FORM_EPS4 {
        return 2;
    }
    if q2.abs() <= FORM_EPS4 && q1.abs() > FORM_EPS4 {
        return 3;
    }
    0
}

fn branch_form(k: &Coeff) -> i32 {
    if k.form == 0 {
        computed_form(k)
    } else {
        k.form
    }
}

fn is_ellipse(k: &Coeff) -> bool {
    branch_form(k) == 1
}

fn is_parabola(k: &Coeff) -> bool {
    branch_form(k) == 3
}

struct Def {
    xcen: f64,
    ycen: f64,
    xax: f64,
    yax: f64,
    rmin: f64,
    rmax: f64,
}

/// `IGESGeom_ConicArc::ComputedDefinition` (`cxx:241-364`).
fn computed_definition(k: &Coeff) -> Result<Def, String> {
    let a = k.a;
    let b = k.b / 2.0;
    let c = k.c;
    let d = k.d / 2.0;
    let e = k.e / 2.0;
    let f = k.f;
    let def = if is_parabola(k) {
        if a.abs() <= FORM_EPS && b.abs() <= FORM_EPS {
            let focal = -d / c;
            Def {
                xcen: (f * c - e * e) / c / d / 2.0,
                ycen: e / c,
                xax: if focal >= 0.0 { 1.0 } else { -1.0 },
                yax: 0.0,
                rmin: focal.abs(),
                rmax: focal.abs(),
            }
        } else {
            let ss = a + c;
            let cc = -(a * d + b * e) / ss;
            let dd = d + (c * d - b * e) / ss;
            let fc = (a * e - b * d) / ss;
            let ee = e + fc;
            let dn = a * ee - dd * b;
            let mut teta = PI / 2.0;
            if b.abs() > FORM_EPS {
                teta = (-a / b).atan();
            }
            if fc < 0.0 {
                teta += PI;
            }
            Def {
                xcen: (cc * ee + f * b) / dn,
                ycen: (-cc * dd - f * a) / dn,
                xax: teta.cos(),
                yax: teta.sin(),
                rmin: fc.abs() / (a * a + b * b).sqrt() / 2.0,
                rmax: fc.abs() / (a * a + b * b).sqrt() / 2.0,
            }
        }
    } else {
        let gdet = a * c * f + 2.0 * b * d * e - c * d * d - a * e * e - b * b * f;
        let pdet = a * c - b * b;
        let xcen = (b * e - c * d) / pdet;
        let ycen = (b * d - a * e) / pdet;
        let term1 = a - c;
        let term2 = 2.0 * b;
        let (cos2t, auxil) = if term1.abs() < REAL_SMALL {
            (1.0, term2)
        } else {
            let t2d = term2 / term1;
            (1.0 / (1.0 + t2d * t2d).sqrt(), (term1 * term1 + term2 * term2).sqrt())
        };
        let cost = ((1.0 + cos2t) / 2.0).sqrt();
        let sint = ((1.0 - cos2t) / 2.0).sqrt();
        let aprim = (a + c + auxil) / 2.0;
        let cprim = (a + c - auxil) / 2.0;
        let t1 = -gdet / (aprim * pdet);
        let t2 = -gdet / (cprim * pdet);
        if is_ellipse(k) {
            let (rmin, rmax) = if t2.sqrt() < t1.sqrt() {
                (t2.sqrt(), t1.sqrt())
            } else {
                (t1.sqrt(), t2.sqrt())
            };
            Def {
                xcen,
                ycen,
                xax: cost,
                yax: sint,
                rmin,
                rmax,
            }
        } else if t1 <= FORM_EPS {
            Def {
                xcen,
                ycen,
                xax: -sint,
                yax: cost,
                rmin: (-t1).sqrt(),
                rmax: t2.sqrt(),
            }
        } else {
            Def {
                xcen,
                ycen,
                xax: cost,
                yax: sint,
                rmin: (-t2).sqrt(),
                rmax: t1.sqrt(),
            }
        }
    };
    if !def.xcen.is_finite() || !def.ycen.is_finite() || !def.rmin.is_finite() || !def.rmax.is_finite()
    {
        return Err("IGES 104 definition is not finite: IGESGeom_ConicArc.cxx:241".into());
    }
    Ok(def)
}
