//! The three `Update` statics of `ChFi3d_Builder_C1.cxx` and their caller
//! `IntersUpdateOnSame`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C1.cxx:118-418`
//! (OCCT 8.0.0).
//!
//! All four functions re-locate a tangency point (`ChFiDS_CommonPoint` on an
//! arc + the parameter of the `ChFiDS_FaceInterference` pcurve on the fillet
//! surface) after one of the two faces was extended or its pcurve re-read.
//! Each answers `false` when no intersection is found, exactly like the OCCT
//! `bool` return; the port's `Err` arms only cover OCCT's throws / null
//! dereferences.

use std::sync::Arc;

use occt_core::elib::clib::in_period;
use occt_core::gp::GpPnt2d;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::INFINITE;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::extrema_cc::ExtremaLocateExtCC;
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::chfi3d::tangent_faces::is_tangent_faces;
use crate::chfi_ds::{CommonPoint, FaceInterference};
use crate::intcurvesurface::perform_curve_surface;
use crate::shape::{Edge, Face, TopoShape};
use crate::topopebrep::tool::adaptor::GeomAdaptorCurve;

/// `recadre(p, ref, isfirst, first, last)` (`cxx:118-138`): adds or subtracts
/// one period and keeps the candidate that is closer to `ref`.
fn recadre(p: f64, reference: f64, is_first: bool, first: f64, last: f64) -> f64 {
    let mut pp = p;
    if is_first {
        pp -= last - first;
    } else {
        pp += last - first;
    }
    if (pp - reference).abs() < (p - reference).abs() {
        return pp;
    }
    p
}

/// `Update(fb, pcfb, surf, fi, cp, p2dbout, isfirst, pared, wop, tol)`
/// (`cxx:147-181`).
///
/// `Extrema_LocateExtCC` between the 3D curve `(pcfb, fb)` and the 3D curve
/// `(fi.PCurveOnSurf(), surf)`; on a hit inside `tol` it writes the pcurve point
/// at the first curve's parameter into `p2dbout`, stores the two parameters,
/// moves the interference parameter and rebuilds the common point.
///
/// Divergence: `pc = fi.PCurveOnSurf()` can be a null handle, which
/// `Geom2dAdaptor_Curve(pc)` would accept but whose `Value` would fault; the
/// port answers `false` for it.
pub fn update_pcurve_on_pcurve(
    fb: &Arc<dyn Surface>,
    pcfb: &Arc<dyn Curve2d>,
    surf: &Arc<dyn Surface>,
    fi: &mut FaceInterference,
    cp: &mut CommonPoint,
    p2dbout: &mut GpPnt2d,
    is_first: bool,
    pared: &mut f64,
    wop: &mut f64,
    tol: f64,
) -> bool {
    // `Adaptor3d_CurveOnSurface c1(pcfb, fb);`
    let c1 = CurveOnSurface::new(Arc::clone(pcfb), Arc::clone(fb));
    // `occ::handle<Geom2d_Curve> pc = fi.PCurveOnSurf();`
    // `occ::handle<Geom2dAdaptor_Curve> hpc = new Geom2dAdaptor_Curve(pc);`
    let Some(pc) = fi.p_curve_on_surf().cloned() else {
        return false;
    };
    // `Adaptor3d_CurveOnSurface c2(hpc, surf);`
    let c2 = CurveOnSurface::new(pc, Arc::clone(surf));

    // `Extrema_LocateExtCC ext(c1, c2, pared, wop);`
    let ext = ExtremaLocateExtCC::new(&c1, &c2, *pared, *wop);
    if ext.is_done() {
        let dist2 = ext.square_distance();
        if dist2 < tol * tol {
            let (ponc1, ponc2) = ext.point();
            let parfb = ponc1.parameter();
            // `p2dbout = pcfb->Value(parfb);`
            *p2dbout = pcfb.d0(parfb);
            *pared = ponc1.parameter();
            *wop = ponc2.parameter();
            fi.set_parameter(*wop, is_first);
            cp.reset();
            cp.set_point(&ponc1.value());
            return true;
        }
    }
    false
}

/// `Update(fb, ct, fi, cp, p2dbout, isfirst, wop)` (`cxx:191-278`).
///
/// `IntCurveSurface_HInter::Perform(ct, fb)` then the two passes of `cxx:213-241`
/// (nearest to `wop` inside the curve range; for a periodic curve a second pass
/// preferring a point within `0.01` of either range end), and the branch that
/// reads back `Point(isol)` / `Point(isolbis)`.
///
/// Faithfulness notes kept from the source:
/// - the first pass folds every candidate with `recadre`, the second compares
///   the **raw** `w` but stores the recadre'd one,
/// - the `!recadrebis` arm folds with `ElCLib::InPeriod` only when the curve is
///   periodic, while the `dist <= distbis` arm of `recadrebis` folds
///   **unconditionally** (`cxx:269`),
/// - `dist` / `distbis` start at `Precision::Infinite()`.
pub fn update_surface_curve(
    fb: &Arc<dyn Surface>,
    ct: &dyn Curve,
    fi: &mut FaceInterference,
    cp: &mut CommonPoint,
    p2dbout: &mut GpPnt2d,
    is_first: bool,
    wop: &mut f64,
) -> bool {
    // `double w, uf = ct->FirstParameter(), ul = ct->LastParameter();`
    let uf = ct.first_parameter();
    let ul = ct.last_parameter();
    // `bool isperiodic = ct->IsPeriodic(), recadrebis = false;`
    let is_periodic = ct.is_periodic();
    let mut recadrebis = false;

    // `Intersection.Perform(ct, fb);` - the port's `HInter` needs the surface
    // window, which `GeomAdaptor_Surface` answers with the surface's own bounds.
    let (u0, u1) = fb.u_range();
    let (v0, v1) = fb.v_range();
    let Ok(intersection) = perform_curve_surface(ct, fb.as_ref(), (uf, ul), (u0, u1, v0, v1))
    else {
        return false;
    };

    let nbp = intersection.nb_points();
    let mut isol = 0usize;
    let mut isolbis = 0usize;
    let mut dist = INFINITE;
    let mut distbis = INFINITE;
    let mut wbis = 0.0f64;

    // `for (i = 1; i <= nbp; i++)` first pass (`cxx:214-227`).
    for i in 1..=nbp {
        let mut w = intersection.point(i - 1).param;
        if is_periodic {
            w = recadre(w, *wop, is_first, uf, ul);
        }
        if uf <= w && ul >= w && (w - *wop).abs() < dist {
            isol = i;
            dist = (w - *wop).abs();
        }
    }
    if is_periodic {
        // Second pass (`cxx:228-241`).
        for i in 1..=nbp {
            let w = intersection.point(i - 1).param;
            if uf <= w
                && ul >= w
                && (w - *wop).abs() < distbis
                && ((w - ul).abs() <= 0.01 || (w - uf).abs() <= 0.01)
            {
                isolbis = i;
                wbis = recadre(w, *wop, is_first, uf, ul);
                distbis = (wbis - *wop).abs();
                recadrebis = true;
            }
        }
    }
    if isol == 0 && isolbis == 0 {
        return false;
    }

    let w = if !recadrebis {
        // `IntCurveSurface_IntersectionPoint pint = Intersection.Point(isol);`
        let pint = intersection.point(isol - 1);
        p2dbout.set_coord(pint.u, pint.v);
        let mut w = pint.param;
        if is_periodic {
            w = in_period(w, uf, ul);
        }
        w
    } else if dist > distbis {
        let pint = intersection.point(isolbis - 1);
        p2dbout.set_coord(pint.u, pint.v);
        wbis
    } else {
        let pint = intersection.point(isol - 1);
        p2dbout.set_coord(pint.u, pint.v);
        let w = pint.param;
        // `w = ElCLib::InPeriod(w, uf, ul);` - unconditional in the source.
        in_period(w, uf, ul)
    };

    fi.set_parameter(w, is_first);
    cp.reset();
    cp.set_point(&ct.d0(w));
    *wop = w;
    true
}

/// `Update(face, edonface, surf, fi, cp, isfirst)` (`cxx:349-417`): the extrema
/// between `(edonface, face)` and the fillet pcurve `(fi.PCurveOnSurf(), surf)`
/// started at the common point's arc parameter, then the `parltg` clamp against
/// the **pcurve on the face's** range (`cxx:381-410`).
///
/// `Geom2dAdaptor_Curve(pc, f, l)` (`cxx:369`) is a pcurve window and becomes
/// [`Geom2dTrimmedCurve`] with the same `[f, l]`.
pub fn update_curve_on_surface(
    face: &Arc<dyn Surface>,
    edonface: &Arc<dyn Curve2d>,
    surf: &Arc<dyn Surface>,
    fi: &mut FaceInterference,
    cp: &mut CommonPoint,
    is_first: bool,
) -> bool {
    if !cp.is_on_arc() {
        return false;
    }
    // `Adaptor3d_CurveOnSurface c1(edonface, face);`
    let c1 = CurveOnSurface::new(Arc::clone(edonface), Arc::clone(face));
    // `double pared = cp.ParameterOnArc();` - `is_on_arc` means `SetArc` ran,
    // so the member is set; a missing one answers `false` instead of reading
    // OCCT's uninitialised `prmarc`.
    let Some(mut pared) = cp.parameter_on_arc() else {
        return false;
    };
    let mut parltg = fi.parameter(is_first);
    let Some(pc) = fi.p_curve_on_surf().cloned() else {
        return false;
    };
    let mut f = fi.first_parameter();
    let mut l = fi.last_parameter();
    let delta = 0.1 * (l - f);
    f = (f - delta).max(pc.first_parameter());
    l = (l + delta).min(pc.last_parameter());
    // `new Geom2dAdaptor_Curve(pc, f, l)`.
    let hpc: Arc<dyn Curve2d> = Arc::new(Geom2dTrimmedCurve::new(Arc::clone(&pc), f, l));
    let c2 = CurveOnSurface::new(hpc, Arc::clone(surf));

    // `Extrema_LocateExtCC ext(c1, c2, pared, parltg);`
    let ext = ExtremaLocateExtCC::new(&c1, &c2, pared, parltg);
    if !ext.is_done() {
        return false;
    }
    let (ponc1, ponc2) = ext.point();
    pared = ponc1.parameter();
    parltg = ponc2.parameter();
    if !(parltg > f && parltg < l) {
        return false;
    }

    // `if (!PConF.IsNull()) { aTrCurve = down_cast<Geom2d_TrimmedCurve>(PConF);
    // if (!aTrCurve.IsNull()) PConF = aTrCurve->BasisCurve();
    // if (!PConF->IsPeriodic()) { ... clamp parltg ... } }` (`cxx:382-410`).
    if let Some(pc_face) = fi.p_curve_on_face() {
        let base = pc_face.trimmed_basis().unwrap_or(pc_face.as_ref());
        if !base.is_periodic() {
            if is_first {
                let fpar = base.first_parameter();
                if parltg < fpar {
                    parltg = fpar;
                }
            } else {
                let lpar = base.last_parameter();
                if parltg > lpar {
                    parltg = lpar;
                }
            }
        }
    }

    fi.set_parameter(parltg, is_first);
    // `cp.SetArc(cp.Tolerance(), cp.Arc(), pared, cp.TransitionOnArc());`
    let tol = cp.tolerance();
    let Some(arc) = cp.arc().cloned() else {
        return false;
    };
    // OCCT leaves `myTraArc` uninitialised in the constructor and `SetArc`
    // always writes it, so the only way to get here without it is a
    // `CommonPoint` built outside `SetArc`; answer `false` instead of reading
    // an indeterminate value.
    let Some(traarc) = cp.transition_on_arc() else {
        return false;
    };
    cp.set_arc(tol, &arc, pared, traarc);
    true
}

/// `IntersUpdateOnSame(HGs, HBs, c3dFI, Fop, Fprol, Eprol, Vtx, isFirst, Tol,
/// FIop, CPop, FprolUV, c3dU)` (`cxx:292-340`).
///
/// Tries the 3D curve of the `Fop` interference against the extended `Fprol`
/// surface first; when that fails and the two faces are tangent along `Eprol`,
/// retries with `Eprol`'s pcurve on `Fprol` against the fillet surface `HGs`,
/// started at the vertex parameter.
pub fn inter_update_on_same(
    hgs: &Arc<dyn Surface>,
    hbs: &Arc<dyn Surface>,
    c3d_fi: &Arc<dyn Curve>,
    fop: &Face,
    fprol: &Face,
    eprol: &Edge,
    vtx: &TopoShape,
    is_first: bool,
    tol: f64,
    fiop: &mut FaceInterference,
    cpop: &mut CommonPoint,
    fprol_uv: &mut GpPnt2d,
    c3d_u: &mut f64,
) -> Result<bool, String> {
    let uf = fiop.first_parameter();
    let ul = fiop.last_parameter();
    // `if (c3dFI->IsPeriodic()) Hc3df = new GeomAdaptor_Curve(c3dFI);
    //  else Hc3df = new GeomAdaptor_Curve(c3dFI, uf, ul);`
    let hc3df: Arc<dyn Curve> = if c3d_fi.is_periodic() {
        Arc::new(GeomAdaptorCurve::of_curve(Arc::clone(c3d_fi)))
    } else {
        Arc::new(GeomAdaptorCurve::with_range(Arc::clone(c3d_fi), uf, ul))
    };

    if update_surface_curve(
        hbs,
        hc3df.as_ref(),
        fiop,
        cpop,
        fprol_uv,
        is_first,
        c3d_u,
    ) {
        return Ok(true);
    }

    // `if (!ChFi3d::IsTangentFaces(Eprol, Fprol, Fop)) return false;`
    // (`hxx:44-48` default `theOrder = GeomAbs_G1`).
    if !is_tangent_faces(eprol, fprol, fop, Shape::G1)? {
        return Ok(false);
    }

    // `gpcprol = BRep_Tool::CurveOnSurface(Eprol, Fprol, uf, ul);`
    let Some((pcprol, _, _)) = curve_on_surface_range(eprol, fprol) else {
        return Err("Failed to get p-curve of edge".to_string());
    };
    // `double partemp = BRep_Tool::Parameter(Vtx, Eprol);`
    let mut partemp = BRepTool::parameter_of_vertex(&eprol.0, vtx);

    Ok(update_pcurve_on_pcurve(
        hbs, &pcprol, hgs, fiop, cpop, fprol_uv, is_first, &mut partemp, c3d_u, tol,
    ))
}
