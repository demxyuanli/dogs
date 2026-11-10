//! `ComputeInternalPointsOnRstr` (`Contap_Contour.cxx:1120-1284`) and
//! `ComputeInternalPoints` (`cxx:1286-1539`): the points of a contour line
//! where the tangent is parallel to the view direction.

use occt_core::gp::{GpVec, GpVec2d};
use occt_core::precision::CONFUSION;
use occt_math::function_set_root::MathFunctionSetWithDerivatives;
use occt_math::{MathMatrix, MathVector};

use occt_topo::int_tools_wline::PntOn2S;
use occt_topo::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::super::line::ContapLine;
use super::super::point::ContapPoint;
use super::super::surf_function::ContapSurfFunction;
use super::super::surf_props;
use super::super::tools::{hcurve2d_d1, hcurve2d_resolution, nb_samples_on_arc};
use super::super::types::{ContapIType, ContapTFunction};

/// `gvec.XYZ() - other.XYZ()`.
fn sub_xyz(a: &occt_core::GpXyz, b: &occt_core::GpXyz) -> occt_core::GpXyz {
    occt_core::GpXyz::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

/// `gp::Resolution()` (`gp.cxx`).
const GP_RESOLUTION: f64 = 1.0e-15;

/// `ComputeInternalPointsOnRstr(Line, Paramf, Paraml, SFunc)`
/// (`cxx:1120-1284`).
pub(crate) fn compute_internal_points_on_rstr(
    line: &mut ContapLine,
    paramf: f64,
    paraml: f64,
    sfunc: &mut ContapSurfFunction,
) {
    if line.type_contour() != ContapIType::Restriction {
        return;
    }
    let Ok(thearc) = line.arc().cloned() else {
        return;
    };
    let Some(surf) = sfunc.surface().cloned() else {
        return;
    };
    let type_func = sfunc.function_type();

    let nbpnts = nb_samples_on_arc(&thearc);
    let mut indexinf = 1;
    let mut vecregard = GpVec::from_xyz(sfunc.direction().xyz());
    let toler = hcurve2d_resolution(&thearc, CONFUSION);
    let mut found = false;
    let mut vecref = GpVec::default();
    let mut paraminf;
    let mut pcour;

    // `cxx:1154-1184`.
    loop {
        paraminf = (f64::from(nbpnts - indexinf) * paramf
            + f64::from(indexinf - 1) * paraml)
            / f64::from(nbpnts - 1);
        let (p2d, d2d) = hcurve2d_d1(&thearc, paraminf);
        let Some((p, d1u, d1v, _n)) = surf_props::deriv_and_norm(&surf, p2d.x(), p2d.y()) else {
            return;
        };
        pcour = p;
        let tgt = combine(&d2d, &d1u, &d1v);

        if tgt.magnitude() > GP_RESOLUTION {
            if type_func == ContapTFunction::ContourPrs || type_func == ContapTFunction::DraftPrs {
                vecregard = GpVec::from_xyz(&sub_xyz(pcour.xyz(), sfunc.eye().xyz()));
            }
            let vr = vecregard.crossed(&tgt);
            if vr.magnitude() <= GP_RESOLUTION {
                indexinf += 1;
            } else {
                vecref = vr;
                found = true;
            }
        } else {
            indexinf += 1;
        }
        if found || indexinf > nbpnts {
            break;
        }
    }

    let mut indexsup = indexinf + 1;
    let mut paramsup;
    let mut toutvu = indexsup > nbpnts;
    // `ok` deliberately lives across the whole scan (`cxx:1136`).
    let mut ok = false;
    let mut paramp = 0.0;

    while !toutvu {
        paramsup = (f64::from(nbpnts - indexsup) * paramf + f64::from(indexsup - 1) * paraml)
            / f64::from(nbpnts - 1);
        let (p2d, d2d) = hcurve2d_d1(&thearc, paramsup);
        let Some((p, d1u, d1v, _n)) = surf_props::deriv_and_norm(&surf, p2d.x(), p2d.y()) else {
            return;
        };
        pcour = p;
        let tgt = combine(&d2d, &d1u, &d1v);

        let vectest = if tgt.magnitude() > GP_RESOLUTION {
            if type_func == ContapTFunction::ContourPrs || type_func == ContapTFunction::DraftPrs {
                vecregard = GpVec::from_xyz(&sub_xyz(pcour.xyz(), sfunc.eye().xyz()));
            }
            vecregard.crossed(&tgt)
        } else {
            GpVec::default()
        };

        if vectest.magnitude() <= GP_RESOLUTION {
            // Look for a genuine sign change.
            indexsup += 1;
        } else {
            if vectest.dot(&vecref) < 0.0 {
                // `cxx:1215-1254`: converge by halving.
                let mut solution = false;
                while !solution {
                    paramp = (paraminf + paramsup) / 2.0;
                    let (p2d, d2d) = hcurve2d_d1(&thearc, paramp);
                    let Some((p, d1u, d1v, _n)) =
                        surf_props::deriv_and_norm(&surf, p2d.x(), p2d.y())
                    else {
                        return;
                    };
                    pcour = p;
                    let tgt = combine(&d2d, &d1u, &d1v);
                    let vtestb = if tgt.magnitude() > GP_RESOLUTION {
                        if type_func == ContapTFunction::ContourPrs
                            || type_func == ContapTFunction::DraftPrs
                        {
                            vecregard = GpVec::from_xyz(&sub_xyz(pcour.xyz(), sfunc.eye().xyz()));
                        }
                        vecregard.crossed(&tgt)
                    } else {
                        GpVec::default()
                    };

                    if vtestb.magnitude() <= GP_RESOLUTION
                        || (paramp - paraminf).abs() <= toler
                        || (paramp - paramsup).abs() <= toler
                    {
                        solution = true;
                        ok = true;
                    } else if vtestb.dot(&vecref) < 0.0 {
                        paramsup = paramp;
                    } else {
                        paraminf = paramp;
                    }
                }

                if ok {
                    // `cxx:1256-1277`: reject a parameter already carried by a
                    // vertex, otherwise add the internal point.
                    let nbvtx = line.nb_vertex();
                    for i in 1..=nbvtx {
                        let has = line
                            .vertex(i)
                            .map(|v| (v.parameter_on_line() - paramp).abs() <= toler)
                            .unwrap_or(false);
                        if has {
                            if let Some(v) = line.vertex_mut(i) {
                                v.set_internal();
                            }
                            ok = false;
                        }
                    }
                    if ok {
                        let p2d = hcurve2d_d1(&thearc, paramp).0;
                        let mut internalp = ContapPoint::with_uv(pcour, p2d.x(), p2d.y());
                        internalp.set_parameter(paramp);
                        internalp.set_internal();
                        line.add_vertex(internalp);
                    }
                }
                // `cxx:1269` also recomputes `paramsup` here, but `cxx:1181`
                // overwrites it at the top of the next iteration before any read.
            }
            vecref = vectest;
            // `cxx:1271` also stores `indexinf = indexsup`; `indexinf` is never
            // read again in this loop, so that store is dead in OCCT as well.
            indexsup += 1;
        }
        toutvu = indexsup > nbpnts;
    }
}

/// `ComputeInternalPoints(Line, SFunc, ureso, vreso)` (`cxx:1286-1539`).
pub(crate) fn compute_internal_points(
    line: &mut ContapLine,
    sfunc: &mut ContapSurfFunction,
    ureso: f64,
    vreso: f64,
) {
    if line.type_contour() != ContapIType::Walking {
        return;
    }
    let Ok(mut nbpnts) = line.nb_pnts() else {
        return;
    };
    let Some(surf) = sfunc.surface().cloned() else {
        return;
    };
    let type_func = sfunc.function_type();

    let mut toler = MathVector::new(1, 2);
    toler.set_value(1, ureso);
    toler.set_value(2, vreso);
    let mut infb = MathVector::new(1, 2);
    let mut supb = MathVector::new(1, 2);
    let (fu, fv) = surface_first_parameters(&surf);
    let (lu, lv) = surface_last_parameters(&surf);
    infb.set_value(1, fu);
    infb.set_value(2, fv);
    supb.set_value(1, lu);
    supb.set_value(2, lv);

    let mut rsnld = occt_math::MathFunctionSetRoot::new(&*sfunc, &toler, 30);

    let mut indexinf = 1;
    let mut vecregard = GpVec::from_xyz(sfunc.direction().xyz());
    let mut found = false;
    let mut vecref = GpVec::default();

    let mut x_inf = MathVector::new(1, 2);
    let mut x_sup = MathVector::new(1, 2);
    let mut x = MathVector::new(1, 2);
    let mut f = MathVector::new(1, 1);
    let mut df = MathMatrix::new(1, 1, 1, 2);

    // `cxx:1330-1356`.
    loop {
        let (u, v) = match line.point(indexinf as usize) {
            Ok(p) => p.parameters_on_s2(),
            Err(_) => return,
        };
        x_inf.set_value(1, u);
        x_inf.set_value(2, v);
        sfunc.values(&x_inf, &mut f, &mut df);
        if !sfunc.is_tangent() {
            let tgt = match sfunc.direction3d() {
                Ok(t) => *t,
                Err(_) => return,
            };
            if type_func == ContapTFunction::ContourPrs || type_func == ContapTFunction::DraftPrs {
                let pv = line.point(indexinf as usize).map(|p| p.value());
                if let Ok(pv) = pv {
                    vecregard = GpVec::from_xyz(&sub_xyz(pv.xyz(), sfunc.eye().xyz()));
                }
            }
            let vr = vecregard.crossed(&tgt);
            if vr.magnitude() <= GP_RESOLUTION {
                indexinf += 1;
            } else {
                vecref = vr;
                found = true;
            }
        } else {
            indexinf += 1;
        }
        if found || indexinf > nbpnts {
            break;
        }
    }

    let mut indexsup = indexinf + 1;
    let mut toutvu = indexsup > nbpnts;
    let mut ok = false;

    while !toutvu {
        let (u, v) = match line.point(indexsup as usize) {
            Ok(p) => p.parameters_on_s2(),
            Err(_) => return,
        };
        x_sup.set_value(1, u);
        x_sup.set_value(2, v);
        sfunc.values(&x_sup, &mut f, &mut df);

        let vectest = if !sfunc.is_tangent() {
            let tgt = match sfunc.direction3d() {
                Ok(t) => *t,
                Err(_) => return,
            };
            if type_func == ContapTFunction::ContourPrs || type_func == ContapTFunction::DraftPrs {
                let pv = line.point(indexsup as usize).map(|p| p.value());
                if let Ok(pv) = pv {
                    vecregard = GpVec::from_xyz(&sub_xyz(pv.xyz(), sfunc.eye().xyz()));
                }
            }
            vecregard.crossed(&tgt)
        } else {
            GpVec::default()
        };

        if vectest.magnitude() <= GP_RESOLUTION {
            indexsup += 1;
        } else {
            if vectest.dot(&vecref) < 0.0 {
                // `cxx:1389-1452`.
                let mut solution = false;
                while !solution {
                    // `cxx:1394`: a shifted point instead of the middle.
                    let acoef = 2.0 / 3.0;
                    x.set_value(1, x_inf.value(1) + acoef * (x_sup.value(1) - x_inf.value(1)));
                    x.set_value(2, x_inf.value(2) + acoef * (x_sup.value(2) - x_inf.value(2)));
                    rsnld.perform_with_bounds(&mut *sfunc, &x, &infb, &supb, false);

                    if !rsnld.is_done() {
                        // `cxx:1402`: `std::cout << "Echec recherche internal points"`.
                        solution = true;
                        ok = false;
                    } else {
                        x = rsnld.root();
                        sfunc.values(&x, &mut f, &mut df);
                        if f.value(1).abs() <= sfunc.tolerance() {
                            let vtestb = if !sfunc.is_tangent() {
                                let tgt = match sfunc.direction3d() {
                                    Ok(t) => *t,
                                    Err(_) => return,
                                };
                                if type_func == ContapTFunction::ContourPrs
                                    || type_func == ContapTFunction::DraftPrs
                                {
                                    vecregard = GpVec::from_xyz(&sub_xyz(sfunc.point().xyz(), sfunc.eye().xyz()));
                                }
                                vecregard.crossed(&tgt)
                            } else {
                                GpVec::default()
                            };
                            if vtestb.magnitude() <= GP_RESOLUTION
                                || ((x.value(1) - x_inf.value(1)).abs() <= toler.value(1)
                                    && (x.value(2) - x_inf.value(2)).abs() <= toler.value(2))
                                || ((x.value(1) - x_sup.value(1)).abs() <= toler.value(1)
                                    && (x.value(2) - x_sup.value(2)).abs() <= toler.value(2))
                            {
                                solution = true;
                                ok = true;
                            } else if vtestb.dot(&vecref) < 0.0 {
                                x_sup = x.clone();
                            } else {
                                x_inf = x.clone();
                            }
                        } else {
                            // `cxx:1447`: `std::cout << "Echec recherche internal points"`.
                            solution = true;
                            ok = false;
                        }
                    }
                }

                if ok {
                    // `cxx:1454-1529`.
                    let mut newpoint = false;
                    let (u, v) = line
                        .point(indexinf as usize)
                        .map(|p| p.parameters_on_s2())
                        .unwrap_or((0.0, 0.0));
                    let mut vinf = GpVec2d::new(x.value(1) - u, x.value(2) - v);
                    let mut paramp;
                    if vinf.x().abs() <= toler.value(1) && vinf.y().abs() <= toler.value(2) {
                        paramp = f64::from(indexinf);
                    } else {
                        paramp = 0.0;
                        let mut index = indexinf + 1;
                        while index <= indexsup {
                            let (u, v) = line
                                .point(index as usize)
                                .map(|p| p.parameters_on_s2())
                                .unwrap_or((0.0, 0.0));
                            let vsup = GpVec2d::new(x.value(1) - u, x.value(2) - v);
                            if vsup.x().abs() <= toler.value(1)
                                && vsup.y().abs() <= toler.value(2)
                            {
                                paramp = f64::from(index);
                                break;
                            } else if vinf.dot(&vsup) < 0.0 {
                                // Between the two points.
                                paramp = f64::from(index);
                                let mut pt2s = PntOn2S { p: occt_core::gp::GpPnt::default(), u1: 0.0, v1: 0.0, u2: 0.0, v2: 0.0 };
                                pt2s.set_point(*sfunc.point());
                                pt2s.set_uv_on(false, x.value(1), x.value(2));
                                if let Some(wl) = line.line_on_2s_mut() {
                                    wl.insert_before(index, pt2s);
                                }
                                // Shift the vertices after the inserted point.
                                let nbvtx = line.nb_vertex();
                                for vi in 1..=nbvtx {
                                    let shift = line
                                        .vertex(vi)
                                        .map(|v| v.parameter_on_line() >= f64::from(index))
                                        .unwrap_or(false);
                                    if shift {
                                        if let Some(v) = line.vertex_mut(vi) {
                                            let p = v.parameter_on_line() + 1.0;
                                            v.set_parameter(p);
                                        }
                                    }
                                }
                                nbpnts += 1;
                                indexsup += 1;
                                newpoint = true;
                                break;
                            } else {
                                vinf = vsup;
                            }
                            index += 1;
                        }
                    }

                    if !newpoint {
                        // On a walking point: match an existing vertex.
                        newpoint = true;
                        let nbvtx = line.nb_vertex();
                        for v in 1..=nbvtx {
                            let same = line
                                .vertex(v)
                                .map(|vp| vp.parameter_on_line() == paramp)
                                .unwrap_or(false);
                            if same {
                                if let Some(vp) = line.vertex_mut(v) {
                                    vp.set_internal();
                                }
                                newpoint = false;
                            }
                        }
                    }

                    if newpoint && paramp > 1.0 && paramp < f64::from(nbpnts) {
                        let mut internalp =
                            ContapPoint::with_uv(*sfunc.point(), x.value(1), x.value(2));
                        internalp.set_parameter(paramp);
                        internalp.set_internal();
                        line.add_vertex(internalp);
                    }
                }

                // `cxx:1531`: re-read the upper bound point.
                let (u, v) = match line.point(indexsup as usize) {
                    Ok(p) => p.parameters_on_s2(),
                    Err(_) => return,
                };
                x_sup.set_value(1, u);
                x_sup.set_value(2, v);
            }
            vecref = vectest;
            indexinf = indexsup;
            indexsup += 1;
            x_inf = x_sup.clone();
        }
        toutvu = indexsup > nbpnts;
    }
}

/// `tgt.SetLinearForm(d2d.X(), d1u, d2d.Y(), d1v)`.
fn combine(d2d: &GpVec2d, d1u: &GpVec, d1v: &GpVec) -> GpVec {
    GpVec::new(
        d2d.x() * d1u.x() + d2d.y() * d1v.x(),
        d2d.x() * d1u.y() + d2d.y() * d1v.y(),
        d2d.x() * d1u.z() + d2d.y() * d1v.z(),
    )
}

/// `Adaptor3d_HSurfaceTool::FirstUParameter/FirstVParameter` on the surface's
/// own bounds (`Adaptor3d_HSurfaceTool.cxx:151-167`).
fn surface_first_parameters(s: &SurfaceAdaptor) -> (f64, f64) {
    (s.first_u_parameter(), s.first_v_parameter())
}

/// `Adaptor3d_HSurfaceTool::LastUParameter/LastVParameter`.
fn surface_last_parameters(s: &SurfaceAdaptor) -> (f64, f64) {
    (s.last_u_parameter(), s.last_v_parameter())
}
