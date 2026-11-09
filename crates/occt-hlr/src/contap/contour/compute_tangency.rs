//! `ComputeTangency` (`Contap_Contour.cxx:556-891`): turns the points found on
//! the restrictions into `IntSurf_PathPoint` start points for the walker, and
//! drops the ones that a `TopTrans_CurveTransition` rejects.

use occt_core::gp::{GpDir, GpDir2d, GpVec, GpVec2d};
use occt_core::precision::CONFUSION;
use occt_math::function_set_root::MathFunctionSetWithDerivatives;
use occt_math::{MathMatrix, MathVector};

use occt_topo::abs::{Orientation, State};
use occt_topo::brep_top_adaptor_topol_tool::TopolTool;
use occt_topo::top_trans_curve_transition::CurveTransition;

use super::super::iwalking::types::PathPoint;
use super::super::surf_function::ContapSurfFunction;
use super::super::surf_props;
use super::super::the_search::ContapTheSearch;
use super::super::tools::{hcont_project, hcurve2d_d1, hcurve2d_value};
use super::TOLE;

/// `ComputeTangency(solrst, Domain, Func, seqpdep, Destination)`
/// (`cxx:556-891`).
///
/// `destination` is OCCT's `NCollection_Array1<int> Destination(1, NbPointRst+1)`:
/// index `i` (1-based) is written by this routine and read by
/// `Contap_Contour::Perform`; index 0 is unused.
pub(crate) fn compute_tangency(
    solrst: &ContapTheSearch,
    _domain: &TopolTool,
    func: &mut ContapSurfFunction,
    seqpdep: &mut Vec<PathPoint>,
    destination: &mut [i32],
) {
    let Some(surf) = func.surface().cloned() else {
        return;
    };
    let nb_points = solrst.nb_points().unwrap_or(0) as i32;
    let mut seqlength: i32 = 0;

    let mut x = MathVector::new(1, 2);
    let mut f = MathVector::new(1, 1);
    let mut d = MathMatrix::new(1, 1, 1, 2);

    for i in 1..=nb_points {
        if destination[i as usize] != 0 {
            continue;
        }
        let Some(pstart) = solrst.point(i as usize) else {
            continue;
        };
        let thearc = pstart.arc();
        let mut theparam = pstart.parameter();
        let ptoproj = hcurve2d_value(thearc, theparam);

        // `cxx:596-615`: drop the points already present on a solution
        // restriction.
        let mut sur_une_restriction_solution = false;
        let nba = solrst.nb_segments().unwrap_or(0);
        for restriction in 1..=nba {
            if sur_une_restriction_solution {
                break;
            }
            let Some(seg) = solrst.segment(restriction) else {
                continue;
            };
            let thearcsol = seg.curve();
            if let Some((_paramproj, pproj)) = hcont_project(thearcsol, &ptoproj) {
                let pprojete = surf.value(pproj.x(), pproj.y()).unwrap_or_default();
                if pstart.value().distance(&pprojete) <= CONFUSION {
                    sur_une_restriction_solution = true;
                }
            }
        }
        if sur_une_restriction_solution {
            continue;
        }

        // `cxx:618-620`.
        let mut arcorien = TopolTool::orientation_of_curve(thearc);
        let mut ispassing = arcorien == Orientation::Internal || arcorien == Orientation::External;

        let (pt2d, mut tg2drst) = hcurve2d_d1(thearc, theparam);
        x.set_value(1, pt2d.x());
        x.set_value(2, pt2d.y());
        let mut ppoint = PathPoint::with_uv(*pstart.value(), x.value(1), x.value(2));

        func.values(&x, &mut f, &mut d);
        if func.is_tangent() {
            // `cxx:628-664`.
            ppoint.set_tangency(true);
            destination[i as usize] = seqlength + 1;
            if !pstart.is_new() {
                if let Some(vtx) = pstart.vertex() {
                    for k in (i + 1)..=nb_points {
                        if destination[k as usize] != 0 {
                            continue;
                        }
                        let Some(pstart2) = solrst.point(k as usize) else {
                            continue;
                        };
                        if pstart2.is_new() {
                            continue;
                        }
                        let Some(vtx2) = pstart2.vertex() else {
                            continue;
                        };
                        if !TopolTool::identical(vtx, vtx2) {
                            continue;
                        }
                        let thearc2 = pstart2.arc();
                        theparam = pstart2.parameter();
                        arcorien = TopolTool::orientation_of_curve(thearc2);
                        ispassing = ispassing
                            && (arcorien == Orientation::Internal
                                || arcorien == Orientation::External);

                        let pt2d = hcurve2d_value(thearc2, theparam);
                        x.set_value(1, pt2d.x());
                        x.set_value(2, pt2d.y());
                        ppoint.add_uv(x.value(1), x.value(2));
                        destination[k as usize] = seqlength + 1;
                    }
                }
            }
            ppoint.set_passing(ispassing);
            seqpdep.push(ppoint);
            seqlength += 1;
        } else {
            // `cxx:666-888`: a potential start point.
            let mut vectg = func.direction3d().map(|v| *v).unwrap_or_default();
            let mut dirtg = func
                .direction2d()
                .unwrap_or_else(|_| GpDir2d::new(1.0, 0.0).expect("+X"));

            let Some((_ptbid, v1, v2, normale)) = surf_props::deriv_and_norm(&surf, x.value(1), x.value(2))
            else {
                continue;
            };
            let mut tg3drst = d1u_plus_d1v(&tg2drst, &v1, &v2);
            if normale.square_magnitude() < occt_core::precision::REAL_EPSILON {
                // `cxx:676-679`: a null normal only prints a diagnostic.
            } else {
                let mut test = vectg.dot(&normale.crossed(&tg3drst));

                if pstart.is_new() {
                    // `cxx:685-707`.
                    let tbis = vectg.normalized().dot(&tg3drst.normalized());
                    if tbis.abs() < 1.0 - TOLE {
                        if (test < 0.0 && arcorien == Orientation::Forward)
                            || (test > 0.0 && arcorien == Orientation::Reversed)
                        {
                            vectg.reverse();
                            dirtg.reverse();
                        }
                        ppoint.set_directions(vectg, dirtg);
                    } else {
                        // keep the point as a tangency stop point
                        ppoint.set_tangency(true);
                    }
                    ppoint.set_passing(ispassing);
                    destination[i as usize] = seqlength + 1;
                    seqpdep.push(ppoint);
                    seqlength += 1;
                } else {
                    // `cxx:708-888`: the complex transition.
                    let bidnorm = GpDir::new(1.0, 1.0, 1.0).expect("(1,1,1)");
                    let mut tobeverified = false;
                    let mut loc_trans;
                    let mut comptrans = CurveTransition::new();
                    comptrans.reset(
                        &GpDir::from_vec(&vectg).unwrap_or(bidnorm),
                        &bidnorm,
                        0.0,
                    );
                    if arcorien != Orientation::Internal && arcorien != Orientation::External {
                        let vtx = pstart.vertex();
                        let vtxorien = vtx
                            .map(|v| v.orientation())
                            .unwrap_or(Orientation::Forward);
                        test = test / vectg.magnitude();
                        test = test / normale.crossed(&tg3drst).magnitude();

                        if test.abs() <= TOLE {
                            tobeverified = true;
                            loc_trans = Orientation::External;
                        } else {
                            if (test > 0.0 && arcorien == Orientation::Forward)
                                || (test < 0.0 && arcorien == Orientation::Reversed)
                            {
                                loc_trans = Orientation::Forward;
                            } else {
                                loc_trans = Orientation::Reversed;
                            }
                            if arcorien == Orientation::Reversed {
                                tg3drst.reverse();
                            }
                        }
                        comptrans.compare(
                            TOLE,
                            &GpDir::from_vec(&tg3drst).unwrap_or(bidnorm),
                            &bidnorm,
                            0.0,
                            loc_trans,
                            vtxorien,
                        );
                    }
                    destination[i as usize] = seqlength + 1;
                    for k in (i + 1)..=nb_points {
                        if destination[k as usize] != 0 {
                            continue;
                        }
                        let Some(pstart2) = solrst.point(k as usize) else {
                            continue;
                        };
                        if pstart2.is_new() {
                            continue;
                        }
                        let (Some(vtx), Some(vtx2)) = (pstart.vertex(), pstart2.vertex()) else {
                            continue;
                        };
                        if !TopolTool::identical(vtx, vtx2) {
                            continue;
                        }
                        let thearc2 = pstart2.arc();
                        theparam = pstart2.parameter();
                        arcorien = TopolTool::orientation_of_curve(thearc2);

                        let (_pt2d, tg2) = hcurve2d_d1(thearc2, theparam);
                        tg2drst = tg2;
                        x.set_value(1, pt2d.x());
                        x.set_value(2, pt2d.y());
                        ppoint.add_uv(x.value(1), x.value(2));

                        if arcorien != Orientation::Internal
                            && arcorien != Orientation::External
                        {
                            ispassing = false;
                            tg3drst = d1u_plus_d1v(&tg2drst, &v1, &v2);
                            test = vectg.dot(&normale.crossed(&tg3drst));
                            test = test / vectg.magnitude();
                            test = test / normale.crossed(&tg3drst).magnitude();

                            let vtxorien = vtx2.orientation();
                            if test.abs() <= TOLE {
                                tobeverified = true;
                                loc_trans = Orientation::External;
                            } else {
                                if (test > 0.0 && arcorien == Orientation::Forward)
                                    || (test < 0.0 && arcorien == Orientation::Reversed)
                                {
                                    loc_trans = Orientation::Forward;
                                } else {
                                    loc_trans = Orientation::Reversed;
                                }
                                if arcorien == Orientation::Reversed {
                                    tg3drst.reverse();
                                }
                            }
                            comptrans.compare(
                                TOLE,
                                &GpDir::from_vec(&tg3drst).unwrap_or(bidnorm),
                                &bidnorm,
                                0.0,
                                loc_trans,
                                vtxorien,
                            );
                        }
                        destination[k as usize] = seqlength + 1;
                    }

                    // `cxx:805-833`.
                    let mut fairpt = true;
                    if !ispassing {
                        let before = comptrans.state_before();
                        let after = comptrans.state_after();
                        if before == State::Unknown || after == State::Unknown {
                            fairpt = false;
                        } else if before == State::In {
                            if after == State::In {
                                ispassing = true;
                            } else {
                                vectg.reverse();
                                dirtg.reverse();
                            }
                        } else if after != State::In {
                            fairpt = false;
                        }
                    }

                    // `cxx:836-867`: do not start along a solution restriction.
                    if fairpt && tobeverified {
                        for k in i..=nb_points {
                            if destination[k as usize] != seqlength + 1 {
                                continue;
                            }
                            let Some(pk) = solrst.point(k as usize) else {
                                continue;
                            };
                            theparam = pk.parameter();
                            let thearc2 = pk.arc();
                            arcorien = TopolTool::orientation_of_curve(thearc2);

                            if arcorien == Orientation::Forward
                                || arcorien == Orientation::Reversed
                            {
                                let (_pt2d, tg2) = hcurve2d_d1(thearc2, theparam);
                                tg3drst = d1u_plus_d1v(&tg2, &v1, &v2);
                                let vtxorien = pk
                                    .vertex()
                                    .map(|v| v.orientation())
                                    .unwrap_or(Orientation::Forward);
                                if (arcorien == Orientation::Forward
                                    && vtxorien == Orientation::Reversed)
                                    || (arcorien == Orientation::Reversed
                                        && vtxorien == Orientation::Forward)
                                {
                                    tg3drst.reverse();
                                }
                                test = vectg.normalized().dot(&tg3drst.normalized());
                                if test >= 1.0 - TOLE {
                                    fairpt = false;
                                    break;
                                }                            }
                        }
                    }

                    // `cxx:868-888`.
                    if fairpt {
                        ppoint.set_directions(vectg, dirtg);
                        ppoint.set_passing(ispassing);
                        seqpdep.push(ppoint);
                        seqlength += 1;
                    } else {
                        for k in i..=nb_points {
                            if destination[k as usize] == seqlength + 1 {
                                destination[k as usize] = -destination[k as usize];
                            }
                        }
                    }
                }
            }
        }
    }
}

/// `tg3drst = tg2drst.X() * v1 + tg2drst.Y() * v2` (`SetLinearForm`).
fn d1u_plus_d1v(tg2drst: &GpVec2d, v1: &GpVec, v2: &GpVec) -> GpVec {
    GpVec::new(
        tg2drst.x() * v1.x() + tg2drst.y() * v2.x(),
        tg2drst.x() * v1.y() + tg2drst.y() * v2.y(),
        tg2drst.x() * v1.z() + tg2drst.y() * v2.z(),
    )
}
