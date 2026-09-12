//! Restriction-line branch of `GeomInt_LineConstructor::Perform`.

use occt_core::gp::GpPnt2d;
use occt_core::precision::PCONFUSION;

use crate::abs::Orientation;
use crate::fclass2d::FaceState;
use crate::geom_int::line_tool;
use crate::geom_int::param_ori::ParameterAndOrientation;
use crate::geom_int::types::GeomIntLine;
use crate::geom_int::LineConstructor;
use crate::int_tools_wline::TransType;

pub(crate) fn perform_restriction(ctor: &mut LineConstructor, line: &GeomIntLine, tol: f64) {
    ctor.set_done(false);
    ctor.seqp_mut().clear();
    let nbvtx = line_tool::nb_vertex(line);
    if nbvtx == 0 {
        ctor.seqp_mut().push(line_tool::first_parameter(line));
        ctor.seqp_mut().push(line_tool::last_parameter(line));
        ctor.set_done(true);
        return;
    }
    let mut seqpss: Vec<ParameterAndOrientation> = Vec::new();
    for i in 1..=nbvtx {
        let thevtx = line_tool::vertex(line, i);
        let prm = thevtx.parameter_on_line();
        let or1 = trans_to_ori(thevtx.on_dom_s1, thevtx.trans1);
        let or2 = trans_to_ori(thevtx.on_dom_s2, thevtx.trans2);
        let mut inserted = false;
        for j in 0..seqpss.len() {
            if (prm - seqpss[j].parameter()).abs() <= tol {
                accumulate(&mut seqpss[j], or1, or2);
                inserted = true;
                break;
            }
            if prm < seqpss[j].parameter() - tol {
                seqpss.insert(j, ParameterAndOrientation::with(prm, or1, or2));
                inserted = true;
                break;
            }
        }
        if !inserted {
            seqpss.push(ParameterAndOrientation::with(prm, or1, or2));
        }
    }

    let mut trim = false;
    let mut dans_s1 = false;
    let mut dans_s2 = false;
    let nb = seqpss.len();
    let mut i_found = nb + 1;
    for (i, item) in seqpss.iter().enumerate() {
        if item.orientation1() != Orientation::Internal {
            trim = true;
            dans_s1 = item.orientation1() != Orientation::Forward;
            i_found = i;
            break;
        }
    }
    if i_found > nb {
        let Some((d1, _)) = ctor.domain_clones() else {
            ctor.set_done(false);
            return;
        };
        for i in 1..=line_tool::nb_vertex(line) {
            let v = line_tool::vertex(line, i);
            if !v.on_dom_s1 {
                if d1.classify(GpPnt2d::new(v.u1, v.v1), tol) == FaceState::Out {
                    ctor.set_done(true);
                    return;
                }
                break;
            }
        }
        dans_s1 = true;
    }
    i_found = nb + 1;
    for (i, item) in seqpss.iter().enumerate() {
        if item.orientation2() != Orientation::Internal {
            trim = true;
            dans_s2 = item.orientation2() != Orientation::Forward;
            i_found = i;
            break;
        }
    }
    if i_found > nb {
        let Some((_, d2)) = ctor.domain_clones() else {
            ctor.set_done(false);
            return;
        };
        for i in 1..=line_tool::nb_vertex(line) {
            let v = line_tool::vertex(line, i);
            if !v.on_dom_s2 {
                if d2.classify(GpPnt2d::new(v.u2, v.v2), tol) == FaceState::Out {
                    ctor.set_done(true);
                    return;
                }
                break;
            }
        }
        dans_s2 = true;
    }
    if !trim {
        ctor.seqp_mut().push(line_tool::first_parameter(line));
        ctor.seqp_mut().push(line_tool::last_parameter(line));
        ctor.set_done(true);
        return;
    }
    let thefirst = line_tool::first_parameter(line);
    let thelast = line_tool::last_parameter(line);
    let mut firstp = thefirst;
    for item in &seqpss {
        let or1 = item.orientation1();
        let or2 = item.orientation2();
        if dans_s1 && dans_s2 {
            if or1 == Orientation::Reversed {
                dans_s1 = false;
            }
            if or2 == Orientation::Reversed {
                dans_s2 = false;
            }
            if !dans_s1 || !dans_s2 {
                let lastp = item.parameter();
                let stofirst = firstp.max(thefirst);
                let stolast = lastp.min(thelast);
                if stolast > stofirst {
                    ctor.seqp_mut().push(stofirst);
                    ctor.seqp_mut().push(stolast);
                }
                if lastp > thelast {
                    break;
                }
            }
        } else {
            if dans_s1 {
                if or1 == Orientation::Reversed {
                    dans_s1 = false;
                }
            } else if or1 == Orientation::Forward {
                dans_s1 = true;
            }
            if dans_s2 {
                if or2 == Orientation::Reversed {
                    dans_s2 = false;
                }
            } else if or2 == Orientation::Forward {
                dans_s2 = true;
            }
            if dans_s1 && dans_s2 {
                firstp = item.parameter();
            }
        }
    }
    if dans_s1 && dans_s2 {
        let lastp = thelast;
        firstp = firstp.max(thefirst);
        if lastp > firstp {
            ctor.seqp_mut().push(firstp);
            ctor.seqp_mut().push(lastp);
        }
    }
    let _ = PCONFUSION;
    ctor.set_done(true);
}

fn trans_to_ori(on_dom: bool, t: TransType) -> Orientation {
    if !on_dom {
        return Orientation::Internal;
    }
    match t {
        TransType::In => Orientation::Forward,
        TransType::Out => Orientation::Reversed,
        TransType::Touch | TransType::Undecided => Orientation::Internal,
    }
}

fn accumulate(valj: &mut ParameterAndOrientation, or1: Orientation, or2: Orientation) {
    if or1 != Orientation::Internal {
        if valj.orientation1() != Orientation::Internal {
            if or1 != valj.orientation1() {
                valj.set_orientation1(Orientation::Internal);
            }
        } else {
            valj.set_orientation1(or1);
        }
    }
    if or2 != Orientation::Internal {
        if valj.orientation2() != Orientation::Internal {
            if or2 != valj.orientation2() {
                valj.set_orientation2(Orientation::Internal);
            }
        } else {
            valj.set_orientation2(or2);
        }
    }
}
