//! `IntImpParGen` (`IntImpParGen.cxx:26-251`): domain normalisation and the
//! position / transition determination used by `IntImpParGen_Intersector`.
//!
//! `IntImpParGen.cxx` and `IntImpParGen_Tool.cxx` define functions with the same
//! names but different bodies (the former compares curvatures against
//! `TOLERANCE_ANGULAIRE` = 1e-8 and normalises the period with a two-condition
//! `while`, the latter against `gp::Resolution()` and a one-condition `while`).
//! Only the `IntImpParGen` bodies are needed by the intersector, so only they
//! are ported here; `IntImpParGen_Tool` stays UNPORTED until a caller needs it.

use crate::gp::{GpPnt2d, GpVec2d};
use crate::intres2d::{
    IntRes2dDomain, IntRes2dPosition, IntRes2dSituation, IntRes2dTransition, IntRes2dTypeTrans,
};

/// `TOLERANCE_ANGULAIRE` (`IntImpParGen.cxx:23`).
const TOLERANCE_ANGULAIRE: f64 = 0.00000001;
/// `DERIVEE_PREMIERE_NULLE` (`IntImpParGen.cxx:24`).
const DERIVEE_PREMIERE_NULLE: f64 = 0.000000000001;

/// `IntImpParGen::NormalizeOnDomain` (`IntImpParGen.cxx:28-45`). OCCT takes
/// `Param` by non-const reference but never writes it; the returned value is
/// the only effect.
pub fn normalize_on_domain(param: f64, the_domain: &IntRes2dDomain) -> f64 {
    let mut mod_param = param;
    if the_domain.is_closed() {
        let (t, p) = the_domain.equivalent_parameters();
        let periode = p - t;
        while mod_param < the_domain.first_parameter()
            && mod_param + periode < the_domain.last_parameter()
        {
            mod_param += periode;
        }
        while mod_param > the_domain.last_parameter()
            && mod_param - periode > the_domain.first_parameter()
        {
            mod_param -= periode;
        }
    }
    mod_param
}

/// `IntImpParGen::DeterminePosition` (`IntImpParGen.cxx:48-82`).
pub fn determine_position(
    the_domain: &IntRes2dDomain,
    pnt1: &GpPnt2d,
    param1: f64,
) -> IntRes2dPosition {
    let mut pos1 = IntRes2dPosition::Middle;

    if the_domain.has_first_point() {
        if pnt1.distance(the_domain.first_point()) <= the_domain.first_tolerance() {
            pos1 = IntRes2dPosition::Head;
        }
    }

    if the_domain.has_last_point() {
        if pnt1.distance(the_domain.last_point()) <= the_domain.last_tolerance() {
            if pos1 == IntRes2dPosition::Head {
                if (param1 - the_domain.last_parameter()).abs()
                    < (param1 - the_domain.first_parameter()).abs()
                {
                    pos1 = IntRes2dPosition::End;
                }
            } else {
                pos1 = IntRes2dPosition::End;
            }
        }
    }

    pos1
}

/// `IntImpParGen::DetermineTransition(Pos1, Tan1, Norm1, T1, Pos2, Tan2,
/// Norm2, T2, Tol)` (`IntImpParGen.cxx:85-205`). `Tan1` / `Tan2` may be
/// overwritten with the normal when the first derivative vanishes. The
/// tolerance parameter is unnamed in OCCT.
#[allow(clippy::too_many_arguments)]
pub fn determine_transition_touch(
    pos1: IntRes2dPosition,
    tan1: &mut GpVec2d,
    norm1: &GpVec2d,
    t1: &mut IntRes2dTransition,
    pos2: IntRes2dPosition,
    tan2: &mut GpVec2d,
    norm2: &GpVec2d,
    t2: &mut IntRes2dTransition,
    _tolerance: f64,
) {
    let mut courbure1 = true;
    let mut courbure2 = true;
    let mut decide = true;

    t1.set_position(pos1);
    t2.set_position(pos2);

    if tan1.square_magnitude() <= DERIVEE_PREMIERE_NULLE {
        *tan1 = *norm1;
        courbure1 = false;
        if tan1.square_magnitude() <= DERIVEE_PREMIERE_NULLE {
            // transition undecided
            decide = false;
        }
    }

    if tan2.square_magnitude() <= DERIVEE_PREMIERE_NULLE {
        *tan2 = *norm2;
        courbure2 = false;
        if tan2.square_magnitude() <= DERIVEE_PREMIERE_NULLE {
            // transition undecided
            decide = false;
        }
    }

    if !decide {
        t1.set_undecided(pos1);
        t2.set_undecided(pos2);
    } else {
        let sgn = tan1.crossed(tan2);
        let norm = tan1.magnitude() * tan2.magnitude();

        if sgn.abs() <= TOLERANCE_ANGULAIRE * norm {
            // Transition TOUCH
            let opos = tan1.dot(tan2) < 0.0;
            if !(courbure1 || courbure2) {
                t1.set_touch(true, pos1, IntRes2dSituation::Unknown, opos);
                t2.set_touch(true, pos2, IntRes2dSituation::Unknown, opos);
            } else {
                let norm_v = GpVec2d::new(-tan1.y(), tan1.x());
                let val1 = if !courbure1 { 0.0 } else { norm_v.dot(norm1) };
                let val2 = if !courbure2 { 0.0 } else { norm_v.dot(norm2) };

                if (val1 - val2).abs() <= TOLERANCE_ANGULAIRE {
                    t1.set_touch(true, pos1, IntRes2dSituation::Unknown, opos);
                    t2.set_touch(true, pos2, IntRes2dSituation::Unknown, opos);
                } else if val2 > val1 {
                    t2.set_touch(true, pos2, IntRes2dSituation::Inside, opos);
                    if opos {
                        t1.set_touch(true, pos1, IntRes2dSituation::Inside, opos);
                    } else {
                        t1.set_touch(true, pos1, IntRes2dSituation::Outside, opos);
                    }
                } else {
                    // Val1 > Val2
                    t2.set_touch(true, pos2, IntRes2dSituation::Outside, opos);
                    if opos {
                        t1.set_touch(true, pos1, IntRes2dSituation::Outside, opos);
                    } else {
                        t1.set_touch(true, pos1, IntRes2dSituation::Inside, opos);
                    }
                }
            }
        } else if sgn < 0.0 {
            t1.set_in_out(false, pos1, IntRes2dTypeTrans::In);
            t2.set_in_out(false, pos2, IntRes2dTypeTrans::Out);
        } else {
            // sgn > 0
            t1.set_in_out(false, pos1, IntRes2dTypeTrans::Out);
            t2.set_in_out(false, pos2, IntRes2dTypeTrans::In);
        }
    }
}

/// `IntImpParGen::DetermineTransition(Pos1, Tan1, T1, Pos2, Tan2, T2, Tol)`
/// (`IntImpParGen.cxx:208-251`). The `Tan1` / `Tan2` parameters are non-const
/// references in OCCT but are only read here. The tolerance is unnamed.
pub fn determine_transition_simple(
    pos1: IntRes2dPosition,
    tan1: &GpVec2d,
    t1: &mut IntRes2dTransition,
    pos2: IntRes2dPosition,
    tan2: &GpVec2d,
    t2: &mut IntRes2dTransition,
    _tolerance: f64,
) -> bool {
    t1.set_position(pos1);
    t2.set_position(pos2);

    let tan1_magnitude = tan1.magnitude();
    if tan1_magnitude <= DERIVEE_PREMIERE_NULLE {
        return false;
    }

    let tan2_magnitude = tan2.magnitude();
    if tan2_magnitude <= DERIVEE_PREMIERE_NULLE {
        return false;
    }

    let sgn = tan1.crossed(tan2);
    let norm = tan1_magnitude * tan2_magnitude;

    if sgn.abs() <= TOLERANCE_ANGULAIRE * norm {
        // Transition TOUCH
        return false;
    } else if sgn < 0.0 {
        t1.set_in_out(false, pos1, IntRes2dTypeTrans::In);
        t2.set_in_out(false, pos2, IntRes2dTypeTrans::Out);
    } else {
        // sgn > 0
        t1.set_in_out(false, pos1, IntRes2dTypeTrans::Out);
        t2.set_in_out(false, pos2, IntRes2dTypeTrans::In);
    }
    true
}
