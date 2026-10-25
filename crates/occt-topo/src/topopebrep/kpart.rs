//! `TopOpeBRep_kpart.cxx` — the VPoint keep/drop filters the curve/point
//! interferences go through.
//!
//! Source (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRep`):
//! - `FUNBREP_Periodize` (`TopOpeBRep_kpart.cxx:49-92`)
//! - `FUNBREP_HasSameGPoint` (`:94-115`)
//! - `FUNBREP_SameUV` (`:117-138`)
//! - `FUNBREP_topowalki` (`:261-458`)
//! - `FUNBREP_topogline` (`:609-790`)
//! - `TopoParameter` (`:796-828`)
//! - `FUNBREP_topokpart` (`:831-916`)
//!
//! `STATIC_lastVPind` (`cxx:41`) is a file-scope `int`: the index of the VPoint
//! the last CPI was computed on, `0` before any CPI. It is a `thread_local`
//! here, like the other OCCT file statics of the port (`ds/connex.rs`,
//! `ds/face_interference_tool.rs`).
//!
//! UNPORTED, in this file:
//! - `FUNBREP_topowalki_new` (`cxx:140-255`) and `FUNBREP_topogline_new`
//!   (`:460-607`) — both are reachable only through the `#ifdef OCCT_DEBUG`
//!   gate `TopOpeBRep_GetcontextNEWKP()` (`:285-293`, `:634-642`); the release
//!   control flow of `FUNBREP_topokpart` never calls them.
//! - `TopOpeBRep::KPart` and the rest of `TopOpeBRep_kpart.cxx`
//!   (`:918-2670`), which belong to `TopOpeBRepBuild`'s `Builder`.

use std::cell::Cell;

use occt_core::precision::PCONFUSION;

use crate::abs::{Orientation, ShapeType, State};
use crate::shape::TopoShape;
use crate::topopebrep::ds::data_structure::DataStructure;
use crate::topopebrep::ds::interference::InterferenceHandle;
use crate::topopebrep::ds::interference_tool::set_parameter;
use crate::topopebrep::ds::point::Point;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::ds::InterferenceList;
use crate::topopebrep::interf::point_geom_tool::PointGeomTool;
use crate::topopebrep::interf::v_point_inter::VPointInter;
use crate::topopebrep::line_inter::LineInter;
use crate::topopebrep::type_line_curve::TypeLineCurve;

thread_local! {
    /// `static int STATIC_lastVPind` (`cxx:41`).
    static LAST_VP_IND: Cell<i32> = const { Cell::new(0) };
}

/// `FUNBREP_Periodize(L, Ifound, PIfound, parline)` (`cxx:49-92`).
///
/// Rewrites the two parameters so that a VPoint sitting exactly on the start or
/// the end of a periodic line is described at the opposite end of the other
/// one, then stores the result on `Ifound`.
fn fun_brep_periodize(
    l: &mut LineInter,
    ifound: &InterferenceHandle,
    p_ifound: &mut f64,
    parline: &mut f64,
) -> Result<(), String> {
    let o_t_ifound = {
        let i = ifound.read().expect("interference lock");
        i.transition().orientation(State::In, ShapeType::Face)
    };
    if l.is_periodic() {
        let (f, last) = l.bounds();
        let onf = (*parline - f).abs() < PCONFUSION;
        let onl = (*parline - last).abs() < PCONFUSION;
        let onfl = onf || onl;
        if onfl {
            if o_t_ifound == Orientation::Forward {
                *p_ifound = f;
                *parline = last;
            } else if o_t_ifound == Orientation::Reversed {
                *p_ifound = last;
                *parline = f;
            } else {
                return Err("FUNBREP_Periodize".to_string());
            }
        } else {
            let m = parline.min(*p_ifound);
            *parline = m;
            *p_ifound = m;
        }
    } else {
        let m = parline.min(*p_ifound);
        *parline = m;
        *p_ifound = m;
    }
    Ok(())
}

/// `FUNBREP_HasSameGPoint(DSP, I, BDS)` (`cxx:94-115`): the DS point `DSP`
/// stands for the same geometric point as `I`'s geometry.
fn fun_brep_has_same_g_point(
    dsp: &Point,
    i: &crate::topopebrep::ds::interference::Interference,
    bds: &DataStructure,
) -> Result<bool, String> {
    use crate::topopebrep::ds::kind::Kind;
    match i.geometry_type() {
        Kind::Point => {
            let p = bds
                .point(i.geometry())
                .ok_or_else(|| "FUNBREP_HasSameGPoint : no point".to_string())?;
            Ok(dsp.is_equal(p))
        }
        Kind::Vertex => {
            let s = bds
                .shape(i.geometry(), true)
                .ok_or_else(|| "FUNBREP_HasSameGPoint : no vertex".to_string())?;
            let p = Point::of_shape(&s)?;
            Ok(dsp.is_equal(&p))
        }
        _ => Ok(false),
    }
}

/// `FUNBREP_SameUV(VP1, VP2, sind, toluv)` (`cxx:117-138`).
fn fun_brep_same_uv(
    vp1: &VPointInter,
    vp2: &VPointInter,
    sind: i32,
    toluv: f64,
) -> Result<bool, String> {
    let (p1, p2) = if sind == 1 {
        (vp1.parameters_on_s1()?, vp2.parameters_on_s1()?)
    } else if sind == 2 {
        (vp1.parameters_on_s2()?, vp2.parameters_on_s2()?)
    } else {
        ((0.0, 0.0), (0.0, 0.0))
    };
    let d = ((p1.0 - p2.0).powi(2) + (p1.1 - p2.1).powi(2)).sqrt();
    Ok(d < toluv)
}

/// `TopoParameter(L, I, parline, closingedge)` (`cxx:796-828`): the current
/// VPoint falls on `I`'s geometry at an equal parameter.
fn topo_parameter(
    l: &LineInter,
    i: &InterferenceHandle,
    parline: f64,
    closingedge: bool,
) -> bool {
    let p_cpi = {
        let g = i.read().expect("interference lock");
        crate::topopebrep::ds::export::fds_parameter(&g)
    };
    let Some(p_cpi) = p_cpi else {
        // `TopOpeBRepDS_InterferenceTool::Parameter` raises for the types
        // without one; `I` is always a curve/point interference here.
        return false;
    };
    if !closingedge {
        (parline - p_cpi).abs() < PCONFUSION
    } else if l.is_periodic() {
        // `cxx:812-821`: leave the (parameter, transition) pairs of the
        // coincident points on the sewing edge to `FUNBREP_topokpart`.
        true
    } else {
        (parline - p_cpi).abs() < PCONFUSION
    }
}

/// `FUNBREP_topowalki(Ifound, DSCIL, L, VP, lasttransLine, BDS, E, F, toluv,
/// CPIfound, samepar, couture, parline, transLine)` (`cxx:261-458`).
///
/// `BDS`, `F`, `toluv` and `parline` are unused by the C++ body (they are
/// commented out in the signature at `cxx:266-279`); the port does not take
/// them.
#[allow(clippy::too_many_arguments)]
fn fun_brep_topowalki(
    ifound: Option<&InterferenceHandle>,
    dscil: &InterferenceList,
    l: &LineInter,
    vp: &VPointInter,
    lasttrans_line: &Transition,
    e: &TopoShape,
    cpi_found: bool,
    samepar: bool,
    couture: bool,
    trans_line: &mut Transition,
) -> Result<bool, String> {
    if l.type_line_curve() != TypeLineCurve::Walking {
        return Err("FUNBREP_topowalki : line is not a walking".to_string());
    }

    let eori = e.orientation();
    let i_vp = vp.index();
    let (i_inon1, i_inonn, _n_inon) = l
        .clone()
        .vp_bounds()
        .ok_or_else(|| "FUNBREP_topowalki : VPBounds".to_string())?;
    let fermee = l.is_v_closed();
    let mut onsort = trans_line.orientation(State::Out, ShapeType::Face) == Orientation::Forward;
    let lastdefinie = !lasttrans_line.is_unknown();
    let mut lastonsort = false;
    let mut lastinin = false;

    if lastdefinie {
        lastonsort = lasttrans_line.orientation(State::Out, ShapeType::Face) == Orientation::Forward;
        lastinin = lasttrans_line.before_state() == State::In;
        lastinin = lastinin && (lasttrans_line.after_state() == State::In);
    }
    let _ = lastinin; // computed at `cxx:307-311`, never read.

    let mut litdefinie = false;
    let mut litonsort = false;
    let mut lit = Transition::default();
    let nointerf = dscil.is_empty();
    if !nointerf {
        let i = dscil.last().expect("non-empty DSCIL").clone();
        let g = i.read().expect("interference lock");
        lit = g.transition();
        litdefinie = !lit.is_unknown();
        if litdefinie {
            litonsort = lit.orientation(State::Out, ShapeType::Face) == Orientation::Forward;
        }
    }

    let mut keep = true;
    let efr = eori == Orientation::Forward || eori == Orientation::Reversed;

    if !litdefinie {
        if onsort {
            let mut cond_cts19305 = !cpi_found;
            cond_cts19305 = cond_cts19305 && !fermee;
            cond_cts19305 = cond_cts19305 && !lastdefinie;
            cond_cts19305 = cond_cts19305 && (i_vp == 1) && (i_vp == i_inon1);
            if cond_cts19305 {
                keep = true;
                *trans_line = Transition::new(State::Out, State::In, ShapeType::Face, ShapeType::Face);
            } else {
                keep = false;
            }
        }
        if i_vp == i_inonn {
            keep = false;
        }
        if !efr && i_vp == i_inon1 {
            *trans_line = Transition::new(State::Out, State::In, ShapeType::Face, ShapeType::Face);
        }
    }

    if keep && litdefinie {
        let lprec_is_entrant = !litonsort;
        let entrantsortant = lprec_is_entrant && onsort;

        if !efr {
            if i_vp == i_inonn {
                *trans_line =
                    Transition::new(State::In, State::Out, ShapeType::Face, ShapeType::Face);
            }
        } else {
            let mut cond_cts19305 = !cpi_found;
            cond_cts19305 = cond_cts19305 && !fermee;
            cond_cts19305 = cond_cts19305 && (lastdefinie && !lastonsort);
            cond_cts19305 = cond_cts19305 && (litdefinie && !litonsort);
            cond_cts19305 = cond_cts19305 && (i_vp == i_inonn);

            if cond_cts19305 {
                keep = true;
                *trans_line =
                    Transition::new(State::In, State::Out, ShapeType::Face, ShapeType::Face);
            } else if !entrantsortant {
                keep = false;
            }
        }

        let samevponcouture = samepar && couture;
        if keep && samevponcouture {
            keep = false;
        }

        if keep && cpi_found && samepar {
            let p_ifound = ifound
                .map(|i| {
                    crate::topopebrep::ds::export::fds_parameter_or_fail(
                        &i.read().expect("interference lock"),
                    )
                })
                .transpose()?;
            let Some(p_ifound) = p_ifound else {
                return Err("FUNBREP_topowalki : null Ifound".to_string());
            };

            let mut cond1 = lastdefinie && !lastonsort;
            cond1 = cond1 && onsort;
            cond1 = cond1 && fermee;
            cond1 = cond1 && (p_ifound == 1.0);
            cond1 = cond1 && (i_vp == i_inon1 + 1) && (i_vp != i_inonn);
            if cond1 {
                return Ok(false);
            }

            let mut cond2 = lastdefinie && !lastonsort;
            // `cxx:439` writes `cond1`, not `cond2`; the assignment is dead
            // because `cond2` is rebuilt from `cond2` below. Kept verbatim.
            cond2 = cond2 && !fermee;
            cond2 = cond2 && (p_ifound == 1.0);
            cond2 = cond2 && (i_vp == i_inon1 + 1) && (i_vp != i_inonn);
            if cond2 {
                return Ok(false);
            }
        }
    }

    if keep {
        LAST_VP_IND.with(|c| c.set(i_vp));
    }
    Ok(keep)
}

/// `FUNBREP_topogline(Ifound, DSCIL, L, VP, BDS, E, F, toluv, CPIfound,
/// samepar, couture, parline, transLine)` (`cxx:609-790`).
#[allow(clippy::too_many_arguments)]
fn fun_brep_topogline(
    ifound: Option<&InterferenceHandle>,
    dscil: &InterferenceList,
    l: &LineInter,
    vp: &VPointInter,
    bds: &DataStructure,
    e: &TopoShape,
    toluv: f64,
    samepar: bool,
    couture: bool,
    parline: &mut f64,
    trans_line: &mut Transition,
) -> Result<bool, String> {
    if l.type_line_curve() == TypeLineCurve::Walking {
        return Err("FUNBREP_topogline : line is not a GLine".to_string());
    }

    let eori = e.orientation();
    let efr = eori == Orientation::Forward || eori == Orientation::Reversed;
    let i_vp = vp.index();
    let (i_inon1, i_inonn, _n_inon) = l
        .clone()
        .vp_bounds()
        .ok_or_else(|| "FUNBREP_topogline : VPBounds".to_string())?;
    let fermee = l.is_v_closed();

    if !efr {
        if i_vp == i_inon1 {
            *trans_line = Transition::new(State::Out, State::In, ShapeType::Face, ShapeType::Face);
        }
        if i_vp == i_inonn {
            *trans_line = Transition::new(State::In, State::Out, ShapeType::Face, ShapeType::Face);
        }
    }
    let mut onsort = trans_line.orientation(State::Out, ShapeType::Face) == Orientation::Forward;

    let mut litdefinie = false;
    let mut litonsort = false;
    let mut lit = Transition::default();
    let dscilempty = dscil.is_empty();

    // `cxx:677-683`: `VP<iVP>` on 3 is described once, on the first `VP`.
    let notkept = !dscilempty && (i_vp == 1);
    if notkept {
        return Ok(false);
    }

    if !dscilempty {
        let i = dscil.last().expect("non-empty DSCIL").clone();
        let g = i.read().expect("interference lock");
        lit = g.transition();
        litdefinie = !lit.is_unknown();
        if litdefinie {
            litonsort = lit.orientation(State::Out, ShapeType::Face) == Orientation::Forward;
        }
    }

    // `cxx:696-718`: `VPmin` with an unknown transition becomes `OUT/IN`,
    // `VPmax` with an unknown transition becomes the complement of the last
    // known one.
    let trliunk = trans_line.is_unknown();
    let isfirstvp = i_vp == i_inon1;
    let islastvp = i_vp == i_inonn;
    if trliunk {
        if isfirstvp {
            *trans_line =
                Transition::new(State::Out, State::In, ShapeType::Face, ShapeType::Face);
            onsort = false;
        }
        if islastvp && litdefinie && !litonsort {
            *trans_line = lit.complement();
            onsort = true;
        }
    }

    let hasfp = l.has_first_point()?;
    let haslp = l.has_last_point()?;
    let hasfol = hasfp || haslp;

    let mut keep = true;

    if !litdefinie {
        if onsort {
            keep = false;
        }
        if i_vp == i_inonn {
            keep = false;
        }
    }

    if litdefinie {
        let lprec_is_entrant = !litonsort;
        let entrantsortant = lprec_is_entrant && onsort;
        if !entrantsortant {
            keep = false;
        }

        let sameparoncouture = samepar && couture;
        if sameparoncouture && hasfol {
            if keep {
                let pnt_vp = PointGeomTool::make_point_of_v_point_inter(vp)?;
                let last_index = LAST_VP_IND.with(|c| c.get());
                let last = dscil.last().expect("non-empty DSCIL").clone();
                let mut samegp = {
                    let g = last.read().expect("interference lock");
                    fun_brep_has_same_g_point(&pnt_vp, &g, bds)?
                };
                // `cxx:756`: `L.VPoint(STATIC_lastVPind)`.
                let last_vp = l
                    .v_point(last_index)
                    .ok_or_else(|| "FUNBREP_topogline : VPoint".to_string())?;
                if fermee {
                    // Two VPoints describing a closing line share the 3d point
                    // and the parameter but not the uv point (cobo121).
                    let same_uv_on1 = fun_brep_same_uv(vp, last_vp, 1, toluv)?;
                    let same_uv_on2 = fun_brep_same_uv(vp, last_vp, 2, toluv)?;
                    samegp = samegp && same_uv_on1 && same_uv_on2;
                }
                if samegp {
                    keep = false;
                }
            }
        }
        if sameparoncouture && !hasfol {
            let Some(ifound) = ifound else {
                return Err("FUNBREP_topogline : null Ifound".to_string());
            };
            let mut p_ifound = {
                let i = ifound.read().expect("interference lock");
                crate::topopebrep::ds::export::fds_parameter_or_fail(&i)?
            };
            let mut lc = l.clone();
            fun_brep_periodize(&mut lc, ifound, &mut p_ifound, parline)?;
            ifound
                .write()
                .expect("interference lock")
                .set_parameter(p_ifound);
            *trans_line = lit.complement();
        }
    }

    if keep {
        LAST_VP_IND.with(|c| c.set(i_vp));
    }
    Ok(keep)
}

/// `FUNBREP_topokpart(Ifound, DSCIL, L, VP, BDS, E, F, toluv, parline,
/// transLine)` (`cxx:831-916`): whether the CPI of the current VPoint must be
/// kept, and the line transition to store with it.
///
/// `parline` is an in/out parameter (it is rewritten by `FUNBREP_Periodize`).
#[allow(clippy::too_many_arguments)]
pub fn fun_brep_topokpart(
    ifound: Option<&InterferenceHandle>,
    dscil: &InterferenceList,
    l: &mut LineInter,
    vp: &VPointInter,
    bds: &DataStructure,
    e: &TopoShape,
    f: &TopoShape,
    toluv: f64,
    parline: &mut f64,
    trans_line: &mut Transition,
) -> Result<bool, String> {
    let cpi_found = ifound.is_some();
    let couture = crate::topopebrep::tool::shape_tool::closed(e, f);
    let samepar = if !cpi_found {
        false
    } else {
        topo_parameter(l, ifound.expect("CPI found"), *parline, couture)
    };

    let mut lasttrans_line = Transition::default();
    if let Some(i) = dscil.last() {
        lasttrans_line = i.read().expect("interference lock").transition();
    }

    if l.type_line_curve() == TypeLineCurve::Walking {
        fun_brep_topowalki(
            ifound,
            dscil,
            l,
            vp,
            &lasttrans_line,
            e,
            cpi_found,
            samepar,
            couture,
            trans_line,
        )
    } else {
        fun_brep_topogline(
            ifound, dscil, l, vp, bds, e, toluv, samepar, couture, parline, trans_line,
        )
    }
}

/// `TopOpeBRepDS_InterferenceTool::Parameter(I, Par)` (`cxx:779`), exposed so
/// `FacesFiller` does not need the `ds::interference_tool` path.
pub fn set_interference_parameter(i: &InterferenceHandle, par: f64) {
    set_parameter(&mut i.write().expect("interference lock"), par);
}
