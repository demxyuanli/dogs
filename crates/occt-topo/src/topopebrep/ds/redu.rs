//! `TopOpeBRepDS_redu` — the 1d/2d interference reduction helpers.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_redu.cxx:35-353`.

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};
use crate::topopebrep::ds::data_structure::DataStructure;
use crate::topopebrep::ds::export::{fds_assign, fds_copy, fds_idata, fds_parameter};
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::{
    fun_select_sk_interference, fun_select_traori_interference, fun_select_trasha_interference,
};
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::tool::closed::is_closing_e_in_shape;
use crate::topopebrep::tool::make_transition::MakeTransition;
use crate::topopebrep::tool::project_curve::fun_tool_par_e_tol;
use crate::topopebrep::tool::topology::{fun_tool_in_s, fun_tool_par_on_e_f_tol};

use super::complete_ds::fun_ds_has_i2d;

/// The eight out-parameters of `FUN_scanloi` (`redu.cxx:35-54`).
pub struct ScanLoi {
    pub l_for: InterferenceList,
    pub n_for: i32,
    pub l_rev: InterferenceList,
    pub n_rev: i32,
    pub l_int: InterferenceList,
    pub n_int: i32,
    pub l_ext: InterferenceList,
    pub n_ext: i32,
}

/// `FUN_scanloi(lI, lFOR, FOR, lREV, REV, lINT, INT, lEXT, EXT)`
/// (`redu.cxx:35-54`): dispatch `lI` by transition orientation.
pub fn fun_scanloi(li: &InterferenceList) -> ScanLoi {
    let mut l_ext: InterferenceList = Vec::new();
    let mut l_for: InterferenceList = Vec::new();
    let mut l_rev: InterferenceList = Vec::new();
    let mut l_int: InterferenceList = Vec::new();
    fds_assign(li, &mut l_ext);
    let n_for = fun_select_traori_interference(&mut l_ext, Orientation::Forward, &mut l_for);
    let n_rev = fun_select_traori_interference(&mut l_ext, Orientation::Reversed, &mut l_rev);
    let n_int = fun_select_traori_interference(&mut l_ext, Orientation::Internal, &mut l_int);
    let n_ext = l_ext.len() as i32;
    ScanLoi {
        l_for,
        n_for,
        l_rev,
        n_rev,
        l_int,
        n_int,
        l_ext,
        n_ext,
    }
}

/// `FUN_ds_redu2d1d(BDS, ISE, I2d, l1d, newT2d)` (`redu.cxx:56-179`).
///
/// `None` is the OCCT `false` return (including the `NYIRAISE` of
/// `redu.cxx:157-158`).
pub fn fun_ds_redu2d1d(
    bds: &DataStructure,
    ise: i32,
    i2d: &crate::topopebrep::ds::interference::InterferenceHandle,
    l1d: &InterferenceList,
) -> Option<Transition> {
    // `redu.cxx:68-81`
    let (sb2, ib2, _sa2, ia2, gt2, _g2, _st2, s2) = fds_idata(&i2d.read().unwrap());
    let t2d = i2d.read().unwrap().transition();
    let o2 = t2d.orientation(crate::abs::State::In, ShapeType::Face);
    let mut new_t2d = Transition::default();
    new_t2d.set_index(ib2);
    new_t2d.set_orientation(o2);
    let ok2 = ib2 == ia2 && sb2 == ShapeType::Face && gt2 == Kind::Vertex;
    if !ok2 {
        return None;
    }

    // `redu.cxx:83-90`
    let se = bds.shape(ise, true)?;
    let f = bds.shape(ib2, true)?;
    let tol_f = BRepTool::face_tolerance(&Face(f.clone())) * 1.0e2; // nyitol
    let e = bds.shape(s2, true)?;
    let tol_e = BRepTool::edge_tolerance(&Edge(e.clone())) * 1.0e2; // nyitol

    if is_closing_e_in_shape(&e, &f, &f) {
        // `redu.cxx:91-128`
        let stb = t2d.before_state();
        let sta = t2d.after_state();
        if stb != sta {
            let pbef = 0.0;
            let paft = 0.0;
            let factor = 1.0e-4;
            let par_se = fds_parameter(&i2d.read().unwrap())?;
            let par_e = fun_tool_par_e_tol(&se, par_se, &e, tol_e)?;
            let uv = fun_tool_par_on_e_f_tol(&e, par_e, &f, tol_f)?;
            let mut mkt = MakeTransition::new();
            let mut stb1 = crate::abs::State::Unknown;
            let mut sta1 = crate::abs::State::Unknown;
            let mut ok = mkt.initialize(&se, pbef, paft, par_se, &f, &uv, factor);
            if ok {
                ok = mkt.set_rest(&e, par_e);
            }
            if ok {
                ok = mkt.mk_ton_e(&mut stb1, &mut sta1);
            }
            if ok {
                new_t2d.set_state_before(stb1);
                new_t2d.set_state_after(sta1);
            }
            return if ok { Some(new_t2d) } else { None };
        }
        return None;
    }

    // `redu.cxx:130-178`
    let mut before_in_1d = false;
    let mut after_in_1d = false;
    for i1d in l1d {
        let (sb1, ib1, sa1, ia1, _gt1, _g1, _st1, _s1) = fds_idata(&i1d.read().unwrap());
        let _ = (sb1, sa1);
        if ib1 != ia1 {
            continue;
        }
        let o1 = i1d
            .read()
            .unwrap()
            .transition()
            .orientation(crate::abs::State::In, ShapeType::Face);

        let esd = bds.shape(ib1, true)?;
        if !fun_tool_in_s(&esd, &f) {
            continue;
        }

        let b_in = o1 == Orientation::Internal || o1 == Orientation::Reversed;
        let a_in = o1 == Orientation::Internal || o1 == Orientation::Forward;
        if b_in && a_in {
            // NYIRAISE : I1d INTERNAL -> NO I2d
            return None;
        }
        if b_in {
            before_in_1d = true;
        }
        if a_in {
            after_in_1d = true;
        }
    }

    if before_in_1d {
        new_t2d.set_state_before(crate::abs::State::In);
    }
    if after_in_1d {
        new_t2d.set_state_after(crate::abs::State::In);
    }
    Some(new_t2d)
}

/// The six out-parameters of `FUN_ds_GetTr` (`redu.cxx:181-353`).
pub struct GetTr {
    pub stb: crate::abs::State,
    pub isb: i32,
    pub bdim: i32,
    pub sta: crate::abs::State,
    pub isa: i32,
    pub adim: i32,
}

/// `FUN_ds_GetTr(BDS, ISE, G, LIG, stb, isb, bdim, sta, isa, adim)`
/// (`redu.cxx:181-353`): the state of the edge `ISE` on both sides of `G`,
/// assuming `IN1d > I2d > I3d` and `OUT3d > OUT2d > OUT1d`.
///
/// OCCT returns `true` unconditionally (`redu.cxx:352`).
pub fn fun_ds_get_tr(ise: i32, lig: &InterferenceList) -> GetTr {
    // `redu.cxx:199-209`
    let mut ligcopy: InterferenceList = Vec::new();
    fds_copy(lig, &mut ligcopy);
    let mut l3d: InterferenceList = Vec::new();
    fds_assign(lig, &mut ligcopy);
    fun_select_sk_interference(&mut ligcopy, Kind::Face, &mut l3d);
    fds_assign(lig, &mut ligcopy);
    let l2d = fun_ds_has_i2d(&ligcopy);
    let mut l1d: InterferenceList = Vec::new();
    fds_assign(lig, &mut ligcopy);
    fun_select_trasha_interference(&mut ligcopy, ShapeType::Edge, &mut l1d);

    // `redu.cxx:211-232`
    let s1 = fun_scanloi(&l1d);
    let s2 = fun_scanloi(&l2d);
    let s3 = fun_scanloi(&l3d);

    let before_in_1d = s1.n_rev + s1.n_int > 0;
    let before_ou_1d = s1.n_for + s1.n_ext != 0;
    let after_in_1d = s1.n_for + s1.n_int > 0;
    let after_ou_1d = s1.n_rev + s1.n_ext != 0;
    let before_in_2d = s2.n_rev + s2.n_int > 0;
    let before_ou_2d = s2.n_for + s2.n_ext != 0;
    let after_in_2d = s2.n_for + s2.n_int > 0;
    let after_ou_2d = s2.n_rev + s2.n_ext != 0;
    let before_in_3d = s3.n_rev + s3.n_int > 0;
    let before_ou_3d = s3.n_for + s3.n_ext != 0;
    let after_in_3d = s3.n_for + s3.n_int > 0;
    let after_ou_3d = s3.n_rev + s3.n_ext != 0;

    let _ = ise;

    // `redu.cxx:234-292` : state before
    let mut stb = crate::abs::State::Unknown;
    let mut isb = 0;
    let mut bdim = 0;
    if before_in_1d {
        stb = crate::abs::State::In;
        bdim = 1;
        isb = first_index_before(&s1.l_rev, &s1.l_int);
    } else if before_in_2d {
        stb = crate::abs::State::In;
        bdim = 2;
        isb = first_index_before(&s2.l_rev, &s2.l_int);
    } else if before_in_3d {
        stb = crate::abs::State::In;
        bdim = 3;
        isb = first_index_before(&s3.l_rev, &s3.l_int);
    } else if before_ou_3d {
        stb = crate::abs::State::Out;
        bdim = 3;
        isb = first_index_before(&s3.l_for, &s3.l_ext);
    } else if before_ou_2d {
        stb = crate::abs::State::Out;
        bdim = 2;
        isb = first_index_before(&s2.l_for, &s2.l_ext);
    } else if before_ou_1d {
        stb = crate::abs::State::Out;
        bdim = 1;
        isb = first_index_before(&s1.l_for, &s1.l_ext);
    }

    // `redu.cxx:294-351` : state after
    let mut sta = crate::abs::State::Unknown;
    let mut isa = 0;
    let mut adim = 0;
    if after_in_1d {
        sta = crate::abs::State::In;
        adim = 1;
        isa = first_index_after(&s1.l_for, &s1.l_int);
    } else if after_in_2d {
        sta = crate::abs::State::In;
        adim = 2;
        isa = first_index_after(&s2.l_for, &s2.l_int);
    } else if after_in_3d {
        sta = crate::abs::State::In;
        adim = 3;
        isa = first_index_after(&s3.l_for, &s3.l_int);
    } else if after_ou_3d {
        sta = crate::abs::State::Out;
        adim = 3;
        isa = first_index_after(&s3.l_rev, &s3.l_ext);
    } else if after_ou_2d {
        sta = crate::abs::State::Out;
        adim = 2;
        isa = first_index_after(&s2.l_rev, &s2.l_ext);
    } else if after_ou_1d {
        sta = crate::abs::State::Out;
        adim = 1;
        isa = first_index_after(&s1.l_rev, &s1.l_ext);
    }

    GetTr {
        stb,
        isb,
        bdim,
        sta,
        isa,
        adim,
    }
}

/// `FDS_copy(a, l) ; FDS_copy(b, l) ; l.First()->Transition().IndexBefore()`
/// (`redu.cxx:244-248` and siblings): the two lists are concatenated and the
/// first element answers. `l1dREV` is filled before `l1dINT` in every arm.
fn first_index_before(a: &InterferenceList, b: &InterferenceList) -> i32 {
    if let Some(h) = a.first().or_else(|| b.first()) {
        return h.read().unwrap().transition().index_before();
    }
    0
}

/// Same as [`first_index_before`] for `IndexAfter()` (`redu.cxx:300-350`).
fn first_index_after(a: &InterferenceList, b: &InterferenceList) -> i32 {
    if let Some(h) = a.first().or_else(|| b.first()) {
        return h.read().unwrap().transition().index_after();
    }
    0
}
