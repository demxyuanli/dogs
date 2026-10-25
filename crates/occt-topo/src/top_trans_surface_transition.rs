//! `TopTrans_SurfaceTransition` — transition of a surface of revolution
//! crossing a surfacic boundary along an intersection curve.
//! Source: `TKGeomAlgo/TopTrans/TopTrans_SurfaceTransition.cxx` (OCCT 8.0.0,
//! 672 lines) with the declaration in `TopTrans_SurfaceTransition.hxx:49-139`.
//!
//! The algorithm keeps a 2x2 table indexed by `(i, j)` where `i` encodes the
//! sign of `cos` and `j` the sign of `sin` of the angle between `beafter`
//! (`myNorm ^ myTgt`) and the normal of the added face inside its matter
//! (`FUN_nCinsideS`, `cxx:22-35`):
//!
//! ```text
//! i = 0,1,2 : cos = 0, > 0, < 0
//! j = 0,1,2 : sin = 0, > 0, < 0
//! ```
//!
//! `myAng` is initialised to `100.` (the "undefined" marker, `cxx:339-345`,
//! `:359-365`) and compared against `Angref == 100.` in `FUN_refnearest`
//! (`cxx:126`, `:176`). The port keeps that sentinel as the literal
//! [`ANG_UNDEF`] constant rather than `NaN`, so the equality test stays exact
//! like the C++ one.
//!
//! `M_Unknown` has two definitions in OCCT (`cxx:112-116`): `(-100)` under
//! `OCCT_DEBUG`, `(-100.)` otherwise; both compare equal to an `int`, so the
//! port uses the single constant `-100`.

use occt_core::gp::GpDir;
use occt_core::precision::{ANGULAR, CONFUSION};

use crate::abs::{Orientation, State};

/// The `100.` "not yet defined" marker of `myAng` (`cxx:341`, `:362`,
/// `:560-561`).
const ANG_UNDEF: f64 = 100.0;

/// `M_Unknown (-100)` (`cxx:112-116`).
const M_UNKNOWN: i32 = -100;
/// `M_noupdate (0)` (`cxx:117`).
const M_NOUPDATE: i32 = 0;
/// `M_updateREF (1)` (`cxx:118`).
const M_UPDATE_REF: i32 = 1;
/// `M_Ointernal (10)` (`cxx:119`).
const M_OINTERNAL: i32 = 10;

/// `BEFORE (2)` / `AFTER (1)` (`cxx:551-552`).
const BEFORE: i32 = 2;
const AFTER: i32 = 1;

/// `FUN_nCinsideS(tgC, ngS)` (`cxx:22-35`): `ngS ^ tgC`, the normal to the
/// curve `C` that is tangent to `S` and oriented inside `S`.
fn n_c_inside_s(tg_c: &GpDir, ng_s: &GpDir) -> GpDir {
    ng_s.crossed(tg_c).unwrap_or_else(|_| GpDir::default_dir())
}

/// `FUN_OO(i)` (`cxx:40-53`): `1 <-> 2`, anything else `0`.
fn oo(i: i32) -> i32 {
    if i == 1 {
        return 2;
    }
    if i == 2 {
        return 1;
    }
    0
}

/// `FUN_Ang(Normref, beafter, TgC, Norm, O)` (`cxx:55-68`). The first
/// parameter is unused in OCCT (`cxx:55`, commented out).
fn fun_ang(beafter: &GpDir, tg_c: &GpDir, norm: &GpDir, o: Orientation) -> f64 {
    let mut diron_f = n_c_inside_s(tg_c, norm);
    if o == Orientation::Reversed {
        diron_f.reverse();
    }
    beafter.angle_with_ref(&diron_f, tg_c)
}

/// `FUN_getSTA(Ang, tola, i, j)` (`cxx:70-93`): the `(cos, sin)` sign code.
fn fun_get_sta(ang: f64, tola: f64) -> (i32, i32) {
    let cos = ang.cos();
    let sin = ang.sin();
    let nullcos = cos.abs() < tola;
    let nullsin = sin.abs() < tola;
    let i = if nullcos {
        0
    } else if cos > 0.0 {
        1
    } else {
        2
    };
    let j = if nullsin {
        0
    } else if sin > 0.0 {
        1
    } else {
        2
    };
    (i, j)
}

/// `FUN_refnearest(Angref, Oriref, Ang, Ori, tola)` (`cxx:120-151`) — the
/// plane/cylinder overload.
fn ref_nearest_plain(
    ang_ref: f64,
    ori_ref: Orientation,
    ang: f64,
    ori: Orientation,
    tola: f64,
) -> i32 {
    if ang_ref == ANG_UNDEF {
        return M_UPDATE_REF;
    }
    let cos_ref = ang_ref.cos();
    let cos = ang.cos();
    let dcos = cos_ref.abs() - cos.abs();
    if dcos.abs() < tola {
        // Tangent case: two boundary faces with the same side and tangent
        // normals. Opposite orientations give the INTERNAL complex
        // transition; EXTERNAL refers to no logical case (`cxx:136-148`).
        if ori.complement() == ori_ref {
            return M_OINTERNAL;
        }
        return M_UNKNOWN;
    }
    if dcos > 0.0 {
        M_NOUPDATE
    } else {
        M_UPDATE_REF
    }
}

/// `FUN_refnearest(i, j, CurvSref, Angref, Oriref, Curvref, Ang, Ori, Curv,
/// tola, TouchFlag)` (`cxx:153-271`) — the curvilinear overload.
#[allow(clippy::too_many_arguments)]
fn ref_nearest_full(
    i: i32,
    j: i32,
    curv_sref: f64,
    ang_ref: f64,
    ori_ref: Orientation,
    curv_ref: f64,
    ang: f64,
    ori: Orientation,
    curv: f64,
    tola: f64,
    touch_flag: &mut bool,
) -> i32 {
    let iisj = i == j;
    let abscos = ang.cos().abs();
    let i0 = (1.0 - abscos).abs() < tola;
    let j0 = abscos < tola;
    let nullcurv = curv == 0.0;
    let curvpos = curv > tola;
    let curvneg = curv < -tola;
    let nullcsref = curv_sref == 0.0;

    if ang_ref == ANG_UNDEF {
        if i0 {
            if iisj && curvneg {
                return M_NOUPDATE;
            }
            if !iisj && curvpos {
                return M_NOUPDATE;
            }
        }
        if j0 {
            if !nullcsref && j == 1 && iisj && (curvpos || nullcurv) {
                return M_UPDATE_REF;
            }
            if !nullcsref && j == 1 && !iisj && (curvneg || nullcurv) {
                return M_UPDATE_REF;
            }
        }
        if iisj && curvpos {
            return M_NOUPDATE;
        }
        if !iisj && curvneg {
            return M_NOUPDATE;
        }
        return M_UPDATE_REF;
    }

    let cos_ref = ang_ref.cos();
    let cos = ang.cos();
    let dcos = cos_ref.abs() - cos.abs();
    let samecos = dcos.abs() < tola;
    if samecos {
        // Analysis for tangent cases: two boundary faces with the same side
        // and a small `dironF` (`cxx:217-267`).
        if (curv_ref - curv).abs() < 1e-4 {
            if ori.complement() == ori_ref {
                return M_OINTERNAL;
            }
            return M_UNKNOWN;
        }

        let mut noupdate = false;
        if iisj && curv_ref > curv {
            noupdate = true;
        }
        if !iisj && curv_ref < curv {
            noupdate = true;
        }
        let mut updateref = if noupdate { M_NOUPDATE } else { M_UPDATE_REF };
        if !j0 {
            return updateref;
        }

        if !noupdate && !nullcsref {
            // (j==1) the face is ABOVE Sref, (j==2) BELOW (`cxx:250-260`).
            if j == 2 && curv.abs() < curv_sref {
                updateref = M_NOUPDATE;
            }
            if j == 1 && curv.abs() > curv_sref {
                updateref = M_NOUPDATE;
            }
        }
        return updateref;
    }

    let updateref = if dcos > 0.0 { M_NOUPDATE } else { M_UPDATE_REF };
    if ori_ref != ori {
        *touch_flag = true; // eap Mar 25 2002
    }
    updateref
}

/// The `NCollection_Array2<double> myAng` / `myCurv` pair plus
/// `NCollection_Array2<TopAbs_Orientation> myOri`, all bounded `(1,2,1,2)`
/// (`cxx:279-281`).
#[derive(Clone, Debug)]
struct Table2 {
    ang: [[f64; 2]; 2],
    curv: [[f64; 2]; 2],
    ori: [[Orientation; 2]; 2],
}

impl Table2 {
    fn new() -> Self {
        Self {
            ang: [[ANG_UNDEF; 2]; 2],
            curv: [[0.0; 2]; 2],
            ori: [[Orientation::Forward; 2]; 2],
        }
    }

    /// OCCT's `myAng(i, j)` with the `NCollection_Array2` lower bound `1`.
    fn ang(&self, i: i32, j: i32) -> f64 {
        self.ang[(i - 1) as usize][(j - 1) as usize]
    }

    fn set_ang(&mut self, i: i32, j: i32, v: f64) {
        self.ang[(i - 1) as usize][(j - 1) as usize] = v;
    }

    fn curv(&self, i: i32, j: i32) -> f64 {
        self.curv[(i - 1) as usize][(j - 1) as usize]
    }

    fn set_curv(&mut self, i: i32, j: i32, v: f64) {
        self.curv[(i - 1) as usize][(j - 1) as usize] = v;
    }

    fn ori(&self, i: i32, j: i32) -> Orientation {
        self.ori[(i - 1) as usize][(j - 1) as usize]
    }

    fn set_ori(&mut self, i: i32, j: i32, v: Orientation) {
        self.ori[(i - 1) as usize][(j - 1) as usize] = v;
    }

    /// The `for (i = 1; i <= 2; i++) for (j = 1; j <= 2; j++) myAng(i, j) =
    /// 100.` reset of both `Reset` overloads (`cxx:339-345`, `:359-365`).
    /// `myCurv` and `myOri` are *not* touched by OCCT, and are not here.
    fn reset_ang(&mut self) {
        for i in 1..=2 {
            for j in 1..=2 {
                self.set_ang(i, j, ANG_UNDEF);
            }
        }
    }
}

/// `TopTrans_SurfaceTransition` (`TopTrans_SurfaceTransition.hxx:49-139`).
#[derive(Clone, Debug)]
pub struct SurfaceTransition {
    tgt: GpDir,
    norm: GpDir,
    beafter: GpDir,
    curv_ref: f64,
    table: Table2,
    is_defined: bool,
    touch_flag: bool,
}

impl SurfaceTransition {
    /// `TopTrans_SurfaceTransition()` (`cxx:277-285`): `myCurvRef = 0`,
    /// `myIsDefined = myTouchFlag = false`. The three arrays are allocated
    /// `(1,2,1,2)` but never filled by the constructor.
    pub fn new() -> Self {
        Self {
            tgt: GpDir::default_dir(),
            norm: GpDir::default_dir(),
            beafter: GpDir::default_dir(),
            curv_ref: 0.0,
            table: Table2::new(),
            is_defined: false,
            touch_flag: false,
        }
    }

    /// `Reset(Tgt, Norm, MaxD, MinD, MaxCurv, MinCurv)`
    /// (`cxx:287-349`).
    ///
    /// The four early exits are reproduced as in OCCT: a plane reference
    /// (`|MaxCurv|, |MinCurv| < tola`, `cxx:306-309`) delegates to
    /// [`SurfaceTransition::reset_line`]; when neither principal direction is
    /// perpendicular to `Tgt` the function returns with `myIsDefined` still
    /// `false` and `myCurvRef` untouched (`cxx:311-318`, the OCCT `NYIxpu`
    /// comment at `:313-315`).
    #[allow(clippy::too_many_arguments)]
    pub fn reset(
        &mut self,
        tgt: &GpDir,
        norm: &GpDir,
        max_d: &GpDir,
        min_d: &GpDir,
        max_curv: f64,
        min_curv: f64,
    ) {
        self.is_defined = false;

        self.norm = *norm;
        self.tgt = *tgt;
        self.beafter = norm.crossed(tgt).unwrap_or_else(|_| GpDir::default_dir());

        let tola = ANGULAR;
        let curismax = max_d.dot(tgt).abs() < tola;
        let curismin = min_d.dot(tgt).abs() < tola;

        if max_curv.abs() < tola && min_curv.abs() < tola {
            self.reset_line(tgt, norm);
            return;
        }

        if !curismax && !curismin {
            // In the plane normal to `myTgt` the boundary face is seen as a
            // boundary curve; the curvature of that curve is not computed
            // (`cxx:311-318`).
            return;
        }

        if curismax {
            self.curv_ref = max_curv.abs();
        }
        if curismin {
            self.curv_ref = min_curv.abs();
        }
        if self.curv_ref < tola {
            self.curv_ref = 0.0;
        }

        self.table.reset_ang();
        self.touch_flag = false; // eap Mar 25 2002
        self.is_defined = true;
    }

    /// `Reset(Tgt, Norm)` (`cxx:351-370`): a straight line, `myCurvRef = 0`.
    pub fn reset_line(&mut self, tgt: &GpDir, norm: &GpDir) {
        self.is_defined = false;

        self.norm = *norm;
        self.tgt = *tgt;
        self.beafter = norm.crossed(tgt).unwrap_or_else(|_| GpDir::default_dir());
        self.table.reset_ang();

        self.curv_ref = 0.0;
        self.touch_flag = false; // eap Mar 25 2002
        self.is_defined = true;
    }

    /// `Compare(Tole, Norm, MaxD, MinD, MaxCurv, MinCurv, S, O)`
    /// (`cxx:371-490`).
    ///
    /// `S` is the topological orientation of the face (how the curve crosses
    /// it), `O` the curve's position on the face (before/after/in).
    #[allow(clippy::too_many_arguments)]
    pub fn compare_curved(
        &mut self,
        tole: f64,
        norm: &GpDir,
        max_d: &GpDir,
        min_d: &GpDir,
        max_curv: f64,
        min_curv: f64,
        s: Orientation,
        o: Orientation,
    ) {
        if !self.is_defined {
            return;
        }

        let mut curv = 0.0;
        let tola = if tole > 0.0 { tole } else { ANGULAR };
        let curismax = max_d.dot(&self.tgt).abs() < tola;
        let curismin = min_d.dot(&self.tgt).abs() < tola;
        if !curismax && !curismin {
            // Same `NYIxpu` gap as in `Reset` (`cxx:392-401`).
            self.is_defined = false;
            return;
        }
        if curismax {
            curv = max_curv.abs();
        }
        if curismin {
            curv = min_curv.abs();
        }
        if self.curv_ref < tola {
            curv = 0.0;
        }
        let diron_f = n_c_inside_s(&self.tgt, norm);
        let prod = diron_f.crossed(norm).map(|d| d.dot(&self.tgt)).unwrap_or(0.0);
        if prod < 0.0 {
            curv = -curv;
        }

        let ang = fun_ang(&self.beafter, &self.tgt, norm, o);

        let (mut i, mut j) = fun_get_sta(ang, tola);

        let kmax = if o == Orientation::Internal { 2 } else { 1 };
        for k in 1..=kmax {
            if k == 2 {
                // Get the opposite Ang.
                i = oo(i);
                j = oo(j);
            }
            let i0 = i == 0;
            let j0 = j == 0;
            let nmax = if i0 || j0 { 2 } else { 1 };
            for n in 1..=nmax {
                if i0 {
                    i = n;
                }
                if j0 {
                    j = n;
                }

                let ang_ref = self.table.ang(i, j);
                let ori_ref = self.table.ori(i, j);
                let curv_ref = self.table.curv(i, j);
                let refn = ref_nearest_full(
                    i,
                    j,
                    self.curv_ref,
                    ang_ref,
                    ori_ref,
                    curv_ref,
                    ang,
                    s,
                    curv,
                    tola,
                    &mut self.touch_flag,
                );
                if refn == M_UNKNOWN {
                    self.is_defined = false;
                    return;
                }
                if refn > 0 {
                    self.table.set_ang(i, j, ang);
                    self.table
                        .set_ori(i, j, if refn == M_OINTERNAL { Orientation::Internal } else { s });
                    self.table.set_curv(i, j, curv);
                }
            }
        }
    }

    /// `Compare(Tole, Norm, S, O)` (`cxx:491-549`) — the plane or cylindric
    /// face overload.
    pub fn compare_plane(
        &mut self,
        tole: f64,
        norm: &GpDir,
        s: Orientation,
        o: Orientation,
    ) {
        if !self.is_defined {
            return;
        }

        let ang = fun_ang(&self.beafter, &self.tgt, norm, o);
        let tola = if tole > 0.0 { tole } else { ANGULAR };

        let (mut i, mut j) = fun_get_sta(ang, tola);

        let kmax = if o == Orientation::Internal { 2 } else { 1 };
        for k in 1..=kmax {
            if k == 2 {
                i = oo(i);
                j = oo(j);
            }

            let i0 = i == 0;
            let j0 = j == 0;
            let nmax = if i0 || j0 { 2 } else { 1 };
            for n in 1..=nmax {
                if i0 {
                    i = n;
                }
                if j0 {
                    j = n;
                }

                let refn = ref_nearest_plain(self.table.ang(i, j), self.table.ori(i, j), ang, s, tola);
                if refn == M_UNKNOWN {
                    self.is_defined = false;
                    return;
                }

                if refn > 0 {
                    self.table.set_ang(i, j, ang);
                    self.table
                        .set_ori(i, j, if refn == M_OINTERNAL { Orientation::Internal } else { s });
                }
            }
        }
    }

    /// `StateBefore()` (`cxx:588-616`).
    pub fn state_before(&self) -> State {
        if !self.is_defined {
            return State::Unknown;
        }
        let mut before = get_state(self, BEFORE, BEFORE);
        if before == State::Unknown {
            // Looking back in `before` for defined states.
            before = get_state(self, AFTER, BEFORE);
            // eap Mar 25 2002
            if self.touch_flag {
                if before == State::Out {
                    before = State::In;
                } else if before == State::In {
                    before = State::Out;
                }
            }
        }
        before
    }

    /// `StateAfter()` (`cxx:618-644`).
    pub fn state_after(&self) -> State {
        if !self.is_defined {
            return State::Unknown;
        }
        let mut after = get_state(self, AFTER, AFTER);
        if after == State::Unknown {
            // Looking back in `before` for defined states.
            after = get_state(self, BEFORE, AFTER);
            // eap Mar 25 2002
            if self.touch_flag {
                if after == State::Out {
                    after = State::In;
                } else if after == State::In {
                    after = State::Out;
                }
            }
        }
        after
    }

    /// `GetBefore(Tran)` (`cxx:646-658`).
    pub fn get_before(tran: Orientation) -> State {
        match tran {
            Orientation::Forward | Orientation::External => State::Out,
            Orientation::Reversed | Orientation::Internal => State::In,
        }
    }

    /// `GetAfter(Tran)` (`cxx:660-672`).
    pub fn get_after(tran: Orientation) -> State {
        match tran {
            Orientation::Forward | Orientation::Internal => State::In,
            Orientation::Reversed | Orientation::External => State::Out,
        }
    }
}

impl Default for SurfaceTransition {
    fn default() -> Self {
        Self::new()
    }
}

/// `FUN_getstate(Ang, Ori, iSTA, iINDEX)` (`cxx:554-586`).
fn get_state(t: &SurfaceTransition, i_sta: i32, i_index: i32) -> State {
    let a1 = t.table.ang(i_sta, 1);
    let a2 = t.table.ang(i_sta, 2);
    let undef1 = a1 == ANG_UNDEF;
    let undef2 = a2 == ANG_UNDEF;
    if undef1 && undef2 {
        return State::Unknown;
    }

    if undef1 || undef2 {
        let jok = if undef1 { 2 } else { 1 };
        let o = t.table.ori(i_sta, jok);
        return if i_index == BEFORE {
            SurfaceTransition::get_before(o)
        } else {
            SurfaceTransition::get_after(o)
        };
    }

    let o1 = t.table.ori(i_sta, 1);
    let o2 = t.table.ori(i_sta, 2);
    let st1 = if i_index == BEFORE {
        SurfaceTransition::get_before(o1)
    } else {
        SurfaceTransition::get_after(o1)
    };
    let st2 = if i_index == BEFORE {
        SurfaceTransition::get_before(o2)
    } else {
        SurfaceTransition::get_after(o2)
    };
    if st1 != st2 {
        return State::Unknown; // Incoherent data
    }
    st1
}

/// `Precision::Confusion()` is included by `TopTrans_SurfaceTransition.cxx`
/// (`cxx:19`) but not read by any of its branches; the import is kept so the
/// dependency list matches the C++ translation unit.
#[allow(dead_code)]
const _UNUSED_CONFUSION: f64 = CONFUSION;
