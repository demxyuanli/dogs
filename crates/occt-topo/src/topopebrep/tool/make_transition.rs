//! `TopOpeBRepTool_makeTransition` — the transition of an edge against a face
//! at a point where the two are tangent.
//!
//! Source: `TKBool/TopOpeBRepTool/TopOpeBRepTool_makeTransition.cxx:64-877`
//! and `.hxx:27-74`.
//!
//! `MkT2donE` (`cxx:351-459`) handles the case where the edge is *in* the face
//! in 2-D: a non-tangent crossing only needs the sign of `tgE . xxES`; a
//! tangent crossing between two quadrics compares their curvatures
//! (`FUN_mkT2dquad`, `cxx:155-190`, `:253-347`). `MkT3onE` (`cxx:728-807`)
//! handles a 3-D tangent edge through the apparent contour of the face, and
//! `MkT3dproj` (`cxx:700-710`) classifies the two nearby points by projection.
//!
//! `MkTonE` (`cxx:811-877`) is the entry point: it tries `MkT3onE`, falls back
//! to `MkT3dproj`, and then repeatedly re-projects while either state is still
//! `ON`/`UNKNOWN`, growing `myfactor` by `(1 - factor) / 5` each round — five
//! rounds at most.
//!
//! `FUN_tolang` (`cxx:57-60`) is `Precision::Angular() * 1.e6`.

use occt_core::gp::{GpDir, GpDir2d, GpPnt2d, GpVec, GpVec2d};

use crate::abs::State;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};

use super::adaptor::{CurveAdaptor, SurfaceAdaptor};
use super::closed::{ori_in_sor, CLOSING, INTERNAL};
use super::curvature::{
    curv_e, curv_f, get_stp_3d_f, getduv, is_quad_e, is_quad_f, nt, stuv_f, xx, INFFIRST, SUPLAST,
};
use super::project_curve::{fun_tool_bounds, fun_tool_value};
use super::tangent::tg_geom_e;

/// `BEFORE (1)` / `AFTER (2)` (`cxx:28-29`).
const BEFORE: i32 = 1;
const AFTER: i32 = 2;

/// `FUN_mkT2dquad`'s answer codes (`cxx:31-37`).
const IS_IN_IF_H1: i32 = 1;
const IS_IN_IF_H2: i32 = 2;
const IS_ON2: i32 = 21;
const IS_ON2_IF_SS: i32 = 10;
const IS_IN2_IF_SS: i32 = 11;
const IS_OU2_IF_SS: i32 = 12;

/// `FUN_tolang()` (`cxx:57-60`): `Precision::Angular() * 1.e6`.
fn fun_tolang() -> f64 {
    occt_core::precision::ANGULAR * 1.0e6
}

/// `FUN_nullcurv(curv)` (`cxx:147-151`): `curv < Precision::Confusion() * 1.e2`.
fn fun_nullcurv(curv: f64) -> bool {
    curv < occt_core::precision::CONFUSION * 1.0e2
}

/// `FUN_mkT2dquad(curvC1, curvC2)` (`cxx:155-190`).
fn mkt2dquad_code(curv_c1: f64, curv_c2: f64) -> i32 {
    let nullc1 = fun_nullcurv(curv_c1);
    let nullc2 = fun_nullcurv(curv_c2);
    if nullc2 && nullc1 {
        return IS_ON2;
    }
    if nullc2 {
        return IS_IN_IF_H1;
    }
    if nullc1 {
        return IS_IN_IF_H2;
    }
    let samec = (curv_c2 - curv_c1).abs() < 1.0e-2; // NYITOLXPU kpartkoletge
    if samec {
        return IS_ON2_IF_SS;
    }
    if curv_c1 > curv_c2 {
        IS_IN2_IF_SS
    } else {
        IS_OU2_IF_SS
    }
}

/// `FUN_getnearpar(e, par, f, l, factor, sta, nearpar)` (`cxx:192-215`).
///
/// `nearpar = par - factor * (l - f)` before the point, `par + factor * (l - f)`
/// after it; `None` is the C++ `false` when the point already sits on the
/// range end the step would move past.
fn get_near_par(
    e: &TopoShape,
    par: f64,
    f: f64,
    l: f64,
    factor: f64,
    sta: i32,
) -> Option<f64> {
    let bc = CurveAdaptor::of_edge(&Edge(e.clone()));
    let tol1d = bc.resolution(bc.tolerance());
    let onf = (par - f).abs() < tol1d;
    let onl = (par - l).abs() < tol1d;
    if onf && sta == BEFORE {
        return None;
    }
    if onl && sta == AFTER {
        return None;
    }
    Some(if sta == BEFORE {
        par - factor * (l - f)
    } else {
        par + factor * (l - f)
    })
}

/// `FUN_tg(e, par, f, l, factor, tg, st)` (`cxx:217-250`): the edge tangent at
/// the point before (`st = BEFORE`) or after the interference point, whichever
/// can be evaluated first.
fn fun_tg(e: &TopoShape, par: f64, f: f64, l: f64, factor: f64) -> Option<(GpDir, i32)> {
    let mut st = BEFORE;
    for nite in 1..=2 {
        if nite == 2 {
            st = AFTER;
        }
        let Some(pn) = get_near_par(e, par, f, l, factor, st) else {
            continue;
        };
        let Some(tmp) = tg_geom_e(pn, e) else {
            continue;
        };
        let tg = GpDir::from_vec(&tmp).ok()?;
        return Some((tg, st));
    }
    None
}

/// `FUN_getsta(mkt, tga1, tga2, xx2, sta)` (`cxx:253-300`).
fn get_sta(mkt: i32, tga1: &GpDir, tga2: &GpDir, xx2: &GpDir) -> Option<State> {
    if mkt == IS_IN_IF_H1 {
        let dot = tga1.dot(xx2);
        return Some(if dot > 0.0 { State::In } else { State::Out });
    }
    if mkt == IS_IN_IF_H2 {
        let dot = tga2.dot(xx2);
        return Some(if dot < 0.0 { State::In } else { State::Out });
    }
    if mkt == IS_ON2_IF_SS {
        let ssided = tga1.dot(tga2) > 0.0;
        return Some(if ssided { State::On } else { State::In });
    }
    if mkt == IS_IN2_IF_SS {
        let dot = tga1.dot(xx2);
        return Some(if dot < 0.0 { State::Out } else { State::In });
    }
    if mkt == IS_OU2_IF_SS {
        let dot = tga2.dot(xx2);
        return Some(if dot < 0.0 { State::In } else { State::Out });
    }
    Some(State::On) // mkt == IS_ON2
}

/// `FUN_mkT2dquad(e1, par1, f1, l1, e2, par2, mkt, xx2, factor, sta)`
/// (`cxx:303-347`): the two near tangents, oriented after the point, fed to
/// [`get_sta`].
fn mkt2dquad(
    e1: &TopoShape,
    par1: f64,
    f1: f64,
    l1: f64,
    e2: &TopoShape,
    par2: f64,
    mkt: i32,
    xx2: &GpDir,
    factor: f64,
) -> Option<State> {
    let mut tga1: Option<GpDir> = None;
    let mut tga2: Option<GpDir> = None;
    let mk1 = mkt == IS_IN_IF_H1 || mkt == IS_ON2_IF_SS || mkt == IS_IN2_IF_SS;
    if mk1 {
        let (tgnear1, st1) = fun_tg(e1, par1, f1, l1, factor)?;
        tga1 = Some(if st1 == AFTER { tgnear1 } else { tgnear1.reversed() });
    }
    let mk2 = mkt == IS_IN_IF_H2 || mkt == IS_ON2_IF_SS || mkt == IS_OU2_IF_SS;
    if mk2 {
        let (f2, l2) = fun_tool_bounds(e2);
        let (tgnear2, st2) = fun_tg(e2, par2, f2, l2, factor)?;
        tga2 = Some(if st2 == AFTER { tgnear2 } else { tgnear2.reversed() });
    }
    let zero = GpDir::default_dir();
    get_sta(mkt, tga1.as_ref().unwrap_or(&zero), tga2.as_ref().unwrap_or(&zero), xx2)
}

/// `FUN_getnearuv(f, uv, factor, sta, duv, nearuv)` (`cxx:461-525`): the uv
/// step from `uv`, unwrapped across a closed direction.
fn get_near_uv(f: &Face, uv: &GpPnt2d, factor: f64, sta: i32, duv: &GpDir2d) -> Option<GpPnt2d> {
    let bs = SurfaceAdaptor::of_face(f);
    let mut xuv = GpVec2d::new(duv.x() * factor, duv.y() * factor);
    if sta == BEFORE {
        xuv.reverse();
    }
    let nearuv = uv.translated_vec(&xuv);

    let (onu, onv) = stuv_f(uv, f);
    let uok = onu == 0;
    let vok = onv == 0;
    if uok && vok {
        return Some(nearuv);
    }

    let mut nearu = nearuv.x();
    let mut nearv = nearuv.y();
    let s = bs.surface()?;

    if onu == INFFIRST || onu == SUPLAST {
        if !s.is_u_closed() {
            return None;
        }
        let uper = s.u_period();
        if onu == INFFIRST {
            nearu += uper;
        } else {
            nearu -= uper;
        }
    }
    if onv == INFFIRST || onv == SUPLAST {
        if !s.is_v_closed() {
            return None;
        }
        let vper = s.v_period();
        if onv == INFFIRST {
            nearv += vper;
        } else {
            nearv -= vper;
        }
    }
    Some(GpPnt2d::new(nearu, nearv))
}

/// `FUN_tgef(f, uv, duv, factor, tge, tg0, tgef, st)` (`cxx:527-564`).
///
/// Because `ntf ^ tge = tg0`, the near tangent of the projected contour is
/// assumed to satisfy `nt(nearuv) ^ tgef(nearuv) = tg0`, i.e.
/// `tgef = tg0 ^ nt`.
fn fun_tgef(
    f: &Face,
    uv: &GpPnt2d,
    duv: &GpDir2d,
    factor: f64,
    tg0: &GpDir,
) -> Option<(GpDir, i32)> {
    let mut st = BEFORE;
    for nite in 1..=2 {
        if nite == 2 {
            st = AFTER;
        }
        let Some(nearuv) = get_near_uv(f, uv, factor, st, duv) else {
            continue;
        };
        let Some(ntf) = nt(&nearuv, f) else {
            return None;
        };
        let tgef = tg0.crossed(&ntf).ok()?;
        return Some((tgef, st));
    }
    None
}

/// `FUN_mkT3dquad(e, pf, pl, par, f, uv, tge, ntf, mkt, factor, sta)`
/// (`cxx:566-619`).
fn mkt3dquad(
    e: &TopoShape,
    pf: f64,
    pl: f64,
    par: f64,
    f: &Face,
    uv: &GpPnt2d,
    tge: &GpDir,
    ntf: &GpDir,
    mkt: i32,
    factor: f64,
) -> Option<State> {
    let xxef = ntf.reversed();
    let tg0 = ntf.crossed(tge).ok()?;

    let mut tgae: Option<GpDir> = None;
    let mut tgaef: Option<GpDir> = None;
    let mke = mkt == IS_IN_IF_H1 || mkt == IS_ON2_IF_SS || mkt == IS_IN2_IF_SS;
    if mke {
        let (tgnear, st) = fun_tg(e, par, pf, pl, factor)?;
        tgae = Some(if st == AFTER { tgnear } else { tgnear.reversed() });
    }
    let mkef = mkt == IS_IN_IF_H2 || mkt == IS_ON2_IF_SS || mkt == IS_OU2_IF_SS;
    if mkef {
        // choice : tgdir(ef,uv) = tgdir(e,pare)
        let fac3d = 0.01; // 0.12345; not only planar faces
        let duv = getduv(f, uv, &GpVec::new(tge.x(), tge.y(), tge.z()), fac3d)?;
        let (tgnear, st) = fun_tgef(f, uv, &duv, factor, &tg0)?;
        tgaef = Some(if st == AFTER { tgnear } else { tgnear.reversed() });
    }
    let zero = GpDir::default_dir();
    get_sta(
        mkt,
        tgae.as_ref().unwrap_or(&zero),
        tgaef.as_ref().unwrap_or(&zero),
        &xxef,
    )
}

/// `FUN_stawithES(tgE, xxES, st)` (`cxx:621-646`).
fn stawith_es(tg_e: &GpDir, xx_es: &GpDir, st: i32) -> State {
    let prod = tg_e.dot(xx_es);
    let tola = fun_tolang();
    if prod.abs() < tola {
        return State::Unknown;
    }
    let positive = prod > 0.0;
    if positive {
        if st == BEFORE {
            State::Out
        } else {
            State::In
        }
    } else if st == BEFORE {
        State::In
    } else {
        State::Out
    }
}

/// `FUN_stawithES(tgE, xxES, st, stt)` (`cxx:648-666`): keep `stt` when the
/// point stays `IN`, downgrade to `OUT` when it leaves.
fn stawith_es_stt(tg_e: &GpDir, xx_es: &GpDir, st: i32, stt: State) -> State {
    if stt == State::Unknown {
        return State::Unknown;
    }
    let st_es = stawith_es(tg_e, xx_es, st);
    if stt == State::In || stt == State::On {
        if st_es == State::In {
            stt
        } else {
            State::Out
        }
    } else {
        State::Unknown
    }
}

/// `FUN_staproj(e, pf, pl, pe, factor, st, f, sta)` (`cxx:668-695`): the state
/// of the point at `pe` stepped `factor` towards `pf`/`pl`.
fn staproj(
    e: &TopoShape,
    pf: f64,
    pl: f64,
    pe: f64,
    factor: f64,
    st: i32,
    f: &Face,
) -> Option<State> {
    let par = get_near_par(e, pe, pf, pl, factor, st)?;
    let pt = fun_tool_value(par, e)?;
    let (_uv, sta) = get_stp_3d_f(&pt, f)?;
    Some(sta)
}

/// `TopOpeBRepTool_makeTransition` (`hxx:27-74`).
#[derive(Debug, Clone)]
pub struct MakeTransition {
    e: TopoShape,
    pb: f64,
    pa: f64,
    pe: f64,
    fs: TopoShape,
    uv: GpPnt2d,
    has_es: bool,
    es: TopoShape,
    pes: f64,
    is_t2d: bool,
    factor: f64,
}

impl Default for MakeTransition {
    fn default() -> Self {
        Self::new()
    }
}

impl MakeTransition {
    /// `TopOpeBRepTool_makeTransition::TopOpeBRepTool_makeTransition()`
    /// (`cxx:64`): the default constructor leaves the geometry defaulted; the
    /// C++ members are uninitialised, the port starts from neutral values.
    pub fn new() -> Self {
        Self {
            e: TopoShape::new(crate::abs::ShapeType::Edge),
            pb: 0.0,
            pa: 0.0,
            pe: 0.0,
            fs: TopoShape::new(crate::abs::ShapeType::Face),
            uv: GpPnt2d::new(0.0, 0.0),
            has_es: false,
            es: TopoShape::new(crate::abs::ShapeType::Edge),
            pes: 0.0,
            is_t2d: false,
            factor: 0.0,
        }
    }

    /// `Initialize(E, pbef, paft, parE, FS, uv, factor)` (`cxx:68-98`).
    ///
    /// `false` when `E` is degenerated (`cxx:75-79`), when `factor` is outside
    /// `[0, 1]` (`cxx:90-94`) or when `EdgeONFace` fails.
    pub fn initialize(
        &mut self,
        e: &TopoShape,
        pbef: f64,
        paft: f64,
        par_e: f64,
        fs: &TopoShape,
        uv: &GpPnt2d,
        factor: f64,
    ) -> bool {
        if BRepTool::is_degenerated(&Edge(e.clone())) {
            return false;
        }
        self.e = e.clone();
        self.pb = pbef;
        self.pa = paft;
        self.pe = par_e;
        self.fs = fs.clone();
        self.uv = *uv;
        self.factor = factor;
        self.has_es = false;

        if !(0.0..=1.0).contains(&factor) {
            return false;
        }
        match super::curvature::edge_on_face(par_e, e, uv, &Face(fs.clone())) {
            Some(ison_fa) => {
                self.is_t2d = ison_fa;
                true
            }
            None => false,
        }
    }

    /// `Setfactor(factor)` (`cxx:103-106`).
    pub fn set_factor(&mut self, factor: f64) {
        self.factor = factor;
    }

    /// `Getfactor()` (`cxx:110-113`).
    pub fn factor(&self) -> f64 {
        self.factor
    }

    /// `IsT2d()` (`cxx:117-120`).
    pub fn is_t2d(&self) -> bool {
        self.is_t2d
    }

    /// `SetRest(ES, parES)` (`cxx:124-137`): `false` for a degenerated
    /// restriction.
    pub fn set_rest(&mut self, es: &TopoShape, par_es: f64) -> bool {
        if BRepTool::is_degenerated(&Edge(es.clone())) {
            return false;
        }
        self.has_es = true;
        self.es = es.clone();
        self.pes = par_es;
        // nyi check <ES> is edge of <myFS>
        true
    }

    /// `HasRest()` (`cxx:142-145`).
    pub fn has_rest(&self) -> bool {
        self.has_es
    }

    /// `MkT2donE(Stb, Sta)` (`cxx:351-459`).
    pub fn mk_t2don_e(&self, stb: &mut State, sta: &mut State) -> bool {
        if !self.is_t2d {
            return false;
        }
        // E is IN 2d(FS), meets no restriction at given point :
        if !self.has_es {
            *stb = State::In;
            *sta = State::In;
            return true;
        }
        // E is IN 2d(FS), meets restriction ES at given point :
        let ori_es_fs = ori_in_sor(&self.es, &self.fs, true);
        if ori_es_fs == 0 {
            return false;
        }
        // ES is closing edge for FS, or ES is INTERNAL in FS :
        if ori_es_fs == INTERNAL || ori_es_fs == CLOSING {
            *stb = State::In;
            *sta = State::In;
            return true;
        }

        let Some(tmp) = tg_geom_e(self.pe, &self.e) else {
            return false;
        };
        let tg_e = GpDir::from_vec(&tmp).ok().unwrap_or_default();
        let face = Face(self.fs.clone());
        let Some(xx_es) = xx(&self.uv, &self.fs, self.pes, &self.es) else {
            return false;
        };

        let tola = fun_tolang();
        let dot = tg_e.dot(&xx_es);

        // E and ES are not tangent at interference point :
        let tgts = dot.abs() < tola;
        if !tgts {
            let dotpos = dot > 0.0;
            if dotpos {
                *stb = State::Out;
                *sta = State::In;
            } else {
                *stb = State::In;
                *sta = State::Out;
            }
            return true;
        }

        // E and ES are tangent, curves are quadratic : transition is
        // INTERNAL/EXTERNAL, elsewhere : transition is FOR/REV/INT/EXT
        // we then use curvatures to compute transition T :
        // xpu090299 PRO13455(E=e7, ES=e9, FS=f11)
        let Some(nt_fs) = nt(&self.uv, &face) else {
            return false;
        };
        let Some(curv_e_v) = curv_e(&self.e, self.pe, &nt_fs) else {
            return false;
        };
        let Some(curv_es_v) = curv_e(&self.es, self.pes, &nt_fs) else {
            return false;
        };

        let quad_e = is_quad_e(&self.e);
        let quad_es = is_quad_e(&self.es);
        if quad_e && quad_es {
            // should return INT/EXT
            let mkt = mkt2dquad_code(curv_e_v, curv_es_v);
            if let Some(s) = mkt2dquad(
                &self.e,
                self.pb,
                self.pa,
                self.pe,
                &self.es,
                self.pes,
                mkt,
                &xx_es,
                self.factor,
            ) {
                *stb = s;
                *sta = s;
                return true;
            }
        }

        // !!!NYI: general case : NYIKPART quadquad, only one state
        false
    }

    /// `MkT3dproj(Stb, Sta)` (`cxx:700-710`): `Stb`/`Sta` are `UNKNOWN` when
    /// the corresponding projection fails, as in the C++ (which returns `okb`
    /// and leaves `Sta` at `UNKNOWN`).
    pub fn mk_t3dproj(&self, stb: &mut State, sta: &mut State) -> bool {
        *stb = State::Unknown;
        *sta = State::Unknown;
        let face = Face(self.fs.clone());
        match staproj(&self.e, self.pb, self.pa, self.pe, self.factor, BEFORE, &face) {
            Some(s) => *stb = s,
            None => return false,
        }
        match staproj(&self.e, self.pb, self.pa, self.pe, self.factor, AFTER, &face) {
            Some(s) => *sta = s,
            None => return false,
        }
        true
    }

    /// `MkT3onE(Stb, Sta)` (`cxx:728-807`): the edge is tangent to the face, so
    /// the state comes from the apparent contour of the face; the general
    /// (non-quadric) case is `NYI` in the C++ and returns `false`.
    pub fn mk_t3on_e(&self, stb: &mut State, sta: &mut State) -> bool {
        if self.is_t2d {
            return false;
        }
        let face = Face(self.fs.clone());
        let Some(tmp) = tg_geom_e(self.pe, &self.e) else {
            return false;
        };
        let tg_e = GpDir::from_vec(&tmp).ok().unwrap_or_default();
        let Some(nt_fs) = nt(&self.uv, &face) else {
            return false;
        };

        let tola = fun_tolang();
        let dot = tg_e.dot(&nt_fs);
        if dot.abs() > tola {
            if dot > 0.0 {
                *stb = State::In;
                *sta = State::Out;
            } else {
                *stb = State::Out;
                *sta = State::In;
            }
            return true;
        }

        // E is tangent to FS at interference point
        let Some(tg0) = nt_fs.crossed(&tg_e).ok() else {
            return false;
        };
        let Some(cur_e) = curv_e(&self.e, self.pe, &tg0) else {
            return false;
        };
        let Some((cur_fs, _direct)) = curv_f(&self.fs, &self.uv, &tg0) else {
            return false;
        };

        let quad_e = is_quad_e(&self.e);
        let quad_fs = is_quad_f(&face);
        if quad_e && quad_fs {
            // should return INT/EXT
            let mkt = mkt2dquad_code(cur_e, cur_fs);
            if let Some(sta_v) = mkt3dquad(
                &self.e,
                self.pb,
                self.pa,
                self.pe,
                &face,
                &self.uv,
                &tg_e,
                &nt_fs,
                mkt,
                self.factor,
            ) {
                if self.has_es {
                    let Some(xx_es) = xx(&self.uv, &self.fs, self.pes, &self.es) else {
                        return false;
                    };
                    *stb = stawith_es_stt(&tg_e, &xx_es, BEFORE, sta_v);
                    *sta = stawith_es_stt(&tg_e, &xx_es, AFTER, sta_v);
                } else {
                    *stb = sta_v;
                    *sta = sta_v;
                }
            }
            return true;
        }

        // NYIGENERALCASE; NYIKPART quadquad, only one state
        false
    }

    /// `MkTonE(Stb, Sta)` (`cxx:811-877`).
    ///
    /// The C++ ignores the failure of the second projection (`cxx:822`,
    /// commented out) and keeps re-projecting while either state is `ON` or
    /// `UNKNOWN`, up to five rounds.
    pub fn mk_ton_e(&mut self, stb: &mut State, sta: &mut State) -> bool {
        *stb = State::Unknown;
        *sta = State::Unknown;
        if self.is_t2d {
            return self.mk_t2don_e(stb, sta);
        }

        if !self.mk_t3on_e(stb, sta) {
            self.mk_t3dproj(stb, sta);
        }

        let Some(tmp) = tg_geom_e(self.pe, &self.e) else {
            return false;
        };
        let tg_e = GpDir::from_vec(&tmp).ok().unwrap_or_default();
        let mut xx_es: Option<GpDir> = None;
        if self.has_es {
            let Some(v) = xx(&self.uv, &self.fs, self.pes, &self.es) else {
                return false;
            };
            xx_es = Some(v);
        }

        let delta = (1.0 - self.factor) / 5.0;
        let face = Face(self.fs.clone());
        for _nite in 1..=5 {
            let kob = *stb == State::On || *stb == State::Unknown;
            let koa = *sta == State::On || *sta == State::Unknown;
            if !koa && !kob {
                return true;
            }

            if kob {
                let okb = staproj(&self.e, self.pb, self.pa, self.pe, self.factor, BEFORE, &face);
                if let Some(stb_v) = okb {
                    *stb = stb_v;
                    if let Some(xx_es_v) = &xx_es {
                        *stb = stawith_es_stt(&tg_e, xx_es_v, BEFORE, stb_v);
                    }
                }
            }
            if koa {
                let oka = staproj(&self.e, self.pb, self.pa, self.pe, self.factor, AFTER, &face);
                if let Some(sta_v) = oka {
                    *sta = sta_v;
                    if let Some(xx_es_v) = &xx_es {
                        *sta = stawith_es_stt(&tg_e, xx_es_v, AFTER, sta_v);
                    }
                }
            }

            self.factor += delta;
        }
        false
    }
}
