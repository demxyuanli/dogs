//! `TopOpeBRepDS_Edge3dInterferenceTool` — 3d complex transition on an edge.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_Edge3dInterferenceTool.hxx:26-70`,
//! `TopOpeBRepDS_Edge3dInterferenceTool.cxx:1-341` (OCCT 8.0.0).
//!
//! The tool reduces the interferences `I = (T on F, G = point/vertex,
//! S = <E>)` laid on a reference edge `Eref` into a transition on `Eref`,
//! using `TopTrans_SurfaceTransition` (`myTool`) fed with the tangents of
//! `Eref` and of `E` at the common point (`cxx:270-271`, `:327`).
//!
//! Translation notes:
//! - `myIsVertex` is never initialized by the constructor
//!   (`cxx:158-162` initializes only `myFaceOriented` and `myrefdef`); the
//!   port defaults it to `POINT` and the only writer is [`InitPointVertex`],
//!   which every caller in `TopOpeBRepDS_EIR` runs before `Init`.
//! - `FUN_Raise()` (`cxx:35-38`) throws `Standard_ProgramError`; every call
//!   here is an `Err` on the same branch.
//! - `BRepAdaptor_Curve2d(EE, FF).Value(pOO)` (`cxx:302-305`) is
//!   `BRep_Tool::CurveOnSurface(E, F, f, l)` (`BRepAdaptor_Curve2d.cxx:55-64`
//!   never builds); [`crate::boptools_2d::curve_on_surface_range`] is that
//!   read. OCCT would evaluate a null adaptor curve (undefined); the port
//!   answers "nothing to compare" instead.
//! - `FUN_tool_nggeomF(parE, E, F, dd)` (`cxx:140`) and the `(uv, F)` overload
//!   (`cxx:327`) are [`fun_tool_ng_geom_f_par`] / [`fun_tool_ng_geom_f`].
//! - `FUN_tool_parVonE` / `FUN_tool_projPonE` (`cxx:97`, `:104`) are
//!   [`fun_tool_par_von_e`] / [`fun_tool_proj_pon_e`].

use occt_core::gp::{GpDir, GpPnt, GpPnt2d};

use crate::abs::{Orientation, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};

use super::interference::{Interference, InterferenceData};
use super::kind::Kind;
use crate::top_trans_surface_transition::SurfaceTransition;
use crate::topopebrep::tool::{
    fun_tool_ng_geom_f, fun_tool_ng_geom_f_par, fun_tool_orient_ein_f_forward,
    fun_tool_par_von_e, fun_tool_value, fun_tool_proj_pon_e,
};

/// `#define POINT (0)` (`cxx:39`).
pub const POINT: i32 = 0;
/// `#define VERTEXonref (1)` (`cxx:40`).
pub const VERTEX_ON_REF: i32 = 1;
/// `#define VERTEXonOO (2)` (`cxx:41`).
pub const VERTEX_ON_OO: i32 = 2;
/// `#define VERTEXonOref (3)` (`cxx:42`).
pub const VERTEX_ON_O_REF: i32 = 3;

/// `M_FORWARD(st)` (`cxx:31`).
fn m_forward(s: Orientation) -> bool {
    s == Orientation::Forward
}

/// `M_REVERSED(st)` (`cxx:32`).
fn m_reversed(s: Orientation) -> bool {
    s == Orientation::Reversed
}

/// `static bool FUN_hasparam(I, paronE)` (`cxx:58-86`).
///
/// The C++ returns `true` untouched when the geometry type is neither
/// `POINT` nor `VERTEX` (`cxx:84`); that arm is kept.
fn fun_hasparam(i: &Interference, paron_e: &mut f64) -> bool {
    match i.geometry_type() {
        Kind::Point => match i.data() {
            InterferenceData::CurvePoint { .. } => match i.parameter() {
                Some(p) => {
                    *paron_e = p;
                    true
                }
                None => false,
            },
            _ => false,
        },
        Kind::Vertex => match i.data() {
            InterferenceData::EdgeVertex(_) => match i.parameter() {
                Some(p) => {
                    *paron_e = p;
                    true
                }
                None => false,
            },
            _ => false,
        },
        _ => true,
    }
}

/// `static bool FUN_paronOOE(OOE, IsVertex, VonOO, P3d, paronOOE)`
/// (`cxx:88-116`): the parameter of the point/vertex `G` on `E`.
fn fun_par_on_ooe(
    oo_e: &TopoShape,
    is_vertex: i32,
    von_oo: Option<&TopoShape>,
    p3d: &GpPnt,
    paron_oo_e: &mut f64,
) -> bool {
    let has_von_oo = is_vertex > 1;
    if has_von_oo {
        return match von_oo {
            Some(v) => match fun_tool_par_von_e(v, oo_e) {
                Some(p) => {
                    *paron_oo_e = p;
                    true
                }
                None => false,
            },
            None => false,
        };
    }
    let (par, dist) = match fun_tool_proj_pon_e(p3d, oo_e) {
        Some(x) => x,
        None => return false,
    };
    *paron_oo_e = par;
    let tol1 = BRepTool::edge_tolerance(&Edge(oo_e.clone()));
    let tol2 = tol1 * 1.0e3;
    let mut tol = tol2;
    if tol > 1.0e-2 {
        tol = 1.0e-2;
    }
    dist <= tol
}

/// `static bool FUN_keepIonF(tgref, parE, E, F, tola)` (`cxx:118-149`):
/// `true` when the interference must be added to the 3d transition on `E`.
fn fun_keep_ion_f(tgref: &GpDir, par_e: f64, e: &TopoShape, f: &TopoShape, tola: f64) -> bool {
    let tmp = match crate::topopebrep::tool::tg_geom_e(par_e, e) {
        Some(v) => v,
        None => return false,
    };
    let tg_e = match GpDir::new(tmp.x(), tmp.y(), tmp.z()) {
        Ok(d) => d,
        Err(_) => return false,
    };
    let prod = tgref.dot(&tg_e).abs();
    if (1.0 - prod).abs() < tola {
        // <Eref> & <E> are tangent edges
        return false;
    }
    let dd = match fun_tool_ng_geom_f_par(par_e, e, f) {
        Some(d) => d,
        None => return false,
    };
    let prod = match tgref.crossed(&tg_e) {
        Ok(c) => c.dot(&dd).abs(),
        Err(_) => return false,
    };
    (1.0 - prod).abs() >= tola
}

/// `TopOpeBRepDS_Edge3dInterferenceTool` (`hxx:26-70`).
#[derive(Debug, Clone)]
pub struct Edge3dInterferenceTool {
    /// `myIsVertex`
    is_vertex: i32,
    /// `myVonOO` — only for `VERTEXonOO || VERTEXonOref`.
    von_oo: Option<TopoShape>,
    /// `myP3d` — only for `POINT || VERTEXonref`.
    p3d: GpPnt,
    /// `myTole` = `Precision::Angular()` (`cxx:199`).
    tole: f64,
    /// `myFaceOriented`
    face_oriented: i32,
    /// `myTgtref`
    tgtref: GpDir,
    /// `myTool`
    tool: SurfaceTransition,
    /// `myrefdef`
    refdef: bool,
}

impl Default for Edge3dInterferenceTool {
    fn default() -> Self {
        Self::new()
    }
}

impl Edge3dInterferenceTool {
    /// `TopOpeBRepDS_Edge3dInterferenceTool()` (`cxx:156-163`):
    /// `myFaceOriented = 0`, `myrefdef = false`; see the module note on
    /// `myIsVertex`.
    pub fn new() -> Self {
        Self {
            is_vertex: POINT,
            von_oo: None,
            p3d: GpPnt::new(0.0, 0.0, 0.0),
            tole: 0.0,
            face_oriented: 0,
            tgtref: GpDir::new(1.0, 0.0, 0.0).unwrap_or_else(|_| GpDir::default_dir()),
            tool: SurfaceTransition::new(),
            refdef: false,
        }
    }

    /// `InitPointVertex(IsVertex, VonOO)` (`cxx:171-179`).
    pub fn init_point_vertex(&mut self, is_vertex: i32, von_oo: Option<TopoShape>) {
        self.is_vertex = is_vertex;
        if is_vertex > 1 {
            self.von_oo = von_oo;
        }
    }

    /// `Init(Eref, E, F, I)` (`cxx:189-271`).
    ///
    /// `myrefdef` stays `false` on every early return; `Add` then re-runs
    /// `Init` on the next interference (`cxx:282-291`).
    pub fn init(
        &mut self,
        e_ref: &TopoShape,
        e: &TopoShape,
        f: &TopoShape,
        i: &Interference,
    ) -> Result<(), String> {
        self.refdef = false;

        self.tole = occt_core::precision::ANGULAR; // NYI

        let mut pref = 0.0;
        if !fun_hasparam(i, &mut pref) {
            return Err("Edge3dInterferenceTool".to_string());
        }
        // <myP3d>
        self.p3d = match fun_tool_value(pref, e_ref) {
            Some(p) => p,
            None => return Err("Edge3dInterferenceTool".to_string()),
        };
        let tmp = match crate::topopebrep::tool::tg_geom_e(pref, e_ref) {
            Some(v) => v,
            None => return Err("Edge3dInterferenceTool".to_string()),
        };
        let tgref = match GpDir::new(tmp.x(), tmp.y(), tmp.z()) {
            Ok(d) => d,
            Err(_) => return Err("Edge3dInterferenceTool".to_string()),
        };

        let mut p_oo = 0.0;
        if !fun_par_on_ooe(e, self.is_vertex, self.von_oo.as_ref(), &self.p3d, &mut p_oo) {
            return Err("Edge3dInterferenceTool".to_string());
        }
        let tmp = match crate::topopebrep::tool::tg_geom_e(p_oo, e) {
            Some(v) => v,
            None => return Err("Edge3dInterferenceTool".to_string()),
        };
        let tg_oo = match GpDir::new(tmp.x(), tmp.y(), tmp.z()) {
            Ok(d) => d,
            Err(_) => return Err("Edge3dInterferenceTool".to_string()),
        };

        let mut dot = tgref.dot(&tg_oo);
        dot = 1.0 - dot.abs();
        let tola = occt_core::precision::CONFUSION;
        let esdm = dot.abs() < tola;
        if esdm {
            return Ok(());
        }
        // NYI : il faut rejeter les interf I = (T,G,S=E) / E sdm with Eref

        // <myrefdef>
        if !fun_keep_ion_f(&tgref, p_oo, e, f, self.tole) {
            // <Eref> is tangent to <F>,
            // If the transition is FORWARD or REVERSED, it describes a 2d
            // transition (while crossing <E> on <F>), we do not keep it.
            let o = i
                .transition()
                .orientation(State::In, crate::abs::ShapeType::Face);
            let is2d = m_forward(o) || m_reversed(o);
            if is2d {
                return Ok(());
            }
        }
        self.refdef = true;

        // <myFaceOriented>
        self.face_oriented = i.transition().index()?;

        // <myTgtref>
        self.tgtref = tgref;

        let norm = tg_oo.crossed(&tgref).unwrap_or_else(|_| {
            GpDir::new(0.0, 0.0, 1.0).unwrap_or_else(|_| GpDir::default_dir())
        });
        self.tool.reset_line(&tg_oo, &norm);
        Ok(())
    }

    /// `Add(Eref, E, F, I)` (`cxx:277-328`).
    pub fn add(
        &mut self,
        e_ref: &TopoShape,
        e: &TopoShape,
        f: &TopoShape,
        i: &Interference,
    ) -> Result<(), String> {
        if !self.refdef {
            self.init(e_ref, e, f, i)?;
        }
        if !self.refdef {
            return Ok(());
        }

        let edge = Edge(e.clone());
        let face = Face(f.clone());

        let mut p_oo = 0.0;
        if !fun_par_on_ooe(e, self.is_vertex, self.von_oo.as_ref(), &self.p3d, &mut p_oo) {
            return Ok(());
        }
        let uv: GpPnt2d = match crate::boptools_2d::curve_on_surface_range(&edge, &face) {
            Some((pc, _first, _last)) => {
                let p = pc.value(p_oo);
                GpPnt2d::new(p.x(), p.y())
            }
            None => return Ok(()),
        };

        if !fun_keep_ion_f(&self.tgtref, p_oo, e, f, self.tole) {
            let o = i
                .transition()
                .orientation(State::In, crate::abs::ShapeType::Face);
            let is2d = m_forward(o) || m_reversed(o);
            if is2d {
                return Ok(());
            }
        }

        let oriloc = i
            .transition()
            .orientation(State::In, crate::abs::ShapeType::Face);
        // xpu : 30/12/97
        let oritan = match fun_tool_orient_ein_f_forward(e, f) {
            Some(o) => o,
            None => return Ok(()),
        };

        let norm = match fun_tool_ng_geom_f(&uv, &face) {
            Some(n) => n,
            None => return Ok(()),
        };
        self.tool.compare_plane(self.tole, &norm, oriloc, oritan);
        Ok(())
    }

    /// `Transition(I)` (`cxx:332-341`).
    pub fn transition(&self, i: &mut Interference) {
        let (stb, sta) = (self.tool.state_before(), self.tool.state_after());
        i.set_support(self.face_oriented);
        let t = i.change_transition();
        t.set(
            stb,
            sta,
            crate::abs::ShapeType::Face,
            crate::abs::ShapeType::Face,
        );
    }
}
