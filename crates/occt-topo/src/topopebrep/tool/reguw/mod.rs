//! `TopOpeBRepTool_REGUW` — regularize the wires of a face.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_REGUW.cxx`
//! (OCCT 8.0.0, 1325 lines) and `TopOpeBRepTool_REGUW.hxx` (123 lines).
//!
//! Layout of the port:
//! * `mod.rs` — the state (`hxx:103-123`), the constructor, `Init` / `InitStep`
//!   (`cxx:290-318`), and the accessors (`cxx:134-186`, `:1198-1325`).
//! * [`map_s`] — `MapS` (`cxx:322-480`) and `SplitEds` (`cxx:187-276`).
//! * [`block`] — `InitBlock` (`cxx:526-674`), `NearestE` (`cxx:678-786`),
//!   `NextinBlock` (`cxx:788-921`) and the `FUN_nextdata` helper (`cxx:482-522`).
//! * [`regu`] — `REGU(istep, Scur, Splits)` (`cxx:927-1132`) and `REGU()`
//!   (`cxx:1136-1197`).
//!
//! `Standard_Failure("TopOpeBRepTool_REGUW : NO INIT")` (`cxx:159`, `:180`,
//! `:190`, `:281`, ...) becomes an assertion: the state machine of the port
//! always calls `Init` first, and a silent `false` here would hide the very
//! misuse OCCT reports.

use std::collections::{HashMap, HashSet};

use occt_core::gp::{GpPnt2d, GpVec2d};

use crate::shape::TopoShape;
use crate::topo_tools::{shape_map_key, ShapeMapKey};

use super::closed::{FORWARD, REVERSED};
use super::connexity::Connexity;
use super::corr_iso::CorrIso;

mod block;
mod map_s;
mod regu;

/// `class TopOpeBRepTool_REGUW` (`hxx:39-123`).
pub struct Reguw {
    /// `myS` (`hxx:104`); `None` is a null `TopoDS_Shape`.
    my_s: Option<TopoShape>,
    /// `myCORRISO` (`hxx:105`).
    my_corr_iso: CorrIso,
    /// `hasnewsplits` (`hxx:106`).
    hasnew_splits: bool,
    /// `myEsplits` (`hxx:107-108`): edge -> its splits.
    my_e_splits: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myOwNw` (`hxx:109`): old wire -> its new wires.
    my_ow_nw: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `mymapvEds` (`hxx:110-111`): vertex -> its connexity. OCCT's
    /// `NCollection_IndexedDataMap` keeps insertion order, so the order is
    /// kept separately (`MapS` fills it, `InitBlock` walks it).
    my_map_v_eds: HashMap<ShapeMapKey, Connexity>,
    /// Insertion order of `my_map_v_eds`.
    my_v_order: Vec<ShapeMapKey>,
    /// `mymapvmultiple` (`hxx:112`).
    my_map_v_multiple: HashSet<ShapeMapKey>,
    /// `myListVmultiple` (`hxx:113`).
    my_list_v_multiple: Vec<TopoShape>,
    /// `iStep` (`hxx:114`).
    i_step: i32,
    /// `mytol2d` (`hxx:115`); set by `REGU` (`cxx:983`).
    my_tol2d: f64,
    /// `isinit0` (`hxx:116`); set by `REGU` (`cxx:984`).
    is_init0: bool,
    /// `myv0` (`hxx:117`).
    my_v0: Option<TopoShape>,
    /// `myp2d0` (`hxx:118`).
    my_p2d0: GpPnt2d,
    /// `myv` (`hxx:119`).
    my_v: Option<TopoShape>,
    /// `myed` (`hxx:120`).
    my_ed: Option<TopoShape>,
    /// `myp2d` (`hxx:121`).
    my_p2d: GpPnt2d,
    /// `mytg2d` (`hxx:122`); OCCT stores a `gp_Dir2d`, the port the equivalent
    /// normalised `gp_Vec2d` (`gp_Vec2d(gp_Dir2d)` is how it reaches
    /// `Matter`).
    my_tg2d: GpVec2d,
}

impl Reguw {
    /// `TopOpeBRepTool_REGUW(Fref)` (`cxx:119-134`).
    ///
    /// The members `mytol2d`, `isinit0`, `myp2d0`, `myp2d`, `mytg2d` and the
    /// two vertex members are left uninitialised by OCCT's constructor; `REGU`
    /// (`cxx:983-984`, `:1001`) and `InitBlock` (`cxx:537-538`, `:558`) always
    /// set them before a read, so the port starts them at zero.
    pub fn new(f_ref: &TopoShape) -> Self {
        Self {
            my_s: None,
            my_corr_iso: CorrIso::new(f_ref),
            hasnew_splits: false,
            my_e_splits: HashMap::new(),
            my_ow_nw: HashMap::new(),
            my_map_v_eds: HashMap::new(),
            my_v_order: Vec::new(),
            my_map_v_multiple: HashSet::new(),
            my_list_v_multiple: Vec::new(),
            i_step: 0,
            my_tol2d: 0.0,
            is_init0: false,
            my_v0: None,
            my_p2d0: GpPnt2d::new(0.0, 0.0),
            my_v: None,
            my_ed: None,
            my_p2d: GpPnt2d::new(0.0, 0.0),
            my_tg2d: GpVec2d::new(0.0, 0.0),
        }
    }

    /// `Fref()` (`cxx:138-141`).
    pub fn f_ref(&self) -> &TopoShape {
        self.my_corr_iso.f_ref()
    }

    /// `SetEsplits(Esplits)` (`cxx:145-150`).
    pub fn set_e_splits(&mut self, e_splits: &HashMap<ShapeMapKey, Vec<TopoShape>>) {
        self.my_e_splits = e_splits.clone();
    }

    /// `GetEsplits(Esplits)` (`cxx:154-163`).
    pub fn e_splits(&self) -> &HashMap<ShapeMapKey, Vec<TopoShape>> {
        assert!(self.has_init(), "TopOpeBRepTool_REGUW : NO INIT");
        &self.my_e_splits
    }

    /// `SetOwNw(OwNw)` (`cxx:167-171`).
    pub fn set_ow_nw(&mut self, ow_nw: &HashMap<ShapeMapKey, Vec<TopoShape>>) {
        self.my_ow_nw = ow_nw.clone();
    }

    /// `GetOwNw(OwNw)` (`cxx:175-185`).
    pub fn ow_nw(&self) -> &HashMap<ShapeMapKey, Vec<TopoShape>> {
        assert!(self.has_init(), "TopOpeBRepTool_REGUW : NO INIT");
        &self.my_ow_nw
    }

    /// `S()` (`cxx:280-287`).
    pub fn s(&self) -> Option<&TopoShape> {
        assert!(self.has_init(), "TopOpeBRepTool_REGUW : NO INIT");
        self.my_s.as_ref()
    }

    /// `Init(S)` (`cxx:290-294`).
    pub fn init(&mut self, s: &TopoShape) {
        self.my_s = Some(s.clone());
        self.init_step(s);
    }

    /// `HasInit()` (`cxx:298-301`).
    pub fn has_init(&self) -> bool {
        self.my_s.is_some()
    }

    /// `InitStep(S)` (`cxx:305-318`).
    ///
    /// `cxx:311-312` calls `myCORRISO.Init(null)` and then `myCORRISO.Init(S)`.
    /// The null call only clears `myERep2d` / `myEds` / `myVEds`
    /// (`TopOpeBRepTool_CORRISO.cxx:130-137`), which the second call clears
    /// again (`cxx:131-133`), so the port issues the surviving call only.
    pub fn init_step(&mut self, s: &TopoShape) {
        assert!(self.has_init(), "TopOpeBRepTool_REGUW : NO INIT");
        self.my_corr_iso.init(s);

        self.my_map_v_eds.clear();
        self.my_v_order.clear();
        self.my_map_v_multiple.clear();
        self.my_list_v_multiple.clear();
    }

    /// `GetSplits(Splits)` (`cxx:1200-1214`).
    pub fn get_splits(&self, splits: &mut Vec<TopoShape>) -> bool {
        assert!(self.has_init(), "TopOpeBRepTool_REGUW : NO INIT");
        let Some(s) = self.my_s.as_ref() else {
            return false;
        };
        match self.my_ow_nw.get(&shape_map_key(s)) {
            Some(l) => {
                *splits = l.clone();
                true
            }
            None => false,
        }
    }

    /// `Connexity(v, co)` (`cxx:1218-1231`).
    pub fn connexity(&self, v: &TopoShape) -> Option<Connexity> {
        assert!(self.has_init(), "TopOpeBRepTool_REGUW : NO INIT");
        self.my_map_v_eds.get(&shape_map_key(v)).cloned()
    }

    /// `AddNewConnexity(v, OriKey, e)` (`cxx:1235-1276`).
    pub fn add_new_connexity(&mut self, v: &TopoShape, ori_key: i32, e: &TopoShape) -> bool {
        assert!(self.has_init(), "TopOpeBRepTool_REGUW : NO INIT");
        let key = shape_map_key(v);
        if !self.my_map_v_eds.contains_key(&key) {
            return false;
        }

        if !self.my_corr_iso.add_new_connexity(e) {
            return false;
        }

        if let Some(co) = self.my_map_v_eds.get_mut(&key) {
            co.add_item(ori_key, e);
        }
        true
    }

    /// `RemoveOldConnexity(v, OriKey, e)` (`cxx:1281-1300`).
    ///
    /// `cxx:1296` comments out the `!ok` early return of the CORRISO call; only
    /// the connexity removal decides the answer.
    pub fn remove_old_connexity(&mut self, v: &TopoShape, ori_key: i32, e: &TopoShape) -> bool {
        assert!(self.has_init(), "TopOpeBRepTool_REGUW : NO INIT");
        let key = shape_map_key(v);
        if !self.my_map_v_eds.contains_key(&key) {
            return false;
        }

        let _ = self.my_corr_iso.remove_old_connexity(e);

        match self.my_map_v_eds.get_mut(&key) {
            Some(co) => co.remove_item_ori(ori_key, e),
            None => false,
        }
    }

    /// `UpdateMultiple(v)` (`cxx:1304-1325`).
    pub fn update_multiple(&mut self, v: &TopoShape) -> bool {
        assert!(self.has_init(), "TopOpeBRepTool_REGUW : NO INIT");
        let key = shape_map_key(v);
        let Some(co) = self.my_map_v_eds.get(&key) else {
            return false;
        };
        if co.is_multiple() {
            if self.my_map_v_multiple.insert(key) {
                self.my_list_v_multiple.push(v.clone());
            }
        }
        true
    }

    /// `(iStep == 1) ? REVERSED : FORWARD` (`cxx:505`, `:539`, `:691`, `:805`).
    fn iv0e1(&self) -> i32 {
        if self.i_step == 1 {
            REVERSED
        } else {
            FORWARD
        }
    }
}
