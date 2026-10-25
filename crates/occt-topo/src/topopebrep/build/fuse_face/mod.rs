//! `TopOpeBRepBuild_FuseFace` — fuse co-surfacic faces connected across a
//! shared edge.
//!
//! Source: `TKBool/TopOpeBRepBuild/TopOpeBRepBuild_FuseFace.hxx:29-85`,
//! `.cxx:49-1398`, `.lxx:20-96` (OCCT 8.0.0). Split by subject:
//!
//! | OCCT | lines | here |
//! |---|---|---|
//! | constructors / `Init` / accessors | `hxx:33-68`, `lxx:20-96`, `cxx:73-118` | this file |
//! | `PerformFace` | `cxx:126-582` | [`perform_face`] |
//! | `PerformEdge` | `cxx:590-681` | [`perform_edge`] |
//! | `ClearEdge` | `cxx:688-822` | [`perform_edge`] |
//! | `ClearVertex` | `cxx:829-841` | [`perform_edge`] |
//! | `GroupShape` | `cxx:845-971` | [`group_shape`] |
//! | `GroupEdge` | `cxx:975-1074` | [`group_edge`] |
//! | `MakeEdge` | `cxx:1078-1177` | [`make_edge`] |
//! | `SameSupport` | `cxx:1181-1398` | [`same_support`] |
//!
//! The class is the fusion step `GFillFaceSFS` runs on `GLOBAL_lfr1`
//! (`ffsfs.cxx:618-667`): the splits `GSplitFaceSFS` produced for one
//! same-domain family are fused back into one face per co-surfacic family, and
//! the survivors become start elements of the `ShellFaceSet`.
//!
//! UNPORTED in this module:
//! * every `OCCT_DEBUG` block — the 32 `#ifdef OCCT_DEBUG` lines of
//!   `FuseFace.cxx` (53, 77, 92, 128, 138, 152, 188, 205, 223, 236, 276, 316,
//!   346, 362, 438, 539, 564, 578, 592, 677, 690, 743, 766, 788, 803, 818, 831,
//!   837, 1082, 1141, 1183, 1234) — tracing only, gated by
//!   `TopOpeBRepBuild_GettraceFUFA()`.

use crate::shape::TopoShape;

pub mod group_edge;
pub mod group_shape;
pub mod make_edge;
pub mod perform_edge;
pub mod perform_face;
pub mod same_support;
pub mod shape_list_map;

/// `class TopOpeBRepBuild_FuseFace` (`hxx:28-85`).
#[derive(Debug, Clone)]
pub struct FuseFace {
    /// `myLIF` (`hxx:79`) — the list of internal faces, appended by no method
    /// of the class.
    my_lif: Vec<TopoShape>,
    /// `myLRF` (`hxx:80`) — the reconstructed faces to fuse.
    my_lrf: Vec<TopoShape>,
    /// `myLFF` (`hxx:81`) — the fused faces.
    my_lff: Vec<TopoShape>,
    /// `myLIE` (`hxx:71`) — the internal edges collected by `ClearEdge`.
    my_lie: Vec<TopoShape>,
    /// `myLEE` (`hxx:72`) — the external edges collected by `ClearEdge`.
    my_lee: Vec<TopoShape>,
    /// `myLME` (`hxx:73`) — the modified edges: an edge seen both FORWARD and
    /// REVERSED (`cxx:394-397`), collected only when `myInternal` is set.
    my_lme: Vec<TopoShape>,
    /// `myLIV` (`hxx:74`) — never written in OCCT.
    my_liv: Vec<TopoShape>,
    /// `myLEV` (`hxx:75`) — never written in OCCT.
    my_lev: Vec<TopoShape>,
    /// `myLMV` (`hxx:76`) — never written in OCCT.
    my_lmv: Vec<TopoShape>,
    /// `myInternal` (`hxx:82`) — `CXM == 2` in `Init`: keep the internal
    /// connections.
    my_internal: bool,
    /// `myModified` (`hxx:83`).
    my_modified: bool,
    /// `myDone` (`hxx:84`).
    my_done: bool,
}

impl Default for FuseFace {
    fn default() -> Self {
        Self::new()
    }
}

impl FuseFace {
    /// `TopOpeBRepBuild_FuseFace()` (`lxx:21-26`).
    pub fn new() -> Self {
        Self {
            my_lif: Vec::new(),
            my_lrf: Vec::new(),
            my_lff: Vec::new(),
            my_lie: Vec::new(),
            my_lee: Vec::new(),
            my_lme: Vec::new(),
            my_liv: Vec::new(),
            my_lev: Vec::new(),
            my_lmv: Vec::new(),
            my_internal: false,
            my_modified: false,
            my_done: false,
        }
    }

    /// `TopOpeBRepBuild_FuseFace(LIF, LRF, CXM)` (`lxx:30-35`).
    pub fn with_lists(lif: &[TopoShape], lrf: &[TopoShape], cxm: i32) -> Self {
        let mut f = Self::new();
        f.init(lif, lrf, cxm);
        f
    }

    /// `Init(LIF, LRF, CXM)` (`cxx:72-117`).
    pub fn init(&mut self, lif: &[TopoShape], lrf: &[TopoShape], cxm: i32) {
        self.my_lif = lif.to_vec();
        self.my_lrf = lrf.to_vec();
        // `cxx:84-92`: only 1 and 2 assign `myInternal`; any other value keeps
        // the default constructor's `false`.
        if cxm == 1 {
            self.my_internal = false;
        } else if cxm == 2 {
            self.my_internal = true;
        }

        // `cxx:107-117`.
        self.my_lff.clear();
        self.my_lie.clear();
        self.my_lee.clear();
        self.my_lme.clear();
        self.my_liv.clear();
        self.my_lev.clear();
        self.my_lmv.clear();
        self.my_modified = false;
        self.my_done = false;
    }

    /// `IsDone()` (`lxx:38-41`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `IsModified()` (`lxx:45-48`).
    pub fn is_modified(&self) -> bool {
        self.my_modified
    }

    /// `LFuseFace()` (`lxx:52-56`).
    pub fn l_fuse_face(&self) -> &[TopoShape] {
        &self.my_lff
    }

    /// `LInternEdge()` (`lxx:60-64`).
    pub fn l_intern_edge(&self) -> &[TopoShape] {
        &self.my_lie
    }

    /// `LExternEdge()` (`lxx:68-71`).
    pub fn l_extern_edge(&self) -> &[TopoShape] {
        &self.my_lee
    }

    /// `LModifEdge()` (`lxx:75-78`).
    pub fn l_modif_edge(&self) -> &[TopoShape] {
        &self.my_lme
    }

    /// `LInternVertex()` (`lxx:82-85`).
    pub fn l_intern_vertex(&self) -> &[TopoShape] {
        &self.my_liv
    }

    /// `LExternVertex()` (`lxx:89-91`).
    pub fn l_extern_vertex(&self) -> &[TopoShape] {
        &self.my_lev
    }

    /// `LModifVertex()` (`lxx:94-96`).
    pub fn l_modif_vertex(&self) -> &[TopoShape] {
        &self.my_lmv
    }

    /// The `myModified = false; myDone = true; myLFF = myLRF; return;` epilogue
    /// repeated by `PerformFace` (`myLFF = myLRF` at `cxx:145`, `:159`, `:195`,
    /// `:211`, `:229`, `:242`, `:282`, `:323`, `:353`, `:369`, `:445`, `:546`,
    /// `:571`) and by `ClearEdge` (`cxx:750`, `:772`, `:794`, `:810`).
    fn no_fuse(&mut self) {
        self.my_modified = false;
        self.my_done = true;
        self.my_lff = self.my_lrf.clone();
    }
}
