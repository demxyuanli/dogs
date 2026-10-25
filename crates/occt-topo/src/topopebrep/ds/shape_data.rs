//! `TopOpeBRepDS_ShapeData` — the per-shape record of `TopOpeBRepDS_DataStructure`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_ShapeData.hxx:24-62`,
//! `TopOpeBRepDS_ShapeData.cxx:20-57`.
//!
//! OCCT's fields are private and only `TopOpeBRepDS_DataStructure` (a `friend
//! class`) touches them besides `Interferences`/`Keep`; the Rust fields are
//! therefore `pub(super)` to `topopebrep::ds`.

use crate::abs::Orientation;
use crate::shape::TopoShape;

use super::interference::InterferenceHandle;
use super::kind::Config;

/// Source: `TopOpeBRepDS_ShapeData.hxx:48-60`
#[derive(Debug, Clone)]
pub struct ShapeData {
    /// `myInterferences`
    pub(super) interferences: Vec<InterferenceHandle>,
    /// `mySameDomain`
    pub(super) same_domain: Vec<TopoShape>,
    /// `mySameDomainRef`
    pub(super) same_domain_ref: i32,
    /// `mySameDomainOri`
    pub(super) same_domain_ori: Config,
    /// `mySameDomainInd`
    pub(super) same_domain_ind: i32,
    /// `myOrientation`
    pub(super) orientation: Orientation,
    /// `myOrientationDef`
    pub(super) orientation_def: bool,
    /// `myAncestorRank`
    pub(super) ancestor_rank: i32,
    /// `myKeep`
    pub(super) keep: bool,
}

impl Default for ShapeData {
    /// `TopOpeBRepDS_ShapeData()` (`cxx:22-31`).
    fn default() -> Self {
        Self {
            interferences: Vec::new(),
            same_domain: Vec::new(),
            same_domain_ref: 0,
            same_domain_ori: Config::UnshGeometry,
            same_domain_ind: 0,
            orientation: Orientation::Forward,
            orientation_def: false,
            ancestor_rank: 0,
            keep: true,
        }
    }
}

impl ShapeData {
    /// `TopOpeBRepDS_ShapeData()` (`cxx:22-31`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Interferences()` (`cxx:35-39`).
    pub fn interferences(&self) -> &[InterferenceHandle] {
        &self.interferences
    }

    /// `ChangeInterferences()` (`cxx:43-47`).
    pub fn change_interferences(&mut self) -> &mut Vec<InterferenceHandle> {
        &mut self.interferences
    }

    /// `Keep()` (`cxx:51-54`).
    pub fn keep(&self) -> bool {
        self.keep
    }

    /// `ChangeKeep(B)` (`cxx:56-59`).
    pub fn change_keep(&mut self, b: bool) {
        self.keep = b;
    }

    /// `mySameDomain` / `mySameDomainRef` / `mySameDomainOri` /
    /// `mySameDomainInd` / `myOrientation` / `myOrientationDef` /
    /// `myAncestorRank` have no OCCT accessors: `TopOpeBRepDS_DataStructure`
    /// (a `friend class`, `hxx:311`) reads and writes the fields directly.
    /// These accessors are that friendship, one field per pair.
    pub fn same_domain(&self) -> &[TopoShape] {
        &self.same_domain
    }

    pub fn change_same_domain(&mut self) -> &mut Vec<TopoShape> {
        &mut self.same_domain
    }

    pub fn same_domain_ref(&self) -> i32 {
        self.same_domain_ref
    }

    pub fn set_same_domain_ref(&mut self, r: i32) {
        self.same_domain_ref = r;
    }

    pub fn same_domain_ori(&self) -> Config {
        self.same_domain_ori
    }

    pub fn set_same_domain_ori(&mut self, o: Config) {
        self.same_domain_ori = o;
    }

    pub fn same_domain_ind(&self) -> i32 {
        self.same_domain_ind
    }

    pub fn set_same_domain_ind(&mut self, i: i32) {
        self.same_domain_ind = i;
    }

    pub fn orientation(&self) -> Orientation {
        self.orientation
    }

    pub fn set_orientation(&mut self, o: Orientation) {
        self.orientation = o;
    }

    pub fn orientation_def(&self) -> bool {
        self.orientation_def
    }

    pub fn set_orientation_def(&mut self, b: bool) {
        self.orientation_def = b;
    }

    pub fn ancestor_rank(&self) -> i32 {
        self.ancestor_rank
    }

    pub fn set_ancestor_rank(&mut self, i: i32) {
        self.ancestor_rank = i;
    }
}
