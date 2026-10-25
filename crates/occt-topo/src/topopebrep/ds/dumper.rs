//! `TopOpeBRepDS_Dumper` — printable names for the shapes of a
//! `TopOpeBRepDS_HDataStructure`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_Dumper.hxx:35-58`,
//! `.cxx:35-150`.
//!
//! | OCCT | lines |
//! |---|---|
//! | constructor | `cxx:35-38` |
//! | `SDumpRefOri(K, I)` | `cxx:42-61` |
//! | `SDumpRefOri(S)` | `cxx:65-73` |
//! | `SPrintShape(I)` | `cxx:77-87` |
//! | `SPrintShape(S)` | `cxx:92-114` |
//! | `SPrintShapeRefOri(S, astr)` | `cxx:116-123` |
//! | `SPrintShapeRefOri(L, astr)` | `cxx:125-150` |
//!
//! `TCollection_AsciiString` becomes `String`; every method returns a fresh
//! one, so the OCCT `const` methods stay `&self`.
//!
//! The reference/orientation pair prints as `(rank, CONFIG)` (`cxx:59`), where
//! `rank` is `SameDomainReference(S)` and `CONFIG` is `SPrint(Config)` of
//! `SameDomainOrientation(S)`. `IsTopology(K)` gates the whole thing
//! (`cxx:47-50`): a geometry kind yields the empty string.
//!
//! `SPrintShape(S)` (`cxx:92-114`) marks the ancestor rank with `*` when the
//! shape is kept and `~` when it is not, in front of the rank-1 and rank-2
//! slots of `TopOpeBRepDS::SPrint(K, I, s1, s2)`.
//!
//! OCCT indexes the data structure through accessors whose `FindKeep`
//! argument defaults to **true** (`TopOpeBRepDS_DataStructure.hxx`), so
//! `BDS.Shape(IS)` / `DS().Shape(S)` / `BDS.KeepShape(IS)` here read
//! `find_keep = true`, while the two `SDumpRefOri` entry points pass the
//! explicit `fk = false` of `cxx:47` / `cxx:69`.

use crate::shape::TopoShape;

use super::h_data_structure::HDataStructure;
use super::kind::{
    config_s_print, is_topology, kind_index_s_print, kind_to_shape, shape_to_kind, Config, Kind,
};

/// `TopOpeBRepDS_Dumper` (`hxx:35-58`).
#[derive(Clone, Copy)]
pub struct Dumper<'a> {
    /// `myHDS` (`hxx:56`) — `None` is OCCT's null handle.
    hds: Option<&'a HDataStructure>,
}

impl<'a> Dumper<'a> {
    /// Constructor (`cxx:35-38`).
    pub fn new(hds: &'a HDataStructure) -> Self {
        Self { hds: Some(hds) }
    }

    /// A dumper over the null handle (`cxx:37` assigns without checking).
    pub fn null() -> Self {
        Self { hds: None }
    }

    /// `SDumpRefOri(K, I)` (`cxx:42-61`).
    pub fn sdump_ref_ori_kind(&self, k: Kind, i: i32) -> String {
        let Some(hds) = self.hds else {
            return String::new();
        };
        // `cxx:47`
        if !is_topology(k) {
            return String::new();
        }
        let ds = hds.ds();
        let fk = false;
        // `cxx:51-55`
        let t = kind_to_shape(k);
        match ds.shape(i, fk) {
            Some(s) if s.shape_type() == t => {}
            _ => return String::new(),
        }
        // `cxx:56-59`
        let Some(s) = hds.shape(i, fk) else {
            return String::new();
        };
        let r = hds.same_domain_reference(&s);
        let o: Config = hds.same_domain_orientation(&s);
        format!("({},{})", self.sprint_shape_index(r), config_s_print(o))
    }

    /// `SDumpRefOri(S)` (`cxx:65-73`).
    ///
    /// `ShapeToKind` raises `Standard_ProgramError` for a `TopAbs_SHAPE`
    /// (`TopOpeBRepDS.cxx:152-183`); the empty string is that case here, the
    /// same value [`Self::sdump_ref_ori_kind`] returns for a non-topological
    /// kind.
    pub fn sdump_ref_ori(&self, s: &TopoShape) -> String {
        let Some(k) = shape_to_kind(s.shape_type()).ok() else {
            return String::new();
        };
        let fk = false;
        let i = match self.hds {
            Some(hds) => hds.shape_index(s, fk),
            None => 0,
        };
        self.sdump_ref_ori_kind(k, i)
    }

    /// `SPrintShape(I)` (`cxx:77-87`): empty for `I` outside `[1, NbShapes]`.
    pub fn sprint_shape_index(&self, is: i32) -> String {
        let Some(hds) = self.hds else {
            return String::new();
        };
        let bds = hds.ds();
        // `cxx:82-85`
        if is < 1 || is > bds.nb_shapes() {
            return String::new();
        }
        // `cxx:86`: `BDS.Shape(IS)` — `FindKeep` defaults to true.
        match bds.shape(is, true) {
            Some(s) => self.sprint_shape(&s),
            None => String::new(),
        }
    }

    /// `SPrintShape(S)` (`cxx:92-114`).
    pub fn sprint_shape(&self, s: &TopoShape) -> String {
        let Some(hds) = self.hds else {
            return String::new();
        };
        let bds = hds.ds();
        // `cxx:93-94`: `DS().Shape(S)` — `FindKeep` defaults to true.
        let is = bds.shape_index(s, true);
        let rank_is = bds.ancestor_rank(is);

        // `cxx:99-108`
        let (s1, s2) = if bds.keep_shape(is, true) {
            (
                if rank_is == 1 { "*" } else { "" },
                if rank_is == 2 { "*" } else { "" },
            )
        } else {
            (
                if rank_is == 1 { "~" } else { "" },
                if rank_is == 2 { "~" } else { "" },
            )
        };

        // `cxx:109-111`
        let Some(k) = shape_to_kind(s.shape_type()).ok() else {
            return String::new();
        };
        kind_index_s_print(k, is, s1, s2)
    }

    /// `SPrintShapeRefOri(S, astr)` (`cxx:116-123`).
    pub fn sprint_shape_ref_ori(&self, s: &TopoShape, astr: &str) -> String {
        format!("{}{} {}", astr, self.sprint_shape(s), self.sdump_ref_ori(s))
    }

    /// `SPrintShapeRefOri(L, astr)` (`cxx:125-150`): one line per shape, the
    /// continuation lines indented by `astr.Length()` blanks.
    pub fn sprint_shape_ref_ori_list(&self, list: &[TopoShape], astr: &str) -> String {
        // `cxx:131-135`
        if list.is_empty() {
            return String::new();
        }

        // `cxx:136-137`
        let mut ss = astr.to_string();
        let bst = " ".repeat(astr.len());

        // `cxx:138-147`
        for (il, s) in list.iter().enumerate() {
            let mut v = self.sprint_shape_ref_ori(s, "");
            if il != 0 {
                v = format!("{}{}", bst, v);
            }
            ss = format!("{}{}\n", ss, v);
        }
        ss
    }
}
