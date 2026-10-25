//! `TopOpeBRepDS_TKI` — the "table of kind-indexed interference lists" that
//! groups a list of interferences by `(Kind, index)` and iterates the
//! non-empty groups.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_TKI.hxx:33-122`,
//! `TopOpeBRepDS_TKI.cxx:1-432`.
//!
//! `myT` is an `HArray1<DataMap<int, List>>` indexed by `Kind + mydelta`, i.e.
//! one table per `Kind` in `[POINT, SOLID]` (`cxx:43-62`). The port keeps the
//! same 1-based table index in `ti` and stores the tables in a `Vec` of
//! `BTreeMap`; OCCT's `NCollection_DataMap` iteration order follows its hash
//! buckets, so the port picks the only deterministic order available (ascending
//! key) and documents that instead of reproducing the hash order.
//!
//! `myITM` is a `DataMap::Iterator` over the current table. The port replaces
//! it with a snapshot of the table's keys plus a cursor: the iteration never
//! inserts into the map (`ChangeInterferences` only reaches an existing
//! binding, see `cxx:139-150`), so the snapshot stays valid.

use std::collections::BTreeMap;

use super::interference::InterferenceHandle;
use super::interference_iterator::InterferenceList;
use super::kind::Kind;

/// `TopOpeBRepDS_TKI` (`hxx:33-122`).
#[derive(Clone, Debug)]
pub struct TKI {
    /// `mydelta` (`cxx:56`): `1 - POINT`, so `Kind + mydelta` is `>= 1`.
    delta: i32,
    /// `myT`: one `DataMap<int, List>` per `Kind` in `[POINT, SOLID]`.
    tables: Vec<BTreeMap<i32, InterferenceList>>,
    /// `myTI` (`cxx:114`): 1-based table index.
    ti: i32,
    /// `myG` (`cxx:115`).
    g: i32,
    /// `myK` (`cxx:119`).
    k: Kind,
    /// `myITM` position: the keys of `tables[ti]` in iteration order.
    keys: Vec<i32>,
    key_pos: usize,
    /// `myEmptyLOI` (`cxx:121`): the list returned when the query misses.
    empty: InterferenceList,
}

impl Default for TKI {
    fn default() -> Self {
        Self::new()
    }
}

impl TKI {
    /// `TopOpeBRepDS_TKI()` (`cxx:37-40`) + `Reset()` (`cxx:44-64`).
    pub fn new() -> Self {
        let mut t = TKI {
            delta: 0,
            tables: Vec::new(),
            ti: 0,
            g: 0,
            k: Kind::Unknown,
            keys: Vec::new(),
            key_pos: 0,
            empty: InterferenceList::new(),
        };
        t.reset();
        t
    }

    /// `Reset()` (`cxx:44-64`): `f = 1`, `l = 1 + (SOLID - POINT)`,
    /// `mydelta = f - POINT`; the tables are allocated once.
    fn reset(&mut self) {
        let ip = Kind::Point as i32;
        let is = Kind::Solid as i32;
        // cxx:49-52: the enumeration must be ordered POINT .. SOLID.
        debug_assert!(ip <= is, "TopOpeBRepDS_TKI : enumeration badly ordered");
        let f = 1;
        let l = f + (is - ip);
        self.delta = f - ip;
        if self.tables.len() != (l - f + 1) as usize {
            self.tables = vec![BTreeMap::new(); (l - f + 1) as usize];
        }
        self.clear();
        self.k = Kind::Unknown;
        self.g = 0;
    }

    /// `Clear()` (`cxx:68-75`).
    pub fn clear(&mut self) {
        for t in self.tables.iter_mut() {
            t.clear();
        }
    }

    /// `FillOnGeometry(L)` (`cxx:78-91`): group by `(GeometryType, Geometry)`.
    pub fn fill_on_geometry(&mut self, l: &InterferenceList) {
        for i in l {
            let h = i.read().expect("interference lock");
            let (gt, g) = (h.geometry_type(), h.geometry());
            drop(h);
            // FDS_data leaves the outputs untouched on a null handle; a
            // non-null handle always has both fields.
            let _ = self.add(gt, g);
            self.change_interferences(gt, g).map(|v| v.push(i.clone()));
        }
    }

    /// `FillOnSupport(L)` (`cxx:94-108`): group by `(SupportType, Support)`.
    pub fn fill_on_support(&mut self, l: &InterferenceList) {
        for i in l {
            let h = i.read().expect("interference lock");
            let (st, s) = (h.support_type(), h.support());
            drop(h);
            let _ = self.add(st, s);
            self.change_interferences(st, s).map(|v| v.push(i.clone()));
        }
    }

    /// `IsBound(K, G)` (`cxx:112-120`).
    pub fn is_bound(&self, k: Kind, g: i32) -> bool {
        if !self.is_valid_kg(k, g) {
            return false;
        }
        let ti = self.kind_to_table_index(k);
        self.table(ti).is_some_and(|t| t.contains_key(&g))
    }

    /// `Interferences(K, G)` (`cxx:124-135`).
    pub fn interferences(&self, k: Kind, g: i32) -> &InterferenceList {
        if !self.is_bound(k, g) {
            return &self.empty;
        }
        let ti = self.kind_to_table_index(k);
        match self.table(ti).and_then(|t| t.get(&g)) {
            Some(v) => v,
            None => &self.empty,
        }
    }

    /// `ChangeInterferences(K, G)` (`cxx:139-150`).
    ///
    /// OCCT returns the member `myEmptyLOI` when the pair is not bound, so a
    /// write through the result is discarded; `None` is that same discard.
    pub fn change_interferences(&mut self, k: Kind, g: i32) -> Option<&mut InterferenceList> {
        if !self.is_bound(k, g) {
            return None;
        }
        let ti = self.kind_to_table_index(k);
        self.table_mut(ti).and_then(|t| t.get_mut(&g))
    }

    /// `HasInterferences(K, G)` (`cxx:154-164`).
    pub fn has_interferences(&self, k: Kind, g: i32) -> bool {
        if !self.is_bound(k, g) {
            return false;
        }
        !self.interferences(k, g).is_empty()
    }

    /// `Add(K, G)` (`cxx:168-183`): binds an empty list.
    pub fn add(&mut self, k: Kind, g: i32) -> Result<(), String> {
        if !self.is_valid_kg(k, g) {
            return Err("TopOpeBRepDS_TKI : Add K G".into());
        }
        let ti = self.kind_to_table_index(k);
        if let Some(t) = self.table_mut(ti) {
            if !t.contains_key(&g) {
                t.insert(g, InterferenceList::new());
            }
        }
        Ok(())
    }

    /// `Add(K, G, HI)` (`cxx:187-201`).
    pub fn add_interference(
        &mut self,
        k: Kind,
        g: i32,
        hi: &InterferenceHandle,
    ) -> Result<(), String> {
        if !self.is_valid_kg(k, g) {
            return Err("TopOpeBRepDS_TKI : Add K G HI".into());
        }
        self.add(k, g)?;
        if let Some(v) = self.change_interferences(k, g) {
            v.push(hi.clone());
        }
        Ok(())
    }

    /// `DumpTKIIterator(s1, s2)` (`cxx:205-218`): the C++ body only prints the
    /// bounds and walks the iterator; the port keeps the walk so the call is
    /// not a no-op, writing the same `K`/`G` pairs to stderr.
    pub fn dump_tki_iterator(&mut self, s1: &str, s2: &str) {
        eprint!("{s1}");
        self.init();
        while self.more() {
            let (k, g) = (self.k, self.g);
            eprintln!("{:?} {}", k, g);
            self.next();
        }
        eprint!("{s2}");
    }

    /// `Init()` (`cxx:222-233`).
    pub fn init(&mut self) {
        self.k = Kind::Unknown;
        self.g = 0;
        self.ti = 1;
        self.k = self.table_index_to_kind(self.ti);
        self.init_itm();
        self.find();
    }

    /// `More()` (`cxx:237-241`).
    pub fn more(&self) -> bool {
        self.is_valid_kg(self.k, self.g)
    }

    /// `Next()` (`cxx:245-261`).
    pub fn next(&mut self) {
        if self.more_itm() {
            self.next_itm();
            self.find();
        } else if self.more_ti() {
            self.next_ti();
            if self.more_ti() {
                self.init_itm();
            }
            self.find();
        }
    }

    /// `Value(K, G)` (`cxx:265-274`).
    pub fn value(&self) -> Option<(Kind, i32, &InterferenceList)> {
        if !self.more() {
            return None;
        }
        Some((self.k, self.g, self.interferences(self.k, self.g)))
    }

    /// `ChangeValue(K, G)` (`cxx:283-295`).
    pub fn change_value(&mut self) -> Option<(Kind, i32, &mut InterferenceList)> {
        if !self.more() {
            return None;
        }
        let (k, g) = (self.k, self.g);
        self.change_interferences(k, g).map(|v| (k, g, v))
    }

    /// `KindToTableIndex(K)` (`cxx:383-388`).
    fn kind_to_table_index(&self, k: Kind) -> i32 {
        k as i32 + self.delta
    }

    /// `TableIndexToKind(TI)` (`cxx:392-397`).
    fn table_index_to_kind(&self, ti: i32) -> Kind {
        kind_from_index(ti - self.delta)
    }

    /// `IsValidTI(TI)` (`cxx:402-409`).
    fn is_valid_ti(&self, ti: i32) -> bool {
        let f = 1;
        let l = self.tables.len() as i32;
        ti >= f && ti <= l
    }

    /// `IsValidK(K)` (`cxx:412-417`): `POINT ..= SOLID` in enum order.
    fn is_valid_k(&self, k: Kind) -> bool {
        let (lo, hi) = (Kind::Point as i32, Kind::Solid as i32);
        let v = k as i32;
        v >= lo && v <= hi
    }

    /// `IsValidG(G)` (`cxx:421-425`).
    fn is_valid_g(&self, g: i32) -> bool {
        g > 0
    }

    /// `IsValidKG(K, G)` (`cxx:428-432`).
    fn is_valid_kg(&self, k: Kind, g: i32) -> bool {
        self.is_valid_k(k) && self.is_valid_g(g)
    }

    /// `myT->Value(TI)`.
    fn table(&self, ti: i32) -> Option<&BTreeMap<i32, InterferenceList>> {
        if !self.is_valid_ti(ti) {
            return None;
        }
        self.tables.get((ti - 1) as usize)
    }

    /// `myT->ChangeValue(TI)`.
    fn table_mut(&mut self, ti: i32) -> Option<&mut BTreeMap<i32, InterferenceList>> {
        if !self.is_valid_ti(ti) {
            return None;
        }
        self.tables.get_mut((ti - 1) as usize)
    }

    /// `myITM.Initialize(myT->Value(myTI))`.
    fn init_itm(&mut self) {
        self.keys = match self.table(self.ti) {
            Some(t) => t.keys().copied().collect(),
            None => Vec::new(),
        };
        self.key_pos = 0;
    }

    /// `MoreTI()` (`cxx:297-300`).
    fn more_ti(&self) -> bool {
        self.is_valid_ti(self.ti)
    }

    /// `NextTI()` (`cxx:305-308`).
    fn next_ti(&mut self) {
        self.ti += 1;
        self.k = self.table_index_to_kind(self.ti);
    }

    /// `MoreITM()` (`cxx:312-316`).
    fn more_itm(&self) -> bool {
        self.key_pos < self.keys.len()
    }

    /// `FindITM()` (`cxx:321-337`).
    fn find_itm(&mut self) {
        while self.more_itm() {
            self.g = self.keys[self.key_pos];
            if self.has_interferences(self.k, self.g) {
                break;
            } else {
                self.key_pos += 1;
            }
        }
    }

    /// `NextITM()` (`cxx:342-348`).
    fn next_itm(&mut self) {
        if !self.more_itm() {
            return;
        }
        self.key_pos += 1;
        self.find_itm();
    }

    /// `Find()` (`cxx:353-379`): advance to the first non-empty group,
    /// walking the tables in ascending index order.
    fn find(&mut self) {
        let mut f = false;
        while self.more_ti() {
            while self.more_itm() {
                self.find_itm();
                f = self.has_interferences(self.k, self.g);
                if f {
                    break;
                }
            }
            if f {
                break;
            } else {
                self.next_ti();
                if self.more_ti() {
                    self.init_itm();
                }
            }
        }
    }
}

/// `(TopOpeBRepDS_Kind)v` restricted to the ordered enumeration
/// `POINT ..= SOLID`; anything else is `UNKNOWN`.
fn kind_from_index(v: i32) -> Kind {
    match v {
        0 => Kind::Point,
        1 => Kind::Curve,
        2 => Kind::Surface,
        3 => Kind::Vertex,
        4 => Kind::Edge,
        5 => Kind::Wire,
        6 => Kind::Face,
        7 => Kind::Shell,
        8 => Kind::Solid,
        _ => Kind::Unknown,
    }
}
