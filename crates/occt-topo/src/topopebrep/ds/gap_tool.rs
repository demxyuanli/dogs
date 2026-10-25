//! `TopOpeBRepDS_GapTool` — the point-geometry index `TopOpeBRepDS_GapFiller`
//! queries while pairing interferences.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_GapTool.hxx:34-96`,
//! `TopOpeBRepDS_GapTool.cxx:22-275` (OCCT 8.0.0).
//!
//! The tool builds two reverse indexes in `Init` (`cxx:52-84`):
//! `myInterToShape` (interference -> the shape whose list holds it) and
//! `myGToI` (point index -> every interference on that point). Both are keyed
//! on pointer identity in OCCT; the port keys on [`std::sync::Arc::ptr_eq`] and
//! on the point index, keeping insertion order because the `NCollection_List`
//! answers of `Interferences` / `SameInterferences` are walked in order.
//!
//! Two OCCT failure modes become `Option`s here, and both are unreachable from
//! the `GapFiller` call sites:
//!
//! * `myGToI(G)` (`cxx:102`, `:157`, `:209`, `:246`, `:267`) is
//!   `NCollection_DataMap::operator()`, which throws when `G` is unbound.
//!   [`GapTool::same_interferences`] and friends answer the empty list.
//! * `occ::down_cast<TopOpeBRepDS_CurvePointInterference>(I)->Parameter()`
//!   (`cxx:205`, `:219`, `:241`, `:253`) dereferences a null handle when the
//!   interference is another subclass. The port checks
//!   [`Interference::is_curve_point`] first and answers `None` when it fails.

use crate::shape::TopoShape;
use crate::topo_tools::is_same;

use super::curve::Curve;
use super::h_data_structure::HDataStructure;
use super::interference::{Interference, InterferenceHandle};
use super::interference_iterator::InterferenceList;
use super::kind::Kind;

/// `TopOpeBRepDS_GapTool` (`hxx:34-96`).
#[derive(Debug, Default)]
pub struct GapTool {
    /// `myInterToShape` (`hxx:93`): `DataMap<handle, TopoDS_Shape>` on handle
    /// identity.
    inter_to_shape: Vec<(InterferenceHandle, TopoShape)>,
    /// `myGToI` (`hxx:94`): `DataMap<int, List<handle>>`.
    g_to_i: Vec<(i32, InterferenceList)>,
}

impl GapTool {
    /// `TopOpeBRepDS_GapTool()` (`cxx:26`): `= default`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `TopOpeBRepDS_GapTool(HDS)` (`cxx:31-34`).
    pub fn with_hds(hds: &HDataStructure) -> Self {
        let mut t = Self::new();
        t.init(hds);
        t
    }

    /// `StoreGToI(GToI, I)` (`cxx:35-47`): bind an empty list for a new
    /// geometry index, then append.
    fn store_g_to_i(&mut self, i: &InterferenceHandle) {
        let g = match i.read() {
            Ok(itf) => itf.geometry(),
            Err(_) => return,
        };
        match self.g_to_i.iter_mut().find(|(k, _)| *k == g) {
            Some((_, list)) => list.push(i.clone()),
            None => self.g_to_i.push((g, vec![i.clone()])),
        }
    }

    /// `Init(HDS)` (`cxx:52-84`).
    pub fn init(&mut self, hds: &HDataStructure) {
        self.inter_to_shape.clear();
        self.g_to_i.clear();

        // `cxx:55-70`: shape interferences of POINT geometry.
        let nb = hds.ds().nb_shapes();
        for i in 1..=nb {
            let Some(s) = hds.ds().shape(i, true) else {
                continue;
            };
            let li = hds.ds().shape_interferences(&s, true).to_vec();
            for handle in li {
                let is_point = match handle.read() {
                    Ok(itf) => itf.geometry_type() == Kind::Point,
                    Err(_) => false,
                };
                if !is_point {
                    continue;
                }
                // `NCollection_DataMap::Bind` throws on a second bind for the
                // same handle, and each interference lives in one list, so the
                // first association wins here.
                if !self
                    .inter_to_shape
                    .iter()
                    .any(|(k, _)| std::sync::Arc::ptr_eq(k, &handle))
                {
                    self.inter_to_shape.push((handle.clone(), s.clone()));
                }
                self.store_g_to_i(&handle);
            }
        }

        // `cxx:71-83`: curve interferences of POINT geometry (no shape bind).
        let nb_curves = hds.ds().nb_curves();
        for i in 1..=nb_curves {
            let li = hds.ds().curve_interferences(i).to_vec();
            for handle in li {
                let is_point = match handle.read() {
                    Ok(itf) => itf.geometry_type() == Kind::Point,
                    Err(_) => false,
                };
                if !is_point {
                    continue;
                }
                self.store_g_to_i(&handle);
            }
        }
    }

    /// `Curve(I, C)` (`cxx:89-115`): the curve `I` sits on, looked up on the
    /// support first and then on every other interference of the same point.
    pub fn curve(&self, i: &InterferenceHandle, hds: &HDataStructure) -> Option<Curve> {
        let (is_point, g, sk, s) = match i.read() {
            Ok(itf) => (
                itf.geometry_type() == Kind::Point,
                itf.geometry(),
                itf.support_type(),
                itf.support(),
            ),
            Err(_) => return None,
        };
        if !is_point {
            return None;
        }
        // `cxx:97-102`
        if sk == Kind::Curve {
            return Some(hds.ds().curve(s).clone());
        }
        // `cxx:103-113`
        for other in self.same_interferences_at(g) {
            let Ok(itf) = other.read() else {
                continue;
            };
            if itf.support_type() == Kind::Curve {
                return Some(hds.ds().curve(itf.support()).clone());
            }
        }
        None
    }

    /// `Interferences(IP)` (`cxx:121-124`).
    pub fn interferences(&self, ip: i32) -> &[InterferenceHandle] {
        self.same_interferences_at(ip)
    }

    /// `myGToI(G)` read access; unbound answers the empty list, see the module
    /// note.
    fn same_interferences_at(&self, g: i32) -> &[InterferenceHandle] {
        match self.g_to_i.iter().find(|(k, _)| *k == g) {
            Some((_, list)) => list.as_slice(),
            None => &[],
        }
    }

    /// `SameInterferences(I)` (`cxx:128-131`).
    pub fn same_interferences(&self, i: &InterferenceHandle) -> &[InterferenceHandle] {
        match i.read() {
            Ok(itf) => self.same_interferences_at(itf.geometry()),
            Err(_) => &[],
        }
    }

    /// `ChangeSameInterferences(I)` (`cxx:135-140`).
    pub fn change_same_interferences(
        &mut self,
        i: &InterferenceHandle,
    ) -> Option<&mut InterferenceList> {
        let g = match i.read() {
            Ok(itf) => itf.geometry(),
            Err(_) => return None,
        };
        match self.g_to_i.iter_mut().find(|(k, _)| *k == g) {
            Some((_, list)) => Some(list),
            None => None,
        }
    }

    /// `EdgeSupport(I, E)` (`cxx:144-175`). `HDS` is unused in OCCT.
    pub fn edge_support(&self, i: &InterferenceHandle) -> Option<TopoShape> {
        let (is_point, g) = match i.read() {
            Ok(itf) => (itf.geometry_type() == Kind::Point, itf.geometry()),
            Err(_) => return None,
        };
        if !is_point {
            return None;
        }
        // `cxx:149-156`
        if let Some(s) = self.shape_of(i) {
            if s.shape_type() == crate::abs::ShapeType::Edge {
                return Some(s.clone());
            }
        }
        // `cxx:157-172`
        for other in self.same_interferences_at(g) {
            if let Some(s) = self.shape_of(other) {
                if s.shape_type() == crate::abs::ShapeType::Edge {
                    return Some(s.clone());
                }
            }
        }
        None
    }
    /// `FacesSupport(I, F1, F2)` (`cxx:179-191`).
    pub fn faces_support(
        &self,
        i: &InterferenceHandle,
        hds: &HDataStructure,
    ) -> Option<(TopoShape, TopoShape)> {
        let c = self.curve(i, hds)?;
        // `C.GetShapes(F1, F2)` (`cxx:187`) copies whatever the curve holds; a
        // null shape stays null in OCCT and every caller then compares it with
        // `IsSame`, so `None` short-circuits the same way.
        let (s1, s2) = c.get_shapes();
        Some((s1.cloned()?, s2.cloned()?))
    }

    /// `ParameterOnEdge(I, E, U)` (`cxx:193-227`).
    pub fn parameter_on_edge(
        &self,
        i: &InterferenceHandle,
        e: &TopoShape,
        _hds: &HDataStructure,
    ) -> Option<f64> {
        let (is_point, g) = match i.read() {
            Ok(itf) => (itf.geometry_type() == Kind::Point, itf.geometry()),
            Err(_) => return None,
        };
        if !is_point {
            return None;
        }
        // `cxx:200-208`
        if let Some(s) = self.shape_of(i) {
            if is_same(s, e) {
                return curve_point_parameter(i);
            }
        }
        // `cxx:209-224`
        for other in self.same_interferences_at(g) {
            if let Some(s) = self.shape_of(other) {
                if is_same(s, e) {
                    return curve_point_parameter(other);
                }
            }
        }
        None
    }

    /// `SetParameterOnEdge(I, E, U)` (`cxx:230-259`).
    pub fn set_parameter_on_edge(&self, i: &InterferenceHandle, e: &TopoShape, u: f64) {
        let (is_point, g) = match i.read() {
            Ok(itf) => (itf.geometry_type() == Kind::Point, itf.geometry()),
            Err(_) => return,
        };
        if !is_point {
            return;
        }
        // `cxx:235-243`
        if let Some(s) = self.shape_of(i) {
            if is_same(s, e) {
                set_curve_point_parameter(i, u);
            }
        }
        // `cxx:244-256`
        for other in self.same_interferences_at(g) {
            if let Some(s) = self.shape_of(other) {
                if is_same(s, e) {
                    set_curve_point_parameter(other, u);
                }
            }
        }
    }

    /// `SetPoint(I, IP)` (`cxx:263-275`): move every interference sharing the
    /// point of `I` onto `ip` and re-index them.
    pub fn set_point(&mut self, i: &InterferenceHandle, ip: i32) {
        let g = match i.read() {
            Ok(itf) => itf.geometry(),
            Err(_) => return,
        };
        if ip == g {
            return;
        }
        // `LI` is a reference to the list bound to the *old* geometry; the
        // driver appends to the new geometry's list while walking it in OCCT,
        // which cannot alias, so the list is snapshotted here.
        let li = self.same_interferences_at(g).to_vec();
        for ii in li {
            if let Ok(mut itf) = ii.write() {
                itf.set_geometry(ip);
            }
            self.store_g_to_i(&ii);
        }
    }

    /// `myInterToShape.Find(I)` (the `IsBound` / `operator()` pair of
    /// `cxx:149-151`, `:201-202`, `:236-237`).
    fn shape_of(&self, i: &InterferenceHandle) -> Option<&TopoShape> {
        self.inter_to_shape
            .iter()
            .find(|(k, _)| std::sync::Arc::ptr_eq(k, i))
            .map(|(_, s)| s)
    }
}

/// `occ::down_cast<TopOpeBRepDS_CurvePointInterference>(I)->Parameter()`.
fn curve_point_parameter(i: &InterferenceHandle) -> Option<f64> {
    let itf: std::sync::RwLockReadGuard<'_, Interference> = i.read().ok()?;
    if !itf.is_curve_point() {
        return None;
    }
    itf.parameter()
}

/// `occ::down_cast<TopOpeBRepDS_CurvePointInterference>(I)->Parameter(U)`.
fn set_curve_point_parameter(i: &InterferenceHandle, u: f64) {
    if let Ok(mut itf) = i.write() {
        if itf.is_curve_point() {
            itf.set_parameter(u);
        }
    }
}
