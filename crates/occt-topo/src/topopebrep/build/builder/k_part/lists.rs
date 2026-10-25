//! `TopOpeBRepBuild_Builder` — the `KP*` list and classification helpers.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_KPart.cxx`
//! (OCCT 8.0.0):
//! - `KPlhg` (`:1463-1494`, both overloads), `KPlhsd` (`:1504-1536`),
//!   `KPls` (`:2203-2231`), `KPContains` (`:2350-2362`),
//!   `KPSameDomain` (`:1720-1779`),
//! - `KPclasSS` (`:1549-1610`, all three overloads), `KPclassF`
//!   (`:2243-2271`), `KPclassFF` (`:2276-2297`), `KPiskoleFF` (`:2307-2344`).
//!
//! `KPlhg` / `KPlhsd` / `KPls` are `const` in OCCT and stay `&self` here; the
//! `KPclas*` family writes `myShapeClassifier`, so it needs `&mut self`.

use crate::abs::State;
use crate::shape::TopoShape;
use crate::topo_tools::is_same;
use crate::topexp::Explorer;

use super::Builder;

impl Builder {
    /// `KPlhg(S, T)` (`KPart.cxx:1463-1468`): how many `T` sub-shapes of `S`
    /// the DS gives a geometry to.
    pub fn kp_lhg(&self, s: &TopoShape, t: crate::abs::ShapeType) -> i32 {
        let mut l = Vec::new();
        self.kp_lhg_into(s, t, &mut l)
    }

    /// `KPlhg(S, T, L)` (`KPart.cxx:1475-1494`).
    pub fn kp_lhg_into(
        &self,
        s: &TopoShape,
        t: crate::abs::ShapeType,
        l: &mut Vec<TopoShape>,
    ) -> i32 {
        let mut n = 0;
        l.clear();
        let mut ex = Explorer::new(s, t);
        while ex.more() {
            let ss = ex.current().clone();
            if self
                .my_data_structure
                .as_ref()
                .map(|hds| hds.has_geometry(&ss))
                .unwrap_or(false)
            {
                n += 1;
                l.push(ss);
            }
            ex.next();
        }
        n
    }

    /// `KPlhsd(S, T)` (`KPart.cxx:1504-1509`): how many `T` sub-shapes of `S`
    /// are `HasSameDomain`.
    pub fn kp_lhsd(&self, s: &TopoShape, t: crate::abs::ShapeType) -> i32 {
        let mut l = Vec::new();
        self.kp_lhsd_into(s, t, &mut l)
    }

    /// `KPlhsd(S, T, L)` (`KPart.cxx:1517-1536`).
    pub fn kp_lhsd_into(
        &self,
        s: &TopoShape,
        t: crate::abs::ShapeType,
        l: &mut Vec<TopoShape>,
    ) -> i32 {
        let mut n = 0;
        l.clear();
        let mut ex = Explorer::new(s, t);
        while ex.more() {
            let ss = ex.current().clone();
            if self
                .my_data_structure
                .as_ref()
                .map(|hds| hds.has_same_domain(&ss, true))
                .unwrap_or(false)
            {
                n += 1;
                l.push(ss);
            }
            ex.next();
        }
        n
    }

    /// `KPls(S, T)` (`KPart.cxx:2203-2208`): the number of `T` sub-shapes.
    pub fn kp_ls(&self, s: &TopoShape, t: crate::abs::ShapeType) -> i32 {
        let mut l = Vec::new();
        self.kp_ls_into(s, t, &mut l)
    }

    /// `KPls(S, T, L)` (`KPart.cxx:2216-2231`).
    pub fn kp_ls_into(
        &self,
        s: &TopoShape,
        t: crate::abs::ShapeType,
        l: &mut Vec<TopoShape>,
    ) -> i32 {
        let mut n = 0;
        l.clear();
        let mut ex = Explorer::new(s, t);
        while ex.more() {
            n += 1;
            l.push(ex.current().clone());
            ex.next();
        }
        n
    }

    /// `KPContains(S, L)` (`KPart.cxx:2350-2362`): `S` is in `L` by `IsSame`.
    pub fn kp_contains(&self, s: &TopoShape, l: &[TopoShape]) -> bool {
        l.iter().any(|sl| is_same(sl, s))
    }

    /// `KPSameDomain(L1, L2)` (`KPart.cxx:1720-1779`).
    ///
    /// Each pass walks only the items prepended by the previous pass — the
    /// count OCCT keeps in `nl1` / `nl2` — and prepends, so the walk order is
    /// the new-items-first order of `NCollection_List::Prepend`.
    pub fn kp_same_domain(&self, l1: &mut Vec<TopoShape>, l2: &mut Vec<TopoShape>) {
        let mut nl1 = l1.len();
        let mut nl2 = l2.len();

        while nl1 > 0 || nl2 > 0 {
            // cxx:1729-1752
            for i in 0..nl1 {
                let s1 = l1[i].clone();
                let sd = self
                    .my_data_structure
                    .as_ref()
                    .map(|hds| hds.same_domain(&s1).to_vec())
                    .unwrap_or_default();
                for s2 in sd {
                    if !self.kp_contains(&s2, l2) {
                        l2.insert(0, s2);
                        nl2 += 1;
                    }
                }
            }
            nl1 = 0;

            // cxx:1754-1778
            for i in 0..nl2 {
                let s2 = l2[i].clone();
                let sd = self
                    .my_data_structure
                    .as_ref()
                    .map(|hds| hds.same_domain(&s2).to_vec())
                    .unwrap_or_default();
                for s1 in sd {
                    if !self.kp_contains(&s1, l1) {
                        l1.insert(0, s1);
                        nl1 += 1;
                    }
                }
            }
            nl2 = 0;
        }
    }

    /// `KPclasSS(S1, exLS1, S2)` (`KPart.cxx:1549-1567`).
    pub fn kp_clas_ss_ls(&mut self, s1: &TopoShape, ex_ls1: &[TopoShape], s2: &TopoShape) -> State {
        let st = self
            .my_shape_classifier
            .borrow_mut()
            .state_shape_shape_av_ls(s1, ex_ls1, s2);
        kp_state(st)
    }

    /// `KPclasSS(S1, exS1, S2)` (`KPart.cxx:1578-1596`).
    pub fn kp_clas_ss_av(&mut self, s1: &TopoShape, ex_s1: Option<&TopoShape>, s2: &TopoShape) -> State {
        let st = match ex_s1 {
            Some(ex_s1) => {
                self.my_shape_classifier
                    .borrow_mut()
                    .state_shape_shape_av_s(s1, ex_s1, s2)
            }
            None => self
                .my_shape_classifier
                .borrow_mut()
                .state_shape_shape_null_avoid(s1, s2),
        };
        kp_state(st)
    }

    /// `KPclasSS(S1, S2)` (`KPart.cxx:1605-1610`): no avoidance.
    pub fn kp_clas_ss(&mut self, s1: &TopoShape, s2: &TopoShape) -> State {
        self.kp_clas_ss_av(s1, None, s2)
    }

    /// `KPclassF(F1, F2)` (`KPart.cxx:2243-2271`): `F1`'s first edge
    /// classified against `F2` in same-domain mode.
    pub fn kp_class_f(&mut self, f1: &TopoShape, f2: &TopoShape) -> State {
        if f1.is_null() || f2.is_null() {
            return State::Unknown;
        }
        let f1f = f1.oriented(crate::abs::Orientation::Forward);
        let f2f = f2.oriented(crate::abs::Orientation::Forward);

        let mut le1 = Vec::new();
        let ne1 = self.kp_ls_into(&f1f, crate::abs::ShapeType::Edge, &mut le1);
        if ne1 == 0 {
            return State::Unknown;
        }
        let e1 = le1[0].clone();

        // cxx:2267-2269: `StateShapeShape(e1, F2F, isamdom = 1)`.
        let st = self
            .my_shape_classifier
            .borrow_mut()
            .state_shape_shape(&e1, &f2f, 1);
        kp_state(st)
    }

    /// `KPclassFF(F1, F2, St1, St2)` (`KPart.cxx:2276-2297`).
    pub fn kp_class_ff(&mut self, f1: &TopoShape, f2: &TopoShape) -> (State, State) {
        let st1 = self.kp_class_f(f1, f2);
        let st2 = self.kp_class_f(f2, f1);
        (st1, st2)
    }

    /// `KPiskoleFF(F1, F2, St1, St2)` (`KPart.cxx:2307-2344`): true when both
    /// states are `IN`/`OUT` and differ, the `iskole` face configuration.
    pub fn kp_iskole_ff(&mut self, f1: &TopoShape, f2: &TopoShape) -> bool {
        matches!(self.kp_iskole_ff_states(f1, f2), Some((true, _, _)))
    }
    /// `KPiskoleFF(F1, F2, St1, St2)` (`KPart.cxx:2307-2344`): `classok` plus
    /// the two states. `None` stands for the `Standard_ConstructionError`
    /// `KPclassFF` raises on a null shape, which `KPiskole` never meets.
    pub fn kp_iskole_ff_states(
        &mut self,
        f1: &TopoShape,
        f2: &TopoShape,
    ) -> Option<(bool, State, State)> {
        if f1.is_null() || f2.is_null() {
            return None;
        }
        // cxx:2329
        let (st1, st2) = self.kp_class_ff(f1, f2);
        // cxx:2330-2343
        let st1ok = st1 == State::Out || st1 == State::In;
        let st2ok = st2 == State::Out || st2 == State::In;
        if !st1ok || !st2ok {
            return Some((false, st1, st2));
        }
        Some((st1 != st2, st1, st2))
    }

    /// `myKPMAPf1f2.Bind(f, {})` when `f` is unbound — the
    /// `if (!myKPMAPf1f2.IsBound(f)) myKPMAPf1f2.Bind(f, los)` idiom of
    /// `KPiskole` (`kpkole.cxx:665-668`, `:674-678`, `:683-686`, `:691-695`).
    pub fn kp_map_f1f2_bind(&mut self, f: &TopoShape) {
        if self.kp_map_f1f2_contains(f) {
            return;
        }
        self.my_kp_map_f1f2.push((f.clone(), Vec::new()));
    }

    /// `myKPMAPf1f2.ChangeFind(f).Append(g)` (`kpkole.cxx:670-671`,
    /// `:688-689`): `f` must be bound by [`Builder::kp_map_f1f2_bind`].
    pub fn kp_map_f1f2_append(&mut self, f: &TopoShape, g: &TopoShape) {
        if let Some((_, los)) = self
            .my_kp_map_f1f2
            .iter_mut()
            .find(|(key, _)| is_same(key, f))
        {
            los.push(g.clone());
        }
    }

    /// `myKPMAPf1f2.IsBound(f)` (`kpkole.cxx:427`, `:448`, `:664`, ...).
    pub fn kp_map_f1f2_contains(&self, f: &TopoShape) -> bool {
        self.my_kp_map_f1f2
            .iter()
            .any(|(key, _)| is_same(key, f))
    }

    /// `myKPMAPf1f2` in `NCollection_DataMap` iteration order. The port keeps
    /// the insertion order in a `Vec` because `MergeKPartiskole` reads the
    /// first non-empty entry as `pfOU` (`kpkole.cxx:135-152`).
    pub fn kp_map_f1f2(&self) -> &[(TopoShape, Vec<TopoShape>)] {
        &self.my_kp_map_f1f2
    }

    /// The `KPiskole` tail (`kpkole.cxx:722-729`): `DS.InitSectionEdges()`
    /// then `DS.AddSectionEdge(E)` for every collected edge.
    pub fn kp_set_section_edges(&mut self, edges: &[TopoShape]) {
        let Some(hds) = self.my_data_structure.as_mut() else {
            return;
        };
        let ds = hds.change_ds();
        ds.init_section_edges();
        for e in edges {
            ds.add_section_edge(e);
        }
    }
}

/// `StateShapeShape` throws `Standard_ConstructionError` on a null reference;
/// the port's classifier reports `Err`, which `Perform` cannot reach with the
/// states `KPclasSS` passes, and the previous state (`UNKNOWN`) stands.
fn kp_state(st: Result<State, String>) -> State {
    st.unwrap_or(State::Unknown)
}
