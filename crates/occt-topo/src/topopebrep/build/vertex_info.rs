//! `TopOpeBRepBuild_VertexInfo` — the edges met around one vertex while the
//! `Builder1` pass walks it.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_VertexInfo.hxx:31-73`,
//! `.cxx:16-153` (OCCT 8.0.0).
//!
//! `myEdgesIn` / `myEdgesOut` / `myLocalEdgesOut` are
//! `NCollection_IndexedMap<TopoDS_Shape>` ([`IndexedShapeSet`], index order =
//! insertion order) and `myEdgesPassed` is a `NCollection_List<TopoDS_Shape>`
//! prepended to ([`Vec`] in reverse use order).

use crate::shape::TopoShape;
use crate::topopebrep::ds::data_structure::IndexedShapeSet;

/// `class TopOpeBRepBuild_VertexInfo` (`hxx:38`).
#[derive(Default)]
pub struct VertexInfo {
    /// `myVertex` (`hxx:65`).
    my_vertex: Option<TopoShape>,
    /// `myCurrent` (`hxx:66`).
    my_current: Option<TopoShape>,
    /// `myCurrentIn` (`hxx:67`).
    my_current_in: Option<TopoShape>,
    /// `mySmart` (`hxx:68`).
    my_smart: bool,
    /// `myEdgesIn` (`hxx:69`).
    my_edges_in: IndexedShapeSet,
    /// `myEdgesOut` (`hxx:70`).
    my_edges_out: IndexedShapeSet,
    /// `myLocalEdgesOut` (`hxx:71`).
    my_local_edges_out: IndexedShapeSet,
    /// `myEdgesPassed` (`hxx:72`).
    my_edges_passed: Vec<TopoShape>,
    /// `myFoundOut` (`hxx:73`).
    my_found_out: i32,
}

impl VertexInfo {
    /// `TopOpeBRepBuild_VertexInfo()` (`cxx:16-21`).
    pub fn new() -> Self {
        Self {
            my_vertex: None,
            my_current: None,
            my_current_in: None,
            my_smart: false,
            my_edges_in: IndexedShapeSet::new(),
            my_edges_out: IndexedShapeSet::new(),
            my_local_edges_out: IndexedShapeSet::new(),
            my_edges_passed: Vec::new(),
            my_found_out: 0,
        }
    }

    /// `SetVertex(aV)` (`cxx:25-28`).
    pub fn set_vertex(&mut self, v: &TopoShape) {
        self.my_vertex = Some(v.clone());
    }

    /// `Vertex()` (`cxx:32-35`).
    pub fn vertex(&self) -> Option<&TopoShape> {
        self.my_vertex.as_ref()
    }

    /// `SetSmart(aFlag)` (`cxx:39-42`).
    pub fn set_smart(&mut self, flag: bool) {
        self.my_smart = flag;
    }

    /// `Smart()` (`cxx:46-49`).
    pub fn smart(&self) -> bool {
        self.my_smart
    }

    /// `NbCases()` (`cxx:53-56`): the extent of `myLocalEdgesOut`.
    pub fn nb_cases(&self) -> i32 {
        self.my_local_edges_out.extent()
    }

    /// `FoundOut()` (`cxx:60-63`).
    pub fn found_out(&self) -> i32 {
        self.my_found_out
    }

    /// `AddIn(anE)` (`cxx:67-70`).
    pub fn add_in(&mut self, e: &TopoShape) {
        self.my_edges_in.add(e);
    }

    /// `AddOut(anE)` (`cxx:74-77`).
    pub fn add_out(&mut self, e: &TopoShape) {
        self.my_edges_out.add(e);
    }

    /// `SetCurrentIn(anE)` (`cxx:81-84`). OCCT stores a `TopoDS_Edge` that may
    /// be null (`Tools2d::Path` starts the walk with `aNullEdge`, `cxx:57-58`).
    pub fn set_current_in(&mut self, e: Option<&TopoShape>) {
        self.my_current_in = e.cloned();
    }

    /// `EdgesIn()` (`cxx:88-91`).
    pub fn edges_in(&self) -> &IndexedShapeSet {
        &self.my_edges_in
    }

    /// `EdgesOut()` (`cxx:95-98`).
    pub fn edges_out(&self) -> &IndexedShapeSet {
        &self.my_edges_out
    }

    /// `ChangeEdgesOut()` (`cxx:102-105`).
    pub fn change_edges_out(&mut self) -> &mut IndexedShapeSet {
        &mut self.my_edges_out
    }

    /// `Dump()` (`cxx:109-117`): `printf` has no port and the method is debug
    /// output only.
    pub fn dump(&self) {}

    /// `AppendPassed(anE)` (`cxx:121-124`): `myEdgesPassed.Prepend(anE)`.
    pub fn append_passed(&mut self, e: &TopoShape) {
        self.my_edges_passed.insert(0, e.clone());
    }

    /// `RemovePassed()` (`cxx:128-131`): `myEdgesPassed.RemoveFirst()`.
    pub fn remove_passed(&mut self) {
        if !self.my_edges_passed.is_empty() {
            self.my_edges_passed.remove(0);
        }
    }

    /// `ListPassed()` (`cxx:135-138`).
    pub fn list_passed(&self) -> &[TopoShape] {
        &self.my_edges_passed
    }

    /// `Prepare(aL)` (`cxx:142-163`).
    pub fn prepare(&mut self, a_l: &[TopoShape]) {
        self.my_local_edges_out.clear();

        let mut tmp_map = IndexedShapeSet::new();
        for e in a_l {
            tmp_map.add(e);
        }

        let nb = self.my_edges_out.extent();
        for i in 1..=nb {
            let Some(a_e) = self.my_edges_out.find_key(i).cloned() else {
                continue;
            };
            if !tmp_map.contains(&a_e) {
                self.my_local_edges_out.add(&a_e);
            }
        }

        tmp_map.clear();
    }

    /// `CurrentOut()` (`cxx:147-173`).
    ///
    /// Returns the first edge of `myLocalEdgesOut` not in `myEdgesPassed` and
    /// sets `myFoundOut = 1`; when every one was passed it clears `myCurrent`,
    /// sets `myFoundOut = 0` and answers null.
    pub fn current_out(&mut self) -> Option<TopoShape> {
        let a_nb_out = self.my_local_edges_out.extent();

        let mut a_map_passed = IndexedShapeSet::new();
        for e in &self.my_edges_passed {
            a_map_passed.add(e);
        }

        for i in 1..=a_nb_out {
            let Some(e) = self.my_local_edges_out.find_key(i).cloned() else {
                continue;
            };
            if !a_map_passed.contains(&e) {
                self.my_current = Some(e);
                self.my_found_out = 1;
                return self.my_current.clone();
            }
        }
        self.my_found_out = 0;
        self.my_current = None;
        None
    }
}
