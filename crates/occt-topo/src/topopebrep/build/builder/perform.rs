//! `TopOpeBRepBuild_Builder` — `Perform` and the vertex build pass.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Builder.cxx:116-145`
//! (both `Perform` overloads), `TopOpeBRepBuild_BuildVertices.cxx:25-34`
//! (`BuildVertices`), `TopOpeBRepBuild_KPart.cxx:169-172` (`IsKPart`)
//! (OCCT 8.0.0).
//!
//! `Perform` is the entry point of the class: it clears the split / merged
//! maps, keeps the DS, builds one vertex per DS point, then builds the edges of
//! every DS curve and the faces of every DS surface, splits the section edges
//! and finally filters and reduces the face interferences of `mySplitON`.
//!
//! UNPORTED, with their OCCT location — the step of `perform_with_shapes` that
//! cannot run yet:
//! - `FindIsKPart` (`TopOpeBRepBuild_KPart.cxx:111-166`) and `MergeKPart`
//!   (`:176-220`): need the five `KP*` predicates of `KPart.cxx`.
//!
//! Every step of `Perform` itself is ported: `BuildEdges` runs at `cxx:126`
//! (`build/builder/build_edges.rs`) and fills `myNewEdges`, which
//! `BuildFaces` at `cxx:127` consumes.
//!
//! `Filter::ProcessFaceInterferences` and `Reducer::ProcessFaceInterferences`
//! *are* ported (`ds/filter/mod.rs`, `ds/reducer.rs`), so the last two
//! statements of `Perform` run here verbatim.

use std::rc::Rc;

use occt_core::gp::GpPnt;

use crate::topopebrep::ds::filter::Filter;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::reducer::Reducer;

use super::Builder;

impl Builder {
    /// `Perform(HDS)` (`Builder.cxx:116-133`).
    pub fn perform(&mut self, hds: HDataStructure) {
        self.clear();
        self.my_data_structure = Some(hds);

        // cxx:124
        self.build_vertices();

        // cxx:125 `SplitEvisoONperiodicF()`
        self.split_eviso_on_periodic_f();

        // cxx:126 `BuildEdges(HDS)`
        self.build_edges_all();

        // cxx:127
        self.build_faces_all();

        // cxx:128-129
        self.my_is_k_part = 0;
        self.init_section();

        // cxx:130
        self.split_section_edges();

        // cxx:131-133
        let psc = Some(Rc::clone(&self.my_shape_classifier));
        if let Some(hds) = self.my_data_structure.as_mut() {
            let mut f = Filter::new(hds, psc);
            f.process_face_interferences(&self.my_split_on);
        }
        if let Some(hds) = self.my_data_structure.as_mut() {
            let mut r = Reducer::new(hds);
            let _ = r.process_face_interferences(&self.my_split_on);
        }
    }

    /// `Perform(HDS, S1, S2)` (`Builder.cxx:138-145`).
    pub fn perform_with_shapes(
        &mut self,
        hds: HDataStructure,
        s1: Option<&crate::shape::TopoShape>,
        s2: Option<&crate::shape::TopoShape>,
    ) {
        self.perform(hds);
        self.my_shape1 = s1.cloned();
        self.my_shape2 = s2.cloned();
        // cxx:144
        self.my_is_k_part = self.find_is_k_part();
    }

    /// `BuildVertices(HDS)` (`BuildVertices.cxx:25-34`).
    ///
    /// `myNewVertices = new HArray1<TopoDS_Shape>(0, n)` becomes a `n`-long
    /// 1-based array: `ChangeNewVertex(iP)` writes index `iP` for
    /// `iP = 1..n`, and index `0` of the OCCT array is never read.
    pub fn build_vertices(&mut self) {
        let n = self
            .my_data_structure
            .as_ref()
            .map(|h| h.nb_points())
            .unwrap_or(0);

        // `HDS->Point(iP)` returns a const reference and throws when the point
        // is absent; the points are copied out first so `myNewVertices` can be
        // written while the DS is only read.
        let mut points: Vec<Option<(GpPnt, f64)>> = Vec::with_capacity(n as usize);
        if let Some(hds) = self.my_data_structure.as_ref() {
            for ip in 1..=n {
                points.push(hds.point(ip).map(|p| (*p.point(), p.tolerance())));
            }
        }

        self.set_new_vertices_len(n as usize);

        for (k, item) in points.into_iter().enumerate() {
            let ip = (k + 1) as i32;
            let Some((p, tol)) = item else {
                continue;
            };
            let v = self.my_build_tool.make_vertex(p, tol);
            if let Some(slot) = self.change_new_vertex(ip) {
                *slot = Some(v);
            }
        }
    }

    /// `IsKPart()` (`KPart.cxx:169-172`).
    pub fn is_k_part(&self) -> i32 {
        self.my_is_k_part
    }
}
