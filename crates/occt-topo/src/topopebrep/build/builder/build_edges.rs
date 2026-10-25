//! `TopOpeBRepBuild_Builder` — `BuildEdges`, both overloads.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_BuildEdges.cxx:38-218`
//! (OCCT 8.0.0).
//!
//! `BuildEdges(iC, HDS)` (`cxx:38-100`) rebuilds the edge of DS curve `iC`
//! from its point interferences: a `PaveSet` on a fresh `MakeEdge(anEdge,
//! curC, BDS)` edge, a `PaveClassifier`, `MakeEdges` to split the edge at the
//! pave vertices, then `RecomputeCurves` on every new edge. When
//! `RecomputeCurves` created a new DS curve, the old curve `iC` is removed and
//! the new edges are registered under the new curve index; otherwise the new
//! edges only inherit the geometry of `anEdge` (`UpdateEdge`).
//!
//! `BuildEdges(HDS)` (`cxx:104-218`) drives the per-curve overload over every
//! curve with `Mother() == 0`, after dropping the curves a previous build left
//! behind (`cxx:110-132`), and finally removes the DS points that no remaining
//! interference refers to (`cxx:144-217`).
//!
//! The two passes that walk the DS twice (`cxx:110-142`) take a snapshot of
//! the curve indices first: OCCT holds a `CurveExplorer` over `BDS` across
//! `BDS.RemoveCurve` / `MakeCurve` calls, which the port's borrow rules forbid,
//! and the explorer's `myMax` is fixed at `Init` in both languages, so a
//! snapshot visits exactly the same curves.
//!
//! `BDS.Shape(is)`, `BDS.ShapeInterferences(is)` and the point-counting
//! `cex.Init(BDS)` all take their default `FindKeep` of `true`
//! (`TopOpeBRepDS_DataStructure.hxx:254`, `:149`,
//! `TopOpeBRepDS_CurveExplorer.hxx:37`); the two passes before them pass
//! `false` explicitly (`cxx:111`, `:134`).

use crate::abs::{ShapeType, State};
use crate::shape::TopoShape;
use crate::topopebrep::build::{EdgeBuilder, LoopSetSource, PaveClassifier, PaveSet};
use crate::topopebrep::ds::curve_explorer::CurveExplorer;
use crate::topopebrep::ds::kind::Kind;

use super::Builder;

impl Builder {
    /// `BuildEdges(iC, HDS)` (`BuildEdges.cxx:38-100`).
    pub fn build_edges(&mut self, i_c: i32) {
        // cxx:47-52: `C3D`, `I1`, `I2` all null => nothing to build.
        let cur_c = {
            let Some(hds) = self.my_data_structure.as_ref() else {
                return;
            };
            let c = hds.curve(i_c);
            if c.curve().is_none() && c.get_sci1().is_none() && c.get_sci2().is_none() {
                return;
            }
            c.clone()
        };

        // cxx:55-57: `MakeEdge(anEdge, curC, HDS->DS())`.
        let mut an_edge = TopoShape::new(ShapeType::Edge);
        {
            let Some(hds) = self.my_data_structure.as_ref() else {
                return;
            };
            self.my_build_tool
                .make_edge_with_ds(&mut an_edge, &cur_c, hds.ds());
        }

        // cxx:58-61: the pave set of the edge and the DS points on the curve.
        let mut pvs = PaveSet::new(&an_edge);
        let mut epit = match self.my_data_structure.as_ref() {
            Some(hds) => hds.curve_points(i_c),
            None => return,
        };
        self.fill_vertex_set(&mut epit, State::In, &mut pvs);

        // cxx:62-68
        let mut vcl = PaveClassifier::new(&an_edge);
        if pvs.has_equal_parameters() {
            let p = pvs.equal_parameters();
            vcl.set_first_parameter(p);
        }
        let closvert = pvs.closed_vertices();
        vcl.set_closed_vertices(closvert);

        // cxx:69-73
        pvs.init_loop();
        if !pvs.more_loop() {
            return;
        }

        // cxx:74-78: `TopOpeBRepBuild_EdgeBuilder EDBU(PVS, VCL)` then
        // `MakeEdges(anEdge, EDBU, EL)` into `ChangeNewEdges(iC)`; the list is
        // taken out of the map because `MakeEdges` needs `&mut self` and the
        // loop below appends to another entry of the same map.
        let mut edbu = EdgeBuilder::with_loop_set(&mut pvs, &mut vcl, false);
        let mut el = std::mem::take(self.change_new_edges(i_c));
        self.make_edges(&an_edge, &mut edbu, &mut el);

        // cxx:79-87
        let mut inew_c = -1;
        for new_edge in &el {
            let bt = self.my_build_tool;
            let Some(hds) = self.my_data_structure.as_mut() else {
                break;
            };
            bt.recompute_curves(&cur_c, &an_edge, &mut inew_c, hds);
            if inew_c != -1 {
                self.change_new_edges(inew_c).push(new_edge.clone());
            }
        }

        // cxx:88-99
        if inew_c != -1 {
            if let Some(hds) = self.my_data_structure.as_mut() {
                hds.remove_curve(i_c);
            }
        } else {
            for new_edge in &el {
                self.my_build_tool.update_edge(&an_edge, new_edge);
            }
        }

        // `EL` is the live `ChangeNewEdges(iC)` list in OCCT; put it back.
        *self.change_new_edges(i_c) = el;
    }

    /// `BuildEdges(HDS)` (`BuildEdges.cxx:104-218`).
    pub fn build_edges_all(&mut self) {
        if self.my_data_structure.is_none() {
            return;
        }

        // cxx:108
        self.my_new_edges.clear();

        // cxx:110-123: the first kept curve that a `MakeCurve` created.
        let ick = {
            let hds = self
                .my_data_structure
                .as_ref()
                .expect("my_data_structure checked above");
            let mut cex = CurveExplorer::with_ds(hds.ds(), false);
            let mut ick = 0;
            while cex.more() {
                let ic = cex.index();
                let ck = cex.is_curve_keep(ic);
                let im = cex.curve_at(ic).map(|c| c.mother()).unwrap_or(0);
                if ck && im != 0 {
                    ick = ic;
                    break;
                }
                cex.next();
            }
            ick
        };

        // cxx:124-132: drop every kept curve, then forget every curve from
        // `ick` on. `RemoveCurve` only clears `Keep` in the port, so the
        // truncation is what actually hides the daughters, exactly as
        // `ChangeNbCurves(ick - 1)` does in OCCT.
        if ick != 0 {
            let to_remove: Vec<i32> = {
                let hds = self
                    .my_data_structure
                    .as_ref()
                    .expect("my_data_structure checked above");
                let mut cex = CurveExplorer::with_ds(hds.ds(), true);
                let mut v = Vec::new();
                while cex.more() {
                    v.push(cex.index());
                    cex.next();
                }
                v
            };
            if let Some(hds) = self.my_data_structure.as_mut() {
                for ic in to_remove {
                    hds.remove_curve(ic);
                }
                hds.change_ds().change_nb_curves(ick - 1);
            }
        }

        // cxx:134-142: build the edges of every curve that is not a daughter.
        let curve_indices: Vec<i32> = {
            let hds = self
                .my_data_structure
                .as_ref()
                .expect("my_data_structure checked above");
            let mut cex = CurveExplorer::with_ds(hds.ds(), false);
            let mut v = Vec::new();
            while cex.more() {
                v.push(cex.index());
                cex.next();
            }
            v
        };
        for ic in curve_indices {
            let im = self
                .my_data_structure
                .as_ref()
                .map(|h| h.curve(ic).mother())
                .unwrap_or(0);
            if im != 0 {
                continue;
            }
            self.build_edges(ic);
        }

        // cxx:144-217: count the references to every DS point from the
        // interferences of the kept curves and of the DS edges, then drop the
        // points nothing refers to.
        let np = {
            let hds = self
                .my_data_structure
                .as_ref()
                .expect("my_data_structure checked above");
            hds.nb_points()
        };
        // `NCollection_HArray1<int> tp(0, np, 0)`: index 0 is never read.
        let mut tp = vec![0i32; (np + 1).max(1) as usize];

        {
            // cxx:147: `cex.Init(BDS)` => `FindOnlyKeep = true`.
            let hds = self
                .my_data_structure
                .as_ref()
                .expect("my_data_structure checked above");
            let ds = hds.ds();
            let mut cex = CurveExplorer::with_ds(ds, true);
            while cex.more() {
                let ic = cex.index();
                for i in ds.curve_interferences(ic) {
                    let guard = i.read().expect("interference lock");
                    count_point_reference(&mut tp, np, guard.geometry_type(), guard.geometry());
                    count_point_reference(&mut tp, np, guard.support_type(), guard.support());
                }
                cex.next();
            }
        }

        {
            // cxx:175-211
            let hds = self
                .my_data_structure
                .as_ref()
                .expect("my_data_structure checked above");
            let ds = hds.ds();
            let ns = ds.nb_shapes();
            for is in 1..=ns {
                let Some(s) = ds.shape(is, true) else {
                    continue;
                };
                if s.shape_type() != ShapeType::Edge {
                    continue;
                }
                for i in ds.shape_interferences_at(is, true) {
                    let guard = i.read().expect("interference lock");
                    count_point_reference(&mut tp, np, guard.geometry_type(), guard.geometry());
                    count_point_reference(&mut tp, np, guard.support_type(), guard.support());
                }
            }
        }

        // cxx:212-217
        if let Some(hds) = self.my_data_structure.as_mut() {
            let ds = hds.change_ds();
            for ip in 1..=np {
                if tp[ip as usize] == 0 {
                    ds.remove_point(ip);
                }
            }
        }
    }
}

/// One `tp.ChangeValue(I) = tp.Value(I) + 1` of `BuildEdges.cxx:157-172` /
/// `:193-208`.
///
/// The `POINT` geometry arm is bounded by `ig <= np` in OCCT (`cxx:160`);
/// the two `POINT` support arms are not (`cxx:169-171`, `:202-204`), and an
/// out-of-range index would raise `Standard_RangeError` there. The port skips
/// such an index instead of panicking.
fn count_point_reference(tp: &mut [i32], np: i32, kind: Kind, index: i32) {
    if kind != Kind::Point {
        return;
    }
    if index < 0 || index > np {
        return;
    }
    tp[index as usize] += 1;
}
