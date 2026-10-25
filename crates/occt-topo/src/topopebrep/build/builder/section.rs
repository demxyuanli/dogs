//! `TopOpeBRepBuild_Builder` — the section edge lists.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Section.cxx:1234-1392`
//! (`Section`, `SectionCurves`, `SectionEdges`, `FillSecEdgeAncestorMap`)
//! (OCCT 8.0.0).
//!
//! `Section()` memoises `mySection`: the new edges of every DS curve, then the
//! section edges of the DS (their `Split(ON)` / `Split(IN)` lists, or the edge
//! itself when it carries neither geometry nor a same-domain partner).
//!
//! `SectionEdges` runs [`Builder::split_section_edges`]
//! (`Section.cxx:304-915`) first, exactly as `cxx:1294` does; a builder that
//! already ran `Perform` (`Builder.cxx:130`) or a previous `SectionEdges` call
//! sees `mySplitSectionEdgesDone` and skips it.

use std::collections::{HashMap, HashSet};

use crate::abs::State;
use crate::shape::TopoShape;
use crate::topo_tools::{shape_map_key, ShapeMapKey};
use crate::topopebrep::ds::curve_explorer::CurveExplorer;

use super::Builder;

impl Builder {
    /// `Section()` (`Section.cxx:1234-1244`).
    ///
    /// The OCCT signature returns `const List&`; the port returns the slice
    /// because `mySection` lives in the builder.
    pub fn section(&mut self) -> &[TopoShape] {
        if !self.my_section_done {
            self.my_section_done = true;
            self.section_curves();
            self.section_edges();
        }
        &self.my_section
    }

    /// `Section(L)` (`Section.cxx:1248-1251`): `L = Section()`.
    pub fn section_list(&mut self, l: &mut Vec<TopoShape>) {
        *l = self.section().to_vec();
    }

    /// `SectionCurves(LSE)` (`Section.cxx:1255-1268`): every new edge of every
    /// DS curve, in `CurveExplorer` order, appended to `mySection`.
    pub fn section_curves(&mut self) {
        let mut collected: Vec<TopoShape> = Vec::new();
        self.section_curves_into(&mut collected);
        self.my_section.extend(collected);
    }

    /// `SectionCurves(LSE)` (`Section.cxx:1255-1268`) into a caller-owned list
    /// — the signature `TopOpeBRepBuild_HBuilder::InitSection(1)` calls
    /// (`HBuilder.cxx:241`).
    pub fn section_curves_into(&mut self, lse: &mut Vec<TopoShape>) {
        let mut collected: Vec<TopoShape> = Vec::new();
        if let Some(hds) = self.my_data_structure.as_ref() {
            let mut cex = CurveExplorer::with_ds(hds.ds(), false);
            while cex.more() {
                let ic = cex.index();
                collected.extend(self.new_edges(ic).iter().cloned());
                cex.next();
            }
        }
        lse.extend(collected);
    }

    /// `SectionEdges(LSE)` (`Section.cxx:1271-1351`), appended to `mySection`.
    pub fn section_edges(&mut self) {
        let mut collected: Vec<TopoShape> = Vec::new();
        self.section_edges_into(&mut collected);
        self.my_section.extend(collected);
    }

    /// `SectionEdges(LSE)` (`Section.cxx:1271-1351`) into a caller-owned list —
    /// the signature `TopOpeBRepBuild_HBuilder::InitSection(2)` calls
    /// (`HBuilder.cxx:245`).
    pub fn section_edges_into(&mut self, lse: &mut Vec<TopoShape>) {
        let mut out: Vec<TopoShape> = Vec::new();
        let mut mos: HashSet<ShapeMapKey> = HashSet::new();

        let nes = match self.my_data_structure.as_ref() {
            Some(hds) => hds.ds().nb_section_edges(),
            None => 0,
        };

        // `cxx:1275-1292`: the iskole case keeps the section edges as they are.
        if self.is_k_part() == 1 {
            for i in 1..=nes {
                let Some(es) = self.section_edge(i) else {
                    continue;
                };
                out.push(es);
            }
            lse.extend(out);
            return;
        }

        // `cxx:1294` `SplitSectionEdges()`. It is a no-op once `Perform`
        // (`Builder.cxx:130`) or a second `SectionEdges` call ran it.
        self.split_section_edges();

        for i in 1..=nes {
            let Some(es) = self.section_edge(i) else {
                continue;
            };
            // cxx:1308-1310
            let issplit_in = self.is_split(&es, State::In);
            let issplit_on = self.is_split(&es, State::On);
            let mut staspl = State::Unknown;
            if issplit_on {
                staspl = State::On;
            } else if issplit_in {
                staspl = State::In;
            }

            // cxx:1325-1334
            if staspl != State::Unknown {
                for s in self.splits(&es, staspl) {
                    if mos.insert(shape_map_key(s)) {
                        out.push(s.clone());
                    }
                }
            } else {
                // cxx:1337-1348
                let (hasgeom, hassame) = match self.my_data_structure.as_ref() {
                    Some(hds) => (hds.has_geometry(&es), hds.has_same_domain(&es, true)),
                    None => (false, false),
                };
                let take = !(hasgeom || hassame);
                if take && mos.insert(shape_map_key(&es)) {
                    out.push(es);
                }
            }
        }

        lse.extend(out);
    }

    /// `FillSecEdgeAncestorMap(aShapeRank, aMapON, anAncMap)` (`Section.cxx:1355-1392`).
    ///
    /// `anAncMap` is a `DataMap<TopoDS_Shape, TopoDS_Shape>`; `Bind` replaces
    /// the value of an already bound key, which a `HashMap` reproduces. The
    /// OCCT traversal order of the map is unspecified, so the ordering of the
    /// port's map is not part of the contract.
    pub fn fill_sec_edge_ancestor_map(
        &self,
        a_shape_rank: i32,
        a_map_on: &HashSet<ShapeMapKey>,
        an_anc_map: &mut HashMap<ShapeMapKey, TopoShape>,
    ) {
        let nse = match self.my_data_structure.as_ref() {
            Some(hds) => hds.ds().nb_section_edges(),
            None => 0,
        };
        for i in 1..=nse {
            let Some(es) = self.section_edge(i) else {
                continue;
            };
            // cxx:1367
            if self.shape_rank(&es) != a_shape_rank {
                continue;
            }
            // cxx:1370-1374
            if a_map_on.contains(&shape_map_key(&es)) {
                an_anc_map.insert(shape_map_key(&es), es);
                continue;
            }
            // cxx:1376-1391
            for st in [State::In, State::On, State::Out] {
                if !self.is_split(&es, st) {
                    continue;
                }
                for a_s in self.splits(&es, st) {
                    if a_map_on.contains(&shape_map_key(a_s)) {
                        an_anc_map.insert(shape_map_key(a_s), es.clone());
                    }
                }
            }
        }
    }

    /// `BDS.SectionEdge(i)` (`TopOpeBRepDS_DataStructure.hxx:234-238`), with the
    /// null test every caller of the loop performs (`cxx:1284`, `:1304`).
    fn section_edge(&self, i: i32) -> Option<TopoShape> {
        self.my_data_structure
            .as_ref()
            .and_then(|hds| hds.ds().section_edge(i, true))
    }
}
