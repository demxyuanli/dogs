//! `TopOpeBRepBuild_Builder::GFillCurveTopologyWES` — the edges carved out by
//! a face's section curves become start elements of a `WireEdgeSet`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_fctwes.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `GFillCurveTopologyWES(F1, G1, WES)` | `60-136` |
//! | `GFillCurveTopologyWES(FCit, G1, WES)` | `138-268` |
//!
//! The first overload walks the face's section curves, skips the ones whose
//! interference is not on the shape the `GTopo` asks for and the ones whose
//! support is a same-domain partner of `F1` (`cxx:92-130`), then hands each
//! surviving curve to the second overload.
//!
//! The second overload turns the transition of the curve into an edge
//! orientation, flips it when the filled face and the receiving face are not
//! geometrically the same (`cxx:177-180`), skips the `Common` of a `Cut`
//! (`cxx:172-176`, `:194-197`) and adds every non-degenerated new edge of the
//! curve after making sure both faces carry a p-curve for it
//! (`cxx:207-264`).
//!
//! UNPORTED: the `OCCT_DEBUG` block (`cxx:70-86`, `:160-191`, `:239-250`,
//! `:252-264`) — trace and dump only.

use crate::abs::{Orientation, ShapeType, State};
use crate::boptools_2d;
use crate::shape::{Edge, Face, TopoShape};
use crate::topopebrep::ds::curve_iterator::CurveIterator;

use super::super::g_topo::GTopo;
use super::super::wire_edge_set::WireEdgeSet;
use super::grid::fun_aresamegeom;
use super::Builder;

impl Builder {
    /// `GFillCurveTopologyWES(F1, G1, WES)` (`cxx:60-136`).
    pub fn g_fill_curve_topology_wes(
        &mut self,
        f1: &TopoShape,
        g1: &GTopo,
        wes: &mut WireEdgeSet,
    ) {
        // `cxx:64-68`
        let (_tb1, _tb2) = g1.states_on();
        let (t1, _t2) = g1.type_();
        let shape_interf = t1;

        // `cxx:88-89`
        let fc_it = match self.my_data_structure.as_ref() {
            Some(hds) => hds.face_curves_of_shape(f1),
            None => return,
        };
        self.my_face_reference = Some(f1.clone());
        self.my_face_to_fill = Some(f1.clone());

        // `cxx:92-108`: the same-domain partners of `F1` cannot contribute a
        // new intersection edge of their own.
        let (hsd, sd_list) = match self.my_data_structure.as_ref() {
            Some(hds) => (
                hds.has_same_domain(f1, true),
                hds.same_domain(f1).to_vec(),
            ),
            None => (false, Vec::new()),
        };
        let mut a_sd_map = crate::topopebrep::ds::data_structure::IndexedShapeSet::new();
        if hsd {
            for sdf in &sd_list {
                a_sd_map.add(sdf);
            }
        }

        self.fill_curve_topology_wes_filtered(fc_it, g1, wes, shape_interf, &a_sd_map);
    }

    /// The `for` loop of the first overload (`cxx:109-133`), split out only so
    /// the iterator borrow and the `&mut self` call of the second overload do
    /// not overlap.
    fn fill_curve_topology_wes_filtered(
        &self,
        mut fc_it: CurveIterator,
        g1: &GTopo,
        wes: &mut WireEdgeSet,
        shape_interf: ShapeType,
        a_sd_map: &crate::topopebrep::ds::data_structure::IndexedShapeSet,
    ) {
        while fc_it.more() {
            let skip = if shape_interf != ShapeType::Shape {
                // `cxx:111-120`
                let (shape_before, shape_after, si) = {
                    let cur = fc_it.current();
                    let t = cur.transition();
                    (t.shape_before(), t.shape_after(), cur.support())
                };
                if shape_after != shape_interf || shape_before != shape_interf {
                    true
                } else {
                    // `cxx:121-129`
                    let ss = match self.my_data_structure.as_ref() {
                        Some(hds) => hds.shape(si, true),
                        None => None,
                    };
                    match ss {
                        Some(ss) => a_sd_map.contains(&ss),
                        None => false,
                    }
                }
            } else {
                false
            };
            if !skip {
                self.g_fill_curve_topology_wes_one(&fc_it, g1, wes);
            }
            fc_it.next();
        }
    }

    /// `GFillCurveTopologyWES(FCit, G1, WES)` (`cxx:138-268`).
    pub fn g_fill_curve_topology_wes_one(
        &self,
        fc_it: &CurveIterator,
        g1: &GTopo,
        wes: &mut WireEdgeSet,
    ) {
        // `cxx:144-147`
        if !fc_it.more() {
            return;
        }

        // `cxx:149-156`
        let (tb1, tb2) = g1.states_on();
        let conf = g1.config1();
        // `cxx:157-162`: the `DIFFORIENTED` arm is commented out in OCCT, so
        // `TB` stays `TB1` whichever configuration the `GTopo` carries.
        let tb = tb1;
        let _ = conf;

        // `cxx:158-159`
        let wesf: Face = wes.face().clone();
        let ftf = match self.my_face_to_fill.clone() {
            Some(f) => f,
            None => return,
        };

        // `cxx:172-176`
        let ope_cut = self.opec12() || self.opec21();
        let com_of_cut = ope_cut && tb1 == tb2 && tb1 == State::In;

        // `cxx:178-180`
        let transition = fc_it.current().transition();
        let mut new_ori_e = transition.orientation(tb, ShapeType::Face);
        let same_geom = fun_aresamegeom(&ftf, &wesf.0);
        if !same_geom {
            new_ori_e = complement_orientation(new_ori_e);
        }

        // `cxx:194-197`
        if com_of_cut {
            return;
        }

        // `cxx:199-203`
        let i_g = fc_it.geometry();
        let l_new_e: Vec<TopoShape> = self.new_edges(i_g).to_vec();
        for ee in l_new_e {
            // `cxx:208-212`
            if crate::brep_tool::BRepTool::is_degenerated(&Edge(ee.clone())) {
                continue;
            }
            let mut e = ee;
            e.set_orientation(new_ori_e);

            // `cxx:217-228`: the filled face needs a p-curve for the new edge.
            if !boptools_2d::has_curve_on_surface(&Edge(e.clone()), &Face(ftf.clone())) {
                let pc = fc_it.p_curve();
                let Some(pc) = pc else {
                    continue;
                };
                self.my_build_tool.p_curve(&ftf, &e, &pc);
            }

            // `cxx:229-250`: so does the face the `WireEdgeSet` belongs to.
            if !boptools_2d::has_curve_on_surface(&Edge(e.clone()), &wesf) {
                // `FC2D_CurveOnSurface` builds the p-curve when it is missing
                // (`boptools_2d.rs:57-64`); OCCT answers a null handle only
                // for a shape that is not an edge or a face, and throws
                // `Standard_ProgramError` there.
                let Ok(c2d) = boptools_2d::make_2d(&Edge(e.clone()), &wesf) else {
                    continue;
                };
                self.my_build_tool.p_curve(&wesf.0, &e, &c2d);
            }

            // `cxx:264`
            wes.add_start_element(&e);
        }
    }
}

/// `TopAbs::Complement` (`TopAbs.cxx`).
fn complement_orientation(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::External,
        Orientation::External => Orientation::Internal,
    }
}
