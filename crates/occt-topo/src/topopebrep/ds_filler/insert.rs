//! `TopOpeBRep_DSFiller` insert pass.
//!
//! Source: `TKBool/TopOpeBRep/TopOpeBRep_DSFiller.cxx:427-1140` (OCCT 8.0.0):
//! `Insert` (`:427-439`), `InsertIntersection` (`:444-719`), `Insert2d`
//! (`:892-905`), `InsertIntersection2d` (`:907-996`) and `Insert1d`
//! (`:1084-1133`).
//!
//! The 3d walk (`cxx:566-671`) dispatches each intersecting pair by type:
//! FACE/FACE to `myFacesFiller`, EDGE/EDGE to `myEdgesFiller`, FACE/EDGE and
//! EDGE/FACE to `myFaceEdgeFiller`. A face/face same-domain couple is followed
//! by its edge/edge couples (`cxx:613-648`), and any later pair that is not
//! edge/edge closes that same-domain couple (`cxx:588-599`).
//!
//! The 2d walk (`cxx:930-985`) is the same domain, restricted to the
//! same-domain face couples `TopOpeBRep_ShapeIntersector2d` reports.

use crate::abs::Orientation;
use crate::shape::{Face, TopoShape};
use crate::topexp::Explorer;
use crate::topopebrep::ds::h_data_structure::HDataStructure;

use super::ds_tools::{
    brep_correct_g_bound, brep_merge_pds, brep_sort_on_parameter, brep_unfill_same_domain,
};
use super::DsFiller;

impl DsFiller {
    /// `Insert(aS1, aS2, HDS, orientFORWARD = true)` (`cxx:427-439`).
    pub fn insert(
        &mut self,
        s1: &TopoShape,
        s2: &TopoShape,
        hds: &mut HDataStructure,
        orient_forward: bool,
    ) -> Result<(), String> {
        self.insert_intersection(s1, s2, hds, orient_forward)?;
        self.complete(hds);
        Ok(())
    }

    /// `InsertIntersection(aS1, aS2, HDS, orientFORWARD)` (`cxx:444-719`).
    pub fn insert_intersection(
        &mut self,
        s1: &TopoShape,
        s2: &TopoShape,
        hds: &mut HDataStructure,
        orient_forward: bool,
    ) -> Result<(), String> {
        // `cxx:448-449`: `FBOX_Prepare()` (`TopOpeBRepTool_box.cxx:22`) and
        // `FC2D_Prepare(aS1, aS2)` (`TopOpeBRepTool_2d.cxx:96`) are UNPORTED
        // process-wide preparations, see the module note.

        // `cxx:451-454` / `cxx:910-914`
        let psc = self.ensure_p_shape_classifier();
        self.faces_filler.set_pshape_classifier(Some(psc));

        if !self.check_insert(s1, s2) {
            // `cxx:470-475`
            hds.add_ancestors(s1);
            hds.add_ancestors(s2);
            self.fdscnx_prepare(s1, s2, hds);
            return Ok(());
        }

        let mut s1 = s1.clone();
        let mut s2 = s2.clone();
        if orient_forward {
            // `cxx:479-489`
            if s1.orientation() == Orientation::Reversed {
                s1 = s1.oriented(Orientation::Forward);
            }
            if s2.orientation() == Orientation::Reversed {
                s2 = s2.oriented(Orientation::Forward);
            }
        }

        // `cxx:491-504`
        hds.change_ds().add_shape_with_ancestor(&s1, 1);
        hds.change_ds().add_shape_with_ancestor(&s2, 2);
        let mut ex1 = Explorer::new(&s1, crate::abs::ShapeType::Solid);
        while ex1.more() {
            let so1 = ex1.current().clone();
            let mut ex2 = Explorer::new(&s2, crate::abs::ShapeType::Solid);
            while ex2.more() {
                let so2 = ex2.current().clone();
                hds.change_ds()
                    .fill_shapes_same_domain(&so1, &so2, true)?;
                ex2.next();
            }
            ex1.next();
        }

        // `cxx:509-519`
        let mut ff_samdom = true;
        let mut isl_ff_samdom = false;
        let mut is_ff_samdom = false;
        let mut is_ee = false;
        let mut lff1: Option<TopoShape> = None;
        let mut lff2: Option<TopoShape> = None;

        // `cxx:520-528`: the rejected-shape maps are cleared before the search.
        hds.change_ds().change_map_of_rejected_shapes_obj().clear();
        hds.change_ds().change_map_of_rejected_shapes_tool().clear();

        // `cxx:530-547`: rejected faces of the object and of the tool.
        let mut list_obj: Vec<TopoShape> = Vec::new();
        let mut list_tool: Vec<TopoShape> = Vec::new();
        if has_solids(&s2) {
            self.shape_intersector.rejected_faces(&s1, &s2, &mut list_obj)?;
        }
        if has_solids(&s1) {
            self.shape_intersector.rejected_faces(&s2, &s1, &mut list_tool)?;
        }

        // `cxx:550-563`
        for s in &list_obj {
            hds.change_ds().change_map_of_rejected_shapes_obj().add(s);
        }
        for s in &list_tool {
            hds.change_ds().change_map_of_rejected_shapes_tool().add(s);
        }

        self.shape_intersector.init_intersection(&s1, &s2)?;
        while self.shape_intersector.more_intersection() {
            // `cxx:580-586`
            let gs1 = self
                .shape_intersector
                .current_geom_shape(1)
                .ok_or_else(|| "CurrentGeomShape : no intersection".to_string())?;
            let gs2 = self
                .shape_intersector
                .current_geom_shape(2)
                .ok_or_else(|| "CurrentGeomShape : no intersection".to_string())?;
            let t1 = gs1.shape_type();
            let t2 = gs2.shape_type();

            // `cxx:588-599`
            let is_ff = t1 == crate::abs::ShapeType::Face && t2 == crate::abs::ShapeType::Face;
            is_ff_samdom = false;
            is_ee = t1 == crate::abs::ShapeType::Edge && t2 == crate::abs::ShapeType::Edge;
            let mut unfill = !is_ee && isl_ff_samdom;

            if unfill {
                let (Some(f1), Some(f2)) = (lff1.clone(), lff2.clone()) else {
                    return Err("BREP_UnfillSameDomain : null face".to_string());
                };
                let psc = self.ensure_p_shape_classifier();
                let mut classifier = psc.borrow_mut();
                // `cxx:597`: the answer is not read again in OCCT either.
                let _ = brep_unfill_same_domain(&f1, &f2, hds, &mut classifier)?;
            }

            if is_ff {
                // `cxx:601-627`
                let ff = self.shape_intersector.change_faces_intersector();
                is_ff_samdom = ff.same_domain()?;
                ff_samdom = ff_samdom && is_ff_samdom;
                if is_ff_samdom {
                    lff1 = Some(ff.face(1)?.0.clone());
                    lff2 = Some(ff.face(2)?.0.clone());
                }
                let ff = self.shape_intersector.change_faces_intersector();
                self.faces_filler.insert(&gs1, &gs2, ff, hds)?;
            } else if t1 == crate::abs::ShapeType::Edge && t2 == crate::abs::ShapeType::Edge {
                // `cxx:629-648`
                let ee = self.shape_intersector.change_edges_intersector();
                ee.set_dimension(2);
                let _ee_samedomain = ee.same_domain();
                if isl_ff_samdom {
                    let (Some(f1), Some(f2)) = (lff1.clone(), lff2.clone()) else {
                        return Err("EdgesFiller : null face".to_string());
                    };
                    self.edges_filler.set_face(1, &f1)?;
                    self.edges_filler.set_face(2, &f2)?;
                }
                let ee = self.shape_intersector.change_edges_intersector();
                self.edges_filler.insert(&gs1, &gs2, ee, hds)?;
            } else if t1 == crate::abs::ShapeType::Face && t2 == crate::abs::ShapeType::Edge {
                // `cxx:650-658`
                let fe = self.shape_intersector.change_face_edge_intersector();
                self.face_edge_filler.insert(&gs1, &gs2, fe, hds)?;
            } else if t1 == crate::abs::ShapeType::Edge && t2 == crate::abs::ShapeType::Face {
                // `cxx:660-667`
                let fe = self.shape_intersector.change_face_edge_intersector();
                self.face_edge_filler.insert(&gs2, &gs1, fe, hds)?;
            }

            // `cxx:669-670`
            isl_ff_samdom = is_ff_samdom;
            self.shape_intersector.next_intersection()?;
        }

        // `cxx:674-688`
        if ff_samdom {
            hds.change_ds();
            crate::topopebrep::ds::export::fun_ds_fill_sdm_faces(hds)?;
        } else {
            crate::topopebrep::ds::export::fun_ds_add_sesdm1d(hds)?;
        }
        if ff_samdom {
            hds.change_ds().set_isfafa(true);
        }

        // `cxx:690-697`
        if isl_ff_samdom && !is_ee {
            let (Some(f1), Some(f2)) = (lff1.clone(), lff2.clone()) else {
                return Err("BREP_UnfillSameDomain : null face".to_string());
            };
            let psc = self.ensure_p_shape_classifier();
            let mut classifier = psc.borrow_mut();
            brep_unfill_same_domain(&f1, &f2, hds, &mut classifier)?;
        }

        // `cxx:699-702`
        brep_sort_on_parameter(hds);
        brep_correct_g_bound(hds);
        brep_merge_pds(hds)?;

        // `cxx:704-706`
        hds.add_ancestors(&s1);
        hds.add_ancestors(&s2);
        self.fdscnx_prepare(&s1, &s2, hds);
        self.fdssdm_prepare(hds);
        Ok(())
    }

    /// `Insert2d(aS1, aS2, HDS)` (`cxx:892-905`).
    pub fn insert_2d(
        &mut self,
        s1: &TopoShape,
        s2: &TopoShape,
        hds: &mut HDataStructure,
    ) -> Result<(), String> {
        self.insert_intersection_2d(s1, s2, hds)?;
        self.complete_ds_2d(hds);
        Ok(())
    }

    /// `InsertIntersection2d(aS1, aS2, HDS)` (`cxx:907-996`).
    pub fn insert_intersection_2d(
        &mut self,
        s1: &TopoShape,
        s2: &TopoShape,
        hds: &mut HDataStructure,
    ) -> Result<(), String> {
        // `cxx:910-914`
        // `cxx:451-454` / `cxx:910-914`
        let psc = self.ensure_p_shape_classifier();
        self.faces_filler.set_pshape_classifier(Some(psc));

        // `cxx:917-920`
        if !self.clear_shape_same_domain(s1, s2, hds) {
            return Ok(());
        }

        // `cxx:922-926`
        hds.change_ds().add_shape_with_ancestor(s1, 1);
        hds.change_ds().add_shape_with_ancestor(s2, 2);

        let mut lff1: Option<TopoShape> = None;
        let mut lff2: Option<TopoShape> = None;
        let mut is_ff_samdom = false;

        self.shape_intersector2d.init_intersection(s1, s2)?;
        while self.shape_intersector2d.more_intersection() {
            // `cxx:935-938`
            let gs1 = self
                .shape_intersector2d
                .current_geom_shape(1)
                .ok_or_else(|| "CurrentGeomShape : no intersection 2d".to_string())?;
            let gs2 = self
                .shape_intersector2d
                .current_geom_shape(2)
                .ok_or_else(|| "CurrentGeomShape : no intersection 2d".to_string())?;
            let t1 = gs1.shape_type();
            let t2 = gs2.shape_type();

            // `cxx:942-967`
            let is_ff = t1 == crate::abs::ShapeType::Face && t2 == crate::abs::ShapeType::Face;
            let is_ee = t1 == crate::abs::ShapeType::Edge && t2 == crate::abs::ShapeType::Edge;
            let mut unfill = !is_ee && is_ff_samdom;
            if unfill {
                let (Some(f1), Some(f2)) = (lff1.clone(), lff2.clone()) else {
                    return Err("BREP_UnfillSameDomain : null face".to_string());
                };
                let psc = self.ensure_p_shape_classifier();
                let mut classifier = psc.borrow_mut();
                let st1 = classifier.state_shape_shape(&f1, &f2, 1)?;
                let st2 = classifier.state_shape_shape(&f2, &f1, 1)?;
                unfill = (st1 == crate::abs::State::Out && st2 == crate::abs::State::Out)
                    || (st1 == crate::abs::State::Unknown && st2 == crate::abs::State::Unknown);
                drop(classifier);
                if unfill {
                    hds.change_ds().unfill_shapes_same_domain(&f1, &f2);
                }
            }

            if is_ff {
                // `cxx:969-976`
                is_ff_samdom = true;
                {
                    let ee = self.shape_intersector2d.change_edges_intersector();
                    ee.set_faces(&gs1, &gs2);
                }
                lff1 = Some(gs1.clone());
                lff2 = Some(gs2.clone());
                hds.change_ds()
                    .fill_shapes_same_domain(&gs1, &gs2, true)?;
            } else if t1 == crate::abs::ShapeType::Edge && t2 == crate::abs::ShapeType::Edge {
                // `cxx:978-984`
                let ee = self.shape_intersector2d.change_edges_intersector();
                ee.set_dimension(2);
                let (Some(f1), Some(f2)) = (lff1.clone(), lff2.clone()) else {
                    return Err("EdgesFiller : null face".to_string());
                };
                self.edges_filler.set_face(1, &f1)?;
                self.edges_filler.set_face(2, &f2)?;
                let ee = self.shape_intersector2d.change_edges_intersector();
                self.edges_filler.insert(&gs1, &gs2, ee, hds)?;
            }
            self.shape_intersector2d.next_intersection()?;
        }

        // `cxx:988-990`
        brep_sort_on_parameter(hds);
        brep_correct_g_bound(hds);
        brep_merge_pds(hds)?;
        Ok(())
    }

    /// `Insert1d(aS1, aS2, aF1, aF2, HDS, orientFORWARD = false)`
    /// (`cxx:1084-1133`).
    pub fn insert_1d(
        &mut self,
        s1: &TopoShape,
        s2: &TopoShape,
        f1: &Face,
        f2: &Face,
        hds: &mut HDataStructure,
        orient_forward: bool,
    ) -> Result<(), String> {
        // `cxx:1090-1093`
        if !self.check_insert(s1, s2) {
            return Ok(());
        }

        let mut s1 = s1.clone();
        let mut s2 = s2.clone();
        if orient_forward {
            // `cxx:1095-1105`
            if s1.orientation() == Orientation::Reversed {
                s1 = s1.oriented(Orientation::Forward);
            }
            if s2.orientation() == Orientation::Reversed {
                s2 = s2.oriented(Orientation::Forward);
            }
        }

        // `cxx:1108-1110`
        hds.change_ds().add_shape_with_ancestor(&s1, 1);
        hds.change_ds().add_shape_with_ancestor(&s2, 2);

        // `cxx:1114-1127`
        self.shape_intersector
            .init_intersection_with_faces(&s1, &s2, f1, f2)?;
        while self.shape_intersector.more_intersection() {
            let gs1 = self
                .shape_intersector
                .current_geom_shape(1)
                .ok_or_else(|| "CurrentGeomShape : no intersection".to_string())?;
            let gs2 = self
                .shape_intersector
                .current_geom_shape(2)
                .ok_or_else(|| "CurrentGeomShape : no intersection".to_string())?;
            let t1 = gs1.shape_type();
            let t2 = gs2.shape_type();
            if t1 == crate::abs::ShapeType::Edge && t2 == crate::abs::ShapeType::Edge {
                let ee = self.shape_intersector.change_edges_intersector();
                ee.set_dimension(1);
                let ee = self.shape_intersector.change_edges_intersector();
                self.edges_filler.insert(&gs1, &gs2, ee, hds)?;
            }
            self.shape_intersector.next_intersection()?;
        }

        // `cxx:1129-1130`: update wires
        self.complete_ds(hds);
        Ok(())
    }
}

/// `TopExp::MapShapes(S, TopAbs_SOLID)` emptiness test (`cxx:533-547`).
fn has_solids(s: &TopoShape) -> bool {
    let mut ex = Explorer::new(s, crate::abs::ShapeType::Solid);
    ex.more()
}
