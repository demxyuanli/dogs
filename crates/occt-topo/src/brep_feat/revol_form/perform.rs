//! `BRepFeat_MakeRevolutionForm::Perform` (`cxx:1168-1438`).
//!
//! The revolution of `myPbase` around `myAxe`, then its trimming by the two
//! planes parallel to `myPln` at `myHeight1` / `myHeight2` from the extreme
//! profile vertex `Pt`, the descendant maps (`MajMap`), the glued faces
//! (`SetGluedFaces`) and the topological reconstruction of the base class
//! (`LFPerform`).

use std::sync::Arc;

use occt_core::gp::{GpDir, GpLin, GpPnt, GpVec};
use occt_geom::curve::Curve;
use occt_geom::line::GeomLine;

use crate::abs::ShapeType;
use crate::brep_feat::algo_api;
use crate::brep_feat::perf_selection::PerfSelection;
use crate::brep_feat::rib_slot::normal;
use crate::brep_feat::statics::{face_until, tool};
use crate::brep_feat::status_error::StatusError;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::locope::cs_intersector::CSIntersector;
use crate::locope::shape_maps::ShapeValueMap;
use crate::locope::RevolutionForm;
use crate::shape::{Face, TopoShape, Vertex};
use crate::topo_tools::{explorer_shapes, is_same};
use crate::transform;

use super::{maj_map, set_glued_faces, MakeRevolutionForm};

impl MakeRevolutionForm {
    /// `Perform()` (`cxx:1168-1438`).
    ///
    /// `Err` is one of the reference's raises: `Standard_NullObject` for the
    /// null `myPln` (`cxx:1260`), `Standard_ConstructionError` for the null
    /// normal used to place the trimming planes (`cxx:1263-1264`),
    /// `StdFail_NotDone` for a `LocOpe_RevolutionForm` / `BRepAlgoAPI_Cut`
    /// that did not build.
    pub fn perform(&mut self) -> Result<(), String> {
        // `cxx:1174-1184`.
        if self.base.my_sbase.is_none()
            || self.base.my_skface.is_none()
            || self.base.my_pbase.is_none()
        {
            self.base.my_status_error = StatusError::NotInitialized;
            self.base.not_done();
            return Ok(());
        }

        let mut pbase = self.base.my_pbase.clone().expect("checked above");
        let mut skface = self.base.my_skface.clone().expect("checked above");

        // `cxx:1185-1198`: `Pt` is the first vertex of `myPbase`. A base with
        // no vertex keeps the default `gp_Pnt`, as the reference's default
        // constructed `gp_Pnt Pt` does.
        let pt = explorer_shapes(&pbase.0, ShapeType::Vertex)
            .first()
            .map(|v| BRepTool::vertex_point_world(&Vertex(v.clone())))
            .unwrap_or_else(GpPnt::default);

        if self.my_angle2 != 0.0 {
            // `cxx:1203-1208`.
            let pbase_rot = transform::rotated(&pbase.0, &self.my_axe, -self.my_angle2);
            let old_edges = explorer_shapes(&pbase.0, ShapeType::Edge);
            let new_edges = explorer_shapes(&pbase_rot, ShapeType::Edge);

            // `cxx:1208-1222`: `myLFMap` values are edges of the profile; each
            // one is replaced by the edge at the same position of the rotated
            // profile. `ex1` walks the *old* profile, so `IsSame` matches.
            let lf_pairs: Vec<(TopoShape, TopoShape)> = self
                .base
                .my_lf_map
                .iter()
                .filter_map(|(k, v)| v.first().map(|e1| (k.clone(), e1.clone())))
                .collect();
            for (key, e1) in lf_pairs {
                for (old, new) in old_edges.iter().zip(new_edges.iter()) {
                    if is_same(old, &e1) {
                        self.base.my_lf_map.bind(&key, vec![new.clone()]);
                        break;
                    }
                }
            }

            // `cxx:1224-1242`: the same walk over `mySlface`'s first edge.
            let sf_pairs: Vec<(TopoShape, TopoShape)> = self
                .my_slface
                .iter()
                .filter_map(|(k, v)| v.first().map(|e1| (k.clone(), e1.clone())))
                .collect();
            for (key, e1) in sf_pairs {
                for (old, new) in old_edges.iter().zip(new_edges.iter()) {
                    if is_same(old, &e1) {
                        self.my_slface.bind(&key, vec![new.clone()]);
                        break;
                    }
                }
            }

            // `cxx:1243-1248`.
            pbase = Face(pbase_rot);
            skface = Face(transform::rotated(&skface.0, &self.my_axe, -self.my_angle2));
            self.base.my_pbase = Some(pbase.clone());
            self.base.my_skface = Some(skface);
        }

        // `cxx:1250-1256`.
        let mut the_form = RevolutionForm::new();
        the_form.perform(&pbase.0, &self.my_axe, self.my_angle1 + self.my_angle2);
        let mut vrai_form = the_form
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone: LocOpe_RevolutionForm".to_string())?;

        maj_map(
            &pbase.0,
            &the_form,
            &mut self.base.my_map,
            &mut self.base.my_f_shape,
            &mut self.base.my_l_shape,
        );

        // `cxx:1258`.
        self.base.my_glued_f.clear();

        // `cxx:1260-1275`: the two trimming planes.
        let pln0 = self
            .my_pln
            .clone()
            .ok_or_else(|| "Standard_NullObject: myPln".to_string())?;
        let builder = TopoBuilder::new();
        let f_plane = builder.make_face_plane(&pln0);
        let nrm = normal(&f_plane, &pt);
        let nvec = GpVec::from_xyz(nrm.xyz());
        let vec1 = nvec.multiplied_scalar(self.my_height1);
        let vec2 = nvec.multiplied_scalar(-self.my_height2);
        let pln1 = pln0.translated_vec(&vec1);
        let pln2 = pln0.translated_vec(&vec2);
        let f1_plane = builder.make_face_plane(&pln1);
        let f2_plane = builder.make_face_plane(&pln2);

        // `cxx:1276-1277`: `FaceUntil` extends both trimming faces onto the
        // basis. The reference would keep a null face and answer
        // `BadIntersect` from the failed sections; the port cannot build a
        // section on a null shape, so it answers the same status here.
        let sbase = self
            .base
            .my_sbase
            .as_ref()
            .ok_or_else(|| "Standard_NullObject: mySbase".to_string())?
            .clone();
        let (Some(f1), Some(f2)) = (
            face_until(&sbase, &f1_plane),
            face_until(&sbase, &f2_plane),
        ) else {
            self.base.my_status_error = StatusError::BadIntersect;
            self.base.not_done();
            return Ok(());
        };

        // `cxx:1279-1297`.
        let mut asi1 = CSIntersector::with_shape(&f1.0);
        let mut asi2 = CSIntersector::with_shape(&f2.0);
        let dir = GpDir::from_vec(&vec1).map_err(|e| format!("Standard_ConstructionError: {e}"))?;
        let normale: Arc<dyn Curve> = Arc::new(GeomLine::new(GpLin::from_pnt_dir(pt, dir)));
        let _ = asi1.perform_curves(&[Some(normale.clone())]);
        let _ = asi2.perform_curves(&[Some(normale)]);
        let ok = asi1.is_done()
            && asi2.is_done()
            && asi1.nb_points(1).unwrap_or(0) == 1
            && asi2.nb_points(1).unwrap_or(0) == 1;
        if !ok {
            self.base.my_status_error = StatusError::BadIntersect;
            self.base.not_done();
            return Ok(());
        }

        // `cxx:1298-1300`.
        let ori1 = asi1.point(1, 1)?.orientation();
        let ori2 = asi2.point(1, 1)?.orientation().reversed();
        let ff1 = asi1.point(1, 1)?.face().clone();
        let ff2 = asi2.point(1, 1)?.face().clone();

        // `cxx:1302-1316`: the two half-space solids that trim the form.
        let mut comp_parts: Vec<TopoShape> = Vec::new();
        if let Some(s1) = tool(&f1.0, &ff1, ori1) {
            comp_parts.push(s1.0);
        }
        if let Some(s2) = tool(&f2.0, &ff2, ori2) {
            comp_parts.push(s2.0);
        }
        let comp = builder.make_compound_of(&comp_parts).0;

        // `cxx:1318`: `BRepAlgoAPI_Cut trP(VraiForm, Comp)`.
        let mut tr_p = algo_api::cut(&vrai_form, &comp);
        tr_p.build()?;

        // `cxx:1322-1366`: the descendant management. The reference walks
        // `myMap` with an iterator while rebinding entries of the same map;
        // the port snapshots the keys and re-reads each value at its turn, so
        // a rebind of `myFShape` / `myLShape` is visible to the later entries
        // just as it is in the reference.
        let keys: Vec<TopoShape> = self.base.my_map.iter().map(|(k, _)| k.clone()).collect();
        let f_shape = self.base.my_f_shape.clone();
        let l_shape = self.base.my_l_shape.clone();
        let vrai_faces = explorer_shapes(&vrai_form, ShapeType::Face);
        for orig in keys {
            let Some(sh) = self
                .base
                .my_map
                .find(&orig)
                .and_then(|v| v.first().cloned())
            else {
                continue;
            };
            for fac in &vrai_faces {
                let thew = explorer_shapes(fac, ShapeType::Wire).first().cloned();
                if let (Some(thew), Some(fs)) = (thew.as_ref(), f_shape.as_ref()) {
                    if is_same(thew, fs) {
                        let desfaces = tr_p.modified(&f2.0);
                        self.base.my_map.bind(fs, desfaces);
                        continue;
                    }
                }
                if let (Some(thew), Some(ls)) = (thew.as_ref(), l_shape.as_ref()) {
                    if is_same(thew, ls) {
                        let desfaces = tr_p.modified(&f1.0);
                        self.base.my_map.bind(ls, desfaces);
                        continue;
                    }
                }
                if is_same(fac, &sh) && !tr_p.is_deleted(fac) {
                    let desfaces = tr_p.modified(fac);
                    if !desfaces.is_empty() {
                        self.base.my_map.bind(&orig, desfaces);
                        break;
                    }
                }
            }
        }

        // `cxx:1368-1388`: the face-to-faces map the glued-face pass reads.
        let mut sliding_map: ShapeValueMap<Vec<TopoShape>> = ShapeValueMap::new();
        for fac in &vrai_faces {
            sliding_map.bind(fac, Vec::new());
            if tr_p.is_deleted(fac) {
                continue;
            }
            let desfaces = tr_p.modified(fac);
            if !desfaces.is_empty() {
                sliding_map.bind(fac, desfaces);
            } else if let Some(list) = sliding_map.change(fac) {
                list.push(fac.clone());
            }
        }

        // `cxx:1390-1391`.
        set_glued_faces(
            &self.my_slface,
            &the_form,
            &sliding_map,
            &mut self.base.my_glued_f,
        );

        // `cxx:1393`.
        vrai_form = tr_p
            .shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone: BRepAlgoAPI_Cut".to_string())?;

        // `cxx:1395-1402`.
        if !self.base.my_glued_f.is_empty() {
            self.base.my_perf_selection = PerfSelection::NoSelection;
        } else {
            self.base.my_perf_selection = PerfSelection::SelectionSh;
        }

        // `cxx:1404-1418`: every profile edge has to survive in the map.
        for e in explorer_shapes(&pbase.0, ShapeType::Edge) {
            if !self.base.my_map.is_bound(&e) {
                self.base.my_status_error = StatusError::IncSlidFace;
                self.base.not_done();
                return Ok(());
            }
        }

        // `cxx:1420`.
        self.base.my_gshape = Some(vrai_form);

        // `cxx:1422-1436`.
        if !self.base.my_glued_f.is_empty() && self.base.my_suntil.is_some() {
            self.base.my_status_error = StatusError::InvShape;
            self.base.not_done();
            return Ok(());
        }

        // `cxx:1438`.
        self.base.lf_perform()?;
        Ok(())
    }
}
