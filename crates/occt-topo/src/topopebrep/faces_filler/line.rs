//! `TopOpeBRep_FacesFiller_1.cxx` — the per-line filling of the DS:
//! `ProcessLine`, `FillLine`, `AddShapesLine` and their geometry helpers.
//!
//! Sources (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRep`):
//! - `ProcessLine` (`TopOpeBRep_FacesFiller_1.cxx:460-487`)
//! - `ResetDSC` (`_1.cxx:493-497`)
//! - `ProcessVPnotonR` (`_1.cxx:507-577`)
//! - `FillLineVPonR` (`_1.cxx:1157-1217`)
//! - `FillLine` (`_1.cxx:1222-1287`)
//! - `AddShapesLine` (`_1.cxx:1292-1421`)
//! - `StoreCurveInterference` (`_1.cxx:1425-1436`)
//! - `GetGeometry` (`_1.cxx:1438-1448`)
//! - `MakeGeometry` (`_1.cxx:1452-1472`)
//! - `GetFFGeometry` (`_1.cxx:1476-1507`)
//!
//! UNPORTED, called from `FillLineVPonR`:
//! - `FUN_GetdgData` (`TopOpeBRep_vprdeg.cxx:179-379`) and `FUN_FillVof12`
//!   (`vprdeg.cxx:105-127`), the degenerated-edge data collection. The port
//!   leaves `data_for_degen_ed` empty, which is exactly the arm OCCT takes when
//!   a line touches no degenerated edge (`vprdeg.cxx:676-680`).
//! - `ProcessRLine` (`_1.cxx:647-1125`), in [`super::rline`].
//! - `ProcessVPR` (`_1.cxx:581-645`) is ported in [`super::vpr`]; it reaches
//!   `ProcessVPondgE` (`vprdeg.cxx:633-801`) through `ProcessVPonR`, which is
//!   the remaining gate.

use occt_geom2d::curve::Curve2d;
use std::sync::Arc;

use occt_core::precision::PCONFUSION;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::shape::{Face, TopoShape};
use crate::topopebrep::ds::curve::Curve as DsCurve;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference::InterferenceHandle;
use crate::topopebrep::ds::interference_tool::make_curve_interference;
use crate::topopebrep::ds::interference_tool::make_face_curve_interference;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::point_iterator::PointIterator;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::faces_intersector::FacesIntersector;
use crate::topopebrep::geom_tool::GeomTool;
use crate::topopebrep::interf::point_geom_tool::PointGeomTool;
use crate::topopebrep::interf::v_point_inter::VPointInter;
use crate::topopebrep::type_line_curve::TypeLineCurve;

/// `FUN_haslastvpon0(L)` (`_1.cxx:1129-1153`): a walking line whose last
/// in/on point has shape index 0.
fn fun_has_last_vp_on0(l: &crate::topopebrep::line_inter::LineInter) -> bool {
    if l.type_line_curve() != TypeLineCurve::Walking {
        return false;
    }
    let nb = l.nb_v_point();
    let mut last = 0;
    for k in 1..=nb {
        if let Some(vp) = l.v_point(k) {
            if vp.keep() {
                last = k;
            }
        }
    }
    let _ = last;
    // `VPBounds(iINON1, iINONn, nINON)` (`cxx:1137-1138`), then the loop looks
    // for `iVP == iINONn && absindex == 0`.
    let mut bounds = l.clone();
    let Some((_f, i_inonn, _n)) = bounds.vp_bounds() else {
        return false;
    };
    match l.v_point(i_inonn) {
        Some(vp) => vp.index() == i_inonn && vp.shape_index() == 0,
        None => false,
    }
}

impl super::FacesFiller {
    /// `ResetDSC()` (`_1.cxx:493-497`).
    fn reset_dsc(&mut self) {
        self.dsc_index = 0;
        self.dscil.clear();
    }

    /// `ProcessLine()` (`_1.cxx:460-487`): fill the DS curve for the current
    /// line, then hand its interferences to `AddShapesLine`.
    pub fn process_line(
        &mut self,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
    ) -> Result<(), String> {
        let reject = !self.line_ok || self.line_index <= 0;
        if reject {
            return Ok(());
        }
        self.reset_dsc();

        let (has_v_pon_r, inl) = {
            let l = facint.change_line(self.line_index);
            (l.has_v_pon_r(), l.inl())
        };
        if has_v_pon_r {
            self.fill_line_vp_on_r(facint, hds)?;
        } else {
            self.fill_line(facint, hds)?;
        }
        if inl {
            return Ok(());
        }

        hds.sort_on_parameter_in_place(&mut self.dscil)?;
        self.add_shapes_line(facint, hds)
    }

    /// `FillLine()` (`_1.cxx:1222-1287`): one curve/point interference per kept
    /// VPoint of the line.
    fn fill_line(
        &mut self,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
    ) -> Result<(), String> {
        let (i_inon1, i_inonn, n_inon) = {
            let l = facint.change_line(self.line_index);
            l.vp_bounds().unwrap_or((0, 0, 0))
        };
        if n_inon == 0 {
            return Ok(());
        }

        let nb = facint.change_line(self.line_index).nb_v_point();
        for k in 1..=nb {
            let vp = match facint.change_line(self.line_index).v_point(k) {
                Some(vp) => vp.clone(),
                None => continue,
            };
            if !vp.keep() {
                continue;
            }

            let mut it = PointIterator::new(self.dscil.clone());
            let cpifound = self.get_geometry(hds, &mut it, &vp)?.is_some();
            let (pv_index, pv_kind) = if cpifound {
                let i = it.value().read().expect("interference lock");
                (i.geometry(), i.geometry_type())
            } else {
                match self.get_ff_geometry(hds.ds(), &vp)? {
                    Some((g, kd)) => (g, kd),
                    None => Self::make_geometry(hds.change_ds(), &vp, 0)?,
                }
            };

            let trans_line = if !cpifound {
                if k == i_inon1 {
                    Transition::from_orientation(Orientation::Forward)
                } else if k == i_inonn {
                    Transition::from_orientation(Orientation::Reversed)
                } else {
                    Transition::default()
                }
            } else {
                it.current().transition().complement()
            };

            let parline = vp.parameter_on_line()?;
            let cpi =
                make_curve_interference(&trans_line, Kind::Curve, 0, pv_kind, pv_index, parline);
            self.store_curve_interference(hds, &cpi);
        }
        Ok(())
    }

    /// `FillLineVPonR()` (`_1.cxx:1157-1217`).
    ///
    /// The `FUN_GetdgData` / `FUN_FillVof12` prologue (`_1.cxx:1166-1167`) is
    /// [`super::degen_data`]; `ProcessVPondgE` reads the map it fills and
    /// returns `false` while it is empty (`vprdeg.cxx:676-680`).
    fn fill_line_vp_on_r(
        &mut self,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
    ) -> Result<(), String> {
        self.kept_vp_nbr = 0;

        // `cxx:1166-1167`
        let (Some(f1), Some(f2)) = (self.f1.clone(), self.f2.clone()) else {
            return Err("TopOpeBRep_FacesFiller::FillLineVPonR : no faces".to_string());
        };
        super::degen_data::fun_get_dg_data(
            hds,
            facint,
            self.line_index,
            &f1.0,
            &f2.0,
            &mut self.data_for_degen_ed,
        )?;
        super::degen_data::fun_fill_v_of12(hds, facint, self.line_index)?;

        let lrest = facint.change_line(self.line_index).type_line_curve()
            == TypeLineCurve::Restriction;
        if lrest {
            self.process_r_line(facint, hds)?;
            return Ok(());
        }

        let (i_inon1, i_inonn, n_inon) = self.line_vp_bounds(facint);
        let _ = (i_inon1, i_inonn);
        if n_inon == 0 {
            return Ok(());
        }

        // `cxx:1189`: `myLineIsonEdge = LSameDomainERL(*myLine, myERL)`.
        {
            let erl = self.erl.clone();
            let global_tol_ff = facint.global_tol_ff();
            let mut l = facint.change_line(self.line_index).clone();
            self.line_is_on_edge =
                super::process_gr::lsame_domain_erl(&mut l, &erl, global_tol_ff)?;
        }
        // `cxx:1192`: `myLastVPison0 = ::FUN_haslastvpon0(*myLine)`.
        self.last_vp_is_on0 = {
            let l = facint.change_line(self.line_index).clone();
            super::line::fun_has_last_vp_on0(&l)
        };

        let nb = facint.change_line(self.line_index).nb_v_point();
        for k in 1..=nb {
            let Some(vp) = facint.change_line(self.line_index).v_point(k).cloned() else {
                continue;
            };
            self.process_vpr(facint, hds, &vp)?;
        }

        if self.line_is_on_edge && !self.dscil.is_empty() {
            self.dscil.clear();
        }
        Ok(())
    }

    /// `AddShapesLine()` (`_1.cxx:1292-1421`): turns `myDSCIL` into a DS curve
    /// with its two face-curve interferences.
    fn add_shapes_line(
        &mut self,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
    ) -> Result<(), String> {
        if self.dscil.is_empty() {
            return Ok(());
        }
        if facint.change_line(self.line_index).inl() {
            return Ok(());
        }

        let f1 = self
            .f1
            .as_ref()
            .ok_or_else(|| "TopOpeBRep_FacesFiller::AddShapesLine : no F1".to_string())?
            .clone();
        let f2 = self
            .f2
            .as_ref()
            .ok_or_else(|| "TopOpeBRep_FacesFiller::AddShapesLine : no F2".to_string())?
            .clone();

        let dsc_index = self.dsc_index;
        let i_f1 = hds.change_ds().add_shape_with_ancestor(&f1.0, 1);
        let i_f2 = hds.change_ds().add_shape_with_ancestor(&f2.0, 2);

        let (pmin, pmax) = hds
            .min_max_on_parameter(&self.dscil)?
            .ok_or_else(|| "TopOpeBRep_FacesFiller::AddShapesLine : empty DSCIL".to_string())?;

        let d = (pmin - pmax).abs();
        let isper = facint.change_line(self.line_index).is_periodic();
        let mut id = (d <= PCONFUSION) && !isper;

        let wline =
            facint.change_line(self.line_index).type_line_curve() == TypeLineCurve::Walking;
        let vclosed = facint.change_line(self.line_index).is_v_closed();
        if wline && !isper && vclosed {
            // xpu240399 : USA60298 : avoid creating curve
            let first = self
                .dscil
                .first()
                .ok_or_else(|| "TopOpeBRep_FacesFiller::AddShapesLine : empty DSCIL".to_string())?;
            let (ipf, kpf) = {
                let i = first.read().expect("interference lock");
                (i.geometry(), i.geometry_type())
            };
            let (ptf, tolf) = if kpf == Kind::Point {
                let pf = hds
                    .ds()
                    .point(ipf)
                    .ok_or_else(|| "TopOpeBRep_FacesFiller::AddShapesLine : no point".to_string())?;
                (*pf.point(), pf.tolerance())
            } else {
                let vf = hds
                    .ds()
                    .shape(ipf, true)
                    .ok_or_else(|| "TopOpeBRep_FacesFiller::AddShapesLine : no vertex".to_string())?;
                let vf = crate::shape::Vertex(vf);
                (BRepTool::vertex_point(&vf), BRepTool::vertex_tolerance(&vf))
            };

            let last = self
                .dscil
                .last()
                .ok_or_else(|| "TopOpeBRep_FacesFiller::AddShapesLine : empty DSCIL".to_string())?;
            let (ipl, kpl) = {
                let i = last.read().expect("interference lock");
                (i.geometry(), i.geometry_type())
            };
            let toll = if kpl == Kind::Point {
                hds.ds()
                    .point(ipl)
                    .map(|pl| pl.tolerance())
                    .ok_or_else(|| {
                        "TopOpeBRep_FacesFiller::AddShapesLine : no point".to_string()
                    })?
            } else {
                let vl = hds.ds().shape(ipl, true).ok_or_else(|| {
                    "TopOpeBRep_FacesFiller::AddShapesLine : no vertex".to_string()
                })?;
                BRepTool::vertex_tolerance(&crate::shape::Vertex(vl))
            };

            let tol = tolf.max(toll);
            let mut onsampt = true;
            let nbwp = facint.change_line(self.line_index).nb_w_point();
            for ii in 1..=nbwp {
                let pp = facint.change_line(self.line_index).w_point(ii).value();
                if !pp.is_equal_with_tolerance(&ptf, tol) {
                    onsampt = false;
                    break;
                }
            }
            if onsampt {
                id = true;
            }
        }

        if id {
            hds.change_ds().change_curve(dsc_index).change_keep(false);
            return Ok(());
        }

        // TopOpeBRep_GeomTool::MakeCurves(pmin, pmax, (*myLine), myF1, myF2, DSC,
        //                                 PC1, PC2) (`cxx:1381`).
        let mut pc1: Option<Arc<dyn Curve2d>> = None;
        let mut pc2: Option<Arc<dyn Curve2d>> = None;
        {
            let mut dsc = hds.change_ds().change_curve(dsc_index);
            let l = facint.change_line(self.line_index);
            GeomTool::make_curves(pmin, pmax, l, &mut dsc, &mut pc1, &mut pc2)?;
        }

        // Patch: avoid making too small edges (`cxx:1384-1398`).
        {
            let dsc = hds.change_ds().change_curve(dsc_index);
            let (fpar, lpar) = dsc
                .range()
                .ok_or_else(|| "TopOpeBRep_FacesFiller::AddShapesLine : no range".to_string())?;
            let curve = dsc
                .curve()
                .cloned()
                .ok_or_else(|| "TopOpeBRep_FacesFiller::AddShapesLine : no curve".to_string())?;
            let the_box = occt_geom::geom_bnd_lib::box_curve(curve.as_ref(), fpar, lpar, 0.0);
            let (xmin, ymin, zmin, xmax, ymax, zmax) = the_box.get().ok_or_else(|| {
                "TopOpeBRep_FacesFiller::AddShapesLine : void box".to_string()
            })?;
            let max_side = (xmax - xmin).max(ymax - ymin).max(zmax - zmin);
            let min_tol = BRepTool::face_tolerance(&f1).min(BRepTool::face_tolerance(&f2));
            if max_side < min_tol {
                hds.change_ds().change_curve(dsc_index).change_keep(false);
                return Ok(());
            }
        }

        hds.change_ds().change_curve(dsc_index).set_tolerance(1.0e-8);

        let t1 = self.face_face_transition(facint, 1)?.clone();
        let t2 = self.face_face_transition(facint, 2)?.clone();

        hds.change_ds()
            .change_curve_interferences(dsc_index)
            .extend(self.dscil.iter().cloned());

        let mut t1i = t1;
        t1i.set_index(i_f2);
        let fci1 = make_face_curve_interference(&t1i, i_f2, dsc_index, pc1);
        hds.store_interference_on_shape(&fci1, &f1.0)?;

        let mut t2i = t2;
        t2i.set_index(i_f1);
        let fci2 = make_face_curve_interference(&t2i, i_f1, dsc_index, pc2);
        hds.store_interference_on_shape(&fci2, &f2.0)?;

        let dsc = hds.change_ds().change_curve(dsc_index);
        dsc.set_shapes(Some(f1.0.clone()), Some(f2.0.clone()));
        dsc.set_sci(Some(fci1), Some(fci2));
        Ok(())
    }

    /// `StoreCurveInterference(I)` (`_1.cxx:1425-1436`).
    pub(super) fn store_curve_interference(&mut self, hds: &mut HDataStructure, i: &InterferenceHandle) {
        if self.dsc_index == 0 {
            let mut dsc = DsCurve::new();
            self.dsc_index = hds.change_ds().add_curve(&mut dsc);
        }
        i.write().expect("interference lock").set_support(self.dsc_index);
        hds.store_interference_in_list(i, &mut self.dscil);
    }

    /// `GetGeometry(IT, VP, G, K)` (`_1.cxx:1438-1448`).
    pub(super) fn get_geometry(
        &self,
        hds: &HDataStructure,
        it: &mut PointIterator,
        vp: &VPointInter,
    ) -> Result<Option<(i32, Kind)>, String> {
        let dsp = PointGeomTool::make_point_of_v_point_inter(vp)?;
        // `myHDS->GetGeometry(IT, DSP, G, K)` leaves the iterator on the hit.
        hds.get_geometry(it, &dsp)
    }

    /// `MakeGeometry(VP, ShapeIndex, K)` (`_1.cxx:1452-1472`).
    pub(super) fn make_geometry(
        ds: &mut crate::topopebrep::ds::data_structure::DataStructure,
        vp: &VPointInter,
        shape_index: i32,
    ) -> Result<(i32, Kind), String> {
        if vp.is_vertex(shape_index)? {
            let s = vp
                .vertex(shape_index)
                .ok_or_else(|| "TopOpeBRep_FacesFiller::MakeGeometry : no vertex".to_string())?;
            let g = ds.add_shape_with_ancestor(s, shape_index);
            Ok((g, Kind::Vertex))
        } else {
            let p = PointGeomTool::make_point_of_v_point_inter(vp)?;
            let g = ds.add_point(&p);
            Ok((g, Kind::Point))
        }
    }

    /// `GetFFGeometry(VP, K, G)` (`_1.cxx:1500-1507`), through
    /// `GetFFGeometry(DSP, K, G)` (`_1.cxx:1476-1497`).
    pub(super) fn get_ff_geometry(
        &self,
        ds: &crate::topopebrep::ds::data_structure::DataStructure,
        vp: &VPointInter,
    ) -> Result<Option<(i32, Kind)>, String> {
        let dsp = PointGeomTool::make_point_of_v_point_inter(vp)?;
        let n = ds.nb_points();
        let mut i = self.ff_first_dsp;
        while i <= n {
            if let Some(oodsp) = ds.point(i) {
                if PointGeomTool::is_equal(&dsp, oodsp) {
                    return Ok(Some((i, Kind::Point)));
                }
            }
            i += 1;
        }
        Ok(None)
    }
}
