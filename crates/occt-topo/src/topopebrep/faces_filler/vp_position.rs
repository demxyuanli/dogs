//! `TopOpeBRep_FacesFiller`'s VPoint positioning passes.
//!
//! Sources (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRep`):
//! - `VP_Position(FACINT)` (`TopOpeBRep_FacesFiller.cxx:377-401`)
//! - `VP_Position(L)` (`cxx:405-418`)
//! - `VP_PositionOnL(L)` (`cxx:426-459`)
//! - `VP_PositionOnR(L)` (`cxx:467-560`)
//! - `VP_Position(VP)` (`cxx:565-688`)
//! - `PequalVPonR(...)` (`cxx:600-660`)
//! - `Lminmax(L, pmin, pmax)` (`TopOpeBRep_ProcessGR.cxx:233-255`)
//! - `IsVPtransLok(L, iVP, SI12, T)` (`TopOpeBRep_ProcessGR.cxx:320-339`)
//! - `TransvpOK(L, ivp, SI, isINOUT)` (`ProcessGR.cxx:341-366`)
//! - `FUN_projPonL(P, L, FF, paramL)` (`TopOpeBRep_FacesFiller.cxx:45-74`)

use occt_core::gp::GpPnt;
use occt_core::precision::PCONFUSION;

use super::FacesFiller;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::faces_intersector::FacesIntersector;
use crate::topopebrep::ff_transition_tool::FFTransitionTool;
use crate::topopebrep::interf::v_point_inter::VPointInter;
use crate::topopebrep::interf::v_point_inter_classifier::VPointInterClassifier;
use crate::topopebrep::line_inter::LineInter;
use crate::topopebrep::type_line_curve::TypeLineCurve;

impl FacesFiller {
    /// `VP_Position(FACINT)` (`cxx:377-401`): restriction lines first, then the
    /// others.
    pub fn vp_position_all(&mut self, facint: &mut FacesIntersector) -> Result<(), String> {
        facint.init_line();
        while facint.more_line() {
            let i = facint.current_line_index();
            if facint.change_line(i).type_line_curve() == TypeLineCurve::Restriction {
                self.vp_position(facint, i)?;
            }
            facint.next_line();
        }

        facint.init_line();
        while facint.more_line() {
            let i = facint.current_line_index();
            if facint.change_line(i).type_line_curve() != TypeLineCurve::Restriction {
                self.vp_position(facint, i)?;
            }
            facint.next_line();
        }
        Ok(())
    }

    /// `VP_Position(L)` (`cxx:407-421`).
    pub fn vp_position(
        &mut self,
        facint: &mut FacesIntersector,
        index: i32,
    ) -> Result<(), String> {
        self.line_index = index;
        let isrest =
            facint.change_line(index).type_line_curve() == TypeLineCurve::Restriction;
        if !isrest {
            self.vp_position_on_l(facint, index)?;
        } else {
            self.vp_position_on_r(facint, index)?;
        }
        facint.change_line(index).set_vp_bounds();
        Ok(())
    }

    /// `VP_PositionOnL(L)` (`cxx:426-459`).
    ///
    /// `PequalVPonR` may write the current VPoint from the other line, so the
    /// two lines are borrowed through
    /// [`FacesIntersector::two_lines_mut`].
    fn vp_position_on_l(
        &mut self,
        facint: &mut FacesIntersector,
        index: i32,
    ) -> Result<(), String> {
        let mut vpc = VPointInterClassifier::new();
        let nb = facint.change_line(index).nb_v_point();
        let nlines = facint.nb_lines();

        for k in 1..=nb {
            let (p3d, vpsi) = {
                let l = facint.change_line(index);
                match l.v_point(k) {
                    Some(vp) => (*vp.value()?, vp.shape_index()),
                    None => continue,
                }
            };

            let mut vp_equal_vp_on_restriction = false;
            for iol in 1..=nlines {
                if iol == index {
                    continue;
                }
                let (Some(l), Some(ol)) = facint.two_lines_mut(index, iol) else {
                    continue;
                };
                let Some(vp) = l.change_v_point(k) else { continue };
                vp_equal_vp_on_restriction =
                    self.pequal_vp_on_r(&p3d, vpsi, vp, ol)?;
                if vp_equal_vp_on_restriction {
                    break;
                }
            }

            if !vp_equal_vp_on_restriction {
                let (tol1, tol2) = facint.get_tolerances();
                let tol = if tol1 > tol2 { tol1 } else { tol2 };
                let assume_in_on = facint.change_line(index).type_line_curve()
                    != TypeLineCurve::Restriction;
                let l = facint.change_line(index);
                let Some(vp) = l.change_v_point(k) else { continue };
                self.vp_position_vp(vp, &mut vpc, assume_in_on, tol)?;
            }
        }
        Ok(())
    }

    /// `VP_PositionOnR(L)` (`cxx:466-555`).
    ///
    /// The `BRepAdaptor_Curve2d` arm (`cxx:479-491`) and the 2D projection of
    /// `FUN_projPonL` (`TopOpeBRep_FacesFiller.cxx:59-70`) are UNPORTED: they
    /// need `FC2D_HasC3D`'s p-curve build plus `FUN_tool_projPonC2D`
    /// (`TopOpeBRepTool_PROJECT.cxx:146-224`). A restriction arc without a 3D
    /// curve is therefore treated as `isline = false` and its points are
    /// unkept, which is the `VP.ChangeKeep(false)` arm of `cxx:532-535`.
    fn vp_position_on_r(
        &mut self,
        facint: &mut FacesIntersector,
        index: i32,
    ) -> Result<(), String> {
        let mut vpc = VPointInterClassifier::new();

        let (e_si, oo_e_si) = {
            let l = facint.change_line(index);
            if l.arc_is_edge(1) {
                (1, 2)
            } else {
                (2, 1)
            }
        };
        let earc = facint
            .change_line(index)
            .arc()
            .cloned()
            .ok_or_else(|| "TopOpeBRep_FacesFiller::VP_PositionOnR : no arc".to_string())?;

        let isline = match crate::brep_tool::BRepTool::edge_curve(&crate::shape::Edge(
            earc.clone(),
        )) {
            Some(c) => crate::topopebrep::tool::geometry::fun_tool_line(c.as_ref()),
            // UNPORTED: `BRepAdaptor_Curve2d::GetType() == GeomAbs_Line`
            // (`cxx:480-491`); see the note above.
            None => false,
        };

        let nb = facint.change_line(index).nb_v_point();
        for k in 1..=nb {
            let is_vertex = facint
                .change_line(index)
                .v_point(k)
                .map(|vp| vp.is_vertex(e_si).unwrap_or(false))
                .unwrap_or(false);
            if is_vertex {
                if !isline {
                    let (tol1, tol2) = facint.get_tolerances();
                    let tol = if tol1 > tol2 { tol1 } else { tol2 };
                    let l = facint.change_line(index);
                    let Some(vp) = l.change_v_point(k) else { continue };
                    self.vp_position_vp(vp, &mut vpc, false, tol)?;
                }
                continue;
            }

            let oo_is_vertex = facint
                .change_line(index)
                .v_point(k)
                .map(|vp| vp.is_vertex(oo_e_si).unwrap_or(false))
                .unwrap_or(false);
            if oo_is_vertex {
                if !isline {
                    let (tol1, tol2) = facint.get_tolerances();
                    let tol = if tol1 > tol2 { tol1 } else { tol2 };
                    let l = facint.change_line(index);
                    let Some(vp) = l.change_v_point(k) else { continue };
                    self.vp_position_vp(vp, &mut vpc, false, tol)?;
                }
                continue;
            }

            let p = match facint
                .change_line(index)
                .v_point(k)
                .map(|vp| vp.value().cloned())
            {
                Some(Ok(p)) => p,
                _ => continue,
            };
            let (arc_is_e, arc_is_oo_e) = {
                let l = facint.change_line(index);
                (l.arc_is_edge(e_si), l.arc_is_edge(oo_e_si))
            };

            if arc_is_e {
                // `FUN_projPonL(P, L, (*this), paramC)` (`cxx:522`).
                let param_c = fun_proj_pon_l(&p, facint, index)?;
                let l = facint.change_line(index);
                let Some(vp) = l.change_v_point(k) else { continue };
                match param_c {
                    Some(param_c) => {
                        vp.set_state(crate::abs::State::On, e_si)?;
                        vp.set_edge_on(&earc, param_c, e_si);
                    }
                    None => vp.change_keep(false), // xpu051198
                }
            }

            if arc_is_oo_e {
                let param_c = fun_proj_pon_l(&p, facint, index)?;
                let l = facint.change_line(index);
                let Some(vp) = l.change_v_point(k) else { continue };
                match param_c {
                    Some(param_c) => {
                        vp.set_state(crate::abs::State::On, oo_e_si)?;
                        vp.set_edge_on(&earc, param_c, oo_e_si);
                    }
                    None => vp.change_keep(false), // xpu051198
                }
            }
        }
        Ok(())
    }

    /// `VP_Position(VP, VPC)` (`cxx:558-607`).
    ///
    /// `cxx:593-595` reads the two tolerances off `myFacesIntersector`; the
    /// caller passes their maximum so the line borrow is not held across the
    /// classifier call.
    fn vp_position_vp(
        &mut self,
        vp: &mut VPointInter,
        vpc: &mut VPointInterClassifier,
        assume_in_on: bool,
        tol: f64,
    ) -> Result<(), String> {
        let si = vp.shape_index();
        let (c1, c2) = match si {
            0 => (true, true),
            1 => (false, true),
            2 => (true, false),
            3 => (true, true),
            _ => (false, false),
        };

        if c1 {
            let f1 = self
                .f1
                .clone()
                .ok_or_else(|| "TopOpeBRep_FacesFiller::VP_Position : no F1".to_string())?;
            vpc.v_point_position(
                &f1,
                vp,
                1,
                &mut self.point_classifier,
                assume_in_on,
                tol,
            )?;
        }
        if c2 {
            let f2 = self
                .f2
                .clone()
                .ok_or_else(|| "TopOpeBRep_FacesFiller::VP_Position : no F2".to_string())?;
            vpc.v_point_position(
                &f2,
                vp,
                2,
                &mut self.point_classifier,
                assume_in_on,
                tol,
            )?;
        }
        Ok(())
    }

    /// `PequalVPonR(P3D, VPsi, VP, Lrest)` (`cxx:611-666`).
    fn pequal_vp_on_r(
        &self,
        p3d: &GpPnt,
        vpsi: i32,
        vp: &mut VPointInter,
        lrest: &mut LineInter,
    ) -> Result<bool, String> {
        if lrest.type_line_curve() != TypeLineCurve::Restriction {
            return Ok(false);
        }

        let eol = lrest
            .arc()
            .cloned()
            .ok_or_else(|| "TopOpeBRep_FacesFiller::PequalVPonR : no arc".to_string())?;
        let eol_si = if lrest.arc_is_edge(1) { 1 } else { 2 };

        let nb = lrest.nb_v_point();
        for k in 1..=nb {
            let vpol = match lrest.v_point(k) {
                Some(v) => v,
                None => continue,
            };
            if vpol.is_vertex(1).unwrap_or(false) {
                continue;
            }
            if vpol.shape_index() != vpsi {
                continue;
            }
            if vpol.state(eol_si)? != crate::abs::State::On {
                continue;
            }
            let p3d_ol = *vpol.value()?;
            let tol_e = crate::brep_tool::BRepTool::edge_tolerance(&crate::shape::Edge(
                eol.clone(),
            ));
            if p3d_ol.is_equal_with_tolerance(p3d, tol_e) {
                let param_col = vpol.edge_on_parameter(eol_si)?;
                vp.set_state(crate::abs::State::On, eol_si)?;
                vp.set_edge_on(&eol, param_col, eol_si);
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// `Lminmax(L, pmin, pmax)` (`TopOpeBRep_ProcessGR.cxx:233-255`).
    pub fn lminmax(l: &mut LineInter) -> (f64, f64) {
        let mut pmin = f64::MAX; // RealLast()
        let mut pmax = -f64::MAX; // RealFirst()
        let nb = l.nb_v_point();
        for k in 1..=nb {
            if let Some(p) = l.v_point(k).and_then(|vp| vp.parameter_on_line().ok()) {
                if p < pmin {
                    pmin = p;
                }
                if p > pmax {
                    pmax = p;
                }
            }
        }
        let d = (pmin - pmax).abs();
        let id = d <= PCONFUSION;
        let isper = l.is_periodic();
        let n = l.nb_v_point();
        if id && isper && n >= 2 {
            let per = l.period();
            pmax = pmin + per;
        }
        (pmin, pmax)
    }

    /// `IsVPtransLok(L, iVP, SI12, T)` (`ProcessGR.cxx:320-339`): `false` when
    /// the point is not on the domain of `SI12` or when the line transition is
    /// unknown.
    pub fn is_vp_trans_lok(
        l: &LineInter,
        ivp: i32,
        si12: i32,
    ) -> Result<Option<Transition>, String> {
        let vp = l
            .v_point(ivp)
            .ok_or_else(|| "TopOpeBRep_FacesFiller::IsVPtransLok".to_string())?;
        let is1 = si12 == 1;
        let vp_on_ed = if is1 {
            vp.is_on_dom_s1()?
        } else {
            vp.is_on_dom_s2()?
        };
        if !vp_on_ed {
            return Ok(None);
        }

        let e = vp
            .edge(si12)
            .ok_or_else(|| "TopOpeBRep_FacesFiller::IsVPtransLok".to_string())?;
        let o = e.orientation();
        let t = FFTransitionTool::process_line_transition_with_index(vp, si12, o);
        if t.is_unknown() {
            Ok(None)
        } else {
            Ok(Some(t))
        }
    }

    /// `TransvpOK(L, ivp, SI, isINOUT)` (`ProcessGR.cxx:341-366`).
    pub fn transvp_ok(
        l: &LineInter,
        ivp: i32,
        si: i32,
        is_in_out: bool,
    ) -> Result<bool, String> {
        let Some(t) = Self::is_vp_trans_lok(l, ivp, si)? else {
            return Ok(false);
        };
        let stb = t.before_state();
        let sta = t.after_state();
        use crate::abs::State;
        if is_in_out {
            Ok(stb == State::In && sta == State::Out)
        } else {
            Ok(stb == State::Out && sta == State::In)
        }
    }
}

/// `FUN_projPonL(P, L, FF, paramL)` (`TopOpeBRep_FacesFiller.cxx:45-74`).
///
/// The 3D arm projects on `BRepAdaptor_Curve(E)`; the 2D arm
/// (`BRepAdaptor_Curve2d` + `FUN_tool_projPonC2D`) is UNPORTED — see
/// [`FacesFiller::vp_position_on_r`] — and answers `None`, i.e. the
/// `paramLdef = false` arm of `cxx:56`.
fn fun_proj_pon_l(
    p: &GpPnt,
    facint: &mut FacesIntersector,
    index: i32,
) -> Result<Option<f64>, String> {
    let e = facint
        .change_line(index)
        .arc()
        .cloned()
        .ok_or_else(|| "TopOpeBRep_FacesFiller::FUN_projPonL : no arc".to_string())?;
    let e = crate::shape::Edge(e);
    if crate::brep_tool::BRepTool::edge_curve(&e).is_none() {
        // UNPORTED: `FUN_tool_projPonC2D` (`TopOpeBRepTool_PROJECT.cxx:146-224`).
        return Ok(None);
    }
    let bac = crate::topopebrep::tool::adaptor::CurveAdaptor::of_edge(&e);
    Ok(crate::topopebrep::tool::project_curve::fun_tool_proj_pon_c(p, &bac)
        .map(|(par, _dist)| par))
}
