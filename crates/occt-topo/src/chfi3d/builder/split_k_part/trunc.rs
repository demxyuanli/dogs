//! `ChFi3d_Builder::Trunc` (`ChFi3d_Builder_SpKP.cxx:372-506`).
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_SpKP.cxx`
//! (OCCT 8.0.0).

use std::sync::Arc;

use occt_geom::surface::Surface;

use crate::chfi_ds::SurfData;
use crate::chfi_ds::SpineHandle;
use crate::topo_tools::edge_vertices;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::super::ChFi3dBuilder;
use super::spine_params::{
    elementary_spine, k_particular, param_on_spine, spine_edge, ya_un_voisin,
};

impl ChFi3dBuilder {
    /// `Trunc(SD, Spine, S1, S2, iedge, isfirst, cntlFiOnS)` (`cxx:372-506`).
    ///
    /// Cuts `SD` at the beginning (`isfirst`) or the end of the elementary
    /// spine `iedge` when the corner point falls outside the edge domain, and
    /// when the deviation is large enough to be a genuine truncation.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn trunc(
        &mut self,
        sd: &mut SurfData,
        spine: &SpineHandle,
        s1: &SurfaceAdaptor,
        s2: &SurfaceAdaptor,
        iedge: i32,
        is_first: bool,
        cntl_fi_on_s: i32,
    ) -> Result<(), String> {
        // Return points and tangents on edge and spine.
        let wtg = sd.interference_on_s1().parameter(is_first);
        let (wsp, _bid) = param_on_spine(
            self.my_ds.change_ds(),
            wtg,
            sd,
            spine,
            iedge,
            false,
            false,
            self.tolesp,
        )?;

        // Modif against Vertex isolated on spine.
        let support = spine_edge(spine, iedge)?;
        let (b1, b2) = edge_vertices(&support);
        let mut bout1 = b1.map(|v| v.0);
        let mut bout2 = b2.map(|v| v.0);
        if support.orientation().is_reversed() {
            std::mem::swap(&mut bout1, &mut bout2);
        }
        if !is_first {
            bout1 = bout2;
        }
        // finmodif

        let bc = elementary_spine(spine, iedge)?;
        let edf = bc.first_parameter();
        let edl = bc.last_parameter();
        let edglen = edl - edf;

        let (ped, ded) = if support.orientation().is_forward() {
            bc.d1(wtg + edf)
                .ok_or_else(|| "ChFi3d_Builder::Trunc: no tangent on the edge".to_string())?
        } else {
            let (p, d) = bc
                .d1(-wtg + edl)
                .ok_or_else(|| "ChFi3d_Builder::Trunc: no tangent on the edge".to_string())?;
            (p, d.reversed())
        };
        let (psp, dsp) = spine
            .write()
            .map_err(|_| "ChFi3d_Builder::Trunc: poisoned spine".to_string())?
            .d1(wsp)
            .ok_or_else(|| "ChFi3d_Builder::Trunc: no tangent on the spine".to_string())?;

        let surf: Arc<dyn Surface> = self
            .my_ds
            .change_ds()
            .surface(sd.surf())
            .surface()
            .cloned()
            .ok_or_else(|| "ChFi3d_Builder::Trunc: null surface".to_string())?;
        let pp1 = sd
            .interference_on_s1()
            .p_curve_on_surf()
            .cloned()
            .ok_or_else(|| "ChFi3d_Builder::Trunc: null pcurve".to_string())?
            .value(wtg);
        let pp2 = sd
            .interference_on_s2()
            .p_curve_on_surf()
            .cloned()
            .ok_or_else(|| "ChFi3d_Builder::Trunc: null pcurve".to_string())?
            .value(wtg);
        let p1 = surf.value(pp1.x(), pp1.y());
        let p2 = surf.value(pp2.x(), pp2.y());

        let mut tron = false;
        let ang = dsp.angle(&ded);
        let dis1 = psp.distance(&ped);
        let dis2 = p1.distance(&p2);
        if ang > std::f64::consts::PI / 18. {
            tron = true;
        }
        if dis1 >= 0.1 * dis2 {
            tron = true;
        }

        if !tron {
            if let Some(ivois) = ya_un_voisin(spine, iedge, is_first)? {
                if let (Some(f1), Some(f2)) = (s1.face(), s2.face()) {
                    let mut fbid = None;
                    let cp1 = sd.vertex(is_first, 1);
                    let cp2 = sd.vertex(is_first, 2);
                    let found1 = cp1.is_on_arc()
                        && self.search_face(Some(spine), cp1, Some(f1), &mut fbid);
                    let found2 = cp2.is_on_arc()
                        && self.search_face(Some(spine), cp2, Some(f2), &mut fbid);
                    if !(found1 || found2) {
                        tron = k_particular(spine, ivois, s1, s2);
                    }
                }
            }
        }

        // modification of lvt against isolated vertex
        if !tron && ya_un_voisin(spine, iedge, is_first)?.is_some() {
            let mut nbed: i32 = -2;
            if let Some(bout1) = bout1.as_ref() {
                for _ in self.my_ve_map.get_from_key(bout1) {
                    nbed += 1;
                }
            }
            if nbed < 3 {
                tron = true;
            }
        }
        // finmodif

        if tron {
            let mut par = 0.;
            let mut d_par = 0.;
            if !is_first {
                par = edglen;
            }
            if cntl_fi_on_s != 0 {
                // Detect the case where FaceInterference ends before the place
                // we are going to truncate SD (eap, occ354).
                let fi_par = sd.interference(cntl_fi_on_s).parameter(!is_first);
                let is_the_case = if is_first { par > fi_par } else { par < fi_par };
                if is_the_case {
                    d_par = par - fi_par;
                    par = fi_par;
                }
            }
            for i in 1..=2 {
                sd.interference_mut(i).set_parameter(par, is_first);
                let pc = sd
                    .interference(i)
                    .p_curve_on_surf()
                    .cloned()
                    .ok_or_else(|| "ChFi3d_Builder::Trunc: null pcurve".to_string())?;
                let p = pc.value(par);
                sd.vertex_mut(is_first, i).reset();
                sd.vertex_mut(is_first, i)
                    .set_point(&surf.value(p.x(), p.y()));
                if is_first {
                    let v = spine
                        .read()
                        .map_err(|_| "ChFi3d_Builder::Trunc: poisoned spine".to_string())?
                        .first_parameter_of(iedge)
                        .ok_or_else(|| {
                            "ChFi3d_Builder::Trunc: FirstParameter".to_string()
                        })?;
                    sd.set_first_spine_param(v - d_par);
                } else {
                    let v = spine
                        .read()
                        .map_err(|_| "ChFi3d_Builder::Trunc: poisoned spine".to_string())?
                        .last_parameter_of(iedge)
                        .ok_or_else(|| {
                            "ChFi3d_Builder::Trunc: LastParameter".to_string()
                        })?;
                    sd.set_last_spine_param(v - d_par);
                }
            }
        }
        Ok(())
    }
}
