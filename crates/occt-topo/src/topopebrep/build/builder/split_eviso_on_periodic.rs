//! `TopOpeBRepBuild_Builder::SplitEvisoONperiodicF`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_GridSS.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `FUN_getG(P, LI, HDS, iEinterf)` | `89-116` |
//! | `FUN_EPIforEvisoONperiodicF(E, F, EPIlist, HDS, loCPI)` | `138-190` |
//! | `FUN_getEPI(LI, EPI)` | `237-258` |
//! | `FUN_getEPIonEds(FOR, HDS, EPI)` | `261-279` |
//! | `FUN_SplitEvisoONperiodicF(HDS, FF)` | `1184-1294` |
//! | `SplitEvisoONperiodicF()` | `1314-1346` |
//!
//! The pass exists for two k-part situations (`cxx:1300-1311`): cylinders
//! same-domain on their cylindrical face with closing edges that are not
//! same-domain, and a cylinder + sphere interfering on the circular edge of
//! the cylinder. In both, the p-curve of the section edge on the closed face is
//! a `Viso` line whose `u` is not bounded in `[0, 2Pi]`; the pass adds a
//! `TopOpeBRepDS_CurvePointInterference` that splits the edge at `u = 2Pi`.
//!
//! `FUN_getEPIonEds(FS, HDS, EPIlist)` is called with the *face* in
//! `FUN_SplitEvisoONperiodicF` (`cxx:1286`) even though its parameter is named
//! `FOR`: the explorer walks the edges of what it is given, so the same call
//! collects the point/face interferences of `FS`'s edges.

use std::sync::Arc;

use occt_core::gp::GpPnt;

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::explorer_shapes;
use crate::topopebrep::build::tools::fc2d_curve_on_surface;
use crate::topopebrep::ds::export::fds_data;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference_tool::make_edge_interference;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::{
    fun_select_gk_interference, fun_select_sk_interference, fun_select_trasha_interference,
};
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::ds::{InterferenceHandle, InterferenceList};
use crate::topopebrep::tool::corr_iso::{CorrIso, SPLITEDGE};
use crate::topopebrep::tool::geometry::fun_tool_is_closed_s;
use crate::topopebrep::tool::project_curve::fun_tool_value;
use crate::topopebrep::tool::topology::fun_tool_pcurveon_f;
use crate::topopebrep::tool::uv_iso_curve;

use super::Builder;

/// `FUN_getG(P, LI, HDS, iEinterf)` (`cxx:89-116`): the geometry index of the
/// first `CurvePointInterference` of `LI` whose data-structure point equals `P`
/// within that point's own tolerance, together with its support index.
///
/// OCCT writes `iEinterf` on every visited interference and reads it only when
/// a match was found (`cxx:169-176`), so the two answers travel together here.
fn fun_get_g(
    p: &GpPnt,
    li: &[InterferenceHandle],
    hds: &HDataStructure,
) -> Option<(i32, i32)> {
    for i in li {
        let interference = i.read().expect("interference lock");
        if !interference.is_curve_point() {
            continue;
        }
        let gi = interference.geometry();
        let i_einterf = interference.support();
        let Some(dsp) = hds.point(gi) else {
            continue;
        };
        if dsp.point().is_equal_with_tolerance(p, dsp.tolerance()) {
            return Some((gi, i_einterf));
        }
    }
    None
}

/// `FUN_getEPI(LI, EPI)` (`cxx:237-258`): keep the `CurvePointInterference`s
/// whose geometry is a POINT and whose support is a FACE.
fn fun_get_epi(li: &[InterferenceHandle], epi: &mut InterferenceList) {
    for i in li {
        let (is_curve_point, gt, st) = {
            let interference = i.read().expect("interference lock");
            (
                interference.is_curve_point(),
                interference.geometry_type(),
                interference.support_type(),
            )
        };
        if !is_curve_point {
            continue;
        }
        if gt != Kind::Point || st != Kind::Face {
            continue;
        }
        epi.push(Arc::clone(i));
    }
}

/// `FUN_getEPIonEds(FOR, HDS, EPI)` (`cxx:261-279`): every edge of `for_` known
/// to the data structure hands its interferences to `fun_get_epi`.
fn fun_get_epi_on_eds(for_: &TopoShape, hds: &HDataStructure, epi: &mut InterferenceList) {
    for e in explorer_shapes(for_, ShapeType::Edge) {
        if !hds.has_shape(&e, true) {
            continue;
        }
        let lii = hds.ds().shape_interferences(&e, true).to_vec();
        fun_get_epi(&lii, epi);
    }
}

/// `FUN_EPIforEvisoONperiodicF(E, F, EPIlist, HDS, loCPI)`
/// (`cxx:138-190`): the `CurvePointInterference` splitting `E` at the point
/// whose p-curve parameter on `F` is `u = 2Pi`, if `E`'s p-curve leaves the
/// `[0, 2Pi]` window at all.
fn fun_epi_for_eviso_on_periodic_f(
    e: &TopoShape,
    f: &TopoShape,
    epi_list: &[InterferenceHandle],
    hds: &HDataStructure,
) -> Option<InterferenceHandle> {
    let face = Face::wrap(f.clone())?;

    // `cxx:146-148`: `CORRISO(F); CORRISO.Init(F);` then `Refclosed(1, uper)`.
    let mut corr_iso = CorrIso::new(f);
    corr_iso.init(f);
    let on_u = corr_iso.ref_closed(1).is_some();
    let tol_f = BRepTool::face_tolerance(&face);
    let tol_u = corr_iso.tol(1, tol_f);
    let (recadre, parone) = corr_iso.edge_out_of_bounds_uv(e, on_u, tol_u);

    // `cxx:155-158`
    if recadre != SPLITEDGE {
        return None;
    }

    // `cxx:159-164`
    let p3d = fun_tool_value(parone, e)?;

    // `cxx:166-176`
    let (i_g, i_einterf) = fun_get_g(&p3d, epi_list, hds)?;
    let i_einterf = match hds.shape(i_einterf, true) {
        Some(s) if s.shape_type() == ShapeType::Edge => i_einterf,
        _ => 0,
    };

    // `cxx:179-188`
    let i_s = hds.shape_index(e, true);
    let mut t = Transition::new(State::In, State::In, ShapeType::Edge, ShapeType::Edge);
    t.set_index(i_s);
    Some(make_edge_interference(
        &t,
        Kind::Edge,
        i_einterf,
        Kind::Point,
        i_g,
        parone,
    ))
}

/// `FUN_SplitEvisoONperiodicF(HDS, FF)` (`cxx:1184-1294`). Always answers
/// `true` in OCCT: every early exit returns `true`, and the only failure path
/// is the `Standard_ProgramError` thrown at `cxx:1259`.
fn fun_split_eviso_on_periodic_f(hds: &mut HDataStructure, ff: &TopoShape) -> bool {
    // `cxx:1186-1192`
    let lli = hds.ds().shape_interferences(ff, true).to_vec();
    if lli.is_empty() {
        return true;
    }

    // `cxx:1193-1198`: `LI` is a copy of the face's interference list, the
    // three selectors below *move* entries out of it.
    let mut li: InterferenceList = lli.clone();

    // `cxx:1200-1212`: `LI1` = FACE-transition interferences on an EDGE
    // geometry, `LI2` = those whose support is a FACE.
    let mut li1: InterferenceList = Vec::new();
    if fun_select_gk_interference(&mut li, Kind::Edge, &mut li1) < 1 {
        return true;
    }
    let mut li2: InterferenceList = Vec::new();
    if fun_select_sk_interference(&mut li1, Kind::Face, &mut li2) < 1 {
        return true;
    }
    let mut li3: InterferenceList = Vec::new();
    if fun_select_trasha_interference(&mut li2, ShapeType::Face, &mut li3) < 1 {
        return true;
    }

    // `cxx:1221-1291`
    for ssi in li3 {
        let (_, gi, _, si) = {
            let interference = ssi.read().expect("interference lock");
            fds_data(&interference)
        };

        let Some(fs) = hds.shape(si, true) else {
            continue;
        };
        let Some(fs_face) = Face::wrap(fs.clone()) else {
            continue;
        };
        // `cxx:1232-1237`
        if !fun_tool_is_closed_s(&fs_face) {
            continue;
        }

        let Some(eg) = hds.shape(gi, true) else {
            continue;
        };
        let Some(eg_edge) = Edge::wrap(eg.clone()) else {
            continue;
        };
        // `cxx:1241-1246`
        if !hds.ds().is_section_edge(&eg, true) {
            continue;
        }

        // `cxx:1251-1266`: the p-curve of `EG` on `FS`, projected on demand.
        let pc = match fc2d_curve_on_surface(&eg_edge, &fs_face) {
            Some((pc, _pf, _pl)) => pc,
            None => {
                if !fun_tool_pcurveon_f(&fs, &eg) {
                    panic!("_Builder::SplitONVisolineonCyl");
                }
                match fc2d_curve_on_surface(&eg_edge, &fs_face) {
                    Some((pc, _pf, _pl)) => pc,
                    // OCCT keeps going with a null handle and dereferences it
                    // in `UVISO` (`cxx:1268`); the port declines the edge.
                    None => continue,
                }
            }
        };

        // `cxx:1268-1273`
        let Some(iso) = uv_iso_curve(pc.as_ref()) else {
            continue;
        };
        if !iso.iso_v {
            continue;
        }

        // `cxx:1286-1291`
        let mut epi_list: InterferenceList = Vec::new();
        fun_get_epi_on_eds(&fs, hds, &mut epi_list);
        let Some(cpi) = fun_epi_for_eviso_on_periodic_f(&eg, &fs, &epi_list, hds) else {
            continue;
        };
        hds.change_ds()
            .change_shape_interferences(&eg)
            .push(cpi);
    }

    true
}

impl Builder {
    /// `SplitEvisoONperiodicF()` (`cxx:1314-1346`): run the pass on every
    /// periodic face of the data structure, in index order.
    pub fn split_eviso_on_periodic_f(&mut self) {
        let Some(hds) = self.my_data_structure.as_mut() else {
            return;
        };
        let sha = hds.nb_shapes();
        for i in 1..=sha {
            let Some(for_) = hds.shape(i, true) else {
                continue;
            };
            if for_.shape_type() != ShapeType::Face {
                continue;
            }
            // `cxx:1328-1334`: `S->IsUPeriodic() || S->IsVPeriodic()`.
            let Some(face) = Face::wrap(for_.clone()) else {
                continue;
            };
            let periodic = match BRepTool::face_surface_world(&face) {
                Some(s) => s.is_u_periodic() || s.is_v_periodic(),
                None => false,
            };
            if !periodic {
                continue;
            }

            // `cxx:1336-1342`
            let mut ff = for_.clone();
            ff.set_orientation(Orientation::Forward);
            if !fun_split_eviso_on_periodic_f(hds, &ff) {
                panic!("_Builder::SplitONVisolineonCyl");
            }
        }
    }
}
