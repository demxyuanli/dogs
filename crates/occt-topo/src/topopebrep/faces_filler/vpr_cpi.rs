//! `TopOpeBRep_vpr.cxx` — the curve/point interference of a VPoint on a
//! restriction, filtered by `FUNBREP_topokpart`.
//!
//! Source: `TopOpeBRep_vpr.cxx:782-853` (`FUN_processCPI`).

use crate::topopebrep::ds::export::fds_parameter;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::interference::InterferenceHandle;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::ds::process_interferences_tool::make_cpv_interference;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::faces_filler::FacesFiller;
use crate::topopebrep::faces_intersector::FacesIntersector;
use crate::topopebrep::interf::v_point_inter::VPointInter;
use crate::topopebrep::kpart::fun_brep_topokpart;

/// `FUN_processCPI(FF, VP, F, ShapeIndex, L, pDS, transLine, DSCIL, Ifound,
/// Gfound, PVKind, PVIndex, keptVPnbr)` (`vpr.cxx:782-853`).
///
/// `keptVPnbr` is an `int&` into `FF.mykeptVPnbr`; the caller writes it back.
#[allow(clippy::too_many_arguments)]
pub fn fun_process_cpi(
    ff: &mut FacesFiller,
    facint: &mut FacesIntersector,
    hds: &mut HDataStructure,
    vp: &VPointInter,
    f: &crate::shape::TopoShape,
    shape_index: i32,
    trans_line: &Transition,
    ifound: Option<&InterferenceHandle>,
    gfound: bool,
    pv_kind: Kind,
    pv_index: i32,
    kept_vp_nbr: &mut i32,
) -> Result<(), String> {
    let oo_shape_index = if shape_index == 1 { 2 } else { 1 };
    let mut t_trans_line = trans_line.clone(); // `cxx:799`
    // prequesitory : current line is not on edge.
    let mut parline = vp.parameter_on_line()?;
    let si_is_vertex = vp.is_vertex(shape_index)?;
    let oo_is_vertex = vp.is_vertex(oo_shape_index)?;
    let e = vp
        .edge(shape_index)
        .ok_or_else(|| "FUN_processCPI : no edge".to_string())?
        .clone();

    // xpu010299 : we do not keep interferences with same parameters on curve
    //             PRO16120(f3,f4 -> null c1)
    if let Some(last) = ff.dscil.last() {
        let par = {
            let i = last.read().expect("interference lock");
            fds_parameter(&i).unwrap_or(0.0)
        };
        let dd = (par - parline).abs(); // en fait, ce sont des entiers
        if dd == 0.0 {
            return Ok(());
        }
    }

    // dist(p2d1,p2d2) < toluv => p2d1, p2d2 are considered equal.
    // NYI : compute uvtol with the original faces. By default, we set
    // toluv = TolClass
    let toluv = 1.0e-8;

    let dscil = ff.dscil.clone();
    let l = facint.change_line(ff.line_index);
    let keep = fun_brep_topokpart(
        ifound,
        &dscil,
        l,
        vp,
        hds.ds(),
        &e,
        f,
        toluv,
        &mut parline,
        &mut t_trans_line,
    )?;

    let mut keep = keep;
    if keep {
        *kept_vp_nbr += 1;
        if *kept_vp_nbr > 2 {
            keep = false;
        }
    }
    if !keep {
        return Ok(());
    }

    let gk_cpv = if gfound {
        pv_kind
    } else if si_is_vertex || oo_is_vertex {
        Kind::Vertex
    } else {
        Kind::Point
    };

    let cpi = make_cpv_interference(&t_trans_line, 0, pv_index, parline, gk_cpv);
    ff.store_curve_interference(hds, &cpi);
    Ok(())
}
