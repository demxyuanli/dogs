//! Statics of `ChFi3d_Builder_SpKP.cxx` that only touch hatching domains:
//! `CompTra` (`:70-82`), `CompCommonPoint` (`:89-109`), `AdjustParam`
//! (`:161-198`), `Tri` (`:547-622`) and `FillSD` (`:627-744`).
//!
//! Sources: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_SpKP.cxx`
//! (OCCT 8.0.0).

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::elib::clib::in_period;
use occt_core::precision::Precision;
use occt_geom::surface::Surface;

use crate::abs::Orientation;
use crate::chfi_ds::CommonPoint;
use crate::chfi_ds::SurfData;
use crate::geom2d_hatch::Hatcher;
use crate::hatch_gen::{Domain, PointOnElement};
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{edge_vertices, is_same};
use crate::topopebrep::ds::DataStructure;

/// `CompTra(O1, O2, isfirst)` (`cxx:70-82`).
pub(super) fn comp_tra(
    o1: Orientation,
    o2: Orientation,
    is_first: bool,
) -> Orientation {
    let composed = Orientation::compose(o1, o2);
    if is_first {
        composed.reversed()
    } else {
        composed
    }
}

/// `CompCommonPoint(FilPoint, arc, PE, Or)` (`cxx:89-109`).
pub(super) fn comp_common_point(
    fil_point: &mut CommonPoint,
    arc: &Edge,
    pe: &PointOnElement,
    or: Orientation,
) {
    // `TopAbs_Orientation pos = PE.Position();`
    let pos = pe.position();
    let (first, last) = edge_vertices(arc);
    // `if (pos == TopAbs_FORWARD) V = TopExp::FirstVertex(arc); else V =
    // TopExp::LastVertex(arc);` — `CumOri` defaults to `false`, i.e. the
    // vertices as stored on the edge.
    let v = if pos == Orientation::Forward {
        first
    } else {
        last
    };
    if let Some(v) = v {
        fil_point.set_vertex(&v.0);
    }
    fil_point.set_arc(
        Precision::PINTERSECTION,
        arc.shape(),
        pe.parameter(),
        Orientation::compose(arc.orientation(), or),
    );
}

/// `AdjustParam(Dom, f, l, wref, period, pitol)` (`cxx:161-198`).
///
/// Answers `(f, l, acheval)`.
pub(super) fn adjust_param(
    dom: &Domain,
    wref: f64,
    period: f64,
    pitol: f64,
) -> (f64, f64, bool) {
    let mut f = if dom.has_first_point() {
        dom.first_point().parameter()
    } else {
        0.
    };
    let mut l = if dom.has_second_point() {
        dom.second_point().parameter()
    } else {
        period
    };
    if period == 0. {
        return (f, l, false);
    }
    f = in_period(f, wref - pitol, wref + period - pitol);
    l = in_period(l, wref + pitol, wref + period + pitol);
    if l < f {
        f -= period;
        return (f, l, true);
    }
    (f, l, false)
}

/// `Tri(H, iH, Ind, wref, period, pitol, Nbdom)` (`cxx:547-622`).
///
/// `ind` is OCCT's 1-based `Ind` array (`ind[i - 1]`), `nb_dom` its `Nbdom`.
/// Sorts the domain indices by increasing first parameter and folds the domain
/// that has no first point onto the one that has no second point.
pub(super) fn tri(
    h: &mut Hatcher,
    ih: i32,
    ind: &mut [usize],
    nb_dom: &mut usize,
    wref: f64,
    period: f64,
    pitol: f64,
) -> bool {
    for i in 1..=*nb_dom {
        ind[i - 1] = i;
    }

    let mut invert = true;
    while invert {
        invert = false;
        for i in 1..*nb_dom {
            let (f1, _, _) =
                adjust_param(h.domain(ih, ind[i - 1] as i32), wref, period, pitol);
            let (f2, _, _) =
                adjust_param(h.domain(ih, ind[i] as i32), wref, period, pitol);
            if f2 < f1 {
                ind.swap(i - 1, i);
                invert = true;
            }
        }
    }

    let mut i_sans_first = 0usize;
    let mut i_sans_last = 0usize;
    if *nb_dom != 1 {
        for i in 1..=*nb_dom {
            if !h.domain(ih, ind[i - 1] as i32).has_first_point() {
                i_sans_first = i;
            }
            if !h.domain(ih, ind[i - 1] as i32).has_second_point() {
                i_sans_last = i;
            }
        }
    }
    if i_sans_first != 0 {
        if i_sans_last == 0 {
            return false;
        }
        // `HatchGen_Domain* Dom = &H.Domain(iH, Ind(iSansFirst)); HatchGen_
        // PointOnHatching* PH = &H.Domain(iH, Ind(iSansLast)).FirstPoint();
        // double NewPar = H.HatchingCurve(iH).FirstParameter() - period
        //                 + H.Domain(iH, Ind(iSansLast)).FirstPoint().Parameter();
        // PH->SetParameter(NewPar); Dom->SetFirstPoint(*PH);`
        let new_par = h.hatching_curve(ih).first_parameter() - period
            + h.domain(ih, ind[i_sans_last - 1] as i32)
                .first_point()
                .parameter();
        let mut ph = h
            .domain(ih, ind[i_sans_last - 1] as i32)
            .first_point()
            .clone();
        ph.set_parameter(new_par);
        let dom_index = ind[i_sans_first - 1];
        h.hatching_mut(ih)
            .change_domain(dom_index)
            .set_first_point(&ph);

        for k in i_sans_last..*nb_dom {
            ind[k - 1] = ind[k];
        }
        *nb_dom -= 1;
    }
    true
}

/// `FillSD(DStr, CD, M, Dom, ponH, isFirst, ons, pitol, bout)`
/// (`cxx:627-744`).
///
/// `m` is OCCT's `NCollection_DataMap<int, handle<Adaptor2d_Curve2d>>`: the
/// element index to the arc that produced it. Only the `BRepAdaptor_Curve2d`
/// down-cast branch of `FillSD` uses it, and that branch needs `Edge()` alone,
/// so an entry is `None` for a `Geom2dAdaptor_Curve` element.
#[allow(clippy::too_many_arguments)]
pub(super) fn fill_sd(
    ds: &mut DataStructure,
    cd: &mut SurfData,
    m: &HashMap<i32, Option<Edge>>,
    dom: &Domain,
    pon_h: f64,
    is_first: bool,
    ons: i32,
    pitol: f64,
    bout: &TopoShape,
) -> Result<(), String> {
    let opp = 3 - ons;
    let surf: Arc<dyn Surface> = ds
        .surface(cd.surf())
        .surface()
        .cloned()
        .ok_or_else(|| "ChFi3d_Builder::SplitKPart: null surface".to_string())?;

    let p_p_h = if is_first && dom.has_first_point() {
        Some(dom.first_point().clone())
    } else if !is_first && dom.has_second_point() {
        Some(dom.second_point().clone())
    } else {
        None
    };

    let Some(p_p_h) = p_p_h else {
        // `CD->ChangeInterference(ons).SetParameter(ponH, isFirst);
        //  pcons->Value(ponH).Coord(x, y);
        //  CD->ChangeVertex(isFirst, ons).SetPoint(Surf->Value(x, y));`
        cd.interference_mut(ons).set_parameter(pon_h, is_first);
        let pcons = cd
            .interference(ons)
            .p_curve_on_surf()
            .cloned()
            .ok_or_else(|| "ChFi3d_Builder::SplitKPart: null pcurve".to_string())?;
        let p = pcons.value(pon_h);
        cd.vertex_mut(is_first, ons).set_point(&surf.value(p.x(), p.y()));
        // Trailing `if (!Popp.IsOnArc())` block.
        if !cd.vertex(is_first, opp).is_on_arc() {
            cd.interference_mut(opp).set_parameter(pon_h, is_first);
            let pcopp = cd
                .interference(opp)
                .p_curve_on_surf()
                .cloned()
                .ok_or_else(|| "ChFi3d_Builder::SplitKPart: null pcurve".to_string())?;
            let q = pcopp.value(pon_h);
            cd.vertex_mut(is_first, opp)
                .set_point(&surf.value(q.x(), q.y()));
        }
        return Ok(());
    };

    // Modification to find already existing vertexes.
    let mut le_type = 1usize;
    let nb_int = p_p_h.nb_points();
    if nb_int > 1 {
        let mut v1: Option<TopoShape> = None;
        let mut v2: Option<TopoShape> = None;
        let mut suite = true;
        // `for (; trouve;)` (`cxx:672`): every exit of OCCT's body is a
        // `break`, so the flag becomes a plain loop here.
        loop {
            let pe = p_p_h.point(le_type);
            let ie = pe.index();
            match m.get(&ie).and_then(|e| e.as_ref()) {
                Some(he) => {
                    let (a, b) = edge_vertices(he);
                    v1 = a.map(|v| v.0);
                    v2 = b.map(|v| v.0);
                }
                None => suite = false,
            }
            let same = |v: &Option<TopoShape>| {
                v.as_ref().is_some_and(|v| is_same(v, bout))
            };
            if (same(&v1) || same(&v2)) && suite {
                break;
            }
            suite = true;
            le_type += 1;
            if le_type > nb_int {
                le_type = 1;
                break;
            }
        }
    }

    let pe = p_p_h.point(le_type);
    let ie = pe.index();
    let Some(Some(he)) = m.get(&ie) else {
        return Ok(());
    };
    let e = he;

    if pe.position() != Orientation::Internal {
        let mut o = cd.interference(ons).transition();
        if is_first {
            o = o.reversed();
        }
        comp_common_point(cd.vertex_mut(is_first, ons), e, pe, o);
    } else {
        let o = comp_tra(
            cd.interference(ons).transition(),
            e.orientation(),
            is_first,
        );
        cd.vertex_mut(is_first, ons)
            .set_arc(pitol, e.shape(), pe.parameter(), o);
    }

    cd.interference_mut(ons).set_parameter(pon_h, is_first);
    let pcadj = cd
        .interference(ons)
        .p_curve_on_surf()
        .cloned()
        .ok_or_else(|| "ChFi3d_Builder::SplitKPart: null pcurve".to_string())?;
    let p = pcadj.value(pon_h);
    cd.vertex_mut(is_first, ons)
        .set_point(&surf.value(p.x(), p.y()));

    if !cd.vertex(is_first, opp).is_on_arc() {
        cd.interference_mut(opp).set_parameter(pon_h, is_first);
        let pcopp = cd
            .interference(opp)
            .p_curve_on_surf()
            .cloned()
            .ok_or_else(|| "ChFi3d_Builder::SplitKPart: null pcurve".to_string())?;
        let q = pcopp.value(pon_h);
        cd.vertex_mut(is_first, opp)
            .set_point(&surf.value(q.x(), q.y()));
    }
    Ok(())
}