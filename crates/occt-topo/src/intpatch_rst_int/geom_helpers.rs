//! File-statics of `IntPatch_RstInt.cxx` (`:48-429`).
//!
//! Ported here, in source order:
//! * `Recadre` (`cxx:49-129`).
//! * the two `Tol3d` inlines (`cxx:133-145`).
//! * `CoincideOnArc` (`cxx:147-176`).
//! * `VerifyTgline` (`cxx:178-218`).
//! * `GetLinePoint2d` (`cxx:220-273`).
//! * `FindParameter` (`cxx:275-421`).
//!
//! `ArePnt2dEqual` (`cxx:423-429`) is UNPORTED: the function is defined and
//! never called in `IntPatch_RstInt.cxx`.
//!
//! Translation notes:
//! * `occ::handle<Adaptor3d_HVertex>` / `<Adaptor2d_Curve2d>` are the port's
//!   [`HVertex`] / [`Curve2dAdaptor`]; `Domain->Tol3d` is
//!   [`TopolTool::tol3d_of_curve`] / [`TopolTool::tol3d_of_vertex`].
//! * `IntPatch_HInterTool::Project` answers `Option` in the port; OCCT's call
//!   site (`cxx:311`, `:316`) ignores its `bool`, which leaves `Param` /
//!   `p2d` unwritten on failure. The port maps that failure to "not found".
//! * `RealLast()` (`cxx:152`, `:174`, `:871`) is `f64::MAX`.

use occt_core::gp::{GpPnt, GpPnt2d, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::Surface;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_surface::SurfaceKind;
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::geom_int::{GeomIntLine, IntPatchIType};
use crate::int_tools_wline::WLine;
use crate::intpatch_h_inter_tool as hinter;

/// `Recadre(typeS1, typeS2, wlin, Param, U1, V1, U2, V2)` (`cxx:49-129`).
///
/// Shifts the UV of a point onto the period of the walking-line sample it is
/// attached to. `Param` is clamped to the line's point range (`cxx:57-65`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn recadre(
    type_s1: SurfaceKind,
    type_s2: SurfaceKind,
    wlin: &WLine,
    param: i32,
    u1: &mut f64,
    v1: &mut f64,
    u2: &mut f64,
    v2: &mut f64,
) {
    let mut param = param;
    let nbpnts = wlin.nb_pnts();
    if param < 1 {
        param = 1;
    } else if param > nbpnts {
        param = nbpnts;
    }
    let (u1p, v1p, u2p, v2p) = wlin.point(param).parameters();
    let two_pi = std::f64::consts::PI + std::f64::consts::PI;

    match type_s1 {
        // `cxx:80-84`: `case GeomAbs_Cylinder/Cone/Sphere/Torus:` then `break`.
        SurfaceKind::Cylinder | SurfaceKind::Cone | SurfaceKind::Sphere | SurfaceKind::Torus => {
            while *u1 < (u1p - 1.5 * std::f64::consts::PI) {
                *u1 += two_pi;
            }
            while *u1 > (u1p + 1.5 * std::f64::consts::PI) {
                *u1 -= two_pi;
            }
        }
        _ => {}
    }
    if type_s1 == SurfaceKind::Torus {
        while *v1 < (v1p - 1.5 * std::f64::consts::PI) {
            *v1 += two_pi;
        }
        while *v1 > (v1p + 1.5 * std::f64::consts::PI) {
            *v1 -= two_pi;
        }
    }

    match type_s2 {
        SurfaceKind::Cylinder | SurfaceKind::Cone | SurfaceKind::Sphere | SurfaceKind::Torus => {
            while *u2 < (u2p - 1.5 * std::f64::consts::PI) {
                *u2 += two_pi;
            }
            while *u2 > (u2p + 1.5 * std::f64::consts::PI) {
                *u2 -= two_pi;
            }
        }
        _ => {}
    }
    // `cxx:118-126`: the V test of the second surface uses `V1p`, not `V2p`.
    // That is what the source says; it is kept verbatim.
    if type_s2 == SurfaceKind::Torus {
        while *v2 < (v1p - 1.5 * std::f64::consts::PI) {
            *v2 += two_pi;
        }
        while *v2 > (v2p + 1.5 * std::f64::consts::PI) {
            *v2 -= two_pi;
        }
    }
}

/// `Tol3d(arc, Domain, tolDef = 0.)` (`cxx:140-145`).
pub(crate) fn tol_3d_arc(arc: &Curve2dAdaptor, domain: &TopolTool, tol_def: f64) -> f64 {
    if domain.has3d() {
        TopolTool::tol3d_of_curve(arc)
    } else if tol_def < CONFUSION {
        CONFUSION
    } else {
        tol_def
    }
}

/// `Tol3d(vtx, Domain, tolDef = 0.)` (`cxx:133-138`).
pub(crate) fn tol_3d_vertex(vtx: &HVertex, domain: &TopolTool, tol_def: f64) -> f64 {
    if domain.has3d() {
        TopolTool::tol3d_of_vertex(vtx)
    } else if tol_def < CONFUSION {
        CONFUSION
    } else {
        tol_def
    }
}

/// `CoincideOnArc(Ptsommet, A, Surf, Toler, Domain, Vtx)` (`cxx:147-176`).
///
/// Returns whether a domain vertex of `arc` coincides with `ptsommet` and, if
/// so, that vertex.
pub(crate) fn coincide_on_arc(
    ptsommet: &GpPnt,
    arc: &Curve2dAdaptor,
    surf: &dyn Surface,
    toler: f64,
    domain: &mut TopolTool,
) -> (bool, Option<HVertex>) {
    let mut distmin = f64::MAX;
    let tolarc = toler.max(tol_3d_arc(arc, domain, 0.0));
    let mut best = None;

    domain.initialize_curve(arc);
    domain.init_vertex_iterator();
    while domain.more_vertex() {
        let vtx1 = match domain.vertex() {
            Some(v) => v,
            None => {
                domain.next_vertex();
                continue;
            }
        };
        let prm = hinter::parameter(&vtx1, arc);
        let p2d = match arc.d0(prm) {
            Some(p) => p,
            None => {
                domain.next_vertex();
                continue;
            }
        };
        let point = surf.d0(p2d.x(), p2d.y());
        let dist = point.distance(ptsommet);
        let tol = tolarc.max(tol_3d_vertex(&vtx1, domain, 0.0));

        if dist <= tol && dist <= distmin {
            // the best coincidence
            distmin = dist;
            best = Some(vtx1);
        }
        domain.next_vertex();
    }
    (distmin < f64::MAX, best)
}

/// `VerifyTgline(wlin, param, Tol, Tgl)` (`cxx:178-218`).
///
/// When the tangent is degenerate, walks the line forwards / backwards to the
/// first sample far enough away and uses that chord instead.
pub(crate) fn verify_tgline(wlin: &WLine, param: i32, tol: f64, tgl: &mut GpVec) {
    if tgl.x().abs() < tol && tgl.y().abs() < tol && tgl.z().abs() < tol {
        // Build a longer tangent, to avoid points that are very close
        // together inside the walking step (`cxx:186-188`).
        let nbpt = wlin.nb_pnts();
        let mut forward = (nbpt - param) >= (param - 1);
        for _n in (1..=2).rev() {
            if forward {
                for i in (param + 1)..=nbpt {
                    let t = GpVec::from_pnts(&wlin.point(param).value(), &wlin.point(i).value());
                    if t.x().abs() >= tol || t.y().abs() >= tol || t.z().abs() >= tol {
                        *tgl = t;
                        return;
                    }
                }
            } else {
                for i in (1..=param - 1).rev() {
                    let t = GpVec::from_pnts(&wlin.point(i).value(), &wlin.point(param).value());
                    if t.x().abs() >= tol || t.y().abs() >= tol || t.z().abs() >= tol {
                        *tgl = t;
                        return;
                    }
                }
            }
            forward = !forward;
        }
    }
}

/// `GetLinePoint2d(L, param, OnFirst, U, V)` (`cxx:220-273`).
///
/// Linear 2d interpolation of the line's samples at the (possibly fractional)
/// vertex parameter `param`.
pub(crate) fn get_line_point2d(
    line: &GeomIntLine,
    param: f64,
    on_first: bool,
) -> Result<(f64, f64), String> {
    let typ_l = line.arc_type();
    let nb_pnts = match typ_l {
        IntPatchIType::Walking => line.as_wline().map(|w| w.nb_pnts()),
        IntPatchIType::Restriction => line.as_rline().map(|r| r.nb_pnts()),
        _ => None,
    }
    .ok_or_else(|| "Standard_DomainError".to_string())?;

    let whole = param.trunc();
    let mut irang = whole as i32;
    let par = if irang == nb_pnts {
        irang -= 1;
        1.0
    } else {
        (param - whole).abs()
    };

    let params = |i: i32| -> (f64, f64) {
        let p = match typ_l {
            IntPatchIType::Walking => line.as_wline().map(|w| *w.point(i)),
            _ => line.as_rline().map(|r| *r.point(i)),
        };
        let p = p.expect("index checked against NbPnts");
        if on_first {
            p.parameters_on_s1()
        } else {
            p.parameters_on_s2()
        }
    };
    let (us1, vs1) = params(irang);
    let (us2, vs2) = params(irang + 1);
    Ok(((1.0 - par) * us1 + par * us2, (1.0 - par) * vs1 + par * vs2))
}

/// `FindParameter(L, OtherSurf, Tol, Ptsom, Ptsom2d, Param, Tgl, ParamApproche,
/// OnFirst)` (`cxx:275-421`).
///
/// Returns the parameter on the line of the point closest to `ptsom` and the
/// tangent there.
#[allow(clippy::too_many_arguments)]
pub(crate) fn find_parameter(
    line: &GeomIntLine,
    other_surf: &dyn Surface,
    tol: f64,
    ptsom: &GpPnt,
    ptsom2d: &GpPnt2d,
    param_approche: i32,
    on_first: bool,
) -> Result<Option<(f64, GpVec)>, String> {
    let tol2 = tol * tol;
    let typ_l = line.arc_type();
    let mut tgl = GpVec::zero();

    if typ_l == IntPatchIType::Restriction {
        let rlin = line
            .as_rline()
            .ok_or_else(|| "Standard_DomainError".to_string())?;
        let use_arc = (!on_first && rlin.arc_on_s1) || (on_first && rlin.arc_on_s2);
        if !use_arc {
            return Ok(None);
        }
        let c = match rlin.c2d.as_ref() {
            Some(c) => c.as_ref(),
            None => return Ok(None),
        };
        let projected = hinter::project(c, ptsom2d);
        let (param, p2d) = match projected {
            Some(v) => v,
            None => return Ok(None),
        };
        let (_, d2d) = c.d1(param);
        let (ptbid, d1u, d1v) = other_surf.d1(p2d.x(), p2d.y());
        if ptbid.square_distance(ptsom) > tol2 {
            return Ok(None);
        }
        tgl = d1u
            .multiplied_scalar(d2d.x())
            .added(&d1v.multiplied_scalar(d2d.y()));
        return Ok(Some((param, tgl)));
    }

    if typ_l != IntPatchIType::Walking {
        return Err("Standard_DomainError".to_string());
    }

    let wlin = line
        .as_wline()
        .ok_or_else(|| "Standard_DomainError".to_string())?;
    let nbpt = wlin.nb_pnts();
    let mut param_search_inf = 1;
    let mut param_search_sup = nbpt;

    if (param_approche - 2) > param_search_inf {
        param_search_inf = param_approche - 2;
    }
    if (param_approche + 2) < param_search_sup {
        param_search_sup = param_approche + 2;
    }

    // first search inside a close bounding around ParamApproche;
    // then search to the nearest end of line;
    // and then search to the farthest end of line.
    let mut inf = [0i32; 3];
    let mut sup = [0i32; 3];
    inf[0] = param_search_inf;
    sup[0] = param_search_sup;
    if param_search_inf - 1 < nbpt - param_search_sup {
        inf[1] = 1;
        sup[1] = param_search_inf;
        inf[2] = param_search_sup;
        sup[2] = nbpt;
    } else {
        inf[1] = param_search_sup;
        sup[1] = nbpt;
        inf[2] = 1;
        sup[2] = param_search_inf;
    }

    for is in 0..3 {
        let mut p1 = wlin.point(inf[is]).value();
        let mut v1 = GpVec::from_pnts(ptsom, &p1);
        let mut norm1 = v1.square_magnitude();
        let mut normmin = tol2;
        let mut ibest = 0;
        let mut param = 0.0;
        if norm1 <= normmin {
            normmin = norm1;
            ibest = inf[is];
        }
        let mut found = false;
        for i in (inf[is] + 1)..=sup[is] {
            let p2 = wlin.point(i).value();
            let v2 = GpVec::from_pnts(ptsom, &p2);
            let norm2 = v2.square_magnitude();
            if v1.dot(&v2) < 0.0 {
                param = (i - 1) as f64 + 1.0 / (1.0 + (norm2 / norm1).sqrt());
                tgl = GpVec::from_pnts(&p1, &p2);
                found = true;
                break;
            } else if norm2 < normmin {
                normmin = norm2;
                ibest = i;
            }
            v1 = v2;
            p1 = p2;
            norm1 = norm2;
        }
        if !found && ibest != 0 {
            param = ibest as f64;
            found = true;
        }
        if found {
            return Ok(Some((param, tgl)));
        }
    }
    Ok(None)
}
