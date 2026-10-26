//! The vertex write-back arms of `IntPatch_RstInt::PutVertexOnLine`.
//!
//! Source: `IntPatch_RstInt.cxx:978-1227`.
//!
//! * `append_new_vertex` is the `if (!ivtx)` arm (`cxx:980-1017`): a fresh
//!   `IntPatch_Point` carrying the refined point, its arc and its transitions.
//! * `update_existing_vertex` is the `else` arm (`cxx:1018-1226`): the refined
//!   point coincides with a vertex already on the line, so only the missing
//!   domain / arc / vertex slots are filled in, and the `SetArc` transplant is
//!   propagated to every occurrence of the same domain vertex.
//!
//! Translation notes:
//! * `occ::handle<...>` comparisons (`ptline.ArcOnS1() == arc`,
//!   `Domain->Identical`) are the port's arc identity token
//!   ([`crate::int_tools_wline::ArcId`]) and `TopolTool::identical`.
//! * OCCT's `IntPatch_Point` keeps both transitions (`traline1` / `tra1`);
//!   the port's [`PatchPoint`] collapses them onto one `trans1` / `trans2`
//!   slot (see `int_tools_wline::PatchPoint`), so `SetArc` writes the line
//!   transition only.
//! * `Nbvtx` is passed by value and returned: two of the arms append to the
//!   vertex sequence while iterating over it (`cxx:1121`, `:1206`).

use occt_core::gp::GpPnt;

use crate::brep_top_adaptor_hvertex::HVertex;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::geom_int::GeomIntLine;
use crate::int_tools_wline::{ArcId, PatchPoint, TransType};

/// Everything the two write-back arms read from the refined point.
pub(crate) struct VertexUpdate {
    /// `ptsommet` (`cxx:467`).
    pub ptsommet: GpPnt,
    /// `vtxTol`.
    pub vtx_tol: f64,
    /// Parameters on S1 / S2, already ordered by `OnFirst`.
    pub u1: f64,
    pub v1: f64,
    pub u2: f64,
    pub v2: f64,
    /// `paramline`.
    pub param_line: f64,
    /// `paramarc`.
    pub param_arc: f64,
    /// `transline` (`cxx:475`).
    pub trans_line: TransType,
    /// `vtxarc` / `VtxOnArc` (`cxx:473-475`).
    pub vtx_arc: Option<HVertex>,
    pub vtx_on_arc: bool,
    /// Identity of the restriction arc `arc` (`cxx:475`).
    pub arc_id: ArcId,
    /// `OnFirst` of the call.
    pub on_first: bool,
}

fn line_add_vertex(line: &mut GeomIntLine, p: PatchPoint) -> Result<(), String> {
    if let Some(w) = line.as_wline_mut() {
        w.vertices.push(p);
        return Ok(());
    }
    if let Some(r) = line.as_rline_mut() {
        r.add_vertex(p);
        return Ok(());
    }
    Err("Standard_DomainError".to_string())
}

fn line_replace(line: &mut GeomIntLine, index: i32, p: PatchPoint) -> Result<(), String> {
    if let Some(w) = line.as_wline_mut() {
        w.replace(index, p);
        return Ok(());
    }
    if let Some(r) = line.as_rline_mut() {
        r.replace(index, p);
        return Ok(());
    }
    Err("Standard_DomainError".to_string())
}

fn line_vertex(line: &GeomIntLine, index: i32) -> Result<PatchPoint, String> {
    line.vertex(index)
        .cloned()
        .ok_or_else(|| "Standard_DomainError".to_string())
}

/// The `if (!ivtx)` arm (`cxx:980-1017`).
pub(crate) fn append_new_vertex(line: &mut GeomIntLine, u: &VertexUpdate) -> Result<(), String> {
    let mut sommet = PatchPoint::new(GpPnt::zero(), 0.0, 0.0, 0.0, 0.0, 0.0);
    // `Sommet.SetValue(ptsommet, vtxTol, false)` ("pour tangence", `cxx:982`).
    sommet.set_value(u.ptsommet, u.vtx_tol, false);
    if u.on_first {
        sommet.set_parameters(u.u1, u.v1, u.u2, u.v2);
    } else {
        sommet.set_parameters(u.u2, u.v2, u.u1, u.v1);
    }
    if u.vtx_on_arc {
        if let Some(v) = &u.vtx_arc {
            sommet.set_vertex(u.on_first, v.clone());
        }
    }
    sommet.set_parameter(u.param_line);
    sommet.set_arc(u.on_first, u.arc_id, u.param_arc, u.trans_line);
    line_add_vertex(line, sommet)
}

/// The `else` arm (`cxx:1018-1226`), i.e. the refined point coincides with the
/// vertex `ivtx` (`ptline` is that vertex).
///
/// `_domain` is kept because OCCT passes `Domain` into the arm; the only use it
/// would have, `Adaptor3d_TopolTool::Identical(V1, V2)`, is
/// `TopolTool::identical`, a static of the port.
pub(crate) fn update_existing_vertex(
    line: &mut GeomIntLine,
    u: &VertexUpdate,
    ptline: PatchPoint,
    ivtx: i32,
    mut nbvtx: i32,
    _domain: &TopolTool,
) -> Result<i32, String> {
    let on_first = u.on_first;
    let mut ptline = ptline;
    let on_different_rst = if on_first {
        ptline.is_on_dom_s1() && ptline.arc_id(true) != u.arc_id
    } else {
        ptline.is_on_dom_s2() && ptline.arc_id(false) != u.arc_id
    };
    ptline.set_tolerance(u.vtx_tol);

    let first_branch = (!ptline.is_vertex_on_s1() && on_first)
        || (!ptline.is_vertex_on_s2() && !on_first)
        || on_different_rst;

    if first_branch {
        let missing_domain = (!ptline.is_on_dom_s2() && !on_first)
            || (!ptline.is_on_dom_s1() && on_first)
            || on_different_rst;
        if missing_domain {
            ptline.set_arc(on_first, u.arc_id, u.param_arc, u.trans_line);
            if u.vtx_on_arc {
                if let Some(v) = &u.vtx_arc {
                    ptline.set_vertex(on_first, v.clone());
                }
            }
            if on_different_rst {
                line_add_vertex(line, ptline)?;
                nbvtx += 1;
            } else {
                line_replace(line, ivtx, ptline)?;
            }
        } else if (on_first && ptline.is_vertex_on_s2())
            || (!on_first && ptline.is_vertex_on_s1())
        {
            let mut sommet = ptline.clone();
            sommet.set_arc(on_first, u.arc_id, u.param_arc, u.trans_line);
            if u.vtx_on_arc {
                if let Some(v) = &u.vtx_arc {
                    sommet.set_vertex(on_first, v.clone());
                }
            }
            line_add_vertex(line, sommet)?;
            nbvtx += 1;
        }
        // else: the source prints "pb dans RstInt Type 1" and does nothing
        // (`cxx:1097-1101`).
        return Ok(nbvtx);
    }

    let vtxref = if on_first {
        ptline.vertex_handle(true)
    } else {
        ptline.vertex_handle(false)
    }
    .cloned()
    .ok_or_else(|| "Standard_DomainError".to_string())?;

    let other_domain_missing =
        (on_first && !ptline.is_on_dom_s2()) || (!on_first && !ptline.is_on_dom_s1());

    if other_domain_missing {
        ptline.set_arc(on_first, u.arc_id, u.param_arc, u.trans_line);
        if u.vtx_on_arc {
            if let Some(v) = &u.vtx_arc {
                ptline.set_vertex(on_first, v.clone());
            }
        }
        line_replace(line, ivtx, ptline)?;

        for k in 1..=nbvtx {
            if k == ivtx {
                continue;
            }
            let mut p = line_vertex(line, k)?;
            let is_vtx = if on_first {
                p.is_vertex_on_s1()
            } else {
                p.is_vertex_on_s2()
            };
            if !is_vtx {
                continue;
            }
            let same = {
                let other = if on_first {
                    p.vertex_handle(true)
                } else {
                    p.vertex_handle(false)
                };
                match other {
                    Some(other) => TopolTool::identical(&vtxref, other),
                    None => false,
                }
            };
            if same {
                if p.tolerance < u.vtx_tol {
                    p.set_tolerance(u.vtx_tol);
                }
                p.set_arc(on_first, u.arc_id, u.param_arc, u.trans_line);
                if u.vtx_on_arc {
                    if let Some(v) = &u.vtx_arc {
                        p.set_vertex(on_first, v.clone());
                    }
                }
                line_replace(line, k, p)?;
            }
        }
    } else if (on_first && ptline.is_vertex_on_s2())
        || (!on_first && ptline.is_vertex_on_s1())
    {
        // "on doit avoir vtxons2 = vtxarc... pas de verif..." (`cxx:1162`):
        // the source does not call `SetVertex` on this copy.
        let mut sommet = ptline.clone();
        sommet.set_arc(on_first, u.arc_id, u.param_arc, u.trans_line);
        line_add_vertex(line, sommet)?;
        nbvtx += 1;

        let mut k = 1;
        while k <= nbvtx {
            if k != ivtx {
                let mut p = line_vertex(line, k)?;
                let is_vtx = if on_first {
                    p.is_vertex_on_s1()
                } else {
                    p.is_vertex_on_s2()
                };
                if is_vtx {
                    let same = {
                        let other = if on_first {
                            p.vertex_handle(true)
                        } else {
                            p.vertex_handle(false)
                        };
                        match other {
                            Some(other) => TopolTool::identical(&vtxref, other),
                            None => false,
                        }
                    };
                    if same {
                        if p.tolerance < u.vtx_tol {
                            p.set_tolerance(u.vtx_tol);
                        }
                        let mut sommet = p.clone();
                        sommet.set_arc(on_first, u.arc_id, u.param_arc, u.trans_line);
                        line_replace(line, k, p)?;
                        line_add_vertex(line, sommet)?;
                        nbvtx += 1;
                    }
                }
            }
            k += 1;
        }
    }
    // else: the source prints "pb dans RstInt Type 2" and does nothing
    // (`cxx:1220-1224`).

    Ok(nbvtx)
}
