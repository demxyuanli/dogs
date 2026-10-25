//! `FUN_ReducerEdge3d` — the 3d complex-transition reducer of
//! `TopOpeBRepDS_EIR`.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_EIR.cxx:107-375` (OCCT 8.0.0).
//!
//! `LI = {I3d = (TonFACE, G=POINT/VERTEX, S=EDGE)}` is folded pairwise with
//! [`Edge3dInterferenceTool`]: the first pair of a geometry group primes the
//! tool (`InitPointVertex` / `Init` / `Add`, `cxx:301-324`), every later pair is
//! `Add`ed (`cxx:326`), and the resulting transition is written into a new
//! `(TonFACE, G, S=FACE)` interference (`cxx:331-364`).
//!
//! Translation notes:
//! - `S2` is redeclared at every `it2` iteration (`cxx:217`) and read only
//!   before `S1 = S2 = BDS.Shape(Eshared)` (`cxx:286`), so the store on `S2` is
//!   dead in OCCT; `S1` persists across the `it2` iterations and is kept in the
//!   port.
//! - `return` inside the loops (`cxx:209`, `:277`, `:281`) leaves the function
//!   with the reductions done so far; the port returns `Ok(())` at the same
//!   points.
//! - `newI` stays null when `I1` is neither an `EdgeVertexInterference` nor a
//!   `CurvePointInterference`; OCCT then keeps `I1` (`cxx:366-368`).

use crate::abs::{Orientation, ShapeType, State};
use crate::shape::TopoShape;

use super::super::data_structure::DataStructure;
use super::super::edge3d_interference_tool::Edge3dInterferenceTool;
use super::super::export::{fds_data, fds_tdata, fun_ds_get_vsdm, fun_ds_sdm};
use super::super::interference::{Interference, InterferenceData};
use super::super::interference_iterator::InterferenceList;
use super::super::interference_tool::make_handle;
use super::super::kind::{Config, Kind};
use super::check_interference_is_valid;
use crate::topopebrep::tool::{fun_tool_eshared, fun_tool_orient_ein_f_forward};

/// `BDS.Shape(I)` for a face index (`cxx:265-266`).
fn shape_at(bds: &DataStructure, i: i32, who: &str) -> Result<TopoShape, String> {
    bds.shape(i, true)
        .ok_or_else(|| format!("FUN_ReducerEdge3d: BDS.Shape({who})"))
}

/// `FUN_ReducerEdge3d(SIX, BDS, LI, reducedLI)` (`cxx:107-375`).
pub fn fun_reducer_edge3d(
    six: i32,
    bds: &DataStructure,
    li: &mut InterferenceList,
    reduced_li: &mut InterferenceList,
) -> Result<(), String> {
    reduced_li.clear();
    if li.len() <= 1 {
        return Ok(());
    }

    let e = shape_at(bds, six, "SIX")?;
    let rank_e = bds.ancestor_rank_of(&e);
    let mut oov: Option<TopoShape> = None;
    let mut gsta: i32 = 0;

    // <-- jyl 231198 : cts21797
    if !li.is_empty() {
        let mut i = 0;
        while i < li.len() {
            let (st1, s1, ib1) = {
                let h = li[i].read().expect("interference lock");
                let (_gt, _g, st, s) = fds_data(&h);
                let (_sb, ib, _sa, _ia) = fds_tdata(&h.transition());
                (st, s, ib)
            };
            let mut torem = false;
            if st1 == Kind::Edge {
                let ee = shape_at(bds, s1, "S1")?;
                let ff = shape_at(bds, ib1, "IB1")?;
                let o = fun_tool_orient_ein_f_forward(&ee, &ff);
                if let Some(o) = o {
                    if o == Orientation::External {
                        torem = true;
                    }
                }
            }
            if torem {
                li.remove(i);
            } else {
                i += 1;
            }
        }
        if li.len() <= 1 {
            return Ok(());
        }
    }
    // --> jyl 231198 : cts21797

    let mut i1 = 0;
    while i1 < li.len() {
        let mut is_complex = false;
        let mut efi_tool = Edge3dInterferenceTool::new();

        let h1 = li[i1].clone();
        let (gt1, g1, st1, s1, ib1) = {
            let h = h1.read().expect("interference lock");
            let (gt, g, st, s) = fds_data(&h);
            let (_sb, ib, _sa, _ia) = fds_tdata(&h.transition());
            (gt, g, st, s, ib)
        };
        let mut s1 = s1;

        // modified by NIZHNY-MKK  Mon Apr  2 15:35:58 2001.BEGIN
        let mut a_vertex: Option<TopoShape> = None;
        if gt1 == Kind::Vertex && g1 != 0 {
            a_vertex = Some(shape_at(bds, g1, "G1")?);
        }
        let e_s1 = shape_at(bds, s1, "S1")?;
        let valid1 = {
            let h = h1.read().expect("interference lock");
            check_interference_is_valid(&h, &e, &e_s1, a_vertex.as_ref())?
        };
        if !valid1 {
            li.remove(i1);
            continue;
        }
        // modified by NIZHNY-MKK  Mon Apr  2 15:36:02 2001.END

        if i1 + 1 >= li.len() {
            return Ok(());
        }

        let mut i2 = i1 + 1;
        while i2 < li.len() {
            let h2 = li[i2].clone();
            let (gt2, g2, st2, s2, ib2) = {
                let h = h2.read().expect("interference lock");
                let (gt, g, st, s) = fds_data(&h);
                let (_sb, ib, _sa, _ia) = fds_tdata(&h.transition());
                (gt, g, st, s, ib)
            };
            let mut s2 = s2;

            // modified by NIZHNY-MKK  Mon Apr  2 15:36:42 2001.BEGIN
            a_vertex = None;
            if gt2 == Kind::Vertex && g2 != 0 {
                a_vertex = Some(shape_at(bds, g2, "G2")?);
            }
            let e_s2 = shape_at(bds, s2, "S2")?;
            let valid2 = {
                let h = h2.read().expect("interference lock");
                check_interference_is_valid(&h, &e, &e_s2, a_vertex.as_ref())?
            };
            if !valid2 {
                li.remove(i2);
                continue;
            }
            // modified by NIZHNY-MKK  Mon Apr  2 15:36:45 2001.END

            let same_g = gt2 == gt1 && g2 == g1;
            if !same_g {
                break;
            }

            // <Gsta>, <OOv>
            if gt1 == Kind::Vertex {
                let vg1 = shape_at(bds, g1, "G1")?;
                let rank_vg1 = bds.ancestor_rank_of(&vg1);
                let g1hsd = fun_ds_get_vsdm(bds, g1);
                if rank_vg1 != rank_e {
                    // vG1 not on E
                    oov = Some(vg1);
                    gsta = if g1hsd.is_some() { 3 } else { 2 };
                } else {
                    // vG on E
                    if let Some(sd) = g1hsd {
                        oov = bds.shape(sd, true);
                    }
                    gsta = if g1hsd.is_some() { 3 } else { 1 };
                }
            }

            let f1 = shape_at(bds, ib1, "IB1")?;
            let f2 = shape_at(bds, ib2, "IB2")?; // F2 != F1

            let same_s = st2 == st1 && s2 == s1;
            if !same_s {
                let foundsh = match oov.as_ref() {
                    Some(v) => fun_tool_eshared(v, &f1, &f2),
                    None => None,
                };
                let Some(eshared) = foundsh else {
                    return Ok(());
                };

                // modified by NIZHNY-MKK  Mon Apr  2 15:37:12 2001.BEGIN
                if !bds.has_shape(&eshared, true) {
                    return Ok(());
                }
                // modified by NIZHNY-MKK  Mon Apr  2 15:37:15 2001.END

                s1 = bds.shape_index(&eshared, true);
                s2 = s1;
            }

            let e1 = shape_at(bds, s1, "S1")?;
            let e2 = shape_at(bds, s2, "S2")?; // E2 == E1

            let sdm = fun_ds_sdm(bds, &e, &e1);
            if sdm {
                i2 += 1;
                continue;
            }

            let mut init = !is_complex;
            let is_vertex = gt1 == Kind::Vertex;
            init = init || is_vertex; // !!!KK a revoir!!!!

            if init {
                // xpu : 04-03-98 : !!!KK a revoir!!!!
                if is_complex {
                    let mut ibid = Interference::default();
                    efi_tool.transition(&mut ibid);
                    let o = ibid
                        .transition()
                        .orientation(State::In, ShapeType::Face);
                    h1.write()
                        .expect("interference lock")
                        .change_transition()
                        .set_orientation(o);
                }
                if !is_complex {
                    efi_tool.init_point_vertex(gsta, oov.clone());
                }
                is_complex = true;
                {
                    let h = h1.read().expect("interference lock");
                    efi_tool.init(&e, &e1, &f1, &h)?;
                    efi_tool.add(&e, &e1, &f1, &h)?;
                }
            } // !isComplex

            {
                let h = h2.read().expect("interference lock");
                efi_tool.add(&e, &e2, &f2, &h)?;
            }
            li.remove(i2);
        } // it2

        if is_complex {
            let mut ibid = Interference::default();
            efi_tool.transition(&mut ibid);
            let mut t = ibid.transition();
            t.set_index(ib1);

            let (isevi, iscpi) = {
                let h = h1.read().expect("interference lock");
                (
                    matches!(h.data(), InterferenceData::EdgeVertex(_)),
                    matches!(h.data(), InterferenceData::CurvePoint { .. }),
                )
            };

            let mut new_i: Option<super::super::interference::InterferenceHandle> = None;
            if isevi {
                let (gbound, param) = {
                    let h = h1.read().expect("interference lock");
                    match h.as_edge_vertex() {
                        Some(d) => (d.shape_shape.g_bound, d.param),
                        None => (false, 0.0),
                    }
                };
                new_i = Some(make_handle(Interference::edge_vertex(
                    t.clone(),
                    Kind::Face,
                    ib1,
                    g1,
                    gbound,
                    Config::UnshGeometry,
                    param,
                )));
            }
            if iscpi {
                let param = h1
                    .read()
                    .expect("interference lock")
                    .parameter()
                    .ok_or_else(|| "FUN_ReducerEdge3d: FDS_Parameter".to_string())?;
                new_i = Some(make_handle(Interference::curve_point(
                    t.clone(),
                    Kind::Face,
                    ib1,
                    Kind::Point,
                    g1,
                    param,
                )));
            }

            if let Some(newi) = new_i {
                reduced_li.push(newi);
                li.remove(i1);
            }
        } else {
            i1 += 1;
        }
    } // it1
    Ok(())
}
