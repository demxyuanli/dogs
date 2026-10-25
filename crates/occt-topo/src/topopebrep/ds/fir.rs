//! `TopOpeBRepDS_FIR` — the face-interference reducer.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_FIR.hxx:30-50`,
//! `TopOpeBRepDS_FIR.cxx:36-571` (OCCT 8.0.0).
//!
//! The translation unit carries four `Standard_EXPORT` free functions that are
//! re-declared by the other reducers; only `FUN_Parameters` (`FaceInterferenceTool.cxx:46-68`)
//! and `FUN_edgeofface` (`ProcessFaceInterferences.cxx:146-...`) are declared here
//! and both are defined elsewhere, so this file does not redefine them.
//!
//! Translation notes:
//! - OCCT passes `NCollection_List&` and mutates it while iterating with list
//!   iterators. The port indexes a `Vec`, so `List::Remove(it)` becomes
//!   `remove(j)` without advancing `j` (`cxx:307`), and `it.Next()` becomes
//!   `i += 1`.
//! - The `Standard_ProgramError` / `Standard_NoSuchObject` throws of
//!   `TopoDS::Face(BDS.Shape(I))` are the `Err` arms of the port.
//! - `FUN_stateedgeface` (`cxx:344-361`) classifies the UV of a point in the
//!   face with `TopOpeBRepTool_SC::FSC_GetPSC(F)` +
//!   `TopOpeBRepTool_ShapeClassifier::StateP2DReference` (`ShapeClassifier.cxx:585-637`).
//!   That method has two arms: the `mymren == 1` one (`:596-621`) and
//!   `BRepTopAdaptor_FClass2d(F, 1e-8).Perform(P2D)` (`:622-632`). Only the
//!   second is ported here (`FClass2d::perform_tab_orien`, the same
//!   `BRepTopAdaptor_FClass2d::Perform` body at `BRepTopAdaptor_FClass2d.cxx:525-665`);
//!   the `mymren == 1` arm is UNPORTED because `mymren` is the edge count of the
//!   referenced face *conditioned on the classifier singleton's `mySameDomain`
//!   being 1* (`ShapeClassifier.cxx:118-139`), and `TopOpeBRepTool_ShapeClassifier`
//!   itself is not ported (it also needs `BRepClass3d_SolidClassifier`).
//! - The `FUN_reduceEDGEgeometry1` "reverse the transitions" blocks
//!   (`cxx:231-292`, `:300-305`) are commented out in OCCT and stay commented
//!   out here; `FUN_reduceEDGEgeometry1` therefore starts from
//!   `STATIC_TOREVERSE = false` inside `FaceInterferenceTool` exactly as
//!   `FaceInterferenceTool::new` does.

use occt_core::gp::{GpPnt, GpPnt2d};

use crate::abs::{ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::fclass2d::{FaceState, FClass2d};
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{shape_map_key, ShapeMapKey};
use crate::topopebrep::tool::adaptor::CurveAdaptor;
use crate::topopebrep::tool::topology::{
    fun_tool_find_p_in_bac, fun_tool_find_p_in_e, fun_tool_orient_ein_f_forward,
};

use super::data_structure::DataStructure;
use super::export::{fds_data, fds_data_is_unknown};
use super::face_interference_tool::{fun_parameters, FaceInterferenceTool};
use super::h_data_structure::HDataStructure;
use super::interference::{Interference, InterferenceData};
use super::interference_iterator::InterferenceList;
use super::kind::Kind;
use super::process_interferences_tool::{
    fun_reducedoublons, fun_select_gk_interference, fun_select_sk_interference,
    fun_select_trasha_interference,
};
use super::tool::DShapeOn1StateMap;

/// `#define MDSke TopOpeBRepDS_EDGE` (`cxx:30`).
const MDS_KE: Kind = Kind::Edge;
/// `#define MDSkf TopOpeBRepDS_FACE` (`cxx:31`).
const MDS_KF: Kind = Kind::Face;

/// `BDS.Shape(S1)` followed by `TopoDS::Face` (`cxx:52`, `:268`, `:287`, `:498`).
fn face_at(bds: &DataStructure, i: i32) -> Result<Face, String> {
    let s = bds
        .shape(i, true)
        .ok_or_else(|| "TopOpeBRepDS_FIR: no such shape".to_string())?;
    if s.shape_type() != ShapeType::Face {
        return Err("TopOpeBRepDS_FIR: TopoDS::Face".to_string());
    }
    Ok(Face(s))
}

/// `TopoDS::Vertex(BDS.Shape(G1))` (`cxx:115`, `:128`).
fn vertex_at(bds: &DataStructure, i: i32) -> Result<Vertex, String> {
    let s = bds
        .shape(i, true)
        .ok_or_else(|| "TopOpeBRepDS_FIR: no such shape".to_string())?;
    if s.shape_type() != ShapeType::Vertex {
        return Err("TopOpeBRepDS_FIR: TopoDS::Vertex".to_string());
    }
    Ok(Vertex(s))
}

/// `FUN_isPonF(LIF, P, BDS, E)` (`cxx:38-70`): `P` is on every face of `LIF`.
///
/// The loop `break`s on the first face `P` is not on, and `Pok` is the answer
/// of the last test otherwise (`cxx:54-67`).
pub fn fun_is_pon_f(
    lif: &InterferenceList,
    p: &GpPnt,
    bds: &DataStructure,
    e: &TopoShape,
) -> Result<bool, String> {
    let mut pok = true;
    for h in lif {
        let (_, _, _, s1) = {
            let i = h.read().expect("interference lock");
            fds_data(&i)
        };
        let f = face_at(bds, s1)?;
        if fun_tool_orient_ein_f_forward(e, &f.0).is_some() {
            pok = true;
        } else {
            // P est NOK pour une face de LIF : arret
            let ok = fun_parameters(p, &f).is_some();
            pok = ok;
            if !ok {
                break;
            }
        }
    }
    Ok(pok)
}

/// `FUN_findPonF(E, BDS, LIF, P, par)` (`cxx:72-156`): the point of `E` that
/// lies on every face of `LIF`, with its parameter on `E`.
///
/// The per-interference loop reads the whole `ShapeInterferences(E)` list, so a
/// `FDS_data` on an element with no parameter falls through to `continue`
/// (`cxx:137-142`).
pub fn fun_find_pon_f(
    e: &TopoShape,
    bds: &DataStructure,
    lif: &InterferenceList,
) -> Result<(GpPnt, f64, bool), String> {
    let bc = CurveAdaptor::of_edge(&Edge(e.clone()));
    let lie = bds.shape_interferences(e, true).to_vec();

    if lie.is_empty() {
        let (p, par) =
            fun_tool_find_p_in_bac(&bc).ok_or_else(|| "FUN_tool_findPinBAC".to_string())?;
        let pok = fun_is_pon_f(lif, &p, bds, e)?;
        return Ok((p, par, pok));
    }

    let mut p = GpPnt::zero();
    let mut par = 0.0;
    let mut pok = false;
    for h in &lie {
        // `bool pardef = false;` (`cxx:96`): every path below assigns before read.
        let pardef;
        let (gt1, g1, _, _) = {
            let i = h.read().expect("interference lock");
            fds_data(&i)
        };
        let is_curve_point = {
            let i = h.read().expect("interference lock");
            matches!(i.data(), InterferenceData::CurvePoint { .. })
        };
        let ssi_gb = {
            let i = h.read().expect("interference lock");
            i.shape_shape_data().map(|d| d.g_bound)
        };
        if is_curve_point {
            par = {
                let i = h.read().expect("interference lock");
                i.parameter()
                    .ok_or_else(|| "FDS_Parameter".to_string())?
            };
            pardef = true;
        } else if let Some(gb) = ssi_gb {
            if gb {
                let v = vertex_at(bds, g1)?;
                p = BRepTool::vertex_point_world(&v);
                par = BRepTool::parameter_of_vertex(e, &v.0);
                pardef = true;
            } else {
                // `pardef = false;` (`cxx:120`) then `if (pardef)` at `cxx:130`
                // tests that very value, so `FUN_tool_projPonC` (`cxx:132-135`)
                // never runs. Kept verbatim.
                pardef = false;
                if gt1 == Kind::Point {
                    let pd = bds
                        .point(g1)
                        .ok_or_else(|| "TopOpeBRepDS_DataStructure::Point".to_string())?;
                    p = pd.point().clone();
                } else if gt1 == Kind::Vertex {
                    let v = vertex_at(bds, g1)?;
                    p = BRepTool::vertex_point_world(&v);
                }
                if pardef {
                    // UNPORTED (never reached): `FUN_tool_projPonC(P, BAC, par, dist)`
                    // at `TopOpeBRepDS_FIR.cxx:132-135`.
                }
            }
        } else {
            continue;
        }

        if !pardef {
            continue;
        }

        p = bc
            .d0(par)
            .ok_or_else(|| "BRepAdaptor_Curve::D0".to_string())?;
        pok = fun_is_pon_f(lif, &p, bds, e)?;
        // P est OK pour toutes les faces de LIF : on arrete de chercher
        if pok {
            break;
        }
    }
    Ok((p, par, pok))
}

/// `FDS_ADDEDGE(TRCF, S, iFI, FITool, FI, F, Ecpx, isEGsp, I)`
/// (`cxx:158-171`): the three leading arguments are debug-only and forwarded
/// to `FITool.Add` as the same `(FI, F, Ecpx, isEGsp, I)`.
fn fds_add_edge(
    fitool: &mut FaceInterferenceTool<'_>,
    fi: &TopoShape,
    f: &TopoShape,
    ecpx: &TopoShape,
    is_egsp: bool,
    i: &Interference,
) -> Result<(), String> {
    fitool.add(fi, f, ecpx, is_egsp, i)
}

/// `FUN_reduceEDGEgeometry1(LI, BDS, iFI, iEG, EGsp, MEsp)` (`cxx:173-315`):
/// the complex transition of every pair of interferences of `LI` sharing the
/// same edge geometry, resolved on the edge `Ecpx`.
///
/// `MEsp` is unused by OCCT (`cxx:178-184`) and is not part of the port's
/// signature.
pub fn fun_reduce_edge_geometry1(
    li: &mut InterferenceList,
    bds: &DataStructure,
    ifi: i32,
    ieg: i32,
    egsp: Option<&TopoShape>,
) -> Result<(), String> {
    if li.is_empty() {
        return Ok(());
    }

    // choix de l'arete Ecpx, lieu de resolution de la transition complexe
    let fi_shape = bds
        .shape(ifi, true)
        .ok_or_else(|| "TopOpeBRepDS_FIR: no such shape".to_string())?;
    if fi_shape.shape_type() != ShapeType::Face {
        return Err("TopOpeBRepDS_FIR: TopoDS::Face".to_string());
    }
    let fi = Face(fi_shape);
    let is_egsp = egsp.is_some();
    let ecpx: TopoShape = match egsp {
        Some(s) => {
            if s.shape_type() != ShapeType::Edge {
                return Err("TopOpeBRepDS_FIR: TopoDS::Edge".to_string());
            }
            s.clone()
        }
        None => {
            let s = bds
                .shape(ieg, true)
                .ok_or_else(|| "TopOpeBRepDS_FIR: no such shape".to_string())?;
            if s.shape_type() != ShapeType::Edge {
                return Err("TopOpeBRepDS_FIR: TopoDS::Edge".to_string());
            }
            s
        }
    };

    let mut fitool = FaceInterferenceTool::new(bds);
    if li.len() >= 2 {
        let (pok, par_pok) = if is_egsp {
            fun_tool_find_p_in_e(&ecpx).ok_or_else(|| "FUN_tool_findPinE".to_string())?
        } else {
            let (p, par, ok) = fun_find_pon_f(&ecpx, bds, li)?;
            if !ok {
                // `if (!isPok) { LI.Clear(); return; }` (`cxx:220-224`)
                li.clear();
                return Ok(());
            }
            (p, par)
        };
        fitool.set_edge_pnt_par(pok, par_pok);
    }

    // xpu :090498 : the `toreverse` blocks of `cxx:231-243` and the two
    // `Complement()` blocks of `cxx:292-305` / `:301-304` are commented out in
    // OCCT.

    let mut i1 = 0;
    while i1 < li.len() {
        let unknown1 = {
            let h = li[i1].read().expect("interference lock");
            fds_data_is_unknown(&h)
        };
        if unknown1 {
            i1 += 1;
            continue;
        }
        let (gt1, g1, st1, s1) = {
            let h = li[i1].read().expect("interference lock");
            fds_data(&h)
        };
        if gt1 != MDS_KE {
            i1 += 1;
            continue;
        }

        let mut is_complex = false; // True if at least two interfs on the same edge
        let f1 = face_at(bds, s1)?;

        let mut i2 = i1 + 1;
        while i2 < li.len() {
            let unknown2 = {
                let h = li[i2].read().expect("interference lock");
                fds_data_is_unknown(&h)
            };
            if unknown2 {
                i2 += 1;
                continue;
            }
            let (gt2, g2, st2, s2) = {
                let h = li[i2].read().expect("interference lock");
                fds_data(&h)
            };
            if !(gt2 == gt1 && g2 == g1 && st2 == st1) {
                i2 += 1;
                continue;
            }
            let f2 = face_at(bds, s2)?;
            if !is_complex {
                is_complex = true;
                let h = li[i1].read().expect("interference lock");
                fitool.init(&fi.0, &ecpx, is_egsp, &h)?;
                fds_add_edge(&mut fitool, &fi.0, &f1.0, &ecpx, is_egsp, &h)?;
            }
            {
                let h = li[i2].read().expect("interference lock");
                fds_add_edge(&mut fitool, &fi.0, &f2.0, &ecpx, is_egsp, &h)?;
            }

            li.remove(i2);
        }
        if is_complex {
            let mut h = li[i1].write().expect("interference lock");
            fitool.transition(&mut h);
        }
        i1 += 1;
    }
    Ok(())
}

/// `FUN_GmapS(LI, BDS, mosd)` (`cxx:317-346`): group the `LI` elements whose
/// geometry is an edge and whose support is a face, keyed by the geometry
/// shape, in first-seen order (OCCT's `IndexedDataMap` order).
fn fun_gmap_s(
    li: &InterferenceList,
    bds: &DataStructure,
) -> (Vec<TopoShape>, Vec<InterferenceList>) {
    let mut keys: Vec<TopoShape> = Vec::new();
    let mut index: std::collections::HashMap<ShapeMapKey, usize> = std::collections::HashMap::new();
    let mut lists: Vec<InterferenceList> = Vec::new();
    for h in li {
        let (gt1, g1, st1, _) = {
            let i = h.read().expect("interference lock");
            fds_data(&i)
        };
        if gt1 != MDS_KE || st1 != MDS_KF {
            continue;
        }
        let Some(sg1) = bds.shape(g1, true) else {
            // `BDS.Shape(G1)` would raise; a valid DS always has it.
            continue;
        };
        let pos = match index.get(&shape_map_key(&sg1)) {
            Some(p) => *p,
            None => {
                keys.push(sg1.clone());
                lists.push(InterferenceList::new());
                let p = keys.len() - 1;
                index.insert(shape_map_key(&sg1), p);
                p
            }
        };
        lists[pos].push(h.clone());
    }
    (keys, lists)
}

/// `FUN_stateedgeface(E, F, P)` (`cxx:344-361`): the state of the point of `E`
/// at `FUN_tool_findPinE`'s parameter in the face `F`.
///
/// `TopOpeBRepTool_ShapeClassifier::StateP2DReference` (`ShapeClassifier.cxx:585-637`)
/// is reproduced by the `BRepTopAdaptor_FClass2d(F, 1e-8).Perform(P2D)` arm;
/// the `mymren == 1` arm (`:596-621`) is UNPORTED (see the module note).
pub fn fun_state_edge_face(e: &TopoShape, f: &TopoShape) -> Result<State, String> {
    let (p, _par) = fun_tool_find_p_in_e(e).ok_or_else(|| "FUN_tool_findPinE".to_string())?;
    let face = face_at_shape(f)?;
    let Some((u, v)) = fun_parameters(&p, &face) else {
        // `bool Pok = FUN_Parameters(P, F, u, v); if (Pok) {...}` (`cxx:354-359`)
        return Ok(State::Unknown);
    };
    // UNPORTED: the `mymren == 1` arm of
    // `TopOpeBRepTool_ShapeClassifier::StateP2DReference` (`ShapeClassifier.cxx:596-621`).
    let classifier = FClass2d::new(&face, 1e-8)?;
    Ok(match classifier.perform_tab_orien(GpPnt2d::new(u, v)) {
        FaceState::In => State::In,
        FaceState::Out => State::Out,
        FaceState::On => State::On,
        FaceState::Unknown => State::Unknown,
    })
}

/// `TopoDS::Face(F)` on an already-held shape (`cxx:351`).
fn face_at_shape(f: &TopoShape) -> Result<Face, String> {
    if f.shape_type() != ShapeType::Face {
        return Err("TopOpeBRepDS_FIR: TopoDS::Face".to_string());
    }
    Ok(Face(f.clone()))
}

/// `FUN_reduceEDGEgeometry(LI, BDS, iFI, MEsp)` (`cxx:368-513`).
///
/// The interferences of `LI` are grouped by geometry edge; a group of one is
/// forwarded as is, a group of two or more is resolved either per split of the
/// geometry edge (`isEGsp`) or as a whole (`cxx:498-508`).
pub fn fun_reduce_edge_geometry(
    li: &mut InterferenceList,
    bds: &DataStructure,
    ifi: i32,
    me_sp: &DShapeOn1StateMap,
) -> Result<(), String> {
    if li.is_empty() {
        return Ok(());
    }

    let (keys, mut lists) = fun_gmap_s(li, bds);

    let mut li_out = InterferenceList::new();
    for idx in 0..keys.len() {
        let eg = keys[idx].clone();
        let ieg = bds.shape_index(&eg, true);

        // donnees samedomain attachees a l'arete iEG
        let esdeg = bds.shape_same_domain_at(ieg).to_vec();
        let egissect = if eg.shape_type() == ShapeType::Edge {
            bds.is_section_edge(&eg, true)
        } else {
            // `TopoDS::Edge(EG)` (`cxx:397`) would raise.
            return Err("TopOpeBRepDS_FIR: TopoDS::Edge".to_string());
        };
        let eghasesd = !esdeg.is_empty();

        let n_ext = lists[idx].len();
        // LIEG = toutes les interferences dont le Support() est une
        // face possedant une interference dont la Geometry() est EG.
        if n_ext == 0 {
            continue;
        }
        if n_ext == 1 {
            li_out.extend(lists[idx].iter().cloned());
        } else if n_ext >= 2 {
            let mut is_egsp = me_sp.contains_key(&shape_map_key(&eg));
            if is_egsp {
                let los1 = me_sp
                    .get(&shape_map_key(&eg))
                    .ok_or_else(|| "FUN_reduceEDGEgeometry: no MEsp entry".to_string())?;
                is_egsp = los1.is_split();
            }

            if is_egsp {
                let los: Vec<TopoShape> = me_sp
                    .get(&shape_map_key(&eg))
                    .map(|l| l.list_on_state().to_vec())
                    .unwrap_or_default();
                for egsp in &los {
                    // LISFIN = liste des interferences de LI dont le Support()
                    // est une face contenant geometriquement l'arete EGsp
                    let mut lisfin = InterferenceList::new();
                    for ilieg in lists[idx].iter() {
                        let (is_, ks) = {
                            let i = ilieg.read().expect("interference lock");
                            (i.support(), i.support_type())
                        };
                        if ks == MDS_KF {
                            let Some(sfilieg) = bds.shape(is_, true) else {
                                continue;
                            };
                            let staef = fun_state_edge_face(egsp, &sfilieg)?;

                            let mut pok = staef == State::In;
                            if eghasesd || egissect {
                                pok = pok || staef == State::On;
                            }

                            if pok {
                                lisfin.push(ilieg.clone());
                            }
                        }
                    } // itLIEG.More

                    if lisfin.len() >= 2 {
                        let gb = {
                            let h = lisfin[0].read().expect("interference lock");
                            h.shape_shape_data().map(|d| d.g_bound).unwrap_or(false)
                        };

                        if gb {
                            // modified by NIZNHY-PKV Thu Mar 16 10:40:57 2000 f
                            //  we have to rest at least one  Interference on the face.
                            //  To kill all of them is too bravely.
                            let an_interference = lisfin[0].clone();
                            lisfin.clear();
                            lisfin.push(an_interference);
                            // modified by NIZNHY-PKV Thu Mar 16 10:41:01 2000 t
                        } else {
                            fun_reduce_edge_geometry1(
                                &mut lisfin,
                                bds,
                                ifi,
                                ieg,
                                Some(egsp),
                            )?;
                        }
                    }

                    if !lisfin.is_empty() {
                        li_out.extend(lisfin);
                    }
                }
            } // isEGsp
            else {
                // iFI = face de reference (indice)
                // E = arete geometrie d'interference (shape), iEG (indice)
                // LIEG = liste d'interferences de geometrie EG
                //        et dont les Support() sont a transitionner complexe
                fun_reduce_edge_geometry1(&mut lists[idx], bds, ifi, ieg, None)?;
                li_out.extend(lists[idx].iter().cloned());
            }
        }
    }

    li.clear();
    li.extend(li_out);
    Ok(())
}

/// `TopOpeBRepDS_FIR` (`hxx:30-50`): the reducer, holding the HDS it rewrites.
pub struct Fir<'a> {
    /// `myHDS`
    hds: &'a mut HDataStructure,
}

impl<'a> Fir<'a> {
    /// `TopOpeBRepDS_FIR(HDS)` (`cxx:517-520`).
    pub fn new(hds: &'a mut HDataStructure) -> Self {
        Self { hds }
    }

    /// `ProcessFaceInterferences(M)` (`cxx:524-544`): every face shape of the
    /// DS, in index order.
    pub fn process_all_face_interferences(
        &mut self,
        m: &DShapeOn1StateMap,
    ) -> Result<(), String> {
        let nshape = self.hds.ds().nb_shapes();
        for i in 1..=nshape {
            let s = self.hds.ds().shape(i, true);
            let Some(s) = s else {
                continue;
            };
            if s.shape_type() == ShapeType::Face {
                self.process_face_interferences(i, m)?;
            }
        }
        Ok(())
    }

    /// `ProcessFaceInterferences(SIX, MEsp)` (`cxx:548-571`).
    ///
    /// `lw` keeps every interference the two `FUN_selectTRASHAinterference`
    /// passes did not move; OCCT drops it (`cxx:563-568` rebuilds `LI` from
    /// `lF`, `lFE`, `lFEF`, `lE`), so the port drops it too.
    pub fn process_face_interferences(
        &mut self,
        six: i32,
        me_sp: &DShapeOn1StateMap,
    ) -> Result<(), String> {
        let li_orig: InterferenceList = self
            .hds
            .ds()
            .shape_interferences_at(six, true)
            .to_vec();
        let mut lw = li_orig;
        let mut lf = InterferenceList::new();
        let mut lfe = InterferenceList::new();
        let mut lfef = InterferenceList::new();
        let mut le = InterferenceList::new();

        fun_select_trasha_interference(&mut lw, ShapeType::Face, &mut lf);
        fun_select_gk_interference(&mut lf, MDS_KE, &mut lfe);
        fun_select_sk_interference(&mut lfe, MDS_KF, &mut lfef);
        fun_select_trasha_interference(&mut lw, ShapeType::Edge, &mut le);
        fun_reduce_edge_geometry(&mut lfef, self.hds.ds(), six, me_sp)?;

        // `LI.Clear(); LI.Append(lF); ... LI.Append(lE);` (`cxx:563-568`)
        let mut li: InterferenceList = lf;
        li.extend(lfe);
        li.extend(lfef);
        li.extend(le);
        // MSV: filter duplicates (`cxx:570`: `FUN_reducedoublons(LI, BDS, SIX)`)
        fun_reducedoublons(&mut li, self.hds.ds(), six);

        let ds_li = self
            .hds
            .change_ds()
            .change_shape_interferences_at(six);
        ds_li.clear();
        ds_li.extend(li);
        Ok(())
    }
}