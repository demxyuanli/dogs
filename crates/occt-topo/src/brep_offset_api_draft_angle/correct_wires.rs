//! `BRepOffsetAPI_DraftAngle::CorrectWires()`
//! (`BRepOffsetAPI_DraftAngle.cxx:272-909`, OCCT 8.0.0).
//!
//! The draft rebuilds every seam edge into several pieces: the p-curves of the
//! tuned face are cut where the *other* wires of the face cross them, and the
//! seam edges are split at those intersections so each piece has two distinct
//! vertices. `CorrectWires` performs exactly that split, in the reference's
//! eleven phases:
//!
//! | OCCT | lines | phase |
//! |---|---|---|
//! | collect the really-closed edges of every wire | `cxx:282-308` | 1 |
//! | 2D intersection of each seam p-curve with the other wires | `cxx:310-443` | 2 |
//! | sort `ParsNonSeam` (and `Seam` / `ParsSeam` with it) | `cxx:445-464` | 3 |
//! | sort `Emap` | `cxx:466-483` | 4 |
//! | seed `EPmap` / `EVmap` / `EWmap` | `cxx:485-499` | 5 |
//! | split the non-seam edges | `cxx:501-592` | 6 |
//! | sort and dedupe `EPmap` / `EVmap` / `EWmap` | `cxx:594-656` | 7 |
//! | split the seam edges (`VEmap`) | `cxx:658-777` | 8 |
//! | drop edges touching a missing seam extremity | `cxx:779-815` | 9 |
//! | rebuild `myShape` from the substitution | `cxx:817-826` | 10 |
//! | rebuild the wires, shifting p-curves by `2*pi` | `cxx:828-908` | 11 |
//!
//! # Translation notes
//!
//! * `TopoDS_Iterator` on a face or a wire (`cxx:286`, `:339`, `:793`,
//!   `:849`) composes the orientation, which is `explorer_shapes` /
//!   [`edges_of_wire`]; `TopoDS_Iterator(aWire, false)` (`cxx:895`) does not,
//!   which is `children_keep_orientation`.
//! * `NCollection_Map<TopoDS_Shape>` over a wire's edges (`cxx:290-295`) is
//!   a `HashSet` over [`shape_map_key`]; the map's hash order is replaced by
//!   the wire's edge order, which is one of the orders OCCT's map may produce
//!   and is the only place the two differ (the phases that follow sort every
//!   sequence they build).
//! * `BRepAdaptor_Curve2d(E, F)` is [`Curve2dAdaptor::of_edge_face`]; the
//!   reference also builds one over `E.Reversed()` (`cxx:331`), which is the
//!   same p-curve read on the reversed occurrence of the edge.
//! * `Geom2dInt_GInter::Perform(AC1, AC2, Tol, Tol)` (`cxx:355`, `:359`) takes
//!   the adaptors, i.e. the p-curve narrowed to the adaptor's own window, so
//!   the port computes both `IntRes2d_Domain`s from
//!   [`Curve2dAdaptor::first_parameter`] / `last_parameter` through
//!   [`Geom2dIntGInter::domain_over_range`] before `perform`.
//!   `IsEmpty()` is `result().is_empty()`.
//! * `BRep_Builder::Range(E, First, Last)` is
//!   [`GeometryRegistry::set_edge_range`]; `UpdateVertex(V, Par, E, Tol)` is
//!   [`update_vertex_on_edge`]; `MakeVertex(P)` is [`TopoBuilder::make_vertex`];
//!   `Add(parent, child)` / `Remove(parent, child)` are [`TopoBuilder::add`] /
//!   [`TopoBuilder::remove`].
//! * `anEdge.EmptyCopied()` (`cxx:538`, `:578`, `:705`, `:731`, `:764`) is
//!   `TopoShape::empty_copied`, which in this port copies only the topology;
//!   `BRep_TEdge::EmptyCopy` (`BRep_TEdge.cxx:105-133`) copies the
//!   representations too, so [`empty_copied_with_geometry`] re-registers the
//!   same 3D curve, range and p-curves. Without that, `BRepAdaptor_Curve(NewE)`
//!   (`cxx:559`) would have no curve to evaluate.
//! * `BRep_TEdge` / `BRep_GCurve` of `cxx:872-891` are the registry's p-curve
//!   slots: `IsCurveOnSurface(theSurf, L)` is "the slot of that face", and
//!   `PCurve(PC->Translated(offset))` rewrites the *first* p-curve of that slot
//!   (the forward one of a seam), then `Remove`+`Append` moves the
//!   representation to the end of the list, i.e. the seam's two p-curves swap
//!   order (`break` at `cxx:890`). [`shift_pcurve_period`] reproduces both
//!   effects.
//! * `Emap(CurEdge)` / `EPmap(Key)` (`cxx:416`, `:529`) are `NCollection_DataMap`
//!   `operator()`, which inserts an empty sequence when the key is absent; the
//!   port keeps insertion-ordered `(key, value)` vectors with the same
//!   behaviour.
//! * `TopTools_DataMapOfShapeShape VEmap` (`cxx:657`) is keyed by `Vfirst` and
//!   valued by the seam edge; it is an ordered vector keyed by `IsSame`.
//! * Two reference loops index past their sequence when an `IsSame` / parameter
//!   lookup fails: the face of a seam edge (`cxx:516-523`, `j` stays
//!   `Eseq.Length()+1`) and the vertex of a split piece (`cxx:739-745`, `j`
//!   stays `Seq2.Length()+1`). OCCT reads out of bounds there, so the port
//!   skips that piece instead of inventing a rule.

use std::collections::HashSet;
use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpTrsf2d, GpVec2d};
use occt_core::precision::CONFUSION;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;

use crate::abs::{Orientation, ShapeType};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_builder_update_vertex::update_vertex_on_edge;
use crate::brep_tool::BRepTool;
use crate::brep_tools::is_really_closed;
use crate::brep_tools_substitution::BRepToolsSubstitution;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::tgeometry::{EdgeGeom, GeometryRegistry};
use crate::topo_tools::{
    children_keep_orientation, edge_vertices, edges_of_wire, explorer_shapes, is_same,
    shape_map_key,
};

use super::BRepOffsetApiDraftAngle;

/// `CorrectWires()` (`cxx:272-909`).
impl BRepOffsetApiDraftAngle {
    pub(super) fn correct_wires(&mut self) {
        // `cxx:275`
        let tol_inter = 1.0e-7;
        // `cxx:280`: `BRep_Builder BB` is created once and reused.
        let bb = TopoBuilder::new();

        let Some(shape) = self.my_shape.clone() else {
            return;
        };

        // ---------------------------------------------------------------
        // `cxx:282-308`: the really-closed edges of every wire.
        // ---------------------------------------------------------------
        let mut eseq: Vec<TopoShape> = Vec::new();
        let mut wseq: Vec<TopoShape> = Vec::new();
        let mut fseq: Vec<TopoShape> = Vec::new();
        for cur_face in explorer_shapes(&shape, ShapeType::Face) {
            let Some(face) = Face::wrap(cur_face.clone()) else {
                continue;
            };
            for cur_wire in explorer_shapes(&cur_face, ShapeType::Wire) {
                let mut seen: HashSet<_> = HashSet::new();
                for edge in edges_of_wire(&Wire(cur_wire.clone())) {
                    // `cxx:290-295`: the map keys on `TopTools_ShapeMapHasher`,
                    // so an edge stored twice is visited once.
                    if !seen.insert(shape_map_key(&edge.0)) {
                        continue;
                    }
                    // `cxx:300-306`
                    if is_really_closed(&edge, &face) {
                        eseq.push(edge.0.clone());
                        wseq.push(cur_wire.clone());
                        fseq.push(cur_face.clone());
                    }
                }
            }
        }

        // `cxx:310-319`: the maps the phases below fill.
        let mut emap: Vec<(TopoShape, Vec<f64>)> = Vec::new();
        let mut non_seam: Vec<TopoShape> = Vec::new();
        let mut non_seam_wires: Vec<TopoShape> = Vec::new();
        let mut pars_non_seam: Vec<Vec<f64>> = Vec::new();
        let mut seam: Vec<Vec<TopoShape>> = Vec::new();
        let mut pars_seam: Vec<Vec<f64>> = Vec::new();
        let mut wfmap: Vec<(TopoShape, TopoShape)> = Vec::new();
        let mut wwmap: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();

        // ---------------------------------------------------------------
        // `cxx:320-443`
        // ---------------------------------------------------------------
        for i in 0..eseq.len() {
            let cur_edge = eseq[i].clone();
            let cur_wire = wseq[i].clone();
            let cur_face = fseq[i].clone();
            let (Some(face), Some(edge)) = (
                Face::wrap(cur_face.clone()),
                Edge::wrap(cur_edge.clone()),
            ) else {
                continue;
            };
            // `cxx:330-331`
            let bac2d1 = Curve2dAdaptor::of_edge_face(&edge, &face);
            let bac2d1r = Curve2dAdaptor::of_edge_face(
                &Edge(cur_edge.clone().oriented(Orientation::Reversed)),
                &face,
            );
            // `cxx:333-335`: the surface the 3D intersection point is taken on.
            let Some(surf) = BRepTool::face_surface(&face) else {
                continue;
            };
            // `cxx:337`
            let tol_cur_e = BRepTool::edge_tolerance(&edge);
            for wire in explorer_shapes(&cur_face, ShapeType::Wire) {
                // `cxx:341-343`
                if is_same(&wire, &cur_wire) {
                    continue;
                }
                let mut pts: Vec<GpPnt> = Vec::new();
                let mut wadd = false;
                for an_edge in edges_of_wire(&Wire(wire.clone())) {
                    // `cxx:352-353`
                    let bac2d2 = Curve2dAdaptor::of_edge_face(&an_edge, &face);
                    // `cxx:355-356`
                    let mut ginter = intersect_adaptors(&bac2d1, &bac2d2, tol_inter);
                    if !ginter.is_done() || ginter.result().is_empty() {
                        // `cxx:359-361`: retry with the reversed seam edge.
                        ginter = intersect_adaptors(&bac2d1r, &bac2d2, tol_inter);
                        if !ginter.is_done() || ginter.result().is_empty() {
                            continue;
                        }
                    }
                    // `cxx:366`
                    wadd = true;
                    // `cxx:367-370`
                    if !wfmap.iter().any(|(k, _)| is_same(k, &wire)) {
                        wfmap.push((wire.clone(), cur_face.clone()));
                    }
                    // `cxx:371-390`: the index of `anEdge` inside `NonSeam`.
                    let mut ind = 0usize;
                    for (j, candidate) in non_seam.iter().enumerate() {
                        if is_same(candidate, &an_edge.0) {
                            ind = j + 1;
                            break;
                        }
                    }
                    if ind == 0 {
                        non_seam.push(an_edge.0.clone());
                        non_seam_wires.push(wire.clone());
                        pars_non_seam.push(Vec::new());
                        seam.push(Vec::new());
                        pars_seam.push(Vec::new());
                        ind = non_seam.len();
                    }
                    // `cxx:391-394`
                    map_slot(&mut emap, &cur_edge);
                    // `cxx:397-400`
                    let tol_cmp = tol_cur_e.max(BRepTool::edge_tolerance(&an_edge));
                    // `cxx:403-429`
                    for k in 1..=ginter.nb_points() {
                        let p2d_int = ginter.point(k);
                        let a_p2d: &GpPnt2d = p2d_int.value();
                        let a_p3d = surf.value(a_p2d.x(), a_p2d.y());
                        // `cxx:411-419`: is that 3D point already recorded?
                        let mut ied = 0usize;
                        for (j, p) in pts.iter().enumerate() {
                            if a_p3d.is_equal_with_tolerance(p, tol_cmp) {
                                ied = j + 1;
                                break;
                            }
                        }
                        if ied == 0 {
                            pts.push(a_p3d);
                            map_slot(&mut emap, &cur_edge).push(p2d_int.param_on_first());
                            pars_non_seam[ind - 1].push(p2d_int.param_on_second());
                            seam[ind - 1].push(cur_edge.clone());
                            pars_seam[ind - 1].push(p2d_int.param_on_first());
                        }
                    }
                }
                // `cxx:430-438`
                if wadd {
                    map_slot(&mut wwmap, &cur_wire).push(wire.clone());
                }
            }
        }

        // ---------------------------------------------------------------
        // `cxx:446-464`: sort `ParsNonSeam`, keeping `Seam` / `ParsSeam`
        // aligned. The reference's `j` / `k` are 1-based.
        // ---------------------------------------------------------------
        for i in 0..non_seam.len() {
            let n = pars_non_seam[i].len();
            for j in 1..n {
                for k in (j + 1)..=n {
                    if pars_non_seam[i][k - 1] < pars_non_seam[i][j - 1] {
                        pars_non_seam[i].swap(j - 1, k - 1);
                        seam[i].swap(j - 1, k - 1);
                        pars_seam[i].swap(j - 1, k - 1);
                    }
                }
            }
        }

        // ---------------------------------------------------------------
        // `cxx:466-483`: sort every `Emap` sequence.
        // ---------------------------------------------------------------
        for (_, seq) in emap.iter_mut() {
            let n = seq.len();
            for i in 1..n {
                for j in (i + 1)..=n {
                    if seq[j - 1] < seq[i - 1] {
                        seq.swap(i - 1, j - 1);
                    }
                }
            }
        }

        // ---------------------------------------------------------------
        // `cxx:485-499`: seed the three maps for every `Emap` key.
        // ---------------------------------------------------------------
        let mut epmap: Vec<(TopoShape, Vec<f64>)> = Vec::new();
        let mut evmap: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
        let mut ewmap: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
        for (key, _) in emap.iter() {
            epmap.push((key.clone(), Vec::new()));
            evmap.push((key.clone(), Vec::new()));
            ewmap.push((key.clone(), Vec::new()));
        }

        // `cxx:501-503`
        let mut sub = BRepToolsSubstitution::new();

        // ---------------------------------------------------------------
        // `cxx:504-592`: split the non-seam edges.
        // ---------------------------------------------------------------
        for i in 0..non_seam.len() {
            let Some(an_edge) = Edge::wrap(non_seam[i].clone()) else {
                continue;
            };
            // `cxx:512` indexes `ParsNonSeam(i)(1)`; the sequence always holds
            // at least the point that created the entry (`cxx:416-425`).
            if pars_non_seam[i].is_empty() {
                continue;
            }
            // `cxx:508-511`
            let (Some(vfirst), Some(vlast)) = edge_vertices(&an_edge) else {
                continue;
            };
            let (first_par, last_par) = BRepTool::edge_parameters(&an_edge);
            let mut firstind = 1usize;
            let par = pars_non_seam[i][0];
            let mut seam_edge = seam[i][0].clone();
            // `cxx:515-523`: the face the first seam edge belongs to.
            let Some(j) = eseq.iter().position(|s| is_same(s, &seam_edge)) else {
                // `cxx:516-523` leaves `j == Eseq.Length()+1` and reads
                // `Fseq(j)` out of bounds; the port drops the edge.
                continue;
            };
            let _the_face = fseq[j].clone();
            let mut prevpar = first_par;
            let mut prev_v = vfirst.0.clone();
            // `cxx:526-534`
            if (par - first_par).abs() <= CONFUSION {
                let _ = update_vertex_on_edge(
                    &vfirst.0,
                    pars_seam[i][0],
                    &seam_edge,
                    BRepTool::vertex_tolerance(&vfirst),
                );
                map_slot(&mut epmap, &seam_edge).push(pars_seam[i][0]);
                map_slot(&mut evmap, &seam_edge).push(vfirst.0.clone());
                map_slot(&mut ewmap, &seam_edge).push(non_seam_wires[i].clone());
                firstind = 2;
            }
            // `cxx:537-576`
            let mut new_edges: Vec<TopoShape> = Vec::new();
            let len = pars_non_seam[i].len();
            let mut jj = firstind;
            while jj <= len {
                // `cxx:538-539`
                let mut new_e = empty_copied_with_geometry(&an_edge.0);
                let par = pars_non_seam[i][jj - 1];
                GeometryRegistry::global().set_edge_range(&new_e, prevpar, par);
                seam_edge = seam[i][jj - 1].clone();
                let new_v: TopoShape;
                // `cxx:545-547`
                if jj == len && (par - last_par).abs() <= CONFUSION {
                    new_v = vlast.0.clone();
                    if firstind == 2 && jj == 2 {
                        // `cxx:548-556`
                        let _ = update_vertex_on_edge(
                            &vlast.0,
                            pars_seam[i][jj - 1],
                            &seam_edge,
                            BRepTool::vertex_tolerance(&vlast),
                        );
                        map_slot(&mut epmap, &seam_edge).push(pars_seam[i][jj - 1]);
                        map_slot(&mut evmap, &seam_edge).push(vlast.0.clone());
                        map_slot(&mut ewmap, &seam_edge).push(non_seam_wires[i].clone());
                        break;
                    }
                } else {
                    // `cxx:558-562`: `BRepAdaptor_Curve(NewE).Value(par)` plus
                    // `BRepLib_MakeVertex(Point)`. An edge with no 3D curve has
                    // no `Value` in the port either, so that piece is left out.
                    let Some(point) = BRepTool::edge_curve(&Edge(new_e.clone()))
                        .map(|c| c.d0(par))
                    else {
                        jj += 1;
                        continue;
                    };
                    let vertex = bb.make_vertex(point, 0.0);
                    new_v = vertex.0.clone();
                    let _ = update_vertex_on_edge(&new_v, par, &new_e, 10.0 * CONFUSION);
                }
                // `cxx:563-566`
                let _ = update_vertex_on_edge(
                    &new_v,
                    pars_seam[i][jj - 1],
                    &seam_edge,
                    10.0 * CONFUSION,
                );
                // `cxx:567-569`
                new_e.set_orientation(Orientation::Forward);
                bb.add(&mut new_e, &prev_v.oriented(Orientation::Forward));
                bb.add(&mut new_e, &new_v.oriented(Orientation::Reversed));
                new_edges.push(new_e);
                // `cxx:570-575`
                map_slot(&mut epmap, &seam_edge).push(pars_seam[i][jj - 1]);
                map_slot(&mut evmap, &seam_edge).push(new_v.clone());
                map_slot(&mut ewmap, &seam_edge).push(non_seam_wires[i].clone());
                prevpar = par;
                prev_v = new_v;
                jj += 1;
            }
            // `cxx:577-589`: the last piece, up to `LastPar`.
            let mut new_e = empty_copied_with_geometry(&an_edge.0);
            let par = last_par;
            if (prevpar - par).abs() > CONFUSION {
                GeometryRegistry::global().set_edge_range(&new_e, prevpar, par);
                new_e.set_orientation(Orientation::Forward);
                bb.add(&mut new_e, &prev_v.oriented(Orientation::Forward));
                bb.add(&mut new_e, &vlast.0.oriented(Orientation::Reversed));
                new_edges.push(new_e);
            }
            // `cxx:592`
            sub.substitute(&an_edge.0, &new_edges);
        }

        // ---------------------------------------------------------------
        // `cxx:596-626`: sort `EPmap`, keeping `EVmap` / `EWmap` in step.
        // ---------------------------------------------------------------
        for idx in 0..epmap.len() {
            let n = epmap[idx].1.len();
            for i in 1..n {
                for j in (i + 1)..=n {
                    let (a, b) = (i - 1, j - 1);
                    if epmap[idx].1[b] < epmap[idx].1[a] {
                        epmap[idx].1.swap(a, b);
                        evmap[idx].1.swap(a, b);
                        ewmap[idx].1.swap(a, b);
                    }
                }
            }
        }

        // ---------------------------------------------------------------
        // `cxx:628-656`: drop consecutive duplicates.
        // ---------------------------------------------------------------
        for idx in 0..epmap.len() {
            let mut remove = true;
            while remove {
                remove = false;
                let mut i = 0usize;
                while i + 1 < epmap[idx].1.len() {
                    if (epmap[idx].1[i] - epmap[idx].1[i + 1]).abs() <= CONFUSION {
                        epmap[idx].1.remove(i + 1);
                        evmap[idx].1.remove(i + 1);
                        ewmap[idx].1.remove(i + 1);
                        remove = true;
                    }
                    i += 1;
                }
            }
        }

        // ---------------------------------------------------------------
        // `cxx:658-777`: split the seam edges.
        // ---------------------------------------------------------------
        let mut vemap: Vec<(TopoShape, TopoShape)> = Vec::new();
        for (key, seq) in emap.iter() {
            let Some(an_edge) = Edge::wrap(key.clone()) else {
                continue;
            };
            if seq.is_empty() {
                continue;
            }
            // `cxx:667-671`
            let seq2 = map_get(&epmap, key);
            let seq_ver = map_get(&evmap, key);
            let _seq_wire = map_get(&ewmap, key);
            let mut onepoint = false;
            // `cxx:673-676`
            let (Some(vfirst), Some(vlast)) = edge_vertices(&an_edge) else {
                continue;
            };
            let (first_par, last_par) = BRepTool::edge_parameters(&an_edge);
            // `cxx:678-680`
            let mut fpar = first_par;
            let mut lpar = seq[0];
            let mut firstind = 1usize;
            // `cxx:681-702`
            if (fpar - lpar).abs() <= CONFUSION {
                // `cxx:685` reads `Seq(2)`; a single-point sequence at the
                // start parameter is out of bounds in the reference.
                if seq.len() < 2 {
                    continue;
                }
                firstind = 2;
                fpar = seq[0];
                lpar = seq[1];
            } else if seq.len() % 2 != 0 {
                // `cxx:690-692`
                vemap.push((vfirst.0.clone(), an_edge.0.clone()));
                firstind = 2;
                fpar = seq[0];
                if seq.len() > 2 {
                    lpar = seq[1];
                } else {
                    onepoint = true;
                }
            }
            // `cxx:704-729`
            let mut new_edges: Vec<TopoShape> = Vec::new();
            if !onepoint {
                let mut new_e = empty_copied_with_geometry(&an_edge.0);
                GeometryRegistry::global().set_edge_range(&new_e, fpar, lpar);
                new_e.set_orientation(Orientation::Forward);
                if firstind == 1 {
                    // `cxx:711-716`
                    bb.add(&mut new_e, &vfirst.0.oriented(Orientation::Forward));
                    bb.add(&mut new_e, &seq_ver[0].oriented(Orientation::Reversed));
                } else {
                    // `cxx:718-726`: `SeqVer(2)` is out of bounds in the
                    // reference when the seam collected a single vertex.
                    if seq_ver.len() < 2 {
                        continue;
                    }
                    bb.add(&mut new_e, &seq_ver[0].oriented(Orientation::Forward));
                    bb.add(&mut new_e, &seq_ver[1].oriented(Orientation::Reversed));
                }
                new_edges.push(new_e);
                firstind += 1;
                // `cxx:730-755`
                let mut i = firstind;
                while i < seq.len() {
                    let mut new_e = empty_copied_with_geometry(&an_edge.0);
                    fpar = seq[i - 1];
                    lpar = seq[i];
                    GeometryRegistry::global().set_edge_range(&new_e, fpar, lpar);
                    // `cxx:738-745`: the vertex whose parameter matches `fpar`.
                    let Some(jf) = seq2
                        .iter()
                        .position(|p| (fpar - p).abs() <= CONFUSION)
                        .map(|p| p + 1)
                    else {
                        // The reference reads `SeqVer(j)` with `j` past the
                        // end here; the port leaves the piece out.
                        i += 2;
                        continue;
                    };
                    if jf + 1 > seq_ver.len() {
                        i += 2;
                        continue;
                    }
                    // `cxx:746-751`
                    new_e.set_orientation(Orientation::Forward);
                    bb.add(&mut new_e, &seq_ver[jf - 1].oriented(Orientation::Forward));
                    bb.add(&mut new_e, &seq_ver[jf].oriented(Orientation::Reversed));
                    new_edges.push(new_e);
                    i += 2;
                }
            }
            // `cxx:757-762`
            let last = seq.len();
            fpar = seq[last - 1];
            lpar = last_par;
            if (fpar - lpar).abs() <= CONFUSION {
                // `cxx:760-762`: `continue` skips the substitution entirely.
                continue;
            }
            // `cxx:763-773`
            let mut new_e = empty_copied_with_geometry(&an_edge.0);
            GeometryRegistry::global().set_edge_range(&new_e, fpar, lpar);
            new_e.set_orientation(Orientation::Forward);
            let Some(last_ver) = seq_ver.last() else {
                continue;
            };
            bb.add(&mut new_e, &last_ver.oriented(Orientation::Forward));
            bb.add(&mut new_e, &vlast.0.oriented(Orientation::Reversed));
            new_edges.push(new_e);
            // `cxx:776`
            sub.substitute(&an_edge.0, &new_edges);
        }

        // ---------------------------------------------------------------
        // `cxx:779-815`: drop edges touching a missing seam extremity.
        // ---------------------------------------------------------------
        for (v, e) in vemap.iter() {
            let Some(w) = eseq
                .iter()
                .position(|s| is_same(s, e))
                .map(|i| wseq[i].clone())
            else {
                continue;
            };
            let mut e_to_remove: Option<TopoShape> = None;
            for cur_e in edges_of_wire(&Wire(w.clone())) {
                // `cxx:799-802`
                if is_same(&cur_e.0, e) {
                    continue;
                }
                let (Some(vfirst), Some(vlast)) = edge_vertices(&cur_e) else {
                    continue;
                };
                if is_same(&vfirst.0, v) || is_same(&vlast.0, v) {
                    e_to_remove = Some(cur_e.0.clone());
                    break;
                }
            }
            // `cxx:812-815`
            if let Some(etr) = e_to_remove {
                let mut w_mut = w.clone();
                w_mut.set_free(true);
                bb.remove(&mut w_mut, &etr);
            }
        }

        // ---------------------------------------------------------------
        // `cxx:817-826`
        // ---------------------------------------------------------------
        sub.build(&shape);
        let mut my_shape = shape;
        // `aSub.IsCopied(myShape)` then `myShape = aSub.Copy(myShape).First()`.
        if let Some(list) = sub.copy(&my_shape) {
            if let Some(first) = list.first() {
                my_shape = first.clone();
            }
        }

        // ---------------------------------------------------------------
        // `cxx:828-908`: rebuild the wires.
        // ---------------------------------------------------------------
        for (cur_wire_key, wires) in wwmap.iter() {
            // `cxx:836-838`
            let Some(cur_wire_copy) = sub.copy(cur_wire_key).and_then(|l| l.first()).cloned()
            else {
                continue;
            };
            let mut cur_wire = cur_wire_copy;
            cur_wire.set_free(true);
            for a_wire_key in wires {
                // `cxx:841-843`
                let Some(the_face) = wfmap
                    .iter()
                    .find(|(k, _)| is_same(k, a_wire_key))
                    .map(|(_, v)| v.clone())
                else {
                    continue;
                };
                let Some(a_wire_copy) = sub.copy(a_wire_key).and_then(|l| l.first()).cloned()
                else {
                    continue;
                };
                let a_wire = a_wire_copy;
                // `cxx:847-892`: adjust the period of the p-curves.
                if let Some(face) = Face::wrap(the_face.clone()) {
                    let face_key = GeometryRegistry::shape_key(&face.0);
                    for an_edge in edges_of_wire(&Wire(a_wire.clone())) {
                        // `cxx:851-857`
                        let (pfirst, plast) = BRepTool::uv_points(&an_edge, &face);
                        let bc2d = Curve2dAdaptor::of_edge_face(&an_edge, &face);
                        let pmid = bc2d.d0(0.5 * (bc2d.first_parameter() + bc2d.last_parameter()));
                        let Some(pmid) = pmid else {
                            continue;
                        };
                        let two_pi = 2.0 * std::f64::consts::PI;
                        let mut offset = GpVec2d::zero();
                        let mut translate = false;
                        // `cxx:858-862`
                        if pfirst.x() - two_pi > CONFUSION
                            || plast.x() - two_pi > CONFUSION
                            || pmid.x() - two_pi > CONFUSION
                        {
                            offset = GpVec2d::new(-two_pi, 0.0);
                            translate = true;
                        }
                        // `cxx:863-870`
                        if pfirst.x() < -CONFUSION
                            || plast.x() < -CONFUSION
                            || pmid.x() < -CONFUSION
                        {
                            offset = GpVec2d::new(two_pi, 0.0);
                            translate = true;
                        }
                        // `cxx:871-891`
                        if translate {
                            shift_pcurve_period(&an_edge.0, face_key, &offset);
                        }
                    }
                }
                // `cxx:894-899`: `TopoDS_Iterator(aWire, false)`.
                for edge in children_keep_orientation(&a_wire) {
                    let mut wire_mut = cur_wire.clone();
                    bb.add(&mut wire_mut, &edge);
                    cur_wire = wire_mut;
                }
                // `cxx:901-906`
                let mut the_face_mut = match sub.copy(&the_face).and_then(|l| l.first()) {
                    Some(first) => first.clone(),
                    None => the_face.clone(),
                };
                the_face_mut.set_free(true);
                bb.remove(&mut the_face_mut, &a_wire);
            }
        }

        self.my_shape = Some(my_shape);
    }
}

/// `Geom2dInt_GInter::Perform(AC1, AC2, Tol, Tol)` (`cxx:355`, `:359`) over two
/// `BRepAdaptor_Curve2d` windows.
fn intersect_adaptors(a: &Curve2dAdaptor, b: &Curve2dAdaptor, tol: f64) -> Geom2dIntGInter {
    let mut ginter = Geom2dIntGInter::new();
    let (Some(c1), Some(c2)) = (a.curve().cloned(), b.curve().cloned()) else {
        return ginter;
    };
    let d1 = Geom2dIntGInter::domain_over_range(
        c1.as_ref(),
        a.first_parameter(),
        a.last_parameter(),
        tol,
    );
    let d2 = Geom2dIntGInter::domain_over_range(
        c2.as_ref(),
        b.first_parameter(),
        b.last_parameter(),
        tol,
    );
    ginter.perform(c1.as_ref(), &d1, c2.as_ref(), &d2, tol, tol);
    ginter
}

/// `anEdge.EmptyCopied()` together with the representations that
/// `BRep_TEdge::EmptyCopy` (`BRep_TEdge.cxx:105-133`) keeps: the 3D curve, its
/// range and the p-curve slots. `TopoShape::empty_copied` only copies the
/// topology, so the geometry side-table entry is re-registered here.
fn empty_copied_with_geometry(src: &TopoShape) -> TopoShape {
    let dst = src.empty_copied();
    let reg = GeometryRegistry::global();
    if let Some(g) = reg.edge_geom(src) {
        reg.set_edge(
            &dst,
            EdgeGeom {
                curve: g.curve,
                curve3d_null: g.curve3d_null,
                first: g.first,
                last: g.last,
                tolerance: g.tolerance,
                same_parameter: g.same_parameter,
                same_range: g.same_range,
                degenerated: g.degenerated,
                pcurves: g.pcurves,
                pcurve_ranges: g.pcurve_ranges,
            },
        );
    }
    dst
}

/// `NCollection_DataMap::operator()`: the sequence bound to `key`, inserted as
/// an empty one when absent.
fn map_slot<'a, T: Default>(map: &'a mut Vec<(TopoShape, T)>, key: &TopoShape) -> &'a mut T {
    match map.iter().position(|(k, _)| is_same(k, key)) {
        Some(i) => &mut map[i].1,
        None => {
            map.push((key.clone(), T::default()));
            let last = map.len() - 1;
            &mut map[last].1
        }
    }
}

/// `NCollection_DataMap::Find` on the same ordered map, cloned.
fn map_get<T: Clone + Default>(map: &[(TopoShape, T)], key: &TopoShape) -> T {
    map.iter()
        .find(|(k, _)| is_same(k, key))
        .map(|(_, v)| v.clone())
        .unwrap_or_default()
}

/// `cxx:871-891`: translate the first p-curve of `edge` on the face identified
/// by `face_key`, then move that representation to the end of the list.
///
/// The reference walks `BRep_TEdge::ChangeCurves()` and stops at the first
/// `BRep_GCurve` with `IsCurveOnSurface(theSurf, L)`; for a seam edge that is
/// the forward p-curve, and `Remove` + `Append` on the representation list puts
/// it after the reversed one. The registry stores the two p-curves of a seam
/// in forward-then-reversed order under one face key, so the same result is a
/// rotation of that pair.
fn shift_pcurve_period(edge: &TopoShape, face_key: usize, offset: &GpVec2d) {
    let reg = GeometryRegistry::global();
    let mut curves: Vec<Arc<dyn Curve2d>> = reg.edge_pcurves(edge, face_key);
    let Some(first) = curves.first() else {
        return;
    };
    let mut trsf = GpTrsf2d::identity();
    trsf.set_translation_vec(offset);
    let shifted: Arc<dyn Curve2d> = Arc::from(first.transformed(&trsf));
    curves[0] = shifted;
    if curves.len() >= 2 {
        // `TE->ChangeCurves().Remove(itcr); TE->ChangeCurves().Append(GC)`.
        let last = curves.len() - 1;
        curves.swap(0, last);
    }
    reg.set_edge_pcurves(edge, face_key, curves);
}
