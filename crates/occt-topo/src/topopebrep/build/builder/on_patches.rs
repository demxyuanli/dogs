//! `TopOpeBRepBuild_Builder::FillOnPatches` and the two helpers it needs.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_GridSS.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `AreFacesCoincideInArea(theBaseFace, theFace, theEdge, allEdges, isSameOri)` | `1597-1780` |
//! | `FillOnPatches(anEdgesON, aBaseFace, avoidMap)` | `1784-2036` |
//! | `FindFacesTouchingEdge(aFace, anEdge, aShRank, aFaces)` | `2038-2085` |
//!
//! `FillOnPatches` builds the faces made only of ON matter that
//! `GFillFaceNotSameDomSFS` could not reach: the ON edges are kept when the
//! face of the other operand that supports them coincides with `aBaseFace` in
//! the area the edge points into, and the coincidence test is
//! [`are_faces_coincide_in_area`].
//!
//! UNPORTED: nothing in `FillOnPatches` / `FindFacesTouchingEdge`. In
//! `AreFacesCoincideInArea`:
//! - `GeomAPI_ProjectPointOnSurf(aPnt, aSurf, umin, umax, vmin, vmax)`
//!   (`cxx:1750`) is the unbounded
//!   [`occt_geom::geom_api::project_point_on_surface`]; the port's helper has no
//!   range-limited overload.
//! - the second `aSurf->D1(u, v, aPnt, d1u, d1v)` (`cxx:1755`) writes into the
//!   same `aPnt` OCCT reuses; the port keeps the projected point.

use std::collections::HashMap;

use occt_core::gp::dir2d::GpDir2d;
use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::vec2d::GpVec2d;
use occt_core::precision::INFINITE;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface;
use crate::brep_class_intersector::{BRepClassEdge, BRepClassIntersector};
use crate::brep_tool::BRepTool;
use crate::geom_bnd_lib_elclib2d::line_value_2d;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{is_equal, is_same, shape_map_key, ShapeMapKey};
use crate::topopebrep::build::tools::{fc2d_curve_on_surface, write_edge_pcurve};
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;
use crate::topopebrep::ds::connex::Connexity;
use crate::topopebrep::ds::data_structure::IndexedShapeSet;
use crate::topopebrep::ds::export::fds_data;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::Builder;

/// `AreFacesCoincideInArea` (`cxx:1597-1780`).
///
/// `Some(is_same_ori)` is OCCT's `true`, `None` its `false`.
pub(crate) fn are_faces_coincide_in_area(
    the_base_face: &TopoShape,
    the_face: &TopoShape,
    the_edge: &TopoShape,
    all_edges: &[TopoShape],
) -> Option<bool> {
    let an_edge_ori = the_edge.orientation();
    if an_edge_ori != Orientation::Forward && an_edge_ori != Orientation::Reversed {
        return None;
    }
    let reverse = an_edge_ori == Orientation::Reversed;

    let a_base_face = Face::wrap(the_base_face.clone())?;
    let a_face = Face::wrap(the_face.clone())?;
    let an_edge = Edge::wrap(the_edge.clone())?;

    // `cxx:1626-1636`: the p-curve on the base face, stored or made.
    let (pcref, pf, pl) = match curve_on_surface(&an_edge, &a_base_face)
        .and_then(|_| crate::boptools_2d::curve_on_surface_range(&an_edge, &a_base_face))
    {
        Some(v) => v,
        None => {
            let (c, f, l) = fc2d_curve_on_surface(&an_edge, &a_base_face)?;
            let tol = BRepTool::edge_tolerance(&an_edge);
            write_edge_pcurve(&an_edge, c.clone(), &a_base_face, tol);
            (c, f, l)
        }
    };

    let t = 0.456789;
    let pm = (1.0 - t) * pf + t * pl;
    let (pt, mut d1) = pcref.d1(pm);
    if d1.magnitude() < occt_core::precision::RESOLUTION {
        return None;
    }
    if reverse {
        d1.reverse();
    }
    // `gp_Vec2d vecInside(-d1.Y(), d1.X())` then `gp_Lin2d aLin(pt, vecInside)`.
    let vec_inside = GpVec2d::new(-d1.y(), d1.x());
    let a_dir = GpDir2d::from_vec2d(&vec_inside).ok()?;
    let a_lin = GpLin2d::from_pnt_dir(pt, a_dir);

    // `cxx:1655-1724`: the nearest intersection of the ray with the other
    // boundary edges.
    let mut has_int = false;
    let mut p_lin_min = f64::MAX;
    let tol2d = occt_core::precision::PCONFUSION;
    let mut max_dist = BRepTool::face_tolerance(&a_base_face)
        .max(BRepTool::face_tolerance(&a_face));
    let mut is_error = false;

    for e in all_edges.iter() {
        let a_e = Edge::wrap(e.clone())?;
        let tol_e = BRepTool::edge_tolerance(&a_e);
        if tol_e > max_dist {
            max_dist = tol_e;
        }
        // `cxx:1674-1677`: skip the ray's own p-curve.
        if is_equal(&a_e.0, &an_edge.0)
            || (a_e.0.orientation() != Orientation::Forward
                && a_e.0.orientation() != Orientation::Reversed
                && is_same(&a_e.0, &an_edge.0))
        {
            continue;
        }

        match curve_on_surface(&a_e, &a_base_face) {
            Some(_) => {}
            None => {
                let Some((c, _f, _l)) = fc2d_curve_on_surface(&a_e, &a_base_face) else {
                    is_error = true;
                    break;
                };
                write_edge_pcurve(&a_e, c, &a_base_face, tol_e);
            }
        }

        let bce = BRepClassEdge::new(a_e.clone(), a_base_face.clone());
        let mut an_inter = BRepClassIntersector::new();
        an_inter.perform(&a_lin, p_lin_min, tol2d, &bce);
        if an_inter.is_done() {
            for i in 1..=an_inter.nb_points() {
                let p_lin = an_inter.point(i).map(|p| p.param_on_first()).unwrap_or(0.0);
                if p_lin > tol2d && p_lin < p_lin_min {
                    p_lin_min = p_lin;
                    has_int = true;
                }
            }
            for i in 1..=an_inter.nb_segments() {
                let Some(a_is) = an_inter.segment(i) else {
                    continue;
                };
                let p_lin_f = if a_is.has_first_point() {
                    a_is.first_point().param_on_first()
                } else {
                    -INFINITE
                };
                let p_lin_l = if a_is.has_last_point() {
                    a_is.last_point().param_on_first()
                } else {
                    INFINITE
                };
                if p_lin_f < tol2d && p_lin_l > -tol2d {
                    is_error = true;
                } else if p_lin_f > tol2d && p_lin_f < p_lin_min {
                    p_lin_min = p_lin_f;
                    has_int = true;
                }
            }
        }
    }

    if is_error || !has_int {
        return None;
    }

    // `cxx:1728-1745`: the normal at a point of the area.
    let a_p2d = line_value_2d(p_lin_min * t, &a_lin);
    let a_bas = SurfaceAdaptor::of_face(&a_base_face);
    let (a_pnt, d1u, d1v) = a_bas.d1(a_p2d.x(), a_p2d.y())?;
    let mut a_norm_base = d1u.crossed(&d1v);
    let mut mag = a_norm_base.magnitude();
    if mag < occt_core::precision::RESOLUTION {
        return None;
    }
    if a_base_face.0.orientation() == Orientation::Reversed {
        mag = -mag;
    }
    a_norm_base = a_norm_base.divided(mag);

    // `cxx:1746-1768`: project the point on the second face and compare normals.
    let a_surf = BRepTool::face_surface(&a_face)?;
    let a_proj = occt_geom::geom_api::project_point_on_surface(
        a_surf.as_ref(),
        &a_pnt,
        occt_core::precision::CONFUSION,
    )?;
    if a_proj.distance > max_dist {
        return None;
    }
    let (_, d1u, d1v) = a_surf.d1(a_proj.u, a_proj.v);
    let mut a_norm = d1u.crossed(&d1v);
    let mut mag = a_norm.magnitude();
    if mag < occt_core::precision::RESOLUTION {
        return None;
    }
    if a_face.0.orientation() == Orientation::Reversed {
        mag = -mag;
    }
    a_norm = a_norm.divided(mag);

    // `cxx:1771-1778`.
    let dot = a_norm_base.dot(&a_norm);
    const MIN_DOT: f64 = 0.9999;
    if dot.abs() < MIN_DOT {
        return None;
    }
    Some(dot > 0.0)
}

impl Builder {
    /// `FillOnPatches(anEdgesON, aBaseFace, avoidMap)` (`cxx:1784-2036`).
    pub fn fill_on_patches(        &mut self,
        an_edges_on: &[TopoShape],
        a_base_face: &TopoShape,
        avoid_map: &IndexedShapeSet,
    ) {
        let mut ff = a_base_face.clone();
        ff.set_orientation(Orientation::Forward);

        let rank_bf = self.shape_rank(a_base_face);
        let rank_opp = if rank_bf == 1 {
            2
        } else if rank_bf == 2 {
            1
        } else {
            return;
        };

        let mut wes = WireEdgeSet::new(&ff);
        let mut a_map_on = IndexedShapeSet::new();
        let mut a_map_on1 = IndexedShapeSet::new();
        let mut a_map_f_state: HashMap<ShapeMapKey, i32> = HashMap::new();

        let mut all_edges: Vec<TopoShape> = Vec::new();
        // `cxx:1812-1821`: ON edges that are not part of `aBaseFace`'s own
        // boundary and were not already used.
        for a_e in an_edges_on.iter() {
            let mut a_e_rev = a_e.clone();
            a_e_rev.reverse();
            if !self.my_on_elem_map.contains(a_e)
                && !self.my_on_elem_map.contains(&a_e_rev)
                && !avoid_map.contains(a_e)
            {
                all_edges.push(a_e.clone());
                a_map_on.add(a_e);
            }
        }

        // `cxx:1822-1843`: the wire case adds its children, the edge case
        // collects the ON boundary edges.
        let mut has_wires = false;
        let extent = self.my_on_elem_map.extent();
        for i in 1..=extent {
            let Some(a_e) = self.my_on_elem_map.find_key(i).cloned() else {
                continue;
            };
            if a_e.shape_type() == ShapeType::Wire {
                for child in crate::topo_tools::explorer_shapes(&a_e, ShapeType::Shape) {
                    if !avoid_map.contains(&child) {
                        all_edges.push(child);
                        has_wires = true;
                    }
                }
            } else if !avoid_map.contains(&a_e) {
                all_edges.push(a_e.clone());
                a_map_on1.add(&a_e);
            }
        }

        // `cxx:1846-1922`: ON parts whose opposite ancestor supports a face
        // coinciding with `aBaseFace`.
        let mut an_anc_map: HashMap<ShapeMapKey, TopoShape> = HashMap::new();
        if !a_map_on.is_empty() {
            let keys: std::collections::HashSet<ShapeMapKey> =
                (1..=a_map_on.extent())
                    .filter_map(|i| a_map_on.find_key(i))
                    .map(shape_map_key)
                    .collect();
            self.fill_sec_edge_ancestor_map(rank_opp, &keys, &mut an_anc_map);
        }

        if !an_anc_map.is_empty() {
            let mut cx = Connexity::new();
            let hds_snapshot = self.data_structure().map(|hds| hds.clone());
            if let Some(hds) = hds_snapshot.as_ref() {
                cx.prepare(hds);
                for a_e in an_edges_on.iter() {
                    let Some(an_anc_e) = an_anc_map.get(&shape_map_key(a_e)) else {
                        continue;
                    };
                    if avoid_map.contains(a_e) {
                        continue;
                    }
                    let a_faces = cx.edge_connexity_shape_index(an_anc_e, hds, rank_opp);
                    let mut is_on_face = false;
                    let mut a_faces_to_check: Vec<TopoShape> = Vec::new();
                    for a_f in a_faces.iter() {
                        if let Some(state) = a_map_f_state.get(&shape_map_key(a_f)) {
                            if *state != 0 {
                                is_on_face = true;
                            }
                        } else {
                            a_faces_to_check.push(a_f.clone());
                        }
                        if is_on_face {
                            break;
                        }
                    }
                    for a_f in a_faces_to_check.iter() {
                        if is_on_face {
                            break;
                        }
                        let same_ori = if a_e.orientation() != Orientation::Forward
                            && a_e.orientation() != Orientation::Reversed
                        {
                            let mut fwd = a_e.clone();
                            fwd.set_orientation(Orientation::Forward);
                            let mut rev = a_e.clone();
                            rev.set_orientation(Orientation::Reversed);
                            are_faces_coincide_in_area(a_base_face, a_f, &fwd, &all_edges)
                                .or_else(|| {
                                    are_faces_coincide_in_area(
                                        a_base_face,
                                        a_f,
                                        &rev,
                                        &all_edges,
                                    )
                                })
                        } else {
                            are_faces_coincide_in_area(a_base_face, a_f, a_e, &all_edges)
                        };
                        match same_ori {
                            Some(true) => {
                                a_map_f_state.insert(shape_map_key(a_f), 1);
                                is_on_face = true;
                            }
                            _ => {
                                a_map_f_state.insert(shape_map_key(a_f), 0);
                            }
                        }
                    }
                    if is_on_face {
                        wes.add_start_element(a_e);
                    }
                }
            }
        }

        // `cxx:1924-2018`: the same test for the ON matter that is part of
        // `aBaseFace`'s own boundary.
        let mut an_anc_map: HashMap<ShapeMapKey, TopoShape> = HashMap::new();
        if !a_map_on1.is_empty() {
            let keys: std::collections::HashSet<ShapeMapKey> =
                (1..=a_map_on1.extent())
                    .filter_map(|i| a_map_on1.find_key(i))
                    .map(shape_map_key)
                    .collect();
            self.fill_sec_edge_ancestor_map(rank_bf, &keys, &mut an_anc_map);
        }

        if has_wires || !an_anc_map.is_empty() {
            let cx = Connexity::new();
            let hds_snapshot = self.data_structure().map(|hds| hds.clone());
            let extent = self.my_on_elem_map.extent();
            for i in 1..=extent {
                let Some(a_e) = self.my_on_elem_map.find_key(i).cloned() else {
                    continue;
                };
                let mut an_edge: Option<TopoShape> = None;
                let mut an_anc_e: Option<TopoShape> = None;
                if a_e.shape_type() == ShapeType::Wire {
                    // `cxx:1937-1951`: one non-degenerated edge of the wire.
                    for child in crate::topo_tools::explorer_shapes(&a_e, ShapeType::Edge) {
                        if avoid_map.contains(&child) {
                            break;
                        }
                        if !BRepTool::is_degenerated(&Edge(child.clone())) {
                            an_edge = Some(child.clone());
                            an_anc_e = Some(child);
                        }
                        if an_anc_e.is_some() {
                            break;
                        }
                    }
                } else if let Some(anc) = an_anc_map.get(&shape_map_key(&a_e)) {
                    if !avoid_map.contains(&a_e) {
                        an_edge = Some(a_e.clone());
                        an_anc_e = Some(anc.clone());
                    }
                }

                let (Some(an_edge), Some(an_anc_e)) = (an_edge, an_anc_e) else {
                    continue;
                };

                // `cxx:1960-1971`: the faces of the other operand touching it.
                let Some(hds) = hds_snapshot.as_ref() else {
                    break;
                };
                let mut a_face_cx: Vec<TopoShape> = Vec::new();
                cx.face_edge_connex_faces(a_base_face, &an_anc_e, hds, &mut a_face_cx);
                if a_face_cx.is_empty() {
                    continue;
                }
                let a_cnx_f = a_face_cx[0].clone();
                a_face_cx.clear();
                let mut a_faces: Vec<TopoShape> = Vec::new();
                self.find_faces_touching_edge(&a_cnx_f, &an_anc_e, rank_opp, &mut a_faces);

                let mut is_on_face = false;
                let mut a_faces_to_check: Vec<TopoShape> = Vec::new();
                for a_f in a_faces.iter() {
                    if let Some(state) = a_map_f_state.get(&shape_map_key(a_f)) {
                        if *state != 0 {
                            is_on_face = true;
                        }
                    } else {
                        a_faces_to_check.push(a_f.clone());
                    }
                    if is_on_face {
                        break;
                    }
                }
                for a_f in a_faces_to_check.iter() {
                    if is_on_face {
                        break;
                    }
                    match are_faces_coincide_in_area(a_base_face, a_f, &an_edge, &all_edges) {
                        Some(true) => {
                            a_map_f_state.insert(shape_map_key(a_f), 1);
                            is_on_face = true;
                        }
                        _ => {
                            a_map_f_state.insert(shape_map_key(a_f), 0);
                        }
                    }
                }
                if is_on_face {
                    if a_e.shape_type() == ShapeType::Wire {
                        wes.add_shape(&a_e);
                    } else {
                        wes.add_start_element(&a_e);
                    }
                }
            }
        }

        // `cxx:2021-2033`: build the patches and keep them as ON faces.
        wes.init_shapes();
        wes.init_start_elements();
        if wes.more_shapes() || wes.more_start_elements() {
            let mut lof: Vec<TopoShape> = Vec::new();
            self.gwes_make_faces(&ff, &mut wes, &mut lof);
            for a_f in lof.iter() {
                if let Some(a_face) = Face::wrap(a_f.clone()) {
                    self.my_on_faces_map.add(&a_face.0, a_base_face.clone());
                }
            }
        }
    }

    /// `FindFacesTouchingEdge(aFace, anEdge, aShRank, aFaces)` (`cxx:2038-2085`).
    pub fn find_faces_touching_edge(
        &self,
        a_face: &TopoShape,
        an_edge: &TopoShape,
        a_sh_rank: i32,
        a_faces: &mut Vec<TopoShape>,
    ) {
        let Some(hds) = self.data_structure() else {
            return;
        };
        let bds = hds.ds();
        let an_edge_ind = bds.shape_index(an_edge, false);
        if an_edge_ind == 0 {
            return;
        }

        let li = bds.shape_interferences(a_face, false);
        for i_handle in li.iter() {
            let Ok(i) = i_handle.read() else {
                continue;
            };
            if !i.is_shape_shape() {
                continue;
            }
            let (gt, gi, st, si) = fds_data(&i);
            if gt != Kind::Edge || st != Kind::Face {
                continue;
            }
            if gi != an_edge_ind {
                continue;
            }
            let tfe = i.transition();
            if tfe.shape_before() != ShapeType::Face || tfe.shape_after() != ShapeType::Face {
                continue;
            }
            let Some(fs) = bds.shape(si, false) else {
                continue;
            };
            if self.shape_rank(&fs) != a_sh_rank {
                continue;
            }
            a_faces.push(fs);
        }
    }
}

