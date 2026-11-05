//! `BRepFill_OffsetWire::PerformWithBiLo` (`cxx:680-1184`): the body of the
//! offset, driven by the bisecting locus.
//!
//! The phases, in order:
//!
//! * the offset curves and circles of the spine (`MakeCircle` / `MakeOffset`),
//!   stored in `myMap` keyed by the spine vertex / edge they come from
//!   (`cxx:772-797`),
//! * the vertices the bisectors put on those parallels, collected in `MapBis`
//!   and `MapVerPar` (`cxx:829-1024`),
//! * the trimming of every parallel by those vertices (`TrimEdge`,
//!   `cxx:1035-1085`), the confusions `MapVV` reports (`cxx:1086-1113`),
//! * `MakeWires` (`cxx:1116`), the lift of every vertex onto the reference
//!   plane (`cxx:1122-1131`), `BRepLib::BuildCurves3d` with the endpoint
//!   reconciliation (`cxx:1133-1176`) and `FixHoles` (`cxx:1178-1181`).
//!
//! Deviations (handle / exception correspondence only, no missing branch):
//! * a null `mySpine` is `Option::is_none`, a null `TopoDS_Shape` member is
//!   `my_shape: Option<TopoShape>`,
//! * `Standard_Failure` raised by `Trim.IntersectWith` / `AddOrConfuse`
//!   (`cxx:917-957`) answers the state `Perform` leaves after its `catch`
//!   (`myShape.Nullify(); myIsDone = false;`),
//! * `NCollection_DataMap` / `NCollection_Sequence` / `NCollection_IndexedDataMap`
//!   are `HashMap` / `Vec`,
//! * `MAT_Node::Distance()` of a null node handle (OCCT dereferences it) is
//!   `f64::MAX` here, which keeps the `MapBis` bind of `cxx:1015-1023` open,
//! * `BRepFill_TrimEdgeTool` takes the bisector curve instead of a
//!   `Bisector_Bisec` (`brep_fill_trim_edge_tool/mod.rs:33-39`), so the
//!   `GeomBis` value is cloned into an `Arc<dyn Curve2d>`.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::{GpPln, GpPnt};
use occt_core::kernel::geomabs::JoinType;
use occt_core::precision::CONFUSION;
use occt_geom::plane::GeomPlane;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_fill_offset_wire::curves::{make_circle, make_offset};
use crate::brep_fill_offset_wire::maps::ShapeToShapeMap;
use crate::brep_fill_offset_wire::perform::{compute_offset_shape, is_altitude_significant};
use crate::brep_fill_offset_wire::trim_edge::trim_edge;
use crate::brep_fill_offset_wire::vertices::{
    set_vertex_point, update_vertex_tolerance, update_vertex_point, vertex_from_node,
};
use crate::brep_fill_offset_wire::BRepFillOffsetWire;
use crate::brep_fill_trim_edge_tool::{Geometry2d, TrimEdgeTool};
use crate::brep_lib_build_curves3d::build_curves3d;
use crate::brep_mat2d_bisecting_locus::BRepMat2dBisectingLocus;
use crate::brep_mat2d_link_topo_bilo::BRepMat2dLinkTopoBilo;
use crate::brep_tool::BRepTool;
use crate::brep_tools_substitution::BRepToolsSubstitution;
use crate::builder::TopoBuilder;
use crate::mat::NodeId;
use crate::mat2d::Mat2dItem;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{
    edge_vertices, edges_of, is_same, shape_map_key, vertices_of, wire_vertices, wires_of_face,
    ShapeMapKey,
};

impl BRepFillOffsetWire {
    /// `PerformWithBiLo(Spine, Offset, Locus, Link, Join, Alt)` (`cxx:680-1184`).
    ///
    /// `Link` is taken by shared reference: the port's
    /// `BRepMat2dLinkTopoBilo::GeneratingShape` only reads the link
    /// (`brep_mat2d_link_topo_bilo.rs:135`).
    pub fn perform_with_bilo(
        &mut self,
        spine: &Face,
        offset: f64,
        locus: &BRepMat2dBisectingLocus,
        link: &BRepMat2dLinkTopoBilo,
        join: JoinType,
        alt: f64,
    ) {
        // `cxx:688-693`.
        self.my_is_done = false;
        let forward = spine.0.oriented(Orientation::Forward);
        self.my_work_spine = Some(Face(forward.clone()));
        self.my_join_type = join;
        self.my_offset = offset;
        self.my_shape = None;

        // `cxx:695-701`.
        if self.my_spine.is_none() {
            self.my_spine = Some(Face(forward));
        }
        self.my_map.clear();

        // `cxx:704-708`: the zero offset is the spine itself.
        if offset.abs() < CONFUSION {
            let my_spine = self.my_spine.clone().expect("mySpine");
            let shape = compute_offset_shape(&my_spine, &mut self.my_map, alt);
            self.my_shape = Some(shape);
            self.my_is_done = true;
            return;
        }

        // `cxx:713-716`.
        if self.kpart_circle(offset, alt) {
            return;
        }

        let work_spine = self.my_work_spine.clone().expect("myWorkSpine");

        // `cxx:744-752`: the reference plane, translated by `Alt` along its
        // normal.
        let Some(ref_plane) = BRepTool::face_surface(&work_spine)
            .and_then(|surface| reference_plane(&surface, &work_spine, alt))
        else {
            self.my_shape = None;
            self.my_is_done = false;
            return;
        };

        // `cxx:754-764`: the two ends of the spine wire for an open result.
        let mut ends: [Option<Vertex>; 2] = [None, None];
        if self.my_is_open_result {
            let my_spine = self.my_spine.clone().expect("mySpine");
            if let Some(wire) = wires_of_face(&my_spine).into_iter().next() {
                let (v1, v2) = wire_vertices(&wire);
                ends = [v1, v2];
            }
        }

        // `cxx:767-770`.
        if locus.number_of_contours() == 0 {
            return;
        }

        // `cxx:772-797`: one parallel per spine edge, one circle per spine
        // vertex.
        for ic in 1..=locus.number_of_contours() {
            let mut pe = link.generating_shape(locus.basic_elt(ic, locus.number_of_elts(ic)));
            for ie in 1..=locus.number_of_elts(ic) {
                let se = link.generating_shape(locus.basic_elt(ic, ie));
                if se.shape_type() == ShapeType::Vertex {
                    // `cxx:779-784`: the ends of an open spine get no circle.
                    let is_end = ends
                        .iter()
                        .flatten()
                        .any(|v| is_same(&v.0, &se));
                    if !is_end {
                        make_circle(
                            &Edge(pe.clone()),
                            &Vertex(se.clone()),
                            &work_spine,
                            self.my_offset,
                            &mut self.my_map,
                            &ref_plane,
                        );
                    }
                } else {
                    // `cxx:786-793`.
                    make_offset(
                        &Edge(se.clone()),
                        &work_spine,
                        self.my_offset,
                        &mut self.my_map,
                        &ref_plane,
                        self.my_is_open_result,
                        self.my_join_type,
                        &ends,
                    );
                    pe = se.clone();
                }
            }
        }

        // `cxx:799-812`: hanging arcs on the ends of an open result.
        if self.my_is_open_result && self.my_join_type == JoinType::Arc {
            if !self.my_map.is_empty()
                && self
                    .my_map
                    .find_key(1)
                    .map_or(false, |k| k.shape_type() == ShapeType::Vertex)
            {
                self.my_map.remove_from_index(1);
            }
            let extent = self.my_map.extent();
            if !self.my_map.is_empty()
                && self
                    .my_map
                    .find_key(extent)
                    .map_or(false, |k| k.shape_type() == ShapeType::Vertex)
            {
                self.my_map.remove_last();
            }
        }

        // `cxx:818-1024`: the vertices the bisectors put on the parallels.
        let mut detromp: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        let mut map_node_vertex: HashMap<NodeId, TopoShape> = HashMap::new();
        let mut map_bis: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        let mut map_ver_par: HashMap<ShapeMapKey, Vec<f64>> = HashMap::new();

        let graph = locus.graph();
        for i in 1..=graph.number_of_arcs() {
            // `cxx:834-835`.
            let arc_id = graph.arc(i);
            let arc = graph.arc_by_id(arc_id);
            let mut reverse = false;
            let bisec = locus.geom_bis(arc, &mut reverse);

            // `cxx:840-841`.
            let (Some(first_elt), Some(second_elt)) = (arc.first_element(), arc.second_element())
            else {
                continue;
            };
            let s = [
                link.generating_shape(first_elt),
                link.generating_shape(second_elt),
            ];

            // `cxx:855-861`: the parallels of both ends of the arc.
            if !(self.my_map.contains(&s[0]) && self.my_map.contains(&s[1])) {
                continue;
            }
            let e0 = Edge(self.my_map.find_from_key(&s[0]).expect("FindFromKey")[0].clone());
            let e1 = Edge(self.my_map.find_from_key(&s[1]).expect("FindFromKey")[0].clone());
            let e = [e0, e1];

            // `cxx:870-883`: the nodes that limit the arc, in the direction the
            // bisector curve is stored.
            let (node1, node2) = if reverse {
                (arc.second_node(), arc.first_node())
            } else {
                (arc.first_node(), arc.second_node())
            };

            // `cxx:887-902`: the arc extremities that are exactly on the offset.
            let mut start_on_edge = false;
            let mut end_on_edge = false;
            let mut vs: Option<Vertex> = None;
            let mut ve: Option<Vertex> = None;
            if let Some(n1) = node1 {
                let node = graph.node_by_id(n1);
                if !node.infinite() {
                    let p = locus.geom_elt_of_node(n1);
                    vs = vertex_from_node(
                        n1,
                        node.distance(),
                        node.infinite(),
                        self.my_offset,
                        &p,
                        &mut map_node_vertex,
                    );
                    start_on_edge = vs.is_some();
                }
            }
            if let Some(n2) = node2 {
                let node = graph.node_by_id(n2);
                if !node.infinite() {
                    let p = locus.geom_elt_of_node(n2);
                    ve = vertex_from_node(
                        n2,
                        node.distance(),
                        node.infinite(),
                        self.my_offset,
                        &p,
                        &mut map_node_vertex,
                    );
                    end_on_edge = ve.is_some();
                }
            }
            // `cxx:900-903`.
            if self.my_join_type == JoinType::Intersection {
                start_on_edge = false;
                end_on_edge = false;
            }

            // `cxx:909-913`.
            let g1 = geometry2d(locus.geom_elt(graph.basic_elt_by_id(first_elt)));
            let g2 = geometry2d(locus.geom_elt(graph.basic_elt_by_id(second_elt)));
            let bis_value = bisec.value().expect("Bisector_Bisec::Value");
            let bis_curve: Arc<dyn Curve2d> = Arc::new(bis_value.clone());
            let trim = TrimEdgeTool::new(bis_curve.clone(), &g1, &g2, self.my_offset);

            // `cxx:917-919`.
            let Ok(mut params) = trim.intersect_with(
                &e[0],
                &e[1],
                &s[0],
                &s[1],
                ends[0].as_ref(),
                ends[1].as_ref(),
                self.my_join_type,
                self.my_is_open_result,
            ) else {
                self.my_shape = None;
                self.my_is_done = false;
                return;
            };

            // `cxx:921-930`: the vertices of the intersections, on the plane
            // `Z = 0`.
            let builder = TopoBuilder::new();
            let mut vertices: Vec<Vertex> = Vec::new();
            for p in &params {
                let p2 = bis_curve.d0(p.x());
                vertices.push(builder.make_vertex(GpPnt::new(p2.x(), p2.y(), 0.0), CONFUSION));
            }

            // `cxx:931-944`.
            if start_on_edge {
                if trim
                    .add_or_confuse(true, &e[0], &e[1], &mut params)
                    .is_err()
                {
                    self.my_shape = None;
                    self.my_is_done = false;
                    return;
                }
                if params.len() == vertices.len() && !params.is_empty() {
                    vertices[0] = vs.clone().expect("VS");
                } else {
                    vertices.insert(0, vs.clone().expect("VS"));
                }
            }
            // `cxx:945-958`.
            if end_on_edge {
                if trim
                    .add_or_confuse(false, &e[0], &e[1], &mut params)
                    .is_err()
                {
                    self.my_shape = None;
                    self.my_is_done = false;
                    return;
                }
                if params.len() == vertices.len() && !params.is_empty() {
                    let last = params.len() - 1;
                    vertices[last] = ve.clone().expect("VE");
                } else {
                    vertices.push(ve.clone().expect("VE"));
                }
            }

            // `cxx:965-977`.
            detromp.entry(shape_map_key(&s[0])).or_default();
            detromp.entry(shape_map_key(&s[1])).or_default();
            self.update_detromp(
                &mut detromp,
                &s[0],
                &s[1],
                &vertices,
                &params,
                &bisec,
                start_on_edge,
                end_on_edge,
                &trim,
            );

            // `cxx:982-1024`.
            if !vertices.is_empty() && params.len() == vertices.len() {
                for k in 0..2 {
                    let key = shape_map_key(&e[k].0);
                    let bis = map_bis.entry(key.clone()).or_default();
                    let par = map_ver_par.entry(key).or_default();
                    for (ii, v) in vertices.iter().enumerate() {
                        bis.push(v.0.clone());
                        par.push(if k == 0 { params[ii].y() } else { params[ii].z() });
                    }
                }
            } else {
                // `cxx:1013-1024`: a complete circle's parallel can sit wholly
                // in the zone between the spine and the offset.
                let d1 = node1.map_or(f64::MAX, |n| graph.node_by_id(n).distance());
                let d2 = node2.map_or(f64::MAX, |n| graph.node_by_id(n).distance());
                for k in 0..2 {
                    let key = shape_map_key(&e[k].0);
                    if !map_bis.contains_key(&key) && d1 > self.my_offset && d2 > self.my_offset {
                        map_bis.insert(key.clone(), Vec::new());
                        map_ver_par.insert(key, Vec::new());
                    }
                }
            }
        }

        // `cxx:1035-1085`: trim every parallel by the vertices collected on it.
        let mut map_vv = ShapeToShapeMap::new();
        let my_spine = self.my_spine.clone().expect("mySpine");
        let extent = self.my_map.extent();
        for j in 1..=extent {
            let current_spine = self.my_map.find_key(j).expect("FindKey").clone();
            let Some(first) = self.my_map.value(j).and_then(|l| l.first()).cloned() else {
                continue;
            };
            let current_edge = Edge(first);
            if let Some(slot) = self.my_map.value_mut(j) {
                slot.clear();
            }

            let key = shape_map_key(&current_edge.0);
            let Some(bis_shapes) = map_bis.get(&key).cloned() else {
                continue;
            };
            if bis_shapes.is_empty() {
                // `cxx:1077-1083`: a complete circle keeps its parallel whole.
                if let Some(slot) = self.my_map.value_mut(j) {
                    slot.push(current_edge.0.clone());
                }
                continue;
            }

            // `cxx:1050-1061`.
            let mut ind_of_e = 0i32;
            if self.my_is_open_result {
                if j == 1 {
                    ind_of_e = 1;
                } else if j == extent {
                    ind_of_e = -1;
                }
            }

            let mut the_ver: Vec<Vertex> = bis_shapes.iter().map(|s| Vertex(s.clone())).collect();
            let mut the_par: Vec<f64> = map_ver_par.get(&key).cloned().unwrap_or_default();
            let detromp_list = detromp
                .get(&shape_map_key(&current_spine))
                .cloned()
                .unwrap_or_default();
            let mut a_seq: Vec<TopoShape> = Vec::new();
            trim_edge(
                &current_edge,
                &current_spine,
                &my_spine,
                &detromp_list,
                &mut the_ver,
                &mut the_par,
                &mut a_seq,
                &mut map_vv,
                ind_of_e,
            );
            if let Some(slot) = self.my_map.value_mut(j) {
                for s in a_seq {
                    slot.push(s);
                }
            }
        }

        // `cxx:1086-1113`: the confused vertices `TrimEdge` reported are
        // replaced by the surviving one everywhere in `myMap`.
        for ind in 1..=map_vv.extent() {
            let old_v = Vertex(map_vv.find_key(ind).expect("FindKey").clone());
            let new_v = Vertex(map_vv.value(ind).expect("MapVV").clone());
            let p1 = BRepTool::vertex_point(&old_v);
            let p2 = BRepTool::vertex_point(&new_v);
            update_vertex_tolerance(&new_v, p1.distance(&p2));

            let mut subst = BRepToolsSubstitution::new();
            subst.substitute(&old_v.0, &[new_v.0.oriented(Orientation::Forward)]);

            let extent = self.my_map.extent();
            for j in 1..=extent {
                let Some(list) = self.my_map.value_mut(j) else {
                    continue;
                };
                for item in list.iter_mut() {
                    subst.build(item);
                    if subst.is_copied(item) {
                        let saved = item.orientation();
                        let first = subst.copy(item).expect("IsCopied")[0].clone();
                        *item = first.oriented(saved);
                    }
                }
            }
        }

        // `cxx:1116`.
        self.make_wires();

        // `cxx:1122-1131`: every vertex was built on the plane `Z = 0`; lift it
        // onto the reference plane.
        let mut map_vertex: Vec<ShapeMapKey> = Vec::new();
        if let Some(shape) = self.my_shape.clone() {
            for v in vertices_of(&shape) {
                let key = shape_map_key(&v.0);
                if !map_vertex.contains(&key) {
                    map_vertex.push(key);
                    let p = BRepTool::vertex_point(&v);
                    let q = ref_plane.d0(p.x(), p.y());
                    update_vertex_point(&v, &q, CONFUSION);
                }
            }
        }

        // `cxx:1133-1134`.
        let Some(shape) = self.my_shape.clone() else {
            self.my_is_done = true;
            return;
        };
        build_curves3d(&shape);

        // `cxx:1135-1176`: reconcile each vertex with the 3D curve ends.
        map_vertex.clear();
        for an_edge in edges_of(&shape) {
            let (v1, v2) = edge_vertices(&an_edge);
            let (Some(v1), Some(v2)) = (v1, v2) else {
                continue;
            };
            let (f, l) = BRepTool::edge_parameters(&an_edge);
            let Some(curve) = BRepTool::edge_curve_world(&an_edge) else {
                continue;
            };
            let f3d = curve.d0(f);
            let l3d = curve.d0(l);
            let dist1 = f3d.distance(&BRepTool::vertex_point(&v1));
            let dist2 = l3d.distance(&BRepTool::vertex_point(&v2));

            for (v, p, dist) in [(v1, f3d, dist1), (v2, l3d, dist2)] {
                let key = shape_map_key(&v.0);
                if !map_vertex.contains(&key) {
                    set_vertex_point(&v, &p);
                    map_vertex.push(key);
                } else {
                    update_vertex_tolerance(&v, 1.5 * dist);
                }
            }
        }

        // `cxx:1178-1181`.
        if !self.my_is_open_result {
            self.fix_holes();
        }

        // `cxx:1183`.
        self.my_is_done = true;
    }
}

/// `BRep_Tool::Surface` down-cast to `Geom_Plane` plus
/// `RefPlane->Translated(anAlt * gp_Vec(RefPlane->Axis().Direction()))`
/// (`cxx:744-752`): the reference plane of the result, at altitude `Alt` above
/// the work spine.
///
/// `None` when the spine surface is not a plane, which is the `Geom_Plane`
/// down-cast failure OCCT reaches with a null `RefPlane`.
fn reference_plane(surface: &Arc<dyn Surface>, work_spine: &Face, alt: f64) -> Option<Arc<dyn Surface>> {
    let pln = surface.gp_pln()?;
    if !is_altitude_significant(alt) {
        return Some(surface.clone());
    }
    // `cxx:746-750`: the altitude follows the work spine's orientation.
    let mut an_alt = alt;
    if work_spine.0.orientation() == Orientation::Reversed {
        an_alt = -alt;
    }
    let axis = pln.axis();
    let dir = axis.direction();
    let mut pos = pln.position();
    let loc = pos.location();
    pos.set_location(GpPnt::new(
        loc.x() + an_alt * dir.x(),
        loc.y() + an_alt * dir.y(),
        loc.z() + an_alt * dir.z(),
    ));
    Some(Arc::new(GeomPlane::new(GpPln::new(pos))))
}

/// `Geom2d_Geometry` of a `Mat2dItem` as `BRepFill_TrimEdgeTool` takes it.
fn geometry2d(item: &Mat2dItem) -> Geometry2d {
    match item {
        Mat2dItem::Curve(c) => Geometry2d::Curve(Arc::new(c.clone())),
        Mat2dItem::Point(p) => Geometry2d::Point(*p),
    }
}
