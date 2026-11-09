//! `LocOpe_WiresOnShape` (`LocOpe_WiresOnShape.hxx` / `.cxx` / `.lxx`).
//!
//! Binds splitting wires / edges onto the faces (or edges) of a shape, and
//! projects their p-curves so `LocOpe_SplitShape` can split the shape.
//!
//! The class is a plain state holder: `myMapEF` (edge -> face) records the
//! bindings, `myMap` (shape -> shape) records what each wire vertex was
//! projected onto, and the iterator walks `myMapEF`. The heavy lifting lives in
//! the file-local helpers, each ported in its own submodule:
//!
//! * [`projection`] - the three `Project` overloads (`cxx:491-718`).
//! * [`put_pcurve`] - `PutPCurve` (`cxx:722-910`).
//! * [`put_pcurves`] - `PutPCurves` (`cxx:913-1288`).
//! * [`internal_intersections`] - `FindInternalIntersections` (`cxx:1291-1511`).
//! * [`add_edges`] - `Add` (`cxx:1515-1624`).

mod add_edges;
mod internal_intersections;
mod projection;
mod put_pcurve;
mod put_pcurves;

use occt_core::gp::GpPnt2d;
use occt_core::precision::{INFINITE, PCONFUSION};

use crate::abs::{Orientation, ShapeType};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_tool::BRepTool;
use crate::brep_tools::compare_vertices;
use crate::shape::{Compound, Edge, Face, TopoShape, Vertex, Wire};
use crate::topo_tools::{edges_of, is_same, vertices_of};

use super::shape_maps::{ShapeDataMap, ShapeDoubleMap, ShapeIndexedMap, ShapeSet};

/// `LocOpe_WiresOnShape` (`hxx:37-119`).
#[derive(Debug, Clone)]
pub struct WiresOnShape {
    /// `myShape` (`hxx:109`).
    my_shape: TopoShape,
    /// `myMapEF` (`hxx:110`): `IndexedDataMap<edge, face>`, so the iterator can
    /// walk it by 1-based index.
    my_map_ef: ShapeIndexedMap<TopoShape>,
    /// `myFacesWithSection` (`hxx:111`).
    my_faces_with_section: ShapeSet,
    /// `myCheckInterior` (`hxx:113`).
    my_check_interior: bool,
    /// `myMap` (`hxx:114`): wire vertex -> projected vertex, or wire vertex ->
    /// edge when the projection fell inside an edge, and edge -> face for the
    /// edges found by `BindAll` (`cxx:353-359`).
    my_map: ShapeDataMap,
    /// `myDone` (`hxx:115`).
    my_done: bool,
    /// `myIndex` (`hxx:116`); `-1` before `InitEdgeIterator`, as in the
    /// reference's constructor (`cxx:99-104`).
    my_index: isize,
}

impl WiresOnShape {
    /// `LocOpe_WiresOnShape(const TopoDS_Shape&)` (`cxx:97-104`).
    pub fn new(shape: TopoShape) -> Self {
        Self {
            my_shape: shape,
            my_map_ef: ShapeIndexedMap::new(),
            my_faces_with_section: ShapeSet::new(),
            my_check_interior: true,
            my_map: ShapeDataMap::new(),
            my_done: false,
            my_index: -1,
        }
    }

    /// `Init(S)` (`cxx:108-115`).
    pub fn init(&mut self, shape: TopoShape) {
        self.my_shape = shape;
        self.my_check_interior = true;
        self.my_done = false;
        // `cxx:113-114`: `myMap.Clear(); myMapEF.Clear();` - the face set is left
        // alone by the reference too.
        self.my_map.clear();
        self.my_map_ef.clear();
    }

    /// `SetCheckInterior(ToCheckInterior)` (`lxx:19-22`).
    pub fn set_check_interior(&mut self, to_check_interior: bool) {
        self.my_check_interior = to_check_interior;
    }

    /// `IsDone()` (`lxx:26-29`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `IsFaceWithSection(aFace)` (`lxx:33-36`).
    pub fn is_face_with_section(&self, face: &TopoShape) -> bool {
        self.my_faces_with_section.contains(face)
    }

    /// `Bind(const TopoDS_Wire&, const TopoDS_Face&)` (`cxx:118-124`).
    pub fn bind_wire(&mut self, w: &Wire, f: &Face) -> Result<(), String> {
        for e in edges_of(&w.0.clone().oriented(Orientation::Forward)) {
            self.bind_edge(&e, f)?;
        }
        Ok(())
    }

    /// `Bind(const TopoDS_Compound&, const TopoDS_Face&)` (`cxx:127-136`).
    pub fn bind_compound(&mut self, comp: &Compound, f: &Face) -> Result<(), String> {
        for e in edges_of(&comp.0.clone().oriented(Orientation::Forward)) {
            self.bind_edge(&e, f)?;
        }
        // `cxx:135`.
        self.my_faces_with_section.add(&f.0);
        Ok(())
    }

    /// `Bind(const TopoDS_Edge&, const TopoDS_Face&)` (`cxx:140-163`).
    pub fn bind_edge(&mut self, e: &Edge, f: &Face) -> Result<(), String> {
        if !self.my_map_ef.contains(&e.0) {
            // `cxx:146-156`: the edge is only recorded when it is *not* already
            // an edge of the face.
            let already = edges_of(&f.0)
                .iter()
                .any(|cur| is_same(&cur.0, &e.0));
            if !already {
                self.my_map_ef.add(&e.0, f.0.clone());
            }
            Ok(())
        } else {
            // `cxx:161-162`: `throw Standard_ConstructionError()`.
            Err("LocOpe_WiresOnShape::Bind: Standard_ConstructionError".to_string())
        }
    }

    /// `Bind(const TopoDS_Edge& Ewir, const TopoDS_Edge& Efac)` (`cxx:168-171`).
    pub fn bind_edge_to_edge(&mut self, e_wir: &Edge, e_fac: &Edge) {
        self.my_map.bind(&e_wir.0, &e_fac.0);
    }

    /// `BindAll()` (`cxx:175-363`).
    ///
    /// Answers `Err` where `PutPCurves` throws `Standard_ConstructionError`
    /// (`cxx:932`, `cxx:1105`); OCCT lets the exception escape `BindAll` too.
    pub fn bind_all(&mut self) -> Result<(), String> {
        if self.my_done {
            return Ok(());
        }
        let mut the_map = ShapeSet::new();
        let mut map_v = ShapeDataMap::new();

        // `cxx:186-219`: for each bound (edge, face)-pair, project its p-curves
        // and decide, per wire vertex, whether it lands on a vertex of the image
        // edge or on the whole edge.
        let bound_pairs: Vec<(TopoShape, TopoShape)> = self.my_map.iter().cloned().collect();
        for (k, v) in &bound_pairs {
            let eref = Edge(k.clone());
            let eimg = Edge(v.clone());

            put_pcurves::put_pcurves(&eref, &eimg, &self.my_shape)?;

            for vtx in vertices_of(&eref.0) {
                if the_map.contains(&vtx.0) {
                    continue; // `cxx:198`: already processed
                }
                let mut matched = false;
                for vtx2 in vertices_of(&eimg.0) {
                    if is_same(&vtx2.0, &vtx.0) {
                        matched = true;
                        break;
                    } else if compare_vertices(&vtx, &vtx2) {
                        map_v.bind(&vtx.0, &vtx2.0);
                        matched = true;
                        break;
                    }
                }
                if !matched {
                    // `cxx:213-216`.
                    map_v.bind(&vtx.0, &eimg.0);
                }
                the_map.add(&vtx.0);
            }
        }

        // `cxx:221-225`.
        let vertex_bindings: Vec<(TopoShape, TopoShape)> = map_v.iter().cloned().collect();
        for (k, v) in &vertex_bindings {
            self.my_map.bind(k, v);
        }

        // `cxx:227-268`: interior intersections, and the split edges they add.
        let mut splits: ShapeIndexedMap<Vec<TopoShape>> = ShapeIndexedMap::new();
        let mut overlapped_edges = ShapeSet::new();
        for ind in 1..=self.my_map_ef.extent() {
            let edg = Edge(self.my_map_ef.find_key(ind).clone());
            let fac = Face(self.my_map_ef.find_from_index(ind).clone());

            put_pcurve::put_pcurve(&edg, &fac);
            // `cxx:238-242`.
            if crate::boptools_2d::curve_on_surface(&edg, &fac).is_none() {
                continue;
            }

            if self.my_check_interior {
                let mut is_overlapped = false;
                internal_intersections::find_internal_intersections(
                    &edg,
                    &fac,
                    &mut splits,
                    &mut is_overlapped,
                );
                if is_overlapped {
                    overlapped_edges.add(&edg.0);
                }
            }
        }

        // `cxx:254-269`.
        for ind in 1..=splits.extent() {
            let an_edge = splits.find_key(ind).clone();
            if overlapped_edges.contains(&an_edge) {
                continue;
            }
            let a_face = self.my_map_ef.find_from_key(&an_edge).cloned();
            self.my_map_ef.remove_key(&an_edge);
            let new_edges = splits.find_from_index(ind).clone();
            if let Some(a_face) = a_face {
                for e in new_edges {
                    self.my_map_ef.add(&e, a_face.clone());
                }
            }
        }

        // `cxx:271-351`: project every wire vertex of every bound edge onto the
        // face and record where it landed.
        let mut a_vert_param = ShapeDoubleMap::new();
        for ind in 1..=self.my_map_ef.extent() {
            let edg = Edge(self.my_map_ef.find_key(ind).clone());
            let fac = Face(self.my_map_ef.find_from_index(ind).clone());

            for vtx in vertices_of(&edg.0) {
                if the_map.contains(&vtx.0) {
                    continue;
                }
                // `cxx:289-296`.
                let vtx_param = BRepTool::parameter_of_vertex(&edg.0, &vtx.0);
                let ba_curve2d = Curve2dAdaptor::of_edge_face(&edg, &fac);
                let p2d = match ba_curve2d.curve() {
                    Some(_) => ba_curve2d
                        .d0(vtx_param)
                        .unwrap_or(GpPnt2d::new(INFINITE, INFINITE)),
                    None => GpPnt2d::new(INFINITE, INFINITE),
                };

                let mut epro: Option<Edge> = None;
                let mut prm = INFINITE;
                let is_projected = self.my_map.is_bound(&vtx.0);
                // `cxx:300-316`.
                if is_projected {
                    let Some(a_sh) = self.my_map.find(&vtx.0).cloned() else {
                        continue;
                    };
                    if a_sh.shape_type() != ShapeType::Edge {
                        continue;
                    }
                    epro = Some(Edge(a_sh));
                    if let Some(p) = a_vert_param.find(&vtx.0) {
                        prm = p;
                    }
                }

                let ok = projection::project_vertex_on_face(
                    &vtx, &p2d, &fac, &mut epro, &mut prm,
                );
                if !ok || is_projected {
                    continue;
                }
                let epro = epro.expect("Project answers true only with an edge");

                // `cxx:320-347`.
                let mut bound = false;
                for vtx2 in vertices_of(&epro.0) {
                    if is_same(&vtx2.0, &vtx.0) {
                        self.my_map.bind(&vtx.0, &vtx2.0);
                        the_map.add(&vtx.0);
                        bound = true;
                        break;
                    } else if compare_vertices(&vtx, &vtx2) {
                        // `cxx:330-331`: `BRep_Tool::Range(Epro, fac, aF1, aL1)`.
                        if let Some((_pc, a_f1, a_l1)) =
                            crate::boptools_2d::curve_on_surface_range(&epro, &fac)
                        {
                            // `cxx:332-334`.
                            if !BRepTool::is_degenerated(&epro)
                                && ((prm - a_f1).abs() <= PCONFUSION
                                    || (prm - a_l1).abs() <= PCONFUSION)
                            {
                                self.my_map.bind(&vtx.0, &vtx2.0);
                                the_map.add(&vtx.0);
                                bound = true;
                                break;
                            }
                        }
                    }
                }
                if !bound {
                    // `cxx:344-348`.
                    self.my_map.bind(&vtx.0, &epro.0);
                    a_vert_param.bind(&vtx.0, prm);
                }
            }
        }

        // `cxx:353-359`: every edge entry of `myMap` becomes a bound edge too.
        let edge_bindings: Vec<(TopoShape, TopoShape)> = self.my_map.iter().cloned().collect();
        for (k, v) in &edge_bindings {
            if k.shape_type() == ShapeType::Edge {
                self.my_map_ef.add(k, v.clone());
            }
        }

        self.my_done = true;
        Ok(())
    }

    /// `InitEdgeIterator()` (`cxx:366-371`).
    pub fn init_edge_iterator(&mut self) -> Result<(), String> {
        self.bind_all()?;
        self.my_index = 1;
        Ok(())
    }

    /// `MoreEdge()` (`cxx:376-380`).
    pub fn more_edge(&self) -> bool {
        self.my_index <= self.my_map_ef.extent() as isize
    }

    /// The 1-based pair the reference reads through `myMapEF.FindKey(myIndex)`
    /// and `myMapEF(myIndex)`; a caller outside `1..=Extent()` is the reference's
    /// `NCollection_IndexedDataMap::OutOfRange` (`cxx:387`, `cxx:396`).
    fn current(&self) -> (&TopoShape, &TopoShape) {
        assert!(
            self.my_index >= 1 && self.my_index <= self.my_map_ef.extent() as isize,
            "NCollection_IndexedDataMap::OutOfRange (myIndex = {})",
            self.my_index
        );
        let i = self.my_index as usize;
        (self.my_map_ef.find_key(i), self.my_map_ef.find_from_index(i))
    }

    /// `Edge()` (`cxx:385-388`).
    pub fn edge(&self) -> Edge {
        Edge(self.current().0.clone())
    }

    /// `OnFace()` (`cxx:393-397`).
    pub fn on_face(&self) -> Face {
        Face(self.current().1.clone())
    }

    /// `OnEdge(TopoDS_Edge& E)` (`cxx:402-411`): `Some` where the reference
    /// answers `true`.
    pub fn on_edge(&self) -> Option<Edge> {
        let key = self.current().0;
        self.my_map.find(key).map(|e| Edge(e.clone()))
    }

    /// `NextEdge()` (`cxx:415-418`).
    pub fn next_edge(&mut self) {
        self.my_index += 1;
    }

    /// `OnVertex(Vw, Vs)` (`cxx:423-435`).
    pub fn on_vertex(&self, v_wire: &Vertex) -> Option<Vertex> {
        let v = self.my_map.find(&v_wire.0)?;
        if v.shape_type() == ShapeType::Vertex {
            Some(Vertex(v.clone()))
        } else {
            None
        }
    }

    /// `OnEdge(V, Ed, prm)` (`cxx:438-447`).
    pub fn on_edge_of_vertex(&self, v: &Vertex) -> Option<(Edge, f64)> {
        let bound = self.my_map.find(&v.0)?;
        if bound.shape_type() == ShapeType::Vertex {
            return None;
        }
        let ed = Edge(bound.clone());
        let prm = projection::project_vertex_on_edge(v, &ed);
        Some((ed, prm))
    }

    /// `OnEdge(V, EdgeFrom, Ed, prm)` (`cxx:451-487`).
    pub fn on_edge_of_vertex_from(&self, v: &Vertex, edge_from: &Edge) -> Option<(Edge, f64)> {
        let bound = self.my_map.find(&v.0)?;
        if bound.shape_type() == ShapeType::Vertex {
            return None;
        }
        let ed = Edge(bound.clone());
        if !self.my_map_ef.contains(&edge_from.0) {
            return None;
        }
        let a_shape = self.my_map_ef.find_from_key(&edge_from.0).cloned();

        // `cxx:469-486`.
        let prm = if BRepTool::edge_curve_world(&ed).is_none()
            && a_shape
                .as_ref()
                .map(|s| s.shape_type() == ShapeType::Face)
                .unwrap_or(false)
        {
            let a_face = Face(a_shape.expect("checked above"));
            let vtx_param = BRepTool::parameter_of_vertex(&edge_from.0, &v.0);
            let ba_curve2d = Curve2dAdaptor::of_edge_face(edge_from, &a_face);
            // `cxx:478`: `BAcurve2d.Value(vtx_param)` - the reference does not
            // check the p-curve here.
            let p2d = ba_curve2d
                .d0(vtx_param)
                .unwrap_or(GpPnt2d::new(INFINITE, INFINITE));
            projection::project_vertex_2d_on_edge(&p2d, &ed, &a_face)
        } else {
            projection::project_vertex_on_edge(v, &ed)
        };
        Some((ed, prm))
    }
}
