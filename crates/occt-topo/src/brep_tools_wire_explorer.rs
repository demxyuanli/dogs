//! `BRepTools_WireExplorer` — walks the edges of a wire in connection order.
//!
//! Source: `ModelingData/TKBRep/BRepTools/BRepTools_WireExplorer.hxx:36-110`
//! and `BRepTools_WireExplorer.cxx:57-872` (OCCT 8.0.0). The static helpers
//! `SelectDouble` (`cxx:758-775`), `SelectDegenerated` (`cxx:777-795`) and
//! `GetNextParamOnPC` (`cxx:797-872`) are the private functions below.
//!
//! Translation notes:
//! - `TopoDS_Edge` / `TopoDS_Vertex` null handles become `Option`; `More()`
//!   is [`WireExplorer::more`].
//! - `NCollection_DataMap<Shape, List<Shape>>` (`myMap`, `cxx:101`) becomes a
//!   `HashMap` keyed by [`ShapeMapKey`] (TShape + location, i.e. `IsSame`
//!   without the orientation, which is what `TopTools_ShapeMapHasher`
//!   compares), and a `Vec` with `NCollection_List` order: append at the end,
//!   `First()` / `Remove(it)` by index.
//! - `NCollection_IndexedMap vmap` (`cxx:223`) only answers "is the wire
//!   open" and "which kept vertex is FORWARD oriented", so its add/remove
//!   dance (`cxx:252-269`) becomes [`toggle_map_key`] on a `Vec` that keeps
//!   `IndexedMap`'s insertion order.
//! - `anInfEmap` (`cxx:225`) is an `NCollection_Map` in OCCT, i.e. hashed
//!   order; the port keeps insertion order. The two differ only when a closed
//!   wire carries more than one infinite edge, which OCCT itself leaves
//!   undefined (the first match wins).
//! - `NCollection_Map` (`myDoubles`, `cxx:104`; `emap`, `cxx:293`) becomes a
//!   `HashSet<ShapeMapKey>`.
//! - `TopExp::Vertices(E, V1, V2, CumOri)` / `FirstVertex` / `LastVertex`
//!   (`TopExp.cxx:183-251`) are [`edge_vertices_cum_ori`] /
//!   [`first_vertex_cum_ori`] / [`last_vertex_cum_ori`], built on
//!   [`cumulated_children`], the port of `TopoDS_Iterator(S, CumOri)`.
//! - `W.IsNull()` (`cxx:98`, `cxx:135`) has no port counterpart (`TopoShape`
//!   is never null); an empty wire reaches the same `return` through the empty
//!   `vmap` path.
//! - `TopoDS::Edge` / `TopoDS::Vertex` down-casts raise `Standard_TypeMismatch`
//!   in OCCT; the port skips a child of the wrong type, which can only happen
//!   on a wire that violates the shape contract.
//! - `throw Standard_NoSuchObject` in `Orientation()` (`cxx:734`) becomes the
//!   `Err` arm of [`WireExplorer::orientation`].
//! - `myReverse` (`cxx:145`, `cxx:220`) is written and never read in OCCT; the
//!   field is kept for the member layout.
//! - `myTolU` / `myTolV` stay `0.` when the face has no registered surface;
//!   OCCT would dereference the null `Handle(Geom_Surface)` in
//!   `GeomAdaptor_Surface` instead.
//! - `Next()` detaches the connection-vertex list from `myMap` for the length
//!   of the selection and writes it back (`cxx:408` returns a reference into
//!   the map; Rust cannot alias `&mut self` with it). The list content and
//!   order are unchanged.

use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::{CONFUSION, INFINITE, PCONFUSION, REAL_SMALL};
use occt_geom::approx_same_parameter::{u_resolution, v_resolution};
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::iterator::cumulated_children;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::topo_tools::{is_same, shape_map_key, vertices_of, ShapeMapKey};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `TopoDS_Iterator(S, false).Value()` sequence: the stored children.
fn stored_children(s: &TopoShape) -> Vec<TopoShape> {
    s.tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .clone()
}

/// `TopExp::Vertices(E, V1, V2, CumOri = true)` (`TopExp.cxx:214-251`): the
/// first child stored `FORWARD` and the last stored `REVERSED`, with the edge
/// orientation folded in.
pub(crate) fn edge_vertices_cum_ori(edge: &Edge) -> (Option<Vertex>, Option<Vertex>) {
    let mut first = None;
    let mut last = None;
    for child in cumulated_children(&edge.0) {
        if child.shape_type() != ShapeType::Vertex {
            continue;
        }
        match child.orientation() {
            Orientation::Forward => first = Some(Vertex(child)),
            Orientation::Reversed => last = Some(Vertex(child)),
            _ => {}
        }
    }
    (first, last)
}

/// `TopExp::FirstVertex(E, CumOri = true)` (`TopExp.cxx:183-194`).
pub(crate) fn first_vertex_cum_ori(edge: &Edge) -> Option<Vertex> {
    edge_vertices_cum_ori(edge).0
}

/// `TopExp::LastVertex(E, CumOri = true)` (`TopExp.cxx:198-209`).
pub(crate) fn last_vertex_cum_ori(edge: &Edge) -> Option<Vertex> {
    edge_vertices_cum_ori(edge).1
}

/// `NCollection_IndexedMap::Add` + `RemoveKey` (`cxx:252-258`): a key already
/// in the map is removed, a missing key is appended. `IndexedMap::Add` answers
/// the existing index for a bound key (`currsize >= ind`), which is what makes
/// OCCT remove it.
fn toggle_map_key(map: &mut Vec<TopoShape>, v: &TopoShape) {
    let key = shape_map_key(v);
    match map.iter().position(|k| shape_map_key(k) == key) {
        Some(pos) => {
            map.remove(pos);
        }
        None => map.push(v.clone()),
    }
}

/// `SelectDouble` (`cxx:758-775`).
fn select_double(doubles: &HashSet<ShapeMapKey>, list: &mut Vec<TopoShape>, e: &mut Edge) -> bool {
    let mut i = 0;
    while i < list.len() {
        let candidate = &list[i];
        if doubles.contains(&shape_map_key(candidate)) && !is_same(candidate, &e.0) {
            *e = Edge(candidate.clone());
            list.remove(i);
            return true;
        }
        i += 1;
    }
    false
}

/// `SelectDegenerated` (`cxx:777-795`). `E` is overwritten while the list is
/// scanned, exactly like the C++ loop, so the `IsSame` test sees the previous
/// candidate.
fn select_degenerated(list: &mut Vec<TopoShape>, e: &mut Edge) -> bool {
    let mut i = 0;
    while i < list.len() {
        if !is_same(&list[i], &e.0) {
            *e = Edge(list[i].clone());
            if BRepTool::is_degenerated(e) {
                list.remove(i);
                return true;
            }
        }
        i += 1;
    }
    false
}

/// `GetNextParamOnPC` (`cxx:797-872`): the parameter of the first p-curve point
/// outside the `tol_u` / `tol_v` box around `p_ref`, walking in `1/1000` of the
/// range steps.
fn next_param_on_pc(
    pc: &dyn Curve2d,
    p_ref: &GpPnt2d,
    first: f64,
    last: f64,
    tol_u: f64,
    tol_v: f64,
    reverse: bool,
) -> f64 {
    let mut result = if reverse { first } else { last };
    let step = (last - first).abs() / 1000.;
    if reverse {
        let mut start = first;
        let mut found = false;
        while !found && start < last {
            start += step;
            let pnt = pc.d0(start);
            if (p_ref.x() - pnt.x()).abs() < tol_u && (p_ref.y() - pnt.y()).abs() < tol_v {
                continue;
            }
            result = start;
            found = true;
        }
        if !found {
            result = last;
        }
        if result > last {
            result = last;
        }
    } else {
        let mut start = last;
        let mut found = false;
        while !found && start > first {
            start -= step;
            let pnt = pc.d0(start);
            if (p_ref.x() - pnt.x()).abs() < tol_u && (p_ref.y() - pnt.y()).abs() < tol_v {
                continue;
            }
            result = start;
            found = true;
        }
        if !found {
            result = first;
        }
        if result < first {
            result = first;
        }
    }
    result
}

/// `BRepTools_WireExplorer` (`hxx:36-110`).
pub struct WireExplorer {
    /// `myMap` (`cxx:101`): vertex occurrence -> the edges leaving it, in wire
    /// order.
    map: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myEdge` (`cxx:102`).
    edge: Option<Edge>,
    /// `myVertex` (`cxx:103`).
    vertex: Option<Vertex>,
    /// `myFace` (`cxx:104`).
    face: Option<Face>,
    /// `myDoubles` (`cxx:105`).
    doubles: HashSet<ShapeMapKey>,
    /// `myReverse` (`cxx:106`): written at `cxx:145` / `cxx:220`, never read by
    /// OCCT.
    #[allow(dead_code)]
    reverse: bool,
    /// `myTolU` (`cxx:107`).
    tol_u: f64,
    /// `myTolV` (`cxx:108`).
    tol_v: f64,
}

impl Default for WireExplorer {
    fn default() -> Self {
        Self::new()
    }
}

impl WireExplorer {
    /// `BRepTools_WireExplorer()` (`cxx:57-63`).
    pub fn new() -> Self {
        Self {
            map: HashMap::new(),
            edge: None,
            vertex: None,
            face: None,
            doubles: HashSet::new(),
            reverse: false,
            tol_u: 0.0,
            tol_v: 0.0,
        }
    }

    /// `BRepTools_WireExplorer(W)` (`cxx:67-72`).
    pub fn of_wire(w: &Wire) -> Self {
        let mut explorer = Self::new();
        explorer.init(w, None);
        explorer
    }

    /// `BRepTools_WireExplorer(W, F)` (`cxx:76-79`).
    pub fn of_wire_on_face(w: &Wire, f: &Face) -> Self {
        let mut explorer = Self::new();
        explorer.init(w, Some(f));
        explorer
    }

    /// `Init(W)` / `Init(W, F)` (`cxx:83-116`).
    pub fn init(&mut self, w: &Wire, f: Option<&Face>) {
        let (mut umin, mut umax, mut vmin, mut vmax) = (0.0, 0.0, 0.0, 0.0);
        if let Some(face) = f {
            // For the faces based on Cone, BSpline and Bezier compute the UV
            // bounds to precise the UV tolerance values (`cxx:106-113`).
            let kind = SurfaceAdaptor::of_face(face).get_type();
            if kind == GeomAbsSurfaceType::Cone
                || kind == GeomAbsSurfaceType::BSplineSurface
                || kind == GeomAbsSurfaceType::BezierSurface
            {
                let bounds = BRepTool::uv_bounds(face);
                umin = bounds.0;
                umax = bounds.1;
                vmin = bounds.2;
                vmax = bounds.3;
            }
        }
        self.init_with_bounds(w, f, umin, umax, vmin, vmax);
    }

    /// `Init(W, F, UMin, UMax, VMin, VMax)` (`cxx:122-382`).
    pub fn init_with_bounds(
        &mut self,
        w: &Wire,
        f: Option<&Face>,
        umin: f64,
        umax: f64,
        vmin: f64,
        vmax: f64,
    ) {
        self.edge = None;
        self.vertex = None;
        self.map.clear();
        self.doubles.clear();

        self.face = f.cloned();
        let mut df_vert_toler = 0.;
        self.reverse = false;

        if let Some(face) = f {
            let gas = SurfaceAdaptor::of_face(face);
            for v in vertices_of(&w.0) {
                df_vert_toler = BRepTool::vertex_tolerance(&v).max(df_vert_toler);
            }
            if df_vert_toler < CONFUSION {
                // Use tolerance of edges (`cxx:153-159`): `TopoDS_Iterator it(W)`
                // keeps the default `cumOri = true` (`TopoDS_Iterator.hxx:50`),
                // so a REVERSED wire flips the edge orientations it yields.
                for child in cumulated_children(&w.0) {
                    if child.shape_type() != ShapeType::Edge {
                        continue;
                    }
                    df_vert_toler = BRepTool::edge_tolerance(&Edge(child)).max(df_vert_toler);
                }
                if df_vert_toler < CONFUSION {
                    // empty wire (`cxx:161-164`).
                    return;
                }
            }

            let kind = gas.get_type();
            if let Some(surf) = gas.surface() {
                self.tol_u = 2. * u_resolution(surf.as_ref(), df_vert_toler);
                self.tol_v = 2. * v_resolution(surf.as_ref(), df_vert_toler);
            }

            // uresolution for cone with infinite vmin vmax is too small
            // (`cxx:171-201`).
            if kind == GeomAbsSurfaceType::Cone {
                let maxtol = 0.0005 * (umax - umin);
                let tol1 = match gas.d1(umin, vmin) {
                    Some((_, d1u, _)) => {
                        let a = d1u.magnitude();
                        if a <= CONFUSION {
                            maxtol
                        } else {
                            maxtol.min(df_vert_toler / a)
                        }
                    }
                    None => maxtol,
                };
                let tol2 = match gas.d1(umin, vmax) {
                    Some((_, d1u, _)) => {
                        let a = d1u.magnitude();
                        if a <= CONFUSION {
                            maxtol
                        } else {
                            maxtol.min(df_vert_toler / a)
                        }
                    }
                    None => maxtol,
                };
                self.tol_u = 2. * tol1.max(tol2);
            }

            if kind == GeomAbsSurfaceType::BSplineSurface
                || kind == GeomAbsSurfaceType::BezierSurface
            {
                let mut max_tol = self.tol_u.max(self.tol_v);
                if let Some((_, du, dv)) = gas.d1((umax - umin) / 2., (vmax - vmin) / 2.) {
                    let md = (du.x() * du.x()
                        + du.y() * du.y()
                        + du.z() * du.z()
                        + dv.x() * dv.x()
                        + dv.y() * dv.y()
                        + dv.z() * dv.z())
                    .sqrt();
                    if md > REAL_SMALL {
                        if md * max_tol / df_vert_toler < 1.5 {
                            max_tol = 1.5 * df_vert_toler / md;
                        }
                        self.tol_u = max_tol;
                        self.tol_v = max_tol;
                    }
                }
            }

            self.reverse = face.0.orientation() == Orientation::Reversed;
        }

        // Map of vertices to know if the wire is open (`cxx:222-223`) and the
        // map of infinite edges (`cxx:225`).
        let mut vmap: Vec<TopoShape> = Vec::new();
        let mut an_inf_emap: Vec<Edge> = Vec::new();
        let mut v1: Option<Vertex> = None;

        // `TopoDS_Iterator it(W)` (`cxx:232`): `cumOri` default true, so the
        // edge views below carry the composed wire orientation.
        let children = cumulated_children(&w.0);
        for child in &children {
            if child.shape_type() != ShapeType::Edge {
                continue;
            }
            let edge = Edge(child.clone());
            let e_ori = edge.0.orientation();
            if e_ori == Orientation::Internal || e_ori == Orientation::External {
                continue;
            }
            let (a_v1, a_v2) = edge_vertices_cum_ori(&edge);

            if let Some(v) = &a_v1 {
                self.map
                    .entry(shape_map_key(&v.0))
                    .or_default()
                    .push(edge.0.clone());
                // add or remove in the vertex map (`cxx:252-258`).
                let mut forward = v.0.clone();
                forward.set_orientation(Orientation::Forward);
                toggle_map_key(&mut vmap, &forward);
            }
            if let Some(v) = &a_v2 {
                let mut reversed = v.0.clone();
                reversed.set_orientation(Orientation::Reversed);
                toggle_map_key(&mut vmap, &reversed);
            }

            if a_v1.is_none() || a_v2.is_none() {
                let (a_f, a_l) = BRepTool::edge_parameters(&edge);
                if e_ori == Orientation::Forward {
                    if a_f == -INFINITE {
                        an_inf_emap.push(edge.clone());
                    }
                } else if a_l == INFINITE {
                    an_inf_emap.push(edge.clone());
                }
            }
        }

        // Construction of the set of double edges (`cxx:294-304`).
        let mut emap: HashSet<ShapeMapKey> = HashSet::new();
        for child in &children {
            let key = shape_map_key(child);
            if !emap.insert(key.clone()) {
                self.doubles.insert(key);
            }
        }

        // If vmap is not empty the wire is open, let us find the first vertex
        // (`cxx:306-325`).
        if !vmap.is_empty() {
            for key in &vmap {
                if key.orientation() == Orientation::Forward {
                    v1 = Some(Vertex(key.clone()));
                    break;
                }
            }
        } else {
            // The wire is infinite: try to find the first vertex, it may be
            // NULL (`cxx:327-348`).
            for edge in &an_inf_emap {
                let an_ori = edge.0.orientation();
                let (a_f, a_l) = BRepTool::edge_parameters(edge);
                if (an_ori == Orientation::Forward && a_f == -INFINITE)
                    || (an_ori == Orientation::Reversed && a_l == INFINITE)
                {
                    self.edge = Some(edge.clone());
                    self.vertex = None;
                    return;
                }
            }

            // Use the first vertex in iterator order (`cxx:351-365`).
            for child in &children {
                if child.shape_type() != ShapeType::Edge {
                    continue;
                }
                let edge = Edge(child.clone());
                let e_ori = edge.0.orientation();
                if e_ori == Orientation::Internal || e_ori == Orientation::External {
                    continue;
                }
                v1 = edge_vertices_cum_ori(&edge).0;
                break;
            }
        }

        let Some(v1) = v1 else {
            return;
        };
        let Some(list) = self.map.get_mut(&shape_map_key(&v1.0)) else {
            return;
        };
        if list.is_empty() {
            return;
        }
        let edge = Edge(list.remove(0));
        self.edge = Some(edge.clone());
        self.vertex = first_vertex_cum_ori(&edge);
    }

    /// `More()` (`cxx:385-388`).
    pub fn more(&self) -> bool {
        self.edge.is_some()
    }

    /// `Next()` (`cxx:392-706`).
    pub fn next(&mut self) {
        let Some(edge) = self.edge.clone() else {
            return;
        };
        let Some(my_vertex) = last_vertex_cum_ori(&edge) else {
            self.vertex = None;
            self.edge = None;
            return;
        };
        self.vertex = Some(my_vertex.clone());
        let key = shape_map_key(&my_vertex.0);
        let Some(mut list) = self.map.remove(&key) else {
            self.edge = None;
            return;
        };
        let next_edge = self.next_from_list(&edge, &my_vertex, &mut list);
        self.map.insert(key, list);
        self.edge = next_edge;
    }

    /// The body of `Next()` once the connection-vertex list is known
    /// (`cxx:408-705`). `None` is OCCT's null `myEdge`.
    fn next_from_list(
        &mut self,
        edge: &Edge,
        my_vertex: &Vertex,
        list: &mut Vec<TopoShape>,
    ) -> Option<Edge> {
        if list.is_empty() {
            return None;
        }
        if list.len() == 1 {
            // OCC325 (`cxx:416-477`).
            let next_edge = Edge(list[0].clone());
            let (a_v1, a_v2) = edge_vertices_cum_ori(&next_edge);
            let Some(a_v1) = a_v1 else {
                return None;
            };
            if !is_same(&a_v1.0, &my_vertex.0) {
                return None;
            }
            if let Some(face) = self.face.clone() {
                if let Some(a_v2) = &a_v2 {
                    if is_same(&a_v1.0, &a_v2.0) {
                        let prev = curve_on_surface_range(edge, &face);
                        let next = curve_on_surface_range(&next_edge, &face);
                        let (Some((a_prev_pc, p11, p12)), Some((a_next_pc, p21, p22))) =
                            (prev, next)
                        else {
                            return None;
                        };
                        let prev_par = if edge.0.orientation() == Orientation::Forward {
                            p12
                        } else {
                            p11
                        };
                        let (next_f_par, next_l_par) =
                            if next_edge.0.orientation() == Orientation::Forward {
                                (p21, p22)
                            } else {
                                (p22, p21)
                            };
                        let p_prev = a_prev_pc.d0(prev_par);
                        let p_next_f = a_next_pc.d0(next_f_par);
                        let p_next_l = a_next_pc.d0(next_l_par);
                        if p_prev.square_distance(&p_next_f) > p_prev.square_distance(&p_next_l) {
                            return None;
                        }
                    }
                }
            }
            list.clear();
            return Some(next_edge);
        }

        if self.face.is_none() {
            // Without a face, try to return edges as logically as possible
            // (`cxx:483-521`).
            let mut candidate = edge.clone();
            if select_degenerated(list, &mut candidate) {
                return Some(candidate);
            }
            let mut candidate = edge.clone();
            if select_double(&self.doubles, list, &mut candidate) {
                return Some(candidate);
            }
            for i in 0..list.len() {
                if !is_same(&list[i], &edge.0) {
                    return Some(Edge(list.remove(i)));
                }
            }
            return None;
        }

        self.select_by_pcurve_angle(edge, list)
    }

    /// The `else` arm of `Next()` (`cxx:521-703`): more than one edge leaves the
    /// connection vertex and a face is known, so the next edge is the one whose
    /// p-curve leaves the vertex with the smallest angle to the current edge's
    /// p-curve direction. Closed (degenerated) edges are considered in a first
    /// pass.
    fn select_by_pcurve_angle(
        &mut self,
        edge: &Edge,
        list: &mut Vec<TopoShape>,
    ) -> Option<Edge> {
        let face = self.face.clone()?;
        let (a_pcurve, df_f_par, df_l_par) = curve_on_surface_range(edge, &face)?;
        // Get the 2D point equal to <myVertex> in 2D for the current edge
        // (`cxx:539-548`).
        let p_ref = if edge.0.orientation() == Orientation::Reversed {
            a_pcurve.d0(df_f_par)
        } else {
            a_pcurve.d0(df_l_par)
        };

        // Get the next 2D point from the current edge's p-curve (`cxx:549-556`).
        let is_reverse = edge.0.orientation() == Orientation::Reversed;
        let df_m_par = next_param_on_pc(
            a_pcurve.as_ref(),
            &p_ref,
            df_f_par,
            df_l_par,
            self.tol_u,
            self.tol_v,
            is_reverse,
        );
        let p_ref_m = a_pcurve.d0(df_m_par);
        let e_ref_dir = GpVec2d::new(p_ref_m.x() - p_ref.x(), p_ref_m.y() - p_ref.y());
        if e_ref_dir.coord.square_modulus() < REAL_SMALL {
            return None;
        }

        let mut k = 1usize;
        let mut k_min = 0usize;
        let mut dmin = f64::MAX;
        let mut df_min_angle = 3. * PI;
        let mut df_cur_angle = 3. * PI;
        let mut is_degenerated = true;

        for _pass in 0..2 {
            for candidate in list.iter() {
                let e = Edge(candidate.clone());
                if is_same(candidate, &edge.0) {
                    k += 1;
                    continue;
                }
                let (a_vert1, a_vert2) = edge_vertices_cum_ori(&e);
                let (Some(a_vert1), Some(a_vert2)) = (a_vert1, a_vert2) else {
                    k += 1;
                    continue;
                };
                let Some((a_pcurve2, f_par, l_par)) = curve_on_surface_range(&e, &face) else {
                    k += 1;
                    continue;
                };
                if is_same(&a_vert1.0, &a_vert2.0) == is_degenerated {
                    let a_p_eb = if e.0.orientation() == Orientation::Reversed {
                        a_pcurve2.d0(l_par)
                    } else {
                        a_pcurve2.d0(f_par)
                    };
                    if (l_par - f_par).abs() > PCONFUSION {
                        let reverse = !(e.0.orientation() == Orientation::Reversed);
                        let a_e_pm = next_param_on_pc(
                            a_pcurve2.as_ref(),
                            &a_p_eb,
                            f_par,
                            l_par,
                            self.tol_u,
                            self.tol_v,
                            reverse,
                        );
                        let mut a_p_ee = a_pcurve2.d0(a_e_pm);
                        if a_p_eb.square_distance(&a_p_ee) <= REAL_SMALL {
                            // seems to be very short curve (`cxx:625-645`).
                            let (_, a_d) = a_pcurve2.d1(a_e_pm);
                            a_p_ee = if e.0.orientation() == Orientation::Reversed {
                                GpPnt2d::new(a_p_eb.x() - a_d.x(), a_p_eb.y() - a_d.y())
                            } else {
                                GpPnt2d::new(a_p_eb.x() + a_d.x(), a_p_eb.y() + a_d.y())
                            };
                            if a_p_eb.square_distance(&a_p_ee) <= REAL_SMALL {
                                k += 1;
                                continue;
                            }
                        }
                        let e_dir =
                            GpVec2d::new(a_p_ee.x() - a_p_eb.x(), a_p_ee.y() - a_p_eb.y());
                        df_cur_angle = e_dir.angle(&e_ref_dir).abs();
                    }

                    if df_cur_angle <= df_min_angle {
                        let mut d = p_ref.square_distance(&a_p_eb);
                        if d <= PCONFUSION {
                            d = 0.;
                        }
                        if (a_p_eb.x() - p_ref.x()).abs() < self.tol_u
                            && (a_p_eb.y() - p_ref.y()).abs() < self.tol_v
                            && d <= dmin
                        {
                            df_min_angle = df_cur_angle;
                            k_min = k;
                            dmin = d;
                        }
                    }
                }
                k += 1;
            }

            if k_min == 0 {
                is_degenerated = false;
                k = 1;
                dmin = f64::MAX;
            } else {
                break;
            }
        }

        if k_min == 0 {
            // probably unclosed in 2d space wire (`cxx:681-686`).
            return None;
        }
        Some(Edge(list.remove(k_min - 1)))
    }

    /// `Current()` (`cxx:710-713`): `None` is the null edge.
    pub fn current(&self) -> Option<&Edge> {
        self.edge.as_ref()
    }

    /// `Current()` as an owned clone, for the loops that keep the shape.
    pub fn current_edge(&self) -> Option<Edge> {
        self.edge.clone()
    }

    /// `Orientation()` (`cxx:717-735`): the orientation of the vertex
    /// occurrence equal to `CurrentVertex()` inside the current edge; an
    /// infinite edge answers `FORWARD`.
    pub fn orientation(&self) -> Result<Orientation, String> {
        let Some(edge) = &self.edge else {
            return Err("Standard_NoSuchObject: BRepTools_WireExplorer::Orientation".to_string());
        };
        let Some(vertex) = &self.vertex else {
            // infinite edge
            return Ok(Orientation::Forward);
        };
        for child in stored_children(&edge.0) {
            if is_same(&child, &vertex.0) {
                return Ok(child.orientation());
            }
        }
        Err("Standard_NoSuchObject: BRepTools_WireExplorer::Orientation".to_string())
    }

    /// `CurrentVertex()` (`cxx:740-743`).
    pub fn current_vertex(&self) -> Option<&Vertex> {
        self.vertex.as_ref()
    }

    /// `Clear()` (`cxx:747-753`).
    pub fn clear(&mut self) {
        self.map.clear();
        self.doubles.clear();
        self.edge = None;
        self.face = None;
        self.vertex = None;
    }
}
