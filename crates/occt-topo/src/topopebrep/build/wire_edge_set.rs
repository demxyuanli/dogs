//! `TopOpeBRepBuild_WireEdgeSet` — the edge set of a face: a bound is a wire,
//! a bound element is an edge.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_WireEdgeSet.hxx:36-139`,
//! `.cxx:51-575` (OCCT 8.0.0).
//!
//! The `ShapeSet` stores the wires (shapes), the edges to start
//! reconstructions (start elements) and, for every vertex, the edges incident
//! to it (`hxx:36-40`). `FindNeighbours` (`cxx:142-167`) walks the vertices of
//! the current edge and keeps the first non-empty incident list;
//! `MakeNeighboursList` (`cxx:173-288`) narrows the list down when the vertex
//! is shared with closing edges of the periodic face.

use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};
use occt_core::precision::ANGULAR;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{explorer_shapes, is_equal, is_same};
use crate::topopebrep::tool::adaptor::CurveAdaptor;
use crate::topopebrep::tool::project::fun_tool_proj_pon_f;
use crate::topopebrep::tool::shape_explorer::ShapeExplorer;

use super::shape_set::{ShapeSet, ShapeSetSource};

/// `class TopOpeBRepBuild_WireEdgeSet` (`hxx:41`).
pub struct WireEdgeSet {
    /// The `TopOpeBRepBuild_ShapeSet` base (`TopAbs_VERTEX`, `cxx:53`).
    base: ShapeSet,
    /// `myFace` (`hxx:138`).
    my_face: Face,
}

impl WireEdgeSet {
    /// `TopOpeBRepBuild_WireEdgeSet(F, Addr)` (`cxx:51-56`).
    pub fn new(f: &TopoShape) -> Self {
        let face = Face::wrap(f.clone()).unwrap_or_else(Face::new);
        Self {
            base: ShapeSet::new(ShapeType::Vertex, false),
            my_face: face,
        }
    }

    /// `Face()` (`cxx:124-128`).
    pub fn face(&self) -> &Face {
        &self.my_face
    }

    /// `AddShape(S)` (`cxx:58-84`).
    ///
    /// A wire of a `U`- or `V`-closed face is added without `CheckShape`
    /// (`cxx:63-74`).
    pub fn add_shape(&mut self, s: &TopoShape) {
        let mut tocheck = true;
        let iswire = s.shape_type() == ShapeType::Wire;
        if iswire {
            if let Some(surface) = BRepTool::face_surface(&self.my_face) {
                let uc = surface.is_u_closed();
                let vc = surface.is_v_closed();
                if uc || vc {
                    tocheck = false;
                }
            }
        }
        let mut chk = true;
        if tocheck {
            chk = self.base.check_shape_of(s, false);
        }
        if !chk {
            return;
        }
        self.base.process_add_shape(s);
    }

    /// `AddStartElement(S)` (`cxx:88-116`).
    ///
    /// A B-spline or Bezier edge is added without `CheckShape`
    /// (`cxx:96-102`).
    pub fn add_start_element(&mut self, s: &TopoShape) {
        let mut tocheck = true;
        let isedge = s.shape_type() == ShapeType::Edge;
        if isedge {
            if let Some(edge) = Edge::wrap(s.clone()) {
                let cac = CurveAdaptor::of_edge(&edge);
                let t = cac.get_type();
                let b = t == crate::geom_abs_curve_type::GeomAbsCurveType::BSplineCurve
                    || t == crate::geom_abs_curve_type::GeomAbsCurveType::BezierCurve;
                tocheck = !b;
            }
        }
        let mut chk = true;
        if tocheck {
            chk = self.base.check_shape_of(s, false);
        }
        if !chk {
            return;
        }
        self.base.process_add_start_element(s);
    }

    /// `AddElement(S)` (`cxx:118-122`): the base implementation.
    pub fn add_element(&mut self, s: &TopoShape) {
        ShapeSet::add_element(&mut self.base, s);
    }

    /// `StartElements()` (`cxx:150-153`): forwarded to the `ShapeSet` base,
    /// which `TopOpeBRepBuild_WireEdgeSet` does not override.
    pub fn start_elements(&self) -> &[TopoShape] {
        self.base.start_elements()
    }

    /// `InitShapes()` — inherited from `TopOpeBRepBuild_ShapeSet`.
    pub fn init_shapes(&mut self) {
        self.base.init_shapes();
    }

    /// `MoreShapes()` — inherited from `TopOpeBRepBuild_ShapeSet`.
    pub fn more_shapes(&self) -> bool {
        self.base.more_shapes()
    }

    /// `InitStartElements()` — inherited from `TopOpeBRepBuild_ShapeSet`.
    pub fn init_start_elements(&mut self) {
        self.base.init_start_elements();
    }

    /// `MoreStartElements()` — inherited from `TopOpeBRepBuild_ShapeSet`.
    pub fn more_start_elements(&self) -> bool {
        self.base.more_start_elements()
    }

    /// `InitNeighbours(E)` (`cxx:131-137`).
    pub fn init_neighbours(&mut self, e: &TopoShape) {
        self.base
            .my_sub_shape_explorer
            .init(e, self.base.my_sub_shape_type, ShapeType::Shape);
        self.base.my_current_shape = Some(e.clone());
        self.find_neighbours();
    }

    /// `FindNeighbours()` (`cxx:142-167`).
    pub fn find_neighbours(&mut self) {
        while self.base.my_sub_shape_explorer.more() {
            let v = match self.base.my_sub_shape_explorer.current() {
                Some(c) => c.clone(),
                None => break,
            };
            let e = self.base.my_current_shape.clone();
            // `myIncidentShapesIter.Initialize(l)`: the port keeps the list in
            // `myCurrentShapeNeighbours` and resets its cursor.
            let l = self.make_neighbours_list(e.as_ref(), &v);
            self.base.my_current_shape_neighbours = l;
            self.base.my_incident_shapes_iter = 0;
            if self.base.more_neighbours() {
                break;
            } else {
                self.base.my_sub_shape_explorer.next();
            }
        }
    }

    /// `MakeNeighboursList(Earg, Varg)` (`cxx:173-288`).
    ///
    /// Without a closing shape among the incident edges the incident list is
    /// returned as is (`cxx:282-286`); otherwise the list is rebuilt from the
    /// edges that connect through `Varg` and, when at least two closing edges
    /// remain, only the one whose tangent crosses `E`'s positively is kept
    /// (`cxx:200-280`).
    pub fn make_neighbours_list(
        &mut self,
        e_arg: Option<&TopoShape>,
        v_arg: &TopoShape,
    ) -> Vec<TopoShape> {
        let Some(e_arg) = e_arg else {
            return Vec::new();
        };
        let e = match Edge::wrap(e_arg.clone()) {
            Some(e) => e,
            None => return Vec::new(),
        };
        let v = match Vertex::wrap(v_arg.clone()) {
            Some(v) => v,
            None => return Vec::new(),
        };
        let l: Vec<TopoShape> = self
            .base
            .my_sub_shape_map
            .find_from_key(v_arg)
            .cloned()
            .unwrap_or_default();

        let nclosing = self.nb_closing_shapes(&l);
        if nclosing == 0 {
            return l;
        }

        // build myCurrentShapeNeighbours = edge list made of connected shapes
        // to Earg through Varg
        let mut neighbours: Vec<TopoShape> = Vec::new();
        for curn in &l {
            let k = self.vertex_connects_edges_closing(v_arg, e_arg, curn);
            if k {
                neighbours.push(curn.clone());
            }
        }

        let newn = self.nb_closing_shapes(&neighbours);
        if newn >= 2 {
            let f = self.my_face.clone();
            let (p_e, mut d1_e) = self.d1_on_face(&e, &v);
            let _ = p_e;
            if e.orientation() == Orientation::Reversed {
                d1_e.reverse();
            }

            let mut kept: Vec<TopoShape> = Vec::new();
            for cur in &neighbours {
                if !self.is_closed(cur) {
                    kept.push(cur.clone());
                    continue;
                }
                let Some(ee) = Edge::wrap(cur.clone()) else {
                    kept.push(cur.clone());
                    continue;
                };
                let (_p_ee, mut d1_ee) = self.d1_on_face(&ee, &v);
                if ee.orientation() == Orientation::Reversed {
                    d1_ee.reverse();
                }

                let cross = d1_e.crossed(&d1_ee);
                let o_ve = self
                    .vertex_connects_edges(v_arg, e_arg, cur)
                    .map(|(o1, _)| o1)
                    .unwrap_or(Orientation::Forward);

                // cxx:268 — keep the closing edge when the relative
                // orientation of V is REVERSED for a positive cross (or
                // FORWARD for a negative one), drop it otherwise.
                let t2 = (cross > 0.0 && o_ve == Orientation::Reversed)
                    || (cross < 0.0 && o_ve == Orientation::Forward);
                if t2 {
                    kept.push(cur.clone());
                }
                let _ = &f;
            }
            return kept;
        }

        neighbours
    }

    /// `VertexConnectsEdges(V, E1, E2, O1, O2)` (`cxx:292-317`): the
    /// orientations of `V` on `E1` and `E2`, `Some` when they differ.
    pub fn vertex_connects_edges(
        &self,
        v: &TopoShape,
        e1: &TopoShape,
        e2: &TopoShape,
    ) -> Option<(Orientation, Orientation)> {
        let mut ex1 = ShapeExplorer::with_shape(e1, ShapeType::Vertex, ShapeType::Shape);
        while ex1.more() {
            let c1 = match ex1.current() {
                Some(c) => c.clone(),
                None => break,
            };
            if is_same(v, &c1) {
                let mut ex2 = ShapeExplorer::with_shape(e2, ShapeType::Vertex, ShapeType::Shape);
                while ex2.more() {
                    let c2 = match ex2.current() {
                        Some(c) => c.clone(),
                        None => break,
                    };
                    if is_same(v, &c2) {
                        let o1 = c1.orientation();
                        let o2 = c2.orientation();
                        if o1 != o2 {
                            return Some((o1, o2));
                        }
                    }
                    ex2.next();
                }
            }
            ex1.next();
        }
        None
    }

    /// `VertexConnectsEdgesClosing(V, E1, E2)` (`cxx:322-388`).
    pub fn vertex_connects_edges_closing(
        &self,
        v: &TopoShape,
        e1: &TopoShape,
        e2: &TopoShape,
    ) -> bool {
        let c1 = self.is_closed(e1);
        let c2 = self.is_closed(e2);

        let mut testconnect = c1 || c2;
        let mut resu = false;

        if c1 && c2 {
            let u1 = if c1 { self.is_u_closed(e1) } else { false };
            let v1 = if c1 { self.is_v_closed(e1) } else { false };
            let u2 = if c2 { self.is_u_closed(e2) } else { false };
            let v2 = if c2 { self.is_v_closed(e2) } else { false };
            let uvdiff = (u1 && v2) || (u2 && v1);
            testconnect = uvdiff;
        }

        if testconnect {
            resu = self.vertex_connects_edges(v, e1, e2).is_some();
        } else {
            // cto 012 O2 arete de couture de face cylindrique
            let oe1 = e1.orientation();
            let oe2 = e2.orientation();
            let iseq = is_equal(e1, e2);
            if c1 && c2 && oe1 == oe2 && !iseq {
                resu = self.vertex_connects_edges(v, e1, e2).is_some();
            }
        }
        resu
    }

    /// `NbClosingShapes(L)` (`cxx:392-405`).
    pub fn nb_closing_shapes(&self, l: &[TopoShape]) -> usize {
        let mut n = 0;
        for s in l {
            if self.is_closed(s) {
                n += 1;
            }
        }
        n
    }

    /// `LocalD1(F, E, V, pE, d1E)` (`cxx:407-439`): the UV point of `V` on `F`
    /// and the 3D tangent of `E` at `V` expressed in the surface's `(d1u, d1v)`
    /// basis.
    fn d1_on_face(&self, e: &Edge, v: &Vertex) -> (GpPnt2d, GpVec2d) {
        let par_e = BRepTool::parameter_of_vertex(&e.0, &v.0);
        let Some(ce) = BRepTool::edge_curve_world(e) else {
            return (GpPnt2d::new(0.0, 0.0), GpVec2d::from_xy(occt_core::gp::GpXY::new(0.0, 0.0)));
        };
        let (p3d, d3d) = ce.d1(par_e);
        let Some(surface) = BRepTool::face_surface(&self.my_face) else {
            return (GpPnt2d::new(0.0, 0.0), GpVec2d::from_xy(occt_core::gp::GpXY::new(0.0, 0.0)));
        };
        let Some((u, vpar, _dist)) = fun_tool_proj_pon_f(&p3d, &self.my_face) else {
            return (GpPnt2d::new(0.0, 0.0), GpVec2d::from_xy(occt_core::gp::GpXY::new(0.0, 0.0)));
        };
        let p_e = GpPnt2d::new(u, vpar);
        let (_p, d1u, d1v) = surface.d1(u, vpar);
        let du = d3d.dot(&d1u);
        let dv = d3d.dot(&d1v);
        (p_e, GpVec2d::from_xy(occt_core::gp::GpXY::new(du, dv)))
    }

    /// `IsClosed(E)` (`cxx:442-454`): `BRep_Tool::IsClosed(E, myFace)`.
    pub fn is_closed(&self, e: &TopoShape) -> bool {
        let Some(ee) = Edge::wrap(e.clone()) else {
            return false;
        };
        BRepTool::is_closed_edge_face(&ee, &self.my_face)
    }

    /// `IsUVISO(E, F, uiso, viso)` (`cxx:458-490`).
    ///
    /// A null 2D curve raises `Standard_ProgramError` (`cxx:468-471`); only a
    /// `Geom2d_Line` p-curve can be `U`- or `V`-iso (`cxx:474-489`).
    pub fn is_uv_iso(e: &Edge, f: &Face) -> (bool, bool) {
        let mut uiso = false;
        let mut viso = false;
        let pc = if let Some(r) = boptools_2d::curve_on_surface_range(e, f) {
            Some(r)
        } else {
            boptools_2d::build_pcurve_for_edge_on_face(e, f)
                .ok()
                .map(|pc| {
                    let (a, b) = BRepTool::edge_parameters(e);
                    (pc, a, b)
                })
        };
        let Some((pc, _first, _last)) = pc else {
            panic!("TopOpeBRepBuild_WireEdgeSet::IsUVISO");
        };

        if pc.is_geom2d_line() {
            let (_, d) = pc.d1(pc.first_parameter());
            if let Ok(dir) = GpDir2d::from_vec2d(&d) {
                let tol = ANGULAR;
                if let Ok(axis_y) = GpDir2d::new(0.0, 1.0) {
                    if dir.is_parallel(&axis_y, tol) {
                        uiso = true;
                    } else if let Ok(axis_x) = GpDir2d::new(1.0, 0.0) {
                        if dir.is_parallel(&axis_x, tol) {
                            viso = true;
                        }
                    }
                }
            }
        }
        (uiso, viso)
    }

    /// `IsUClosed(E)` (`cxx:493-500`).
    pub fn is_u_closed(&self, e: &TopoShape) -> bool {
        let Some(ee) = Edge::wrap(e.clone()) else {
            return false;
        };
        let (closed, _) = Self::is_uv_iso(&ee, &self.my_face);
        closed
    }

    /// `IsVClosed(E)` (`cxx:503-510`).
    pub fn is_v_closed(&self, e: &TopoShape) -> bool {
        let Some(ee) = Edge::wrap(e.clone()) else {
            return false;
        };
        let (_, closed) = Self::is_uv_iso(&ee, &self.my_face);
        closed
    }

    /// `SName(S, sb, sa)` (`cxx:539-545`).
    pub fn s_name(&self, _s: &TopoShape, sb: &str, _sa: &str) -> String {
        sb.to_string()
    }

    /// `SNameori(S, sb, sa)` (`cxx:549-555`).
    pub fn s_name_ori(&self, _s: &TopoShape, sb: &str, _sa: &str) -> String {
        sb.to_string()
    }

    /// `SName(L, sb, sa)` (`cxx:559-565`).
    pub fn s_name_list(&self, _l: &[TopoShape], _sb: &str, _sa: &str) -> String {
        String::new()
    }

    /// `SNameori(L, sb, sa)` (`cxx:569-575`).
    pub fn s_name_ori_list(&self, _l: &[TopoShape], _sb: &str, _sa: &str) -> String {
        String::new()
    }

    /// `SNameVEE(V, E1, E2)` (`cxx:516-520`).
    pub fn s_name_vee(&self, _v: &TopoShape, _e1: &TopoShape, _e2: &TopoShape) -> String {
        String::new()
    }

    /// `SNameVEL(V, E, L)` (`cxx:524-530`).
    pub fn s_name_vel(&self, _v: &TopoShape, _e: &TopoShape, _l: &[TopoShape]) -> String {
        String::new()
    }

    /// `DumpSS()` (`cxx:534`).
    pub fn dump_ss(&self) {}

    /// The edges of `S` as `TopExp_Explorer(S, TopAbs_EDGE)` yields them.
    pub fn edges_of(s: &TopoShape) -> Vec<TopoShape> {
        explorer_shapes(s, ShapeType::Edge)
    }
}

impl ShapeSetSource for WireEdgeSet {
    fn base(&self) -> &ShapeSet {
        &self.base
    }

    fn base_mut(&mut self) -> &mut ShapeSet {
        &mut self.base
    }

    fn add_shape(&mut self, s: &TopoShape) {
        WireEdgeSet::add_shape(self, s);
    }

    fn add_start_element(&mut self, s: &TopoShape) {
        WireEdgeSet::add_start_element(self, s);
    }

    fn add_element(&mut self, s: &TopoShape) {
        WireEdgeSet::add_element(self, s);
    }

    fn init_neighbours(&mut self, s: &TopoShape) {
        WireEdgeSet::init_neighbours(self, s);
    }

    fn find_neighbours(&mut self) {
        WireEdgeSet::find_neighbours(self);
    }

    fn make_neighbours_list(&mut self, e: Option<&TopoShape>, v: &TopoShape) -> Vec<TopoShape> {
        WireEdgeSet::make_neighbours_list(self, e, v)
    }
}
