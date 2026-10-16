//! `BRepSweep_Rotation` geometric hooks and `BRepSweep_Trsf::Process`.
//!
//! Sources: `BRepSweep_Rotation.cxx`, `BRepSweep_Trsf.cxx`.
//! `SetContinuity` (`BRepSweep_Trsf.cxx:97-201`) writes edge continuity through
//! `BRepLProp::Continuity`. `SplitShell` calls `BRepTools_Quilt`
//! (`BRepSweep_Rotation.cxx:797-801`).

use std::sync::Arc;

use occt_core::elib::{clib2d, slib};
use occt_core::gp::{
    GpAx1, GpAx2, GpAx22d, GpCirc, GpDir, GpDir2d, GpLin, GpPln, GpPnt, GpPnt2d, GpVec,
};
use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION};
use occt_core::toploc::TopLocLocation;
use occt_geom::curve::Curve;
use occt_geom::{GeomCircle, GeomPlane, Surface};
use occt_geom2d::{Curve2d, Geom2dCircle, Geom2dLine};

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::iterator::cumulated_children;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topexp::Explorer;
use crate::topo_tools_full::is_same;

use super::grid::{SweepOps, SweepTables};
use super::quilt::Quilt;
use super::revolution_surface::revolution_surface;
use super::num_shape::{NumShape, NumShapeIterator};

const PI: f64 = std::f64::consts::PI;

pub struct Rotation {
    pub tables: SweepTables,
    pub location: TopLocLocation,
    pub copy: bool,
    pub angle: f64,
    pub axe: GpAx1,
}

impl Rotation {
    pub fn new(
        shape: TopoShape,
        dir: NumShape,
        location: TopLocLocation,
        axe: GpAx1,
        angle: f64,
        copy: bool,
    ) -> Result<Self, String> {
        if angle < ANGULAR {
            return Err("BRepSweep_Rotation::Constructor".into());
        }
        let rot = Self {
            tables: SweepTables::new(shape, dir),
            location,
            copy,
            angle,
            axe,
        };
        rot.init();
        if let Some(err) = rot.tables.error.borrow().clone() {
            return Err(err);
        }
        Ok(rot)
    }

    /// `BRepSweep_Trsf::Init` (`BRepSweep_Trsf.cxx:47-57`).
    fn init(&self) {
        if self.copy {
            return;
        }
        let mut it = NumShapeIterator::new(self.tables.dir);
        while it.more() {
            let _ = self.process(&self.tables.gen.clone(), it.value());
            it.next();
        }
    }

    /// `BRepSweep_Trsf::Process` (`BRepSweep_Trsf.cxx:59-93`).
    fn process(&self, gen: &TopoShape, dir_v: NumShape) -> bool {
        let dotrsf = dir_v.index == 2 && !self.tables.dir.closed;
        let i_d = self.tables.dir_tool.index(dir_v) as usize;
        let i_g = self.tables.gen_tool.index(gen);
        if i_g == 0 || i_d == 0 {
            return false;
        }
        if self.is_invariant(gen) {
            self.tables.shapes.borrow_mut()[i_g][i_d] = Some(gen.clone());
            self.tables.built.borrow_mut()[i_g][i_d] = true;
            return true;
        }
        let mut touch = false;
        for child in cumulated_children(gen) {
            if self.process(&child, dir_v) {
                touch = true;
            }
        }
        if !touch || !dotrsf {
            let mut new_shape = gen.clone();
            if dotrsf {
                new_shape.move_location(&self.location);
            }
            self.tables.shapes.borrow_mut()[i_g][i_d] = Some(new_shape);
            self.tables.built.borrow_mut()[i_g][i_d] = true;
        }
        touch
    }

    pub fn shape(&self) -> Option<TopoShape> {
        super::grid::shape_result(self)
    }
}

impl SweepOps for Rotation {
    fn tables(&self) -> &SweepTables {
        &self.tables
    }
    fn make_empty_vertex(&self, gen: &TopoShape, dir: NumShape) -> TopoShape {
        if !self.copy {
            self.tables.fail("BRepSweep_Rotation::MakeEmptyVertex");
            return TopoShape::new(ShapeType::Vertex);
        }
        let mut p = BRepTool::vertex_point_world(&Vertex(gen.clone()));
        if dir.index == 2 {
            p = p.transformed(&self.location.transformation());
        }
        let tol = BRepTool::vertex_tolerance(&Vertex(gen.clone()));
        let v = TopoBuilder::new().make_vertex(p, tol);
        if dir.index == 1 && self.is_invariant(gen) && self.tables.dir_tool.nb_shapes() == 3 {
            let ig = self.tables.gen_tool.index(gen);
            self.tables.built.borrow_mut()[ig][3] = true;
            self.tables.shapes.borrow_mut()[ig][3] = Some(v.0.clone());
        }
        v.0
    }
    fn make_empty_directing_edge(&self, gen: &TopoShape, _dir: NumShape) -> TopoShape {
        let p = BRepTool::vertex_point_world(&Vertex(gen.clone()));
        let dirz = *self.axe.direction();
        let vdir = GpVec::new(dirz.x(), dirz.y(), dirz.z());
        let mut origin = *self.axe.location();
        let along = vdir.dot(&GpVec::from_pnts(&origin, &p));
        origin = GpPnt::new(
            origin.x() + along * vdir.x(),
            origin.y() + along * vdir.y(),
            origin.z() + along * vdir.z(),
        );
        let tol = BRepTool::vertex_tolerance(&Vertex(gen.clone()));
        let builder = TopoBuilder::new();
        if origin.distance(&p) <= CONFUSION {
            let axis = GpAx2::new(origin, dirz, GpDir::from_axis(occt_core::gp::dir::DirAxis::X))
                .unwrap_or_else(|_| GpAx2::default());
            let edge = builder.make_edge(
                Arc::new(GeomCircle::new(GpCirc::new(axis, 0.0))),
                0.0,
                2.0 * PI,
            );
            GeometryRegistry::global().set_degenerated(&edge.0, true);
            edge.0
        } else {
            let xdir = GpDir::from_vec(&GpVec::from_pnts(&origin, &p)).unwrap_or(dirz);
            let axis = GpAx2::new(origin, dirz, xdir).unwrap_or_else(|_| GpAx2::default());
            let radius = origin.distance(&p);
            builder
                .make_edge(Arc::new(GeomCircle::new(GpCirc::new(axis, radius))), 0.0, 2.0 * PI)
                .0
        }
    }
    fn make_empty_generating_edge(&self, gen: &TopoShape, dir: NumShape) -> TopoShape {
        let src = Edge(gen.clone());
        let builder = TopoBuilder::new();
        let edge = if BRepTool::is_degenerated(&src) {
            let tol = BRepTool::edge_tolerance(&src);
            let (first, last) = BRepTool::edge_parameters(&src);
            let e = TopoShape::new(ShapeType::Edge);
            {
                let mut ts = e.tshape.write().expect("tshape");
                let c = ts.edge_core_mut();
                c.first = first;
                c.last = last;
                c.tolerance = tol;
                c.degenerated = true;
                c.same_parameter = true;
                c.same_range = true;
            }
            e
        } else {
            let (first, last) = BRepTool::edge_parameters(&src);
            let tol = BRepTool::edge_tolerance(&src);
            let curve = BRepTool::edge_curve(&src).map(|c| {
                let tr = gen.location().transformation();
                let moved: Arc<dyn Curve> = Arc::from(c.transformed(&tr));
                if dir.index == 2 {
                    let sweep = self.location.transformation();
                    Arc::from(moved.transformed(&sweep))
                } else {
                    moved
                }
            });
            if let Some(curve) = curve {
                let e = builder.make_edge(curve, first, last);
                GeometryRegistry::global().raise_edge_tolerance(&e.0, tol);
                e.0
            } else {
                TopoShape::new(ShapeType::Edge)
            }
        };
        if dir.index == 1 && self.is_invariant(gen) && self.tables.dir_tool.nb_shapes() == 3 {
            let ig = self.tables.gen_tool.index(gen);
            self.tables.built.borrow_mut()[ig][3] = true;
            self.tables.shapes.borrow_mut()[ig][3] = Some(edge.clone());
        }
        edge
    }
    fn make_empty_face(&self, gen: &TopoShape, dir: NumShape) -> TopoShape {
        let (surface, tol) = if gen.shape_type() == ShapeType::Edge {
            let edge = Edge(gen.clone());
            let tol = BRepTool::edge_tolerance(&edge);
            let (first, last) = BRepTool::edge_parameters(&edge);
            let curve = BRepTool::edge_curve_world(&edge);
            let surface = curve
                .map(|c| revolution_surface(c, first, last, &self.axe))
                .unwrap_or_else(|| Arc::new(GeomPlane::new(GpPln::default())) as Arc<dyn Surface>);
            (surface, tol)
        } else {
            let face = Face(gen.clone());
            let tol = BRepTool::face_tolerance(&face);
            let mut surface = BRepTool::face_surface_world(&face)
                .unwrap_or_else(|| Arc::new(GeomPlane::new(GpPln::default())) as Arc<dyn Surface>);
            if dir.index == 2 {
                let sweep = self.location.transformation();
                surface = Arc::from(surface.transformed(&sweep));
            }
            (surface, tol)
        };
        let face = TopoBuilder::new().make_face(surface, &[]);
        GeometryRegistry::global().set_face_tolerance(&face.0, tol);
        face.0
    }
    fn set_parameters(
        &self,
        new_face: &TopoShape,
        new_vertex: &TopoShape,
        gen_f: &TopoShape,
        gen_v: &TopoShape,
        _dir: NumShape,
    ) {
        // UV of the generating vertex on the generating face.
        // `BRep_Tool::Parameters` is the stored PointOnSurface; absent storage
        // projects. Cap faces are copies, so the generating UV is reused.
        // `BRep_Tool::Parameters` (`BRepSweep_Rotation.cxx:261-274`) reads a
        // stored PointOnSurface. This port has no such slot, so the cap-face
        // UV is not written.
        let _ = (new_face, new_vertex, gen_f, gen_v);
    }
    fn set_directing_parameter(
        &self,
        new_edge: &TopoShape,
        new_vertex: &TopoShape,
        _gen_v: &TopoShape,
        _dir_e: NumShape,
        dir_v: NumShape,
    ) {
        let mut param = 0.0;
        let mut ori = Orientation::Forward;
        if dir_v.index == 2 {
            param = self.angle;
            ori = Orientation::Reversed;
        }
        let mut v = new_vertex.clone();
        v.set_orientation(ori);
        set_vertex_param(new_edge, &v, param);
    }
    fn set_generating_parameter(
        &self,
        new_edge: &TopoShape,
        new_vertex: &TopoShape,
        gen_e: &TopoShape,
        gen_v: &TopoShape,
        _dir: NumShape,
    ) {
        let mut v = new_vertex.clone();
        v.set_orientation(gen_v.orientation());
        let param = vertex_parameter(gen_v, gen_e);
        set_vertex_param(new_edge, &v, param);
    }
    fn set_pcurve(
        &self,
        new_face: &TopoShape,
        new_edge: &TopoShape,
        gen_f: &TopoShape,
        gen_e: &TopoShape,
        _dir: NumShape,
        orien: Orientation,
    ) {
        let face_key = GeometryRegistry::shape_key(gen_f);
        if let Some(pc) = GeometryRegistry::global().edge_pcurve(gen_e, face_key) {
            set_the_pcurve(new_edge, new_face, orien, pc);
        }
    }
    fn set_generating_pcurve(
        &self,
        new_face: &TopoShape,
        new_edge: &TopoShape,
        _gen_e: &TopoShape,
        _dir_e: NumShape,
        dir_v: NumShape,
        orien: Orientation,
    ) {
        let line = generating_line(self, new_face, new_edge, dir_v);
        set_the_pcurve(new_edge, new_face, orien, Arc::new(line));
    }
    fn set_directing_pcurve(
        &self,
        new_face: &TopoShape,
        new_edge: &TopoShape,
        gen_e: &TopoShape,
        gen_v: &TopoShape,
        _dir_e: NumShape,
        orien: Orientation,
    ) {
        let pc = directing_pcurve(self, new_face, gen_e, gen_v);
        set_the_pcurve(new_edge, new_face, orien, pc);
    }
    fn direct_solid(&self, gen: &TopoShape, _dir: NumShape) -> Orientation {
        direct_solid(gen, &self.axe)
    }
    fn ggd_is_to_add(
        &self,
        new_shape: &TopoShape,
        new_sub: &TopoShape,
        gen: &TopoShape,
        sub_gen: &TopoShape,
        dir: NumShape,
    ) -> bool {
        if new_shape.shape_type() == ShapeType::Face
            && new_sub.shape_type() == ShapeType::Edge
            && gen.shape_type() == ShapeType::Edge
            && sub_gen.shape_type() == ShapeType::Vertex
            && dir.ty == ShapeType::Edge
        {
            if face_is_plane(new_shape) {
                return !self.is_invariant(sub_gen);
            }
        }
        true
    }
    fn gdd_is_to_add(
        &self,
        new_shape: &TopoShape,
        new_sub: &TopoShape,
        gen: &TopoShape,
        dir: NumShape,
        sub_dir: NumShape,
    ) -> bool {
        let full = (self.angle - 2.0 * PI).abs() > ANGULAR;
        if new_shape.shape_type() == ShapeType::Solid
            && new_sub.shape_type() == ShapeType::Face
            && gen.shape_type() == ShapeType::Face
            && dir.ty == ShapeType::Edge
            && sub_dir.ty == ShapeType::Vertex
        {
            return full;
        }
        if new_shape.shape_type() == ShapeType::Face
            && new_sub.shape_type() == ShapeType::Edge
            && gen.shape_type() == ShapeType::Edge
            && dir.ty == ShapeType::Edge
            && sub_dir.ty == ShapeType::Vertex
            && face_is_plane(new_shape)
        {
            return full;
        }
        true
    }
    fn separated_wires(
        &self,
        new_shape: &TopoShape,
        new_sub: &TopoShape,
        gen: &TopoShape,
        sub_gen: &TopoShape,
        dir: NumShape,
    ) -> bool {
        if new_shape.shape_type() == ShapeType::Face
            && new_sub.shape_type() == ShapeType::Edge
            && gen.shape_type() == ShapeType::Edge
            && sub_gen.shape_type() == ShapeType::Vertex
            && dir.ty == ShapeType::Edge
            && face_is_plane(new_shape)
        {
            return (self.angle - 2.0 * PI).abs() <= ANGULAR;
        }
        false
    }
    fn split_shell(&self, shell: &TopoShape) -> TopoShape {
        // `BRepSweep_Rotation::SplitShell` (`cxx:797-801`).
        let mut quilt = Quilt::new();
        quilt.add(shell);
        quilt.shells()
    }
    fn set_continuity(&self, gen: &TopoShape, dir: NumShape) {
        super::continuity::set_continuity(self, gen, dir);
    }
    fn has_shape(&self, gen: &TopoShape, dir: NumShape) -> bool {
        if dir.ty == ShapeType::Edge && gen.shape_type() == ShapeType::Edge {
            let edge = Edge(gen.clone());
            if BRepTool::is_degenerated(&edge) {
                return false;
            }
            if BRepTool::edge_curve(&edge).is_none() {
                return false;
            }
            if self.is_invariant(gen) {
                return false;
            }
            let mut exp = Explorer::new(&self.tables.gen, ShapeType::Face);
            while exp.more() {
                if is_really_closed(gen, exp.current()) {
                    return false;
                }
                exp.next();
            }
        }
        true
    }
    fn is_invariant(&self, gen: &TopoShape) -> bool {
        is_invariant(gen, &self.axe)
    }
}

fn is_invariant(gen: &TopoShape, axe: &GpAx1) -> bool {
    let lin = GpLin::from_pnt_dir(*axe.location(), *axe.direction());
    if gen.shape_type() == ShapeType::Edge {
        let edge = Edge(gen.clone());
        let Some(curve) = BRepTool::edge_curve_world(&edge) else {
            return false;
        };
        let line = curve.gp_line().is_some();
        let spline = curve.bspline_poles().is_some();
        let bezier = curve.bezier_poles().is_some();
        if !(line || spline || bezier) {
            return false;
        }
        let verts = cumulated_children(gen);
        if verts.len() < 2 {
            return false;
        }
        if !is_invariant(&verts[0], axe) || !is_invariant(&verts[1], axe) {
            return false;
        }
        if line {
            return true;
        }
        let tol = BRepTool::vertex_tolerance(&Vertex(verts[0].clone()))
            .max(BRepTool::vertex_tolerance(&Vertex(verts[1].clone())));
        let poles: Vec<GpPnt> = curve
            .bspline_poles()
            .map(|p| p.to_vec())
            .or_else(|| curve.bezier_poles().map(|p| p.to_vec()))
            .unwrap_or_default();
        poles.iter().all(|p| lin.distance(p) <= tol)
    } else if gen.shape_type() == ShapeType::Vertex {
        let p = BRepTool::vertex_point_world(&Vertex(gen.clone()));
        lin.distance(&p) <= BRepTool::vertex_tolerance(&Vertex(gen.clone()))
    } else {
        false
    }
}

pub(super) fn is_really_closed(edge: &TopoShape, face: &TopoShape) -> bool {
    let mut fwd = false;
    let mut rev = false;
    for wire in cumulated_children(face) {
        for e in cumulated_children(&wire) {
            if is_same(&e, edge) {
                if e.orientation() == Orientation::Forward {
                    fwd = true;
                }
                if e.orientation() == Orientation::Reversed {
                    rev = true;
                }
            }
        }
    }
    fwd && rev
}

fn face_is_plane(face: &TopoShape) -> bool {
    BRepTool::face_surface(&Face(face.clone()))
        .and_then(|s| s.gp_pln())
        .is_some()
}

pub(super) fn set_vertex_param(edge: &TopoShape, vertex: &TopoShape, param: f64) {
    let key = GeometryRegistry::shape_key(vertex);
    let mut ts = edge.tshape.write().expect("tshape");
    let core = ts.edge_core_mut();
    if let Some(slot) = core.vertex_params.iter_mut().find(|(k, _)| *k == key) {
        slot.1 = param;
    } else {
        core.vertex_params.push((key, param));
    }
}

pub(super) fn vertex_parameter(vertex: &TopoShape, edge: &TopoShape) -> f64 {
    let key = GeometryRegistry::shape_key(vertex);
    if let Ok(ts) = edge.tshape.read() {
        if let Some(core) = ts.edge_core() {
            if let Some((_, p)) = core.vertex_params.iter().find(|(k, _)| *k == key) {
                return *p;
            }
        }
    }
    let (first, last) = BRepTool::edge_parameters(&Edge(edge.clone()));
    if vertex.orientation() == Orientation::Reversed {
        last
    } else {
        first
    }
}

pub(super) fn set_the_pcurve(edge: &TopoShape, face: &TopoShape, orien: Orientation, curve: Arc<dyn Curve2d>) {
    if BRepTool::is_degenerated(&Edge(edge.clone())) {
        GeometryRegistry::global().raise_edge_tolerance(edge, BRepTool::edge_tolerance(&Edge(edge.clone())));
    }
    let tol = pcurve_tolerance(edge, face, curve.as_ref());
    let face_key = GeometryRegistry::shape_key(face);
    let plane = face_is_plane(face);
    let existing = if plane {
        None
    } else {
        GeometryRegistry::global().edge_pcurve(edge, face_key)
    };
    let curves = match existing {
        None => vec![curve],
        Some(old) if orien == Orientation::Reversed => vec![old, curve],
        Some(old) => vec![curve, old],
    };
    GeometryRegistry::global().set_edge_pcurves(edge, face_key, curves);
    GeometryRegistry::global().raise_edge_tolerance(edge, tol);
}

fn pcurve_tolerance(edge: &TopoShape, face: &TopoShape, pc: &dyn Curve2d) -> f64 {
    let e = Edge(edge.clone());
    if BRepTool::is_degenerated(&e) {
        return BRepTool::edge_tolerance(&e);
    }
    let Some(c3d) = BRepTool::edge_curve_world(&e) else {
        return BRepTool::edge_tolerance(&e);
    };
    let Some(surf) = BRepTool::face_surface_world(&Face(face.clone())) else {
        return BRepTool::edge_tolerance(&e);
    };
    let (first, last) = BRepTool::edge_parameters(&e);
    let mut d2: f64 = 0.0;
    let nn = 23.0;
    for i in 0..=23 {
        let t = (i as f64) / nn;
        let u = first * (1.0 - t) + last * t;
        let pc3d = c3d.value(u);
        let uv = pc.d0(u);
        let cons = surf.d0(uv.x(), uv.y());
        if !cons.x().is_finite() || !cons.y().is_finite() || !cons.z().is_finite() {
            return f64::MAX;
        }
        d2 = d2.max(pc3d.square_distance(&cons));
    }
    let mut d = 1.5 * d2.sqrt();
    if d < 1.0e-7 {
        d = 1.0e-7;
    }
    d
}

fn generating_line(rot: &Rotation, face: &TopoShape, edge: &TopoShape, dir_v: NumShape) -> Geom2dLine {
    let surf = BRepTool::face_surface_world(&Face(face.clone()));
    let dx = GpDir2d::new(1.0, 0.0).unwrap_or_default();
    let dy = GpDir2d::new(0.0, 1.0).unwrap_or_default();
    let fallback = |ang: f64| Geom2dLine::from_pnt_dir(GpPnt2d::new(ang, 0.0), dy);
    let Some(surf) = surf else {
        return fallback(0.0);
    };
    if let Some(pln) = surf.gp_pln() {
        let Some(curve) = BRepTool::edge_curve_world(&Edge(edge.clone())) else {
            rot.tables.fail("BRepSweep_Rotation::SetGeneratingPCurve");
            return fallback(0.0);
        };
        let Some(gl) = curve.gp_line() else {
            rot.tables.fail("BRepSweep_Rotation::SetGeneratingPCurve");
            return fallback(0.0);
        };
        let loc = gl.location();
        let (u, v) = slib::plane_parameters(&pln.position(), &loc);
        let dir = gl.direction();
        let ax3 = pln.position();
        let d2 = GpDir2d::new(dir.dot(ax3.x_direction()), dir.dot(ax3.y_direction()))
            .unwrap_or(dx);
        return Geom2dLine::from_pnt_dir(GpPnt2d::new(u, v), d2);
    }
    if let Some(tor) = surf.gp_torus() {
        let Some(curve) = BRepTool::edge_curve_world(&Edge(edge.clone())) else {
            return fallback(0.0);
        };
        let u0 = curve.first_parameter();
        let point = curve.value(u0);
        let (mut u, mut v) = if point.distance(&tor.location()) < CONFUSION {
            (0.0, PI)
        } else {
            slib::torus_parameters(&tor.position(), tor.major_radius(), tor.minor_radius(), &point)
        };
        v = clib2d::in_period(v, 0.0, 2.0 * PI);
        if (2.0 * PI - v) <= PCONFUSION {
            v -= 2.0 * PI;
        }
        if dir_v.index == 2 {
            let mut u_left = u - rot.angle;
            clib2d::adjust_periodic(-PI, PI, PCONFUSION, &mut u_left, &mut u);
        } else {
            let mut u_right = u + rot.angle;
            clib2d::adjust_periodic(-PI, PI, PCONFUSION, &mut u, &mut u_right);
        }
        return Geom2dLine::from_pnt_dir(GpPnt2d::new(u, v - u0), dy);
    }
    if let Some(sph) = surf.gp_sphere() {
        let Some(curve) = BRepTool::edge_curve_world(&Edge(edge.clone())) else {
            return fallback(0.0);
        };
        let u0 = curve.first_parameter();
        let point = curve.value(u0);
        let (_, v) = slib::sphere_parameters(&sph.position(), &point);
        let u = if dir_v.index == 2 { rot.angle } else { 0.0 };
        return Geom2dLine::from_pnt_dir(GpPnt2d::new(u, v - u0), dy);
    }
    let ang = if dir_v.index == 2 { rot.angle } else { 0.0 };
    fallback(ang)
}

fn directing_pcurve(rot: &Rotation, face: &TopoShape, gen_e: &TopoShape, gen_v: &TopoShape) -> Arc<dyn Curve2d> {
    let dx = GpDir2d::new(1.0, 0.0).unwrap_or_default();
    let dy = GpDir2d::new(0.0, 1.0).unwrap_or_default();
    let par = vertex_parameter(gen_v, gen_e);
    let p2 = BRepTool::vertex_point_world(&Vertex(gen_v.clone()));
    let Some(surf) = BRepTool::face_surface_world(&Face(face.clone())) else {
        return Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, par), dx));
    };
    if let Some(pln) = surf.gp_pln() {
        let p1 = pln.location();
        let r = p1.distance(&p2);
        let (u, v) = slib::plane_parameters(&pln.position(), &p2);
        let d = GpDir2d::new(u, v).unwrap_or(dx);
        if let Ok(ax) = GpAx22d::new(GpPnt2d::new(0.0, 0.0), d, dy) {
            return Arc::new(Geom2dCircle::new(occt_core::gp::GpCirc2d::new(ax, r)));
        }
    } else if let Some(cone) = surf.gp_cone() {
        let (_, v) = slib::cone_parameters(&cone.position(), cone.radius(), cone.semi_angle(), &p2);
        return Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, v), dx));
    } else if let Some(sph) = surf.gp_sphere() {
        let (_, v) = slib::sphere_parameters(&sph.position(), &p2);
        return Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, v), dx));
    } else if let Some(tor) = surf.gp_torus() {
        let Some(curve) = BRepTool::edge_curve_world(&Edge(gen_e.clone())) else {
            return Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, par), dx));
        };
        let p1 = curve.value(curve.first_parameter());
        let (u1, v1) = if p1.distance(&tor.location()) < CONFUSION {
            (0.0, PI)
        } else {
            slib::torus_parameters(&tor.position(), tor.major_radius(), tor.minor_radius(), &p1)
        };
        let p2c = curve.value(curve.last_parameter());
        let (_, v2) = if p2c.distance(&tor.location()) < CONFUSION {
            (0.0, PI)
        } else {
            slib::torus_parameters(&tor.position(), tor.major_radius(), tor.minor_radius(), &p2c)
        };
        let mut v1 = v1;
        let mut v2 = v2;
        clib2d::adjust_periodic(0.0, 2.0 * PI, PCONFUSION, &mut v1, &mut v2);
        let mut u1 = u1;
        let mut u2 = u1 + rot.angle;
        clib2d::adjust_periodic(-PI, PI, PCONFUSION, &mut u1, &mut u2);
        let vv = if gen_v.orientation() == Orientation::Forward { v1 } else { v2 };
        return Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(u1, vv), dx));
    }
    Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, par), dx))
}

fn direct_solid(gen: &TopoShape, axe: &GpAx1) -> Orientation {
    let face = Face(gen.clone());
    let Some(surf) = BRepTool::face_surface_world(&face) else {
        return Orientation::Forward;
    };
    // `BRepSweep_Rotation::DirectSolid` samples `BRepAdaptor_Surface(TopoDS::Face(aGenS))`
    // (`BRepSweep_Rotation.cxx:641`, `:649-657`), whose `FirstUParameter` /
    // `LastUParameter` / `FirstVParameter` / `LastVParameter` are the **face's**
    // UV box, not the surface's natural range: its constructor defaults to
    // `Restriction = true` and runs `BRepTools::UVBounds(F, ...)`
    // (`BRepAdaptor_Surface.cxx:57-81`). A plane's natural range is infinite,
    // so sampling the surface range puts the midpoint at a non-point of the
    // face and the sign of `aX` degenerates.
    let (u1, u2, v1, v2) = crate::brep_surface::face_uv_bounds(&face);
    let mut tx = 0.5;
    let mut ux = tx * (u1 + u2);
    let mut vx = tx * (v1 + v2);
    let (mut p, mut du, mut dv) = surf.d1(ux, vx);
    let loc = *axe.location();
    let dir = GpVec::new(axe.direction().x(), axe.direction().y(), axe.direction().z());
    let mut radial = GpVec::from_pnts(&loc, &p);
    radial = radial.crossed(&dir);
    let tol2 = CONFUSION * CONFUSION;
    if radial.square_magnitude() < tol2 {
        tx = 0.43213918;
        ux = u1 * (1.0 - tx) + u2 * tx;
        vx = v1 * (1.0 - tx) + v2 * tx;
        let sample = surf.d1(ux, vx);
        p = sample.0;
        du = sample.1;
        dv = sample.2;
        radial = GpVec::from_pnts(&loc, &p).crossed(&dir);
    }
    let x = radial.dot(&du.crossed(&dv));
    if x > 0.0 {
        Orientation::Forward
    } else {
        Orientation::Reversed
    }
}

