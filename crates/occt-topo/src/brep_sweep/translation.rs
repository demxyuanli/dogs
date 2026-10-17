//! `BRepSweep_Translation` geometric hooks and `BRepPrimAPI_MakePrism`.
//!
//! Sources: `BRepSweep_Translation.cxx`, `BRepSweep_Prism.cxx`,
//! `BRepPrimAPI_MakePrism.cxx`. `Copy` is false and `Canonize` is true for
//! `BRepPrimAPI_MakePrism(S, V)` (`BRepPrimAPI_MakePrism.hxx:56-58`).
//! The plane arms of `SetGeneratingPCurve` and `SetDirectingPCurve` are
//! commented out (`BRepSweep_Translation.cxx:336-364`, `:396-411`).

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpDir2d, GpLin, GpPnt, GpPnt2d, GpTrsf, GpVec};
use occt_core::precision::{Precision, ANGULAR, CONFUSION};
use occt_core::toploc::TopLocLocation;
use occt_geom::curve::Curve;
use occt_geom::{GeomCylinder, GeomLine, GeomPlane, GeomSurfaceOfLinearExtrusion, Surface};

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::iterator::cumulated_children;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topexp::Explorer;

use super::grid::{SweepOps, SweepTables};
use super::num_shape::{NumShape, NumShapeIterator};
use super::quilt::Quilt;
use super::rotation::{is_really_closed, set_the_pcurve, set_vertex_param, vertex_parameter};

const PI: f64 = std::f64::consts::PI;

pub struct Translation {
    pub tables: SweepTables,
    pub location: TopLocLocation,
    pub copy: bool,
    pub canonize: bool,
    pub vec: GpVec,
}

impl Translation {
    pub fn new(
        shape: TopoShape,
        dir: NumShape,
        location: TopLocLocation,
        vec: GpVec,
        copy: bool,
        canonize: bool,
    ) -> Result<Self, String> {
        if vec.magnitude() < CONFUSION {
            return Err("BRepSweep_Translation::Constructor".into());
        }
        let sweep = Self {
            tables: SweepTables::new(shape, dir),
            location,
            copy,
            canonize,
            vec,
        };
        sweep.init();
        if let Some(err) = sweep.tables.error.borrow().clone() {
            return Err(err);
        }
        Ok(sweep)
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

/// `BRepPrimAPI_MakePrism(S, V)` (`BRepPrimAPI_MakePrism.cxx:36-42`).
/// `Copy` is false and `Canonize` is true.
pub fn make_prism(shape: &TopoShape, vec: &GpVec) -> Result<TopoShape, String> {
    let _ = crate::brep_lib_build_curves3d::build_curves3d(shape);
    if vec.magnitude() <= CONFUSION {
        return Err("BRepSweep_Prism::Constructor".into());
    }
    let mut num = NumShape::dummy();
    num.init(2, ShapeType::Edge, false, false, false);
    let mut gpt = GpTrsf::identity();
    gpt.set_translation_vec(vec);
    let location = TopLocLocation::composed(&gpt, TopLocLocation::identity());
    let sweep = Translation::new(shape.clone(), num, location, *vec, false, true)?;
    sweep
        .shape()
        .ok_or_else(|| "BRepSweep_Prism::Shape".to_string())
}

impl SweepOps for Translation {
    fn tables(&self) -> &SweepTables {
        &self.tables
    }
    fn make_empty_vertex(&self, gen: &TopoShape, dir: NumShape) -> TopoShape {
        // `cxx:104-119`. Called only when `myCopy` is true.
        if !self.copy {
            self.tables.fail("BRepSweep_Translation::MakeEmptyVertex");
            return TopoShape::new(ShapeType::Vertex);
        }
        let mut p = BRepTool::vertex_point_world(&Vertex(gen.clone()));
        if dir.index == 2 {
            p = p.transformed(&self.location.transformation());
        }
        let tol = BRepTool::vertex_tolerance(&Vertex(gen.clone()));
        TopoBuilder::new().make_vertex(p, tol).0
    }
    fn make_empty_directing_edge(&self, gen: &TopoShape, _dir: NumShape) -> TopoShape {
        // `cxx:124-132`: `gp_Lin(P, myVec)`.
        let p = BRepTool::vertex_point_world(&Vertex(gen.clone()));
        let Ok(dir) = GpDir::from_vec(&self.vec) else {
            self.tables.fail("BRepSweep_Translation::MakeEmptyDirectingEdge");
            return TopoShape::new(ShapeType::Edge);
        };
        let tol = BRepTool::vertex_tolerance(&Vertex(gen.clone()));
        let edge = TopoBuilder::new().make_edge(
            Arc::new(GeomLine::new(GpLin::from_pnt_dir(p, dir))),
            0.0,
            self.vec.magnitude(),
        );
        GeometryRegistry::global().raise_edge_tolerance(&edge.0, tol);
        edge.0
    }
    fn make_empty_generating_edge(&self, gen: &TopoShape, dir: NumShape) -> TopoShape {
        // `cxx:137-165`. Called only when `myCopy` is true.
        if !self.copy {
            self.tables.fail("BRepSweep_Translation::MakeEmptyVertex");
            return TopoShape::new(ShapeType::Edge);
        }
        let src = Edge(gen.clone());
        let builder = TopoBuilder::new();
        if BRepTool::is_degenerated(&src) {
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
            return e;
        }
        let (first, last) = BRepTool::edge_parameters(&src);
        let tol = BRepTool::edge_tolerance(&src);
        let curve = BRepTool::edge_curve(&src).map(|c| {
            let tr = gen.location().transformation();
            let moved: Arc<dyn Curve> = Arc::from(c.transformed(&tr));
            if dir.index == 2 {
                Arc::from(moved.transformed(&self.location.transformation()))
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
    }
    fn make_empty_face(&self, gen: &TopoShape, dir: NumShape) -> TopoShape {
        // `cxx:222-278`.
        let (surface, tol) = if dir.ty == ShapeType::Edge {
            let edge = Edge(gen.clone());
            let tol = BRepTool::edge_tolerance(&edge);
            let (first, last) = BRepTool::edge_parameters(&edge);
            let Some(curve) = BRepTool::edge_curve(&edge) else {
                self.tables.fail("BRepSweep_Translation::MakeEmptyFace");
                return TopoShape::new(ShapeType::Face);
            };
            let tr = gen.location().transformation();
            let curve: Arc<dyn Curve> = Arc::from(curve.transformed(&tr));
            let Ok(mut direction) = GpDir::from_vec(&self.vec) else {
                self.tables.fail("BRepSweep_Translation::MakeEmptyFace");
                return TopoShape::new(ShapeType::Face);
            };
            direction.reverse();
            let surface: Arc<dyn Surface> = if self.canonize {
                canon_extrusion(curve.as_ref(), direction, first, last)
            } else {
                Arc::new(GeomSurfaceOfLinearExtrusion::new(curve, direction))
            };
            (surface, tol)
        } else {
            let face = Face(gen.clone());
            let tol = BRepTool::face_tolerance(&face);
            let mut surface = BRepTool::face_surface_world(&face)
                .unwrap_or_else(|| Arc::new(GeomPlane::new(occt_core::gp::GpPln::default())));
            if dir.index == 2 {
                let mut shift = GpTrsf::identity();
                shift.set_translation_vec(&self.vec);
                surface = Arc::from(surface.transformed(&shift));
            }
            (surface, tol)
        };
        let face = TopoBuilder::new().make_face(surface, &[]);
        GeometryRegistry::global().set_face_tolerance(&face.0, tol);
        if std::env::var("SWEEP_DBG").is_ok() {
            eprintln!(
                "[mef] dir idx={} ty={:?} gen_loc={:?} face_loc={:?} surf={}",
                dir.index,
                dir.ty,
                gen.location().transformation(),
                face.0.location().transformation(),
                face.0.tshape.read().unwrap().face_core().is_some(),
            );
        }
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
        // `BRep_Tool::Parameters` (`cxx:177`) reads a stored PointOnSurface.
        // This port has no such slot, so the cap-face UV is not written.
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
        // `cxx:187-201`.
        let param = if dir_v.index == 2 {
            self.vec.magnitude()
        } else {
            0.0
        };
        set_vertex_param(new_edge, new_vertex, param);
    }
    fn set_generating_parameter(
        &self,
        new_edge: &TopoShape,
        new_vertex: &TopoShape,
        gen_e: &TopoShape,
        gen_v: &TopoShape,
        _dir: NumShape,
    ) {
        // `cxx:206-217`.
        let mut v = new_vertex.clone();
        v.set_orientation(gen_v.orientation());
        set_vertex_param(new_edge, &v, vertex_parameter(gen_v, gen_e));
    }
    fn set_pcurve(
        &self,
        new_face: &TopoShape,
        new_edge: &TopoShape,
        gen_f: &TopoShape,
        gen_e: &TopoShape,
        _dir: NumShape,
        _orien: Orientation,
    ) {
        // `cxx:283-316`. A closed edge copies both pcurves.
        let reg = GeometryRegistry::global();
        let face_key = GeometryRegistry::shape_key(gen_f);
        let closed = BRepTool::is_closed_edge_face(&Edge(gen_e.clone()), &Face(gen_f.clone()));
        if closed {
            let curves = reg.edge_pcurves(gen_e, face_key);
            if curves.len() >= 2 {
                reg.set_edge_pcurves(new_edge, GeometryRegistry::shape_key(new_face), curves);
                return;
            }
        }
        if let Some(pc) = reg.edge_pcurve(gen_e, face_key) {
            set_the_pcurve(new_edge, new_face, Orientation::Forward, pc);
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
        // `cxx:321-377`. The plane arm is commented out.
        if face_is_plane(new_face) {
            return;
        }
        let v = if dir_v.index == 2 {
            -self.vec.magnitude()
        } else {
            0.0
        };
        let Ok(dx) = GpDir2d::new(1.0, 0.0) else {
            return;
        };
        let line = occt_geom2d::Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, v), dx);
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
        // `cxx:382-415`. The plane arm is commented out.
        if face_is_plane(new_face) {
            return;
        }
        let Ok(ny) = GpDir2d::new(0.0, -1.0) else {
            return;
        };
        let par = vertex_parameter(gen_v, gen_e);
        let line = occt_geom2d::Geom2dLine::from_pnt_dir(GpPnt2d::new(par, 0.0), ny);
        set_the_pcurve(new_edge, new_face, orien, Arc::new(line));
    }
    fn direct_solid(&self, gen: &TopoShape, _dir: NumShape) -> Orientation {
        // `cxx:420-435`.
        let Some(surf) = BRepTool::face_surface_world(&Face(gen.clone())) else {
            return Orientation::Forward;
        };
        let (u1, u2) = surf.u_range();
        let (v1, v2) = surf.v_range();
        let (_, du, dv) = surf.d1((u1 + u2) / 2.0, (v1 + v2) / 2.0);
        if self.vec.dot(&du.crossed(&dv)) > 0.0 {
            Orientation::Reversed
        } else {
            Orientation::Forward
        }
    }
    fn ggd_is_to_add(
        &self,
        _new_shape: &TopoShape,
        _new_sub: &TopoShape,
        _gen: &TopoShape,
        _sub_gen: &TopoShape,
        _dir: NumShape,
    ) -> bool {
        true
    }
    fn gdd_is_to_add(
        &self,
        _new_shape: &TopoShape,
        _new_sub: &TopoShape,
        _gen: &TopoShape,
        _dir: NumShape,
        _sub_dir: NumShape,
    ) -> bool {
        true
    }
    fn separated_wires(
        &self,
        _new_shape: &TopoShape,
        _new_sub: &TopoShape,
        _gen: &TopoShape,
        _sub_gen: &TopoShape,
        _dir: NumShape,
    ) -> bool {
        false
    }
    fn split_shell(&self, shell: &TopoShape) -> TopoShape {
        let mut quilt = Quilt::new();
        quilt.add(shell);
        quilt.shells()
    }
    fn set_continuity(&self, gen: &TopoShape, dir: NumShape) {
        super::continuity::set_continuity(self, gen, dir);
    }
    fn has_shape(&self, gen: &TopoShape, dir: NumShape) -> bool {
        // `cxx:473-517`. The sewing-edge list walk is commented out; the
        // replacement is `BRepTools::IsReallyClosed`.
        if dir.ty == ShapeType::Edge && gen.shape_type() == ShapeType::Edge {
            let edge = Edge(gen.clone());
            if BRepTool::is_degenerated(&edge) {
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
    fn is_invariant(&self, _gen: &TopoShape) -> bool {
        false
    }
}

fn face_is_plane(face: &TopoShape) -> bool {
    BRepTool::face_surface(&Face(face.clone()))
        .and_then(|s| s.gp_pln())
        .is_some()
}

/// `GeomAdaptor_SurfaceOfLinearExtrusion::GetType` (`cxx:269-328`), `Plane`
/// (`cxx:332-370`) and `Cylinder` (`cxx:375-386`). Any other type stays
/// `Geom_SurfaceOfLinearExtrusion`.
fn canon_extrusion(
    curve: &dyn Curve,
    direction: GpDir,
    first: f64,
    last: f64,
) -> Arc<dyn Surface> {
    let extrusion = Arc::new(GeomSurfaceOfLinearExtrusion::new(
        Arc::from(curve.clone_dyn()),
        direction,
    ));
    if let Some(line) = curve.gp_line() {
        let line_dir = line.direction();
        if !direction.is_parallel_tol(&line_dir, ANGULAR) {
            if let Some(pln) = extrusion_plane(curve, direction, first, last) {
                return Arc::new(GeomPlane::new(pln));
            }
        }
        return extrusion;
    }
    if let Some(circ) = curve.gp_circ() {
        let axis = *circ.axis().direction();
        if direction.is_parallel_tol(&axis, ANGULAR) {
            let pos = circ.position();
            let mut ax3 = GpAx3::from_ax2(&pos);
            if direction.dot(&axis) < 0.0 {
                ax3.z_reverse();
            }
            if let Ok(cyl) = occt_core::gp::GpCylinder::new(ax3, circ.radius()) {
                return Arc::new(GeomCylinder::new(cyl));
            }
        } else if is_normal_tol(&direction, &axis, ANGULAR) {
            if let Some(pln) = extrusion_plane(curve, direction, first, last) {
                return Arc::new(GeomPlane::new(pln));
            }
        }
        return extrusion;
    }
    let axis = curve
        .gp_ellipse()
        .map(|e| *e.axis().direction())
        .or_else(|| curve.gp_parabola().map(|e| *e.axis().direction()))
        .or_else(|| curve.gp_hyperbola().map(|e| *e.axis().direction()));
    if let Some(axis) = axis {
        if is_normal_tol(&direction, &axis, ANGULAR) {
            if let Some(pln) = extrusion_plane(curve, direction, first, last) {
                return Arc::new(GeomPlane::new(pln));
            }
        }
    }
    extrusion
}

fn is_normal_tol(a: &GpDir, b: &GpDir, tol: f64) -> bool {
    (a.angle_tol(b) - PI / 2.0).abs() <= tol
}

fn extrusion_plane(curve: &dyn Curve, direction: GpDir, mut first: f64, mut last: f64) -> Option<occt_core::gp::GpPln> {
    // `Plane` (`GeomAdaptor_SurfaceOfLinearExtrusion.cxx:332-370`).
    if Precision::is_negative_infinite(first) && Precision::is_positive_infinite(last) {
        first = -100.0;
        last = 100.0;
    } else if Precision::is_negative_infinite(first) {
        first = last - 200.0;
    } else if Precision::is_positive_infinite(last) {
        last = first + 200.0;
    }
    let deltau = (last - first) / 20.0;
    let mut point = GpPnt::new(0.0, 0.0, 0.0);
    let mut d1 = GpVec::zero();
    let mut new_z = GpVec::zero();
    for i in 1..=21 {
        let prm = first + (i - 1) as f64 * deltau;
        let sample = curve.d1(prm);
        point = sample.0;
        d1 = sample.1;
        if d1.magnitude() <= CONFUSION {
            continue;
        }
        new_z = d1.normalized().crossed(&GpVec::new(direction.x(), direction.y(), direction.z()));
        if new_z.magnitude() > 1.0e-12 {
            break;
        }
    }
    if new_z.magnitude() <= 1.0e-12 {
        return None;
    }
    let zdir = GpDir::from_vec(&new_z).ok()?;
    let xdir = GpDir::from_vec(&d1).ok()?;
    let mut ax3 = GpAx3::new(point, zdir, &xdir).ok()?;
    let y = ax3.y_direction();
    if direction.dot(y) < 0.0 {
        ax3.y_reverse();
    }
    Some(occt_core::gp::GpPln::new(ax3))
}
