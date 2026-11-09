//! `BRepOffset_Offset::Init(Vertex, LEdge, Offset, Polynomial, TolApp, Conti)`
//! (`BRepOffset_Offset.cxx:1559-1696`): the sphere a vertex of the source
//! shape turns into when its incident edges are all offset with an arc join.
//!
//! The axis of the sphere is built the way the `.cxx` builds it:
//! `:1590` the vertex point is the origin, `:1594-1598` the incident edges are
//! assembled into a wire and run through `ShapeFix_Shape`, `:1601-1603`
//! `BRepGProp::LinearProperties` gives the wire's centre of mass (the X
//! direction points from it to the vertex), and `:1605-1610`
//! `GetFarestCorner` + `gce_MakePln(Origin, BaryCenter, FarestCorner)` give the
//! Z direction of the `gp_Ax3`.
//!
//! UNPORTED, with the OCCT lines that keep it out:
//!
//! * the `Polynomial` arm (`:1622-1633`): `GeomConvert_ApproxSurface` has no
//!   counterpart, so `SS` stays the exact `Geom_SphericalSurface`. Every
//!   `BRepOffset_MakeOffset` call site passes the default `Polynomial = false`.
//! * `OCCT_DEBUG` (`:1577-1587`): debug output only.
//! * `BRepTools::Update(myFace)` (`:1695`): the port keeps no cached UV points
//!   nor triangulation.

use std::sync::Arc;

use occt_core::gp::{GpAx2d, GpAx3, GpDir, GpDir2d, GpPnt, GpPnt2d, GpSphere, GpTrsf2d, GpVec, GpVec2d};
use occt_core::gprop::gprops::PEquation;
use occt_core::kernel::geomabs::Shape;
use occt_geom::geom_proj_lib;
use occt_geom::{Curve, GeomSphere, GeomTrimmedCurveBasis, Surface};
use occt_geom2d::Curve2d;

use crate::abs::Orientation;
use crate::brep_gprop_full::linear_properties;
use crate::brep_lib_build_curves3d::build_curve3d;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;

use crate::shape::{Edge, Vertex, Wire};
use crate::shhealing::ShapeFixShape;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::vertices_of;

use super::compute_curve3d::update_edge_pcurve;
use super::{BRepOffsetOffset, BRepOffsetStatus};

/// `gp::DX2d()` (`gp.hxx`).
const DX2D: GpDir2d = GpDir2d { x: 1.0, y: 0.0 };

/// `GetFarestCorner(const TopoDS_Wire&)` (`cxx:87-111`): the point of the pair
/// of wire vertices with the largest squared distance. `thePoint` stays the
/// default-constructed `gp_Pnt` when every distance is zero, so this answers
/// the origin then, exactly like the `.cxx`.
fn farest_corner(w: &Wire) -> GpPnt {
    let verts = vertices_of(&w.0);
    let mut max_dist = 0.0;
    let mut the_point = GpPnt::zero();
    for v1 in &verts {
        for v2 in &verts {
            let p1 = BRepTool::vertex_point_world(v1);
            let p2 = BRepTool::vertex_point_world(v2);
            let d = p1.square_distance(&p2);
            if d > max_dist {
                max_dist = d;
                the_point = p1;
            }
        }
    }
    the_point
}

impl BRepOffsetOffset {
    /// `BRepOffset_Offset::Init(Vertex, LEdge, Offset, Polynomial, TolApp,
    /// Conti)` (`cxx:1559-1696`).
    ///
    /// `l_edge` is OCCT's `LEdge` list: the images of the edges incident to
    /// `vertex` on the adjacent offset faces, already reversed the way
    /// `BRepOffset_MakeOffset::BuildOffsetByArc` collects them (`:2036`).
    pub fn init_vertex(
        &mut self,
        vertex: &Vertex,
        l_edge: &[Edge],
        offset: f64,
        polynomial: bool,
        tol_app: f64,
        conti: Shape,
    ) -> Result<(), String> {
        self.my_status = BRepOffsetStatus::Good;
        self.my_shape = Some(vertex.0.clone());

        // `cxx:1589`: `Origin = BRep_Tool::Pnt(Vertex)`.
        let origin = BRepTool::vertex_point_world(vertex);

        // `cxx:1594-1598`: `BRepLib_MakeWire` on the edge list, then
        // `ShapeFix_Shape`.
        let mut mw = MakeWire::new();
        for e in l_edge {
            let _ = mw.add(e);
        }
        let raw_wire = mw.wire();
        let mut fixer = ShapeFixShape::new(&raw_wire.0);
        fixer.perform();
        let the_wire = Wire(fixer.shape().clone());

        // `cxx:1601-1605`: `BRepGProp::LinearProperties(theWire, GlobalProps)`
        // and `Xdir = gp_Vec(BaryCenter, Origin)`.
        let props = linear_properties(&the_wire.0)?;
        let bary_center = props.center();
        let xdir = GpVec::from_xyz(&origin.coord.subtracted(&bary_center.coord));

        // `cxx:1606-1610`: the plane through the three points; its normal is
        // the Z direction of the sphere axis.
        let farest = farest_corner(&the_wire);
        let plane = PEquation::from_points(&origin, &bary_center, &farest)?;
        let vdir: GpDir = plane.to_plane().position().direction();
        let axis = GpAx3::new(origin, vdir, &xdir_dir(&xdir))?;

        // `cxx:1612`: `Geom_SphericalSurface(Axis, |Offset|)`.
        let sphere = GpSphere::new(axis, offset.abs())
            .map_err(|e| format!("Geom_SphericalSurface: {e}"))?;
        let surf: Arc<dyn Surface> = Arc::new(GeomSphere::new(sphere));

        // `cxx:1614-1633`: `SS = S`; the polynomial approximation arm has no
        // counterpart (see the module note), so `SS` always stays exact.
        let _ = (polynomial, tol_app, conti);

        let tol = BRepTool::vertex_tolerance(vertex);
        let builder = TopoBuilder::new();
        let mut my_face = builder.make_face_with_tolerance(surf.clone(), &[], tol);

        // `cxx:1636-1689`: project every incident edge onto the sphere and
        // accumulate them into the wire.
        let mut w = Wire::new();
        for e in l_edge {
            // `cxx:1642-1649`: a missing 3D curve is built first.
            if BRepTool::edge_curve(e).is_none() {
                let _ = build_curve3d(e, BRepTool::edge_tolerance(e), Shape::C1, 14, 0);
            }
            let (f, l) = BRepTool::edge_parameters(e);
            let Some(c) = BRepTool::edge_curve(e) else {
                return Err("BRepOffset_Offset::Init(Vertex): an edge has no curve".into());
            };
            let trimmed: Arc<dyn Curve> = Arc::new(GeomTrimmedCurveBasis::new(c, f, l));
            let c_world: Arc<dyn Curve> = if e.0.location().is_identity() {
                trimmed
            } else {
                Arc::from(trimmed.transformed(&e.0.location().transformation()))
            };

            // `cxx:1651`: `GeomProjLib::Curve2d(C, S)`.
            let pcurve = geom_proj_lib::curve_2d(&c_world, &surf)
                .ok_or("GeomProjLib::Curve2d: no projection")?;
            let pcurve = adjust_to_sphere_limits(pcurve, f);

            // `cxx:1680-1682`: `UpdateEdge` + `Range` + `Add`.
            update_edge_pcurve(e, &pcurve, &my_face, tol);
            GeometryRegistry::global().set_pcurve_range(
                &e.0,
                GeometryRegistry::shape_key(&my_face.0),
                f,
                l,
            );
            builder.add_edge(&mut w, e);
        }

        // `cxx:1684-1692`: a negative offset reverses the wire, or the face.
        if offset < 0.0 {
            builder.add_wire(&mut my_face, &Wire(w.0.oriented(Orientation::Reversed)));
            my_face.0.reverse();
        } else {
            builder.add_wire(&mut my_face, &w);
        }

        self.my_face = Some(my_face);
        Ok(())
    }
}

/// `gp_Dir(Xdir)` (`cxx:1610`): the X direction of the sphere axis, `gp::DX()`
/// when `Xdir` has no modulus (the origin coincides with the barycentre).
fn xdir_dir(v: &GpVec) -> GpDir {
    GpDir::from_xyz(&v.coord).unwrap_or(GpDir::from_axis(occt_core::gp::dir::DirAxis::X))
}

/// `cxx:1652-1680`: move the projected curve into the sphere's canonical
/// parameter range. `U' = U + PI + 2kPI`, `V' = +/-PI + 2k'PI`: the mirror is
/// about the horizontal line `V = +/- PI/2`, then the curve is translated by
/// `+/- PI` in `U` so its first point falls back into `[0, 2PI]`.
fn adjust_to_sphere_limits(curve: Arc<dyn Curve2d>, f: f64) -> Arc<dyn Curve2d> {
    let p2d: GpPnt2d = curve.value(f);
    let mirror_y = if p2d.y() < -std::f64::consts::FRAC_PI_2 {
        -std::f64::consts::FRAC_PI_2
    } else if p2d.y() > std::f64::consts::FRAC_PI_2 {
        std::f64::consts::FRAC_PI_2
    } else {
        return curve;
    };

    let mut c = curve.clone_dyn();
    let mut t = GpTrsf2d::identity();
    let ax = GpAx2d::new(GpPnt2d::new(0.0, mirror_y), DX2D);
    t.set_mirror_ax2d(&ax);
    c.transform(&t);

    // `cxx:1669-1677`: `Tr(M_PI, 0)`, reversed when the first U is above PI.
    let mut tr = GpVec2d::new(std::f64::consts::PI, 0.0);
    if p2d.x() > std::f64::consts::PI {
        tr.reverse();
    }
    let mut t2 = GpTrsf2d::identity();
    t2.set_translation_vec(&tr);
    c.transform(&t2);
    Arc::from(c)
}
