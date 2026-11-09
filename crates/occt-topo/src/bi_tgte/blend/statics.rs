//! The file-local statics of `BiTgte_Blend.cxx:96-706`, in reference order:
//!
//! | function | lines |
//! |---|---|
//! | `IsOnRestriction` | `96-136` |
//! | `Add` | `140-194` |
//! | `IsInFace` | `198-211` |
//! | `KPartCurve3d` | `213-377` |
//! | `Touched` | `477-497` |
//! | `FindVertex` | `499-543` |
//! | `MakeDegeneratedEdge` | `545-584` |
//! | `Orientation` | `586-607` |
//! | `FindCreatedEdge` | `611-703` |
//!
//! Translation notes:
//! - `NCollection_Map<TopoDS_Shape>` (a set) is `HashSet<ShapeMapKey>` keyed by
//!   [`shape_map_key`], i.e. `TopTools_ShapeMapHasher` identity (same `TShape`
//!   and location, orientation ignored).
//! - `NCollection_List<TopoDS_Shape>` is `Vec<TopoShape>`; `L.First()` is
//!   `[0]` and `L.Last()` is `.last()`.
//! - `Map.Add(S)` on an `IndexedMap` is [`ShapeIndexedMap::add`]; a repeated key
//!   keeps its first index (`NCollection_IndexedMap::Add`).
//! - `BRep_Tool::CurveOnSurface(E, F, f, l)` is
//!   [`crate::boptools_2d::curve_on_surface_range`]; `BRep_Tool::Parameter(V, E,
//!   F)` is [`BRepTool::parameter_of_vertex_on_face`].
//! - `Geom2dAPI_ProjectPointOnCurve::Init(P, PC, f, l)` + `NbPoints() > 0` +
//!   `LowerDistance() < Tol` becomes the minimum distance of `P` to the p-curve
//!   restricted to `[f, l]`: the interior minimum is
//!   [`project_point_on_curve_2d`]'s answer and the two interval ends are
//!   evaluated explicitly, because the port's projector is unrestricted.
//! - `BRep_Builder::UpdateEdge(E, C, Loc, Tol)` is
//!   [`crate::brep_offset::update_edge_curve3d`] and `Degenerated(E, true)` is
//!   [`GeometryRegistry::set_degenerated`].
//! - `ElSLib::*Iso` is the `occt_core::elib::slib` family; `gp::DX2d()` /
//!   `gp::DY2d()` are built from `GpDir2d::new`.
//! - `Touched` has an empty body in the reference (`cxx:483-496` is a comment
//!   block), so the port is a no-op with the same signature.

use std::collections::HashSet;
use std::sync::Arc;

use occt_core::elib::slib;
use occt_core::gp::{GpAx1, GpCirc, GpDir2d, GpLin, GpPnt, GpPnt2d, GpVec};
use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION};
use occt_geom::{Curve, GeomCircle, GeomLine, Surface};
use occt_geom2d::Curve2d;
use occt_geom2d::geom2d_api::{project_point_on_curve_2d, Project2d};

use crate::abs::{Orientation, ShapeType};
use crate::bi_tgte::CurveOnEdge;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_offset::{ShapeIndexedMap};
use crate::brep_offset::{update_edge_curve3d, Analyse, BRepOffsetOffset};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::geom_abs_curve_type::{curve2d_kind, GeomAbsCurveType};
use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, explorer_shapes, is_equal, is_same, shape_map_key, ShapeMapKey};

use super::offset_map::OffsetMap;

/// `static bool IsOnRestriction` (`cxx:96-136`): is `v`, a vertex of `cur_e`,
/// on a restriction of `f`? Answers that restriction edge in `e`.
pub(super) fn is_on_restriction(
    v: &Vertex,
    cur_e: &Edge,
    f: &Face,
    e: &mut Option<Edge>,
) -> bool {
    // `cxx:112-114`: `CurC`, the parameter `U` of `V` on `(CurE, F)`, and
    // `P = CurC->Value(U)`. OCCT dereferences `CurC` unconditionally; the port
    // answers `false` for a missing p-curve instead of panicking.
    let Some((cur_c, _, _)) = curve_on_surface_range(cur_e, f) else {
        return false;
    };
    let u = match BRepTool::parameter_of_vertex_on_face_checked(&v.0, cur_e, f) {
        Some(u) => u,
        None => return false,
    };
    let p = cur_c.d0(u);

    // `cxx:120`: `Tol = BRep_Tool::Tolerance(V)` (`cxx:116-119` keeps the
    // exaggerated tolerance commented out).
    let tol = BRepTool::vertex_tolerance(v);

    // `cxx:121-134`.
    for occ in explorer_shapes(&f.0, ShapeType::Edge) {
        let candidate = Edge(occ.clone());
        *e = Some(candidate.clone());
        let Some((pc, pf, pl)) = curve_on_surface_range(&candidate, f) else {
            continue;
        };
        if restricted_distance(&pc, &p, pf, pl) < tol {
            return true;
        }
    }
    false
}

/// `Geom2dAPI_ProjectPointOnCurve::Init(P, PC, f, l).LowerDistance()`: the
/// minimum distance of `p` to `pc` over `[f, l]`.
fn restricted_distance(pc: &Arc<dyn Curve2d>, p: &GpPnt2d, f: f64, l: f64) -> f64 {
    let mut best = f64::INFINITY;
    if let Some(Project2d { distance, .. }) = project_point_on_curve_2d(pc.as_ref(), p, PCONFUSION)
    {
        best = best.min(distance);
    }
    for t in [f, l] {
        if t.is_finite() {
            best = best.min(pc.d0(t).distance(p));
        }
    }
    best
}

/// `static void Add` (`cxx:140-194`): record the restrictions of `s` that the
/// offset `of` maps onto `e`.
pub(super) fn add(
    e: &Edge,
    map: &mut ShapeIndexedMap,
    s: &TopoShape,
    of: &BRepOffsetOffset,
    analyse: &Analyse,
    warning_sur_bord_libre: bool,
) {
    let shape_type = s.shape_type();
    if shape_type == ShapeType::Face {
        // `cxx:151-173`.
        for occ in explorer_shapes(s, ShapeType::Edge) {
            let ori_e = Edge(occ.clone());
            let Some(ie) = of.generated(&ori_e.0) else {
                continue;
            };
            let ie = Edge(ie);
            if !is_equal(&e.0, &ie.0) {
                continue;
            }
            if warning_sur_bord_libre {
                // `cxx:161-168`: nothing is done on a free border.
                let l = analyse.ancestors(&ori_e.0);
                if l.len() == 1 {
                    break;
                }
            }
            map.add(&occ);
            break;
        }
    } else if shape_type == ShapeType::Edge {
        // `cxx:176-193`.
        for occ in explorer_shapes(s, ShapeType::Vertex) {
            let Some(ie) = of.generated(&occ) else {
                continue;
            };
            let ie = Edge(ie);
            if !is_equal(&e.0, &ie.0) {
                continue;
            }
            for anc in analyse.ancestors(&occ) {
                map.add(anc);
            }
            break;
        }
    }
}

/// `static bool IsInFace` (`cxx:198-211`): `e` is an edge of `f` by `IsSame`.
pub(super) fn is_in_face(e: &Edge, f: &Face) -> bool {
    explorer_shapes(&f.0, ShapeType::Edge)
        .iter()
        .any(|occ| is_same(occ, &e.0))
}

/// `static void KPartCurve3d` (`cxx:213-377`): give `edge` the analytic 3D
/// curve matching its p-curve `curve` on `surf` when `surf` is an elementary
/// surface, or the planar 3D curve on a plane.
pub(super) fn kpart_curve3d(edge: &Edge, curve: &Arc<dyn Curve2d>, surf: &Arc<dyn Surface>) {
    let cty = curve2d_kind(curve.as_ref());
    let sty = surface_kind(surf.as_ref());
    let tol = CONFUSION;
    let dx2d = GpDir2d::new(1.0, 0.0).unwrap();
    let dy2d = GpDir2d::new(0.0, 1.0).unwrap();

    // `cxx:231`: a plane goes through `GeomAPI::To3d` below; every other
    // surface only gets an iso when the p-curve is a line.
    if sty != GeomAbsSurfaceType::Plane {
        if cty != GeomAbsCurveType::Line {
            return;
        }
        let Some(line) = curve.gp_lin2d() else {
            return;
        };
        let d = *line.direction();
        if d.is_parallel(&dx2d, ANGULAR) {
            // `cxx:235-309`: Iso V.
            match sty {
                GeomAbsSurfaceType::Sphere => {
                    let Some(sph) = surf.gp_sphere() else {
                        return;
                    };
                    let p = line.location();
                    if (p.y().abs() - std::f64::consts::FRAC_PI_2).abs() < PCONFUSION {
                        GeometryRegistry::global().set_degenerated(&edge.0, true);
                    } else {
                        let axis = *sph.position();
                        let mut ci = slib::sphere_v_iso(&axis, sph.radius(), p.y());
                        rotate_iso(&mut ci, &axis, p.x());
                        update_with_circle(edge, &ci, &d, &dx2d, tol);
                    }
                }
                GeomAbsSurfaceType::Cylinder => {
                    let Some(cyl) = surf.gp_cylinder() else {
                        return;
                    };
                    let p = line.location();
                    let axis = cyl.position();
                    let mut ci = slib::cylinder_v_iso(&axis, cyl.radius(), p.y());
                    rotate_iso(&mut ci, &axis, p.x());
                    update_with_circle(edge, &ci, &d, &dx2d, tol);
                }
                GeomAbsSurfaceType::Cone => {
                    let Some(cone) = surf.gp_cone() else {
                        return;
                    };
                    let p = line.location();
                    let axis = cone.position();
                    let mut ci =
                        slib::cone_v_iso(&axis, cone.radius(), cone.semi_angle(), p.y());
                    rotate_iso(&mut ci, &axis, p.x());
                    update_with_circle(edge, &ci, &d, &dx2d, tol);
                }
                GeomAbsSurfaceType::Torus => {
                    let Some(tore) = surf.gp_torus() else {
                        return;
                    };
                    let p = line.location();
                    let axis = *tore.position();
                    let mut ci = slib::torus_v_iso(
                        &axis,
                        tore.major_radius(),
                        tore.minor_radius(),
                        p.y(),
                    );
                    rotate_iso(&mut ci, &axis, p.x());
                    update_with_circle(edge, &ci, &d, &dx2d, tol);
                }
                _ => {}
            }
        } else if d.is_parallel(&dy2d, ANGULAR) {
            // `cxx:310-371`: Iso U.
            match sty {
                GeomAbsSurfaceType::Sphere => {
                    let Some(sph) = surf.gp_sphere() else {
                        return;
                    };
                    let p = line.location();
                    let axis = *sph.position();
                    // `cxx:319`: iso 0 then the two rotations of `:321-328`.
                    let mut ci = slib::sphere_u_iso(&axis, sph.radius(), 0.0);
                    let mut d_rev = axis.x_direction().crossed(&axis.y_direction());
                    if let Ok(rev) = axis.x_direction().crossed(&axis.direction()) {
                        d_rev = Ok(rev);
                    }
                    if let Ok(d_rev) = d_rev {
                        let axe_rev = GpAx1::new(axis.location(), d_rev);
                        ci.rotate(&axe_rev, p.y());
                    }
                    if let Ok(d_rev) = axis.x_direction().crossed(&axis.y_direction()) {
                        let axe_rev = GpAx1::new(axis.location(), d_rev);
                        ci.rotate(&axe_rev, p.x());
                    }
                    update_with_circle(edge, &ci, &d, &dy2d, tol);
                }
                GeomAbsSurfaceType::Cylinder => {
                    let Some(cyl) = surf.gp_cylinder() else {
                        return;
                    };
                    let p = line.location();
                    let mut l = slib::cylinder_u_iso(&cyl.position(), cyl.radius(), p.x());
                    translate_line(&mut l, p.y());
                    update_with_line(edge, &l, &d, &dy2d, tol);
                }
                GeomAbsSurfaceType::Cone => {
                    let Some(cone) = surf.gp_cone() else {
                        return;
                    };
                    let p = line.location();
                    let mut l = slib::cone_u_iso(
                        &cone.position(),
                        cone.radius(),
                        cone.semi_angle(),
                        p.x(),
                    );
                    translate_line(&mut l, p.y());
                    update_with_line(edge, &l, &d, &dy2d, tol);
                }
                // `cxx:367-369`: a torus iso U is not built.
                GeomAbsSurfaceType::Torus => {}
                _ => {}
            }
        }
    } else {
        // `cxx:373-376`: `GeomAPI::To3d(Curve, Plane)`.
        let Some(pln) = surf.gp_pln() else {
            return;
        };
        let axes = pln.position().ax2();
        let Some(c3d) = occt_geom::geom_lib::to_3d(&axes, curve.as_ref()) else {
            return;
        };
        update_edge_curve3d(edge, &c3d, tol);
    }
}

/// `gp_Circ Ci` rotated by `angle` around the axis normal of `axis` through its
/// location (`cxx:250-251`, `:268-269`, `:286-287`, `:300-301`).
fn rotate_iso(ci: &mut GpCirc, axis: &occt_core::gp::GpAx3, angle: f64) {
    if let Ok(rev) = axis.x_direction().crossed(&axis.y_direction()) {
        let axe_rev = GpAx1::new(axis.location(), rev);
        ci.rotate(&axe_rev, angle);
    }
}

/// `TheBuilder.UpdateEdge(Edge, Circle, Loc, Tol)` with the optional
/// `Circle->Reverse()` of `cxx:253-257` / `:271-275` / `:288-292` / `:302-306`.
fn update_with_circle(edge: &Edge, ci: &GpCirc, d: &GpDir2d, iso_dir: &GpDir2d, tol: f64) {
    let mut circle = GeomCircle::new(ci.clone());
    if d.is_opposite(iso_dir, ANGULAR) {
        circle.reverse();
    }
    let c: Arc<dyn Curve> = Arc::new(circle);
    update_edge_curve3d(edge, &c, tol);
}

/// `TheBuilder.UpdateEdge(Edge, Line, Loc, Tol)` with `Line->Reverse()` when
/// `D.IsOpposite(gp::DY2d())` (`cxx:345-350`, `:356-361`).
fn update_with_line(edge: &Edge, l: &GpLin, d: &GpDir2d, iso_dir: &GpDir2d, tol: f64) {
    let mut line = GeomLine::new(*l);
    if d.is_opposite(iso_dir, ANGULAR) {
        line.reverse();
    }
    let c: Arc<dyn Curve> = Arc::new(line);
    update_edge_curve3d(edge, &c, tol);
}

/// `gp_Vec Tr(L.Direction()); Tr.Multiply(P.Y()); L.Translate(Tr)`
/// (`cxx:342-344`, `:357-359`).
fn translate_line(l: &mut GpLin, amount: f64) {
    let dir = l.direction();
    let v = GpVec::from_xyz(&dir.xyz()).multiply_scalar(amount);
    l.translate_vec(&v);
}

/// `static void Touched` (`cxx:477-497`): the reference body is entirely
/// commented out (`cxx:483-496`), so nothing is done.
pub(super) fn touched(
    _analyse: &Analyse,
    _stop_faces: &HashSet<ShapeMapKey>,
    _shape: &TopoShape,
    _touched_by_cork: &mut HashSet<ShapeMapKey>,
) {
}

/// `static TopoDS_Vertex FindVertex` (`cxx:499-543`): the vertex of an edge of
/// `map` that stands for `p` within its own tolerance, else within `tol` (and
/// then its tolerance is raised to `tol`).
pub(super) fn find_vertex(p: &GpPnt, map: &EdgeSet, tol: f64) -> Option<Vertex> {
    let tol_carre = tol * tol;
    for e in map.iter() {
        let (v1, v2) = edge_vertices(e);
        for vv in [v1, v2].into_iter().flatten() {
            // `cxx:520-528`: the vertex's own tolerance.
            let tol2 = BRepTool::vertex_tolerance(&vv);
            let tol2 = tol2 * tol2;
            let p1 = BRepTool::vertex_point(&vv);
            let dist = p.square_distance(&p1);
            if dist <= tol2 {
                return Some(vv);
            }
            // `cxx:529-536`: otherwise with the required tolerance.
            if tol_carre > tol2 && dist <= tol_carre {
                GeometryRegistry::global().set_vertex_tolerance(&vv.0, tol);
                return Some(vv);
            }
        }
    }
    None
}

/// `NCollection_Map<TopoDS_Shape> MapOnV` (`cxx:1930`, `:1943`): the `IsSame`
/// identity set of the edges created around one vertex of the intersection, in
/// insertion order (the order `FindVertex` and the `ComputeSurfaces`
/// re-scan of `cxx:2059-2090` walk).
#[derive(Default, Clone)]
pub(super) struct EdgeSet {
    edges: Vec<Edge>,
    keys: HashSet<ShapeMapKey>,
}

impl EdgeSet {
    pub(super) fn new() -> Self {
        Self::default()
    }

    /// `MapOnV.Add(E)`: a repeated key keeps the first entry
    /// (`NCollection_Map::Add`).
    pub(super) fn add(&mut self, e: Edge) {
        if e.0.shape_type() != ShapeType::Edge {
            return;
        }
        if self.keys.insert(shape_map_key(&e.0)) {
            self.edges.push(e);
        }
    }

    pub(super) fn iter(&self) -> impl Iterator<Item = &Edge> {
        self.edges.iter()
    }
}

/// `static TopoDS_Edge MakeDegeneratedEdge(CC, VfOnE)` (`cxx:545-584`).
pub(super) fn make_degenerated_edge(cc: &Arc<dyn Curve>, vf_on_e: Option<&Vertex>) -> Edge {
    let b = TopoBuilder::new();
    let tol = CONFUSION;
    // `cxx:554-559`: kill the trimmed curves.
    let mut basis = cc.clone();
    while let Some(inner) = basis.trim_basis() {
        basis = inner;
    }

    // `cxx:561-574`.
    let (v1, v2) = match vf_on_e {
        None => {
            let p = basis.d0(basis.first_parameter());
            let v = b.make_vertex(p, tol);
            (v.clone(), v)
        }
        Some(v) => (v.clone(), v.clone()),
    };

    // `cxx:576-582`: `B.MakeEdge(E, C, Tol)`, `B.Add(E, V1)`, `B.Add(E, V2)`
    // (FORWARD / REVERSED) and then the range of the *trimmed* curve.
    let mut e = b.make_edge(basis, cc.first_parameter(), cc.last_parameter());
    b.add_edge_vertices(&mut e, &v1, &v2);
    e
}

/// `static TopAbs_Orientation Orientation(E, F, L)` (`cxx:586-607`): the
/// orientation of `e` in the descendant list `l`, complemented when `f` is
/// REVERSED.
pub(super) fn orientation(e: &Edge, f: &Face, l: &[TopoShape]) -> Orientation {
    let mut orien = Orientation::Forward;
    for s in l {
        if is_same(s, &e.0) {
            orien = s.orientation();
            break;
        }
    }
    if f.0.orientation() == Orientation::Reversed {
        orien = orien.reversed();
    }
    orien
}

/// `static TopoDS_Edge FindCreatedEdge(V1, E, MapSF, MapOnV, CenterAnalyse,
/// Radius, Tol)` (`cxx:611-703`).
///
/// `map_on_v` is the `NCollection_Map` the reference fills; the port keeps the
/// `IsSame` identity set plus the matching edge list so that the later
/// [`find_vertex_in`] scan can read the vertices (an `NCollection_Map` stores
/// the edges themselves).
pub(super) fn find_created_edge(
    v1: &Vertex,
    e: &Edge,
    map_sf: &OffsetMap,
    map_on_v: &mut EdgeSet,
    center_analyse: &Analyse,
    radius: f64,
    tol: f64,
) -> Option<Edge> {
    // `cxx:620-624`.
    if !center_analyse.has_ancestor(&v1.0) {
        return None;
    }

    let mut tang_e: Vec<TopoShape> = Vec::new();
    center_analyse.tangent_edges(e, v1, &mut tang_e);

    // `cxx:629-683`.
    let mut e1: Option<Edge> = None;
    let mut find = false;
    let mut index = 0;
    while index < tang_e.len() && !find {
        let et = Edge(tang_e[index].clone());
        index += 1;
        if let Some(of) = map_sf.get(&et.0) {
            if let Some(gen) = of.generated(&v1.0) {
                let edge = Edge(gen);
                map_on_v.add(edge.clone());
                e1 = Some(edge);
            }
            find = true;
        } else {
            // `cxx:645-682`: the sharing of vertices in case of tangent
            // consecutive 3 edges, the second of which degenerates the tube.
            let Some(cet) = BRepTool::edge_curve(&et) else {
                continue;
            };
            let mut cet = cet;
            while let Some(inner) = cet.trim_basis() {
                cet = inner;
            }
            let Some(circ) = cet.gp_circ() else {
                continue;
            };
            if (circ.radius() - radius.abs()).abs() > tol {
                continue;
            }
            let (u1, u2) = edge_vertices(&et);
            let mut u1 = u1;
            if let (Some(a), Some(b)) = (u1.as_ref(), u2.as_ref()) {
                if is_same(&a.0, &v1.0) {
                    u1 = Some(b.clone());
                }
            }
            let Some(u1) = u1 else {
                continue;
            };
            let mut tang2: Vec<TopoShape> = Vec::new();
            center_analyse.tangent_edges(&et, &u1, &mut tang2);
            for et2 in tang2 {
                let et2 = Edge(et2);
                if let Some(of) = map_sf.get(&et2.0) {
                    if let Some(gen) = of.generated(&u1.0) {
                        map_on_v.add(Edge(gen));
                    }
                }
            }
        }
    }

    // `cxx:684-701`.
    if !find && center_analyse.has_ancestor(&v1.0) {
        for occ in center_analyse.ancestors(&v1.0) {
            if let Some(of) = map_sf.get(occ) {
                if let Some(gen) = of.generated(&v1.0) {
                    map_on_v.add(Edge(gen));
                }
            }
        }
    }

    e1
}

/// `new BiTgte_CurveOnEdge(CurCutE, OE)` used by `ComputeSurfaces`
/// (`cxx:1944`, `:1992`).
pub(super) fn curve_on_edge(eon_f: &Edge, edge: &Edge) -> CurveOnEdge {
    CurveOnEdge::new(eon_f, edge)
}




