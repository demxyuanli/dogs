//! The curve side of `BRepFill_OffsetWire`: `KPartCircle` (`cxx:192-305`) with
//! `MakeCircle` (`cxx:2130-2166`) and `MakeOffset` (`cxx:2169-2306`).
//!
//! `KPartCircle` is the shortcut taken before the bisecting locus is used: when
//! the work spine has a single edge and that edge is a closed circle (or the
//! result is open), the offset is that one edge offset along its p-curve and no
//! map is needed. `MakeOffset` does the same offsetting per spine edge for the
//! general case, and `MakeCircle` builds the arc joining two offset edges around
//! a spine vertex.
//!
//! `Adaptor2d_OffsetCurve` has no port; its two accessors are reproduced here as
//! [`offset_line`] and [`offset_circle`] (`Adaptor2d_OffsetCurve.cxx`,
//! `::Line` / `::Circle`). The offset itself is
//! [`occt_geom2d::offset::Geom2dOffsetCurve`], the faithful `Geom2d_OffsetCurve`
//! port whose normal is `(d1.y, -d1.x)` (`Geom2d_OffsetCurveUtils.pxx:30-44`);
//! `offset2d::OffsetCurve2d` is *not* used because its normal is
//! `(-d1.y, d1.x)`, i.e. the opposite side.
//!
//! Deviation: `BRepLib_MakeEdge(C2d, Plane)` attaches the p-curve to the
//! reference plane handle. The port attaches it to the work-spine face, whose
//! UV parameterization is the same plane; the two differ only when
//! `PerformWithBiLo` translates the reference plane by `Alt` (`cxx:743-751`),
//! and there the 2D coordinates are unchanged.

use std::sync::Arc;

use occt_core::gp::{GpAx2, GpAx2d, GpAx22d, GpCirc2d, GpDir2d, GpLin2d, GpTrsf, GpVec};
use occt_core::kernel::geomabs::{CurveType, JoinType};
use occt_core::precision::{CONFUSION, RESOLUTION};
use occt_core::toploc::TopLocLocation;
use occt_geom::geom_lib;
use occt_geom::Surface;
use occt_geom2d::adaptor_curve::Geom2dAdaptorCurve;
use occt_geom2d::circle::Geom2dCircle;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::offset::Geom2dOffsetCurve;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_fill_offset_wire::{BRepFillOffsetWire, IndexedShapeMap};
use crate::brep_lib_make_edge;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, edges_of, is_same, shape_map_key};

impl BRepFillOffsetWire {
    /// `KPartCircle(myWorkSpine, Offset, myIsOpenResult, Alt, myShape, myMap,
    /// myIsDone)` (`cxx:192-305`): `true` when the shortcut applied.
    pub(crate) fn kpart_circle(&mut self, offset: f64, alt: f64) -> bool {
        let Some(work_spine) = self.my_work_spine.clone() else {
            return false;
        };

        // `cxx:203-212`: exactly one edge, otherwise the shortcut is skipped.
        let edges = edges_of(&work_spine.0);
        if edges.len() != 1 {
            return false;
        }
        let e = edges[0].clone();

        // `cxx:214-222`: `BRep_Tool::Curve(E, L, f, l)`, trimmed curve unwrapped.
        let location = e.0.location().clone();
        let Some(mut c) = BRepTool::edge_curve(&e) else {
            return false;
        };
        if let Some(basis) = c.trim_basis() {
            c = basis;
        }

        // `cxx:229-230`: a closed circle, or an open result.
        let closed_circle = c.gp_circ().is_some() && BRepTool::is_closed_edge(&e);
        if !(closed_circle || self.my_is_open_result) {
            return false;
        }

        // `cxx:234-235`: the p-curve range replaces the 3D range.
        let Some((a_pcurve, first, last)) = curve_on_surface_range(&e, &work_spine) else {
            panic!("BRep_Tool::CurveOnSurface(E, mySpine) raised");
        };
        let ahc = Geom2dAdaptorCurve::with_range(a_pcurve.clone(), first, last);

        let mut an_offset = offset;
        let oc: Arc<dyn Curve2d> = match ahc.get_type() {
            Some(CurveType::Line) => {
                // `cxx:239-245`.
                if e.0.orientation() == Orientation::Forward {
                    an_offset *= -1.0;
                }
                Arc::new(Geom2dLine::new(*offset_line(&ahc, an_offset).position()))
            }
            Some(CurveType::Circle) => {
                // `cxx:249-264`.
                if e.0.orientation() == Orientation::Forward {
                    an_offset *= -1.0;
                }
                if !BRepTool::is_closed_edge(&e) {
                    an_offset *= -1.0;
                }
                let the_circ = ahc
                    .circle()
                    .expect("Geom2dAdaptor_Curve::Circle on a circle");
                if an_offset > 0.0 || an_offset.abs() < the_circ.radius() {
                    Arc::new(Geom2dCircle::new(GpCirc2d::new(
                        *the_circ.position(),
                        the_circ.radius() + an_offset,
                    )))
                } else {
                    // `cxx:261-264`.
                    self.my_is_done = false;
                    return false;
                }
            }
            _ => {
                // `cxx:267-275`.
                if e.0.orientation() == Orientation::Forward {
                    an_offset *= -1.0;
                }
                let g2dt = Geom2dTrimmedCurve::new(a_pcurve.clone(), first, last);
                Arc::new(Geom2dOffsetCurve::new(Arc::new(g2dt), an_offset))
            }
        };

        // `cxx:277-279`: `BRepLib_MakeEdge(OC, aPlane, f, l)` + `BuildCurve3d`.
        let builder = TopoBuilder::new();
        let surface = BRepTool::face_surface(&work_spine).expect("BRep_Tool::Surface");
        let mut shape = make_edge_2d(&builder, &oc, &surface, &work_spine, first, last)
            .expect("BRepLib_MakeEdge(OC, aPlane, f, l)");

        // `cxx:281-282`.
        shape.0.set_orientation(e.0.orientation());
        shape.0.set_location(&location);

        // `cxx:283-293`: the altitude translation along the plane normal.
        if alt.abs() > RESOLUTION {
            let pln = surface.gp_pln().expect("BRepAdaptor_Surface.Plane()");
            let axis = pln.axis();
            let dir = axis.direction();
            let trans = GpVec::new(alt * dir.x(), alt * dir.y(), alt * dir.z());
            let mut trsf = GpTrsf::identity();
            trsf.set_translation_vec(&trans);
            shape
                .0
                .set_location(&TopLocLocation::composed(&trsf, TopLocLocation::identity()));
        }

        // `cxx:295-300`.
        self.my_map.add(&e.0, vec![shape.0.clone()]);
        let wire = builder.make_wire(&[Edge(shape.0.clone())]);
        self.my_shape = Some(wire.0);
        self.my_is_done = true;
        true
    }
}

/// `MakeCircle(E, V, F, Offset, Map, RefPlane)` (`cxx:2130-2166`): the arc of
/// radius `|Offset|` joining the two offset edges at the spine vertex `v`.
pub(crate) fn make_circle(
    e: &Edge,
    v: &Vertex,
    f: &Face,
    offset: f64,
    map: &mut IndexedShapeMap,
    ref_plane: &Arc<dyn Surface>,
) {
    // `cxx:2139-2150`: `D1` at the far end of the p-curve, reversed on a
    // FORWARD edge.
    let mut gc = curve_on_surface_range(e, f)
        .map(|(c, _, _)| c)
        .unwrap_or_else(|| panic!("BRep_Tool::CurveOnSurface(E, F) raised"));
    if let Some(basis) = gc.trimmed_basis() {
        gc = Arc::from(basis.clone_dyn());
    }
    let (pf, pl) = BRepTool::edge_parameters(e);
    let (p, mut dx) = if e.0.orientation() == Orientation::Forward {
        gc.d1(pl)
    } else {
        gc.d1(pf)
    };
    if e.0.orientation() == Orientation::Forward {
        dx = occt_core::gp::GpVec2d::new(-dx.x(), -dx.y());
    }

    // `cxx:2152-2153`: `Geom2d_Circle(gp_Ax2d(P, gp_Dir2d(DX)), |Offset|,
    // Offset < 0)`.
    let axis = GpAx2d::new(p, GpDir2d::new(dx.x(), dx.y()).expect("gp_Dir2d(DX)"));
    let circ = Geom2dCircle::new(GpCirc2d::new(
        circle_axes(&axis, offset < 0.0),
        offset.abs(),
    ));
    let circ: Arc<dyn Curve2d> = Arc::new(circ);
    let (first, last) = (circ.first_parameter(), circ.last_parameter());

    // `cxx:2156-2161`.
    let builder = TopoBuilder::new();
    let oe = make_edge_2d(&builder, &circ, ref_plane, f, first, last)
        .expect("BRepLib_MakeEdge(Circ, RefPlane)");
    map.add(&v.0, vec![oe.0.clone()]);
}

/// `MakeOffset(E, F, Offset, Map, RefPlane, IsOpenResult, theJoinType, Ends)`
/// (`cxx:2169-2306`): the parallel of one spine edge.
pub(crate) fn make_offset(
    e: &Edge,
    f: &Face,
    offset: f64,
    map: &mut IndexedShapeMap,
    ref_plane: &Arc<dyn Surface>,
    is_open_result: bool,
    the_join_type: JoinType,
    ends: &[Option<Vertex>; 2],
) {
    // `cxx:2178-2183`.
    let mut an_offset = offset;
    if e.0.orientation() == Orientation::Forward {
        an_offset *= -1.0;
    }

    let (g2d, mut first, mut last) = curve_on_surface_range(e, f)
        .unwrap_or_else(|| panic!("BRep_Tool::CurveOnSurface(E, F) raised"));

    // `cxx:2189-2203`: an open result leaves the two end edges unextended.
    let mut to_extend_first = true;
    let mut to_extend_last = true;
    if is_open_result {
        let (v1, v2) = edge_vertices(e);
        let same = |v: &Option<Vertex>, w: &Option<Vertex>| match (v, w) {
            (Some(a), Some(b)) => is_same(&a.0, &b.0),
            _ => false,
        };
        if same(&v1, &ends[0]) || same(&v1, &ends[1]) {
            to_extend_first = false;
        }
        if same(&v2, &ends[0]) || same(&v2, &ends[1]) {
            to_extend_last = false;
        }
    }

    let ac = Geom2dAdaptorCurve::with_range(g2d.clone(), first, last);
    let g2d_oc: Option<Arc<dyn Curve2d>> = match ac.get_type() {
        Some(CurveType::Circle) => {
            // `cxx:2207-2256`.
            let c1 = ac.circle().expect("Geom2dAdaptor_Curve::Circle");
            let axes = c1.position();
            let crossed = axes.vxdir.x * axes.vydir.y - axes.vxdir.y * axes.vydir.x;
            let signe = if crossed > 0.0 { -1.0 } else { 1.0 };

            if an_offset * signe < c1.radius() - CONFUSION {
                let cc =
                    offset_circle(&ac, an_offset).expect("Adaptor2d_OffsetCurve::Circle raised");
                let delta = 2.0 * std::f64::consts::PI - last + first;
                if the_join_type == JoinType::Arc {
                    if to_extend_first {
                        first -= 0.2 * delta;
                    }
                    if to_extend_last {
                        last += 0.2 * delta;
                    }
                } else if to_extend_first && to_extend_last {
                    let old_l = last;
                    first = old_l + delta / 2.0;
                    last = first + 2.0 * std::f64::consts::PI;
                } else if to_extend_first {
                    first = last;
                    last = first + 2.0 * std::f64::consts::PI;
                } else if to_extend_last {
                    last = first + 2.0 * std::f64::consts::PI;
                }
                Some(Arc::new(Geom2dTrimmedCurve::new(
                    Arc::new(Geom2dCircle::new(cc)),
                    first,
                    last,
                )))
            } else {
                // `cxx:2209`: the offset reaches past the radius, no curve.
                None
            }
        }
        Some(CurveType::Line) => {
            // `cxx:2258-2285`.
            let line = offset_line(&ac, an_offset);
            let delta = last - first;
            if to_extend_first {
                if the_join_type == JoinType::Arc {
                    first -= delta;
                } else {
                    first = f64::NEG_INFINITY;
                }
            }
            if to_extend_last {
                if the_join_type == JoinType::Arc {
                    last += delta;
                } else {
                    last = f64::INFINITY;
                }
            }
            Some(Arc::new(Geom2dTrimmedCurve::new(
                Arc::new(Geom2dLine::new(*line.position())),
                first,
                last,
            )))
        }
        _ => {
            // `cxx:2287-2293`.
            let g2dt = Geom2dTrimmedCurve::new(g2d.clone(), first, last);
            Some(Arc::new(Geom2dOffsetCurve::new(
                Arc::new(g2dt),
                an_offset,
            )))
        }
    };

    // `cxx:2295-2303`.
    if let Some(g2d_oc) = g2d_oc {
        let builder = TopoBuilder::new();
        let (first, last) = (g2d_oc.first_parameter(), g2d_oc.last_parameter());
        let mut oe = make_edge_2d(&builder, &g2d_oc, ref_plane, f, first, last)
            .expect("BRepLib_MakeEdge(G2dOC, RefPlane)");
        oe.0.set_orientation(e.0.orientation());
        map.add(&e.0, vec![oe.0.clone()]);
    }
}

/// `gp_Ax22d(P, D, Sense)` (`gp_Ax22d.hxx`): `Y = ±rot90(X)`, `Sense = false`
/// reversing it.
fn circle_axes(axis: &GpAx2d, sense: bool) -> GpAx22d {
    let vx = axis.vdir;
    let mut vy = GpDir2d::new(-vx.y, vx.x).expect("gp_Ax22d Y direction");
    if !sense {
        vy = vy.reversed();
    }
    GpAx22d::new(axis.loc, vx, vy).expect("gp_Ax22d perpendicular directions")
}

/// `Adaptor2d_OffsetCurve::Line()`: `D1(0)` of the offset curve as a line.
fn offset_line(ac: &Geom2dAdaptorCurve, offset: f64) -> GpLin2d {
    let basis = ac.curve().expect("empty adaptor curve").clone();
    let off = Geom2dOffsetCurve::new(basis, offset);
    let (p, v) = off.d1(0.0);
    GpLin2d::from_pnt_dir(p, GpDir2d::new(v.x(), v.y()).expect("gp_Lin2d(P, V)"))
}

/// `Adaptor2d_OffsetCurve::Circle()`: the concentric circle, axis flipped when
/// the radius would go negative.
fn offset_circle(ac: &Geom2dAdaptorCurve, offset: f64) -> Option<GpCirc2d> {
    let c1 = ac.circle()?;
    if offset == 0.0 {
        return Some(c1);
    }
    let axes = c1.position();
    let mut radius = c1.radius();
    let crossed = axes.vxdir.x * axes.vydir.y - axes.vxdir.y * axes.vydir.x;
    let signe = if crossed > 0.0 { 1.0 } else { -1.0 };
    radius += signe * offset;
    if radius > 0.0 {
        Some(GpCirc2d::new(*axes, radius))
    } else if radius < 0.0 {
        let ax = GpAx22d { point: axes.point, vxdir: axes.vxdir.reversed(), vydir: axes.vydir };
        Some(GpCirc2d::new(ax, -radius))
    } else {
        None
    }
}

/// `BRepLib_MakeEdge(C2d, Plane, f, l)` (`BRepLib_MakeEdge.cxx:905-1094`)
/// followed by `BRepLib::BuildCurve3d`: the edge carries the 3D lift of the
/// p-curve (`GeomLib::To3d`) and the p-curve itself on `face`.
fn make_edge_2d(
    b: &TopoBuilder,
    g2d: &Arc<dyn Curve2d>,
    surface: &Arc<dyn Surface>,
    face: &Face,
    first: f64,
    last: f64,
) -> Option<Edge> {
    let pln = surface.gp_pln()?;
    let axes: GpAx2 = pln.position().ax2();
    let c3d = geom_lib::to_3d(&axes, g2d.as_ref())?;
    let mut edge = b.make_edge(c3d, first, last);
    let fk = shape_map_key(&face.0).0;
    let reg = GeometryRegistry::global();
    reg.set_edge_pcurve(&edge.0, fk, g2d.clone());
    reg.set_pcurve_range(&edge.0, fk, first, last);
    if let Some(build) = brep_lib_make_edge::init_2d(g2d, surface, first, last) {
        brep_lib_make_edge::add_endpoint_vertices(b, &mut edge, &build);
    }
    Some(edge)
}
