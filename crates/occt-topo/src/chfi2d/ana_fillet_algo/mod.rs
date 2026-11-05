//! `ChFi2d_AnaFilletAlgo` — an analytical fillet between two planar edges.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi2d/ChFi2d_AnaFilletAlgo.hxx:24-142`
//! / `.cxx:1-1101` (OCCT 8.0.0). Segments and arcs of circle only.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `ChFi2d_AnaFilletAlgo()` (`cxx:82-104`) | 23 | [`AnaFilletAlgo::new`] |
//! | `ChFi2d_AnaFilletAlgo(wire, plane)` (`cxx:110-131`) | 22 | [`AnaFilletAlgo::from_wire`] |
//! | `ChFi2d_AnaFilletAlgo(e1, e2, plane)` (`cxx:137-161`) | 25 | [`AnaFilletAlgo::from_edges`] |
//! | `Init(wire, plane)` (`cxx:164-255`) | 92 | [`AnaFilletAlgo::init_wire`] |
//! | `Init(e1, e2, plane)` (`cxx:257-325`) | 69 | [`AnaFilletAlgo::init_edges`] |
//! | `Perform(radius)` (`cxx:326-571`) | 246 | [`AnaFilletAlgo::perform`] |
//! | `Result(e1, e2)` (`cxx:573-579`) | 7 | [`AnaFilletAlgo::result`] |
//!
//! The four analytical cases live in [`fillet_cases`], `isCW` / `IsEqual` /
//! `Cut` in [`fillet_support`].
//!
//! Translation notes:
//! - `BRepAdaptor_Curve` is the stored curve through
//!   [`crate::brep_lib_make_edge::adaptor_curve`] (untrimmed basis plus the
//!   interval in basis parameters); `GetType()` is [`Curve::is_line`] /
//!   [`Curve::gp_circ`],
//! - `TopExp::Vertices(E, v1, v2, cumOri = true)` is
//!   [`crate::shhealing::first_vertex`] / [`crate::shhealing::last_vertex`],
//! - `ProjLib::Project(gp_Pln, gp_Pnt)` is
//!   [`occt_geom::projlib::project_pln_pnt`], `ElSLib::Value` is
//!   [`occt_core::elib::slib::plane_value`],
//! - `BRepBuilderAPI_MakeWire` is [`MakeWire`], `BRepBuilderAPI_MakeFace(plane)`
//!   is [`TopoBuilder::make_face_plane`], `ShapeAnalysis_Wire::
//!   CheckSelfIntersection` is [`crate::shhealing::check_self_intersection`],
//! - a null result edge (`fillet` / `shrinke1` / `shrinke2`) is `None`; OCCT
//!   sets all three or none of them, so [`AnaFilletAlgo::result`] answers one
//!   `Option`,
//! - the `Standard_TypeMismatch` / `Standard_Failure` arms of the two `Init`
//!   are the `Err(String)` of this port, with the OCCT message text.

mod fillet_cases;
mod fillet_support;

use std::sync::Arc;

use occt_core::elib::slib;
use occt_core::gp::{GpAx2, GpCirc, GpPln};
use occt_core::precision::CONFUSION;
use occt_geom::circle::GeomCircle;
use occt_geom::curve::Curve;

use crate::abs::Orientation;
use crate::brep_lib_make_edge;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::iterator::ShapeIterator;
use crate::shape::{Edge, Wire};
use crate::shhealing::{check_self_intersection, first_vertex, last_vertex};

/// `myLeft neighbour` (`hxx:112-123`) and `myRight neighbour`
/// (`hxx:126-136`): `e1, segment1, x11, y11, x12, y12, xc1, yc1, radius1, cw1`.
#[derive(Clone, Default)]
struct Neighbour {
    /// `TopoDS_Edge e1;` — null until `Init`.
    edge: Option<Edge>,
    /// `bool segment1;` — `true` for a line, `false` for an arc of circle.
    segment: bool,
    /// `x11, y11` — the projected first vertex (`cxx:206-207`).
    x1: f64,
    y1: f64,
    /// `x12, y12` — the projected last vertex (`cxx:207-208`), overwritten by
    /// `Cut` (`cxx:1094-1097`).
    x2: f64,
    y2: f64,
    /// `xc1, yc1` — the projected circle centre (`cxx:213-214`), `0` on a line.
    xc: f64,
    yc: f64,
    /// `radius1` (`cxx:216`), `0` on a line.
    radius: f64,
    /// `cw1 = isCW(AC1)` (`cxx:217`).
    cw: bool,
}

impl Neighbour {
    /// The left (`cxx:189-219`) and right (`cxx:221-253`) neighbour blocks.
    fn init(edge: Edge, plane: &GpPln) -> Result<Self, String> {
        // `BRepAdaptor_Curve AC(e);` — a null curve handle is dereferenced by
        // the `GetType()` call in OCCT.
        let stored = BRepTool::edge_curve(&edge)
            .ok_or_else(|| "ChFi2d_AnaFilletAlgo: the edge has no 3D curve".to_string())?;
        if !stored.is_line() && stored.gp_circ().is_none() {
            return Err("A segment or an arc of circle is expected.".to_string());
        }
        let (curve, first, last) = brep_lib_make_edge::adaptor_curve(&stored);

        // `TopExp::Vertices(e, v1, v2, true); if (v1.IsNull() || v2.IsNull())`
        // (`cxx:193-199`).
        let (Some(v1), Some(v2)) = (first_vertex(&edge), last_vertex(&edge)) else {
            return Err("An infinite edge.".to_string());
        };
        // `cxx:201-204`: both vertices projected on the plane.
        let p1 = BRepTool::vertex_point(&v1);
        let p2 = BRepTool::vertex_point(&v2);
        let q1 = occt_geom::projlib::project_pln_pnt(plane, &p1);
        let q2 = occt_geom::projlib::project_pln_pnt(plane, &p2);

        let mut neighbour = Neighbour {
            edge: Some(edge),
            segment: true,
            x1: q1.x(),
            y1: q1.y(),
            x2: q2.x(),
            y2: q2.y(),
            ..Neighbour::default()
        };

        // `cxx:209-218`: `if (AC.GetType() == GeomAbs_Circle)`.
        if let Some(circle) = stored.gp_circ() {
            neighbour.segment = false;
            let loc = occt_geom::projlib::project_pln_pnt(plane, &circle.location());
            neighbour.xc = loc.x();
            neighbour.yc = loc.y();
            neighbour.radius = circle.radius();
            neighbour.cw = fillet_support::is_cw(&curve, first, last);
        }
        Ok(neighbour)
    }

    /// `BRepAdaptor_Curve(edge)` on the edge as it stands now — re-read after
    /// `Cut` replaced it (`cxx:371`).
    fn curve_range(&self) -> Option<(Arc<dyn Curve>, f64, f64)> {
        fillet_support::curve_range(self.edge.as_ref()?)
    }
}

/// The five out-parameters of the four analytical cases (`cxx:338-343`):
/// `xc, yc`, `start`, `end`, `cw`, plus the two optional points
/// `(xstart, ystart)` / `(xend, yend)` that `DBL_MAX` marks as "not computed".
#[derive(Clone)]
struct FilletDef {
    xc: f64,
    yc: f64,
    cw: bool,
    /// Parameter on the left neighbour.
    start: f64,
    /// Parameter on the right neighbour.
    end: f64,
    xstart: f64,
    ystart: f64,
    xend: f64,
    yend: f64,
}

impl Default for FilletDef {
    fn default() -> Self {
        Self {
            xc: 0.0,
            yc: 0.0,
            cw: false,
            start: 0.0,
            end: 0.0,
            xstart: f64::MAX,
            ystart: f64::MAX,
            xend: f64::MAX,
            yend: f64::MAX,
        }
    }
}

/// The three edges [`AnaFilletAlgo::result`] answers: `fillet`, `shrinke1`,
/// `shrinke2` (`hxx:139-141`).
#[derive(Clone)]
pub struct AnaFilletResult {
    /// `fillet` — the arc of the fillet.
    pub fillet: Edge,
    /// `shrinke1` — the left neighbour limited to the fillet start.
    pub shrinke1: Edge,
    /// `shrinke2` — the right neighbour limited to the fillet end.
    pub shrinke2: Edge,
}

/// `ChFi2d_AnaFilletAlgo` (`hxx:24-142`).
#[derive(Clone, Default)]
pub struct AnaFilletAlgo {
    /// `gp_Pln plane;` (`hxx:111`).
    plane: GpPln,
    /// `TopoDS_Edge e1; ...` (`hxx:114-123`).
    left: Neighbour,
    /// `TopoDS_Edge e2; ...` (`hxx:127-136`).
    right: Neighbour,
    /// `TopoDS_Edge fillet;` (`hxx:139`).
    fillet: Option<Edge>,
    /// `TopoDS_Edge shrinke1;` (`hxx:140`).
    shrinke1: Option<Edge>,
    /// `TopoDS_Edge shrinke2;` (`hxx:141`).
    shrinke2: Option<Edge>,
}

impl AnaFilletAlgo {
    /// `ChFi2d_AnaFilletAlgo()` (`cxx:82-104`): every neighbour field at zero,
    /// both edges and the three results null.
    pub fn new() -> Self {
        Self::default()
    }

    /// `ChFi2d_AnaFilletAlgo(theWire, thePlane)` (`cxx:110-131`).
    pub fn from_wire(wire: &Wire, plane: &GpPln) -> Result<Self, String> {
        let mut algo = Self::new();
        algo.init_wire(wire, plane)?;
        Ok(algo)
    }

    /// `ChFi2d_AnaFilletAlgo(theEdge1, theEdge2, thePlane)` (`cxx:137-161`).
    pub fn from_edges(edge1: &Edge, edge2: &Edge, plane: &GpPln) -> Result<Self, String> {
        let mut algo = Self::new();
        algo.init_edges(edge1, edge2, plane)?;
        Ok(algo)
    }

    /// `Init(theWire, thePlane)` (`cxx:164-255`): the first two wire children
    /// (only into a still-null member, exactly like the `e1.IsNull()` chain),
    /// then the two neighbour blocks.
    pub fn init_wire(&mut self, wire: &Wire, plane: &GpPln) -> Result<(), String> {
        self.plane = plane.clone();
        // `TopoDS_Iterator itr(theWire); for (; itr.More(); itr.Next())`
        // (`cxx:169-180`): `TopoDS::Edge` on a non-edge child raises
        // `Standard_TypeMismatch`.
        for child in ShapeIterator::of_shape(&wire.0) {
            let edge = Edge::wrap(child)
                .ok_or_else(|| "ChFi2d_AnaFilletAlgo::Init: wire child is not an edge".to_string())?;
            if self.left.edge.is_none() {
                self.left.edge = Some(edge);
            } else if self.right.edge.is_none() {
                self.right.edge = Some(edge);
            }
        }
        // `cxx:181-185`.
        if self.left.edge.is_none() || self.right.edge.is_none() {
            return Err(
                "The algorithm expects a wire consisting of two linear or circular edges."
                    .to_string(),
            );
        }

        // `cxx:188-219`: the left neighbour.
        let left = Neighbour::init(self.left.edge.clone().unwrap(), &self.plane)?;
        self.left = left;
        // `cxx:221-252`: the right neighbour.
        let right = Neighbour::init(self.right.edge.clone().unwrap(), &self.plane)?;
        self.right = right;
        Ok(())
    }

    /// `Init(theEdge1, theEdge2, thePlane)` (`cxx:257-325`): finds the common
    /// point, orients both edges into a wire and hands that wire to
    /// [`AnaFilletAlgo::init_wire`].
    pub fn init_edges(&mut self, edge1: &Edge, edge2: &Edge, plane: &GpPln) -> Result<(), String> {
        // `cxx:263-269`.
        let (v11, v12) = (first_vertex(edge1), last_vertex(edge1));
        let (v21, v22) = (first_vertex(edge2), last_vertex(edge2));
        let (Some(v11), Some(v12), Some(v21), Some(v22)) = (v11, v12, v21, v22) else {
            return Err("An infinite edge.".to_string());
        };
        let p11 = BRepTool::vertex_point(&v11);
        let p12 = BRepTool::vertex_point(&v12);
        let p21 = BRepTool::vertex_point(&v21);
        let p22 = BRepTool::vertex_point(&v22);

        // `cxx:277-288`.
        let pcommon = if fillet_support::is_equal(&p11, &p21) || fillet_support::is_equal(&p11, &p22)
        {
            p11
        } else if fillet_support::is_equal(&p12, &p21) || fillet_support::is_equal(&p12, &p22) {
            p12
        } else {
            return Err("The edges have no common point.".to_string());
        };

        // `cxx:291-299`.
        let mut is1st_reversed = false;
        let mut is2nd_reversed = false;
        if fillet_support::is_equal(&pcommon, &p11) {
            is1st_reversed = true;
        } else if fillet_support::is_equal(&pcommon, &p22) {
            is2nd_reversed = true;
        }

        // `cxx:302-323`: `BRepBuilderAPI_MakeWire mkWire;` with the two `Add`
        // calls, then `IsDone()`.
        let mut mkwire = MakeWire::new();
        let first = if is1st_reversed {
            reversed(edge1)
        } else {
            edge1.clone()
        };
        let _ = mkwire.add(&first);
        let second = if is2nd_reversed {
            reversed(edge2)
        } else {
            edge2.clone()
        };
        let _ = mkwire.add(&second);
        if !mkwire.is_done() {
            return Err("Can't make a wire.".to_string());
        }

        let wire = mkwire.wire();
        self.init_wire(&wire, plane)
    }

    /// `Perform(theRadius)` (`cxx:326-571`): the self-intersection cut, the
    /// analytical case, the left-handed-plane flip, the "neighbour almost
    /// disappears" retry, then the fillet edge and the two shrunk neighbours.
    pub fn perform(&mut self, radius: f64) -> bool {
        // `cxx:332-335`.
        if self.left.edge.is_none() || self.right.edge.is_none() || radius < CONFUSION {
            return false;
        }

        // `cxx:338-343`.
        let mut def = FilletDef::default();
        let mut is_cut = false;

        // `cxx:347-369`: the analytical algorithm works for non-intersecting
        // arcs only.
        if !self.left.segment || !self.right.segment {
            // `BRepBuilderAPI_MakeWire mkWire(e1, e2);` (`cxx:349`).
            let mut mkwire = MakeWire::new();
            let _ = mkwire.add(self.left.edge.as_ref().unwrap());
            let _ = mkwire.add(self.right.edge.as_ref().unwrap());
            if mkwire.is_done() {
                let wire = mkwire.wire();
                // `BRepBuilderAPI_MakeFace mkFace(plane);` (`cxx:352`): the whole
                // plane, unbounded, and always done.
                let face = TopoBuilder::new().make_face_plane(&self.plane);
                // `ShapeAnalysis_Wire analyzer(W, F, Precision::Confusion());`
                // then `CheckSelfIntersection()` (`cxx:358-359`).
                if check_self_intersection(&wire, &face) {
                    // `cxx:360-366`.
                    is_cut = true;
                    if !self.cut() {
                        return false;
                    }
                }
            }
        }

        // `BRepAdaptor_Curve AC1(e1), AC2(e2);` (`cxx:371`).
        let (Some((c1, f1, l1)), Some((c2, f2, l2))) =
            (self.left.curve_range(), self.right.curve_range())
        else {
            return false;
        };

        // `cxx:372-388`: choose the case.
        let mut b_ret = match (self.left.segment, self.right.segment) {
            (true, true) => {
                fillet_cases::segment_fillet_segment(&self.left, &self.right, radius, &mut def)
            }
            (true, false) => fillet_cases::segment_fillet_arc(&self.left, &self.right, radius, &mut def),
            (false, true) => fillet_cases::arc_fillet_segment(&self.left, &self.right, radius, &mut def),
            (false, false) => fillet_cases::arc_fillet_arc(&self.left, &self.right, radius, &mut def),
        };
        // `cxx:390-393`.
        if !b_ret {
            return false;
        }

        // `cxx:396-399`: invert the fillet for a left-handed plane.
        if !self.plane.is_direct() {
            def.cw = !def.cw;
        }

        // `cxx:402-406`: `gp_Ax2(center, cw ? -normal : normal)`.
        let mut center = slib::plane_value(&self.plane, def.xc, def.yc);
        let normal = self.plane.pos.direction();
        let z_dir = if def.cw { normal.reversed() } else { normal };
        let mut circ = GpCirc::new(GpAx2::from_axis(center, z_dir), radius);

        // `cxx:409-455`: a fillet may only shrink a neighbour edge, it can't
        // prolongate it.
        let delta1 = l1 - f1;
        let delta2 = l2 - f2;
        if !is_cut && (def.start > delta1 || def.end > delta2) {
            // `cxx:412-418`: try to reduce the fillet radius a little
            // (1.e-5 mm) when a neighbour edge almost disappears.
            let little = 100.0 * CONFUSION;
            let d1 = (def.start - delta1).abs();
            let d2 = (def.end - delta2).abs();
            if d1 < little || d2 < little {
                let reduced = radius - little;
                b_ret = match (self.left.segment, self.right.segment) {
                    (true, true) => fillet_cases::segment_fillet_segment(
                        &self.left,
                        &self.right,
                        reduced,
                        &mut def,
                    ),
                    (true, false) => fillet_cases::segment_fillet_arc(
                        &self.left,
                        &self.right,
                        reduced,
                        &mut def,
                    ),
                    (false, true) => fillet_cases::arc_fillet_segment(
                        &self.left,
                        &self.right,
                        reduced,
                        &mut def,
                    ),
                    (false, false) => {
                        fillet_cases::arc_fillet_arc(&self.left, &self.right, reduced, &mut def)
                    }
                };
                if b_ret {
                    // `cxx:437-442`.
                    if !self.plane.is_direct() {
                        def.cw = !def.cw;
                    }
                    // `cxx:444-447`.
                    center = slib::plane_value(&self.plane, def.xc, def.yc);
                    circ.set_location(&center);
                    let _ = circ.set_radius(reduced);
                } else {
                    // `cxx:448-451`.
                    return false;
                }
            } else {
                // `cxx:453-455`.
                return false;
            }
        }

        if b_ret {
            let e1 = self.left.edge.clone().unwrap();
            let e2 = self.right.edge.clone().unwrap();

            // `cxx:461-477`: the fillet start.
            let pstart = if def.xstart != f64::MAX {
                slib::plane_value(&self.plane, def.xstart, def.ystart)
            } else if e1.0.orientation() == Orientation::Forward {
                c1.d0(l1 - def.start)
            } else {
                c1.d0(f1 + def.start)
            };
            // `cxx:478-493`: the fillet end.
            let pend = if def.xend != f64::MAX {
                slib::plane_value(&self.plane, def.xend, def.yend)
            } else if e2.0.orientation() == Orientation::Forward {
                c2.d0(f2 + def.end)
            } else {
                c2.d0(l2 - def.end)
            };

            // `cxx:496-498`: `BRepBuilderAPI_MakeEdge mkEdge(circ, pstart, pend);`
            let fillet_curve: Arc<dyn Curve> = Arc::new(GeomCircle::new(circ));
            match brep_lib_make_edge::init_curve_points(&fillet_curve, &pstart, &pend) {
                Some(build) => {
                    let builder = TopoBuilder::new();
                    let mut fillet = builder.make_edge(fillet_curve.clone(), build.first, build.last);
                    brep_lib_make_edge::add_endpoint_vertices(&builder, &mut fillet, &build);
                    self.fillet = Some(fillet);

                    // `cxx:503-534`: limit the left neighbour.
                    self.shrinke1 = None;
                    let (p1, p2) = if e1.0.orientation() == Orientation::Forward {
                        (c1.d0(f1), pstart)
                    } else {
                        (pstart, c1.d0(l1))
                    };
                    self.shrinke1 = fillet_support::make_edge_points(&c1, &p1, &p2);

                    // `cxx:536-562`: limit the right neighbour. OCCT tests `e1`'s
                    // orientation here, not `e2`'s (`cxx:538`).
                    self.shrinke2 = None;
                    let (p1, p2) = if e1.0.orientation() == Orientation::Forward {
                        (pend, c2.d0(l2))
                    } else {
                        (c2.d0(f2), pend)
                    };
                    self.shrinke2 = fillet_support::make_edge_points(&c2, &p1, &p2);

                    // `cxx:566`.
                    b_ret = self.shrinke1.is_some() && self.shrinke2.is_some();
                }
                // `cxx:498-499`: `bRet = mkEdge.IsDone();`.
                None => b_ret = false,
            }
        }

        // `cxx:570`.
        b_ret
    }

    /// `Result(theE1, theE2)` (`cxx:573-579`) with the three null handles as
    /// `None`: OCCT fills all three or none, so one `Option` carries them.
    pub fn result(&self) -> Option<AnaFilletResult> {
        Some(AnaFilletResult {
            fillet: self.fillet.clone()?,
            shrinke1: self.shrinke1.clone()?,
            shrinke2: self.shrinke2.clone()?,
        })
    }

    /// `if (!Cut(plane, e1, e2)) return false;` (`cxx:362-365`): the two member
    /// edges are replaced and `(x12, y12)` / `(x21, y21)` move to the
    /// intersection point (`cxx:1094-1097`).
    fn cut(&mut self) -> bool {
        let mut e1 = self.left.edge.clone().unwrap();
        let mut e2 = self.right.edge.clone().unwrap();
        match fillet_support::cut(&self.plane, &mut e1, &mut e2) {
            Some(p) => {
                self.left.edge = Some(e1);
                self.right.edge = Some(e2);
                self.left.x2 = p.x();
                self.left.y2 = p.y();
                self.right.x1 = p.x();
                self.right.y1 = p.y();
                true
            }
            None => false,
        }
    }
}

/// `TopoDS_Edge(E.Reversed())` (`cxx:306`, `:313`): a copy sharing the `TShape`
/// with the opposite orientation.
fn reversed(edge: &Edge) -> Edge {
    let mut e = edge.clone();
    e.0.reverse();
    e
}
