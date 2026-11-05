//! `ChFi2d_ChamferAPI` — a chamfer between two linear edges of a plane.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi2d/ChFi2d_ChamferAPI.hxx:25-57` /
//! `.cxx:21-160` (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `ChFi2d_ChamferAPI()` | `cxx:21-30` |
//! | `ChFi2d_ChamferAPI(theWire)` | `cxx:33-43` |
//! | `ChFi2d_ChamferAPI(theEdge1, theEdge2)` | `cxx:46-56` |
//! | `Init(theWire)` | `cxx:60-79` |
//! | `Init(theEdge1, theEdge2)` | `cxx:83-87` |
//! | `Perform()` | `cxx:92-126` |
//! | `Result(E1, E2, L1, L2)` | `cxx:130-159` |
//!
//! Translation notes:
//! - `BRep_Tool::Curve(E, first, last)` is
//!   [`crate::brep_tool::BRepTool::edge_curve`] plus
//!   [`crate::brep_tool::BRepTool::edge_parameters`]; a null handle OCCT would
//!   dereference (`myCurve1->Value(...)`, `cxx:96`) is the `Err` here,
//! - `BRepBuilderAPI_MakeEdge(C, U1, U2)` is `BRepLib_MakeEdge::Init` with null
//!   vertices (`BRepLib_MakeEdge.cxx:512-517`), i.e.
//!   [`crate::brep_lib_make_edge::init_3d`] plus
//!   [`crate::brep_lib_make_edge::add_endpoint_vertices`]; its
//!   `!IsDone()` (`ParameterOutOfRange` `cxx:653-657`, `LineThroughIdenticPoints`
//!   `cxx:660-663`) makes `BRepBuilderAPI_MakeEdge::Edge()` throw, so the port
//!   answers `Err`,
//! - `GC_MakeLine(P1, P2)` + `BRepBuilderAPI_MakeEdge(line, P1, P2)` (`cxx:145-148`)
//!   is the `gp_Lin(P1, dir(P1, P2))` edge of `BRepLib_MakeEdge::Init(C, P1, P2)`
//!   (`cxx:522-541`): the vertices are built on the two points with
//!   `BRepLib::Precision()` (`= Precision::Confusion()`), `V2` aliases `V1` when
//!   the two points are closer than that, and the range is `[0, |P2 - P1|]`,
//!   because `gce_MakeLin(P1, P2)` puts the origin on `P1`,
//! - the null edge `Result` answers when a length does not fit on its edge
//!   (`cxx:136-143`) is `Ok(None)`; `theEdge1` / `theEdge2` keep their values
//!   there, exactly like the two C++ early returns,
//! - `TopoDS::Edge(itr.Value())` (`cxx:67`, `:71`) on a wire child that is not
//!   an edge raises `Standard_TypeMismatch` in OCCT; that arm is the `Err` of
//!   [`ChamferApi::init_wire`].

use std::sync::Arc;

use occt_core::gp::{GpDir, GpLin, GpPnt, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::curve::Curve;
use occt_geom::line::GeomLine;

use crate::brep_lib_make_edge;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::iterator::ShapeIterator;
use crate::shape::{Edge, Wire};

/// `ChFi2d_ChamferAPI` (`hxx:25-57`).
#[derive(Clone, Default)]
pub struct ChamferApi {
    /// `TopoDS_Edge myEdge1, myEdge2;` (`hxx:52`). OCCT leaves them null until
    /// `Init`; `None` is that null handle.
    edge1: Option<Edge>,
    edge2: Option<Edge>,
    /// `occ::handle<Geom_Curve> myCurve1, myCurve2;` (`hxx:53`), filled by
    /// `Perform` and null before it.
    curve1: Option<Arc<dyn Curve>>,
    curve2: Option<Arc<dyn Curve>>,
    /// `double myStart1, myEnd1, myStart2, myEnd2;` (`hxx:54`).
    start1: f64,
    end1: f64,
    start2: f64,
    end2: f64,
    /// `bool myCommonStart1, myCommonStart2;` (`hxx:55`).
    common_start1: bool,
    common_start2: bool,
}

impl ChamferApi {
    /// `ChFi2d_ChamferAPI()` (`cxx:21-30`): every length at `0.0`, both flags
    /// `false`, both edges null.
    pub fn new() -> Self {
        Self::default()
    }

    /// `ChFi2d_ChamferAPI(theWire)` (`cxx:33-43`): the two flags / lengths of
    /// the empty constructor, then `Init(theWire)`.
    pub fn from_wire(wire: &Wire) -> Result<Self, String> {
        let mut api = Self::new();
        api.init_wire(wire)?;
        Ok(api)
    }

    /// `ChFi2d_ChamferAPI(theEdge1, theEdge2)` (`cxx:46-56`).
    pub fn from_edges(edge1: &Edge, edge2: &Edge) -> Self {
        Self {
            edge1: Some(edge1.clone()),
            edge2: Some(edge2.clone()),
            ..Self::new()
        }
    }

    /// `Init(theWire)` (`cxx:60-79`): the first two wire children, in wire
    /// order, then `Init(E1, E2)`.
    pub fn init_wire(&mut self, wire: &Wire) -> Result<(), String> {
        // `TopoDS_Iterator itr(theWire); for (; itr.More(); itr.Next())`
        // (`cxx:62-77`).
        let mut e1: Option<Edge> = None;
        let mut e2: Option<Edge> = None;
        for child in ShapeIterator::of_shape(&wire.0) {
            let edge = Edge::wrap(child)
                .ok_or_else(|| "ChFi2d_ChamferAPI::Init: wire child is not an edge".to_string())?;
            if e1.is_none() {
                e1 = Some(edge);
            } else if e2.is_none() {
                e2 = Some(edge);
            } else {
                break;
            }
        }
        // `Init(E1, E2)` (`cxx:78`). A wire with fewer than two children leaves
        // one of the two handles null, exactly like OCCT.
        self.edge1 = e1;
        self.edge2 = e2;
        Ok(())
    }

    /// `Init(theEdge1, theEdge2)` (`cxx:83-87`).
    pub fn init_edges(&mut self, edge1: &Edge, edge2: &Edge) {
        self.edge1 = Some(edge1.clone());
        self.edge2 = Some(edge2.clone());
    }

    /// `Perform()` (`cxx:92-126`): the two curves and their ranges, then the
    /// four endpoint equality tests that set `myCommonStart1` /
    /// `myCommonStart2`.
    pub fn perform(&mut self) -> Result<bool, String> {
        // `myCurve1 = BRep_Tool::Curve(myEdge1, myStart1, myEnd1);` (`cxx:94`).
        let edge1 = self
            .edge1
            .as_ref()
            .ok_or_else(|| "ChFi2d_ChamferAPI::Perform: myEdge1 is null".to_string())?;
        let edge2 = self
            .edge2
            .as_ref()
            .ok_or_else(|| "ChFi2d_ChamferAPI::Perform: myEdge2 is null".to_string())?;
        self.curve1 = Some(
            BRepTool::edge_curve(edge1)
                .ok_or_else(|| "ChFi2d_ChamferAPI::Perform: myEdge1 has no curve".to_string())?,
        );
        let (s1, e1) = BRepTool::edge_parameters(edge1);
        self.start1 = s1;
        self.end1 = e1;
        self.curve2 = Some(
            BRepTool::edge_curve(edge2)
                .ok_or_else(|| "ChFi2d_ChamferAPI::Perform: myEdge2 has no curve".to_string())?,
        );
        let (s2, e2) = BRepTool::edge_parameters(edge2);
        self.start2 = s2;
        self.end2 = e2;

        let c1 = self.curve1.as_ref().unwrap();
        let c2 = self.curve2.as_ref().unwrap();
        // `if (myCurve1->Value(myStart1).IsEqual(myCurve2->Value(myEnd2),
        //  Precision::Confusion()))` (`cxx:96-100`).
        let p1_start = c1.d0(self.start1);
        let p1_end = c1.d0(self.end1);
        let p2_start = c2.d0(self.start2);
        let p2_end = c2.d0(self.end2);
        if p1_start.is_equal_with_tolerance(&p2_end, CONFUSION) {
            self.common_start1 = true;
            self.common_start2 = false;
        } else if p1_end.is_equal_with_tolerance(&p2_start, CONFUSION) {
            // `cxx:103-107`.
            self.common_start1 = false;
            self.common_start2 = true;
        } else if p1_end.is_equal_with_tolerance(&p2_end, CONFUSION) {
            // `cxx:110-114`.
            self.common_start1 = false;
            self.common_start2 = false;
        } else {
            // `cxx:117-121`.
            self.common_start1 = true;
            self.common_start2 = true;
        }
        // `return true;` (`cxx:125`).
        Ok(true)
    }

    /// `Result(theEdge1, theEdge2, theLength1, theLength2)` (`cxx:130-159`):
    /// the chamfer edge plus the two divided input edges, written back through
    /// `the_edge1` / `the_edge2`.
    ///
    /// `Ok(None)` is the null `aResult` of the two length tests
    /// (`cxx:136-143`), the only arm that leaves `the_edge1` / `the_edge2`
    /// untouched.
    pub fn result(
        &self,
        the_edge1: &mut Edge,
        the_edge2: &mut Edge,
        the_length1: f64,
        the_length2: f64,
    ) -> Result<Option<Edge>, String> {
        // `if (std::abs(myEnd1 - myStart1) < theLength1) return aResult;`
        // (`cxx:136-139`).
        if (self.end1 - self.start1).abs() < the_length1 {
            return Ok(None);
        }
        // `if (std::abs(myEnd2 - myStart2) < theLength2) return aResult;`
        // (`cxx:140-143`).
        if (self.end2 - self.start2).abs() < the_length2 {
            return Ok(None);
        }

        // `cxx:145-149`.
        let a_common1 = (if self.common_start1 { self.start1 } else { self.end1 })
            + if (self.start1 > self.end1) ^ self.common_start1 {
                the_length1
            } else {
                -the_length1
            };
        let a_common2 = (if self.common_start2 { self.start2 } else { self.end2 })
            + if (self.start2 > self.end2) ^ self.common_start2 {
                the_length2
            } else {
                -the_length2
            };

        let curve1 = self
            .curve1
            .as_ref()
            .ok_or_else(|| "ChFi2d_ChamferAPI::Result: myCurve1 is null".to_string())?;
        let curve2 = self
            .curve2
            .as_ref()
            .ok_or_else(|| "ChFi2d_ChamferAPI::Result: myCurve2 is null".to_string())?;
        let p1 = curve1.d0(a_common1);
        let p2 = curve2.d0(a_common2);
        let builder = TopoBuilder::new();

        // `GC_MakeLine aML(myCurve1->Value(aCommon1), myCurve2->Value(aCommon2));`
        // then `BRepBuilderAPI_MakeEdge(aML.Value(), p1, p2)` (`cxx:145-148`).
        let a_result = make_line_edge(&builder, &p1, &p2)?;

        // `BRepBuilderAPI_MakeEdge aDivider1(myCurve1, aCommon1,
        //  (myCommonStart1 ? myEnd1 : myStart1));` (`cxx:150-151`).
        *the_edge1 = make_edge_over(
            &builder,
            curve1,
            a_common1,
            if self.common_start1 { self.end1 } else { self.start1 },
        )?;
        // `BRepBuilderAPI_MakeEdge aDivider2(myCurve2, aCommon2,
        //  (myCommonStart2 ? myEnd2 : myStart2));` (`cxx:153-154`).
        *the_edge2 = make_edge_over(
            &builder,
            curve2,
            a_common2,
            if self.common_start2 { self.end2 } else { self.start2 },
        )?;

        // `return aResult;` (`cxx:158`).
        Ok(Some(a_result))
    }
}

/// `BRepLib_MakeEdge::Init(C, p1, p2)` (`cxx:512-517`) — the
/// `BRepBuilderAPI_MakeEdge(C, U1, U2)` form: null vertices, so `Init` builds
/// them on the curve's two endpoints (`BRepLib_MakeEdge.cxx:740-775`) and
/// reorders the range.
fn make_edge_over(
    builder: &TopoBuilder,
    curve: &Arc<dyn Curve>,
    first: f64,
    last: f64,
) -> Result<Edge, String> {
    let build = brep_lib_make_edge::init_3d(curve, first, last).ok_or_else(|| {
        "ChFi2d_ChamferAPI::Result: BRepLib_MakeEdge not done (range / identity)".to_string()
    })?;
    let mut e = builder.make_edge(curve.clone(), build.first, build.last);
    brep_lib_make_edge::add_endpoint_vertices(builder, &mut e, &build);
    Ok(e)
}

/// `GC_MakeLine(P1, P2)` + `BRepBuilderAPI_MakeEdge(line, P1, P2)`
/// (`ChFi2d_ChamferAPI.cxx:145-148`).
///
/// `gce_MakeLin(P1, P2)` (`gce_MakeLin.cxx:26-40`) fails on confused points and
/// then `aML.Value()` is null, which `BRepBuilderAPI_MakeEdge` dereferences —
/// the `Err` here. `BRepLib_MakeEdge::Init(C, P1, P2)`
/// (`BRepLib_MakeEdge.cxx:522-541`) is `gp_Lin(P1, dir(P2 - P1))` with the two
/// endpoint vertices at `BRepLib::Precision()` and `V2 = V1` when the points
/// are closer than that.
fn make_line_edge(builder: &TopoBuilder, p1: &GpPnt, p2: &GpPnt) -> Result<Edge, String> {
    let dir = GpDir::from_vec(&GpVec::from_pnts(p1, p2))
        .map_err(|_| "ChFi2d_ChamferAPI::Result: GC_MakeLine on confused points".to_string())?;
    let dist = p1.distance(p2);
    let lin = GpLin::from_pnt_dir(*p1, dir);
    let mut e = builder.make_edge(Arc::new(GeomLine::new(lin)), 0.0, dist);
    let v1 = builder.make_vertex(*p1, CONFUSION);
    let v2 = if dist < CONFUSION {
        // `BRepLib_MakeEdge.cxx:530-537`: `if (P1.Distance(P2) < Tol) V2 = V1;`.
        e.set_closed(true);
        v1.clone()
    } else {
        builder.make_vertex(*p2, CONFUSION)
    };
    builder.add_edge_vertices(&mut e, &v1, &v2);
    Ok(e)
}
