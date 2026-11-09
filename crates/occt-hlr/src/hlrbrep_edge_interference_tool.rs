//! `HLRBRep_EdgeInterferenceTool`
//! (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_EdgeInterferenceTool.{hxx,cxx,lxx}`,
//! `hxx:35-89`, `cxx:27-100`, `lxx:19-75`): the tool object the
//! `HLRAlgo_Interference` list builder drives while it classifies the
//! interferences of one edge.
//!
//! `LoadEdge` (`cxx:34-47`) fills the two `inter[]` slots with the **bounds** of
//! the current edge - `(parameter, tolerance, vertex index)` of the start and
//! of the end - so `SameVertexAndInterference` can recognise an interference
//! sitting exactly on one of them. `EdgeGeometry` /`InterferenceBoundaryGeometry`
//! (`cxx:51-100`) are the two 2D frames the transition computation needs: the
//! tangent/normal/curvature of the edge at a parameter, and of the face edge a
//! boundary interference names. Both go through
//! `HLRBRep_Data::LocalLEGeometry2D` / `LocalFEGeometry2D` and lift the 2D
//! directions into `gp_Dir(x, y, 0)`.
//!
//! `SameInterferences` (`cxx:64-74`) only answers `true` when both interferences
//! carry a non-zero vertex index and the two indices are equal; index `0` means
//! "not a vertex" and never matches.
//!
//! Deviations:
//! * `myDS` (`hxx:86`) is a `handle<HLRBRep_Data>` in OCCT; the port takes
//!   `&HlrBRepData` on the three methods that read it (`LoadEdge`,
//!   `EdgeGeometry`, `InterferenceBoundaryGeometry`), so the tool carries no
//!   self-reference.
//! * `EdgeGeometry` / `InterferenceBoundaryGeometry` write `gp_Dir` by
//!   reference in OCCT; the port returns `Option<(GpDir, GpDir, f64)>`, `None`
//!   where OCCT would have read an undefined tangent (its `Local*Geometry2D`
//!   throws before writing).
//! * `inter[]` is `HLRAlgo_Intersection inter[2]` (`hxx:87`) - the two edge
//!   bounds, not the interferences. `CurrentVertex()` (`lxx:44-47`) returns
//!   `inter[cur]` directly, so `HLRBRep_VertexList::Current` sees an
//!   `HLRAlgo_Intersection` for boundary vertices and the interference
//!   `.Intersection()` for interference vertices.

use occt_core::gp::GpDir;

use occt_topo::abs::Orientation;
use crate::hlr_algo::{HlrAlgoInterference, HlrAlgoIntersection};
use crate::hlrbrep_data::HlrBRepData;

/// `HLRBRep_EdgeInterferenceTool` (`hxx:35-89`).
#[derive(Debug, Clone)]
pub struct HlrBRepEdgeInterferenceTool {
    /// `HLRAlgo_Intersection inter[2]` (`hxx:87`).
    inter: [HlrAlgoIntersection; 2],
    /// `int cur` (`hxx:88`): the vertex cursor, `0` = start, `1` = end.
    cur: usize,
}

impl Default for HlrBRepEdgeInterferenceTool {
    fn default() -> Self {
        Self::new()
    }
}

impl HlrBRepEdgeInterferenceTool {
    /// `HLRBRep_EdgeInterferenceTool(DS)` (`cxx:27-30`): the constructor only
    /// stores the data handle, which the port passes per call instead.
    pub fn new() -> Self {
        Self {
            inter: [HlrAlgoIntersection::new(); 2],
            cur: 0,
        }
    }

    /// `LoadEdge()` (`cxx:34-47`): `(p1, t1, VSta)` into `inter[0]` and
    /// `(p2, t2, VEnd)` into `inter[1]`, from the status bounds of
    /// `myDS->Edge()`.
    pub fn load_edge(&mut self, ds: &HlrBRepData) {
        let ed = ds.e_data(ds.edge());
        let (p1, t1, p2, t2) = ed.status().bounds();
        self.inter[0].set_parameter(p1);
        self.inter[0].set_tolerance(t1);
        self.inter[0].set_index(ed.v_sta());
        self.inter[1].set_parameter(p2);
        self.inter[1].set_tolerance(t2);
        self.inter[1].set_index(ed.v_end());
    }

    /// `InitVertices()` (`lxx:23-26`).
    pub fn init_vertices(&mut self) {
        self.cur = 0;
    }

    /// `MoreVertices()` (`lxx:30-33`).
    pub fn more_vertices(&self) -> bool {
        self.cur < 2
    }

    /// `NextVertex()` (`lxx:37-40`).
    pub fn next_vertex(&mut self) {
        self.cur += 1;
    }

    /// `CurrentVertex()` (`lxx:44-47`): `inter[cur]`.
    pub fn current_vertex(&self) -> &HlrAlgoIntersection {
        &self.inter[self.cur]
    }

    /// `CurrentOrientation()` (`lxx:51-59`): `FORWARD` at the start,
    /// `REVERSED` at the end.
    pub fn current_orientation(&self) -> Orientation {
        if self.cur == 0 {
            Orientation::Forward
        } else {
            Orientation::Reversed
        }
    }

    /// `CurrentParameter()` (`lxx:63-66`).
    pub fn current_parameter(&self) -> f64 {
        self.inter[self.cur].parameter()
    }

    /// `IsPeriodic()` (`lxx:70-73`): always `false`.
    pub fn is_periodic(&self) -> bool {
        false
    }

    /// `ParameterOfInterference(I)` (`lxx:72-75`).
    pub fn parameter_of_interference(&self, i: &HlrAlgoInterference) -> f64 {
        i.intersection().parameter()
    }

    /// `EdgeGeometry(Param, Tgt, Nrm, CrLE)` (`cxx:51-59`): the 2D frame of the
    /// current edge lifted into the `z = 0` plane.
    pub fn edge_geometry(
        &self,
        ds: &HlrBRepData,
        param: f64,
    ) -> Option<(GpDir, GpDir, f64)> {
        let (tg, nm, cr) = ds.local_le_geometry_2d(param)?;
        Some((
            GpDir::new(tg.x, tg.y, 0.0).unwrap_or_else(|_| GpDir::default_dir()),
            GpDir::new(nm.x, nm.y, 0.0).unwrap_or_else(|_| GpDir::default_dir()),
            cr,
        ))
    }

    /// `SameInterferences(I1, I2)` (`cxx:64-74`): equal non-zero vertex indices.
    pub fn same_interferences(&self, i1: &HlrAlgoInterference, i2: &HlrAlgoInterference) -> bool {
        let ind1 = i1.intersection().index();
        let ind2 = i2.intersection().index();
        if ind1 != 0 && ind2 != 0 {
            return ind1 == ind2;
        }
        false
    }

    /// `SameVertexAndInterference(I)` (`cxx:78-85`): the interference sits on
    /// the current vertex, or its orientation matches the current one.
    pub fn same_vertex_and_interference(&self, i: &HlrAlgoInterference) -> bool {
        if i.intersection().index() == self.inter[self.cur].index() {
            return true;
        }
        i.intersection().orientation() == self.current_orientation()
    }

    /// `InterferenceBoundaryGeometry(I, Tang, Norm, CrFE)` (`cxx:90-100`): the
    /// 2D frame of the face edge named by the interference boundary.
    pub fn interference_boundary_geometry(
        &self,
        ds: &HlrBRepData,
        i: &HlrAlgoInterference,
    ) -> Option<(GpDir, GpDir, f64)> {
        let (fe, param) = i.boundary().value_2d();
        let (tg, nm, cr) = ds.local_fe_geometry_2d(fe, param)?;
        Some((
            GpDir::new(tg.x, tg.y, 0.0).unwrap_or_else(|_| GpDir::default_dir()),
            GpDir::new(nm.x, nm.y, 0.0).unwrap_or_else(|_| GpDir::default_dir()),
            cr,
        ))
    }
}
