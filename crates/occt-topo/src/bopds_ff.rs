//! F/F interference payload of the DS: `BOPDS_Curve`, `BOPDS_Point`, `BOPDS_InterfFF`.
//!
//! Source: `BOPDS_Curve.hxx/.lxx`, `BOPDS_Point.hxx/.lxx`, `BOPDS_Interf.hxx`
//! (`BOPDS_InterfFF`). `BOPAlgo_PaveFiller::PerformFF` (`BOPAlgo_PaveFiller_6.cxx`)
//! writes these records; `MakeBlocks` / `PostTreatFF` in the same file read them
//! to build section edges. The 3-D curve plus pcurves are the payload of
//! `IntTools_Curve`, held here as `Arc<dyn Curve>` / `Curve2d`.

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::{GpPnt, GpPnt2d};
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;

use crate::bopds::BopdsPaveBlock;

/// Intersection point of two faces. Source: `BOPDS_Point`.
#[derive(Debug, Clone)]
pub struct BopdsPoint {
    pnt: GpPnt,
    pnt2d1: GpPnt2d,
    pnt2d2: GpPnt2d,
    index: i64,
}

impl BopdsPoint {
    /// Empty constructor (`myPnt` / `myPnt2D*` at 99, `myIndex` -1).
    pub fn new() -> Self {
        Self {
            pnt: GpPnt::new(99.0, 99.0, 99.0),
            pnt2d1: GpPnt2d::new(99.0, 99.0),
            pnt2d2: GpPnt2d::new(99.0, 99.0),
            index: -1,
        }
    }

    pub fn set_pnt(&mut self, p: GpPnt) {
        self.pnt = p;
    }

    pub fn pnt(&self) -> &GpPnt {
        &self.pnt
    }

    pub fn set_pnt2d1(&mut self, p: GpPnt2d) {
        self.pnt2d1 = p;
    }

    pub fn pnt2d1(&self) -> &GpPnt2d {
        &self.pnt2d1
    }

    pub fn set_pnt2d2(&mut self, p: GpPnt2d) {
        self.pnt2d2 = p;
    }

    pub fn pnt2d2(&self) -> &GpPnt2d {
        &self.pnt2d2
    }

    pub fn set_index(&mut self, index: usize) {
        self.index = index as i64;
    }

    pub fn index(&self) -> Option<usize> {
        if self.index >= 0 {
            Some(self.index as usize)
        } else {
            None
        }
    }
}

impl Default for BopdsPoint {
    fn default() -> Self {
        Self::new()
    }
}

/// Intersection curve of two faces. Source: `BOPDS_Curve`.
#[derive(Clone)]
pub struct BopdsCurve {
    curve: Option<Arc<dyn Curve>>,
    pcurve1: Option<Arc<dyn Curve2d>>,
    pcurve2: Option<Arc<dyn Curve2d>>,
    first: f64,
    last: f64,
    pave_blocks: Vec<BopdsPaveBlock>,
    techno_vertices: Vec<usize>,
    bbox: BndBox,
    tolerance: f64,
    tangential_tolerance: f64,
}

impl BopdsCurve {
    /// Empty constructor.
    pub fn new() -> Self {
        Self {
            curve: None,
            pcurve1: None,
            pcurve2: None,
            first: 0.0,
            last: 0.0,
            pave_blocks: Vec::new(),
            techno_vertices: Vec::new(),
            bbox: BndBox::new(),
            tolerance: 0.0,
            tangential_tolerance: 0.0,
        }
    }

    pub fn set_curve(&mut self, c: Arc<dyn Curve>) {
        self.curve = Some(c);
    }

    pub fn curve(&self) -> Option<&Arc<dyn Curve>> {
        self.curve.as_ref()
    }

    pub fn set_pcurves(&mut self, c1: Option<Arc<dyn Curve2d>>, c2: Option<Arc<dyn Curve2d>>) {
        self.pcurve1 = c1;
        self.pcurve2 = c2;
    }

    pub fn pcurve1(&self) -> Option<&Arc<dyn Curve2d>> {
        self.pcurve1.as_ref()
    }

    pub fn pcurve2(&self) -> Option<&Arc<dyn Curve2d>> {
        self.pcurve2.as_ref()
    }

    pub fn set_range(&mut self, first: f64, last: f64) {
        self.first = first;
        self.last = last;
    }

    pub fn range(&self) -> (f64, f64) {
        (self.first, self.last)
    }

    pub fn set_box(&mut self, bbox: BndBox) {
        self.bbox = bbox;
    }

    pub fn bounding_box(&self) -> &BndBox {
        &self.bbox
    }

    pub fn change_box(&mut self) -> &mut BndBox {
        &mut self.bbox
    }

    pub fn set_pave_blocks(&mut self, pbs: Vec<BopdsPaveBlock>) {
        self.pave_blocks = pbs;
    }

    pub fn pave_blocks(&self) -> &[BopdsPaveBlock] {
        &self.pave_blocks
    }

    pub fn change_pave_blocks(&mut self) -> &mut Vec<BopdsPaveBlock> {
        &mut self.pave_blocks
    }

    /// Creates the initial pave block of the curve (`BOPDS_Curve::InitPaveBlock1`).
    pub fn init_pave_block1(&mut self) {
        if self.pave_blocks.is_empty() {
            self.pave_blocks.push(BopdsPaveBlock::new());
        }
    }

    pub fn change_pave_block1(&mut self) -> &mut BopdsPaveBlock {
        self.init_pave_block1();
        &mut self.pave_blocks[0]
    }

    pub fn techno_vertices(&self) -> &[usize] {
        &self.techno_vertices
    }

    pub fn change_techno_vertices(&mut self) -> &mut Vec<usize> {
        &mut self.techno_vertices
    }

    pub fn has_edge(&self) -> bool {
        self.pave_blocks.iter().any(|pb| pb.has_edge() && pb.edge() != 0)
    }

    pub fn set_tolerance(&mut self, tol: f64) {
        self.tolerance = tol;
    }

    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// `BOPDS_Curve::SetTangentialTolerance` is not a setter on the DS curve;
    /// the value lives on `IntTools_Curve`. Stored here so `CorrectToleranceOfSE`
    /// can read `TangentialTolerance()`.
    pub fn set_tangential_tolerance(&mut self, tol: f64) {
        self.tangential_tolerance = tol;
    }

    /// `BOPDS_Curve::TangentialTolerance` (`myCurve.TangentialTolerance()`).
    pub fn tangential_tolerance(&self) -> f64 {
        self.tangential_tolerance
    }
}

impl Default for BopdsCurve {
    fn default() -> Self {
        Self::new()
    }
}

impl std::fmt::Debug for BopdsCurve {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BopdsCurve")
            .field("has_curve", &self.curve.is_some())
            .field("range", &(self.first, self.last))
            .field("pave_blocks", &self.pave_blocks.len())
            .field("tolerance", &self.tolerance)
            .field("tangential_tolerance", &self.tangential_tolerance)
            .finish()
    }
}

/// Face/face interference. Source: `BOPDS_InterfFF`.
#[derive(Debug, Clone)]
pub struct BopdsInterfFf {
    index1: usize,
    index2: usize,
    tangent_faces: bool,
    curves: Vec<BopdsCurve>,
    points: Vec<BopdsPoint>,
}

impl BopdsInterfFf {
    /// Empty constructor with the two face indices.
    pub fn new(index1: usize, index2: usize) -> Self {
        Self {
            index1,
            index2,
            tangent_faces: false,
            curves: Vec::new(),
            points: Vec::new(),
        }
    }

    pub fn set_indices(&mut self, i1: usize, i2: usize) {
        self.index1 = i1;
        self.index2 = i2;
    }

    pub fn indices(&self) -> (usize, usize) {
        (self.index1, self.index2)
    }

    pub fn index1(&self) -> usize {
        self.index1
    }

    pub fn index2(&self) -> usize {
        self.index2
    }

    /// Reserve curve/point capacity (`BOPDS_InterfFF::Init`).
    pub fn init(&mut self, nb_curves: usize, nb_points: usize) {
        if nb_curves > 0 {
            self.curves.reserve(nb_curves);
        }
        if nb_points > 0 {
            self.points.reserve(nb_points);
        }
    }

    pub fn set_tangent_faces(&mut self, flag: bool) {
        self.tangent_faces = flag;
    }

    pub fn tangent_faces(&self) -> bool {
        self.tangent_faces
    }

    pub fn curves(&self) -> &[BopdsCurve] {
        &self.curves
    }

    pub fn change_curves(&mut self) -> &mut Vec<BopdsCurve> {
        &mut self.curves
    }

    pub fn points(&self) -> &[BopdsPoint] {
        &self.points
    }

    pub fn change_points(&mut self) -> &mut Vec<BopdsPoint> {
        &mut self.points
    }
}
