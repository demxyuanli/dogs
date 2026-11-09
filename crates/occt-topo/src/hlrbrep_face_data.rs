//! `HLRBRep_FaceData` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_FaceData.hxx:35-186`,
//! `HLRBRep_FaceData.cxx:21-90`, `HLRBRep_FaceData.lxx:19-303`): one face of the
//! exact HLR data structure - the surface, the wire/edge blocks and the
//! visibility flags the hidden-line removal toggles.
//!
//! The OCCT `Wires()` accessor returns the `occ::handle<HLRAlgo_WiresBlock>&`
//! so callers mutate the shared block through the handle. The port stores
//! `Arc<HlrAlgoWiresBlock>` and hands out `&mut` through `Arc::make_mut`
//! (the block is `Clone`), which keeps the same "one face owns its wire block"
//! semantics without interior mutability.

use std::sync::Arc;

use crate::abs::Orientation;
use crate::hlr_algo::{HlrAlgoEdgesBlock, HlrAlgoWiresBlock};
use crate::hlrbrep_surface::HlrBRepSurface;
use crate::shape::Face;
use crate::topopebrep::tool::shape_tool;

/// `HLRBRep_FaceData::EMaskFlags` (`hxx:152-172`).
const EMASK_ORIENT: i32 = 15;
const FMASK_SELECTED: i32 = 16;
const FMASK_BACK: i32 = 32;
const FMASK_SIDE: i32 = 64;
const FMASK_CLOSED: i32 = 128;
const FMASK_HIDING: i32 = 256;
const FMASK_SIMPLE: i32 = 512;
const FMASK_CUT: i32 = 1024;
const FMASK_WITH_OUT_L: i32 = 2048;
const FMASK_PLANE: i32 = 4096;
const FMASK_CYLINDER: i32 = 8192;
const FMASK_CONE: i32 = 16384;
const FMASK_SPHERE: i32 = 32768;
const FMASK_TORUS: i32 = 65536;

/// `HLRBRep_FaceData`.
#[derive(Clone)]
pub struct HlrBRepFaceData {
    /// `int myFlags` (`hxx:176`).
    my_flags: i32,
    /// `occ::handle<HLRAlgo_WiresBlock> myWires` (`hxx:177`).
    my_wires: Option<Arc<HlrAlgoWiresBlock>>,
    /// `HLRBRep_Surface myGeometry` (`hxx:178`).
    my_geometry: HlrBRepSurface,
    /// `double mySize` (`hxx:179`).
    my_size: f64,
    /// `float myTolerance` (`hxx:180`).
    my_tolerance: f32,
}

impl Default for HlrBRepFaceData {
    /// `HLRBRep_FaceData()` (`cxx:32-37`): `myFlags = 0`, `mySize = 0`, then
    /// `Selected(true)`.
    fn default() -> Self {
        let mut self_ = Self {
            my_flags: 0,
            my_wires: None,
            my_geometry: HlrBRepSurface::new(),
            my_size: 0.0,
            my_tolerance: 0.0,
        };
        self_.set_selected(true);
        self_
    }
}

impl HlrBRepFaceData {
    /// `HLRBRep_FaceData()` (`cxx:32-37`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Set(FG, Or, Cl, NW)` (`cxx:41-51`).
    pub fn set(&mut self, face: &Face, or: Orientation, closed: bool, nb_wires: usize) {
        self.set_closed(closed);
        self.my_geometry.set_surface(face);
        self.my_tolerance = shape_tool::tolerance_of(face.shape())
            .map(|t| t as f32)
            .unwrap_or(0.0);
        self.set_orientation(or);
        self.my_wires = Some(Arc::new(HlrAlgoWiresBlock::new(nb_wires)));
    }

    /// `SetWire(WI, NE)` (`cxx:55-58`).
    pub fn set_wire(&mut self, wi: usize, nb_edges: usize) {
        if let Some(wires) = self.my_wires.as_mut() {
            Arc::make_mut(wires).set(wi, Arc::new(HlrAlgoEdgesBlock::new(nb_edges)));
        }
    }

    /// `SetWEdge(WI, EWI, EI, Or, OutL, Inte, Dble, IsoL)` (`cxx:62-90`).
    #[allow(clippy::too_many_arguments)]
    pub fn set_wedge(
        &mut self,
        wi: usize,
        ewi: usize,
        ei: i32,
        or: Orientation,
        out_l: bool,
        inte: bool,
        dble: bool,
        iso_l: bool,
    ) {
        let Some(wires) = self.my_wires.as_mut() else {
            return;
        };
        let Some(wire) = Arc::make_mut(wires).wire_mut(wi) else {
            return;
        };
        wire.set_edge(ewi, ei);
        wire.set_orientation(ewi, or);
        wire.set_out_line(ewi, out_l);
        wire.set_internal(ewi, inte);
        wire.set_double(ewi, dble);
        wire.set_iso_line(ewi, iso_l);
    }

    fn mask(&self, m: i32) -> bool {
        (self.my_flags & m) != 0
    }

    fn set_mask(&mut self, m: i32, b: bool) {
        if b {
            self.my_flags |= m;
        } else {
            self.my_flags &= !m;
        }
    }

    /// `Selected()` / `Selected(B)` (`lxx:19-31`).
    pub fn selected(&self) -> bool {
        self.mask(FMASK_SELECTED)
    }

    pub fn set_selected(&mut self, b: bool) {
        self.set_mask(FMASK_SELECTED, b);
    }

    /// `Back()` / `Back(B)` (`lxx:35-47`).
    pub fn back(&self) -> bool {
        self.mask(FMASK_BACK)
    }

    pub fn set_back(&mut self, b: bool) {
        self.set_mask(FMASK_BACK, b);
    }

    /// `Side()` / `Side(B)` (`lxx:51-63`).
    pub fn side(&self) -> bool {
        self.mask(FMASK_SIDE)
    }

    pub fn set_side(&mut self, b: bool) {
        self.set_mask(FMASK_SIDE, b);
    }

    /// `Closed()` / `Closed(B)` (`lxx:67-79`).
    pub fn closed(&self) -> bool {
        self.mask(FMASK_CLOSED)
    }

    pub fn set_closed(&mut self, b: bool) {
        self.set_mask(FMASK_CLOSED, b);
    }

    /// `Hiding()` / `Hiding(B)` (`lxx:83-95`).
    pub fn hiding(&self) -> bool {
        self.mask(FMASK_HIDING)
    }

    pub fn set_hiding(&mut self, b: bool) {
        self.set_mask(FMASK_HIDING, b);
    }

    /// `Simple()` / `Simple(B)` (`lxx:99-111`).
    pub fn simple(&self) -> bool {
        self.mask(FMASK_SIMPLE)
    }

    pub fn set_simple(&mut self, b: bool) {
        self.set_mask(FMASK_SIMPLE, b);
    }

    /// `Cut()` / `Cut(B)` (`lxx:115-127`).
    pub fn cut(&self) -> bool {
        self.mask(FMASK_CUT)
    }

    pub fn set_cut(&mut self, b: bool) {
        self.set_mask(FMASK_CUT, b);
    }

    /// `WithOutL()` / `WithOutL(B)` (`lxx:131-143`).
    pub fn with_out_l(&self) -> bool {
        self.mask(FMASK_WITH_OUT_L)
    }

    pub fn set_with_out_l(&mut self, b: bool) {
        self.set_mask(FMASK_WITH_OUT_L, b);
    }

    /// `Plane()` / `Plane(B)` (`lxx:147-159`).
    pub fn plane(&self) -> bool {
        self.mask(FMASK_PLANE)
    }

    pub fn set_plane(&mut self, b: bool) {
        self.set_mask(FMASK_PLANE, b);
    }

    /// `Cylinder()` / `Cylinder(B)` (`lxx:163-175`).
    pub fn cylinder(&self) -> bool {
        self.mask(FMASK_CYLINDER)
    }

    pub fn set_cylinder(&mut self, b: bool) {
        self.set_mask(FMASK_CYLINDER, b);
    }

    /// `Cone()` / `Cone(B)` (`lxx:179-191`).
    pub fn cone(&self) -> bool {
        self.mask(FMASK_CONE)
    }

    pub fn set_cone(&mut self, b: bool) {
        self.set_mask(FMASK_CONE, b);
    }

    /// `Sphere()` / `Sphere(B)` (`lxx:195-207`).
    pub fn sphere(&self) -> bool {
        self.mask(FMASK_SPHERE)
    }

    pub fn set_sphere(&mut self, b: bool) {
        self.set_mask(FMASK_SPHERE, b);
    }

    /// `Torus()` / `Torus(B)` (`lxx:211-223`).
    pub fn torus(&self) -> bool {
        self.mask(FMASK_TORUS)
    }

    pub fn set_torus(&mut self, b: bool) {
        self.set_mask(FMASK_TORUS, b);
    }

    /// `Size()` / `Size(S)` (`lxx:231-241`).
    pub fn size(&self) -> f64 {
        self.my_size
    }

    pub fn set_size(&mut self, s: f64) {
        self.my_size = s;
    }

    /// `Orientation()` (`lxx:245-248`): the low four bits of `myFlags`.
    pub fn orientation(&self) -> Orientation {
        match self.my_flags & EMASK_ORIENT {
            1 => Orientation::Reversed,
            2 => Orientation::Internal,
            3 => Orientation::External,
            _ => Orientation::Forward,
        }
    }

    /// `Orientation(O)` (`lxx:252-256`).
    pub fn set_orientation(&mut self, o: Orientation) {
        let bits = match o {
            Orientation::Forward => 0,
            Orientation::Reversed => 1,
            Orientation::Internal => 2,
            Orientation::External => 3,
        };
        self.my_flags &= !EMASK_ORIENT;
        self.my_flags |= bits & EMASK_ORIENT;
    }

    /// `Wires()` (`lxx:260-263`).
    pub fn wires(&self) -> Option<&Arc<HlrAlgoWiresBlock>> {
        self.my_wires.as_ref()
    }

    /// `Wires()` mutable overload.
    pub fn wires_mut(&mut self) -> Option<&mut HlrAlgoWiresBlock> {
        self.my_wires.as_mut().map(Arc::make_mut)
    }

    /// `Geometry()` (`lxx:267-270`).
    pub fn geometry(&self) -> &HlrBRepSurface {
        &self.my_geometry
    }

    /// `Geometry()` mutable overload.
    pub fn geometry_mut(&mut self) -> &mut HlrBRepSurface {
        &mut self.my_geometry
    }

    /// `Tolerance()` (`lxx:299-302`).
    pub fn tolerance(&self) -> f32 {
        self.my_tolerance
    }
}
