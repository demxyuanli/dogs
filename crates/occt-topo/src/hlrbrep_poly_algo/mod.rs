//! Hidden-line removal on triangulations. Source: `HLRBRep_PolyAlgo`
//! (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_PolyAlgo.hxx:93-638`,
//! `HLRBRep_PolyAlgo.cxx:73-4452`).
//!
//! Module layout (by OCCT responsibility, not by line count):
//!
//! - this file: framework (`cxx:73-287`, `:4331-4452`), the shape / edge / face
//!   maps and the projection matrix seams;
//! - [`node_adjacency`]: the node and adjacency layer (`cxx:3513-3996`);
//! - [`bi_points`]: the bi-point assembly cluster (`cxx:878-1377` + `:1379-1886`
//!   + `:1888-2978`);
//! - [`outlines`]: the outline maintenance passes (`cxx:2980-3161`,
//!   `:3163-3511`, `:3998-4153`);
//! - [`finalize`]: `UpdateEdgesBiPoints` / `UpdatePolyData` (`cxx:4155-4329`).
//!
//! Mirrored OCCT members that are only written by `StoreShell` (`TMat`/`TLoc`,
//! `TTMa`/`TTLo`, `TIMa`/`TILo`) exist here and stay at their OCCT initial value
//! (identity translation, zero matrix) until that batch lands; the two helpers
//! that fill them are placed at the same points as the OCCT assignments
//! (`cxx:302-341` and `cxx:429-446`).
//!
//! UNPORTED list: `Update` (`cxx:149-201`), `StoreShell` (`cxx:289-623`),
//! `Normal` (`cxx:626-715`).
//!
//! `StoreShell` / `Normal` need no new geometry primitive: `BRep_Tool::Surface`
//! is [`GeometryRegistry::face_surface`] (the local surface, exactly what
//! `BRep_Tool::Surface(theFace, theLoc)` returns), `GeomAdaptor_Surface`'s
//! `D1`/`D2` are [`occt_geom::Surface::d1`]/[`occt_geom::Surface::d2`] (the five
//! elementary surfaces carry the analytic `ElSLib` second derivatives) and
//! `GeomAdaptor_Surface::GetType() == GeomAbs_Plane` is
//! [`crate::brep_surface::is_plane_type`].
//!
//! The bi-point cluster reads its polygons through
//! [`crate::brep_tool::BRepTool::polygon_on_triangulation`]; nothing commits a
//! `BRep_PolygonOnTriangulation` representation yet (`BRepMesh_ModelPostProcessor`'s
//! `PolygonCommitter` is not ported), so the connex arms currently answer exactly
//! what OCCT answers for an edge that was never meshed. See [`bi_points`] for the
//! per-line deviations.

mod bi_points;
mod finalize;
mod node_adjacency;
mod outlines;

use std::collections::{HashMap, HashSet};

use occt_core::gp::{GpTrsf, GpXyz};

use crate::abs::ShapeType;
use crate::hlr_algo::{
    BiPointPoints, HlrAlgoEdgeStatus, HlrAlgoPolyAlgo, HlrAlgoPolyInternalData,
    HlrAlgoPolyInternalSegment,
};
use crate::shape::TopoShape;
use crate::topo_tools::{
    explorer_shapes, explorer_shapes_avoiding, shape_map_key, ShapeMapKey,
};

/// `NCollection_Array1<occ::handle<HLRAlgo_PolyInternalData>>` (`cxx:289`): the
/// per-face internal data, 1-based with slot 0 unused and `None` for the null
/// slots OCCT leaves behind (a face with no triangulation).
pub type PolyInternalDataArray = Vec<Option<HlrAlgoPolyInternalData>>;

/// `HLRBRep_PolyAlgo` private node flags (`cxx:57-66`).
pub mod nmask {
    /// `NMsk_Vert` (`cxx:59`).
    pub const VERT: i32 = 1;
    pub const OUT_L: i32 = 2;
    pub const NORM: i32 = 4;
    pub const FUCK: i32 = 8;
    pub const EDGE: i32 = 16;
    pub const MOVE: i32 = 32;
}

/// `NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>`: 1-based
/// lookup by [`ShapeMapKey`] with insertion order preserved.
#[derive(Debug, Clone, Default)]
pub struct ShapeIndexedMap {
    items: Vec<TopoShape>,
    keys: HashMap<ShapeMapKey, usize>,
}

impl ShapeIndexedMap {
    pub fn new() -> Self {
        Self::default()
    }

    /// `NCollection_IndexedMap::Add` — returns the 1-based index of `s`,
    /// inserting it when absent.
    pub fn add(&mut self, s: &TopoShape) -> usize {
        let key = shape_map_key(s);
        if let Some(&i) = self.keys.get(&key) {
            return i;
        }
        self.items.push(s.clone());
        let i = self.items.len();
        self.keys.insert(key, i);
        i
    }

    /// `NCollection_IndexedMap::FindIndex` — 0 when `s` is absent.
    pub fn find_index(&self, s: &TopoShape) -> usize {
        self.keys.get(&shape_map_key(s)).copied().unwrap_or(0)
    }

    /// `NCollection_IndexedMap::Extent` (`FindKey` range upper bound).
    pub fn extent(&self) -> usize {
        self.items.len()
    }

    /// `NCollection_IndexedMap::IsEmpty`.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// `NCollection_IndexedMap::Clear`.
    pub fn clear(&mut self) {
        self.items.clear();
        self.keys.clear();
    }

    /// `NCollection_IndexedMap::FindKey(I)` (`hxx` accessor), 1-based.
    pub fn find_key(&self, i: usize) -> Option<&TopoShape> {
        if i >= 1 {
            self.items.get(i - 1)
        } else {
            None
        }
    }
}

/// The four `HLRBRep_PolyAlgo::Hide` out-parameters together with the `PointsT`
/// it returns by reference (`cxx:4370-4383`).
#[derive(Debug, Clone)]
pub struct PolyAlgoSegment {
    pub points: BiPointPoints,
    pub shape: TopoShape,
    pub reg1: bool,
    pub regn: bool,
    pub outl: bool,
    pub intl: bool,
}

/// Calculates the visible and hidden parts of a triangulated shape in a
/// projection. Source: `HLRBRep_PolyAlgo`.
#[derive(Debug, Clone)]
pub struct HlrBRepPolyAlgo {
    proj: crate::hlr_algo::HlrAlgoProjector,
    t_mat: [[f64; 3]; 3],
    t_loc: [f64; 3],
    tt_ma: [[f64; 3]; 3],
    tt_lo: [f64; 3],
    ti_ma: [[f64; 3]; 3],
    ti_lo: [f64; 3],
    shapes: Vec<TopoShape>,
    e_map: ShapeIndexedMap,
    f_map: ShapeIndexedMap,
    algo: HlrAlgoPolyAlgo,
    debug: bool,
    tol_sta: f64,
    tol_end: f64,
    tol_angular: f64,
}

impl Default for HlrBRepPolyAlgo {
    /// `HLRBRep_PolyAlgo()` (`cxx:73-81`).
    fn default() -> Self {
        Self {
            proj: crate::hlr_algo::HlrAlgoProjector::default(),
            t_mat: [[0.0; 3]; 3],
            t_loc: [0.0; 3],
            tt_ma: [[0.0; 3]; 3],
            tt_lo: [0.0; 3],
            ti_ma: [[0.0; 3]; 3],
            ti_lo: [0.0; 3],
            shapes: Vec::new(),
            e_map: ShapeIndexedMap::new(),
            f_map: ShapeIndexedMap::new(),
            algo: HlrAlgoPolyAlgo::new(),
            debug: false,
            tol_sta: 0.1,
            tol_end: 0.9,
            tol_angular: 0.001,
        }
    }
}

impl HlrBRepPolyAlgo {
    /// `HLRBRep_PolyAlgo()` (`cxx:73-81`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `HLRBRep_PolyAlgo(const TopoDS_Shape&)` (`cxx:102-111`).
    pub fn from_shape(shape: TopoShape) -> Self {
        let mut a = Self::default();
        a.shapes.push(shape);
        a
    }

    /// `NbShapes()` (`hxx:109`).
    pub fn nb_shapes(&self) -> usize {
        self.shapes.len()
    }

    /// `Shape(I)` (`cxx:115-118`), 1-based.
    pub fn shape(&self, i: usize) -> Option<&TopoShape> {
        if i >= 1 {
            self.shapes.get(i - 1)
        } else {
            None
        }
    }

    /// `Shape(I)` mutable form (`cxx:115-118`).
    pub fn change_shape(&mut self, i: usize) -> Option<&mut TopoShape> {
        if i >= 1 {
            self.shapes.get_mut(i - 1)
        } else {
            None
        }
    }

    /// `Remove(I)` (`cxx:121-129`). The `Standard_OutOfRange` raise for an
    /// unknown index is `Err` here.
    pub fn remove(&mut self, i: usize) -> Result<(), &'static str> {
        if i == 0 || i > self.shapes.len() {
            return Err("HLRBRep_PolyAlgo::Remove : unknown Shape");
        }
        self.shapes.remove(i - 1);
        self.algo.clear();
        self.e_map.clear();
        self.f_map.clear();
        Ok(())
    }

    /// `Index(S)` (`cxx:133-147`): 1-based, 0 when absent.
    pub fn index(&self, s: &TopoShape) -> usize {
        let key = shape_map_key(s);
        self.shapes
            .iter()
            .position(|x| shape_map_key(x) == key)
            .map(|p| p + 1)
            .unwrap_or(0)
    }

    /// `Load(S)` (`hxx:121-124`).
    pub fn load(&mut self, shape: TopoShape) {
        self.shapes.push(shape);
    }

    /// `Algo()` (`hxx:126`).
    pub fn algo(&self) -> &HlrAlgoPolyAlgo {
        &self.algo
    }

    /// `Algo()` mutable form.
    pub fn change_algo(&mut self) -> &mut HlrAlgoPolyAlgo {
        &mut self.algo
    }

    /// `Projector()` (`hxx:131`).
    pub fn projector(&self) -> &crate::hlr_algo::HlrAlgoProjector {
        &self.proj
    }

    /// `Projector(theProj)` (`hxx:133`).
    pub fn set_projector(&mut self, proj: crate::hlr_algo::HlrAlgoProjector) {
        self.proj = proj;
    }

    /// `TolAngular()` (`hxx:135`).
    pub fn tol_angular(&self) -> f64 {
        self.tol_angular
    }

    /// `TolAngular(theTol)` (`hxx:137`).
    pub fn set_tol_angular(&mut self, tol: f64) {
        self.tol_angular = tol;
    }

    /// `TolCoef()` (`hxx:139`) — returns `myTolSta`.
    pub fn tol_coef(&self) -> f64 {
        self.tol_sta
    }

    /// `TolCoef(theTol)` (`hxx:141-144`): also stores `1 - theTol` in
    /// `myTolEnd`.
    pub fn set_tol_coef(&mut self, tol: f64) {
        self.tol_sta = tol;
        self.tol_end = 1.0 - tol;
    }

    /// `Debug()` (`hxx:179`).
    pub fn debug(&self) -> bool {
        self.debug
    }

    /// `Debug(theDebug)` (`hxx:181`).
    pub fn set_debug(&mut self, d: bool) {
        self.debug = d;
    }

    /// `TMultiply(theX, theY, theZ, theVPO)` (`cxx:4331-4340`): apply the
    /// projection matrix, adding `TLoc` unless `theVPO`.
    pub fn tmultiply(&self, x: &mut f64, y: &mut f64, z: &mut f64, vec_part_only: bool) {
        let add = if vec_part_only { 0.0 } else { 1.0 };
        let t = self.t_mat;
        let l = self.t_loc;
        let xt = t[0][0] * *x + t[0][1] * *y + t[0][2] * *z + add * l[0];
        let yt = t[1][0] * *x + t[1][1] * *y + t[1][2] * *z + add * l[1];
        *z = t[2][0] * *x + t[2][1] * *y + t[2][2] * *z + add * l[2];
        *x = xt;
        *y = yt;
    }

    /// `TMultiply(gp_XYZ&, theVPO)` (`hxx:581-588`).
    pub fn tmultiply_xyz(&self, p: &mut GpXyz, vec_part_only: bool) {
        let (mut x, mut y, mut z) = (p.x(), p.y(), p.z());
        self.tmultiply(&mut x, &mut y, &mut z, vec_part_only);
        *p = GpXyz::new(x, y, z);
    }

    /// `TTMultiply(theX, theY, theZ, theVPO)` (`cxx:4344-4352`).
    pub fn ttmultiply(&self, x: &mut f64, y: &mut f64, z: &mut f64, vec_part_only: bool) {
        let add = if vec_part_only { 0.0 } else { 1.0 };
        let t = self.tt_ma;
        let l = self.tt_lo;
        let xt = t[0][0] * *x + t[0][1] * *y + t[0][2] * *z + add * l[0];
        let yt = t[1][0] * *x + t[1][1] * *y + t[1][2] * *z + add * l[1];
        *z = t[2][0] * *x + t[2][1] * *y + t[2][2] * *z + add * l[2];
        *x = xt;
        *y = yt;
    }

    /// `TTMultiply(gp_XYZ&, theVPO)` (`hxx:592-598`).
    pub fn ttmultiply_xyz(&self, p: &mut GpXyz, vec_part_only: bool) {
        let (mut x, mut y, mut z) = (p.x(), p.y(), p.z());
        self.ttmultiply(&mut x, &mut y, &mut z, vec_part_only);
        *p = GpXyz::new(x, y, z);
    }

    /// `TIMultiply(theX, theY, theZ, theVPO)` (`cxx:4357-4365`).
    pub fn timultiply(&self, x: &mut f64, y: &mut f64, z: &mut f64, vec_part_only: bool) {
        let add = if vec_part_only { 0.0 } else { 1.0 };
        let t = self.ti_ma;
        let l = self.ti_lo;
        let xt = t[0][0] * *x + t[0][1] * *y + t[0][2] * *z + add * l[0];
        let yt = t[1][0] * *x + t[1][1] * *y + t[1][2] * *z + add * l[1];
        *z = t[2][0] * *x + t[2][1] * *y + t[2][2] * *z + add * l[2];
        *x = xt;
        *y = yt;
    }

    /// `TIMultiply(gp_XYZ&, theVPO)` (`hxx:604-612`).
    pub fn timultiply_xyz(&self, p: &mut GpXyz, vec_part_only: bool) {
        let (mut x, mut y, mut z) = (p.x(), p.y(), p.z());
        self.timultiply(&mut x, &mut y, &mut z, vec_part_only);
        *p = GpXyz::new(x, y, z);
    }

    /// `cxx:302-341`: read the projection transform and its inverse into
    /// `TMat`/`TLoc` and `TIMa`/`TILo`. In OCCT this is the prologue of
    /// `StoreShell`; the port keeps it as a named seam until `StoreShell`
    /// (UNPORTED) calls it.
    #[allow(dead_code)]
    fn read_projection_matrices(&mut self) {
        copy_trsf(self.proj.transformation(), &mut self.t_mat, &mut self.t_loc);
        copy_trsf(
            self.proj.inverted_transformation(),
            &mut self.ti_ma,
            &mut self.ti_lo,
        );
    }

    /// `cxx:429-446`: `aTT = aLoc.Transformation(); aTT.PreMultiply(aProjTrsf)`
    /// read into `TTMa`/`TTLo`. Called by `StoreShell`'s face loop and by
    /// `InsertOnOutLine`'s prologue.
    fn read_face_matrices(&mut self, loc: &GpTrsf) {
        // `gp_Trsf::PreMultiply` (`gp_Trsf.cxx`): `this = aT * this`, i.e. the
        // location applied after the projection.
        let mut a_tt = loc.clone();
        a_tt.pre_multiply(self.proj.transformation());
        copy_trsf(&a_tt, &mut self.tt_ma, &mut self.tt_lo);
    }

    /// `MakeShape()` (`cxx:204-220`): a compound of every loaded shape, `None`
    /// when nothing was loaded.
    pub fn make_shape(&self) -> Option<TopoShape> {
        if self.shapes.is_empty() {
            return None;
        }
        let mut comp = TopoShape::new(ShapeType::Compound);
        let builder = crate::builder::TopoBuilder::new();
        for s in &self.shapes {
            builder.add(&mut comp, s);
        }
        Some(comp)
    }

    /// `InitShape(theShape, theIsoledF, theIsoledE)` (`cxx:224-287`): count the
    /// shells that carry a triangulation, detect isolated faces and isolated
    /// edges, and size the poly-shell array. Returns
    /// `(nb_shell, isoled_f, isoled_e)`.
    pub fn init_shape(&mut self, shape: &TopoShape) -> (i32, bool, bool) {
        let mut shape_map0: HashSet<ShapeMapKey> = HashSet::new();
        let mut nb_shell = 0i32;

        // `cxx:234-251`.
        for shell in explorer_shapes(shape, ShapeType::Shell) {
            let mut has_trian = false;
            for face in explorer_shapes(&shell, ShapeType::Face) {
                if has_triangulation(&face) && shape_map0.insert(shape_map_key(&face)) {
                    has_trian = true;
                }
            }
            if has_trian {
                nb_shell += 1;
            }
        }

        // `cxx:254-267`: `More() && !theIsoledF`.
        let mut isoled_f = false;
        for face in explorer_shapes_avoiding(shape, ShapeType::Face, ShapeType::Shell) {
            if isoled_f {
                break;
            }
            if has_triangulation(&face) && shape_map0.insert(shape_map_key(&face)) {
                isoled_f = true;
            }
        }
        if isoled_f {
            nb_shell += 1;
        }

        // `cxx:270-279`: the first `TopAbs_EDGE` occurrence outside a face.
        let isoled_e =
            !explorer_shapes_avoiding(shape, ShapeType::Edge, ShapeType::Face).is_empty();
        if isoled_e {
            nb_shell += 1;
        }

        if nb_shell > 0 {
            self.algo.init(nb_shell as usize);
        }
        (nb_shell, isoled_f, isoled_e)
    }

    /// `InitHide()` (`hxx:150`).
    pub fn init_hide(&mut self) {
        self.algo.init_hide();
    }

    /// `MoreHide()` (`hxx:152`).
    pub fn more_hide(&self) -> bool {
        self.algo.more_hide()
    }

    /// `NextHide()` (`hxx:154`).
    pub fn next_hide(&mut self) {
        self.algo.next_hide();
    }

    /// `Hide(status, S, reg1, regn, outl, intl)` (`cxx:4370-4383`): the shape is
    /// `myFMap(anIndex)` for an internal outline and `myEMap(anIndex)` otherwise.
    pub fn hide(&mut self, status: &mut HlrAlgoEdgeStatus) -> Option<PolyAlgoSegment> {
        let info = self.algo.hide(status);
        let map = if info.intl { &self.f_map } else { &self.e_map };
        let shape = map.find_key(info.index as usize)?.clone();
        Some(PolyAlgoSegment {
            points: info.points,
            shape,
            reg1: info.reg1,
            regn: info.regn,
            outl: info.outl,
            intl: info.intl,
        })
    }

    /// `InitShow()` (`hxx:164`).
    pub fn init_show(&mut self) {
        self.algo.init_show();
    }

    /// `MoreShow()` (`hxx:166`).
    pub fn more_show(&self) -> bool {
        self.algo.more_show()
    }

    /// `NextShow()` (`hxx:168`).
    pub fn next_show(&mut self) {
        self.algo.next_show();
    }

    /// `Show(S, reg1, regn, outl, intl)` (`cxx:4386-4397`).
    pub fn show(&self) -> Option<PolyAlgoSegment> {
        let info = self.algo.show();
        let map = if info.intl { &self.f_map } else { &self.e_map };
        let shape = map.find_key(info.index as usize)?.clone();
        Some(PolyAlgoSegment {
            points: info.points,
            shape,
            reg1: info.reg1,
            regn: info.regn,
            outl: info.outl,
            intl: info.intl,
        })
    }
}

/// `BRep_Tool::Triangulation(aFace, aLoc).IsNull()` negated
/// (`BRep_Tool.cxx:112-122`): whether the face carries any triangulation. The
/// storage accessor is `TShape::face_triangulation` (`tshape.rs:133`).
pub(crate) fn has_triangulation(shape: &TopoShape) -> bool {
    shape.tshape.read().unwrap().face_triangulation().is_some()
}

/// `gp_Trsf` -> `(matrix, translation)` for the `TMat`/`TLoc` trios
/// (`cxx:303-341`). `gp_Trsf::VectorialPart` folds the scale in.
fn copy_trsf(t: &GpTrsf, mat: &mut [[f64; 3]; 3], loc: &mut [f64; 3]) {
    let v = t.translation_part();
    loc[0] = v.x();
    loc[1] = v.y();
    loc[2] = v.z();
    let m = t.vectorial_part();
    for (r, row) in mat.iter_mut().enumerate() {
        for (c, cell) in row.iter_mut().enumerate() {
            *cell = m.value(r + 1, c + 1);
        }
    }
}

/// `aCoef2 * uv1 + aCoef1 * uv2` (`cxx:3579`, `:3605`).
pub(crate) fn xy_scaled_sum(
    a: &occt_core::gp::GpXY,
    ka: f64,
    b: &occt_core::gp::GpXY,
    kb: f64,
) -> occt_core::gp::GpXY {
    occt_core::gp::GpXY::new(ka * a.x() + kb * b.x(), ka * a.y() + kb * b.y())
}

/// The `aNorm > 0` / `(1 / aNorm) * aXYZ` / `(1, 0, 0)` arm of
/// `HLRAlgo_PolyInternalData::AddNode` (`cxx:462-479`) and `ChangeNode`
/// (`cxx:3588-3624`).
pub(crate) fn blended_normal(a: &GpXyz, ka: f64, b: &GpXyz, kb: f64) -> GpXyz {
    let x = GpXyz::new(
        ka * a.x() + kb * b.x(),
        ka * a.y() + kb * b.y(),
        ka * a.z() + kb * b.z(),
    );
    let norm = x.modulus();
    if norm > 0.0 {
        x.multiplied(1.0 / norm)
    } else {
        GpXyz::new(1.0, 0.0, 0.0)
    }
}
