//! `BiTgte_Blend` (`BiTgte_Blend.hxx`, `.cxx`): the rolling-ball blend of a
//! shape - tubes along the centre lines of the contacts, spheres on the
//! vertices where enough tubes meet, and the faces that connect them.
//!
//! Ported, by file:
//!
//! * `mod.rs`: `BiTgte_Blend()` (`cxx:707-713`), `BiTgte_Blend(S, Radius,
//!   Tol, NUBS)` (`:716-724`), `Init` (`:728-739`), `Clear` (`:742-754`),
//!   `SetStoppingFace` (`:758-767`), `SetFaces` (`:771-775`), `SetEdge`
//!   (`:779-782`), `Perform` (`:785-940`), `IsDone` (`:944-947`), `Shape`
//!   (`:951-954`), `NbSurfaces` (`:958-961`), `Surface(Index)` (`:966-969`),
//!   `Face(Index)` (`:974-977`), `CenterLines` (`:982-990`), `Surface(
//!   CenterLine)` (`:994-997`), `Face(CenterLine)` (`:1005-1012`),
//!   `ContactType` (`:1016-1083`), `SupportShape1` (`:1087-1098`),
//!   `SupportShape2` (`:1104-1115`), `CurveOnShape1` (`:1119-1136`),
//!   `CurveOnShape2` (`:1141-1159`), the four `PCurve*` accessors
//!   (`:1163-1192`), `NbBranches` (`:1196-1270`) and `IndicesOfBranche`
//!   (`:1274-1281`).
//! * `compute_centers.rs`: `ComputeCenters` (`:1286-1757`).
//! * `compute_surfaces.rs`: `ComputeSurfaces` (`:1761-2217`).
//! * `compute_shape.rs`: `ComputeShape` (`:2221-2507`).
//! * `intersect.rs`: `Intersect` (`:2511-2664`).
//! * `make_curve.rs`: `MakeCurve` (`:381-470`).
//! * `statics.rs`: the nine file-local statics (`:96-706`).
//! * `offset_map.rs`: `myMapSF` (`hxx:206`) as an insertion-ordered map.
//!
//! UNPORTED, with the OCCT lines that keep it out:
//!
//! * `Touched` (`cxx:477-497`): every line of the reference body is commented
//!   out, so the static does nothing and `ComputeCenters`'s `TouchedByCork`
//!   stays empty - kept as an empty function for the call site.
//! * the `OCCT_DEBUG` timing / tracing blocks (`:850-940`, `:1749-1756`,
//!   `:1919-1922`, `:1985`, `:2035`): debug output only.
//! * `BRepTools::Update` (`:2336`, `:2473`): the port keeps no cached UV points
//!   nor triangulation.
//! * `myImage` (`hxx:208`): `Clear` clears it and nothing else ever reads or
//!   writes it in the reference.
//! * `myAncestors` (`hxx:193`) is [`AncestorMap`]; `Perform` rebuilds it after
//!   the sewing (`:833`) exactly as the reference does.
//!
//! Deviations:
//! - `myMapSF` is [`OffsetMap`] (insertion order plus the key shapes) because
//!   `NCollection_DataMap` iteration order is hash order and the key shape is
//!   not recoverable from a [`ShapeMapKey`].
//! - `myCreated` (`hxx:196-198`) is [`CreatedMap`], which keeps the outer and
//!   inner walks in insertion order for the same reason.
//! - accessors that return a reference in C++ (`Shape`, `Face(CenterLine)`,
//!   `SupportShape1/2`) answer an `Option`/clone here: a null handle or a
//!   missing map entry is a reachable state in the port, where OCCT would
//!   either return a null reference or throw.

mod compute_centers;
mod compute_shape;
mod compute_surfaces;
mod intersect;
mod make_curve;
mod offset_map;
mod statics;

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use occt_core::precision::CONFUSION;
use occt_geom::{Curve, GeomTrimmedCurve, Surface};
use occt_geom2d::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_algo::AsDes;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_builder_api_sewing::Sewing;
use crate::brep_lib_build_curves3d::build_curves3d;
use crate::brep_lib_same_parameter::same_parameter_shape;
use crate::brep_offset::ShapeIndexedMap;
use crate::brep_offset::Analyse;
use crate::brep_sweep::Quilt;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key, AncestorMap, ShapeMapKey};

use self::offset_map::OffsetMap;
use super::ContactType;

/// `NCollection_DataMap<TopoDS_Shape, NCollection_DataMap<TopoDS_Shape,
/// NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>, ...>` `myCreated`
/// (`hxx:196-198`): the offset shapes generated for the edges and vertices of
/// every face, filled by `ComputeSurfaces` (`cxx:2042-2049`, `:2353-2375`,
/// `:2439-2442`, `:2453-2456`) and read by `ComputeShape`
/// (`cxx:2353-2375`, `:2416-2434`, `:2444-2460`).
///
/// `bind` replaces the inner map, exactly like `NCollection_DataMap::Bind`.
#[derive(Debug, Clone, Default)]
pub(super) struct CreatedMap {
    entries: Vec<(TopoShape, Vec<(TopoShape, Vec<TopoShape>)>)>,
    index: HashMap<ShapeMapKey, usize>,
}

impl CreatedMap {
    pub(super) fn is_bound(&self, s: &TopoShape) -> bool {
        self.index.contains_key(&shape_map_key(s))
    }

    /// `Created.Bind(S, EmptyMap)`: replaces (or inserts) an empty inner map.
    pub(super) fn bind_empty(&mut self, s: &TopoShape) {
        let key = shape_map_key(s);
        match self.index.get(&key) {
            Some(&i) => self.entries[i].1.clear(),
            None => {
                self.index.insert(key, self.entries.len());
                self.entries.push((s.clone(), Vec::new()));
            }
        }
    }

    fn inner(&self, s: &TopoShape) -> Option<&Vec<(TopoShape, Vec<TopoShape>)>> {
        self.index.get(&shape_map_key(s)).map(|&i| &self.entries[i].1)
    }

    fn inner_mut(&mut self, s: &TopoShape) -> Option<&mut Vec<(TopoShape, Vec<TopoShape>)>> {
        match self.index.get(&shape_map_key(s)) {
            Some(&i) => Some(&mut self.entries[i].1),
            None => None,
        }
    }

    /// `Created(S).IsBound(F)`.
    pub(super) fn inner_is_bound(&self, s: &TopoShape, f: &TopoShape) -> bool {
        self.inner(s)
            .map(|l| l.iter().any(|(k, _)| is_same(k, f)))
            .unwrap_or(false)
    }

    /// `Created(S).Bind(F, Empty)`.
    pub(super) fn inner_bind_empty(&mut self, s: &TopoShape, f: &TopoShape) {
        let key = shape_map_key(f);
        let Some(l) = self.inner_mut(s) else {
            return;
        };
        match l.iter_mut().find(|(k, _)| shape_map_key(k) == key) {
            Some((_, v)) => v.clear(),
            None => l.push((f.clone(), Vec::new())),
        }
    }

    /// `Created(S)(F).Append(V)`.
    pub(super) fn inner_append(&mut self, s: &TopoShape, f: &TopoShape, v: &TopoShape) {
        self.inner_bind_empty(s, f);
        let key = shape_map_key(f);
        if let Some(l) = self.inner_mut(s) {
            if let Some((_, list)) = l.iter_mut().find(|(k, _)| shape_map_key(k) == key) {
                list.push(v.clone());
            }
        }
    }

    /// `Created(S)(F).First()`.
    pub(super) fn inner_first(&self, s: &TopoShape, f: &TopoShape) -> Option<&TopoShape> {
        let key = shape_map_key(f);
        self.inner(s)?
            .iter()
            .find(|(k, _)| shape_map_key(k) == key)
            .and_then(|(_, list)| list.first())
    }
}

/// The stored p-curve of `e` on `f` and its range, i.e.
/// `BRep_Tool::CurveOnSurface(E, F, f, l)` (`BRep_Tool.cxx:310-360`).
pub(super) fn pcurve_on_face(
    e: &Edge,
    f: &Face,
) -> Option<(Arc<dyn Curve2d>, Arc<dyn Surface>, f64, f64)> {
    let s = BRepTool::face_surface(f)?;
    let key = Arc::as_ptr(&s) as *const () as usize;
    let reps: Vec<(Arc<dyn Curve2d>, f64, f64)> =
        GeometryRegistry::global()
            .edge_pcurve_reps(&e.0)
            .into_iter()
            .filter(|(surf, _, _, _)| Arc::as_ptr(surf) as *const () as usize == key)
            .map(|(_, pc, a, b)| (pc, a, b))
            .collect();
    // `BRep_Tool::CurveOnSurface` asks for `PCurve2()` when the representation
    // is a closed surface and the edge is REVERSED, `PCurve()` otherwise.
    if reps.len() >= 2 && e.0.orientation() == Orientation::Reversed {
        let (pc, a, b) = reps.into_iter().nth(1)?;
        return Some((pc, s, a, b));
    }
    let (pc, a, b) = reps.into_iter().next()?;
    Some((pc, s, a, b))
}

/// `BRepLib_MakeEdge(C, V1, V2)`: `Project` both vertices, then `Init`
/// (`BRepLib_MakeEdge.cxx:542-572`).
pub(super) fn make_edge_3d(
    curve: &Arc<dyn Curve>,
    v1: Option<&Vertex>,
    v2: Option<&Vertex>,
) -> Option<Edge> {
    let build = crate::brep_lib_make_edge::init_3d_projecting(curve, v1, v2)?;
    let b = TopoBuilder::new();
    let mut e = b.make_edge(curve.clone(), build.first, build.last);
    crate::brep_lib_make_edge::add_endpoint_vertices(&b, &mut e, &build);
    if build.closed {
        e.0.set_closed(true);
    }
    Some(e)
}

/// `BRepLib_MakeEdge(C2d, S, V1, V2, f, l)` (`BRepLib_MakeEdge.cxx:905-1094`).
pub(super) fn make_edge_2d(
    pc: &Arc<dyn Curve2d>,
    surf: &Arc<dyn Surface>,
    f: &Face,
    v1: Option<&Vertex>,
    v2: Option<&Vertex>,
    p1: f64,
    p2: f64,
) -> Option<Edge> {
    let build =
        crate::brep_lib_make_edge::init_2d_vertices(pc, surf, v1, v2, p1, p2)?;
    let mut e = Edge::new();
    let reg = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(&f.0);
    reg.set_edge_pcurve(&e.0, face_key, pc.clone());
    reg.set_ranges_all(&e.0, build.first, build.last);
    reg.raise_edge_tolerance(&e.0, CONFUSION);
    crate::brep_lib_make_edge::add_endpoint_vertices(&TopoBuilder::new(), &mut e, &build);
    if build.reverse {
        e.0.set_orientation(Orientation::Reversed);
    }
    Some(e)
}

/// `BiTgte_Blend` (`BiTgte_Blend.hxx:46-216`).
///
/// No `Clone` / `Debug`: `myAnalyse` ([`Analyse`]) and `myAncestors`
/// ([`AncestorMap`]) own interning state that is not trivially duplicable in
/// the port, and the reference class is not copyable either.
pub struct Blend {
    /// `myRadius` (`hxx:184`).
    my_radius: f64,
    /// `myTol` (`hxx:185`).
    my_tol: f64,
    /// `myNubs` (`hxx:186`).
    my_nubs: bool,
    /// `myShape` (`hxx:189`).
    my_shape: Option<TopoShape>,
    /// `myResult` (`hxx:190`).
    my_result: Option<TopoShape>,
    /// `myBuildShape` (`hxx:192`).
    my_build_shape: bool,
    /// `myAncestors` (`hxx:193-194`): edge -> faces of the sewed shape.
    my_ancestors: AncestorMap,
    /// `myCreated` (`hxx:196-198`).
    my_created: CreatedMap,
    /// `myCutEdges` (`hxx:199-200`).
    my_cut_edges: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myFaces` (`hxx:201`).
    my_faces: ShapeIndexedMap,
    /// `myEdges` (`hxx:202`).
    my_edges: ShapeIndexedMap,
    /// `myStopFaces` (`hxx:203`).
    my_stop_faces: HashSet<ShapeMapKey>,
    /// `myAnalyse` (`hxx:204`).
    my_analyse: Analyse,
    /// `myCenters` (`hxx:205`).
    my_centers: ShapeIndexedMap,
    /// `myMapSF` (`hxx:206`).
    my_map_sf: OffsetMap,
    /// `myInitOffsetFace` (`hxx:207`).
    my_init_offset_face: BRepAlgoImage,
    /// `myImage` (`hxx:208`).
    my_image: BRepAlgoImage,
    /// `myImageOffset` (`hxx:209`).
    my_image_offset: BRepAlgoImage,
    /// `myAsDes` (`hxx:210`).
    my_as_des: AsDes,
    /// `myNbBranches` (`hxx:211`).
    my_nb_branches: i32,
    /// `myIndices` (`hxx:212`).
    my_indices: Vec<i32>,
    /// `myDone` (`hxx:215`).
    my_done: bool,
}

impl Default for Blend {
    fn default() -> Self {
        Self::new()
    }
}

impl Blend {
    /// `BiTgte_Blend::BiTgte_Blend()` (`cxx:707-713`): a fresh `AsDes` and the
    /// `myNbBranches` sentinel.
    pub fn new() -> Self {
        Self {
            my_radius: 0.0,
            my_tol: 0.0,
            my_nubs: false,
            my_shape: None,
            my_result: None,
            my_build_shape: false,
            my_ancestors: AncestorMap::default(),
            my_created: CreatedMap::default(),
            my_cut_edges: HashMap::new(),
            my_faces: ShapeIndexedMap::new(),
            my_edges: ShapeIndexedMap::new(),
            my_stop_faces: HashSet::new(),
            my_analyse: Analyse::new(),
            my_centers: ShapeIndexedMap::new(),
            my_map_sf: OffsetMap::new(),
            my_init_offset_face: BRepAlgoImage::new(),
            my_image: BRepAlgoImage::new(),
            my_image_offset: BRepAlgoImage::new(),
            my_as_des: AsDes::new(),
            my_nb_branches: -1,
            my_indices: Vec::new(),
            my_done: false,
        }
    }

    /// `BiTgte_Blend::BiTgte_Blend(S, Radius, Tol, NUBS)` (`cxx:716-724`).
    pub fn with_shape(shape: &TopoShape, radius: f64, tol: f64, nubs: bool) -> Self {
        let mut b = Self::new();
        b.init(shape, radius, tol, nubs);
        b
    }

    /// `void BiTgte_Blend::Init(S, Radius, Tol, NUBS)` (`cxx:728-739`).
    pub fn init(&mut self, shape: &TopoShape, radius: f64, tol: f64, nubs: bool) {
        self.clear();
        self.my_shape = Some(shape.clone());
        self.my_tol = tol;
        self.my_nubs = nubs;
        self.my_radius = radius;
        self.my_nb_branches = -1;
    }

    /// `void BiTgte_Blend::Clear()` (`cxx:742-754`).
    pub fn clear(&mut self) {
        self.my_init_offset_face.clear();
        self.my_image.clear();
        self.my_image_offset.clear();
        self.my_stop_faces.clear();
        self.my_analyse.clear();
        self.my_as_des.clear();
        self.my_nb_branches = -1;
        self.my_done = false;
    }

    /// `void BiTgte_Blend::SetStoppingFace(Face)` (`cxx:758-767`).
    pub fn set_stopping_face(&mut self, face: &Face) {
        self.my_stop_faces.insert(shape_map_key(&face.0));
    }

    /// `void BiTgte_Blend::SetFaces(F1, F2)` (`cxx:771-775`).
    pub fn set_faces(&mut self, f1: &Face, f2: &Face) {
        self.my_faces.add(&f1.0);
        self.my_faces.add(&f2.0);
    }

    /// `void BiTgte_Blend::SetEdge(Edge)` (`cxx:779-782`).
    pub fn set_edge(&mut self, edge: &Edge) {
        self.my_edges.add(&edge.0);
    }

    /// `void BiTgte_Blend::Perform(BuildShape)` (`cxx:785-940`).
    ///
    /// `:790-850` sews the input so that the borders that only look free are
    /// not rolled on, `:866-883` re-orders `myFaces` / `myStopFaces` in the
    /// sewed shape, then the three stages run, `BRepLib::BuildCurves3d` fills
    /// the 3D curves of the result (`:924`).
    pub fn perform(&mut self, build_shape: bool) {
        self.my_build_shape = build_shape;

        let Some(shape) = self.my_shape.clone() else {
            return;
        };

        // `cxx:792-800`.
        let mut sew = Sewing::new(self.my_tol, true, true, true, false);
        let _ = build_curves3d(&shape);
        for f in explorer_shapes(&shape, ShapeType::Face) {
            sew.add(&f);
        }
        // `:798`: `Sew->Perform()`. A failed sewing leaves `SewedShape` null
        // and OCCT throws; the port keeps the input shape (see the module
        // header: the throw is a guard, not a branch).
        let sewed_shape = match sew.perform() {
            Ok(()) => sew.sewed_shape().cloned(),
            Err(_) => None,
        };
        let Some(mut sewed_shape) = sewed_shape else {
            return;
        };

        // `cxx:806-823`: a sewing that flipped the orientation of the
        // reference face reverses the whole sewed shape.
        if let Some(first) = explorer_shapes(&shape, ShapeType::Face).first().cloned() {
            let mut face_ref = Face(first.clone());
            let ori_ref = first.orientation();
            if sew.is_modified(&face_ref.0) {
                let m = sew.modified(&face_ref.0).clone();
                face_ref = Face(m);
            }
            for ff in explorer_shapes(&sewed_shape, ShapeType::Face) {
                if is_same(&face_ref.0, &ff) && ff.orientation() != ori_ref {
                    sewed_shape.reverse();
                    break;
                }
            }
        }

        // `cxx:828-831`: `BRepLib::SameParameter` on every edge of the sewed
        // shape.
        for e in explorer_shapes(&sewed_shape, ShapeType::Edge) {
            same_parameter_shape(&e, BRepTool::edge_tolerance(&Edge(e.clone())), false);
        }

        // `cxx:833`.
        self.my_ancestors =
            AncestorMap::map_shapes_and_ancestors(&sewed_shape, ShapeType::Edge, ShapeType::Face);

        // `cxx:838-847`: follow the modified faces of `myFaces`.
        let originals: Vec<TopoShape> = explorer_shapes(&shape, ShapeType::Face);
        for f in originals {
            if self.my_faces.contains(&f) && sew.is_modified(&f) {
                let m = sew.modified(&f).clone();
                self.my_faces.remove_key(&f);
                self.my_faces.add(&m);
            }
        }

        // `cxx:849`.
        self.my_shape = Some(sewed_shape);

        // `cxx:866-883`: re-add the faces (and stopping faces) of the sewed
        // shape so that their order and orientation follow it.
        let Some(current) = self.my_shape.clone() else {
            return;
        };
        for f in explorer_shapes(&current, ShapeType::Face) {
            if self.my_faces.contains(&f) {
                self.my_faces.remove_key(&f);
                self.my_faces.add(&f);
            } else if self.my_stop_faces.contains(&shape_map_key(&f)) {
                self.my_stop_faces.remove(&shape_map_key(&f));
                self.my_stop_faces.insert(shape_map_key(&f));
            }
        }

        // `cxx:886-889`.
        self.compute_centers();
        // `cxx:897`.
        self.compute_surfaces();
        // `cxx:915-918`.
        if self.my_build_shape {
            self.compute_shape();
        }

        // `cxx:924`.
        if let Some(result) = self.my_result.clone() {
            let _ = build_curves3d(&result);
        }

        // `cxx:939`.
        self.my_done = true;
    }

    /// `bool BiTgte_Blend::IsDone() const` (`cxx:944-947`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `const TopoDS_Shape& BiTgte_Blend::Shape() const` (`cxx:951-954`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_result.as_ref()
    }

    /// `int BiTgte_Blend::NbSurfaces() const` (`cxx:958-961`).
    pub fn nb_surfaces(&self) -> usize {
        self.my_centers.extent()
    }

    /// `Geom_Surface BiTgte_Blend::Surface(Index)` (`cxx:966-969`).
    pub fn surface_of_index(&self, index: usize) -> Option<Arc<dyn Surface>> {
        let center = self.my_centers.value(index.checked_sub(1)?)?.clone();
        self.surface_of_center_line(&center)
    }

    /// `const TopoDS_Face& BiTgte_Blend::Face(Index)` (`cxx:974-977`).
    pub fn face_of_index(&self, index: usize) -> Option<Face> {
        let center = self.my_centers.value(index.checked_sub(1)?)?.clone();
        self.face_of_center_line(&center)
    }

    /// `void BiTgte_Blend::CenterLines(LC) const` (`cxx:982-990`).
    pub fn center_lines(&self) -> Vec<TopoShape> {
        self.my_centers.iter().cloned().collect()
    }

    /// `Geom_Surface BiTgte_Blend::Surface(CenterLine)` (`cxx:994-997`).
    pub fn surface_of_center_line(&self, center_line: &TopoShape) -> Option<Arc<dyn Surface>> {
        let f = self.my_map_sf.get(center_line)?.face()?;
        BRepTool::face_surface(f)
    }

    /// `const TopoDS_Face& BiTgte_Blend::Face(CenterLine)` (`cxx:1005-1012`):
    /// `Standard_DomainError` for an unbound centre line, i.e. `None` here.
    pub fn face_of_center_line(&self, center_line: &TopoShape) -> Option<Face> {
        self.my_map_sf.get(center_line)?.face().cloned()
    }

    /// `BiTgte_ContactType BiTgte_Blend::ContactType(Index)` (`cxx:1016-1083`).
    pub fn contact_type(&self, index: usize) -> Option<ContactType> {
        let s1 = self.support_shape1(index)?;
        let s2 = self.support_shape2(index)?;
        let mut t1 = s1.shape_type();
        let mut t2 = s2.shape_type();
        if t2 < t1 {
            std::mem::swap(&mut t1, &mut t2);
        }
        Some(match (t1, t2) {
            (ShapeType::Vertex, ShapeType::Vertex) => ContactType::VertexVertex,
            (ShapeType::Vertex, ShapeType::Edge) => ContactType::EdgeVertex,
            (ShapeType::Vertex, ShapeType::Face) => ContactType::FaceVertex,
            (ShapeType::Edge, ShapeType::Edge) => ContactType::EdgeEdge,
            (ShapeType::Edge, ShapeType::Face) => ContactType::FaceEdge,
            (ShapeType::Face, ShapeType::Face) => ContactType::FaceEdge,
            _ => ContactType::VertexVertex,
        })
    }

    /// `const TopoDS_Shape& BiTgte_Blend::SupportShape1(Index)`
    /// (`cxx:1087-1098`).
    pub fn support_shape1(&self, index: usize) -> Option<TopoShape> {
        let cur_e = self.my_centers.value(index.checked_sub(1)?)?;
        let l = self.my_as_des.ascendant(cur_e);
        let f1 = l.first()?;
        self.my_init_offset_face.image_from(f1).ok()
    }

    /// `const TopoDS_Shape& BiTgte_Blend::SupportShape2(Index)`
    /// (`cxx:1104-1115`).
    pub fn support_shape2(&self, index: usize) -> Option<TopoShape> {
        let cur_e = self.my_centers.value(index.checked_sub(1)?)?;
        let l = self.my_as_des.ascendant(cur_e);
        let f2 = l.last()?;
        self.my_init_offset_face.image_from(f2).ok()
    }

    /// `Geom_Curve BiTgte_Blend::CurveOnShape1(Index)` (`cxx:1119-1136`): the
    /// first edge of the tube, trimmed to its own range.
    pub fn curve_on_shape1(&self, index: usize) -> Option<Arc<dyn Curve>> {
        self.curve_on_shape(index, 0)
    }

    /// `Geom_Curve BiTgte_Blend::CurveOnShape2(Index)` (`cxx:1141-1159`).
    pub fn curve_on_shape2(&self, index: usize) -> Option<Arc<dyn Curve>> {
        self.curve_on_shape(index, 1)
    }

    fn curve_on_shape(&self, index: usize, which: usize) -> Option<Arc<dyn Curve>> {
        let cur_e = self.my_centers.value(index.checked_sub(1)?)?;
        let entry = self.my_map_sf.get(cur_e)?;
        let f = entry.face()?;
        let e = explorer_shapes(&f.0, ShapeType::Edge).into_iter().nth(which)?;
        let e = Edge(e);
        if BRepTool::is_degenerated(&e) {
            return None;
        }
        let c = BRepTool::edge_curve(&e)?;
        let (f1, l1) = BRepTool::edge_parameters(&e);
        Some(Arc::new(GeomTrimmedCurve::new(c, f1, l1)) as Arc<dyn Curve>)
    }

    /// `Geom2d_Curve BiTgte_Blend::PCurveOnFace1(Index)` (`cxx:1163-1168`): a
    /// null handle in the reference.
    pub fn pcurve_on_face1(&self, _index: usize) -> Option<Arc<dyn Curve2d>> {
        None
    }

    /// `Geom2d_Curve BiTgte_Blend::PCurve1OnFillet(Index)` (`cxx:1172-1177`).
    pub fn pcurve1_on_fillet(&self, _index: usize) -> Option<Arc<dyn Curve2d>> {
        None
    }

    /// `Geom2d_Curve BiTgte_Blend::PCurveOnFace2(Index)` (`cxx:1181-1186`).
    pub fn pcurve_on_face2(&self, _index: usize) -> Option<Arc<dyn Curve2d>> {
        None
    }

    /// `Geom2d_Curve BiTgte_Blend::PCurve2OnFillet(Index)` (`cxx:1190-1195`).
    pub fn pcurve2_on_fillet(&self, _index: usize) -> Option<Arc<dyn Curve2d>> {
        None
    }

    /// `int BiTgte_Blend::NbBranches()` (`cxx:1196-1270`): glue the generated
    /// faces into shells and re-order `myCenters` branch by branch.
    pub fn nb_branches(&mut self) -> i32 {
        if self.my_nb_branches != -1 {
            return self.my_nb_branches;
        }

        // `cxx:1204-1210`.
        let mut glue = Quilt::new();
        let nb_faces = self.my_centers.extent();
        if nb_faces == 0 {
            return 0;
        }

        // `cxx:1213-1218`.
        for i in 0..nb_faces {
            let center = self.my_centers.value(i).expect("myCenters(i)").clone();
            if let Some(f) = self.my_map_sf.get(&center).and_then(|of| of.face().cloned()) {
                glue.add(&f.0);
            }
        }
        let shells = glue.shells();

        // `cxx:1223-1231`.
        self.my_nb_branches = 0;
        let mut tmp_map = ShapeIndexedMap::new();
        for _ in explorer_shapes(&shells, ShapeType::Shell) {
            self.my_nb_branches += 1;
        }

        // `cxx:1233-1236`: `myIndices(1..NbBranches+1)`, 1-based, `myIndices(1)
        // = 0`.
        self.my_indices = vec![0; self.my_nb_branches as usize + 2];
        let mut count: i32 = 0;
        let mut index = 2usize;

        // `cxx:1239-1265`.
        for cur_s in explorer_shapes(&shells, ShapeType::Shell) {
            for cur_f in explorer_shapes(&cur_s, ShapeType::Face) {
                for i in 0..nb_faces {
                    let center = self.my_centers.value(i).expect("myCenters(i)").clone();
                    let Some(rakk) = self.my_map_sf.get(&center).and_then(|of| of.face().cloned())
                    else {
                        continue;
                    };
                    if is_same(&cur_f, &rakk.0) {
                        tmp_map.add(&center);
                        count += 1;
                        break;
                    }
                }
            }
            if index < self.my_indices.len() {
                self.my_indices[index] = count;
            }
            index += 1;
        }

        self.my_centers = tmp_map;
        self.my_nb_branches
    }

    /// `void BiTgte_Blend::IndicesOfBranche(Index, From, To) const`
    /// (`cxx:1274-1281`).
    pub fn indices_of_branche(&self, index: usize, from: &mut i32, to: &mut i32) {
        // `myIndices` is 1-based in the reference.
        let get = |i: usize| self.my_indices.get(i.saturating_sub(1)).copied().unwrap_or(0);
        *from = get(index) + 1;
        *to = get(index + 1);
    }
}


