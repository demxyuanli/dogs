//! Construction and `Init` (`ChFi2d_Builder.cxx:66-239`), including the
//! `IsIssuedFrom` helper (`cxx:208-239`).

use std::collections::HashSet;

use crate::abs::ShapeType;
use crate::brep_lib_build_curves3d::build_curves3d;
use crate::brep_tool::BRepTool;
use crate::chfi2d::ConstructionError;
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{explorer_shapes, map_shapes, shape_map_key, ShapeMapKey};

use super::ChFi2dBuilder;

impl Default for ChFi2dBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl ChFi2dBuilder {
    /// `ChFi2d_Builder()` (`cxx:66-68`): `status = ChFi2d_NotPlanar`.
    pub fn new() -> Self {
        Self {
            status: ConstructionError::NotPlanar,
            ref_face: None,
            new_face: None,
            fillets: Vec::new(),
            chamfers: Vec::new(),
            history: Vec::new(),
        }
    }

    /// `ChFi2d_Builder(const TopoDS_Face& F)` (`cxx:73-94`). The constructor
    /// clears nothing (the three containers are empty in a fresh object), so it
    /// shares the body of [`ChFi2dBuilder::init`].
    pub fn with_face(f: &Face) -> Self {
        let mut builder = Self::new();
        builder.init_face(f);
        builder
    }

    /// `Init(const TopoDS_Face& F)` (`cxx:99-122`).
    pub fn init(&mut self, f: &Face) {
        // `cxx:105-107`.
        self.fillets.clear();
        self.chamfers.clear();
        self.history.clear();
        self.init_face(f);
    }

    /// `cxx:81-94` / `cxx:110-121`: the body shared by the constructor and
    /// `Init(F)`.
    fn init_face(&mut self, f: &Face) {
        if surface_is_plane(f) {
            // `newFace = refFace = F; newFace.Orientation(FORWARD);`
            // `BRepLib::BuildCurves3d(newFace);` (`cxx:83-88`, `:112-115`).
            self.ref_face = Some(f.clone());
            let new_face = Face(f.0.oriented(crate::abs::Orientation::Forward));
            build_curves3d(&new_face.0);
            self.new_face = Some(new_face);
            self.status = ConstructionError::Ready;
        } else {
            self.status = ConstructionError::NotPlanar;
        }
    }

    /// `Init(const TopoDS_Face& RefFace, const TopoDS_Face& ModFace)`
    /// (`cxx:125-203`): adopt `ModFace` as the working face and classify every
    /// edge of `ModFace` that is absent from `RefFace` — a trimmed copy of a
    /// reference edge gives a history binding, otherwise the edge is a new
    /// fillet (circle) or chamfer (line).
    pub fn init_from_faces(&mut self, ref_face: &Face, mod_face: &Face) {
        // `cxx:130-136`.
        self.fillets.clear();
        self.chamfers.clear();
        self.history.clear();
        // `cxx:140-144`: `!BRep_Tool::Surface(RefFace)->IsKind(Geom_Plane)`.
        if !surface_is_plane(ref_face) {
            self.status = ConstructionError::NotPlanar;
            return;
        }
        // `cxx:146-149`.
        self.ref_face = Some(ref_face.clone());
        let new_face = Face(mod_face.0.oriented(crate::abs::Orientation::Forward));
        self.new_face = Some(new_face.clone());
        self.status = ConstructionError::Ready;

        // `cxx:154-155`: `TopExp::MapShapes(refFace, TopAbs_EDGE, refEdgesMap)`.
        let ref_edges: Vec<Edge> = map_shapes(&ref_face.0, &[ShapeType::Edge])
            .remove(&ShapeType::Edge)
            .unwrap_or_default()
            .into_iter()
            .filter_map(Edge::wrap)
            .collect();
        let ref_keys: HashSet<ShapeMapKey> = ref_edges.iter().map(|e| shape_map_key(&e.0)).collect();

        // `cxx:157-166`: the edges of `newFace` that are not in `refEdgesMap`.
        let mut new_edges: Vec<Edge> = Vec::new();
        for child in explorer_shapes(&new_face.0, ShapeType::Edge) {
            let Some(edge) = Edge::wrap(child) else { continue };
            if !ref_keys.contains(&shape_map_key(&edge.0)) {
                new_edges.push(edge);
            }
        }

        // `cxx:169-201`.
        for edge in &new_edges {
            if let Some(basis) = is_issued_from(edge, &ref_edges, &ref_keys) {
                self.bind_history(&basis.0, &edge.0);
                continue;
            }
            // `cxx:182-195`: the 3D curve decides — `Geom_Circle` is a fillet,
            // `Geom_Line` a chamfer, anything else is
            // `ChFi2d_InitialisationError`.
            let Some(curve) = BRepTool::edge_curve(edge) else {
                self.status = ConstructionError::InitialisationError;
                return;
            };
            if curve.gp_circ().is_some() {
                self.fillets.push(edge.0.clone());
            } else if curve.is_line() {
                self.chamfers.push(edge.0.clone());
            } else {
                self.status = ConstructionError::InitialisationError;
                return;
            }
        }
    }
}

/// `BRep_Tool::Surface(F, Loc)->IsKind(STANDARD_TYPE(Geom_Plane))`
/// (`ChFi2d_Builder.cxx:85`, `:113`, `:140`).
///
/// `IsKind` is a subclass test: a `Geom_RectangularTrimmedSurface` over a plane
/// answers `false` in OCCT. `gp_pln()` reports the exact `Geom_Plane` only, so
/// it is the same predicate.
fn surface_is_plane(f: &Face) -> bool {
    GeometryRegistry::global()
        .face_surface(&f.0)
        .map(|s| s.gp_pln().is_some())
        .unwrap_or(false)
}

/// `IsIssuedFrom(E, Map, BasisEdge)` (`cxx:208-239`): the reference edge whose
/// 3D curve is the **same handle** as `E`'s (`c1 == c2`, `cxx:228`) and whose
/// parameter range overlaps `E`'s (`cxx:229-231`).
fn is_issued_from(e: &Edge, candidates: &[Edge], map: &HashSet<ShapeMapKey>) -> Option<Edge> {
    // `BRep_Tool::Curve(E, loc1, f1, L1)` (`cxx:219`).
    let c1 = BRepTool::edge_curve(e)?;
    let (f1, l1) = BRepTool::edge_parameters(e);
    // `for (int i = 1; i <= Map.Extent(); i++)` (`cxx:222`).
    for candidate in candidates {
        if !map.contains(&shape_map_key(&candidate.0)) {
            continue;
        }
        // `BRep_Tool::Curve(currentEdge, loc2, f2, L2)` (`cxx:227`).
        let Some(c2) = BRepTool::edge_curve(candidate) else { continue };
        if !std::sync::Arc::ptr_eq(&c1, &c2) {
            continue;
        }
        let (f2, l2) = BRepTool::edge_parameters(candidate);
        // `cxx:229-231`.
        if ((f1 > f2 && f1 < l2) || (l1 > f2 && l1 < l2))
            || ((f1 > l2 && f1 < f2) || (l1 > l2 && l1 < f2))
        {
            // `BasisEdge = currentEdge; BasisEdge.Orientation(E.Orientation())`
            // (`cxx:233-235`).
            return Some(Edge(candidate.0.oriented(e.0.orientation())));
        }
    }
    None
}
