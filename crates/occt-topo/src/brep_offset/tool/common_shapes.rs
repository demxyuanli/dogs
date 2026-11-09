//! `BRepOffset_Tool` helpers that map shapes of one face / vertex / edge onto
//! another: `FindCommonShapes` (`BRepOffset_Tool.cxx:572-618`),
//! `MapVertexEdges` (`:3685-3716`) and `BuildNeighbour` (`:3720-3789`).

use std::collections::{HashMap, HashSet};

use crate::abs::{Orientation, ShapeType};
use crate::brep_offset::tool::BRepOffsetTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools::{
    edge_vertices, explorer_shapes, is_same_occurrence, shape_map_key, ShapeMapKey,
};

impl BRepOffsetTool {
    /// `BRepOffset_Tool::FindCommonShapes(theF1, theF2, theLE, theLV)`
    /// (`cxx:572-580`): the shared edges and shared vertices of two faces.
    ///
    /// Like OCCT, the two output lists are *appended* to (the type overload
    /// below is the one that clears), and the answer is whether either list
    /// ended up non-empty.
    pub fn find_common_shapes(
        f1: &Face,
        f2: &Face,
        the_le: &mut Vec<TopoShape>,
        the_lv: &mut Vec<TopoShape>,
    ) -> bool {
        let b_found_edges =
            Self::find_common_shapes_of_type(&f1.0, &f2.0, ShapeType::Edge, the_le);
        let b_found_vertices =
            Self::find_common_shapes_of_type(&f1.0, &f2.0, ShapeType::Vertex, the_lv);
        b_found_edges || b_found_vertices
    }

    /// `BRepOffset_Tool::FindCommonShapes(theS1, theS2, theType, theLSC)`
    /// (`cxx:584-618`): the shapes of type `the_type` present in both `s1` and
    /// `s2`, in `s2`'s exploration order and without repetitions.
    ///
    /// OCCT's fences are `NCollection_Map<TopoDS_Shape, ShapeMapHasher>`, i.e.
    /// `IsSame` equality - [`ShapeMapKey`] is exactly that key, so the
    /// orientation of a shared sub-shape plays no part.
    pub fn find_common_shapes_of_type(
        s1: &TopoShape,
        s2: &TopoShape,
        the_type: ShapeType,
        the_lsc: &mut Vec<TopoShape>,
    ) -> bool {
        the_lsc.clear();
        // `cxx:592-598`.
        let a_ms: HashSet<ShapeMapKey> =
            explorer_shapes(s1, the_type).iter().map(shape_map_key).collect();
        if a_ms.is_empty() {
            return false;
        }
        // `cxx:604-614`: `aMFence.Add(aS2)` both deduplicates and reports the
        // first sighting.
        let mut a_m_fence: HashSet<ShapeMapKey> = HashSet::new();
        for a_s2 in explorer_shapes(s2, the_type) {
            let key = shape_map_key(&a_s2);
            if a_ms.contains(&key) && a_m_fence.insert(key) {
                the_lsc.push(a_s2);
            }
        }
        // `cxx:617`.
        !the_lsc.is_empty()
    }

    /// `BRepOffset_Tool::MapVertexEdges(S, MVE)` (`cxx:3685-3716`): for every
    /// vertex of `S`, the edges incident to it, each stored once.
    ///
    /// The map is keyed by [`ShapeMapKey`] (OCCT's
    /// `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
    /// ShapeMapHasher>`, `IsSame` semantics). `TopExp_Explorer(S.Oriented(
    /// FORWARD), EDGE)` + the `DejaVu` fence is reproduced by the `seen` set:
    /// an edge occurrence met twice is appended only on its first sighting.
    pub fn map_vertex_edges(s: &TopoShape) -> HashMap<ShapeMapKey, Vec<TopoShape>> {
        let mut mev: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        let mut deja_vu: HashSet<ShapeMapKey> = HashSet::new();
        for shape in explorer_shapes(&s.oriented(Orientation::Forward), ShapeType::Edge) {
            let e = Edge(shape);
            if !deja_vu.insert(shape_map_key(&e.0)) {
                continue;
            }
            let (v1, v2) = edge_vertices(&e);
            // `cxx:3699-3705`: `Bind(V1, empty)` then `Append(E)`.
            let (Some(v1), Some(v2)) = (v1, v2) else {
                // Port guard: OCCT would call `MEV(V1)` on a null vertex.
                continue;
            };
            mev.entry(shape_map_key(&v1.0)).or_default().push(e.0.clone());
            if !is_same_occurrence(&v1.0, &v2.0) {
                mev.entry(shape_map_key(&v2.0)).or_default().push(e.0.clone());
            }
        }
        mev
    }

    /// `BRepOffset_Tool::BuildNeighbour(W, F, NOnV1, NOnV2)` (`cxx:3720-3789`):
    /// for every edge of the wire, its neighbour through the first (`NOnV1`)
    /// and last (`NOnV2`) vertex.
    ///
    /// Both maps are OCCT `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape>`
    /// with `ShapeMapHasher` keys, so a later `Bind` of the same edge
    /// overwrites the earlier one - `HashMap::insert` has the same effect.
    ///
    /// `None` stands for an empty wire: `wexp.Current()` at `cxx:3734` faults
    /// on OCCT's empty explorer, and the `while (wexp.More())` body would then
    /// read a null edge.
    pub fn build_neighbour(
        w: &Wire,
        f: &Face,
        n_on_v1: &mut HashMap<ShapeMapKey, TopoShape>,
        n_on_v2: &mut HashMap<ShapeMapKey, TopoShape>,
    ) -> Option<()> {
        // `cxx:3730-3733`: both shapes oriented FORWARD.
        let a_local_face = f.0.oriented(Orientation::Forward);
        let a_local_wire = w.0.oriented(Orientation::Forward);
        let a_face = Face::wrap(a_local_face)?;
        let a_wire = Wire::wrap(a_local_wire)?;
        let mut wexp = WireExplorer::of_wire_on_face(&a_wire, &a_face);

        // `cxx:3734-3738`.
        let first_e = wexp.current_edge()?;
        let mut cur_e = first_e.clone();
        let mut prec_e = first_e.clone();
        let (v1, v2) = edge_vertices(&cur_e);
        let (mut v1, mut v2) = (v1?, v2?);
        let fv1 = v1.clone();
        let fv2 = v2.clone();
        let mut vp1 = v1.clone();
        let mut vp2 = v2.clone();
        wexp.next();

        // `cxx:3740-3768`.
        while wexp.more() {
            cur_e = wexp.current_edge()?;
            let (nv1, nv2) = edge_vertices(&cur_e);
            v1 = nv1?;
            v2 = nv2?;
            if is_same_occurrence(&v1.0, &vp1.0) {
                n_on_v1.insert(shape_map_key(&prec_e.0), cur_e.0.clone());
                n_on_v1.insert(shape_map_key(&cur_e.0), prec_e.0.clone());
            }
            if is_same_occurrence(&v1.0, &vp2.0) {
                n_on_v2.insert(shape_map_key(&prec_e.0), cur_e.0.clone());
                n_on_v1.insert(shape_map_key(&cur_e.0), prec_e.0.clone());
            }
            if is_same_occurrence(&v2.0, &vp1.0) {
                n_on_v1.insert(shape_map_key(&prec_e.0), cur_e.0.clone());
                n_on_v2.insert(shape_map_key(&cur_e.0), prec_e.0.clone());
            }
            if is_same_occurrence(&v2.0, &vp2.0) {
                n_on_v2.insert(shape_map_key(&prec_e.0), cur_e.0.clone());
                n_on_v2.insert(shape_map_key(&cur_e.0), prec_e.0.clone());
            }
            prec_e = cur_e.clone();
            vp1 = v1.clone();
            vp2 = v2.clone();
            wexp.next();
        }

        // `cxx:3769-3788`: the closing edge pair.
        if is_same_occurrence(&v1.0, &fv1.0) {
            n_on_v1.insert(shape_map_key(&first_e.0), cur_e.0.clone());
            n_on_v1.insert(shape_map_key(&cur_e.0), first_e.0.clone());
        }
        if is_same_occurrence(&v1.0, &fv2.0) {
            n_on_v2.insert(shape_map_key(&first_e.0), cur_e.0.clone());
            n_on_v1.insert(shape_map_key(&cur_e.0), first_e.0.clone());
        }
        if is_same_occurrence(&v2.0, &fv1.0) {
            n_on_v1.insert(shape_map_key(&first_e.0), cur_e.0.clone());
            n_on_v2.insert(shape_map_key(&cur_e.0), first_e.0.clone());
        }
        if is_same_occurrence(&v2.0, &fv2.0) {
            n_on_v2.insert(shape_map_key(&first_e.0), cur_e.0.clone());
            n_on_v2.insert(shape_map_key(&cur_e.0), first_e.0.clone());
        }
        Some(())
    }
}
