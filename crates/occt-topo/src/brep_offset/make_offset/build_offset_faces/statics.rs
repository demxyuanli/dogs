//! File-local helpers of `BRepOffset_MakeOffset_1.cxx` (`cxx:74-521`).
//!
//! These are the anonymous-namespace functions the file-local
//! [`super::BRepOffsetBuildOffsetFaces`] class drives. Each one keeps the
//! reference's branch structure; only the container types differ (the OCCT
//! collections are carried by their Rust counterparts, see each item).

use std::collections::HashMap;

use occt_core::gp::GpVec;

use crate::abs::{Orientation, ShapeType};
use crate::bop_builder2::BopBuilder;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::brep_offset::{Analyse, Dmvv};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use crate::pave_ff_se::is_micro_edge;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{
    edge_vertices, explorer_shapes, is_same, shape_map_key, shapes_of, vertex_position, ShapeMapKey,
};

/// `AddToContainer` overload set (`cxx:78-114`).
///
/// The reference has one overload per destination container:
///
/// | OCCT overload | `cxx` | answer |
/// |---------------|-------|--------|
/// | `NCollection_List` | `:79-82` | void - always "added" |
/// | `NCollection_Map` | `:84-87` | whether the shape was new |
/// | `NCollection_IndexedMap` | `:89-95` | whether the shape was new |
/// | `TopoDS_Shape&` (compound) | `:97-100` | void |
///
/// The Rust trait folds the four into one `bool`: the list and compound arms
/// answer `true`, so a caller that tests the answer only skips work for the
/// map arms, exactly like the reference's `!theMFence || Add(...)` condition.
/// The `NCollection_Map` and `NCollection_IndexedMap` overloads both map onto
/// [`ShapeIndexedMap`], whose insertion order is extra state that does not
/// change membership.
pub trait AddToContainer {
    /// `AddToContainer(theS, theOut)`: adds the shape, answering whether the
    /// destination treats it as new.
    fn add_to_container(&mut self, s: &TopoShape) -> bool;
}

/// `AddToContainer(theS, NCollection_List<TopoDS_Shape>&)` (`cxx:79-82`).
impl AddToContainer for Vec<TopoShape> {
    fn add_to_container(&mut self, s: &TopoShape) -> bool {
        self.push(s.clone());
        true
    }
}

/// `AddToContainer(theS, NCollection_Map&)` (`cxx:84-87`) and
/// `AddToContainer(theS, NCollection_IndexedMap&)` (`cxx:89-95`).
impl AddToContainer for ShapeIndexedMap {
    fn add_to_container(&mut self, s: &TopoShape) -> bool {
        let a_nb = self.extent();
        self.add(s) > a_nb
    }
}

/// `AddToContainer(theS, TopoDS_Shape&)` (`cxx:97-100`):
/// `BRep_Builder().Add(theSOut, theS)`.
impl AddToContainer for TopoShape {
    fn add_to_container(&mut self, s: &TopoShape) -> bool {
        TopoBuilder::new().add(self, s);
        true
    }
}

/// `AddToContainer(theKey, theValue, NCollection_DataMap<Shape, List>&)`
/// (`cxx:102-114`): binds an empty list on the first use of the key, then
/// appends.
pub fn add_to_container_of_map(
    the_map: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    the_key: &TopoShape,
    the_value: &TopoShape,
) {
    the_map
        .entry(shape_map_key(the_key))
        .or_default()
        .push(the_value.clone());
}

/// `static bool TakeModified(...)` (`cxx:124-164`), both overloads.
///
/// * `cxx:124-154` - the overload that takes `theMFence`;
/// * `cxx:156-164` - the overload that passes a null fence, carried here as
///   `fence: None`.
///
/// Answers `true` when `theS` had images in `theImages`: every image (or the
/// shape itself, when unbound) is put into `theContainer`, filtered by the
/// fence when one is given.
pub fn take_modified<C: AddToContainer>(
    the_s: &TopoShape,
    the_images: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    the_container: &mut C,
    mut fence: Option<&mut ShapeIndexedMap>,
) -> bool {
    match the_images.get(&shape_map_key(the_s)) {
        Some(p_ls_im) => {
            for a_s_im in p_ls_im {
                let a_add = match fence.as_deref_mut() {
                    Some(the_m_fence) => the_m_fence.add_to_container(a_s_im),
                    None => true,
                };
                if a_add {
                    the_container.add_to_container(a_s_im);
                }
            }
            true
        }
        None => {
            let a_add = match fence.as_deref_mut() {
                Some(the_m_fence) => the_m_fence.add_to_container(the_s),
                None => true,
            };
            if a_add {
                the_container.add_to_container(the_s);
            }
            false
        }
    }
}

/// The `IsEmpty()` that `hasData` tests on its container.
pub trait HasData {
    /// `Container::IsEmpty()`.
    fn is_empty_container(&self) -> bool;
}

impl<T> HasData for Vec<T> {
    fn is_empty_container(&self) -> bool {
        self.is_empty()
    }
}

impl<K, V> HasData for HashMap<K, V> {
    fn is_empty_container(&self) -> bool {
        self.is_empty()
    }
}

impl HasData for ShapeIndexedMap {
    fn is_empty_container(&self) -> bool {
        self.extent() == 0
    }
}

impl HasData for Dmvv {
    fn is_empty_container(&self) -> bool {
        self.extent() == 0
    }
}

/// `static bool hasData(const Container* theData)` (`cxx:170-175`).
///
/// The reference also tests the pointer for null. Every destination is an
/// owned container in this port (see the module note of
/// [`super::BRepOffsetBuildOffsetFaces`]), so an unset destination is empty and
/// the null arm cannot be reached.
pub fn has_data<C: HasData>(the_data: &C) -> bool {
    !the_data.is_empty_container()
}

/// `static void AppendToList(List&, const TopoDS_Shape&)` (`cxx:180-191`):
/// appends only when no `IsSame` element is present.
pub fn append_to_list(the_list: &mut Vec<TopoShape>, the_shape: &TopoShape) {
    for a_s in the_list.iter() {
        if is_same(a_s, the_shape) {
            return;
        }
    }
    the_list.push(the_shape.clone());
}

/// `static bool ProcessMicroEdge(const TopoDS_Edge&, const IntTools_Context&)`
/// (`cxx:197-216`).
///
/// `BOPTools_AlgoTools::IsMicroEdge(theEdge, theCtx)` maps onto
/// [`crate::pave_ff_se::is_micro_edge`], which builds its own
/// `IntTools_ShrunkRange`; the reference's context argument therefore has no
/// counterpart and the signature drops it. The line-edge tolerance update keeps
/// `BRep_Builder::UpdateVertex(V, Tol)` semantics - `max(current, Tol)`
/// (`BRep_TVertex.hxx:49-52`).
pub fn process_micro_edge(the_edge: &Edge) -> bool {
    // `cxx:201-205`: both vertices must be present.
    let (v1, v2) = edge_vertices(the_edge);
    let (Some(a_v1), Some(a_v2)) = (v1, v2) else {
        return false;
    };

    // `cxx:208`: `BOPTools_AlgoTools::IsMicroEdge` with `theCheckSplittable`
    // defaulting to `true`.
    let b_micro = is_micro_edge(the_edge, true);
    if b_micro {
        // `cxx:209`: `BRepAdaptor_Curve(theEdge).GetType() == GeomAbs_Line`.
        let is_line = match BRepTool::edge_curve(the_edge) {
            Some(a_c) => curve_kind(a_c.as_ref()) == GeomAbsCurveType::Line,
            None => false,
        };
        if is_line {
            // `cxx:210-213`: `aLen = Distance(PntV1, PntV2)`, then
            // `UpdateVertex(V, aLen / 2)` on both vertices.
            let a_len = vertex_position(&a_v1).distance(&vertex_position(&a_v2));
            let a_tol = a_len / 2.;
            a_v1.set_tolerance(BRepTool::vertex_tolerance(&a_v1).max(a_tol));
            a_v2.set_tolerance(BRepTool::vertex_tolerance(&a_v2).max(a_tol));
        }
    }

    b_micro
}

/// `static void UpdateOrigins(const List&, DataMap<Shape, List>&,
/// BOPAlgo_Builder&)` (`cxx:220-261`).
///
/// `theGF.Modified(aS)` maps onto `BopHistory::modified` through
/// [`BopBuilder::history`], the port's images/origins carrier.
pub fn update_origins(
    the_la: &[TopoShape],
    the_origins: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    the_gf: &BopBuilder,
) {
    for a_s in the_la {
        // `cxx:230-234`: skip the shapes the general fuse did not modify.
        let a_ls_im = match the_gf.history().modified(a_s) {
            Some(a_list) if !a_list.is_empty() => a_list,
            _ => continue,
        };

        // `cxx:237-242`: the stored list when the source is already a key,
        // otherwise a one-element list seeded with the source itself.
        let p_ls: Vec<TopoShape> = match the_origins.get(&shape_map_key(a_s)) {
            Some(a_list) => a_list.clone(),
            None => vec![a_s.clone()],
        };

        for a_s_im in a_ls_im {
            match the_origins.get_mut(&shape_map_key(a_s_im)) {
                // `cxx:248-256`: merge the two lists.
                Some(p_ls_or) => {
                    for a_value in &p_ls {
                        append_to_list(p_ls_or, a_value);
                    }
                }
                // `cxx:257-260`: bind a copy of the source list.
                None => {
                    the_origins.insert(shape_map_key(a_s_im), p_ls.clone());
                }
            }
        }
    }
}

/// `static void UpdateImages(const List&, DataMap<Shape, List>&,
/// BOPAlgo_Builder&, Map&)` (`cxx:265-305`).
///
/// `theGF.Images()` maps onto `BopHistory::images` through
/// [`BopBuilder::history`].
pub fn update_images(
    the_la: &[TopoShape],
    the_images: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    the_gf: &BopBuilder,
    the_modified: &mut ShapeIndexedMap,
) {
    for a_s in the_la {
        let a_key = shape_map_key(a_s);
        // `cxx:276-287`: no list under the key yet - bind the modified images.
        let Some(p_ls_im) = the_images.get(&a_key).cloned() else {
            let a_ls_im: &[TopoShape] = match the_gf.history().modified(a_s) {
                Some(a_list) => a_list,
                None => &[],
            };
            if !a_ls_im.is_empty() {
                the_images.insert(a_key, a_ls_im.to_vec());
                the_modified.add_to_container(a_s);
            }
            continue;
        };

        // `cxx:288-298`: re-run TakeModified over the existing images.
        let mut a_m_fence = ShapeIndexedMap::new();
        let mut a_ls_im_new: Vec<TopoShape> = Vec::new();
        let mut b_modified = false;
        for a_s_im in &p_ls_im {
            b_modified |= take_modified(
                a_s_im,
                the_gf.history().images(),
                &mut a_ls_im_new,
                Some(&mut a_m_fence),
            );
        }

        // `cxx:300-304`.
        if b_modified {
            the_images.insert(a_key, a_ls_im_new);
            the_modified.add_to_container(a_s);
        }
    }
}

/// `static void UpdateImages(const List& theLA, DataMap& theImages,
/// const BOPAlgo_Builder& theGF, IndexedMap& theModified)`
/// (`cxx:265-305`) over the `IndexedDataMap` form of the images:
/// `BRepOffset_BuildOffsetFaces::myOEImages` (`hxx:230`). The binding arm of
/// `cxx:276-287` uses `Bind`, which replaces the key's list, so the key order
/// of the `IndexedDataMap` is preserved.
pub fn update_images_dmvv(
    the_la: &[TopoShape],
    the_images: &mut Dmvv,
    the_gf: &BopBuilder,
    the_modified: &mut ShapeIndexedMap,
) {
    for a_s in the_la {
        let a_key = shape_map_key(a_s);
        // `cxx:276-287`: no list under the key yet - bind the modified images.
        let Some(p_ls_im) = the_images.seek(a_s).map(|l| l.to_vec()) else {
            let a_ls_im: &[TopoShape] = match the_gf.history().modified(a_s) {
                Some(a_list) => a_list,
                None => &[],
            };
            if !a_ls_im.is_empty() {
                the_images.bind_list(a_s, a_ls_im);
                the_modified.add_to_container(a_s);
            }
            continue;
        };

        // `cxx:288-298`: re-run TakeModified over the existing images.
        let mut a_m_fence = ShapeIndexedMap::new();
        let mut a_ls_im_new: Vec<TopoShape> = Vec::new();
        let mut b_modified = false;
        for a_s_im in &p_ls_im {
            b_modified |= take_modified(
                a_s_im,
                the_gf.history().images(),
                &mut a_ls_im_new,
                Some(&mut a_m_fence),
            );
        }

        // `cxx:300-304`.
        if b_modified {
            the_images.bind_list(a_s, &a_ls_im_new);
            the_modified.add_to_container(a_s);
        }
        let _ = a_key;
    }
}

/// `static bool checkConnectionsOfFace(const TopoDS_Shape&, const List&)`
/// (`cxx:1301-1341`).
///
/// Counts the vertices the face shares with the other faces of the list and
/// answers `true` as soon as a second one is found. The reference builds the
/// two vertex sets with `TopExp::MapShapes` (distinct `IsSame` keys), which
/// [`shapes_of`] already answers.
pub fn check_connections_of_face(the_face: &TopoShape, the_lf: &[TopoShape]) -> bool {
    // `cxx:1305-1315`.
    let mut a_shape_vert = ShapeIndexedMap::new();
    for a_shape in the_lf {
        if is_same(a_shape, the_face) {
            continue;
        }
        for a_v in shapes_of(a_shape, ShapeType::Vertex) {
            a_shape_vert.add(&a_v);
        }
    }

    // `cxx:1316-1334`.
    let mut a_nb_connections = 0;
    for a_v in shapes_of(the_face, ShapeType::Vertex) {
        if a_shape_vert.contains(&a_v) {
            a_nb_connections += 1;
        }
        if a_nb_connections > 1 {
            return true;
        }
    }
    false
}

/// `static void FindCommonParts(const List&, const List&, List&,
/// const TopAbs_ShapeEnum = TopAbs_EDGE)` (`cxx:312-346`).
pub fn find_common_parts(
    the_ls1: &[TopoShape],
    the_ls2: &[TopoShape],
    the_lsc: &mut Vec<TopoShape>,
    the_type: ShapeType,
) {
    // `cxx:318-327`: map the sub-shapes of the first list.
    let mut a_ms1 = ShapeIndexedMap::new();
    for a_s in the_ls1 {
        for a_sub in explorer_shapes(a_s, the_type) {
            a_ms1.add(&a_sub);
        }
    }
    if a_ms1.extent() == 0 {
        return;
    }

    // `cxx:330-345`: append the matching sub-shapes of the second list once.
    let mut a_m_fence = ShapeIndexedMap::new();
    for a_s in the_ls2 {
        for a_st in explorer_shapes(a_s, the_type) {
            if a_ms1.contains(&a_st) && a_m_fence.add_to_container(&a_st) {
                the_lsc.push(a_st);
            }
        }
    }
}

/// `static int NbPoints(const TopoDS_Edge&)` (`cxx:351-362`): the number of
/// sample points used for the average direction of an edge.
pub fn nb_points(the_edge: &Edge) -> usize {
    match BRepTool::edge_curve(the_edge) {
        // `cxx:356-357`.
        Some(a_c) if curve_kind(a_c.as_ref()) == GeomAbsCurveType::Line => 1,
        // `cxx:358-359`.
        _ => 11,
    }
}

/// `static bool FindShape(const TopoDS_Shape&, const TopoDS_Shape&,
/// const BRepOffset_Analyse*, TopoDS_Shape&)` (`cxx:367-405`).
///
/// `theAnalyse->Descendants(theSWhere)` maps onto
/// [`Analyse::descendants`] with `update = false` (the reference's default
/// argument).
pub fn find_shape(
    the_s_what: &TopoShape,
    the_s_where: &TopoShape,
    the_analyse: Option<&Analyse>,
    the_res: &mut TopoShape,
) -> bool {
    let mut b_found = false;
    let a_type = the_s_what.shape_type();

    // `cxx:372-385`.
    for a_s in explorer_shapes(the_s_where, a_type) {
        if is_same(&a_s, the_s_what) {
            *the_res = a_s;
            b_found = true;
            break;
        }
    }

    // `cxx:387-402`.
    if !b_found {
        if let Some(a_analyse) = the_analyse {
            if let Some(a_ld) = a_analyse.descendants(the_s_where, false) {
                for a_s in a_ld {
                    if is_same(&a_s, the_s_what) {
                        *the_res = a_s;
                        b_found = true;
                        break;
                    }
                }
            }
        }
    }

    b_found
}

/// `static gp_Vec GetAverageTangent(const TopoDS_Shape&, const int)`
/// (`cxx:491-521`).
///
/// `BRep_Tool::Curve(aE, aT1, aT2)` maps onto
/// [`BRepTool::edge_curve`] plus [`BRepTool::edge_parameters`].
pub fn get_average_tangent(the_s: &TopoShape, the_nb_p: usize) -> GpVec {
    let mut a_va = GpVec::zero();

    // `cxx:494-495`.
    for a_shape in explorer_shapes(the_s, ShapeType::Edge) {
        let Some(a_e) = Edge::wrap(a_shape.clone()) else {
            continue;
        };
        // `cxx:499-501`.
        let Some(a_c) = BRepTool::edge_curve(&a_e) else {
            continue;
        };
        let (a_t1, a_t2) = BRepTool::edge_parameters(&a_e);

        // `cxx:503-511`.
        let mut a_v_sum = GpVec::zero();
        let a_dt = (a_t2 - a_t1) / the_nb_p as f64;
        let mut a_t = a_t1;
        while a_t <= a_t2 {
            let (_, a_v) = a_c.d1(a_t);
            a_v_sum = a_v_sum.added(&a_v.normalized());
            a_t += a_dt;
        }

        // `cxx:513-516`.
        if a_shape.orientation() == Orientation::Reversed {
            a_v_sum.reverse();
        }

        // `cxx:518`.
        a_va = a_va.added(&a_v_sum);
    }

    a_va
}
