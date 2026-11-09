//! `HLRBRep_HLRToShape` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_HLRToShape.hxx:60-175`,
//! `.cxx:33-305`, `.lxx:19-183`): the filter that turns the exact HLR results
//! held by an [`HlrBRepAlgo`] into 2D compounds.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `HLRBRep_HLRToShape(A)` | `cxx:33-35` | [`HlrBRepHlrToShape::new`] |
//! | `InternalCompound(typ, visible, S, In3d)` | `cxx:43-150` | the private `internal_compound` |
//! | `DrawFace(...)` | `cxx:154-226` | the private `draw_face` |
//! | `DrawEdge(...)` | `cxx:232-305` | the private `draw_edge` |
//! | the 13 entry points | `lxx:21-183` | the `*_compound` / `*_compound_of` pairs |
//!
//! # Translation notes
//!
//! * `myAlgo` is an `occ::handle<HLRBRep_Algo>`; the port owns an
//!   [`HlrBRepAlgo`] value. Every extraction call reads and rewrites the same
//!   `myDS` through it, which is the sharing the handle gives in OCCT.
//! * `DS->Projector().Scaled(true)` (`cxx:51`) and `Scaled(false)`
//!   (`cxx:143`) write into the data structure's projector. The port's
//!   `HlrBRepData::my_proj` is an `Arc` also handed to every
//!   `HLRBRep_Curve` (`hlrbrep_data/update.rs:42`), so `Arc::make_mut` clones
//!   it when those curves are alive: the toggle reaches the data structure's
//!   own `Projector()` but **not** the curve evaluations. OCCT keeps one
//!   projector object and the curves see the toggle too. This only shows up
//!   through `HLRBRep_Curve::Value` in the fallback arm of
//!   [`crate::hlrbrep::make_edge`] and only when the projector carries a scale
//!   factor or a translation (see the batch report).
//! * `NCollection_IndexedMap::FindIndex` (`cxx:96`, `:108`) is an `IsSame`
//!   lookup; the port builds one `shape_map_key` map per call, first
//!   occurrence wins, and the entry order is the `EdgeMap()` / `FaceMap()`
//!   insertion order.
//! * `TopExp_Explorer(S, TopAbs_EDGE, TopAbs_FACE)` (`cxx:106`) is
//!   [`crate::topo_tools::explorer_shapes_avoiding`].
//! * `bool& added` (`cxx:47`) is a `&mut bool`, `TopoDS_Shape& Result`
//!   (`cxx:48`) a `&mut TopoShape`; a null result is `None`.

use std::collections::HashMap;
use std::sync::Arc;

use crate::abs::ShapeType;
use crate::builder::TopoBuilder;
use crate::hlr_algo::HlrAlgoEdgeIterator;
use crate::hlrbrep;
use crate::hlrbrep_algo::HlrBRepAlgo;
use crate::hlrbrep_data::HlrBRepData;
use crate::hlrbrep_edge_data::HlrBRepEdgeData;
use crate::hlrbrep_face_iterator::HlrBRepFaceIterator;
use crate::hlrbrep_type_of_resulting_edge::HlrBRepTypeOfResultingEdge;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, explorer_shapes_avoiding, shape_map_key, ShapeMapKey};

/// `NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>::FindIndex`
/// of one stored map: 1-based position, first occurrence of a key wins.
fn index_of(map: &[TopoShape]) -> HashMap<ShapeMapKey, i32> {
    let mut out = HashMap::with_capacity(map.len());
    for (i, s) in map.iter().enumerate() {
        out.entry(shape_map_key(s)).or_insert(i as i32 + 1);
    }
    out
}

/// `HLRBRep_HLRToShape`.
pub struct HlrBRepHlrToShape {
    /// `occ::handle<HLRBRep_Algo> myAlgo` (`hxx:173`).
    my_algo: HlrBRepAlgo,
}

impl HlrBRepHlrToShape {
    /// `HLRBRep_HLRToShape(const occ::handle<HLRBRep_Algo>& A)` (`cxx:33-35`).
    pub fn new(a: HlrBRepAlgo) -> Self {
        Self { my_algo: a }
    }

    /// The algorithm the filter reads.
    pub fn algo(&self) -> &HlrBRepAlgo {
        &self.my_algo
    }

    /// `InternalCompound(typ, visible, S, In3d)` (`cxx:43-150`).
    ///
    /// `None` is the null `TopoDS_Shape` of `cxx:145-148`: either the data
    /// structure is null (`cxx:49`) or nothing was added (`cxx:145`).
    fn internal_compound(
        &mut self,
        typ: i32,
        visible: bool,
        s: Option<&TopoShape>,
        in3d: bool,
    ) -> Option<TopoShape> {
        let mut added = false;

        // `cxx:48-55`: the default range covers the whole data structure.
        // `DS->NbEdges()` / `NbFaces()` are read before the optional shape
        // narrows them (`cxx:52-54`), so the null check of `cxx:49` comes first.
        let (nb_edges, nb_faces) = {
            let ds = self.my_algo.base().data_structure()?;
            (ds.nb_edges(), ds.nb_faces())
        };
        let mut e1 = 1;
        let mut e2 = nb_edges;
        let mut f1 = 1;
        let mut f2 = nb_faces;

        // `cxx:57-70`.
        let mut explor = false;
        if let Some(s) = s {
            let index = self.my_algo.index(s);
            if index == 0 {
                // `cxx:63-65`: the shape is not loaded, walk its own faces.
                explor = true;
            } else if let Some(sb) = self.my_algo.base().shape_bounds(index) {
                // `cxx:68`: `myAlgo->ShapeBounds(index).Bounds(v1, v2, e1, e2, f1, f2)`.
                let (_v1, _v2, a, b, c, d) = sb.bounds();
                e1 = a;
                e2 = b;
                f1 = c;
                f2 = d;
            }
        }

        let ds = Arc::make_mut(self.my_algo.base_mut().data_structure_mut()?);

        // `cxx:51`: `DS->Projector().Scaled(true)`.
        if let Ok(()) = Arc::make_mut(ds.projector_mut()).scaled(true) {}

        // `cxx:72-87`: the selection reset. `ed = &ChangeValue(e1 - 1); ed++`
        // is the `e1..=e2` walk over the 1-based slots.
        for ie in e1..=e2 {
            let ed = ds.e_data_mut(ie);
            if ed.selected() && !ed.vertical() {
                ed.set_used(false);
                ed.set_hide_count(0);
            } else {
                ed.set_used(true);
            }
        }

        // `cxx:71-72`: `BRep_Builder B; B.MakeCompound(TopoDS::Compound(Result))`.
        let builder = TopoBuilder::new();
        let mut result = builder.make_shape(ShapeType::Compound);

        if explor {
            // `cxx:88-119`.
            let s = s.expect("explor implies a shape");
            let edges = index_of(ds.edge_map());
            let faces = index_of(ds.face_map());

            // `cxx:94-101`: every FACE occurrence of S.
            for f in explorer_shapes(s, ShapeType::Face) {
                if let Some(&iface) = faces.get(&shape_map_key(&f)) {
                    Self::draw_face_impl(visible, typ, iface, ds, &mut result, &mut added, in3d);
                }
            }
            if typ >= 3 {
                // `cxx:106-117`: every EDGE of S that is not inside one of its
                // faces.
                for e in explorer_shapes_avoiding(s, ShapeType::Edge, ShapeType::Face) {
                    let Some(&ie) = edges.get(&shape_map_key(&e)) else {
                        continue;
                    };
                    let ed = ds.e_data_mut(ie);
                    if !ed.used() {
                        Self::draw_edge_impl(visible, false, typ, ed, &mut result, &mut added, in3d);
                        ed.set_used(true);
                    }
                }
            }
        } else {
            // `cxx:120-141`.
            for iface in f1..=f2 {
                Self::draw_face_impl(visible, typ, iface, ds, &mut result, &mut added, in3d);
            }
            if typ >= 3 {
                for ie in e1..=e2 {
                    let ed = ds.e_data_mut(ie);
                    if !ed.used() {
                        Self::draw_edge_impl(visible, false, typ, ed, &mut result, &mut added, in3d);
                        ed.set_used(true);
                    }
                }
            }
        }

        // `cxx:143`: `DS->Projector().Scaled(false)`.
        if let Ok(()) = Arc::make_mut(ds.projector_mut()).scaled(false) {}

        // `cxx:145-148`.
        if !added {
            return None;
        }
        Some(result)
    }

    /// `DrawFace(visible, typ, iface, DS, Result, added, In3d)`
    /// (`cxx:154-226`): walk the face's edge block and draw the edges the
    /// `typ` filter keeps.
    fn draw_face_impl(
        visible: bool,
        typ: i32,
        iface: i32,
        ds: &mut HlrBRepData,
        result: &mut TopoShape,
        added: &mut bool,
        in3d: bool,
    ) {
        // `cxx:162-163`: `HLRBRep_FaceIterator Itf; Itf.InitEdge(FDataArray(iface))`.
        let mut itf = HlrBRepFaceIterator::new();
        itf.init_edge(ds.f_data(iface));

        while itf.more_edge() {
            // `cxx:166-168`.
            let ie = itf.edge();
            let edf = ds.e_data_mut(ie);
            if !edf.used() {
                // `cxx:170-198`: the `todraw` arms.
                let todraw = if typ == 1 {
                    itf.iso_line()
                } else if typ == 2 {
                    // outlines
                    if in3d {
                        itf.internal() || itf.out_line()
                    } else {
                        itf.internal()
                    }
                } else if typ == 3 {
                    edf.rg1_line() && !edf.rgn_line() && !itf.out_line()
                } else if typ == 4 {
                    edf.rgn_line() && !itf.out_line()
                } else {
                    !itf.iso_line() && !itf.internal() && (!edf.rg1_line() || itf.out_line())
                };

                if todraw {
                    // `cxx:199-203`.
                    Self::draw_edge_impl(visible, true, typ, edf, result, added, in3d);
                    edf.set_used(true);
                } else if (typ > 4 || typ == 2) && edf.rg1_line() && !itf.out_line() {
                    // `cxx:206-219`: a sharp/outline miss on an Rg1 edge is
                    // retried against the next face.
                    let hc = edf.hide_count();
                    if hc > 0 {
                        edf.set_used(true);
                    } else {
                        edf.set_hide_count(hc + 1);
                    }
                } else {
                    // `cxx:221-223`.
                    edf.set_used(true);
                }
            }
            itf.next_edge();
        }
    }

    /// `DrawEdge(visible, inFace, typ, ed, Result, added, In3d)`
    /// (`cxx:232-305`): the visible (or hidden) parts of one edge as 2D (or 3D)
    /// edges added to `result`.
    fn draw_edge_impl(
        visible: bool,
        in_face: bool,
        typ: i32,
        ed: &HlrBRepEdgeData,
        result: &mut TopoShape,
        added: &mut bool,
        in3d: bool,
    ) {
        // `cxx:239-254`: the `todraw` arms.
        let todraw = if in_face {
            true
        } else if typ == 3 {
            ed.rg1_line() && !ed.rgn_line()
        } else if typ == 4 {
            ed.rgn_line()
        } else {
            !ed.rg1_line()
        };
        if !todraw {
            return;
        }

        // `cxx:262-264`.
        let builder = TopoBuilder::new();
        let mut it = HlrAlgoEdgeIterator::new();

        if visible {
            // `cxx:266-283`.
            it.init_visible(ed.status());
            while it.more_visible() {
                let (sta, _tol_sta, end, _tol_end) = it.visible();
                if let Some(e) = make_result_edge(ed, sta, end, in3d) {
                    // `cxx:279-281`: `B.Add(Result, E); added = true;`.
                    builder.add(result, &e.0);
                    *added = true;
                }
                it.next_visible();
            }
        } else {
            // `cxx:286-302`.
            it.init_hidden(ed.status());
            while it.more_hidden() {
                let (sta, _tol_sta, end, _tol_end) = it.hidden();
                if let Some(e) = make_result_edge(ed, sta, end, in3d) {
                    builder.add(result, &e.0);
                    *added = true;
                }
                it.next_hidden();
            }
        }
    }

    // ---- the entry points of `lxx:21-183` ----

    /// `VCompound()` (`lxx:21-24`): visible sharp edges.
    pub fn v_compound(&mut self) -> Option<TopoShape> {
        self.internal_compound(5, true, None, false)
    }

    /// `VCompound(S)` (`lxx:28-31`).
    pub fn v_compound_of(&mut self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(5, true, Some(s), false)
    }

    /// `Rg1LineVCompound()` (`lxx:35-38`).
    pub fn rg1_line_v_compound(&mut self) -> Option<TopoShape> {
        self.internal_compound(3, true, None, false)
    }

    /// `Rg1LineVCompound(S)` (`lxx:42-45`).
    pub fn rg1_line_v_compound_of(&mut self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(3, true, Some(s), false)
    }

    /// `RgNLineVCompound()` (`lxx:49-52`).
    pub fn rgn_line_v_compound(&mut self) -> Option<TopoShape> {
        self.internal_compound(4, true, None, false)
    }

    /// `RgNLineVCompound(S)` (`lxx:56-59`).
    pub fn rgn_line_v_compound_of(&mut self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(4, true, Some(s), false)
    }

    /// `OutLineVCompound()` (`lxx:63-66`).
    pub fn out_line_v_compound(&mut self) -> Option<TopoShape> {
        self.internal_compound(2, true, None, false)
    }

    /// `OutLineVCompound3d()` (`lxx:70-73`).
    pub fn out_line_v_compound_3d(&mut self) -> Option<TopoShape> {
        self.internal_compound(2, true, None, true)
    }

    /// `OutLineVCompound(S)` (`lxx:77-80`).
    pub fn out_line_v_compound_of(&mut self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(2, true, Some(s), false)
    }

    /// `IsoLineVCompound()` (`lxx:84-87`).
    pub fn iso_line_v_compound(&mut self) -> Option<TopoShape> {
        self.internal_compound(1, true, None, false)
    }

    /// `IsoLineVCompound(S)` (`lxx:91-94`).
    pub fn iso_line_v_compound_of(&mut self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(1, true, Some(s), false)
    }

    /// `HCompound()` (`lxx:98-101`): hidden sharp edges.
    pub fn h_compound(&mut self) -> Option<TopoShape> {
        self.internal_compound(5, false, None, false)
    }

    /// `HCompound(S)` (`lxx:105-108`).
    pub fn h_compound_of(&mut self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(5, false, Some(s), false)
    }

    /// `Rg1LineHCompound()` (`lxx:112-115`).
    pub fn rg1_line_h_compound(&mut self) -> Option<TopoShape> {
        self.internal_compound(3, false, None, false)
    }

    /// `Rg1LineHCompound(S)` (`lxx:119-122`).
    pub fn rg1_line_h_compound_of(&mut self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(3, false, Some(s), false)
    }

    /// `RgNLineHCompound()` (`lxx:126-129`).
    pub fn rgn_line_h_compound(&mut self) -> Option<TopoShape> {
        self.internal_compound(4, false, None, false)
    }

    /// `RgNLineHCompound(S)` (`lxx:133-136`).
    pub fn rgn_line_h_compound_of(&mut self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(4, false, Some(s), false)
    }

    /// `OutLineHCompound()` (`lxx:140-143`).
    pub fn out_line_h_compound(&mut self) -> Option<TopoShape> {
        self.internal_compound(2, false, None, false)
    }

    /// `OutLineHCompound(S)` (`lxx:148-151`).
    pub fn out_line_h_compound_of(&mut self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(2, false, Some(s), false)
    }

    /// `IsoLineHCompound()` (`lxx:155-158`).
    pub fn iso_line_h_compound(&mut self) -> Option<TopoShape> {
        self.internal_compound(1, false, None, false)
    }

    /// `IsoLineHCompound(S)` (`lxx:162-165`).
    pub fn iso_line_h_compound_of(&mut self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(1, false, Some(s), false)
    }

    /// `CompoundOfEdges(type, visible, In3d)` (`lxx:169-174`).
    pub fn compound_of_edges(
        &mut self,
        r#type: HlrBRepTypeOfResultingEdge,
        visible: bool,
        in3d: bool,
    ) -> Option<TopoShape> {
        self.internal_compound(r#type.as_int(), visible, None, in3d)
    }

    /// `CompoundOfEdges(S, type, visible, In3d)` (`lxx:178-183`).
    pub fn compound_of_edges_of(
        &mut self,
        s: &TopoShape,
        r#type: HlrBRepTypeOfResultingEdge,
        visible: bool,
        in3d: bool,
    ) -> Option<TopoShape> {
        self.internal_compound(r#type.as_int(), visible, Some(s), in3d)
    }
}

/// The two `HLRBRep::MakeEdge` / `MakeEdge3d` calls of `cxx:269-274` and
/// `:289-294`, sharing the `if (!E.IsNull())` test of `cxx:276` / `:296`.
fn make_result_edge(
    ed: &HlrBRepEdgeData,
    sta: f64,
    end: f64,
    in3d: bool,
) -> Option<crate::shape::Edge> {
    if in3d {
        hlrbrep::make_edge_3d(ed.geometry(), sta, end)
    } else {
        hlrbrep::make_edge(ed.geometry(), sta, end)
    }
}
