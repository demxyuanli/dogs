//! `HLRBRep_InternalAlgo` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_InternalAlgo.hxx:34-139`,
//! `.cxx:50-1018`): the driver that stacks several shapes into one
//! `HLRBRep_Data` and then selects / shows / hides the edges of each.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `HLRBRep_InternalAlgo()` | `cxx:50-53` | [`HlrBRepInternalAlgo::new`] |
//! | `HLRBRep_InternalAlgo(A)` | `cxx:57-64` | [`HlrBRepInternalAlgo::from`] |
//! | `Projector(P)` / `Projector()` | `cxx:68-76` | [`HlrBRepInternalAlgo::set_projector`] / [`HlrBRepInternalAlgo::projector`] |
//! | `Load(S, SData, nbIso)` | `cxx:187-194` | [`HlrBRepInternalAlgo::load_with_data`] |
//! | `Load(S, nbIso)` | `cxx:198-202` | [`HlrBRepInternalAlgo::load`] |
//! | `Index(S)` | `cxx:206-220` | [`HlrBRepInternalAlgo::index`] |
//! | `Remove(I)` | `cxx:224-232` | [`HlrBRepInternalAlgo::remove`] |
//! | `ShapeData(I, SData)` | `cxx:236-242` | [`HlrBRepInternalAlgo::set_shape_data`] |
//! | `SeqOfShapeBounds()` | `cxx:246-249` | [`HlrBRepInternalAlgo::seq_of_shape_bounds`] |
//! | `NbShapes()` | `cxx:253-256` | [`HlrBRepInternalAlgo::nb_shapes`] |
//! | `ShapeBounds(I)` | `cxx:260-266` | [`HlrBRepInternalAlgo::shape_bounds`] / `_mut` |
//! | `InitEdgeStatus()` | `cxx:270-348` | [`HlrBRepInternalAlgo::init_edge_status`] |
//! | `Select()` | `cxx:352-373` | [`HlrBRepInternalAlgo::select`] |
//! | `Select(I)` | `cxx:377-405` | [`HlrBRepInternalAlgo::select_shape`] |
//! | `SelectEdge(I)` | `cxx:409-427` | [`HlrBRepInternalAlgo::select_edge`] |
//! | `SelectFace(I)` | `cxx:431-450` | [`HlrBRepInternalAlgo::select_face`] |
//! | `ShowAll()` / `ShowAll(I)` | `cxx:454-492` | [`HlrBRepInternalAlgo::show_all`] / [`HlrBRepInternalAlgo::show_all_of`] |
//! | `HideAll()` / `HideAll(I)` | `cxx:496-534` | [`HlrBRepInternalAlgo::hide_all`] / [`HlrBRepInternalAlgo::hide_all_of`] |
//! | `PartialHide()` | `cxx:538-556` | [`HlrBRepInternalAlgo::partial_hide`] |
//! | `Hide()` | `cxx:560-588` | [`HlrBRepInternalAlgo::hide`] |
//! | `Hide(I)` | `cxx:592-608` | [`HlrBRepInternalAlgo::hide_of`] |
//! | `Hide(I, J)` | `cxx:613-653` | [`HlrBRepInternalAlgo::hide_by`] |
//! | `HideSelected(I, SideFace)` | `cxx:658-1004` | [`HlrBRepInternalAlgo::hide_selected`] |
//! | `Debug` | `cxx:1006-1018` | [`HlrBRepInternalAlgo::set_debug`] / [`HlrBRepInternalAlgo::debug`] |
//! | `DataStructure()` | `hxx:122` | [`HlrBRepInternalAlgo::data_structure`] |
//!
//! # UNPORTED
//!
//! * The `myDebug` `std::cout` transcripts (`cxx:545-548`, `:567-570`,
//!   `:601-604`, `:643-647`, `:689-748`, `:766-790`, `:937-963`, `:968-996`)
//!   are `OCCT_DEBUG` progress output; `debug` only stores the flag as OCCT
//!   does.
//!
//! # Translation notes
//!
//! * `NCollection_Sequence<HLRBRep_ShapeBounds> myShapes` (`hxx:135`) is a
//!   `Vec` with the OCCT 1-based index translated by [`HlrBRepInternalAlgo::shape_bounds`];
//!   the raisers of `cxx:227`, `:239`, `:263` and `:380` become `None`.
//! * `myMapOfShapeTool` (`hxx:136`) is the domain-tool cache
//!   [`ShapeToTopAdaptorToolMap`] keyed by `shape_map_key`; `Update` fills it
//!   through `HLRBRep_ShapeToHLR::Load`.
//! * `myDS` (`hxx:133`) is `Option<Arc<HlrBRepData>>`; `Nullify()` is
//!   `self.my_ds = None`.

use std::sync::Arc;

use crate::abs::Orientation;
use crate::hlr_algo::{HlrAlgoProjector, MinMaxIndices};
use crate::hlrbrep_data::{HlrBRepData, MINMAX_SIGN_MASK};
use crate::hlrbrep_face_iterator::HlrBRepFaceIterator;
use crate::hlrbrep_hider::HlrBRepHider;
use crate::hlrbrep_shape_bounds::HlrBRepShapeBounds;
use crate::hlrtopobrep::out_liner::{HlrTopoBRepOutLiner, ShapeToTopAdaptorToolMap};
use std::any::Any;

/// `HLRBRep_InternalAlgo`.
pub struct HlrBRepInternalAlgo {
    /// `occ::handle<HLRBRep_Data> myDS` (`hxx:133`).
    my_ds: Option<Arc<HlrBRepData>>,
    /// `HLRAlgo_Projector myProj` (`hxx:134`).
    my_proj: Arc<HlrAlgoProjector>,
    /// `NCollection_Sequence<HLRBRep_ShapeBounds> myShapes` (`hxx:135`).
    my_shapes: Vec<HlrBRepShapeBounds>,
    /// `NCollection_DataMap<TopoDS_Shape, BRepTopAdaptor_Tool,
    /// TopTools_ShapeMapHasher> myMapOfShapeTool` (`hxx:136`): the domain tools
    /// `HLRBRep_ShapeToHLR::Load` fills (`cxx:95`).
    my_map_of_shape_tool: ShapeToTopAdaptorToolMap,
    /// `bool myDebug` (`hxx:137`).
    my_debug: bool,
}

impl HlrBRepInternalAlgo {
    /// `HLRBRep_InternalAlgo()` (`cxx:50-53`): `myDebug(false)`, the projector
    /// and the shape sequence left as their C++ default members leave them.
    pub fn new() -> Self {
        Self {
            my_ds: None,
            my_proj: Arc::new(HlrAlgoProjector::new()),
            my_shapes: Vec::new(),
            my_map_of_shape_tool: ShapeToTopAdaptorToolMap::new(),
            my_debug: false,
        }
    }

    /// `HLRBRep_InternalAlgo(const handle(HLRBRep_InternalAlgo)& A)`
    /// (`cxx:57-64`): shares the data structure and the shapes of `a`.
    pub fn from(a: &HlrBRepInternalAlgo) -> Self {
        Self {
            my_ds: a.my_ds.clone(),
            my_proj: Arc::clone(&a.my_proj),
            my_shapes: a.my_shapes.clone(),
            my_map_of_shape_tool: a.my_map_of_shape_tool.clone(),
            my_debug: a.my_debug,
        }
    }

    /// `Projector(const HLRAlgo_Projector& P)` (`cxx:68-71`).
    pub fn set_projector(&mut self, p: Arc<HlrAlgoProjector>) {
        self.my_proj = p;
    }

    /// `HLRAlgo_Projector& Projector()` (`cxx:74-77`).
    pub fn projector(&self) -> &Arc<HlrAlgoProjector> {
        &self.my_proj
    }

    /// The non-const `Projector()` of `cxx:74-77`.
    pub fn projector_mut(&mut self) -> &mut Arc<HlrAlgoProjector> {
        &mut self.my_proj
    }

    /// `Update()` (`cxx:82-186`).
    ///
    /// For every `HLRBRep_ShapeBounds` it runs
    /// [`crate::hlrbrep_shape_to_hlr::load`] (`cxx:95`), rewrites the bounds to
    /// `(1, dv, 1, de, 1, df)` (`cxx:117`) and accumulates the three sizes.
    /// A single shape becomes `myDS`; several are spliced through
    /// `HLRBRep_Data::Write` (`cxx:134-144`). Then `myDS->Update(myProj)`
    /// (`cxx:148`) and the per-shape `MinMax` fold of `cxx:150-184`.
    ///
    /// The `catch (Standard_Failure)` arm of `cxx:98-114` becomes an `Err` from
    /// `load`, which stores an empty `HLRBRep_Data(0, 0, 0)` exactly as OCCT
    /// does.
    pub fn update(&mut self) -> Result<(), &'static str> {
        if self.my_shapes.is_empty() {
            return Ok(()); // `if (!myShapes.IsEmpty())` (`cxx:84`)
        }
        let n = self.my_shapes.len();
        let proj = Arc::clone(&self.my_proj);
        let mut loaded: Vec<HlrBRepData> = Vec::with_capacity(n);
        let (mut nv, mut ne, mut nf) = (0, 0, 0);

        // `cxx:89-121`.
        let mut shapes = std::mem::take(&mut self.my_shapes);
        for sb in shapes.iter_mut() {
            let nb_iso = sb.nb_of_iso();
            let (dv, de, df) = match sb.take_shape() {
                Some(mut arc) => {
                    // OCCT fills through `handle<HLRTopoBRep_OutLiner>`; the port
                    // mutates the shared outliner in place, cloning only when
                    // another owner still holds it.
                    let fill_result = {
                        let s = Arc::make_mut(&mut arc);
                        crate::hlrbrep_shape_to_hlr::load(
                            s,
                            &proj,
                            &mut self.my_map_of_shape_tool,
                            nb_iso,
                        )
                    };
                    sb.set_shape(arc);
                    match fill_result {
                        Ok(ds) => {
                            let sizes = (ds.nb_vertices(), ds.nb_edges(), ds.nb_faces());
                            loaded.push(ds);
                            sizes
                        }
                        Err(_) => {
                            loaded.push(HlrBRepData::new(0, 0, 0));
                            (0, 0, 0)
                        }
                    }
                }
                None => {
                    loaded.push(HlrBRepData::new(0, 0, 0));
                    (0, 0, 0)
                }
            };
            // `SB = HLRBRep_ShapeBounds(SB.Shape(), SB.ShapeData(),
            // SB.NbOfIso(), 1, dv, 1, de, 1, df);` (`cxx:117`).
            let shape = sb.shape().cloned();
            let s_data = sb.shape_data().cloned();
            *sb = match shape {
                Some(shape) => HlrBRepShapeBounds::with_shape_data(
                    shape, s_data, nb_iso, 1, dv, 1, de, 1, df,
                ),
                None => HlrBRepShapeBounds::new_bounded(
                    Arc::new(HlrTopoBRepOutLiner::new()),
                    nb_iso,
                    1,
                    dv,
                    1,
                    de,
                    1,
                    df,
                ),
            };
            nv += dv;
            ne += de;
            nf += df;
        }

        // `cxx:124-145`.
        let merged = if n == 1 {
            loaded.pop().expect("one loaded data structure")
        } else {
            let mut ds = HlrBRepData::new(nv, ne, nf);
            let (mut v, mut e, mut f) = (0, 0, 0);
            for (i, sb) in shapes.iter_mut().enumerate() {
                let (dv, de, df) = sb.sizes();
                sb.translate(v, e, f);
                ds.write(&loaded[i], v, e, f);
                v += dv;
                e += de;
                f += df;
            }
            ds
        };
        self.my_shapes = shapes;
        self.my_ds = Some(Arc::new(merged));

        // `myDS->Update(myProj);` (`cxx:148`).
        if let Some(ds) = self.my_ds.as_mut() {
            let proj = Arc::clone(&proj);
            Arc::make_mut(ds).update(&proj)?;
        }

        // `cxx:150-184`: the per-shape `MinMax` fold.
        let Some(ds) = self.my_ds.as_ref() else {
            return Ok(());
        };
        for sb in self.my_shapes.iter_mut() {
            let (v1, v2, e1, e2, f1, f2) = sb.bounds();
            let _ = (v1, v2);
            let mut min_max = MinMaxIndices::default();
            crate::hlrbrep_shape_to_hlr::update_min_max_of_shape(
                ds,
                e1,
                e2,
                f1,
                f2,
                &mut min_max,
            );
            sb.update_min_max(&min_max);
        }
        Ok(())
    }

    /// `Load(S, SData, nbIso)` (`cxx:187-194`).
    pub fn load_with_data(
        &mut self,
        s: Arc<HlrTopoBRepOutLiner>,
        s_data: Option<Arc<dyn Any + Send + Sync>>,
        nb_iso: i32,
    ) {
        self.my_shapes
            .push(HlrBRepShapeBounds::with_shape_data(s, s_data, nb_iso, 0, 0, 0, 0, 0, 0));
        self.my_ds = None; // `myDS.Nullify()` (`cxx:193`)
    }

    /// `Load(S, nbIso)` (`cxx:198-202`).
    pub fn load(&mut self, s: Arc<HlrTopoBRepOutLiner>, nb_iso: i32) {
        self.my_shapes
            .push(HlrBRepShapeBounds::new_bounded(s, nb_iso, 0, 0, 0, 0, 0, 0));
        self.my_ds = None;
    }

    /// `Index(S)` (`cxx:206-220`): the 1-based index of `s`, `0` if absent.
    pub fn index(&self, s: &Arc<HlrTopoBRepOutLiner>) -> i32 {
        for (i, sb) in self.my_shapes.iter().enumerate() {
            if let Some(shape) = sb.shape() {
                if Arc::ptr_eq(shape, s) {
                    return i as i32 + 1;
                }
            }
        }
        0 // `cxx:219`
    }

    /// `Remove(I)` (`cxx:224-232`). OCCT raises `Standard_OutOfRange` for an
    /// index outside `1..=Length()`; here the call is ignored instead.
    pub fn remove(&mut self, i: i32) {
        if i == 0 || i as usize > self.my_shapes.len() {
            return; // `Standard_OutOfRange_Raise_if` (`cxx:227-228`)
        }
        self.my_shapes.remove(i as usize - 1);
        self.my_ds = None;
    }

    /// `ShapeData(I, SData)` (`cxx:236-242`).
    pub fn set_shape_data(&mut self, i: i32, s_data: Option<Arc<dyn Any + Send + Sync>>) {
        if i == 0 || i as usize > self.my_shapes.len() {
            return; // `Standard_OutOfRange_Raise_if` (`cxx:239-240`)
        }
        self.my_shapes[i as usize - 1].set_shape_data(s_data);
    }

    /// `SeqOfShapeBounds()` (`cxx:246-249`).
    pub fn seq_of_shape_bounds(&self) -> &[HlrBRepShapeBounds] {
        &self.my_shapes
    }

    /// The mutable `SeqOfShapeBounds()` of `cxx:246-249`.
    pub fn seq_of_shape_bounds_mut(&mut self) -> &mut Vec<HlrBRepShapeBounds> {
        &mut self.my_shapes
    }

    /// `NbShapes()` (`cxx:253-256`).
    pub fn nb_shapes(&self) -> i32 {
        self.my_shapes.len() as i32
    }

    /// `ShapeBounds(I)` (`cxx:260-266`), `None` where OCCT raises.
    pub fn shape_bounds(&self, i: i32) -> Option<&HlrBRepShapeBounds> {
        if i < 1 {
            return None;
        }
        self.my_shapes.get(i as usize - 1)
    }

    /// The mutable `ShapeBounds(I)` of `cxx:260-266`.
    pub fn shape_bounds_mut(&mut self, i: i32) -> Option<&mut HlrBRepShapeBounds> {
        if i < 1 {
            return None;
        }
        self.my_shapes.get_mut(i as usize - 1)
    }

    /// `InitEdgeStatus()` (`cxx:270-348`): initialise the visibility status of
    /// the selected edges from the back faces of a closed shell.
    ///
    /// OCCT dereferences `myDS` without a test; here a missing data structure
    /// is the `if (!myDS.IsNull())` guard the other entries carry
    /// (`cxx:354`, `:464`, ...).
    pub fn init_edge_status(&mut self) {
        let Some(ds) = self.my_ds.as_mut() else {
            return;
        };
        let data = Arc::make_mut(ds);

        let ne = data.nb_edges();
        let nf = data.nb_faces();

        // `cxx:280-287`: every selected edge starts fully visible.
        for e in 1..=ne {
            let ed = data.e_data_mut(e);
            if ed.selected() {
                ed.status_mut().show_all();
            }
        }

        let mut face_itr = HlrBRepFaceIterator::new();

        // `cxx:289-306`: the edges of a selected face start hidden.
        for f in 1..=nf {
            let fd = data.f_data(f);
            if fd.selected() {
                face_itr.init_edge(fd);
                while face_itr.more_edge() {
                    let e = face_itr.edge();
                    let edf = data.e_data_mut(e);
                    if edf.selected() {
                        edf.status_mut().hide_all();
                    }
                    face_itr.next_edge();
                }
            }
        }

        // `cxx:309-347`: re-show the edges of the faces that stay visible.
        for f in 1..=nf {
            let fd = data.f_data(f);
            let mut visible = true;
            if fd.selected() && fd.closed() {
                if fd.side() {
                    visible = false;
                } else if !fd.with_out_l() {
                    match fd.orientation() {
                        Orientation::Reversed => visible = fd.back(),
                        Orientation::Forward => visible = !fd.back(),
                        // `TopAbs_EXTERNAL` / `TopAbs_INTERNAL` (`cxx:327-331`).
                        Orientation::External | Orientation::Internal => visible = true,
                    }
                }
            }
            if visible {
                face_itr.init_edge(fd);
                while face_itr.more_edge() {
                    let e = face_itr.edge();
                    let edf = data.e_data_mut(e);
                    if edf.selected() && !edf.vertical() {
                        edf.status_mut().show_all();
                    }
                    face_itr.next_edge();
                }
            }
        }
    }

    /// `Select()` (`cxx:352-373`): select the whole data structure.
    pub fn select(&mut self) {
        let Some(ds) = self.my_ds.as_mut() else {
            return;
        };
        let data = Arc::make_mut(ds);
        let ne = data.nb_edges();
        let nf = data.nb_faces();
        for e in 1..=ne {
            data.e_data_mut(e).set_selected(true);
        }
        for f in 1..=nf {
            data.f_data_mut(f).set_selected(true);
        }
    }

    /// `Select(I)` (`cxx:377-405`): select only the edges and faces of shape
    /// `i`.
    pub fn select_shape(&mut self, i: i32) {
        let Some((e1, e2, f1, f2)) = self.bounds_of(i) else {
            return;
        };
        let Some(ds) = self.my_ds.as_mut() else {
            return;
        };
        let data = Arc::make_mut(ds);
        let ne = data.nb_edges();
        let nf = data.nb_faces();
        for e in 1..=ne {
            data.e_data_mut(e).set_selected(e >= e1 && e <= e2);
        }
        for f in 1..=nf {
            data.f_data_mut(f).set_selected(f >= f1 && f <= f2);
        }
    }

    /// `SelectEdge(I)` (`cxx:409-427`): select only the edges of shape `i`.
    pub fn select_edge(&mut self, i: i32) {
        let Some((e1, e2, _f1, _f2)) = self.bounds_of(i) else {
            return;
        };
        let Some(ds) = self.my_ds.as_mut() else {
            return;
        };
        let data = Arc::make_mut(ds);
        let ne = data.nb_edges();
        for e in 1..=ne {
            data.e_data_mut(e).set_selected(e >= e1 && e <= e2);
        }
    }

    /// `SelectFace(I)` (`cxx:431-450`): select only the faces of shape `i`.
    pub fn select_face(&mut self, i: i32) {
        let Some((_e1, _e2, f1, f2)) = self.bounds_of(i) else {
            return;
        };
        let Some(ds) = self.my_ds.as_mut() else {
            return;
        };
        let data = Arc::make_mut(ds);
        let nf = data.nb_faces();
        for f in 1..=nf {
            data.f_data_mut(f).set_selected(f >= f1 && f <= f2);
        }
    }

    /// `ShowAll()` (`cxx:454-468`).
    pub fn show_all(&mut self) {
        let Some(ds) = self.my_ds.as_mut() else {
            return;
        };
        let data = Arc::make_mut(ds);
        let ne = data.nb_edges();
        for e in 1..=ne {
            data.e_data_mut(e).status_mut().show_all();
        }
    }

    /// `ShowAll(I)` (`cxx:472-492`).
    pub fn show_all_of(&mut self, i: i32) {
        if self.my_ds.is_none() {
            return;
        }
        if i == 0 || i as usize > self.my_shapes.len() {
            return; // `cxx:477-478`
        }
        self.select_shape(i);
        let Some(ds) = self.my_ds.as_mut() else {
            return;
        };
        let data = Arc::make_mut(ds);
        let ne = data.nb_edges();
        for e in 1..=ne {
            let ed = data.e_data_mut(e);
            if ed.selected() {
                ed.status_mut().show_all();
            }
        }
    }

    /// `HideAll()` (`cxx:496-508`).
    pub fn hide_all(&mut self) {
        let Some(ds) = self.my_ds.as_mut() else {
            return;
        };
        let data = Arc::make_mut(ds);
        let ne = data.nb_edges();
        for e in 1..=ne {
            data.e_data_mut(e).status_mut().hide_all();
        }
    }

    /// `HideAll(I)` (`cxx:512-534`).
    pub fn hide_all_of(&mut self, i: i32) {
        if self.my_ds.is_none() {
            return;
        }
        if i == 0 || i as usize > self.my_shapes.len() {
            return; // `cxx:517-518`
        }
        self.select_shape(i);
        let Some(ds) = self.my_ds.as_mut() else {
            return;
        };
        let data = Arc::make_mut(ds);
        let ne = data.nb_edges();
        for e in 1..=ne {
            let ed = data.e_data_mut(e);
            if ed.selected() {
                ed.status_mut().hide_all();
            }
        }
    }

    /// `PartialHide()` (`cxx:538-556`): `Hide(i)` for every shape, then
    /// `Select()`.
    pub fn partial_hide(&mut self) {
        if self.my_ds.is_none() {
            return;
        }
        let n = self.my_shapes.len() as i32;
        for i in 1..=n {
            self.hide_of(i);
        }
        self.select();
    }

    /// `Hide()` (`cxx:560-588`): every shape by itself, then every ordered pair
    /// `i != j`, then `Select()`.
    pub fn hide(&mut self) {
        if self.my_ds.is_none() {
            return;
        }
        let n = self.my_shapes.len() as i32;
        for i in 1..=n {
            self.hide_of(i);
        }
        for i in 1..=n {
            for j in 1..=n {
                if i != j {
                    self.hide_by(i, j);
                }
            }
        }
        self.select();
    }

    /// `Hide(I)` (`cxx:592-608`): `Select(I)`, `InitEdgeStatus()`,
    /// `HideSelected(I, true)`.
    pub fn hide_of(&mut self, i: i32) {
        if self.my_ds.is_none() {
            return;
        }
        if i == 0 || i as usize > self.my_shapes.len() {
            return; // `cxx:595-596`
        }
        self.select_shape(i);
        self.init_edge_status();
        self.hide_selected(i, true);
    }

    /// `Hide(I, J)` (`cxx:613-653`): hide shape `i` by shape `j`.
    ///
    /// The twelve `& 0x80008000` pair tests (`cxx:627-641`) reject the pair
    /// before the selection; the accepted arm selects the edges of `i` and the
    /// faces of `j` and hides with `SideFace = false` (`cxx:649-651`).
    pub fn hide_by(&mut self, i: i32, j: i32) {
        if self.my_ds.is_none() {
            return;
        }
        if i == 0 || i as usize > self.my_shapes.len() || j == 0 || j as usize > self.my_shapes.len()
        {
            return; // `cxx:617-618`
        }
        if i == j {
            return self.hide_of(i); // `cxx:621-623`
        }

        let (Some(sbi), Some(sbj)) = (self.shape_bounds(i), self.shape_bounds(j)) else {
            return;
        };
        let min_max_shbi: MinMaxIndices = *sbi.min_max();
        let min_max_shbj: MinMaxIndices = *sbj.min_max();

        // `cxx:627-642`: the sixteen `& 0x80008000` pair tests.
        let mut overlaps = true;
        for k in 0..8 {
            if (min_max_shbj.max[k] - min_max_shbi.min[k]) & MINMAX_SIGN_MASK != 0
                || (min_max_shbi.max[k] - min_max_shbj.min[k]) & MINMAX_SIGN_MASK != 0
            {
                overlaps = false;
                break;
            }
        }
        if !overlaps {
            return;
        }

        self.select_edge(i);
        self.select_face(j);
        self.hide_selected(i, false);
    }

    /// `HideSelected(I, SideFace)` (`cxx:658-1004`): sort the selected hiding
    /// faces of shape `i` by surface type then size (`cxx:795-927`) and let
    /// [`HlrBRepHider::hide`] run over them.
    ///
    /// The `OCCT_DEBUG` counters and transcripts (`cxx:660-674`, `:689-748`,
    /// `:766-790`, `:937-963`, `:968-996`) are not ported.
    pub fn hide_selected(&mut self, i: i32, side_face: bool) {
        // `myDS` is moved out so the hider can mutate it while
        // `my_map_of_shape_tool` is borrowed mutably at the same time.
        let Some(ds) = self.my_ds.take() else {
            return;
        };
        let mut ds = ds;
        {
            let data = Arc::make_mut(&mut ds);
            self.hide_selected_in(data, i, side_face);
        }
        self.my_ds = Some(ds);
    }

    /// The body of `HideSelected` over an owned data structure.
    fn hide_selected_in(&mut self, data: &mut HlrBRepData, i: i32, side_face: bool) {
        let Some((e1, e2, _f1, _f2)) = self.bounds_of(i) else {
            return; // `Standard_OutOfRange_Raise_if` of `cxx:678`
        };
        if e2 < e1 {
            return;
        }

        {
            let Some(sb) = self.my_shapes.get(i as usize - 1) else {
                return;
            };
            let min_max = *sb.min_max();
            data.init_bound_sort(&min_max, e1, e2);
        }

        let mut hider = HlrBRepHider::new();
        let nf = data.nb_faces();
        if nf == 0 {
            return;
        }

        if side_face {
            for f in 1..=nf {
                let fd = data.f_data(f);
                if fd.selected() && fd.side() {
                    hider.own_hiding(f);
                }
            }
        }

        // `Val` / `Size` of `cxx:795-844`, 1-based like the OCCT arrays.
        let n = nf as usize;
        let mut val = vec![0i32; n + 1];
        let mut size = vec![0.0f64; n + 1];
        let mut index = vec![0i32; n + 1];
        for f in 1..=nf {
            let fd = data.f_data(f);
            let mut v = if fd.plane() {
                10
            } else if fd.cylinder() {
                9
            } else if fd.cone() {
                8
            } else if fd.sphere() {
                7
            } else if fd.torus() {
                6
            } else {
                0
            };
            if fd.cut() {
                v -= 10;
            }
            if fd.side() {
                v -= 100;
            }
            if fd.with_out_l() {
                v -= 20;
            }
            val[f as usize] = v;
            size[f as usize] = fd.size();
        }
        for (tt, slot) in index.iter_mut().enumerate().take(nf as usize + 1).skip(1) {
            *slot = tt as i32;
        }

        // `cxx:869-927`: the heapsort, verbatim with 1-based indices.
        if nf > 2 {
            let mut l = (nf >> 1) + 1;
            let mut ir = nf;
            let mut rra: i32;
            loop {
                if l > 1 {
                    l -= 1;
                    rra = index[l as usize];
                } else {
                    rra = index[ir as usize];
                    index[ir as usize] = index[1];
                    ir -= 1;
                    if ir == 1 {
                        index[1] = rra;
                        break;
                    }
                }
                let mut ii = l;
                let mut k = l + l;
                while k <= ir {
                    if k < ir {
                        let a = index[k as usize] as usize;
                        let b = index[(k + 1) as usize] as usize;
                        if val[a] > val[b] {
                            k += 1;
                        } else if val[a] == val[b] && size[a] > size[b] {
                            k += 1;
                        }
                    }
                    let kk = index[k as usize] as usize;
                    let rr = rra as usize;
                    let take = val[rr] > val[kk] || (val[rr] == val[kk] && size[rr] > size[kk]);
                    if take {
                        index[ii as usize] = index[k as usize];
                        ii = k;
                        k <<= 1;
                    } else {
                        k = ir + 1;
                    }
                }
                index[ii as usize] = rra;
            }
        }

        for f in 1..=nf {
            let fi = index[f as usize];
            let fd = data.f_data(fi);
            if fd.selected() && fd.hiding() {
                hider.hide(data, fi, &mut self.my_map_of_shape_tool);
            }
        }
    }

    /// `Debug(bool)` (`cxx:1006-1009`).
    pub fn set_debug(&mut self, deb: bool) {
        self.my_debug = deb;
    }

    /// `Debug()` (`cxx:1013-1016`).
    pub fn debug(&self) -> bool {
        self.my_debug
    }

    /// `DataStructure()` (`hxx:122`).
    pub fn data_structure(&self) -> Option<&Arc<HlrBRepData>> {
        self.my_ds.as_ref()
    }

    /// The mutable `myDS`, for the driver layer that fills it.
    pub fn set_data_structure(&mut self, ds: Option<Arc<HlrBRepData>>) {
        self.my_ds = ds;
    }

    /// The mutable `myDS` handle itself, for the callers that rewrite the data
    /// structure in place: `HLRBRep_HLRToShape` mutates `EDataArray` `Used` /
    /// `HideCount` and `Projector` through `myAlgo->DataStructure()`
    /// (`HLRBRep_HLRToShape.cxx:76-143`).
    pub fn data_structure_mut(&mut self) -> Option<&mut Arc<HlrBRepData>> {
        self.my_ds.as_mut()
    }

    /// `myShapes(I).Bounds(v1, v2, e1, e2, f1, f2)` as `(e1, e2, f1, f2)`.
    /// Callers test `myDS` first, exactly as OCCT does, and the range check is
    /// the `Standard_OutOfRange_Raise_if` of each entry.
    fn bounds_of(&self, i: i32) -> Option<(i32, i32, i32, i32)> {
        if i == 0 || i as usize > self.my_shapes.len() {
            return None;
        }
        let (_v1, _v2, e1, e2, f1, f2) = self.my_shapes[i as usize - 1].bounds();
        Some((e1, e2, f1, f2))
    }
}

impl Default for HlrBRepInternalAlgo {
    fn default() -> Self {
        Self::new()
    }
}
