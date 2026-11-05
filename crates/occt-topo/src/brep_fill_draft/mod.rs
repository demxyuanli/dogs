//! `BRepFill_Draft` — the skin a wire sweeps out at a constant draft angle
//! (`BRepFill_Draft.cxx:66-975`, `.hxx:16-101`, OCCT 8.0.0).
//!
//! Port status against the reference:
//!
//! * Ported: the constructor's three input arms (`cxx:195-262`), the
//!   "closed but not declared" wire fix-up (`cxx:266-276`), `SetOptions`
//!   (`cxx:286-292`), `SetDraft` (`cxx:294-302`), the three `Perform`
//!   overloads (`cxx:306-407`), `Init` (`cxx:412-460`), `BuildShell`
//!   (`cxx:467-531`), `Fuse` (`cxx:538-824`), `Sewing` (`cxx:828-916`) and the
//!   four accessors (`cxx:924-975`).
//! * The vertex arm of `Generated` (`cxx:929-938`) is unreachable:
//!   `TopoDS::Edge(S)` raises `Standard_TypeMismatch` for a shape that is not
//!   an edge (`TopoDS.hxx:92-96`), so `E` is never null and the loop can never
//!   run. The port reports the mismatch instead.
//!
//! Translation notes:
//! - `handle(BRepFill_DraftLaw)` is `Arc<DraftLaw>` so `Init` / `GoodOrientation`
//!   / `Generated` reach the concrete class; `handle(BRepFill_SectionLaw)` is
//!   `Arc<dyn SectionLaw>`.
//! - `Bnd_Box` is [`BndBox`]. `Bnd_Box::Get(Xmin, Ymin, Zmin, Xmax, Ymax,
//!   Zmax)` is [`BndBox::get`], whose tuple is `(xmin, xmax, ymin, ymax, zmin,
//!   zmax)` — see [`geometry::bounds`].
//! - `TopLoc_Location(Tf)` applied to the wire (`cxx:92-95`) is
//!   [`crate::shape_ops::transformed_copy`], which bakes the transformation
//!   into fresh geometry; `BndLib` reads the same curves either way.
//! - `NCollection_HArray2<TopoDS_Shape>` is
//!   [`crate::brep_fill_sweep::helpers::ShapeTable`]. The sweep's tables are
//!   copied out of it because `Fuse` / `Sewing` rewrite the cells in place
//!   (`myFaces->SetValue`, `cxx:810`, `cxx:881`).
//! - A null `TopoDS_Wire` / `TopoDS_Shell` / `TopoDS_Solid` / `TopoDS_Shape` is
//!   `None`; `myTop` is never null because the constructor rejected every
//!   usable shape type.
//!
//! Two OCCT call sequences have no one-to-one counterpart in the BOP port and
//! are kept equivalent rather than approximated:
//!
//! - `Fuse` performs one `BOPAlgo_PaveFiller` and then hands it to
//!   `BRepAlgoAPI_Section` (`cxx:599`) and to `PerformWithFiller`
//!   (`cxx:694`). The port runs the same intersection phase inside each of
//!   those algorithms ([`crate::bop_builder2::BopBuilder::perform_section`] and
//!   [`crate::bop_builder2::BopBuilder::perform`]); `PaveFiller::perform` is a
//!   function of the arguments and options alone, so every run fills the same
//!   `BopdsDS`.
//! - `aHistory->Merge(aLO, aBuilder)` restricts the cut history to the first
//!   argument (`BRepTools_History.hxx:209-214`). The port merges the whole cut
//!   history ([`BopHistory::merge_algorithm`]); the only shapes read back
//!   afterwards are the section/face tables, which live under that argument.

mod build_shell;
mod fuse;
mod geometry;
mod init;
mod sewing;

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::{GpAx3, GpDir, GpPln, GpPnt};
use occt_core::kernel::geomabs::Shape;
use occt_geom::rectangular_trimmed::GeomRectangularTrimmedSurface;

use crate::abs::ShapeType;
use crate::brep_fill_draft_law::DraftLaw;
use crate::brep_fill_section_law::SectionLaw;
use crate::brep_fill_sweep::helpers::ShapeTable;
use crate::brep_fill_sweep::TransitionStyle;
use crate::brep_tool::BRepTool;
use crate::shape::{Face, Shell, TopoShape, Wire};
use crate::topo_tools::{is_same, shape_map_key, wire_vertices, ShapeMapKey};

/// `Standard_ConstructionError("BRepFill_Draft")` (`cxx:248`, `cxx:255`,
/// `cxx:261`).
pub const CONSTRUCTION_ERROR: &str = "Standard_ConstructionError: BRepFill_Draft";

/// `BRepFill_Draft` (`BRepFill_Draft.hxx:37-101`): the class data.
pub struct BrepFillDraft {
    /// `myDir` (`hxx:82`). Reversed in place by [`geometry::longueur`]
    /// (`cxx:124`) and read by every later step.
    my_dir: GpDir,
    /// `myAngle` (`hxx:83`). Negated by `Init` when the law runs the other way
    /// (`cxx:423`).
    my_angle: f64,
    /// `angmin` (`hxx:84`): `SetOptions`' `Min`.
    angmin: f64,
    /// `angmax` (`hxx:85`): `SetOptions`' `Max`.
    angmax: f64,
    /// `myTol` (`hxx:86`): `1.e-4` (`cxx:279`).
    my_tol: f64,
    /// `myLoc` (`hxx:87`): the draft location law built by `Init`.
    my_loc: Option<Arc<DraftLaw>>,
    /// `mySec` (`hxx:88`): the generating section law.
    my_sec: Option<Arc<dyn SectionLaw>>,
    /// `mySections` (`hxx:89`): `1 x (NbLaw + 1)` section shapes.
    my_sections: ShapeTable,
    /// `myFaces` (`hxx:90`): `1 x NbLaw` swept faces.
    my_faces: ShapeTable,
    /// `myGenerated` (`hxx:91`).
    my_generated: Vec<TopoShape>,
    /// `myShape` (`hxx:92`): null is `None`.
    my_shape: Option<TopoShape>,
    /// `myTop` (`hxx:93`): the input shape.
    my_top: TopoShape,
    /// `myShell` (`hxx:94`): the skin without the stopping face.
    my_shell: Option<Shell>,
    /// `myWire` (`hxx:95`): the generating wire.
    my_wire: Wire,
    /// `myCont` (`hxx:96`): `GeomAbs_C1` (`cxx:280`).
    my_cont: Shape,
    /// `myStyle` (`hxx:97`): `BRepFill_Right` (`hxx:38`).
    my_style: TransitionStyle,
    /// `IsInternal` (`hxx:98`).
    is_internal: bool,
    /// `myDone` (`hxx:99`).
    my_done: bool,
}

impl BrepFillDraft {
    /// `BRepFill_Draft(const TopoDS_Shape& S, const gp_Dir& Dir, const double
    /// Angle)` (`cxx:193-284`).
    ///
    /// The `throw Standard_ConstructionError` arms of the reference are the
    /// `Err` here.
    pub fn new(shape: &TopoShape, dir: &GpDir, angle: f64) -> Result<Self, String> {
        // `cxx:196-198`: every handle starts null.
        let mut this = Self {
            my_dir: *dir,
            my_angle: 0.0,
            angmin: 0.0,
            angmax: 0.0,
            my_tol: 1.0e-4,
            my_loc: None,
            my_sec: None,
            my_sections: ShapeTable::default(),
            my_faces: ShapeTable::default(),
            my_generated: Vec::new(),
            my_shape: None,
            my_top: shape.clone(),
            my_shell: None,
            my_wire: Wire::new(),
            my_cont: Shape::C1,
            my_style: TransitionStyle::Right,
            is_internal: false,
            my_done: false,
        };

        // `cxx:200-262`.
        this.my_wire = match shape.shape_type() {
            ShapeType::Wire => Wire::wrap(shape.clone()).ok_or(CONSTRUCTION_ERROR)?,
            ShapeType::Face => {
                // `TopoDS_Iterator Exp(S); myWire = TopoDS::Wire(Exp.Value())`
                // (`cxx:207-208`): the face's first boundary wire.
                let face = Face::wrap(shape.clone()).ok_or(CONSTRUCTION_ERROR)?;
                crate::topo_tools::wires_of_face(&face)
                    .into_iter()
                    .next()
                    .ok_or(CONSTRUCTION_ERROR)?
            }
            ShapeType::Shell => free_border_wire(shape)?,
            _ => return Err(CONSTRUCTION_ERROR.to_string()),
        };

        // `cxx:266-276`: "Attention to closed non declared wires".
        if !this.my_wire.0.is_closed() {
            let (vf, vl) = wire_vertices(&this.my_wire);
            if let (Some(vf), Some(vl)) = (vf, vl) {
                if is_same(&vf.0, &vl.0) {
                    this.my_wire.0.set_closed(true);
                }
            }
        }

        // `cxx:278-282`.
        this.my_angle = angle.abs();
        this.my_dir = *dir;
        this.my_top = shape.clone();
        this.my_done = false;
        this.my_tol = 1.0e-4;
        this.my_cont = Shape::C1;
        // `SetOptions()` with the `hxx:38-40` defaults.
        this.set_options(TransitionStyle::Right, 0.01, 3.0);
        this.set_draft(false);
        Ok(this)
    }

    /// `SetOptions(Style, Min, Max)` (`cxx:286-292`).
    pub fn set_options(&mut self, style: TransitionStyle, min: f64, max: f64) {
        self.my_style = style;
        self.angmin = min;
        self.angmax = max;
    }

    /// `SetDraft(Internal)` (`cxx:294-302`).
    pub fn set_draft(&mut self, internal: bool) {
        self.is_internal = internal;
    }

    /// `Perform(LengthMax)` (`cxx:306-317`): the skin alone, `LengthMax` long.
    pub fn perform(&mut self, length_max: f64) -> Result<(), String> {
        // `cxx:309-310`: the stopping surface stays null.
        let (w_box, _trsf) = geometry::compute_trsf(&self.my_wire, &self.my_dir)?;
        self.init(length_max, &w_box)?;
        self.build_shell(None, false)?;
        self.sewing()?;
        Ok(())
    }

    /// `Perform(Surface, KeepInsideSurface)` (`cxx:319-349`): the skin stopped
    /// by a surface.
    pub fn perform_on_surface(
        &mut self,
        surface: &Arc<dyn occt_geom::Surface>,
        keep_inside_surface: bool,
    ) -> Result<(), String> {
        // `cxx:331-336`.
        let (w_box, trsf) = geometry::compute_trsf(&self.my_wire, &self.my_dir)?;
        let surf: Arc<dyn occt_geom::Surface> = Arc::from(surface.transformed(&trsf));

        // `cxx:337-339`: the box of the stop surface.
        let mut s_box = BndBox::new();
        s_box.set_void();
        crate::geom_bnd_lib_surface3d::add_surface_with_bounds(&surf, 0.1, &mut s_box);

        // `cxx:342-343`.
        let mut pt = GpPnt::new(0.0, 0.0, 0.0);
        let mut len = geometry::longueur(&w_box, &s_box, &mut self.my_dir, &mut pt);
        len /= self.my_angle.cos().abs();

        // `cxx:346-348`: `Init` receives the given (untransformed) surface,
        // `BuildShell` the transformed one (`cxx:333-335`).
        self.init(len, &w_box)?;
        self.build_shell(Some(&surf), !keep_inside_surface)?;
        self.sewing()?;
        Ok(())
    }

    /// `Perform(StopShape, KeepOutSide)` (`cxx:351-407`): the skin stopped by a
    /// shape.
    pub fn perform_on_shape(
        &mut self,
        stop_shape: &TopoShape,
        keep_out_side: bool,
    ) -> Result<(), String> {
        // `cxx:363`.
        let (w_box, trsf) = geometry::compute_trsf(&self.my_wire, &self.my_dir)?;

        // `cxx:365-390`: the grouped box of every face of the stop shape.
        // `BSurf` is declared outside the loop and never reset, so
        // `SBox.Add(BSurf)` adds the accumulated box; adding every face box
        // into `SBox` is the same set of extents.
        let mut s_box = BndBox::new();
        s_box.set_void();
        let mut b_surf = BndBox::new();
        for face in crate::topo_tools::faces_of(stop_shape) {
            let (u_min, u_max, v_min, v_max) = crate::brep_tools::uv_bounds(&face);
            let Some(surface) = BRepTool::face_surface(&face) else {
                continue;
            };
            // `cxx:381-382`: `BRep_Tool::Surface(F)->Transformed(Trsf)`.
            let surf: Arc<dyn occt_geom::Surface> = Arc::from(surface.transformed(&trsf));
            crate::geom_bnd_lib_surface3d::add_surface(
                &surf, u_min, u_max, v_min, v_max, 0.1, &mut b_surf,
            );
            s_box.add_box(&b_surf);
        }

        // `cxx:392-393`.
        let mut pt = GpPnt::new(0.0, 0.0, 0.0);
        let mut len = geometry::longueur(&w_box, &s_box, &mut self.my_dir, &mut pt);
        len /= self.my_angle.cos().abs();

        // `cxx:397-400`: back to absolute space, then the stopping plane.
        let inv = trsf.inverted().map_err(|e| format!("gp_Trsf::Inverted: {e}"))?;
        pt.transform(&inv);
        let plan = Arc::new(occt_geom::plane::GeomPlane::new(GpPln::new(
            GpAx3::from_point_dir(pt, self.my_dir).map_err(|e| format!("gp_Ax3: {e}"))?,
        )));
        // `cxx:400`: `Geom_RectangularTrimmedSurface(Plan, -L, L, -L, L)`.
        let surf: Arc<dyn occt_geom::Surface> = Arc::new(
            GeomRectangularTrimmedSurface::uv(plan.clone(), -len, len, -len, len)
                .map_err(|e| format!("Geom_RectangularTrimmedSurface: {e}"))?,
        );

        // `cxx:403-406`.
        self.init(len * 1.01, &w_box)?;
        self.build_shell(Some(&surf), false)?;
        self.fuse(stop_shape, keep_out_side)?;
        self.sewing()?;
        Ok(())
    }

    /// `IsDone()` (`cxx:970-974`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `Shell()` (`cxx:963-968`): a null shell is `None`.
    pub fn shell(&self) -> Option<&Shell> {
        self.my_shell.as_ref()
    }

    /// `Shape()` (`cxx:954-958`): a null shape is `None`.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `Generated(S)` (`cxx:922-952`): the shape the sweep made from the edge
    /// `S`, or an empty list when `S` is not one of the path's edges.
    ///
    /// The reference's vertex arm (`cxx:929-938`) is unreachable — see the
    /// module note — so a shape that is not an edge is the
    /// `Standard_TypeMismatch` its cast raises.
    pub fn generated(&mut self, s: &TopoShape) -> Result<&[TopoShape], String> {
        self.my_generated.clear();
        // `cxx:928`: `E = TopoDS::Edge(S)`.
        let e = crate::shape::Edge::wrap(s.clone())
            .ok_or_else(|| "Standard_TypeMismatch: TopoDS::Edge".to_string())?;
        // `cxx:941-948`.
        let loc = self
            .my_loc
            .as_ref()
            .ok_or_else(|| "Standard_NullObject: BRepFill_Draft::myLoc".to_string())?;
        let nb_law = loc.base().nb_law();
        for ii in 1..=nb_law {
            let is_path_edge = match loc.base().edge(ii) {
                Some(lei) => is_same(&e.0, &lei.0),
                None => false,
            };
            if is_path_edge {
                let face = self.my_faces.value(1, ii).cloned();
                if let Some(f) = face {
                    self.my_generated.push(f);
                }
                break;
            }
        }
        Ok(&self.my_generated)
    }

    /// `myLoc` (`hxx:87`), `None` before `Init` ran.
    pub fn location_law(&self) -> Option<&Arc<DraftLaw>> {
        self.my_loc.as_ref()
    }

    /// `IsInternal` (`hxx:98`).
    pub fn is_internal(&self) -> bool {
        self.is_internal
    }

    /// `NbLaw()` on the stored location law — the `myFaces` / `mySections`
    /// column count `Fuse` and `Sewing` walk.
    fn nb_law(&self) -> usize {
        self.my_loc.as_ref().map(|l| l.base().nb_law()).unwrap_or(0)
    }
}

/// The `TopoDS_SHELL` arm of the constructor (`cxx:211-258`): the shell's free
/// borders — every non-degenerated edge with exactly one face ancestor —
/// assembled into one wire.
///
/// `TopExp::MapShapesAndAncestors` builds an `NCollection_IndexedDataMap`, and
/// the reference walks it by `FindKey(i)`, i.e. in insertion order. The port
/// keeps that order explicitly instead of walking a `HashMap`.
fn free_border_wire(shape: &TopoShape) -> Result<Wire, String> {
    let mut order: Vec<ShapeMapKey> = Vec::new();
    let mut edges: std::collections::HashMap<ShapeMapKey, TopoShape> =
        std::collections::HashMap::new();
    let mut ancestors: std::collections::HashMap<ShapeMapKey, Vec<TopoShape>> =
        std::collections::HashMap::new();

    // `cxx:217`: the `(EDGE, FACE)` pass of `MapShapesAndAncestors`.
    for face in crate::topo_tools::explorer_shapes(shape, ShapeType::Face) {
        for e in crate::topo_tools::explorer_shapes(&face, ShapeType::Edge) {
            let k = shape_map_key(&e);
            if !edges.contains_key(&k) {
                order.push(k.clone());
                edges.insert(k.clone(), e.clone());
            }
            ancestors.entry(k).or_default().push(face.clone());
        }
    }
    // `cxx:110-119`: the trailing `(EDGE, FACE)` pass, which only creates
    // empty entries (their faces count is 0, so they never enter the list).
    for e in crate::topo_tools::explorer_shapes_avoiding(shape, ShapeType::Edge, ShapeType::Face) {
        let k = shape_map_key(&e);
        if !edges.contains_key(&k) {
            order.push(k.clone());
            edges.insert(k.clone(), e);
            ancestors.insert(k, Vec::new());
        }
    }

    // `cxx:220-231`.
    let mut list: Vec<crate::shape::Edge> = Vec::new();
    for k in order {
        let (Some(edge), Some(faces)) = (edges.get(&k), ancestors.get(&k)) else {
            continue;
        };
        let Some(edge) = crate::shape::Edge::wrap(edge.clone()) else {
            return Err(CONSTRUCTION_ERROR.to_string());
        };
        // `cxx:224`: skip degenerated edges; `cxx:226-229`: keep the free ones.
        if !BRepTool::is_degenerated(&edge) && faces.len() == 1 {
            list.push(edge);
        }
    }

    // `cxx:234-257`.
    if list.is_empty() {
        return Err(CONSTRUCTION_ERROR.to_string());
    }
    let mut mw = crate::brep_lib_make_wire::MakeWire::new();
    for e in &list {
        mw.add(e).map_err(|_| CONSTRUCTION_ERROR.to_string())?;
    }
    if mw.error() != crate::brep_lib_make_wire::WireError::WireDone {
        return Err(CONSTRUCTION_ERROR.to_string());
    }
    Ok(mw.wire())
}
