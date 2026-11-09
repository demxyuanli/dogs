//! `Draft_Modification` - the draft-angle `BRepTools_Modification`
//! (`Draft_Modification.cxx:54-510`) plus the private overloads of
//! `NewSurface` / `NewCurve` and the propagation that fills the records
//! (`Draft_Modification_1.cxx:116-708`, `:1780-2044`, `:2403-2453`).
//!
//! | OCCT | file | status |
//! |---|---|---|
//! | constructor / `Clear` / `Init` / `Add` / `Remove` / `IsDone` / `Error` / `ProblematicShape` / `ConnectedFaces` / `ModifiedFaces` (`cxx:56-224`) | `mod.rs` | ported |
//! | `NewSurface(F, ...)` (`cxx:228-256`) | `mod.rs` | ported |
//! | `NewCurve(E, ...)` (`cxx:261-286`) | `mod.rs` | ported |
//! | `NewPoint` (`cxx:291-305`) | `mod.rs` | ported |
//! | `NewCurve2d` (`cxx:310-427`) | `mod.rs` | ported |
//! | `NewParameter` (`cxx:432-495`) | `mod.rs` | ported |
//! | `Continuity` (`cxx:502-510`) | `mod.rs` | ported |
//! | `InternalAdd` (`_1.cxx:116-420`) | `internal_add.rs` | ported |
//! | `Propagate` (`_1.cxx:424-705`) | `propagate.rs` | ported |
//! | private `NewSurface` + `FindRotation` (`_1.cxx:1780-1960`, `:2403-2453`) | `new_surface.rs` | ported |
//! | private `NewCurve` (`_1.cxx:1963-2044`) | `new_curve.rs` | ported |
//! | `Perform` (`_1.cxx:708-1777`) | `perform/` | ported |
//! | `Choose` / `Parameter` / `SmartParameter` (`_1.cxx:2047-2380`) | `perform/` | ported |
//!
//! `Perform` is the entry point: it drives the phase order of the reference
//! (`_1.cxx:708-725`) and every phase lives in its own file under `perform/`.
//!
//! # Translation notes
//!
//! * The three `NCollection_IndexedDataMap` members are [`ShapeDataMap`], an
//!   insertion-ordered 1-based map keyed by [`shape_map_key`] - the
//!   `TopTools_ShapeMapHasher` equality `TopoDS_Shape::IsSame` tests, so
//!   orientation is ignored (`hxx:203-207`).
//! * `badShape` is `Option<TopoShape>`: `TopoShape` has no null state in this
//!   port, so `IsNull()` is `is_none()` and `Nullify()` is `None`.
//! * `Draft_Modification` inherits `BRepTools_Modification`; the Rust port
//!   keeps the same six call shapes (`NewSurface`, `NewCurve`, `NewPoint`,
//!   `NewCurve2d`, `NewParameter`, `Continuity`) as inherent methods that
//!   answer `Result<bool, String>`, the `Err` being the reference's
//!   `Standard_DomainError` / `Standard_NoSuchObject` / `StdFail_NotDone`
//!   throws.
//! * `NewCurve2d` needs the "recadre" block (`cxx:368-422`) because the port's
//!   `Geom_TrimmedCurve` reparameterizes the trim onto `[0, 1]`
//!   (`trimmed.rs:17-38`) where OCCT keeps the basis parameters; the projection
//!   is therefore run over the trimmed curve's own range, which is the same
//!   geometric domain as `cxx:362-365`.

mod internal_add;
mod new_curve;
mod new_surface;
mod perform;
mod propagate;
mod shape_data_map;

use std::sync::Arc;

use occt_core::gp::{GpDir, GpPln, GpTrsf2d, GpVec2d};
use occt_core::precision::PCONFUSION;
use occt_geom::curve::Curve;
use occt_geom::geom_proj_lib;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_lib_continuity::stored_continuity;
use crate::brep_tool::BRepTool;
use crate::brep_tools::eval_and_update_tol;
use crate::builder::TopoBuilder;
use crate::draft::edge_info::DraftEdgeInfo;
use crate::draft::error_status::DraftErrorStatus;
use crate::draft::face_info::DraftFaceInfo;
use crate::draft::vertex_info::DraftVertexInfo;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, explorer_shapes, is_equal, is_same, AncestorMap};

use shape_data_map::ShapeDataMap;

/// `class Draft_Modification` (`Draft_Modification.hxx:47-217`).
pub struct DraftModification {
    /// `myFMap` (`hxx:203`).
    my_f_map: ShapeDataMap<DraftFaceInfo>,
    /// `myEMap` (`hxx:204`).
    my_e_map: ShapeDataMap<DraftEdgeInfo>,
    /// `myVMap` (`hxx:205`).
    my_v_map: ShapeDataMap<DraftVertexInfo>,
    /// `myComp` (`hxx:208`).
    my_comp: bool,
    /// `myShape` (`hxx:209`).
    my_shape: TopoShape,
    /// `badShape` (`hxx:210`).
    bad_shape: Option<TopoShape>,
    /// `errStat` (`hxx:211`).
    err_stat: DraftErrorStatus,
    /// `curFace` (`hxx:212`).
    cur_face: Face,
    /// `conneF` (`hxx:213`).
    conne_f: Vec<TopoShape>,
    /// `myEFMap` (`hxx:214-215`).
    my_ef_map: AncestorMap,
}

impl DraftModification {
    /// `Draft_Modification(S)` (`cxx:54-61`).
    pub fn new(s: &TopoShape) -> Self {
        Self {
            my_f_map: ShapeDataMap::new(),
            my_e_map: ShapeDataMap::new(),
            my_v_map: ShapeDataMap::new(),
            my_comp: false,
            my_shape: s.clone(),
            bad_shape: None,
            err_stat: DraftErrorStatus::NoError,
            cur_face: Face::new(),
            conne_f: Vec::new(),
            my_ef_map: AncestorMap::map_shapes_and_ancestors(s, ShapeType::Edge, ShapeType::Face),
        }
    }

    /// `Clear()` (`cxx:64-74`).
    pub fn clear(&mut self) {
        self.my_comp = false;
        self.my_f_map.clear();
        self.my_e_map.clear();
        self.my_v_map.clear();
        self.my_ef_map = AncestorMap::default();
        self.bad_shape = None;
        self.err_stat = DraftErrorStatus::NoError;
    }

    /// `Init(S)` (`cxx:78-83`).
    pub fn init(&mut self, s: &TopoShape) {
        self.my_shape = s.clone();
        self.clear();
        self.my_ef_map =
            AncestorMap::map_shapes_and_ancestors(s, ShapeType::Edge, ShapeType::Face);
    }

    /// `Add(F, Direction, Angle, NeutralPlane, Flag = true)` (`cxx:87-105`).
    ///
    /// `Err` is the reference's `Standard_ConstructionError` thrown when an
    /// error is already pending (`cxx:94-97`).
    pub fn add(
        &mut self,
        f: &Face,
        direction: &GpDir,
        angle: f64,
        neutral_plane: &GpPln,
        flag: bool,
    ) -> Result<bool, String> {
        // `cxx:94-97`
        if self.bad_shape.is_some() {
            return Err("Standard_ConstructionError: Draft_Modification::Add".to_string());
        }
        // `cxx:98-101`: `Clear()` also empties `myEFMap` (`cxx:71`), so a
        // second `Add` after `Perform` runs `InternalAdd` with no ancestor
        // map, exactly like the reference (whose `FindFromKey` then throws).
        if self.my_comp {
            self.clear();
        }
        // `cxx:102`
        self.cur_face = f.clone();
        Ok(self.internal_add(f, direction, angle, neutral_plane, flag))
    }

    /// `Remove(F)` (`cxx:108-153`).
    pub fn remove(&mut self, f: &Face) -> Result<(), String> {
        // `cxx:110-113`
        if !self.my_f_map.contains(&f.0) || self.my_comp {
            return Err("Standard_NoSuchObject: Draft_Modification::Remove".to_string());
        }
        // `cxx:117-132`
        self.conne_f.clear();
        self.cur_face = self
            .my_f_map
            .find_from_key(&f.0)
            .and_then(|fi| fi.root_face_value())
            .cloned()
            .unwrap_or_else(Face::new);
        for i in 1..=self.my_f_map.extent() {
            let the_f = self.my_f_map.find_key(i).clone();
            let same_root = self
                .my_f_map
                .find_from_key(&the_f)
                .and_then(|fi| fi.root_face_value())
                .map(|r| is_same(&r.0, &self.cur_face.0))
                .unwrap_or(false);
            if same_root {
                self.conne_f.push(the_f.clone());
                if self.bad_shape.as_ref().map(|b| is_same(b, &the_f)).unwrap_or(false) {
                    self.bad_shape = None;
                }
            }
        }
        // `cxx:133-137`
        for s in self.conne_f.clone() {
            self.my_f_map.remove_key(&s);
        }
        // `cxx:139-153`
        self.conne_f.clear();
        for i in 1..=self.my_e_map.extent() {
            let the_e = self.my_e_map.find_key(i).clone();
            let same_root = self
                .my_e_map
                .find_from_key(&the_e)
                .and_then(|ei| ei.root_face_value())
                .map(|r| is_same(&r.0, &self.cur_face.0))
                .unwrap_or(false);
            if same_root {
                self.conne_f.push(the_e);
            }
        }
        for s in self.conne_f.clone() {
            self.my_e_map.remove_key(&s);
        }
        Ok(())
    }

    /// `IsDone()` (`cxx:158-161`).
    pub fn is_done(&self) -> bool {
        self.my_comp && self.bad_shape.is_none()
    }

    /// `Error()` (`cxx:165-168`).
    pub fn error(&self) -> DraftErrorStatus {
        self.err_stat
    }

    /// `ProblematicShape()` (`cxx:172-175`).
    pub fn problematic_shape(&self) -> Option<&TopoShape> {
        self.bad_shape.as_ref()
    }

    /// `ConnectedFaces(F)` (`cxx:179-202`): every face sharing `F`'s root face.
    ///
    /// `Err` carries the reference's `Standard_NoSuchObject` /
    /// `StdFail_NotDone` throws.
    pub fn connected_faces(&mut self, f: &Face) -> Result<&[TopoShape], String> {
        // `cxx:181-184`
        if !self.my_f_map.contains(&f.0) {
            return Err("Standard_NoSuchObject: Draft_Modification::ConnectedFaces".to_string());
        }
        // `cxx:185-188`
        if !self.is_done() {
            return Err("StdFail_NotDone: Draft_Modification::ConnectedFaces".to_string());
        }
        // `cxx:189-199`
        self.conne_f.clear();
        self.cur_face = self
            .my_f_map
            .find_from_key(&f.0)
            .and_then(|fi| fi.root_face_value())
            .cloned()
            .unwrap_or_else(Face::new);
        for i in 1..=self.my_f_map.extent() {
            let the_f = self.my_f_map.find_key(i).clone();
            let same_root = self
                .my_f_map
                .find_from_key(&the_f)
                .and_then(|fi| fi.root_face_value())
                .map(|r| is_same(&r.0, &self.cur_face.0))
                .unwrap_or(false);
            if same_root {
                self.conne_f.push(the_f);
            }
        }
        Ok(&self.conne_f)
    }

    /// `ModifiedFaces()` (`cxx:206-224`): every face with a root face.
    pub fn modified_faces(&mut self) -> Result<&[TopoShape], String> {
        // `cxx:208-212`
        if self.bad_shape.is_some() {
            return Err("StdFail_NotDone: Draft_Modification::ModifiedFaces".to_string());
        }
        // `cxx:213-222`
        self.conne_f.clear();
        for i in 1..=self.my_f_map.extent() {
            let the_f = self.my_f_map.find_key(i).clone();
            if self
                .my_f_map
                .find_from_key(&the_f)
                .and_then(|fi| fi.root_face_value())
                .is_some()
            {
                self.conne_f.push(the_f);
            }
        }
        Ok(&self.conne_f)
    }

    /// `NewSurface(F, S, L, Tol, RevWires, RevFace)`
    /// (`Draft_Modification.cxx:228-256`).
    pub fn new_surface_of(
        &self,
        f: &Face,
        s: &mut Option<Arc<dyn Surface>>,
        tol: &mut f64,
        rev_wires: &mut bool,
        rev_face: &mut bool,
    ) -> Result<bool, String> {
        // `cxx:233-236`
        if !self.is_done() {
            return Err("Standard_DomainError: Draft_Modification::NewSurface".to_string());
        }
        // `cxx:240-243`
        let Some(fi) = self.my_f_map.find_from_key(&f.0) else {
            return Ok(false);
        };
        if !fi.new_geometry() {
            return Ok(false);
        }
        // `cxx:245-253`
        *rev_wires = false;
        *rev_face = false;
        *tol = BRepTool::face_tolerance(f);
        // `cxx:249-252`: `S = BRep_Tool::Surface(F, L); L.Identity();
        // S = myFMap.FindFromKey(F).Geometry();` - the location is reset, the
        // geometry of the record replaces the looked-up surface.
        *s = fi.geometry().cloned();
        Ok(true)
    }

    /// `NewCurve(E, C, L, Tol)` (`Draft_Modification.cxx:260-287`).
    pub fn new_curve_of(
        &self,
        e: &Edge,
        c: &mut Option<Arc<dyn Curve>>,
        tol: &mut f64,
    ) -> Result<bool, String> {
        // `cxx:265-268`
        if !self.is_done() {
            return Err("Standard_DomainError: Draft_Modification::NewCurve".to_string());
        }
        // `cxx:270-273`
        let Some(einf) = self.my_e_map.find_from_key(&e.0) else {
            return Ok(false);
        };
        // `cxx:276-279`
        if !einf.new_geometry() {
            return Ok(false);
        }
        // `cxx:281-285`
        *tol = einf.tolerance().max(BRepTool::edge_tolerance(e));
        *c = einf.geometry().cloned();
        Ok(true)
    }

    /// `NewPoint(V, P, Tol)` (`Draft_Modification.cxx:291-306`).
    pub fn new_point_of(
        &self,
        v: &Vertex,
        p: &mut occt_core::gp::GpPnt,
        tol: &mut f64,
    ) -> Result<bool, String> {
        // `cxx:295-298`
        if !self.is_done() {
            return Err("Standard_DomainError: Draft_Modification::NewPoint".to_string());
        }
        // `cxx:299-302`
        let Some(vinf) = self.my_v_map.find_from_key(&v.0) else {
            return Ok(false);
        };
        // `cxx:304-305`
        *tol = BRepTool::vertex_tolerance(v);
        *p = *vinf.geometry();
        Ok(true)
    }

    /// `NewCurve2d(E, F, NewE, NewF, C, Tol)`
    /// (`Draft_Modification.cxx:310-427`).
    pub fn new_curve2d_of(
        &self,
        e: &Edge,
        f: &Face,
        new_e: &Edge,
        c: &mut Option<Arc<dyn Curve2d>>,
        tol: &mut f64,
    ) -> Result<bool, String> {
        // `cxx:319-322`
        if !self.is_done() {
            return Err("Standard_DomainError: Draft_Modification::NewCurve2d".to_string());
        }
        // `cxx:324-327`
        if !self.my_e_map.contains(&e.0) {
            return Ok(false);
        }
        // `cxx:329-330`
        let (mut fp, mut lp) = BRepTool::edge_parameters(new_e);
        // `cxx:332`: `myFMap.FindFromKey(F).Geometry()`.
        let Some(sb) = self.my_f_map.find_from_key(&f.0).and_then(|fi| fi.geometry()).cloned()
        else {
            return Err("Standard_NoSuchObject: Draft_Modification::NewCurve2d".to_string());
        };
        // `cxx:334`
        *tol = BRepTool::edge_tolerance(e);

        // `cxx:336-345`: the recorded pcurves win.
        let einf = self.my_e_map.find_from_key(&e.0).unwrap();
        let first_hit = einf
            .first_face()
            .map(|ff| is_same(&ff.0, &f.0))
            .unwrap_or(false);
        let second_hit = einf
            .second_face()
            .map(|sf| is_same(&sf.0, &f.0))
            .unwrap_or(false);
        let recorded = if first_hit {
            einf.first_pc().cloned()
        } else if second_hit {
            einf.second_pc().cloned()
        } else {
            None
        };
        let mut new_c2d = recorded;
        if new_c2d.is_none() {
            // `cxx:346-355`
            if !einf.new_geometry() {
                let (fpi, lpi) = BRepTool::edge_parameters(e);
                if fpi <= fp && fp <= lpi && fpi <= lp && lp <= lpi {
                    return Ok(false);
                }
            }
            // `cxx:357-365`: `Geom_TrimmedCurve` over the new edge's range, then
            // the projection on that same (now ordered) range.
            let (a, b) = BRepTool::edge_parameters(new_e);
            let (nfirst, nlast) = (a.min(b), a.max(b));
            let Some(geom) = einf.geometry().cloned() else {
                return Err("Standard_NoSuchObject: Draft_Modification::NewCurve2d".to_string());
            };
            let tc: Arc<dyn Curve> =
                Arc::new(occt_geom::trimmed::GeomTrimmedCurve::new(geom, nfirst, nlast));
            fp = tc.first_parameter();
            lp = tc.last_parameter();
            GeometryRegistry::global().set_edge_range(&new_e.0, fp, lp);
            let mut ptol = *tol;
            let Some(pc) = geom_proj_lib::curve_2d_range(&tc, fp, lp, &sb, &mut ptol) else {
                return Ok(false);
            };
            *tol = ptol;
            new_c2d = Some(pc);
        }

        // `cxx:368-388`: is the support one the 2*pi shift has to reconcile?
        let mut sbasis = sb.clone();
        let mut typs = surface_kind(&sbasis);
        if typs == SurfaceKind::RectangularTrimmed {
            sbasis = sbasis.rectangular_trimmed_basis().unwrap();
            typs = surface_kind(&sbasis);
        }
        let mut je_recadre = false;
        if typs == SurfaceKind::LinearExtrusion {
            if let Some(a_c) = sbasis.extrusion_basis_curve() {
                if a_c.gp_circ().is_some() {
                    je_recadre = true;
                }
            }
        }
        je_recadre = je_recadre
            || typs == SurfaceKind::Cylinder
            || typs == SurfaceKind::Sphere
            || typs == SurfaceKind::Cone;

        // `cxx:390-422`
        if je_recadre {
            let mut b_translate = false;
            let mut a_v2dt = GpVec2d::new(0.0, 0.0);
            let vectra = GpVec2d::new(2.0 * std::f64::consts::PI, 0.0);
            let existing = crate::boptools_2d::curve_on_surface_oriented(e, f, false)
                .map(|(pc, t1, t2)| (pc, t1, t2));
            let (pf, confidence) = match existing {
                Some((pc, t1, t2)) => (pc.d0(0.5 * (t1 + t2)), pc),
                None => {
                    return Err(
                        "Standard_NoSuchObject: Draft_Modification::NewCurve2d".to_string()
                    )
                }
            };
            let _ = confidence;
            let new_pf = new_c2d.as_ref().unwrap().d0(0.5 * (fp + lp));
            let a_d2 = new_pf.square_distance(&pf);
            if new_pf.translated_vec(&vectra).square_distance(&pf) < a_d2 {
                a_v2dt = vectra;
                b_translate = !b_translate;
            } else if new_pf.translated_vec(&vectra.multiplied_scalar(-1.0)).square_distance(&pf) < a_d2
            {
                a_v2dt = vectra.multiplied_scalar(-1.0);
                b_translate = !b_translate;
            }
            if b_translate {
                let mut trsf = GpTrsf2d::identity();
                trsf.set_translation_vec(&a_v2dt);
                let mut shifted: Box<dyn Curve2d> = new_c2d.as_ref().unwrap().clone_dyn();
                shifted.transform(&trsf);
                new_c2d = Some(Arc::from(shifted));
            }
        }

        // `cxx:423-424`
        let Some(c3d) = BRepTool::edge_curve(new_e) else {
            return Ok(false);
        };
        let (rf, rl) = BRepTool::edge_parameters(new_e);
        *tol = eval_and_update_tol(new_e, &c3d, new_c2d.as_ref().unwrap(), &sbasis, rf, rl);
        *c = new_c2d;
        Ok(true)
    }

    /// `NewParameter(V, E, P, Tol)` (`Draft_Modification.cxx:431-498`).
    pub fn new_parameter_of(
        &self,
        v: &Vertex,
        e: &Edge,
        p: &mut f64,
        tol: &mut f64,
    ) -> Result<bool, String> {
        // `cxx:438-441`
        if !self.is_done() {
            return Err("Standard_DomainError: Draft_Modification::NewParameter".to_string());
        }
        // `cxx:443-446`
        if !self.my_v_map.contains(&v.0) {
            return Ok(false);
        }
        // `cxx:448-449`
        *p = self
            .my_v_map
            .find_from_key(&v.0)
            .map(|vi| vi.parameter(e))
            .unwrap_or(0.0);
        // `cxx:450-454`
        let Some(mut gc) = self.my_e_map.find_from_key(&e.0).and_then(|ei| ei.geometry()).cloned()
        else {
            return Err("Standard_NoSuchObject: Draft_Modification::NewParameter".to_string());
        };
        if gc.is_geom_trimmed() {
            if let Some(basis) = gc.trim_basis() {
                gc = basis;
            }
        }
        // `cxx:456-492`
        if gc.is_closed() {
            let (fv, _lv) = edge_vertices(e);
            let mut paramf = match &fv {
                Some(fv) => {
                    let has = self.my_v_map.contains(&fv.0);
                    if has {
                        self.my_v_map
                            .find_from_key(&fv.0)
                            .map(|vi| vi.parameter(e))
                            .unwrap_or(0.0)
                    } else {
                        BRepTool::parameter_of_vertex(&e.0, &fv.0)
                    }
                }
                None => 0.0,
            };
            // `cxx:469-487`
            let first_par = gc.first_parameter();
            let last_par = gc.last_parameter();
            // `cxx:476-483`
            if let Some(fv) = &fv {
                if (paramf - last_par).abs() <= PCONFUSION {
                    paramf = first_par;
                    let fv_oriented = fv.0.oriented(e.0.orientation());
                    if is_equal(&v.0, &fv_oriented) {
                        *p = paramf;
                    }
                }
                // `cxx:484-491`
                let fv_oriented = fv.0.oriented(e.0.orientation());
                if !is_equal(&v.0, &fv_oriented) && *p <= paramf {
                    if gc.is_periodic() {
                        *p += gc.period();
                    } else {
                        *p = gc.last_parameter();
                    }
                }
            }
        }
        // `cxx:494`
        *tol = BRepTool::vertex_tolerance(v).max(BRepTool::edge_tolerance(e));
        Ok(true)
    }

    /// `Continuity(E, F1, F2, NewE, NewF1, NewF2)` (`cxx:502-510`).
    pub fn continuity_of(&self, e: &Edge, f1: &Face, f2: &Face) -> u8 {
        stored_continuity(e, f1, f2)
    }
}

/// `Draft_Modification : BRepTools_Modification` (`Draft_Modification.hxx:47`).
///
/// The six overloads are the ones `BRepTools_Modifier::Perform` calls while it
/// rebuilds the shape; they delegate to the inherent out-parameter methods.
/// The inherent methods answer `Err` only for the reference's
/// `Standard_DomainError` / `Standard_NoSuchObject` throws, and OCCT's modifier
/// would let those propagate out of `Perform`; `DoModif` only runs after
/// `IsDone()` is true, so the `Err` arm is unreachable from
/// `BRepOffsetAPI_DraftAngle::Build` and is reported as "not modified" here.
impl crate::brep_tools_modifier::Modification for DraftModification {
    /// `NewSurface` (`Draft_Modification.cxx:228-256`).
    fn new_surface(&self, f: &TopoShape) -> Option<crate::brep_tools_modifier::NewSurfaceInfo> {
        let face = Face::wrap(f.clone())?;
        let mut surface = None;
        let mut tol = 0.0;
        let (mut rev_wires, mut rev_face) = (false, false);
        let ok = self
            .new_surface_of(&face, &mut surface, &mut tol, &mut rev_wires, &mut rev_face)
            .ok()?;
        if !ok {
            return None;
        }
        Some(crate::brep_tools_modifier::NewSurfaceInfo {
            surface: surface?,
            // `L.Identity()` (`cxx:250`): the record replaces the looked-up
            // surface, so the location is the identity.
            location: occt_core::toploc::TopLocLocation::identity(),
            tolerance: tol,
            rev_wires,
            rev_face,
        })
    }

    /// `NewCurve` (`Draft_Modification.cxx:260-287`).
    fn new_curve(&self, e: &TopoShape) -> Option<crate::brep_tools_modifier::NewCurveInfo> {
        let edge = Edge::wrap(e.clone())?;
        let mut curve = None;
        let mut tol = 0.0;
        let ok = self.new_curve_of(&edge, &mut curve, &mut tol).ok()?;
        if !ok {
            return None;
        }
        Some(crate::brep_tools_modifier::NewCurveInfo {
            curve,
            // `L.Identity()` (`cxx:284`).
            location: occt_core::toploc::TopLocLocation::identity(),
            tolerance: tol,
        })
    }

    /// `NewPoint` (`Draft_Modification.cxx:291-306`).
    fn new_point(&self, v: &TopoShape) -> Option<crate::brep_tools_modifier::NewPointInfo> {
        let vertex = Vertex::wrap(v.clone())?;
        let mut p = occt_core::gp::GpPnt::new(0.0, 0.0, 0.0);
        let mut tol = 0.0;
        let ok = self.new_point_of(&vertex, &mut p, &mut tol).ok()?;
        if !ok {
            return None;
        }
        Some(crate::brep_tools_modifier::NewPointInfo {
            point: p,
            tolerance: tol,
        })
    }

    /// `NewCurve2d` (`Draft_Modification.cxx:310-427`).
    fn new_curve2d(
        &self,
        e: &TopoShape,
        f: &TopoShape,
        new_e: &TopoShape,
        _new_f: &TopoShape,
    ) -> Option<Arc<dyn Curve2d>> {
        let (edge, face, new_edge) = (
            Edge::wrap(e.clone())?,
            Face::wrap(f.clone())?,
            Edge::wrap(new_e.clone())?,
        );
        let mut curve = None;
        let mut tol = 0.0;
        let ok = self
            .new_curve2d_of(&edge, &face, &new_edge, &mut curve, &mut tol)
            .ok()?;
        if !ok {
            return None;
        }
        curve
    }

    /// `NewParameter` (`Draft_Modification.cxx:431-498`).
    fn new_parameter(
        &self,
        v: &TopoShape,
        e: &TopoShape,
    ) -> Option<crate::brep_tools_modifier::NewParameterInfo> {
        let (vertex, edge) = (Vertex::wrap(v.clone())?, Edge::wrap(e.clone())?);
        let mut p = 0.0;
        let mut tol = 0.0;
        let ok = self.new_parameter_of(&vertex, &edge, &mut p, &mut tol).ok()?;
        if !ok {
            return None;
        }
        Some(crate::brep_tools_modifier::NewParameterInfo {
            parameter: p,
            tolerance: tol,
        })
    }

    /// `Continuity(E, F1, F2, NewE, NewF1, NewF2)` (`cxx:502-510`).
    fn continuity(
        &self,
        e: &TopoShape,
        f1: &TopoShape,
        f2: &TopoShape,
        _new_e: &TopoShape,
        _new_f1: &TopoShape,
        _new_f2: &TopoShape,
    ) -> occt_core::kernel::geomabs::Shape {
        let (Some(edge), Some(face1), Some(face2)) = (
            Edge::wrap(e.clone()),
            Face::wrap(f1.clone()),
            Face::wrap(f2.clone()),
        ) else {
            return occt_core::kernel::geomabs::Shape::C0;
        };
        crate::brep_tools_modifier::geom_abs_shape_from_u8(self.continuity_of(&edge, &face1, &face2))
    }
}

/// The surface categories `NewCurve2d` discriminates on.
#[derive(PartialEq, Eq, Clone, Copy)]
enum SurfaceKind {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    LinearExtrusion,
    RectangularTrimmed,
    Other,
}

fn surface_kind(s: &Arc<dyn Surface>) -> SurfaceKind {
    if s.gp_pln().is_some() {
        SurfaceKind::Plane
    } else if s.gp_cylinder().is_some() {
        SurfaceKind::Cylinder
    } else if s.gp_cone().is_some() {
        SurfaceKind::Cone
    } else if s.gp_sphere().is_some() {
        SurfaceKind::Sphere
    } else if s.is_surface_of_linear_extrusion() {
        SurfaceKind::LinearExtrusion
    } else if s.rectangular_trimmed_basis().is_some() {
        SurfaceKind::RectangularTrimmed
    } else {
        SurfaceKind::Other
    }
}

/// `static TopAbs_Orientation Orientation(S, F)`
/// (`Draft_Modification_1.cxx:2383-2398`): the occurrence orientation of the
/// first face of `s` that is the same as `f`, `FORWARD` when there is none.
pub(super) fn shape_face_orientation(s: &TopoShape, f: &Face) -> Orientation {
    for face in explorer_shapes(s, ShapeType::Face) {
        if is_same(&face, &f.0) {
            return face.orientation();
        }
    }
    Orientation::Forward
}

/// The builder `internal_add` uses for the free-border face of `Propagate`
/// (`_1.cxx:688-692`: `BRep_Builder B; B.MakeFace(TheNewFace, S2,
/// Precision::Confusion())`).
pub(super) fn make_face(surface: Arc<dyn Surface>, tolerance: f64) -> Face {
    TopoBuilder::new().make_face_with_tolerance(surface, &[], tolerance)
}
