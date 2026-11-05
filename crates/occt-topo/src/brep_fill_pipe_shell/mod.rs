//! `BRepFill_PipeShell` — a swept shell from a spine, a set of sections and a
//! trihedron mode (`BRepFill_PipeShell.cxx:1-1681`, `.hxx:47-261`).
//!
//! Port status, keeping the OCCT control flow:
//!
//! * Ported: the data layer, the constructor (`cxx:223-252`), `Set(Frenet)`
//!   (`cxx:256-273`), `SetDiscrete` (`cxx:279-289`), `Set(gp_Ax2)`
//!   (`cxx:292-304`), `Set(gp_Dir)` (`cxx:309-317`), `Set(SpineSupport)`
//!   (`cxx:323-337`), `Set(AuxiliarySpine, CurvilinearEquivalence,
//!   KeepContact)` (`cxx:343-431`), `SetMaxDegree` / `SetMaxSegments` /
//!   `SetForceApproxC1`, both `Add` overloads (`cxx:462-548`, including the
//!   automatic-homothety law `cxx:479-540`), both `SetLaw` overloads
//!   (`cxx:554-582`), `DeleteProfile` (`cxx:586-606`), `IsReady`, `GetStatus`,
//!   `SetTolerance`, `SetTransition` (`cxx:637-649`), `Simulate`
//!   (`cxx:651-701`), `Build` (`cxx:706-846`), `MakeSolid` (`cxx:849-914`), the
//!   five accessors, `Generated` (`cxx:946-956`), `Prepare`
//!   (`cxx:958-1236`, including the multi-section branch and the contact
//!   retiming), `Place` (`cxx:1242-1258`), `ResetLoc` (`cxx:1264-1275`),
//!   `BuildHistory` (`cxx:1282-1615`), `BuildBoundaries` (`cxx:1618-1681`) and
//!   the statics `PerformTransition` (`cxx:130-148`), `PerformPlan`
//!   (`cxx:152-183`) and `IsSameOriented` (`cxx:187-220`).
//! * `ResetLoc` and the contact retiming both need
//!   `occ::down_cast<GeomFill_LocationGuide>`, which the port spells
//!   `LocationLaw::as_location_guide` (no RTTI on trait objects).
//!
//! Translation notes:
//! - `handle(BRepFill_SectionLaw)` / `handle(BRepFill_LocationLaw)` are
//!   `Arc<dyn SectionLaw>` / `Arc<dyn LocationLaw>`; a null handle is `None`.
//! - A null `TopoDS_Shape` is `None`. `myShape` after a failed `Build` is the
//!   empty shell OCCT stores, so it is `Some` and not null.
//! - `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>` is a
//!   [`HashMap`] keyed by [`ShapeMapKey`].

pub mod build;
pub mod history;
pub mod prepare;
pub mod simulate;

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::{GpAx2, GpDir, GpPnt, GpVec};

use occt_geom::curve::Curve;
use occt_geom::geom_fill::plane_curve_inter::plane_intersections;
use occt_geom::geom_fill::{
    ConstantBiNormal, CorrectedFrenet, CurveAndTrihedron, DiscreteTrihedron, Fixed, Frenet,
    GeomFillPipeError, GeomFillTrihedron, GuideTrihedronAC, GuideTrihedronPlan, LocationGuide,
    LocationLaw as GeomFillLocationLaw,
};
use occt_geom::law::LawInterpol;

use crate::brep_adaptor_comp_curve::CompCurveAdaptor;
use crate::brep_fill::search_origin;
use crate::brep_fill_acr_law::AcrLaw;
use crate::brep_fill_compatible_wires::CompatibleWires;
use crate::brep_fill_edge3d_law::Edge3DLaw;
use crate::brep_fill_edge_on_surf_law::EdgeOnSurfLaw;
use crate::brep_fill_location_law::LocationLaw;
use crate::brep_fill_section::Section;
use crate::brep_fill_section_law::SectionLaw;
use crate::brep_fill_sweep::TransitionStyle;
use crate::brep_gprop_full::linear_properties;
use crate::brep_lib_make_face::make_face_from_wire;
use crate::brep_tool::BRepTool;
use crate::shape::{Shell, TopoShape, Vertex, Wire};
use crate::topo_tools::{is_same, shape_map_key, wire_vertices, ShapeMapKey};

pub use self::build::{build_boundaries, is_same_oriented, perform_plan};
pub use self::prepare::{compute_section, perform_transition};
pub use self::simulate::simulate;

/// `BRepFill_TypeOfContact` (`BRepFill_TypeOfContact.hxx`): how the auxiliary
/// spine constrains the swept section.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TypeOfContact {
    /// `BRepFill_NoContact`: the normal is `PQ`, nothing touches the guide.
    NoContact,
    /// `BRepFill_Contact`: the section stays in contact with the guide.
    Contact,
    /// `BRepFill_ContactOnBorder`: the guide becomes a boundary of the sweep.
    ContactOnBorder,
}

/// The smallest V tolerance the pipe shell hands to `BRepFill_Sweep`
/// (`cxx:760`: the literal `1.e-5` for the 2D tolerance).
const DEFAULT_TOL_2D: f64 = 1.0e-5;

/// `BRepFill_PipeShell` (`BRepFill_PipeShell.hxx:47-261`): the class data.
pub struct PipeShell {
    /// `mySpine` (`hxx:231`).
    my_spine: Wire,
    /// `myFirst` (`hxx:232`), a null shape is `None`.
    my_first: Option<TopoShape>,
    /// `myLast` (`hxx:233`).
    my_last: Option<TopoShape>,
    /// `myShape` (`hxx:234`).
    my_shape: Option<TopoShape>,
    /// `mySeq` (`hxx:235`), 1-based in OCCT.
    my_seq: Vec<Section>,
    /// `WSeq` (`hxx:236`).
    w_seq: Vec<Wire>,
    /// `myIndOfSec` (`hxx:237`).
    my_ind_of_sec: Vec<usize>,
    /// `myEdgeNewEdges` (`hxx:238-239`).
    my_edge_new_edges: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myGenMap` (`hxx:240-241`).
    my_gen_map: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myTol3d` (`hxx:242`).
    my_tol_3d: f64,
    /// `myBoundTol` (`hxx:243`).
    my_bound_tol: f64,
    /// `myTolAngular` (`hxx:244`).
    my_tol_angular: f64,
    /// `angmin` (`hxx:245`).
    ang_min: f64,
    /// `angmax` (`hxx:246`).
    ang_max: f64,
    /// `myMaxDegree` (`hxx:247`).
    my_max_degree: i32,
    /// `myMaxSegments` (`hxx:248`).
    my_max_segments: i32,
    /// `myForceApproxC1` (`hxx:249`).
    my_force_approx_c1: bool,
    /// `myLaw` (`hxx:250`).
    my_law: Option<Arc<dyn occt_geom::law::LawFunction>>,
    /// `myIsAutomaticLaw` (`hxx:251`).
    my_is_automatic_law: bool,
    /// `myLocation` (`hxx:252`).
    my_location: Option<Arc<dyn LocationLaw>>,
    /// `mySection` (`hxx:253`).
    my_section: Option<Arc<dyn SectionLaw>>,
    /// `myTrihedron` (`hxx:256`).
    my_trihedron: GeomFillTrihedron,
    /// `myTransition` (`hxx:257`).
    my_transition: TransitionStyle,
    /// `myStatus` (`hxx:258`).
    my_status: GeomFillPipeError,
    /// `myErrorOnSurf` (`hxx:259`).
    my_error_on_surf: f64,
    /// `myIsBuildHistory` (`hxx:260`).
    my_is_build_history: bool,
}

impl PipeShell {
    /// `BRepFill_PipeShell(Spine)` (`cxx:223-252`).
    ///
    /// `angmin` / `angmax` are not initialised by the OCCT constructor (only
    /// `SetTransition` writes them); the port starts them at the
    /// `SetTransition` defaults so the value handed to `SetAngularControl` is
    /// defined.
    pub fn new(spine: &Wire) -> Self {
        // `cxx:244-251`: "Attention to closed non-declared wire !".
        if !spine.is_closed() {
            let (vf, vl) = crate::topo_tools::wire_vertices(spine);
            if let (Some(vf), Some(vl)) = (vf, vl) {
                if is_same(&vf.0, &vl.0) {
                    spine.set_closed(true);
                }
            }
        }
        let mut this = Self {
            my_spine: spine.clone(),
            my_first: None,
            my_last: None,
            my_shape: None,
            my_seq: Vec::new(),
            w_seq: Vec::new(),
            my_ind_of_sec: Vec::new(),
            my_edge_new_edges: HashMap::new(),
            my_gen_map: HashMap::new(),
            my_tol_3d: 0.0,
            my_bound_tol: 0.0,
            my_tol_angular: 0.0,
            ang_min: 0.01,
            ang_max: 6.0,
            my_max_degree: 11,
            my_max_segments: 100,
            my_force_approx_c1: false,
            my_law: None,
            my_is_automatic_law: false,
            my_location: None,
            my_section: None,
            my_trihedron: GeomFillTrihedron::CorrectedFrenet,
            my_transition: TransitionStyle::Modified,
            my_status: GeomFillPipeError::PipeOk,
            my_error_on_surf: 0.0,
            my_is_build_history: true,
        };
        this.set_tolerance(1.0e-4, 1.0e-4, 1.0e-2);
        this
    }

    /// `Set(IsFrenet)` (`cxx:256-273`).
    pub fn set_frenet(&mut self, is_frenet: bool) -> Result<(), String> {
        let trihedron: Arc<dyn occt_geom::geom_fill::TrihedronLaw> = if is_frenet {
            Arc::new(Frenet::new())
        } else {
            Arc::new(CorrectedFrenet::new())
        };
        self.my_trihedron = GeomFillTrihedron::Frenet;
        let loc: Arc<dyn GeomFillLocationLaw> =
            Arc::new(CurveAndTrihedron::new(trihedron));
        self.my_location = Some(Arc::new(Edge3DLaw::new(&self.my_spine, &loc)));
        // "It is required to relocalize sections."
        self.my_section = None;
        Ok(())
    }

    /// `SetDiscrete()` (`cxx:279-289`).
    pub fn set_discrete(&mut self) -> Result<(), String> {
        self.my_trihedron = GeomFillTrihedron::DiscreteTrihedron;
        let trihedron: Arc<dyn occt_geom::geom_fill::TrihedronLaw> =
            Arc::new(DiscreteTrihedron::new());
        let loc: Arc<dyn GeomFillLocationLaw> =
            Arc::new(CurveAndTrihedron::new(trihedron));
        self.my_location = Some(Arc::new(Edge3DLaw::new(&self.my_spine, &loc)));
        self.my_section = None;
        Ok(())
    }

    /// `Set(const gp_Ax2&)` (`cxx:292-304`).
    pub fn set_ax2(&mut self, axe: &GpAx2) -> Result<(), String> {
        self.my_trihedron = GeomFillTrihedron::Fixed;
        let v1 = GpVec::from_xyz(&axe.direction().xyz());
        let v2 = GpVec::from_xyz(&axe.x_direction().xyz());
        let trihedron: Arc<dyn occt_geom::geom_fill::TrihedronLaw> =
            Arc::new(Fixed::new(&v1, &v2)?);
        let loc: Arc<dyn GeomFillLocationLaw> =
            Arc::new(CurveAndTrihedron::new(trihedron));
        self.my_location = Some(Arc::new(Edge3DLaw::new(&self.my_spine, &loc)));
        self.my_section = None;
        Ok(())
    }

    /// `Set(const gp_Dir&)` (`cxx:309-317`).
    pub fn set_bi_normal(&mut self, bi_normal: &GpDir) -> Result<(), String> {
        self.my_trihedron = GeomFillTrihedron::ConstantNormal;
        let trihedron: Arc<dyn occt_geom::geom_fill::TrihedronLaw> =
            Arc::new(ConstantBiNormal::new(bi_normal));
        let loc: Arc<dyn GeomFillLocationLaw> =
            Arc::new(CurveAndTrihedron::new(trihedron));
        self.my_location = Some(Arc::new(Edge3DLaw::new(&self.my_spine, &loc)));
        self.my_section = None;
        Ok(())
    }

    /// `Set(const TopoDS_Shape& SpineSupport)` (`cxx:323-337`): build the
    /// Darboux location law from the pcurves of the spine on `spine_support`.
    ///
    /// OCCT answers `false` (and leaves `myLocation` / `myTrihedron` /
    /// `mySection` untouched) when the law has no result.
    pub fn set_spine_support(&mut self, spine_support: &TopoShape) -> Result<bool, String> {
        // `BRepFill_EdgeOnSurfLaw loc(mySpine, SpineSupport)` (`cxx:328`).
        let loc = EdgeOnSurfLaw::new(&self.my_spine, spine_support);
        // `B = loc->HasResult()` (`cxx:329`).
        let b = loc.has_result();
        if b {
            // `myLocation = loc; myTrihedron = GeomFill_IsDarboux;
            //  mySection.Nullify()` (`cxx:331-335`).
            self.my_location = Some(Arc::new(loc));
            self.my_trihedron = GeomFillTrihedron::Darboux;
            self.my_section = None; // It is required to relocalize the sections.
        }
        Ok(b)
    }

    /// `Set(AuxiliarySpine, CurvilinearEquivalence, KeepContact)`
    /// (`cxx:343-431`): put a guide law in place of the spine trihedron.
    ///
    /// `cxx:356-396` reorients the guide first — an open guide against an open
    /// spine goes through `BRepFill_CompatibleWires`, a *closed* guide goes
    /// through `BRepFill::SearchOrigin` on the spine's start point and start
    /// tangent. `cxx:398-430` then builds either a curvilinear-equivalence law
    /// (`GeomFill_GuideTrihedronAC` + `BRepFill_ACRLaw`) or a plane law
    /// (`GeomFill_GuideTrihedronPlan` + `BRepFill_Edge3DLaw`).
    pub fn set_auxiliary_spine(
        &mut self,
        auxiliary_spine: &Wire,
        curvilinear_equivalence: bool,
        keep_contact: TypeOfContact,
    ) -> Result<(), String> {
        // `cxx:346-350`: `TheGuide = AuxiliarySpine; SpClose = mySpine.Closed();
        // GuideClose = AuxiliarySpine.Closed()`.
        let mut the_guide = auxiliary_spine.clone();
        let sp_close = self.my_spine.is_closed();
        let guide_close = auxiliary_spine.is_closed();

        // `cxx:352-355`.
        if keep_contact == TypeOfContact::ContactOnBorder {
            self.my_is_automatic_law = true;
        }

        if !sp_close && !guide_close {
            // `cxx:356-372`: open reorientation of the guide.
            let seq = [self.my_spine.0.clone(), the_guide.0.clone()];
            let mut cw = CompatibleWires::with_sections(&seq);
            cw.set_percent(0.1);
            cw.perform(false)?;
            if !cw.is_done() {
                return Err("StdFail_NotDone: Uncompatible wires".to_string());
            }
            // `TopoDS::Wire(CW.Shape().Value(2))`.
            let second = cw
                .shape()
                .get(1)
                .cloned()
                .ok_or_else(|| "Standard_OutOfRange: CW.Shape().Value(2)".to_string())?;
            the_guide = Wire::wrap(second).ok_or_else(|| {
                "Standard_TypeMismatch: BRepFill_PipeShell::Set".to_string()
            })?;
        } else if guide_close {
            // `cxx:373-396`: closed guide, determine the origin and reorient.
            let (sp_or, dir) = if !sp_close {
                let (vf, vl) = wire_vertices(&self.my_spine);
                let (Some(vf), Some(vl)) = (vf, vl) else {
                    return Err(
                        "BRepFill_PipeShell::Set: the spine has no boundary vertex".to_string()
                    );
                };
                let sp_or = BRepTool::vertex_point(&vf);
                let p = BRepTool::vertex_point(&vl);
                // `gp_Vec V(P, SpOr)`.
                let v = GpVec::from_pnts(&p, &sp_or);
                // `SpOr.BaryCenter(0.5, P, 0.5)`.
                let mid = GpPnt::from_xyz(
                    &sp_or.xyz().multiplied(0.5).added(&p.xyz().multiplied(0.5)),
                );
                (mid, v)
            } else {
                // `BRepAdaptor_CompCurve BC(mySpine); BC.D1(0, SpOr, Dir)`.
                let bc = CompCurveAdaptor::of_wire(&self.my_spine, false);
                bc.d1(0.0)
            };
            search_origin(&mut the_guide, &sp_or, &dir, 100.0 * self.my_tol_3d)?;
        }

        // `cxx:398-399`: `Guide = new BRepAdaptor_CompCurve(TheGuide)`; the
        // default `KnotByCurvilinearAbcissa` is `false`.
        let guide: Arc<dyn Curve> =
            Arc::new(CompCurveAdaptor::of_wire(&the_guide, false));

        if curvilinear_equivalence {
            // `cxx:401-414`.
            if keep_contact == TypeOfContact::Contact
                || keep_contact == TypeOfContact::ContactOnBorder
            {
                self.my_trihedron = GeomFillTrihedron::GuideACWithContact;
            } else {
                self.my_trihedron = GeomFillTrihedron::GuideAC;
            }
            let t_law = Arc::new(GuideTrihedronAC::new(Arc::clone(&guide)));
            let loc = LocationGuide::new(t_law)?;
            let law = AcrLaw::new(&self.my_spine, &loc)?;
            self.my_location = Some(Arc::new(law));
        } else {
            // `cxx:415-429`.
            if keep_contact == TypeOfContact::Contact
                || keep_contact == TypeOfContact::ContactOnBorder
            {
                self.my_trihedron = GeomFillTrihedron::GuidePlanWithContact;
            } else {
                self.my_trihedron = GeomFillTrihedron::GuidePlan;
            }
            let t_law = Arc::new(GuideTrihedronPlan::new(Arc::clone(&guide)));
            let loc: Arc<dyn GeomFillLocationLaw> = Arc::new(LocationGuide::new(t_law)?);
            let law = Edge3DLaw::new(&self.my_spine, &loc);
            self.my_location = Some(Arc::new(law));
        }

        // `cxx:430`: `mySection.Nullify()` — the sections must be relocalized.
        self.my_section = None;
        Ok(())
    }

    /// `SetMaxDegree(NewMaxDegree)` (`cxx:435-438`).
    pub fn set_max_degree(&mut self, new_max_degree: i32) {
        self.my_max_degree = new_max_degree;
    }

    /// `SetMaxSegments(NewMaxSegments)` (`cxx:442-445`).
    pub fn set_max_segments(&mut self, new_max_segments: i32) {
        self.my_max_segments = new_max_segments;
    }

    /// `SetForceApproxC1(ForceApproxC1)` (`cxx:453-456`).
    pub fn set_force_approx_c1(&mut self, force_approx_c1: bool) {
        self.my_force_approx_c1 = force_approx_c1;
    }

    /// `Add(Profile, WithContact, WithCorrection)` (`cxx:462-471`).
    pub fn add(
        &mut self,
        profile: &TopoShape,
        with_contact: bool,
        with_correction: bool,
    ) -> Result<(), String> {
        self.add_at(profile, &Vertex::new(), with_contact, with_correction)?;
        self.reset_loc();
        Ok(())
    }

    /// `Add(Profile, Location, WithContact, WithCorrection)` (`cxx:475-548`).
    pub fn add_at(
        &mut self,
        profile: &TopoShape,
        location: &Vertex,
        with_contact: bool,
        with_correction: bool,
    ) -> Result<(), String> {
        let _ = self.delete_profile(profile); // No duplication
        if self.my_is_automatic_law {
            // `cxx:479-540`: the automatic homothety law of a guide sweep.
            self.add_automatic_law(profile, location, with_contact, with_correction)?;
        } else {
            let a_section =
                Section::with_profile(profile, location, with_contact, with_correction)?;
            self.my_seq.push(a_section);
            self.my_section = None;
            self.reset_loc();
        }
        Ok(())
    }

    /// The `myIsAutomaticLaw` arm of `Add` (`cxx:479-540`): rebuild the single
    /// section, ask the guide law for the width/abscissa samples, normalise
    /// them by the initial section width, and interpolate them into `myLaw`.
    fn add_automatic_law(
        &mut self,
        profile: &TopoShape,
        location: &Vertex,
        with_contact: bool,
        with_correction: bool,
    ) -> Result<(), String> {
        // `cxx:481-487`.
        self.my_seq.clear();
        let mut a_section = Section::with_profile(profile, location, with_contact, with_correction)?;
        a_section.set(true);
        self.my_seq.push(a_section);
        self.my_section = None;
        self.reset_loc();

        // `cxx:489-491`: `Loc = down_cast<GeomFill_LocationGuide>(myLocation->Law(1))`.
        let law_handle = self
            .my_location
            .as_ref()
            .and_then(|loc| loc.base().law(1))
            .ok_or_else(|| {
                "Standard_NoSuchObject: BRepFill_PipeShell::Add: myLocation->Law(1)".to_string()
            })?;
        let guide_law = law_handle.as_location_guide().ok_or_else(|| {
            "Standard_TypeMismatch: BRepFill_PipeShell::Add automatic law".to_string()
        })?;

        // `cxx:492-494`: `Loc->ComputeAutomaticLaw(ParAndRad)`.
        let (status, mut par_and_rad) = guide_law.compute_automatic_law();
        if status != GeomFillPipeError::PipeOk || par_and_rad.is_empty() {
            return Err(
                "GeomFill_LocationGuide::ComputeAutomaticLaw: no law computed".to_string(),
            );
        }

        // `cxx:496-499`: the barycentre of the profile, `BRepGProp::LinearProperties`.
        let global_props = linear_properties(profile)?;
        let barycenter = global_props.center();

        // `cxx:501-503`: the plane of the profile face, `BRepLib_MakeFace(W, true)`.
        let profile_wire = Wire::wrap(profile.clone())
            .ok_or_else(|| "Standard_TypeMismatch: Add(Profile): not a wire".to_string())?;
        let profile_face = make_face_from_wire(&profile_wire, true)
            .ok_or_else(|| "StdFail_NotDone: BRepLib_MakeFace".to_string())?;
        let the_plane = BRepTool::face_surface(&profile_face)
            .ok_or_else(|| "Standard_NullObject: BRep_Tool::Surface".to_string())?;
        let plane = the_plane
            .gp_pln()
            .ok_or_else(|| "Standard_TypeMismatch: BRepLib_MakeFace: not a plane".to_string())?;
        let plane_pnt = plane.location();
        // `gp_Pln::Axis().Direction()` is a unit `gp_Dir`; the intersector takes
        // the plane normal as a `GpVec` (`IntSurf_Quadric::Gradient`).
        let plane_normal = GpVec::from_xyz(plane.axis().direction().xyz());

        // `cxx:504-529`: for the location curve and the guide, the intersection
        // point nearest to the barycentre of the profile.
        let curves: [Option<Arc<dyn Curve>>; 2] =
            [guide_law.get_curve(), guide_law.guide()];
        let mut points_on_spines: [Option<GpPnt>; 2] = [None, None];
        for (i, curve) in curves.iter().enumerate() {
            let curve = curve.as_ref().ok_or_else(|| {
                "Standard_NullObject: GeomFill_LocationGuide::GetCurve".to_string()
            })?;
            let mut min_dist = f64::MAX;
            for (_, pnt) in plane_intersections(curve.as_ref(), &plane_pnt, &plane_normal) {
                let dist = barycenter.distance(&pnt);
                if dist < min_dist {
                    min_dist = dist;
                    points_on_spines[i] = Some(pnt);
                }
            }
        }
        let (Some(p0), Some(p1)) = (points_on_spines[0], points_on_spines[1]) else {
            return Err(
                "BRepFill_PipeShell::Add: the plane misses a guide curve".to_string()
            );
        };

        // `cxx:531-538`: `InitialWidth`, then `Y /= InitialWidth`.
        let initial_width = p0.distance(&p1);
        let nb_par_rad = par_and_rad.len();
        for slot in par_and_rad.iter_mut() {
            slot.set_y(slot.y() / initial_width);
        }

        // `cxx:540-547`: `Law_Interpol::Set(ParAndRad->Array1(), IsPeriodic)`.
        let is_periodic =
            (par_and_rad[0].y() - par_and_rad[nb_par_rad - 1].y()).abs()
                < occt_core::precision::CONFUSION;
        let pairs: Vec<(f64, f64)> =
            par_and_rad.iter().map(|p| (p.x(), p.y())).collect();
        let mut law = LawInterpol::new();
        law.set(&pairs, is_periodic)?;
        self.my_law = Some(Arc::new(law));
        Ok(())
    }

    /// `SetLaw(Profile, L, WithContact, WithCorrection)` (`cxx:554-564`).
    pub fn set_law(
        &mut self,
        profile: &TopoShape,
        l: Arc<dyn occt_geom::law::LawFunction>,
        with_contact: bool,
        with_correction: bool,
    ) -> Result<(), String> {
        self.set_law_at(profile, l, &Vertex::new(), with_contact, with_correction)?;
        self.reset_loc();
        Ok(())
    }

    /// `SetLaw(Profile, L, Location, WithContact, WithCorrection)`
    /// (`cxx:569-582`).
    pub fn set_law_at(
        &mut self,
        profile: &TopoShape,
        l: Arc<dyn occt_geom::law::LawFunction>,
        location: &Vertex,
        with_contact: bool,
        with_correction: bool,
    ) -> Result<(), String> {
        self.my_seq.clear();
        let mut a_section =
            Section::with_profile(profile, location, with_contact, with_correction)?;
        a_section.set(true);
        self.my_seq.push(a_section);
        self.my_law = Some(l);
        self.my_section = None;
        self.reset_loc();
        Ok(())
    }

    /// `DeleteProfile(Profile)` (`cxx:586-606`).
    pub fn delete_profile(&mut self, profile: &TopoShape) -> Result<(), String> {
        let mut trouve = false;
        let mut ii = 0usize;
        while ii < self.my_seq.len() && !trouve {
            let same = self.my_seq[ii]
                .original_shape()
                .map(|shape| is_same(shape, profile))
                .unwrap_or(false);
            if same {
                trouve = true;
                self.my_seq.remove(ii);
            } else {
                ii += 1;
            }
        }
        if trouve {
            self.my_section = None;
        }
        self.reset_loc();
        Ok(())
    }

    /// `IsReady()` (`cxx:609-612`).
    pub fn is_ready(&self) -> bool {
        !self.my_seq.is_empty()
    }

    /// `GetStatus()` (`cxx:616-619`).
    pub fn get_status(&self) -> GeomFillPipeError {
        self.my_status
    }

    /// `SetTolerance(Tol3d, BoundTol, TolAngular)` (`cxx:623-630`).
    pub fn set_tolerance(&mut self, tol_3d: f64, bound_tol: f64, tol_angular: f64) {
        self.my_tol_3d = tol_3d;
        self.my_bound_tol = bound_tol;
        self.my_tol_angular = tol_angular;
    }

    /// `SetTransition(Mode, Angmin, Angmax)` (`cxx:637-649`).
    pub fn set_transition(&mut self, mode: TransitionStyle, ang_min: f64, ang_max: f64) {
        if self.my_transition != mode {
            self.my_section = None; // It is required to relocalize the sections.
        }
        self.my_transition = mode;
        self.ang_min = ang_min;
        self.ang_max = ang_max;
    }

    /// `SetIsBuildHistory` (`hxx:125-128`).
    pub fn set_is_build_history(&mut self, is_build_history: bool) {
        self.my_is_build_history = is_build_history;
    }

    /// `IsBuildHistory` (`hxx:131`).
    pub fn is_build_history(&self) -> bool {
        self.my_is_build_history
    }

    /// `Shape()` (`cxx:917-920`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `ErrorOnSurface()` (`cxx:924-927`).
    pub fn error_on_surface(&self) -> f64 {
        self.my_error_on_surf
    }

    /// `FirstShape()` (`cxx:932-935`).
    pub fn first_shape(&self) -> Option<&TopoShape> {
        self.my_first.as_ref()
    }

    /// `LastShape()` (`cxx:938-941`).
    pub fn last_shape(&self) -> Option<&TopoShape> {
        self.my_last.as_ref()
    }

    /// `Profiles(Profiles)` (`hxx:203-207`).
    pub fn profiles(&self) -> Vec<TopoShape> {
        self.my_seq
            .iter()
            .filter_map(|sec| sec.original_shape().cloned())
            .collect()
    }

    /// `Spine()` (`hxx:211-214`).
    pub fn spine(&self) -> &Wire {
        &self.my_spine
    }

    /// `Generated(S, L)` (`cxx:946-956`).
    pub fn generated(&self, the_shape: &TopoShape) -> Vec<TopoShape> {
        self.my_gen_map
            .get(&shape_map_key(the_shape))
            .cloned()
            .unwrap_or_default()
    }

    // ---- accessors used by the sub-modules ----

    pub(super) fn spine_ref(&self) -> &Wire {
        &self.my_spine
    }

    pub(super) fn seq(&self) -> &[Section] {
        &self.my_seq
    }

    pub(super) fn location(&self) -> Option<&Arc<dyn LocationLaw>> {
        self.my_location.as_ref()
    }

    pub(super) fn location_mut(&mut self) -> &mut Option<Arc<dyn LocationLaw>> {
        &mut self.my_location
    }

    pub(super) fn section(&self) -> Option<&Arc<dyn SectionLaw>> {
        self.my_section.as_ref()
    }

    pub(super) fn section_mut(&mut self) -> &mut Option<Arc<dyn SectionLaw>> {
        &mut self.my_section
    }

    pub(super) fn tol_3d(&self) -> f64 {
        self.my_tol_3d
    }

    pub(super) fn bound_tol(&self) -> f64 {
        self.my_bound_tol
    }

    pub(super) fn tol_angular(&self) -> f64 {
        self.my_tol_angular
    }

    pub(super) fn tol_2d(&self) -> f64 {
        DEFAULT_TOL_2D
    }

    pub(super) fn ang_min(&self) -> f64 {
        self.ang_min
    }

    pub(super) fn ang_max(&self) -> f64 {
        self.ang_max
    }

    pub(super) fn max_degree(&self) -> i32 {
        self.my_max_degree
    }

    pub(super) fn max_segments(&self) -> i32 {
        self.my_max_segments
    }

    pub(super) fn force_approx_c1(&self) -> bool {
        self.my_force_approx_c1
    }

    pub(super) fn trihedron(&self) -> GeomFillTrihedron {
        self.my_trihedron
    }

    pub(super) fn transition(&self) -> TransitionStyle {
        self.my_transition
    }

    pub(super) fn set_status(&mut self, status: GeomFillPipeError) {
        self.my_status = status;
    }

    pub(super) fn set_shape(&mut self, shape: Option<TopoShape>) {
        self.my_shape = shape;
    }

    pub(super) fn set_first(&mut self, shape: Option<TopoShape>) {
        self.my_first = shape;
    }

    pub(super) fn set_last(&mut self, shape: Option<TopoShape>) {
        self.my_last = shape;
    }

    pub(super) fn first_ref(&self) -> Option<&TopoShape> {
        self.my_first.as_ref()
    }

    pub(super) fn first_mut(&mut self) -> &mut Option<TopoShape> {
        &mut self.my_first
    }

    pub(super) fn last_mut(&mut self) -> &mut Option<TopoShape> {
        &mut self.my_last
    }

    pub(super) fn shape_mut(&mut self) -> Option<&mut TopoShape> {
        self.my_shape.as_mut()
    }

    pub(super) fn last_ref(&self) -> Option<&TopoShape> {
        self.my_last.as_ref()
    }

    pub(super) fn set_error_on_surf(&mut self, error: f64) {
        self.my_error_on_surf = error;
    }

    pub(super) fn w_seq_mut(&mut self) -> &mut Vec<Wire> {
        &mut self.w_seq
    }

    pub(super) fn ind_of_sec_mut(&mut self) -> &mut Vec<usize> {
        &mut self.my_ind_of_sec
    }

    pub(super) fn edge_new_edges_mut(&mut self) -> &mut HashMap<ShapeMapKey, Vec<TopoShape>> {
        &mut self.my_edge_new_edges
    }

    pub(super) fn law_ref(&self) -> Option<&Arc<dyn occt_geom::law::LawFunction>> {
        self.my_law.as_ref()
    }

    /// The empty shell OCCT stores in `myShape` when `Build` fails
    /// (`cxx:714-720`, `cxx:833-838`).
    pub(super) fn empty_shell() -> Option<TopoShape> {
        Some(Shell::new().0)
    }
}
