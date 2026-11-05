//! `FilletSurf_InternalBuilder` — the private worker of `FilletSurf_Builder`.
//!
//! Source: `ModelingAlgorithms/TKFillet/FilletSurf/FilletSurf_InternalBuilder.cxx`
//! (OCCT 8.0.0), `hxx:40-235`.
//!
//! OCCT derives the class from `ChFi3d_FilBuilder` and overrides the four
//! `PerformSurf`; the overridden half lives in [`ops`], the contour surgery and
//! the result accessors here. `Done()` answers the inherited `ChFi3d_Builder::
//! done` member, which `Perform` never writes (it bypasses `Compute`
//! deliberately, `cxx:320-338`) — see [`FilletSurfInternalBuilder::done`].
//!
//! Divergences that are not heuristics:
//! - `Add` answers `Result<i32, String>`: the six numbered codes are the OCCT
//!   return values (`cxx:171-179`), the `Err` is the `Standard_Failure`
//!   `ChFi3d_FilBuilder::Add` and `Spine::D1` may raise — OCCT lets it escape
//!   the caller's constructor (`FilletSurf_Builder.cxx:41`),
//! - `TopoDS::Edge(S)` of a shape that is not an edge raises
//!   `Standard_TypeMismatch` in OCCT; the port answers the same `4` the null
//!   edge gets (`cxx:190-193`),
//! - an empty `myListStripe` makes `First()` raise in OCCT; the accessors
//!   answer `None` / `0` instead.

pub mod edge_parameter;
pub mod ops;

use std::sync::Arc;

use occt_core::gp::GpCirc;
use occt_geom::{Curve as GeomCurve, GeomCircle, GeomTrimmedCurve, Surface as GeomSurface};
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::brep_lib_continuity::stored_continuity;
use crate::chfi3d::builder::ChFi3dFilBuilder;
use crate::chfi3d::builder_0_helpers::ch_fi3d_perform_el_spine;
use crate::chfi3d::fillet_shape::FilletShape;
use crate::chfi_ds::{ElSpine, Spine, StripeHandle, SurfData};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::is_same;

use self::edge_parameter::compute_edge_parameter;
use self::ops::FilletSurfOps;
use super::status::StatusType;

/// `FilletSurf_InternalBuilder` (`FilletSurf_InternalBuilder.hxx:40-235`).
pub struct FilletSurfInternalBuilder {
    /// The `ChFi3d_FilBuilder` base class sub-object.
    base: ChFi3dFilBuilder,
}

impl FilletSurfInternalBuilder {
    /// `FilletSurf_InternalBuilder(S, FShape, Ta, Tapp3d, Tapp2d)`
    /// (`cxx:155-166`): `ChFi3d_FilBuilder(S, FShape, Ta)`, then
    /// `SetParams(Ta, Tapp3d, Tapp2d, Tapp3d, Tapp2d, 1.e-3)` and
    /// `SetContinuity(GeomAbs_C2, Ta)`.
    pub fn new(shape: &TopoShape, fshape: FilletShape, ta: f64, tapp3d: f64, tapp2d: f64) -> Self {
        let mut this = Self {
            base: ChFi3dFilBuilder::new(shape, fshape, ta),
        };
        this.base
            .base_mut()
            .set_params(ta, tapp3d, tapp2d, tapp3d, tapp2d, 1.0e-3);
        this.base
            .base_mut()
            .set_continuity(occt_core::kernel::geomabs::Shape::C2, ta);
        this
    }

    /// The `ChFi3d_FilBuilder` sub-object.
    pub fn base(&self) -> &ChFi3dFilBuilder {
        &self.base
    }

    /// The `ChFi3d_Builder` half of the sub-object.
    pub fn builder(&self) -> &crate::chfi3d::builder::ChFi3dBuilder {
        self.base.base()
    }

    /// `Add(E, R)` (`cxx:180-318`): creates the contour on the set of edges.
    ///
    /// The six answer codes are `cxx:171-179`:
    /// `0` no problem, `1` empty list, `2` the edges are not G1, `3` two
    /// connected faces on a same support are not G1, `4` the edge is not on the
    /// shape, `5` the edge is not sharp.
    pub fn add(&mut self, edges: &[TopoShape], radius: f64) -> Result<i32, String> {
        // `if (E.IsEmpty()) { return 1; }` (`cxx:182-185`).
        if edges.is_empty() {
            return Ok(1);
        }

        // The sharpness pre-check (`cxx:187-227`).
        for value in edges {
            let Some(cured) = Edge::wrap(value.clone()) else {
                return Ok(4);
            };
            // `if (!myEFMap.Contains(cured)) { return 4; }` (`cxx:194-197`).
            if !self.base.base().ef_map().contains(&cured.0) {
                return Ok(4);
            }
            // `for (It.Initialize(myEFMap(cured)); ...)` (`cxx:199-215`).
            let mut ff1: Option<TopoShape> = None;
            let mut ff2: Option<TopoShape> = None;
            for value in self.base.base().ef_map().get_from_key(&cured.0) {
                match ff1.as_ref() {
                    None => ff1 = Some(value.clone()),
                    Some(first) => {
                        ff2 = Some(value.clone());
                        if !is_same(&value.clone(), first) {
                            break;
                        }
                    }
                }
            }
            let (Some(ff1), Some(ff2)) = (ff1, ff2) else {
                return Ok(5);
            };
            // `if (ff1.IsSame(ff2)) { return 5; }` (`cxx:217-221`).
            if is_same(&ff1, &ff2) {
                return Ok(5);
            }
            // `if (BRep_Tool::Continuity(cured, ff1, ff2) != GeomAbs_C0)`
            // (`cxx:222-226`). `stored_continuity` answers `C0` (0) for an
            // edge that was never encoded, exactly like `BRep_Tool`.
            let cont = stored_continuity(&cured, &Face(ff1), &Face(ff2));
            if cont != 0 {
                return Ok(5);
            }
        }

        // `TopoDS_Edge ed = TopoDS::Edge(E.First()); ed.Orientation(
        // TopAbs_FORWARD); ChFi3d_FilBuilder::Add(R, ed);` (`cxx:229-230`).
        let Some(mut ed) = Edge::wrap(edges[0].clone()) else {
            return Ok(4);
        };
        ed.0.set_orientation(Orientation::Forward);
        self.base.add_with_radius(radius, &ed)?;

        // `occ::handle<ChFiDS_Stripe> st = myListStripe.First();
        // occ::handle<ChFiDS_Spine>& sp = st->ChangeSpine();` (`cxx:233-235`).
        let stripe = self.first_stripe().ok_or_else(|| {
            "FilletSurf_InternalBuilder::Add: the contour was not created".to_string()
        })?;
        let mut spine = {
            let guard = stripe
                .read()
                .map_err(|_| "FilletSurf_InternalBuilder::Add: poisoned stripe".to_string())?;
            guard.spine().cloned()
        }
        .ok_or_else(|| "FilletSurf_InternalBuilder::Add: the stripe has no spine".to_string())?;
        let periodic = spine
            .read()
            .map_err(|_| "FilletSurf_InternalBuilder::Add: poisoned spine".to_string())?
            .is_periodic();

        // `for (It.Initialize(E); ...) { if (!Contains(cured)) { return 2; } }`
        // (`cxx:238-246`).
        for value in edges {
            let Some(cured) = Edge::wrap(value.clone()) else {
                return Ok(2);
            };
            if self.base.base().contains(&cured) == 0 {
                return Ok(2);
            }
        }

        // `occ::handle<ChFiDS_FilSpine> newsp = new ChFiDS_FilSpine(); bool
        // debut = false; int premierquinyestpas = 0; int yatrou = 0;`
        // (`cxx:247-251`).
        let mut newsp = Spine::new_fil_spine();
        let mut debut = false;
        let mut premier_qui_n_y_est_pas = 0usize;
        let mut yatrou = 0usize;
        // `for (int i = 1; i <= sp->NbEdges(); i++)` (`cxx:252-269`).
        let nb_edges = spine
            .read()
            .map_err(|_| "FilletSurf_InternalBuilder::Add: poisoned spine".to_string())?
            .nb_edges();
        for i in 1..=nb_edges {
            let cured = {
                let guard = spine
                    .read()
                    .map_err(|_| "FilletSurf_InternalBuilder::Add: poisoned spine".to_string())?;
                guard.edges(i).cloned()
            };
            let Some(cured) = cured else {
                continue;
            };
            if is_in_list(&cured, edges) {
                debut = true;
                if premier_qui_n_y_est_pas != 0 {
                    yatrou = 1;
                    break;
                }
                newsp.set_edges(&cured);
            } else if debut && premier_qui_n_y_est_pas == 0 {
                premier_qui_n_y_est_pas = i;
            }
        }

        // `if (!periodic && yatrou) { return 2; }` (`cxx:270-273`).
        if !periodic && yatrou != 0 {
            return Ok(2);
        }
        // The periodic wraparound (`cxx:274-293`).
        if periodic && yatrou != 0 {
            let mut vrai_trou = false;
            let mut a_local_debut = false;
            let mut i = nb_edges;
            while i > yatrou {
                let cured = {
                    let guard = spine.read().map_err(|_| {
                        "FilletSurf_InternalBuilder::Add: poisoned spine".to_string()
                    })?;
                    guard.edges(i).cloned()
                };
                if let Some(cured) = cured {
                    if is_in_list(&cured, edges) {
                        if vrai_trou {
                            return Ok(2);
                        }
                        newsp.put_in_first(&cured);
                    } else if a_local_debut {
                        vrai_trou = true;
                    }
                }
                a_local_debut = true;
                i -= 1;
            }
        }

        // `if (newsp->NbEdges() != sp->NbEdges()) { newsp->Load();
        // newsp->SetRadius(R); sp = newsp; }` (`cxx:295-301`).
        if newsp.nb_edges() != nb_edges {
            newsp.load();
            newsp
                .set_radius(radius)
                .ok_or_else(|| "FilletSurf_InternalBuilder::Add: SetRadius failed".to_string())?;
            spine = Arc::new(std::sync::RwLock::new(newsp));
            stripe
                .write()
                .map_err(|_| "FilletSurf_InternalBuilder::Add: poisoned stripe".to_string())?
                .set_spine(Arc::clone(&spine));
        }

        // `occ::handle<ChFiDS_ElSpine> hels = new ChFiDS_ElSpine(); gp_Vec
        // TFirst, TLast; gp_Pnt PFirst, PLast; sp->D1(sp->FirstParameter(),
        // PFirst, TFirst); sp->D1(sp->LastParameter(), PLast, TLast);`
        // (`cxx:304-308`).
        let mut hels = Arc::new(ElSpine::new());
        let (p_first, t_first, p_last, t_last, first_param, last_param) = {
            let mut guard = spine
                .write()
                .map_err(|_| "FilletSurf_InternalBuilder::Add: poisoned spine".to_string())?;
            let first_param = guard.first_parameter();
            let Some((p_first, t_first)) = guard.d1(first_param) else {
                return Err("FilletSurf_InternalBuilder::Add: Spine::D1 failed".to_string());
            };
            let Some(last_param) = guard.last_parameter() else {
                return Err("FilletSurf_InternalBuilder::Add: the spine is not loaded".to_string());
            };
            let Some((p_last, t_last)) = guard.d1(last_param) else {
                return Err("FilletSurf_InternalBuilder::Add: Spine::D1 failed".to_string());
            };
            (p_first, t_first, p_last, t_last, first_param, last_param)
        };
        // `hels->FirstParameter(sp->FirstParameter());
        // hels->SetFirstPointAndTgt(PFirst, TFirst);
        // hels->LastParameter(sp->LastParameter());
        // hels->SetLastPointAndTgt(PLast, TLast);` (`cxx:309-313`).
        {
            let es = Arc::make_mut(&mut hels);
            es.set_first_parameter(first_param);
            es.set_first_point_and_tgt(p_first, t_first);
            es.set_last_parameter(last_param);
            es.set_last_point_and_tgt(p_last, t_last);
        }
        // `ChFi3d_PerformElSpine(hels, sp, myConti, tolesp);
        // sp->AppendElSpine(hels); sp->SplitDone(true);` (`cxx:314-316`).
        ch_fi3d_perform_el_spine(
            &mut hels,
            &spine,
            self.base.base().continuity(),
            self.base.base().tolesp(),
            false,
        )?;
        {
            let mut guard = spine
                .write()
                .map_err(|_| "FilletSurf_InternalBuilder::Add: poisoned spine".to_string())?;
            guard.append_el_spine(hels);
            guard.set_split_done(true);
        }
        Ok(0)
    }

    /// `Perform()` (`cxx:322-338`): `PerformSetOfKGen` on the single contour.
    pub fn perform(&mut self) -> Result<(), String> {
        self.perform_set_of_k_gen(false)
    }

    /// `Simulate()` (`cxx:833-848`): the same walk in simulation mode.
    pub fn simulate(&mut self) -> Result<(), String> {
        self.perform_set_of_k_gen(true)
    }

    /// The shared body of `Perform` / `Simulate`.
    fn perform_set_of_k_gen(&mut self, simul: bool) -> Result<(), String> {
        // `occ::handle<ChFiDS_Stripe> Stripe = myListStripe.First();
        // occ::handle<NCollection_HSequence<...>>& HData =
        // Stripe->ChangeSetOfSurfData();
        // HData = new NCollection_HSequence<...>();` (`cxx:326-329`,
        // `:836-839`).
        let stripe = self
            .first_stripe()
            .ok_or_else(|| "FilletSurf_InternalBuilder: no contour".to_string())?;
        stripe
            .write()
            .map_err(|_| "FilletSurf_InternalBuilder: poisoned stripe".to_string())?
            .change_set_of_surf_data()
            .clear();

        // `occ::handle<ChFiDS_Spine>& Spine = Stripe->ChangeSpine();
        // StripeOrientations(Spine, RefOr1, RefOr2, RefChoix);
        // Stripe->OrientationOnFace1(RefOr1); Stripe->OrientationOnFace2(RefOr2);
        // Stripe->Choix(RefChoix);` (`cxx:329-336`, `:839-846`).
        let spine = {
            let guard = stripe
                .read()
                .map_err(|_| "FilletSurf_InternalBuilder: poisoned stripe".to_string())?;
            guard.spine().cloned()
        }
        .ok_or_else(|| "FilletSurf_InternalBuilder: the stripe has no spine".to_string())?;
        let (ref_or1, ref_or2, ref_choix) = self.base.base().stripe_orientations(&spine)?;
        {
            let mut guard = stripe
                .write()
                .map_err(|_| "FilletSurf_InternalBuilder: poisoned stripe".to_string())?;
            guard.set_orientation_on_face1(ref_or1);
            guard.set_orientation_on_face2(ref_or2);
            guard.set_choix(ref_choix);
        }

        // `PerformSetOfKGen(Stripe, false);` (`cxx:337`).
        let mut surf_ops = FilletSurfOps::new(self.base.fillet_shape());
        self.base
            .base_mut()
            .perform_set_of_k_gen(&mut surf_ops, &stripe, simul)
    }

    /// `Done()` (`cxx:598-601`): the inherited `ChFi3d_Builder::done` member.
    ///
    /// `ChFi3d_Builder::done` is written by `Compute` (`ChFi3d_Builder.cxx:233`
    /// and the tail of the same member) and by `Reset`; `FilletSurf::
    /// InternalBuilder::Perform` runs `PerformSetOfKGen` directly and never
    /// touches it, so the answer stays `false` unless a caller ran `Compute`
    /// on the same object. The port keeps that behaviour.
    pub fn done(&self) -> bool {
        self.base.base().is_done()
    }

    /// `NbSurface()` (`cxx:609-612`).
    pub fn nb_surface(&self) -> usize {
        match self.first_stripe().and_then(|s| s.read().ok().map(|g| g.set_of_surf_data().len())) {
            Some(len) => len,
            None => 0,
        }
    }

    /// `SurfaceFillet(Index)` (`cxx:619-625`): the NUBS surface of index
    /// `Index` (`1`-based).
    pub fn surface_fillet(&self, index: usize) -> Option<Arc<dyn GeomSurface>> {
        let isurf = self.surf_data(index)?.surf();
        self.base.base().ds().surface(isurf).surface().cloned()
    }

    /// `TolApp3d(Index)` (`cxx:632-637`).
    pub fn tol_app3d(&self, index: usize) -> Option<f64> {
        let isurf = self.surf_data(index)?.surf();
        Some(self.base.base().ds().surface(isurf).tolerance())
    }

    /// `SupportFace1(Index)` (`cxx:642-648`).
    pub fn support_face1(&self, index: usize) -> Option<Face> {
        let isurf = self.surf_data(index)?.index_of_s1();
        self.base
            .base()
            .ds()
            .shape(isurf, false)
            .and_then(Face::wrap)
    }

    /// `SupportFace2(Index)` (`cxx:653-659`).
    pub fn support_face2(&self, index: usize) -> Option<Face> {
        let isurf = self.surf_data(index)?.index_of_s2();
        self.base
            .base()
            .ds()
            .shape(isurf, false)
            .and_then(Face::wrap)
    }

    /// `CurveOnFace1(Index)` (`cxx:666-670`).
    pub fn curve_on_face1(&self, index: usize) -> Option<Arc<dyn GeomCurve>> {
        let icurv = self.surf_data(index)?.interference_on_s1().line_index();
        self.base.base().ds().curve(icurv).curve().cloned()
    }

    /// `CurveOnFace2(Index)` (`cxx:676-680`).
    pub fn curve_on_face2(&self, index: usize) -> Option<Arc<dyn GeomCurve>> {
        let icurv = self.surf_data(index)?.interference_on_s2().line_index();
        self.base.base().ds().curve(icurv).curve().cloned()
    }

    /// `PCurveOnFace1(Index)` (`cxx:687-690`).
    pub fn p_curve_on_face1(&self, index: usize) -> Option<Arc<dyn Curve2d>> {
        self.surf_data(index)?
            .interference_on_s1()
            .p_curve_on_face()
            .cloned()
    }

    /// `PCurve1OnFillet(Index)` (`cxx:696-699`).
    pub fn p_curve1_on_fillet(&self, index: usize) -> Option<Arc<dyn Curve2d>> {
        self.surf_data(index)?
            .interference_on_s1()
            .p_curve_on_surf()
            .cloned()
    }

    /// `PCurveOnFace2(Index)` (`cxx:705-708`).
    pub fn p_curve_on_face2(&self, index: usize) -> Option<Arc<dyn Curve2d>> {
        self.surf_data(index)?
            .interference_on_s2()
            .p_curve_on_face()
            .cloned()
    }

    /// `PCurve2OnFillet(Index)` (`cxx:714-717`).
    pub fn p_curve2_on_fillet(&self, index: usize) -> Option<Arc<dyn Curve2d>> {
        self.surf_data(index)?
            .interference_on_s2()
            .p_curve_on_surf()
            .cloned()
    }

    /// `FirstParameter()` (`cxx:725-742`): the parameter of the fillet on the
    /// first edge.
    pub fn first_parameter(&self) -> f64 {
        let Some(stripe) = self.first_stripe() else {
            return 0.0;
        };
        let (spine, p) = {
            let Ok(guard) = stripe.read() else {
                return 0.0;
            };
            let Some(spine) = guard.spine().cloned() else {
                return 0.0;
            };
            let Some(sd) = guard.set_of_surf_data().first() else {
                return 0.0;
            };
            (spine, sd.first_spine_param())
        };
        let periodic = match spine.read() {
            Ok(guard) => guard.is_periodic(),
            Err(_) => return 0.0,
        };
        // `ind = 1; if (sp->IsPeriodic()) { ind = sp->Index(p); }`
        // (`cxx:730-734`) — `Index(W)` has `Forward = true` by default.
        let ind = if periodic {
            match spine.read() {
                Ok(guard) => match guard.index(p, true) {
                    Some(ind) => ind,
                    None => return 0.0,
                },
                Err(_) => return 0.0,
            }
        } else {
            1
        };
        let mut ep = 0.0f64;
        if compute_edge_parameter(
            &spine,
            ind,
            p,
            &mut ep,
            self.base.base().tolapp3d(),
        ) {
            return ep;
        }
        0.0
    }

    /// `LastParameter()` (`cxx:748-764`): the parameter of the fillet on the
    /// last edge.
    pub fn last_parameter(&self) -> f64 {
        let Some(stripe) = self.first_stripe() else {
            return 0.0;
        };
        let nb_surface = self.nb_surface();
        let (spine, p, nb_edges) = {
            let Ok(guard) = stripe.read() else {
                return 0.0;
            };
            let Some(spine) = guard.spine().cloned() else {
                return 0.0;
            };
            let Some(sd) = nb_surface
                .checked_sub(1)
                .and_then(|k| guard.set_of_surf_data().get(k))
            else {
                return 0.0;
            };
            let nb_edges = match spine.read() {
                Ok(s) => s.nb_edges(),
                Err(_) => return 0.0,
            };
            (spine, sd.last_spine_param(), nb_edges)
        };
        let periodic = match spine.read() {
            Ok(guard) => guard.is_periodic(),
            Err(_) => return 0.0,
        };
        // `ind = sp->NbEdges(); if (sp->IsPeriodic()) { ind = sp->Index(p); }`
        // (`cxx:753-757`).
        let ind = if periodic {
            match spine.read() {
                Ok(guard) => match guard.index(p, true) {
                    Some(ind) => ind,
                    None => return 0.0,
                },
                Err(_) => return 0.0,
            }
        } else {
            nb_edges as i32
        };
        let mut ep = 0.0f64;
        if compute_edge_parameter(
            &spine,
            ind,
            p,
            &mut ep,
            self.base.base().tolapp3d(),
        ) {
            return ep;
        }
        0.0
    }

    /// `StartSectionStatus()` (`cxx:778-795`).
    pub fn start_section_status(&self) -> StatusType {
        let (ison_edge1, ison_edge2) = match self.surf_data(1) {
            Some(sd) => (
                sd.vertex_first_on_s1().is_on_arc(),
                sd.vertex_first_on_s2().is_on_arc(),
            ),
            None => (false, false),
        };
        status_of(ison_edge1, ison_edge2)
    }

    /// `EndSectionStatus()` (`cxx:807-826`).
    pub fn end_section_status(&self) -> StatusType {
        let nb = self.nb_surface();
        let (ison_edge1, ison_edge2) = match self.surf_data(nb) {
            Some(sd) => (
                sd.vertex_last_on_s1().is_on_arc(),
                sd.vertex_last_on_s2().is_on_arc(),
            ),
            None => (false, false),
        };
        status_of(ison_edge1, ison_edge2)
    }

    /// `NbSection(IndexSurf)` (`cxx:854-857`).
    pub fn nb_section(&self, index_surf: usize) -> usize {
        match self.base.sect(1, index_surf) {
            Some(sec) => sec.len(),
            None => 0,
        }
    }

    /// `Section(IndexSurf, IndexSec)` (`cxx:866-872`): the trimmed circular
    /// arc of the `IndexSec`-th section of `SurfaceFillet(IndexSurf)`.
    pub fn section(&self, index_surf: usize, index_sec: usize) -> Option<GeomTrimmedCurve> {
        let sec = self.base.sect(1, index_surf)?;
        let (circ, deb, fin): (GpCirc, f64, f64) = sec.value(index_sec as i32)?.get_circ();
        let circle: Arc<dyn GeomCurve> = Arc::new(GeomCircle::new(circ));
        Some(GeomTrimmedCurve::new(circle, deb, fin))
    }

    /// `myListStripe.First()`, with the null answer of an empty list.
    fn first_stripe(&self) -> Option<StripeHandle> {
        self.base.base().stripes().first().cloned()
    }

    /// `myListStripe.First()->SetOfSurfData()->Value(Index)`, `Index`
    /// `1`-based.
    fn surf_data(&self, index: usize) -> Option<Arc<SurfData>> {
        if index == 0 {
            return None;
        }
        let stripe = self.first_stripe()?;
        let guard = stripe.read().ok()?;
        guard.set_of_surf_data().get(index - 1).cloned()
    }
}

/// `static bool isinlist(const TopoDS_Shape& E,
/// const NCollection_List<TopoDS_Shape>& L)` (`cxx:59-70`).
fn is_in_list(e: &TopoShape, list: &[TopoShape]) -> bool {
    list.iter().any(|value| is_same(value, e))
}

/// The common tail of `StartSectionStatus` / `EndSectionStatus`
/// (`cxx:782-794`, `:809-825`).
fn status_of(ison_edge1: bool, ison_edge2: bool) -> StatusType {
    if ison_edge1 && ison_edge2 {
        StatusType::TwoExtremityOnEdge
    } else if !ison_edge1 && !ison_edge2 {
        StatusType::NoExtremityOnEdge
    } else {
        StatusType::OneExtremityOnEdge
    }
}
