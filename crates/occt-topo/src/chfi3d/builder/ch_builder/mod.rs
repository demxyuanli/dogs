//! `ChFi3d_ChBuilder` — the chamfer subclass of `ChFi3d_Builder`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_ChBuilder.cxx`,
//! `ChFi3d_ChBuilder_C2.cxx`, `ChFi3d_ChBuilder_C3.cxx` (OCCT 8.0.0). This
//! module is the **data layer**: the constructor, the three `Add` overloads,
//! the distance / distance-angle / mode API, `Simulate` / `NbSurf` / `Sect`,
//! `ResetContour` and the `ChFi3dCornerOps` entry points.
//!
//! Module layout:
//! ```text
//! statics        SearchCommonFaces (`cxx:67-93`), ExtentSpineOnCommonFace
//!                (`cxx:105-186`), ConexFaces (`cxx:2283-2318`)
//! extent_corner  ExtentOne/Two/ThreeCorner (`cxx:1953-2212`)
//! set_regul      SetRegul (`cxx:2218-2274`)
//! simul_k_part   SimulKPart (`cxx:655-716`)
//! simul_surf     `SimulSurf` on two faces (`cxx:717-1222`)
//! perform_surf   `PerformSurf` on two faces (`cxx:1618-1853`)
//! perform_first_section `PerformFirstSection` (`cxx:1319-1614`)
//! perform_three_corner `PerformThreeCorner` (`ChFi3d_ChBuilder_C3.cxx:202-1142`)
//! surf_ops       the `ChFi3dSurfOps` impl of `ChSurfOps`
//! ```
//!
//! The surface half of the subclass is complete: `PerformSurf` / `SimulSurf` /
//! `SimulKPart` / `PerformFirstSection` and `PerformThreeCorner` all have a
//! body, and the corner side (`Extent*Corner`, `SetRegul`,
//! `PerformTwoCorner`) is ported as well.
//!
//! `PerformThreeCorner` (`ChFi3d_ChBuilder_C3.cxx:202-1142`) lives in
//! [`perform_three_corner`], together with its three file statics (`CoPlanar`
//! `:64-88`, `BoundSurf` `:96-114`, `ComputeIntersection` `:131-197`).
//!
//! The three `SimulSurf` and the three `PerformSurf` edge-contact arms are
//! `throw Standard_Failure("... Not Implemented")` in OCCT itself
//! (`cxx:1223-1310`, `:1854-1940`), so they keep the trait defaults that throw
//! the same way. [`ChSurfOps`] answers the other entries; see [`surf_ops`].
//!
//! Divergences that are not heuristics:
//! - OCCT derives `ChFi3d_ChBuilder` from `ChFi3d_Builder`; the port keeps the
//!   base as a field (`base`) and the only subclass member (`myMode`,
//!   `hxx:409`) next to it. The corner bodies need `myEFMap` read-only, so the
//!   three `Extent*Corner` trait methods take `&ChFi3dBuilder`, like
//!   `SetRegul` already does,
//! - `occ::down_cast<ChFiDS_ChamfSpine>(Value(IC))` is
//!   [`ChFi3dChBuilder::chamf_spine`], which checks `SpineClass::Chamf`,
//! - `occ::handle<NCollection_HArray1<ChFiDS_CircSection>>` is
//!   [`CircSectionArray`], and a null handle is `None`,
//! - `throw Standard_DomainError("the face is not common to any of edges of the
//!   contour")` (`cxx:304`, `:419`, `:523`) is the `Err` of `set_dist` /
//!   `set_dists` / `set_dist_angle`,
//! - `ChFiDS_ChamfSpine::GetDist` / `Dists` / `GetDistAngle` raise on a spine
//!   of another kind; the port answers `None`, reported as `Err`.

mod extent_corner;
mod perform_first_section;
mod perform_surf;
mod perform_three_corner;
mod set_regul;
mod simul_k_part;
mod simul_surf;
mod statics;
mod surf_ops;

use std::sync::{Arc, RwLock};

use crate::abs::Orientation;
use crate::chfi3d::side::concave_side;
use crate::chfi_ds::{
    ChamfMethod, ChamfMode, CircSectionArray, Spine, SpineClass, SpineHandle, Stripe, StripeHandle,
};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::is_same;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::{ChFi3dBuilder, ChFi3dCornerOps, OneCornerExit};

/// `ChFi3d_ChBuilder` (`ChFi3d_ChBuilder.hxx:46-410`).
///
/// The base class half lives in `base`; the only member OCCT adds is `myMode`
/// (`hxx:409`).
pub struct ChFi3dChBuilder {
    /// `ChFi3d_Builder` — the base class sub-object.
    base: ChFi3dBuilder,
    /// `ChFiDS_ChamfMode myMode;` (`hxx:409`).
    mode: ChamfMode,
}

/// The corner-side carrier of the chamfer subclass.
///
/// OCCT's corner bodies are members of the same object as the base, so they
/// read `myEFMap` through `this`; the port hands the base in explicitly (see
/// the module note). It carries no state of its own: the chamfer's
/// `Extent*Corner` and `SetRegul` read only base members.
pub struct ChSurfOps;

impl ChFi3dChBuilder {
    /// `ChFi3d_ChBuilder(S, Ta)` (`cxx:188-193`): `ChFi3d_Builder(S, Ta)` then
    /// `myMode = ChFiDS_ClassicChamfer`.
    pub fn new(shape: &TopoShape, ta: f64) -> Self {
        Self {
            base: ChFi3dBuilder::new(shape, ta),
            mode: ChamfMode::ClassicChamfer,
        }
    }

    /// The `ChFi3d_Builder` sub-object.
    pub fn base(&self) -> &ChFi3dBuilder {
        &self.base
    }

    /// The mutable `ChFi3d_Builder` sub-object.
    pub fn base_mut(&mut self) -> &mut ChFi3dBuilder {
        &mut self.base
    }

    /// `SetMode(theMode)` (`cxx:539-542`).
    pub fn set_mode(&mut self, the_mode: ChamfMode) {
        self.mode = the_mode;
    }

    /// `Mode()` (`cxx:557-560`).
    pub fn mode(&self) -> ChamfMode {
        self.mode
    }

    /// `IsChamfer(IC)` (`cxx:546-554`): the kind of section of the `IC`-th
    /// contour.
    pub fn is_chamfer(&self, ic: usize) -> Result<ChamfMethod, String> {
        let spine = self.chamf_spine(ic)?;
        let guard = spine
            .read()
            .map_err(|_| "ChFi3d_ChBuilder: poisoned spine".to_string())?;
        Ok(guard.is_chamfer())
    }

    /// `Add(E)` (`cxx:201-222`): a new contour on `E` with a `ChFiDS_ChamfSpine`
    /// of the default distances.
    pub fn add(&mut self, e: &Edge) -> Result<(), String> {
        // `if (!Contains(E) && myEFMap.Contains(E))` (`cxx:205`).
        if self.base.contains(e) != 0 || !self.base.my_ef_map.contains(&e.0) {
            return Ok(());
        }
        // `TopoDS_Face dummy; ... Sp = new ChFiDS_ChamfSpine(tolesp);`
        // (`cxx:203-209`).
        let tolesp = self.base.tolesp;
        let spine: SpineHandle = Arc::new(RwLock::new(Spine::with_tolerance_chamf(tolesp)));
        let stripe: StripeHandle = Arc::new(RwLock::new(Stripe::new()));
        stripe
            .write()
            .map_err(|_| "Add: poisoned stripe".to_string())?
            .set_spine(spine.clone());

        // `TopoDS_Edge E_wnt = E; E_wnt.Orientation(TopAbs_FORWARD);
        // Spine->SetEdges(E_wnt);` (`cxx:212-214`).
        let mut e_wnt = e.0.clone();
        e_wnt.set_orientation(Orientation::Forward);
        spine
            .write()
            .map_err(|_| "Add: poisoned spine".to_string())?
            .set_edges(&e_wnt);

        // `if (PerformElement(Spine, -1, dummy)) { PerformExtremity(Spine);
        // Spine->Load(); myListStripe.Append(Stripe); }` (`cxx:215-221`).
        if self.base.perform_element(&spine, -1.0, None)? {
            self.base.perform_extremity(&spine)?;
            spine
                .write()
                .map_err(|_| "Add: poisoned spine".to_string())?
                .load();
            self.base.my_list_stripe.push(stripe);
        }
        Ok(())
    }

    /// `Add(Dis, E)` (`cxx:232-260`): `Add(E)` with `myMode`, then the
    /// symmetric distance.
    ///
    /// The statement order is OCCT's: `Load`, `myListStripe.Append`,
    /// `SetDist`, `PerformExtremity` — different from `Add(E)` above, which
    /// extends before appending.
    pub fn add_with_dist(&mut self, dis: f64, e: &Edge) -> Result<(), String> {
        if self.base.contains(e) != 0 || !self.base.my_ef_map.contains(&e.0) {
            return Ok(());
        }
        let tolesp = self.base.tolesp;
        let spine: SpineHandle = Arc::new(RwLock::new(Spine::with_tolerance_chamf(tolesp)));
        let stripe: StripeHandle = Arc::new(RwLock::new(Stripe::new()));
        stripe
            .write()
            .map_err(|_| "Add: poisoned stripe".to_string())?
            .set_spine(spine.clone());

        let mut e_wnt = e.0.clone();
        e_wnt.set_orientation(Orientation::Forward);
        // `Spine->SetMode(myMode); Spine->SetEdges(E_wnt);` (`cxx:247-249`).
        {
            let mut guard = spine
                .write()
                .map_err(|_| "Add: poisoned spine".to_string())?;
            guard.set_mode(self.mode);
            guard.set_edges(&e_wnt);
        }
        // `if (PerformElement(Spine, -1, dummy)) { Spine->Load();
        // myListStripe.Append(Stripe); Spine->SetDist(Dis);
        // PerformExtremity(Spine); }` (`cxx:250-258`).
        if self.base.perform_element(&spine, -1.0, None)? {
            spine
                .write()
                .map_err(|_| "Add: poisoned spine".to_string())?
                .load();
            self.base.my_list_stripe.push(stripe);
            spine
                .write()
                .map_err(|_| "Add: poisoned spine".to_string())?
                .set_dist(dis);
            self.base.perform_extremity(&spine)?;
        }
        Ok(())
    }

    /// `Add(Dis1, Dis2, E, F)` (`cxx:326-361`): two distances on the common
    /// face `F`.
    ///
    /// `Offset` is `-1` except for
    /// `ChFiDS_ConstThroatWithPenetrationChamfer`, where it is
    /// `min(Dis1, Dis2)` (`cxx:343-347`).
    pub fn add_with_two_dists(
        &mut self,
        dis1: f64,
        dis2: f64,
        e: &Edge,
        f: &Face,
    ) -> Result<(), String> {
        if self.base.contains(e) != 0 || !self.base.my_ef_map.contains(&e.0) {
            return Ok(());
        }
        let tolesp = self.base.tolesp;
        let spine: SpineHandle = Arc::new(RwLock::new(Spine::with_tolerance_chamf(tolesp)));
        let stripe: StripeHandle = Arc::new(RwLock::new(Stripe::new()));
        stripe
            .write()
            .map_err(|_| "Add: poisoned stripe".to_string())?
            .set_spine(spine.clone());

        let mut e_wnt = e.0.clone();
        e_wnt.set_orientation(Orientation::Forward);
        // `Spine->SetMode(myMode); double Offset = -1; if (myMode ==
        // ChFiDS_ConstThroatWithPenetrationChamfer) Offset = min(Dis1, Dis2);
        // Spine->SetEdges(E_wnt);` (`cxx:342-349`).
        let offset = if self.mode == ChamfMode::ConstThroatWithPenetrationChamfer {
            dis1.min(dis2)
        } else {
            -1.0
        };
        {
            let mut guard = spine
                .write()
                .map_err(|_| "Add: poisoned spine".to_string())?;
            guard.set_mode(self.mode);
            guard.set_edges(&e_wnt);
        }
        // `if (PerformElement(Spine, Offset, F)) { Spine->Load();
        // myListStripe.Append(Stripe); Spine->SetDists(Dis1, Dis2);
        // PerformExtremity(Spine); }` (`cxx:350-359`).
        if self.base.perform_element(&spine, offset, Some(f))? {
            spine
                .write()
                .map_err(|_| "Add: poisoned spine".to_string())?
                .load();
            self.base.my_list_stripe.push(stripe);
            spine
                .write()
                .map_err(|_| "Add: poisoned spine".to_string())?
                .set_dists(dis1, dis2);
            self.base.perform_extremity(&spine)?;
        }
        Ok(())
    }

    /// `AddDA(Dis1, Angle, E, F)` (`cxx:445-471`): one distance and an angle.
    ///
    /// Unlike the two other `Add` overloads this body does *not* call
    /// `Spine->SetMode(myMode)` (`cxx:453-461`).
    pub fn add_with_dist_angle(
        &mut self,
        dis: f64,
        angle: f64,
        e: &Edge,
        f: &Face,
    ) -> Result<(), String> {
        if self.base.contains(e) != 0 || !self.base.my_ef_map.contains(&e.0) {
            return Ok(());
        }
        let tolesp = self.base.tolesp;
        let spine: SpineHandle = Arc::new(RwLock::new(Spine::with_tolerance_chamf(tolesp)));
        let stripe: StripeHandle = Arc::new(RwLock::new(Stripe::new()));
        stripe
            .write()
            .map_err(|_| "Add: poisoned stripe".to_string())?
            .set_spine(spine.clone());

        let mut e_wnt = e.0.clone();
        e_wnt.set_orientation(Orientation::Forward);
        spine
            .write()
            .map_err(|_| "Add: poisoned spine".to_string())?
            .set_edges(&e_wnt);
        // `if (PerformElement(Spine, -1, F)) { Spine->Load();
        // myListStripe.Append(Stripe); Spine->SetDistAngle(Dis1, Angle);
        // PerformExtremity(Spine); }` (`cxx:461-469`).
        if self.base.perform_element(&spine, -1.0, Some(f))? {
            spine
                .write()
                .map_err(|_| "Add: poisoned spine".to_string())?
                .load();
            self.base.my_list_stripe.push(stripe);
            spine
                .write()
                .map_err(|_| "Add: poisoned spine".to_string())?
                .set_dist_angle(dis, angle);
            self.base.perform_extremity(&spine)?;
        }
        Ok(())
    }

    /// `SetDist(Dis, IC, F)` (`cxx:268-309`).
    ///
    /// The search loop walks the edges of the contour until one of its two
    /// common faces is `F`; `F1` / `F2` are swapped when `F` is the second one
    /// (`cxx:295-300`), which only affects the dead `BRepAdaptor_Surface`
    /// locals — `FirstF1` / `FirstF2` are written at `cxx:287-291` and never
    /// read by this body, and `Sb1` / `Sb2` are declared (`cxx:279`) and never
    /// initialised. Both are kept as noted dead code.
    pub fn set_dist(&mut self, dis: f64, ic: usize, f: &Face) -> Result<(), String> {
        if ic > self.base.nb_elements() {
            return Ok(());
        }
        let spine = self.chamf_spine(ic)?;
        let nb_edges = spine
            .read()
            .map_err(|_| "SetDist: poisoned spine".to_string())?
            .nb_edges();
        // `while ((i <= csp->NbEdges()) && (!Found))` (`cxx:282-292`).
        let mut found = false;
        let mut i = 1usize;
        while i <= nb_edges && !found {
            let edge = spine
                .read()
                .map_err(|_| "SetDist: poisoned spine".to_string())?
                .edges(i)
                .cloned()
                .ok_or_else(|| format!("SetDist: ChFiDS_Spine::Edges({i}) out of range"))?;
            let (f1, f2) = statics::search_common_faces(&self.base.my_ef_map, &edge)?;
            found = f1.as_ref().is_some_and(|x| is_same(&x.0, &f.0))
                || f2.as_ref().is_some_and(|x| is_same(&x.0, &f.0));
            i += 1;
        }
        if !found {
            // `throw Standard_DomainError("the face is not common to any of
            // edges of the contour");` (`cxx:304`).
            return Err(
                "SetDist: the face is not common to any of edges of the contour".to_string(),
            );
        }
        // `csp->SetDist(Dis);` (`cxx:301`).
        spine
            .write()
            .map_err(|_| "SetDist: poisoned spine".to_string())?
            .set_dist(dis);
        Ok(())
    }

    /// `GetDist(IC)` (`cxx:312-316`).
    pub fn get_dist(&self, ic: usize) -> Result<f64, String> {
        let spine = self.chamf_spine(ic)?;
        let guard = spine
            .read()
            .map_err(|_| "GetDist: poisoned spine".to_string())?;
        guard
            .get_dist()
            .ok_or_else(|| "GetDist: ChFiDS_ChamfSpine::GetDist on another kind".to_string())
    }

    /// `SetDists(Dis1, Dis2, IC, F)` (`cxx:368-425`).
    ///
    /// The `ConcaveSide` of the found edge and the one of `Edges(1)` decide the
    /// order of the two distances (`cxx:403-415`).
    pub fn set_dists(
        &mut self,
        dis1: f64,
        dis2: f64,
        ic: usize,
        f: &Face,
    ) -> Result<(), String> {
        if ic > self.base.nb_elements() {
            return Ok(());
        }
        let spine = self.chamf_spine(ic)?;
        let nb_edges = spine
            .read()
            .map_err(|_| "SetDists: poisoned spine".to_string())?
            .nb_edges();
        // `while ((i <= csp->NbEdges()) && (!Found)) { SearchCommonFaces(...);
        // if (i == 1) { FirstF1 = F1; FirstF2 = F2; } Found = ...; i++; }`
        // (`cxx:385-396`).
        let mut first_f: Option<(Face, Face)> = None;
        let mut found: Option<(Face, Face, TopoShape, usize)> = None;
        let mut i = 1usize;
        while i <= nb_edges && found.is_none() {
            let edge = spine
                .read()
                .map_err(|_| "SetDists: poisoned spine".to_string())?
                .edges(i)
                .cloned()
                .ok_or_else(|| format!("SetDists: ChFiDS_Spine::Edges({i}) out of range"))?;
            let (f1, f2) = statics::search_common_faces(&self.base.my_ef_map, &edge)?;
            let (Some(f1), Some(f2)) = (f1, f2) else {
                // OCCT leaves a null face here and hands it to
                // `BRepAdaptor_Surface::Initialize`; the port stops.
                return Err("SetDists: a spine edge has no two common faces".to_string());
            };
            if i == 1 {
                first_f = Some((f1.clone(), f2.clone()));
            }
            if is_same(&f1.0, &f.0) || is_same(&f2.0, &f.0) {
                found = Some((f1, f2, edge, i));
            }
            i += 1;
        }

        // `if (Found) { if (F2.IsSame(F)) { F2 = F1; F1 = F; } ... }`
        // (`cxx:394-416`).
        let Some((mut f1, mut f2, edge, index)) = found else {
            return Err(
                "SetDists: the face is not common to any of edges of the contour".to_string(),
            );
        };
        if is_same(&f2.0, &f.0) {
            f2 = f1;
            f1 = f.clone();
        }
        // `Sb1.Initialize(F1); Sb2.Initialize(F2); Choix =
        // ChFi3d::ConcaveSide(Sb1, Sb2, csp->Edges(i - 1), Or1, Or2);`
        // (`cxx:401-405`).
        let mut or1 = Orientation::Forward;
        let mut or2 = Orientation::Forward;
        let choix = concave_side(
            &SurfaceAdaptor::of_face(&f1),
            &SurfaceAdaptor::of_face(&f2),
            &Edge(edge),
            &mut or1,
            &mut or2,
        )?;
        // `Sb1.Initialize(FirstF1); Sb2.Initialize(FirstF2); ChoixConge =
        // ChFi3d::ConcaveSide(Sb1, Sb2, csp->Edges(1), Or1, Or2);`
        // (`cxx:406-409`).
        let (first_f1, first_f2) = first_f.ok_or_else(|| {
            format!("SetDists: the contour of index {index} has no first edge")
        })?;
        let first_edge = spine
            .read()
            .map_err(|_| "SetDists: poisoned spine".to_string())?
            .edges(1)
            .cloned()
            .ok_or_else(|| "SetDists: ChFiDS_Spine::Edges(1) out of range".to_string())?;
        let choix_conge = concave_side(
            &SurfaceAdaptor::of_face(&first_f1),
            &SurfaceAdaptor::of_face(&first_f2),
            &Edge(first_edge),
            &mut or1,
            &mut or2,
        )?;
        // `if (ChoixConge % 2 != Choix % 2) csp->SetDists(Dis2, Dis1); else
        // csp->SetDists(Dis1, Dis2);` (`cxx:410-416`).
        let mut guard = spine
            .write()
            .map_err(|_| "SetDists: poisoned spine".to_string())?;
        if choix_conge % 2 != choix % 2 {
            guard.set_dists(dis2, dis1);
        } else {
            guard.set_dists(dis1, dis2);
        }
        Ok(())
    }

    /// `Dists(IC)` (`cxx:428-435`).
    pub fn dists(&self, ic: usize) -> Result<(f64, f64), String> {
        let spine = self.chamf_spine(ic)?;
        let guard = spine
            .read()
            .map_err(|_| "Dists: poisoned spine".to_string())?;
        guard
            .dists()
            .ok_or_else(|| "Dists: ChFiDS_ChamfSpine::Dists on another kind".to_string())
    }

    /// `SetDistAngle(Dis, Angle, IC, F)` (`cxx:480-526`).
    ///
    /// `Sb1` / `Sb2` are initialised twice (`cxx:514-518`) and never read — the
    /// body only locates `F` and calls `SetDistAngle`. Kept as noted dead code.
    pub fn set_dist_angle(
        &mut self,
        dis: f64,
        angle: f64,
        ic: usize,
        f: &Face,
    ) -> Result<(), String> {
        if ic > self.base.nb_elements() {
            return Ok(());
        }
        let spine = self.chamf_spine(ic)?;
        let nb_edges = spine
            .read()
            .map_err(|_| "SetDistAngle: poisoned spine".to_string())?
            .nb_edges();
        // `while ((i <= csp->NbEdges()) && (!Found)) { SearchCommonFaces(...);
        // if (i == 1) { FirstF1 = F1; FirstF2 = F2; } Found = ...; i++; }`
        // (`cxx:496-507`).
        let mut found = false;
        let mut i = 1usize;
        while i <= nb_edges && !found {
            let edge = spine
                .read()
                .map_err(|_| "SetDistAngle: poisoned spine".to_string())?
                .edges(i)
                .cloned()
                .ok_or_else(|| format!("SetDistAngle: ChFiDS_Spine::Edges({i}) out of range"))?;
            let (f1, f2) = statics::search_common_faces(&self.base.my_ef_map, &edge)?;
            found = f1.as_ref().is_some_and(|x| is_same(&x.0, &f.0))
                || f2.as_ref().is_some_and(|x| is_same(&x.0, &f.0));
            i += 1;
        }
        if !found {
            // `throw Standard_DomainError("the face is not common to any edges
            // of the contour");` (`cxx:523`).
            return Err(
                "SetDistAngle: the face is not common to any edges of the contour".to_string(),
            );
        }
        // `csp->SetDistAngle(Dis, Angle);` (`cxx:519`).
        spine
            .write()
            .map_err(|_| "SetDistAngle: poisoned spine".to_string())?
            .set_dist_angle(dis, angle);
        Ok(())
    }

    /// `GetDistAngle(IC)` (`cxx:530-535`).
    pub fn get_dist_angle(&self, ic: usize) -> Result<(f64, f64), String> {
        let spine = self.chamf_spine(ic)?;
        let guard = spine
            .read()
            .map_err(|_| "GetDistAngle: poisoned spine".to_string())?;
        guard.get_dist_angle().ok_or_else(|| {
            "GetDistAngle: ChFiDS_ChamfSpine::GetDistAngle on another kind".to_string()
        })
    }

    /// `ResetContour(IC)` (`cxx:564-571`): `chsp->Reset(true)`.
    pub fn reset_contour(&mut self, ic: usize) -> Result<(), String> {
        if ic > self.base.nb_elements() {
            return Ok(());
        }
        let spine = self.chamf_spine(ic)?;
        let mut guard = spine
            .write()
            .map_err(|_| "ResetContour: poisoned spine".to_string())?;
        guard
            .reset(true)
            .ok_or_else(|| "ResetContour: ChFiDS_ChamfSpine::Reset".to_string())
    }

    /// `Simulate(IC)` (`cxx:576-608`): `PerformSetOfSurf(stripe, true)` on the
    /// `IC`-th stripe, dispatched to the [`ChSurfOps`] carrier — OCCT reaches
    /// the same object through `this`.
    ///
    /// The `#ifdef OCCT_DEBUG` chronometer block (`cxx:579-586`, `:597-606`)
    /// is not translated.
    pub fn simulate(&mut self, ic: usize) -> Result<(), String> {
        let some = ic
            .checked_sub(1)
            .and_then(|k| self.base.my_list_stripe.get(k).cloned());
        if let Some(stripe) = some {
            let mut ops = ChSurfOps;
            self.base.perform_set_of_surf(&mut ops, &stripe, true)?;
        }
        Ok(())
    }

    /// `NbSurf(IC)` (`cxx:614-625`): the number of `SurfData` of the `IC`-th
    /// stripe, `0` when there is no such stripe.
    pub fn nb_surf(&self, ic: usize) -> usize {
        let some = ic
            .checked_sub(1)
            .and_then(|k| self.base.my_list_stripe.get(k))
            .and_then(|s| s.read().ok());
        match some {
            Some(guard) => guard.set_of_surf_data().len(),
            None => 0,
        }
    }

    /// `Sect(IC, IS)` (`cxx:631-646`): the simulated section of the `IS`-th
    /// `SurfData` of the `IC`-th stripe; `None` is OCCT's null handle.
    pub fn sect(&self, ic: usize, is: usize) -> Option<Arc<CircSectionArray>> {
        let stripe = ic
            .checked_sub(1)
            .and_then(|k| self.base.my_list_stripe.get(k))?;
        let guard = stripe.read().ok()?;
        let sd = guard.set_of_surf_data().get(is.checked_sub(1)?)?;
        sd.simul().cloned()
    }

    /// `ChFi3d_Builder::Compute()` reached on this subclass
    /// (`BRepFilletAPI_MakeChamfer.cxx:277`: `myBuilder.Compute()`).
    ///
    /// Translation boundary: the port's base entry point
    /// [`ChFi3dBuilder::compute`] takes the two groups of `ChFi3d_Builder`
    /// virtuals as `&mut dyn` parameters, because Rust cannot alias one object
    /// into two `&mut` references. [`ChSurfOps`] carries no state, so the two
    /// carriers built here are interchangeable with the single `this` OCCT
    /// passes implicitly.
    pub fn compute(&mut self) -> Result<(), String> {
        let mut corner_ops = ChSurfOps;
        let mut surf_ops = ChSurfOps;
        self.base.compute(&mut corner_ops, &mut surf_ops)
    }

    /// The spine of the `IC`-th stripe as a `ChFiDS_ChamfSpine`, the
    /// `occ::down_cast<ChFiDS_ChamfSpine>(Value(IC))` of the data layer.
    fn chamf_spine(&self, ic: usize) -> Result<SpineHandle, String> {
        let spine = self
            .base
            .value(ic)
            .ok_or_else(|| format!("ChFi3d_ChBuilder: Value({ic}) has no spine"))?;
        let is_chamf = spine
            .read()
            .map_err(|_| "ChFi3d_ChBuilder: poisoned spine".to_string())?
            .class()
            == SpineClass::Chamf;
        if !is_chamf {
            return Err(format!(
                "ChFi3d_ChBuilder: Value({ic}) is not a ChFiDS_ChamfSpine"
            ));
        }
        Ok(spine)
    }
}

impl ChFi3dCornerOps for ChSurfOps {
    /// `SetRegul()` (`cxx:2218-2274`) — see [`set_regul`].
    fn set_regul(&mut self, base: &ChFi3dBuilder) -> Result<(), String> {
        set_regul::set_regul(base)
    }

    /// `PerformTwoCorner(Index)` (`ChFi3d_ChBuilder_C2.cxx:22-27`): the body is
    /// `PerformTwoCornerbyInter(Index);`, the base member.
    fn perform_two_corner(
        &mut self,
        base: &mut ChFi3dBuilder,
        index: usize,
    ) -> Result<(), String> {
        base.perform_two_corner_by_inter(index)
    }

    /// `PerformThreeCorner(Index)` (`ChFi3d_ChBuilder_C3.cxx:202-1142`) — see
    /// [`perform_three_corner::perform_three_corner`].
    fn perform_three_corner(
        &mut self,
        base: &mut ChFi3dBuilder,
        index: usize,
    ) -> Result<(), String> {
        perform_three_corner::perform_three_corner(base, index)
    }

    /// `ExtentOneCorner(V, S)` (`cxx:1953-1976`) — see [`extent_corner`].
    fn extent_one_corner(
        &mut self,
        _base: &ChFi3dBuilder,
        v: &TopoShape,
        s: StripeHandle,
    ) -> Result<(), String> {
        extent_corner::extent_one_corner(v, &s)
    }

    /// `ExtentTwoCorner(V, LS)` (`cxx:1982-2107`).
    fn extent_two_corner(
        &mut self,
        base: &ChFi3dBuilder,
        v: &TopoShape,
        ls: &[StripeHandle],
    ) -> Result<(), String> {
        extent_corner::extent_two_corner(&base.my_ef_map, v, ls)
    }

    /// `ExtentThreeCorner(V, LS)` (`cxx:2111-2212`).
    fn extent_three_corner(
        &mut self,
        base: &ChFi3dBuilder,
        v: &TopoShape,
        ls: &[StripeHandle],
    ) -> Result<(), String> {
        extent_corner::extent_three_corner(&base.my_ef_map, v, ls)
    }

    /// `PerformOneCorner(Index, PrepareOnSame)` (`ChFi3d_Builder_C1.cxx:611-1647`),
    /// the base member the chamfer subclass does not override.
    fn perform_one_corner(
        &mut self,
        base: &mut ChFi3dBuilder,
        index: usize,
        prepare_on_same: bool,
    ) -> Result<(), String> {
        match base.perform_one_corner(index, prepare_on_same)? {
            OneCornerExit::Done => Ok(()),
            OneCornerExit::IntersectionAtEnd => base.dispatch_end_intersection_exit(index),
        }
    }

    /// `PerformMoreSurfdata(Index)` (`ChFi3d_Builder_C1.cxx:3771-4531`), the
    /// base member the chamfer subclass does not override.
    fn perform_more_surfdata(
        &mut self,
        base: &mut ChFi3dBuilder,
        index: usize,
    ) -> Result<(), String> {
        base.perform_more_surfdata(index)
    }
}
