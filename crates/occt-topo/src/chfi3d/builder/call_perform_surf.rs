//! `ChFi3d_Builder::CallPerformSurf` and `ChFi3d_Builder::StripeOrientations`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx` (OCCT
//! 8.0.0): `CallPerformSurf` (`:635-818`), `StripeOrientations` (`:826-857`).
//!
//! `CallPerformSurf` is the single place where the base class decides *which*
//! pair of faces and which `Choix` the subclass' `PerformSurf` / `SimulSurf`
//! receives, and where it retries once on the second-choice face. Both virtuals
//! are pure in OCCT, so the port drives them through
//! [`ChFi3dSurfOps`](super::surf_ops::ChFi3dSurfOps).
//!
//! Divergences that are not heuristics (each is also stated at the call site):
//! - `occ::handle<BRepAdaptor_Surface>` becomes `SurfaceAdaptor`; the "null
//!   handle" of `HS3` / `HS4` is `Option<&SurfaceAdaptor>`, and the two output
//!   handles `Surf1` / `Surf2` are `&mut Option<SurfaceAdaptor>` because OCCT
//!   declares them null and only assigns them at the end,
//! - `It1->Initialize((const occ::handle<Adaptor3d_Surface>&)HSon1)` becomes
//!   `TopolTool::initialize(surface, face)`: the port's adaptor carries its
//!   face, which is what `BRepAdaptor_Surface::Surface()` / `::Face()` pair up
//!   in OCCT. A null adaptor surface answers `Err` where OCCT would raise,
//! - `NCollection_Sequence<occ::handle<ChFiDS_SurfData>>&` is
//!   `&mut Vec<Arc<SurfData>>` (the shape of `ChFiDS_Stripe::SetOfSurfData`),
//! - `surf1` / `surf2` are written even when both attempts failed, as in OCCT.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_math::MathVector;

use crate::abs::Orientation;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi3d::builder_helpers::conex_faces;
use crate::chfi3d::side::{concave_side, next_side};
use crate::chfi_ds::{ElSpine, SpineHandle, StripeHandle, SurfData};
use crate::shape::{Edge, Face};
use crate::topo_tools::is_same;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::surf_ops::ChFi3dSurfOps;
use super::ChFi3dBuilder;

/// `It->Initialize((const occ::handle<Adaptor3d_Surface>&)HS)` on the adaptor's
/// own surface and face.
pub(super) fn initialize_topol_tool(
    it: &Rc<RefCell<TopolTool>>,
    hs: &SurfaceAdaptor,
) -> Result<(), String> {
    let surface = hs
        .surface()
        .cloned()
        .ok_or_else(|| "CallPerformSurf: null surface on the adaptor".to_string())?;
    let face = hs
        .face()
        .cloned()
        .ok_or_else(|| "CallPerformSurf: null face on the adaptor".to_string())?;
    it.borrow_mut().initialize(surface, &face)
}

/// `HS->Face().Orientation()`.
fn adaptor_orientation(hs: &SurfaceAdaptor) -> Result<Orientation, String> {
    hs.face()
        .map(|f| f.0.orientation())
        .ok_or_else(|| "CallPerformSurf: null face on the adaptor".to_string())
}

impl ChFi3dBuilder {
    /// `CallPerformSurf(Stripe, Simul, SeqSD, SD, HGuide, Spine, HS1, HS3, pp1,
    /// pp3, It1, HS2, HS4, pp2, pp4, It2, MaxStep, Fleche, TolGuide, First,
    /// Last, Inside, Appro, forward, RecOnS1, RecOnS2, Soldep, intf, intl,
    /// Surf1, Surf2)` (`cxx:635-818`).
    ///
    /// `ops` is the port of the two protected pure virtuals the call dispatches
    /// to ([`ChFi3dSurfOps`]). `Err` covers the failures OCCT raises
    /// (`TopLoc`-free `TopoDS::Face`, a null adaptor surface, `ChFi3d::NextSide`
    /// on an internal / external orientation).
    #[allow(clippy::too_many_arguments)]
    pub fn call_perform_surf(
        &mut self,
        ops: &mut dyn ChFi3dSurfOps,
        stripe: &StripeHandle,
        simul: bool,
        seq_sd: &mut Vec<Arc<SurfData>>,
        sd: &mut SurfData,
        hguide: &mut Arc<ElSpine>,
        spine: Option<&SpineHandle>,
        hs1: &SurfaceAdaptor,
        hs3: Option<&SurfaceAdaptor>,
        pp1: &GpPnt2d,
        pp3: &GpPnt2d,
        it1: &Rc<RefCell<TopolTool>>,
        hs2: &SurfaceAdaptor,
        hs4: Option<&SurfaceAdaptor>,
        pp2: &GpPnt2d,
        pp4: &GpPnt2d,
        it2: &Rc<RefCell<TopolTool>>,
        max_step: f64,
        fleche: f64,
        // `const double /*TolGuide*/` and `const bool /*Appro*/` are commented
        // out in OCCT (`cxx:653`, `:657`): the two virtuals receive `tolesp`
        // and `Inside` instead. Kept positionally, like the C++ declaration.
        _tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        _approx: bool,
        forward: bool,
        rec_on_s1: bool,
        rec_on_s2: bool,
        soldep: &mut MathVector,
        intf: &mut i32,
        intl: &mut i32,
        surf1: &mut Option<SurfaceAdaptor>,
        surf2: &mut Option<SurfaceAdaptor>,
    ) -> Result<(), String> {
        let (stripe_or1, stripe_or2, stripe_choix) = {
            let guard = stripe
                .read()
                .map_err(|_| "CallPerformSurf: poisoned stripe".to_string())?;
            (
                guard.orientation_on_face1(),
                guard.orientation_on_face2(),
                guard.choix(),
            )
        };

        let mut hson1 = hs1.clone();
        let mut hson2 = hs2.clone();
        // Definition of the domain of path It1, It2
        initialize_topol_tool(it1, &hson1)?;
        initialize_topol_tool(it2, &hson2)?;

        let mut or1 = adaptor_orientation(hs1)?;
        let mut or2 = adaptor_orientation(hs2)?;
        let mut choix = next_side(
            &mut or1,
            &mut or2,
            stripe_or1,
            stripe_or2,
            stripe_choix,
        );
        soldep.set_value(1, pp1.x());
        soldep.set_value(2, pp1.y());
        soldep.set_value(3, pp2.x());
        soldep.set_value(4, pp2.y());

        let thef = *first;
        let thel = *last;

        // OCCT passes `Inside` for the `Inside` *and* the `Appro` parameter of
        // both virtuals (`cxx:689-736`); the `/*Appro*/` of this function's own
        // parameter list is unused. Kept verbatim.
        let mut isdone = if simul {
            ops.simul_surf(
                sd,
                hguide,
                spine,
                choix,
                &hson1,
                it1,
                &hson2,
                it2,
                self.tolesp,
                first,
                last,
                inside,
                inside,
                forward,
                rec_on_s1,
                rec_on_s2,
                soldep,
                intf,
                intl,
                &self.ctx(),
            )
        } else {
            ops.perform_surf(
                self,
                seq_sd,
                hguide,
                spine,
                choix,
                &hson1,
                it1,
                &hson2,
                it2,
                max_step,
                fleche,
                self.tolesp,
                first,
                last,
                inside,
                inside,
                forward,
                rec_on_s1,
                rec_on_s2,
                soldep,
                intf,
                intl,
            )
        };

        // Case of error
        if !isdone {
            *first = thef;
            *last = thel;
            let mut reprise = false;
            if let Some(hs3) = hs3 {
                hson1 = hs3.clone();
                initialize_topol_tool(it1, &hson1)?;
                or1 = adaptor_orientation(hs3)?;
                soldep.set_value(1, pp3.x());
                soldep.set_value(2, pp3.y());
                reprise = true;
            } else if let Some(hs4) = hs4 {
                hson2 = hs4.clone();
                initialize_topol_tool(it2, &hson2)?;
                or2 = adaptor_orientation(hs4)?;
                soldep.set_value(3, pp4.x());
                soldep.set_value(4, pp4.y());
                reprise = true;
            }

            if reprise {
                choix = next_side(
                    &mut or1,
                    &mut or2,
                    stripe_or1,
                    stripe_or2,
                    stripe_choix,
                );
                isdone = if simul {
                    ops.simul_surf(
                        sd,
                        hguide,
                        spine,
                        choix,
                        &hson1,
                        it1,
                        &hson2,
                        it2,
                        self.tolesp,
                        first,
                        last,
                        inside,
                        inside,
                        forward,
                        rec_on_s1,
                        rec_on_s2,
                        soldep,
                        intf,
                        intl,
                        &self.ctx(),
                    )
                } else {
                    ops.perform_surf(
                        self,
                        seq_sd,
                        hguide,
                        spine,
                        choix,
                        &hson1,
                        it1,
                        &hson2,
                        it2,
                        max_step,
                        fleche,
                        self.tolesp,
                        first,
                        last,
                        inside,
                        inside,
                        forward,
                        rec_on_s1,
                        rec_on_s2,
                        soldep,
                        intf,
                        intl,
                    )
                };
            }
        }
        // `void` in OCCT: the final `isdone` is dropped.
        let _ = isdone;
        *surf1 = Some(hson1);
        *surf2 = Some(hson2);
        Ok(())
    }

    /// `StripeOrientations(Spine, Or1, Or2, ChoixConge)` (`cxx:826-857`): the
    /// reference orientation that fixes the concave face of the fillet, and the
    /// `Choix` of the first band.
    ///
    /// OCCT answers `bool` and always returns `true`; `Err` carries the failures
    /// it would raise (`Spine->Edges(1)` out of range, `ChFi3d::ConcaveSide`).
    pub fn stripe_orientations(
        &self,
        spine: &SpineHandle,
    ) -> Result<(Orientation, Orientation, i32), String> {
        let an_edge = {
            let guard = spine
                .read()
                .map_err(|_| "StripeOrientations: poisoned spine".to_string())?;
            guard.edges(1).cloned().ok_or_else(|| {
                "StripeOrientations: ChFiDS_Spine::Edges(1) out of range".to_string()
            })?
        };
        // `TopoDS::Face(myEdgeFirstFace(anEdge))`: a missing binding is the null
        // face, which `IsSame` then answers `false` for.
        let first_face = self.edge_first_face(&an_edge).and_then(|f| Face::wrap(f.clone()));
        let (ff1, ff2) = conex_faces(&Edge(an_edge.clone()), &self.my_ef_map);
        let (Some(mut ff1), Some(mut ff2)) = (ff1, ff2) else {
            return Err("StripeOrientations: ChFi3d_conexfaces found no face pair".to_string());
        };
        if let Some(first_face) = &first_face {
            if is_same(&ff2.0, &first_face.0) {
                std::mem::swap(&mut ff1, &mut ff2);
            }
        }
        let of1 = ff1.0.orientation();
        let f1 = Face(ff1.0.oriented(Orientation::Forward));
        let of2 = ff2.0.orientation();
        let f2 = Face(ff2.0.oriented(Orientation::Forward));
        let sb1 = SurfaceAdaptor::of_face(&f1);
        let sb2 = SurfaceAdaptor::of_face(&f2);

        let mut or1 = Orientation::Forward;
        let mut or2 = Orientation::Forward;
        let choix_conge = concave_side(&sb1, &sb2, &Edge(an_edge), &mut or1, &mut or2)?;
        // `TopAbs::Compose(Or1, Of1)`: `Of1` is the reference orientation, i.e.
        // the port's first argument.
        Ok((
            Orientation::compose(of1, or1),
            Orientation::compose(of2, or2),
            choix_conge,
        ))
    }
}
