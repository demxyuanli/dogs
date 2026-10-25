//! The three obstacle arms of `PerformSetOfSurfOnElSpine`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:2567-2809`
//! (OCCT 8.0.0).
//!
//! When `StartSol` left one of the two contacts on an arc of the neighbouring
//! face (`obstacleon1` / `obstacleon2`), the section is a surface / curve or
//! curve / curve blend, i.e. one of the four *edge-contact* overloads of
//! `PerformSurf` / `SimulSurf` instead of the face / face pair. The `Choix` is
//! adjusted per case exactly as OCCT does:
//!
//! | case | `Choix` source | adjustment |
//! |---|---|---|
//! | both (`cxx:2568-2658`) | `NextSide(Or1, Or2, ..., Choix)` | `4` / `1` |
//! | `S1` (`cxx:2659-2738`) | `NextSide(..., -Choix)` | `+1` / `-1` |
//! | `S2` (`cxx:2739-2808`) | `NextSide(Or1, Or2, ..., Choix)` | none |
//!
//! `SoldepCC` carries the two guideline parameters (`w1`, `w2`); `SoldepCS`
//! carries the guideline parameter of the contact and the two parameters of
//! the other face (`cxx:2585-2586`, `:2676-2678`, `:2747-2749`).
//!
//! Divergences that are not heuristics (each is also stated at the call site):
//! - `occ::handle<BRepAdaptor_Curve2d>` is `Option<&Curve2dAdaptor>` on the
//!   trait, so the null `HC1` / `HC2` stay `None`,
//! - the `PerformSurf` / `SimulSurf` overloads are the throwing defaults of
//!   `ChFi3d_Builder_NotImp.cxx` (`surf_ops.rs`), so every arm propagates their
//!   `Err("... Not Implemented")` until the subclasses are ported,
//! - the `SD` handle of `cxx:2564` is `seq_sd[0]`; the `SimulSurf` arms reach
//!   it through `Arc::make_mut` (`SimulSurf` fills the one `SurfData`,
//!   `PerformSurf` appends to the sequence).

use std::sync::Arc;

use crate::chfi3d::side::next_side;

use super::Walk;

impl<'a> Walk<'a> {
    /// `obstacleon1 && obstacleon2` (`cxx:2568-2658`).
    pub(super) fn both_obstacles(
        &mut self,
        seq_sd: &mut Vec<Arc<crate::chfi_ds::SurfData>>,
    ) -> Result<(), String> {
        // `TopAbs_Orientation Or1 = HSref1->Face().Orientation(); Or2 =
        // HSref2->Face().Orientation();` (`cxx:2569-2570`).
        let mut or1 = self
            .hsref1
            .face()
            .map(|face| face.0.orientation())
            .ok_or_else(|| "PerformSetOfSurfOnElSpine: HSref1 without face".to_string())?;
        let mut or2 = self
            .hsref2
            .face()
            .map(|face| face.0.orientation())
            .ok_or_else(|| "PerformSetOfSurfOnElSpine: HSref2 without face".to_string())?;
        // `int Choix = ChFi3d::NextSide(Or1, Or2,
        // Stripe->OrientationOnFace1(), Stripe->OrientationOnFace2(),
        // Stripe->Choix());` (`cxx:2571-2576`).
        let (stripe_or1, stripe_or2, stripe_choix) = {
            let guard = self
                .stripe
                .read()
                .map_err(|_| "PerformSetOfSurfOnElSpine: poisoned stripe".to_string())?;
            (
                guard.orientation_on_face1(),
                guard.orientation_on_face2(),
                guard.choix(),
            )
        };
        let mut choix = next_side(&mut or1, &mut or2, stripe_or1, stripe_or2, stripe_choix);
        // `if (Choix % 2 == 0) Choix = 4; else Choix = 1;`
        // (`cxx:2579-2585`).
        choix = if choix % 2 == 0 { 4 } else { 1 };

        // `SoldepCC(1) = w1; SoldepCC(2) = w2;` (`cxx:2585-2586`).
        self.soldep_cc.set_value(1, self.w1);
        self.soldep_cc.set_value(2, self.w2);

        // `if (Simul) SimulSurf(SD, ...) else PerformSurf(SeqSD, ...)`
        // (`cxx:2587-2656`).
        let tolesp = self.builder.tolesp;
        if self.simul {
            self.ops.simul_surf_contact(
                Arc::make_mut(&mut seq_sd[0]),
                &*self.guide,
                Some(&self.spine),
                choix,
                &self.hs1,
                self.it1,
                self.hc1.as_ref(),
                &self.hsref1,
                &self.hcref1,
                &mut self.decroch1,
                or1,
                &self.hs2,
                self.it2,
                self.hc2.as_ref(),
                &self.hsref2,
                &self.hcref2,
                &mut self.decroch2,
                or2,
                self.locfleche,
                tolesp,
                &mut self.first,
                &mut self.last,
                self.inside,
                self.inside,
                self.forward,
                self.rec_p1,
                self.rec_rst1,
                self.rec_p2,
                self.rec_rst2,
                &self.soldep_cc,
            )?;
        } else {
            self.ops.perform_surf_contact(
                seq_sd,
                &*self.guide,
                Some(&self.spine),
                choix,
                &self.hs1,
                self.it1,
                self.hc1.as_ref(),
                &self.hsref1,
                &self.hcref1,
                &mut self.decroch1,
                or1,
                &self.hs2,
                self.it2,
                self.hc2.as_ref(),
                &self.hsref2,
                &self.hcref2,
                &mut self.decroch2,
                or2,
                self.max_step,
                self.locfleche,
                tolesp,
                &mut self.first,
                &mut self.last,
                self.inside,
                self.inside,
                self.forward,
                self.rec_p1,
                self.rec_rst1,
                self.rec_p2,
                self.rec_rst2,
                &self.soldep_cc,
            )?;
        }

        // `SD->ChangeIndexOfS1(DStr.AddShape(HS1->Face()));
        // SD->ChangeIndexOfS2(DStr.AddShape(HS2->Face()));`
        // (`cxx:2657-2658`).
        self.set_indices(Arc::make_mut(&mut seq_sd[0]))
    }

    /// `obstacleon1` (`cxx:2659-2738`).
    pub(super) fn obstacle_on_s1(
        &mut self,
        seq_sd: &mut Vec<Arc<crate::chfi_ds::SurfData>>,
    ) -> Result<(), String> {
        // `TopAbs_Orientation Or1 = HSref1->Face().Orientation(); Or2 =
        // HS2->Face().Orientation();` (`cxx:2660-2661`).
        let mut or1 = self
            .hsref1
            .face()
            .map(|face| face.0.orientation())
            .ok_or_else(|| "PerformSetOfSurfOnElSpine: HSref1 without face".to_string())?;
        let mut or2 = self
            .hs2
            .face()
            .map(|face| face.0.orientation())
            .ok_or_else(|| "PerformSetOfSurfOnElSpine: HS2 without face".to_string())?;
        // `int Choix = ChFi3d::NextSide(Or1, Or2,
        // Stripe->OrientationOnFace1(), Stripe->OrientationOnFace2(),
        // -Stripe->Choix());` (`cxx:2662-2667`).
        let (stripe_or1, stripe_or2, stripe_choix) = {
            let guard = self
                .stripe
                .read()
                .map_err(|_| "PerformSetOfSurfOnElSpine: poisoned stripe".to_string())?;
            (
                guard.orientation_on_face1(),
                guard.orientation_on_face2(),
                guard.choix(),
            )
        };
        let mut choix = next_side(&mut or1, &mut or2, stripe_or1, stripe_or2, -stripe_choix);
        // `if (Choix % 2 == 1) Choix++; else Choix--;`
        // (`cxx:2669-2675`).
        choix = if choix % 2 == 1 { choix + 1 } else { choix - 1 };

        // `SoldepCS(3) = w1; SoldepCS(1) = pp2.X(); SoldepCS(2) = pp2.Y();`
        // (`cxx:2676-2678`).
        self.soldep_cs.set_value(3, self.w1);
        self.soldep_cs.set_value(1, self.pp2.x());
        self.soldep_cs.set_value(2, self.pp2.y());

        // `if (Simul) SimulSurf(SD, ...) else PerformSurf(SeqSD, ...)`
        // (`cxx:2679-2736`).
        let tolesp = self.builder.tolesp;
        if self.simul {
            self.ops.simul_surf_on_s1(
                Arc::make_mut(&mut seq_sd[0]),
                &*self.guide,
                Some(&self.spine),
                choix,
                &self.hs1,
                self.it1,
                self.hc1.as_ref(),
                &self.hsref1,
                &self.hcref1,
                &mut self.decroch1,
                &self.hs2,
                self.it2,
                or2,
                self.locfleche,
                tolesp,
                &mut self.first,
                &mut self.last,
                self.inside,
                self.inside,
                self.forward,
                self.rec_p1,
                self.rec_s2,
                self.rec_rst1,
                &self.soldep_cs,
            )?;
        } else {
            self.ops.perform_surf_on_s1(
                seq_sd,
                &*self.guide,
                Some(&self.spine),
                choix,
                &self.hs1,
                self.it1,
                self.hc1.as_ref(),
                &self.hsref1,
                &self.hcref1,
                &mut self.decroch1,
                &self.hs2,
                self.it2,
                or2,
                self.max_step,
                self.locfleche,
                tolesp,
                &mut self.first,
                &mut self.last,
                self.inside,
                self.inside,
                self.forward,
                self.rec_p1,
                self.rec_s2,
                self.rec_rst1,
                &self.soldep_cs,
            )?;
        }

        // `SD->ChangeIndexOfS1(DStr.AddShape(HS1->Face()));
        // SD->ChangeIndexOfS2(DStr.AddShape(HS2->Face())); decroch2 = false;`
        // (`cxx:2735-2737`).
        self.set_indices(Arc::make_mut(&mut seq_sd[0]))?;
        self.decroch2 = false;
        Ok(())
    }

    /// `obstacleon2` (`cxx:2739-2808`).
    pub(super) fn obstacle_on_s2(
        &mut self,
        seq_sd: &mut Vec<Arc<crate::chfi_ds::SurfData>>,
    ) -> Result<(), String> {
        // `TopAbs_Orientation Or1 = HS1->Face().Orientation(); Or2 =
        // HSref2->Face().Orientation();` (`cxx:2740-2741`).
        let mut or1 = self
            .hs1
            .face()
            .map(|face| face.0.orientation())
            .ok_or_else(|| "PerformSetOfSurfOnElSpine: HS1 without face".to_string())?;
        let mut or2 = self
            .hsref2
            .face()
            .map(|face| face.0.orientation())
            .ok_or_else(|| "PerformSetOfSurfOnElSpine: HSref2 without face".to_string())?;
        // `int Choix = ChFi3d::NextSide(Or1, Or2,
        // Stripe->OrientationOnFace1(), Stripe->OrientationOnFace2(),
        // Stripe->Choix());` (`cxx:2742-2747`).
        let (stripe_or1, stripe_or2, stripe_choix) = {
            let guard = self
                .stripe
                .read()
                .map_err(|_| "PerformSetOfSurfOnElSpine: poisoned stripe".to_string())?;
            (
                guard.orientation_on_face1(),
                guard.orientation_on_face2(),
                guard.choix(),
            )
        };
        let choix = next_side(&mut or1, &mut or2, stripe_or1, stripe_or2, stripe_choix);

        // `SoldepCS(3) = w2; SoldepCS(1) = pp1.X(); SoldepCS(2) = pp1.Y();`
        // (`cxx:2747-2749`).
        self.soldep_cs.set_value(3, self.w2);
        self.soldep_cs.set_value(1, self.pp1.x());
        self.soldep_cs.set_value(2, self.pp1.y());

        // `if (Simul) SimulSurf(SD, ...) else PerformSurf(SeqSD, ...)`
        // (`cxx:2750-2806`).
        let tolesp = self.builder.tolesp;
        if self.simul {
            self.ops.simul_surf_on_s2(
                Arc::make_mut(&mut seq_sd[0]),
                &*self.guide,
                Some(&self.spine),
                choix,
                &self.hs1,
                self.it1,
                or1,
                &self.hs2,
                self.it2,
                self.hc2.as_ref(),
                &self.hsref2,
                &self.hcref2,
                &mut self.decroch2,
                self.locfleche,
                tolesp,
                &mut self.first,
                &mut self.last,
                self.inside,
                self.inside,
                self.forward,
                self.rec_p2,
                self.rec_s1,
                self.rec_rst2,
                &self.soldep_cs,
            )?;
        } else {
            self.ops.perform_surf_on_s2(
                seq_sd,
                &*self.guide,
                Some(&self.spine),
                choix,
                &self.hs1,
                self.it1,
                or1,
                &self.hs2,
                self.it2,
                self.hc2.as_ref(),
                &self.hsref2,
                &self.hcref2,
                &mut self.decroch2,
                self.max_step,
                self.locfleche,
                tolesp,
                &mut self.first,
                &mut self.last,
                self.inside,
                self.inside,
                self.forward,
                self.rec_p2,
                self.rec_s1,
                self.rec_rst2,
                &self.soldep_cs,
            )?;
        }

        // `SD->ChangeIndexOfS1(DStr.AddShape(HS1->Face()));
        // SD->ChangeIndexOfS2(DStr.AddShape(HS2->Face())); decroch1 = false;`
        // (`cxx:2807-2809`).
        self.set_indices(Arc::make_mut(&mut seq_sd[0]))?;
        self.decroch1 = false;
        Ok(())
    }
}
