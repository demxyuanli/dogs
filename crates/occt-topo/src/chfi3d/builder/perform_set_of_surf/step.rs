//! One turn of the `while (!fini)` loop of `PerformSetOfSurfOnElSpine`.
//!
//! Source: `ChFi3d_Builder_2.cxx:2416-2988` (OCCT 8.0.0).

use super::*;

impl<'a> Walk<'a> {
    /// One turn of the `while (!fini)` loop (`cxx:2416-2988`).
    #[allow(unused_assignments)]
    pub(super) fn step(&mut self) -> Result<StepOutcome, String> {
        // `Ok1 = true, Ok2 = true;` (`cxx:2419`).
        self.ok1 = true;
        self.ok2 = true;

        // `if (!Spine->IsPeriodic()) { ... }` (`cxx:2420-2464`).
        let (spine_periodic, first_on_same, first_tangent, last_on_same, last_tangent) = {
            let guard = self
                .spine
                .read()
                .map_err(|_| "PerformSetOfSurfOnElSpine: poisoned spine".to_string())?;
            (
                guard.is_periodic(),
                guard.first_status() == State::OnSame,
                guard.is_tangency_extremity(true),
                guard.last_status() == State::OnSame,
                guard.is_tangency_extremity(false),
            )
        };
        if !spine_periodic {
            if self.wf < self.builder.tolesp && self.complete == self.inside {
                self.intf = if first_on_same { 2 } else { 1 };
            }
            if first_tangent {
                self.intf = 4;
                self.set_guide_first(self.wf_sav);
            }
            if self.wl - self.lastedlastp > -self.builder.tolesp {
                self.intl = if last_on_same { 2 } else { 1 };
            }
            if last_tangent {
                self.intl = 4;
                self.set_guide_last(self.wl_sav);
            }
        }

        // `if (intf && !forward) Vref = Spine->FirstVertex(); if (intl &&
        // forward) Vref = Spine->LastVertex();` (`cxx:2465-2470`).
        if self.intf != 0 && !self.forward {
            self.vref = self
                .spine
                .read()
                .map_err(|_| "PerformSetOfSurfOnElSpine: poisoned spine".to_string())?
                .first_vertex();
        }
        if self.intl != 0 && self.forward {
            self.vref = self
                .spine
                .read()
                .map_err(|_| "PerformSetOfSurfOnElSpine: poisoned spine".to_string())?
                .last_vertex();
        }

        // `if (!ref.IsNull()) { Ok1 = StartSol(... 1 ...); Ok2 =
        // StartSol(... 2 ...); }` (`cxx:2471-2509`).
        if let Some(ref_sd) = self.ref_sd.clone() {
            self.ok1 = self.builder.start_sol_with_ref(
                Some(&self.spine),
                &mut self.hs1,
                &mut self.pp1,
                &mut self.hc1,
                &mut self.w1,
                &ref_sd,
                !self.forward,
                1,
                &mut self.hsref1,
                &mut self.hcref1,
                &mut self.rec_p1,
                &mut self.rec_s1,
                &mut self.rec_rst1,
                &mut self.obstacleon1,
                &mut self.hs3,
                &mut self.pp3,
                self.decroch1,
                self.vref.as_ref(),
            )?;
            self.ok2 = self.builder.start_sol_with_ref(
                Some(&self.spine),
                &mut self.hs2,
                &mut self.pp2,
                &mut self.hc2,
                &mut self.w2,
                &ref_sd,
                !self.forward,
                2,
                &mut self.hsref2,
                &mut self.hcref2,
                &mut self.rec_p2,
                &mut self.rec_s2,
                &mut self.rec_rst2,
                &mut self.obstacleon2,
                &mut self.hs4,
                &mut self.pp4,
                self.decroch2,
                self.vref.as_ref(),
            )?;
        }

        // No more connected faces: build the tangent plane to continue the path
        // till the output on the other face (`cxx:2511-2553`).
        if (!self.ok1 && self.hc1.is_none()) || (!self.ok2 && self.hc2.is_none()) {
            if (self.intf != 0 && !self.forward) || (self.intl != 0 && self.forward) {
                let ref_sd = self.ref_sd.as_ref().ok_or_else(|| {
                    "PerformSetOfSurfOnElSpine: no reference SurfData".to_string()
                })?;
                if !self.ok1 {
                    crate::chfi3d::builder_2_helpers::ch_fi3d_build_plane(
                        &self.builder.my_ds,
                        &mut self.hs1,
                        &mut self.pp1,
                        ref_sd,
                        !self.forward,
                        1,
                    )?;
                }
                if !self.ok2 {
                    crate::chfi3d::builder_2_helpers::ch_fi3d_build_plane(
                        &self.builder.my_ds,
                        &mut self.hs2,
                        &mut self.pp2,
                        ref_sd,
                        !self.forward,
                        2,
                    )?;
                }
                if self.intf != 0 {
                    self.intf = 5;
                } else if self.intl != 0 {
                    self.intl = 5;
                }
                if self.forward {
                    self.set_guide_first(self.wf);
                } else {
                    self.set_guide_last(self.wl);
                }
            } else {
                return Err("PerformSetOfSurfOnElSpine : Chaining is impossible.".into());
            }
        }

        // `It1->Initialize(HSon1); It2->Initialize(HSon2);`
        // (`cxx:2556-2560`).
        initialize_topol_tool(self.it1, &self.hs1)?;
        initialize_topol_tool(self.it2, &self.hs2)?;

        // `SD = new ChFiDS_SurfData(); NCollection_Sequence<handle> SeqSD;
        // SeqSD.Append(SD);` (`cxx:2563-2565`).
        let mut seq_sd: Vec<Arc<SurfData>> = vec![Arc::new(SurfData::new())];

        // `if (obstacleon1 && obstacleon2) { ... } else if (obstacleon1) { ...
        // } else if (obstacleon2) { ... } else { ... CallPerformSurf ... }`
        // (`cxx:2567-2843`).
        if self.obstacleon1 && self.obstacleon2 {
            self.both_obstacles(&mut seq_sd)?;
        } else if self.obstacleon1 {
            self.obstacle_on_s1(&mut seq_sd)?;
        } else if self.obstacleon2 {
            self.obstacle_on_s2(&mut seq_sd)?;
        } else {
            // `const handle(Adaptor3d_TopolTool)& aTT1 = It1;` — the same tool.
            let mut sd = SurfData::new();
            let mut surf1 = Some(self.hs1.clone());
            let mut surf2 = Some(self.hs2.clone());
            self.builder.call_perform_surf(
                &mut *self.ops,
                self.stripe,
                self.simul,
                &mut seq_sd,
                &mut sd,
                &mut *self.guide,
                Some(&self.spine),
                &self.hs1,
                self.hs3.as_ref(),
                &self.pp1,
                &self.pp3,
                self.it1,
                &self.hs2,
                self.hs4.as_ref(),
                &self.pp2,
                &self.pp4,
                self.it2,
                self.max_step,
                self.locfleche,
                self.builder.tolesp,
                &mut self.first,
                &mut self.last,
                self.inside,
                self.inside,
                self.forward,
                self.rec_s1,
                self.rec_s2,
                &mut self.soldep,
                &mut self.intf,
                &mut self.intl,
                &mut surf1,
                &mut surf2,
            )?;
            // `SimulSurf` fills `SD`; `PerformSurf` works on `SeqSD`
            // (`cxx:2571-2588`, `:2811-2841`).
            if self.simul {
                seq_sd[0] = Arc::new(sd);
            }
            // `HS1` / `HS2` are the `Surf1` / `Surf2` outputs of
            // `CallPerformSurf` (`cxx:818`).
            if let Some(surf1) = surf1 {
                self.hs1 = surf1;
            }
            if let Some(surf2) = surf2 {
                self.hs2 = surf2;
            }
            self.decroch1 = false;
            self.decroch2 = false;
        }

        // `if (!done) { ... } else { ... }` (`cxx:2845-2928`).
        if !self.builder.done {
            // A failure in a part of extension is not serious: stop there.
            if (!self.ok1 && !self.obstacleon1) || (!self.ok2 && !self.obstacleon2) {
                self.builder.done = true;
                self.inside = false;
                if self.forward {
                    self.intl = 1;
                } else {
                    self.intf = 1;
                }
            } else {
                self.set_spine_error(crate::chfi_ds::ErrorStatus::WalkingFailure)?;
                return Err("CallPerformSurf : Path failed!".into());
            }
        } else {
            self.splice(&mut seq_sd)?;
        }

        // `if (Inside) { Inside = false; Firstsov = First; if
        // (Guide.IsPeriodic()) { complete = false; wf =
        // Guide.FirstParameter(); wl = Guide.LastParameter(); } }`
        // (`cxx:2929-2939`).
        if self.inside {
            self.inside = false;
            self.firstsov = self.first;
            if self.guide.is_periodic() {
                self.complete = false;
                self.wf = self.guide.first_parameter();
                self.wl = self.guide.last_parameter();
            }
        }

        // `if (forward) { ... }` (`cxx:2940-2981`).
        let mut fini = false;
        if self.forward {
            fini = ((self.wl - self.last) <= 10.0 * self.builder.tolesp)
                || (self.intl != 0 && !(self.obstacleon1 || self.obstacleon2));
            if !fini && self.guide.is_periodic() {
                let period = self
                    .guide
                    .period()
                    .ok_or_else(|| "PerformSetOfSurfOnElSpine: no period".to_string())?;
                if (self.wl - self.last) < period * 1.0e-3 {
                    // It is tested if reframing of extremes is done at the same
                    // edge, which closes the loop (`cxx:2946-2965`).
                    let (thefirst, thelast) = {
                        let guard = self.stripe.read().map_err(|_| {
                            "PerformSetOfSurfOnElSpine: poisoned stripe".to_string()
                        })?;
                        let seq = guard.set_of_surf_data();
                        (seq.first().cloned(), seq.last().cloned())
                    };
                    let thefirst = thefirst.ok_or_else(|| {
                        "PerformSetOfSurfOnElSpine: empty SurfData sequence".to_string()
                    })?;
                    let thelast = thelast.ok_or_else(|| {
                        "PerformSetOfSurfOnElSpine: empty SurfData sequence".to_string()
                    })?;
                    if thefirst.vertex_first_on_s1().is_on_arc()
                        && thelast.vertex_last_on_s1().is_on_arc()
                    {
                        fini = match (
                            thefirst.vertex_first_on_s1().arc(),
                            thelast.vertex_last_on_s1().arc(),
                        ) {
                            (Some(a), Some(b)) => crate::topo_tools::is_same(a, b),
                            _ => false,
                        };
                    }
                    if !fini
                        && thefirst.vertex_first_on_s2().is_on_arc()
                        && thelast.vertex_last_on_s2().is_on_arc()
                    {
                        fini = match (
                            thefirst.vertex_first_on_s2().arc(),
                            thelast.vertex_last_on_s2().arc(),
                        ) {
                            (Some(a), Some(b)) => crate::topo_tools::is_same(a, b),
                            _ => false,
                        };
                    }
                    if fini {
                        // It is ended! (`cxx:2962-2965`): OCCT returns without
                        // restoring the guideline range.
                        return Ok(StepOutcome::EarlyReturn);
                    }
                }
            }
            if fini && self.complete {
                // Restart in the opposite direction (`cxx:2968-2975`).
                self.ref_sd = {
                    let guard = self.stripe.read().map_err(|_| {
                        "PerformSetOfSurfOnElSpine: poisoned stripe".to_string()
                    })?;
                    guard.set_of_surf_data().first().cloned()
                };
                if self.ref_sd.is_none() {
                    // `Sequence().First()` on an empty sequence throws in OCCT.
                    return Err(
                        "PerformSetOfSurfOnElSpine: empty SurfData sequence".to_string()
                    );
                }
                self.forward = false;
                fini = false;
                self.first = self.firstsov;
            } else {
                self.first = self.last;
                self.last = self.wl;
            }
        }

        // `if (!forward) { fini = ((First - wf) <= 10. * tolesp ||
        // (intf && !(obstacleon1 || obstacleon2))); complete = false;
        // Last = wf; }` (`cxx:2982-2987`).
        if !self.forward {
            fini = ((self.first - self.wf) <= 10.0 * self.builder.tolesp)
                || (self.intf != 0 && !(self.obstacleon1 || self.obstacleon2));
            self.complete = false;
            self.last = self.wf;
        }

        Ok(if fini {
            StepOutcome::Done
        } else {
            StepOutcome::Continue
        })
    }
}
