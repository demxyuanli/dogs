//! The starting solution of `PerformSetOfSurfOnElSpine`.
//!
//! Source: `ChFi3d_Builder_2.cxx:2283-2412` (OCCT 8.0.0).

use super::*;

impl<'a> Walk<'a> {
    /// `cxx:2283-2385`: the starting solution and the section bounds.
    pub(super) fn start(&mut self) -> Result<(), String> {
        // `bool forward = true; bool Inside = false; double First = wf;
        // double Last = wl; bool Ok1 = true, Ok2 = true;
        // TopoDS_Vertex Vref;` (`cxx:2277-2284`) — already in the struct.
        if self.ref_sd.is_none() && self.raf.is_none() {
            // `Inside = true; StartSol(Stripe, HGuide, HS1, HS2, It1, It2,
            // pp1, pp2, First);` (`cxx:2288-2290`).
            self.inside = true;
            self.builder.start_sol_el_spine(
                &mut *self.ops,
                self.stripe,
                &*self.guide,
                &mut self.hs1,
                &mut self.hs2,
                self.it1,
                self.it2,
                &mut self.pp1,
                &mut self.pp2,
                &mut self.first,
            )?;
            self.last = self.wf;
            if self.guide.is_periodic() {
                // `Last = First - Guide.Period(); ...` (`cxx:2292-2306`).
                let period = self
                    .guide
                    .period()
                    .ok_or_else(|| "PerformSetOfSurfOnElSpine: no period".to_string())?;
                self.last = self.first - period;
                Arc::make_mut(self.guide).save_first_parameter();
                Arc::make_mut(self.guide).set_first_parameter(self.last);
                Arc::make_mut(self.guide).save_last_parameter();
                Arc::make_mut(self.guide).set_last_parameter(self.first * 1.1);
                if let Some(offset) = self.offset.as_mut() {
                    Arc::make_mut(offset).save_first_parameter();
                    Arc::make_mut(offset).set_first_parameter(self.last);
                    Arc::make_mut(offset).save_last_parameter();
                    Arc::make_mut(offset).set_last_parameter(self.first * 1.1);
                    self.flush_offset();
                }
            }
        } else {
            // `if (!Spine->IsPeriodic() && (wl - lastedlastp > -tolesp))
            // Vref = Spine->LastVertex();` (`cxx:2311-2314`).
            let (spine_periodic, last_vertex) = {
                let guard = self
                    .spine
                    .read()
                    .map_err(|_| "PerformSetOfSurfOnElSpine: poisoned spine".to_string())?;
                (guard.is_periodic(), guard.last_vertex())
            };
            if !spine_periodic && (self.wl - self.lastedlastp > -self.builder.tolesp) {
                self.vref = last_vertex;
            }
            if self.ref_sd.is_none() {
                // `if (!Spine->IsPeriodic() && (wf < tolesp))
                // Vref = Spine->FirstVertex(); ref = raf; forward = false;
                // First = wl; Last = Guide.FirstParameter();`
                // (`cxx:2316-2323`).
                let (spine_periodic, first_vertex) = {
                    let guard = self
                        .spine
                        .read()
                        .map_err(|_| "PerformSetOfSurfOnElSpine: poisoned spine".to_string())?;
                    (guard.is_periodic(), guard.first_vertex())
                };
                if !spine_periodic && (self.wf < self.builder.tolesp) {
                    self.vref = first_vertex;
                }
                self.ref_sd = self.raf.clone();
                self.forward = false;
                self.first = self.wl;
                self.last = self.guide.first_parameter();
            }

            // `Ok1 = StartSol(Spine, HS1, pp1, HC1, w1, ref, !forward, 1,
            // HSref1, HCref1, RecP1, RecS1, RecRst1, obstacleon1, HS3, pp3,
            // decroch1, Vref);` and the `2` twin (`cxx:2326-2360`).
            let ref_sd = self
                .ref_sd
                .clone()
                .ok_or_else(|| "PerformSetOfSurfOnElSpine: no reference SurfData".to_string())?;
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
            // `HC1.Nullify(); HC2.Nullify();` (`cxx:2361-2362`).
            self.hc1 = None;
            self.hc2 = None;

            if self.ok1 && self.ok2 {
                if self.forward {
                    self.set_guide_first(self.wf);
                } else {
                    self.set_guide_last(self.wl);
                }
            }
        }

        // `bool fini = false; bool complete = Inside;`
        // (`cxx:2386-2387`).
        self.complete = self.inside;

        // `if (!Guide.IsPeriodic()) { int indf = Spine->Index(wf); int indl =
        // Spine->Index(wl, false); if (Spine->IsPeriodic() && (indl < indf))
        // indl += nbed; nbed = indl - indf + 1; }` (`cxx:2388-2396`).
        if !self.guide.is_periodic() {
            let guard = self
                .spine
                .read()
                .map_err(|_| "PerformSetOfSurfOnElSpine: poisoned spine".to_string())?;
            let indf = guard
                .index(self.wf, true)
                .ok_or_else(|| "PerformSetOfSurfOnElSpine: Index(wf)".to_string())?;
            let mut indl = guard
                .index(self.wl, false)
                .ok_or_else(|| "PerformSetOfSurfOnElSpine: Index(wl)".to_string())?;
            if guard.is_periodic() && indl < indf {
                indl += self.nbed;
            }
            self.nbed = indl - indf + 1;
        }

        // `double bidf = wf, bidl = wl; if (!Spine->IsPeriodic()) { ... }`
        // (`cxx:2402-2410`).
        let mut bidf = self.wf;
        let mut bidl = self.wl;
        let (spine_periodic, spine_last_param) = {
            let guard = self
                .spine
                .read()
                .map_err(|_| "PerformSetOfSurfOnElSpine: poisoned spine".to_string())?;
            (
                guard.is_periodic(),
                guard.last_parameter_of(guard.nb_edges() as i32),
            )
        };
        if !spine_periodic {
            let last_param = spine_last_param
                .ok_or_else(|| "PerformSetOfSurfOnElSpine: empty spine".to_string())?;
            bidf = 0.0f64.max(self.wf);
            bidl = self.wl.min(last_param);
            // PMN 20/07/98 : attention in case if there is only extension.
            if (bidl - bidf) < 0.01 * last_param {
                bidf = self.wf;
                bidl = self.wl;
            }
        }
        // `double MaxStep = (bidl - bidf) * 0.05 / nbed;`
        // (`cxx:2411`).
        self.max_step = (bidl - bidf) * 0.05 / f64::from(self.nbed);
        Ok(())
    }
}
