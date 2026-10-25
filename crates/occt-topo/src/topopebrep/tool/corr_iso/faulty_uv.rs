//! `EdgeOUTofBoundsUV` / `EdgeWithFaultyUV` of `TopOpeBRepTool_CORRISO`.
//!
//! Ported from `TopOpeBRepTool_CORRISO.cxx:507-1157`.

use super::*;

impl CorrIso {
    /// `EdgeOUTofBoundsUV(E, onU, tolx, parspE)` (`cxx:753-855`): does the
    /// pcurve of `E` leave the `Fref` window along `x`, and how far?
    ///
    /// Answers `(recadre, parspE)`. `parspE` starts at `-1.e7` (`cxx:759`) and
    /// is only set by the `SPLITEDGE` arm of an `x`-crossing non-iso pcurve
    /// (`cxx:826-828`).
    pub fn edge_out_of_bounds_uv(&self, e: &TopoShape, on_u: bool, tolx: f64) -> (i32, f64) {
        let mut parsp_e = -1.0e7;
        let mut recadre = INIT;

        let Some(c2df) = self.uv_rep(e) else {
            // `cxx:761-765`: `return false` from an `int` function.
            return (SPLITEDGE, parsp_e);
        };
        let (pc, f, l, tol) = c2df.pc();
        let Some(pc) = pc else {
            // Unreachable: `Init` refuses an edge whose editable p-curve is
            // null (`cxx:162-165`).
            return (INIT, parsp_e);
        };

        let xfirst = if on_u {
            self.gas.first_u_parameter()
        } else {
            self.gas.first_v_parameter()
        };
        let xlast = if on_u {
            self.gas.last_u_parameter()
        } else {
            self.gas.last_v_parameter()
        };
        let xperiod = if on_u { self.u_per } else { self.v_per };

        if let Some(iso) = uv_iso_curve(pc.as_ref()) {
            // `cxx:781-830`: the p-curve is an iso.
            let in_x = (on_u && iso.iso_u) || ((!on_u) && iso.iso_v);
            if in_x {
                let xpar = if on_u { iso.origin.x() } else { iso.origin.y() };
                let toosmall = xpar < xfirst - tolx;
                let tobig = xpar > xfirst + xperiod + tolx;
                if toosmall {
                    recadre = INCREASE;
                }
                if tobig {
                    recadre = DECREASE;
                }
                return (recadre, parsp_e);
            }
            let in_y = (on_u && iso.iso_v) || ((!on_u) && iso.iso_u);
            if in_y {
                // `cxx:810-828`: the p-curve describes `[minx, maxx]` across
                // the `x` window.
                let d2ddir = if on_u { iso.dir.y } else { iso.dir.x };
                let reverse = d2ddir < 0.0;
                let xfactor = if reverse { -1.0 } else { 1.0 };
                let max = if reverse { f } else { l };
                let min = if reverse { l } else { f };
                let maxuv = pc.d0(max);
                let minuv = pc.d0(min);

                let maxx = if on_u { maxuv.x() } else { maxuv.y() };
                let minx = if on_u { minuv.x() } else { minuv.y() };

                let (rec, maxsup) = fun_tool_recadre(minx, maxx, xfirst, xlast, tolx);
                recadre = rec;
                if recadre == SPLITEDGE {
                    let xbound = if maxsup { xperiod } else { 0.0 };
                    parsp_e = max - xfactor * (maxx - xbound);
                }
                return (recadre, parsp_e);
            }
        } else {
            // `cxx:831-853`: the p-curve is not an iso.
            let mut bn2d = BndBox2d::new();
            let tol_e = BRepTool::edge_tolerance(&Edge(e.clone()));
            let toladd = tol_e.max(tol);
            bnd_lib_add2d::add_geom2d_range(pc.as_ref(), f, l, toladd, &mut bn2d);
            let (umin, vmin, umax, vmax) = bn2d
                .get()
                .unwrap_or((1.0e7, 1.0e7, -1.0e7, -1.0e7));
            let xmin = if on_u { umin } else { vmin };
            let xmax = if on_u { umax } else { vmax };

            let (rec, _maxsup) = fun_tool_recadre(xmin, xmax, xfirst, xlast, tolx);
            recadre = rec;
            if recadre == SPLITEDGE {
                // `cxx:845-851` (`NYIxpu271198`): the intersection of the
                // p-curve with the `x`-iso seam is not computed, the arm
                // answers the `INIT` sentinel instead.
                return (INIT, parsp_e);
            }
            return (recadre, parsp_e);
        }
        (recadre, parsp_e)
    }

    /// `EdgesOUTofBoundsUV(EdsToCheck, onU, tolx, FyEds)` (`cxx:859-885`):
    /// `(edge, 1)` for an edge that must move up, `(edge, -1)` for down.
    pub fn edges_out_of_bounds_uv(
        &self,
        eds_to_check: &[TopoShape],
        on_u: bool,
        tolx: f64,
    ) -> Vec<(TopoShape, i32)> {
        let mut fy_eds: Vec<(TopoShape, i32)> = Vec::new();
        for e in eds_to_check {
            let (recadre, _sspar) = self.edge_out_of_bounds_uv(e, on_u, tolx);
            if recadre == SPLITEDGE {
                // `cxx:871-874`: `FUN_Raise()` prints in an `OCCT_DEBUG` build
                // and leaves the edge out of `FyEds`.
                continue;
            }
            if recadre == INCREASE {
                bind_fy(&mut fy_eds, e, INCREASE);
            }
            if recadre == DECREASE {
                bind_fy(&mut fy_eds, e, DECREASE);
            }
        }
        fy_eds
    }

    /// `EdgeWithFaultyUV(E, Ivfaulty)` (`cxx:889-1069`): `0` for a valid edge,
    /// `1` / `2` for the faulty bound, `3` for both.
    ///
    /// A bound is valid when another edge of the same vertex carries the same
    /// UV geometry at the *opposite* parameter index.
    pub fn edge_with_faulty_uv(&self, e: &TopoShape) -> i32 {
        let mut ivfaulty = 0;
        let mut tttol: f64 = 1.0e-8;
        let tttol_f = BRepTool::face_tolerance(&Face(self.f_ref.clone()));
        let tttuv_f = self.tol(1, tttol_f).max(self.tol(2, tttol_f));
        let tttol_e = BRepTool::edge_tolerance(&Edge(e.clone()));
        let tttuv_e = self.tol(1, tttol_e).max(self.tol(2, tttol_e));

        let o_e = e.orientation();
        if o_e == Orientation::Internal || o_e == Orientation::External {
            return 0;
        }

        let (ves1, ves2) = vertices(e);
        let closed = match (&ves1, &ves2) {
            (Some(a), Some(b)) => is_same_occurrence(a, b),
            _ => false,
        };
        if closed {
            // `cxx:915-922`: a closed edge is assumed valid.
            return 0;
        }

        let mut nfyv = 0;
        for iv_e in 1..=2 {
            let ve = if iv_e == 1 { ves1.clone() } else { ves2.clone() };
            let Some(ve) = ve else {
                continue;
            };

            let Some(parve) = par_e(iv_e, e) else {
                continue;
            };
            let Some(c2df) = self.uv_rep(e) else {
                return 0;
            };
            let Some(uvve) = uvf(parve, &c2df) else {
                return 0;
            };

            let tttol_ve = BRepTool::vertex_tolerance(&crate::shape::Vertex(ve.clone()));
            let tttuv_ve = self.tol(1, tttol_ve).max(self.tol(2, tttol_ve));

            let Some(loe) = self.v_eds.get(&shape_map_key(&ve)) else {
                return 0;
            };

            let mut ve_ok = false;
            for other in loe.iter() {
                let oe = other.orientation();
                if is_same_occurrence(other, e) {
                    continue;
                }
                if oe == Orientation::Internal || oe == Orientation::External {
                    continue;
                }
                let Some(a_c2df) = self.e_rep_2d.get(&shape_map_key(other)) else {
                    return 0;
                };

                let (o_ves1, o_ves2) = vertices(other);
                for ive in 1..=2 {
                    let o_ve = if ive == 1 { o_ves1.clone() } else { o_ves2.clone() };
                    let Some(o_ve) = o_ve else {
                        continue;
                    };
                    if !is_same_occurrence(&o_ve, &ve) {
                        continue;
                    }
                    let Some(pare) = par_e(ive, other) else {
                        continue;
                    };
                    let Some(uvve2) = uvf(pare, a_c2df) else {
                        continue;
                    };
                    if ive == iv_e {
                        continue;
                    }
                    let tttol_ve2 = BRepTool::vertex_tolerance(&crate::shape::Vertex(o_ve.clone()));
                    let tttuv_ve2 = self.tol(1, tttol_ve2).max(self.tol(2, tttol_ve2));

                    tttol = tttol.max(tttuv_f.max(
                        tttuv_e.max(tttuv_e.max(tttuv_ve2.max(tttuv_ve))),
                    ));
                    if uvve.is_equal(&uvve2, tttol) {
                        ve_ok = true;
                        break;
                    }
                }
                if ve_ok {
                    break;
                }
            }

            if !ve_ok {
                nfyv += 1;
                ivfaulty = iv_e;
            }
        }
        if nfyv == 2 {
            ivfaulty = 3;
        }
        ivfaulty
    }

    /// `EdgesWithFaultyUV(EdsToCheck, nfybounds, FyEds, stopatfirst)`
    /// (`cxx:1073-1138`).
    ///
    /// `nfybounds` filters on the number of faulty bounds: `1` and `2` require
    /// exactly that many, `3` requires any (`cxx:1112-1124`).
    pub fn edges_with_faulty_uv(
        &self,
        eds_to_check: &[TopoShape],
        nfybounds: i32,
        stopatfirst: bool,
    ) -> Vec<(TopoShape, i32)> {
        let mut fy_eds: Vec<(TopoShape, i32)> = Vec::new();
        for echk in eds_to_check {
            let ivfaulty = self.edge_with_faulty_uv(echk);
            if ivfaulty == 0 {
                continue;
            }
            let nfyv = if ivfaulty == 3 { 2 } else { 1 };

            let found = if nfybounds == 1 || nfybounds == 2 {
                nfyv == nfybounds
            } else if nfybounds == 3 {
                nfyv > 0
            } else {
                false
            };

            if found {
                bind_fy(&mut fy_eds, echk, ivfaulty);
                if stopatfirst {
                    return fy_eds;
                }
            }
        }
        fy_eds
    }

    /// `EdgeWithFaultyUV(EdsToCheck, nfybounds, fyE, Ifaulty)`
    /// (`cxx:1142-1157`): the first faulty edge of the list.
    pub fn edge_with_faulty_uv_first(
        &self,
        eds_to_check: &[TopoShape],
        nfybounds: i32,
    ) -> Option<(TopoShape, i32)> {
        let fy_eds = self.edges_with_faulty_uv(eds_to_check, nfybounds, true);
        fy_eds.into_iter().next()
    }

}
