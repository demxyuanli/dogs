//! The common-point edge limiting of `PerformMoreThreeCorner` and the free
//! border parameter.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:2198-2366`
//! (OCCT 8.0.0).
//!
//! Two blocks:
//! - `cxx:2198-2315` walks the sharp slots, evaluates the edge point at
//!   `p(ic, icplus)`, and pulls the two neighbouring common points in
//!   (`ChangeVertex(isfirst, 3 - jf)` for the `icplus` side,
//!   `ChangeVertex(isfirst, jf)` for the `icmoins` side). When the edge point
//!   is nearer to the vertex than both common points (`dS < d1 || dS < d2`) and
//!   the two common points are not the same point, the nearer one that lies on
//!   an arc steps the parameters back to `ParameterOnArc()`; otherwise the
//!   common point is only used when it is within its own tolerance of the edge
//!   point (`cp.Point().Distance(PE) < cp.Tolerance()`), when both neighbours
//!   share it, or when there is a single fillet (`nconges == 1`), and then the
//!   slot is flagged `ponctuel`. Both arms record the largest accepted
//!   tolerance below `1` in `tolcp`.
//! - `cxx:2318-2366` re-reads the common point of the slot next to a free
//!   border and, when its arc is that border, sets both `p` cells of the free
//!   edge slot to `ParameterOnArc()`; the second interference on the surface is
//!   tried when the first one does not match.
//!
//! Translation notes:
//! - `cp1` / `cp2` are `ChFiDS_CommonPoint` copies in OCCT (`cxx:2215`),
//!   assigned from the `ChangeVertex(...)` reference; the port reads the same
//!   cell through [`MoreThreeCorner::common_point_at`] and keeps the copy,
//! - `jfp` is `PerformMoreThreeCorner`'s function-level local again: the
//!   `icplus` side writes `3 - jf(icplus)` (`cxx:2218`, `:2248`, `:2284`),
//!   the `icmoins` side reads `jf(icmoins)` without writing it (`cxx:2230`,
//!   `:2263`, `:2303`),
//! - `tolapp = 1.e-4` is the function-level constant of `cxx:1313`, not a port
//!   side threshold,
//! - `trouve` is only a local switch in the second block (`cxx:2342`), and
//!   `CP1` there is re-read from the surface interference `2` after the first
//!   one missed.

use crate::brep_tool::BRepTool;
use crate::shape::{TopoShape, Vertex};
use crate::topo_tools::is_same;
use crate::topopebrep::tool::adaptor::CurveAdaptor;

use super::calcul::calcul_param;
use super::surf_index::indices;
use super::MoreThreeCorner;

/// `tolapp` (`cxx:1313`).
const TOLAPP: f64 = 1.0e-4;

impl MoreThreeCorner {
    /// `cxx:2198-2366`: the `tolcp` limiting and the free-border parameter.
    /// `nconges` is `PerformMoreThreeCorner`'s own argument (`cxx:1237`), read
    /// by the `nconges == 1` arms of `cxx:2276` and `:2292`.
    pub(super) fn limit_by_common_points(
        &mut self,
        vtx: &TopoShape,
        nconges: usize,
    ) -> Result<(), String> {
        let nedge = self.nedge as usize;
        let vertex = Vertex::wrap(vtx.clone())
            .ok_or_else(|| "PerformMoreThreeCorner: V1 is not a vertex".to_string())?;
        let sommet = BRepTool::vertex_point(&vertex);

        // `double tolcp = 0;` (`cxx:2201`) — function scope, read by the
        // `GeomPlate` section at `cxx:2875` / `:2889`.
        let mut tolcp = 0.0f64;

        // `cxx:2203-2316`.
        if !self.deuxconges {
            for ic in 0..nedge {
                if !self.sharp[ic] {
                    continue;
                }
                let (icplus, icmoins) = indices(nedge, ic);
                // `BRepAdaptor_Curve C(TopoDS::Edge(Evive.Value(ic)));
                //  PE = C.Value(p.Value(ic, icplus));` (`cxx:2210-2211`).
                let edge = self.edge_at(ic)?;
                let curve = CurveAdaptor::of_edge(&edge);
                let pe = curve.d0(self.p[ic][icplus]).ok_or_else(|| {
                    "PerformMoreThreeCorner: Evive has no 3d curve".to_string()
                })?;
                let d_s = pe.distance(&sommet);

                // `cxx:2214-2232`.
                let mut d1 = 0.0f64;
                let mut d2 = 0.0f64;
                let cp1 = if !self.sharp[icplus] {
                    let isfirst = self.sens[icplus] == 1;
                    self.jfp = 3 - self.jf[icplus];
                    let cp = self.common_point_at(icplus, isfirst, self.jfp)?;
                    d1 = cp.point().distance(&sommet);
                    Some(cp)
                } else {
                    None
                };
                let cp2 = if !self.sharp[icmoins] {
                    let isfirst = self.sens[icmoins] == 1;
                    let cp = self.common_point_at(icmoins, isfirst, self.jf[icmoins])?;
                    d2 = cp.point().distance(&sommet);
                    Some(cp)
                } else {
                    None
                };

                // `bool samecompoint = false; if (!sharp(icmoins) &&
                //  !sharp(icplus)) samecompoint = cp1.Point().Distance(
                //  cp2.Point()) < tolapp;` (`cxx:2232-2237`).
                let mut samecompoint = false;
                if let (Some(a), Some(b)) = (cp1.as_ref(), cp2.as_ref()) {
                    samecompoint = a.point().distance(b.point()) < TOLAPP;
                }

                if (d_s < d1 || d_s < d2) && !samecompoint {
                    // step back till Common Points.
                    if d2 < d1 && cp1.as_ref().is_some_and(|c| c.is_on_arc()) {
                        // cp1 is chosen (`cxx:2242-2256`).
                        let cp = cp1.as_ref().unwrap();
                        let param = cp.parameter_on_arc().ok_or_else(|| {
                            "PerformMoreThreeCorner: the common point has no arc parameter"
                                .to_string()
                        })?;
                        self.p[ic][icmoins] = param;
                        self.p[ic][icplus] = self.p[ic][icmoins];
                        let isfirst = self.sens[icplus] == 1;
                        self.jfp = 3 - self.jf[icplus];
                        let para = self.param_at(icplus, self.jfp, self.i[icplus][ic], isfirst)?;
                        self.p[icplus][ic] = para;
                        let tol = cp.tolerance();
                        if tol > tolcp && tol < 1.0 {
                            tolcp = tol;
                        }
                    } else if cp2.as_ref().is_some_and(|c| c.is_on_arc()) {
                        // cp2 is chosen (`cxx:2257-2268`).
                        let cp = cp2.as_ref().unwrap();
                        let param = cp.parameter_on_arc().ok_or_else(|| {
                            "PerformMoreThreeCorner: the common point has no arc parameter"
                                .to_string()
                        })?;
                        self.p[ic][icmoins] = param;
                        self.p[ic][icplus] = self.p[ic][icmoins];
                        let isfirst = self.sens[icmoins] == 1;
                        let para =
                            self.param_at(icmoins, self.jf[icmoins], self.i[icmoins][ic], isfirst)?;
                        self.p[icmoins][ic] = para;
                        let tol = cp.tolerance();
                        if tol > tolcp && tol < 1.0 {
                            tolcp = tol;
                        }
                    }
                } else {
                    // step back till Common Point only if it is very close.
                    // `cxx:2272-2292`.
                    if let Some(cp) = cp1.as_ref() {
                        if (cp.point().distance(&pe) < cp.tolerance()
                            || samecompoint
                            || nconges == 1)
                            && cp.is_on_arc()
                        {
                            let param = cp.parameter_on_arc().ok_or_else(|| {
                                "PerformMoreThreeCorner: the common point has no arc parameter"
                                    .to_string()
                            })?;
                            self.p[ic][icmoins] = param;
                            self.ponctuel[ic] = true;
                            self.p[ic][icplus] = self.p[ic][icmoins];
                            let isfirst = self.sens[icplus] == 1;
                            self.jfp = 3 - self.jf[icplus];
                            let para =
                                self.param_at(icplus, self.jfp, self.i[icplus][ic], isfirst)?;
                            self.p[icplus][ic] = para;
                            let tol = cp.tolerance();
                            if tol > tolcp && tol < 1.0 {
                                tolcp = tol;
                            }
                        }
                    }
                    // `cxx:2293-2312`.
                    if let Some(cp) = cp2.as_ref() {
                        if (cp.point().distance(&pe) < cp.tolerance()
                            || samecompoint
                            || nconges == 1)
                            && cp.is_on_arc()
                        {
                            let param = cp.parameter_on_arc().ok_or_else(|| {
                                "PerformMoreThreeCorner: the common point has no arc parameter"
                                    .to_string()
                            })?;
                            self.ponctuel[icmoins] = true;
                            self.p[ic][icmoins] = param;
                            self.p[ic][icplus] = self.p[ic][icmoins];
                            let isfirst = self.sens[icmoins] == 1;
                            let para = self.param_at(
                                icmoins,
                                self.jf[icmoins],
                                self.i[icmoins][ic],
                                isfirst,
                            )?;
                            self.p[icmoins][ic] = para;
                            let tol = cp.tolerance();
                            if tol > tolcp && tol < 1.0 {
                                tolcp = tol;
                            }
                        }
                    }
                }
            }
        }
        self.tolcp = tolcp;

        // `cxx:2318-2366`: in case of a free border the parameter
        // corresponding to the common point on the free edge is chosen.
        for ic in 0..nedge {
            let edge = self.edge_at(ic)?;
            if !self.is_free_border(&edge) {
                continue;
            }
            let (icplus, icmoins) = indices(nedge, ic);
            // `indic = libre(ic) ? icmoins : icplus;` (`cxx:2330-2337`).
            let indic = if self.libre[ic] { icmoins } else { icplus };
            if self.sharp[indic] {
                continue;
            }
            let isfirst = self.sens[indic] == 1;
            // `ChangeVertex(isfirst, 1)` (`cxx:2341`).
            let cp1 = self.common_point_at(indic, isfirst, 1)?;
            let mut trouve = false;
            if cp1.is_on_arc() {
                if let Some(arc) = cp1.arc() {
                    if is_same(arc, &edge.0) {
                        let param = cp1.parameter_on_arc().ok_or_else(|| {
                            "PerformMoreThreeCorner: the common point has no arc parameter"
                                .to_string()
                        })?;
                        self.p[ic][icmoins] = param;
                        self.p[ic][icplus] = param;
                        trouve = true;
                    }
                }
            }
            if !trouve {
                // `ChangeVertex(isfirst, 2)` (`cxx:2355`).
                let cp2 = self.common_point_at(indic, isfirst, 2)?;
                if cp2.is_on_arc() {
                    if let Some(arc) = cp2.arc() {
                        if is_same(arc, &edge.0) {
                            let param = cp2.parameter_on_arc().ok_or_else(|| {
                                "PerformMoreThreeCorner: the common point has no arc parameter"
                                    .to_string()
                            })?;
                            self.p[ic][icmoins] = param;
                            self.p[ic][icplus] = param;
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// `Calcul_Param(CD.Value(slot), jfposit, ind, isfirst, para)`
    /// (`cxx:2249`, `:2263`, `:2284`, `:2303`).
    fn param_at(
        &self,
        slot: usize,
        jfposit: i32,
        ind: i32,
        isfirst: bool,
    ) -> Result<f64, String> {
        let handle = self.cd[slot]
            .clone()
            .ok_or_else(|| "PerformMoreThreeCorner: CD is null".to_string())?;
        let guard = handle
            .read()
            .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
        calcul_param(&guard, jfposit, ind, isfirst)
    }
}
