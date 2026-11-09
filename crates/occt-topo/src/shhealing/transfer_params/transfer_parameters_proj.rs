//! `ShapeAnalysis_TransferParametersProj`.

use super::prelude::*;
use super::curve_range::*;
use super::TransferParameters;

/// `ShapeAnalysis_TransferParametersProj`
/// (`ShapeAnalysis_TransferParametersProj.cxx`).
#[derive(Clone)]
pub struct TransferParametersProj {
    base: TransferParameters,
    /// `myCurve` (`hxx:97`).
    curve: Option<Arc<dyn Curve>>,
    /// `myCurve2d` (`hxx:98`).
    curve2d: Option<Arc<dyn Curve2d>>,
    /// The surface behind `myAC3d` (`Proj.cxx:97`).
    surface: Option<Arc<dyn Surface>>,
    /// `myAC3d`'s window, `(f2d, l2d)` (`Proj.cxx:94`, `:99`).
    ac_first: f64,
    ac_last: f64,
    /// `myPrecision` (`hxx:102`) - `BRep_Tool::Tolerance(E)` (`Proj.cxx:74`).
    precision: f64,
    /// `myForceProj` (`hxx:104`).
    force_proj: bool,
    /// `myInitOK` (`hxx:105`).
    init_ok: bool,
}

impl Default for TransferParametersProj {
    fn default() -> Self {
        Self::new()
    }
}

impl TransferParametersProj {
    /// `ShapeAnalysis_TransferParametersProj()` (`Proj.cxx:49-54`):
    /// `myPrecision = 0`, `myMaxTolerance = 1`, `myForceProj = false`,
    /// `myInitOK = false`.
    pub fn new() -> Self {
        Self {
            base: TransferParameters::new(),
            curve: None,
            curve2d: None,
            surface: None,
            ac_first: 0.0,
            ac_last: 0.0,
            precision: 0.0,
            force_proj: false,
            init_ok: false,
        }
    }

    /// `ShapeAnalysis_TransferParametersProj(E, F)` (`Proj.cxx:59-66`).
    pub fn with_edge_face(e: &Edge, f: &Face) -> Self {
        let mut t = Self::new();
        t.init(e, f);
        t
    }

    /// `Init` (`Proj.cxx:69-103`).
    pub fn init(&mut self, e: &Edge, f: &Face) {
        self.init_ok = false; // `Proj.cxx:71`
        self.base.init(e, Some(f)); // `Proj.cxx:72`
        // `myEdge = E` is already done by the base `Init` (`cxx:52`).
        self.precision = BRepTool::edge_tolerance(e); // `Proj.cxx:74`

        // `myCurve = BRep_Tool::Curve(E, myFirst, myLast)` (`Proj.cxx:77`); the
        // base already read the same range at `cxx:54-56`.
        self.curve = BRepTool::edge_curve(e);
        if self.curve.is_none() {
            self.base.first = 0.0; // `Proj.cxx:79-82`
            self.base.last = 1.0;
            return;
        }
        if self.base.face.is_none() {
            return; // `Proj.cxx:84-87` (`F.IsNull()`)
        }
        // `sae.PCurve(E, F, myCurve2d, f2d, l2d, false)` (`Proj.cxx:91-92`).
        if let Some((c2d, f2d, l2d)) = curve_on_surface_oriented(e, f, false) {
            self.surface = BRepTool::face_surface(f); // `Proj.cxx:97`
            self.curve2d = Some(c2d);
            self.ac_first = f2d;
            self.ac_last = l2d;
            self.init_ok = true; // `Proj.cxx:102`
        }
    }

    /// `SetMaxTolerance` (`cxx:79-82`).
    pub fn set_max_tolerance(&mut self, maxtol: f64) {
        self.base.set_max_tolerance(maxtol);
    }

    /// `ForceProjection` (`Proj.cxx:555-558`).
    pub fn force_projection(&mut self) -> &mut bool {
        &mut self.force_proj
    }

    /// The `!myInitOK || (!myForceProj && myPrecision < myMaxTolerance &&
    /// BRep_Tool::SameParameter(myEdge))` gate repeated at `Proj.cxx:112-116`,
    /// `:186-189`, `:226-228`, `:292-294`, `:542-545`.
    fn use_linear(&self) -> bool {
        !self.init_ok
            || (!self.force_proj
                && self.precision < self.base.max_tolerance
                && BRepTool::same_parameter(&self.base.edge))
    }

    /// `myAC3d` (`Proj.cxx:95-101`): `surface(pcurve(t))` over `[f2d, l2d]`,
    /// `f2d` / `l2d` coming from the COS representation.
    fn ac3d(&self) -> Option<CurveOnSurface> {
        Some(CurveOnSurface::new(
            self.curve2d.clone()?,
            self.surface.clone()?,
            self.ac_first,
            self.ac_last,
        ))
    }

    /// `Adaptor3d_CurveOnSurface Ad1(AC2d(myCurve2d, First, Last), AdS)`
    /// (`Proj.cxx:200-201`, `:462-464`).
    fn ac3d_range(&self, first: f64, last: f64) -> Option<CurveOnSurface> {
        Some(CurveOnSurface::new(
            self.curve2d.clone()?,
            self.surface.clone()?,
            first,
            last,
        ))
    }

    /// `(FirstParameter(), LastParameter())` of `myAC3d` (`Proj.cxx:122-123`,
    /// `:232`, `:241-242`, `:332-336`).
    fn ac_range(&self) -> (f64, f64) {
        (self.ac_first, self.ac_last)
    }

    /// `myCurve->IsClosed()` (`Proj.cxx:145`). The port's substitute is the
    /// coincident-end test (see `shape_analysis_curve`), so a trimmed arc of a
    /// periodic basis is not treated as closed.
    fn curve_is_closed(&self) -> bool {
        let Some(c) = self.curve.as_ref() else {
            return false;
        };
        let (f, l) = (c.first_parameter(), c.last_parameter());
        f.is_finite() && l.is_finite() && c.d0(f).distance(&c.d0(l)) <= CONFUSION
    }

    /// `Perform(Knots, To2d)` (`Proj.cxx:107-176`).
    pub fn perform(&self, knots: &[f64], to2d: bool) -> Vec<f64> {
        if self.use_linear() {
            return self.base.perform(knots, to2d); // `Proj.cxx:112-116`
        }
        let len = knots.len();
        let preci = 2.0 * PCONFUSION; // `Proj.cxx:121`
        let (first, last) = if to2d {
            self.ac_range()
        } else {
            (self.base.first, self.base.last)
        }; // `Proj.cxx:122-123`
        let mut max_par = first; // `Proj.cxx:124`
        let last_par = last;
        let mut prev_par = max_par;
        let mut res_knots = Vec::with_capacity(len);
        for &knot in knots {
            // `Proj.cxx:130-142`
            let par = self.preform_segment(knot, to2d, prev_par, last_par);
            prev_par = par;
            if prev_par > last_par {
                prev_par -= preci; // `Proj.cxx:135-138`
            }
            res_knots.push(par);
            if par > max_par {
                max_par = par;
            }
        }
        // `Proj.cxx:145-160`: pdn correcting on periodic.
        if self.curve_is_closed() {
            for i in (0..len).rev() {
                if res_knots[i] < max_par {
                    let base_last = if to2d {
                        self.ac_last
                    } else {
                        self.curve.as_ref().map(|c| c.last_parameter()).unwrap_or(last)
                    };
                    res_knots[i] = base_last - (len - i - 1) as f64 * preci;
                } else {
                    break;
                }
            }
        }
        // `Proj.cxx:162-174`: pdn correction on range.
        for k in res_knots.iter_mut() {
            *k = k.max(first).min(last);
        }
        res_knots
    }

    /// `Perform(Knot, To2d)` (`Proj.cxx:222-252`).
    pub fn perform_param(&self, knot: f64, to2d: bool) -> f64 {
        if self.use_linear() {
            return self.base.perform_param(knot, to2d); // `Proj.cxx:226-228`
        }
        let (first, last) = if to2d {
            self.ac_range() // `Proj.cxx:232`
        } else {
            (self.base.first, self.base.last) // `Proj.cxx:236`
        };
        let res = self.preform_segment(knot, to2d, first, last);
        // `Proj.cxx:240-251`: pdn correction on range.
        res.max(first).min(last)
    }

    /// `PreformSegment` (`Proj.cxx:180-217`).
    fn preform_segment(&self, param: f64, to2d: bool, first: f64, last: f64) -> f64 {
        let lin_par = self.base.perform_param(param, to2d); // `Proj.cxx:185`
        if self.use_linear() {
            return lin_par; // `Proj.cxx:186-190`
        }
        let (proj_param, lin_dev, proj_dev);
        if to2d {
            // `Proj.cxx:198-203`.
            let (Some(curve), Some(ad1)) = (self.curve.as_ref(), self.ac3d_range(first, last))
            else {
                return lin_par;
            };
            let p1 = curve.d0(param); // `myCurve->Value(Param)`
            let projected = project_adaptor(&ad1, &p1, self.precision, true);
            proj_param = projected.param;
            proj_dev = projected.distance;
            lin_dev = p1.distance(&ad1.d0(lin_par));
        } else {
            // `Proj.cxx:205-210`.
            let (Some(curve), Some(ad1)) =
                (self.curve.as_ref(), self.ac3d_range(self.ac_first, self.ac_last))
            else {
                return lin_par;
            };
            let p1 = ad1.d0(param); // `myAC3d.Value(Param)`
            let projected = project_range(curve.as_ref(), &p1, self.precision, first, last, false);
            proj_param = projected.param;
            proj_dev = projected.distance;
            lin_dev = p1.distance(&curve.d0(lin_par));
        }
        // `Proj.cxx:212-215`.
        if lin_dev <= proj_dev || (lin_dev < self.precision && lin_dev <= 2.0 * proj_dev) {
            return lin_par;
        }
        proj_param
    }

    /// `TransferRange(newEdge, prevPar, currPar, Is2d)` (`Proj.cxx:285-536`).
    pub fn transfer_range(&self, new_edge: &mut Edge, prev_par: f64, curr_par: f64, is2d: bool) {
        if self.use_linear() {
            // `Proj.cxx:291-295`
            self.base.transfer_range(new_edge, prev_par, curr_par, is2d);
            return;
        }
        let reg = GeometryRegistry::global();
        let mut samerange = true; // `Proj.cxx:297`
        copy_ranges(new_edge, &self.base.edge, 0.0, 1.0); // `Proj.cxx:300`
        let preci = PCONFUSION; // `Proj.cxx:305`
        let (first_par, last_par) = if prev_par < curr_par {
            (prev_par, curr_par)
        } else {
            (curr_par, prev_par)
        }; // `Proj.cxx:307-314`

        let (p1, p2, alpha, beta);
        if is2d {
            // `Proj.cxx:317-337`.
            let Some(ad) = self.ac3d() else {
                return;
            };
            p1 = ad.d0(first_par); // `myAC3d.Value(firstPar)`
            if point_is_infinite(&p1) {
                reg.set_same_range(&new_edge.0, false);
                return;
            }
            p2 = ad.d0(last_par);
            if point_is_infinite(&p2) {
                reg.set_same_range(&new_edge.0, false);
                return;
            }
            let fact = self.ac_last - self.ac_first; // `Proj.cxx:332`
            let (mut a, mut b) = (0.0, 1.0);
            if fact > epsilon(self.ac_last) {
                a = (first_par - self.ac_first) / fact; // `Proj.cxx:335`
                b = (last_par - self.ac_first) / fact; // `Proj.cxx:336`
            }
            alpha = a;
            beta = b;
        } else {
            // `Proj.cxx:339-361`.
            let Some(curve) = self.curve.as_ref() else {
                return;
            };
            p1 = curve.d0(first_par); // `myCurve->Value(firstPar)`
            if point_is_infinite(&p1) {
                reg.set_same_range(&new_edge.0, false);
                return;
            }
            p2 = curve.d0(last_par);
            if point_is_infinite(&p2) {
                reg.set_same_range(&new_edge.0, false);
                return;
            }
            let fact = self.base.last - self.base.first; // `Proj.cxx:356`
            let (mut a, mut b) = (0.0, 1.0);
            if fact > epsilon(self.base.last) {
                a = (first_par - self.base.first) / fact; // `Proj.cxx:359`
                b = (last_par - self.base.first) / fact; // `Proj.cxx:360`
            }
            alpha = a;
            beta = b;
        }
        let use_linear_first = alpha < preci; // `Proj.cxx:362`
        let use_linear_last = 1.0 - beta < preci; // `Proj.cxx:363`

        // `Proj.cxx:367-376`: walk every representation of `newEdge`.
        let Some(mut to_geom) = reg.edge_geom(&new_edge.0) else {
            return;
        };

        // The `BRep_Curve3D` representation (`Proj.cxx:381-451`).
        {
            let (mut ppar1, mut ppar2) = if !is2d {
                (first_par, last_par) // `Proj.cxx:383-387`
            } else {
                let c3d = to_geom.curve.clone(); // `toGC->Curve3D()` (`Proj.cxx:390`)
                let first = to_geom.first; // `Proj.cxx:393-395`
                let last = to_geom.last;
                let len = last - first;
                let lin_first = first + alpha * len; // `Proj.cxx:405`
                let lin_last = first + beta * len; // `Proj.cxx:406`
                let gac = CurveRange::new(c3d.clone(), first, last); // `Proj.cxx:400`
                let proj1 = next_project(&gac, lin_first, &p1, self.precision); // `Proj.cxx:407`
                let proj2 = next_project(&gac, lin_last, &p2, self.precision); // `Proj.cxx:408`
                let use_linear = (proj1.param - proj2.param).abs() < preci; // `Proj.cxx:404`
                let pos1 = c3d.d0(lin_first); // `Proj.cxx:409-410`
                let pos2 = c3d.d0(lin_last);
                let d01 = pos1.distance(&p1); // `Proj.cxx:411-412`
                let d02 = pos2.distance(&p2);
                let (mut a, mut b) = (proj1.param, proj2.param);
                // `Proj.cxx:413-421`
                if use_linear_first
                    || use_linear
                    || d01 <= proj1.distance
                    || (d01 < self.precision && d01 <= 2.0 * proj1.distance)
                {
                    a = lin_first;
                }
                if use_linear_last
                    || use_linear
                    || d02 <= proj2.distance
                    || (d02 < self.precision && d02 <= 2.0 * proj2.distance)
                {
                    b = lin_last;
                }
                (a, b)
            };
            // `Proj.cxx:423-443`: order and nudge.
            if ppar1 > ppar2 {
                std::mem::swap(&mut ppar1, &mut ppar2);
            }
            if ppar2 - ppar1 < preci {
                if ppar1 - to_geom.first < preci {
                    ppar2 += 2.0 * preci;
                } else if to_geom.last - ppar2 < preci {
                    ppar1 -= 2.0 * preci;
                } else {
                    ppar1 -= preci;
                    ppar2 += preci;
                }
            }
            to_geom.first = ppar1; // `Proj.cxx:444` (`toGC->SetRange`)
            to_geom.last = ppar2;
            if ppar1 != first_par || ppar2 != last_par {
                samerange = false; // `Proj.cxx:447-450`
            }
        }

        // The `BRep_CurveOnSurface` representations (`Proj.cxx:452-531`).
        let keys: Vec<usize> = to_geom.pcurves.keys().copied().collect();
        for key in keys {
            let Some(c2d) = to_geom.pcurves.get(&key).and_then(|v| v.first().cloned()) else {
                continue;
            };
            let Some(surface) = self.surface.clone() else {
                continue;
            };
            let (first, last) = to_geom
                .pcurve_ranges
                .get(&key)
                .copied()
                .unwrap_or((to_geom.first, to_geom.last)); // `Proj.cxx:457-459`
            let len = last - first;
            let ad1 = CurveOnSurface::new(c2d.clone(), surface, first, last); // `Proj.cxx:462-464`
            let lin_first = first + alpha * len; // `Proj.cxx:471`
            let lin_last = first + beta * len; // `Proj.cxx:472`
            let proj1 = next_project(&ad1, lin_first, &p1, self.precision); // `Proj.cxx:473`
            let proj2 = next_project(&ad1, lin_last, &p2, self.precision); // `Proj.cxx:474`
            // `Proj.cxx:476-478`
            let is_first_on_end = (proj1.param - first) / len < PCONFUSION;
            let is_last_on_end = (last - proj2.param) / len < PCONFUSION;
            let use_linear = (proj1.param - proj2.param).abs() < PCONFUSION;
            let mut local_linear_first = use_linear_first;
            let mut local_linear_last = use_linear_last;
            if is_first_on_end && !local_linear_first {
                local_linear_first = true; // `Proj.cxx:479-482`
            }
            if is_last_on_end && !local_linear_last {
                local_linear_last = true; // `Proj.cxx:483-486`
            }
            let pos1 = ad1.d0(lin_first); // `Proj.cxx:488-489`
            let pos2 = ad1.d0(lin_last);
            let d01 = pos1.distance(&p1); // `Proj.cxx:490-491`
            let d02 = pos2.distance(&p2);
            let mut ppar1 = proj1.param;
            let mut ppar2 = proj2.param;
            // `Proj.cxx:492-499`
            if local_linear_first
                || use_linear
                || d01 <= proj1.distance
                || (d01 < self.precision && d01 <= 2.0 * proj1.distance)
            {
                ppar1 = lin_first;
            }
            if local_linear_last
                || use_linear
                || d02 <= proj2.distance
                || (d02 < self.precision && d02 <= 2.0 * proj2.distance)
            {
                ppar2 = lin_last;
            }
            if ppar1 > ppar2 {
                std::mem::swap(&mut ppar1, &mut ppar2); // `Proj.cxx:501-508`
            }
            ppar1 = correct_parameter(&c2d, ppar1); // `Proj.cxx:509-510`
            ppar2 = correct_parameter(&c2d, ppar2);
            if ppar2 - ppar1 < preci {
                // `Proj.cxx:511-524`
                if ppar1 - first < preci {
                    ppar2 += 2.0 * preci;
                } else if last - ppar2 < preci {
                    ppar1 -= 2.0 * preci;
                } else {
                    ppar1 -= preci;
                    ppar2 += preci;
                }
            }
            to_geom.pcurve_ranges.insert(key, (ppar1, ppar2)); // `Proj.cxx:527`
            if ppar1 != first_par || ppar2 != last_par {
                samerange = false; // `Proj.cxx:528-531`
            }
        }
        reg.set_edge(&new_edge.0, to_geom);
        reg.set_same_range(&new_edge.0, samerange); // `Proj.cxx:534`
    }

    /// `IsSameRange` (`Proj.cxx:539-551`).
    pub fn is_same_range(&self) -> bool {
        if self.use_linear() {
            self.base.is_same_range()
        } else {
            false // `Proj.cxx:549`
        }
    }
}
