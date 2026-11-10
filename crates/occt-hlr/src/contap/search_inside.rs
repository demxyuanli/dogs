//! `Contap_TheSearchInside` - `IntStart_SearchInside.gxx` instantiated with
//! `Contap_SurfFunction` / `Adaptor3d_Surface` / `Adaptor3d_TopolTool`
//! (`Contap_TheSearchInside_0.cxx:1-40`).
//!
//! Sources: `IntStart_SearchInside.gxx:27-312`, `IntStart_SearchInside.lxx:1-40`.
//!
//! Two `Perform` overloads:
//! - `gxx:42-265` samples the domain (`TheTopolTool::NbSamples(U/V)`,
//!   `SamplePoint`) and runs `math_FunctionSetRoot` around each sample plus the
//!   twelve corner points OCCT adds for tangent corners.
//! - `gxx:271-312` converges from one given `(UStart, VStart)`.
//!
//! Both append only interior points that are strictly inside the domain
//! (`TopAbs_IN`).

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::CONFUSION;
use occt_math::function_set_root::MathFunctionSetWithDerivatives;
use occt_math::{MathMatrix, MathVector};

use occt_topo::brep_top_adaptor_topol_tool::TopolTool;
use occt_topo::fclass2d::FaceState;
use occt_topo::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::iwalking::types::InteriorPoint;
use super::surf_function::ContapSurfFunction;

/// `IntStart_SearchInside` (`Contap_TheSearchInside.hxx:29-60`).
#[derive(Clone)]
pub struct ContapTheSearchInside {
    /// `done` (`hxx:56`).
    done: bool,
    /// `list` (`hxx:57`): the interior points, in insertion order.
    list: Vec<InteriorPoint>,
}

impl Default for ContapTheSearchInside {
    /// `IntStart_SearchInside()` (`gxx:27-30`).
    fn default() -> Self {
        Self::new()
    }
}

impl ContapTheSearchInside {
    /// `IntStart_SearchInside()` (`gxx:27-30`).
    pub fn new() -> Self {
        Self {
            done: false,
            list: Vec::new(),
        }
    }

    /// `Perform(Func, PS, T, Epsilon)` (`gxx:42-265`).
    pub fn perform(
        &mut self,
        func: &mut ContapSurfFunction,
        ps: &SurfaceAdaptor,
        topol: &mut TopolTool,
        epsilon: f64,
    ) {
        self.done = false;
        self.list.clear();

        // `gxx:50-62`: `Binf/Bsup` are also the bounds used by
        // `math_FunctionSetRoot::Perform`.
        let (umin, vmin) = (ps.first_u_parameter(), ps.first_v_parameter());
        let (umax, vmax) = (ps.last_u_parameter(), ps.last_v_parameter());
        let mut binf = MathVector::new(1, 2);
        let mut bsup = MathVector::new(1, 2);
        binf.set_value(1, umin);
        binf.set_value(2, vmin);
        bsup.set_value(1, umax);
        bsup.set_value(2, vmax);

        let nb_sample_u = topol.nb_samples_u();
        let nb_sample_v = topol.nb_samples_v();
        let nb_sample = topol.nb_samples();

        let mut du = bsup.value(1) - binf.value(1);
        let mut dv = bsup.value(2) - binf.value(2);
        du /= f64::from(nb_sample_u) * 0.5;
        dv /= f64::from(nb_sample_v) * 0.5;

        // `gxx:71-73`.
        let toler1 = ps.u_resolution(CONFUSION);
        let toler2 = ps.v_resolution(CONFUSION);

        // `gxx:77-85`: keep the interior points away from the restrictions.
        let mut max_toler1_toler2 = toler1;
        if toler2 > max_toler1_toler2 {
            max_toler1_toler2 = toler2;
        }
        max_toler1_toler2 *= 1000.0;
        if max_toler1_toler2 > du * 0.001 {
            max_toler1_toler2 = du * 0.001;
        }
        if max_toler1_toler2 > dv * 0.001 {
            max_toler1_toler2 = dv * 0.001;
        }

        func.set_surface(ps.clone());
        let tol = func.tolerance();

        let mut toler = MathVector::new(1, 2);
        toler.set_value(1, toler1);
        toler.set_value(2, toler2);
        // `gxx:89`: `math_FunctionSetRoot Rsnld(Func, toler)` (100 iterations).
        let mut rsnld = occt_math::MathFunctionSetRoot::new(&*func, &toler, 100);

        // `gxx:92-95`.
        let umin = umin + du * 0.01;
        let vmin = vmin + dv * 0.01;
        let umax = umax - du * 0.01;
        let vmax = vmax - dv * 0.01;

        let mut uvap = MathVector::new(1, 2);

        // `gxx:118`: `for (i = 1; i <= Nbsample + 12; i++)`.
        for i in 1..=(nb_sample + 12) {
            let mut s2d = GpPnt2d::new(0.0, 0.0);
            let mut nepastester = false;
            if i <= nb_sample {
                // `gxx:124-126`: `T->SamplePoint(i, s2d, s3d)`.
                if let Some((p2d, _p3d)) = topol.sample_point(i) {
                    s2d = p2d;
                }
                uvap.set_value(1, s2d.x());
                uvap.set_value(2, s2d.y());

                // `gxx:129-133`.
                let u1 = umin.max(uvap.value(1) - du);
                let v1 = vmin.max(uvap.value(2) - dv);
                let u2 = umax.min(uvap.value(1) + du);
                let v2 = vmax.min(uvap.value(2) + dv);
                binf.set_value(1, u1);
                binf.set_value(2, v1);
                bsup.set_value(1, u2);
                bsup.set_value(2, v2);

                // `gxx:136-143`: skip a sample whose function value is far
                // from the diagonal of its box.
                let pextrm1 = ps.value(u1, v1).unwrap_or_default();
                let pextrm2 = ps.value(u2, v2).unwrap_or_default();
                let mut valf = MathVector::new(1, 1);
                let mut df = MathMatrix::new(1, 1, 1, 2);
                func.values(&uvap, &mut valf, &mut df);
                let rvalf = valf.value(1);
                let dist_pp = pextrm1.square_distance(&pextrm2);
                if rvalf * rvalf > 3.0 * dist_pp {
                    nepastester = true;
                }
            } else {
                // `gxx:150-198`: the twelve extra corner points.
                let (su, sv) = match i - nb_sample {
                    1 => (umin + du * 0.02, vmin + dv * 0.02),
                    2 => (umax - du * 0.02, vmin + dv * 0.02),
                    3 => (umin + du * 0.02, vmax - dv * 0.02),
                    4 => (umax - du * 0.02, vmax - dv * 0.02),
                    5 => (umin + du * 0.02, vmin + dv * 0.02),
                    6 => (umax - du * 0.02, vmin + dv * 0.02),
                    7 => (umin + du * 0.02, vmax - dv * 0.02),
                    8 => (umax - du * 0.02, vmax - dv * 0.02),
                    9 => (umin + du * 0.005, vmin + dv * 0.005),
                    10 => (umax - du * 0.005, vmin + dv * 0.005),
                    11 => (umin + du * 0.005, vmax - dv * 0.005),
                    _ => (umax - du * 0.005, vmax - dv * 0.005),
                };
                s2d = GpPnt2d::new(su, sv);
                uvap.set_value(1, s2d.x());
                uvap.set_value(2, s2d.y());

                // `gxx:202-207`.
                binf.set_value(1, umin.max(uvap.value(1) - du));
                binf.set_value(2, vmin.max(uvap.value(2) - dv));
                bsup.set_value(1, umax.min(uvap.value(1) + du));
                bsup.set_value(2, vmax.min(uvap.value(2) + dv));
            }

            if nepastester {
                continue;
            }

            // `gxx:211-261`.
            rsnld.perform_with_bounds(&mut *func, &uvap, &binf, &bsup, false);
            if !rsnld.is_done() {
                continue;
            }
            if func.root().abs() > tol {
                continue;
            }
            if func.is_tangent() {
                continue;
            }
            let psol = *func.point();
            let root = rsnld.root();
            uvap.set_value(1, root.value(1));
            uvap.set_value(2, root.value(2));

            // `gxx:222-246`: is the point new?
            let mut testpnt = !self.list.is_empty();
            let mut j = 0usize;
            while testpnt {
                let ipj = &self.list[j];
                let pj = ipj.value();
                if (pj.x() - psol.x()).abs() <= epsilon
                    && (pj.y() - psol.y()).abs() <= epsilon
                    && (pj.z() - psol.z()).abs() <= epsilon
                    && (uvap.value(1) - ipj.u_parameter()).abs() <= toler1
                    && (uvap.value(2) - ipj.v_parameter()).abs() <= toler2
                {
                    testpnt = false;
                } else {
                    j += 1;
                    testpnt = j < self.list.len();
                }
            }
            let is_new = j >= self.list.len();

            if is_new {
                // `gxx:248-256`.
                let situ = topol.classify(
                    GpPnt2d::new(uvap.value(1), uvap.value(2)),
                    max_toler1_toler2,
                    false,
                );
                if situ == FaceState::In {
                    let d3d = func.direction3d().map(|v| *v).unwrap_or_else(|_| {
                        occt_core::gp::GpVec::new(0.0, 0.0, 0.0)
                    });
                    let d2d = func
                        .direction2d()
                        .map(|d| GpVec2d::new(d.x(), d.y()))
                        .unwrap_or_else(|_| GpVec2d::new(0.0, 0.0));
                    self.list
                        .push(InteriorPoint::new(psol, uvap.value(1), uvap.value(2), d3d, d2d));
                }
            }
        }

        self.done = true;
    }

    /// `Perform(Func, PS, UStart, VStart)` (`gxx:271-312`).
    pub fn perform_at(
        &mut self,
        func: &mut ContapSurfFunction,
        ps: &SurfaceAdaptor,
        u_start: f64,
        v_start: f64,
    ) {
        self.done = false;
        self.list.clear();

        let mut binf = MathVector::new(1, 2);
        let mut bsup = MathVector::new(1, 2);
        let mut toler = MathVector::new(1, 2);
        binf.set_value(1, ps.first_u_parameter());
        binf.set_value(2, ps.first_v_parameter());
        bsup.set_value(1, ps.last_u_parameter());
        bsup.set_value(2, ps.last_v_parameter());
        toler.set_value(1, ps.u_resolution(CONFUSION));
        toler.set_value(2, ps.v_resolution(CONFUSION));

        // `gxx:288-289`.
        if u_start - binf.value(1) > -toler.value(1)
            && u_start - bsup.value(1) < toler.value(1)
            && v_start - binf.value(2) > -toler.value(2)
            && v_start - bsup.value(2) < toler.value(2)
        {
            func.set_surface(ps.clone());
            let mut uvap = MathVector::new(1, 2);
            uvap.set_value(1, u_start);
            uvap.set_value(2, v_start);

            let mut rsnld = occt_math::MathFunctionSetRoot::new(&*func, &toler, 100);
            rsnld.perform_with_bounds(&mut *func, &uvap, &binf, &bsup, false);
            if rsnld.is_done() {
                let tol = func.tolerance();
                let valf = func.root();
                if valf.abs() <= tol && !func.is_tangent() {
                    let psol = *func.point();
                    let root = rsnld.root();
                    uvap.set_value(1, root.value(1));
                    uvap.set_value(2, root.value(2));
                    let d3d = func.direction3d().map(|v| *v).unwrap_or_else(|_| {
                        occt_core::gp::GpVec::new(0.0, 0.0, 0.0)
                    });
                    let d2d = func
                        .direction2d()
                        .map(|d| GpVec2d::new(d.x(), d.y()))
                        .unwrap_or_else(|_| GpVec2d::new(0.0, 0.0));
                    self.list
                        .push(InteriorPoint::new(psol, uvap.value(1), uvap.value(2), d3d, d2d));
                }
            }
        }

        self.done = true;
    }

    /// `IsDone()` (`lxx:19-22`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `NbPoints()` (`lxx:24-31`): `0` where OCCT raises `StdFail_NotDone`.
    pub fn nb_points(&self) -> usize {
        if !self.done {
            return 0;
        }
        self.list.len()
    }

    /// `Value(Index)` (`lxx:33-40`), 1-based; `None` where OCCT raises
    /// `Standard_OutOfRange`.
    pub fn value(&self, index: usize) -> Option<&InteriorPoint> {
        index.checked_sub(1).and_then(|i| self.list.get(i))
    }

    /// The whole sequence, 1-based indices dropped.
    pub fn points(&self) -> &[InteriorPoint] {
        &self.list
    }
}
