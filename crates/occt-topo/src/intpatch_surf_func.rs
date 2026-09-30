//! `IntPatch_TheSurfFunction` (`IntImp_ZerImpFunc.gxx`).
//! F(u,v) = `IntSurf_Quadric::Distance`(S(u,v)).
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::gp::{GpDir2d, GpPnt, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::Surface;

use crate::intpatch::impimp::{distance, gradient, ImplicitQuad};

const EPS_ANG2: f64 = 1.0e-16;
const TOL_PETIT: f64 = 1.0e-16;

/// `IntPatch_TheSurfFunction`.
pub(crate) struct SurfFunction<'a> {
    surf: &'a dyn Surface,
    quad: &'a ImplicitQuad,
    tol: f64,
    u: f64,
    v: f64,
    pntsol: GpPnt,
    valf: f64,
    computed: bool,
    tangent: bool,
    derived: bool,
    tgdu: f64,
    tgdv: f64,
    gradient: GpVec,
    d1u: GpVec,
    d1v: GpVec,
    d3d: GpVec,
    d2d: GpDir2d,
}

impl<'a> SurfFunction<'a> {
    pub(crate) fn new(surf: &'a dyn Surface, quad: &'a ImplicitQuad, tol: f64) -> Self {
        Self {
            surf,
            quad,
            tol,
            u: 0.0,
            v: 0.0,
            pntsol: GpPnt::zero(),
            valf: 0.0,
            computed: false,
            tangent: false,
            derived: false,
            tgdu: 0.0,
            tgdv: 0.0,
            gradient: GpVec::zero(),
            d1u: GpVec::zero(),
            d1v: GpVec::zero(),
            d3d: GpVec::zero(),
            d2d: GpDir2d::default(),
        }
    }

    pub(crate) fn tolerance(&self) -> f64 {
        self.tol
    }

    pub(crate) fn root(&self) -> f64 {
        self.valf
    }

    pub(crate) fn point(&self) -> GpPnt {
        self.pntsol
    }

    pub(crate) fn value(&mut self, u: f64, v: f64) -> f64 {
        self.u = u;
        self.v = v;
        self.pntsol = self.surf.d0(u, v);
        self.valf = distance(self.quad, &self.pntsol);
        self.computed = false;
        self.derived = false;
        self.valf
    }

    pub(crate) fn derivatives(&mut self, u: f64, v: f64) -> (f64, f64) {
        self.u = u;
        self.v = v;
        let (p, d1u, d1v) = self.surf.d1(u, v);
        self.pntsol = p;
        self.d1u = d1u;
        self.d1v = d1v;
        self.gradient = gradient(self.quad, &self.pntsol);
        self.computed = false;
        self.derived = true;
        (self.d1u.dot(&self.gradient), self.d1v.dot(&self.gradient))
    }

    pub(crate) fn values(&mut self, u: f64, v: f64) -> (f64, f64, f64) {
        self.u = u;
        self.v = v;
        let (p, d1u, d1v) = self.surf.d1(u, v);
        self.pntsol = p;
        self.d1u = d1u;
        self.d1v = d1v;
        self.valf = distance(self.quad, &self.pntsol);
        self.gradient = gradient(self.quad, &self.pntsol);
        self.computed = false;
        self.derived = true;
        (
            self.valf,
            self.d1u.dot(&self.gradient),
            self.d1v.dot(&self.gradient),
        )
    }

    pub(crate) fn is_tangent(&mut self) -> bool {
        if !self.computed {
            self.computed = true;
            if !self.derived {
                let (p, d1u, d1v) = self.surf.d1(self.u, self.v);
                self.pntsol = p;
                self.d1u = d1u;
                self.d1v = d1v;
                self.derived = true;
            }
            self.tgdu = self.gradient.dot(&self.d1v);
            self.tgdv = -self.gradient.dot(&self.d1u);
            let n2grad = self.gradient.square_magnitude();
            let n2grad_eps = n2grad * EPS_ANG2;
            let n2d1u = self.d1u.square_magnitude();
            let n2d1v = self.d1v.square_magnitude();
            self.tangent = (self.tgdu * self.tgdu <= n2grad_eps * n2d1v)
                && (self.tgdv * self.tgdv <= n2grad_eps * n2d1u);
            if !self.tangent {
                self.d3d = self
                    .d1u
                    .multiplied_scalar(self.tgdu)
                    .added(&self.d1v.multiplied_scalar(self.tgdv));
                match GpDir2d::new(self.tgdu, self.tgdv) {
                    Ok(d) => self.d2d = d,
                    Err(_) => self.tangent = true,
                }
                if self.d3d.magnitude() <= TOL_PETIT {
                    self.tangent = true;
                }
            }
        }
        self.tangent
    }

    pub(crate) fn direction3d(&mut self) -> GpVec {
        let _ = self.is_tangent();
        self.d3d
    }

    pub(crate) fn direction2d(&mut self) -> GpDir2d {
        let _ = self.is_tangent();
        self.d2d
    }
}

/// `math_FunctionSetRoot` for one equation and two variables (minimum-norm Newton).
pub(crate) fn function_set_root(
    func: &mut SurfFunction<'_>,
    uv0: (f64, f64),
    binf: (f64, f64),
    bsup: (f64, f64),
    tol_uv: (f64, f64),
) -> Option<(f64, f64)> {
    let mut u = uv0.0.clamp(binf.0, bsup.0);
    let mut v = uv0.1.clamp(binf.1, bsup.1);
    for _ in 0..100 {
        let (f, ju, jv) = func.values(u, v);
        if f.abs() <= func.tolerance() {
            return Some((u, v));
        }
        let j2 = ju * ju + jv * jv;
        if j2 <= CONFUSION * CONFUSION {
            return if f.abs() <= func.tolerance() {
                Some((u, v))
            } else {
                None
            };
        }
        let du = -f * ju / j2;
        let dv = -f * jv / j2;
        u = (u + du).clamp(binf.0, bsup.0);
        v = (v + dv).clamp(binf.1, bsup.1);
        if du.abs() <= tol_uv.0 && dv.abs() <= tol_uv.1 {
            let _ = func.values(u, v);
            return if func.root().abs() <= func.tolerance() * 100.0 {
                Some((u, v))
            } else {
                None
            };
        }
    }
    let _ = func.values(u, v);
    if func.root().abs() <= func.tolerance() {
        Some((u, v))
    } else {
        None
    }
}

/// Keep a UV box non-empty for Newton.
pub(crate) fn clamp_box(
    mut lo: (f64, f64),
    mut hi: (f64, f64),
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
) -> ((f64, f64), (f64, f64)) {
    lo.0 = lo.0.max(umin);
    lo.1 = lo.1.max(vmin);
    hi.0 = hi.0.min(umax);
    hi.1 = hi.1.min(vmax);
    if hi.0 < lo.0 {
        hi.0 = lo.0;
    }
    if hi.1 < lo.1 {
        hi.1 = lo.1;
    }
    (lo, hi)
}
