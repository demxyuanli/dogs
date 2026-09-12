//! `IntWalk_TheInt2S` (`IntImp_Int2S.gxx`).

use occt_core::gp::{GpDir, GpDir2d};
use occt_core::precision::CONFUSION;
use occt_geom::Surface;

use crate::int_tools_wline::{u_resolution, v_resolution, PntOn2S};

use super::func::{function_set_root_3, ZerParFunc};
use super::iso::{compute_tangence, ConstIso};

/// Intersection of two parametric surfaces at a close UV guess.
pub struct TheInt2S<'a> {
    s1: &'a dyn Surface,
    s2: &'a dyn Surface,
    func: ZerParFunc,
    done: bool,
    empty: bool,
    pint: PntOn2S,
    tangent: bool,
    d3d: GpDir,
    d2d1: GpDir2d,
    d2d2: GpDir2d,
    tol: f64,
    ua0: f64,
    va0: f64,
    ua1: f64,
    va1: f64,
    ub0: f64,
    vb0: f64,
    ub1: f64,
    vb1: f64,
}

impl<'a> TheInt2S<'a> {
    pub fn new(s1: &'a dyn Surface, s2: &'a dyn Surface, tol_tangency: f64) -> Self {
        let func = ZerParFunc::new(s1, s2);
        let (ua0, ua1) = s1.u_range();
        let (va0, va1) = s1.v_range();
        let (ub0, ub1) = s2.u_range();
        let (vb0, vb1) = s2.v_range();
        Self {
            s1,
            s2,
            func,
            done: true,
            empty: true,
            pint: PntOn2S {
                p: occt_core::gp::GpPnt::zero(),
                u1: 0.0,
                v1: 0.0,
                u2: 0.0,
                v2: 0.0,
            },
            tangent: false,
            d3d: GpDir::default_dir(),
            d2d1: GpDir2d::default(),
            d2d2: GpDir2d::default(),
            tol: tol_tangency * tol_tangency,
            ua0,
            va0,
            ua1,
            va1,
            ub0,
            vb0,
            ub1,
            vb1,
        }
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn is_empty(&self) -> bool {
        self.empty
    }

    pub fn is_tangent(&self) -> bool {
        self.tangent
    }

    pub fn point(&self) -> PntOn2S {
        self.pint
    }

    pub fn change_point(&mut self) -> &mut PntOn2S {
        &mut self.pint
    }

    pub fn direction(&self) -> GpDir {
        self.d3d
    }

    pub fn direction_on_s1(&self) -> GpDir2d {
        self.d2d1
    }

    pub fn direction_on_s2(&self) -> GpDir2d {
        self.d2d2
    }

    /// `Perform(Param, Rsnld, ChoixIso)`.
    pub fn perform_iso(&mut self, param: [f64; 4], choix: ConstIso) -> ConstIso {
        self.done = true;
        let (uvap, binf, bsup, tolerance) = self.func.compute_parameters(choix, param);
        let Some(root) = function_set_root_3(
            &mut self.func,
            self.s1,
            self.s2,
            uvap,
            binf,
            bsup,
            tolerance,
        ) else {
            self.empty = true;
            return choix;
        };
        let mut best = choix;
        if self.func.root().abs() <= self.tol {
            let mut uvres = [0.0; 4];
            let (tang, bc) = self.func.is_tangent(root, &mut uvres);
            self.empty = false;
            self.tangent = tang;
            best = bc;
            let p = self.func.point();
            self.pint = PntOn2S {
                p,
                u1: uvres[0],
                v1: uvres[1],
                u2: uvres[2],
                v2: uvres[3],
            };
            if !self.tangent {
                if let Some(d) = self.func.direction() {
                    self.d3d = d;
                }
                self.d2d1 = self.func.direction_on_s1();
                self.d2d2 = self.func.direction_on_s2();
            }
        } else {
            self.empty = true;
        }
        best
    }

    /// `Perform(Param, Rsnld)` — pick the best iso via `ComputeTangence`.
    pub fn perform(&mut self, param: [f64; 4]) -> ConstIso {
        let (p1, du1, dv1) = self.s1.d1(param[0], param[1]);
        let (p2, du2, dv2) = self.s2.d1(param[2], param[3]);
        let _ = (p1, p2);
        let dpuv = [du1, dv1, du2, dv2];
        let eps = [
            u_resolution(self.s1, CONFUSION),
            v_resolution(self.s1, CONFUSION),
            u_resolution(self.s2, CONFUSION),
            v_resolution(self.s2, CONFUSION),
        ];
        let mut tg = [0.0; 4];
        let mut choix = [ConstIso::UOnS1; 4];
        self.empty = true;
        if compute_tangence(&dpuv, &eps, &mut tg, &mut choix) {
            return ConstIso::UOnS1;
        }
        let mut best = choix[0];
        let mut i = 0;
        while self.empty && i <= 3 {
            let current = self.perform_iso(param, choix[i]);
            if !self.empty {
                best = current;
            }
            i += 1;
        }
        if !self.empty {
            self.snap_to_bounds(eps);
        }
        best
    }

    fn snap_to_bounds(&mut self, eps: [f64; 4]) {
        let mut duv = [self.pint.u1, self.pint.v1, self.pint.u2, self.pint.v2];
        let uvd = [self.ua0, self.va0, self.ub0, self.vb0];
        let uvf = [self.ua1, self.va1, self.ub1, self.vb1];
        for i in 0..4 {
            if uvd[i].is_finite() && duv[i] <= uvd[i] - eps[i] {
                duv[i] = uvd[i];
            } else if uvf[i].is_finite() && duv[i] >= uvf[i] + eps[i] {
                duv[i] = uvf[i];
            }
        }
        self.pint.u1 = duv[0];
        self.pint.v1 = duv[1];
        self.pint.u2 = duv[2];
        self.pint.v2 = duv[3];
    }
}
