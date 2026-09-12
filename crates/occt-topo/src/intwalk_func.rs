//! `IntWalk_TheFunctionOfTheInt2S` / `IntImp_ZerParFunc`.

use occt_core::gp::{GpDir, GpDir2d, GpPnt, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::Surface;

use crate::int_tools_wline::{u_resolution, v_resolution};

use super::iso::{compute_tangence, ConstIso};

/// Three-equation zero of `S1 - S2` with one UV held as iso.
pub(crate) struct ZerParFunc {
    pub ua0: f64,
    pub va0: f64,
    pub ua1: f64,
    pub va1: f64,
    pub ub0: f64,
    pub vb0: f64,
    pub ub1: f64,
    pub vb1: f64,
    pub ures1: f64,
    pub vres1: f64,
    pub ures2: f64,
    pub vres2: f64,
    pub chx_iso: ConstIso,
    pub param_const: f64,
    pub pntsol1: GpPnt,
    pub pntsol2: GpPnt,
    pub dpuv: [GpVec; 4],
    pub f: [f64; 3],
    pub tgduv: [f64; 4],
    pub tangent: bool,
}

impl ZerParFunc {
    pub(crate) fn new(s1: &dyn Surface, s2: &dyn Surface) -> Self {
        let (ua0, ua1) = s1.u_range();
        let (va0, va1) = s1.v_range();
        let (ub0, ub1) = s2.u_range();
        let (vb0, vb1) = s2.v_range();
        Self {
            ua0,
            va0,
            ua1,
            va1,
            ub0,
            vb0,
            ub1,
            vb1,
            ures1: u_resolution(s1, CONFUSION),
            vres1: v_resolution(s1, CONFUSION),
            ures2: u_resolution(s2, CONFUSION),
            vres2: v_resolution(s2, CONFUSION),
            chx_iso: ConstIso::UOnS1,
            param_const: 0.0,
            pntsol1: GpPnt::zero(),
            pntsol2: GpPnt::zero(),
            dpuv: [GpVec::zero(); 4],
            f: [0.0; 3],
            tgduv: [0.0; 4],
            tangent: false,
        }
    }

    pub(crate) fn root(&self) -> f64 {
        self.f[0] * self.f[0] + self.f[1] * self.f[1] + self.f[2] * self.f[2]
    }

    pub(crate) fn point(&self) -> GpPnt {
        GpPnt::new(
            0.5 * (self.pntsol1.x() + self.pntsol2.x()),
            0.5 * (self.pntsol1.y() + self.pntsol2.y()),
            0.5 * (self.pntsol1.z() + self.pntsol2.z()),
        )
    }

    pub(crate) fn direction(&self) -> Option<GpDir> {
        if self.tangent {
            return None;
        }
        let v = self.dpuv[0]
            .multiplied_scalar(self.tgduv[0])
            .added(&self.dpuv[1].multiplied_scalar(self.tgduv[1]));
        GpDir::from_vec(&v).ok()
    }

    pub(crate) fn direction_on_s1(&self) -> GpDir2d {
        GpDir2d::new(self.tgduv[0], self.tgduv[1]).unwrap_or_default()
    }

    pub(crate) fn direction_on_s2(&self) -> GpDir2d {
        GpDir2d::new(self.tgduv[2], self.tgduv[3]).unwrap_or_default()
    }

    pub(crate) fn values(
        &mut self,
        s1: &dyn Surface,
        s2: &dyn Surface,
        x: [f64; 3],
    ) -> ([f64; 3], [[f64; 3]; 3]) {
        let d = self.derivatives(s1, s2, x);
        self.f = [
            self.pntsol1.x() - self.pntsol2.x(),
            self.pntsol1.y() - self.pntsol2.y(),
            self.pntsol1.z() - self.pntsol2.z(),
        ];
        (self.f, d)
    }

    fn derivatives(&mut self, s1: &dyn Surface, s2: &dyn Surface, x: [f64; 3]) -> [[f64; 3]; 3] {
        match self.chx_iso {
            ConstIso::UOnS1 => {
                let (p1, du1, dv1) = s1.d1(self.param_const, x[0]);
                let (p2, du2, dv2) = s2.d1(x[1], x[2]);
                self.pntsol1 = p1;
                self.pntsol2 = p2;
                self.dpuv = [du1, dv1, du2, dv2];
                jac_row(dv1, du2, dv2, true)
            }
            ConstIso::VOnS1 => {
                let (p1, du1, dv1) = s1.d1(x[0], self.param_const);
                let (p2, du2, dv2) = s2.d1(x[1], x[2]);
                self.pntsol1 = p1;
                self.pntsol2 = p2;
                self.dpuv = [du1, dv1, du2, dv2];
                jac_row(du1, du2, dv2, true)
            }
            ConstIso::UOnS2 => {
                let (p1, du1, dv1) = s1.d1(x[0], x[1]);
                let (p2, du2, dv2) = s2.d1(self.param_const, x[2]);
                self.pntsol1 = p1;
                self.pntsol2 = p2;
                self.dpuv = [du1, dv1, du2, dv2];
                [
                    [du1.x(), dv1.x(), -dv2.x()],
                    [du1.y(), dv1.y(), -dv2.y()],
                    [du1.z(), dv1.z(), -dv2.z()],
                ]
            }
            ConstIso::VOnS2 => {
                let (p1, du1, dv1) = s1.d1(x[0], x[1]);
                let (p2, du2, dv2) = s2.d1(x[2], self.param_const);
                self.pntsol1 = p1;
                self.pntsol2 = p2;
                self.dpuv = [du1, dv1, du2, dv2];
                [
                    [du1.x(), dv1.x(), -du2.x()],
                    [du1.y(), dv1.y(), -du2.y()],
                    [du1.z(), dv1.z(), -du2.z()],
                ]
            }
        }
    }

    pub(crate) fn compute_parameters(
        &mut self,
        choix: ConstIso,
        param: [f64; 4],
    ) -> ([f64; 3], [f64; 3], [f64; 3], [f64; 3]) {
        self.chx_iso = choix;
        let (uvap, mut binf, mut bsup, tol) = match choix {
            ConstIso::UOnS1 => {
                self.param_const = param[0];
                (
                    [param[1], param[2], param[3]],
                    [self.va0, self.ub0, self.vb0],
                    [self.va1, self.ub1, self.vb1],
                    [self.vres1, self.ures2, self.vres2],
                )
            }
            ConstIso::VOnS1 => {
                self.param_const = param[1];
                (
                    [param[0], param[2], param[3]],
                    [self.ua0, self.ub0, self.vb0],
                    [self.ua1, self.ub1, self.vb1],
                    [self.ures1, self.ures2, self.vres2],
                )
            }
            ConstIso::UOnS2 => {
                self.param_const = param[2];
                (
                    [param[0], param[1], param[3]],
                    [self.ua0, self.va0, self.vb0],
                    [self.ua1, self.va1, self.vb1],
                    [self.ures1, self.vres1, self.vres2],
                )
            }
            ConstIso::VOnS2 => {
                self.param_const = param[3];
                (
                    [param[0], param[1], param[2]],
                    [self.ua0, self.va0, self.ub0],
                    [self.ua1, self.va1, self.ub1],
                    [self.ures1, self.vres1, self.ures2],
                )
            }
        };
        for i in 0..3 {
            if binf[i].is_finite() && bsup[i].is_finite() {
                let incr = (bsup[i] - binf[i]) * 0.01;
                binf[i] -= incr;
                bsup[i] += incr;
            }
        }
        (uvap, binf, bsup, tol)
    }

    pub(crate) fn is_tangent(&mut self, uvap: [f64; 3], param: &mut [f64; 4]) -> (bool, ConstIso) {
        match self.chx_iso {
            ConstIso::UOnS1 => {
                param[0] = self.param_const;
                param[1] = uvap[0];
                param[2] = uvap[1];
                param[3] = uvap[2];
            }
            ConstIso::VOnS1 => {
                param[1] = self.param_const;
                param[0] = uvap[0];
                param[2] = uvap[1];
                param[3] = uvap[2];
            }
            ConstIso::UOnS2 => {
                param[2] = self.param_const;
                param[0] = uvap[0];
                param[1] = uvap[1];
                param[3] = uvap[2];
            }
            ConstIso::VOnS2 => {
                param[3] = self.param_const;
                param[0] = uvap[0];
                param[1] = uvap[1];
                param[2] = uvap[2];
            }
        }
        let eps = [self.ures1, self.vres1, self.ures2, self.vres2];
        let mut tab = [ConstIso::UOnS1; 4];
        self.tangent = compute_tangence(&self.dpuv, &eps, &mut self.tgduv, &mut tab);
        if !self.tangent {
            self.chx_iso = tab[0];
        }
        (self.tangent, self.chx_iso)
    }
}

fn jac_row(a: GpVec, b: GpVec, c: GpVec, neg_bc: bool) -> [[f64; 3]; 3] {
    let sb = if neg_bc { -1.0 } else { 1.0 };
    [
        [a.x(), sb * b.x(), sb * c.x()],
        [a.y(), sb * b.y(), sb * c.y()],
        [a.z(), sb * b.z(), sb * c.z()],
    ]
}

/// `math_FunctionSetRoot` for 3 equations / 3 variables.
pub(crate) fn function_set_root_3(
    func: &mut ZerParFunc,
    s1: &dyn Surface,
    s2: &dyn Surface,
    uv0: [f64; 3],
    binf: [f64; 3],
    bsup: [f64; 3],
    tol: [f64; 3],
) -> Option<[f64; 3]> {
    let mut x = [
        clamp_finite(uv0[0], binf[0], bsup[0]),
        clamp_finite(uv0[1], binf[1], bsup[1]),
        clamp_finite(uv0[2], binf[2], bsup[2]),
    ];
    for _ in 0..100 {
        let (f, j) = func.values(s1, s2, x);
        let n2 = f[0] * f[0] + f[1] * f[1] + f[2] * f[2];
        if n2.sqrt() <= tol[0].abs().max(tol[1].abs()).max(tol[2].abs()).max(1.0e-14) {
            return Some(x);
        }
        let Some(dx) = solve3(j, [-f[0], -f[1], -f[2]]) else {
            return if n2 <= 1.0e-20 { Some(x) } else { None };
        };
        let mut alpha = 1.0;
        let mut accepted = x;
        let mut best = n2;
        for _ in 0..20 {
            let cand = [
                clamp_finite(x[0] + alpha * dx[0], binf[0], bsup[0]),
                clamp_finite(x[1] + alpha * dx[1], binf[1], bsup[1]),
                clamp_finite(x[2] + alpha * dx[2], binf[2], bsup[2]),
            ];
            let (fc, _) = func.values(s1, s2, cand);
            let nc = fc[0] * fc[0] + fc[1] * fc[1] + fc[2] * fc[2];
            if nc < best {
                best = nc;
                accepted = cand;
                break;
            }
            alpha *= 0.5;
        }
        let du = [
            (accepted[0] - x[0]).abs(),
            (accepted[1] - x[1]).abs(),
            (accepted[2] - x[2]).abs(),
        ];
        x = accepted;
        if du[0] <= tol[0] && du[1] <= tol[1] && du[2] <= tol[2] {
            let _ = func.values(s1, s2, x);
            return Some(x);
        }
    }
    let _ = func.values(s1, s2, x);
    Some(x)
}

fn clamp_finite(x: f64, lo: f64, hi: f64) -> f64 {
    let mut y = x;
    if lo.is_finite() {
        y = y.max(lo);
    }
    if hi.is_finite() {
        y = y.min(hi);
    }
    y
}

fn solve3(a: [[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let mut m = [
        [a[0][0], a[0][1], a[0][2], b[0]],
        [a[1][0], a[1][1], a[1][2], b[1]],
        [a[2][0], a[2][1], a[2][2], b[2]],
    ];
    for col in 0..3 {
        let mut piv = col;
        let mut best = m[col][col].abs();
        for row in (col + 1)..3 {
            let v = m[row][col].abs();
            if v > best {
                best = v;
                piv = row;
            }
        }
        if best < 1.0e-30 {
            return None;
        }
        if piv != col {
            m.swap(col, piv);
        }
        let diag = m[col][col];
        for j in col..4 {
            m[col][j] /= diag;
        }
        for row in 0..3 {
            if row == col {
                continue;
            }
            let f = m[row][col];
            for j in col..4 {
                m[row][j] -= f * m[col][j];
            }
        }
    }
    Some([m[0][3], m[1][3], m[2][3]])
}
