//! `Plate_Plate::SolEm` — the `(iu, iv)`-th derivative of the fundamental
//! solution of the Laplacian at the power order.
//!
//! Source: `Plate_Plate.cxx:1373-1883` (`SolEm`), `Plate_Plate.lxx:20-47`
//! (`Polm`).
//!
//! The `(u <-> v)` swap of `cxx:1384-1398`, the `Uold` / `Vold` memoisation of
//! `cxx:1401-1420` (OCCT reaches it through a `const_cast`; the port uses
//! `Cell`s) and the hand-unrolled `R^expo` of `cxx:1423-1448` are reproduced
//! verbatim. Every `(IU, IV)` pair that OCCT leaves with `DUV = 0` (its
//! `default:` arms) falls into the final `_ => 0.0` arm.

use occt_core::gp::GpXY;

use super::Plate;

impl Plate {
    /// `SolEm(point2d, iu, iv)` (`cxx:1375-1883`).
    pub(super) fn sol_em(&self, point2d: GpXY, iu: i32, iv: i32) -> f64 {
        let (i_u, i_v, u, v) = if iv > iu {
            // SolEm is symmetric in (u <-> v): swap u and v if iv > iu
            // to avoid some code
            (iv, iu, point2d.y() * self.ddv[1], point2d.x() * self.ddu[1])
        } else {
            (iu, iv, point2d.x() * self.ddu[1], point2d.y() * self.ddv[1])
        };

        if u == self.uold.get() && v == self.vold.get() {
            if self.r_cache.get() < 1.0e-20 {
                return 0.0;
            }
        } else {
            self.uold.set(u);
            self.vold.set(v);
            self.u2.set(u * u);
            let r = self.u2.get() + v * v;
            self.r_cache.set(r);
            if r < 1.0e-20 {
                return 0.0;
            }
            self.l.set(r.ln());
        }

        let big_r = self.r_cache.get();
        let big_l = self.l.get();
        let big_u2 = self.u2.get();
        let m = self.order as f64;
        let mm1 = (self.order - 1) as f64;

        // `int expo = mm1 - IU - IV;` with the unrolled power of `R`.
        let expo = (self.order - 1) - i_u - i_v;
        let pr = if expo == 0 {
            1.0
        } else if expo < 0 {
            let mut pr = big_r;
            for _ in 1..(-expo) {
                pr *= big_r;
            }
            1.0 / pr
        } else {
            let mut pr = big_r;
            for _ in 1..expo {
                pr *= big_r;
            }
            pr
        };

        let duv = match (i_u, i_v) {
            (0, 0) => pr * big_l,

            (1, 0) => 2.0 * pr * u * (1.0 + big_l * mm1),

            (1, 1) => {
                let m2 = m * m;
                4.0 * pr * u * v * ((2.0 * m - 3.0) + (m2 - 3.0 * m + 2.0) * big_l)
            }

            (2, 0) => {
                let m2 = m * m;
                2.0 * pr
                    * (big_r - big_l * big_r + big_l * m * big_r - 6.0 * big_u2
                        + 4.0 * big_l * big_u2
                        + 4.0 * m * big_u2
                        - 6.0 * big_l * m * big_u2
                        + 2.0 * big_l * m2 * big_u2)
            }

            (2, 1) => {
                let m2 = m * m;
                let m3 = m2 * m;
                let mut duv = -3.0 * big_r
                    + 2.0 * big_l * big_r
                    + 2.0 * m * big_r
                    - 3.0 * big_l * m * big_r
                    + big_l * m2 * big_r
                    + 22.0 * big_u2
                    - 12.0 * big_l * big_u2
                    - 24.0 * m * big_u2
                    + 22.0 * big_l * m * big_u2
                    + 6.0 * m2 * big_u2
                    - 12.0 * big_l * m2 * big_u2
                    + 2.0 * big_l * m3 * big_u2;
                duv *= 4.0 * pr * v;
                duv
            }

            (2, 2) => {
                let m2 = m * m;
                let m3 = m2 * m;
                let m4 = m2 * m2;
                let v2 = v * v;
                let r2 = big_r * big_r;
                let mut duv = -3.0 * r2
                    + 2.0 * big_l * r2
                    + 2.0 * m * r2
                    - 3.0 * big_l * m * r2
                    + big_l * m2 * r2
                    + 22.0 * big_r * big_u2
                    - 12.0 * big_l * big_r * big_u2
                    - 24.0 * m * big_r * big_u2
                    + 22.0 * big_l * m * big_r * big_u2
                    + 6.0 * m2 * big_r * big_u2
                    - 12.0 * big_l * m2 * big_r * big_u2;
                duv += 2.0 * big_l * m3 * big_r * big_u2 + 22.0 * big_r * v2
                    - 12.0 * big_l * big_r * v2
                    - 24.0 * m * big_r * v2
                    + 22.0 * big_l * m * big_r * v2
                    + 6.0 * m2 * big_r * v2
                    - 12.0 * big_l * m2 * big_r * v2
                    + 2.0 * big_l * m3 * big_r * v2
                    - 200.0 * big_u2 * v2
                    + 96.0 * big_l * big_u2 * v2;
                duv += 280.0 * m * big_u2 * v2
                    - 200.0 * big_l * m * big_u2 * v2
                    - 120.0 * m2 * big_u2 * v2
                    + 140.0 * big_l * m2 * big_u2 * v2
                    + 16.0 * m3 * big_u2 * v2
                    - 40.0 * big_l * m3 * big_u2 * v2
                    + 4.0 * big_l * m4 * big_u2 * v2;
                duv * 4.0 * pr
            }

            (3, 0) => {
                let m2 = m * m;
                let m3 = m2 * m;
                let mut duv = -9.0 * big_r
                    + 6.0 * big_l * big_r
                    + 6.0 * m * big_r
                    - 9.0 * big_l * m * big_r
                    + 3.0 * big_l * m2 * big_r
                    + 22.0 * big_u2
                    - 12.0 * big_l * big_u2
                    - 24.0 * m * big_u2
                    + 22.0 * big_l * m * big_u2
                    + 6.0 * m2 * big_u2
                    - 12.0 * big_l * m2 * big_u2
                    + 2.0 * big_l * m3 * big_u2;
                duv *= 4.0 * pr * u;
                duv
            }

            (3, 1) => {
                let m2 = m * m;
                let m3 = m2 * m;
                let m4 = m2 * m2;
                let mut duv = 33.0 * big_r
                    - 18.0 * big_l * big_r
                    - 36.0 * m * big_r
                    + 33.0 * big_l * m * big_r
                    + 9.0 * m2 * big_r
                    - 18.0 * big_l * m2 * big_r
                    + 3.0 * big_l * m3 * big_r
                    - 100.0 * big_u2
                    + 48.0 * big_l * big_u2
                    + 140.0 * m * big_u2
                    - 100.0 * big_l * m * big_u2
                    - 60.0 * m2 * big_u2
                    + 70.0 * big_l * m2 * big_u2;
                duv += 8.0 * m3 * big_u2 - 20.0 * big_l * m3 * big_u2
                    + 2.0 * big_l * m4 * big_u2;
                duv *= 8.0 * pr * u * v;
                duv
            }

            (3, 2) => {
                let m2 = m * m;
                let m3 = m2 * m;
                let m4 = m2 * m2;
                let m5 = m4 * m;
                let ru2 = big_r * big_u2;
                let v2 = v * v;
                let rv2 = big_r * v2;
                let u2v2 = v2 * big_u2;
                let r2 = big_r * big_r;

                let mut duv = -100.0 * ru2
                    + 48.0 * big_l * ru2
                    + 140.0 * m * ru2
                    - 100.0 * big_l * m * ru2
                    - 60.0 * m2 * ru2
                    + 70.0 * big_l * m2 * ru2
                    + 8.0 * m3 * ru2
                    - 20.0 * big_l * m3 * ru2
                    + 2.0 * big_l * m4 * ru2
                    - 300.0 * rv2
                    + 144.0 * big_l * rv2
                    + 420.0 * m * rv2
                    - 300.0 * big_l * m * rv2
                    - 180.0 * m2 * rv2
                    + 210.0 * big_l * m2 * rv2
                    + 24.0 * m3 * rv2
                    - 60.0 * big_l * m3 * rv2
                    + 6.0 * big_l * m4 * rv2
                    + 33.0 * r2
                    - 18.0 * big_l * r2
                    - 36.0 * m * r2
                    + 33.0 * big_l * m * r2
                    + 9.0 * m2 * r2
                    - 18.0 * big_l * m2 * r2
                    + 3.0 * big_l * m3 * r2
                    + 1096.0 * u2v2
                    - 480.0 * big_l * u2v2
                    - 1800.0 * m * u2v2
                    + 1096.0 * big_l * m * u2v2
                    + 1020.0 * m2 * u2v2
                    - 900.0 * big_l * m2 * u2v2
                    - 240.0 * m3 * u2v2
                    + 340.0 * big_l * m3 * u2v2
                    + 20.0 * m4 * u2v2
                    - 60.0 * big_l * m4 * u2v2
                    + 4.0 * big_l * m5 * u2v2;
                duv *= 8.0 * pr * u;
                duv
            }

            (3, 3) => {
                let m2 = m * m;
                let m3 = m2 * m;
                let m4 = m2 * m2;
                let m5 = m3 * m2;
                let m6 = m3 * m3;
                let ru2 = big_r * big_u2;
                let v2 = v * v;
                let rv2 = big_r * v2;
                let u2v2 = v2 * big_u2;
                let r2 = big_r * big_r;

                let mut duv = 1644.0 * ru2
                    - 720.0 * big_l * ru2
                    - 2700.0 * m * ru2
                    + 1644.0 * big_l * m * ru2
                    + 1530.0 * m2 * ru2
                    - 1350.0 * big_l * m2 * ru2
                    - 360.0 * m3 * ru2
                    + 510.0 * big_l * m3 * ru2
                    + 30.0 * m4 * ru2
                    - 90.0 * big_l * m4 * ru2
                    + 6.0 * big_l * m5 * ru2
                    + 1644.0 * rv2
                    - 720.0 * big_l * rv2
                    - 2700.0 * m * rv2
                    + 1644.0 * big_l * m * rv2
                    + 1530.0 * m2 * rv2
                    - 1350.0 * big_l * m2 * rv2
                    - 360.0 * m3 * rv2
                    + 510.0 * big_l * m3 * rv2
                    + 30.0 * m4 * rv2
                    - 90.0 * big_l * m4 * rv2
                    + 6.0 * big_l * m5 * rv2
                    - 450.0 * r2
                    + 216.0 * big_l * r2
                    + 630.0 * m * r2
                    - 450.0 * big_l * m * r2
                    - 270.0 * m2 * r2
                    + 315.0 * big_l * m2 * r2
                    + 36.0 * m3 * r2
                    - 90.0 * big_l * m3 * r2
                    + 9.0 * big_l * m4 * r2
                    - 7056.0 * u2v2
                    + 2880.0 * big_l * u2v2
                    + 12992.0 * m * u2v2
                    - 7056.0 * big_l * m * u2v2
                    - 8820.0 * m2 * u2v2
                    + 6496.0 * big_l * m2 * u2v2
                    + 2800.0 * m3 * u2v2
                    - 2940.0 * big_l * m3 * u2v2
                    - 420.0 * m4 * u2v2
                    + 700.0 * big_l * m4 * u2v2
                    + 24.0 * m5 * u2v2
                    - 84.0 * big_l * m5 * u2v2
                    + 4.0 * big_l * m6 * u2v2;
                duv *= 16.0 * pr * u * v;
                duv
            }

            (4, 0) => {
                let m2 = m * m;
                let m3 = m2 * m;
                let m4 = m2 * m2;
                let u4 = big_u2 * big_u2;
                let r2 = big_r * big_r;
                let mut duv = -9.0 * r2
                    + 6.0 * big_l * r2
                    + 6.0 * m * r2
                    - 9.0 * big_l * m * r2
                    + 3.0 * big_l * m2 * r2
                    + 132.0 * big_r * big_u2
                    - 72.0 * big_l * big_r * big_u2
                    - 144.0 * m * big_r * big_u2
                    + 132.0 * big_l * m * big_r * big_u2
                    + 36.0 * m2 * big_r * big_u2
                    - 72.0 * big_l * m2 * big_r * big_u2;
                duv += 12.0 * big_l * m3 * big_r * big_u2 - 200.0 * u4
                    + 96.0 * big_l * u4
                    + 280.0 * m * u4
                    - 200.0 * big_l * m * u4
                    - 120.0 * m2 * u4
                    + 140.0 * big_l * m2 * u4
                    + 16.0 * m3 * u4
                    - 40.0 * big_l * m3 * u4
                    + 4.0 * big_l * m4 * u4;
                duv * 4.0 * pr
            }

            (4, 1) => {
                let m2 = m * m;
                let m3 = m2 * m;
                let m4 = m2 * m2;
                let m5 = m2 * m3;
                let u4 = big_u2 * big_u2;
                let ru2 = big_r * big_u2;
                let r2 = big_r * big_r;

                let mut duv = -600.0 * ru2
                    + 288.0 * big_l * ru2
                    + 840.0 * m * ru2
                    - 600.0 * big_l * m * ru2
                    - 360.0 * m2 * ru2
                    + 420.0 * big_l * m2 * ru2
                    + 48.0 * m3 * ru2
                    - 120.0 * big_l * m3 * ru2
                    + 12.0 * big_l * m4 * ru2
                    + 33.0 * r2
                    - 18.0 * big_l * r2
                    - 36.0 * m * r2
                    + 33.0 * big_l * m * r2
                    + 9.0 * m2 * r2
                    - 18.0 * big_l * m2 * r2
                    + 3.0 * big_l * m3 * r2
                    + 1096.0 * u4
                    - 480.0 * big_l * u4
                    - 1800.0 * m * u4
                    + 1096.0 * big_l * m * u4
                    + 1020.0 * m2 * u4
                    - 900.0 * big_l * m2 * u4
                    - 240.0 * m3 * u4
                    + 340.0 * big_l * m3 * u4
                    + 20.0 * m4 * u4
                    - 60.0 * big_l * m4 * u4
                    + 4.0 * big_l * m5 * u4;
                duv *= 8.0 * pr * v;
                duv
            }

            (4, 2) => {
                let m2 = m * m;
                let m3 = m2 * m;
                let m4 = m2 * m2;
                let m5 = m2 * m3;
                let m6 = m3 * m3;
                let u4 = big_u2 * big_u2;
                let r2 = big_r * big_r;
                let r3 = r2 * big_r;
                let v2 = v * v;
                let u2v2 = v2 * big_u2;
                let ru2v2 = big_r * u2v2;
                let u4v2 = u4 * v2;
                let r2u2 = r2 * big_u2;
                let ru4 = big_r * u4;
                let r2v2 = r2 * v2;

                let mut duv = 6576.0 * ru2v2
                    - 2880.0 * big_l * ru2v2
                    - 10800.0 * m * ru2v2
                    + 6576.0 * big_l * m * ru2v2
                    + 6120.0 * m2 * ru2v2
                    - 5400.0 * big_l * m2 * ru2v2
                    - 1440.0 * m3 * ru2v2
                    + 2040.0 * big_l * m3 * ru2v2
                    + 120.0 * m4 * ru2v2
                    - 360.0 * big_l * m4 * ru2v2
                    + 24.0 * big_l * m5 * ru2v2
                    + 1096.0 * ru4
                    - 480.0 * big_l * ru4
                    - 1800.0 * m * ru4
                    + 1096.0 * big_l * m * ru4
                    + 1020.0 * m2 * ru4
                    - 900.0 * big_l * m2 * ru4
                    - 240.0 * m3 * ru4
                    + 340.0 * big_l * m3 * ru4
                    + 20.0 * m4 * ru4
                    - 60.0 * big_l * m4 * ru4
                    + 4.0 * big_l * m5 * ru4
                    - 600.0 * r2u2
                    + 288.0 * big_l * r2u2
                    + 840.0 * m * r2u2
                    - 600.0 * big_l * m * r2u2
                    - 360.0 * m2 * r2u2
                    + 420.0 * big_l * m2 * r2u2
                    + 48.0 * m3 * r2u2
                    - 120.0 * big_l * m3 * r2u2
                    + 12.0 * big_l * m4 * r2u2
                    - 300.0 * r2v2
                    + 144.0 * big_l * r2v2
                    + 420.0 * m * r2v2
                    - 300.0 * big_l * m * r2v2
                    - 180.0 * m2 * r2v2
                    + 210.0 * big_l * m2 * r2v2
                    + 24.0 * m3 * r2v2
                    - 60.0 * big_l * m3 * r2v2
                    + 6.0 * big_l * m4 * r2v2
                    + 33.0 * r3
                    - 18.0 * big_l * r3
                    - 36.0 * m * r3
                    + 33.0 * big_l * m * r3
                    + 9.0 * m2 * r3
                    - 18.0 * big_l * m2 * r3
                    + 3.0 * big_l * m3 * r3
                    - 14112.0 * u4v2
                    + 5760.0 * big_l * u4v2
                    + 25984.0 * m * u4v2
                    - 14112.0 * big_l * m * u4v2
                    - 17640.0 * m2 * u4v2
                    + 12992.0 * big_l * m2 * u4v2
                    + 5600.0 * m3 * u4v2
                    - 5880.0 * big_l * m3 * u4v2
                    - 840.0 * m4 * u4v2
                    + 1400.0 * big_l * m4 * u4v2
                    + 48.0 * m5 * u4v2
                    - 168.0 * big_l * m5 * u4v2
                    + 8.0 * big_l * m6 * u4v2;
                duv *= 8.0 * pr;
                duv
            }

            (4, 3) => {
                let m2 = m * m;
                let m3 = m2 * m;
                let m4 = m2 * m2;
                let m5 = m2 * m3;
                let m6 = m3 * m3;
                let m7 = m3 * m4;
                let u4 = big_u2 * big_u2;
                let r2 = big_r * big_r;
                let r3 = r2 * big_r;
                let v2 = v * v;
                let u2v2 = v2 * big_u2;
                let ru2v2 = big_r * u2v2;
                let u4v2 = u4 * v2;
                let r2u2 = r2 * big_u2;
                let r2v2 = r2 * v2;
                let ru4 = big_r * u4;

                let mut duv = -42336.0 * ru2v2
                    + 17280.0 * big_l * ru2v2
                    + 77952.0 * m * ru2v2
                    - 42336.0 * big_l * m * ru2v2
                    - 52920.0 * m2 * ru2v2
                    + 38976.0 * big_l * m2 * ru2v2
                    + 16800.0 * m3 * ru2v2
                    - 17640.0 * big_l * m3 * ru2v2
                    - 2520.0 * m4 * ru2v2
                    + 4200.0 * big_l * m4 * ru2v2
                    + 144.0 * m5 * ru2v2
                    - 504.0 * big_l * m5 * ru2v2
                    + 24.0 * big_l * m6 * ru2v2
                    - 21168.0 * ru4
                    + 8640.0 * big_l * ru4
                    + 38976.0 * m * ru4
                    - 21168.0 * big_l * m * ru4
                    - 26460.0 * m2 * ru4
                    + 19488.0 * big_l * m2 * ru4
                    + 8400.0 * m3 * ru4
                    - 8820.0 * big_l * m3 * ru4
                    - 1260.0 * m4 * ru4
                    + 2100.0 * big_l * m4 * ru4
                    + 72.0 * m5 * ru4
                    - 252.0 * big_l * m5 * ru4
                    + 12.0 * big_l * m6 * ru4
                    + 9864.0 * r2u2
                    - 4320.0 * big_l * r2u2
                    - 16200.0 * m * r2u2
                    + 9864.0 * big_l * m * r2u2
                    + 9180.0 * m2 * r2u2
                    - 8100.0 * big_l * m2 * r2u2
                    - 2160.0 * m3 * r2u2
                    + 3060.0 * big_l * m3 * r2u2
                    + 180.0 * m4 * r2u2
                    - 540.0 * big_l * m4 * r2u2
                    + 36.0 * big_l * m5 * r2u2
                    + 1644.0 * r2v2
                    - 720.0 * big_l * r2v2
                    - 2700.0 * m * r2v2
                    + 1644.0 * big_l * m * r2v2
                    + 1530.0 * m2 * r2v2
                    - 1350.0 * big_l * m2 * r2v2
                    - 360.0 * m3 * r2v2
                    + 510.0 * big_l * m3 * r2v2
                    + 30.0 * m4 * r2v2
                    - 90.0 * big_l * m4 * r2v2
                    + 6.0 * big_l * m5 * r2v2
                    - 450.0 * r3
                    + 216.0 * big_l * r3
                    + 630.0 * m * r3
                    - 450.0 * big_l * m * r3
                    - 270.0 * m2 * r3
                    + 315.0 * big_l * m2 * r3
                    + 36.0 * m3 * r3
                    - 90.0 * big_l * m3 * r3
                    + 9.0 * big_l * m4 * r3
                    + 104544.0 * u4v2
                    - 40320.0 * big_l * u4v2
                    - 210112.0 * m * u4v2
                    + 104544.0 * big_l * m * u4v2
                    + 162456.0 * m2 * u4v2
                    - 105056.0 * big_l * m2 * u4v2
                    - 62720.0 * m3 * u4v2
                    + 54152.0 * big_l * m3 * u4v2
                    + 12880.0 * m4 * u4v2
                    - 15680.0 * big_l * m4 * u4v2
                    - 1344.0 * m5 * u4v2
                    + 2576.0 * big_l * m5 * u4v2
                    + 56.0 * m6 * u4v2
                    - 224.0 * big_l * m6 * u4v2
                    + 8.0 * big_l * m7 * u4v2;
                duv *= 16.0 * pr * v;
                duv
            }

            (5, 0) => {
                let m2 = m * m;
                let m3 = m2 * m;
                let m4 = m2 * m2;
                let m5 = m2 * m3;
                let u4 = big_u2 * big_u2;
                let r2 = big_r * big_r;
                let ru2 = big_r * big_u2;

                let mut duv = -1000.0 * ru2
                    + 480.0 * big_l * ru2
                    + 1400.0 * m * ru2
                    - 1000.0 * big_l * m * ru2
                    - 600.0 * m2 * ru2
                    + 700.0 * big_l * m2 * ru2
                    + 80.0 * m3 * ru2
                    - 200.0 * big_l * m3 * ru2
                    + 20.0 * big_l * m4 * ru2
                    + 165.0 * r2
                    - 90.0 * big_l * r2
                    - 180.0 * m * r2
                    + 165.0 * big_l * m * r2
                    + 45.0 * m2 * r2
                    - 90.0 * big_l * m2 * r2
                    + 15.0 * big_l * m3 * r2
                    + 1096.0 * u4
                    - 480.0 * big_l * u4
                    - 1800.0 * m * u4
                    + 1096.0 * big_l * m * u4
                    + 1020.0 * m2 * u4
                    - 900.0 * big_l * m2 * u4
                    - 240.0 * m3 * u4
                    + 340.0 * big_l * m3 * u4
                    + 20.0 * m4 * u4
                    - 60.0 * big_l * m4 * u4
                    + 4.0 * big_l * m5 * u4;
                duv *= 8.0 * pr * u;
                duv
            }

            (5, 1) => {
                let m2 = m * m;
                let m3 = m2 * m;
                let m4 = m2 * m2;
                let m5 = m2 * m3;
                let m6 = m3 * m3;
                let u4 = big_u2 * big_u2;
                let ru2 = big_r * big_u2;
                let r2 = big_r * big_r;

                let mut duv = 5480.0 * ru2
                    - 2400.0 * big_l * ru2
                    - 9000.0 * m * ru2
                    + 5480.0 * big_l * m * ru2
                    + 5100.0 * m2 * ru2
                    - 4500.0 * big_l * m2 * ru2
                    - 1200.0 * m3 * ru2
                    + 1700.0 * big_l * m3 * ru2
                    + 100.0 * m4 * ru2
                    - 300.0 * big_l * m4 * ru2
                    + 20.0 * big_l * m5 * ru2
                    - 750.0 * r2
                    + 360.0 * big_l * r2
                    + 1050.0 * m * r2
                    - 750.0 * big_l * m * r2
                    - 450.0 * m2 * r2
                    + 525.0 * big_l * m2 * r2
                    + 60.0 * m3 * r2
                    - 150.0 * big_l * m3 * r2
                    + 15.0 * big_l * m4 * r2
                    - 7056.0 * u4
                    + 2880.0 * big_l * u4
                    + 12992.0 * m * u4
                    - 7056.0 * big_l * m * u4
                    - 8820.0 * m2 * u4
                    + 6496.0 * big_l * m2 * u4
                    + 2800.0 * m3 * u4
                    - 2940.0 * big_l * m3 * u4
                    - 420.0 * m4 * u4
                    + 700.0 * big_l * m4 * u4
                    + 24.0 * m5 * u4
                    - 84.0 * big_l * m5 * u4
                    + 4.0 * big_l * m6 * u4;
                duv *= 16.0 * pr * u * v;
                duv
            }

            (5, 2) => {
                let m2 = m * m;
                let m3 = m2 * m;
                let m4 = m2 * m2;
                let m5 = m2 * m3;
                let m6 = m3 * m3;
                let m7 = m3 * m4;
                let u4 = big_u2 * big_u2;
                let r2 = big_r * big_r;
                let r3 = r2 * big_r;
                let v2 = v * v;
                let u2v2 = v2 * big_u2;
                let ru2v2 = big_r * u2v2;
                let u4v2 = u4 * v2;
                let r2u2 = r2 * big_u2;
                let r2v2 = r2 * v2;
                let ru4 = big_r * u4;

                let mut duv = -70560.0 * ru2v2
                    + 28800.0 * big_l * ru2v2
                    + 129920.0 * m * ru2v2
                    - 70560.0 * big_l * m * ru2v2
                    - 88200.0 * m2 * ru2v2
                    + 64960.0 * big_l * m2 * ru2v2
                    + 28000.0 * m3 * ru2v2
                    - 29400.0 * big_l * m3 * ru2v2
                    - 4200.0 * m4 * ru2v2
                    + 7000.0 * big_l * m4 * ru2v2
                    + 240.0 * m5 * ru2v2
                    - 840.0 * big_l * m5 * ru2v2
                    + 40.0 * big_l * m6 * ru2v2
                    - 7056.0 * ru4
                    + 2880.0 * big_l * ru4
                    + 12992.0 * m * ru4
                    - 7056.0 * big_l * m * ru4
                    - 8820.0 * m2 * ru4
                    + 6496.0 * big_l * m2 * ru4
                    + 2800.0 * m3 * ru4
                    - 2940.0 * big_l * m3 * ru4
                    - 420.0 * m4 * ru4
                    + 700.0 * big_l * m4 * ru4
                    + 24.0 * m5 * ru4
                    - 84.0 * big_l * m5 * ru4
                    + 4.0 * big_l * m6 * ru4
                    + 5480.0 * r2u2
                    - 2400.0 * big_l * r2u2
                    - 9000.0 * m * r2u2
                    + 5480.0 * big_l * m * r2u2
                    + 5100.0 * m2 * r2u2
                    - 4500.0 * big_l * m2 * r2u2
                    - 1200.0 * m3 * r2u2
                    + 1700.0 * big_l * m3 * r2u2
                    + 100.0 * m4 * r2u2
                    - 300.0 * big_l * m4 * r2u2
                    + 20.0 * big_l * m5 * r2u2
                    + 8220.0 * r2v2
                    - 3600.0 * big_l * r2v2
                    - 13500.0 * m * r2v2
                    + 8220.0 * big_l * m * r2v2
                    + 7650.0 * m2 * r2v2
                    - 6750.0 * big_l * m2 * r2v2
                    - 1800.0 * m3 * r2v2
                    + 2550.0 * big_l * m3 * r2v2
                    + 150.0 * m4 * r2v2
                    - 450.0 * big_l * m4 * r2v2
                    + 30.0 * big_l * m5 * r2v2
                    - 750.0 * r3
                    + 360.0 * big_l * r3
                    + 1050.0 * m * r3
                    - 750.0 * big_l * m * r3
                    - 450.0 * m2 * r3
                    + 525.0 * big_l * m2 * r3
                    + 60.0 * m3 * r3
                    - 150.0 * big_l * m3 * r3
                    + 15.0 * big_l * m4 * r3
                    + 104544.0 * u4v2
                    - 40320.0 * big_l * u4v2
                    - 210112.0 * m * u4v2
                    + 104544.0 * big_l * m * u4v2
                    + 162456.0 * m2 * u4v2
                    - 105056.0 * big_l * m2 * u4v2
                    - 62720.0 * m3 * u4v2
                    + 54152.0 * big_l * m3 * u4v2
                    + 12880.0 * m4 * u4v2
                    - 15680.0 * big_l * m4 * u4v2
                    - 1344.0 * m5 * u4v2
                    + 2576.0 * big_l * m5 * u4v2
                    + 56.0 * m6 * u4v2
                    - 224.0 * big_l * m6 * u4v2
                    + 8.0 * big_l * m7 * u4v2;
                duv *= 16.0 * pr * u;
                duv
            }

            (6, 0) => {
                let m2 = m * m;
                let m3 = m2 * m;
                let m4 = m2 * m2;
                let m5 = m2 * m3;
                let m6 = m3 * m3;
                let u4 = big_u2 * big_u2;
                let u6 = big_u2 * u4;
                let r2 = big_r * big_r;
                let r3 = r2 * big_r;
                let r2u2 = r2 * big_u2;
                let ru4 = big_r * u4;

                let mut duv = 16440.0 * ru4
                    - 7200.0 * big_l * ru4
                    - 27000.0 * m * ru4
                    + 16440.0 * big_l * m * ru4
                    + 15300.0 * m2 * ru4
                    - 13500.0 * big_l * m2 * ru4
                    - 3600.0 * m3 * ru4
                    + 5100.0 * big_l * m3 * ru4
                    + 300.0 * m4 * ru4
                    - 900.0 * big_l * m4 * ru4
                    + 60.0 * big_l * m5 * ru4
                    - 4500.0 * r2u2
                    + 2160.0 * big_l * r2u2
                    + 6300.0 * m * r2u2
                    - 4500.0 * big_l * m * r2u2
                    - 2700.0 * m2 * r2u2
                    + 3150.0 * big_l * m2 * r2u2
                    + 360.0 * m3 * r2u2
                    - 900.0 * big_l * m3 * r2u2
                    + 90.0 * big_l * m4 * r2u2
                    + 165.0 * r3
                    - 90.0 * big_l * r3
                    - 180.0 * m * r3
                    + 165.0 * big_l * m * r3
                    + 45.0 * m2 * r3
                    - 90.0 * big_l * m2 * r3
                    + 15.0 * big_l * m3 * r3
                    - 14112.0 * u6
                    + 5760.0 * big_l * u6
                    + 25984.0 * m * u6
                    - 14112.0 * big_l * m * u6
                    - 17640.0 * m2 * u6
                    + 12992.0 * big_l * m2 * u6
                    + 5600.0 * m3 * u6
                    - 5880.0 * big_l * m3 * u6
                    - 840.0 * m4 * u6
                    + 1400.0 * big_l * m4 * u6
                    + 48.0 * m5 * u6
                    - 168.0 * big_l * m5 * u6
                    + 8.0 * big_l * m6 * u6;
                duv *= 8.0 * pr;
                duv
            }

            _ => 0.0,
        };

        duv * self.ddu[iu as usize] * self.ddv[iv as usize]
    }
}
