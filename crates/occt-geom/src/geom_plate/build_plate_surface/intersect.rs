//! `Intersect` (`GeomPlate_BuildPlateSurface.cxx:1909-2145`): finds the 2D
//! intersections between the constraint p-curves (each curve also against
//! itself) and books, per curve, the parameter intervals where a constraint
//! point must be dropped to avoid duplicating the same 3D constraint
//! (`PntInter`) plus the incompatible G1-G1 / G0-G1 zones (`PntG1G1`).
//!
//! Translation notes:
//! - `Ci` / `Cj` are `Geom2dAdaptor_Curve`s loaded with
//!   `myLinCont->Value(i)->Curve2dOnSurf()`; the port keeps the `Curve2d`
//!   handles (`:1926-1935`). A null p-curve is an error here where OCCT would
//!   dereference a null handle.
//! - The `Perform(Ci, myTol2d*10, myTol2d*10)` self-intersection overload
//!   (`:1936`) is `Geom2dIntGInter::perform_curve`; the two-curve overload
//!   (`:1938-1940`) is `perform_curves`. Both `TolConf` / `Tol` are
//!   `myTol2d * 10`.
//! - `Ci.Resolution(Tol)` is `Curve2d::resolution`.
//! - `std::abs(std::abs(A1) - M_PI) < myTolAng` (`:2019`, `:2079`) is copied
//!   as written, even though `A1` was already folded into `[0, pi/2]` by the
//!   `if (A1 > M_PI/2) A1 = M_PI - A1` just above; the port does not repair
//!   OCCT's own arithmetic.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_geom2d::Curve2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;

use super::BuildPlateSurface;

impl BuildPlateSurface {
    /// `myLinCont->Value(i)->Curve2dOnSurf()` with the `1..=N` to `0..N` shift
    /// and a null-handle error.
    fn p_curve(&self, i: usize) -> Result<Arc<dyn Curve2d>, String> {
        self.lin(i)
            .curve2d_on_surf()?
            .ok_or_else(|| "GeomPlate_BuildPlateSurface::Intersect: null p-curve".to_string())
    }

    /// `Intersect(PntInter, PntG1G1)` (`cxx:1911-2145`).
    ///
    /// `pnt_inter` / `pnt_g1g1` are the caller's `NCollection_HArray1<
    /// NCollection_Sequence<double>>(1, NTLinCont)` in 0-based form.
    pub(super) fn intersect(
        &self,
        pnt_inter: &mut [Vec<f64>],
        pnt_g1g1: &mut [Vec<f64>],
    ) -> Result<(), String> {
        let nt_lin_cont = self.lin_cont.len();
        let mut intersection = Geom2dIntGInter::new();
        let tol_conf = self.tol2d * 10.0;
        for i in 0..nt_lin_cont {
            let ci = self.p_curve(i)?;
            // `cxx:1930`: `for (j = i; j <= NTLinCont; j++)`, the curve
            // against itself included.
            for j in i..nt_lin_cont {
                let cj = self.p_curve(j)?;
                if i == j {
                    // `cxx:1936`: `Perform(Ci, myTol2d*10, myTol2d*10)` -- the
                    // single-curve overload, whose domain is
                    // `ComputeDomain(Ci, max(TolConf, Tol))`.
                    let domain = intersection.compute_domain(ci.as_ref(), tol_conf);
                    intersection.perform_curve(ci.as_ref(), &domain, tol_conf, tol_conf);
                } else {
                    // `cxx:1938-1940`.
                    intersection.perform_curves(ci.as_ref(), cj.as_ref(), tol_conf, tol_conf);
                }
                // `cxx:1943`.
                if intersection.result().is_empty() {
                    continue;
                }
                let nbpt = intersection.nb_points();
                for k in 0..nbpt {
                    // `cxx:1948`.
                    let int2d = intersection.point(k);
                    let par1 = int2d.param_on_first();
                    let par2 = int2d.param_on_second();
                    let p1 = self.lin(i).d0(par1);
                    let p2 = self.lin(j).d0(par2);
                    if p1.distance(&p2) < self.tol3d {
                        // `cxx:1959-1987`: the length of the interval, in 2D,
                        // corresponding in 3D to `myTol3d`.
                        let mut tolint = ci.resolution(self.tol3d);
                        let (_, v2d) = ci.d1(par1);
                        let mut aux = v2d.magnitude();
                        if aux > 1.0e-7 {
                            aux = self.tol3d / aux;
                            if aux > 100.0 * tolint {
                                tolint *= 100.0;
                            } else {
                                tolint = aux;
                            }
                        } else {
                            tolint *= 100.0;
                        }
                        pnt_inter[i].push(par1 - tolint);
                        pnt_inter[i].push(par1 + tolint);
                        // `cxx:1989`: G1-G1.
                        if self.lin(i).order() == 1 && self.lin(j).order() == 1 {
                            self.g1_g1_zone(i, j, par1, par2, pnt_g1g1)?;
                        }
                        // `cxx:2042`: G0-G1.
                        let oi = self.lin(i).order();
                        let oj = self.lin(j).order();
                        if (oi == 0 && oj == 1) || (oi == 1 && oj == 0) {
                            self.g0_g1_zone(i, j, par1, par2, pnt_g1g1)?;
                        }
                    } else {
                        // `cxx:2113-2137`: the 2D intersection maps to two 3D
                        // points further apart than `myTol3d`.
                        let dist = p1.distance(&p2);
                        let tolint = ci.resolution(dist);
                        pnt_inter[i].push(par1 - tolint);
                        pnt_inter[i].push(par1 + tolint);
                        if j != i {
                            let tolint = cj.resolution(dist);
                            pnt_inter[j].push(par2 - tolint);
                            pnt_inter[j].push(par2 + tolint);
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// `cxx:1989-2041`: the G1-G1 incompatibility test and the interval it
    /// books on both curves.
    fn g1_g1_zone(
        &self,
        i: usize,
        j: usize,
        par1: f64,
        par2: f64,
        pnt_g1g1: &mut [Vec<f64>],
    ) -> Result<(), String> {
        // `cxx:1992-1993`: the `D2` of the 3D curves behind the constraints.
        let (_, v11, v12, _v13, _v14, v15) = self.lin(i).d2(par1)?;
        let (_, v21, v22, _v23, _v24, v25) = self.lin(j).d2(par2)?;
        // `cxx:1994-1998`.
        let v16 = v11.crossed(&v12);
        let v26 = v21.crossed(&v22);
        let mut ant = v16.angle(&v26);
        if ant > PI / 2.0 {
            ant = PI - ant;
        }
        // `cxx:2001`.
        if (v16.dot(&v15) - v16.dot(&v25)).abs() > self.tol3d / 1000.0
            || ant.abs() > self.tol3d / 1000.0
        {
            // `cxx:2006-2019`: `Tol` and the 2D angle between the p-curves.
            let (tol, _) = self.g1_g1_tolerance(i, j, par1, par2, ant)?;
            // `cxx:2028-2038`.
            let coin = self.p_curve(i)?.resolution(tol);
            pnt_g1g1[i].push(par1 - coin);
            pnt_g1g1[i].push(par1 + coin);
            let coin = self.p_curve(j)?.resolution(tol);
            pnt_g1g1[j].push(par2 - coin);
            pnt_g1g1[j].push(par2 + coin);
        }
        Ok(())
    }

    /// `cxx:2007-2026`: `Tol = 100*myTol3d`, raised to `100000*myTol3d` when
    /// the two p-curves meet at `pi` within `myTolAng`. The second value is
    /// `A1` (`cxx:2012`), kept for the port's own bookkeeping.
    fn g1_g1_tolerance(
        &self,
        i: usize,
        j: usize,
        par1: f64,
        par2: f64,
        _ant: f64,
    ) -> Result<(f64, f64), String> {
        let mut tol = 100.0 * self.tol3d;
        let (_, v1) = self.p_curve(i)?.d1(par1);
        let (_, v2) = self.p_curve(j)?.d1(par2);
        let mut a1 = v1.angle(&v2);
        if a1 > PI / 2.0 {
            a1 = PI - a1;
        }
        if (a1.abs() - PI).abs() < self.tol_ang {
            tol = 100000.0 * self.tol3d;
        }
        Ok((tol, a1))
    }

    /// `cxx:2042-2112`: the G0-G1 incompatibility test and the interval booked
    /// on the G1 constraint only.
    fn g0_g1_zone(
        &self,
        i: usize,
        j: usize,
        par1: f64,
        par2: f64,
        pnt_g1g1: &mut [Vec<f64>],
    ) -> Result<(), String> {
        let oi = self.lin(i).order();
        // `cxx:2047-2057`: the G0 constraint's 3D curve tangent against the
        // normal of the G1 constraint's surface.
        let (vec, _p, vec_u, vec_v): (GpVec, GpPnt, GpVec, GpVec) = if oi == 0 {
            let curve = self
                .lin(i)
                .curve3d()
                .ok_or_else(|| "GeomPlate_BuildPlateSurface::Intersect: no 3D curve".to_string())?;
            let (p, v) = curve.d1(par1);
            let (_, u, vv) = self.lin(j).d1(par2)?;
            (v, p, u, vv)
        } else {
            let curve = self
                .lin(j)
                .curve3d()
                .ok_or_else(|| "GeomPlate_BuildPlateSurface::Intersect: no 3D curve".to_string())?;
            let (p, v) = curve.d1(par2);
            let (_, u, vv) = self.lin(i).d1(par1)?;
            (v, p, u, vv)
        };
        // `cxx:2058-2061`.
        let n = vec_u.crossed(&vec_v);
        let angle = (PI / 2.0 - vec.angle(&n)).abs();
        if angle > self.tol_ang / 10.0 {
            // `cxx:2065-2080`.
            let (tol, _) = self.g1_g1_tolerance(i, j, par1, par2, angle)?;
            // `cxx:2086-2111`.
            if oi == 1 {
                let mut coin = self.p_curve(i)?.resolution(tol);
                coin *= angle / self.tol_ang * 10.0;
                pnt_g1g1[i].push(par1 - coin);
                pnt_g1g1[i].push(par1 + coin);
            } else {
                let mut coin = self.p_curve(j)?.resolution(tol);
                coin *= angle / self.tol_ang * 10.0;
                pnt_g1g1[j].push(par2 - coin);
                pnt_g1g1[j].push(par2 + coin);
            }
        }
        Ok(())
    }
}
