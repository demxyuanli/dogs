//! `HLRBRep_Curve` parameter map and projection update:
//! `Parameter2d` (`cxx:50-76`), `Parameter3d` (`cxx:78-105`), `Update`
//! (`cxx:107-224`), `UpdateMinMax` (`cxx:228-287`) and `Z` (`cxx:291-297`).

use occt_core::gp::{GpDir, GpLin, GpPnt, GpVec};
use occt_core::precision::ANGULAR;

use occt_topo::geom_abs_curve_type::GeomAbsCurveType;
use crate::hlr_algo::update_min_max as hlr_algo_update_min_max;

use super::HlrBRepCurve;

impl HlrBRepCurve {
    /// `Parameter2d(P3d)` (`cxx:50-76`): the parameter on the 2D curve from
    /// the parameter on the 3D curve. For lines the closed formula in the
    /// comment at `cxx:52-56` is used; for ellipses the offset is `myOX`;
    /// everything else is the identity.
    pub fn parameter2d(&self, p3d: f64) -> f64 {
        match self.type_ {
            GeomAbsCurveType::Line => {
                if self.perspective() {
                    let f_moz = self.of - self.oz;
                    return self.of * p3d * (self.vx * f_moz + self.ox * self.vz)
                        / (f_moz * (f_moz - p3d * self.vz));
                }
                p3d * self.vx
            }
            GeomAbsCurveType::Ellipse => p3d + self.ox,
            _ => p3d,
        }
    }

    /// `Parameter3d(P2d)` (`cxx:78-105`): the inverse of
    /// [`HlrBRepCurve::parameter2d`].
    pub fn parameter3d(&self, p2d: f64) -> f64 {
        if self.type_ == GeomAbsCurveType::Line {
            if self.perspective() {
                let f_moz = self.of - self.oz;
                return p2d * f_moz * f_moz
                    / (f_moz * (self.of * self.vx + p2d * self.vz) + self.of * self.ox * self.vz);
            }
            return if self.vx <= occt_core::precision::RESOLUTION {
                p2d
            } else {
                p2d / self.vx
            };
        } else if self.type_ == GeomAbsCurveType::Ellipse {
            return p2d - self.ox;
        }
        p2d
    }

    /// `Project(P3D, X, Y, Z)` of the projector (`HLRAlgo_Projector.cxx:262-317`)
    /// with the identity fallback used when no projector is attached.
    fn project_coords(&self, p: &GpPnt) -> (f64, f64, f64) {
        match &self.proj {
            Some(proj) => proj.project_coords(p),
            None => (p.x(), p.y(), p.z()),
        }
    }

    /// `Z(U)` (`cxx:291-297`): the Z coordinate of the curve point in the
    /// viewing coordinate system. OCCT applies `Transformation()`, not the
    /// scaled `FullTransformation()`, and does *not* apply the perspective
    /// divide.
    pub fn z(&self, u: f64) -> f64 {
        let mut p3d = self.value3d(u);
        p3d.transform(&self.transformation());
        p3d.z()
    }

    /// `Value3D(U)` (`lxx:44-50`).
    pub fn value3d(&self, u: f64) -> GpPnt {
        match &self.curve {
            Some(c) => c.d0(u),
            None => GpPnt::new(0.0, 0.0, 0.0),
        }
    }

    /// `Update(TotMin, TotMax)` (`cxx:107-224`): classify the projected curve
    /// and, when it degenerates to a line, derive the `myOX` / `myVX` /
    /// `myVZ` / `myOZ` / `myOF` scalars of the parameter map.
    pub fn update(&mut self, tot_min: &mut [f64; 16], tot_max: &mut [f64; 16]) -> f64 {
        let typ = self.curve_kind_3d();
        self.type_ = GeomAbsCurveType::OtherCurve;

        match typ {
            GeomAbsCurveType::Line => {
                self.type_ = typ;
            }
            GeomAbsCurveType::Circle => {
                if !self.perspective() {
                    let trsf = self.transformation();
                    let dz = GpDir::new(0.0, 0.0, 1.0).expect("gp::DZ()");
                    let circle = match self.curve.as_ref().and_then(|c| c.gp_circ()) {
                        Some(c) => c,
                        None => return self.update_min_max(tot_min, tot_max),
                    };
                    let mut d1 = *circle.axis().direction();
                    d1.transform(&trsf);
                    if d1.is_parallel_tol(&dz, ANGULAR) {
                        // `cxx:125-127`
                        self.type_ = GeomAbsCurveType::Circle;
                    } else if d1.dot(&dz).abs() < ANGULAR * 10.0 {
                        // `cxx:129-133`: "*10: The minor radius of ellipse
                        // should not be too small."
                        self.type_ = GeomAbsCurveType::OtherCurve;
                    } else {
                        // `cxx:134-141`: compute the angle offset.
                        self.type_ = GeomAbsCurveType::Ellipse;
                        let d3 = match d1.crossed(&dz) {
                            Ok(d) => d,
                            Err(_) => return self.update_min_max(tot_min, tot_max),
                        };
                        let mut d2 = *circle.x_axis().direction();
                        d2.transform(&trsf);
                        self.ox = d3.angle_with_ref(&d2, &d1);
                    }
                }
            }
            GeomAbsCurveType::Ellipse => {
                if !self.perspective() {
                    // `cxx:145-156`
                    let trsf = self.transformation();
                    let dz = GpDir::new(0.0, 0.0, 1.0).expect("gp::DZ()");
                    let ellipse = match self.curve.as_ref().and_then(|c| c.gp_ellipse()) {
                        Some(e) => e,
                        None => return self.update_min_max(tot_min, tot_max),
                    };
                    let mut d1 = *ellipse.axis().direction();
                    d1.transform(&trsf);
                    if d1.is_parallel_tol(&dz, ANGULAR) {
                        self.ox = 0.0; // no offset on the angle
                        self.type_ = GeomAbsCurveType::Ellipse;
                    }
                }
            }
            GeomAbsCurveType::BezierCurve => {
                // `cxx:158-168`
                let degree = match self.curve.as_ref().and_then(|c| c.bezier_poles()) {
                    Some(p) => p.len() as i32 - 1,
                    None => 0,
                };
                if degree == 1 {
                    self.type_ = GeomAbsCurveType::Line;
                } else if !self.perspective() {
                    self.type_ = typ;
                }
            }
            GeomAbsCurveType::BSplineCurve => {
                // `cxx:169-174`
                if !self.perspective() {
                    self.type_ = typ;
                }
            }
            _ => {}
        }

        if self.type_ == GeomAbsCurveType::Line {
            // `cxx:180-223`: compute the values for a line.
            let mut l3d = 1.0; // length of the 3d bezier curve
            let lin: GpLin = if typ == GeomAbsCurveType::Line {
                match self.curve.as_ref().and_then(|c| c.gp_line()) {
                    Some(l) => l,
                    None => return self.update_min_max(tot_min, tot_max),
                }
            } else {
                // bezier of degree 1 (`cxx:188-197`)
                let (pl, vl) = match &self.curve {
                    Some(c) => c.d1(0.0),
                    None => return self.update_min_max(tot_min, tot_max),
                };
                let dir = match GpDir::from_vec(&vl) {
                    Ok(d) => d,
                    Err(_) => return self.update_min_max(tot_min, tot_max),
                };
                l3d = pl.distance(&match &self.curve {
                    Some(c) => c.value(1.0),
                    None => pl,
                });
                GpLin::from_pnt_dir(pl, dir)
            };
            let mut p = lin.location();
            let mut v = GpVec::new(lin.direction().x(), lin.direction().y(), lin.direction().z());
            let trsf = self.transformation();
            p.transform(&trsf);
            v.transform(&trsf);
            if self.perspective() {
                // `cxx:199-215`
                let (f, mut vfx) = self.d1(0.0);
                let _ = vfx.normalize();
                self.vx = (vfx.x() * v.x() + vfx.y() * v.y()) * l3d;
                let l = -(vfx.x() * f.x() + vfx.y() * f.y());
                let f = occt_core::gp::GpPnt2d::new(f.x() + vfx.x() * l, f.y() + vfx.y() * l);
                self.ox = vfx.x() * (p.x() - f.x()) + vfx.y() * (p.y() - f.y());
                let focus = self.focus().unwrap_or(0.0);
                let vfz = GpVec::new(-f.x(), -f.y(), focus);
                self.of = vfz.magnitude();
                let vfz = vfz.divide(self.of);
                self.vz = vfz.dot(&v) * l3d;
                self.oz = vfz.dot(&GpVec::new(p.x() - f.x(), p.y() - f.y(), p.z()));
            } else {
                // `cxx:219-221`
                self.vx = (v.x() * v.x() + v.y() * v.y()).sqrt() * l3d;
            }
        }
        self.update_min_max(tot_min, tot_max)
    }

    /// `UpdateMinMax(TotMin, TotMax)` (`cxx:228-287`): accumulate the 3D
    /// projection of the curve into the 16-direction min/max arrays and
    /// return the chordal tolerance (`tolMinMax`).
    ///
    /// The three points `xa` / `xb` / `x` lag by one and two samples exactly
    /// as in the C++ loop: `xa <- xb`, `xb <- x`, then `x` is refreshed, which
    /// makes `dx2` the previous chord and `dx1` the two-step secant.
    #[allow(unused_assignments)]
    pub fn update_min_max(&mut self, tot_min: &mut [f64; 16], tot_max: &mut [f64; 16]) -> f64 {
        let mut a = self.first_parameter;
        let b = self.last_parameter;
        let mut tol_min_max = 0.0;

        let (mut x, mut y, mut z) = self.project_coords(&self.value3d(a));
        hlr_algo_update_min_max(x, y, z, tot_min, tot_max);

        if self.type_ != GeomAbsCurveType::Line {
            let nb_pnt: i32 = 30;
            let step = (b - a) / (nb_pnt as f64 + 1.0);
            let (mut xa, mut ya, mut za) = (0.0, 0.0, 0.0);
            let (mut xb, mut yb, mut zb) = (0.0, 0.0, 0.0);

            for i in 1..=nb_pnt {
                a += step;
                xa = xb;
                ya = yb;
                za = zb;
                xb = x;
                yb = y;
                zb = z;
                let (nx, ny, nz) = self.project_coords(&self.value3d(a));
                x = nx;
                y = ny;
                z = nz;
                hlr_algo_update_min_max(x, y, z, tot_min, tot_max);
                if i >= 2 {
                    let dx1 = x - xa;
                    let dy1 = y - ya;
                    let dz1 = z - za;
                    let dd1 = (dx1 * dx1 + dy1 * dy1 + dz1 * dz1).sqrt();
                    if dd1 > 0.0 {
                        let dx2 = xb - xa;
                        let dy2 = yb - ya;
                        let dz2 = zb - za;
                        let dd2 = (dx2 * dx2 + dy2 * dy2 + dz2 * dz2).sqrt();
                        if dd2 > 0.0 {
                            let p = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (dd1 * dd2);
                            let ex = xa + p * dx1 - xb;
                            let ey = ya + p * dy1 - yb;
                            let ez = za + p * dz1 - zb;
                            let de = (ex * ex + ey * ey + ez * ez).sqrt();
                            if de > tol_min_max {
                                tol_min_max = de;
                            }
                        }
                    }
                }
            }
        }
        let (nx, ny, nz) = self.project_coords(&self.value3d(b));
        hlr_algo_update_min_max(nx, ny, nz, tot_min, tot_max);
        tol_min_max
    }
}
