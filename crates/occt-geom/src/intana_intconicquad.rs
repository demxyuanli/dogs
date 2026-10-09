//! `IntAna_IntConicQuad` -- analytic intersection between a `gp` conic
//! (line, circle, ellipse, parabola, hyperbola) and an implicit quadric
//! (`IntAna_Quadric`), with the conic/plane cases routed through
//! `IntAna_QuadQuadGeo` + `IntAna2d_AnaIntersection`.
//! Source: `IntAna_IntConicQuad.cxx` / `.hxx` / `.lxx`.

use occt_core::gp::{GpAx2d, GpAx22d, GpCirc2d, GpDir2d, GpLin2d, GpPnt2d, GpVec2d};
use occt_core::intana2d::IntAna2dAnaIntersection;
use occt_core::math_direct_poly_roots::DirectPolynomialRoots;
use occt_core::math_trig_roots::TrigonometricFunctionRoots;

use super::analytic_intersections::{plane_coeffs, quadric_quadric_planes, QuadricIntersection};
use super::prelude::*;
use super::IntAnaQuadric;

/// `PIpPI = M_PI + M_PI` (`IntAna_IntConicQuad.cxx:44`).
const PI_PI: f64 = 2.0 * std::f64::consts::PI;

/// `RealEpsilon()` (`Standard_Real.hxx:161-164`).
const REAL_EPSILON: f64 = f64::EPSILON;

/// `IntAna_IntConicQuad` (`IntAna_IntConicQuad.hxx:29-196`, `.cxx:29-575`).
#[derive(Debug, Clone, Copy)]
pub struct IntAnaIntConicQuad {
    done: bool,
    parallel: bool,
    inquadric: bool,
    nbpts: i32,
    pnts: [GpPnt; 4],
    paramonc: [f64; 4],
}

/// Signed distance from `P` to the plane (unit normal), i.e. `gp_Pln::Distance`.
fn plane_distance(p: &GpPln, pt: &GpPnt) -> f64 {
    let (a, b, c, d) = plane_coeffs(p);
    a * pt.x() + b * pt.y() + c * pt.z() + d
}

/// `Axex.Dot(V)` / `Axey.Dot(V)` for a `gp_Dir` against a `gp_Vec`.
fn dir_dot_vec(d: &GpDir, v: &GpVec) -> f64 {
    d.x() * v.x() + d.y() * v.y() + d.z() * v.z()
}

impl Default for IntAnaIntConicQuad {
    fn default() -> Self {
        Self::new()
    }
}

impl IntAnaIntConicQuad {
    /// Empty ctor (`cxx:50-57`).
    pub fn new() -> Self {
        Self {
            done: false,
            parallel: false,
            inquadric: false,
            nbpts: 0,
            pnts: [GpPnt::default(); 4],
            paramonc: [0.0; 4],
        }
    }

    // -----------------------------------------------------------------------
    // Line - quadric (`cxx:62-127`)
    // -----------------------------------------------------------------------

    /// Ctor `IntAna_IntConicQuad(const gp_Lin&, const IntAna_Quadric&)`.
    pub fn line_quadric(l: &GpLin, quad: &IntAnaQuadric) -> Self {
        let mut r = Self::new();
        r.perform_line_quadric(l, quad);
        r
    }

    /// `Perform(const gp_Lin&, const IntAna_Quadric&)` (`cxx:66-127`).
    pub fn perform_line_quadric(&mut self, l: &GpLin, quad: &IntAnaQuadric) {
        self.done = false;
        self.inquadric = false;
        self.parallel = false;
        self.nbpts = 0;

        let c = quad.coefficients();
        let (qxx, qyy, qzz) = (c[0], c[1], c[2]);
        let (qxy, qxz, qyz) = (c[3], c[4], c[5]);
        let (qx, qy, qz, qcte) = (c[6], c[7], c[8], c[9]);

        let dir = l.direction();
        let (lx, ly, lz) = (dir.x(), dir.y(), dir.z());
        let loc = l.location();
        let (lx0, ly0, lz0) = (loc.x(), loc.y(), loc.z());

        let a0 = qcte
            + qxx * lx0 * lx0
            + qyy * ly0 * ly0
            + qzz * lz0 * lz0
            + 2.0
                * (lx0 * (qx + qxy * ly0 + qxz * lz0) + ly0 * (qy + qyz * lz0) + qz * lz0);

        let a1 = 2.0
            * (lx * (qx + qxx * lx0 + qxy * ly0 + qxz * lz0)
                + ly * (qy + qxy * lx0 + qyy * ly0 + qyz * lz0)
                + lz * (qz + qxz * lx0 + qyz * ly0 + qzz * lz0));

        let a2 = qxx * lx * lx
            + qyy * ly * ly
            + qzz * lz * lz
            + 2.0 * (lx * (qxy * ly + qxz * lz) + qyz * ly * lz);

        let pol = DirectPolynomialRoots::new3(a2, a1, a0);
        if pol.is_done() {
            self.done = true;
            if pol.infinite_roots() {
                self.inquadric = true;
            } else {
                self.nbpts = pol.nb_solutions();
                for i in 1..=self.nbpts {
                    let t = pol.value(i);
                    self.paramonc[(i - 1) as usize] = t;
                    self.pnts[(i - 1) as usize] =
                        GpPnt::new(lx0 + lx * t, ly0 + ly * t, lz0 + lz * t);
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Circle - quadric (`cxx:132-196`)
    // -----------------------------------------------------------------------

    /// Ctor `IntAna_IntConicQuad(const gp_Circ&, const IntAna_Quadric&)`.
    pub fn circle_quadric(c: &GpCirc, quad: &IntAnaQuadric) -> Self {
        let mut r = Self::new();
        r.perform_circle_quadric(c, quad);
        r
    }

    /// `Perform(const gp_Circ&, const IntAna_Quadric&)` (`cxx:136-196`).
    pub fn perform_circle_quadric(&mut self, c: &GpCirc, quad: &IntAnaQuadric) {
        self.done = false;
        self.inquadric = false;
        self.parallel = false;

        let nc = quad.new_coefficients(&GpAx3::from_ax2(&c.position()));
        let (qxx, qyy) = (nc[0], nc[1]);
        let qxy = nc[3];
        let (qx, qy, qcte) = (nc[6], nc[7], nc[9]);

        let r = c.radius();
        let rr = r * r;

        let p_coscos = rr * qxx;
        let p_sinsin = rr * qyy;
        let p_sin = r * qy;
        let p_cos = r * qx;
        let p_cossin = rr * qxy;
        let p_cte = qcte;

        let pol = TrigonometricFunctionRoots::new_abcde(
            p_coscos - p_sinsin,
            p_cossin,
            p_cos + p_cos,
            p_sin + p_sin,
            p_cte + p_sinsin,
            0.0,
            PI_PI,
        );

        if pol.is_done() {
            self.done = true;
            if pol.infinite_roots() {
                self.inquadric = true;
            } else {
                self.nbpts = pol.nb_solutions();
                for i in 1..=self.nbpts {
                    let t = pol.value(i);
                    self.paramonc[(i - 1) as usize] = t;
                    self.pnts[(i - 1) as usize] = clib::circle_value(c, t);
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Elips - quadric (`cxx:200-263`)
    // -----------------------------------------------------------------------

    /// Ctor `IntAna_IntConicQuad(const gp_Elips&, const IntAna_Quadric&)`.
    pub fn ellipse_quadric(e: &GpElips, quad: &IntAnaQuadric) -> Self {
        let mut r = Self::new();
        r.perform_ellipse_quadric(e, quad);
        r
    }

    /// `Perform(const gp_Elips&, const IntAna_Quadric&)` (`cxx:204-263`).
    pub fn perform_ellipse_quadric(&mut self, e: &GpElips, quad: &IntAnaQuadric) {
        self.done = false;
        self.inquadric = false;
        self.parallel = false;

        let nc = quad.new_coefficients(&GpAx3::from_ax2(e.position()));
        let (qxx, qyy) = (nc[0], nc[1]);
        let qxy = nc[3];
        let (qx, qy, qcte) = (nc[6], nc[7], nc[9]);

        let r = e.major_radius();
        let rmin = e.minor_radius();

        let p_coscos = r * r * qxx;
        let p_sinsin = rmin * rmin * qyy;
        let p_sin = rmin * qy;
        let p_cos = r * qx;
        let p_cossin = r * rmin * qxy;
        let p_cte = qcte;

        let pol = TrigonometricFunctionRoots::new_abcde(
            p_coscos - p_sinsin,
            p_cossin,
            p_cos + p_cos,
            p_sin + p_sin,
            p_cte + p_sinsin,
            0.0,
            PI_PI,
        );

        if pol.is_done() {
            self.done = true;
            if pol.infinite_roots() {
                self.inquadric = true;
            } else {
                self.nbpts = pol.nb_solutions();
                for i in 1..=self.nbpts {
                    let t = pol.value(i);
                    self.paramonc[(i - 1) as usize] = t;
                    self.pnts[(i - 1) as usize] = clib::ellipse_value(e, t);
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Parab - quadric (`cxx:268-326`)
    // -----------------------------------------------------------------------

    /// Ctor `IntAna_IntConicQuad(const gp_Parab&, const IntAna_Quadric&)`.
    pub fn parabola_quadric(p: &GpParab, quad: &IntAnaQuadric) -> Self {
        let mut r = Self::new();
        r.perform_parabola_quadric(p, quad);
        r
    }

    /// `Perform(const gp_Parab&, const IntAna_Quadric&)` (`cxx:272-326`).
    pub fn perform_parabola_quadric(&mut self, p: &GpParab, quad: &IntAnaQuadric) {
        self.done = false;
        self.inquadric = false;
        self.parallel = false;

        let nc = quad.new_coefficients(&GpAx3::from_ax2(p.position()));
        let (qxx, qyy) = (nc[0], nc[1]);
        let qxy = nc[3];
        let (qx, qy, qcte) = (nc[6], nc[7], nc[9]);

        let f = p.focal;
        let un_sur_2p = 0.25 / f;

        let a4 = qxx * un_sur_2p * un_sur_2p;
        let a3 = (qxy + qxy) * un_sur_2p;
        let a2 = qyy + (qx + qx) * un_sur_2p;
        let a1 = qy + qy;
        let a0 = qcte;

        let pol = DirectPolynomialRoots::new5(a4, a3, a2, a1, a0);
        if pol.is_done() {
            self.done = true;
            if pol.infinite_roots() {
                self.inquadric = true;
            } else {
                self.nbpts = pol.nb_solutions();
                for i in 1..=self.nbpts {
                    let t = pol.value(i);
                    self.paramonc[(i - 1) as usize] = t;
                    self.pnts[(i - 1) as usize] = clib::parabola_value(p, t);
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Hypr - quadric (`cxx:330-397`)
    // -----------------------------------------------------------------------

    /// Ctor `IntAna_IntConicQuad(const gp_Hypr&, const IntAna_Quadric&)`.
    pub fn hyperbola_quadric(h: &GpHypr, quad: &IntAnaQuadric) -> Self {
        let mut r = Self::new();
        r.perform_hyperbola_quadric(h, quad);
        r
    }

    /// `Perform(const gp_Hypr&, const IntAna_Quadric&)` (`cxx:334-397`).
    pub fn perform_hyperbola_quadric(&mut self, h: &GpHypr, quad: &IntAnaQuadric) {
        self.done = false;
        self.inquadric = false;
        self.parallel = false;

        let nc = quad.new_coefficients(&GpAx3::from_ax2(h.position()));
        let (qxx, qyy) = (nc[0], nc[1]);
        let qxy = nc[3];
        let (qx, qy, qcte) = (nc[6], nc[7], nc[9]);

        let r = h.major_radius;
        let rmin = h.minor_radius;
        let rr = r * r;
        let rminrmin = rmin * rmin;
        let rr_prod = r * rmin;

        let a4 = rr * qxx + rr_prod * (qxy + qxy) + rminrmin * qyy;
        let a3 = 4.0 * (r * qx + rmin * qy);
        let a2 = 2.0 * ((qcte + qcte) + qxx * rr - qyy * rminrmin);
        let a1 = 4.0 * (r * qx - rmin * qy);
        let a0 = qxx * rr - rr_prod * (qxy + qxy) + qyy * rminrmin;

        let pol = DirectPolynomialRoots::new5(a4, a3, a2, a1, a0);
        if pol.is_done() {
            self.done = true;
            if pol.infinite_roots() {
                self.inquadric = true;
            } else {
                self.nbpts = pol.nb_solutions();
                let mut bonnesolutions = 0i32;
                for i in 1..=self.nbpts {
                    let t = pol.value(i);
                    if t >= REAL_EPSILON {
                        let lnt = t.ln();
                        self.paramonc[bonnesolutions as usize] = lnt;
                        self.pnts[bonnesolutions as usize] = clib::hyperbola_value(h, lnt);
                        bonnesolutions += 1;
                    }
                }
                self.nbpts = bonnesolutions;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Line - plane (`cxx:436-492`)
    // -----------------------------------------------------------------------

    /// Ctor `IntAna_IntConicQuad(const gp_Lin&, const gp_Pln&, Tolang, Tol, Len)`.
    pub fn line_plane(l: &GpLin, p: &GpPln, tol_ang: f64, tol: f64, len: f64) -> Self {
        let mut r = Self::new();
        r.perform_line_plane(l, p, tol_ang, tol, len);
        r
    }

    /// `Perform(const gp_Lin&, const gp_Pln&, Tolang, Tol, Len)` (`cxx:436-492`).
    pub fn perform_line_plane(
        &mut self,
        l: &GpLin,
        p: &GpPln,
        tol_ang: f64,
        tol: f64,
        len: f64,
    ) {
        self.done = false;

        let (a, b, c, d) = plane_coeffs(p);
        let orig = l.location();
        let dir = l.direction();
        let (al, bl, cl) = (dir.x(), dir.y(), dir.z());

        let direc = a * al + b * bl + c * cl;
        let dis = a * orig.x() + b * orig.y() + c * orig.z() + d;

        self.parallel = false;
        if direc.abs() < tol_ang {
            self.parallel = true;
            if len != 0.0 && direc != 0.0 {
                // Check the distance from the bounding point of the line to the plane.
                let p1 = GpPnt::new(orig.x() - dis * a, orig.y() - dis * b, orig.z() - dis * c);
                let p2 = GpPnt::new(p1.x() + len * al, p1.y() + len * bl, p1.z() + len * cl);
                if plane_distance(p, &p2) > tol {
                    self.parallel = false;
                }
            }
        }
        if self.parallel {
            self.inquadric = dis.abs() < tol_ang;
        } else {
            self.parallel = false;
            self.inquadric = false;
            self.nbpts = 1;
            self.paramonc[0] = -dis / direc;
            self.pnts[0] = GpPnt::new(
                orig.x() + self.paramonc[0] * al,
                orig.y() + self.paramonc[0] * bl,
                orig.z() + self.paramonc[0] * cl,
            );
        }
        self.done = true;
    }

    // -----------------------------------------------------------------------
    // Circle - plane (`cxx:494-560`)
    // -----------------------------------------------------------------------

    /// Ctor `IntAna_IntConicQuad(const gp_Circ&, const gp_Pln&, Tolang, Tol)`.
    pub fn circle_plane(c: &GpCirc, p: &GpPln, tol_ang: f64, tol: f64) -> Self {
        let mut r = Self::new();
        r.perform_circle_plane(c, p, tol_ang, tol);
        r
    }

    /// `Perform(const gp_Circ&, const gp_Pln&, Tolang, Tol)` (`cxx:494-560`).
    pub fn perform_circle_plane(&mut self, c: &GpCirc, p: &GpPln, tol_ang: f64, tol: f64) {
        self.done = false;

        let pos = c.position();
        let plconic = GpPln::new(GpAx3::from_ax2(&pos));
        let intp = quadric_quadric_planes(&plconic, p, tol_ang, tol);
        let ligsol = match intp {
            QuadricIntersection::None => {
                self.parallel = true;
                let distmax = plane_distance(p, &c.location()) + c.radius() * tol_ang;
                self.inquadric = distmax < tol;
                self.done = true;
                return;
            }
            QuadricIntersection::Same => {
                self.inquadric = true;
                self.done = true;
                return;
            }
            QuadricIntersection::Line(l) => l,
            _ => return,
        };

        self.inquadric = false;
        self.parallel = false;

        let v0 = GpVec::from_pnts(&plconic.location(), &ligsol.location());
        let axex = *plconic.position().x_direction();
        let axey = *plconic.position().y_direction();

        let orig = GpPnt2d::new(dir_dot_vec(&axex, &v0), dir_dot_vec(&axey, &v0));
        let ligsol_dir = ligsol.direction();
        let ligsol_vec = GpVec::new(ligsol_dir.x(), ligsol_dir.y(), ligsol_dir.z());
        let dire = GpVec2d::new(
            dir_dot_vec(&axex, &ligsol_vec),
            dir_dot_vec(&axey, &ligsol_vec),
        );

        let Some(dir2d) = GpDir2d::new(dire.x(), dire.y()).ok() else {
            return;
        };
        let ligs = GpLin2d::new(GpAx2d::new(orig, dir2d));

        let ax2d_bid = GpAx2d::new(GpPnt2d::zero(), GpDir2d::default());
        let cir = GpCirc2d::new(GpAx22d::from_xdir(ax2d_bid.loc, ax2d_bid.vdir), c.radius());

        let mut int2d = IntAna2dAnaIntersection::new();
        int2d.perform_lin_circ(&ligs, &cir);
        if !int2d.is_done() {
            return;
        }

        let n = int2d.nb_points().min(4);
        self.nbpts = n;
        let loc = plconic.location();
        for i in 1..=n {
            let resul = *int2d.point(i).value();
            let (x, y) = (resul.x(), resul.y());
            self.pnts[(i - 1) as usize] = GpPnt::new(
                loc.x() + x * axex.x() + y * axey.x(),
                loc.y() + x * axex.y() + y * axey.y(),
                loc.z() + x * axex.z() + y * axey.z(),
            );
            self.paramonc[(i - 1) as usize] = int2d.point(i).param_on_second();
        }
        self.done = true;
    }

    // -----------------------------------------------------------------------
    // Accessors (`IntAna_IntConicQuad.lxx:19-89`)
    // -----------------------------------------------------------------------

    /// `IsDone` (`lxx:19-22`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `IsInQuadric` (`lxx:24-31`).
    pub fn is_in_quadric(&self) -> bool {
        assert!(self.done, "StdFail_NotDone in IntAna_IntConicQuad::IsInQuadric");
        self.inquadric
    }

    /// `IsParallel` (`lxx:33-40`).
    pub fn is_parallel(&self) -> bool {
        assert!(self.done, "StdFail_NotDone in IntAna_IntConicQuad::IsParallel");
        self.parallel
    }

    /// `NbPoints` (`lxx:42-56`).
    pub fn nb_points(&self) -> i32 {
        assert!(self.done, "StdFail_NotDone in IntAna_IntConicQuad::NbPoints");
        assert!(
            !(self.parallel || self.inquadric),
            "Standard_DomainError in IntAna_IntConicQuad::NbPoints"
        );
        self.nbpts
    }

    /// `Point(N)` (`lxx:58-71`), 1-based.
    pub fn point(&self, index: i32) -> &GpPnt {
        assert!(self.done, "StdFail_NotDone in IntAna_IntConicQuad::Point");
        assert!(
            !(self.parallel || self.inquadric),
            "Standard_DomainError in IntAna_IntConicQuad::Point"
        );
        assert!(
            index > 0 && index <= self.nbpts,
            "Standard_OutOfRange in IntAna_IntConicQuad::Point"
        );
        &self.pnts[(index - 1) as usize]
    }

    /// `ParamOnConic(N)` (`lxx:73-89`), 1-based.
    pub fn param_on_conic(&self, index: i32) -> f64 {
        assert!(self.done, "StdFail_NotDone in IntAna_IntConicQuad::ParamOnConic");
        assert!(
            !(self.parallel || self.inquadric),
            "Standard_DomainError in IntAna_IntConicQuad::ParamOnConic"
        );
        assert!(
            index > 0 && index <= self.nbpts,
            "Standard_OutOfRange in IntAna_IntConicQuad::ParamOnConic"
        );
        self.paramonc[(index - 1) as usize]
    }
}
