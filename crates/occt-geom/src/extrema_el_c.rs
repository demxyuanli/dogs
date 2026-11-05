//! `Extrema_ExtElC` — distances between two elementary curves.
//!
//! Source: `ModelingData/TKGeomBase/Extrema/Extrema_ExtElC.hxx:36-92`,
//! `.cxx:256-1168` (OCCT 8.0.0).
//!
//! Ported constructors:
//!
//! | pair | OCCT | port |
//! |---|---|---|
//! | `Lin` / `Lin` | `.cxx:268-358` | [`ExtremaExtElC::lines`] |
//! | `Lin` / `Circ` | `.cxx:471-625` | [`ExtremaExtElC::line_circle`] |
//! | `Circ` / `Circ` | `.cxx:955-1114` | [`ExtremaExtElC::circles`] |
//!
//! The remaining three (line/ellipse `.cxx:627-754`, line/hyperbola
//! `.cxx:757-859`, line/parabola `.cxx:862-952`) are **UNPORTED**; nothing in
//! the port calls them yet.
//!
//! `Lin`/`Circ` needs the internal `ExtremaExtElC_TrigonometricRoots`
//! (`.cxx:41-234`), the `IntAna_IntQuadQuad` helper duplicate, and
//! `PlanarLineCircleExtrema` (`.cxx:361-469`), which solves the coplanar case
//! with `Extrema_ExtElC2d` + `IntAna2d_AnaIntersection`. `Circ`/`Circ` also uses
//! the file-local `RefineDir` (`.cxx:1171-1196`).
//!
//! `Extrema_POnCurv` (`Extrema_POnCurv.hxx:24-60`) is defined here because this
//! is its only consumer: a parameter on a curve plus the point at that
//! parameter.

use occt_core::elib::clib;
use occt_core::gp::{
    GpAx22d, GpAx3, GpCirc, GpCirc2d, GpDir, GpDir2d, GpLin, GpLin2d, GpPln, GpPnt, GpPnt2d, GpXyz,
};
use occt_core::intana2d::IntAna2dAnaIntersection;
use occt_core::math_trig_roots::TrigonometricFunctionRoots;
use occt_core::precision::{Precision, ANGULAR, CONFUSION, REAL_SMALL, RESOLUTION};
use occt_geom2d::extrema2d::ExtremaExtElC2d;

/// `RealLast()` (`Standard_Real.hxx:110-113`) — the initial `mySqDist` value
/// (`cxx:261-264`).
const REAL_LAST: f64 = f64::MAX;

/// `RealEpsilon()` (`Standard_Real.hxx:124-128`) — `DBL_EPSILON`.
const REAL_EPSILON: f64 = f64::EPSILON;

/// `gp_Pln::SquareDistance(theP)` (`gp_Pln.hxx:159-163`).
fn pln_square_distance_pnt(pl: &GpPln, p: &GpPnt) -> f64 {
    let d = pl.signed_distance_pnt(p);
    d * d
}

/// `static void RefineDir(gp_Dir&)` (`cxx:1171-1196`): snap a direction whose
/// component is within `RealEpsilon()` of `+-1` onto that axis.
fn refine_dir(dir: &GpDir) -> GpDir {
    let mut cx = [dir.x(), dir.y(), dir.z()];
    for i in 0..3 {
        let one = if cx[i] > 0.0 { 1.0 } else { -1.0 };
        if cx[i] > one - REAL_EPSILON && cx[i] < one + REAL_EPSILON {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            cx[i] = one;
            cx[j] = 0.0;
            cx[k] = 0.0;
            return GpDir::new(cx[0], cx[1], cx[2]).expect("RefineDir: unit direction");
        }
    }
    *dir
}

/// `class ExtremaExtElC_TrigonometricRoots` (`cxx:41-234`) — the
/// `IntAna_IntQuadQuad` helper duplicate: sorted roots of
/// `CC*cos^2 + SC*cos*sin + C*cos + S*sin + Cte`, with the near-zero-coefficient
/// retry of `cxx:214-233` and the `aPrecision` root filter of `cxx:158-181`.
///
/// `IsARoot` (`cxx:70-90`) has no caller in `Extrema_ExtElC.cxx` and is not
/// ported.
#[derive(Debug, Clone)]
struct ExtremaExtElCTrigRoots {
    /// `double Roots[4]` (`cxx:49`).
    roots: [f64; 4],
    /// `bool done` (`cxx:50`).
    done: bool,
    /// `int NbRoots` (`cxx:51`).
    nb_roots: usize,
    /// `bool infinite_roots` (`cxx:52`).
    infinite_roots: bool,
}

impl ExtremaExtElCTrigRoots {
    /// `ExtremaExtElC_TrigonometricRoots(CC, SC, C, S, Cte, Binf, Bsup)`
    /// (`cxx:121-234`).
    fn new(cc: f64, sc: f64, c: f64, s: f64, cte: f64, binf: f64, bsup: f64) -> Self {
        let mut out = Self {
            roots: [0.0; 4],
            done: false,
            nb_roots: 0,
            infinite_roots: false,
        };
        let mut lcc = cc;
        let mut lsc = sc;
        let mut lc = c;
        let mut ls = s;
        let mut lcte = cte;
        let mut nbessai = 1;
        while nbessai <= 2 && !out.done {
            let mtfr = TrigonometricFunctionRoots::new_abcde(lcc, lsc, lc, ls, lcte, binf, bsup);
            if mtfr.is_done() {
                out.done = true;
                if mtfr.infinite_roots() {
                    out.infinite_roots = true;
                } else {
                    // `cxx:146-181`.
                    let two_pi = std::f64::consts::PI + std::f64::consts::PI;
                    out.nb_roots = mtfr.nb_solutions() as usize;
                    for i in 0..out.nb_roots {
                        let mut v = mtfr.value((i + 1) as i32);
                        if v < 0.0 {
                            v += two_pi;
                        }
                        if v > two_pi {
                            v -= two_pi;
                        }
                        out.roots[i] = v;
                    }
                    // `cxx:153-157`: `std::max` of the raw coefficients, not
                    // their magnitudes.
                    let a_max_coef = f64::max(f64::max(f64::max(cc, sc), f64::max(c, s)), cte);
                    let a_precision = f64::max(1.0e-8, 1.0e-12 * a_max_coef);

                    let sv_nb_roots = out.nb_roots;
                    for i in 0..sv_nb_roots {
                        let co = out.roots[i].cos();
                        let si = out.roots[i].sin();
                        let y = co * (cc * co + (sc + sc) * si + c) + s * si + cte;
                        if y.abs() > a_precision {
                            out.nb_roots -= 1;
                            out.roots[i] = 1000.0;
                        }
                    }
                    // `cxx:183-197`: bubble sort over the untouched slot count.
                    loop {
                        let mut triee = true;
                        for i in 1..sv_nb_roots {
                            if out.roots[i] < out.roots[i - 1] {
                                triee = false;
                                out.roots.swap(i, i - 1);
                            }
                        }
                        if triee {
                            break;
                        }
                    }
                    // `cxx:198-208`.
                    out.infinite_roots = false;
                    if out.nb_roots == 0
                        && (cc.abs() + sc.abs() + c.abs() + s.abs()) < 1e-10
                        && cte.abs() < 1e-10
                    {
                        out.infinite_roots = true;
                    }
                }
            } else {
                // `cxx:214-233`: retry once with near-zero coefficients zeroed.
                if cc.abs() < 1e-10 {
                    lcc = 0.0;
                }
                if sc.abs() < 1e-10 {
                    lsc = 0.0;
                }
                if c.abs() < 1e-10 {
                    lc = 0.0;
                }
                if s.abs() < 1e-10 {
                    ls = 0.0;
                }
                if cte.abs() < 1e-10 {
                    lcte = 0.0;
                }
                nbessai += 1;
            }
        }
        out
    }

    /// `IsDone()` (`cxx:63`).
    fn is_done(&self) -> bool {
        self.done
    }

    /// `InfiniteRoots()` (`cxx:95-102`).
    fn infinite_roots(&self) -> bool {
        assert!(self.done, "StdFail_NotDone: ExtremaExtElC_TrigonometricRoots::InfiniteRoots");
        self.infinite_roots
    }

    /// `NbSolutions()` (`cxx:83-93`).
    fn nb_solutions(&self) -> usize {
        assert!(self.done, "StdFail_NotDone: ExtremaExtElC_TrigonometricRoots::NbSolutions");
        self.nb_roots
    }

    /// `Value(n)` (`cxx:105-113`), 1-based.
    fn value(&self, n: usize) -> f64 {
        assert!(
            self.done && n <= self.nb_roots,
            "StdFail_NotDone: ExtremaExtElC_TrigonometricRoots::Value"
        );
        self.roots[n - 1]
    }
}

/// `Extrema_POnCurv` (`Extrema_POnCurv.hxx:24-60`): `myU` plus `myP`, plus the
/// 1-based readers `Value()` / `Parameter()`.
#[derive(Debug, Clone, Copy, Default)]
pub struct ExtremaPOnCurv {
    /// `myU` (`hxx:58`).
    pub u: f64,
    /// `myP` (`hxx:59`).
    pub p: GpPnt,
}

impl ExtremaPOnCurv {
    /// `Extrema_POnCurv(theU, theP)` (`hxx:38-42`).
    pub fn new(u: f64, p: GpPnt) -> Self {
        Self { u, p }
    }

    /// `Parameter()` (`hxx:56-57`).
    pub fn parameter(&self) -> f64 {
        self.u
    }

    /// `Value()` (`hxx:52-53`).
    pub fn value(&self) -> &GpPnt {
        &self.p
    }

    /// `SetValues(theU, theP)` (`hxx:44-48`).
    pub fn set_values(&mut self, u: f64, p: GpPnt) {
        self.u = u;
        self.p = p;
    }
}

/// `class Extrema_ExtElC` (`hxx:36-92`): the six `myPoint` slots of
/// `hxx:92` are the largest array OCCT declares for any of its constructors;
/// `Lin`/`Circ` fills up to four and `Circ`/`Circ` up to six.
#[derive(Debug, Clone)]
pub struct ExtremaExtElC {
    /// `myDone` (`hxx:88`).
    done: bool,
    /// `myIsPar` (`hxx:89`).
    is_par: bool,
    /// `myNbExt` (`hxx:90`).
    nb_ext: usize,
    /// `mySqDist[6]` (`hxx:91`).
    sq_dist: [f64; 6],
    /// `myPoint[6][2]` (`hxx:92`).
    point: [[ExtremaPOnCurv; 2]; 6],
}

impl ExtremaExtElC {
    /// `Extrema_ExtElC()` (`cxx:256-265`): `mySqDist` starts at `RealLast()`.
    pub fn new() -> Self {
        Self {
            done: false,
            is_par: false,
            nb_ext: 0,
            sq_dist: [REAL_LAST; 6],
            point: [[ExtremaPOnCurv::default(); 2]; 6],
        }
    }

    /// `Extrema_ExtElC(const gp_Lin& C1, const gp_Lin& C2, const double)`
    /// (`cxx:268-358`). The tolerance argument is unnamed and unused in OCCT
    /// (the parallel test uses `gp::Resolution()` and `Precision::Angular()`),
    /// so the port takes no tolerance.
    pub fn lines(c1: &GpLin, c2: &GpLin) -> Self {
        let mut this = Self::new();

        // `const gp_Dir &aD1 = theC1.Position().Direction(), &aD2 = ...;
        // const double aCosA = aD1.Dot(aD2); const double aSqSinA = 1.0 -
        // aCosA * aCosA;` (`cxx:324-327`).
        let d1 = c1.direction();
        let d2 = c2.direction();
        let cos_a = d1.dot(&d2);
        let sq_sin_a = 1.0 - cos_a * cos_a;
        let mut u1 = 0.0f64;
        let mut u2 = 0.0f64;

        // `if (aSqSinA < gp::Resolution() || aD1.IsParallel(aD2,
        // Precision::Angular())) { myIsPar = true; } else { ... }`
        // (`cxx:328-339`).
        if sq_sin_a < REAL_SMALL || d1.is_parallel_tol(&d2, ANGULAR) {
            this.is_par = true;
        } else {
            let l = c2.location().xyz().subtracted(c1.location().xyz());
            let d1l = d1.xyz().dot(&l);
            let d2l = d2.xyz().dot(&l);
            u1 = (d1l - cos_a * d2l) / sq_sin_a;
            u2 = (cos_a * d1l - d2l) / sq_sin_a;
            this.is_par = Precision::is_infinite(u1) || Precision::is_infinite(u2);
        }

        // `if (myIsPar) { mySqDist[0] = theC2.SquareDistance(theC1.Location());
        // myNbExt = 1; myDone = true; return; }` (`cxx:341-348`).
        if this.is_par {
            this.sq_dist[0] = c2.square_distance(&c1.location());
            this.nb_ext = 1;
            this.done = true;
            return this;
        }

        // `const gp_Pnt aP1(ElCLib::Value(aU1, theC1)), aP2(ElCLib::Value(aU2,
        // theC2)); mySqDist[myNbExt] = aP1.SquareDistance(aP2);
        // myPoint[myNbExt][0] = Extrema_POnCurv(aU1, aP1); myPoint[myNbExt][1]
        // = Extrema_POnCurv(aU2, aP2); myNbExt = 1; myDone = true;`
        // (`cxx:352-357`).
        let p1 = clib::line_value(c1, u1);
        let p2 = clib::line_value(c2, u2);
        this.sq_dist[0] = p1.square_distance(&p2);
        this.point[0] = [ExtremaPOnCurv::new(u1, p1), ExtremaPOnCurv::new(u2, p2)];
        this.nb_ext = 1;
        this.done = true;
        this
    }

    /// `PlanarLineCircleExtrema(theLin, theCirc)` (`cxx:361-469`): the coplanar
    /// line/circle case, solved in the circle's 2D frame. Returns `true` when the
    /// caller must return immediately, exactly like the OCCT `bool`.
    fn planar_line_circle_extrema(&mut self, the_lin: &GpLin, the_circ: &GpCirc) -> bool {
        let pos = the_circ.position();
        let a_dir_c = pos.direction();
        let a_dir_l = the_lin.direction();
        if a_dir_c.dot(&a_dir_l).abs() > ANGULAR {
            return false;
        }

        // `cxx:383-407`: project the line into the circle frame.
        let a_c_loc = *the_circ.location().xyz();
        let a_d_cx = *pos.x_direction().xyz();
        let a_d_cy = *pos.y_direction().xyz();
        let a_l_loc = *the_lin.location().xyz();
        let a_l_dir = *the_lin.direction().xyz();
        let a_vec_cl = a_l_loc.subtracted(&a_c_loc);

        let a_pc = GpPnt2d::new(0.0, 0.0);
        let a_circ_axis = GpAx22d::new(
            a_pc,
            GpDir2d::new(1.0, 0.0).expect("gp::DX2d"),
            GpDir2d::new(0.0, 1.0).expect("gp::DY2d"),
        )
        .expect("gp_Ax22d(gp::DX2d, gp::DY2d)");
        let a_circ2d = GpCirc2d::new(a_circ_axis, the_circ.radius());
        let a_pl = GpPnt2d::new(a_vec_cl.dot(&a_d_cx), a_vec_cl.dot(&a_d_cy));
        let a_dl = GpDir2d::new(a_l_dir.dot(&a_d_cx), a_l_dir.dot(&a_d_cy))
            .expect("gp_Dir2d: the line direction cannot vanish in its own plane");
        let a_lin2d = GpLin2d::from_pnt_dir(a_pl, a_dl);

        // `cxx:410-414`.
        let an_ext2d = ExtremaExtElC2d::new_line_circle(&a_lin2d, &a_circ2d);
        let mut an_inters = IntAna2dAnaIntersection::new();
        an_inters.perform_lin_circ(&a_lin2d, &a_circ2d);

        self.done = an_ext2d.is_done() || an_inters.is_done();
        if !self.done {
            return true;
        }

        let a_nb_extr = an_ext2d.nb_ext();
        let a_nb_sol = an_inters.nb_points() as usize;

        // `cxx:425-464`: the 2D extrema first, then the intersection points.
        for an_extr_id in 1..=(a_nb_extr + a_nb_sol) {
            let (a_lin_par, a_circ_par) = if an_extr_id <= a_nb_extr {
                let (p_lin2d, p_circ2d) = an_ext2d.points(an_extr_id);
                (p_lin2d.parameter(), p_circ2d.parameter())
            } else {
                let ip = an_inters.point((an_extr_id - a_nb_extr) as i32);
                (ip.param_on_first(), ip.param_on_second())
            };

            let a_p_on_l = clib::line_value(the_lin, a_lin_par);
            let a_p_on_c = clib::circle_value(the_circ, a_circ_par);

            let i = self.nb_ext;
            self.sq_dist[i] = a_p_on_l.square_distance(&a_p_on_c);
            self.point[i] = [
                ExtremaPOnCurv::new(a_lin_par, a_p_on_l),
                ExtremaPOnCurv::new(a_circ_par, a_p_on_c),
            ];
            self.nb_ext += 1;
        }

        true
    }

    /// `Extrema_ExtElC(const gp_Lin& C1, const gp_Circ& C2, const double)`
    /// (`cxx:471-625`). The tolerance argument is unnamed and unused in OCCT.
    pub fn line_circle(c1: &GpLin, c2: &GpCirc) -> Self {
        let mut this = Self::new();

        // `cxx:483-486`.
        if this.planar_line_circle_extrema(c1, c2) {
            return this;
        }

        // `cxx:488-498`: the line direction in the reference of the circle.
        let d1 = c1.direction();
        let pos2 = c2.position();
        let x2 = *pos2.x_direction();
        let y2 = *pos2.y_direction();
        let z2 = pos2.direction();
        let dx = d1.dot(&x2);
        let dy = d1.dot(&y2);
        let dz = d1.dot(&z2);
        let d = refine_dir(&GpDir::new(dx, dy, dz).expect("gp_Dir: the line direction is unit"));
        let (dx, dy, _dz) = (d.x(), d.y(), d.z());

        // `cxx:501-521`: `O2O1` in the reference frame of the circle.
        let o1 = c1.location();
        let o2 = c2.location();
        let mut o2o1 = o1.xyz().subtracted(o2.xyz());
        let a_tol_ro2o1 = RESOLUTION;
        let a_ro2o1 = o2o1.modulus();
        if a_ro2o1 > a_tol_ro2o1 {
            o2o1 = o2o1.multiplied(1.0 / a_ro2o1);
            let a_d_o2o1 = GpDir::new(
                o2o1.dot(x2.xyz()),
                o2o1.dot(y2.xyz()),
                o2o1.dot(z2.xyz()),
            )
            .expect("gp_Dir: O2O1 is unit");
            o2o1 = refine_dir(&a_d_o2o1).xyz().multiplied(a_ro2o1);
        } else {
            o2o1 = GpXyz::new(
                o2o1.dot(x2.xyz()),
                o2o1.dot(y2.xyz()),
                o2o1.dot(z2.xyz()),
            );
        }

        // `cxx:523`: `Vxyz = (D.XYZ() * (O2O1.Dot(D))) - O2O1.XYZ()`.
        let vxyz = d.xyz().multiplied(o2o1.dot(d.xyz())).subtracted(&o2o1);

        // `cxx:559-580`: coefficients of the equation in `cos` and `sin`,
        // divided by `R`, with the `1.e-12` zeroing of `cxx:561-578`.
        let a_tol = 1.0e-12;
        let r = c2.radius();
        let a5 = r * dx * dy;
        let mut a1 = -2.0 * a5;
        let mut a2 = 0.5 * r * (dx * dx - dy * dy);
        let mut a3 = vxyz.y();
        let mut a4 = -vxyz.x();
        let mut a5 = a5;
        if a1 >= -a_tol && a1 <= a_tol {
            a1 = 0.0;
        }
        if a2 >= -a_tol && a2 <= a_tol {
            a2 = 0.0;
        }
        if a3 >= -a_tol && a3 <= a_tol {
            a3 = 0.0;
        }
        if a4 >= -a_tol && a4 <= a_tol {
            a4 = 0.0;
        }
        if a5 >= -a_tol && a5 <= a_tol {
            a5 = 0.0;
        }

        // `cxx:588-601`.
        let two_pi = std::f64::consts::PI + std::f64::consts::PI;
        let sol = ExtremaExtElCTrigRoots::new(a1, a2, a3, a4, a5, 0.0, two_pi);
        if !sol.is_done() {
            return this;
        }
        if sol.infinite_roots() {
            this.is_par = true;
            this.sq_dist[0] = r * r;
            this.nb_ext = 1;
            this.done = true;
            return this;
        }

        // `cxx:603-621`.
        let nb_sol = sol.nb_solutions();
        for no_sol in 1..=nb_sol {
            let u2 = sol.value(no_sol);
            let p2 = clib::circle_value(c2, u2);
            let u1 = p2.xyz().subtracted(o1.xyz()).dot(d1.xyz());
            let p1 = clib::line_value(c1, u1);
            let i = this.nb_ext;
            this.sq_dist[i] = p1.square_distance(&p2);
            this.point[i] = [ExtremaPOnCurv::new(u1, p1), ExtremaPOnCurv::new(u2, p2)];
            this.nb_ext += 1;
        }
        this.done = true;
        this
    }

    /// `Extrema_ExtElC(const gp_Circ& C1, const gp_Circ& C2)` (`cxx:955-1114`).
    pub fn circles(c1: &GpCirc, c2: &GpCirc) -> Self {
        let mut this = Self::new();

        let a_tol_a = ANGULAR;
        let a_tol_d = CONFUSION;
        let a_tol_d2 = a_tol_d * a_tol_d;

        // `cxx:976-989`.
        let a_pc1 = c1.location();
        let a_dc1 = *c1.axis().direction();
        let a_pc2 = c2.location();
        let a_dc2 = *c2.axis().direction();
        let a_plc1 = GpPln::new(GpAx3::from_point_dir(a_pc1, a_dc1).expect("gp_Ax3(aPc1, aDc1)"));

        let a_d2 = pln_square_distance_pnt(&a_plc1, &a_pc2);
        let b_is_same_plane = a_dc1.is_parallel_tol(&a_dc2, a_tol_a) && a_d2 < a_tol_d2;
        if !b_is_same_plane {
            return this;
        }

        // `cxx:995-1007`: the same axis.
        let a_dc2_sq = a_pc1.square_distance(&a_pc2);
        if a_dc2_sq < a_tol_d2 {
            this.is_par = true;
            this.nb_ext = 1;
            this.done = true;
            let a_dr = c1.radius() - c2.radius();
            this.sq_dist[0] = a_dr * a_dr;
            return this;
        }

        this.done = true;

        // `cxx:1013-1033`: `aC1` is the larger circle.
        let (a_c1, a_c2, j1, j2): (&GpCirc, &GpCirc, usize, usize) = if c2.radius() > c1.radius() {
            (c2, c1, 1, 0)
        } else {
            (c1, c2, 0, 1)
        };
        let a_r1 = a_c1.radius();
        let a_r2 = a_c2.radius();
        let a_pc1 = a_c1.location();
        let a_pc2 = a_c2.location();

        let a_d12 = a_pc1.distance(&a_pc2);
        let a_dir12 = GpDir::new(
            a_pc2.x() - a_pc1.x(),
            a_pc2.y() - a_pc1.y(),
            a_pc2.z() - a_pc1.z(),
        )
        .expect("gp_Dir(aVec12): the centres are apart by more than a confusion");

        // `cxx:1037-1072`: four common solutions.
        this.nb_ext = 4;
        let a_p11 = GpPnt::from_xyz(&a_pc1.xyz().subtracted(&a_dir12.xyz().multiplied(a_r1)));
        let a_p12 = GpPnt::from_xyz(&a_pc1.xyz().added(&a_dir12.xyz().multiplied(a_r1)));
        let a_p21 = GpPnt::from_xyz(&a_pc2.xyz().subtracted(&a_dir12.xyz().multiplied(a_r2)));
        let a_p22 = GpPnt::from_xyz(&a_pc2.xyz().added(&a_dir12.xyz().multiplied(a_r2)));

        let a_t11 = clib::circle_parameter(&a_c1.position(), &a_p11);
        let a_t12 = clib::circle_parameter(&a_c1.position(), &a_p12);
        let a_t21 = clib::circle_parameter(&a_c2.position(), &a_p21);
        let a_t22 = clib::circle_parameter(&a_c2.position(), &a_p22);

        this.point[0][j1] = ExtremaPOnCurv::new(a_t11, a_p11);
        this.point[0][j2] = ExtremaPOnCurv::new(a_t21, a_p21);
        this.sq_dist[0] = a_p11.square_distance(&a_p21);

        this.point[1][j1] = ExtremaPOnCurv::new(a_t11, a_p11);
        this.point[1][j2] = ExtremaPOnCurv::new(a_t22, a_p22);
        this.sq_dist[1] = a_p11.square_distance(&a_p22);

        this.point[2][j1] = ExtremaPOnCurv::new(a_t12, a_p12);
        this.point[2][j2] = ExtremaPOnCurv::new(a_t21, a_p21);
        this.sq_dist[2] = a_p12.square_distance(&a_p21);

        this.point[3][j1] = ExtremaPOnCurv::new(a_t12, a_p12);
        this.point[3][j2] = ExtremaPOnCurv::new(a_t22, a_p22);
        this.sq_dist[3] = a_p12.square_distance(&a_p22);

        // `cxx:1076-1111`: the intersection points, when the circles overlap.
        let b_out = a_d12 > (a_r1 + a_r2 + a_tol_d);
        let b_in = a_d12 < (a_r1 - a_r2 - a_tol_d);
        if !b_out && !b_in {
            let a_alpha = 0.5 * (a_r1 * a_r1 - a_r2 * a_r2 + a_d12 * a_d12) / a_d12;
            let mut a_val = a_r1 * a_r1 - a_alpha * a_alpha;
            if a_val < 0.0 {
                a_val = -a_val;
            }
            let a_beta = a_val.sqrt();
            let a_pt = a_pc1.xyz().added(&a_dir12.xyz().multiplied(a_alpha));

            // `aDLt = aDc1 ^ aDir12` uses the original `C1` axis direction.
            let a_dlt = a_dc1
                .crossed(&a_dir12)
                .expect("gp_Dir(aDc1 ^ aDir12): the circles are not coaxial");
            let a_pl1 = GpPnt::from_xyz(&a_pt.added(&a_dlt.xyz().multiplied(a_beta)));
            let a_pl2 = GpPnt::from_xyz(&a_pt.subtracted(&a_dlt.xyz().multiplied(a_beta)));

            let a_dist2 = a_pl1.square_distance(&a_pl2);
            let b_nb_ext6 = a_dist2 > a_tol_d2;

            let mut a_t = [0.0f64; 2];
            this.nb_ext = 5;
            a_t[j1] = clib::circle_parameter(&a_c1.position(), &a_pl1);
            a_t[j2] = clib::circle_parameter(&a_c2.position(), &a_pl1);
            this.point[4][j1] = ExtremaPOnCurv::new(a_t[j1], a_pl1);
            this.point[4][j2] = ExtremaPOnCurv::new(a_t[j2], a_pl1);
            this.sq_dist[4] = 0.0;

            if b_nb_ext6 {
                this.nb_ext = 6;
                a_t[j1] = clib::circle_parameter(&a_c1.position(), &a_pl2);
                a_t[j2] = clib::circle_parameter(&a_c2.position(), &a_pl2);
                this.point[5][j1] = ExtremaPOnCurv::new(a_t[j1], a_pl2);
                this.point[5][j2] = ExtremaPOnCurv::new(a_t[j2], a_pl2);
                this.sq_dist[5] = 0.0;
            }
        }

        this
    }

    /// `IsDone()` (`cxx:1117-1120`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `IsParallel()` (`cxx:1124-1131`): OCCT throws `StdFail_NotDone`.
    pub fn is_parallel(&self) -> Result<bool, String> {
        if !self.is_done() {
            return Err("StdFail_NotDone: Extrema_ExtElC::IsParallel".into());
        }
        Ok(self.is_par)
    }

    /// `NbExt()` (`cxx:1135-1142`): OCCT throws `StdFail_NotDone`.
    pub fn nb_ext(&self) -> Result<usize, String> {
        if !self.is_done() {
            return Err("StdFail_NotDone: Extrema_ExtElC::NbExt".into());
        }
        Ok(self.nb_ext)
    }

    /// `SquareDistance(N)` (`cxx:1146-1154`), 1-based: OCCT throws
    /// `Standard_OutOfRange` outside `1..=NbExt()`.
    pub fn square_distance(&self, n: usize) -> Result<f64, String> {
        let nb = self.nb_ext()?;
        if n < 1 || n > nb {
            return Err("Standard_OutOfRange: Extrema_ExtElC::SquareDistance".into());
        }
        Ok(self.sq_dist[n - 1])
    }

    /// `Points(N, P1, P2)` (`cxx:1158-1166`), 1-based, same range check.
    pub fn points(&self, n: usize) -> Result<(ExtremaPOnCurv, ExtremaPOnCurv), String> {
        let nb = self.nb_ext()?;
        if n < 1 || n > nb {
            return Err("Standard_OutOfRange: Extrema_ExtElC::Points".into());
        }
        Ok((self.point[n - 1][0], self.point[n - 1][1]))
    }
}

impl Default for ExtremaExtElC {
    fn default() -> Self {
        Self::new()
    }
}
