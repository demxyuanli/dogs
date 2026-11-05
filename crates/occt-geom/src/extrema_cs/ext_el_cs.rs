//! `Extrema_ExtElCS` — extrema between an elementary curve and an elementary
//! surface. Source: `Extrema_ExtElCS.cxx:45-806`, `.hxx:45-140`.
//!
//! | pair | OCCT | port |
//! |---|---|---|
//! | `Lin` / `Pln` | `cxx:55-69` | [`ExtElCS::perform_line_plane`] |
//! | `Lin` / `Cylinder` | `cxx:76-184` | [`ExtElCS::perform_line_cylinder`] |
//! | `Lin` / `Cone` | `cxx:191-195` | `Standard_NotImplemented` |
//! | `Lin` / `Sphere` | `cxx:203-273` | [`ExtElCS::perform_line_sphere`] |
//! | `Lin` / `Torus` | `cxx:280-284` | `Standard_NotImplemented` |
//! | `Circ` / `Pln` | `cxx:291-390` | [`ExtElCS::perform_circle_plane`] |
//! | `Circ` / `Cylinder` | `cxx:396-544` | [`ExtElCS::perform_circle_cylinder`] |
//! | `Circ` / `Cone` | `cxx:546-554` | `Standard_NotImplemented` |
//! | `Circ` / `Sphere` | `cxx:558-697` | [`ExtElCS::perform_circle_sphere`] |
//! | `Circ` / `Torus` | `cxx:699-706` | `Standard_NotImplemented` |
//! | `Hypr` / `Pln` | `cxx:708-764` | [`ExtElCS::perform_hyperbola_plane`] |
//!
//! `Circ` / `Cylinder` (`cxx:400-544`) starts from `Extrema_ExtElC anExtC(anAxis,
//! C, 0.)` (`cxx:410`), i.e. `Extrema_ExtElC(gp_Lin, gp_Circ)`
//! (`Extrema_ExtElC.cxx:471-625`). `Circ` / `Sphere` (`cxx:558-697`) uses
//! `Extrema_ExtElC(C, aCInt)` (`cxx:631`) — `Extrema_ExtElC(gp_Circ, gp_Circ)`
//! (`Extrema_ExtElC.cxx:955-1114`) — in the intersection-is-a-circle branch
//! (`cxx:628-695`). Both live in [`crate::extrema_el_c`].

use crate::extrema_pc::point_curve_extrema_all;
use crate::extrema_surf::{point_cylinder_extrema, point_sphere_extrema};
use crate::intana::{quadric_quadric_plane_sphere, IntAnaIntConicQuad, IntAnaQuadric, QuadricIntersection};
use crate::line::GeomLine;
use occt_core::gp::GpCirc;
use occt_core::precision::RESOLUTION;

use super::prelude::*;

/// `gp_Pln::SignedDistance(theLin)` (`gp_Pln.hxx:355-363`): zero unless the
/// line direction is parallel to the plane (`gp::Resolution()`), otherwise the
/// signed distance of the line's location.
fn pln_signed_distance_lin(pl: &GpPln, l: &GpLin) -> f64 {
    let dir = l.direction();
    if !pl.pos.direction().is_normal_tol(&dir, RESOLUTION) {
        return 0.0;
    }
    pl.signed_distance_pnt(&l.location())
}

/// `gp_Pln::SquareDistance(theP)` (`gp_Pln.hxx:159-163`).
fn pln_square_distance_pnt(pl: &GpPln, p: &GpPnt) -> f64 {
    let d = pl.signed_distance_pnt(p);
    d * d
}

/// `Extrema_ExtElCS` (`Extrema_ExtElCS.hxx:45-140`).
pub struct ExtElCS {
    /// `myDone` (`hxx:129`).
    done: bool,
    /// `myIsPar` (`hxx:130`).
    is_par: bool,
    /// `mySqDist` (`hxx:125`).
    sq_dist: Vec<f64>,
    /// `myPoint1` (`hxx:126`).
    point1: Vec<ExtremaPOnCurv>,
    /// `myPoint2` (`hxx:127`).
    point2: Vec<ExtremaPOnSurf>,
}

impl ExtElCS {
    /// `Extrema_ExtElCS()` (`cxx:45-50`).
    pub fn new() -> Self {
        ExtElCS {
            done: false,
            is_par: false,
            sq_dist: Vec::new(),
            point1: Vec::new(),
            point2: Vec::new(),
        }
    }

    /// Reinitialise the way every `Perform` overload starts.
    fn reset(&mut self, done: bool) {
        self.done = done;
        self.is_par = false;
        self.sq_dist.clear();
        self.point1.clear();
        self.point2.clear();
    }

    fn push(&mut self, sq_dist: f64, p1: ExtremaPOnCurv, p2: ExtremaPOnSurf) {
        self.sq_dist.push(sq_dist);
        self.point1.push(p1);
        self.point2.push(p2);
    }

    /// `Perform(gp_Lin, gp_Pln)` (`cxx:55-69`).
    pub fn perform_line_plane(&mut self, c: &GpLin, s: &GpPln) {
        self.reset(true);
        if c.direction().is_normal_tol(&s.axis().direction(), ANGULAR) {
            let d = pln_signed_distance_lin(s, c);
            self.sq_dist.push(d * d);
            self.is_par = true;
        }
    }

    /// `Perform(gp_Lin, gp_Cylinder)` (`cxx:76-184`).
    pub fn perform_line_cylinder(&mut self, c: &GpLin, s: &GpCylinder) {
        self.reset(false);
        let pos = s.position();
        let radius = s.radius();
        let axis = GpLin::new(*pos.axis());
        let extrem = crate::extrema_el_c::ExtremaExtElC::lines(&axis, c);

        let mut is_parallel = extrem.is_parallel().unwrap_or(false);
        if !is_parallel {
            // `cxx:97-99`.
            let (p_on_c1, p_on_c2) = extrem
                .points(1)
                .expect("Extrema_ExtElC::Points: not done");
            let pon_axis = *p_on_c1.value();
            let pc = *p_on_c2.value();

            // `cxx:102-133`: the line crosses the cylinder.
            if radius - pon_axis.distance(&pc) > PCONFUSION {
                let quadric = IntAnaQuadric::from_cylinder(s);
                let inter = IntAnaIntConicQuad::line_quadric(c, &quadric);
                if inter.is_done() && inter.is_in_quadric() {
                    is_parallel = true;
                } else if inter.is_done() {
                    let nb = inter.nb_points();
                    for i in 1..=nb {
                        let p_int = *inter.point(i);
                        let w = inter.param_on_conic(i);
                        let (u, v) = slib::cylinder_parameters(&pos, &p_int);
                        self.push(
                            0.0,
                            ExtremaPOnCurv::new(w, p_int),
                            ExtremaPOnSurf::new(u, v, p_int),
                        );
                    }
                }
            } else {
                // `cxx:135-152`: tangent or outside — the distances are the
                // ones from the perpendicular foot `PC` on the line.
                let exps = point_cylinder_extrema(s, &pc);
                if !exps.is_empty() {
                    for pair in &exps {
                        let q = pair.p2;
                        let (u, v) = slib::cylinder_parameters(&pos, &q);
                        self.push(
                            pc.square_distance(&q),
                            p_on_c2,
                            ExtremaPOnSurf::new(u, v, q),
                        );
                    }
                }
            }
            self.done = true;
        }

        if is_parallel {
            // `cxx:155-183`.
            let mut a_dist = extrem
                .square_distance(1)
                .expect("Extrema_ExtElC::SquareDistance: not done");
            let a_nb_ext = extrem.nb_ext().unwrap_or(1);
            for i in 2..=a_nb_ext {
                if let Ok(d) = extrem.square_distance(i) {
                    if d < a_dist {
                        a_dist = d;
                    }
                }
            }
            let d = a_dist.sqrt() - radius;
            self.sq_dist.push(d * d);
            self.done = true;
            self.is_par = true;
        }
    }

    /// `Perform(gp_Lin, gp_Cone)` (`cxx:191-195`): `Standard_NotImplemented`.
    pub fn perform_line_cone(&mut self, _c: &GpLin, _s: &GpCone) {
        panic!("Standard_NotImplemented: Extrema_ExtElCS::Perform(gp_Lin, gp_Cone)");
    }

    /// `Perform(gp_Lin, gp_Sphere)` (`cxx:203-273`).
    pub fn perform_line_sphere(&mut self, c: &GpLin, s: &GpSphere) {
        self.reset(false);
        let a_center = s.location();

        // `cxx:217-244`: `Extrema_ExtPElC(aCenter, C, Angular(),
        // RealFirst(), RealLast())` — the line arm of `Extrema_ExtPElC`
        // (`Extrema_ExtPElC.cxx:58-81`) is one solution, the perpendicular foot.
        let line = GeomLine::new(c.clone());
        let extrem = point_curve_extrema_all(&line, &a_center);

        if let Some(first) = extrem.first() {
            let my_p_on_c1 = ExtremaPOnCurv::new(first.u1, first.p2);
            // `cxx:222-246`: the line meets the sphere.
            if my_p_on_c1.value().distance(&a_center) <= s.radius() {
                let quadric = IntAnaQuadric::from_sphere(s);
                let inter = IntAnaIntConicQuad::line_quadric(c, &quadric);
                if inter.is_done() {
                    let nb = inter.nb_points();
                    for i in 1..=nb {
                        let p = *inter.point(i);
                        let w = inter.param_on_conic(i);
                        let (u, v) = slib::sphere_parameters(&s.position(), &p);
                        self.push(
                            0.0,
                            ExtremaPOnCurv::new(w, p),
                            ExtremaPOnSurf::new(u, v, p),
                        );
                    }
                }
            }

            // `cxx:249-270`: the extra perpendiculars from the foot.
            let exps = point_sphere_extrema(s, my_p_on_c1.value());
            for pair in &exps {
                let q = pair.p2;
                let v = pair.v1.expect("point_sphere_extrema: v1");
                self.push(
                    my_p_on_c1.value().square_distance(&q),
                    my_p_on_c1,
                    ExtremaPOnSurf::new(pair.u1, v, q),
                );
            }
        }
        // `cxx:271`: `myDone = true;` is unconditional.
        self.done = true;
    }

    /// `Perform(gp_Lin, gp_Torus)` (`cxx:280-284`).
    pub fn perform_line_torus(&mut self, _c: &GpLin, _s: &GpTorus) {
        panic!("Standard_NotImplemented: Extrema_ExtElCS::Perform(gp_Lin, gp_Torus)");
    }

    /// `Perform(gp_Circ, gp_Pln)` (`cxx:291-390`).
    pub fn perform_circle_plane(&mut self, c: &GpCirc, s: &GpPln) {
        self.reset(true);
        let pos = c.position();
        let n_circ = pos.direction();
        let pln_axis = s.axis();
        let n_pln = pln_axis.direction();

        let mut is_parallel = n_circ.is_parallel_tol(&n_pln, ANGULAR);
        let mut inter = IntAnaIntConicQuad::new();
        let mut t = [0.0f64; 2];
        if !is_parallel {
            // `cxx:312-325`.
            let ext_line = n_circ
                .crossed(&n_pln)
                .expect("standard_construction_error")
                .crossed(&n_circ)
                .expect("standard_construction_error");
            let x_dir = *pos.x_direction();
            t[0] = x_dir.angle_with_ref(&ext_line, &n_circ);
            if t[0] < 0.0 {
                t[0] += std::f64::consts::PI;
            }
            t[1] = t[0] + std::f64::consts::PI;
            // `cxx:326-341`. `myNbExt = 2 + NbPoints()` only sizes the arrays in
            // OCCT; the pushes below keep the same order.
            inter = IntAnaIntConicQuad::circle_plane(c, s, ANGULAR, CONFUSION);
            if inter.is_done() && inter.is_in_quadric() {
                is_parallel = true;
            }
        }

        if !is_parallel {
            // `cxx:343-381`: the two planes of the minor circle.
            for i in 0..2 {
                let pc = clib::circle_value(c, t[i]);
                let (u, v) = slib::plane_parameters(&s.position(), &pc);
                let pp = slib::plane_value(s, u, v);
                self.push(
                    pc.square_distance(&pp),
                    ExtremaPOnCurv::new(t[i], pc),
                    ExtremaPOnSurf::new(u, v, pp),
                );
            }
            // `cxx:366-380`: append the intersection points.
            if inter.is_done() && !inter.is_in_quadric() {
                let nb = inter.nb_points();
                if nb > 1 {
                    for i in 1..=nb {
                        let t_int = inter.param_on_conic(i);
                        let pc = clib::circle_value(c, t_int);
                        let (u, v) = slib::plane_parameters(&s.position(), &pc);
                        let pp = slib::plane_value(s, u, v);
                        self.push(
                            pc.square_distance(&pp),
                            ExtremaPOnCurv::new(t_int, pc),
                            ExtremaPOnSurf::new(u, v, pp),
                        );
                    }
                }
            }
        }

        if is_parallel {
            // `cxx:383-390`.
            self.sq_dist.push(pln_square_distance_pnt(s, &c.location()));
            self.is_par = true;
        }
    }

    /// `Perform(gp_Circ, gp_Cylinder)` (`cxx:400-544`).
    pub fn perform_circle_cylinder(&mut self, c: &GpCirc, s: &GpCylinder) {
        self.reset(false);

        // `cxx:407-415`.
        let an_axis = GpLin::new(s.axis());
        let an_ext_c = crate::extrema_el_c::ExtremaExtElC::line_circle(&an_axis, c);
        if !an_ext_c.is_done() {
            return;
        }

        let mut is_parallel = false;
        if an_ext_c.is_parallel().unwrap_or(false) {
            is_parallel = true;
        } else {
            let a_nb_ext = an_ext_c.nb_ext().unwrap_or(0);
            let a_cyl_rad = s.radius();
            let a_tol_conf = CONFUSION;

            // `cxx:429-438`: the intersection of the circle with the cylinder
            // also decides the parallel status.
            let a_cyl_quad = IntAnaQuadric::from_cylinder(s);
            let a_circ_cyl_inter = IntAnaIntConicQuad::circle_quadric(c, &a_cyl_quad);
            let mut a_nb_inter = 0;
            if a_circ_cyl_inter.is_done() && a_circ_cyl_inter.is_in_quadric() {
                is_parallel = true;
            } else if a_circ_cyl_inter.is_done() {
                a_nb_inter = a_circ_cyl_inter.nb_points();
            }

            if !is_parallel {
                // `cxx:446-483`: two shifted points along the axis/circle
                // direction for every line/circle extremum.
                for i in 1..=a_nb_ext {
                    let a_sq_dist = an_ext_c
                        .square_distance(i)
                        .expect("Extrema_ExtElC::SquareDistance: not done");
                    let a_dist = a_sq_dist.sqrt();
                    let (a_p_on_axis, a_p_on_circ) = an_ext_c
                        .points(i)
                        .expect("Extrema_ExtElC::Points: not done");

                    if a_sq_dist <= a_tol_conf * a_tol_conf {
                        continue;
                    }

                    let a_dir = GpDir::new(
                        a_p_on_axis.value().x() - a_p_on_circ.value().x(),
                        a_p_on_axis.value().y() - a_p_on_circ.value().y(),
                        a_p_on_axis.value().z() - a_p_on_circ.value().z(),
                    )
                    .expect("gp_Dir: the two extrema points are more than a confusion apart");
                    let a_shift = [a_dist + a_cyl_rad, a_dist - a_cyl_rad];

                    for shift in a_shift {
                        let a_pnt_on_cyl = GpPnt::from_xyz(
                            &a_p_on_circ.value().xyz().added(&a_dir.xyz().multiplied(shift)),
                        );
                        let (a_u, a_v) = slib::cylinder_parameters(&s.position(), &a_pnt_on_cyl);
                        self.push(
                            shift * shift,
                            a_p_on_circ,
                            ExtremaPOnSurf::new(a_u, a_v, a_pnt_on_cyl),
                        );
                    }
                }

                // `cxx:485-506`: the intersection points are zero-distance
                // solutions.
                for i in 1..=a_nb_inter {
                    let a_inter_pnt = *a_circ_cyl_inter.point(i);
                    let a_u = clib::circle_parameter(&c.position(), &a_inter_pnt);
                    let (a_uu, a_vv) = slib::cylinder_parameters(&s.position(), &a_inter_pnt);
                    self.push(
                        0.0,
                        ExtremaPOnCurv::new(a_u, a_inter_pnt),
                        ExtremaPOnSurf::new(a_uu, a_vv, a_inter_pnt),
                    );
                }
            }
        }

        // `cxx:509-511`.
        self.done = true;

        if is_parallel {
            // `cxx:513-541`.
            self.is_par = true;
            let mut a_dist = an_ext_c
                .square_distance(1)
                .expect("Extrema_ExtElC::SquareDistance: not done");
            let a_nb_ext = an_ext_c.nb_ext().unwrap_or(1);
            for i in 2..=a_nb_ext {
                let a_d = an_ext_c
                    .square_distance(i)
                    .expect("Extrema_ExtElC::SquareDistance: not done");
                if a_d < a_dist {
                    a_dist = a_d;
                }
            }
            let d = a_dist.sqrt() - s.radius();
            self.sq_dist.push(d * d);
        }
    }

    /// `Perform(gp_Circ, gp_Cone)` (`cxx:546-554`).
    pub fn perform_circle_cone(&mut self, _c: &GpCirc, _s: &GpCone) {
        panic!("Standard_NotImplemented: Extrema_ExtElCS::Perform(gp_Circ, gp_Cone)");
    }

    /// `Perform(gp_Circ, gp_Sphere)` (`cxx:558-697`).
    pub fn perform_circle_sphere(&mut self, c: &GpCirc, s: &GpSphere) {
        self.reset(false);

        // `cxx:572-587`: the circle's axis passes through the centre.
        let circ_axis = GpLin::new(c.axis());
        if circ_axis.square_distance(&s.location()) < SQUARE_CONFUSION {
            let a_sq_dist_loc = c.location().square_distance(&s.location());
            let a_sq_dist = a_sq_dist_loc + c.radius() * c.radius();
            let a_dist = a_sq_dist.sqrt() - s.radius();
            self.sq_dist.push(a_dist * a_dist);
            self.is_par = true;
            self.done = true;
            return;
        }

        // `cxx:589-594`: intersect the sphere with the circle's plane.
        let circ_ax1 = c.axis();
        let c_pln = GpPln::new(
            GpAx3::from_point_dir(c.location(), *circ_ax1.direction())
                .expect("gp_Ax3::from_point_dir"),
        );
        let an_inter = quadric_quadric_plane_sphere(&c_pln, s);
        if let QuadricIntersection::Circle(a_c_int) = an_inter {
            // `cxx:628-695`: the intersection is a circle. Both circles are
            // coplanar by construction, so `Extrema_ExtElC(C, aCInt)` resolves
            // the pairs of extremal points in that plane.
            let an_ext_c = crate::extrema_el_c::ExtremaExtElC::circles(c, &a_c_int);
            let a_nb_ext_c = an_ext_c.nb_ext().unwrap_or(0);
            let is_valid = an_ext_c.is_done()
                && !an_ext_c.is_parallel().unwrap_or(true)
                && a_nb_ext_c > 0;
            if !is_valid {
                // `cxx:639-643`: the parallel case was already handled above.
                return;
            }

            self.done = true;

            // `cxx:647-655`: the minimal distance.
            let mut a_min_sq_dist = f64::MAX;
            for i in 1..=a_nb_ext_c {
                let a_sq_dist = an_ext_c
                    .square_distance(i)
                    .expect("Extrema_ExtElC::SquareDistance: not done");
                if a_sq_dist < a_min_sq_dist {
                    a_min_sq_dist = a_sq_dist;
                }
            }

            // `cxx:657-666`: every solution within a square confusion.
            let mut a_sols = Vec::new();
            for i in 1..=a_nb_ext_c {
                let a_sq_dist = an_ext_c
                    .square_distance(i)
                    .expect("Extrema_ExtElC::SquareDistance: not done");
                if a_sq_dist - a_min_sq_dist < SQUARE_CONFUSION {
                    a_sols.push(i);
                }
            }

            // `cxx:668-695`: the distance is measured to the point on the
            // sphere rebuilt from the parameters.
            for idx in a_sols {
                let (p1, _p2) = an_ext_c.points(idx).expect("Extrema_ExtElC::Points: not done");
                let (a_u, a_v) = slib::sphere_parameters(&s.position(), p1.value());
                let a_p_on_s = slib::sphere_value(s, a_u, a_v);
                self.push(
                    p1.value().square_distance(&a_p_on_s),
                    p1,
                    ExtremaPOnSurf::new(a_u, a_v, a_p_on_s),
                );
            }
            return;
        }

        // `cxx:596-626`: empty, a point: one minimal solution.
        let a_t = clib::circle_parameter(&c.position(), &s.location());
        let a_p_on_c = clib::circle_value(c, a_t);
        let (a_u, a_v) = slib::sphere_parameters(&s.position(), &a_p_on_c);
        let a_p_on_s = slib::sphere_value(s, a_u, a_v);
        self.push(
            a_p_on_c.square_distance(&a_p_on_s),
            ExtremaPOnCurv::new(a_t, a_p_on_c),
            ExtremaPOnSurf::new(a_u, a_v, a_p_on_s),
        );
        self.done = true;
    }

    /// `Perform(gp_Circ, gp_Torus)` (`cxx:699-706`).
    pub fn perform_circle_torus(&mut self, _c: &GpCirc, _s: &GpTorus) {
        panic!("Standard_NotImplemented: Extrema_ExtElCS::Perform(gp_Circ, gp_Torus)");
    }

    /// `Perform(gp_Hypr, gp_Pln)` (`cxx:708-764`).
    pub fn perform_hyperbola_plane(&mut self, c: &GpHypr, s: &GpPln) {
        self.reset(true);
        let pos = c.position();
        let n_hypr = pos.direction();
        let pln_axis = s.axis();
        let n_pln = pln_axis.direction();

        if n_hypr.is_parallel_tol(&n_pln, ANGULAR) {
            self.sq_dist.push(pln_square_distance_pnt(s, &c.location()));
            self.is_par = true;
            return;
        }

        let x_dir = *pos.x_direction();
        let y_dir = *pos.y_direction();
        let a = c.minor_radius * n_pln.dot(&y_dir);
        let b = c.major_radius * n_pln.dot(&x_dir);
        if b.abs() > a.abs() {
            let t = -0.5 * ((a + b) / (b - a)).ln();
            let ph = clib::hyperbola_value(c, t);
            let sq_dist = pln_square_distance_pnt(s, &ph);
            let (u, v) = slib::plane_parameters(&s.position(), &ph);
            let pp = slib::plane_value(s, u, v);
            self.push(
                sq_dist,
                ExtremaPOnCurv::new(t, ph),
                ExtremaPOnSurf::new(u, v, pp),
            );
        }
    }

    /// `IsDone()` (`cxx:766-769`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `IsParallel()` (`cxx:799-806`): OCCT throws `StdFail_NotDone`.
    pub fn is_parallel(&self) -> bool {
        assert!(self.is_done(), "Extrema_ExtElCS::IsParallel: StdFail_NotDone");
        self.is_par
    }

    /// `NbExt()` (`cxx:771-777`).
    pub fn nb_ext(&self) -> usize {
        assert!(self.is_done(), "Extrema_ExtElCS::NbExt: StdFail_NotDone");
        self.sq_dist.len()
    }

    /// `SquareDistance(N)` (`cxx:779-787`), 1-based.
    pub fn square_distance(&self, n: usize) -> f64 {
        assert!(
            n >= 1 && n <= self.nb_ext(),
            "Extrema_ExtElCS::SquareDistance: Standard_OutOfRange"
        );
        self.sq_dist[n - 1]
    }

    /// `Points(N, P1, P2)` (`cxx:789-797`), 1-based.
    pub fn points(&self, n: usize) -> (ExtremaPOnCurv, ExtremaPOnSurf) {
        assert!(
            n >= 1 && n <= self.nb_ext(),
            "Extrema_ExtElCS::Points: Standard_OutOfRange"
        );
        (self.point1[n - 1], self.point2[n - 1])
    }
}

impl Default for ExtElCS {
    fn default() -> Self {
        ExtElCS::new()
    }
}
