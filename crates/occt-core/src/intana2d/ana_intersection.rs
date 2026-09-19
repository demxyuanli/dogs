//! 2D analytical intersection of `gp` elements and conics.
//! Source: `IntAna2d_AnaIntersection.hxx/.lxx/.cxx` and the eight
//! `IntAna2d_AnaIntersection_<n>.cxx` `Perform` overloads.
use std::f64::consts::PI;

use crate::gp::{GpAx2d, GpCirc2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpVec2d};

use super::conic::IntAna2dConic;
use super::int_point::IntAna2dIntPoint;
use super::outils::{
    coord_ancien_repere, traitement_points_confondus, MyDirectPolynomialRoots,
};

#[derive(Debug, Clone, Copy)]
pub struct IntAna2dAnaIntersection {
    done: bool,
    para: bool,
    iden: bool,
    empt: bool,
    nbp: i32,
    lpnt: [IntAna2dIntPoint; 4],
}

/// `RealEpsilon()` (`Standard_Real.hxx:161-164`).
const REAL_EPSILON: f64 = f64::EPSILON;

impl IntAna2dAnaIntersection {
    /// Empty ctor (`cxx:22-31`).
    pub fn new() -> Self {
        Self {
            done: false,
            para: false,
            iden: false,
            empt: true,
            nbp: 0,
            lpnt: [IntAna2dIntPoint::new(); 4],
        }
    }

    /// `IsDone` (`lxx:17-20`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `IsEmpty` (`lxx:22-31`).
    pub fn is_empty(&self) -> bool {
        assert!(self.done, "StdFail_NotDone in IntAna2d_AnaIntersection");
        self.nbp == 0 && !self.iden
    }

    /// `IdenticalElements` (`lxx:33-41`).
    pub fn identical_elements(&self) -> bool {
        assert!(self.done, "StdFail_NotDone in IntAna2d_AnaIntersection");
        self.iden
    }

    /// `ParallelElements` (`lxx:43-51`).
    pub fn parallel_elements(&self) -> bool {
        assert!(self.done, "StdFail_NotDone in IntAna2d_AnaIntersection");
        self.para
    }

    /// `NbPoints` (`lxx:53-61`).
    pub fn nb_points(&self) -> i32 {
        assert!(self.done, "StdFail_NotDone in IntAna2d_AnaIntersection");
        self.nbp
    }

    /// `Point(N)` (`lxx:63-83`), 1-based.
    pub fn point(&self, n: i32) -> &IntAna2dIntPoint {
        assert!(self.done, "StdFail_NotDone in IntAna2d_AnaIntersection");
        assert!(
            n > 0 && n <= self.nbp,
            "Standard_OutOfRange in IntAna2d_AnaIntersection::Point"
        );
        &self.lpnt[(n - 1) as usize]
    }

    /// `Perform(const gp_Lin2d&, const gp_Lin2d&)` (`_1.cxx:23-129`).
    pub fn perform_lin_lin(&mut self, l1: &GpLin2d, l2: &GpLin2d) {
        self.done = false;

        let (a1, b1, c1) = lin_coefficients(l1);
        let (a2, b2, c2) = lin_coefficients(l2);

        let det = a1
            .abs()
            .max(a2.abs().max(b1.abs().max(b2.abs())));

        let (al1, be1, ga1, al2, be2, ga2);
        if a1.abs() == det {
            al1 = a1;
            be1 = b1;
            ga1 = c1;
            al2 = a2;
            be2 = b2;
            ga2 = c2;
        } else if b1.abs() == det {
            al1 = b1;
            be1 = a1;
            ga1 = c1;
            al2 = b2;
            be2 = a2;
            ga2 = c2;
        } else if a2.abs() == det {
            al1 = a2;
            be1 = b2;
            ga1 = c2;
            al2 = a1;
            be2 = b1;
            ga2 = c1;
        } else {
            al1 = b2;
            be1 = a2;
            ga1 = c2;
            al2 = b1;
            be2 = a1;
            ga2 = c1;
        }

        let rap = al2 / al1;
        let denom = be2 - rap * be1;

        if denom.abs() <= REAL_EPSILON {
            // Coincident directions.
            self.para = true;
            self.nbp = 0;
            if (ga2 - rap * ga1).abs() <= REAL_EPSILON {
                // Coincident lines.
                self.iden = true;
                self.empt = false;
            } else {
                // Parallel lines.
                self.iden = false;
                self.empt = true;
            }
        } else {
            self.para = false;
            self.iden = false;
            self.empt = false;
            self.nbp = 1;
            let mut xs = (be1 * ga2 / al1 - be2 * ga1 / al1) / denom;
            let mut ys = (rap * ga1 - ga2) / denom;

            if ((a1.abs() != det) && (b1.abs() == det))
                || ((a1.abs() != det) && (b1.abs() != det) && (a2.abs() != det))
            {
                std::mem::swap(&mut xs, &mut ys);
            }

            let la;
            let mu;
            if a1.abs() >= b1.abs() {
                la = (ys - l1.location().y()) / a1;
            } else {
                la = (l1.location().x() - xs) / b1;
            }
            if a2.abs() >= b2.abs() {
                mu = (ys - l2.location().y()) / a2;
            } else {
                mu = (l2.location().x() - xs) / b2;
            }
            self.lpnt[0].set_value4(xs, ys, la, mu);
        }
        self.done = true;
    }

    /// `Perform(const gp_Circ2d&, const gp_Circ2d&)` (`_2.cxx:23-211`).
    pub fn perform_circ_circ(&mut self, c1: &GpCirc2d, c2: &GpCirc2d) {
        self.done = false;
        let d = c1.location().distance(&c2.location());
        let r1 = c1.radius;
        let r2 = c2.radius;
        let sum = r1 + r2;
        let dif = (r1 - r2).abs();

        if d <= REAL_EPSILON {
            // Concentric circles.
            self.para = true;
            self.nbp = 0;
            if dif <= REAL_EPSILON {
                // Coincident circles.
                self.empt = false;
                self.iden = true;
            } else {
                // Parallel circles.
                self.empt = true;
                self.iden = false;
            }
        } else if (d - sum) > crate::precision::epsilon(sum) {
            // Circles exterior to each other, no solution.
            self.empt = true;
            self.para = false;
            self.iden = false;
            self.nbp = 0;
        } else if (d - sum).abs() <= crate::precision::epsilon(sum) {
            // Externally tangent circles.
            self.empt = false;
            self.para = false;
            self.iden = false;
            self.nbp = 1;
            let ax = vec_between(&c1.location(), &c2.location());
            let ox1 = GpVec2d::from_dir2d(&c1.x_axis().vdir);
            let ox2 = GpVec2d::from_dir2d(&c2.x_axis().vdir);

            let xs = (c1.location().x() * r2 + c2.location().x() * r1) / sum;
            let ys = (c1.location().y() * r2 + c2.location().y() * r1) / sum;
            let mut ang1 = ox1.angle(&ax); // Result between -PI and +PI
            let ang2 = ox2.angle(&ax) + PI;
            if ang1 < 0.0 {
                ang1 = 2.0 * PI + ang1; // Normalize to range [0, 2PI]
            }
            self.lpnt[0].set_value4(xs, ys, ang1, ang2);
        } else if ((sum - d) > crate::precision::epsilon(sum))
            && ((d - dif) > crate::precision::epsilon(d + dif))
        {
            self.empt = false;
            self.para = false;
            self.iden = false;
            self.nbp = 2;
            let ax = vec_between(&c1.location(), &c2.location());
            let ox1 = GpVec2d::from_dir2d(&c1.x_axis().vdir);
            let ox2 = GpVec2d::from_dir2d(&c2.x_axis().vdir);
            let ref1 = ox1.angle(&ax); // Result between -PI and +PI
            let ref2 = ox2.angle(&ax); // Result between -PI and +PI

            let mut l1 = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
            if r1 * r1 - l1 * l1 < 0.0 {
                l1 = if l1 > 0.0 { r1 } else { -r1 };
            }
            let h = (r1 * r1 - l1 * l1).sqrt();

            let xs1 = c1.location().x() + l1 * ax.x() / d - h * ax.y() / d;
            let ys1 = c1.location().y() + l1 * ax.y() / d + h * ax.x() / d;

            let xs2 = c1.location().x() + l1 * ax.x() / d + h * ax.y() / d;
            let ys2 = c1.location().y() + l1 * ax.y() / d - h * ax.x() / d;

            let sint1 = h / r1;
            let cost1 = l1 / r1;

            let sint2 = h / r2;
            let cost2 = (l1 - d) / r2;

            // ang1 and ang2 correspond to the solutions with positive sine when
            // the reference axis is the axis of centers C1C2. Arccos is used
            // between pi/2 and 3pi/2, arcsin otherwise.
            let ang1 = if cost1.abs() <= 0.707 {
                cost1.acos()
            } else {
                let mut a = sint1.asin();
                if cost1 < 0.0 {
                    a = PI - a;
                }
                a
            };
            let ang2 = if cost2.abs() <= 0.707 {
                cost2.acos()
            } else {
                let mut a = sint2.asin();
                if cost2 < 0.0 {
                    a = PI - a;
                }
                a
            };
            let mut ang11 = ref1 + ang1;
            let mut ang21 = ref2 + ang2;
            let mut ang12 = ref1 - ang1;
            let mut ang22 = ref2 - ang2;
            if ang11 < 0.0 {
                ang11 = 2.0 * PI + ang11;
            } else if ang11 >= 2.0 * PI {
                ang11 -= 2.0 * PI;
            }
            if ang21 < 0.0 {
                ang21 = 2.0 * PI + ang21;
            } else if ang21 >= 2.0 * PI {
                ang21 -= 2.0 * PI;
            }
            if ang12 < 0.0 {
                ang12 = 2.0 * PI + ang12;
            } else if ang12 >= 2.0 * PI {
                ang12 -= 2.0 * PI;
            }
            if ang22 < 0.0 {
                ang22 = 2.0 * PI + ang22;
            } else if ang22 >= 2.0 * PI {
                ang22 -= 2.0 * PI;
            }
            self.lpnt[0].set_value4(xs1, ys1, ang11, ang21);
            self.lpnt[1].set_value4(xs2, ys2, ang12, ang22);
        } else if (d - dif).abs() <= crate::precision::epsilon(sum) {
            // Internally tangent circles.
            self.empt = false;
            self.para = false;
            self.iden = false;
            self.nbp = 1;
            let mut ax = vec_between(&c1.location(), &c2.location());
            if c1.radius < c2.radius {
                ax.reverse();
            }

            let ox1 = GpVec2d::from_dir2d(&c1.x_axis().vdir);
            let ox2 = GpVec2d::from_dir2d(&c2.x_axis().vdir);
            let mut ang1 = ox1.angle(&ax); // Result between -PI and +PI
            let mut ang2 = ox2.angle(&ax);
            if ang1 < 0.0 {
                ang1 = 2.0 * PI + ang1; // Normalize to range [0, 2PI]
            }
            if ang2 < 0.0 {
                ang2 = 2.0 * PI + ang2; // Normalize to range [0, 2PI]
            }
            let xs = (c1.location().x() * r2 - c2.location().x() * r1) / (r2 - r1);
            let ys = (c1.location().y() * r2 - c2.location().y() * r1) / (r2 - r1);
            self.lpnt[0].set_value4(xs, ys, ang1, ang2);
        } else {
            // d < dif - Resolution and d != 0: one circle inside the other.
            self.empt = true;
            self.para = false;
            self.iden = false;
            self.nbp = 0;
        }
        self.done = true;
    }

    /// `Perform(const gp_Lin2d&, const gp_Circ2d&)` (`_3.cxx:27-107`).
    pub fn perform_lin_circ(&mut self, l: &GpLin2d, c: &GpCirc2d) {
        self.done = false;
        self.iden = false;
        self.para = false;

        let (a, b, c0) = lin_coefficients(l);
        let d = a * c.location().x() + b * c.location().y() + c0;

        if d.abs() - c.radius > crate::precision::epsilon(c.radius) {
            self.empt = true;
            self.nbp = 0;
        } else {
            // At least 1 solution.
            self.empt = false;
            if (d.abs() - c.radius).abs() <= crate::precision::epsilon(c.radius) {
                // Tangency case.
                self.nbp = 1;
                let xs = c.location().x() - d * a;
                let ys = c.location().y() - d * b;
                let p = crate::gp::GpPnt2d::new(xs, ys);
                let u = crate::elib::clib::line2d_parameter(l, &p);
                let ang = crate::elib::clib::circle2d_parameter(c.position(), &p);
                self.lpnt[0].set_value4(xs, ys, u, ang);
            } else {
                // 2 intersection points.
                self.nbp = 2;
                let h = (c.radius * c.radius - d * d).sqrt();
                let xs1 = c.location().x() - d * a - h * b;
                let ys1 = c.location().y() - d * b + h * a;
                let xs2 = c.location().x() - d * a + h * b;
                let ys2 = c.location().y() - d * b - h * a;

                let p1 = crate::gp::GpPnt2d::new(xs1, ys1);
                let p2 = crate::gp::GpPnt2d::new(xs2, ys2);
                let u1 = crate::elib::clib::line2d_parameter(l, &p1);
                let u2 = crate::elib::clib::line2d_parameter(l, &p2);
                let ang1 = crate::elib::clib::circle2d_parameter(c.position(), &p1);
                let ang2 = crate::elib::clib::circle2d_parameter(c.position(), &p2);

                self.lpnt[0].set_value4(xs1, ys1, u1, ang1);
                self.lpnt[1].set_value4(xs2, ys2, u2, ang2);
            }
        }
        self.done = true;
    }

    /// `Perform(const gp_Lin2d&, const IntAna2d_Conic&)` (`_4.cxx:26-76`).
    pub fn perform_lin_conic(&mut self, l: &GpLin2d, conic: &IntAna2dConic) {
        self.done = false;
        self.nbp = 0;
        self.para = false;
        self.iden = false;

        let (a, b, c, d, e, f) = conic.coefficients();
        let (dr_a, dr_b, _dr_c) = lin_coefficients(l);
        let x0 = l.location().x();
        let y0 = l.location().y();

        // Parameter L: X = Xo - L*DR_B and Y = Yo + L*DR_A
        let px0 = f + x0 * (d + d + a * x0 + 2.0 * c * y0) + y0 * (e + e + b * y0);
        let px1 = 2.0
            * (e * dr_a - d * dr_b + x0 * (c * dr_a - a * dr_b) + y0 * (b * dr_a - c * dr_b));
        let px2 = dr_a * (b * dr_a - 2.0 * c * dr_b) + a * (dr_b * dr_b);

        let sol = MyDirectPolynomialRoots::new3(px2, px1, px0);

        if !sol.is_done() {
            self.done = false;
            return;
        }
        if sol.infinite_roots() {
            self.iden = true;
            self.done = true;
            return;
        }
        self.nbp = sol.nb_solutions();
        for i in 1..=self.nbp {
            let s = sol.value(i);
            let tx = x0 - s * dr_b;
            let ty = y0 + s * dr_a;
            self.lpnt[(i - 1) as usize].set_value3(tx, ty, s);
        }
        traitement_points_confondus(&mut self.nbp, &mut self.lpnt);
        self.done = true;
    }

    /// `Perform(const gp_Circ2d&, const IntAna2d_Conic&)` (`_5.cxx:26-86`).
    pub fn perform_circ_conic(&mut self, circle: &GpCirc2d, conic: &IntAna2dConic) {
        let c_is_direct = circle.position().vxdir.crossed(&circle.position().vydir) >= 0.0;
        let radius = circle.radius;
        let radius_p2 = radius * radius;

        self.done = false;
        self.nbp = 0;
        self.para = false;
        self.empt = false;
        self.iden = false;

        let axe_rep = GpAx2d::new(circle.location(), circle.position().vxdir);

        let (a, b, c, d, e, f) = conic.new_coefficients(&axe_rep);

        // Parameter a with x = Radius*cos(a) and y = Radius*sin(a).
        let pss = b * radius_p2;
        let pcc = a * radius_p2 - pss; // COS ^2
        let p2sc = c * radius_p2; // 2 SIN COS
        let pc = 2.0 * d * radius; // COS
        let ps = 2.0 * e * radius; // SIN
        let pcte = f + pss; // 1

        let sol = crate::math_trig_roots::TrigonometricFunctionRoots::new_abcde(
            pcc,
            p2sc,
            pc,
            ps,
            pcte,
            0.0,
            2.0 * PI,
        );

        if !sol.is_done() {
            self.done = false;
            return;
        }
        if sol.infinite_roots() {
            self.iden = true;
            self.done = true;
            return;
        }
        self.nbp = sol.nb_solutions();
        for i in 1..=self.nbp {
            let mut s = sol.value(i);
            let mut tx = radius * s.cos();
            let mut ty = radius * s.sin();
            coord_ancien_repere(&mut tx, &mut ty, &axe_rep);
            if !c_is_direct {
                s = PI + PI - s;
            }
            self.lpnt[(i - 1) as usize].set_value3(tx, ty, s);
        }
        traitement_points_confondus(&mut self.nbp, &mut self.lpnt);
        self.done = true;
    }

    /// `Perform(const gp_Elips2d&, const IntAna2d_Conic&)` (`_6.cxx:28-84`).
    pub fn perform_elips_conic(&mut self, elips: &GpElips2d, conic: &IntAna2dConic) {
        let e_is_direct = elips.pos.vxdir.crossed(&elips.pos.vydir) >= 0.0;
        let minor_radius = elips.minor_radius;
        let major_radius = elips.major_radius;

        self.done = false;
        self.nbp = 0;
        self.para = false;
        self.iden = false;
        self.empt = false;

        let axe_rep = GpAx2d::new(elips.pos.point, elips.pos.vxdir);

        let (a, b, c, d, e, f) = conic.new_coefficients(&axe_rep);

        // Parameter: a with x = MajorRadius*cos(a) and y = MinorRadius*sin(a).
        let pss = b * minor_radius * minor_radius; // SIN ^2
        let pcc = a * major_radius * major_radius - pss; // COS ^2
        let p2sc = c * major_radius * minor_radius; // 2 SIN COS
        let pc = 2.0 * d * major_radius; // COS
        let ps = 2.0 * e * minor_radius; // SIN
        let pcte = f + pss; // 1

        let sol = crate::math_trig_roots::TrigonometricFunctionRoots::new_abcde(
            pcc,
            p2sc,
            pc,
            ps,
            pcte,
            0.0,
            2.0 * PI,
        );

        if !sol.is_done() {
            self.done = false;
            return;
        }
        if sol.infinite_roots() {
            self.iden = true;
            self.done = true;
            return;
        }
        self.nbp = sol.nb_solutions();
        for i in 1..=self.nbp {
            let mut s = sol.value(i);
            let mut tx = major_radius * s.cos();
            let mut ty = minor_radius * s.sin();
            coord_ancien_repere(&mut tx, &mut ty, &axe_rep);
            if !e_is_direct {
                s = PI + PI - s;
            }
            self.lpnt[(i - 1) as usize].set_value3(tx, ty, s);
        }
        traitement_points_confondus(&mut self.nbp, &mut self.lpnt);
        self.done = true;
    }

    /// `Perform(const gp_Parab2d&, const IntAna2d_Conic&)` (`_7.cxx:27-82`).
    ///
    /// The OCCT body does not return after `!Sol.IsDone()`, so `done` ends up
    /// `true` on every path; that quirk is reproduced here.
    pub fn perform_parab_conic(&mut self, p: &GpParab2d, conic: &IntAna2dConic) {
        let p_is_direct = p.pos.vxdir.crossed(&p.pos.vydir) >= 0.0;
        let un_sur_2p = 0.5 / (2.0 * p.focal);

        self.done = false;
        self.nbp = 0;
        self.para = false;
        self.empt = false;
        self.iden = false;

        let axe_rep = GpAx2d::new(p.pos.point, p.pos.vxdir);

        let (a, b, c, d, e, f) = conic.new_coefficients(&axe_rep);

        // Parameter y: x = y^2 / (2p).
        let px0 = f;
        let px1 = e + e;
        let px2 = b + un_sur_2p * (d + d);
        let px3 = (c + c) * un_sur_2p;
        let px4 = a * (un_sur_2p * un_sur_2p);

        let sol = MyDirectPolynomialRoots::new5(px4, px3, px2, px1, px0);

        if !sol.is_done() {
            self.done = false;
        } else {
            if sol.infinite_roots() {
                self.iden = true;
                self.done = true;
            }
            self.nbp = sol.nb_solutions();
            for i in 1..=self.nbp {
                let mut s = sol.value(i);
                let mut tx = un_sur_2p * s * s;
                let mut ty = s;
                coord_ancien_repere(&mut tx, &mut ty, &axe_rep);
                if !p_is_direct {
                    s = -s;
                }
                self.lpnt[(i - 1) as usize].set_value3(tx, ty, s);
            }
            traitement_points_confondus(&mut self.nbp, &mut self.lpnt);
        }
        self.done = true;
    }

    /// `Perform(const gp_Hypr2d&, const IntAna2d_Conic&)` (`_8.cxx:42-123`).
    pub fn perform_hypr_conic(&mut self, h: &GpHypr2d, conic: &IntAna2dConic) {
        let h_is_direct = h.pos.vxdir.crossed(&h.pos.vydir) >= 0.0;
        let minor_radius = h.minor_radius;
        let major_radius = h.major_radius;

        self.done = false;
        self.nbp = 0;
        self.para = false;
        self.iden = false;
        self.empt = false;

        let axe_rep = GpAx2d::new(h.pos.point, h.pos.vxdir);
        let (a, b, c, d, e, f) = conic.new_coefficients(&axe_rep);

        let a_major_radius_p2 = a * major_radius * major_radius;
        let b_minor_radius_p2 = b * minor_radius * minor_radius;
        let c_2_major_minor_radius = c * 2.0 * major_radius * minor_radius;

        // Parameter t with x = MajorRadius*Ch(t), y = MinorRadius*Sh(t). The
        // polynomial is rewritten in Exp(t); the coefficients below are those
        // of P multiplied by 4*Exp(t)^2.
        let px0 = a_major_radius_p2 - c_2_major_minor_radius + b_minor_radius_p2;
        let px1 = 4.0 * (d * major_radius - e * minor_radius);
        let px2 = 2.0 * (a_major_radius_p2 + 2.0 * f - b_minor_radius_p2);
        let px3 = 4.0 * (d * major_radius + e * minor_radius);
        let px4 = a_major_radius_p2 + c_2_major_minor_radius + b_minor_radius_p2;

        let sol = MyDirectPolynomialRoots::new5(px4, px3, px2, px1, px0);

        if !sol.is_done() {
            self.done = false;
            return;
        }
        if sol.infinite_roots() {
            self.iden = true;
            self.done = true;
            return;
        }
        // We have X = (CosH(t)*major_radius)/2, Y = (SinH(t)*minor_radius)/2;
        // the resolution is in S = Exp(t).
        self.nbp = sol.nb_solutions();
        let mut nb_sol_valides = 0;
        for i in 1..=self.nbp {
            let mut s = sol.value(i);
            if s > REAL_EPSILON {
                let mut tx = 0.5 * major_radius * (s + 1.0 / s);
                let mut ty = 0.5 * minor_radius * (s - 1.0 / s);

                //--- Are we on the correct branch of the hyperbola? We assume
                //--- the deviation on curve1 is zero (the point comes from the
                //--- parameterisation).
                nb_sol_valides += 1;
                coord_ancien_repere(&mut tx, &mut ty, &axe_rep);
                s = s.ln();
                if !h_is_direct {
                    s = -s;
                }
                self.lpnt[(nb_sol_valides - 1) as usize].set_value3(tx, ty, s);
            }
        }
        self.nbp = nb_sol_valides;
        traitement_points_confondus(&mut self.nbp, &mut self.lpnt);
        self.done = true;
    }
}

impl Default for IntAna2dAnaIntersection {
    fn default() -> Self {
        Self::new()
    }
}

/// `gp_Lin2d::Coefficients` (`gp_Lin2d.hxx:91-96`).
fn lin_coefficients(l: &GpLin2d) -> (f64, f64, f64) {
    let a = l.pos.vdir.y;
    let b = -l.pos.vdir.x;
    let c = -(a * l.pos.loc.x() + b * l.pos.loc.y());
    (a, b, c)
}

/// `gp_Vec2d(const gp_Pnt2d&, const gp_Pnt2d&)` (`gp_Vec2d.hxx`).
fn vec_between(p1: &crate::gp::GpPnt2d, p2: &crate::gp::GpPnt2d) -> GpVec2d {
    GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y())
}
