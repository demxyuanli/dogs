//! Triangle/triangle probing: `maxSR` / `minSR` / `project6` / `TriContact`
//! (`IntPolyh_MaillageAffinage.cxx:1389-1669`), `TestNbPoints`
//! (`cxx:1677-1757`), `StartingPointsResearch` (`cxx:1767-1897`),
//! `NextStartingPointsResearch` (`cxx:1902-2048`) and `TriangleCompare`
//! (`cxx:3141-3185`).

use super::super::edge_couple::IntPolyhCouple;
use super::super::point::IntPolyhPoint;
use super::super::start_point::IntPolyhStartPoint;
use super::{MaillageAffinage, SQUARE_MY_CONFUSION_PRECISION};

/// `maxSR(a, b, c)` (`cxx:1393-1406`).
pub(super) fn max_sr(a: f64, b: f64, c: f64) -> f64 {
    let mut t = a;
    if b > t {
        t = b;
    }
    if c > t {
        t = c;
    }
    t
}

/// `minSR(a, b, c)` (`cxx:1413-1425`).
pub(super) fn min_sr(a: f64, b: f64, c: f64) -> f64 {
    let mut t = a;
    if b < t {
        t = b;
    }
    if c < t {
        t = c;
    }
    t
}

/// `project6(ax, p1, p2, p3, q1, q2, q3)` (`cxx:1429-1457`).
pub(super) fn project6(
    ax: &IntPolyhPoint,
    p1: &IntPolyhPoint,
    p2: &IntPolyhPoint,
    p3: &IntPolyhPoint,
    q1: &IntPolyhPoint,
    q2: &IntPolyhPoint,
    q3: &IntPolyhPoint,
) -> bool {
    let pp1 = ax.dot(p1);
    let pp2 = ax.dot(p2);
    let pp3 = ax.dot(p3);
    let qq1 = ax.dot(q1);
    let qq2 = ax.dot(q2);
    let qq3 = ax.dot(q3);

    let mx1 = max_sr(pp1, pp2, pp3);
    let mn1 = min_sr(pp1, pp2, pp3);
    let mx2 = max_sr(qq1, qq2, qq3);
    let mn2 = min_sr(qq1, qq2, qq3);

    if mn1 > mx2 {
        return false;
    }
    if mn2 > mx1 {
        return false;
    }
    true
}

/// `TriContact(P1, P2, P3, Q1, Q2, Q3, Angle)` (`cxx:1466-1669`): `(contact,
/// Angle)`. `Angle` keeps its incoming value when a normal is too small, as in
/// the `.cxx` where it is only written inside the guard.
pub fn tri_contact(
    p1: &IntPolyhPoint,
    p2: &IntPolyhPoint,
    p3: &IntPolyhPoint,
    q1: &IntPolyhPoint,
    q2: &IntPolyhPoint,
    q3: &IntPolyhPoint,
    angle: &mut f64,
) -> bool {
    if max_sr(p1.x, p2.x, p3.x) < min_sr(q1.x, q2.x, q3.x) {
        return false;
    }
    if max_sr(p1.y, p2.y, p3.y) < min_sr(q1.y, q2.y, q3.y) {
        return false;
    }
    if max_sr(p1.z, p2.z, p3.z) < min_sr(q1.z, q2.z, q3.z) {
        return false;
    }

    if min_sr(p1.x, p2.x, p3.x) > max_sr(q1.x, q2.x, q3.x) {
        return false;
    }
    if min_sr(p1.y, p2.y, p3.y) > max_sr(q1.y, q2.y, q3.y) {
        return false;
    }
    if min_sr(p1.z, p2.z, p3.z) > max_sr(q1.z, q2.z, q3.z) {
        return false;
    }

    // `cxx:1521-1547`: the triangles are translated so that P1 is the origin.
    let pp1 = IntPolyhPoint::new(0.0, 0.0, 0.0, 0.0, 0.0);
    let pp2 = sub_xyz(p2, p1);
    let pp3 = sub_xyz(p3, p1);
    let qq1 = sub_xyz(q1, p1);
    let qq2 = sub_xyz(q2, p1);
    let qq3 = sub_xyz(q3, p1);

    let e1 = sub_xyz(&pp2, &pp1);
    let e2 = sub_xyz(&pp3, &pp2);
    let e3 = sub_xyz(&pp1, &pp3);

    let f1 = sub_xyz(&qq2, &qq1);
    let f2 = sub_xyz(&qq3, &qq2);
    let f3 = sub_xyz(&qq1, &qq3);

    let mut n1 = IntPolyhPoint::default();
    n1.cross(&e1, &e2); // normal to the first triangle
    let mut m1 = IntPolyhPoint::default();
    m1.cross(&f1, &f2); // normal to the second triangle

    let mut g1 = IntPolyhPoint::default();
    let mut g2 = IntPolyhPoint::default();
    let mut g3 = IntPolyhPoint::default();
    let mut h1 = IntPolyhPoint::default();
    let mut h2 = IntPolyhPoint::default();
    let mut h3 = IntPolyhPoint::default();
    g1.cross(&e1, &n1);
    g2.cross(&e2, &n1);
    g3.cross(&e3, &n1);
    h1.cross(&f1, &m1);
    h2.cross(&f2, &m1);
    h3.cross(&f3, &m1);

    let mut ef11 = IntPolyhPoint::default();
    let mut ef12 = IntPolyhPoint::default();
    let mut ef13 = IntPolyhPoint::default();
    let mut ef21 = IntPolyhPoint::default();
    let mut ef22 = IntPolyhPoint::default();
    let mut ef23 = IntPolyhPoint::default();
    let mut ef31 = IntPolyhPoint::default();
    let mut ef32 = IntPolyhPoint::default();
    let mut ef33 = IntPolyhPoint::default();
    ef11.cross(&e1, &f1);
    ef12.cross(&e1, &f2);
    ef13.cross(&e1, &f3);
    ef21.cross(&e2, &f1);
    ef22.cross(&e2, &f2);
    ef23.cross(&e2, &f3);
    ef31.cross(&e3, &f1);
    ef32.cross(&e3, &f2);
    ef33.cross(&e3, &f3);

    // `cxx:1585-1657`: the separating-axis tests.
    for ax in [
        &n1, &m1, &ef11, &ef12, &ef13, &ef21, &ef22, &ef23, &ef31, &ef32, &ef33, &g1, &g2, &g3,
        &h1, &h2, &h3,
    ] {
        if !project6(ax, &pp1, &pp2, &pp3, &qq1, &qq2, &qq3) {
            return false;
        }
    }

    // `cxx:1658-1667`: cosinus of the angle between the normals.
    let mut sq_modn1 = -1.0;
    let mut sq_modm1 = -1.0;
    sq_modn1 = n1.square_modulus();
    if sq_modn1 > SQUARE_MY_CONFUSION_PRECISION {
        sq_modm1 = m1.square_modulus();
    }
    if sq_modm1 > SQUARE_MY_CONFUSION_PRECISION {
        *angle = n1.dot(&m1) / (sq_modn1.sqrt() * sq_modm1.sqrt());
    }
    true
}

/// `P - p1` with the UV slots left at zero, matching the `.cxx` translation
/// which only writes `X`/`Y`/`Z`.
fn sub_xyz(p: &IntPolyhPoint, base: &IntPolyhPoint) -> IntPolyhPoint {
    IntPolyhPoint::new(p.x - base.x, p.y - base.y, p.z - base.z, 0.0, 0.0)
}

/// `TestNbPoints(TriSurfID, NbPoints, NbPointsTotal, Pt1, Pt2, SP1, SP2)`
/// (`cxx:1677-1757`).
pub(super) fn test_nb_points(
    _tri_surf_id: i32,
    nb_points: i32,
    nb_points_total: &mut i32,
    pt1: &IntPolyhStartPoint,
    pt2: &IntPolyhStartPoint,
    sp1: &mut IntPolyhStartPoint,
    sp2: &mut IntPolyhStartPoint,
) {
    if nb_points > 2 {
        return;
    }

    if nb_points == 1 && *nb_points_total == 0 {
        *sp1 = *pt1;
        *nb_points_total = 1;
    } else if nb_points == 1 && *nb_points_total == 1 {
        if pt1.check_same_sp(sp1) != 1 {
            *sp2 = *pt1;
            *nb_points_total = 2;
        }
    } else if nb_points == 1 && *nb_points_total == 2 {
        if sp1.check_same_sp(pt1) != 0 || sp2.check_same_sp(pt1) != 0 {
            *nb_points_total = 2;
        } else {
            *nb_points_total = 3;
        }
    } else if nb_points == 2 && *nb_points_total == 0 {
        *sp1 = *pt1;
        *sp2 = *pt2;
        *nb_points_total = 2;
    } else if nb_points == 2 && *nb_points_total == 1 {
        // There is also Pt1 != Pt2.
        if sp1.check_same_sp(pt1) != 0 {
            *sp2 = *pt2;
            *nb_points_total = 2;
        } else if sp1.check_same_sp(pt2) != 0 {
            *sp2 = *pt1;
            *nb_points_total = 2;
        } else {
            *nb_points_total = 3; // case SP1 != Pt1 && SP1 != Pt2
        }
    } else if nb_points == 2 && *nb_points_total == 2 {
        // There is also SP1 != SP2.
        if sp1.check_same_sp(pt1) != 0 || sp1.check_same_sp(pt2) != 0 {
            if sp2.check_same_sp(pt1) != 0 || sp2.check_same_sp(pt2) != 0 {
                *nb_points_total = 2;
            } else {
                *nb_points_total = 3;
            }
        } else {
            *nb_points_total = 3;
        }
    }
}

impl MaillageAffinage<'_> {
    /// `TriContact(P1, P2, P3, Q1, Q2, Q3, Angle)` (`cxx:1466`).
    pub fn tri_contact(
        &self,
        p1: &IntPolyhPoint,
        p2: &IntPolyhPoint,
        p3: &IntPolyhPoint,
        q1: &IntPolyhPoint,
        q2: &IntPolyhPoint,
        q3: &IntPolyhPoint,
        angle: &mut f64,
    ) -> bool {
        tri_contact(p1, p2, p3, q1, q2, q3, angle)
    }

    /// `StartingPointsResearch(T1, T2, SP1, SP2)` (`cxx:1767-1897`): number of
    /// start points found (0, 1, 2 or 3 when the triangles are coplanar).
    pub fn starting_points_research(
        &self,
        t1: i32,
        t2: i32,
        sp1: &mut IntPolyhStartPoint,
        sp2: &mut IntPolyhStartPoint,
    ) -> i32 {
        let tri1 = self.triangles1.get(t1);
        let tri2 = self.triangles2.get(t2);

        let (p1, p2, p3) = super::triangle_points(&tri1, &self.points1);
        let (q1, q2, q3) = super::triangle_points(&tri2, &self.points2);

        let e1 = p2.sub(&p1);
        let e2 = p3.sub(&p2);
        let e3 = p1.sub(&p3);

        let f1 = q2.sub(&q1);
        let f2 = q3.sub(&q2);
        let f3 = q1.sub(&q3);

        let mut nn1 = IntPolyhPoint::default();
        nn1.cross(&e1, &e2); // normal to the first triangle
        let mut mm1 = IntPolyhPoint::default();
        mm1.cross(&f1, &f2); // normal to the second triangle

        let nn1modulus = nn1.square_modulus().sqrt();
        let mm1modulus = mm1.square_modulus().sqrt();

        let mut nb_points_total = 0;

        // Check T1 normal.
        if nn1modulus.abs() < super::MY_CONFUSION_PRECISION {
            // `cxx:1809-1811`: nothing to do.
        } else {
            let n1 = nn1.divide(nn1modulus);
            if nb_points_total < 3 {
                let (mut pt1, mut pt2) = (
                    IntPolyhStartPoint::default(),
                    IntPolyhStartPoint::default(),
                );
                let nb_points = self.triangle_edge_contact(
                    1, 1, &tri1, &tri2, &p1, &p2, &p3, &e1, &e2, &e3, &q1, &q2, &f1, &n1, &mut pt1,
                    &mut pt2,
                );
                test_nb_points(1, nb_points, &mut nb_points_total, &pt1, &pt2, sp1, sp2);
            }
            if nb_points_total < 3 {
                let (mut pt1, mut pt2) = (
                    IntPolyhStartPoint::default(),
                    IntPolyhStartPoint::default(),
                );
                let nb_points = self.triangle_edge_contact(
                    1, 2, &tri1, &tri2, &p1, &p2, &p3, &e1, &e2, &e3, &q2, &q3, &f2, &n1, &mut pt1,
                    &mut pt2,
                );
                test_nb_points(1, nb_points, &mut nb_points_total, &pt1, &pt2, sp1, sp2);
            }
            if nb_points_total < 3 {
                let (mut pt1, mut pt2) = (
                    IntPolyhStartPoint::default(),
                    IntPolyhStartPoint::default(),
                );
                let nb_points = self.triangle_edge_contact(
                    1, 3, &tri1, &tri2, &p1, &p2, &p3, &e1, &e2, &e3, &q3, &q1, &f3, &n1, &mut pt1,
                    &mut pt2,
                );
                test_nb_points(1, nb_points, &mut nb_points_total, &pt1, &pt2, sp1, sp2);
            }
        }

        // Check T2 normal.
        if mm1modulus.abs() < super::MY_CONFUSION_PRECISION {
            // `cxx:1844-1846`: nothing to do.
        } else {
            let m1 = mm1.divide(mm1modulus);
            if nb_points_total < 3 {
                let (mut pt1, mut pt2) = (
                    IntPolyhStartPoint::default(),
                    IntPolyhStartPoint::default(),
                );
                let nb_points = self.triangle_edge_contact(
                    2, 1, &tri1, &tri2, &q1, &q2, &q3, &f1, &f2, &f3, &p1, &p2, &e1, &m1, &mut pt1,
                    &mut pt2,
                );
                test_nb_points(2, nb_points, &mut nb_points_total, &pt1, &pt2, sp1, sp2);
            }
            if nb_points_total < 3 {
                let (mut pt1, mut pt2) = (
                    IntPolyhStartPoint::default(),
                    IntPolyhStartPoint::default(),
                );
                let nb_points = self.triangle_edge_contact(
                    2, 2, &tri1, &tri2, &q1, &q2, &q3, &f1, &f2, &f3, &p2, &p3, &e2, &m1, &mut pt1,
                    &mut pt2,
                );
                test_nb_points(2, nb_points, &mut nb_points_total, &pt1, &pt2, sp1, sp2);
            }
            if nb_points_total < 3 {
                let (mut pt1, mut pt2) = (
                    IntPolyhStartPoint::default(),
                    IntPolyhStartPoint::default(),
                );
                let nb_points = self.triangle_edge_contact(
                    2, 3, &tri1, &tri2, &q1, &q2, &q3, &f1, &f2, &f3, &p3, &p1, &e3, &m1, &mut pt1,
                    &mut pt2,
                );
                test_nb_points(2, nb_points, &mut nb_points_total, &pt1, &pt2, sp1, sp2);
            }
        }

        // `cxx:1878-1892`.
        if nb_points_total == 2 {
            sp1.set_couple_value(t1, t2);
            sp2.set_couple_value(t1, t2);
        } else if nb_points_total == 1 || nb_points_total == 3 {
            sp1.set_couple_value(t1, t2);
        }

        nb_points_total
    }

    /// `NextStartingPointsResearch(T1, T2, SPInit, SPNext)` (`cxx:1902-2048`).
    pub fn next_starting_points_research(
        &self,
        t1: i32,
        t2: i32,
        sp_init: &IntPolyhStartPoint,
        sp_next: &mut IntPolyhStartPoint,
    ) -> i32 {
        let mut nb_points_total = 0;
        let edge_init1 = sp_init.e1;
        let edge_init2 = sp_init.e2;

        if t1 < 0 || t2 < 0 {
            nb_points_total = 0;
        } else {
            let tri1 = self.triangles1.get(t1);
            let tri2 = self.triangles2.get(t2);

            let (p1, p2, p3) = super::triangle_points(&tri1, &self.points1);
            let (q1, q2, q3) = super::triangle_points(&tri2, &self.points2);

            let e1 = p2.sub(&p1);
            let e2 = p3.sub(&p2);
            let e3 = p1.sub(&p3);

            let f1 = q2.sub(&q1);
            let f2 = q3.sub(&q2);
            let f3 = q1.sub(&q3);

            let mut nn1 = IntPolyhPoint::default();
            nn1.cross(&e1, &e2); // normal to the first triangle
            let mut mm1 = IntPolyhPoint::default();
            mm1.cross(&f1, &f2); // normal to the second triangle

            let nn1modulus = nn1.square_modulus().sqrt();
            let mm1modulus = mm1.square_modulus().sqrt();

            let mut sp1 = IntPolyhStartPoint::default();
            let mut sp2 = IntPolyhStartPoint::default();

            // Check T1 normal.
            if nn1modulus.abs() < super::MY_CONFUSION_PRECISION {
                // `cxx:1956-1958`: nothing to do.
            } else {
                let n1 = nn1.divide(nn1modulus);
                if nb_points_total < 3 && edge_init2 != tri2.first_edge() {
                    let (mut pt1, mut pt2) = (
                        IntPolyhStartPoint::default(),
                        IntPolyhStartPoint::default(),
                    );
                    let nb_points = self.triangle_edge_contact(
                        1, 1, &tri1, &tri2, &p1, &p2, &p3, &e1, &e2, &e3, &q1, &q2, &f1, &n1,
                        &mut pt1, &mut pt2,
                    );
                    test_nb_points(
                        1,
                        nb_points,
                        &mut nb_points_total,
                        &pt1,
                        &pt2,
                        &mut sp1,
                        &mut sp2,
                    );
                }
                if nb_points_total < 3 && edge_init2 != tri2.second_edge() {
                    let (mut pt1, mut pt2) = (
                        IntPolyhStartPoint::default(),
                        IntPolyhStartPoint::default(),
                    );
                    let nb_points = self.triangle_edge_contact(
                        1, 2, &tri1, &tri2, &p1, &p2, &p3, &e1, &e2, &e3, &q2, &q3, &f2, &n1,
                        &mut pt1, &mut pt2,
                    );
                    test_nb_points(
                        1,
                        nb_points,
                        &mut nb_points_total,
                        &pt1,
                        &pt2,
                        &mut sp1,
                        &mut sp2,
                    );
                }
                if nb_points_total < 3 && edge_init2 != tri2.third_edge() {
                    let (mut pt1, mut pt2) = (
                        IntPolyhStartPoint::default(),
                        IntPolyhStartPoint::default(),
                    );
                    let nb_points = self.triangle_edge_contact(
                        1, 3, &tri1, &tri2, &p1, &p2, &p3, &e1, &e2, &e3, &q3, &q1, &f3, &n1,
                        &mut pt1, &mut pt2,
                    );
                    test_nb_points(
                        1,
                        nb_points,
                        &mut nb_points_total,
                        &pt1,
                        &pt2,
                        &mut sp1,
                        &mut sp2,
                    );
                }
            }

            // Check T2 normal.
            if mm1modulus.abs() < super::MY_CONFUSION_PRECISION {
                // `cxx:1988-1990`: nothing to do.
            } else {
                let m1 = mm1.divide(mm1modulus);
                if nb_points_total < 3 && edge_init1 != tri1.first_edge() {
                    let (mut pt1, mut pt2) = (
                        IntPolyhStartPoint::default(),
                        IntPolyhStartPoint::default(),
                    );
                    let nb_points = self.triangle_edge_contact(
                        2, 1, &tri1, &tri2, &q1, &q2, &q3, &f1, &f2, &f3, &p1, &p2, &e1, &m1,
                        &mut pt1, &mut pt2,
                    );
                    test_nb_points(
                        2,
                        nb_points,
                        &mut nb_points_total,
                        &pt1,
                        &pt2,
                        &mut sp1,
                        &mut sp2,
                    );
                }
                if nb_points_total < 3 && edge_init1 != tri1.second_edge() {
                    let (mut pt1, mut pt2) = (
                        IntPolyhStartPoint::default(),
                        IntPolyhStartPoint::default(),
                    );
                    let nb_points = self.triangle_edge_contact(
                        2, 2, &tri1, &tri2, &q1, &q2, &q3, &f1, &f2, &f3, &p2, &p3, &e2, &m1,
                        &mut pt1, &mut pt2,
                    );
                    test_nb_points(
                        2,
                        nb_points,
                        &mut nb_points_total,
                        &pt1,
                        &pt2,
                        &mut sp1,
                        &mut sp2,
                    );
                }
                if nb_points_total < 3 && edge_init1 != tri1.third_edge() {
                    let (mut pt1, mut pt2) = (
                        IntPolyhStartPoint::default(),
                        IntPolyhStartPoint::default(),
                    );
                    let nb_points = self.triangle_edge_contact(
                        2, 3, &tri1, &tri2, &q1, &q2, &q3, &f1, &f2, &f3, &p3, &p1, &e3, &m1,
                        &mut pt1, &mut pt2,
                    );
                    test_nb_points(
                        2,
                        nb_points,
                        &mut nb_points_total,
                        &pt1,
                        &pt2,
                        &mut sp1,
                        &mut sp2,
                    );
                }
            }

            // `cxx:2019-2043`.
            if nb_points_total == 1 {
                if sp1.check_same_sp(sp_init) != 0 {
                    nb_points_total = 0;
                } else {
                    *sp_next = sp1;
                }
            } else if nb_points_total == 2 && sp1.check_same_sp(sp_init) != 0 {
                nb_points_total = 1; // SP1 and SPInit are identical
                *sp_next = sp2;
            } else if nb_points_total == 2 && sp2.check_same_sp(sp_init) != 0 {
                nb_points_total = 1; // SP2 and SPInit are identical
                *sp_next = sp1;
            }
        }

        sp_next.set_couple_value(t1, t2);
        nb_points_total
    }

    /// `TriangleCompare()` (`cxx:3141-3185`): fills the couples of triangles
    /// that interfere and returns their number.
    pub fn triangle_compare(&mut self) -> i32 {
        // Find couples with interfering bounding boxes.
        let interfering = super::refine::get_interfering_triangles(
            &mut self.triangles1,
            &self.points1,
            &mut self.triangles2,
            &self.points2,
        );
        if interfering.is_empty() {
            return 0;
        }

        let mut couple_angle = -2.0;
        for (i_s1, list) in interfering {
            let (p1, p2, p3) = super::triangle_points(&self.triangles1.get(i_s1), &self.points1);

            for i_s2 in list {
                let (q1, q2, q3) =
                    super::triangle_points(&self.triangles2.get(i_s2), &self.points2);

                if tri_contact(&p1, &p2, &p3, &q1, &q2, &q3, &mut couple_angle) {
                    self.contacts
                        .push(IntPolyhCouple::new(i_s1, i_s2, couple_angle));
                    self.triangles1.get_mut(i_s1).has_intersection = true;
                    self.triangles2.get_mut(i_s2).has_intersection = true;
                }
            }
        }

        self.contacts.len() as i32
    }
}
