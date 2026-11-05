//! `MAT2d_Tool2d::Projection` (`MAT2d_Tool2d.cxx:560-652`),
//! `IsSameDistance` (`cxx:655-785`), `IntersectBisector` (`cxx:787-1101`) and
//! `CheckEnds` (`cxx:1385-1415`).
//!
//! `IntersectBisector` is the query `MAT2d_Mat2d` uses to expand its tree: the
//! two bisectors are intersected, then a candidate is kept only when it really
//! is equidistant from the four items the two bisectors separate
//! (`IsSameDistance`) and when it is not the start of either bisector.
//!
//! `MAT_Bisector` stores its edges as graph ids, not handles, so the edge
//! numbers `cxx:816-819` / `:1027-1030` read with `FirstEdge()->EdgeNumber()`
//! are resolved by the caller and passed in.

use std::f64::consts::PI;

use occt_core::elib::clib::adjust_periodic;
use occt_core::gp::{GpDir2d, GpPnt2d};
use occt_core::kernel::geomabs::JoinType;
use occt_core::precision::{CONFUSION, INFINITE};
use occt_geom2d::bisector::{bisector_parameter, BisectorInter};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::extrema2d::ExtremaExtPC2d;

use super::trim::domain;
use super::{are_neighbours, is_infinite, Mat2dTool2d, TOLCONF};
use crate::mat::MatBisector;
use crate::mat2d::item::Mat2dItem;

/// `CheckEnds(Elt, PCom, Distance, Tol)` (`cxx:1385-1415`): true when one of the
/// ends of the item lies at `Distance` from `PCom`, i.e. when `PCom` projects
/// onto the end of the item rather than onto its interior.
pub(crate) fn check_ends(elt: &Mat2dItem, p_com: &GpPnt2d, distance: f64, tol: f64) -> bool {
    match elt {
        // `cxx:1394-1397`: a point item has no end.
        Mat2dItem::Point(_) => false,
        // `cxx:1398-1413`.
        Mat2dItem::Curve(_) => {
            let df = p_com.distance(&elt.start_point());
            let dl = p_com.distance(&elt.end_point());
            (df - distance).abs() <= tol || (dl - distance).abs() <= tol
        }
    }
}

impl Mat2dTool2d {
    /// `Projection(IEdge, PCom, Distance)` (`cxx:560-652`): the distance of
    /// `PCom` to the item, restricted to the part of the item between its
    /// connexions.
    ///
    /// The OCCT out-parameter becomes the second half of the pair; `false` (the
    /// `NbExt() < 1` arm, `cxx:630-633`) leaves `Precision::Infinite()` there.
    pub fn projection(&self, i_edge: i32, p_com: &GpPnt2d) -> (bool, f64) {
        let eps = TOLCONF;
        let circuit = self.circuit_rc();

        // `cxx:570-576`: a point item is its own projection.
        let Some(curve) = circuit.value(i_edge).curve() else {
            return (
                true,
                p_com.distance(&circuit.value(i_edge).point().expect("point item")),
            );
        };

        // `cxx:583-599`: the range of the item between its connexions.
        let mut param_min = curve.first_parameter();
        let mut param_max = curve.last_parameter();
        if circuit.connexion_on(i_edge) {
            param_min = circuit.connexion(i_edge).borrow().parameter_on_second();
        }
        let nb_items = circuit.number_of_items();
        let i_next = if i_edge == nb_items { 1 } else { i_edge + 1 };
        if circuit.connexion_on(i_next) {
            param_max = circuit.connexion(i_next).borrow().parameter_on_first();
            if curve.basis_curve().is_periodic() {
                adjust_periodic(0.0, 2.0 * PI, eps, &mut param_min, &mut param_max);
            }
        }

        // `cxx:602-622`: widen the range, in arc length on a circle.
        let circle = curve.gp_circ2d();
        if let Some(circle) = &circle {
            let radius = circle.radius();
            let mut eps_circ = 100.0 * eps;
            if radius < 1.0 {
                eps_circ = eps / radius;
            }
            if param_max - param_min + 2.0 * eps_circ < 2.0 * PI {
                param_max += eps_circ;
                param_min -= eps_circ;
            }
        } else {
            param_max += eps;
            param_min -= eps;
        }

        // `cxx:626-648`.
        let extremas = ExtremaExtPC2d::new_range(curve, p_com, param_min, param_max);
        if extremas.is_done() {
            if extremas.nb_ext() < 1 {
                return (false, INFINITE);
            }
            let mut distance = INFINITE;
            for i in 1..=extremas.nb_ext() {
                if extremas.square_distance(i) < distance {
                    distance = extremas.square_distance(i);
                }
            }
            (true, distance.sqrt())
        } else if let Some(circle) = circle {
            (true, circle.radius())
        } else {
            (true, INFINITE)
        }
    }

    /// `IsSameDistance(BisectorOne, BisectorTwo, PCom, Distance)`
    /// (`cxx:655-785`): true when `PCom` is at the same distance from the four
    /// items the two bisectors separate.
    ///
    /// `edges_one` / `edges_two` are the `FirstEdge` / `SecondEdge` numbers of
    /// the two bisectors (`cxx:664-667`).
    pub fn is_same_distance(
        &self,
        edges_one: [i32; 2],
        edges_two: [i32; 2],
        p_com: &GpPnt2d,
    ) -> (bool, f64) {
        let eps = 1.0e-7;
        let [i_edge1, i_edge2] = edges_one;
        let [i_edge3, i_edge4] = edges_two;
        let circuit = self.circuit_rc();

        // `cxx:670-696`: a projection that failed because `PCom` sees the item
        // end-on takes the distance of the projection that succeeded.
        let (is_done1, dist1) = self.projection(i_edge1, p_com);
        let (is_done2, dist2) = self.projection(i_edge2, p_com);
        let mut dist = [dist1, dist2, 0.0, 0.0];
        if is_done1 {
            if !is_done2 {
                let tol = CONFUSION.max(eps * dist[0]);
                if check_ends(circuit.value(i_edge2), p_com, dist[0], tol) {
                    dist[1] = dist[0];
                }
            }
        } else if is_done2 {
            let tol = CONFUSION.max(eps * dist[1]);
            if check_ends(circuit.value(i_edge1), p_com, dist[1], tol) {
                dist[0] = dist[1];
            }
        }

        // `cxx:698-724`: the two edges of the second bisector, reusing the
        // distances already known for the first one.
        let (is_done3, dist3) = if i_edge3 == i_edge1 {
            (true, dist[0])
        } else if i_edge3 == i_edge2 {
            (true, dist[1])
        } else {
            self.projection(i_edge3, p_com)
        };
        dist[2] = dist3;
        let (is_done4, dist4) = if i_edge4 == i_edge1 {
            (true, dist[0])
        } else if i_edge4 == i_edge2 {
            (true, dist[1])
        } else {
            self.projection(i_edge4, p_com)
        };
        dist[3] = dist4;

        // `cxx:726-748`.
        if is_done3 {
            if !is_done4 {
                let tol = CONFUSION.max(eps * dist[2]);
                if check_ends(circuit.value(i_edge4), p_com, dist[2], tol) {
                    dist[3] = dist[2];
                }
            }
        } else if is_done4 {
            let tol = CONFUSION.max(eps * dist[3]);
            if check_ends(circuit.value(i_edge3), p_com, dist[3], tol) {
                dist[2] = dist[3];
            }
        }

        // `cxx:757-784`.
        let eps_dist = TOLCONF * 300.0;
        let mut distance = dist[0];
        let intersection_join = self.join_type() == JoinType::Intersection;
        if intersection_join && is_infinite(distance) {
            for value in dist.iter().skip(1) {
                if !is_infinite(*value) {
                    distance = *value;
                    break;
                }
            }
        }
        for value in dist.iter() {
            if intersection_join && is_infinite(*value) {
                continue;
            }
            if (*value - distance).abs() > eps_dist {
                return (false, INFINITE);
            }
        }

        (true, distance)
    }

    /// `IntersectBisector(BisectorOne, BisectorTwo, IntPnt)` (`cxx:787-1101`):
    /// the distance of the retained intersection point, `Precision::Infinite()`
    /// when the two bisectors do not meet in a valid point.
    ///
    /// `IntPnt` (`cxx:788`) becomes the second half of the pair: the index of the
    /// point stored in `theGeomPnts`, `None` when the intersection is rejected.
    pub fn intersect_bisector(
        &mut self,
        bisector_one: &mut MatBisector,
        edges_one: [i32; 2],
        bisector_two: &mut MatBisector,
        edges_two: [i32; 2],
    ) -> (f64, Option<i32>) {
        let tolerance = TOLCONF;
        let [if1, is1] = edges_one;
        let [if2, is2] = edges_two;

        // `cxx:798-808`: `ChangeValue()` only checks that the handles are not
        // null; the geometry is read, never written, so the clones below are
        // equivalent to the OCCT handles.
        let bisec1 = self.geom_bis(bisector_one.bisector_number()).clone();
        let bisec2 = self.geom_bis(bisector_two.bisector_number()).clone();
        let bisector1 = bisec1
            .value()
            .expect("MAT2d_Tool2d::IntersectBisector: null thebisector")
            .clone();
        let bisector2 = bisec2
            .value()
            .expect("MAT2d_Tool2d::IntersectBisector: null thebisector")
            .clone();

        // `cxx:810-834`: two bisectors issued from the two halves of the same
        // connexion meet at its middle, which is not a valid intersection.
        let nb_items = self.number_of_items();
        if are_neighbours(if1, is1, nb_items) && are_neighbours(if2, is2, nb_items) {
            let circuit = self.circuit_rc();
            if circuit.connexion_on(is2) && circuit.connexion_on(is1) {
                let c1 = circuit.connexion(is1).borrow();
                let c2 = circuit.connexion(is2).borrow();
                if c2.index_first_line() == c1.index_second_line()
                    && c1.index_first_line() == c2.index_second_line()
                {
                    return (INFINITE, None);
                }
            }
        }

        // `cxx:836-847`.
        let domain1 = domain(&bisector1, tolerance);
        let domain2 = domain(&bisector2, tolerance);
        if domain1.last_parameter() - domain1.first_parameter() < tolerance {
            return (INFINITE, None);
        }
        if domain2.last_parameter() - domain2.first_parameter() < tolerance {
            return (INFINITE, None);
        }

        // `cxx:874-886`.
        let mut intersect = BisectorInter::new();
        intersect.perform(
            &bisec1,
            &domain1,
            &bisec2,
            &domain2,
            tolerance,
            tolerance,
            true,
        );

        // `cxx:892-904`.
        if !intersect.is_done() || intersect.is_empty() {
            return (INFINITE, None);
        }

        let mut distance_mini = INFINITE;
        let mut param1 = INFINITE;
        let mut param2 = INFINITE;
        let mut solution_valid = false;
        let mut point_solution = GpPnt2d::new(0.0, 0.0);

        // `cxx:907-971`: the small segments are scanned for a point that really
        // is equidistant from the four items.
        let mut distance = 0.0;
        if intersect.nb_segments() >= 1 {
            let max_segment_length = 10.0 * tolerance;
            for i in 1..=intersect.nb_segments() {
                let segment = intersect.segment(i);
                let mut point_retained = false;
                let mut point_on_segment = GpPnt2d::new(0.0, 0.0);

                if segment.has_first_point() && segment.has_last_point() {
                    let p1 = *segment.first_point().value();
                    let p2 = *segment.last_point().value();
                    let segment_length = p1.distance(&p2);
                    if segment_length <= tolerance {
                        point_on_segment = p1;
                        let same = self.is_same_distance(edges_one, edges_two, &point_on_segment);
                        distance = same.1;
                        if same.0 {
                            point_retained = true;
                        }
                    } else if segment_length <= max_segment_length {
                        let dir = GpDir2d::new(p2.x() - p1.x(), p2.y() - p1.y())
                            .expect("MAT2d_Tool2d::IntersectBisector: null direction");
                        let mut dist = 0.0;
                        while dist <= segment_length + tolerance {
                            point_on_segment =
                                GpPnt2d::new(p1.x() + dist * dir.x(), p1.y() + dist * dir.y());
                            let same =
                                self.is_same_distance(edges_one, edges_two, &point_on_segment);
                            distance = same.1;
                            if same.0 {
                                point_retained = true;
                                break;
                            }
                            dist += tolerance;
                        }
                    }
                }

                // `cxx:953-968`: keep the candidate with the smallest parameters
                // on both bisectors.
                if point_retained {
                    let param_a =
                        bisector_parameter(&**bisector1.basis_curve(), &point_on_segment);
                    let param_b =
                        bisector_parameter(&**bisector2.basis_curve(), &point_on_segment);
                    if param_a < param1 && param_b < param2 {
                        param1 = param_a;
                        param2 = param_b;
                        distance_mini = distance;
                        point_solution = point_on_segment;
                        solution_valid = true;
                    }
                }
            }
        }

        // `cxx:972-998`.
        if intersect.nb_points() != 1 {
            for i in 1..=intersect.nb_points() {
                let candidate = intersect.point(i);
                let same = self.is_same_distance(edges_one, edges_two, candidate.value());
                if same.0 && same.1 > tolerance {
                    let param_a = candidate.param_on_first();
                    let param_b = candidate.param_on_second();
                    if param_a < param1 && param_b < param2 {
                        param1 = param_a;
                        param2 = param_b;
                        distance_mini = same.1;
                        point_solution = *candidate.value();
                        solution_valid = true;
                    }
                }
            }
        } else {
            let candidate = intersect.point(1);
            point_solution = *candidate.value();
            param1 = candidate.param_on_first();
            param2 = candidate.param_on_second();
            let same = self.is_same_distance(edges_one, edges_two, &point_solution);
            distance_mini = same.1;
            solution_valid = same.0;
        }

        // `cxx:1000-1007`.
        if !solution_valid {
            return (INFINITE, None);
        }
        self.the_number_of_pnts += 1;
        self.the_geom_pnts
            .insert(self.the_number_of_pnts, point_solution);
        let int_pnt = self.the_number_of_pnts;

        // `cxx:1009-1075`: a solution sitting on the start of one of the
        // bisectors is rejected, unless that bisector is issued from a connexion
        // between neighbouring items.
        let mut extremite_controle = true;
        if self.circuit_rc().connexion_on(is1)
            && are_neighbours(if1, is1, nb_items)
            && are_neighbours(if2, is2, nb_items)
            && is1 == if2
        {
            extremite_controle = false;
            param1 += tolerance;
        }
        if self.circuit_rc().connexion_on(is2)
            && are_neighbours(if1, is1, nb_items)
            && are_neighbours(if2, is2, nb_items)
            && is1 == if2
        {
            extremite_controle = false;
            param2 += tolerance;
        }
        if extremite_controle {
            let start1 = bisector1.d0(bisector1.first_parameter());
            if start1.distance(&point_solution) < tolerance {
                return (INFINITE, None);
            }
            let start2 = bisector2.d0(bisector2.first_parameter());
            if start2.distance(&point_solution) < tolerance {
                return (INFINITE, None);
            }
        }

        // `cxx:1077-1087`.
        if bisector_one.second_parameter() < INFINITE
            && bisector_one.second_parameter() < param1 * (1.0 - tolerance)
        {
            return (INFINITE, None);
        }
        if bisector_two.first_parameter() < INFINITE
            && bisector_two.first_parameter() < param2 * (1.0 - tolerance)
        {
            return (INFINITE, None);
        }

        // `cxx:1089-1090`.
        bisector_one.set_second_parameter(param1);
        bisector_two.set_first_parameter(param2);

        (distance_mini, Some(int_pnt))
    }
}
