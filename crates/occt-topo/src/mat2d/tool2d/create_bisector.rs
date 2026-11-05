//! `MAT2d_Tool2d::CreateBisector` (`MAT2d_Tool2d.cxx:277-410`) and
//! `TrimBisec` (`cxx:419-487`).
//!
//! The bisector between two items of the circuit is built by `Bisector_Bisec`
//! from the issue point and the two tangents the `MAT` nodes already carry,
//! then restricted by the half-lines that start at the connexions when one of
//! the items is a circle (the restriction the propagation of the algorithm
//! needs to stay inside its area).

use std::sync::Arc;

use occt_core::gp::GpDir2d;
use occt_core::precision::INFINITE;
use occt_geom2d::bisector::BisectorBisec;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use super::trim::set_trim;
use super::{are_neighbours, Mat2dTool2d, TOLCONF};
use crate::mat::MatBisector;
use crate::mat2d::item::Mat2dItem;

impl Mat2dTool2d {
    /// `CreateBisector(abisector)` (`cxx:277-410`) as `MAT2d_Mat2d` calls it
    /// (`MAT2d_Mat2d.cxx:250` / `:981`): the issue point and the two vectors are
    /// read from the bisector, and the resulting bisector number and
    /// `Sense(1.)` are written back (`cxx:392-395`).
    ///
    /// `edge1` / `edge2` are the `EdgeNumber()` of the bisector's two edges,
    /// which the caller resolves out of its own edge arena.
    pub fn create_bisector_of(
        &mut self,
        bisector: &mut MatBisector,
        edge1: i32,
        edge2: i32,
    ) -> i32 {
        let issue_point = bisector.issue_point();
        let first_vector = bisector.first_vector();
        let second_vector = bisector.second_vector();
        let number =
            self.create_bisector(edge1, edge2, issue_point, first_vector, second_vector);
        bisector.set_bisector_number(number);
        bisector.set_sense(1.0);
        number
    }

    /// `CreateBisector(abisector)` (`cxx:277-410`).
    ///
    /// `edge1` / `edge2` are `abisector->FirstEdge()->EdgeNumber()` and
    /// `abisector->SecondEdge()->EdgeNumber()`, `issue_point` / `first_vector`
    /// / `second_vector` the indices `MAT2d_Mat2d` stored in the bisector. The
    /// returned index is the one `abisector->BisectorNumber()` takes
    /// (`cxx:393`).
    pub fn create_bisector(
        &mut self,
        edge1: i32,
        edge2: i32,
        issue_point: i32,
        first_vector: i32,
        second_vector: i32,
    ) -> i32 {
        thread_local! {
            static DEPTH: std::cell::Cell<u32> = const { std::cell::Cell::new(0) };
        }
        struct Guard;
        impl Drop for Guard {
            fn drop(&mut self) {
                DEPTH.with(|d| d.set(d.get() - 1));
            }
        }
        DEPTH.with(|d| d.set(d.get() + 1));
        let _guard = Guard;
        if DEPTH.with(|d| d.get()) > 40 {
            panic!("create_bisector recursion depth edge1={edge1} edge2={edge2}");
        }
        let circuit = self.circuit_rc();
        let nb_items = circuit.number_of_items();
        let mut ontheline = are_neighbours(edge1, edge2, nb_items);
        let initial_neighbour = ontheline;

        // `cxx:288-291`.
        if circuit.connexion_on(edge2) {
            ontheline = false;
        }

        let elt1 = circuit.value(edge1);
        let elt2 = circuit.value(edge2);

        let issue = self.geom_pnt(issue_point);
        let vector1 = self.geom_vec(first_vector);
        let vector2 = self.geom_vec(second_vector);
        let direction = self.direction();
        let join_type = self.join_type();
        let tolerance = TOLCONF;

        let mut bisector = BisectorBisec::new();
        match (elt1, elt2) {
            // `cxx:332-343`.
            (Mat2dItem::Curve(item1), Mat2dItem::Curve(item2)) => bisector
                .perform_curve_curve(
                    item1, item2, &issue, &vector1, &vector2, direction, join_type, tolerance,
                    ontheline,
                )
                .expect("MAT2d_Tool2d::CreateBisector: StdFail_NotDone"),
            // `cxx:344-358`.
            (Mat2dItem::Point(point1), Mat2dItem::Point(point2)) => bisector
                .perform_point_point(
                    point1, point2, &issue, &vector1, &vector2, direction, tolerance, ontheline,
                )
                .expect("MAT2d_Tool2d::CreateBisector: StdFail_NotDone"),
            // `cxx:359-370`.
            (Mat2dItem::Point(point1), Mat2dItem::Curve(item2)) => bisector
                .perform_point_curve(
                    point1, item2, &issue, &vector1, &vector2, direction, tolerance, ontheline,
                )
                .expect("MAT2d_Tool2d::CreateBisector: StdFail_NotDone"),
            // `cxx:371-383`.
            (Mat2dItem::Curve(item1), Mat2dItem::Point(point2)) => bisector
                .perform_curve_point(
                    item1, point2, &issue, &vector1, &vector2, direction, tolerance, ontheline,
                )
                .expect("MAT2d_Tool2d::CreateBisector: StdFail_NotDone"),
        }

        // `cxx:385-388`: restriction of the bisector.
        self.trim_bisec(&mut bisector, edge1, initial_neighbour, 1);
        self.trim_bisec(&mut bisector, edge2, initial_neighbour, 2);

        self.the_number_of_bisectors += 1;
        self.the_geom_bisectors
            .insert(self.the_number_of_bisectors, bisector);
        self.the_number_of_bisectors
    }

    /// `TrimBisec(B1, IndexEdge, InitialNeighbour, StartOrEnd)`
    /// (`cxx:419-487`): the bisector is cut by the half-line that starts at the
    /// connexion of a circular item.
    fn trim_bisec(
        &self,
        b1: &mut BisectorBisec,
        index_edge: i32,
        initial_neighbour: bool,
        start_or_end: i32,
    ) {
        let circuit = self.circuit_rc();
        let nb_items = circuit.number_of_items();
        let i_next = if index_edge == nb_items {
            1
        } else {
            index_edge + 1
        };

        // `cxx:435-438`: a point item has no basis curve.
        let Some(edge) = circuit.value(index_edge).curve() else {
            return;
        };
        if initial_neighbour {
            return;
        }

        // `cxx:440-439`: the basis of the item, which must be an exact circle.
        let basis = edge.basis_curve();
        let Some(circle) = basis.gp_circ2d() else {
            return;
        };
        let origin = circle.location();

        let mut line1: Option<Geom2dLine> = None;
        let mut line2: Option<Geom2dLine> = None;

        // `cxx:445-454`.
        if circuit.connexion_on(index_edge) && start_or_end == 1 {
            let p2 = circuit.connexion(index_edge).borrow().point_on_first();
            let dir = GpDir2d::new(p2.x() - origin.x(), p2.y() - origin.y())
                .expect("MAT2d_Tool2d::TrimBisec: null direction");
            line1 = Some(Geom2dLine::from_pnt_dir(origin, dir));
        }
        // `cxx:456-464`.
        if circuit.connexion_on(i_next) && start_or_end == 2 {
            let p2 = circuit.connexion(i_next).borrow().point_on_second();
            let dir = GpDir2d::new(p2.x() - origin.x(), p2.y() - origin.y())
                .expect("MAT2d_Tool2d::TrimBisec: null direction");
            line2 = Some(Geom2dLine::from_pnt_dir(origin, dir));
        }

        // `cxx:466-469`.
        if line1.is_none() && line2.is_none() {
            return;
        }

        // `cxx:474-484`.
        for line in [line1, line2].into_iter().flatten() {
            let basis: Arc<dyn Curve2d> = Arc::new(line);
            let trimmed = Geom2dTrimmedCurve::new(basis, 0.0, INFINITE);
            set_trim(b1, &trimmed);
        }
    }
}
