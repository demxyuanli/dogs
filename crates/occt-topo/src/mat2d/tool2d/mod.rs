//! `MAT2d_Tool2d` (`MAT2d_Tool2d.cxx:82-1481`, `.hxx:41-177`): the geometry of
//! the bisecting locus that `MAT2d_Mat2d` walks.
//!
//! The tool wraps the circuit built by `MAT2d_Circuit` and owns the three
//! geometric tables the algorithm indexes into: `theGeomPnts` (points),
//! `theGeomVecs` (tangents) and `theGeomBisectors` (the `Bisector_Bisec` of
//! every `MAT_Bisector`). The indices handed out by `FirstPoint`,
//! `TangentBefore`, `TangentAfter` and `Tangent` are the ones stored in the
//! `MAT_Bisector` / `MAT_Node` fields of the `MAT` graph.
//!
//! `MAT_Bisector` stores its edges as graph ids, not handles, so the edge
//! numbers `MAT2d_Tool2d.cxx` reads with `FirstEdge()->EdgeNumber()` are
//! resolved by the caller and passed in: see [`Mat2dTool2d::create_bisector`]
//! and [`Mat2dTool2d::intersect_bisector`].

mod create_bisector;
mod intersect;
mod trim;

use std::collections::BTreeMap;
use std::rc::Rc;

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::kernel::geomabs::JoinType;
use occt_geom2d::bisector::BisectorBisec;
use occt_geom2d::curve::Curve2d;

use crate::mat::MatSide;
use crate::mat2d::circuit::Mat2dCircuit;
use crate::mat2d::item::Mat2dItem;

/// `MAT2d_TOLCONF` (`MAT2d_Tool2d.cxx:79`).
pub(crate) const TOLCONF: f64 = 1.0e-7;

/// `class MAT2d_Tool2d` (`MAT2d_Tool2d.hxx:41-177`).
#[derive(Clone)]
pub struct Mat2dTool2d {
    /// `theDirection` (`hxx:167`).
    the_direction: f64,
    /// `theJoinType` (`hxx:168`).
    the_join_type: JoinType,
    /// `theNumberOfBisectors` (`hxx:169`).
    the_number_of_bisectors: i32,
    /// `theNumberOfPnts` (`hxx:170`).
    the_number_of_pnts: i32,
    /// `theNumberOfVecs` (`hxx:171`).
    the_number_of_vecs: i32,
    /// `theCircuit` (`hxx:172`): shared with the caller, which reads it back
    /// with `Circuit()` (`BRepMAT2d_BisectingLocus.cxx:154`).
    the_circuit: Option<Rc<Mat2dCircuit>>,
    /// `theGeomBisectors` (`hxx:173`).
    the_geom_bisectors: BTreeMap<i32, BisectorBisec>,
    /// `theGeomPnts` (`hxx:174`).
    the_geom_pnts: BTreeMap<i32, GpPnt2d>,
    /// `theGeomVecs` (`hxx:175`).
    the_geom_vecs: BTreeMap<i32, GpVec2d>,
    /// `theLinesLength` (`hxx:176`).
    the_lines_length: Vec<i32>,
}

impl Default for Mat2dTool2d {
    /// `MAT2d_Tool2d()` (`cxx:84-91`).
    fn default() -> Self {
        Self::new()
    }
}

impl std::fmt::Debug for Mat2dTool2d {
    /// `Bisector_Bisec` holds a trimmed curve, which is not `Debug`; the dump of
    /// the tool only reports the size of the three tables.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Mat2dTool2d")
            .field("theDirection", &self.the_direction)
            .field("theJoinType", &self.the_join_type)
            .field("theNumberOfBisectors", &self.the_number_of_bisectors)
            .field("theNumberOfPnts", &self.the_number_of_pnts)
            .field("theNumberOfVecs", &self.the_number_of_vecs)
            .field("theGeomBisectors", &self.the_geom_bisectors.len())
            .field("theGeomPnts", &self.the_geom_pnts.len())
            .field("theGeomVecs", &self.the_geom_vecs.len())
            .field("theLinesLength", &self.the_lines_length.len())
            .finish()
    }
}

impl Mat2dTool2d {
    /// `MAT2d_Tool2d()` (`cxx:84-91`): the direction is `1.` and the join type
    /// `GeomAbs_Arc`.
    pub fn new() -> Self {
        Self {
            the_direction: 1.0,
            the_join_type: JoinType::Arc,
            the_number_of_bisectors: 0,
            the_number_of_pnts: 0,
            the_number_of_vecs: 0,
            the_circuit: None,
            the_geom_bisectors: BTreeMap::new(),
            the_geom_pnts: BTreeMap::new(),
            the_geom_vecs: BTreeMap::new(),
            the_lines_length: Vec::new(),
        }
    }

    /// `InitItems(EquiCircuit)` (`cxx:94-107`): every table is cleared and the
    /// circuit is kept.
    pub fn init_items(&mut self, equi_circuit: Rc<Mat2dCircuit>) {
        self.the_geom_bisectors.clear();
        self.the_geom_pnts.clear();
        self.the_geom_vecs.clear();
        self.the_lines_length.clear();
        self.the_number_of_bisectors = 0;
        self.the_number_of_vecs = 0;
        self.the_number_of_pnts = 0;
        self.the_circuit = Some(equi_circuit);
    }

    /// `Sense(aSide)` (`cxx:109-119`).
    pub fn sense(&mut self, a_side: MatSide) {
        if a_side == MatSide::Left {
            self.the_direction = 1.0;
        } else {
            self.the_direction = -1.0;
        }
    }

    /// `SetJoinType(aJoinType)` (`cxx:123-126`).
    pub fn set_join_type(&mut self, a_join_type: JoinType) {
        self.the_join_type = a_join_type;
    }

    /// `NumberOfItems()` (`cxx:131-134`).
    pub fn number_of_items(&self) -> i32 {
        self.circuit().number_of_items()
    }

    /// `ToleranceOfConfusion()` (`cxx:138-141`).
    pub fn tolerance_of_confusion(&self) -> f64 {
        2.0 * TOLCONF
    }

    /// `FirstPoint(anitem, dist)` (`cxx:145-174`): the middle of the connexion
    /// that starts `anitem`, or the start of the item itself.
    pub fn first_point(&mut self, anitem: i32, dist: &mut f64) -> i32 {
        self.the_number_of_pnts += 1;

        let circuit = self.circuit_rc();
        if circuit.connexion_on(anitem) {
            let connexion = circuit.connexion(anitem).borrow();
            let p1 = connexion.point_on_first();
            let p2 = connexion.point_on_second();
            let middle = GpPnt2d::new((p1.x() + p2.x()) * 0.5, (p1.y() + p2.y()) * 0.5);
            *dist = p1.distance(&p2) * 0.5;
            self.the_geom_pnts.insert(self.the_number_of_pnts, middle);
            return self.the_number_of_pnts;
        }

        *dist = 0.0;
        let point = match circuit.value(anitem) {
            Mat2dItem::Curve(c) => c.d0(c.first_parameter()),
            Mat2dItem::Point(p) => *p,
        };
        self.the_geom_pnts.insert(self.the_number_of_pnts, point);
        self.the_number_of_pnts
    }

    /// `TangentBefore(anitem, IsOpenResult)` (`cxx:178-217`): the tangent at
    /// the end of `anitem`, i.e. the vector of the connexion that follows it.
    pub fn tangent_before(&mut self, anitem: i32, is_open_result: bool) -> i32 {
        self.the_number_of_vecs += 1;

        let circuit = self.circuit_rc();
        let nb_items = circuit.number_of_items();
        let item = if !is_open_result {
            if anitem == nb_items {
                1
            } else {
                anitem + 1
            }
        } else if anitem == nb_items {
            anitem - 1
        } else {
            anitem + 1
        };

        let vector = if circuit.connexion_on(item) {
            let c = circuit.connexion(item).borrow();
            let p1 = c.point_on_first();
            let p2 = c.point_on_second();
            GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y())
        } else {
            match circuit.value(anitem) {
                Mat2dItem::Curve(c) => c.eval_dn(c.last_parameter(), 1),
                Mat2dItem::Point(_) => {
                    let c = circuit.value(item).expect_curve();
                    let parameter = if is_open_result && anitem == nb_items {
                        c.last_parameter()
                    } else {
                        c.first_parameter()
                    };
                    c.eval_dn(parameter, 1)
                }
            }
        };

        self.the_geom_vecs.insert(self.the_number_of_vecs, vector);
        self.the_number_of_vecs
    }

    /// `TangentAfter(anitem, IsOpenResult)` (`cxx:221-262`): the reversed
    /// tangent at the origin of `anitem`.
    pub fn tangent_after(&mut self, anitem: i32, is_open_result: bool) -> i32 {
        self.the_number_of_vecs += 1;

        let circuit = self.circuit_rc();
        if circuit.connexion_on(anitem) {
            let c = circuit.connexion(anitem).borrow();
            let p1 = c.point_on_first();
            let p2 = c.point_on_second();
            self.the_geom_vecs
                .insert(self.the_number_of_vecs, GpVec2d::new(p1.x() - p2.x(), p1.y() - p2.y()));
            return self.the_number_of_vecs;
        }

        let nb_items = circuit.number_of_items();
        let vector = match circuit.value(anitem) {
            Mat2dItem::Curve(c) => c.eval_dn(c.first_parameter(), 1),
            Mat2dItem::Point(_) => {
                let item = if !is_open_result {
                    if anitem == 1 {
                        nb_items
                    } else {
                        anitem - 1
                    }
                } else if anitem == 1 {
                    2
                } else {
                    anitem - 1
                };
                let c = circuit.value(item).expect_curve();
                let parameter = if is_open_result && anitem == 1 {
                    c.first_parameter()
                } else {
                    c.last_parameter()
                };
                c.eval_dn(parameter, 1)
            }
        };

        self.the_geom_vecs
            .insert(self.the_number_of_vecs, vector.reversed());
        self.the_number_of_vecs
    }

    /// `Tangent(bisector)` (`cxx:266-272`): the tangent at the end of the
    /// bisector.
    pub fn tangent(&mut self, bisector: i32) -> i32 {
        self.the_number_of_vecs += 1;
        let vector = {
            let b = self
                .geom_bis(bisector)
                .value()
                .expect("MAT2d_Tool2d::Tangent: null thebisector");
            b.eval_dn(b.last_parameter(), 1)
        };
        self.the_geom_vecs.insert(self.the_number_of_vecs, vector);
        self.the_number_of_vecs
    }

    /// `Distance(abisector, param1, param2)` (`cxx:1106-1120`).
    pub fn distance(&self, abisector: &crate::mat::MatBisector, param1: f64, param2: f64) -> f64 {
        let mut dist = infinite();
        if !is_infinite(param1) && !is_infinite(param2) {
            let b = self
                .geom_bis(abisector.bisector_number())
                .value()
                .expect("MAT2d_Tool2d::Distance: null thebisector");
            let p1 = b.d0(param1);
            let p2 = b.d0(param2);
            dist = p1.distance(&p2);
        }
        dist
    }

    /// `GeomBis(Index)` (`cxx:1146-1149`).
    pub fn geom_bis(&self, index: i32) -> &BisectorBisec {
        self.the_geom_bisectors
            .get(&index)
            .unwrap_or_else(|| panic!("MAT2d_Tool2d::GeomBis: {index} is not bound"))
    }

    /// `ChangeGeomBis(Index)` (`cxx:1153-1156`).
    pub fn change_geom_bis(&mut self, index: i32) -> &mut BisectorBisec {
        self.the_geom_bisectors
            .get_mut(&index)
            .unwrap_or_else(|| panic!("MAT2d_Tool2d::ChangeGeomBis: {index} is not bound"))
    }

    /// `GeomElt(Index)` (`cxx:1160-1163`): `theCircuit->Value(Index)`.
    pub fn geom_elt(&self, index: i32) -> &Mat2dItem {
        self.circuit().value(index)
    }

    /// `GeomPnt(Index)` (`cxx:1167-1170`).
    pub fn geom_pnt(&self, index: i32) -> GpPnt2d {
        *self
            .the_geom_pnts
            .get(&index)
            .unwrap_or_else(|| panic!("MAT2d_Tool2d::GeomPnt: {index} is not bound"))
    }

    /// `GeomVec(Index)` (`cxx:1174-1177`).
    pub fn geom_vec(&self, index: i32) -> GpVec2d {
        *self
            .the_geom_vecs
            .get(&index)
            .unwrap_or_else(|| panic!("MAT2d_Tool2d::GeomVec: {index} is not bound"))
    }

    /// `Circuit()` (`cxx:1179-1182`).
    pub fn circuit(&self) -> &Mat2dCircuit {
        self.the_circuit
            .as_deref()
            .expect("MAT2d_Tool2d: no circuit (InitItems was not called)")
    }

    /// `theDirection` (`hxx:167`).
    pub(crate) fn direction(&self) -> f64 {
        self.the_direction
    }

    /// `theJoinType` (`hxx:168`).
    pub(crate) fn join_type(&self) -> JoinType {
        self.the_join_type
    }

    /// `theCircuit` as a shared handle, so that the tables above can be filled
    /// while the circuit is read.
    pub(crate) fn circuit_rc(&self) -> Rc<Mat2dCircuit> {
        Rc::clone(
            self.the_circuit
                .as_ref()
                .expect("MAT2d_Tool2d: no circuit (InitItems was not called)"),
        )
    }
}

/// `AreNeighbours(IEdge1, IEdge2, NbEdge)` (`cxx:1271-1283`): true when the two
/// items are consecutive on a closed contour of `NbEdge` items.
pub(crate) fn are_neighbours(edge1: i32, edge2: i32, nb_edge: i32) -> bool {
    (edge1 - edge2).abs() == 1 || (edge1 - edge2).abs() == nb_edge - 1
}

/// `Precision::IsInfinite(value)`.
pub(crate) fn is_infinite(value: f64) -> bool {
    occt_core::precision::Precision::is_infinite(value)
}

/// `Precision::Infinite()`.
pub(crate) fn infinite() -> f64 {
    occt_core::precision::INFINITE
}
