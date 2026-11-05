//! `Geom2dHatch_Element` and `Geom2dHatch_Elements`.
//!
//! Sources: `ModelingAlgorithms/TKGeomAlgo/Geom2dHatch/Geom2dHatch_Element.hxx:25-53`,
//! `.cxx:20-72`, `Geom2dHatch_Elements.hxx:34-93`, `.cxx:27-263` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `Geom2dAdaptor_Curve` has no port class; the port's [`Curve2d`] trait
//!   already carries `D0` / `D1` / `FirstParameter` / `LastParameter`, so the
//!   elements store an `Arc<dyn Curve2d>` (the same substitution
//!   [`crate::brep_adaptor_curve2d`] documents).
//! - `NCollection_DataMap<int, Geom2dHatch_Element>` becomes a `BTreeMap`, and
//!   the `Iterator Iter` member becomes a sorted key snapshot plus a cursor
//!   (`init_edges` / `more_edges` / `next_edge` / `current_edge`). OCCT walks
//!   the map in bucket order, the port walks it in ascending key order; element
//!   indices are handed out in ascending order by
//!   `Geom2dHatch_Hatcher::AddElement` (`Hatcher.cxx:157-168`), so the two
//!   orders coincide for every map the hatcher builds.
//! - `Geom2dHatch_Elements(const Geom2dHatch_Elements&)` (`cxx:30-39`) does not
//!   copy the map at all, it only resets the four cursors; that constructor has
//!   no caller, so the port derives `Clone` instead and records the difference
//!   here.

use std::collections::BTreeMap;
use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpLin2d, GpPnt2d, GpVec2d};
use occt_core::precision::{Precision, PCONFUSION};
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;

/// `Probing_Start` (`cxx:26`).
const PROBING_START: f64 = 0.123;
/// `Probing_End` (`cxx:27`).
const PROBING_END: f64 = 0.8;
/// `Probing_Step` (`cxx:28`).
const PROBING_STEP: f64 = 0.2111;
/// `Precision::SquarePConfusion()` (`Precision.hxx:176`).
const SQUARE_PCONFUSION: f64 = PCONFUSION * PCONFUSION;
/// `RealLast()` (`Standard_Real.hxx:179-182`).
const REAL_LAST: f64 = f64::MAX;

/// `Geom2dHatch_Element` (`hxx:25-53`).
#[derive(Clone)]
pub struct Element {
    /// `myCurve`
    my_curve: Arc<dyn Curve2d>,
    /// `myOrientation`
    my_orientation: Orientation,
}

impl Element {
    /// `Geom2dHatch_Element(Curve, Orientation = TopAbs_FORWARD)` (`cxx:24-32`).
    pub fn new(curve: Arc<dyn Curve2d>, orientation: Orientation) -> Self {
        Self {
            my_curve: curve,
            my_orientation: orientation,
        }
    }

    /// `Curve()` (`cxx:38-41`).
    pub fn curve(&self) -> &Arc<dyn Curve2d> {
        &self.my_curve
    }

    /// `ChangeCurve()` (`cxx:47-50`).
    pub fn change_curve(&mut self) -> &mut Arc<dyn Curve2d> {
        &mut self.my_curve
    }

    /// `Orientation(Orientation)` (`cxx:57-60`).
    pub fn set_orientation(&mut self, orientation: Orientation) {
        self.my_orientation = orientation;
    }

    /// `Orientation()` (`cxx:67-70`).
    pub fn orientation(&self) -> Orientation {
        self.my_orientation
    }
}

/// `Geom2dHatch_Elements` (`hxx:34-93`): the boundary of the 2D region.
#[derive(Clone)]
pub struct Elements {
    /// `myMap`
    my_map: BTreeMap<i32, Element>,
    /// The sorted key snapshot standing in for `Iter`.
    iter_keys: Vec<i32>,
    /// The current position in the snapshot standing in for `Iter`.
    iter_pos: usize,
    /// `NumWire`
    num_wire: i32,
    /// `NumEdge` — written by `InitEdges` (`cxx:218`), never read, as in OCCT.
    #[allow(dead_code)]
    num_edge: i32,
    /// `myCurEdge`
    my_cur_edge: i32,
    /// `myCurEdgePar`
    my_cur_edge_par: f64,
}

impl Default for Elements {
    fn default() -> Self {
        Self::new()
    }
}

impl Elements {
    /// `Geom2dHatch_Elements()` (`cxx:41-48`).
    pub fn new() -> Self {
        Self {
            my_map: BTreeMap::new(),
            iter_keys: Vec::new(),
            iter_pos: 0,
            num_wire: 0,
            num_edge: 0,
            my_cur_edge: 1,
            my_cur_edge_par: PROBING_START,
        }
    }

    /// `Clear()` (`cxx:50-53`).
    pub fn clear(&mut self) {
        self.my_map.clear();
        self.iter_keys.clear();
        self.iter_pos = 0;
    }

    /// `Bind(K, I)` (`cxx:66-69`): answers whether the key was free.
    pub fn bind(&mut self, k: i32, i: Element) -> bool {
        self.my_map.insert(k, i).is_none()
    }

    /// `IsBound(K)` (`cxx:55-58`).
    pub fn is_bound(&self, k: i32) -> bool {
        self.my_map.contains_key(&k)
    }

    /// `UnBind(K)` (`cxx:60-63`).
    pub fn unbind(&mut self, k: i32) -> bool {
        self.my_map.remove(&k).is_some()
    }

    /// `Find(K)` (`cxx:71-74`). OCCT raises `OutOfRange` for a free key, the
    /// port panics through the indexing operation.
    pub fn find(&self, k: i32) -> &Element {
        &self.my_map[&k]
    }

    /// `ChangeFind(K)` (`cxx:76-79`).
    pub fn change_find(&mut self, k: i32) -> &mut Element {
        self.my_map
            .get_mut(&k)
            .expect("Geom2dHatch_Elements::ChangeFind")
    }

    /// `CheckPoint(P)` (`cxx:82-86`): the probing point is always accepted.
    pub fn check_point(&mut self, _p: &GpPnt2d) -> bool {
        true
    }

    /// `Reject(P)` (`cxx:91-95`): no early rejection.
    pub fn reject(&self, _p: &GpPnt2d) -> bool {
        false
    }

    /// `Segment(P, L, Par)` (`cxx:97-101`).
    pub fn segment(&mut self, p: &GpPnt2d, l: &mut GpLin2d, par: &mut f64) -> bool {
        self.my_cur_edge = 1;
        self.my_cur_edge_par = PROBING_START;
        self.other_segment(p, l, par)
    }

    /// `OtherSegment(P, L, Par)` (`cxx:103-197`).
    ///
    /// Walks the elements from `myCurEdge` on and probes each one at the
    /// relative parameters `myCurEdgePar + k * Probing_Step` below
    /// `Probing_End`; the first probe whose line is not collinear with the
    /// element and whose endpoints are off that line is accepted. The returned
    /// parameter is the squared distance from `P` to the probed point, except
    /// in the endpoint branch (`cxx:181`) where it is the plain distance.
    pub fn other_segment(&mut self, p: &GpPnt2d, l: &mut GpLin2d, par: &mut f64) -> bool {
        let keys: Vec<i32> = self.my_map.keys().copied().collect();
        for (index, key) in keys.iter().enumerate() {
            let i = (index + 1) as i32;
            if i < self.my_cur_edge {
                continue; // `cxx:111-114`
            }
            // `Geom2dHatch_Element& Item = myMap.ChangeFind(...)` (`cxx:118`);
            // the port clones the `Arc` so that the probing cursors below stay
            // mutable while the curve is read.
            let (orientation, e) = {
                let item = self
                    .my_map
                    .get(key)
                    .expect("Geom2dHatch_Elements::OtherSegment");
                (item.orientation(), Arc::clone(item.curve()))
            };
            if orientation == Orientation::Forward || orientation == Orientation::Reversed {
                let mut a_f_par = e.first_parameter();
                let mut a_l_par = e.last_parameter();
                if Precision::is_negative_infinite(a_f_par) {
                    if Precision::is_positive_infinite(a_l_par) {
                        a_f_par = -1.;
                        a_l_par = 1.;
                    } else {
                        a_f_par = a_l_par - 1.;
                    }
                } else if Precision::is_positive_infinite(a_l_par) {
                    a_l_par = a_f_par + 1.;
                }

                while self.my_cur_edge_par < PROBING_END {
                    let a_param =
                        self.my_cur_edge_par * a_f_par + (1. - self.my_cur_edge_par) * a_l_par;
                    let (a_p_on_c, a_tan_vec) = e.d1(a_param);
                    let a_lin_vec = GpVec2d::new(p.x() - a_p_on_c.x(), p.y() - a_p_on_c.y());
                    *par = a_lin_vec.square_magnitude();

                    let mut skipped = false;
                    if *par > SQUARE_PCONFUSION {
                        let a_lin_dir = GpDir2d::from_vec2d(&a_lin_vec)
                            .expect("Geom2dHatch_Elements::OtherSegment: null direction");
                        let a_tan_mod = a_tan_vec.square_magnitude();
                        if a_tan_mod < SQUARE_PCONFUSION {
                            skipped = true; // `cxx:154-157`
                        } else {
                            let a_tan_unit = a_tan_vec.divided(a_tan_mod.sqrt());
                            let a_sin_a = a_tan_unit.crossed(&GpVec2d::from_dir2d(&a_lin_dir));
                            if a_sin_a.abs() < 0.001 {
                                // Angle too small: line and edge would look
                                // tangent, which breaks the classifier.
                                if self.my_cur_edge_par + PROBING_STEP < PROBING_END {
                                    skipped = true; // `cxx:159-165`
                                }
                            }
                        }
                        if !skipped {
                            *l = GpLin2d::from_pnt_dir(*p, a_lin_dir);
                            let mut a_p_on_c2 = e.d0(a_f_par);
                            if square_distance(l, &a_p_on_c2) > SQUARE_PCONFUSION {
                                a_p_on_c2 = e.d0(a_l_par);
                                if square_distance(l, &a_p_on_c2) > SQUARE_PCONFUSION {
                                    self.my_cur_edge_par += PROBING_STEP;
                                    if self.my_cur_edge_par >= PROBING_END {
                                        self.my_cur_edge += 1;
                                        self.my_cur_edge_par = PROBING_START;
                                    }
                                    *par = par.sqrt();
                                    return true; // `cxx:166-184`
                                }
                            }
                        }
                    }
                    self.my_cur_edge_par += PROBING_STEP;
                }
            }
            self.my_cur_edge += 1;
            self.my_cur_edge_par = PROBING_START;
        }

        *par = REAL_LAST;
        *l = GpLin2d::from_pnt_dir(*p, GpDir2d::new(1., 0.).expect("unit X direction"));
        false // `cxx:194-196`
    }

    /// `InitWires()` (`cxx:202-205`).
    pub fn init_wires(&mut self) {
        self.num_wire = 0;
    }

    /// `MoreWires()` (`cxx:242-245`): the whole element set is one wire.
    pub fn more_wires(&self) -> bool {
        self.num_wire == 0
    }

    /// `NextWire()` (`cxx:248-251`).
    pub fn next_wire(&mut self) {
        self.num_wire += 1;
    }

    /// `RejectWire(L, Par)` (`cxx:207-211`): no early rejection.
    pub fn reject_wire(&self, _l: &GpLin2d, _par: f64) -> bool {
        false
    }

    /// `InitEdges()` (`cxx:216-220`).
    pub fn init_edges(&mut self) {
        self.num_edge = 0;
        self.iter_keys = self.my_map.keys().copied().collect();
        self.iter_pos = 0;
    }

    /// `MoreEdges()` (`cxx:255-258`).
    pub fn more_edges(&self) -> bool {
        self.iter_pos < self.iter_keys.len()
    }

    /// `NextEdge()` (`cxx:263-266`).
    pub fn next_edge(&mut self) {
        self.iter_pos += 1;
    }

    /// `RejectEdge(L, Par)` (`cxx:225-228`): no early rejection.
    pub fn reject_edge(&self, _l: &GpLin2d, _par: f64) -> bool {
        false
    }

    /// `CurrentEdge(E, Or)` (`cxx:233-237`).
    pub fn current_edge(&self) -> (Arc<dyn Curve2d>, Orientation) {
        let key = self.iter_keys[self.iter_pos];
        let item = self.find(key);
        (Arc::clone(item.curve()), item.orientation())
    }
}

/// `gp_Lin2d::SquareDistance(P)` (`gp_Lin2d.hxx:118-128`).
fn square_distance(l: &GpLin2d, p: &GpPnt2d) -> f64 {
    let d = l.distance(p);
    d * d
}
