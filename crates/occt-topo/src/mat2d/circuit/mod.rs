//! `MAT2d_Circuit` (`MAT2d_Circuit.hxx:35-113`, `MAT2d_Circuit.cxx:57-775`):
//! the closed circuit walked by the medial-axis algorithm, built from the
//! figure's lines plus the connexions of a [`Mat2dMiniPath`].
//!
//! `Perform` (`cxx:66-226`) decides which lines are open, closes each line on
//! itself (`InsertCorner`, `InitOpen`), doubles the open ones (`DoubleLine`),
//! runs the mini path and lays the items of the reference line plus the
//! portions between consecutive connexions into `geomElements`
//! (`ConstructCircuit`, `cxx:360-486`).
//!
//! Items are `MAT2dItem`s and the connexions are shared
//! `Rc<RefCell<Mat2dConnexion>>` handles, because `DoubleLine` rewrites the
//! connexion parameters of the mini path in place (`cxx:545-613`) and
//! `ConstructCircuit` then reads them back.

mod construct;
mod corners;
mod double_line;

use std::collections::HashMap;

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::kernel::geomabs::JoinType;
use occt_core::precision::CONFUSION;

use crate::mat2d::bi_int::Mat2dBiInt;
use crate::mat2d::connexion::Mat2dConnexion;
use crate::mat2d::item::{Mat2dFigure, Mat2dItem};
use crate::mat2d::mini_path::{ConnexionHandle, Mat2dMiniPath};

/// `class MAT2d_Circuit : public Standard_Transient` (`MAT2d_Circuit.hxx:39-115`).
#[derive(Debug, Clone)]
pub struct Mat2dCircuit {
    /// `direction` (`hxx:108`).
    direction: f64,
    /// `geomElements` (`hxx:109`).
    geom_elements: Vec<Mat2dItem>,
    /// `connexionMap` (`hxx:110`).
    connexion_map: HashMap<i32, ConnexionHandle>,
    /// `linkRefEqui` (`hxx:111`): `(line, curve) -> items of this circuit`.
    link_ref_equi: HashMap<Mat2dBiInt, Vec<i32>>,
    /// `linesLength` (`hxx:112`).
    lines_length: Vec<i32>,
    /// `myJoinType` (`hxx:113`).
    my_join_type: JoinType,
    /// `myIsOpenResult` (`hxx:114`).
    my_is_open_result: bool,
}

impl Default for Mat2dCircuit {
    /// `MAT2d_Circuit(GeomAbs_Arc, false)` (`cxx:57-63`).
    fn default() -> Self {
        Self::new(JoinType::Arc, false)
    }
}

impl Mat2dCircuit {
    /// `MAT2d_Circuit(aJoinType, IsOpenResult)` (`cxx:57-63`): `direction`
    /// starts at 0.
    pub fn new(a_join_type: JoinType, is_open_result: bool) -> Self {
        Self {
            direction: 0.0,
            geom_elements: Vec::new(),
            connexion_map: HashMap::new(),
            link_ref_equi: HashMap::new(),
            lines_length: Vec::new(),
            my_join_type: a_join_type,
            my_is_open_result: is_open_result,
        }
    }

    /// `Perform(FigItem, IsClosed, IndRefLine, Trigo)` (`cxx:66-226`).
    pub fn perform(
        &mut self,
        fig_item: &mut Mat2dFigure,
        is_closed: &[bool],
        ind_ref_line: i32,
        trigo: bool,
    ) {
        let nb_lines = fig_item.len() as i32;

        if trigo {
            self.direction = 1.0;
        } else {
            self.direction = -1.0;
        }

        // `cxx:90-95`: reinitialisation.
        self.geom_elements.clear();
        self.connexion_map.clear();
        self.link_ref_equi.clear();
        self.lines_length.clear();

        // `cxx:99-123`: detection of the open lines.
        let mut open = vec![false; nb_lines as usize + 1];
        for i in 1..=nb_lines {
            let line = &fig_item[(i - 1) as usize];
            let p1 = line.first().expect("MAT2d_Circuit: empty line").start_point();
            let p2 = line.last().expect("MAT2d_Circuit: empty line").end_point();
            if is_closed[(i - 1) as usize] {
                open[i as usize] = false;
            } else if p1.is_equal(&p2, CONFUSION) {
                open[i as usize] = false;
            } else {
                open[i as usize] = true;
            }
        }

        // `cxx:128-143`: sharp corners are inserted, or the ends of an open
        // line are added.
        for i in 1..=nb_lines {
            let line = &mut fig_item[(i - 1) as usize];
            if open[i as usize] {
                corners::init_open(line);
                self.lines_length.push(line.len() as i32);
            } else {
                self.insert_corner(line, self.direction);
                self.lines_length.push(line.len() as i32);
            }
        }

        // `cxx:146-157`: a single line needs no path.
        if nb_lines == 1 {
            if open[1] {
                let mut empty = Vec::new();
                self.double_line(&mut fig_item[0], &mut empty, None, self.direction);
                self.lines_length[0] = fig_item[0].len() as i32;
            }
            self.geom_elements = fig_item[0].clone();
            let length = self.geom_elements.len() as i32;
            self.update_link(1, 1, 1, length);
            let line_len = fig_item[0].len() as i32;
            self.lines_length.push(line_len);
            return;
        }

        // `cxx:164`: the mini path.
        let mut road = Mat2dMiniPath::new();
        road.perform(fig_item, ind_ref_line, trigo);

        // `cxx:172-198`: closure of the open lines.
        for i in 1..=nb_lines {
            if open[i as usize] {
                let cf = if road.is_root(i) {
                    None
                } else {
                    Some(road.father(i).clone())
                };
                if road.is_connexions_from(i) {
                    let mut from = std::mem::take(
                        road.connexions_from_mut(i)
                            .expect("MAT2d_Circuit: connexions from an open line"),
                    );
                    self.double_line(&mut fig_item[(i - 1) as usize], &mut from, cf, self.direction);
                    *road
                        .connexions_from_mut(i)
                        .expect("MAT2d_Circuit: connexions from an open line") = from;
                } else {
                    let mut empty = Vec::new();
                    self.double_line(
                        &mut fig_item[(i - 1) as usize],
                        &mut empty,
                        cf,
                        self.direction,
                    );
                }
                self.lines_length[(i - 1) as usize] = fig_item[(i - 1) as usize].len() as i32;
            }
        }

        // `cxx:199`: the path is rebuilt after the connexions were rewritten.
        road.run_on_connexions();

        // `cxx:216-217`.
        self.construct_circuit(fig_item, ind_ref_line, &road);
    }

    /// `PassByLast(C1, C2)` (`cxx:795-825`).
    pub fn pass_by_last(&self, c1: &Mat2dConnexion, c2: &Mat2dConnexion) -> bool {
        if c2.index_first_line() == c1.index_second_line() {
            if c2.index_item_on_first() < c1.index_item_on_second() {
                return true;
            } else if c2.index_item_on_first() == c1.index_item_on_second() {
                if c1.index_first_line() == c2.index_second_line() {
                    return true;
                }
                if c2.parameter_on_first() == c1.parameter_on_second() {
                    let v1 = vec_between(&c1.point_on_second(), &c1.point_on_first());
                    let v2 = vec_between(&c2.point_on_first(), &c2.point_on_second());
                    if v1.crossed(&v2) * self.direction > 0.0 {
                        return true;
                    }
                } else if c2.parameter_on_first() < c1.parameter_on_second() {
                    return true;
                }
            }
        }
        false
    }

    /// `Side(C1, Line)` (`cxx:774-791`).
    pub fn side(&self, c1: &Mat2dConnexion, line: &[Mat2dItem]) -> f64 {
        let v1 = vec_between(&c1.point_on_first(), &c1.point_on_second());
        let item = &line[(c1.index_item_on_first() - 1) as usize];
        let v2 = item
            .dn(c1.parameter_on_first(), 1)
            .expect("MAT2d_Circuit::Side: item is not a curve");
        if v1.crossed(&v2) > 0.0 {
            -1.0
        } else {
            1.0
        }
    }

    /// `UpDateLink(IFirst, ILine, ICurveFirst, ICurveLast)` (`cxx:829-854`).
    pub fn update_link(&mut self, first: i32, line: i32, curve_first: i32, curve_last: i32) {
        let mut i_equi = first;
        for i in curve_first..=curve_last {
            let key = Mat2dBiInt::new(line, i);
            self.link_ref_equi.entry(key).or_default().push(i_equi);
            i_equi += 1;
        }
    }

    /// `NumberOfItems()` (`cxx:707-710`).
    pub fn number_of_items(&self) -> i32 {
        self.geom_elements.len() as i32
    }

    /// `LineLength(IndexLine)` (`cxx:714-717`): `linesLength(I)`.
    pub fn line_length(&self, index_line: i32) -> i32 {
        self.lines_length[(index_line - 1) as usize]
    }

    /// `Value(Index)` (`cxx:721-724`).
    pub fn value(&self, index: i32) -> &Mat2dItem {
        &self.geom_elements[(index - 1) as usize]
    }

    /// `geomElements` (`hxx:109`).
    pub fn geom_elements(&self) -> &[Mat2dItem] {
        &self.geom_elements
    }

    /// `RefToEqui(IndLine, IndCurve)` (`cxx:728-732`):
    /// `linkRefEqui(MAT2d_BiInt(IndLine, IndCurve))`.
    pub fn ref_to_equi(&self, ind_line: i32, ind_curve: i32) -> &[i32] {
        self.link_ref_equi
            .get(&Mat2dBiInt::new(ind_line, ind_curve))
            .unwrap_or_else(|| {
                panic!(
                    "MAT2d_Circuit::RefToEqui: ({ind_line}, {ind_curve}) is not bound"
                )
            })
            .as_slice()
    }

    /// `Connexion(I)` (`cxx:759-762`): `connexionMap(I)`.
    pub fn connexion(&self, index: i32) -> &ConnexionHandle {
        self.connexion_map
            .get(&index)
            .unwrap_or_else(|| panic!("MAT2d_Circuit::Connexion: item {index} is not bound"))
    }

    /// `ConnexionOn(I)` (`cxx:766-769`).
    pub fn connexion_on(&self, index: i32) -> bool {
        self.connexion_map.contains_key(&index)
    }

    /// `myJoinType` (`hxx:113`).
    pub fn join_type(&self) -> JoinType {
        self.my_join_type
    }

    /// `myIsOpenResult` (`hxx:114`).
    pub fn is_open_result(&self) -> bool {
        self.my_is_open_result
    }

    /// `direction` (`hxx:108`).
    pub fn direction(&self) -> f64 {
        self.direction
    }

    /// `connexionMap.Bind(Index, C)` (`cxx:379`, `:417`, `:431`, `:442`).
    pub(crate) fn bind_connexion(&mut self, index: i32, c: ConnexionHandle) {
        self.connexion_map.insert(index, c);
    }

    /// `connexionMap.UnBind(Index)` (`cxx:446`).
    pub(crate) fn unbind_connexion(&mut self, index: i32) {
        self.connexion_map.remove(&index);
    }

    /// `linkRefEqui.ChangeFind(aBiInt)` (`cxx:739`).
    pub(crate) fn link_ref_equi_mut(&mut self, key: &Mat2dBiInt) -> &mut Vec<i32> {
        self.link_ref_equi.entry(*key).or_default()
    }

    /// `linkRefEqui` keys, the domain of the `Iterator` of `cxx:456-463`.
    pub(crate) fn link_ref_equi_keys(&self) -> Vec<Mat2dBiInt> {
        self.link_ref_equi.keys().copied().collect()
    }

    /// `linkRefEqui(Key).Length()` (`cxx:461`).
    pub(crate) fn link_ref_equi_len(&self, key: &Mat2dBiInt) -> usize {
        self.link_ref_equi.get(key).map(|v| v.len()).unwrap_or(0)
    }
}

/// `FigItem.ChangeValue(i)`: 1-based access to a line of the figure.
pub(crate) fn line_mut(fig_item: &mut Mat2dFigure, index: i32) -> &mut Vec<Mat2dItem> {
    &mut fig_item[(index - 1) as usize]
}

/// `gp_Vec2d(P1, P2)` (`gp_Vec2d.hxx`): the translation `P2 - P1`.
pub(crate) fn vec_between(p1: &GpPnt2d, p2: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y())
}
