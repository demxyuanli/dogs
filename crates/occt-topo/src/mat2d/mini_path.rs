//! `MAT2d_MiniPath` (`MAT2d_MiniPath.hxx:35-115`, `MAT2d_MiniPath.cxx:19-435`):
//! the shortest tree of `MAT2d_Connexion` that links every line of a
//! `MAT2d` figure.
//!
//! `Perform` (`cxx:45-131`) computes the minimum-distance connexion between
//! every pair of lines with `MinimumL1L2` (`cxx:294-435`), grows the tree from
//! line `IndStart` by repeatedly taking the shortest connexion between the
//! reached set and the rest, and stores it through `Append` (`cxx:143-181`),
//! which keeps the connexions that start on a line ordered by `IsAfter`.
//! `RunOnConnexions` (`cxx:228-242`) / `ExploSons` (`cxx:247-286`) then flatten
//! that tree into `Path` (`cxx:186-190`), the sequence the circuit walks.
//!
//! OCCT shares the connexions as `occ::handle`s: `MAT2d_Circuit::DoubleLine`
//! rewrites `IndexItemOnFirst` / `ParameterOnFirst` of the connexions it is
//! handed, and `ConstructCircuit` reads those same objects back out of `Path`.
//! The port keeps that aliasing with `Rc<RefCell<...>>` handles.

use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::extrema2d::{ExtremaExtCC2d, ExtremaExtPC2d};

use super::connexion::Mat2dConnexion;
use super::item::{Mat2dFigure, Mat2dItem};

/// `occ::handle(MAT2d_Connexion)`.
pub type ConnexionHandle = Rc<RefCell<Mat2dConnexion>>;

/// `class MAT2d_MiniPath` (`MAT2d_MiniPath.hxx:45-116`).
#[derive(Debug, Clone)]
pub struct Mat2dMiniPath {
    /// `theConnexions` (`hxx:110`): the connexions issued from a line, ordered
    /// by `MAT2d_Connexion::IsAfter`.
    the_connexions: HashMap<i32, Vec<ConnexionHandle>>,
    /// `theFather` (`hxx:111`): the connexion that ends on a line.
    the_father: HashMap<i32, ConnexionHandle>,
    /// `thePath` (`hxx:112`).
    the_path: Vec<ConnexionHandle>,
    /// `theDirection` (`hxx:113`).
    the_direction: f64,
    /// `indStart` (`hxx:114`).
    ind_start: i32,
}

impl Default for Mat2dMiniPath {
    /// `MAT2d_MiniPath()` (`cxx:31-35`).
    fn default() -> Self {
        Self::new()
    }
}

impl Mat2dMiniPath {
    /// `MAT2d_MiniPath()` (`cxx:31-35`): `theDirection = 1.0`, `indStart = 0`.
    pub fn new() -> Self {
        Self {
            the_connexions: HashMap::new(),
            the_father: HashMap::new(),
            the_path: Vec::new(),
            the_direction: 1.0,
            ind_start: 0,
        }
    }

    /// `Perform(Figure, IndStart, Sense)` (`cxx:45-131`).
    pub fn perform(&mut self, figure: &Mat2dFigure, ind_start: i32, sense: bool) {
        let nb_lines = figure.len() as i32;
        if std::env::var("MAT2D_DEBUG").is_ok() {
            let lens: Vec<usize> = figure.iter().map(|l| l.len()).collect();
            eprintln!(
                "MINIPATH lines={} lens={:?} ind_start={} sense={}",
                nb_lines, lens, ind_start, sense
            );
        }
        // `cxx:51`: `NCollection_Array2<handle> Connexion(1, NbLines, 1, NbLines)`.
        let mut connexion: Vec<Vec<Option<ConnexionHandle>>> =
            vec![vec![None; nb_lines as usize + 1]; nb_lines as usize + 1];

        self.ind_start = ind_start;
        self.the_direction = 1.0;
        if sense {
            self.the_direction = -1.0;
        }

        // `cxx:65-73`: the minimum-distance connexion of every pair of lines,
        // plus its reverse.
        for i in 1..nb_lines {
            for j in i + 1..=nb_lines {
                let c = Rc::new(RefCell::new(self.minimum_l1_l2(figure, i, j)));
                connexion[j as usize][i as usize] =
                    Some(Rc::new(RefCell::new(c.borrow().reversed())));
                connexion[i as usize][j as usize] = Some(c);
            }
        }

        let mut set1: Vec<i32> = Vec::new();
        let mut set2: Vec<i32> = Vec::new();
        let mut isuiv: i32 = 0;
        let mut min_on_set1: i32 = 0;
        let mut min_on_set2: i32 = 0;

        // `cxx:86-93`: `Set1` starts with the reference line, `Set2` holds the
        // others.
        set1.push(ind_start);
        for i in 1..=nb_lines {
            if i != ind_start {
                set2.push(i);
            }
        }

        // `cxx:105-125`: add the shortest connexion between `Set1` and `Set2`,
        // move its line of `Set2` into `Set1`, repeat.
        while !set2.is_empty() {
            let mut dist_s1_s2 = f64::MAX;
            for &indice_line1 in &set1 {
                for (j, &indice_line2) in set2.iter().enumerate() {
                    let c = connexion[indice_line1 as usize][indice_line2 as usize]
                        .as_ref()
                        .expect("MAT2d_MiniPath::Perform: connexion not computed");
                    let distance = c.borrow().distance();
                    if distance < dist_s1_s2 {
                        isuiv = j as i32 + 1;
                        dist_s1_s2 = distance;
                        min_on_set1 = indice_line1;
                        min_on_set2 = indice_line2;
                    }
                }
            }
            // `cxx:121-123`.
            set1.push(set2[(isuiv - 1) as usize]);
            set2.remove((isuiv - 1) as usize);
            let c = connexion[min_on_set1 as usize][min_on_set2 as usize]
                .as_ref()
                .expect("MAT2d_MiniPath::Perform: connexion not computed")
                .clone();
            self.append(&c);
        }

        self.run_on_connexions();
    }

    /// `Append(C)` (`cxx:143-181`): insert `C` among the connexions issued from
    /// its first line, ordered by `IsAfter(C, theDirection)`, and bind it as the
    /// father of its second line.
    fn append(&mut self, c: &ConnexionHandle) {
        let first_line = c.borrow().index_first_line();
        let second_line = c.borrow().index_second_line();

        let Some(seq) = self.the_connexions.get_mut(&first_line) else {
            // `cxx:147-152`.
            self.the_connexions.insert(first_line, vec![c.clone()]);
            self.the_father.insert(second_line, c.clone());
            return;
        };

        // `cxx:155-167`: the insertion index, 1-based in OCCT.
        let mut index_after = 0usize;
        let direction = self.the_direction;
        {
            let c_ref = c.borrow();
            for (i, cc) in seq.iter().enumerate() {
                if cc.borrow().is_after(&c_ref, direction) {
                    index_after = i + 1;
                    break;
                }
            }
        }
        if index_after == 0 {
            // `cxx:174-177`.
            seq.push(c.clone());
        } else {
            // `cxx:179`: `InsertBefore`.
            seq.insert(index_after - 1, c.clone());
        }
        self.the_father.insert(second_line, c.clone());
    }

    /// `Path()` (`cxx:186-190`).
    pub fn path(&self) -> &[ConnexionHandle] {
        &self.the_path
    }

    /// `IsConnexionsFrom(i)` (`cxx:193-196`).
    pub fn is_connexions_from(&self, index: i32) -> bool {
        self.the_connexions.contains_key(&index)
    }

    /// `ConnexionsFrom(i)` (`cxx:204-207`): `theConnexions.ChangeFind(i)`.
    pub fn connexions_from(&self, index: i32) -> Option<&Vec<ConnexionHandle>> {
        self.the_connexions.get(&index)
    }

    /// `ConnexionsFrom(i)` (`cxx:204-207`) as the mutable sequence
    /// `MAT2d_Circuit::DoubleLine` re-orders (`cxx:511-513`).
    pub fn connexions_from_mut(&mut self, index: i32) -> Option<&mut Vec<ConnexionHandle>> {
        self.the_connexions.get_mut(&index)
    }

    /// `IsRoot(ILine)` (`cxx:211-214`).
    pub fn is_root(&self, line: i32) -> bool {
        line == self.ind_start
    }

    /// `Father(ILine)` (`cxx:219-222`): `theFather.ChangeFind(ILine)`.
    pub fn father(&self, line: i32) -> &ConnexionHandle {
        self.the_father
            .get(&line)
            .unwrap_or_else(|| panic!("MAT2d_MiniPath::Father: line {line} is not bound"))
    }

    /// `theDirection` (`hxx:113`), read by `MAT2d_Circuit`.
    pub fn direction(&self) -> f64 {
        self.the_direction
    }

    /// `RunOnConnexions()` (`cxx:228-242`): rebuild `thePath` from the current
    /// `theConnexions`.
    pub fn run_on_connexions(&mut self) {
        let Some(sc) = self.the_connexions.get(&self.ind_start).cloned() else {
            panic!("MAT2d_MiniPath::RunOnConnexions: root has no connexion");
        };

        self.the_path.clear();
        for c in sc {
            self.the_path.push(c.clone());
            self.explo_sons(&c);
            self.the_path
                .push(Rc::new(RefCell::new(c.borrow().reversed())));
        }
    }

    /// `ExploSons(CResult, CRef)` (`cxx:247-286`).
    fn explo_sons(&mut self, c_ref: &ConnexionHandle) {
        let index = c_ref.borrow().index_second_line();

        let Some(sc) = self.the_connexions.get(&index).cloned() else {
            return;
        };
        let crr = Rc::new(RefCell::new(c_ref.borrow().reversed()));
        let direction = self.the_direction;

        // `cxx:262-272`: the children after `CRef->Reverse()`.
        for c in &sc {
            if c.borrow().is_after(&crr.borrow(), direction) {
                self.the_path.push(c.clone());
                self.explo_sons(c);
                self.the_path
                    .push(Rc::new(RefCell::new(c.borrow().reversed())));
            }
        }

        // `cxx:274-287`: the children up to `CRef->Reverse()`.
        for c in &sc {
            if !c.borrow().is_after(&crr.borrow(), direction) {
                self.the_path.push(c.clone());
                self.explo_sons(c);
                self.the_path
                    .push(Rc::new(RefCell::new(c.borrow().reversed())));
            } else {
                break;
            }
        }
    }

    /// `MinimumL1L2(Figure, IL1, IL2)` (`cxx:294-435`): the connexion of
    /// minimum distance between the items of line `IL1` and those of line
    /// `IL2`, oriented from `IL1` to `IL2`.
    fn minimum_l1_l2(&self, figure: &Mat2dFigure, il1: i32, il2: i32) -> Mat2dConnexion {
        let l1 = &figure[(il1 - 1) as usize];
        let l2 = &figure[(il2 - 1) as usize];

        let mut dist_l1_l2_2 = f64::MAX;
        let mut imin_c1 = 0;
        let mut imin_c2 = 0;
        let mut parameter_on_c1 = 0.0;
        let mut parameter_on_c2 = 0.0;
        let mut point1 = GpPnt2d::zero();
        let mut point2 = GpPnt2d::zero();

        for (ic1, item1) in l1.iter().enumerate() {
            for (ic2, item2) in l2.iter().enumerate() {
                match (item1, item2) {
                    (Mat2dItem::Point(p1), Mat2dItem::Point(p2)) => {
                        // `cxx:345-357`.
                        let dist_p1_p2_2 = p1.square_distance(p2);
                        if dist_p1_p2_2 <= dist_l1_l2_2 {
                            dist_l1_l2_2 = dist_p1_p2_2;
                            imin_c1 = ic1 as i32 + 1;
                            imin_c2 = ic2 as i32 + 1;
                            point1 = *p1;
                            point2 = *p2;
                            parameter_on_c1 = 0.0;
                            parameter_on_c2 = 0.0;
                        }
                    }
                    (Mat2dItem::Point(p1), Mat2dItem::Curve(_)) => {
                        // `cxx:358-379`.
                        let c2 = item2.expect_curve();
                        let extremas = ExtremaExtPC2d::new(c2, p1);
                        if extremas.is_done() {
                            for i in 1..=extremas.nb_ext() {
                                if extremas.square_distance(i) < dist_l1_l2_2 {
                                    dist_l1_l2_2 = extremas.square_distance(i);
                                    imin_c1 = ic1 as i32 + 1;
                                    imin_c2 = ic2 as i32 + 1;
                                    let (u2, p2) = extremas.point(i);
                                    parameter_on_c1 = 0.0;
                                    parameter_on_c2 = u2;
                                    point1 = *p1;
                                    point2 = p2;
                                }
                            }
                        }
                    }
                    (Mat2dItem::Curve(_), Mat2dItem::Point(p2)) => {
                        // `cxx:380-402`.
                        let c1 = item1.expect_curve();
                        let extremas = ExtremaExtPC2d::new(c1, p2);
                        if extremas.is_done() {
                            for i in 1..=extremas.nb_ext() {
                                if extremas.square_distance(i) < dist_l1_l2_2 {
                                    dist_l1_l2_2 = extremas.square_distance(i);
                                    imin_c1 = ic1 as i32 + 1;
                                    imin_c2 = ic2 as i32 + 1;
                                    let (u1, p1) = extremas.point(i);
                                    parameter_on_c2 = 0.0;
                                    parameter_on_c1 = u1;
                                    point1 = p1;
                                    point2 = *p2;
                                }
                            }
                        }
                    }
                    (Mat2dItem::Curve(_), Mat2dItem::Curve(_)) => {
                        // `cxx:403-424`.
                        let c1 = item1.expect_curve();
                        let c2 = item2.expect_curve();
                        let extremas = ExtremaExtCC2d::new(c1, c2);
                        if std::env::var("MAT2D_DEBUG").is_ok() {
                            let a = c1.d0(c1.first_parameter());
                            let b = c1.d0(c1.last_parameter());
                            let d = c2.d0(c2.first_parameter());
                            let e = c2.d0(c2.last_parameter());
                            let l1 = occt_geom2d::adaptor_curve::basis2d(c1).gp_lin2d();
                            let l2 = occt_geom2d::adaptor_curve::basis2d(c2).gp_lin2d();
                            eprintln!(
                                "CCPAIR {}:{} types {:?}/{:?} done={} par={} nbext={}",
                                ic1 + 1,
                                ic2 + 1,
                                occt_geom2d::adaptor_curve::type2d(c1),
                                occt_geom2d::adaptor_curve::type2d(c2),
                                extremas.is_done(),
                                extremas.is_parallel(),
                                extremas.nb_ext()
                            );
                            eprintln!(
                                "  C1 trim=[{:.3},{:.3}] seg ({:.3},{:.3})->({:.3},{:.3}) line loc=({:.3},{:.3}) dir=({:.3},{:.3})",
                                c1.first_parameter(), c1.last_parameter(),
                                a.x(), a.y(), b.x(), b.y(),
                                l1.map(|l| l.pos.loc.x()).unwrap_or(f64::NAN),
                                l1.map(|l| l.pos.loc.y()).unwrap_or(f64::NAN),
                                l1.map(|l| l.pos.vdir.x).unwrap_or(f64::NAN),
                                l1.map(|l| l.pos.vdir.y).unwrap_or(f64::NAN),
                            );
                            eprintln!(
                                "  C2 trim=[{:.3},{:.3}] seg ({:.3},{:.3})->({:.3},{:.3}) line loc=({:.3},{:.3}) dir=({:.3},{:.3})",
                                c2.first_parameter(), c2.last_parameter(),
                                d.x(), d.y(), e.x(), e.y(),
                                l2.map(|l| l.pos.loc.x()).unwrap_or(f64::NAN),
                                l2.map(|l| l.pos.loc.y()).unwrap_or(f64::NAN),
                                l2.map(|l| l.pos.vdir.x).unwrap_or(f64::NAN),
                                l2.map(|l| l.pos.vdir.y).unwrap_or(f64::NAN),
                            );
                        }
                        if !extremas.is_parallel() && extremas.is_done() {
                            for i in 1..=extremas.nb_ext() {
                                if extremas.square_distance(i) < dist_l1_l2_2 {
                                    dist_l1_l2_2 = extremas.square_distance(i);
                                    imin_c1 = ic1 as i32 + 1;
                                    imin_c2 = ic2 as i32 + 1;
                                    let (u1, p1, u2, p2) = extremas.points(i);
                                    parameter_on_c1 = u1;
                                    parameter_on_c2 = u2;
                                    point1 = p1;
                                    point2 = p2;
                                }
                            }
                        }
                    }
                }
            }
        }

        // `cxx:427-434`.
        if std::env::var("MAT2D_DEBUG").is_ok() && (imin_c1 == 0 || imin_c2 == 0) {
            eprintln!(
                "MINL1L2 ZERO il1={} il2={} len1={} len2={} d2={:?}",
                il1,
                il2,
                l1.len(),
                l2.len(),
                dist_l1_l2_2
            );
            for (k, it) in l1.iter().enumerate() {
                eprintln!("  L1[{}] {:?}", k + 1, item_kind(it));
            }
            for (k, it) in l2.iter().enumerate() {
                eprintln!("  L2[{}] {:?}", k + 1, item_kind(it));
            }
        }
        Mat2dConnexion::with_values(
            il1,
            il2,
            imin_c1,
            imin_c2,
            dist_l1_l2_2.sqrt(),
            parameter_on_c1,
            parameter_on_c2,
            point1,
            point2,
        )
    }
}

/// `gp_Vec2d(P1, P2)` (`gp_Vec2d.hxx`): the translation `P2 - P1`.
#[allow(dead_code)]
fn vec_between(p1: &GpPnt2d, p2: &GpPnt2d) -> GpVec2d {
    GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y())
}

/// Debug helper: `P` for a point item, `C` for a curve item.
#[allow(dead_code)]
fn item_kind(it: &Mat2dItem) -> &'static str {
    match it {
        Mat2dItem::Point(_) => "P",
        Mat2dItem::Curve(_) => "C",
    }
}
