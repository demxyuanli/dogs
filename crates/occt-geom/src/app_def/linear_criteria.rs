//! `AppDef_LinearCriteria` — the smoothing functional
//! `W_quad J_LS + W_qual (p1 J1 + p2 J2 + p3 J3)` and the assembly data the
//! variational optimizer drives.
//!
//! Source: `ModelingData/TKGeomBase/AppDef/AppDef_LinearCriteria.cxx:42-852` and
//! `AppDef_LinearCriteria.hxx` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `myCriteria[3]` become a `[Box<dyn Criterion>; 3]`; `Set(Coeff)` clones the
//!   coefficient array for each criterion because the ported criteria own their
//!   `myCoeff`.
//! - `myPntWeight` and `myParameters` are `HArray1<double>` with OCCT bounds;
//!   they are [`MathVector`] here, so `Lower()` / `Value(i)` keep their
//!   semantics.
//! - `myCache` is a 1-based `HArray1<double>` in OCCT; the port keeps a 0-based
//!   `Vec<f64>` and subtracts one at every `Value(ii)`.
//! - `AssemblyTable()` returns the `[dimension][element]` global index lists
//!   directly, i.e. the `NCollection_HArray2` content `FEmTool_Assembly`
//!   consumes.

use std::cell::RefCell;
use std::rc::Rc;

use occt_core::kernel::geomabs::Shape;

use occt_math::fem_tool::array2::Array2;
use occt_math::fem_tool::criteria::{Criterion, LinearFlexion, LinearJerk, LinearTension};
use occt_math::fem_tool::curve::Curve;
use occt_math::matrix::MathMatrix;
use occt_math::vector::MathVector;

use crate::app_def::multi_line::MultiLine;
use crate::app_def::my_line_tool::MyLineTool;
use crate::app_def::smooth_criterion::SmoothCriterion;

/// `order(B)` (`cxx:38-41`): `PLib_HermitJacobi::NivConstr()`.
fn order(base: &occt_math::plib_hermit_jacobi::HermitJacobi) -> i32 {
    base.niv_constr()
}

impl SmoothCriterion for LinearCriteria {
    fn set_parameters(&mut self, parameters: &MathVector) {
        LinearCriteria::set_parameters(self, parameters);
    }

    fn set_curve(&mut self, c: &Rc<RefCell<Curve>>) {
        LinearCriteria::set_curve(self, c);
    }

    fn get_curve(&self) -> Option<Rc<RefCell<Curve>>> {
        LinearCriteria::get_curve(self)
    }

    fn set_estimation(&mut self, e1: f64, e2: f64, e3: f64) {
        LinearCriteria::set_estimation(self, e1, e2, e3);
    }

    fn est_length(&mut self) -> &mut f64 {
        LinearCriteria::est_length(self)
    }

    fn estimation(&self) -> (f64, f64, f64) {
        let e = LinearCriteria::estimation(self);
        (e[0], e[1], e[2])
    }

    fn assembly_table(&self) -> Vec<Vec<Vec<i32>>> {
        LinearCriteria::assembly_table(self)
    }

    fn dependence_table(&self) -> Vec<Vec<i32>> {
        LinearCriteria::dependence_table(self)
    }

    fn quality_values(
        &mut self,
        j1min: f64,
        j2min: f64,
        j3min: f64,
        j1: &mut f64,
        j2: &mut f64,
        j3: &mut f64,
    ) -> i32 {
        LinearCriteria::quality_values(self, j1min, j2min, j3min, j1, j2, j3)
    }

    fn error_values(
        &mut self,
        max_error: &mut f64,
        quadratic_error: &mut f64,
        average_error: &mut f64,
    ) {
        LinearCriteria::error_values(self, max_error, quadratic_error, average_error);
    }

    fn hessian(&mut self, element: usize, dim1: i32, dim2: i32, h: &mut MathMatrix) {
        LinearCriteria::hessian(self, element, dim1, dim2, h);
    }

    fn gradient(&mut self, element: usize, dim: i32, g: &mut MathVector) {
        LinearCriteria::gradient(self, element, dim, g);
    }

    fn input_vector(&mut self, x: &MathVector, ass_table: &[Vec<Vec<i32>>]) {
        LinearCriteria::input_vector(self, x, ass_table);
    }

    fn set_weight(&mut self, quadratic_weight: f64, quality_weight: f64, percent: [f64; 3]) {
        LinearCriteria::set_weight(self, quadratic_weight, quality_weight, percent);
    }

    fn weight(&self) -> (f64, f64) {
        LinearCriteria::weight(self)
    }

    fn set_pnt_weight(&mut self, weight: &MathVector) {
        LinearCriteria::set_pnt_weight(self, weight);
    }
}

/// `GeomAbs_Shape` from the constraint order (`cxx:84-95`).
fn constraint_order(order: i32) -> Shape {
    match order {
        0 => Shape::C0,
        1 => Shape::C1,
        2 => Shape::C2,
        _ => Shape::C0,
    }
}

/// `AppDef_LinearCriteria` (`AppDef_LinearCriteria.hxx`).
pub struct LinearCriteria {
    ssp: MultiLine,
    quadratic_weight: f64,
    quality_weight: f64,
    pnt_weight: MathVector,
    length: f64,
    estimation: [f64; 3],
    percent: [f64; 3],
    criteria: Vec<Box<dyn Criterion>>,
    curve: Option<Rc<RefCell<Curve>>>,
    parameters: Option<MathVector>,
    /// `myE` — element the cache was built for.
    e: i32,
    /// `myCache` (`cxx:839`).
    cache: Vec<f64>,
    /// `IF` / `IL` (`cxx:817-846`).
    if_: i32,
    il: i32,
}

impl LinearCriteria {
    /// `AppDef_LinearCriteria(SSP, FirstPoint, LastPoint)` (`cxx:46-60`).
    pub fn new(ssp: &MultiLine, first_point: i32, last_point: i32) -> Self {
        let mut pnt_weight = MathVector::new(first_point as usize, last_point as usize);
        pnt_weight.init(1.0);
        Self {
            ssp: ssp.clone(),
            quadratic_weight: 0.0,
            quality_weight: 0.0,
            pnt_weight,
            length: 0.0,
            estimation: [0.0; 3],
            percent: [0.0; 3],
            criteria: Vec::new(),
            curve: None,
            parameters: None,
            e: 0,
            cache: Vec::new(),
            if_: 0,
            il: 0,
        }
    }

    /// `SetParameters(Parameters)` (`cxx:64-68`).
    pub fn set_parameters(&mut self, parameters: &MathVector) {
        self.parameters = Some(parameters.clone());
        self.e = 0; // Cache becomes invalid.
    }

    /// `SetCurve(C)` (`cxx:73-156`).
    ///
    /// `myCurve` is a `Handle(FEmTool_Curve)`, so the criteria and the caller
    /// share one mutable curve; the port keeps that sharing with
    /// `Rc<RefCell<Curve>>` and uses `Rc::ptr_eq` for `myCurve != C`
    /// (`cxx:108`).
    pub fn set_curve(&mut self, c: &Rc<RefCell<Curve>>) {
        match self.curve.as_ref() {
            None => {
                let (mx_deg, nb_dim, ord) = {
                    let cu = c.borrow();
                    (cu.base().work_degree(), cu.dimension() as i32, order(cu.base()))
                };
                self.rebuild_criteria(mx_deg, nb_dim, ord);
                self.curve = Some(Rc::clone(c));
            }
            Some(old) if !Rc::ptr_eq(old, c) => {
                let (old_mx_deg, old_nb_dim, old_order) = {
                    let cu = old.borrow();
                    (cu.base().work_degree(), cu.dimension() as i32, order(cu.base()))
                };

                self.curve = Some(Rc::clone(c));

                let (mx_deg, nb_dim, ord) = {
                    let cu = c.borrow();
                    (cu.base().work_degree(), cu.dimension() as i32, order(cu.base()))
                };

                if mx_deg != old_mx_deg || ord != old_order {
                    self.rebuild_criteria(mx_deg, nb_dim, ord);
                } else if nb_dim != old_nb_dim {
                    let coeff = Array2::new(0, 0, 1, nb_dim, 0.0);
                    for crit in self.criteria.iter_mut() {
                        crit.set_coeff(coeff.clone());
                    }
                }
            }
            Some(_) => {}
        }
    }

    /// `cxx:96-106` / `cxx:131-154`.
    fn rebuild_criteria(&mut self, mx_deg: i32, nb_dim: i32, ord: i32) {
        let co = constraint_order(ord);
        self.criteria = vec![
            Box::new(LinearTension::new(mx_deg, co).expect("LinearTension")),
            Box::new(LinearFlexion::new(mx_deg, co).expect("LinearFlexion")),
            Box::new(LinearJerk::new(mx_deg, co).expect("LinearJerk")),
        ];
        let coeff = Array2::new(0, 0, 1, nb_dim, 0.0);
        for crit in self.criteria.iter_mut() {
            crit.set_coeff(coeff.clone());
        }
    }

    /// `GetCurve(C)` (`cxx:160-163`).
    pub fn get_curve(&self) -> Option<Rc<RefCell<Curve>>> {
        self.curve.as_ref().map(Rc::clone)
    }

    /// `SetEstimation(E1, E2, E3)` (`cxx:166-173`).
    pub fn set_estimation(&mut self, e1: f64, e2: f64, e3: f64) {
        self.estimation = [e1, e2, e3];
    }

    /// `EstLength()` (`cxx:175-178`).
    pub fn est_length(&mut self) -> &mut f64 {
        &mut self.length
    }

    /// `GetEstimation(E1, E2, E3)` (`cxx:181-188`).
    pub fn estimation(&self) -> [f64; 3] {
        self.estimation
    }

    /// `AssemblyTable()` (`cxx:190-273`): `[dimension][element]` global index
    /// lists, indexed in element order (OCCT arrays are 1-based).
    pub fn assembly_table(&self) -> Vec<Vec<Vec<i32>>> {
        let curve = self.curve.as_ref().expect("AppDef_LinearCriteria::AssemblyTable").clone();
        let cu = curve.borrow();
        let nb_dim = cu.dimension();
        let nb_elm = cu.nb_elements();
        let nc1 = order(cu.base()) + 1;
        let mx_deg = cu.base().work_degree();

        let mut ass_table = vec![vec![Vec::new(); nb_elm]; nb_dim];
        let mut nb_glob_var = 0i32;
        let mut gi0;

        // dim = 1, el = 1 (`cxx:210-229`).
        let mut glob_index = vec![0i32; (mx_deg + 1) as usize];
        for i in 0..nc1 {
            nb_glob_var += 1;
            glob_index[i as usize] = nb_glob_var;
        }
        gi0 = mx_deg - 2 * nc1 + 1;
        for i in nc1..(2 * nc1) {
            nb_glob_var += 1;
            glob_index[i as usize] = nb_glob_var + gi0;
        }
        for i in (2 * nc1)..=mx_deg {
            nb_glob_var += 1;
            glob_index[i as usize] = nb_glob_var - nc1;
        }
        gi0 = nb_glob_var - nc1 + 1;
        ass_table[0][0] = glob_index;

        // dim = 1, el = 2..NbElm (`cxx:232-252`).
        for el in 2..=nb_elm {
            let mut glob_index = vec![0i32; (mx_deg + 1) as usize];
            for i in 0..nc1 {
                glob_index[i as usize] = gi0 + i;
            }
            gi0 = mx_deg - 2 * nc1 + 1;
            for i in nc1..(2 * nc1) {
                nb_glob_var += 1;
                glob_index[i as usize] = nb_glob_var + gi0;
            }
            for i in (2 * nc1)..=mx_deg {
                nb_glob_var += 1;
                glob_index[i as usize] = nb_glob_var - nc1;
            }
            gi0 = nb_glob_var - nc1 + 1;
            ass_table[0][el - 1] = glob_index;
        }

        // Other dimensions (`cxx:256-270`).
        gi0 = nb_glob_var;
        for dim in 2..=nb_dim {
            for el in 1..=nb_elm {
                let aux = ass_table[0][el - 1].clone();
                let glob_index: Vec<i32> =
                    aux.iter().map(|v| v + nb_glob_var).collect();
                ass_table[dim - 1][el - 1] = glob_index;
            }
            nb_glob_var += gi0;
        }

        ass_table
    }

    /// `DependenceTable()` (`cxx:277-294`).
    pub fn dependence_table(&self) -> Vec<Vec<i32>> {
        let curve = self.curve.as_ref().expect("dependence_table").clone();
        let dim = curve.borrow().dimension();
        let mut dep = vec![vec![0i32; dim]; dim];
        for i in 0..dim {
            dep[i][i] = 1;
        }
        dep
    }

    /// `QualityValues(J1min, J2min, J3min, J1, J2, J3)` (`cxx:298-440`).
    pub fn quality_values(
        &mut self,
        j1min: f64,
        j2min: f64,
        j3min: f64,
        j1: &mut f64,
        j2: &mut f64,
        j3: &mut f64,
    ) -> i32 {
        let (nb_dim, nb_elm) = {
            let curve = self.curve.as_ref().expect("quality_values").clone();
            let cu = curve.borrow();
            (cu.dimension() as i32, cu.nb_elements())
        };

        let mut deg = 0i32;
        let mut coeff = Array2::new(0, 0, 1, nb_dim, 0.0);

        *j1 = 0.0;
        *j2 = 0.0;
        *j3 = 0.0;
        for el in 1..=nb_elm {
            let curve = self.curve.as_ref().unwrap().clone();
            let cu = curve.borrow();
            let curdeg = cu.degree(el);
            if deg != curdeg {
                deg = curdeg;
                coeff = Array2::new(0, deg, 1, nb_dim, 0.0);
            }

            cu.get_element(el, &mut coeff);

            let knots = cu.knots();
            let (u_first, u_last) = (knots[el], knots[el + 1]);
            drop(cu);

            self.criteria[0].set_coeff(coeff.clone());
            self.criteria[0].set_interval(u_first, u_last);
            *j1 += self.criteria[0].value();

            self.criteria[1].set_coeff(coeff.clone());
            self.criteria[1].set_interval(u_first, u_last);
            *j2 += self.criteria[1].value();

            self.criteria[2].set_coeff(coeff.clone());
            self.criteria[2].set_interval(u_first, u_last);
            *j3 += self.criteria[2].value();
        }

        let jes_min = [j1min, j2min, j3min];
        let val_cri = [*j1, *j2, *j3];

        let mut icdana = 0;

        // (2) Test the improvement of estimates (`cxx:363-378`).
        for i in 0..=2 {
            if val_cri[i] < 0.8 * self.estimation[i] && self.estimation[i] > jes_min[i] {
                if icdana < 1 {
                    icdana = 1;
                }
                if val_cri[i] < 0.1 * self.estimation[i] {
                    icdana = 2;
                }
                self.estimation[i] = (1.05 * val_cri[i]).max(jes_min[i]);
            }
        }

        // (3) Update the estimates (`cxx:381-438`).
        let upd = |icdana: &mut i32, val: f64, est: &mut f64, factor: f64, poor: f64| {
            if val > *est * factor {
                *est += val * poor;
                if *icdana == 0 {
                    if val > *est * 10.0 {
                        *icdana = 2;
                    } else {
                        *icdana = 1;
                    }
                } else {
                    *icdana = 2;
                }
            }
        };
        upd(&mut icdana, val_cri[0], &mut self.estimation[0], 2.0, 0.1);
        {
            // `cxx:400-418` — the threshold is 20 and the poor-conditioning
            // ratio 100 for J2 and J3.
            if val_cri[1] > self.estimation[1] * 20.0 {
                self.estimation[1] += val_cri[1] * 0.1;
                if icdana == 0 {
                    if val_cri[1] > self.estimation[1] * 100.0 {
                        icdana = 2;
                    } else {
                        icdana = 1;
                    }
                } else {
                    icdana = 2;
                }
            }
            if val_cri[2] > self.estimation[2] * 20.0 {
                self.estimation[2] += val_cri[2] * 0.05;
                if icdana == 0 {
                    if val_cri[2] > self.estimation[2] * 100.0 {
                        icdana = 2;
                    } else {
                        icdana = 1;
                    }
                } else {
                    icdana = 2;
                }
            }
        }

        icdana
    }

    /// `ErrorValues(MaxError, QuadraticError, AverageError)` (`cxx:443-510`).
    pub fn error_values(
        &mut self,
        max_error: &mut f64,
        quadratic_error: &mut f64,
        average_error: &mut f64,
    ) {
        let nb_dim = {
            let curve = self.curve.as_ref().expect("error_values").clone();
            let d = curve.borrow().dimension() as i32;
            d
        };
        let nb_p2d = MyLineTool::nb_p2d(&self.ssp);
        let nb_p3d = MyLineTool::nb_p3d(&self.ssp);

        if nb_dim != 2 * nb_p2d + 3 * nb_p3d {
            return; // `Standard_DomainError` (`cxx:457-461`).
        }

        let params = self.parameters.clone().expect("error_values parameters");
        let mut base_point = vec![0.0f64; nb_dim as usize];

        *max_error = 0.0;
        *quadratic_error = 0.0;
        *average_error = 0.0;

        for i in params.lower()..=params.upper() {
            {
                let curve = self.curve.as_ref().unwrap().clone();
                let mut cu = curve.borrow_mut();
                cu.d0(params.value(i), &mut base_point);
            }

            let mut c0 = 0usize;
            let (tab_p3d, tab_p2d) = MyLineTool::value(&self.ssp, i as i32);
            for ipnt in 0..nb_p3d as usize {
                let p3d = &tab_p3d[ipnt];
                let dx = base_point[c0] - p3d.x();
                let dy = base_point[c0 + 1] - p3d.y();
                let dz = base_point[c0 + 2] - p3d.z();
                let sqr = dx * dx + dy * dy + dz * dz;
                let dist = sqr.sqrt();
                *max_error = max_error.max(dist);
                *quadratic_error += sqr;
                *average_error += dist;
                c0 += 3;
            }

            for ipnt in 0..nb_p2d as usize {
                let p2d = &tab_p2d[ipnt];
                let dx = base_point[c0] - p2d.x();
                let dy = base_point[c0 + 1] - p2d.y();
                let sqr = dx * dx + dy * dy;
                let dist = sqr.sqrt();
                *max_error = max_error.max(dist);
                *quadratic_error += sqr;
                *average_error += dist;
                c0 += 2;
            }
        }
    }

    /// `Hessian(Element, Dimension1, Dimension2, H)` (`cxx:514-619`).
    pub fn hessian(&mut self, element: usize, dim1: i32, dim2: i32, h: &mut MathMatrix) {
        let (mx_deg, ord) = {
            let curve = self.curve.as_ref().expect("hessian").clone();
            let cu = curve.borrow();
            (cu.base().work_degree(), order(cu.base()))
        };
        if dim1 < 1
            || dim2 < 1
            || self.dependence_table()[dim1 as usize - 1][dim2 as usize - 1] == 0
        {
            return; // `Standard_DomainError` (`cxx:524-527`).
        }

        let mut aux_h = MathMatrix::with_init(
            0,
            h.row_count() - 1,
            0,
            h.col_count() - 1,
            0.0,
        );

        let (u_first, u_last) = {
            let curve = self.curve.as_ref().unwrap().clone();
            let cu = curve.borrow();
            let knots = cu.knots();
            (knots[element], knots[element + 1])
        };

        // Quality criterion part of Hessian (`cxx:546-556`).
        h.init(0.0);
        for icrit in 0..=2 {
            self.criteria[icrit].set_interval(u_first, u_last);
            self.criteria[icrit].hessian(dim1, dim2, &mut aux_h);
            let w = self.quality_weight * self.percent[icrit] / self.estimation[icrit];
            // `H += w * AuxH` (`cxx:553`).
            let (lr, lc) = (h.lower_row(), h.lower_col());
            for r in 0..h.row_count() {
                for c in 0..h.col_count() {
                    let v = h.value(lr + r, lc + c) + w * aux_h.value(r, c);
                    h.set_value(lr + r, lc + c, v);
                }
            }
        }

        // Least square part of Hessian (`cxx:558-607`).
        aux_h.init(0.0);

        let coeff = (u_last - u_first) / 2.0;
        let deg_h = 2 * ord + 1;
        let params = self.parameters.clone().expect("hessian parameters");
        let di = self.pnt_weight.lower() as i32 - params.lower() as i32;

        if self.e != element as i32 {
            self.build_cache(element);
        }

        let mx_deg_u = mx_deg as usize;
        let deg_h_u = deg_h as usize;
        let mut ii = 1usize; // 1-based index into myCache
        let mut ipnt = self.if_;
        while ipnt <= self.il {
            let poid = self.pnt_weight.value((di + ipnt) as usize) * 2.0;
            let bv: Vec<f64> = (0..=mx_deg_u)
                .map(|j| self.cache[ii - 1 + j])
                .collect();

            for i in 0..=deg_h_u {
                let k1 = if i as i32 <= ord { i } else { i - ord as usize - 1 };
                let curcoeff = coeff.powi(k1 as i32) * poid;
                for j in i..=deg_h_u {
                    let k2 = if j as i32 <= ord { j } else { j - ord as usize - 1 };
                    let v = aux_h.value(i, j) + curcoeff * coeff.powi(k2 as i32) * bv[j];
                    aux_h.set_value(i, j, v);
                }
                for j in (deg_h_u + 1)..=mx_deg_u {
                    let v = aux_h.value(i, j) + curcoeff * bv[j];
                    aux_h.set_value(i, j, v);
                }
            }

            for i in (deg_h_u + 1)..=mx_deg_u {
                let curcoeff = bv[i] * poid;
                for j in i..=mx_deg_u {
                    let v = aux_h.value(i, j) + curcoeff * bv[j];
                    aux_h.set_value(i, j, v);
                }
            }

            ipnt += 1;
            ii += mx_deg_u + 1;
        }

        let i0 = h.lower_row();
        let j0 = h.lower_col();
        let mut i1 = i0;
        for i in 0..=mx_deg_u {
            let mut j1 = j0 + i;
            for j in i..=mx_deg_u {
                let v = h.value(i1, j1) + self.quadratic_weight * aux_h.value(i, j);
                h.set_value(i1, j1, v);
                h.set_value(j1, i1, v);
                j1 += 1;
            }
            i1 += 1;
        }
    }

    /// `Gradient(Element, Dimension, G)` (`cxx:623-736`).
    pub fn gradient(&mut self, element: usize, dim: i32, g: &mut MathVector) {
        let nb_p2d = MyLineTool::nb_p2d(&self.ssp);
        let nb_p3d = MyLineTool::nb_p3d(&self.ssp);

        if dim > 2 * nb_p2d + 3 * nb_p3d {
            return; // `Standard_DomainError` (`cxx:634-637`).
        }

        let (in_3d, ind_pnt, ind_crd) = if dim <= 3 * nb_p3d {
            let mut crd = dim % 3;
            let mut pnt = dim / 3;
            if crd == 0 {
                crd = 3;
            } else {
                pnt += 1;
            }
            (true, pnt, crd)
        } else {
            let mut crd = (dim - 3 * nb_p3d) % 2;
            let mut pnt = (dim - 3 * nb_p3d) / 2;
            if crd == 0 {
                crd = 2;
            } else {
                pnt += 1;
            }
            (false, pnt, crd)
        };

        let (u_first, u_last) = {
            let curve = self.curve.as_ref().unwrap().clone();
            let cu = curve.borrow();
            let knots = cu.knots();
            (knots[element], knots[element + 1])
        };
        let coeff = (u_last - u_first) / 2.0;

        let (ord, mx_deg) = {
            let curve = self.curve.as_ref().expect("gradient").clone();
            let cu = curve.borrow();
            (order(cu.base()), cu.base().work_degree())
        };
        let deg_h = 2 * ord + 1;
        let i0 = g.lower();
        let params = self.parameters.clone().expect("gradient parameters");
        let di = self.pnt_weight.lower() as i32 - params.lower() as i32;

        if self.e != element as i32 {
            self.build_cache(element);
        }

        g.init(0.0);

        let mx_deg_u = mx_deg as usize;
        let mut ii = 0usize; // 0-based index into the cache, OCCT `ii` starts at 1
        let mut ipnt = self.if_;
        while ipnt <= self.il {
            let (tab_p3d, tab_p2d) = MyLineTool::value(&self.ssp, ipnt);
            let pnt = if in_3d {
                let p = &tab_p3d[(ind_pnt - 1) as usize];
                match ind_crd {
                    1 => p.x(),
                    2 => p.y(),
                    _ => p.z(),
                }
            } else {
                let p = &tab_p2d[(ind_pnt - 1) as usize];
                if ind_crd == 1 {
                    p.x()
                } else {
                    p.y()
                }
            };

            let curcoeff = pnt * self.pnt_weight.value((di + ipnt) as usize);
            for i in 0..=mx_deg_u {
                let idx = i0 + i;
                let v = g.value(idx) + self.cache[ii + i] * curcoeff;
                g.set_value(idx, v);
            }
            ii += mx_deg_u + 1;
            ipnt += 1;
        }

        g.multiply_scalar(2.0 * self.quadratic_weight);

        for i in 0..=deg_h as usize {
            let k = if (i as i32) <= ord { i } else { i - ord as usize - 1 };
            let curcoeff = coeff.powi(k as i32);
            let v = g.value(i0 + i) * curcoeff;
            g.set_value(i0 + i, v);
        }
    }

    /// `InputVector(X, AssTable)` (`cxx:738-763`).
    pub fn input_vector(&mut self, x: &MathVector, ass_table: &[Vec<Vec<i32>>]) {
        let (nb_dim, nb_elm, mx_deg) = {
            let curve = self.curve.as_ref().expect("input_vector").clone();
            let cu = curve.borrow();
            (
                cu.dimension() as i32,
                cu.nb_elements(),
                cu.base().work_degree(),
            )
        };
        let mut coeff_el = Array2::new(0, mx_deg, 1, nb_dim, 0.0);
        let i0 = x.lower() as i32 - 1;

        for el in 1..=nb_elm {
            for dim in 1..=nb_dim {
                let glob_index = &ass_table[dim as usize - 1][el - 1];
                for i in 0..=mx_deg {
                    let v = x.value((i0 + glob_index[i as usize]) as usize);
                    coeff_el.set(i as usize, dim as usize, v);
                }
            }
            let curve = self.curve.as_ref().unwrap().clone();
            let mut cu = curve.borrow_mut();
            cu.set_degree(el, mx_deg);
            cu.set_element(el, &coeff_el);
        }
    }

    /// `SetWeight(QuadraticWeight, QualityWeight, percentJ1..J3)`
    /// (`cxx:766-789`).
    pub fn set_weight(&mut self, quadratic_weight: f64, quality_weight: f64, p: [f64; 3]) {
        if quadratic_weight < 0.0 || quality_weight < 0.0 || p.iter().any(|v| *v < 0.0) {
            return; // `Standard_DomainError` (`cxx:772-781`).
        }

        self.quadratic_weight = quadratic_weight;
        self.quality_weight = quality_weight;

        let total = p[0] + p[1] + p[2];
        self.percent[0] = p[0] / total;
        self.percent[1] = p[1] / total;
        self.percent[2] = p[2] / total;
    }

    /// `GetWeight(QuadraticWeight, QualityWeight)` (`cxx:792-797`).
    pub fn weight(&self) -> (f64, f64) {
        (self.quadratic_weight, self.quality_weight)
    }

    /// `SetWeight(Weight)` (`cxx:800-803`).
    pub fn set_pnt_weight(&mut self, weight: &MathVector) {
        self.pnt_weight = weight.clone();
    }

    /// `BuildCache(Element)` (`cxx:806-852`).
    fn build_cache(&mut self, element: usize) {
        let curve = self.curve.as_ref().unwrap().clone();
        let (u_first, u_last, order_, base) = {
            let cu = curve.borrow();
            let knots = cu.knots();
            (
                knots[element],
                knots[element + 1],
                cu.base().work_degree() + 1,
                cu.base().clone(),
            )
        };

        let params = self.parameters.clone().expect("build_cache parameters");
        let p_lower = params.lower() as i32;
        let p_upper = params.upper() as i32;

        self.if_ = 0;
        let mut ipnt = p_lower;
        while ipnt <= p_upper {
            let t = params.value(ipnt as usize);
            if (t > u_first && t <= u_last) || (element == 1 && t == u_first) {
                if self.if_ == 0 {
                    self.if_ = ipnt;
                }
                self.il = ipnt;
            } else if t > u_last {
                break;
            }
            ipnt += 1;
        }

        if self.if_ != 0 {
            let npnt = (self.il - self.if_ + 1) as usize;
            self.cache = vec![0.0f64; npnt * order_ as usize];

            let mut ipnt = self.if_;
            let mut ii = 0usize;
            while ipnt <= self.il {
                let t = params.value(ipnt as usize);
                let coeff = 2.0 / (u_last - u_first);
                let c0 = -(u_last + u_first) / 2.0;
                let s = (t + c0) * coeff;
                let mut basic_value = vec![0.0f64; order_ as usize];
                base.d0(s, &mut basic_value);
                self.cache[ii..ii + order_ as usize].copy_from_slice(&basic_value);
                ipnt += 1;
                ii += order_ as usize;
            }
        } else {
            // No points in the interval (`cxx:851-854`).
            self.if_ = self.il;
            self.il -= 1;
        }
        self.e = element as i32;
    }
}
