//! `FEmTool_ElementaryCriterion` and its three linear criteria.
//!
//! Source: `ModelingData/TKGeomBase/FEmTool/FEmTool_ElementaryCriterion.{hxx,cxx}`,
//! `FEmTool_LinearTension.cxx:27-241`, `FEmTool_LinearFlexion.cxx:37-233`,
//! `FEmTool_LinearJerk.cxx:35-244` (OCCT 8.0.0).
//!
//! Translation notes:
//! - The reference matrix is built from the degree-14 Hermite/Jacobi base and
//!   cached per (order, derivative order), mirroring the OCCT `static int Order`
//!   / `static math_Vector MatrixElemts` pair (`Tension.cxx:43-66`,
//!   `Flexion.cxx:43-64`, `Jerk.cxx:39-66`). The `WorkDegree > 14`
//!   `Standard_ConstructionError` is checked only when the cache has to be
//!   rebuilt, exactly like OCCT.
//! - `myCoeff` is the `(coefficient, dimension)` array handed over by
//!   `AppDef_LinearCriteria::QualityValues` (`LinearCriteria.cxx:327-330`), so
//!   `RowLength()` is the dimension count and `ColLength()` the coefficient
//!   count — see [`crate::fem_tool::Array2::row_length`].
//! - `math_Vector::Multiply(H, X)` (`math_VectorBase.lxx:436-459`) is
//!   lower-bound relative; the port reads the matrix through
//!   `lower_row()/lower_col()` offsets.

use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

use occt_core::kernel::geomabs::Shape;

use crate::fem_tool::array2::Array2;
use crate::fem_tool::gauss_set::gauss_set_integration;
use crate::fem_tool::ref_matrix::ElementsOfRefMatrix;
use crate::matrix::MathMatrix;
use crate::plib_hermit_jacobi::HermitJacobi;
use crate::vector::MathVector;

/// `math::GaussPointsMax()` (`math.cxx:24-27`).
const GAUSS_POINTS_MAX: i32 = 61;

/// `static int WDeg = 14` (`LinearTension.cxx:43`, `LinearFlexion.cxx:43`,
/// `LinearJerk.cxx:39`).
const W_DEG: i32 = 14;

/// `FEmTool_ElementaryCriterion` (`hxx:31-65`): the interface the assembly
/// drives.
pub trait Criterion {
    /// `Set(Coeff)` (`hxx:39`).
    fn set_coeff(&mut self, coeff: Array2<f64>);

    /// `Set(FirstKnot, LastKnot)` (`cxx:24-28`).
    fn set_interval(&mut self, first_knot: f64, last_knot: f64);

    /// `DependenceTable()` (`hxx:45`).
    fn dependence_table(&self) -> Option<Array2<i32>>;

    /// `Value()` (`hxx:48`).
    fn value(&mut self) -> f64;

    /// `Hessian(Dim1, Dim2, H)` (`hxx:54`).
    fn hessian(&mut self, dim1: i32, dim2: i32, h: &mut MathMatrix);

    /// `Gradient(Dim, G)` (`hxx:58`).
    fn gradient(&mut self, dim: i32, g: &mut MathVector);
}

/// The `myCoeff` / `myFirst` / `myLast` members of
/// `FEmTool_ElementaryCriterion` (`hxx:63-65`).
struct ElementaryCriteria;

impl ElementaryCriteria {
    /// `Set(FirstKnot, LastKnot)` (`FEmTool_ElementaryCriterion.cxx:24-28`).
    fn set_interval(first: &mut f64, last: &mut f64, first_knot: f64, last_knot: f64) {
        *first = first_knot;
        *last = last_knot;
    }
}

/// `PLib::NivConstr(ConstraintOrder)` (`FEmTool_LinearTension.cxx:45`).
fn niv_constr(order: Shape) -> i32 {
    occt_core::bspl::plib_jacobi::niv_constr(order)
}

/// The shared `MatrixElemts` cache: `static int Order` / `static math_Vector
/// MatrixElemts` per criterion class, keyed here by `(der_order, order)`.
fn ref_matrix_elements(
    work_degree: i32,
    constraint_order: Shape,
    der_order: i32,
    order: i32,
) -> Result<Vec<f64>, ()> {
    static CACHE: OnceLock<Mutex<HashMap<(i32, i32), Vec<f64>>>> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    let key = (der_order, order);

    let mut guard = cache.lock().unwrap();
    if let Some(elements) = guard.get(&key) {
        return Ok(elements.clone());
    }

    // `if (WorkDegree > WDeg) throw Standard_ConstructionError("Degree too
    // high")` (`Tension.cxx:47-51`), reached only on a cache miss.
    if work_degree > W_DEG {
        return Err(());
    }

    let the_base = HermitJacobi::new(W_DEG, constraint_order).map_err(|_| ())?;
    let elem = ElementsOfRefMatrix::new(&the_base, der_order).map_err(|_| ())?;

    let max_degree = W_DEG + 1;
    let an_order = std::cmp::min(4 * (max_degree / 2 + 1), GAUSS_POINTS_MAX) as usize;
    let elements = gauss_set_integration(elem.nb_equations(), an_order, -1.0, 1.0, |u, f| {
        elem.value(u, f)
    })
    .ok_or(())?;

    guard.insert(key, elements.clone());
    Ok(elements)
}

/// Fills `RefMatrix(0..WorkDegree, 0..WorkDegree)` from the cached elements
/// (`Tension.cxx:67-76`, `Flexion.cxx:65-74`, `Jerk.cxx:67-77`).
fn build_ref_matrix(work_degree: i32, elements: &[f64]) -> MathMatrix {
    let n = work_degree.max(0) as usize;
    let mut ref_matrix = MathMatrix::new(0, n, 0, n);
    let mut ii = 0usize;
    for i in 0..=n {
        ref_matrix.set_value(i, i, elements[ii]);
        let mut jj = ii + 1;
        for j in (i + 1)..=n {
            let v = elements[jj];
            ref_matrix.set_value(j, i, v);
            ref_matrix.set_value(i, j, v);
            jj += 1;
        }
        ii += (W_DEG + 1 - i as i32).max(0) as usize;
    }
    ref_matrix
}

/// `G.Multiply(H, X)` (`math_VectorBase.lxx:436-459`).
fn multiply_mat_vec(g: &mut MathVector, h: &MathMatrix, x: &MathVector) {
    let n_rows = h.row_count();
    let n_cols = h.col_count();
    let lr = h.lower_row();
    let lc = h.lower_col();
    for r in 0..n_rows {
        let mut sum = 0.0;
        for c in 0..n_cols {
            sum += h.value(lr + r, lc + c) * x.value(x.lower() + c);
        }
        g.set_value(g.lower() + r, sum);
    }
}

/// `FEmTool_LinearTension` (`hxx:31-51`).
pub struct LinearTension {
    coeff: Option<Array2<f64>>,
    first: f64,
    last: f64,
    order: i32,
    ref_matrix: MathMatrix,
}

impl LinearTension {
    /// `FEmTool_LinearTension(WorkDegree, ConstraintOrder)`
    /// (`Tension.cxx:34-77`).
    pub fn new(work_degree: i32, constraint_order: Shape) -> Result<Self, ()> {
        let order = niv_constr(constraint_order);
        let elements = ref_matrix_elements(work_degree, constraint_order, 1, order)?;
        Ok(Self {
            coeff: None,
            first: 0.0,
            last: 0.0,
            order,
            ref_matrix: build_ref_matrix(work_degree, &elements),
        })
    }
}

impl Criterion for LinearTension {
    fn set_coeff(&mut self, coeff: Array2<f64>) {
        self.coeff = Some(coeff);
    }

    fn set_interval(&mut self, first_knot: f64, last_knot: f64) {
        ElementaryCriteria::set_interval(&mut self.first, &mut self.last, first_knot, last_knot);
    }

    /// `DependenceTable()` (`Tension.cxx:79-97`).
    fn dependence_table(&self) -> Option<Array2<i32>> {
        let coeff = self.coeff.as_ref()?;
        let lc = coeff.lower_col();
        let uc = coeff.upper_col();
        let mut dep_tab = Array2::new(lc as i32, uc as i32, lc as i32, uc as i32, 0);
        // `for (i = 1; i <= myCoeff->RowLength(); i++) DepTab->SetValue(i, i, 1);`
        for i in 1..=coeff.row_length() {
            dep_tab.set_value(i, i, 1);
        }
        Some(dep_tab)
    }

    /// `Value()` (`Tension.cxx:99-148`).
    fn value(&mut self) -> f64 {
        let coeff = match self.coeff.as_ref() {
            Some(c) => c,
            None => return 0.0,
        };
        let deg = std::cmp::min(coeff.col_length() as i32 - 1, self.ref_matrix.upper_row() as i32);
        let j0 = coeff.lower_row();
        let deg_h = std::cmp::min(2 * self.order + 1, deg);
        let nb_dim = coeff.row_length();

        let mut new_coeff = Array2::new(1, nb_dim as i32, 0, deg, 0.0);

        let c = (self.last - self.first) / 2.0;
        let cteh3 = 2.0 / c;

        for i in 0..=deg {
            let k1 = if i <= self.order { i } else { i - self.order - 1 };
            let mfact = c.powf(k1 as f64);
            for dim in 1..=nb_dim {
                let v = coeff.at(j0 + i as usize, dim) * mfact;
                new_coeff.set(dim, i as usize, v);
            }
        }

        for i in (deg_h + 1)..=deg {
            for dim in 1..=nb_dim {
                let v = coeff.at(j0 + i as usize, dim);
                new_coeff.set(dim, i as usize, v);
            }
        }

        let mut j = 0.0;
        for dim in 1..=nb_dim {
            for i in 0..=deg {
                let mut jline = 0.5 * self.ref_matrix.value(i as usize, i as usize)
                    * new_coeff.at(dim, i as usize);
                for jj in 0..i {
                    jline += self.ref_matrix.value(i as usize, jj as usize)
                        * new_coeff.at(dim, jj as usize);
                }
                j += jline * new_coeff.at(dim, i as usize);
            }
        }

        cteh3 * j
    }

    /// `Hessian(Dimension1, Dimension2, H)` (`Tension.cxx:151-219`).
    fn hessian(&mut self, dim1: i32, dim2: i32, h: &mut MathMatrix) {
        let dep_tab = match self.dependence_table() {
            Some(t) => t,
            None => return,
        };

        if (dim1 as usize) < dep_tab.lower_row()
            || (dim1 as usize) > dep_tab.upper_row()
            || (dim2 as usize) < dep_tab.lower_col()
            || (dim2 as usize) > dep_tab.upper_col()
        {
            return; // `Standard_OutOfRange` (`Tension.cxx:157-160`).
        }
        if *dep_tab.value(dim1 as usize, dim2 as usize) == 0 {
            return; // `Standard_DomainError` (`Tension.cxx:162-165`).
        }

        let deg = std::cmp::min(
            self.ref_matrix.upper_row() as i32,
            h.row_count() as i32 - 1,
        );
        let deg_h = std::cmp::min(2 * self.order + 1, deg);

        let c = (self.last - self.first) / 2.0;
        let cteh3 = 2.0 / c;
        let i0 = h.lower_row();
        let j0 = h.lower_col();

        h.init(0.0);

        let mut i1 = i0;
        for i in 0..=deg {
            let k1 = if i <= self.order { i } else { i - self.order - 1 };
            let mfact = c.powf(k1 as f64) * cteh3;
            // Hermite*Hermite part of matrix.
            let mut j1 = j0 + i as usize;
            for j in i..=deg_h {
                let k2 = if j <= self.order { j } else { j - self.order - 1 };
                let v = mfact
                    * c.powf(k2 as f64)
                    * self.ref_matrix.value(i as usize, j as usize);
                h.set_value(i1, j1, v);
                if i != j {
                    h.set_value(j1, i1, v);
                }
                j1 += 1;
            }
            // Hermite*Jacobi part of matrix.
            let mut j1 = j0 + deg_h as usize + 1;
            for j in (deg_h + 1)..=deg {
                let v = mfact * self.ref_matrix.value(i as usize, j as usize);
                h.set_value(i1, j1, v);
                h.set_value(j1, i1, v);
                j1 += 1;
            }
            i1 += 1;
        }

        // Jacoby*Jacobi part of matrix.
        let mut i1 = i0 + deg_h as usize + 1;
        for i in (deg_h + 1)..=deg {
            let mut j1 = j0 + i as usize;
            for j in i..=deg {
                let v = cteh3 * self.ref_matrix.value(i as usize, j as usize);
                h.set_value(i1, j1, v);
                if i != j {
                    h.set_value(j1, i1, v);
                }
                j1 += 1;
            }
            i1 += 1;
        }
    }

    /// `Gradient(Dimension, G)` (`Tension.cxx:222-241`).
    fn gradient(&mut self, dim: i32, g: &mut MathVector) {
        let coeff = match self.coeff.as_ref() {
            Some(c) => c,
            None => return,
        };
        if (dim as usize) < coeff.lower_col() || (dim as usize) > coeff.upper_col() {
            return; // `Standard_OutOfRange` (`Tension.cxx:224-227`).
        }

        let deg = std::cmp::min(g.len() as i32 - 1, coeff.col_length() as i32 - 1);

        let mut x = MathVector::new(0, deg as usize);
        let i1 = coeff.lower_row();
        for i in 0..=deg {
            x.set_value(i as usize, coeff.at(i1 + i as usize, dim as usize));
        }

        let mut h = MathMatrix::new(0, deg as usize, 0, deg as usize);
        self.hessian(dim, dim, &mut h);

        multiply_mat_vec(g, &h, &x);
    }
}

/// `FEmTool_LinearFlexion` (`hxx:31-51`).
pub struct LinearFlexion {
    coeff: Option<Array2<f64>>,
    first: f64,
    last: f64,
    order: i32,
    ref_matrix: MathMatrix,
}

impl LinearFlexion {
    /// `FEmTool_LinearFlexion(WorkDegree, ConstraintOrder)`
    /// (`Flexion.cxx:37-75`).
    pub fn new(work_degree: i32, constraint_order: Shape) -> Result<Self, ()> {
        let order = niv_constr(constraint_order);
        let elements = ref_matrix_elements(work_degree, constraint_order, 2, order)?;
        Ok(Self {
            coeff: None,
            first: 0.0,
            last: 0.0,
            order,
            ref_matrix: build_ref_matrix(work_degree, &elements),
        })
    }
}

impl Criterion for LinearFlexion {
    fn set_coeff(&mut self, coeff: Array2<f64>) {
        self.coeff = Some(coeff);
    }

    fn set_interval(&mut self, first_knot: f64, last_knot: f64) {
        ElementaryCriteria::set_interval(&mut self.first, &mut self.last, first_knot, last_knot);
    }

    /// `DependenceTable()` (`Flexion.cxx:80-99`).
    fn dependence_table(&self) -> Option<Array2<i32>> {
        let coeff = self.coeff.as_ref()?;
        let lc = coeff.lower_col();
        let uc = coeff.upper_col();
        let mut dep_tab = Array2::new(lc as i32, uc as i32, lc as i32, uc as i32, 0);
        for i in lc..=uc {
            dep_tab.set_value(i, i, 1);
        }
        Some(dep_tab)
    }

    /// `Value()` (`Flexion.cxx:104-155`).
    fn value(&mut self) -> f64 {
        let coeff = match self.coeff.as_ref() {
            Some(c) => c,
            None => return 0.0,
        };
        let deg = std::cmp::min(coeff.col_length() as i32 - 1, self.ref_matrix.upper_row() as i32);
        let j0 = coeff.lower_row();
        let deg_h = std::cmp::min(2 * self.order + 1, deg);
        let nb_dim = coeff.row_length();

        let mut new_coeff = Array2::new(1, nb_dim as i32, 0, deg, 0.0);

        let c = (self.last - self.first) / 2.0;
        let cteh3 = 2.0 / c.powf(3.0);

        for i in 0..=deg {
            let k1 = if i <= self.order { i } else { i - self.order - 1 };
            let mfact = c.powf(k1 as f64);
            for dim in 1..=nb_dim {
                let v = coeff.at(j0 + i as usize, dim) * mfact;
                new_coeff.set(dim, i as usize, v);
            }
        }

        for i in (deg_h + 1)..=deg {
            for dim in 1..=nb_dim {
                let v = coeff.at(j0 + i as usize, dim);
                new_coeff.set(dim, i as usize, v);
            }
        }

        let mut j = 0.0;
        for dim in 1..=nb_dim {
            for i in 0..=deg {
                let mut jline = 0.5 * self.ref_matrix.value(i as usize, i as usize)
                    * new_coeff.at(dim, i as usize);
                for jj in 0..i {
                    jline += self.ref_matrix.value(i as usize, jj as usize)
                        * new_coeff.at(dim, jj as usize);
                }
                j += jline * new_coeff.at(dim, i as usize);
            }
        }

        if j < 0.0 {
            j = 0.0;
        }
        cteh3 * j
    }

    /// `Hessian(Dimension1, Dimension2, H)` (`Flexion.cxx:158-213`): unlike
    /// tension / jerk this body addresses `H(i, j)` directly, so `H` must have
    /// lower bounds `0`.
    fn hessian(&mut self, dim1: i32, dim2: i32, h: &mut MathMatrix) {
        let dep_tab = match self.dependence_table() {
            Some(t) => t,
            None => return,
        };

        if (dim1 as usize) < dep_tab.lower_row()
            || (dim1 as usize) > dep_tab.upper_row()
            || (dim2 as usize) < dep_tab.lower_col()
            || (dim2 as usize) > dep_tab.upper_col()
        {
            return; // `Standard_OutOfRange` (`Flexion.cxx:164-167`).
        }
        if *dep_tab.value(dim1 as usize, dim2 as usize) == 0 {
            return; // `Standard_DomainError` (`Flexion.cxx:169-172`).
        }

        let deg = std::cmp::min(
            self.ref_matrix.upper_row() as i32,
            h.row_count() as i32 - 1,
        );
        let deg_h = std::cmp::min(2 * self.order + 1, deg);

        let c = (self.last - self.first) / 2.0;
        let cteh3 = 2.0 / c.powf(3.0);

        h.init(0.0);

        for i in 0..=deg_h {
            let k1 = if i <= self.order { i } else { i - self.order - 1 };
            let mfact = c.powf(k1 as f64) * cteh3;
            // Hermite*Hermite part of matrix.
            for j in i..=deg_h {
                let k2 = if j <= self.order { j } else { j - self.order - 1 };
                let v = mfact
                    * c.powf(k2 as f64)
                    * self.ref_matrix.value(i as usize, j as usize);
                h.set_value(i as usize, j as usize, v);
                if i != j {
                    h.set_value(j as usize, i as usize, v);
                }
            }
            // Hermite*Jacobi part of matrix.
            for j in (deg_h + 1)..=deg {
                let v = mfact * self.ref_matrix.value(i as usize, j as usize);
                h.set_value(i as usize, j as usize, v);
                h.set_value(j as usize, i as usize, v);
            }
        }

        // Jacoby*Jacobi part of matrix.
        for i in (deg_h + 1)..=deg {
            for j in i..=deg {
                let v = cteh3 * self.ref_matrix.value(i as usize, j as usize);
                h.set_value(i as usize, j as usize, v);
                if i != j {
                    h.set_value(j as usize, i as usize, v);
                }
            }
        }
    }

    /// `Gradient(Dimension, G)` (`Flexion.cxx:218-236`).
    fn gradient(&mut self, dim: i32, g: &mut MathVector) {
        let coeff = match self.coeff.as_ref() {
            Some(c) => c,
            None => return,
        };
        if (dim as usize) < coeff.lower_col() || (dim as usize) > coeff.upper_col() {
            return; // `Standard_OutOfRange` (`Flexion.cxx:220-223`).
        }

        let deg = std::cmp::min(g.len() as i32 - 1, coeff.col_length() as i32 - 1);

        let mut x = MathVector::new(0, deg as usize);
        let mut h = MathMatrix::new(0, deg as usize, 0, deg as usize);
        let i1 = coeff.lower_row();
        for i in 0..=deg {
            x.set_value(i as usize, coeff.at(i1 + i as usize, dim as usize));
        }

        self.hessian(dim, dim, &mut h);

        multiply_mat_vec(g, &h, &x);
    }
}

/// `FEmTool_LinearJerk` (`hxx:31-51`).
pub struct LinearJerk {
    coeff: Option<Array2<f64>>,
    first: f64,
    last: f64,
    order: i32,
    ref_matrix: MathMatrix,
}

impl LinearJerk {
    /// `FEmTool_LinearJerk(WorkDegree, ConstraintOrder)` (`Jerk.cxx:35-78`).
    pub fn new(work_degree: i32, constraint_order: Shape) -> Result<Self, ()> {
        let order = niv_constr(constraint_order);
        let elements = ref_matrix_elements(work_degree, constraint_order, 3, order)?;
        Ok(Self {
            coeff: None,
            first: 0.0,
            last: 0.0,
            order,
            ref_matrix: build_ref_matrix(work_degree, &elements),
        })
    }
}

impl Criterion for LinearJerk {
    fn set_coeff(&mut self, coeff: Array2<f64>) {
        self.coeff = Some(coeff);
    }

    fn set_interval(&mut self, first_knot: f64, last_knot: f64) {
        ElementaryCriteria::set_interval(&mut self.first, &mut self.last, first_knot, last_knot);
    }

    /// `DependenceTable()` (`Jerk.cxx:80-99`).
    fn dependence_table(&self) -> Option<Array2<i32>> {
        let coeff = self.coeff.as_ref()?;
        let lc = coeff.lower_col();
        let uc = coeff.upper_col();
        let mut dep_tab = Array2::new(lc as i32, uc as i32, lc as i32, uc as i32, 0);
        for i in lc..=uc {
            dep_tab.set_value(i, i, 1);
        }
        Some(dep_tab)
    }

    /// `Value()` (`Jerk.cxx:101-156`).
    fn value(&mut self) -> f64 {
        let coeff = match self.coeff.as_ref() {
            Some(c) => c,
            None => return 0.0,
        };
        let deg = std::cmp::min(coeff.col_length() as i32 - 1, self.ref_matrix.upper_row() as i32);
        let j0 = coeff.lower_row();
        let deg_h = std::cmp::min(2 * self.order + 1, deg);
        let nb_dim = coeff.row_length();

        let mut new_coeff = Array2::new(1, nb_dim as i32, 0, deg, 0.0);

        let c = (self.last - self.first) / 2.0;
        let cteh3 = 2.0 / c.powf(5.0);

        for i in 0..=deg_h {
            let k1 = if i <= self.order { i } else { i - self.order - 1 };
            let mfact = c.powf(k1 as f64);
            for dim in 1..=nb_dim {
                let v = coeff.at(j0 + i as usize, dim) * mfact;
                new_coeff.set(dim, i as usize, v);
            }
        }

        for i in (deg_h + 1)..=deg {
            for dim in 1..=nb_dim {
                let v = coeff.at(j0 + i as usize, dim);
                new_coeff.set(dim, i as usize, v);
            }
        }

        let mut j = 0.0;
        for dim in 1..=nb_dim {
            for i in 0..=deg {
                let mut jline = 0.5 * self.ref_matrix.value(i as usize, i as usize)
                    * new_coeff.at(dim, i as usize);
                for jj in 0..i {
                    jline += self.ref_matrix.value(i as usize, jj as usize)
                        * new_coeff.at(dim, jj as usize);
                }
                j += jline * new_coeff.at(dim, i as usize);
                if j < 0.0 {
                    j = 0.0;
                }
            }
        }

        cteh3 * j
    }

    /// `Hessian(Dimension1, Dimension2, H)` (`Jerk.cxx:158-223`).
    fn hessian(&mut self, dim1: i32, dim2: i32, h: &mut MathMatrix) {
        let dep_tab = match self.dependence_table() {
            Some(t) => t,
            None => return,
        };

        if (dim1 as usize) < dep_tab.lower_row()
            || (dim1 as usize) > dep_tab.upper_row()
            || (dim2 as usize) < dep_tab.lower_col()
            || (dim2 as usize) > dep_tab.upper_col()
        {
            return; // `Standard_OutOfRange` (`Jerk.cxx:164-167`).
        }
        if *dep_tab.value(dim1 as usize, dim2 as usize) == 0 {
            return; // `Standard_DomainError` (`Jerk.cxx:169-172`).
        }

        let deg = std::cmp::min(
            self.ref_matrix.upper_row() as i32,
            h.row_count() as i32 - 1,
        );
        let deg_h = std::cmp::min(2 * self.order + 1, deg);

        let c = (self.last - self.first) / 2.0;
        let cteh3 = 2.0 / c.powf(5.0);
        let i0 = h.lower_row();
        let j0 = h.lower_col();

        h.init(0.0);

        let mut i1 = i0;
        for i in 0..=deg_h {
            let k1 = if i <= self.order { i } else { i - self.order - 1 };
            let mfact = c.powf(k1 as f64) * cteh3;
            // Hermite*Hermite part of matrix.
            let mut j1 = j0 + i as usize;
            for j in i..=deg_h {
                let k2 = if j <= self.order { j } else { j - self.order - 1 };
                let v = mfact
                    * c.powf(k2 as f64)
                    * self.ref_matrix.value(i as usize, j as usize);
                h.set_value(i1, j1, v);
                if i != j {
                    h.set_value(j1, i1, v);
                }
                j1 += 1;
            }
            // Hermite*Jacobi part of matrix.
            let mut j1 = j0 + deg_h as usize + 1;
            for j in (deg_h + 1)..=deg {
                let v = mfact * self.ref_matrix.value(i as usize, j as usize);
                h.set_value(i1, j1, v);
                h.set_value(j1, i1, v);
                j1 += 1;
            }
            i1 += 1;
        }

        // Jacoby*Jacobi part of matrix.
        let mut i1 = i0 + deg_h as usize + 1;
        for i in (deg_h + 1)..=deg {
            let mut j1 = j0 + i as usize;
            for j in i..=deg {
                let v = cteh3 * self.ref_matrix.value(i as usize, j as usize);
                h.set_value(i1, j1, v);
                if i != j {
                    h.set_value(j1, i1, v);
                }
                j1 += 1;
            }
            i1 += 1;
        }
    }

    /// `Gradient(Dimension, G)` (`Jerk.cxx:225-244`).
    fn gradient(&mut self, dim: i32, g: &mut MathVector) {
        let coeff = match self.coeff.as_ref() {
            Some(c) => c,
            None => return,
        };
        if (dim as usize) < coeff.lower_col() || (dim as usize) > coeff.upper_col() {
            return; // `Standard_OutOfRange` (`Jerk.cxx:227-230`).
        }

        let deg = std::cmp::min(g.len() as i32 - 1, coeff.col_length() as i32 - 1);

        let mut x = MathVector::new(0, deg as usize);
        let i1 = coeff.lower_row();
        for i in 0..=deg {
            x.set_value(i as usize, coeff.at(i1 + i as usize, dim as usize));
        }

        let mut h = MathMatrix::new(0, deg as usize, 0, deg as usize);
        self.hessian(dim, dim, &mut h);

        multiply_mat_vec(g, &h, &x);
    }
}
