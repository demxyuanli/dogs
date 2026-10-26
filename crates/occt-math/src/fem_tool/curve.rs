//! `FEmTool_Curve` — a piecewise polynomial curve defined by Hermite/Jacobi
//! elements.
//!
//! Source: `ModelingData/TKGeomBase/FEmTool/FEmTool_Curve.{hxx,cxx:30-580}`
//! (OCCT 8.0.0).
//!
//! Translation notes:
//! - All member arrays are stored with a leading unused slot so the Rust
//!   indices are the OCCT 1-based ones (`myCoeff(i)` → `coeff[i]`). The OCCT
//!   pointer arithmetic over the element blocks (`(Index - 1) * (WorkDegree + 1)
//!   * Dimension`) is kept verbatim.
//! - `PLib::NoDerivativeEvalPolynomial(S, Degree, Dimension, Degree * Dimension,
//!   myPoly(Ptr), …)` (`cxx:207-213`) addresses the **highest** coefficient
//!   block through `PolynomialCoeff + DegreeDimension` (`PLib.cxx:1036`); the
//!   port passes the whole element block to
//!   [`occt_core::bspl::coefs::no_derivative_eval_polynomial`], whose layout is
//!   lowest coefficient first.
//! - `PLib_HermitJacobi::ToCoefficients` returns the canonical coefficients
//!   instead of filling the out-parameter (`plib_hermit_jacobi.rs:122-177`).
//! - `Standard_OutOfRange` (`cxx:67`, `:109`, `:497`, `:503`) is reported as
//!   `false` / no-op, matching the callers in `AppDef_Variational`, which never
//!   index outside the curve.

use occt_core::bspl::coefs::no_derivative_eval_polynomial;

use crate::fem_tool::array2::Array2;
use crate::fem_tool::eval_length::eval_length;
use crate::plib_hermit_jacobi::HermitJacobi;

/// `FEmTool_Curve` (`hxx:35-146`).
#[derive(Clone)]
pub struct Curve {
    nb_elements: usize,
    dimension: usize,
    base: HermitJacobi,
    /// `myKnots` (`HArray1<double>(1, NbElements + 1)`); index 0 unused.
    knots: Vec<f64>,
    /// `myDegree` (`Array1<int>(1, NbElements)`); index 0 unused.
    degree: Vec<i32>,
    /// `myCoeff` (`Array1<double>(1, Dimension * NbElements * (WorkDegree + 1))`).
    coeff: Vec<f64>,
    /// `myPoly` — the same elements in the canonical basis.
    poly: Vec<f64>,
    /// `myDeri` (`Array1<double>(1, Dimension * NbElements * WorkDegree)`).
    deri: Vec<f64>,
    /// `myDsecn` (`Array1<double>(1, Dimension * NbElements * (WorkDegree - 1))`).
    dsecn: Vec<f64>,
    has_poly: Vec<i32>,
    has_deri: Vec<i32>,
    has_secn: Vec<i32>,
    /// `myLength` (`Array1<double>(1, NbElements)`), `-1` when not computed.
    length: Vec<f64>,
    uf: f64,
    ul: f64,
    denom: f64,
    u_sum: f64,
    index_cache: usize,
    ptr: usize,
}

impl Curve {
    /// `FEmTool_Curve(Dimension, NbElements, TheBase, Tolerance)`
    /// (`cxx:30-55`).
    pub fn new(dimension: i32, nb_elements: i32, base: &HermitJacobi, tolerance: f64) -> Self {
        let _ = tolerance; // `const double` parameter, unused in the body.
        let nb_el = nb_elements.max(0) as usize;
        let dim = dimension.max(0) as usize;
        let work_degree = base.work_degree();
        let wd = work_degree.max(0) as usize;

        let mut knots = vec![0.0f64; nb_el + 2];
        knots[1] = 0.0;
        let mut degree = vec![0i32; nb_el + 1];
        for i in 1..=nb_el {
            degree[i] = work_degree;
        }
        let mut s = Self {
            nb_elements: nb_el,
            dimension: dim,
            base: base.clone(),
            knots,
            degree,
            coeff: vec![0.0; dim * nb_el * (wd + 1) + 1],
            poly: vec![0.0; dim * nb_el * (wd + 1) + 1],
            deri: vec![0.0; dim * nb_el * wd + 1],
            dsecn: vec![0.0; dim * nb_el * wd.saturating_sub(1) + 1],
            has_poly: vec![0; nb_el + 1],
            has_deri: vec![0; nb_el + 1],
            has_secn: vec![0; nb_el + 1],
            length: vec![-1.0; nb_el + 1],
            uf: 0.0,
            ul: 0.0,
            denom: 0.0,
            u_sum: 0.0,
            index_cache: 0,
            ptr: 0,
        };
        s.knots[0] = 0.0;
        s
    }

    /// `Knots()` (`cxx:57-60`): OCCT 1-based indices, index 0 is unused.
    pub fn knots(&self) -> &[f64] {
        &self.knots
    }

    /// `Knots()` as a mutable reference (`cxx:57-60`).
    pub fn knots_mut(&mut self) -> &mut [f64] {
        &mut self.knots
    }

    /// `NbElements()` (`cxx:470-473`).
    pub fn nb_elements(&self) -> usize {
        self.nb_elements
    }

    /// `Dimension()` (`cxx:475-478`).
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// `Base()` (`cxx:480-483`).
    pub fn base(&self) -> &HermitJacobi {
        &self.base
    }

    /// `Degree(IndexOfElement)` (`cxx:485-488`).
    pub fn degree(&self, index_of_element: usize) -> i32 {
        self.degree[index_of_element]
    }

    /// `SetDegree(IndexOfElement, Degree)` (`cxx:492-506`).
    pub fn set_degree(&mut self, index_of_element: usize, degree: i32) {
        if degree <= self.base.work_degree() {
            self.degree[index_of_element] = degree;
            self.has_poly[index_of_element] = 0;
            self.has_deri[index_of_element] = 0;
            self.has_secn[index_of_element] = 0;
            self.length[index_of_element] = -1.0;
        }
        // `Degree > WorkDegree` raises `Standard_OutOfRange` (`cxx:502-505`).
    }

    /// `SetElement(IndexOfElement, Coeffs)` (`cxx:64-102`).
    pub fn set_element(&mut self, index_of_element: usize, coeffs: &Array2<f64>) {
        if index_of_element > self.nb_elements || index_of_element < 1 {
            return; // `Standard_OutOfRange` (`cxx:66-69`).
        }
        let deg_base = self.base.work_degree();
        let deg = self.degree[index_of_element] as usize;
        let dim = self.dimension;

        let i_base = (index_of_element - 1) * (deg_base.max(0) as usize + 1) * dim;
        let mut i1 = i_base as isize - dim as isize;
        let mut i2 = coeffs.lower_row() - 1;
        let j1 = coeffs.lower_col() - 1;
        for _i in 1..=deg + 1 {
            i1 += dim as isize;
            i2 += 1;
            for j in 1..=dim {
                self.coeff[(i1 + j as isize) as usize] = coeffs.at(i2, j1 + j);
            }
        }

        let stenor = (self.knots[index_of_element + 1] - self.knots[index_of_element]) / 2.0;
        let niv_constr = self.base.niv_constr().max(0) as usize;
        let mut i1 = i_base;
        let mut i2 = i_base + (niv_constr + 1) * dim;
        for i in 1..=niv_constr {
            i1 += dim;
            i2 += dim;
            let mfact = stenor.powf(i as f64);
            for j in 1..=dim {
                self.coeff[i1 + j] *= mfact;
                self.coeff[i2 + j] *= mfact;
            }
        }

        self.has_poly[index_of_element] = 0;
        self.has_deri[index_of_element] = 0;
        self.has_secn[index_of_element] = 0;
        self.length[index_of_element] = -1.0;
    }

    /// `GetElement(IndexOfElement, Coeffs)` (`cxx:107-141`).
    pub fn get_element(&self, index_of_element: usize, coeffs: &mut Array2<f64>) {
        if index_of_element > self.nb_elements || index_of_element < 1 {
            return; // `Standard_OutOfRange` (`cxx:109-112`).
        }
        let deg_base = self.base.work_degree();
        let deg = self.degree[index_of_element] as usize;
        let dim = self.dimension;

        let i_base = (index_of_element - 1) * (deg_base.max(0) as usize + 1) * dim;
        let mut i1 = i_base as isize - dim as isize;
        let mut i2 = coeffs.lower_row() - 1;
        let j1 = coeffs.lower_col() - 1;
        for _i in 1..=deg + 1 {
            i1 += dim as isize;
            i2 += 1;
            for j in 1..=dim {
                coeffs.set(i2, j1 + j, self.coeff[(i1 + j as isize) as usize]);
            }
        }

        let stenor = 2.0 / (self.knots[index_of_element + 1] - self.knots[index_of_element]);
        let i2 = coeffs.lower_row();
        let i3 = i2 + self.base.niv_constr().max(0) as usize + 1;
        for i in 1..=self.base.niv_constr().max(0) as usize {
            let mfact = stenor.powf(i as f64);
            for j in (j1 + 1)..=dim {
                *coeffs.change_value(i2 + i, j) *= mfact;
                *coeffs.change_value(i3 + i, j) *= mfact;
            }
        }
    }

    /// `GetPolynom(Coeffs)` (`cxx:145-163`): writes the canonical coefficients
    /// of every element, `Coeffs(i - 1) = myPoly(i)`.
    pub fn get_polynom(&mut self, coeffs: &mut [f64]) {
        for index_of_element in 1..=self.nb_elements {
            if self.has_poly[index_of_element] == 0 {
                self.update(index_of_element, 0);
            }
        }
        let n = self.poly.len().min(coeffs.len() + 1);
        for i in 1..n {
            coeffs[i - 1] = self.poly[i];
        }
    }

    /// `D0(U, Pnt)` (`cxx:167-213`).
    pub fn d0(&mut self, u: f64, pnt: &mut [f64]) {
        self.search_span(u);
        let idx = self.index_cache;
        let deg = self.degree[idx];
        if self.has_poly[idx] == 0 {
            self.update(idx, 0);
        }

        // Parameter normalization: S [-1, 1].
        let s = (2.0 * u - self.u_sum) * self.denom;
        let r = no_derivative_eval_polynomial(s, deg, self.dimension as i32, &self.poly[self.ptr..]);
        for (k, v) in r.iter().enumerate() {
            if k < pnt.len() {
                pnt[k] = *v;
            }
        }
    }

    /// `D1(U, Vec)` (`cxx:217-273`).
    pub fn d1(&mut self, u: f64, vec: &mut [f64]) {
        self.search_span(u);
        let idx = self.index_cache;
        let deg = self.degree[idx];
        if self.has_deri[idx] == 0 {
            self.update(idx, 1);
        }

        let s = (2.0 * u - self.u_sum) * self.denom;
        let base = 1 + (idx - 1) * self.base.work_degree().max(0) as usize * self.dimension;
        let r = no_derivative_eval_polynomial(
            s,
            deg - 1,
            self.dimension as i32,
            &self.deri[base..],
        );
        let scale = 2.0 * self.denom;
        for (k, v) in r.iter().enumerate() {
            if k < vec.len() {
                vec[k] = *v * scale;
            }
        }
    }

    /// `D2(U, Vec)` (`cxx:277-333`).
    pub fn d2(&mut self, u: f64, vec: &mut [f64]) {
        self.search_span(u);
        let idx = self.index_cache;
        let deg = self.degree[idx];
        if self.has_secn[idx] == 0 {
            self.update(idx, 2);
        }

        let s = (2.0 * u - self.u_sum) * self.denom;
        let stride = self.base.work_degree().max(0) as usize - 1;
        let base = 1 + (idx - 1) * stride * self.dimension;
        let r = no_derivative_eval_polynomial(
            s,
            deg - 2,
            self.dimension as i32,
            &self.dsecn[base..],
        );
        let scale = 4.0 * self.denom * self.denom;
        for (k, v) in r.iter().enumerate() {
            if k < vec.len() {
                vec[k] = *v * scale;
            }
        }
    }

    /// `Length(FirstU, LastU, Length)` (`cxx:337-468`).
    pub fn length_of(&mut self, first_u: f64, last_u: f64) -> f64 {
        if first_u > last_u {
            return 0.0; // `Standard_OutOfRange` (`cxx:339-342`).
        }
        let n = self.nb_elements;

        let low = if self.knots[1] > first_u {
            1
        } else {
            let mut low = 1usize;
            while low <= n {
                if first_u >= self.knots[low] && first_u <= self.knots[low + 1] {
                    break;
                }
                low += 1;
            }
            low
        };
        let mut low = low;
        if low > n {
            low = n;
        }

        let high = if self.knots[1] > last_u {
            1
        } else {
            let mut high = low;
            while high <= n {
                if last_u >= self.knots[high] && last_u <= self.knots[high + 1] {
                    break;
                }
                high += 1;
            }
            high
        };
        let mut high = high;
        if self.knots[n + 1] < last_u {
            high = n;
        }

        let deg_base = self.base.work_degree().max(0) as usize;
        let mut length = 0.0;

        let first_s = (2.0 * first_u - self.knots[low] - self.knots[low + 1])
            / (self.knots[low + 1] - self.knots[low]);
        let last_s = (2.0 * last_u - self.knots[high] - self.knots[high + 1])
            / (self.knots[high + 1] - self.knots[high]);

        if low == high {
            let ptr = (low - 1) * (deg_base + 1) * self.dimension + 1;
            let deg = self.degree[low];
            if self.has_poly[low] == 0 {
                self.update(low, 0);
            }
            return eval_length(
                deg,
                self.dimension as i32,
                &self.poly[ptr..],
                first_s,
                last_s,
            );
        }

        let deg = self.degree[low];
        let ptr = (low - 1) * (deg_base + 1) * self.dimension + 1;
        if self.has_poly[low] == 0 {
            self.update(low, 0);
        }
        if first_s < -1.0 {
            let li = eval_length(deg, self.dimension as i32, &self.poly[ptr..], first_s, -1.0);
            length += li;
            if self.length[low] < 0.0 {
                let li = eval_length(deg, self.dimension as i32, &self.poly[ptr..], -1.0, 1.0);
                self.length[low] = li;
            }
            length += self.length[low];
        } else {
            let li = eval_length(deg, self.dimension as i32, &self.poly[ptr..], first_s, 1.0);
            length += li;
        }

        let deg = self.degree[high];
        let ptr = (high - 1) * (deg_base + 1) * self.dimension + 1;
        if self.has_poly[high] == 0 {
            self.update(high, 0);
        }
        if last_s > 1.0 {
            let li = eval_length(deg, self.dimension as i32, &self.poly[ptr..], 1.0, last_s);
            length += li;
            if self.length[high] < 0.0 {
                let li = eval_length(deg, self.dimension as i32, &self.poly[ptr..], -1.0, 1.0);
                self.length[high] = li;
            }
            length += self.length[high];
        } else {
            let li = eval_length(deg, self.dimension as i32, &self.poly[ptr..], -1.0, last_s);
            length += li;
        }

        for i in (low + 1)..high {
            if self.length[i] < 0.0 {
                let ptr = (i - 1) * (deg_base + 1) * self.dimension + 1;
                let deg = self.degree[i];
                if self.has_poly[i] == 0 {
                    self.update(i, 0);
                }
                let li = eval_length(deg, self.dimension as i32, &self.poly[ptr..], -1.0, 1.0);
                self.length[i] = li;
            }
            length += self.length[i];
        }

        length
    }

    /// `ReduceDegree(IndexOfElement, Tol, NewDegree, MaxError)`
    /// (`cxx:509-527`).
    pub fn reduce_degree(&mut self, index_of_element: usize, tol: f64) -> (i32, f64) {
        let deg = self.degree[index_of_element];
        let ptr = (index_of_element - 1) * (self.base.work_degree().max(0) as usize + 1)
            * self.dimension
            + 1;

        let (mut new_degree, max_error) =
            self.base
                .reduce_degree(self.dimension as i32, deg, tol, &self.coeff[ptr..]);
        new_degree = std::cmp::max(new_degree, 2 * self.base.niv_constr() + 1);

        if new_degree < deg {
            self.degree[index_of_element] = new_degree;
            self.has_poly[index_of_element] = 0;
            self.has_deri[index_of_element] = 0;
            self.has_secn[index_of_element] = 0;
            self.length[index_of_element] = -1.0;
        }
        (new_degree, max_error)
    }

    /// The span search shared by `D0` / `D1` / `D2` (`cxx:172-199`).
    fn search_span(&mut self, u: f64) {
        let needs_search = self.index_cache == 0
            || u < self.uf
            || u > self.ul
            || self.knots[self.index_cache] != self.uf
            || self.knots[self.index_cache + 1] != self.ul;

        if needs_search {
            if u <= self.knots[2] {
                self.index_cache = 1;
            } else {
                let mut index = 2usize;
                while index <= self.nb_elements {
                    if u >= self.knots[index] && u <= self.knots[index + 1] {
                        break;
                    }
                    index += 1;
                }
                if index > self.nb_elements {
                    index = self.nb_elements;
                }
                self.index_cache = index;
            }
            let idx = self.index_cache;
            self.uf = self.knots[idx];
            self.ul = self.knots[idx + 1];
            self.denom = 1.0 / (self.ul - self.uf);
            self.u_sum = self.uf + self.ul;
            self.ptr = (idx - 1) * (self.base.work_degree().max(0) as usize + 1) * self.dimension + 1;
        }
    }

    /// `Update(Index, Order)` (`cxx:531-579`).
    fn update(&mut self, index: usize, order: i32) {
        let deg_base = self.base.work_degree().max(0) as usize;
        let deg = self.degree[index].max(0) as usize;
        let dim = self.dimension;

        if self.has_poly[index] == 0 {
            let ptr = (index - 1) * (deg_base + 1) * dim + 1;
            let n = dim * (deg + 1);
            let coeffs = self.base.to_coefficients(
                dim as i32,
                deg as i32,
                &self.coeff[ptr..ptr + n],
            );
            self.poly[ptr..ptr + n].copy_from_slice(&coeffs);
            self.has_poly[index] = 1;
        }

        if order >= 1 {
            let mut i1 = (index - 1) * deg_base * dim - dim;
            let mut i2 = (index - 1) * (deg_base + 1) * dim;
            if self.has_deri[index] == 0 {
                for i in 1..=deg {
                    i1 += dim;
                    i2 += dim;
                    for j in 1..=dim {
                        self.deri[i1 + j] = i as f64 * self.poly[i2 + j];
                    }
                }
                self.has_deri[index] = 1;
            }
            if order >= 2 && self.has_secn[index] == 0 {
                let stride = deg_base.saturating_sub(1);
                let mut i1 = (index - 1) * stride * dim - dim;
                let mut i2 = (index - 1) * deg_base * dim;
                for i in 1..deg {
                    i1 += dim;
                    i2 += dim;
                    for j in 1..=dim {
                        self.dsecn[i1 + j] = i as f64 * self.deri[i2 + j];
                    }
                }
                self.has_secn[index] = 1;
            }
        }
    }
}
