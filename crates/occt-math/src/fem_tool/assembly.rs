//! `FEmTool_Assembly` — global finite-element assembly of the elementary
//! matrices / vectors, the linear constraints `G X = C`, and the condensed
//! solve `H X = B - Gt l`.
//!
//! Source: `ModelingData/TKGeomBase/FEmTool/FEmTool_Assembly.cxx:28-698` and
//! `FEmTool_Assembly.hxx:30-105` (OCCT 8.0.0).
//!
//! Translation notes:
//! - Global variable numbers are kept 1:1 with OCCT (`B` spans
//!   `MinIndex(Table)..=MaxIndex(Table)`); nothing is renumbered.
//! - `NCollection_List<HArray1<double>>` for each `G` line becomes an ordered
//!   `Vec<Array1d>`; `Array1d` reproduces the bounds-preserving `HArray1`.
//! - `OCCT_DEBUG` dumps (`cxx:213-218`) and the commented-out retry loops
//!   (`cxx:219-232`, `cxx:385-392`) are not translated.

use crate::fem_tool::array1::Array1d;
use crate::fem_tool::profile_matrix::{ProfileMatrix, SparseMatrix};
use crate::matrix::MathMatrix;
use crate::vector::MathVector;

/// One line of the constraint matrix `G` together with its `C` value.
struct Constraint {
    terms: Vec<Array1d>,
    value: f64,
}

/// `FEmTool_Assembly` (`hxx:37-105`).
pub struct Assembly {
    /// `myDepTable` — 1 when dimensions `d1`, `d2` are coupled.
    dep_table: Vec<Vec<i32>>,
    /// `myRefTable` — `[dimension][element]` global variable numbers.
    ref_table: Vec<Vec<Vec<i32>>>,
    is_solved: bool,
    h: ProfileMatrix,
    b: MathVector,
    g_h_gt: Option<ProfileMatrix>,
    g: Vec<Constraint>,
}

impl Assembly {
    /// `FEmTool_Assembly(Dependence, Table)` (`cxx:90-137`).
    ///
    /// `dependence[dim1][dim2]` and `ref_table[dim][element]` are 0-based
    /// slices of the OCCT arrays, whose own bounds are 1-based.
    pub fn new(dependence: &[Vec<i32>], ref_table: &[Vec<Vec<i32>>]) -> Self {
        let mut imin = i32::MAX;
        let mut imax = i32::MIN;
        for dim in ref_table {
            for el in dim {
                for &v in el {
                    imin = imin.min(v);
                    imax = imax.max(v);
                }
            }
        }
        if imin == i32::MAX {
            imin = 1;
            imax = 1;
        }

        let b = MathVector::new(imin as usize, imax as usize);

        // FirstIndexes(1..B.Length()) initialized to B.Length() (`cxx:101-103`).
        let n = b.len();
        let mut first_indexes = vec![n as i32; n];

        for dim in ref_table {
            for el in dim {
                let i_min = el.iter().copied().min().unwrap_or(imin);
                for &v in el {
                    let i = (v - b_lower(&b)) as usize;
                    first_indexes[i] = first_indexes[i].min(i_min);
                }
            }
        }

        let mut dep_table = dependence.to_vec();
        // `myDepTable(1, Dependence.ColLength(), 1, Dependence.RowLength())`
        // copied from `Dependence` is a transposed raw copy; `Dependence` is
        // square and symmetric in `AppDef_LinearCriteria::DependenceTable`
        // (`cxx:286-293`), so the identity layout is kept as given.
        for row in dep_table.iter_mut() {
            for v in row.iter_mut() {
                if *v != 0 {
                    *v = 1;
                }
            }
        }

        let mut asm = Self {
            dep_table,
            ref_table: ref_table.to_vec(),
            is_solved: false,
            h: ProfileMatrix::new(&first_indexes),
            b,
            g_h_gt: None,
            g: Vec::new(),
        };
        asm.nullify_matrix();
        asm.nullify_vector();
        asm
    }

    /// `NullifyMatrix()` (`cxx:139-143`).
    pub fn nullify_matrix(&mut self) {
        self.h.init(0.0);
        self.is_solved = false;
    }

    /// `AddMatrix(Element, Dimension1, Dimension2, Mat)` (`cxx:146-180`).
    pub fn add_matrix(&mut self, element: usize, dim1: usize, dim2: usize, mat: &MathMatrix) {
        if self.dep_table[dim1 - 1][dim2 - 1] == 0 {
            return; // `Standard_DomainError` (`cxx:153-156`).
        }

        let t1 = &self.ref_table[dim1 - 1][element - 1];
        let t2 = &self.ref_table[dim2 - 1][element - 1];

        let nvaru = std::cmp::min(t1.len(), mat.row_count());
        let i0 = 1 - b_lower(&self.b);
        let mat_r0 = mat.lower_row() as i32;
        let mat_c0 = mat.lower_col() as i32;

        for i in 0..nvaru {
            let ii = mat_r0 + (i as i32 + 1);
            let i_glob = t1[i] + i0;
            for j in 0..=i {
                let j_glob = t2[j] + i0;
                let col = mat_c0 + (j as i32 + 1);
                let v = mat.value(ii as usize, col as usize);
                *self.h.change_value(i_glob as usize, j_glob as usize) += v;
            }
        }

        self.is_solved = false;
    }

    /// `NullifyVector()` (`cxx:184-187`).
    pub fn nullify_vector(&mut self) {
        self.b.init(0.0);
    }

    /// `AddVector(Element, Dimension, Vec)` (`cxx:189-203`).
    pub fn add_vector(&mut self, element: usize, dim: usize, vec: &MathVector) {
        let t = &self.ref_table[dim - 1][element - 1];
        let nvaru = std::cmp::min(t.len(), vec.len());
        let i0 = vec.lower() as i32 - 1;
        let b_low = b_lower(&self.b);

        for i in 0..nvaru {
            let gi = (t[i] - b_low) as usize;
            let v = self.b.value(b_low as usize + gi) + vec.value((i0 + i as i32 + 1) as usize);
            self.b.set_value(b_low as usize + gi, v);
        }
    }

    /// Number of global variables (`cxx:466-471`).
    pub fn nb_glob_var(&self) -> usize {
        self.b.len()
    }

    /// `AssemblyTable()` (`hxx:85-88`).
    pub fn assembly_table(&self) -> &[Vec<Vec<i32>>] {
        &self.ref_table
    }

    /// `ResetConstraint()` (`cxx:478-483`).
    pub fn reset_constraint(&mut self) {
        self.g.clear();
    }

    /// `NullifyConstraint()` (`cxx:485-499`).
    pub fn nullify_constraint(&mut self) {
        for c in self.g.iter_mut() {
            c.value = 0.0;
            for t in c.terms.iter_mut() {
                t.init(0.0);
            }
        }
    }

    /// `AddConstraint(IndexofConstraint, Element, Dimension, LinearForm, Value)`
    /// (`cxx:504-697`).
    pub fn add_constraint(
        &mut self,
        index: usize,
        element: usize,
        dim: usize,
        linear_form: &MathVector,
        value: f64,
    ) {
        while self.g.len() < index {
            self.g.push(Constraint { terms: Vec::new(), value: 0.0 });
        }

        let indexes = self.ref_table[dim - 1][element - 1].clone();
        let mut i_min = self.nb_glob_var() as i32;
        let mut i_max = 0i32;
        for &v in &indexes {
            i_min = i_min.min(v);
            i_max = i_max.max(v);
        }

        let terms = &mut self.g[index - 1].terms;

        // Position `target` of the coefficient array that receives the
        // contribution (`Coeff` in `cxx:530-688`).
        let target: usize;

        if terms.is_empty() {
            let mut coeff = Array1d::zeros(i_min, i_max);
            coeff.init(0.0);
            terms.push(coeff);
            target = 0;
        } else {
            let mut s1 = 0usize;
            let mut s2 = 0usize;
            let mut aux1: Option<Array1d> = None;
            let mut aux2: Option<Array1d> = None;

            for (idx, a) in terms.iter().enumerate() {
                let i = idx + 1;
                if i_min >= a.lower() {
                    s1 = i;
                    aux1 = Some(a.clone());
                    if i_max <= a.upper() {
                        s2 = s1;
                        break;
                    }
                }
                if i_max <= a.upper() {
                    s2 = i;
                    aux2 = Some(a.clone());
                }
            }

            if s1 != s2 {
                let a1 = aux1.clone();
                let a2 = aux2.clone();
                if s1 == 0 {
                    let a2 = a2.expect("add_constraint: s1 == 0 without Aux2");
                    if i_max < a2.lower() {
                        // Inserting before the first segment.
                        let mut coeff = Array1d::zeros(i_min, i_max);
                        coeff.init(0.0);
                        terms.insert(0, coeff);
                        target = 0;
                    } else {
                        // Merge the new and the first segment.
                        let mut coeff = Array1d::zeros(i_min, a2.upper());
                        for i in i_min..a2.lower() {
                            coeff.set_value(i, 0.0);
                        }
                        for i in a2.lower()..=a2.upper() {
                            coeff.set_value(i, a2.value(i));
                        }
                        terms[0] = coeff;
                        target = 0;
                    }
                } else if s2 == 0 {
                    let a1 = a1.expect("add_constraint: s2 == 0 without Aux1");
                    if i_min > a1.upper() {
                        let mut coeff = Array1d::zeros(i_min, i_max);
                        coeff.init(0.0);
                        terms.push(coeff);
                        target = terms.len() - 1;
                    } else {
                        // Merge the new and the last segment (`L.Last() = Coeff`).
                        let mut coeff = Array1d::zeros(a1.lower(), i_max);
                        for i in a1.lower()..=a1.upper() {
                            coeff.set_value(i, a1.value(i));
                        }
                        for i in (a1.upper() + 1)..=i_max {
                            coeff.set_value(i, 0.0);
                        }
                        let last = terms.len() - 1;
                        terms[last] = coeff;
                        target = last;
                    }
                } else {
                    let a1 = a1.expect("add_constraint: s1 set without Aux1");
                    let a2 = a2.expect("add_constraint: s2 set without Aux2");
                    if i_min <= a1.upper() && i_max < a2.lower() {
                        // Merge s1 and the new.
                        let mut coeff = Array1d::zeros(a1.lower(), i_max);
                        for i in a1.lower()..=a1.upper() {
                            coeff.set_value(i, a1.value(i));
                        }
                        for i in (a1.upper() + 1)..=i_max {
                            coeff.set_value(i, 0.0);
                        }
                        terms[s1 - 1] = coeff;
                        target = s1 - 1;
                    } else if i_min > a1.upper() && i_max >= a2.lower() {
                        // Merge the new and the first segment.
                        let mut coeff = Array1d::zeros(i_min, a2.upper());
                        for i in i_min..a2.lower() {
                            coeff.set_value(i, 0.0);
                        }
                        for i in a2.lower()..=a2.upper() {
                            coeff.set_value(i, a2.value(i));
                        }
                        terms[s2 - 1] = coeff;
                        target = s2 - 1;
                    } else if i_min > a1.upper() && i_max < a2.lower() {
                        // Inserting the new between s1 and s2.
                        let mut coeff = Array1d::zeros(i_min, i_max);
                        coeff.init(0.0);
                        terms.insert(s1, coeff);
                        target = s1;
                    } else {
                        // Merge s1, the new and s2, then remove s2.
                        let mut coeff = Array1d::zeros(a1.lower(), a2.upper());
                        for i in a1.lower()..=a1.upper() {
                            coeff.set_value(i, a1.value(i));
                        }
                        for i in (a1.upper() + 1)..a2.lower() {
                            coeff.set_value(i, 0.0);
                        }
                        for i in a2.lower()..=a2.upper() {
                            coeff.set_value(i, a2.value(i));
                        }
                        terms[s1 - 1] = coeff;
                        terms.remove(s2 - 1);
                        target = s1 - 1;
                    }
                }
            } else {
                // `Coeff = Iter.Value()` — the existing segment at `s1`.
                target = s1 - 1;
            }
        }

        let linear_upper = linear_form.upper() as i32;
        let linear_lower = linear_form.lower() as i32;
        let mut j = linear_lower;
        for &gv in &indexes {
            let v = linear_form.value(j as usize);
            self.g[index - 1].terms[target].add_at(gv, v);
            j += 1;
            if j > linear_upper {
                break;
            }
        }

        self.g[index - 1].value += value;
    }
}

#[inline]
fn b_lower(b: &MathVector) -> i32 {
    b.lower() as i32
}

impl Assembly {
    /// `Solve()` (`cxx:207-395`).
    pub fn solve(&mut self) -> bool {
        self.is_solved = self.h.decompose();
        if self.is_solved && !self.g.is_empty() {
            self.solve_with_constraints();
        }
        self.is_solved
    }

    fn solve_with_constraints(&mut self) {
        let nvar = self.h.row_number();
        let n_low = b_lower(&self.b);

        if self.g_h_gt.is_none() || self.g_h_gt.as_ref().map(|g| g.row_number()) != Some(self.g.len()) {
            // H1 — the sparsity pattern of H with the not-in-profile blocks
            // zeroed (`cxx:240-275`).
            let mut h1 = vec![1u8; (nvar + 1) * (nvar + 1)];
            let at = |i: usize, j: usize| i * (nvar + 1) + j;
            let mut block_beg = 1usize;
            let mut i = 2usize;
            while i <= nvar {
                let block_end = i - 1;
                if !self.h.is_in_profile(i, block_end) {
                    let mut block = true;
                    let mut j = i + 1;
                    while j <= nvar {
                        if self.h.is_in_profile(j, block_end) {
                            block = false;
                            break;
                        }
                        j += 1;
                    }
                    if block {
                        for jj in i..=nvar {
                            for k in block_beg..=block_end {
                                h1[at(jj, k)] = 0;
                                h1[at(k, jj)] = 0;
                            }
                        }
                        block_beg = block_end + 1;
                    } else {
                        i = j;
                    }
                }
                i += 1;
            }

            let ng = self.g.len();
            let mut first_indexes = vec![1i32; ng];
            for i in 1..=ng {
                for j in (1..=i).rev() {
                    let mut zero = true;
                    'outer: for a in &self.g[i - 1].terms {
                        for k in a.lower()..=a.upper() {
                            for b in &self.g[j - 1].terms {
                                for l in b.lower()..=b.upper() {
                                    if h1[at(k as usize, l as usize)] != 0 {
                                        zero = false;
                                        break 'outer;
                                    }
                                }
                            }
                        }
                    }
                    if !zero {
                        first_indexes[i - 1] = j as i32;
                        break;
                    }
                }
            }
            self.g_h_gt = Some(ProfileMatrix::new(&first_indexes));
        }

        let g_h_gt = self.g_h_gt.as_mut().expect("GHGt assigned above");
        g_h_gt.init(0.0);

        let ng = self.g.len();
        for i in 1..=ng {
            let mut gi = MathVector::new(n_low as usize, (nvar as i32 + n_low - 1) as usize);
            for a in &self.g[i - 1].terms {
                for j in a.lower()..=a.upper() {
                    gi.set_value(j as usize, a.value(j));
                }
            }
            let mut qi = MathVector::new(n_low as usize, (nvar as i32 + n_low - 1) as usize);
            self.h.solve(&gi, &mut qi);

            for k in i..=ng {
                if !g_h_gt.is_in_profile(k, i) {
                    continue;
                }
                let mut m = 0.0;
                for a in &self.g[k - 1].terms {
                    for j in a.lower()..=a.upper() {
                        m += qi.value(j as usize) * a.value(j);
                    }
                }
                *g_h_gt.change_value(k, i) = m;
            }
        }

        self.is_solved = g_h_gt.decompose();
    }

    /// `Solution(Solution)` (`cxx:397-463`).
    pub fn solution(&self, solution: &mut MathVector) {
        let n_low = b_lower(&self.b);
        let n_high = self.b.upper() as i32;

        if self.g.is_empty() {
            self.h.solve(&self.b, solution);
            return;
        }

        let mut v1 = MathVector::new(n_low as usize, n_high as usize);
        self.h.solve(&self.b, &mut v1);

        let ng = self.g.len();
        let mut v2 = MathVector::new(1, ng);
        for i in 1..=ng {
            let mut m = 0.0;
            for a in &self.g[i - 1].terms {
                for j in a.lower()..=a.upper() {
                    m += v1.value(j as usize) * a.value(j);
                }
            }
            v2.set_value(i, m - self.g[i - 1].value);
        }

        let mut l = MathVector::new(1, ng);
        let g_h_gt = self.g_h_gt.as_ref().expect("Solution requires a decomposed GHGt");
        g_h_gt.solve(&v2, &mut l);

        v1 = self.b.clone();
        for i in 1..=ng {
            for a in &self.g[i - 1].terms {
                for j in a.lower()..=a.upper() {
                    let nv = v1.value(j as usize) - l.value(i) * a.value(j);
                    v1.set_value(j as usize, nv);
                }
            }
        }

        self.h.solve(&v1, solution);
    }
}

/// `HArray1<double>` from a plain slice, 1-based.
pub fn harray_from_slice(values: &[f64]) -> Array1d {
    Array1d::from_slice(1, values)
}
