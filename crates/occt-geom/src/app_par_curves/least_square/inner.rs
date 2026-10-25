//! The internal contact-system builders of `AppParCurves_LeastSquare`:
//! `SearchIndex` (`gxx:1814-1860`) and the three `MakeTAA` overloads
//! (`gxx:1511-1704`, `:1706-1764`, `:1766-1812`).
//!
//! C++ overloads on the second parameter's type become three names here:
//! [`LeastSquare::make_taa_vector`] (`MakeTAA(math_Vector&, math_Vector&)`),
//! [`LeastSquare::make_taa_matrix`] (`MakeTAA(math_Vector&, math_Matrix&)`) and
//! [`LeastSquare::make_taa_packed`] (`MakeTAA(math_Vector&)`).

use super::*;

impl LeastSquare {
    /// `SearchIndex(Index)` (`gxx:1814-1860`): the packed positions of the
    /// diagonal of the contact system.
    pub(crate) fn search_index(&self, index: &mut MathIntVector) {
        let nincx = self.resfin - self.resinit + 1;
        index.set_value(1, 1);

        if self.myknots.is_none() {
            if self.resinit <= self.resfin {
                let mut l = 1;
                for i in 2..=nincx {
                    l += 1;
                    let v = index.value(l - 1) + i;
                    index.set_value(l, v);
                }
            }
        } else {
            let mut iinit = self.resinit;
            let mut jinit = self.resinit;
            let mut ifin = (self.deg + 1).min(self.resfin);
            let len = self.myknots.as_ref().unwrap().len() as i32;
            let mults = self.mymults.clone().unwrap_or_default();
            let mut i2 = 1;

            for k in 2..=len {
                for i in iinit..=ifin {
                    for _j in jinit..=i {
                        if i2 != 1 {
                            let v = index.value(i2 - 1) + i - jinit + 1;
                            index.set_value(i2 as usize, v);
                        }
                    }
                    i2 += 1;
                }
                iinit = ifin + 1;
                let d = ifin + mults[(k - 1) as usize];
                ifin = d.min(self.resfin);
                jinit = (d - self.deg).max(self.resinit);
            }
        }
    }

    /// `MakeTAA(TheA, myTAB)` (`gxx:1511-1704`): the tangent-constrained
    /// contact system, linearised into the packed vectors `TheA` / `myTAB`.
    pub(crate) fn make_taa_vector(&self, the_a: &mut MathVector, my_tab: &mut MathVector) {
        let nincx = self.resfin - self.resinit + 1;
        let neq = self.last_p - self.first_p + 1;

        let ninc1 = if at_least_tangency(self.first_constraint)
            && at_least_tangency(self.last_constraint)
        {
            self.ninc - 1
        } else {
            self.ninc
        };

        let myfirst = self.a.lower_row() as i32;
        let mylast = myfirst + self.nlignes - 1;
        let mut taf1 = 0.0;
        let mut taf2 = 0.0;
        let mut taf3 = 0.0;
        let mut tab1 = 0.0;
        let mut tab2 = 0.0;
        let na1 = self.na - 1;
        let mut my_b = MathVector::new(myfirst as usize, mylast as usize);
        let mut my_v1 = MathVector::new(myfirst as usize, mylast as usize);
        let mut my_v2 = MathVector::new(myfirst as usize, mylast as usize);
        let mut the_v1 = MathVector::new(1, self.ninc as usize);
        let mut the_v2 = MathVector::new(1, self.ninc as usize);

        let first_t = at_least_tangency(self.first_constraint);
        let last_t = at_least_tangency(self.last_constraint);
        let first_p_pass = at_least_pass(self.first_constraint);
        let last_p_pass = at_least_pass(self.last_constraint);

        for i in self.first_p..=self.last_p {
            let ai2 = self.a.value(i as usize, 2);
            let aid = self.a.value(i as usize, self.nbpoles as usize);
            let mut xx = 0.0;
            let mut yy = 0.0;
            if first_p_pass {
                xx = self.a.value(i as usize, 1);
            }
            if first_t {
                xx += ai2;
            }
            if last_p_pass {
                yy = self.a.value(i as usize, self.nbpoles as usize);
            }
            if last_t {
                yy += aid;
            }
            let mut i2 = 1;
            let mut nrow = myfirst - self.first_p;
            for _ci in 1..=self.nbp {
                let i21 = i2 + 1;
                let i22 = i2 + 2;
                let ix = i + nrow;
                let iy = ix + neq;
                let iz = iy + neq;
                if first_t {
                    my_v1.set_value(ix as usize, ai2 * self.vec1t.value(i2 as usize));
                    my_v1.set_value(iy as usize, ai2 * self.vec1t.value(i21 as usize));
                    my_v1.set_value(iz as usize, ai2 * self.vec1t.value(i22 as usize));
                }
                if last_t {
                    my_v2.set_value(ix as usize, -aid * self.vec2t.value(i2 as usize));
                    my_v2.set_value(iy as usize, -aid * self.vec2t.value(i21 as usize));
                    my_v2.set_value(iz as usize, -aid * self.vec2t.value(i22 as usize));
                }
                let b_ix = self.mypoints.value(i as usize, i2 as usize)
                    - xx * self.mypoints.value(self.myfirstp as usize, i2 as usize)
                    - yy * self.mypoints.value(self.mylastp as usize, i2 as usize);
                let b_iy = self.mypoints.value(i as usize, i21 as usize)
                    - xx * self.mypoints.value(self.myfirstp as usize, i21 as usize)
                    - yy * self.mypoints.value(self.mylastp as usize, i21 as usize);
                let b_iz = self.mypoints.value(i as usize, i22 as usize)
                    - xx * self.mypoints.value(self.myfirstp as usize, i22 as usize)
                    - yy * self.mypoints.value(self.mylastp as usize, i22 as usize);
                my_b.set_value(ix as usize, b_ix);
                my_b.set_value(iy as usize, b_iy);
                my_b.set_value(iz as usize, b_iz);
                i2 += 3;
                nrow += 3 * neq;
            }

            for _ci in 1..=self.nbp2d {
                let i21 = i2 + 1;
                let ix = i + nrow;
                let iy = ix + neq;
                if first_t {
                    my_v1.set_value(ix as usize, ai2 * self.vec1t.value(i2 as usize));
                    my_v1.set_value(iy as usize, ai2 * self.vec1t.value(i21 as usize));
                }
                if last_t {
                    my_v2.set_value(ix as usize, -aid * self.vec2t.value(i2 as usize));
                    my_v2.set_value(iy as usize, -aid * self.vec2t.value(i21 as usize));
                }
                let b_ix = self.mypoints.value(i as usize, i2 as usize)
                    - xx * self.mypoints.value(self.myfirstp as usize, i2 as usize)
                    - yy * self.mypoints.value(self.mylastp as usize, i2 as usize);
                let b_iy = self.mypoints.value(i as usize, i21 as usize)
                    - xx * self.mypoints.value(self.myfirstp as usize, i21 as usize)
                    - yy * self.mypoints.value(self.mylastp as usize, i21 as usize);
                my_b.set_value(ix as usize, b_ix);
                my_b.set_value(iy as usize, b_iy);
                nrow += 2 * neq;
                i2 += 2;
            }
        }

        // Construction de TA*A et TA*B:
        for k in self.first_p..=self.last_p {
            let indexdeb = self.myindex.value(k as usize) + 1;
            let indexfin = indexdeb + self.deg;
            let jinit = indexdeb.max(self.resinit);
            let jfin = indexfin.min(self.resfin);
            let k1 = k + myfirst - self.first_p;
            for i in 0..=na1 {
                let nb = (i * neq + k1) as usize;
                let v1 = if first_t { my_v1.value(nb) } else { 0.0 };
                let v2 = if last_t { my_v2.value(nb) } else { 0.0 };
                let b = my_b.value(nb);
                let inc = i * nincx - self.resinit + 1;
                for j in jinit..=jfin {
                    let akj = self.a.value(k as usize, j as usize);
                    let u = (j + inc) as usize;
                    if first_t {
                        let cur = the_v1.value(u);
                        the_v1.set_value(u, cur + akj * v1);
                    }
                    if last_t {
                        let cur = the_v2.value(u);
                        the_v2.set_value(u, cur + akj * v2);
                    }
                    let cur = my_tab.value(u);
                    my_tab.set_value(u, cur + akj * b);
                }
                if first_t {
                    taf1 += v1 * v1;
                    tab1 += v1 * b;
                }
                if last_t {
                    taf2 += v2 * v2;
                    tab2 += v2 * b;
                }
                if first_t && last_t {
                    taf3 += v1 * v2;
                }
            }
        }

        if first_t {
            the_v1.set_value(ninc1 as usize, taf1);
            my_tab.set_value(ninc1 as usize, tab1);
        }
        if last_t {
            the_v2.set_value(self.ninc as usize, taf2);
            my_tab.set_value(self.ninc as usize, tab2);
        }
        if first_t && last_t {
            the_v2.set_value(ninc1 as usize, taf3);
        }

        if self.resinit <= self.resfin {
            let mut index = MathIntVector::new(1, nincx as usize);
            self.search_index(&mut index);
            let mut aa = MathVector::new(1, index.value(nincx as usize) as usize);
            self.make_taa_packed(&mut aa);

            let mut kk = 1;
            for _k in 1..=self.na {
                for i in 1..=aa.len() {
                    the_a.set_value(kk, aa.value(i));
                    kk += 1;
                }
            }
        }

        let length = the_a.len() as i32;

        if first_t && last_t {
            for j in 1..=ninc1 {
                the_a.set_value(
                    (length - 2 * self.ninc + j + 1) as usize,
                    the_v1.value(j as usize),
                );
            }
            for j in 1..=self.ninc {
                the_a.set_value(
                    (length - self.ninc + j) as usize,
                    the_v2.value(j as usize),
                );
            }
        } else if first_t {
            for j in 1..=self.ninc {
                the_a.set_value(
                    (length - self.ninc + j) as usize,
                    the_v1.value(j as usize),
                );
            }
        } else if last_t {
            for j in 1..=self.ninc {
                the_a.set_value(
                    (length - self.ninc + j) as usize,
                    the_v2.value(j as usize),
                );
            }
        }
    }

    /// `MakeTAA(AA, myTAB)` (`gxx:1706-1764`): the unconstrained block and the
    /// right-hand side.
    pub(crate) fn make_taa_matrix(&self, aa: &mut MathVector, my_tab: &mut MathMatrix) {
        let mut the_a = MathMatrix::new(
            self.resinit.max(0) as usize,
            self.resfin.max(0) as usize,
            self.resinit.max(0) as usize,
            self.resfin.max(0) as usize,
        );
        the_a.init(0.0);

        for k in self.first_p..=self.last_p {
            let indexdeb = self.myindex.value(k as usize) + 1;
            let indexfin = indexdeb + self.deg;
            let jinit = indexdeb.max(self.resinit);
            let jfin = indexfin.min(self.resfin);
            for i in jinit..=jfin {
                let akj = self.a.value(k as usize, i as usize);
                for j in jinit..=i {
                    let cur = the_a.value(i as usize, j as usize);
                    the_a.set_value(
                        i as usize,
                        j as usize,
                        cur + self.a.value(k as usize, j as usize) * akj,
                    );
                }
                for j in 1..=self.b2.col_count() as i32 {
                    let cur = my_tab.value(i as usize, j as usize);
                    let b2 = self.b2.value(k as usize, j as usize);
                    my_tab.set_value(i as usize, j as usize, cur + akj * b2);
                }
            }
        }

        let len = match &self.myknots {
            Some(k) => k.len() as i32,
            None => 2,
        };
        let mults = self.mymults.clone().unwrap_or_default();
        let mut iinit = self.resinit;
        let mut jinit = self.resinit;
        let mut ifin = (self.deg + 1).min(self.resfin);
        let mut i2 = 1;

        for k in 2..=len {
            for i in iinit..=ifin {
                for j in jinit..=i {
                    aa.set_value(i2, the_a.value(i as usize, j as usize));
                    i2 += 1;
                }
            }
            if self.mymults.is_some() {
                iinit = ifin + 1;
                let d = ifin + mults[(k - 1) as usize];
                ifin = d.min(self.resfin);
                jinit = (d - self.deg).max(self.resinit);
            }
        }
    }

    /// `MakeTAA(AA)` (`gxx:1766-1812`): the packed unconstrained block alone.
    pub(crate) fn make_taa_packed(&self, aa: &mut MathVector) {
        let mut the_a = MathMatrix::new(
            self.resinit.max(0) as usize,
            self.resfin.max(0) as usize,
            self.resinit.max(0) as usize,
            self.resfin.max(0) as usize,
        );

        for k in self.first_p..=self.last_p {
            let indexdeb = self.myindex.value(k as usize) + 1;
            let indexfin = indexdeb + self.deg;
            let jinit = indexdeb.max(self.resinit);
            let jfin = indexfin.min(self.resfin);
            for i in jinit..=jfin {
                let akj = self.a.value(k as usize, i as usize);
                for j in jinit..=i {
                    let cur = the_a.value(i as usize, j as usize);
                    the_a.set_value(
                        i as usize,
                        j as usize,
                        cur + self.a.value(k as usize, j as usize) * akj,
                    );
                }
            }
        }

        let len = match &self.myknots {
            Some(k) => k.len() as i32,
            None => 2,
        };
        let mults = self.mymults.clone().unwrap_or_default();
        let mut i2 = 1;
        let mut iinit = self.resinit;
        let mut jinit = self.resinit;
        let mut ifin = (self.deg + 1).min(self.resfin);

        for k in 2..=len {
            for i in iinit..=ifin {
                for j in jinit..=i {
                    aa.set_value(i2, the_a.value(i as usize, j as usize));
                    i2 += 1;
                }
            }
            if self.mymults.is_some() {
                iinit = ifin + 1;
                let d = ifin + mults[(k - 1) as usize];
                ifin = d.min(self.resfin);
                jinit = (d - self.deg).max(self.resinit);
            }
        }
    }
}

