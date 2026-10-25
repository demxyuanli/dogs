//! `IsSolution` / `Decroch` of `BRepBlend_SurfRstConstRad`.
//!
//! Source: `ModelingAlgorithms/TKFillet/BRepBlend/BRepBlend_SurfRstConstRad.cxx:311-427`
//! (`IsSolution`), `:523-575` (`Decroch`) (OCCT 8.0.0).

use occt_core::gp::{GpPnt, GpVec};
use occt_geom::curve::Curve;
use occt_math::{svd_padded, svd_solve_occt, Gauss, MathMatrix, MathVector};

use super::equations::pnt_vec;
use super::SurfRstConstRad;

impl SurfRstConstRad {
    /// `IsSolution(Sol, Tol)` (`cxx:311-427`): the three values are within
    /// tolerance, then the tangents are solved for from the jacobian
    /// (`math_Gauss`, falling back to `math_SVD`), and `maxang` / `minang` /
    /// `distmin` are updated from the section angle.
    pub(super) fn check_solution(&mut self, sol: &MathVector, tol: f64) -> bool {
        let mut valsol = MathVector::new(1, 3);
        let mut gradsol = MathMatrix::new(1, 3, 1, 3);
        self.compute_values(sol, &mut valsol, &mut gradsol);

        if !(valsol.value(1).abs() <= tol
            && valsol.value(2).abs() <= tol
            && valsol.value(3).abs() <= 2.0 * tol * self.ray.abs())
        {
            self.istangent = true;
            return false;
        }

        // Calculation of tangents (`cxx:325-389`).
        self.pt2ds = occt_core::gp::GpPnt2d::new(sol.value(1), sol.value(2));
        self.prmrst = sol.value(3);
        self.pt2drst = self.rst.d0(self.prmrst);
        let (pts, d1u1, d1v1) = self.surf.d1(sol.value(1), sol.value(2));
        self.pts = pts;
        let (ptrst, d1) = self.cons.d1(sol.value(3));
        self.ptrst = ptrst;

        let dnplan = GpVec::zero().set_linear_form_2(
            1.0 / self.normtg,
            &self.d2gui,
            -1.0 / self.normtg * self.nplan.dot(&self.d2gui),
            &self.nplan,
        );

        let mut secmember = MathVector::new(1, 3);
        let temp = pnt_vec(&self.pts).subtracted(&pnt_vec(&self.ptgui));
        secmember.set_value(1, self.normtg - dnplan.dot(&temp));
        let temp = pnt_vec(&self.ptrst).subtracted(&pnt_vec(&self.ptgui));
        secmember.set_value(2, self.normtg - dnplan.dot(&temp));

        let ns = d1u1.crossed(&d1v1);
        let ncrossns = self.nplan.crossed(&ns);
        let ndotns = self.nplan.dot(&ns);
        let norm = ncrossns.magnitude();

        let grosterme = ncrossns.dot(&dnplan.crossed(&ns)) / norm / norm;
        let temp = GpVec::zero().set_linear_form_3(
            self.ray / norm * (dnplan.dot(&ns) - grosterme * ndotns),
            &self.nplan,
            self.ray * ndotns / norm,
            &dnplan,
            self.ray * grosterme / norm,
            &ns,
        );

        let ns = GpVec::zero().set_linear_form_2(
            ndotns / norm,
            &self.nplan,
            -1.0 / norm,
            &ns,
        );
        let resul = GpVec::zero().set_linear_form_2(
            self.ray,
            &ns,
            1.0,
            &GpVec::from_pnts(&self.ptrst, &self.pts),
        );
        secmember.set_value(3, -2.0 * temp.dot(&resul));

        let resol = Gauss::new(&gradsol, super::GAUSS_DEFAULT_EPS);
        if resol.is_done() {
            if let Some(x) = resol.solve(&secmember) {
                secmember = x;
            }
            self.istangent = false;
        } else {
            let sing_rs = svd_padded(&gradsol);
            if sing_rs.is_done {
                let dedt = secmember.clone();
                svd_solve_occt(&sing_rs, &dedt, &mut secmember, super::SVD_EPS);
                self.istangent = false;
            } else {
                self.istangent = true;
            }
        }

        if !self.istangent {
            self.tgs = GpVec::zero().set_linear_form_2(
                secmember.value(1),
                &d1u1,
                secmember.value(2),
                &d1v1,
            );
            self.tgrst = d1.multiplied_scalar(secmember.value(3));
            self.tg2ds = occt_core::gp::GpVec2d::new(secmember.value(1), secmember.value(2));
            let (_, d1urst, d1vrst) = self.surfrst.d1(self.pt2drst.x(), self.pt2drst.y());
            let mut a = 0.0;
            let mut b = 0.0;
            super::equations::t3dto2d(&mut a, &mut b, &self.tgrst, &d1urst, &d1vrst);
            self.tg2drst = occt_core::gp::GpVec2d::new(a, b);
        }

        // Update of maxang (`cxx:392-421`). `ns` may be reversed, exactly as
        // OCCT reverses its local copy before building `ns2`.
        let mut ns = ns;
        if self.ray > 0.0 {
            ns.reverse();
        }
        let ns2 = resul.normalized().reversed();

        let cosa = ns.dot(&ns2);
        let mut sina = self.nplan.dot(&ns.crossed(&ns2));
        if self.choix % 2 != 0 {
            // `// nplan is changed to -nplan` (`cxx:405`).
            sina = -sina;
        }

        let mut angle = cosa.acos();
        if sina < 0.0 {
            angle = 2.0 * std::f64::consts::PI - angle;
        }

        if angle > self.maxang {
            self.maxang = angle;
        }
        if angle < self.minang {
            self.minang = angle;
        }
        self.distmin = self.distmin.min(self.pts.distance(&self.ptrst));

        true
    }

    /// `Decroch(Sol, NS, TgS)` (`cxx:523-575`): a section leaves the
    /// restriction when the restriction normal inside the section plane has a
    /// positive projection onto the restriction tangent (`dot < 1e-10`).
    ///
    /// `surfref` / `rstref` are the face and the pcurve `Set(SurfRef, RstRef)`
    /// stored. OCCT dereferences those two handles without a null check
    /// (`cxx:547-549`); the port's `Option`s are unwrapped, so a missing
    /// `Set` panics exactly where OCCT would fault.
    pub(super) fn compute_decroch(&self, sol: &MathVector, ns_out: &mut GpVec, tgs_out: &mut GpVec) -> bool {
        let (bid, d1u, d1v) = self.surf.d1(sol.value(1), sol.value(2));
        let ns = d1u.crossed(&d1v);
        *ns_out = ns;

        let norm = self.nplan.crossed(&ns).magnitude();
        let unsurnorm = 1.0 / norm;
        let mut ns_in_plane =
            GpVec::zero().set_linear_form_2(self.nplan.dot(&ns) * unsurnorm, &self.nplan, -unsurnorm, &ns);

        let center = GpPnt::from_xyz(
            &bid.coord
                .added(&ns_in_plane.coord.multiply(self.ray)),
        );
        if self.choix > 2 {
            ns_in_plane.reverse();
        }
        let mut tgs = self.nplan.crossed(&GpVec::from_pnts(&center, &bid));
        if self.choix % 2 == 1 {
            tgs.reverse();
        }
        *tgs_out = tgs;

        let rstref = self
            .rstref
            .as_ref()
            .expect("BRepBlend_SurfRstConstRad::Decroch: Set(SurfRef, RstRef) not called");
        let surfref = self
            .surfref
            .as_ref()
            .expect("BRepBlend_SurfRstConstRad::Decroch: Set(SurfRef, RstRef) not called");
        let uv = rstref.d0(sol.value(3));
        let (bid, d1u, d1v) = surfref.d1(uv.x(), uv.y());

        let nrst = d1u.crossed(&d1v);
        let norm = self.nplan.crossed(&nrst).magnitude();
        let unsurnorm = 1.0 / norm;
        let mut nrst_in_plane =
            GpVec::zero().set_linear_form_2(self.nplan.dot(&nrst) * unsurnorm, &self.nplan, -unsurnorm, &nrst);
        let centptrst = GpVec::from_pnts(&center, &bid);
        if centptrst.dot(&nrst_in_plane) < 0.0 {
            nrst_in_plane.reverse();
        }
        let mut tgrst = self.nplan.crossed(&centptrst);
        if self.choix % 2 == 1 {
            tgrst.reverse();
        }

        let mut nt = nrst_in_plane.magnitude();
        nt *= tgrst.magnitude();
        if nt.abs() < super::DECROCH_NT_EPS {
            // Singularity or Incoherence (`cxx:566-569`).
            return false;
        }
        let mut dot = nrst_in_plane.dot(&tgrst);
        dot /= nt;

        dot < super::DECROCH_DOT_EPS
    }
}
