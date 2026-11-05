//! `BlendFunc_ConstRad::ComputeValues` — the `Order == 2` arm.
//!
//! Source: `ModelingAlgorithms/TKFillet/BlendFunc/BlendFunc_ConstRad.cxx:434-774`
//! (OCCT 8.0.0). The prologue of `ComputeValues` and its order-0 / order-1
//! arms are in [`super::compute_values`]; this file fills `D2EDX2`, and (under
//! `ByParam`) `D2EDXDT` and `D2EDT2`.

use occt_core::gp::GpVec;

use super::compute_values::SectionPlaneTerms;
use super::ConstRad;

impl ConstRad {
    /// `cxx:434-774`: the second derivatives of the four equations with respect
    /// to `X`, then (when `by_param`) with respect to `(X, t)` and to `t`.
    pub(super) fn compute_second_order(&mut self, terms: &SectionPlaneTerms, by_param: bool) {
        // `D2EDX2.Init(0)` (`cxx:440`).
        self.d2edx2.init(0.0);

        // ================ the second derivative of equation 1 (`cxx:442-448`).
        self.d2edx2
            .set_value(1, 1, 1, self.nplan.dot(&self.d2u1) / 2.0);
        let mix = self.nplan.dot(&self.d2uv1) / 2.0;
        self.d2edx2.set_value(1, 2, 1, mix);
        self.d2edx2.set_value(1, 1, 2, mix);
        self.d2edx2
            .set_value(1, 2, 2, self.nplan.dot(&self.d2v1) / 2.0);

        self.d2edx2
            .set_value(1, 3, 3, self.nplan.dot(&self.d2u2) / 2.0);
        let mix = self.nplan.dot(&self.d2uv2) / 2.0;
        self.d2edx2.set_value(1, 4, 3, mix);
        self.d2edx2.set_value(1, 3, 4, mix);
        self.d2edx2
            .set_value(1, 4, 4, self.nplan.dot(&self.d2v2) / 2.0);

        // ================
        // ==  Surface 1 ==
        // ================
        let mut carre = terms.invnorm1 * terms.invnorm1;
        let mut cube = carre * terms.invnorm1;
        // Derived double compared to u1 (`cxx:453-477`).
        let d2ns1u1 = GpVec::zero().set_linear_form_3(
            1.0,
            &self.d3u1.crossed(&self.d1v1),
            2.0,
            &self.d2u1.crossed(&self.d2uv1),
            1.0,
            &self.d1u1.crossed(&self.d3uuv1),
        );
        let dprim = terms
            .ncrossns1
            .dot(&self.nplan.crossed(&self.dns1u1));
        let smallterm = -2.0 * dprim * cube;
        let dsecn = terms
            .ncrossns1
            .dot(&self.nplan.crossed(&d2ns1u1))
            + self.nplan.crossed(&self.dns1u1).square_magnitude();
        let grosterme = (3.0 * dprim * dprim * carre - dsecn) * cube;

        let temp =
            GpVec::zero().set_linear_form_2(grosterme * terms.ndotns1, &self.nplan, -grosterme, &self.nsurf1);
        let p1 = self.nplan.dot(&self.dns1u1);
        let p2 = self.nplan.dot(&d2ns1u1);
        self.d2ndu1 = GpVec::zero().set_linear_form_3(
            terms.invnorm1 * p2 + smallterm * p1,
            &self.nplan,
            -smallterm,
            &self.dns1u1,
            -terms.invnorm1,
            &d2ns1u1,
        );
        self.d2ndu1 = self.d2ndu1.added(&temp);
        let resul = self.d2ndu1.multiplied_scalar(self.ray1).added(&self.d2u1);
        self.d2edx2.set_value(2, 1, 1, resul.x());
        self.d2edx2.set_value(3, 1, 1, resul.y());
        self.d2edx2.set_value(4, 1, 1, resul.z());

        // Derived double compared to u1, v1 (`cxx:478-503`).
        let d2ns1uv1 = self
            .d3uuv1
            .crossed(&self.d1v1)
            .added(&self.d2u1.crossed(&self.d2v1))
            .added(&self.d1u1.crossed(&self.d3uvv1));
        let mut uterm = terms
            .ncrossns1
            .dot(&self.nplan.crossed(&self.dns1u1));
        let mut vterm = terms
            .ncrossns1
            .dot(&self.nplan.crossed(&self.dns1v1));
        let dsecn = self
            .nplan
            .crossed(&self.dns1v1)
            .dot(&self.nplan.crossed(&self.dns1u1))
            + terms
                .ncrossns1
                .dot(&self.nplan.crossed(&d2ns1uv1));
        let grosterme = (3.0 * uterm * vterm * carre - dsecn) * cube;
        uterm *= -cube; // and only now (`cxx:486`)
        vterm *= -cube;

        let p1 = self.nplan.dot(&self.dns1u1);
        let p2 = self.nplan.dot(&self.dns1v1);
        let temp = GpVec::zero().set_linear_form_3(
            grosterme * terms.ndotns1,
            &self.nplan,
            -grosterme,
            &self.nsurf1,
            -terms.invnorm1,
            &d2ns1uv1,
        );
        self.d2nduv1 = GpVec::zero().set_linear_form_3(
            terms.invnorm1 * self.nplan.dot(&d2ns1uv1) + uterm * p2 + vterm * p1,
            &self.nplan,
            -uterm,
            &self.dns1v1,
            -vterm,
            &self.dns1u1,
        );
        self.d2nduv1 = self.d2nduv1.added(&temp);
        let resul = self
            .d2nduv1
            .multiplied_scalar(self.ray1)
            .added(&self.d2uv1);
        self.d2edx2.set_value(2, 2, 1, resul.x());
        self.d2edx2.set_value(2, 1, 2, resul.x());
        self.d2edx2.set_value(3, 2, 1, resul.y());
        self.d2edx2.set_value(3, 1, 2, resul.y());
        self.d2edx2.set_value(4, 2, 1, resul.z());
        self.d2edx2.set_value(4, 1, 2, resul.z());

        // Derived double compared to v1 (`cxx:505-526`).
        let d2ns1v1 = GpVec::zero().set_linear_form_3(
            1.0,
            &self.d1u1.crossed(&self.d3v1),
            2.0,
            &self.d2uv1.crossed(&self.d2v1),
            1.0,
            &self.d3uvv1.crossed(&self.d1v1),
        );
        let dprim = terms
            .ncrossns1
            .dot(&self.nplan.crossed(&self.dns1v1));
        let smallterm = -2.0 * dprim * cube;
        let dsecn = terms
            .ncrossns1
            .dot(&self.nplan.crossed(&d2ns1v1))
            + self.nplan.crossed(&self.dns1v1).square_magnitude();
        let grosterme = (3.0 * dprim * dprim * carre - dsecn) * cube;

        let p1 = self.nplan.dot(&self.dns1v1);
        let p2 = self.nplan.dot(&d2ns1v1);
        let temp =
            GpVec::zero().set_linear_form_2(grosterme * terms.ndotns1, &self.nplan, -grosterme, &self.nsurf1);
        self.d2ndv1 = GpVec::zero().set_linear_form_3(
            terms.invnorm1 * p2 + smallterm * p1,
            &self.nplan,
            -smallterm,
            &self.dns1v1,
            -terms.invnorm1,
            &d2ns1v1,
        );
        self.d2ndv1 = self.d2ndv1.added(&temp);
        let resul = self.d2ndv1.multiplied_scalar(self.ray1).added(&self.d2v1);
        self.d2edx2.set_value(2, 2, 2, resul.x());
        self.d2edx2.set_value(3, 2, 2, resul.y());
        self.d2edx2.set_value(4, 2, 2, resul.z());

        // ================
        // ==  Surface 2 ==
        // ================
        carre = terms.invnorm2 * terms.invnorm2;
        cube = carre * terms.invnorm2;
        // Derived double compared to u2 (`cxx:534-554`).
        let d2ns1u2 = GpVec::zero().set_linear_form_3(
            1.0,
            &self.d3u2.crossed(&self.d1v2),
            2.0,
            &self.d2u2.crossed(&self.d2uv2),
            1.0,
            &self.d1u2.crossed(&self.d3uuv2),
        );
        let dprim = terms
            .ncrossns2
            .dot(&self.nplan.crossed(&self.dns1u2));
        let smallterm = -2.0 * dprim * cube;
        let dsecn = terms
            .ncrossns2
            .dot(&self.nplan.crossed(&d2ns1u2))
            + self.nplan.crossed(&self.dns1u2).square_magnitude();
        let grosterme = (3.0 * dprim * dprim * carre - dsecn) * cube;

        let temp =
            GpVec::zero().set_linear_form_2(grosterme * terms.ndotns2, &self.nplan, -grosterme, &self.nsurf2);
        let p1 = self.nplan.dot(&self.dns1u2);
        let p2 = self.nplan.dot(&d2ns1u2);
        self.d2ndu2 = GpVec::zero().set_linear_form_3(
            terms.invnorm2 * p2 + smallterm * p1,
            &self.nplan,
            -smallterm,
            &self.dns1u2,
            -terms.invnorm2,
            &d2ns1u2,
        );
        self.d2ndu2 = self.d2ndu2.added(&temp);
        let resul = GpVec::zero().set_linear_form_2(-self.ray2, &self.d2ndu2, -1.0, &self.d2u2);
        self.d2edx2.set_value(2, 3, 3, resul.x());
        self.d2edx2.set_value(3, 3, 3, resul.y());
        self.d2edx2.set_value(4, 3, 3, resul.z());

        // Derived double compared to u2, v2 (`cxx:556-583`).
        let d2ns1uv2 = self
            .d3uuv2
            .crossed(&self.d1v2)
            .added(&self.d2u2.crossed(&self.d2v2))
            .added(&self.d1u2.crossed(&self.d3uvv2));
        let mut uterm = terms
            .ncrossns2
            .dot(&self.nplan.crossed(&self.dns1u2));
        let mut vterm = terms
            .ncrossns2
            .dot(&self.nplan.crossed(&self.dns1v2));
        let dsecn = self
            .nplan
            .crossed(&self.dns1v2)
            .dot(&self.nplan.crossed(&self.dns1u2))
            + terms
                .ncrossns2
                .dot(&self.nplan.crossed(&d2ns1uv2));
        let grosterme = (3.0 * uterm * vterm * carre - dsecn) * cube;
        uterm *= -cube; // and only now (`cxx:566`)
        vterm *= -cube;

        let p1 = self.nplan.dot(&self.dns1u2);
        let p2 = self.nplan.dot(&self.dns1v2);
        let temp = GpVec::zero().set_linear_form_3(
            grosterme * terms.ndotns2,
            &self.nplan,
            -grosterme,
            &self.nsurf2,
            -terms.invnorm2,
            &d2ns1uv2,
        );
        self.d2nduv2 = GpVec::zero().set_linear_form_3(
            terms.invnorm2 * self.nplan.dot(&d2ns1uv2) + uterm * p2 + vterm * p1,
            &self.nplan,
            -uterm,
            &self.dns1v2,
            -vterm,
            &self.dns1u2,
        );
        self.d2nduv2 = self.d2nduv2.added(&temp);
        let resul = GpVec::zero().set_linear_form_2(-self.ray2, &self.d2nduv2, -1.0, &self.d2uv2);
        self.d2edx2.set_value(2, 4, 3, resul.x());
        self.d2edx2.set_value(2, 3, 4, resul.x());
        self.d2edx2.set_value(3, 4, 3, resul.y());
        self.d2edx2.set_value(3, 3, 4, resul.y());
        self.d2edx2.set_value(4, 4, 3, resul.z());
        self.d2edx2.set_value(4, 3, 4, resul.z());

        // Derived double compared to v2 (`cxx:585-606`).
        let d2ns1v2 = GpVec::zero().set_linear_form_3(
            1.0,
            &self.d1u2.crossed(&self.d3v2),
            2.0,
            &self.d2uv2.crossed(&self.d2v2),
            1.0,
            &self.d3uvv2.crossed(&self.d1v2),
        );
        let dprim = terms
            .ncrossns2
            .dot(&self.nplan.crossed(&self.dns1v2));
        let smallterm = -2.0 * dprim * cube;
        let dsecn = terms
            .ncrossns2
            .dot(&self.nplan.crossed(&d2ns1v2))
            + self.nplan.crossed(&self.dns1v2).square_magnitude();
        let grosterme = (3.0 * dprim * dprim * carre - dsecn) * cube;

        let p1 = self.nplan.dot(&self.dns1v2);
        let p2 = self.nplan.dot(&d2ns1v2);
        let temp =
            GpVec::zero().set_linear_form_2(grosterme * terms.ndotns2, &self.nplan, -grosterme, &self.nsurf2);
        self.d2ndv2 = GpVec::zero().set_linear_form_3(
            terms.invnorm2 * p2 + smallterm * p1,
            &self.nplan,
            -smallterm,
            &self.dns1v2,
            -terms.invnorm2,
            &d2ns1v2,
        );
        self.d2ndv2 = self.d2ndv2.added(&temp);
        let resul = GpVec::zero().set_linear_form_2(-self.ray2, &self.d2ndv2, -1.0, &self.d2v2);
        self.d2edx2.set_value(2, 4, 4, resul.x());
        self.d2edx2.set_value(3, 4, 4, resul.y());
        self.d2edx2.set_value(4, 4, 4, resul.z());

        if by_param {
            self.compute_second_order_in_t(terms);
        }
    }

    /// `cxx:608-774`: `D2EDXDT` and `D2EDT2` when the guide parameter is
    /// explicit.
    fn compute_second_order_in_t(&mut self, terms: &SectionPlaneTerms) {
        // ---------- Derivation double in t, X ------------------------------.
        self.d2edxdt
            .set_value(1, 1, self.dnplan.dot(&self.d1u1) / 2.0);
        self.d2edxdt
            .set_value(1, 2, self.dnplan.dot(&self.d1v1) / 2.0);
        self.d2edxdt
            .set_value(1, 3, self.dnplan.dot(&self.d1u2) / 2.0);
        self.d2edxdt
            .set_value(1, 4, self.dnplan.dot(&self.d1v2) / 2.0);

        let mut carre = terms.invnorm1 * terms.invnorm1;
        let mut cube = carre * terms.invnorm1;
        // --> Derived compared to u1 and t (`cxx:619-641`).
        let tterm = terms
            .ncrossns1
            .dot(&self.dnplan.crossed(&self.nsurf1));
        let smallterm = -tterm * cube;
        // Derived from the norm (`cxx:624`).
        let mut uterm = terms
            .ncrossns1
            .dot(&self.nplan.crossed(&self.dns1u1));
        let dsecn = self
            .nplan
            .crossed(&self.dns1u1)
            .dot(&self.dnplan.crossed(&self.nsurf1))
            + terms
                .ncrossns1
                .dot(&self.dnplan.crossed(&self.dns1u1));
        let grosterme = (3.0 * uterm * tterm * carre - dsecn) * cube;
        uterm *= -cube;

        let p1 = self.dnplan.dot(&self.nsurf1);
        let p2 = self.nplan.dot(&self.dns1u1);
        let p12 = self.dnplan.dot(&self.dns1u1);

        self.d2ndtu1 = GpVec::zero().set_linear_form_3(
            terms.invnorm1 * p12 + smallterm * p2 + uterm * p1 + grosterme * terms.ndotns1,
            &self.nplan,
            terms.invnorm1 * p2 + uterm * terms.ndotns1,
            &self.dnplan,
            -smallterm,
            &self.dns1u1,
        );
        self.d2ndtu1 = self
            .d2ndtu1
            .subtracted(&self.nsurf1.multiplied_scalar(grosterme));

        let resul = self.d2ndtu1.multiplied_scalar(self.ray1);
        self.d2edxdt.set_value(2, 1, resul.x());
        self.d2edxdt.set_value(3, 1, resul.y());
        self.d2edxdt.set_value(4, 1, resul.z());

        // --> Derived compared to v1 and t (`cxx:644-668`).
        // Derived from the norm (`cxx:648`).
        uterm = terms
            .ncrossns1
            .dot(&self.nplan.crossed(&self.dns1v1));
        let dsecn = self
            .nplan
            .crossed(&self.dns1v1)
            .dot(&self.dnplan.crossed(&self.nsurf1))
            + terms
                .ncrossns1
                .dot(&self.dnplan.crossed(&self.dns1v1));
        let grosterme = (3.0 * uterm * tterm * carre - dsecn) * cube;
        uterm *= -cube;

        let p1 = self.dnplan.dot(&self.nsurf1);
        let p2 = self.nplan.dot(&self.dns1v1);
        let p12 = self.dnplan.dot(&self.dns1v1);
        self.d2ndtv1 = GpVec::zero().set_linear_form_3(
            terms.invnorm1 * p12 + uterm * p1 + smallterm * p2 + grosterme * terms.ndotns1,
            &self.nplan,
            terms.invnorm1 * p2 + uterm * terms.ndotns1,
            &self.dnplan,
            -smallterm,
            &self.dns1v1,
        );
        self.d2ndtv1 = self
            .d2ndtv1
            .subtracted(&self.nsurf1.multiplied_scalar(grosterme));

        let resul = self.d2ndtv1.multiplied_scalar(self.ray1);
        self.d2edxdt.set_value(2, 2, resul.x());
        self.d2edxdt.set_value(3, 2, resul.y());
        self.d2edxdt.set_value(4, 2, resul.z());

        carre = terms.invnorm2 * terms.invnorm2;
        cube = carre * terms.invnorm2;
        // --> Derived compared to u2 and t (`cxx:672-697`).
        let tterm = terms
            .ncrossns2
            .dot(&self.dnplan.crossed(&self.nsurf2));
        let smallterm = -tterm * cube;
        // Derived from the norm (`cxx:677`).
        let mut uterm = terms
            .ncrossns2
            .dot(&self.nplan.crossed(&self.dns1u2));
        let dsecn = self
            .nplan
            .crossed(&self.dns1u2)
            .dot(&self.dnplan.crossed(&self.nsurf2))
            + terms
                .ncrossns2
                .dot(&self.dnplan.crossed(&self.dns1u2));
        let grosterme = (3.0 * uterm * tterm * carre - dsecn) * cube;
        uterm *= -cube;

        let p1 = self.dnplan.dot(&self.nsurf2);
        let p2 = self.nplan.dot(&self.dns1u2);
        let p12 = self.dnplan.dot(&self.dns1u2);

        self.d2ndtu2 = GpVec::zero().set_linear_form_3(
            terms.invnorm2 * p12 + smallterm * p2 + uterm * p1 + grosterme * terms.ndotns2,
            &self.nplan,
            terms.invnorm2 * p2 + uterm * terms.ndotns2,
            &self.dnplan,
            -smallterm,
            &self.dns1u2,
        );
        self.d2ndtu2 = self
            .d2ndtu2
            .subtracted(&self.nsurf2.multiplied_scalar(grosterme));

        let resul = self.d2ndtu2.multiplied_scalar(-self.ray2);
        self.d2edxdt.set_value(2, 3, resul.x());
        self.d2edxdt.set_value(3, 3, resul.y());
        self.d2edxdt.set_value(4, 3, resul.z());

        // --> Derived compared to v2 and t (`cxx:700-723`).
        // Derived from the norm (`cxx:702`).
        uterm = terms
            .ncrossns2
            .dot(&self.nplan.crossed(&self.dns1v2));
        let dsecn = self
            .nplan
            .crossed(&self.dns1v2)
            .dot(&self.dnplan.crossed(&self.nsurf2))
            + terms
                .ncrossns2
                .dot(&self.dnplan.crossed(&self.dns1v2));
        let grosterme = (3.0 * uterm * tterm * carre - dsecn) * cube;
        uterm *= -cube;

        let p1 = self.dnplan.dot(&self.nsurf2);
        let p2 = self.nplan.dot(&self.dns1v2);
        let p12 = self.dnplan.dot(&self.dns1v2);

        self.d2ndtv2 = GpVec::zero().set_linear_form_3(
            terms.invnorm2 * p12 + smallterm * p2 + uterm * p1 + grosterme * terms.ndotns2,
            &self.nplan,
            terms.invnorm2 * p2 + uterm * terms.ndotns2,
            &self.dnplan,
            -smallterm,
            &self.dns1v2,
        );
        self.d2ndtv2 = self
            .d2ndtv2
            .subtracted(&self.nsurf2.multiplied_scalar(grosterme));

        let resul = self.d2ndtv2.multiplied_scalar(-self.ray2);
        self.d2edxdt.set_value(2, 4, resul.x());
        self.d2edxdt.set_value(3, 4, resul.y());
        self.d2edxdt.set_value(4, 4, resul.z());

        // ---------- Derivation double in t ----------------------------------.
        // Derived from n1 compared to w (`cxx:726-744`).
        carre = terms.invnorm1 * terms.invnorm1;
        cube = carre * terms.invnorm1;
        // Derived from the norm (`cxx:729`).
        let dprim = terms
            .ncrossns1
            .dot(&self.dnplan.crossed(&self.nsurf1));
        let smallterm = -2.0 * dprim * cube;
        let dsecn = self.dnplan.crossed(&self.nsurf1).square_magnitude()
            + terms
                .ncrossns1
                .dot(&self.d2nplan.crossed(&self.nsurf1));
        let grosterme = (3.0 * dprim * dprim * carre - dsecn) * cube;

        let p1 = self.dnplan.dot(&self.nsurf1);
        let p2 = self.d2nplan.dot(&self.nsurf1);

        let temp =
            GpVec::zero().set_linear_form_2(grosterme * terms.ndotns1, &self.nplan, -grosterme, &self.nsurf1);
        self.d2n1w = GpVec::zero().set_linear_form_3(
            smallterm * p1 + terms.invnorm1 * p2,
            &self.nplan,
            smallterm * terms.ndotns1 + 2.0 * terms.invnorm1 * p1,
            &self.dnplan,
            terms.ndotns1 * terms.invnorm1,
            &self.d2nplan,
        );
        self.d2n1w = self.d2n1w.added(&temp);

        // Derived from n2 compared to w (`cxx:746-764`).
        carre = terms.invnorm2 * terms.invnorm2;
        cube = carre * terms.invnorm2;
        // Derived from the norm (`cxx:749`).
        let dprim = terms
            .ncrossns2
            .dot(&self.dnplan.crossed(&self.nsurf2));
        let smallterm = -2.0 * dprim * cube;
        let dsecn = self.dnplan.crossed(&self.nsurf2).square_magnitude()
            + terms
                .ncrossns2
                .dot(&self.d2nplan.crossed(&self.nsurf2));
        let grosterme = (3.0 * dprim * dprim * carre - dsecn) * cube;

        let p1 = self.dnplan.dot(&self.nsurf2);
        let p2 = self.d2nplan.dot(&self.nsurf2);

        let temp =
            GpVec::zero().set_linear_form_2(grosterme * terms.ndotns2, &self.nplan, -grosterme, &self.nsurf2);
        self.d2n2w = GpVec::zero().set_linear_form_3(
            smallterm * p1 + terms.invnorm2 * p2,
            &self.nplan,
            smallterm * terms.ndotns2 + 2.0 * terms.invnorm2 * p1,
            &self.dnplan,
            terms.ndotns2 * terms.invnorm2,
            &self.d2nplan,
        );
        self.d2n2w = self.d2n2w.added(&temp);

        let temp = GpVec::from_xyz(
            &self
                .pts1
                .coord
                .added(&self.pts2.coord)
                .multiply(0.5)
                .subtracted(&self.ptgui.coord),
        );
        self.d2edt2.set_value(
            1,
            self.d2nplan.dot(&temp) - 2.0 * self.dnplan.dot(&self.d1gui) - self.nplan.dot(&self.d2gui),
        );
        self.d2edt2
            .set_value(2, self.ray1 * self.d2n1w.x() - self.ray2 * self.d2n2w.x());
        self.d2edt2
            .set_value(3, self.ray1 * self.d2n1w.y() - self.ray2 * self.d2n2w.y());
        self.d2edt2
            .set_value(4, self.ray1 * self.d2n1w.z() - self.ray2 * self.d2n2w.z());
    }
}
