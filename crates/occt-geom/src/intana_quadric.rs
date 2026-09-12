//! `IntAna_Quadric` — implicit quadric by 10 polynomial coefficients.
//! Source: `IntAna_Quadric.cxx` / `.hxx`.

use occt_core::elib::slib;
use occt_core::gp::{GpAx3, GpCone, GpCylinder, GpPln, GpPnt, GpSphere, GpTrsf};

/// Ten coefficients of
/// `CXX x^2 + CYY y^2 + CZZ z^2 + 2(CXY xy + CXZ xz + CYZ yz) + 2(CX x + CY y + CZ z) + CCte`.
#[derive(Debug, Clone)]
pub struct IntAnaQuadric {
    cxx: f64,
    cyy: f64,
    czz: f64,
    cxy: f64,
    cxz: f64,
    cyz: f64,
    cx: f64,
    cy: f64,
    cz: f64,
    ccte: f64,
    special_points: Vec<GpPnt>,
}

impl Default for IntAnaQuadric {
    fn default() -> Self {
        Self {
            cxx: 0.0,
            cyy: 0.0,
            czz: 0.0,
            cxy: 0.0,
            cxz: 0.0,
            cyz: 0.0,
            cx: 0.0,
            cy: 0.0,
            cz: 0.0,
            ccte: 1.0,
            special_points: Vec::new(),
        }
    }
}

impl IntAnaQuadric {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_plane(p: &GpPln) -> Self {
        let mut q = Self::new();
        q.set_plane(p);
        q
    }

    pub fn from_sphere(sph: &GpSphere) -> Self {
        let mut q = Self::new();
        q.set_sphere(sph);
        q
    }

    pub fn from_cylinder(cyl: &GpCylinder) -> Self {
        let mut q = Self::new();
        q.set_cylinder(cyl);
        q
    }

    pub fn from_cone(cone: &GpCone) -> Self {
        let mut q = Self::new();
        q.set_cone(cone);
        q
    }

    pub fn set_plane(&mut self, p: &GpPln) {
        let (cx, cy, cz, ccte) = p.coefficients();
        self.cx = 0.5 * cx;
        self.cy = 0.5 * cy;
        self.cz = 0.5 * cz;
        self.ccte = ccte;
        self.cxx = 0.0;
        self.cyy = 0.0;
        self.czz = 0.0;
        self.cxy = 0.0;
        self.cxz = 0.0;
        self.cyz = 0.0;
        self.special_points.clear();
    }

    pub fn set_cylinder(&mut self, cyl: &GpCylinder) {
        let c = cyl.coefficients();
        self.set_ten(c);
        self.special_points.clear();
    }

    pub fn set_cone(&mut self, cone: &GpCone) {
        let c = cone.coefficients();
        self.set_ten(c);
        self.special_points.clear();
        let a_v = -cone.radius() / cone.semi_angle().sin();
        self.special_points.push(slib::cone_value(cone, 0.0, a_v));
    }

    pub fn set_sphere(&mut self, sph: &GpSphere) {
        let c = sph.coefficients();
        self.set_ten(c);
        self.special_points.clear();
        self.special_points
            .push(slib::sphere_value(sph, 0.0, -std::f64::consts::FRAC_PI_2));
        self.special_points
            .push(slib::sphere_value(sph, 0.0, std::f64::consts::FRAC_PI_2));
    }

    fn set_ten(&mut self, c: [f64; 10]) {
        self.cxx = c[0];
        self.cyy = c[1];
        self.czz = c[2];
        self.cxy = c[3];
        self.cxz = c[4];
        self.cyz = c[5];
        self.cx = c[6];
        self.cy = c[7];
        self.cz = c[8];
        self.ccte = c[9];
    }

    pub fn coefficients(&self) -> [f64; 10] {
        [
            self.cxx, self.cyy, self.czz, self.cxy, self.cxz, self.cyz, self.cx, self.cy, self.cz,
            self.ccte,
        ]
    }

    pub fn special_points(&self) -> &[GpPnt] {
        &self.special_points
    }

    /// `IntAna_Quadric::NewCoefficients` in the local frame of `axis`.
    pub fn new_coefficients(&self, axis: &GpAx3) -> [f64; 10] {
        let mut trans = GpTrsf::identity();
        trans.set_transformation(axis);
        let trans = trans.inverted().unwrap_or(trans);
        let t11 = trans.value(1, 1);
        let t12 = trans.value(1, 2);
        let t13 = trans.value(1, 3);
        let t14 = trans.value(1, 4);
        let t21 = trans.value(2, 1);
        let t22 = trans.value(2, 2);
        let t23 = trans.value(2, 3);
        let t24 = trans.value(2, 4);
        let t31 = trans.value(3, 1);
        let t32 = trans.value(3, 2);
        let t33 = trans.value(3, 3);
        let t34 = trans.value(3, 4);
        let cxx = self.cxx;
        let cyy = self.cyy;
        let czz = self.czz;
        let cxy = self.cxy;
        let cxz = self.cxz;
        let cyz = self.cyz;
        let cx = self.cx;
        let cy = self.cy;
        let cz = self.cz;
        let ccte = self.ccte;
        let t11_p2 = t11 * t11;
        let t21_p2 = t21 * t21;
        let t31_p2 = t31 * t31;
        let t12_p2 = t12 * t12;
        let t22_p2 = t22 * t22;
        let t32_p2 = t32 * t32;
        let t13_p2 = t13 * t13;
        let t23_p2 = t23 * t23;
        let t33_p2 = t33 * t33;
        let t14_p2 = t14 * t14;
        let t24_p2 = t24 * t24;
        let t34_p2 = t34 * t34;
        let ccte_n = ccte
            + t14_p2 * cxx
            + t24_p2 * cyy
            + t34_p2 * czz
            + 2.0
                * (t14 * (cx + t24 * cxy + t34 * cxz)
                    + t24 * (cy + t34 * cyz)
                    + t34 * cz);
        let cxx_n = t11_p2 * cxx
            + t21_p2 * cyy
            + t31_p2 * czz
            + 2.0 * (t11 * (t21 * cxy + t31 * cxz) + t21 * t31 * cyz);
        let cyy_n = t12_p2 * cxx
            + t22_p2 * cyy
            + t32_p2 * czz
            + 2.0 * (t12 * (t22 * cxy + t32 * cxz) + t22 * t32 * cyz);
        let czz_n = t13_p2 * cxx
            + t33_p2 * czz
            + t23_p2 * cyy
            + 2.0 * (t13 * (t23 * cxy + t33 * cxz) + t23 * t33 * cyz);
        let cz_n = t13 * cx
            + t13 * (t14 * cxx + t24 * cxy + t34 * cxz)
            + t14 * (t23 * cxy + t33 * cxz)
            + t23 * (cy + t24 * cyy + t34 * cyz)
            + t33 * (t24 * cyz + cz + t34 * czz);
        let cx_n = t11 * (cx + t14 * cxx + t24 * cxy + t34 * cxz)
            + t14 * (t21 * cxy + t31 * cxz)
            + t21 * (cy + t24 * cyy + t34 * cyz)
            + t31 * (t24 * cyz + cz + t34 * czz);
        let cxy_n = t11 * (t12 * cxx + t22 * cxy + t32 * cxz)
            + t12 * (t21 * cxy + t31 * cxz)
            + t21 * (t22 * cyy + t32 * cyz)
            + t31 * (t22 * cyz + t32 * czz);
        let cxz_n = t11 * (t13 * cxx + t23 * cxy + t33 * cxz)
            + t13 * (t21 * cxy + t31 * cxz)
            + t21 * (t23 * cyy + t33 * cyz)
            + t31 * (t23 * cyz + t33 * czz);
        let cy_n = t12 * (cx + t14 * cxx + t24 * cxy + t34 * cxz)
            + t14 * (t22 * cxy + t32 * cxz)
            + t22 * (cy + t24 * cyy + t34 * cyz)
            + t32 * (cz + t24 * cyz + t34 * czz);
        let cyz_n = t12 * (t13 * cxx + t23 * cxy + t33 * cxz)
            + t13 * (t22 * cxy + t32 * cxz)
            + t22 * (t23 * cyy + t33 * cyz)
            + t32 * (t23 * cyz + t33 * czz);
        [cxx_n, cyy_n, czz_n, cxy_n, cxz_n, cyz_n, cx_n, cy_n, cz_n, ccte_n]
    }
}
