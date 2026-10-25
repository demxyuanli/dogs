//! `GeomFill_DraftTrihedron` — the trihedron of a deported (draft-angle)
//! sweep: a constant binormal plus a normal tilted by `Angle` from the
//! surface's `v` direction.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_DraftTrihedron.hxx` and
//! `GeomFill_DraftTrihedron.cxx:32-371` (OCCT 8.0.0): `DDeriv` `cxx:32-45`,
//! constructor `cxx:47-52`, `SetAngle` `cxx:56-60`, `D0` `cxx:62-100`,
//! `D1` `cxx:106-157`, `D2` `cxx:161-236`, `Copy` `cxx:238-243`,
//! `NbIntervals` `cxx:245-268`, `Intervals` `cxx:272-293`,
//! `GetAverageLaw` `cxx:297-323`, `IsConstant` `cxx:327-330`,
//! `IsOnlyBy3dCurve` `cxx:334-371`.
//!
//! Translation notes:
//! - `B` (`hxx`) is a `gp_Vec` normalized by the constructor; it is kept a
//!   unit [`GpVec`].
//! - `myAngle` / `myCos` are `M_PI/2 + Angle` and its cosine (`cxx:58-60`);
//!   `Copy` rebuilds the law from `myAngle - M_PI/2` (`cxx:241`).
//! - `myCurve->GetType()` is [`crate::extrema_pc::extrema_curve_type`], the
//!   same `GeomAdaptor_Curve::GetType()` the rest of the port uses;
//!   `myCurve->Line().Direction()` / `Circle().Axis()` etc. are the
//!   `gp_line` / `gp_circ` / `gp_ellipse` / `gp_hyperbola` / `gp_parabola`
//!   accessors of [`crate::curve::Curve`].
//! - `Standard_OutOfRange` from the two switches is a panic, OCCT's
//!   exception counterpart (same convention as `GeomFill_Frenet`).

use std::sync::Arc;

use occt_core::gp::GpVec;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::ANGULAR;

use crate::curve::Curve;
use crate::extrema_pc::{extrema_curve_type, GeomAbsCurveType};
use crate::geom_fill::frenet::deriv::dderiv;
use crate::geom_fill::trihedron_law::{TrihedronD1, TrihedronD2, TrihedronLaw, TrihedronLawCurve};

/// `GeomFill_DraftTrihedron` (`GeomFill_DraftTrihedron.hxx:22-65`).
pub struct DraftTrihedron {
    /// `myCurve` / `myTrimmed` (`GeomFill_TrihedronLaw.hxx:111-112`).
    data: TrihedronLawCurve,
    /// `B` (`hxx`): the constant bi-normal, normalized by the constructor.
    b: GpVec,
    /// `myAngle` (`hxx`): `M_PI/2 + Angle`.
    angle: f64,
    /// `myCos` (`hxx`): `cos(myAngle)`.
    cos: f64,
}

impl DraftTrihedron {
    /// `GeomFill_DraftTrihedron(BiNormal, Angle)` (`cxx:47-52`): `B` is
    /// normalized and `SetAngle` applied.
    pub fn new(bi_normal: &GpVec, angle: f64) -> Self {
        let mut me = Self {
            data: TrihedronLawCurve::new(),
            b: bi_normal.normalized(),
            angle: 0.0,
            cos: 0.0,
        };
        me.set_angle(angle);
        me
    }

    /// `GeomFill_DraftTrihedron::SetAngle` (`cxx:56-60`).
    pub fn set_angle(&mut self, angle: f64) {
        self.angle = std::f64::consts::FRAC_PI_2 + angle;
        self.cos = self.angle.cos();
    }

    /// `B` (`hxx`).
    pub fn bi_normal(&self) -> GpVec {
        self.b
    }

    /// `myAngle` (`hxx`); `Copy` feeds `myAngle - M_PI/2` back in.
    pub fn angle(&self) -> f64 {
        self.angle
    }

    /// `myTrimmed` (`hxx`, inherited from the base).
    fn trimmed(&self) -> Result<Arc<dyn Curve>, String> {
        self.data.curve_or_err("GeomFill_DraftTrihedron")
    }

    /// The `NbIntervals` / `Intervals` switch (`cxx:250-267`, `cxx:274-291`).
    fn raise_shape(s: u8) -> Shape {
        match s {
            0 => Shape::C2,
            2 => Shape::C3,
            4 | 5 | 6 => Shape::CN,
            _ => panic!("Standard_OutOfRange: GeomFill_DraftTrihedron::Intervals"),
        }
    }
}

impl TrihedronLaw for DraftTrihedron {
    /// `GeomFill_TrihedronLaw::SetCurve` (`GeomFill_TrihedronLaw.cxx:25-30`).
    fn set_curve(&self, curve: Arc<dyn Curve>) -> bool {
        self.data.set_curve(curve)
    }

    /// `GeomFill_DraftTrihedron::Copy` (`cxx:238-243`): rebuilt from
    /// `(B, myAngle - M_PI/2)` and the same curve.
    fn copy(&self) -> Arc<dyn TrihedronLaw> {
        let opy = DraftTrihedron::new(&self.b, self.angle - std::f64::consts::FRAC_PI_2);
        if let Some(curve) = self.data.curve() {
            TrihedronLaw::set_curve(&opy, curve);
        }
        Arc::new(opy)
    }

    /// `GeomFill_DraftTrihedron::D0` (`cxx:62-100`).
    fn d0(&self, param: f64) -> Result<(GpVec, GpVec, GpVec), String> {
        let curve = self.trimmed()?;
        let (_p, mut tangent) = curve.d1(param);
        tangent = tangent.normalized();

        let mut b = tangent.crossed(&self.b);
        let normb = b.magnitude();
        // `cxx:78-81`: OCCT divides first, then reports the failure.
        b = b.divided(normb);
        if normb < 1.0e-12 {
            return Err("GeomFill_DraftTrihedron::D0: degenerate bi-normal".into());
        }

        let v = b.crossed(&tangent);
        let mu = self.cos;
        let normal = b
            .multiplied_scalar((1.0 - mu * mu).sqrt())
            .added(&v.multiplied_scalar(mu));

        let tangent = normal.crossed(&self.b).normalized();
        let bi_normal = tangent.crossed(&normal);
        Ok((tangent, normal, bi_normal))
    }

    /// `GeomFill_DraftTrihedron::D1` (`cxx:106-157`).
    fn d1(&self, param: f64) -> Result<TrihedronD1, String> {
        let curve = self.trimmed()?;
        let (_p, mut tangent, aux) = curve.d2(param);

        let norm_t = tangent.magnitude();
        tangent = tangent.divided(norm_t);
        let mut d_tangent = tangent
            .multiplied_scalar(-tangent.dot(&aux))
            .added(&aux);
        d_tangent = d_tangent.divided(norm_t);

        let mut b = tangent.crossed(&self.b);
        let normb = b.magnitude();
        if normb < 1.0e-12 {
            return Err("GeomFill_DraftTrihedron::D1: degenerate bi-normal".into());
        }
        b = b.divided(normb);
        let aux = d_tangent.crossed(&self.b);
        let mut db = b.multiplied_scalar(-b.dot(&aux)).added(&aux);
        db = db.divided(normb);

        let v = b.crossed(&tangent);
        let dv = db.crossed(&tangent).added(&b.crossed(&d_tangent));

        let mu = self.cos;
        let rac = (1.0 - mu * mu).sqrt();
        let normal = b.multiplied_scalar(rac).added(&v.multiplied_scalar(mu));
        let d_normal = db.multiplied_scalar(rac).added(&dv.multiplied_scalar(mu));

        let tangent = normal.crossed(&self.b);
        let norm_t = tangent.magnitude();
        let tangent = tangent.divided(norm_t);
        let aux = d_normal.crossed(&self.b);
        let d_tangent = tangent
            .multiplied_scalar(-tangent.dot(&aux))
            .added(&aux)
            .divided(norm_t);

        let bi_normal = tangent.crossed(&normal);
        let d_bi_normal = d_tangent
            .crossed(&normal)
            .added(&tangent.crossed(&d_normal));

        Ok(TrihedronD1 {
            tangent,
            dtangent: d_tangent,
            normal,
            dnormal: d_normal,
            binormal: bi_normal,
            dbinormal: d_bi_normal,
        })
    }

    /// `GeomFill_DraftTrihedron::D2` (`cxx:161-236`).
    fn d2(&self, param: f64) -> Result<TrihedronD2, String> {
        let curve = self.trimmed()?;
        let (_p, t, aux, aux2) = curve.d3(param);
        let d2t = dderiv(&t, &aux, &aux2);

        let norm_t = t.magnitude();
        let t = t.divided(norm_t);
        let dot = t.dot(&aux);
        let d_tangent = t.multiplied_scalar(-dot).added(&aux).divided(norm_t);

        let b = t.crossed(&self.b);
        let normb = b.magnitude();
        if normb < 1.0e-12 {
            return Err("GeomFill_DraftTrihedron::D2: degenerate bi-normal".into());
        }

        let aux = d_tangent.crossed(&self.b);
        let aux2 = d2t.crossed(&self.b);
        let d2b = dderiv(&b, &aux, &aux2);
        let b = b.divided(normb);
        let dot = b.dot(&aux);
        let db = b.multiplied_scalar(-dot).added(&aux).divided(normb);

        let v = b.crossed(&t);
        let dv = db.crossed(&t).added(&b.crossed(&d_tangent));
        let d2v = d2b
            .crossed(&t)
            .added(&db.crossed(&d_tangent).multiplied_scalar(2.0))
            .added(&b.crossed(&d2t));

        let mu = self.cos;
        let rac = (1.0 - mu * mu).sqrt();
        let normal = b.multiplied_scalar(rac).added(&v.multiplied_scalar(mu));
        let d_normal = db.multiplied_scalar(rac).added(&dv.multiplied_scalar(mu));
        let d2_normal = d2b.multiplied_scalar(rac).added(&d2v.multiplied_scalar(mu));

        let tangent = normal.crossed(&self.b);
        let norm_t = tangent.magnitude();

        let aux = d_normal.crossed(&self.b);
        let aux2 = d2_normal.crossed(&self.b);
        let d2_tangent = dderiv(&tangent, &aux, &aux2);
        let tangent = tangent.divided(norm_t);
        let dot = tangent.dot(&aux);
        let d_tangent = tangent
            .multiplied_scalar(-dot)
            .added(&aux)
            .divided(norm_t);

        let bi_normal = tangent.crossed(&normal);
        let d_bi_normal = d_tangent
            .crossed(&normal)
            .added(&tangent.crossed(&d_normal));
        let d2_bi_normal = d2_tangent
            .crossed(&normal)
            .added(&d_tangent.crossed(&d_normal).multiplied_scalar(2.0))
            .added(&tangent.crossed(&d2_normal));

        Ok(TrihedronD2 {
            tangent,
            dtangent: d_tangent,
            d2tangent: d2_tangent,
            normal,
            dnormal: d_normal,
            d2normal: d2_normal,
            binormal: bi_normal,
            dbinormal: d_bi_normal,
            d2binormal: d2_bi_normal,
        })
    }

    /// `GeomFill_DraftTrihedron::NbIntervals` (`cxx:245-268`).
    fn nb_intervals(&self, s: u8) -> i32 {
        let tmp_s = Self::raise_shape(s);
        match self.data.curve() {
            Some(c) => c.nb_intervals(tmp_s as u8),
            None => 1,
        }
    }

    /// `GeomFill_DraftTrihedron::Intervals` (`cxx:272-293`).
    fn intervals(&self, s: u8) -> Vec<f64> {
        let tmp_s = Self::raise_shape(s);
        match self.data.curve() {
            Some(c) => c.parameter_intervals(tmp_s as u8),
            None => vec![0.0, 0.0],
        }
    }

    /// `GeomFill_TrihedronLaw::SetInterval` (`GeomFill_TrihedronLaw.cxx:57-60`).
    fn set_interval(&self, first: f64, last: f64) -> Result<(), String> {
        self.data.set_interval(first, last)
    }

    /// `GeomFill_TrihedronLaw::GetInterval` (`GeomFill_TrihedronLaw.cxx:62-67`).
    fn get_interval(&self) -> (f64, f64) {
        self.data.interval()
    }

    /// `GeomFill_DraftTrihedron::GetAverageLaw` (`cxx:297-323`).
    fn get_average_law(&self) -> (GpVec, GpVec, GpVec) {
        const NUM: usize = 20;
        let mut a_tangent = GpVec::zero();
        let mut a_normal = GpVec::zero();
        let mut a_bi_normal = GpVec::zero();
        let (first, last) = self.data.interval();
        let step = (last - first) / NUM as f64;
        for i in 0..=NUM {
            let mut param = first + i as f64 * step;
            if param > last {
                param = last;
            }
            if let Ok((t, n, bn)) = self.d0(param) {
                a_tangent = a_tangent.added(&t);
                a_normal = a_normal.added(&n);
                a_bi_normal = a_bi_normal.added(&bn);
            }
        }
        // `cxx:320-322`: note the odd order -- `ANormal` and `ABiNormal` are
        // divided before `ATangent`, and no renormalization follows.
        a_normal = a_normal.divided((NUM + 1) as f64);
        a_bi_normal = a_bi_normal.divided((NUM + 1) as f64);
        a_tangent = a_tangent.divided((NUM + 1) as f64);
        (a_tangent, a_normal, a_bi_normal)
    }

    /// `GeomFill_DraftTrihedron::IsConstant` (`cxx:327-330`).
    fn is_constant(&self) -> bool {
        match self.data.curve() {
            Some(c) => matches!(extrema_curve_type(c.as_ref()), GeomAbsCurveType::Line),
            None => false,
        }
    }

    /// `GeomFill_DraftTrihedron::IsOnlyBy3dCurve` (`cxx:334-371`): the curve's
    /// plane normal must be parallel to `B`.
    fn is_only_by_3d_curve(&self) -> bool {
        let Some(curve) = self.data.curve() else {
            return false;
        };
        match extrema_curve_type(curve.as_ref()) {
            GeomAbsCurveType::Circle => curve
                .gp_circ()
                .map(|c| c.axis())
                .map(|axis| {
                    GpVec::from_xyz(axis.direction().xyz()).is_parallel_ang(&self.b, ANGULAR)
                })
                .unwrap_or(false),
            GeomAbsCurveType::Ellipse => curve
                .gp_ellipse()
                .map(|c| *c.axis())
                .map(|axis| {
                    GpVec::from_xyz(axis.direction().xyz()).is_parallel_ang(&self.b, ANGULAR)
                })
                .unwrap_or(false),
            GeomAbsCurveType::Hyperbola => curve
                .gp_hyperbola()
                .map(|c| *c.axis())
                .map(|axis| {
                    GpVec::from_xyz(axis.direction().xyz()).is_parallel_ang(&self.b, ANGULAR)
                })
                .unwrap_or(false),
            GeomAbsCurveType::Parabola => curve
                .gp_parabola()
                .map(|c| *c.axis())
                .map(|axis| {
                    GpVec::from_xyz(axis.direction().xyz()).is_parallel_ang(&self.b, ANGULAR)
                })
                .unwrap_or(false),
            // `cxx:358-361`: for a line the curve plane normal is the line
            // direction itself.
            GeomAbsCurveType::Line => curve
                .gp_line()
                .map(|l| {
                    GpVec::from_xyz(l.direction().xyz()).is_parallel_ang(&self.b, ANGULAR)
                })
                .unwrap_or(false),
            // `cxx:363-364`: "pas de risques".
            _ => false,
        }
    }
}
