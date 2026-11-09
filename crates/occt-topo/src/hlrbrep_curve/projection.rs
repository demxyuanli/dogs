//! `HLRBRep_Curve` 2D evaluation and analytic reconstruction:
//! `Value3D` / `D0` / `D1` (`lxx:44-62`), `Value` / `D0` / `D1` / `D2` / `D3`
//! / `DN` (`cxx:315-412`, `lxx:121-127`), `Tangent` (`cxx:301-312`), `Line`
//! (`cxx:415-423`), `Circle` (`cxx:426-431`), `Ellipse` (`cxx:434-459`),
//! `Hyperbola` (`cxx:462-465`) and `Parabola` (`cxx:468-471`).

use occt_core::gp::{
    GpAx22d, GpCirc2d, GpDir2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpPnt, GpPnt2d, GpVec,
    GpVec2d,
};

use super::HlrBRepCurve;

/// `RealLast()` / `RealFirst()` (`Standard_Real.hxx`), as used by the
/// `LProp_CurveUtils` chord fallback.
const REAL_LAST: f64 = f64::MAX;
const REAL_FIRST: f64 = -f64::MAX;

impl HlrBRepCurve {
    /// `D0(U, gp_Pnt&)` (`lxx:55-58`): the 3D curve point.
    pub fn d0_3d(&self, u: f64) -> GpPnt {
        self.value3d(u)
    }

    /// `D1(U, gp_Pnt&, gp_Vec&)` (`lxx:62-65`).
    pub fn d1_3d(&self, u: f64) -> (GpPnt, GpVec) {
        match &self.curve {
            Some(c) => c.d1(u),
            None => (GpPnt::new(0.0, 0.0, 0.0), GpVec::new(0.0, 0.0, 0.0)),
        }
    }

    /// `Value(U)` (`lxx:121-127`).
    pub fn value(&self, u: f64) -> GpPnt2d {
        self.d0(u)
    }

    /// `myProj->Transform(P)` of `Poles(C, TP)` / `PolesAndWeights(C, TP, TW)`
    /// (`HLRBRep_Curve.cxx:524-527`, `:596-599`): the projector's plain
    /// transformation, then the `X` / `Y` pair. No projector leaves the point
    /// untouched, which is the `X` / `Y` of a point already in the view plane.
    pub fn project_to_2d(&self, p: &GpPnt) -> GpPnt2d {
        let mut q = *p;
        if let Some(proj) = &self.proj {
            q.transform(proj.transformation());
        }
        GpPnt2d::new(q.x(), q.y())
    }

    /// `D0(U, gp_Pnt2d&)` (`cxx:315-331`): transform the 3D point and project.
    pub fn d0(&self, u: f64) -> GpPnt2d {
        let p3d = self.value3d(u);
        match &self.proj {
            Some(proj) => {
                let mut pout = GpPnt2d::new(0.0, 0.0);
                proj.project(&p3d, &mut pout);
                pout
            }
            None => GpPnt2d::new(p3d.x(), p3d.y()),
        }
    }

    /// `D1(U, gp_Pnt2d&, gp_Vec2d&)` (`cxx:333-362`).
    ///
    /// The perspective arm uses the closed formula printed at `cxx:337-343`;
    /// the axonometric arm goes through `HLRAlgo_Projector::Project` with the
    /// derivative (`cxx:360-361`, the `OCC155` fix).
    pub fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d) {
        let (p3d, v13d) = self.d1_3d(u);
        if self.perspective() {
            let trsf = self.transformation();
            let mut p = p3d;
            let mut v = v13d;
            p.transform(&trsf);
            v.transform(&trsf);
            let f = self.focus().unwrap_or(0.0);
            let r = 1.0 - p.z() / f;
            let e = v.z() / (f * r * r);
            (
                GpPnt2d::new(p.x() / r, p.y() / r),
                GpVec2d::new(v.x() / r + p.x() * e, v.y() / r + p.y() * e),
            )
        } else {
            match &self.proj {
                Some(proj) => proj.project_with_derivative(&p3d, &v13d),
                None => (
                    GpPnt2d::new(p3d.x(), p3d.y()),
                    GpVec2d::new(v13d.x(), v13d.y()),
                ),
            }
        }
    }

    /// `D2(U, gp_Pnt2d&, gp_Vec2d&, gp_Vec2d&)` (`cxx:366-402`).
    pub fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        let (p3d, v13d, v23d) = match &self.curve {
            Some(c) => c.d2(u),
            None => (
                GpPnt::new(0.0, 0.0, 0.0),
                GpVec::new(0.0, 0.0, 0.0),
                GpVec::new(0.0, 0.0, 0.0),
            ),
        };
        let trsf = self.transformation();
        let mut p = p3d;
        let mut v1 = v13d;
        let mut v2 = v23d;
        p.transform(&trsf);
        v1.transform(&trsf);
        v2.transform(&trsf);
        if self.perspective() {
            // `cxx:386-395`
            let f = self.focus().unwrap_or(0.0);
            let r = 1.0 - p.z() / f;
            let q = f * r * r;
            let e = v1.z() / q;
            let c = e * v1.z() / (f * r);
            (
                GpPnt2d::new(p.x() / r, p.y() / r),
                GpVec2d::new(v1.x() / r + p.x() * e, v1.y() / r + p.y() * e),
                GpVec2d::new(
                    v2.x() / r + 2.0 * v1.x() * e + p.x() * v2.z() / q + 2.0 * p.x() * c,
                    v2.y() / r + 2.0 * v1.y() * e + p.y() * v2.z() / q + 2.0 * p.y() * c,
                ),
            )
        } else {
            // `cxx:397-402`
            (
                GpPnt2d::new(p.x(), p.y()),
                GpVec2d::new(v1.x(), v1.y()),
                GpVec2d::new(v2.x(), v2.y()),
            )
        }
    }

    /// `D3(...)` (`cxx:406`): an empty body in OCCT, so the out-parameters are
    /// left untouched. The port answers zeroes; no caller reads them.
    pub fn d3(&self, _u: f64) -> (GpPnt2d, GpVec2d, GpVec2d, GpVec2d) {
        (
            GpPnt2d::new(0.0, 0.0),
            GpVec2d::new(0.0, 0.0),
            GpVec2d::new(0.0, 0.0),
            GpVec2d::new(0.0, 0.0),
        )
    }

    /// `DN(U, N)` (`cxx:410-412`): `return gp_Vec2d();`.
    pub fn dn(&self, _u: f64, _n: i32) -> GpVec2d {
        GpVec2d::new(0.0, 0.0)
    }

    /// `Tangent(AtStart, P, D)` (`cxx:301-312`).
    ///
    /// OCCT builds `HLRBRep_CLProps CLP(2, Epsilon(1.))` on the 2D projection
    /// and asks `IsTangentDefined` then `Tangent`. With derivative order 2 the
    /// significance test walks order 1 then order 2
    /// (`LProp_CurveUtils.hxx:333-381`), and the order-2 arm orients the
    /// second derivative by the chord `P(u-delta)...P(u+delta)`
    /// (`LProp_CurveUtils.hxx:179-217`). Both are reproduced here; the
    /// parameter passed to the properties is the *3D* `U`, exactly as
    /// `cxx:305-310` does.
    pub fn tangent(&self, at_start: bool) -> Option<(GpPnt2d, GpDir2d)> {
        let u = if at_start {
            self.first_parameter
        } else {
            self.last_parameter
        };
        let p = self.d0(u);
        let eps = occt_core::precision::REAL_EPSILON;

        let (_, v1) = self.d1(u);
        if v1.x() * v1.x() + v1.y() * v1.y() > eps * eps {
            return Some((p, GpDir2d::from_vec2d(&v1).ok()?));
        }
        let (_, _, mut av) = self.d2(u);
        if av.x() * av.x() + av.y() * av.y() <= eps * eps {
            return None;
        }
        let u_inf = self.first_parameter();
        let u_sup = self.last_parameter();
        let du = if u_sup >= REAL_LAST || u_inf <= REAL_FIRST {
            0.0
        } else {
            u_sup - u_inf
        };
        let delta = (du * 1.0e-3).max(1.0e-7);
        let other_u = if u - u_inf < delta { u + delta } else { u - delta };
        let p1 = self.d0(u.min(other_u));
        let p2 = self.d0(u.max(other_u));
        let chord = GpVec2d::new(p2.x() - p1.x(), p2.y() - p1.y());
        if av.dot(&chord) < 0.0 {
            av.multiply_scalar(-1.0);
        }
        Some((p, GpDir2d::from_vec2d(&av).ok()?))
    }

    /// `Line()` (`cxx:415-423`): `gp_Lin2d(P, V)` from `D1(0)`.
    pub fn line(&self) -> Option<GpLin2d> {
        let (p, v) = self.d1(0.0);
        Some(GpLin2d::from_pnt_dir(p, GpDir2d::from_vec2d(&v).ok()?))
    }

    /// `Circle()` (`cxx:426-431`): transform the 3D circle into the view frame
    /// then `ProjLib::Project(gp_Pln(gp::XOY()), C)`.
    ///
    /// `ProjLib_Plane::Project(gp_Circ)` (`ProjLib_Plane.cxx:110-124`) maps the
    /// centre and both axis directions through `EvalPnt2d` / `EvalDir2d`
    /// (`ProjLib_Plane.cxx:87-96`); for `gp::XOY()` those are the identity on
    /// `X` / `Y`, which is inlined here instead of porting the whole
    /// `ProjLib_Plane` class.
    pub fn circle(&self) -> Option<GpCirc2d> {
        let mut c = self.curve.as_ref()?.gp_circ()?;
        c.transform(&self.transformation());
        let pos = c.position();
        let p2d = GpPnt2d::new(c.location().x(), c.location().y());
        let x2d = GpDir2d::new(pos.x_direction().x(), pos.x_direction().y()).ok()?;
        let y2d = GpDir2d::new(pos.y_direction().x(), pos.y_direction().y()).ok()?;
        let ax = GpAx22d::new(p2d, x2d, y2d).ok()?;
        Some(GpCirc2d::new(ax, c.radius()))
    }

    /// `Ellipse()` (`cxx:434-459`): a projected 3D ellipse, or the flattened
    /// form of a projected circle.
    pub fn ellipse(&self) -> Option<GpElips2d> {
        let curve = self.curve.as_ref()?;
        let trsf = self.transformation();
        if self.curve_kind_3d() == crate::geom_abs_curve_type::GeomAbsCurveType::Ellipse {
            // `cxx:438-443`
            let mut e = curve.gp_ellipse()?;
            e.transform(&trsf);
            let pos = e.position();
            let p2d = GpPnt2d::new(e.location().x(), e.location().y());
            let x2d = GpDir2d::new(pos.x_direction().x(), pos.x_direction().y()).ok()?;
            let y2d = GpDir2d::new(pos.y_direction().x(), pos.y_direction().y()).ok()?;
            let ax = GpAx22d::new(p2d, x2d, y2d).ok()?;
            return Some(GpElips2d::new(ax, e.major_radius(), e.minor_radius()));
        }
        // `cxx:444-459`: this is a circle.
        let mut c = curve.gp_circ()?;
        c.transform(&trsf);
        let circ_axis = c.axis();
        let d1 = *circ_axis.direction();
        let dz = occt_core::gp::GpDir::new(0.0, 0.0, 1.0).ok()?;
        let d3 = d1.crossed(&dz).ok()?;
        let d2 = d1.crossed(&d3).ok()?;
        let rap = (d2.x() * d2.x() + d2.y() * d2.y()).sqrt();
        let d = GpDir2d::new(d1.y(), -d1.x()).ok()?;
        let p = GpPnt2d::new(c.location().x(), c.location().y());
        let ax = GpAx22d::from_xdir(p, d);
        let mut el = GpElips2d::new(ax, c.radius(), c.radius() * rap);
        if d1.z() < 0.0 {
            // `gp_Elips2d::Reverse()` reverses both axis directions.
            el.pos.x_reverse();
            el.pos.y_reverse();
        }
        Some(el)
    }

    /// `Hyperbola()` (`cxx:462-465`): `return gp_Hypr2d();`.
    pub fn hyperbola(&self) -> GpHypr2d {
        GpHypr2d::default()
    }

    /// `Parabola()` (`cxx:468-471`): `return gp_Parab2d();`.
    pub fn parabola(&self) -> GpParab2d {
        GpParab2d::default()
    }
}
