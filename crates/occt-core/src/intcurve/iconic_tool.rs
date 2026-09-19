//! `IntCurve_IConicTool` (TKGeomAlgo, `IntCurve/IntCurve_IConicTool.hxx`,
//! `.cxx`).
//!
//! Can be plugged into `IntImpParGen` / `IntCurve_IntImpConicParConic` in place
//! of a tool. It carries the same conic data as `IntCurve_PConic` plus an
//! `Abs_To_Object` transform, and it is the `Abs_To_Object` that makes
//! `Distance` / `GradDistance` implicit-curve functions: the signed distance of
//! an *absolute* point is evaluated after mapping the point into the conic's own
//! frame, and the gradient is mapped back.
use crate::elib::clib2d;
use crate::gp::{GpAx22d, GpCirc2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpPnt2d, GpTrsf2d, GpVec2d};
use crate::kernel::geomabs::CurveType;

/// `IntCurve_IConicTool` (`IntCurve_IConicTool.hxx:36-83`). `prm1`/`prm2`/`prm3`
/// are the three scalars; which one means what depends on `type`, exactly as
/// the OCCT `#define` block `cxx:27-51` records:
///
/// | `type`      | `prm1`  | `prm2`       | `prm3` |
/// |-------------|---------|--------------|--------|
/// | `Line`      | `a`     | `b`          | `c`    |
/// | `Circle`    | `r`     | `x0`         | `y0`   |
/// | `Ellipse`   | `a`     | `b`          | `c`    |
/// | `Parabola`  | `f`     | `2p` = `4f`  | -      |
/// | `Hyperbola` | `a`     | `b`          | -      |
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntCurveIConicTool {
    prm1: f64,
    prm2: f64,
    prm3: f64,
    axis: GpAx22d,
    type_curve: CurveType,
    abs_to_object: GpTrsf2d,
}

impl Default for IntCurveIConicTool {
    /// `IntCurve_IConicTool()` (`cxx:54-60`): all scalars zero, `Axis` default
    /// constructed, `Abs_To_Object` default constructed, type
    /// `GeomAbs_OtherCurve`. The OCCT comment asks for the type to be left
    /// undefined; `GeomAbs_OtherCurve` is what the initialiser list uses.
    fn default() -> Self {
        Self {
            prm1: 0.0,
            prm2: 0.0,
            prm3: 0.0,
            axis: GpAx22d::default(),
            type_curve: CurveType::OtherCurve,
            abs_to_object: GpTrsf2d::default(),
        }
    }
}

impl IntCurveIConicTool {
    /// `IntCurve_IConicTool(const gp_Lin2d&)` (`cxx:73-82`). `prm1..prm3` are
    /// `Line.Coefficients(a, b, c)`; the axis is
    /// `gp_Ax22d(Line.Position(), true)`, i.e. right-handed.
    pub fn from_lin2d(line: &GpLin2d) -> Self {
        let (a, b, c) = lin2d_coefficients(line);
        Self {
            prm1: a,
            prm2: b,
            prm3: c,
            axis: GpAx22d::from_xdir(line.pos.loc, line.pos.vdir),
            type_curve: CurveType::Line,
            abs_to_object: GpTrsf2d::default(),
        }
    }

    /// `IntCurve_IConicTool(const gp_Elips2d&)` (`cxx:85-97`). `prm3` is the
    /// linear eccentricity `c = sqrt(a^2 - b^2)`; `Abs_To_Object` is the change
    /// of basis `gp::OX2d() -> Axis.XAxis()`.
    pub fn from_elips2d(elips: &GpElips2d) -> Self {
        let a = elips.major_radius;
        let b = elips.minor_radius;
        let mut abs_to_object = GpTrsf2d::identity();
        abs_to_object.set_transformation(&crate::gp::ox2d(), &elips.pos.x_axis());
        Self {
            prm1: a,
            prm2: b,
            prm3: (a * a - b * b).sqrt(),
            axis: elips.pos,
            type_curve: CurveType::Ellipse,
            abs_to_object,
        }
    }

    /// `IntCurve_IConicTool(const gp_Circ2d&)` (`cxx:100-111`): `prm1 = r`,
    /// `prm2 = x0`, `prm3 = y0` (the centre in *absolute* coordinates).
    pub fn from_circ2d(c: &GpCirc2d) -> Self {
        let mut abs_to_object = GpTrsf2d::identity();
        abs_to_object.set_transformation(&crate::gp::ox2d(), &c.pos.x_axis());
        Self {
            prm1: c.radius,
            prm2: c.pos.point.x(),
            prm3: c.pos.point.y(),
            axis: c.pos,
            type_curve: CurveType::Circle,
            abs_to_object,
        }
    }

    /// `IntCurve_IConicTool(const gp_Parab2d&)` (`cxx:114-124`):
    /// `prm1 = Focal`, `prm2 = 4 * Focal`.
    pub fn from_parab2d(p: &GpParab2d) -> Self {
        let mut abs_to_object = GpTrsf2d::identity();
        abs_to_object.set_transformation(&crate::gp::ox2d(), &p.pos.x_axis());
        Self {
            prm1: p.focal,
            prm2: 4.0 * p.focal,
            prm3: 0.0,
            axis: p.pos,
            type_curve: CurveType::Parabola,
            abs_to_object,
        }
    }

    /// `IntCurve_IConicTool(const gp_Hypr2d&)` (`cxx:128-137`).
    pub fn from_hypr2d(h: &GpHypr2d) -> Self {
        let mut abs_to_object = GpTrsf2d::identity();
        abs_to_object.set_transformation(&crate::gp::ox2d(), &h.pos.x_axis());
        Self {
            prm1: h.major_radius,
            prm2: h.minor_radius,
            prm3: 0.0,
            axis: h.pos,
            type_curve: CurveType::Hyperbola,
            abs_to_object,
        }
    }

    /// `IntCurve_IConicTool::Value` (`cxx:140-158`). OCCT prints
    /// `"### Erreur sur le  type de la courbe ###"` and returns `(0, 0)` for an
    /// unset type; the Rust port returns the same value without the print.
    pub fn value(&self, x: f64) -> GpPnt2d {
        match self.type_curve {
            CurveType::Line => clib2d::line_value_ax2d(x, &self.axis.x_axis()),
            CurveType::Ellipse => {
                clib2d::ellipse_value_ax22d(x, &self.axis, self.prm1, self.prm2)
            }
            CurveType::Circle => clib2d::circle_value_ax22d(x, &self.axis, self.prm1),
            CurveType::Parabola => clib2d::parabola_value_ax22d(x, &self.axis, self.prm1),
            CurveType::Hyperbola => {
                clib2d::hyperbola_value_ax22d(x, &self.axis, self.prm1, self.prm2)
            }
            _ => GpPnt2d::new(0.0, 0.0),
        }
    }

    /// `IntCurve_IConicTool::D1` (`cxx:162-186`).
    pub fn d1(&self, x: f64) -> (GpPnt2d, GpVec2d) {
        match self.type_curve {
            CurveType::Line => clib2d::line_d1_ax2d(x, &self.axis.x_axis()),
            CurveType::Ellipse => {
                clib2d::ellipse_d1_ax22d(x, &self.axis, self.prm1, self.prm2)
            }
            CurveType::Circle => clib2d::circle_d1_ax22d(x, &self.axis, self.prm1),
            CurveType::Parabola => clib2d::parabola_d1_ax22d(x, &self.axis, self.prm1),
            CurveType::Hyperbola => {
                clib2d::hyperbola_d1_ax22d(x, &self.axis, self.prm1, self.prm2)
            }
            // `cxx:181-184`: no output is written for an unknown type. Returning
            // the uninitialised values OCCT leaves behind is not reproducible,
            // so zeros stand in.
            _ => (GpPnt2d::new(0.0, 0.0), GpVec2d::new(0.0, 0.0)),
        }
    }

    /// `IntCurve_IConicTool::D2` (`cxx:189-214`). The line arm normalises the
    /// second derivative to zero, per `cxx:194`.
    pub fn d2(&self, x: f64) -> (GpPnt2d, GpVec2d, GpVec2d) {
        match self.type_curve {
            CurveType::Line => {
                let (p, tan) = clib2d::line_d1_ax2d(x, &self.axis.x_axis());
                (p, tan, GpVec2d::new(0.0, 0.0))
            }
            CurveType::Ellipse => {
                clib2d::ellipse_d2_ax22d(x, &self.axis, self.prm1, self.prm2)
            }
            CurveType::Circle => clib2d::circle_d2_ax22d(x, &self.axis, self.prm1),
            CurveType::Parabola => clib2d::parabola_d2_ax22d(x, &self.axis, self.prm1),
            CurveType::Hyperbola => {
                clib2d::hyperbola_d2_ax22d(x, &self.axis, self.prm1, self.prm2)
            }
            _ => (
                GpPnt2d::new(0.0, 0.0),
                GpVec2d::new(0.0, 0.0),
                GpVec2d::new(0.0, 0.0),
            ),
        }
    }

    /// `IntCurve_IConicTool::Distance` (`cxx:220-279`) — the signed distance
    /// between the point and the *implicit* curve. The `AN_ELIPS` alternative
    /// body (`cxx:230-236`) is disabled upstream and therefore not ported.
    pub fn distance(&self, the_point: &GpPnt2d) -> f64 {
        match self.type_curve {
            CurveType::Line => self.prm1 * the_point.x() + self.prm2 * the_point.y() + self.prm3,
            CurveType::Ellipse => {
                let p = the_point.transformed(&self.abs_to_object);
                let x = p.x();
                let y = p.y() * (self.prm1 / self.prm2);
                (x * x + y * y).sqrt() - self.prm1
            }
            CurveType::Circle => {
                let dx = self.prm2 - the_point.x();
                let dy = self.prm3 - the_point.y();
                (dx * dx + dy * dy).sqrt() - self.prm1
            }
            CurveType::Parabola => {
                // `Distance(X, Y) = Y**2 - 2 P X` (`cxx:253`).
                let p = the_point.transformed(&self.abs_to_object);
                p.y() * p.y() - self.prm2 * p.x()
            }
            CurveType::Hyperbola => {
                // `Distance(X, Y) = (X/a)**2 - (Y/b)**2 - 1` for `X > 0`
                // (`cxx:257-272`). The negative branch negates `X/a` so the
                // gradient of `x -> |x|` keeps pushing back to the `+X` branch.
                let p = the_point.transformed(&self.abs_to_object);
                let aa = self.prm1 * self.prm1;
                let bb = self.prm2 * self.prm2;
                if p.x() > 0.0 {
                    (p.x() * p.x()) / aa - (p.y() * p.y()) / bb - 1.0
                } else {
                    (-p.x() * p.x()) / aa - (p.y() * p.y()) / bb - 1.0
                }
            }
            _ => 0.0,
        }
    }

    /// `IntCurve_IConicTool::GradDistance` (`cxx:282-370`). The `AN_ELIPS`
    /// alternative body (`cxx:307-330`) is disabled upstream and not ported.
    /// `Object_To_Abs` is `Abs_To_Object.Inverted()` (`cxx:26`); a singular
    /// `Abs_To_Object` cannot occur because `SetTransformation` is always given
    /// an orthonormal basis pair, so the inversion cannot fail.
    pub fn grad_distance(&self, the_point: &GpPnt2d) -> GpVec2d {
        let object_to_abs = self
            .abs_to_object
            .inverted()
            .expect("IntCurve_IConicTool::GradDistance: Abs_To_Object is never singular");
        match self.type_curve {
            CurveType::Line => GpVec2d::new(self.prm1, self.prm2),
            CurveType::Circle => {
                let p = the_point.transformed(&self.abs_to_object);
                let mut gradx = 0.0;
                let mut grady = 0.0;
                let temp1 = (p.y() * p.y() + p.x() * p.x()).sqrt();
                if temp1 != 0.0 {
                    gradx = p.x() / temp1;
                    grady = p.y() / temp1;
                }
                GpVec2d::new(gradx, grady).transformed(&object_to_abs)
            }
            CurveType::Ellipse => {
                let p = the_point.transformed(&self.abs_to_object);
                let mut gradx = 0.0;
                let mut grady = 0.0;
                let x = p.x();
                let y = p.y() * (self.prm1 / self.prm2);
                let temp1 = (y * y + x * x).sqrt();
                if temp1 != 0.0 {
                    gradx = x / temp1;
                    grady = (y * (self.prm1 / self.prm2)) / temp1;
                }
                GpVec2d::new(gradx, grady).transformed(&object_to_abs)
            }
            CurveType::Parabola => {
                // `Distance(X, Y) = Y**2 - 2 P X` (`cxx:350`).
                let p = the_point.transformed(&self.abs_to_object);
                GpVec2d::new(-self.prm2, p.y() + p.y()).transformed(&object_to_abs)
            }
            CurveType::Hyperbola => {
                // `Distance(X, Y) = (X/a)**2 - (Y/b)**2 - 1` (`cxx:357-363`).
                let p = the_point.transformed(&self.abs_to_object);
                let aa = self.prm1 * self.prm1;
                let bb = self.prm2 * self.prm2;
                GpVec2d::new(2.0 * p.x().abs() / aa, -2.0 * p.y() / bb).transformed(&object_to_abs)
            }
            _ => GpVec2d::new(0.0, 0.0),
        }
    }

    /// `IntCurve_IConicTool::FindParameter` (`cxx:372-414`). Circle and ellipse
    /// parameters are wrapped into `[0, 2*PI]` before being returned.
    pub fn find_parameter(&self, p: &GpPnt2d) -> f64 {
        let two_pi = 2.0 * std::f64::consts::PI;
        match self.type_curve {
            CurveType::Line => clib2d::line_parameter_ax2d(&self.axis.x_axis(), p),
            CurveType::Circle => {
                let mut param = clib2d::circle_parameter_ax22d(&self.axis, p);
                if param < 0.0 {
                    param += two_pi;
                }
                param
            }
            CurveType::Ellipse => {
                let mut param =
                    clib2d::ellipse_parameter_ax22d(&self.axis, self.prm1, self.prm2, p);
                if param < 0.0 {
                    param += two_pi;
                }
                param
            }
            CurveType::Parabola => clib2d::parabola_parameter_ax22d(&self.axis, p),
            CurveType::Hyperbola => {
                clib2d::hyperbola_parameter_ax22d(&self.axis, self.prm1, self.prm2, p)
            }
            _ => 0.0,
        }
    }
}

/// `gp_Lin2d::Coefficients(A, B, C)` (`gp_Lin2d.hxx:91-96`). The line is
/// `A*X + B*Y + C = 0` with `(A, B)` the Y direction of the line's own
/// placement (`A = Yd.Y`, `B = -Yd.X`); because that direction is already unit
/// the triple needs no normalisation. This is why the distance in
/// `IntCurve_IConicTool::Distance` is a *signed* distance.
fn lin2d_coefficients(l: &GpLin2d) -> (f64, f64, f64) {
    let dir = l.pos.vdir;
    let a = dir.y;
    let b = -dir.x;
    let c = -(a * l.pos.loc.x() + b * l.pos.loc.y());
    (a, b, c)
}
