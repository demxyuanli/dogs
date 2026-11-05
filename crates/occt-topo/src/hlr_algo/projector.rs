//! Projector object used by hidden line removal. Source: `HLRAlgo_Projector`
//! (`ModelingAlgorithms/TKHLR/HLRAlgo/HLRAlgo_Projector.hxx`,
//! `HLRAlgo_Projector.lxx`, `HLRAlgo_Projector.cxx:38-396`).

use occt_core::gp::{GpAx2, GpAx3, GpDir, GpLin, GpPnt, GpPnt2d, GpTrsf, GpVec, GpVec2d};
use occt_core::precision::ANGULAR;

/// `TrsfType` (`cxx:113-150`): the four hard coded standard projections
/// (`0` axonometric, `1` top, `2` front, `3` axo) and `-1` for any other
/// transformation. The tolerances are the literal `1e-15` of the source.
fn trsf_type(trsf: &GpTrsf) -> i32 {
    const EPS: f64 = 1e-15;
    const C07: f64 = 0.7071067811865476;
    let m = trsf.vectorial_part();
    let v = |r: usize, c: usize| m.value(r, c);

    if (v(1, 1) - 1.0).abs() < EPS && (v(2, 2) - 1.0).abs() < EPS && (v(3, 3) - 1.0).abs() < EPS {
        return 1; // top
    } else if (v(1, 1) - C07).abs() < EPS
        && (v(1, 2) + 0.5).abs() < EPS
        && (v(1, 3) - 0.5).abs() < EPS
        && (v(2, 1) - C07).abs() < EPS
        && (v(2, 2) - 0.5).abs() < EPS
        && (v(2, 3) + 0.5).abs() < EPS
        && v(3, 1).abs() < EPS
        && (v(3, 2) - C07).abs() < EPS
        && (v(3, 3) - C07).abs() < EPS
    {
        return 0;
    } else if (v(1, 1) - 1.0).abs() < EPS
        && (v(2, 3) - 1.0).abs() < EPS
        && (v(3, 2) + 1.0).abs() < EPS
    {
        return 2; // front
    } else if (v(1, 1) - C07).abs() < EPS
        && (v(1, 2) - C07).abs() < EPS
        && v(1, 3).abs() < EPS
        && (v(2, 1) + 0.5).abs() < EPS
        && (v(2, 2) - 0.5).abs() < EPS
        && (v(2, 3) - C07).abs() < EPS
        && (v(3, 1) - 0.5).abs() < EPS
        && (v(3, 2) + 0.5).abs() < EPS
        && (v(3, 3) - C07).abs() < EPS
    {
        return 3; // axo
    }
    -1
}

/// Implements a projector object: transforms and projects points and vectors
/// for the removal of hidden lines. Source: `HLRAlgo_Projector`.
#[derive(Debug, Clone)]
pub struct HlrAlgoProjector {
    my_type: i32,
    persp: bool,
    focus: f64,
    scaled_trsf: GpTrsf,
    trsf: GpTrsf,
    inv_trsf: GpTrsf,
    d1: GpVec2d,
    d2: GpVec2d,
    d3: GpVec2d,
}

impl Default for HlrAlgoProjector {
    /// `HLRAlgo_Projector()` (`cxx:38-43`): no perspective, identity
    /// `myScaledTrsf`, then `Scaled()` (the `Result` cannot fail for the
    /// identity transform).
    fn default() -> Self {
        let mut p = Self {
            my_type: -1,
            persp: false,
            focus: 0.0,
            scaled_trsf: GpTrsf::identity(),
            trsf: GpTrsf::identity(),
            inv_trsf: GpTrsf::identity(),
            d1: GpVec2d::zero(),
            d2: GpVec2d::zero(),
            d3: GpVec2d::zero(),
        };
        let _ = p.scaled(false);
        p
    }
}

impl HlrAlgoProjector {
    /// `HLRAlgo_Projector()` (`cxx:38-43`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `HLRAlgo_Projector(const gp_Ax2&)` (`cxx:47-56`): axonometric projector
    /// with the viewing coordinate system `cs`. OCCT passes `gp_Ax2` where a
    /// `gp_Ax3` is expected, i.e. the implicit `gp_Ax3(gp_Ax2)` conversion.
    pub fn from_ax2(cs: &GpAx2) -> Result<Self, &'static str> {
        let mut p = Self {
            my_type: -1,
            persp: false,
            focus: 0.0,
            scaled_trsf: GpTrsf::identity(),
            trsf: GpTrsf::identity(),
            inv_trsf: GpTrsf::identity(),
            d1: GpVec2d::zero(),
            d2: GpVec2d::zero(),
            d3: GpVec2d::zero(),
        };
        p.scaled_trsf.set_transformation(&GpAx3::from_ax2(cs));
        p.scaled(false)?;
        p.set_direction();
        Ok(p)
    }

    /// `HLRAlgo_Projector(const gp_Ax2&, const double Focus)` (`cxx:58-66`):
    /// perspective projector.
    pub fn from_ax2_with_focus(cs: &GpAx2, focus: f64) -> Result<Self, &'static str> {
        let mut p = Self {
            my_type: -1,
            persp: true,
            focus,
            scaled_trsf: GpTrsf::identity(),
            trsf: GpTrsf::identity(),
            inv_trsf: GpTrsf::identity(),
            d1: GpVec2d::zero(),
            d2: GpVec2d::zero(),
            d3: GpVec2d::zero(),
        };
        p.scaled_trsf.set_transformation(&GpAx3::from_ax2(cs));
        p.scaled(false)?;
        p.set_direction();
        Ok(p)
    }

    /// `HLRAlgo_Projector(const gp_Trsf&, bool, double)` (`cxx:70-78`).
    pub fn from_trsf(t: &GpTrsf, persp: bool, focus: f64) -> Result<Self, &'static str> {
        let mut p = Self {
            my_type: -1,
            persp,
            focus,
            scaled_trsf: t.clone(),
            trsf: GpTrsf::identity(),
            inv_trsf: GpTrsf::identity(),
            d1: GpVec2d::zero(),
            d2: GpVec2d::zero(),
            d3: GpVec2d::zero(),
        };
        p.scaled(false)?;
        p.set_direction();
        Ok(p)
    }

    /// `HLRAlgo_Projector(const gp_Trsf&, bool, double, const gp_Vec2d&,
    /// const gp_Vec2d&, const gp_Vec2d&)` (`cxx:82-94`): the projector with
    /// given minmax directions.
    pub fn from_trsf_directions(
        t: &GpTrsf,
        persp: bool,
        focus: f64,
        v1: &GpVec2d,
        v2: &GpVec2d,
        v3: &GpVec2d,
    ) -> Result<Self, &'static str> {
        let mut p = Self {
            my_type: -1,
            persp,
            focus,
            scaled_trsf: t.clone(),
            trsf: GpTrsf::identity(),
            inv_trsf: GpTrsf::identity(),
            d1: *v1,
            d2: *v2,
            d3: *v3,
        };
        p.scaled(false)?;
        Ok(p)
    }

    /// `Set` (`cxx:98-106`).
    pub fn set(&mut self, t: &GpTrsf, persp: bool, focus: f64) -> Result<(), &'static str> {
        self.persp = persp;
        self.focus = focus;
        self.scaled_trsf = t.clone();
        self.scaled(false)?;
        self.set_direction();
        Ok(())
    }

    /// `Scaled` (`cxx:152-168`).
    pub fn scaled(&mut self, on: bool) -> Result<(), &'static str> {
        self.my_type = -1;
        self.trsf = self.scaled_trsf.clone();
        if !on {
            self.trsf.set_scale_factor(1.0)?;
            if !self.persp {
                self.trsf.set_translation_part(&GpVec::new(0.0, 0.0, 0.0));
                self.my_type = trsf_type(&self.trsf);
            }
        }
        self.inv_trsf = self.trsf.clone();
        self.inv_trsf.invert()?;
        Ok(())
    }

    /// `SetDirection` (`cxx:360-388`).
    fn set_direction(&mut self) {
        let mut v1 = GpVec::new(1.0, 0.0, 0.0);
        v1.transform(&self.trsf);
        if (v1.x().abs() + v1.y().abs()) < ANGULAR {
            v1 = GpVec::new(1.0, 1.0, 0.0);
        }
        let d1 = GpVec2d::new(v1.x(), v1.y());
        self.d1.set_coord(-d1.y(), d1.x());

        let mut v2 = GpVec::new(0.0, 1.0, 0.0);
        v2.transform(&self.trsf);
        if (v2.x().abs() + v2.y().abs()) < ANGULAR {
            v2 = GpVec::new(1.0, 1.0, 0.0);
        }
        let d2 = GpVec2d::new(v2.x(), v2.y());
        self.d2.set_coord(-d2.y(), d2.x());

        let mut v3 = GpVec::new(0.0, 0.0, 1.0);
        v3.transform(&self.trsf);
        if (v3.x().abs() + v3.y().abs()) < ANGULAR {
            v3 = GpVec::new(1.0, 1.0, 0.0);
        }
        let d3 = GpVec2d::new(v3.x(), v3.y());
        self.d3.set_coord(-d3.y(), d3.x());
    }

    /// `Directions` (`lxx:24-29`): `(D1, D2, D3)`.
    pub fn directions(&self) -> (GpVec2d, GpVec2d, GpVec2d) {
        (self.d1, self.d2, self.d3)
    }

    /// `Perspective` (`lxx:33-36`).
    pub fn perspective(&self) -> bool {
        self.persp
    }

    /// `Transformation` (`cxx:392-396`).
    pub fn transformation(&self) -> &GpTrsf {
        &self.trsf
    }

    /// `InvertedTransformation` (`lxx:46-49`).
    pub fn inverted_transformation(&self) -> &GpTrsf {
        &self.inv_trsf
    }

    /// `FullTransformation` (`lxx:53-56`).
    pub fn full_transformation(&self) -> &GpTrsf {
        &self.scaled_trsf
    }

    /// `Focus` (`lxx:58-63`): OCCT raises `Standard_NoSuchObject` when the
    /// projector is not a perspective one.
    pub fn focus(&self) -> Result<f64, &'static str> {
        if !self.persp {
            return Err("HlrAlgoProjector::Not a Perspective");
        }
        Ok(self.focus)
    }

    /// `Transform(gp_Vec&)` (`lxx:66-69`).
    pub fn transform_vec(&self, d: &mut GpVec) {
        d.transform(&self.trsf);
    }

    /// `Transform(gp_Pnt&)` (`lxx:72-75`).
    pub fn transform_pnt(&self, p: &mut GpPnt) {
        p.transform(&self.trsf);
    }

    /// `Project(const gp_Pnt&, gp_Pnt2d&)` (`cxx:172-238`): transform and
    /// apply perspective if needed.
    pub fn project(&self, p: &GpPnt, pout: &mut GpPnt2d) {
        match self.my_type {
            0 => {
                // axono standard
                let x07 = p.x() * 0.7071067811865475;
                let y05 = p.y() * 0.5;
                let z05 = p.z() * 0.5;
                pout.set_coord(x07 - y05 + z05, x07 + y05 - z05);
            }
            1 => {
                // top
                pout.set_coord(p.x(), p.y());
            }
            2 => {
                pout.set_coord(p.x(), p.z());
            }
            3 => {
                let xmy05 = (p.x() - p.y()) * 0.5;
                let z07 = p.z() * 0.7071067811865476;
                let x = 0.7071067811865476 * (p.x() + p.y());
                let y = -xmy05 + z07;
                pout.set_coord(x, y);
            }
            _ => {
                let mut p2 = *p;
                self.transform_pnt(&mut p2);
                if self.persp {
                    let r = 1.0 - p2.z() / self.focus;
                    pout.set_coord(p2.x() / r, p2.y() / r);
                } else {
                    pout.set_coord(p2.x(), p2.y());
                }
            }
        }
    }

    /// `Project(const gp_Pnt&, double&, double&, double&)` (`cxx:262-317`):
    /// returns `(X, Y, Z)`.
    pub fn project_coords(&self, p: &GpPnt) -> (f64, f64, f64) {
        match self.my_type {
            0 => {
                // axono standard
                let x07 = p.x() * 0.7071067811865475;
                let y05 = p.y() * 0.5;
                let z05 = p.z() * 0.5;
                (
                    x07 - y05 + z05,
                    x07 + y05 - z05,
                    0.7071067811865475 * (p.y() + p.z()),
                )
            }
            1 => (p.x(), p.y(), p.z()),
            2 => (p.x(), p.z(), -p.y()),
            3 => {
                let xmy05 = (p.x() - p.y()) * 0.5;
                let z07 = p.z() * 0.7071067811865476;
                (
                    0.7071067811865476 * (p.x() + p.y()),
                    -xmy05 + z07,
                    xmy05 + z07,
                )
            }
            _ => {
                let mut p2 = *p;
                self.transform_pnt(&mut p2);
                let (mut x, mut y, z) = (p2.x(), p2.y(), p2.z());
                if self.persp {
                    let r = 1.0 - z / self.focus;
                    x /= r;
                    y /= r;
                }
                (x, y, z)
            }
        }
    }

    /// `Project(const gp_Pnt&, const gp_Vec&, gp_Pnt2d&, gp_Vec2d&)`
    /// (`cxx:322-343`): transform and apply perspective, with the
    /// derivative `D1out`. Returns `(Pout, D1out)`.
    pub fn project_with_derivative(&self, p: &GpPnt, d1: &GpVec) -> (GpPnt2d, GpVec2d) {
        let mut pp = *p;
        pp.transform(&self.trsf);
        let mut dd1 = *d1;
        dd1.transform(&self.trsf);
        if self.persp {
            let r = 1.0 - pp.z() / self.focus;
            let pout = GpPnt2d::new(pp.x() / r, pp.y() / r);
            let d1out = GpVec2d::new(
                dd1.x() / r + pp.x() * dd1.z() / (self.focus * r * r),
                dd1.y() / r + pp.y() * dd1.z() / (self.focus * r * r),
            );
            (pout, d1out)
        } else {
            (GpPnt2d::new(pp.x(), pp.y()), GpVec2d::new(dd1.x(), dd1.y()))
        }
    }

    /// `Shoot` (`cxx:347-358`): the line through the eye towards the 2d point
    /// `(X, Y)`. `gp_Dir` raises on a null direction, reported as `Err`.
    pub fn shoot(&self, x: f64, y: f64) -> Result<GpLin, &'static str> {
        let mut l = if self.persp {
            GpLin::from_pnt_dir(
                GpPnt::new(0.0, 0.0, self.focus),
                GpDir::new(x, y, -self.focus)?,
            )
        } else {
            // `gp_Dir::D::NZ` is `(0, 0, -1)`.
            GpLin::from_pnt_dir(GpPnt::new(x, y, 0.0), GpDir::new(0.0, 0.0, -1.0)?)
        };
        l.transform(&self.inv_trsf);
        Ok(l)
    }
}

