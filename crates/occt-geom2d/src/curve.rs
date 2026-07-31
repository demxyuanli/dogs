//! Abstract 2D parametric curve trait. Source: `Geom2d_Curve.hxx`
use occt_core::gp::{GpPnt2d, GpVec2d, GpTrsf2d};

/// Parametric 2D curve. Replaces OCCT Geom2d_Curve.
pub trait Curve2d: Send + Sync {
    fn d0(&self, u: f64) -> GpPnt2d;
    fn d1(&self, u: f64) -> (GpPnt2d, GpVec2d);
    fn d2(&self, u: f64) -> (GpPnt2d, GpVec2d, GpVec2d);
    fn value(&self, u: f64) -> GpPnt2d { self.d0(u) }
    fn first_parameter(&self) -> f64;
    fn last_parameter(&self) -> f64;
    fn is_periodic(&self) -> bool { false }
    fn period(&self) -> f64 { 0.0 }
    fn continuity(&self) -> u8;
    fn transform(&mut self, t: &GpTrsf2d);
    fn reverse(&mut self);

    fn clone_dyn(&self) -> Box<dyn Curve2d>;

    fn transformed(&self, t: &GpTrsf2d) -> Box<dyn Curve2d> {
        let mut c = self.clone_dyn();
        c.transform(t);
        c
    }
    fn reversed(&self) -> Box<dyn Curve2d> {
        let mut c = self.clone_dyn();
        c.reverse();
        c
    }
}

// ponytail: Arc<dyn Curve2d> works for cloning without custom Clone impls
