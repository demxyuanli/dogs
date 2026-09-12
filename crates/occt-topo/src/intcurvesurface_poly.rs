//! `IntCurveSurface_ThePolygonOfHInter` / `ThePolyhedronOfHInter`.
//!
//! Source: `IntCurveSurface_ThePolygonOfHInter.cxx`,
//! `IntCurveSurface_ThePolyhedronOfHInter.cxx`. These are the bounding
//! polylines / UV grids `HInter` builds before the exact CS solve. IntCurvesFace
//! uses the polyhedron AABB as `RejectFace`.

use occt_core::bnd::BndBox;
use occt_core::gp::{GpLin, GpPnt};
use occt_geom::{Curve, Surface};

/// Curve polyline (`IntCurveSurface_ThePolygonOfHInter`).
#[derive(Clone)]
pub struct ThePolygon {
    pub params: Vec<f64>,
    pub points: Vec<GpPnt>,
    pub box_: BndBox,
}

impl ThePolygon {
    /// Sample `curve` on `[t0, t1]` with `nb` segments (`NbPointsOnCurve`).
    pub fn of_curve(curve: &dyn Curve, t0: f64, t1: f64, nb: usize) -> Self {
        let n = nb.max(2);
        let mut params = Vec::with_capacity(n);
        let mut points = Vec::with_capacity(n);
        let mut box_ = BndBox::new();
        for i in 0..n {
            let t = t0 + (t1 - t0) * i as f64 / (n - 1) as f64;
            let p = curve.d0(t);
            params.push(t);
            box_.add_point(&p);
            points.push(p);
        }
        Self {
            params,
            points,
            box_,
        }
    }

    /// Line segment as a two-point polygon (`ThePolygonOfHInter` line ctor).
    pub fn of_line(lin: &GpLin, t0: f64, t1: f64) -> Self {
        let p0 = occt_core::elib::clib::line_value(lin, t0);
        let p1 = occt_core::elib::clib::line_value(lin, t1);
        let mut box_ = BndBox::new();
        box_.add_point(&p0);
        box_.add_point(&p1);
        Self {
            params: vec![t0, t1],
            points: vec![p0, p1],
            box_,
        }
    }

    pub fn nb_points(&self) -> usize {
        self.points.len()
    }
}

/// Surface UV grid (`IntCurveSurface_ThePolyhedronOfHInter`).
#[derive(Clone)]
pub struct ThePolyhedron {
    pub u: Vec<f64>,
    pub v: Vec<f64>,
    pub box_: BndBox,
}

impl ThePolyhedron {
    /// Uniform UV sample (`ThePolyhedronOfHInter(S, nbU, nbV, U1, V1, U2, V2)`).
    pub fn of_surface(
        surface: &dyn Surface,
        nbu: usize,
        nbv: usize,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) -> Self {
        let nu = nbu.max(2);
        let nv = nbv.max(2);
        let mut u = Vec::with_capacity(nu);
        let mut v = Vec::with_capacity(nv);
        for i in 0..nu {
            u.push(u1 + (u2 - u1) * i as f64 / (nu - 1) as f64);
        }
        for j in 0..nv {
            v.push(v1 + (v2 - v1) * j as f64 / (nv - 1) as f64);
        }
        let mut box_ = BndBox::new();
        for &ui in &u {
            for &vj in &v {
                box_.add_point(&surface.d0(ui, vj));
            }
        }
        Self { u, v, box_ }
    }

    /// Whether the curve polygon AABB is out of this polyhedron (`IsOut`).
    pub fn is_out_polygon(&self, poly: &ThePolygon) -> bool {
        if self.box_.is_void() || poly.box_.is_void() {
            return true;
        }
        self.box_.is_out_box(&poly.box_)
    }
}
