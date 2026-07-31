//! Bidirectional conversions between gp primitives and Geom curves/surfaces.
//! Source: `Geom_Line.hxx`, `Geom_Plane.hxx`, `Geom_Circle.hxx`.

use std::sync::Arc;
use crate::curve::Curve;
use crate::surface::Surface;
use crate::{GeomCircle, GeomLine, GeomPlane};
use occt_core::gp::{GpCirc, GpLin, GpPln, GpPnt, GpXyz};

/// Wrap a gp line as a `Geom_Line`.
pub fn lin_to_geom(l: &GpLin) -> GeomLine {
    GeomLine::new(*l)
}

/// Unwrap the underlying gp line from a `Geom_Line`.
pub fn geom_to_lin(c: &GeomLine) -> GpLin {
    *c.lin()
}

/// Wrap a gp plane as a `Geom_Plane`.
pub fn pln_to_geom(p: &GpPln) -> GeomPlane {
    GeomPlane::new(p.clone())
}

/// Unwrap the underlying gp plane from a `Geom_Plane`.
pub fn geom_to_pln(s: &GeomPlane) -> GpPln {
    s.pln().clone()
}

/// Wrap a gp circle as a `Geom_Circle`.
pub fn circ_to_geom(c: &GpCirc) -> GeomCircle {
    GeomCircle::new(c.clone())
}

/// Unwrap the underlying gp circle from a `Geom_Circle`.
pub fn geom_to_circ(c: &GeomCircle) -> GpCirc {
    c.circ().clone()
}

/// Extract the coordinate triple from a point.
pub fn pnt_xyz(p: &GpPnt) -> GpXyz {
    p.coord
}

/// Build a point from a coordinate triple.
pub fn xyz_pnt(v: &GpXyz) -> GpPnt {
    GpPnt::from_xyz(v)
}

/// Box -> Arc for curves (the crate's handle type).
pub fn curve_as_handle(c: Box<dyn Curve>) -> Arc<dyn Curve> {
    Arc::from(c)
}

/// Box -> Arc for surfaces (the crate's handle type).
pub fn surface_as_handle(s: Box<dyn Surface>) -> Arc<dyn Surface> {
    Arc::from(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use occt_core::gp::{GpAx2, GpDir};

    #[test]
    fn lin_round_trip() {
        let l = GpLin::from_pnt_dir(GpPnt::new(1., 2., 3.), GpDir::new(0., 1., 0.).unwrap());
        let l2 = geom_to_lin(&lin_to_geom(&l));
        assert_eq!(l, l2);
    }

    #[test]
    fn circ_round_trip() {
        let c = GpCirc::new(GpAx2::standard(), 2.5);
        let c2 = geom_to_circ(&circ_to_geom(&c));
        assert_eq!(c.pos, c2.pos);
        assert_eq!(c.radius, c2.radius);
    }

    #[test]
    fn pln_round_trip() {
        let p = GpPln::default();
        let p2 = geom_to_pln(&pln_to_geom(&p));
        assert_eq!(p.pos.location(), p2.pos.location());
        assert_eq!(p.pos.direction(), p2.pos.direction());
        assert_eq!(p.pos.x_direction(), p2.pos.x_direction());
        assert_eq!(p.pos.y_direction(), p2.pos.y_direction());
    }

    #[test]
    fn pnt_xyz_round_trip() {
        let p = GpPnt::new(1., 2., 3.);
        let v = pnt_xyz(&p);
        assert_eq!(v, GpXyz::new(1., 2., 3.));
        let p2 = xyz_pnt(&v);
        assert_eq!(p, p2);
    }

    #[test]
    fn handles() {
        let l = GpLin::from_pnt_dir(GpPnt::new(0., 0., 0.), GpDir::new(1., 0., 0.).unwrap());
        let h = curve_as_handle(Box::new(lin_to_geom(&l)));
        assert!(h.d0(1.0).x() > 0.9);
    }
}
