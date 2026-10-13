//! `gp_Ax3` frame and circle helpers shared by the surface and revolution emitters (chain to the `GeomToIGES_GeomSurface` builder frames).
use super::prelude::*;

/// Finite sampling bounds of a surface (unbounded ranges clamp to ±1).
pub(super) fn surf_bounds(s: &dyn Surface) -> (f64, f64, f64, f64) {
    let (u0, u1) = s.u_range();
    let (v0, v1) = s.v_range();
    let clamp = |a: f64, b: f64| {
        if a.is_finite() && b.is_finite() && b > a {
            (a, b)
        } else {
            (-1.0, 1.0)
        }
    };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

/// `Precision::IsNegativeInfinite` / `IsPositiveInfinite` clamp of a parameter
/// range before it is used as a curve parameter
/// (`GeomToIGES_GeomSurface.cxx:806-813` cone, `:1058-1065` tabulated
/// cylinder, `:1144-1150` revolution).
pub(super) fn clamp_infinite_range(v0: f64, v1: f64) -> (f64, f64) {
    let a = if occt_core::precision::Precision::is_negative_infinite(v0) {
        -occt_core::precision::INFINITE
    } else {
        v0
    };
    let b = if occt_core::precision::Precision::is_positive_infinite(v1) {
        occt_core::precision::INFINITE
    } else {
        v1
    };
    (a, b)
}

/// The `gp_Ax3(gp_Ax2(...))` frames the elementary-surface transfers build their
/// generatrix circles on:
///
/// - sphere: `gp_Ax2(gp::Origin(), -gp::DY(), gp::DX())`
///   (`GeomToIGES_GeomSurface.cxx:893`);
/// - torus: `gp_Ax2(gp_Pnt(MajorRadius, 0., 0.), -gp::DY(), gp::DX())`
///   (`cxx:960`), i.e. the same frame shifted to the tube centre.
///
/// The circle is written in this frame (entity 100 + its own 124, see
/// [`IgesWriter::emit_circular_arc`]) and referenced by the 120.
pub(super) fn revolution_generatrix_frame(location: GpPnt) -> GpAx3 {
    // `-gp::DY()` is `gp::DY()` reversed (`gp_Dir::Reverse`), so its zero
    // components carry a negative sign; building `(0, -1, 0)` from literals
    // loses them and the 124's signed zeros no longer match OCCT. The 3-arg
    // `gp_Ax2(N, Vx)` ctor (`gp_Ax2.hxx:73-80`) sets `vxdir = N ^ (Vx ^ N)`
    // (`gp_Dir::CrossCross`) and `vydir = N ^ vxdir`; `GpAx3::new` computes the
    // same `vydir` from the orthogonalised `vxdir`.
    let n = GpDir::from_axis(occt_core::gp::dir::DirAxis::Y).reversed();
    let vx = GpDir::from_axis(occt_core::gp::dir::DirAxis::X);
    let inner = vx.xyz().crossed(n.xyz());
    let xdir = GpDir::from_xyz(&n.xyz().crossed(&inner)).unwrap_or(vx);
    GpAx3::new(location, n, &xdir).expect("axis2 placement (Origin, -DY, DX)")
}

/// `ElCLib::CircleValue` on `frame` (`clib::circle_value`).
pub(super) fn circle_point(frame: &GpAx3, radius: f64, u: f64) -> GpPnt {
    let p = frame
        .location()
        .coord
        .added(&frame.x_direction().xyz().multiplied(radius * u.cos()))
        .added(&frame.y_direction().xyz().multiplied(radius * u.sin()));
    GpPnt::from_xyz(&p)
}

/// Inverse of `IGESConvGeom_GeomBuilder::SetPosition(gp_Ax3)`
/// (`GeomBuilder.cxx:137-143`, `gp_Trsf::SetTransformation` then
/// `gp_Trsf::Invert` at `gp_Trsf.cxx:418-425`).
pub(super) fn builder_local_trsf(frame: &GpAx3) -> GpTrsf {
    let mut ps = GpTrsf::identity();
    ps.set_transformation_from_to(frame, &xoy_frame());
    ps.scale = 1.0 / ps.scale;
    ps.matrix = ps.matrix.transpose();
    ps.loc = ps.matrix.multiplied(&ps.loc).multiply(-ps.scale);
    ps
}

pub(super) fn eval_xyz(local: &GpTrsf, p: &GpPnt) -> (f64, f64, f64) {
    let q = p.transformed(local);
    (q.x(), q.y(), q.z())
}

/// `gp::XOY()` — the standard `gp_Ax3` frame (origin, `+Z` main, `+X`).
pub(super) fn xoy_frame() -> GpAx3 {
    GpAx3::new(
        GpPnt::zero(),
        GpDir::from_axis(occt_core::gp::dir::DirAxis::Z),
        &GpDir::from_axis(occt_core::gp::dir::DirAxis::X),
    )
    .expect("gp::XOY()")
}

/// 3×3 determinant.
pub(super) fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Circumcenter of three non-collinear 3D points (Cramer's rule on the
/// perpendicular-bisector system).
pub(super) fn circle_center3(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
    let d1 = GpVec::from_pnts(a, b);
    let d2 = GpVec::from_pnts(a, c);
    let n = d1.xyz().crossed(d2.xyz());
    if n.square_modulus() < 1e-30 {
        return None;
    }
    let n2 = |p: &GpPnt| p.x() * p.x() + p.y() * p.y() + p.z() * p.z();
    let mat = [
        [d1.x(), d1.y(), d1.z()],
        [d2.x(), d2.y(), d2.z()],
        [n.x, n.y, n.z],
    ];
    let rhs = [
        0.5 * (n2(b) - n2(a)),
        0.5 * (n2(c) - n2(a)),
        a.x() * n.x + a.y() * n.y + a.z() * n.z,
    ];
    let d = det3(&mat);
    if d.abs() < 1e-30 {
        return None;
    }
    let mut o = [0.0; 3];
    for k in 0..3 {
        let mut m = mat;
        for i in 0..3 {
            m[i][k] = rhs[i];
        }
        o[k] = det3(&m) / d;
    }
    Some(GpPnt::new(o[0], o[1], o[2]))
}

