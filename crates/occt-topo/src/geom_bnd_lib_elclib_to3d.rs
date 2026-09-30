//! Lift 2D elementary geometry into a 3D plane.
//!
//! Source: `ElCLib.cxx` `To3d` overloads (1339-1421). A `gp_Ax2` is the
//! embedding plane: `P3 = Loc + X * Xd + Y * Yd`.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::gp::{
    GpAx1, GpAx2, GpAx2d, GpAx22d, GpCirc, GpCirc2d, GpDir, GpDir2d, GpElips, GpElips2d, GpHypr,
    GpHypr2d, GpLin, GpLin2d, GpParab, GpParab2d, GpPnt, GpPnt2d, GpVec, GpVec2d,
};

/// `ElCLib::To3d(Pos, gp_Pnt2d)`.
pub fn pnt_to3d(pos: &GpAx2, p: &GpPnt2d) -> GpPnt {
    let xd = pos.x_direction();
    let yd = pos.y_direction();
    let loc = pos.location();
    GpPnt::new(
        p.x() * xd.x() + p.y() * yd.x() + loc.x(),
        p.x() * xd.y() + p.y() * yd.y() + loc.y(),
        p.x() * xd.z() + p.y() * yd.z() + loc.z(),
    )
}

/// `ElCLib::To3d(Pos, gp_Dir2d)` — `Xd * V.X + Yd * V.Y`, then renormalise.
pub fn dir_to3d(pos: &GpAx2, v: &GpDir2d) -> GpDir {
    let vec = vec_to3d(pos, &GpVec2d::from_dir2d(v));
    GpDir::from_vec(&vec).unwrap_or_else(|_| *pos.x_direction())
}

/// `ElCLib::To3d(Pos, gp_Vec2d)`.
pub fn vec_to3d(pos: &GpAx2, v: &GpVec2d) -> GpVec {
    let xd = pos.x_direction();
    let yd = pos.y_direction();
    let mut vx = GpVec::new(xd.x(), xd.y(), xd.z()).multiply_scalar(v.x());
    let vy = GpVec::new(yd.x(), yd.y(), yd.z()).multiply_scalar(v.y());
    vx = vx.add(&vy);
    vx
}

/// `ElCLib::To3d(Pos, gp_Ax2d)`.
pub fn ax2d_to3d(pos: &GpAx2, a: &GpAx2d) -> GpAx1 {
    let p = pnt_to3d(pos, a.location());
    let v = dir_to3d(pos, a.direction());
    GpAx1::new(p, v)
}

/// `ElCLib::To3d(Pos, gp_Ax22d)` — `Ax2(P, VX × VY, VX)`.
pub fn ax22d_to3d(pos: &GpAx2, a: &GpAx22d) -> GpAx2 {
    let p = pnt_to3d(pos, a.location());
    let vx = vec_to3d(pos, &GpVec2d::from_dir2d(a.x_direction()));
    let vy = vec_to3d(pos, &GpVec2d::from_dir2d(a.y_direction()));
    let z = vx.crossed(&vy);
    let zdir = GpDir::from_vec(&z).unwrap_or_else(|_| pos.direction());
    let xdir = GpDir::from_vec(&vx).unwrap_or_else(|_| *pos.x_direction());
    GpAx2::new(p, zdir, xdir).unwrap_or_else(|_| *pos)
}

/// `ElCLib::To3d(Pos, gp_Lin2d)`.
pub fn lin_to3d(pos: &GpAx2, l: &GpLin2d) -> GpLin {
    GpLin::new(ax2d_to3d(pos, l.position()))
}

/// `ElCLib::To3d(Pos, gp_Circ2d)`.
pub fn circ_to3d(pos: &GpAx2, c: &GpCirc2d) -> GpCirc {
    GpCirc::new(ax22d_to3d(pos, c.position()), c.radius())
}

/// `ElCLib::To3d(Pos, gp_Elips2d)`.
pub fn elips_to3d(pos: &GpAx2, e: &GpElips2d) -> GpElips {
    GpElips::new(ax22d_to3d(pos, &e.pos), e.major_radius, e.minor_radius)
}

/// `ElCLib::To3d(Pos, gp_Hypr2d)`.
pub fn hypr_to3d(pos: &GpAx2, h: &GpHypr2d) -> GpHypr {
    GpHypr {
        pos: ax22d_to3d(pos, &h.pos),
        major_radius: h.major_radius,
        minor_radius: h.minor_radius,
    }
}

/// `ElCLib::To3d(Pos, gp_Parab2d)`.
pub fn parab_to3d(pos: &GpAx2, prb: &GpParab2d) -> GpParab {
    GpParab {
        pos: ax22d_to3d(pos, &prb.pos),
        focal: prb.focal,
    }
}
