//! `gce_MakeCirc2d` (2D circumscribed circle).
//!
//! Source: `ModelingData/TKGeomBase/gce/gce_MakeCirc2d.cxx:31-113`, the
//! `(P1, P2, P3)` constructor. It is the only overload this port needs:
//! `GC_MakeArcOfCircle2d(P1, P2, P3)` (`GC_MakeArcOfCircle2d.cxx:31-47`) is the
//! single caller, from `Geom2dConvert_ApproxArcsSegments::makeCircle`
//! (`Geom2dConvert_ApproxArcsSegments.cxx:119-137`).
//!
//! Deviations, all confined to the error paths:
//! - `gp_Dir2d(X, Y)` raises `Standard_ConstructionError` on a zero vector in
//!   OCCT; here it is an `Err`.
//! - The `(P1, P2, P3)` constructor leaves `TheError` untouched when
//!   `IntAna2d_AnaIntersection` is not done or is empty (`cxx:97-111`), so
//!   `IsDone()` reads an uninitialised member. Those two paths answer `Err`
//!   here instead of a garbage status.

use occt_core::elib::clib2d::parameter_circ2d;
use occt_core::gp::{GpAx2d, GpAx22d, GpCirc2d, GpDir2d, GpLin2d, GpPnt2d};
use occt_core::intana2d::IntAna2dAnaIntersection;
use occt_core::precision::REAL_SMALL;

/// `gce_MakeCirc2d(const gp_Pnt2d& P1, const gp_Pnt2d& P2, const gp_Pnt2d& P3)`
/// (`gce_MakeCirc2d.cxx:31-113`).
pub fn gce_make_circ2d_3p(
    p1: &GpPnt2d,
    p2: &GpPnt2d,
    p3: &GpPnt2d,
) -> Result<GpCirc2d, &'static str> {
    let (x1, y1) = (p1.x(), p1.y());
    let (x2, y2) = (p2.x(), p2.y());
    let (x3, y3) = (p3.x(), p3.y());

    // `cxx:49-61`: all three points coincident -> zero-radius circle at `P1`
    // whose axis is `gp_Ax2d(P1, gp::DX2d())`.
    let d12 = p1.distance(p2);
    let d13 = p1.distance(p3);
    let d23 = p2.distance(p3);
    if d12 < REAL_SMALL && d13 < REAL_SMALL && d23 < REAL_SMALL {
        return Ok(GpCirc2d::new(
            GpAx22d::from_xdir(*p1, GpDir2d::default()),
            0.0,
        ));
    }

    // `cxx:62-90`: the two perpendicular bisectors, with the degenerate
    // replacements for a coincident pair.
    let dir12 = GpDir2d::new(y1 - y2, x2 - x1)?;
    let dir13 = GpDir2d::new(y1 - y3, x3 - x1)?;
    let mid12 = GpPnt2d::new(0.5 * (x1 + x2), 0.5 * (y1 + y2));
    let mid13 = GpPnt2d::new(0.5 * (x1 + x3), 0.5 * (y1 + y3));

    let mut l1 = GpLin2d::new(GpAx2d::new(GpPnt2d::zero(), GpDir2d::default()));
    let mut l2 = GpLin2d::new(GpAx2d::new(GpPnt2d::zero(), GpDir2d::default()));
    if d12 > REAL_SMALL {
        l1 = GpLin2d::new(GpAx2d::new(mid12, dir12));
    }
    if d13 > REAL_SMALL {
        l2 = GpLin2d::new(GpAx2d::new(mid13, dir13));
    }
    if d13 <= REAL_SMALL {
        l2 = GpLin2d::new(GpAx2d::new(*p1, dir12));
    } else if d12 <= REAL_SMALL {
        l1 = GpLin2d::new(GpAx2d::new(*p1, dir13));
    } else if d23 <= REAL_SMALL {
        l2 = GpLin2d::new(GpAx2d::new(*p1, dir12));
    }

    let mut intp = IntAna2dAnaIntersection::new();
    intp.perform_lin_lin(&l1, &l2);
    if !intp.is_done() || intp.is_empty() {
        return Err("gce_MakeCirc2d: perpendicular bisectors do not intersect");
    }

    // `cxx:91-110`.
    let p_int = *intp.point(1).value();
    let da = p1.distance(&p_int);
    let db = p2.distance(&p_int);
    let dc = p3.distance(&p_int);
    let (xc, yc) = (p_int.x(), p_int.y());
    let d1 = GpDir2d::new(x1 - xc, y1 - yc)?;
    let d2 = GpDir2d::new(xc - x3, yc - y3)?;
    // `gp_Ax22d(P, Vx, Vy)`: the X direction is `Vx`; the Y direction is rebuilt
    // from `Vx` with the side chosen by the sign of `Vx ^ Vy` (`gp_Ax22d.hxx:60-73`).
    let vydir = if d1.crossed(&d2) >= 0.0 {
        GpDir2d { x: -d1.y, y: d1.x }
    } else {
        GpDir2d { x: d1.y, y: -d1.x }
    };
    let mut circ = GpCirc2d::new(
        GpAx22d {
            point: p_int,
            vxdir: d1,
            vydir,
        },
        (da + db + dc) / 3.0,
    );

    // `cxx:102-107`: orient the circle so that the three parameters grow.
    let a1 = parameter_circ2d(&circ, p1);
    let a2 = parameter_circ2d(&circ, p2);
    let a3 = parameter_circ2d(&circ, p3);
    if a1 > a2 || a2 > a3 {
        circ = circ.reversed();
    }
    Ok(circ)
}
