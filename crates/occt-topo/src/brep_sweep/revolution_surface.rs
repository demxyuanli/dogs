//! GeomAdaptor_SurfaceOfRevolution::Load and GetType.
//! Source: GeomAdaptor_SurfaceOfRevolution.cxx:100-199 and :370-555.

use std::sync::Arc;

use occt_core::gp::{
    GpAx1, GpAx2, GpAx3, GpCone, GpCylinder, GpDir, GpLin, GpPln, GpPnt, GpSphere, GpTorus, GpVec,
};
use occt_core::precision::{ANGULAR, CONFUSION};
use occt_geom::curve::Curve;
use occt_geom::{
    GeomCone, GeomCylinder, GeomPlane, GeomSphere, GeomSurfaceOfRevolution, GeomTorus, Surface,
};
/// `GeomAdaptor_SurfaceOfRevolution::Load` + `GetType` (`GeomAdaptor_SurfaceOfRevolution.cxx:100-199`, `:370-467`).
pub(super) fn revolution_surface(curve: Arc<dyn Curve>, first: f64, last: f64, axis: &GpAx1) -> Arc<dyn Surface> {
    let axe = revolution_ax3(curve.as_ref(), first, last, axis);
    match classify(curve.as_ref(), first, last, axis, &axe) {
        RevKind::Plane => {
            let pln = revolution_plane(curve.as_ref(), &axe, axis);
            Arc::new(GeomPlane::new(pln))
        }
        RevKind::Cylinder => {
            let p = curve.value(0.0);
            let r = GpVec::from_pnts(&axe.location(), &p).dot(&GpVec::new(
                axe.x_direction().x(),
                axe.x_direction().y(),
                axe.x_direction().z(),
            ));
            match GpCylinder::new(axe, r) {
                Ok(c) => Arc::new(GeomCylinder::new(c)),
                Err(_) => Arc::new(GeomSurfaceOfRevolution::new(curve, *axis)),
            }
        }
        RevKind::Sphere => {
            if let Some(circ) = curve.gp_circ() {
                let mut ax = axe;
                ax.set_location(circ.location());
                match GpSphere::new(ax, circ.radius()) {
                    Ok(s) => Arc::new(GeomSphere::new(s)),
                    Err(_) => Arc::new(GeomSurfaceOfRevolution::new(curve, *axis)),
                }
            } else {
                Arc::new(GeomSurfaceOfRevolution::new(curve, *axis))
            }
        }
        RevKind::Cone => {
            if let Some(line) = curve.gp_line() {
                let mut angle = axe.direction().angle(&line.direction());
                let p0 = curve.value(0.0);
                let r = axe.location().distance(&p0);
                if r >= CONFUSION {
                    let mut op0 = GpVec::from_pnts(&axe.location(), &p0);
                    let ldir = line.direction();
                    let denom = ldir.dot(axe.x_direction());
                    if denom.abs() > 0.0 {
                        let t = op0.dot(&GpVec::new(
                            axe.x_direction().x(),
                            axe.x_direction().y(),
                            axe.x_direction().z(),
                        )) / denom;
                        let ld = GpVec::new(ldir.x(), ldir.y(), ldir.z());
                        op0 = op0.subtracted(&ld.multiplied_scalar(t));
                        let ad = GpVec::new(axe.direction().x(), axe.direction().y(), axe.direction().z());
                        if op0.dot(&ad) > 0.0 {
                            angle = -angle;
                        }
                    }
                }
                let cone = GpCone {
                    pos: axe,
                    radius: r,
                    semi_angle: angle,
                };
                Arc::new(GeomCone::new(cone))
            } else {
                Arc::new(GeomSurfaceOfRevolution::new(curve, *axis))
            }
        }
        RevKind::Torus => {
            if let Some(circ) = curve.gp_circ() {
                let major = GpLin::from_pnt_dir(*axis.location(), *axis.direction()).distance(&circ.location());
                match GpTorus::new(axe, major, circ.radius()) {
                    Ok(t) => Arc::new(GeomTorus::new(t)),
                    Err(_) => Arc::new(GeomSurfaceOfRevolution::new(curve, *axis)),
                }
            } else {
                Arc::new(GeomSurfaceOfRevolution::new(curve, *axis))
            }
        }
        RevKind::Revolution => Arc::new(GeomSurfaceOfRevolution::new(curve, *axis)),
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum RevKind {
    Plane,
    Cylinder,
    Sphere,
    Cone,
    Torus,
    Revolution,
}

fn revolution_ax3(curve: &dyn Curve, first: f64, last: f64, axis: &GpAx1) -> GpAx3 {
    let mut oz = *axis.direction();
    let mut yrev = false;
    if let Some(line) = curve.gp_line() {
        if line.direction().dot(&oz) < 0.0 {
            yrev = true;
            oz.reverse();
        }
    }
    let p = curve.value(0.0);
    let q = if curve.gp_circ().is_some() {
        curve.gp_circ().unwrap().location()
    } else if classify(curve, first, last, axis, &GpAx3::default()) == RevKind::Cone {
        if let Some(line) = curve.gp_line() {
            if GpLin::from_pnt_dir(*axis.location(), *axis.direction()).distance(&p) <= CONFUSION {
                let loc = line.location();
                let d = line.direction();
                GpPnt::new(loc.x() + d.x(), loc.y() + d.y(), loc.z() + d.z())
            } else {
                p
            }
        } else {
            p
        }
    } else if !first.is_finite() {
        p
    } else {
        curve.value(first)
    };
    let dz = *axis.direction();
    let o0 = *axis.location();
    let shift = GpVec::from_pnts(&o0, &p).dot(&GpVec::new(dz.x(), dz.y(), dz.z()));
    let origin = GpPnt::new(
        o0.x() + shift * dz.x(),
        o0.y() + shift * dz.y(),
        o0.z() + shift * dz.z(),
    );
    let ox = if GpLin::from_pnt_dir(*axis.location(), *axis.direction()).distance(&q) > CONFUSION {
        GpDir::from_vec(&GpVec::from_pnts(&origin, &q)).unwrap_or(oz)
    } else {
        let mut ratio = 1.0;
        let mut pp = p;
        let mut dist = 0.0;
        while dist < CONFUSION && ratio < 100.0 {
            let t = first + (last - first) / ratio;
            pp = curve.value(t);
            dist = GpLin::from_pnt_dir(*axis.location(), *axis.direction()).distance(&pp);
            ratio += 1.0;
        }
        let ozv = GpVec::new(oz.x(), oz.y(), oz.z());
        let crossed = ozv.crossed(&GpVec::from_pnts(&origin, &pp));
        let oxv = crossed.crossed(&ozv);
        GpDir::from_vec(&oxv).unwrap_or(oz)
    };
    let mut axe = GpAx3::new(origin, oz, &ox).unwrap_or_else(|_| GpAx3::default());
    if yrev {
        axe.y_reverse();
    } else if let Some(circ) = curve.gp_circ() {
        if let Ok(cross) = ox.crossed(&oz) {
            if cross.dot(&circ.axis().direction()) < 0.0 {
                axe.z_reverse();
            }
        }
    }
    axe
}

fn classify(curve: &dyn Curve, first: f64, last: f64, axis: &GpAx1, axe: &GpAx3) -> RevKind {
    if let Some(line) = curve.gp_line() {
        let axe_dir = line.position().direction();
        if axis.direction().is_parallel_tol(&axe_dir, ANGULAR) {
            let p = curve.value(0.0);
            let r = GpVec::from_pnts(&axe.location(), &p).dot(&GpVec::new(
                axe.x_direction().x(),
                axe.x_direction().y(),
                axe.x_direction().z(),
            ));
            if r > CONFUSION {
                return RevKind::Cylinder;
            }
        } else if dir_is_normal(axis.direction(), &axe_dir, ANGULAR) {
            return RevKind::Plane;
        } else {
            let trimmed = first.is_finite() && last.is_finite();
            if trimmed {
                let pf = curve.value(first);
                let pl = curve.value(last);
                let len = pf.distance(&pl);
                let vlin = GpVec::from_pnts(&pf, &pl);
                let vaxe = GpVec::new(axis.direction().x(), axis.direction().y(), axis.direction().z());
                let projlen = vaxe.dot(&vlin).abs();
                if (len - projlen) <= CONFUSION {
                    let p = curve.value(0.0);
                    let r = GpVec::from_pnts(&axe.location(), &p).dot(&GpVec::new(
                        axe.x_direction().x(),
                        axe.x_direction().y(),
                        axe.x_direction().z(),
                    ));
                    if r > CONFUSION {
                        return RevKind::Cylinder;
                    }
                } else if projlen <= CONFUSION {
                    return RevKind::Plane;
                }
            }
            let line_loc = line.location();
            let v = GpVec::from_pnts(axis.location(), &line_loc);
            let w = GpVec::new(axe_dir.x(), axe_dir.y(), axe_dir.z());
            let axis_dir = GpVec::new(axis.direction().x(), axis.direction().y(), axis.direction().z());
            let proj = w.dot(&axis_dir).abs();
            if v.dot(&axis_dir.crossed(&w)).abs() <= CONFUSION
                && proj >= CONFUSION
                && proj <= 1.0 - CONFUSION
            {
                return RevKind::Cone;
            }
        }
    } else if let Some(circ) = curve.gp_circ() {
        if !ax2_coplanar_ax1(&circ.position(), axis, CONFUSION, ANGULAR) {
            return RevKind::Revolution;
        }
        let lin = GpLin::from_pnt_dir(*axis.location(), *axis.direction());
        if lin.distance(&circ.location()) <= CONFUSION {
            return RevKind::Sphere;
        }
        if lin.distance(&circ.location()) > circ.radius() {
            return RevKind::Torus;
        }
    }
    RevKind::Revolution
}

fn revolution_plane(curve: &dyn Curve, axe: &GpAx3, axis: &GpAx1) -> GpPln {
    let mut ax = axe.clone();
    let pon = curve.value(0.0);
    let dot = GpVec::from_pnts(axis.location(), &pon).dot(&GpVec::new(
        axis.direction().x(),
        axis.direction().y(),
        axis.direction().z(),
    ));
    let p = GpPnt::new(
        axis.location().x() + dot * axis.direction().x(),
        axis.location().y() + dot * axis.direction().y(),
        axis.location().z() + dot * axis.direction().z(),
    );
    ax.set_location(p);
    if let Some(line) = curve.gp_line() {
        if ax.x_direction().dot(&line.direction()) >= -CONFUSION {
            ax.x_reverse();
        }
    }
    GpPln::new(ax)
}

fn ax2_coplanar_ax1(ax2: &GpAx2, ax1: &GpAx1, lin_tol: f64, ang_tol: f64) -> bool {
    let dd = ax2.direction();
    let pp = ax2.location();
    let ap = ax1.location();
    let d1 = (dd.x() * (ap.x() - pp.x()) + dd.y() * (ap.y() - pp.y()) + dd.z() * (ap.z() - pp.z())).abs();
    d1 <= lin_tol && dir_is_normal(&dd, ax1.direction(), ang_tol)
}

fn dir_is_normal(a: &GpDir, b: &GpDir, ang_tol: f64) -> bool {
    (std::f64::consts::FRAC_PI_2 - a.angle(b).abs()).abs() <= ang_tol
}
