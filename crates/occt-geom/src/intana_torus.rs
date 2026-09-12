//! `IntAna_QuadQuadGeo` torus pairs.
//! Source: `IntAna_QuadQuadGeo.cxx` Perform(Pln/Cyl/Cone/Sphere/Torus, Torus).

use std::f64::consts::{FRAC_PI_2, PI};

use occt_core::gp::{
    GpAx1, GpAx2, GpCirc, GpCone, GpCylinder, GpDir, GpLin, GpPln, GpPnt, GpSphere, GpTorus, GpVec,
    GpXyz,
};
use occt_core::precision::ANGULAR;

const EPS_DISTANCE: f64 = 1.0e-14;
const EPS_CYL_DELTA_RADIUS: f64 = 1.0e-13;

/// Closed-form torus intersection. `Fail` is OCCT `IntAna_NoGeometricSolution`.
#[derive(Debug, Clone)]
pub enum TorusIntersection {
    Fail,
    Empty,
    Same,
    Circles(Vec<GpCirc>),
}

/// `DirToAx2` — Ax2 whose Z is `d`.
fn dir_to_ax2(p: GpPnt, d: &GpDir) -> GpAx2 {
    let (x, y, z) = (d.x(), d.y(), d.z());
    let (ax, ay, az) = (x.abs(), y.abs(), z.abs());
    let xdir = if ax == 0.0 || (ax < ay && ax < az) {
        GpDir::new(0.0, -z, y)
    } else if ay == 0.0 || (ay < ax && ay < az) {
        GpDir::new(-z, 0.0, x)
    } else {
        GpDir::new(-y, x, 0.0)
    };
    let xdir = xdir.unwrap_or_else(|_| GpDir::new(1.0, 0.0, 0.0).unwrap());
    GpAx2::new(p, *d, xdir).unwrap_or_else(|_| GpAx2::standard())
}

fn circle_at(center: GpPnt, dir: &GpDir, radius: f64) -> GpCirc {
    GpCirc::new(dir_to_ax2(center, dir), radius)
}

fn dirs_parallel(a: &GpDir, b: &GpDir, ang: f64) -> bool {
    let a0 = a.angle(b);
    a0 <= ang || (PI - a0) <= ang
}

fn dirs_normal(a: &GpDir, b: &GpDir, ang: f64) -> bool {
    (FRAC_PI_2 - a.angle(b)).abs() <= ang
}

fn plane_coeffs(p: &GpPln) -> (f64, f64, f64, f64) {
    let n = p.pos.direction();
    let loc = p.location();
    let (a, b, c) = (n.x(), n.y(), n.z());
    let d = -(a * loc.x() + b * loc.y() + c * loc.z());
    (a, b, c, d)
}

fn plane_distance(pln: &GpPln, p: &GpPnt) -> f64 {
    let (a, b, c, d) = plane_coeffs(pln);
    (a * p.x() + b * p.y() + c * p.z() + d).abs()
}

fn rotate_line(lin: &GpLin, ax: &GpAx1, angle: f64) -> GpLin {
    let loc = lin.location().rotated(ax, angle);
    let tip = lin
        .location()
        .translated_vec(&GpVec::from_xyz(lin.direction().xyz()))
        .rotated(ax, angle);
    let dir = GpDir::from_vec(&GpVec::from_pnts(&loc, &tip)).unwrap_or(lin.direction());
    GpLin::from_pnt_dir(loc, dir)
}

fn pnt_from_xyz(xyz: &GpXyz) -> GpPnt {
    GpPnt::from_xyz(xyz)
}

/// Plane ∩ torus. Port of `IntAna_QuadQuadGeo::Perform(gp_Pln, gp_Torus)`.
pub fn quadric_quadric_plane_torus(pln: &GpPln, tor: &GpTorus, tol: f64) -> TorusIntersection {
    let r_min = tor.minor_radius();
    let r_maj = tor.major_radius();
    if r_min >= r_maj {
        return TorusIntersection::Fail;
    }
    let pln_ax = pln.axis();
    let tor_ax = tor.axis();
    let parallel = dirs_parallel(pln_ax.direction(), tor_ax.direction(), ANGULAR);
    let normal = if parallel {
        false
    } else {
        dirs_normal(pln_ax.direction(), tor_ax.direction(), ANGULAR)
    };
    if !normal && !parallel {
        return TorusIntersection::Fail;
    }
    let tor_loc = *tor_ax.location();
    if parallel {
        let (a, b, c, d) = plane_coeffs(pln);
        let dist = a * tor_loc.x() + b * tor_loc.y() + c * tor_loc.z() + d;
        let a_dr = dist.abs() - r_min;
        if a_dr > EPS_CYL_DELTA_RADIUS {
            return TorusIntersection::Empty;
        }
        let dist = if a_dr.abs() < EPS_CYL_DELTA_RADIUS {
            if dist < 0.0 {
                -r_min
            } else {
                r_min
            }
        } else {
            dist
        };
        let a_dt = (r_min * r_min - dist * dist).abs().sqrt();
        let pt1 = GpPnt::new(
            tor_loc.x() - dist * a,
            tor_loc.y() - dist * b,
            tor_loc.z() - dist * c,
        );
        let dir1 = *tor_ax.direction();
        let mut circs = vec![circle_at(pt1, &dir1, r_maj + a_dt)];
        if a_dr < -EPS_CYL_DELTA_RADIUS && a_dt > tol {
            circs.push(circle_at(pt1, &dir1, r_maj - a_dt));
        }
        TorusIntersection::Circles(circs)
    } else {
        if plane_distance(pln, &tor_loc) > EPS_DISTANCE {
            return TorusIntersection::Fail;
        }
        let dir1 = *pln_ax.direction();
        let a_dir = tor_ax.direction().crossed(&dir1);
        let Ok(a_dir) = a_dir else {
            return TorusIntersection::Fail;
        };
        let offset = GpVec::from_xyz(a_dir.xyz()).multiplied_scalar(r_maj);
        let pt1 = tor_loc.translated_vec(&offset);
        let pt2 = tor_loc.translated_vec(&offset.reversed());
        TorusIntersection::Circles(vec![
            circle_at(pt1, &dir1, r_min),
            circle_at(pt2, &dir1, r_min),
        ])
    }
}

/// Cylinder ∩ torus. Port of `IntAna_QuadQuadGeo::Perform(gp_Cylinder, gp_Torus)`.
pub fn quadric_quadric_cylinder_torus(
    cyl: &GpCylinder,
    tor: &GpTorus,
    tol: f64,
) -> TorusIntersection {
    let r_min = tor.minor_radius();
    let r_maj = tor.major_radius();
    if r_min >= r_maj {
        return TorusIntersection::Fail;
    }
    let cyl_ax = cyl.axis();
    let tor_ax = tor.axis();
    let lin = GpLin::new(*tor_ax);
    if !dirs_parallel(tor_ax.direction(), cyl_ax.direction(), ANGULAR)
        || lin.distance(&cyl.location()) > EPS_DISTANCE
    {
        return TorusIntersection::Fail;
    }
    let r_cyl = cyl.radius();
    if (r_cyl + tol) < (r_maj - r_min) || (r_cyl - tol) > (r_maj + r_min) {
        return TorusIntersection::Empty;
    }
    let a_dist = (r_min * r_min - (r_cyl - r_maj) * (r_cyl - r_maj)).abs().sqrt();
    let dir1 = *tor_ax.direction();
    let tor_loc = tor_ax.location().coord;
    let pt1 = pnt_from_xyz(&tor_loc.added(&dir1.xyz().multiplied(a_dist)));
    let mut circs = vec![circle_at(pt1, &dir1, r_cyl)];
    if a_dist > tol && r_cyl > (r_maj - r_min) && r_cyl < (r_maj + r_min) {
        let pt2 = pnt_from_xyz(&tor_loc.added(&dir1.xyz().multiplied(-a_dist)));
        circs.push(circle_at(pt2, &dir1, r_cyl));
    }
    TorusIntersection::Circles(circs)
}

/// Cone ∩ torus. Port of `IntAna_QuadQuadGeo::Perform(gp_Cone, gp_Torus)`.
pub fn quadric_quadric_cone_torus(cone: &GpCone, tor: &GpTorus, tol: f64) -> TorusIntersection {
    let r_min = tor.minor_radius();
    let r_maj = tor.major_radius();
    if r_min >= r_maj {
        return TorusIntersection::Fail;
    }
    let con_ax = cone.axis();
    let tor_ax = tor.axis();
    let lin = GpLin::new(*tor_ax);
    let apex = cone.apex();
    if !dirs_parallel(tor_ax.direction(), con_ax.direction(), ANGULAR)
        || lin.distance(&apex) > EPS_DISTANCE
    {
        return TorusIntersection::Fail;
    }
    let angle = cone.semi_angle();
    let tor_loc = *tor_ax.location();
    let pn = tor_loc.translated_vec(
        &GpVec::from_xyz(tor.y_axis().direction().xyz()).multiplied_scalar(r_maj),
    );
    let Ok(dn) = GpDir::from_vec(&GpVec::from_pnts(&tor_loc, &pn)) else {
        return TorusIntersection::Fail;
    };
    let ax_rot = GpAx1::new(apex, dn);
    let con_l = rotate_line(&lin, &ax_rot, angle);
    let dl = con_l.direction();
    let mut x_dir = *tor.x_axis().direction();
    let mut pts = Vec::new();
    let mut params = Vec::new();
    for i in 0..2 {
        if i == 1 {
            x_dir = x_dir.reversed();
        }
        let pct = tor_loc.translated_vec(&GpVec::from_xyz(x_dir.xyz()).multiplied_scalar(r_maj));
        let a_dist = con_l.distance(&pct);
        if a_dist > r_min + tol {
            continue;
        }
        let ph = pct
            .coord
            .subtracted(&con_l.normal(&pct).direction().xyz().multiplied(a_dist));
        let a_dt = (r_min * r_min - a_dist * a_dist).abs().sqrt();
        let d_val = dl.xyz().multiplied(a_dt);
        let p = pnt_from_xyz(&ph.added(&d_val));
        params.push(lin.distance(&p));
        pts.push(pnt_from_xyz(
            &p.coord.subtracted(&x_dir.xyz().multiplied(*params.last().unwrap())),
        ));
        if a_dist < r_min && a_dt > tol {
            let p = pnt_from_xyz(&ph.subtracted(&d_val));
            params.push(lin.distance(&p));
            pts.push(pnt_from_xyz(
                &p.coord
                    .subtracted(&x_dir.xyz().multiplied(*params.last().unwrap())),
            ));
        }
    }
    if pts.is_empty() {
        return TorusIntersection::Empty;
    }
    let dir = *tor_ax.direction();
    let circs = pts
        .into_iter()
        .zip(params)
        .map(|(pt, r)| circle_at(pt, &dir, r))
        .collect();
    TorusIntersection::Circles(circs)
}

/// Sphere ∩ torus. Port of `IntAna_QuadQuadGeo::Perform(gp_Sphere, gp_Torus)`.
pub fn quadric_quadric_sphere_torus(
    sph: &GpSphere,
    tor: &GpTorus,
    tol: f64,
) -> TorusIntersection {
    let r_min = tor.minor_radius();
    let r_maj = tor.major_radius();
    if r_min >= r_maj {
        return TorusIntersection::Fail;
    }
    let tor_ax = tor.axis();
    let lin = GpLin::new(*tor_ax);
    let sph_loc = sph.location();
    if lin.distance(&sph_loc) > EPS_DISTANCE {
        return TorusIntersection::Fail;
    }
    let x_dir = *tor.x_axis().direction();
    let tor_loc_tube =
        pnt_from_xyz(&tor_ax.location().coord.added(&x_dir.xyz().multiplied(r_maj)));
    let r_sph = sph.radius();
    let vec12 = GpVec::from_pnts(&tor_loc_tube, &sph_loc);
    let a_dist = vec12.magnitude();
    if (a_dist - tol) > (r_min + r_sph) || (a_dist + tol) < (r_min - r_sph).abs() {
        return TorusIntersection::Empty;
    }
    if a_dist <= f64::EPSILON {
        if (r_min - r_sph).abs() > tol {
            return TorusIntersection::Empty;
        }
        let dir = *tor_ax.direction();
        return TorusIntersection::Circles(vec![circle_at(tor_loc_tube, &dir, r_min)]);
    }
    let an_alpha = 0.5 * (r_min * r_min - r_sph * r_sph + a_dist * a_dist) / a_dist;
    let a_beta = (r_min * r_min - an_alpha * an_alpha).abs().sqrt();
    let Ok(dir12) = GpDir::from_vec(&vec12) else {
        return TorusIntersection::Empty;
    };
    let ph = tor_loc_tube
        .coord
        .added(&dir12.xyz().multiplied(an_alpha));
    let Ok(dc) = tor.y_axis().direction().crossed(&dir12) else {
        return TorusIntersection::Fail;
    };
    let d_val = dc.xyz().multiplied(a_beta);
    let p = pnt_from_xyz(&ph.added(&d_val));
    let param1 = lin.distance(&p);
    let pt1 = pnt_from_xyz(&p.coord.subtracted(&x_dir.xyz().multiplied(param1)));
    let dir1 = *tor_ax.direction();
    let mut circs = vec![circle_at(pt1, &dir1, param1)];
    if a_dist < (r_sph + r_min) && a_dist > (r_sph - r_min).abs() && d_val.modulus() > tol {
        let p = pnt_from_xyz(&ph.subtracted(&d_val));
        let param2 = lin.distance(&p);
        let pt2 = pnt_from_xyz(&p.coord.subtracted(&x_dir.xyz().multiplied(param2)));
        circs.push(circle_at(pt2, &dir1, param2));
    }
    TorusIntersection::Circles(circs)
}

/// Torus ∩ torus. Port of `IntAna_QuadQuadGeo::Perform(gp_Torus, gp_Torus)`.
pub fn quadric_quadric_torus_torus(t1: &GpTorus, t2: &GpTorus, tol: f64) -> TorusIntersection {
    let r_min1 = t1.minor_radius();
    let r_maj1 = t1.major_radius();
    let r_min2 = t2.minor_radius();
    let r_maj2 = t2.major_radius();
    let ax1 = t1.axis();
    let ax2 = t2.axis();
    let loc1 = *ax1.location();
    let loc2 = *ax2.location();
    let l1 = GpLin::new(*ax1);
    if !dirs_parallel(ax1.direction(), ax2.direction(), ANGULAR) || l1.distance(&loc2) > EPS_DISTANCE
    {
        return TorusIntersection::Fail;
    }
    if loc1.distance(&loc2) <= tol
        && (r_min1 - r_min2).abs() <= tol
        && (r_maj1 - r_maj2).abs() <= tol
    {
        return TorusIntersection::Same;
    }
    if r_min1 >= r_maj1 || r_min2 >= r_maj2 {
        return TorusIntersection::Fail;
    }
    let x_dir1 = *t1.x_axis().direction();
    let p1 = pnt_from_xyz(&loc1.coord.added(&x_dir1.xyz().multiplied(r_maj1)));
    let p2 = pnt_from_xyz(&loc2.coord.added(&x_dir1.xyz().multiplied(r_maj2)));
    let v12 = GpVec::from_pnts(&p1, &p2);
    let a_dist = v12.magnitude();
    if (a_dist - tol) > (r_min1 + r_min2) || (a_dist + tol) < (r_min1 - r_min2).abs() {
        return TorusIntersection::Empty;
    }
    if a_dist <= f64::EPSILON {
        if (r_min1 - r_min2).abs() > tol {
            return TorusIntersection::Empty;
        }
        let dir = *ax1.direction();
        return TorusIntersection::Circles(vec![circle_at(p1, &dir, r_min1)]);
    }
    let an_alpha = 0.5 * (r_min1 * r_min1 - r_min2 * r_min2 + a_dist * a_dist) / a_dist;
    let a_beta = (r_min1 * r_min1 - an_alpha * an_alpha).abs().sqrt();
    let Ok(dir12) = GpDir::from_vec(&v12) else {
        return TorusIntersection::Empty;
    };
    let ph = p1.coord.added(&dir12.xyz().multiplied(an_alpha));
    let Ok(dc) = t1.y_axis().direction().crossed(&dir12) else {
        return TorusIntersection::Fail;
    };
    let d_val = dc.xyz().multiplied(a_beta);
    let p = pnt_from_xyz(&ph.added(&d_val));
    let param1 = l1.distance(&p);
    let pt1 = pnt_from_xyz(&p.coord.subtracted(&x_dir1.xyz().multiplied(param1)));
    let dir1 = *ax1.direction();
    let mut circs = vec![circle_at(pt1, &dir1, param1)];
    if a_dist < (r_min1 + r_min2) && a_dist > (r_min1 - r_min2).abs() && d_val.modulus() > tol {
        let p = pnt_from_xyz(&ph.subtracted(&d_val));
        let param2 = l1.distance(&p);
        let pt2 = pnt_from_xyz(&p.coord.subtracted(&x_dir1.xyz().multiplied(param2)));
        circs.push(circle_at(pt2, &dir1, param2));
    }
    TorusIntersection::Circles(circs)
}
