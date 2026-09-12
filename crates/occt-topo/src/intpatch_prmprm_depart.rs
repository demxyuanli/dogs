//! `IntPatch_PrmPrmIntersection::PointDepart`.

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;
use occt_geom::Surface;

use crate::int_tools_wline::PntOn2S;

use super::prmprm_t3bits::{
    code_reject, dans_grille, grille_integer, integer_grille, nb_points_grille, remplit, T3Bits,
};

const XNBI: i32 = 200;
const IC15: i32 = 15;
const LIM: i32 = 3;

/// Sample both surfaces and emit approximate start points (`PointDepart`).
pub(crate) fn point_depart(
    s1: &dyn Surface,
    su1: i32,
    sv1: i32,
    s2: &dyn Surface,
    su2: i32,
    sv2: i32,
) -> Vec<PntOn2S> {
    let mut su1 = (IC15 * su1).min(XNBI).max(2);
    let mut sv1 = (IC15 * sv1).min(XNBI).max(2);
    let mut su2 = (IC15 * su2).min(XNBI).max(2);
    let mut sv2 = (IC15 * sv2).min(XNBI).max(2);

    let (u0_1, u1_1) = finite_uv(s1.u_range());
    let (v0_1, v1_1) = finite_uv(s1.v_range());
    let du1 = (u1_1 - u0_1) / (su1 - 1) as f64;
    let dv1 = (v1_1 - v0_1) / (sv1 - 1) as f64;
    let mut p1 = vec![GpPnt::zero(); (su1 * sv1) as usize];
    let mut box1 = BndBox::new();
    let mut dmax_on1: f64 = 0.0;
    for i in 0..su1 {
        let u = u0_1 + i as f64 * du1;
        for j in 0..sv1 {
            let v = v0_1 + j as f64 * dv1;
            let p = s1.d0(u, v);
            p1[idx(i, j, sv1)] = p;
            box1.add_point(&p);
            if i > 0 && j > 0 {
                let q = p1[idx(i - 1, j - 1, sv1)];
                let d = (q.x() - p.x()).abs() + (q.y() - p.y()).abs() + (q.z() - p.z()).abs();
                dmax_on1 = dmax_on1.max(d);
            }
        }
    }
    box1.enlarge(1.0e-8);

    let (u0_2, u1_2) = finite_uv(s2.u_range());
    let (v0_2, v1_2) = finite_uv(s2.v_range());
    let du2 = (u1_2 - u0_2) / (su2 - 1) as f64;
    let dv2 = (v1_2 - v0_2) / (sv2 - 1) as f64;
    let mut p2 = vec![GpPnt::zero(); (su2 * sv2) as usize];
    let mut box2 = BndBox::new();
    let mut dmax_on2: f64 = 0.0;
    for i in 0..su2 {
        let u = u0_2 + i as f64 * du2;
        for j in 0..sv2 {
            let v = v0_2 + j as f64 * dv2;
            let p = s2.d0(u, v);
            p2[idx(i, j, sv2)] = p;
            box2.add_point(&p);
            if i > 0 && j > 0 {
                let q = p2[idx(i - 1, j - 1, sv2)];
                let d = (q.x() - p.x()).abs() + (q.y() - p.y()).abs() + (q.z() - p.z()).abs();
                dmax_on2 = dmax_on2.max(d);
            }
        }
    }
    box2.enlarge(1.0e-8);
    if box1.is_out_box(&box2) {
        return Vec::new();
    }
    let Some((x10, y10, z10, x11, y11, z11)) = box1.get() else {
        return Vec::new();
    };
    let Some((x20, y20, z20, x21, y21, z21)) = box2.get() else {
        return Vec::new();
    };
    let mut x0 = x10.max(x20);
    let mut y0 = y10.max(y20);
    let mut z0 = z10.max(z20);
    let mut x1 = x11.min(x21);
    let mut y1 = y11.min(y21);
    let mut z1 = z11.min(z21);
    let mut dmax = dmax_on1.max(dmax_on2);
    dmax += dmax;
    x0 -= dmax;
    y0 -= dmax;
    z0 -= dmax;
    x1 += dmax;
    y1 += dmax;
    z1 += dmax;
    let mut bx10 = x10 - dmax;
    let mut by10 = y10 - dmax;
    let mut bz10 = z10 - dmax;
    let mut bx11 = x11 + dmax;
    let mut by11 = y11 + dmax;
    let mut bz11 = z11 + dmax;
    let mut bx20 = x20 - dmax;
    let mut by20 = y20 - dmax;
    let mut bz20 = z20 - dmax;
    let mut bx21 = x21 + dmax;
    let mut by21 = y21 + dmax;
    let mut bz21 = z21 + dmax;

    let nbg = nb_points_grille() as f64;
    let mut dx = (x1 - x0) / nbg;
    let mut dy = (y1 - y0) / nbg;
    let mut dz = (z1 - z0) / nbg;
    let dmx = dx.max(dy).max(dz);
    if dx < dmx * 0.01 {
        dx = dmx * 0.01;
    }
    if dy < dmx * 0.01 {
        dy = dmx * 0.01;
    }
    if dz < dmx * 0.01 {
        dz = dmx * 0.01;
    }
    let dx2 = dx * 0.5;
    let dy2 = dy * 0.5;
    let dz2 = dz * 0.5;

    let mut ip1 = vec![-1i32; p1.len()];
    let mut ds2 = vec![0i32; p1.len()];
    for i in 0..su1 {
        for j in 0..sv1 {
            let p = p1[idx(i, j, sv1)];
            ds2[idx(i, j, sv1)] = code_reject(bx20, by20, bz20, bx21, by21, bz21, p.x(), p.y(), p.z());
            let ix = ((p.x() - x0 + dx2) / dx) as i32;
            if dans_grille(ix) {
                let iy = ((p.y() - y0 + dy2) / dy) as i32;
                if dans_grille(iy) {
                    let iz = ((p.z() - z0 + dz2) / dz) as i32;
                    if dans_grille(iz) {
                        ip1[idx(i, j, sv1)] = grille_integer(ix, iy, iz);
                    }
                }
            }
        }
    }
    let mut ip2 = vec![-1i32; p2.len()];
    let mut ds1 = vec![0i32; p2.len()];
    for i in 0..su2 {
        for j in 0..sv2 {
            let p = p2[idx(i, j, sv2)];
            ds1[idx(i, j, sv2)] = code_reject(bx10, by10, bz10, bx11, by11, bz11, p.x(), p.y(), p.z());
            let ix = ((p.x() - x0 + dx2) / dx) as i32;
            if dans_grille(ix) {
                let iy = ((p.y() - y0 + dy2) / dy) as i32;
                if dans_grille(iy) {
                    let iz = ((p.z() - z0 + dz2) / dz) as i32;
                    if dans_grille(iz) {
                        ip2[idx(i, j, sv2)] = grille_integer(ix, iy, iz);
                    }
                }
            }
        }
    }

    let mut m1 = T3Bits::new(nb_points_grille());
    let mut m2 = T3Bits::new(nb_points_grille());
    for i in 0..su1 - 1 {
        for j in 0..sv1 - 1 {
            let a = ds2[idx(i, j, sv1)];
            if a & ds2[idx(i + 1, j, sv1)] == 0 || a & ds2[idx(i + 1, j + 1, sv1)] == 0 {
                remplit(
                    ip1[idx(i, j, sv1)],
                    ip1[idx(i + 1, j, sv1)],
                    ip1[idx(i + 1, j + 1, sv1)],
                    &mut m1,
                );
            }
            if a & ds2[idx(i, j + 1, sv1)] == 0 || a & ds2[idx(i + 1, j + 1, sv1)] == 0 {
                remplit(
                    ip1[idx(i, j, sv1)],
                    ip1[idx(i, j + 1, sv1)],
                    ip1[idx(i + 1, j + 1, sv1)],
                    &mut m1,
                );
            }
        }
    }
    for i in 0..su2 - 1 {
        for j in 0..sv2 - 1 {
            let a = ds1[idx(i, j, sv2)];
            if a & ds1[idx(i + 1, j, sv2)] == 0 || a & ds1[idx(i + 1, j + 1, sv2)] == 0 {
                remplit(
                    ip2[idx(i, j, sv2)],
                    ip2[idx(i + 1, j, sv2)],
                    ip2[idx(i + 1, j + 1, sv2)],
                    &mut m2,
                );
            }
            if a & ds1[idx(i, j + 1, sv2)] == 0 || a & ds1[idx(i + 1, j + 1, sv2)] == 0 {
                remplit(
                    ip2[idx(i, j, sv2)],
                    ip2[idx(i, j + 1, sv2)],
                    ip2[idx(i + 1, j + 1, sv2)],
                    &mut m2,
                );
            }
        }
    }

    let mut out = Vec::new();
    let mut newind = 0i32;
    let mut indice = 0i32;
    while m1.and_next(&mut m2, &mut newind) {
        indice -= 1;
        let (i, j, k) = integer_grille(newind);
        if dans_grille(i - 1)
            && dans_grille(j - 1)
            && dans_grille(k - 1)
            && dans_grille(i + 1)
            && dans_grille(j + 1)
            && dans_grille(k + 1)
        {
            let mut nb = 0;
            'lim: for si in -1..=1 {
                for sj in -1..=1 {
                    for sk in -1..=1 {
                        let lu = grille_integer(i + si, j + sj, k + sk);
                        if m1.val(lu) != 0 && m2.val(lu) != 0 {
                            nb += 1;
                            if nb >= LIM {
                                break 'lim;
                            }
                        }
                    }
                }
            }
            if nb >= LIM {
                for si in -1..=1 {
                    for sj in -1..=1 {
                        for sk in -1..=1 {
                            if si != 0 || sj != 0 || sk != 0 {
                                let lu = grille_integer(i + si, j + sj, k + sk);
                                m1.raz(lu);
                            }
                        }
                    }
                }
            }
        }
        let p = GpPnt::new(dx * i as f64 + x0, dy * j as f64 + y0, dz * k as f64 + z0);
        let mut nu1 = -1i32;
        let mut nv1 = 0i32;
        for nu in 0..su1 {
            for nv in 0..sv1 {
                if ip1[idx(nu, nv, sv1)] == newind {
                    ip1[idx(nu, nv, sv1)] = indice;
                    nu1 = nu;
                    nv1 = nv;
                    break;
                }
            }
            if nu1 >= 0 {
                break;
            }
        }
        let mut nu2 = -1i32;
        let mut nv2 = 0i32;
        if nu1 >= 0 {
            for nu in 0..su2 {
                for nv in 0..sv2 {
                    if ip2[idx(nu, nv, sv2)] == newind {
                        ip2[idx(nu, nv, sv2)] = indice;
                        nu2 = nu;
                        nv2 = nv;
                        break;
                    }
                }
                if nu2 >= 0 {
                    break;
                }
            }
        }
        if nu1 >= 0 && nu2 >= 0 {
            out.push(PntOn2S {
                p,
                u1: u0_1 + nu1 as f64 * du1,
                v1: v0_1 + nv1 as f64 * dv1,
                u2: u0_2 + nu2 as f64 * du2,
                v2: v0_2 + nv2 as f64 * dv2,
            });
        } else {
            let (u13, v13) = nearest3(&p1, su1, sv1, u0_1, v0_1, du1, dv1, &p);
            let (u23, v23) = nearest3(&p2, su2, sv2, u0_2, v0_2, du2, dv2, &p);
            out.push(PntOn2S {
                p,
                u1: u13,
                v1: v13,
                u2: u23,
                v2: v23,
            });
        }
    }
    let _ = (bx10, by10, bz10, bx11, by11, bz11, bx20, by20, bz20, bx21, by21, bz21);
    out
}

fn idx(i: i32, j: i32, sv: i32) -> usize {
    (i * sv + j) as usize
}

fn finite_uv(r: (f64, f64)) -> (f64, f64) {
    let (a, b) = r;
    match (a.is_finite(), b.is_finite()) {
        (true, true) => (a, b),
        (true, false) => (a, a + 100.0),
        (false, true) => (b - 100.0, b),
        (false, false) => (-100.0, 100.0),
    }
}

fn nearest3(
    pts: &[GpPnt],
    su: i32,
    sv: i32,
    u0: f64,
    v0: f64,
    du: f64,
    dv: f64,
    p: &GpPnt,
) -> (f64, f64) {
    let mut dist = [f64::MAX; 3];
    let mut uu = [0.0; 3];
    let mut vv = [0.0; 3];
    for i in 0..su {
        let u = u0 + i as f64 * du;
        for j in 0..sv {
            let v = v0 + j as f64 * dv;
            let t = pts[idx(i, j, sv)].square_distance(p);
            if dist[0] < dist[1] {
                dist.swap(0, 1);
                uu.swap(0, 1);
                vv.swap(0, 1);
            }
            if dist[1] < dist[2] {
                dist.swap(1, 2);
                uu.swap(1, 2);
                vv.swap(1, 2);
            }
            if dist[0] < dist[1] {
                dist.swap(0, 1);
                uu.swap(0, 1);
                vv.swap(0, 1);
            }
            if t < dist[0] {
                dist[0] = t;
                uu[0] = u;
                vv[0] = v;
            }
        }
    }
    ((uu[0] + uu[1] + uu[2]) / 3.0, (vv[0] + vv[1] + vv[2]) / 3.0)
}
