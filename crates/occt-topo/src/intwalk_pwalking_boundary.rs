//! `IntWalk_PWalking` boundary handling: snapping a walking line onto the
//! surface boundaries.
//!
//! Source: `IntWalk/IntWalk_PWalking.cxx`:
//! `HandleSingleSingularPoint` (`:2587-2710`),
//! `SeekPointOnBoundary` (`:2716-2947`),
//! `PutToBoundary` (`:2951-3156`),
//! `SeekAdditionalPoints` (`:3159-3294`).

use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, PCONFUSION, SQUARE_CONFUSION};
use occt_geom::Surface;

use crate::int_tools_wline::{u_resolution, v_resolution, PntOn2S};

use super::int2s::TheInt2S;
use super::iso::ConstIso;
use super::pwalking::PWalking;
use super::pwalking_dist::adjust_to_domain;

/// `IntImp_ConstIsoparametric` of the locked parameter of each of the four
/// UV coordinates, in u1, v1, u2, v2 order
/// (`IntWalk_PWalking.cxx:2602-2605`).
const LOCKED_DIR: [ConstIso; 4] = [
    ConstIso::UOnS1,
    ConstIso::VOnS1,
    ConstIso::UOnS2,
    ConstIso::VOnS2,
];

/// `IntWalk_PWalking.cxx:2786-2792`: midpoint of the two 3D values.
fn midpoint(p1: &GpPnt, p2: &GpPnt) -> GpPnt {
    GpPnt::new(
        0.5 * (p1.x() + p2.x()),
        0.5 * (p1.y() + p2.y()),
        0.5 * (p1.z() + p2.z()),
    )
}

impl<'a> PWalking<'a> {
    /// `HandleSingleSingularPoint` (`IntWalk_PWalking.cxx:2587-2710`).
    ///
    /// `pnt` is in u1, v1, u2, v2 order and is updated in place.
    pub(crate) fn handle_single_singular_point(&self, tol3d: f64, pnt: &mut [f64; 4]) -> bool {
        let (u1a, u1b) = self.s1.u_range();
        let (v1a, v1b) = self.s1.v_range();
        let (u2a, u2b) = self.s2.u_range();
        let (v2a, v2b) = self.s2.v_range();

        // u1, v1, u2, v2 order is used.
        let low = [u1a, v1a, u2a, v2a];
        let upp = [u1b, v1b, u2b, v2b];

        // Create new intersector with new tolerance.
        let mut int = TheInt2S::new(self.s1, self.s2, tol3d);

        for i in 0..4 {
            if (pnt[i] - low[i]).abs() < PCONFUSION || (pnt[i] - upp[i]).abs() < PCONFUSION {
                int.perform_iso(*pnt, LOCKED_DIR[i]);

                if !int.is_done() {
                    continue;
                }

                if int.is_empty() {
                    continue;
                }

                let p = int.point();
                let mut pars = [p.u1, p.v1, p.u2, p.v2];

                let surfs: [&dyn Surface; 2] = [self.s1, self.s2];
                // Local resolutions
                let tol2 = tol3d * tol3d;
                for k in 0..2 {
                    let iu = 2 * k;
                    let iv = iu + 1;
                    let (p_int, du, dv) = surfs[k].d1(pars[iu], pars[iv]);
                    let mut modu = du.magnitude();
                    if modu > CONFUSION {
                        let tol_u = tol3d / modu;
                        if (low[iu] - pars[iu]).abs() < tol_u {
                            let p = surfs[k].value(low[iu], pars[iv]);
                            if p_int.square_distance(&p) < tol2 {
                                pars[iu] = low[iu];
                            }
                        } else if (upp[iu] - pars[iu]).abs() < tol_u {
                            let p = surfs[k].value(upp[iu], pars[iv]);
                            if p_int.square_distance(&p) < tol2 {
                                pars[iu] = upp[iu];
                            }
                        }
                    }
                    modu = dv.magnitude();
                    if modu > CONFUSION {
                        let tol_v = tol3d / modu;
                        if (low[iv] - pars[iv]).abs() < tol_v {
                            let p = surfs[k].value(pars[iu], low[iv]);
                            if p_int.square_distance(&p) < tol2 {
                                pars[iv] = low[iv];
                            }
                        } else if (upp[iv] - pars[iv]).abs() < tol_v {
                            let p = surfs[k].value(pars[iu], upp[iv]);
                            if p_int.square_distance(&p) < tol2 {
                                pars[iv] = upp[iv];
                            }
                        }
                    }
                }

                for j in 0..4 {
                    pnt[j] = pars[j];
                }
                let mut is_in_domain = true;
                for j in 0..4 {
                    if (pnt[j] - low[j] + PCONFUSION) * (pnt[j] - upp[j] - PCONFUSION) > 0.0 {
                        is_in_domain = false;
                        break;
                    }
                }

                if is_in_domain {
                    return true;
                }
            }
        }

        false
    }

    /// `SeekPointOnBoundary` (`IntWalk_PWalking.cxx:2716-2947`).
    ///
    /// Minimizes the 3D distance between the two surfaces starting from the
    /// candidate parameters and, when the result is a true intersection, inserts
    /// it at the head (`is_the_first`) or the tail of the walking line, deleting
    /// the points whose insertion would create a hairpin bend.
    pub(crate) fn seek_point_on_boundary(
        &mut self,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        is_the_first: bool,
    ) -> bool {
        let mut is_ok = false;

        let (u1a, u1b) = self.s1.u_range();
        let (v1a, v1b) = self.s1.v_range();
        let (u2a, u2b) = self.s2.u_range();
        let (v2a, v2b) = self.s2.v_range();

        // u1, v1, u2, v2 order is used.
        let low = [u1a, v1a, u2a, v2a];
        let upp = [u1b, v1b, u2b, v2b];

        // Tune solution tolerance according with object size.
        let res1 = (PCONFUSION / u_resolution(self.s1, 1.0))
            .max(PCONFUSION / v_resolution(self.s1, 1.0));
        let res2 = (PCONFUSION / u_resolution(self.s2, 1.0))
            .max(PCONFUSION / v_resolution(self.s2, 1.0));
        let tol3d = res1.max(res2);
        let tol = CONFUSION.max(tol3d);

        // u1, v1, u2, v2 order is used.
        let mut pnt = [u1, v1, u2, v2];
        let mut singular_pnt = pnt;

        let mut nb_iter = 20;
        let mut status;
        loop {
            nb_iter -= 1;
            status = self.distance_minimize_by_gradient(&mut pnt, None);
            if status && !adjust_to_domain(&mut pnt[0..4], &low, &upp) {
                break;
            }

            let p2 = self.s2.value(pnt[2], pnt[3]);
            status =
                self.distance_minimize_by_extrema(self.s2, &p2, &mut pnt[0..2], None);
            if status && !adjust_to_domain(&mut pnt[0..2], &low[0..2], &upp[0..2])
            {
                break;
            }

            let p1 = self.s1.value(pnt[0], pnt[1]);
            status =
                self.distance_minimize_by_extrema(self.s1, &p1, &mut pnt[2..4], None);
            if status && !adjust_to_domain(&mut pnt[2..4], &low[2..4], &upp[2..4])
            {
                break;
            }

            if status || nb_iter <= 0 {
                break;
            }
        }

        // Handle singular points.
        let singular_status = self.handle_single_singular_point(tol, &mut singular_pnt);
        if singular_status {
            pnt = singular_pnt;
        }

        if !status && !singular_status {
            return is_ok;
        }

        let p1 = self.s1.value(pnt[0], pnt[1]);
        let p2 = self.s2.value(pnt[2], pnt[3]);
        let p_int = midpoint(&p1, &p2);

        let sq_dist = p_int.square_distance(&p1);
        if sq_dist > tol * tol {
            return is_ok;
        }

        // Found point is true intersection point
        let ip = PntOn2S {
            p: p_int,
            u1: pnt[0],
            v1: pnt[1],
            u2: pnt[2],
            v2: pnt[3],
        };

        // The main idea of checks below is to define if insertion of
        // addition point (on the boundary) does not lead to invalid
        // intersection curve (e.g. having a loop).
        //
        // Loops are detected with rotation angle of the Walking-line (WL).
        // If there is hairpin bend then insertion is forbidden.
        if is_the_first {
            while self.line.len() > 1 {
                let nb_pnts = self.line.len();

                let mut p_ind = 1usize;
                let mut a_p1 = self.line[0].value();
                while p_ind <= nb_pnts {
                    a_p1 = self.line[p_ind - 1].value();
                    if a_p1.square_distance(&p_int) > SQUARE_CONFUSION {
                        break;
                    }
                    p_ind += 1;
                }

                let mut a_p2 = a_p1;
                p_ind += 1;
                while p_ind <= nb_pnts {
                    a_p2 = self.line[p_ind - 1].value();
                    if a_p1.square_distance(&a_p2) > SQUARE_CONFUSION {
                        break;
                    }
                    p_ind += 1;
                }

                if p_ind > nb_pnts {
                    return is_ok;
                }

                let dir01 = a_p1.xyz().subtracted(p_int.xyz());
                let dir12 = a_p2.xyz().subtracted(a_p1.xyz());

                if dir01.dot(&dir12) > 0.0 {
                    break;
                }

                self.line.remove(0);
            }

            let a_p1 = self.line[0].value();
            if a_p1.square_distance(&p_int) <= SQUARE_CONFUSION {
                self.line.remove(0);
            }

            self.line.insert(0, ip);
            is_ok = true;
        } else {
            while self.line.len() > 1 {
                let nb_pnts = self.line.len();

                let mut a_p_curr = self.line[nb_pnts - 1].value();
                let mut p_ind = nb_pnts;
                while p_ind > 0 {
                    a_p_curr = self.line[p_ind - 1].value();
                    if a_p_curr.square_distance(&p_int) > SQUARE_CONFUSION {
                        break;
                    }
                    p_ind -= 1;
                }

                let mut a_p_prev = a_p_curr;
                p_ind -= 1;
                while p_ind > 0 {
                    a_p_prev = self.line[p_ind - 1].value();
                    if a_p_curr.square_distance(&a_p_prev) > SQUARE_CONFUSION {
                        break;
                    }
                    p_ind -= 1;
                }

                if p_ind < 1 {
                    return is_ok;
                }

                let dir_pc = a_p_curr.xyz().subtracted(a_p_prev.xyz());
                let dir_cn = p_int.xyz().subtracted(a_p_curr.xyz());

                if dir_pc.dot(&dir_cn) > 0.0 {
                    break;
                }

                self.line.pop();
            }

            let nb_pnts = self.line.len();
            let a_p1 = self.line[nb_pnts - 1].value();
            if a_p1.square_distance(&p_int) <= SQUARE_CONFUSION {
                self.line.pop();
            }
            self.line.push(ip);

            is_ok = true;
        }

        is_ok
    }

    /// `PutToBoundary` (`IntWalk_PWalking.cxx:2951-3156`): snap the first and
    /// the last walking-line points onto a boundary lying closer than the
    /// parallelism tolerance.
    pub(crate) fn put_to_boundary(&mut self) -> bool {
        const TOL_MIN: f64 = CONFUSION;

        let mut has_been_added = false;

        let (u1b_first, u1b_last) = self.s1.u_range();
        let (u2b_first, u2b_last) = self.s2.u_range();
        let (v1b_first, v1b_last) = self.s1.v_range();
        let (v2b_first, v2b_last) = self.s2.v_range();

        let mut tol: f64 = 1.0;
        tol = tol.min(u1b_last - u1b_first);
        tol = tol.min(u2b_last - u2b_first);
        tol = tol.min(v1b_last - v1b_first);
        tol = tol.min(v2b_last - v2b_first) * 1.0e-3;

        if tol <= 2.0 * TOL_MIN {
            return has_been_added;
        }

        let mut is_need_adding = false;

        let (is_u1_parallel, is_v1_parallel) = self.is_parallel(true, tol);
        let (is_u2_parallel, is_v2_parallel) = self.is_parallel(false, tol);

        let (mut u1, mut v1, mut u2, mut v2) = self.line[0].parameters();
        let mut delta;

        if !is_v1_parallel {
            delta = u1 - u1b_first;
            if TOL_MIN < delta && delta < tol {
                u1 = u1b_first;
                is_need_adding = true;
            } else {
                delta = u1b_last - u1;
                if TOL_MIN < delta && delta < tol {
                    u1 = u1b_last;
                    is_need_adding = true;
                }
            }
        }

        if !is_v2_parallel {
            delta = u2 - u2b_first;
            if TOL_MIN < delta && delta < tol {
                u2 = u2b_first;
                is_need_adding = true;
            } else {
                delta = u2b_last - u2;
                if TOL_MIN < delta && delta < tol {
                    u2 = u2b_last;
                    is_need_adding = true;
                }
            }
        }

        if !is_u1_parallel {
            delta = v1 - v1b_first;
            if TOL_MIN < delta && delta < tol {
                v1 = v1b_first;
                is_need_adding = true;
            } else {
                delta = v1b_last - v1;
                if TOL_MIN < delta && delta < tol {
                    v1 = v1b_last;
                    is_need_adding = true;
                }
            }
        }

        if !is_u2_parallel {
            delta = v2 - v2b_first;
            if TOL_MIN < delta && delta < tol {
                v2 = v2b_first;
                is_need_adding = true;
            } else {
                delta = v2b_last - v2;
                if TOL_MIN < delta && delta < tol {
                    v2 = v2b_last;
                    is_need_adding = true;
                }
            }
        }

        if is_need_adding {
            has_been_added = self.seek_point_on_boundary(u1, v1, u2, v2, true);
        }

        let nb_pnts = self.line.len();
        is_need_adding = false;
        let (mut u1, mut v1, mut u2, mut v2) = self.line[nb_pnts - 1].parameters();

        if !is_v1_parallel {
            delta = u1 - u1b_first;
            if TOL_MIN < delta && delta < tol {
                u1 = u1b_first;
                is_need_adding = true;
            } else {
                delta = u1b_last - u1;
                if TOL_MIN < delta && delta < tol {
                    u1 = u1b_last;
                    is_need_adding = true;
                }
            }
        }

        if !is_v2_parallel {
            delta = u2 - u2b_first;
            if TOL_MIN < delta && delta < tol {
                u2 = u2b_first;
                is_need_adding = true;
            } else {
                delta = u2b_last - u2;
                if TOL_MIN < delta && delta < tol {
                    u2 = u2b_last;
                    is_need_adding = true;
                }
            }
        }

        if !is_u1_parallel {
            delta = v1 - v1b_first;
            if TOL_MIN < delta && delta < tol {
                v1 = v1b_first;
                is_need_adding = true;
            } else {
                delta = v1b_last - v1;
                if TOL_MIN < delta && delta < tol {
                    v1 = v1b_last;
                    is_need_adding = true;
                }
            }
        }

        if !is_u2_parallel {
            delta = v2 - v2b_first;
            if TOL_MIN < delta && delta < tol {
                v2 = v2b_first;
                is_need_adding = true;
            } else {
                delta = v2b_last - v2;
                if TOL_MIN < delta && delta < tol {
                    v2 = v2b_last;
                    is_need_adding = true;
                }
            }
        }

        if is_need_adding {
            has_been_added = self.seek_point_on_boundary(u1, v1, u2, v2, false);
        }

        has_been_added
    }

    /// `SeekAdditionalPoints` (`IntWalk_PWalking.cxx:3159-3294`): insert
    /// intermediate points until the line has at least `the_min_nb_points`.
    pub(crate) fn seek_additional_points(&mut self, the_min_nb_points: i32) -> bool {
        const TOL: f64 = 1.0e-14;
        let mut nb_points = self.line.len() as i32;
        if nb_points > the_min_nb_points {
            return true;
        }

        let (u1b_first, u1b_last) = self.s1.u_range();
        let (u2b_first, u2b_last) = self.s2.u_range();
        let (v1b_first, v1b_last) = self.s1.v_range();
        let (v2b_first, v2b_last) = self.s2.v_range();

        let mut is_precise = false;

        let mut pnt = [0.0f64; 4];

        let mut nb_points_prev = 0;
        while nb_points < the_min_nb_points && nb_points != nb_points_prev {
            nb_points_prev = nb_points;
            let mut fp = 1i32;
            let mut lp;
            while fp < nb_points {
                lp = fp + 1;
                let (u1f, v1f, u2f, v2f) = self.line[(fp - 1) as usize].parameters();
                let (u1l, v1l, u2l, v2l) = self.line[(lp - 1) as usize].parameters();

                pnt[0] = 0.5 * (u1f + u1l);
                if pnt[0] < u1b_first {
                    pnt[0] = u1b_first;
                }
                if pnt[0] > u1b_last {
                    pnt[0] = u1b_last;
                }

                pnt[1] = 0.5 * (v1f + v1l);
                if pnt[1] < v1b_first {
                    pnt[1] = v1b_first;
                }
                if pnt[1] > v1b_last {
                    pnt[1] = v1b_last;
                }

                pnt[2] = 0.5 * (u2f + u2l);
                if pnt[2] < u2b_first {
                    pnt[2] = u2b_first;
                }
                if pnt[2] > u2b_last {
                    pnt[2] = u2b_last;
                }

                pnt[3] = 0.5 * (v2f + v2l);
                if pnt[3] < v2b_first {
                    pnt[3] = v2b_first;
                }
                if pnt[3] > v2b_last {
                    pnt[3] = v2b_last;
                }

                let mut status;
                let mut nb_iter = 5;
                loop {
                    status = self.distance_minimize_by_gradient(&mut pnt, None);
                    if status {
                        break;
                    }

                    let p2 = self.s2.value(pnt[2], pnt[3]);
                    status =
                        self.distance_minimize_by_extrema(self.s2, &p2, &mut pnt[0..2], None);
                    if status {
                        break;
                    }

                    let p1 = self.s1.value(pnt[0], pnt[1]);
                    status =
                        self.distance_minimize_by_extrema(self.s1, &p1, &mut pnt[2..4], None);
                    if status {
                        break;
                    }

                    nb_iter -= 1;
                    if nb_iter <= 0 {
                        break;
                    }
                }

                if status {
                    let a_p1 = self.s1.value(pnt[0], pnt[1]);
                    let a_p2 = self.s2.value(pnt[2], pnt[3]);
                    let p_int = midpoint(&a_p1, &a_p2);

                    let sq_dist1 = p_int.square_distance(&a_p1);
                    let sq_dist2 = p_int.square_distance(&a_p2);

                    if sq_dist1 < TOL && sq_dist2 < TOL {
                        let ip = PntOn2S {
                            p: p_int,
                            u1: pnt[0],
                            v1: pnt[1],
                            u2: pnt[2],
                            v2: pnt[3],
                        };
                        self.line.insert((lp - 1) as usize, ip);

                        is_precise = true;

                        nb_points += 1;
                        if nb_points >= the_min_nb_points {
                            break;
                        }
                    } else {
                        lp -= 1;
                    }
                }

                fp = lp + 1;
            }
        }

        is_precise
    }
}
