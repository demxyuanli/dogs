//! Line bookkeeping of the walking engine: `MakeWalkingPoint`
//! (`IntWalk_IWalking.gxx:2904-2940`), `AddPointInCurrentLine`
//! (`gxx:2894-2903`), `OpenLine` (`gxx:2942-2987`), `IsValidEndPoint`
//! (`gxx:2989-3001`), `RemoveTwoEndPoints` (`gxx:3003-3016`) and the two
//! `IsPointOnLine` overloads (`gxx:3018-3140`).

use occt_core::gp::{GpPnt, GpPnt2d, GpVec};
use occt_core::precision::REAL_EPSILON;
use occt_math::function_set_root::MathFunctionSetWithDerivatives;
use occt_math::{MathFunctionSetRoot, MathMatrix, MathVector};

use occt_topo::int_tools_wline::PntOn2S;

use super::types::{math_vec2, pnt2d_on, set_uv, uv_on, PathPoint};
use super::{ContapIWLineOfTheIWalking, ContapTheIWalking};

/// `MakeWalkingPoint(Case, U, V, sp, Psol)` (`gxx:2904-2940`).
///
/// `reversed` is the walking flag: `SetValue(P, reversed, U, V)` writes the
/// parameters on surface 1 when it is set.
pub(super) fn make_walking_point(
    reversed: bool,
    case: i32,
    u: f64,
    v: f64,
    sp: &mut super::super::surf_function::ContapSurfFunction,
    psol: &mut PntOn2S,
) {
    if case == 1 || case == 2 {
        psol.set_point(*sp.point());
        set_uv(psol, !reversed, u, v);
    } else if case == 11 || case == 12 {
        let mut ff = MathVector::new(1, 1);
        let mut dd = MathMatrix::new(1, 1, 1, 2);
        let _ = sp.values(&math_vec2(u, v), &mut ff, &mut dd);
        make_walking_point(reversed, case - 10, u, v, sp, psol);
    } else {
        panic!("MakeWalkingPoint: Standard_ConstructionError");
    }
}

impl ContapTheIWalking {
    /// `AddPointInCurrentLine(N, PathPnt, CurrentLine)` (`gxx:2894-2903`).
    pub(super) fn add_point_in_current_line(
        &self,
        n: i32,
        path_pnt: &PathPoint,
        line: &mut ContapIWLineOfTheIWalking,
    ) {
        // `IntSurf_PntOn2S Psol;` - a fresh point, so only the parameters of
        // the current surface are written.
        let mut psol = PntOn2S {
            p: *path_pnt.value(),
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
        };
        set_uv(
            &mut psol,
            !self.reversed,
            self.wd1[n as usize].ustart,
            self.wd1[n as usize].vstart,
        );
        line.add_point(psol);
    }

    /// `OpenLine(N, Psol, Pnts1, sp, Line)` (`gxx:2942-2987`): reverse the
    /// line and restart in the other direction.
    pub(super) fn open_line(
        &mut self,
        n: i32,
        psol: &PntOn2S,
        pnts1: &[PathPoint],
        sp: &mut super::super::surf_function::ContapSurfFunction,
        line: &mut ContapIWLineOfTheIWalking,
    ) {
        self.previous_point = *line.value(1);
        let (u, v) = uv_on(&self.previous_point, self.reversed);
        let mut ff = MathVector::new(1, 1);
        let mut dd = MathMatrix::new(1, 1, 1, 2);
        let _ = sp.values(&math_vec2(u, v), &mut ff, &mut dd);
        if let Ok(d3d) = sp.direction3d() {
            self.previous_d3d = *d3d;
        }
        if let Ok(d2d) = sp.direction2d() {
            self.previous_d2d = d2d;
        }

        if n > 0 {
            let path_pnt = pnts1[(n - 1) as usize].clone();
            line.add_status_first_with_point(false, true, n, path_pnt.clone());
            self.add_point_in_current_line(n, &path_pnt, line);
        } else {
            if n < 0 {
                line.add_point(*psol);
            }
            line.add_status_first(false, false);
        }
        line.reverse();
        let nb = line.nb_points();
        line.set_tangent_vector(self.previous_d3d.reversed(), nb);
    }

    /// `IsValidEndPoint(IndOfPoint, IndOfLine)` (`gxx:2989-3001`).
    pub(super) fn is_valid_end_point(&mut self, ind_of_point: i32, ind_of_line: i32) -> bool {
        if self.point_line_line.is_empty() {
            return true;
        }
        let Some(list) = self.point_line_line.get_mut(&ind_of_point) else {
            return false;
        };
        if let Some(pos) = list.iter().position(|&v| v == ind_of_line) {
            list.remove(pos);
            return true;
        }
        false
    }

    /// `RemoveTwoEndPoints(IndOfPoint)` (`gxx:3003-3016`).
    pub(super) fn remove_two_end_points(&mut self, ind_of_point: i32) {
        if let Some(list) = self.point_line_line.get(&ind_of_point) {
            if list.is_empty() {
                return;
            }
            let line1 = list[0];
            let line2 = *list.last().unwrap_or(&line1);
            self.seq_alone
                .retain(|&v| v != line1 && v != line2);
        }
    }

    /// `IsPointOnLine(theP2d, Irang)` (`gxx:3018-3037`).
    pub(super) fn is_point_on_line_2d(&self, the_p2d: &GpPnt2d, irang: i32) -> bool {
        let Some(a_line) = self
            .lines
            .get((irang.unsigned_abs() as usize).saturating_sub(1))
        else {
            return false;
        };
        let nbp = a_line.nb_points();
        for i in 1..=nbp {
            let p2d1 = pnt2d_on(a_line.value(i), self.reversed);
            if (p2d1.coord.x - the_p2d.coord.x).abs() <= self.vec_tolerance[0]
                && (p2d1.coord.y - the_p2d.coord.y).abs() <= self.vec_tolerance[1]
            {
                return true;
            }
            if i < nbp {
                let p2d2 = pnt2d_on(a_line.value(i + 1), self.reversed);
                let pp1 = occt_core::gp::GpVec2d::new(
                    p2d1.coord.x - the_p2d.coord.x,
                    p2d1.coord.y - the_p2d.coord.y,
                );
                let pp2 = occt_core::gp::GpVec2d::new(
                    p2d2.coord.x - the_p2d.coord.x,
                    p2d2.coord.y - the_p2d.coord.y,
                );
                if pp1.dot(&pp2) < 0.0 {
                    return true;
                }
            }
        }
        false
    }

    /// `IsPointOnLine(thePOn2S, InfBounds, SupBounds, Solver, Func)`
    /// (`gxx:3041-3140`): project the point on the nearest segment of every
    /// already computed line and refine the projection with `Solver`.
    pub(super) fn is_point_on_line(
        &self,
        the_p_on2s: &PntOn2S,
        inf_bounds: &MathVector,
        sup_bounds: &MathVector,
        solver: &mut MathFunctionSetRoot,
        the_func: &mut super::super::surf_function::ContapSurfFunction,
    ) -> bool {
        let a_p3d: GpPnt = the_p_on2s.p;

        for a_line in self.lines.iter() {
            if is_out_box(a_line, &a_p3d) {
                continue;
            }

            let mut a_u_min = 0.0;
            let mut a_v_min = 0.0;
            let mut a_min_sq_dist = f64::MAX;
            let nbp = a_line.nb_points();
            let mut pt_idx = 1;
            while pt_idx < nbp {
                let a_p1 = a_line.value(pt_idx).p;
                let a_p2 = a_line.value(pt_idx + 1).p;

                let a_p1p = GpVec::from_pnts(&a_p1, &a_p3d);
                let a_p1p2 = GpVec::from_pnts(&a_p1, &a_p2);

                let a_sq12 = a_p1p2.square_magnitude();
                if a_sq12 < occt_core::precision::RESOLUTION {
                    pt_idx += 1;
                    continue;
                }

                let a_dp = a_p1p.dot(&a_p1p2);
                if a_dp < 0.0 || a_dp > a_sq12 {
                    pt_idx += 1;
                    continue;
                }

                let a_sq_d = a_p1p.cross_square_magnitude(&a_p1p2) / a_sq12;
                if a_sq_d < a_min_sq_dist {
                    a_min_sq_dist = a_sq_d;

                    let a_l1 = a_dp / a_sq12;
                    let a_l2 = 1.0 - a_l1;
                    if a_l1 < REAL_EPSILON || a_l2 < REAL_EPSILON {
                        return true;
                    }

                    let (a_u1, a_v1) = uv_on(a_line.value(pt_idx), self.reversed);
                    let (a_u2, a_v2) = uv_on(a_line.value(pt_idx + 1), self.reversed);
                    a_u_min = a_l1 * a_u2 + a_l2 * a_u1;
                    a_v_min = a_l1 * a_v2 + a_l2 * a_v1;
                }
                pt_idx += 1;
            }

            if a_min_sq_dist > occt_core::precision::INFINITE {
                continue;
            }

            let mut a_vec_prms = math_vec2(a_u_min, a_v_min);
            solver.perform_with_bounds(the_func, &a_vec_prms, inf_bounds, sup_bounds, false);
            if !solver.is_done() {
                continue;
            }

            a_vec_prms = solver.root();

            let surface = match the_func.surface() {
                Some(s) => s,
                None => continue,
            };
            let a_pa = match surface.value(a_u_min, a_v_min) {
                Some(p) => p,
                None => continue,
            };
            let a_pb = match surface.value(a_vec_prms.value(1), a_vec_prms.value(2)) {
                Some(p) => p,
                None => continue,
            };
            let a_sq_d1 = a_pb.square_distance(&a_p3d);
            let a_sq_d2 = a_pa.square_distance(&a_pb);

            if a_sq_d1 < 4.0 * a_sq_d2 {
                return true;
            }
        }

        false
    }
}

/// `IntSurf_LineOn2S::IsOutBox(Pxyz)` (`IntSurf_LineOn2S.cxx`): the line's 3D
/// box enlarged by 1% of its largest extent.
fn is_out_box(line: &ContapIWLineOfTheIWalking, p: &GpPnt) -> bool {
    if line.nb_points() == 0 {
        return false;
    }
    let mut min = [f64::MAX; 3];
    let mut max = [-f64::MAX; 3];
    for i in 1..=line.nb_points() {
        let q = line.value(i).p;
        let c = [q.x(), q.y(), q.z()];
        for k in 0..3 {
            min[k] = min[k].min(c[k]);
            max[k] = max[k].max(c[k]);
        }
    }
    let dx = max[0] - min[0];
    let dy = max[1] - min[1];
    let dz = max[2] - min[2];
    let largest = if dx > dy { dx.max(dz) } else { dy.max(dz) };
    let gap = largest * 0.01;

    let c = [p.x(), p.y(), p.z()];
    for k in 0..3 {
        if c[k] < min[k] - gap || c[k] > max[k] + gap {
            return true;
        }
    }
    false
}
