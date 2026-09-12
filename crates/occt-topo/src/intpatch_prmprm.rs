//! `IntPatch_PrmPrmIntersection` — two bi-parametric surfaces.
//! Source: `IntPatch_PrmPrmIntersection.cxx` Perform without a start point:
//! `PointDepart` then `IntWalk_PWalking`.

use occt_core::gp::{GpPnt, GpVec};
use occt_geom::Surface;

use crate::geom_int::{GeomIntLine, TopolTool};
use crate::int_tools_wline::{PntOn2S, WLine, WLineWay};
use crate::intpatch::{intersect_general_surfaces, SurfaceIntersection};
use crate::intwalk::PWalking;

#[path = "intpatch_prmprm_t3bits.rs"]
mod prmprm_t3bits;
#[path = "intpatch_prmprm_depart.rs"]
mod prmprm_depart;

use prmprm_depart::point_depart;

/// Parametric × parametric intersection.
#[derive(Clone)]
pub struct PrmPrmIntersection {
    done: bool,
    empty: bool,
    slin: Vec<GeomIntLine>,
}

impl PrmPrmIntersection {
    pub fn new() -> Self {
        Self {
            done: true,
            empty: true,
            slin: Vec::new(),
        }
    }

    /// `Perform(S1, D1, S2, D2, TolTangency, Epsilon, Deflection, Increment)`.
    pub fn perform(
        &mut self,
        s1: &dyn Surface,
        d1: &TopolTool,
        s2: &dyn Surface,
        d2: &TopolTool,
        tol_tangency: f64,
        epsilon: f64,
        deflection: f64,
        increment: f64,
    ) {
        self.done = true;
        self.empty = true;
        self.slin.clear();

        let nbu1 = d1.nb_samples_u(s1).max(2);
        let nbv1 = d1.nb_samples_v(s1).max(2);
        let nbu2 = d2.nb_samples_u(s2).max(2);
        let nbv2 = d2.nb_samples_v(s2).max(2);
        let starts = point_depart(s1, nbu1, nbv1, s2, nbu2, nbv2);
        let n_inc = increment;
        let mut pw = PWalking::new(s1, s2, tol_tangency, epsilon, deflection, n_inc);
        let seuil = 15.0 * increment * increment;

        for st in &starts {
            let start = [st.u1, st.v1, st.u2, st.v2];
            let Some(first) = pw.perform_first_point(start) else {
                continue;
            };
            if on_existing_line(&first, &self.slin, deflection) {
                continue;
            }
            pw.perform(start);
            if !pw.is_done() || pw.nb_points() <= 2 {
                continue;
            }
            let last = pw.value(pw.nb_points());
            let debut = pw.value(1).p;
            let fin = last.p;
            let mut rejet = false;
            for line in &self.slin {
                if let Some(wl) = line.as_wline() {
                    if is_point_on_line(&last, wl, deflection) {
                        rejet = true;
                        break;
                    }
                    if wl.nb_pnts() >= 1 {
                        let a = wl.point(1).p;
                        let b = wl.point(wl.nb_pnts()).p;
                        if debut.distance(&a) < tol_tangency || fin.distance(&b) < tol_tangency {
                            rejet = true;
                            break;
                        }
                    }
                }
            }
            if rejet {
                continue;
            }
            let mut wl = WLine::new();
            wl.set_creating_way(WLineWay::PrmPrm);
            for p in pw.line() {
                wl.add(*p);
            }
            wl.ensure_end_vertices();
            if let Some(purged) = crate::intpatch::wline_tool::compute_purged_wline(
                &wl, s1, s2, d1, d2,
            ) {
                if purged.nb_pnts() >= 2 {
                    self.slin.push(GeomIntLine::Walking(purged));
                }
            } else if wl.nb_pnts() >= 2 {
                self.slin.push(GeomIntLine::Walking(wl));
            }
        }

        if self.slin.is_empty() {
            match intersect_general_surfaces(s1, s2, tol_tangency) {
                SurfaceIntersection::Curves(ics) => {
                    for ic in &ics {
                        self.slin
                            .push(GeomIntLine::Walking(WLine::from_intersection_curve_way(
                                ic,
                                WLineWay::PrmPrm,
                            )));
                    }
                }
                SurfaceIntersection::Coincident | SurfaceIntersection::None => {}
            }
        }
        self.empty = self.slin.is_empty();
        let _ = seuil;
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn is_empty(&self) -> bool {
        self.empty
    }

    pub fn lines(&self) -> &[GeomIntLine] {
        &self.slin
    }
}

impl Default for PrmPrmIntersection {
    fn default() -> Self {
        Self::new()
    }
}

fn on_existing_line(p: &PntOn2S, lines: &[GeomIntLine], deflection: f64) -> bool {
    lines.iter().any(|l| {
        l.as_wline()
            .map(|wl| is_point_on_line(p, wl, deflection))
            .unwrap_or(false)
    })
}

/// `IsPointOnLine` (segment projection within `Deflection`).
fn is_point_on_line(p: &PntOn2S, wl: &WLine, deflection: f64) -> bool {
    let n = wl.nb_pnts();
    if n < 2 {
        return false;
    }
    let def2 = deflection * deflection;
    let start = p.p;
    for ll in 1..n {
        let pa = wl.point(ll).p;
        let pb = wl.point(ll + 1).p;
        if dist_point_segment_sq(&start, &pa, &pb) <= def2 {
            return true;
        }
    }
    false
}

fn dist_point_segment_sq(p: &GpPnt, a: &GpPnt, b: &GpPnt) -> f64 {
    let ab = GpVec::from_pnts(a, b);
    let ap = GpVec::from_pnts(a, p);
    let l2 = ab.square_magnitude();
    if l2 <= f64::EPSILON {
        return p.square_distance(a);
    }
    let t = (ap.dot(&ab) / l2).clamp(0.0, 1.0);
    let q = GpPnt::new(
        a.x() + t * (b.x() - a.x()),
        a.y() + t * (b.y() - a.y()),
        a.z() + t * (b.z() - a.z()),
    );
    p.square_distance(&q)
}
