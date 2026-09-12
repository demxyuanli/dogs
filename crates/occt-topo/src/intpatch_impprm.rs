//! `IntPatch_ImpPrmIntersection` — implicit quadric × parametric surface.
//! Source: `IntPatch_ImpPrmIntersection.cxx` Perform: TheSOnBounds,
//! ComputeTangency, SearchInside, TheIWalking.

use occt_core::gp::{GpDir2d, GpVec};
use occt_geom::Surface;

use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::geom_int::{GeomIntLine, RestrictionArc, TopolTool};
use crate::int_tools_wline::{u_resolution, v_resolution, WLine};
use crate::intpatch::impimp::{
    distance, process_rline, process_segments, quadric_tolerance, search_on_bounds, try_set_quad,
    ImplicitQuad, PathPoint,
};

#[path = "intpatch_surf_func.rs"]
mod surf_func;
#[path = "intpatch_search_inside.rs"]
mod search_inside;
#[path = "intpatch_iwalking.rs"]
mod iwalking;

use iwalking::{iwalking_perform, lines_to_wlines, WalkStart};
use search_inside::search_inside;
use surf_func::SurfFunction;

#[path = "intpatch_impprm_ends.rs"]
mod ends;
use ends::attach_wline_ends;

/// Analytic types accepted as the implicit side (`ImpPrm` Quadric switch).
fn is_impprm_quadric(s: &dyn Surface) -> bool {
    matches!(
        classify_surface(s),
        SurfaceKind::Plane | SurfaceKind::Cylinder | SurfaceKind::Sphere | SurfaceKind::Cone
    )
}

/// Implicit × parametric intersection.
#[derive(Clone)]
pub struct ImpPrmIntersection {
    done: bool,
    empty: bool,
    slin: Vec<GeomIntLine>,
}

impl ImpPrmIntersection {
    pub fn new() -> Self {
        Self {
            done: false,
            empty: true,
            slin: Vec::new(),
        }
    }

    /// `Perform(S1, D1, S2, D2, TolArc, TolTang, Fleche, Pas)`.
    pub fn perform(
        &mut self,
        s1: &dyn Surface,
        d1: &TopolTool,
        s2: &dyn Surface,
        d2: &TopolTool,
        tol_arc: f64,
        tol_tang: f64,
        fleche: f64,
        pas: f64,
    ) {
        self.done = false;
        self.empty = true;
        self.slin.clear();

        let reversed = !is_impprm_quadric(s1);
        let (quad, prm, d_prm, quad_surf) = if !reversed {
            let Some(q) = try_set_quad(s1) else {
                return;
            };
            if !matches!(
                q,
                ImplicitQuad::Plane(_)
                    | ImplicitQuad::Cylinder(_)
                    | ImplicitQuad::Cone(_)
                    | ImplicitQuad::Sphere(_)
            ) {
                return;
            }
            (q, s2, d2, s1)
        } else {
            let Some(q) = try_set_quad(s2) else {
                return;
            };
            if !matches!(
                q,
                ImplicitQuad::Plane(_)
                    | ImplicitQuad::Cylinder(_)
                    | ImplicitQuad::Cone(_)
                    | ImplicitQuad::Sphere(_)
            ) {
                return;
            }
            (q, s1, d1, s2)
        };

        let a_local_pas = get_local_step(prm, pas);
        let mut func = SurfFunction::new(prm, &quad, quadric_tolerance(&quad));
        let solrst = search_on_bounds(prm, &quad, d_prm, tol_arc, tol_tang);
        if !solrst.done {
            return;
        }

        let (seqpdep, dest) = compute_tangency(&solrst.points, &mut func, prm);
        let mut search_ins = true;
        if matches!(quad, ImplicitQuad::Plane(_)) && !solrst.segments.is_empty() {
            search_ins = plane_needs_interior(&mut func, prm, d_prm);
        }

        let seqpins = if search_ins {
            search_inside(&mut func, prm, d_prm, tol_tang)
        } else {
            Vec::new()
        };

        // Walking when there are PathPoints or InteriorPoints
        // (`IntPatch_ImpPrmIntersection.cxx` ~847). Segments are handled below
        // even when walking is empty (`NbSegm` block ~1405).
        if !seqpdep.is_empty() || !seqpins.is_empty() {
            let raw = iwalking_perform(
                &seqpdep,
                &seqpins,
                &mut func,
                prm,
                reversed,
                a_local_pas,
                fleche,
            );
            let converted = lines_to_wlines(&raw, &quad, quad_surf, prm, reversed);
            let f_tol = quadric_tolerance(&quad).max(tol_tang);
            let kind1 = classify_surface(s1);
            let kind2 = classify_surface(s2);
            for (mut w, iw) in converted {
                if wline_on_quadric(&w, &quad, f_tol) {
                    attach_wline_ends(
                        &mut w,
                        &iw,
                        &seqpdep,
                        &solrst.points,
                        &dest,
                        reversed,
                        kind1,
                        kind2,
                        &quad,
                        prm,
                        quad_surf,
                        tol_arc,
                    );
                    self.slin.push(GeomIntLine::Walking(w));
                }
            }
        }

        // Whole-arc restriction solutions → RLine (`on_first` = parametric is S1).
        let on_first = reversed;
        if !solrst.segments.is_empty() {
            process_segments(&solrst.segments, &mut self.slin, s1, s2, on_first, tol_arc);
            process_rline(&mut self.slin, s1, s2, tol_arc, true);
        }

        self.empty = self.slin.is_empty();
        self.done = true;
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

impl Default for ImpPrmIntersection {
    fn default() -> Self {
        Self::new()
    }
}

/// Keep a walking line only when every 3D point lies on the implicit quadric
/// (`IntSurf_Quadric::Distance` vs `TolTang` / quadric tolerance).
fn wline_on_quadric(w: &WLine, quad: &ImplicitQuad, f_tol: f64) -> bool {
    let n = w.nb_pnts();
    if n < 2 {
        return false;
    }
    for i in 1..=n {
        if distance(quad, &w.point(i).p).abs() > f_tol {
            return false;
        }
    }
    true
}

/// `GetLocalStep` (BSpline/Bezier resolution branch).
fn get_local_step(surf: &dyn Surface, step: f64) -> f64 {
    let mut local = step;
    match classify_surface(surf) {
        SurfaceKind::Other => {
            let min_res = u_resolution(surf, occt_core::precision::CONFUSION)
                .min(v_resolution(surf, occt_core::precision::CONFUSION));
            if min_res < 1.0e-10 {
                local = 0.0001;
            }
        }
        _ => {}
    }
    local.min(step)
}

/// Plane + restriction segments: skip SearchInside unless F changes sign.
fn plane_needs_interior(func: &mut SurfFunction<'_>, prm: &dyn Surface, d: &TopolTool) -> bool {
    let n = d.nb_samples(prm);
    if n < 1 {
        return false;
    }
    let (s2d, _) = d.sample_point(prm, 1);
    let v0 = func.value(s2d.x(), s2d.y());
    let sign0 = v0.signum();
    for i in 2..=n {
        let (s2d, _) = d.sample_point(prm, i);
        let v = func.value(s2d.x(), s2d.y());
        if sign0 * v < 0.0 {
            return true;
        }
    }
    false
}

/// `ComputeTangency` for UV-box restriction zeros (`IsNew` branch).
/// Vertex/`TopTrans_CurveTransition` merge (`cxx` 329–465) is not-ported:
/// UV-box `TopolTool` has no `HVertex`.
fn compute_tangency(
    points: &[PathPoint],
    func: &mut SurfFunction<'_>,
    prm: &dyn Surface,
) -> (Vec<WalkStart>, Vec<i32>) {
    let mut out = Vec::new();
    let mut dest = vec![0i32; points.len()];
    for (i, pt) in points.iter().enumerate() {
        let _ = func.values(pt.u, pt.v);
        if func.is_tangent() {
            dest[i] = (out.len() + 1) as i32;
            out.push(WalkStart {
                p: pt.p,
                u: pt.u,
                v: pt.v,
                d3d: GpVec::zero(),
                d2d: GpDir2d::default(),
                tangent: true,
            });
            continue;
        }
        let mut vectg = func.direction3d();
        let mut dirtg = func.direction2d();
        let (_, d1u, d1v) = prm.d1(pt.u, pt.v);
        let (_p2d, d2d) = RestrictionArc::d1(&pt.arc, pt.param_on_arc);
        let v2 = d1u
            .multiplied_scalar(d2d.x())
            .added(&d1v.multiplied_scalar(d2d.y()));
        let v1 = d1u.crossed(&d1v);
        let test = vectg.dot(&v1.crossed(&v2));
        // `Adaptor3d_TopolTool::Orientation(Curve2d)` is always FORWARD.
        if test < 0.0 {
            vectg.reverse();
            dirtg.reverse();
        }
        dest[i] = (out.len() + 1) as i32;
        out.push(WalkStart {
            p: pt.p,
            u: pt.u,
            v: pt.v,
            d3d: vectg,
            d2d: dirtg,
            tangent: false,
        });
    }
    (out, dest)
}
