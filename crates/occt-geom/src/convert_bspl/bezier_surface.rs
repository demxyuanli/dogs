// ---------------------------------------------------------------------------
// BSpline surface → Bezier surface patches
// (GeomConvert_BSplineSurfaceToBezierSurface).
// ---------------------------------------------------------------------------

use occt_core::bspl::bezier::boehm_insert;
use occt_core::bspl::knots::{hunt, insert_knot, multiplicity};
use occt_core::gp::GpPnt;

use crate::bspline_surface::GeomBSplineSurface;

/// A single Bezier patch of a split B-spline surface, over the parameter
/// rectangle `[u0, u1] × [v0, v1]` with a `(deg_u+1) × (deg_v+1)` pole grid.
#[derive(Debug, Clone)]
pub struct BezierSurfacePatch {
    pub poles: Vec<Vec<GpPnt>>,
    pub weights: Option<Vec<Vec<f64>>>,
    pub u0: f64,
    pub u1: f64,
    pub v0: f64,
    pub v1: f64,
}

fn insert_surface_u(s: &GeomBSplineSurface, u: f64, mult: usize) -> GeomBSplineSurface {
    let nu = s.poles.len();
    let nv = s.poles[0].len();
    let deg = s.deg_u;
    let new_nu = nu + mult;
    let mut new_poles = vec![vec![GpPnt::zero(); nv]; new_nu];
    let mut new_weights = s.weights.as_ref().map(|_| vec![vec![0.0; nv]; new_nu]);
    for j in 0..nv {
        let mut col: Vec<GpPnt> = (0..nu).map(|i| s.poles[i][j]).collect();
        let mut wcol = s.weights.as_ref().map(|w| (0..nu).map(|i| w[i][j]).collect::<Vec<_>>());
        let idx = hunt(&s.knots_u, u);
        for _ in 0..mult {
            boehm_insert(&mut col, &s.knots_u, idx, u, deg, wcol.as_mut());
        }
        for i in 0..new_nu {
            new_poles[i][j] = col[i];
            if let (Some(nw), Some(wc)) = (new_weights.as_mut(), wcol.as_ref()) {
                nw[i][j] = wc[i];
            }
        }
    }
    let mut knots_u = s.knots_u.clone();
    for _ in 0..mult {
        knots_u = insert_knot(&knots_u, u, 1);
    }
    GeomBSplineSurface::from_parts(
        new_poles,
        new_weights,
        knots_u,
        s.knots_v().to_vec(),
        s.deg_u(),
        s.deg_v(),
        s.u_periodic(),
        s.v_periodic(),
    )
}

fn insert_surface_v(s: &GeomBSplineSurface, v: f64, mult: usize) -> GeomBSplineSurface {
    let nu = s.poles.len();
    let nv = s.poles[0].len();
    let deg = s.deg_v;
    let new_nv = nv + mult;
    let mut new_poles = vec![vec![GpPnt::zero(); new_nv]; nu];
    let mut new_weights = s.weights.as_ref().map(|_| vec![vec![0.0; new_nv]; nu]);
    for i in 0..nu {
        let mut row: Vec<GpPnt> = s.poles[i].clone();
        let mut wrow = s.weights.as_ref().map(|w| w[i].clone());
        let idx = hunt(&s.knots_v, v);
        for _ in 0..mult {
            boehm_insert(&mut row, &s.knots_v, idx, v, deg, wrow.as_mut());
        }
        for j in 0..new_nv {
            new_poles[i][j] = row[j];
            if let (Some(nw), Some(wr)) = (new_weights.as_mut(), wrow.as_ref()) {
                nw[i][j] = wr[j];
            }
        }
    }
    let mut knots_v = s.knots_v.clone();
    for _ in 0..mult {
        knots_v = insert_knot(&knots_v, v, 1);
    }
    GeomBSplineSurface::from_parts(
        new_poles,
        new_weights,
        s.knots_u().to_vec(),
        knots_v,
        s.deg_u(),
        s.deg_v(),
        s.u_periodic(),
        s.v_periodic(),
    )
}

fn distinct_in_range(knots: &[f64], first: f64, last: f64) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    for &k in knots {
        if k < first - 1e-15 || k > last + 1e-15 {
            continue;
        }
        if out.last().map_or(true, |&b| (b - k).abs() > 1e-15) {
            out.push(k);
        }
    }
    out
}

/// Convert a clamped B-spline surface into its Bezier patches by knot
/// insertion: every interior knot is raised to multiplicity `degree`, after
/// which each span is a Bezier patch. Mirrors
/// `GeomConvert_BSplineSurfaceToBezierSurface`.
pub fn bspline_surface_to_bezier_surface(s: &GeomBSplineSurface) -> Vec<BezierSurfacePatch> {
    let mut s = s.clone();
    let du = s.deg_u;
    let dv = s.deg_v;
    let (fu, lu) = (s.knots_u[du], s.knots_u[s.knots_u.len() - 1 - du]);
    let (fv, lv) = (s.knots_v[dv], s.knots_v[s.knots_v.len() - 1 - dv]);

    let interior_u: Vec<f64> = s
        .knots_u
        .iter()
        .copied()
        .filter(|&k| k > fu + 1e-15 && k < lu - 1e-15)
        .collect::<Vec<_>>()
        .into_iter()
        .fold(Vec::new(), |mut acc, k| {
            if acc.last().map_or(true, |&l| (l - k).abs() > 1e-15) {
                acc.push(k);
            }
            acc
        });
    for u in interior_u {
        let m = multiplicity(&s.knots_u, u);
        if m < du {
            s = insert_surface_u(&s, u, du - m);
        }
    }
    let interior_v: Vec<f64> = s
        .knots_v
        .iter()
        .copied()
        .filter(|&k| k > fv + 1e-15 && k < lv - 1e-15)
        .collect::<Vec<_>>()
        .into_iter()
        .fold(Vec::new(), |mut acc, k| {
            if acc.last().map_or(true, |&l| (l - k).abs() > 1e-15) {
                acc.push(k);
            }
            acc
        });
    for v in interior_v {
        let m = multiplicity(&s.knots_v, v);
        if m < dv {
            s = insert_surface_v(&s, v, dv - m);
        }
    }

    let fu = s.knots_u[du];
    let lu = s.knots_u[s.knots_u.len() - 1 - du];
    let fv = s.knots_v[dv];
    let lv = s.knots_v[s.knots_v.len() - 1 - dv];
    let bpu = distinct_in_range(&s.knots_u, fu, lu);
    let bpv = distinct_in_range(&s.knots_v, fv, lv);
    if bpu.len() < 2 || bpv.len() < 2 {
        return Vec::new();
    }

    let mut patches = Vec::with_capacity((bpu.len() - 1) * (bpv.len() - 1));
    for iu in 0..bpu.len() - 1 {
        for iv in 0..bpv.len() - 1 {
            let u_val = iu * du;
            let v_val = iv * dv;
            let mut poles = Vec::with_capacity(du + 1);
            let mut weights = s.weights.as_ref().map(|_| Vec::with_capacity(du + 1));
            for i in u_val..=u_val + du {
                poles.push(s.poles[i][v_val..=v_val + dv].to_vec());
                if let (Some(w), Some(sw)) = (weights.as_mut(), s.weights.as_ref()) {
                    w.push(sw[i][v_val..=v_val + dv].to_vec());
                }
            }
            patches.push(BezierSurfacePatch {
                poles,
                weights,
                u0: bpu[iu],
                u1: bpu[iu + 1],
                v0: bpv[iv],
                v1: bpv[iv + 1],
            });
        }
    }
    patches
}
