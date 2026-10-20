use occt_core::gp::GpPnt;
use occt_core::precision::epsilon;
use occt_geom::extrema_pc::extrema_ext_pc_range;
use occt_geom::Curve;

// ---------------------------------------------------------------------------
// Projection onto a restricted range (`GeomAPI_ProjectPointOnCurve`)
// ---------------------------------------------------------------------------

/// `GeomAPI_ProjectPointOnCurve` restricted to `[t1, t2]`
/// (`GeomAPI_ProjectPointOnCurve.cxx:135-156`): `NbPoints() == 0` is the
/// "cannot project" answer, otherwise the *smallest* squared distance over all
/// extrema (no `myIsMin` filter) and its parameter.
pub(in crate::edge_edge) struct ProjPc<'a> {
    curve: &'a dyn Curve,
    t1: f64,
    t2: f64,
}

impl<'a> ProjPc<'a> {
    pub(in crate::edge_edge) fn new(curve: &'a dyn Curve, t1: f64, t2: f64) -> Self {
        Self { curve, t1, t2 }
    }

    /// `Perform` + `NbPoints()`: `None` when nothing projects.
    pub(in crate::edge_edge) fn project(&self, p: &GpPnt) -> Option<(f64, f64)> {
        let sols = extrema_ext_pc_range(self.curve, p, self.t1, self.t2);
        let best = sols.iter().reduce(|a, b| if b.sq_dist < a.sq_dist { b } else { a })?;
        Some((best.sq_dist.sqrt(), best.u))
    }
}

/// Running extremum search state shared by `DistPC` / `FindDistPC`
/// (`IntTools_EdgeEdge.cxx:1210-1362`); `d_max` is a *maximum* for `i_c = 1`
/// and a *minimum* for `i_c = -1`, exactly as in OCCT.
#[derive(Debug, Clone, Copy)]
pub(in crate::edge_edge) struct DistState {
    pub(in crate::edge_edge) d_max: f64,
    pub(in crate::edge_edge) t1_max: f64,
    pub(in crate::edge_edge) t2_max: f64,
}

/// The 7-argument `DistPC` (`IntTools_EdgeEdge.cxx:1332-1362`): the error code
/// plus the distance and the second parameter.
pub(in crate::edge_edge) fn dist_pc(
    proj: &ProjPc<'_>,
    c1: &dyn Curve,
    t1: f64,
    criteria: f64,
    i_c: i32,
) -> (i32, f64, f64) {
    let p1 = c1.d0(t1);
    let Some((d, t2)) = proj.project(&p1) else {
        // `iErr = 1`: the point of C1 cannot be projected on C2.
        return (1, 0.0, 0.0);
    };
    let i_err = if (i_c as f64) * (d - criteria) > 0.0 {
        2 // the distance is too big or too small
    } else {
        0
    };
    (i_err, d, t2)
}

/// The 9-argument `DistPC` (`IntTools_EdgeEdge.cxx:1301-1328`): as above plus
/// the running extremum update. Returns `(iErr, aD)`.
fn dist_pc_max(
    proj: &ProjPc<'_>,
    c1: &dyn Curve,
    t1: f64,
    criteria: f64,
    i_c: i32,
    st: &mut DistState,
) -> (i32, f64) {
    let (i_err, d, t2) = dist_pc(proj, c1, t1, criteria, i_c);
    if i_err == 1 {
        return (i_err, d);
    }
    if (i_c as f64) * (d - st.d_max) > 0.0 {
        st.d_max = d;
        st.t1_max = t1;
        st.t2_max = t2;
    }
    (i_err, d)
}

/// `FindDistPC` (`IntTools_EdgeEdge.cxx:1210-1297`) — golden-section search of
/// the extremal distance of C1 on `[t1a, t1b]` against the projector.
pub(in crate::edge_edge) fn find_dist_pc(
    proj: &ProjPc<'_>,
    c1: &dyn Curve,
    t1a: f64,
    t1b: f64,
    criteria: f64,
    eps: f64,
    st: &mut DistState,
    b_max_dist: bool,
) -> i32 {
    let i_c = if b_max_dist { 1 } else { -1 };
    let gs = 0.618_033_988_749_894_8; // = 0.5*(1+sqrt(5))-1
    let mut a = t1a;
    let mut b = t1b;
    st.t1_max = 0.0;
    st.t2_max = 0.0;

    // check bounds
    let (i_err, _) = dist_pc_max(proj, c1, a, criteria, i_c, st);
    if i_err == 2 {
        return i_err;
    }
    let (i_err, _) = dist_pc_max(proj, c1, b, criteria, i_c, st);
    if i_err == 2 {
        return i_err;
    }
    let mut xp = a + (b - a) * gs;
    let mut xl = b - (b - a) * gs;
    let (i_err, mut yp) = dist_pc_max(proj, c1, xp, criteria, i_c, st);
    if i_err != 0 {
        return i_err;
    }
    let (i_err, mut yl) = dist_pc_max(proj, c1, xl, criteria, i_c, st);
    if i_err != 0 {
        return i_err;
    }

    let an_eps = eps.max(epsilon(a.abs().max(b.abs())) * 10.0);
    loop {
        let i_err;
        if (i_c as f64) * (yp - yl) > 0.0 {
            a = xl;
            xl = xp;
            yl = yp;
            xp = a + (b - a) * gs;
            let (e, d) = dist_pc_max(proj, c1, xp, criteria, i_c, st);
            i_err = e;
            yp = d;
        } else {
            b = xp;
            xp = xl;
            yp = yl;
            xl = b - (b - a) * gs;
            let (e, d) = dist_pc_max(proj, c1, xl, criteria, i_c, st);
            i_err = e;
            yl = d;
        }
        if i_err != 0 {
            if i_err == 2 && !b_max_dist {
                let xp_mid = (a + b) * 0.5;
                let _ = dist_pc_max(proj, c1, xp_mid, criteria, i_c, st);
            }
            return i_err;
        }
        if (b - a) < an_eps {
            break;
        }
    }
    0
}
