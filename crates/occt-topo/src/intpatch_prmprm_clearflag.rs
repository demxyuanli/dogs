//! `IntPatch_PrmPrmIntersection::Perform(Surf1, D1, Surf2, D2, TolTangency,
//! Epsilon, Deflection, Increment, ClearFlag)` (`cxx:2488-3402`).
//!
//! The overload first takes the sample parameters of both domains
//! (`Adaptor3d_TopolTool::SamplePnts`, `cxx:2495-2515`) and picks between two
//! strategies on the sample count (`Limit = 2500`, `cxx:2497`, `:2524`):
//!
//! * `NbU1*NbV1 <= Limit && NbU2*NbV2 <= Limit`: build an
//!   `IntPolyh_Intersection` on the same sampling (`:2533-2543`) and walk from
//!   its section lines (`:2575-2940`) and tangent zones (`:2950-3186`), then
//!   `AdjustOnPeriodic` (`:2944`). Split between
//!   `intpatch_prmprm_clearflag_lines.rs` and
//!   `intpatch_prmprm_clearflag_zones.rs`.
//! * otherwise: `PointDepart` (`:3200-3201`) over the `15 * NbSamples` lattice
//!   and walk from those approximate points (`:3239-3402`), in
//!   `intpatch_prmprm_clearflag_depart.rs`.

use occt_geom::Surface;

use crate::adaptor3d_topol_tool_sample::{is_uniform_sampling, sample_pnts};
use crate::geom_int::{GeomIntLine, TopolTool};
use crate::intpolyh::IntPolyhIntersection;

use super::clearflag_depart::walk_depart;
use super::clearflag_lines::walk_section_lines;
use super::clearflag_zones::walk_tangent_zones;
use super::periodic::adjust_on_periodic;

/// `Limit` (`cxx:2497`).
const LIMIT: i32 = 2500;

/// `NbU1 = NbV1 = NbU2 = NbV2 = 10` (`cxx:2498`, passed to `SamplePnts` at
/// `:2500-2501` as `theNUmin` / `theNVmin`).
const SAMPLE_MIN: i32 = 10;

/// The `Periods[4]` array of `cxx:2517-2522`.
pub(crate) fn periods_of(s1: &dyn Surface, s2: &dyn Surface) -> [f64; 4] {
    [
        if s1.is_u_periodic() { s1.u_period() } else { 0.0 },
        if s1.is_v_periodic() { s1.v_period() } else { 0.0 },
        if s2.is_u_periodic() { s2.u_period() } else { 0.0 },
        if s2.is_v_periodic() { s2.v_period() } else { 0.0 },
    ]
}

/// `Perform(Surf1, D1, Surf2, D2, TolTangency, Epsilon, Deflection, Increment,
/// ClearFlag)` (`cxx:2488-3402`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform_clear_flag(
    s1: &dyn Surface,
    d1: &TopolTool,
    s2: &dyn Surface,
    d2: &TopolTool,
    tol_tangency: f64,
    epsilon: f64,
    deflection: f64,
    increment: f64,
    clear_flag: bool,
    slin: &mut Vec<GeomIntLine>,
    done: &mut bool,
    empty: &mut bool,
) {
    // `D1->SamplePnts(Deflection, NbU1, NbV1)` / `D2->SamplePnts(...)`
    // (`cxx:2500-2501`).
    let sp1 = sample_pnts(s1, deflection, SAMPLE_MIN, SAMPLE_MIN);
    let sp2 = sample_pnts(s2, deflection, SAMPLE_MIN, SAMPLE_MIN);

    // `NbU1 = D1->NbSamplesU()` .. `NbV2 = D2->NbSamplesV()` (`cxx:2503-2506`).
    let nb_u1 = sp1.nb_u();
    let nb_v1 = sp1.nb_v();
    let nb_u2 = sp2.nb_u();
    let nb_v2 = sp2.nb_v();

    if nb_u1 * nb_v1 <= LIMIT && nb_v2 * nb_u2 <= LIMIT {
        *empty = true;
        if clear_flag {
            slin.clear();
        }

        // `new IntPolyh_Intersection(Surf1, NbU1, NbV1, Surf2, NbU2, NbV2)` at
        // `cxx:2537` when either domain samples uniformly, the parameter-array
        // constructor at `:2541` otherwise.
        let inter = if is_uniform_sampling(s1) || is_uniform_sampling(s2) {
            IntPolyhIntersection::with_sampling(s1, nb_u1, nb_v1, s2, nb_u2, nb_v2)
        } else {
            IntPolyhIntersection::with_parameters(
                s1,
                &sp1.u_pars,
                &sp1.v_pars,
                s2,
                &sp2.u_pars,
                &sp2.v_pars,
            )
        };

        // `done = Interference.IsDone()` (`cxx:2551`).
        *done = inter.is_done();
        if !*done {
            return;
        }

        let periods = periods_of(s1, s2);
        walk_section_lines(
            &inter,
            s1,
            d1,
            s2,
            d2,
            tol_tangency,
            epsilon,
            deflection,
            increment,
            &periods,
            slin,
            empty,
        );

        // `AdjustOnPeriodic(Surf1, Surf2, SLin)` (`cxx:2944`).
        adjust_on_periodic(s1, s2, slin);

        walk_tangent_zones(
            &inter,
            s1,
            d1,
            s2,
            d2,
            tol_tangency,
            epsilon,
            deflection,
            increment,
            slin,
            empty,
        );
        return;
    }

    // `PointDepart(LOn2S, Surf1, NbU1, NbV1, Surf2, NbU2, NbV2)`
    // (`cxx:3200-3201`).
    let points = super::prmprm_depart::point_depart(s1, nb_u1, nb_v1, s2, nb_u2, nb_v2);
    *empty = true;
    *done = true;

    walk_depart(
        &points,
        s1,
        d1,
        s2,
        d2,
        tol_tangency,
        epsilon,
        deflection,
        increment,
        slin,
        empty,
    );
}
