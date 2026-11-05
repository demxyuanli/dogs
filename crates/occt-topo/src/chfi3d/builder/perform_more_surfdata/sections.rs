//! `PerformMoreSurfdata` — the two surface sections and the fillet arc.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C1.cxx:3918-4156`
//! (OCCT 8.0.0).
//!
//! The block runs `ChFi3d_ComputeArete` over the two end common points of
//! `SeqFil(anInd)`, intersects the end face's surface with the fillet surface
//! of `anInd` (`cxx:3986`), then — when the far end of that section does not
//! land on either arc point — trims the fillet curve to the projection of that
//! far end and intersects a second time against the surface of `anIndPrev`
//! (`cxx:4067-4156`).
//!
//! Two details of the source control flow are preserved because they decide
//! which line survives the loops:
//!
//! - both `IntSS` loops stop on the *first* line whose two ends come within
//!   `aTol3d` of the arc point, and the line handles / `aTolex` keep that last
//!   examined line (`cxx:4005-4025`, `:4130-4148`),
//! - the `Extrema_ExtPC` loop has `&& !isFound` in its condition but never sets
//!   `isFound`, so every minimum is scanned and `aPar` keeps the *last*
//!   qualifying one (`cxx:4061-4076`).

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom::curve::Curve;
use occt_geom::extrema_pc::extrema_ext_pc_range_done;
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::chfi3d::builder_0_helpers::ch_fi3d_compute_arete;
use crate::chfi3d::builder_0_helpers::ch_fi3d_eval_tol_reached;
use crate::chfi3d::builder_0_helpers::ComputeArete;
use crate::geom_int::IntSS;
use crate::topopebrep::tool::adaptor::GeomAdaptorCurve;

use super::super::ChFi3dBuilder;
use super::MoreSurfdata;

/// `const double aTol3d = 1.e-4;` (`cxx:3777`).
pub(super) const TOL3D: f64 = 1.0e-4;

/// The locals `cxx:3918-4156` produce for the storage block.
pub(super) struct Section {
    /// `aPCurv1` — the p-curve of `aCracc` on the surface of `anInd`.
    ///
    /// `aCracc`, `aParf` and `aParl` of `cxx:3922-3957` are not carried over:
    /// the source reads `aCracc` only through the trimming of `cxx:4097-4155`
    /// and never reads `aParf` / `aParl` at all.
    pub pcurv1: Option<Arc<dyn Curve2d>>,
    /// `aTolReached`
    pub tol_reached: f64,
    /// `aCint1` — the `aSurfPrev` / `aSurf(anInd)` section through `aCP1onArc`.
    pub cint1: Arc<dyn Curve>,
    /// `aPCint11` — its p-curve on `aSurfPrev`.
    pub pcint11: Arc<dyn Curve2d>,
    /// `aPCint12` — its p-curve on `aSurf(anInd)`.
    pub pcint12: Arc<dyn Curve2d>,
    /// `aTolex1` — the evaluated tolerance of `aCint1`.
    pub tolex1: f64,
    /// `aTrCracc` — `aCracc` trimmed to the range the storage block uses.
    pub trcracc: Arc<dyn Curve>,
    /// `aCint2` — only when `is_do_second_section`.
    pub cint2: Option<Arc<dyn Curve>>,
    /// `aPCint21`
    pub pcint21: Option<Arc<dyn Curve2d>>,
    /// `aPCint22`
    pub pcint22: Option<Arc<dyn Curve2d>>,
    /// `aTolex2`
    pub tolex2: f64,
    /// `aPext` — the section end that is not on either arc point.
    pub pext: GpPnt,
    /// `isDoSecondSection`
    pub is_do_second_section: bool,
    /// `anOrSD1` — `Orientation()` of `SeqFil(anInd)` before the `Remove`.
    pub or_sd1: Orientation,
    /// `anOrSD2` — `Orientation()` of `SeqFil(anIndPrev)`.
    pub or_sd2: Orientation,
}

/// `cxx:3918-4156`. `Ok(None)` is one of the four bare `return`s.
pub(super) fn compute(b: &mut ChFi3dBuilder, st: &MoreSurfdata) -> Result<Option<Section>, String> {
    let tol3d = TOL3D;
    // `aSurfData = aSeqSurfData.Value(anInd);` (`cxx:3920`).
    let fd = super::fd_of(st)?;

    // `ChFi3d_ComputeArete(aSurfData->VertexLastOnS1(),
    //  aSurfData->InterferenceOnS1().PCurveOnSurf()->Value(...LastParameter()),
    //  aSurfData->VertexLastOnS2(), ...FirstParameter..., DStr.Surface(aSurfData
    //  ->Surf()).Surface(), aCracc, aPCurv1, aParf, aParl, aTol3d, tol2d,
    //  aTolReached, 0)` (`cxx:3922-3957`).
    let (cv1, uv1, cv2, uv2) = if st.is_first {
        (
            fd.vertex_last_on_s1().clone(),
            last_pcurve_point(fd.interference_on_s1())?,
            fd.vertex_last_on_s2().clone(),
            last_pcurve_point(fd.interference_on_s2())?,
        )
    } else {
        (
            fd.vertex_first_on_s1().clone(),
            first_pcurve_point(fd.interference_on_s1())?,
            fd.vertex_first_on_s2().clone(),
            first_pcurve_point(fd.interference_on_s2())?,
        )
    };
    let a_surf = {
        let ds = b.my_ds.change_ds();
        ds.surface(fd.surf()).surface().cloned()
    }
    .ok_or_else(|| "PerformMoreSurfdata: DStr.Surface(Fd->Surf()) is null".to_string())?;
    let mut arete = ComputeArete::default();
    ch_fi3d_compute_arete(
        &cv1,
        &uv1,
        &cv2,
        &uv2,
        &a_surf,
        &mut arete,
        tol3d,
        b.tol2d,
        0,
    )?;
    let cracc = arete
        .c3d
        .clone()
        .ok_or_else(|| "PerformMoreSurfdata: ChFi3d_ComputeArete returned no C3d".to_string())?;

    // `aSurf = DStr.Surface(aSurfData->Surf()).Surface();` (`cxx:3984`) — the
    // same surface the arete used.
    // `GeomInt_IntSS anInterSS(aSurfPrev, aSurf, 1.e-7, true, true, true);`
    // (`cxx:3986-3987`).
    let mut inter_ss = IntSS::new();
    inter_ss.perform(
        Arc::clone(&st.surf_prev),
        Arc::clone(&a_surf),
        1.0e-7,
        true,
        true,
        true,
    );
    // `if (!anInterSS.IsDone()) return;` (`cxx:4000-4003`).
    if !inter_ss.is_done() {
        return Ok(None);
    }

    // `for (i = 1; i <= anInterSS.NbLines() && !isFound; i++)` (`cxx:4007`).
    let mut cint1 = None;
    let mut pcint11 = None;
    let mut pcint12 = None;
    let mut tolex1 = 0.0;
    let mut pext1 = GpPnt::new(0.0, 0.0, 0.0);
    let mut pext2 = GpPnt::new(0.0, 0.0, 0.0);
    let mut is_found = false;
    let mut i = 1;
    while i <= inter_ss.nb_lines() && !is_found {
        let line = inter_ss
            .line(i)
            .ok_or_else(|| "PerformMoreSurfdata: IntSS line is missing".to_string())?;
        let c = Arc::clone(&line.curve);
        let p1 = line
            .pcurve1
            .clone()
            .ok_or_else(|| "PerformMoreSurfdata: IntSS has no p-curve on S1".to_string())?;
        let p2 = line
            .pcurve2
            .clone()
            .ok_or_else(|| "PerformMoreSurfdata: IntSS has no p-curve on S2".to_string())?;
        // `aTolex1 = ChFi3d_EvalTolReached(H1, aPCint11, H2, aPCint12, aCint1);`
        // (`cxx:4012`).
        tolex1 = ch_fi3d_eval_tol_reached(
            st.surf_prev.as_ref(),
            p1.as_ref(),
            a_surf.as_ref(),
            p2.as_ref(),
            c.as_ref(),
        );
        // `aCint1->D0(First/LastParameter(), aPext1 / aPext2);` (`cxx:4014-4015`).
        pext1 = c.d0(c.first_parameter());
        pext2 = c.d0(c.last_parameter());
        // `if (aPext1.Distance(aCP1onArc.Point()) <= aTol3d || aPext2.Distance(
        //  aCP1onArc.Point()) <= aTol3d) isFound = true;` (`cxx:4020-4025`).
        if pext1.distance(st.cp1_on_arc.point()) <= tol3d
            || pext2.distance(st.cp1_on_arc.point()) <= tol3d
        {
            is_found = true;
        }
        cint1 = Some(c);
        pcint11 = Some(p1);
        pcint12 = Some(p2);
        i += 1;
    }
    // `if (!isFound) return;` (`cxx:4028-4031`).
    if !is_found {
        return Ok(None);
    }
    let cint1 = cint1
        .ok_or_else(|| "PerformMoreSurfdata: IntSS produced no line".to_string())?;
    let pcint11 =
        pcint11.ok_or_else(|| "PerformMoreSurfdata: IntSS produced no line".to_string())?;
    let pcint12 =
        pcint12.ok_or_else(|| "PerformMoreSurfdata: IntSS produced no line".to_string())?;

    // `cxx:4033-4048`: the section end that is on neither arc point.
    let mut pext = GpPnt::new(0.0, 0.0, 0.0);
    let mut is_pext_found = false;
    if pext1.distance(st.cp2_on_arc.point()) > tol3d
        && pext1.distance(st.cp1_on_arc.point()) > tol3d
    {
        pext = pext1;
        is_pext_found = true;
    } else if pext2.distance(st.cp2_on_arc.point()) > tol3d
        && pext2.distance(st.cp1_on_arc.point()) > tol3d
    {
        pext = pext2;
        is_pext_found = true;
    }

    // `bool isDoSecondSection = false; double aPar = 0.;` (`cxx:4049-4050`).
    let mut is_do_second_section = false;
    let mut par = 0.0;
    if is_pext_found {
        // `GeomAdaptor_Curve aCad(aCracc); Extrema_ExtPC anExt(aPext, aCad,
        //  aTol3d);` (`cxx:4055-4056`).
        let (sols, done) = extrema_ext_pc_range_done(
            cracc.as_ref(),
            &pext,
            cracc.first_parameter(),
            cracc.last_parameter(),
        );
        // `if (!anExt.IsDone()) return;` (`cxx:4058-4061`).
        if !done {
            return Ok(None);
        }
        // `cxx:4062-4076`: see the module note — no early exit despite
        // `&& !isFound`.
        for sol in &sols {
            if sol.is_min && pext.distance(&sol.point) <= tol3d {
                par = sol.u;
                is_do_second_section = true;
            }
        }
    }

    // `anOrSD1 = aSurfData->Orientation(); aSurfData = aSeqSurfData.Value(
    //  anIndPrev); anOrSD2 = aSurfData->Orientation(); aSurf = DStr.Surface(
    //  aSurfData->Surf()).Surface();` (`cxx:4083-4087`).
    let or_sd1 = fd.orientation();
    let sd_prev = super::fd_at(st, st.ind_prev)?;
    let or_sd2 = sd_prev.orientation();
    let a_surf_prev_fillet = {
        let ds = b.my_ds.change_ds();
        ds.surface(sd_prev.surf()).surface().cloned()
    }
    .ok_or_else(|| "PerformMoreSurfdata: DStr.Surface(SDprev->Surf()) is null".to_string())?;

    let mut cint2 = None;
    let mut pcint21 = None;
    let mut pcint22 = None;
    let mut tolex2 = 0.0;
    let trcracc: Arc<dyn Curve>;
    if is_do_second_section {
        // `aCracc->D0(aCracc->FirstParameter(), aPext1);` (`cxx:4097`).
        let pext1_first = cracc.d0(cracc.first_parameter());
        let par1 = if pext1_first.distance(st.cp2_not_on_arc.point()) <= tol3d {
            cracc.first_parameter()
        } else {
            cracc.last_parameter()
        };
        // `aTrCracc = new Geom_TrimmedCurve(aCracc, aPar1, aPar / aPar, aPar1);`
        // (`cxx:4109-4115`).
        trcracc = if par1 < par {
            Arc::new(GeomAdaptorCurve::with_range(Arc::clone(&cracc), par1, par))
        } else {
            Arc::new(GeomAdaptorCurve::with_range(Arc::clone(&cracc), par, par1))
        };

        // Second section (`cxx:4118`).
        let mut inter_ss2 = IntSS::new();
        inter_ss2.perform(
            Arc::clone(&st.surf_prev),
            Arc::clone(&a_surf_prev_fillet),
            1.0e-7,
            true,
            true,
            true,
        );
        // `if (!anInterSS2.IsDone()) return;` (`cxx:4121-4124`).
        if !inter_ss2.is_done() {
            return Ok(None);
        }
        // `for (i = 1; i <= anInterSS2.NbLines() && !isFound; i++)`
        // (`cxx:4132-4148`).
        let mut is_found2 = false;
        let mut i = 1;
        while i <= inter_ss2.nb_lines() && !is_found2 {
            let line = inter_ss2
                .line(i)
                .ok_or_else(|| "PerformMoreSurfdata: IntSS2 line is missing".to_string())?;
            let c = Arc::clone(&line.curve);
            let p1 = line
                .pcurve1
                .clone()
                .ok_or_else(|| "PerformMoreSurfdata: IntSS2 has no p-curve on S1".to_string())?;
            let p2 = line
                .pcurve2
                .clone()
                .ok_or_else(|| "PerformMoreSurfdata: IntSS2 has no p-curve on S2".to_string())?;
            tolex2 = ch_fi3d_eval_tol_reached(
                st.surf_prev.as_ref(),
                p1.as_ref(),
                a_surf_prev_fillet.as_ref(),
                p2.as_ref(),
                c.as_ref(),
            );
            let e1 = c.d0(c.first_parameter());
            let e2 = c.d0(c.last_parameter());
            if e1.distance(st.cp2_on_arc.point()) <= tol3d
                || e2.distance(st.cp2_on_arc.point()) <= tol3d
            {
                is_found2 = true;
            }
            cint2 = Some(c);
            pcint21 = Some(p1);
            pcint22 = Some(p2);
            i += 1;
        }
        // `if (!isFound) return;` (`cxx:4147-4151`).
        if !is_found2 {
            return Ok(None);
        }
    } else {
        // `aTrCracc = new Geom_TrimmedCurve(aCracc, FirstParameter(),
        //  LastParameter());` (`cxx:4154-4155`).
        trcracc = Arc::new(GeomAdaptorCurve::with_range(
            Arc::clone(&cracc),
            cracc.first_parameter(),
            cracc.last_parameter(),
        ));
    }

    Ok(Some(Section {
        pcurv1: arete.p_curve.clone(),
        tol_reached: arete.tolreached,
        cint1,
        pcint11,
        pcint12,
        tolex1,
        trcracc,
        cint2,
        pcint21,
        pcint22,
        tolex2,
        pext,
        is_do_second_section,
        or_sd1,
        or_sd2,
    }))
}

/// `aSurfData->InterferenceOnS1().PCurveOnSurf()->Value(LastParameter())`.
fn last_pcurve_point(fi: &crate::chfi_ds::FaceInterference) -> Result<occt_core::gp::GpPnt2d, String> {
    fi.p_curve_on_surf()
        .map(|c| c.value(fi.last_parameter()))
        .ok_or_else(|| "PerformMoreSurfdata: InterferenceOnS* has no p-curve".to_string())
}

/// `aSurfData->InterferenceOnS1().PCurveOnSurf()->Value(FirstParameter())`.
fn first_pcurve_point(
    fi: &crate::chfi_ds::FaceInterference,
) -> Result<occt_core::gp::GpPnt2d, String> {
    fi.p_curve_on_surf()
        .map(|c| c.value(fi.first_parameter()))
        .ok_or_else(|| "PerformMoreSurfdata: InterferenceOnS* has no p-curve".to_string())
}
