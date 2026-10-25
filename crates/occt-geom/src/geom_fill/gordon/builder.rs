//! `GeomFill_GordonBuilder` — the Gordon surface builder.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/GeomFill/GeomFill_GordonBuilder.cxx`
//! (OCCT 8.0.0).
//!
//! A Gordon surface interpolates a network of profiles and guides: the profiles
//! are skinned into `S_profiles`, the guides into `S_guides`, the profile/guide
//! intersection points into a tensor surface `S_tensor`, and the result is the
//! boolean sum `S_profiles + S_guides - S_tensor`.
//!
//! Translation notes:
//! - OCCT keeps every surface as a handle and `Surface()` /
//!   `ProfileSurface()` / `GuideSurface()` / `TensorSurface()` throw
//!   `StdFail_NotDone` before `Perform` (`cxx:152-191`); the port returns
//!   `None`.
//! - `myToUseParallel` (`cxx:264-299`) only selects the `OSD_Parallel::For`
//!   policy; the port runs the loop serially, since the approximation code has
//!   no parallel executor. The results do not depend on the policy — every pole
//!   column is interpolated independently.
//! - `weave` of the boolean sum uses `WeightsArray()`: OCCT stores unit weights
//!   for a non-rational surface, so [`super::super::bspline_surface::GeomBSplineSurface::weight`]
//!   returns 1.0 there.
//! - `StdFail_NotDone` and `Standard_ConstructionError` from the interpolation
//!   and knot operations become the `Option`/`Err` returns documented below.

use crate::bspline_curve::GeomBSplineCurve;
use crate::bspline_surface::GeomBSplineSurface;
use crate::curve::Curve;
use crate::geom_api_interpolate::GeomAPIInterpolate;
use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, PCONFUSION};

/// `GeomFill_GordonBuilder` (`GeomFill_GordonBuilder.hxx`).
pub struct GordonBuilder {
    /// `myProfiles`.
    profiles: Vec<GeomBSplineCurve>,
    /// `myGuides`.
    guides: Vec<GeomBSplineCurve>,
    /// `myProfileParams`.
    profile_params: Vec<f64>,
    /// `myGuideParams`.
    guide_params: Vec<f64>,
    /// `myTolerance`.
    tolerance: f64,
    /// `myIsUClosed`.
    is_u_closed: bool,
    /// `myIsVClosed`.
    is_v_closed: bool,
    /// `myToUseParallel`.
    to_use_parallel: bool,
    /// `myIsDone`.
    is_done: bool,
    /// `mySurface`.
    surface: Option<GeomBSplineSurface>,
    /// `myProfileSurface`.
    profile_surface: Option<GeomBSplineSurface>,
    /// `myGuideSurface`.
    guide_surface: Option<GeomBSplineSurface>,
    /// `myTensorSurface`.
    tensor_surface: Option<GeomBSplineSurface>,
}

impl Default for GordonBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl GordonBuilder {
    /// `GeomFill_GordonBuilder::GeomFill_GordonBuilder()`
    /// (`GeomFill_GordonBuilder.cxx:28`).
    pub fn new() -> Self {
        Self {
            profiles: Vec::new(),
            guides: Vec::new(),
            profile_params: Vec::new(),
            guide_params: Vec::new(),
            tolerance: 0.0,
            is_u_closed: false,
            is_v_closed: false,
            to_use_parallel: false,
            is_done: false,
            surface: None,
            profile_surface: None,
            guide_surface: None,
            tensor_surface: None,
        }
    }

    /// `GeomFill_GordonBuilder::Init` (`GeomFill_GordonBuilder.cxx:32-53`).
    pub fn init(
        &mut self,
        profiles: Vec<GeomBSplineCurve>,
        guides: Vec<GeomBSplineCurve>,
        profile_params: Vec<f64>,
        guide_params: Vec<f64>,
        tolerance: f64,
        is_u_closed: bool,
        is_v_closed: bool,
    ) {
        self.profiles = profiles;
        self.guides = guides;
        self.profile_params = profile_params;
        self.guide_params = guide_params;
        self.tolerance = tolerance;
        self.is_u_closed = is_u_closed;
        self.is_v_closed = is_v_closed;
        self.is_done = false;
        self.surface = None;
        self.profile_surface = None;
        self.guide_surface = None;
        self.tensor_surface = None;
    }

    /// `GeomFill_GordonBuilder::SetParallelMode`
    /// (`GeomFill_GordonBuilder.hxx:74-76`). `OSD_Parallel::For`
    /// (`cxx:264-299`) only chooses the executor; the port keeps the per-column
    /// interpolation serial, so the flag is stored and reported but has no
    /// effect on the result.
    pub fn set_parallel_mode(&mut self, to_use_parallel: bool) {
        self.to_use_parallel = to_use_parallel;
    }

    /// `GeomFill_GordonBuilder::IsParallelMode`
    /// (`GeomFill_GordonBuilder.hxx:78-79`).
    pub fn is_parallel_mode(&self) -> bool {
        self.to_use_parallel
    }

    /// `GeomFill_GordonBuilder::Perform` (`GeomFill_GordonBuilder.cxx:57-146`).
    pub fn perform(&mut self) {
        self.is_done = false;

        if self.profiles.len() < 2 || self.guides.len() < 2 {
            return;
        }
        if self.profile_params.len() != self.profiles.len()
            || self.guide_params.len() != self.guides.len()
        {
            return;
        }
        // `cxx:74-86`: rational sections are rejected.
        if self.profiles.iter().any(|c| c.is_rational()) || self.guides.iter().any(|c| c.is_rational())
        {
            return;
        }

        // `cxx:92-105`: `intersection_pnts(guide, profile) =
        // profile.Value(guideParam)`.
        let nb_guides = self.guides.len();
        let nb_profiles = self.profiles.len();
        let mut intersections = vec![vec![GpPnt::zero(); nb_profiles]; nb_guides];
        for (p, profile) in self.profiles.iter().enumerate() {
            for (g, _guide) in self.guides.iter().enumerate() {
                let param = self.guide_params[g];
                intersections[g][p] = profile.d0(param);
            }
        }

        // Step 1: skin the profiles (V direction).
        let Some(profile_surface) =
            self.build_skin_surface(&self.profiles, &self.profile_params, self.is_v_closed)
        else {
            return;
        };
        self.profile_surface = Some(profile_surface);

        // Step 2: skin the guides (U direction), then swap U and V.
        let Some(mut guide_surface) =
            self.build_skin_surface(&self.guides, &self.guide_params, self.is_u_closed)
        else {
            return;
        };
        guide_surface.exchange_uv();
        self.guide_surface = Some(guide_surface);

        // Step 3: the tensor surface through the intersection grid.
        let Some(tensor_surface) = self.build_tensor_surface(
            &intersections,
            &self.guide_params,
            &self.profile_params,
            self.is_u_closed,
            self.is_v_closed,
        ) else {
            return;
        };
        self.tensor_surface = Some(tensor_surface);

        // Step 4: bring all three to common degrees and knots.
        let mut profile_surface = self.profile_surface.take().unwrap();
        let mut guide_surface = self.guide_surface.take().unwrap();
        let mut tensor_surface = self.tensor_surface.take().unwrap();
        let unified = unify_surfaces(&mut profile_surface, &mut guide_surface).is_ok()
            && unify_surfaces(&mut profile_surface, &mut tensor_surface).is_ok()
            && unify_surfaces(&mut guide_surface, &mut tensor_surface).is_ok();
        // The three surfaces stay reachable through the getters even when a
        // later step fails, as the OCCT handles do.
        self.profile_surface = Some(profile_surface);
        self.guide_surface = Some(guide_surface);
        self.tensor_surface = Some(tensor_surface);
        if !unified {
            return;
        }

        // Step 5: boolean sum.
        let profile_surface = self.profile_surface.as_ref().unwrap();
        let guide_surface = self.guide_surface.as_ref().unwrap();
        let tensor_surface = self.tensor_surface.as_ref().unwrap();
        let Some(surface) = compute_boolean_sum(profile_surface, guide_surface, tensor_surface)
        else {
            return;
        };
        self.surface = Some(surface);
        self.is_done = true;
    }

    /// `GeomFill_GordonBuilder::Surface()` (`GeomFill_GordonBuilder.cxx:152-158`);
    /// `None` where OCCT throws `StdFail_NotDone`.
    pub fn surface(&self) -> Option<&GeomBSplineSurface> {
        if self.is_done {
            self.surface.as_ref()
        } else {
            None
        }
    }

    /// `GeomFill_GordonBuilder::ProfileSurface()`
    /// (`GeomFill_GordonBuilder.cxx:162-168`).
    pub fn profile_surface(&self) -> Option<&GeomBSplineSurface> {
        if self.is_done {
            self.profile_surface.as_ref()
        } else {
            None
        }
    }

    /// `GeomFill_GordonBuilder::GuideSurface()`
    /// (`GeomFill_GordonBuilder.cxx:172-178`).
    pub fn guide_surface(&self) -> Option<&GeomBSplineSurface> {
        if self.is_done {
            self.guide_surface.as_ref()
        } else {
            None
        }
    }

    /// `GeomFill_GordonBuilder::TensorSurface()`
    /// (`GeomFill_GordonBuilder.cxx:182-191`).
    pub fn tensor_surface(&self) -> Option<&GeomBSplineSurface> {
        if self.is_done {
            self.tensor_surface.as_ref()
        } else {
            None
        }
    }

    /// `IsDone()` — the port's stand-in for a successful `Perform`.
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `GeomFill_GordonBuilder::buildSkinSurface`
    /// (`GeomFill_GordonBuilder.cxx:195-317`): interpolate every pole column of
    /// the sections along the section parameters and assemble the tensor grid.
    fn build_skin_surface(
        &self,
        sections: &[GeomBSplineCurve],
        section_params: &[f64],
        is_closed: bool,
    ) -> Option<GeomBSplineSurface> {
        if sections.len() < 2 {
            return None;
        }
        let nb_poles_u = sections[0].nb_poles();
        let nb_sections = sections.len();

        // `cxx:208-218`: closedness only when the first and last sections
        // coincide.
        let mut make_closed = is_closed;
        if make_closed && !sections[0].is_equal(&sections[nb_sections - 1], CONFUSION) {
            make_closed = false;
        }

        // `cxx:220-232`: periodic interpolation uses N-1 points with N
        // parameters (the wrap-around parameter closes the loop).
        let nb_interp_pts = if make_closed {
            nb_sections - 1
        } else {
            nb_sections
        };
        let nb_param_pts = nb_sections;
        let params = section_params[..nb_param_pts].to_vec();

        // `cxx:234-252`: interpolate the first pole column to fix the V knots.
        let first_points: Vec<GpPnt> = (0..nb_interp_pts)
            .map(|i| sections[i].poles.first().copied())
            .collect::<Option<Vec<_>>>()?;
        let first_curve = interpolate(&first_points, &params, make_closed)?;
        let nb_poles_v = first_curve.nb_poles();

        let mut poles = vec![vec![GpPnt::zero(); nb_poles_v]; nb_poles_u];
        for j in 0..nb_poles_v {
            poles[0][j] = *first_curve.pole(j);
        }

        // `cxx:262-299`: every remaining pole column. OCCT parallelizes this
        // loop; the port keeps the same per-column interpolation serially.
        for u in 1..nb_poles_u {
            let points: Vec<GpPnt> = (0..nb_interp_pts)
                .map(|i| sections[i].poles.get(u).copied())
                .collect::<Option<Vec<_>>>()?;
            let curve = interpolate(&points, &params, make_closed)?;
            for j in 0..nb_poles_v {
                poles[u][j] = *curve.pole(j);
            }
        }

        // `cxx:309-317`: U knots from the first section, V knots from the first
        // interpolation curve.
        GeomBSplineSurface::from_poles_knots_mults(
            poles,
            sections[0].distinct_knots().to_vec(),
            first_curve.distinct_knots().to_vec(),
            sections[0].multiplicities().to_vec(),
            first_curve.multiplicities().to_vec(),
            sections[0].degree(),
            first_curve.degree(),
        )
        .ok()
    }

    /// `GeomFill_GordonBuilder::buildTensorSurface`
    /// (`GeomFill_GordonBuilder.cxx:321-413`).
    fn build_tensor_surface(
        &self,
        points: &[Vec<GpPnt>],
        u_params: &[f64],
        v_params: &[f64],
        is_u_closed: bool,
        is_v_closed: bool,
    ) -> Option<GeomBSplineSurface> {
        let nb_u = points.len();
        let nb_v = points.first().map_or(0, |r| r.len());
        if nb_u < 2 || nb_v < 2 {
            return None;
        }

        // `cxx:338-351`: periodic U interpolation needs N-1 points with N
        // parameters, and only when the first and last rows coincide.
        let mut make_u_closed = is_u_closed && nb_u >= 3;
        if make_u_closed {
            for v in 0..nb_v {
                // `gp_Pnt::IsEqual(other, Precision::Confusion())`.
                if points[0][v].distance(&points[nb_u - 1][v]) > CONFUSION {
                    make_u_closed = false;
                    break;
                }
            }
        }
        let nb_interp_u = if make_u_closed { nb_u - 1 } else { nb_u };
        let u_param_values = u_params[..nb_u].to_vec();

        // `cxx:365-387`: interpolate each row to get a curve in U.
        let mut row_curves: Vec<GeomBSplineCurve> = Vec::with_capacity(nb_v);
        for v in 0..nb_v {
            let row_points: Vec<GpPnt> = (0..nb_interp_u).map(|u| points[u][v]).collect();
            row_curves.push(interpolate(&row_points, &u_param_values, make_u_closed)?);
        }

        // `cxx:389-395`: raise every row curve to the highest degree.
        let mut max_degree = 0usize;
        for c in &row_curves {
            max_degree = max_degree.max(c.degree());
        }
        for c in row_curves.iter_mut() {
            if c.increase_degree(max_degree).is_err() {
                return None;
            }
        }

        // `cxx:398-411`: unify the knot vectors before skinning. OCCT snapshots
        // each curve's `Knots()` / `Multiplicities()` (copies) and inserts them
        // into the first curve, then inserts the first curve's unified arrays
        // into all the others.
        let other_knots: Vec<(Vec<f64>, Vec<i32>)> = row_curves[1..]
            .iter()
            .map(|c| (c.distinct_knots().to_vec(), c.multiplicities().to_vec()))
            .collect();
        for (knots, mults) in &other_knots {
            if row_curves[0]
                .insert_knots(knots, Some(mults), PCONFUSION, false)
                .is_err()
            {
                return None;
            }
        }
        let unified_knots = row_curves[0].distinct_knots().to_vec();
        let unified_mults = row_curves[0].multiplicities().to_vec();
        for c in row_curves.iter_mut().skip(1) {
            if c.insert_knots(&unified_knots, Some(&unified_mults), PCONFUSION, false)
                .is_err()
            {
                return None;
            }
        }

        // `cxx:411`: skin the row curves in V.
        self.build_skin_surface(&row_curves, &v_params[..nb_v], is_v_closed)
    }
}

/// `GeomAPI_Interpolate` through `points` with explicit `params`
/// (`GeomFill_GordonBuilder.cxx:241-249`, `:282-287`, `:378-381`).
fn interpolate(
    points: &[GpPnt],
    params: &[f64],
    periodic: bool,
) -> Option<GeomBSplineCurve> {
    let mut interp = GeomAPIInterpolate::with_parameters(
        points.to_vec(),
        params.to_vec(),
        periodic,
        CONFUSION,
    )
    .ok()?;
    interp.perform();
    interp.curve().cloned()
}

/// `GeomFill_GordonBuilder::unifySurfaces` (`GeomFill_GordonBuilder.cxx:417-453`).
fn unify_surfaces(
    surf1: &mut GeomBSplineSurface,
    surf2: &mut GeomBSplineSurface,
) -> Result<(), String> {
    let deg_u = surf1.deg_u().max(surf2.deg_u());
    let deg_v = surf1.deg_v().max(surf2.deg_v());
    surf1.increase_degree(deg_u, deg_v)?;
    surf2.increase_degree(deg_u, deg_v)?;

    // `cxx:430-440`: U knots, using snapshots of both vectors because the
    // insertion rewrites the target's own arrays.
    let knots1 = surf1.uknots().to_vec();
    let mults1 = surf1.umults().to_vec();
    let knots2 = surf2.uknots().to_vec();
    let mults2 = surf2.umults().to_vec();
    surf1.insert_u_knots(&knots2, &mults2, PCONFUSION, false)?;
    surf2.insert_u_knots(&knots1, &mults1, PCONFUSION, false)?;

    // `cxx:443-452`: V knots, same two-pass approach.
    let knots1 = surf1.vknots().to_vec();
    let mults1 = surf1.vmults().to_vec();
    let knots2 = surf2.vknots().to_vec();
    let mults2 = surf2.vmults().to_vec();
    surf1.insert_v_knots(&knots2, &mults2, PCONFUSION, false)?;
    surf2.insert_v_knots(&knots1, &mults1, PCONFUSION, false)?;
    Ok(())
}

/// `GeomFill_GordonBuilder::computeBooleanSum`
/// (`GeomFill_GordonBuilder.cxx:457-502`).
///
/// The formula works in weighted (homogeneous) coordinates:
/// `w_gordon = w_prof + w_guide - w_tensor` and
/// `w_gordon * P_gordon = w_prof * P_prof + w_guide * P_guide - w_tensor * P_tensor`.
/// `None` where OCCT returns a null handle (a vanishing weight, `cxx:486-489`).
fn compute_boolean_sum(
    profile: &GeomBSplineSurface,
    guide: &GeomBSplineSurface,
    tensor: &GeomBSplineSurface,
) -> Option<GeomBSplineSurface> {
    let mut result = profile.clone();
    let nb_u = result.nb_poles_u();
    let nb_v = result.nb_poles_v();
    for u in 0..nb_u {
        for v in 0..nb_v {
            let w_prof = profile.weight(u, v);
            let w_guide = guide.weight(u, v);
            let w_tensor = tensor.weight(u, v);
            let w_gordon = w_prof + w_guide - w_tensor;
            if w_gordon.abs() < CONFUSION {
                return None;
            }
            let pp = profile.pole(u, v);
            let pg = guide.pole(u, v);
            let pt = tensor.pole(u, v);
            let x = (pp.x() * w_prof + pg.x() * w_guide - pt.x() * w_tensor) / w_gordon;
            let y = (pp.y() * w_prof + pg.y() * w_guide - pt.y() * w_tensor) / w_gordon;
            let z = (pp.z() * w_prof + pg.z() * w_guide - pt.z() * w_tensor) / w_gordon;
            result.set_pole(u, v, GpPnt::new(x, y, z));
            result.set_weight(u, v, w_gordon).ok()?;
        }
    }
    Some(result)
}
