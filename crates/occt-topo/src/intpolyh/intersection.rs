//! `IntPolyh_Intersection` (`IntPolyh_Intersection.hxx:45-119`,
//! `.cxx:33-568`): API class of the surface/surface intersection by
//! triangulation.
//!
//! The intersection runs inside the constructors, as in the `.cxx`: a caller
//! checks [`IntPolyhIntersection::is_done`] before reading the section lines.

use occt_geom::Surface;

use super::arrays::{ArrayOfPointNormal, ArrayOfSectionLines, ArrayOfTangentZones};
use super::edge_couple::IntPolyhCouple;
use super::maillage_affinage::MaillageAffinage;
use super::start_point::IntPolyhStartPoint;

/// ~~cos of 5 deg~~, the `anEps` of `IsAdvRequired` (`cxx:488`) and the `eps` of
/// `AnalyzeIntersection` (`cxx:547`).
const ANGLE_EPS: f64 = 0.996;

/// The eight values `GetLinePoint` hands back (`cxx:95-117`).
#[derive(Debug, Clone, Copy, Default)]
pub struct IntPolyhLinePoint {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub u1: f64,
    pub v1: f64,
    pub u2: f64,
    pub v2: f64,
    /// `SP.GetAngle()`.
    pub incidence: f64,
}

/// The seven values `GetTangentZonePoint` hands back (`cxx:122-139`).
#[derive(Debug, Clone, Copy, Default)]
pub struct IntPolyhTangentPoint {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub u1: f64,
    pub v1: f64,
    pub u2: f64,
    pub v2: f64,
}

/// `IntPolyh_Intersection` (`hxx:45`).
pub struct IntPolyhIntersection<'a> {
    /// `mySurf1` / `mySurf2` (`hxx:143-144`).
    surface1: &'a dyn Surface,
    surface2: &'a dyn Surface,
    /// `myNbSU1` / `myNbSV1` / `myNbSU2` / `myNbSV2` (`hxx:145-148`).
    nb_su1: i32,
    nb_sv1: i32,
    nb_su2: i32,
    nb_sv2: i32,
    /// `myIsDone` / `myIsParallel` (`hxx:149-150`).
    is_done: bool,
    is_parallel: bool,
    /// `mySectionLines` / `myTangentZones` (`hxx:151-152`).
    section_lines: ArrayOfSectionLines,
    tangent_zones: ArrayOfTangentZones,
}

impl<'a> IntPolyhIntersection<'a> {
    /// `IntPolyh_Intersection(theS1, theS2)` (`cxx:35-48`): 10x10 sampling on
    /// both surfaces, intersection performed.
    pub fn new(surface1: &'a dyn Surface, surface2: &'a dyn Surface) -> Self {
        Self::with_sampling(surface1, 10, 10, surface2, 10, 10)
    }

    /// `IntPolyh_Intersection(theS1, theNbSU1, theNbSV1, theS2, theNbSU2,
    /// theNbSV2)` (`cxx:52-70`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_sampling(
        surface1: &'a dyn Surface,
        nb_su1: i32,
        nb_sv1: i32,
        surface2: &'a dyn Surface,
        nb_su2: i32,
        nb_sv2: i32,
    ) -> Self {
        let mut this = Self::unperformed(surface1, nb_su1, nb_sv1, surface2, nb_su2, nb_sv2);
        this.perform();
        this
    }

    /// `IntPolyh_Intersection(theS1, theUPars1, theVPars1, theS2, theUPars2,
    /// theVPars2)` (`cxx:74-92`): the sampling is given by the caller.
    #[allow(clippy::too_many_arguments)]
    pub fn with_parameters(
        surface1: &'a dyn Surface,
        u_pars1: &[f64],
        v_pars1: &[f64],
        surface2: &'a dyn Surface,
        u_pars2: &[f64],
        v_pars2: &[f64],
    ) -> Self {
        let mut this = Self::unperformed(
            surface1,
            u_pars1.len() as i32,
            v_pars1.len() as i32,
            surface2,
            u_pars2.len() as i32,
            v_pars2.len() as i32,
        );
        this.perform_with_params(u_pars1, v_pars1, u_pars2, v_pars2);
        this
    }

    /// The state the three `.cxx` constructors establish before `Perform`
    /// (`cxx:37-47`, `:57-69`, `:81-91`).
    fn unperformed(
        surface1: &'a dyn Surface,
        nb_su1: i32,
        nb_sv1: i32,
        surface2: &'a dyn Surface,
        nb_su2: i32,
        nb_sv2: i32,
    ) -> Self {
        let mut section_lines = ArrayOfSectionLines::default();
        section_lines.init(1000);
        let mut tangent_zones = ArrayOfTangentZones::default();
        tangent_zones.init(10000);
        Self {
            surface1,
            surface2,
            nb_su1,
            nb_sv1,
            nb_su2,
            nb_sv2,
            is_done: false,
            is_parallel: false,
            section_lines,
            tangent_zones,
        }
    }

    /// `IsDone()` (`hxx:79`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `IsParallel()` (`hxx:82`).
    pub fn is_parallel(&self) -> bool {
        self.is_parallel
    }

    /// `NbSectionLines()` (`hxx:85`).
    pub fn nb_section_lines(&self) -> i32 {
        self.section_lines.nb_items()
    }

    /// `NbPointsInLine(IndexLine)` (`hxx:87-90`).
    pub fn nb_points_in_line(&self, index_line: i32) -> i32 {
        self.section_lines.get(index_line - 1).nb_start_points()
    }

    /// `NbTangentZones()` (`hxx:93`).
    pub fn nb_tangent_zones(&self) -> i32 {
        self.tangent_zones.nb_items()
    }

    /// `NbPointsInTangentZone()` (`hxx:96`): always 1 in the `.cxx`.
    pub fn nb_points_in_tangent_zone(&self) -> i32 {
        1
    }

    /// `GetLinePoint(IndexLine, IndexPoint, ...)` (`cxx:95-117`).
    pub fn get_line_point(&self, index_line: i32, index_point: i32) -> IntPolyhLinePoint {
        let sp = self.section_lines.get(index_line - 1).value(index_point - 1);
        IntPolyhLinePoint {
            x: sp.x,
            y: sp.y,
            z: sp.z,
            u1: sp.u1,
            v1: sp.v1,
            u2: sp.u2,
            v2: sp.v2,
            incidence: sp.angle,
        }
    }

    /// `GetTangentZonePoint(IndexLine, IndexPoint, ...)` (`cxx:122-139`): the
    /// `IndexPoint` argument is ignored, as in the `.cxx`.
    pub fn get_tangent_zone_point(&self, index_line: i32) -> IntPolyhTangentPoint {
        let sp: IntPolyhStartPoint = self.tangent_zones.get(index_line - 1);
        IntPolyhTangentPoint {
            x: sp.x,
            y: sp.y,
            z: sp.z,
            u1: sp.u1,
            v1: sp.v1,
            u2: sp.u2,
            v2: sp.v2,
        }
    }

    /// The section lines, for the callers that walk them directly.
    pub fn section_lines(&self) -> &ArrayOfSectionLines {
        &self.section_lines
    }

    /// The tangent zones, for the callers that walk them directly.
    pub fn tangent_zones(&self) -> &ArrayOfTangentZones {
        &self.tangent_zones
    }

    /// `Perform()` (`cxx:144-150`).
    fn perform(&mut self) {
        let (u_pars1, v_pars1) =
            super::tools::make_sampling(self.surface1, self.nb_su1, self.nb_sv1, false);
        let (u_pars2, v_pars2) =
            super::tools::make_sampling(self.surface2, self.nb_su2, self.nb_sv2, false);
        self.perform_with_params(&u_pars1, &v_pars1, &u_pars2, &v_pars2);
    }

    /// `Perform(theUPars1, theVPars1, theUPars2, theVPars2)` (`cxx:155-240`).
    fn perform_with_params(
        &mut self,
        u_pars1: &[f64],
        v_pars1: &[f64],
        u_pars2: &[f64],
        v_pars2: &[f64],
    ) {
        self.is_done = true;

        let defl_tol1 = super::tools::compute_deflection(self.surface1, u_pars1, v_pars1);
        let defl_tol2 = super::tools::compute_deflection(self.surface2, u_pars2, v_pars2);

        // Standard intersection (`cxx:168-178`).
        let (is_std_done, mut maillage_std) =
            self.perform_maillage(u_pars1, v_pars1, u_pars2, v_pars2, defl_tol1, defl_tol2);
        if !is_std_done {
            self.is_done = false;
            return;
        }
        let nb_couples_std = maillage_std.get_couples().len() as i32;

        if !self.is_adv_required(&maillage_std) {
            maillage_std.start_points_chain(&mut self.section_lines, &mut self.tangent_zones);
            return;
        }

        // Advanced intersection (`cxx:191-238`).
        let mut points1 = ArrayOfPointNormal::default();
        let mut points2 = ArrayOfPointNormal::default();
        super::tools::fill_array_of_point_normal(self.surface1, u_pars1, v_pars1, &mut points1);
        super::tools::fill_array_of_point_normal(self.surface2, u_pars2, v_pars2, &mut points2);

        let mut maillage_fr = None;
        let mut maillage_rf = None;
        let mut maillage_ff = None;
        let mut maillage_rr = None;
        let mut is_adv_done = true;

        // `cxx:291-328`: the shifts are evaluated FR, RF, FF, RR and the chain
        // stops at the first failure.
        if is_adv_done {
            let (ok, m) = self.perform_maillage_normals(
                u_pars1, v_pars1, u_pars2, v_pars2, defl_tol1, defl_tol2, &points1, &points2,
                true, false,
            );
            is_adv_done = ok;
            if ok {
                maillage_fr = Some(m);
            }
        }
        if is_adv_done {
            let (ok, m) = self.perform_maillage_normals(
                u_pars1, v_pars1, u_pars2, v_pars2, defl_tol1, defl_tol2, &points1, &points2,
                false, true,
            );
            is_adv_done = ok;
            if ok {
                maillage_rf = Some(m);
            }
        }
        if is_adv_done {
            let (ok, m) = self.perform_maillage_normals(
                u_pars1, v_pars1, u_pars2, v_pars2, defl_tol1, defl_tol2, &points1, &points2,
                true, true,
            );
            is_adv_done = ok;
            if ok {
                maillage_ff = Some(m);
            }
        }
        if is_adv_done {
            let (ok, m) = self.perform_maillage_normals(
                u_pars1, v_pars1, u_pars2, v_pars2, defl_tol1, defl_tol2, &points1, &points2,
                false, false,
            );
            is_adv_done = ok;
            if ok {
                maillage_rr = Some(m);
            }
        }

        let nb_couples_adv;
        if is_adv_done {
            let (mut ff, mut fr, mut rf, mut rr) = (
                maillage_ff.take().unwrap(),
                maillage_fr.take().unwrap(),
                maillage_rf.take().unwrap(),
                maillage_rr.take().unwrap(),
            );
            nb_couples_adv = ff.get_couples().len() as i32
                + fr.get_couples().len() as i32
                + rf.get_couples().len() as i32
                + rr.get_couples().len() as i32;

            if nb_couples_adv > 0 {
                Self::merge_couples(
                    ff.get_couples_mut(),
                    fr.get_couples_mut(),
                    rf.get_couples_mut(),
                    rr.get_couples_mut(),
                );
            }

            if nb_couples_adv > 0 {
                // `cxx:218-221`: the chain order is FF, FR, RF, RR.
                ff.start_points_chain(&mut self.section_lines, &mut self.tangent_zones);
                fr.start_points_chain(&mut self.section_lines, &mut self.tangent_zones);
                rf.start_points_chain(&mut self.section_lines, &mut self.tangent_zones);
                rr.start_points_chain(&mut self.section_lines, &mut self.tangent_zones);
                return;
            }
        }

        // `cxx:224-228`: fall back to the standard intersection.
        if nb_couples_std > 0 {
            maillage_std.start_points_chain(&mut self.section_lines, &mut self.tangent_zones);
        }
    }

    /// `PerformMaillage(...)` (`cxx:351-394`), the standard (unshifted) one.
    #[allow(clippy::too_many_arguments)]
    fn perform_maillage(
        &mut self,
        u_pars1: &[f64],
        v_pars1: &[f64],
        u_pars2: &[f64],
        v_pars2: &[f64],
        defl_tol1: f64,
        defl_tol2: f64,
    ) -> (bool, MaillageAffinage<'a>) {
        let mut maillage = MaillageAffinage::with_sampling(
            self.surface1,
            u_pars1.len() as i32,
            v_pars1.len() as i32,
            self.surface2,
            u_pars2.len() as i32,
            v_pars2.len() as i32,
        );

        maillage.fill_array_of_pnt_with_params(1, u_pars1, v_pars1, Some(defl_tol1));
        maillage.fill_array_of_pnt_with_params(2, u_pars2, v_pars2, Some(defl_tol2));

        let fin_ttc = compute_intersection(&mut maillage);

        // If no intersecting triangles are found, try the enlarged surfaces
        // (`cxx:373-390`).
        if fin_ttc == 0 {
            let (enlarge_u1, enlarge_v1) = super::tools::is_enlarge_possible(self.surface1);
            let (enlarge_u2, enlarge_v2) = super::tools::is_enlarge_possible(self.surface2);
            if enlarge_u1 || enlarge_v1 || enlarge_u2 || enlarge_v2 {
                maillage.set_enlarge_zone(true);
                maillage.fill_array_of_pnt(1);
                maillage.fill_array_of_pnt(2);
                compute_intersection(&mut maillage);
                maillage.set_enlarge_zone(false);
            }
        }

        let is_done = self.analyze_intersection(&mut maillage);
        (is_done, maillage)
    }

    /// `PerformMaillage(..., thePoints1, thePoints2, theIsFirstFwd,
    /// theIsSecondFwd, ...)` (`cxx:398-422`), the shifted one.
    #[allow(clippy::too_many_arguments)]
    fn perform_maillage_normals(
        &mut self,
        u_pars1: &[f64],
        v_pars1: &[f64],
        u_pars2: &[f64],
        v_pars2: &[f64],
        defl_tol1: f64,
        defl_tol2: f64,
        points1: &ArrayOfPointNormal,
        points2: &ArrayOfPointNormal,
        is_first_fwd: bool,
        is_second_fwd: bool,
    ) -> (bool, MaillageAffinage<'a>) {
        let mut maillage = MaillageAffinage::with_sampling(
            self.surface1,
            u_pars1.len() as i32,
            v_pars1.len() as i32,
            self.surface2,
            u_pars2.len() as i32,
            v_pars2.len() as i32,
        );

        maillage.fill_array_of_pnt_with_normals(1, is_first_fwd, points1, u_pars1, v_pars1, defl_tol1);
        maillage.fill_array_of_pnt_with_normals(2, is_second_fwd, points2, u_pars2, v_pars2, defl_tol2);

        compute_intersection(&mut maillage);

        let is_done = self.analyze_intersection(&mut maillage);
        (is_done, maillage)
    }

    /// `MergeCouples(anArrayFF, anArrayFR, anArrayRF, anArrayRR)`
    /// (`cxx:433-453`): a couple is kept in the first list that holds it.
    fn merge_couples(
        array_ff: &mut Vec<IntPolyhCouple>,
        array_fr: &mut Vec<IntPolyhCouple>,
        array_rf: &mut Vec<IntPolyhCouple>,
        array_rr: &mut Vec<IntPolyhCouple>,
    ) {
        let mut fence: std::collections::HashSet<(i32, i32)> = std::collections::HashSet::new();
        for list in [array_ff, array_fr, array_rf, array_rr] {
            list.retain(|couple| fence.insert(couple.hash_key()));
        }
    }

    /// `IsAdvRequired(theMaillage)` (`cxx:461-502`).
    fn is_adv_required(&self, maillage: &MaillageAffinage<'a>) -> bool {
        let couples = maillage.get_couples();
        let nb_couples = couples.len() as i32;

        let mut is_adv_req = nb_couples == 0 && !self.is_parallel;
        if is_adv_req {
            return is_adv_req;
        }

        if nb_couples > 10 {
            return is_adv_req;
        }

        for couple in couples {
            if couple.angle.abs() > ANGLE_EPS {
                is_adv_req = true;
                break;
            }
        }
        is_adv_req
    }

    /// `AnalyzeIntersection(theMaillage)` (`cxx:537-566`).
    fn analyze_intersection(&mut self, maillage: &mut MaillageAffinage<'a>) -> bool {
        let fin_ttc = maillage.get_couples().len() as i32;
        if fin_ttc > 200 {
            let npara = maillage
                .get_couples()
                .iter()
                .filter(|couple| couple.angle.abs() > ANGLE_EPS)
                .count() as i32;
            if npara >= maillage.get_array_of_triangles(1).nb_items()
                || npara >= maillage.get_array_of_triangles(2).nb_items()
            {
                maillage.get_couples_mut().clear();
                self.is_parallel = true;
                return true;
            }
        }
        true
    }
}

/// `ComputeIntersection(theMaillage)` (`cxx:509-530`): box, triangles, edges,
/// deflection refinement, then the triangle/couple search.
fn compute_intersection(maillage: &mut MaillageAffinage<'_>) -> i32 {
    maillage.common_box();
    maillage.fill_array_of_triangles(1);
    maillage.fill_array_of_triangles(2);
    maillage.fill_array_of_edges(1);
    maillage.fill_array_of_edges(2);
    maillage.triangles_deflections_refinement_bsb();
    maillage.triangle_compare()
}
