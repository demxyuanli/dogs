//! `GeomPlate_BuildPlateSurface` (`GeomPlate_BuildPlateSurface.hxx:44-278`,
//! `GeomPlate_BuildPlateSurface.cxx:80-2809`): builds a plate surface that
//! deforms an initial surface so that it satisfies curve and point constraints.
//!
//! | OCCT | here |
//! |---|---|
//! | ctor `(NPoints, TabCurve, Tang, Degree, ...)` (`cxx:80-142`) | [`BuildPlateSurface::from_curves`] |
//! | ctor `(Surf, Degree, NbPtsOnCur, ...)` (`cxx:148-182`) | [`BuildPlateSurface::from_surface`] |
//! | ctor `(Degree, NbPtsOnCur, ...)` (`cxx:186-219`) | [`BuildPlateSurface::new`] |
//! | `Init` / `LoadInitSurface` / `Add` / `SetNbBounds` (`cxx:388-432`) | the same names |
//! | `Perform` (`cxx:438-747`) | [`BuildPlateSurface::perform`] in `perform.rs` |
//! | `EcartContraintesMil` (`cxx:751-865`) | `verify::ecart_contraintes_mil` |
//! | `ProjectPoint` / `ProjectCurve` / `ProjectedCurve` (`cxx:253-376`) | `projection.rs` |
//! | `ComputeSurfInit` (`cxx:1455-1903`) | `perform::compute_surf_init` |
//! | `Intersect` (`cxx:1909-2145`) | `intersect.rs` |
//! | `Discretise` / `CalculNbPtsInit` (`cxx:2156-2400`) | `discretise.rs` |
//! | `LoadCurve` / `LoadPoint` (`cxx:2404-2584`) | `load.rs` |
//! | `VerifSurface` / `VerifPoints` (`cxx:2587-2781`) | `verify.rs` |
//! | `CourbeJointive` (`cxx:1349-1444`) | `verify::courbe_jointive` |
//! | `Disc2dContour` / `Disc3dContour` (`cxx:869-1163`) | `contours.rs` |
//! | `Sense` / `Curves2d` / `Order` (`cxx:1191-1233`) | the same names |
//! | the `G0Error` / `G1Error` / `G2Error` family (`cxx:1236-1323`) | the same names |
//!
//! Translation notes:
//! - Every `occ::handle<NCollection_HArray1<T>>` becomes a `Vec<T>` with 0-based
//!   indices; the OCCT `1..=N` loops become `0..N` with the same arithmetic.
//! - `myLinCont` / `myPntCont` are `NCollection_HSequence<handle>` that the
//!   builder mutates in place (`SetValue`, `ChangeValue`, `Append`, `Clear`); the
//!   port stores the constraints by value in a `Vec` and the `CurveConstraint(order)`
//!   / `PointConstraint(order)` accessors return a borrow instead of a handle.
//!   No caller in OCCT writes through those accessors.
//! - `myProj` (`Extrema_ExtPS`) is a member that `Perform` re-initialises on the
//!   current `mySurfInit`; here the surface is a per-call argument of
//!   [`crate::extrema_surf::ExtPs`], so `project_point` builds the engine with
//!   the same `(Uinf, Usup, Vinf, Vsup, TolU, TolV)` window (`cxx:516-517`,
//!   `:542-543`, `:565-566`).
//! - `Standard_ConstructionError` (`cxx:106-131`, `:170-177`, `:205-212`) becomes
//!   `Err`, and `myNbIter < 1` / `myDegree < 2` are still the two checks.
//! - `GeomPlate_MakeApprox` calls (`cxx:537`, `:1890`) use `CritOrder == -1`, the
//!   only arm [`crate::geom_plate::MakeApprox`] ports.

mod contours;
mod discretise;
mod intersect;
mod load;
mod perform;
mod projection;
mod verify;

use std::sync::Arc;

use occt_core::gp::GpPnt;

use crate::Surface;
use crate::geom_plate::{BoundaryCurve, CurveConstraint, PlateSurface, PointConstraint};
use crate::plate::Plate;

pub use perform::PerformOutcome;

/// `GeomPlate_BuildPlateSurface` (`hxx:44-278`).
pub struct BuildPlateSurface {
    /// `myLinCont` (`hxx:248`).
    lin_cont: Vec<CurveConstraint>,
    /// `myParCont` (`hxx:249`): the parameters of the constraint points.
    par_cont: Vec<Vec<f64>>,
    /// `myPlateCont` (`hxx:250`): the parameters actually loaded into the plate.
    plate_cont: Vec<Vec<f64>>,
    /// `myPntCont` (`hxx:251`).
    pnt_cont: Vec<PointConstraint>,
    /// `mySurfInit` (`hxx:252`).
    surf_init: Option<Arc<dyn Surface>>,
    /// `myPlanarSurfInit` (`hxx:253`).
    planar_surf_init: Option<Arc<dyn Surface>>,
    /// `myGeomPlateSurface` (`hxx:254`).
    geom_plate_surface: Option<Arc<PlateSurface>>,
    /// `myPlate` (`hxx:255`).
    plate: Plate,
    /// `myPrevPlate` (`hxx:256`).
    prev_plate: Plate,
    /// `myAnisotropie` (`hxx:257`).
    anisotropie: bool,
    /// `mySense` (`hxx:258`).
    sense: Vec<i32>,
    /// `myDegree` (`hxx:259`).
    degree: i32,
    /// `myInitOrder` (`hxx:260`).
    init_order: Vec<i32>,
    /// `myG0Error` (`hxx:261`).
    g0_error: f64,
    /// `myG1Error` (`hxx:262`).
    g1_error: f64,
    /// `myG2Error` (`hxx:263`).
    g2_error: f64,
    /// `myNbPtsOnCur` (`hxx:264`).
    nb_pts_on_cur: i32,
    /// `mySurfInitIsGive` (`hxx:265`).
    surf_init_is_give: bool,
    /// `myNbIter` (`hxx:266`).
    nb_iter: i32,
    /// `myTol2d` (`hxx:267`). `myProj` is not stored, see the module note.
    tol2d: f64,
    /// `myTol3d` (`hxx:268`).
    tol3d: f64,
    /// `myTolAng` (`hxx:269`).
    tol_ang: f64,
    /// `myTolU` (`hxx:270`).
    tol_u: f64,
    /// `myTolV` (`hxx:271`).
    tol_v: f64,
    /// `myNbBounds` (`hxx:272`).
    nb_bounds: i32,
    /// `myIsLinear` (`hxx:273`).
    is_linear: bool,
    /// `myFree` (`hxx:274`).
    free: bool,
}

impl BuildPlateSurface {
    /// `GeomPlate_BuildPlateSurface(NPoints, TabCurve, Tang, Degree, NbIter,
    /// Tol2d, Tol3d, TolAng, TolCurv, Anisotropie)` (`cxx:80-142`): the old
    /// constructor that takes an array of curves.
    ///
    /// `n_points`, `tab_curve` and `tang` are the `1..=NTCurve` contents of the
    /// three `NCollection_HArray1` arguments.
    pub fn from_curves(
        n_points: &[i32],
        tab_curve: &[Arc<dyn crate::curve::Curve>],
        tang: &[i32],
        degree: i32,
        nb_iter: i32,
        tol2d: f64,
        tol3d: f64,
        tol_ang: f64,
        anisotropie: bool,
    ) -> Result<Self, String> {
        let mut this = Self {
            lin_cont: Vec::new(),
            par_cont: Vec::new(),
            plate_cont: Vec::new(),
            pnt_cont: Vec::new(),
            surf_init: None,
            planar_surf_init: None,
            geom_plate_surface: None,
            plate: Plate::default(),
            prev_plate: Plate::default(),
            anisotropie,
            sense: Vec::new(),
            degree,
            init_order: Vec::new(),
            g0_error: 0.0,
            g1_error: 0.0,
            g2_error: 0.0,
            nb_pts_on_cur: 0,
            surf_init_is_give: false,
            nb_iter,
            tol2d,
            tol3d,
            tol_ang,
            tol_u: 0.0,
            tol_v: 0.0,
            nb_bounds: 0,
            is_linear: true,
            free: false,
        };
        let nt_curve = tab_curve.len();
        // `cxx:105-108`.
        if this.nb_iter < 1 {
            return Err("GeomPlate :  Number of iteration must be >= 1".to_string());
        }
        // `cxx:109-112`.
        if nt_curve == 0 {
            return Err("GeomPlate : the bounds Array is null".to_string());
        }
        // `cxx:113-116`.
        if tang.is_empty() {
            return Err("GeomPlate : the constraints Array is null".to_string());
        }
        // `cxx:117-126`.
        let mut nbp = 0;
        for i in 0..nt_curve {
            nbp += n_points[i];
        }
        if nbp == 0 {
            return Err(
                "GeomPlate : the resolution is impossible if the number of constraints points is 0"
                    .to_string(),
            );
        }
        // `cxx:127-130`.
        if this.degree < 2 {
            return Err("GeomPlate ; the degree resolution must be upper of 2".to_string());
        }
        // `cxx:132-138`: the old constructor fills `myLinCont` through
        // `GeomPlate_CurveConstraint(Boundary, Order, NPt)` (`hxx:139-142` of
        // `GeomPlate_CurveConstraint.hxx`: TolDist 1.e-4, TolAng 1.e-2,
        // TolCurv 0.1).
        for i in 0..nt_curve {
            let cont = CurveConstraint::from_boundary(
                BoundaryCurve::Curve3d(tab_curve[i].clone()),
                tang[i],
                n_points[i],
                1.0e-4,
                1.0e-2,
                0.1,
            )?;
            this.lin_cont.push(cont);
        }
        Ok(this)
    }

    /// `GeomPlate_BuildPlateSurface(Surf, Degree, NbPtsOnCur, NbIter, Tol2d,
    /// Tol3d, TolAng, TolCurv, Anisotropie)` (`cxx:148-182`).
    ///
    /// `TolCurv` is unused in OCCT (`cxx:157` takes it unnamed); the port drops
    /// the parameter for the same reason.
    pub fn from_surface(
        surf: Arc<dyn Surface>,
        degree: i32,
        nb_pts_on_cur: i32,
        nb_iter: i32,
        tol2d: f64,
        tol3d: f64,
        tol_ang: f64,
        anisotropie: bool,
    ) -> Result<Self, String> {
        let mut this = Self::blank(degree, nb_pts_on_cur, nb_iter, tol2d, tol3d, tol_ang, anisotropie);
        this.surf_init = Some(surf);
        // `cxx:167-170`.
        if this.nb_iter < 1 {
            return Err("GeomPlate :  Number of iteration must be >= 1".to_string());
        }
        // `cxx:171-174`.
        if this.degree < 2 {
            return Err("GeomPlate : the degree must be above 2".to_string());
        }
        // `cxx:179`.
        this.surf_init_is_give = true;
        Ok(this)
    }

    /// `GeomPlate_BuildPlateSurface(Degree, NbPtsOnCur, NbIter, Tol2d, Tol3d,
    /// TolAng, TolCurv, Anisotropie)` (`cxx:186-219`).
    pub fn new(
        degree: i32,
        nb_pts_on_cur: i32,
        nb_iter: i32,
        tol2d: f64,
        tol3d: f64,
        tol_ang: f64,
        anisotropie: bool,
    ) -> Result<Self, String> {
        let this = Self::blank(degree, nb_pts_on_cur, nb_iter, tol2d, tol3d, tol_ang, anisotropie);
        // `cxx:206-209`.
        if this.nb_iter < 1 {
            return Err("GeomPlate :  Number of iteration must be >= 1".to_string());
        }
        // `cxx:210-213`.
        if this.degree < 2 {
            return Err("GeomPlate : the degree resolution must be upper of 2".to_string());
        }
        Ok(this)
    }

    /// The field block shared by the two constructors that take
    /// `(Degree, NbPtsOnCur, ...)` (`cxx:196-218`, `:158-180`).
    fn blank(
        degree: i32,
        nb_pts_on_cur: i32,
        nb_iter: i32,
        tol2d: f64,
        tol3d: f64,
        tol_ang: f64,
        anisotropie: bool,
    ) -> Self {
        Self {
            lin_cont: Vec::new(),
            par_cont: Vec::new(),
            plate_cont: Vec::new(),
            pnt_cont: Vec::new(),
            surf_init: None,
            planar_surf_init: None,
            geom_plate_surface: None,
            plate: Plate::default(),
            prev_plate: Plate::default(),
            anisotropie,
            sense: Vec::new(),
            degree,
            init_order: Vec::new(),
            g0_error: 0.0,
            g1_error: 0.0,
            g2_error: 0.0,
            nb_pts_on_cur,
            surf_init_is_give: false,
            nb_iter,
            tol2d,
            tol3d,
            tol_ang,
            tol_u: 0.0,
            tol_v: 0.0,
            nb_bounds: 0,
            is_linear: true,
            free: false,
        }
    }

    /// `Init` (`cxx:388-394`): drops every constraint.
    pub fn init(&mut self) {
        self.lin_cont.clear();
        self.pnt_cont.clear();
    }

    /// `LoadInitSurface` (`cxx:400-404`).
    pub fn load_init_surface(&mut self, surf: Arc<dyn Surface>) {
        self.surf_init = Some(surf);
        self.surf_init_is_give = true;
    }

    /// `Add(const handle(GeomPlate_CurveConstraint)&)` (`cxx:414-417`).
    pub fn add_curve_constraint(&mut self, cont: CurveConstraint) {
        self.lin_cont.push(cont);
    }

    /// `SetNbBounds` (`cxx:419-422`).
    pub fn set_nb_bounds(&mut self, nb_bounds: i32) {
        self.nb_bounds = nb_bounds;
    }

    /// `Add(const handle(GeomPlate_PointConstraint)&)` (`cxx:429-432`).
    pub fn add_point_constraint(&mut self, cont: PointConstraint) {
        self.pnt_cont.push(cont);
    }

    /// `IsDone` (`cxx:1167-1170`).
    pub fn is_done(&self) -> bool {
        self.plate.is_done()
    }

    /// `Surface` (`cxx:1175-1178`).
    pub fn surface(&self) -> Option<&Arc<PlateSurface>> {
        self.geom_plate_surface.as_ref()
    }

    /// `SurfInit` (`cxx:1184-1187`).
    pub fn surf_init(&self) -> Option<&Arc<dyn Surface>> {
        self.surf_init.as_ref()
    }

    /// `Sense` (`cxx:1192-1201`): the orientations of the curves in the order
    /// [`Self::curves2d`] returns them.
    pub fn sense(&self) -> Result<Vec<i32>, String> {
        let mut sens = Vec::with_capacity(self.lin_cont.len());
        for _ in 0..self.lin_cont.len() {
            sens.push(0);
        }
        for i in 0..self.lin_cont.len() {
            // `mySense->Value(myInitOrder->Value(i))` (`cxx:1199`).
            let at = self.order_at(i)?;
            sens[i] = self.sense_at(at)?;
        }
        Ok(sens)
    }

    /// `Curves2d` (`cxx:1207-1218`).
    pub fn curves2d(&self) -> Result<Vec<Option<Arc<dyn occt_geom2d::Curve2d>>>, String> {
        let mut out = Vec::with_capacity(self.lin_cont.len());
        for i in 0..self.lin_cont.len() {
            let at = self.order_at(i)?;
            out.push(self.lin_cont[at].curve2d_on_surf()?);
        }
        Ok(out)
    }

    /// `Order` (`cxx:1223-1233`).
    pub fn order(&self) -> Result<Vec<i32>, String> {
        if self.init_order.len() < self.lin_cont.len() {
            return Err("GeomPlate_BuildPlateSurface::Order: no initial order".to_string());
        }
        let mut result = vec![0; self.lin_cont.len()];
        for i in 0..self.lin_cont.len() {
            // `result->SetValue(myInitOrder->Value(i), i)` (`cxx:1229`).
            let at = self.order_at(i)?;
            result[at] = (i + 1) as i32;
        }
        Ok(result)
    }

    /// `myInitOrder->Value(i)` with the `1..=N` to `0..N` shift and the
    /// null-handle guard.
    fn order_at(&self, i: usize) -> Result<usize, String> {
        let value = self
            .init_order
            .get(i)
            .ok_or_else(|| "GeomPlate_BuildPlateSurface: no initial order".to_string())?;
        if *value < 1 {
            return Err("GeomPlate_BuildPlateSurface: initial order out of range".to_string());
        }
        Ok((*value - 1) as usize)
    }

    /// `mySense->Value(j)` with the `1..=N` to `0..N` shift.
    fn sense_at(&self, j: usize) -> Result<i32, String> {
        self.sense
            .get(j)
            .copied()
            .ok_or_else(|| "GeomPlate_BuildPlateSurface: no sense".to_string())
    }

    /// `G0Error` (`cxx:1237-1240`).
    pub fn g0_error(&self) -> f64 {
        self.g0_error
    }

    /// `G1Error` (`cxx:1245-1248`).
    pub fn g1_error(&self) -> f64 {
        self.g1_error
    }

    /// `G2Error` (`cxx:1253-1256`).
    pub fn g2_error(&self) -> f64 {
        self.g2_error
    }

    /// `G0Error(Index)` (`cxx:1262-1278`).
    pub fn g0_error_of(&mut self, index: usize) -> Result<f64, String> {
        let (tdistance, tangle, tcurvature) = self.ecart_arrays(index)?;
        let mut max_distance = 0.0;
        for i in 0..self.nb_pts_on_cur.max(0) as usize {
            if tdistance[i] > max_distance {
                max_distance = tdistance[i];
            }
        }
        let _ = (tangle, tcurvature);
        Ok(max_distance)
    }

    /// `G1Error(Index)` (`cxx:1283-1300`).
    pub fn g1_error_of(&mut self, index: usize) -> Result<f64, String> {
        let (_, tangle, _) = self.ecart_arrays(index)?;
        let mut max_angle = 0.0;
        for i in 0..self.nb_pts_on_cur.max(0) as usize {
            if tangle[i] > max_angle {
                max_angle = tangle[i];
            }
        }
        Ok(max_angle)
    }

    /// `G2Error(Index)` (`cxx:1305-1322`).
    pub fn g2_error_of(&mut self, index: usize) -> Result<f64, String> {
        let (_, _, tcurvature) = self.ecart_arrays(index)?;
        let mut max_curvature = 0.0;
        for i in 0..self.nb_pts_on_cur.max(0) as usize {
            if tcurvature[i] > max_curvature {
                max_curvature = tcurvature[i];
            }
        }
        Ok(max_curvature)
    }

    /// The three `NCollection_HArray1<double>(1, myNbPtsOnCur)` work arrays of
    /// `G0Error(Index)` / `G1Error(Index)` / `G2Error(Index)`
    /// (`cxx:1262-1266`, `:1285-1289`, `:1308-1312`).
    fn ecart_arrays(
        &mut self,
        index: usize,
    ) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), String> {
        let n = self.nb_pts_on_cur.max(0) as usize;
        let mut tdistance = vec![0.0; n];
        let mut tangle = vec![0.0; n];
        let mut tcurvature = vec![0.0; n];
        self.ecart_contraintes_mil(index, &mut tdistance, &mut tangle, &mut tcurvature)?;
        Ok((tdistance, tangle, tcurvature))
    }

    /// `CurveConstraint(order)` (`cxx:1325-1328`).
    ///
    /// OCCT returns the handle at the `1..=N` index; the port returns a borrow
    /// and treats `order == 0` (a null handle in OCCT) as out of range.
    pub fn curve_constraint(&self, order: usize) -> Result<&CurveConstraint, String> {
        if order == 0 {
            return Err("GeomPlate_BuildPlateSurface: null curve constraint".to_string());
        }
        self.lin_cont
            .get(order - 1)
            .ok_or_else(|| "GeomPlate_BuildPlateSurface: curve constraint out of range".to_string())
    }

    /// `PointConstraint(order)` (`cxx:1330-1333`).
    pub fn point_constraint(&self, order: usize) -> Result<&PointConstraint, String> {
        if order == 0 {
            return Err("GeomPlate_BuildPlateSurface: null point constraint".to_string());
        }
        self.pnt_cont
            .get(order - 1)
            .ok_or_else(|| "GeomPlate_BuildPlateSurface: point constraint out of range".to_string())
    }

    /// `myLinCont` by index, 0-based, for the submodules.
    fn lin(&self, i: usize) -> &CurveConstraint {
        &self.lin_cont[i]
    }

    /// `myLinCont->ChangeValue(i)` by index, 0-based.
    fn lin_mut(&mut self, i: usize) -> &mut CurveConstraint {
        &mut self.lin_cont[i]
    }

    /// `myPntCont` by index, 0-based.
    fn pnt(&self, i: usize) -> &PointConstraint {
        &self.pnt_cont[i]
    }

    /// `myPntCont->ChangeValue(i)` by index, 0-based.
    fn pnt_mut(&mut self, i: usize) -> &mut PointConstraint {
        &mut self.pnt_cont[i]
    }

    /// `myParCont->Value(i)` by index, 0-based.
    fn par(&self, i: usize) -> &Vec<f64> {
        &self.par_cont[i]
    }

    /// `myPlateCont->Value(i)` by index, 0-based.
    fn plate_at(&self, i: usize) -> &Vec<f64> {
        &self.plate_cont[i]
    }

    /// `ComputeAnisotropie` (`cxx:2783-2796`): both arms return `1.0` in OCCT.
    pub fn compute_anisotropie(&self) -> f64 {
        if self.anisotropie {
            // `cxx:2787-2789`: "Temporary".
            1.0
        } else {
            1.0
        }
    }

    /// `IsOrderG1` (`cxx:2798-2808`).
    pub fn is_order_g1(&self) -> bool {
        for i in 0..self.lin_cont.len() {
            if self.lin(i).order() < 1 {
                return false;
            }
        }
        true
    }

    /// `TrierTab` (`cxx:239-249`): `tab[tab[i] - 1] = i + 1` for every `i`.
    fn trier_tab(tab: &mut [i32]) {
        let nb = tab.len();
        let mut tri = vec![0; nb];
        for i in 0..nb {
            let at = (tab[i] - 1) as usize;
            tri[at] = (i + 1) as i32;
        }
        tab.copy_from_slice(&tri);
    }

    /// `myParCont->Value(c).Length()` with the null-handle guard.
    fn par_len(&self, c: usize) -> usize {
        self.par_cont.get(c).map(|v| v.len()).unwrap_or(0)
    }

    /// `myPlateCont->Value(c).Length()` with the null-handle guard.
    fn plate_len(&self, c: usize) -> usize {
        self.plate_cont.get(c).map(|v| v.len()).unwrap_or(0)
    }
}

/// The `gp_Pnt` a constraint evaluates to, used by the four `D0` call sites of
/// `Discretise` / `LoadCurve` / `LoadPoint` (`cxx:2285-2289`, `:2327-2332`,
/// `:2425`, `:2534`).
fn constraint_pnt(p: &GpPnt) -> GpPnt {
    *p
}
