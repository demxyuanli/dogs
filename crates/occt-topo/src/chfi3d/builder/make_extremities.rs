//! `ChFi3d_MakeExtremities`, `ChFi3d_SingularExtremity` and `IsFree`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx`
//! (OCCT 8.0.0):
//!
//! - `ChFi3d_SingularExtremity` (`cxx:1635-1743`),
//! - `IsFree` (`cxx:1744-1765`),
//! - `ChFi3d_MakeExtremities` (`cxx:1766-2046`).
//!
//! Translation notes:
//! - `DStr` is the builder's [`HDataStructure`] and `EFMap` its
//!   [`ChFiDSMap`], so the two statics become methods,
//! - `Stripe->ChangeFirstPCurve()` is written directly in OCCT; the port keeps
//!   the computed p-curve local and stores it on the stripe once the DS is no
//!   longer borrowed,
//! - `ChFi3d_CoupeParPlan` fills `C3d` / `pcurve` / `tolreached` / `Pardeb` /
//!   `Parfin` through references; the port answers them in
//!   [`CoupeParPlan`](crate::chfi3d::builder_2_helpers::CoupeParPlan),
//! - `sp->SetFirstStatus(ChFiDS_FreeBoundary)` / `SetLastStatus` run on the
//!   shared-mutable spine handle, as everywhere else in the crate.

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt2d;
use occt_geom::Surface;

use crate::chfi_ds::{CommonPoint, State, StripeHandle, SurfData};
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_bound_surf, ch_fi3d_compute_arete, ch_fi3d_compute_pcurv,
    ch_fi3d_enlarge_box_edge, ch_fi3d_enlarge_box_surf_data, ch_fi3d_index_of_surf_data,
    ch_fi3d_index_point_in_ds, ch_fi3d_set_point_tolerance, ComputeArete,
};
use crate::chfi3d::builder_2_helpers::{ch_fi3d_coupe_par_plan, CoupeParPlan};
use crate::shape::TopoShape;
use crate::topopebrep::ds::Curve as DsCurve;

use super::ChFi3dBuilder;

/// `IsFree(E, EFMap)` (`cxx:1744-1765`): every face of `EFMap(E)` is the same
/// face, and `E` is a key of the map.
pub fn is_free(e: &TopoShape, ef_map: &crate::chfi_ds::ChFiDSMap) -> bool {
    // `if (!EFMap.Contains(E)) return false;` (`cxx:1747-1750`).
    if !ef_map.contains(e) {
        return false;
    }
    // `for (It.Initialize(EFMap(E)); It.More(); It.Next())` (`cxx:1751-1762`).
    let mut fref: Option<&TopoShape> = None;
    for value in ef_map.get_from_key(e) {
        match fref {
            None => fref = Some(value),
            Some(first) => {
                if !crate::topo_tools::is_same(first, value) {
                    return false;
                }
            }
        }
    }
    true
}

/// `Surface(Fd)` — `DStr.Surface(Fd->Surf()).Surface()`.
fn surf_data_surface(
    ds: &crate::topopebrep::ds::DataStructure,
    fd: &SurfData,
) -> Result<Arc<dyn Surface>, String> {
    ds.surface(fd.surf())
        .surface()
        .cloned()
        .ok_or_else(|| "ChFi3d_MakeExtremities: DStr.Surface(Surf()).Surface() is null".to_string())
}

impl ChFi3dBuilder {
    /// `ChFi3d_SingularExtremity(stripe, DStr, Vtx, tol3d, tol2d)`
    /// (`cxx:1635-1743`).
    pub fn singular_extremity(
        &mut self,
        stripe: &StripeHandle,
        vtx: &TopoShape,
        tol3d: f64,
        tol2d: f64,
    ) -> Result<(), String> {
        // `if (stripe->Spine()->IsPeriodic()) { isfirst = true; Fd =
        // ...Sequence().First(); } else { num = ChFi3d_IndexOfSurfData(Vtx,
        // stripe, sens); Fd = ...Value(num); isfirst = (sens == 1); }`
        // (`cxx:1652-1666`).
        let (fd, isfirst) = {
            let guard = stripe
                .read()
                .map_err(|_| "ChFi3d_SingularExtremity: poisoned stripe".to_string())?;
            let spine = guard
                .spine()
                .cloned()
                .ok_or_else(|| "ChFi3d_SingularExtremity: the stripe has no spine".to_string())?;
            let periodic = spine
                .read()
                .map_err(|_| "ChFi3d_SingularExtremity: poisoned spine".to_string())?
                .is_periodic();
            if periodic {
                let fd = guard
                    .set_of_surf_data()
                    .first()
                    .cloned()
                    .ok_or_else(|| {
                        "ChFi3d_SingularExtremity: empty SurfData sequence".to_string()
                    })?;
                (fd, true)
            } else {
                drop(guard);
                let mut sens = 0;
                let num = ch_fi3d_index_of_surf_data(vtx, stripe, &mut sens)?;
                let guard = stripe
                    .read()
                    .map_err(|_| "ChFi3d_SingularExtremity: poisoned stripe".to_string())?;
                let fd = usize::try_from(num)
                    .ok()
                    .and_then(|n| n.checked_sub(1))
                    .and_then(|k| guard.set_of_surf_data().get(k))
                    .cloned()
                    .ok_or_else(|| {
                        "ChFi3d_SingularExtremity: SurfData index out of range".to_string()
                    })?;
                (fd, sens == 1)
            }
        };

        // `const ChFiDS_CommonPoint& CV1 = Fd->Vertex(isfirst, 1); CV2 =
        // Fd->Vertex(isfirst, 2);` (`cxx:1668-1669`).
        let cv1 = fd.vertex(isfirst, 1);
        let cv2 = fd.vertex(isfirst, 2);
        // `if (CV1.Point().IsEqual(CV2.Point(), 0))` (`cxx:1670-1672`).
        if !cv1.point().is_equal(cv2.point()) {
            return Ok(());
        }

        let is_periodic = {
            let guard = stripe
                .read()
                .map_err(|_| "ChFi3d_SingularExtremity: poisoned stripe".to_string())?;
            let spine = guard
                .spine()
                .cloned()
                .ok_or_else(|| "ChFi3d_SingularExtremity: the stripe has no spine".to_string())?;
            let sid = spine
                .read()
                .map_err(|_| "ChFi3d_SingularExtremity: poisoned spine".to_string())?;
            let periodic = sid.is_periodic();
            periodic
        };

        // `Ivtx = ChFi3d_IndexPointInDS(CV1, DStr);` (`cxx:1671`).
        let cv1 = cv1.clone();
        let cv2 = cv2.clone();
        let ivtx = ch_fi3d_index_point_in_ds(&cv1, self.my_ds.change_ds());

        // The first pass, on `Fd` (`cxx:1673-1704`).
        self.singular_pass(stripe, &fd, &cv1, &cv2, ivtx, isfirst, tol3d, tol2d)?;

        // `if (stripe->Spine()->IsPeriodic())` (`cxx:1705-1738`): the operation
        // is renewed on the last SurfData with its last parameters.
        if is_periodic {
            let fd_last = {
                let guard = stripe
                    .read()
                    .map_err(|_| "ChFi3d_SingularExtremity: poisoned stripe".to_string())?;
                guard.set_of_surf_data().last().cloned().ok_or_else(|| {
                    "ChFi3d_SingularExtremity: empty SurfData sequence".to_string()
                })?
            };
            self.singular_pass(stripe, &fd_last, &cv1, &cv2, ivtx, false, tol3d, tol2d)?;
        }
        Ok(())
    }

    /// One `ChFi3d_ComputeArete` + store block of
    /// [`Self::singular_extremity`] (`cxx:1673-1704`, `cxx:1714-1737`).
    #[allow(clippy::too_many_arguments)]
    fn singular_pass(
        &mut self,
        stripe: &StripeHandle,
        fd: &Arc<SurfData>,
        cv1: &CommonPoint,
        cv2: &CommonPoint,
        ivtx: i32,
        isfirst: bool,
        tol3d: f64,
        tol2d: f64,
    ) -> Result<(), String> {
        // `VOnS1 = Fd->InterferenceOnS1().PCurveOnSurf()->Value(Fd->
        // InterferenceOnS1().FirstParameter());` — `Last` when `!isfirst`
        // (`cxx:1675-1681`).
        let fi1 = fd.interference_on_s1();
        let fi2 = fd.interference_on_s2();
        let pc1 = fi1.p_curve_on_surf().ok_or_else(|| {
            "ChFi3d_SingularExtremity: InterferenceOnS1 has no p-curve".to_string()
        })?;
        let pc2 = fi2.p_curve_on_surf().ok_or_else(|| {
            "ChFi3d_SingularExtremity: InterferenceOnS2 has no p-curve".to_string()
        })?;
        let (p1, p2) = if isfirst {
            (fi1.first_parameter(), fi2.first_parameter())
        } else {
            (fi1.last_parameter(), fi2.last_parameter())
        };
        let von_s1: GpPnt2d = pc1.value(p1);
        let von_s2: GpPnt2d = pc2.value(p2);

        // `ChFi3d_ComputeArete(CV1, VOnS1, CV2, VOnS2,
        // DStr.Surface(Fd->Surf()).Surface(), C3d, PCurv, ...)` (`cxx:1683-1695`).
        let surface = surf_data_surface(self.my_ds.ds(), fd)?;
        let mut arete = ComputeArete::default();
        ch_fi3d_compute_arete(cv1, &von_s1, cv2, &von_s2, &surface, &mut arete, tol3d, tol2d, 0)?;

        // `Crv = TopOpeBRepDS_Curve(C3d, tolreached); Icurv =
        // DStr.AddCurve(Crv);` (`cxx:1696-1697`).
        let mut crv = DsCurve::with_curve(arete.c3d.clone(), arete.tolreached, false);
        let icurv = self.my_ds.change_ds().add_curve(&mut crv);

        // `stripe->SetCurve(Icurv, isfirst); SetParameters(isfirst, Pardeb,
        // Parfin); ChangePCurve(isfirst) = PCurv; SetIndexPoint(Ivtx,
        // isfirst, 1); SetIndexPoint(Ivtx, isfirst, 2);` (`cxx:1699-1703`).
        let mut guard = stripe
            .write()
            .map_err(|_| "ChFi3d_SingularExtremity: poisoned stripe".to_string())?;
        guard.set_curve(icurv, isfirst);
        guard.set_parameters(isfirst, arete.pardeb, arete.parfin);
        if isfirst {
            *guard.change_first_p_curve() = arete.p_curve;
        } else {
            *guard.change_last_p_curve() = arete.p_curve;
        }
        guard.set_index_point(ivtx, isfirst, 1);
        guard.set_index_point(ivtx, isfirst, 2);
        Ok(())
    }

    /// `ChFi3d_MakeExtremities(Stripe, DStr, EFMap, tol3d, tol2d)`
    /// (`cxx:1766-2046`): the 3D curve and p-curve of each extremity, in the
    /// periodic and the free-boundary cases.
    pub fn make_extremities(
        &mut self,
        stripe: &StripeHandle,
        tol3d: f64,
        tol2d: f64,
    ) -> Result<(), String> {
        let spine = {
            let guard = stripe
                .read()
                .map_err(|_| "ChFi3d_MakeExtremities: poisoned stripe".to_string())?;
            guard
                .spine()
                .cloned()
                .ok_or_else(|| "ChFi3d_MakeExtremities: the stripe has no spine".to_string())?
        };
        let periodic = spine
            .read()
            .map_err(|_| "ChFi3d_MakeExtremities: poisoned spine".to_string())?
            .is_periodic();

        if periodic {
            let (sdf, sdl) = {
                let guard = stripe
                    .read()
                    .map_err(|_| "ChFi3d_MakeExtremities: poisoned stripe".to_string())?;
                let sds = guard.set_of_surf_data();
                (
                    sds.first().cloned(),
                    sds.last().cloned().or_else(|| sds.first().cloned()),
                )
            };
            let (Some(sdf), Some(sdl)) = (sdf, sdl) else {
                return Ok(());
            };
            // `CV1 = SDF->VertexFirstOnS1(); CV2 = SDF->VertexFirstOnS2();`
            // (`cxx:1780-1781`).
            let cv1 = sdf.vertex(true, 1).clone();
            let cv2 = sdf.vertex(true, 2).clone();
            if cv1.point().is_equal(cv2.point()) {
                // `if (CV1.IsVertex()) ChFi3d_SingularExtremity(...)`
                // (`cxx:1845-1849`).
                if let Some(v) = cv1.vertex() {
                    let v = v.clone();
                    return self.singular_extremity(stripe, &v, tol3d, tol2d);
                }
                return Ok(());
            }

            // `cxx:1783-1807`: the arete on the first SurfData.
            let fi1 = sdf.interference_on_s1();
            let fi2 = sdf.interference_on_s2();
            let uv1 = fi1
                .p_curve_on_surf()
                .ok_or_else(|| "ChFi3d_MakeExtremities: SDF/1 has no p-curve".to_string())?
                .value(fi1.first_parameter());
            let uv2 = fi2
                .p_curve_on_surf()
                .ok_or_else(|| "ChFi3d_MakeExtremities: SDF/2 has no p-curve".to_string())?
                .value(fi2.first_parameter());
            let surface = surf_data_surface(self.my_ds.ds(), &sdf)?;
            let mut arete = ComputeArete::default();
            ch_fi3d_compute_arete(
                &cv1, &uv1, &cv2, &uv2, &surface, &mut arete, tol3d, tol2d, 0,
            )?;
            let c3d = arete.c3d.clone();
            let mut crv = DsCurve::with_curve(c3d.clone(), arete.tolreached, false);
            let icurv = self.my_ds.change_ds().add_curve(&mut crv);
            let i_first_s1 = ch_fi3d_index_point_in_ds(&cv1, self.my_ds.change_ds());
            let i_first_s2 = ch_fi3d_index_point_in_ds(&cv2, self.my_ds.change_ds());

            // `cxx:1799-1807`.
            {
                let mut guard = stripe
                    .write()
                    .map_err(|_| "ChFi3d_MakeExtremities: poisoned stripe".to_string())?;
                guard.set_curve(icurv, true);
                guard.set_parameters(true, arete.pardeb, arete.parfin);
                *guard.change_first_p_curve() = arete.p_curve.clone();
                guard.change_index_first_point_on_s1(i_first_s1);
                guard.change_index_first_point_on_s2(i_first_s2);
                guard.change_last_parameters(arete.pardeb, arete.parfin);
                let first_curve = guard.first_curve();
                guard.change_last_curve(first_curve);
                let first_s1 = guard.index_first_point_on_s1();
                let first_s2 = guard.index_first_point_on_s2();
                guard.change_index_last_point_on_s1(first_s1);
                guard.change_index_last_point_on_s2(first_s2);
            }

            // `cxx:1809-1821`: the p-curve on the last SurfData, with its own
            // last parameters.
            let fi1 = sdl.interference_on_s1();
            let fi2 = sdl.interference_on_s2();
            let uv1 = fi1
                .p_curve_on_surf()
                .ok_or_else(|| "ChFi3d_MakeExtremities: SDL/1 has no p-curve".to_string())?
                .value(fi1.last_parameter());
            let uv2 = fi2
                .p_curve_on_surf()
                .ok_or_else(|| "ChFi3d_MakeExtremities: SDL/2 has no p-curve".to_string())?
                .value(fi2.last_parameter());
            let sdl_surface = surf_data_surface(self.my_ds.ds(), &sdl)?;
            let c3d = c3d.ok_or_else(|| "ChFi3d_MakeExtremities: no 3D curve".to_string())?;
            let mut last_pcurve = None;
            let tol_pc = ch_fi3d_compute_pcurv(
                c3d.as_ref(),
                &uv1,
                &uv2,
                &mut last_pcurve,
                sdl_surface.as_ref(),
                arete.pardeb,
                arete.parfin,
                tol3d,
                false,
            )?;

            // `double oldtol = DStr.ChangeCurve(ICurv).Tolerance();
            // DStr.ChangeCurve(ICurv).Tolerance(max(oldtol, tolreached));`
            // (`cxx:1822-1824`).
            {
                let curve = self.my_ds.change_ds().change_curve(icurv);
                let oldtol = curve.tolerance();
                curve.set_tolerance(oldtol.max(tol_pc));
            }

            // `cxx:1825-1842`.
            let mut b1 = BndBox::new();
            let mut b2 = BndBox::new();
            if let (Some(arc), Some(param)) = (cv1.arc(), cv1.parameter_on_arc()) {
                ch_fi3d_enlarge_box_edge(arc, self.my_ef_map.get_from_key(arc), param, &mut b1);
            }
            if let (Some(arc), Some(param)) = (cv2.arc(), cv2.parameter_on_arc()) {
                ch_fi3d_enlarge_box_edge(arc, self.my_ef_map.get_from_key(arc), param, &mut b2);
            }
            {
                let guard = stripe
                    .read()
                    .map_err(|_| "ChFi3d_MakeExtremities: poisoned stripe".to_string())?;
                ch_fi3d_enlarge_box_surf_data(
                    self.my_ds.ds(),
                    Some(&guard),
                    &sdf,
                    &mut b1,
                    &mut b2,
                    true,
                )?;
                ch_fi3d_enlarge_box_surf_data(
                    self.my_ds.ds(),
                    Some(&guard),
                    &sdl,
                    &mut b1,
                    &mut b2,
                    false,
                )?;
            }
            let (i_first_s1, i_first_s2) = {
                let guard = stripe
                    .read()
                    .map_err(|_| "ChFi3d_MakeExtremities: poisoned stripe".to_string())?;
                (
                    guard.index_first_point_on_s1(),
                    guard.index_first_point_on_s2(),
                )
            };
            if !cv1.is_vertex() {
                ch_fi3d_set_point_tolerance(self.my_ds.change_ds(), &b1, i_first_s1);
            }
            if !cv2.is_vertex() {
                ch_fi3d_set_point_tolerance(self.my_ds.change_ds(), &b2, i_first_s2);
            }
            // `stripe->ChangeLastPCurve()` was not written yet: store it.
            if last_pcurve.is_some() {
                let mut guard = stripe
                    .write()
                    .map_err(|_| "ChFi3d_MakeExtremities: poisoned stripe".to_string())?;
                *guard.change_last_p_curve() = last_pcurve;
            }
            return Ok(());
        }

        // `cxx:1860-1951`: the start extremity.
        let sddeb = {
            let guard = stripe
                .read()
                .map_err(|_| "ChFi3d_MakeExtremities: poisoned stripe".to_string())?;
            guard.set_of_surf_data().first().cloned()
        };
        let Some(sddeb) = sddeb else {
            return Ok(());
        };
        let cpdeb1 = sddeb.vertex(true, 1).clone();
        let cpdeb2 = sddeb.vertex(true, 2).clone();
        let mut freedeb = spine
            .read()
            .map_err(|_| "ChFi3d_MakeExtremities: poisoned spine".to_string())?
            .first_status()
            == State::FreeBoundary;
        if !freedeb && cpdeb1.is_on_arc() && cpdeb2.is_on_arc() {
            // `freedeb = IsFree(cpdeb1.Arc(), EFMap) && IsFree(cpdeb2.Arc(),
            // EFMap);` (`cxx:1865-1868`).
            let (Some(arc1), Some(arc2)) = (cpdeb1.arc(), cpdeb2.arc()) else {
                return Ok(());
            };
            freedeb = is_free(arc1, &self.my_ef_map) && is_free(arc2, &self.my_ef_map);
        }
        if freedeb {
            self.free_extremity(stripe, &spine, &sddeb, &cpdeb1, &cpdeb2, true, tol3d, tol2d)?;
        }

        // `cxx:1952-2042`: the end extremity.
        let sdlin = {
            let guard = stripe
                .read()
                .map_err(|_| "ChFi3d_MakeExtremities: poisoned stripe".to_string())?;
            guard.set_of_surf_data().last().cloned()
        };
        let Some(sdlin) = sdlin else {
            return Ok(());
        };
        let cpfin1 = sdlin.vertex(false, 1).clone();
        let cpfin2 = sdlin.vertex(false, 2).clone();
        let mut freefin = spine
            .read()
            .map_err(|_| "ChFi3d_MakeExtremities: poisoned spine".to_string())?
            .last_status()
            == State::FreeBoundary;
        if !freefin && cpfin1.is_on_arc() && cpfin2.is_on_arc() {
            let (Some(arc1), Some(arc2)) = (cpfin1.arc(), cpfin2.arc()) else {
                return Ok(());
            };
            freefin = is_free(arc1, &self.my_ef_map) && is_free(arc2, &self.my_ef_map);
        }
        if freefin {
            self.free_extremity(stripe, &spine, &sdlin, &cpfin1, &cpfin2, false, tol3d, tol2d)?;
        }
        Ok(())
    }

    /// One free-boundary extremity (`cxx:1869-1950` for `isfirst`, the mirror
    /// `cxx:1961-2041` otherwise).
    #[allow(clippy::too_many_arguments)]
    fn free_extremity(
        &mut self,
        stripe: &StripeHandle,
        spine: &crate::chfi_ds::SpineHandle,
        sd: &Arc<SurfData>,
        cp1: &CommonPoint,
        cp2: &CommonPoint,
        isfirst: bool,
        tol3d: f64,
        tol2d: f64,
    ) -> Result<(), String> {
        // `sp->SetFirstStatus(ChFiDS_FreeBoundary);` (`cxx:1871`,
        // `cxx:1962`).
        spine
            .write()
            .map_err(|_| "ChFi3d_MakeExtremities: poisoned spine".to_string())?
            .set_status(State::FreeBoundary, isfirst);

        if !cp1.point().is_equal(cp2.point()) {
            let fi1 = sd.interference_on_s1();
            let fi2 = sd.interference_on_s2();
            let param1 = if isfirst {
                fi1.first_parameter()
            } else {
                fi1.last_parameter()
            };
            let param2 = if isfirst {
                fi2.first_parameter()
            } else {
                fi2.last_parameter()
            };
            let uv1 = fi1
                .p_curve_on_surf()
                .ok_or_else(|| "ChFi3d_MakeExtremities: InterferenceOnS1 has no p-curve".to_string())?
                .value(param1);
            let uv2 = fi2
                .p_curve_on_surf()
                .ok_or_else(|| "ChFi3d_MakeExtremities: InterferenceOnS2 has no p-curve".to_string())?
                .value(param2);

            // `HConge = ChFi3d_BoundSurf(DStr, SD, 1, 2);` (`cxx:1886`).
            let hconge = ch_fi3d_bound_surf(self.my_ds.ds(), sd, 1, 2)?;
            // `ChFi3d_CoupeParPlan(...)` (`cxx:1887-1898`).
            let mut coupe = CoupeParPlan::default();
            ch_fi3d_coupe_par_plan(
                cp1, cp2, &hconge, &uv1, &uv2, tol3d, tol2d, &mut coupe,
            );
            let (c3d, pcurve, tolreached, pardeb, parfin) = if coupe.plane {
                (
                    coupe.c3d.clone(),
                    coupe.pcurve.clone(),
                    coupe.tolreached,
                    coupe.pardeb,
                    coupe.parfin,
                )
            } else {
                // `ChFi3d_ComputeArete(...)` (`cxx:1898-1915`).
                let surface = surf_data_surface(self.my_ds.ds(), sd)?;
                let mut arete = ComputeArete::default();
                ch_fi3d_compute_arete(
                    cp1, &uv1, cp2, &uv2, &surface, &mut arete, tol3d, tol2d, 0,
                )?;
                (
                    arete.c3d,
                    arete.p_curve,
                    arete.tolreached,
                    arete.pardeb,
                    arete.parfin,
                )
            };

            // `Crv = TopOpeBRepDS_Curve(C3d, tolreached);
            // Stripe->ChangeFirstCurve(DStr.AddCurve(Crv));` (`cxx:1916-1918`).
            let mut crv = DsCurve::with_curve(c3d, tolreached, false);
            let icurv = self.my_ds.change_ds().add_curve(&mut crv);
            let i1 = ch_fi3d_index_point_in_ds(cp1, self.my_ds.change_ds());
            let i2 = ch_fi3d_index_point_in_ds(cp2, self.my_ds.change_ds());
            {
                let mut guard = stripe
                    .write()
                    .map_err(|_| "ChFi3d_MakeExtremities: poisoned stripe".to_string())?;
                if isfirst {
                    guard.change_first_curve(icurv);
                    guard.change_first_parameters(pardeb, parfin);
                    guard.change_index_first_point_on_s1(i1);
                    guard.change_index_first_point_on_s2(i2);
                    *guard.change_first_p_curve() = pcurve;
                } else {
                    guard.change_last_curve(icurv);
                    guard.change_last_parameters(pardeb, parfin);
                    guard.change_index_last_point_on_s1(i1);
                    guard.change_index_last_point_on_s2(i2);
                    *guard.change_last_p_curve() = pcurve;
                }
            }

            // `cxx:1919-1937`.
            let mut b1 = BndBox::new();
            let mut b2 = BndBox::new();
            if let (Some(arc), Some(param)) = (cp1.arc(), cp1.parameter_on_arc()) {
                ch_fi3d_enlarge_box_edge(arc, self.my_ef_map.get_from_key(arc), param, &mut b1);
            }
            if let (Some(arc), Some(param)) = (cp2.arc(), cp2.parameter_on_arc()) {
                ch_fi3d_enlarge_box_edge(arc, self.my_ef_map.get_from_key(arc), param, &mut b2);
            }
            {
                let guard = stripe
                    .read()
                    .map_err(|_| "ChFi3d_MakeExtremities: poisoned stripe".to_string())?;
                ch_fi3d_enlarge_box_surf_data(
                    self.my_ds.ds(),
                    Some(&guard),
                    sd,
                    &mut b1,
                    &mut b2,
                    isfirst,
                )?;
            }
            let (idx1, idx2) = {
                let guard = stripe
                    .read()
                    .map_err(|_| "ChFi3d_MakeExtremities: poisoned stripe".to_string())?;
                if isfirst {
                    (
                        guard.index_first_point_on_s1(),
                        guard.index_first_point_on_s2(),
                    )
                } else {
                    (
                        guard.index_last_point_on_s1(),
                        guard.index_last_point_on_s2(),
                    )
                }
            };
            if !cp1.is_vertex() {
                ch_fi3d_set_point_tolerance(self.my_ds.change_ds(), &b1, idx1);
            }
            if !cp2.is_vertex() {
                ch_fi3d_set_point_tolerance(self.my_ds.change_ds(), &b2, idx2);
            }
            return Ok(());
        }

        // `else if (cp1.IsVertex()) ChFi3d_SingularExtremity(Stripe, DStr,
        // cp1.Vertex(), tol3d, tol2d);` (`cxx:1939-1943`, `:2030-2034`).
        if let Some(v) = cp1.vertex() {
            let v = v.clone();
            return self.singular_extremity(stripe, &v, tol3d, tol2d);
        }
        Ok(())
    }
}
