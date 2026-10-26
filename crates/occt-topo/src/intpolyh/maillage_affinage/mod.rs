//! `IntPolyh_MaillageAffinage` (`IntPolyh_MaillageAffinage.hxx:38-258`,
//! `.cxx:251-3997`): intersection of two surfaces through the intersection of
//! their triangulations.
//!
//! Layout of the port:
//! * this file - construction, sampling, the point/edge/triangle arrays, the
//!   common box and the deflection bookkeeping;
//! * `refine` - `GetInterferingTriangles` plus the deflection-driven
//!   refinement (`cxx:1119-1387`);
//! * `contact` - `TriContact` / `project6` / `StartingPointsResearch` /
//!   `NextStartingPointsResearch` / `TriangleCompare` (`cxx:1429-2048`,
//!   `cxx:3141-3185`);
//! * `edge_contact` - `TriangleEdgeContact` and
//!   `CalculPtsInterTriEdgeCoplanaires` (`cxx:2051-3138`);
//! * `chain` - `StartPointsChain` and its helpers (`cxx:3188-3781`);
//! * `degenerated` - `IsDegenerated` / `DegeneratedIndex` (`cxx:3863-3996`).

pub mod chain;
pub mod common_box;
pub mod contact;
pub mod degenerated;
pub mod edge_contact;
pub mod refine;

use occt_core::bnd::BndBox;
use occt_geom::Surface;

use super::arrays::{ArrayOfEdges, ArrayOfPoints, ArrayOfTriangles};
use super::edge_couple::IntPolyhCouple;
use super::point::IntPolyhPoint;
use super::tools;
use super::triangle::{middle_refinement_at, IntPolyhTriangle};

pub use degenerated::{degenerated_index, is_degenerated};

/// `MyTolerance` of `IntPolyh_MaillageAffinage.cxx:51`.
pub const MY_TOLERANCE: f64 = 10.0e-7;
/// `MyConfusionPrecision` of `IntPolyh_MaillageAffinage.cxx:52`.
pub const MY_CONFUSION_PRECISION: f64 = 10.0e-12;
/// `SquareMyConfusionPrecision` of `IntPolyh_MaillageAffinage.cxx:53`.
pub const SQUARE_MY_CONFUSION_PRECISION: f64 = 10.0e-24;
/// `RealLast()` as used at `IntPolyh_MaillageAffinage.cxx:1096-1097`,
/// `:2108`, `:2984`.
pub const REAL_LAST: f64 = f64::MAX;

/// `IntPolyh_MaillageAffinage` (`hxx:38`).
pub struct MaillageAffinage<'a> {
    /// `MaSurface1` / `MaSurface2` (`hxx:232-233`).
    pub(crate) surface1: &'a dyn Surface,
    pub(crate) surface2: &'a dyn Surface,
    /// `MyBox1` / `MyBox2` (`hxx:234-235`).
    pub(crate) box1: BndBox,
    pub(crate) box2: BndBox,
    /// `NbSamplesU1` / `NbSamplesU2` / `NbSamplesV1` / `NbSamplesV2`
    /// (`hxx:236-239`).
    pub(crate) nb_samples_u1: i32,
    pub(crate) nb_samples_u2: i32,
    pub(crate) nb_samples_v1: i32,
    pub(crate) nb_samples_v2: i32,
    /// `FlecheMax1` / `FlecheMax2` / `FlecheMin1` / `FlecheMin2`
    /// (`hxx:240-243`).
    pub(crate) fleche_max1: f64,
    pub(crate) fleche_max2: f64,
    pub(crate) fleche_min1: f64,
    pub(crate) fleche_min2: f64,
    /// `TPoints1` / `TPoints2` (`hxx:246-247`).
    pub(crate) points1: ArrayOfPoints,
    pub(crate) points2: ArrayOfPoints,
    /// `TEdges1` / `TEdges2` (`hxx:248-249`).
    pub(crate) edges1: ArrayOfEdges,
    pub(crate) edges2: ArrayOfEdges,
    /// `TTriangles1` / `TTriangles2` (`hxx:250-251`).
    pub(crate) triangles1: ArrayOfTriangles,
    pub(crate) triangles2: ArrayOfTriangles,
    /// `TTrianglesContacts` (`hxx:254`). OCCT uses an `NCollection_List`; only
    /// the insertion order is read by the algorithm.
    pub(crate) contacts: Vec<IntPolyhCouple>,
    /// `myEnlargeZone` (`hxx:257`).
    pub(crate) enlarge_zone: bool,
}

impl<'a> MaillageAffinage<'a> {
    /// `IntPolyh_MaillageAffinage(S1, S2, PRINT)` (`cxx:270-291`): the
    /// `PRINT` argument is unused, as in the `.cxx`.
    pub fn new(surface1: &'a dyn Surface, surface2: &'a dyn Surface) -> Self {
        Self::with_sampling(surface1, 10, 10, surface2, 10, 10)
    }

    /// `IntPolyh_MaillageAffinage(S1, NbSU1, NbSV1, S2, NbSU2, NbSV2, PRINT)`
    /// (`cxx:251-268`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_sampling(
        surface1: &'a dyn Surface,
        nb_su1: i32,
        nb_sv1: i32,
        surface2: &'a dyn Surface,
        nb_su2: i32,
        nb_sv2: i32,
    ) -> Self {
        Self {
            surface1,
            surface2,
            box1: BndBox::new(),
            box2: BndBox::new(),
            nb_samples_u1: nb_su1,
            nb_samples_u2: nb_su2,
            nb_samples_v1: nb_sv1,
            nb_samples_v2: nb_sv2,
            fleche_max1: 0.0,
            fleche_max2: 0.0,
            fleche_min1: 0.0,
            fleche_min2: 0.0,
            points1: ArrayOfPoints::default(),
            points2: ArrayOfPoints::default(),
            edges1: ArrayOfEdges::default(),
            edges2: ArrayOfEdges::default(),
            triangles1: ArrayOfTriangles::default(),
            triangles2: ArrayOfTriangles::default(),
            contacts: Vec::new(),
            enlarge_zone: false,
        }
    }

    /// The surface selected by `SurfID`.
    pub(crate) fn surface(&self, surf_id: i32) -> &'a dyn Surface {
        if surf_id == 1 {
            self.surface1
        } else {
            self.surface2
        }
    }

    /// `MakeSampling(SurfID, theUPars, theVPars)` (`cxx:293-315`).
    pub fn make_sampling(&self, surf_id: i32) -> (Vec<f64>, Vec<f64>) {
        let (nb_su, nb_sv) = if surf_id == 1 {
            (self.nb_samples_u1, self.nb_samples_v1)
        } else {
            (self.nb_samples_u2, self.nb_samples_v2)
        };
        tools::make_sampling(self.surface(surf_id), nb_su, nb_sv, self.enlarge_zone)
    }

    /// `FillArrayOfPnt(SurfID)` (`cxx:321-328`).
    pub fn fill_array_of_pnt(&mut self, surf_id: i32) {
        let (u_pars, v_pars) = self.make_sampling(surf_id);
        self.fill_array_of_pnt_with_params(surf_id, &u_pars, &v_pars, None);
    }

    /// `FillArrayOfPnt(SurfID, isShiftFwd)` (`cxx:335-341`).
    pub fn fill_array_of_pnt_shifted(&mut self, surf_id: i32, is_shift_fwd: bool) {
        let (u_pars, v_pars) = self.make_sampling(surf_id);
        self.fill_array_of_pnt_shifted_with_params(surf_id, is_shift_fwd, &u_pars, &v_pars, None);
    }

    /// `FillArrayOfPnt(SurfID, Upars, Vpars, theDeflTol)` (`cxx:348-408`).
    pub fn fill_array_of_pnt_with_params(
        &mut self,
        surf_id: i32,
        u_pars: &[f64],
        v_pars: &[f64],
        defl_tol: Option<f64>,
    ) {
        let surface = self.surface(surf_id);
        let (nb_u, nb_v) = if surf_id == 1 {
            (self.nb_samples_u1, self.nb_samples_v1)
        } else {
            (self.nb_samples_u2, self.nb_samples_v2)
        };

        let mut jd1 = 0;
        let mut jd2 = 0;
        let mut id1 = 0;
        let mut id2 = 0;
        degenerated_index(v_pars, nb_v, surface, 1, &mut jd1, &mut jd2);
        if !(jd1 != 0 || jd2 != 0) {
            degenerated_index(u_pars, nb_u, surface, 2, &mut id1, &mut id2);
        }

        let (points, a_box) = if surf_id == 1 {
            (&mut self.points1, &mut self.box1)
        } else {
            (&mut self.points2, &mut self.box2)
        };

        points.init(nb_u * nb_v);
        let mut i_cnt: i32 = 0;
        for i in 1..=nb_u {
            let deg_i = id1 == i || id2 == i;
            let a_u = u_pars[(i - 1) as usize];
            for j in 1..=nb_v {
                let a_v = v_pars[(j - 1) as usize];
                let p = surface.value(a_u, a_v);
                let ip = points.get_mut(i_cnt);
                ip.set(p.x(), p.y(), p.z(), a_u, a_v, 1);
                let deg = deg_i || jd1 == j || jd2 == j;
                if deg {
                    ip.degenerated = true;
                }
                i_cnt += 1;
                a_box.add_point(&p);
            }
        }
        points.set_nb_items(i_cnt);

        // `cxx:400-404`.
        let mut tol = match defl_tol {
            Some(t) => t,
            None => tools::compute_deflection(surface, u_pars, v_pars),
        };
        tol *= 1.2;

        // `cxx:407-408`.
        if let Some((a1, a2, a3, b1, b2, b3)) = a_box.get() {
            a_box.update(a1 - tol, a2 - tol, a3 - tol, b1 + tol, b2 + tol, b3 + tol);
        }
        a_box.enlarge(MY_TOLERANCE);
    }

    /// `FillArrayOfPnt(SurfID, isShiftFwd, thePointsNorm, theUPars, theVPars,
    /// theDeflTol)` (`cxx:413-478`).
    pub fn fill_array_of_pnt_with_normals(
        &mut self,
        surf_id: i32,
        is_shift_fwd: bool,
        points_norm: &super::arrays::ArrayOfPointNormal,
        u_pars: &[f64],
        v_pars: &[f64],
        defl_tol: f64,
    ) {
        let surface = self.surface(surf_id);
        let (nb_u, nb_v) = if surf_id == 1 {
            (self.nb_samples_u1, self.nb_samples_v1)
        } else {
            (self.nb_samples_u2, self.nb_samples_v2)
        };

        let mut jd1 = 0;
        let mut jd2 = 0;
        let mut id1 = 0;
        let mut id2 = 0;
        degenerated_index(v_pars, nb_v, surface, 1, &mut jd1, &mut jd2);
        if !(jd1 != 0 || jd2 != 0) {
            degenerated_index(u_pars, nb_u, surface, 2, &mut id1, &mut id2);
        }

        let (points, a_box) = if surf_id == 1 {
            (&mut self.points1, &mut self.box1)
        } else {
            (&mut self.points2, &mut self.box2)
        };

        points.init(points_norm.nb_items());

        let mut i_cnt: i32 = 0;
        for i in 1..=nb_u {
            let a_u = u_pars[(i - 1) as usize];
            let deg_i = id1 == i || id2 == i;
            for j in 1..=nb_v {
                let a_v = v_pars[(j - 1) as usize];

                let pn = points_norm.get(i_cnt);
                let mut norm = pn.normal.multiplied_scalar(1.5 * defl_tol);
                if !is_shift_fwd {
                    norm.reverse();
                }
                let p = occt_core::gp::GpPnt::new(
                    pn.point.x() + norm.x(),
                    pn.point.y() + norm.y(),
                    pn.point.z() + norm.z(),
                );

                let ip = points.get_mut(i_cnt);
                ip.set(p.x(), p.y(), p.z(), a_u, a_v, 1);
                let deg = deg_i || jd1 == j || jd2 == j;
                if deg {
                    ip.degenerated = true;
                }

                i_cnt += 1;
                a_box.add_point(&p);
            }
        }
        points.set_nb_items(i_cnt);

        // `cxx:472-477`.
        let tol = defl_tol * 1.2;
        if let Some((a1, a2, a3, b1, b2, b3)) = a_box.get() {
            a_box.update(a1 - tol, a2 - tol, a3 - tol, b1 + tol, b2 + tol, b3 + tol);
        }
        a_box.enlarge(MY_TOLERANCE);
    }

    /// `FillArrayOfPnt(SurfID, isShiftFwd, Upars, Vpars, theDeflTol)`
    /// (`cxx:483-499`).
    ///
    /// The `.cxx:498` calls `FillArrayOfPnt(1, ...)` with the surface id
    /// hard-coded to `1`; the port keeps that, it is not a typo.
    pub fn fill_array_of_pnt_shifted_with_params(
        &mut self,
        surf_id: i32,
        is_shift_fwd: bool,
        u_pars: &[f64],
        v_pars: &[f64],
        defl_tol: Option<f64>,
    ) {
        let surface = self.surface(surf_id);
        let a_tol = match defl_tol {
            Some(t) => t,
            None => tools::compute_deflection(surface, u_pars, v_pars),
        };

        let mut a_points = super::arrays::ArrayOfPointNormal::default();
        tools::fill_array_of_point_normal(surface, u_pars, v_pars, &mut a_points);

        // `cxx:498`: literally `FillArrayOfPnt(1, ...)`.
        self.fill_array_of_pnt_with_normals(1, is_shift_fwd, &a_points, u_pars, v_pars, a_tol);
    }


    /// `FillArrayOfEdges(SurfID)` (`cxx:770-951`).
    pub fn fill_array_of_edges(&mut self, surf_id: i32) {
        let (nb_su, nb_sv) = if surf_id == 1 {
            (self.nb_samples_u1, self.nb_samples_v1)
        } else {
            (self.nb_samples_u2, self.nb_samples_v2)
        };
        let (edges, triangles) = if surf_id == 1 {
            (&mut self.edges1, &mut self.triangles1)
        } else {
            (&mut self.edges2, &mut self.triangles2)
        };

        let nb_edges = 3 * nb_su * nb_sv - 2 * (nb_su + nb_sv) + 1;
        edges.init(nb_edges);

        let mut cpt: i32 = 0;

        // `cxx:788-807`: the first cell.
        {
            let e = edges.get_mut(cpt);
            e.point1 = 0;
            e.point2 = 1;
            e.triangle2 = 0;
            let ec = edges.get(cpt);
            triangles.get_mut(0).set_edge_and_orientation(&ec, cpt);
            cpt += 1;

            let e = edges.get_mut(cpt);
            e.point1 = 0;
            e.point2 = nb_sv;
            e.triangle1 = 1;
            let ec = edges.get(cpt);
            triangles.get_mut(1).set_edge_and_orientation(&ec, cpt);
            cpt += 1;

            let e = edges.get_mut(cpt);
            e.point1 = 0;
            e.point2 = nb_sv + 1;
            e.triangle1 = 0;
            e.triangle2 = 1;
            let ec = edges.get(cpt);
            triangles.get_mut(0).set_edge_and_orientation(&ec, cpt);
            triangles.get_mut(1).set_edge_and_orientation(&ec, cpt);
            cpt += 1;
        }

        // `cxx:809-836`: mesh along U = u0.
        let mut pnt_init = 1;
        for bmv in 1..nb_sv - 1 {
            {
                let e = edges.get_mut(cpt);
                e.point1 = pnt_init;
                e.point2 = pnt_init + 1;
                e.triangle2 = bmv * 2;
                let ec = edges.get(cpt);
                triangles.get_mut(bmv * 2).set_edge_and_orientation(&ec, cpt);
                cpt += 1;

                let e = edges.get_mut(cpt);
                e.point1 = pnt_init;
                e.point2 = pnt_init + nb_sv + 1;
                e.triangle1 = bmv * 2;
                e.triangle2 = bmv * 2 + 1;
                let ec = edges.get(cpt);
                triangles.get_mut(bmv * 2).set_edge_and_orientation(&ec, cpt);
                triangles.get_mut(bmv * 2 + 1).set_edge_and_orientation(&ec, cpt);
                cpt += 1;

                let e = edges.get_mut(cpt);
                e.point1 = pnt_init;
                e.point2 = pnt_init + nb_sv;
                e.triangle1 = bmv * 2 + 1;
                e.triangle2 = bmv * 2 - 2;
                let ec = edges.get(cpt);
                triangles.get_mut(bmv * 2 + 1).set_edge_and_orientation(&ec, cpt);
                triangles.get_mut(bmv * 2 - 2).set_edge_and_orientation(&ec, cpt);
                cpt += 1;
            }
            pnt_init += 1;
        }

        // `cxx:838-869`: mesh along V = v0.
        pnt_init = nb_sv;
        for bmu in 1..nb_su - 1 {
            {
                let e = edges.get_mut(cpt);
                e.point1 = pnt_init;
                e.point2 = pnt_init + 1;
                e.triangle1 = (bmu - 1) * (nb_sv - 1) * 2 + 1;
                e.triangle2 = bmu * (nb_sv - 1) * 2;
                let ec = edges.get(cpt);
                triangles
                    .get_mut((bmu - 1) * (nb_sv - 1) * 2 + 1)
                    .set_edge_and_orientation(&ec, cpt);
                triangles
                    .get_mut(bmu * (nb_sv - 1) * 2)
                    .set_edge_and_orientation(&ec, cpt);
                cpt += 1;

                let e = edges.get_mut(cpt);
                e.point1 = pnt_init;
                e.point2 = pnt_init + nb_sv + 1;
                e.triangle1 = bmu * (nb_sv - 1) * 2;
                e.triangle2 = bmu * (nb_sv - 1) * 2 + 1;
                let ec = edges.get(cpt);
                triangles
                    .get_mut(bmu * (nb_sv - 1) * 2)
                    .set_edge_and_orientation(&ec, cpt);
                triangles
                    .get_mut(bmu * (nb_sv - 1) * 2 + 1)
                    .set_edge_and_orientation(&ec, cpt);
                cpt += 1;

                let e = edges.get_mut(cpt);
                e.point1 = pnt_init;
                e.point2 = pnt_init + nb_sv;
                e.triangle1 = bmu * (nb_sv - 1) * 2 + 1;
                let ec = edges.get(cpt);
                triangles
                    .get_mut(bmu * (nb_sv - 1) * 2 + 1)
                    .set_edge_and_orientation(&ec, cpt);
                cpt += 1;
            }
            pnt_init += nb_sv;
        }

        // `cxx:871-921`: interior points, three edges per point.
        pnt_init = nb_sv + 1;
        for bmu in 1..nb_su - 1 {
            for bmv in 1..nb_sv - 1 {
                {
                    let e = edges.get_mut(cpt);
                    e.point1 = pnt_init;
                    e.point2 = pnt_init + 1;
                    e.triangle1 = (nb_sv - 1) * 2 * (bmu - 1) + bmv * 2 + 1;
                    e.triangle2 = (nb_sv - 1) * 2 * bmu + bmv * 2;
                    let ec = edges.get(cpt);
                    triangles
                        .get_mut((nb_sv - 1) * 2 * (bmu - 1) + bmv * 2 + 1)
                        .set_edge_and_orientation(&ec, cpt);
                    triangles
                        .get_mut((nb_sv - 1) * 2 * bmu + bmv * 2)
                        .set_edge_and_orientation(&ec, cpt);
                    cpt += 1;

                    let e = edges.get_mut(cpt);
                    e.point1 = pnt_init;
                    e.point2 = pnt_init + nb_sv + 1;
                    e.triangle1 = (nb_sv - 1) * 2 * bmu + bmv * 2;
                    e.triangle2 = (nb_sv - 1) * 2 * bmu + bmv * 2 + 1;
                    let ec = edges.get(cpt);
                    triangles
                        .get_mut((nb_sv - 1) * 2 * bmu + bmv * 2)
                        .set_edge_and_orientation(&ec, cpt);
                    triangles
                        .get_mut((nb_sv - 1) * 2 * bmu + bmv * 2 + 1)
                        .set_edge_and_orientation(&ec, cpt);
                    cpt += 1;

                    let e = edges.get_mut(cpt);
                    e.point1 = pnt_init;
                    e.point2 = pnt_init + nb_sv;
                    e.triangle1 = (nb_sv - 1) * 2 * bmu + bmv * 2 + 1;
                    e.triangle2 = (nb_sv - 1) * 2 * bmu + bmv * 2 - 2;
                    let ec = edges.get(cpt);
                    triangles
                        .get_mut((nb_sv - 1) * 2 * bmu + bmv * 2 + 1)
                        .set_edge_and_orientation(&ec, cpt);
                    triangles
                        .get_mut((nb_sv - 1) * 2 * bmu + bmv * 2 - 2)
                        .set_edge_and_orientation(&ec, cpt);
                    cpt += 1;
                }
                pnt_init += 1;
            }
            pnt_init += 1;
            pnt_init += 1;
        }

        // `cxx:923-935`: close the mesh at U = u1.
        pnt_init = (nb_su - 1) * nb_sv;
        for bmv in 0..nb_sv - 1 {
            let e = edges.get_mut(cpt);
            e.point1 = pnt_init;
            e.point2 = pnt_init + 1;
            e.triangle1 = (nb_su - 2) * (nb_sv - 1) * 2 + bmv * 2 + 1;
            let ec = edges.get(cpt);
            triangles
                .get_mut((nb_su - 2) * (nb_sv - 1) * 2 + bmv * 2 + 1)
                .set_edge_and_orientation(&ec, cpt);
            cpt += 1;
            pnt_init += 1;
        }

        // `cxx:937-949`: close the mesh at V = v1.
        for bmv in 0..nb_su - 1 {
            let e = edges.get_mut(cpt);
            e.point1 = nb_sv - 1 + bmv * nb_sv;
            e.point2 = nb_sv - 1 + (bmv + 1) * nb_sv;
            e.triangle2 = bmv * 2 * (nb_sv - 1) + (nb_sv - 2) * 2;
            let ec = edges.get(cpt);
            triangles
                .get_mut(bmv * 2 * (nb_sv - 1) + (nb_sv - 2) * 2)
                .set_edge_and_orientation(&ec, cpt);
            cpt += 1;
        }
        edges.set_nb_items(cpt);
    }

    /// `FillArrayOfTriangles(SurfID)` (`cxx:960-1019`).
    pub fn fill_array_of_triangles(&mut self, surf_id: i32) {
        let (nb_su, nb_sv) = if surf_id == 1 {
            (self.nb_samples_u1, self.nb_samples_v1)
        } else {
            (self.nb_samples_u2, self.nb_samples_v2)
        };
        let (points, triangles) = if surf_id == 1 {
            (&self.points1, &mut self.triangles1)
        } else {
            (&self.points2, &mut self.triangles2)
        };

        let mut cpt: i32 = 0;
        let mut pnt_init: i32 = 0;

        triangles.init(2 * (nb_su - 1) * (nb_sv - 1));
        for _bmu in 0..nb_su - 1 {
            for _bmv in 0..nb_sv - 1 {
                // FIRST TRIANGLE (`cxx:977-993`).
                {
                    let t = triangles.get_mut(cpt);
                    t.points = [pnt_init, pnt_init + 1, pnt_init + nb_sv + 1];
                    let p0 = points.get(pnt_init).part_of_common;
                    let p1 = points.get(pnt_init + 1).part_of_common;
                    let p2 = points.get(pnt_init + nb_sv + 1).part_of_common;
                    if (p0 & p1) != 0 && (p1 & p2) != 0 && (p2 & p0) != 0 {
                        t.is_intersection_possible = false;
                    }
                }
                cpt += 1;

                // SECOND TRIANGLE (`cxx:995-1009`).
                {
                    let t = triangles.get_mut(cpt);
                    t.points = [pnt_init, pnt_init + nb_sv + 1, pnt_init + nb_sv];
                    let p0 = points.get(pnt_init).part_of_common;
                    let p1 = points.get(pnt_init + nb_sv + 1).part_of_common;
                    let p2 = points.get(pnt_init + nb_sv).part_of_common;
                    if (p0 & p1) != 0 && (p1 & p2) != 0 && (p2 & p0) != 0 {
                        t.is_intersection_possible = false;
                    }
                }
                cpt += 1;

                pnt_init += 1;
            }
            pnt_init += 1;
        }
        triangles.set_nb_items(cpt);
    }

    /// `CommonPartRefinement()` (`cxx:1026-1046`).
    pub fn common_part_refinement(&mut self) {
        let fin_init1 = self.triangles1.nb_items();
        for i in 0..fin_init1 {
            if self.triangles1.get(i).is_intersection_possible {
                middle_refinement_at(
                    &mut self.triangles1,
                    i,
                    self.surface1,
                    &mut self.points1,
                    &mut self.edges1,
                );
            }
        }

        let fin_init2 = self.triangles2.nb_items();
        for i in 0..fin_init2 {
            if self.triangles2.get(i).is_intersection_possible {
                middle_refinement_at(
                    &mut self.triangles2,
                    i,
                    self.surface2,
                    &mut self.points2,
                    &mut self.edges2,
                );
            }
        }
    }

    /// `LocalSurfaceRefinement(SurfID)` (`cxx:1051-1077`).
    pub fn local_surface_refinement(&mut self, surf_id: i32) {
        if surf_id == 1 {
            let fin_init1 = self.triangles1.nb_items();
            for i in 0..fin_init1 {
                if self.triangles1.get(i).is_intersection_possible {
                    middle_refinement_at(
                        &mut self.triangles1,
                        i,
                        self.surface1,
                        &mut self.points1,
                        &mut self.edges1,
                    );
                }
            }
        }
        if surf_id == 2 {
            let fin_init2 = self.triangles2.nb_items();
            for i in 0..fin_init2 {
                if self.triangles2.get(i).is_intersection_possible {
                    middle_refinement_at(
                        &mut self.triangles2,
                        i,
                        self.surface2,
                        &mut self.points2,
                        &mut self.edges2,
                    );
                }
            }
        }
    }

    /// `ComputeDeflections(SurfID)` (`cxx:1088-1116`).
    pub fn compute_deflections(&mut self, surf_id: i32) {
        let (surface, points, triangles) = if surf_id == 1 {
            (self.surface1, &mut self.points1, &mut self.triangles1)
        } else {
            (self.surface2, &mut self.points2, &mut self.triangles2)
        };

        let mut fleche_max = -REAL_LAST;
        let mut fleche_min = REAL_LAST;
        let fin_tt = triangles.nb_items();
        for i in 0..fin_tt {
            let fleche = triangles.get_mut(i).compute_deflection(surface, points);
            if fleche > fleche_max {
                fleche_max = fleche;
            }
            if fleche < fleche_min {
                fleche_min = fleche;
            }
        }

        if surf_id == 1 {
            self.fleche_max1 = fleche_max;
            self.fleche_min1 = fleche_min;
        } else {
            self.fleche_max2 = fleche_max;
            self.fleche_min2 = fleche_min;
        }
    }

    /// `GetArrayOfPoints(SurfID)` (`cxx:3783-3790`).
    pub fn get_array_of_points(&self, surf_id: i32) -> &ArrayOfPoints {
        if surf_id == 1 {
            &self.points1
        } else {
            &self.points2
        }
    }

    /// `GetArrayOfEdges(SurfID)` (`cxx:3794-3802`).
    pub fn get_array_of_edges(&self, surf_id: i32) -> &ArrayOfEdges {
        if surf_id == 1 {
            &self.edges1
        } else {
            &self.edges2
        }
    }

    /// `GetArrayOfTriangles(SurfID)` (`cxx:3805-3814`).
    pub fn get_array_of_triangles(&self, surf_id: i32) -> &ArrayOfTriangles {
        if surf_id == 1 {
            &self.triangles1
        } else {
            &self.triangles2
        }
    }

    /// `GetFinTE(SurfID)`.
    pub fn get_fin_te(&self, surf_id: i32) -> i32 {
        if surf_id == 1 {
            self.edges1.nb_items()
        } else {
            self.edges2.nb_items()
        }
    }

    /// `GetFinTT(SurfID)`.
    pub fn get_fin_tt(&self, surf_id: i32) -> i32 {
        if surf_id == 1 {
            self.triangles1.nb_items()
        } else {
            self.triangles2.nb_items()
        }
    }

    /// `GetBox(SurfID)` (`cxx:3817-3826`).
    pub fn get_box(&self, surf_id: i32) -> BndBox {
        if surf_id == 1 {
            self.box1.clone()
        } else {
            self.box2.clone()
        }
    }

    /// `GetCouples()` (`cxx:3828-3831`).
    pub fn get_couples(&self) -> &[IntPolyhCouple] {
        &self.contacts
    }

    /// `GetCouples()` as a mutable list, for the couple bookkeeping of
    /// `StartPointsChain`.
    pub fn get_couples_mut(&mut self) -> &mut Vec<IntPolyhCouple> {
        &mut self.contacts
    }

    /// `SetEnlargeZone(EnlargeZone)` (`cxx:3835-3838`).
    pub fn set_enlarge_zone(&mut self, enlarge_zone: bool) {
        self.enlarge_zone = enlarge_zone;
    }

    /// `GetEnlargeZone()` (`cxx:3842-3845`).
    pub fn get_enlarge_zone(&self) -> bool {
        self.enlarge_zone
    }

    /// `GetMinDeflection(SurfID)` (`cxx:3849-3852`).
    pub fn get_min_deflection(&self, surf_id: i32) -> f64 {
        if surf_id == 1 {
            self.fleche_min1
        } else {
            self.fleche_min2
        }
    }

    /// `GetMaxDeflection(SurfID)` (`cxx:3856-3859`).
    pub fn get_max_deflection(&self, surf_id: i32) -> f64 {
        if surf_id == 1 {
            self.fleche_max1
        } else {
            self.fleche_max2
        }
    }
}

/// The three points of a triangle, resolved through an array of points.
pub(crate) fn triangle_points(
    triangle: &IntPolyhTriangle,
    points: &ArrayOfPoints,
) -> (IntPolyhPoint, IntPolyhPoint, IntPolyhPoint) {
    let (p1, p2, p3) = triangle.points3();
    (points.get(p1), points.get(p2), points.get(p3))
}
