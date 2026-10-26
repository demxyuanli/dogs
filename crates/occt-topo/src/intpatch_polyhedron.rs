//! `IntPatch_Polyhedron` and `IntPatch_PolyhedronTool`.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_Polyhedron.hxx` /
//! `.cxx` and `IntPatch_PolyhedronTool.hxx` / `.lxx`.
//!
//! The polyhedron is the `nbdeltaU x nbdeltaV` UV grid of surface points, two
//! triangles per cell, carrying the `(u, v)` of every node and the plane of
//! every triangle. It is the geometry `IntPatch_InterferencePolyhedron`
//! intersects and the one `IntPatch_RstInt` walks.
//!
//! Ported members (OCCT 8.0.0):
//! * `NbPOnU` / `NbPOnV` (`cxx:38-53`): `IntPatch_HInterTool::NbSamplesU/V`
//!   capped at `NBMAXUV` (30).
//! * both ctors (`cxx:70-136`, `cxx:138-206`).
//! * `Parameters` (`cxx:208-212`), `DeflectionOverEstimation` get/set
//!   (`cxx:216-236`), `Bounding` (`cxx:239-242`), `FillBounding` (`cxx:246-275`),
//!   `ComponentsBounding` (`cxx:278-281`), `NbTriangles` (`cxx:285-288`),
//!   `NbPoints` (`cxx:292-295`), `TriConnex` (`cxx:299-569`),
//!   `PlaneEquation` (`cxx:572-613`), `Contain` (`cxx:618-630`),
//!   `Size` (`cxx:638-642`), `Triangle` (`cxx:646-663`), the four `Point`
//!   overloads (`cxx:668-707`).
//! * `IntPatch_PolyhedronTool` (`hxx:34-78` / `lxx:17-64`) is the
//!   [`IntfPolyhedronTool`] implementation below - in OCCT it is a class of
//!   statics that forward to the very same members.
//!
//! Translation notes:
//! * the three `void*` arrays `C_MyPnts` / `C_MyU` / `C_MyV` (`hxx:142-144`)
//!   are the typeless storage OCCT shares with `IntCurveSurface_Polyhedron`
//!   (`cxx:57-66` frees them through `Destroy`). The port has no such sharing,
//!   so they become typed `Vec`s indexed 1-based with a dummy element 0, and
//!   `Destroy` is `Drop`.
//! * `Dump` (`cxx:634`) has an empty body in OCCT and is not translated.

use occt_core::bnd::BndBox;
use occt_core::gp::{GpPnt, GpVec};
use occt_core::intf::IntfPolyhedronTool;

use occt_geom::Surface;

use crate::brep_surface::classify_surface;
use crate::intpatch_h_inter_tool::{nb_samples_u, nb_samples_v, NB_MAX_UV};

/// `LONGUEUR_MINI_EDGE_TRIANGLE` (`cxx:31`).
const LONGUEUR_MINI_EDGE_TRIANGLE: f64 = 1e-14;

/// `DEFLECTION_COEFF` (`cxx:32`).
const DEFLECTION_COEFF: f64 = 1.1;

/// `IntPatch_Polyhedron` (`hxx:26-148`).
#[derive(Clone)]
pub struct IntPatchPolyhedron {
    /// `TheBnd` (`hxx:136`).
    the_bnd: BndBox,
    /// `TheComponentsBnd` (`hxx:137`), 1-based.
    the_components_bnd: Vec<BndBox>,
    /// `TheDeflection` (`hxx:138`).
    the_deflection: f64,
    /// `nbdeltaU` (`hxx:139`).
    nb_delta_u: usize,
    /// `nbdeltaV` (`hxx:140`).
    nb_delta_v: usize,
    /// `C_MyPnts` (`hxx:142`), 1-based.
    points: Vec<GpPnt>,
    /// `C_MyU` (`hxx:143`), 1-based.
    us: Vec<f64>,
    /// `C_MyV` (`hxx:144`), 1-based.
    vs: Vec<f64>,
    /// `UMinSingular` (`hxx:145`).
    umin_singular: bool,
    /// `UMaxSingular` (`hxx:146`).
    umax_singular: bool,
    /// `VMinSingular` (`hxx:147`).
    vmin_singular: bool,
    /// `VMaxSingular` (`hxx:148`).
    vmax_singular: bool,
}

/// `static int NbPOnU(S)` (`cxx:38-44`).
fn nb_p_on_u(s: &dyn Surface) -> usize {
    let nbpu = nb_samples_u(s);
    (if nbpu > NB_MAX_UV { NB_MAX_UV } else { nbpu }).max(1) as usize
}

/// `static int NbPOnV(S)` (`cxx:47-53`).
fn nb_p_on_v(s: &dyn Surface) -> usize {
    let nbpv = nb_samples_v(s);
    (if nbpv > NB_MAX_UV { NB_MAX_UV } else { nbpv }).max(1) as usize
}

impl IntPatchPolyhedron {
    /// `IntPatch_Polyhedron(Surface)` (`cxx:70-136`).
    pub fn new(surface: &dyn Surface) -> Self {
        Self::build(surface, nb_p_on_u(surface), nb_p_on_v(surface))
    }

    /// `IntPatch_Polyhedron(Surface, nbu, nbv)` (`cxx:138-206`).
    pub fn with_sizes(surface: &dyn Surface, nbu: i32, nbv: i32) -> Self {
        Self::build(surface, nbu.max(1) as usize, nbv.max(1) as usize)
    }

    /// The shared body of both ctors (`cxx:72-93` = `cxx:141-162`).
    fn build(surface: &dyn Surface, nb_delta_u: usize, nb_delta_v: usize) -> Self {
        let mut this = Self {
            the_bnd: BndBox::new(),
            the_components_bnd: Vec::new(),
            the_deflection: occt_core::precision::epsilon(100.0),
            nb_delta_u,
            nb_delta_v,
            points: Vec::new(),
            us: Vec::new(),
            vs: Vec::new(),
            umin_singular: crate::intpatch_h_inter_tool::singular_on_u_min(surface),
            umax_singular: crate::intpatch_h_inter_tool::singular_on_u_max(surface),
            vmin_singular: crate::intpatch_h_inter_tool::singular_on_v_min(surface),
            vmax_singular: crate::intpatch_h_inter_tool::singular_on_v_max(surface),
        };

        // `const int t = (nbdeltaU + 1) * (nbdeltaV + 1) + 1;`
        let t = (nb_delta_u + 1) * (nb_delta_v + 1) + 1;
        this.points = vec![GpPnt::new(0.0, 0.0, 0.0); t];
        this.us = vec![0.0; t];
        this.vs = vec![0.0; t];

        let (u0, u1) = surface.u_range();
        let (v0, v1) = surface.v_range();

        let du = (u1 - u0) / nb_delta_u as f64;
        let dv = (v1 - v0) / nb_delta_v as f64;

        let mut index = 1;
        for i1 in 0..=nb_delta_u {
            for i2 in 0..=nb_delta_v {
                let u = u0 + i1 as f64 * du;
                let v = v0 + i2 as f64 * dv;
                let p = surface.d0(u, v);
                this.points[index] = p;
                this.us[index] = u;
                this.vs[index] = v;
                this.the_bnd.add_point(&p);
                index += 1;
            }
        }

        // `double tol = PolyUtils::ComputeMaxDeflection(*Surface, *this,
        // NbTriangles()); tol *= DEFLECTION_COEFF;` (`cxx:128-131`)
        let mut tol = this.compute_max_deflection(surface);
        tol *= DEFLECTION_COEFF;

        this.set_deflection_over_estimation(tol);
        this.fill_bounding();
        this
    }

    /// `IntCurveSurface_PolyhedronUtils::ComputeMaxDeflection`
    /// (`PolyhedronUtils.pxx:658-694`) over this polyhedron's own triangles.
    fn compute_max_deflection(&self, surface: &dyn Surface) -> f64 {
        let ntri = self.nb_triangles();
        if ntri == 0 {
            return 0.0;
        }
        let mut tol: f64 = 0.0;
        for i in 1..=ntri {
            let (i1, i2, i3) = self.triangle(i);
            let (p1, u1, v1) = self.point_uv(i1);
            let (p2, u2, v2) = self.point_uv(i2);
            let (p3, u3, v3) = self.point_uv(i3);
            let u_center = (u1 + u2 + u3) / 3.0;
            let v_center = (v1 + v2 + v3) / 3.0;
            let a_center = surface.d0(u_center, v_center);
            let tol1 = Self::compute_deflection_with_center(&p1, &p2, &p3, &a_center);
            if tol1 > tol {
                tol = tol1;
            }
        }
        tol
    }

    /// `PolyhedronUtils::ComputeDeflectionWithCenter` (`pxx:615-653`).
    fn compute_deflection_with_center(p1: &GpPnt, p2: &GpPnt, p3: &GpPnt, center: &GpPnt) -> f64 {
        if p1.square_distance(p2) <= LONGUEUR_MINI_EDGE_TRIANGLE
            || p1.square_distance(p3) <= LONGUEUR_MINI_EDGE_TRIANGLE
            || p2.square_distance(p3) <= LONGUEUR_MINI_EDGE_TRIANGLE
        {
            return 0.0;
        }
        let xyz1 = GpVec::from_pnts(p1, p2);
        let xyz2 = GpVec::from_pnts(p2, p3);
        let xyz3 = GpVec::from_pnts(p3, p1);
        let mut normal = xyz1
            .crossed(&xyz2)
            .added(&xyz2.crossed(&xyz3))
            .added(&xyz3.crossed(&xyz1));
        let norm_len = normal.magnitude();
        if norm_len < f64::EPSILON {
            return 0.0;
        }
        normal = normal.divided(norm_len);
        let p1p = GpVec::from_pnts(p1, center);
        p1p.dot(&normal).abs()
    }

    /// `Parameters(Index, U, V)` (`cxx:208-212`).
    pub fn parameters(&self, index: usize) -> (f64, f64) {
        (self.us[index], self.vs[index])
    }

    /// `DeflectionOverEstimation(flec)` (`cxx:216-229`).
    pub fn set_deflection_over_estimation(&mut self, flec: f64) {
        if flec < 0.0001 {
            self.the_deflection = 0.0001;
            self.the_bnd.enlarge(0.0001);
        } else {
            self.the_deflection = flec;
            self.the_bnd.enlarge(flec);
        }
    }

    /// `DeflectionOverEstimation()` (`cxx:232-236`).
    pub fn deflection_over_estimation(&self) -> f64 {
        self.the_deflection
    }

    /// `Bounding()` (`cxx:239-242`).
    pub fn bounding(&self) -> &BndBox {
        &self.the_bnd
    }

    /// `FillBounding()` (`cxx:246-275`).
    pub fn fill_bounding(&mut self) {
        let ntri = self.nb_triangles();
        self.the_components_bnd = vec![BndBox::new(); ntri + 1];
        for i_tri in 1..=ntri {
            let (p1, p2, p3) = self.triangle(i_tri);
            let pp1 = self.points[p1 as usize];
            let pp2 = self.points[p2 as usize];
            let pp3 = self.points[p3 as usize];
            let mut boite = BndBox::new();
            boite.set_void();
            if pp1.square_distance(&pp2) > LONGUEUR_MINI_EDGE_TRIANGLE
                && pp1.square_distance(&pp3) > LONGUEUR_MINI_EDGE_TRIANGLE
                && pp2.square_distance(&pp3) > LONGUEUR_MINI_EDGE_TRIANGLE
            {
                boite.add_point(&pp1);
                boite.add_point(&pp2);
                boite.add_point(&pp3);
            }
            boite.enlarge(self.the_deflection);
            self.the_components_bnd[i_tri] = boite;
        }
    }

    /// `ComponentsBounding()` (`cxx:278-281`), 1-based.
    pub fn components_bounding(&self) -> &[BndBox] {
        &self.the_components_bnd
    }

    /// `NbTriangles()` (`cxx:285-288`): `nbdeltaU * nbdeltaV * 2`.
    pub fn nb_triangles(&self) -> usize {
        self.nb_delta_u * self.nb_delta_v * 2
    }

    /// `NbPoints()` (`cxx:292-295`): `(nbdeltaU + 1) * (nbdeltaV + 1)`.
    pub fn nb_points(&self) -> usize {
        (self.nb_delta_u + 1) * (self.nb_delta_v + 1)
    }

    /// `TriConnex(Triang, Pivot, Pedge, TriCon, OtherP)` (`cxx:299-569`).
    /// Returns `(return_value, TriCon, OtherP)`.
    pub fn tri_connex(&self, triang: i32, pivot: i32, pedge: i32) -> (i32, i32, i32) {
        let nb_delta_u = self.nb_delta_u as i32;
        let nb_delta_v = self.nb_delta_v as i32;

        let pivotm1 = pivot - 1;
        let nb_delta_vp1 = nb_delta_v + 1;
        let nb_delta_vm2 = nb_delta_v + nb_delta_v;

        // Pivot position in the MaTriangle
        let lig_p = pivotm1 / nb_delta_vp1;
        let col_p = pivotm1 - lig_p * nb_delta_vp1;

        // Point sur Edge position in the MaTriangle and edge typ
        let mut lig_e = 0i32;
        let mut col_e = 0i32;
        let typ_e;
        if pedge != 0 {
            lig_e = (pedge - 1) / nb_delta_vp1;
            col_e = (pedge - 1) - (lig_e * nb_delta_vp1);
            typ_e = if lig_p == lig_e {
                1
            } else if col_p == col_e {
                2
            } else {
                3
            };
        } else {
            typ_e = 0;
        }

        // Triangle position, general case
        let mut lin_t = 0i32;
        let mut col_t = 0i32;
        let mut lin_o = 0i32;
        let mut col_o = 0i32;
        let mut typ_e = typ_e;
        if triang != 0 {
            let t = (triang - 1) / nb_delta_vm2;
            let tt = (triang - 1) - t * nb_delta_vm2;
            lin_t = 1 + t;
            col_t = 1 + tt;
            if typ_e == 0 {
                if lig_p == lin_t {
                    lig_e = lig_p - 1;
                    col_e = col_p - 1;
                    typ_e = 3;
                } else if col_t == lig_p + lig_p {
                    lig_e = lig_p;
                    col_e = col_p - 1;
                    typ_e = 1;
                } else {
                    lig_e = lig_p + 1;
                    col_e = col_p + 1;
                    typ_e = 3;
                }
            }
            match typ_e {
                1 => {
                    // Horizontal
                    if lin_t == lig_p {
                        lin_t += 1;
                        lin_o = lig_p + 1;
                        col_o = col_p.max(col_e);
                    } else {
                        lin_t -= 1;
                        lin_o = lig_p - 1;
                        col_o = col_p.min(col_e);
                    }
                }
                2 => {
                    // Vertical
                    if col_t == col_p + col_p {
                        col_t += 1;
                        lin_o = lig_p.max(lig_e);
                        col_o = col_p + 1;
                    } else {
                        col_t -= 1;
                        lin_o = lig_p.min(lig_e);
                        col_o = col_p - 1;
                    }
                }
                _ => {
                    // Oblique
                    if (col_t & 1) == 0 {
                        col_t -= 1;
                        lin_o = lig_p.max(lig_e);
                        col_o = col_p.min(col_e);
                    } else {
                        col_t += 1;
                        lin_o = lig_p.min(lig_e);
                        col_o = col_p.max(col_e);
                    }
                }
            }
        } else {
            // Unknown Triangle position
            if pedge == 0 {
                // Unknown edge
                lin_t = 1.max(lig_p);
                col_t = 1.max(col_p + col_p);
                if lig_p == 0 {
                    lin_o = lig_p + 1;
                } else {
                    lin_o = lig_p - 1;
                }
                col_o = col_p;
            } else {
                // Known edge: take the left or down connectivity
                match typ_e {
                    1 => {
                        // Horizontal
                        lin_t = lig_p + 1;
                        col_t = col_p.max(col_e);
                        col_t += col_t;
                        lin_o = lig_p + 1;
                        col_o = col_p.max(col_e);
                    }
                    2 => {
                        // Vertical
                        lin_t = lig_p.max(lig_e);
                        col_t = col_p + col_p;
                        lin_o = lig_p.min(lig_e);
                        col_o = col_p - 1;
                    }
                    _ => {
                        // Oblique
                        lin_t = lig_p.max(lig_e);
                        col_t = col_p + col_e;
                        lin_o = lig_p.max(lig_e);
                        col_o = col_p.min(col_e);
                    }
                }
            }
        }

        let mut tri_con = (lin_t - 1) * nb_delta_vm2 + col_t;

        if lin_t < 1 {
            lin_o = 0;
            col_o = col_p + col_p - col_e;
            if col_o < 0 {
                col_o = 0;
                lin_o = 1;
            } else if col_o > nb_delta_v {
                col_o = nb_delta_v;
                lin_o = 1;
            }
            tri_con = 0;
        } else if lin_t > nb_delta_u {
            lin_o = nb_delta_u;
            col_o = col_p + col_p - col_e;
            if col_o < 0 {
                col_o = 0;
                lin_o = nb_delta_u - 1;
            } else if col_o > nb_delta_v {
                col_o = nb_delta_v;
                lin_o = nb_delta_u - 1;
            }
            tri_con = 0;
        }

        if col_t < 1 {
            col_o = 0;
            lin_o = lig_p + lig_p - lig_e;
            if lin_o < 0 {
                lin_o = 0;
                col_o = 1;
            } else if lin_o > nb_delta_u {
                lin_o = nb_delta_u;
                col_o = 1;
            }
            tri_con = 0;
        } else if col_t > nb_delta_v {
            col_o = nb_delta_v;
            lin_o = lig_p + lig_p - lig_e;
            if lin_o < 0 {
                lin_o = 0;
                col_o = nb_delta_v - 1;
            } else if lin_o > nb_delta_u {
                lin_o = nb_delta_u;
                col_o = nb_delta_v - 1;
            }
            tri_con = 0;
        }

        let mut other_p = lin_o * nb_delta_vp1 + col_o + 1;

        // Detection of the cases where the returned triangle is invalid
        // (`cxx:536-565`).
        if pedge != 0
            && self.point(pivot).square_distance(&self.point(pedge)) <= LONGUEUR_MINI_EDGE_TRIANGLE
        {
            other_p = 0;
            tri_con = triang;
            return (tri_con, tri_con, other_p);
        }
        if pedge != 0
            && self.point(other_p).square_distance(&self.point(pedge))
                <= LONGUEUR_MINI_EDGE_TRIANGLE
        {
            return (0, tri_con, other_p);
        }
        (tri_con, tri_con, other_p)
    }

    /// `PlaneEquation(Triang, NormalVector, PolarDistance)` (`cxx:572-613`).
    pub fn plane_equation(&self, triang: usize) -> (GpVec, f64) {
        let (i1, i2, i3) = self.triangle(triang);
        let p1 = self.points[i1 as usize].xyz();
        let p2 = self.points[i2 as usize].xyz();
        let p3 = self.points[i3 as usize].xyz();

        let v1 = p2.subtracted(&p1);
        let v2 = p3.subtracted(&p2);
        let v3 = p1.subtracted(&p3);

        if v1.square_modulus() <= LONGUEUR_MINI_EDGE_TRIANGLE
            || v2.square_modulus() <= LONGUEUR_MINI_EDGE_TRIANGLE
            || v3.square_modulus() <= LONGUEUR_MINI_EDGE_TRIANGLE
        {
            return (GpVec::new(1.0, 0.0, 0.0), 0.0);
        }

        let mut normal = v1
            .crossed(&v2)
            .added(&v2.crossed(&v3))
            .added(&v3.crossed(&v1));
        let norm_len = normal.modulus();
        if norm_len < f64::EPSILON {
            (GpVec::from_xyz(&normal), 0.0)
        } else {
            normal = normal.divided(norm_len);
            let polar = normal.dot(&p1);
            (GpVec::from_xyz(&normal), polar)
        }
    }

    /// `Contain(Triang, ThePnt)` (`cxx:618-630`).
    pub fn contain(&self, triang: usize, the_pnt: &GpPnt) -> bool {
        let (i1, i2, i3) = self.triangle(triang);
        let p1 = self.points[i1 as usize].xyz();
        let p2 = self.points[i2 as usize].xyz();
        let p3 = self.points[i3 as usize].xyz();

        let v1 = p2.subtracted(&p1).crossed(&the_pnt.xyz().subtracted(&p1));
        let v2 = p3.subtracted(&p2).crossed(&the_pnt.xyz().subtracted(&p2));
        let v3 = p1.subtracted(&p3).crossed(&the_pnt.xyz().subtracted(&p3));
        v1.dot(&v2) >= 0.0 && v2.dot(&v3) >= 0.0 && v3.dot(&v1) >= 0.0
    }

    /// `Size(nbdu, nbdv)` (`cxx:638-642`).
    pub fn size(&self) -> (usize, usize) {
        (self.nb_delta_u, self.nb_delta_v)
    }

    /// `Triangle(Index, P1, P2, P3)` (`cxx:646-663`), 1-based.
    pub fn triangle(&self, index: usize) -> (i32, i32, i32) {
        let nb_delta_v = self.nb_delta_v as i32;
        let line = 1 + ((index as i32 - 1) / (nb_delta_v * 2));
        let colon = 1 + ((index as i32 - 1) % (nb_delta_v * 2));
        let colpnt = (colon + 1) / 2;

        let p1 = (line - 1) * (nb_delta_v + 1) + colpnt;
        let p2 = line * (nb_delta_v + 1) + colpnt + ((colon - 1) % 2);
        let p3 = (line - 1 + (colon % 2)) * (nb_delta_v + 1) + colpnt + 1;
        (p1, p2, p3)
    }

    /// `Point(Index, U, V)` (`cxx:668-675`).
    pub fn point_uv(&self, index: i32) -> (GpPnt, f64, f64) {
        let i = index as usize;
        (self.points[i], self.us[i], self.vs[i])
    }

    /// `Point(Index)` (`cxx:681-685`).
    pub fn point(&self, index: i32) -> GpPnt {
        self.points[index as usize]
    }

    /// `Point(Index, gp_Pnt&)` (`cxx:702-706`).
    pub fn point_out(&self, index: i32) -> GpPnt {
        self.points[index as usize]
    }

    /// `HasUMinSingularity()` / `HasUMaxSingularity()` /
    /// `HasVMinSingularity()` / `HasVMaxSingularity()` (`hxx:114-121`).
    pub fn has_umin_singularity(&self) -> bool {
        self.umin_singular
    }
    pub fn has_umax_singularity(&self) -> bool {
        self.umax_singular
    }
    pub fn has_vmin_singularity(&self) -> bool {
        self.vmin_singular
    }
    pub fn has_vmax_singularity(&self) -> bool {
        self.vmax_singular
    }
}

impl IntfPolyhedronTool for IntPatchPolyhedron {
    type Polyhedron = IntPatchPolyhedron;

    /// `IntPatch_PolyhedronTool::Bounding` (`lxx:17-21`).
    fn bounding(polyh: &Self::Polyhedron) -> &BndBox {
        polyh.bounding()
    }

    /// `ComponentsBounding` (`lxx:23-27`), 1-based.
    fn components_bounding(polyh: &Self::Polyhedron) -> &[BndBox] {
        polyh.components_bounding()
    }

    /// `DeflectionOverEstimation` (`lxx:29-33`).
    fn deflection_over_estimation(polyh: &Self::Polyhedron) -> f64 {
        polyh.deflection_over_estimation()
    }

    /// `NbTriangles` (`lxx:35-39`).
    fn nb_triangles(polyh: &Self::Polyhedron) -> usize {
        polyh.nb_triangles()
    }

    /// `Triangle` (`lxx:41-47`), 1-based.
    fn triangle(polyh: &Self::Polyhedron, index: i32) -> (i32, i32, i32) {
        polyh.triangle(index as usize)
    }

    /// `Point` (`lxx:49-53`).
    fn point(polyh: &Self::Polyhedron, index: i32) -> GpPnt {
        polyh.point(index)
    }

    /// `TriConnex` (`lxx:55-64`).
    fn tri_connex(
        polyh: &Self::Polyhedron,
        triang: i32,
        pivot: i32,
        pedge: i32,
    ) -> (i32, i32, i32) {
        polyh.tri_connex(triang, pivot, pedge)
    }

    /// `IntPatch_PolyhedronTool` has no `IsOnBound` / `GetBorderDeflection`
    /// (`hxx:34-78`), unlike `IntCurveSurface_ThePolyhedronToolOfHInter`. The
    /// `Intf` generic engine never reaches those two for a polyhedron built
    /// this way, so the conservative answers are `false` / `0`.
    fn is_on_bound(_polyh: &Self::Polyhedron, _index1: i32, _index2: i32) -> bool {
        false
    }

    /// See [`Self::is_on_bound`].
    fn border_deflection(_polyh: &Self::Polyhedron) -> f64 {
        0.0
    }
}

impl IntPatchPolyhedron {
    /// The `classify_surface` helper kept next to the ctor for callers that
    /// need the same `GeomAdaptor_Surface` type split (unused by OCCT's
    /// `IntPatch_Polyhedron`, which only asks the sample counts).
    pub fn surface_kind(s: &dyn Surface) -> crate::brep_surface::SurfaceKind {
        classify_surface(s)
    }
}
