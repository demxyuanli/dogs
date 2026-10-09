//! `IntCurveSurface_ThePolyhedronOfHInter` and
//! `IntCurveSurface_ThePolyhedronToolOfHInter`.
//!
//! Source: `IntCurveSurface_ThePolyhedronOfHInter.hxx/.cxx` and
//! `IntCurveSurface_PolyhedronUtils.pxx` (TKGeomAlgo). The polyhedron is the
//! `nbdeltaU x nbdeltaV` UV grid of surface points (two triangles per cell) that
//! carries the `(u, v)` of every node, so
//! [`super::section_point_params::section_point_to_parameters`] can invert a
//! section point back to surface parameters without any projection.
//!
//! Arrays are 1-based, exactly as OCCT indexes them (`C_MyPnts[Index]` with
//! `Index` from 1), so the `Triangle` / `TriConnex` / `IsOnBound` index
//! arithmetic is a literal transcription and can be passed straight to the
//! interference engine.

use occt_core::bnd::BndBox;
use occt_core::gp::{GpPnt, GpVec};
use occt_core::intf::IntfPolyhedronTool;
use occt_core::precision::epsilon;
use occt_geom::Surface;

/// `THE_MIN_EDGE_LENGTH_SQUARED` (`PolyhedronUtils.pxx:35`).
const THE_MIN_EDGE_LENGTH_SQUARED: f64 = 1e-15;

/// `THE_MIN_DEFLECTION` (`PolyhedronUtils.pxx:838`).
const THE_MIN_DEFLECTION: f64 = 0.0001;

/// `RealFirst()` (`Standard_Real.hxx`).
const REAL_FIRST: f64 = f64::MIN;

/// `IntCurveSurface_ThePolyhedronOfHInter` (`...hxx:27-157`).
#[derive(Clone)]
pub struct ThePolyhedronOfHInter {
    /// `nbdeltaU` (`...hxx:139`).
    nb_delta_u: usize,
    /// `nbdeltaV` (`...hxx:140`).
    nb_delta_v: usize,
    /// `TheBnd` (`...hxx:141`).
    the_bnd: BndBox,
    /// `TheComponentsBnd` (`...hxx:142`), 1-based: element 0 is unused.
    the_components_bnd: Vec<BndBox>,
    /// `TheDeflection` (`...hxx:143`).
    the_deflection: f64,
    /// `C_MyPnts` (`...hxx:146`), 1-based.
    pnts: Vec<GpPnt>,
    /// `C_MyU` (`...hxx:147`), 1-based.
    us: Vec<f64>,
    /// `C_MyV` (`...hxx:148`), 1-based.
    vs: Vec<f64>,
    /// `C_MyIsOnBounds` (`...hxx:154`), 1-based.
    is_on_bounds: Vec<bool>,
    /// `UMinSingular` (`...hxx:149`).
    umin_singular: bool,
    /// `UMaxSingular` (`...hxx:150`).
    umax_singular: bool,
    /// `VMinSingular` (`...hxx:151`).
    vmin_singular: bool,
    /// `VMaxSingular` (`...hxx:152`).
    vmax_singular: bool,
    /// `TheBorderDeflection` (`...hxx:153`).
    the_border_deflection: f64,
}

impl ThePolyhedronOfHInter {
    /// `ThePolyhedronOfHInter(Surface, nbdU, nbdV, u1, v1, u2, v2)`
    /// (`...cxx:46-68`).
    pub fn new(
        surface: &dyn Surface,
        nbd_u: usize,
        nbd_v: usize,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) -> Self {
        let nb_delta_u = if nbd_u < 3 { 3 } else { nbd_u };
        let nb_delta_v = if nbd_v < 3 { 3 } else { nbd_v };
        let mut me = Self::with_shape(nb_delta_u, nb_delta_v);
        me.init_uniform(surface, u1, v1, u2, v2);
        me
    }

    /// `ThePolyhedronOfHInter(Surface, Upars, Vpars)` (`...cxx:72-96`).
    pub fn with_params(surface: &dyn Surface, upars: &[f64], vpars: &[f64]) -> Self {
        assert!(upars.len() >= 2, "Upars must contain at least two values");
        assert!(vpars.len() >= 2, "Vpars must contain at least two values");
        let mut me = Self::with_shape(upars.len() - 1, vpars.len() - 1);
        me.init_with_params(surface, upars, vpars);
        me
    }

    /// `AllocateArrays` (`PolyhedronUtils.pxx:816-828`) plus the field
    /// initializers shared by both ctors (`...cxx:55-65`, `:82-92`).
    fn with_shape(nb_delta_u: usize, nb_delta_v: usize) -> Self {
        let n = (nb_delta_u + 1) * (nb_delta_v + 1) + 1;
        Self {
            nb_delta_u,
            nb_delta_v,
            the_bnd: BndBox::new(),
            the_components_bnd: Vec::new(),
            the_deflection: epsilon(100.0),
            pnts: vec![GpPnt::new(0.0, 0.0, 0.0); n],
            us: vec![0.0; n],
            vs: vec![0.0; n],
            is_on_bounds: vec![false; n],
            umin_singular: false,
            umax_singular: false,
            vmin_singular: false,
            vmax_singular: false,
            the_border_deflection: 0.0,
        }
    }

    /// `Init(Surface, U0, V0, U1, V1)` (`...cxx:107-137`) =
    /// `PolyhedronUtils::InitUniform` (`PolyhedronUtils.pxx:50-100`).
    fn init_uniform(&mut self, surface: &dyn Surface, u0: f64, v0: f64, u1: f64, v1: f64) {
        let du = (u1 - u0) / self.nb_delta_u as f64;
        let dv = (v1 - v0) / self.nb_delta_v as f64;
        let nb_u = self.nb_delta_u + 1;
        let nb_v = self.nb_delta_v + 1;
        let mut index = 1usize;
        for i1 in 0..nb_u {
            for i2 in 0..nb_v {
                let u = u0 + i1 as f64 * du;
                let v = v0 + i2 as f64 * dv;
                let p = surface.d0(u, v);
                self.pnts[index] = p;
                self.us[index] = u;
                self.vs[index] = v;
                self.is_on_bounds[index] =
                    i1 == 0 || i1 == self.nb_delta_u || i2 == 0 || i2 == self.nb_delta_v;
                self.the_bnd.add_point(&p);
                index += 1;
            }
        }

        let ntri = self.nb_triangles();
        let tol = self.compute_max_deflection(surface, ntri);
        self.set_deflection_over_estimation(tol * 1.2);
        self.fill_bounding();
        self.the_border_deflection =
            Self::compute_max_border_deflection(surface, u0, v0, u1, v1, self.nb_delta_u, self.nb_delta_v);
    }

    /// `Init(Surface, Upars, Vpars)` (`...cxx:139-168`) =
    /// `PolyhedronUtils::InitWithParams` (`PolyhedronUtils.pxx:114-145`).
    fn init_with_params(&mut self, surface: &dyn Surface, upars: &[f64], vpars: &[f64]) {
        for i1 in 0..=self.nb_delta_u {
            for i2 in 0..=self.nb_delta_v {
                let u = upars[i1];
                let v = vpars[i2];
                let index = i1 * (self.nb_delta_v + 1) + i2 + 1;
                let p = surface.d0(u, v);
                self.pnts[index] = p;
                self.us[index] = u;
                self.vs[index] = v;
                self.is_on_bounds[index] =
                    i1 == 0 || i1 == self.nb_delta_u || i2 == 0 || i2 == self.nb_delta_v;
                self.the_bnd.add_point(&p);
            }
        }

        let ntri = self.nb_triangles();
        let tol = self.compute_max_deflection(surface, ntri);
        self.set_deflection_over_estimation(tol * 1.2);
        self.fill_bounding();
        self.the_border_deflection = Self::compute_max_border_deflection(
            surface,
            upars[0],
            vpars[0],
            upars[upars.len() - 1],
            vpars[vpars.len() - 1],
            self.nb_delta_u,
            self.nb_delta_v,
        );
    }

    /// `NbTriangles()` (`...cxx:220-223` = `PolyhedronUtils.pxx:233-236`).
    pub fn nb_triangles(&self) -> usize {
        Self::triangles_of(self.nb_delta_u, self.nb_delta_v)
    }

    /// `NbPoints()` (`...cxx:227-230` = `PolyhedronUtils.pxx:242-245`).
    pub fn nb_points(&self) -> usize {
        Self::points_of(self.nb_delta_u, self.nb_delta_v)
    }

    /// `PolyUtils::NbTriangles` (`PolyhedronUtils.pxx:233-236`).
    pub fn triangles_of(nb_delta_u: usize, nb_delta_v: usize) -> usize {
        nb_delta_u * nb_delta_v * 2
    }

    /// `PolyUtils::NbPoints` (`PolyhedronUtils.pxx:242-245`).
    pub fn points_of(nb_delta_u: usize, nb_delta_v: usize) -> usize {
        (nb_delta_u + 1) * (nb_delta_v + 1)
    }

    /// `PolyUtils::Triangle` (`PolyhedronUtils.pxx:252-262`), 1-based.
    pub fn triangle_indices(&self, index: usize) -> (i32, i32, i32) {
        Self::triangle_of(index, self.nb_delta_v)
    }

    /// `PolyUtils::Triangle` (`PolyhedronUtils.pxx:252-262`), 1-based.
    pub fn triangle_of(index: usize, nb_delta_v: usize) -> (i32, i32, i32) {
        let line = 1 + (index - 1) / (nb_delta_v * 2);
        let colon = 1 + (index - 1) % (nb_delta_v * 2);
        let colpnt = (colon + 1) / 2;

        let p1 = (line - 1) * (nb_delta_v + 1) + colpnt;
        let p2 = line * (nb_delta_v + 1) + colpnt + (colon - 1) % 2;
        let p3 = (line - 1 + colon % 2) * (nb_delta_v + 1) + colpnt + 1;
        (p1 as i32, p2 as i32, p3 as i32)
    }

    /// `Parameters(Index, U, V)` (`...cxx:172-176` =
    /// `PolyhedronUtils.pxx:855-868`).
    pub fn parameters(&self, index: i32) -> (f64, f64) {
        (self.us[index as usize], self.vs[index as usize])
    }

    /// `Point(Index)` (`...cxx:296-299`).
    pub fn point(&self, index: i32) -> GpPnt {
        self.pnts[index as usize]
    }

    /// `Point(Index, U, V)` (`...cxx:286-293`).
    pub fn point_uv(&self, index: i32) -> (GpPnt, f64, f64) {
        (self.pnts[index as usize], self.us[index as usize], self.vs[index as usize])
    }

    /// `Bounding()` (`...cxx:192-196`).
    pub fn bounding(&self) -> &BndBox {
        &self.the_bnd
    }

    /// `ComponentsBounding()` (`...cxx:210-214`).
    pub fn components_bounding(&self) -> &[BndBox] {
        &self.the_components_bnd[1..]
    }

    /// `DeflectionOverEstimation()` (`...cxx:186-189`).
    pub fn deflection_over_estimation(&self) -> f64 {
        self.the_deflection
    }

    /// `DeflectionOverEstimation(flec)` (`...cxx:179-183` =
    /// `PolyhedronUtils.pxx:836-848`).
    pub fn set_deflection_over_estimation(&mut self, flec: f64) {
        if flec < THE_MIN_DEFLECTION {
            self.the_deflection = THE_MIN_DEFLECTION;
            self.the_bnd.enlarge(THE_MIN_DEFLECTION);
        } else {
            self.the_deflection = flec;
            self.the_bnd.enlarge(flec);
        }
    }

    /// `GetBorderDeflection()` (`...hxx:151`).
    pub fn get_border_deflection(&self) -> f64 {
        self.the_border_deflection
    }

    /// `FillBounding()` (`...cxx:200-208` =
    /// `PolyhedronUtils.pxx:578-608`).
    pub fn fill_bounding(&mut self) {
        let ntri = self.nb_triangles();
        self.the_components_bnd = vec![BndBox::new(); ntri + 1];
        for tri in 1..=ntri {
            let (n1, n2, n3) = self.triangle_indices(tri);
            let p1 = self.pnts[n1 as usize];
            let p2 = self.pnts[n2 as usize];
            let p3 = self.pnts[n3 as usize];
            let mut boite = BndBox::new();
            boite.set_void();
            if p1.square_distance(&p2) > THE_MIN_EDGE_LENGTH_SQUARED
                && p1.square_distance(&p3) > THE_MIN_EDGE_LENGTH_SQUARED
                && p2.square_distance(&p3) > THE_MIN_EDGE_LENGTH_SQUARED
            {
                boite.add_point(&p1);
                boite.add_point(&p2);
                boite.add_point(&p3);
                boite.enlarge(self.the_deflection);
            }
            boite.enlarge(self.the_deflection);
            self.the_components_bnd[tri] = boite;
        }
    }

    /// `Contain(Triang, ThePnt)` (`...cxx:256-262` =
    /// `PolyhedronUtils.pxx:556-567`).
    pub fn contain(&self, triang: usize, the_pnt: &GpPnt) -> bool {
        let (n1, n2, n3) = self.triangle_indices(triang);
        let p1 = self.pnts[n1 as usize];
        let p2 = self.pnts[n2 as usize];
        let p3 = self.pnts[n3 as usize];
        Self::contain_point(&p1, &p2, &p3, the_pnt)
    }

    /// `PolyUtils::Contain` (`PolyhedronUtils.pxx:556-567`).
    pub fn contain_point(p1: &GpPnt, p2: &GpPnt, p3: &GpPnt, the_pnt: &GpPnt) -> bool {
        let v1 = GpVec::from_pnts(p1, p2).crossed(&GpVec::from_pnts(p1, the_pnt));
        let v2 = GpVec::from_pnts(p2, p3).crossed(&GpVec::from_pnts(p2, the_pnt));
        let v3 = GpVec::from_pnts(p3, p1).crossed(&GpVec::from_pnts(p3, the_pnt));
        v1.dot(&v2) >= 0.0 && v2.dot(&v3) >= 0.0 && v3.dot(&v1) >= 0.0
    }

    /// `PlaneEquation(Triang, NormalVector, PolarDistance)` (`...cxx:245-254` =
    /// `PolyhedronUtils.pxx:508-549`). The `(1, 0, 0) / 0` degenerate answer is
    /// the OCCT one.
    pub fn plane_equation(&self, triang: usize) -> (GpVec, f64) {
        let (n1, n2, n3) = self.triangle_indices(triang);
        let p1 = self.pnts[n1 as usize];
        let p2 = self.pnts[n2 as usize];
        let p3 = self.pnts[n3 as usize];
        Self::plane_equation_of(&p1, &p2, &p3)
    }

    /// `PolyUtils::PlaneEquation` (`PolyhedronUtils.pxx:508-549`).
    pub fn plane_equation_of(p1: &GpPnt, p2: &GpPnt, p3: &GpPnt) -> (GpVec, f64) {
        let v1 = GpVec::from_pnts(p1, p2);
        let v2 = GpVec::from_pnts(p2, p3);
        let v3 = GpVec::from_pnts(p3, p1);

        if v1.square_magnitude() <= THE_MIN_EDGE_LENGTH_SQUARED
            || v2.square_magnitude() <= THE_MIN_EDGE_LENGTH_SQUARED
            || v3.square_magnitude() <= THE_MIN_EDGE_LENGTH_SQUARED
        {
            return (GpVec::new(1.0, 0.0, 0.0), 0.0);
        }

        let mut normal = v1.crossed(&v2).added(&v2.crossed(&v3)).added(&v3.crossed(&v1));
        let norm_len = normal.magnitude();
        if norm_len < f64::EPSILON {
            (normal, 0.0)
        } else {
            normal = normal.divided(norm_len);
            let polar = normal.dot(&GpVec::from_xyz(p1.xyz()));
            (normal, polar)
        }
    }

    /// `TriConnex(Triang, Pivot, Pedge, TriCon, OtherP)` (`...cxx:234-238` =
    /// `PolyhedronUtils.pxx:275-501`). Returns `(theTriCon, TriCon, OtherP)`.
    pub fn tri_connex(&self, triang: i32, pivot: i32, pedge: i32) -> (i32, i32, i32) {
        let nb_delta_v = self.nb_delta_v as i32;
        let nb_delta_u = self.nb_delta_u as i32;

        let pivot_m1 = pivot - 1;
        let nb_delta_vp1 = nb_delta_v + 1;
        let nb_delta_vm2 = nb_delta_v + nb_delta_v;

        let lig_p = pivot_m1 / nb_delta_vp1;
        let col_p = pivot_m1 - lig_p * nb_delta_vp1;

        let mut lig_e = 0i32;
        let mut col_e = 0i32;
        let mut typ_e = 0i32;
        if pedge != 0 {
            lig_e = (pedge - 1) / nb_delta_vp1;
            col_e = (pedge - 1) - (lig_e * nb_delta_vp1);
            if lig_p == lig_e {
                typ_e = 1;
            } else if col_p == col_e {
                typ_e = 2;
            } else {
                typ_e = 3;
            }
        }

        let mut lin_t = 0i32;
        let mut col_t = 0i32;
        let mut lin_o = 0i32;
        let mut col_o = 0i32;

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
                3 => {
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
                _ => {}
            }
        } else if pedge == 0 {
            lin_t = if 1 > lig_p { 1 } else { lig_p };
            col_t = if 1 > col_p + col_p { 1 } else { col_p + col_p };
            lin_o = if lig_p == 0 { lig_p + 1 } else { lig_p - 1 };
            col_o = col_p;
        } else {
            match typ_e {
                1 => {
                    lin_t = lig_p + 1;
                    col_t = col_p.max(col_e);
                    col_t += col_t;
                    lin_o = lig_p + 1;
                    col_o = col_p.max(col_e);
                }
                2 => {
                    lin_t = lig_p.max(lig_e);
                    col_t = col_p + col_p;
                    lin_o = lig_p.min(lig_e);
                    col_o = col_p - 1;
                }
                3 => {
                    lin_t = lig_p.max(lig_e);
                    col_t = col_p + col_e;
                    lin_o = lig_p.max(lig_e);
                    col_o = col_p.min(col_e);
                }
                _ => {}
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

        let other_p = lin_o * nb_delta_vp1 + col_o + 1;
        (tri_con, tri_con, other_p)
    }

    /// `IsOnBound(Index1, Index2)` (`...cxx:320-326` =
    /// `PolyhedronUtils.pxx:753-780`).
    pub fn is_on_bound(&self, index1: i32, index2: i32) -> bool {
        let nb_delta_u = self.nb_delta_u as i32;
        let nb_delta_v = self.nb_delta_v as i32;

        let diff = (index1 - index2).abs();
        if diff != 1 && diff != nb_delta_v + 1 {
            return false;
        }
        for i in 0..=nb_delta_u {
            if index1 == 1 + i * (nb_delta_v + 1) && index2 == index1 - 1 {
                return false;
            }
            if index1 == (1 + i) * (nb_delta_v + 1) && index2 == index1 + 1 {
                return false;
            }
        }
        self.is_on_bounds[index1 as usize] && self.is_on_bounds[index2 as usize]
    }

    /// `UMinSingularity(Sing)` (`...cxx:331-334`).
    pub fn set_umin_singularity(&mut self, sing: bool) {
        self.umin_singular = sing;
    }

    /// `UMaxSingularity(Sing)` (`...cxx:338-341`).
    pub fn set_umax_singularity(&mut self, sing: bool) {
        self.umax_singular = sing;
    }

    /// `VMinSingularity(Sing)` (`...cxx:345-348`).
    pub fn set_vmin_singularity(&mut self, sing: bool) {
        self.vmin_singular = sing;
    }

    /// `VMaxSingularity(Sing)` (`...cxx:352-355`).
    pub fn set_vmax_singularity(&mut self, sing: bool) {
        self.vmax_singular = sing;
    }

    /// `HasUMinSingularity()` (`...cxx:359-362`).
    pub fn has_umin_singularity(&self) -> bool {
        self.umin_singular
    }

    /// `HasUMaxSingularity()` (`...cxx:366-369`).
    pub fn has_umax_singularity(&self) -> bool {
        self.umax_singular
    }

    /// `HasVMinSingularity()` (`...cxx:373-376`).
    pub fn has_vmin_singularity(&self) -> bool {
        self.vmin_singular
    }

    /// `HasVMaxSingularity()` (`...cxx:380-383`).
    pub fn has_vmax_singularity(&self) -> bool {
        self.vmax_singular
    }

    /// `PolyUtils::ComputeMaxDeflection` (`PolyhedronUtils.pxx:664-694`).
    fn compute_max_deflection(&self, surface: &dyn Surface, nb_triangles: usize) -> f64 {
        if nb_triangles == 0 {
            return 0.0;
        }
        let mut tol = 0.0;
        for i in 1..=nb_triangles {
            let (n1, n2, n3) = self.triangle_indices(i);
            let (p1, u1, v1) = self.point_uv(n1);
            let (p2, u2, v2) = self.point_uv(n2);
            let (p3, u3, v3) = self.point_uv(n3);
            let u_center = (u1 + u2 + u3) / 3.0;
            let v_center = (v1 + v2 + v3) / 3.0;
            let center = surface.d0(u_center, v_center);
            let tol1 = Self::deflection_with_center(&p1, &p2, &p3, &center);
            if tol1 > tol {
                tol = tol1;
            }
        }
        tol
    }

    /// `PolyUtils::ComputeDeflectionWithCenter` (`PolyhedronUtils.pxx:616-653`).
    fn deflection_with_center(p1: &GpPnt, p2: &GpPnt, p3: &GpPnt, center: &GpPnt) -> f64 {
        if p1.square_distance(p2) <= THE_MIN_EDGE_LENGTH_SQUARED
            || p1.square_distance(p3) <= THE_MIN_EDGE_LENGTH_SQUARED
            || p2.square_distance(p3) <= THE_MIN_EDGE_LENGTH_SQUARED
        {
            return 0.0;
        }
        let xyz1 = GpVec::from_pnts(p1, p2);
        let xyz2 = GpVec::from_pnts(p2, p3);
        let xyz3 = GpVec::from_pnts(p3, p1);
        let mut normal = xyz1.crossed(&xyz2).added(&xyz2.crossed(&xyz3)).added(&xyz3.crossed(&xyz1));
        let norm_len = normal.magnitude();
        if norm_len < f64::EPSILON {
            return 0.0;
        }
        normal = normal.divided(norm_len);
        normal.dot(&GpVec::from_pnts(p1, center)).abs()
    }

    /// `PolyUtils::ComputeMaxBorderDeflection` (`PolyhedronUtils.pxx:707-745`).
    fn compute_max_border_deflection(
        surface: &dyn Surface,
        u0: f64,
        v0: f64,
        u1: f64,
        v1: f64,
        nb_delta_u: usize,
        nb_delta_v: usize,
    ) -> f64 {
        let mut max_deflection = REAL_FIRST;
        for d in [
            Self::compute_border_deflection(surface, u0, v0, v1, true, nb_delta_v),
            Self::compute_border_deflection(surface, u1, v0, v1, true, nb_delta_v),
            Self::compute_border_deflection(surface, v0, u0, u1, false, nb_delta_u),
            Self::compute_border_deflection(surface, v1, u0, u1, false, nb_delta_u),
        ] {
            if d > max_deflection {
                max_deflection = d;
            }
        }
        max_deflection
    }

    /// `PolyUtils::ComputeBorderDeflection` (`PolyhedronUtils.pxx:156-226`).
    /// `EvaluateGrid` is just `D0` at each parameter, so the grid is evaluated
    /// directly.
    fn compute_border_deflection(
        surface: &dyn Surface,
        parameter: f64,
        pmin: f64,
        pmax: f64,
        is_u_iso: bool,
        nb_samples: usize,
    ) -> f64 {
        if nb_samples == 0 {
            return 0.0;
        }
        let delta = (pmax - pmin) / nb_samples as f64;
        let eval = |varying: f64| {
            if is_u_iso {
                surface.d0(parameter, varying)
            } else {
                surface.d0(varying, parameter)
            }
        };
        let mut deflection = REAL_FIRST;
        for i in 0..nb_samples {
            let p1 = eval(pmin + i as f64 * delta);
            let p2 = eval(pmin + (i + 1) as f64 * delta);
            let par_mid = eval(pmin + (i as f64 + 0.5) * delta);
            let mid = GpPnt::new(
                0.5 * (p2.x() + p1.x()),
                0.5 * (p2.y() + p1.y()),
                0.5 * (p2.z() + p1.z()),
            );
            let dist = GpVec::from_pnts(&mid, &par_mid).magnitude();
            if dist > deflection {
                deflection = dist;
            }
        }
        deflection
    }

    /// Port-only: whether a curve polygon's AABB is disjoint from this
    /// polyhedron's. The `IntCurvesFace` quick reject
    /// (`Intf_Tool::PolyhedronBox`), not part of
    /// `IntCurveSurface_ThePolyhedronOfHInter`.
    pub fn is_out_polygon(&self, poly: &super::polygon::ThePolygonOfHInter) -> bool {
        if self.the_bnd.is_void() || poly.bounding().is_void() {
            return true;
        }
        self.the_bnd.is_out_box(poly.bounding())
    }
}

/// `IntCurveSurface_ThePolyhedronToolOfHInter` (`...hxx:27-112`).
pub struct ThePolyhedronToolOfHInter;

impl IntfPolyhedronTool for ThePolyhedronToolOfHInter {
    type Polyhedron = ThePolyhedronOfHInter;

    fn bounding(the_polyh: &Self::Polyhedron) -> &BndBox {
        the_polyh.bounding()
    }

    fn components_bounding(the_polyh: &Self::Polyhedron) -> &[BndBox] {
        the_polyh.components_bounding()
    }

    fn deflection_over_estimation(the_polyh: &Self::Polyhedron) -> f64 {
        the_polyh.deflection_over_estimation()
    }

    fn nb_triangles(the_polyh: &Self::Polyhedron) -> usize {
        the_polyh.nb_triangles()
    }

    fn triangle(the_polyh: &Self::Polyhedron, index: i32) -> (i32, i32, i32) {
        the_polyh.triangle_indices(index as usize)
    }

    fn point(the_polyh: &Self::Polyhedron, index: i32) -> GpPnt {
        the_polyh.point(index)
    }

    fn tri_connex(
        the_polyh: &Self::Polyhedron,
        triang: i32,
        pivot: i32,
        pedge: i32,
    ) -> (i32, i32, i32) {
        the_polyh.tri_connex(triang, pivot, pedge)
    }

    fn is_on_bound(the_polyh: &Self::Polyhedron, index1: i32, index2: i32) -> bool {
        the_polyh.is_on_bound(index1, index2)
    }

    fn get_border_deflection(the_polyh: &Self::Polyhedron) -> f64 {
        the_polyh.get_border_deflection()
    }
}
