//! Port of ShapeExtend_CompositeSurface
//! (ShapeExtend_CompositeSurface.hxx/.cxx, 758 lines).
//!
//! A 2D grid of Geom_Surface patches with joint parameter values, used by
//! ShapeFix_ComposeShell as myGrid. OCCT's NCollection_HArray2 grid is
//! represented here as patches[i][j] where the first index is the U direction
//! (OCCT ColLength) and the second the V direction (OCCT RowLength), both
//! 1-based in the API.

use std::sync::Arc;

use occt_core::gp::{GpPnt2d, GpTrsf2d, GpVec2d, TrsfForm};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::Surface;

/// ShapeExtend_Parametrisation (ShapeExtend_Parametrisation.hxx).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Parametrisation {
    Natural,
    Uniform,
    Unitary,
}

/// ShapeExtend_CompositeSurface.
#[derive(Clone, Default)]
pub struct CompositeSurface {
    /// myPatches: patches[u - 1][v - 1].
    patches: Vec<Vec<Arc<dyn Surface>>>,
    /// myUJointValues (length NbUPatches + 1).
    u_joint_values: Vec<f64>,
    /// myVJointValues (length NbVPatches + 1).
    v_joint_values: Vec<f64>,
    /// myUClosed / myVClosed.
    u_closed: bool,
    v_closed: bool,
}

/// LimitValue (cxx:657-660).
fn limit_value(par: f64) -> f64 {
    if par.is_infinite() {
        if par < 0.0 {
            -10000.0
        } else {
            10000.0
        }
    } else {
        par
    }
}

impl CompositeSurface {
    /// ShapeExtend_CompositeSurface() (cxx:32).
    pub fn new() -> Self {
        Self::default()
    }

    /// ShapeExtend_CompositeSurface(GridSurf, param) (cxx:36-41).
    pub fn with_grid(patches: Vec<Vec<Arc<dyn Surface>>>, param: Parametrisation) -> Self {
        let mut s = Self::new();
        s.init(patches, param);
        s
    }

    /// Init(GridSurf, param) (cxx:55-66).
    pub fn init(&mut self, patches: Vec<Vec<Arc<dyn Surface>>>, param: Parametrisation) -> bool {
        if patches.is_empty() {
            return false;
        }
        self.patches = patches;
        self.compute_joint_values(param);
        self.check_connectivity(CONFUSION)
    }

    /// Init(GridSurf, UJoints, VJoints) (cxx:70-92).
    pub fn init_joints(
        &mut self,
        patches: Vec<Vec<Arc<dyn Surface>>>,
        u_joints: &[f64],
        v_joints: &[f64],
    ) -> bool {
        if patches.is_empty() {
            return false;
        }
        self.patches = patches;
        let mut ok = true;
        if !self.set_u_joint_values(u_joints) || !self.set_v_joint_values(v_joints) {
            ok = false;
            self.compute_joint_values(Parametrisation::Natural);
        }
        if !self.check_connectivity(CONFUSION) {
            return false;
        }
        ok
    }

    /// NbUPatches() (cxx:96-99) = ColLength.
    pub fn nb_u_patches(&self) -> usize {
        self.patches.len()
    }

    /// NbVPatches() (cxx:103-106) = RowLength.
    pub fn nb_v_patches(&self) -> usize {
        self.patches.first().map_or(0, |r| r.len())
    }

    /// Patch(i, j) (cxx:110-113), 1-based.
    pub fn patch(&self, i: usize, j: usize) -> Option<&Arc<dyn Surface>> {
        self.patches.get(i.checked_sub(1)?)?.get(j.checked_sub(1)?)
    }

    /// UJointValues() (cxx:125-128).
    pub fn u_joint_values(&self) -> &[f64] {
        &self.u_joint_values
    }

    /// VJointValues() (cxx:132-135).
    pub fn v_joint_values(&self) -> &[f64] {
        &self.v_joint_values
    }

    /// UJointValue(i) (cxx:139-142), 1-based.
    pub fn u_joint_value(&self, i: usize) -> f64 {
        self.u_joint_values.get(i - 1).copied().unwrap_or(0.0)
    }

    /// VJointValue(i) (cxx:146-149), 1-based.
    pub fn v_joint_value(&self, i: usize) -> f64 {
        self.v_joint_values.get(i - 1).copied().unwrap_or(0.0)
    }

    /// SetUJointValues (cxx:153-173): length must be NbU + 1 and strictly
    /// increasing above PConfusion.
    pub fn set_u_joint_values(&mut self, joints: &[f64]) -> bool {
        let nb_u = self.nb_u_patches();
        if joints.len() != nb_u + 1 {
            return false;
        }
        for i in 1..joints.len() {
            if joints[i] - joints[i - 1] < PCONFUSION {
                return false;
            }
        }
        self.u_joint_values = joints.to_vec();
        true
    }

    /// SetVJointValues (cxx:177-197).
    pub fn set_v_joint_values(&mut self, joints: &[f64]) -> bool {
        let nb_v = self.nb_v_patches();
        if joints.len() != nb_v + 1 {
            return false;
        }
        for i in 1..joints.len() {
            if joints[i] - joints[i - 1] < PCONFUSION {
                return false;
            }
        }
        self.v_joint_values = joints.to_vec();
        true
    }

    /// SetUFirstValue (cxx:201-214).
    pub fn set_u_first_value(&mut self, u_first: f64) {
        if self.u_joint_values.is_empty() {
            return;
        }
        let shift = u_first - self.u_joint_values[0];
        for v in self.u_joint_values.iter_mut() {
            *v += shift;
        }
    }

    /// SetVFirstValue (cxx:218-231).
    pub fn set_v_first_value(&mut self, v_first: f64) {
        if self.v_joint_values.is_empty() {
            return;
        }
        let shift = v_first - self.v_joint_values[0];
        for v in self.v_joint_values.iter_mut() {
            *v += shift;
        }
    }

    /// LocateUParameter (cxx:235-246).
    pub fn locate_u_parameter(&self, u: f64) -> usize {
        let nb = self.nb_u_patches();
        for i in 2..=nb {
            if u < self.u_joint_value(i) {
                return i - 1;
            }
        }
        nb
    }

    /// LocateVParameter (cxx:250-261).
    pub fn locate_v_parameter(&self, v: f64) -> usize {
        let nb = self.nb_v_patches();
        for i in 2..=nb {
            if v < self.v_joint_value(i) {
                return i - 1;
            }
        }
        nb
    }

    /// LocateUVPoint (cxx:265-269).
    pub fn locate_uv_point(&self, pnt: &GpPnt2d) -> (usize, usize) {
        (
            self.locate_u_parameter(pnt.x()),
            self.locate_v_parameter(pnt.y()),
        )
    }

    /// Patch(U, V) (cxx:273-277).
    pub fn patch_uv(&self, u: f64, v: f64) -> Option<&Arc<dyn Surface>> {
        self.patch(self.locate_u_parameter(u), self.locate_v_parameter(v))
    }

    /// Patch(pnt) (cxx:281-284).
    pub fn patch_pnt(&self, pnt: &GpPnt2d) -> Option<&Arc<dyn Surface>> {
        let (i, j) = self.locate_uv_point(pnt);
        self.patch(i, j)
    }

    /// ULocalToGlobal (cxx:288-296).
    pub fn u_local_to_global(&self, i: usize, j: usize, u: f64) -> f64 {
        let Some(p) = self.patch(i, j) else { return u };
        let (u1, u2, _, _) = (p.u_range().0, p.u_range().1, 0.0, 0.0);
        let scale = (self.u_joint_value(i + 1) - self.u_joint_value(i)) / (u2 - u1);
        u * scale + (self.u_joint_value(i) - u1 * scale)
    }

    /// VLocalToGlobal (cxx:300-308).
    pub fn v_local_to_global(&self, i: usize, j: usize, v: f64) -> f64 {
        let Some(p) = self.patch(i, j) else { return v };
        let (v1, v2) = (p.v_range().0, p.v_range().1);
        let scale = (self.v_joint_value(j + 1) - self.v_joint_value(j)) / (v2 - v1);
        v * scale + (self.v_joint_value(j) - v1 * scale)
    }

    /// LocalToGlobal (cxx:312-324).
    pub fn local_to_global(&self, i: usize, j: usize, uv: &GpPnt2d) -> GpPnt2d {
        GpPnt2d::new(
            self.u_local_to_global(i, j, uv.x()),
            self.v_local_to_global(i, j, uv.y()),
        )
    }

    /// UGlobalToLocal (cxx:328-336).
    pub fn u_global_to_local(&self, i: usize, j: usize, u: f64) -> f64 {
        let Some(p) = self.patch(i, j) else { return u };
        let (u1, u2) = (p.u_range().0, p.u_range().1);
        let scale = (u2 - u1) / (self.u_joint_value(i + 1) - self.u_joint_value(i));
        u * scale + (u1 - self.u_joint_value(i) * scale)
    }

    /// VGlobalToLocal (cxx:340-348).
    pub fn v_global_to_local(&self, i: usize, j: usize, v: f64) -> f64 {
        let Some(p) = self.patch(i, j) else { return v };
        let (v1, v2) = (p.v_range().0, p.v_range().1);
        let scale = (v2 - v1) / (self.v_joint_value(j + 1) - self.v_joint_value(j));
        v * scale + (v1 - self.v_joint_value(j) * scale)
    }

    /// GlobalToLocal (cxx:352-365).
    pub fn global_to_local(&self, i: usize, j: usize, uv: &GpPnt2d) -> GpPnt2d {
        GpPnt2d::new(
            self.u_global_to_local(i, j, uv.x()),
            self.v_global_to_local(i, j, uv.y()),
        )
    }

    /// GlobalToLocalTransformation (cxx:369-393): returns (uFact, Trsf,
    /// needT) where needT is the OCCT return value
    /// (uFact != 1 || Trsf.Form() != gp_Identity).
    pub fn global_to_local_transformation(
        &self,
        i: usize,
        j: usize,
    ) -> (f64, GpTrsf2d, bool) {
        let Some(p) = self.patch(i, j) else {
            return (1.0, GpTrsf2d::identity(), false);
        };
        let (u1, u2) = (p.u_range().0, p.u_range().1);
        let (v1, v2) = (p.v_range().0, p.v_range().1);
        let scaleu = (u2 - u1) / (self.u_joint_value(i + 1) - self.u_joint_value(i));
        let scalev = (v2 - v1) / (self.v_joint_value(j + 1) - self.v_joint_value(j));
        let shift = GpVec2d::new(
            u1 / scaleu - self.u_joint_value(i),
            v1 / scalev - self.v_joint_value(j),
        );
        let u_fact = scaleu / scalev;
        let mut trsf = GpTrsf2d::identity();
        let mut val = GpTrsf2d::identity();
        if shift.x() != 0.0 || shift.y() != 0.0 {
            val.set_translation_vec(&shift);
            trsf = val.multiplied(&trsf);
        }
        if scalev != 1.0 {
            let mut sc = GpTrsf2d::identity();
            let _ = sc.set_scale(&GpPnt2d::new(0.0, 0.0), scalev);
            trsf = sc.multiplied(&trsf);
        }
        let need_t = u_fact != 1.0 || trsf.form() != TrsfForm::Identity;
        (u_fact, trsf, need_t)
    }

    /// Bounds (cxx:463-469).
    pub fn bounds(&self) -> (f64, f64, f64, f64) {
        (
            self.u_joint_value(1),
            self.u_joint_value(self.nb_u_patches() + 1),
            self.v_joint_value(1),
            self.v_joint_value(self.nb_v_patches() + 1),
        )
    }

    /// IsUPeriodic (cxx:473-476).
    pub fn is_u_periodic(&self) -> bool {
        false
    }

    /// IsVPeriodic (cxx:480-483).
    pub fn is_v_periodic(&self) -> bool {
        false
    }

    /// IsUClosed (cxx:524-527).
    pub fn is_u_closed(&self) -> bool {
        self.u_closed
    }

    /// IsVClosed (cxx:531-534).
    pub fn is_v_closed(&self) -> bool {
        self.v_closed
    }

    /// Value(pnt) (cxx:591-599).
    pub fn value_pnt(&self, pnt: &GpPnt2d) -> occt_core::gp::GpPnt {
        let (i, j) = self.locate_uv_point(pnt);
        let uv = self.global_to_local(i, j, pnt);
        match self.patch(i, j) {
            Some(p) => p.d0(uv.x(), uv.y()),
            None => occt_core::gp::GpPnt::zero(),
        }
    }

    /// Value(U, V) (EvalD0 / cxx:538-544).
    pub fn value_uv(&self, u: f64, v: f64) -> occt_core::gp::GpPnt {
        self.value_pnt(&GpPnt2d::new(u, v))
    }

    /// ComputeJointValues (cxx:603-653).
    pub fn compute_joint_values(&mut self, param: Parametrisation) {
        let nb_u = self.nb_u_patches();
        let nb_v = self.nb_v_patches();
        self.u_joint_values = vec![0.0; nb_u + 1];
        self.v_joint_values = vec![0.0; nb_v + 1];
        match param {
            Parametrisation::Natural => {
                let mut u = 0.0;
                for i in 1..=nb_u {
                    let r = self.patch(i, 1).map(|p| p.u_range()).unwrap_or((0.0, 1.0));
                    if i == 1 {
                        u = r.0;
                        self.u_joint_values[0] = u;
                    }
                    u += r.1 - r.0;
                    self.u_joint_values[i] = u;
                }
                let mut v = 0.0;
                for i in 1..=nb_v {
                    let r = self.patch(1, i).map(|p| p.v_range()).unwrap_or((0.0, 1.0));
                    if i == 1 {
                        v = r.0;
                        self.v_joint_values[0] = v;
                    }
                    v += r.1 - r.0;
                    self.v_joint_values[i] = v;
                }
            }
            _ => {
                let (mut stepu, mut stepv) = (1.0, 1.0);
                if param == Parametrisation::Unitary {
                    stepu /= nb_u.max(1) as f64;
                    stepv /= nb_v.max(1) as f64;
                }
                for i in 0..=nb_u {
                    self.u_joint_values[i] = i as f64 * stepu;
                }
                for i in 0..=nb_v {
                    self.v_joint_values[i] = i as f64 * stepv;
                }
            }
        }
    }

    /// CheckConnectivity (cxx:675-758): 23 samples along each shared boundary
    /// and myUClosed / myVClosed set from the first wrap pair.
    pub fn check_connectivity(&mut self, prec: f64) -> bool {
        const NPOINTS: usize = 23;
        let nb_u = self.nb_u_patches();
        let nb_v = self.nb_v_patches();
        let mut ok = true;
        if nb_u == 0 || nb_v == 0 {
            return false;
        }
        // u direction
        let mut j = nb_u;
        for i in 1..=nb_u {
            let mut maxdist2 = 0.0f64;
            for k in 1..=nb_v {
                let (Some(sj), Some(si)) = (self.patch(j, k), self.patch(i, k)) else {
                    continue;
                };
                let (_, uj2, vj1, vj2) = (
                    limit_value(sj.u_range().0),
                    limit_value(sj.u_range().1),
                    limit_value(sj.v_range().0),
                    limit_value(sj.v_range().1),
                );
                let (ui1, _, vi1, vi2) = (
                    limit_value(si.u_range().0),
                    limit_value(si.u_range().1),
                    limit_value(si.v_range().0),
                    limit_value(si.v_range().1),
                );
                let stepj = (vj2 - vj1) / (NPOINTS as f64 - 1.0);
                let stepi = (vi2 - vi1) / (NPOINTS as f64 - 1.0);
                for s in 0..NPOINTS {
                    let parj = vj1 + stepj * s as f64;
                    let pari = vi1 + stepi * s as f64;
                    let d2 = sj.d0(uj2, parj).distance(&si.d0(ui1, pari));
                    let d2 = d2 * d2;
                    if maxdist2 < d2 {
                        maxdist2 = d2;
                    }
                }
            }
            if i == 1 {
                self.u_closed = maxdist2 <= prec * prec;
            } else if maxdist2 > prec * prec {
                ok = false;
            }
            j = i;
        }
        // v direction
        let mut j = nb_v;
        for i in 1..=nb_v {
            let mut maxdist2 = 0.0f64;
            for k in 1..=nb_u {
                let (Some(sj), Some(si)) = (self.patch(k, j), self.patch(k, i)) else {
                    continue;
                };
                let (uj1, uj2, _, vj2) = (
                    limit_value(sj.u_range().0),
                    limit_value(sj.u_range().1),
                    limit_value(sj.v_range().0),
                    limit_value(sj.v_range().1),
                );
                let (ui1, ui2, vi1, _) = (
                    limit_value(si.u_range().0),
                    limit_value(si.u_range().1),
                    limit_value(si.v_range().0),
                    limit_value(si.v_range().1),
                );
                let stepj = (uj2 - uj1) / (NPOINTS as f64 - 1.0);
                let stepi = (ui2 - ui1) / (NPOINTS as f64 - 1.0);
                for s in 0..NPOINTS {
                    let parj = uj1 + stepj * s as f64;
                    let pari = ui1 + stepi * s as f64;
                    let d2 = sj.d0(parj, vj2).distance(&si.d0(pari, vi1));
                    let d2 = d2 * d2;
                    if maxdist2 < d2 {
                        maxdist2 = d2;
                    }
                }
            }
            if i == 1 {
                self.v_closed = maxdist2 <= prec * prec;
            } else if maxdist2 > prec * prec {
                ok = false;
            }
            j = i;
        }
        ok
    }
}
