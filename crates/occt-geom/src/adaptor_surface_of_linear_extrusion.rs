//! `GeomAdaptor_SurfaceOfLinearExtrusion`, the part `Draft_Modification`
//! consumes: the surface-type classification and the `Plane` / `Cylinder`
//! conversions.
//!
//! Source: `ModelingData/TKG3d/GeomAdaptor/GeomAdaptor_SurfaceOfLinearExtrusion.cxx`
//! - `GetType` (`cxx:269-328`), `Plane` (`cxx:332-371`) and `Cylinder`
//! (`cxx:375-387`).
//!
//! # Translation notes
//!
//! * `GetType` tests `myBasisCurve->GetType()`, i.e. the `GeomAbs_CurveType`
//!   of the basis adaptor, not its dynamic type; the port uses the
//!   [`Curve`] accessors that answer the same classification
//!   (`gp_line`, `gp_circ`, `gp_ellipse`, `gp_parabola`, `gp_hyperbola`).
//! * `gp_Ax3 Ax3(P, gp_Dir(newZ), gp_Dir(D1u))` is
//!   [`GpAx3::from_point_dir_x`], which raises on a null vector the way
//!   `gp_Dir` does; both `Plane` raise sites become `None`.
//! * `Plane` probes at most 21 points whenever the basis range is infinite
//!   (`cxx:341-354`); the loop breaks on the first non-degenerate
//!   `D1u.Normalized().Crossed(myDirection)` (`cxx:356-363`).

use occt_core::gp::{GpAx3, GpCylinder, GpDir, GpPln, GpPnt, GpVec};
use occt_core::kernel::geomabs::SurfaceType;
use occt_core::precision::{Precision, ANGULAR};

use crate::curve::Curve;

/// `GeomAdaptor_SurfaceOfLinearExtrusion`: a basis curve swept along a
/// direction.
pub struct GeomAdaptorSurfaceOfLinearExtrusion {
    basis: std::sync::Arc<dyn Curve>,
    direction: GpDir,
}

impl GeomAdaptorSurfaceOfLinearExtrusion {
    /// `GeomAdaptor_SurfaceOfLinearExtrusion(C, D)`.
    pub fn new(basis: std::sync::Arc<dyn Curve>, direction: GpDir) -> Self {
        Self { basis, direction }
    }

    /// `Load(C, D)`.
    pub fn load(&mut self, basis: std::sync::Arc<dyn Curve>, direction: GpDir) {
        self.basis = basis;
        self.direction = direction;
    }

    /// `BasisCurve()`'s `Adaptor3d_Curve`.
    pub fn curve(&self) -> &std::sync::Arc<dyn Curve> {
        &self.basis
    }

    /// `GetType()` (`cxx:269-328`): `GeomAbs_Plane` / `GeomAbs_Cylinder` when
    /// the swept basis degenerates into one, `GeomAbs_SurfaceOfExtrusion`
    /// otherwise.
    pub fn get_type(&self) -> SurfaceType {
        // The port's `Curve` accessors stand in for `Adaptor3d_Curve::GetType`.
        if let Some(line) = self.basis.gp_line() {
            // `cxx:273-281`
            if !self.direction.is_parallel_tol(&line.direction(), ANGULAR) {
                return SurfaceType::Plane;
            }
        } else if let Some(circ) = self.basis.gp_circ() {
            // `cxx:283-294`
            let circ_axis = circ.axis();
            let d = circ_axis.direction();
            if self.direction.is_parallel_tol(d, ANGULAR) {
                return SurfaceType::Cylinder;
            } else if self.direction.is_normal_tol(d, ANGULAR) {
                return SurfaceType::Plane;
            }
        } else if let Some(elips) = self.basis.gp_ellipse() {
            // `cxx:296-303`
            if self.direction.is_normal_tol(&elips.position().direction(), ANGULAR) {
                return SurfaceType::Plane;
            }
        } else if let Some(parab) = self.basis.gp_parabola() {
            // `cxx:305-312`
            if self.direction.is_normal_tol(&parab.position().direction(), ANGULAR) {
                return SurfaceType::Plane;
            }
        } else if let Some(hypr) = self.basis.gp_hyperbola() {
            // `cxx:314-321`
            if self.direction.is_normal_tol(&hypr.position().direction(), ANGULAR) {
                return SurfaceType::Plane;
            }
        }
        // `cxx:327`: `GeomAbs_SurfaceOfExtrusion`.
        SurfaceType::SurfaceOfExtrusion
    }

    /// `Plane()` (`cxx:332-371`): `None` where the reference raises
    /// `Standard_NoSuchObject`.
    pub fn plane(&self) -> Option<GpPln> {
        if self.get_type() != SurfaceType::Plane {
            return None;
        }
        let mut u_first = self.basis.first_parameter();
        let mut u_last = self.basis.last_parameter();
        // `cxx:341-354`
        if Precision::is_negative_infinite(u_first) && Precision::is_positive_infinite(u_last) {
            u_first = -100.0;
            u_last = 100.0;
        } else if Precision::is_negative_infinite(u_first) {
            u_first = u_last - 200.0;
        } else if Precision::is_positive_infinite(u_last) {
            u_last = u_first + 200.0;
        }
        let deltau = (u_last - u_first) / 20.0;
        let mut p = GpPnt::new(0.0, 0.0, 0.0);
        let mut d1u = GpVec::zero();
        let mut new_z = GpVec::zero();
        // `cxx:355-364`
        for i in 1..=21 {
            let prm = u_first + (i - 1) as f64 * deltau;
            let (pt, du) = self.basis.d1(prm);
            p = pt;
            d1u = du;
            new_z = du.normalized().crossed(&GpVec::from_xyz(self.direction.xyz()));
            if new_z.magnitude() > 1.0e-12 {
                break;
            }
        }
        // `cxx:365-370`
        let z = GpDir::from_vec(&new_z).ok()?;
        let vx = GpDir::from_vec(&d1u).ok()?;
        let mut ax3 = GpAx3::from_point_dir_x(p, z, &vx).ok()?;
        if self.direction.dot(ax3.y_direction()) < 0.0 {
            ax3.y_reverse();
        }
        Some(GpPln::new(ax3))
    }

    /// `Cylinder()` (`cxx:375-387`): `None` where the reference raises
    /// `Standard_NoSuchObject`.
    pub fn cylinder(&self) -> Option<GpCylinder> {
        if self.get_type() != SurfaceType::Cylinder {
            return None;
        }
        let circ = self.basis.gp_circ()?;
        // `cxx:382-385`
        let mut ax3 = GpAx3::from_ax2(&circ.position());
        let circ_axis = circ.axis();
        if self.direction.dot(circ_axis.direction()) < 0.0 {
            ax3.z_reverse();
        }
        GpCylinder::new(ax3, circ.radius()).ok()
    }
}
