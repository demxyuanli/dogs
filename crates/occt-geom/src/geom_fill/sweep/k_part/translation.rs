//! `GeomFill_Sweep::BuildKPart`, branch (1): a rectilinear trajectory
//! (`GeomFill_Sweep.cxx:594-798`, OCCT 8.0.0).
//!
//! `(1.1)` constant section (`cxx:633-715`) recognises a plane (`1.1.a`) or a
//! cylinder (`1.1.b`) and otherwise emits a surface of linear extrusion
//! (`1.1.c`); `(1.2)` a conical law (`cxx:717-798`) recognises a cone.
//!
//! Translation notes:
//! - `Section->Transform(Tf2)` / `Section->Translate(...)` (`cxx:705`,
//!   `cxx:729-730`, `cxx:737`, `cxx:746`) mutate the handle
//!   `ConstantSection` / `CirclSection` return, which in OCCT is the law's own
//!   curve. The port copies it first (`clone_dyn`) so the law is not mutated;
//!   the produced surface is the same.
//! - `AC.GetType()` / `AC.FirstParameter()` / `AC.LastParameter()` are the
//!   `GeomAdaptor_Curve` of that (untransformed) section:
//!   [`adaptor_curve_type`], `first_parameter`, `last_parameter`.
//! - `AC.D1(0, P2, dsection)` (`cxx:735`) is `d1(0)`.
//! - `gp_Ax2(P, DS ^ DP, DS)` relies on the implicit `gp_Dir(const gp_Vec&)`
//!   (`gp_Dir.hxx:338-343`), hence [`dir_of`].

use std::sync::Arc;

use occt_core::gp::{GpAx2, GpAx3, GpCone, GpCylinder, GpPln, GpTrsf, GpVec};
use occt_core::kernel::geomabs::CurveType;

use super::super::adaptor_curve_type;
use super::state::{ax3_from_point_dir_x, dir_of, occt_infinite, trsf_from_matrix, KPartState};
use super::{is_sweep_parallel_spine, Sweep};
use crate::cone::GeomCone;
use crate::cylinder::GeomCylinder;
use crate::geom_fill::section_law::SectionLaw;
use crate::plane::GeomPlane;
use crate::surface_of_linear_extrusion::GeomSurfaceOfLinearExtrusion;

/// `TolProd` of the cylindrical test (`GeomFill_Sweep.cxx:670`).
const TOL_PROD: f64 = 1.0e-6;

impl Sweep {
    /// `GeomFill_Sweep::BuildKPart` branch (1) (`cxx:594-798`).
    ///
    /// `None` means "carry on with branch (2)"; `Some(v)` is the early
    /// `return v` of `BuildKPart`.
    pub(super) fn kpart_translation(&mut self, st: &mut KPartState) -> Option<bool> {
        let sec = match self.sec.clone() {
            Some(s) => s,
            None => return Some(false),
        };

        // `cxx:596-601`: the translation of the location law.
        let (m, ds) = match self.loc.d0(1.0) {
            Ok(v) => v,
            Err(_) => return Some(false),
        };
        let (_, v) = match self.loc.d0(0.0) {
            Ok(v) => v,
            Err(_) => return Some(false),
        };
        let mut dp = ds.subtracted(&v);
        dp.normalize();

        // `cxx:602-628`: the affine transform of the law, rejected when
        // `gp_Trsf::SetValues` raises (`IsTrsf = false; return false`).
        let tf2 = match trsf_from_matrix(&m, &v) {
            Ok(t) => t,
            Err(_) => return Some(false),
        };

        let mut error = st.error;

        if sec.is_constant(&mut error) {
            // `(1.1)` `cxx:633-715`: extrusion.
            st.error = error;
            let section = match sec.constant_section() {
                Ok(c) => c,
                Err(_) => return Some(false),
            };
            let section_type = adaptor_curve_type(&section);
            st.u_first = occt_infinite(section.first_parameter());
            st.u_last = occt_infinite(section.last_parameter());

            // `(1.1.a)` `cxx:641-665`: plane.
            if section_type == CurveType::Line {
                // `cxx:644-648` (OCC5073): the spine has to be parallel to the
                // first-to-last section translation.
                if !is_sweep_parallel_spine(&self.loc, sec.as_ref(), st.tol) {
                    return Some(false);
                }
                let mut line = match section.gp_line() {
                    Some(l) => l,
                    None => return Some(false),
                };
                line.transform(&tf2);
                let mut spine_dir = GpVec::from_xyz(line.position().direction().xyz());
                spine_dir.normalize();
                let levier = spine_dir.dot(&dp).abs();
                st.s_error = error + levier * (self.last - self.first).abs();
                if st.s_error <= st.tol {
                    st.ok = true;
                    let normal = match dir_of(&spine_dir.crossed(&dp)) {
                        Ok(d) => d,
                        Err(_) => return Some(false),
                    };
                    let plane_x = match dir_of(&spine_dir) {
                        Ok(d) => d,
                        Err(_) => return Some(false),
                    };
                    let axis_of_plane = match GpAx2::new(line.location(), normal, plane_x) {
                        Ok(a) => a,
                        Err(_) => return Some(false),
                    };
                    st.s = Some(Arc::new(GeomPlane::new(GpPln::new(
                        GpAx3::from_ax2(&axis_of_plane),
                    ))));
                } else {
                    st.s_error = 0.0;
                }
            }

            // `(1.1.b)` `cxx:668-699`: cylinder.
            if section_type == CurveType::Circle {
                let mut circ = match section.gp_circ() {
                    Some(c) => c,
                    None => return Some(false),
                };
                circ.transform(&tf2);

                let mut axis_dir = GpVec::from_xyz(circ.position().direction().xyz());
                axis_dir.normalize();
                let levier = axis_dir.cross_magnitude(&dp) * circ.radius();
                st.s_error = levier * (self.last - self.first).abs();
                if st.s_error <= TOL_PROD {
                    st.ok = true;
                    // `gp_Ax3 axe(C.Location(), DP, C.Position().XDirection())`
                    // (`cxx:681`).
                    let dp_dir = match dir_of(&dp) {
                        Ok(d) => d,
                        Err(_) => return Some(false),
                    };
                    let axe = match ax3_from_point_dir_x(
                        circ.location(),
                        dp_dir,
                        circ.position().x_direction(),
                    ) {
                        Ok(a) => a,
                        Err(_) => return Some(false),
                    };
                    let cylinder = match GpCylinder::new(axe, circ.radius()) {
                        Ok(c) => c,
                        Err(_) => return Some(false),
                    };
                    st.s = Some(Arc::new(GeomCylinder::new(cylinder)));
                    if circ
                        .position()
                        .direction()
                        .is_opposite_ang(&axe.direction(), 0.1)
                    {
                        // `cxx:686-694`: the parametric orientation is reversed.
                        let l = 2.0 * std::f64::consts::PI - st.u_first;
                        let f = 2.0 * std::f64::consts::PI - st.u_last;
                        st.u_first = f;
                        st.u_last = l;
                        st.u_reversed = true;
                    }
                } else {
                    st.s_error = 0.0;
                }
            }

            // `(1.1.c)` `cxx:701-715`: a plain extrusion. The `IsTrsf = false`
            // arm of `cxx:711-714` is unreachable here because `(1.1)` already
            // returned.
            if !st.ok {
                let mut transformed = section.clone_dyn();
                transformed.transform(&tf2);
                let dp_dir = match dir_of(&dp) {
                    Ok(d) => d,
                    Err(_) => return Some(false),
                };
                st.s = Some(Arc::new(GeomSurfaceOfLinearExtrusion::new(
                    Arc::from(transformed),
                    dp_dir,
                )));
                st.s_error = 0.0;
                st.ok = true;
            }
        } else {
            // `(1.2)` `cxx:717-798`: conical law.
            let mut cone_error = st.error;
            if sec.is_conical_law(&mut cone_error) {
                if let Some(result) = self.kpart_cone(st, &sec, &tf2, &dp, cone_error) {
                    return Some(result);
                }
            }
        }

        None
    }

    /// `(1.2)` `GeomFill_Sweep.cxx:720-797`: the cone of a conical section law.
    fn kpart_cone(
        &mut self,
        st: &mut KPartState,
        sec: &Arc<dyn SectionLaw>,
        tf2: &GpTrsf,
        dp: &GpVec,
        error: f64,
    ) -> Option<bool> {
        // `cxx:728-736`: the last circle. `Section` is copied before the
        // `Transform` / `Translate` calls of `cxx:729-730`.
        let mut section = match sec.circl_section(self.s_last) {
            Ok(c) => c.clone_dyn(),
            Err(_) => return Some(false),
        };
        section.transform(tf2);
        let mut shift = GpTrsf::identity();
        shift.set_translation_vec(&dp.multiplied_scalar(self.last));
        section.transform(&shift);
        let circ_last = match section.gp_circ() {
            Some(c) => c,
            None => return Some(false),
        };
        let centre2 = circ_last.location();
        let (p2, dsection) = section.d1(0.0);
        let r2 = circ_last.radius();

        // `cxx:737-743`: the first circle.
        let mut section = match sec.circl_section(self.s_first) {
            Ok(c) => c.clone_dyn(),
            Err(_) => return Some(false),
        };
        section.transform(tf2);
        let mut shift = GpTrsf::identity();
        shift.set_translation_vec(&dp.multiplied_scalar(self.first));
        section.transform(&shift);
        let circ_first = match section.gp_circ() {
            Some(c) => c,
            None => return Some(false),
        };
        let centre1 = circ_first.location();
        let p1 = section.d0(0.0);
        let r1 = circ_first.radius();

        // `cxx:745-750`: the circle at the apex parameter, untransformed.
        let apex_param =
            self.s_first - self.first * (self.s_last - self.s_first) / (self.last - self.first);
        let section = match sec.circl_section(apex_param) {
            Ok(c) => c,
            Err(_) => return Some(false),
        };
        let circ_apex = match section.gp_circ() {
            Some(c) => c,
            None => return Some(false),
        };
        let centre0 = circ_apex.location();

        // `cxx:752-769`: the half angle is only accepted strictly inside
        // `(0.01, PI / 2 - 0.01)`.
        let mut axis_vec = GpVec::from_pnts(&centre1, &p1);
        if axis_vec.magnitude() < 1.0e-9 {
            axis_vec = GpVec::from_pnts(&centre2, &p2);
        }
        let line_vec = GpVec::from_pnts(&p1, &p2);
        let dir_vec = GpVec::from_pnts(&centre1, &centre2);

        let mut angle = line_vec.angle(&dir_vec);
        if !(angle > 0.01 && angle < std::f64::consts::FRAC_PI_2 - 0.01) {
            return None;
        }
        if r2 < r1 {
            angle = -angle;
        }
        st.s_error = error;
        let dir = match dir_of(&dir_vec) {
            Ok(d) => d,
            Err(_) => return Some(false),
        };
        let normal = match dir_of(&axis_vec) {
            Ok(d) => d,
            Err(_) => return Some(false),
        };
        // `gp_Ax3 Axis(Centre0, Dir, N)` (`cxx:769`): the surface axis follows
        // the centres, the cone `ReferencePlane` follows `N`.
        let axis = match ax3_from_point_dir_x(centre0, dir, &normal) {
            Ok(a) => a,
            Err(_) => return Some(false),
        };
        let cone = match GpCone::new(axis, circ_apex.radius(), angle) {
            Ok(c) => c,
            Err(_) => return Some(false),
        };
        st.s = Some(Arc::new(GeomCone::new(cone)));

        // `cxx:771-778`: the parametric slide and the U bounds.
        st.v_first = self.first / angle.cos();
        st.v_last = self.last / angle.cos();
        st.u_first = occt_infinite(section.first_parameter());
        st.u_last = occt_infinite(section.last_parameter());

        // `cxx:779-786`: `S->VIso(VLast)->D1(0, pbis, diso)`.
        if let Some(iso) = st.s.as_ref().and_then(|s| s.v_iso_curve(st.v_last)) {
            let (_, diso) = iso.d1(0.0);
            if diso.magnitude() > 1.0e-9 && dsection.magnitude() > 1.0e-9 {
                st.u_reversed = diso.is_opposite_ang(&dsection, 0.1);
            }
        }
        if st.u_reversed {
            // `cxx:787-793`: the parametric orientation is reversed.
            let l = 2.0 * std::f64::consts::PI - st.u_first;
            let f = 2.0 * std::f64::consts::PI - st.u_last;
            st.u_first = f;
            st.u_last = l;
        }

        st.ok = true;
        None
    }
}
