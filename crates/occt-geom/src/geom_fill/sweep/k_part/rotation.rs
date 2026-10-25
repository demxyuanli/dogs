//! `GeomFill_Sweep::BuildKPart`, branch (2): a circular trajectory
//! (`GeomFill_Sweep.cxx:802-1112`, OCCT 8.0.0).
//!
//! `(2.1)` a circular section (`cxx:867-1001`) recognises a sphere (`2.1.a`)
//! or a torus (`2.1.b`); `(2.2)` a linear section (`cxx:1003-1098`) recognises
//! a cylinder (`2.2.a`) or a cone (`2.2.b`); `(2.3)` otherwise emits a surface
//! of revolution (`cxx:1100-1111`).
//!
//! Translation notes:
//! - The commented-out `try { Tf2.SetValues(...) } catch (...)` of `cxx:836-853`
//!   leaves a `Standard_ConstructionError` from `gp_Trsf::SetValues`
//!   uncaught in OCCT, i.e. `BuildKPart` unwinds. Rust has no exceptions, so
//!   the port returns "no KPart" instead of aborting; this is only reachable
//!   for a degenerate (null determinant) location matrix.
//! - `Section->Transform(Tf2)` mutates the law's own curve in OCCT
//!   (`cxx:869`, `cxx:1105`); the port copies the handle first.
//! - `ElSLib::Parameters(theSphere, P, U, V)` (`cxx:920-921`) is
//!   `slib::sphere_parameters(sphere.position(), &p)`.

use std::sync::Arc;

use occt_core::elib::slib;
use occt_core::gp::{GpAx1, GpCone, GpCylinder, GpPnt, GpSphere, GpTorus, GpVec};
use occt_core::kernel::geomabs::CurveType;

use super::super::adaptor_curve_type;
use super::state::{
    ax3_from_point_dir_x, dir_of, occt_infinite, trsf_from_matrix, KPartState, GP_RESOLUTION,
    P_CONFUSION,
};
use super::Sweep;
use crate::cone::GeomCone;
use crate::curve::Curve;
use crate::cylinder::GeomCylinder;
use crate::line::GeomLine;
use crate::sphere::GeomSphere;
use crate::surface::Surface;
use crate::surface_of_revolution::GeomSurfaceOfRevolution;
use crate::torus::GeomTorus;
use crate::trimmed::GeomTrimmedCurve;

impl Sweep {
    /// `GeomFill_Sweep::BuildKPart` branch (2) (`cxx:802-1112`).
    ///
    /// `None` means "carry on"; `Some(v)` is the early `return v` of
    /// `BuildKPart` (`cxx:826` returns `false`, `cxx:945` returns `Ok`).
    pub(super) fn kpart_rotation(&mut self, st: &mut KPartState) -> Option<bool> {
        let sec = match self.sec.clone() {
            Some(s) => s,
            None => return Some(false),
        };
        let mut error = st.error;
        if !sec.is_constant(&mut error) {
            return None;
        }
        st.error = error;

        // `cxx:808-810`: the trajectory.
        st.is_v_periodic = (self.last - self.first - 2.0 * std::f64::consts::PI).abs() < 1.0e-15;
        let (m, ds_param) = match self.loc.d0(0.1) {
            Ok(v) => v,
            Err(_) => return Some(false),
        };
        let (_, v) = match self.loc.d0(0.0) {
            Ok(v) => v,
            Err(_) => return Some(false),
        };
        let mut centre = match self.loc.rotation() {
            Ok(c) => c,
            Err(_) => return Some(false),
        };

        let mut dp = ds_param.subtracted(&v);
        let mut ds = GpVec::from_xyz(&v.xyz().subtracted(centre.xyz()));
        let rot_radius = ds.magnitude();
        if rot_radius > 1.0e-15 {
            ds.normalize();
        } else {
            // `cxx:822-826`: a degenerate rotation is not a KPart.
            return Some(false);
        }
        // `cxx:828-830`.
        let mut dn = ds.crossed(&dp);
        dn.normalize();
        dp = dn.crossed(&ds);
        dp.normalize();

        let tf2 = match trsf_from_matrix(&m, &v) {
            Ok(t) => t,
            Err(_) => return Some(false),
        };

        // `cxx:854-862`: the section.
        let section = match sec.constant_section() {
            Ok(c) => c,
            Err(_) => return Some(false),
        };
        let section_type = adaptor_curve_type(&section);
        st.u_first = occt_infinite(section.first_parameter());
        st.u_last = occt_infinite(section.last_parameter());

        // `(2.1)` `cxx:867-1001`: torus or sphere.
        if section_type == CurveType::Circle {
            let mut circ = match section.gp_circ() {
                Some(c) => c,
                None => return Some(false),
            };
            let mut is_good_side = true;
            circ.transform(&tf2);

            // `cxx:871-886`: cut the section plane through the rotation centre.
            let mut dc = GpVec::from_xyz(&circ.location().xyz().subtracted(centre.xyz()));
            let shift_along_dn = dn.multiplied_scalar(dc.dot(&dn));
            centre.translate_vec(&shift_along_dn);
            dc = GpVec::from_xyz(&circ.location().xyz().subtracted(centre.xyz()));
            let radius = dc.magnitude();
            if radius > st.tol && dc.dot(&ds) < 0.0 {
                is_good_side = false;
            }
            if radius < st.tol / 100.0 {
                dc = ds;
            }

            // `cxx:887-900`: the error of the section plane rotation.
            let mut nc = GpVec::from_xyz(circ.position().direction().xyz());
            nc.normalize();
            error = nc.dot(&dn).abs();
            error += nc.dot(&ds).abs();
            error *= circ.radius();
            if error <= st.tol {
                st.s_error = error;
                error += radius;
                if radius <= st.tol {
                    // `(2.1.a)` `cxx:901-957`: sphere.
                    //
                    // `cxx:902`: `f` / `l` start at the section bounds and are
                    // then overwritten by the V parameters of the trimmed
                    // section on the sphere (both arms below assign).
                    let f;
                    let l;
                    st.s_error = error;
                    // `cxx:905`: `Centre.BaryCenter(1.0, C.Location(), 1.0)`.
                    let centre_sum = GpPnt::new(
                        (centre.x() + circ.location().x()) / 2.0,
                        (centre.y() + circ.location().y()) / 2.0,
                        (centre.z() + circ.location().z()) / 2.0,
                    );
                    let dn_dir = match dir_of(&dn) {
                        Ok(d) => d,
                        Err(_) => return Some(false),
                    };
                    let ds_dir = match dir_of(&ds) {
                        Ok(d) => d,
                        Err(_) => return Some(false),
                    };
                    // `gp_Ax3 AxisOfSphere(Centre, DN, DS)` (`cxx:906`).
                    let axis = match ax3_from_point_dir_x(centre_sum, dn_dir, &ds_dir) {
                        Ok(a) => a,
                        Err(_) => return Some(false),
                    };
                    let a_radius = circ.radius();
                    let sphere = match GpSphere::new(axis, a_radius) {
                        Ok(s) => s,
                        Err(_) => return Some(false),
                    };
                    st.s = Some(Arc::new(GeomSphere::new(sphere)));

                    // `cxx:912-933`: the U window follows the trimmed section
                    // because V cannot be controlled (`myExchUV = true`).
                    let fpar = section.first_parameter();
                    let lpar = section.last_parameter();
                    let mut the_section =
                        GeomTrimmedCurve::new(section.clone(), fpar, lpar);
                    the_section.transform(&tf2);
                    let first_point = the_section.d0(the_section.first_parameter());
                    let last_point = the_section.d0(the_section.last_parameter());
                    let (u_first_on_sec, v_first_on_sec) =
                        slib::sphere_parameters(sphere.position(), &first_point);
                    let (u_last_on_sec, v_last_on_sec) =
                        slib::sphere_parameters(sphere.position(), &last_point);
                    if v_first_on_sec < v_last_on_sec {
                        f = v_first_on_sec;
                        l = v_last_on_sec;
                    } else {
                        // `cxx:929-932`: the parametric orientation is reversed.
                        f = v_last_on_sec;
                        l = v_first_on_sec;
                        st.u_reversed = true;
                    }

                    if (l - f).abs() <= P_CONFUSION
                        || (u_last_on_sec - u_first_on_sec).abs()
                            > std::f64::consts::FRAC_PI_2
                    {
                        // `cxx:937-944`: a "degenerated" or "twisted" patch that
                        // a trimmed sphere cannot represent: `return Ok`.
                        st.u_reversed = false;
                        return Some(st.ok);
                    }

                    if f >= -std::f64::consts::FRAC_PI_2
                        && l <= std::f64::consts::FRAC_PI_2
                    {
                        st.ok = true;
                        st.exch_uv = true;
                        st.u_first = f;
                        st.u_last = l;
                    } else {
                        // `cxx:953-956`: restore what is needed.
                        st.u_reversed = false;
                    }
                } else if is_good_side {
                    // `(2.1.b)` `cxx:957-999`: torus.
                    let dn_dir = match dir_of(&dn) {
                        Ok(d) => d,
                        Err(_) => return Some(false),
                    };
                    let dc_dir = match dir_of(&dc) {
                        Ok(d) => d,
                        Err(_) => return Some(false),
                    };
                    // `gp_Ax3 AxisOfTore(Centre, DN, DC)` (`cxx:960`).
                    let axis = match ax3_from_point_dir_x(centre, dn_dir, &dc_dir) {
                        Ok(a) => a,
                        Err(_) => return Some(false),
                    };
                    let torus = match GpTorus::new(axis, radius, circ.radius()) {
                        Ok(t) => t,
                        Err(_) => return Some(false),
                    };
                    let surface: Arc<dyn Surface> = Arc::new(GeomTorus::new(torus));
                    st.s = Some(surface.clone());

                    // `cxx:966-970`: the isoparametric circle at U = 0.
                    let iso = match surface.u_iso_curve(0.0).and_then(|c| c.gp_circ()) {
                        Some(c) => c,
                        None => return Some(false),
                    };
                    let axeiso = iso.position();
                    if circ
                        .position()
                        .direction()
                        .is_opposite_ang(&axeiso.direction(), 0.1)
                    {
                        // `cxx:971-980`: the parametric orientation is reversed.
                        let l = 2.0 * std::f64::consts::PI - st.u_first;
                        let f = 2.0 * std::f64::consts::PI - st.u_last;
                        st.u_first = f;
                        st.u_last = l;
                        st.u_reversed = true;
                    }

                    // `cxx:981-986`: the angular slide.
                    let rot = circ.position().x_direction().angle_with_ref(
                        axeiso.x_direction(),
                        &axeiso.direction(),
                    );
                    st.u_first -= rot;
                    st.u_last -= rot;

                    st.exch_uv = true;
                    // `cxx:990-993`: the seam may not be where the periodic
                    // window expects it.
                    if st.is_u_periodic && st.u_first.abs() > P_CONFUSION {
                        st.is_u_periodic = false;
                    }
                    st.ok = true;
                }
            } else {
                st.s_error = 0.0;
            }
        }

        // `(2.2)` `cxx:1003-1098`: cylinder or cone.
        if section_type == CurveType::Line {
            let mut line = match section.gp_line() {
                Some(l) => l,
                None => return Some(false),
            };
            line.transform(&tf2);
            let dl = GpVec::from_xyz(line.direction().xyz());
            let levier = section
                .first_parameter()
                .abs()
                .max(section.last_parameter());
            st.s_error = error + levier * dl.dot(&dp).abs();
            if st.s_error <= st.tol {
                // `cxx:1015-1018`: the parametrisation direction.
                let reverse = dl.dot(&dn) < 0.0;
                let dn_dir = match dir_of(&dn) {
                    Ok(d) => d,
                    Err(_) => return Some(false),
                };
                // `gp_Lin Dir(Centre, DN)` (`cxx:1016`).
                let dir_line = occt_core::gp::GpLin::new(GpAx1::new(centre, dn_dir));
                // `cxx:1020-1026`: the centre of the surface on the axis.
                let o1o2 = GpVec::from_pnts(&centre, &line.location());
                let trans = dn.multiplied_scalar(dn.dot(&o1o2));
                let centre_of_surf = centre.translated_vec(&trans);
                let ds_on_axis = GpVec::from_xyz(
                    &line.location().xyz().subtracted(centre_of_surf.xyz()),
                );

                error = st.s_error;
                error += dl.cross_magnitude(&dn) * levier;
                if error <= st.tol {
                    // `(2.2.a)` `cxx:1031-1045`: cylinder.
                    st.s_error = error;
                    let mut axis = match occt_core::gp::GpAx3::from_point_dir(
                        centre_of_surf,
                        dn_dir,
                    ) {
                        Ok(a) => a,
                        Err(_) => return Some(false),
                    };
                    if ds_on_axis.square_magnitude() > GP_RESOLUTION {
                        let x_dir = match dir_of(&ds_on_axis) {
                            Ok(d) => d,
                            Err(_) => return Some(false),
                        };
                        axis.set_x_direction_projected(&x_dir);
                    }
                    let cylinder = match GpCylinder::new(axis, line.distance(&centre_of_surf)) {
                        Ok(c) => c,
                        Err(_) => return Some(false),
                    };
                    st.s = Some(Arc::new(GeomCylinder::new(cylinder)));
                    st.ok = true;
                    st.exch_uv = true;
                } else {
                    // `cxx:1047-1075`: the half angle of the cone.
                    let mut angle = dir_line.angle(&line).abs();
                    if angle > std::f64::consts::FRAC_PI_2 {
                        angle = std::f64::consts::PI - angle;
                    }
                    if reverse {
                        angle = -angle;
                    }
                    let aux = ds_on_axis.dot(&dl);
                    if aux < 0.0 {
                        angle = -angle;
                    }
                    if (angle.abs() - std::f64::consts::FRAC_PI_2).abs() > 0.01 {
                        // `(2.2.b)` `cxx:1062-1070`: cone.
                        let radius = centre_of_surf.distance(&line.location());
                        let ds_dir = match dir_of(&ds_on_axis) {
                            Ok(d) => d,
                            Err(_) => return Some(false),
                        };
                        // `gp_Ax3 Axis(CentreOfSurf, Dir.Direction(), DS)`
                        // (`cxx:1067`).
                        let axis = match ax3_from_point_dir_x(centre_of_surf, dn_dir, &ds_dir) {
                            Ok(a) => a,
                            Err(_) => return Some(false),
                        };
                        let cone = match GpCone::new(axis, radius, angle) {
                            Ok(c) => c,
                            Err(_) => return Some(false),
                        };
                        st.s = Some(Arc::new(GeomCone::new(cone)));
                        st.exch_uv = true;
                        st.ok = true;
                    } else {
                        // `cxx:1072-1075`: undecided, the error goes back to 0.
                        st.s_error = 0.0;
                    }
                }

                if st.ok && reverse {
                    // `cxx:1079-1092`: reverse the U parameter.
                    let cl = GeomLine::new(line);
                    let uf = cl.reversed_parameter(st.u_last);
                    let ul = cl.reversed_parameter(st.u_first);
                    st.u_first = ul;
                    st.u_last = uf;
                }
            } else {
                st.s_error = 0.0;
            }
        }

        // `(2.3)` `cxx:1100-1111`: revolution.
        if !st.ok {
            let mut transformed = section.clone_dyn();
            transformed.transform(&tf2);
            let dn_dir = match dir_of(&dn) {
                Ok(d) => d,
                Err(_) => return Some(false),
            };
            let axis = GpAx1::new(centre, dn_dir);
            st.s = Some(Arc::new(GeomSurfaceOfRevolution::new(
                Arc::from(transformed),
                axis,
            )));
            st.exch_uv = true;
            st.s_error = 0.0;
            st.ok = true;
        }

        None
    }
}
