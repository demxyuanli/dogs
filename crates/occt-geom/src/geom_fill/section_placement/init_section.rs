//! Frame and `Gabarit` initialization of `GeomFill_SectionPlacement`
//! (`GeomFill_SectionPlacement.cxx:157-380`), plus the `Geom_BSplineCurve::
//! LocateU` query the mean-plane sampling needs (`cxx:285-286`).

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::bspl::locate::hunt_occt;
use occt_core::gp::{GpDir, GpPnt};
use occt_core::precision::{
    COMPUTATIONAL, CONFUSION, PCONFUSION,
};

use crate::curve::Curve;
use crate::extrema_pc::{extrema_ext_pc_range_done, ExtPcSolution};
use crate::geom_bnd_lib::box_curve;
use crate::geom_lib::geom_lib_axe_of_inertia;

use super::math_helpers::tangente;
use super::SectionPlacement;

impl SectionPlacement {
    /// Constructor tail: `Gabarit` from the section box (`cxx:182-204`),
    /// `TheAxe` from the section tangent (`cxx:206-212`) or from the conic axis
    /// (`cxx:214-240`), then the mean plane when the section is not planar
    /// (`cxx:252-371`).
    pub(super) fn init_section(&mut self) {
        // Boite d'encombrement de la section pour en deduire le gabarit
        let mut box_ = BndBox::new();
        if self.is_point {
            box_.add_point(&self.point);
        } else {
            let section = self
                .section_curve
                .as_ref()
                .expect("GeomFill_SectionPlacement: null section");
            // `BndLib_Add3dCurve::Add(myAdpSection, 1.0e-4, box)` (`cxx:192`).
            box_.add_box(&box_curve(
                section.as_ref(),
                section.first_parameter(),
                section.last_parameter(),
                1.0e-4,
            ));
        }
        // `box.Get(...)` (`cxx:194`): `Bnd_Box::Get` throws for a void box
        // (`Bnd_Box.cxx:214-217`). An open bound is `Precision::Infinite()`
        // there and +/-inf in `occt_core::bnd::BndBox::get`; only the fact that
        // the diagonal is infinite reaches `Gabarit`, so the value is the same.
        let (xmin, xmax, ymin, ymax, zmin, zmax) =
            box_.get().expect("Bnd_Box is void");

        let dx = xmax - xmin;
        let dy = ymax - ymin;
        let dz = zmax - zmin;
        self.gabarit = (dx * dx + dy * dy + dz * dz).sqrt() / 2.0;

        // Cas des toute petite
        self.gabarit += CONFUSION;

        // Initialisation de TheAxe pour les cas singulier
        if !self.is_point {
            let section = self
                .section_curve
                .clone()
                .expect("GeomFill_SectionPlacement: null section");

            let mid = (section.first_parameter() + section.last_parameter()) / 2.0;
            let (p, v) = tangente(section.as_ref(), mid);
            self.the_axe.set_location(p);
            self.the_axe.set_direction(
                // `gp_Dir::SetXYZ` on a null vector is OCCT's
                // `gp_VectorWithNullMagnitude` throw.
                GpDir::new(v.x(), v.y(), v.z()).expect("gp_Dir: null magnitude"),
            );

            // y a t'il un Plan moyen ?
            // `int NbPoles = 0;` (`cxx:182`): the conic branches leave it 0.
            let mut nb_poles: i32 = 0;
            if let Some(circ) = section.gp_circ() {
                self.is_plan = true;
                self.the_axe = circ.axis();
            } else if let Some(elips) = section.gp_ellipse() {
                self.is_plan = true;
                self.the_axe = *elips.axis();
            } else if let Some(hypr) = section.gp_hyperbola() {
                self.is_plan = true;
                self.the_axe = *hypr.axis();
            } else if let Some(parab) = section.gp_parabola() {
                self.is_plan = true;
                self.the_axe = *parab.axis();
            } else if section.is_line() {
                // Pas de Plan !!
                nb_poles = 0;
            } else if let Some(poles) = section.bezier_poles() {
                nb_poles = poles.len() as i32;
            } else if let Some(poles) = section.bspline_poles() {
                nb_poles = poles.len() as i32;
            } else {
                nb_poles = 21;
            }

            if !self.is_plan && nb_poles > 2 {
                let points = mean_plane_points(section.as_ref(), nb_poles);
                // `GeomLib::AxeOfInertia(Pnts, axe, issing, Precision::Confusion())`
                let (axe, is_singular) = geom_lib_axe_of_inertia(&points, CONFUSION);
                if !is_singular {
                    self.is_plan = true;
                    self.the_axe.set_location(axe.location());
                    self.the_axe.set_direction(axe.direction());
                }
            }
        }
    }
}

/// `Extrema_ExtPC::Perform` on the section (`cxx:373-376` + `cxx:426-431`).
pub(super) fn extrema_on_section(
    section: &dyn Curve,
    point: &GpPnt,
) -> (Vec<ExtPcSolution>, bool) {
    extrema_ext_pc_range_done(
        section,
        point,
        section.first_parameter(),
        section.last_parameter(),
    )
}

/// The point cloud `GeomLib::AxeOfInertia` is called on (`cxx:252-364`).
fn mean_plane_points(section: &dyn Curve, nb_poles: i32) -> Vec<GpPnt> {
    let first = section.first_parameter();
    let last = section.last_parameter();

    // `myAdpSection.IsClosed()` (`cxx:302`, `cxx:347`, `cxx:359`): for a
    // B-spline `Geom_BSplineCurve::IsClosed` compares the square distance with
    // `Precision::Computational()` (`Geom_BSplineCurve_1.cxx:146-149`), every
    // other curve uses the endpoint `Precision::Confusion()` form.
    let is_closed = if section.bspline_knots().is_some() {
        section
            .value(first)
            .square_distance(&section.value(last))
            <= COMPUTATIONAL
    } else {
        section.value(first).distance(&section.value(last)) <= CONFUSION
    };

    let mut first = first;
    let mut last = last;
    if section.is_periodic() {
        // Correct boundaries to avoid mistake of LocateU
        let a_curve: Arc<dyn Curve> = if section.is_geom_trimmed() {
            section
                .trim_basis()
                .unwrap_or_else(|| Arc::from(section.clone_dyn()))
        } else {
            Arc::from(section.clone_dyn())
        };
        let u_first = a_curve.first_parameter();
        let a_period = a_curve.period();
        let u1 = u_first + ((first - u_first) / a_period).floor() * a_period;
        let u2 = u1 + a_period;
        if (first - u1).abs() <= PCONFUSION {
            first = u1;
        }
        if (last - u2).abs() <= PCONFUSION {
            last = u2;
        }
    }

    if section.bspline_knots().is_some() {
        let knots = unique_knots(section.bspline_knots().unwrap());
        let (i1, i2) = bspline_locate_u(section, &knots, first, CONFUSION);
        let (i3, i4) = bspline_locate_u(section, &knots, last, CONFUSION);
        let nb_knots = i3 - i2 + 1;

        let nb_local_pnts = 10i32;
        let mut nb_pnts = (nb_knots - 1) * nb_local_pnts;
        if nb_pnts < 0 {
            nb_pnts = 0;
        }
        if i1 != i2 {
            nb_pnts += nb_local_pnts;
        }
        if i3 != i4 && first < knot_at(&knots, i3) {
            nb_pnts += nb_local_pnts;
        }
        if !is_closed {
            nb_pnts += 1;
        }

        let mut pnts = Vec::with_capacity(nb_pnts.max(0) as usize);
        if i1 != i2 {
            let locallast = if knot_at(&knots, i2) < last {
                knot_at(&knots, i2)
            } else {
                last
            };
            let delta = (locallast - first) / nb_local_pnts as f64;
            for j in 0..nb_local_pnts {
                let t = first + j as f64 * delta;
                pnts.push(section.value(t));
            }
        }
        for i in i2..i3 {
            let mut t = knot_at(&knots, i);
            let delta = (knot_at(&knots, i + 1) - t) / nb_local_pnts as f64;
            for _ in 0..nb_local_pnts {
                pnts.push(section.value(t));
                t += delta;
            }
        }
        if i3 != i4 && first < knot_at(&knots, i3) {
            let mut t = knot_at(&knots, i3);
            let delta = (last - t) / nb_local_pnts as f64;
            for _ in 0..nb_local_pnts {
                pnts.push(section.value(t));
                t += delta;
            }
        }
        if !is_closed {
            pnts.push(section.value(last));
        }
        return pnts;
    }

    // other type
    let mut nb_pnts = nb_poles - 1;
    if !is_closed {
        nb_pnts += 1;
    }
    let mut pnts = Vec::with_capacity(nb_pnts.max(0) as usize);
    let delta = (last - first) / (nb_poles - 1) as f64;
    for i in 0..(nb_poles - 1) {
        let t = first + i as f64 * delta;
        pnts.push(section.value(t));
    }
    if !is_closed {
        pnts.push(section.value(last));
    }
    pnts
}

/// `knots(I)` with OCCT's 1-based `Array1` indexing.
fn knot_at(knots: &[f64], i: i32) -> f64 {
    knots[(i - 1) as usize]
}

/// Distinct values of a flat knot vector, i.e. `Geom_BSplineCurve::Knots()`
/// (the non-flat array `LocateU` indexes by default).
fn unique_knots(flat: &[f64]) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    for &k in flat {
        if out.last().map_or(true, |&prev| prev != k) {
            out.push(k);
        }
    }
    out
}

/// `Geom_BSplineCurve::LocateU(U, ParametricTolerance, I1, I2)`
/// (`Geom_BSplineCurve_1.cxx:702-752`) on the curve's distinct knots
/// (`WithKnotRepetition = false`).
fn bspline_locate_u(
    section: &dyn Curve,
    knots: &[f64],
    u: f64,
    parametric_tolerance: f64,
) -> (i32, i32) {
    let len = knots.len() as i32;
    let mut new_u = u;

    // `PeriodicNormalization` (`Geom_BSplineCurve_1.cxx:711`): the flat knot
    // window is `myFlatKnots(1)` .. `myFlatKnots(Upper - Degree)`.
    if section.is_periodic() {
        let flat = section
            .bspline_knots()
            .expect("GeomFill_SectionPlacement: B-spline knots");
        let degree = section.nurbs_degree().unwrap_or(0);
        let p_lo = flat[0];
        let p_hi = flat[flat.len() - 1 - degree];
        let period = section.period();
        while new_u > p_hi {
            new_u -= period;
        }
        while new_u < p_lo {
            new_u += period;
        }
    }

    let u_first = knots[0];
    let u_last = knots[(len - 1) as usize];
    let tol = parametric_tolerance.abs();

    if (new_u - u_first).abs() <= tol {
        (1, 1)
    } else if (new_u - u_last).abs() <= tol {
        (len, len)
    } else if new_u < u_first {
        (0, 1)
    } else if new_u > u_last {
        (len, len + 1)
    } else {
        let mut i1 = hunt_occt(knots, new_u).clamp(1, len);
        while i1 + 1 <= len && (knot_at(knots, i1 + 1) - new_u).abs() <= tol {
            i1 += 1;
        }
        if (knot_at(knots, i1) - new_u).abs() <= tol {
            (i1, i1)
        } else {
            (i1, i1 + 1)
        }
    }
}
