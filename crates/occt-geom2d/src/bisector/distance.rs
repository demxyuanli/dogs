//! `Bisector_BisecAna::Distance` (`Bisector_BisecAna.cxx:64-221`): the distance
//! from a point to a candidate bisector plus the sector test that decides
//! whether the candidate is oriented the right way.

use occt_core::elib::clib2d;
use occt_core::gp::{GpDir2d, GpPnt2d, GpVec2d};

use crate::gcc_int::GccIntBisec;

use super::BisectorBisecAna;

/// Out-parameters of [`BisectorBisecAna::distance`]: OCCT passes `aparameter`,
/// `asense` and `astatus` by reference from `Bisector_BisecAna.cxx:66-67`.
#[derive(Debug, Clone, Copy)]
pub struct DistanceResult {
    /// The returned distance from `apoint` to the evaluated point of `abisector`.
    pub distance: f64,
    /// `aparameter`: the parameter of `apoint` on the bisector arc.
    pub parameter: f64,
    /// `asense`: the orientation of the bisector arc.
    pub sense: bool,
    /// `astatus`: false when the candidate is not in the proper sector.
    pub status: bool,
}

impl BisectorBisecAna {
    /// `Bisector_BisecAna::Distance` (`cxx:64-221`).
    ///
    /// The arc `abisector` is evaluated at the parameter of `apoint`; the signed
    /// areas of `(V1, tangent)` and `(tangent, V2)` then decide whether the
    /// candidate lies in the sector delimited by `afirstvector`/`asecondvector`.
    /// `is_bisec_of_two_lines` selects the two-line variant, which skips the
    /// whole second block (`cxx:163-215`).
    pub fn distance(
        apoint: &GpPnt2d,
        abisector: &GccIntBisec,
        afirstvector: &GpVec2d,
        asecondvector: &GpVec2d,
        vec_ref: &GpVec2d,
        adirection: f64,
        mut asense: bool,
        is_bisec_of_two_lines: bool,
    ) -> DistanceResult {
        let mut status = true;
        // `gp_Pnt2d point; gp_Vec2d tangent;` (`cxx:88-89`) - both keep their
        // default value when `ArcType()` has no branch, which is what OCCT does
        // for `GccInt_Pnt`.
        let mut point = GpPnt2d::zero();
        let mut tangent = GpVec2d::zero();
        let mut aparameter = 0.0;

        match abisector {
            GccIntBisec::Line(l) => {
                aparameter = clib2d::parameter_lin2d(l, apoint);
                let (p, t) = clib2d::line_d1_ax2d(aparameter, l.position());
                point = p;
                tangent = t;
            }
            GccIntBisec::Circle(c) => {
                aparameter = clib2d::parameter_circ2d(c, apoint);
                let (p, t) = clib2d::circle_d1_ax22d(aparameter, c.position(), c.radius());
                point = p;
                tangent = t;
            }
            GccIntBisec::Hyperbola(h) => {
                aparameter = clib2d::parameter_hypr2d(h, apoint);
                let (p, t) = clib2d::hyperbola_d1_ax22d(
                    aparameter,
                    &h.pos,
                    h.major_radius,
                    h.minor_radius,
                );
                point = p;
                tangent = t;
            }
            GccIntBisec::Parabola(pb) => {
                aparameter = clib2d::parameter_parab2d(pb, apoint);
                let (p, t) = clib2d::parabola_d1_ax22d(aparameter, &pb.pos, pb.focal);
                point = p;
                tangent = t;
            }
            GccIntBisec::Ellipse(e) => {
                aparameter = clib2d::parameter_elips2d(e, apoint);
                let (p, t) =
                    clib2d::ellipse_d1_ax22d(aparameter, &e.pos, e.major_radius, e.minor_radius);
                point = p;
                tangent = t;
            }
            GccIntBisec::Point(_) => {}
        }

        let distance = apoint.distance(&point);

        let afirstdir = GpDir2d::from_vec2d(afirstvector).unwrap_or_default();
        let aseconddir = GpDir2d::from_vec2d(asecondvector).unwrap_or_default();
        let tangdir = GpDir2d::from_vec2d(&tangent).unwrap_or_default();
        let secdirrev = aseconddir.reversed();

        if asense {
            // `cxx:136-161`: the status is only determined when the tangent to
            // the bisector is the bisector of the two vectors.
            let sin_plat = 1.0e-3;
            if afirstdir.crossed(&aseconddir).abs() < sin_plat {
                if afirstdir.dot(&aseconddir) >= 0.0 {
                    // Tangent mixed: the scalar product must be close to 1.
                    if tangdir.dot(&afirstdir).abs() < 0.5 {
                        status = false;
                    }
                } else {
                    // Opposed tangents: the scalar product must be close to 0.
                    if tangdir.dot(&afirstdir).abs() > 0.5 {
                        status = false;
                    }
                }
            } else if (afirstdir.crossed(&tangdir)) * (tangdir.crossed(&aseconddir)) < -1.0e-8 {
                status = false;
            }
        } else {
            asense = true;

            if !is_bisec_of_two_lines {
                // `cxx:163-215`, including the `1.e-4` tolerance introduced by
                // the 2002 Airbus fix (`cxx:166-172`).
                let a_tol = 1.0e-4;
                let side = (afirstdir.crossed(&secdirrev)) * adirection;

                if side < -0.1 {
                    // Input.
                    if (afirstdir.crossed(&tangdir)) * adirection < a_tol
                        && (secdirrev.crossed(&tangdir)) * adirection < a_tol
                    {
                        asense = false;
                    }
                } else if side > 0.1 {
                    // Output.
                    if (afirstdir.crossed(&tangdir)) * adirection < a_tol
                        || (secdirrev.crossed(&tangdir)) * adirection < a_tol
                    {
                        asense = false;
                    }
                } else if afirstdir.dot(&secdirrev) > 0.0 {
                    // Flat, tangents together.
                    if (afirstdir.crossed(&tangdir)) * adirection < 0.0 {
                        asense = false;
                    }
                } else {
                    // Flat, tangents turned back (`cxx:191-203`).
                    if afirstdir.dot(&tangdir) < 0.0 {
                        asense = false;
                    }
                }

                // `cxx:206-213`: the OCC26185 fix - the reference vector pins
                // the sense when it is not null.
                if vec_ref.square_magnitude() != 0.0 {
                    let dir_ref = GpDir2d::from_vec2d(vec_ref).unwrap_or_default();
                    if tangdir.dot(&dir_ref) < 0.0 {
                        asense = false;
                    }
                }
            }
        }

        DistanceResult {
            distance,
            parameter: aparameter,
            sense: asense,
            status,
        }
    }

    /// `Bisector_BisecAna::Parameter` (`cxx:1659-1692`): the parameter of `P` on
    /// the **basis** curve, dispatched on the exact elementary type.
    ///
    /// OCCT switches on `STANDARD_TYPE(Geom2d_Line/Circle/Hyperbola/Parabola/
    /// Ellipse)`; the `gp_*2d` accessors used here also answer for a
    /// `Geom2d_OffsetCurve` of that basis, which cannot arise because the basis
    /// of a `Bisector_BisecAna` is always one of the `GccAna_*Bisec` arcs.
    pub fn parameter(&self, p: &GpPnt2d) -> f64 {
        let Some(basis) = self.geom2d_curve() else {
            return 0.0;
        };
        if let Some(l) = basis.gp_lin2d() {
            return clib2d::parameter_lin2d(&l, p);
        }
        if let Some(c) = basis.gp_circ2d() {
            return clib2d::parameter_circ2d(&c, p);
        }
        if let Some(h) = basis.gp_hypr2d() {
            return clib2d::parameter_hypr2d(&h, p);
        }
        if let Some(pb) = basis.gp_parab2d() {
            return clib2d::parameter_parab2d(&pb, p);
        }
        if let Some(e) = basis.gp_elips2d() {
            return clib2d::parameter_elips2d(&e, p);
        }
        0.0
    }
}
