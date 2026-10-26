//! `AppBlend_AppSurf::Perform` — the three overloads.
//!
//! Source: `TKGeomAlgo/AppBlend/AppBlend_AppSurf.gxx:170-1048` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `InternalPerform(Lin, F, SpApprox, UseSmoothing)` (`gxx:189-592`): the
//!   `!SpApprox` arm runs `AppDef_Compute`
//!   ([`ComputeLine`](crate::approx::compute_line::ComputeLine)); the `SpApprox`
//!   arm runs `AppDef_BSplineCompute`
//!   ([`BsplComputeLine`](crate::approx::bspl_compute_line::BsplComputeLine)).
//!   `UseSmoothing` (`gxx:472-541`) drives `AppDef_Variational`
//!   ([`Variational`](crate::app_def::variational::Variational)).
//! - `Perform(Lin, SecGen, NbMaxP)` (`gxx:597-1048`) is the troncon split: it
//!   scales every row into `[0, 1]`, approximates each troncon with
//!   `AppDef_Compute`, and concatenates the pieces with
//!   `Approx_MCurvesToBSpCurve` ([`MCurvesToBSpCurve`]).
//! - `std::cout` diagnostics (`gxx:412-415`, `:466-468`, `:532-537`, `:973-977`,
//!   `:998-1001`) are `OCCT_DEBUG`-only and not ported; the values they print
//!   are still computed.
//! - `static bool scal = 1` (`gxx:36`) is never written in the file, so the
//!   `if (!scal)` reset block (`gxx:759-773`) is dead code kept for parity.

use occt_core::approx::ApproxParametrizationType;
use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::REAL_SMALL;
use occt_math::MathVector;

use crate::app_def::multi_line::MultiLine;
use crate::app_def::multi_point_constraint::MultiPointConstraint;
use crate::app_def::variational::Variational;
use crate::app_par_curves::constraint::{Constraint, ConstraintCouple};
use crate::app_par_curves::multi_bspline_curve::MultiBSpCurve;
use crate::approx::bspl_compute_line::BsplComputeLine;
use crate::approx::compute_line::ComputeLine;
use crate::approx::m_curves_to_bsp_curve::MCurvesToBSpCurve;
use crate::geom_fill::app_surf::debug::app_blend_get_context_approx_with_no_tgt;
use crate::geom_fill::app_surf::row_scratch::{RowScratch, SCAL};
use crate::geom_fill::app_surf::{AppBlendSectionGenerator, AppSurf};
use crate::geom_fill::line::Line;

impl AppSurf {
    /// `Perform(Lin, SecGen, SpApprox)` (`gxx:170-176`).
    ///
    /// `SpApprox == true` selects the `AppDef_BSplineCompute` branch
    /// (`gxx:420-471`, [`BsplComputeLine`](crate::approx::bspl_compute_line::BsplComputeLine)).
    pub fn perform<G: AppBlendSectionGenerator>(
        &mut self,
        lin: &Line,
        f: &G,
        sp_approx: bool,
    ) -> Result<(), String> {
        self.internal_perform(lin, f, sp_approx, false)
    }

    /// `PerformSmoothing(Lin, SecGen)` (`gxx:180-183`): the variational branch
    /// (`gxx:472-541`).
    pub fn perform_smoothing<G: AppBlendSectionGenerator>(
        &mut self,
        lin: &Line,
        f: &G,
    ) -> Result<(), String> {
        self.internal_perform(lin, f, true, true)
    }

    /// `Perform(Lin, SecGen, NbMaxP)` (`gxx:597-1048`).
    pub fn perform_nb_max_p<G: AppBlendSectionGenerator>(
        &mut self,
        lin: &Line,
        f: &G,
        nb_max_p: i32,
    ) -> Result<(), String> {
        self.done = false;
        // `Lin.IsNull()` (`gxx:602-605`): `Line` is a value type here.
        let nb_point_tot = lin.nb_points();

        let (nb_upoles, _nb_uknots, udeg, nb_poles2d) = f.get_shape()?;
        self.udeg = udeg;
        self.tab_u_knots = f.knots()?;
        self.tab_u_mults = f.mults()?;
        let nb_upoles = nb_upoles.max(0) as usize;
        let nb_poles2d = nb_poles2d.max(0) as usize;
        let mut s = RowScratch::new(nb_upoles, nb_poles2d);

        self.seq_poles2d.clear();

        let mut cfirst = Constraint::NoConstraint;
        let mut clast = Constraint::NoConstraint;
        let mut mytol3d = 0.0f64;
        let mut mytol2d = 0.0f64;

        // `X = Y = Z = RealLast(); DX = DY = DZ = RealFirst();` (`gxx:628-630`).
        let mut x = f64::MAX;
        let mut y = f64::MAX;
        let mut z = f64::MAX;
        let mut dx = -f64::MAX;
        let mut dy = -f64::MAX;
        let mut dz = -f64::MAX;

        // `X2d` / `Y2d` start at `RealLast`, `DX2d` / `DY2d` at `RealFirst`
        // (`gxx:635-641`).
        let n2 = nb_poles2d.max(1);
        let mut x2d = vec![f64::MAX; n2];
        let mut y2d = vec![f64::MAX; n2];
        let mut dx2d = vec![-f64::MAX; n2];
        let mut dy2d = vec![-f64::MAX; n2];

        // On calcule les boites de chaque ligne (`gxx:648-702`).
        for i in 1..=nb_point_tot {
            let (values, _with_deriv) = f.section_with_derivatives(lin.point(i))?;
            if values.poles.len() != nb_upoles {
                return Err("Standard_OutOfRange: AppBlend_AppSurf::Section".into());
            }
            for j in 1..=nb_upoles {
                let p = values.poles[j - 1];
                let (px, py, pz) = (p.x(), p.y(), p.z());
                if px < x {
                    x = px;
                }
                if px > dx {
                    dx = px;
                }
                if py < y {
                    y = py;
                }
                if py > dy {
                    dy = py;
                }
                if pz < z {
                    z = pz;
                }
                if pz > dz {
                    dz = pz;
                }
            }
            for j in 1..=nb_poles2d {
                let p = values.poles2d[j - 1];
                let (px, py) = (p.x(), p.y());
                if px < x2d[j - 1] {
                    x2d[j - 1] = px;
                }
                if px > dx2d[j - 1] {
                    dx2d[j - 1] = px;
                }
                if py < y2d[j - 1] {
                    y2d[j - 1] = py;
                }
                if py > dy2d[j - 1] {
                    dy2d[j - 1] = py;
                }
            }
        }

        // On calcule pour chaque ligne la transformation vers 0 1.
        let seuil = 1000.0 * self.tol3d;
        let seuil2d = 1000.0 * self.tol2d;
        if (dx - x) < seuil {
            dx = 1.0;
            x = 0.0;
        } else {
            dx = 1.0 / (dx - x);
            x *= -dx;
        }
        if (dy - y) < seuil {
            dy = 1.0;
            y = 0.0;
        } else {
            dy = 1.0 / (dy - y);
            y *= -dy;
        }
        if (dz - z) < seuil {
            dz = 1.0;
            z = 0.0;
        } else {
            dz = 1.0 / (dz - z);
            z *= -dz;
        }
        for j in 1..=nb_poles2d {
            if (dx2d[j - 1] - x2d[j - 1]) < seuil2d {
                dx2d[j - 1] = 1.0;
                x2d[j - 1] = 0.0;
            } else {
                dx2d[j - 1] = 1.0 / (dx2d[j - 1] - x2d[j - 1]);
                x2d[j - 1] *= -dx2d[j - 1];
            }
            if (dy2d[j - 1] - y2d[j - 1]) < seuil2d {
                dy2d[j - 1] = 1.0;
                y2d[j - 1] = 0.0;
            } else {
                dy2d[j - 1] = 1.0 / (dy2d[j - 1] - y2d[j - 1]);
                y2d[j - 1] *= -dy2d[j - 1];
            }
        }
        if !SCAL {
            dx = 1.0;
            x = 0.0;
            dy = 1.0;
            y = 0.0;
            dz = 1.0;
            z = 0.0;
            for j in 1..=nb_poles2d {
                dx2d[j - 1] = 1.0;
                x2d[j - 1] = 0.0;
                dy2d[j - 1] = 1.0;
                y2d[j - 1] = 0.0;
            }
        }

        // Keep "inter-troncons" parameters, not only first and last
        // (`gxx:774-784`).
        let mut a_param_seq: Vec<f64> = Vec::new();
        if self.knownp {
            a_param_seq.push(f.parameter(lin.point(1))?);
        }

        let mut concat = MCurvesToBSpCurve::new();

        // On calcule le nombre de troncons (`gxx:788-808`).
        let mut nbtronc = nb_point_tot / nb_max_p;
        let mut reste = nb_point_tot - (nbtronc * nb_max_p);
        let mut nmax = nb_max_p;
        if nbtronc > 0 && reste > 0 {
            nmax = nb_point_tot / (nbtronc + 1);
            if nmax > (2 * nb_max_p) / 3 {
                nbtronc += 1;
                reste = nb_point_tot - (nbtronc * nmax);
            } else {
                nmax = nb_max_p;
            }
        } else if nbtronc == 0 {
            nbtronc = 1;
            nmax = reste;
            reste = 0;
        }

        let mut troncsize = vec![0i32; nbtronc.max(0) as usize];
        let mut troncstart = vec![0i32; nbtronc.max(0) as usize];
        let rab = reste / nbtronc + 1;
        let mut start = 1;
        for itronc in 1..=nbtronc {
            troncstart[(itronc - 1) as usize] = start;
            let rabrab = rab.min(reste);
            if reste > 0 {
                reste -= rabrab;
            }
            troncsize[(itronc - 1) as usize] = nmax + rabrab + 1;
            start += nmax + rabrab;
        }
        troncsize[(nbtronc - 1) as usize] -= 1;

        for itronc in 1..=nbtronc {
            let nb_point = troncsize[(itronc - 1) as usize];
            let st_point = troncstart[(itronc - 1) as usize];
            let mut mult_l = MultiLine::with_size(nb_point.max(0) as usize);

            for i in 1..=nb_point {
                let i_lin = st_point + i - 1;
                let (values, mut with_deriv) = f.section_with_derivatives(lin.point(i_lin))?;
                if app_blend_get_context_approx_with_no_tgt() {
                    with_deriv = false;
                }
                s.load(&values, with_deriv)?;

                for j in 1..=nb_poles2d {
                    let p = s.tab_p2d[j - 1];
                    let (px, py) = (p.x(), p.y());
                    s.tab_app_p2d[j - 1] = GpPnt2d::new(
                        dx2d[j - 1] * px + x2d[j - 1],
                        dy2d[j - 1] * py + y2d[j - 1],
                    );
                    if with_deriv {
                        let v = s.tab_v2d[j - 1];
                        s.tab_app_v2d[j - 1] =
                            GpVec2d::new(dx2d[j - 1] * v.x(), dy2d[j - 1] * v.y());
                    }
                }
                for j in 1..=nb_upoles {
                    let w = s.tab_w[j - 1];
                    if with_deriv {
                        let dw = s.tab_dw[j - 1];
                        s.tab_app_v2d[nb_poles2d + j - 1] = GpVec2d::new(dw, 0.0);
                        let mut new_dv = GpPnt::new(0.0, 0.0, 0.0).xyz().clone();
                        new_dv.set_linear_form(
                            dw,
                            s.tab_app_p[j - 1].xyz(),
                            w,
                            &s.tab_app_v[j - 1].coord,
                        );
                        s.tab_app_v[j - 1] = GpVec::new(dx * new_dv.x(), dy * new_dv.y(), dz * new_dv.z());
                    }
                    s.tab_app_p[j - 1] = GpPnt::from_xyz(&s.tab_app_p[j - 1].xyz().multiply(w));
                    s.tab_app_p2d[nb_poles2d + j - 1] = GpPnt2d::new(w, 0.0);
                    let p = s.tab_app_p[j - 1];
                    s.tab_app_p[j - 1] = GpPnt::new(
                        dx * p.x() + x,
                        dy * p.y() + y,
                        dz * p.z() + z,
                    );
                }

                let (mult_p, cons) = s.multi_point(with_deriv);
                if i == 1 {
                    cfirst = cons;
                } else if i == nb_point {
                    clast = cons;
                }
                mult_l.set_value(i, &mult_p);
            }

            // IFV 04.06.07 occ13904
            if nb_point == 2 {
                self.dmin = 1;
                if cfirst == Constraint::PassPoint && clast == Constraint::PassPoint {
                    self.dmax = 1;
                }
            }

            let mut ufloc = 0.0f64;
            let mut ulloc = 0.0f64;
            let mut the_approx = ComputeLine::with_settings(
                self.dmin,
                self.dmax,
                self.tol3d,
                self.tol2d,
                self.nbit,
                true,
                ApproxParametrizationType::ChordLength,
                false,
            );
            if self.knownp {
                let mut the_params = MathVector::new(1, nb_point as usize);
                // On recale les parametres entre 0 et 1.
                ufloc = f.parameter(lin.point(st_point))?;
                ulloc = f.parameter(lin.point(st_point + nb_point - 1))?;
                for i in 1..=nb_point {
                    let i_lin = st_point + i - 1;
                    the_params.set_value(
                        i as usize,
                        (f.parameter(lin.point(i_lin))? - ufloc) / (ulloc - ufloc),
                    );
                }
                the_approx = ComputeLine::from_parameters(
                    &the_params,
                    self.dmin,
                    self.dmax,
                    self.tol3d,
                    self.tol2d,
                    self.nbit,
                    true,
                    true,
                );
            }
            the_approx.set_constraints(cfirst, clast);
            the_approx.perform(&mult_l);

            // To know internal parameters if multicurve is approximated by
            // several Bezier's (`gxx:923-925`).
            let mut a_pole_dist_seq: Vec<f64> = Vec::new();
            let mut a_whole_dist = 0.0f64;

            for index in 1..=the_approx.nb_multi_curves() as i32 {
                let (the_tol3d, the_tol2d) = the_approx.error(index);
                mytol3d = mytol3d.max(the_tol3d / dx);
                mytol3d = mytol3d.max(the_tol3d / dy);
                mytol3d = mytol3d.max(the_tol3d / dz);
                {
                    let mucu = the_approx.change_value(index);
                    for j in 1..=nb_upoles {
                        mucu.transform(
                            j as i32,
                            -x / dx,
                            1.0 / dx,
                            -y / dy,
                            1.0 / dy,
                            -z / dz,
                            1.0 / dz,
                        );
                    }
                    for j in 1..=nb_poles2d {
                        mucu.transform2d(
                            (j + nb_upoles) as i32,
                            -x2d[j - 1] / dx2d[j - 1],
                            1.0 / dx2d[j - 1],
                            -y2d[j - 1] / dy2d[j - 1],
                            1.0 / dy2d[j - 1],
                        );
                        mytol2d = mytol2d.max(the_tol2d / dx2d[j - 1]);
                        mytol2d = mytol2d.max(the_tol2d / dy2d[j - 1]);
                    }
                }
                let mucu = the_approx.value(index).clone();
                concat.append(&mucu);

                if self.knownp && the_approx.nb_multi_curves() > 1 {
                    // `mucu.Pole(Index, 1)` / `mucu.Pole(Index,
                    // mucu.NbPoles())` (`gxx:953-954`): OCCT passes the
                    // multi-curve index as `CuIndex`; kept verbatim.
                    let a_first_pole = mucu.pole(index, 1).clone();
                    let a_last_pole = mucu.pole(index, mucu.nb_poles() as i32).clone();
                    a_pole_dist_seq.push(a_first_pole.distance(&a_last_pole));
                    a_whole_dist += a_pole_dist_seq[a_pole_dist_seq.len() - 1];
                }
            }

            if self.knownp {
                let mut i_u = ufloc;
                for i_dist in 0..a_pole_dist_seq.len().saturating_sub(1) {
                    i_u += a_pole_dist_seq[i_dist] / a_whole_dist * (ulloc - ufloc);
                    a_param_seq.push(i_u);
                }
                a_param_seq.push(ulloc);
            }
        }

        self.tol3d_reached = mytol3d;
        self.tol2d_reached = mytol2d;

        concat.perform();
        let mult_c = concat.value().clone();
        self.store_knots_and_mults(&mult_c);

        if self.knownp {
            if a_param_seq.len() != self.tab_v_knots.len() {
                let u1 = f.parameter(lin.point(1))?;
                let u2 = f.parameter(lin.point(lin.nb_points()))?;
                occt_core::bspl::knots::reparametrize(u1, u2, &mut self.tab_v_knots);
            } else {
                for (i_knot, value) in a_param_seq.iter().enumerate() {
                    self.tab_v_knots[i_knot] = *value;
                }
            }
        }

        if !self.unpack_poles(&mult_c, nb_upoles, nb_poles2d)? {
            return Ok(());
        }
        self.done = true;
        Ok(())
    }

    /// `InternalPerform(Lin, SecGen, SpApprox, UseSmoothing)` (`gxx:189-592`).
    fn internal_perform<G: AppBlendSectionGenerator>(
        &mut self,
        lin: &Line,
        f: &G,
        sp_approx: bool,
        use_smoothing: bool,
    ) -> Result<(), String> {
        self.done = false;
        // `Lin.IsNull()` (`gxx:193-196`): `Line` is a value type here.

        let nb_point = lin.nb_points();
        let mut mult_l = MultiLine::with_size(nb_point.max(0) as usize);

        let (nb_upoles, _nb_uknots, udeg, nb_poles2d) = f.get_shape()?;
        self.udeg = udeg;
        self.tab_u_knots = f.knots()?;
        self.tab_u_mults = f.mults()?;
        let nb_upoles = nb_upoles.max(0) as usize;
        let nb_poles2d = nb_poles2d.max(0) as usize;
        let mut s = RowScratch::new(nb_upoles, nb_poles2d);

        self.seq_poles2d.clear();

        let (mult_p, cfirst) = self.row_constraint(f, lin.point(1), &mut s)?;
        mult_l.set_value(1, &mult_p);

        for i in 2..nb_point {
            let (mult_p, _cons) = if sp_approx {
                // `gxx:277-292`: the `SpApprox` overload writes no derivatives.
                self.row_constraint_without_derivatives(f, lin.point(i), &mut s)?
            } else {
                self.row_constraint(f, lin.point(i), &mut s)?
            };
            mult_l.set_value(i, &mult_p);
        }

        let (mult_p, clast) = self.row_constraint(f, lin.point(nb_point), &mut s)?;
        mult_l.set_value(nb_point, &mult_p);

        // IFV 04.06.07 occ13904
        if nb_point == 2 {
            self.dmin = 1;
            if cfirst == Constraint::PassPoint && clast == Constraint::PassPoint {
                self.dmax = 1;
            }
        }

        let mult_c = if !sp_approx {
            let mut theapprox = if self.knownp {
                let mut the_params = MathVector::new(1, nb_point as usize);
                // On recale les parametres entre 0 et 1.
                the_params.set_value(1, 0.0);
                the_params.set_value(nb_point as usize, 1.0);
                let uf = f.parameter(lin.point(1))?;
                let ul = f.parameter(lin.point(nb_point))? - uf;
                for i in 2..nb_point {
                    the_params.set_value(i as usize, (f.parameter(lin.point(i))? - uf) / ul);
                }
                ComputeLine::from_parameters(
                    &the_params,
                    self.dmin,
                    self.dmax,
                    self.tol3d,
                    self.tol2d,
                    self.nbit,
                    true,
                    true,
                )
            } else {
                ComputeLine::with_settings(
                    self.dmin,
                    self.dmax,
                    self.tol3d,
                    self.tol2d,
                    self.nbit,
                    true,
                    self.paramtype,
                    false,
                )
            };
            theapprox.set_constraints(cfirst, clast);
            theapprox.perform(&mult_l);

            let mut mytol3d = 0.0f64;
            let mut mytol2d = 0.0f64;
            for index in 1..=theapprox.nb_multi_curves() as i32 {
                let (the_tol3d, the_tol2d) = theapprox.error(index);
                mytol3d = mytol3d.max(the_tol3d);
                mytol2d = mytol2d.max(the_tol2d);
            }
            // `gxx:401-418` leaves `tol3dreached` untouched on this arm; the
            // two maxima only feed the `OCCT_DEBUG` print.
            let _ = (mytol3d, mytol2d);
            theapprox.spline_value().clone()
        } else if !use_smoothing {
            // `gxx:421-471`: `AppDef_BSplineCompute`.
            let use_squares = self.nbit == 0;
            let mut theapprox = BsplComputeLine::with_settings(
                self.dmin,
                self.dmax,
                self.tol3d,
                self.tol2d,
                self.nbit,
                true,
                self.paramtype,
                use_squares,
            );
            // `gxx:428-443`, kept verbatim: a `C0` request sets 0 and is then
            // overwritten by the trailing `else`, which is what OCCT does.
            if self.continuity == Shape::C0 {
                theapprox.set_continuity(0);
            }
            if self.continuity == Shape::C1 {
                theapprox.set_continuity(1);
            } else if self.continuity == Shape::C2 {
                theapprox.set_continuity(2);
            } else {
                theapprox.set_continuity(3);
            }

            theapprox.set_constraints(cfirst, clast);

            if self.knownp {
                let mut the_params = MathVector::new(1, nb_point as usize);
                // On recale les parametres entre 0 et 1.
                the_params.set_value(1, 0.0);
                the_params.set_value(nb_point as usize, 1.0);
                let uf = f.parameter(lin.point(1))?;
                let ul = f.parameter(lin.point(nb_point))? - uf;
                for i in 2..nb_point {
                    the_params.set_value(i as usize, (f.parameter(lin.point(i))? - uf) / ul);
                }
                theapprox.init(
                    self.dmin,
                    self.dmax,
                    self.tol3d,
                    self.tol2d,
                    self.nbit,
                    true,
                    ApproxParametrizationType::IsoParametric,
                    true,
                );
                theapprox.set_parameters(&the_params);
            }
            theapprox.perform(&mult_l);

            let (mytol3d, mytol2d) = theapprox.error();
            self.tol3d_reached = mytol3d;
            self.tol2d_reached = mytol2d;
            theapprox.value().clone()
        } else {
            // `gxx:472-541`: the variational algorithm.
            //
            // Every `return;` of this arm leaves `done` at `false`, which the
            // port spells `return Ok(())`.
            let mut tab_of_cc = vec![ConstraintCouple::default(); nb_point.max(1) as usize]; // `gxx:476-485`
            for i in 1..=nb_point {
                tab_of_cc[(i - 1) as usize] =
                    ConstraintCouple::from_parts(i, Constraint::NoConstraint);
            }
            tab_of_cc[0].set_constraint(cfirst);
            tab_of_cc[(nb_point - 1) as usize].set_constraint(clast);

            // `AppDef_Variational Variation(multL, 1, NbPoint, TABofCC)`
            // (`gxx:489`), whose defaults come from
            // `AppDef_Variational.hxx:62-74`.
            let mut variation = match Variational::with_defaults(&mult_l, 1, nb_point, &tab_of_cc) {
                Ok(v) => v,
                Err(()) => return Ok(()),
            };

            // `gxx:492-505`.
            let the_max_segments = 1000;
            let the_with_min_max = false;
            let the_with_cutting = true;

            if variation.set_max_degree(self.dmax).is_err() {
                return Ok(());
            }
            if variation.set_continuity(self.continuity).is_err() {
                return Ok(());
            }
            let _ = variation.set_max_segment(the_max_segments);
            if variation.set_tolerance(self.tol3d).is_err() {
                return Ok(());
            }
            if variation.set_with_min_max(the_with_min_max).is_err() {
                return Ok(());
            }
            let _ = variation.set_with_cutting(the_with_cutting);
            variation.set_nb_iterations(self.nbit);
            if variation
                .set_criterium_weight3(self.critweights[0], self.critweights[1], self.critweights[2])
                .is_err()
            {
                return Ok(());
            }

            if !variation.is_created() {
                return Ok(()); // `gxx:508-511`.
            }
            if variation.is_over_constrained() {
                return Ok(()); // `gxx:513-516`.
            }
            if variation.approximate().is_err() {
                return Ok(()); // `gxx:518-525`: `catch (Standard_Failure const&)`.
            }
            if !variation.is_done() {
                return Ok(()); // `gxx:527-530`.
            }

            let mytol3d = match variation.max_error() {
                Ok(e) => e,
                Err(()) => return Ok(()),
            };
            let mytol2d = 0.0;
            self.tol3d_reached = mytol3d;
            self.tol2d_reached = mytol2d;
            match variation.value() {
                Ok(v) => v,
                Err(()) => return Ok(()),
            }
        };

        self.store_knots_and_mults(&mult_c);

        if self.knownp && !use_smoothing {
            // `BSplCLib::Reparametrize(F.Parameter(Lin->Point(1)),
            // F.Parameter(Lin->Point(NbPoint)), tabVKnots)` (`gxx:552-557`).
            let u1 = f.parameter(lin.point(1))?;
            let u2 = f.parameter(lin.point(nb_point))?;
            occt_core::bspl::knots::reparametrize(u1, u2, &mut self.tab_v_knots);
        }

        if !self.unpack_poles(&mult_c, nb_upoles, nb_poles2d)? {
            return Ok(());
        }
        self.done = true;
        Ok(())
    }

    /// `F.Section(P, Poles, Poles2d, Weigths)` plus the weight scaling and the
    /// `AppDef_MultiPointConstraint` it feeds (`gxx:279-291`): the `SpApprox`
    /// interior-row overload, which writes no derivatives.
    fn row_constraint_without_derivatives<G: AppBlendSectionGenerator>(
        &self,
        f: &G,
        p: i32,
        s: &mut RowScratch,
    ) -> Result<(MultiPointConstraint, Constraint), String> {
        let values = f.section(p)?;
        s.load(&values, false)?;
        s.weight_poles(false);
        Ok(s.multi_point(false))
    }

    /// One `F.Section(...)` plus the weight scaling and the
    /// `AppDef_MultiPointConstraint` it feeds (`gxx:236-271`).
    fn row_constraint<G: AppBlendSectionGenerator>(
        &self,
        f: &G,
        p: i32,
        s: &mut RowScratch,
    ) -> Result<(MultiPointConstraint, Constraint), String> {
        let (values, mut with_deriv) = f.section_with_derivatives(p)?;
        if app_blend_get_context_approx_with_no_tgt() {
            with_deriv = false;
        }
        s.load(&values, with_deriv)?;
        s.weight_poles(with_deriv);
        Ok(s.multi_point(with_deriv))
    }

    /// `vdeg = multC.Degree()` and the two knot / mult arrays
    /// (`gxx:543-561`, `:982-1017`).
    fn store_knots_and_mults(&mut self, mult_c: &MultiBSpCurve) {
        self.vdeg = mult_c.degree();
        self.tab_v_knots = mult_c.knots().to_vec();
        self.tab_v_mults = mult_c.multiplicities().to_vec();
    }

    /// The pole / weight unpacking of `gxx:563-591` (also `:1019-1048`):
    /// divide each pole by its weight and store the weight separately.
    ///
    /// OCCT answers `done = false` and returns when a weight is below
    /// `gp::Resolution()`; the port returns `Ok(false)` for the same early
    /// exit.
    fn unpack_poles(
        &mut self,
        mult_c: &MultiBSpCurve,
        nb_upoles: usize,
        nb_poles2d: usize,
    ) -> Result<bool, String> {
        let nb_vpoles = mult_c.base().nb_poles();
        self.tab_poles = vec![vec![GpPnt::new(0.0, 0.0, 0.0); nb_vpoles]; nb_upoles];
        self.tab_weights = vec![vec![0.0f64; nb_vpoles]; nb_upoles];

        for j in 1..=nb_upoles {
            let newtab_p = mult_c.base().curve3d(j as i32);
            let newtab_p2d = mult_c
                .base()
                .curve2d((j + nb_upoles + nb_poles2d) as i32);
            for k in 1..=nb_vpoles {
                // pour les courbes rationnelles il faut maintenant diviser les
                // poles par leurs poids respectifs
                let a_weight = newtab_p2d[k - 1].x();
                let p = newtab_p[k - 1];
                self.tab_poles[j - 1][k - 1] =
                    GpPnt::new(p.x() / a_weight, p.y() / a_weight, p.z() / a_weight);
                if a_weight < REAL_SMALL {
                    self.done = false;
                    return Ok(false);
                }
                self.tab_weights[j - 1][k - 1] = a_weight;
            }
        }

        for j in 1..=nb_poles2d {
            let poles = mult_c.base().curve2d((nb_upoles + j) as i32);
            self.seq_poles2d.push(poles);
        }
        Ok(true)
    }
}


