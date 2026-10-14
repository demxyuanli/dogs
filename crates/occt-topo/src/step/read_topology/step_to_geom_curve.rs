//! `StepToGeom::MakeCurve`: the 3D curve arm dispatch into lines, conics, B-splines and the composite / trimmed forms.
use super::*;

impl<'a> Resolver<'a> {
    pub(in crate::step) fn resolve_curve(&self, id: usize) -> Result<Arc<dyn Curve>, String> {
        if let Some(c) = self.curve_cache.borrow().get(&id) {
            return Ok(c.clone());
        }
        let rec = self.record(id)?;
        let curve: Arc<dyn Curve> = match rec.type_name.as_str() {
            "LINE" => {
                let pnt = parse_ref(&rec.args[1]).ok_or("LINE: bad point ref")?;
                let vec = parse_ref(&rec.args[2]).ok_or("LINE: bad vector ref")?;
                let p = self.resolve_point(pnt)?;
                let v = self.resolve_vector(vec)?;
                let d = GpDir::from_vec(&v).map_err(|e| format!("LINE: {e}"))?;
                Arc::new(GeomLine::new(GpLin::from_pnt_dir(p, d)))
            }
            "CIRCLE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("CIRCLE: bad axis ref")?;
                // `StepToGeom::MakeCircle` (`StepToGeom.cxx:1212-1225`):
                // `SC->Radius() * LengthFactor()` (`cxx:1222`).
                let r = parse_f64(&rec.args[2])? * self.length_factor;
                Arc::new(GeomCircle::new(GpCirc::new(self.resolve_axis2(ax)?, r)))
            }
            "ELLIPSE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("ELLIPSE: bad axis ref")?;
                // `StepToGeom::MakeEllipse` (`StepToGeom.cxx:1536-1566`):
                // `majorR = SemiAxis1 * LF`, `minorR = SemiAxis2 * LF`
                // (`cxx:1549-1550`); when `majorR - minorR < 0` OCCT turns the X
                // direction by `A.XDirection() ^ A.Direction()` and swaps the radii
                // (`cxx:1557-1561`).
                let maj = parse_f64(&rec.args[2])? * self.length_factor;
                let min = parse_f64(&rec.args[3])? * self.length_factor;
                let mut ax2 = self.resolve_axis2(ax)?;
                let (maj, min) = if maj - min >= 0.0 {
                    (maj, min)
                } else {
                    if let Ok(xd) = ax2.x_direction().crossed(&ax2.direction()) {
                        ax2.set_x_direction(xd);
                    }
                    (min, maj)
                };
                Arc::new(GeomEllipse::new(GpElips::new(ax2, maj, min)))
            }
            "HYPERBOLA" => {
                let ax = parse_ref(&rec.args[1]).ok_or("HYPERBOLA: bad axis ref")?;
                // `StepToGeom::MakeHyperbola` (`StepToGeom.cxx:1604-1622`):
                // `SemiAxis() * LF`, `SemiImagAxis() * LF` (`cxx:1616-1617`).
                let maj = parse_f64(&rec.args[2])? * self.length_factor;
                let min = parse_f64(&rec.args[3])? * self.length_factor;
                Arc::new(GeomHyperbola::new(GpHypr::new(
                    self.resolve_axis2(ax)?,
                    maj,
                    min,
                )))
            }
            "PARABOLA" => {
                let ax = parse_ref(&rec.args[1]).ok_or("PARABOLA: bad axis ref")?;
                // `StepToGeom::MakeParabola` (`StepToGeom.cxx:1695-1710`):
                // `SC->FocalDist() * LengthFactor()` (`cxx:1706`).
                let f = parse_f64(&rec.args[2])? * self.length_factor;
                Arc::new(GeomParabola::new(GpParab::new(self.resolve_axis2(ax)?, f)))
            }
            "SURFACE_CURVE" | "SEAM_CURVE" => {
                // SURFACE_CURVE/SEAM_CURVE(name, curve_3d, pcurves, master_rep):
                // the 3D curve is the second argument; the pcurve list
                // (associated_geometry) is carried for the face-level pcurve
                // association. SEAM_CURVE is the seam of a closed surface.
                let c3d = parse_ref(&rec.args[1]).ok_or("SURFACE_CURVE: bad 3D curve ref")?;
                self.surface_curve_pcurves
                    .borrow_mut()
                    .insert(id, parse_ref_list(&rec.args[2]));
                self.resolve_curve(c3d)?
            }
            "BEZIER_CURVE" | "UNIFORM_CURVE" | "QUASI_UNIFORM_CURVE" => {
                // `StepToGeom::MakeBSplineCurve` (`StepToGeom.cxx:295-450`): the STEP
                // Bezier / uniform / quasi-uniform curve is converted into a
                // `BSplineCurveWithKnots` before being mapped:
                //   Bezier       -> knots {0,1}, multiplicities degree+1 (`:310-317`)
                //   Uniform      -> n_poles + degree + 1 knots i-1, mults 1 (`:338-347`)
                //   QuasiUniform -> n_poles - degree + 1 knots i-1, mults 1 except
                //                   both ends = degree+1 (`:362-384`)
                // The `_AND_RATIONAL_B_SPLINE_CURVE` complex forms add the weights
                // (`:385-421`, `:422-458`); a rational **Bezier** curve is the
                // exception — OCCT has no Bezier+rational arm, and its plain
                // `BezierCurve` arm (`:295-320`) never reads weights.
                let degree = parse_f64(&rec.args[1])? as usize;
                let poles: Vec<GpPnt> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point(r))
                    .collect::<Result<Vec<_>, _>>()?;
                // `MakeBSplineCurveCommon` synthesises the knot / multiplicity
                // lists this family implies (`StepToGeom.cxx:310-317`,
                // `:338-347`, `:362-384`) and then runs the shared periodic /
                // closed-curve logic.
                let n = poles.len();
                let (mults, knots): (Vec<usize>, Vec<f64>) = match rec.type_name.as_str() {
                    "BEZIER_CURVE" => {
                        (vec![degree + 1, degree + 1], vec![0.0, 1.0])
                    }
                    "UNIFORM_CURVE" => {
                        let nb = n + degree + 1;
                        (vec![1; nb], (0..nb).map(|i| i as f64).collect())
                    }
                    _ => {
                        let nb = n.saturating_sub(degree) + 1;
                        let mut m = vec![1usize; nb];
                        m[0] = degree + 1;
                        m[nb - 1] = degree + 1;
                        (m, (0..nb).map(|i| i as f64).collect())
                    }
                };
                // A merged rational complex carries the weights at index 3; a plain
                // entity has `curve_form` there (a marker like `.UNSPECIFIED.`).
                let weights_arg = rec.args.get(3).map(|s| s.trim().to_string());
                let weights = match weights_arg.as_deref() {
                    Some(w) if w.starts_with('(') && rec.type_name != "BEZIER_CURVE" => {
                        Some(parse_real_list(w))
                    }
                    _ => None,
                };
                if let Some(w) = weights.as_ref() {
                    if w.len() != poles.len() {
                        return Err(format!("{}: weight count mismatch", rec.type_name));
                    }
                }
                let curve = make_bspline_curve_with_knots(
                    degree,
                    poles,
                    weights,
                    &mults,
                    &knots,
                    curve_record_closed(rec),
                )?;
                Arc::new(curve)
            }
            "CURVE_REPLICA" => {
                // `StepToGeom::MakeCurve` CurveReplica arm (`StepToGeom.cxx:1351-1371`):
                // `C1 = MakeCurve(ParentCurve)`, then `C1->Transform(T1)` with
                // `T1 = MakeTransformation3d(Transformation)`. The guard
                // `!T.IsNull() && PC != SC` (`cxx:1358`) rejects a cyclic replica.
                let parent_ref = parse_ref(&rec.args[1]).ok_or("CURVE_REPLICA: bad parent ref")?;
                let trsf_ref =
                    parse_ref(&rec.args[2]).ok_or("CURVE_REPLICA: bad transformation ref")?;
                if parent_ref == id {
                    return Err(format!("CURVE_REPLICA: cyclic parent (#{id})"));
                }
                let parent = self.resolve_curve(parent_ref)?;
                let t = self.make_transformation3d(trsf_ref)?;
                Arc::from(parent.transformed(&t))
            }
            "B_SPLINE_CURVE_WITH_KNOTS" => {
                // Layout (10 args): name, degree, control_points, weights|SELF,
                // curve_form, closed, self_intersect, knots, multiplicities, knot_spec.
                let degree = parse_f64(&rec.args[1])? as usize;
                let poles: Vec<GpPnt> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point(r))
                    .collect::<Result<Vec<_>, _>>()?;
                let weights_arg = rec.args.get(3).map(|s| s.trim().to_string());
                // B_SPLINE_CURVE_WITH_KNOTS layout:
                // (name, degree, control_points, curve_form|weights, curve_form,
                //  closed, self_intersect, knot_multiplicities, knots, knot_spec).
                // `MakeBSplineCurveCommon` (`StepToGeom.cxx:776-927`) merges
                // duplicate knots, clamps the multiplicities, trims the poles and
                // then decides the periodic / closed-curve representation.
                let knot_mults = parse_usize_list(rec.args.get(6).map(|s| s.as_str()).unwrap_or("()"));
                let knot_values = parse_real_list(rec.args.get(7).map(|s| s.as_str()).unwrap_or("()"));
                let weights = match weights_arg.as_deref() {
                    // No weights: non-rational curve (curve_form is UNSPECIFIED
                    // or a non-rational flag like CIRCULAR/LINEAR).
                    None | Some("SELF") | Some(".UNSPECIFIED.") | Some(".CIRCULAR.")
                    | Some(".LINEAR.") => None,
                    Some(w) if w.starts_with('.') => None,
                    Some(w) => {
                        let weights = parse_real_list(w);
                        if weights.len() != poles.len() {
                            return Err("B_SPLINE_CURVE: weight count mismatch".into());
                        }
                        Some(weights)
                    }
                };
                let curve = make_bspline_curve_with_knots(
                    degree,
                    poles,
                    weights,
                    &knot_mults,
                    &knot_values,
                    curve_record_closed(rec),
                )?;
                Arc::new(curve)
            }
            "POLYLINE" => {
                // A connected sequence of CARTESIAN_POINTs; represent it as a
                // degree-1 B-spline that passes through every vertex.
                let pts: Vec<GpPnt> = parse_ref_list(&rec.args[1])
                    .into_iter()
                    .map(|r| self.resolve_point(r))
                    .collect::<Result<Vec<_>, _>>()?;
                if pts.len() < 2 {
                    return Err("POLYLINE: need at least 2 points".into());
                }
                let knots = uniform_knots_for(pts.len(), 1);
                Arc::new(
                    GeomBSplineCurve::new(pts, knots, 1)
                        .map_err(|e| format!("POLYLINE: {e}"))?,
                )
            }
            "B_SPLINE_CURVE" => {
                // Plain B-spline (no explicit knots): the knot vector is the
                // clamped uniform one implied by the pole count and degree.
                let degree = parse_f64(&rec.args[1])? as usize;
                let poles: Vec<GpPnt> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point(r))
                    .collect::<Result<Vec<_>, _>>()?;
                let knots = uniform_knots_for(poles.len(), degree);
                Arc::new(
                    GeomBSplineCurve::new(poles, knots, degree)
                        .map_err(|e| format!("B_SPLINE_CURVE: {e}"))?,
                )
            }
            "TRIMMED_CURVE" => {
                // `StepToGeom::MakeTrimmedCurve` (`StepToGeom.cxx:2323-2496`) with
                // `ExtractParameter` (`:2221-2317`).
                // Layout: name, basis_curve, trim_1, trim_2, sense_agreement,
                // master_representation.
                let basis_ref = parse_ref(&rec.args[1]).ok_or("TRIMMED_CURVE: bad basis ref")?;
                let basis_rec = self.record(basis_ref)?;
                let basis_kind = basis_rec.type_name.clone();
                let basis = self.resolve_curve(basis_ref)?;
                let (params1, points1) = parse_trimming_select(&rec.args[2]);
                let (params2, points2) = parse_trimming_select(&rec.args[3]);
                let sense = !matches!(
                    rec.args.get(4).map(|s| s.trim()),
                    Some(".F.") | Some("F") | Some("false")
                );
                // `cxx:2339-2350`: `.CARTESIAN.` = 1, `.PARAMETER.` = 2, else 0.
                let master_rep: i32 = match rec.args.get(5).map(|s| s.trim()) {
                    Some(".CARTESIAN.") => 1,
                    Some(".PARAMETER.") => 2,
                    _ => 0,
                };
                // `cxx:2352-2373`: with an unspecified master representation (or a
                // parameter one carrying two selects on each side), a Cartesian
                // point on both trims means the trims are points, not parameters.
                let n1 = params1.len() + points1.len();
                let n2 = params2.len() + points2.len();
                let is_point = (master_rep == 0 || (master_rep == 2 && n1 > 1 && n2 > 1))
                    && !points1.is_empty()
                    && !points2.is_empty();
                // `cxx:2375-2392`: parameter scaling per basis type.
                let (fact, shift) = match basis_kind.as_str() {
                    "LINE" => {
                        let vec_ref =
                            parse_ref(&basis_rec.args[2]).ok_or("LINE: bad vector ref")?;
                        // `cxx:2380`: `Dir()->Magnitude() * LengthFactor` — the port's
                        // `resolve_vector` already applies `LengthFactor`
                        // (`StepToGeom::MakeVectorWithMagnitude`, `cxx:2577`).
                        (self.resolve_vector(vec_ref)?.xyz().modulus(), 0.0)
                    }
                    "CIRCLE" | "ELLIPSE" => {
                        // `cxx:2386`: `PlaneAngleFactor`; `cxx:2387-2392`: a π/2 shift
                        // for an ellipse whose `SemiAxis1 - SemiAxis2 < 0`.
                        let mut shift = 0.0;
                        if basis_kind == "ELLIPSE" {
                            let s1 = parse_f64(&basis_rec.args[2])?;
                            let s2 = parse_f64(&basis_rec.args[3])?;
                            if s1 - s2 < 0.0 {
                                shift = 0.5 * PI;
                            }
                        }
                        (self.plane_angle_factor, shift)
                    }
                    _ => (1.0, 0.0),
                };
                // `cxx:2394-2425`: a conic whose placement has no reference direction
                // cannot be trimmed by parameters; OCCT returns the full period with
                // the sense of `SenseAgreement` (the port's `GeomTrimmedCurve` carries
                // no sense flag, so the full-period span is what is modelled).
                if matches!(basis_kind.as_str(), "CIRCLE" | "ELLIPSE") && !is_point && master_rep != 1 {
                    let ax_ref = parse_ref(&basis_rec.args[1]).ok_or("conic: bad axis ref")?;
                    if !self.axis2_has_ref_direction(ax_ref)? {
                        return Ok(Arc::new(GeomTrimmedCurve::new(basis, 0.0, 2.0 * PI)));
                    }
                }
                let t1 = self.extract_parameter(&*basis, &params1, &points1, master_rep, fact, shift)?;
                let t2 = self.extract_parameter(&*basis, &params2, &points2, master_rep, fact, shift)?;
                let (Some(mut trim1), Some(mut trim2)) = (t1, t2) else {
                    // `cxx:2495`: no parameter on either side → OCCT returns null.
                    return Err(format!("TRIMMED_CURVE: no trimming parameters (#{id})"));
                };
                let cf = basis.first_parameter();
                let cl = basis.last_parameter();
                // `cxx:2438-2459`: clamp into the basis range when it is not periodic.
                if !basis.is_periodic() {
                    trim1 = trim1.clamp(cf, cl);
                    trim2 = trim2.clamp(cf, cl);
                }
                if (trim1 - trim2).abs() < occt_core::precision::PCONFUSION {
                    if basis.is_periodic() {
                        // `cxx:2462-2465`: `ElCLib::AdjustPeriodic(cf, cl, PConfusion, …)`.
                        occt_core::elib::clib2d::adjust_periodic(
                            cf,
                            cl,
                            occt_core::precision::PCONFUSION,
                            &mut trim1,
                            &mut trim2,
                        );
                    } else {
                        // `cxx:2466-2480`: the closed (non-periodic) basis arm is
                        // UNPORTED — the `Curve` trait exposes no `IsClosed`, so the
                        // null-result branch OCCT takes otherwise is returned here.
                        return Err(format!("TRIMMED_CURVE: degenerate trim on closed basis (#{id})"));
                    }
                }
                // `cxx:2486-2493`: `SenseAgreement` selects the argument order. The
                // port's `GeomTrimmedCurve` stores only the span (no `Sense` flags,
                // `trimmed.rs`), so an agreed trim is written as-is and an opposite one
                // is written reversed — the reversal flags themselves are UNPORTED.
                if sense {
                    Arc::new(GeomTrimmedCurve::new(basis, trim1, trim2))
                } else {
                    Arc::new(GeomTrimmedCurve::new(basis, trim2, trim1))
                }
            }
            "OFFSET_CURVE_3D" => {
                // Layout: name, basis_curve, direction, distance, self_intersect,
                // curve_form.
                let basis_ref = parse_ref(&rec.args[1]).ok_or("OFFSET_CURVE_3D: bad curve ref")?;
                let dir_ref = parse_ref(&rec.args[2]).ok_or("OFFSET_CURVE_3D: bad dir ref")?;
                let offset = parse_f64(&rec.args[3])?;
                let basis = self.resolve_curve(basis_ref)?;
                let dir = self.resolve_direction(dir_ref)?;
                Arc::new(GeomOffsetCurve::new(basis, offset, dir))
            }
            other => {
                self.warn(format!("unsupported curve entity {other} (#{id})"));
                return Err(format!("unsupported curve entity {other} (#{id})"));
            }
        };
        self.curve_cache.borrow_mut().insert(id, curve.clone());
        Ok(curve)
    }

}
