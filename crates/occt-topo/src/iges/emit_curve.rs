//! Curve entity emitters. Source: `GeomToIGES_GeomCurve` (line 110, point 116, circle 100, B-spline/Bezier 126).
use super::prelude::*;
use super::*;

impl IgesWriter {

    // ---- geometry entities ----

    pub(super) fn emit_point(&mut self, p: &GpPnt) -> usize {
        self.emit(
            116,
            0,
            format!("116,{},{},{};", num(p.x()), num(p.y()), num(p.z())),
        )
    }

    pub(super) fn emit_line(&mut self, a: &GpPnt, b: &GpPnt) -> usize {
        self.emit(
            110,
            0,
            format!(
                "110,{},{},{},{},{},{};",
                num(a.x()),
                num(a.y()),
                num(a.z()),
                num(b.x()),
                num(b.y()),
                num(b.z())
            ),
        )
    }

    /// Entity 108 (`GeomToIGES_GeomSurface::TransferSurface(Geom_Plane)`,
    /// `GeomToIGES_GeomSurface.cxx:609-688`, the `write.iges.plane.mode == 0`
    /// arm - that static is initialised to 0, `IGESData.cxx:182-186`; written by
    /// `IGESGeom_ToolPlane::WriteOwnParams`, `IGESGeom_ToolPlane.cxx:129-145`):
    ///
    /// `108, A, B, C, D, bounding_curve, attach_x, attach_y, attach_z, symbol_size;`
    ///
    /// - `A,B,C,D` are `start->Coefficients(...)` with `D` negated (`cxx:625-627`)
    ///   - the plane as `A*x + B*y + C*z = D`, all divided by the file unit,
    ///   which is what the port's `d = n . location` already is;
    /// - the bounding curve is null and the symbol size is 0 (`cxx:630` passes
    ///   `nullptr` and `0`), so both write `0`;
    /// - the attachment point is the plane's location / unit (`cxx:629`);
    /// - the form stays 0: `TransferPlaneSurface` never calls `SetFormNumber`,
    ///   and `IGESGeom_Plane::Init` only re-applies the form already stored
    ///   (`IGESGeom_Plane.cxx:44`), which is 0 on a fresh entity.
    pub(super) fn emit_plane(&mut self, pln: &GpPln) -> usize {
        // `GeomToIGES_GeomSurface.cxx:625-629`: `start->Coefficients(A,B,C,D)`
        // then `D = -D` "because of difference in Geom_Plane class and Type 108".
        // `gp_Pln::Coefficients` (`gp_Pln.hxx`) negates the stored main direction
        // when the placement is indirect, which is exactly what
        // `Geom_Plane::UReverse` (`gp_Pln.hxx:101`, `XReverse`) produces for a
        // REVERSED face.
        let (a, b, c, cte) = pln.coefficients();
        let d = -cte;
        let loc = pln.location().coord;
        self.emit(
            108,
            0,
            format!(
                "108,{},{},{},{},0,{},{},{},0.;",
                num(a),
                num(b),
                num(c),
                num(d),
                num(loc.x),
                num(loc.y),
                num(loc.z)
            ),
        )
    }

    /// Entity 100 (`GeomToIGES_GeomCurve::TransferCurve(Geom_Circle)`,
    /// `GeomToIGES_GeomCurve.cxx:533-603`, written by
    /// `IGESGeom_ToolCircularArc::WriteOwnParams`,
    /// `IGESGeom_ToolCircularArc.cxx:79-89`):
    ///
    /// `100, ZT, Xc, Yc, Xs, Ys, Xe, Ye;`
    ///
    /// Every coordinate is **2-D, in the arc's own plane**: `TransferCircle`
    /// places the builder on `gp_Ax3(arc.Position())` (`cxx:560-563`) and writes
    /// each point through `Build.EvalXYZ` (`cxx:580-583`). The placing frame
    /// itself rides on the entity's transformation matrix (`cxx:591-597`).
    pub(super) fn emit_circular_arc(&mut self, frame: &GpAx3, radius: f64, u1: f64, u2: f64) -> usize {
        // `TransferCurve` (`GeomToIGES_GeomCurve.cxx:569-586`) evaluates the
        // centre and the two ends with `Geom_Circle::D0`, then
        // `IGESConvGeom_GeomBuilder::EvalXYZ` (`GeomBuilder.cxx:212-216`):
        // `thepos.Inverted().Transforms` after `SetPosition(gp_Ax3)`
        // (`:137-143`). `gp_Trsf::Invert` transposes (`gp_Trsf.cxx:418-425`).
        // Writing `R*(cos U, sin U)` skips that round trip; a frame whose axes
        // carry a sub-resolution component cancels the libm residual to a real
        // 0 only after the inverse is applied. `GetUnit()` is 1 for this writer
        // (the 124 below is emitted with unit 1).
        let local = builder_local_trsf(frame);
        let (xc, yc, zc) = eval_xyz(&local, &frame.location());
        let (xs, ys, _) = eval_xyz(&local, &circle_point(frame, radius, u1));
        let (xe, ye, _) = eval_xyz(&local, &circle_point(frame, radius, u2));
        let de = self.emit(
            100,
            0,
            format!(
                "100,{},{},{},{},{},{},{};",
                num(zc),
                num(xc),
                num(yc),
                num(xs),
                num(ys),
                num(xe),
                num(ye)
            ),
        );
        if !frame_is_identity(frame) {
            let trsf = self.emit_transformation_matrix(frame, 1.0);
            self.set_trsf(de, trsf);
        }
        de
    }

    /// `GeomToIGES_GeomCurve::TransferCurve(Geom_Circle, Udeb, Ufin)`
    /// (`GeomToIGES_GeomCurve.cxx:533-603`): before the arc is built a `Udeb`
    /// within `gp::Resolution()` snaps to 0 (`:549-552`), and a request spanning
    /// a whole period writes the end point coincident with the start
    /// (`:571-576`, `gka BUG 6542`).
    pub(super) fn emit_circle_transfer(&mut self, frame: &GpAx3, radius: f64, u1: f64, u2: f64) -> usize {
        let first = if u1.abs() <= occt_core::precision::REAL_SMALL {
            0.0
        } else {
            u1
        };
        let last = if (u2 - u1 - 2.0 * std::f64::consts::PI).abs()
            <= occt_core::precision::PCONFUSION
        {
            first
        } else {
            u2
        };
        self.emit_circular_arc(frame, radius, first, last)
    }

    // ---- topology entities ----

    /// `TransferCurve(Geom_OffsetCurve)` under the default
    /// `write.iges.offset.mode = 0` (`IGESData.cxx:193`,
    /// `GeomToIGES_GeomCurve.cxx:916-919`):
    /// `TransferCurve(GeomConvert::CurveToBSplineCurve(start), U1, U2)`.
    ///
    /// `TransferCurve(Geom_TrimmedCurve)` forwards to the immediate basis with
    /// the same range (`cxx:467-475`; the nested-trim assignment at `:472` is
    /// overwritten, so only the basis that recursion reaches is written). The
    /// port's `GeomTrimmedCurve` evaluates on `[0, 1]` while
    /// `untrimmed_basis` reports the basis interval (`trimmed.rs`), so a
    /// remapped trim's request is mapped back before the spline is segmented.
    /// The mode `1` arm (entity 130) is not the default and is not taken.
    /// A failed approximation is the `Standard_ConstructionError` OCCT throws
    /// (`GeomConvert.cxx:350-352`); the caller then keeps its null-handle path.
    pub(super) fn emit_offset_curve_default(
        &mut self,
        curve: &dyn Curve,
        a: f64,
        b: f64,
    ) -> Option<usize> {
        let (target, u0, u1) = offset_curve_request(curve, a, b)?;
        let bs = occt_geom::convert_bspl::curve_to_bspline_curve(
            target.as_ref(),
            occt_core::convert::ParameterisationType::TgtThetaOver2,
        )
        .ok()?;
        self.emit_bspline_curve(&bs, u0, u1)
    }

    /// `GeomToIGES_GeomCurve::TransferCurve(Geom_Curve, Udeb, Ufin)`
    /// (`GeomToIGES_GeomCurve.cxx:94-126` dispatching, `:133-161` on the bounded
    /// family, `:481-528` on the conics): a `Geom_Line` becomes entity 110, a
    /// `Geom_Circle` entity 100, a conic entity 104, a B-spline/Bezier entity 126.
    /// Returns `None` where OCCT's transfer returns a null handle (a narrow
    /// request on a Bezier needs the unported `Geom_BezierCurve::Segment`; any
    /// other curve type has no branch at all).
    ///
    /// A `Geom_OffsetCurve` (and a trimmed curve whose basis is one,
    /// `cxx:467-475`) is taken first, matching the default
    /// `write.iges.offset.mode = 0` arm (`cxx:916-919`).
    pub(super) fn emit_curve_range(&mut self, curve: &dyn Curve, a: f64, b: f64) -> Option<usize> {
        if let Some(idx) = self.emit_offset_curve_default(curve, a, b) {
            return Some(idx);
        }
        match iges_curve_kind(curve) {
            IgCurveKind::Line => {
                let p1 = curve.d0(a);
                let p2 = curve.d0(b);
                Some(self.emit_line(&p1, &p2))
            }
            // `GeomToIGES_GeomCurve::TransferCircle`
            // (`GeomToIGES_GeomCurve.cxx:533-603`): the centre, the axis and the
            // radius come from the `Geom_Circle`'s own `gp_Circ` - never from a
            // three-point fit on sampled points - and the start/end parameters are
            // used as they are, apart from the two special cases in
            // `emit_circle_transfer`.
            IgCurveKind::Circle => {
                let c = curve.gp_circ()?;
                Some(self.emit_circle_transfer(&c.pos.to_ax3(), c.radius, a, b))
            }
            // `TransferCurve(Geom_BSplineCurve)` (`cxx:279-423`): a periodic curve
            // is written through a `SetNotPeriodic` copy, and a requested range
            // narrower than the curve's own through a `Segment`-ed copy, both
            // inside `emit_bspline_curve` (`cxx:294-357`). The remaining
            // UNPORTED case is a narrow request on a **Bezier** curve, whose
            // `Geom_BezierCurve::Segment` is not ported, see there.
            IgCurveKind::Bounded => self.emit_bspline_curve(curve, a, b),
            // `TransferCurve(Geom_Conic)` (`GeomToIGES_GeomCurve.cxx:483-528`:
            // the circle is `:533-603`, the ellipse / hyperbola / parabola
            // transfers are `:608-700`, `:707-773`, `:780-845`) writes entity 104.
            // Any other curve type has no branch in `TransferCurve` and yields a
            // null handle; the caller then writes the chord-line stand-in
            // (`emit_edge_curve`).
            IgCurveKind::Conic => self.emit_conic_arc(curve, a, b),
            IgCurveKind::Other => None,
        }
    }

    /// `BRepToIGES_BRWire::TransferEdge(edge, originMap, false)`
    /// (`BRepToIGES_BRWire.cxx:266-335`): the edge's 3-D curve entity.
    ///
    /// Returns `None` where OCCT leaves `ICurve` null and the caller skips the
    /// edge (`TransferWire`, `cxx:715-723`): the edge has no 3-D curve at all
    /// (a degenerated edge - `BRep_Builder::Degenerated` drops the curve,
    /// `BRep_Builder.cxx:1073-1085`, so `BRep_Tool::Curve` at `cxx:282` returns
    /// null; the port's `EdgeGeom::curve` is non-nullable and carries the
    /// `degenerated` flag instead, see `builder.rs:326-331`).
    pub(super) fn emit_edge_curve(&mut self, e: &Edge) -> Option<usize> {
        let Some(curve) = BRepTool::edge_curve(e) else {
            let (v1, v2) = edge_vertices(e);
            let p1 = v1.map(|v| BRepTool::vertex_point(&v)).unwrap_or_default();
            let p2 = v2.map(|v| BRepTool::vertex_point(&v)).unwrap_or_default();
            return Some(self.emit_line(&p1, &p2));
        };
        if BRepTool::is_degenerated(e) {
            return None;
        }
        let (first, last) = BRepTool::edge_parameters(e);
        // `BRepToIGES_BRWire::TransferEdge` (`cxx:302-314`): every call from the
        // face path passes `theIsBRepMode = false`, so a REVERSED edge is written
        // through `Curve3d->Reverse()` with the range mapped by
        // `ReversedParameter(Last)` / `ReversedParameter(First)`.
        let (curve, a, b): (Arc<dyn Curve>, f64, f64) =
            if e.0.orientation() == crate::abs::Orientation::Reversed {
                let u1 = curve.reversed_parameter(last);
                let u2 = curve.reversed_parameter(first);
                (Arc::from(curve.reversed()), u1, u2)
            } else {
                (curve, first, last)
            };
        match self.emit_curve_range(curve.as_ref(), a, b) {
            Some(idx) => Some(idx),
            None => {
                // The transfer above returned a null handle: OCCT writes no curve
                // at all for those types, the port keeps the chord-line stand-in
                // (the chord-line stand-in, see `emit_curve_range`).
                let p1 = curve.d0(a);
                let p2 = curve.d0(b);
                Some(self.emit_line(&p1, &p2))
            }
        }
    }

    /// Entity 126 (`GeomToIGES_GeomCurve::TransferCurve(Geom_BSplineCurve)`,
    /// `GeomToIGES_GeomCurve.cxx:279-423`, written by
    /// `IGESGeom_ToolBSplineCurve::WriteOwnParams`,
    /// `IGESGeom_ToolBSplineCurve.cxx:64-105`):
    ///
    /// `126, index, degree, planar, closed, polynomial, periodic,
    ///   knots[-degree .. index+1], weights[0 .. index], poles[0 .. index],
    ///   UMin, UMax, Normal;`
    ///
    /// with `index = nb_poles - 1` and the knot array the curve's flattened knot
    /// sequence. The planar flag and the normal come from `IsPlanar(start, Norm)`
    /// (`cxx:414` -> `:236-272`), which runs on the curve **as received**, not on
    /// the unperiodized / `Segment`-ed copy whose poles are written out: for a
    /// B-spline that is `ArePolesPlanar(Poles())` (`cxx:191-228`),
    /// `P(n) x P(1) + sum P(i) x P(i+1)`, normalised (or `(0,0,1)` with
    /// `planar = false` when the sum is shorter than `Precision::Confusion()`),
    /// then every pole must sit at the same distance from the plane through
    /// `P(1)`. Returns `None` for the cases this port cannot express.
    pub(super) fn emit_bspline_curve(&mut self, start: &dyn Curve, u_deb: f64, u_fin: f64) -> Option<usize> {
        // `cxx:294-307`: a periodic curve is written through a non-periodic copy
        // (`SetNotPeriodic`). `IPerio` is captured from the curve as received
        // (`cxx:297`) and passed to `Init` (`cxx:419`) unchanged, so the 126
        // parameter stays 1 after the copy is opened.
        let unperiodized;
        let curve: &dyn Curve = if start.is_periodic() {
            let poles = start.bspline_poles()?.to_vec();
            let knots = start.bspline_knots()?.to_vec();
            let degree = start.nurbs_degree()?;
            let mut c = occt_geom::bspline_curve::GeomBSplineCurve {
                poles,
                weights: start.bspline_weights().map(|w| w.to_vec()),
                knots,
                degree,
                periodic: true,
            };
            c.set_not_periodic();
            unperiodized = c;
            &unperiodized
        } else {
            start
        };

        // `cxx:309-318`: an infinite bound is replaced by `±Precision::Infinite()`.
        let mut umin = if occt_core::precision::Precision::is_negative_infinite(u_deb) {
            -occt_core::precision::INFINITE
        } else {
            u_deb
        };
        let mut umax = if occt_core::precision::Precision::is_positive_infinite(u_fin) {
            occt_core::precision::INFINITE
        } else {
            u_fin
        };

        // `cxx:320-331`: protect against exceptions in `Segment()`; the range is
        // pulled onto the curve's own bounds (one-sided, as OCCT writes the
        // clamped values into the entity's `UMin`/`UMax`).
        let (curve_first, curve_last) = (curve.first_parameter(), curve.last_parameter());
        if umin - curve_first < occt_core::precision::PCONFUSION {
            umin = curve_first;
        }
        if curve_last - umax < occt_core::precision::PCONFUSION {
            umax = curve_last;
        }

        // `cxx:332-356`: cut the curve for E3 — a copy is `Segment`-ed and
        // replaces the original. OCCT catches the failure and keeps the
        // untrimmed copy; a null copy likewise leaves the curve unchanged.
        let trimmed;
        let curve: &dyn Curve = if umin - curve_first > occt_core::precision::PCONFUSION
            || curve_last - umax > occt_core::precision::PCONFUSION
        {
            // UNPORTED: `TransferCurve(Geom_BezierCurve)` (`cxx:441-448`) puts
            // the Bezier in a `Geom_TrimmedCurve` and converts it through
            // `GeomConvert::CurveToBSplineCurve`, which needs
            // `Geom_BezierCurve::Segment` (not ported) — a narrow request on a
            // Bezier keeps returning `None` (the caller writes a chord line).
            let bs = bspline_copy_of(curve)?;
            let mut bs = bs;
            if (umax - umin).abs() > occt_core::precision::PCONFUSION {
                // default `theTolerance = Precision::PConfusion()`.
                let _ = bs.segment(umin, umax, occt_core::precision::PCONFUSION);
            }
            trimmed = bs;
            &trimmed
        } else {
            curve
        };

        let (poles, knots, degree) = if let (Some(p), Some(k), Some(d)) = (
            curve.bspline_poles(),
            curve.bspline_knots(),
            curve.nurbs_degree(),
        ) {
            (p.to_vec(), k.to_vec(), d)
        } else if let Some(p) = curve.bezier_poles() {
            // `TransferCurve(Geom_BezierCurve)` (`cxx:430-450`) goes through
            // `GeomConvert::CurveToBSplineCurve`: one span whose end knots carry
            // multiplicity degree+1.
            let d = p.len().checked_sub(1)?;
            let mut k = vec![curve_first; d + 1];
            k.extend(std::iter::repeat(curve_last).take(d + 1));
            (p.to_vec(), k, d)
        } else {
            return None;
        };
        if poles.len() < 2 || degree == 0 || knots.len() != poles.len() + degree + 1 {
            return None;
        }
        let index = poles.len() - 1;


        let (planar, normal) = curve_planar_and_normal(start);
        // `cxx:415-418`: the normal is flipped when it points down.
        let normal = if normal.z() < 0.0 {
            GpVec::new(-normal.x(), -normal.y(), -normal.z())
        } else {
            normal
        };

        // `cxx:360`: `IPolyn = !mycurve->IsRational()`. `IsRational` is the
        // `myRational` flag (`Geom_BSplineCurve_1.cxx:683-686`), set by
        // `Rational(Weights)` (`Geom_BSplineCurve.cxx:98-108`, `:206-208`).
        // A stored weight array whose values are all equal is not rational, and
        // the entity then carries `BSplCLib::UnitWeights` (`:220-222`).
        //
        // `WriteOwnParams` does not send that flag. It sends
        // `IGESGeom_BSplineCurve::IsPolynomial()` (`IGESGeom_BSplineCurve.cxx:105-124`,
        // called from `IGESGeom_ToolBSplineCurve.cxx:242`), which reports
        // polynomial when every weight is within `1e-10` of the first one.
        let stored = curve.bspline_weights().filter(|w| w.len() == poles.len());
        let geom_rational = stored.as_ref().is_some_and(|w| weights_are_rational(w));
        let weights: Vec<f64> = if geom_rational {
            stored.unwrap().to_vec()
        } else {
            vec![1.0; poles.len()]
        };
        let polynomial = !geom_rational
            || weights
                .first()
                .is_none_or(|w0| weights.iter().all(|wi| (wi - w0).abs() <= 1.0e-10));
        let closed = poles
            .first()
            .zip(poles.last())
            .map(|(f, l)| f.distance(l) <= occt_core::precision::CONFUSION)
            .unwrap_or(false);

        let mut s = format!(
            "126,{index},{degree},{},{},{},{}",
            u8::from(planar),
            u8::from(closed),
            u8::from(polynomial),
            u8::from(start.is_periodic())
        );
        for k in &knots {
            s.push(',');
            s.push_str(&num(*k));
        }
        for w in &weights {
            s.push(',');
            s.push_str(&num(*w));
        }
        for p in &poles {
            s.push_str(&format!(",{},{},{}", num(p.x()), num(p.y()), num(p.z())));
        }
        // `cxx:419`: the entity carries the clamped/segmented `Umin`/`Umax`,
        // not the raw requested range.
        s.push_str(&format!(
            ",{},{},{},{},{};",
            num(umin),
            num(umax),
            num(normal.x()),
            num(normal.y()),
            num(normal.z())
        ));
        Some(self.emit(126, 0, s))
    }

}
