//! Conic-arc (104), whole-period ellipse and offset-surface emitters, plus the face -> basis-surface bridge. Source: `GeomToIGES_GeomCurve` conic arms and `GeomToIGES_GeomSurface`.
use super::prelude::*;
use super::*;

impl IgesWriter {

    /// `GeomToIGES_GeomCurve::TransferCurve(Geom_Ellipse)` full-period arm
    /// (`GeomToIGES_GeomCurve.cxx:620-645`): the ellipse copy is rotated by
    /// `Udeb` (`:626-628`, direction chosen by `gp_Ax3(pos).Direct()`), converted
    /// to a B-spline, re-parameterised onto `[Udeb, Udeb + 2*PI]` (`:641-643`)
    /// and transferred as entity 126.
    ///
    /// `GeomConvert_ApproxCurve(aCopy, Precision::Approximation(), GeomAbs_C1,
    /// 100, 6)` is tried first (`:632-636`); `GeomConvert::CurveToBSplineCurve
    /// (copystart, Convert_QuasiAngular)` (`:637-640`) is only the fallback when
    /// that produces no curve.
    pub(super) fn emit_whole_period_ellipse(&mut self, e: &occt_core::gp::GpElips, a: f64, b: f64) -> Option<usize> {
        use std::f64::consts::PI;
        let pos = *e.position();
        // `copystart->SetPosition(pos.Rotated(pos.Axis(), gp_Ax3(pos).Direct() ? Udeb : 2*PI - Udeb))`
        // (`cxx:626-628`).
        let angle = if pos.to_ax3().is_direct() { a } else { 2.0 * PI - a };
        let mut copy = *e;
        copy.set_position(pos.rotated(&pos.axis().clone(), angle));
        let rotated = occt_geom::ellipse::GeomEllipse::new(copy);
        // `if (approx.HasResult()) Bspline = approx.Curve(); if (Bspline.IsNull())
        //  Bspline = GeomConvert::CurveToBSplineCurve(copystart, Convert_QuasiAngular);`
        // (`cxx:632-640`).
        let approx = occt_geom::GeomConvertApproxCurve::new(
            &rotated,
            occt_core::precision::APPROXIMATION,
            occt_core::kernel::geomabs::Shape::C1,
            100,
            6,
        );
        let mut bs = match approx.curve() {
            Some(curve) => curve.clone(),
            // `GeomConvert::CurveToBSplineCurve(copystart, Convert_QuasiAngular)` (`cxx:639`).
            None => occt_geom::convert_bspl::curve_to_bspline_curve(
                &rotated,
                occt_core::convert::ParameterisationType::QuasiAngular,
            )
            .ok()?,
        };
        // `Knots = Bspline->Knots(); BSplCLib::Reparametrize(Udeb, Udeb + 2*PI, Knots);
        //  Bspline->SetKnots(Knots);` (`cxx:641-643`).
        let (mut uknots, _) = bs.distinct_knots_and_mults();
        occt_core::bspl::knots::reparametrize(a, a + 2.0 * PI, &mut uknots);
        bs.set_knots(&uknots).ok()?;
        // `return TransferCurve(Bspline, Udeb, Ufin)` (`cxx:644`). The periodic
        // flag is read inside that transfer, before `SetNotPeriodic`.
        self.emit_bspline_curve(&bs, a, b)
    }

    /// Entity 104 (`GeomToIGES_GeomCurve::TransferCurve(Geom_Ellipse)`,
    /// `GeomToIGES_GeomCurve.cxx:608-700`; `(Geom_Hyperbola)`, `:707-773`;
    /// `(Geom_Parabola)`, `:780-845`; written by
    /// `IGESGeom_ToolConicArc::WriteOwnParams`, `IGESGeom_ToolConicArc.cxx:103-125`):
    ///
    /// `104, A, B, C, D, E, F, ZT, start.x, start.y, end.x, end.y;`
    ///
    /// The coefficients come from the 2d conic built on the identity frame (see
    /// [`gp_elips2d_coefficients`]) with the radii divided by the file unit, the
    /// DE form number from [`conic_form_number`] (`IGESGeom_ConicArc::Init`,
    /// `IGESGeom_ConicArc.cxx:33-53`), `ZT = 0` and the arc's end points in the
    /// conic's own frame (`Build.EvalXYZ`, `:669-670`). A frame other than the
    /// absolute one is recorded as entity 124 on the DE card
    /// (`:692-697`). A full-period ellipse is routed to
    /// [`emit_whole_period_ellipse`](IgesWriter::emit_whole_period_ellipse)
    /// instead, exactly as OCCT does at `:620-645`.
    pub(super) fn emit_conic_arc(&mut self, curve: &dyn Curve, a: f64, b: f64) -> Option<usize> {
        use std::f64::consts::PI;
        let unit = IGES_UNIT;
        let (abc, pos, u1, u2) = if let Some(e) = curve.gp_ellipse() {
            if (b - a - 2.0 * PI).abs() <= occt_core::precision::PCONFUSION {
                // `GeomToIGES_GeomCurve.cxx:620-645`: a trimmed **full-period**
                // ellipse (IGES 104 cannot carry the whole-period semantics) is
                // re-parameterised onto `[Udeb, Udeb + 2*PI]` and written as a
                // B-spline (entity 126) instead.
                return self.emit_whole_period_ellipse(&e, a, b);
            }
            // `cxx:649-654`: `|Udeb| <= gp::Resolution()` is snapped to 0.
            let u1 = if a.abs() <= occt_core::precision::REAL_SMALL {
                0.0
            } else {
                a
            };
            let (fa, fb, fc, fd, fe, ff) =
                gp_elips2d_coefficients(e.major_radius() / unit, e.minor_radius() / unit);
            // `E2d.Coefficients(A, C, B, D, E, F)` (`cxx:675`) then
            // `Init(A, 2*B, C, 2*D, 2*E, F)` (`cxx:678-686`).
            ([fa, 2.0 * fc, fb, 2.0 * fd, 2.0 * fe, ff], *e.position(), u1, b)
        } else if let Some(h) = curve.gp_hyperbola() {
            let (u1, u2) = occt_infinite_range(a, b);
            let (fa, fb, fc, fd, fe, ff) =
                gp_hypr2d_coefficients(h.major_radius / unit, h.minor_radius / unit);
            // `H2d.Coefficients(A, C, B, D, E, F)` (`cxx:749`) then
            // `Init(A, B, C, D, E, F)` (`cxx:751-759`; no doubling here).
            ([fa, fc, fb, fd, fe, ff], *h.position(), u1, u2)
        } else if let Some(p) = curve.gp_parabola() {
            let (u1, u2) = occt_infinite_range(a, b);
            let (fa, fb, fc, fd, fe, ff) = gp_parab2d_coefficients(p.focal / unit);
            // `P2d.Coefficients(A, C, B, D, E, F)` (`cxx:821`) then
            // `Init(A, B, C, D, E, F)` (`cxx:823-831`; no doubling here).
            ([fa, fc, fb, fd, fe, ff], *p.position(), u1, u2)
        } else {
            return None;
        };
        let frame = pos.to_ax3();
        // `Build.EvalXYZ((start->Value(U1)).XYZ(), ...)` (`cxx:669-670`) followed
        // by `gp_XY(Xs / GetUnit(), Ys / GetUnit())` in `Init` (`cxx:678-686`).
        let (sx, sy, _) = frame_local_point(&frame, &curve.d0(u1));
        let (sx, sy) = (sx / unit, sy / unit);
        let (ex, ey, _) = frame_local_point(&frame, &curve.d0(u2));
        let (ex, ey) = (ex / unit, ey / unit);
        let [fa, fb, fc, fd, fe, ff] = abc;
        let form = conic_form_number(fa, fb, fc, fd, fe, ff);
        let de = self.emit(
            104,
            form,
            format!(
                "104,{},{},{},{},{},{},{},{},{},{},{};",
                num(fa),
                num(fb),
                num(fc),
                num(fd),
                num(fe),
                num(ff),
                num(0.0),
                num(sx),
                num(sy),
                num(ex),
                num(ey)
            ),
        );
        // `cxx:688-697`: the transformation matrix is created **after** the conic
        // arc and recorded on it (`Conic->InitTransf(TMat)`).
        if !frame_is_identity(&frame) {
            let t = self.emit_transformation_matrix(&frame, unit);
            self.set_trsf(de, t);
        }
        Some(de)
    }

    /// Base surface entity for a face.
    ///
    /// `BRepToIGES_BRShell::TransferFace` (`BRepToIGES_BRShell.cxx:246-262`)
    /// transfers the face's surface with
    /// `GeomToIGES_GeomSurface::TransferSurface(Surf, U1, U2, V1, V2)` and, when
    /// that returns a null handle, logs `"the basic surface is a null entity"`
    /// and **drops the face** (`:257-261`). There is no plane stand-in.
    ///
    /// `reversal` is `TransferFace`'s REVERSED branch: OCCT already replaced the
    /// face surface with the U-reversed, trim-stripped one (`:126-146`), mirrored
    /// every p-curve about `U = aCenter` and recomputed the UV bounds on the new
    /// face (`:253`). The port transfers the original face and mirrors the
    /// p-curves on the fly (`transfer_edge_uv`), so the box it hands to
    /// `TransferSurface` has to be the mirror of the original one -
    /// `[2*aCenter - U1, 2*aCenter - U0]`, V untouched - exactly the box OCCT
    /// reads off the mirrored copy. Handing the unmirrored box over instead turns
    /// every elementary arm's `2*PI - U2, 2*PI - U1` into a whole-period shift of
    /// the revolution range (a quarter cylinder at `[0, PI/2]` came out as
    /// `[3*PI/2, 2*PI]`).
    pub(super) fn emit_face_surface(&mut self, f: &Face, reversal: Option<&FaceReversal>) -> Option<usize> {
        let surf = match reversal {
            Some(r) => r.surf.clone(),
            None => BRepTool::face_surface(f)?,
        };
        let (u0, u1, v0, v1) = face_uv_bounds_finite(f);
        let (u0, u1) = match reversal {
            Some(r) => (r.two_center - u1, r.two_center - u0),
            None => (u0, u1),
        };
        self.emit_surface_entity(surf.as_ref(), u0, u1, v0, v1)
    }

    /// `GeomToIGES_GeomSurface::TransferSurface(const Geom_Surface&, Udeb, Ufin,
    /// Vdeb, Vfin)` (`GeomToIGES_GeomSurface.cxx:112-147`): the surface type
    /// dispatcher.
    ///
    /// The dispatch order is OCCT's: the `Geom_RectangularTrimmedSurface` wrapper
    /// is stripped first (`Geom_BoundedSurface` family, `:125-128` -> `:492-515`,
    /// which recurses on `BasisSurface()`), then the exact `Geom_BoundedSurface` /
    /// `Geom_ElementarySurface` / `Geom_SweptSurface` / `Geom_OffsetSurface` types
    /// are tried in that order. No sampling is involved; a surface with no
    /// matching arm has no branch in OCCT either and returns a null handle there
    /// and `None` here.
    pub(super) fn emit_surface_entity(
        &mut self,
        s: &dyn Surface,
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
    ) -> Option<usize> {
        // `TransferSurface(Geom_BoundedSurface)` (`cxx:492-515`): a rectangular
        // trim restarts the whole dispatch on its basis.
        if let Some(b) = s.rectangular_trimmed_basis() {
            return self.emit_surface_entity(b.as_ref(), u0, u1, v0, v1);
        }
        if let Some(pln) = s.gp_pln() {
            return Some(self.emit_plane(&pln));
        }
        if let Some(cy) = s.gp_cylinder() {
            // `cxx:815-817`: the generatrix is `Geom_Line(gp_Pnt(R,0,0), +Z)`, so
            // `Value(t) = (R, 0, t)`; `V1`/`V2` are the clamped `UVBounds` of the
            // face (`cxx:806-813`).
            let (v0, v1) = clamp_infinite_range(v0, v1);
            let r = cy.radius();
            let generatrix = self.emit_line(&GpPnt::new(r, 0.0, v0), &GpPnt::new(r, 0.0, v1));
            return Some(self.emit_local_revolution_surface(
                generatrix,
                &cy.position(),
                u0,
                u1,
            ));
        }
        if let Some(co) = s.gp_cone() {
            // `cxx:816-818`: `Geom_Line(gp_Pnt(RefRadius,0,0),
            // gp_Dir(sin a, 0, cos a))`, so `Value(t) = (RefRadius + t*sin a, 0,
            // t*cos a)`. The semi-angle is used as it is stored, sign included.
            let (v0, v1) = clamp_infinite_range(v0, v1);
            let (semi_sin, semi_cos) = co.semi_angle().sin_cos();
            let ref_radius = co.radius();
            let generatrix = self.emit_line(
                &GpPnt::new(ref_radius + v0 * semi_sin, 0.0, v0 * semi_cos),
                &GpPnt::new(ref_radius + v1 * semi_sin, 0.0, v1 * semi_cos),
            );
            return Some(self.emit_local_revolution_surface(
                generatrix,
                &co.position(),
                u0,
                u1,
            ));
        }
        if let Some(sp) = s.gp_sphere() {
            // `cxx:893-894`: the generatrix is the half circle
            // `Geom_Circle(gp_Ax2(gp::Origin(), -gp::DY(), gp::DX()), R)`
            // travelled over the face's own V range (`Su(-PI/2) .. Sv(+PI/2)`).
            let generatrix = self.emit_circle_transfer(
                &revolution_generatrix_frame(GpPnt::zero()),
                sp.radius(),
                v0,
                v1,
            );
            return Some(self.emit_local_revolution_surface(
                generatrix,
                &sp.position(),
                u0,
                u1,
            ));
        }
        if let Some(t) = s.gp_torus() {
            // `cxx:960-962`: `Geom_Circle(gp_Ax2(gp_Pnt(Major,0,0), -gp::DY(),
            // gp::DX()), Minor)`.
            let generatrix = self.emit_circle_transfer(
                &revolution_generatrix_frame(GpPnt::new(t.major_radius(), 0.0, 0.0)),
                t.minor_radius(),
                v0,
                v1,
            );
            return Some(self.emit_local_revolution_surface(
                generatrix,
                &t.position(),
                u0,
                u1,
            ));
        }
        // `GeomToIGES_GeomSurface::TransferSurface` (`cxx:154-188`) on the
        // Bounded family: the port writes entity 128.
        if let Some(idx) = self.emit_bspline_surface(s, u0, u1, v0, v1) {
            return Some(idx);
        }
        // `cxx:1000-1025` on the swept family - reached after the bounded family
        // (B-spline / Bezier / trimmed), which is the order `cxx:112-147`
        // dispatches in.
        match swept_surface_kind(s) {
            Some(SweptKind::Extrusion) => {
                if let Some(idx) = self.emit_tabulated_cylinder(s, v0, v1) {
                    return Some(idx);
                }
            }
            Some(SweptKind::Revolution) => {
                if let Some(idx) = self.emit_surface_of_revolution(s, u0, u1, v0, v1) {
                    return Some(idx);
                }
            }
            None => {}
        }
        // `cxx:1195-1237` on the offset family (entity 140 `OffsetSurface`). The
        // basis surface is put through this same dispatcher, so the 140 nests
        // whatever entity the basis maps to.
        if let (Some(basis), Some(distance)) = (s.offset_basis_surface(), s.offset_distance()) {
            // `cxx:1208-1214`: the indicator is read at the mid point of the
            // offset's own `Bounds`, which `Geom_OffsetSurface::Bounds` delegates
            // to its basis.
            let (bu0, bu1) = s.u_range();
            let (bv0, bv1) = s.v_range();
            let (um, vm) = (bounds_mid_point(bu0, bu1), bounds_mid_point(bv0, bv1));
            let (_, d1u, d1v) = basis.d1(um, vm);
            // `cxx:1217-1219`: `GeomLProp_SLProps(TheSurf, Um, Vm, 1,
            // Precision::Confusion()).Normal()`. OCCT's `Normal()` throws
            // `LProp_NotDefined` when `CSLib::Normal` does not report
            // `CSLib_Done`; the port reports the same condition as an unhandled
            // surface and lets the face be dropped.
            let normal = surface_normal_from_d1(&d1u, &d1v, occt_core::precision::CONFUSION)?;
            // `cxx:1215`: the basis is transferred with the *caller's* parameter
            // range, not with its own bounds.
            let surface = self.emit_surface_entity(basis.as_ref(), u0, u1, v0, v1);
            return Some(self.emit_offset_surface(&normal, distance / IGES_UNIT, surface));
        }
        // Everything else has no branch in `TransferSurface` and returns a null
        // handle: the face is dropped, never replaced by a plane.
        None
    }

    /// Entity 140 (`GeomToIGES_GeomSurface::TransferSurface(Geom_OffsetSurface)`,
    /// `GeomToIGES_GeomSurface.cxx:1195-1237`, written by
    /// `IGESGeom_ToolOffsetSurface::WriteOwnParams`,
    /// `IGESGeom_ToolOffsetSurface.cxx:105-113`):
    ///
    /// `140, indicator.x, indicator.y, indicator.z, distance, surface;`
    ///
    /// - the indicator is the basis surface's unit normal at the mid point of the
    ///   offset's bounds, every component divided by the file unit
    ///   (`cxx:1217-1226`);
    /// - the distance is the signed `Offset()` divided by the file unit
    ///   (`cxx:1216`);
    /// - the surface is the recursively transferred basis (`cxx:1215`). When that
    ///   transfer yields a null handle OCCT still builds the 140 and
    ///   `IGESData_IGESWriter::Send` writes the reference as `0`, so the port
    ///   writes `0` too instead of dropping the entity.
    pub(super) fn emit_offset_surface(
        &mut self,
        indicator: &GpVec,
        distance: f64,
        surface: Option<usize>,
    ) -> usize {
        let params = format!(
            "140,{},{},{},{},",
            num(indicator.x() / IGES_UNIT),
            num(indicator.y() / IGES_UNIT),
            num(indicator.z() / IGES_UNIT),
            num(distance)
        );
        match surface {
            Some(idx) => self.emit_refs(140, 0, format!("{params}#0;"), &[idx]),
            None => self.emit_refs(140, 0, format!("{params}0;"), &[]),
        }
    }

}
