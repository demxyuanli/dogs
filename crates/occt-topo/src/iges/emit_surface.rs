//! Surface, revolution and transformation-matrix entity emitters. Source: `GeomToIGES_GeomSurface` with `IGESConvGeom_GeomBuilder` and `IGESGeom_ToolTransformationMatrix`.
use super::prelude::*;
use super::*;

impl IgesWriter {

    /// Entity 123 Direction (`IGESGeom_ToolDirection::WriteOwnParams`,
    /// `IGESGeom_ToolDirection.cxx:64-74`): the three components of the unit
    /// vector.
    pub(super) fn emit_direction(&mut self, d: &GpDir) -> usize {
        self.emit(
            123,
            0,
            format!("123,{},{},{};", num(d.x()), num(d.y()), num(d.z())),
        )
    }

    /// Entity 196 (`GeomToIGES_GeomSurface::TransferSphericalSurface`,
    /// `GeomToIGES_GeomSurface.cxx:1366-1400`, written by
    /// `IGESSolid_ToolSphericalSurface::WriteOwnParams`,
    /// `IGESSolid_ToolSphericalSurface.cxx:70-80`):
    /// `196, center_point, radius, axis_direction, reference_direction;` — the
    /// parametrised form, which is the one OCCT's transfer produces.
    ///
    /// **Not wired**: this is the BRep-mode arm (`myBRepMode && myAnalytic`,
    /// `cxx:579-586`). OCCT's default `write.iges.brep.mode = 0`
    /// (`IGESData.cxx:90-94`) keeps the writer in Faces mode, where a sphere
    /// becomes a 120 revolution surface instead (see
    /// [`IgesWriter::emit_local_revolution_surface`]).
    pub(super) fn emit_spherical_surface(
        &mut self,
        center: &GpPnt,
        radius: f64,
        axis: &GpDir,
        x_dir: &GpDir,
    ) -> usize {
        let c = self.emit_point(center);
        let a = self.emit_direction(axis);
        let r = self.emit_direction(x_dir);
        self.emit_refs(
            196,
            0,
            format!("196,#0,{},#1,#2;", num(radius)),
            &[c, a, r],
        )
    }

    /// Entity 192 (`GeomToIGES_GeomSurface::TransferCylindricalSurface`,
    /// `GeomToIGES_GeomSurface.cxx:1280-1314`, written by
    /// `IGESSolid_ToolCylindricalSurface::WriteOwnParams`):
    /// `192, location_point, axis_direction, radius, reference_direction;`.
    ///
    /// **Not wired**: BRep-mode only (`myBRepMode && myAnalytic`, `cxx:555-562`);
    /// in Faces mode a cylinder is a 120.
    pub(super) fn emit_cylindrical_surface(
        &mut self,
        location: &GpPnt,
        axis: &GpDir,
        radius: f64,
        x_dir: &GpDir,
    ) -> usize {
        let l = self.emit_point(location);
        let a = self.emit_direction(axis);
        let r = self.emit_direction(x_dir);
        self.emit_refs(
            192,
            0,
            format!("192,#0,#1,{},#2;", num(radius)),
            &[l, a, r],
        )
    }

    /// Entity 194 (`TransferConicalSurface`, `cxx:1318-1362`, written by
    /// `IGESSolid_ToolConicalSurface::WriteOwnParams`):
    /// `194, location_point, axis_direction, ref_radius, semi_angle_deg,
    /// reference_direction;`. A negative semi-angle is written by mirroring the
    /// reference point through the apex, negating the angle and reversing the
    /// reference direction (`cxx:1344-1350`).
    ///
    /// **Not wired**: BRep-mode only (`myBRepMode && myAnalytic`, `cxx:567-574`);
    /// in Faces mode a cone is a 120.
    pub(super) fn emit_conical_surface(
        &mut self,
        location: &GpPnt,
        apex: &GpPnt,
        axis: &GpDir,
        ref_radius: f64,
        semi_angle: f64,
        x_dir: &GpDir,
    ) -> usize {
        let (loc, angle, xd) = if semi_angle < 0.0 {
            (
                GpPnt::new(
                    2.0 * apex.x() - location.x(),
                    2.0 * apex.y() - location.y(),
                    2.0 * apex.z() - location.z(),
                ),
                -semi_angle,
                GpDir::new(-x_dir.x(), -x_dir.y(), -x_dir.z())
                    .unwrap_or_else(|_| *x_dir),
            )
        } else {
            (*location, semi_angle, *x_dir)
        };
        let l = self.emit_point(&loc);
        let a = self.emit_direction(axis);
        let r = self.emit_direction(&xd);
        self.emit_refs(
            194,
            0,
            format!(
                "194,#0,#1,{},{},#2;",
                num(ref_radius),
                num(angle * 180.0 / std::f64::consts::PI)
            ),
            &[l, a, r],
        )
    }

    /// Entity 198 (`TransferToroidalSurface`, `cxx:1402-1435`, written by
    /// `IGESSolid_ToolToroidalSurface::WriteOwnParams`):
    /// `198, center_point, axis_direction, major_radius, minor_radius,
    /// reference_direction;`.
    ///
    /// **Not wired**: BRep-mode only (`myBRepMode && myAnalytic`, `cxx:591-598`);
    /// in Faces mode a torus is a 120.
    pub(super) fn emit_toroidal_surface(
        &mut self,
        center: &GpPnt,
        axis: &GpDir,
        major: f64,
        minor: f64,
        x_dir: &GpDir,
    ) -> usize {
        let c = self.emit_point(center);
        let a = self.emit_direction(axis);
        let r = self.emit_direction(x_dir);
        self.emit_refs(
            198,
            0,
            format!("198,#0,#1,{},{},#2;", num(major), num(minor)),
            &[c, a, r],
        )
    }

    /// Entity 128 (`GeomToIGES_GeomSurface::TransferBSplineSurface`,
    /// `GeomToIGES_GeomSurface.cxx:191-460`, written by
    /// `IGESGeom_ToolBSplineSurface::WriteOwnParams`,
    /// `IGESGeom_ToolBSplineSurface.cxx:64-125`):
    ///
    /// `128, indU, indV, degU, degV, closedU, closedV, polynomial, periodicU,
    ///   periodicV, knotU[-degU .. indU+1], knotV[-degV .. indV+1],
    ///   weights[0..indU][0..indV], poles[0..indU][0..indV],
    ///   UMin, UMax, VMin, VMax;`
    ///
    /// with `indU = nb_poles_u - 1`, `indV = nb_poles_v - 1` and the flattened
    /// knot vectors. A Bezier surface becomes the equivalent single-span
    /// B-spline (end knots with multiplicity degree+1) as
    /// `GeomConvert::SurfaceToBSplineSurface` would produce.
    ///
    /// `closedU/V` come from `IsUClosed`/`IsVClosed` (`cxx:347-348`, see
    /// [`iso_rows_equal`]) of the **unperiodized** surface, `polynomial` from
    /// `Polynom = !(RationU || RationV)` (`cxx:355`; the port's B-spline surface
    /// carries one rational flag for both directions), and `periodicU/V` from the
    /// original `IsUPeriodic`/`IsVPeriodic` (`cxx:235-236`, `:448-449`). The
    /// bounds fix implements both arms (`cxx:244-302`) and a periodic direction is
    /// unperiodized before its knots and poles are read (`cxx:303-343`,
    /// `SetUNotPeriodic`/`SetVNotPeriodic` → `BSplSLib::Unperiodize`,
    /// [`occt_core::bspl::unperiodize_direction`]). The bounds fix and the
    /// re-origin implement both arms of `cxx:244-343`, including the
    /// `SetUOrigin`/`SetVOrigin` call (`cxx:310-320`, `:330-340`, issue 26138)
    /// that OCCT performs when the written range's period shift differs from the
    /// range start's.
    pub(super) fn emit_bspline_surface(
        &mut self,
        surf: &dyn Surface,
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
    ) -> Option<usize> {
        let (mut poles, mut knots_u, mut knots_v, deg_u, deg_v) = match (
            surf.bspline_surface_poles(),
            surf.bspline_surface_uknots(),
            surf.bspline_surface_vknots(),
        ) {
            (Some(p), Some(ku), Some(kv)) => (
                p.to_vec(),
                ku.to_vec(),
                kv.to_vec(),
                surf.u_degree().max(1) as usize,
                surf.v_degree().max(1) as usize,
            ),
            _ => {
                if std::env::var("IGES_TRACE").is_ok() {
                    eprintln!("TRACE bspline missing data: is_bspline={}", surf.is_bspline_surface());
                }
                return None;
            }
        };
        if poles.is_empty() || poles[0].is_empty() {
            return None;
        }
        let (nu0, nv0) = (poles.len(), poles[0].len());

        // `TransferSurface(Geom_RectangularTrimmedSurface)` recurses on
        // `BasisSurface()` before this branch (`cxx:492-515`), so every quantity
        // read below - poles, knots, degree, periodicity and `Bounds` - belongs to
        // the untrimmed basis (`st` there plays the role of `start` here).
        let basis = surf.rectangular_trimmed_basis();
        let bs: &dyn Surface = match basis.as_deref() {
            Some(b) => b,
            None => surf,
        };

        let stored_weights = bs.bspline_surface_weights();
        let rational = stored_weights.is_some();
        let mut weights: Vec<Vec<f64>> = match stored_weights {
            Some(w) if w.len() == nu0 && w.iter().all(|r| r.len() == nv0) => w.to_vec(),
            _ => vec![vec![1.0; nv0]; nu0],
        };
        let polynomial = !rational;

        // `GeomToIGES_GeomSurface.cxx:235-236`: `PeriodicU/V` are read from the
        // **original** surface and written verbatim (`:448-449`).
        let period_u = bs.is_u_periodic();
        let period_v = bs.is_v_periodic();

        // An **open** direction's flat knot vector and pole grid must describe
        // the same representation: `nb_poles + degree + 1`. A periodic direction
        // is *not* trimmed to the open length -
        // `Geom_BSplineSurface::SetUPeriodic` (`Geom_BSplineSurface_1.cxx:940-981`)
        // trims the poles to `BSplCLib::NbPoles(degree, true, mults)` while
        // `updateUKnots` rebuilds the flat sequence with the periodic extension
        // (`BSplCLib::KnotSequenceLength`, `BSplCLib.cxx:455-474`) - and the base
        // multiplicities that would let the length be predicted are no longer
        // recoverable from the extended sequence (the leading `degree + 1 -
        // mults(1)` entries are the wrapped tail, shifted by one period). Such a
        // direction is therefore checked after the unperiodization below, on the
        // open representation OCCT's writer actually emits.
        let open_flat_ok = |flat: &[f64], nb: usize, deg: usize, periodic: bool| -> bool {
            periodic || flat.len() == nb + deg + 1
        };
        if !open_flat_ok(&knots_u, nu0, deg_u, period_u)
            || !open_flat_ok(&knots_v, nv0, deg_v, period_v)
        {
            if std::env::var("IGES_TRACE").is_ok() {
                eprintln!(
                    "TRACE bspline layout mismatch: nu={nu0} nv={nv0} deg=({deg_u},{deg_v}) ku={} kv={} | nbuv=({},{}) uvdeg=({},{}) pu={period_u} pv={period_v} trim={} bs={}",
                    knots_u.len(),
                    knots_v.len(),
                    surf.nb_u_poles(),
                    surf.nb_v_poles(),
                    surf.u_degree(),
                    surf.v_degree(),
                    surf.rectangular_trimmed_basis().is_some(),
                    surf.is_bspline_surface(),
                );
            }
            return None;
        }

        // `cxx:239-301`: the written range is clamped to the surface's own
        // bounds; the periodic arm snaps an end that already sits on a bound and
        // otherwise shifts the range into the period
        // (`ShapeAnalysis::AdjustToPeriod`), truncating it to one period. This
        // runs **before** the re-origin below, which needs `uShift`/`vShift`.
        let (su0, su1) = bs.u_range();
        let (sv0, sv1) = bs.v_range();
        let (u1_arg, v1_arg) = (u1, v1);
        let (mut u0, mut u1) = (u0, u1);
        let mut u_shift = 0.0;
        if period_u {
            if (u0 - su0).abs() < occt_core::precision::PCONFUSION {
                u0 = su0;
            }
            if (u1 - su1).abs() < occt_core::precision::PCONFUSION {
                u1 = su1;
            }
            u_shift = crate::pcurve_full::adjust_to_period(u0, su0, su1);
            u0 += u_shift;
            u1 += u_shift;
            if u1 - u0 > su1 - su0 {
                u1 = u0 + (su1 - su0);
            }
        } else {
            u0 = u0.max(su0);
            u1 = u1.min(su1);
        }
        let (mut v0, mut v1) = (v0, v1);
        let mut v_shift = 0.0;
        if period_v {
            if (v0 - sv0).abs() < occt_core::precision::PCONFUSION {
                v0 = sv0;
            }
            if (v1 - sv1).abs() < occt_core::precision::PCONFUSION {
                v1 = sv1;
            }
            v_shift = crate::pcurve_full::adjust_to_period(v0, sv0, sv1);
            v0 += v_shift;
            v1 += v_shift;
            if v1 - v0 > sv1 - sv0 {
                v1 = v0 + (sv1 - sv0);
            }
        } else {
            v0 = v0.max(sv0);
            v1 = v1.min(sv1);
        }

        // `cxx:303-343`: `SetUOrigin`/`SetVOrigin` re-origin a periodic
        // B-spline direction when the written range's own period shift differs
        // from the range's start one (issue 26138: synchronize the p-curve
        // ranges with the surface bounds). `LocateU(Vmin)` gives the window
        // index the new origin starts at; the re-origin only permutes the
        // periodic representation, so the geometry is unchanged.
        let u_reorigin = period_u
            && (u_shift - crate::pcurve_full::adjust_to_period(u1_arg, su0, su1)).abs()
                > occt_core::precision::PCONFUSION;
        let v_reorigin = period_v
            && (v_shift - crate::pcurve_full::adjust_to_period(v1_arg, sv0, sv1)).abs()
                > occt_core::precision::PCONFUSION;
        if u_reorigin || v_reorigin {
            let mut ms = occt_geom::bspline_surface::GeomBSplineSurface {
                poles: poles.clone(),
                knots_u: knots_u.clone(),
                knots_v: knots_v.clone(),
                deg_u,
                deg_v,
                weights: if rational { Some(weights.clone()) } else { None },
                u_periodic: period_u,
                v_periodic: period_v,
                res_latch: occt_geom::bspline_surface::ResolutionLatch::new(),
            };
            let mut changed = false;
            if u_reorigin {
                let (left, _right) = ms.locate_u(u0, occt_core::precision::PCONFUSION);
                if ms.set_u_origin(left).is_ok() {
                    changed = true;
                }
            }
            if v_reorigin {
                let (left, _right) = ms.locate_v(v0, occt_core::precision::PCONFUSION);
                if ms.set_v_origin(left).is_ok() {
                    changed = true;
                }
            }
            if changed {
                poles = ms.poles;
                if let Some(w) = ms.weights {
                    weights = w;
                }
                knots_u = ms.knots_u;
                knots_v = ms.knots_v;
            }
        }

        // `cxx:303-343`: a periodic B-spline surface is **unperiodized** before
        // its knots and poles are read (`SetUNotPeriodic` / `SetVNotPeriodic` →
        // `BSplSLib::Unperiodize`, `Geom_BSplineSurface_1.cxx:1238-1300`), so
        // the written 128 has the open knot vector IGES expects. The pole
        // extension is the cyclic wrap of `BSplCLib::Unperiodize`
        // (`BSplCLib.cxx:3076-3079`, modulus `Poles.Length()` =
        // `BSplCLib::NbPoles(degree, true, mults)`); weights follow the same wrap
        // because `BSplSLib::Unperiodize` works on homogeneous poles.
        if period_u {
            let (nf, map) =
                occt_core::bspl::unperiodize_direction(deg_u as i32, &knots_u, nu0);
            knots_u = nf;
            poles = map.iter().map(|k| poles[*k].clone()).collect();
            weights = map.iter().map(|k| weights[*k].clone()).collect();
        }
        if period_v {
            let (nf, map) =
                occt_core::bspl::unperiodize_direction(deg_v as i32, &knots_v, nv0);
            knots_v = nf;
            for row in poles.iter_mut() {
                *row = map.iter().map(|k| row[*k].clone()).collect();
            }
            for row in weights.iter_mut() {
                *row = map.iter().map(|k| row[*k]).collect();
            }
        }

        let (nu, nv) = (poles.len(), poles[0].len());
        if knots_u.len() != nu + deg_u + 1 || knots_v.len() != nv + deg_v + 1 {
            if std::env::var("IGES_TRACE").is_ok() {
                eprintln!(
                    "TRACE bspline post-unperiodize mismatch: nu={nu} nv={nv} deg=({deg_u},{deg_v}) ku={} kv={} pu={period_u} pv={period_v}",
                    knots_u.len(),
                    knots_v.len()
                );
            }
            return None;
        }
        let (ind_u, ind_v) = (nu - 1, nv - 1);

        // `Geom_BSplineSurface_1.cxx:1026-1130` `SetUOrigin`/`SetVOrigin` are
        // applied **before** the unperiodization when the written range's own
        // period shift differs from the range start's (`cxx:303-343`), see the
        // `u_reorigin`/`v_reorigin` block above; the flags below are read from
        // the resulting surface, exactly like OCCT reads them after its own
        // re-origin.
        //
        // `CloseU = mysurface->IsUClosed()` / `CloseV` (`cxx:347-348`), read from
        // the **unperiodized** surface: a periodic direction whose first and last
        // pole rows coincide is closed, which the pole comparison below detects
        // after the cyclic extension.
        let closed_u = iso_rows_equal(&poles[0], &poles[nu - 1])
            && (!rational || weights[0].iter().zip(&weights[nu - 1]).all(|(a, b)| weight_equal(*a, *b)));
        let closed_v = (0..nu).all(|i| {
            iso_rows_equal(&poles[i][0..1], &poles[i][nv - 1..nv])
                && (!rational || weight_equal(weights[i][0], weights[i][nv - 1]))
        });

        let mut s = format!(
            "128,{ind_u},{ind_v},{deg_u},{deg_v},{},{},{},{},{}",
            u8::from(closed_u),
            u8::from(closed_v),
            u8::from(polynomial),
            u8::from(period_u),
            u8::from(period_v)
        );
        for k in &knots_u {
            s.push(',');
            s.push_str(&num(*k));
        }
        for k in &knots_v {
            s.push(',');
            s.push_str(&num(*k));
        }
        for j in 0..nv {
            for i in 0..nu {
                s.push(',');
                s.push_str(&num(weights[i][j]));
            }
        }
        for j in 0..nv {
            for i in 0..nu {
                let p = &poles[i][j];
                s.push_str(&format!(",{},{},{}", num(p.x()), num(p.y()), num(p.z())));
            }
        }
        s.push_str(&format!(
            ",{},{},{},{};",
            num(u0),
            num(u1),
            num(v0),
            num(v1)
        ));
        Some(self.emit(128, 0, s))
    }

    /// Entity 120 (`GeomToIGES_GeomSurface::TransferSurface(
    /// Geom_SurfaceOfRevolution)`, `GeomToIGES_GeomSurface.cxx:1111-1188`, written
    /// by `IGESGeom_ToolSurfaceOfRevolution::WriteOwnParams`, `:119-127`):
    ///
    /// `120, axis_line, generatrix, start_angle, end_angle;`
    ///
    /// - the generatrix is `GC.TransferCurve(BasisCurve, V1, V2)` (`cxx:1155`);
    /// - the axis is an `IGESGeom_Line` (`cxx:1163-1183`) whose start point is
    ///   the axis location and whose second point is `Location - Direction` -
    ///   the CAS.CADE axis reversed (`#30 rln`, `#36 BUC60328 face 7`), the same
    ///   convention the reader inverts (`IGESToBRep_TopoSurface.cxx:768`);
    /// - the angles come from `Surf->Init(Axis, Generatrix, 2*M_PI - U2,
    ///   2*M_PI - U1)` (`cxx:1185`), the CAS.CADE/IGES phase convention.
    pub(super) fn emit_surface_of_revolution(
        &mut self,
        surf: &dyn Surface,
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
    ) -> Option<usize> {
        let basis = surf.revolution_basis_curve()?;
        let axis = surf.revolution_axis()?;
        // `GC.TransferCurve(Curve, V1, V2)` (`cxx:1155`); a generatrix type with
        // no IGES writer in `TransferCurve` keeps the port's chord-line stand-in,
        // the same fallback `emit_edge_curve` uses.
        let generatrix = match self.emit_curve_range(basis.as_ref(), v0, v1) {
            Some(i) => i,
            None => {
                let p0 = basis.d0(v0);
                let p1 = basis.d0(v1);
                self.emit_line(&p0, &p1)
            }
        };
        let loc = axis.location();
        let d = axis.direction();
        let axis_line = self.emit_line(
            &loc,
            &GpPnt::new(loc.x() - d.x(), loc.y() - d.y(), loc.z() - d.z()),
        );
        let tau = 2.0 * std::f64::consts::PI;
        Some(self.emit_refs(
            120,
            0,
            format!(
                "120,#0,#1,{},{};",
                num(tau - u1),
                num(tau - u0)
            ),
            &[axis_line, generatrix],
        ))
    }

    /// `GeomToIGES_GeomSurface::TransferSurface(Geom_CylindricalSurface /
    /// Geom_ConicalSurface / Geom_SphericalSurface / Geom_ToroidalSurface)`
    /// (`GeomToIGES_GeomSurface.cxx:696-768`, `:777-852`, `:859-925`,
    /// `:933-993`).
    ///
    /// These four overloads are reached only from the `myBRepMode && myAnalytic`
    /// arm of the `Geom_ElementarySurface` branch (`cxx:552-598`). Under OCCT's
    /// default `write.iges.brep.mode = 0` (`IGESData.cxx:90-94`) the writer runs
    /// in Faces mode, so they are **not** taken and the elementary surfaces are
    /// written through the plain `TransferSurface(Geom_Surface)` overload as
    /// `IGESGeom_SurfaceOfRevolution` (120). This helper is that Faces-mode 120:
    /// the generatrix and the axis are built in the surface's **local** frame and
    /// the frame itself rides on a 124.
    ///
    /// - the axis is always the local Z axis, written as an `IGESGeom_Line` from
    ///   `(0,0,1)` to `(0,0,0)` (`cxx:744`, `:830`, `:902`, `:969`; the
    ///   `#30 rln 19.10.98` "IGES axis = reversed CAS.CADE axis");
    /// - `generatrix` is `GC.TransferCurve(<local generatrix>, V1, V2)` already
    ///   written by the caller (`cxx:819`, `:820`, `:896`, `:964`);
    /// - the angles are `Init(Axis, Generatrix, 2*PI - U2, 2*PI - U1)`
    ///   (`cxx:831`, `:907`, `:971`);
    /// - `IGESConvGeom_GeomBuilder::SetPosition(<surface>.Position())` plus
    ///   `IsIdentity()` decide whether the 124 is written
    ///   (`cxx:841-849`, `:915-921`, `:981-988`;
    ///   `IGESConvGeom_GeomBuilder.cxx:137-198`).
    pub(super) fn emit_local_revolution_surface(
        &mut self,
        generatrix: usize,
        frame: &GpAx3,
        u0: f64,
        u1: f64,
    ) -> usize {
        let axis = self.emit_line(&GpPnt::new(0.0, 0.0, 1.0), &GpPnt::new(0.0, 0.0, 0.0));
        let tau = 2.0 * std::f64::consts::PI;
        let de = self.emit_refs(
            120,
            0,
            format!("120,#0,#1,{},{};", num(tau - u1), num(tau - u0)),
            &[axis, generatrix],
        );
        if !frame_is_identity(frame) {
            let trsf = self.emit_transformation_matrix(frame, 1.0);
            self.set_trsf(de, trsf);
        }
        de
    }

    /// Entity 122 (`GeomToIGES_GeomSurface::TransferSurface(Geom_SurfaceOfLinearExtrusion)`,
    /// `GeomToIGES_GeomSurface.cxx:1032-1104`, written by
    /// `IGESGeom_ToolTabulatedCylinder::WriteOwnParams`, `:90-98`):
    ///
    /// `122, directrix, end_point.x, end_point.y, end_point.z;`
    ///
    /// - the U range is re-read from the surface's own `Bounds`
    ///   (`cxx:1067-1071`, the `OCC9490` fix), not from the caller's range;
    /// - the V range keeps OCCT's infinite handling (`cxx:1058-1065`);
    /// - the directrix is the basis curve translated so that its origin coincides
    ///   with the directrix origin (`cxx:1075-1096`): when
    ///   `|V1| > Precision::Confusion()` it is shifted by
    ///   `Value(U1,V1) - Value(U1,0)`, i.e. by `V1 * Direction`;
    /// - the terminate point is `start->Value(U1, V2)` (`cxx:1078`).
    pub(super) fn emit_tabulated_cylinder(&mut self, surf: &dyn Surface, v0: f64, v1: f64) -> Option<usize> {
        let basis = surf.extrusion_basis_curve()?;
        let (u1, u2) = surf.u_range();
        let vv1 = if occt_core::precision::Precision::is_negative_infinite(v0) {
            -occt_core::precision::INFINITE
        } else {
            v0
        };
        let vv2 = if occt_core::precision::Precision::is_positive_infinite(v1) {
            occt_core::precision::INFINITE
        } else {
            v1
        };
        let end = surf.value(u1, vv2);
        let copy = if vv1.abs() > occt_core::precision::CONFUSION {
            let shift = GpVec::from_pnts(&surf.value(u1, 0.0), &surf.value(u1, vv1));
            basis.translated(&shift)
        } else {
            basis.clone_dyn()
        };
        let directrix = match self.emit_curve_range(copy.as_ref(), u1, u2) {
            Some(i) => i,
            None => {
                let p0 = copy.d0(u1);
                let p1 = copy.d0(u2);
                self.emit_line(&p0, &p1)
            }
        };
        Some(self.emit_refs(
            122,
            0,
            format!(
                "122,#0,{},{},{};",
                num(end.x()),
                num(end.y()),
                num(end.z())
            ),
            &[directrix],
        ))
    }

    /// Entity 124 (`IGESConvGeom_GeomBuilder::MakeTransformation`,
    /// `IGESConvGeom_GeomBuilder.cxx:218-237`, written by
    /// `IGESGeom_ToolTransformationMatrix::WriteOwnParams`,
    /// `IGESGeom_ToolTransformationMatrix.cxx:90-104`): the frame's 3x4 matrix
    /// `R11 R12 R13 T1 R21 R22 R23 T2 R31 R32 R33 T3`, the translation column
    /// divided by the file unit; form 1 when the frame is left-handed
    /// (`rs->SetFormNumber(1)` when `thepos.IsNegative()`).
    pub(super) fn emit_transformation_matrix(&mut self, frame: &GpAx3, unit: f64) -> usize {
        // `IGESConvGeom_GeomBuilder::SetPosition(const gp_Ax3&)`
        // (`GeomBuilder.cxx:137-143`) runs `ps.SetTransformation(pos, gp::XOY())`,
        // then `MakeTransformation` (`GeomBuilder.cxx:218-237`) reads
        // `thepos.Value(i, j)` and divides column 4 by the file unit. Going
        // through `gp_Trsf::SetTransformation` (`gp_Trsf.cxx:172-194`) matters:
        // the trailing `matrix.Multiply(MA1)` (identity times the frame's basis
        // columns) normalises the signed zeros of `Value(i, j)`, and `IsNegative`
        // (`gp_Trsf::IsNegative`, `gp_Trsf.cxx`) is the form-number test.
        let mut ps = GpTrsf::identity();
        ps.set_transformation_from_to(frame, &xoy_frame());
        let form = if ps.is_negative() { 1 } else { 0 };
        let mut s = String::from("124");
        for i in 1..=3 {
            for j in 1..=3 {
                s.push(',');
                s.push_str(&num(ps.value(i, j)));
            }
            s.push(',');
            s.push_str(&num(ps.value(i, 4) / unit));
        }
        s.push(';');
        self.emit(124, form, s)
    }

}
