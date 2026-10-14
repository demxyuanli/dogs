//! `StepToTopoDS_TranslateFace`: `ADVANCED_FACE` resolution - surface, bounds, natural bound, vertex-loop detection. Source: `StepToTopoDS_TranslateFace.cxx`.
use super::*;

impl<'a> Resolver<'a> {
    pub(in crate::step) fn resolve_face(&self, rec: &'a Record) -> Result<TopoShape, String> {
        // Two argument layouts occur in the wild:
        //  * STEP-214 (ISO standard, written by OCCT/FreeCAD and by this port's
        //    writer since the cylinder-validation fix): surface is the third
        //    argument — `ADVANCED_FACE(name, bounds, surface, same_sense)`.
        //  * legacy files from this port's earlier writer: surface is the
        //    second argument — `ADVANCED_FACE('', #surface, (bounds), .T.)`.
        // Detect by checking which argument holds a surface reference.
        let surf_ref = if let Some(r) = parse_ref(&rec.args[2]) {
            r
        } else {
            parse_ref(&rec.args[1]).ok_or("ADVANCED_FACE: bad surface ref")?
        };
        let surface = self.resolve_surface(surf_ref)?;
        // StepToTopoDS_TranslateFace.cxx:557-567: a B-spline surface STEP face
        // is forced periodic when geometrically closed (ShapeCustom_Surface::
        // ConvertToPeriodic). Without it a pcurve running once around the closed
        // direction leaves the parameter window and is extrapolated.
        let is_bspline_surface_entity = self
            .record(surf_ref)
            .map(|r| {
                r.type_name == "B_SPLINE_SURFACE"
                    || r.type_name == "B_SPLINE_SURFACE_WITH_KNOTS"
            })
            .unwrap_or(false);
        let surface = if is_bspline_surface_entity {
            crate::shape_custom_surface::convert_to_periodic(&surface, self.precision)
                .unwrap_or(surface)
        } else {
            surface
        };
        // Bounds are the other of the two leading arguments.
        let bounds = if rec.args.get(2).map(|s| s.starts_with('(')).unwrap_or(false) {
            parse_ref_list(&rec.args[2])
        } else {
            parse_ref_list(&rec.args[1])
        };
        // `StepToTopoDS_TranslateFace.cxx:572-573` / `716-752`.
        let same_sense = parse_logical(rec.args.last().map(String::as_str), true);
        let a_same_sense = if self.step_surface_is_reversed(surf_ref) {
            !same_sense
        } else {
            same_sense
        };
        // `StepToTopoDS_TranslateFace.cxx:608-634` — a lone VertexLoop on a
        // sphere / BSpline / revolution is the whole closed surface; add
        // `BRepLib_MakeFace` natural bounds and skip the vertex loop itself.
        // `cxx:629` passes `Precision()`, the translator's `myPrecision`
        // (`STEPControl_ActorRead.cxx:2370-2384`) = file uncertainty x
        // `LengthFactor` (`step_precision`), not `Precision::Confusion()`.
        if bounds.len() == 1
            && self.bound_loop_is_vertex_loop(bounds[0])
            && (surface.gp_sphere().is_some()
                || surface.is_bspline_surface()
                || surface.is_surface_of_revolution())
        {
            let mut face = crate::brep_lib_make_face::make_face_from_surface(
                surface,
                self.precision,
            );
            face.0.set_orientation(if same_sense {
                Orientation::Forward
            } else {
                Orientation::Reversed
            });
            return Ok(face.0);
        }
        let mut wires = Vec::with_capacity(bounds.len());
        // `TranslateEdgeLoop.cxx:272`: `ForwardWire = FaceBound->Orientation()`
        // for SelectForwardSeam Compose — captured before sameSense reverse.
        let mut wire_bound_oris = Vec::with_capacity(bounds.len());
        for &b in &bounds {
            let mut s = self.resolve_shape(b)?;
            if !s.is_wire() {
                return Err(format!("#{b}: expected wire in face bounds"));
            }
            wire_bound_oris.push(s.orientation());
            // Bound orientation is already on the wire from `resolve_outer_bound`.
            // A reversed `Face_Surface` (or negative-major torus) flips it again
            // so the stored wire matches CAS.CADE (`cxx:714-723`).
            if !a_same_sense {
                s.reverse();
            }
            wires.push(Wire(s));
        }
        let mut face = self.b.make_face(surface, &wires);
        face.0.set_orientation(if same_sense {
            Orientation::Forward
        } else {
            Orientation::Reversed
        });
        let face_key = GeometryRegistry::shape_key(&face.0);
        // Associate each wire edge's SURFACE_CURVE pcurve with this face's
        // surface (`BRep_Builder::UpdateEdge(edge, pcurve, face, tol)`).
        let face_ori = face.0.orientation();
        for (w, &wire_o) in wires.iter_mut().zip(wire_bound_oris.iter()) {
            // EdgeO is ORIENTED_EDGE storage (`cxx:700`), not wire-composed.
            // Iterate EDGE_LOOP order so last seam UpdateEdge matches cxx.
            let loop_edges = edges_stored_on_wire(w);
            // `StepToTopoDS_GeometricTool.cxx:104-116`: the seam test counts how
            // many oriented edges of the loop reference the same step edge.
            let loop_refs: Vec<usize> = loop_edges
                .iter()
                .map(|x| {
                    self.edge_curve_ref
                        .borrow()
                        .get(&(Arc::as_ptr(&x.0.tshape) as usize))
                        .copied()
                        .unwrap_or(0)
                })
                .collect();
            for e in &loop_edges {
                let key = Arc::as_ptr(&e.0.tshape) as usize;
                let curve_ref = self.edge_curve_ref.borrow().get(&key).copied().unwrap_or(0);
                let nb_oe = loop_refs.iter().filter(|&&r| r == curve_ref).count();
                self.associate_edge_pcurve(e, surf_ref, face_key, face_ori, wire_o, nb_oe)?;
            }
            // `TranslateEdgeLoop.cxx:844-868` EdgeProjAux before CheckPCurves.
            // cxx:815 is only `B.Add(W,E)`; after CheckPCurves (cxx:875) the
            // function returns — no OrientEdgesOnWire / wire child reorder.
            // `TranslateEdgeLoop.cxx:236` `preci = Precision()`, set by
            // `STEPControl_ActorRead` from the file's
            // `GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT` uncertainty or from
            // `read.precision.val` (default 1.e-03); see `step_precision`.
            let preci = self.precision;
            // `StepToTopoDS_TranslateEdgeLoop.cxx:875` -> `CheckPCurves`
            // (`:105-177`) -> `XSAlgo_ShapeProcessor::CheckPCurve`
            // (`XSAlgo_ShapeProcessor.cxx:344-401`): the reader's "advanced check"
            // drops a pcurve whose U/V span wraps the surface or whose end
            // points disagree with the 3D curve beyond the read precision; the
            // dropped pcurve is re-projected by `ShapeFix_Edge::FixAddPCurve`
            // in the `ShapeFix` pass below (OCCT's `XSAlgo` pipeline).
            for e in &loop_edges {
                crate::shhealing::xsalgo_check_pcurve(e, &face, self.precision);
            }
            crate::shhealing::project_wire_pcurve_ranges(w, &face, preci);
            match self.current_heal_context() {
                // `ShapeFix_Face::Perform` runs its wire round on
                // `S = Context()->Apply(myFace)` (`ShapeFix_Face.cxx:379-380`),
                // i.e. after the shared `ReShape` has substituted every edge a
                // neighbouring face's `ShapeFix_ComposeShell` recorded. The
                // `FixAddPCurve` therefore has to land on the *substituted* edge
                // (`ShapeFix_ComposeShell.cxx:3451`), otherwise the pcurve is
                // attached to the pre-substitution edge that `Replace` drops and
                // the seam face later reads a pcurve-less wire.
                Some(mut ctx) => {
                    let applied = ctx.apply(&face.0);
                    let perf_face = if applied.is_face() {
                        Face(applied)
                    } else {
                        face.clone()
                    };
                    // `TopoDS_Iterator iter(S, false)` (`ShapeFix_Face.cxx:400`):
                    // the wire is handed to `ShapeFix_Wire::Load` with its stored
                    // orientation, not the one composed through the face.
                    for mut aw in wires_stored_on_face(&perf_face) {
                        crate::shhealing::check_pcurves_and_shift(
                            &mut aw, &perf_face, preci, false,
                        );
                    }
                }
                None => {
                    crate::shhealing::check_pcurves_and_shift(w, &face, preci, false);
                }
            }
        }
        // `BRepLib_MakeFace.cxx:860-866` forced SameParameter is already in
        // `make_face_uv` for natural-bound Offset faces. Wiring it on STEP
        // Offset faces (after EdgeProjAux) densifies Shape 6141 -> 6213 vs
        // occ 6150; TranslateFace.cxx has no equivalent call.
        //
        // `FromSTEP.FixShape` (`ShapeProcess_OperLibrary.cxx:830`) runs
        // `ShapeFix_Face::Perform`, whose `FixMissingSeam` step is
        // `ShapeFix_Face.cxx:492-498`. Measured this round on `data/occ/a3n00.stp`
        // by toggling the operator in the OCCT probe (`--faceids` counts
        // `TopExp_Explorer(face, TopAbs_WIRE)`): ShapeProcess ON gives
        // `{1:208, 2:9, 4:1, 6:4, 10:4}`, `--nofix` gives
        // `{1:163, 2:53, 4:2, 6:4, 10:4}` — the latter is the file's own
        // `FACE_BOUND` histogram and the port's, so the operator demonstrably
        // merges 45 two-bound faces into one seam-closed wire during import
        // (specs/_a3n00_gap_analysis.md §9.368.6). §9.219's link-map reading
        // ("`Perform` is not called on this read path") is superseded by that
        // measurement; see D27/§9.290.
        //
        // The seam construction and the `w2 != null` merge this needs are
        // `ShapeFix_Face.cxx:1899-2330`, ported in
        // `shhealing/shape_fix_face.rs::fix_missing_seam` (verified per face:
        // a3n00 f=192 turns 2 wires / 1 edge each into 1 wire / 5 edges, the same
        // structure the OCCT model has). `check_pcurves_and_shift` then repairs
        // the pcurves of the rebuilt wires (`cxx:365-480` first wire round).
        let mut sff = match self.current_heal_context() {
            Some(ctx) => crate::shhealing::ShapeFixFace::with_face_and_context(&face, ctx),
            None => crate::shhealing::ShapeFixFace::with_face(&face),
        };
        sff.result = Some(face.0.clone());
        // `ShapeFix_Face::Perform`'s tail (`ShapeFix_Face.cxx:482-498`):
        // `myResult = myFace`, then `FixPeriodicDegenerated` (`cxx:486-489`)
        // — the conic degenerate-apex edge — then `FixMissingSeam`
        // (`cxx:492-498`), which picks up the 2-wire cone face through
        // `Context()->Apply` (`cxx:1737-1741`).
        sff.fix_periodic_degenerated();
        let seam_ret = sff.fix_missing_seam();
        if seam_ret {
            if let Some(res) = sff.result.clone() {
                // `ShapeFix_Face::FixMissingSeam` returns `CompShell.Result()`
                // (`cxx:2266`), which is a face **or** a shell when
                // `MakeFacesOnPatch` splits the periodic face into several
                // patches (`ShapeFix_ComposeShell.cxx:2978-3271`). OCCT
                // stores it with `Context()->Replace(myFace, myResult)`
                // (`cxx:2268`) without a type test; `ShapeFix_Shape`
                // materialises it through `Context()->Apply(S)`
                // (`ShapeFix_Shape.cxx:257/294`), which re-adds a
                // replacement's children of the original type
                // (`ShapeBuild_ReShape.cxx:282-299`). The shell is expanded
                // by the caller (`resolve_shell`).
                let sub_faces = if res.shape_type() == crate::abs::ShapeType::Face {
                    vec![crate::shape::Face(res.clone())]
                } else if res.shape_type() == crate::abs::ShapeType::Shell {
                    crate::topo_tools_full::faces_of(&res)
                } else {
                    Vec::new()
                };
                if !sub_faces.is_empty() {
                    for rf in &sub_faces {
                        for mut w in wires_stored_on_face(rf) {
                            crate::shhealing::check_pcurves_and_shift(
                                &mut w,
                                rf,
                                self.precision,
                                false,
                            );
                        }
                    }
                    // `ShapeFix_Face::Perform`'s second wire round
                    // (`ShapeFix_Face.cxx:583-598`): `FixLoopWire` splits a
                    // boundary wire that revisits a vertex into one wire per
                    // loop. It runs on `myResult`'s wires, i.e. after the
                    // seam step, for every face (not only the healed ones).
                    for rf in &sub_faces {
                        crate::shhealing::split_loop_wires(rf);
                    }
                    // `ShapeFix_Face::Perform`'s second-part loop
                    // (`ShapeFix_Face.cxx:676-716`) runs on every face of
                    // `myResult`, substituting the rebuilt faces through the
                    // shared context (`cxx:1639-1641`, `cxx:917`, `cxx:1098`).
                    let mut healed: Vec<crate::shape::TopoShape> =
                        Vec::with_capacity(sub_faces.len());
                    for rf in &sub_faces {
                        healed.push(fix_face_second_part(&mut sff, rf));
                    }
                    let split_happened = healed
                        .iter()
                        .any(|s| s.shape_type() != crate::abs::ShapeType::Face);
                    let res = if split_happened {
                        // `FixSplitFace` replaced a face by a compound
                        // (`ShapeFix_Face.cxx:2988-2999`).
                        self.b.make_compound_of(&healed).0
                    } else if res.shape_type() == crate::abs::ShapeType::Shell {
                        let faces: Vec<crate::shape::Face> = healed
                            .iter()
                            .map(|s| crate::shape::Face(s.clone()))
                            .collect();
                        self.b.make_shell(&faces).0
                    } else {
                        healed.remove(0)
                    };
                    return Ok(res);
                }
            }
        }
        crate::shhealing::split_loop_wires(&face);
        // `ShapeFix_Face::Perform`'s second-part loop (`ShapeFix_Face.cxx:676-716`):
        // the wire steps (`cxx:676-685`), `FixOrientation` (`cxx:687-697`),
        // `FixAddNaturalBound` (`cxx:702-708`) and `FixSplitFace` (`cxx:710-716`).
        let res = fix_face_second_part(&mut sff, &face);
        Ok(res)
    }

    /// True when `bound_id` is a `FACE_BOUND` / `FACE_OUTER_BOUND` whose loop
    /// is a `VERTEX_LOOP` (`StepToTopoDS_TranslateFace.cxx:608`).
    pub(super) fn bound_loop_is_vertex_loop(&self, bound_id: usize) -> bool {
        let Ok(rec) = self.record(bound_id) else {
            return false;
        };
        match rec.type_name.as_str() {
            "FACE_BOUND" | "FACE_OUTER_BOUND" => parse_ref(&rec.args[1])
                .and_then(|id| self.record(id).ok())
                .map(|r| r.type_name == "VERTEX_LOOP")
                .unwrap_or(false),
            "VERTEX_LOOP" => true,
            _ => false,
        }
    }

}
