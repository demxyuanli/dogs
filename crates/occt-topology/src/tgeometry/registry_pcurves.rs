use super::*;

impl GeometryRegistry {
    /// The pcurve of edge `s` on the face identified by `face_key` (see
    /// `shape_key`), if one has been attached. For a seam edge this is the
    /// forward pcurve; use [`GeometryRegistry::edge_pcurves`] to get both.
    ///
    /// When no entry exists for `face_key`, a pcurve attached for another face
    /// that shares the same surface handle is returned (`BRep_Tool::CurveOnSurface`
    /// keys representations by `Geom_Surface`, not by face TShape).
    pub fn edge_pcurve(&self, s: &TopoShape, face_key: usize) -> Option<Arc<dyn Curve2d>> {
        self.edge_pcurves(s, face_key).into_iter().next()
    }

    /// All pcurves of edge `s` on the face identified by `face_key`
    /// (forward-then-reversed for a seam edge, one for a normal edge).
    pub fn edge_pcurves(&self, s: &TopoShape, face_key: usize) -> Vec<Arc<dyn Curve2d>> {
        let key = self.repr_key(face_key);
        let ts = s.tshape.read().unwrap();
        ts.edge_pcurves()
            .map(|p| p.pcurves(key))
            .unwrap_or_default()
    }

    /// `BRep_Tool::CurveOnSurface(E, F)`: match both the surface handle and
    /// `F.Location().Predivided(E.Location())`.
    pub fn edge_pcurves_on_face(
        &self,
        edge: &TopoShape,
        face: &TopoShape,
    ) -> Vec<Arc<dyn Curve2d>> {
        let key = self.repr_key(Self::shape_key(face));
        let wanted = face.location().predivided(edge.location());
        let ts = edge.tshape.read().unwrap();
        let Some(pcurves) = ts.edge_pcurves() else {
            return Vec::new();
        };
        if !pcurves.location(key).is_equal(&wanted) {
            return Vec::new();
        }
        pcurves.pcurves(key)
    }

    /// Attach a pcurve to edge `s` for the face identified by `face_key`.
    /// Mirrors `BRep_Builder::UpdateEdge(edge, curve2d, face, tol)`.
    pub fn set_edge_pcurve(&self, s: &TopoShape, face_key: usize, curve: Arc<dyn Curve2d>) {
        if crate::debug_flags::pc_set() {
            let info = curve
                .bspline_copy2d()
                .map(|b| format!("bs np={} nk={}", b.xs.len(), b.knots.len()))
                .unwrap_or_else(|| "nobs".to_string());
            eprintln!(
                "PCSET1 edge={:p} face={face_key} [{:.17e},{:.17e}] {info}",
                Arc::as_ptr(&s.tshape),
                curve.first_parameter(),
                curve.last_parameter()
            );
        }
        let (first, last) = self.edge_parameters(s);
        let key = self.repr_key(face_key);
        s.tshape
            .write()
            .unwrap()
            .edge_pcurves_mut()
            .set_pcurve(key, curve, first, last);
    }

    /// Replace the pcurves of edge `s` on the face identified by `face_key`.
    /// Mirrors the seam overload `BRep_Builder::UpdateEdge(edge, c1, c2, face)`.
    pub fn set_edge_pcurves(&self, s: &TopoShape, face_key: usize, curves: Vec<Arc<dyn Curve2d>>) {
        if crate::debug_flags::pc_set() {
            for c in &curves {
                let info = c
                    .bspline_copy2d()
                    .map(|b| format!("bs np={} nk={}", b.xs.len(), b.knots.len()))
                    .unwrap_or_else(|| "nobs".to_string());
                eprintln!(
                    "PCSET edge={:p} face={face_key} [{:.17e},{:.17e}] {info}",
                    Arc::as_ptr(&s.tshape),
                    c.first_parameter(),
                    c.last_parameter()
                );
            }
        }
        if crate::debug_flags::pc_s2d() && curves.len() >= 2 {
            let (e0, e1) = self.edge_parameters(s);
            let mut line = format!(
                "PCS2 edge={:p} face={face_key} n={}",
                Arc::as_ptr(&s.tshape),
                curves.len()
            );
            for (q, c) in curves.iter().enumerate() {
                let (mut a, mut b) = (c.first_parameter(), c.last_parameter());
                if !(a.is_finite() && b.is_finite()) {
                    a = e0;
                    b = e1;
                }
                let p0 = c.d0(a);
                let p1 = c.d0(b);
                line += &format!(
                    " q{q}=({:.6},{:.6})->({:.6},{:.6})",
                    p0.x(),
                    p0.y(),
                    p1.x(),
                    p1.y()
                );
            }
            let bt = std::backtrace::Backtrace::force_capture();
            let text = format!("{bt}");
            let mut names: Vec<String> = Vec::new();
            for l in text.lines() {
                if let Some(pos) = l.find("occt_topo::") {
                    names.push(l[pos..].trim().to_string());
                    if names.len() >= 8 {
                        break;
                    }
                }
            }
            eprintln!("{line} {}", names.join(" <- "));
        }
        let (first, last) = self.edge_parameters(s);
        let key = self.repr_key(face_key);
        s.tshape
            .write()
            .unwrap()
            .edge_pcurves_mut()
            .set_pcurves(key, curves, first, last);
    }

    /// `BRep_Builder::UpdateEdge(E, C1, C2, F, Tol)`: store the
    /// CurveOnSurface representation location relative to the edge.
    pub fn set_edge_pcurves_on_face(
        &self,
        edge: &TopoShape,
        face: &TopoShape,
        curves: Vec<Arc<dyn Curve2d>>,
    ) {
        let (first, last) = self.edge_parameters(edge);
        let key = self.repr_key(Self::shape_key(face));
        let location = face.location().predivided(edge.location());
        edge.tshape
            .write()
            .unwrap()
            .edge_pcurves_mut()
            .set_pcurves_at(key, location, curves, first, last);
    }

    /// `ShapeBuild_Edge::CopyPCurves` (`ShapeBuild_Edge.cxx:360-413`) at the
    /// `BRep_TEdge` level: copy every CurveOnSurface representation (pcurve
    /// list + range) of `from` onto `to`, keyed by the surface pointer.
    ///
    /// Deliberately independent of `edge_geom`: OCCT's CurveRepresentation
    /// list exists for edges that carry **only** pcurves.
    /// `ShapeFix_ComposeShell::SplitByLine` (`ShapeFix_ComposeShell.cxx:2061-2070`)
    /// creates exactly such an edge (no 3D curve, two pcurves + range), and
    /// `DispatchWires` copies it with `ShapeBuild_Edge::Copy`; routing that
    /// through `edge_geom` dropped the pcurves and the resulting face had
    /// mesh-less edges.
    pub fn copy_edge_pcurve_slots(&self, to: &TopoShape, from: &TopoShape) {
        if Arc::ptr_eq(&from.tshape, &to.tshape) {
            return;
        }
        let (curves, ranges, locations, continuities) = {
            let fts = from.tshape.read().expect("poisoned TShape lock");
            match fts.edge_pcurves() {
                Some(g) => (
                    g.curves.clone(),
                    g.ranges.clone(),
                    g.locations.clone(),
                    g.continuities.clone(),
                ),
                None => return,
            }
        };
        let mut tts = to.tshape.write().expect("poisoned TShape lock");
        let slot = tts.edge_pcurves_mut();
        // `ShapeBuild_Edge::CopyPCurves` (`cxx:364-412`) walks the source
        // representation list and appends each one the target does not carry
        // yet, so the target's list order follows the source's.
        for (k, v) in curves.iter() {
            slot.curves.insert(*k, v.clone());
        }
        for (k, v) in ranges.iter() {
            slot.ranges.insert(*k, *v);
        }
        for (k, v) in locations.iter() {
            slot.locations.insert(*k, v.clone());
        }
        for (k, v) in continuities.iter() {
            slot.continuities.insert(*k, *v);
        }
    }

    /// `BRep_Builder::Range(edge, face, first, last)` (`BRep_Builder.cxx:1121`).
    pub fn set_pcurve_range(&self, s: &TopoShape, face_key: usize, first: f64, last: f64) {
        let key = self.repr_key(face_key);
        s.tshape
            .write()
            .unwrap()
            .edge_pcurves_mut()
            .ranges
            .insert(key, (first, last));
    }

    /// COS representation `[First, Last]` for `face_key`, with the same
    /// same-surface fallback as [`GeometryRegistry::edge_pcurve`].
    pub fn pcurve_range(&self, s: &TopoShape, face_key: usize) -> Option<(f64, f64)> {
        let key = self.repr_key(face_key);
        let ts = s.tshape.read().unwrap();
        ts.edge_pcurves()?.ranges.get(&key).copied()
    }

    pub fn pcurve_range_on_face(&self, edge: &TopoShape, face: &TopoShape) -> Option<(f64, f64)> {
        let key = self.repr_key(Self::shape_key(face));
        let wanted = face.location().predivided(edge.location());
        let ts = edge.tshape.read().unwrap();
        let pcurves = ts.edge_pcurves()?;
        if !pcurves.location(key).is_equal(&wanted) {
            return None;
        }
        pcurves.ranges.get(&key).copied()
    }

    /// `BRep_Builder::Range(edge, face, first, last)` (`BRep_Builder.cxx:1121`)
    /// addressed by the `Geom_Surface` handle instead of a face TShape: the
    /// reference call sites that hold only the surface
    /// (`BRepFill_Sweep.cxx:1642`) cannot form a face key.
    pub fn set_pcurve_range_on_surface(
        &self,
        s: &TopoShape,
        surface: &Arc<dyn Surface>,
        first: f64,
        last: f64,
    ) {
        let key = Arc::as_ptr(surface) as *const () as usize;
        s.tshape
            .write()
            .unwrap()
            .edge_pcurves_mut()
            .ranges
            .insert(key, (first, last));
    }

    /// `BRep_Builder::SameRange`.
    pub fn set_same_range(&self, s: &TopoShape, value: bool) {
        s.tshape.write().unwrap().edge_core_mut().same_range = value;
    }

    /// `BRep_Builder::SameParameter`.
    pub fn set_same_parameter(&self, s: &TopoShape, value: bool) {
        s.tshape.write().unwrap().edge_core_mut().same_parameter = value;
    }

    /// `BRep_Builder::UpdateEdge` tolerance write.
    pub fn set_edge_tolerance(&self, s: &TopoShape, tol: f64) {
        s.tshape.write().unwrap().edge_core_mut().tolerance = tol;
    }

    /// CurveOnSurface representations: surface, pcurve, COS `[first, last]`.
    pub fn edge_pcurve_reps(
        &self,
        s: &TopoShape,
    ) -> Vec<(Arc<dyn Surface>, Arc<dyn Curve2d>, f64, f64)> {
        let by_ptr = self.surface_by_ptr.read().unwrap();
        // T-25 end phase: pcurves, their COS ranges and the 3D range all come
        // off the edge's own `TShape`.
        let (slot, three_d) = {
            let ts = s.tshape.read().unwrap();
            let Some(core) = ts.edge_core() else {
                return Vec::new();
            };
            (ts.edge_pcurves().cloned(), (core.first, core.last))
        };
        let Some(g) = slot else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for (&fk, pcs) in g.curves.iter() {
            let Some(surf) = by_ptr.get(&fk) else {
                continue;
            };
            let (a, b) = g.ranges.get(&fk).copied().unwrap_or(three_d);
            for pc in pcs {
                out.push((surf.clone(), pc.clone(), a, b));
            }
        }
        out
    }

    /// CurveOnSurface representations including
    /// `BRep_CurveRepresentation::Location`.
    pub fn edge_pcurve_reps_with_location(
        &self,
        s: &TopoShape,
    ) -> Vec<(Arc<dyn Surface>, Arc<dyn Curve2d>, f64, f64, TopLocLocation)> {
        let by_ptr = self.surface_by_ptr.read().unwrap();
        let (slot, three_d) = {
            let ts = s.tshape.read().unwrap();
            let Some(core) = ts.edge_core() else {
                return Vec::new();
            };
            (ts.edge_pcurves().cloned(), (core.first, core.last))
        };
        let Some(g) = slot else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for (&key, curves) in g.curves.iter() {
            let Some(surface) = by_ptr.get(&key) else {
                continue;
            };
            let (first, last) = g.ranges.get(&key).copied().unwrap_or(three_d);
            let location = g.location(key);
            for pcurve in curves {
                out.push((
                    surface.clone(),
                    pcurve.clone(),
                    first,
                    last,
                    location.clone(),
                ));
            }
        }
        out
    }

    /// `BRep_Builder::UpdateEdge(E, C2, S, L, Tol)` (`BRep_Builder.cxx:655-672`
    /// through `UpdateCurves` `cxx:104-167` and `cxx:149-164`): the pcurve of
    /// `edge` on the `Geom_Surface` `surface` is replaced, keeping the 3D
    /// representation's range; a null `curve` removes that representation
    /// instead (`cxx:149`: only a non-null `C` is appended). The edge
    /// tolerance is raised to `tol` without ever being lowered
    /// (`BRep_TEdge::UpdateTolerance`, `BRep_TEdge.lxx:33-37`).
    ///
    /// Unlike [`GeometryRegistry::set_edge_pcurve`] the key is the surface
    /// pointer itself, so a surface that no face has registered yet still
    /// anchors the representation.
    pub fn update_edge_pcurve_on_surface(
        &self,
        edge: &TopoShape,
        surface: &Arc<dyn Surface>,
        curve: Option<Arc<dyn Curve2d>>,
        tol: f64,
    ) {
        let key = Arc::as_ptr(surface) as *const () as usize;
        let (first, last) = self.edge_parameters(edge);
        let mut ts = edge.tshape.write().unwrap();
        match curve {
            Some(c) => ts.edge_pcurves_mut().set_pcurve(key, c, first, last),
            None => {
                let g = ts.edge_pcurves_mut();
                g.curves.remove(&key);
                g.ranges.remove(&key);
                g.locations.remove(&key);
            }
        }
        let c = ts.edge_core_mut();
        if c.tolerance < tol {
            c.tolerance = tol;
        }
    }

    /// The seam overload `BRep_Builder::UpdateEdge(E, C1, C2, S, L, Tol)`
    /// (`BRep_Builder.cxx:1194-1204`): both pcurves of a seam edge attached to
    /// the same `Geom_Surface` handle, keyed by that surface's pointer.
    pub fn set_edge_pcurves_on_surface(
        &self,
        edge: &TopoShape,
        surface: &Arc<dyn Surface>,
        curves: Vec<Arc<dyn Curve2d>>,
        tol: f64,
    ) {
        let key = Arc::as_ptr(surface) as *const () as usize;
        let (first, last) = self.edge_parameters(edge);
        let mut ts = edge.tshape.write().unwrap();
        ts.edge_pcurves_mut().set_pcurves(key, curves, first, last);
        let c = ts.edge_core_mut();
        if c.tolerance < tol {
            c.tolerance = tol;
        }
    }

    /// `ShapeBuild_Edge::RemovePCurve` — drop every pcurve of `edge` on the
    /// same `Geom_Surface` handle as `face`.
    pub fn remove_pcurves_on_surface(&self, edge: &TopoShape, face: &TopoShape) {
        let Some(want) = self.face_surface(face) else {
            return;
        };
        let want_key = Arc::as_ptr(&want) as *const () as usize;
        let mut ts = edge.tshape.write().unwrap();
        if let Some(g) = ts.edge_pcurves.as_mut() {
            g.curves.remove(&want_key);
            g.ranges.remove(&want_key);
            g.locations.remove(&want_key);
        }
    }

    /// Stable pointer-identity key for a shape, usable as a `HashMap` key.
    /// Same value as the internal `key` used to index the side-table.
    pub fn shape_key(s: &TopoShape) -> usize {
        key(s)
    }
}
