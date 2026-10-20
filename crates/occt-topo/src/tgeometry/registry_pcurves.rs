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

    /// Attach a pcurve to edge `s` for the face identified by `face_key`.
    /// Mirrors `BRep_Builder::UpdateEdge(edge, curve2d, face, tol)`.
    pub fn set_edge_pcurve(&self, s: &TopoShape, face_key: usize, curve: Arc<dyn Curve2d>) {
        if std::env::var("PCSET").is_ok() {
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
        if std::env::var("PCBTRACE").is_ok()
            && curve.bspline_copy2d().map(|b| b.xs.len()) == Some(25)
            && curve.first_parameter() < -0.2
        {
            let bt = std::backtrace::Backtrace::force_capture();
            let text = format!("{bt}");
            let mut names: Vec<String> = Vec::new();
            for line in text.lines() {
                if let Some(pos) = line.find("occt_topo::") {
                    let f = line[pos..].trim().to_string();
                    names.push(f);
                    if names.len() >= 8 {
                        break;
                    }
                }
            }
            eprintln!("PCBTRACE edge={:p} {}", Arc::as_ptr(&s.tshape), names.join(" <- "));
        }
        let (first, last) = self.edge_parameters(s);
        let key = self.repr_key(face_key);
        s.tshape.write().unwrap().edge_pcurves_mut().set_pcurve(key, curve, first, last);
    }

    /// Replace the pcurves of edge `s` on the face identified by `face_key`.
    /// Mirrors the seam overload `BRep_Builder::UpdateEdge(edge, c1, c2, face)`.
    pub fn set_edge_pcurves(&self, s: &TopoShape, face_key: usize, curves: Vec<Arc<dyn Curve2d>>) {
        if std::env::var("PCSET").is_ok() {
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
        if std::env::var("PCS2DBG").is_ok() && curves.len() >= 2 {
            let (e0, e1) = self.edge_parameters(s);
            let mut line = format!("PCS2 edge={:p} face={face_key} n={}", Arc::as_ptr(&s.tshape), curves.len());
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
        if std::env::var("PCBTRACE2").is_ok() {
            let hit = curves.iter().any(|c| {
                c.bspline_copy2d().map(|b| b.xs.len()) == Some(25)
                    && (c.first_parameter() - 1.1414320265529101).abs() < 1e-12
            });
            if hit {
                let bt = std::backtrace::Backtrace::force_capture();
                let text = format!("{bt}");
                let mut names: Vec<String> = Vec::new();
                for line in text.lines() {
                    if let Some(pos) = line.find("occt_topo::") {
                        names.push(line[pos..].trim().to_string());
                        if names.len() >= 12 {
                            break;
                        }
                    }
                }
                for (q, c) in curves.iter().enumerate() {
                    let (f, l) = (c.first_parameter(), c.last_parameter());
                    match c.bspline_copy2d() {
                        Some(bs) => {
                            let ys: Vec<String> =
                                bs.ys.iter().map(|y| format!("{y:.17e}")).collect();
                            let xs: Vec<String> =
                                bs.xs.iter().take(3).map(|x| format!("{x:.17e}")).collect();
                            eprintln!(
                                "PCB2 edge={:p} face={face_key} q={q} f={f:.17e} l={l:.17e} \
deg={} np={} nk={} x[:3]=[{}] y=[{}] {}",
                                Arc::as_ptr(&s.tshape),
                                bs.degree(),
                                bs.xs.len(),
                                bs.knots.len(),
                                xs.join(","),
                                ys.join(","),
                                names.join(" <- ")
                            );
                        }
                        None => eprintln!(
                            "PCB2 edge={:p} face={face_key} q={q} f={f:.17e} l={l:.17e} nobs {}",
                            Arc::as_ptr(&s.tshape),
                            names.join(" <- ")
                        ),
                    }
                }
            }
        }
        if std::env::var("PCBTRACE").is_ok() {
            for c in &curves {
                if c.bspline_copy2d().map(|b| b.xs.len()) == Some(25)
                    && (c.first_parameter() - 1.1414320265529101).abs() < 1e-12
                {
                    let bt = std::backtrace::Backtrace::force_capture();
                    let text = format!("{bt}");
                    let mut names: Vec<String> = Vec::new();
                    for line in text.lines() {
                        if let Some(pos) = line.find("occt_topo::") {
                            names.push(line[pos..].trim().to_string());
                            if names.len() >= 10 {
                                break;
                            }
                        }
                    }
                    let mut info = String::new();
                    for c in &curves {
                        let (f, l) = (c.first_parameter(), c.last_parameter());
                        let a = c.d0(f);
                        let m = c.d0(0.5 * (f + l));
                        let b = c.d0(l);
                        let poles = c
                            .bspline_copy2d()
                            .map(|bs| {
                                let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
                                for y in &bs.ys {
                                    lo = lo.min(*y);
                                    hi = hi.max(*y);
                                }
                                let head: Vec<String> =
                                    bs.ys.iter().take(4).map(|y| format!("{y:.17e}")).collect();
                                format!(" pymin={lo:.17e} pymax={hi:.17e} py=[{}]", head.join(","))
                            })
                            .unwrap_or_default();
                        info.push_str(&format!(
                            " [f={f:.17e} y0={:.17e} ym={:.17e} y1={:.17e} x0={:.17e}{poles}]",
                            a.y(),
                            m.y(),
                            b.y(),
                            a.x()
                        ));
                    }
                    eprintln!(
                        "PBTRACE edge={:p} n={} face={face_key}{info} {}",
                        Arc::as_ptr(&s.tshape),
                        curves.len(),
                        names.join(" <- ")
                    );
                }
            }
        }
        let (first, last) = self.edge_parameters(s);
        let key = self.repr_key(face_key);
        s.tshape.write().unwrap().edge_pcurves_mut().set_pcurves(key, curves, first, last);
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
        let (curves, ranges) = {
            let fts = from.tshape.read().expect("poisoned TShape lock");
            match fts.edge_pcurves() {
                Some(g) => (g.curves.clone(), g.ranges.clone()),
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
        }
    }

    /// Stable pointer-identity key for a shape, usable as a `HashMap` key.
    /// Same value as the internal `key` used to index the side-table.
    pub fn shape_key(s: &TopoShape) -> usize {
        key(s)
    }

}
