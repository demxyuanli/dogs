//! `TrslUV` / `GetnewS` / `AddNewConnexity` / `RemoveOldConnexity` of
//! `TopOpeBRepTool_CORRISO`.
//!
//! Ported from `TopOpeBRepTool_CORRISO.cxx:1161-1399`.

use super::*;

impl CorrIso {
    /// `TrslUV(onU, FyEds)` (`cxx:1161-1215`): shift every faulty pcurve by one
    /// period, in the direction `FyEds` asks for.
    pub fn trsl_uv(&mut self, on_u: bool, fy_eds: &[(TopoShape, i32)]) -> bool {
        let tt2d = if on_u {
            match self.ref_closed(1) {
                Some(uper) => GpVec2d::new(uper, 0.0),
                None => return false,
            }
        } else {
            match self.ref_closed(2) {
                Some(vper) => GpVec2d::new(0.0, vper),
                None => return false,
            }
        };

        for (e, itt) in fy_eds.iter() {
            let Some(mut c2df) = self.uv_rep(e) else {
                return false;
            };
            let v = match *itt {
                SPLITEDGE => return false,
                INCREASE => tt2d,
                DECREASE => tt2d.multiplied_scalar(-1.0),
                _ => return false,
            };
            trsl_uv_c2df(&v, &mut c2df);
            self.set_uv_rep(e, &c2df);
        }
        true
    }

    /// `GetnewS(newS)` (`cxx:1236-1287`): write the stored 2-D representations
    /// back onto `S` as trimmed curves of the face.
    ///
    /// `BB.UpdateEdge(E, cu, curr, newS, tol)` (`cxx:1279-1281`) is the
    /// two-pcurve seam form, `BB.UpdateEdge(E, cu, newS, tol)`
    /// (`cxx:1284`) the ordinary one.
    pub fn get_new_s(&self) -> Option<TopoShape> {
        let s = self.s.as_ref()?;
        if s.shape_type() != ShapeType::Face {
            return None;
        }
        let new_s = s.clone();
        let reg = GeometryRegistry::global();
        let face_key = GeometryRegistry::shape_key(&new_s);

        for e in self.eds.iter() {
            let ori_e = e.orientation();
            let c2df = self.uv_rep(e)?;
            let (pc, f, l, _tol) = c2df.pc();
            let pc = pc?;
            let cu: Arc<dyn Curve2d> = Arc::new(Geom2dTrimmedCurve::new(pc, f, l));

            let err = e.clone().oriented(ori_e.complement());
            if let Some(c2df_rr) = self.uv_rep(&err) {
                let (pcrr, frr, lrr, _tolrr) = c2df_rr.pc();
                if let Some(pcrr) = pcrr {
                    let curr: Arc<dyn Curve2d> =
                        Arc::new(Geom2dTrimmedCurve::new(pcrr, frr, lrr));
                    if ori_e.is_forward() {
                        reg.set_edge_pcurves(e, face_key, vec![cu, curr]);
                    }
                }
            } else {
                reg.set_edge_pcurve(e, face_key, cu);
            }
        }
        Some(new_s)
    }

    /// `AddNewConnexity(V, E)` (`cxx:1293-1340`): register `E` (and its pcurve)
    /// if the tool does not know it yet. `V` is unused in OCCT.
    pub fn add_new_connexity(&mut self, e: &TopoShape) -> bool {
        let isb = self.e_rep_2d.contains_key(&shape_map_key(e));

        if !isb {
            let has_old = fc2d_has_old_curve_on_surface(e, &self.f_ref).is_some();
            let Some(c2df) = fc2d_editable_curve_on_surface(e, &self.f_ref, false) else {
                return false;
            };
            if !has_old {
                if let Some(pc) = c2df.curve.clone() {
                    fc2d_add_new_curve_on_surface(
                        &pc, e, &self.f_ref, c2df.first, c2df.last, c2df.tol,
                    );
                }
            }
            self.e_rep_2d.insert(shape_map_key(e), c2df);
        }

        if !isb {
            self.eds.push(e.clone());
        }

        for v in explorer_shapes(e, ShapeType::Vertex) {
            self.v_eds
                .entry(shape_map_key(&v))
                .or_default()
                .push(e.clone());
        }
        true
    }

    /// `RemoveOldConnexity(V, E)` (`cxx:1344-1399`): forget `E`'s pcurve, its
    /// place in `myEds` and its entry under each of its vertices. `V` is unused
    /// in OCCT.
    ///
    /// OCCT answers the `done` flag of the last removal it attempted — `true`
    /// when the edge was found under any vertex, even after it was already
    /// removed from the other structures.
    pub fn remove_old_connexity(&mut self, e: &TopoShape) -> bool {
        let isb = self.e_rep_2d.contains_key(&shape_map_key(e));
        if isb {
            self.e_rep_2d.remove(&shape_map_key(e));
        }

        if isb {
            if let Some(pos) = self.eds.iter().position(|x| is_equal(x, e)) {
                self.eds.remove(pos);
            }
        }

        let mut done = false;
        for v in explorer_shapes(e, ShapeType::Vertex) {
            let key = shape_map_key(&v);
            let Some(loe) = self.v_eds.get_mut(&key) else {
                return false;
            };
            if let Some(pos) = loe.iter().position(|x| is_equal(x, e)) {
                done = true;
                loe.remove(pos);
            }
        }
        done
    }
}

/// `TopoDS_Shape::IsEqual` (`TopoDS_Shape.hxx:270-275`): same `TShape`, same
/// location **and** same orientation.
fn is_equal(a: &TopoShape, b: &TopoShape) -> bool {
    crate::topo_tools::is_equal(a, b)
}
