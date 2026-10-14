//! Top-level `StepToTopoDS_TranslateTool::Bind` dispatch over the STEP entity type, plus `StepToTopoDS_TranslateVertex`.
use super::*;

impl<'a> Resolver<'a> {
    pub(in crate::step) fn resolve_shape(&self, id: usize) -> Result<TopoShape, String> {
        if let Some(s) = self.shape_cache.borrow().get(&id) {
            return Ok(s.clone());
        }
        let rec = self.record(id)?;
        // `StepToTopoDS_TranslateTool::Bind` — a `VERTEX_POINT` another record was
        // bound onto (`StepToTopoDS_TranslateEdgeLoop.cxx:384-396`, `:466-477`)
        // resolves to that vertex instead of a fresh one.
        if rec.type_name == "VERTEX_POINT" {
            if let Some(v) = self.vertex_bind.borrow().get(&id) {
                return Ok(v.clone());
            }
        }
        if !self.resolving.borrow_mut().insert(id) {
            return Err(format!("cyclic entity reference #{id}"));
        }
        let shape = match rec.type_name.as_str() {
            "VERTEX_POINT" => self.resolve_vertex(rec),
            "EDGE_CURVE" => self.resolve_edge(rec),
            "ORIENTED_EDGE" => self.resolve_oriented_edge(rec),
            "EDGE_LOOP" => self.resolve_loop(rec),
            "VERTEX_LOOP" => self.resolve_vertex_loop(rec),
            "FACE_OUTER_BOUND" => self.resolve_outer_bound(rec),
            "FACE_BOUND" => self.resolve_outer_bound(rec),
            "ADVANCED_FACE" => self.resolve_face(rec),
            "CLOSED_SHELL" => self.resolve_shell(rec),
            "MANIFOLD_SOLID_BREP" => self.resolve_solid(rec),
            "BREP_WITH_VOIDS" => self.resolve_brep_with_voids(rec),
            "MAPPED_ITEM" => self.resolve_mapped_item(rec),
            other => {
                self.warn(format!("unsupported topological entity {other} (#{id})"));
                Err(format!("unsupported entity {other} (#{id})"))
            }
        };
        self.resolving.borrow_mut().remove(&id);
        match shape {
            Ok(s) => {
                self.shape_cache.borrow_mut().insert(id, s.clone());
                Ok(s)
            }
            Err(e) => Err(e),
        }
    }

    pub(in crate::step) fn resolve_vertex(&self, rec: &'a Record) -> Result<TopoShape, String> {
        if rec.args.len() < 2 {
            return Err("VERTEX_POINT: bad args".into());
        }
        let pid = parse_ref(&rec.args[1]).ok_or("VERTEX_POINT: bad point ref")?;
        let p = self.resolve_point(pid)?;
        Ok(self.b.make_vertex(p, 0.0).0)
    }

}
