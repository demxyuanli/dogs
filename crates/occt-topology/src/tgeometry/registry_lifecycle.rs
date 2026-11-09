use super::*;

impl GeometryRegistry {

    /// Drop all geometry entries owned by `s`. Call when a shape is discarded
    /// to keep the side-table from growing without bound.
    pub fn clear_shape(&self, s: &TopoShape) {
        let k = key(s);

        self.face_surfaces.write().unwrap().remove(&k);
        // T-25: the edge geometry lives on the shape itself now.
        {
            let mut ts = s.tshape.write().unwrap();
            ts.edge_pcurves = None;
            ts.edge_core = None;
            ts.vertex_core = None;
            ts.face_core = None;
        }
    }

    /// Remove every geometry entry keyed by the raw `TShape` address. Called
    /// from `TShape::drop` so entries die with their shape — this prevents a
    /// stale entry from leaking into a future shape that reuses the address.
    pub fn remove_by_ptr(&self, ptr: usize, id: u64) {
        {
            let mut ids = self.ids.write().unwrap();
            if ids.get(&ptr).copied() != Some(id) {
                // The address has been reused by a later shape: its geometry is
                // not ours to erase.
                return;
            }
            ids.remove(&ptr);
        }

        self.face_surfaces.write().unwrap().remove(&ptr);
    }

    /// Number of live entries (vertices + edges + faces).
    pub fn len(&self) -> usize {
        self.ids.read().unwrap().len()
    }

    pub fn is_empty(&self) -> bool { self.len() == 0 }

    /// Remove every entry (for tests / teardown).
    pub fn clear_all(&self) {
        self.ids.write().unwrap().clear();

        self.face_surfaces.write().unwrap().clear();
    }
}
