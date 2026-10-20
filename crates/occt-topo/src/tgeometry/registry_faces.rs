use super::*;

impl GeometryRegistry {

    pub fn set_face(&self, s: &TopoShape, geom: FaceGeom) {
        let k = key(s);
        self.ids.write().unwrap().insert(k, shape_id(s));
        // T-25 batch 4: the face geometry lives on the shape itself now.
        {
            let mut ts = s.tshape.write().unwrap();
            let c = ts.face_core_mut();
            c.surface = Some(geom.surface.clone());
            c.tolerance = geom.tolerance;
            c.natural_restriction = geom.natural_restriction;
        }
        self.face_surfaces.write().unwrap().insert(k, geom.surface.clone());
        self.surface_by_ptr
            .write()
            .unwrap()
            .insert(Arc::as_ptr(&geom.surface) as *const () as usize, geom.surface);
    }

    /// Canonical `BRep_GCurve` representation key: OCCT keys a CurveOnSurface
    /// representation by the `Geom_Surface` handle (plus location); the port
    /// has identity locations, so the registered surface's data pointer is the
    /// key. Faces that share a surface therefore share one representation, as
    /// in OCCT.
    pub(super) fn repr_key(&self, face_key: usize) -> usize {
        self.face_surfaces
            .read()
            .unwrap()
            .get(&face_key)
            .map(|s| Arc::as_ptr(s) as *const () as usize)
            .unwrap_or(face_key)
    }

    /// Debug accessor for `repr_key`, for the in-tree pcurve diagnostics.
    pub fn debug_repr_key(&self, face_key: usize) -> usize {
        self.repr_key(face_key)
    }

    /// (face_key, surface_ptr) pairs currently registered.
    pub fn debug_face_surfaces(&self) -> Vec<(usize, usize)> {
        self.face_surfaces
            .read()
            .unwrap()
            .iter()
            .map(|(k, s)| (*k, Arc::as_ptr(s) as *const () as usize))
            .collect()
    }

    /// Kind + u/v range of the surface registered under `ptr`.
    pub fn debug_surface_desc(&self, ptr: usize) -> String {
        let by_ptr = self.surface_by_ptr.read().unwrap();
        match by_ptr.get(&ptr) {
            Some(s) => {
                let (u0, u1) = s.u_range();
                let (v0, v1) = s.v_range();
                format!(
                    "{:?} u=[{u0:.6},{u1:.6}] v=[{v0:.6},{v1:.6}]",
                    crate::pcurve_full::classify_surface_kind(s.as_ref())
                )
            }
            None => format!("ptr{ptr}:unregistered"),
        }
    }

    pub fn face_geom(&self, s: &TopoShape) -> Option<FaceGeom> {
        // T-25 batch 4: read off the face's own `TShape`; the side table is the
        // safety net while other slots are still being lifted.
        // T-25 end phase: the face geometry comes off its own `TShape`.
        let ts = s.tshape.read().unwrap();
        let c = ts.face_core()?;
        let surf = c.surface.clone()?;
        Some(FaceGeom {
            surface: surf,
            tolerance: c.tolerance,
            natural_restriction: c.natural_restriction,
        })
    }

    /// The underlying surface handle (clone of the Arc).
    pub fn face_surface(&self, s: &TopoShape) -> Option<Arc<dyn Surface>> {
        self.face_geom(s).map(|g| g.surface)
    }

    pub fn face_tolerance(&self, s: &TopoShape) -> f64 {
        // `BRep_Tool::Tolerance(face)` -> `BRep_TFace::Tolerance()`, whose
        // initial value is `RealEpsilon()` (`BRep_TFace.cxx:29`).
        self.face_geom(s)
            .map(|g| g.tolerance)
            .unwrap_or(occt_core::precision::REAL_EPSILON)
    }

    pub fn natural_restriction(&self, s: &TopoShape) -> bool {
        self.face_geom(s).map(|g| g.natural_restriction).unwrap_or(true)
    }

    /// `BRep_Builder::NaturalRestriction(face, flag)`.
    pub fn set_natural_restriction(&self, s: &TopoShape, flag: bool) {
        s.tshape.write().unwrap().face_core_mut().natural_restriction = flag;
    }

    /// `BRep_Builder::UpdateFace` tolerance write.
    pub fn set_face_tolerance(&self, s: &TopoShape, tol: f64) {
        s.tshape.write().unwrap().face_core_mut().tolerance = tol;
    }

    /// Register the surface side-table entry for a face whose `TFace` surface
    /// was written directly (a rebuild through `TopoDS_Shape::EmptyCopied`,
    /// `BRep_TFace.cxx:37-43`). `repr_key` resolves pcurves through
    /// `face_surfaces`, so a rebuilt face whose entry is missing would fall back
    /// to the face pointer and lose every pcurve: OCCT instead finds a
    /// CurveOnSurface representation by the face's own surface
    /// (`BRep_Tool::CurveOnSurface`, `BRep_Tool.cxx:347-357`), which survives a
    /// face-only rebuild because the surface handle is carried over.
    pub fn register_face_surface(&self, s: &TopoShape, surface: &Arc<dyn Surface>) {
        let k = key(s);
        self.face_surfaces
            .write()
            .unwrap()
            .insert(k, surface.clone());
    }

}
