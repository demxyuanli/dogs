use super::prelude::*;
use super::*;

impl RangeSplitter for UndefinedRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner.inner.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner.inner.inner
    }

    fn adjust_range(&mut self) {
        self.inner.adjust_range();
    }

    fn generate_surface_nodes(&self, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        generate_nurbs_grid(self, params)
    }

    fn parameters_u(&self) -> Option<&ParamSet> {
        Some(&self.inner.inner.u_params)
    }
    fn parameters_v(&self) -> Option<&ParamSet> {
        Some(&self.inner.inner.v_params)
    }
    fn parameters_u_mut(&mut self) -> Option<&mut ParamSet> {
        Some(&mut self.inner.inner.u_params)
    }
    fn parameters_v_mut(&mut self) -> Option<&mut ParamSet> {
        Some(&mut self.inner.inner.v_params)
    }

    fn get_undefined_interval_nb(&self, _is_u: bool, _continuity: u8) -> i32 {
        1
    }
}

/// Splitter that seeds the U/V parameter sets from boundary points only. Source:
/// `BRepMesh_BoundaryParamsRangeSplitter`.
pub struct BoundaryParamsRangeSplitter {
    pub(super) inner: NURBSRangeSplitter,
}

impl BoundaryParamsRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: NURBSRangeSplitter::new(),
        }
    }
}

impl RangeSplitter for BoundaryParamsRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner.inner.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner.inner.inner
    }

    fn reset(&mut self, dface: &MeshFace, params: &MeshParameters) {
        self.reset_base(dface, params);
        // The derived boundary/torus splitters seed their own UV maps from
        // AddPoint; a fresh Reset must clear them (the base reset below only
        // resets the bottom DefaultRangeSplitter).
        self.inner.inner.u_params.clear();
        self.inner.inner.v_params.clear();
    }

    fn add_point(&mut self, point: GpPnt2d) {
        self.add_point_base(point);
        self.inner.inner.u_params.insert(point.x());
        self.inner.inner.v_params.insert(point.y());
    }

    fn adjust_range(&mut self) {
        self.inner.adjust_range();
    }

    fn generate_surface_nodes(&self, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        generate_nurbs_grid(self, params)
    }

    fn parameters_u(&self) -> Option<&ParamSet> {
        Some(&self.inner.inner.u_params)
    }
    fn parameters_v(&self) -> Option<&ParamSet> {
        Some(&self.inner.inner.v_params)
    }
    fn parameters_u_mut(&mut self) -> Option<&mut ParamSet> {
        Some(&mut self.inner.inner.u_params)
    }
    fn parameters_v_mut(&mut self) -> Option<&mut ParamSet> {
        Some(&mut self.inner.inner.v_params)
    }

    /// Source: `BRepMesh_BoundaryParamsRangeSplitter::initParameters`.
    fn init_parameters(&self) -> bool {
        true
    }

    fn seeds_cn_intervals(&self) -> bool {
        false
    }
}

/// Splitter for extrusion surfaces — the interval count follows the basis
/// curve, simplified here to a single interval. Source:
/// `BRepMesh_ExtrusionRangeSplitter`.
pub struct ExtrusionRangeSplitter {
    pub(super) inner: NURBSRangeSplitter,
}

impl ExtrusionRangeSplitter {
    /// Creates an empty splitter.
    pub fn new() -> Self {
        Self {
            inner: NURBSRangeSplitter::new(),
        }
    }
}

impl RangeSplitter for ExtrusionRangeSplitter {
    fn base(&self) -> &DefaultRangeSplitter {
        &self.inner.inner.inner
    }
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter {
        &mut self.inner.inner.inner
    }

    fn adjust_range(&mut self) {
        self.inner.adjust_range();
    }

    fn generate_surface_nodes(&self, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        generate_nurbs_grid(self, params)
    }

    fn parameters_u(&self) -> Option<&ParamSet> {
        Some(&self.inner.inner.u_params)
    }
    fn parameters_v(&self) -> Option<&ParamSet> {
        Some(&self.inner.inner.v_params)
    }
    fn parameters_u_mut(&mut self) -> Option<&mut ParamSet> {
        Some(&mut self.inner.inner.u_params)
    }
    fn parameters_v_mut(&mut self) -> Option<&mut ParamSet> {
        Some(&mut self.inner.inner.v_params)
    }

    fn get_undefined_interval_nb(&self, _is_u: bool, continuity: u8) -> i32 {
        // `BRepMesh_ExtrusionRangeSplitter.cxx:24-47`: BasisCurve NbIntervals,
        // then NbPoles-1 for Bezier/BSpline when the curve is a single span.
        let Some(surf) = self.surface() else {
            return 1;
        };
        let curve = surf.extrusion_basis_curve().or_else(|| {
            let (v0, _) = surf.v_range();
            surf.u_iso_curve(v0)
        });
        let Some(curve) = curve else {
            return 1;
        };
        let mut n = curve.nb_intervals(continuity);
        if n == 1 {
            if let Some(poles) = curve.bspline_poles().or_else(|| curve.bezier_poles()) {
                n = (poles.len() as i32) - 1;
            }
        }
        n.max(1)
    }
}

/// Creates the range splitter matching a surface's analytic type. Source:
/// `BRepMesh_MeshAlgoFactory` splitter selection.
pub fn create_range_splitter(surface: &dyn Surface) -> Box<dyn RangeSplitter> {
    match classify_surface(surface) {
        SurfaceType::Plane => Box::new(DefaultRangeSplitter::new()),
        SurfaceType::Sphere => Box::new(SphereRangeSplitter::new()),
        SurfaceType::Cylinder => Box::new(CylinderRangeSplitter::new()),
        SurfaceType::Cone => Box::new(ConeRangeSplitter::new()),
        SurfaceType::Torus => Box::new(TorusRangeSplitter::new()),
        SurfaceType::SurfaceOfRevolution => Box::new(BoundaryParamsRangeSplitter::new()),
        SurfaceType::SurfaceOfExtrusion => Box::new(ExtrusionRangeSplitter::new()),
        SurfaceType::BezierSurface | SurfaceType::BSplineSurface => Box::new(NURBSRangeSplitter::new()),
        SurfaceType::OffsetSurface | SurfaceType::OtherSurface => {
            Box::new(UndefinedRangeSplitter::new())
        }
    }
}

/// Whether `BRepMesh_MeshAlgoFactory::GetAlgo` wraps the splitter in
/// `DelaunayDeflectionControlMeshAlgo` (`cxx:60-117`).
///
/// Revolution / extrusion / Bezier / BSpline / Offset / Other always use it.
/// Plane / sphere / cylinder / cone / torus use it only when
/// `EnableControlSurfaceDeflectionAllSurfaces` is set.
pub fn factory_uses_deflection_control(surface: &dyn Surface, params: &MeshParameters) -> bool {
    if params.enable_control_surface_deflection_all_surfaces {
        return true;
    }
    match classify_surface(surface) {
        SurfaceType::SurfaceOfRevolution
        | SurfaceType::SurfaceOfExtrusion
        | SurfaceType::BezierSurface
        | SurfaceType::BSplineSurface
        | SurfaceType::OffsetSurface
        | SurfaceType::OtherSurface => true,
        SurfaceType::Plane
        | SurfaceType::Sphere
        | SurfaceType::Cylinder
        | SurfaceType::Cone
        | SurfaceType::Torus => false,
    }
}

/// `BRepMesh_NURBSRangeSplitter::grabParamsOfEdges(Edge_Internal, Param_U|Param_V)`
/// (`cxx:492-536`). Frontier samples are not added here — that arm is commented
/// out in OCCT. `BoundaryParamsRangeSplitter::initParameters` skips this.
pub fn grab_params_of_internal_edges(
    model: &MeshModel,
    face_index: usize,
    splitter: &mut dyn RangeSplitter,
) {
    if !splitter.seeds_cn_intervals() {
        return;
    }
    let Ok(face) = model.face(face_index) else {
        return;
    };
    for &wire_index in face.wires() {
        let Ok(wire) = model.wire(wire_index) else {
            continue;
        };
        for j in 0..wire.edges_nb() {
            let Ok(edge_index) = wire.edge(j) else {
                continue;
            };
            let Ok(edge) = model.edge(edge_index) else {
                continue;
            };
            for p in 0..edge.pcurves_nb() {
                let Ok(pc) = edge.pcurve(p) else {
                    continue;
                };
                if pc.face() != face_index || !pc.is_internal() {
                    continue;
                }
                for pt in pc.points() {
                    if let Some(u) = splitter.parameters_u_mut() {
                        u.insert(pt.x());
                    }
                    if let Some(v) = splitter.parameters_v_mut() {
                        v.insert(pt.y());
                    }
                }
            }
        }
    }
}
