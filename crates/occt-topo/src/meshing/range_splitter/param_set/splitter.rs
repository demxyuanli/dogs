use super::intervals::update_range;
use super::set::ParamSet;
use super::super::prelude::*;

/// Range splitter — computes the discrete UV range of a face and, for the
/// analytical / NURBS splitters, generates interior (surface) nodes.
/// Source: `BRepMesh_DefaultRangeSplitter`.
pub trait RangeSplitter {
    /// Access to the embedded base splitter state.
    fn base(&self) -> &DefaultRangeSplitter;
    /// Mutable access to the embedded base splitter state.
    fn base_mut(&mut self) -> &mut DefaultRangeSplitter;

    // ---- lifecycle ----

    /// Resets the splitter for the given discrete face. Must be called before
    /// first use. Source: `Reset`.
    fn reset(&mut self, dface: &MeshFace, params: &MeshParameters) {
        self.reset_base(dface, params);
    }

    /// Base implementation of `Reset` (used by derived splitters).
    fn reset_base(&mut self, dface: &MeshFace, params: &MeshParameters) {
        self.base_mut().reset(dface, params);
    }

    /// Registers a boundary point. Source: `AddPoint`.
    fn add_point(&mut self, point: GpPnt2d) {
        self.add_point_base(point);
    }

    /// Base implementation of `AddPoint`.
    fn add_point_base(&mut self, point: GpPnt2d) {
        self.base_mut().add_point(point);
    }

    /// Computes the discrete range from the geometric range and the registered
    /// boundary points. Source: `AdjustRange`.
    fn adjust_range(&mut self) {
        self.adjust_range_base();
    }

    /// Base implementation of `AdjustRange` (used by derived splitters).
    fn adjust_range_base(&mut self) {
        let (mut du_first, mut du_second) = self.base().range_u();
        let (mut dv_first, mut dv_second) = self.base().range_v();
        let surface = self.surface().cloned();
        let Some(surf) = surface else {
            self.base_mut().set_valid(false);
            return;
        };

        // `BRepMesh_DefaultRangeSplitter::GetSurface()` (`hxx:81`) is the
        // discrete face's own adaptor, which `IMeshData_Face` builds
        // UNRESTRICTED (`IMeshData_Face.hxx:69`
        // `new BRepAdaptor_Surface(GetFace(), false)`). Its
        // `First/Last U|VParameter` therefore come from `Geom_Surface::Bounds`
        // (`BRepAdaptor_Surface.cxx:79` `Load(aSurface, trsf)`), NOT from
        // `BRepTools::UVBounds`. `AdjustRange` (`cxx:45-81`) must use that plain
        // surface range; clipping the discrete range to the face's pcurve box is
        // the `Restriction=true` adaptor's behaviour and collapses faces whose
        // stored pcurves lie off the surface domain (e.g. a plane face whose
        // pcurves sit in a shifted parameterization), which leaves
        // `computeLengthV` == 0 and `IsValid` == false.
        let (gu0, gu1) = surf.u_range();
        let (gv0, gv1) = surf.v_range();
        update_range(gu0, gu1, surf.is_u_periodic(), &mut du_first, &mut du_second);
        if du_second < du_first {
            self.base_mut().set_valid(false);
            return;
        }

        update_range(gv0, gv1, surf.is_v_periodic(), &mut dv_first, &mut dv_second);
        if dv_second < dv_first {
            self.base_mut().set_valid(false);
            return;
        }

        self.base_mut().set_range_u((du_first, du_second));
        self.base_mut().set_range_v((dv_first, dv_second));

        let len_u = self.base().compute_length_u(surf.as_ref());
        let len_v = self.base().compute_length_v(surf.as_ref());
        let valid = len_u > PCONFUSION && len_v > PCONFUSION;
        self.base_mut().set_valid(valid);
        if valid {
            self.compute_tolerance(len_u, len_v);
            self.compute_delta(len_u, len_v);
        }
        // TEMPORARY diagnostic (SPLITDBG): raw splitter inputs/outputs so they can
        // be diffed against `occt_probe --uvsum` with `OCCT_TOPO_TRACE_SPLITTER=1`.
        if std::env::var("SPLITDBG").is_ok() {
            let (tu, tv) = self.base().tolerance_uv();
            let (da, db) = self.base().delta();
            let fkey = self
                .base()
                .dface
                .as_ref()
                .map(|df| crate::tgeometry::GeometryRegistry::shape_key(&df.face().0))
                .unwrap_or(0);
            eprintln!(
                "SPLITDBG key={fkey} ru=[{:.17e},{:.17e}] rv=[{:.17e},{:.17e}] len_u={:.17e} len_v={:.17e} tol=[{:.17e},{:.17e}] delta=[{:.17e},{:.17e}] defl={:.17e} valid={}",
                self.base().range_u().0,
                self.base().range_u().1,
                self.base().range_v().0,
                self.base().range_v().1,
                len_u,
                len_v,
                tu,
                tv,
                da,
                db,
                self.base().deflection(),
                valid as i32
            );
        }
    }

    /// True when the computed range is valid. Source: `IsValid`.
    fn is_valid(&self) -> bool {
        self.base().is_valid()
    }

    /// Scales a point between real parametric space and the face basis.
    /// Source: `Scale`.
    fn scale(&self, point: GpPnt2d, to_face_basis: bool) -> GpPnt2d {
        self.base().scale(point, to_face_basis)
    }

    /// Returns the interior nodes generated from the surface data, or `None`
    /// when the splitter generates no nodes. Source: `GenerateSurfaceNodes`.
    fn generate_surface_nodes(&self, params: &MeshParameters) -> Option<Vec<GpPnt2d>> {
        self.base().generate_surface_nodes(params)
    }

    /// 3D point corresponding to the given UV parameter. Source: `Point`.
    fn point(&self, point2d: GpPnt2d) -> GpPnt {
        self.base().point(point2d)
    }

    // ---- accessors ----

    /// Discrete U range.
    fn range_u(&self) -> (f64, f64) {
        self.base().range_u()
    }
    /// Discrete V range.
    fn range_v(&self) -> (f64, f64) {
        self.base().range_v()
    }
    /// Scale factors per direction.
    fn delta(&self) -> (f64, f64) {
        self.base().delta()
    }
    /// Parametric tolerances per direction.
    fn tolerance_uv(&self) -> (f64, f64) {
        self.base().tolerance_uv()
    }
    /// The discrete face model.
    fn dface(&self) -> Option<&MeshFace> {
        self.base().dface()
    }
    /// The face surface.
    fn surface(&self) -> Option<&Arc<dyn Surface>> {
        self.base().surface()
    }
    /// The face deflection.
    fn deflection(&self) -> f64 {
        self.base().deflection()
    }

    // ---- UV-parameter maps (UVParamRangeSplitter and derivatives) ----

    /// Sorted U parameters collected from boundary points, when the splitter
    /// tracks them.
    fn parameters_u(&self) -> Option<&ParamSet> {
        None
    }
    /// Sorted V parameters collected from boundary points, when the splitter
    /// tracks them.
    fn parameters_v(&self) -> Option<&ParamSet> {
        None
    }
    /// Mutable U parameter map. Used by `grabParamsOfEdges`.
    fn parameters_u_mut(&mut self) -> Option<&mut ParamSet> {
        None
    }
    /// Mutable V parameter map. Used by `grabParamsOfEdges`.
    fn parameters_v_mut(&mut self) -> Option<&mut ParamSet> {
        None
    }

    // ---- virtual hooks overridden by derived splitters ----

    /// Computes the parametric tolerances. Source: `computeTolerance`.
    fn compute_tolerance(&mut self, len_u: f64, len_v: f64) {
        self.base_mut().compute_tolerance(len_u, len_v);
    }

    /// Computes the scale factors. Source: `computeDelta`.
    fn compute_delta(&mut self, len_u: f64, len_v: f64) {
        self.base_mut().compute_delta(len_u, len_v);
    }

    /// Number of intervals to subdivide an undefined range into. Source:
    /// `getUndefinedIntervalNb` — `(isU ? NbUPoles : NbVPoles) - 1` for NURBS.
    fn undefined_interval_nb(&self, _is_u: bool, continuity: u8) -> i32 {
        ((continuity as i32) + 1).clamp(2, 32)
    }

    /// Whether the surface parameter sets can be initialized. Source:
    /// `initParameters`.
    fn init_parameters(&self) -> bool {
        self.surface().is_some()
    }

    /// True when `initParameters` seeds U/V from CN intervals.
    /// `BRepMesh_BoundaryParamsRangeSplitter` returns false.
    fn seeds_cn_intervals(&self) -> bool {
        true
    }
}

/// Base range splitter. Source: `BRepMesh_DefaultRangeSplitter`.
pub struct DefaultRangeSplitter {
    pub(in crate::meshing::range_splitter) dface: Option<MeshFace>,
    pub(in crate::meshing::range_splitter) surface: Option<Arc<dyn Surface>>,
    pub(in crate::meshing::range_splitter) deflection: f64,
    pub(in crate::meshing::range_splitter) range_u: (f64, f64),
    pub(in crate::meshing::range_splitter) range_v: (f64, f64),
    pub(in crate::meshing::range_splitter) delta: (f64, f64),
    pub(in crate::meshing::range_splitter) tolerance: (f64, f64),
    pub(in crate::meshing::range_splitter) is_valid: bool,
}
