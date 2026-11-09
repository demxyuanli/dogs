#[derive(Clone, Debug)]
pub(in crate::iges_read) enum EntityData {
    /// 402 form 1. `IGESBasic_ToolGroup::ReadOwnParams` (`cxx:41-74`).
    Group { entities: Vec<usize> },
    /// 402 form 9. `IGESBasic_ToolSingleParent::ReadOwnParams` (`cxx:42-105`),
    /// selected for `CaseIGES(402, 9) = 14`
    /// (`IGESBasic_ReadWriteModule.cxx:60-72`, `:263-270`). `NbParentEntities`
    /// is consumed but not stored: `OwnCheck`/`OwnCorrect`
    /// (`IGESBasic_ToolSingleParent.cxx:148-162`, `:177-191`) only report and
    /// rewrite it and the port has no check/correct stage.
    SingleParent {
        parent: Option<usize>,
        children: Vec<Option<usize>>,
    },
    /// 308. `IGESBasic_ToolSubfigureDef::ReadOwnParams` (`cxx:42-79`).
    /// Depth and the name are stored by `Init` and are not used by transfer.
    SubfigureDef {
        #[allow(dead_code)]
        depth: i32,
        #[allow(dead_code)]
        name: String,
        entities: Vec<Option<usize>>,
    },
    /// 408. `IGESBasic_ToolSingularSubfigure::ReadOwnParams` (`cxx:40-112`).
    SingularSubfigure {
        subfigure: usize,
        translation: [f64; 3],
        has_scale: bool,
        scale: f64,
    },
    /// 144. `IGESGeom_ToolTrimmedSurface::ReadOwnParams` (`cxx:42-191`).
    Trimmed {
        surface: usize,
        outer: Option<usize>,
        inners: Vec<Option<usize>>,
    },
    /// 108. `IGESGeom_ToolPlane::ReadOwnParams` (`cxx:48-124`).
    /// IGES equation is `ax+by+cz = d`.
    Plane {
        a: f64,
        b: f64,
        c: f64,
        d: f64,
        bounding: Option<usize>,
    },
    /// 142. `IGESGeom_ToolCurveOnSurface::ReadOwnParams` (`cxx:41-139`).
    CurveOnSurface {
        mode: i32,
        surface: Option<usize>,
        curve_uv: Option<usize>,
        curve_3d: Option<usize>,
        preference: i32,
    },
    /// 102. `IGESGeom_ToolCompositeCurve::ReadOwnParams` (`cxx:38-67`).
    Composite { curves: Vec<usize> },
    /// 110. `IGESGeom_ToolLine::ReadOwnParams` (`cxx:39-61`).
    Line { start: [f64; 3], end: [f64; 3] },
    /// 100. `IGESGeom_ToolCircularArc::ReadOwnParams` (`cxx:42-77`).
    /// `ZT`, center XY, start XY, end XY.
    CircularArc {
        zt: f64,
        center: [f64; 2],
        start: [f64; 2],
        end: [f64; 2],
    },
    /// 104. `IGESGeom_ToolConicArc::ReadOwnParams` (`cxx:48-98`).
    /// `A..F`, `ZT`, start XY, end XY. Form 0 is recomputed by `ComputedFormNumber`.
    /// 106. `IGESGeom_ToolCopiousData::ReadOwnParams` (`cxx:43-109`).
    /// Data type 1 is XY plus `ZPlane`, type 2 is XYZ, type 3 is XYZ plus a vector.
    Copious {
        data_type: i32,
        z_plane: f64,
        tuples: Vec<f64>,
    },
    /// 112. `IGESGeom_ToolSplineCurve::ReadOwnParams` (`cxx:42-180`).
    /// `continuity_degree` and the terminate values are stored by `Init` and
    /// are not read by `SplineCurveFromIGES`.
    SplineCurve {
        spline_type: i32,
        #[allow(dead_code)]
        continuity_degree: i32,
        nb_dimensions: i32,
        breaks: Vec<f64>,
        x_poly: Vec<[f64; 4]>,
        y_poly: Vec<[f64; 4]>,
        z_poly: Vec<[f64; 4]>,
        #[allow(dead_code)]
        x_values: [f64; 4],
        #[allow(dead_code)]
        y_values: [f64; 4],
        #[allow(dead_code)]
        z_values: [f64; 4],
    },
    /// 114. `IGESGeom_ToolSplineSurface::ReadOwnParams` (`cxx:42-227`).
    /// `patch_type` is stored by `Init` and is not read by `SplineSurfaceFromIGES`.
    SplineSurface {
        boundary_type: i32,
        #[allow(dead_code)]
        patch_type: i32,
        u_breaks: Vec<f64>,
        v_breaks: Vec<f64>,
        x_poly: Vec<Vec<[f64; 16]>>,
        y_poly: Vec<Vec<[f64; 16]>>,
        z_poly: Vec<Vec<[f64; 16]>>,
    },
    ConicArc {
        a: f64,
        b: f64,
        c: f64,
        d: f64,
        e: f64,
        f: f64,
        zt: f64,
        start: [f64; 2],
        end: [f64; 2],
    },
    /// 118. `IGESGeom_ToolRuledSurface::ReadOwnParams` (`cxx:48-107`).
    /// `IsRuledByParameter` is `FormNumber() != 0` (`IGESGeom_RuledSurface.cxx:64-67`).
    RuledSurface {
        first: usize,
        second: usize,
        direction: i32,
        #[allow(dead_code)]
        dev_flag: i32,
    },
    /// 122. `IGESGeom_ToolTabulatedCylinder::ReadOwnParams` (`cxx:52-87`).
    Tabulated {
        directrix: usize,
        end: [f64; 3],
    },
    /// 141. `IGESGeom_ToolBoundary::ReadOwnParams` (`cxx:50-212`).
    Boundary {
        boundary_type: i32,
        preference: i32,
        surface: Option<usize>,
        curves: Vec<BoundaryCurve>,
    },
    /// 143. `IGESGeom_ToolBoundedSurface::ReadOwnParams` (`cxx:46-98`).
    Bounded {
        representation: i32,
        surface: usize,
        boundaries: Vec<usize>,
    },
    /// 140. `IGESGeom_ToolOffsetSurface::ReadOwnParams` (`cxx:51-98`).
    /// Transfer uses `Distance` and `Surface` only (`TopoSurface.cxx:1059`, `1121`).
    OffsetSurface {
        #[allow(dead_code)]
        indicator: [f64; 3],
        distance: f64,
        surface: usize,
    },
    /// 120. `IGESGeom_ToolSurfaceOfRevolution::ReadOwnParams` (`cxx:38-117`).
    Revolution {
        axis: usize,
        generatrix: usize,
        start_angle: f64,
        end_angle: f64,
    },
    /// 124. `IGESGeom_ToolTransformationMatrix::ReadOwnParams` (`cxx:42-73`).
    /// 3x4 row-major.
    TransformMatrix { rows: [[f64; 4]; 3] },
    /// 116. `IGESGeom_ToolPoint::ReadOwnParams` (`cxx:42-90`). The display symbol is skipped.
    Point { xyz: [f64; 3] },
    /// 130. `IGESGeom_ToolOffsetCurve::ReadOwnParams` (`cxx:47-209`).
    /// `TransferOffsetCurve` uses the type, the first distance, the normal, and
    /// the parameter window. The function, taper, and second distance are stored.
    OffsetCurve {
        base: usize,
        offset_type: i32,
        #[allow(dead_code)]
        function: Option<usize>,
        #[allow(dead_code)]
        function_coord: i32,
        #[allow(dead_code)]
        tapered_type: i32,
        distance1: f64,
        #[allow(dead_code)]
        arc_length1: f64,
        #[allow(dead_code)]
        distance2: f64,
        #[allow(dead_code)]
        arc_length2: f64,
        normal: [f64; 3],
        start_param: f64,
        end_param: f64,
    },
    /// 123. `IGESGeom_ToolDirection::ReadOwnParams` (`cxx:38-68`).
    Direction { xyz: [f64; 3] },
    /// 128. `IGESGeom_ToolBSplineSurface::ReadOwnParams` (`cxx:46-268`).
    /// Parameter bounds are read as Umin, Vmin, Umax, Vmax (`cxx:196-205`).
    BSplineSurface {
        degree_u: i32,
        degree_v: i32,
        closed_u: bool,
        closed_v: bool,
        polynomial: bool,
        periodic_u: bool,
        periodic_v: bool,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
        /// Row-major in U, then V. File order is V-major (`cxx:144-187`).
        weights: Vec<Vec<f64>>,
        poles: Vec<Vec<[f64; 3]>>,
        #[allow(dead_code)]
        u_min: f64,
        #[allow(dead_code)]
        u_max: f64,
        #[allow(dead_code)]
        v_min: f64,
        #[allow(dead_code)]
        v_max: f64,
    },
    /// 190. `IGESSolid_ToolPlaneSurface::ReadOwnParams` (`cxx:43-151`).
    PlaneSurface {
        location: usize,
        normal: usize,
        reference: Option<usize>,
    },
    /// 192. `IGESSolid_ToolCylindricalSurface::ReadOwnParams` (`cxx:36-71`).
    CylindricalSurface {
        location: usize,
        axis: usize,
        radius: f64,
        reference: Option<usize>,
    },
    /// 194. `IGESSolid_ToolConicalSurface::ReadOwnParams` (`cxx:36-72`).
    ConicalSurface {
        location: usize,
        axis: usize,
        radius: f64,
        semi_angle: f64,
        reference: Option<usize>,
    },
    /// 196. `IGESSolid_ToolSphericalSurface::ReadOwnParams` (`cxx:36-71`).
    SphericalSurface {
        center: usize,
        radius: f64,
        axis: Option<usize>,
        reference: Option<usize>,
    },
    /// 198. `IGESSolid_ToolToroidalSurface::ReadOwnParams` (`cxx:36-68`).
    ToroidalSurface {
        center: usize,
        axis: usize,
        major: f64,
        minor: f64,
        reference: Option<usize>,
    },
    /// 126. `IGESGeom_ToolBSplineCurve::ReadOwnParams` (`cxx:45-227`).
    BSplineCurve {
        #[allow(dead_code)]
        upper: i32,
        degree: i32,
        #[allow(dead_code)]
        planar: bool,
        #[allow(dead_code)]
        closed: bool,
        polynomial: bool,
        #[allow(dead_code)]
        periodic: bool,
        knots: Vec<f64>,
        weights: Vec<f64>,
        poles: Vec<[f64; 3]>,
        u_min: f64,
        u_max: f64,
    },
    /// Types outside the transferred family. Transfer cites the OCCT arm.
    Raw { type_number: i32 },
    /// 502. `IGESSolid_ToolVertexList::ReadOwnParams` (`cxx:46-91`).
    VertexList { vertices: Vec<[f64; 3]> },
    /// 504. `IGESSolid_ToolEdgeList::ReadOwnParams` (`cxx:47-297`).
    EdgeList {
        edges: Vec<super::super::igessolid::SolidEdgeItem>,
    },
    /// 508. `IGESSolid_ToolLoop::ReadOwnParams` (`cxx:53-232`).
    SolidLoop {
        edges: Vec<super::super::igessolid::SolidLoopItem>,
    },
    /// 510. `IGESSolid_ToolFace::ReadOwnParams` (`cxx:47-140`).
    SolidFace {
        surface: Option<usize>,
        loops: Vec<Option<usize>>,
        /// `IGESSolid_Face::OuterLoopFlag` is stored by `Init`; the transfer
        /// picks the outer loop from `Loop.Orientation`, not from this.
        #[allow(dead_code)]
        outer_loop: bool,
    },
    /// 514. `IGESSolid_ToolShell::ReadOwnParams` (`cxx:47-124`).
    SolidShell {
        faces: Vec<Option<usize>>,
        orientations: Vec<bool>,
    },
    /// 186. `IGESSolid_ToolManifoldSolid::ReadOwnParams` (`cxx:43-149`).
    ManifoldSolid {
        shell: Option<usize>,
        oriented: bool,
        void_shells: Vec<Option<usize>>,
        /// `VoidOrientationFlag` is read by the tool but
        /// `TransferManifoldSolid` reuses the single `isoriented` flag.
        #[allow(dead_code)]
        void_flags: Vec<bool>,
    },
}

/// One model-space curve of entity 141, with its sense and parameter curves.
#[derive(Clone, Debug)]
pub(in crate::iges_read) struct BoundaryCurve {
    pub model: usize,
    pub sense: i32,
    pub parameters: Vec<usize>,
}
