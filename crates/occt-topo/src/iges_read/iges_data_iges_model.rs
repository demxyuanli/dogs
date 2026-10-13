//! `IGESData_IGESModel` directory and own parameters for the planar family.
//!
//! Field order follows each tool's `ReadOwnParams`.

use super::iges_data_param_reader::{IgesParam, ParamReader};

/// `IGESData_BasicEditor::UnitFlagValue` (`IGESData_BasicEditor.cxx:492-520`)
/// divided by the default cascade unit 1 (`IGESData_GlobalSection::UnitValue`).
pub(in crate::iges_read) fn unit_flag_value(flag: i32) -> f64 {
    match flag {
        1 => 25.4,
        2 => 1.0,
        3 => 1.0,
        4 => 304.8,
        5 => 1_609_344.0,
        6 => 1000.0,
        7 => 1_000_000.0,
        8 => 0.0254,
        9 => 0.001,
        10 => 10.0,
        11 => 0.0000254,
        _ => 0.0,
    }
}

/// Directory sequence (1-based, odd for the first line of a pair) to the
/// 1-based entity rank. `IGESData_IGESReaderTool` uses `(1 + pointer) / 2`.
pub(in crate::iges_read) fn de_pointer_to_entity(de: i32) -> Result<Option<usize>, String> {
    if de == 0 {
        return Ok(None);
    }
    if de < 0 {
        return Err(format!(
            "negative IGES DE pointer {de} is not ported (IGESData_IGESReaderTool.cxx)"
        ));
    }
    Ok(Some(((de as usize) + 1) / 2))
}

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
        edges: Vec<super::igessolid::SolidEdgeItem>,
    },
    /// 508. `IGESSolid_ToolLoop::ReadOwnParams` (`cxx:53-232`).
    SolidLoop {
        edges: Vec<super::igessolid::SolidLoopItem>,
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

#[derive(Clone, Debug)]
pub(in crate::iges_read) struct IgesEntity {
    pub index: usize,
    pub type_number: i32,
    pub form_number: i32,
    /// `IGESData_IGESEntity::BlankStatus` (`cxx:256-259`), low 4 bits of the
    /// two-digit blank field. `00` visible, `01` blanked.
    pub blank_status: i32,
    /// Directory field 7. Zero means no 124 matrix.
    pub transform: i32,
    pub data: EntityData,
    pub shared: bool,
}

#[derive(Clone, Debug)]
pub(in crate::iges_read) struct IgesModel {
    pub entities: Vec<IgesEntity>,
    /// `IGESData_GlobalSection::UnitValue`.
    pub unit_factor: f64,
    /// Global parameter 19, minimum resolution.
    pub resolution: f64,
}

impl IgesModel {
    pub(in crate::iges_read) fn entity(&self, index: usize) -> Result<&IgesEntity, String> {
        self.entities
            .get(index)
            .filter(|e| e.index == index)
            .ok_or_else(|| format!("IGES entity {index} is out of range"))
    }
}

/// One model-space curve of entity 141, with its sense and parameter curves.
#[derive(Clone, Debug)]
pub(in crate::iges_read) struct BoundaryCurve {
    pub model: usize,
    pub sense: i32,
    pub parameters: Vec<usize>,
}

pub(in crate::iges_read) struct DirLine {
    pub type_number: i32,
    pub form_number: i32,
    pub blank_status: i32,
    pub transform: i32,
    pub de_seq: i32,
}

/// Directory fields filled by `iges_Dsect` (`analiges.c:75-112`).

fn parse_boundary(pr: &mut ParamReader<'_>) -> Result<EntityData, String> {
    // `IGESGeom_ToolBoundary.cxx:66-212`.
    let boundary_type = pr.read_integer()?;
    let preference = pr.read_integer()?;
    let surface = pr.read_entity(true)?;
    let num = pr.read_integer()?;
    if num <= 0 {
        return Err("XTSEP_126: IGESGeom_ToolBoundary.cxx:106-108".into());
    }
    let mut curves = Vec::with_capacity(num as usize);
    for _ in 0..num {
        let model = pr
            .read_entity(false)?
            .ok_or_else(|| "IGES 141 null model curve".to_string())?;
        let sense = pr.read_integer()?;
        let count = pr.read_integer()?;
        if count < 0 {
            return Err("XTSEP_129: IGESGeom_ToolBoundary.cxx:202-205".into());
        }
        let mut parameters = Vec::with_capacity(count as usize);
        for _ in 0..count {
            if let Some(curve) = pr.read_entity(true)? {
                parameters.push(curve);
            }
        }
        curves.push(BoundaryCurve {
            model,
            sense,
            parameters,
        });
    }
    Ok(EntityData::Boundary {
        boundary_type,
        preference,
        surface,
        curves,
    })
}

fn parse_bounded(pr: &mut ParamReader<'_>) -> Result<EntityData, String> {
    // `IGESGeom_ToolBoundedSurface.cxx:60-98`.
    let representation = pr.read_integer()?;
    let surface = pr
        .read_entity(false)?
        .ok_or_else(|| "IGES 143 null basis surface".to_string())?;
    let num = pr.read_integer()?;
    if num <= 0 {
        return Err("XTSEP_167: IGESGeom_ToolBoundedSurface.cxx:100-103".into());
    }
    let mut boundaries = Vec::with_capacity(num as usize);
    for _ in 0..num {
        let boundary = pr
            .read_entity(false)?
            .ok_or_else(|| "IGES 143 null boundary".to_string())?;
        boundaries.push(boundary);
    }
    Ok(EntityData::Bounded {
        representation,
        surface,
        boundaries,
    })
}

/// `IGESBasic_ToolSingleParent::ReadOwnParams` (`cxx:42-105`), the 402 form 9
/// layout: parent count, child count, the single parent, then the children.
/// `ReadEnts` only runs for `nbval > 0` (`cxx:98-100`), so a non-positive child
/// count leaves the child array null and `NbChildren()` at 0. A null pointer is
/// kept as `None` (`ReadEntity` / `ReadEnts` report it and continue,
/// `cxx:75-96`).
fn parse_single_parent(pr: &mut ParamReader<'_>) -> Result<EntityData, String> {
    let _nb_parents = pr.read_integer()?;
    // `OwnCheck` (`cxx:186-190`) only reports `XSTEP_204` when the parent count
    // is not 1 and transfers anyway; the port has no check stage.
    let nb = pr.read_integer()?;
    let parent = pr.read_entity(true)?;
    let mut children = Vec::new();
    if nb > 0 {
        children.reserve(nb as usize);
        for _ in 0..nb {
            children.push(pr.read_entity(true)?);
        }
    }
    Ok(EntityData::SingleParent { parent, children })
}

fn parse_own(type_number: i32, form_number: i32, tokens: &[IgesParam]) -> Result<EntityData, String> {
    if tokens.is_empty() {
        return Err(format!("IGES entity {type_number} has no parameters"));
    }
    match &tokens[0] {
        IgesParam::Int(n) if *n == type_number => {}
        IgesParam::Real(r) if *r == type_number as f64 => {}
        _ => {
            return Err(format!(
                "IGES parameter type does not match directory type {type_number}"
            ));
        }
    }
    let mut pr = ParamReader::new(tokens);
    match type_number {
        402 => {
            // `IGESBasic_ReadWriteModule::CaseIGES` (`cxx:60-72`): form 9 is
            // `IGESBasic_SingleParent`, forms 1/7/14/15 are the `IGESBasic_Group`
            // family (one entity list each).
            if form_number == 9 {
                return parse_single_parent(&mut pr);
            }
            let nb = pr.read_integer()?;
            if nb < 0 {
                return Err("IGES 402 entity count is negative".into());
            }
            let mut entities = Vec::with_capacity(nb as usize);
            for _ in 0..nb {
                let ix = pr.read_entity(false)?;
                entities.push(ix.ok_or_else(|| "IGES 402 null entity".to_string())?);
            }
            Ok(EntityData::Group { entities })
        }
        308 => super::iges_to_brep_curve_and_surface::parse_subfigure_def(&mut pr),
        408 => super::iges_to_brep_curve_and_surface::parse_singular_subfigure(&mut pr),
        144 => {
            let surface = pr
                .read_entity(false)?
                .ok_or_else(|| "IGES 144 null surface".to_string())?;
            let _flag = pr.read_integer()?;
            let count = pr.read_integer()?;
            if count < 0 {
                return Err("IGES 144 inner contour count is negative".into());
            }
            // Outer is read before the inner list (`IGESGeom_ToolTrimmedSurface.cxx:111-149`).
            let outer = pr.read_entity(true)?;
            let mut inners = Vec::with_capacity(count as usize);
            for _ in 0..count {
                inners.push(pr.read_entity(false)?);
            }
            Ok(EntityData::Trimmed {
                surface,
                outer,
                inners,
            })
        }
        108 => {
            let a = pr.read_real()?;
            let b = pr.read_real()?;
            let c = pr.read_real()?;
            let d = pr.read_real()?;
            let bounding = if pr.is_defined() {
                pr.read_entity(true)?
            } else {
                None
            };
            // Display symbol is not used by `TransferPlane`. Consume it when present.
            if pr.is_defined() {
                let _attach = pr.read_xyz()?;
                if pr.is_defined() {
                    let _size = pr.read_real()?;
                }
            }
            Ok(EntityData::Plane {
                a,
                b,
                c,
                d,
                bounding,
            })
        }
        142 => {
            let mode = pr.read_integer()?;
            let surface = pr.read_entity(false)?;
            let curve_uv = pr.read_entity(true)?;
            let curve_3d = pr.read_entity(true)?;
            let preference = pr.read_integer()?;
            Ok(EntityData::CurveOnSurface {
                mode,
                surface,
                curve_uv,
                curve_3d,
                preference,
            })
        }
        102 => {
            let num = pr.read_integer()?;
            if num <= 0 {
                return Err("IGES 102 curve count is not positive".into());
            }
            let mut curves = Vec::with_capacity(num as usize);
            for _ in 0..num {
                let ix = pr.read_entity(false)?;
                curves.push(ix.ok_or_else(|| "IGES 102 null curve".to_string())?);
            }
            Ok(EntityData::Composite { curves })
        }
        110 => {
            let start = pr.read_xyz()?;
            let end = pr.read_xyz()?;
            Ok(EntityData::Line { start, end })
        }
        106 => {
            // `IGESGeom_ToolCopiousData.cxx:58-104`.
            let data_type = pr.read_integer()?;
            let nb = pr.read_integer()?;
            if nb <= 0 {
                return Err("IGES 106 tuple count is not positive".into());
            }
            let z_plane = if data_type == 1 { pr.read_real()? } else { 0.0 };
            let upper = if data_type == 1 {
                2 * nb
            } else if data_type == 2 {
                3 * nb
            } else {
                6 * nb
            };
            let mut tuples = Vec::with_capacity(upper as usize);
            for _ in 0..upper {
                tuples.push(pr.read_real()?);
            }
            Ok(EntityData::Copious {
                data_type,
                z_plane,
                tuples,
            })
        }
        112 => super::spline_curve::parse_spline_curve(&mut pr),
        114 => super::spline_surface::parse_spline_surface(&mut pr),
        104 => {
            // `IGESGeom_ToolConicArc.cxx:70-84`.
            let a = pr.read_real()?;
            let b = pr.read_real()?;
            let c = pr.read_real()?;
            let d = pr.read_real()?;
            let e = pr.read_real()?;
            let f = pr.read_real()?;
            let zt = pr.read_real()?;
            let start = pr.read_xy()?;
            let end = pr.read_xy()?;
            Ok(EntityData::ConicArc {
                a,
                b,
                c,
                d,
                e,
                f,
                zt,
                start,
                end,
            })
        }
        100 => {
            // `IGESGeom_ToolCircularArc.cxx:59-66`.
            let zt = pr.read_real()?;
            let center = pr.read_xy()?;
            let start = pr.read_xy()?;
            let end = pr.read_xy()?;
            Ok(EntityData::CircularArc {
                zt,
                center,
                start,
                end,
            })
        }
        118 => {
            // `IGESGeom_ToolRuledSurface.cxx:48-107`.
            let first = pr
                .read_entity(false)?
                .ok_or_else(|| "IGES 118 null first curve".to_string())?;
            let second = pr
                .read_entity(false)?
                .ok_or_else(|| "IGES 118 null second curve".to_string())?;
            let direction = pr.read_integer()?;
            let dev_flag = pr.read_integer()?;
            Ok(EntityData::RuledSurface {
                first,
                second,
                direction,
                dev_flag,
            })
        }
        122 => {
            // `IGESGeom_ToolTabulatedCylinder.cxx:59-83`.
            let directrix = pr
                .read_entity(false)?
                .ok_or_else(|| "IGES 122 null directrix".to_string())?;
            let end = pr.read_xyz()?;
            Ok(EntityData::Tabulated { directrix, end })
        }
        141 => parse_boundary(&mut pr),
        143 => parse_bounded(&mut pr),
        140 => {
            // `IGESGeom_ToolOffsetSurface.cxx:57-98`.
            let indicator = pr.read_xyz()?;
            let distance = pr.read_real()?;
            let surface = pr
                .read_entity(false)?
                .ok_or_else(|| "IGES 140 null basis surface".to_string())?;
            Ok(EntityData::OffsetSurface {
                indicator,
                distance,
                surface,
            })
        }
        120 => {
            // `IGESGeom_ToolSurfaceOfRevolution.cxx:51-108`.
            let axis = pr
                .read_entity(false)?
                .ok_or_else(|| "IGES 120 null axis".to_string())?;
            let generatrix = pr
                .read_entity(false)?
                .ok_or_else(|| "IGES 120 null generatrix".to_string())?;
            let start_angle = pr.read_real()?;
            let end_angle = pr.read_real()?;
            Ok(EntityData::Revolution {
                axis,
                generatrix,
                start_angle,
                end_angle,
            })
        }
        124 => {
            // `IGESGeom_ToolTransformationMatrix.cxx:53-68`, row-major 3x4.
            let mut rows = [[0.0; 4]; 3];
            for row in &mut rows {
                for cell in row.iter_mut() {
                    *cell = pr.read_real()?;
                }
            }
            Ok(EntityData::TransformMatrix { rows })
        }
        126 => parse_bspline(&mut pr),
        116 => {
            let xyz = pr.read_xyz()?;
            if pr.defined_else_skip() {
                let _symbol = pr.read_entity(true)?;
            }
            Ok(EntityData::Point { xyz })
        }
        130 => super::offset_curve::parse_offset_curve(&mut pr),
        123 => {
            let xy = pr.read_xy()?;
            let z = if pr.defined_else_skip() {
                pr.read_real()?
            } else {
                0.0
            };
            Ok(EntityData::Direction {
                xyz: [xy[0], xy[1], z],
            })
        }
        128 => super::iges_surface_params::parse_bspline_surface(&mut pr),
        190 => super::iges_surface_params::parse_plane_surface(&mut pr, form_number),
        192 => super::iges_surface_params::parse_cylindrical_surface(&mut pr, form_number),
        194 => super::iges_surface_params::parse_conical_surface(&mut pr, form_number),
        196 => super::iges_surface_params::parse_spherical_surface(&mut pr, form_number),
        198 => super::iges_surface_params::parse_toroidal_surface(&mut pr, form_number),
        502 => super::igessolid::parse_vertex_list(&mut pr),
        504 => super::igessolid::parse_edge_list(&mut pr),
        508 => super::igessolid::parse_loop(&mut pr),
        510 => super::igessolid::parse_face(&mut pr),
        514 => super::igessolid::parse_shell(&mut pr),
        186 => super::igessolid::parse_manifold_solid(&mut pr),
        _ => Ok(EntityData::Raw { type_number }),
    }
}

fn parse_bspline(pr: &mut ParamReader<'_>) -> Result<EntityData, String> {
    // `IGESGeom_ToolBSplineCurve.cxx:83-161`.
    let mut upper = pr.read_integer()?;
    if upper < 0 {
        upper = 0;
    }
    let degree = match pr.read_integer() {
        Ok(n) => n,
        Err(_) => 0,
    };
    let planar = pr.read_integer()? != 0;
    let closed = pr.read_integer()? != 0;
    let polynomial = pr.read_integer()? != 0;
    let periodic = pr.read_integer()? != 0;
    let nb_knots = upper + degree + 2;
    if nb_knots < 0 {
        return Err("IGES 126 knot count is negative".into());
    }
    let mut knots = Vec::with_capacity(nb_knots as usize);
    for _ in 0..nb_knots {
        knots.push(pr.read_real()?);
    }
    let n_poles = upper + 1;
    if n_poles < 0 {
        return Err("IGES 126 pole count is negative".into());
    }
    let mut weights = Vec::with_capacity(n_poles as usize);
    for _ in 0..n_poles {
        weights.push(pr.read_real()?);
    }
    let mut poles = Vec::with_capacity(n_poles as usize);
    for _ in 0..n_poles {
        poles.push(pr.read_xyz()?);
    }
    let u_min = pr.read_real()?;
    let u_max = pr.read_real()?;
    // Unit normal is optional (`cxx:167-212`) and unused by TransferBSplineCurve.
    if pr.defined_else_skip() {
        let _nx = pr.read_real()?;
    }
    if pr.defined_else_skip() {
        let _ny = pr.read_real()?;
    }
    if pr.defined_else_skip() {
        let _nz = pr.read_real()?;
    }
    Ok(EntityData::BSplineCurve {
        upper,
        degree,
        planar,
        closed,
        polynomial,
        periodic,
        knots,
        weights,
        poles,
        u_min,
        u_max,
    })
}

fn pointer_targets(data: &EntityData) -> Vec<usize> {
    let mut v = Vec::new();
    match data {
        EntityData::Group { entities } => v.extend(entities.iter().copied()),
        EntityData::SingleParent { parent, children } => {
            if let Some(p) = parent {
                v.push(*p);
            }
            for c in children.iter().flatten() {
                v.push(*c);
            }
        }
        EntityData::SubfigureDef { entities, .. } => {
            for entity in entities {
                if let Some(index) = entity {
                    v.push(*index);
                }
            }
        }
        EntityData::SingularSubfigure { subfigure, .. } => v.push(*subfigure),
        EntityData::Trimmed {
            surface,
            outer,
            inners,
        } => {
            v.push(*surface);
            if let Some(o) = outer {
                v.push(*o);
            }
            for i in inners.iter().flatten() {
                v.push(*i);
            }
        }
        EntityData::Plane { bounding, .. } => {
            if let Some(b) = bounding {
                v.push(*b);
            }
        }
        EntityData::CurveOnSurface {
            surface,
            curve_uv,
            curve_3d,
            ..
        } => {
            if let Some(s) = surface {
                v.push(*s);
            }
            if let Some(u) = curve_uv {
                v.push(*u);
            }
            if let Some(c) = curve_3d {
                v.push(*c);
            }
        }
        EntityData::Composite { curves } => v.extend(curves.iter().copied()),
        EntityData::Tabulated { directrix, .. } => v.push(*directrix),
        EntityData::OffsetSurface { surface, .. } => v.push(*surface),
        EntityData::Bounded {
            surface, boundaries, ..
        } => {
            v.push(*surface);
            v.extend(boundaries.iter().copied());
        }
        EntityData::Boundary {
            surface, curves, ..
        } => {
            if let Some(s) = surface {
                v.push(*s);
            }
            for curve in curves {
                v.push(curve.model);
                v.extend(curve.parameters.iter().copied());
            }
        }
        EntityData::RuledSurface { first, second, .. } => {
            v.push(*first);
            v.push(*second);
        }
        EntityData::Revolution { axis, generatrix, .. } => {
            v.push(*axis);
            v.push(*generatrix);
        }
        EntityData::PlaneSurface {
            location,
            normal,
            reference,
        } => {
            v.push(*location);
            v.push(*normal);
            if let Some(r) = reference {
                v.push(*r);
            }
        }
        EntityData::CylindricalSurface {
            location,
            axis,
            reference,
            ..
        }
        | EntityData::ConicalSurface {
            location,
            axis,
            reference,
            ..
        } => {
            v.push(*location);
            v.push(*axis);
            if let Some(r) = reference {
                v.push(*r);
            }
        }
        EntityData::SphericalSurface {
            center,
            axis,
            reference,
            ..
        } => {
            v.push(*center);
            if let Some(a) = axis {
                v.push(*a);
            }
            if let Some(r) = reference {
                v.push(*r);
            }
        }
        EntityData::ToroidalSurface {
            center,
            axis,
            reference,
            ..
        } => {
            v.push(*center);
            v.push(*axis);
            if let Some(r) = reference {
                v.push(*r);
            }
        }
        EntityData::Line { .. }
        | EntityData::CircularArc { .. }
        | EntityData::ConicArc { .. }
        | EntityData::Copious { .. }
        | EntityData::SplineCurve { .. }
        | EntityData::SplineSurface { .. }
        | EntityData::TransformMatrix { .. }
        | EntityData::BSplineCurve { .. }
        | EntityData::BSplineSurface { .. } => {}
        EntityData::OffsetCurve { base, function, .. } => {
            v.push(*base);
            if let Some(f) = function {
                v.push(*f);
            }
        }
        EntityData::Point { .. }
        | EntityData::Direction { .. }
        | EntityData::Raw { .. } => {}
        // `IGESSolid_*::OwnShared` (`IGESSolid_ToolVertexList.cxx:...`,
        // `EdgeList.cxx:...`, `Loop.cxx:...`, `Face.cxx:...`, `Shell.cxx:...`,
        // `ManifoldSolid.cxx:165-179`).
        EntityData::VertexList { .. } => {}
        EntityData::EdgeList { edges } => {
            for edge in edges {
                if let Some(c) = edge.curvature {
                    v.push(c);
                }
                if let Some(l) = edge.start_vertex_list {
                    v.push(l);
                }
                if let Some(l) = edge.end_vertex_list {
                    v.push(l);
                }
            }
        }
        EntityData::SolidLoop { edges } => {
            for edge in edges {
                if let Some(e) = edge.edge {
                    v.push(e);
                }
                for c in edge.param_curves.iter().flatten() {
                    v.push(*c);
                }
            }
        }
        EntityData::SolidFace { surface, loops, .. } => {
            if let Some(s) = surface {
                v.push(*s);
            }
            for l in loops.iter().flatten() {
                v.push(*l);
            }
        }
        EntityData::SolidShell { faces, .. } => {
            for f in faces.iter().flatten() {
                v.push(*f);
            }
        }
        EntityData::ManifoldSolid {
            shell, void_shells, ..
        } => {
            if let Some(s) = shell {
                v.push(*s);
            }
            for s in void_shells.iter().flatten() {
                v.push(*s);
            }
        }
    }
    v
}

pub(in crate::iges_read) fn build_model(
    global: &[IgesParam],
    directory: &[DirLine],
    params_by_de: &std::collections::HashMap<i32, Vec<IgesParam>>,
) -> Result<IgesModel, String> {
    // Global parameters are 1-based and include the two delimiter fields.
    // Parameter 14 is the unit flag, 19 is the minimum resolution.
    let unit_flag = match global.get(13) {
        Some(IgesParam::Int(n)) => *n,
        Some(IgesParam::Real(r)) => *r as i32,
        _ => 2,
    };
    let resolution = match global.get(18) {
        Some(IgesParam::Real(r)) => *r,
        Some(IgesParam::Int(n)) => *n as f64,
        _ => 1.0e-7,
    };
    let unit_factor = unit_flag_value(unit_flag);
    if unit_factor == 0.0 {
        return Err(format!("IGES unit flag {unit_flag} has no UnitFlagValue"));
    }

    let mut entities = Vec::with_capacity(directory.len() + 1);
    entities.push(IgesEntity {
        index: 0,
        type_number: 0,
        form_number: 0,
        blank_status: 0,
        transform: 0,
        data: EntityData::Raw { type_number: 0 },
        shared: false,
    });
    for (i, dir) in directory.iter().enumerate() {
        let index = i + 1;
        let tokens = params_by_de.get(&dir.de_seq).ok_or_else(|| {
            format!(
                "IGES entity {index} (DE {}) has no parameter data",
                dir.de_seq
            )
        })?;
        let data = parse_own(dir.type_number, dir.form_number, tokens)?;
        entities.push(IgesEntity {
            index,
            type_number: dir.type_number,
            form_number: dir.form_number,
            blank_status: dir.blank_status,
            transform: dir.transform,
            data,
            shared: false,
        });
    }

    let mut shared = vec![false; entities.len()];
    for ent in entities.iter().skip(1) {
        for p in pointer_targets(&ent.data) {
            if p < shared.len() {
                shared[p] = true;
            } else {
                return Err(format!("IGES pointer to entity {p} is out of range"));
            }
        }
    }
    for (i, flag) in shared.iter().enumerate().skip(1) {
        entities[i].shared = *flag;
    }

    Ok(IgesModel {
        entities,
        unit_factor,
        resolution,
    })
}
