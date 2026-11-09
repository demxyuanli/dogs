use crate::iges_read::iges_data_param_reader::{IgesParam, ParamReader};

use super::entities::{BoundaryCurve, EntityData};

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

pub(super) fn parse_own(
    type_number: i32,
    form_number: i32,
    tokens: &[IgesParam],
) -> Result<EntityData, String> {
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
        308 => crate::iges_read::iges_to_brep_curve_and_surface::parse_subfigure_def(&mut pr),
        408 => crate::iges_read::iges_to_brep_curve_and_surface::parse_singular_subfigure(&mut pr),
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
        112 => crate::iges_read::spline_curve::parse_spline_curve(&mut pr),
        114 => crate::iges_read::spline_surface::parse_spline_surface(&mut pr),
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
        130 => crate::iges_read::offset_curve::parse_offset_curve(&mut pr),
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
        128 => crate::iges_read::iges_surface_params::parse_bspline_surface(&mut pr),
        190 => crate::iges_read::iges_surface_params::parse_plane_surface(&mut pr, form_number),
        192 => crate::iges_read::iges_surface_params::parse_cylindrical_surface(&mut pr, form_number),
        194 => crate::iges_read::iges_surface_params::parse_conical_surface(&mut pr, form_number),
        196 => crate::iges_read::iges_surface_params::parse_spherical_surface(&mut pr, form_number),
        198 => crate::iges_read::iges_surface_params::parse_toroidal_surface(&mut pr, form_number),
        502 => crate::iges_read::igessolid::parse_vertex_list(&mut pr),
        504 => crate::iges_read::igessolid::parse_edge_list(&mut pr),
        508 => crate::iges_read::igessolid::parse_loop(&mut pr),
        510 => crate::iges_read::igessolid::parse_face(&mut pr),
        514 => crate::iges_read::igessolid::parse_shell(&mut pr),
        186 => crate::iges_read::igessolid::parse_manifold_solid(&mut pr),
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
