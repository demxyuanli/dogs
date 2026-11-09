//! Own parameters for IGES basic surfaces 128 and 190-198.
//!
//! `IGESGeom_ToolBSplineSurface::ReadOwnParams` (`cxx:46-268`) and the
//! `IGESSolid_Tool*Surface::ReadOwnParams` readers.

use occt_core::precision::PCONFUSION;

use super::iges_data_iges_model::EntityData;
use super::iges_data_param_reader::ParamReader;

pub(in crate::iges_read) fn parse_bspline_surface(
    pr: &mut ParamReader<'_>,
) -> Result<EntityData, String> {
    // `IGESGeom_ToolBSplineSurface.cxx:73-268`.
    let index_u = pr.read_integer()?;
    let index_v = pr.read_integer()?;
    let degree_u = pr.read_integer()?;
    let degree_v = pr.read_integer()?;
    let closed_u = parse_flag(pr)?;
    let closed_v = parse_flag(pr)?;
    let polynomial = parse_flag(pr)?;
    let periodic_u = parse_flag(pr)?;
    let periodic_v = parse_flag(pr)?;
    let n_knots_u = index_u + degree_u + 2;
    let n_knots_v = index_v + degree_v + 2;
    if n_knots_u < 0 || n_knots_v < 0 || index_u < 0 || index_v < 0 {
        return Err("IGES 128 index is negative".into());
    }
    let mut knots_u = Vec::with_capacity(n_knots_u as usize);
    for _ in 0..n_knots_u {
        knots_u.push(pr.read_real()?);
    }
    let mut knots_v = Vec::with_capacity(n_knots_v as usize);
    for _ in 0..n_knots_v {
        knots_v.push(pr.read_real()?);
    }
    let nu = (index_u + 1) as usize;
    let nv = (index_v + 1) as usize;
    let mut weights = vec![vec![0.0; nv]; nu];
    let mut bad_weight = false;
    for j in 0..nv {
        for i in 0..nu {
            let w = pr.read_real()?;
            if w < PCONFUSION {
                bad_weight = true;
            }
            weights[i][j] = w;
        }
    }
    if bad_weight {
        for row in &mut weights {
            for w in row {
                *w = 1.0;
            }
        }
    }
    let mut poles = vec![vec![[0.0; 3]; nv]; nu];
    for j in 0..nv {
        for i in 0..nu {
            poles[i][j] = pr.read_xyz()?;
        }
    }
    let u_min = pr.read_real()?;
    let v_min = pr.read_real()?;
    let u_max = pr.read_real()?;
    let v_max = pr.read_real()?;
    Ok(EntityData::BSplineSurface {
        degree_u,
        degree_v,
        closed_u,
        closed_v,
        polynomial,
        periodic_u,
        periodic_v,
        knots_u,
        knots_v,
        weights,
        poles,
        u_min,
        u_max,
        v_min,
        v_max,
    })
}

pub(in crate::iges_read) fn parse_plane_surface(
    pr: &mut ParamReader<'_>,
    form: i32,
) -> Result<EntityData, String> {
    let location = require_entity(pr, "190 location")?;
    let normal = require_entity(pr, "190 normal")?;
    let reference = if form == 1 {
        Some(require_entity(pr, "190 reference")?)
    } else {
        None
    };
    Ok(EntityData::PlaneSurface {
        location,
        normal,
        reference,
    })
}

pub(in crate::iges_read) fn parse_cylindrical_surface(
    pr: &mut ParamReader<'_>,
    form: i32,
) -> Result<EntityData, String> {
    let location = require_entity(pr, "192 location")?;
    let axis = require_entity(pr, "192 axis")?;
    let radius = pr.read_real()?;
    let reference = if form == 1 {
        Some(require_entity(pr, "192 reference")?)
    } else {
        None
    };
    Ok(EntityData::CylindricalSurface {
        location,
        axis,
        radius,
        reference,
    })
}

pub(in crate::iges_read) fn parse_conical_surface(
    pr: &mut ParamReader<'_>,
    form: i32,
) -> Result<EntityData, String> {
    let location = require_entity(pr, "194 location")?;
    let axis = require_entity(pr, "194 axis")?;
    let radius = pr.read_real()?;
    let semi_angle = pr.read_real()?;
    let reference = if form == 1 {
        Some(require_entity(pr, "194 reference")?)
    } else {
        None
    };
    Ok(EntityData::ConicalSurface {
        location,
        axis,
        radius,
        semi_angle,
        reference,
    })
}

pub(in crate::iges_read) fn parse_spherical_surface(
    pr: &mut ParamReader<'_>,
    form: i32,
) -> Result<EntityData, String> {
    let center = require_entity(pr, "196 center")?;
    let radius = pr.read_real()?;
    let (axis, reference) = if form == 1 {
        (
            Some(require_entity(pr, "196 axis")?),
            Some(require_entity(pr, "196 reference")?),
        )
    } else {
        (None, None)
    };
    Ok(EntityData::SphericalSurface {
        center,
        radius,
        axis,
        reference,
    })
}

pub(in crate::iges_read) fn parse_toroidal_surface(
    pr: &mut ParamReader<'_>,
    form: i32,
) -> Result<EntityData, String> {
    let center = require_entity(pr, "198 center")?;
    let axis = require_entity(pr, "198 axis")?;
    let major = pr.read_real()?;
    let minor = pr.read_real()?;
    let reference = if form == 1 {
        Some(require_entity(pr, "198 reference")?)
    } else {
        None
    };
    Ok(EntityData::ToroidalSurface {
        center,
        axis,
        major,
        minor,
        reference,
    })
}

fn parse_flag(pr: &mut ParamReader<'_>) -> Result<bool, String> {
    Ok(pr.read_integer()? != 0)
}

fn require_entity(pr: &mut ParamReader<'_>, what: &str) -> Result<usize, String> {
    pr.read_entity(false)?
        .ok_or_else(|| format!("IGES {what} pointer is null"))
}
