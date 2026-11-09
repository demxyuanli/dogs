//! `BRepOffset_Offset::Init`'s degenerated-edge singularity handling
//! (`BRepOffset_Offset.cxx:526-807`).
//!
//! An offset face can close onto a singular point (a cone apex, the centre of
//! an offset of a surface of revolution). The block finds the degenerated
//! edges of the source face, decides along which iso direction the
//! singularity runs, and either trims the offset surface back to the apex
//! (cone arm, `cxx:616-634`) or extends it until the generatrix meets the
//! apex line (offset-surface arm, `cxx:635-805`). `HasSingularity` and the two
//! apex points are then used by the wire rebuild (`init_face.rs`).
//!
//! Port guard (never a branch): `occ::down_cast<Geom_OffsetSurface>(TheSurf)`
//! (`cxx:642`, `:683`, `:724`, `:767`) answers a null handle once an earlier
//! block of the same arm trimmed `TheSurf` (a `Geom_RectangularTrimmedSurface`
//! is not a `Geom_OffsetSurface`), and the following `->BasisSurface()` would
//! dereference it. The port skips that block; `offset_basis_surface` returns
//! `None` for the same objects, so the detection is the same.
//!
//! `GeomLib::ExtendSurfByLength` (`cxx:661`, `:702`, `:743`, `:786`) makes the
//! surface *longer* along one iso direction; the port's
//! `geom_lib_extend_surf::extend_surf_by_length` answers the adjusted surface,
//! and a `None` there is OCCT's `Standard_ConstructionError`.

use std::sync::Arc;

use occt_core::elib::slib;
use occt_core::gp::{GpDir, GpLin, GpPnt, GpVec};
use occt_core::kernel::geomabs::JoinType;
use occt_core::precision::{CONFUSION, PCONFUSION, SQUARE_CONFUSION};
use occt_geom::{GeomLine, GeomRectangularTrimmedSurface, Surface};

use crate::abs::ShapeType;
use crate::boptools_2d;
use crate::brep_offset::BRepOffsetTool;
use crate::brep_tool::BRepTool;
use crate::brep_tools;
use crate::shape::{Edge, Face};
use crate::topo_tools::explorer_shapes;

/// `TolApex = 1.e-5` (`cxx:567`).
const TOL_APEX: f64 = 1.0e-5;

/// What `cxx:526-807` leaves behind for the wire rebuild (`cxx:829-868`,
/// `:914-956`).
#[derive(Debug, Clone, Copy)]
pub(crate) struct SingularityInfo {
    /// `HasSingularity` (`cxx:533`).
    pub has: bool,
    /// `MinApex` (`cxx:532`).
    pub min_apex: GpPnt,
    /// `MaxApex` (`cxx:532`).
    pub max_apex: GpPnt,
}

impl SingularityInfo {
    fn none() -> Self {
        Self {
            has: false,
            min_apex: GpPnt::zero(),
            max_apex: GpPnt::zero(),
        }
    }
}

/// `TheSurf->DynamicType() == STANDARD_TYPE(Geom_ConicalSurface)`.
///
/// `gp_cone()` answers through a `Geom_RectangularTrimmedSurface` (it mirrors
/// `GeomAdaptor_Surface::GetType()`), so the exact-type test also requires
/// that the surface is not trimmed.
fn is_exact_cone(s: &Arc<dyn Surface>) -> bool {
    s.rectangular_trimmed_spec().is_none() && s.gp_cone().is_some()
}

/// `BRepOffset_Offset.cxx:526-807`.
///
/// `src_surf` is the source surface `S` of `cxx:480-522` (used at `cxx:553-554`
/// to test the degenerated edges), `the_surf` starts as
/// `BRepOffset::Surface`'s result (`cxx:524`) and is trimmed or extended in
/// place.
pub(crate) fn process_singularities(
    face: &Face,
    src_surf: &Arc<dyn Surface>,
    the_surf: &mut Arc<dyn Surface>,
    offset_outside: bool,
    join_type: JoinType,
) -> Result<SingularityInfo, String> {
    let mut info = SingularityInfo::none();

    // `cxx:535`.
    let (uf1, uf2, vf1, vf2) = brep_tools::uv_bounds(face);

    // `cxx:536-538`.
    if (offset_outside && join_type == JoinType::Arc)
        || !(is_exact_cone(the_surf) || the_surf.is_offset_surface())
    {
        return Ok(info);
    }

    // `cxx:540-564`.
    let mut deg_edges: Vec<Edge> = Vec::new();
    for shape in explorer_shapes(&face.0, ShapeType::Edge) {
        let edge = Edge(shape);
        if !BRepTool::is_degenerated(&edge) {
            continue;
        }
        let Some((pc, a_f, a_l)) = boptools_2d::curve_on_surface_range(&edge, face) else {
            continue;
        };
        let a_fp = pc.value(a_f);
        let a_lp = pc.value(a_l);
        let a_f_pnt = src_surf.d0(a_fp.x(), a_fp.y());
        let a_l_pnt = src_surf.d0(a_lp.x(), a_lp.y());
        if a_f_pnt.square_distance(&a_l_pnt) < SQUARE_CONFUSION {
            deg_edges.push(edge);
        }
    }
    if deg_edges.is_empty() {
        return Ok(info);
    }

    // `cxx:567-615`: which iso direction carries the singularity.
    let Some((first_pc, fpar, lpar)) = boptools_2d::curve_on_surface_range(&deg_edges[0], face)
    else {
        return Ok(info);
    };
    let fp2d = first_pc.value(fpar);
    let lp2d = first_pc.value(lpar);
    let u_iso_degen = (fp2d.x() - lp2d.x()).abs() <= PCONFUSION;

    let (mut u_min, mut u_max, mut v_min, mut v_max) = (false, false, false, false);
    if deg_edges.len() == 2 {
        // `cxx:578-590`.
        if u_iso_degen {
            u_min = true;
            u_max = true;
        } else {
            v_min = true;
            v_max = true;
        }
    } else if u_iso_degen {
        // `cxx:593-603`.
        if (fp2d.x() - uf1).abs() <= CONFUSION {
            u_min = true;
        } else {
            u_max = true;
        }
    } else if (fp2d.y() - vf1).abs() <= CONFUSION {
        // `cxx:604-614`.
        v_min = true;
    } else {
        v_max = true;
    }

    if is_exact_cone(the_surf) {
        // `cxx:616-634`.
        let Some(cone) = the_surf.gp_cone() else {
            return Ok(info);
        };
        let apex = cone.apex();
        let (_, v_apex) = slib::cone_parameters(
            &cone.position(),
            cone.radius(),
            cone.semi_angle(),
            &apex,
        );
        if v_min {
            let trimmed =
                GeomRectangularTrimmedSurface::uv(the_surf.clone(), uf1, uf2, v_apex, vf2)?;
            *the_surf = Arc::new(trimmed);
            info.min_apex = apex;
            info.has = true;
        } else if v_max {
            let trimmed =
                GeomRectangularTrimmedSurface::uv(the_surf.clone(), uf1, uf2, vf1, v_apex)?;
            *the_surf = Arc::new(trimmed);
            info.max_apex = apex;
            info.has = true;
        }
        return Ok(info);
    }

    // `cxx:635-805`: the `Geom_OffsetSurface` arm, four independent blocks in
    // OCCT's order, each working on the surface the previous one produced.

    // `cxx:637-677`: UminDegen.
    if u_min {
        if let Some(iso) = the_surf.u_iso_curve(uf1) {
            if BRepOffsetTool::gabarit(iso.as_ref()) > TOL_APEX {
                let p_first = the_surf.d0(uf1, vf1);
                let papex = basis_first(the_surf, uf1, vf1)?;
                let p_quart = the_surf.d0(uf1, 0.75 * vf1 + 0.25 * vf2);
                let p_mid = the_surf.d0(uf1, 0.5 * (vf1 + vf2));
                let dn = basis_dn(the_surf, uf1, vf1, 1, 0)?;
                let (length, pint1) = apex_length(papex, p_first, p_quart, p_mid, dn)?;
                if offset_outside {
                    let a_surf: Arc<dyn Surface> = Arc::new(GeomRectangularTrimmedSurface::uv(
                        the_surf.clone(),
                        uf1,
                        uf2,
                        vf1,
                        vf2,
                    )?);
                    let extended =
                        occt_geom::geom_lib_extend_surf::extend_surf_by_length(
                            &a_surf, length, 1, true, false,
                        )
                        .ok_or("GeomLib::ExtendSurfByLength failed")?;
                    let u1 = extended.u_range().0;
                    info.min_apex = extended.d0(u1, vf1);
                    *the_surf = extended;
                } else {
                    let Some(viso) = the_surf.v_iso_curve(vf1) else {
                        return Err("BRepOffset_Offset: no VIso at UminDegen".into());
                    };
                    let Some(proj) =
                        occt_geom::geom_api::project_point_on_curve(viso.as_ref(), &pint1, CONFUSION)
                    else {
                        return Err("BRepOffset_Offset: projection failed at UminDegen".into());
                    };
                    let new_first_u = proj.parameter;
                    let trimmed = GeomRectangularTrimmedSurface::uv(
                        the_surf.clone(),
                        new_first_u,
                        uf2,
                        vf1,
                        vf2,
                    )?;
                    info.min_apex = trimmed.d0(new_first_u, vf1);
                    *the_surf = Arc::new(trimmed);
                }
                info.has = true;
            }
        }
    }

    // `cxx:678-718`: UmaxDegen.
    if u_max {
        if let Some(iso) = the_surf.u_iso_curve(uf2) {
            if BRepOffsetTool::gabarit(iso.as_ref()) > TOL_APEX {
                let p_first = the_surf.d0(uf2, vf1);
                let papex = basis_first(the_surf, uf2, vf1)?;
                let p_quart = the_surf.d0(uf2, 0.75 * vf1 + 0.25 * vf2);
                let p_mid = the_surf.d0(uf2, 0.5 * (vf1 + vf2));
                let dn = basis_dn(the_surf, uf2, vf1, 1, 0)?;
                let (length, pint1) = apex_length(papex, p_first, p_quart, p_mid, dn)?;
                if offset_outside {
                    let a_surf: Arc<dyn Surface> = Arc::new(GeomRectangularTrimmedSurface::uv(
                        the_surf.clone(),
                        uf1,
                        uf2,
                        vf1,
                        vf2,
                    )?);
                    let extended =
                        occt_geom::geom_lib_extend_surf::extend_surf_by_length(
                            &a_surf, length, 1, true, true,
                        )
                        .ok_or("GeomLib::ExtendSurfByLength failed")?;
                    let u2 = extended.u_range().1;
                    info.max_apex = extended.d0(u2, vf1);
                    *the_surf = extended;
                } else {
                    let Some(viso) = the_surf.v_iso_curve(vf1) else {
                        return Err("BRepOffset_Offset: no VIso at UmaxDegen".into());
                    };
                    let Some(proj) =
                        occt_geom::geom_api::project_point_on_curve(viso.as_ref(), &pint1, CONFUSION)
                    else {
                        return Err("BRepOffset_Offset: projection failed at UmaxDegen".into());
                    };
                    let new_last_u = proj.parameter;
                    let trimmed = GeomRectangularTrimmedSurface::uv(
                        the_surf.clone(),
                        uf1,
                        new_last_u,
                        vf1,
                        vf2,
                    )?;
                    info.max_apex = trimmed.d0(new_last_u, vf1);
                    *the_surf = Arc::new(trimmed);
                }
                info.has = true;
            }
        }
    }

    // `cxx:719-761`: VminDegen.
    if v_min {
        if let Some(iso) = the_surf.v_iso_curve(vf1) {
            if BRepOffsetTool::gabarit(iso.as_ref()) > TOL_APEX {
                let p_first = the_surf.d0(uf1, vf1);
                let papex = basis_first(the_surf, uf1, vf1)?;
                let p_quart = the_surf.d0(0.75 * uf1 + 0.25 * uf2, vf1);
                let p_mid = the_surf.d0(0.5 * (uf1 + uf2), vf1);
                let dn = basis_dn(the_surf, uf1, vf1, 0, 1)?;
                let (length, pint1) = apex_length(papex, p_first, p_quart, p_mid, dn)?;
                if offset_outside {
                    let a_surf: Arc<dyn Surface> = Arc::new(GeomRectangularTrimmedSurface::uv(
                        the_surf.clone(),
                        uf1,
                        uf2,
                        vf1,
                        vf2,
                    )?);
                    let extended =
                        occt_geom::geom_lib_extend_surf::extend_surf_by_length(
                            &a_surf, length, 1, false, false,
                        )
                        .ok_or("GeomLib::ExtendSurfByLength failed")?;
                    let v1 = extended.v_range().0;
                    info.min_apex = extended.d0(uf1, v1);
                    *the_surf = extended;
                } else {
                    let Some(uiso) = the_surf.u_iso_curve(uf1) else {
                        return Err("BRepOffset_Offset: no UIso at VminDegen".into());
                    };
                    let Some(proj) =
                        occt_geom::geom_api::project_point_on_curve(uiso.as_ref(), &pint1, CONFUSION)
                    else {
                        return Err("BRepOffset_Offset: projection failed at VminDegen".into());
                    };
                    let new_first_v = proj.parameter;
                    let trimmed = GeomRectangularTrimmedSurface::uv(
                        the_surf.clone(),
                        uf1,
                        uf2,
                        new_first_v,
                        vf2,
                    )?;
                    info.min_apex = trimmed.d0(uf1, new_first_v);
                    *the_surf = Arc::new(trimmed);
                }
                info.has = true;
            }
        }
    }

    // `cxx:762-804`: VmaxDegen.
    if v_max {
        if let Some(iso) = the_surf.v_iso_curve(vf2) {
            if BRepOffsetTool::gabarit(iso.as_ref()) > TOL_APEX {
                let p_first = the_surf.d0(uf1, vf2);
                let papex = basis_first(the_surf, uf1, vf2)?;
                let p_quart = the_surf.d0(0.75 * uf1 + 0.25 * uf2, vf2);
                let p_mid = the_surf.d0(0.5 * (uf1 + uf2), vf2);
                let dn = basis_dn(the_surf, uf1, vf2, 0, 1)?;
                let (length, pint1) = apex_length(papex, p_first, p_quart, p_mid, dn)?;
                if offset_outside {
                    let a_surf: Arc<dyn Surface> = Arc::new(GeomRectangularTrimmedSurface::uv(
                        the_surf.clone(),
                        uf1,
                        uf2,
                        vf1,
                        vf2,
                    )?);
                    let extended =
                        occt_geom::geom_lib_extend_surf::extend_surf_by_length(
                            &a_surf, length, 1, false, true,
                        )
                        .ok_or("GeomLib::ExtendSurfByLength failed")?;
                    let v2 = extended.v_range().1;
                    info.max_apex = extended.d0(uf1, v2);
                    *the_surf = extended;
                } else {
                    let Some(uiso) = the_surf.u_iso_curve(uf1) else {
                        return Err("BRepOffset_Offset: no UIso at VmaxDegen".into());
                    };
                    let Some(proj) =
                        occt_geom::geom_api::project_point_on_curve(uiso.as_ref(), &pint1, CONFUSION)
                    else {
                        return Err("BRepOffset_Offset: projection failed at VmaxDegen".into());
                    };
                    let new_last_v = proj.parameter;
                    let trimmed = GeomRectangularTrimmedSurface::uv(
                        the_surf.clone(),
                        uf1,
                        uf2,
                        vf1,
                        new_last_v,
                    )?;
                    info.max_apex = trimmed.d0(uf1, new_last_v);
                    *the_surf = Arc::new(trimmed);
                }
                info.has = true;
            }
        }
    }

    Ok(info)
}

/// `BasisSurf->Value(u, v)` for `BasisSurf = down_cast<Geom_OffsetSurface>(
/// TheSurf)->BasisSurface()` (`cxx:642-645`, `:683-686`, `:724-727`,
/// `:767-770`). `Err` is the null down-cast guard described in the module
/// header.
fn basis_first(the_surf: &Arc<dyn Surface>, u: f64, v: f64) -> Result<GpPnt, String> {
    let basis = the_surf
        .offset_basis_surface()
        .ok_or("BRepOffset_Offset: null Geom_OffsetSurface down-cast")?;
    Ok(basis.d0(u, v))
}

/// `BasisSurf->DN(u, v, nu, nv)` (`cxx:651`, `:692`, `:733`, `:776`).
fn basis_dn(
    the_surf: &Arc<dyn Surface>,
    u: f64,
    v: f64,
    nu: i32,
    nv: i32,
) -> Result<GpVec, String> {
    let basis = the_surf
        .offset_basis_surface()
        .ok_or("BRepOffset_Offset: null Geom_OffsetSurface down-cast")?;
    match basis.try_eval_dn(u, v, nu, nv) {
        Ok(vec) if vec.magnitude() > 0.0 => Ok(vec),
        _ => Err("BRepOffset_Offset: undefined surface derivative".into()),
    }
}

/// `cxx:649-656`: the apex line meets the generatrix line; the answer is
/// `length = Pfirst.Distance(Pint1)` together with `Pint1`, which the
/// projection arm needs (`cxx:670-671`, `:711-712`, `:752-753`, `:795-796`).
///
/// `None` is the null `gp_Dir` of `Geom_Line` (`cxx:650`) or an empty
/// `GeomAPI_ExtremaCurveCurve` (`cxx:655`), both of which OCCT lets throw.
fn apex_length(
    papex: GpPnt,
    pfirst: GpPnt,
    pquart: GpPnt,
    pmid: GpPnt,
    dn: GpVec,
) -> Result<(f64, GpPnt), String> {
    let dir_apex = GpVec::from_pnts(&pfirst, &pquart).crossed(&GpVec::from_pnts(&pfirst, &pmid));
    let apex_dir = GpDir::from_vec(&dir_apex)
        .map_err(|_| "BRepOffset_Offset: null apex line direction".to_string())?;
    let gen_dir = GpDir::from_vec(&dn)
        .map_err(|_| "BRepOffset_Offset: null generatrix direction".to_string())?;
    let line_apex = GeomLine::new(GpLin::from_pnt_dir(papex, apex_dir));
    let line_gen = GeomLine::new(GpLin::from_pnt_dir(pfirst, gen_dir));
    let extrema = occt_geom::extrema::curve_curve_extrema_with_samples(&line_gen, &line_apex, 0);
    let pair = extrema
        .first()
        .ok_or("BRepOffset_Offset: empty GeomAPI_ExtremaCurveCurve")?;
    Ok((pfirst.distance(&pair.p1), pair.p1))
}
