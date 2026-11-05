//! `totalsurf` — the file static of `BRepFill_NSections.cxx:136-358`.
//!
//! Concatenates every section into one `Geom_BSplineCurve` (per column of the
//! section table), profiles the lot through `GeomFill_SectionGenerator`, and
//! approximates the result into the single `Geom_BSplineSurface` all the
//! `GeomFill_NSections` laws of one `BRepFill_NSections` share.
//!
//! Translation notes:
//! - `NCollection_Array2<TopoDS_Shape>` with `Value(i, j)` is
//!   `Vec<Vec<Option<TopoShape>>>` indexed `[i - 1][j - 1]`; a `None` slot is
//!   the null shape OCCT down-casts with `TopoDS::Edge` (which raises
//!   `Standard_NullObject`), answered as an `Err` here.
//! - `GeomFill_AppSurf` is [`AppSurf`] built with exactly the arguments of
//!   `cxx:261`: `(2, 6, myPres3d, myPres3d, 0, true)`. Unlike
//!   `GeomFill_NSections::ComputeSurface` (`GeomFill_NSections.cxx:573-579`)
//!   there is **no** `SetContinuity` call, so the `AppBlend_AppSurf` default
//!   `GeomAbs_C2` stays in place.
//! - `Geom_BSplineSurface(Poles, Weights, UKnots, VKnots, UMults, VMults,
//!   UDeg, VDeg, UPeriodic, VPeriodic)` (`Geom_BSplineSurface.cxx:230-292`) is
//!   [`GeomBSplineSurface::from_parts`] over the flat knot vectors
//!   `BSplCLib::KnotSequence` produces (`Geom_BSplineSurface.cxx:1160-1195`).
//!   The 12-argument constructor keeps the weights it is handed (it only runs
//!   `Rational()` to set the two flags), so a weight grid whose entries all
//!   agree is dropped to `None` — the non-rational surface the `cxx:299`
//!   comment asks for; equal weights cancel out of the rational evaluation,
//!   so the geometry is unchanged.

use occt_core::gp::GpPnt;
use occt_core::precision::PCONFUSION;
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::bspline_surface::GeomBSplineSurface;
use occt_geom::convert_bspl::CompCurveToBSplineCurve;
use occt_core::convert::ParameterisationType;
use occt_geom::geom_fill::app_surf::AppSurf;
use occt_geom::geom_fill::line::Line;
use occt_geom::geom_fill::section_generator::SectionGenerator;

use crate::brep_fill_n_sections::edge_to_bspline::{edge_to_bspline, point_curve_of_edge};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::common_vertex;

/// `totalsurf(shapes, NbSects, NbEdges, params, w1Point, w2Point, uClosed,
/// vClosed, myPres3d)` (`cxx:136-358`).
#[allow(clippy::too_many_arguments)]
pub fn total_surf(
    shapes: &[Vec<Option<TopoShape>>],
    nb_sects: usize,
    nb_edges: usize,
    params: &[f64],
    w1_point: bool,
    w2_point: bool,
    u_closed: bool,
    v_closed: bool,
    pres3d: f64,
) -> Result<GeomBSplineSurface, String> {
    if nb_sects == 0 || nb_edges == 0 {
        return Err("Standard_OutOfRange: BRepFill_NSections::totalsurf".to_string());
    }

    let mut section = SectionGenerator::new();
    let mut bs1: Option<GeomBSplineCurve> = None;
    let mut jdeb = 1usize;
    let mut jfin = nb_sects;

    // `cxx:152-169`: a punctual first section is profiled as a point curve.
    if w1_point {
        jdeb += 1;
        let edge = shape_edge(shapes, 1, 1)?;
        section.add_curve(&point_curve_of_edge(&edge)?)?;
    }

    // `cxx:171-174`.
    if w2_point {
        jfin -= 1;
    }

    // `cxx:176-232`.
    let mut j = jdeb;
    while j <= jfin {
        // `cxx:179-182`: the closing column of a looping set of sections
        // reuses the first one instead of being read again.
        if j == jfin && v_closed {
            let closed = bs1
                .clone()
                .ok_or_else(|| "Standard_NullObject: BRepFill_NSections::totalsurf".to_string())?;
            section.add_curve(&closed)?;
            j += 1;
            continue;
        }

        // `cxx:186-190`: the first edge initialises the concatenator.
        let mut prev_edge = shape_edge(shapes, 1, j)?;
        let curv_bs = edge_to_bspline(&prev_edge)?;
        let mut comp = CompCurveToBSplineCurve::from_bspline(
            curv_bs,
            ParameterisationType::QuasiAngular,
        );

        for i in 2..=nb_edges {
            let next_edge = shape_edge(shapes, i, j)?;
            let curv_bs = edge_to_bspline(&next_edge)?;

            // `cxx:200-209`: the tolerance is the common vertex's, or
            // `Precision::Confusion()` when the edges do not share one.
            let mut eps_v = occt_core::precision::CONFUSION;
            if let Some(v) = common_vertex(&prev_edge, &next_edge) {
                eps_v = BRepTool::vertex_tolerance(&Vertex(v));
            }
            // `cxx:211-216`: concatenation, retried with `200 * epsV`.
            let mut bof = comp
                .add(&curv_bs, eps_v, true, false, 1)
                .unwrap_or(false);
            if !bof {
                bof = comp
                    .add(&curv_bs, 200.0 * eps_v, true, false, 1)
                    .unwrap_or(false);
            }
            prev_edge = next_edge;
        }

        // `cxx:221-229`: the concatenated section, remembered when the first
        // column of a looping set is done.
        let bs = comp
            .bspline_curve()
            .cloned()
            .ok_or_else(|| "Standard_ConstructionError: GeomConvert_CompCurveToBSplineCurve".to_string())?;
        section.add_curve(&bs)?;
        if j == jdeb && v_closed {
            bs1 = Some(bs);
        }
        j += 1;
    }

    // `cxx:234-248`: a punctual last section is profiled as a point curve too.
    if w2_point {
        let edge = shape_edge(shapes, nb_edges, nb_sects)?;
        section.add_curve(&point_curve_of_edge(&edge)?)?;
    }

    // `cxx:250-257`.
    section.set_param(params);
    section.perform(PCONFUSION)?;

    let line = Line::with_nb_points(nb_sects as i32);
    // `cxx:259-264`: `GeomFill_AppSurf(2, 6, myPres3d, myPres3d, 0, true)`.
    let mut approx = AppSurf::with_settings(2, 6, pres3d, pres3d, 0, true);
    approx.perform(&line, &section, true)?;
    if !approx.is_done() {
        // `SurfPoles()` of a not-done `AppBlend_AppSurf` raises
        // `StdFail_NotDone` (`AppBlend_AppSurf.gxx:1073-1090`).
        return Err("StdFail_NotDone: GeomFill_AppSurf::Perform".to_string());
    }

    // `cxx:265-287`: the periodic drops of the seam mults and pole column.
    let u_periodic = u_closed;
    let v_periodic = v_closed;
    let poles_all = approx.surf_poles().to_vec();
    let weights_all = approx.surf_weights().to_vec();
    let u_knots = approx.surf_u_knots().to_vec();
    let v_knots = approx.surf_v_knots().to_vec();
    let mut u_mults = approx.surf_u_mults().to_vec();
    let mut v_mults = approx.surf_v_mults().to_vec();
    let deg_u = approx.u_degree() as usize;
    let deg_v = approx.v_degree() as usize;

    let mut nup = poles_all.len();
    let mut nvp = poles_all.first().map_or(0, |row| row.len());

    if u_periodic {
        let nbuk = u_knots.len();
        if nbuk > 0 {
            u_mults[0] -= 1;
            u_mults[nbuk - 1] -= 1;
        }
        nup = nup.saturating_sub(1);
    }
    // `cxx:280-286`.
    if v_periodic {
        let nbvk = v_knots.len();
        if nbvk > 0 {
            v_mults[0] -= 1;
            v_mults[nbvk - 1] -= 1;
        }
        nvp = nvp.saturating_sub(1);
    }

    // `cxx:288-297`.
    if nup == 0 || nvp == 0 {
        return Err("StdFail_NotDone: GeomFill_AppSurf::Perform".to_string());
    }
    let mut poles: Vec<Vec<GpPnt>> = vec![vec![GpPnt::zero(); nvp]; nup];
    let mut weights: Vec<Vec<f64>> = vec![vec![1.0; nvp]; nup];
    for jj in 0..nvp {
        for ii in 0..nup {
            poles[ii][jj] = poles_all[ii][jj];
            weights[ii][jj] = weights_all[ii][jj];
        }
    }

    // `cxx:299-334`: is the surface really rational? `Vrational` walks the U
    // index inside one V column, `Urational` walks the V index inside one U row.
    const TOL_EPS: f64 = 1.0e-13;
    let mut v_rational = false;
    for jj in 0..nvp {
        if v_rational {
            break;
        }
        for ii in 0..nup.saturating_sub(1) {
            if (weights[ii][jj] - weights[ii + 1][jj]).abs() > TOL_EPS {
                v_rational = true;
                break;
            }
        }
    }
    let mut u_rational = false;
    for ii in 0..nup {
        if u_rational {
            break;
        }
        for jj in 0..nvp.saturating_sub(1) {
            if (weights[ii][jj] - weights[ii][jj + 1]).abs() > TOL_EPS {
                u_rational = true;
                break;
            }
        }
    }

    // `cxx:335-345`: a surface whose weights all agree keeps a single weight.
    if !v_rational && !u_rational {
        let the_weight = weights[0][0];
        for row in weights.iter_mut() {
            for w in row.iter_mut() {
                *w = the_weight;
            }
        }
    }

    // `cxx:347-357`.
    let flat_u = if u_periodic {
        occt_core::bspl::knots::knot_sequence_periodic(&u_knots, &u_mults, deg_u as i32)
    } else {
        occt_core::bspl::banded_interp::knot_sequence(&u_knots, &u_mults, deg_u as i32)
    };
    let flat_v = if v_periodic {
        occt_core::bspl::knots::knot_sequence_periodic(&v_knots, &v_mults, deg_v as i32)
    } else {
        occt_core::bspl::banded_interp::knot_sequence(&v_knots, &v_mults, deg_v as i32)
    };
    let weight_grid = if u_rational || v_rational {
        Some(weights)
    } else {
        None
    };
    Ok(GeomBSplineSurface::from_parts(
        poles,
        weight_grid,
        flat_u,
        flat_v,
        deg_u,
        deg_v,
        u_periodic,
        v_periodic,
    ))
}

/// `TopoDS::Edge(shapes.Value(i, j))` with the 1-based `(i, j)` of the OCCT
/// table; a null slot raises `Standard_NullObject` in OCCT.
fn shape_edge(shapes: &[Vec<Option<TopoShape>>], i: usize, j: usize) -> Result<Edge, String> {
    shapes
        .get(i - 1)
        .and_then(|row| row.get(j - 1))
        .and_then(|slot| slot.clone())
        .and_then(Edge::wrap)
        .ok_or_else(|| "Standard_NullObject: BRepFill_NSections::totalsurf".to_string())
}
