//! `BRepAlgo::ConcatenateWireC0(W)` (`BRepAlgo.cxx:404-802`).
//!
//! The wire is first healed by `ShapeFix_Shape` at `Precision::Confusion()`
//! (`cxx:411-417`), then walked with `BRepTools_WireExplorer` into per-curve
//! parameter sequences (`cxx:419-425`). Two consecutive edges whose *curves*
//! coincide - the same handle (`cxx:469-473`) or, for the same curve type, an
//! equal line / circle / ellipse / hyperbola / parabola within `LinTol` and
//! `AngTol` (`cxx:476-571`) - are merged into one parameter range instead of a
//! new entry, with the periodic adjustment of `cxx:578-647`. The collected
//! ranges are then trimmed onto their bases (`cxx:697-710`), C0 knots are
//! fixed, and the B-splines are joined by `GeomConvert::ConcatC1`
//! (`cxx:748-751`) and, when more than one remains,
//! `GeomConvert_CompCurveToBSplineCurve` (`cxx:766-771`). The result is a
//! single edge (`cxx:779-786`), or a copy of the only curve's edge
//! (`cxx:787-793`), reversed when the wire started REVERSED (`cxx:795-799`).
//!
//! `Ok(None)` is the early null `ResEdge` of `cxx:682-685` (no edge carried a
//! 3D curve). A `BRepLib_MakeEdge` that did not complete cannot be stored as
//! an `Edge`, so it is `Err`.

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::elib::clib;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{ANGULAR, CONFUSION, PCONFUSION};
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::convert_bspl::{curve_to_bspline_curve, CompCurveToBSplineCurve};
use occt_geom::curve::Curve;
use occt_geom::geom_convert_concat::{c0_bspline_to_c1_bspline_curve, concat_c1};
use occt_geom::geom_lprop::continuity;
use occt_geom::trimmed::GeomTrimmedCurveBasis;

use crate::abs::Orientation;
use crate::brep_lib_build_curves3d;
use crate::brep_lib_make_edge;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use crate::shape::{Edge, Vertex, Wire};
use crate::shhealing::ShapeFixShape;
use crate::topo_tools::is_same;

pub fn concatenate_wire_c0(a_wire: &Wire) -> Result<Option<Edge>, String> {
    // `cxx:406-407`.
    let lin_tol = CONFUSION;
    let ang_tol = ANGULAR;

    // `cxx:410-411`.
    let mut the_wire = a_wire.clone();
    brep_lib_build_curves3d::build_curves3d(&the_wire.0);

    // `cxx:412-417`: `SetPrecision(LinTol)`, `SetMaxTolerance(LinTol)`,
    // `Perform()`. `SetMinTolerance` is not called, so the wire tool keeps the
    // `ShapeFix_Root` default `Precision::Confusion()`.
    let mut fixer = ShapeFixShape::new(&the_wire.0);
    fixer.set_precision(lin_tol);
    fixer.set_max_tolerance(lin_tol);
    fixer.perform();
    the_wire = Wire(fixer.shape().clone());

    // `cxx:419-425`.
    let mut curve_seq: Vec<Arc<dyn Curve>> = Vec::new();
    let mut fpar_seq: Vec<f64> = Vec::new();
    let mut lpar_seq: Vec<f64> = Vec::new();
    let mut tol_seq: Vec<f64> = Vec::new();
    let mut is_fwd_seq: Vec<bool> = Vec::new();
    let mut cur_type = GeomAbsCurveType::OtherCurve;
    let mut first_vertex: Option<Vertex> = None;

    let mut wexp = WireExplorer::of_wire(&the_wire);
    while wexp.more() {
        let edge = wexp
            .current_edge()
            .ok_or_else(|| "BRepAlgo ConcatenateWireC0: empty edge".to_string())?;
        // `cxx:432-438`: `BRep_Tool::Curve(anEdge, fpar, lpar)`, skipping an
        // edge whose 3D curve is null.
        let Some(a_curve) = BRepTool::edge_curve(&edge) else {
            wexp.next();
            continue;
        };
        let (mut fpar, mut lpar) = BRepTool::edge_parameters(&edge);

        // `cxx:440-442`: `GeomAdaptor_Curve aGACurve(aCurve)`; `GetType()` and
        // `Curve()` are the adaptor's - the trim is unwrapped.
        let (basis_curve, _, _) = brep_lib_make_edge::adaptor_curve(&a_curve);
        let a_type = curve_kind(basis_curve.as_ref());
        let is_fwd = wexp.orientation()? != Orientation::Reversed;

        // `cxx:444-451`.
        if basis_curve.is_periodic() {
            clib::adjust_periodic(
                basis_curve.first_parameter(),
                basis_curve.last_parameter(),
                PCONFUSION,
                &mut fpar,
                &mut lpar,
            );
        }

        if curve_seq.is_empty() {
            // `cxx:453-461`.
            curve_seq.push(a_curve);
            fpar_seq.push(fpar);
            lpar_seq.push(lpar);
            is_fwd_seq.push(is_fwd);
            cur_type = a_type;
            first_vertex = wexp.current_vertex().cloned();
        } else {
            // `cxx:465-467`: `RealFirst()` / `RealLast()`.
            let mut is_same_curve = false;
            let mut new_fpar = -1.0e100;
            let mut new_lpar = 1.0e100;
            let prev = curve_seq
                .last()
                .expect("curve sequence is not empty")
                .clone();
            let (prev_basis, _, _) = brep_lib_make_edge::adaptor_curve(&prev);

            if Arc::ptr_eq(&a_curve, &prev) {
                // `cxx:469-474`: the same `Geom_Curve` handle.
                new_fpar = fpar;
                new_lpar = lpar;
                is_same_curve = true;
            } else if a_type == cur_type {
                // `cxx:476-571`: the type-specific coincidence tests.
                match a_type {
                    GeomAbsCurveType::Line => {
                        if let (Some(a_line), Some(prev_line)) =
                            (basis_curve.gp_line(), prev_basis.gp_line())
                        {
                            if a_line.contains(&prev_line.location(), lin_tol)
                                && a_line
                                    .direction()
                                    .is_parallel_tol(&prev_line.direction(), ang_tol)
                            {
                                let p1 = clib::line_value(&a_line, fpar);
                                let p2 = clib::line_value(&a_line, lpar);
                                new_fpar = clib::line_parameter(&prev_line, &p1);
                                new_lpar = clib::line_parameter(&prev_line, &p2);
                                is_same_curve = true;
                            }
                        }
                    }
                    GeomAbsCurveType::Circle => {
                        if let (Some(a_circle), Some(prev_circle)) =
                            (basis_curve.gp_circ(), prev_basis.gp_circ())
                        {
                            if a_circle.location().distance(&prev_circle.location()) <= lin_tol
                                && (a_circle.radius() - prev_circle.radius()).abs() <= lin_tol
                                && a_circle
                                    .axis()
                                    .direction()
                                    .is_parallel_tol(&prev_circle.axis().direction(), ang_tol)
                            {
                                let p1 = clib::circle_value(&a_circle, fpar);
                                let p2 = clib::circle_value(&a_circle, lpar);
                                let prev_pos = prev_circle.position();
                                new_fpar = clib::circle_parameter(&prev_pos, &p1);
                                new_lpar = clib::circle_parameter(&prev_pos, &p2);
                                is_same_curve = true;
                            }
                        }
                    }
                    GeomAbsCurveType::Ellipse => {
                        if let (Some(an_ellipse), Some(prev_ellipse)) =
                            (basis_curve.gp_ellipse(), prev_basis.gp_ellipse())
                        {
                            if an_ellipse.focus1().distance(&prev_ellipse.focus1()) <= lin_tol
                                && an_ellipse.focus2().distance(&prev_ellipse.focus2()) <= lin_tol
                                && (an_ellipse.major_radius() - prev_ellipse.major_radius()).abs()
                                    <= lin_tol
                                && (an_ellipse.minor_radius() - prev_ellipse.minor_radius()).abs()
                                    <= lin_tol
                                && an_ellipse
                                    .axis()
                                    .direction()
                                    .is_parallel_tol(&prev_ellipse.axis().direction(), ang_tol)
                            {
                                let p1 = clib::ellipse_value(&an_ellipse, fpar);
                                let p2 = clib::ellipse_value(&an_ellipse, lpar);
                                new_fpar = clib::parameter_elips(&prev_ellipse, &p1);
                                new_lpar = clib::parameter_elips(&prev_ellipse, &p2);
                                is_same_curve = true;
                            }
                        }
                    }
                    GeomAbsCurveType::Hyperbola => {
                        if let (Some(a_hypr), Some(prev_hypr)) =
                            (basis_curve.gp_hyperbola(), prev_basis.gp_hyperbola())
                        {
                            if a_hypr.focus1().distance(&prev_hypr.focus1()) <= lin_tol
                                && a_hypr.focus2().distance(&prev_hypr.focus2()) <= lin_tol
                                && (a_hypr.major_radius - prev_hypr.major_radius).abs() <= lin_tol
                                && (a_hypr.minor_radius - prev_hypr.minor_radius).abs() <= lin_tol
                                && a_hypr
                                    .axis()
                                    .direction()
                                    .is_parallel_tol(&prev_hypr.axis().direction(), ang_tol)
                            {
                                let p1 = clib::hyperbola_value(&a_hypr, fpar);
                                let p2 = clib::hyperbola_value(&a_hypr, lpar);
                                new_fpar = clib::parameter_hypr(&prev_hypr, &p1);
                                new_lpar = clib::parameter_hypr(&prev_hypr, &p2);
                                is_same_curve = true;
                            }
                        }
                    }
                    GeomAbsCurveType::Parabola => {
                        if let (Some(a_parab), Some(prev_parab)) =
                            (basis_curve.gp_parabola(), prev_basis.gp_parabola())
                        {
                            if a_parab.location().distance(&prev_parab.location()) <= lin_tol
                                && a_parab.focus().distance(&prev_parab.focus()) <= lin_tol
                                && (a_parab.focal - prev_parab.focal).abs() <= lin_tol
                                && a_parab
                                    .axis()
                                    .direction()
                                    .is_parallel_tol(&prev_parab.axis().direction(), ang_tol)
                            {
                                let p1 = clib::parabola_value(&a_parab, fpar);
                                let p2 = clib::parabola_value(&a_parab, lpar);
                                new_fpar = clib::parameter_parab(&prev_parab, &p1);
                                new_lpar = clib::parameter_parab(&prev_parab, &p2);
                                is_same_curve = true;
                            }
                        }
                    }
                    _ => {}
                }
            }

            if is_same_curve {
                // `cxx:575-643`.
                let is_same_dir = is_fwd
                    == *is_fwd_seq
                        .last()
                        .expect("the sequence has an entry per curve");
                let prev_is_fwd = *is_fwd_seq
                    .last()
                    .expect("the sequence has an entry per curve");

                if basis_curve.is_periodic() {
                    // `cxx:578-626`.
                    let a_period = basis_curve.period();
                    if is_same_dir {
                        // `cxx:583-589`.
                        while new_fpar > new_lpar {
                            new_fpar -= a_period;
                        }
                    } else {
                        // `cxx:590-601`.
                        while new_lpar > new_fpar {
                            new_lpar -= a_period;
                        }
                        std::mem::swap(&mut new_lpar, &mut new_fpar);
                    }

                    // `cxx:604-626`.
                    if prev_is_fwd {
                        // The current curve comes after the previous one.
                        let l_par = *lpar_seq
                            .last()
                            .expect("the sequence has an entry per curve");
                        clib::adjust_periodic(
                            l_par,
                            l_par + a_period,
                            PCONFUSION,
                            &mut new_fpar,
                            &mut new_lpar,
                        );
                    } else {
                        // The current curve comes before the previous one.
                        let f_par = *fpar_seq
                            .last()
                            .expect("the sequence has an entry per curve");
                        clib::adjust_periodic(
                            f_par - a_period,
                            f_par,
                            PCONFUSION,
                            &mut new_fpar,
                            &mut new_lpar,
                        );
                    }
                } else if !is_same_dir {
                    // `cxx:627-632`: non-periodic curves with opposite
                    // directions.
                    std::mem::swap(&mut new_lpar, &mut new_fpar);
                }

                // `cxx:634-642`: the merge writes into the *previous* entry.
                if prev_is_fwd {
                    *lpar_seq
                        .last_mut()
                        .expect("the sequence has an entry per curve") = new_lpar;
                } else {
                    *fpar_seq
                        .last_mut()
                        .expect("the sequence has an entry per curve") = new_fpar;
                }
            } else {
                // `cxx:644-655`: add the new curve.
                curve_seq.push(a_curve);
                fpar_seq.push(fpar);
                lpar_seq.push(lpar);
                is_fwd_seq.push(is_fwd);
                // `cxx:653`: `BRep_Tool::Tolerance(wexp.CurrentVertex())`.
                let cv = wexp
                    .current_vertex()
                    .ok_or_else(|| "BRepAlgo ConcatenateWireC0: null vertex".to_string())?;
                tol_seq.push(BRepTool::vertex_tolerance(cv));
                cur_type = a_type;
            }
        }
        wexp.next();
    }

    // `cxx:658-659`.
    let last_vertex = wexp.current_vertex().cloned();
    tol_seq.push(
        last_vertex
            .as_ref()
            .map(BRepTool::vertex_tolerance)
            .ok_or_else(|| "BRepAlgo ConcatenateWireC0: null last vertex".to_string())?,
    );

    // `cxx:661-667`.
    let is_reverse = is_fwd_seq.first().map(|f| !f).unwrap_or(false);

    // `cxx:669-681`.
    let (first_final, last_final) = if is_reverse {
        (last_vertex.clone(), first_vertex.clone())
    } else {
        (first_vertex.clone(), last_vertex.clone())
    };
    let first_final = first_final
        .ok_or_else(|| "BRepAlgo ConcatenateWireC0: null first vertex".to_string())?
        .0
        .oriented(Orientation::Forward);
    let last_final = last_final
        .ok_or_else(|| "BRepAlgo ConcatenateWireC0: null last vertex".to_string())?
        .0
        .oriented(Orientation::Reversed);

    // `cxx:682-685`.
    if curve_seq.is_empty() {
        return Ok(None);
    }

    // `cxx:687-688`.
    let nb_curve = curve_seq.len();
    let mut tab: Vec<GeomBSplineCurve> = Vec::with_capacity(nb_curve);
    let mut tab_tol_vertex: Vec<f64> = vec![0.0; nb_curve];

    let res_edge: Option<Edge>;

    if nb_curve > 1 {
        for i in 1..=nb_curve {
            // `cxx:699-702`: kill a `Geom_TrimmedCurve` basis first.
            if curve_seq[i - 1].is_geom_trimmed() {
                let (basis, _, _) = brep_lib_make_edge::adaptor_curve(&curve_seq[i - 1]);
                curve_seq[i - 1] = basis;
            }

            // `cxx:704-705`: `Geom_TrimmedCurve(aCurve, FparSeq(i), LparSeq(i))`
            // then `GeomConvert::CurveToBSplineCurve`.
            let trimmed = GeomTrimmedCurveBasis::with_sense(
                Arc::clone(&curve_seq[i - 1]),
                fpar_seq[i - 1],
                lpar_seq[i - 1],
                true,
            )
            .map_err(|e| format!("BRepAlgo ConcatenateWireC0: Geom_TrimmedCurve ({e})"))?;
            let bs = curve_to_bspline_curve(&trimmed, ParameterisationType::TgtThetaOver2)
                .map_err(|e| format!("BRepAlgo ConcatenateWireC0: CurveToBSplineCurve ({e:?})"))?;
            // `cxx:706`: `GeomConvert::C0BSplineToC1BSplineCurve(.., Confusion)`.
            let mut bs = c0_bspline_to_c1_bspline_curve(&bs, CONFUSION, ANGULAR)
                .map_err(|e| format!("BRepAlgo ConcatenateWireC0: C0BSplineToC1 ({e})"))?;

            // `cxx:708-711`.
            if !is_fwd_seq[i - 1] {
                bs.reverse();
            }

            // `cxx:719-722`: the junction tolerance of every entry but the
            // first, i.e. `TolSeq(i - 1) * 5`.
            if i > 1 {
                tab_tol_vertex[i - 2] = tol_seq[i - 2] * 5.0;
            }
            tab.push(bs);
        }
        // `cxx:723`.
        tab_tol_vertex[nb_curve - 1] = *tol_seq
            .last()
            .expect("one tolerance is appended after the walk")
            * 5.0;

        // `cxx:725-741`: the closed flag, gated on the two end vertices being
        // the same and the two end curves meeting at `GeomAbs_G1` or better.
        let mut closed_flag = false;
        let mut closed_tolerance = 0.0f64;
        let same_ends = match (&first_vertex, &last_vertex) {
            (Some(a), Some(b)) => is_same(&a.0, &b.0),
            _ => false,
        };
        if same_ends {
            let cont = continuity(
                &tab[0],
                &tab[nb_curve - 1],
                tab[0].first_parameter(),
                tab[nb_curve - 1].last_parameter(),
                false,
                false,
                lin_tol,
                ang_tol,
            )
            .map_err(|e| format!("BRepAlgo ConcatenateWireC0: {e}"))?;
            if cont >= Shape::G1 {
                closed_flag = true;
                closed_tolerance = BRepTool::vertex_tolerance(
                    &first_vertex
                        .clone()
                        .expect("the shared end vertex exists"),
                );
            }
        }

        // `cxx:748-751`.
        let (_indices, mut concatenated) =
            concat_c1(&mut tab, &tab_tol_vertex, &mut closed_flag, closed_tolerance)
                .map_err(|e| format!("BRepAlgo ConcatenateWireC0: ConcatC1 ({e})"))?;

        // `cxx:753-773`: more than one concatenated curve is folded back into
        // one by `GeomConvert_CompCurveToBSplineCurve` with
        // `MaxTolVer = max(TolSeq) * 5`.
        if concatenated.len() > 1 {
            let max_tol_ver = tol_seq.iter().cloned().fold(lin_tol, f64::max) * 5.0;
            let first_curve = concatenated
                .first()
                .ok_or_else(|| "BRepAlgo ConcatenateWireC0: empty ConcatC1 result".to_string())?;
            let mut concat = CompCurveToBSplineCurve::from_curve(
                first_curve,
                ParameterisationType::TgtThetaOver2,
            )
            .map_err(|e| format!("BRepAlgo ConcatenateWireC0: CompCurve ({e:?})"))?;
            for curve in concatenated.iter().skip(1) {
                // `cxx:769`: the `.cxx` ignores the `Add` result.
                let _ = concat.add(curve, max_tol_ver, true, true, 0);
            }
            concatenated[0] = concat
                .into_curve()
                .ok_or_else(|| "BRepAlgo ConcatenateWireC0: null composed curve".to_string())?;
        }

        // `cxx:775-778`.
        if is_reverse {
            concatenated
                .first_mut()
                .ok_or_else(|| "BRepAlgo ConcatenateWireC0: empty ConcatC1 result".to_string())?
                .reverse();
        }

        // `cxx:779-786`.
        let curve: Arc<dyn Curve> = Arc::new(
            concatenated
                .first()
                .ok_or_else(|| "BRepAlgo ConcatenateWireC0: empty ConcatC1 result".to_string())?
                .clone(),
        );
        let build = brep_lib_make_edge::init_3d_vertices(
            &curve,
            Some(&Vertex(first_final.clone())),
            Some(&Vertex(last_final.clone())),
            curve.first_parameter(),
            curve.last_parameter(),
        )
        .ok_or_else(|| "BRepAlgo ConcatenateWireC0: BRepLib_MakeEdge not done".to_string())?;
        let builder = TopoBuilder::new();
        let mut edge = builder.make_edge(Arc::clone(&curve), build.first, build.last);
        brep_lib_make_edge::add_endpoint_vertices(&builder, &mut edge, &build);
        res_edge = Some(edge);
    } else {
        // `cxx:787-793`.
        if curve_seq[0].is_geom_trimmed() {
            let (basis, _, _) = brep_lib_make_edge::adaptor_curve(&curve_seq[0]);
            curve_seq[0] = basis;
        }
        // `cxx:791`: `CurveSeq(1)->Copy()`.
        let a_copy_curve: Arc<dyn Curve> = Arc::from(curve_seq[0].clone_dyn());
        let build = brep_lib_make_edge::init_3d_vertices(
            &a_copy_curve,
            Some(&Vertex(first_final.clone())),
            Some(&Vertex(last_final.clone())),
            fpar_seq[0],
            lpar_seq[0],
        )
        .ok_or_else(|| "BRepAlgo ConcatenateWireC0: BRepLib_MakeEdge not done".to_string())?;
        let builder = TopoBuilder::new();
        let mut edge = builder.make_edge(Arc::clone(&a_copy_curve), build.first, build.last);
        brep_lib_make_edge::add_endpoint_vertices(&builder, &mut edge, &build);
        res_edge = Some(edge);
    }

    // `cxx:795-800`.
    let mut res_edge = res_edge
        .ok_or_else(|| "BRepAlgo ConcatenateWireC0: empty result edge".to_string())?;
    if is_reverse {
        res_edge.0.reverse();
    }
    Ok(Some(res_edge))
}
