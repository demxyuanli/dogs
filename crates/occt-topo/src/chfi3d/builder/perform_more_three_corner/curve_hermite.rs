//! `CurveHermite` - the Hermite bridge between two stripes of a corner.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx:283-527`
//! (OCCT 8.0.0).
//!
//! The function evaluates the two end curves at the corner parameters, builds
//! the cubic Hermite segment between them (`PLib::HermiteCoefficients` +
//! `PLib::CoefficientsPoles`), makes a `Geom_BezierCurve` edge out of it and
//! projects that edge on every face that separates the two stripes
//! (`BRepAlgo_NormalProjection`, `SetParams(1.e-4, 1.e-4, GeomAbs_C1, 14, 16)`,
//! **without** `SetLimit(false)`).
//!
//! Translation notes:
//! - the `NCollection_Sequence` out-parameters (`proj2d`, `cproj`, `Eproj`,
//!   `param`) are plain vectors, appended in the same order; `Eproj` holds
//!   `Option<Edge>` because the `MapE1.Extent() == 0` arms append the
//!   loop-local null `E1` (`cxx:459-471`),
//! - `BRep_Tool::Curve(E, u1, u2)` is [`BRepTool::edge_curve`] plus
//!   [`BRepTool::edge_parameters`]: the OCCT wrappers of `cxx:445-455` are
//!   `Geom2d_TrimmedCurve(pcurve, up1, up2)` (the p-curve range) and
//!   `Geom_TrimmedCurve(curve, up1, up2)` where `up1`/`up2` have just been
//!   **overwritten** by the edge range, so the two wrappers do not share a
//!   range,
//! - `Extrema_ExtCC(C, L)` is [`GGenExtCC`]: `L` is the Bezier, so OCCT's
//!   analytic `Extrema_ExtElC` dispatch (`Extrema_ExtCC.cxx:251-305`) never
//!   applies and the general `Extrema_ECC` path is the faithful one.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpVec};
use occt_core::kernel::geomabs::Shape as GeomAbsShape;
use occt_core::precision::CONFUSION;
use occt_geom::bezier_curve::GeomBezierCurve;
use occt_geom::bspl_tang_extend::coefficients_poles_dim;
use occt_geom::extrema_cc::GGenExtCC;
use occt_geom::extrema_pc::extrema_ext_pc_range;
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_math::hermite_coefficients;

use crate::brep_lib_make_edge;
use crate::brep_lib_make_face::make_face_from_surface;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::boptools_2d::curve_on_surface_range;
use crate::normal_projection::NormalProjection;
use crate::shape::{Edge, Face};
use crate::topo_tools::edges_of;
use crate::topopebrep::ds::DataStructure;
use crate::topopebrep::tool::adaptor::GeomAdaptorCurve;

/// One end of the Hermite bridge: the argument block of `CurveHermite` for
/// either side (`cxx:285-297` / `:290-296`).
pub(super) struct HermiteSide<'a> {
    /// `CDicmoins` / `CDicplus`
    pub stripe: &'a crate::chfi_ds::StripeHandle,
    /// `jficmoins` / `3 - jficplus` - the `InterferenceOnS1`/`S2` selector of
    /// `cxx:318-328` (left) and `cxx:336-345` (right, which flips the side).
    pub jfp: i32,
    /// `icmoins` / `icplus` - a 1-based `SetOfSurfData` index.
    pub index: i32,
    /// `picmoins` / `picplus` - the parameter on the end curve.
    pub par: f64,
    /// `sensicmoins` / `sensicplus`
    pub sens: i32,
    /// `sharpicmoins` / `sharpicplus`
    pub sharp: bool,
    /// `Eviveicmoins` / `Eviveicplus`
    pub evive: &'a Edge,
}

/// `CurveHermite(...)` (`cxx:283-527`).
#[allow(clippy::too_many_arguments)]
pub(super) fn curve_hermite(
    ds: &DataStructure,
    left: &HermiteSide,
    right: &HermiteSide,
    nbface: usize,
    ecom: &[Edge],
    face: &[Option<Face>],
    proj2d: &mut Vec<Option<Arc<dyn Curve2d>>>,
    cproj: &mut Vec<Option<Arc<dyn Curve>>>,
    eproj: &mut Vec<Option<Edge>>,
    param: &mut Vec<f64>,
    error: &mut f64,
) -> Result<(), String> {
    // `cxx:315-347`: the 3D curve of each end - the edge's own curve for a
    // sharp end, the DS curve of the interference otherwise.
    let c1 = end_curve(ds, left)?;
    let c2 = end_curve(ds, right)?;

    // `cxx:354-355`.
    let (p01, d11) = c1.d1(left.par);
    let (p02, d12) = c2.d1(right.par);

    // `cxx:356-359`.
    let size = 4;
    let mat = hermite_coefficients(0.0, 1.0, 1, 1)
        .ok_or_else(|| "CurveHermite: PLib::HermiteCoefficients failed".to_string())?;

    // `cxx:360-381`: the four Hermite coefficients. `Cont(2)` is negated when
    // the end points *into* the bridge, `Cont(4)` when it points *away*.
    let l1 = p01.distance(&p02);
    let mut cont = [GpVec::zero(); 4];
    let lambda = 1.0 / (d11.magnitude() / l1).max(1.0e-6);
    cont[0] = GpVec::from_xyz(p01.xyz());
    cont[1] = if left.sens == 1 {
        d11.multiplied_scalar(-lambda)
    } else {
        d11.multiplied_scalar(lambda)
    };
    let lambda = 1.0 / (d12.magnitude() / l1).max(1.0e-6);
    cont[2] = GpVec::from_xyz(p02.xyz());
    cont[3] = if right.sens == 1 {
        d12.multiplied_scalar(lambda)
    } else {
        d12.multiplied_scalar(-lambda)
    };

    // `cxx:384-391`: `ExtraCoeffs(jj) = sum_ii MatCoefs(ii, jj) * Cont(ii)`.
    let mut extra = [GpVec::zero(); 4];
    for ii in 0..size {
        for jj in 0..size {
            extra[jj] = extra[jj].added(&cont[ii].multiplied_scalar(mat[ii][jj]));
        }
    }

    // `cxx:392`: `PLib::CoefficientsPoles(ExtraCoeffs, PLib::NoWeights(),
    // ExtrapPoles, PLib::NoWeights())`.
    let mut flat_coefs = Vec::with_capacity(size * 3);
    for v in extra.iter() {
        let q = v.xyz();
        flat_coefs.push(q.x);
        flat_coefs.push(q.y);
        flat_coefs.push(q.z);
    }
    let mut flat_poles = vec![0.0; size * 3];
    coefficients_poles_dim(3, &flat_coefs, &mut flat_poles);
    let poles = (0..size)
        .map(|i| {
            GpPnt::new(
                flat_poles[i * 3],
                flat_poles[i * 3 + 1],
                flat_poles[i * 3 + 2],
            )
        })
        .collect::<Vec<_>>();

    // `cxx:393-395`: `new Geom_BezierCurve(ExtrapPoles)` +
    // `BRepLib_MakeEdge Bedge(Bezier)`.
    let bezier = GeomBezierCurve::new(poles).map_err(|e| format!("CurveHermite: {e}"))?;
    let bezier: Arc<dyn Curve> = Arc::new(bezier);
    let builder = TopoBuilder::new();
    let mut edg = builder.make_edge(
        bezier.clone(),
        bezier.first_parameter(),
        bezier.last_parameter(),
    );
    if let Some(build) = brep_lib_make_edge::init_3d(
        &bezier,
        bezier.first_parameter(),
        bezier.last_parameter(),
    ) {
        brep_lib_make_edge::add_endpoint_vertices(&builder, &mut edg, &build);
    }

    // `cxx:398-474`: one projection per separating face.
    *error = 1.0e-30;
    for nb in 1..=nbface {
        // A null `Face.Value(nb)` makes OCCT's `Init` / `Build` fail, which
        // appends nothing to any of the four sequences - the same effect as
        // skipping the face here.
        let Some(mut f) = face.get(nb - 1).cloned().flatten() else {
            continue;
        };
        // `cxx:404-410`.
        let mut proj = NormalProjection::new();
        proj.init(f.0.clone());
        proj.add(&edg.0);
        proj.set_params(1.0e-4, 1.0e-4, GeomAbsShape::C1 as u8, 14, 16);
        let _ = proj.build();

        if proj.is_done() {
            let mut map_e1 = edges_of(proj.projection());
            if !map_e1.is_empty() {
                // `cxx:415-428`: more than one projected edge - rebuild the
                // face on its own surface and project again.
                if map_e1.len() != 1 {
                    if let Some(surface) = BRepTool::face_surface(&f) {
                        let bface = make_face_from_surface(surface, CONFUSION);
                        f = bface;
                    }
                    proj.init(f.0.clone());
                    let _ = proj.build();
                    map_e1 = if proj.is_done() {
                        edges_of(proj.projection())
                    } else {
                        Vec::new()
                    };
                }
                if !map_e1.is_empty() {
                    // `cxx:430-439`: the first non-degenerated edge. `E1`
                    // keeps the first entry when every edge is degenerated
                    // (`trouve` never set, the loop still assigned it).
                    let mut e1 = map_e1[0].clone();
                    let mut trouve = false;
                    for candidate in map_e1.iter() {
                        e1 = candidate.clone();
                        if !BRepTool::is_degenerated(&e1) {
                            trouve = true;
                            break;
                        }
                    }
                    let _ = trouve;

                    // `cxx:441`.
                    eproj.push(Some(e1.clone()));

                    // `cxx:442-447`.
                    let (proj1, up1, up2) = curve_on_surface_range(&e1, &f).ok_or_else(|| {
                        "Failed to get p-curve of edge".to_string()
                    })?;
                    proj2d.push(Some(Arc::new(Geom2dTrimmedCurve::new(
                        proj1, up1, up2,
                    )) as Arc<dyn Curve2d>));

                    // `cxx:448-455`: `BRep_Tool::Curve(E1, up1, up2)` overwrites
                    // the p-curve range with the edge's own, so the trimmed 3D
                    // curve spans the edge range, not the p-curve one.
                    let proj1c = BRepTool::edge_curve(&e1)
                        .ok_or_else(|| "Failed to get 3D curve of edge".to_string())?;
                    let (c_first, c_last) = BRepTool::edge_parameters(&e1);
                    cproj.push(Some(Arc::new(GeomAdaptorCurve::with_range(
                        proj1c, c_first, c_last,
                    )) as Arc<dyn Curve>));
                    let tol = BRepTool::edge_tolerance(&e1);
                    if *error > tol {
                        *error = tol;
                    }
                } else {
                    // `cxx:459-465`.
                    eproj.push(None);
                    proj2d.push(None);
                    cproj.push(None);
                }
            } else {
                // `cxx:466-471`.
                eproj.push(None);
                proj2d.push(None);
                cproj.push(None);
            }
        }
    }

    // `cxx:475-526`: move the shared parameters onto the Bezier. `p01` is the
    // function-scope point of `cxx:308` (written by `c1->D1` above) and is
    // carried across the iterations, exactly as in the source.
    let mut p01 = p01;
    for nb in 1..=nbface.saturating_sub(1) {
        let edge = ecom
            .get(nb - 1)
            .ok_or_else(|| format!("CurveHermite: Ecom({nb}) is missing"))?;
        let c_curve = BRepTool::edge_curve(edge)
            .ok_or_else(|| "Failed to get 3D curve of edge".to_string())?;
        let (c_first, c_last) = BRepTool::edge_parameters(edge);
        let c = GeomAdaptorCurve::with_range(c_curve, c_first, c_last);

        // `cxx:477-478`.
        let par_nb = *param
            .get(nb - 1)
            .ok_or_else(|| format!("CurveHermite: param({nb}) is missing"))?;
        let p02 = c.d0(par_nb);

        // `cxx:479-481`.
        let mut ext = GGenExtCC::new(
            &c as &dyn Curve,
            bezier.as_ref(),
            (c_first, c_last),
            (bezier.first_parameter(), bezier.last_parameter()),
        );
        let _ = ext.perform();

        // `cxx:482-503`.
        if ext.is_done() && !ext.is_parallel() && ext.nb_ext() != 0 {
            let (u, p_on_c, _v, p_on_l) = ext.points(1);
            if p_on_c.distance(&p_on_l) < CONFUSION {
                param[nb - 1] = u;
            } else if let Some(p) = projection_end_point(cproj, nb) {
                p01 = p;
            }
        }

        // `cxx:504-524`. OCCT re-runs this arm whenever the extrema are
        // missing *or* empty, even after the block above.
        if !ext.is_done() || ext.nb_ext() == 0 {
            if let Some(p) = projection_end_point(cproj, nb) {
                p01 = p;
            }
            if p01.distance(&p02) > 1.0e-4 {
                let sols =
                    extrema_ext_pc_range(&c, &p01, c.first_parameter(), c.last_parameter());
                if let Some(s) = sols.first() {
                    param[nb - 1] = s.u;
                }
            }
        }
    }

    Ok(())
}

/// `c1` / `c2` of `cxx:315-347`.
fn end_curve(
    ds: &DataStructure,
    side: &HermiteSide,
) -> Result<Arc<dyn Curve>, String> {
    if side.sharp {
        // `c1 = BRep_Tool::Curve(Eviveicmoins, up1, up2);`
        return BRepTool::edge_curve(side.evive)
            .ok_or_else(|| "Failed to get 3D curve of edge".to_string());
    }
    // `jficmoins == 1 ? InterferenceOnS1() : InterferenceOnS2()` then
    // `DStr.Curve(ilin).Curve()` (`cxx:318-328`, `:333-345`).
    let jf = side.jfp;
    let guard = side
        .stripe
        .read()
        .map_err(|_| "CurveHermite: poisoned stripe".to_string())?;
    let surf_data = guard
        .set_of_surf_data()
        .get(side.index as usize - 1)
        .ok_or_else(|| {
            format!("CurveHermite: SurfData {} is out of range", side.index)
        })?;
    let ilin = if jf == 1 {
        surf_data.interference_on_s1().line_index()
    } else {
        surf_data.interference_on_s2().line_index()
    };
    ds.curve(ilin)
        .curve()
        .cloned()
        .ok_or_else(|| "Failed to get 3D curve of edge".to_string())
}

/// The `cproj.Value(nb)` / `cproj.Value(nb + 1)` end points of
/// `cxx:490-501` and `:505-513`: the last parameter of the previous
/// projection, else the first parameter of the next one.
fn projection_end_point(cproj: &[Option<Arc<dyn Curve>>], nb: usize) -> Option<GpPnt> {
    if let Some(Some(c)) = cproj.get(nb - 1) {
        return Some(c.d0(c.last_parameter()));
    }
    if let Some(Some(c)) = cproj.get(nb) {
        return Some(c.d0(c.first_parameter()));
    }
    None
}
