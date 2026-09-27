//! `ShapeAnalysis_Wire::CheckSmallArea` (`ShapeAnalysis_Wire.cxx:2004-2098`).

use occt_core::gp::{GpPnt2d, GpVec};
use occt_core::precision::CONFUSION;

use crate::boptools_2d;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Face, Wire};
use crate::topo_tools_full::edges_of_wire;

/// `ShapeAnalysis_Wire::CheckSmallArea(theWire)` (`ShapeAnalysis_Wire.cxx:2004-2098`).
///
/// Returns `true` when the wire encloses a null/small area on `face`:
/// the sampled 3D cross-product area estimator is below `length * Precision`,
/// and the exact `BRepGProp::SurfaceProperties` area of a temporary face
/// carrying only this wire is below `0.5 * LinearProperties.Mass() * Precision`.
///
/// `myPrecision` is `ShapeFix_Root::Precision()` = `Precision::Confusion()`
/// (`ShapeFix_Root.cxx`); the port uses `CONFUSION` for it.
pub fn check_small_area(wire: &Wire, face: &Face) -> bool {
    let Some(surf) = BRepTool::face_surface(face) else {
        return false;
    };
    let edges = edges_of_wire(wire);
    let nb_edges = edges.len();
    if nb_edges < 1 {
        return false;
    }
    const A_NB_CONTROL: usize = 23;
    let an_inv = 1.0 / (A_NB_CONTROL - 1) as f64;

    // cxx:2019-2036: mid point for the closed contour, averaged over every
    // edge's 22 interior samples.
    let mut cx = 0.0f64;
    let mut cy = 0.0f64;
    for e in &edges {
        let Some((c2d, f, l)) = boptools_2d::curve_on_surface_range(e, face) else {
            return false; // cxx:2024-2028 FAIL2
        };
        for i in 1..A_NB_CONTROL {
            let v = an_inv * ((A_NB_CONTROL - 1 - i) as f64 * f + i as f64 * l);
            let p = c2d.d0(v);
            cx += p.x();
            cy += p.y();
        }
    }
    let inv = 1.0 / (nb_edges * (A_NB_CONTROL - 1)) as f64;
    let center2d = GpPnt2d::new(cx * inv, cy * inv);

    // cxx:2038-2077: approximated area in 3D.
    let center = surf.d0(center2d.x(), center2d.y());
    let mut cross = GpVec::new(0.0, 0.0, 0.0);
    let mut length = 0.0f64;
    // cxx:2039-2040: aPnt3d / aPrev3d carry over from one edge to the next.
    let mut pnt3d = occt_core::gp::GpPnt::new(0.0, 0.0, 0.0);
    let mut prev3d = GpVec::new(0.0, 0.0, 0.0);
    for (j, e) in edges.iter().enumerate() {
        let Some(c3d) = BRepTool::edge_curve(e) else {
            return false; // cxx:2046-2051 FAIL2
        };
        let (f, l) = BRepTool::edge_parameters(e);
        if !(f.is_finite() && l.is_finite()) {
            continue; // cxx:2052-2055
        }
        let begin = if j == 0 {
            let p = c3d.d0(f);
            pnt3d = p;
            prev3d = GpVec::from_pnts(&center, &p);
            1
        } else {
            0
        };
        for i in begin..A_NB_CONTROL {
            let u = an_inv * ((A_NB_CONTROL - 1 - i) as f64 * f + i as f64 * l);
            let pnt = c3d.d0(u);
            let v = GpVec::from_pnts(&center, &pnt);
            cross = cross.added(&prev3d.crossed(&v));
            length += pnt3d.distance(&pnt);
            pnt3d = pnt;
            prev3d = v;
        }
    }

    // cxx:2079-2096: exact check when the estimator says the area is small.
    let tolerance = length * CONFUSION;
    if cross.magnitude() < tolerance {
        let builder = TopoBuilder::new();
        let a_face = builder.make_face(surf, &[wire.clone()]); // cxx:2085-2086
        if let (Ok((sp, _)), Ok(lp)) = (
            crate::brep_gprop_full::surface_properties(&a_face.0),
            crate::brep_gprop_full::linear_properties(&a_face.0),
        ) {
            let new_tolerance = lp.mass() * CONFUSION;
            if sp.mass().abs() < 0.5 * new_tolerance {
                return true; // cxx:2091-2095 DONE1
            }
        }
    }
    false
}
