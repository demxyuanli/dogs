//! The three `Perform` overloads of `LocOpe_CSIntersector`
//! (`LocOpe_CSIntersector.cxx:69-177`).

use std::sync::Arc;

use occt_core::gp::{GpCirc, GpLin};
use occt_core::precision::PCONFUSION;
use occt_geom::{Curve, GeomCircle};

use crate::abs::ShapeType;
use crate::int_curves_face::FaceIntersector;
use crate::shape::Face;
use crate::topexp::Explorer;

use super::add_points::add_points;
use super::CSIntersector;

/// The faces of `myShape`, in `TopExp_Explorer(S, TopAbs_FACE)` order
/// (`cxx:84`).
fn faces(m: &CSIntersector) -> Vec<Face> {
    let mut out = Vec::new();
    let mut ex = Explorer::new(&m.my_shape, ShapeType::Face);
    while ex.more() {
        if let Some(f) = Face::wrap(ex.current().clone()) {
            out.push(f);
        }
        ex.next();
    }
    out
}

/// `Perform(Slin)` (`cxx:69-98`).
pub(super) fn perform_lines(m: &mut CSIntersector, lines: &[GpLin]) -> Result<(), String> {
    if m.my_shape.is_null() || lines.is_empty() {
        return Err("Standard_ConstructionError".to_string());
    }
    m.my_done = false;
    m.my_points = vec![Vec::new(); lines.len()];

    // `cxx:84-85`: `binf = RealFirst()`, `bsup = RealLast()`.
    let binf = -f64::MAX;
    let bsup = f64::MAX;

    for theface in faces(m) {
        let mut the_int = FaceIntersector::new(theface.clone(), PCONFUSION, true, true);
        for (i, line) in lines.iter().enumerate() {
            the_int.perform(line, binf, bsup);
            if the_int.is_done() {
                add_points(&the_int, &mut m.my_points[i], &theface);
            }
        }
    }
    m.my_done = true;
    Ok(())
}

/// `Perform(Scir)` (`cxx:102-137`): every circle is loaded into one
/// `GeomAdaptor_Curve` over `[0, 2*pi]`, and the face intersector is built with
/// tolerance `0.` (`cxx:125`).
pub(super) fn perform_circles(m: &mut CSIntersector, circles: &[GpCirc]) -> Result<(), String> {
    if m.my_shape.is_null() || circles.is_empty() {
        return Err("Standard_ConstructionError".to_string());
    }
    m.my_done = false;
    m.my_points = vec![Vec::new(); circles.len()];

    let binf = 0.0f64;
    let bsup = 2.0 * std::f64::consts::PI;

    for theface in faces(m) {
        let mut the_int = FaceIntersector::new(theface.clone(), 0.0, true, true);
        for (i, circ) in circles.iter().enumerate() {
            let curve: Arc<dyn Curve> = Arc::new(GeomCircle::new(circ.clone()));
            the_int.perform_curve(curve.as_ref(), binf, bsup);
            if the_int.is_done() {
                add_points(&the_int, &mut m.my_points[i], &theface);
            }
        }
    }
    m.my_done = true;
    Ok(())
}

/// `Perform(Scur)` (`cxx:140-177`): the range of every element is its adaptor's
/// own `FirstParameter`/`LastParameter` (`cxx:168-169`); a null handle is
/// skipped.
pub(super) fn perform_curves(
    m: &mut CSIntersector,
    curves: &[Option<Arc<dyn Curve>>],
) -> Result<(), String> {
    if m.my_shape.is_null() || curves.is_empty() {
        return Err("Standard_ConstructionError".to_string());
    }
    m.my_done = false;
    m.my_points = vec![Vec::new(); curves.len()];

    for theface in faces(m) {
        let mut the_int = FaceIntersector::new(theface.clone(), 0.0, true, true);
        for (i, curve) in curves.iter().enumerate() {
            let Some(curve) = curve else {
                continue;
            };
            let binf = curve.first_parameter();
            let bsup = curve.last_parameter();
            the_int.perform_curve(curve.as_ref(), binf, bsup);
            if the_int.is_done() {
                add_points(&the_int, &mut m.my_points[i], &theface);
            }
        }
    }
    m.my_done = true;
    Ok(())
}
