//! Curve helpers shared by face and edge unification.
//!
//! Source: `ShapeUpgrade_UnifySameDomain.cxx:268-303`.

use occt_core::gp::{GpAx2d, GpDir2d, GpPnt2d, GpVec2d};
use occt_core::precision::CONFUSION;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

/// `TryMakeLine`: recognizes an exactly unit-parameterized line by the same
/// endpoint-length and ten-sample tests used by OCCT.
pub(super) fn try_make_line(pcurve: &dyn Curve2d, first: f64, last: f64) -> Option<Geom2dLine> {
    let first_point = pcurve.d0(first);
    let last_point = pcurve.d0(last);
    let vector = GpVec2d::new(
        last_point.x() - first_point.x(),
        last_point.y() - first_point.y(),
    );
    let squared_length = vector.square_magnitude();
    let squared_parameter_length = (last - first) * (last - first);
    if (squared_length - squared_parameter_length).abs() > CONFUSION {
        return None;
    }

    let direction = GpDir2d::from_vec2d(&vector).ok()?;
    let origin = GpPnt2d::new(
        first_point.x() - direction.x() * first,
        first_point.y() - direction.y() * first,
    );
    let line = occt_core::gp::GpLin2d::from_pnt_dir(origin, direction);
    let delta = (last - first) / 10.0;
    for sample in 1..10 {
        if line.distance(&pcurve.d0(first + sample as f64 * delta)) > CONFUSION {
            return None;
        }
    }
    Some(Geom2dLine::new(GpAx2d::new(origin, direction)))
}
