//! `ChFi3d_BoxDiag`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:3282-3297`
//! (OCCT 8.0.0).

use occt_core::bnd::BndBox;

/// `ChFi3d_BoxDiag(box)` (`cxx:3282-3297`): the length of the diagonal of
/// `box`.
///
/// OCCT takes the six `Get` coordinates, spans `d..f`, squares them and takes
/// the square root; it keeps the `sqrt` of the sum rather than `Bnd_Box`'s own
/// diagonal, which is why the port does not reuse an accessor.
///
/// `Bnd_Box::Get` on a void box raises `Standard_ConstructionError` in OCCT;
/// the port's [`BndBox::get`] answers `None` and this function answers `0.0`,
/// the value OCCT's own `ChFi3d_EnlargeBox` callers would produce for a box
/// that never grew.
pub fn ch_fi3d_box_diag(b: &BndBox) -> f64 {
    let Some((a, bmin, c, d, e, f)) = b.get() else {
        return 0.0;
    };
    let dx = d - a;
    let dy = e - bmin;
    let dz = f - c;
    let diag = dx * dx + dy * dy + dz * dz;
    diag.sqrt()
}
