//! Port of `Intf_Polygon2d`
//! (`src/ModelingAlgorithms/TKGeomAlgo/Intf/Intf_Polygon2d.hxx/.cxx/.lxx`).
//!
//! OCCT's `Intf_Polygon2d` is an abstract base carrying the bounding box
//! (`myBox`) plus the segmentation query its intersector needs. The port makes
//! it a trait; implementors store their own `BndBox2d` and expose it through
//! [`IntfPolygon2d::bounding`].

use crate::bnd::box2d::BndBox2d;
use crate::gp::GpPnt2d;

/// `Intf_Polygon2d` (`Intf_Polygon2d.hxx:29-56`).
pub trait IntfPolygon2d {
    /// `Bounding()` (`Intf_Polygon2d.lxx:19-22`).
    fn bounding(&self) -> &BndBox2d;

    /// `Closed()` (`Intf_Polygon2d.cxx:20-23`). Default `false`.
    fn closed(&self) -> bool {
        false
    }

    /// `DeflectionOverEstimation()` (`Intf_Polygon2d.hxx:45`).
    fn deflection_over_estimation(&self) -> f64;

    /// `NbSegments()` (`Intf_Polygon2d.hxx:47`).
    fn nb_segments(&self) -> usize;

    /// `Segment(Index, theBegin, theEnd)` (`Intf_Polygon2d.hxx:49-52`).
    /// OCCT uses 1-based indexing.
    fn segment(&self, index: usize, the_begin: &mut GpPnt2d, the_end: &mut GpPnt2d);
}
