//! `TopOpeBRep_PointGeomTool` — the filler-side factory turning the three
//! intersection-point representations into a `TopOpeBRepDS_Point`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_PointGeomTool.hxx:30-44`,
//! `.cxx:29-98`.
//!
//! Ported: `MakePoint(const TopOpeBRep_Point2d&)` (`cxx:73-77`),
//! `MakePoint(const TopoDS_Shape&)` (`cxx:87-90`),
//! `MakePoint(const TopOpeBRep_FaceEdgeIntersector&)` (`cxx:79-84`),
//! `MakePoint(const TopOpeBRep_VPointInter&)` (`cxx:29-67`) and
//! `IsEqual(DSP1, DSP2)` (`cxx:94-98`).

use crate::shape::{Edge, TopoShape};
use crate::topopebrep::ds::Point;
use crate::topopebrep::face_edge_intersector::FaceEdgeIntersector;
use crate::topopebrep::interf::v_point_inter::VPointInter;

use super::point2d::Point2d;

/// `class TopOpeBRep_PointGeomTool` (`hxx:30-44`).
pub struct PointGeomTool;

impl PointGeomTool {
    /// `MakePoint(const TopOpeBRep_VPointInter& IP)` (`cxx:29-67`): the
    /// tolerance is widened by the edge tolerances the point lies on.
    ///
    /// The `si == 3` arm (`cxx:47-64`) caps `t1` at `0.9` (the second cap is
    /// a copy of the first, kept as OCCT writes it) and then uses
    /// `2.5 * (t1 + t2)`.
    pub fn make_point_of_v_point_inter(ip: &VPointInter) -> Result<Point, String> {
        let tolip = ip.tolerance()?;
        let mut tolout = tolip;
        let si = ip.shape_index();
        if si == 0 {
            tolout = tolip;
        } else if si == 1 || si == 2 {
            let s = ip.edge(si).ok_or_else(|| {
                "TopOpeBRep_PointGeomTool::MakePoint : no edge".to_string()
            })?;
            tolout = crate::brep_tool::BRepTool::edge_tolerance(&Edge(s.clone()));
        } else if si == 3 {
            let e1 = ip.edge(1).ok_or_else(|| {
                "TopOpeBRep_PointGeomTool::MakePoint : no edge".to_string()
            })?;
            let e2 = ip.edge(2).ok_or_else(|| {
                "TopOpeBRep_PointGeomTool::MakePoint : no edge".to_string()
            })?;
            let mut t1 = crate::brep_tool::BRepTool::edge_tolerance(&Edge(e1.clone()));
            let t2 = crate::brep_tool::BRepTool::edge_tolerance(&Edge(e2.clone()));
            if t1 > 0.9 {
                t1 = 0.9;
            }
            if t1 > 0.9 {
                t1 = 0.9;
            }
            tolout = 2.5 * (t1 + t2);
        }
        tolout = if tolout > tolip { tolout } else { tolip };
        Ok(Point::with_tolerance(*ip.value()?, tolout))
    }

    /// `MakePoint(const TopOpeBRep_Point2d& P2D)` (`cxx:73-77`).
    pub fn make_point_of_point2d(p2d: &Point2d) -> Point {
        Point::with_tolerance(*p2d.value(), p2d.tolerance())
    }

    /// `MakePoint(const TopOpeBRep_FaceEdgeIntersector&)` (`cxx:78-82`).
    pub fn make_point_of_face_edge_int(feint: &FaceEdgeIntersector) -> Point {
        Point::with_tolerance(feint.value(), feint.tolerance())
    }

    /// `MakePoint(const TopoDS_Shape& S)` (`cxx:87-90`).
    pub fn make_point_of_shape(s: &TopoShape) -> Result<Point, String> {
        Point::of_shape(s)
    }

    /// `IsEqual(DSP1, DSP2)` (`cxx:94-98`).
    pub fn is_equal(dsp1: &Point, dsp2: &Point) -> bool {
        dsp1.is_equal(dsp2)
    }
}
