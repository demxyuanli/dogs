//! `ShapeFix_Face` mode fields and the `ShapeFix_Root::NeedFix` gate they feed.
//!
//! Source: `ShapeFix_Face.hxx:285-295` (fields), `ShapeFix_Face.cxx:130-145`
//! (`ClearModes`), `ShapeFix_Root.hxx:112` / `ShapeFix_Root.lxx:101-104`
//! (`NeedFix(flag, def = true)`).
//!
//! Both read paths leave every one of these keys at its `ClearModes` value:
//!
//! * `ShapeProcess_OperLibrary::fixshape` (`cxx:818-836`) copies each key from
//!   `ShapeProcess` (`FixWireMode` … `FixSplitFaceMode`), whose `DESTEP` /
//!   `DEIGES` defaults are `DE_ShapeFixParameters.hxx:29` `FixMode::FixOrNot`
//!   (`= -1`) — the same as `ClearModes`.
//! * `myAutoCorrectPrecisionMode` and `myFixPeriodicDegenerated` are *not*
//!   among the keys `fixshape` sets, so they keep `ClearModes`' `1` / `-1`
//!   (`ShapeFix_Face.cxx:143-144`). The DE struct's `AutoCorrectPrecisionMode`
//!   (`DE_ShapeFixParameters.hxx:56` `FixMode::Fix`) agrees with that `1`.
//!
//! So [`ShapeFixFaceModes::default`] and [`ShapeFixFaceModes::reader`] are the
//! same set; `reader` exists only to state where the values come from.

use crate::shhealing::wire_fix::{need_fix, FIX, FIX_OR_NOT};

/// The `ShapeFix_Face` mode set (`ShapeFix_Face.hxx:285-295`).
#[derive(Debug, Clone, PartialEq)]
pub struct ShapeFixFaceModes {
    /// `myFixWireMode` (`ShapeFix_Face.cxx:134`): `Perform`'s two wire rounds
    /// (`cxx:365`, `cxx:509`).
    pub fix_wire_mode: i32,
    /// `myFixOrientationMode` (`cxx:135`): `FixOrientation` (`cxx:691`).
    pub fix_orientation_mode: i32,
    /// `myFixAddNaturalBoundMode` (`cxx:136`): read by `isNeedAddNaturalBound`
    /// (`cxx:1125`) under `NeedFix`, so `-1` resolves to false here — the
    /// natural-bound repair is only reached when the mode is explicitly `1`
    /// *and* the surface is U/V periodic.
    pub fix_add_natural_bound_mode: i32,
    /// `myFixMissingSeamMode` (`cxx:137`): `FixMissingSeam` (`cxx:492-498`).
    pub fix_missing_seam_mode: i32,
    /// `myFixSmallAreaWireMode` (`cxx:138`): `FixSmallAreaWire` (`cxx:731`).
    pub fix_small_area_wire_mode: i32,
    /// `myRemoveSmallAreaFaceMode` (`cxx:139`): the `isRemoveSmallFace`
    /// argument of `FixSmallAreaWire` (`cxx:733`).
    pub remove_small_area_face_mode: i32,
    /// `myFixLoopWiresMode` (`cxx:141`): `FixLoopWire` inside the second wire
    /// round (`cxx:583`).
    pub fix_loop_wires_mode: i32,
    /// `myFixIntersectingWiresMode` (`cxx:140`): `FixIntersectingWires`
    /// (`cxx:680`).
    pub fix_intersecting_wires_mode: i32,
    /// `myFixSplitFaceMode` (`cxx:142`): `FixSplitFace` (`cxx:711`).
    pub fix_split_face_mode: i32,
    /// `myAutoCorrectPrecisionMode` (`cxx:143`): the `1.00001` precision
    /// reduction at the head of the first wire round (`cxx:386-397`).
    pub auto_correct_precision_mode: i32,
    /// `myFixPeriodicDegenerated` (`cxx:144`): `FixPeriodicDegenerated`
    /// (`cxx:486`).
    pub fix_periodic_degenerated: i32,
}

impl Default for ShapeFixFaceModes {
    /// `ShapeFix_Face::ClearModes` (`ShapeFix_Face.cxx:130-145`).
    fn default() -> Self {
        ShapeFixFaceModes {
            fix_wire_mode: FIX_OR_NOT,
            fix_orientation_mode: FIX_OR_NOT,
            fix_add_natural_bound_mode: FIX_OR_NOT,
            fix_missing_seam_mode: FIX_OR_NOT,
            fix_small_area_wire_mode: FIX_OR_NOT,
            remove_small_area_face_mode: FIX_OR_NOT,
            fix_loop_wires_mode: FIX_OR_NOT,
            fix_intersecting_wires_mode: FIX_OR_NOT,
            fix_split_face_mode: FIX_OR_NOT,
            auto_correct_precision_mode: FIX,
            fix_periodic_degenerated: FIX_OR_NOT,
        }
    }
}

impl ShapeFixFaceModes {
    /// The mode set `ShapeProcess_OperLibrary::fixshape` produces on the
    /// `FromSTEP.FixShape` / `FromIGES.FixShape` operators
    /// (`ShapeProcess_OperLibrary.cxx:818-836`).
    pub fn reader() -> Self {
        Self::default()
    }

    /// `NeedFix(myFixWireMode)` (`cxx:365`, `cxx:509`).
    pub fn need_fix_wire(&self) -> bool {
        need_fix(self.fix_wire_mode, true)
    }
    /// `NeedFix(myFixOrientationMode)` (`cxx:691`).
    pub fn need_fix_orientation(&self) -> bool {
        need_fix(self.fix_orientation_mode, true)
    }
    /// `NeedFix(myFixAddNaturalBoundMode)` (`cxx:1125`).
    pub fn need_fix_add_natural_bound(&self) -> bool {
        need_fix(self.fix_add_natural_bound_mode, true)
    }
    /// `NeedFix(myFixMissingSeamMode)` — `FixMissingSeam`'s own head is
    /// `if (!NeedFix(myFixMissingSeamMode)) return false;` in effect: the
    /// caller (`cxx:492-498`) calls it unconditionally, and the tool's
    /// `ClearModes` `-1` resolves to true.
    pub fn need_fix_missing_seam(&self) -> bool {
        need_fix(self.fix_missing_seam_mode, true)
    }
    /// `NeedFix(myFixSmallAreaWireMode, false)` (`cxx:731`).
    pub fn need_fix_small_area_wire(&self) -> bool {
        need_fix(self.fix_small_area_wire_mode, false)
    }
    /// `NeedFix(myRemoveSmallAreaFaceMode, false)` (`cxx:733`).
    pub fn remove_small_area_face(&self) -> bool {
        need_fix(self.remove_small_area_face_mode, false)
    }
    /// `NeedFix(myFixLoopWiresMode)` (`cxx:583`).
    pub fn need_fix_loop_wires(&self) -> bool {
        need_fix(self.fix_loop_wires_mode, true)
    }
    /// `NeedFix(myFixIntersectingWiresMode)` (`cxx:680`).
    pub fn need_fix_intersecting_wires(&self) -> bool {
        need_fix(self.fix_intersecting_wires_mode, true)
    }
    /// `NeedFix(myFixSplitFaceMode)` (`cxx:711`).
    pub fn need_fix_split_face(&self) -> bool {
        need_fix(self.fix_split_face_mode, true)
    }
    /// `if (myAutoCorrectPrecisionMode)` (`cxx:387`) — a plain truth test, not
    /// a `NeedFix` call.
    pub fn auto_correct_precision(&self) -> bool {
        self.auto_correct_precision_mode != 0
    }
    /// `NeedFix(myFixPeriodicDegenerated)` (`cxx:486`).
    pub fn need_fix_periodic_degenerated(&self) -> bool {
        need_fix(self.fix_periodic_degenerated, true)
    }

    /// `FixSmallAreaWire`'s `theIsRemoveSmallFace` (`cxx:733`).
    pub fn is_remove_small_face(&self) -> bool {
        self.remove_small_area_face()
    }
}

/// `ShapeFix::LeastEdgeSize` (`ShapeFix.cxx:793-815`): the square root of the
/// smallest 3-point bounding-box diagonal over every edge's 3D curve, sampled
/// at `first`, `(first + last) / 2` and `last`.
///
/// `RealLast()` seeds the running minimum, so a shape with no 3D curve reports
/// `sqrt(RealLast())`, as in OCCT.
pub fn least_edge_size(shape: &crate::shape::TopoShape) -> f64 {
    let mut res = f64::MAX; // `RealLast()`
    for edge in crate::topo_tools::edges_of(shape) {
        let Some(c3d) = crate::brep_tool::BRepTool::edge_curve(&edge) else {
            continue; // cxx:799-803
        };
        let (first, last) = crate::brep_tool::BRepTool::edge_parameters(&edge);
        let mut bb = occt_core::bnd::BndBox::new();
        bb.add_point(&c3d.d0(first));
        bb.add_point(&c3d.d0(last));
        bb.add_point(&c3d.d0((last + first) / 2.0));
        let lo = bb.corner_min();
        let hi = bb.corner_max();
        let size = (hi.x() - lo.x()) * (hi.x() - lo.x())
            + (hi.y() - lo.y()) * (hi.y() - lo.y())
            + (hi.z() - lo.z()) * (hi.z() - lo.z());
        if size < res {
            res = size;
        }
    }
    res.sqrt()
}
