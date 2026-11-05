//! `Init` (`BRepFill_OffsetWire.cxx:327-366`) and the bisecting-locus setup it
//! performs.
//!
//! `Init` runs once per offset wire: it copies the spine as a FORWARD face,
//! replaces it by the explorer's rebuilt face (`BRepMAT2d_Explorer::Perform`,
//! which re-adds the wires with the p-curves cut at their curvature extrema),
//! prepares the work spine ([`BRepFillOffsetWire::prepare_spine`]), takes the
//! `KPartCircle` shortcut when it applies (`cxx:352-355`) and otherwise computes
//! the left bisecting locus plus the topology link (`cxx:362-365`).

use occt_core::kernel::geomabs::JoinType;

use crate::abs::Orientation;
use crate::brep_fill_offset_wire::BRepFillOffsetWire;
use crate::brep_mat2d_explorer::BRepMat2dExplorer;
use crate::mat::MatSide;
use crate::shape::Face;

impl BRepFillOffsetWire {
    /// `BRepFill_OffsetWire(Spine, Join, IsOpenResult)` (`cxx:318-322`).
    pub fn with_face(spine: &Face, join: JoinType, is_open_result: bool) -> Self {
        let mut oper = Self::new();
        oper.init(spine, join, is_open_result);
        oper
    }

    /// `Init(Spine, Join, IsOpenResult)` (`cxx:327-366`).
    pub fn init(&mut self, spine: &Face, join: JoinType, is_open_result: bool) {
        self.my_is_done = false;
        let local = spine.0.oriented(Orientation::Forward);
        self.my_spine = Some(Face(local.clone()));
        self.my_join_type = join;
        self.my_is_open_result = is_open_result;

        self.my_map.clear();
        self.my_map_spine.clear();

        // `cxx:345-348`: the spine is replaced by the explorer's rebuilt face.
        let mut exp = BRepMat2dExplorer::new();
        exp.perform(&Face(local.clone()));
        let modified = exp.modified_shape(&local);
        self.my_spine = Some(Face(modified));

        // `cxx:348`.
        self.prepare_spine();

        // `cxx:350-356`: the single-closed-circle shortcut.
        if self.kpart_circle(1.0, 0.0) {
            return;
        }

        // `cxx:362-365`.
        let work = self
            .my_work_spine
            .clone()
            .expect("BRepFill_OffsetWire::Init: null work spine");
        exp.perform(&work);
        self.my_bilo.compute(
            &mut exp,
            1,
            MatSide::Left,
            self.my_join_type,
            self.my_is_open_result,
        );
        self.my_link.perform(&exp, &self.my_bilo);
    }
}

