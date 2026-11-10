//! `Intf_Tool::NbSegments` / `BeginParam` / `EndParam`
//! (`Intf_Tool.cxx:1634-1653`).

use super::structure::IntfTool;

impl IntfTool {
    /// `Intf_Tool::NbSegments()` (`Intf_Tool.cxx:1634-1637`).
    pub fn nb_segments(&self) -> usize {
        self.nb_seg
    }

    /// `Intf_Tool::BeginParam(SegmentNum)` (`Intf_Tool.cxx:1641-1645`).
    ///
    /// The argument is 1-based; `Standard_OutOfRange` (`Intf_Tool.cxx:1643`)
    /// is raised as a panic.
    pub fn begin_param(&self, segment_num: usize) -> f64 {
        assert!(
            segment_num >= 1 && segment_num <= self.nb_seg,
            "Intf_Tool::BeginParam"
        );
        self.begin_on_curve[segment_num - 1]
    }

    /// `Intf_Tool::EndParam(SegmentNum)` (`Intf_Tool.cxx:1649-1653`).
    ///
    /// The argument is 1-based; `Standard_OutOfRange` (`Intf_Tool.cxx:1651`)
    /// is raised as a panic.
    pub fn end_param(&self, segment_num: usize) -> f64 {
        assert!(
            segment_num >= 1 && segment_num <= self.nb_seg,
            "Intf_Tool::EndParam"
        );
        self.end_on_curve[segment_num - 1]
    }
}
