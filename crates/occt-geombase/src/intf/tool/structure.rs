//! `Intf_Tool` state (`Intf_Tool.hxx:35-85`) and the constructor
//! (`Intf_Tool.cxx:38-48`).
/// `Intf_Tool` (`Intf_Tool.hxx:35-85`).
#[derive(Clone, Debug)]
pub struct IntfTool {
    /// `nbSeg` (`Intf_Tool.hxx:77`).
    pub(super) nb_seg: usize,
    /// `beginOnCurve[6]` (`Intf_Tool.hxx:78`).
    pub(super) begin_on_curve: [f64; 6],
    /// `endOnCurve[6]` (`Intf_Tool.hxx:79`).
    pub(super) end_on_curve: [f64; 6],
    /// `bord[12]` (`Intf_Tool.hxx:80`).
    pub(super) bord: [i32; 12],
    /// `xint[12]` (`Intf_Tool.hxx:81`).
    pub(super) xint: [f64; 12],
    /// `yint[12]` (`Intf_Tool.hxx:82`).
    pub(super) yint: [f64; 12],
    /// `zint[12]` (`Intf_Tool.hxx:83`). Written by the 3D conic paths
    /// only, which are UNPORTED (see the module header).
    #[allow(dead_code)]
    pub(super) zint: [f64; 12],
    /// `parint[12]` (`Intf_Tool.hxx:84`).
    pub(super) parint: [f64; 12],
}

impl Default for IntfTool {
    fn default() -> Self {
        Self::new()
    }
}

impl IntfTool {
    /// `Intf_Tool::Intf_Tool()` (`Intf_Tool.cxx:38-48`).
    pub fn new() -> Self {
        Self {
            nb_seg: 0,
            begin_on_curve: [0.0; 6],
            end_on_curve: [0.0; 6],
            bord: [0; 12],
            xint: [0.0; 12],
            yint: [0.0; 12],
            zint: [0.0; 12],
            parint: [0.0; 12],
        }
    }
}
