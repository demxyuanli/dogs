//! `TopOpeBRepBuild_BuilderON::GFillONPartsWES2` — the `yap4` branch.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_BuilderON.cxx:3013-3043`
//! (OCCT 8.0.0).
//!
//! `yap4` is the unconditional fallback of `!FFinSD(so|do)`: `EspON` is added
//! with the orientation of the transition of `FOR` against `TB1`, without any
//! further test. `yap3b` and the `FFinSD` false arms can both reach it
//! (`cxx:1769-1771`, `:1757-1760`).
//!
//! The `#ifdef OCCT_DEBUG` traces (`cxx:3015-3027`, `:3032-3038`) are not
//! translated.

use crate::abs::ShapeType;
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;

use super::super::Builder;
use super::wes2::Wes2Ctx;
use super::wes2_fcx::FcxCtx;

impl Wes2Ctx {
    /// `yap4` (`cxx:3013-3043`).
    pub(super) fn branch_yap4(&self, _fcx: &FcxCtx, pb: &Builder, wes: &mut WireEdgeSet) {
        // `cxx:3029-3031`
        let otfe = self.tfe.orientation(self.tb1, ShapeType::Face);
        let newori_e = otfe;

        // `cxx:3031-3041`
        let mut new_e = self.esp_on.clone();
        pb.my_build_tool.set_orientation(&mut new_e, newori_e);
        wes.add_start_element(&new_e);
    }
}
