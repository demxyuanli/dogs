//! `TopOpeBRepBuild_BuilderON::GFillONPartsWES2` — the `yap7` branch.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_BuilderON.cxx:3046-3121`
//! (OCCT 8.0.0).
//!
//! `yap7` is the `hsdFOR` fallback of `!FFinSD(so|do)`. Inside a cut that
//! builds the same state on both sides, the split is only kept when `EG` bounds
//! `FOR` or one of its same-domain faces (`cxx:3062-3089`), and the start
//! orientation of `EspON` is complemented when `FOR` is not the reference face
//! with the same geometric orientation (`cxx:3091-3107`).
//!
//! The `#ifdef OCCT_DEBUG` traces (`cxx:3049-3061`, `:3112-3117`) are not
//! translated.

use crate::abs::ShapeType;
use crate::topexp::Explorer;
use crate::topopebrep::build::wire_edge_set::WireEdgeSet;
use crate::topopebrep::ds::kind::Config;
use crate::topo_tools::is_same;

use super::super::Builder;
use super::wes2::Wes2Ctx;
use super::wes2_fcx::FcxCtx;

impl Wes2Ctx {
    /// `yap7` (`cxx:3046-3121`).
    pub(super) fn branch_yap7(&self, _fcx: &FcxCtx, pb: &Builder, wes: &mut WireEdgeSet) {
        let Some(hds) = pb.data_structure() else {
            return;
        };
        let bds = hds.ds();

        // `cxx:3062-3077` (xpu290598 CTS20212): whether `EG` bounds `FOR` or one
        // of its same-domain faces.
        let mut isbound = false;
        for sd in bds.shape_same_domain(&self.for_) {
            let mut ex = Explorer::new(sd, ShapeType::Edge);
            while ex.more() {
                if is_same(ex.current(), &self.eg) {
                    isbound = true;
                    break;
                }
                ex.next();
            }
        }

        // `cxx:3078-3079`
        let mut newori_e = self.tfe.orientation(self.tb1, ShapeType::Face);

        // `cxx:3081-3089` (CTS20212: tspIN(f7), tspON(e4))
        if self.com_of_cut && !self.eg_bound_for && !isbound {
            return;
        }

        // `cxx:3091-3107` (CTS20212: tspIN(f7), tspON(e5))
        let samegeomref = if self.for_is_ref {
            true
        } else {
            let mut same = self.c_for == Config::SameOriented;
            if self.o_for != self.oref_for {
                same = !same;
            }
            same
        };
        if !samegeomref {
            newori_e = newori_e.complement();
        }

        // `cxx:3109-3118`
        let mut new_e = self.esp_on.clone();
        pb.my_build_tool.set_orientation(&mut new_e, newori_e);
        wes.add_start_element(&new_e);
    }
}
