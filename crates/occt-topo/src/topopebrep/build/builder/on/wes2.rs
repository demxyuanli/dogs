//! `TopOpeBRepBuild_BuilderON::GFillONPartsWES2` — the ON-part router.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_BuilderON.cxx:554-3123`
//! (OCCT 8.0.0).
//!
//! The C++ function reads the whole interference context once (`cxx:556-790`),
//! then runs a chain of mutually exclusive `if (yapXX) { ...; return; }`
//! branches. The port keeps that shape: [`Wes2Ctx`] is the context block and
//! [`Builder::g_fill_on_parts_wes2`] is the chain. Each branch is a method of
//! `Wes2Ctx` in `wes2_yap0.rs` / `wes2_yapc.rs`, so no single file carries the
//! ~2,500 lines of the C++ function.
//!
//! UNPORTED: nothing left in the `yapXX` chain — `yap00`, `yap0`, `yapc1`,
//! `yapc2`, `yapc3`, `yap1`, `yap2`, `yap6`, `yap3`, `yap5`, `yap4` and `yap7`
//! are all ported. `GFillONPartsWES2` is called only by `GFillONPartsWES1`
//! (`cxx:280-289`), which is called by `GFillONPartsWES` (`on.cxx:28-49`) and
//! hence by `GFillFacesWESK` and `GSplitFace` — see [`super`].
//!
//! The `#ifdef OCCT_DEBUG` blocks (`cxx:602-604`, `:648-653`, `:713-731`,
//! `:788-808`, and the `debaddpwes` / `FUN_cout` calls) only print traces and
//! are not translated.

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topopebrep::ds::connex::Connexity;
use crate::topopebrep::ds::export::{
    fds_data, fds_has_same_domain3d, fun_ds_on_esd,
};
use crate::topopebrep::ds::interference::Interference;
use crate::topopebrep::ds::kind::{Config, Kind};
use crate::topopebrep::ds::process_interferences_tool::fun_select_gi_interference;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::tool::closed::is_closing_e_in_shape;
use crate::topopebrep::tool::project::fun_tool_proj_pon_bounded_f;
use crate::topopebrep::tool::tangent::curves_same_oriented;
use crate::topopebrep::tool::{
    fun_tool_bounds, fun_tool_maxtol, fun_tool_value,
};

use super::super::super::g_topo::GTopo;
use super::super::super::wire_edge_set::WireEdgeSet;
use super::super::Builder;
use super::helpers::{fun_build_tb, fun_find_e_sd};

/// The context block of `GFillONPartsWES2` (`cxx:556-790`).
///
/// The fields are the locals of that block, kept 1:1 with the C++. Several of
/// them are read only by the branches that are still UNPORTED
/// (`cxx:1605-3123`), hence the `dead_code` allowance.
#[allow(dead_code)]
pub(super) struct Wes2Ctx {
    pub(super) for_: TopoShape,
    pub(super) esp_on: TopoShape,
    pub(super) eg: TopoShape,
    pub(super) fs: TopoShape,
    pub(super) i_for: i32,
    pub(super) gi: i32,
    pub(super) i_fs: i32,
    pub(super) tb1: State,
    pub(super) tb2: State,
    pub(super) tfe: Transition,
    pub(super) eg_bound_for: bool,
    pub(super) o_for: Orientation,
    pub(super) o_fs: Orientation,
    pub(super) is_closed_ff: bool,
    pub(super) is_closed_fs: bool,
    pub(super) is_closed: bool,
    pub(super) is_rest: bool,
    pub(super) rank_fs: i32,
    pub(super) rank_eg: i32,
    pub(super) rank_for: i32,
    pub(super) otfe: Orientation,
    pub(super) tfe_bef: State,
    pub(super) tfe_aft: State,
    pub(super) eghassd: bool,
    pub(super) eghassd3d: bool,
    /// `le3` (`cxx:605-606`): the 3-D same-domain partners of `EG` collected by
    /// `FDS_HasSameDomain3d`; only the `FCX` tail reads it.
    pub(super) le3: Vec<TopoShape>,
    pub(super) ie3d: i32,
    pub(super) oe3d: Orientation,
    pub(super) ie2d: i32,
    pub(super) o2d_ori: Orientation,
    pub(super) lfcx: Vec<TopoShape>,
    pub(super) hsd_for: bool,
    pub(super) c_for: Config,
    pub(super) iref_for: i32,
    pub(super) oref_for: Orientation,
    pub(super) for_is_ref: bool,
    pub(super) ope_fus: bool,
    pub(super) opec12: bool,
    pub(super) ope_cut: bool,
    pub(super) ope_com: bool,
    pub(super) com_of_cut: bool,
    pub(super) tb_for: State,
    pub(super) tol_a: f64,
    pub(super) tol_eg: f64,
    pub(super) tol_for: f64,
    pub(super) tol_fs: f64,
    /// `factor` (`cxx:652`): the `1.e2` multiplier of the provisional
    /// tolerances; `yap6` also scales the tolerance of `Esd` with it
    /// (`cxx:2521`).
    pub(super) factor: f64,
    pub(super) par_on: f64,
    pub(super) esp_on_esd: bool,
    pub(super) iesd: i32,
    pub(super) eon_fs: bool,
    /// `eONsoEsd` (`cxx:684-686`): whether `EspON` and its same-domain partner
    /// are same-oriented.
    pub(super) e_on_so_esd: bool,
    pub(super) eghassd_on: bool,
    pub(super) eghassd3d_on: bool,
    pub(super) yap00: bool,
    pub(super) yap0: bool,
    pub(super) yap0bis: bool,
    pub(super) yapc1: bool,
    pub(super) yapc2a: bool,
    pub(super) yapc2b: bool,
    pub(super) yapc2c: bool,
    pub(super) yapc2: bool,
    pub(super) yapc3: bool,
    pub(super) hassd3d_on: bool,
    pub(super) e3closed_for: bool,
}

impl Wes2Ctx {
    /// `cxx:556-790`: the context block `GFillONPartsWES2` computes before its
    /// first branch. `None` covers the early `return`s of the block
    /// (`cxx:679-681`, `:690-692`, `:696-698`) — the C++ leaves the function
    /// without adding a start element, the same as a `false` branch.
    fn new(
        b: &Builder,
        for_: &TopoShape,
        g: &GTopo,
        i: &Interference,
        esp_on: &TopoShape,
    ) -> Option<Self> {
        let hds = b.data_structure()?;
        let bds = hds.ds();

        let (tb1, tb2) = g.states_on();
        // `cxx:559-568`
        let (gt, gi, st, si) = fds_data(i);
        if gt != Kind::Edge || st != Kind::Face {
            return None;
        }
        let tfe = i.transition();
        let eg_bound_for = i.g_bound().unwrap_or(false);

        let e_on = esp_on;
        let eg = bds.shape(gi, true)?;
        let fs = bds.shape(si, true)?;
        let i_for = bds.shape_index(for_, true);
        let i_fs = bds.shape_index(&fs, true);

        // `cxx:581` `tola = Precision::Angular() * 1.e3`
        let tol_a = occt_core::precision::ANGULAR * 1.0e3;

        let o_for = bds.shape(i_for, true)?.orientation();
        let o_fs = bds.shape(i_fs, true)?.orientation();

        // `cxx:587-589`
        let is_closed_ff = is_closing_e_in_shape(&eg, for_, for_);
        let is_closed_fs = is_closing_e_in_shape(&eg, &fs, &fs);
        let is_closed = is_closed_ff || is_closed_fs;
        let is_rest = bds.is_section_edge(&eg, false);

        // `cxx:593-595`
        let rank_fs = b.shape_rank(&fs);
        let rank_eg = b.shape_rank(&eg);
        let rank_for = b.shape_rank(for_);
        let otfe = tfe.orientation(State::In, ShapeType::Face);
        let tfe_bef = tfe.before_state();
        let tfe_aft = tfe.after_state();

        // `cxx:604-606`
        let eghassd = hds.has_same_domain(&eg, true);
        let (eghassd3d, le3) = fds_has_same_domain3d(bds, &eg, true);
        let mut ie3d = 0;
        let mut oe3d = Orientation::External;
        let ie2d = 0;
        let o2d_ori = Orientation::External;

        // `cxx:612-614`
        let mut lfcx: Vec<TopoShape> = Vec::new();
        {
            let mut connex = Connexity::new();
            connex.prepare(hds);
            connex.face_edge_connex_faces(&fs, &eg, hds, &mut lfcx);
        }

        // `cxx:616-631`
        let hsd_for = hds.has_same_domain(for_, true);
        let mut c_for = Config::UnshGeometry;
        let mut iref_for = 0;
        let mut oref_for = Orientation::External;
        let mut for_is_ref = false;
        if hsd_for {
            c_for = bds.same_domain_ori(i_for);
            iref_for = bds.same_domain_ref_of(for_);
            if iref_for != 0 {
                oref_for = bds.shape(iref_for, true)?.orientation();
            }
            for_is_ref = iref_for == i_for;
        }

        // `cxx:633-641`
        let ope_fus = b.opefus();
        let opec12 = b.opec12();
        let ope_cut = b.opec12() || b.opec21();
        let ope_com = b.opecom();
        let com_of_cut = ope_cut && (tb1 == tb2) && (tb1 == State::In);

        let tb_for = fun_build_tb(b, rank_for);

        // `cxx:652-658`
        let factor = 1.0e2;
        let tol_eg = BRepTool::edge_tolerance(&Edge(eg.clone())) * factor;
        let tol_for = BRepTool::face_tolerance(&Face(for_.clone())) * factor;
        let tol_fs = BRepTool::face_tolerance(&Face(fs.clone())) * factor;

        // `cxx:660-662`
        let (f, l) = fun_tool_bounds(e_on);
        let x = 0.456189;
        let par_on = (1.0 - x) * f + x * l;

        // `cxx:665-687`: `FUN_ds_ONesd` then the orientation comparison.
        let mut esp_on_esd = false;
        let mut iesd = 0;
        if eghassd {
            esp_on_esd = fun_ds_on_esd(bds, gi, esp_on, &mut iesd);
        }
        let mut eon_so_esd = false;
        if eghassd && iesd != 0 {
            let esd = bds.shape(iesd, true)?;
            eon_so_esd = curves_same_oriented(e_on, &esd)?;
        }

        // `cxx:688-702`
        let mut eon_fs = true;
        if eghassd && !esp_on_esd {
            let p_on = fun_tool_value(par_on, e_on)?;
            let (_u, _v, d) = fun_tool_proj_pon_bounded_f(&p_on, &Face(fs.clone()))?;
            let tol_f = fun_tool_maxtol(&fs);
            eon_fs = d < tol_f;
        }

        let eghassd_on = eghassd && esp_on_esd;
        let eghassd3d_on = eghassd3d && esp_on_esd;

        // `cxx:709-726`: `yap00`.
        let mut yap00 = true;
        yap00 = yap00 && !is_closed;
        yap00 = yap00 && eg_bound_for;
        yap00 = yap00 && esp_on_esd;
        yap00 = yap00 && !hsd_for;
        if yap00 {
            let mut l_i_copy = bds.shape_interferences(for_, true).to_vec();
            let mut l_i_gesd = Vec::new();
            let n_gesd = fun_select_gi_interference(&mut l_i_copy, iesd, &mut l_i_gesd);
            if n_gesd != 0 {
                yap00 = false;
            }
        }

        // `cxx:728-742`: `yap0` and `yap0bis`.
        let mut yap0 = true;
        yap0 = yap0 && !is_closed;
        yap0 = yap0 && eg_bound_for;
        yap0 = yap0 && !eghassd;

        let mut yap0bis = true;
        yap0bis = yap0bis && !is_closed;
        yap0bis = yap0bis && eg_bound_for;
        yap0bis = yap0bis && eghassd;
        yap0bis = yap0bis && !eghassd_on;
        yap0bis = yap0bis && eon_fs;
        yap0bis = yap0bis && !com_of_cut;
        if yap0bis {
            yap0 = true;
        }

        // `cxx:745-750`
        let mut hassd3d_on = false;
        if !eg_bound_for && eghassd3d {
            ie3d = fun_find_e_sd(bds, e_on, &eg, for_, &mut oe3d, 3);
            hassd3d_on = ie3d != 0;
        }

        // `cxx:752-789`
        let yapc1 = is_closed && !eg_bound_for && eghassd3d && hassd3d_on;
        let yapc2a = is_closed && !eghassd3d;
        let yapc2b = is_closed && !eg_bound_for && eghassd3d && !hassd3d_on;
        let yapc2c = is_closed && eg_bound_for && eghassd3d;
        let yapc2 = yapc2a || yapc2b || yapc2c;

        let mut yapc3 = !is_closed && eghassd3d && hassd3d_on;
        let mut e3closed_for = false;
        if hassd3d_on {
            let e3d = bds.shape(ie3d, true)?;
            e3closed_for = BRepTool::is_closed_edge_face(&Edge(e3d), &Face(for_.clone()));
        }
        yapc3 = yapc3 && e3closed_for;

        Some(Self {
            for_: for_.clone(),
            esp_on: esp_on.clone(),
            eg,
            fs,
            i_for,
            gi,
            i_fs,
            tb1,
            tb2,
            tfe,
            eg_bound_for,
            o_for,
            o_fs,
            is_closed_ff,
            is_closed_fs,
            is_closed,
            is_rest,
            rank_fs,
            rank_eg,
            rank_for,
            otfe,
            tfe_bef,
            tfe_aft,
            eghassd,
            eghassd3d,
            le3,
            ie3d,
            oe3d,
            ie2d,
            o2d_ori,
            lfcx,
            hsd_for,
            c_for,
            iref_for,
            oref_for,
            for_is_ref,
            ope_fus,
            opec12,
            ope_cut,
            ope_com,
            com_of_cut,
            tb_for,
            tol_a,
            tol_eg,
            tol_for,
            tol_fs,
            factor,
            par_on,
            esp_on_esd,
            iesd,
            eon_fs,
            e_on_so_esd: eon_so_esd,
            eghassd_on,
            eghassd3d_on,
            yap00,
            yap0,
            yap0bis,
            yapc1,
            yapc2a,
            yapc2b,
            yapc2c,
            yapc2,
            yapc3,
            hassd3d_on,
            e3closed_for,
        })
    }
}

impl Builder {
    /// `GFillONPartsWES2(I, EspON)` (`cxx:554-3123`): the branch chain.
    pub(super) fn g_fill_on_parts_wes2(
        &self,
        for_: &TopoShape,
        g: &GTopo,
        i: &Interference,
        esp_on: &TopoShape,
        wes: &mut WireEdgeSet,
    ) {
        let Some(c) = Wes2Ctx::new(self, for_, g, i, esp_on) else {
            return;
        };

        // `cxx:791-877`
        if c.yap00 {
            c.branch_yap00(self, wes);
            return;
        }
        // `cxx:879-1199`
        if c.yap0 {
            c.branch_yap0(self, wes);
            return;
        }
        // `cxx:1201-1297`
        if c.yapc1 {
            c.branch_yapc1(self, wes);
            return;
        }
        // `cxx:1299-1424`
        if c.yapc2 {
            c.branch_yapc2(self, wes);
            return;
        }
        // `cxx:1426-1603`
        if c.yapc3 {
            c.branch_yapc3(self, wes);
            return;
        }
        // `cxx:1605-1790`: the `FCX` block and the `yap1..yap7` predicates.
        let Some(fcx) = c.fcx_block(self) else {
            return;
        };
        // `cxx:1790-2071`
        if fcx.yap1 {
            c.branch_yap1(&fcx, self, wes);
            return;
        }
        // `cxx:2074-2439`
        if fcx.yap2 {
            c.branch_yap2(&fcx, self, wes);
            return;
        }
        // `cxx:2442-2808`
        if fcx.yap6 {
            c.branch_yap6(&fcx, self, wes);
            return;
        }
        // `cxx:2811-2914`
        if fcx.yap3 {
            c.branch_yap3(&fcx, self, wes);
            return;
        }
        // `cxx:2917-3010`
        if fcx.yap5 {
            c.branch_yap5(&fcx, self, wes);
            return;
        }
        // `cxx:3013-3043`
        if fcx.yap4 {
            c.branch_yap4(&fcx, self, wes);
            return;
        }
        // `cxx:3046-3121`
        if fcx.yap7 {
            c.branch_yap7(&fcx, self, wes);
            return;
        }
    }
}
