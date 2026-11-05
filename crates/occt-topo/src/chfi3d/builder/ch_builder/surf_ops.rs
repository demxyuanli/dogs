//! `ChFi3dSurfOps` for [`ChSurfOps`] — the surface half of
//! `ChFi3d_ChBuilder`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_ChBuilder.cxx` (OCCT
//! 8.0.0). The trait carries the eight `ChFi3d_Builder` virtuals that
//! `PerformSetOfSurf` / `StartSol` dispatch to; this file states which of them
//! the chamfer subclass answers with a real body:
//!
//! - `SimulSurf` on two faces (`cxx:717-1222`) — [`super::simul_surf`],
//! - `SimulKPart` (`cxx:655-716`) — [`super::simul_k_part`],
//! - `PerformSurf` on two faces (`cxx:1618-1853`) — [`super::perform_surf`],
//! - `PerformFirstSection` (`cxx:1319-1614`) —
//!   [`super::perform_first_section`],
//!
//! The six edge-contact overloads are **not** overridden here: OCCT's chamfer
//! bodies for all six are `throw Standard_Failure("... Not Implemented")`
//! (`cxx:1223-1310` for the three `SimulSurf`, `:1854-1940` for the three
//! `PerformSurf`), which is exactly the default body of the trait
//! (`ChFi3d_Builder_NotImp.cxx:29-220`).
//!
//! `ChSurfOps` carries no state: every chamfer body reads the section kind
//! (`ChFiDS_ChamfSpine`) and the mode off the spine handle, never off
//! `this->myMode`.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_math::MathVector;

use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi3d::builder::{ChFi3dBuilder, ChFi3dCtx, ChFi3dSurfOps};
use crate::chfi_ds::{ElSpine, SpineHandle, SurfData};
use crate::topopebrep::ds::HDataStructure;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::perform_surf::perform_surf_face_face;
use super::simul_k_part::simul_k_part_impl;
use super::simul_surf::simul_surf_face_face;
use super::ChSurfOps;

impl ChFi3dSurfOps for ChSurfOps {
    /// `PerformSurf(Data, Guide, Spine, Choix, S1, I1, S2, I2, MaxStep,
    /// Fleche, TolGuide, First, Last, Inside, Appro, Forward, RecOnS1,
    /// RecOnS2, Soldep, Intf, Intl)` (`cxx:1618-1853`) — see
    /// [`super::perform_surf::perform_surf_face_face`].
    #[allow(clippy::too_many_arguments)]
    fn perform_surf(
        &mut self,
        base: &mut ChFi3dBuilder,
        data: &mut Vec<Arc<SurfData>>,
        guide: &mut Arc<ElSpine>,
        spine: Option<&SpineHandle>,
        choix: i32,
        s1: &SurfaceAdaptor,
        i1: &Rc<RefCell<TopolTool>>,
        s2: &SurfaceAdaptor,
        i2: &Rc<RefCell<TopolTool>>,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        approx: bool,
        forward: bool,
        rec_on_s1: bool,
        rec_on_s2: bool,
        soldep: &MathVector,
        intf: &mut i32,
        intl: &mut i32,
    ) -> bool {
        perform_surf_face_face(
            base, data, guide, spine, choix, s1, i1, s2, i2, max_step, fleche, tol_guide, first,
            last, inside, approx, forward, rec_on_s1, rec_on_s2, soldep, intf, intl,
        )
    }

    /// `SimulSurf(...)` on two faces (`cxx:717-1222`) — see
    /// [`super::simul_surf::simul_surf_face_face`].
    #[allow(clippy::too_many_arguments)]
    fn simul_surf(
        &mut self,
        data: &mut SurfData,
        guide: &mut Arc<ElSpine>,
        spine: Option<&SpineHandle>,
        choix: i32,
        s1: &SurfaceAdaptor,
        i1: &Rc<RefCell<TopolTool>>,
        s2: &SurfaceAdaptor,
        i2: &Rc<RefCell<TopolTool>>,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        approx: bool,
        forward: bool,
        rec_on_s1: bool,
        rec_on_s2: bool,
        soldep: &MathVector,
        intf: &mut i32,
        intl: &mut i32,
        ctx: &ChFi3dCtx<'_>,
    ) -> bool {
        simul_surf_face_face(
            data, guide, spine, choix, s1, i1, s2, i2, tol_guide, first, last, inside, approx,
            forward, rec_on_s1, rec_on_s2, soldep, intf, intl, ctx,
        )
        // `throw Standard_ConstructionError("SimulSurf : this is not the
        // spine of a chamfer")` (`cxx:739-743`) is the `Err`, which the caller
        // sees as the same "did not simulate".
        .unwrap_or(false)
    }

    /// `PerformFirstSection(Spine, HGuide, Choix, S1, S2, I1, I2, Par,
    /// SolDep, Pos1, Pos2)` (`cxx:1319-1614`) — see
    /// [`super::perform_first_section::perform_first_section`].
    ///
    /// `throw Standard_ConstructionError("PerformSurf : this is not the spine
    /// of a chamfer")` (`cxx:1332-1336`) is the `Err`, which the caller sees
    /// as the same "section not found".
    #[allow(clippy::too_many_arguments)]
    fn perform_first_section(
        &mut self,
        spine: &SpineHandle,
        guide: &Arc<ElSpine>,
        choix: i32,
        s1: &SurfaceAdaptor,
        s2: &SurfaceAdaptor,
        i1: &Rc<RefCell<TopolTool>>,
        i2: &Rc<RefCell<TopolTool>>,
        par: f64,
        sol_dep: &mut MathVector,
        pos1: &mut crate::abs::State,
        pos2: &mut crate::abs::State,
        ctx: &ChFi3dCtx<'_>,
    ) -> bool {
        super::perform_first_section::perform_first_section(
            spine, guide, choix, s1, s2, i1, i2, par, sol_dep, pos1, pos2, ctx,
        )
        .unwrap_or(false)
    }

    /// `SimulKPart(SD)` (`cxx:655-716`) — see
    /// [`super::simul_k_part::simul_k_part_impl`].
    fn simul_k_part(&mut self, sd: &mut SurfData, ds: &HDataStructure) -> Result<(), String> {
        simul_k_part_impl(sd, ds)
    }
}
