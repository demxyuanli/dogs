//! `ChFi3d_Builder::PerformOneCorner` — one fillet stripe at one corner.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C1.cxx:611-1647`
//! (OCCT 8.0.0). The body is split by OCCT line ranges:
//!
//! - [`setup`] — `cxx:611-931`: stripe selection, the `ChFiDS_OnSame` /
//!   `FindFace` face resolution, `IntersUpdateOnSame` / `Update`, and the
//!   `thePrepareOnSame` early return (`cxx:899-904`),
//! - [`intersection`] — `cxx:933-1098`: `ChFi3d_BoundSurf` + `ChFi3d_BoundFac`,
//!   `ChFi3d_ComputeCurves`, and the sewing-edge split (`ChFi3d_Couture` +
//!   `Extrema_ExtCC`),
//! - [`ds_write`] — `cxx:1100-1393`: the `ChFi3d_IndexPointInDS` bookkeeping,
//!   the straight (`!intcouture`) DS write with the `OCC119` overlap guard, and
//!   the split (`intcouture`) DS write,
//! - [`tail`] — `cxx:1395-1647`: `ChFi3d_EnlargeBox`, the `onsame && inters`
//!   "missing end of curve" block (`VARIANT 1`, `cxx:1403`), the final
//!   `ChFi3d_EnlargeBox` calls and the point tolerances.
//!
//! Translation notes:
//! - `PerformIntersectionAtEnd(Index)` is a corner entry point the port injects
//!   through [`super::ChFi3dCornerOps`]; the two OCCT spots that call it and
//!   return (`cxx:821-825`, `cxx:910-914`) answer
//!   [`OneCornerExit::IntersectionAtEnd`] instead of recursing.
//! - `Fd` (`SeqFil.ChangeValue(num)`) cannot be held as a `&mut SurfData` into
//!   the stripe while the body also reads the DS, so the port copies the entry
//!   into [`OneCorner::fd`] and writes it back with [`OneCorner::flush`]; the
//!   surplus `SurfData` removals (`cxx:637-652`) are applied to the real
//!   sequence immediately, before any copy is taken.
//! - `ChFi3d_TrimmedCurve`-valued locals (`curv1` / `curv2`) are
//!   `GeomAdaptor_Curve` over `(first, last)` in the port: `Geom_TrimmedCurve`
//!   in this port reparameterises to `[0, 1]`, while every DS interference the
//!   body writes carries OCCT basis parameters (`Udeb`, `par2`, `Ufin`).
//! - `ComputeCurve2d` (`cxx:456-483`) is ported in
//!   `builder_c1/compute_curve2d.rs` and wired at the `intcouture` DS write;
//!   the `c2d1` / `c2d2` of `cxx:1338-1342` / `:1569-1575` come from it
//!   (`ds_write::write_split`, `tail`).

mod ds_write;
mod intersection;
pub(super) mod setup;
mod tail;

use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt2d;
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::chfi_ds::{StripeHandle, SurfData};
use crate::shape::{Edge, Face, TopoShape};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::ChFi3dBuilder;

/// The early exits of `PerformOneCorner` that hand the corner over to
/// `PerformIntersectionAtEnd`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OneCornerExit {
    /// `PerformOneCorner` ran to its end (`cxx:1647`).
    Done,
    /// `PerformIntersectionAtEnd(Index)` must run next: either the edge to
    /// prolong was not found (`cxx:821-825`) or `FindFace` failed
    /// (`cxx:910-914`). `ChFi3d_Builder` cannot dispatch to it itself — it is
    /// one of the injected corner entry points — so the caller does.
    IntersectionAtEnd,
}

/// The locals of `PerformOneCorner` that outlive one phase.
pub(super) struct OneCorner {
    /// `Vtx`
    pub vtx: TopoShape,
    /// `stripe`
    pub stripe: StripeHandle,
    /// `SeqFil` — the stripe's surf-data sequence, kept in sync with the stripe.
    pub seq: Vec<Arc<SurfData>>,
    /// `num`, 1-based like OCCT.
    pub num: usize,
    /// `isfirst`
    pub is_first: bool,
    /// `onsame`
    pub onsame: bool,
    /// `Fv`
    pub fv: Option<Face>,
    /// `Fop`
    pub fop: Option<Face>,
    /// `Fad`
    pub fad: Option<Face>,
    /// `Arcpiv`
    pub arcpiv: Option<TopoShape>,
    /// `Arcspine` (`cxx:675-681`) — needed by the `VARIANT 1` block (`cxx:1405`).
    pub arcspine: Option<TopoShape>,
    /// `Arcprol`
    pub arcprol: Option<TopoShape>,
    /// `OArcprolv`
    pub o_arcprol_v: Orientation,
    /// `OArcprolop`
    pub o_arcprol_op: Orientation,
    /// `IFadArc`
    pub i_fad_arc: i32,
    /// `IFopArc`
    pub i_fop_arc: i32,
    /// `Bs`
    pub hbs: Option<SurfaceAdaptor>,
    /// `Bad`
    pub hbad: Option<SurfaceAdaptor>,
    /// `Bop`
    pub hbop: Option<SurfaceAdaptor>,
    /// `HGs` before the `ChFi3d_BoundSurf` reassignment (`cxx:867-868`)
    pub hgs: Option<Arc<dyn Surface>>,
    /// `Fd`
    pub fd: SurfData,
    /// `p2dbout`
    pub p2dbout: GpPnt2d,
    /// `wop`
    pub wop: f64,
    /// `inters`
    pub inters: bool,
    /// `Cc`
    pub cc: Option<Arc<dyn Curve>>,
    /// `Ps` — the p-curve of `Cc` on the fillet surface.
    pub ps: Option<Arc<dyn Curve2d>>,
    /// `Pc` — the p-curve of `Cc` on `Fv`.
    pub pc: Option<Arc<dyn Curve2d>>,
    /// `Udeb`
    pub udeb: f64,
    /// `Ufin`
    pub ufin: f64,
    /// `tolreached`
    pub tolreached: f64,
    /// `couture`
    pub couture: bool,
    /// `edgecouture`
    pub edgecouture: Option<Edge>,
    /// `intcouture`
    pub intcouture: bool,
    /// `par1`
    pub par1: f64,
    /// `par2`
    pub par2: f64,
    /// `indpt`
    pub indpt: i32,
    /// `Icurv1`
    pub icurv1: i32,
    /// `Icurv2`
    pub icurv2: i32,
    /// `curv1`
    pub curv1: Option<Arc<dyn Curve>>,
    /// `curv2`
    pub curv2: Option<Arc<dyn Curve>>,
    /// `box1`
    pub box1: BndBox,
    /// `box2`
    pub box2: BndBox,
    /// Set when the body must stop and answer an early exit.
    pub exit: Option<OneCornerExit>,
}

impl OneCorner {
    /// Writes the copied `Fd` back into the stripe's sequence.
    pub(super) fn flush(&mut self) {
        if let Some(slot) = self.num.checked_sub(1).and_then(|i| self.seq.get_mut(i)) {
            *Arc::make_mut(slot) = self.fd.clone();
        }
        if let Ok(mut stripe) = self.stripe.write() {
            *stripe.change_set_of_surf_data() = self.seq.clone();
        }
    }
}

impl ChFi3dBuilder {
    /// `PerformOneCorner(Index, thePrepareOnSame = false)`
    /// (`ChFi3d_Builder_C1.cxx:611-1647`).
    ///
    /// `Err` carries OCCT's throws of this body; the two
    /// `PerformIntersectionAtEnd` early returns answer
    /// [`OneCornerExit::IntersectionAtEnd`].
    pub fn perform_one_corner(
        &mut self,
        index: usize,
        prepare_on_same: bool,
    ) -> Result<OneCornerExit, String> {
        // `ChFi3d_SelectStripe` answering false returns without touching
        // anything (`cxx:624-627`).
        let Some(mut oc) = setup::setup(self, index, prepare_on_same)? else {
            return Ok(OneCornerExit::Done);
        };
        if let Some(exit) = oc.exit {
            oc.flush();
            return Ok(exit);
        }
        // `else { throw Standard_NotImplemented("OneCorner : bouchon non
        // ecrit"); }` (`cxx:1063-1065`).
        if !oc.inters {
            return Err("OneCorner : bouchon non ecrit".to_string());
        }
        intersection::intersect(self, &mut oc)?;
        ds_write::write(self, &mut oc)?;
        tail::extend(self, &mut oc)?;
        oc.flush();
        Ok(OneCornerExit::Done)
    }
}
