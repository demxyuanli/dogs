//! `ChFi3d_Builder::IntersectMoreCorner` — the corner intersection taken when
//! `PerformIntersectionAtEnd` hands the cap over.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C1.cxx:4708-5445`
//! (OCCT 8.0.0). The body is split by OCCT line ranges:
//!
//! - [`setup`] — `cxx:4708-4988`: the first stripe of the vertex (no
//!   `ChFi3d_SelectStripe`), the surplus `SurfData` removals, the cap-side face
//!   resolution (`IFadArc` from the two arcs, `Fv` from `myEFMap`, `Arcprol`
//!   from `myVEMap`, `cherche_face`, `ChFi3d_ExtendSurface`) and the three
//!   `Update` calls that re-locate the corner points,
//! - [`intersection`] — `cxx:5004-5123`: `ChFi3d_BoundSurf` + `ChFi3d_BoundFac`,
//!   `ChFi3d_ComputeCurves`, and the sewing-edge split (`ChFi3d_Couture` +
//!   `Extrema_ExtCC`),
//! - [`ds_write`] — `cxx:5129-5241`: the `IShape` / `Et` orientation, the
//!   `ChFi3d_IndexPointInDS` bookkeeping, and the straight (`!intcouture`) /
//!   split (`intcouture`) DS writes,
//! - [`tail`] — `cxx:5243-5445`: `ChFi3d_EnlargeBox`, the "missing end of
//!   curve" `IZob` between `Fop` and `Fv` (`Arcprolbis`), the final
//!   `ChFi3d_EnlargeBox` calls and the point tolerances.
//!
//! Translation notes:
//! - OCCT keeps `Fd` as a `ChFiDS_SurfData&` into the stripe's sequence; the
//!   port copies it into [`MoreCorner::fd`] and writes it back with
//!   [`MoreCorner::flush`]. The surplus removals of the prologue are applied to
//!   the real sequence at once, because `ChFi3d_IndexOfSurfData` re-reads it.
//! - `Update` / `IntersUpdateOnSame` are the same file statics the port already
//!   uses in `perform_one_corner`; the first two calls of `cxx:4967-4983` are
//!   exactly `IntersUpdateOnSame` with `Fprol = Fv` and `Eprol = Arcprol`, so
//!   they are dispatched through
//!   [`crate::chfi3d::builder_c1::inter_update_on_same`].
//! - `Geom_TrimmedCurve`-valued locals (`curv1` / `curv2`) are
//!   `GeomAdaptor_Curve` over `(first, last)` in the port, for the reason stated
//!   in `perform_one_corner`.
//! - `ComputeCurve2d` is ported (`builder_c1/compute_curve2d.rs`), so the
//!   `intcouture` write takes its `c2d1` / `c2d2` from the real projection
//!   (`ds_write`, `cxx:5188` / `:5192`).
//! - `else { throw Standard_NotImplemented("OneCorner : cap not written"); }`
//!   (`cxx:5125-5128`) is kept as the `Err` of [`ChFi3dBuilder::intersect_more_corner`].
//!
//! Divergences that are not heuristics: a null handle OCCT dereferences (a
//! missing surface, p-curve or arc) becomes an `Err` with the message of the
//! nearest OCCT throw where one exists.

mod ds_write;
mod intersection;
mod setup;
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

/// The locals of `IntersectMoreCorner` that outlive one phase.
pub(super) struct MoreCorner {
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
    /// `Fv`
    pub fv: Option<Face>,
    /// `Fad` — the arc-side face. OCCT assigns it at `cxx:4825` / `cxx:4853` and
    /// never reads it later in this body, so the field documents the local only.
    #[allow(dead_code)]
    pub fad: Option<Face>,
    /// `Fop`
    pub fop: Option<Face>,
    /// `Fopbis` — the fourth face, used by the `OArcprolbis` choice
    /// (`cxx:5326-5340`).
    pub fopbis: Option<Face>,
    /// `Arcpiv`
    pub arcpiv: Option<TopoShape>,
    /// `Arcspine`
    pub arcspine: TopoShape,
    /// `Arcprol`
    pub arcprol: TopoShape,
    /// `Arcprolbis`
    pub arcprolbis: Option<TopoShape>,
    /// `OArcprolv`
    pub o_arcprol_v: Orientation,
    /// `OArcprolop`
    pub o_arcprol_op: Orientation,
    /// `OArcprolbis`
    pub o_arcprolbis: Orientation,
    /// `IFadArc`
    pub i_fad_arc: i32,
    /// `IFopArc`
    pub i_fop_arc: i32,
    /// `Bs`
    pub hbs: Option<SurfaceAdaptor>,
    /// `Bad` — initialized at `cxx:4931`, never read by OCCT afterwards.
    #[allow(dead_code)]
    pub hbad: Option<SurfaceAdaptor>,
    /// `Bop`
    pub hbop: Option<SurfaceAdaptor>,
    /// `HGs` — the fillet surface handed to the `cxx:4984` `Update`. The prologue
    /// already consumed it, so the field only documents the local.
    #[allow(dead_code)]
    pub hgs: Option<Arc<dyn Surface>>,
    /// `Fd`
    pub fd: SurfData,
    /// `p2dbout`
    pub p2dbout: GpPnt2d,
    /// `wop` — `FiadArc.Parameter(isfirst)` (`cxx:4945`), read only by the two
    /// `Update` calls of `cxx:4966` / `cxx:4984`.
    #[allow(dead_code)]
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
}

impl MoreCorner {
    /// Writes the copied `Fd` back into the stripe's sequence. OCCT holds `Fd`
    /// by reference, so its edits survive a throw; the port flushes on both
    /// paths.
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
    /// `IntersectMoreCorner(Index)` (`ChFi3d_Builder_C1.cxx:4708-5445`).
    ///
    /// `Err` carries OCCT's throws of this body, including the
    /// `Standard_NotImplemented("OneCorner : cap not written")` of the cap arm.
    pub fn intersect_more_corner(&mut self, index: usize) -> Result<(), String> {
        let mut mc = setup::setup(self, index)?;
        let result = (|| {
            // `if (inters) { ... } else { throw Standard_NotImplemented("OneCorner :
            // cap not written"); }` (`cxx:5004-5128`).
            if !mc.inters {
                return Err("OneCorner : cap not written".to_string());
            }
            intersection::intersect(self, &mut mc)?;
            ds_write::write(self, &mut mc)?;
            tail::extend(self, &mut mc)
        })();
        mc.flush();
        result
    }
}
