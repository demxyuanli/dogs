//! `ChFi3d_Builder::PerformIntersectionAtEnd` — intersection at the end of a
//! fillet with at least two faces.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C1.cxx:1828-3770`
//! (OCCT 8.0.0). The body is split by OCCT line ranges:
//!
//! - [`prologue`] — `cxx:1828-2136`: the face / edge counts (`ChFi3d_nbface`,
//!   `ChFi3d_NbNotDegeneratedEdges`, `ChFi3d_ChercheBordsLibres`), the sewing
//!   edge of the vertex (`ChFi3d_CoutureOnVertex`), the regular-adjacent-edge
//!   analysis (`ChFi3d_edge_common_faces` + `ChFi3d::IsTangentFaces`), the two
//!   `PerformMoreSurfdata` / `PerformMoreThreeCorner` dispatches, the bounded
//!   fillet surface (`ChFi3d_BoundSurf`), the orientation of the end curve and
//!   the `checkShrink` detection.
//!
//! Translation notes:
//! - `PerformMoreSurfdata(Index)` (`cxx:2000-2004`) is UNPORTED
//!   (`ChFi3d_Builder_C1.cxx:3771-4531`), and `PerformMoreThreeCorner(Index, 1)`
//!   (`cxx:2008`, `:2163`) is one of the injected corner entry points, so both
//!   answer an [`EndIntersectionExit`] instead of recursing; the caller
//!   dispatches them like `PerformOneCorner` does for
//!   `PerformIntersectionAtEnd`.
//! - `SeqFil` (`stripe->ChangeSetOfSurfData()->ChangeSequence()`) cannot be
//!   held as a `&mut` into the stripe while the body also reads the DS, so the
//!   port keeps [`EndIntersection::seq`] as a copy and writes it back with
//!   [`EndIntersection::flush`]. The one in-place edit of the prologue,
//!   `SeqFil.Remove(num)` (`cxx:1965`), is applied to the real sequence
//!   immediately because `ChFi3d_IndexOfSurfData` re-reads it right after.

mod prologue;

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::chfi_ds::{FaceInterference, SpineHandle, State, StripeHandle, SurfData};
use crate::shape::{Edge, Face, TopoShape};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::ChFi3dBuilder;

/// The ways `PerformIntersectionAtEnd` hands its corner over to another body.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EndIntersectionExit {
    /// `cxx:2000-2004`: `if (MoreSurfdata(Index)) { PerformMoreSurfdata(Index);
    /// return; }`. The callee is UNPORTED (`cxx:3771-4531`).
    MoreSurfdata,
    /// `cxx:2006-2010` / `cxx:2161-2165`: `PerformMoreThreeCorner(Index, 1);
    /// return;` with the `nbcourb` argument OCCT passes.
    MoreThreeCorner(usize),
    /// The port stops after the prologue: `cxx:2136-3770` is UNPORTED.
    RestUnported,
}

/// The locals of `PerformIntersectionAtEnd` that the prologue produces.
///
/// Fields that only `cxx:2136-3770` reads are still written here, exactly as
/// the still-unported stages expect them.
#[allow(dead_code)]
pub(super) struct EndIntersection {
    /// `Vtx`
    pub vtx: TopoShape,
    /// `stripe`
    pub stripe: StripeHandle,
    /// `spine`
    pub spine: SpineHandle,
    /// `SeqFil` — the stripe's surf-data sequence, kept in sync with the stripe.
    pub seq: Vec<Arc<SurfData>>,
    /// `nface` (`cxx:1858`), possibly incremented by `couture || bordlibre`.
    pub nface: usize,
    /// `nbarete` (`cxx:1861-1870`)
    pub nbarete: i32,
    /// `bordlibre` (`cxx:1862`)
    pub bord_libre: bool,
    /// `edgelibre1` (`cxx:1862`)
    pub edge_libre1: Option<Edge>,
    /// `edgelibre2` (`cxx:1862`)
    pub edge_libre2: Option<Edge>,
    /// `couture` (`cxx:1877-1885`)
    pub couture: bool,
    /// `facecouture`
    pub face_couture: Option<Face>,
    /// `edgecouture`
    pub edge_couture: Option<Edge>,
    /// `nbsurf` (`cxx:1894`)
    pub nbsurf: usize,
    /// `nbedge` (`cxx:1895`)
    pub nbedge: usize,
    /// `num` — 1-based like OCCT.
    pub num: usize,
    /// `num1` — `num + 1` or `num - 1`; outside the sequence when out of range.
    pub num1: i64,
    /// `isfirst` (`cxx:1897`)
    pub is_first: bool,
    /// `state` (`cxx:1898-1910`)
    pub state: State,
    /// `F1` — first face of `EdgeSpine`, later `DStr.Shape(Fd->IndexOfS1())`.
    pub f1: Option<Face>,
    /// `F2` — second face of `EdgeSpine`, later `DStr.Shape(Fd->IndexOfS2())`.
    pub f2: Option<Face>,
    /// `Eadj1`
    pub eadj1: Option<Edge>,
    /// `Eadj2`
    pub eadj2: Option<Edge>,
    /// `reg1`
    pub reg1: bool,
    /// `reg2`
    pub reg2: bool,
    /// `face` — the `FindFace` answer (`cxx:1982-1997`).
    pub face: Option<Face>,
    /// `onecorner` (`cxx:1986-1997`)
    pub one_corner: bool,
    /// `HGs` (`cxx:2011`)
    pub hgs: SurfaceAdaptor,
    /// `Fi1` (`cxx:2012`)
    pub fi1: FaceInterference,
    /// `Fi2` (`cxx:2013`)
    pub fi2: FaceInterference,
    /// `Cc`, `Pc`, `Ps` (`cxx:2017-2018`) — filled by the unported stages.
    pub cc: Option<Arc<dyn Curve>>,
    /// `Pc`
    pub pc: Option<Arc<dyn Curve2d>>,
    /// `Ps`
    pub ps: Option<Arc<dyn Curve2d>>,
    /// `orsurfdata` (`cxx:2021`)
    pub or_surf_data: Orientation,
    /// `Isurf` (`cxx:2022`)
    pub isurf: i32,
    /// `IsurfPrev` (`cxx:2022-2028`)
    pub isurf_prev: i32,
    /// `SDprev` (`cxx:2023-2028`)
    pub sdprev: Option<Arc<SurfData>>,
    /// `tolpt` (`cxx:2031`)
    pub tolpt: f64,
    /// `tolreached` (`cxx:2032`) — filled by the unported stages.
    pub tol_reached: f64,
    /// `orcourbe` (`cxx:2071-2078`)
    pub or_courbe: Orientation,
    /// `Ishape1` (`cxx:2041`)
    pub ishape1: i32,
    /// `trafil1` (`cxx:2042-2058`)
    pub trafil1: Orientation,
    /// `isOnSame1`, `isOnSame2` (`cxx:2082-2089`)
    pub is_on_same1: bool,
    /// `isOnSame2`
    pub is_on_same2: bool,
    /// `checkShrink` (`cxx:2088-2128`)
    pub check_shrink: bool,
    /// `isShrink` (`cxx:2088-2089`)
    pub is_shrink: bool,
    /// `isUShrink` (`cxx:2088-2128`)
    pub is_u_shrink: bool,
    /// `checkShrParam` (`cxx:2090`)
    pub check_shr_param: f64,
    /// `prevSDParam` (`cxx:2090`)
    pub prev_sd_param: f64,
    /// `midP2d` (`cxx:2091`)
    pub mid_p2d: GpPnt2d,
    /// `midIpoint` (`cxx:2092`)
    pub mid_ipoint: i32,
}

impl EndIntersection {
    /// Writes the copied `SeqFil` back into the stripe.
    pub(super) fn flush(&mut self) {
        if let Ok(mut stripe) = self.stripe.write() {
            *stripe.change_set_of_surf_data() = self.seq.clone();
        }
    }
}

impl ChFi3dBuilder {
    /// `PerformIntersectionAtEnd(Index)` (`ChFi3d_Builder_C1.cxx:1828-3770`).
    ///
    /// `Err` carries the throws of the ported prologue; the corner hand-overs
    /// answer an [`EndIntersectionExit`].
    pub fn perform_intersection_at_end(
        &mut self,
        index: usize,
    ) -> Result<EndIntersectionExit, String> {
        let mut st = match prologue::prologue(self, index)? {
            prologue::PrologueOut::Continue(st) => st,
            prologue::PrologueOut::Exit(exit) => return Ok(exit),
        };
        // `cxx:2136-3770` is UNPORTED: "find faces intersecting with the fillet
        // and edges limiting intersections" (`:2130-2596`), the per-face
        // intersection and surface extension (`:2597-3127`), the DS storage
        // (`:3128-3250`), the edge limitations (`:3251-3537`) and the shrink /
        // tolerance tail (`:3538-3769`).
        st.flush();
        Ok(EndIntersectionExit::RestUnported)
    }
}
