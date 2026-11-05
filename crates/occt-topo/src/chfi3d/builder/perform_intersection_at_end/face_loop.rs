//! `PerformIntersectionAtEnd`'s per-face intersection loop —
//! `ChFi3d_Builder_C1.cxx:2582-3145`.
//!
//! Ported blocks, in OCCT order:
//! - `cxx:2582-2596` — the loop's surviving locals (`Pardeb` / `Parfin`, the
//!   four 2d points, `Hc1` / `Hc2`, `proledge` / `prolface` / `shrink` /
//!   `faceprol` / `indcurve` / `indpoint1` / `indpoint2` / `paredge1` /
//!   `paredge2` / `tolex` / `extend` / `intersOnSameFailed`),
//! - `cxx:2600-2625` — the array clears, `faceprol[1] = facesau`, the
//!   `!isOnSame1 && !isOnSame2` shrink cancel and the `nb = isOnSame1 ? 0 : 1`
//!   start,
//! - `cxx:2627-2657` — the loop head: `E2 = Edge[nb]`, `F = nb ? Face[nb-1] :
//!   F1`, the `faceprol` fill, the null-face throw and `Sfacemoins1 =
//!   BRep_Tool::Surface(F)`,
//! - `cxx:2658-2759` — the `nb == 1` arm: `Hc1 = CurveOnSurface(Edge[0],
//!   Face[0])`, the `isOnSame1` `IntersUpdateOnSame` / point / `SetArc` block,
//!   the `Hc1.IsNull()` extension block with its `Extrema_ExtPS` re-projection,
//!   and the `pfac1` / `paredge1` / `pfil1` / `pfildeb` writes,
//! - `cxx:2760-2766` — the `nb != 1` carry (`pfil1 = pfil2`, `paredge1 =
//!   paredge2`, `pfac1 = pint`),
//! - `cxx:2993-2996` — the `!nb` continue,
//! - `cxx:3054-3126` — `Bs` on face `nb - 1`, `Pardeb` / `Parfin`,
//!   `ChFi3d_Boite` + `ChFi3d_BoundFac`,
//! - `cxx:3131-3145` — `ChFi3d_ComputeCurves`, whose `false` is
//!   `PerformMoreThreeCorner(Index, 1); return;`.
//!
//! `cxx:2767-2952` (the `E2` curve × fillet intersection with the
//! `GeomLib::ExtendSurfByLength` and `GeomInt_IntSS` fallbacks),
//! `cxx:2953-2992` (the `Hc2` arm) and `cxx:2998-3052` (the `isOnSame2`
//! `IntersUpdateOnSame`) are [`super::edge_cut`].
//!
//! Translation notes:
//! - `Bs.Initialize(F, Restriction)` is [`SurfaceAdaptor::of_face`] /
//!   [`SurfaceAdaptor::of_face_restricted`]; `Bs.Initialize(faceprol[i])` is
//!   the same on the face `BRep_Builder::MakeFace(faceprol[i], S, L, Tol)`
//!   builds around the extended surface ([`make_prol_face`]).
//! - `faceprol` stays a `Face` per slot, exactly as in OCCT, even where OCCT
//!   writes it and the loop never reads it back (`cxx:2610`, `:2640`).
//! - `Extrema_ExtPS(P, S, TolU, TolV, Extrema_ExtFlag_MIN)` is read through
//!   `Point(1).Parameter()`: the MIN flag keeps one extremum, which is the
//!   global minimum [`point_surface_extrema`] answers.
//! - `C2dint1` / `C2dint2` (`cxx:2651`) are local to one iteration and are
//!   only written by the `GeomInt_IntSS` fallback, so they live in
//!   [`FaceLoop`] as `c2d_int1` / `c2d_int2`.

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_geom::extrema_surf::point_surface_extrema;
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_boite, ch_fi3d_bound_fac, ch_fi3d_compute_curves, ComputeCurves,
};
use crate::chfi3d::builder_c1::{ch_fi3d_extend_surface, inter_update_on_same};
use crate::chfi_ds::{CommonPoint, FaceInterference};
use crate::shape::Face;
use crate::topopebrep::ds::{InterferenceHandle, Point};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::super::ChFi3dBuilder;
use super::corner_faces::{edge_at, face_at};
use super::{end_points, EndIntersection, EndIntersectionExit, NN};

/// The locals of `cxx:2582-2596` that survive between loop iterations.
pub(super) struct FaceLoop {
    /// `Pardeb(1..4)`
    pub pardeb: [f64; 4],
    /// `Parfin(1..4)`
    pub parfin: [f64; 4],
    /// `pfil1`
    pub pfil1: GpPnt2d,
    /// `pfac1`
    pub pfac1: GpPnt2d,
    /// `pfil2`
    pub pfil2: GpPnt2d,
    /// `pfac2`
    pub pfac2: GpPnt2d,
    /// `pint`
    pub pint: GpPnt2d,
    /// `pfildeb`
    pub pfildeb: GpPnt2d,
    /// `Hc1`
    pub hc1: Option<Arc<dyn Curve2d>>,
    /// `Hc2`
    pub hc2: Option<Arc<dyn Curve2d>>,
    /// `C2dint1` (`cxx:2651`)
    pub c2d_int1: Option<Arc<dyn Curve2d>>,
    /// `C2dint2` (`cxx:2651`)
    pub c2d_int2: Option<Arc<dyn Curve2d>>,
    /// `proledge[nn]`
    pub prol_edge: Vec<i32>,
    /// `prolface[nn + 1]`
    pub prol_face: Vec<i32>,
    /// `shrink[nn]`
    pub shrink: Vec<i32>,
    /// `faceprol[nn]`
    pub face_prol: Vec<Option<Face>>,
    /// `indcurve[nn]` (`cxx:2588`) — the DS index of each face's `Cc`.
    pub ind_curve: Vec<i32>,
    /// `Interfedge[nn]` (`cxx:2590`) — may stay null on the `extend` arm.
    pub interf_edge: Vec<Option<InterferenceHandle>>,
    /// `InterfPC[nn]` (`cxx:2591`)
    pub interf_pc: Vec<Option<InterferenceHandle>>,
    /// `InterfPS[nn]` (`cxx:2591`) — may stay null under `shrink[nb - 1]`.
    pub interf_ps: Vec<Option<InterferenceHandle>>,
    /// `indpoint1`
    pub ind_point1: i32,
    /// `indpoint2`
    pub ind_point2: i32,
    /// `paredge1`
    pub paredge1: f64,
    /// `paredge2` — starts at `0.` (`cxx:2593`).
    pub paredge2: f64,
    /// `tolex` — starts at `1.e-4` (`cxx:2593`).
    pub tolex: f64,
    /// `extend`
    pub extend: bool,
    /// `intersOnSameFailed`
    pub inters_on_same_failed: bool,
    /// `Sfacemoins1` (`cxx:2595`) — persists across iterations like OCCT's.
    pub sface_moins1: Option<Arc<dyn Surface>>,
    /// `Sface` (`cxx:2595`) — written by the `GeomInt_IntSS` fallback only, and
    /// read back by the `extend` arm of `cxx:3477` / `:3520`, so it persists.
    pub sface: Option<Arc<dyn Surface>>,
    /// `cint` (`cxx:2650`) — per-iteration, so cleared at the loop head.
    pub cint: Option<Arc<dyn Curve>>,
    /// `Cc` (`cxx:2017`) — the `ChFi3d_ComputeCurves` answer of the iteration.
    pub cc: Option<Arc<dyn Curve>>,
    /// `Ps` (`cxx:2018`)
    pub ps: Option<Arc<dyn Curve2d>>,
    /// `Pc` (`cxx:2018`)
    pub pc: Option<Arc<dyn Curve2d>>,
    /// `tolreached` (`cxx:2032`)
    pub tol_reached: f64,
    /// `p1` (`cxx:3149`) — `Cc->FirstParameter()`.
    pub p1: f64,
    /// `p2` (`cxx:3150`) — `Cc->LastParameter()`.
    pub p2: f64,
    /// `p2d1` (`cxx:3152`) — `Pc` then `Ps` at `p1` (`cxx:3159-3164`).
    pub p2d1: GpPnt2d,
    /// `p2d2` (`cxx:3152`)
    pub p2d2: GpPnt2d,
    /// `HBs` (`cxx:3058`) — the bounded face adaptor.
    pub hbs: Option<SurfaceAdaptor>,
    /// `orface` (`cxx:3237`)
    pub orface: Orientation,
    /// `ori` (`cxx:3237`)
    pub ori: Orientation,
    /// `orien` (`cxx:3240-3245`)
    pub orien: Orientation,
}

impl FaceLoop {
    /// `cxx:2582-2596`.
    pub(super) fn new() -> Self {
        Self {
            pardeb: [0.0; 4],
            parfin: [0.0; 4],
            pfil1: GpPnt2d::zero(),
            pfac1: GpPnt2d::zero(),
            pfil2: GpPnt2d::zero(),
            pfac2: GpPnt2d::zero(),
            pint: GpPnt2d::zero(),
            pfildeb: GpPnt2d::zero(),
            hc1: None,
            hc2: None,
            c2d_int1: None,
            c2d_int2: None,
            prol_edge: vec![0; NN],
            prol_face: vec![0; NN + 1],
            shrink: vec![0; NN],
            face_prol: vec![None; NN],
            ind_curve: vec![0; NN],
            interf_edge: vec![None; NN],
            interf_pc: vec![None; NN],
            interf_ps: vec![None; NN],
            ind_point1: 0,
            ind_point2: 0,
            paredge1: 0.0,
            paredge2: 0.0,
            tolex: 1.0e-4,
            extend: false,
            inters_on_same_failed: false,
            sface_moins1: None,
            sface: None,
            cint: None,
            cc: None,
            ps: None,
            pc: None,
            tol_reached: 0.0,
            p1: 0.0,
            p2: 0.0,
            p2d1: GpPnt2d::zero(),
            p2d2: GpPnt2d::zero(),
            hbs: None,
            orface: Orientation::Forward,
            ori: Orientation::Forward,
            orien: Orientation::Forward,
        }
    }

    /// `Bs.Initialize(prolface[i] ? faceprol[i] : Face[i], Restriction)`.
    ///
    /// Every OCCT call site pairs the `prolface` slot with the same `Face[]`
    /// index (`cxx:2664-2671`, `:3001-3008`, `:3054-3061`), so `i` is both.
    pub(super) fn face_adaptor(
        &self,
        st: &EndIntersection,
        i: usize,
        restriction: bool,
    ) -> Result<SurfaceAdaptor, String> {
        let face = if self.prol_face[i] != 0 {
            self.face_prol[i].clone()
        } else {
            Some(face_at(st, i, 0)?)
        };
        let face = face.ok_or_else(|| {
            "PerformIntersectionAtEnd: Bs.Initialize has a null face".to_string()
        })?;
        Ok(if restriction {
            SurfaceAdaptor::of_face_restricted(&face)
        } else {
            SurfaceAdaptor::of_face(&face)
        })
    }
}

/// `BRep_Builder::MakeFace(faceprol[i], S, L, Tol)`.
pub(super) fn make_prol_face(surface: Arc<dyn Surface>, tolerance: f64) -> Face {
    TopoBuilder::new().make_face_with_tolerance(surface, &[], tolerance)
}

/// The `cxx:2582-3145` loop. `Ok(Some(exit))` is one of the two
/// `PerformMoreThreeCorner(Index, 1); return;` hand-overs inside it.
pub(super) fn intersect_faces(
    b: &mut ChFi3dBuilder,
    st: &mut EndIntersection,
) -> Result<Option<EndIntersectionExit>, String> {
    let mut lp = FaceLoop::new();
    let nbface = st.nbface;
    if nbface < 0 {
        return Err("PerformIntersectionAtEnd: negative nbface".to_string());
    }
    // `for (nb = 1; nb <= nbface; nb++) { prolface[nb - 1] = 0; proledge[nb - 1]
    //  = 0; shrink[nb - 1] = 0; }` (`cxx:2600-2605`).
    for i in 0..nbface as usize {
        lp.prol_face[i] = 0;
        lp.prol_edge[i] = 0;
        lp.shrink[i] = 0;
    }
    // `proledge[nbface] = 0; prolface[nn] = 0;` (`cxx:2606-2607`).
    lp.prol_edge[nbface as usize] = 0;
    lp.prol_face[NN] = 0;
    // `if (oneintersection1 || oneintersection2) faceprol[1] = facesau;`
    // (`cxx:2608-2611`).
    if st.one_intersection1 || st.one_intersection2 {
        if let Some(saved) = st.face_sau.clone() {
            lp.face_prol[1] = Some(saved);
        }
    }
    // `if (!isOnSame1 && !isOnSame2) checkShrink = false;` (`cxx:2612-2615`).
    if !st.is_on_same1 && !st.is_on_same2 {
        st.check_shrink = false;
    }
    // `if (isOnSame1) nb = 0; else nb = 1;` (`cxx:2617-2624`).
    let mut nb: i32 = if st.is_on_same1 { 0 } else { 1 };

    // `for (; nb <= nbface; nb++)` (`cxx:2627`).
    while nb <= nbface {
        let (mut cv1, mut cv2) = end_points(st)?;
        let mut fi1 = st.fi1.clone();
        let mut fi2 = st.fi2.clone();
        // `extend = false; E2 = Edge[nb];` (`cxx:2629-2630`).
        lp.extend = false;
        lp.c2d_int1 = None;
        lp.c2d_int2 = None;
        lp.cint = None;
        let e2 = edge_at(st, nb as usize, 2630)?;
        // `F = nb ? Face[nb - 1] : F1; if (nb && !prolface[nb - 1]) faceprol[nb
        //  - 1] = F;` (`cxx:2631-2642`).
        let f = if nb == 0 {
            st.f1
                .clone()
                .ok_or_else(|| "PerformIntersectionAtEnd: F1 is null".to_string())?
        } else {
            let f = face_at(st, (nb - 1) as usize, 2637)?;
            if lp.prol_face[(nb - 1) as usize] == 0 {
                lp.face_prol[(nb - 1) as usize] = Some(f.clone());
            }
            f
        };
        // `if (F.IsNull()) throw Standard_NullObject(...)` (`cxx:2644-2647`).
        // `face_at` already answered a null face as an error.
        // `Sfacemoins1 = BRep_Tool::Surface(F);` (`cxx:2649`).
        let mut sface_moins1 = BRepTool::face_surface(&f)
            .ok_or_else(|| "IntersectionAtEnd: Trying to intersect with NULL face".to_string())?;

        if nb == 1 {
            first_side(b, st, &mut lp, &mut cv1, &mut fi1, &f, &mut sface_moins1)?;
        } else {
            // `pfil1 = pfil2; paredge1 = paredge2; pfac1 = pint;`
            // (`cxx:2762-2764`).
            lp.pfil1 = lp.pfil2;
            lp.paredge1 = lp.paredge2;
            lp.pfac1 = lp.pint;
        }

        // `if (nb != nbface || isOnSame2)` (`cxx:2767`), else `cxx:2953-2992`.
        let exit = if nb != nbface || st.is_on_same2 {
            super::edge_cut::cut(b, st, &mut lp, &e2, &f, &mut sface_moins1, nb)?
        } else {
            super::edge_cut::last_edge_pcurve(st, &mut lp, &mut cv2, &mut fi2, &e2, &mut sface_moins1)?;
            None
        };
        // `if (!nb) continue;` (`cxx:2993-2996`). The `exit` of `cut` returns
        // through the same `PerformMoreThreeCorner(Index, 1); return;` the
        // `!nb` path never reaches.
        if nb == 0 || exit.is_some() {
            store_all(st, &cv1, &cv2, &fi1, &fi2);
            if let Some(exit) = exit {
                return Ok(Some(exit));
            }
            nb += 1;
            continue;
        }
        // `if (nb == nbface && isOnSame2)` (`cxx:2998-3052`).
        if nb == nbface && st.is_on_same2 {
            super::edge_cut::same_end_update(b, st, &mut lp, &mut cv2, &mut fi2, &e2, nb)?;
        }
        // `Sfacemoins1` (`cxx:2595`) is the loop-persistent handle the storage
        // of `cxx:3429-3450` / `:3472-3485` reads back.
        lp.sface_moins1 = Some(Arc::clone(&sface_moins1));
        // `Bs.Initialize(prolface[nb - 1] ? faceprol[nb - 1] : Face[nb - 1]);`
        // (`cxx:3054-3061`).
        let bs = lp.face_adaptor(st, (nb - 1) as usize, true)?;

        // `Pardeb(1..4) = ...; Parfin(1..4) = ...;` (`cxx:3114-3121`).
        lp.pardeb = [lp.pfil1.x(), lp.pfil1.y(), lp.pfac1.x(), lp.pfac1.y()];
        lp.parfin = [lp.pfil2.x(), lp.pfil2.y(), lp.pfac2.x(), lp.pfac2.y()];
        // `ChFi3d_Boite(pfac1, pfac2, uu1, uu2, vv1, vv2);
        //  ChFi3d_BoundFac(Bs, uu1, uu2, vv1, vv2);` (`cxx:3123-3125`).
        let boite = ch_fi3d_boite(&lp.pfac1, &lp.pfac2);
        let hbs = ch_fi3d_bound_fac(&bs, boite.u_min, boite.u_max, boite.v_min, boite.v_max);
        // `if (!ChFi3d_ComputeCurves(HGs, HBs, Pardeb, Parfin, Cc, Ps, Pc,
        //  tolapp3d, tol2d, tolreached, nbface == 1)) {
        //  PerformMoreThreeCorner(Index, 1); return; }` (`cxx:3131-3145`).
        let mut out = ComputeCurves::default();
        let computed = ch_fi3d_compute_curves(
            &st.hgs,
            &hbs,
            lp.pardeb,
            lp.parfin,
            &mut out,
            b.tolapp3d,
            b.tol2d,
            nbface == 1,
        );
        store_all(st, &cv1, &cv2, &fi1, &fi2);
        if !computed {
            return Ok(Some(EndIntersectionExit::MoreThreeCorner(1)));
        }
        // `HBs` is the `ChFi3d_BoundFac` adaptor the storage samples
        // (`cxx:3156-3158`).
        lp.hbs = Some(hbs);
        lp.cc = out.c3d.clone();
        lp.ps = out.pc1.clone();
        lp.pc = out.pc2.clone();
        lp.tol_reached = out.tolreached;
        // `cxx:3146-3264`: the point / curve storage and the two end-edge
        // interferences.
        super::ds_store::store(b, st, &mut lp, &cv1, &cv2, nb)?;
        // `cxx:3262-3532`: the limitation of `Edge[nb]`.
        super::edge_limits::limit(b, st, &mut lp, &cv1, &cv2, &e2, nb)?;
        // `cxx:3534-3572`: the shrunk part of the end.
        super::shrink::shrink(b, st, &mut lp, nb)?;
        // `cxx:3576-3585`: `InterfPC` / `InterfPS` and the `indpoint1 =
        //  indpoint2` carry.
        super::ds_store::store_face(b, st, &mut lp, nb)?;
        nb += 1;
    }
    // `cxx:3587-3610`: the `CV1` / `CV2` resets and the appends onto the faces
    // and edges.
    super::ds_store::finish(b, st, &mut lp)?;
    // `cxx:3611-3766`: `stripe->InDS(isfirst)` or the shrink tail.
    super::shrink::tail(b, st, &mut lp)?;
    Ok(None)
}

/// `cxx:2658-2759`: the `nb == 1` arm.
fn first_side(
    b: &mut ChFi3dBuilder,
    st: &EndIntersection,
    lp: &mut FaceLoop,
    cv1: &mut CommonPoint,
    fi1: &mut FaceInterference,
    f: &Face,
    sface_moins1: &mut Arc<dyn Surface>,
) -> Result<(), String> {
    // `Hc1 = BRep_Tool::CurveOnSurface(Edge[0], Face[0], Ubid, Ubid);`
    // (`cxx:2660`).
    let e0 = edge_at(st, 0, 2660)?;
    let f0 = face_at(st, 0, 2660)?;
    lp.hc1 = curve_on_surface_range(&e0, &f0).map(|(c, _, _)| c);
    if st.is_on_same1 {
        // `Bs.Initialize(prolface[0] ? faceprol[0] : Face[0], false);`
        // (`cxx:2664-2671`).
        let bs = lp.face_adaptor(st, 0, false)?;
        let hgs = st
            .hgs
            .surface()
            .cloned()
            .ok_or_else(|| "PerformIntersectionAtEnd: HGs has no surface".to_string())?;
        let hbs = bs
            .surface()
            .cloned()
            .ok_or_else(|| "PerformIntersectionAtEnd: Bs has no surface".to_string())?;
        // `const handle<Geom_Curve>& c3df = DStr.Curve(Fi1.LineIndex()).Curve();`
        // (`cxx:2672`).
        let c3d_fi = {
            let ds = b.my_ds.change_ds();
            let line = fi1.line_index();
            if line <= 0 {
                return Err(
                    "PerformIntersectionAtEnd: Fi1.LineIndex() is not in the DS".to_string(),
                );
            }
            ds.curve(line).curve().cloned().ok_or_else(|| {
                "PerformIntersectionAtEnd: DStr.Curve(Fi1.LineIndex()) is null".to_string()
            })?
        };
        // `double Ufi = Fi2.Parameter(isfirst); ... IntersUpdateOnSame(HGs, HBs,
        //  c3df, F1, Face[0], Edge[0], Vtx, isfirst, 10 * tolapp3d, Fi, CV1,
        //  pfac1, Ufi)` (`cxx:2673-2690`).
        let f1 = st
            .f1
            .clone()
            .ok_or_else(|| "PerformIntersectionAtEnd: F1 is null".to_string())?;
        let mut ufi = st.fi2.parameter(st.is_first);
        let ok = inter_update_on_same(
            &hgs,
            &hbs,
            &c3d_fi,
            &f1,
            &f0,
            &e0,
            &st.vtx,
            st.is_first,
            10.0 * b.tolapp3d,
            fi1,
            cv1,
            &mut lp.pfac1,
            &mut ufi,
        )?;
        if !ok {
            return Err("IntersectionAtEnd: pb intersection Face - Fi".to_string());
        }
        // `if (intersOnSameFailed) { proj.Init(pfac1, C2dint2.IsNull() ? Hc1 :
        //  C2dint2); paredge2 = proj.LowerDistanceParameter(); }`
        // (`cxx:2692-2705`). `C2dint2` is null on this arm.
        if lp.inters_on_same_failed {
            let hc1 = lp
                .hc1
                .clone()
                .ok_or_else(|| "IntersectionAtEnd: Hc1 is null".to_string())?;
            lp.paredge2 = project_2d(&hc1, &lp.pfac1)?;
        }
        // `TopOpeBRepDS_Point tpoint(CV1.Point(), tolapp3d); indpoint1 =
        //  DStr.AddPoint(tpoint); stripe->SetIndexPoint(indpoint1, isfirst,
        //  1);` (`cxx:2707-2709`).
        let p = *cv1.point();
        lp.ind_point1 = add_point(b, p, b.tolapp3d);
        if let Ok(mut stripe) = st.stripe.write() {
            stripe.set_index_point(lp.ind_point1, st.is_first, 1);
        }
        // `TopExp::Vertices(Edge[0], vert1, vert2); arcOri = Vtx.IsSame(vert1) ?
        //  FORWARD : REVERSED; CV1.SetArc(tolapp3d, Edge[0], paredge2,
        //  arcOri);` (`cxx:2711-2714`).
        let ori = super::edge_cut::arc_orientation(&e0, &st.vtx);
        cv1.set_arc(b.tolapp3d, &e0.0, lp.paredge2, ori);
    } else if lp.hc1.is_none() {
        // `ChFi3d_ExtendSurface(Sfacemoins1, prolface[0]);` (`cxx:2722`).
        *sface_moins1 = ch_fi3d_extend_surface(sface_moins1, &mut lp.prol_face[0]);
        if lp.prol_face[0] != 0 {
            // `extend = true; BRE.MakeFace(faceprol[0], Sfacemoins1, F.Location(),
            //  tol);` (`cxx:2725-2728`).
            lp.extend = true;
            let tol = BRepTool::face_tolerance(f);
            lp.face_prol[0] = Some(make_prol_face(Arc::clone(sface_moins1), tol));
            if !st.is_on_same1 {
                // `Extrema_ExtPS ext(CV1.Point(), Asurf(Sfacemoins1), tol, tol,
                //  Extrema_ExtFlag_MIN); if (ext.IsDone()) pfac1 =
                //  ext.Point(1).Parameter(uc1, vc1);` (`cxx:2731-2741`).
                let p = *cv1.point();
                if let Some(pair) = point_surface_extrema(sface_moins1.as_ref(), &p) {
                    lp.pfac1 = GpPnt2d::new(pair.u1, pair.v1.unwrap_or(0.0));
                }
            }
        }
    } else {
        // `pfac1 = Hc1->Value(CV1.ParameterOnArc());` (`cxx:2746`).
        let hc1 = lp.hc1.as_ref().expect("checked above");
        let par = cv1
            .parameter_on_arc()
            .ok_or_else(|| "IntersectionAtEnd: CV1 is not on an arc".to_string())?;
        lp.pfac1 = hc1.value(par);
    }
    // `paredge1 = CV1.ParameterOnArc();` (`cxx:2749`).
    lp.paredge1 = cv1
        .parameter_on_arc()
        .ok_or_else(|| "IntersectionAtEnd: CV1 is not on an arc".to_string())?;
    // `pfil1 = Fi1.PCurveOnSurf()->Value(Fi1.Parameter(...)); pfildeb = pfil1;`
    // (`cxx:2750-2758`).
    lp.pfil1 = fi_param_on_surf(fi1, st.is_first)?;
    lp.pfildeb = lp.pfil1;
    Ok(())
}

/// The `isfirst` / `!isfirst` parameter pick of `Fi.PCurveOnSurf()->Value(...)`
/// (`cxx:2750-2756`, `:2984-2990`).
pub(super) fn fi_param_on_surf(
    fi: &FaceInterference,
    is_first: bool,
) -> Result<GpPnt2d, String> {
    let pc = fi
        .p_curve_on_surf()
        .ok_or_else(|| "PerformIntersectionAtEnd: Fi has no p-curve on surf".to_string())?;
    let par = if fi.line_index() != 0 {
        fi.parameter(is_first)
    } else {
        fi.parameter(!is_first)
    };
    Ok(pc.value(par))
}

/// `Geom2dAPI_ProjectPointOnCurve::Init(P, C)` + `LowerDistanceParameter()`.
pub(super) fn project_2d(c: &Arc<dyn Curve2d>, p: &GpPnt2d) -> Result<f64, String> {
    occt_geom2d::geom2d_api::project_point_on_curve_2d(c.as_ref(), p, occt_core::precision::PCONFUSION)
        .map(|proj| proj.parameter)
        .ok_or_else(|| "IntersectionAtEnd: Geom2dAPI_ProjectPointOnCurve found nothing".to_string())
}

/// `DStr.AddPoint(TopOpeBRepDS_Point(P, tol))`.
pub(super) fn add_point(b: &mut ChFi3dBuilder, p: GpPnt, tol: f64) -> i32 {
    b.my_ds
        .change_ds()
        .add_point(&Point::with_tolerance(p, tol))
}

/// Writes the four `Change*` references back into `SeqFil(num)`.
fn store_all(
    st: &mut EndIntersection,
    cv1: &CommonPoint,
    cv2: &CommonPoint,
    fi1: &FaceInterference,
    fi2: &FaceInterference,
) {
    let Some(slot) = st.num.checked_sub(1).and_then(|i| st.seq.get_mut(i)) else {
        return;
    };
    let fd = Arc::make_mut(slot);
    *fd.vertex_mut(st.is_first, 1) = cv1.clone();
    *fd.vertex_mut(st.is_first, 2) = cv2.clone();
    *fd.interference_on_s1_mut() = fi1.clone();
    *fd.interference_on_s2_mut() = fi2.clone();
}
