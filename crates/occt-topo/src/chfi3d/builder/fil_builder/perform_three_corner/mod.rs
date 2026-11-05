//! `ChFi3d_FilBuilder::PerformThreeCorner` - the corner shared by three fillet
//! stripes of a vertex.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder_C3.cxx:241-1331`
//! (OCCT 8.0.0). The file statics of the same `.cxx` are ported next to the
//! body: `SearchPivot` (`cxx:78-105`) already lives in
//! [`crate::chfi3d::builder_helpers::search_pivot`], `SearchFD` (`cxx:109-211`)
//! in [`crate::chfi3d::builder_0_helpers::ch_fi3d_search_fd`], and
//! `ToricCorner` (`cxx:219-234`) is [`toric_corner`] here.
//!
//! The body is split by OCCT line ranges:
//!
//! - this file - `cxx:241-739`: the three stripes and their `ChFi3d_IndexOfSurfData`
//!   indices, the three `SearchFD` calls, the pivot analysis of `cxx:345-391`,
//!   the four `ChFi3d_ExtrSpineCarac` probes and the torus / sphere / `c1pointu`
//!   criteria of `cxx:412-537`, then the surface frame, the `Choix` of
//!   `cxx:610-646`, the two face p-curves and the direct `ComputeCorner` calls,
//! - [`fill`] - `cxx:739-1015`: the guideline and `ComputeData` path, then the
//!   `GeomFill_ConstrainedFilling` path,
//! - [`ds_update`] - `cxx:1016-1331`: the DS update of the four stripes and the
//!   truncation of the three `SurfData` sequences.
//!
//! Translation notes:
//! - OCCT's `done` is [`ChFi3dBuilder::done`](super::super::ChFi3dBuilder),
//! - `CD[i]` are `ChFiDS_Stripe` handles, so `CD[0] == CD[1]` of `cxx:280-295`
//!   is [`Arc::ptr_eq`],
//! - the `SurfData` OCCT holds by reference (`fddeb` / `fdfin` / `fdpiv`) are
//!   written in place through [`with_sd`], which locks the stripe and takes
//!   `Arc::make_mut` on the sequence slot: the three slots alias exactly when
//!   the C++ handles do,
//! - the four `Bnd_Box` of `cxx:1046-1054` are one array with an alias map,
//!   because `pbl1 = pbf1` (`cxx:1048-1051`) makes two of the names the same
//!   box when `c1pointu`,
//! - a null handle OCCT dereferences (a missing stripe, an out-of-range
//!   `SurfData`, a p-curve that is not there, a face index that is not a face)
//!   is the `Err` of this port.

mod ds_update;
mod fill;

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_core::precision;

use crate::abs::Orientation;
use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_extr_spine_carac, ch_fi3d_index_of_surf_data, ch_fi3d_search_fd,
};
use crate::chfi3d::builder_helpers::search_pivot;
use crate::chfi_ds::{Stripe, StripeHandle, SurfData};
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::shape::{Face, TopoShape};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::super::corner_sd::{sd, sd_len, with_sd};
use super::super::ChFi3dBuilder;

/// `ToricCorner(F, rd, rf, v)` (`cxx:219-234`): the corner is a torus, or a
/// sphere limited by isos, when the two radii agree, the pivot face is a plane
/// and `v` lies in that plane's normal direction.
fn toric_corner(face: &Face, rd: f64, rf: f64, v: &GpVec) -> bool {
    if (rd - rf).abs() > precision::CONFUSION {
        return false;
    }
    let bs = SurfaceAdaptor::of_face(face);
    if bs.get_type() != GeomAbsSurfaceType::Plane {
        return false;
    }
    let Some(pl) = bs.plane() else {
        return false;
    };
    // `double scal1 = std::abs(bs.Plane().Position().XDirection().Dot(v));`
    let pos = pl.position();
    let xd = pos.x_direction();
    let yd = pos.y_direction();
    let scal1 = (xd.x() * v.x() + xd.y() * v.y() + xd.z() * v.z()).abs();
    let scal2 = (yd.x() * v.x() + yd.y() * v.y() + yd.z() * v.z()).abs();
    scal1 <= precision::CONFUSION && scal2 <= precision::CONFUSION
}

/// The state the three OCCT ranges share.
pub(super) struct ThreeCorner {
    /// `CD[0..3]` (`cxx:257`).
    pub stripes: [StripeHandle; 3],
    /// `Index[0..3]` (`cxx:252`).
    pub index: [i32; 3],
    /// `sens[0..3]` (`cxx:268`).
    pub sens: [i32; 3],
    /// `i[0..3][0..3]` (`cxx:266`).
    pub i: [[i32; 3]; 3],
    /// `jf[0..3][0..3]` (`cxx:260`).
    pub jf: [[i32; 3]; 3],
    /// `p[0..3][0..3]` (`cxx:268`).
    pub p: [[f64; 3]; 3],
    /// `oksea[0..3]` (`cxx:261`).
    pub oksea: [bool; 3],
    /// `sameside[0..3]` (`cxx:261`).
    pub sameside: [bool; 3],
    /// `face[0..3]` (`cxx:259`).
    pub face: [Option<Face>; 3],
    /// `pivot`, `deb`, `fin` (`cxx:252`).
    pub pivot: usize,
    /// `deb` (`cxx:252`).
    pub deb: usize,
    /// `fin` (`cxx:252`).
    pub fin: usize,
    /// `pivdif` (`cxx:254`).
    pub pivdif: bool,
    /// `c1pointu` (`cxx:255`).
    pub c1pointu: bool,
    /// `c1toric` (`cxx:256`).
    pub c1toric: bool,
    /// `c1spheric` (`cxx:257`).
    pub c1spheric: bool,
    /// `filling` (`cxx:270`).
    pub filling: bool,
    /// `Vtx` (`cxx:250`).
    pub vtx: TopoShape,
}

/// The locals the [`fill`] range needs from the head.
pub(super) struct CornerFrame {
    /// `OFac` (`cxx:551`): the pivot face orientation.
    pub o_fac: Orientation,
    /// `Fac` (`cxx:552`): `BRepAdaptor_Surface(face[pivot])`.
    pub fac: SurfaceAdaptor,
    /// `Surf` (`cxx:595`): `GeomAdaptor_Surface(gasurf)`.
    pub surf: SurfaceAdaptor,
    /// `IFac` (`cxx:600`).
    pub ifac: std::rc::Rc<std::cell::RefCell<crate::brep_top_adaptor_topol_tool::TopolTool>>,
    /// `ISurf` (`cxx:603`).
    pub isurf: std::rc::Rc<std::cell::RefCell<crate::brep_top_adaptor_topol_tool::TopolTool>>,
    /// `oo1` (`cxx:612`).
    pub oo1: Orientation,
    /// `oo2` (`cxx:613`).
    pub oo2: Orientation,
    /// `o1` (`cxx:611`).
    pub o1: Orientation,
    /// `o2` (`cxx:612`).
    pub o2: Orientation,
    /// `choix` (`cxx:610`).
    pub choix: i32,
    /// `Rdeb` (`cxx:412`).
    pub rdeb: f64,
    /// `Rfin` (`cxx:412`).
    pub rfin: f64,
    /// `Rdp` (`cxx:412`).
    pub rdp: f64,
    /// `Rfp` (`cxx:412`).
    pub rfp: f64,
    /// `Pdeb` (`cxx:413`).
    pub pdeb: GpPnt,
    /// `Vdeb` (`cxx:414`).
    pub vdeb: GpVec,
    /// `Pfin` (`cxx:413`).
    pub pfin: GpPnt,
    /// `Vfin` (`cxx:414`).
    pub vfin: GpVec,
    /// `pfac1` (`cxx:647`).
    pub pfac1: GpPnt2d,
    /// `vfac1` (`cxx:648`).
    pub vfac1: GpVec2d,
    /// `pfac2` (`cxx:647`).
    pub pfac2: GpPnt2d,
    /// `vfac2` (`cxx:648`).
    pub vfac2: GpVec2d,
    /// `psurf1` (`cxx:647`).
    pub psurf1: GpPnt2d,
    /// `psurf2` (`cxx:647`).
    pub psurf2: GpPnt2d,
}

/// `ChFi3d_ExtrSpineCarac(DStr, CD[k], i, p, jf, sens[k], ...)` for one stripe.
fn extr(
    base: &ChFi3dBuilder,
    c: &ThreeCorner,
    k: usize,
    j: usize,
    jf: i32,
) -> Result<(GpPnt, GpVec, f64), String> {
    let guard = c.stripes[k]
        .read()
        .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
    ch_fi3d_extr_spine_carac(
        &base.my_ds,
        &guard,
        c.i[k][j],
        c.p[k][j],
        jf,
        c.sens[k],
    )
}

/// `SearchFD(DStr, cd1, cd2, sens1, sens2, ..., ind1, ind2, ...)` (`cxx:296-339`).
#[allow(clippy::too_many_arguments)]
fn search_fd(
    base: &ChFi3dBuilder,
    c: &ThreeCorner,
    a: usize,
    b: usize,
) -> Result<crate::chfi3d::builder_0_helpers::SearchFdAnswer, String> {
    let ga = c.stripes[a]
        .read()
        .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
    let gb = c.stripes[b]
        .read()
        .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
    ch_fi3d_search_fd(
        &base.my_ds,
        &ga,
        &gb,
        c.sens[a],
        c.sens[b],
        c.index[a],
        c.index[b],
    )
}

/// `PerformThreeCorner(Index)` (`cxx:241-1331`).
///
/// `shape` is the `myShape` of `ChFi3d_FilBuilder` (`hxx:377`), read by the two
/// blend functions of [`fill`].
pub(super) fn perform_three_corner(
    base: &mut ChFi3dBuilder,
    shape: crate::blend::blend_func::BlendFuncSectionShape,
    index: usize,
) -> Result<(), String> {
    // `done = false; const TopoDS_Vertex& Vtx = myVDataMap.FindKey(Jndex);`
    // (`cxx:250-251`).
    base.done = false;
    let vtx = base
        .my_v_data_map
        .find_key(index)
        .cloned()
        .ok_or_else(|| format!("PerformThreeCorner: myVDataMap has no key {index}"))?;

    // `for (It.Initialize(myVDataMap(Jndex)), ii = 0; It.More() && ii < 3;
    // It.Next(), ii++) { Index[ii] = ChFi3d_IndexOfSurfData(Vtx, It.Value(),
    // sens[ii]); CD[ii] = It.Value(); }` (`cxx:272-277`).
    let stripes: Vec<StripeHandle> = base
        .my_v_data_map
        .get_from_index(index)
        .map(|s| s.to_vec())
        .ok_or_else(|| format!("PerformThreeCorner: myVDataMap({index}) is out of range"))?;
    if stripes.len() < 3 {
        return Err(format!(
            "PerformThreeCorner: the vertex carries {} stripe(s), three are required",
            stripes.len()
        ));
    }
    let mut c = ThreeCorner {
        stripes: [stripes[0].clone(), stripes[1].clone(), stripes[2].clone()],
        index: [0; 3],
        sens: [0; 3],
        i: [[0; 3]; 3],
        jf: [[0; 3]; 3],
        p: [[0.0; 3]; 3],
        oksea: [false; 3],
        sameside: [false; 3],
        face: [None, None, None],
        pivot: 0,
        deb: 0,
        fin: 0,
        pivdif: true,
        c1pointu: false,
        c1toric: false,
        c1spheric: false,
        filling: false,
        vtx: vtx.clone(),
    };
    for ii in 0..3 {
        let mut sens = 0;
        c.index[ii] = ch_fi3d_index_of_surf_data(&vtx, &c.stripes[ii], &mut sens)?;
        c.sens[ii] = sens;
    }

    // `if (CD[0] == CD[1]) { Index[1] = ...Length(); sens[1] = -1; } else if
    // (CD[1] == CD[2]) { ... } else if (CD[0] == CD[2]) { ... }`
    // (`cxx:280-295`).
    if Arc::ptr_eq(&c.stripes[0], &c.stripes[1]) {
        c.index[1] = sd_len(&c.stripes[1])?;
        c.sens[1] = -1;
    } else if Arc::ptr_eq(&c.stripes[1], &c.stripes[2]) {
        c.index[2] = sd_len(&c.stripes[2])?;
        c.sens[2] = -1;
    } else if Arc::ptr_eq(&c.stripes[0], &c.stripes[2]) {
        c.index[2] = sd_len(&c.stripes[2])?;
        c.sens[2] = -1;
    }

    // `oksea[2] = SearchFD(DStr, CD[0], CD[1], ...)` (`cxx:296-309`).
    let a = search_fd(base, &c, 0, 1)?;
    c.oksea[2] = a.found;
    c.i[0][1] = a.i1;
    c.i[1][0] = a.i2;
    c.p[0][1] = a.p1;
    c.p[1][0] = a.p2;
    c.face[2] = a.face.clone();
    c.sameside[2] = a.sameside;
    c.jf[0][1] = a.jf1;
    c.jf[1][0] = a.jf2;
    // `oksea[1] = SearchFD(DStr, CD[0], CD[2], ...)` (`cxx:310-327`).
    let a = search_fd(base, &c, 0, 2)?;
    c.oksea[1] = a.found;
    c.i[0][2] = a.i1;
    c.i[2][0] = a.i2;
    c.p[0][2] = a.p1;
    c.p[2][0] = a.p2;
    c.face[1] = a.face.clone();
    c.sameside[1] = a.sameside;
    c.jf[0][2] = a.jf1;
    c.jf[2][0] = a.jf2;
    // `oksea[0] = SearchFD(DStr, CD[1], CD[2], ...)` (`cxx:328-339`).
    let a = search_fd(base, &c, 1, 2)?;
    c.oksea[0] = a.found;
    c.i[1][2] = a.i1;
    c.i[2][1] = a.i2;
    c.p[1][2] = a.p1;
    c.p[2][1] = a.p2;
    c.face[0] = a.face.clone();
    c.sameside[0] = a.sameside;
    c.jf[1][2] = a.jf1;
    c.jf[2][1] = a.jf2;

    // `// Analyze concavities of 3 fillets` (`cxx:340-391`).
    if c.oksea[2] && c.oksea[1] && !c.sameside[2] && !c.sameside[1] {
        c.pivot = 0;
        c.deb = 1;
        c.fin = 2;
    } else if c.oksea[2] && c.oksea[0] && !c.sameside[2] && !c.sameside[0] {
        c.pivot = 1;
        c.deb = 2;
        c.fin = 0;
    } else if c.oksea[1] && c.oksea[0] && !c.sameside[1] && !c.sameside[0] {
        c.pivot = 2;
        c.deb = 0;
        c.fin = 1;
    } else if c.oksea[0] && c.oksea[1] && c.oksea[2] {
        // 3 concavities identic (`cxx:363-387`).
        let pivot = search_pivot(c.sens, c.p, base.tol2d);
        if pivot < 0 {
            base.perform_more_three_corner(index, 3)?;
            return Ok(());
        }
        c.pivot = pivot as usize;
        c.deb = (c.pivot + 1) % 3;
        c.fin = (c.pivot + 2) % 3;
        c.pivdif = false;
        if (c.p[0][1] - c.p[0][2]).abs() <= base.tol2d
            && (c.p[1][0] - c.p[1][2]).abs() <= base.tol2d
            && (c.p[2][0] - c.p[2][1]).abs() <= base.tol2d
        {
            c.c1pointu = true;
        }
    } else {
        base.perform_more_three_corner(index, 3)?;
        return Ok(());
    }

    // `int ifacdeb, ifacfin; ifacdeb = CD[deb]->ChangeSetOfSurfData()->Value(
    // i[deb][pivot])->Index(3 - jf[deb][pivot]); ifacfin = ...; if (ifacfin !=
    // ifacdeb) { PerformMoreThreeCorner(Jndex, 3); return; }` (`cxx:393-402`).
    let ifacdeb = sd(&c.stripes[c.deb], c.i[c.deb][c.pivot])?
        .index(3 - c.jf[c.deb][c.pivot]);
    let ifacfin = sd(&c.stripes[c.fin], c.i[c.fin][c.pivot])?
        .index(3 - c.jf[c.fin][c.pivot]);
    if ifacfin != ifacdeb {
        base.perform_more_three_corner(index, 3)?;
        return Ok(());
    }
    // `if (i[pivot][deb] != i[pivot][fin]) { PerformMoreThreeCorner(Jndex, 3);
    // return; }` (`cxx:403-411`).
    if c.i[c.pivot][c.deb] != c.i[c.pivot][c.fin] {
        base.perform_more_three_corner(index, 3)?;
        return Ok(());
    }

    // `double Rdeb, Rfin, Rdp, Rfp; gp_Pnt Pdeb, Pfin, Pdp, Pfp; gp_Vec Vdeb,
    // Vfin, Vdp, Vfp;` (`cxx:412-415`). OCCT leaves `Pdp` / `Pfp` / `Vfp`
    // unused afterwards, and so does the port.

    // `if (c1pointu) { gp_Pnt pbid; gp_Vec qv[3]; double qr[3]; for (ii = 0; ii
    // <= 2; ii++) { jj = (ii + 1) % 3; kk = (ii + 2) % 3;
    // ChFi3d_ExtrSpineCarac(DStr, CD[jj], i[jj][ii], p[jj][ii], 1, sens[jj],
    // pbid, qv[jj], qr[jj]); } ... }` (`cxx:416-443`).
    if c.c1pointu {
        let mut qv = [GpVec::new(0.0, 0.0, 0.0); 3];
        let mut qr = [0.0f64; 3];
        for ii in 0..3 {
            let jj = (ii + 1) % 3;
            let (_, v, r) = extr(base, &c, jj, ii, 1)?;
            qv[jj] = v;
            qr[jj] = r;
        }
        // `for (ii = 0; ii <= 2 && !c1toric; ii++) { jj = (ii + 1) % 3; kk =
        // (ii + 2) % 3; if (ToricCorner(face[ii], qr[jj], qr[kk], qv[ii])) {
        // c1toric = true; pivot = ii; deb = jj; fin = kk; } }`
        // (`cxx:426-436`).
        for ii in 0..3 {
            if c.c1toric {
                break;
            }
            let jj = (ii + 1) % 3;
            let kk = (ii + 2) % 3;
            let Some(face) = c.face[ii].clone() else {
                continue;
            };
            if toric_corner(&face, qr[jj], qr[kk], &qv[ii]) {
                c.c1toric = true;
                c.pivot = ii;
                c.deb = jj;
                c.fin = kk;
            }
        }
        // `if (!c1toric) { c1spheric = (std::abs(qr[0] - qr[1]) < tolapp3d &&
        // std::abs(qr[0] - qr[2]) < tolapp3d); }` (`cxx:437-441`).
        if !c.c1toric {
            c.c1spheric =
                (qr[0] - qr[1]).abs() < base.tolapp3d && (qr[0] - qr[2]).abs() < base.tolapp3d;
        }
    }

    // `int jjjd = jf[deb][pivot], jjjf = jf[fin][pivot];` (`cxx:449`).
    let jjjd = c.jf[c.deb][c.pivot];
    let jjjf = c.jf[c.fin][c.pivot];
    // The four `ChFi3d_ExtrSpineCarac` probes (`cxx:450-485`). `Pdp` / `Pfp`
    // are not read afterwards in OCCT.
    let (pdeb, vdeb, rdeb) = extr(base, &c, c.deb, c.pivot, jjjd)?;
    let (pfin, vfin, rfin) = extr(base, &c, c.fin, c.pivot, jjjf)?;
    let (_, vdp, rdp) = extr(base, &c, c.pivot, c.deb, 0)?;
    let (_, _vfp, rfp) = extr(base, &c, c.pivot, c.fin, 0)?;

    // `if (!pivdif) { ... double distest = ptestdeb.Distance(ptestfin); if
    // (distest < (Rdp + Rfp) * 0.05) { filling = true; } if (distest < (Rdp +
    // Rfp) * 0.005) { c1pointu = true; } }` (`cxx:486-521`).
    if !c.pivdif {
        let (ptestdeb, _, _) = extr(base, &c, c.deb, c.pivot, c.jf[c.deb][c.fin])?;
        let (ptestfin, _, _) = extr(base, &c, c.fin, c.pivot, c.jf[c.fin][c.deb])?;
        let distest = ptestdeb.distance(&ptestfin);
        if distest < (rdp + rfp) * 0.05 {
            c.filling = true;
        }
        if distest < (rdp + rfp) * 0.005 {
            c.c1pointu = true;
        }
    }

    // `if (!c1pointu) { if (!pivdif) { c1pointu = (...) } if (std::abs(
    // p[pivot][deb] - p[pivot][fin]) <= tol2d) { c1toric =
    // ToricCorner(face[pivot], Rdeb, Rfin, Vdp); } }` (`cxx:523-534`).
    if !c.c1pointu {
        if !c.pivdif {
            c.c1pointu = (c.p[c.deb][c.pivot] - c.p[c.deb][c.fin]).abs() <= base.tol2d
                && (c.p[c.fin][c.pivot] - c.p[c.fin][c.deb]).abs() <= base.tol2d;
        }
        if (c.p[c.pivot][c.deb] - c.p[c.pivot][c.fin]).abs() <= base.tol2d {
            if let Some(face) = c.face[c.pivot].clone() {
                c.c1toric = toric_corner(&face, rdeb, rfin, &vdp);
            }
        }
    }

    // `occ::handle<ChFiDS_SurfData>& fddeb = ...; fdfin = ...; fdpiv =
    // CD[pivot]->ChangeSetOfSurfData()->ChangeValue(i[pivot][deb]);`
    // (`cxx:542-547`). The three handles are reached through
    // [`with_sd`] / [`sd`] below.

    // `TopAbs_Orientation OFac = face[pivot].Orientation(); occ::handle<
    // BRepAdaptor_Surface> Fac = new BRepAdaptor_Surface(face[pivot]);`
    // (`cxx:551-552`).
    let pivot_face = c.face[c.pivot]
        .clone()
        .ok_or_else(|| "PerformThreeCorner: the pivot has no face".to_string())?;
    let o_fac = pivot_face.0.orientation();
    let fac = SurfaceAdaptor::of_face(&pivot_face);

    // `const ChFiDS_FaceInterference& bid1 = CD[pivot]->SetOfSurfData()
    // ->Value(i[pivot][deb])->InterferenceOnS1(); ppp1 =
    // bid1.PCurveOnSurf()->Value(bid1.FirstParameter());` and `bid2` with
    // `InterferenceOnS2` / `LastParameter()` (`cxx:555-561`).
    let sd_piv = sd(&c.stripes[c.pivot], c.i[c.pivot][c.deb])?;
    let ppp1 = sd_piv
        .interference_on_s1()
        .p_curve_on_surf()
        .ok_or_else(|| "PerformThreeCorner: the pivot has no S1 p-curve".to_string())?
        .value(sd_piv.interference_on_s1().first_parameter());
    let ppp2 = sd_piv
        .interference_on_s2()
        .p_curve_on_surf()
        .ok_or_else(|| "PerformThreeCorner: the pivot has no S2 p-curve".to_string())?
        .value(sd_piv.interference_on_s2().last_parameter());
    let mut uu1 = ppp1.x();
    let mut uu2 = ppp2.x();
    let mut vv1 = ppp1.y();
    let mut vv2 = ppp2.y();

    // `GeomAdaptor_Surface gasurf(surface, uu1, uu2, vv1, vv2); GeomAbs_SurfaceType
    // styp = gasurf.GetType();` (`cxx:562-568`).
    let piv_surface = base
        .my_ds
        .ds()
        .surface(sd_piv.surf())
        .surface()
        .cloned()
        .ok_or_else(|| "PerformThreeCorner: DStr.Surface(Surf()).Surface() is null".to_string())?;
    let mut gasurf = SurfaceAdaptor::of_surface_with_range(
        Arc::clone(&piv_surface),
        uu1,
        uu2,
        vv1,
        vv2,
    );
    let styp = gasurf.get_type();
    if styp == GeomAbsSurfaceType::Cylinder {
        // `double h = vv2 - vv1; vv1 -= 0.5 * h; vv2 += 0.5 * h; gasurf.Load(...)`
        // (`cxx:569-578`).
        let h = vv2 - vv1;
        vv1 -= 0.5 * h;
        vv2 += 0.5 * h;
        gasurf = SurfaceAdaptor::of_surface_with_range(Arc::clone(&piv_surface), uu1, uu2, vv1, vv2);
    } else if styp == GeomAbsSurfaceType::Torus {
        // `double h = uu2 - uu1; uu1 -= 0.1 * h; uu2 += 0.1 * h; gasurf.Load(...)`
        // (`cxx:579-588`).
        let h = uu2 - uu1;
        uu1 -= 0.1 * h;
        uu2 += 0.1 * h;
        gasurf = SurfaceAdaptor::of_surface_with_range(Arc::clone(&piv_surface), uu1, uu2, vv1, vv2);
    } else if styp == GeomAbsSurfaceType::BezierSurface
        || styp == GeomAbsSurfaceType::BSplineSurface
    {
        // `gasurf.Load(surface);` (`cxx:590-593`).
        gasurf = SurfaceAdaptor::of_surface(Arc::clone(&piv_surface));
    }

    // `occ::handle<GeomAdaptor_Surface> Surf = new GeomAdaptor_Surface(gasurf);
    // ... bidsurf = new GeomAdaptor_Surface(Fac->GeomSurfaceOriginal());
    // IFac = new Adaptor3d_TopolTool(bidsurf); ISurf = new
    // Adaptor3d_TopolTool(Surf);` (`cxx:595-603`).
    let fac_surface = Arc::clone(
        fac.surface()
            .ok_or_else(|| "PerformThreeCorner: Fac has no surface".to_string())?,
    );
    let mut ifac = crate::brep_top_adaptor_topol_tool::TopolTool::new();
    ifac.initialize(fac_surface, &pivot_face)?;
    let mut isurf = crate::brep_top_adaptor_topol_tool::TopolTool::new();
    isurf.initialize_surface_only(Arc::clone(
        gasurf
            .surface()
            .ok_or_else(|| "PerformThreeCorner: Surf has no surface".to_string())?,
    ));
    let ifac = std::rc::Rc::new(std::cell::RefCell::new(ifac));
    let isurf = std::rc::Rc::new(std::cell::RefCell::new(isurf));

    // `occ::handle<ChFiDS_Stripe> corner = new ChFiDS_Stripe(); cornerset =
    // corner->ChangeSetOfSurfData(); cornerset = new ...HSequence; coin = new
    // ChFiDS_SurfData(); cornerset->Append(coin);` (`cxx:604-608`). OCCT
    // appends first and fills in place; nothing reads the sequence before the
    // end of this function, so the port keeps `coin` as a local.
    let mut corner = Stripe::new();
    let mut coin = SurfData::new();

    // `int choix = CD[deb]->Choix(); o1 = face[pivot].Orientation(); o2 =
    // fdpiv->Orientation(); oo1 = o1; oo2 = o2; os1 =
    // CD[deb]->OrientationOnFace1(); os2 = CD[deb]->OrientationOnFace2();`
    // (`cxx:610-617`).
    let mut choix = c.stripes[c.deb]
        .read()
        .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?
        .choix();
    let mut o1 = o_fac;
    let mut o2 = sd(&c.stripes[c.pivot], c.i[c.pivot][c.deb])?.orientation();
    let oo1 = o1;
    let oo2 = o2;
    let (os1, os2) = {
        let guard = c.stripes[c.deb]
            .read()
            .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
        (guard.orientation_on_face1(), guard.orientation_on_face2())
    };
    if c.jf[c.deb][c.fin] == 1 {
        // `choix = ChFi3d::NextSide(o1, o2, os1, os2, choix); if (sens[deb] == 1)
        // { if (choix % 2 == 1) choix++; else choix--; }` (`cxx:618-631`).
        choix = crate::chfi3d::side::next_side(&mut o1, &mut o2, os1, os2, choix);
        if c.sens[c.deb] == 1 {
            if choix % 2 == 1 {
                choix += 1;
            } else {
                choix -= 1;
            }
        }
    } else {
        // `choix = ChFi3d::NextSide(o2, o1, os1, os2, -choix); if (sens[deb] ==
        // -1) { ... }` (`cxx:633-645`).
        choix = crate::chfi3d::side::next_side(&mut o2, &mut o1, os1, os2, -choix);
        if c.sens[c.deb] == -1 {
            if choix % 2 == 1 {
                choix += 1;
            } else {
                choix -= 1;
            }
        }
    }

    // `CD[deb]->...->Interference(jf[deb][fin]).PCurveOnFace()->D1(
    // p[deb][pivot], pfac1, vfac1); CD[fin]->...->Interference(jf[fin][deb])
    // .PCurveOnFace()->D1(p[fin][pivot], pfac2, vfac2);` (`cxx:650-662`).
    let (pfac1, vfac1) = {
        let s = sd(&c.stripes[c.deb], c.i[c.deb][c.pivot])?;
        let pc = s
            .interference(c.jf[c.deb][c.fin])
            .p_curve_on_face()
            .ok_or_else(|| "PerformThreeCorner: CD[deb] has no face p-curve".to_string())?;
        pc.d1(c.p[c.deb][c.pivot])
    };
    let (pfac2, vfac2) = {
        let s = sd(&c.stripes[c.fin], c.i[c.fin][c.pivot])?;
        let pc = s
            .interference(c.jf[c.fin][c.deb])
            .p_curve_on_face()
            .ok_or_else(|| "PerformThreeCorner: CD[fin] has no face p-curve".to_string())?;
        pc.d1(c.p[c.fin][c.pivot])
    };
    // `psurf1 = CD[pivot]->...->Interference(jf[pivot][deb]).PCurveOnSurf()
    // ->Value(p[pivot][deb]); psurf2 = ...Value(p[pivot][fin]);`
    // (`cxx:663-672`).
    let psurf1 = {
        let s = sd(&c.stripes[c.pivot], c.i[c.pivot][c.deb])?;
        s.interference(c.jf[c.pivot][c.deb])
            .p_curve_on_surf()
            .ok_or_else(|| "PerformThreeCorner: the pivot has no S1 surface p-curve".to_string())?
            .value(c.p[c.pivot][c.deb])
    };
    let psurf2 = {
        let s = sd(&c.stripes[c.pivot], c.i[c.pivot][c.fin])?;
        s.interference(c.jf[c.pivot][c.fin])
            .p_curve_on_surf()
            .ok_or_else(|| "PerformThreeCorner: the pivot has no S2 surface p-curve".to_string())?
            .value(c.p[c.pivot][c.fin])
    };

    // `done = false;` (`cxx:675`).
    base.done = false;

    let frame = CornerFrame {
        o_fac,
        fac,
        surf: gasurf,
        ifac,
        isurf,
        oo1,
        oo2,
        o1,
        o2,
        choix,
        rdeb,
        rfin,
        rdp,
        rfp,
        pdeb,
        vdeb,
        pfin,
        vfin,
        pfac1,
        vfac1,
        pfac2,
        vfac2,
        psurf1,
        psurf2,
    };

    // `if (c1toric) { done = ChFiKPart_ComputeData::ComputeCorner(DStr, coin,
    // Fac, Surf, oo1, oo2, o1, o2, Rdeb, Rdp, pfac1, pfac2, psurf1, psurf2); }`
    // (`cxx:686-714`).
    if c.c1toric {
        base.done = crate::chfi_k_part::ch_fi_kpart_compute_corner_plane_cylinder(
            base.my_ds.change_ds(),
            &mut coin,
            &frame.fac,
            &frame.surf,
            frame.oo1,
            frame.oo2,
            frame.o1,
            frame.o2,
            frame.rdeb,
            frame.rdp,
            frame.pfac1,
            frame.pfac2,
            frame.psurf1,
            frame.psurf2,
        )?;
    } else if c.c1spheric {
        // `done = ChFiKPart_ComputeData::ComputeCorner(DStr, coin, Fac, Surf,
        // oo1, oo2, o1, o2, Rdp, pfac1, psurf1, psurf2);` (`cxx:715-735`).
        base.done = crate::chfi_k_part::ch_fi_kpart_compute_corner_sphere(
            base.my_ds.change_ds(),
            &mut coin,
            &frame.fac,
            &frame.surf,
            frame.oo1,
            frame.oo2,
            frame.o1,
            frame.o2,
            frame.rdp,
            frame.pfac1,
            frame.psurf1,
            frame.psurf2,
        )?;
    } else if c.c1pointu {
        // `else if (c1pointu) { filling = true; }` (`cxx:736-740`).
        c.filling = true;
    }

    // `if (!done) { ... }` (`cxx:741-1016`).
    if !base.done {
        base.done = fill::build_corner(base, shape, &mut c, &frame, &mut coin)?;
    }

    // `if (!c1pointu) { p[deb][fin] = p[deb][pivot]; p[fin][deb] =
    // p[fin][pivot]; }` (`cxx:1016-1021`).
    if !c.c1pointu {
        c.p[c.deb][c.fin] = c.p[c.deb][c.pivot];
        c.p[c.fin][c.deb] = c.p[c.fin][c.pivot];
    }

    // `if (done) { // Update of 4 Stripes and the DS ... }` (`cxx:1022-1289`).
    if base.done {
        ds_update::write_ds(base, &mut c, &mut corner, &mut coin)?;
    }

    // `// The data corners are truncated and index is updated.`
    // (`cxx:1289-1331`).
    ds_update::trim_and_append(base, &mut c, &mut corner, &coin)?;
    Ok(())
}

