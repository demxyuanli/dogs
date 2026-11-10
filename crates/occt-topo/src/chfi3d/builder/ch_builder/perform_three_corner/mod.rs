//! `ChFi3d_ChBuilder::PerformThreeCorner` - the corner shared by three chamfer
//! stripes of a vertex.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_ChBuilder_C3.cxx:202-1142`
//! (OCCT 8.0.0). The file statics of the same `.cxx` are [`statics`]:
//! `CoPlanar` (`cxx:64-88`), `BoundSurf` (`cxx:96-114`) and
//! `ComputeIntersection` (`cxx:131-197`). `SearchPivot` (`cxx:78-105` of
//! `ChFi3d_FilBuilder_C3.cxx`) already lives in
//! [`crate::chfi3d::builder_helpers::search_pivot`], `SearchFD` in
//! [`crate::chfi3d::builder_0_helpers::ch_fi3d_search_fd`], `IsInFront` in
//! [`crate::chfi3d::builder_0_helpers::ch_fi3d_is_in_front`].
//!
//! The body is split by OCCT line ranges:
//!
//! - this file - `cxx:202-572`: `issmooth`, the three stripes and their
//!   `ChFi3d_IndexOfSurfData` indices, the three `SearchFD` calls, the
//!   concavity analysis of `cxx:304-395` (with its `IsInFront` retries), the
//!   three `ChFi3d_ComputesIntPC` probes, the `c1triangle` widening of
//!   `cxx:529-547` and the `c1plan` test,
//! - [`plan`] - `cxx:573-796`: the corner carried by a plane through the first
//!   three points, with its three `ComputeIntersection` calls,
//! - [`fill`] - `cxx:798-886`: the `GeomFill_ConstrainedFilling` corner,
//! - [`ds_update`] - `cxx:888-1142`: the four stripes, the DS, and the
//!   truncation of the three `SurfData` sequences.
//!
//! Translation notes:
//! - OCCT's `done` is [`ChFi3dBuilder::done`](super::super::ChFi3dBuilder),
//! - `issmooth` is declared `bool issmooth = false;` (`cxx:207`) and never
//!   assigned, so every `if (issmooth)` arm is dead; both arms are translated
//!   and the flag stays a local, exactly as in the `.cxx` (`Tgpiv` / `Tg3` are
//!   default-constructed `gp_Vec2d`, which is the zero vector),
//! - `IFac` / `ISurf` (`cxx:437-442`) are built and never read by the body, so
//!   they are not constructed,
//! - `CD[i]` are `ChFiDS_Stripe` handles, so `CD[0] == CD[1]` of `cxx:233-249`
//!   is [`Arc::ptr_eq`],
//! - the `SurfData` OCCT holds by reference (`fddeb` / `fdfin` / `fdpiv`) are
//!   read through [`sd`](super::super::corner_sd::sd),
//! - `Standard_Failure("FD en vis a vis non trouvees")` (`cxx:399`),
//!   `"Echec intersection PCurves OnCommonFace"` (`cxx:402`),
//!   `Standard_NotImplemented("coin mutant non programme")` (`cxx:415`) and
//!   `"coin sur plusieurs faces non programme"` (`cxx:424`) are the `Err` of
//!   this port,
//! - `visavis` is an out-parameter `IsInFront` writes and the body never reads
//!   (`cxx:302`, `:324`, `:353`, `:378`); [`IsInFrontState::visavis`] carries
//!   it and is dropped.

mod ds_update;
mod fill;
mod plan;
mod statics;

use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_geom2d::curve::Curve2d;

use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_bound_surf, ch_fi3d_computes_int_pc_point, ch_fi3d_index_of_surf_data,
    ch_fi3d_is_in_front, ch_fi3d_search_fd, IsInFrontState,
};
use crate::chfi3d::builder_helpers::search_pivot;
use crate::chfi_ds::{Stripe, StripeHandle, SurfData};
use crate::shape::{Face, TopoShape};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::super::corner_sd::{data_surface, sd, sd_len};
use super::super::ChFi3dBuilder;

/// The state the four OCCT ranges share.
pub(super) struct ThreeChamfer {
    /// `CD[0..3]` (`cxx:212`).
    pub stripes: [StripeHandle; 3],
    /// `Index[0..3]` (`cxx:213`).
    pub index: [i32; 3],
    /// `sens[0..3]` (`cxx:222`).
    pub sens: [i32; 3],
    /// `i[0..3][0..3]` (`cxx:200`).
    pub i: [[i32; 3]; 3],
    /// `jf[0..3][0..3]` (`cxx:215`).
    pub jf: [[i32; 3]; 3],
    /// `p[0..3][0..3]` (`cxx:223`).
    pub p: [[f64; 3]; 3],
    /// `oksea[0..3]` (`cxx:216`).
    pub oksea: [bool; 3],
    /// `sameside[0..3]` (`cxx:216`).
    pub sameside: [bool; 3],
    /// `face[0..3]` (`cxx:215`).
    pub face: [Option<Face>; 3],
    /// `pivot` (`cxx:213`).
    pub pivot: usize,
    /// `deb` (`cxx:213`).
    pub deb: usize,
    /// `fin` (`cxx:213`).
    pub fin: usize,
    /// `CornerAllSame` (`cxx:300`).
    pub corner_all_same: bool,
    /// `c1triangle` (`cxx:224`).
    pub c1triangle: bool,
    /// `Vtx` (`cxx:209`).
    pub vtx: TopoShape,
}

/// The locals of `cxx:486-560` the two corner arms read.
pub(super) struct CornerProbe {
    /// `Fac` (`cxx:435`): `BRepAdaptor_Surface(face[pivot])`.
    pub fac: SurfaceAdaptor,
    /// `Surf` (`cxx:439-440`): `ChFi3d_BoundSurf(DStr, fdpiv, jf[pivot][deb],
    /// jf[pivot][fin])`.
    pub surf: SurfaceAdaptor,
    /// `p2d[0..4]` (`cxx:461`): `deb`, `fin`, `pivot` and the fourth point.
    pub p2d: [GpPnt2d; 4],
    /// `p3d[0..4]` (`cxx:462`): the same four points in space.
    pub p3d: [GpPnt; 4],
    /// `deru` (`cxx:520`) and `derv` (`cxx:520`) are consumed by `norpl`
    /// right where `Fac->D1` produces them, so they are not carried here.
    /// `norpl` (`cxx:523`).
    pub norpl: GpVec,
    /// `Tgpiv` (`cxx:500`): the face tangent at `p2d[pivot]`, only written by
    /// the dead `issmooth` arm.
    pub tgpiv: GpVec2d,
    /// `Tg3` (`cxx:500`): the face tangent at `p2d[3]`, only written by the
    /// dead `issmooth` arm.
    pub tg3: GpVec2d,
}

/// What the `c1plan` arm answers to [`ds_update`].
pub(super) struct PlanCorner {
    /// `Icf` (`cxx:563`).
    pub icf: i32,
    /// `Icl` (`cxx:563`).
    pub icl: i32,
    /// `debpc1` (`cxx:564`).
    pub debpc1: Option<Arc<dyn Curve2d>>,
    /// `finpc1` (`cxx:564`).
    pub finpc1: Option<Arc<dyn Curve2d>>,
}

/// `ChFi3d_IsInFront(DStr, CD[a], CD[b], ..., face[k], sameside[k],
/// jf[a][b], jf[b][a], visavis, Vtx, false, true)` (`cxx:312-329`,
/// `:341-357`, `:366-383`): the retry `PerformThreeCorner` runs when the third
/// `SearchFD` failed.
fn is_in_front(
    base: &ChFi3dBuilder,
    c: &mut ThreeChamfer,
    a: usize,
    b: usize,
    k: usize,
) -> Result<bool, String> {
    let ca = c.stripes[a]
        .read()
        .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
    let cb = c.stripes[b]
        .read()
        .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
    let mut st = IsInFrontState::new(c.p[a][b], c.p[b][a]);
    st.face = c.face[k].clone();
    st.sameside = c.sameside[k];
    st.jf1 = c.jf[a][b];
    st.jf2 = c.jf[b][a];
    ch_fi3d_is_in_front(
        &base.my_ds,
        &ca,
        &cb,
        c.i[a][b],
        c.i[b][a],
        c.sens[a],
        c.sens[b],
        &mut st,
        Some(&c.vtx),
        false,
        true,
    )?;
    drop(ca);
    drop(cb);
    c.p[a][b] = st.p1;
    c.p[b][a] = st.p2;
    c.face[k] = st.face;
    c.sameside[k] = st.sameside;
    c.jf[a][b] = st.jf1;
    c.jf[b][a] = st.jf2;
    Ok(st.ok)
}

/// `SearchFD(DStr, CD[a], CD[b], sens[a], sens[b], ..., Index[a], Index[b],
/// face[k], sameside[k], jf[a][b], jf[b][a])`.
fn search_fd(
    base: &ChFi3dBuilder,
    c: &mut ThreeChamfer,
    a: usize,
    b: usize,
) -> Result<bool, String> {
    let ga = c.stripes[a]
        .read()
        .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
    let gb = c.stripes[b]
        .read()
        .map_err(|_| "PerformThreeCorner: poisoned stripe".to_string())?;
    let answer = ch_fi3d_search_fd(&base.my_ds, &ga, &gb, c.sens[a], c.sens[b], c.index[a], c.index[b])?;
    drop(ga);
    drop(gb);
    c.i[a][b] = answer.i1;
    c.i[b][a] = answer.i2;
    c.p[a][b] = answer.p1;
    c.p[b][a] = answer.p2;
    c.face[3 - a - b] = answer.face;
    c.sameside[3 - a - b] = answer.sameside;
    c.jf[a][b] = answer.jf1;
    c.jf[b][a] = answer.jf2;
    Ok(answer.found)
}

/// `PerformThreeCorner(Index)` (`cxx:202-1142`).
pub(super) fn perform_three_corner(
    base: &mut ChFi3dBuilder,
    index: usize,
) -> Result<(), String> {
    // `bool issmooth = false;` (`cxx:207`).
    let issmooth = false;

    // `const TopoDS_Vertex& Vtx = myVDataMap.FindKey(Jndex);` (`cxx:209`).
    let vtx = base
        .my_v_data_map
        .find_key(index)
        .cloned()
        .ok_or_else(|| format!("PerformThreeCorner: myVDataMap has no key {index}"))?;

    // `for (It.Initialize(myVDataMap(Jndex)), ii = 0; It.More() && ii < 3;
    // It.Next(), ii++) { Index[ii] = ChFi3d_IndexOfSurfData(...); CD[ii] =
    // It.Value(); }` (`cxx:227-231`).
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
    let mut c = ThreeChamfer {
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
        corner_all_same: false,
        c1triangle: false,
        vtx,
    };
    for ii in 0..3 {
        let mut sens = 0;
        c.index[ii] = ch_fi3d_index_of_surf_data(&c.vtx, &c.stripes[ii], &mut sens)?;
        c.sens[ii] = sens;
    }

    // `if (CD[0] == CD[1]) { Index[1] = ...Length(); sens[1] = -1; } else if
    // (CD[1] == CD[2]) { ... } else if (CD[0] == CD[2]) { ... }`
    // (`cxx:233-249`).
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

    // `oksea[2] = ChFi3d_SearchFD(DStr, CD[0], CD[1], ...)` (`cxx:250-264`),
    // `oksea[1] = ... CD[0], CD[2] ...` (`cxx:265-279`), `oksea[0] = ...
    // CD[1], CD[2] ...` (`cxx:280-297`).
    c.oksea[2] = search_fd(base, &mut c, 0, 1)?;
    c.oksea[1] = search_fd(base, &mut c, 0, 2)?;
    c.oksea[0] = search_fd(base, &mut c, 1, 2)?;

    // `bool CornerAllSame = false; bool okinter = true; bool visavis;`
    // (`cxx:300-302`).
    let mut okinter = true;

    // `Analyse des concavites des 3 chanfreins` (`cxx:304-395`).
    if c.oksea[2] && c.oksea[1] && !c.sameside[2] && !c.sameside[1] {
        // `pivot = 0; deb = 1; fin = 2; if (!oksea[0]) { okinter =
        // ChFi3d_IsInFront(DStr, CD[1], CD[2], ...); }` (`cxx:307-329`).
        c.pivot = 0;
        c.deb = 1;
        c.fin = 2;
        if !c.oksea[0] {
            okinter = is_in_front(base, &mut c, 1, 2, 0)?;
        }
    } else if c.oksea[2] && c.oksea[0] && !c.sameside[2] && !c.sameside[0] {
        // `pivot = 1; deb = 2; fin = 0; if (!oksea[1]) { ... }`
        // (`cxx:331-357`).
        c.pivot = 1;
        c.deb = 2;
        c.fin = 0;
        if !c.oksea[1] {
            okinter = is_in_front(base, &mut c, 0, 2, 1)?;
        }
    } else if c.oksea[1] && c.oksea[0] && !c.sameside[1] && !c.sameside[0] {
        // `pivot = 2; deb = 0; fin = 1; if (!oksea[2]) { ... }`
        // (`cxx:358-384`).
        c.pivot = 2;
        c.deb = 0;
        c.fin = 1;
        if !c.oksea[2] {
            okinter = is_in_front(base, &mut c, 0, 1, 2)?;
        }
    } else if c.oksea[0] && c.oksea[1] && c.oksea[2] {
        // `// 3 concavites identiques.` (`cxx:385-395`).
        let pivot = search_pivot(c.sens, c.p, base.tol2d);
        // `if (pivot < 0) { // on prend un pivot au hasard!!!!!!!!!!!!!!!
        // pivot = 0; }` (`cxx:388-391`).
        c.pivot = if pivot < 0 { 0 } else { pivot as usize };
        c.deb = (c.pivot + 1) % 3;
        c.fin = (c.pivot + 2) % 3;
        c.corner_all_same = true;
    } else {
        // `throw Standard_Failure("FD en vis a vis non trouvees");`
        // (`cxx:397-400`).
        return Err("FD en vis a vis non trouvees".to_string());
    }

    // `if (!okinter) { throw Standard_Failure("Echec intersection PCurves
    // OnCommonFace"); }` (`cxx:400-405`).
    if !okinter {
        return Err("Echec intersection PCurves OnCommonFace".to_string());
    }

    // `if (i[pivot][deb] != i[pivot][fin]) { throw
    // Standard_NotImplemented("coin mutant non programme"); }` (`cxx:413-416`).
    if c.i[c.pivot][c.deb] != c.i[c.pivot][c.fin] {
        return Err(
            "coin mutant non programme (ChFi3d_ChBuilder_C3.cxx:413-416)".to_string(),
        );
    }
    // `if (oksea[pivot] && (i[deb][pivot] != i[deb][fin] || i[fin][pivot] !=
    // i[fin][deb])) { throw Standard_NotImplemented("coin sur plusieurs faces
    // non programme"); }` (`cxx:422-425`).
    if c.oksea[c.pivot]
        && (c.i[c.deb][c.pivot] != c.i[c.deb][c.fin]
            || c.i[c.fin][c.pivot] != c.i[c.fin][c.deb])
    {
        return Err(
            "coin sur plusieurs faces non programme (ChFi3d_ChBuilder_C3.cxx:422-425)".to_string(),
        );
    }

    // `occ::handle<ChFiDS_SurfData>& fddeb = CD[deb]->...->ChangeValue(
    // i[deb][pivot]); fdfin = ...; fdpiv = CD[pivot]->...->ChangeValue(
    // i[pivot][deb]);` (`cxx:427-431`). The three handles are read through
    // [`sd`] and written in place by [`ds_update`].
    let fddeb = sd(&c.stripes[c.deb], c.i[c.deb][c.pivot])?;
    let fdfin = sd(&c.stripes[c.fin], c.i[c.fin][c.pivot])?;
    let fdpiv = sd(&c.stripes[c.pivot], c.i[c.pivot][c.deb])?;

    // `occ::handle<BRepAdaptor_Surface> Fac = new BRepAdaptor_Surface(
    // face[pivot]);` (`cxx:435`).
    let pivot_face = c.face[c.pivot]
        .clone()
        .ok_or_else(|| "PerformThreeCorner: the pivot has no face".to_string())?;
    let fac = SurfaceAdaptor::of_face(&pivot_face);
    // `Surf = ChFi3d_BoundSurf(DStr, fdpiv, jf[pivot][deb], jf[pivot][fin]);`
    // (`cxx:439-440`).
    let surf = ch_fi3d_bound_surf(
        base.my_ds.ds(),
        &fdpiv,
        c.jf[c.pivot][c.deb],
        c.jf[c.pivot][c.fin],
    )?;

    // `corner = new ChFiDS_Stripe(); cornerset = new ...HSequence; coin = new
    // ChFiDS_SurfData(); cornerset->Append(coin);` (`cxx:445-449`). The
    // sequence is empty before the append, so the port keeps `coin` local and
    // appends it in [`ds_update::trim_and_append`].
    let mut corner = Stripe::new();
    let mut coin = SurfData::new();

    // `HSdeb = new GeomAdaptor_Surface(GeomAdaptor_Surface(DStr.Surface(
    // fddeb->Surf()).Surface()));` and the same for `HSfin` / `HSpiv`
    // (`cxx:453-459`).
    let hsdeb = SurfaceAdaptor::of_surface(data_surface(base.my_ds.ds(), &fddeb)?);
    let hsfin = SurfaceAdaptor::of_surface(data_surface(base.my_ds.ds(), &fdfin)?);
    let hspiv = SurfaceAdaptor::of_surface(data_surface(base.my_ds.ds(), &fdpiv)?);

    // `gp_Pnt2d p2d[4]; gp_Pnt p3d[4], PSom;` (`cxx:461-462`).
    let mut p2d = [GpPnt2d::new(0.0, 0.0); 4];
    let mut p3d = [GpPnt::new(0.0, 0.0, 0.0); 4];

    // `ChFi3d_ComputesIntPC(fdpiv->Interference(jf[pivot][deb]),
    // fddeb->Interference(jf[deb][pivot]), HSpiv, HSdeb, p[pivot][deb],
    // p[deb][pivot], p3d[fin]);` (`cxx:464-469`). OCCT's third point is the
    // pivot / deb intersection and it lands in `p3d[fin]`; kept as written.
    let mut u1 = c.p[c.pivot][c.deb];
    let mut u2 = c.p[c.deb][c.pivot];
    p3d[c.fin] = ch_fi3d_computes_int_pc_point(
        fdpiv.interference(c.jf[c.pivot][c.deb]),
        fddeb.interference(c.jf[c.deb][c.pivot]),
        &hspiv,
        &hsdeb,
        &mut u1,
        &mut u2,
    )?;
    c.p[c.pivot][c.deb] = u1;
    c.p[c.deb][c.pivot] = u2;

    // `ChFi3d_ComputesIntPC(fdpiv->Interference(jf[pivot][fin]),
    // fdfin->Interference(jf[fin][pivot]), HSpiv, HSfin, p[pivot][fin],
    // p[fin][pivot], p3d[deb]);` (`cxx:471-477`).
    let mut u1 = c.p[c.pivot][c.fin];
    let mut u2 = c.p[c.fin][c.pivot];
    p3d[c.deb] = ch_fi3d_computes_int_pc_point(
        fdpiv.interference(c.jf[c.pivot][c.fin]),
        fdfin.interference(c.jf[c.fin][c.pivot]),
        &hspiv,
        &hsfin,
        &mut u1,
        &mut u2,
    )?;
    c.p[c.pivot][c.fin] = u1;
    c.p[c.fin][c.pivot] = u2;

    // `ChFi3d_ComputesIntPC(fddeb->Interference(jf[deb][fin]),
    // fdfin->Interference(jf[fin][deb]), HSdeb, HSfin, p[deb][fin],
    // p[fin][deb], PSom);` (`cxx:478-484`).
    let mut u1 = c.p[c.deb][c.fin];
    let mut u2 = c.p[c.fin][c.deb];
    let psom = ch_fi3d_computes_int_pc_point(
        fddeb.interference(c.jf[c.deb][c.fin]),
        fdfin.interference(c.jf[c.fin][c.deb]),
        &hsdeb,
        &hsfin,
        &mut u1,
        &mut u2,
    )?;
    c.p[c.deb][c.fin] = u1;
    c.p[c.fin][c.deb] = u2;

    // `gp_Vec2d Tg3, Tgpiv;` (`cxx:500`): `issmooth` is false, so both stay
    // the null vector the C++ default constructor makes.
    let mut tg3 = GpVec2d::new(0.0, 0.0);
    let mut tgpiv = GpVec2d::new(0.0, 0.0);
    // `if (issmooth) { ...->D1(p[deb][pivot], p2d[pivot], Tgpiv); ...
    // ->D1(p[fin][pivot], p2d[3], Tg3); } else { p2d[pivot] = ...->Value(
    // p[deb][pivot]); p2d[3] = ...->Value(p[fin][pivot]); }` (`cxx:504-515`).
    if issmooth {
        let (uv, d1) = fddeb
            .interference(c.jf[c.deb][c.fin])
            .p_curve_on_face()
            .ok_or_else(|| "PerformThreeCorner: CD[deb] has no face p-curve".to_string())?
            .d1(c.p[c.deb][c.pivot]);
        p2d[c.pivot] = uv;
        tgpiv = d1;
        let (uv, d1) = fdfin
            .interference(c.jf[c.fin][c.deb])
            .p_curve_on_face()
            .ok_or_else(|| "PerformThreeCorner: CD[fin] has no face p-curve".to_string())?
            .d1(c.p[c.fin][c.pivot]);
        p2d[3] = uv;
        tg3 = d1;
    } else {
        p2d[c.pivot] = fddeb
            .interference(c.jf[c.deb][c.fin])
            .p_curve_on_face()
            .ok_or_else(|| "PerformThreeCorner: CD[deb] has no face p-curve".to_string())?
            .value(c.p[c.deb][c.pivot]);
        p2d[3] = fdfin
            .interference(c.jf[c.fin][c.deb])
            .p_curve_on_face()
            .ok_or_else(|| "PerformThreeCorner: CD[fin] has no face p-curve".to_string())?
            .value(c.p[c.fin][c.pivot]);
    }

    // `p2d[fin] = fdpiv->Interference(jf[pivot][deb]).PCurveOnSurf()->Value(
    // p[pivot][deb]); p2d[deb] = fdpiv->Interference(jf[pivot][fin])
    // .PCurveOnSurf()->Value(p[pivot][fin]);` (`cxx:516-518`).
    p2d[c.fin] = fdpiv
        .interference(c.jf[c.pivot][c.deb])
        .p_curve_on_surf()
        .ok_or_else(|| "PerformThreeCorner: the pivot has no S1 surface p-curve".to_string())?
        .value(c.p[c.pivot][c.deb]);
    p2d[c.deb] = fdpiv
        .interference(c.jf[c.pivot][c.fin])
        .p_curve_on_surf()
        .ok_or_else(|| "PerformThreeCorner: the pivot has no S2 surface p-curve".to_string())?
        .value(c.p[c.pivot][c.fin]);

    // `Fac->D1(p2d[pivot].X(), p2d[pivot].Y(), p3d[pivot], deru, derv);
    // gp_Vec norpl = deru.Crossed(derv); p3d[3] = Fac->Value(p2d[3].X(),
    // p2d[3].Y());` (`cxx:522-527`).
    let (pt, deru, derv) = fac
        .d1(p2d[c.pivot].x(), p2d[c.pivot].y())
        .ok_or_else(|| "PerformThreeCorner: Fac->D1 is not available".to_string())?;
    p3d[c.pivot] = pt;
    let norpl = deru.crossed(&derv);
    p3d[3] = fac
        .value(p2d[3].x(), p2d[3].y())
        .ok_or_else(|| "PerformThreeCorner: Fac->Value is not available".to_string())?;

    // `double DistMin = (p3d[3]).Distance(p3d[fin]); double DistTmp =
    // (p3d[pivot]).Distance(p3d[deb]); double DistDebFin =
    // (p3d[pivot]).Distance(p3d[3]); if (DistTmp > DistMin) { DistMin =
    // DistTmp; }` (`cxx:529-537`).
    let mut dist_min = p3d[3].distance(&p3d[c.fin]);
    let dist_tmp = p3d[c.pivot].distance(&p3d[c.deb]);
    let dist_deb_fin = p3d[c.pivot].distance(&p3d[3]);
    if dist_tmp > dist_min {
        dist_min = dist_tmp;
    }
    // `if (CornerAllSame) { c1triangle = (DistDebFin < 0.3 * DistMin); }`
    // (`cxx:542-545`).
    if c.corner_all_same {
        c.c1triangle = dist_deb_fin < 0.3 * dist_min;
    }
    // `if (c1triangle) { p3d[pivot] = PSom; }` (`cxx:547-551`).
    if c.c1triangle {
        p3d[c.pivot] = psom;
    }

    // `bool c1plan = c1triangle; gp_Vec v1(p3d[pivot], p3d[deb]); gp_Vec
    // v2(p3d[pivot], p3d[fin]); gp_Vec nor = v1.Crossed(v2); done = false;`
    // (`cxx:556-561`).
    let v1 = GpVec::from_pnts(&p3d[c.pivot], &p3d[c.deb]);
    let v2 = GpVec::from_pnts(&p3d[c.pivot], &p3d[c.fin]);
    let nor = v1.crossed(&v2);
    base.done = false;

    let probe = CornerProbe {
        fac,
        surf,
        p2d,
        p3d,
        norpl,
        tgpiv,
        tg3,
    };

    // `if (!c1triangle) { c1plan = CoPlanar(p3d[0], p3d[1], p3d[2], p3d[3]); }`
    // (`cxx:566-570`).
    let c1plan = if c.c1triangle {
        true
    } else {
        statics::co_planar(&probe.p3d[0], &probe.p3d[1], &probe.p3d[2], &probe.p3d[3])
    };

    let mut plan = None;
    if c1plan {
        // `if (c1plan) { ... done = true; }` (`cxx:573-796`).
        let built = plan::build_plan_corner(base, &c, &mut coin, &probe, &nor, &v1)?;
        base.done = true;
        plan = Some(built);
    } else {
        // `else { // !c1plan ... done = CompleteData(...); }` (`cxx:797-886`).
        base.done = fill::build_corner(base, &c, &probe, &mut coin)?;
    }

    // `if (done) { ... }` (`cxx:889-1099`).
    if base.done {
        ds_update::write_ds(base, &mut c, &mut corner, &mut coin, plan.as_ref())?;
    }

    // `// On tronque les corners data et met a jour les index.`
    // (`cxx:1101-1142`).
    ds_update::trim_and_append(base, &mut c, &mut corner, &coin)?;
    Ok(())
}
