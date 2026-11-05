//! `IntersectMoreCorner`'s prologue — `ChFi3d_Builder_C1.cxx:4708-4988`.
//!
//! The first stripe of the vertex (OCCT does not call `ChFi3d_SelectStripe`
//! here), the surplus `SurfData` removals, the cap-side face resolution, the
//! extension of `Fv`'s surface and the three `Update` calls.

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_geom::Surface;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::chfi_ds::SurfData;
use crate::chfi3d::builder_0_helpers::ch_fi3d_index_of_surf_data;
use crate::chfi3d::builder_c1::{
    ch_fi3d_extend_surface, cherche_face, inter_update_on_same, update_curve_on_surface,
};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{is_same, shapes_of};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::super::perform_one_corner::setup::{ds_face, edge_contains_vertex, face_edges};
use super::super::ChFi3dBuilder;
use super::MoreCorner;

/// `IntersectMoreCorner`'s prologue (`cxx:4708-4988`).
pub(super) fn setup(b: &mut ChFi3dBuilder, index: usize) -> Result<MoreCorner, String> {
    // `StrIt.Initialize(myVDataMap(Index)); stripe = StrIt.Value();`
    // (`cxx:4716-4718`).
    let Some(vtx) = b.my_v_data_map.find_key(index).cloned() else {
        return Err("IntersectMoreCorner: myVDataMap has no key for the index".to_string());
    };
    let stripe = b
        .my_v_data_map
        .get_from_index(index)
        .and_then(|l| l.first().cloned())
        .ok_or_else(|| "IntersectMoreCorner: myVDataMap(Index) is empty".to_string())?;

    // `spine = stripe->Spine(); SeqFil = stripe->ChangeSetOfSurfData()
    //  ->ChangeSequence();` (`cxx:4719-4721`).
    let spine = {
        let guard = stripe
            .read()
            .map_err(|_| "IntersectMoreCorner: poisoned stripe".to_string())?;
        guard
            .spine()
            .cloned()
            .ok_or_else(|| "IntersectMoreCorner: the stripe has no spine".to_string())?
    };
    let mut seq: Vec<Arc<SurfData>> = {
        let guard = stripe
            .read()
            .map_err(|_| "IntersectMoreCorner: poisoned stripe".to_string())?;
        guard.set_of_surf_data().to_vec()
    };

    // `num = ChFi3d_IndexOfSurfData(Vtx, stripe, sens); isfirst = (sens == 1);`
    // (`cxx:4728-4729`).
    let mut sens = 0;
    let mut num = ch_fi3d_index_of_surf_data(&vtx, &stripe, &mut sens)? as usize;
    let is_first = sens == 1;

    // `cxx:4730-4746`: the surplus parts are removed from the real sequence.
    let bad = |sd: &SurfData| sd.index_of_s1() == 0 || sd.index_of_s2() == 0;
    if is_first {
        while num < seq.len() && bad(&seq[num - 1]) {
            seq.remove(num - 1);
        }
    } else {
        while num > 1 && bad(&seq[num - 1]) {
            seq.remove(num - 1);
            num -= 1;
        }
    }
    if let Ok(mut guard) = stripe.write() {
        *guard.change_set_of_surf_data() = seq.clone();
    }

    // `occ::handle<ChFiDS_SurfData>& Fd = SeqFil.ChangeValue(num);` (`cxx:4748`).
    let fd: SurfData = num
        .checked_sub(1)
        .and_then(|i| seq.get(i))
        .map(|sd| sd.as_ref().clone())
        .ok_or_else(|| "IntersectMoreCorner: no SurfData left at the corner".to_string())?;
    // `ChFiDS_CommonPoint& CV1 = Fd->ChangeVertex(isfirst, 1);` and `CV2`
    // (`cxx:4749-4750`).
    let cv1 = fd.vertex(is_first, 1).clone();
    let cv2 = fd.vertex(is_first, 2).clone();

    // `Arcspine = isfirst ? spine->Edges(1) : spine->Edges(spine->NbEdges());`
    // (`cxx:4759-4766`).
    let arcspine = {
        let guard = spine
            .read()
            .map_err(|_| "IntersectMoreCorner: poisoned spine".to_string())?;
        let i = if is_first { 1 } else { guard.nb_edges() };
        guard.edges(i).cloned()
    }
    .ok_or_else(|| "IntersectMoreCorner: the spine has no edge".to_string())?;

    // `bool inters = true; int IFadArc = 1, IFopArc = 2;` (`cxx:4782-4783`).
    // `Fop = TopoDS::Face(DStr.Shape(Fd->Index(IFopArc)));` (`cxx:4784`) is a
    // dead store in OCCT: `cxx:4826` always overwrites it before any read, so
    // the port initializes `fop` once `IFopArc` is settled (`cxx:4825-4826`).
    // `IFopArc = 2` is likewise dead because `cxx:4821` sets it from `IFadArc`
    // before its first read.
    let mut i_fad_arc = 1;
    let mut i_fop_arc;

    // `cxx:4791-4822`: which of the two arcs holds the corner vertex.
    if !cv1.is_on_arc() && !cv2.is_on_arc() {
        return Err("Corner intersmore : no point on arc".to_string());
    } else if cv1.is_on_arc() && cv2.is_on_arc() {
        // The first `TopExp_Explorer` over `CV1.Arc()` (`cxx:4798-4804`) only
        // breaks out of its own loop and is not read again; `sur2` is the
        // `CV2.Arc()` walk (`cxx:4805-4812`).
        let sur2 = cv2
            .arc()
            .map(|a| edge_contains_vertex(a, &vtx))
            .unwrap_or(false);
        if sur2 {
            i_fad_arc = 2;
        }
    } else if cv2.is_on_arc() {
        i_fad_arc = 2;
    }
    i_fop_arc = 3 - i_fad_arc;

    // `Arcpiv = Fd->Vertex(isfirst, IFadArc).Arc(); Fad = Face(DStr.Shape(
    //  Fd->Index(IFadArc))); Fop = Face(DStr.Shape(Fd->Index(IFopArc)));`
    // (`cxx:4824-4826`).
    let mut arcpiv = fd.vertex(is_first, i_fad_arc).arc().cloned();
    let mut fad = ds_face(b.my_ds.change_ds(), fd.index(i_fad_arc));
    let mut fop = ds_face(b.my_ds.change_ds(), fd.index(i_fop_arc));

    // `cxx:4829-4836`: the face at end, first non-`Fad` face of `myEFMap`.
    let mut fv = search_face_at_end(b, arcpiv.as_ref(), fad.as_ref());

    // `cxx:4838-4847`: does the face at end contain the vertex?
    let isinface = fv
        .as_ref()
        .map(|f| {
            shapes_of(&f.0, ShapeType::Vertex)
                .iter()
                .any(|v| is_same(v, &vtx))
        })
        .unwrap_or(false);
    // `cxx:4848-4865`: swap the two arcs and search again. `Fv` is not reset by
    // OCCT, so a failed second search keeps the previous face.
    if !isinface {
        i_fad_arc = 3 - i_fad_arc;
        i_fop_arc = 3 - i_fop_arc;
        arcpiv = fd.vertex(is_first, i_fad_arc).arc().cloned();
        fad = ds_face(b.my_ds.change_ds(), fd.index(i_fad_arc));
        fop = ds_face(b.my_ds.change_ds(), fd.index(i_fop_arc));
        if let Some(found) = search_face_at_end(b, arcpiv.as_ref(), fad.as_ref()) {
            fv = Some(found);
        }
    }
    // `if (Fv.IsNull()) throw StdFail_NotDone("OneCorner : face at end is not
    //  found");` (`cxx:4867-4870`).
    let Some(fv) = fv else {
        return Err("OneCorner : face at end is not found".to_string());
    };

    // `Fv.Orientation(TopAbs_FORWARD); Fad.Orientation(TopAbs_FORWARD);`
    // (`cxx:4872-4873`).
    let fv = Face(fv.0.oriented(Orientation::Forward));
    let fad = fad.map(|f| Face(f.0.oriented(Orientation::Forward)));

    // `cxx:4876-4890`: the edge that will be extended.
    let mut arcprol = None;
    let mut o_arcprol_v = Orientation::Forward;
    for item in b.my_ve_map.get_from_key(&vtx) {
        if let Some(ap) = arcpiv.as_ref() {
            if is_same(item, ap) {
                continue;
            }
        }
        for e in face_edges(&fv) {
            if is_same(item, &e) {
                arcprol = Some(item.clone());
                o_arcprol_v = e.orientation();
                break;
            }
        }
        if arcprol.is_some() {
            break;
        }
    }
    // `if (Arcprol.IsNull()) throw StdFail_NotDone("IntersectMoreCorner: edge
    //  to be extended is not found");` (`cxx:4893-4896`).
    let Some(arcprol) = arcprol else {
        return Err("IntersectMoreCorner: edge to be extended is not found".to_string());
    };

    // `Fopbis.Orientation(TopAbs_FORWARD); cherche_face(myVFMap(Vtx), Arcprol,
    //  Fad, Fv, Fv, Fopbis); Fop.Orientation(TopAbs_FORWARD); for (ex.Init(
    //  Fopbis, TopAbs_EDGE); ...) if (Arcprol.IsSame(ex.Current()))
    //  OArcprolop = ex.Current().Orientation();` (`cxx:4900-4912`).
    let vf: Vec<TopoShape> = b.my_vf_map.get_from_key(&vtx).to_vec();
    let fad_ref = fad
        .clone()
        .ok_or_else(|| "OneCorner : face at end is not found".to_string())?;
    let fopbis = cherche_face(&vf, &Edge(arcprol.clone()), &fad_ref, &fv, &fv)?;
    let mut o_arcprol_op = Orientation::Forward;
    for e in face_edges(&fopbis) {
        if is_same(&arcprol, &e) {
            o_arcprol_op = e.orientation();
            break;
        }
    }

    // `Sface = BRep_Tool::Surface(Fv); ChFi3d_ExtendSurface(Sface, prol); tol =
    //  BRep_Tool::Tolerance(Fv); BRE.MakeFace(FFv, Sface, tol);` (`cxx:4913-4921`).
    let sface = BRepTool::face_surface_world(&fv)
        .ok_or_else(|| "IntersectMoreCorner: Fv has no surface".to_string())?;
    let mut prol = 0;
    let sface = ch_fi3d_extend_surface(&sface, &mut prol);
    // `tol` feeds `MakeFace` only; the port gives `Bs` the extended surface.
    let _tol = BRepTool::face_tolerance(&fv);
    let hbs = if prol != 0 {
        // `Bs.Initialize(FFv, false); DStr.SetNewSurface(Fv, Sface);`
        // (`cxx:4922-4926`).
        b.my_ds.change_ds().set_new_surface(&fv.0, sface.clone());
        SurfaceAdaptor::of_surface(sface)
    } else {
        // `Bs.Initialize(Fv, false);` (`cxx:4929`).
        SurfaceAdaptor::of_face(&fv)
    };
    // `Bad.Initialize(Fad); Bop.Initialize(Fop);` (`cxx:4931-4932`).
    let fop = fop.ok_or_else(|| "OneCorner : face at end is not found".to_string())?;
    let hbad = SurfaceAdaptor::of_face_restricted(&fad_ref);
    let hbop = SurfaceAdaptor::of_face_restricted(&fop);

    // `HGs = new GeomAdaptor_Surface(DStr.Surface(Fd->Surf()).Surface());`
    // (`cxx:4947-4948`).
    let hgs: Arc<dyn Surface> = b
        .my_ds
        .change_ds()
        .surface(fd.surf())
        .surface()
        .cloned()
        .ok_or_else(|| "IntersectMoreCorner: DStr.Surface(Surf()).Surface() is null".to_string())?;

    // `CPopArc` / `FiopArc` / `CPadArc` / `FiadArc` are references into `Fd`
    // (`cxx:4939-4942`); `wop` is `FiadArc.Parameter(isfirst)` (`cxx:4945`).
    let mut c_pop = fd.vertex(is_first, i_fop_arc).clone();
    let mut fi_op = fd.interference(i_fop_arc).clone();
    let mut c_pad = fd.vertex(is_first, i_fad_arc).clone();
    let mut fi_ad = fd.interference(i_fad_arc).clone();
    let mut wop = fi_ad.parameter(is_first);

    // `cxx:4955-4987`: `IntersUpdateOnSame` over the extended face at end
    // (`Update(HBs, Hc3df, ...)` then the tangent-faces retry with `Arcprol`'s
    // pcurve on `Fv`), then `Update(HBs, pced, HGs, FiadArc, CPadArc, isfirst)`.
    let c3df = b
        .my_ds
        .change_ds()
        .curve(fi_op.line_index())
        .curve()
        .cloned()
        .ok_or_else(|| "IntersectMoreCorner: FiopArc.LineIndex() is not in the DS".to_string())?;
    let hbs_surface = hbs
        .surface()
        .cloned()
        .ok_or_else(|| "IntersectMoreCorner: Bs has no surface".to_string())?;
    let mut p2dbout = GpPnt2d::zero();
    let inters = inter_update_on_same(
        &hgs,
        &hbs_surface,
        &c3df,
        &fop,
        &fv,
        &Edge(arcprol.clone()),
        &vtx,
        is_first,
        10.0 * b.tolapp3d,
        &mut fi_op,
        &mut c_pop,
        &mut p2dbout,
        &mut wop,
    )?;

    let arc = c_pad
        .arc()
        .cloned()
        .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
    let (pced, _, _) = curve_on_surface_range(&Edge(arc), &fv)
        .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
    // `Update(HBs, pced, HGs, FiadArc, CPadArc, isfirst);` — the OCCT return
    // value is not read.
    let _ = update_curve_on_surface(&hbs_surface, &pced, &hgs, &mut fi_ad, &mut c_pad, is_first);

    // The `Fd` references are written back (OCCT never copied them out).
    let mut fd = fd;
    *fd.vertex_mut(is_first, i_fop_arc) = c_pop;
    *fd.interference_mut(i_fop_arc) = fi_op;
    *fd.vertex_mut(is_first, i_fad_arc) = c_pad;
    *fd.interference_mut(i_fad_arc) = fi_ad;

    Ok(MoreCorner {
        vtx,
        stripe,
        seq,
        num,
        is_first,
        fv: Some(fv),
        fad: Some(fad_ref),
        fop: Some(fop),
        fopbis: Some(fopbis),
        arcpiv,
        arcspine,
        arcprol,
        arcprolbis: None,
        o_arcprol_v,
        o_arcprol_op,
        o_arcprolbis: Orientation::Forward,
        i_fad_arc,
        i_fop_arc,
        hbs: Some(hbs),
        hbad: Some(hbad),
        hbop: Some(hbop),
        hgs: Some(hgs),
        fd,
        p2dbout,
        wop,
        inters,
        cc: None,
        ps: None,
        pc: None,
        udeb: 0.0,
        ufin: 0.0,
        tolreached: 0.0,
        couture: false,
        edgecouture: None,
        intcouture: false,
        par1: 0.0,
        par2: 0.0,
        indpt: 0,
        icurv1: 0,
        icurv2: 0,
        curv1: None,
        curv2: None,
        box1: occt_core::bnd::BndBox::new(),
        box2: occt_core::bnd::BndBox::new(),
    })
}

/// `cxx:4829-4836` / `cxx:4857-4864`: the first face of `myEFMap(Arcpiv)` that
/// is not `Fad`.
fn search_face_at_end(
    b: &ChFi3dBuilder,
    arcpiv: Option<&TopoShape>,
    fad: Option<&Face>,
) -> Option<Face> {
    let arcpiv = arcpiv?;
    for item in b.my_ef_map.get_from_key(arcpiv) {
        // `if (!Fad.IsSame(It.Value()))` — a null `Fad` is never `IsSame`.
        let same_as_fad = fad.map(|f| is_same(item, &f.0)).unwrap_or(false);
        if !same_as_fad {
            return Face::wrap(item.clone());
        }
    }
    None
}
