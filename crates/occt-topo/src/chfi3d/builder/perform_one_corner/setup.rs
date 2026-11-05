//! `PerformOneCorner` prologue — `ChFi3d_Builder_C1.cxx:611-931`.
//!
//! Stripe selection, the surplus `SurfData` removals, the `ChFiDS_OnSame` face
//! resolution (`Arcpiv` / `Fad` / `Fv` / `Arcprol`), the `FindFace` alternative,
//! and the `Update` / `IntersUpdateOnSame` calls that re-locate the corner
//! points before the surfaces are intersected.

use std::sync::Arc;

use occt_geom::Surface;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::chfi_ds::{State, StripeHandle, SurfData};
use crate::chfi3d::builder_0_helpers::ch_fi3d_index_of_surf_data;
use crate::chfi3d::builder_c1::{
    ch_fi3d_extend_surface, ch_fi3d_select_stripe, inter_update_on_same, update_curve_on_surface,
    StripeIterator,
};
use crate::chfi3d::builder_helpers::edge_state;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{edge_vertices, is_same, shapes_of};
use crate::topopebrep::ds::DataStructure;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;
use crate::boptools_2d::curve_on_surface_range;

use super::super::ChFi3dBuilder;
use super::{OneCorner, OneCornerExit};

/// `DStr.Shape(I)` as a face (`TopoDS::Face(DStr.Shape(I))`).
pub(crate) fn ds_face(ds: &DataStructure, i: i32) -> Option<Face> {
    ds.shape(i, false).map(Face)
}

/// `for (ex.Init(E, TopAbs_VERTEX); ...) if (Vtx.IsSame(ex.Current()))` for an
/// edge: the two vertices `TopExp_Explorer` reports, either of which may be
/// absent on a degenerate edge.
pub(crate) fn edge_contains_vertex(shape: &TopoShape, vtx: &TopoShape) -> bool {
    let e = Edge(shape.clone());
    let (v1, v2) = edge_vertices(&e);
    v1.map(|v| is_same(&v.0, vtx)).unwrap_or(false)
        || v2.map(|v| is_same(&v.0, vtx)).unwrap_or(false)
}

/// `TopExp_Explorer(F.Oriented(TopAbs_FORWARD), TopAbs_EDGE)`: the face's wires
/// with the face orientation composed in, then each wire's edges with the wire
/// orientation composed in. That is the orientation `ex.Current()` reports.
pub(crate) fn face_edges(f: &Face) -> Vec<TopoShape> {
    let mut out = Vec::new();
    for wire in crate::iterator::cumulated_children(&f.0) {
        if wire.shape_type() != ShapeType::Wire {
            continue;
        }
        out.extend(crate::iterator::cumulated_children(&wire));
    }
    out
}

/// `PerformOneCorner`'s prologue. `Ok(None)` when `ChFi3d_SelectStripe`
/// answers false (`cxx:624-627`, nothing was modified); otherwise the locals the
/// rest of the body needs, with [`OneCorner::exit`] set when OCCT hands the
/// corner over to `PerformIntersectionAtEnd`.
pub(super) fn setup(
    b: &mut ChFi3dBuilder,
    index: usize,
    prepare_on_same: bool,
) -> Result<Option<OneCorner>, String> {
    // `const TopoDS_Vertex& Vtx = myVDataMap.FindKey(Index);` (`cxx:619`).
    let Some(vtx) = b.my_v_data_map.find_key(index).cloned() else {
        return Err("PerformOneCorner: myVDataMap has no key for the index".to_string());
    };

    // `StrIt.Initialize(myVDataMap(Index)); if (!ChFi3d_SelectStripe(StrIt, Vtx,
    // thePrepareOnSame)) return;` (`cxx:622-627`).
    let stripes: Vec<StripeHandle> = b
        .my_v_data_map
        .get_from_index(index)
        .unwrap_or(&[])
        .to_vec();
    let mut it = StripeIterator::new(&stripes);
    if !ch_fi3d_select_stripe(&mut it, &vtx, prepare_on_same)? {
        return Ok(None);
    }
    let Some(stripe) = it.value().cloned() else {
        return Ok(None);
    };

    // `const occ::handle<ChFiDS_Spine>& spine = stripe->Spine();` (`cxx:628`).
    let spine = {
        let guard = stripe
            .read()
            .map_err(|_| "PerformOneCorner: poisoned stripe".to_string())?;
        guard
            .spine()
            .cloned()
            .ok_or_else(|| "PerformOneCorner: the stripe has no spine".to_string())?
    };

    // `SeqFil = stripe->ChangeSetOfSurfData()->ChangeSequence();` (`cxx:629-630`).
    let mut seq: Vec<Arc<SurfData>> = {
        let guard = stripe
            .read()
            .map_err(|_| "PerformOneCorner: poisoned stripe".to_string())?;
        guard.set_of_surf_data().to_vec()
    };

    // `num = ChFi3d_IndexOfSurfData(Vtx, stripe, sens); isfirst = (sens == 1);`
    // (`cxx:635-637`). `ChFi3d_IndexOfSurfData` raises in its last arm.
    let mut sens = 0;
    let mut num = ch_fi3d_index_of_surf_data(&vtx, &stripe, &mut sens)? as usize;
    let is_first = sens == 1;

    // `cxx:638-652`: remove the surplus parts at the corner, in place, so the
    // early `PerformIntersectionAtEnd` sees the shortened sequence.
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

    // `occ::handle<ChFiDS_SurfData>& Fd = SeqFil.ChangeValue(num);` (`cxx:654`).
    let Some(fd) = num.checked_sub(1).and_then(|i| seq.get(i)).cloned() else {
        return Err("PerformOneCorner: no SurfData left at the corner".to_string());
    };
    let fd: SurfData = fd.as_ref().clone();

    // `stat = isfirst ? spine->FirstStatus() : spine->LastStatus();
    //  bool onsame = (stat == ChFiDS_OnSame);` (`cxx:663-672`).
    let stat = {
        let guard = spine
            .read()
            .map_err(|_| "PerformOneCorner: poisoned spine".to_string())?;
        if is_first {
            guard.first_status()
        } else {
            guard.last_status()
        }
    };
    let onsame = stat == State::OnSame;

    // `Arcspine = isfirst ? spine->Edges(1) : spine->Edges(spine->NbEdges());`
    // (`cxx:675-681`).
    let arcspine = {
        let guard = spine
            .read()
            .map_err(|_| "PerformOneCorner: poisoned spine".to_string())?;
        let i = if is_first { 1 } else { guard.nb_edges() };
        guard.edges(i).cloned()
    };
    let Some(arcspine) = arcspine else {
        return Err("PerformOneCorner: the spine has no edge".to_string());
    };

    // `Fop = TopoDS::Face(DStr.Shape(Fd->Index(IFopArc)));` (`cxx:699`), with
    // `IFopArc = 2` (`cxx:698`).
    let fop0 = b
        .my_ds
        .change_ds()
        .shape(fd.index(2), false)
        .map(Face);

    let mut oc = OneCorner {
        vtx: vtx.clone(),
        stripe: stripe.clone(),
        seq,
        num,
        is_first,
        onsame,
        fv: None,
        fop: fop0,
        fad: None,
        arcpiv: None,
        arcspine: Some(arcspine.clone()),
        arcprol: None,
        o_arcprol_v: Orientation::Forward,
        o_arcprol_op: Orientation::Forward,
        i_fad_arc: 1,
        i_fop_arc: 2,
        hbs: None,
        hbad: None,
        hbop: None,
        hgs: None,
        p2dbout: occt_core::gp::GpPnt2d::zero(),
        wop: 0.0,
        inters: true,
        fd,
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
        exit: None,
    };

    if onsame {
        // `cxx:707-711`.
        let cv1 = oc.fd.vertex(is_first, 1).clone();
        let cv2 = oc.fd.vertex(is_first, 2).clone();
        if !cv1.is_on_arc() && !cv2.is_on_arc() {
            return Err("Corner OnSame : no point on arc".to_string());
        } else if cv1.is_on_arc() && cv2.is_on_arc() {
            // `cxx:715-731`: does the corner vertex lie on each arc?
            let sur1 = cv1
                .arc()
                .map(|a| edge_contains_vertex(a, &vtx))
                .unwrap_or(false);
            let sur2 = cv2
                .arc()
                .map(|a| edge_contains_vertex(a, &vtx))
                .unwrap_or(false);
            if sur1 && sur2 {
                // `cxx:733-741`.
                let e0 = Edge(cv1.arc().cloned().unwrap_or_else(|| arcspine.clone()));
                let e1 = Edge(cv2.arc().cloned().unwrap_or_else(|| arcspine.clone()));
                let e2 = Edge(arcspine.clone());
                if edge_state([&e0, &e1, &e2], &b.my_ef_map)? != State::OnDiff {
                    oc.i_fad_arc = 2;
                }
            } else if sur2 {
                oc.i_fad_arc = 2;
            }
        } else if cv2.is_on_arc() {
            oc.i_fad_arc = 2;
        }
        // `cxx:751`.
        oc.i_fop_arc = 3 - oc.i_fad_arc;

        // `cxx:753-756`.
        oc.arcpiv = oc.fd.vertex(is_first, oc.i_fad_arc).arc().cloned();
        oc.fad = ds_face(b.my_ds.change_ds(), oc.fd.index(oc.i_fad_arc));
        oc.fop = ds_face(b.my_ds.change_ds(), oc.fd.index(oc.i_fop_arc));

        // `cxx:757-764`: the face at end, last non-`Fad` face of `myEFMap`.
        let mut fv: Option<Face> = None;
        if let (Some(arpiv), Some(fad)) = (oc.arcpiv.clone(), oc.fad.clone()) {
            for item in b.my_ef_map.get_from_key(&arpiv) {
                if !is_same(item, &fad.0) {
                    fv = Face::wrap(item.clone());
                }
            }
        }

        // `cxx:766-776`: does the face at end contain the vertex?
        let mut isinface = false;
        if let Some(f) = fv.as_ref() {
            for v in shapes_of(&f.0, ShapeType::Vertex) {
                if is_same(&v, &vtx) {
                    isinface = true;
                    break;
                }
            }
        }
        // `cxx:777-793`: swap the two arcs when it does not.
        if !isinface
            && oc
                .fd
                .vertex(is_first, 3 - oc.i_fad_arc)
                .is_on_arc()
        {
            oc.i_fad_arc = 3 - oc.i_fad_arc;
            oc.i_fop_arc = 3 - oc.i_fop_arc;
            oc.arcpiv = oc.fd.vertex(is_first, oc.i_fad_arc).arc().cloned();
            oc.fad = ds_face(b.my_ds.change_ds(), oc.fd.index(oc.i_fad_arc));
            oc.fop = ds_face(b.my_ds.change_ds(), oc.fd.index(oc.i_fop_arc));
            let mut again: Option<Face> = None;
            if let (Some(arpiv), Some(fad)) = (oc.arcpiv.clone(), oc.fad.clone()) {
                for item in b.my_ef_map.get_from_key(&arpiv) {
                    if !is_same(item, &fad.0) {
                        again = Face::wrap(item.clone());
                    }
                }
            }
            fv = again;
        }

        // `if (Fv.IsNull()) throw StdFail_NotDone("OneCorner : face at end not
        // found");` (`cxx:795-798`).
        let Some(fv) = fv else {
            return Err("OneCorner : face at end not found".to_string());
        };
        // `Fv.Orientation(TopAbs_FORWARD);` and the same for `Fad`, `Fop`
        // (`cxx:801-803`).
        let fv = Face(fv.0.oriented(Orientation::Forward));
        oc.fad = oc
            .fad
            .clone()
            .map(|f| Face(f.0.oriented(Orientation::Forward)));
        oc.fop = oc
            .fop
            .clone()
            .map(|f| Face(f.0.oriented(Orientation::Forward)));

        // `cxx:805-819`: the edge that will be extended.
        for item in b.my_ve_map.get_from_key(&vtx) {
            if let Some(arpiv) = oc.arcpiv.as_ref() {
                if is_same(item, arpiv) {
                    continue;
                }
            }
            for e in face_edges(&fv) {
                if is_same(item, &e) {
                    oc.arcprol = Some(item.clone());
                    oc.o_arcprol_v = e.orientation();
                    break;
                }
            }
            if oc.arcprol.is_some() {
                break;
            }
        }
        // `if (Arcprol.IsNull()) { PerformIntersectionAtEnd(Index); return; }`
        // (`cxx:821-825`).
        if oc.arcprol.is_none() {
            oc.exit = Some(OneCornerExit::IntersectionAtEnd);
            return Ok(Some(oc));
        }
        // `cxx:826-833`: the orientation of `Arcprol` on `Fop`.
        if let Some(fop) = oc.fop.as_ref() {
            for e in face_edges(fop) {
                if let Some(ap) = oc.arcprol.as_ref() {
                    if is_same(ap, &e) {
                        oc.o_arcprol_op = e.orientation();
                        break;
                    }
                }
            }
        }

        // `Sface = BRep_Tool::Surface(Fv); ChFi3d_ExtendSurface(Sface, prol);
        //  tol = BRep_Tool::Tolerance(Fv); BRE.MakeFace(FFv, Sface, tol);`
        // (`cxx:835-843`).
        let sface = BRepTool::face_surface_world(&fv)
            .ok_or_else(|| "PerformOneCorner: Fv has no surface".to_string())?;
        let mut prol = 0;
        let sface = ch_fi3d_extend_surface(&sface, &mut prol);
        // `tol` only feeds `MakeFace`, which exists to hand `Bs` the extended
        // surface; the port gives `Bs` that surface directly below.
        let _tol = BRepTool::face_tolerance(&fv);
        if prol != 0 {
            // `Bs.Initialize(FFv, false); DStr.SetNewSurface(Fv, Sface);`
            // (`cxx:845-848`).
            oc.hbs = Some(SurfaceAdaptor::of_surface(sface.clone()));
            b.my_ds.change_ds().set_new_surface(&fv.0, sface);
        } else {
            // `Bs.Initialize(Fv, false);` (`cxx:850-851`).
            oc.hbs = Some(SurfaceAdaptor::of_face(&fv));
        }
        // `Bad.Initialize(Fad); Bop.Initialize(Fop);` (`cxx:853-854`).
        let fad = oc
            .fad
            .clone()
            .ok_or_else(|| "PerformOneCorner: IFadArc is not a face".to_string())?;
        let fop = oc
            .fop
            .clone()
            .ok_or_else(|| "PerformOneCorner: IFopArc is not a face".to_string())?;
        oc.hbad = Some(SurfaceAdaptor::of_face_restricted(&fad));
        oc.hbop = Some(SurfaceAdaptor::of_face_restricted(&fop));

        oc.fv = Some(fv);
    }

    // `HGs = new GeomAdaptor_Surface(DStr.Surface(Fd->Surf()).Surface());`
    // (`cxx:867-868`).
    let hgs: Arc<dyn Surface> = {
        let ds = b.my_ds.change_ds();
        ds.surface(oc.fd.surf())
            .surface()
            .cloned()
            .ok_or_else(|| "PerformOneCorner: DStr.Surface(Surf()).Surface() is null".to_string())?
    };
    oc.hgs = Some(Arc::clone(&hgs));

    // `CPopArc` / `FiopArc` / `CPadArc` / `FiadArc` are references into `Fd`
    // (`cxx:859-862`); the port works on copies and writes them back below.
    let mut c_pop = oc.fd.vertex(is_first, oc.i_fop_arc).clone();
    let mut fi_op = oc.fd.interference(oc.i_fop_arc).clone();
    let mut c_pad = oc.fd.vertex(is_first, oc.i_fad_arc).clone();
    let mut fi_ad = oc.fd.interference(oc.i_fad_arc).clone();

    // `double wop = Fd->ChangeInterference(IFadArc).Parameter(isfirst);`
    // (`cxx:866`).
    oc.wop = fi_ad.parameter(is_first);

    if onsame {
        // `cxx:872-898`: `IntersUpdateOnSame` over the extended face at end.
        let save_c_pop = c_pop.clone();
        let c3df = {
            let ds = b.my_ds.change_ds();
            ds.curve(fi_op.line_index()).curve().cloned().ok_or_else(|| {
                "PerformOneCorner: FiopArc.LineIndex() is not in the DS".to_string()
            })?
        };
        let fop = oc
            .fop
            .clone()
            .ok_or_else(|| "PerformOneCorner: IFopArc is not a face".to_string())?;
        let fv = oc
            .fv
            .clone()
            .ok_or_else(|| "OneCorner : face at end not found".to_string())?;
        let arcprol = Edge(
            oc.arcprol
                .clone()
                .ok_or_else(|| "PerformOneCorner: no Arcprol".to_string())?,
        );
        let hbs_surf = oc
            .hbs
            .as_ref()
            .and_then(|a| a.surface())
            .cloned()
            .ok_or_else(|| "PerformOneCorner: Bs has no surface".to_string())?;
        let mut p2dbout = occt_core::gp::GpPnt2d::zero();
        let mut wop = oc.wop;
        let inters = inter_update_on_same(
            &hgs,
            &hbs_surf,
            &c3df,
            &fop,
            &fv,
            &arcprol,
            &vtx,
            is_first,
            10.0 * b.tolapp3d,
            &mut fi_op,
            &mut c_pop,
            &mut p2dbout,
            &mut wop,
        )?;
        oc.p2dbout = p2dbout;
        oc.wop = wop;
        oc.inters = inters;

        // `pced->Initialize(CPadArc.Arc(), Fv); if ((FiadArc.LastParameter() -
        // FiadArc.FirstParameter()) > 10 * tolesp) Update(...)` (`cxx:890-897`).
        let fv = oc.fv.clone().unwrap();
        let arc = c_pad
            .arc()
            .cloned()
            .ok_or_else(|| "PerformOneCorner: CPadArc is not on an arc".to_string())?;
        let (pced, _, _) = curve_on_surface_range(&Edge(arc), &fv)
            .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
        if (fi_ad.last_parameter() - fi_ad.first_parameter()) > 10.0 * b.tolesp {
            update_curve_on_surface(&hbs_surf, &pced, &hgs, &mut fi_ad, &mut c_pad, is_first);
        }

        *oc.fd.vertex_mut(is_first, oc.i_fop_arc) = c_pop.clone();
        *oc.fd.interference_mut(oc.i_fop_arc) = fi_op;
        *oc.fd.vertex_mut(is_first, oc.i_fad_arc) = c_pad;
        *oc.fd.interference_mut(oc.i_fad_arc) = fi_ad;

        // `cxx:899-904`: `saveCPopArc.SetPoint(CPopArc.Point()); CPopArc =
        // saveCPopArc; return;`
        if prepare_on_same {
            let mut save = save_c_pop;
            save.set_point(c_pop.point());
            *oc.fd.vertex_mut(is_first, oc.i_fop_arc) = save;
            return Ok(Some(oc));
        }
    } else {
        // `inters = FindFace(Vtx, CV1, CV2, Fv, Fop); if (!inters) {
        //  PerformIntersectionAtEnd(Index); return; }` (`cxx:909-914`).
        let cv1 = oc.fd.vertex(is_first, 1).clone();
        let cv2 = oc.fd.vertex(is_first, 2).clone();
        let fop = oc
            .fop
            .clone()
            .ok_or_else(|| "OneCorner : face at end not found".to_string())?;
        let Some(fv) = b.find_face(&vtx, &cv1, &cv2, Some(&fop)) else {
            oc.exit = Some(OneCornerExit::IntersectionAtEnd);
            return Ok(Some(oc));
        };
        oc.fv = Some(fv.clone());

        // `Bs.Initialize(Fv);` (`cxx:916`).
        oc.hbs = Some(SurfaceAdaptor::of_face_restricted(&fv));
        let hbs_surf = oc
            .hbs
            .as_ref()
            .and_then(|a| a.surface())
            .cloned()
            .ok_or_else(|| "PerformOneCorner: Bs has no surface".to_string())?;

        // `pced->Initialize(CV1.Arc(), Fv); Update(HBs, pced, HGs,
        //  Fd->ChangeInterferenceOnS1(), CV1, isfirst);` then the S2 pair
        // (`cxx:917-920`).
        for on_s in 1..=2 {
            let arc = oc
                .fd
                .vertex(is_first, on_s)
                .arc()
                .cloned()
                .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
            let (pced, _, _) = curve_on_surface_range(&Edge(arc), &fv)
                .ok_or_else(|| "Failed to get p-curve of edge".to_string())?;
            let mut fi = oc.fd.interference(on_s).clone();
            let mut cp = oc.fd.vertex(is_first, on_s).clone();
            update_curve_on_surface(&hbs_surf, &pced, &hgs, &mut fi, &mut cp, is_first);
            *oc.fd.interference_mut(on_s) = fi;
            *oc.fd.vertex_mut(is_first, on_s) = cp;
        }
    }

    Ok(Some(oc))
}
