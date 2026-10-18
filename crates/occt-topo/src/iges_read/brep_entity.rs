//! `IGESToBRep_BRepEntity` (`IGESToBRep_BRepEntity.cxx`) for the B-rep family
//! reached from `TransferCurveAndSurface` (`IGESToBRep_CurveAndSurface.cxx:184-188`).
//!
//! `TransferBRepEntity` (`cxx:110-136`) dispatches a 510 `IGESSolid_Face`, a
//! 514 `IGESSolid_Shell`, and a 186 `IGESSolid_ManifoldSolid`; the 502 / 504 /
//! 508 list entities fall into the `else` and fail with `IGES_1005` (a null
//! result). The sub-results are cached the way OCCT caches them:
//! `TransferFace` / `TransferShell` / `TransferManifoldSolid` / `TransferLoop`
//! use the single-shape binder (`HasShapeResult` / `SetShapeResult` /
//! `GetShapeResult`), while `TransferVertex` / `TransferEdge` fill and read the
//! 1-based `TransferBRep_ShapeListBinder` of their 502 / 504 entity.
//!
//! Three OCCT calls of `TransferEdge` / `TransferFace` / `TransferShell` have
//! no counterpart call in this module. Each is listed with its OCCT file and
//! line, and with the reason the result is already carried by the port:
//!
//! * `BRep_Builder::UpdateVertex(V, Par, E, 0.)` (`IGESToBRep_BRepEntity.cxx:258-259`,
//!   `:265-266`, `:272-273`) runs after `B.Add(E, V1/V2)`, so the vertex is
//!   found in the edge and the orientation resolves to `FORWARD` / `REVERSED`;
//!   that arm only writes `GC->First(Par)` / `GC->Last(Par)` on every curve
//!   representation (`BRep_Builder.cxx:1220-1312`). Only the `INTERNAL` arm
//!   writes `BRep_PointRepresentation` entries into the vertex, and it cannot
//!   be reached here. The following `B.Range(E, first, last)` writes the same
//!   `First` / `Last` to the same set (`Only3d = false`, `BRep_Builder.cxx:1091-1117`),
//!   so this port calls only `GeometryRegistry::set_ranges_all`
//!   (`rebuild_edge`), and the `Tol = 0.` argument cannot raise a tolerance
//!   either (`BRep_TVertex::UpdateTolerance` keeps the larger value).
//! * `BRepTools::Update(F)` (`cxx:527`) is `UpdateFaceUVPoints`
//!   (`BRepTools.cxx:383-390`, `:488-523`): it caches `myUV1` / `myUV2` of every
//!   `BRep_CurveOnSurface` as `pcurve.D0(rep.First())` / `pcurve.D0(rep.Last())`
//!   (`BRep_CurveOnSurface.cxx:100-114`), because a fresh `EmptyCopy` face has
//!   `Checked() == false`. The port keeps no `UV1` / `UV2` fields: it returns
//!   those two values on demand (`BRep_Tool::UVPoints` equivalent in
//!   `meshing/shape_tool.rs`), so there is nothing to cache.
//! * `F.Orientable(true)` / `S.Orientable(true)` (`cxx:528`, `:600`) set
//!   `TopoDS_TShape`'s orientable flag, which the port leaves at its `true`
//!   default (`TopoDS_TShape`'s constructor ORs in `Bit_Orientable`,
//!   `TopoDS_TShape.hxx:166-170`; `abs.rs` `ShapeFlags::default`).

use occt_core::gp::{GpPnt, GpTrsf2d};
use occt_core::precision::CONFUSION;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex};
use crate::shhealing::{first_vertex, last_vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::is_same;

use super::iges_data_iges_model::EntityData;
use super::iges_to_brep_basic_curve::eps_geom;
use super::iges_to_brep_curve_and_surface::{is_topo_curve, is_topo_surface, Session};
use super::iges_to_brep_topo_curve::{
    boundary_max_tol, edges_of, transfer_pcurve, transfer_topo_curve, WireDataBoundary,
};
use super::iges_to_brep_topo_surface::{transfer_topo_surface, uv_map};

/// `IGESToBRep_BRepEntity::TransferBRepEntity` (`cxx:110-136`).
pub(in crate::iges_read) fn transfer_brep_entity(
    session: &mut Session<'_>,
    index: usize,
) -> Result<Option<TopoShape>, String> {
    let type_number = session.model.entity(index)?.type_number;
    match type_number {
        510 => transfer_face(session, index),
        514 => transfer_shell(session, index),
        186 => transfer_manifold_solid(session, index),
        // `else`: `Message_Msg Msg1005("IGES_1005")` and a null result
        // (`cxx:131-135`). Reachable for 502 / 504 / 508.
        _ => Ok(None),
    }
}

/// `IGESToBRep_BRepEntity::TransferVertex` (`cxx:143-171`).
fn transfer_vertex(
    session: &mut Session<'_>,
    list_index: usize,
    index: usize,
) -> Result<Option<Vertex>, String> {
    // `int nbshapes = NbShapeResult(start); if (nbshapes == 0) { ... fill ... }`
    // (`cxx:148-159`). OCCT fills the whole list on the first call and caches it.
    if session.nb_shape_result(list_index) == 0 {
        let vertices = match &session.model.entity(list_index)?.data {
            EntityData::VertexList { vertices } => vertices.clone(),
            _ => return Err(format!("IGES {list_index} is not a 502 VertexList")),
        };
        let unit = session.model.unit_factor;
        let builder = TopoBuilder::new();
        for xyz in vertices {
            // `point.Scale(gp_Pnt(0,0,0), GetUnitFactor())` (`cxx:152-153`).
            let point = GpPnt::new(xyz[0] * unit, xyz[1] * unit, xyz[2] * unit);
            // pdn 12.03.99 S4135 Constructing vertex with minimal tolerance.
            let v = builder.make_vertex(point, CONFUSION);
            session.add_shape_result(list_index, Some(v.0));
        }
    }

    // `TopoDS_Shape Sh = GetShapeResult(start, index); if (Sh.IsNull()) { ... }`
    // (`cxx:161-169`): the `IGES_1156` warning has no shape effect.
    Ok(session
        .list_shape_result(list_index, index)
        .and_then(Vertex::wrap))
}

/// `IGESToBRep_BRepEntity::TransferEdge` (`cxx:176-310`).
fn transfer_edge(
    session: &mut Session<'_>,
    list_index: usize,
    index: usize,
) -> Result<Option<TopoShape>, String> {
    if session.nb_shape_result(list_index) == 0 {
        let edges = match &session.model.entity(list_index)?.data {
            EntityData::EdgeList { edges } => edges.clone(),
            _ => return Err(format!("IGES {list_index} is not a 504 EdgeList")),
        };
        let builder = TopoBuilder::new();
        for item in &edges {
            // Vertices (`cxx:189-198`).
            let v1 = match item.start_vertex_list {
                Some(list) => transfer_vertex(
                    session,
                    list,
                    item.start_vertex_index.unwrap_or(0).max(0) as usize,
                )?,
                None => None,
            };
            let v2 = match item.end_vertex_list {
                Some(list) => {
                    transfer_vertex(session, list, item.end_vertex_index.unwrap_or(0).max(0) as usize)?
                }
                None => None,
            };

            // Curve (`cxx:202-208`): null, not a topo curve, a 142
            // `IGESGeom_CurveOnSurface`, or a 141 `IGESGeom_Boundary`.
            let curve_ok = match item.curvature {
                Some(c) => {
                    let t = session.model.entity(c)?.type_number;
                    is_topo_curve(t) && t != 142 && t != 141
                }
                None => false,
            };
            if !curve_ok {
                // Message_Msg Msg1306("IGES_1306") (`cxx:205-208`).
                // add null shape to avoid shift of indexing (`cxx:211`).
                session.add_shape_result(list_index, None);
                continue;
            }

            let curve_index = item.curvature.expect("curve checked above");
            let shape = transfer_topo_curve(session, curve_index)?;
            if shape.is_edge() {
                let edge = Edge(shape);
                let e = rebuild_edge(&builder, &edge, &v1, &v2)?;
                session.add_shape_result(list_index, Some(e));
            } else if shape.is_wire() {
                // `Sh.ShapeType() == TopAbs_WIRE` (`cxx:279-286`): "pas traite",
                // `IGES_1325`, then the wire itself is stored.
                session.add_shape_result(list_index, Some(shape));
            }
            // No `else` in OCCT (`cxx:278-287`), so any other shape type adds
            // nothing and shifts the 1-based indexing.
        }
    }

    // `GetShapeResult(start, index)` (`cxx:300-309`): the `IGES_1156` warning
    // when null has no shape effect.
    Ok(session.list_shape_result(list_index, index))
}

/// The edge rewrite of `TransferEdge` (`cxx:217-277`).
///
/// `Vf` / `Vl` are `TopExp::Vertices(edge, Vf, Vl)`. The vertices are added
/// FORWARD / REVERSED and the range is written to every representation
/// (`B.Range(E, first, last)`, `Only3d = false`, `BRep_Builder.cxx:1091-1117`).
fn rebuild_edge(
    builder: &TopoBuilder,
    edge: &Edge,
    v1: &Option<Vertex>,
    v2: &Option<Vertex>,
) -> Result<TopoShape, String> {
    let (Some(v1), Some(v2)) = (v1, v2) else {
        // OCCT dereferences the null vertex (`BRep_Tool::Pnt`, `cxx:243-246`).
        return Err("IGES 504 edge vertex is null: IGESToBRep_BRepEntity.cxx:190-198".into());
    };
    let (Some(vf), Some(vl)) = (first_vertex(edge), last_vertex(edge)) else {
        return Err("IGES 504 edge has no endpoints: IGESToBRep_BRepEntity.cxx:221".into());
    };

    let reg = GeometryRegistry::global();
    let (first, last) = reg.edge_parameters(&edge.0);
    // `BRep_Tool::Curve(edge, loc, first, last)` (`cxx:227`): the port keeps the
    // 3D curve in the edge core and does not model the curve location
    // separately (transfer results carry an identity location here).
    let crv = reg
        .edge_curve(&edge.0)
        .ok_or_else(|| "IGES 504 edge has no 3D curve: IGESToBRep_BRepEntity.cxx:227".to_string())?;
    // In the case of a conic, it is necessary to reverse IGES direction of
    // travel inverse to CASCADE direction of travel (`cxx:228-239`).
    let new_c3d = if crv.is_geom_trimmed() {
        crv.trim_basis().unwrap_or(crv)
    } else {
        crv
    };

    let mut e = builder.make_edge(new_c3d, first, last);

    let p1 = BRepTool::vertex_point(v1);
    let p2 = BRepTool::vertex_point(v2);
    let pf = BRepTool::vertex_point(&vf);
    let pl = BRepTool::vertex_point(&vl);
    let dist1f = p1.distance(&pf);
    let dist2f = p2.distance(&pf);
    let dist1l = p1.distance(&pl);
    let dist2l = p2.distance(&pl);

    if is_same(&v1.0, &v2.0)
        || dist1f + dist2l <= dist1l + dist2f + CONFUSION
    {
        add_edge_vertices(builder, &mut e, v1, v2);
    } else {
        // modification mjm of 13/10/97 : Reverse the edge? (`cxx:268-273`).
        e.0.reverse();
        add_edge_vertices(builder, &mut e, v1, v2);
    }
    // `B.Range(E, first, last)` (`cxx:261`, `:274`).
    reg.set_ranges_all(&e.0, first, last);
    Ok(e.0)
}

fn add_edge_vertices(builder: &TopoBuilder, e: &mut Edge, v1: &Vertex, v2: &Vertex) {
    let mut vf = v1.0.clone();
    vf.set_orientation(Orientation::Forward);
    builder.add(&mut e.0, &vf);
    let mut vl = v2.0.clone();
    vl.set_orientation(Orientation::Reversed);
    builder.add(&mut e.0, &vl);
}

/// `IGESToBRep_BRepEntity::TransferLoop` (`cxx:315-461`).
fn transfer_loop(
    session: &mut Session<'_>,
    loop_index: usize,
    face: &TopoShape,
    trans: &GpTrsf2d,
    u_fact: f64,
) -> Result<Option<TopoShape>, String> {
    if session.has_shape_result(loop_index) {
        return Ok(session.get_shape_result(loop_index));
    }

    let items = match &session.model.entity(loop_index)?.data {
        EntityData::SolidLoop { edges } => edges.clone(),
        _ => return Err(format!("IGES {loop_index} is not a 508 Loop")),
    };

    // `IGESToBRep_IGESBoundary::Init` (`IGESToBRep_IGESBoundary.cxx:41-64`) with
    // `filepreference = 3` (`cxx:326`).
    let precision = eps_geom(session.model.resolution) * session.model.unit_factor;
    let maxtol = boundary_max_tol(session);
    let mut boundary = WireDataBoundary::new(precision, maxtol, 3);
    // `occ::handle<ShapeExtend_WireData> sewd` starts null and is only assigned
    // at `iedge == 1` (`cxx:323`, `:424-427`).
    let mut sewd_started = false;

    for (i0, item) in items.iter().enumerate() {
        let iedge = i0 + 1;
        let itype = item.edge_type.unwrap_or(-1);
        let indexlist = item.list_index.unwrap_or(0).max(0) as usize;
        let orientation = item.orientation;
        let nbparam = item.param_curves.len();

        let Some(edge_index) = item.edge else {
            // `IGES_1365` "Loop : one edge is null" (`cxx:343-347`).
            continue;
        };
        let edge_type = session.model.entity(edge_index)?.type_number;

        // edge (`cxx:351-388`).
        let mut curve3d: Vec<Edge> = Vec::new();
        if itype == 1 && edge_type == 502 {
            let Some(v1) = transfer_vertex(session, edge_index, indexlist)? else {
                continue;
            };
            // `B.MakeEdge(E)` + the V1 FORWARD / V1 REVERSED pair +
            // `B.Degenerated(E, true)` (`cxx:358-367`).
            let p = BRepTool::vertex_point(&v1);
            let mut e = TopoBuilder::new().make_degenerated_edge(&p);
            let mut f = v1.0.clone();
            f.set_orientation(Orientation::Forward);
            TopoBuilder::new().add(&mut e.0, &f);
            let mut r = v1.0.clone();
            r.set_orientation(Orientation::Reversed);
            TopoBuilder::new().add(&mut e.0, &r);
            curve3d.push(e);
        } else if itype == 0 && edge_type == 504 {
            let Some(shape) = transfer_edge(session, edge_index, indexlist)? else {
                // `continue; // skip non-translated edge` (`cxx:374-377`).
                continue;
            };
            curve3d.extend(edges_of(&shape)?);
        } else {
            // `IGES_1365` "Improper type for the edge" (`cxx:381-386`).
            continue;
        }

        if !orientation {
            // `curve3d->Reverse()` (`cxx:389-392`),
            // `ShapeExtend_WireData::Reverse` (`ShapeExtend_WireData.cxx:483-506`).
            curve3d.reverse();
            for e in curve3d.iter_mut() {
                e.0.reverse();
            }
        }

        // processing of 2d curves (`cxx:394-413`). `Curves2d` stays null unless
        // the single-edge / single-param-curve / one-topo-curve test passes.
        let mut curves2d: Vec<Option<usize>> = Vec::new();
        if curve3d.len() == 1 && nbparam == 1 {
            if let Some(&Some(pc)) = item.param_curves.first() {
                if param_curve_count(session, pc)? == 1 {
                    curves2d = item.param_curves.clone();
                }
            }
        }

        let lsewd = boundary.transfer(
            session, face, trans, u_fact, &curve3d, &curves2d, !orientation, iedge,
        )?;
        if iedge == 1 {
            // `sewd = IB->WireData(); // initialization` (`cxx:423-427`).
            sewd_started = true;
        }

        // to keep sharing of edges all geometric representations should be put
        // into the edge from EdgeList (`cxx:428-439`).
        if curve3d.len() == 1 && lsewd.len() == 1 {
            let fromedge = lsewd[0].clone();
            let toedge = curve3d[0].clone();
            if !is_same(&fromedge.0, &toedge.0) {
                let reg = GeometryRegistry::global();
                reg.remove_edge_curve3d(&toedge.0); // `sbe.RemoveCurve3d(toedge)` (`cxx:435`)
                transfer_pcurve(&fromedge, &toedge, face); // `cxx:436`
                // `sewd->Set(toedge, sewd->Index(fromedge))` (`cxx:437`).
                set_wire_edge(&mut boundary.out, &toedge, &fromedge);
            }
        }
    }

    // `if (!sewd.IsNull()) { mywire = sewd->Wire(); }` (`cxx:442-447`).
    let result = if sewd_started {
        Some(TopoBuilder::new().make_wire(&boundary.out).0)
    } else {
        None
    };
    // `SetShapeResult(start, mywire)` (`cxx:448`), even when the wire is null.
    session.set_shape_result(loop_index, result.clone());
    Ok(result)
}

/// `ShapeExtend_WireData::Index` (`ShapeExtend_WireData.cxx:616-626`) then
/// `Set` (`cxx:456-475`): an edge not found at index 0 rewrites the last slot
/// (`myEdges->SetValue((num > 0 ? num : NbEdges()), edge)`). The seam exemption
/// in `Index` cannot trigger here: the accumulated `mysewd` has no seam
/// bookkeeping.
fn set_wire_edge(out: &mut [Edge], to: &Edge, from: &Edge) {
    if out.is_empty() {
        return;
    }
    let pos = out
        .iter()
        .position(|e| is_same(&e.0, &from.0))
        .unwrap_or(out.len() - 1);
    out[pos] = to.clone();
}

/// `IGESToBRep::IGESCurveToSequenceOfIGESCurve` length (`IGESToBRep.cxx:365-388`):
/// composites recurse, a topo curve that is not a 116 point counts as one,
/// anything else counts zero.
fn param_curve_count(session: &Session<'_>, index: usize) -> Result<usize, String> {
    match &session.model.entity(index)?.data {
        EntityData::Composite { curves } => {
            let mut total = 0;
            for child in curves.clone() {
                total += param_curve_count(session, child)?;
            }
            Ok(total)
        }
        _ => {
            let t = session.model.entity(index)?.type_number;
            Ok(usize::from(is_topo_curve(t) && t != 116))
        }
    }
}

/// `IGESToBRep_BRepEntity::TransferFace` (`cxx:464-557`).
fn transfer_face(
    session: &mut Session<'_>,
    index: usize,
) -> Result<Option<TopoShape>, String> {
    if session.has_shape_result(index) {
        return Ok(session.get_shape_result(index));
    }

    let (surface, loops) = match &session.model.entity(index)?.data {
        EntityData::SolidFace { surface, loops, .. } => (*surface, loops.clone()),
        _ => return Err(format!("IGES {index} is not a 510 Face")),
    };

    // `if (surf.IsNull() || !IsTopoSurface(surf) || IsKind(Plane) ||
    // IsKind(BoundedSurface) || IsKind(TrimmedSurface) || IsKind(SingleParent))`
    // (`cxx:480-485`).
    let surf_ok = match surface {
        Some(s) => {
            let ent = session.model.entity(s)?;
            // `IsKind(SingleParent)` is the 402 form 9 class.
            is_topo_surface(session.model, s)
                && !matches!(ent.type_number, 108 | 143 | 144)
                && !matches!(&ent.data, EntityData::SingleParent { .. })
        }
        None => false,
    };
    if !surf_ok {
        // Message_Msg Msg196("XSTEP_196") (`cxx:487-490`).
        session.set_shape_result(index, None);
        return Ok(None);
    }
    let surf = surface.expect("surface checked above");

    // `TopoDS_Shape myshape = TS.ParamSurface(surf, trans, uFact)` (`cxx:498`).
    let myshape = transfer_topo_surface(session, surf)?;
    let Some(src_face) = face_of(&myshape) else {
        // The `myshape.ShapeType() == TopAbs_FACE` gate (`cxx:502`) leaves the
        // result unbound for any other shape (e.g. a multi-face shell).
        return Ok(None);
    };

    let reg = GeometryRegistry::global();
    let Some(surface_geom) = reg.face_surface(&src_face.0) else {
        return Err(format!("IGES {index} face has no surface: IGESToBRep_BRepEntity.cxx:1919"));
    };
    // `F = TopoDS::Face(myshape); F.EmptyCopy();` (`cxx:506-507`): a fresh face
    // over the same surface with no wires, keeping the location and tolerance.
    let mut new_face = TopoBuilder::new().make_face(surface_geom, &[]);
    if !src_face.0.location().is_identity() {
        new_face.0.move_location(&src_face.0.location());
    }
    reg.set_face_tolerance(&new_face.0, reg.face_tolerance(&src_face.0));

    // `if (nbloops == 0) { B.NaturalRestriction(F, true); }` (`cxx:509-512`).
    if loops.is_empty() {
        reg.set_natural_restriction(&new_face.0, true);
    }

    // `gp_Trsf2d trans; double uFact;` filled by `ParamSurface`
    // (`cxx:494-498`); the port reads them back through `uv_map`.
    let (trans, u_fact) = uv_map(session, surf, &src_face.0)?;

    // Loops (`cxx:515-523`).
    for lp in loops.iter().flatten() {
        if let Some(shape) = transfer_loop(session, *lp, &new_face.0, &trans, u_fact)? {
            TopoBuilder::new().add(&mut new_face.0, &shape);
        }
    }

    // `BRepTools::Update(F)` (`cxx:527`) has no port call: it only caches the
    // edge UV points (`BRepTools.cxx:488-523`), which the port evaluates on
    // demand. `F.Orientable(true)` (`cxx:528`) is already the `TopoDS_TShape`
    // default (`TopoDS_TShape.hxx:166-170`).
    session.set_shape_result(index, Some(new_face.0.clone()));
    Ok(Some(new_face.0))
}

/// `TopoDS_Face` selection of `ParamSurface` (`IGESToBRep_TopoSurface.cxx:1812-1844`):
/// a face passes straight through, a shell is accepted only when it holds
/// exactly one face.
fn face_of(shape: &TopoShape) -> Option<Face> {
    if shape.is_face() {
        return Face::wrap(shape.clone());
    }
    if shape.is_shell() {
        let children = shape.tshape.read().ok()?.children.clone();
        let faces: Vec<TopoShape> = children.into_iter().filter(|c| c.is_face()).collect();
        if faces.len() == 1 {
            return Face::wrap(faces.into_iter().next().expect("one face"));
        }
    }
    None
}

/// `IGESToBRep_BRepEntity::TransferShell` (`cxx:562-622`).
fn transfer_shell(
    session: &mut Session<'_>,
    index: usize,
) -> Result<Option<TopoShape>, String> {
    if session.has_shape_result(index) {
        return Ok(session.get_shape_result(index));
    }

    let (faces, orientations) = match &session.model.entity(index)?.data {
        EntityData::SolidShell { faces, orientations } => (faces.clone(), orientations.clone()),
        _ => return Err(format!("IGES {index} is not a 514 Shell")),
    };

    let mut shell = Shell::new();
    let nbfaces = faces.len();
    if nbfaces != 0 {
        let mut closed = true; // : 39
        for (i, face) in faces.iter().enumerate() {
            let orientation = orientations.get(i).copied().unwrap_or(false);
            let shape = match face {
                Some(f) => transfer_face(session, *f)?,
                // `start->Face(iface)` can be null (`cxx:578`); OCCT would
                // dereference it, the port treats it as a missing face.
                None => None,
            };
            let Some(mut shape) = shape else {
                closed = false; // : 39 by abv 15.12.97
                continue;
            };
            if !orientation {
                shape.reverse();
            }
            TopoBuilder::new().add(&mut shell.0, &shape);
        }
        if !closed {
            // Message_Msg Msg1360("IGES_1360") (`cxx:594-597`).
        }
        shell.0.set_closed(closed); // : 39
        // `S.Orientable(true)` (`cxx:600`) is already the `TopoDS_TShape`
        // default (`TopoDS_TShape.hxx:166-170`).
        session.set_shape_result(index, Some(shell.0.clone()));
        Ok(Some(shell.0))
    } else {
        // Message_Msg Msg200("XSTEP_200") (`cxx:604-606`): no result is bound,
        // so the shell stays retryable like in OCCT.
        Ok(None)
    }
}

/// `IGESToBRep_BRepEntity::TransferManifoldSolid` (`cxx:627-698`).
fn transfer_manifold_solid(
    session: &mut Session<'_>,
    index: usize,
) -> Result<Option<TopoShape>, String> {
    if session.has_shape_result(index) {
        return Ok(session.get_shape_result(index));
    }

    let (shell, isoriented, void_shells) = match &session.model.entity(index)?.data {
        EntityData::ManifoldSolid {
            shell,
            oriented,
            void_shells,
            ..
        } => (*shell, *oriented, void_shells.clone()),
        _ => return Err(format!("IGES {index} is not a 186 ManifoldSolid")),
    };

    let mut solid = Solid::new();
    let shape = match shell {
        Some(s) => transfer_shell(session, s)?,
        None => None,
    };
    if let Some(shape) = shape {
        if shape.is_shell() {
            let mut shell = Shell(shape);
            if !isoriented {
                shell.0.reverse();
            }
            TopoBuilder::new().add(&mut solid.0, &shell.0);
        }
        if !void_shells.is_empty() {
            for void in void_shells.iter().flatten() {
                let a_sh = transfer_shell(session, *void)?;
                if let Some(a_sh) = a_sh {
                    if a_sh.is_shell() {
                        let mut shell = Shell(a_sh);
                        if !isoriented {
                            shell.0.reverse();
                        }
                        TopoBuilder::new().add(&mut solid.0, &shell.0);
                    }
                }
                // The null-void else (`cxx:674-680`) adds a null `TopoDS_Shell`,
                // which `TopoDS_Builder::Add` ignores; the port adds nothing.
            }
        }
    }
    // `SetShapeResult(start, S)` runs even when the outer shell is null
    // (`cxx:683`), so an empty solid is bound.
    session.set_shape_result(index, Some(solid.0.clone()));
    Ok(Some(solid.0))
}
