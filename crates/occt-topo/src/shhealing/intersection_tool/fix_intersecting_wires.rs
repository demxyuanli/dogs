//! `ShapeFix_IntersectionTool::FixIntersectingWires`
//! (`ShapeFix_IntersectionTool.cxx:1834-2530`).
//!
//! The *face*-level twin of `FixSelfIntersectWire` (`intersection_tool.rs`):
//! it pairs up every two boundary wires of one face, intersects every
//! non-adjacent edge cross pair, and splits / cuts / welds on a hit. The
//! `ShapeBuild_ReShape` context is `ShapeFix_Face::myContext`; the tool records
//! edge, wire and face substitutions there and returns the rebuilt face through
//! the caller's `myFace` out-parameter (the tail `face = newface`).
//!
//! Driven by `ShapeFix_Face::FixIntersectingWires`
//! (`ShapeFix_Face.cxx:2821-2825`), which builds
//! `ShapeFix_IntersectionTool(Context(), Precision(), MaxTolerance())` and
//! calls this; the port keeps the free function here and gains the
//! `ShapeFixFace::fix_intersecting_wires` wrapper.
//!
//! A modified wire is rebuilt as a new FORWARD wire from its current edge list
//! (`ShapeExtend_WireData::Wire`, `rebuild_forward_wire`) and recorded with
//! `ctx.replace(old_wire, new_wire)` (`cxx:2489-2498`), then the final
//! `myContext->Replace(face, newface)` is recorded as before.

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpPnt;
use occt_core::intres2d::{IntRes2dDomain, IntRes2dPosition};
use occt_core::precision::CONFUSION;

use occt_geom2d::geom2d_int::Geom2dIntGInter;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_oriented;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Face, TopoShape, Vertex, Wire};
use crate::shape_fix_compose_shell::ReShape;
use crate::shhealing::wire_fix::{
    copy_replace_vertices_with, first_vertex, last_vertex, wire_set_edge_composed,
};
use crate::shhealing::wire_self_inter::point_on_edge;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of_wire, is_same, vertex_position, vertices_of};

use super::{
    create_boxes_2d, cut_edge, find_vert_and_split_edge, select_int_pnt, skey, split_edge1,
    split_edge2, union_vertexes, update_vertex_tolerance, wire_edge_at, wire_len, Boxes, TOLINT,
};

/// `ShapeExtend_WireData::Wire()` (`ShapeExtend_WireData.cxx:651-675`): a new
/// FORWARD wire built from the current composed edge list (`edges_of_wire`
/// composes the wire orientation the way `TopoDS_Iterator` with CumOri does).
fn rebuild_forward_wire(wire: &Wire) -> Wire {
    let edges = edges_of_wire(wire);
    let mut rebuilt = TopoBuilder::new().make_wire(&edges);
    rebuilt.0.set_orientation(Orientation::Forward);
    rebuilt
}

/// `ShapeFix_IntersectionTool::FixIntersectingWires(TopoDS_Face& face)`
/// (`cxx:1834-2530`): returns `isDone`. On success `*face` becomes the rebuilt
/// face (`face = newface`, `cxx:2527`) and `ctx` records
/// `Replace(old_face, new_face)`.
pub fn fix_intersecting_wires(face: &mut Face, ctx: &mut dyn ReShape) -> bool {
    // `cxx:1836-1840`: `myContext.IsNull() || face.IsNull()` - the port always
    // hands in both.

    // `cxx:1842-1858`: `TopoDS_Iterator(SF, false)` splits the direct children
    // into the oriented wires and everything else.
    let ori = face.0.orientation();
    let children: Vec<TopoShape> = {
        let ts = face.0.tshape.read().expect("poisoned TShape lock");
        ts.children.clone()
    };
    let mut seq_wir: Vec<Wire> = Vec::new();
    let mut seq_nm: Vec<TopoShape> = Vec::new();
    for child in children {
        if child.shape_type() != ShapeType::Wire
            || (child.orientation() != Orientation::Forward
                && child.orientation() != Orientation::Reversed)
        {
            seq_nm.push(child);
            continue;
        }
        seq_wir.push(Wire(child));
    }
    // `cxx:1859-1862`: `gka 06.09.04`.
    if seq_wir.len() < 2 {
        return false;
    }

    // `cxx:1864-1869`: the face's widest vertex tolerance.
    let mut max_tol_vert = 0.0f64;
    for v in vertices_of(&face.0) {
        max_tol_vert = max_tol_vert.max(BRepTool::vertex_tolerance(&v));
    }

    // `cxx:1870-1872`: `isDone = false`, `ShapeAnalysis_Surface sas`.
    let mut is_done = false;
    let surf = match BRepTool::face_surface(face) {
        Some(s) => s,
        None => return false,
    };

    // `cxx:1874-1888`: precompute the 2D edge boxes and the total box per wire.
    let nb_wir = seq_wir.len();
    let mut seq_boxes: Vec<Boxes> = Vec::with_capacity(nb_wir);
    let mut seq_total: Vec<BndBox2d> = Vec::with_capacity(nb_wir);
    for w in &seq_wir {
        let mut boxes: Boxes = Boxes::new();
        let total = create_boxes_2d(w, face, &mut boxes);
        seq_boxes.push(boxes);
        seq_total.push(total);
    }

    // `cxx:1889-2507`: `for (n1 = 1; n1 <= NbWires-1; n1++)`.
    for n1 in 1..=nb_wir - 1 {
        let mut wire1 = seq_wir[n1 - 1].clone();
        // `boxes1` is the `aSeqWirEdgeBoxes.ChangeValue(n1)` reference, a local
        // copy here written back after every n2 step.
        let mut boxes1 = seq_boxes[n1 - 1].clone();
        // `aBox1` is a plain copy fetched once, never refreshed by the update
        // tail (`cxx:1893`) - keep it stale exactly like the cxx.
        let a_box1 = seq_total[n1 - 1];

        for n2 in n1 + 1..=nb_wir {
            let mut wire2 = seq_wir[n2 - 1].clone();
            let mut boxes2 = seq_boxes[n2 - 1].clone();
            let a_box2 = seq_total[n2 - 1];
            // `cxx:1903-1906`.
            if !a_box1.is_void() && !a_box2.is_void() && a_box1.is_out_box(&a_box2) {
                continue;
            }

            // `cxx:1907-1909`.
            let mut nb_modif = 0i64;
            let mut nb_replaced = 0i64;
            let mut has_modif_wire = false;
            // `cxx:1910`: `for (num1 = 1; num1 <= sewd1->NbEdges() && NbModif < 30; num1++)`.
            let mut num1 = 1i64;
            while num1 <= wire_len(&wire1) && nb_modif < 30 {
                let edge1 = match wire_edge_at(&wire1, num1) {
                    Some(e) => e,
                    None => break,
                };
                // `cxx:1915`: `for (num2 = 1; num2 <= sewd2->NbEdges() && NbModif < 30; num2++)`.
                let mut num2 = 1i64;
                loop {
                    if !(num2 <= wire_len(&wire2) && nb_modif < 30) {
                        break;
                    }
                    // 0 = `num2++` (continue / fall through), 1 = `num2--; continue`
                    // (net stay), 2 = `break` out of the num2 loop. Each
                    // `break 'body` carries the C++ transfer it mirrors.
                    let mut ctl = 0u8;
                    'body: {
                        let edge2 = match wire_edge_at(&wire2, num2) {
                            Some(e) => e,
                            None => break 'body,
                        };
                        // `cxx:1919-1921`.
                        if is_same(&edge1.0, &edge2.0) {
                            break 'body;
                        }
                        // `cxx:1922-1925`.
                        if BRepTool::is_degenerated(&edge1) || BRepTool::is_degenerated(&edge2) {
                            break 'body;
                        }
                        // `cxx:1926-1929`.
                        let (Some(b1_box), Some(b2_box)) = (
                            boxes1.get(&skey(&edge1)).copied(),
                            boxes2.get(&skey(&edge2)).copied(),
                        ) else {
                            break 'body;
                        };
                        // `cxx:1931`.
                        if b1_box.is_out_box(&b2_box) {
                            break 'body;
                        }

                        // `cxx:1937-1945`: `sae.PCurve` failure is `continue`.
                        let Some((crv1, a1, b1)) = curve_on_surface_oriented(&edge1, face, false)
                        else {
                            break 'body;
                        };
                        let Some((crv2, a2, b2)) = curve_on_surface_oriented(&edge2, face, false)
                        else {
                            break 'body;
                        };

                        // `cxx:1946-1954`.
                        let d1 = IntRes2dDomain::bounded(
                            &crv1.d0(a1),
                            a1,
                            TOLINT,
                            &crv1.d0(b1),
                            b1,
                            TOLINT,
                        );
                        let d2 = IntRes2dDomain::bounded(
                            &crv2.d0(a2),
                            a2,
                            TOLINT,
                            &crv2.d0(b2),
                            b2,
                            TOLINT,
                        );
                        let mut inter = Geom2dIntGInter::new();
                        inter.perform(crv1.as_ref(), &d1, crv2.as_ref(), &d2, TOLINT, TOLINT);
                        if !inter.is_done() {
                            break 'body;
                        }

                        // `cxx:1956-2042`: intersection is a point.
                        if inter.nb_points() > 0 && inter.nb_points() < 3 {
                            let ip = select_int_pnt(&inter);
                            let tr1 = *ip.transition_of_first();
                            let tr2 = *ip.transition_of_second();

                            // `cxx:1962-1992`.
                            if tr1.position_on_curve() == IntRes2dPosition::Middle
                                && tr2.position_on_curve() == IntRes2dPosition::Middle
                            {
                                let param1 = ip.param_on_first();
                                let param2 = ip.param_on_second();
                                let pi1 =
                                    point_on_edge(&edge1, surf.as_ref(), crv1.as_ref(), param1);
                                let pi2 =
                                    point_on_edge(&edge2, surf.as_ref(), crv2.as_ref(), param2);
                                let p0 = GpPnt::new(
                                    0.5 * (pi1.x() + pi2.x()),
                                    0.5 * (pi1.y() + pi2.y()),
                                    0.5 * (pi1.z() + pi2.z()),
                                );
                                let tol_v = ((pi1.distance(&pi2) / 2.0) * 1.00001)
                                    .max(CONFUSION);
                                let v = TopoBuilder::new().make_vertex(p0, tol_v);
                                max_tol_vert = max_tol_vert.max(tol_v);
                                let is_split_edge2 = split_edge1(
                                    &mut wire2,
                                    ctx,
                                    face,
                                    num2,
                                    param2,
                                    &v,
                                    tol_v,
                                    &mut boxes2,
                                );
                                if is_split_edge2 {
                                    nb_modif += 1;
                                    num2 -= 1;
                                }
                                if split_edge1(
                                    &mut wire1,
                                    ctx,
                                    face,
                                    num1,
                                    param1,
                                    &v,
                                    tol_v,
                                    &mut boxes1,
                                ) {
                                    nb_modif += 1;
                                    num1 -= 1;
                                    ctl = 2;
                                    break 'body;
                                }
                                if is_split_edge2 {
                                    // `continue` with the `num2--` already applied.
                                    ctl = 1;
                                    break 'body;
                                }
                            }
                            // `cxx:1993-2011`.
                            if tr1.position_on_curve() == IntRes2dPosition::Middle
                                && tr2.position_on_curve() != IntRes2dPosition::Middle
                            {
                                let param1 = ip.param_on_first();
                                if find_vert_and_split_edge(
                                    param1,
                                    &edge1,
                                    &edge2,
                                    crv1.as_ref(),
                                    &mut max_tol_vert,
                                    &mut num1,
                                    &mut wire1,
                                    face,
                                    &mut boxes1,
                                    true,
                                    ctx,
                                ) {
                                    nb_modif += 1;
                                    ctl = 2;
                                    break 'body; // `break`
                                }
                            }
                            // `cxx:2012-2033`.
                            if tr1.position_on_curve() != IntRes2dPosition::Middle
                                && tr2.position_on_curve() == IntRes2dPosition::Middle
                            {
                                let param2 = ip.param_on_second();
                                if find_vert_and_split_edge(
                                    param2,
                                    &edge2,
                                    &edge1,
                                    crv2.as_ref(),
                                    &mut max_tol_vert,
                                    &mut num2,
                                    &mut wire2,
                                    face,
                                    &mut boxes2,
                                    true,
                                    ctx,
                                ) {
                                    nb_modif += 1;
                                    break 'body; // `continue` (ctl stays 0)
                                }
                            }
                            // `cxx:2034-2041`.
                            if tr1.position_on_curve() != IntRes2dPosition::Middle
                                && tr2.position_on_curve() != IntRes2dPosition::Middle
                            {
                                if union_vertexes(
                                    &mut wire2,
                                    ctx,
                                    &edge1,
                                    &edge2,
                                    num2,
                                    &mut boxes2,
                                    &b2_box,
                                ) {
                                    nb_replaced += 1;
                                }
                            }
                        }

                        // `cxx:2045`, `gka 06.09.04`.
                        has_modif_wire = has_modif_wire || nb_modif != 0 || nb_replaced != 0;

                        // `cxx:2047-2482`: intersection is a segment.
                        if inter.nb_segments() == 1 {
                            let seg = inter.segment(1);
                            if seg.has_first_point() && seg.has_last_point() {
                                let ipf = seg.first_point().clone();
                                let ipl = seg.last_point().clone();
                                let p11 = ipf.param_on_first();
                                let p21 = ipf.param_on_second();
                                let p12 = ipl.param_on_first();
                                let p22 = ipl.param_on_second();
                                let pnt11 =
                                    point_on_edge(&edge1, surf.as_ref(), crv1.as_ref(), p11);
                                let pnt12 =
                                    point_on_edge(&edge1, surf.as_ref(), crv1.as_ref(), p12);
                                let pnt21 =
                                    point_on_edge(&edge2, surf.as_ref(), crv2.as_ref(), p21);
                                let pnt22 =
                                    point_on_edge(&edge2, surf.as_ref(), crv2.as_ref(), p22);

                                // `cxx:2067-2113`: analysis for edge1.
                                let (Some(sv1), Some(sv2)) =
                                    (first_vertex(&edge1), last_vertex(&edge1))
                                else {
                                    break 'body;
                                };
                                let pv1 = vertex_position(&sv1);
                                let pv2 = vertex_position(&sv2);
                                let mut is_modified1 = false;
                                let mut new_v: Option<Vertex> = None;
                                let mut newtol = 0.0f64;

                                let dist1 = pnt11.distance(&pv1);
                                let dist2 = pnt12.distance(&pv1);
                                let maxdist = dist1.max(dist2);
                                let pdist = if edge1.0.orientation() == Orientation::Reversed {
                                    (b1 - p11).abs().max((b1 - p12).abs())
                                } else {
                                    (a1 - p11).abs().max((a1 - p12).abs())
                                };
                                if maxdist < max_tol_vert || pdist < (b1 - a1).abs() * 0.01 {
                                    newtol = maxdist;
                                    new_v = Some(sv1.clone());
                                    is_modified1 = true;
                                }
                                let dist1 = pnt11.distance(&pv2);
                                let dist2 = pnt12.distance(&pv2);
                                let maxdist = dist1.max(dist2);
                                let pdist = if edge1.0.orientation() == Orientation::Reversed {
                                    (a1 - p11).abs().max((a1 - p12).abs())
                                } else {
                                    (b1 - p11).abs().max((b1 - p12).abs())
                                };
                                if maxdist < max_tol_vert || pdist < (b1 - a1).abs() * 0.01 {
                                    if (is_modified1 && maxdist < newtol) || !is_modified1 {
                                        newtol = maxdist;
                                        new_v = Some(sv2.clone());
                                        is_modified1 = true;
                                    }
                                }
                                if is_modified1 {
                                    // `cxx:2114-2142`: cut edge1 and widen NewV.
                                    let dista = (a1 - p11).abs() + (a1 - p12).abs();
                                    let distb = (b1 - p11).abs() + (b1 - p12).abs();
                                    let pend = if dista > distb { a1 } else { b1 };
                                    let cut =
                                        if (pend - p11).abs() > (pend - p12).abs() { p12 } else { p11 };
                                    let mut is_cut_line = false;
                                    if !cut_edge(&edge1, pend, cut, face, &mut is_cut_line) {
                                        // `IsModified1 = false; continue;` - the
                                        // local is dead past the `continue`.
                                        break 'body;
                                    }
                                    if let Some(nv) = &new_v {
                                        if newtol > BRepTool::vertex_tolerance(nv) {
                                            update_vertex_tolerance(nv, newtol * 1.00001);
                                        }
                                    }
                                }

                                // `cxx:2144-2218`: analysis for edge2.
                                let (Some(sv12), Some(sv22)) =
                                    (first_vertex(&edge2), last_vertex(&edge2))
                                else {
                                    break 'body;
                                };
                                let pv12 = vertex_position(&sv12);
                                let pv22 = vertex_position(&sv22);
                                let mut is_modified2 = false;

                                let dist1 = pnt21.distance(&pv12);
                                let dist2 = pnt22.distance(&pv12);
                                let maxdist = dist1.max(dist2);
                                let pdist = if edge2.0.orientation() == Orientation::Reversed {
                                    (b2 - p21).abs().max((b2 - p22).abs())
                                } else {
                                    (a2 - p21).abs().max((a2 - p22).abs())
                                };
                                if maxdist < max_tol_vert || pdist < (b2 - a2).abs() * 0.01 {
                                    newtol = maxdist;
                                    new_v = Some(sv12.clone());
                                    is_modified2 = true;
                                }
                                let dist1 = pnt21.distance(&pv22);
                                let dist2 = pnt22.distance(&pv22);
                                let maxdist = dist1.max(dist2);
                                let pdist = if edge2.0.orientation() == Orientation::Reversed {
                                    (a2 - p21).abs().max((a2 - p22).abs())
                                } else {
                                    (b2 - p21).abs().max((b2 - p22).abs())
                                };
                                if maxdist < max_tol_vert || pdist < (b2 - a2).abs() * 0.01 {
                                    if (is_modified2 && maxdist < newtol) || !is_modified2 {
                                        newtol = maxdist;
                                        new_v = Some(sv22.clone());
                                        is_modified2 = true;
                                    }
                                }
                                if is_modified2 {
                                    // `cxx:2185-2218`: cut edge2 and widen NewV.
                                    let dista = (a2 - p21).abs() + (a2 - p22).abs();
                                    let distb = (b2 - p21).abs() + (b2 - p22).abs();
                                    let pend = if dista > distb { a2 } else { b2 };
                                    let cut =
                                        if (pend - p21).abs() > (pend - p22).abs() { p22 } else { p21 };
                                    let mut is_cut_line = false;
                                    if !cut_edge(&edge2, pend, cut, face, &mut is_cut_line) {
                                        // `IsModified2 = false; continue;` - the
                                        // local is dead past the `continue`.
                                        break 'body;
                                    }
                                    if let Some(nv) = &new_v {
                                        if newtol > BRepTool::vertex_tolerance(nv) {
                                            update_vertex_tolerance(nv, newtol * 1.00001);
                                        }
                                    }
                                }

                                // `cxx:2219-2228`.
                                if is_modified1 || is_modified2 {
                                    num2 -= 1;
                                    has_modif_wire = true;
                                    ctl = 1;
                                    break 'body; // `num2--; continue;`
                                }

                                if (p12 - p11).abs() > (b1 - a1).abs() / 2.0
                                    || (p22 - p21).abs() > (b2 - a2).abs() / 2.0
                                {
                                    // `cxx:2232-2449`: the segment is big, split each
                                    // intersecting edge on three edges and drop the
                                    // middle (segment) pair.
                                    let p01 = GpPnt::new(
                                        0.5 * (pnt11.x() + pnt21.x()),
                                        0.5 * (pnt11.y() + pnt21.y()),
                                        0.5 * (pnt11.z() + pnt21.z()),
                                    );
                                    let p02 = GpPnt::new(
                                        0.5 * (pnt12.x() + pnt22.x()),
                                        0.5 * (pnt12.y() + pnt22.y()),
                                        0.5 * (pnt12.z() + pnt22.z()),
                                    );
                                    let mut tol_v1 =
                                        pnt11.distance(&p01).max(pnt21.distance(&p01));
                                    tol_v1 = tol_v1.max(CONFUSION) * 1.00001;
                                    let mut tol_v2 =
                                        pnt12.distance(&p02).max(pnt22.distance(&p02));
                                    tol_v2 = tol_v2.max(CONFUSION) * 1.00001;
                                    if tol_v1 > max_tol_vert || tol_v2 > max_tol_vert {
                                        break 'body; // `continue`
                                    }

                                    has_modif_wire = true; // `cxx:2252`
                                    let mut new_v1: Option<Vertex> = None;
                                    let mut new_v2: Option<Vertex> = None;
                                    let mut akey1 = 0i32;
                                    let mut akey2 = 0i32;
                                    // `cxx:2258-2293`: P01 / P02 against V1 / V2.
                                    if p01.distance(&pv1)
                                        < tol_v1.max(BRepTool::vertex_tolerance(&sv1))
                                    {
                                        new_v1 = Some(sv1.clone());
                                        if tol_v1 > BRepTool::vertex_tolerance(&sv1) {
                                            update_vertex_tolerance(&sv1, tol_v1);
                                        }
                                        akey1 += 1;
                                    }
                                    if p01.distance(&pv2)
                                        < tol_v1.max(BRepTool::vertex_tolerance(&sv2))
                                    {
                                        new_v1 = Some(sv2.clone());
                                        if tol_v1 > BRepTool::vertex_tolerance(&sv2) {
                                            update_vertex_tolerance(&sv2, tol_v1);
                                        }
                                        akey1 += 1;
                                    }
                                    if p02.distance(&pv1)
                                        < tol_v2.max(BRepTool::vertex_tolerance(&sv1))
                                    {
                                        new_v2 = Some(sv1.clone());
                                        if tol_v2 > BRepTool::vertex_tolerance(&sv1) {
                                            update_vertex_tolerance(&sv1, tol_v2);
                                        }
                                        akey2 += 1;
                                    }
                                    if p02.distance(&pv2)
                                        < tol_v2.max(BRepTool::vertex_tolerance(&sv2))
                                    {
                                        new_v2 = Some(sv2.clone());
                                        if tol_v2 > BRepTool::vertex_tolerance(&sv2) {
                                            update_vertex_tolerance(&sv2, tol_v2);
                                        }
                                        akey2 += 1;
                                    }
                                    // `cxx:2294-2297`.
                                    if akey1 > 1 || akey2 > 1 {
                                        break 'body; // `continue`
                                    }
                                    // `cxx:2299-2307`: prepare vertices.
                                    if akey1 == 0 {
                                        new_v1 =
                                            Some(TopoBuilder::new().make_vertex(p01, tol_v1));
                                    }
                                    if akey2 == 0 {
                                        new_v2 =
                                            Some(TopoBuilder::new().make_vertex(p02, tol_v2));
                                    }
                                    let (Some(nv1), Some(nv2)) = (new_v1, new_v2) else {
                                        break 'body;
                                    };

                                    // `cxx:2309-2345`: split edge1.
                                    let mut numseg1 = num1;
                                    if akey1 == 0 && akey2 > 0 {
                                        if split_edge1(
                                            &mut wire1,
                                            ctx,
                                            face,
                                            num1,
                                            p11,
                                            &nv1,
                                            tol_v1,
                                            &mut boxes1,
                                        ) {
                                            nb_modif += 1;
                                            numseg1 = num1 + 1;
                                        }
                                    }
                                    if akey1 > 0 && akey2 == 0 {
                                        if split_edge1(
                                            &mut wire1,
                                            ctx,
                                            face,
                                            num1,
                                            p12,
                                            &nv2,
                                            tol_v2,
                                            &mut boxes1,
                                        ) {
                                            nb_modif += 1;
                                            numseg1 = num1;
                                        }
                                    }
                                    if akey1 == 0 && akey2 == 0 {
                                        // `cxx:2325-2344`: split by p11 then by p12,
                                        // on the next edge when p12 is external.
                                        let mut num1split2 = num1;
                                        if split_edge1(
                                            &mut wire1,
                                            ctx,
                                            face,
                                            num1,
                                            p11,
                                            &nv1,
                                            tol_v1,
                                            &mut boxes1,
                                        ) {
                                            nb_modif += 1;
                                            if let Some(tmp_e) = wire_edge_at(&wire1, num1) {
                                                if let Some((_, a, b)) =
                                                    curve_on_surface_oriented(&tmp_e, face, false)
                                                {
                                                    if (a - p12) * (b - p12) > 0.0 {
                                                        // p12 external for [a, b].
                                                        num1split2 += 1;
                                                    }
                                                }
                                            }
                                        }
                                        if split_edge1(
                                            &mut wire1,
                                            ctx,
                                            face,
                                            num1split2,
                                            p12,
                                            &nv2,
                                            tol_v2,
                                            &mut boxes1,
                                        ) {
                                            nb_modif += 1;
                                            numseg1 = num1 + 1;
                                        }
                                    }
                                    // `cxx:2346`: `SegE = sewd1->Edge(numseg1)`.
                                    let Some(mut seg_e) = wire_edge_at(&wire1, numseg1) else {
                                        break 'body;
                                    };

                                    // `cxx:2351-2399`: replace edge2 vertices when needed.
                                    let mut v12 = sv12.clone();
                                    let mut v22 = sv22.clone();
                                    let mut cur_edge2 = edge2.clone();
                                    let mut akey1 = 0i32;
                                    let mut akey2 = 0i32;
                                    if p01.distance(&pv12) < tol_v1 {
                                        tol_v1 += p01.distance(&pv12);
                                        update_vertex_tolerance(&nv1, tol_v1);
                                        v12 = nv1.clone();
                                        let new_e = copy_replace_vertices_with(
                                            &cur_edge2,
                                            Some(&nv1),
                                            Some(&v22),
                                        );
                                        ctx.replace(&cur_edge2.0, &new_e.0);
                                        wire_set_edge_composed(
                                            &mut wire2,
                                            num2 as usize,
                                            &new_e,
                                        );
                                        boxes2.insert(skey(&new_e), b2_box);
                                        cur_edge2 = new_e;
                                        akey1 = 1;
                                    }
                                    if p01.distance(&pv22) < tol_v1 {
                                        tol_v1 += p01.distance(&pv22);
                                        update_vertex_tolerance(&nv1, tol_v1);
                                        v22 = nv1.clone();
                                        let new_e = copy_replace_vertices_with(
                                            &cur_edge2,
                                            Some(&v12),
                                            Some(&nv1),
                                        );
                                        ctx.replace(&cur_edge2.0, &new_e.0);
                                        wire_set_edge_composed(
                                            &mut wire2,
                                            num2 as usize,
                                            &new_e,
                                        );
                                        boxes2.insert(skey(&new_e), b2_box);
                                        cur_edge2 = new_e;
                                        akey1 = 2;
                                    }
                                    if p02.distance(&pv12) < tol_v2 {
                                        tol_v2 += p02.distance(&pv12);
                                        update_vertex_tolerance(&nv2, tol_v2);
                                        v12 = nv2.clone();
                                        let new_e = copy_replace_vertices_with(
                                            &cur_edge2,
                                            Some(&nv2),
                                            Some(&v22),
                                        );
                                        ctx.replace(&cur_edge2.0, &new_e.0);
                                        wire_set_edge_composed(
                                            &mut wire2,
                                            num2 as usize,
                                            &new_e,
                                        );
                                        boxes2.insert(skey(&new_e), b2_box);
                                        cur_edge2 = new_e;
                                        akey2 = 1;
                                    }
                                    if p02.distance(&pv22) < tol_v2 {
                                        tol_v2 += p02.distance(&pv22);
                                        update_vertex_tolerance(&nv2, tol_v2);
                                        v22 = nv2.clone();
                                        let new_e = copy_replace_vertices_with(
                                            &cur_edge2,
                                            Some(&v12),
                                            Some(&nv2),
                                        );
                                        ctx.replace(&cur_edge2.0, &new_e.0);
                                        wire_set_edge_composed(
                                            &mut wire2,
                                            num2 as usize,
                                            &new_e,
                                        );
                                        boxes2.insert(skey(&new_e), b2_box);
                                        cur_edge2 = new_e;
                                        akey2 = 2;
                                    }

                                    // `cxx:2400-2438`: split edge2.
                                    let mut numseg2 = num2;
                                    if akey1 == 0 && akey2 > 0 {
                                        if split_edge1(
                                            &mut wire2,
                                            ctx,
                                            face,
                                            num2,
                                            p21,
                                            &nv1,
                                            tol_v1,
                                            &mut boxes2,
                                        ) {
                                            nb_modif += 1;
                                            numseg2 = num2 + 1;
                                        }
                                    }
                                    if akey1 > 0 && akey2 == 0 {
                                        if split_edge1(
                                            &mut wire2,
                                            ctx,
                                            face,
                                            num2,
                                            p22,
                                            &nv2,
                                            tol_v2,
                                            &mut boxes2,
                                        ) {
                                            nb_modif += 1;
                                            numseg2 = num2;
                                        }
                                    }
                                    if akey1 == 0 && akey2 == 0 {
                                        let mut num2split2 = num2;
                                        if split_edge1(
                                            &mut wire2,
                                            ctx,
                                            face,
                                            num2,
                                            p21,
                                            &nv1,
                                            tol_v1,
                                            &mut boxes2,
                                        ) {
                                            nb_modif += 1;
                                            numseg2 = num2 + 1;
                                            if let Some(tmp_e) = wire_edge_at(&wire2, num2) {
                                                if let Some((_, a, b)) =
                                                    curve_on_surface_oriented(&tmp_e, face, false)
                                                {
                                                    if (a - p22) * (b - p22) > 0.0 {
                                                        num2split2 += 1;
                                                    }
                                                }
                                            }
                                        }
                                        if split_edge1(
                                            &mut wire2,
                                            ctx,
                                            face,
                                            num2split2,
                                            p22,
                                            &nv2,
                                            tol_v2,
                                            &mut boxes2,
                                        ) {
                                            nb_modif += 1;
                                            numseg2 = num2 + 1;
                                        }
                                    }

                                    // `cxx:2440-2449`.
                                    let Some(tmp_e) = wire_edge_at(&wire2, numseg2) else {
                                        break 'body;
                                    };
                                    if let Some(bx) = boxes1.get(&skey(&seg_e)).copied() {
                                        boxes2.insert(skey(&tmp_e), bx);
                                    }
                                    let same_first = match (
                                        first_vertex(&seg_e),
                                        first_vertex(&tmp_e),
                                    ) {
                                        (Some(a), Some(b)) => is_same(&a.0, &b.0),
                                        _ => false,
                                    };
                                    if !same_first {
                                        seg_e.0.reverse();
                                    }
                                    ctx.replace(&tmp_e.0, &seg_e.0);
                                    wire_set_edge_composed(&mut wire2, numseg2 as usize, &seg_e);
                                    num1 -= 1;
                                    ctl = 2;
                                    break 'body;
                                } else {
                                    // `cxx:2451-2479`: split each intersecting edge on
                                    // two edges.
                                    let p0 = GpPnt::new(
                                        0.5 * (pnt11.x() + pnt12.x()),
                                        0.5 * (pnt11.y() + pnt12.y()),
                                        0.5 * (pnt11.z() + pnt12.z()),
                                    );
                                    let param1 = 0.5 * (p11 + p12);
                                    let param2 = 0.5 * (p21 + p22);
                                    let pnt10 =
                                        point_on_edge(&edge1, surf.as_ref(), crv1.as_ref(), param1);
                                    let pnt20 =
                                        point_on_edge(&edge2, surf.as_ref(), crv2.as_ref(), param2);
                                    let dist1 = pnt11.distance(&p0).max(pnt12.distance(&pnt10));
                                    let dist2 = pnt21.distance(&p0).max(pnt22.distance(&pnt10));
                                    let mut tol_v = dist1.max(dist2);
                                    tol_v = tol_v.max(pnt10.distance(&pnt20)) * 1.00001;
                                    let nv = TopoBuilder::new().make_vertex(pnt10, tol_v);
                                    max_tol_vert = max_tol_vert.max(tol_v);
                                    has_modif_wire = true;
                                    let mut stay = false;
                                    if split_edge2(
                                        &mut wire2,
                                        ctx,
                                        face,
                                        num2,
                                        p21,
                                        p22,
                                        &nv,
                                        tol_v,
                                        &mut boxes2,
                                    ) {
                                        nb_modif += 1;
                                        stay = true;
                                    }
                                    if split_edge2(
                                        &mut wire1,
                                        ctx,
                                        face,
                                        num1,
                                        p11,
                                        p12,
                                        &nv,
                                        tol_v,
                                        &mut boxes1,
                                    ) {
                                        nb_modif += 1;
                                        num1 -= 1;
                                        ctl = 2;
                                        break 'body;
                                    }
                                    if stay {
                                        ctl = 1;
                                        break 'body; // `num2--; continue;`
                                    }
                                }
                            }
                        }
                    }

                    if ctl == 2 {
                        break;
                    }
                    num2 += if ctl == 1 { 0 } else { 1 };
                }
                num1 += 1;
            }

            // `cxx:2483-2505`: persist the modified wires and rethink their boxes.
            if has_modif_wire {
                is_done = true;
                // `cxx:2489-2491`: `SeqWir.SetValue(n1, sewd1->Wire());
                // myContext->Replace(wire1, sewd1->Wire()); wire1 = sewd1->Wire();`
                let new_wire1 = rebuild_forward_wire(&wire1);
                ctx.replace(&wire1.0, &new_wire1.0);
                seq_wir[n1 - 1] = new_wire1.clone();
                wire1 = new_wire1;
                boxes1.clear();
                let a_new_box1 = create_boxes_2d(&wire1, face, &mut boxes1);
                seq_total[n1 - 1] = a_new_box1;
                // `cxx:2496-2498`: the same for wire2.
                let new_wire2 = rebuild_forward_wire(&wire2);
                ctx.replace(&wire2.0, &new_wire2.0);
                seq_wir[n2 - 1] = new_wire2.clone();
                wire2 = new_wire2;
                boxes2.clear();
                let a_new_box2 = create_boxes_2d(&wire2, face, &mut boxes2);
                seq_total[n2 - 1] = a_new_box2;
            }
            seq_boxes[n1 - 1] = boxes1.clone();
            seq_boxes[n2 - 1] = boxes2.clone();
        }
    }

    // `cxx:2508-2529`: rebuild the face over `SeqWir` plus the untouched children.
    if is_done {
        let builder = TopoBuilder::new();
        let mut new_face = builder.make_face(surf.clone(), &seq_wir);
        new_face.0.set_orientation(Orientation::Forward); // cxx:2512
        // `myFace.EmptyCopied()` carries the location and the whole `BRep_TFace`
        // (surface, tolerance, natural restriction).
        new_face.0.set_location(face.0.location());
        if let Some(geom) = GeometryRegistry::global().face_geom(&face.0) {
            GeometryRegistry::global().set_face(&new_face.0, geom);
        }
        for s in &seq_nm {
            builder.add(&mut new_face.0, s); // cxx:2520-2524
        }
        new_face.0.set_orientation(ori); // cxx:2525
        ctx.replace(&face.0, &new_face.0); // cxx:2526
        *face = new_face; // cxx:2527
    }
    is_done
}
