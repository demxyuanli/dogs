//! Bi-point assembly of the outline update. Source: `HLRBRep_PolyAlgo.cxx:878-2978`
//! (`InitBiPointsWithConnexity`, the three `Interpolation` overloads and the
//! three `MoveOrInsertPoint` overloads).
//!
//! The C++ overloads of `Interpolation` / `MoveOrInsertPoint` differ only in how
//! many faces bound the edge run being split: one face (an edge with a single
//! adjacent face), two faces (an edge between two faces, one interpolated point
//! on each), and two faces with *two* interpolated points. The port keeps the
//! same three levels: [`HlrBRepPolyAlgo::interpolation_one`] /
//! [`HlrBRepPolyAlgo::move_or_insert_one`], `interpolation_two` /
//! `move_or_insert_two` and `move_or_insert_two_two`.
//!
//! The long in/out parameter lists (`theX1 .. theU2`) become [`Segment`], the
//! `(theNod11*, theNod12*, theI1p1, theI1p2, theI1, thePid1)` bundle of one side
//! becomes [`FaceSide`], and the `(X3..ZTI3, Coef3, U3, insP3, mP3P1)` out bundle
//! of one interpolation becomes [`Interp`]. OCCT hands out `NodeIndices&` /
//! `NodeData&` references into a reallocating array; the port passes the node
//! index instead and re-reads the slot, which is exactly what OCCT's post-
//! `AddNode` re-fetch (`cxx:2039-2041`, `:2349-2354`, `:2819-2830`) does.
//!
//! UNPORTED / deviations (each with its OCCT line):
//! * `BRep_Tool::PolygonOnTriangulation` needs a `BRep_PolygonOnTriangulation`
//!   representation on the edge (`cxx:900`, `:1070`, `:1079`). The port has the
//!   slot and the reader ([`crate::brep_tool::BRepTool::polygon_on_triangulation`])
//!   but no producer: `BRepMesh_ModelPostProcessor::PolygonCommitter`
//!   (`BRepMesh_ModelPostProcessor.cxx:90-155`) is not ported, so the connex
//!   arms behave exactly like OCCT does for an edge that was never meshed
//!   (null handle -> the `#ifdef OCCT_DEBUG` branch at `cxx:1055-1063`,
//!   `:1301-1309`).
//! * A polygon without parameters (`!HasParameters()`): OCCT dereferences the
//!   null `par` handle at `cxx:908` / `:928` (undefined behaviour); the port
//!   skips the pairing for that edge instead of reproducing a crash.
//! * `BRep_Tool::Continuity(theEdge, aF1, aF2)` (`cxx:1080`) has no port: the
//!   port's `BRep_TEdge` carries no `BRep_Continuity` representation list, so
//!   the `theRg` argument arrives as the already-decoded flag
//!   (`GeomAbs_G1` -> 1, `GeomAbs_G2` -> 2, `cxx:1512-1519`). A caller that has
//!   no continuity record passes 0, which is the `GeomAbs_C0` that
//!   `BRep_Tool::Continuity` answers when nothing matches.
//! * `myBCurv.Value(theU3)` (`cxx:1874`) needs the edge's 3D curve. When the
//!   edge has none, OCCT's `BRepAdaptor_Curve` has no curve to evaluate (the
//!   `Initialize` at `cxx:891` already left it null); the port answers "no
//!   interpolated point", i.e. the same `false` the enclosing `NewNode` failure
//!   arm produces.
//! * `TopoDS_Edge::Reverse()` (`cxx:1076`) only flips the local copy's
//!   orientation: the port threads the flag into the second polygon lookup
//!   instead of mutating a shape.
//! * `HLRAlgo_PolyInternalData::IncPINod` / `IncTData` / `IncPISeg` reallocate
//!   *two* arrays at once in OCCT (`PolyInternalData.cxx:923-1038`); the port's
//!   arrays belong to one face object, so the second out-parameter disappears
//!   (same deviation T-408 recorded).

use std::sync::Arc;

use occt_core::gp::GpXyz;
use occt_core::poly::Triangulation;
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;

use super::{nmask, HlrBRepPolyAlgo, PolyInternalDataArray};
use crate::hlr_algo::{
    HlrAlgoBiPoint, HlrAlgoPolyInternalData, PolyInternalNodeData, PolyInternalNodeIndices,
};
use crate::shape::{Edge, TopoShape};

/// The `theX1 .. theU2` block of `Interpolation` / `MoveOrInsertPoint`.
///
/// `x1/y1/z1` and `x2/y2/z2` are the **projected** endpoints (node space), the
/// `xt*/yt*/zt*` fields the 3D ones; `u1`/`u2` are the edge parameters. Every
/// field is in/out: `MoveOrInsertPoint` overwrites an endpoint when the
/// interpolated point snaps onto it (`cxx:1951-1959`, `:2001-2009`).
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct Segment {
    pub x1: f64,
    pub y1: f64,
    pub z1: f64,
    pub x2: f64,
    pub y2: f64,
    pub z2: f64,
    pub xt1: f64,
    pub yt1: f64,
    pub zt1: f64,
    pub xt2: f64,
    pub yt2: f64,
    pub zt2: f64,
    pub u1: f64,
    pub u2: f64,
}

/// One side of a bi-point pair: `(theNod1*, theNod2*, theI1p1, theI1p2, theI1,
/// thePid1)`. `face` is `theI1` (1-based) and doubles as the
/// [`PolyInternalDataArray`] slot, because `thePid1` is always `thePID(theI1)`
/// (`cxx:902`, `:1081`).
#[derive(Debug, Clone, Copy)]
pub(crate) struct FaceSide {
    pub face: usize,
    pub p1: i32,
    pub p2: i32,
}

/// The `(X3, Y3, Z3, XTI3, YTI3, ZTI3, Coef3, U3, insP3, mP3P1)` out bundle of
/// `Interpolation(theU1, theU2, ...)` (`cxx:1857-1884`).
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct Interp {
    /// Projected point (`TMultiply` of the 3D point, `cxx:1878`).
    pub x: f64,
    pub y: f64,
    pub z: f64,
    /// 3D point (`myBCurv.Value(theU3)`, `cxx:1874-1877`).
    pub xt: f64,
    pub yt: f64,
    pub zt: f64,
    pub coef: f64,
    pub u: f64,
    pub ins: bool,
    /// `mP3P1`: the crossing sits nearer `P1` than `P2`.
    pub m_p3_p1: bool,
}

/// `BRep_Tool::Triangulation(face, loc)` for a bare [`TopoShape`].
fn face_triangulation(face: &TopoShape) -> Option<Arc<Triangulation>> {
    face.tshape.read().unwrap().face_triangulation().cloned()
}

/// `BRepAdaptor_Curve2d(theEdge, theFace)` (`BRepAdaptor_Curve2d.cxx:39-53`):
/// the pcurve `BRep_Tool::CurveOnSurface` answers for that face, with the
/// `IsPeriodic()` / `Period()` the `U2` wrap at `cxx:961-964` needs.
///
/// `reversed` selects the second (reversed) pcurve of a seam edge, the way
/// `BRep_Tool::CurveOnSurface` does when the edge orientation and the face
/// orientation disagree.
fn pcurve_of(edge: &TopoShape, face: &TopoShape, reversed: bool) -> Option<Arc<dyn Curve2d>> {
    let face_key = crate::tgeometry::GeometryRegistry::shape_key(face);
    let pcs = crate::tgeometry::GeometryRegistry::global().edge_pcurves(edge, face_key);
    if pcs.is_empty() {
        return None;
    }
    let pos = if reversed && pcs.len() >= 2 { 1 } else { 0 };
    Some(pcs[pos].clone())
}

fn node_data(data: &HlrAlgoPolyInternalData, i: i32) -> PolyInternalNodeData {
    *data.pi_nod()[i as usize].data()
}

fn node_indices(data: &HlrAlgoPolyInternalData, i: i32) -> PolyInternalNodeIndices {
    *data.pi_nod()[i as usize].indices()
}

/// The `theNod1.RValues.Point = gp_XYZ(theX3, ...)`, `PCu1`/`PCu2`, `Scal = 0`
/// and `Flag |= NMsk_OutL` block that `MoveOrInsertPoint` repeats for every face
/// side (`cxx:1958-1977`, `:2008-2027`, and their two-face twins).
fn move_node_onto_point(
    data: &mut HlrAlgoPolyInternalData,
    node: i32,
    the_i_edge: i32,
    x3: f64,
    y3: f64,
    z3: f64,
    u3: f64,
) {
    let (edg1, edg2) = {
        let ind = node_indices(data, node);
        (ind.edg1, ind.edg2)
    };
    let nd = data.change_pi_nod()[node as usize].change_data();
    nd.point = GpXyz::new(x3, y3, z3);
    // `cxx:1959-1973`: the parameter lands on whichever slot holds this edge.
    if edg1 == the_i_edge {
        nd.pcu1 = u3;
    } else if edg2 == the_i_edge {
        nd.pcu2 = u3;
    }
    nd.scal = 0.0;
    data.change_pi_nod()[node as usize].change_indices().flag |= nmask::OUT_L;
}

impl HlrBRepPolyAlgo {
    /// `Interpolation(theU1, theU2, theNod1RValues, theNod2RValues, theP3,
    /// theX3..theZTI3, theCoef3, theU3, theIsNull, mP3P1)` (`cxx:1857-1884`):
    /// the crossing of the two nodes' `Scal` values, evaluated on the edge's 3D
    /// curve and projected with `TMultiply`.
    fn interpolation_point(
        &self,
        bcurv: &Option<Arc<dyn Curve>>,
        u1: f64,
        u2: f64,
        n1: &PolyInternalNodeData,
        n2: &PolyInternalNodeData,
    ) -> Interp {
        if let Some((coef, m_p3_p1)) = Self::new_node(n1.scal, n2.scal, self.tol_angular) {
            if let Some(c) = bcurv {
                let u3 = u1 + (u2 - u1) * coef;
                let p = c.value(u3);
                let (mut x, mut y, mut z) = (p.x(), p.y(), p.z());
                let (xt, yt, zt) = (x, y, z);
                // `cxx:1878-1883`: `TMultiply` on the 3D point, `XTI3` keeps the
                // untransformed one.
                self.tmultiply(&mut x, &mut y, &mut z, false);
                return Interp {
                    x,
                    y,
                    z,
                    xt,
                    yt,
                    zt,
                    coef,
                    u: u3,
                    ins: true,
                    m_p3_p1,
                };
            }
        }
        Interp::default()
    }

    /// `InitBiPointsWithConnexity(theIEdge, theEdge, theList, thePID, theLS,
    /// theIsConnex)` (`cxx:878-1377`).
    ///
    /// `the_ls` is the `theLS` face list of `TopExp::MapShapesAndAncestors`
    /// (`cxx:582`): the faces incident to the edge, empty when the edge is
    /// isolated (`cxx:588-590`). `the_i_edge` is `theIEdge`, the index the
    /// bi-points carry as `ShapeIndex`.
    pub fn init_bi_points_with_connexity(
        &mut self,
        the_i_edge: usize,
        the_edge: &TopoShape,
        the_list: &mut Vec<HlrAlgoBiPoint>,
        the_pid: &mut PolyInternalDataArray,
        the_ls: &[TopoShape],
        the_is_connex: bool,
        rg_flag: i32,
    ) {
        let edge = Edge(the_edge.clone());
        // `myBCurv.Initialize(theEdge)` (`cxx:891`).
        let bcurv = crate::brep_tool::BRepTool::edge_curve(&edge);
        let edge_reversed = the_edge.orientation().is_reversed();
        let mut seg = Segment::default();

        if the_is_connex {
            let nb_connex = the_ls.len();
            if nb_connex == 1 {
                let f1 = &the_ls[0];
                let i1 = self.f_map.find_index(f1);
                if i1 == 0 {
                    // `thePID.Value(i1)` with `i1 == 0` would be out of range in
                    // OCCT (`cxx:903`); a face that was never mapped cannot be
                    // one of `theLS`'s entries, so this is unreachable.
                    return;
                }
                let Some(tr1) = face_triangulation(f1) else {
                    return;
                };
                let Some(hpol) = crate::brep_tool::BRepTool::polygon_on_triangulation(
                    &edge,
                    &tr1,
                    edge_reversed,
                ) else {
                    return;
                };
                let Some(par) = hpol.parameters.as_ref() else {
                    return;
                };
                if i1 >= the_pid.len() || the_pid[i1].is_none() {
                    return;
                }
                let a_pol1: Vec<i32> = hpol.nodes.iter().map(|&n| n as i32).collect();
                let a_nb_pol = a_pol1.len();
                if a_nb_pol == 0 {
                    return;
                }
                // `myPC.Initialize(theEdge, aF1)` (`cxx:905`).
                let pc = pcurve_of(the_edge, f1, edge_reversed);

                {
                    let data = the_pid[i1].as_mut().unwrap();
                    // `cxx:917-919`: only the two ends get `NMsk_Vert`; for a
                    // closed polygon both writes hit the same node.
                    data.change_pi_nod()[a_pol1[0] as usize]
                        .change_indices()
                        .flag |= nmask::VERT;
                    data.change_pi_nod()[a_pol1[a_nb_pol - 1] as usize]
                        .change_indices()
                        .flag |= nmask::VERT;
                    // `cxx:921-935`.
                    for i_pol in 0..a_nb_pol {
                        let node = a_pol1[i_pol] as usize;
                        let par_i = par[i_pol];
                        let (edg1, _) = {
                            let ind = data.pi_nod()[node].indices();
                            (ind.edg1, ind.edg2)
                        };
                        if edg1 == 0 || edg1 == the_i_edge as i32 {
                            data.change_pi_nod()[node].change_indices().edg1 = the_i_edge as i32;
                            data.change_pi_nod()[node].change_data().pcu1 = par_i;
                        } else {
                            data.change_pi_nod()[node].change_indices().edg2 = the_i_edge as i32;
                            data.change_pi_nod()[node].change_data().pcu2 = par_i;
                        }
                    }
                }

                let mut i1p2 = a_pol1[0];
                // `cxx:938-940`: `aNode12Indices = aNode11Indices` and the same
                // for the values, i.e. node 12 is node 11 = `aPol1(1)`.
                let n12 = node_data(the_pid[i1].as_ref().unwrap(), i1p2);
                let ind12 = node_indices(the_pid[i1].as_ref().unwrap(), i1p2);
                seg.x2 = n12.point.x();
                seg.y2 = n12.point.y();
                seg.z2 = n12.point.z();
                seg.xt2 = seg.x2;
                seg.yt2 = seg.y2;
                seg.zt2 = seg.z2;
                if ind12.edg1 == the_i_edge as i32 {
                    seg.u2 = n12.pcu1;
                } else if ind12.edg2 == the_i_edge as i32 {
                    seg.u2 = n12.pcu2;
                }
                // `cxx:959-1027`.
                the_pid[i1].as_mut().unwrap().change_pi_nod()[i1p2 as usize]
                    .change_indices()
                    .flag |= nmask::EDGE;
                self.timultiply(&mut seg.xt2, &mut seg.yt2, &mut seg.zt2, false);
                if a_pol1[0] == a_pol1[a_nb_pol - 1] {
                    if let Some(pc) = &pc {
                        if pc.is_periodic() {
                            // `cxx:961-964`.
                            seg.u2 -= pc.period();
                        }
                    }
                }

                if a_nb_pol == 2 && crate::brep_tool::BRepTool::is_degenerated(&edge) {
                    // `cxx:966-971`.
                    let data = the_pid[i1].as_mut().unwrap();
                    Self::check_degenerated_segment(
                        data,
                        self.tol_angular,
                        a_pol1[0],
                        a_pol1[a_nb_pol - 1],
                    );
                    self.update_around_node(data, a_pol1[0]);
                    self.update_around_node(data, a_pol1[a_nb_pol - 1]);
                } else {
                    for i_pol in 1..a_nb_pol {
                        let i1p1 = i1p2;
                        i1p2 = a_pol1[i_pol];
                        // `cxx:999-1008`.
                        seg.x1 = seg.x2;
                        seg.y1 = seg.y2;
                        seg.z1 = seg.z2;
                        seg.xt1 = seg.xt2;
                        seg.yt1 = seg.yt2;
                        seg.zt1 = seg.zt2;
                        seg.u1 = seg.u2;
                        let n12 = node_data(the_pid[i1].as_ref().unwrap(), i1p2);
                        let ind12 = node_indices(the_pid[i1].as_ref().unwrap(), i1p2);
                        seg.x2 = n12.point.x();
                        seg.y2 = n12.point.y();
                        seg.z2 = n12.point.z();
                        seg.xt2 = seg.x2;
                        seg.yt2 = seg.y2;
                        seg.zt2 = seg.z2;
                        if ind12.edg1 == the_i_edge as i32 {
                            seg.u2 = n12.pcu1;
                        } else if ind12.edg2 == the_i_edge as i32 {
                            seg.u2 = n12.pcu2;
                        }
                        the_pid[i1].as_mut().unwrap().change_pi_nod()[i1p2 as usize]
                            .change_indices()
                            .flag |= nmask::EDGE;
                        self.timultiply(&mut seg.xt2, &mut seg.yt2, &mut seg.zt2, false);
                        let side = FaceSide {
                            face: i1,
                            p1: i1p1,
                            p2: i1p2,
                        };
                        self.interpolation_one(&bcurv, the_list, &mut seg, the_i_edge, side, the_pid);
                    }
                }
            } else if nb_connex == 2 {
                self.init_bi_points_two_faces(
                    the_i_edge,
                    &edge,
                    the_edge,
                    the_list,
                    the_pid,
                    the_ls,
                    &bcurv,
                    &mut seg,
                    rg_flag,
                );
            }
        } else {
            // `cxx:1311-1374`: no connexity, so the edge contributes its
            // `Poly_Polygon3D` segments directly.
            let Some(a_polyg) = crate::brep_tool::BRepTool::edge_polygon3d(&edge) else {
                return;
            };
            // `cxx:1315-1335`: `aTT = aLoc.Transformation();
            // aTT.PreMultiply(aProjTrsf)` read into `TTMa`/`TTLo`. The port's
            // edge location is the identity, so this is the projection alone.
            self.read_face_matrices(&occt_core::gp::GpTrsf::identity());
            let a_pol = &a_polyg.nodes;
            if a_pol.is_empty() {
                return;
            }
            let a_nb_pol1 = a_pol.len();
            let p1 = a_pol[0];
            seg.x2 = p1.x();
            seg.y2 = p1.y();
            seg.z2 = p1.z();
            self.ttmultiply(&mut seg.x2, &mut seg.y2, &mut seg.z2, false);
            seg.xt2 = seg.x2;
            seg.yt2 = seg.y2;
            seg.zt2 = seg.z2;
            self.timultiply(&mut seg.xt2, &mut seg.yt2, &mut seg.zt2, false);
            for j_pol in 1..a_nb_pol1 {
                seg.x1 = seg.x2;
                seg.y1 = seg.y2;
                seg.z1 = seg.z2;
                seg.xt1 = seg.xt2;
                seg.yt1 = seg.yt2;
                seg.zt1 = seg.zt2;
                let p2 = a_pol[j_pol];
                seg.x2 = p2.x();
                seg.y2 = p2.y();
                seg.z2 = p2.z();
                self.ttmultiply(&mut seg.x2, &mut seg.y2, &mut seg.z2, false);
                seg.xt2 = seg.x2;
                seg.yt2 = seg.y2;
                seg.zt2 = seg.z2;
                self.timultiply(&mut seg.xt2, &mut seg.yt2, &mut seg.zt2, false);
                // `cxx:1363-1364`: `Pnt*` are the 3D points (`XTI*`), `PntP*`
                // the projected ones (`X*`).
                the_list.insert(
                    0,
                    HlrAlgoBiPoint::with_flag(
                        seg.xt1, seg.yt1, seg.zt1, seg.xt2, seg.yt2, seg.zt2, seg.x1, seg.y1,
                        seg.z1, seg.x2, seg.y2, seg.z2, the_i_edge as i32, 0,
                    ),
                );
            }
        }
    }

    /// The `aNbConnex == 2` arm of `InitBiPointsWithConnexity`
    /// (`cxx:1064-1309`).
    #[allow(clippy::too_many_arguments)]
    fn init_bi_points_two_faces(
        &mut self,
        the_i_edge: usize,
        edge: &Edge,
        the_edge: &TopoShape,
        the_list: &mut Vec<HlrAlgoBiPoint>,
        the_pid: &mut PolyInternalDataArray,
        the_ls: &[TopoShape],
        bcurv: &Option<Arc<dyn Curve>>,
        seg: &mut Segment,
        rg_flag: i32,
    ) {
        let f1 = &the_ls[0];
        let i1 = self.f_map.find_index(f1);
        let orig_reversed = the_edge.orientation().is_reversed();
        let f2 = &the_ls[1];
        let i2 = self.f_map.find_index(f2);
        // `cxx:1075-1078`: the same face twice means the edge is its own seam:
        // OCCT reverses the local edge so the second lookup answers the other
        // polygon of the closed triangulation.
        let reversed2 = if i1 == i2 {
            !orig_reversed
        } else {
            orig_reversed
        };
        if i1 == 0 || i2 == 0 {
            return;
        }
        let Some(tr1) = face_triangulation(f1) else {
            return;
        };
        let Some(tr2) = face_triangulation(f2) else {
            return;
        };
        let pol1 = crate::brep_tool::BRepTool::polygon_on_triangulation(edge, &tr1, orig_reversed);
        let pol2 = crate::brep_tool::BRepTool::polygon_on_triangulation(edge, &tr2, reversed2);
        let (Some(hpol1), Some(hpol2)) = (pol1, pol2) else {
            return;
        };
        let Some(par) = hpol1.parameters.as_ref() else {
            return;
        };
        if i1 >= the_pid.len()
            || i2 >= the_pid.len()
            || the_pid[i1].is_none()
            || the_pid[i2].is_none()
        {
            return;
        }
        let a_pol1: Vec<i32> = hpol1.nodes.iter().map(|&n| n as i32).collect();
        let a_pol2: Vec<i32> = hpol2.nodes.iter().map(|&n| n as i32).collect();
        let a_nb_pol1 = a_pol1.len();
        if a_nb_pol1 == 0 || a_pol2.len() < a_nb_pol1 {
            return;
        }
        // `cxx:1085`: `myPC.Initialize(theEdge, aF1)` — after the reverse, so
        // the reversed pcurve of a seam edge is the one that is consulted.
        let pc = pcurve_of(the_edge, f1, orig_reversed);

        // `cxx:1109-1112`: four end nodes get `NMsk_Vert`.
        for (face, pol) in [(i1, &a_pol1), (i2, &a_pol2)] {
            let data = the_pid[face].as_mut().unwrap();
            data.change_pi_nod()[pol[0] as usize]
                .change_indices()
                .flag |= nmask::VERT;
            data.change_pi_nod()[pol[a_nb_pol1 - 1] as usize]
                .change_indices()
                .flag |= nmask::VERT;
        }
        // `cxx:1114-1143`: both polygons' nodes take the edge parameter.
        for i_pol in 0..a_nb_pol1 {
            let par_i = par[i_pol];
            for (face, pol) in [(i1, &a_pol1), (i2, &a_pol2)] {
                let data = the_pid[face].as_mut().unwrap();
                let node = pol[i_pol] as usize;
                let edg1 = data.pi_nod()[node].indices().edg1;
                if edg1 == 0 || edg1 == the_i_edge as i32 {
                    data.change_pi_nod()[node].change_indices().edg1 = the_i_edge as i32;
                    data.change_pi_nod()[node].change_data().pcu1 = par_i;
                } else {
                    data.change_pi_nod()[node].change_indices().edg2 = the_i_edge as i32;
                    data.change_pi_nod()[node].change_data().pcu2 = par_i;
                }
            }
        }

        let mut i1p2 = a_pol1[0];
        let mut i2p2 = a_pol2[0];
        // `cxx:1151-1159`: node 12 is node 11, node 22 is node 21.
        let n12 = node_data(the_pid[i1].as_ref().unwrap(), i1p2);
        let ind12 = node_indices(the_pid[i1].as_ref().unwrap(), i1p2);
        seg.x2 = n12.point.x();
        seg.y2 = n12.point.y();
        seg.z2 = n12.point.z();
        seg.xt2 = seg.x2;
        seg.yt2 = seg.y2;
        seg.zt2 = seg.z2;
        if ind12.edg1 == the_i_edge as i32 {
            seg.u2 = n12.pcu1;
        } else if ind12.edg2 == the_i_edge as i32 {
            seg.u2 = n12.pcu2;
        }
        // `cxx:1169-1171`.
        the_pid[i1].as_mut().unwrap().change_pi_nod()[i1p2 as usize]
            .change_indices()
            .flag |= nmask::EDGE;
        the_pid[i2].as_mut().unwrap().change_pi_nod()[i2p2 as usize]
            .change_indices()
            .flag |= nmask::EDGE;
        self.timultiply(&mut seg.xt2, &mut seg.yt2, &mut seg.zt2, false);
        if a_pol1[0] == a_pol1[a_nb_pol1 - 1] {
            if let Some(pc) = &pc {
                if pc.is_periodic() {
                    // `cxx:1172-1175`.
                    seg.u2 -= pc.period();
                }
            }
        }

        if a_nb_pol1 == 2 && crate::brep_tool::BRepTool::is_degenerated(edge) {
            // `cxx:1177-1186`.
            {
                let data = the_pid[i1].as_mut().unwrap();
                Self::check_degenerated_segment(data, self.tol_angular, a_pol1[0], a_pol1[1]);
                self.update_around_node(data, a_pol1[0]);
                self.update_around_node(data, a_pol1[1]);
            }
            {
                let data = the_pid[i2].as_mut().unwrap();
                Self::check_degenerated_segment(data, self.tol_angular, a_pol2[0], a_pol2[1]);
                self.update_around_node(data, a_pol2[0]);
                self.update_around_node(data, a_pol2[1]);
            }
        } else {
            for i_pol in 1..a_nb_pol1 {
                let i1p1 = i1p2;
                let i2p1 = i2p2;
                i1p2 = a_pol1[i_pol];
                i2p2 = a_pol2[i_pol];
                // `cxx:1231-1239`.
                seg.x1 = seg.x2;
                seg.y1 = seg.y2;
                seg.z1 = seg.z2;
                seg.xt1 = seg.xt2;
                seg.yt1 = seg.yt2;
                seg.zt1 = seg.zt2;
                seg.u1 = seg.u2;
                let n12 = node_data(the_pid[i1].as_ref().unwrap(), i1p2);
                let ind12 = node_indices(the_pid[i1].as_ref().unwrap(), i1p2);
                seg.x2 = n12.point.x();
                seg.y2 = n12.point.y();
                seg.z2 = n12.point.z();
                seg.xt2 = seg.x2;
                seg.yt2 = seg.y2;
                seg.zt2 = seg.z2;
                if ind12.edg1 == the_i_edge as i32 {
                    seg.u2 = n12.pcu1;
                } else if ind12.edg2 == the_i_edge as i32 {
                    seg.u2 = n12.pcu2;
                }
                // `cxx:1256-1257`: the flag is OR-ed twice in OCCT (a duplicated
                // line); one OR is the same state.
                the_pid[i1].as_mut().unwrap().change_pi_nod()[i1p2 as usize]
                    .change_indices()
                    .flag |= nmask::EDGE;
                self.timultiply(&mut seg.xt2, &mut seg.yt2, &mut seg.zt2, false);
                let side1 = FaceSide {
                    face: i1,
                    p1: i1p1,
                    p2: i1p2,
                };
                let side2 = FaceSide {
                    face: i2,
                    p1: i2p1,
                    p2: i2p2,
                };
                self.interpolation_two(
                    bcurv, the_list, seg, the_i_edge, side1, side2, the_pid, rg_flag,
                );
            }
        }
    }

    /// `HLRBRep_PolyAlgo::Interpolation(list, X1..ZTI2, theIEdge, U1, U2,
    /// Nod11, Nod12, I1p1, I1p2, I1, Pid1, ...)` (`cxx:1379-1463`), the
    /// one-face overload.
    fn interpolation_one(
        &mut self,
        bcurv: &Option<Arc<dyn Curve>>,
        list: &mut Vec<HlrAlgoBiPoint>,
        seg: &mut Segment,
        the_i_edge: usize,
        f1: FaceSide,
        the_pid: &mut PolyInternalDataArray,
    ) {
        let n11 = node_data(the_pid[f1.face].as_ref().unwrap(), f1.p1);
        let n12 = node_data(the_pid[f1.face].as_ref().unwrap(), f1.p2);
        let p3 = self.interpolation_point(bcurv, seg.u1, seg.u2, &n11, &n12);
        self.move_or_insert_one(list, seg, the_i_edge as i32, f1, the_pid, &p3, 0);
    }

    /// The one-face `Interpolation` overload of `MoveOrInsertPoint`
    /// (`cxx:1888-2112`).
    #[allow(clippy::too_many_arguments)]
    fn move_or_insert_one(
        &mut self,
        list: &mut Vec<HlrAlgoBiPoint>,
        seg: &mut Segment,
        the_i_edge: i32,
        f1: FaceSide,
        the_pid: &mut PolyInternalDataArray,
        p3: &Interp,
        the_flag: i32,
    ) {
        let mut an_ins3 = p3.ins;
        if an_ins3 && p3.m_p3_p1 {
            // `cxx:1935-1980`: P1 ---> P3.
            let is_vert =
                node_indices(the_pid[f1.face].as_ref().unwrap(), f1.p1).flag & nmask::VERT != 0;
            if !is_vert && p3.coef < self.tol_sta {
                an_ins3 = false;
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    self.change_node(data, f1.p1, f1.p2, p3.coef, p3.x, p3.y, p3.z, true);
                    move_node_onto_point(data, f1.p1, the_i_edge, p3.x, p3.y, p3.z, p3.u);
                    self.update_around_node(data, f1.p1);
                }
                seg.x1 = p3.x;
                seg.y1 = p3.y;
                seg.z1 = p3.z;
                seg.xt1 = p3.xt;
                seg.yt1 = p3.yt;
                seg.zt1 = p3.zt;
                seg.u1 = p3.u;
                if let Some(front) = list.first_mut() {
                    // `cxx:1978-1979`.
                    let pts = front.change_points();
                    pts.pntp2 = GpXyz::new(p3.x, p3.y, p3.z);
                    pts.pnt2 = GpXyz::new(p3.xt, p3.yt, p3.zt);
                }
            }
        }
        if an_ins3 && !p3.m_p3_p1 {
            // `cxx:1981-2028`: P2 ---> P3.
            let is_vert =
                node_indices(the_pid[f1.face].as_ref().unwrap(), f1.p2).flag & nmask::VERT != 0;
            if !is_vert && p3.coef > self.tol_end {
                an_ins3 = false;
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    self.change_node(data, f1.p1, f1.p2, p3.coef, p3.x, p3.y, p3.z, false);
                    move_node_onto_point(data, f1.p2, the_i_edge, p3.x, p3.y, p3.z, p3.u);
                    self.update_around_node(data, f1.p2);
                }
                seg.x2 = p3.x;
                seg.y2 = p3.y;
                seg.z2 = p3.z;
                seg.xt2 = p3.xt;
                seg.yt2 = p3.yt;
                seg.zt2 = p3.zt;
                seg.u2 = p3.u;
            }
        }
        if an_ins3 {
            // `cxx:2029-2090`: p1 i1p3 p2.
            let an_i1p3 = {
                let data = the_pid[f1.face].as_mut().unwrap();
                let n11 = node_data(data, f1.p1);
                let n12 = node_data(data, f1.p2);
                let new_node = data.add_node(&n11, &n12, p3.coef, p3.x, p3.y, p3.z);
                let ind = data.change_pi_nod()[new_node as usize].change_indices();
                ind.edg1 = the_i_edge;
                ind.flag |= nmask::OUT_L | nmask::EDGE;
                let nd = data.change_pi_nod()[new_node as usize].change_data();
                nd.pcu1 = p3.u;
                nd.scal = 0.0;
                data.update_links_split(f1.p1, f1.p2, new_node);
                self.update_around_node(data, new_node);
                new_node
            };
            list.insert(
                0,
                HlrAlgoBiPoint::from_face1_with_flag(
                    seg.xt1, seg.yt1, seg.zt1, p3.xt, p3.yt, p3.zt, seg.x1, seg.y1, seg.z1, p3.x,
                    p3.y, p3.z, the_i_edge, f1.face as i32, f1.p1, an_i1p3, the_flag,
                ),
            );
            list.insert(
                0,
                HlrAlgoBiPoint::from_face1_with_flag(
                    p3.xt, p3.yt, p3.zt, seg.xt2, seg.yt2, seg.zt2, p3.x, p3.y, p3.z, seg.x2, seg.y2,
                    seg.z2, the_i_edge, f1.face as i32, an_i1p3, f1.p2, the_flag,
                ),
            );
        } else {
            // `cxx:2092-2111`: p1 p2.
            list.insert(
                0,
                HlrAlgoBiPoint::from_face1_with_flag(
                    seg.xt1, seg.yt1, seg.zt1, seg.xt2, seg.yt2, seg.zt2, seg.x1, seg.y1, seg.z1,
                    seg.x2, seg.y2, seg.z2, the_i_edge, f1.face as i32, f1.p1, f1.p2, the_flag,
                ),
            );
        }
    }

    /// `Interpolation(..., theRg, Nod11, Nod12, ..., Nod21, Nod22, ...)`
    /// (`cxx:1466-1853`): the two-face overload. One interpolated point per
    /// face, then the crossing nearest `P1` (or both) is inserted.
    #[allow(clippy::too_many_arguments)]
    fn interpolation_two(
        &mut self,
        bcurv: &Option<Arc<dyn Curve>>,
        list: &mut Vec<HlrAlgoBiPoint>,
        seg: &mut Segment,
        the_i_edge: usize,
        f1: FaceSide,
        f2: FaceSide,
        the_pid: &mut PolyInternalDataArray,
        rg_flag: i32,
    ) {
        // `cxx:1512-1519`: `flag += 1` for G1, `+= 2` for G2.
        let flag = rg_flag;
        let n11 = node_data(the_pid[f1.face].as_ref().unwrap(), f1.p1);
        let n12 = node_data(the_pid[f1.face].as_ref().unwrap(), f1.p2);
        let p3 = self.interpolation_point(bcurv, seg.u1, seg.u2, &n11, &n12);
        let n21 = node_data(the_pid[f2.face].as_ref().unwrap(), f2.p1);
        let n22 = node_data(the_pid[f2.face].as_ref().unwrap(), f2.p2);
        let p4 = self.interpolation_point(bcurv, seg.u1, seg.u2, &n21, &n22);
        if p3.ins || p4.ins {
            if !p4.ins {
                // p1 i1p3 p2
                self.move_or_insert_two(
                    list,
                    seg,
                    the_i_edge as i32,
                    f1,
                    f2,
                    the_pid,
                    &p3,
                    flag,
                );
            } else if !p3.ins {
                // p1 i2p4 p2 — the sides swap, so `f2` plays the "first" role.
                self.move_or_insert_two(
                    list,
                    seg,
                    the_i_edge as i32,
                    f2,
                    f1,
                    the_pid,
                    &p4,
                    flag,
                );
            } else if (p4.coef - p3.coef).abs() < self.tol_sta {
                // `cxx:1654`: the two crossings coincide -> keep the face-2 one.
                self.move_or_insert_two(
                    list,
                    seg,
                    the_i_edge as i32,
                    f2,
                    f1,
                    the_pid,
                    &p4,
                    flag,
                );
            } else if p4.coef < p3.coef {
                // p1 i2p4 i1p3 p2
                self.move_or_insert_two_two(
                    list,
                    seg,
                    the_i_edge as i32,
                    f2,
                    f1,
                    the_pid,
                    &p4,
                    &p3,
                    flag,
                );
            } else {
                // p1 i1p3 i2p4 p2
                self.move_or_insert_two_two(
                    list,
                    seg,
                    the_i_edge as i32,
                    f1,
                    f2,
                    the_pid,
                    &p3,
                    &p4,
                    flag,
                );
            }
        } else {
            // `cxx:1830-1851`: p1 p2.
            list.insert(
                0,
                HlrAlgoBiPoint::from_two_faces_with_flag(
                    seg.xt1,
                    seg.yt1,
                    seg.zt1,
                    seg.xt2,
                    seg.yt2,
                    seg.zt2,
                    seg.x1,
                    seg.y1,
                    seg.z1,
                    seg.x2,
                    seg.y2,
                    seg.z2,
                    the_i_edge as i32,
                    f1.face as i32,
                    f1.p1,
                    f1.p2,
                    f2.face as i32,
                    f2.p1,
                    f2.p2,
                    flag,
                ),
            );
        }
    }

    /// The two-face, one-point `MoveOrInsertPoint` (`cxx:2116-2449`).
    ///
    /// `f1`/`f2` are `(theNod1*, theI1*)` and `(theNod2*, theI2*)`; the point is
    /// `theP3` with its bundle.
    #[allow(clippy::too_many_arguments)]
    fn move_or_insert_two(
        &mut self,
        list: &mut Vec<HlrAlgoBiPoint>,
        seg: &mut Segment,
        the_i_edge: i32,
        f1: FaceSide,
        f2: FaceSide,
        the_pid: &mut PolyInternalDataArray,
        p3: &Interp,
        the_flag: i32,
    ) {
        let mut an_ins3 = p3.ins;
        if an_ins3 && p3.m_p3_p1 {
            // `cxx:2168-2249`: P1 ---> P3 on both faces.
            let is_vert =
                node_indices(the_pid[f1.face].as_ref().unwrap(), f1.p1).flag & nmask::VERT != 0;
            if !is_vert && p3.coef < self.tol_sta {
                an_ins3 = false;
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    self.change_node(data, f1.p1, f1.p2, p3.coef, p3.x, p3.y, p3.z, true);
                    move_node_onto_point(data, f1.p1, the_i_edge, p3.x, p3.y, p3.z, p3.u);
                    self.update_around_node(data, f1.p1);
                }
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    self.change_node(data, f2.p1, f2.p2, p3.coef, p3.x, p3.y, p3.z, true);
                    move_node_onto_point(data, f2.p1, the_i_edge, p3.x, p3.y, p3.z, p3.u);
                    self.update_around_node(data, f2.p1);
                }
                seg.x1 = p3.x;
                seg.y1 = p3.y;
                seg.z1 = p3.z;
                seg.xt1 = p3.xt;
                seg.yt1 = p3.yt;
                seg.zt1 = p3.zt;
                seg.u1 = p3.u;
                if let Some(front) = list.first_mut() {
                    let pts = front.change_points();
                    pts.pntp2 = GpXyz::new(p3.x, p3.y, p3.z);
                    pts.pnt2 = GpXyz::new(p3.xt, p3.yt, p3.zt);
                }
            }
        }
        if an_ins3 && !p3.m_p3_p1 {
            // `cxx:2251-2329`: P2 ---> P3 on both faces.
            let is_vert =
                node_indices(the_pid[f1.face].as_ref().unwrap(), f1.p2).flag & nmask::VERT != 0;
            if !is_vert && p3.coef > self.tol_end {
                an_ins3 = false;
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    self.change_node(data, f1.p1, f1.p2, p3.coef, p3.x, p3.y, p3.z, false);
                    move_node_onto_point(data, f1.p2, the_i_edge, p3.x, p3.y, p3.z, p3.u);
                    self.update_around_node(data, f1.p2);
                }
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    self.change_node(data, f2.p1, f2.p2, p3.coef, p3.x, p3.y, p3.z, false);
                    move_node_onto_point(data, f2.p2, the_i_edge, p3.x, p3.y, p3.z, p3.u);
                    self.update_around_node(data, f2.p2);
                }
                seg.x2 = p3.x;
                seg.y2 = p3.y;
                seg.z2 = p3.z;
                seg.xt2 = p3.xt;
                seg.yt2 = p3.yt;
                seg.zt2 = p3.zt;
                seg.u2 = p3.u;
            }
        }
        if an_ins3 {
            // `cxx:2331-2423`: p1 i1p3 p2 with a node on each face.
            let (an_i1p3, an_i2p3) = {
                let data1 = the_pid[f1.face].as_mut().unwrap();
                let n11 = node_data(data1, f1.p1);
                let n12 = node_data(data1, f1.p2);
                let new1 = data1.add_node(&n11, &n12, p3.coef, p3.x, p3.y, p3.z);
                data1.change_pi_nod()[new1 as usize].change_indices().edg1 = the_i_edge;
                data1.change_pi_nod()[new1 as usize].change_indices().flag |=
                    nmask::OUT_L | nmask::EDGE;
                data1.change_pi_nod()[new1 as usize].change_data().pcu1 = p3.u;
                data1.change_pi_nod()[new1 as usize].change_data().scal = 0.0;
                let data2 = the_pid[f2.face].as_mut().unwrap();
                let n21 = node_data(data2, f2.p1);
                let n22 = node_data(data2, f2.p2);
                let new2 = data2.add_node(&n21, &n22, p3.coef, p3.x, p3.y, p3.z);
                data2.change_pi_nod()[new2 as usize].change_indices().edg1 = the_i_edge;
                data2.change_pi_nod()[new2 as usize].change_indices().flag |=
                    nmask::OUT_L | nmask::EDGE;
                data2.change_pi_nod()[new2 as usize].change_data().pcu1 = p3.u;
                data2.change_pi_nod()[new2 as usize].change_data().scal = 0.0;
                (new1, new2)
            };
            {
                let data = the_pid[f1.face].as_mut().unwrap();
                data.update_links_split(f1.p1, f1.p2, an_i1p3);
                self.update_around_node(data, an_i1p3);
            }
            {
                let data = the_pid[f2.face].as_mut().unwrap();
                data.update_links_split(f2.p1, f2.p2, an_i2p3);
                self.update_around_node(data, an_i2p3);
            }
            list.insert(
                0,
                HlrAlgoBiPoint::from_two_faces_with_flag(
                    seg.xt1,
                    seg.yt1,
                    seg.zt1,
                    p3.xt,
                    p3.yt,
                    p3.zt,
                    seg.x1,
                    seg.y1,
                    seg.z1,
                    p3.x,
                    p3.y,
                    p3.z,
                    the_i_edge,
                    f1.face as i32,
                    f1.p1,
                    an_i1p3,
                    f2.face as i32,
                    f2.p1,
                    an_i2p3,
                    the_flag,
                ),
            );
            list.insert(
                0,
                HlrAlgoBiPoint::from_two_faces_with_flag(
                    p3.xt,
                    p3.yt,
                    p3.zt,
                    seg.xt2,
                    seg.yt2,
                    seg.zt2,
                    p3.x,
                    p3.y,
                    p3.z,
                    seg.x2,
                    seg.y2,
                    seg.z2,
                    the_i_edge,
                    f1.face as i32,
                    an_i1p3,
                    f1.p2,
                    f2.face as i32,
                    an_i2p3,
                    f2.p2,
                    the_flag,
                ),
            );
        } else {
            // `cxx:2426-2448`: p1 p2.
            list.insert(
                0,
                HlrAlgoBiPoint::from_two_faces_with_flag(
                    seg.xt1,
                    seg.yt1,
                    seg.zt1,
                    seg.xt2,
                    seg.yt2,
                    seg.zt2,
                    seg.x1,
                    seg.y1,
                    seg.z1,
                    seg.x2,
                    seg.y2,
                    seg.z2,
                    the_i_edge,
                    f1.face as i32,
                    f1.p1,
                    f1.p2,
                    f2.face as i32,
                    f2.p1,
                    f2.p2,
                    the_flag,
                ),
            );
        }
    }

    /// The two-face, two-point `MoveOrInsertPoint` (`cxx:2453-2976`).
    ///
    /// `f1`/`f2` pair with `p3`/`p4`; the snapping arms touch `f1`+`f2` for
    /// `P1 -> P3`, `f1`+`f2` for `P2 -> P4`, and then fall back to the
    /// one-point overloads when only one of the two points survived.
    #[allow(clippy::too_many_arguments)]
    fn move_or_insert_two_two(
        &mut self,
        list: &mut Vec<HlrAlgoBiPoint>,
        seg: &mut Segment,
        the_i_edge: i32,
        f1: FaceSide,
        f2: FaceSide,
        the_pid: &mut PolyInternalDataArray,
        p3: &Interp,
        p4: &Interp,
        the_flag: i32,
    ) {
        let mut an_ins3 = p3.ins;
        let mut an_ins4 = p4.ins;
        if an_ins3 && p3.m_p3_p1 {
            // `cxx:2516-2597`: P1 ---> P3.
            let is_vert =
                node_indices(the_pid[f1.face].as_ref().unwrap(), f1.p1).flag & nmask::VERT != 0;
            if !is_vert && p3.coef < self.tol_sta {
                an_ins3 = false;
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    self.change_node(data, f1.p1, f1.p2, p3.coef, p3.x, p3.y, p3.z, true);
                    move_node_onto_point(data, f1.p1, the_i_edge, p3.x, p3.y, p3.z, p3.u);
                    self.update_around_node(data, f1.p1);
                }
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    self.change_node(data, f2.p1, f2.p2, p3.coef, p3.x, p3.y, p3.z, true);
                    move_node_onto_point(data, f2.p1, the_i_edge, p3.x, p3.y, p3.z, p3.u);
                    self.update_around_node(data, f2.p1);
                }
                seg.x1 = p3.x;
                seg.y1 = p3.y;
                seg.z1 = p3.z;
                seg.xt1 = p3.xt;
                seg.yt1 = p3.yt;
                seg.zt1 = p3.zt;
                seg.u1 = p3.u;
                if let Some(front) = list.first_mut() {
                    let pts = front.change_points();
                    pts.pntp2 = GpXyz::new(p3.x, p3.y, p3.z);
                    pts.pnt2 = GpXyz::new(p3.xt, p3.yt, p3.zt);
                }
            }
        }
        if an_ins4 && !p4.m_p3_p1 {
            // `cxx:2599-2678`: P2 ---> P4.
            let is_vert =
                node_indices(the_pid[f1.face].as_ref().unwrap(), f1.p2).flag & nmask::VERT != 0;
            if !is_vert && p4.coef > self.tol_end {
                an_ins4 = false;
                // `cxx:2603-2630`: `ChangeNode` on the face-2 side first, then on
                // the face-1 side; each already re-orients its own node.
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    self.change_node(data, f2.p1, f2.p2, p4.coef, p4.x, p4.y, p4.z, false);
                }
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    self.change_node(data, f1.p1, f1.p2, p4.coef, p4.x, p4.y, p4.z, false);
                }
                // `cxx:2639-2665`: the parameter writes of the two sides.
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    move_node_onto_point(data, f1.p2, the_i_edge, p4.x, p4.y, p4.z, p4.u);
                    self.update_around_node(data, f1.p2);
                }
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    move_node_onto_point(data, f2.p2, the_i_edge, p4.x, p4.y, p4.z, p4.u);
                    self.update_around_node(data, f2.p2);
                }
                seg.x2 = p4.x;
                seg.y2 = p4.y;
                seg.z2 = p4.z;
                seg.xt2 = p4.xt;
                seg.yt2 = p4.yt;
                seg.zt2 = p4.zt;
                seg.u2 = p4.u;
            }
        }
        if an_ins3 || an_ins4 {
            if !an_ins4 {
                // `cxx:2681-2733`: p1 i1p3 p2.
                self.move_or_insert_two(
                    list, seg, the_i_edge, f1, f2, the_pid, p3, the_flag,
                );
            } else if !an_ins3 {
                // `cxx:2734-2786`: p1 i2p4 p2.
                self.move_or_insert_two(
                    list, seg, the_i_edge, f2, f1, the_pid, p4, the_flag,
                );
            } else {
                // `cxx:2787-2952`: p1 i1p3 i2p4 p2, a node on each face for each
                // of the two points. `AddNode` runs face-1/face-2/face-1/face-2 in
                // that order (`cxx:2787-2815`); the two nodes of one face keep
                // their relative numbering, which is all that matters.
                let an_i1p3 = {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    let n11 = node_data(data, f1.p1);
                    let n12 = node_data(data, f1.p2);
                    data.add_node(&n11, &n12, p3.coef, p3.x, p3.y, p3.z)
                };
                let an_i2p3 = {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    let n21 = node_data(data, f2.p1);
                    let n22 = node_data(data, f2.p2);
                    data.add_node(&n21, &n22, p3.coef, p3.x, p3.y, p3.z)
                };
                let an_i1p4 = {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    let n11 = node_data(data, f1.p1);
                    let n12 = node_data(data, f1.p2);
                    data.add_node(&n11, &n12, p4.coef, p4.x, p4.y, p4.z)
                };
                let an_i2p4 = {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    let n21 = node_data(data, f2.p1);
                    let n22 = node_data(data, f2.p2);
                    data.add_node(&n21, &n22, p4.coef, p4.x, p4.y, p4.z)
                };
                // `cxx:2831-2851`: `Edg1`/`PCu1`/`Scal = 0`/`OutL`/`Edge` per new
                // node, face 1 first.
                for (face, node, u) in [
                    (f1.face, an_i1p3, p3.u),
                    (f2.face, an_i2p3, p3.u),
                    (f1.face, an_i1p4, p4.u),
                    (f2.face, an_i2p4, p4.u),
                ] {
                    let data = the_pid[face].as_mut().unwrap();
                    data.change_pi_nod()[node as usize].change_indices().edg1 = the_i_edge;
                    data.change_pi_nod()[node as usize].change_indices().flag |=
                        nmask::OUT_L | nmask::EDGE;
                    data.change_pi_nod()[node as usize].change_data().pcu1 = u;
                    data.change_pi_nod()[node as usize].change_data().scal = 0.0;
                }
                // `cxx:2853-2888`: the four splits, in OCCT's order.
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    data.update_links_split(f1.p1, f1.p2, an_i1p3);
                }
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    data.update_links_split(f2.p1, f2.p2, an_i2p3);
                }
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    data.update_links_split(an_i2p3, f2.p2, an_i2p4);
                }
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    data.update_links_split(an_i1p3, f1.p2, an_i1p4);
                }
                // `cxx:2888-2890`.
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    self.update_around_node(data, an_i1p3);
                }
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    self.update_around_node(data, an_i2p3);
                }
                {
                    let data = the_pid[f1.face].as_mut().unwrap();
                    self.update_around_node(data, an_i1p4);
                }
                {
                    let data = the_pid[f2.face].as_mut().unwrap();
                    self.update_around_node(data, an_i2p4);
                }
                list.insert(
                    0,
                    HlrAlgoBiPoint::from_two_faces_with_flag(
                        seg.xt1,
                        seg.yt1,
                        seg.zt1,
                        p3.xt,
                        p3.yt,
                        p3.zt,
                        seg.x1,
                        seg.y1,
                        seg.z1,
                        p3.x,
                        p3.y,
                        p3.z,
                        the_i_edge,
                        f1.face as i32,
                        f1.p1,
                        an_i1p3,
                        f2.face as i32,
                        f2.p1,
                        an_i2p3,
                        the_flag,
                    ),
                );
                list.insert(
                    0,
                    HlrAlgoBiPoint::from_two_faces_with_flag(
                        p3.xt,
                        p3.yt,
                        p3.zt,
                        p4.xt,
                        p4.yt,
                        p4.zt,
                        p3.x,
                        p3.y,
                        p3.z,
                        p4.x,
                        p4.y,
                        p4.z,
                        the_i_edge,
                        f1.face as i32,
                        an_i1p3,
                        an_i1p4,
                        f2.face as i32,
                        an_i2p3,
                        an_i2p4,
                        the_flag,
                    ),
                );
                list.insert(
                    0,
                    HlrAlgoBiPoint::from_two_faces_with_flag(
                        p4.xt,
                        p4.yt,
                        p4.zt,
                        seg.xt2,
                        seg.yt2,
                        seg.zt2,
                        p4.x,
                        p4.y,
                        p4.z,
                        seg.x2,
                        seg.y2,
                        seg.z2,
                        the_i_edge,
                        f1.face as i32,
                        an_i1p4,
                        f1.p2,
                        f2.face as i32,
                        an_i2p4,
                        f2.p2,
                        the_flag,
                    ),
                );
            }
        } else {
            // `cxx:2956-2976`: p1 p2.
            list.insert(
                0,
                HlrAlgoBiPoint::from_two_faces_with_flag(
                    seg.xt1,
                    seg.yt1,
                    seg.zt1,
                    seg.xt2,
                    seg.yt2,
                    seg.zt2,
                    seg.x1,
                    seg.y1,
                    seg.z1,
                    seg.x2,
                    seg.y2,
                    seg.z2,
                    the_i_edge,
                    f1.face as i32,
                    f1.p1,
                    f1.p2,
                    f2.face as i32,
                    f2.p1,
                    f2.p2,
                    the_flag,
                ),
            );
        }
    }
}
