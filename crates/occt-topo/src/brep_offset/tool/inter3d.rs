//! `BRepOffset_Tool::Inter3D` (`BRepOffset_Tool.cxx:1441-1941`): the 3D
//! intersection of two faces, answering the section edges oriented on each
//! face (`L1` / `L2`).
//!
//! Control flow, in the reference's order:
//! 1. two untrimmed planes take the direct `PerformPlanes` path
//!    (`cxx:1451-1468`);
//! 2. `BRepLib::BuildCurves3d` on both faces and `UpdateVertexTolerances`
//!    (`cxx:1470-1475`);
//! 3. a `BOPAlgo_PaveFiller` runs on the two faces alone
//!    (`cxx:1477-1481`); with a reference edge, `CheckIntersFF` narrows the
//!    section edges to the nearest connexity block (`cxx:1483-1487`);
//! 4. every section pave block of every `BOPDS_InterfFF` curve is re-read
//!    (`cxx:1504-1591`): the pcurve on a face that has none is attached
//!    after `BOPTools_AlgoTools2D::AdjustPCurveOnFace`, and
//!    `OrientSection` gives the two orientations (`OUT` reverses both);
//! 5. `addPCurve1/2` are cleared for planar faces and `isEl1/2` set for the
//!    elementary ones (`cxx:1594-1623`);
//! 6. with more than one section edge and a non-elementary pair, the edges
//!    are glued (`AssembleEdge`, `cxx:1705-1889`): from `TrueEdges` when the
//!    reference edge was given, otherwise by chaining `L1`'s edges into wires
//!    (`AreConnex` / `AngleWireEdge`, `cxx:1728-1787`) and ordering each
//!    wire's edges from a non-autonomous vertex (`cxx:1793-1867`) or with
//!    `BRepTools_WireExplorer` (`cxx:1868-1873`);
//! 7. when the assembly shortened the list, `L1` / `L2` are rebuilt from it
//!    with `BRepLib::SameParameter`, `ReconstructPCurves` and a second
//!    `SameParameter` for tolerances of at least `1.e-4`
//!    (`cxx:1892-1929`);
//! 8. otherwise every remaining edge only takes `BRepLib::SameParameter`
//!    (`cxx:1931-1941`).
//!
//! Port guards (no OCCT counterpart): a pave-block edge index that does not
//! resolve to an edge, a failed pave filler and an `OrientSection` that cannot
//! read both pcurves are skipped instead of dereferencing a null handle; a
//! closed wire whose edges are all autonomous (the reference reads the first
//! element of an empty `Elist`, `cxx:1835-1841`) is passed through
//! `BRepTools_WireExplorer` instead, and a chain step that finds no matching
//! vertex stops the chain (the reference loops forever).
//!
//! `BRepTools::UVBounds` for the untrimmed-plane test is
//! [`crate::brep_tools_uv::uv_bounds_face`], i.e. the pcurve-based bounds that
//! `BRepAdaptor_Surface::Initialize(F, true)` uses, not the surface's own
//! infinite range - otherwise every planar face would look untrimmed.
//!
//! `BRep_Builder::UpdateEdge(anEdge, aC2d, F, aTolEdge)` is
//! [`GeometryRegistry::set_edge_pcurve`] plus a raised edge tolerance;
//! `AdjustPCurveOnFace`'s two overloads (`BOPTools_AlgoTools2D.cxx:209-244`)
//! both reduce to `AdjustPCurveOnSurf` over `[f, l]`, which is
//! [`AlgoTools2D::adjust_pcurve_on_surf`].

use std::sync::Arc;

use occt_core::gp::GpVec;
use occt_core::precision::CONFUSION;
use occt_geom::extrema_pc::extrema_ext_pc_range_done;
use occt_geom2d::Curve2d;

use crate::abs::State;
use crate::algo_tools::AlgoTools2D;
use crate::boptools_2d;
use crate::brep_lib_build_curves3d::build_curves3d;
use crate::brep_lib_same_parameter::same_parameter_edge_inplace;
use crate::brep_tool::BRepTool;
use crate::brep_tools_uv::uv_bounds_face;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};
use crate::pave_filler::PaveFiller;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edges_of_wire, is_same};

use super::BRepOffsetTool;

/// True for the five `Geom_ElementarySurface` subclasses
/// (`Geom_Plane` / `Geom_CylindricalSurface` / `Geom_ConicalSurface` /
/// `Geom_SphericalSurface` / `Geom_ToroidalSurface`): OCCT's
/// `IsKind(STANDARD_TYPE(Geom_ElementarySurface))` (`cxx:1606`, `:1618`).
fn is_elementary(kind: GeomAbsSurfaceType) -> bool {
    matches!(
        kind,
        GeomAbsSurfaceType::Plane
            | GeomAbsSurfaceType::Cylinder
            | GeomAbsSurfaceType::Cone
            | GeomAbsSurfaceType::Sphere
            | GeomAbsSurfaceType::Torus
    )
}

/// `V1.IsSame(V2)` on two possibly null vertices: two null shapes are `IsSame`
/// (`TopoDS_Shape::IsSame`).
fn same_vertex(a: &Option<Vertex>, b: &Option<Vertex>) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => is_same(&x.0, &y.0),
        (None, None) => true,
        _ => false,
    }
}

/// `IsAutonomVertex(V, PDS)` on a possibly null vertex; an absent vertex is
/// not in the DS, which the two-argument overload answers with `false`
/// (`cxx:839-911`).
fn is_autonom(vertex: Option<&Vertex>, pds: &crate::bopds::BopdsDS) -> bool {
    match vertex {
        Some(v) => BRepOffsetTool::is_autonom_vertex(v, pds),
        None => false,
    }
}

/// `BRep_Builder::UpdateEdge(E, C, F, Tol)`: attach the pcurve and raise the
/// edge tolerance to `tol`.
fn update_edge_pcurve(edge: &Edge, face: &Face, pcurve: Arc<dyn Curve2d>, tol: f64) {
    let reg = GeometryRegistry::global();
    reg.set_edge_pcurve(&edge.0, GeometryRegistry::shape_key(&face.0), pcurve);
    if BRepTool::edge_tolerance(edge) < tol {
        reg.set_edge_tolerance(&edge.0, tol);
    }
}

impl BRepOffsetTool {
    /// `Inter3D(F1, F2, L1, L2, Side, RefEdge, theRefFace1, theRefFace2)`
    /// (`cxx:1441-1941`).
    ///
    /// `l1` and `l2` are cleared first and receive the section edges oriented
    /// on `F1` and on `F2`.
    #[allow(clippy::too_many_arguments)]
    pub fn inter3d(
        f1: &Face,
        f2: &Face,
        l1: &mut Vec<TopoShape>,
        l2: &mut Vec<TopoShape>,
        side: State,
        ref_edge: Option<&Edge>,
        ref_face1: Option<&Face>,
        ref_face2: Option<&Face>,
    ) {
        // `cxx:1451-1468`: planar and not trimmed - `BRepAdaptor_Surface` with
        // `Restriction = false` first (its own bounds), then `Initialize(F,
        // true)` and `IsInf(LastUParameter())` / `IsInf(LastVParameter())`.
        let surf1 = BRepTool::face_surface(f1);
        let surf2 = BRepTool::face_surface(f2);
        if let (Some(s1), Some(s2)) = (&surf1, &surf2) {
            if surface_kind(s1.as_ref()) == GeomAbsSurfaceType::Plane
                && surface_kind(s2.as_ref()) == GeomAbsSurfaceType::Plane
            {
                let (_, umax1, _, vmax1) = uv_bounds_face(f1);
                if Self::is_inf(umax1) && Self::is_inf(vmax1) {
                    let (_, umax2, _, vmax2) = uv_bounds_face(f2);
                    if Self::is_inf(umax2) && Self::is_inf(vmax2) {
                        Self::perform_planes(f1, f2, side, l1, l2);
                        return;
                    }
                }
            }
        }

        // `cxx:1470-1475`.
        build_curves3d(&f1.0);
        build_curves3d(&f2.0);
        Self::update_vertex_tolerances(f1);
        Self::update_vertex_tolerances(f2);

        // `cxx:1477-1481`.
        let mut filler = PaveFiller::new();
        filler.set_arguments(&[f1.0.clone(), f2.0.clone()]);
        if filler.perform().is_err() {
            return;
        }

        // `cxx:1483-1487`.
        let mut true_edges: Vec<Edge> = Vec::new();
        if let Some(re) = ref_edge {
            Self::check_inters_ff(filler.ds(), re, &mut true_edges);
        }

        let pds = filler.ds();

        // `cxx:1496-1497`.
        l1.clear();
        l2.clear();

        // `cxx:1504-1591`.
        for a_ff in pds.interf_ff() {
            for a_bc in a_ff.curves() {
                for a_pb in a_bc.pave_blocks() {
                    let Some(shape) = pds.shape(a_pb.edge()) else {
                        continue;
                    };
                    let Some(an_edge) = Edge::wrap(shape.clone()) else {
                        continue;
                    };
                    // `cxx:1522-1526`.
                    if !true_edges.is_empty() && !true_edges.iter().any(|t| is_same(&t.0, &an_edge.0))
                    {
                        continue;
                    }

                    // `cxx:1529-1539`.
                    let (f, l) = BRepTool::edge_parameters(&an_edge);
                    let a_c3d = BRepTool::edge_curve(&an_edge);
                    let a_tol_edge = BRepTool::edge_tolerance(&an_edge);

                    // `cxx:1540-1559`: the pcurve on F1.
                    if !boptools_2d::has_curve_on_surface(&an_edge, f1) {
                        if let Some(pcurve) = a_bc.pcurve1() {
                            let mut a_c2d = pcurve.clone();
                            // `cxx:1543-1557`: with a 3D curve the pcurve is
                            // shifted by whole surface periods; the periodic
                            // and trimmed-curve overloads both use `[f, l]`.
                            if a_c3d.is_some() {
                                if let Ok(adjusted) =
                                    AlgoTools2D::adjust_pcurve_on_surf(&a_c2d, f1, f, l, a_tol_edge)
                                {
                                    a_c2d = adjusted;
                                }
                            }
                            update_edge_pcurve(&an_edge, f1, a_c2d, a_tol_edge);
                        }
                    }

                    // `cxx:1560-1579`: the pcurve on F2.
                    if !boptools_2d::has_curve_on_surface(&an_edge, f2) {
                        if let Some(pcurve) = a_bc.pcurve2() {
                            let mut a_c2d = pcurve.clone();
                            if a_c3d.is_some() {
                                if let Ok(adjusted) =
                                    AlgoTools2D::adjust_pcurve_on_surf(&a_c2d, f2, f, l, a_tol_edge)
                                {
                                    a_c2d = adjusted;
                                }
                            }
                            update_edge_pcurve(&an_edge, f2, a_c2d, a_tol_edge);
                        }
                    }

                    // `cxx:1580-1588`.
                    let Some((mut o1, mut o2)) = Self::orient_section(&an_edge, f1, f2) else {
                        continue;
                    };
                    if side == State::Out {
                        o1 = o1.reversed();
                        o2 = o2.reversed();
                    }
                    l1.push(an_edge.0.oriented(o1));
                    l2.push(an_edge.0.oriented(o2));
                }
            }
        }

        // `cxx:1593-1623`.
        let mut add_pcurve1 = true;
        let mut add_pcurve2 = true;
        let mut is_el1 = false;
        let mut is_el2 = false;
        if let Some(s) = &surf1 {
            let kind = surface_kind(s.as_ref());
            if kind == GeomAbsSurfaceType::Plane {
                add_pcurve1 = false;
            } else if is_elementary(kind) {
                is_el1 = true;
            }
        }
        if let Some(s) = &surf2 {
            let kind = surface_kind(s.as_ref());
            if kind == GeomAbsSurfaceType::Plane {
                add_pcurve2 = false;
            } else if is_elementary(kind) {
                is_el2 = true;
            }
        }

        // `cxx:1625-1702`: drop the section edges outside the reference edge's
        // range when the reference faces have closed surfaces.
        if l1.len() > 1 && (!is_el1 || !is_el2) {
            if let Some(rf1) = ref_face1 {
                if let (Some(rf2), Some(re)) = (ref_face2, ref_edge) {
                    let (v1, v2) = Self::edge_vertices(re);
                    // `if (!aV1.IsSame(aV2))`: only an open reference edge.
                    if !same_vertex(&v1, &v2) {
                        if let (Some(rs1), Some(rs2)) =
                            (BRepTool::face_surface(rf1), BRepTool::face_surface(rf2))
                        {
                            if rs1.is_u_closed()
                                || rs1.is_v_closed()
                                || rs2.is_u_closed()
                                || rs2.is_v_closed()
                            {
                                Self::remove_out_of_range_edges(l1, l2, re);
                            }
                        }
                    }
                }
            }
        }

        // `cxx:1704-1890`.
        if l1.len() > 1 && (!is_el1 || !is_el2) {
            let mut eseq: Vec<TopoShape> = Vec::new();

            if !true_edges.is_empty() {
                // `cxx:1709-1728`: `TrueEdges` is walked from its last to its
                // first element.
                let edges_for_concat: Vec<Edge> = true_edges.iter().rev().cloned().collect();
                match Self::assemble_edge(
                    pds,
                    f1,
                    f2,
                    add_pcurve1,
                    add_pcurve2,
                    &edges_for_concat,
                ) {
                    Some(assembled) => eseq.push(assembled.0),
                    None => eseq.extend(true_edges.iter().rev().map(|e| e.0.clone())),
                }
            } else {
                // `cxx:1730-1890`.
                let mut edges: Vec<Edge> = l1
                    .iter()
                    .filter_map(|s| Edge::wrap(s.clone()))
                    .collect();
                let mut wseq: Vec<Wire> = Vec::new();
                let builder = TopoBuilder::new();

                while !edges.is_empty() {
                    // `cxx:1739-1742`.
                    let mut an_edge = edges[0].clone();
                    let mut a_wire = builder.make_wire(std::slice::from_ref(&an_edge));
                    let mut candidates: Vec<usize> = Vec::new();
                    let mut matched: Option<usize> = None;
                    for (k, wire) in wseq.iter().enumerate() {
                        if Self::are_connex(wire, &a_wire) {
                            candidates.push(1);
                            matched = Some(k);
                            break;
                        }
                    }
                    let Some(matched) = matched else {
                        // `cxx:1754-1758`.
                        wseq.push(a_wire);
                        edges.remove(0);
                        continue;
                    };

                    // `cxx:1759-1771`.
                    let mut res_wire = wseq[matched].clone();
                    for j in 1..edges.len() {
                        an_edge = edges[j].clone();
                        a_wire = builder.make_wire(std::slice::from_ref(&an_edge));
                        if Self::are_connex(&res_wire, &a_wire) {
                            candidates.push(j + 1);
                        }
                    }

                    // `cxx:1772-1789`.
                    let mut minind = 1usize;
                    if candidates.len() > 1 {
                        let mut min_angle = f64::MAX;
                        for j in 1..=candidates.len() {
                            let edge = &edges[candidates[j - 1] - 1];
                            if let Some(angle) = Self::angle_wire_edge(&res_wire, edge) {
                                if angle < min_angle {
                                    min_angle = angle;
                                    minind = j;
                                }
                            }
                        }
                    }
                    let chosen = candidates[minind - 1] - 1;
                    builder.add_edge(&mut res_wire, &edges[chosen]);
                    wseq[matched] = res_wire;
                    edges.remove(chosen);
                }

                // `cxx:1792-1889`.
                for wire in &wseq {
                    let local = Self::edges_for_concat_of_wire(wire, pds);
                    match Self::assemble_edge(pds, f1, f2, add_pcurve1, add_pcurve2, &local) {
                        Some(assembled) => eseq.push(assembled.0),
                        None => eseq.extend(local.iter().rev().map(|e| e.0.clone())),
                    }
                }
            }

            // `cxx:1892-1929`.
            if eseq.len() < l1.len() {
                l1.clear();
                l2.clear();
                for a_shape in eseq {
                    let Some(an_edge) = Edge::wrap(a_shape) else {
                        continue;
                    };
                    same_parameter_edge_inplace(&an_edge, CONFUSION);
                    let edge_tol = BRepTool::edge_tolerance(&an_edge);
                    if edge_tol > 1.0e-2 {
                        continue;
                    }
                    if edge_tol >= 1.0e-4 {
                        Self::reconstruct_pcurves(&an_edge);
                        same_parameter_edge_inplace(&an_edge, CONFUSION);
                    }
                    let Some((mut o1, mut o2)) = Self::orient_section(&an_edge, f1, f2) else {
                        continue;
                    };
                    if side == State::Out {
                        o1 = o1.reversed();
                        o2 = o2.reversed();
                    }
                    l1.push(an_edge.0.oriented(o1));
                    l2.push(an_edge.0.oriented(o2));
                }
            }
        } else {
            // `cxx:1931-1941`.
            for shape in l1.iter() {
                if let Some(an_edge) = Edge::wrap(shape.clone()) {
                    same_parameter_edge_inplace(&an_edge, CONFUSION);
                }
            }
        }
    }

    /// `cxx:1632-1701`: keep only the section edge with the smallest angle
    /// between the reference edge's midpoint and the edge's own midpoint.
    fn remove_out_of_range_edges(l1: &mut Vec<TopoShape>, l2: &mut Vec<TopoShape>, ref_edge: &Edge) {
        let Some(ref_curve) = BRepTool::edge_curve_world(ref_edge) else {
            return;
        };
        let (ref_first, ref_last) = BRepTool::edge_parameters(ref_edge);
        let ref_pnt = ref_curve.d0(0.5 * (ref_first + ref_last));

        let mut min_angle = f64::MAX;
        let mut min_angle_edge: Option<Edge> = None;

        for shape in l1.iter() {
            let Some(an_edge) = Edge::wrap(shape.clone()) else {
                continue;
            };
            let Some(curve) = BRepTool::edge_curve_world(&an_edge) else {
                continue;
            };
            let (first, last) = BRepTool::edge_parameters(&an_edge);
            let mid_pnt = curve.d0(0.5 * (first + last));
            let ref_to_mid = GpVec::from_pnts(&ref_pnt, &mid_pnt);

            // `Extrema_ExtPC aProjector(aRefPnt, aBAcurve)` (`cxx:1653`).
            let (solutions, done) =
                extrema_ext_pc_range_done(curve.as_ref(), &ref_pnt, first, last);
            if !done {
                continue;
            }
            let mut imin = 0usize;
            let mut min_sq_dist = f64::MAX;
            for (index, solution) in solutions.iter().enumerate() {
                if solution.sq_dist < min_sq_dist {
                    min_sq_dist = solution.sq_dist;
                    imin = index + 1;
                }
            }
            if imin == 0 {
                continue;
            }
            let projection = solutions[imin - 1].point;
            let ref_to_proj = GpVec::from_pnts(&ref_pnt, &projection);
            let an_angle = ref_to_proj.angle(&ref_to_mid);
            if an_angle < min_angle {
                min_angle = an_angle;
                min_angle_edge = Some(an_edge);
            }
        }

        let Some(min_edge) = min_angle_edge else {
            return;
        };
        // `cxx:1683-1698`: the two lists lose the other entries in step.
        let mut kept1: Vec<TopoShape> = Vec::new();
        let mut kept2: Vec<TopoShape> = Vec::new();
        for (a, b) in l1.iter().zip(l2.iter()) {
            if is_same(a, &min_edge.0) {
                kept1.push(a.clone());
                kept2.push(b.clone());
            }
        }
        *l1 = kept1;
        *l2 = kept2;
    }

    /// `cxx:1795-1873`: the edge order of one wire, used as `EdgesForConcat`.
    ///
    /// A closed wire starts from an edge with a non-autonomous endpoint
    /// (`cxx:1799-1841`) and then walks the shared vertices (`cxx:1842-1866`);
    /// an open one is walked with `BRepTools_WireExplorer` (`cxx:1868-1873`).
    fn edges_for_concat_of_wire(wire: &Wire, pds: &crate::bopds::BopdsDS) -> Vec<Edge> {
        let mut local: Vec<Edge> = Vec::new();
        // `if (aWire.Closed())` reads the `TopoDS_TWire` flag; `BB.MakeWire`
        // plus `BB.Add` never raise it, so the sequence below is the OCCT
        // sequence as well (the `else` arm).
        if wire.is_closed() {
            let mut elist: Vec<Edge> = Vec::new();
            let mut start_found = false;
            let mut start_vertex: Option<Vertex> = None;
            let mut start_edge: Option<Edge> = None;

            for an_edge in edges_of_wire(wire) {
                if start_found {
                    elist.push(an_edge);
                    continue;
                }
                let (v1, v2) = Self::edge_vertices(&an_edge);
                if !is_autonom(v1.as_ref(), pds) {
                    start_vertex = v2;
                    start_edge = Some(an_edge);
                    start_found = true;
                } else if !is_autonom(v2.as_ref(), pds) {
                    start_vertex = v1;
                    start_edge = Some(an_edge);
                    start_found = true;
                } else {
                    elist.push(an_edge);
                }
            }

            if !start_found {
                let Some(first) = elist.first().cloned() else {
                    // Port guard: OCCT reads the first element of an empty list.
                    return Self::wire_explorer_edges(wire);
                };
                elist.remove(0);
                let (v1, _) = Self::edge_vertices(&first);
                start_edge = Some(first);
                start_vertex = v1;
            }

            let Some(start_edge) = start_edge else {
                return Self::wire_explorer_edges(wire);
            };
            local.push(start_edge);

            while !elist.is_empty() {
                let mut found = false;
                let mut index = 0usize;
                for (i, an_edge) in elist.iter().enumerate() {
                    let (v1, v2) = Self::edge_vertices(an_edge);
                    if same_vertex(&v1, &start_vertex) {
                        start_vertex = v2;
                        index = i;
                        found = true;
                        break;
                    } else if same_vertex(&v2, &start_vertex) {
                        start_vertex = v1;
                        index = i;
                        found = true;
                        break;
                    }
                }
                if !found {
                    // Port guard: the reference `while` never terminates here.
                    break;
                }
                local.push(elist.remove(index));
            }
        } else {
            local = Self::wire_explorer_edges(wire);
        }
        local
    }

    /// `BRepTools_WireExplorer` over the wire (`cxx:1869-1872`).
    fn wire_explorer_edges(wire: &Wire) -> Vec<Edge> {
        let mut out: Vec<Edge> = Vec::new();
        let mut explorer = WireExplorer::of_wire(wire);
        while explorer.more() {
            if let Some(edge) = explorer.current_edge() {
                out.push(edge);
            }
            explorer.next();
        }
        out
    }
}
