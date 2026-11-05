//! `BRepFill_OffsetWire::Perform` (`cxx:454-638`) and the `CheckBadEdges` /
//! `PerformCurve` / `QuasiFleche` helpers it uses (`cxx:2658-2920`).
//!
//! `Perform` is the driver: the `KPartCircle` shortcut, then `CheckBadEdges`
//! which asks the bisecting locus whether it sees every spine edge (an edge
//! whose curvature radius is below `Offset` is split into its sampled points),
//! and finally `PerformWithBiLo` - either on the locus computed by `Init` or on
//! a locus recomputed after the badly seen edges are subdivided.
//!
//! Deviations (handle / exception correspondence only, no missing branch):
//! * `Standard_Failure` out of the `try` of `cxx:456-616` and the
//!   `Standard_ConstructionError` of `cxx:634` answer the state `Perform`
//!   leaves after its `catch` (`myShape.Nullify(); myIsDone = false;`). The
//!   OCCT exceptions of that body are panics in the port, so [`perform`] runs
//!   `perform_body` under `catch_unwind` and restores exactly that state.
//! * `BRepLib_MakeEdge(V1, V2)` with a null vertex raises in OCCT; the port
//!   keeps that bad edge unsubstituted ([`make_edge_from_vertices`]).
//! * `NCollection_DataMap<...>::Bind` is `HashMap::insert`, and
//!   `NCollection_List` / `NCollection_Sequence` are `Vec`.

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt, GpVec, GpVec2d};
use occt_core::kernel::geomabs::CurveType;
use occt_core::precision::{CONFUSION, REAL_SMALL, RESOLUTION};
use occt_geom::curve::Curve;
use occt_geom2d::cl_props2d::GeomLPropCLProps2d;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_fill_offset_wire::{BRepFillOffsetWire, IndexedShapeMap};
use crate::brep_mat2d_bisecting_locus::BRepMat2dBisectingLocus;
use crate::brep_mat2d_explorer::BRepMat2dExplorer;
use crate::brep_mat2d_link_topo_bilo::BRepMat2dLinkTopoBilo;
use crate::brep_tool::BRepTool;
use crate::brep_tools_substitution::BRepToolsSubstitution;
use crate::builder::TopoBuilder;
use crate::mat::MatSide;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, shape_map_key, wires_of, wires_of_face};

impl BRepFillOffsetWire {
    /// `Perform(Offset, Alt)` (`cxx:454-638`).
    ///
    /// `cxx:456-616` is a `try` / `catch (Standard_Failure const&)` whose catch
    /// leaves `myShape` null and `myIsDone` false. The port raises the OCCT
    /// exceptions of that body as panics, so the body runs under
    /// `catch_unwind` and the catch state is restored here.
    pub fn perform(&mut self, offset: f64, alt: f64) {
        let prev_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {}));
        let land = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            self.perform_body(offset, alt)
        }));
        std::panic::set_hook(prev_hook);
        if land.is_err() {
            // `cxx:606-617`.
            self.my_shape = None;
            self.my_is_done = false;
        }
    }

    /// The `try` block of `Perform` (`cxx:458-605`).
    fn perform_body(&mut self, offset: f64, alt: f64) {
        // `cxx:460`.
        self.my_call_gen = false;

        // `cxx:462-465`: the single closed circle shortcut.
        if self.kpart_circle(offset, alt) {
            return;
        }

        // `cxx:467-469`.
        let Some(work_spine) = self.my_work_spine.clone() else {
            // `CheckBadEdges` on a null spine raises `Standard_NullObject`.
            self.my_shape = None;
            self.my_is_done = false;
            return;
        };
        let mut bad_edges: Vec<TopoShape> = Vec::new();
        if !check_bad_edges(&work_spine, offset, &self.my_bilo, &self.my_link, &mut bad_edges) {
            // `LProp_NotDefined` out of `CheckBadEdges` is a `Standard_Failure`
            // of the `try` of `cxx:456-616`.
            self.my_shape = None;
            self.my_is_done = false;
            return;
        }

        if bad_edges.is_empty() {
            // `cxx:600-602`.
            let bilo = std::mem::replace(&mut self.my_bilo, BRepMat2dBisectingLocus::new());
            let link = std::mem::replace(&mut self.my_link, BRepMat2dLinkTopoBilo::new());
            let join = self.my_join_type;
            self.perform_with_bilo(&work_spine, offset, &bilo, &link, join, alt);
            self.my_bilo = bilo;
            self.my_link = link;
        } else {
            // `cxx:471-599`.
            self.substitute_bad_edges(&work_spine, &bad_edges, offset, alt);
        }

        // `cxx:619-636`: a closed result must have closed wires.
        if self.my_is_done && !self.my_is_open_result {
            let Some(shape) = self.my_shape.clone() else {
                return;
            };
            for wire in wires_of(&shape) {
                if !wire.0.is_closed() {
                    self.my_shape = None;
                    self.my_is_done = false;
                    return;
                }
            }
        }
    }

    /// `cxx:471-599`: rebuild the badly seen spine edges as chains of
    /// `PerformCurve` samples, substitute them in `myWorkSpine`, and recompute
    /// the locus on the substituted face.
    fn substitute_bad_edges(
        &mut self,
        work_spine: &Face,
        bad_edges: &[TopoShape],
        offset: f64,
        alt: f64,
    ) {
        let builder = TopoBuilder::new();
        let old_work_spine = work_spine.clone();
        let a_defl = 0.01 * offset.abs();
        let mut subst = BRepToolsSubstitution::new();

        // `cxx:484-548`.
        for an_e in bad_edges {
            let e = Edge(an_e.clone());
            let (vf, vl) = edge_vertices(&e);
            let (fpar, lpar) = BRepTool::edge_parameters(&e);
            let Some(curve) = BRepTool::edge_curve(&e) else {
                // `GeomAdaptor_Curve` over a null curve raises.
                continue;
            };

            // `cxx:496`: `PerformCurve(..., aDefl, f, l, Confusion(), 2)`.
            let mut parameters: Vec<f64> = Vec::new();
            let mut points: Vec<GpPnt> = Vec::new();
            perform_curve(
                &mut parameters,
                &mut points,
                &curve,
                a_defl,
                fpar,
                lpar,
                CONFUSION,
                2,
            );

            // `cxx:499-524`.
            let mut a_l: Vec<TopoShape> = Vec::new();
            let mut failed = false;
            if points.len() > 2 {
                let mut fv = vf.clone();
                for np in 2..points.len() {
                    let lp = points[np - 1];
                    let lv = builder.make_vertex(lp, CONFUSION);
                    match make_edge_from_vertices(&builder, fv.as_ref(), Some(&lv)) {
                        Some(ne) => a_l.push(ne.0),
                        None => {
                            failed = true;
                            break;
                        }
                    }
                    fv = Some(lv);
                }
                if !failed {
                    match make_edge_from_vertices(&builder, fv.as_ref(), vl.as_ref()) {
                        Some(ne) => a_l.push(ne.0),
                        None => failed = true,
                    }
                }
            } else {
                // `cxx:518-524`.
                match make_edge_from_vertices(&builder, vf.as_ref(), vl.as_ref()) {
                    Some(ne) => a_l.push(ne.0),
                    None => failed = true,
                }
            }
            if failed {
                // OCCT raises out of the `try`, leaving the whole `Perform`
                // with no result.
                self.my_shape = None;
                self.my_is_done = false;
                return;
            }

            // `cxx:526-547`: the new edges inherit the old edge's spine image.
            if let Some((_key_shape, value)) = self.my_map_spine.get(&shape_map_key(&e.0)).cloned() {
                for ne in &a_l {
                    self.my_map_spine.insert(shape_map_key(ne), (ne.clone(), value.clone()));
                    let (nv1, nv2) = edge_vertices(&Edge(ne.clone()));
                    for v in [nv1, nv2].into_iter().flatten() {
                        self.my_map_spine
                            .entry(shape_map_key(&v.0))
                            .or_insert_with(|| (v.0.clone(), value.clone()));
                    }
                }
                self.my_map_spine.remove(&shape_map_key(&e.0));
            }

            // `cxx:548`.
            subst.substitute(an_e, &a_l);
        }

        // `cxx:551-568`: copy the wires, then re-register them on the face.
        let mut wwmap: Vec<(TopoShape, Vec<TopoShape>)> = Vec::new();
        for wire in wires_of_face(work_spine) {
            subst.build(&wire.0);
            if subst.is_copied(&wire.0) {
                let new_wire = subst.copy(&wire.0).expect("IsCopied")[0].clone();
                new_wire.set_closed(wire.0.is_closed());
                wwmap.push((wire.0.clone(), vec![new_wire]));
            }
        }
        subst.clear();
        for (key, value) in &wwmap {
            subst.substitute(key, value);
        }

        // `cxx:575-599`.
        subst.build(&work_spine.0);
        if !subst.is_copied(&work_spine.0) {
            let bilo = std::mem::replace(&mut self.my_bilo, BRepMat2dBisectingLocus::new());
            let link = std::mem::replace(&mut self.my_link, BRepMat2dLinkTopoBilo::new());
            let join = self.my_join_type;
            self.perform_with_bilo(work_spine, offset, &bilo, &link, join, alt);
            self.my_bilo = bilo;
            self.my_link = link;
            return;
        }

        let new_work = Face(subst.copy(&work_spine.0).expect("IsCopied")[0].clone());
        let mut exp = BRepMat2dExplorer::new();
        exp.perform(&new_work);
        let mut bilo = BRepMat2dBisectingLocus::new();
        bilo.compute(
            &mut exp,
            1,
            MatSide::Left,
            self.my_join_type,
            self.my_is_open_result,
        );
        if !bilo.is_done() {
            self.my_shape = None;
            self.my_is_done = false;
            return;
        }
        let mut link = BRepMat2dLinkTopoBilo::new();
        link.perform(&exp, &bilo);
        let join = self.my_join_type;
        self.perform_with_bilo(&new_work, offset, &bilo, &link, join, alt);
        self.my_work_spine = Some(old_work_spine);
    }
}

/// `BRepLib_MakeEdge(V1, V2)` (`BRepLib_MakeEdge.cxx:217-230`): the straight
/// line through the two vertex points with those vertices as ends.
///
/// `None` where OCCT raises: a null vertex (`Standard_NullObject`) or two
/// coincident points (`gp_Dir` on a null vector).
fn make_edge_from_vertices(
    builder: &TopoBuilder,
    v1: Option<&Vertex>,
    v2: Option<&Vertex>,
) -> Option<Edge> {
    let (v1, v2) = (v1?, v2?);
    let p1 = BRepTool::vertex_point(v1);
    let p2 = BRepTool::vertex_point(v2);
    if p1.distance(&p2) <= REAL_SMALL {
        return None;
    }
    Some(builder.make_edge_segment_with_vertices(&p1, &p2, v1, v2))
}

/// `CheckBadEdges(Spine, Offset, Locus, Link, BadEdges)` (`cxx:2658-2778`).
///
/// A spine edge is badly seen when the bisecting locus does not see it: its
/// p-curve is not a line or a circle and its curvature at an end is not below
/// `1 / Offset` on the offset side, or its tangent is not defined there.
///
/// `false` is the `LProp_NotDefined` arm of `CentreOfCurvature`
/// (`LProp_CurveUtils.hxx:461-462`), which OCCT raises out of `Perform`'s
/// `try`; the caller answers the state its `catch` leaves.
pub(crate) fn check_bad_edges(
    spine: &Face,
    offset: f64,
    locus: &BRepMat2dBisectingLocus,
    link: &BRepMat2dLinkTopoBilo,
    bad_edges: &mut Vec<TopoShape>,
) -> bool {
    // `cxx:2665-2667`.
    let f = Face(spine.0.oriented(Orientation::Forward));
    let eps = CONFUSION;
    let lim_curv = 1.0 / offset;
    let mut a_map: Vec<crate::topo_tools::ShapeMapKey> = Vec::new();

    // `cxx:2671-2775`.
    for ic in 1..=locus.number_of_contours() {
        for ie in 1..=locus.number_of_elts(ic) {
            let se = link.generating_shape(locus.basic_elt(ic, ie));
            if se.shape_type() != ShapeType::Edge {
                continue;
            }
            let key = shape_map_key(&se);
            if a_map.contains(&key) {
                continue;
            }
            let e = Edge(se.clone());
            let Some((g2d, fpar, lpar)) = curve_on_surface_range(&e, &f) else {
                // `BRep_Tool::CurveOnSurface` on a null p-curve raises.
                continue;
            };

            // `cxx:2691-2695`.
            let ac = occt_geom2d::adaptor_curve::Geom2dAdaptorCurve::with_range(
                g2d.clone(),
                fpar,
                lpar,
            );
            let Some(curve_type) = ac.get_type() else {
                continue;
            };
            if curve_type == CurveType::Line || curve_type == CurveType::Circle {
                continue;
            }
            // `cxx:2697-2702`.
            let reverse = e.0.orientation() == Orientation::Forward;

            let mut bad = false;
            for u in [fpar, lpar] {
                // `cxx:2706-2741`.
                let mut props = GeomLPropCLProps2d::new(g2d.as_ref(), u, 2, eps);
                if !props.is_tangent_defined() {
                    bad = true;
                    break;
                }
                let p = props.value();
                let crv = props.curvature();
                if crv >= eps {
                    // `cxx:2722-2739`: `(N.Y, -N.X)`, reversed on a FORWARD
                    // edge, then compared with the centre-of-curvature
                    // direction.
                    let t = props.tangent();
                    let mut n = GpDir2d::new(t.y, -t.x).expect("gp_Dir2d(N.Y, -N.X)");
                    if reverse {
                        n = n.reversed();
                    }
                    // `CentreOfCurvature` raises `LProp_NotDefined` when
                    // `|Curvature()| <= eps`, which is the `crv == eps` corner
                    // of the guard above.
                    if crv <= eps {
                        return false;
                    }
                    let pc = props.centre_of_curvature();
                    let dir = GpVec2d::new(pc.x() - p.x(), pc.y() - p.y());
                    if GpVec2d::new(n.x, n.y).dot(&dir) > 0.0 && lim_curv <= crv + eps {
                        bad = true;
                        break;
                    }
                }
            }
            if bad {
                // `cxx:2711-2714` / `:2718-2722`: `BadEdges.Append(SE);
                // aMap.Add(SE); continue;`.
                bad_edges.push(se.clone());
                a_map.push(key);
            }
        }
    }
    true
}

/// `PerformCurve(Parameters, Points, C, Deflection, U1, U2, EPSILON, Nbmin)`
/// (`cxx:2781-2824`): append `U1`, then the samples up to `U2`.
fn perform_curve(
    parameters: &mut Vec<f64>,
    points: &mut Vec<GpPnt>,
    c: &Arc<dyn Curve>,
    deflection: f64,
    u1: f64,
    u2: f64,
    epsilon: f64,
    nbmin: i32,
) -> bool {
    let uu1 = u1.min(u2);
    let uu2 = u1.max(u2);

    // `cxx:2793-2798`.
    let (pdeb, ddeb) = c.d1(uu1);
    parameters.push(uu1);
    points.push(pdeb);

    // `cxx:2800-2821`.
    let (pfin, dfin) = c.d1(uu2);
    let a_delta = uu2 - uu1;
    let a_dist = pdeb.distance(&pfin);
    if (a_delta / a_dist) > 5.0e-14 {
        quasi_fleche(
            c,
            deflection * deflection,
            uu1,
            &pdeb,
            &ddeb,
            uu2,
            &pfin,
            &dfin,
            nbmin,
            epsilon * epsilon,
            parameters,
            points,
        );
    }
    true
}

/// `QuasiFleche(C, Deflection2, Udeb, Pdeb, Vdeb, Ufin, Pfin, Vfin, Nbmin, Eps,
/// Parameters, Points)` (`cxx:2827-2923`).
#[allow(clippy::too_many_arguments)]
fn quasi_fleche(
    c: &Arc<dyn Curve>,
    deflection2: f64,
    udeb: f64,
    pdeb: &GpPnt,
    vdeb: &GpVec,
    ufin: f64,
    pfin: &GpPnt,
    vfin: &GpVec,
    nbmin: i32,
    eps: f64,
    parameters: &mut Vec<f64>,
    points: &mut Vec<GpPnt>,
) {
    // `cxx:2841-2856`.
    let pts_length = points.len();
    let mut u_delta = ufin - udeb;
    let (p_delta, v_delta) = if nbmin > 2 {
        u_delta /= f64::from(nbmin - 1);
        c.d1(udeb + u_delta)
    } else {
        (*pfin, *vfin)
    };

    // `cxx:2858-2885`: the arrow, `IntWalk_IWalking::TestDeflection` style,
    // falling back to the midpoint deviation.
    let norme = GpVec::from_pnts(pdeb, &p_delta).square_magnitude();
    let mut the_fleche = 0.0;
    let mut fleche_ok = false;
    if norme > eps {
        let n1 = vdeb.square_magnitude();
        let n2 = v_delta.square_magnitude();
        if n1 > eps && n2 > eps {
            let va = vdeb.normalized();
            let vb = v_delta.normalized();
            let a = va.xyz();
            let b = vb.xyz();
            let normediff =
                (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
            if normediff > eps {
                the_fleche = normediff * norme / 64.0;
                fleche_ok = true;
            }
        }
    }
    if !fleche_ok {
        let pmid = GpPnt::new(
            (pdeb.x() + p_delta.x()) / 2.0,
            (pdeb.y() + p_delta.y()) / 2.0,
            (pdeb.z() + p_delta.z()) / 2.0,
        );
        let pverif = c.d0(udeb + u_delta / 2.0);
        the_fleche = pmid.square_distance(&pverif);
    }

    // `cxx:2887-2899`.
    if the_fleche < deflection2 {
        parameters.push(udeb + u_delta);
        points.push(p_delta);
    } else {
        quasi_fleche(
            c,
            deflection2,
            udeb,
            pdeb,
            vdeb,
            udeb + u_delta,
            &p_delta,
            &v_delta,
            3,
            eps,
            parameters,
            points,
        );
    }

    // `cxx:2905-2921`.
    if nbmin > 2 {
        let nb = nbmin - (points.len() - pts_length) as i32;
        quasi_fleche(
            c,
            deflection2,
            udeb + u_delta,
            &p_delta,
            &v_delta,
            ufin,
            pfin,
            vfin,
            nb,
            eps,
            parameters,
            points,
        );
    }
}

/// `Compute(Spine, aShape, Map, Alt)` (`cxx:643-676`): the zero-offset case, a
/// compound of the spine wires translated by `Alt` along `Z`.
pub(crate) fn compute_offset_shape(
    spine: &Face,
    map: &mut IndexedShapeMap,
    alt: f64,
) -> TopoShape {
    let builder = TopoBuilder::new();
    let mut shape = builder.make_shape(ShapeType::Compound);

    // `cxx:650-654`.
    let mut a = alt;
    if spine.0.orientation() == Orientation::Reversed {
        a = -alt;
    }
    let mut trsf = occt_core::gp::GpTrsf::identity();
    trsf.set_translation_vec(&GpVec::new(0.0, 0.0, a));
    let location = occt_core::toploc::TopLocLocation::composed(
        &trsf,
        occt_core::toploc::TopLocLocation::identity(),
    );

    // `cxx:656-675`.
    for cur_w in wires_of_face(spine) {
        let new_w = cur_w.0.located(&location);
        builder.add(&mut shape, &new_w);
        let moved_edges = crate::topo_tools::edges_of(&new_w);
        for (e1, e2) in crate::topo_tools::edges_of_wire(&cur_w).iter().zip(moved_edges) {
            map.add(&e1.0, vec![e2.0]);
        }
    }
    shape
}

/// `fabs(Alt) > gp::Resolution()`, the test `PerformWithBiLo` uses to decide
/// whether the reference plane is translated (`cxx:744`).
pub(crate) fn is_altitude_significant(alt: f64) -> bool {
    alt.abs() > RESOLUTION
}
