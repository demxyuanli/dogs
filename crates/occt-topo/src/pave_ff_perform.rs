//! `BOPAlgo_PaveFiller::PerformFF` — full OCCT 8.0.0 translation.
//!
//! Source: `BOPAlgo_PaveFiller_6.cxx:285-622`.
//!
//! Sequence:
//! 1. Update FaceInfo On/In for every FACE/FACE iterator pair and every
//!    already-touched face (`HasReference` / FaceInfo already allocated).
//! 2. If the iterator is empty, return.
//! 3. Glue != Off: append empty `InterfFF` records and skip FaceFace.
//! 4. Glue Off: CheckPlanes for plane/plane; build an EE-new-vertex map;
//!    for non-planar pairs, shift a closed (seam) face onto the exact EE
//!    intersection when the EE vertex sits off the edges; run FaceFace with
//!    `GetEFPnts` start points; `ToleranceFF` plus the shift as `aTolFF`.
//! 5. On FaceFace failure: empty `InterfFF` plus an intersection-failed
//!    warning (`AddIntersectionFailedWarning`).
//! 6. On success: `CheckCurve` each line (`Bnd_Box::IsThin(3*Confusion)`),
//!    expand the box by `aTolFF + max vertex tolerance`, store points.
//!
//! Section edges are still born later in `MakeBlocks` / `PostTreatFF`.
//! `BOPTools_Parallel::Perform` is sequential in this port.

use std::collections::{HashMap, HashSet};

use occt_core::gp::GpVec;
use occt_core::precision::CONFUSION;
use occt_geom::Surface;

use crate::abs::ShapeType;
use crate::bopds_ff::BopdsInterfFf;
use crate::bopds_ff::BopdsCurve;
use crate::brep_surface::{classify_surface, face_is_planar, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::int_face_face::FaceFace;
use crate::int_tools_curve_box;
use crate::int_tools_full::IntToolsContext;
use crate::iterator::ShapeIterator;
use crate::pave_ff_misc::check_planes;
use crate::pave_ff_paves::get_ef_pnts;
use crate::pave_filler::{GlueEnum, PaveFiller};
use crate::pave_intersect::collect_pairs;
use crate::shape::{Edge, Face, Vertex};
use crate::transform::translated;

/// Analytic surface kinds that skip the `5.e-6` floor in [`tolerance_ff`].
fn is_analytic_ff(kind: SurfaceKind) -> bool {
    matches!(
        kind,
        SurfaceKind::Plane
            | SurfaceKind::Cylinder
            | SurfaceKind::Cone
            | SurfaceKind::Sphere
            | SurfaceKind::Torus
    )
}

/// `ToleranceFF(BRepAdaptor_Surface, BRepAdaptor_Surface)` at `_6.cxx:3922`.
///
/// `max(Tol1, Tol2)`, then floored at `5.e-6` when either surface is not
/// plane/cylinder/cone/sphere/torus.
pub fn tolerance_ff(face1: &Face, face2: &Face) -> f64 {
    let a_tol1 = BRepTool::face_tolerance(face1);
    let a_tol2 = BRepTool::face_tolerance(face2);
    let mut a_tol_ff = a_tol1.max(a_tol2);
    let k1 = BRepTool::face_surface(face1)
        .map(|s| classify_surface(s.as_ref()))
        .unwrap_or(SurfaceKind::Other);
    let k2 = BRepTool::face_surface(face2)
        .map(|s| classify_surface(s.as_ref()))
        .unwrap_or(SurfaceKind::Other);
    if !is_analytic_ff(k1) || !is_analytic_ff(k2) {
        a_tol_ff = a_tol_ff.max(5.0e-6);
    }
    a_tol_ff
}

/// `IsPlaneFF` at `_6.cxx:84` — Geom_Plane, or plane under offset / trim.
pub fn is_plane_ff(surf: &dyn Surface) -> bool {
    classify_surface(surf) == SurfaceKind::Plane
}

/// `IsClosedFF` at `_6.cxx:106`.
///
/// For a non-plane surface OCCT walks the edge representations looking for
/// `IsCurveOnClosedSurface`. The port treats a seam as an edge that appears
/// twice on the face, which is the topological counterpart of a closed
/// pcurve pair. For a plane OCCT only checks triangulation; without a
/// triangulation the plane path returns false.
pub fn is_closed_ff(ds: &crate::bopds::BopdsDS, n_f: usize, n_e: usize, is_plane: bool) -> bool {
    if is_plane {
        return false;
    }
    let Some(e_shape) = ds.shape(n_e) else {
        return false;
    };
    if BRepTool::is_closed_edge(&Edge(e_shape.clone())) {
        return true;
    }
    edge_count_on_face(ds, n_f, n_e) >= 2
}

fn edge_count_on_face(ds: &crate::bopds::BopdsDS, n_f: usize, n_e: usize) -> usize {
    let Some(f_shape) = ds.shape(n_f) else {
        return 0;
    };
    let mut count = 0usize;
    for child in ShapeIterator::of_shape(f_shape) {
        if child.shape_type() == ShapeType::Edge {
            if ds.index(&child) == Some(n_e) {
                count += 1;
            }
            continue;
        }
        if child.shape_type() != ShapeType::Wire {
            continue;
        }
        for edge in ShapeIterator::of_shape(&child) {
            if edge.shape_type() != ShapeType::Edge {
                continue;
            }
            if ds.index(&edge) == Some(n_e) {
                count += 1;
            }
        }
    }
    count
}

/// Direct children of a face that are edges, including those nested in wires.
pub fn face_edges(ds: &crate::bopds::BopdsDS, n_f: usize) -> Vec<usize> {
    let Some(f_shape) = ds.shape(n_f) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for child in ShapeIterator::of_shape(f_shape) {
        if child.shape_type() == ShapeType::Edge {
            if let Some(n) = ds.index(&child) {
                out.push(n);
            }
            continue;
        }
        if child.shape_type() != ShapeType::Wire {
            continue;
        }
        for edge in ShapeIterator::of_shape(&child) {
            if edge.shape_type() != ShapeType::Edge {
                continue;
            }
            if let Some(n) = ds.index(&edge) {
                out.push(n);
            }
        }
    }
    out
}

/// Max vertex tolerance of a face (`BRep_Tool::MaxTolerance(face, VERTEX)`).
pub fn max_vertex_tolerance_of_face(ds: &crate::bopds::BopdsDS, n_f: usize) -> f64 {
    let Some(si) = ds.shape_info(n_f) else {
        return 0.0;
    };
    let mut m = 0.0f64;
    for &n in si.sub_shapes() {
        let Some(s) = ds.shape(n) else {
            continue;
        };
        if s.shape_type() != ShapeType::Vertex {
            continue;
        }
        m = m.max(BRepTool::vertex_tolerance(&Vertex(s.clone())));
    }
    m
}

/// EE new-vertex map used by the closed-edge shift (`_6.cxx:335-360`).
fn build_ee_new_vertex_map(ds: &crate::bopds::BopdsDS) -> HashMap<(usize, usize), Vec<usize>> {
    let mut map: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for it in ds.interf_ee() {
        let Some(n_vn) = it.get_index_new() else {
            continue;
        };
        let (n_e1, n_e2) = it.indices();
        let key = (n_e1.min(n_e2), n_e1.max(n_e2));
        map.entry(key).or_default().push(n_vn);
    }
    map
}

fn ee_vertices<'a>(
    map: &'a HashMap<(usize, usize), Vec<usize>>,
    n_e1: usize,
    n_e2: usize,
) -> Option<&'a [usize]> {
    let key = (n_e1.min(n_e2), n_e1.max(n_e2));
    map.get(&key).map(|v| v.as_slice())
}

/// Closed-edge shift of `_6.cxx:399-486`.
///
/// When a pair of edges (at least one closed/seam) produced an EE vertex that
/// is farther from the two edges than the vertex tolerance, translate the
/// closed face so the edges meet at the exact intersection. The shift
/// distance becomes a floor on `aTolFF`.
fn shift_closed_faces(
    ds: &crate::bopds::BopdsDS,
    ctx: &IntToolsContext,
    n_f1: usize,
    n_f2: usize,
    face1: &Face,
    face2: &Face,
    ee_map: &HashMap<(usize, usize), Vec<usize>>,
) -> (Face, Face, f64) {
    let mut shifted1 = face1.clone();
    let mut shifted2 = face2.clone();
    let mut shift_value = 0.0;

    let plane1 = BRepTool::face_surface(face1)
        .map(|s| is_plane_ff(s.as_ref()))
        .unwrap_or(false);
    let plane2 = BRepTool::face_surface(face2)
        .map(|s| is_plane_ff(s.as_ref()))
        .unwrap_or(false);
    if plane1 && plane2 {
        return (shifted1, shifted2, shift_value);
    }

    let edges1 = face_edges(ds, n_f1);
    let edges2 = face_edges(ds, n_f2);
    let mut found = false;
    for &an_edge_index1 in &edges1 {
        if found {
            break;
        }
        let Some(e1s) = ds.shape(an_edge_index1).cloned() else {
            continue;
        };
        let edge1 = Edge(e1s);
        let closed1 = is_closed_ff(ds, n_f1, an_edge_index1, plane1);
        for &an_edge_index2 in &edges2 {
            if found {
                break;
            }
            let Some(e2s) = ds.shape(an_edge_index2).cloned() else {
                continue;
            };
            let edge2 = Edge(e2s);
            let closed2 = is_closed_ff(ds, n_f2, an_edge_index2, plane2);
            if !closed1 && !closed2 {
                continue;
            }
            let Some(verts) = ee_vertices(ee_map, an_edge_index1, an_edge_index2) else {
                continue;
            };
            for &a_vertex_index in verts {
                let Some(vs) = ds.shape(a_vertex_index).cloned() else {
                    continue;
                };
                let vertex = Vertex(vs);
                let a_vertex_point = BRepTool::vertex_point(&vertex);
                let t1 = ctx.project_point_on_edge(&edge1, &a_vertex_point);
                let t2 = ctx.project_point_on_edge(&edge2, &a_vertex_point);
                if t1.is_none() && t2.is_none() {
                    continue;
                }
                let a_p1 = if let (Some(t), Some(c)) = (t1, BRepTool::edge_curve(&edge1)) {
                    c.d0(t)
                } else {
                    a_vertex_point
                };
                let a_p2 = if let (Some(t), Some(c)) = (t2, BRepTool::edge_curve(&edge2)) {
                    c.d0(t)
                } else {
                    a_vertex_point
                };
                let a_shift_dist = a_p1.distance(&a_p2);
                if a_shift_dist > BRepTool::vertex_tolerance(&vertex) {
                    let vec = if closed1 {
                        GpVec::from_pnts(&a_p1, &a_p2)
                    } else {
                        GpVec::from_pnts(&a_p2, &a_p1)
                    };
                    if closed1 {
                        shifted1 = Face(translated(&shifted1.0, &vec));
                    } else {
                        shifted2 = Face(translated(&shifted2.0, &vec));
                    }
                    shift_value = a_shift_dist;
                    found = true;
                    break;
                }
            }
        }
    }
    (shifted1, shifted2, shift_value)
}

fn add_intersection_failed_warning(f: &mut PaveFiller, n_f1: usize, n_f2: usize) {
    f.add_warning(format!(
        "The intersection of the pair of faces {n_f1} and {n_f2} has failed"
    ));
}

struct FaceFaceJob {
    n_f1: usize,
    n_f2: usize,
    face1: Face,
    face2: Face,
    orig1: Face,
    orig2: Face,
    tol_ff: f64,
    starts: Vec<(f64, f64, f64, f64)>,
}

fn plane_angle(f1: &Face, f2: &Face) -> Option<f64> {
    let p1 = crate::brep_surface::face_plane(f1)?;
    let p2 = crate::brep_surface::face_plane(f2)?;
    Some(p1.axis().direction().angle(p2.axis().direction()))
}

/// `BOPAlgo_PaveFiller::PerformFF` (`_6.cxx:285`).
pub fn perform_ff(f: &mut PaveFiller) -> Result<(), String> {
    let pairs = collect_pairs(f.ds(), ShapeType::Face, ShapeType::Face);

    let mut fence: HashSet<usize> = HashSet::new();
    for &(n_f1, n_f2) in &pairs {
        fence.insert(n_f1);
        fence.insert(n_f2);
    }
    let n_src = f.ds().nb_source_shapes();
    for i in 0..n_src {
        let is_face = f
            .ds()
            .shape_info(i)
            .map(|s| s.shape_type() == ShapeType::Face)
            .unwrap_or(false);
        if is_face && f.ds().face_info(i).is_some() {
            fence.insert(i);
        }
    }
    f.ds_mut().update_face_info_on_faces(&fence);
    f.ds_mut().update_face_info_in_faces(&fence);

    if pairs.is_empty() {
        return Ok(());
    }

    let glue = f.glue();
    if glue != GlueEnum::None {
        for (n_f1, n_f2) in pairs {
            let mut rec = BopdsInterfFf::new(n_f1, n_f2);
            rec.set_tangent_faces(false);
            rec.init(0, 0);
            f.ds_mut().append_interf_ff(rec);
        }
        return Ok(());
    }

    let ee_map = build_ee_new_vertex_map(f.ds());
    let fuzzy = f.fuzzy_value();
    let ctx_tools = IntToolsContext::new();

    let mut jobs: Vec<FaceFaceJob> = Vec::new();
    for (n_f1, n_f2) in pairs {
        let Some(f1_shape) = f.ds().shape(n_f1).cloned() else {
            continue;
        };
        let Some(f2_shape) = f.ds().shape(n_f2).cloned() else {
            continue;
        };
        let face1 = Face(f1_shape);
        let face2 = Face(f2_shape);

        if face_is_planar(&face1) && face_is_planar(&face2) && !check_planes(f.ds(), n_f1, n_f2) {
            let mut rec = BopdsInterfFf::new(n_f1, n_f2);
            rec.init(0, 0);
            f.ds_mut().append_interf_ff(rec);
            continue;
        }

        let (shifted1, shifted2, shift_value) =
            shift_closed_faces(f.ds(), &ctx_tools, n_f1, n_f2, &face1, &face2, &ee_map);
        let a_tol_ff = shift_value.max(tolerance_ff(&face1, &face2));
        let starts: Vec<(f64, f64, f64, f64)> = get_ef_pnts(f, n_f1, n_f2)
            .into_iter()
            .map(|p| (p.u1, p.v1, p.u2, p.v2))
            .collect();
        jobs.push(FaceFaceJob {
            n_f1,
            n_f2,
            orig1: face1,
            orig2: face2,
            face1: shifted1,
            face2: shifted2,
            tol_ff: a_tol_ff,
            starts,
        });
    }

    for job in jobs {
        let n_f1 = job.n_f1;
        let n_f2 = job.n_f2;
        let mut ff = FaceFace::new();
        ff.set_face1(job.face1.clone());
        ff.set_face2(job.face2.clone());
        ff.set_tolerance(job.tol_ff.max(fuzzy));
        if !job.starts.is_empty() {
            ff.set_list(job.starts.clone());
        }
        let failed = match ff.perform() {
            Err(_) => true,
            Ok(()) => !ff.is_done(),
        };
        if failed {
            let mut rec = BopdsInterfFf::new(n_f1, n_f2);
            rec.init(0, 0);
            f.ds_mut().append_interf_ff(rec);
            add_intersection_failed_warning(f, n_f1, n_f2);
            continue;
        }

        let b_tangent = ff.tangent_faces();
        let a_tol_ff = job.tol_ff;
        let res = ff.result();
        let a_nb_curves = res.nb_curves();
        let a_nb_points = 0usize;

        if a_nb_curves > 0 || a_nb_points > 0 {
            f.ds_mut().add_interf(n_f1, n_f2);
        }

        let mut rec = BopdsInterfFf::new(n_f1, n_f2);
        rec.set_tangent_faces(b_tangent);
        rec.init(a_nb_curves, a_nb_points);

        let mut a_box_expand = a_tol_ff;
        if a_nb_curves > 0 {
            let max_v = max_vertex_tolerance_of_face(f.ds(), n_f1)
                .max(max_vertex_tolerance_of_face(f.ds(), n_f2));
            a_box_expand += max_v;
        }

        let tang_default = int_tools_curve_box::curve_tangential_tolerance(
            face_is_planar(&job.orig1),
            face_is_planar(&job.orig2),
            BRepTool::face_tolerance(&job.orig1),
            BRepTool::face_tolerance(&job.orig2),
            plane_angle(&job.orig1, &job.orig2),
        );

        for c in res.curves() {
            let r = c.range;
            let mut box_ = occt_core::bnd::BndBox::new();
            let b_valid = int_tools_curve_box::check_curve(
                Some(c.curve.as_ref()),
                r.first,
                r.last,
                a_tol_ff.max(CONFUSION),
                tang_default,
                &mut box_,
            );
            if !b_valid {
                continue;
            }
            int_tools_curve_box::enlarge_box(&mut box_, a_box_expand);
            let mut nc = BopdsCurve::new();
            nc.set_curve(c.curve.clone());
            nc.set_pcurves(c.pcurve1.clone(), c.pcurve2.clone());
            nc.set_range(r.first, r.last);
            nc.set_tolerance(a_tol_ff);
            nc.set_tangential_tolerance(tang_default);
            nc.set_box(box_);
            nc.init_pave_block1();
            rec.change_curves().push(nc);
        }

        f.ds_mut().append_interf_ff(rec);
    }
    Ok(())
}
