use super::prelude::*;
use super::{EdgeInfo, SmartMap};
use super::angles::*;

/// `BOPAlgo_WireSplitter::RefineAngles` (`BOPAlgo_WireSplitter_1.cxx:905-918`,
/// `:925-1029`): for every vertex, when exactly **two** boundary edges meet,
/// re-angle the interior OUT edges that do not already lie inside the boundary
/// wedge — through [`refine_angle_2d`], or, when the vertex carries exactly two
/// interior edges and the 2-D intersection finds nothing, by nudging the angle
/// just inside the wedge (`:996-1000`). The updated angle is looked up by edge
/// (`aDMSR` is a map of edge → angle), so both occurrences of a shared section
/// edge at the vertex get it, the IN one shifted by π (`:1009-1027`).
pub(super) fn refine_angles(face: &Face, map: &mut SmartMap) {
    for slot in 0..map.verts.len() {
        let v = map.verts[slot].clone();
        refine_angles_at_vertex(face, &v, &mut map.infos[slot]);
    }
}

pub(super) fn refine_angles_at_vertex(face: &Face, v: &TopoShape, infos: &mut [EdgeInfo]) {
    let mut a1 = 0.0_f64; // angle of the outgoing boundary edge
    let mut a2 = 0.0_f64; // angle of the incoming boundary edge
    let mut i_cnt_bnd = 0usize;
    let mut i_cnt_int = 0usize;
    for ei in infos.iter() {
        if !ei.is_inside {
            i_cnt_bnd += 1;
            if ei.is_in {
                a2 = ei.angle;
            } else {
                a1 = ei.angle;
            }
        } else {
            i_cnt_int += 1;
        }
    }
    if i_cnt_bnd != 2 {
        return;
    }
    let delta = clockwise_angle(a2, a1);
    let mut refined: HashMap<usize, f64> = HashMap::new();
    for ei in infos.iter() {
        if !ei.is_inside || ei.is_in {
            continue;
        }
        if clockwise_angle(a2, ei.angle) < delta {
            continue; // already inside the wedge
        }
        let ek = GeometryRegistry::shape_key(&ei.edge.0);
        match refine_angle_2d(face, v, &ei.edge, a1, a2, delta) {
            Some(new_val) => {
                refined.insert(ek, new_val);
            }
            None if i_cnt_int == 2 => {
                let new_val = if ei.angle <= a1 {
                    a1 + ANGULAR
                } else {
                    a2 - ANGULAR
                };
                refined.insert(ek, new_val);
            }
            None => {}
        }
    }
    if refined.is_empty() {
        return;
    }
    for ei in infos.iter_mut() {
        if let Some(&a) = refined.get(&GeometryRegistry::shape_key(&ei.edge.0)) {
            ei.angle = if ei.is_in { a + PI } else { a };
        }
    }
}

/// `RefineAngle2D` (`BOPAlgo_WireSplitter_1.cxx:1033-1125`): intersect the
/// edge's p-curve with the ray through the vertex along each boundary angle and
/// take the angle of the point where the curve leaves the wedge. Returns the
/// refined angle, or `None` when neither ray produces a usable point.
/// `BOPAlgo_WireSplitter_1.cxx::RefineAngle2D` (`:1033-1125`).
///
/// Substitutions registered for audit §13 / task T-33:
/// * OCCT intersects through `Geom2dInt_GInter::Perform(aGAC1, aDomain1,
///   aGAC2, aDomain2, aTolInt, aTolInt)` with an explicit `IntRes2d_Domain`
///   (`cxx:1080`); this port calls
///   [`occt_geom2d::geom2d_api::intersect_curves`] and applies the
///   `[t1, t2]` restriction itself (`:1084-1100`). The two agree when the
///   intersector returns every intersection point of that domain.
pub(super) fn refine_angle_2d(
    face: &Face,
    v: &TopoShape,
    e: &Edge,
    a1: f64,
    a2: f64,
    delta: f64,
) -> Option<f64> {
    let cf = 0.01;
    let tol_int = 1.0e-10;
    let c = make_2d(e, face).ok()?;
    let (t1, t2) = BRepTool::edge_parameters(e);
    // `aTV = BRep_Tool::Parameter(aV, aE, myFace)` (`cxx:1060`) — the vertex
    // parameter on the *pcurve*, not on the 3D edge.
    let tv = vertex_parameter_on_face(v, e);
    let pv = c.d0(tv);
    let t_op = if (tv - t1).abs() < (tv - t2).abs() {
        t2
    } else {
        t1
    };
    // `aGAC1.Load(aC2D, aT1, aT2)`: the intersection is restricted to the
    // CurveOnSurface range, which for a bounded p-curve is its own domain.
    let max_dt = 0.3 * (t2 - t1);
    for i in 0..2 {
        let ai = if i == 0 { a1 } else { a2 + PI };
        let dir = GpDir2d::from_vec2d(&GpVec2d::new(ai.cos(), ai.sin())).ok()?;
        let line = Geom2dLine::from_pnt_dir(pv, dir);
        let points = occt_geom2d::geom2d_api::intersect_curves(c.as_ref(), &line, tol_int);
        let mut t1max = tv;
        let mut t2max = -1.0_f64;
        for p in points {
            if p.u1 < t1 - tol_int || p.u1 > t2 + tol_int {
                continue;
            }
            if p.u2 > t2max && (p.u1 - tv).abs() < max_dt {
                t2max = p.u2;
                t1max = p.u1;
            }
        }
        if t2max <= 0.0 {
            continue;
        }
        let d_t = t_op - t1max;
        if d_t.abs() < tol_int {
            continue;
        }
        let t = t1max + cf * d_t;
        let p_val = c.d0(t);
        let vec2 = GpVec2d::new(p_val.x() - pv.x(), p_val.y() - pv.y());
        let Ok(dir2) = GpDir2d::from_vec2d(&vec2) else {
            continue;
        };
        let angle = dir2d_angle(&dir2);
        if clockwise_angle(a2, angle) < delta {
            return Some(angle);
        }
    }
    None
}

