use super::next_param::next_param_on_pc;
use super::pcurve::pcurve_and_range_on_face;
use super::wire_edges::{edge_oriented_ends, same_vertex};

use super::super::prelude::*;

/// `BRepTools_WireExplorer::Next` multi-candidate selection
/// (`BRepTools_WireExplorer.cxx:521-705`): returns the index into `list` of the
/// next edge, or `None` when the walk stops (`myEdge = TopoDS_Edge()`).
pub(super) fn wire_next_by_uv(
    edge: &Edge,
    face: &Face,
    list: &[Edge],
    tol_u: f64,
    tol_v: f64,
) -> Option<usize> {
    let (pcurve, df_f_par, df_l_par) = pcurve_and_range_on_face(edge, face)?;
    // The 2D point of the current vertex on the current edge (`cxx:548-547`):
    // MAX pcurve parameter for a FORWARD edge, MIN for a REVERSED one.
    let rev = edge.0.orientation().is_reversed();
    let p_ref = if rev {
        pcurve.d0(df_f_par)
    } else {
        pcurve.d0(df_l_par)
    };
    // The direction from the current vertex back into the current edge
    // (`cxx:551-556`), used as the angle reference.
    let df_m_par =
        next_param_on_pc(pcurve.as_ref(), &p_ref, df_f_par, df_l_par, tol_u, tol_v, rev);
    let p_refm = pcurve.d0(df_m_par);
    let an_e_ref_dir = GpVec2d::new(p_refm.x() - p_ref.x(), p_refm.y() - p_ref.y());
    if an_e_ref_dir.square_magnitude() < REAL_SMALL {
        return None;
    }

    let mut k_min = 0usize;
    let mut is_degenerated = true;
    let mut dmin = f64::MAX;
    let mut df_min_angle = 3.0 * PI;
    let mut df_cur_angle = 3.0 * PI;
    for _pass in 0..2 {
        let mut k = 1usize;
        for e in list {
            if Arc::ptr_eq(&e.0.tshape, &edge.0.tshape) {
                k += 1;
                continue;
            }
            let (v1, v2) = edge_oriented_ends(e);
            let (Some(v1), Some(v2)) = (v1, v2) else {
                k += 1;
                continue;
            };
            let Some((pc, f, l)) = pcurve_and_range_on_face(e, face) else {
                k += 1;
                continue;
            };
            // `cxx:602`: pass 1 keeps closed edges (`V1.IsSame(V2)`), pass 2 the
            // rest.
            if same_vertex(&v1, &v2) == is_degenerated {
                let peb = if e.0.orientation().is_reversed() {
                    pc.d0(l)
                } else {
                    pc.d0(f)
                };
                if (l - f).abs() > PCONFUSION {
                    let reverse = !e.0.orientation().is_reversed();
                    let epm =
                        next_param_on_pc(pc.as_ref(), &peb, f, l, tol_u, tol_v, reverse);
                    let mut pee = pc.d0(epm);
                    if peb.square_distance(&pee) <= REAL_SMALL {
                        // Very short curve: take its tangent (`cxx:626-643`).
                        let (_, d) = pc.d1(epm);
                        pee = if e.0.orientation().is_reversed() {
                            GpPnt2d::new(peb.x() - d.x(), peb.y() - d.y())
                        } else {
                            GpPnt2d::new(peb.x() + d.x(), peb.y() + d.y())
                        };
                        if peb.square_distance(&pee) <= REAL_SMALL {
                            k += 1;
                            continue;
                        }
                    }
                    let an_e_dir = GpVec2d::new(pee.x() - peb.x(), pee.y() - peb.y());
                    df_cur_angle = an_e_dir.angle(&an_e_ref_dir).abs();
                }
                if df_cur_angle <= df_min_angle {
                    let mut d = p_ref.square_distance(&peb);
                    if d <= PCONFUSION {
                        d = 0.0;
                    }
                    if (peb.x() - p_ref.x()).abs() < tol_u
                        && (peb.y() - p_ref.y()).abs() < tol_v
                        && d <= dmin
                    {
                        df_min_angle = df_cur_angle;
                        k_min = k;
                        dmin = d;
                    }
                }
            }
            k += 1;
        }
        if k_min == 0 {
            is_degenerated = false;
            dmin = f64::MAX;
        } else {
            break;
        }
    }
    if k_min == 0 { None } else { Some(k_min - 1) }
}
