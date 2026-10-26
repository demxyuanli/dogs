//! `CalculPtsInterTriEdgeCoplanaires(...)` (`IntPolyh_MaillageAffinage.cxx:2051-2546`).

use occt_core::gp::GpVec;

use super::super::super::point::IntPolyhPoint;
use super::super::super::start_point::IntPolyhStartPoint;
use super::super::super::triangle::IntPolyhTriangle;
use super::super::{MY_CONFUSION_PRECISION, SQUARE_MY_CONFUSION_PRECISION};

/// `CalculPtsInterTriEdgeCoplanaires(...)` (`cxx:2051-2546`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn calcul_pts_inter_tri_edge_coplanaires(
    tri_surf_id: i32,
    normale_tri: &IntPolyhPoint,
    tri1: &IntPolyhTriangle,
    tri2: &IntPolyhTriangle,
    pe1: &IntPolyhPoint,
    pe2: &IntPolyhPoint,
    edge: &IntPolyhPoint,
    edge_index: i32,
    pt1: &IntPolyhPoint,
    pt2: &IntPolyhPoint,
    cote: &IntPolyhPoint,
    cote_index: i32,
    sp1: &mut IntPolyhStartPoint,
    sp2: &mut IntPolyhStartPoint,
    nb_points: &mut i32,
) {
    // `cxx:2068-2083`: the `.cxx` divides by the *square* magnitude.
    let mut a_ve = GpVec::new(edge.x, edge.y, edge.z);
    let mut a_vc = GpVec::new(cote.x, cote.y, cote.z);
    let a_de = a_ve.square_magnitude();
    let a_dc = a_vc.square_magnitude();
    if a_de > SQUARE_MY_CONFUSION_PRECISION {
        a_ve = a_ve.divided(a_de);
    }
    if a_dc > SQUARE_MY_CONFUSION_PRECISION {
        a_vc = a_vc.divided(a_dc);
    }

    if !a_ve.is_parallel_ang(&a_vc, MY_CONFUSION_PRECISION) {
        // `cxx:2086-2277`: the edge and the side are not parallel.
        let mut per = IntPolyhPoint::default();
        per.cross(normale_tri, cote);
        let p1p = per.dot(pe1);
        let p2p = per.dot(pe2);
        let p0p = per.dot(pt1);

        if (((p1p >= p0p) && (p0p >= p2p)) || ((p1p <= p0p) && (p0p <= p2p)))
            && ((p1p - p2p).abs() > MY_CONFUSION_PRECISION)
        {
            let lambda = (p1p - p0p) / (p1p - p2p);
            let pie = if lambda.abs() < MY_CONFUSION_PRECISION {
                // lambda == 0
                *pe1
            } else if lambda.abs() > 1.0 - MY_CONFUSION_PRECISION {
                // lambda == 1
                *pe2
            } else {
                pe1.add(&edge.multiplication(lambda))
            };

            let mut alpha = super::super::REAL_LAST;
            if cote.x != 0.0 {
                alpha = (pie.x - pt1.x) / cote.x;
            } else if cote.y != 0.0 {
                alpha = (pie.y - pt1.y) / cote.y;
            } else if cote.z != 0.0 {
                alpha = (pie.z - pt1.z) / cote.z;
            }

            if alpha < -MY_CONFUSION_PRECISION {
                // `cxx:2134-2136`: rejected.
            } else if *nb_points == 0 {
                sp1.set_xyz(pie.x, pie.y, pie.z);
                if tri_surf_id == 1 {
                    if alpha.abs() < MY_CONFUSION_PRECISION {
                        // alpha == 0
                        sp1.set_uv1(pt1.u, pt1.v);
                        sp1.set_uv2(pie.u, pie.v);
                        sp1.e1 = -1;
                    }
                    if alpha.abs() > 1.0 - MY_CONFUSION_PRECISION {
                        // alpha == 1
                        sp1.set_uv1(pt2.u, pt2.v);
                        sp1.set_uv2(pie.u, pie.v);
                        sp1.e1 = -1;
                    } else {
                        sp1.set_uv1(pt1.u + cote.u * alpha, pt1.v + cote.v * alpha);
                        sp1.set_uv2(pie.u, pie.v);
                        sp1.e1 = tri1.edge_number(cote_index);
                        if tri1.edge_orientation(cote_index) > 0 {
                            sp1.lambda1 = alpha;
                        } else {
                            sp1.lambda1 = 1.0 - alpha;
                        }
                    }
                    *nb_points += 1;
                } else if tri_surf_id == 2 {
                    if alpha.abs() < MY_CONFUSION_PRECISION {
                        // alpha == 0
                        sp1.set_uv2(pt1.u, pt1.v);
                        sp1.set_uv1(pie.u, pie.v);
                        sp1.e2 = -1;
                    }
                    if alpha.abs() > 1.0 - MY_CONFUSION_PRECISION {
                        // alpha == 1
                        sp1.set_uv2(pt2.u, pt2.v);
                        sp1.set_uv1(pie.u, pie.v);
                        sp1.e2 = -1;
                    } else {
                        sp1.set_uv1(pie.u, pie.v);
                        sp1.set_uv2(pt1.u + cote.u * alpha, pt1.v + cote.v * alpha);
                        sp1.e2 = tri2.edge_number(cote_index);
                        if tri2.edge_orientation(cote_index) > 0 {
                            sp1.lambda2 = alpha;
                        } else {
                            sp1.lambda2 = 1.0 - alpha;
                        }
                    }
                    *nb_points += 1;
                }
            } else if *nb_points == 1 {
                sp2.set_xyz(pie.x, pie.y, pie.z);
                if tri_surf_id == 1 {
                    if alpha.abs() < MY_CONFUSION_PRECISION {
                        // alpha == 0
                        sp2.set_uv1(pt1.u, pt1.v);
                        sp2.set_uv2(pie.u, pie.v);
                        sp2.e1 = -1;
                    }
                    if alpha.abs() > 1.0 - MY_CONFUSION_PRECISION {
                        // alpha == 1
                        sp2.set_uv1(pt2.u, pt2.v);
                        sp2.set_uv2(pie.u, pie.v);
                        sp2.e1 = -1;
                    } else {
                        sp2.set_uv1(pt1.u + cote.u * alpha, pt1.v + cote.v * alpha);
                        sp2.set_uv2(pie.u, pie.v);
                        sp2.e1 = tri1.edge_number(cote_index);
                        if tri1.edge_orientation(cote_index) > 0 {
                            sp2.lambda1 = alpha;
                        } else {
                            sp2.lambda1 = 1.0 - alpha;
                        }
                    }
                    *nb_points += 1;
                } else if tri_surf_id == 2 {
                    if alpha.abs() < MY_CONFUSION_PRECISION {
                        // alpha == 0
                        sp2.set_uv2(pt1.u, pt1.v);
                        sp2.set_uv1(pie.u, pie.v);
                        sp2.e2 = -1;
                    }
                    if alpha.abs() > 1.0 - MY_CONFUSION_PRECISION {
                        // alpha == 1
                        sp2.set_uv2(pt2.u, pt2.v);
                        sp2.set_uv1(pie.u, pie.v);
                        sp2.e2 = -1;
                    } else {
                        sp2.set_uv1(pie.u, pie.v);
                        sp2.set_uv2(pt1.u + cote.u * alpha, pt1.v + cote.v * alpha);
                        sp2.e2 = tri2.edge_number(cote_index);
                        if tri1.edge_orientation(cote_index) > 0 {
                            sp2.lambda2 = alpha;
                        } else {
                            sp2.lambda2 = 1.0 - alpha;
                        }
                    }
                    *nb_points += 1;
                }
            }
        }
    } else {
        // `cxx:2278-2508`: side and edge are parallel; the points are projected
        // on that side.
        let pe1p = cote.dot(pe1);
        let pe2p = cote.dot(pe2);
        let pt1p = cote.dot(pt1);
        let pt2p = cote.dot(pt2);
        let mut lambda1 = 0.0;
        let mut lambda2 = 0.0;
        let mut alpha1 = 0.0;
        let mut alpha2 = 0.0;
        let mut pep1 = IntPolyhPoint::default();
        let mut ptp1 = IntPolyhPoint::default();
        let mut pep2 = IntPolyhPoint::default();
        let mut ptp2 = IntPolyhPoint::default();

        if pe1p > pe2p {
            if (pt1p < pe1p) && (pe1p <= pt2p) {
                lambda1 = 0.0;
                pep1 = *pe1;
                alpha1 = (pe1p - pt1p) / (pt2p - pt1p);
                ptp1 = pt1.add(&cote.multiplication(alpha1));
                *nb_points = 1;
                if pt1p <= pe2p {
                    lambda2 = 1.0;
                    pep2 = *pe2;
                    alpha2 = (pe2p - pt1p) / (pt2p - pt1p);
                    ptp2 = pt1.add(&cote.multiplication(alpha2));
                    *nb_points = 2;
                } else {
                    lambda2 = (pt1p - pe1p) / (pe2p - pe1p);
                    pep2 = pe1.add(&edge.multiplication(lambda2));
                    alpha2 = 0.0;
                    ptp2 = *pt1;
                    *nb_points = 2;
                }
            } else if (pt2p < pe1p) && (pe1p <= pt1p) {
                lambda1 = 0.0;
                pep1 = *pe1;
                alpha1 = (pt1p - pe1p) / (pt1p - pt2p);
                ptp1 = pt1.add(&cote.multiplication(alpha1));
                *nb_points = 1;
                if pt2p <= pe2p {
                    lambda2 = 1.0;
                    pep2 = *pe2;
                    alpha2 = (pe2p - pt1p) / (pt2p - pt1p);
                    ptp2 = pt1.add(&cote.multiplication(alpha2));
                    *nb_points = 2;
                } else {
                    lambda2 = (pt2p - pe1p) / (pe2p - pe1p);
                    pep2 = pe1.add(&edge.multiplication(lambda2));
                    alpha2 = 1.0;
                    ptp2 = *pt2;
                    *nb_points = 2;
                }
            }
        }

        if pe1p < pe2p {
            if (pt1p < pe2p) && (pe2p <= pt2p) {
                lambda1 = 1.0;
                pep1 = *pe2;
                alpha1 = (pe2p - pt1p) / (pt2p - pt1p);
                ptp1 = pt1.add(&cote.multiplication(alpha1));
                *nb_points = 1;
                if pt1p <= pe1p {
                    lambda2 = 0.0;
                    pep2 = *pe1;
                    alpha2 = (pe1p - pt1p) / (pt2p - pt1p);
                    ptp2 = pt1.add(&cote.multiplication(alpha2));
                    *nb_points = 2;
                } else {
                    lambda2 = (pt1p - pe1p) / (pe2p - pe1p);
                    pep2 = pe2.add(&edge.multiplication(lambda2));
                    alpha2 = 0.0;
                    ptp2 = *pt1;
                    *nb_points = 2;
                }
            } else if (pt2p < pe2p) && (pe2p <= pt1p) {
                lambda1 = 1.0;
                pep1 = *pe2;
                alpha1 = (pt1p - pe2p) / (pt1p - pt2p);
                ptp1 = pt1.add(&cote.multiplication(alpha1));
                *nb_points = 1;
                if pt2p <= pe1p {
                    lambda2 = 0.0;
                    pep2 = *pe1;
                    alpha2 = (pe1p - pt1p) / (pt2p - pt1p);
                    ptp2 = pt1.add(&cote.multiplication(alpha2));
                    *nb_points = 2;
                } else {
                    lambda2 = (pt2p - pe1p) / (pe2p - pe1p);
                    pep2 = pe1.add(&edge.multiplication(lambda2));
                    alpha2 = 1.0;
                    ptp2 = *pt2;
                    *nb_points = 2;
                }
            }
        }

        if *nb_points != 0 {
            sp1.set_xyz(pep1.x, pep1.y, pep1.z);
            if tri_surf_id == 1 {
                // cote belongs to Tri1
                sp1.set_uv1(ptp1.u, ptp1.v);
                sp1.set_uv2(pep1.u, pep1.v);
                sp1.e1 = tri1.edge_number(cote_index);
                if tri1.edge_orientation(cote_index) > 0 {
                    sp1.lambda1 = alpha1;
                } else {
                    sp1.lambda1 = 1.0 - alpha1;
                }
                if tri2.edge_orientation(edge_index) > 0 {
                    sp1.lambda2 = lambda1;
                } else {
                    sp1.lambda2 = 1.0 - lambda1;
                }
            }
            if tri_surf_id == 2 {
                // cote belongs to Tri2
                sp1.set_uv1(pep1.u, ptp1.v);
                sp1.set_uv2(ptp1.u, pep1.v);
                sp1.e2 = tri2.edge_number(cote_index);
                if tri2.edge_orientation(cote_index) > 0 {
                    sp1.lambda1 = alpha1;
                } else {
                    sp1.lambda1 = 1.0 - alpha1;
                }
                if tri1.edge_orientation(edge_index) > 0 {
                    sp1.lambda2 = lambda1;
                } else {
                    sp1.lambda2 = 1.0 - lambda1;
                }
            }

            // `cxx:2448-2451`: PEP1 and PEP2 are the same point.
            if *nb_points == 2
                && (pep1.u - pep2.u).abs() < MY_CONFUSION_PRECISION
                && (pep1.v - pep2.v).abs() < MY_CONFUSION_PRECISION
            {
                *nb_points = 1;
            }
            if *nb_points == 2 {
                sp2.set_xyz(pep2.x, pep2.y, pep2.z);
                if tri_surf_id == 1 {
                    sp2.set_uv1(ptp2.u, ptp2.v);
                    sp2.set_uv2(pep2.u, pep2.v);
                    sp2.e1 = tri1.edge_number(cote_index);
                    if tri1.edge_orientation(cote_index) > 0 {
                        sp2.lambda1 = alpha1;
                    } else {
                        sp2.lambda1 = 1.0 - alpha1;
                    }
                    if tri2.edge_orientation(edge_index) > 0 {
                        sp2.lambda2 = lambda1;
                    } else {
                        sp2.lambda2 = 1.0 - lambda1;
                    }
                }
                if tri_surf_id == 2 {
                    sp2.set_uv1(pep2.u, ptp2.v);
                    sp2.set_uv2(ptp2.u, pep2.v);
                    sp2.e2 = tri2.edge_number(cote_index);
                    if tri1.edge_orientation(cote_index) > 0 {
                        sp2.lambda1 = alpha1;
                    } else {
                        sp2.lambda1 = 1.0 - alpha1;
                    }
                    if tri2.edge_orientation(edge_index) > 0 {
                        sp2.lambda2 = lambda1;
                    } else {
                        sp2.lambda2 = 1.0 - lambda1;
                    }
                }
            }
        }
    }

    // `cxx:2508-2545`: a point placed on a top has its edge reset.
    if *nb_points > 0 {
        if sp1.lambda1.abs() < MY_CONFUSION_PRECISION {
            sp1.e1 = -1;
        }
        if (sp1.lambda1 - 1.0).abs() < MY_CONFUSION_PRECISION {
            sp1.e1 = -1;
        }
        if sp1.lambda2.abs() < MY_CONFUSION_PRECISION {
            sp1.e2 = -1;
        }
        if (sp1.lambda2 - 1.0).abs() < MY_CONFUSION_PRECISION {
            sp1.e2 = -1;
        }
    }
    if *nb_points == 2 {
        if sp2.lambda1.abs() < MY_CONFUSION_PRECISION {
            sp2.e1 = -1;
        }
        if (sp2.lambda1 - 1.0).abs() < MY_CONFUSION_PRECISION {
            sp2.e1 = -1;
        }
        if sp2.lambda2.abs() < MY_CONFUSION_PRECISION {
            sp2.e2 = -1;
        }
        if (sp2.lambda2 - 1.0).abs() < MY_CONFUSION_PRECISION {
            sp2.e2 = -1;
        }
    }
}
