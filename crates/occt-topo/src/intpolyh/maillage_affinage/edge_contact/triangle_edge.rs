//! `TriangleEdgeContact` (`IntPolyh_MaillageAffinage.cxx:2550-3138`).


use super::super::super::point::IntPolyhPoint;
use super::super::super::start_point::IntPolyhStartPoint;
use super::super::super::triangle::IntPolyhTriangle;
use super::super::{MaillageAffinage, MY_CONFUSION_PRECISION};
use super::calcul_pts_inter_tri_edge_coplanaires;

impl MaillageAffinage<'_> {
    /// `TriangleEdgeContact(TriSurfID, EdgeIndex, Tri1, Tri2, PT1, PT2, PT3,
    /// Cote12, Cote23, Cote31, PE1, PE2, Edge, NormaleT, SP1, SP2)`
    /// (`cxx:2550-3138`).
    #[allow(clippy::too_many_arguments)]
    pub fn triangle_edge_contact(
        &self,
        tri_surf_id: i32,
        edge_index: i32,
        tri1: &IntPolyhTriangle,
        tri2: &IntPolyhTriangle,
        pt1: &IntPolyhPoint,
        pt2: &IntPolyhPoint,
        pt3: &IntPolyhPoint,
        cote12: &IntPolyhPoint,
        cote23: &IntPolyhPoint,
        cote31: &IntPolyhPoint,
        pe1: &IntPolyhPoint,
        pe2: &IntPolyhPoint,
        edge: &IntPolyhPoint,
        normale_t: &IntPolyhPoint,
        sp1: &mut IntPolyhStartPoint,
        sp2: &mut IntPolyhStartPoint,
    ) -> i32 {
        // `cxx:2570-2584`.
        if tri_surf_id == 1 {
            sp1.e2 = tri2.edge_number(edge_index);
            sp2.e2 = tri2.edge_number(edge_index);
        } else if tri_surf_id == 2 {
            sp1.e1 = tri1.edge_number(edge_index);
            sp2.e1 = tri1.edge_number(edge_index);
        }

        let mut nb_points = 0;

        if normale_t.square_modulus() == 0.0 {
            // `cxx:2589-2591`.
        } else if cote12.square_modulus() == 0.0
            || cote23.square_modulus() == 0.0
            || cote31.square_modulus() == 0.0
        {
            // `cxx:2592-2595`.
        } else if edge.square_modulus() == 0.0 {
            // `cxx:2596-2598`.
        } else {
            let pe1v = normale_t.dot(pe1);
            let pe2v = normale_t.dot(pe2);
            let pt1v = normale_t.dot(pt1);

            if (pe1v - pt1v).abs() < MY_CONFUSION_PRECISION
                && (pe2v - pt1v).abs() < MY_CONFUSION_PRECISION
            {
                // `cxx:2606-2801`: edge and triangle are coplanar (at most two
                // contact points). The tops of the triangle are projected on
                // the perpendicular to the edge.
                let mut perp_edge = IntPolyhPoint::default();
                perp_edge.cross(normale_t, edge);
                let pp1 = perp_edge.dot(pt1);
                let pp2 = perp_edge.dot(pt2);
                let pp3 = perp_edge.dot(pt3);
                let ppe1 = perp_edge.dot(pe1);

                if (pp1 - pp2).abs() < MY_CONFUSION_PRECISION
                    && (pp1 - pp3).abs() < MY_CONFUSION_PRECISION
                {
                    // `cxx:2619-2622`.
                } else {
                    if ((pp1 >= ppe1) && (pp2 <= ppe1) && (pp3 <= ppe1))
                        || ((pp1 <= ppe1) && (pp2 >= ppe1) && (pp3 >= ppe1))
                    {
                        // Two sides with the common top PT1.
                        calcul_pts_inter_tri_edge_coplanaires(
                            tri_surf_id,
                            normale_t,
                            tri1,
                            tri2,
                            pe1,
                            pe2,
                            edge,
                            edge_index,
                            pt1,
                            pt2,
                            cote12,
                            1,
                            sp1,
                            sp2,
                            &mut nb_points,
                        );

                        if nb_points > 1
                            && (sp2.u1 - sp1.u1).abs() < MY_CONFUSION_PRECISION
                            && (sp1.v1 - sp2.v1).abs() < MY_CONFUSION_PRECISION
                        {
                            nb_points = 1;
                        }

                        // Second side.
                        if nb_points < 2 {
                            calcul_pts_inter_tri_edge_coplanaires(
                                tri_surf_id,
                                normale_t,
                                tri1,
                                tri2,
                                pe1,
                                pe2,
                                edge,
                                edge_index,
                                pt3,
                                pt1,
                                cote31,
                                3,
                                sp1,
                                sp2,
                                &mut nb_points,
                            );
                        }
                    }

                    if nb_points > 1
                        && (sp1.u1 - sp2.u1).abs() < MY_CONFUSION_PRECISION
                        && (sp1.v2 - sp2.v1).abs() < MY_CONFUSION_PRECISION
                    {
                        nb_points = 1;
                    }
                    if nb_points >= 2 {
                        return nb_points;
                    } else if (((pp2 >= ppe1) && (pp1 <= ppe1) && (pp3 <= ppe1))
                        || ((pp2 <= ppe1) && (pp1 >= ppe1) && (pp3 >= ppe1)))
                        && nb_points < 2
                    {
                        // Two sides with the common top PT2.
                        calcul_pts_inter_tri_edge_coplanaires(
                            tri_surf_id,
                            normale_t,
                            tri1,
                            tri2,
                            pe1,
                            pe2,
                            edge,
                            edge_index,
                            pt1,
                            pt2,
                            cote12,
                            1,
                            sp1,
                            sp2,
                            &mut nb_points,
                        );

                        if nb_points > 1
                            && (sp2.u1 - sp1.u1).abs() < MY_CONFUSION_PRECISION
                            && (sp1.v1 - sp2.v1).abs() < MY_CONFUSION_PRECISION
                        {
                            nb_points = 1;
                        }

                        // Second side.
                        if nb_points < 2 {
                            calcul_pts_inter_tri_edge_coplanaires(
                                tri_surf_id,
                                normale_t,
                                tri1,
                                tri2,
                                pe1,
                                pe2,
                                edge,
                                edge_index,
                                pt2,
                                pt3,
                                cote23,
                                2,
                                sp1,
                                sp2,
                                &mut nb_points,
                            );
                        }
                    }
                    if nb_points > 1
                        && (sp2.u1 - sp1.u1).abs() < MY_CONFUSION_PRECISION
                        && (sp1.v1 - sp2.v1).abs() < MY_CONFUSION_PRECISION
                    {
                        nb_points = 1;
                    }
                    if nb_points >= 2 {
                        return nb_points;
                    } else if (((pp3 >= ppe1) && (pp1 <= ppe1) && (pp2 <= ppe1))
                        || ((pp3 <= ppe1) && (pp1 >= ppe1) && (pp2 >= ppe1)))
                        && nb_points < 2
                    {
                        // Two sides with the common top PT3.
                        calcul_pts_inter_tri_edge_coplanaires(
                            tri_surf_id,
                            normale_t,
                            tri1,
                            tri2,
                            pe1,
                            pe2,
                            edge,
                            edge_index,
                            pt3,
                            pt1,
                            cote31,
                            3,
                            sp1,
                            sp2,
                            &mut nb_points,
                        );

                        if nb_points > 1
                            && (sp2.u1 - sp1.u1).abs() < MY_CONFUSION_PRECISION
                            && (sp1.v1 - sp2.v1).abs() < MY_CONFUSION_PRECISION
                        {
                            nb_points = 1;
                        }

                        // Second side.
                        if nb_points < 2 {
                            calcul_pts_inter_tri_edge_coplanaires(
                                tri_surf_id,
                                normale_t,
                                tri1,
                                tri2,
                                pe1,
                                pe2,
                                edge,
                                edge_index,
                                pt2,
                                pt3,
                                cote23,
                                2,
                                sp1,
                                sp2,
                                &mut nb_points,
                            );
                        }
                    }
                    if nb_points > 1
                        && (sp2.u1 - sp1.u1).abs() < MY_CONFUSION_PRECISION
                        && (sp2.v1 - sp1.v1).abs() < MY_CONFUSION_PRECISION
                    {
                        nb_points = 1;
                    }
                    if nb_points >= 2 {
                        return nb_points;
                    }
                }
            } else if ((pe1v >= pt1v) && (pt1v >= pe2v))
                || ((pe1v <= pt1v) && (pt1v <= pe2v))
            {
                // `cxx:2803-3130`: non-coplanar edge and triangle (one contact
                // point).
                let lambda = (pe1v - pt1v) / (pe1v - pe2v);
                let mut pi = IntPolyhPoint::default();
                if lambda < -MY_CONFUSION_PRECISION {
                    // `cxx:2809-2811`: rejected, `PI` keeps its default value
                    // and the computation below still runs.
                } else if lambda.abs() < MY_CONFUSION_PRECISION {
                    // lambda == 0
                    pi = *pe1;
                    if tri_surf_id == 1 {
                        sp1.e2 = -1;
                    } else {
                        sp1.e1 = -1;
                    }
                } else if (lambda - 1.0).abs() < MY_CONFUSION_PRECISION {
                    // lambda == 1
                    pi = *pe2;
                    if tri_surf_id == 1 {
                        sp1.e2 = -1;
                    } else {
                        sp1.e1 = -1;
                    }
                } else {
                    pi = pe1.add(&edge.multiplication(lambda));
                    if tri_surf_id == 1 {
                        if tri2.edge_orientation(edge_index) > 0 {
                            sp1.lambda2 = lambda;
                        } else {
                            sp1.lambda2 = 1.0 - lambda;
                        }
                    }
                    if tri_surf_id == 2 {
                        if tri1.edge_orientation(edge_index) > 0 {
                            sp1.lambda1 = lambda;
                        } else {
                            sp1.lambda1 = 1.0 - lambda;
                        }
                    }
                }

                let cote23x = cote23.x;
                let mut d1 = 0.0;
                let alpha;
                let beta;

                // Combination Eq1 Eq2.
                if cote23x.abs() > MY_CONFUSION_PRECISION {
                    d1 = cote12.y - cote12.x * cote23.y / cote23x;
                }
                if d1.abs() > MY_CONFUSION_PRECISION {
                    alpha = (pi.y - pt1.y - (pi.x - pt1.x) * cote23.y / cote23x) / d1;
                    if alpha < -MY_CONFUSION_PRECISION || alpha > 1.0 + MY_CONFUSION_PRECISION {
                        return 0;
                    }
                    beta = (pi.x - pt1.x - alpha * cote12.x) / cote23x;
                } else if cote12.x.abs() > MY_CONFUSION_PRECISION
                    && cote23x.abs() < MY_CONFUSION_PRECISION
                {
                    // There is Cote23.X() == 0.
                    alpha = (pi.x - pt1.x) / cote12.x;
                    if alpha < -MY_CONFUSION_PRECISION || alpha > 1.0 + MY_CONFUSION_PRECISION {
                        return 0;
                    }
                    if cote23.y.abs() > MY_CONFUSION_PRECISION {
                        beta = (pi.y - pt1.y - alpha * cote12.y) / cote23.y;
                    } else if cote23.z.abs() > MY_CONFUSION_PRECISION {
                        beta = (pi.z - pt1.z - alpha * cote12.z) / cote23.z;
                    } else {
                        beta = 0.0;
                    }
                } else if cote23.x.abs() > MY_CONFUSION_PRECISION
                    && (cote12.z - cote12.x * cote23.z / cote23.x).abs()
                        > MY_CONFUSION_PRECISION
                {
                    // Combination Eq1 and Eq3.
                    let d3 = cote12.z - cote12.x * cote23.z / cote23.x;
                    alpha = (pi.z - pt1.z - (pi.x - pt1.x) * cote23.z / cote23.x) / d3;
                    if alpha < -MY_CONFUSION_PRECISION || alpha > 1.0 + MY_CONFUSION_PRECISION {
                        return 0;
                    }
                    beta = (pi.x - pt1.x - alpha * cote12.x) / cote23.x;
                } else if cote23.y.abs() > MY_CONFUSION_PRECISION
                    && (cote12.z - cote12.y * cote23.z / cote23.y).abs()
                        > MY_CONFUSION_PRECISION
                {
                    // Combination Eq2 and Eq3.
                    let d4 = cote12.z - cote12.y * cote23.z / cote23.y;
                    alpha = (pi.z - pt1.z - (pi.y - pt1.y) * cote23.z / cote23.y) / d4;
                    if alpha < -MY_CONFUSION_PRECISION || alpha > 1.0 + MY_CONFUSION_PRECISION {
                        return 0;
                    }
                    beta = (pi.y - pt1.y - alpha * cote12.y) / cote23.y;
                } else if cote12.y.abs() > MY_CONFUSION_PRECISION
                    && cote23.y.abs() < MY_CONFUSION_PRECISION
                {
                    // Combination Eq2 and Eq3 with Cote23.Y() == 0.
                    alpha = (pi.y - pt1.y) / cote12.y;
                    if alpha < -MY_CONFUSION_PRECISION || alpha > 1.0 + MY_CONFUSION_PRECISION {
                        return 0;
                    }
                    if cote23.z.abs() > MY_CONFUSION_PRECISION {
                        beta = (pi.z - pt1.z - alpha * cote12.z) / cote23.z;
                    } else {
                        // `cxx:2959-2961` prints "Cote PT2PT3 nul1" here.
                        beta = 0.0;
                    }
                } else if cote12.z.abs() > MY_CONFUSION_PRECISION
                    && cote23.z.abs() < MY_CONFUSION_PRECISION
                {
                    // Combination Eq1 and Eq3 with Cote23.Z() == 0.
                    alpha = (pi.z - pt1.z) / cote12.z;
                    if alpha < -MY_CONFUSION_PRECISION || alpha > 1.0 + MY_CONFUSION_PRECISION {
                        return 0;
                    }
                    if cote23.x.abs() > MY_CONFUSION_PRECISION {
                        beta = (pi.x - pt1.x - alpha * cote12.x) / cote23.x;
                    } else {
                        beta = 0.0;
                    }
                } else {
                    // Particular case not processed.
                    alpha = super::super::REAL_LAST;
                    beta = super::super::REAL_LAST;
                }

                if beta < -MY_CONFUSION_PRECISION || beta > alpha + MY_CONFUSION_PRECISION {
                    return 0;
                }

                sp1.set_xyz(pi.x, pi.y, pi.z);

                if tri_surf_id == 1 {
                    sp1.set_uv2(pi.u, pi.v);
                    sp1.set_uv1(
                        pt1.u + cote12.u * alpha + cote23.u * beta,
                        pt1.v + cote12.v * alpha + cote23.v * beta,
                    );
                    nb_points += 1;
                    if alpha < MY_CONFUSION_PRECISION {
                        // alpha == 0 -> beta == 0
                        sp1.set_xyz(pt1.x, pt1.y, pt1.z);
                        sp1.set_uv1(pt1.u, pt1.v);
                        sp1.e1 = -1;
                    } else if beta < MY_CONFUSION_PRECISION && (1.0 - alpha).abs() < MY_CONFUSION_PRECISION
                    {
                        // beta == 0, alpha == 1
                        sp1.set_xyz(pt2.x, pt2.y, pt2.z);
                        sp1.set_uv1(pt2.u, pt2.v);
                        sp1.e1 = -1;
                    } else if (beta - 1.0).abs() < MY_CONFUSION_PRECISION
                        && (1.0 - alpha).abs() < MY_CONFUSION_PRECISION
                    {
                        // beta == 1, alpha == 1
                        sp1.set_xyz(pt3.x, pt3.y, pt3.z);
                        sp1.set_uv1(pt3.u, pt3.v);
                        sp1.e1 = -1;
                    } else if beta < MY_CONFUSION_PRECISION {
                        // beta == 0
                        sp1.e1 = tri1.edge_number(1);
                        if tri1.edge_orientation(1) > 0 {
                            sp1.lambda1 = alpha;
                        } else {
                            sp1.lambda1 = 1.0 - alpha;
                        }
                    } else if (beta - alpha).abs() < MY_CONFUSION_PRECISION {
                        // beta == alpha
                        sp1.e1 = tri1.edge_number(3);
                        if tri1.edge_orientation(3) > 0 {
                            sp1.lambda1 = 1.0 - alpha;
                        } else {
                            sp1.lambda1 = alpha;
                        }
                    } else if (alpha - 1.0).abs() < MY_CONFUSION_PRECISION {
                        // alpha == 1
                        sp1.e1 = tri1.edge_number(2);
                        if tri1.edge_orientation(2) > 0 {
                            sp1.lambda1 = beta;
                        } else {
                            sp1.lambda1 = 1.0 - beta;
                        }
                    }
                } else if tri_surf_id == 2 {
                    sp1.set_uv1(pi.u, pi.v);
                    sp1.set_uv2(
                        pt1.u + cote12.u * alpha + cote23.u * beta,
                        pt1.v + cote12.v * alpha + cote23.v * beta,
                    );
                    nb_points += 1;
                    if alpha < MY_CONFUSION_PRECISION {
                        // alpha == 0 -> beta == 0
                        sp1.set_xyz(pt1.x, pt1.y, pt1.z);
                        sp1.set_uv2(pt1.u, pt1.v);
                        sp1.e2 = -1;
                    } else if beta < MY_CONFUSION_PRECISION && (1.0 - alpha).abs() < MY_CONFUSION_PRECISION
                    {
                        // beta == 0, alpha == 1
                        sp1.set_xyz(pt2.x, pt2.y, pt2.z);
                        sp1.set_uv2(pt2.u, pt2.v);
                        sp1.e2 = -1;
                    } else if (beta - 1.0).abs() < MY_CONFUSION_PRECISION
                        && (1.0 - alpha).abs() < MY_CONFUSION_PRECISION
                    {
                        // beta == 1, alpha == 1
                        sp1.set_xyz(pt3.x, pt3.y, pt3.z);
                        sp1.set_uv2(pt3.u, pt3.v);
                        sp1.e2 = -1;
                    } else if beta < MY_CONFUSION_PRECISION {
                        // beta == 0
                        sp1.e2 = tri2.edge_number(1);
                        if tri2.edge_orientation(1) > 0 {
                            sp1.lambda2 = alpha;
                        } else {
                            sp1.lambda2 = 1.0 - alpha;
                        }
                    } else if (beta - alpha).abs() < MY_CONFUSION_PRECISION {
                        // beta == alpha
                        sp1.e2 = tri2.edge_number(3);
                        if tri2.edge_orientation(3) > 0 {
                            sp1.lambda2 = 1.0 - alpha;
                        } else {
                            sp1.lambda2 = alpha;
                        }
                    } else if (alpha - 1.0).abs() < MY_CONFUSION_PRECISION {
                        // alpha == 1
                        sp1.e2 = tri2.edge_number(2);
                        if tri2.edge_orientation(2) > 0 {
                            sp1.lambda2 = alpha;
                        } else {
                            sp1.lambda2 = 1.0 - alpha;
                        }
                    }
                }
            } else {
                return 0;
            }
        }

        nb_points
    }
}
