//! `IntCurveSurface_InterUtils::SectionPointToParameters`.
//!
//! Source: `IntCurveSurface_InterUtils.pxx:736-849` (TKGeomAlgo). This is the
//! only place the numeric `IntCurveSurface_Inter` path turns a section point
//! back into surface `(u, v)`: the parameters come from the **polyhedron** node
//! or triangle the interference engine reported, never from a projection.
//!
//! The `Intf_VERTEX` / `Intf_EDGE` arms interpolate the node parameters; the
//! `Intf_FACE` arm barycentrically blends the triangle's three node parameters
//! (and falls back to the closest edge when the point is degenerate against the
//! triangle).

use occt_core::gp::GpVec;
use occt_core::intf::{IntfPIType, IntfSectionPoint};
use occt_core::precision::COMPUTATIONAL;

use super::polygon::ThePolygonOfHInter;
use super::polyhedron::ThePolyhedronOfHInter;

/// `SectionPointToParameters(Sp, Polyhedron, Polygon, U, V, W)`
/// (`IntCurveSurface_InterUtils.pxx:739-849`). Returns `(U, V, W)`.
pub fn section_point_to_parameters(
    sp: &IntfSectionPoint,
    polyhedron: &ThePolyhedronOfHInter,
    polygon: &ThePolygonOfHInter,
) -> (f64, f64, f64) {
    let p = *sp.pnt();
    let (mut u1, mut v1) = (0.0f64, 0.0f64);

    let (typ, adr1, adr2, param) = sp.info_second_2();
    match typ {
        // `Polyhedron.Parameters(Adr1, u1, v1)` (`:761-765`).
        IntfPIType::Vertex => {
            let (u, v) = polyhedron.parameters(adr1);
            u1 = u;
            v1 = v;
        }
        // `:766-772`.
        IntfPIType::Edge => {
            let (ua, va) = polyhedron.parameters(adr1);
            let (ub, vb) = polyhedron.parameters(adr2);
            u1 = ua + param * (ub - ua);
            v1 = va + param * (vb - va);
        }
        // `:773-834`.
        IntfPIType::Face => {
            let (pt1, pt2, pt3) = polyhedron.triangle_indices(adr1 as usize);
            let pa = polyhedron.point(pt1);
            let pb = polyhedron.point(pt2);
            let pc = polyhedron.point(pt3);
            let (ua, va) = polyhedron.parameters(pt1);
            let (ub, vb) = polyhedron.parameters(pt2);
            let (uc, vc) = polyhedron.parameters(pt3);

            let ab = GpVec::from_pnts(&pa, &pb);
            let bc = GpVec::from_pnts(&pb, &pc);
            let ca_vec = GpVec::from_pnts(&pc, &pa);
            let normale = ab.crossed(&GpVec::from_pnts(&pa, &pc));

            let cc = ab.crossed(&GpVec::from_pnts(&pa, &p)).dot(&normale);
            let ca = bc.crossed(&GpVec::from_pnts(&pb, &p)).dot(&normale);
            let cb = ca_vec.crossed(&GpVec::from_pnts(&pc, &p)).dot(&normale);
            let cabc = ca + cb + cc;

            if cabc.abs() > COMPUTATIONAL {
                // `:788-795`: barycentric blend of the three node parameters.
                let ca = ca / cabc;
                let cb = cb / cabc;
                let cc = cc / cabc;
                u1 = ca * ua + cb * ub + cc * uc;
                v1 = ca * va + cb * vb + cc * vc;
            } else {
                // `:797-832`: degenerate against the triangle plane, project
                // onto the longest edge and interpolate there.
                let sq_ab = ab.square_magnitude();
                let sq_bc = bc.square_magnitude();
                let sq_ca = ca_vec.square_magnitude();
                let sq_conf = occt_core::precision::SQUARE_CONFUSION;

                if sq_ab >= sq_bc && sq_ab >= sq_ca && sq_ab > sq_conf {
                    let t = (GpVec::from_pnts(&pa, &p).dot(&ab) / sq_ab).clamp(0.0, 1.0);
                    u1 = ua + t * (ub - ua);
                    v1 = va + t * (vb - va);
                } else if sq_bc >= sq_ca && sq_bc > sq_conf {
                    let t = (GpVec::from_pnts(&pb, &p).dot(&bc) / sq_bc).clamp(0.0, 1.0);
                    u1 = ub + t * (uc - ub);
                    v1 = vb + t * (vc - vb);
                } else if sq_ca > sq_conf {
                    let t = (GpVec::from_pnts(&pc, &p).dot(&ca_vec) / sq_ca).clamp(0.0, 1.0);
                    u1 = uc + t * (ua - uc);
                    v1 = vc + t * (va - vc);
                } else {
                    u1 = ua;
                    v1 = va;
                }
            }
        }
        // `default: break` (`:835-837`).
        IntfPIType::External => {}
    }

    // `:839-848`: the curve parameter comes from the polygon segment the
    // interference engine reported.
    let (_, seg_index, param) = sp.info_first();
    let w = polygon.approx_param_on_curve(seg_index, param);
    (u1, v1, w)
}
