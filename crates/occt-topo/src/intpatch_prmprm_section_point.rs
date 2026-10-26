//! `SectionPointToParameters` (`IntPatch_PrmPrmIntersection.cxx:1584-1687`).
//!
//! Turns an `Intf_SectionPoint` of the polyhedron-polyhedron interference back
//! into UV parameters on both surfaces: a vertex is looked up, an edge is
//! interpolated, a face gets barycentric weights from the three triangles
//! corners.

use occt_core::gp::{GpPnt, GpVec};
use occt_core::intf::{IntfPIType, IntfSectionPoint};

use crate::intpatch_polyhedron::IntPatchPolyhedron;

/// `SectionPointToParameters(Sp, Poly1, Poly2, u1, v1, u2, v2)`
/// (`cxx:1584-1687`).
///
/// The four parameters are in-out: OCCT's `default:` arm of both `switch`es
/// leaves them untouched (`cxx:1636-1641`, `:1680-1685`), so a point whose
/// type is neither vertex, edge nor face keeps the value it had.
pub fn section_point_to_parameters(
    sp: &IntfSectionPoint,
    poly1: &IntPatchPolyhedron,
    poly2: &IntPatchPolyhedron,
    u1: &mut f64,
    v1: &mut f64,
    u2: &mut f64,
    v2: &mut f64,
) {
    let p = sp.pnt();
    parameters_on(sp.info_first_2(), poly1, p, u1, v1);
    parameters_on(sp.info_second_2(), poly2, p, u2, v2);
}

/// One `switch (typ)` of `cxx:1598-1641` / `:1645-1685`.
fn parameters_on(
    info: (IntfPIType, i32, i32, f64),
    poly: &IntPatchPolyhedron,
    p: &GpPnt,
    u: &mut f64,
    v: &mut f64,
) {
    let (typ, adr1, adr2, param) = info;
    match typ {
        // `Intf_VERTEX`: `Poly.Parameters(Adr1, u, v)`.
        IntfPIType::Vertex => {
            let (a, b) = poly.parameters(adr1 as usize);
            *u = a;
            *v = b;
        }
        // `Intf_EDGE`: interpolate along the edge with `Param`.
        IntfPIType::Edge => {
            let (u1, v1) = poly.parameters(adr1 as usize);
            let (u2, v2) = poly.parameters(adr2 as usize);
            *u = u1 + param * (u2 - u1);
            *v = v1 + param * (v2 - v1);
        }
        // `Intf_FACE`: barycentric weights of `P` inside the triangle.
        IntfPIType::Face => {
            let (pt1, pt2, pt3) = poly.triangle(adr1 as usize);
            let pa = poly.point(pt1);
            let pb = poly.point(pt2);
            let pc = poly.point(pt3);
            let (ua, va) = poly.parameters(pt1 as usize);
            let (ub, vb) = poly.parameters(pt2 as usize);
            let (uc, vc) = poly.parameters(pt3 as usize);
            let normale = GpVec::from_pnts(&pa, &pb).crossed(&GpVec::from_pnts(&pa, &pc));
            let cc = GpVec::from_pnts(&pa, &pb)
                .crossed(&GpVec::from_pnts(&pa, p))
                .dot(&normale);
            let ca = GpVec::from_pnts(&pb, &pc)
                .crossed(&GpVec::from_pnts(&pb, p))
                .dot(&normale);
            let cb = GpVec::from_pnts(&pc, &pa)
                .crossed(&GpVec::from_pnts(&pc, p))
                .dot(&normale);
            let cabc = ca + cb + cc;
            // `cxx:1629-1631` divides without a guard; keep it.
            let (ca, cb, cc) = (ca / cabc, cb / cabc, cc / cabc);
            *u = ca * ua + cb * ub + cc * uc;
            *v = ca * va + cb * vb + cc * vc;
        }
        // `default:` of both `switch`es: nothing is written.
        IntfPIType::External => {}
    }
}
