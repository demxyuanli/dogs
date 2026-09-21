use super::prelude::*;
use super::*;

/// Linear properties of a polyline (`BRepGProp_MeshCinert::Perform`).
pub fn mesh_cinert(nodes: &[GpPnt], loc: &GpPnt) -> GProps {
    let order = 2;
    let (gp, gw) = gauss_legendre(-1.0, 1.0, order);
    let mut inert = Inertia::default();
    for i in 0..nodes.len().saturating_sub(1) {
        let p1 = &nodes[i];
        let p2 = &nodes[i + 1];
        let dir = p2.coord.subtracted(&p1.coord);
        let upper = dir.modulus();
        if upper < 1e-12 {
            continue;
        }
        let d = dir.divided(upper);
        let um = 0.5 * upper;
        for j in 0..order {
            let u = um + um * gp[j];
            let p = GpPnt::from_xyz(&p1.coord.added(&d.multiplied(u)));
            let ds = gw[j];
            let (x, y, z) = (p.x() - loc.x(), p.y() - loc.y(), p.z() - loc.z());
            inert.mass += ds;
            inert.ix += x * ds;
            inert.iy += y * ds;
            inert.iz += z * ds;
            inert.ixy += x * y * ds;
            inert.iyz += y * z * ds;
            inert.ixz += x * z * ds;
            inert.ixx += (y * y + z * z) * ds;
            inert.iyy += (x * x + z * z) * ds;
            inert.izz += (x * x + y * y) * ds;
        }
        inert.mul(um);
    }
    let (mass, g, mat) = convert_s(&inert);
    GProps { dim: mass, loc: *loc, g, inertia: mat }
}
