//! Mass and inertia properties of polyhedra. Source: `GProp_GProps`
use crate::gp::{GpPnt, GpXyz};

/// Mass properties: mass, center of mass, 3x3 inertia tensor.
#[derive(Debug, Clone)]
pub struct InertiaProps {
    pub mass: f64,
    pub center: GpPnt,
    pub inertia: [f64; 9], // row-major 3x3, about global origin
}

impl InertiaProps {
    pub fn new() -> Self { Self { mass: 0.0, center: GpPnt::zero(), inertia: [0.0; 9] } }
}

/// Compute mass properties of a triangle-mesh solid.
/// vertices: mesh vertices. tris: (i,j,k) index triples.
/// Uses tetrahedron decomposition about the origin.
pub fn compute_inertia(vertices: &[GpPnt], tris: &[(usize, usize, usize)], density: f64) -> InertiaProps {
    let mut props = InertiaProps::new();
    for &(i, j, k) in tris {
        let a = vertices[i].coord;
        let b = vertices[j].coord;
        let c = vertices[k].coord;
        let vol6 = a.dot_cross(&b, &c); // 6*volume (signed)
        let vol = vol6 / 6.0;
        let mass = vol * density;
        props.mass += mass;
        // Tetra centroid = (a+b+c)/4
        let centroid = a.added(&b).added(&c).multiplied(0.25);
        let cm = props.center.coord.multiplied(props.mass - mass);
        let cw = centroid.multiplied(mass);
        let new_total = props.mass;
        if new_total.abs() > 1e-30 {
            props.center = GpPnt::from_xyz(&cm.added(&cw).divided(new_total));
        }
        // Inertia: tetrahedron contribution about origin
        let a0 = a.added(&b).added(&c);
        let ab = a.added(&b);
        let bc = b.added(&c);
        let ca = c.added(&a);
        for m in 0..3 {
            for n in 0..3 {
                // Second moment integral over tetra: (vol/20) * Σ_{pairs} (monomial products)
                // Approx via 4-vertex symmetric sum:
                let va = [a.x, a.y, a.z]; let vb = [b.x, b.y, b.z];
                let vc = [c.x, c.y, c.z];
                let pairs = [
                    va[m]*va[n], va[m]*vb[n], va[m]*vc[n],
                    vb[m]*vb[n], vb[m]*vc[n], vc[m]*vc[n],
                ];
                let mut sum = 0.0;
                for p in &pairs { sum += p; }
                props.inertia[m*3+n] += density * vol6 / 120.0 * sum;
            }
        }
    }
    props
}

/// Shift inertia tensor to center of mass (parallel axis theorem).
/// Takes inertia about origin, returns inertia about center.
pub fn inertia_about_center(p: &InertiaProps) -> [f64; 9] {
    let cx = p.center.x(); let cy = p.center.y(); let cz = p.center.z();
    let m = p.mass;
    // I_cm = I_origin - m*(r²I - rrᵀ)
    let r2 = cx*cx + cy*cy + cz*cz;
    let mut out = p.inertia;
    out[0] -= m * (r2 - cx*cx);
    out[1] -= m * (-cx*cy);
    out[2] -= m * (-cx*cz);
    out[3] -= m * (-cy*cx);
    out[4] -= m * (r2 - cy*cy);
    out[5] -= m * (-cy*cz);
    out[6] -= m * (-cz*cx);
    out[7] -= m * (-cz*cy);
    out[8] -= m * (r2 - cz*cz);
    out
}

/// Jacobi eigenvalues of a symmetric 3x3 matrix.
pub fn inertia_eigenvalues(i: &[f64; 9]) -> [f64; 3] {
    let mut a = [[i[0], i[1], i[2]], [i[3], i[4], i[5]], [i[6], i[7], i[8]]];
    let mut eig = [a[0][0], a[1][1], a[2][2]];
    for _ in 0..50 {
        let mut p = 0usize; let mut q = 1usize; let mut max_off = 0.0;
        for x in 0..3 { for y in (x+1)..3 {
            if a[x][y].abs() > max_off { max_off = a[x][y].abs(); p = x; q = y; }
        }}
        if max_off < 1e-15 { break; }
        let theta = 0.5 * (eig[q] - eig[p]) / a[p][q];
        let t = 1.0 / (theta.abs() + (1.0 + theta*theta).sqrt());
        let t = if theta < 0.0 { -t } else { t };
        let c = 1.0 / (1.0 + t*t).sqrt();
        let s = t * c;
        let tau = s / (1.0 + c);
        let h = t * a[p][q];
        eig[p] -= h; eig[q] += h;
        a[p][q] = 0.0;
        for j in 0..3 {
            if j != p && j != q {
                let g = a[p][j]; let h = a[q][j];
                a[p][j] = g - s*(h + g*tau);
                a[q][j] = h + s*(g - h*tau);
            }
        }
    }
    eig.sort_by(|x, y| y.partial_cmp(x).unwrap());
    eig
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cube_mass_one() {
        // Unit cube, 12 triangles (2 per face), density 1 → mass 1
        let vs = [
            GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(1.,1.,0.), GpPnt::new(0.,1.,0.),
            GpPnt::new(0.,0.,1.), GpPnt::new(1.,0.,1.), GpPnt::new(1.,1.,1.), GpPnt::new(0.,1.,1.),
        ];
        let tris = [
            (0,1,5),(0,5,4),(1,2,6),(1,6,5),(2,3,7),(2,7,6),
            (3,0,4),(3,4,7),(0,3,2),(0,2,1),(4,5,6),(4,6,7),
        ];
        let p = compute_inertia(&vs, &tris, 1.0);
        assert!((p.mass - 1.0).abs() < 1e-6, "mass {}", p.mass);
        assert!((p.center.x() - 0.5).abs() < 1e-3);
    }

    #[test]
    fn tetra_volume() {
        let vs = [GpPnt::new(0.,0.,0.), GpPnt::new(1.,0.,0.), GpPnt::new(0.,1.,0.), GpPnt::new(0.,0.,1.)];
        let tris = [(0,1,2),(0,1,3),(0,2,3),(1,2,3)];
        let p = compute_inertia(&vs, &tris, 1.0);
        // Unit right tetra volume = 1/6
        assert!((p.mass - 1.0/6.0).abs() < 1e-8, "mass {}", p.mass);
    }
}
