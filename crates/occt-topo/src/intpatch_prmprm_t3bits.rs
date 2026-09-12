//! `IntPatch_PrmPrmIntersection_T3Bits` and voxel fill (`Remplit*`).

const BASE: i32 = 128;
const BASE_M1: i32 = 127;
const DECAL: i32 = 7;
const DECAL2: i32 = 14;

pub(crate) struct T3Bits {
    p: Vec<i32>,
}

impl T3Bits {
    pub(crate) fn new(size: i32) -> Self {
        let nb = ((size as i64) * (size as i64) * (size as i64) / 32) as usize;
        Self { p: vec![0; nb] }
    }

    pub(crate) fn add(&mut self, t: i32) {
        let idx = (t >> 5) as usize;
        if idx < self.p.len() {
            self.p[idx] |= 1 << ((t as u32) & 31);
        }
    }

    pub(crate) fn val(&self, t: i32) -> i32 {
        let idx = (t >> 5) as usize;
        if idx >= self.p.len() {
            return 0;
        }
        self.p[idx] & (1 << ((t as u32) & 31))
    }

    pub(crate) fn raz(&mut self, t: i32) {
        let idx = (t >> 5) as usize;
        if idx < self.p.len() {
            self.p[idx] &= !(1 << ((t as u32) & 31));
        }
    }

    pub(crate) fn and_next(&mut self, oth: &mut Self, indice: &mut i32) -> bool {
        let mut k = (*indice >> 5) as usize;
        while k < self.p.len() && k < oth.p.len() {
            let mut r = self.p[k] & oth.p[k];
            if r != 0 {
                let mut c = 0u32;
                loop {
                    if r & 1 != 0 {
                        let op = ((k as i32) << 5) | (c as i32);
                        self.raz(op);
                        oth.raz(op);
                        *indice = op;
                        return true;
                    }
                    c += 1;
                    r >>= 1;
                    if c >= 32 {
                        break;
                    }
                }
            }
            k += 1;
        }
        false
    }
}

pub(crate) fn grille_integer(ix: i32, iy: i32, iz: i32) -> i32 {
    ix | (iy << DECAL) | (iz << DECAL2)
}

pub(crate) fn integer_grille(tt: i32) -> (i32, i32, i32) {
    let ix = tt & BASE_M1;
    let t = tt >> DECAL;
    let iy = t & BASE_M1;
    let iz = t >> DECAL;
    (ix, iy, iz)
}

pub(crate) fn dans_grille(t: i32) -> bool {
    t >= 0 && t < BASE
}

pub(crate) fn nb_points_grille() -> i32 {
    BASE
}

pub(crate) fn code_reject(
    x0: f64,
    y0: f64,
    z0: f64,
    x1: f64,
    y1: f64,
    z1: f64,
    x: f64,
    y: f64,
    z: f64,
) -> i32 {
    let mut code = 0;
    if x < x0 {
        code = 1;
    }
    if y < y0 {
        code |= 2;
    }
    if z < z0 {
        code |= 4;
    }
    if x > x1 {
        code |= 8;
    }
    if y > y1 {
        code |= 16;
    }
    if z > z1 {
        code |= 32;
    }
    code
}

pub(crate) fn remplit(a: i32, b: i32, c: i32, map: &mut T3Bits) {
    if a != -1 {
        map.add(a);
    }
    if b != -1 {
        map.add(b);
    }
    if c != -1 {
        map.add(c);
    }
    if a != -1 && b != -1 && c != -1 {
        let (iax, iay, iaz) = integer_grille(a);
        let (ibx, iby, ibz) = integer_grille(b);
        let (icx, icy, icz) = integer_grille(c);
        remplit_tri(iax, iay, iaz, ibx, iby, ibz, icx, icy, icz, map);
    }
}

fn remplit_lin(x1: i32, y1: i32, z1: i32, x2: i32, y2: i32, z2: i32, map: &mut T3Bits) {
    let xg = (x1 - x2).abs();
    let yg = (y1 - y2).abs();
    let zg = (z1 - z2).abs();
    if dans_grille(x1) && dans_grille(y1) && dans_grille(z1) {
        map.add(grille_integer(x1, y1, z1));
    }
    if xg <= 1 && yg <= 1 && zg <= 1 {
        return;
    }
    let mx = (x1 + x2) >> 1;
    let my = (y1 + y2) >> 1;
    let mz = (z1 + z2) >> 1;
    remplit_lin(x1, y1, z1, mx, my, mz, map);
    remplit_lin(x2, y2, z2, mx, my, mz, map);
}

fn remplit_tri(
    x1: i32,
    y1: i32,
    z1: i32,
    x2: i32,
    y2: i32,
    z2: i32,
    x3: i32,
    y3: i32,
    z3: i32,
    map: &mut T3Bits,
) {
    if x1 == x2 && x1 == x3 && y1 == y2 && y1 == y3 && z1 == z2 && z1 == z3 {
        if dans_grille(x1) && dans_grille(y1) && dans_grille(z1) {
            map.add(grille_integer(x1, y1, z1));
        }
        return;
    }
    let xg = (x1 + x2 + x3) / 3;
    let yg = (y1 + y2 + y3) / 3;
    let zg = (z1 + z2 + z3) / 3;
    if xg == x1 && yg == y1 && zg == z1 {
        remplit_lin(x1, y1, z1, x2, y2, z2, map);
        remplit_lin(x1, y1, z1, x3, y3, z3, map);
        return;
    }
    if xg == x2 && yg == y2 && zg == z2 {
        remplit_lin(x2, y2, z2, x1, y1, z1, map);
        remplit_lin(x2, y2, z2, x3, y3, z3, map);
        return;
    }
    if xg == x3 && yg == y3 && zg == z3 {
        remplit_lin(x3, y3, z3, x2, y2, z2, map);
        remplit_lin(x3, y3, z3, x1, y1, z1, map);
        return;
    }
    if dans_grille(xg) && dans_grille(yg) && dans_grille(zg) {
        map.add(grille_integer(xg, yg, zg));
    }
    if xg != x3 || yg != y3 || zg != z3 {
        remplit_tri(x1, y1, z1, x2, y2, z2, xg, yg, zg, map);
    }
    if xg != x1 || yg != y1 || zg != z1 {
        remplit_tri(xg, yg, zg, x2, y2, z2, x3, y3, z3, map);
    }
    if xg != x2 || yg != y2 || zg != z2 {
        remplit_tri(x1, y1, z1, xg, yg, zg, x3, y3, z3, map);
    }
}
