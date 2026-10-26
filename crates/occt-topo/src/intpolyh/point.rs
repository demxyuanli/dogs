//! `IntPolyh_Point` (`ModelingAlgorithms/TKGeomAlgo/IntPolyh/IntPolyh_Point.hxx:23-157`,
//! `.cxx:21-148`).
//!
//! A triangulation node with a 3D position and the UV parameters it was
//! sampled from, plus a "part of common" flag (`myPOC`) and the degenerated
//! flag the deflection code reads.

use occt_core::precision::COMPUTATIONAL;

/// `IntPolyh_Point` (`hxx:23`).
#[derive(Debug, Clone, Copy)]
pub struct IntPolyhPoint {
    /// `myX` / `myY` / `myZ`.
    pub x: f64,
    pub y: f64,
    pub z: f64,
    /// `myU` / `myV`: parameters on the surface.
    pub u: f64,
    pub v: f64,
    /// `myPOC` (`hxx:155`): 0 when the point is not common with the other
    /// surface.
    pub part_of_common: i32,
    /// `myDegenerated` (`hxx:156`).
    pub degenerated: bool,
}

impl Default for IntPolyhPoint {
    /// `IntPolyh_Point()` (`hxx:30-39`): `myPOC = 1`.
    fn default() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            u: 0.0,
            v: 0.0,
            part_of_common: 1,
            degenerated: false,
        }
    }
}

impl IntPolyhPoint {
    /// `IntPolyh_Point(x, y, z, u, v)` (`hxx:43-51`).
    pub fn new(x: f64, y: f64, z: f64, u: f64, v: f64) -> Self {
        Self {
            x,
            y,
            z,
            u,
            v,
            part_of_common: 1,
            degenerated: false,
        }
    }

    /// `Set(x, y, z, u, v, II = 1)` (`hxx:71-87`).
    pub fn set(&mut self, x: f64, y: f64, z: f64, u: f64, v: f64, ii: i32) {
        self.x = x;
        self.y = y;
        self.z = z;
        self.u = u;
        self.v = v;
        self.part_of_common = ii;
    }

    /// `Middle(MySurface, P1, P2)` (`cxx:23-36`): UV midpoint, 3D evaluated on
    /// the surface, so the new node lies on the surface and not on the chord.
    pub fn middle(&mut self, surface: &dyn occt_geom::Surface, p1: &Self, p2: &Self) {
        self.u = (p1.u + p2.u) * 0.5;
        self.v = (p1.v + p2.v) * 0.5;
        let p = surface.value(self.u, self.v);
        self.x = p.x();
        self.y = p.y();
        self.z = p.z();
    }

    /// `Add(P1)` (`cxx:39-49`).
    pub fn add(&self, p1: &Self) -> Self {
        let mut r = Self::default();
        r.x = self.x + p1.x;
        r.y = self.y + p1.y;
        r.z = self.z + p1.z;
        r.u = self.u + p1.u;
        r.v = self.v + p1.v;
        r
    }

    /// `Sub(P1)` (`cxx:53-64`).
    pub fn sub(&self, p1: &Self) -> Self {
        let mut r = Self::default();
        r.x = self.x - p1.x;
        r.y = self.y - p1.y;
        r.z = self.z - p1.z;
        r.u = self.u - p1.u;
        r.v = self.v - p1.v;
        r
    }

    /// `Divide(RR)` (`cxx:68-79`): a near-zero divisor leaves the default
    /// zero point, exactly as the `.cxx` does.
    pub fn divide(&self, rr: f64) -> Self {
        let mut r = Self::default();
        if rr.abs() > COMPUTATIONAL {
            r.x = self.x / rr;
            r.y = self.y / rr;
            r.z = self.z / rr;
            r.u = self.u / rr;
            r.v = self.v / rr;
        }
        r
    }

    /// `Multiplication(RR)` (`cxx:84-94`).
    pub fn multiplication(&self, rr: f64) -> Self {
        let mut r = Self::default();
        r.x = self.x * rr;
        r.y = self.y * rr;
        r.z = self.z * rr;
        r.u = self.u * rr;
        r.v = self.v * rr;
        r
    }

    /// `SquareModulus()` (`cxx:98-101`).
    pub fn square_modulus(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// `SquareDistance(P2)` (`cxx:105-110`).
    pub fn square_distance(&self, p2: &Self) -> f64 {
        (self.x - p2.x) * (self.x - p2.x)
            + (self.y - p2.y) * (self.y - p2.y)
            + (self.z - p2.z) * (self.z - p2.z)
    }

    /// `Dot(b)` (`cxx:114-118`).
    pub fn dot(&self, b: &Self) -> f64 {
        self.x * b.x + self.y * b.y + self.z * b.z
    }

    /// `Cross(a, b)` (`cxx:122-127`): writes the cross product into `self`;
    /// the UV slots are left untouched, as in the `.cxx`.
    pub fn cross(&mut self, a: &Self, b: &Self) {
        self.x = a.y * b.z - a.z * b.y;
        self.y = a.z * b.x - a.x * b.z;
        self.z = a.x * b.y - a.y * b.x;
    }

    /// 3D position as an `occt_core` point.
    pub fn pnt(&self) -> occt_core::gp::GpPnt {
        occt_core::gp::GpPnt::new(self.x, self.y, self.z)
    }

    /// `operator+(P1)` (`hxx:112`).
    pub fn plus(&self, p1: &Self) -> Self {
        self.add(p1)
    }

    /// `operator-(P1)` (`hxx:117`).
    pub fn minus(&self, p1: &Self) -> Self {
        self.sub(p1)
    }

    /// `operator/(rr)` (`hxx:122`).
    pub fn divided(&self, rr: f64) -> Self {
        self.divide(rr)
    }
}
