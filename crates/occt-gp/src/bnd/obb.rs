//! Oriented Bounding Box. Source: `Bnd_OBB.hxx`
use crate::gp::{xyz::GpXyz, pnt::GpPnt, dir::GpDir, ax3::GpAx3};

/// Oriented bounding box: center + 3 orthogonal axes + 3 half-dimensions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BndOBB {
    center: GpXyz,
    axes: [GpXyz; 3],  // X, Y, Z directions (unit vectors)
    hdims: [f64; 3],    // half-sizes along X, Y, Z
    is_aabox: bool,
}

impl BndOBB {
    /// Empty (invalid) OBB. Source: `Bnd_OBB.hxx:54`
    pub fn new() -> Self {
        Self { center: GpXyz::zero(), axes: [GpXyz::zero(); 3], hdims: [-1.0; 3], is_aabox: false }
    }

    /// Constructor with all parameters. Source: `Bnd_OBB.hxx:61`
    pub fn from_params(center: &GpPnt, xdir: &GpDir, ydir: &GpDir, zdir: &GpDir, hx: f64, hy: f64, hz: f64) -> Self {
        Self {
            center: center.coord,
            axes: [*xdir.xyz(), *ydir.xyz(), *zdir.xyz()],
            hdims: [hx, hy, hz],
            is_aabox: false,
        }
    }

    /// From axis-aligned box. Source: `Bnd_OBB.hxx` (from Bnd_Box)
    pub fn from_aabb(box3d: &super::box3d::BndBox) -> Self {
        if let Some((xmin, xmax, ymin, ymax, zmin, zmax)) = box3d.get() {
            let cx = (xmin + xmax) * 0.5;
            let cy = (ymin + ymax) * 0.5;
            let cz = (zmin + zmax) * 0.5;
            let hx = (xmax - xmin) * 0.5;
            let hy = (ymax - ymin) * 0.5;
            let hz = (zmax - zmin) * 0.5;
            Self {
                center: GpXyz::new(cx, cy, cz),
                axes: [GpXyz::new(1.,0.,0.), GpXyz::new(0.,1.,0.), GpXyz::new(0.,0.,1.)],
                hdims: [hx, hy, hz],
                is_aabox: true,
            }
        } else {
            Self::new()
        }
    }

    pub fn is_void(&self) -> bool { self.hdims[0] < 0.0 || self.is_aabox && self.hdims[0] <= 0.0 }
    pub fn is_aabox(&self) -> bool { self.is_aabox }
    pub fn center(&self) -> GpPnt { GpPnt::from_xyz(&self.center) }
    pub fn x_direction(&self) -> GpDir { GpDir::from_xyz(&self.axes[0]).unwrap_or_default() }
    pub fn y_direction(&self) -> GpDir { GpDir::from_xyz(&self.axes[1]).unwrap_or_default() }
    pub fn z_direction(&self) -> GpDir { GpDir::from_xyz(&self.axes[2]).unwrap_or_default() }
    pub fn half_sizes(&self) -> (f64, f64, f64) { (self.hdims[0], self.hdims[1], self.hdims[2]) }

    /// Transform point to local OBB coordinates. Source: `Bnd_OBB.hxx` (IsOut)
    fn to_local(&self, p: &GpPnt) -> GpXyz {
        let d = p.coord.subtracted(&self.center);
        GpXyz::new(d.dot(&self.axes[0]), d.dot(&self.axes[1]), d.dot(&self.axes[2]))
    }

    /// Is point outside OBB?
    pub fn is_out_point(&self, p: &GpPnt) -> bool {
        let local = self.to_local(p);
        local.x.abs() > self.hdims[0] || local.y.abs() > self.hdims[1] || local.z.abs() > self.hdims[2]
    }

    /// Is OBB completely inside another OBB?
    pub fn is_completely_inside(&self, other: &Self) -> bool {
        // Check 8 corners of self are inside other
        for sx in [-1.0, 1.0].iter() {
        for sy in [-1.0, 1.0].iter() {
        for sz in [-1.0, 1.0].iter() {
            let corner = self.center
                .added(&self.axes[0].multiplied(*sx * self.hdims[0]))
                .added(&self.axes[1].multiplied(*sy * self.hdims[1]))
                .added(&self.axes[2].multiplied(*sz * self.hdims[2]));
            if other.is_out_point(&GpPnt::from_xyz(&corner)) { return false; }
        }}}
        true
    }

    /// Position as Ax3. Source: `Bnd_OBB.hxx` (Position)
    pub fn position(&self) -> GpAx3 {
        GpAx3::new(
            GpPnt::from_xyz(&self.center),
            GpDir::from_xyz(&self.axes[2]).unwrap_or_default(),
            &GpDir::from_xyz(&self.axes[0]).unwrap_or_default(),
        ).unwrap_or_default()
    }
}

impl Default for BndOBB { fn default() -> Self { Self::new() } }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn obb_point_test() {
        let obb = BndOBB::from_params(
            &GpPnt::new(0.,0.,0.),
            &GpDir::from_axis(crate::gp::dir::DirAxis::X),
            &GpDir::from_axis(crate::gp::dir::DirAxis::Y),
            &GpDir::from_axis(crate::gp::dir::DirAxis::Z),
            1.0, 2.0, 3.0,
        );
        assert!(!obb.is_out_point(&GpPnt::new(0.5, 0., 0.)));
        assert!(obb.is_out_point(&GpPnt::new(2.0, 0., 0.)));
    }

    #[test]
    fn obb_from_aabb() {
        let mut aabb = crate::bnd::box3d::BndBox::new();
        aabb.add_point(&GpPnt::new(-1., -2., -3.));
        aabb.add_point(&GpPnt::new(1., 2., 3.));
        let obb = BndOBB::from_aabb(&aabb);
        assert!(!obb.is_void());
        assert_eq!(obb.half_sizes(), (1.0, 2.0, 3.0));
    }
}
