//! Sorted bounding box for spatial indexing. Source: `Bnd_BoundSortBox.hxx`
use crate::gp::GpPnt;
use crate::bnd::box3d::BndBox;

/// Element stored in a sorted box grid.
#[derive(Debug, Clone, Copy)]
pub struct SortBoxElement {
    pub index: usize,
    pub box3d: BndBox,
}

/// Sorted box structure for spatial queries.
/// Divides space into a grid of boxes for fast overlap/intersection tests.
#[derive(Debug, Clone)]
pub struct BoundSortBox {
    elements: Vec<SortBoxElement>,
    total_box: BndBox,
    nx: usize, ny: usize, nz: usize,
    delta: [f64; 3],
    origin: [f64; 3],
}

impl BoundSortBox {
    /// Initialize with elements and grid dimensions.
    pub fn new(elements: Vec<SortBoxElement>, nx: usize, ny: usize, nz: usize) -> Self {
        let total_box = elements.iter().fold(BndBox::new(), |mut acc, e| { acc.add_box(&e.box3d); acc });
        let cmin = total_box.corner_min(); let cmax = total_box.corner_max();
        let (ox, oy, oz) = (cmin.x(), cmin.y(), cmin.z());
        let (mx, my, mz) = (cmax.x(), cmax.y(), cmax.z());
        let delta = [(mx - ox) / nx as f64, (my - oy) / ny as f64, (mz - oz) / nz as f64];
        let delta = [if delta[0] == 0.0 { 1.0 } else { delta[0] }, if delta[1] == 0.0 { 1.0 } else { delta[1] }, if delta[2] == 0.0 { 1.0 } else { delta[2] }];
        Self { elements, total_box, nx, ny, nz, delta, origin: [ox, oy, oz] }
    }

    /// Find elements that overlap with query box.
    pub fn find_overlapping(&self, query: &BndBox) -> Vec<usize> {
        let mut result = Vec::new();
        if self.elements.is_empty() || query.is_out_box(&self.total_box) { return result; }
        // Determine grid range that overlaps
        let (qxmin, qxmax, qymin, qymax, qzmin, qzmax) = match self.clip_range(query) {
            Some(v) => v, None => return result,
        };
        for elem in &self.elements {
            if elem.box3d.is_out_box(query) { continue; }
            // Quick grid check
            if self.in_range(elem, qxmin, qxmax, qymin, qymax, qzmin, qzmax) {
                result.push(elem.index);
            }
        }
        result
    }

    fn clip_range(&self, q: &BndBox) -> Option<(usize,usize,usize,usize,usize,usize)> {
        if q.is_out_box(&self.total_box) { return None; }
        let (xmin, xmax, ymin, ymax, zmin, zmax) = q.get()?;
        let ixmin = ((xmin - self.origin[0]) / self.delta[0]).max(0.0) as usize;
        let ixmin = ((xmin - self.origin[0]) / self.delta[0]).max(0.0) as usize;
        let ixmax = ((xmax - self.origin[0]) / self.delta[0]).min(self.nx as f64 - 1.0).max(0.0) as usize;
        let iymin = ((ymin - self.origin[1]) / self.delta[1]).max(0.0) as usize;
        let iymax = ((ymax - self.origin[1]) / self.delta[1]).min(self.ny as f64 - 1.0).max(0.0) as usize;
        let izmin = ((zmin - self.origin[2]) / self.delta[2]).max(0.0) as usize;
        let izmax = ((zmax - self.origin[2]) / self.delta[2]).min(self.nz as f64 - 1.0).max(0.0) as usize;
        Some((ixmin.min(self.nx-1), ixmax.min(self.nx-1), iymin.min(self.ny-1), iymax.min(self.ny-1), izmin.min(self.nz-1), izmax.min(self.nz-1)))
    }

    fn in_range(&self, elem: &SortBoxElement, xmin:usize,xmax:usize,ymin:usize,ymax:usize,zmin:usize,zmax:usize) -> bool {
        let c = elem.box3d.corner_min();
        let ix = ((c.x() - self.origin[0]) / self.delta[0]) as usize;
        let iy = ((c.y() - self.origin[1]) / self.delta[1]) as usize;
        let iz = ((c.z() - self.origin[2]) / self.delta[2]) as usize;
        ix <= xmax && (ix + 1) >= xmin && iy <= ymax && (iy + 1) >= ymin && iz <= zmax && (iz + 1) >= zmin
    }
}
