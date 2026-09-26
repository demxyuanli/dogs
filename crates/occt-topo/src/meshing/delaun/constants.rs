use super::prelude::*;


pub(super) const ANG_DEV_1DEG: f64 = PI / 180.0;
pub(super) const ANG_DEV_90DEG: f64 = 90.0 * ANG_DEV_1DEG;
pub(super) const ANGLE_2PI: f64 = 2.0 * PI;

pub(super) const PREC: f64 = PCONFUSION;
pub(super) const PREC2: f64 = PREC * PREC;

/// Minimal 2-D axis-aligned bounding box. Source: `Bnd_B2d`.
#[derive(Debug, Clone, Copy)]
pub struct BndB2 {
    pub(super) min_x: f64,
    pub(super) min_y: f64,
    pub(super) max_x: f64,
    pub(super) max_y: f64,
    pub(super) is_void: bool,
}

impl BndB2 {
    pub(super) fn void() -> Self {
        Self {
            min_x: f64::INFINITY,
            min_y: f64::INFINITY,
            max_x: f64::NEG_INFINITY,
            max_y: f64::NEG_INFINITY,
            is_void: true,
        }
    }

    pub(super) fn add_pnt(&mut self, p: GpPnt2d) {
        if self.is_void {
            self.min_x = p.x();
            self.min_y = p.y();
            self.max_x = p.x();
            self.max_y = p.y();
            self.is_void = false;
        } else {
            self.min_x = self.min_x.min(p.x());
            self.max_x = self.max_x.max(p.x());
            self.min_y = self.min_y.min(p.y());
            self.max_y = self.max_y.max(p.y());
        }
    }

    pub(super) fn add_xy(&mut self, p: GpXY) {
        self.add_pnt(GpPnt2d::from_xy(p));
    }

    pub(super) fn enlarge(&mut self, tol: f64) {
        if !self.is_void {
            self.min_x -= tol;
            self.min_y -= tol;
            self.max_x += tol;
            self.max_y += tol;
        }
    }

    pub(super) fn is_out(&self, other: &BndB2) -> bool {
        if self.is_void || other.is_void {
            return true;
        }
        self.max_x < other.min_x
            || other.max_x < self.min_x
            || self.max_y < other.min_y
            || other.max_y < self.min_y
    }

    pub(super) fn get(&self) -> (f64, f64, f64, f64) {
        (self.min_x, self.min_y, self.max_x, self.max_y)
    }
}

pub(super) fn update_bnd_box(p1: GpXY, p2: GpXY, b: &mut BndB2) {
    b.add_xy(p1);
    b.add_xy(p2);
    b.enlarge(PREC);
}

/// Cell-filtered store of circumcircles keyed by triangle id.
///
/// Source: `BRepMesh_CircleTool` + `BRepMesh_CircleInspector`. Circles are
/// stored in a `HashMap` grid keyed by 2-D cells; a query point only inspects
/// the circles whose clamped bounding box covers the query's cell.
pub struct CircleTool {
    pub(super) tolerance: f64,
    pub(super) sq_tolerance: f64,
    pub(super) cell_size: GpXY,
    pub(super) face_min: GpXY,
    pub(super) face_max: GpXY,
    pub(super) circles: Vec<DelaunCircle>,
    pub(super) grid: HashMap<(i64, i64), Vec<i32>>,
}

impl CircleTool {
    pub(super) fn new() -> Self {
        Self {
            tolerance: PREC,
            sq_tolerance: PREC * PREC,
            cell_size: GpXY::new(10.0, 10.0),
            face_min: GpXY::zero(),
            face_max: GpXY::zero(),
            circles: Vec::new(),
            grid: HashMap::new(),
        }
    }

    pub(super) fn set_min_max_size(&mut self, min: GpXY, max: GpXY) {
        self.face_min = min;
        self.face_max = max;
    }

    pub(super) fn set_cell_size(&mut self, size_x: f64, size_y: f64) {
        self.cell_size = GpXY::new(if size_x > 0.0 { size_x } else { 1.0 }, if size_y > 0.0 { size_y } else { 1.0 });
        self.grid.clear();
    }

    /// `NCollection_CellFilter::Cell` (`NCollection_CellFilter.hxx:252-266`):
    /// `index = (int)(coord / cellSize)` with origin at 0, truncation toward zero.
    /// `myFaceMin`/`myFaceMax` only clamp the AABB in `bind`, they do not shift cells.
    pub(super) fn cell_of(&self, p: GpXY) -> (i64, i64) {
        let cx = (p.x / self.cell_size.x) as i64;
        let cy = (p.y / self.cell_size.y) as i64;
        (cx, cy)
    }

    pub(super) fn bind(&mut self, index: i32, location: GpXY, radius: f64) {
        if self.circles.len() < index as usize {
            self.circles.resize(index as usize, DelaunCircle::default());
        }
        self.circles[(index - 1) as usize] = DelaunCircle::with_radius(GpPnt2d::from_xy(location), radius);

        let min_x = (location.x - radius).max(self.face_min.x);
        let max_x = (location.x + radius).min(self.face_max.x);
        let min_y = (location.y - radius).max(self.face_min.y);
        let max_y = (location.y + radius).min(self.face_max.y);
        let min_cell = self.cell_of(GpXY::new(min_x, min_y));
        let max_cell = self.cell_of(GpXY::new(max_x, max_y));
        for ci in min_cell.0..=max_cell.0 {
            for cj in min_cell.1..=max_cell.1 {
                // CellFilter prepends (`NCollection_CellFilter.hxx:356-359`);
                // push then Select `.rev()` matches newest-first Inspect.
                self.grid.entry((ci, cj)).or_default().push(index);
            }
        }
    }

    /// Computes the circumcircle of three points. Source: `MakeCircle`.
    pub(super) fn make_circle(p1: GpXY, p2: GpXY, p3: GpXY) -> Option<(GpXY, f64)> {
        let sq_prec = PREC2;
        let link1 = GpXY::new(p3.x - p2.x, p2.y - p3.y);
        if link1.square_modulus() < sq_prec {
            return None;
        }
        let link2 = GpXY::new(p1.x - p3.x, p3.y - p1.y);
        if link2.square_modulus() < sq_prec {
            return None;
        }
        let link3 = GpXY::new(p2.x - p1.x, p1.y - p2.y);
        if link3.square_modulus() < sq_prec {
            return None;
        }
        let d = 2.0 * (p1.x * link1.y + p2.x * link2.y + p3.x * link3.y);
        // `BRepMesh_CircleTool.cxx:112`: `if (std::abs(aD) < gp::Resolution())`
        // with `gp::Resolution()` = `RealSmall()` = `DBL_MIN` (`gp.hxx:60`).
        if d.abs() < REAL_SMALL {
            return None;
        }
        let inv_d = 1.0 / d;
        let sq1 = p1.square_modulus();
        let sq2 = p2.square_modulus();
        let sq3 = p3.square_modulus();
        let loc = GpXY::new(
            (sq1 * link1.y + sq2 * link2.y + sq3 * link3.y) * inv_d,
            (sq1 * link1.x + sq2 * link2.x + sq3 * link3.x) * inv_d,
        );
        let r_sq = (p1.subtracted(&loc).square_modulus())
            .max(p2.subtracted(&loc).square_modulus())
            .max(p3.subtracted(&loc).square_modulus());
        let r = r_sq.sqrt() + 2.0 * f64::EPSILON;
        Some((loc, r))
    }

    /// Binds a circumcircle to the triangle index; returns `false` when the
    /// points are degenerate (no circle can be built). Source: `Bind`.
    pub(super) fn bind_circle(&mut self, index: i32, p1: GpXY, p2: GpXY, p3: GpXY) -> bool {
        match Self::make_circle(p1, p2, p3) {
            Some((loc, r)) => {
                self.bind(index, loc, r);
                true
            }
            None => false,
        }
    }

    /// Binds an implicit zero (invalid) circle. Source: `MocBind`.
    pub(super) fn moc_bind(&mut self, index: i32) {
        if self.circles.len() < index as usize {
            self.circles.resize(index as usize, DelaunCircle::default());
        }
        self.circles[(index - 1) as usize] = DelaunCircle::with_radius(GpPnt2d::zero(), -1.0);
    }

    /// Deletes the circle with the given index. Source: `Delete`.
    pub(super) fn delete(&mut self, index: i32) {
        if let Some(c) = self.circles.get_mut((index - 1) as usize) {
            if c.is_created {
                c.set_radius(-1.0);
            }
        }
    }

    /// Returns indices of all circles shot by the point (containing it within
    /// the circle tolerance). Source: `BRepMesh_CircleTool::Select` (`cxx:167-172`)
    /// + `BRepMesh_CircleInspector::Inspect` (`hxx:83-115`).
    ///
    /// CellFilter Inspect walks the cell list newest-first; deleted circles
    /// (`radius < 0`) are skipped (`CellFilter_Purge`).
    pub fn select(&self, point: GpXY) -> Vec<i32> {
        let mut shot = Vec::new();
        let key = self.cell_of(point);
        if let Some(list) = self.grid.get(&key) {
            for &idx in list.iter().rev() {
                let circle = &self.circles[(idx - 1) as usize];
                if !circle.is_created {
                    continue;
                }
                let dx = point.x - circle.center.x();
                let dy = point.y - circle.center.y();
                if dx * dx + dy * dy - circle.radius_sq <= self.sq_tolerance {
                    shot.push(idx);
                }
            }
        }
        shot
    }
}

/// Replacement mode for `create_and_replace_polygon_link`. Source: `ReplaceFlag`.
#[derive(Clone, Copy)]
pub(super) enum ReplaceFlag {
    Replace,
    InsertAfter,
    InsertBefore,
}

/// Stack of element ranges used by `cleanup_polygon`. Source: `StackOfFrames`.
pub(super) struct StackOfFrames {
    pub(super) frames: Vec<(usize, usize)>,
}

impl StackOfFrames {
    pub(super) fn new() -> Self {
        Self { frames: Vec::new() }
    }

    pub(super) fn push_frame(&mut self, start: usize, end: usize) {
        self.frames.push((start, end));
    }

    pub(super) fn pop_element(&mut self) -> usize {
        let (cur, end) = self.frames.last().copied().expect("pop on empty stack");
        self.frames.last_mut().unwrap().0 = cur + 1;
        if cur + 1 == end {
            self.frames.pop();
        }
        cur
    }

    pub(super) fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }
}

pub(super) fn map_bind(map: &mut BTreeMap<i32, bool>, key: i32, value: bool) -> bool {
    match map.entry(key) {
        std::collections::btree_map::Entry::Occupied(_) => false,
        std::collections::btree_map::Entry::Vacant(v) => {
            v.insert(value);
            true
        }
    }
}

/// Delaunay triangulator over a [`DelaunDataStructure`].
///
/// Mirrors `BRepMesh_Delaun`. The full Watson pipeline is ported: super
/// triangle, incremental insertion with the circle cell filter, constraint
/// (frontier / fixed) edge insertion via polygon meshing, frontier adjustment,
/// cleanup and auxiliary-element removal.
pub struct Delaun {
    pub(super) mesh_data: DelaunDataStructure,
    pub(super) circles: CircleTool,
    pub(super) sup_vert: Vec<i32>,
    pub(super) init_circles: bool,
    pub(super) sup_trian: DelaunTriangle,
    /// Set when `addTriangle` hit OCCT's `Standard_OutOfRange` condition (a link
    /// already carrying two triangles). OCCT lets the exception escape to
    /// `BRepMesh_BaseMeshAlgo.cxx:52-62`, which swallows it and leaves the face
    /// unmeshed; callers here read this flag and abort the polygon.
    pub(super) failed: bool,
}
