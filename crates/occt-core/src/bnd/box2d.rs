//! Axis-aligned 2D bounding box. Source: `Bnd_Box2d.hxx` / `Bnd_Box2d.cxx`.
use crate::gp::{GpDir2d, GpLin2d, GpPnt2d, GpTrsf2d, GpXY};

/// Open-direction Get() sentinel used by `Bnd_Box2d.cxx` (`THE_BND_PRECISION_INFINITE`).
const BND2_PRECISION_INFINITE: f64 = 1e100;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BndBox2d {
    xmin: f64, xmax: f64, ymin: f64, ymax: f64,
    gap: f64,
    flags: u8,
}

const V2_VOID: u8  = 0b0000_0001;
const V2_XO: u8    = 0b0000_0010;
const V2_X1: u8    = 0b0000_0100;
const V2_YO: u8    = 0b0000_1000;
const V2_Y1: u8    = 0b0001_0000;
const V2_WHOLE: u8 = V2_XO | V2_X1 | V2_YO | V2_Y1;

impl BndBox2d {
    pub fn new() -> Self { Self { xmin:0.,xmax:0.,ymin:0.,ymax:0.,gap:0.,flags:V2_VOID } }
    pub fn from_corners(min: &GpPnt2d, max: &GpPnt2d) -> Self {
        Self { xmin:min.x(),xmax:max.x(),ymin:min.y(),ymax:max.y(),gap:0.,flags:0 }
    }

    pub fn is_void(&self) -> bool { self.flags & V2_VOID != 0 }
    pub fn is_whole(&self) -> bool { (self.flags & V2_WHOLE) == V2_WHOLE }
    pub fn is_finite(&self) -> bool { self.flags & (V2_VOID | V2_WHOLE) == 0 }

    pub fn is_open_xmin(&self) -> bool { self.flags & V2_XO != 0 }
    pub fn is_open_xmax(&self) -> bool { self.flags & V2_X1 != 0 }
    pub fn is_open_ymin(&self) -> bool { self.flags & V2_YO != 0 }
    pub fn is_open_ymax(&self) -> bool { self.flags & V2_Y1 != 0 }

    pub fn open_xmin(&mut self) { self.flags |= V2_XO; }
    pub fn open_xmax(&mut self) { self.flags |= V2_X1; }
    pub fn open_ymin(&mut self) { self.flags |= V2_YO; }
    pub fn open_ymax(&mut self) { self.flags |= V2_Y1; }

    pub fn set_void(&mut self) { self.flags = V2_VOID; self.gap = 0.0; }
    pub fn set_whole(&mut self) { self.flags = V2_WHOLE; }
    pub fn set_gap(&mut self, g: f64) { self.gap = g.abs(); }
    pub fn gap(&self) -> f64 { self.gap }
    /// `Bnd_Box2d::Enlarge(theTol)` (`Bnd_Box2d.hxx`):
    /// `Gap = max(Gap, |theTol|)`.
    pub fn enlarge(&mut self, t: f64) { self.gap = self.gap.max(t.abs()); }

    /// `Bnd_Box2d::Update(xmin, ymin, xmax, ymax)`.
    pub fn update(&mut self, xmin: f64, ymin: f64, xmax: f64, ymax: f64) {
        if self.flags & V2_VOID != 0 {
            self.xmin = xmin;
            self.ymin = ymin;
            self.xmax = xmax;
            self.ymax = ymax;
            self.flags &= !V2_VOID;
            return;
        }
        if self.flags & V2_XO == 0 { self.xmin = self.xmin.min(xmin); }
        if self.flags & V2_X1 == 0 { self.xmax = self.xmax.max(xmax); }
        if self.flags & V2_YO == 0 { self.ymin = self.ymin.min(ymin); }
        if self.flags & V2_Y1 == 0 { self.ymax = self.ymax.max(ymax); }
    }

    /// `Bnd_Box2d::Update(X, Y)` — add a single point by coordinates.
    pub fn update_point(&mut self, x: f64, y: f64) {
        self.update(x, y, x, y);
    }

    /// `Bnd_Box2d::GetXMin` (open → `-THE_BND_PRECISION_INFINITE`).
    pub fn xmin(&self) -> f64 {
        if self.flags & V2_XO != 0 { -BND2_PRECISION_INFINITE } else { self.xmin - self.gap }
    }
    pub fn xmax(&self) -> f64 {
        if self.flags & V2_X1 != 0 { BND2_PRECISION_INFINITE } else { self.xmax + self.gap }
    }
    pub fn ymin(&self) -> f64 {
        if self.flags & V2_YO != 0 { -BND2_PRECISION_INFINITE } else { self.ymin - self.gap }
    }
    pub fn ymax(&self) -> f64 {
        if self.flags & V2_Y1 != 0 { BND2_PRECISION_INFINITE } else { self.ymax + self.gap }
    }

    /// `Bnd_Box2d::Get`. Returns `None` when the box is void (OCCT throws).
    pub fn get(&self) -> Option<(f64, f64, f64, f64)> {
        if self.is_void() { return None; }
        Some((self.xmin(), self.ymin(), self.xmax(), self.ymax()))
    }

    /// `Bnd_Box2d::Add(gp_Dir2d)`.
    pub fn add_dir(&mut self, d: &GpDir2d) {
        let eps = f64::EPSILON;
        if d.x() < -eps { self.open_xmin(); }
        else if d.x() > eps { self.open_xmax(); }
        if d.y() < -eps { self.open_ymin(); }
        else if d.y() > eps { self.open_ymax(); }
    }

    /// `Bnd_Box2d::IsOut(Bnd_Box2d)`.
    pub fn is_out_box(&self, other: &Self) -> bool {
        if self.flags == 0 && other.flags == 0 {
            let delta = other.gap + self.gap;
            if self.xmin - other.xmax > delta { return true; }
            if other.xmin - self.xmax > delta { return true; }
            if self.ymin - other.ymax > delta { return true; }
            if other.ymin - self.ymax > delta { return true; }
            return false;
        }
        if self.is_void() || other.is_void() { return true; }
        if self.is_whole() || other.is_whole() { return false; }
        let Some((oxmin, oymin, oxmax, oymax)) = other.get() else { return true; };
        if self.flags & V2_XO == 0 && oxmax < self.xmin - self.gap { return true; }
        if self.flags & V2_X1 == 0 && oxmin > self.xmax + self.gap { return true; }
        if self.flags & V2_YO == 0 && oymax < self.ymin - self.gap { return true; }
        if self.flags & V2_Y1 == 0 && oymin > self.ymax + self.gap { return true; }
        false
    }

    /// `Bnd_Box2d::IsOut(gp_Lin2d)`.
    pub fn is_out_lin(&self, line: &GpLin2d) -> bool {
        if self.is_whole() { return false; }
        if self.is_void() { return true; }
        let Some((xmin, ymin, xmax, ymax)) = self.get() else { return true; };
        let cx = 0.5 * (xmin + xmax);
        let cy = 0.5 * (ymin + ymax);
        let hx = (xmax - cx).abs();
        let hy = (ymax - cy).abs();
        let dx = line.direction().x();
        let dy = line.direction().y();
        let loc = line.location();
        let prod0 = dx * (cy - loc.y()) - dy * (cx - loc.x());
        let prod1 = dx * hy;
        let prod2 = dy * hx;
        prod0.abs() > prod1.abs() + prod2.abs()
    }

    /// `Bnd_Box2d::IsOut(P0, P1)` — segment vs box.
    pub fn is_out_segment(&self, p0: &GpPnt2d, p1: &GpPnt2d) -> bool {
        if self.is_whole() { return false; }
        if self.is_void() { return true; }
        let Some((xmin, ymin, xmax, ymax)) = self.get() else { return true; };
        let mut status = true;
        let seg = GpXY::new(p1.x() - p0.x(), p1.y() - p0.y());
        let cx = 0.5 * (xmin + xmax);
        let cy = 0.5 * (ymin + ymax);
        let hx = (xmax - cx).abs();
        let hy = (ymax - cy).abs();
        let center_minus_p0 = GpXY::new(cx - p0.x(), cy - p0.y());
        let prod0 = seg.crossed(&center_minus_p0);
        let prod1 = seg.x() * hy;
        let prod2 = seg.y() * hx;
        if prod0.abs() <= prod1.abs() + prod2.abs() {
            let hseg_x = 0.5 * seg.x();
            let hseg_y = 0.5 * seg.y();
            let mx = (p0.x() + hseg_x - cx).abs();
            let my = (p0.y() + hseg_y - cy).abs();
            status = mx > hx + hseg_x.abs() || my > hy + hseg_y.abs();
        }
        status
    }

    /// `Bnd_Box2d::Distance`.
    pub fn distance_box(&self, other: &Self) -> f64 {
        if self.is_void() || other.is_void() { return 0.0; }
        let Some((xmin1, ymin1, xmax1, ymax1)) = self.get() else { return 0.0; };
        let Some((xmin2, ymin2, xmax2, ymax2)) = other.get() else { return 0.0; };
        fn dist_axis(min1: f64, max1: f64, min2: f64, max2: f64) -> f64 {
            if min1 > max2 { let d = min1 - max2; return d * d; }
            if min2 > max1 { let d = min2 - max1; return d * d; }
            0.0
        }
        (dist_axis(xmin1, xmax1, xmin2, xmax2) + dist_axis(ymin1, ymax1, ymin2, ymax2)).sqrt()
    }

    pub fn add_point(&mut self, p: &GpPnt2d) {
        if self.is_whole() { return; }
        if self.is_void() { self.xmin=p.x(); self.xmax=p.x(); self.ymin=p.y(); self.ymax=p.y(); self.flags=0; return; }
        if self.flags & V2_XO != 0 || p.x() < self.xmin { self.xmin = p.x(); self.flags &= !V2_XO; }
        if self.flags & V2_X1 != 0 || p.x() > self.xmax { self.xmax = p.x(); self.flags &= !V2_X1; }
        if self.flags & V2_YO != 0 || p.y() < self.ymin { self.ymin = p.y(); self.flags &= !V2_YO; }
        if self.flags & V2_Y1 != 0 || p.y() > self.ymax { self.ymax = p.y(); self.flags &= !V2_Y1; }
    }

    pub fn add_box(&mut self, other: &Self) {
        if other.is_void() || self.is_whole() { return; }
        if other.is_whole() { self.set_whole(); return; }
        if self.is_void() { *self = *other; return; }
        if self.flags & V2_XO != 0 || other.flags & V2_XO != 0 { self.flags |= V2_XO; } else { self.xmin = self.xmin.min(other.xmin); }
        if self.flags & V2_X1 != 0 || other.flags & V2_X1 != 0 { self.flags |= V2_X1; } else { self.xmax = self.xmax.max(other.xmax); }
        if self.flags & V2_YO != 0 || other.flags & V2_YO != 0 { self.flags |= V2_YO; } else { self.ymin = self.ymin.min(other.ymin); }
        if self.flags & V2_Y1 != 0 || other.flags & V2_Y1 != 0 { self.flags |= V2_Y1; } else { self.ymax = self.ymax.max(other.ymax); }
        self.gap = self.gap.max(other.gap);
    }

    pub fn is_out(&self, p: &GpPnt2d) -> bool {
        if self.is_void() { return true; }
        if self.is_whole() { return false; }
        let g = self.gap;
        (self.flags & V2_XO == 0 && p.x() < self.xmin - g)
            || (self.flags & V2_X1 == 0 && p.x() > self.xmax + g)
            || (self.flags & V2_YO == 0 && p.y() < self.ymin - g)
            || (self.flags & V2_Y1 == 0 && p.y() > self.ymax + g)
    }

    pub fn corner_min(&self) -> GpPnt2d { GpPnt2d::new(self.xmin - self.gap, self.ymin - self.gap) }
    pub fn corner_max(&self) -> GpPnt2d { GpPnt2d::new(self.xmax + self.gap, self.ymax + self.gap) }

    /// `Bnd_Box2d::Center`.
    pub fn center(&self) -> Option<GpPnt2d> {
        if self.is_void() {
            return None;
        }
        Some(GpPnt2d::new(
            0.5 * (self.xmin() + self.xmax()),
            0.5 * (self.ymin() + self.ymax()),
        ))
    }

    /// `Bnd_Box2d::Transformed` (`Bnd_Box2d.cxx:155`).
    ///
    /// Identity: copy. Translation: shift closed sides. Otherwise: open
    /// directions become transformed `gp_Dir2d`s via `Add(Dir)`, remaining
    /// finite corners are transformed as points, gap is preserved.
    pub fn transformed(&self, t: &GpTrsf2d) -> Self {
        let a_f = t.form();
        let mut a_new = *self;
        if self.is_void() {
            return a_new;
        }
        if a_f == crate::gp::trsf_form::TrsfForm::Identity {
            return a_new;
        }
        if a_f == crate::gp::trsf_form::TrsfForm::Translation {
            let a_dx = t.translation_part().x();
            let a_dy = t.translation_part().y();
            if self.flags & V2_XO == 0 {
                a_new.xmin += a_dx;
            }
            if self.flags & V2_X1 == 0 {
                a_new.xmax += a_dx;
            }
            if self.flags & V2_YO == 0 {
                a_new.ymin += a_dy;
            }
            if self.flags & V2_Y1 == 0 {
                a_new.ymax += a_dy;
            }
            return a_new;
        }
        let dir_xmin = GpDir2d { x: -1.0, y: 0.0 };
        let dir_xmax = GpDir2d { x: 1.0, y: 0.0 };
        let dir_ymin = GpDir2d { x: 0.0, y: -1.0 };
        let dir_ymax = GpDir2d { x: 0.0, y: 1.0 };
        let mut a_vertex = [true; 4];
        let mut a_d: [GpDir2d; 6] = [dir_xmax; 6];
        let mut a_nb_dirs = 0usize;
        if self.flags & V2_XO != 0 {
            a_d[a_nb_dirs] = dir_xmin;
            a_nb_dirs += 1;
            a_vertex[0] = false;
            a_vertex[2] = false;
        }
        if self.flags & V2_X1 != 0 {
            a_d[a_nb_dirs] = dir_xmax;
            a_nb_dirs += 1;
            a_vertex[1] = false;
            a_vertex[3] = false;
        }
        if self.flags & V2_YO != 0 {
            a_d[a_nb_dirs] = dir_ymin;
            a_nb_dirs += 1;
            a_vertex[0] = false;
            a_vertex[1] = false;
        }
        if self.flags & V2_Y1 != 0 {
            a_d[a_nb_dirs] = dir_ymax;
            a_nb_dirs += 1;
            a_vertex[2] = false;
            a_vertex[3] = false;
        }
        a_new.set_void();
        for i in 0..a_nb_dirs {
            let d = a_d[i].transformed(t);
            a_new.add_dir(&d);
        }
        let a_p = [
            GpPnt2d::new(self.xmin, self.ymin),
            GpPnt2d::new(self.xmax, self.ymin),
            GpPnt2d::new(self.xmin, self.ymax),
            GpPnt2d::new(self.xmax, self.ymax),
        ];
        for i in 0..4 {
            if a_vertex[i] {
                a_new.add_point(&a_p[i].transformed(t));
            }
        }
        a_new.gap = self.gap;
        a_new
    }

    pub fn transform(&mut self, t: &GpTrsf2d) {
        *self = self.transformed(t);
    }
}

impl Default for BndBox2d { fn default() -> Self { Self::new() } }
