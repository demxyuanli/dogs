//! `CommonBox()` and its marking tail (`IntPolyh_MaillageAffinage.cxx:503-762`).
//!
//! `CommonBox` intersects the two surface boxes; `mark_common` then grows the
//! flat box and marks every node of both samplings with the sides it is out of.

use occt_core::bnd::BndBox;

use super::MaillageAffinage;

impl MaillageAffinage<'_> {
    /// `CommonBox()` (`cxx:503-507`).
    pub fn common_box(&mut self) {
        let b1 = self.box1.clone();
        let b2 = self.box2.clone();
        let (x_min, y_min, z_min, x_max, y_max, z_max) = self.common_box_of(&b1, &b2);
        self.mark_common(x_min, y_min, z_min, x_max, y_max, z_max);
    }

    /// `CommonBox(B1, B2, XMin, YMin, ZMin, XMax, YMax, ZMax)` (`cxx:517-613`),
    /// the box part only.
    pub fn common_box_of(
        &self,
        _b1: &BndBox,
        _b2: &BndBox,
    ) -> (f64, f64, f64, f64, f64, f64) {
        let (x10, y10, z10, x11, y11, z11) =
            self.box1.get().unwrap_or((0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        let (x20, y20, z20, x21, y21, z21) =
            self.box2.get().unwrap_or((0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

        let mut x_min = 0.0;
        let mut y_min = 0.0;
        let mut z_min = 0.0;
        let mut x_max = 0.0;
        let mut y_max = 0.0;
        let mut z_max = 0.0;

        if (x10 > x21) || (x20 > x11) || (y10 > y21) || (y20 > y11) || (z10 > z21) || (z20 > z11) {
            // `cxx:538-540`: disjoint boxes keep the null box.
        } else {
            if x11 <= x21 {
                x_max = x11;
            } else if x21 <= x11 {
                x_max = x21;
            }
            if x20 <= x10 {
                x_min = x10;
            } else if x10 <= x20 {
                x_min = x20;
            }
            if y11 <= y21 {
                y_max = y11;
            } else if y21 <= y11 {
                y_max = y21;
            }
            if y20 <= y10 {
                y_min = y10;
            } else if y10 <= y20 {
                y_min = y20;
            }
            if z11 <= z21 {
                z_max = z11;
            } else if z21 <= z11 {
                z_max = z21;
            }
            if z20 <= z10 {
                z_min = z10;
            } else if z10 <= z20 {
                z_min = z20;
            }
        }

        (x_min, y_min, z_min, x_max, y_max, z_max)
    }

    /// The tail of `CommonBox` (`cxx:617-762`): extend the flat box and mark
    /// the points of both surfaces with the code of the sides they are out of.
    pub fn mark_common(
        &mut self,
        x_min: f64,
        y_min: f64,
        z_min: f64,
        x_max: f64,
        y_max: f64,
        z_max: f64,
    ) {
        let mut x = x_max - x_min;
        let mut y = y_max - y_min;
        let mut z = z_max - z_min;

        if x == 0.0 && y != 0.0 {
            x = y * 0.1;
        } else if x == 0.0 && z != 0.0 {
            x = z * 0.1;
        } else {
            x *= 0.1;
        }

        if y == 0.0 && x != 0.0 {
            y = x * 0.1;
        } else if y == 0.0 && z != 0.0 {
            y = z * 0.1;
        } else {
            y *= 0.1;
        }

        if z == 0.0 && x != 0.0 {
            z = x * 0.1;
        } else if z == 0.0 && y != 0.0 {
            z = y * 0.1;
        } else {
            z *= 0.1;
        }

        let x_min = x_min - x;
        let x_max = x_max + x;
        let y_min = y_min - y;
        let y_max = y_max + y;
        let z_min = z_min - z;
        let z_max = z_max + z;

        let fin_tp1 = self.points1.nb_items();
        for i in 0..fin_tp1 {
            let pt = self.points1.get_mut(i);
            let mut r = if pt.x < x_min {
                1
            } else if pt.x > x_max {
                2
            } else {
                0
            };
            if pt.y < y_min {
                r |= 4;
            } else if pt.y > y_max {
                r |= 8;
            }
            if pt.z < z_min {
                r |= 16;
            } else if pt.z > z_max {
                r |= 32;
            }
            pt.part_of_common = r;
        }

        let fin_tp2 = self.points2.nb_items();
        for i in 0..fin_tp2 {
            let pt = self.points2.get_mut(i);
            let mut r = if pt.x < x_min {
                1
            } else if pt.x > x_max {
                2
            } else {
                0
            };
            if pt.y < y_min {
                r |= 4;
            } else if pt.y > y_max {
                r |= 8;
            }
            if pt.z < z_min {
                r |= 16;
            } else if pt.z > z_max {
                r |= 32;
            }
            pt.part_of_common = r;
        }
    }
}
