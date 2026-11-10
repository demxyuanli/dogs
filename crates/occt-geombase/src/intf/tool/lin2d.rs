//! `Intf_Tool::Lin2dBox` (`Intf_Tool.cxx:52-206`).

use occt_foundation::precision::Precision;
use occt_gp::bnd::BndBox2d;
use occt_gp::gp::GpLin2d;

use super::helpers::set2d_pnt_dir;
use super::structure::IntfTool;

impl IntfTool {
    /// `Intf_Tool::Lin2dBox` (`Intf_Tool.cxx:52-206`).
    pub fn lin2d_box(&mut self, l2d: &GpLin2d, domain: &BndBox2d, box_lin: &mut BndBox2d) {
        self.nb_seg = 0;
        box_lin.set_void();
        if domain.is_whole() {
            set2d_pnt_dir(box_lin, &l2d.location(), l2d.direction());
            let rev = l2d.direction().reversed();
            box_lin.add_dir(&rev);
            self.nb_seg = 1;
            self.begin_on_curve[0] = -Precision::INFINITE;
            self.end_on_curve[0] = Precision::INFINITE;
            return;
        } else if domain.is_void() {
            return;
        }

        let Some((xmin, ymin, xmax, ymax)) = domain.get() else {
            return;
        };
        let loc_x = l2d.location().x();
        let loc_y = l2d.location().y();
        let dir_x = l2d.direction().x;
        let dir_y = l2d.direction().y;

        let mut parmin = -Precision::INFINITE;
        let mut parmax = Precision::INFINITE;
        let mut xmin_set = 0.0f64;
        let mut xmax_set = 0.0f64;
        let mut ymin_set = 0.0f64;
        let mut ymax_set = 0.0f64;
        let x_to_set: bool;
        let y_to_set: bool;

        if dir_x > 0.0 {
            if domain.is_open_xmin() {
                parmin = -Precision::INFINITE;
            } else {
                parmin = (xmin - loc_x) / dir_x;
            }
            if domain.is_open_xmax() {
                parmax = Precision::INFINITE;
            } else {
                parmax = (xmax - loc_x) / dir_x;
            }
            x_to_set = true;
        } else if dir_x < 0.0 {
            if domain.is_open_xmax() {
                parmin = -Precision::INFINITE;
            } else {
                parmin = (xmax - loc_x) / dir_x;
            }
            if domain.is_open_xmin() {
                parmax = Precision::INFINITE;
            } else {
                parmax = (xmin - loc_x) / dir_x;
            }
            x_to_set = true;
        } else {
            // Parallel to axis X
            if loc_x < xmin || xmax < loc_x {
                return;
            }
            xmin_set = loc_x;
            xmax_set = loc_x;
            x_to_set = false;
        }

        if dir_y > 0.0 {
            let parcur = if domain.is_open_ymin() {
                -Precision::INFINITE
            } else {
                (ymin - loc_y) / dir_y
            };
            parmin = parmin.max(parcur);
            let parcur = if domain.is_open_ymax() {
                Precision::INFINITE
            } else {
                (ymax - loc_y) / dir_y
            };
            parmax = parmax.min(parcur);
            y_to_set = true;
        } else if dir_y < 0.0 {
            let parcur = if domain.is_open_ymax() {
                -Precision::INFINITE
            } else {
                (ymax - loc_y) / dir_y
            };
            parmin = parmin.max(parcur);
            let parcur = if domain.is_open_ymin() {
                Precision::INFINITE
            } else {
                (ymin - loc_y) / dir_y
            };
            parmax = parmax.min(parcur);
            y_to_set = true;
        } else {
            // Parallel to axis Y
            if loc_y < ymin || ymax < loc_y {
                return;
            }
            ymin_set = loc_y;
            ymax_set = loc_y;
            y_to_set = false;
        }

        self.nb_seg += 1;
        self.begin_on_curve[0] = parmin;
        self.end_on_curve[0] = parmax;

        if x_to_set {
            let par1 = loc_x + parmin * dir_x;
            let par2 = loc_x + parmax * dir_x;
            xmin_set = par1.min(par2);
            xmax_set = par1.max(par2);
        }

        if y_to_set {
            let par1 = loc_y + parmin * dir_y;
            let par2 = loc_y + parmax * dir_y;
            ymin_set = par1.min(par2);
            ymax_set = par1.max(par2);
        }

        box_lin.update(xmin_set, ymin_set, xmax_set, ymax_set);
    }
}
