//! `Intf_Tool::LinBox` (`Intf_Tool.cxx:746-963`).

use occt_foundation::precision::Precision;
use occt_gp::bnd::BndBox;
use occt_gp::gp::GpLin;

use super::helpers::add3d_dir;
use super::structure::IntfTool;

impl IntfTool {
    /// `Intf_Tool::LinBox` (`Intf_Tool.cxx:746-963`).
    pub fn lin_box(&mut self, l: &GpLin, domain: &BndBox, box_lin: &mut BndBox) {
        self.nb_seg = 0;
        box_lin.set_void();
        if domain.is_whole() {
            // boxLin.Set(L.Location(), L.Direction()) + Add(Direction().Reversed())
            box_lin.set_void();
            box_lin.add_point(&l.location());
            let dir = l.direction();
            add3d_dir(box_lin, &dir);
            let rev = dir.reversed();
            add3d_dir(box_lin, &rev);
            self.nb_seg = 1;
            self.begin_on_curve[0] = -Precision::INFINITE;
            self.end_on_curve[0] = Precision::INFINITE;
            return;
        } else if domain.is_void() {
            return;
        }

        let Some((xmin, xmax, ymin, ymax, zmin, zmax)) = domain.get() else {
            return;
        };
        let loc = l.location();
        let dir = l.direction();

        let mut parmin = -Precision::INFINITE;
        let mut parmax = Precision::INFINITE;
        let mut xmin_set = 0.0f64;
        let mut xmax_set = 0.0f64;
        let mut ymin_set = 0.0f64;
        let mut ymax_set = 0.0f64;
        let mut zmin_set = 0.0f64;
        let mut zmax_set = 0.0f64;
        let x_to_set: bool;
        let y_to_set: bool;
        let z_to_set: bool;

        if dir.x() > 0.0 {
            parmin = if domain.is_open_xmin() {
                -Precision::INFINITE
            } else {
                (xmin - loc.x()) / dir.x()
            };
            parmax = if domain.is_open_xmax() {
                Precision::INFINITE
            } else {
                (xmax - loc.x()) / dir.x()
            };
            x_to_set = true;
        } else if dir.x() < 0.0 {
            parmin = if domain.is_open_xmax() {
                -Precision::INFINITE
            } else {
                (xmax - loc.x()) / dir.x()
            };
            parmax = if domain.is_open_xmin() {
                Precision::INFINITE
            } else {
                (xmin - loc.x()) / dir.x()
            };
            x_to_set = true;
        } else {
            // Perpendicular to axis X
            if loc.x() < xmin || xmax < loc.x() {
                return;
            }
            xmin_set = loc.x();
            xmax_set = loc.x();
            x_to_set = false;
        }

        if dir.y() > 0.0 {
            let parcur = if domain.is_open_ymin() {
                -Precision::INFINITE
            } else {
                (ymin - loc.y()) / dir.y()
            };
            parmin = parmin.max(parcur);
            let parcur = if domain.is_open_ymax() {
                Precision::INFINITE
            } else {
                (ymax - loc.y()) / dir.y()
            };
            parmax = parmax.min(parcur);
            y_to_set = true;
        } else if dir.y() < 0.0 {
            let parcur = if domain.is_open_ymax() {
                -Precision::INFINITE
            } else {
                (ymax - loc.y()) / dir.y()
            };
            parmin = parmin.max(parcur);
            let parcur = if domain.is_open_ymin() {
                Precision::INFINITE
            } else {
                (ymin - loc.y()) / dir.y()
            };
            parmax = parmax.min(parcur);
            y_to_set = true;
        } else {
            // Perpendicular to axis Y
            if loc.y() < ymin || ymax < loc.y() {
                return;
            }
            ymin_set = loc.y();
            ymax_set = loc.y();
            y_to_set = false;
        }

        if dir.z() > 0.0 {
            let parcur = if domain.is_open_zmin() {
                -Precision::INFINITE
            } else {
                (zmin - loc.z()) / dir.z()
            };
            parmin = parmin.max(parcur);
            let parcur = if domain.is_open_zmax() {
                Precision::INFINITE
            } else {
                (zmax - loc.z()) / dir.z()
            };
            parmax = parmax.min(parcur);
            z_to_set = true;
        } else if dir.z() < 0.0 {
            let parcur = if domain.is_open_zmax() {
                -Precision::INFINITE
            } else {
                (zmax - loc.z()) / dir.z()
            };
            parmin = parmin.max(parcur);
            let parcur = if domain.is_open_zmin() {
                Precision::INFINITE
            } else {
                (zmin - loc.z()) / dir.z()
            };
            parmax = parmax.min(parcur);
            z_to_set = true;
        } else {
            // Perpendicular to axis Z
            if loc.z() < zmin || zmax < loc.z() {
                return;
            }
            zmin_set = loc.z();
            zmax_set = loc.z();
            z_to_set = false;
        }

        self.nb_seg += 1;
        self.begin_on_curve[0] = parmin;
        self.end_on_curve[0] = parmax;

        if x_to_set {
            let par1 = loc.x() + parmin * dir.x();
            let par2 = loc.x() + parmax * dir.x();
            xmin_set = par1.min(par2);
            xmax_set = par1.max(par2);
        }

        if y_to_set {
            let par1 = loc.y() + parmin * dir.y();
            let par2 = loc.y() + parmax * dir.y();
            ymin_set = par1.min(par2);
            ymax_set = par1.max(par2);
        }

        if z_to_set {
            let par1 = loc.z() + parmin * dir.z();
            let par2 = loc.z() + parmax * dir.z();
            zmin_set = par1.min(par2);
            zmax_set = par1.max(par2);
        }

        box_lin.update(xmin_set, ymin_set, zmin_set, xmax_set, ymax_set, zmax_set);
    }
}
