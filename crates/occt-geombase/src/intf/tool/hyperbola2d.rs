//! `Intf_Tool::Hypr2dBox` (`Intf_Tool.cxx:210-369`) and
//! `Intf_Tool::Inters2d(gp_Hypr2d, Bnd_Box2d)` (`Intf_Tool.cxx:373-473`).

use occt_foundation::precision::Precision;
use occt_gp::bnd::BndBox2d;
use occt_gp::elib::clib2d::{hyperbola_d1_ax22d, hyperbola_value_ax22d};
use occt_gp::gp::{GpHypr2d, GpLin2d, GpPnt2d, GpXY};
use occt_gp::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};

use super::helpers::{dir2d_nx, dir2d_ny, dir2d_x, dir2d_y};
use super::structure::IntfTool;

impl IntfTool {
    /// `Intf_Tool::Hypr2dBox` (`Intf_Tool.cxx:210-369`).
    pub fn hypr2d_box(
        &mut self,
        hypr2d: &GpHypr2d,
        domain: &BndBox2d,
        box_hypr2d: &mut BndBox2d,
    ) {
        self.nb_seg = 0;
        box_hypr2d.set_void();
        if domain.is_whole() {
            box_hypr2d.set_whole();
            self.nb_seg = 1;
            self.begin_on_curve[0] = -Precision::INFINITE;
            self.end_on_curve[0] = Precision::INFINITE;
            return;
        } else if domain.is_void() {
            return;
        }

        let nb_pi = self.inters2d_hypr(hypr2d, domain);

        if nb_pi > 0 {
            let Some((mut xmin, mut ymin, mut xmax, mut ymax)) = domain.get() else {
                return;
            };

            for npi in 0..nb_pi as usize {
                xmin = xmin.min(self.xint[npi]);
                xmax = xmax.max(self.xint[npi]);
                ymin = ymin.min(self.yint[npi]);
                ymax = ymax.max(self.yint[npi]);
            }
            box_hypr2d.update(xmin, ymin, xmax, ymax);

            // Selection sort of (parint, bord) by parint (`Intf_Tool.cxx:245-266`).
            for npi in 0..nb_pi as usize {
                let mut npk = npi;
                for npj in (npi + 1)..nb_pi as usize {
                    if self.parint[npj] < self.parint[npk] {
                        npk = npj;
                    }
                }
                if npk != npi {
                    self.parint.swap(npk, npi);
                    self.bord.swap(npk, npi);
                }
            }

            let mut sinan = 0.0f64;
            let mut out = true;

            for npi in 0..nb_pi as usize {
                let (_, tan) = hyperbola_d1_ax22d(
                    self.parint[npi],
                    &hypr2d.pos,
                    hypr2d.major_radius,
                    hypr2d.minor_radius,
                );
                match self.bord[npi] {
                    1 => sinan = GpXY::new(-1.0, 0.0).crossed(&tan.xy()),
                    2 => sinan = GpXY::new(0.0, -1.0).crossed(&tan.xy()),
                    3 => sinan = GpXY::new(1.0, 0.0).crossed(&tan.xy()),
                    4 => sinan = GpXY::new(0.0, 1.0).crossed(&tan.xy()),
                    _ => {}
                }
                if sinan.abs() > Precision::ANGULAR {
                    if sinan > 0.0 {
                        if self.nb_seg < 6 {
                            out = false;
                            self.begin_on_curve[self.nb_seg] = self.parint[npi];
                            self.nb_seg += 1;
                        }
                    } else {
                        if out && self.nb_seg < 6 {
                            self.begin_on_curve[self.nb_seg] = -Precision::INFINITE;
                            self.nb_seg += 1;
                        }
                        if self.nb_seg > 0 {
                            self.end_on_curve[self.nb_seg - 1] = self.parint[npi];
                        }
                        out = true;

                        let ipmin = if self.begin_on_curve[self.nb_seg - 1] < -10.0 {
                            -10
                        } else {
                            self.begin_on_curve[self.nb_seg - 1] as i32
                        };
                        let ipmax = if self.end_on_curve[self.nb_seg - 1] > 10.0 {
                            10
                        } else {
                            self.end_on_curve[self.nb_seg - 1] as i32
                        };
                        let ipmin = ipmin * 10 + 1;
                        let ipmax = ipmax * 10 - 1;
                        let mut ip = ipmin;
                        while ip <= ipmax {
                            let pas = if ip.abs() <= 10 { 1 } else { 10 };
                            box_hypr2d.add_point(&hyperbola_value_ax22d(
                                f64::from(ip) / 10.0,
                                &hypr2d.pos,
                                hypr2d.major_radius,
                                hypr2d.minor_radius,
                            ));
                            ip += pas;
                        }
                    }
                }
            }
            if !out && self.nb_seg > 0 {
                self.end_on_curve[self.nb_seg - 1] = Precision::INFINITE;
            }
        } else if !domain.is_out(&hyperbola_value_ax22d(
            0.0,
            &hypr2d.pos,
            hypr2d.major_radius,
            hypr2d.minor_radius,
        )) {
            *box_hypr2d = *domain;
            self.begin_on_curve[0] = -Precision::INFINITE;
            self.end_on_curve[0] = Precision::INFINITE;
            self.nb_seg = 1;
        }
    }

    /// `Intf_Tool::Inters2d(gp_Hypr2d, Bnd_Box2d)` (`Intf_Tool.cxx:373-473`).
    fn inters2d_hypr(&mut self, curv: &GpHypr2d, domain: &BndBox2d) -> i32 {
        let mut nbpi: i32 = 0;

        let Some((xmin, ymin, xmax, ymax)) = domain.get() else {
            return 0;
        };

        if !domain.is_open_ymax() {
            let l1 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, ymax), dir2d_nx());
            let mut inters1 = IntAna2dAnaIntersection::new();
            inters1.perform_hypr_conic(curv, &IntAna2dConic::from_lin2d(&l1));
            if inters1.is_done() && !inters1.is_empty() {
                for npi in 1..=inters1.nb_points() {
                    let x = inters1.point(npi).value().x();
                    self.xint[nbpi as usize] = x;
                    if xmin < x && x <= xmax {
                        self.yint[nbpi as usize] = ymax;
                        self.parint[nbpi as usize] = inters1.point(npi).param_on_first();
                        self.bord[nbpi as usize] = 1;
                        nbpi += 1;
                    }
                }
            }
        }

        if !domain.is_open_xmin() {
            let l2 = GpLin2d::from_pnt_dir(GpPnt2d::new(xmin, 0.0), dir2d_ny());
            let mut inters2 = IntAna2dAnaIntersection::new();
            inters2.perform_hypr_conic(curv, &IntAna2dConic::from_lin2d(&l2));
            if inters2.is_done() && !inters2.is_empty() {
                for npi in 1..=inters2.nb_points() {
                    let y = inters2.point(npi).value().y();
                    self.yint[npi as usize] = y;
                    if ymin < y && y <= ymax {
                        self.xint[npi as usize] = xmin;
                        self.parint[npi as usize] = inters2.point(npi).param_on_first();
                        self.bord[npi as usize] = 2;
                        nbpi += 1;
                    }
                }
            }
        }

        if !domain.is_open_ymin() {
            let l3 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, ymin), dir2d_x());
            let mut inters3 = IntAna2dAnaIntersection::new();
            inters3.perform_hypr_conic(curv, &IntAna2dConic::from_lin2d(&l3));
            if inters3.is_done() && !inters3.is_empty() {
                for npi in 1..=inters3.nb_points() {
                    let x = inters3.point(npi).value().x();
                    self.xint[npi as usize] = x;
                    if xmin <= x && x < xmax {
                        self.yint[npi as usize] = ymin;
                        self.parint[npi as usize] = inters3.point(npi).param_on_first();
                        self.bord[npi as usize] = 3;
                        nbpi += 1;
                    }
                }
            }
        }

        if !domain.is_open_xmax() {
            let l4 = GpLin2d::from_pnt_dir(GpPnt2d::new(xmax, 0.0), dir2d_y());
            let mut inters4 = IntAna2dAnaIntersection::new();
            inters4.perform_hypr_conic(curv, &IntAna2dConic::from_lin2d(&l4));
            if inters4.is_done() && !inters4.is_empty() {
                for npi in 1..=inters4.nb_points() {
                    let y = inters4.point(npi).value().y();
                    self.yint[npi as usize] = y;
                    if ymin <= y && y < ymax {
                        self.xint[npi as usize] = xmax;
                        self.parint[npi as usize] = inters4.point(npi).param_on_first();
                        self.bord[npi as usize] = 4;
                        nbpi += 1;
                    }
                }
            }
        }
        nbpi
    }
}
