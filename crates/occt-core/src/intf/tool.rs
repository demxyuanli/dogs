//! Port of `Intf_Tool` (`src/ModelingAlgorithms/TKGeomAlgo/Intf/Intf_Tool.hxx` +
//! `Intf_Tool.cxx`).
//!
//! `Intf_Tool` builds the part of a bounding box cut by an *infinite* line,
//! hyperbola or parabola: for each "context" (a `Bnd_Box` / `Bnd_Box2d` domain)
//! it returns the clipped box plus the parameter ranges of the curve segments
//! that run through that domain (`NbSegments`, `BeginParam`, `EndParam`).
//!
//! Ported members (OCCT 8.0.0, tag V8_0_0):
//! * `Intf_Tool()`                        `Intf_Tool.cxx:38-48`
//! * `Lin2dBox`                           `Intf_Tool.cxx:52-206`
//! * `Hypr2dBox`                          `Intf_Tool.cxx:210-369`
//! * `Inters2d(gp_Hypr2d, Bnd_Box2d)`     `Intf_Tool.cxx:373-473`
//! * `Parab2dBox`                         `Intf_Tool.cxx:477-638`
//! * `Inters2d(gp_Parab2d, Bnd_Box2d)`    `Intf_Tool.cxx:642-742`
//! * `LinBox`                             `Intf_Tool.cxx:746-963`
//! * `NbSegments`                         `Intf_Tool.cxx:1634-1637`
//! * `BeginParam`                         `Intf_Tool.cxx:1641-1645`
//! * `EndParam`                           `Intf_Tool.cxx:1649-1653`
//!
//! UNPORTED members:
//! * `HyprBox`              `Intf_Tool.cxx:967-1119`
//! * `Inters3d(gp_Hypr)`    `Intf_Tool.cxx:1123-1302`
//! * `Inters3d(gp_Parab)`   `Intf_Tool.cxx:1306-1485`
//! * `ParabBox`             `Intf_Tool.cxx:1489-1630`
//!
//! Reason: all four call `IntAna_IntConicQuad` (`Intf_Tool.cxx:1133, 1156,
//! 1179, 1202, 1225, 1248, 1316, 1339, 1362, 1385, 1408, 1431`) together with
//! the 3D `ElCLib::D1(gp_Hypr/gp_Parab, ...)` (`Intf_Tool.cxx:1015,
//! 1534`). Neither `IntAna_IntConicQuad` nor the 3D conic
//! `ElCLib::D1`/`Value` overloads exist in this crate yet (only
//! `IntAna2d` and the 2D `elib::clib2d` overloads are ported), so no
//! faithful body can be written. They are deliberately absent rather than
//! approximated.
//!
//! Note: the surrounding `intf::mod` prose mentions `LinTetra` /
//! `LinSphere`, but neither symbol exists in this class - the real member
//! list is the ten entries above (`Intf_Tool.hxx:42-66`).

use crate::bnd::{BndBox, BndBox2d};
use crate::elib::clib2d::{
    hyperbola_d1_ax22d, hyperbola_value_ax22d, parabola_d1_ax22d, parabola_value_ax22d,
};
use crate::gp::{GpDir, GpDir2d, GpHypr2d, GpLin, GpLin2d, GpParab2d, GpPnt2d, GpXY};
use crate::intana2d::{IntAna2dAnaIntersection, IntAna2dConic};
use crate::precision::Precision;

/// `gp_Dir2d::D::NX`, the direction `(-1, 0)` used by the `Inters2d`
/// boundary lines (`Intf_Tool.cxx:383, 652`).
fn dir2d_nx() -> GpDir2d {
    GpDir2d::new(-1.0, 0.0).expect("gp_Dir2d::D::NX")
}

/// `gp_Dir2d::D::NY`, the direction `(0, -1)` (`Intf_Tool.cxx:406, 675`).
fn dir2d_ny() -> GpDir2d {
    GpDir2d::new(0.0, -1.0).expect("gp_Dir2d::D::NY")
}

/// `gp_Dir2d::D::X`, the direction `(1, 0)` (`Intf_Tool.cxx:429, 698`).
fn dir2d_x() -> GpDir2d {
    GpDir2d::new(1.0, 0.0).expect("gp_Dir2d::D::X")
}

/// `gp_Dir2d::D::Y`, the direction `(0, 1)` (`Intf_Tool.cxx:452, 721`).
fn dir2d_y() -> GpDir2d {
    GpDir2d::new(0.0, 1.0).expect("gp_Dir2d::D::Y")
}

/// `Bnd_Box2d::Set(const gp_Pnt2d&, const gp_Dir2d&)` (`Bnd_Box2d.hxx:97-102`):
/// `SetVoid()` + `Gap = 0` + `Add(P, D)`.
fn set2d_pnt_dir(b: &mut BndBox2d, p: &GpPnt2d, d: &GpDir2d) {
    b.set_void();
    b.add_point(p);
    b.add_dir(d);
}

/// `Bnd_Box::Add(const gp_Dir&)` (`Bnd_Box.cxx:627-655`). `BndBox`
/// in this crate has no equivalent method, so the port spells out the `Open*`
/// calls; the `RealEpsilon()` comparisons are `f64::EPSILON`.
fn add3d_dir(b: &mut BndBox, d: &GpDir) {
    let eps = f64::EPSILON;
    if d.x() < -eps {
        b.open_xmin();
    } else if d.x() > eps {
        b.open_xmax();
    }
    if d.y() < -eps {
        b.open_ymin();
    } else if d.y() > eps {
        b.open_ymax();
    }
    if d.z() < -eps {
        b.open_zmin();
    } else if d.z() > eps {
        b.open_zmax();
    }
}

/// `Intf_Tool` (`Intf_Tool.hxx:35-85`).
#[derive(Clone, Debug)]
pub struct IntfTool {
    /// `nbSeg` (`Intf_Tool.hxx:77`).
    nb_seg: usize,
    /// `beginOnCurve[6]` (`Intf_Tool.hxx:78`).
    begin_on_curve: [f64; 6],
    /// `endOnCurve[6]` (`Intf_Tool.hxx:79`).
    end_on_curve: [f64; 6],
    /// `bord[12]` (`Intf_Tool.hxx:80`).
    bord: [i32; 12],
    /// `xint[12]` (`Intf_Tool.hxx:81`).
    xint: [f64; 12],
    /// `yint[12]` (`Intf_Tool.hxx:82`).
    yint: [f64; 12],
    /// `zint[12]` (`Intf_Tool.hxx:83`). Written by the 3D conic paths
    /// only, which are UNPORTED (see the module header).
    #[allow(dead_code)]
    zint: [f64; 12],
    /// `parint[12]` (`Intf_Tool.hxx:84`).
    parint: [f64; 12],
}

impl Default for IntfTool {
    fn default() -> Self {
        Self::new()
    }
}

impl IntfTool {
    /// `Intf_Tool::Intf_Tool()` (`Intf_Tool.cxx:38-48`).
    pub fn new() -> Self {
        Self {
            nb_seg: 0,
            begin_on_curve: [0.0; 6],
            end_on_curve: [0.0; 6],
            bord: [0; 12],
            xint: [0.0; 12],
            yint: [0.0; 12],
            zint: [0.0; 12],
            parint: [0.0; 12],
        }
    }

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

    /// `Intf_Tool::Hypr2dBox` (`Intf_Tool.cxx:210-369`).
    pub fn hypr2d_box(
        &mut self,
        the_hypr2d: &GpHypr2d,
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

        let nb_pi = self.inters2d_hypr(the_hypr2d, domain);

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
                    &the_hypr2d.pos,
                    the_hypr2d.major_radius,
                    the_hypr2d.minor_radius,
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
                                &the_hypr2d.pos,
                                the_hypr2d.major_radius,
                                the_hypr2d.minor_radius,
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
            &the_hypr2d.pos,
            the_hypr2d.major_radius,
            the_hypr2d.minor_radius,
        )) {
            *box_hypr2d = *domain;
            self.begin_on_curve[0] = -Precision::INFINITE;
            self.end_on_curve[0] = Precision::INFINITE;
            self.nb_seg = 1;
        }
    }

    /// `Intf_Tool::Inters2d(gp_Hypr2d, Bnd_Box2d)` (`Intf_Tool.cxx:373-473`).
    fn inters2d_hypr(&mut self, the_curv: &GpHypr2d, domain: &BndBox2d) -> i32 {
        let mut nbpi: i32 = 0;

        let Some((xmin, ymin, xmax, ymax)) = domain.get() else {
            return 0;
        };

        if !domain.is_open_ymax() {
            let l1 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, ymax), dir2d_nx());
            let mut inters1 = IntAna2dAnaIntersection::new();
            inters1.perform_hypr_conic(the_curv, &IntAna2dConic::from_lin2d(&l1));
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
            inters2.perform_hypr_conic(the_curv, &IntAna2dConic::from_lin2d(&l2));
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
            inters3.perform_hypr_conic(the_curv, &IntAna2dConic::from_lin2d(&l3));
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
            inters4.perform_hypr_conic(the_curv, &IntAna2dConic::from_lin2d(&l4));
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

    /// `Intf_Tool::Parab2dBox` (`Intf_Tool.cxx:477-638`).
    pub fn parab2d_box(
        &mut self,
        the_parab2d: &GpParab2d,
        domain: &BndBox2d,
        box_parab2d: &mut BndBox2d,
    ) {
        self.nb_seg = 0;
        box_parab2d.set_void();
        if domain.is_whole() {
            box_parab2d.set_whole();
            self.nb_seg = 1;
            self.begin_on_curve[0] = -Precision::INFINITE;
            self.end_on_curve[0] = Precision::INFINITE;
            return;
        } else if domain.is_void() {
            return;
        }

        let nb_pi = self.inters2d_parab(the_parab2d, domain);

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
            box_parab2d.update(xmin, ymin, xmax, ymax);

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
                let (_, tan) =
                    parabola_d1_ax22d(self.parint[npi], &the_parab2d.pos, the_parab2d.focal);
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
                            box_parab2d.add_point(&parabola_value_ax22d(
                                f64::from(ip) / 10.0,
                                &the_parab2d.pos,
                                the_parab2d.focal,
                            ));
                            ip += pas;
                        }
                    }
                }
            }
            if !out && self.nb_seg > 0 {
                self.end_on_curve[self.nb_seg - 1] = Precision::INFINITE;
            }
        } else if !domain.is_out(&parabola_value_ax22d(0.0, &the_parab2d.pos, the_parab2d.focal)) {
            *box_parab2d = *domain;
            self.begin_on_curve[0] = -Precision::INFINITE;
            self.end_on_curve[0] = Precision::INFINITE;
            self.nb_seg = 1;
        }
    }

    /// `Intf_Tool::Inters2d(gp_Parab2d, Bnd_Box2d)` (`Intf_Tool.cxx:642-742`).
    fn inters2d_parab(&mut self, the_curv: &GpParab2d, domain: &BndBox2d) -> i32 {
        let mut nbpi: i32 = 0;

        let Some((xmin, ymin, xmax, ymax)) = domain.get() else {
            return 0;
        };

        if !domain.is_open_ymax() {
            let l1 = GpLin2d::from_pnt_dir(GpPnt2d::new(0.0, ymax), dir2d_nx());
            let mut inters1 = IntAna2dAnaIntersection::new();
            inters1.perform_parab_conic(the_curv, &IntAna2dConic::from_lin2d(&l1));
            if inters1.is_done() && !inters1.is_empty() {
                for npi in 1..=inters1.nb_points() {
                    let x = inters1.point(npi).value().x();
                    self.xint[npi as usize] = x;
                    if xmin < x && x <= xmax {
                        self.yint[npi as usize] = ymax;
                        self.parint[npi as usize] = inters1.point(npi).param_on_first();
                        self.bord[npi as usize] = 1;
                        nbpi += 1;
                    }
                }
            }
        }

        if !domain.is_open_xmin() {
            let l2 = GpLin2d::from_pnt_dir(GpPnt2d::new(xmin, 0.0), dir2d_ny());
            let mut inters2 = IntAna2dAnaIntersection::new();
            inters2.perform_parab_conic(the_curv, &IntAna2dConic::from_lin2d(&l2));
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
            inters3.perform_parab_conic(the_curv, &IntAna2dConic::from_lin2d(&l3));
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
            inters4.perform_parab_conic(the_curv, &IntAna2dConic::from_lin2d(&l4));
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

    /// `Intf_Tool::NbSegments()` (`Intf_Tool.cxx:1634-1637`).
    pub fn nb_segments(&self) -> usize {
        self.nb_seg
    }

    /// `Intf_Tool::BeginParam(SegmentNum)` (`Intf_Tool.cxx:1641-1645`).
    ///
    /// The argument is 1-based; `Standard_OutOfRange` (`Intf_Tool.cxx:1643`)
    /// is raised as a panic.
    pub fn begin_param(&self, segment_num: usize) -> f64 {
        assert!(
            segment_num >= 1 && segment_num <= self.nb_seg,
            "Intf_Tool::BeginParam"
        );
        self.begin_on_curve[segment_num - 1]
    }

    /// `Intf_Tool::EndParam(SegmentNum)` (`Intf_Tool.cxx:1649-1653`).
    ///
    /// The argument is 1-based; `Standard_OutOfRange` (`Intf_Tool.cxx:1651`)
    /// is raised as a panic.
    pub fn end_param(&self, segment_num: usize) -> f64 {
        assert!(
            segment_num >= 1 && segment_num <= self.nb_seg,
            "Intf_Tool::EndParam"
        );
        self.end_on_curve[segment_num - 1]
    }
}
