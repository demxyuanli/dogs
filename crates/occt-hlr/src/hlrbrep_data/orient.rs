//! `HLRBRep_Data::OrientOutLine` (`HLRBRep_Data.cxx:1746-1879`) and
//! `HLRBRep_Data::OrientOthEdge` (`HLRBRep_Data.cxx:1881-1991`).
//!
//! `Update` reaches both through the `FD.WithOutL() && !FD.Side()` arm
//! (`cxx:841-844`): once a face carries outline edges, every outline edge is
//! oriented against the face normal (`OrientOutLine`), and every remaining
//! edge of the face gets its orientation flipped when the projected normal
//! points away (`OrientOthEdge`).
//!
//! Neither function touches `HLRBRep_Intersector`: both are driven by
//! `HLRBRep_EdgeFaceTool::UVPoint` (`HLRBRep_EdgeFaceTool.cxx:61-104`), the
//! local surface properties `HLRBRep_SLProps` and, for `OrientOutLine`,
//! `HLRBRep_EdgeFaceTool::CurvatureValue` (`cxx:26-57`) plus the normal
//! curvature sign correction that follows it.
//!
//! `myFEGeom = &(ed1.ChangeGeometry())` and `const HLRBRep_Curve& EC =
//! ed1.Geometry()` both alias the same edge record in OCCT; the port reads the
//! immutable `geometry()` for both.

use occt_core::gp::{GpDir, GpPnt2d, GpVec};
use occt_core::precision::RESOLUTION;

use occt_topo::abs::Orientation;
use crate::hlrbrep_data::HlrBRepData;
use crate::hlrbrep_edge_face_tool::{curvature_value, uv_point};
use crate::hlrbrep_face_data::HlrBRepFaceData;

/// `gp_Dir(gp_Dir::D::NZ)`: the projection direction for a parallel projector.
fn normal_z() -> GpDir {
    GpDir::new(0.0, 0.0, 1.0).expect("gp_Dir(0, 0, 1)")
}

/// The 3D point at which `OrientOutLine` probes the face normal
/// (`HLRBRep_Data.cxx:1780-1795`): the start vertex parameter when only the end
/// vertex exists, the end vertex parameter when only the start exists, the
/// middle parameter when both exist, and `0` when neither does.
fn probe_parameter(ec: &crate::hlrbrep_curve::HlrBRepCurve, vsta: i32, vend: i32) -> f64 {
    if vsta == 0 && vend == 0 {
        0.0
    } else if vsta == 0 {
        ec.parameter3d(ec.last_parameter())
    } else if vend == 0 {
        ec.parameter3d(ec.first_parameter())
    } else {
        ec.parameter3d((ec.last_parameter() + ec.first_parameter()) / 2.0)
    }
}

impl HlrBRepData {
    /// `OrientOutLine(I, FD)` (`HLRBRep_Data.cxx:1746-1879`): orients every
    /// outline / internal edge of `fd` from the sign of the projected face
    /// normal crossed with the edge tangent, and reports whether the first
    /// closed non-cut face ended up inverted.
    pub(crate) fn orient_out_line(&mut self, fd: &mut HlrBRepFaceData) -> bool {
        let t = self.projector().transformation().clone();
        let ti = self.projector().inverted_transformation().clone();
        let perspective = self.projector().perspective();
        let focus = self.projector().focus().ok();

        let mut inverted = false;
        let mut first_inversion = true;

        let nw = fd.wires().map(|w| w.nb_wires()).unwrap_or(0);
        for iw1 in 1..=nw {
            let ne1 = match fd.wires().and_then(|w| w.wire(iw1)) {
                Some(b) => b.nb_edges(),
                None => continue,
            };
            for ie1 in 1..=ne1 {
                let (my_fe, is_double, is_iso, out_l, internal, old_ori) = {
                    let b = match fd.wires().and_then(|w| w.wire(iw1)) {
                        Some(b) => b,
                        None => continue,
                    };
                    (
                        b.edge(ie1),
                        b.is_double(ie1),
                        b.iso_line(ie1),
                        b.out_line(ie1),
                        b.internal(ie1),
                        b.orientation(ie1),
                    )
                };

                // `cxx:1765-1772`.
                let vertical = self.my_e_data[my_fe as usize].vertical();
                self.my_e_data[my_fe as usize].set_used(is_double || is_iso || vertical);

                if (out_l || internal) && !vertical {
                    let (p, uv) = {
                        let ed = &self.my_e_data[my_fe as usize];
                        let ec = ed.geometry();
                        let p = probe_parameter(ec, ed.v_sta(), ed.v_end());
                        (p, uv_point(p, ec, fd.geometry()))
                    };
                    if let Some((pu, pv)) = uv {
                        // `mySLProps.SetParameters(pu, pv); EC.D1(p, Pt, Tg);`
                        // (`cxx:1802-1803`).
                        self.my_sl_props_mut().set_parameters(pu, pv);
                        let (mut pt, mut tg) =
                            self.my_e_data[my_fe as usize].geometry().d1_3d(p);

                        // `gp_Dir V` (`cxx:1804-1815`).
                        let mut v = if perspective {
                            let mut p2d = GpPnt2d::default();
                            self.projector().project(&pt, &mut p2d);
                            GpDir::new(p2d.x(), p2d.y(), -focus.unwrap_or(0.0))
                                .unwrap_or_else(|_| normal_z())
                        } else {
                            normal_z()
                        };
                        v.transform(&ti);

                        if self.my_sl_props_mut().is_normal_defined() {
                            // `CurvatureValue(iFaceGeom, pu, pv, V)` and the
                            // `curv > 0` normal flip (`cxx:1818-1825`).
                            let curv = curvature_value(fd.geometry(), pu, pv, &v);
                            // `gp_Vec Nm = mySLProps.Normal()` (`cxx:1819`);
                            // `Normal()` only throws when `IsNormalDefined()`
                            // is false, which the branch above excluded.
                            let nrm = self
                                .my_sl_props()
                                .normal()
                                .copied()
                                .expect("mySLProps.Normal() after IsNormalDefined()");
                            let mut nm = GpVec::new(nrm.x(), nrm.y(), nrm.z());
                            if curv > 0.0 {
                                nm.reverse();
                            }
                            tg.transform(&t);
                            pt.transform(&t);
                            nm.transform(&t);
                            // `Nm.Cross(Tg)` (`cxx:1829`): in-place cross
                            // product, so `Nm` becomes `Nm ^ Tg`.
                            nm = nm.crossed(&tg);

                            // `cxx:1834-1841`.
                            let r = if perspective {
                                let f = focus.unwrap_or(0.0);
                                nm.z() * f
                                    - (nm.x() * pt.x() + nm.y() * pt.y() + nm.z() * pt.z())
                            } else {
                                nm.z()
                            };
                            let my_fe_ori = if r > 0.0 {
                                Orientation::Forward
                            } else {
                                Orientation::Reversed
                            };

                            // `cxx:1843-1851`.
                            if !fd.cut() && fd.closed() && first_inversion {
                                if (old_ori == my_fe_ori)
                                    != (fd.orientation() == Orientation::Forward)
                                {
                                    first_inversion = false;
                                    inverted = true;
                                }
                            }

                            if let Some(w) = fd.wires_mut() {
                                if let Some(b) = w.wire_mut(iw1) {
                                    b.set_orientation(ie1, my_fe_ori);
                                }
                            }
                        }
                    }
                    // `ed1.Used(true)` (`cxx:1877`).
                    self.my_e_data[my_fe as usize].set_used(true);
                }
            }
        }
        inverted
    }

    /// `OrientOthEdge(I, FD)` (`HLRBRep_Data.cxx:1881-1991`): for every edge of
    /// `fd` that `OrientOutLine` did not touch, flip its orientation when the
    /// projected face normal points away from the viewer.
    pub(crate) fn orient_oth_edge(&mut self, fd: &mut HlrBRepFaceData) {
        let t = self.projector().transformation().clone();
        let perspective = self.projector().perspective();
        let focus = self.projector().focus().ok();

        let nw = fd.wires().map(|w| w.nb_wires()).unwrap_or(0);
        for iw1 in 1..=nw {
            let ne1 = match fd.wires().and_then(|w| w.wire(iw1)) {
                Some(b) => b.nb_edges(),
                None => continue,
            };
            for ie1 in 1..=ne1 {
                let (my_fe, mut my_fe_ori) = {
                    let b = match fd.wires().and_then(|w| w.wire(iw1)) {
                        Some(b) => b,
                        None => continue,
                    };
                    (b.edge(ie1), b.orientation(ie1))
                };

                if self.my_e_data[my_fe as usize].used() {
                    continue;
                }
                // `ed1.Used(true);` (`cxx:1904`).
                self.my_e_data[my_fe as usize].set_used(true);

                let (_p, pt, uv) = {
                    let ec = self.my_e_data[my_fe as usize].geometry();
                    let p = ec.parameter3d((ec.last_parameter() + ec.first_parameter()) / 2.0);
                    (p, ec.value3d(p), uv_point(p, ec, fd.geometry()))
                };
                let Some((pu, pv)) = uv else { continue };

                // `mySLProps.SetParameters(pu, pv)` (`cxx:1909`).
                self.my_sl_props_mut().set_parameters(pu, pv);
                if !self.my_sl_props_mut().is_normal_defined() {
                    continue;
                }
                // `gp_Vec Nm = mySLProps.Normal()` (`cxx:1912`).
                let nrm = self
                    .my_sl_props()
                    .normal()
                    .copied()
                    .expect("mySLProps.Normal() after IsNormalDefined()");
                let mut nm = GpVec::new(nrm.x(), nrm.y(), nrm.z());
                let mut pt = pt;
                pt.transform(&t);
                nm.transform(&t);

                // `cxx:1917-1922`.
                let r = if perspective {
                    let f = focus.unwrap_or(0.0);
                    nm.z() * f - (nm.x() * pt.x() + nm.y() * pt.y() + nm.z() * pt.z())
                } else {
                    nm.z()
                };
                if r < 0.0 {
                    my_fe_ori = my_fe_ori.reversed();
                    if let Some(w) = fd.wires_mut() {
                        if let Some(b) = w.wire_mut(iw1) {
                            b.set_orientation(ie1, my_fe_ori);
                        }
                    }
                }
            }
        }
    }
}

/// `gp::Resolution()` is only read by `OrientOutLine`'s debug print of
/// `Tg.Magnitude()` (`HLRBRep_Data.cxx:1830-1836`); kept here so the constant
/// the `.cxx` names stays referenced.
#[allow(dead_code)]
const TANGENT_RESOLUTION: f64 = RESOLUTION;
