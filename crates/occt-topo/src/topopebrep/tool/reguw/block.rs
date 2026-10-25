//! `TopOpeBRepTool_REGUW::InitBlock` / `NearestE` / `NextinBlock` and the
//! `FUN_nextdata` helper.
//!
//! Sources: `TopOpeBRepTool_REGUW.cxx:482-522` (`FUN_nextdata`),
//! `:526-674` (`InitBlock`), `:678-786` (`NearestE`), `:788-921`
//! (`NextinBlock`).

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::ANGULAR;

use crate::abs::{Orientation, ShapeType};
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key};

use super::super::closed::{
    is_closing_e_in_wire, par_e, vertex, CLOSING,
};
use super::super::geometry::{fun_tool_line_2d, fun_tool_quad_2d, try_tg_2d_app};
use super::super::matter::matter_2d;
use super::Reguw;

/// `FUN_nextdata(iStep, e, pc, v, p2d, tg)` (`cxx:482-522`).
///
/// Answers `(v, p2d, tg)`; the 2-D tangent is reapproximated on a quadric
/// (non-line) p-curve, because the end tangent of such a curve can be
/// degenerate (`cxx:502-518`).
pub(crate) fn fun_nextdata(
    i_step: i32,
    e: &TopoShape,
    pc: &dyn occt_geom2d::curve::Curve2d,
) -> Option<(TopoShape, GpPnt2d, GpVec2d)> {
    let iv1e1 = if i_step == 1 {
        super::super::closed::FORWARD
    } else {
        super::super::closed::REVERSED
    };
    let v = vertex(iv1e1, e)?;
    let par1 = par_e(iv1e1, e)?;

    let line = fun_tool_line_2d(pc);
    let quad = fun_tool_quad_2d(pc);
    let approx = quad && !line;

    let (p2d, tg2d) = if approx {
        let _p2d = pc.d0(par1);
        let iv0e1 = if i_step == 1 {
            super::super::closed::REVERSED
        } else {
            super::super::closed::FORWARD
        };
        let par0 = par_e(iv0e1, e)?;
        // the point near the other end; only the tangent is kept (`cxx:508-517`).
        let x = 0.2345;
        let par = (1.0 - x) * par1 + x * par0;
        let (_pmil, tg2d) = pc.d1(par);
        (_p2d, tg2d)
    } else {
        pc.d1(par1)
    };

    let mut tg = tg2d.normalized().unwrap_or(tg2d);
    if e.orientation() == Orientation::Reversed {
        tg.reverse();
    }
    Some((v, p2d, tg))
}

impl Reguw {
    /// `InitBlock()` (`cxx:526-674`): pick the starting vertex and edge of the
    /// current wire, and the starting 2-D point and tangent.
    pub fn init_block(&mut self) -> bool {
        assert!(self.has_init(), "TopOpeBRepTool_REGUW : NO INIT");

        self.my_v0 = None;
        self.my_p2d0 = GpPnt2d::new(1.0e7, 1.0e7);
        self.my_ed = None;
        let iv0e1 = self.iv0e1();

        // updating <mymapvmultiple> and myListVmultiple (`cxx:545-560`)
        let mut keep: Vec<TopoShape> = Vec::new();
        for vmu in self.my_list_v_multiple.clone() {
            let key = shape_map_key(&vmu);
            let mult = self
                .my_map_v_eds
                .get(&key)
                .map(|c| c.is_multiple())
                .unwrap_or(false);
            if mult {
                keep.push(vmu);
            } else {
                self.my_map_v_multiple.remove(&key);
            }
        }
        self.my_list_v_multiple = keep;

        // myv0 (`cxx:562-587`)
        if self.my_list_v_multiple.is_empty() {
            for key in self.my_v_order.clone() {
                let Some(co) = self.my_map_v_eds.get(&key) else {
                    continue;
                };
                let mut lea: Vec<TopoShape> = Vec::new();
                let nea = co.item(iv0e1, &mut lea);
                let mut leb: Vec<TopoShape> = Vec::new();
                let neb = co.item(CLOSING, &mut leb);
                let ne = nea + neb;
                if ne != 0 {
                    self.my_v0 = co.key().cloned();
                    break;
                }
            }
        } else {
            self.my_v0 = self.my_list_v_multiple.first().cloned();
        }
        let Some(v0) = self.my_v0.clone() else {
            return false;
        };

        // myed (`cxx:589-627`)
        let Some(co) = self.my_map_v_eds.get(&shape_map_key(&v0)) else {
            return false;
        };
        let mut lea: Vec<TopoShape> = Vec::new();
        let nea = co.item(iv0e1, &mut lea);
        let mut leb: Vec<TopoShape> = Vec::new();
        let neb = co.item(CLOSING, &mut leb);
        let mut le: Vec<TopoShape> = Vec::new();
        le.extend(lea.iter().cloned());
        le.extend(leb.iter().cloned());
        let ne = nea + neb;
        if ne == 0 {
            return false;
        }

        if nea > 0 {
            self.my_ed = le.first().cloned();
        } else {
            // `<myv0> CLOSING in <myed>` (`cxx:608-626`). `cxx:613` reads
            // `IsClosingE(myed, ...)` while `myed` is still null, so `iscE` is
            // false and the first element is taken.
            let s = self.my_s.clone();
            let f_ref = self.my_corr_iso.f_ref().clone();
            for eb in &le {
                let isc_e = match (&self.my_ed, &s) {
                    (Some(ed), Some(s)) => is_closing_e_in_wire(ed, s, &f_ref),
                    _ => false,
                };
                if !isc_e {
                    self.my_ed = Some(eb.clone());
                    break;
                }
                let iov0 = super::super::closed::ori_in_sor_closed(&v0, eb);
                if iov0 != iv0e1 {
                    continue;
                }
            }
        }
        let Some(myed) = self.my_ed.clone() else {
            return false;
        };

        // `cxx:629-638`: `myed` leaves the connexity of its vertices.
        for vcur in explorer_shapes(&myed, ShapeType::Vertex) {
            if let Some(cco) = self.my_map_v_eds.get_mut(&shape_map_key(&vcur)) {
                cco.remove_item(&myed);
            }
        }

        // myp2d0 (`cxx:640-650`)
        let Some(c2df) = self.my_corr_iso.uv_rep(&myed) else {
            return false;
        };
        let Some(pc) = c2df.curve.clone() else {
            return false;
        };
        let Some(par0) = par_e(iv0e1, &myed) else {
            return false;
        };
        self.my_p2d0 = pc.d0(par0);

        // myv, myp2d, mytg2d (`cxx:652-653`)
        let Some((v, p2d, tg)) = fun_nextdata(self.i_step, &myed, pc.as_ref()) else {
            return false;
        };
        self.my_v = Some(v);
        self.my_p2d = p2d;
        self.my_tg2d = tg;
        true
    }

    /// `NearestE(loe, efound)` (`cxx:678-786`): the edge of `loe` whose tangent
    /// makes the smallest matter (step 1) / non-matter (step 2) angle with
    /// `mytg2d`. Equal angles are ambiguous and answer `false`.
    pub fn nearest_e(&self, loe: &[TopoShape], efound: &mut Option<TopoShape>) -> bool {
        *efound = None;
        let fac = 0.45678;
        let tola = ANGULAR;
        let iv0e1 = self.iv0e1();

        let mut it = loe.iter();
        let Some(first) = it.next() else {
            return true;
        };
        let mut found = first.clone();

        let Some(c2defound) = self.my_corr_iso.uv_rep(&found) else {
            return false;
        };
        let mut tg2dfound = try_tg_2d_app(iv0e1, &found, &c2defound, fac);
        if found.orientation() == Orientation::Reversed {
            tg2dfound.reverse();
        }
        let mut angfound = if self.i_step == 1 {
            matter_2d(&self.my_tg2d, &tg2dfound)
        } else {
            2.0 * std::f64::consts::PI - matter_2d(&tg2dfound, &self.my_tg2d)
        };

        for ei in it {
            // for an INTERNAL edge eI -> eF+eR
            if let Some(myed) = &self.my_ed {
                if is_same(ei, myed) {
                    continue;
                }
            }
            let Some(c2dei) = self.my_corr_iso.uv_rep(ei) else {
                return false;
            };
            let mut tg2di = try_tg_2d_app(iv0e1, ei, &c2dei, fac);
            if ei.orientation() == Orientation::Reversed {
                tg2di.reverse();
            }
            let angi = if self.i_step == 1 {
                matter_2d(&self.my_tg2d, &tg2di)
            } else {
                2.0 * std::f64::consts::PI - matter_2d(&tg2di, &self.my_tg2d)
            };
            let eq = (angi - angfound).abs() < tola;
            if eq {
                return false;
            }
            if angi < angfound {
                found = ei.clone();
                angfound = angi;
            }
        }
        *efound = Some(found);
        true
    }

    /// `NextinBlock()` (`cxx:788-921`): step to the next edge of the wire.
    pub fn next_in_block(&mut self) -> bool {
        assert!(self.has_init(), "TopOpeBRepTool_REGUW : NO INIT");

        let iv0e1 = self.iv0e1();
        let Some(myv) = self.my_v.clone() else {
            return false;
        };
        let Some(co) = self.my_map_v_eds.get(&shape_map_key(&myv)) else {
            return false;
        };

        // {e} : e is connexed to <myv> && ori(<myv>,e)=iv0e1 (`cxx:807-819`)
        let mut lea: Vec<TopoShape> = Vec::new();
        let nea = co.item(iv0e1, &mut lea);
        let mut leb: Vec<TopoShape> = Vec::new();
        let neb = co.item(CLOSING, &mut leb);
        let mut le: Vec<TopoShape> = Vec::new();
        le.extend(lea.iter().cloned());
        le.extend(leb.iter().cloned());
        let _ = nea + neb;

        // `cxx:821-873`: keep the edges whose bound at `<myv>` falls on
        // `<myp2d>`.
        let mut kept: Vec<TopoShape> = Vec::new();
        for e in le {
            if let Some(myed) = &self.my_ed {
                // xpu240299 FRA60275, fsp4
                if is_same(&e, myed) {
                    continue;
                }
            }
            let Some(c2df) = self.my_corr_iso.uv_rep(&e) else {
                return false;
            };
            let Some(pc) = c2df.curve.clone() else {
                return false;
            };
            let Some(par) = par_e(iv0e1, &e) else {
                return false;
            };
            let p2de = pc.d0(par);
            // p2d(myv,ed)=p2d(myv,e)
            if p2de.is_equal(&self.my_p2d, self.my_tol2d) {
                kept.push(e);
            }
        }
        let ne = kept.len();
        if ne == 0 {
            return false;
        }

        // myed (`cxx:875-891`)
        if ne == 1 {
            self.my_ed = Some(kept[0].clone());
        } else {
            let mut efound: Option<TopoShape> = None;
            if !self.nearest_e(&kept, &mut efound) {
                return false;
            }
            self.my_ed = efound;
        }
        let Some(myed) = self.my_ed.clone() else {
            return false;
        };

        // `cxx:893-902`: `myed` leaves the connexity of its vertices.
        for vcur in explorer_shapes(&myed, ShapeType::Vertex) {
            if let Some(cco) = self.my_map_v_eds.get_mut(&shape_map_key(&vcur)) {
                cco.remove_item(&myed);
            }
        }

        // `cxx:904-907`: `UVRep` is called without checking its answer; a
        // missing entry would dereference a null handle in OCCT.
        let Some(c2df) = self.my_corr_iso.uv_rep(&myed) else {
            return false;
        };
        let Some(pc) = c2df.curve.clone() else {
            return false;
        };

        // `cxx:908`: myv, myp2d, mytg2d
        let Some((v, p2d, tg)) = fun_nextdata(self.i_step, &myed, pc.as_ref()) else {
            return false;
        };
        self.my_v = Some(v);
        self.my_p2d = p2d;
        self.my_tg2d = tg;
        true
    }
}
