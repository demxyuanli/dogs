//! ShapeFix_ComposeShell::SplitByGrid (ShapeFix_ComposeShell.cxx:2131-2275):
//! split every wire segment by the U- and V-seams of the composite surface.

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpDir2d, GpLin2d, GpPnt2d, GpVec2d};

use crate::brep_tools::add_uv_bounds_on_wire;
use crate::builder::TopoBuilder;

use super::helpers::*;
use super::shell::ComposeShell;
use super::wire_segment::WireSegment;

impl ComposeShell {
    /// ShapeFix_ComposeShell::SplitByGrid (cxx:2131-2275).
    pub fn split_by_grid(&mut self, seqw: &mut Vec<WireSegment>) {
        // cxx:2135-2138.
        let face = match self.face() {
            Some(f) => f.clone(),
            None => return,
        };
        let (uf, ul, vf, vl) = crate::brep_tools::uv_bounds(&face);
        let (umin, umax, vmin, vmax) = self.grid.bounds();

        // cxx:2142.
        let pprec = TOLINT;

        if self.closed_mode {
            // cxx:2152-2199. OCCT builds `myFace.EmptyCopied()` with the
            // segment's wire and asks `ShapeAnalysis::GetFaceUVBounds`. The
            // pcurves live on the shared surface, so `BRepTools::AddUVBounds
            // (Face, Wire)` on the real face is the same box.
            for w in seqw.iter_mut() {
                let builder = TopoBuilder::new();
                let wire = builder.make_wire(w.edges());
                let mut box2d = BndBox2d::new();
                add_uv_bounds_on_wire(&face, &wire, &mut box2d);
                let (uf1b, vf1b, ul1b, vl1b) = match box2d.get() {
                    Some(v) => v,
                    None => (0.0, 0.0, 0.0, 0.0),
                };
                let mut uf1 = uf1b;
                let mut ul1 = ul1b;
                let mut vf1 = vf1b;
                let mut vl1 = vl1b;
                // cxx:2165-2174.
                let shift_u = if self.closed_mode && self.u_closed {
                    adjust_to_period(ul1 - pprec, self.grid.u_joint_value(1), self.grid.u_joint_value(2))
                } else {
                    0.0
                };
                let shift_v = if self.closed_mode && self.v_closed {
                    adjust_to_period(vl1 - pprec, self.grid.v_joint_value(1), self.grid.v_joint_value(2))
                } else {
                    0.0
                };
                uf1 += shift_u;
                ul1 += shift_u;
                vf1 += shift_v;
                vl1 += shift_v;
                // cxx:2182-2198.
                let iumin = 0.max(get_patch_index(uf1 + pprec, self.grid.u_joint_values(), self.u_closed));
                let iumax = get_patch_index(ul1 - pprec, self.grid.u_joint_values(), self.u_closed) + 1;
                let ivmin = 0.max(get_patch_index(vf1 + pprec, self.grid.v_joint_values(), self.v_closed));
                let ivmax = get_patch_index(vl1 - pprec, self.grid.v_joint_values(), self.v_closed) + 1;
                for j in 1..=w.nb_edges() {
                    w.define_iu_min(j, iumin);
                    w.define_iu_max(j, iumax);
                    w.define_iv_min(j, ivmin);
                    w.define_iv_max(j, ivmax);
                }
            }
        } else {
            // cxx:2203-2225.
            let iumin = get_patch_index(uf + pprec, self.grid.u_joint_values(), self.u_closed);
            let iumax = get_patch_index(ul - pprec, self.grid.u_joint_values(), self.u_closed) + 1;
            for w in seqw.iter_mut() {
                for j in 1..=w.nb_edges() {
                    w.define_iu_min(j, iumin);
                    w.define_iu_max(j, iumax);
                }
            }
            let ivmin = get_patch_index(vf + pprec, self.grid.v_joint_values(), self.v_closed);
            let ivmax = get_patch_index(vl - pprec, self.grid.v_joint_values(), self.v_closed) + 1;
            for w in seqw.iter_mut() {
                for j in 1..=w.nb_edges() {
                    w.define_iv_min(j, ivmin);
                    w.define_iv_max(j, ivmax);
                }
            }
        }

        // cxx:2227-2251: split by U lines.
        let u_start = if self.u_closed { 1 } else { 2 };
        for i in u_start..=self.grid.nb_u_patches() {
            let pos = GpPnt2d::new(self.grid.u_joint_value(i), 0.0);
            let dir = GpDir2d::new(0.0, 1.0).expect("direction");
            let line = GpLin2d::from_pnt_dir(pos, dir);
            if !self.closed_mode && self.u_closed {
                let period = umax - umin;
                let x = pos.x();
                let mut sh = adjust_to_period(x, uf, uf + period);
                while x + sh <= ul + pprec {
                    let ln = line.translated_vec(&GpVec2d::new(sh, 0.0));
                    let cut_index = get_patch_index(x + sh + pprec, self.grid.u_joint_values(), self.u_closed);
                    self.split_by_line_wires(seqw, &ln, true, cut_index);
                    sh += period;
                }
            } else {
                self.split_by_line_wires(seqw, &line, true, i as i32);
            }
        }

        // cxx:2253-2273: split by V lines.
        let v_start = if self.v_closed { 1 } else { 2 };
        for i in v_start..=self.grid.nb_v_patches() {
            let pos = GpPnt2d::new(0.0, self.grid.v_joint_value(i));
            let dir = GpDir2d::new(1.0, 0.0).expect("direction");
            let line = GpLin2d::from_pnt_dir(pos, dir);
            if !self.closed_mode && self.v_closed {
                let period = vmax - vmin;
                let y = pos.y();
                let mut sh = adjust_to_period(y, vf, vf + period);
                while y + sh <= vl + pprec {
                    let ln = line.translated_vec(&GpVec2d::new(0.0, sh));
                    let cut_index = get_patch_index(y + sh + pprec, self.grid.v_joint_values(), self.v_closed);
                    self.split_by_line_wires(seqw, &ln, false, cut_index);
                    sh += period;
                }
            } else {
                self.split_by_line_wires(seqw, &line, false, i as i32);
            }
        }
    }
}
