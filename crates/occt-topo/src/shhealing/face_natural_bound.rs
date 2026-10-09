//! `ShapeFix_Face::FixAddNaturalBound` (`ShapeFix_Face.cxx:876-1108`).
//!
//! Detects a missing natural boundary on a double-closed surface (sphere,
//! torus, closed BSpline) and adds it, first shifting every existing wire into
//! the free interval left by the others. `ShapeFix_Face::Perform` runs it
//! (`cxx:704-708`) after the wire fixes and before `FixSplitFace` (`cxx:711`),
//! and its return value is what turns `NeedSplit` off.
//!
//! UNPORTED (`cxx:1030-1078`): the sphere arm that merges a hole whose
//! degenerated edge coincides with a degenerated edge of the natural bound. It
//! needs `ShapeExtend_WireData::SetLast` / `Add(sewd, at)` and
//! `ShapeFix_Wire::Load` / `FixConnected` / `FixDegenerated` on a wire data,
//! none of which exist as such in the port. The arm only fires when
//! `GetType() == GeomAbs_Sphere && ws.Length() == nb + 1`, so it is registered
//! rather than approximated.

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::CONFUSION;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::Face;

use super::face_geom_helpers::{
    cut_interval, empty_copied_face, find_best_interval, is_surface_uv_infinite, shift_2d_wire,
    shift_is_negligible, split_face_children,
};
use super::ShapeFixFace;
use crate::shape_fix_compose_shell::ReShape;

impl ShapeFixFace {
    /// `ShapeFix_Face::FixAddNaturalBound()` (`ShapeFix_Face.cxx:876-1108`).
    /// Returns `isAdded`; on success `myFace` is the rebuilt face and the
    /// context records `Replace(old_face, new_face)` (`cxx:917`, `cxx:1098`).
    pub fn fix_add_natural_bound(&mut self) -> bool {
        // cxx:878-881.
        let Some(surf) = self.surf.clone() else {
            return false;
        };
        // cxx:883-887: `myFace = TopoDS::Face(Context()->Apply(myFace))`.
        let mut face = match self.face.clone() {
            Some(f) => f,
            None => return false,
        };
        {
            let applied = self.context.apply(&face.0);
            if applied.is_face() {
                face = Face(applied);
            }
        }
        // cxx:889-904.
        let (mut ws, vs) = split_face_children(&face);
        let builder = TopoBuilder::new();

        // cxx:906-939: an empty face gets a standard natural-bound face.
        if ws.is_empty() && !is_surface_uv_infinite(surf.as_ref()) {
            let mut new_face =
                crate::brep_lib_make_face::make_face_from_surface(surf.clone(), CONFUSION);
            new_face.0.set_orientation(face.0.orientation()); // cxx:913
            self.context.replace(&face.0, &new_face.0); // cxx:917
            // cxx:923-930: `ShapeFix_Edge::FixVertexTolerance(edg, myFace)` on
            // every edge of the new face.
            for e in crate::topo_tools_full::edges_of(&new_face.0) {
                crate::shhealing::fix_vertex_tolerance_edge(&e);
            }
            // cxx:932-937: `BRepTools::Update` is a no-op in the port;
            // `myResult = myFace`.
            self.face = Some(new_face.clone());
            self.result = Some(new_face.0.clone());
            return true;
        }

        // cxx:941-945.
        if !self.is_need_add_natural_bound(&ws) {
            return false;
        }

        // cxx:947-955.
        let (suf, sul) = surf.u_range();
        let (svf, svl) = surf.v_range();
        let mut int_u: Vec<GpPnt2d> = vec![GpPnt2d::new(suf, sul)];
        let mut int_v: Vec<GpPnt2d> = vec![GpPnt2d::new(svf, svl)];
        let nb = ws.len();
        let mut centers: Vec<GpPnt2d> = Vec::with_capacity(nb);
        // `ShapeAnalysis_Surface::IsUClosed()` / `IsVClosed()`
        // (`ShapeAnalysis_Surface.cxx:661` / `:868`): their first test is
        // `mySurf->IsUClosed()` / `IsVClosed()`, which every surface reaching
        // here (`IsSurfaceUVPeriodic` held) answers immediately, so the port
        // reads the surface closure directly.
        let is_u_closed = surf.is_u_closed();
        let is_v_closed = surf.is_v_closed();

        for i in 0..nb {
            // cxx:957-965: `aWireFace = TopoDS::Face(myFace.EmptyCopied())`,
            // `aB.Add(aWireFace, aw)`, `ShapeAnalysis::GetFaceUVBounds`.
            let a_wire_face = builder.make_face(surf.clone(), std::slice::from_ref(&ws[i]));
            let (umin, umax, vmin, vmax) = crate::brep_surface::face_uv_bounds(&a_wire_face);
            if is_u_closed {
                cut_interval(&mut int_u, &GpPnt2d::new(umin, umax), sul - suf);
            }
            if is_v_closed {
                cut_interval(&mut int_v, &GpPnt2d::new(vmin, vmax), svl - svf);
            }
            centers.push(GpPnt2d::new(0.5 * (umin + umax), 0.5 * (vmin + vmax)));
        }

        // cxx:979-987.
        let mut shift = GpPnt2d::new(0.0, 0.0);
        if is_u_closed {
            shift = GpPnt2d::new(find_best_interval(&int_u), shift.y());
        }
        if is_v_closed {
            shift = GpPnt2d::new(shift.x(), find_best_interval(&int_v));
        }

        // cxx:989-1003: move every existing wire into the free interval.
        let center = GpPnt2d::new(
            shift.x() + 0.5 * (sul - suf),
            shift.y() + 0.5 * (svl - svf),
        );
        for (i, wire) in ws.iter().enumerate() {
            let mut sh = GpPnt2d::new(0.0, 0.0);
            if is_u_closed {
                sh = GpPnt2d::new(
                    crate::shhealing::adjust_by_period(centers[i].x(), center.x(), sul - suf),
                    sh.y(),
                );
            }
            if is_v_closed {
                sh = GpPnt2d::new(
                    sh.x(),
                    crate::shhealing::adjust_by_period(centers[i].y(), center.y(), svl - svf),
                );
            }
            // cxx:1002: `Shift2dWire(wire, myFace, sh.XY(), mySurf)`.
            shift_2d_wire(
                wire,
                &face,
                &GpVec2d::new(sh.x(), sh.y()),
                surf.as_ref(),
                false,
            );
        }

        // cxx:1005-1027: the natural bound itself, plus the shift it needs.
        let location = face.0.location().clone();
        let bound_surf = BRepTool::face_surface(&face).unwrap_or_else(|| surf.clone());
        let mut ftmp = crate::brep_lib_make_face::make_face_from_surface(bound_surf, CONFUSION);
        ftmp.0.set_location(&location); // cxx:1011
        let (bound_wires, _) = split_face_children(&ftmp);
        for wire in bound_wires {
            ws.push(wire.clone());
            if shift_is_negligible(&shift) {
                continue; // cxx:1020-1023
            }
            shift_2d_wire(
                &wire,
                &face,
                &GpVec2d::new(shift.x(), shift.y()),
                surf.as_ref(),
                true, // cxx:1024
            );
        }

        // UNPORTED (`cxx:1028-1078`): the `GeomAbs_Sphere && ws.Length() == nb
        // + 1` merge of a hole that shares a degenerated edge with the natural
        // bound. See the module header.

        // cxx:1080-1102: rebuild the face over `ws` then `vs`.
        let mut s = empty_copied_face(&face);
        s.0.set_orientation(Orientation::Forward); // cxx:1083
        for w in &ws {
            builder.add(&mut s.0, &w.0); // cxx:1085-1088
        }
        for v in &vs {
            builder.add(&mut s.0, v); // cxx:1089-1092
        }
        if !self.my_fwd {
            s.0.set_orientation(Orientation::Reversed); // cxx:1093-1096
        }
        self.context.replace(&face.0, &s.0); // cxx:1098
        // cxx:1100-1101: `myFace = TopoDS::Face(S)` + `BRepTools::Update`.
        self.face = Some(s);
        true
    }
}
