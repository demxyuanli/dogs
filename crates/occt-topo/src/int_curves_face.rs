//! `IntCurvesFace_Intersector` / `IntCurvesFace_ShapeIntersector`.
//!
//! Source: `IntCurvesFace_Intersector.cxx` (`InternalCall` at 228,
//! `Perform(Lin)` at ~360). A line is intersected with the face surface by
//! [`crate::intcurvesurface::perform_curve_surface`] (`IntCurveSurface_HInter`);
//! each CS point is classified in UV by [`crate::fclass2d::FClass2d`]. IN/ON
//! points are kept, and the CS transition is reversed when the face is
//! `REVERSED`.
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::bnd::BndBox;
use occt_core::gp::{GpLin, GpPnt, GpPnt2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::GeomLine;

use crate::abs::Orientation;
use crate::bbox_from_geometry::shape_bbox;
use crate::brep_surface::{classify_surface, face_uv_bounds, SurfaceKind};
use crate::brep_tool::BRepTool;
use crate::fclass2d::{FaceState, FClass2d};
use crate::intcurvesurface::{self, IntersectionPoint, State};
use crate::intcurvesurface::{
    ThePolygonOfHInter as ThePolygon, ThePolyhedronOfHInter as ThePolyhedron,
};
use crate::shape::{Face, TopoShape};
use crate::topo_tools_full::faces_of;

/// `IntCurveSurface_TransitionOnCurve`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Transition {
    In,
    Out,
    Tangent,
}

/// One line/face hit kept by the intersector.
#[derive(Debug, Clone, Copy)]
pub struct FaceHit {
    pub w: f64,
    pub u: f64,
    pub v: f64,
    pub pnt: GpPnt,
    pub transition: Transition,
    pub state: FaceState,
}

/// `IntCurvesFace_Intersector`.
pub struct FaceIntersector {
    face: Face,
    tol: f64,
    use_bound_tol: bool,
    done: bool,
    parallel: bool,
    box_: BndBox,
    /// Parameter window of the intersected surface (see [`adaptor_domain`]).
    domain: (f64, f64, f64, f64),
    polyhedron: Option<ThePolyhedron>,
    pub(crate) pnts: Vec<FaceHit>,
}

impl FaceIntersector {
    /// `IntCurvesFace_Intersector(Face, Tol, aRestr, UseBToler)`.
    pub fn new(face: Face, tol: f64, restr: bool, use_bound_tol: bool) -> Self {
        let box_ = shape_bbox(&face.0);
        let domain = adaptor_domain(&face, restr);
        let polyhedron = build_polyhedron(&face, domain);
        Self {
            face,
            tol: if tol > 0.0 { tol } else { CONFUSION },
            use_bound_tol,
            done: false,
            parallel: false,
            box_,
            domain,
            polyhedron,
            pnts: Vec::new(),
        }
    }

    pub fn bounding(&self) -> BndBox {
        self.box_
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn is_parallel(&self) -> bool {
        self.parallel
    }

    pub fn nb_pnt(&self) -> i32 {
        self.pnts.len() as i32
    }

    pub fn w_parameter(&self, i: i32) -> f64 {
        self.pnts[(i as usize).saturating_sub(1)].w
    }

    pub fn u_parameter(&self, i: i32) -> f64 {
        self.pnts[(i as usize).saturating_sub(1)].u
    }

    pub fn v_parameter(&self, i: i32) -> f64 {
        self.pnts[(i as usize).saturating_sub(1)].v
    }

    pub fn pnt(&self, i: i32) -> GpPnt {
        self.pnts[(i as usize).saturating_sub(1)].pnt
    }

    pub fn transition(&self, i: i32) -> Transition {
        self.pnts[(i as usize).saturating_sub(1)].transition
    }

    pub fn state(&self, i: i32) -> FaceState {
        self.pnts[(i as usize).saturating_sub(1)].state
    }

    pub fn classify_uv_point(&self, puv: GpPnt2d) -> FaceState {
        FClass2d::new(&self.face, self.tol)
            .map(|c| c.perform(puv))
            .unwrap_or(FaceState::Unknown)
    }

    /// `Perform(Lin, PInf, PSup)`.
    pub fn perform(&mut self, lin: &GpLin, pinf: f64, psup: f64) {
        self.done = false;
        self.parallel = false;
        self.pnts.clear();
        let (mut lo, mut hi) = if pinf <= psup {
            (pinf, psup)
        } else {
            (psup, pinf)
        };
        // `Intf_Tool::LinBox` (`IntCurvesFace_Intersector.cxx:392-447`): clip
        // the (possibly infinite) line to the face polyhedron box before the
        // polygon out-test. `ThePolygon::of_line` with +/-inf yields a void
        // box and would reject every hit.
        if let Some(ph) = &self.polyhedron {
            match lin_box_clip(lin, ph.bounding(), lo, hi) {
                None => {
                    self.done = true;
                    return;
                }
                Some((a, b)) => {
                    lo = a;
                    hi = b;
                }
            }
            let poly = ThePolygon::of_line(lin, lo, hi);
            if ph.is_out_polygon(&poly) {
                self.done = true;
                return;
            }
        }
        let Some(surf) = BRepTool::face_surface(&self.face) else {
            self.done = true;
            return;
        };
        let curve = GeomLine::new(*lin);
        let (u0, u1, v0, v1) = self.domain;
        let Ok(res) =
            intcurvesurface::perform_curve_surface(&curve, surf.as_ref(), (lo, hi), (u0, u1, v0, v1))
        else {
            self.done = true;
            return;
        };
        if res.nb_segments() > 0 {
            self.parallel = true;
        }
        let cl = match FClass2d::new(&self.face, self.tol) {
            Ok(c) => c,
            Err(_) => {
                for &p in res.points() {
                    self.pnts.push(FaceHit {
                        w: p.param,
                        u: p.u,
                        v: p.v,
                        pnt: p.pnt,
                        transition: match p.state {
                            State::In => Transition::In,
                            State::Out => Transition::Out,
                            State::On | State::Unknown => Transition::Tangent,
                        },
                        state: FaceState::In,
                    });
                }
                self.pnts.sort_by(|a, b| a.w.total_cmp(&b.w));
                self.done = true;
                return;
            }
        };
        for &p in res.points() {
            if p.param + PCONFUSION < lo || p.param - PCONFUSION > hi {
                continue;
            }
            if let Some(hit) = keep_point(&self.face, &cl, p) {
                self.pnts.push(hit);
            }
        }
        self.pnts.sort_by(|a, b| a.w.total_cmp(&b.w));
        self.done = true;
    }
}

fn keep_point(face: &Face, cl: &FClass2d, p: IntersectionPoint) -> Option<FaceHit> {
    // OCCT `IntCurvesFace_Intersector::InternalCall` classifies the intersection
    // point at its own `(U, V)` (`IntCurvesFace_Intersector.cxx:254`, `:293-294`)
    // and never re-projects it. `intcurvesurface` reports those parameters in the
    // surface's own frame (`IntCurveSurface_InterUtils.pxx:901-925`), and
    // `Perform` above feeds it the very `BRepTool::face_surface` the classifier
    // uses, so the re-projection this used to do is redundant.
    let (u, v) = (p.u, p.v);
    let st = cl.perform(GpPnt2d::new(u, v));
    if st != FaceState::In && st != FaceState::On {
        return None;
    }
    if st == FaceState::On {
        if let Some(surf) = BRepTool::face_surface(face) {
            if surf.d0(u, v).distance(&p.pnt) > BRepTool::face_tolerance(face).max(CONFUSION) * 10.0
            {
                return None;
            }
        }
    }
    let mut tran = match p.state {
        State::In => Transition::In,
        State::Out => Transition::Out,
        State::On | State::Unknown => Transition::Tangent,
    };
    if tran != Transition::Tangent && face.0.orientation() == Orientation::Reversed {
        tran = match tran {
            Transition::In => Transition::Out,
            Transition::Out => Transition::In,
            Transition::Tangent => Transition::Tangent,
        };
    }
    Some(FaceHit {
        w: p.param,
        u,
        v,
        pnt: p.pnt,
        transition: tran,
        state: st,
    })
}

fn lin_box_clip(lin: &GpLin, bbox: &BndBox, tmin: f64, tmax: f64) -> Option<(f64, f64)> {
    let (x0, x1, y0, y1, z0, z1) = bbox.get()?;
    let o = lin.location();
    let d = lin.direction();
    let mut t0 = tmin;
    let mut t1 = tmax;
    let slabs = [
        (x0, x1, o.x(), d.x()),
        (y0, y1, o.y(), d.y()),
        (z0, z1, o.z(), d.z()),
    ];
    for (p0, p1, orig, dir) in slabs {
        if dir.abs() <= 1e-30 {
            if orig < p0 || orig > p1 {
                return None;
            }
            continue;
        }
        let mut ta = (p0 - orig) / dir;
        let mut tb = (p1 - orig) / dir;
        if ta > tb {
            std::mem::swap(&mut ta, &mut tb);
        }
        t0 = t0.max(ta);
        t1 = t1.min(tb);
        if t0 > t1 {
            return None;
        }
    }
    let pad = 0.05 * (t1 - t0);
    t0 -= pad;
    t1 += pad;
    if (t1 - t0) < 1e-10 {
        t0 -= 1e-10;
        t1 += 1e-10;
    }
    if t0 > tmax || t1 < tmin {
        return None;
    }
    t0 = t0.max(tmin);
    t1 = t1.min(tmax);
    if t0 > t1 - 1e-9 {
        return None;
    }
    Some((t0, t1))
}

/// `BRepTools::UVBounds` of `face`, made finite for the intersection engines.
///
/// `IntCurvesFace_Intersector` / `IntTools_EdgeFace` both hand
/// `IntCurveSurface_HInter` a surface adaptor; the port's engines need a finite
/// UV window instead, and this is the one the port's `IntCurvesFace_Intersector`
/// uses (unbounded directions get the `[-1e3, 1e3]` proxy window, a degenerate
/// window is widened by 1).
pub(crate) fn finite_uv(face: &Face) -> (f64, f64, f64, f64) {
    let (u0, u1, v0, v1) = face_uv_bounds(face);
    finite_window(u0, u1, v0, v1)
}

/// Parameter window of the `BRepAdaptor_Surface` that
/// `IntCurvesFace_Intersector` builds with `Initialize(Face, aRestr)`
/// (`IntCurvesFace_Intersector.cxx:137`): the face UV box when `restr` is set,
/// otherwise the native range of the underlying surface.
fn adaptor_domain(face: &Face, restr: bool) -> (f64, f64, f64, f64) {
    if restr {
        return finite_uv(face);
    }
    match BRepTool::face_surface(face) {
        Some(surf) => {
            let (u0, u1) = surf.u_range();
            let (v0, v1) = surf.v_range();
            finite_window(u0, u1, v0, v1)
        }
        None => finite_uv(face),
    }
}

fn finite_window(u0: f64, u1: f64, v0: f64, v1: f64) -> (f64, f64, f64, f64) {
    let clamp = |a: f64, b: f64| {
        let mut x = a;
        let mut y = b;
        if !x.is_finite() {
            x = -1.0e3;
        }
        if !y.is_finite() {
            y = 1.0e3;
        }
        if y < x {
            std::mem::swap(&mut x, &mut y);
        }
        if (y - x).abs() < PCONFUSION {
            y = x + 1.0;
        }
        (x, y)
    };
    let (u0, u1) = clamp(u0, u1);
    let (v0, v1) = clamp(v0, v1);
    (u0, u1, v0, v1)
}

fn build_polyhedron(face: &Face, domain: (f64, f64, f64, f64)) -> Option<ThePolyhedron> {
    let surf = BRepTool::face_surface(face)?;
    if matches!(
        classify_surface(surf.as_ref()),
        SurfaceKind::Plane
            | SurfaceKind::Cylinder
            | SurfaceKind::Cone
            | SurfaceKind::Sphere
            | SurfaceKind::Torus
    ) {
        return None;
    }
    let (u0, u1, v0, v1) = domain;
    Some(ThePolyhedron::new(
        surf.as_ref(),
        20,
        20,
        u0,
        v0,
        u1,
        v1,
    ))
}

/// `IntCurvesFace_ShapeIntersector` — one line against every face of `shape`.
pub struct ShapeIntersector {
    faces: Vec<FaceIntersector>,
    done: bool,
    hits: Vec<(usize, FaceHit)>,
}

impl ShapeIntersector {
    pub fn load(shape: &TopoShape, tol: f64) -> Self {
        let faces = faces_of(shape)
            .into_iter()
            .map(|f| FaceIntersector::new(f, tol, true, true))
            .collect();
        Self {
            faces,
            done: false,
            hits: Vec::new(),
        }
    }

    pub fn perform(&mut self, lin: &GpLin, pinf: f64, psup: f64) {
        self.hits.clear();
        for (i, fi) in self.faces.iter_mut().enumerate() {
            fi.perform(lin, pinf, psup);
            if fi.is_done() {
                for k in 1..=fi.nb_pnt() {
                    self.hits.push((i, fi.pnts[(k as usize) - 1]));
                }
            }
        }
        self.hits.sort_by(|a, b| a.1.w.total_cmp(&b.1.w));
        self.done = true;
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn nb_pnt(&self) -> i32 {
        self.hits.len() as i32
    }

    pub fn w_parameter(&self, i: i32) -> f64 {
        self.hits[(i as usize).saturating_sub(1)].1.w
    }

    pub fn pnt(&self, i: i32) -> GpPnt {
        self.hits[(i as usize).saturating_sub(1)].1.pnt
    }

    pub fn face_index(&self, i: i32) -> usize {
        self.hits[(i as usize).saturating_sub(1)].0
    }
}
