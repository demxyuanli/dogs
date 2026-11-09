//! `LocOpe_FindEdgesInFace` (`LocOpe_FindEdgesInFace.hxx:34-61`,
//! `.cxx:40-193`, `.lxx:23-52`).
//!
//! `Set(S, F)` keeps the edges of `S` that lie on the face `F`: an edge of `F`
//! itself is kept unconditionally (`cxx:87-98`), and every other edge is kept
//! only when its curve is a line or a circle **on** the face surface
//! (`cxx:102-158`) and three points of it classify as not-`OUT` inside the
//! face's UV domain (`cxx:161-190`).

use std::sync::Arc;

use occt_core::elib::slib;
use occt_core::gp::{GpCirc, GpLin, GpPln, GpPnt2d};
use occt_core::precision::{CONFUSION, ANGULAR};
use occt_geom::{Curve, Surface};

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face, TopoShape};
use crate::topexp::Explorer;
use crate::topo_tools::is_same;

pub struct FindEdgesInFace {
    my_shape: TopoShape,
    my_face: Face,
    my_list: Vec<Edge>,
    index: usize,
}

impl Default for FindEdgesInFace {
    fn default() -> Self {
        Self::new()
    }
}

impl FindEdgesInFace {
    /// `LocOpe_FindEdgesInFace()` (`.lxx:23`).
    pub fn new() -> Self {
        Self {
            my_shape: TopoShape::new(ShapeType::Compound),
            my_face: Face::wrap(TopoShape::new(ShapeType::Face))
                .expect("a FACE-typed TopoShape wraps as a Face"),
            my_list: Vec::new(),
            index: 0,
        }
    }

    /// `LocOpe_FindEdgesInFace(S, F)` (`.lxx:27-30`).
    pub fn with_shape_and_face(s: &TopoShape, f: &Face) -> Self {
        let mut r = Self::new();
        r.set(s, f);
        r
    }

    /// `Set(S, F)` (`.cxx:42-193`).
    pub fn set(&mut self, s: &TopoShape, f: &Face) {
        self.my_shape = s.clone();
        self.my_face = f.clone();
        self.my_list.clear();

        // `cxx:60-72`: only a plane or a cylinder is handled ("pour le moment");
        // the surface of a `Geom_RectangularTrimmedSurface` is unwrapped once.
        let Some(surf) = BRepTool::face_surface(f) else {
            return;
        };
        let surf = unwrap_trimmed(&surf);
        let plane = surf.gp_pln();
        let cylinder = surf.gp_cylinder();
        if plane.is_none() && cylinder.is_none() {
            return;
        }

        let Some(adaptor) = BRepTool::face_surface(f) else {
            return;
        };
        let mut topol_tool = TopolTool::new();
        if topol_tool.initialize(adaptor, f).is_err() {
            return;
        }

        let mut seen: Vec<TopoShape> = Vec::new();
        let mut ex = Explorer::new(s, ShapeType::Edge);
        while ex.more() {
            let Some(edg) = Edge::wrap(ex.current().clone()) else {
                ex.next();
                continue;
            };
            ex.next();
            // `cxx:81-86`: `if (M.Contains(edg)) continue;`
            if seen.iter().any(|k| is_same(k, &edg.0)) {
                continue;
            }
            // `cxx:87-99`: an edge shared with the face is taken as is.
            let mut shared = false;
            let mut exf = Explorer::new(f.shape(), ShapeType::Edge);
            while exf.more() {
                if is_same(exf.current(), &edg.0) {
                    shared = true;
                    break;
                }
                exf.next();
            }
            if shared {
                seen.push(edg.0.clone());
                self.my_list.push(edg);
                continue;
            }

            let Some(curve) = BRepTool::edge_curve_world(&edg) else {
                continue;
            };
            let (fpar, lpar) = BRepTool::edge_parameters(&edg);
            let curve = unwrap_trimmed_curve(&curve);
            // `cxx:108-112`: only lines and circles are handled.
            let line = curve.gp_line();
            let circle = curve.gp_circ();
            if line.is_none() && circle.is_none() {
                continue;
            }

            let to_add = match (&line, &circle) {
                (Some(li), _) => line_on_surface(li, &plane, &cylinder),
                (None, Some(ci)) => circle_on_surface(ci, &plane, &cylinder),
                (None, None) => false,
            };
            if !to_add {
                continue;
            }

            // `cxx:161-190`: classify three points (start, end, middle).
            let p = [
                curve.d0(fpar),
                curve.d0(lpar),
                curve.d0((fpar + lpar) / 2.0),
            ];
            let mut inside = true;
            for pt in &p {
                let (u, v) = match (&plane, &cylinder) {
                    (Some(pl), _) => slib::plane_parameters(&pl.position(), pt),
                    (None, Some(cy)) => slib::cylinder_parameters(&cy.position(), pt),
                    (None, None) => (0.0, 0.0),
                };
                if topol_tool.classify(GpPnt2d::new(u, v), CONFUSION, true) == FaceState::Out {
                    inside = false;
                    break;
                }
            }
            if inside {
                seen.push(edg.0.clone());
                self.my_list.push(edg);
            }
        }
    }

    /// `Init()` (`.lxx:34-37`).
    pub fn init(&mut self) {
        self.index = 0;
    }

    /// `More()` (`.lxx:41-44`).
    pub fn more(&self) -> bool {
        self.index < self.my_list.len()
    }

    /// `Edge()` (`.lxx:48-51`).
    pub fn edge(&self) -> Option<&Edge> {
        self.my_list.get(self.index)
    }

    /// `Next()` (`.lxx:55-58`).
    pub fn next(&mut self) {
        self.index += 1;
    }
}

/// `cxx:63-68`: `Geom_RectangularTrimmedSurface -> BasisSurface()`.
fn unwrap_trimmed(s: &Arc<dyn Surface>) -> Arc<dyn Surface> {
    match s.rectangular_trimmed_basis() {
        Some(b) => b,
        None => s.clone(),
    }
}

/// `cxx:103-108`: `Geom_TrimmedCurve -> BasisCurve()`.
fn unwrap_trimmed_curve(c: &Arc<dyn Curve>) -> Arc<dyn Curve> {
    match c.trim_basis() {
        Some(b) => b,
        None => c.clone(),
    }
}

/// `cxx:124-139`.
fn line_on_surface(li: &GpLin, plane: &Option<GpPln>, cylinder: &Option<occt_core::gp::GpCylinder>) -> bool {
    if let Some(pl) = plane {
        pl.contains_lin(li, CONFUSION, ANGULAR)
    } else if let Some(cy) = cylinder {
        cy.axis().is_parallel_tol(li.position(), ANGULAR)
            && (li.distance(&cy.location()) - cy.radius()).abs() < CONFUSION
    } else {
        false
    }
}

/// `cxx:141-158`.
fn circle_on_surface(
    ci: &GpCirc,
    plane: &Option<GpPln>,
    cylinder: &Option<occt_core::gp::GpCylinder>,
) -> bool {
    if let Some(pl) = plane {
        // `pl.Position().IsCoplanar(ci.Position(), Tol, TolAng)`: `gp_Ax3`'s
        // overload takes a `gp_Ax3`, and OCCT converts the `gp_Ax2` implicitly
        // (`gp_Ax3(const gp_Ax2&)`).
        pl.position()
            .is_coplanar(&occt_core::gp::GpAx3::from_ax2(&ci.position()), CONFUSION, ANGULAR)
    } else if let Some(cy) = cylinder {
        (cy.radius() - ci.radius()).abs() < CONFUSION
            && cy.axis().is_coaxial_tol(&ci.axis(), ANGULAR, CONFUSION)
    } else {
        false
    }
}
