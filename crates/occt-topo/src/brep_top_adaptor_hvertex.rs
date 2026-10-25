//! `BRepTopAdaptor_HVertex` — an `Adaptor3d_HVertex` over a `TopoDS_Vertex`
//! and the p-curve adaptor it lies on.
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepTopAdaptor/BRepTopAdaptor_HVertex.hxx:31-64`,
//! `.cxx:17-200`.
//!
//! Translation notes:
//! - The `occ::handle<BRepAdaptor_Curve2d> myCurve` member (`hxx:62`) is only
//!   read by `Parameter` / `Resolution` through their `Adaptor2d_Curve2d`
//!   argument, never as a fallback, so the port stores no field for it.
//! - `Value()` returns `gp_Pnt2d(RealFirst(), RealFirst())` by design
//!   (`cxx:36-40`, the commented-out original call) — "do nothing". `RealFirst()`
//!   is `f64::MIN` here (`Standard_Real.hxx`), not a geometric point.
//! - `Parameter(C)` / `Resolution(C)` down-cast their `Adaptor2d_Curve2d`
//!   argument to `BRepAdaptor_Curve2d` (`cxx:42-46`, `:50-53`); the port
//!   receives that adaptor directly, so the down-cast disappears.
//! - `IsSame(Other)` down-casts to `BRepTopAdaptor_HVertex` (`cxx:196-200`);
//!   likewise a direct `&HVertex`.

use occt_core::gp::{GpPnt2d, GpVec2d};
use occt_core::precision::Precision;

use crate::abs::Orientation;
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_tool::BRepTool;
use crate::shape::Vertex;
use crate::topo_tools::is_same;

/// `class BRepTopAdaptor_HVertex : public Adaptor3d_HVertex` (`hxx:31-64`).
#[derive(Debug, Clone)]
pub struct HVertex {
    /// `myVtx` (`hxx:61`).
    vtx: Vertex,
}

impl HVertex {
    /// `BRepTopAdaptor_HVertex(Vtx, Curve)` (`cxx:21-26`).
    pub fn new(v: &Vertex) -> Self {
        Self { vtx: v.clone() }
    }

    /// `Vertex()` (`hxx:43`).
    pub fn vertex(&self) -> &Vertex {
        &self.vtx
    }

    /// `ChangeVertex()` (`hxx:45`).
    pub fn change_vertex(&mut self) -> &mut Vertex {
        &mut self.vtx
    }

    /// `Value()` (`cxx:31-35`): the "do nothing" arm, `(RealFirst(), RealFirst())`.
    pub fn value(&self) -> GpPnt2d {
        GpPnt2d::new(f64::MIN, f64::MIN)
    }

    /// `Parameter(C)` (`cxx:42-46`): `BRep_Tool::Parameter(myVtx, C->Edge(),
    /// C->Face())`.
    pub fn parameter(&self, c: &Curve2dAdaptor) -> f64 {
        BRepTool::parameter_of_vertex_on_face(&self.vtx.0, c.edge(), c.face())
    }

    /// `Resolution(C)` (`cxx:48-190`).
    ///
    /// Verbatim control flow: the UV resolution is the `max` of the two curve
    /// resolutions at `UResolution(tv)` / `VResolution(tv)`; the parametric
    /// estimate `tv / |DC|` is capped at `1.0e7` when `tv > 1.0e7 * |DC|`; the
    /// estimate is walked one tolerance along the curve in the vertex's
    /// orientation, clamped to the curve range, and a distance check refines it
    /// twice. The result is `min(ResUV, ResUV1)`.
    pub fn resolution(&self, c: &Curve2dAdaptor) -> f64 {
        let Some(surf) = BRepTool::face_surface(c.face()) else {
            return 0.0;
        };
        let tv = BRepTool::vertex_tolerance(&self.vtx);
        let p = self.parameter(c);
        let or = self.orientation();
        let Some((p2d0, v2d0)) = c.d1(p) else {
            return 0.0;
        };
        let (p0, du, dv) = surf.d1(p2d0.x(), p2d0.y());
        let dc = GpVec2d::new(
            v2d0.x() * du.x() + v2d0.y() * dv.x(),
            v2d0.x() * du.y() + v2d0.y() * dv.y(),
        );
        let mut mag = dc.magnitude();

        let u_res = occt_geom::approx_same_parameter::u_resolution(surf.as_ref(), tv);
        let v_res = occt_geom::approx_same_parameter::v_resolution(surf.as_ref(), tv);
        let t_u_res = c.resolution(u_res).unwrap_or(0.0);
        let t_v_res = c.resolution(v_res).unwrap_or(0.0);
        let res_uv1 = t_u_res.max(t_v_res);

        if mag < 1e-12 {
            return res_uv1;
        }

        let mut res_uv;
        // "for lack of better options limit the parametric solution to
        // 10 million*tolerance of the point" (`cxx:105-108`).
        if tv > 1.0e7 * mag {
            res_uv = 1.0e7;
        } else {
            res_uv = tv / mag;
        }

        let mut pp;
        if or == Orientation::Reversed {
            pp = p + res_uv;
        } else {
            pp = p - res_uv;
        }

        let u_min = c.first_parameter();
        let u_max = c.last_parameter();
        if pp > u_max {
            pp = u_max;
        }
        if pp < u_min {
            pp = u_min;
        }

        let Some(p2d) = c.d0(pp) else {
            return res_uv1;
        };
        let p1 = surf.d0(p2d.x(), p2d.y());
        let mut dist = p0.distance(&p1);
        if dist > 1e-12 && (dist > 1.1 * tv || dist < 0.8 * tv) {
            // Refine if possible (`cxx:135-190`).
            let mut pp2;
            if or == Orientation::Reversed {
                pp2 = p + tv / dist;
            } else {
                pp2 = p - tv / dist;
            }
            if pp2 > u_max {
                pp2 = u_max;
            }
            if pp2 < u_min {
                pp2 = u_min;
            }

            let (p2d_b, v2d_b) = match c.d1(pp2) {
                Some(v) => v,
                None => return res_uv1,
            };
            let (p1_b, du_b, dv_b) = surf.d1(p2d_b.x(), p2d_b.y());
            let dc_b = GpVec2d::new(
                v2d_b.x() * du_b.x() + v2d_b.y() * dv_b.x(),
                v2d_b.x() * du_b.y() + v2d_b.y() * dv_b.y(),
            );
            let dist1_b = p0.distance(&p1_b);
            let mut dist1 = dist1_b;
            if (dist1 - tv).abs() < (dist - tv).abs() {
                // Take the result of interpolation (`cxx:160-164`).
                res_uv = tv / dist;
                dist = dist1;
            }

            mag = dc_b.magnitude();
            if tv > 1.0e7 * mag {
                mag = tv * 1.0e-7;
            }
            let mut pp3;
            if or == Orientation::Reversed {
                pp3 = p + tv / mag;
            } else {
                pp3 = p - tv / mag;
            }
            if pp3 > u_max {
                pp3 = u_max;
            }
            if pp3 < u_min {
                pp3 = u_min;
            }

            let Some(p2d_c) = c.d0(pp3) else {
                return res_uv1;
            };
            let p1_c = surf.d0(p2d_c.x(), p2d_c.y());
            dist1 = p0.distance(&p1_c);
            if (dist1 - tv).abs() < (dist - tv).abs() {
                // Take the new estimation (`cxx:184-188`).
                res_uv = tv / mag;
                dist = dist1;
            }
        }
        let _ = dist;

        res_uv.min(res_uv1)
    }

    /// `Orientation()` (`cxx:194-197`): `myVtx.Orientation()`.
    pub fn orientation(&self) -> Orientation {
        self.vtx.0.orientation()
    }

    /// `IsSame(Other)` (`cxx:198-200`): `myVtx.IsSame(brhv->Vertex())`.
    pub fn is_same(&self, other: &HVertex) -> bool {
        is_same(&self.vtx.0, &other.vtx.0)
    }
}

/// `Adaptor3d_HVertex::Resolution` reads `BRep_Tool::Tolerance(V)` raised to
/// `Precision::Confusion()`; the constant is kept in scope so the port's
/// `vertex_tolerance` mapping stays visible.
#[allow(dead_code)]
const CONFUSION: f64 = Precision::CONFUSION;
