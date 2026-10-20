//! `ParametricDelaun` (face-basis Delaun exposing parametric UV to deflection
//! control) and the `CellsCount` / `ComputeErrFactors` helpers.

use std::collections::HashMap;

use occt_core::gp::{GpPnt, GpPnt2d, GpXY};
use occt_core::precision::SQUARE_CONFUSION;
use occt_geom::geom_api::project_point_on_surface;
use occt_geom::Surface;

use crate::abs::Orientation;

use super::super::data_model::*;
use super::super::delaun::Delaun;
use super::super::delaun_data::DelaunDataStructure;
use super::super::delaun_types::{DelaunLink, DelaunTriangle, DelaunVertex, VertexState};
use super::super::face_discret::Classifier;
#[allow(unused_imports)]
use super::super::mesh_algo::*;
use super::super::face_discret::PointState;
use super::super::deflection_control::{
    DelaunayDeflectionControlMeshAlgo, DeflectionMesh,
};
use super::super::parameters::MeshParameters;
use super::super::range_splitter::{
    mesh_surface_type, create_range_splitter, factory_uses_deflection_control,
    grab_params_of_internal_edges, RangeSplitter, MeshSurfaceType,
};
use super::*;

/// Delaun in face-basis coords, exposing parametric UV to deflection control
/// (`getNodePoint2d` = `Scale(coord, false)`).
pub(super) struct ParametricDelaun<'a> {
    pub(super) delaun: &'a mut Delaun,
    pub(super) ru: (f64, f64),
    pub(super) rv: (f64, f64),
    pub(super) delta: (f64, f64),
}

impl ParametricDelaun<'_> {
    fn to_param(&self, p: GpPnt2d) -> GpPnt2d {
        GpPnt2d::new(
            p.x() * self.delta.0 + self.ru.0,
            p.y() * self.delta.1 + self.rv.0,
        )
    }

    fn to_face(&self, p: GpPnt2d) -> GpPnt2d {
        if self.delta.0.abs() <= 1e-16 || self.delta.1.abs() <= 1e-16 {
            return p;
        }
        GpPnt2d::new(
            (p.x() - self.ru.0) / self.delta.0,
            (p.y() - self.rv.0) / self.delta.1,
        )
    }

    /// `DelaunayNodeInsertionMeshAlgo::insertNodes` (`hxx:104-131`).
    pub(super) fn insert_parametric(
        &mut self,
        algo: &mut DelaunayNodeInsertionMeshAlgo,
        surface: &dyn Surface,
        nodes: &[GpPnt2d],
    ) -> bool {
        // `hxx:108-111`: the early return tests the INPUT list, not the
        // classifier's output.
        if nodes.is_empty() {
            return false;
        }
        let mut idxs = Vec::new();
        let mut out_cls = 0usize;
        for &uv in nodes {
            if !algo.classifier.is_inside(&uv) {
                out_cls += 1;
                continue;
            }
            let p3d = surface.d0(uv.x(), uv.y());
            let vertex = DelaunVertex::new(self.to_face(uv), p3d, 0, VertexState::Free);
            let idx = self.delaun.add_node(vertex);
            if idx as usize > algo.nodes_map.len() {
                algo.nodes_map.push(p3d);
            }
            idxs.push(idx);
        }
        if crate::debug_flags::dc_dbg() {
            eprintln!(
                "DCDBG insert in={} out_cls={} added={}",
                nodes.len(),
                out_cls,
                idxs.len()
            );
        }
        // `hxx:125`: `AddVertices` is called even when every node was
        // classified out, so `ProcessConstraints()` still runs.
        self.delaun.add_vertices(&mut idxs);
        // `hxx:130`: only the RETURN VALUE reflects the classifier's output.
        !idxs.is_empty()
    }
}

impl DeflectionMesh for ParametricDelaun<'_> {
    fn elements_of_domain(&self) -> Vec<i32> {
        self.delaun.elements_of_domain()
    }
    fn triangle(&self, id: i32) -> DelaunTriangle {
        self.delaun.triangle(id)
    }
    fn link(&self, id: i32) -> DelaunLink {
        self.delaun.link(id)
    }
    fn link_movability(&self, id: i32) -> VertexState {
        self.delaun.link_movability(id)
    }
    fn vertex(&self, id: i32) -> DelaunVertex {
        let mut v = self.delaun.vertex(id);
        v.location = self.to_param(v.location);
        v
    }
    fn shot_triangles(&self, p: GpXY) -> Option<Vec<i32>> {
        // `rejectByMinSize` (`hxx:425-427`): `CircleTool::Select(Scale(p, true))`.
        let face = self.to_face(GpPnt2d::from_xy(p));
        Some(self.delaun.circles().select(face.coord))
    }
}

/// Axis-aligned UV bounds of a point set.
pub(super) fn uv_bounds(points: &[GpPnt2d]) -> (f64, f64, f64, f64) {
    let (mut umin, mut umax, mut vmin, mut vmax) =
        (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
    for p in points {
        umin = umin.min(p.x());
        umax = umax.max(p.x());
        vmin = vmin.min(p.y());
        vmax = vmax.max(p.y());
    }
    if !umin.is_finite() {
        return (0.0, 1.0, 0.0, 1.0);
    }
    (umin, umax, vmin, vmax)
}

/// `BRepMesh_GeomTool::CellsCount` + `ComputeErrFactors` + `AdjustCellsCounts`
/// (`BRepMesh_GeomTool.cxx:32-84, 86-144, 465-511`).
pub(super) fn geom_tool_cells_count(
    surface: &dyn Surface,
    vertices_nb: i32,
    deflection: f64,
    range_u: (f64, f64),
    range_v: (f64, f64),
    delta: (f64, f64),
) -> (i32, i32) {
    let kind = mesh_surface_type(surface);
    let (err_u, err_v) = compute_err_factors(surface, kind, deflection);
    let du = range_u.1 - range_u.0;
    let dv = range_v.1 - range_v.0;
    let mut cells_u;
    let mut cells_v;
    if kind == MeshSurfaceType::Torus {
        cells_u = 2.0_f64.powf((du / delta.0).log10()).ceil() as i32;
        cells_v = 2.0_f64.powf((dv / delta.1).log10()).ceil() as i32;
    } else if kind == MeshSurfaceType::Cylinder {
        cells_u = 2.0_f64.powf((du / delta.0 / dv).log10()).ceil() as i32;
        cells_v = 2.0_f64.powf((dv / err_v).log10()).ceil() as i32;
    } else {
        cells_u = 2.0_f64.powf((du / delta.0 / err_u).log10()).ceil() as i32;
        cells_v = 2.0_f64.powf((dv / delta.1 / err_v).log10()).ceil() as i32;
    }
    if kind == MeshSurfaceType::OtherSurface {
        return (-1, -1);
    }
    let logn = 2.0_f64.powf((vertices_nb as f64).log10()).ceil() as i32;
    match kind {
        MeshSurfaceType::Plane => {
            cells_u = logn;
            cells_v = logn;
        }
        MeshSurfaceType::Cylinder | MeshSurfaceType::Cone => {
            cells_v = logn;
        }
        MeshSurfaceType::SurfaceOfExtrusion | MeshSurfaceType::SurfaceOfRevolution => {
            if let Some(c) = basis_curve_for_err(surface, kind) {
                let planar_like = c.is_line()
                    || (c.bspline_poles().is_some()
                        && c.nurbs_degree().map(|d| d < 2).unwrap_or(false));
                if planar_like {
                    if kind == MeshSurfaceType::SurfaceOfExtrusion {
                        cells_u = logn;
                    } else {
                        cells_v = logn;
                    }
                }
            }
            if kind == MeshSurfaceType::SurfaceOfExtrusion {
                cells_v = logn;
            }
        }
        MeshSurfaceType::BezierSurface | MeshSurfaceType::BSplineSurface => {
            if surface.u_degree() < 2 {
                cells_u = logn;
            }
            if surface.v_degree() < 2 {
                cells_v = logn;
            }
        }
        _ => {}
    }
    (cells_u.max(2), cells_v.max(2))
}

/// `ComputeErrFactors` (`BRepMesh_GeomTool.cxx:32-84`).
pub(super) fn compute_err_factors(surface: &dyn Surface, kind: MeshSurfaceType, deflection: f64) -> (f64, f64) {
    let mut err_u = deflection * 10.0;
    let mut err_v = deflection * 10.0;
    match kind {
        MeshSurfaceType::Cylinder | MeshSurfaceType::Cone | MeshSurfaceType::Sphere | MeshSurfaceType::Torus => {}
        MeshSurfaceType::SurfaceOfExtrusion | MeshSurfaceType::SurfaceOfRevolution => {
            if let Some(c) = basis_curve_for_err(surface, kind) {
                if c.bspline_poles().is_some() {
                    if let Some(deg) = c.nurbs_degree() {
                        if deg > 2 {
                            let nk = (c.nb_intervals(0) + 1) as f64;
                            err_v /= deg as f64 * nk;
                        }
                    }
                }
            }
        }
        MeshSurfaceType::BezierSurface => {
            let du = surface.u_degree();
            let dv = surface.v_degree();
            if du > 2 {
                err_u /= du as f64;
            }
            if dv > 2 {
                err_v /= dv as f64;
            }
        }
        MeshSurfaceType::BSplineSurface => {
            let du = surface.u_degree();
            let dv = surface.v_degree();
            if du > 2 {
                err_u /= du as f64 * (surface.nb_u_intervals(0) + 1) as f64;
            }
            if dv > 2 {
                err_v /= dv as f64 * (surface.nb_v_intervals(0) + 1) as f64;
            }
        }
        MeshSurfaceType::Plane | MeshSurfaceType::OffsetSurface | MeshSurfaceType::OtherSurface => {
            err_u = 1.0;
            err_v = 1.0;
        }
    }
    (err_u, err_v)
}

pub(super) fn basis_curve_for_err(
    surface: &dyn Surface,
    kind: MeshSurfaceType,
) -> Option<std::sync::Arc<dyn occt_geom::Curve>> {
    match kind {
        MeshSurfaceType::SurfaceOfExtrusion => surface.extrusion_basis_curve(),
        MeshSurfaceType::SurfaceOfRevolution => surface.revolution_basis_curve(),
        _ => None,
    }
}
