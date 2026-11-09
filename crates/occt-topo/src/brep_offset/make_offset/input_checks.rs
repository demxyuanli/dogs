//! `BRepOffset_MakeOffset::CheckInputData` (`BRepOffset_MakeOffset.cxx:
//! 4384-4498`) and its file-local helper `checkSinglePoint` (`:4583-4622`).
//!
//! `CheckInputData` runs four checks in the reference order:
//!
//! 1. `:4391-4414` non-null offset: `|myOffset| > myTol`, or at least one
//!    `myFaceOffset` value above `myTol`, otherwise
//!    `BRepOffset_NullOffset`.
//! 2. `:4417-4421` connectivity of `myFaceComp` (`IsConnectedShell`), else
//!    `BRepOffset_NotConnectedShell`.
//! 3. `:4446-4452` per distinct `TShape` of `myFaceComp` (`aPresenceMap`),
//!    `Geom_Surface::Continuity() == GeomAbs_C0` is
//!    `BRepOffset_C0Geometry`.
//! 4. `:4455-4497` a 21x21 grid over `BRepTools::UVBounds` plus every vertex'
//!    `BRep_Tool::Parameters` value, each through `checkSinglePoint`, which
//!    reports `BRepOffset_BadNormalsOnGeometry` when a partial derivative
//!    vanishes away from a known degenerated point or when the two partials
//!    are collinear.
//!
//! `myBadShape` is reset at `:4388-4389` and `checkSinglePoint` never binds it
//! in the reference either (the `GetBadShape` doc comment is the only place
//! that claims otherwise), so it stays null.
//!
//! UNPORTED:
//!
//! * the `Message_ProgressScope` guard `!aPS.More()`
//!   (`BRepOffset_MakeOffset.cxx:4430-4434`): no progress machinery in the
//!   port, same as the rest of `BRepOffset`.
//!
//! Port guards, each replacing an OCCT exception:
//!
//! * `BRep_Tool::Surface(aF, L)` (`:4445`) answers a null handle for a face
//!   without a surface, and the following `Continuity()` faults; the port
//!   reports `BRepOffset_C0Geometry` (the only error that check can raise).
//! * `BRep_Tool::Parameters(aV, aF)` (`:4488`) throws `Standard_NoSuchObject`
//!   when the vertex has no parameter on the face; the port skips that vertex.

use std::collections::HashSet;

use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{CONFUSION, SQUARE_CONFUSION};
use occt_geom::Surface;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::Face;
use crate::topo_tools::{edge_vertices, explorer_shapes};

use super::shape_prep::is_connected_shell;
use super::BRepOffsetMakeOffset;
use crate::brep_offset::Error;

/// `checkSinglePoint(theUParam, theVParam, theSurf, theBadPoints)`
/// (`cxx:4583-4622`).
pub fn check_single_point(u: f64, v: f64, surf: &dyn Surface, bad_points: &[GpPnt]) -> Error {
    let (pnt, d1u, d1v) = surf.d1(u, v);

    // `cxx:4592-4615`: a vanishing partial is only allowed at a point of a
    // degenerated edge.
    if d1u.square_magnitude() < SQUARE_CONFUSION || d1v.square_magnitude() < SQUARE_CONFUSION {
        let is_known_bad = bad_points
            .iter()
            .any(|bad| pnt.square_distance(bad) < SQUARE_CONFUSION);
        if is_known_bad {
            return Error::NoError;
        }
        return Error::BadNormalsOnGeometry;
    }

    // `cxx:4617-4621`: collinear iso-lines.
    if d1u.is_parallel_ang(&d1v, CONFUSION) {
        return Error::BadNormalsOnGeometry;
    }

    Error::NoError
}

impl BRepOffsetMakeOffset {
    /// `CheckInputData(theRange)` (`cxx:4384-4498`).
    pub fn check_input_data(&mut self) -> bool {
        // `cxx:4386-4389`.
        self.my_error = Error::NoError;
        self.my_bad_shape = None;

        // `cxx:4391-4414`: the offset must not be null, unless one face has an
        // own non-null offset.
        if self.my_offset.abs() <= self.my_tol {
            let mut is_found = false;
            for (_, value) in self.my_face_offset.iter() {
                if value.abs() > self.my_tol {
                    is_found = true;
                    break;
                }
            }
            if !is_found {
                self.my_error = Error::NullOffset;
                return false;
            }
        }

        // `cxx:4417-4421`.
        if !is_connected_shell(&self.my_face_comp.0) {
            self.my_error = Error::NotConnectedShell;
            return false;
        }

        // `cxx:4423-4425`.
        const A_PNT_PER_DIM: usize = 20;

        // `cxx:4426`: `NCollection_Map<Handle(TopoDS_TShape)>`, i.e. the
        // partner shapes that share the same geometry are visited once.
        let mut a_presence_map: HashSet<usize> = HashSet::new();

        for face_shape in explorer_shapes(&self.my_face_comp.0, ShapeType::Face) {
            let Some(a_f) = Face::wrap(face_shape) else {
                continue;
            };
            if !a_presence_map.insert(a_f.0.shape_addr()) {
                continue;
            }

            // `cxx:4445-4447`: `BRep_Tool::Surface(aF, L)` + `BRepTools::UVBounds`.
            let Some(a_surf) = BRepTool::face_surface(&a_f) else {
                self.my_error = Error::C0Geometry;
                return false;
            };
            let (a_umin, a_umax, a_vmin, a_vmax) = BRepTool::uv_bounds(&a_f);

            // `cxx:4449-4453`.
            if a_surf.continuity() == Shape::C0 as u8 {
                self.my_error = Error::C0Geometry;
                return false;
            }

            // `cxx:4455-4465`: the points of the degenerated edges, excluded
            // from the bad-normal report.
            let mut a_bad3d_pnts: Vec<GpPnt> = Vec::new();
            for edge_shape in explorer_shapes(&a_f.0, ShapeType::Edge) {
                let e = crate::shape::Edge(edge_shape);
                if BRepTool::is_degenerated(&e) {
                    let (v_first, _) = edge_vertices(&e);
                    if let Some(v_first) = v_first {
                        a_bad3d_pnts.push(BRepTool::vertex_point(&v_first));
                    }
                }
            }

            // `cxx:4468-4481`: the geometry grid.
            for i in 0..=A_PNT_PER_DIM {
                let a_u_param =
                    a_umin + (a_umax - a_umin) * i as f64 / A_PNT_PER_DIM as f64;
                for j in 0..=A_PNT_PER_DIM {
                    let a_v_param =
                        a_vmin + (a_vmax - a_vmin) * j as f64 / A_PNT_PER_DIM as f64;
                    self.my_error = check_single_point(
                        a_u_param,
                        a_v_param,
                        a_surf.as_ref(),
                        &a_bad3d_pnts,
                    );
                    if self.my_error != Error::NoError {
                        return false;
                    }
                }
            }

            // `cxx:4484-4496`: the parameters of the face's vertices.
            for v_shape in explorer_shapes(&a_f.0, ShapeType::Vertex) {
                let Ok(a_pnt2d) = BRepTool::parameters_of_vertex_on_face(&v_shape, &a_f) else {
                    continue;
                };
                self.my_error = check_single_point(
                    a_pnt2d.x(),
                    a_pnt2d.y(),
                    a_surf.as_ref(),
                    &a_bad3d_pnts,
                );
                if self.my_error != Error::NoError {
                    return false;
                }
            }
        }

        true
    }
}
