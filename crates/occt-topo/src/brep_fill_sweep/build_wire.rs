//! `BRepFill_Sweep::BuildWire` (`BRepFill_Sweep.cxx:2078-2210`): the sweep of a
//! vertex section, i.e. the generating wire itself.
//!
//! Translation notes:
//! - `GeomFill_Sweep` is [`occt_geom::geom_fill::sweep::Sweep`];
//!   `Sweep.Surface()` answers `Option<Arc<dyn Surface>>` where OCCT leaves a
//!   null handle, which is `false` here.
//! - `NCollection_HArray2<TopoDS_Shape>` is [`ShapeTable`];
//!   `myUEdges(1,2,1,NbPath)` is `ShapeTable::new(2, NbPath)` and
//!   `myVEdges(1,1,1,NbPath+1)` is `ShapeTable::new(1, NbPath + 1)`.
//! - `TopoDS::Vertex(shape)` is an unchecked cast in OCCT. The port only takes
//!   the slot when it really holds a vertex; `FirstShape` / `LastShape` put a
//!   *wire* into that table (`cxx:2163`, `cxx:2184`), and dereferencing it as a
//!   vertex (`BRep_Tool::Pnt`) is not reproducible, so that arm would be
//!   undefined behaviour in the reference as well.

use occt_geom::curve::Curve;
use occt_geom::geom_fill::sweep::Sweep as GeomFillSweep;

use crate::brep_lib_make_edge::{add_endpoint_vertices, init_3d_vertices};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::abs::ShapeType;
use crate::shape::{Vertex, Wire};
use crate::tgeometry::GeometryRegistry;

use super::helpers::{self, ShapeTable};
use super::{BrepFillSweep, TransitionStyle};

/// `IsG1(0, myTol3d)`'s default `AngularTolerance`
/// (`BRepFill_LocationLaw.hxx:102`).
const DEFAULT_ANGULAR_TOL: f64 = 1.0e-4;

impl BrepFillSweep {
    /// `BRepFill_Sweep::BuildWire(Transition)` (`cxx:2078-2210`).
    pub fn build_wire(&mut self, _transition: TransitionStyle) -> bool {
        let reg = GeometryRegistry::global();
        let b = TopoBuilder::new();
        let nb_path = self.my_loc.base().nb_law();
        let vclose = self.my_loc.base().is_closed()
            && self
                .my_loc
                .base()
                .is_g1(0, self.my_tol_3d, DEFAULT_ANGULAR_TOL)
                >= 0;
        self.error = 0.0;

        // `cxx:2103-2106`.
        self.my_faces = ShapeTable::new(1, nb_path);
        self.my_u_edges = ShapeTable::new(2, nb_path);
        self.my_v_edges = ShapeTable::new(1, nb_path + 1);

        // `cxx:2094-2097`: one wire for the whole sweep.
        let mut wire = Wire::new();

        // `cxx:2108-2208`.
        for ipath in 1..=nb_path {
            // `cxx:2109-2112`.
            let Some(law) = self.my_loc.base().law(ipath) else {
                return false;
            };
            let mut sweep = GeomFillSweep::with_k_part(law, self.k_part);
            sweep.set_tolerance(
                self.my_tol_3d,
                self.my_bound_tol,
                self.my_tol_2d,
                self.my_tol_angular,
            );
            sweep.set_force_approx_c1(self.my_force_approx_c1);
            let Some(section) = self.my_sec.base().law(1) else {
                return false;
            };
            sweep.build(
                section,
                self.my_approx_style,
                self.my_continuity,
                self.my_degmax,
                self.my_segmax,
            );
            if !sweep.is_done() {
                return false;
            }
            let Some(s) = sweep.surface() else {
                return false;
            };

            // `cxx:2118-2141`.
            let (iso, first, last): (std::sync::Arc<dyn Curve>, f64, f64) =
                if sweep.exchange_uv() {
                    let (u_first, u_last) = s.u_range();
                    let (v_first, v_last) = s.v_range();
                    let val = if sweep.u_reversed() { v_last } else { v_first };
                    let Some(iso) = s.v_iso_curve(val) else {
                        return false;
                    };
                    (iso, u_first, u_last)
                } else {
                    let (u_first, u_last) = s.u_range();
                    let (v_first, v_last) = s.v_range();
                    let val = if sweep.u_reversed() { u_last } else { u_first };
                    let Some(iso) = s.u_iso_curve(val) else {
                        return false;
                    };
                    (iso, v_first, v_last)
                };

            // `cxx:2143-2170`: the vertex at the path end.
            if ipath < nb_path {
                self.build_vertex_into(&iso, false, first, last, 1, ipath + 1);
            } else if vclose {
                // `cxx:2150-2158`.
                if let Some(v) = vertex_in(&self.my_v_edges, 1, 1) {
                    self.my_v_edges
                        .set_value(1, ipath + 1, Some(v.0.clone()));
                    let p1 = iso.d0(last);
                    let tol = p1.distance(&BRepTool::vertex_point(&v));
                    reg.raise_vertex_tolerance(&v.0, tol);
                }
            } else if self.last_shape.is_some() {
                // `cxx:2161-2164`: the reference stores `FirstShape` (a wire)
                // into the vertex table; see the module note.
                let stored = self.first_shape.as_ref().map(|w| w.0.clone());
                if let Some(shape) = stored {
                    self.my_v_edges.set_value(1, nb_path, Some(shape));
                }
            } else {
                self.build_vertex_into(&iso, false, first, last, 1, nb_path + 1);
            }

            // `cxx:2172-2178`.
            if ipath > 1 {
                let p1 = iso.d0(first);
                if let Some(v) = vertex_in(&self.my_v_edges, 1, ipath) {
                    let tol = p1.distance(&BRepTool::vertex_point(&v));
                    reg.raise_vertex_tolerance(&v.0, tol);
                }
            }
            // `cxx:2179-2188`.
            if ipath == 1 {
                if self.first_shape.is_some() {
                    let stored = self.first_shape.as_ref().map(|w| w.0.clone());
                    if let Some(shape) = stored {
                        self.my_v_edges.set_value(1, 1, Some(shape));
                    }
                } else {
                    self.build_vertex_into(&iso, true, first, last, 1, 1);
                }
            }

            // `cxx:2191-2207`: `BRepLib_MakeEdge::Init(Iso, V1, V2, F, L)`.
            let (Some(v1), Some(v2)) =
                (vertex_in(&self.my_v_edges, 1, ipath), vertex_in(&self.my_v_edges, 1, ipath + 1))
            else {
                // The table slots hold a wire (`cxx:2163`); OCCT's `TopoDS::Vertex`
                // cast on it is not reproducible, see the module note.
                return false;
            };
            let (f, l) = (iso.first_parameter(), iso.last_parameter());
            let Some(build) = init_3d_vertices(&iso, Some(&v1), Some(&v2), f, l) else {
                return false;
            };
            let mut e = b.make_edge(iso.clone(), build.first, build.last);
            add_endpoint_vertices(&b, &mut e, &build);
            reg.raise_edge_tolerance(&e.0, sweep.error_on_surface());

            // `cxx:2205-2207`.
            b.add_edge(&mut wire, &e);
            self.my_faces.set_value(1, ipath, Some(e.0.clone()));
            self.my_shape = Some(wire.0.clone());
        }

        true
    }

    /// `BuildVertex(Iso, isfirst, First, Last, V)` (`cxx:1440-1458`) writing
    /// into `table(row, col)`.
    fn build_vertex_into(
        &mut self,
        iso: &std::sync::Arc<dyn Curve>,
        isfirst: bool,
        first: f64,
        last: f64,
        row: usize,
        col: usize,
    ) {
        let mut v = Vertex::new();
        helpers::build_vertex(iso, isfirst, first, last, &mut v);
        self.my_v_edges.set_value(row, col, Some(v.0));
    }
}

/// `TopoDS::Vertex(shape)` guarded by the actual shape type.
fn vertex_in(table: &ShapeTable, row: usize, col: usize) -> Option<Vertex> {
    let shape = table.value(row, col)?;
    if shape.shape_type() == ShapeType::Vertex {
        Some(Vertex(shape.clone()))
    } else {
        None
    }
}
