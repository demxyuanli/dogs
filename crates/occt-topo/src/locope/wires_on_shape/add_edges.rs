//! `LocOpe_WiresOnShape::Add(theEdges)` (`LocOpe_WiresOnShape.cxx:1515-1624`).
//!
//! Binds every edge of the given splitting shapes to the face of `myShape` it
//! lies on, without an explicit edge -> face pairing: the candidate face is the
//! first one whose `BRepBndLib::Add` box contains the edge box and whose surface
//! extremum from the edge's midpoint lands inside that face.

use std::collections::HashSet;

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt2d;
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::extrema_surf::{ExtPs, ExtremaExtFlag};

use crate::abs::ShapeType;
use crate::brep_bnd_lib;
use crate::brep_tool::BRepTool;
use crate::brep_top_adaptor_fclass2d::BrepTopAdaptorFClass2d;
use crate::fclass2d::FaceState;
use crate::locope::wires_on_shape::WiresOnShape;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{explorer_shapes, explorer_shapes_avoiding};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

impl WiresOnShape {
    /// `Add(const NCollection_Sequence<TopoDS_Shape>& theEdges)` (`cxx:1515-1624`).
    ///
    /// Returns `true` when at least one (edge, face) pair was bound, exactly as
    /// the reference answers `!anUsedEdges.IsEmpty()`.
    ///
    /// `Err` is the reference's `Standard_ConstructionError` from `Bind`
    /// (`cxx:161-162`), reachable when the same edge is offered twice for the
    /// same face.
    pub fn add(&mut self, the_edges: &[TopoShape]) -> Result<bool, String> {
        let mut an_edges: Vec<Edge> = Vec::new();
        let nb = the_edges.len();
        // `cxx:1519`: `NCollection_Array1<Bnd_Box> anEdgeBoxes(1, nb)` - sized by
        // the *number of input shapes*, yet filled from the inner edge loop, so a
        // shape with several edges overwrites the same slot and `anEdges` ends up
        // longer than the array. Both containers are reproduced as-is.
        let mut an_edge_boxes: Vec<BndBox> = vec![BndBox::new(); nb];

        // `cxx:1520-1539`.
        for (item_index, a_cur_split) in the_edges.iter().enumerate() {
            for shape in explorer_shapes_avoiding(a_cur_split, ShapeType::Edge, ShapeType::Shape) {
                let Some(a_cur_e) = Edge::wrap(shape) else {
                    continue;
                };
                // `cxx:1528-1532`.
                let mut a_box_e = BndBox::new();
                brep_bnd_lib::add_close(&a_cur_e.0, &mut a_box_e);
                if a_box_e.is_void() {
                    continue;
                }
                // `cxx:1533-1538`.
                let a_tol_e = BRepTool::edge_tolerance(&a_cur_e);
                a_box_e.set_gap(a_tol_e);
                an_edge_boxes[item_index] = a_box_e;
                an_edges.push(a_cur_e);
            }
        }

        // `cxx:1540-1622`.
        let mut an_used_edges: HashSet<usize> = HashSet::new();
        for face_shape in explorer_shapes(&self.my_shape, ShapeType::Face) {
            let Some(a_cur_f) = Face::wrap(face_shape) else {
                continue;
            };
            // `cxx:1548-1551`.
            let mut a_box_f = BndBox::new();
            brep_bnd_lib::add(&a_cur_f.0, &mut a_box_f, true);
            if a_box_f.is_void() {
                continue;
            }

            // `cxx:1552`: `BRepAdaptor_Surface anAdF(aCurF, false)` - the
            // `Restriction = false` arm, whose parameters are the surface's own.
            let an_ad_f = SurfaceAdaptor::of_face(&a_cur_f);
            let Some(an_ad_f_surf) = an_ad_f.surface() else {
                continue;
            };
            // `cxx:1555`: `NCollection_Handle<BRepTopAdaptor_FClass2d>
            // aCheckStateTool` is built lazily on the first accepted extremum.
            let mut a_check_state_tool: Option<BrepTopAdaptorFClass2d> = None;

            // `cxx:1558-1566`: `Extrema_ExtPS` initialized once per face over the
            // face's parameter window with `TolU = TolV = Precision::Confusion()`.
            let mut an_extr = ExtPs::new();
            an_extr.initialize(
                an_ad_f_surf.as_ref(),
                an_ad_f.first_u_parameter(),
                an_ad_f.last_u_parameter(),
                an_ad_f.first_v_parameter(),
                an_ad_f.last_v_parameter(),
                CONFUSION,
                CONFUSION,
            );
            // `cxx:1567`: `anExtr.SetFlag(Extrema_ExtFlag_MIN)`.
            an_extr.set_flag(ExtremaExtFlag::Min);

            let nb_boxes = an_edge_boxes.len();
            for i in 1..=nb_boxes {
                // `cxx:1570-1573`.
                if an_used_edges.contains(&i) {
                    continue;
                }
                // `cxx:1575-1578`.
                if a_box_f.is_out_box(&an_edge_boxes[i - 1]) {
                    continue;
                }

                // `cxx:1580-1589`: `anEdges(i)` is the reference's 1-based
                // `NCollection_Sequence::Value`, which raises
                // `Standard_OutOfRange` beyond `Length()`; the flat vector below
                // keeps the same indexing so the same input raises it.
                let Some(a_cur_e) = an_edges.get(i - 1).cloned() else {
                    panic!(
                        "NCollection_Sequence::Value: index {} out of range (Length = {})",
                        i,
                        an_edges.len()
                    );
                };
                let Some(a_c) = BRepTool::edge_curve_world(&a_cur_e) else {
                    // `cxx:1585-1589`.
                    an_used_edges.insert(i);
                    continue;
                };
                let (a_f, a_l) = BRepTool::edge_parameters(&a_cur_e);
                let a_p = a_c.d0((a_f + a_l) * 0.5);
                // `cxx:1590-1595`.
                an_extr.perform(&a_p);
                if !an_extr.is_done() || an_extr.nb_ext() == 0 {
                    continue;
                }

                // `cxx:1597-1621`.
                let a_tol_e = BRepTool::edge_tolerance(&a_cur_e);
                let a_tol2 = (a_tol_e + CONFUSION) * (a_tol_e + CONFUSION);
                for n in 1..=an_extr.nb_ext() {
                    let a_dist2 = an_extr.square_distance(n);
                    if a_dist2 > a_tol2 {
                        continue;
                    }
                    let (a_u, a_v, _) = an_extr.point(n);
                    if a_check_state_tool.is_none() {
                        // `cxx:1612`: `BRepTopAdaptor_FClass2d(aCurF, PConfusion)`.
                        a_check_state_tool = Some(BrepTopAdaptorFClass2d::new(&a_cur_f, PCONFUSION));
                    }
                    let tool = a_check_state_tool.as_ref().expect("built just above");
                    if tool.perform(GpPnt2d::new(a_u, a_v), true) == FaceState::In {
                        self.bind_edge(&a_cur_e, &a_cur_f)?;
                        an_used_edges.insert(i);
                    }
                }
            }
        }

        Ok(!an_used_edges.is_empty())
    }
}
