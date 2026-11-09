//! `MapEV` construction: the generated / recovered face behind every edge of
//! the split face (`LocOpe_SplitDrafts.cxx:405-470`).
//!
//! For an edge that lies on the original wire the draft wall is the pipe of that
//! edge along `theLinePipe`; any other edge of the split face borrows the face on
//! its far side. `MapEV` is an indexed map keyed on both edges and (later)
//! vertices, exactly like the reference, whose `FindIndex` doubles as the
//! "already handled" flag.

use std::sync::Arc;

use occt_core::precision::CONFUSION;
use occt_geom::{Curve, Surface};

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::locope::split_drafts::perform::{pipe_surface, PIPE_ERROR};
use crate::locope::split_drafts::state::PerformState;
use crate::locope::split_drafts::SplitDrafts;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same};

impl SplitDrafts {
    /// `cxx:405-470`. `Ok(false)` is the reference's bare `return` on an edge
    /// that does not have exactly two ancestor faces.
    pub(crate) fn build_generated(
        &mut self,
        f: &Face,
        st: &mut PerformState,
    ) -> Result<bool, String> {
        let b = TopoBuilder::new();
        // `cxx:411`: `myMap(F)` - the faces the split produced under `F`.
        let faces: Vec<TopoShape> = self.my_map.find(&f.0).cloned().unwrap_or_default();
        for fac in faces {
            for edg in explorer_shapes(&fac, ShapeType::Edge) {
                // `cxx:417-420`.
                if st.map_ev.find_index(&edg).is_some() {
                    continue;
                }
                if st.map_w.is_bound(&edg) {
                    // `cxx:421-450`: edge of the original wire.
                    let edge = Edge(edg.clone());
                    let mut c: Arc<dyn Curve> = match BRepTool::edge_curve(&edge) {
                        Some(c) => c,
                        None => return Err(PIPE_ERROR.to_string()),
                    };
                    // `cxx:426-430`: unwrap the trim, then place the curve.
                    if let Some(basis) = c.trim_basis() {
                        c = basis;
                    }
                    let loc = edge.0.location();
                    c = Arc::from(c.transformed(&loc.transformation()));

                    let ps = pipe_surface(&st.line_pipe, c)?;
                    // `cxx:441-445`.
                    let ps: Arc<dyn Surface> = match ps.rectangular_trimmed_basis() {
                        Some(basis) => basis,
                        None => ps,
                    };
                    let new_face = b.make_face_with_tolerance(ps, &[], CONFUSION);
                    st.map_ev.add(&edg, new_face.0);
                } else {
                    // `cxx:451-466`: recover the face across the edge.
                    let Some(ancestors) = st.the_map_ef.find_from_key(&edg) else {
                        return Ok(false);
                    };
                    if ancestors.len() != 2 {
                        return Ok(false);
                    }
                    let other = if is_same(&ancestors[0], &fac) {
                        ancestors[1].clone()
                    } else {
                        ancestors[0].clone()
                    };
                    st.map_ev.add(&edg, other);
                }
            }
        }
        Ok(true)
    }
}
