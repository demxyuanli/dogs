//! `BRepFill_Evolved::AddTopAndBottom` (`BRepFill_Evolved.cxx:2078-2225`) and
//! `MakeSolid` (`cxx:2229-2262`).
//!
//! `AddTopAndBottom` closes the two ends of the sweep: for each end vertex of
//! the profile it puts every edge the sweep generated on that vertex onto the
//! plane `z = -Altitud(V)` through a `BRepAlgo_Loop`, orients those edges
//! against the support edge / the support vertex, and adds the resulting faces
//! to the quilt (`myTop` is the first end, `myBottom` the second).
//!
//! `MakeSolid` wraps every shell of `myShape` in a solid, reversing the shells
//! that classify `IN` at the infinite point.

use occt_core::gp::GpVec;

use crate::abs::{Orientation, ShapeType};
use crate::brep_algo_loop::BRepAlgoLoop;
use crate::brep_class3d::SolidClassifier;
use crate::brep_sweep::Quilt;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Shell, TopoShape};
use crate::topo_tools::{edges_of, explorer_shapes, shape_map_key, vertices_of, wire_vertices};
use crate::topopebrep::tool::adaptor::CurveAdaptor;

use super::evolved::BRepFillEvolved;
use super::side::{altitud, complement, FILL_CONFUSION};

impl BRepFillEvolved {
    /// `AddTopAndBottom(Glue)` (`cxx:2078-2225`).
    pub fn add_top_and_bottom(&mut self, glue: &mut Quilt) {
        // `cxx:2080-2087`.
        let Some(profile) = self.my_profile.clone() else {
            return;
        };
        let (v0, v1) = wire_vertices(&profile);
        let (Some(v0), Some(v1)) = (v0, v1) else {
            return;
        };
        if crate::topo_tools::is_same(&v0.0, &v1.0) {
            return;
        }

        let Some(my_spine) = self.my_spine.clone() else {
            return;
        };
        let builder = TopoBuilder::new();
        // `cxx:2089`: `ToReverse` lives across both ends.
        let mut to_reverse = false;

        for i in 0..2 {
            let vi = if i == 0 { &v0 } else { &v1 };

            // `cxx:2093-2097`: the support face `gp_Pln(0, 0, 1, -Altitud(V))`.
            let mut loop_builder = BRepAlgoLoop::new();
            let pln = super::side::plane_at_z(-altitud(vi));
            let face = builder.make_face_plane(&pln);
            loop_builder.init(&face);

            // `cxx:2099-2138`: the edges generated on the spine edges.
            let mut view: Vec<crate::topo_tools::ShapeMapKey> = Vec::new();
            for es in edges_of(&my_spine.0) {
                let generated = self.generated_shapes(&es.0, &vi.0);
                let mut compute_orientation = false;
                for e_shape in generated {
                    let e = Edge(e_shape);
                    if !compute_orientation {
                        let c1 = CurveAdaptor::of_edge(&es);
                        let c2 = CurveAdaptor::of_edge(&e);
                        let (f, l) = BRepTool::edge_parameters(&e);
                        let (fs, ls) = BRepTool::edge_parameters(&es);
                        let u = 0.3 * f + 0.7 * l;
                        let us = 0.3 * fs + 0.7 * ls;
                        let (Some((_, v1_dir)), Some((_, v2_dir))) = (c1.d1(us), c2.d1(u))
                        else {
                            continue;
                        };
                        to_reverse = v1_dir.dot(&v2_dir) < 0.0;
                        compute_orientation = true;
                    }

                    // `cxx:2129-2134`.
                    let mut or = es.0.orientation();
                    if to_reverse {
                        or = complement(or);
                    }
                    loop_builder.add_const_edge(&Edge(e.0.oriented(or)));
                }
            }

            // `cxx:2140-2185`: the edges generated on the spine vertices.
            let pv = BRepTool::vertex_point(vi);
            let is_out = pv.y() < 0.0;
            for es in vertices_of(&my_spine.0) {
                let key = shape_map_key(&es.0);
                if view.contains(&key) {
                    continue;
                }
                view.push(key);

                let generated = self.generated_shapes(&es.0, &vi.0);
                for e_shape in generated {
                    let e = Edge(e_shape);
                    // `cxx:2152`.
                    if BRepTool::is_degenerated(&e) {
                        continue;
                    }
                    let c = CurveAdaptor::of_edge(&e);
                    let (f, l) = BRepTool::edge_parameters(&e);
                    let u = 0.3 * f + 0.7 * l;
                    let p = BRepTool::vertex_point(&es);
                    let Some((pc, vc)) = c.d1(u) else {
                        continue;
                    };
                    let appc = GpVec::new(pc.x() - p.x(), pc.y() - p.y(), pc.z() - p.z());
                    let prod = appc.crossed(&vc);
                    to_reverse = if is_out {
                        prod.z() < 0.0
                    } else {
                        prod.z() > 0.0
                    };
                    let or = if to_reverse {
                        Orientation::Reversed
                    } else {
                        Orientation::Forward
                    };
                    loop_builder.add_const_edge(&Edge(e.0.oriented(or)));
                }
            }

            // `cxx:2187-2188`.
            loop_builder.perform();
            loop_builder.wires_to_faces();
            let new_faces: Vec<TopoShape> = loop_builder.new_faces().to_vec();

            // `cxx:2192-2221`.
            let mut bouchon = builder.make_compound_of(&[]).0;
            let mut j = 0;
            for face_shape in &new_faces {
                j += 1;
                glue.add(face_shape);
                if j == 1 && i == 0 {
                    self.my_top = Some(face_shape.clone());
                }
                if j == 1 && i == 1 {
                    self.my_bottom = Some(face_shape.clone());
                }
                builder.add(&mut bouchon, face_shape);
            }
            if i == 0 && j > 1 {
                self.my_top = Some(bouchon.clone());
            }
            if i == 1 && j > 1 {
                self.my_bottom = Some(bouchon);
            }
        }
    }

    /// `MakeSolid()` (`cxx:2229-2262`).
    pub fn make_solid(&mut self) {
        let builder = TopoBuilder::new();
        let Some(shape) = self.my_shape.clone() else {
            return;
        };

        // `cxx:2231-2238`.
        let mut res = builder.make_compound_of(&[]).0;
        let mut ish = 0;
        let mut sol: Option<TopoShape> = None;

        // `cxx:2239-2256`.
        for sh in explorer_shapes(&shape, ShapeType::Shell) {
            // `cxx:2243-2244`: `B.MakeSolid(Sol); B.Add(Sol, Sh)`.
            let mut current = builder.make_solid(&[Shell(sh.clone())]).0;
            let mut classifier = SolidClassifier::new();
            classifier.load(current.clone());
            classifier.perform_infinite_point(FILL_CONFUSION);
            if classifier.state() == FaceState::In {
                current = builder
                    .make_solid(&[Shell(sh.oriented(complement(sh.orientation())))])
                    .0;
            }
            builder.add(&mut res, &current);
            sol = Some(current);
            ish += 1;
        }

        // `cxx:2256-2261`.
        if ish == 1 {
            if let Some(sol) = sol {
                self.my_shape = Some(sol);
            }
        } else {
            self.my_shape = Some(res);
        }
    }
}
