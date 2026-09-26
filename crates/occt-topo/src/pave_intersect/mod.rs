//! Pairwise intersection of the PaveFiller — port of `BOPAlgo_PaveFiller_{1..4}.cxx`.
//!
//! This module contains the *intersection execution* of the Boolean component's
//! PaveFiller (Phase 19, wave C2b). The [`crate::pave_filler::PaveFiller`] owns
//! the data structure ([`BopdsDS`]) and the pipeline; the heavy per-pair
//! geometric work lives here:
//!
//! | Stage          | Source                           | Port                                                        |
//! |----------------|----------------------------------|-------------------------------------------------------------|
//! | Vertex/Vertex  | `PerformVV`                      | [`perform_vv`] — fuse coincident vertices into SD vertices  |
//! | Vertex/Edge    | `PerformVE` / `IntersectVE`      | [`perform_ve`] — project vertices onto edges, insert paves  |
//! | Edge/Edge      | `PerformEE` + `TreatNewVertices` | [`perform_ee`] — edge/edge intersection via [`EdgeEdge`]     |
//! | Vertex/Face    | `PerformVF`                      | [`perform_vf`] — classify vertices against faces            |
//! | Split          | `SplitPaveBlocks`                | [`split_pave_blocks`] — split pave blocks with extra paves  |
//!
//! The module works entirely on the DS through two handles:
//!
//! * [`BopdsDS`] — the shape registry + pave-block pool (read via
//!   [`crate::pave_filler::PaveFiller::ds`], write via `ds_mut`);
//! * [`FillCtx`] — the per-run options the intersection routines need
//!   (fuzzy value, glue mode, non-destructive flag) plus the error report.
//!
//! Every public function returns `Result<(), String>` and, on a hard failure,
//! records a fatal alert on the filler via
//! [`crate::pave_filler::PaveFiller::add_error`].
mod prelude {

pub(crate) use std::collections::{HashMap, HashSet};


pub(crate) use occt_core::gp::{GpPnt, GpPnt2d, GpVec};
pub(crate) use occt_core::precision::{CONFUSION, PCONFUSION, RESOLUTION};
pub(crate) use occt_geom::Curve;

pub(crate) use crate::abs::ShapeType;
pub(crate) use crate::algo_tools::{AlgoTools, D_TOLERANCE};
pub(crate) use crate::bopds::{
    BopdsCommonBlock, BopdsDS, BopdsFaceInfo, BopdsIterator, BopdsPave, BopdsPaveBlock,
    BopdsShapeInfo,
};
pub(crate) use crate::boptools_2d::intermediate_point;
pub(crate) use crate::bopds_ff::{BopdsCurve, BopdsInterfFf};
pub(crate) use crate::brep_surface::face_is_planar;
pub(crate) use crate::brep_tool::BRepTool;
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::edge_edge::EdgeEdge;
pub(crate) use crate::edge_face::EdgeFace;
pub(crate) use crate::fclass2d::FaceState;
pub(crate) use crate::int_face_face::FaceFace;
pub(crate) use crate::int_tools_full::IntToolsContext;
pub(crate) use crate::inttools_data::{CommonPartType, IntRange};
pub(crate) use crate::pave_filler::{EdgeRangeDistance, GlueEnum, PaveFiller};
pub(crate) use crate::shape::{Edge, Face, TopoShape, Vertex};
pub(crate) use crate::topo_tools_full::vertices_of;

}

use prelude::*;


// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep_extrema::test_box::unit_box;
    use occt_core::gp::{GpAx3, GpDir, GpLin, GpVec};

    /// A unit box with its corner at `(ox, oy, oz)` (fresh shapes, so the DS
    /// sees two distinct arguments even when their corners coincide).
    struct BoxAt {
        solid: crate::shape::Solid,
        vertices: Vec<Vertex>,
        #[allow(dead_code)] // edge/face lists kept for future fixture use
        edges: Vec<Edge>,
        #[allow(dead_code)]
        faces: Vec<Face>,
    }

    /// Build a unit box with the corner at the given offset.
    fn box_at(ox: f64, oy: f64, oz: f64) -> BoxAt {
        let b = TopoBuilder::new();
        let corners = [
            GpPnt::new(ox, oy, oz),
            GpPnt::new(1.0 + ox, oy, oz),
            GpPnt::new(1.0 + ox, 1.0 + oy, oz),
            GpPnt::new(ox, 1.0 + oy, oz),
            GpPnt::new(ox, oy, 1.0 + oz),
            GpPnt::new(1.0 + ox, oy, 1.0 + oz),
            GpPnt::new(1.0 + ox, 1.0 + oy, 1.0 + oz),
            GpPnt::new(ox, 1.0 + oy, 1.0 + oz),
        ];
        let vertices: Vec<Vertex> = corners.iter().map(|p| b.make_vertex(*p, 1e-7)).collect();
        let edge_idx: [(usize, usize); 12] = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        let mut edges = Vec::new();
        for &(i, j) in &edge_idx {
            let p1 = corners[i];
            let p2 = corners[j];
            let dir = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)).expect("distinct");
            let lin = GpLin::from_pnt_dir(p1, dir);
            let mut e = b.make_edge(
                std::sync::Arc::new(occt_geom::GeomLine::new(lin)),
                0.0,
                p1.distance(&p2),
            );
            b.add_edge_vertices(&mut e, &vertices[i], &vertices[j]);
            edges.push(e);
        }
        // Face planes: the same normal/u axes as the unit box, origins shifted.
        let face_planes: [(GpPnt, GpDir, GpDir); 6] = [
            (GpPnt::new(ox, oy, oz), GpDir::new(0.0, 0.0, -1.0).unwrap(), GpDir::new(0.0, 1.0, 0.0).unwrap()),
            (GpPnt::new(ox, oy, 1.0 + oz), GpDir::new(0.0, 0.0, 1.0).unwrap(), GpDir::new(1.0, 0.0, 0.0).unwrap()),
            (GpPnt::new(ox, oy, oz), GpDir::new(0.0, -1.0, 0.0).unwrap(), GpDir::new(1.0, 0.0, 0.0).unwrap()),
            (GpPnt::new(ox, 1.0 + oy, oz), GpDir::new(0.0, 1.0, 0.0).unwrap(), GpDir::new(0.0, 0.0, 1.0).unwrap()),
            (GpPnt::new(ox, oy, oz), GpDir::new(-1.0, 0.0, 0.0).unwrap(), GpDir::new(0.0, 0.0, 1.0).unwrap()),
            (GpPnt::new(1.0 + ox, oy, oz), GpDir::new(1.0, 0.0, 0.0).unwrap(), GpDir::new(0.0, 1.0, 0.0).unwrap()),
        ];
        let face_edge_sets: [[usize; 4]; 6] = [
            [0, 1, 2, 3],
            [4, 5, 6, 7],
            [0, 9, 4, 8],
            [2, 10, 6, 11],
            [3, 11, 7, 8],
            [1, 10, 5, 9],
        ];
        let mut faces = Vec::new();
        for fi in 0..6 {
            let (origin, normal, u_dir) = face_planes[fi];
            let ax3 = GpAx3::new(origin, normal, &u_dir).expect("perpendicular axes");
            let mut face = b.make_face_plane(&occt_core::gp::GpPln::new(ax3));
            let wire = b.make_wire(&face_edge_sets[fi].map(|ei| edges[ei].clone()));
            b.add_wire(&mut face, &wire);
            faces.push(face);
        }
        let shell = b.make_shell(&faces);
        let solid = b.make_solid(&[shell]);
        BoxAt { solid, vertices, edges, faces }
    }

    /// Initialize the pave blocks of every source edge of `ds`.
    fn init_all_edge_blocks(ds: &mut BopdsDS) {
        let n = ds.nb_source_shapes();
        for i in 0..n {
            if ds.shape_info(i).map(|s| s.shape_type()) == Some(ShapeType::Edge) {
                ds.init_pave_blocks_for_edge(i);
            }
        }
    }

    #[test]
    fn merge_vertices_averages_points_and_takes_max_tolerance() {
        let b = TopoBuilder::new();
        let v1 = b.make_vertex(GpPnt::new(0.0, 0.0, 0.0), 1e-7);
        let v2 = b.make_vertex(GpPnt::new(0.0, 0.0, 0.0), 2e-7);
        let merged = merge_vertices(&[v1.0.clone(), v2.0.clone()]).unwrap();
        let mv = Vertex(merged);
        let p = BRepTool::vertex_point(&mv);
        assert!((p.x().abs() < 1e-9) && (p.y().abs() < 1e-9) && (p.z().abs() < 1e-9));
        assert!(BRepTool::vertex_tolerance(&mv) >= 2e-7);
    }

    #[test]
    fn make_sd_vertices_links_component_and_interferes() {
        let b = TopoBuilder::new();
        let v1 = b.make_vertex(GpPnt::new(0.0, 0.0, 0.0), 1e-7);
        let v2 = b.make_vertex(GpPnt::new(0.0, 0.0, 0.0), 1e-7);
        let mut ds = BopdsDS::new();
        let n1 = ds.append(v1.into()).unwrap();
        let n2 = ds.append(v2.into()).unwrap();
        let nv = make_sd_vertices(&mut ds, &[n1, n2], true).unwrap();
        assert_eq!(ds.get_same_domain_index(n1), nv);
        assert_eq!(ds.get_same_domain_index(n2), nv);
        assert!(ds.has_interf_pair(n1, n2));
    }

    #[test]
    fn perform_vv_merges_shared_box_corners() {
        let a = unit_box();
        let bbox = box_at(1.0, 0.0, 0.0);
        let mut f = PaveFiller::new();
        f.set_arguments(&[a.solid.0.clone(), bbox.solid.0.clone()]);
        f.init().unwrap();
        // The 4 corners of the shared face (x = 1) coincide.
        let idx_a = [
            a.vertices[1].0.clone(), // (1,0,0)
            a.vertices[2].0.clone(), // (1,1,0)
            a.vertices[5].0.clone(), // (1,0,1)
            a.vertices[6].0.clone(), // (1,1,1)
        ];
        let idx_b = [
            bbox.vertices[0].0.clone(), // (1,0,0)
            bbox.vertices[3].0.clone(), // (1,1,0)
            bbox.vertices[4].0.clone(), // (1,0,1)
            bbox.vertices[7].0.clone(), // (1,1,1)
        ];
        perform_vv(&mut f).unwrap();
        for i in 0..4 {
            let na = f.ds().index(&idx_a[i]).unwrap();
            let nb = f.ds().index(&idx_b[i]).unwrap();
            assert!(
                f.ds().has_shape_sd(na).is_some() || f.ds().has_shape_sd(nb).is_some(),
                "corner pair {i} merged into an SD vertex"
            );
            assert_eq!(
                f.ds().get_same_domain_index(na),
                f.ds().get_same_domain_index(nb),
                "corner pair {i} share the same SD vertex"
            );
        }
        // A non-coincident corner pair is not merged.
        let n0a = f.ds().index(&a.vertices[0].0).unwrap(); // (0,0,0)
        let n0b = f.ds().index(&bbox.vertices[1].0).unwrap(); // (2,0,0)
        assert_ne!(f.ds().get_same_domain_index(n0a), f.ds().get_same_domain_index(n0b));
    }

    #[test]
    fn vertex_on_edge_finds_interior_parameter() {
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(2.0, 0.0, 0.0));
        let v = b.make_vertex(GpPnt::new(0.5, 0.0, 0.0), 1e-7);
        let ctx = IntToolsContext::new();
        let (t, tol) = vertex_on_edge(&ctx, &v, &e, 1e-7).unwrap().unwrap();
        assert!((t - 0.5).abs() < 1e-6, "t = {t}");
        assert!(tol >= 1e-7);
        // A vertex separated from the edge is not found.
        let far = b.make_vertex(GpPnt::new(0.5, 5.0, 0.0), 1e-7);
        assert!(vertex_on_edge(&ctx, &far, &e, 1e-7).unwrap().is_none());
    }

    #[test]
    fn perform_ve_paves_vertex_on_edge() {
        // A = [0,1]^3, B shifted by (0.5, 0, 0). B's corner (0.5, 0, 0) lies on
        // A's bottom-front edge (0,0,0)-(1,0,0), so perform_ve must insert a
        // pave at t = 0.5 and split A's edge into two blocks.
        let a = unit_box();
        let bbox = box_at(0.5, 0.0, 0.0);
        let mut f = PaveFiller::new();
        f.set_arguments(&[a.solid.0.clone(), bbox.solid.0.clone()]);
        f.init().unwrap();
        init_all_edge_blocks(f.ds_mut());
        perform_ve(&mut f).unwrap();

        let na_e = f.ds().index(&a.edges[0].0).unwrap();
        let nb_v = f.ds().index(&bbox.vertices[0].0).unwrap(); // (0.5, 0, 0)
        // The edge was split into two blocks by the interior pave.
        let blocks = f.ds().pave_blocks(na_e);
        assert_eq!(blocks.len(), 2, "edge split at the vertex projection");
        // The interior block references B's corner vertex.
        let mid = blocks
            .iter()
            .find(|pb| pb.pave1().index() == nb_v || pb.pave2().index() == nb_v)
            .expect("a block is bounded by B's corner vertex");
        let (t1, t2) = mid.range();
        assert!(t1.abs() < 1e-6 || (t2 - 1.0).abs() < 1e-6, "range {t1}..{t2} touches the split");
        // The V/E interference was recorded.
        assert!(f.ds().has_interf_pair(nb_v, na_e));
    }

    #[test]
    fn perform_ee_creates_vertex_at_crossing() {
        // A = [0,1]^3, B shifted by (0.5, 0.5, 0). A's right-bottom edge
        // (1,0,0)-(1,1,0) crosses B's bottom-front edge (0.5,0.5,0)-(1.5,0.5,0)
        // at (1, 0.5, 0), interior to both.
        let a = unit_box();
        let bbox = box_at(0.5, 0.5, 0.0);
        let mut f = PaveFiller::new();
        f.set_arguments(&[a.solid.0.clone(), bbox.solid.0.clone()]);
        f.init().unwrap();
        init_all_edge_blocks(f.ds_mut());
        perform_ee(&mut f).unwrap();

        // A new vertex close to (1, 0.5, 0) must be present among the new
        // shapes appended after the source shapes.
        let n_src = f.ds().nb_source_shapes();
        let mut found = false;
        for i in n_src..f.ds().nb_shapes() {
            if f.ds().shape_info(i).map(|s| s.shape_type()) != Some(ShapeType::Vertex) {
                continue;
            }
            let p = BRepTool::vertex_point(&Vertex(f.ds().shape(i).unwrap().clone()));
            if p.distance(&GpPnt::new(1.0, 0.5, 0.0)) < 1e-6 {
                found = true;
                break;
            }
        }
        assert!(found, "crossing vertex created in the DS");

        // The crossing edge of A is split at the crossing: the block carrying
        // the extra pave was replaced by elementary blocks (OCCT
        // `SplitPaveBlocks` at the end of `IntersectEE`).
        let na_e = f.ds().index(&a.edges[1].0).unwrap();
        let blocks = f.ds().pave_blocks(na_e);
        assert!(blocks.len() >= 2, "A's right-bottom edge is split at the crossing");
        assert!(
            blocks.iter().all(|pb| pb.ext_paves().is_empty()),
            "all extra paves consumed into elementary blocks"
        );
    }

    #[test]
    fn split_pave_blocks_splits_block_with_two_extra_paves() {
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(3.0, 0.0, 0.0));
        let v1 = b.make_vertex(GpPnt::new(1.0, 0.0, 0.0), 1e-7);
        let v2 = b.make_vertex(GpPnt::new(2.0, 0.0, 0.0), 1e-7);
        let mut ds = BopdsDS::new();
        let n_e = ds.append(e.0.clone()).unwrap();
        let n1 = ds.append(v1.into()).unwrap();
        let n2 = ds.append(v2.into()).unwrap();
        ds.init_pave_blocks_for_edge(n_e);
        // Add two extra paves inside the block.
        {
            let pbs = ds.change_pave_blocks_mut(n_e);
            pbs[0].append_ext_pave(BopdsPave::new(n1, 1.0));
            pbs[0].append_ext_pave(BopdsPave::new(n2, 2.0));
        }
        assert_eq!(ds.pave_blocks(n_e).len(), 1);
        assert!(ds.pave_blocks(n_e)[0].is_to_update());

        let mut ctx = FillCtx { fuzzy: 1e-7, glue_full: false, non_destructive: false, errors: Vec::new() };
        split_pave_blocks_impl(&mut ds, &mut ctx, &[n_e]).unwrap();
        let blocks = ds.pave_blocks(n_e);
        assert_eq!(blocks.len(), 3, "2 extra paves split 1 block into 3");
        let ranges: Vec<(f64, f64)> = blocks.iter().map(|pb| pb.range()).collect();
        assert_eq!(ranges[0], (0.0, 1.0));
        assert_eq!(ranges[1], (1.0, 2.0));
        assert_eq!(ranges[2], (2.0, 3.0));
    }

    #[test]
    fn perform_vf_classifies_inside_on_and_out() {
        let box_solid = unit_box();
        let b = TopoBuilder::new();
        // Vertices exactly on the bottom-face plane (z = 0) so the bounding-box
        // candidate filter (which has no tolerance gap) still pairs them.
        let inside = b.make_vertex(GpPnt::new(0.5, 0.5, 0.0), 1e-7);
        let on_edge = b.make_vertex(GpPnt::new(0.5, 0.0, 0.0), 1e-7);
        let mid_box = b.make_vertex(GpPnt::new(0.5, 0.5, 0.5), 1e-7);
        // Package the three vertices as a second argument (compound).
        let comp = b.make_compound_of(&[inside.0.clone(), on_edge.0.clone(), mid_box.0.clone()]);

        let mut f = PaveFiller::new();
        f.set_arguments(&[box_solid.solid.0.clone(), comp.0.clone()]);
        f.init().unwrap();
        perform_vf(&mut f).unwrap();

        // The bottom face of the box is face 0 (z = 0).
        let n_face = f.ds().index(&box_solid.faces[0].0).unwrap();
        let pool = f.ds().face_info_pool();
        let info = pool.iter().find(|fi| fi.face_index == n_face);
        assert!(info.is_some(), "face info initialized");
        let info = info.unwrap();
        let n_inside = f.ds().index(&inside.0).unwrap();
        let n_on = f.ds().index(&on_edge.0).unwrap();
        let n_mid = f.ds().index(&mid_box.0).unwrap();
        // On-face vertices are kept in the `verts` list, never in the section
        // `paves` (which must hold only edge indices).
        assert!(info.paves().is_empty(), "section paves hold no vertices");
        let recorded: Vec<usize> = info.verts().iter().map(|p| p.0).collect();
        assert!(recorded.contains(&n_inside), "inside vertex recorded");
        assert!(recorded.contains(&n_on), "on-boundary vertex recorded");
        assert!(!recorded.contains(&n_mid), "floating vertex not recorded");
        // The interference table records the inside/on pairs.
        assert!(f.ds().has_interf_pair(n_inside, n_face));
        assert!(f.ds().has_interf_pair(n_on, n_face));
    }

    #[test]
    fn perform_vf_glue_full_only_initializes_face_info() {
        let box_solid = unit_box();
        let b = TopoBuilder::new();
        let v = b.make_vertex(GpPnt::new(0.5, 0.5, 0.0), 1e-7);
        let comp = b.make_compound_of(&[v.0.clone()]);
        let mut f = PaveFiller::new();
        f.set_glue(GlueEnum::Full);
        f.set_arguments(&[box_solid.solid.0.clone(), comp.0.clone()]);
        f.init().unwrap();
        perform_vf(&mut f).unwrap();
        let n_face = f.ds().index(&box_solid.faces[0].0).unwrap();
        let pool = f.ds().face_info_pool();
        assert!(
            pool.iter().any(|fi| fi.face_index == n_face),
            "face info initialized in glue-full mode"
        );
        let n_v = f.ds().index(&v.0).unwrap();
        assert!(!f.ds().has_interf_pair(n_v, n_face));
    }

    #[test]
    fn treat_new_vertices_fuses_coincident_seeds() {
        let mut ds = BopdsDS::new();
        let b = TopoBuilder::new();
        let e1 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 0.0, 0.0));
        let e2 = b.make_edge_segment(&GpPnt::new(0.0, 0.0, 0.0), &GpPnt::new(1.0, 1.0, 0.0));
        let ne1 = ds.append(e1.0.clone()).unwrap();
        let ne2 = ds.append(e2.0.clone()).unwrap();
        ds.init_pave_blocks_for_edge(ne1);
        ds.init_pave_blocks_for_edge(ne2);
        let seeds = vec![
            NewVertexSeed {
                point: GpPnt::new(0.5, 0.0, 0.0),
                tol: 1e-7,
                edge_a: ne1,
                t_a: 0.5,
                edge_b: ne2,
                t_b: 0.5,
                range1_first: 0.5,
                range1_last: 0.5,
            },
            NewVertexSeed {
                point: GpPnt::new(0.5, 0.0, 0.0),
                tol: 1e-7,
                edge_a: ne1,
                t_a: 0.5,
                edge_b: ne2,
                t_b: 0.5,
                range1_first: 0.5,
                range1_last: 0.5,
            },
        ];
        let idx = treat_new_vertices(&mut ds, 1e-7, &seeds).unwrap();
        assert_eq!(idx.len(), 1, "coincident seeds fuse into one vertex");
        assert_eq!(ds.pave_blocks(ne1)[0].ext_paves().len(), 1, "single extra pave on edge 1");
        assert_eq!(ds.pave_blocks(ne2)[0].ext_paves().len(), 1, "single extra pave on edge 2");
    }

    // -----------------------------------------------------------------------
    // perform_ef
    // -----------------------------------------------------------------------

    #[test]
    fn perform_ef_creates_vertex_at_piercing_point() {
        // A unit box + an edge piercing its bottom face (z = 0) at (0.5, 0.5, 0).
        let box_solid = unit_box();
        let b = TopoBuilder::new();
        let e = b.make_edge_segment(&GpPnt::new(0.5, 0.5, -1.0), &GpPnt::new(0.5, 0.5, 1.0));
        let comp = b.make_compound_of(&[e.0.clone()]);

        let mut f = PaveFiller::new();
        f.set_arguments(&[box_solid.solid.0.clone(), comp.0.clone()]);
        f.init().unwrap();
        init_all_edge_blocks(f.ds_mut());
        perform_ef(&mut f).unwrap();
        assert!(!f.has_errors(), "errors: {:?}", f.errors());

        // A new vertex close to (0.5, 0.5, 0) must exist among the new shapes.
        let n_src = f.ds().nb_source_shapes();
        let mut found = false;
        for i in n_src..f.ds().nb_shapes() {
            if f.ds().shape_info(i).map(|s| s.shape_type()) != Some(ShapeType::Vertex) {
                continue;
            }
            let p = BRepTool::vertex_point(&Vertex(f.ds().shape(i).unwrap().clone()));
            if p.distance(&GpPnt::new(0.5, 0.5, 0.0)) < 1e-4 {
                found = true;
                break;
            }
        }
        assert!(found, "piercing vertex created in the DS");

        // The E/F interference is recorded and the edge is split at the
        // piercing point (the extra pave was consumed into elementary blocks).
        let n_e = f.ds().index(&e.0).unwrap();
        let n_f = f.ds().index(&box_solid.faces[0].0).unwrap();
        assert!(f.ds().has_interf_pair(n_e, n_f));
        let blocks = f.ds().pave_blocks(n_e);
        assert!(
            blocks.iter().any(|pb| pb.pave1().parameter() == 1.0 || pb.pave2().parameter() == 1.0),
            "piercing point is a bound of an elementary block"
        );
        assert!(
            blocks.iter().all(|pb| pb.ext_paves().is_empty()),
            "all extra paves consumed into elementary blocks"
        );
    }

    #[test]
    fn perform_ef_ignores_face_boundary_edge() {
        // Edge 0 of the box IS a boundary of face 0 — no new vertex must be
        // created for the pair (the EF pass skips boundary sub-shapes).
        let box_solid = unit_box();
        let mut f = PaveFiller::new();
        f.set_arguments(&[box_solid.solid.0.clone()]);
        f.init().unwrap();
        let n_before = f.ds().nb_shapes();
        init_all_edge_blocks(f.ds_mut());
        perform_ef(&mut f).unwrap();
        assert!(!f.has_errors());
        // One argument → no inter-argument pairs; the DS is untouched.
        assert_eq!(f.ds().nb_shapes(), n_before);
    }

    // -----------------------------------------------------------------------
    // perform_ff
    // -----------------------------------------------------------------------

    #[test]
    fn perform_ff_check_planes_skips_unshared_planar_pair() {
        // Box A = [0,1]^3, box B shifted by (0.5, 0.5, 0). A's bottom (z = 0)
        // and B's front (y = 0.5) are planes that do not share On/In vertices
        // yet, so `CheckPlanes` skips FaceFace (`BOPAlgo_PaveFiller_6.cxx`).
        // The pair is recorded as an empty InterfFF; section geometry comes
        // from EF / later stages.
        let a = unit_box();
        let bbox = box_at(0.5, 0.5, 0.0);
        let mut f = PaveFiller::new();
        f.set_arguments(&[a.solid.0.clone(), bbox.solid.0.clone()]);
        f.init().unwrap();
        perform_ff(&mut f).unwrap();
        assert!(!f.has_errors(), "errors: {:?}", f.errors());

        let n_src = f.ds().nb_source_shapes();
        let mut found = false;
        for i in n_src..f.ds().nb_shapes() {
            if f.ds().shape_info(i).map(|s| s.shape_type()) != Some(ShapeType::Edge) {
                continue;
            }
            let e = Edge(f.ds().shape(i).unwrap().clone());
            let Some(curve) = BRepTool::edge_curve(&e) else { continue };
            let (t1, t2) = BRepTool::edge_parameters(&e);
            if !t1.is_finite() || !t2.is_finite() || t2 <= t1 {
                continue;
            }
            let mid = curve.d0(0.5 * (t1 + t2));
            if mid.distance(&GpPnt::new(0.75, 0.5, 0.0)) < 1e-4 {
                found = true;
                break;
            }
        }
        assert!(!found, "CheckPlanes must not emit a FaceFace section edge");

        let n_fa = f.ds().index(&a.faces[0].0).unwrap(); // A bottom (z = 0)
        let n_fb = f.ds().index(&bbox.faces[2].0).unwrap(); // B front (y = 0.5)
        // OCCT CheckPlanes skip: Append InterfFF Init(0,0), no AddInterf.
        assert!(
            f.ds().interf_ff().iter().any(|ff| {
                let (a, b) = ff.indices();
                (a == n_fa && b == n_fb) || (a == n_fb && b == n_fa)
            }),
            "CheckPlanes skip must still record an empty InterfFF"
        );
    }

    #[test]
    fn perform_ff_disjoint_faces_produce_no_edges() {
        // A single box: no inter-argument face pairs, so no section edges.
        let box_solid = unit_box();
        let mut f = PaveFiller::new();
        f.set_arguments(&[box_solid.solid.0.clone()]);
        f.init().unwrap();
        let n_before = f.ds().nb_shapes();
        perform_ff(&mut f).unwrap();
        assert!(!f.has_errors());
        assert_eq!(f.ds().nb_shapes(), n_before);
    }

    #[test]
    fn update_vertex_records_increased_tolerance() {
        // `update_vertex_ds` (port of `BOPAlgo_PaveFiller::UpdateVertex`) records
        // the vertex in `increased_ss` when its tolerance grows, so the
        // repeat-intersection stage can find it.
        let a = unit_box();
        let mut f = PaveFiller::new();
        f.set_arguments(&[a.solid.0.clone()]);
        f.init().unwrap();
        let n_v = f.ds().index(&a.vertices[0].0).unwrap(); // corner (0,0,0)
        let n_new = update_vertex_ds(f.ds_mut(), n_v, 0.5, false);
        assert_eq!(n_new, n_v, "non-destructive off → in-place update");
        assert!(
            f.ds().increased_ss().contains(&n_v),
            "the increased vertex must be recorded for the repeat intersection"
        );
        // The DS box grew to cover the enlarged tolerance sphere.
        let (x0, x1, _, _, _, _) = f.ds().box_of(n_v).unwrap().get().unwrap();
        assert!((x1 - x0 - 1.0).abs() < 1e-9, "corner + 0.5 tolerance → span 1.0, got {}", x1 - x0);
    }

    #[test]
    fn intersect_ext_pairs_finds_vertex_on_edge() {
        // A = [0,1]^3, B shifted by (0.5,0,0). B's corner (0.5,0,0) lies on A's
        // edge (0,0,0)-(1,0,0), so `intersect_ext_pairs` (port of
        // `BOPDS_Iterator::IntersectExt`) must report the V/E pair for that
        // vertex.
        let a = unit_box();
        let b = box_at(0.5, 0.0, 0.0);
        let mut f = PaveFiller::new();
        f.set_arguments(&[a.solid.0.clone(), b.solid.0.clone()]);
        f.init().unwrap();
        let nb_v = f.ds().index(&b.vertices[0].0).unwrap(); // (0.5,0,0)
        let na_e = f.ds().index(&a.edges[0].0).unwrap();
        let map: HashSet<usize> = [nb_v].into_iter().collect();
        let buckets = crate::bopds::intersect_ext_pairs(f.ds(), &map);
        let ve = &buckets[1]; // V/E bucket
        assert!(
            ve.iter().any(|&(v, e)| (v == nb_v || e == nb_v) && (v == na_e || e == na_e)),
            "expected the V/E pair of the on-edge vertex, got {ve:?}"
        );
    }
}

mod fill_ctx;
mod vertex_face;
mod face_face;
mod perform;
pub use fill_ctx::*;
pub use vertex_face::*;

pub use perform::*;
