//! `BRepFill_TrimShellCorner` (`BRepFill_TrimShellCorner.cxx`).
//!
//! The class trims the two corner shells of `BRepFill_Sweep::PerformCorner` so
//! that the two pipe parts end on the bisector plane. Entire class is ported:
//!
//! * constructor and `AddBounds` / `AddUEdges` / `AddVEdges` (`cxx:206-262`);
//! * `Perform` (`cxx:266-410`): the two shells built from the face columns,
//!   the vertex tolerance scan, `BOPAlgo_PaveFiller` over the pair with the
//!   `4 * Precision::Confusion()` fuzzy value, `SplitUEdges`, then the
//!   `InterfFF` walk dispatching to `MakeFacesNonSec` / `MakeFacesSec`;
//! * `IsDone` / `HasSection` / `Modified` (`cxx:414-436`);
//! * `MakeFacesNonSec` (`cxx:443-783`) and `MakeFacesSec` (`cxx:789-1018`);
//! * `ChooseSection` (`cxx:1020-1285`) and its geometry leaves (`cxx:2620-2855`);
//! * `SplitUEdges` (`cxx:1289-1392`);
//! * `FindFreeVertices` (`cxx:1409-1435`);
//! * `FindCommonVertex(BOPDS_PDS, ...)` (`cxx:1440-1508`);
//! * `GetUEdges` (`cxx:1511-1575`);
//! * the section walk `FillGap` (`cxx:1577-1613`), `FindNextEdge`
//!   (`cxx:1616-1703`), `CheckAndOrientEdges` (`cxx:1705-1827`),
//!   `FindVertex` (`cxx:1829-1887`), `FindNextVertex` (`cxx:1889-1931`),
//!   `GetPave` (`cxx:1933-1958`), `FindFromUEdge` (`cxx:1962-2070`) and
//!   `FindFromVEdge` (`cxx:2073-2513`);
//! * `RemoveEdges` (`cxx:2520-2547`) and `FilterSectionEdges` (`cxx:2549-2617`).
//!
//! UNPORTED:
//! * `IntTools_BeanFaceIntersector` arm of `FilterSectionEdges`
//!   (`cxx:2576-2606`) — dead in OCCT because `MakeFacesSec` passes an
//!   always-null `aSecPlane` (`cxx:800` -> `cxx:802`);
//! * `BOPDS_DS::Shape(-1)` in `FindVertex` (`cxx:1867-1881`): the Rust pave
//!   carries `usize::MAX` and the walk then reports failure instead of
//!   dereferencing a null shape.
//!
//! Translation notes:
//! - `NCollection_HArray2<TopoDS_Shape>` is [`ShapeTable`] (1-based);
//!   `NCollection_HArray1<TopoDS_Shape>` is a `Vec<Option<TopoShape>>` with an
//!   unused slot 0, the same convention `shell_tables.rs` documents.
//! - `myHistMap` is a `HashMap<ShapeMapKey, Vec<TopoShape>>`
//!   (`NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
//!   TopTools_ShapeMapHasher>`).
//! - `BRep_Builder::Add(E, V)` + `Range(E, f, l)` is `init_3d_vertices` +
//!   `TopoBuilder::make_edge` + `add_endpoint_vertices`, the same pair the rest
//!   of `brep_fill_sweep` uses to build a bounded edge on an existing curve.
use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx2, GpDir, GpLin, GpPln, GpPnt, GpPnt2d, GpVec};
use occt_core::precision::{CONFUSION, REAL_SMALL, SQUARE_CONFUSION};
use occt_geom::approx_same_parameter::u_resolution;
use occt_geom::curve::Curve;
use occt_geom::gcpnts::uniform_abscissa;
use occt_geom::geom_lib::geom_lib_axe_of_inertia;
use occt_geom::line::GeomLine;

use crate::abs::{Orientation, ShapeType};
use crate::bopds::BopdsDS;
use crate::bopds_pave::BopdsPave;
use crate::bopds_ff::BopdsCurve;
use crate::brep_lib_make_edge::{add_endpoint_vertices, init_3d_vertices};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::inttools_data::CommonPartType;
use crate::pave_filler::PaveFiller;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::shape_fix_compose_shell::{empty_copied_edge, MapReShape, ReShape};
use crate::shhealing::{first_vertex, last_vertex, ShapeFixEdge};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{
    common_vertex, edge_vertices, edges_of_wire, explorer_shapes, is_same,
    map_shapes_and_ancestors, shape_map_key, vertices_of, wire_vertices, ShapeMapKey,
};

use super::helpers::ShapeTable;
use super::TransitionStyle;

/// `BRepFill_TrimShellCorner` (`BRepFill_TrimShellCorner.hxx`).
pub struct TrimShellCorner {
    /// `myTransition`.
    pub(super) transition: TransitionStyle,
    /// `myAxeOfBisPlane`.
    pub(super) axe_of_bis_plane: GpAx2,
    /// `myIntPointCrossDir`.
    pub(super) int_point_cross_dir: GpVec,
    /// `myShape1` / `myShape2`.
    pub(super) shape1: Option<TopoShape>,
    pub(super) shape2: Option<TopoShape>,
    /// `myBounds` — the boundary edges of the sections, one column per part.
    pub(super) bounds: Option<ShapeTable>,
    /// `myUEdges`.
    pub(super) u_edges: Option<ShapeTable>,
    /// `myVEdges` — `NCollection_HArray1<TopoDS_Shape>`, 1-based (slot 0
    /// unused).
    pub(super) v_edges: Vec<Option<TopoShape>>,
    /// `myFaces`.
    pub(super) faces: Option<ShapeTable>,
    /// `myDone`.
    pub(super) done: bool,
    /// `myHasSection`.
    pub(super) has_section: bool,
    /// `myHistMap` — `NCollection_DataMap<TopoDS_Shape,
    /// NCollection_List<TopoDS_Shape>>` keyed like `ShapeMapKey`.
    pub(super) hist_map: HashMap<ShapeMapKey, Vec<TopoShape>>,
}

impl TrimShellCorner {
    /// Constructor (`cxx:206-222`): copies the face table, stores the
    /// transition, the bisector frame and the crossing direction of the two
    /// paths. `myDone` / `myHasSection` start `false`.
    pub fn new(
        the_faces: &ShapeTable,
        the_transition: TransitionStyle,
        the_axe_of_bis_plane: GpAx2,
        the_int_point_cross_dir: GpVec,
    ) -> Self {
        Self {
            transition: the_transition,
            axe_of_bis_plane: the_axe_of_bis_plane,
            int_point_cross_dir: the_int_point_cross_dir,
            shape1: None,
            shape2: None,
            bounds: None,
            u_edges: None,
            v_edges: Vec::new(),
            faces: Some(the_faces.clone()),
            done: false,
            has_section: false,
            hist_map: HashMap::new(),
        }
    }

    /// `AddBounds(Bounds)` (`cxx:228-236`): copies the boundary table.
    pub fn add_bounds(&mut self, the_bounds: &ShapeTable) {
        self.bounds = Some(the_bounds.clone());
    }

    /// `AddUEdges(UEdges)` (`cxx:240-250`): copies the `u`-iso edge table.
    pub fn add_u_edges(&mut self, the_u_edges: &ShapeTable) {
        self.u_edges = Some(the_u_edges.clone());
    }

    /// `AddVEdges(theVEdges, theIndex)` (`cxx:252-262`): takes column
    /// `the_index` of the `v`-iso table into the 1-based `myVEdges` array.
    pub fn add_v_edges(&mut self, the_v_edges: &ShapeTable, the_index: usize) {
        let up = the_v_edges.nb_rows();
        self.v_edges = vec![None; up + 1];
        for i in 1..=up {
            self.v_edges[i] = the_v_edges.value(i, the_index).cloned();
        }
    }

    /// `IsDone()` (`cxx:414-417`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `HasSection()` (`cxx:421-425`).
    pub fn has_section(&self) -> bool {
        self.has_section
    }

    /// `Modified(theShape, theModified)` (`cxx:428-436`): the recorded images
    /// of `theShape`, empty when the map carries no entry.
    pub fn modified(&self, the_shape: &TopoShape) -> Vec<TopoShape> {
        self.hist_map
            .get(&shape_map_key(the_shape))
            .cloned()
            .unwrap_or_default()
    }

    /// `Perform()` (`cxx:266-410`).
    pub fn perform(&mut self) -> Result<(), String> {
        self.done = false;
        self.hist_map.clear();

        let Some(faces) = self.faces.clone() else {
            return Ok(());
        };
        // `cxx:274-277`: only the two-part layout is handled.
        if faces.nb_rows() != 2 {
            return Ok(());
        }
        if faces.nb_cols() < 2 {
            return Ok(());
        }

        // `cxx:281-301`: one shell per part, `Closed(BRep_Tool::IsClosed)`.
        let b = TopoBuilder::new();
        let mut shells: Vec<TopoShape> = Vec::new();
        for jj in 1..=faces.nb_cols() {
            let mut shell_faces: Vec<crate::shape::Face> = Vec::new();
            for ii in 1..=faces.nb_rows() {
                if let Some(f) = faces.value(ii, jj) {
                    shell_faces.push(crate::shape::Face(f.clone()));
                }
            }
            let shell = b.make_shell(&shell_faces);
            shell.set_closed(BRepTool::is_closed(&shell.0));
            shells.push(shell.0);
        }
        let Some(shape1) = shells.first().cloned() else {
            return Ok(());
        };
        let Some(shape2) = shells.get(1).cloned() else {
            return Ok(());
        };
        self.shape1 = Some(shape1.clone());
        self.shape2 = Some(shape2.clone());

        // `cxx:304-314`: the maximum vertex tolerance over both shells.
        let mut max_tol = 0.0f64;
        for s in [&shape1, &shape2] {
            for v in vertices_of(s) {
                max_tol = max_tol.max(BRepTool::vertex_tolerance(&v));
            }
        }

        // `cxx:316-331`: pave the pair; the fuzzy value only rises to
        // `4 * Precision::Confusion()` when both shells are already tight.
        let mut a_fuzzy = 4.0 * CONFUSION;
        let mut pf = PaveFiller::new();
        pf.set_arguments(&[shape1.clone(), shape2.clone()]);
        if max_tol < 1.005 * CONFUSION {
            a_fuzzy = a_fuzzy.max(pf.fuzzy_value());
            pf.set_fuzzy_value(a_fuzzy);
        }
        pf.perform()?;
        if pf.has_errors() {
            return Ok(());
        }

        // `cxx:337-340`: split the `u`-iso edges against the interferences.
        {
            let Some(u_edges) = self.u_edges.clone() else {
                return Ok(());
            };
            if !split_u_edges(&u_edges, pf.ds(), &self.int_point_cross_dir, &mut self.hist_map) {
                return Ok(());
            }
        }

        // `cxx:342-408`: walk the face/face interferences of the pair.
        let nb_ffs = pf.ds().interf_ff().len();
        for ii in 1..=faces.nb_rows() {
            let Some(f1) = faces.value(ii, 1).cloned() else {
                continue;
            };
            let Some(f2) = faces.value(ii, faces.nb_cols()).cloned() else {
                continue;
            };
            let (Some(an_index1), Some(an_index2)) = (pf.ds().index(&f1), pf.ds().index(&f2)) else {
                continue;
            };

            for i in 0..nb_ffs {
                let (n_f1, n_f2) = pf.ds().interf_ff()[i].indices();
                let a_nb_p = pf.ds().interf_ff()[i].points().len();
                let a_nb_c = pf.ds().interf_ff()[i].curves().len();
                if a_nb_p == 0 && a_nb_c == 0 && !pf.ds().has_interf_sub_shapes(n_f1, n_f2) {
                    continue;
                }
                if n_f1 != an_index1 || n_f2 != an_index2 {
                    continue;
                }
                // `cxx:376-388`: a section exists when some intersection curve
                // carries pave blocks.
                let mut has_sec = false;
                for j in 0..a_nb_c {
                    if !pf.ds().interf_ff()[i].curves()[j]
                        .pave_blocks()
                        .is_empty()
                    {
                        has_sec = true;
                        break;
                    }
                }
                if !has_sec {
                    if !self.make_faces_non_sec(ii, pf.ds(), an_index1, an_index2)? {
                        self.hist_map.clear();
                        return Ok(());
                    }
                } else if !self.make_faces_sec(ii, pf.ds(), an_index1, an_index2, i)? {
                    self.hist_map.clear();
                    return Ok(());
                }
                break;
            }
        }
        self.done = true;
        Ok(())
    }

    /// `MakeFacesNonSec` (`cxx:443-783`).
    ///
    /// The bounds of the corner row do not share a section curve: the two face
    /// columns are joined at the common vertex found on their bound edges
    /// (`FindCommonVertex`) and/or between their `u`-iso edge pairs, then the
    /// substituted faces are recorded in `myHistMap`.
    fn make_faces_non_sec(
        &mut self,
        the_index: usize,
        the_ds: &BopdsDS,
        the_face_index1: usize,
        the_face_index2: usize,
    ) -> Result<bool, String> {
        // `cxx:448-452`.
        let mut has_new_edge = false;
        let mut new_edge: Option<Edge> = None;

        let Some(bounds) = self.bounds.clone() else {
            return Ok(false);
        };
        let Some(u_edges) = self.u_edges.clone() else {
            return Ok(false);
        };
        let (Some(a_e1), Some(a_e2)) = (
            bounds.value(the_index, 1).cloned(),
            bounds.value(the_index, 2).cloned(),
        ) else {
            return Ok(false);
        };

        // `cxx:455-464`: the common vertex between the two bounds.
        let an_index1 = the_ds.index(&a_e1).unwrap_or(usize::MAX);
        let an_index2 = the_ds.index(&a_e2).unwrap_or(usize::MAX);
        let mut common_vertex: Option<TopoShape> = None;
        let mut apar1 = 0.0f64;
        let mut apar2 = 0.0f64;
        let vertex_found = find_common_vertex(
            the_ds,
            an_index1,
            an_index2,
            &self.int_point_cross_dir,
            &mut common_vertex,
            &mut apar1,
            &mut apar2,
        );

        let mut substitutor = MapReShape::new();

        // `cxx:467-517`: the common vertices between the `u`-iso edge pairs.
        let mut common_vertices: Vec<TopoShape> = Vec::new();
        let mut common_flag = 0i32;
        let low_col = 1usize;
        let up_col = u_edges.nb_cols();
        for (ueit, eindex) in [(1usize, the_index), (2usize, the_index + 1)] {
            let (Some(a_shape1), Some(a_shape2)) = (
                u_edges.value(eindex, low_col).cloned(),
                u_edges.value(eindex, up_col).cloned(),
            ) else {
                continue;
            };
            // `cxx:479-497`: an already-substituted `u`-iso edge wins.
            let a_ue1 = hist_first(&self.hist_map, &a_shape1).unwrap_or_else(|| a_shape1.clone());
            let a_ue2 = hist_first(&self.hist_map, &a_shape2).unwrap_or_else(|| a_shape2.clone());
            // `cxx:499-506`.
            if !is_same(&a_shape1, &a_ue1) {
                substitutor.replace(&forward(&a_shape1), &forward(&a_ue1));
            }
            if !is_same(&a_shape2, &a_ue2) {
                substitutor.replace(&forward(&a_shape2), &forward(&a_ue2));
            }
            // `cxx:509-517`: `LastVertex(aUE1)` / `FirstVertex(aUE2)`.
            let v1 = last_vertex(&Edge(a_ue1));
            let v2 = first_vertex(&Edge(a_ue2));
            if let (Some(v1), Some(v2)) = (v1, v2) {
                if is_same(&v1.0, &v2.0) {
                    common_flag = if common_flag == 0 { ueit as i32 } else { 3 };
                    common_vertices.push(v1.0);
                }
            }
        }

        // `cxx:520-537`.
        if vertex_found {
            if common_vertices.len() != 1 {
                return Ok(false);
            }
            let Some(cv) = common_vertex.clone() else {
                return Ok(false);
            };
            let cv = Vertex(cv);
            let first_cv = Vertex(common_vertices[0].clone());
            new_edge = if common_flag == 1 {
                make_edge_vv(&first_cv, &cv)
            } else {
                make_edge_vv(&cv, &first_cv)
            };
            has_new_edge = true;
        }

        // `cxx:539-546`.
        if common_vertices.len() == 2 {
            new_edge = make_edge_vv(
                &Vertex(common_vertices[0].clone()),
                &Vertex(common_vertices[1].clone()),
            );
            has_new_edge = true;
        }

        let builder = TopoBuilder::new();
        for fit in 1..=2usize {
            // `cxx:548-577`: assemble the compound whose free vertices are the
            // unmatched ends of the two `u`-iso edges.
            let mut comp_children: Vec<TopoShape> = Vec::new();
            let mut map_v: Vec<TopoShape> = Vec::new();
            for eindex in [the_index, the_index + 1] {
                let Some(a_shape) = u_edges.value(eindex, low_col + fit - 1).cloned() else {
                    continue;
                };
                let a_ue =
                    hist_first(&self.hist_map, &a_shape).unwrap_or_else(|| a_shape.clone());
                let a_v = if fit == 1 {
                    first_vertex(&Edge(a_ue.clone()))
                } else {
                    last_vertex(&Edge(a_ue.clone()))
                };
                if let Some(v) = a_v {
                    map_v.push(v.0);
                }
                comp_children.push(a_ue);
            }
            if has_new_edge {
                if let Some(ne) = new_edge.clone() {
                    comp_children.push(ne.0.clone());
                    store_vedge_in_hist_map(&self.v_edges, the_index, &ne.0, &mut self.hist_map);
                }
            }
            let comp = builder.make_compound_of(&comp_children);

            // `cxx:579-584`.
            let mut alone_vertices: Vec<TopoShape> = Vec::new();
            find_free_vertices(&comp.0, &map_v, &mut alone_vertices);
            if !alone_vertices.is_empty() && alone_vertices.len() != 2 {
                return Ok(false);
            }

            // `cxx:586-591`.
            let a_face_index = if fit == 1 {
                the_face_index1
            } else {
                the_face_index2
            };
            let Some(mut a_face) = the_ds.shape(a_face_index).cloned() else {
                return Ok(false);
            };
            let a_face_ori = a_face.orientation();
            a_face.set_orientation(Orientation::Forward);

            // `cxx:593-618`.
            if has_new_edge {
                if let Some(mut ne) = new_edge.clone() {
                    ne.0.set_orientation(Orientation::Forward);
                    // `cxx:596-612`: reverse the new edge when both of its ends
                    // sit closer to the *other* end of the reference bound edge.
                    let a_e_ref = Edge(if fit == 1 {
                        a_e1.clone()
                    } else {
                        a_e2.clone()
                    });
                    let (v1, v2) = (
                        first_vertex(&ne),
                        last_vertex(&ne),
                    );
                    let (r1, r2) = (
                        first_vertex(&a_e_ref),
                        last_vertex(&a_e_ref),
                    );
                    if let (Some(v1), Some(v2), Some(r1), Some(r2)) = (v1, v2, r1, r2) {
                        let p1 = BRepTool::vertex_point(&v1);
                        let p2 = BRepTool::vertex_point(&v2);
                        let q1 = BRepTool::vertex_point(&r1);
                        let q2 = BRepTool::vertex_point(&r2);
                        let p11 = p1.distance(&q1);
                        let p22 = p2.distance(&q2);
                        let p12 = p1.distance(&q2);
                        let p21 = p2.distance(&q1);
                        if p11 > p12 && p22 > p21 {
                            ne.0.reverse();
                        }
                    }
                    // `cxx:614-618`: add the pcurve on a non-planar face.
                    ShapeFixEdge::new().fix_add_pcurve(&ne, &a_face);
                    new_edge = Some(ne);
                }
            }

            // `cxx:620-726`.
            let mut ordered: Vec<Edge> = Vec::new();
            if alone_vertices.is_empty() {
                if let Some(ne) = new_edge.clone() {
                    ordered.push(ne);
                }
            } else {
                let an_e_index = if fit == 1 { an_index1 } else { an_index2 };
                let mut found1 = false;
                let mut found2 = false;
                let mut param1 = 0.0f64;
                let mut param2 = 0.0f64;
                // `cxx:629-657`: the two paves of the bound edge that carry the
                // free vertices.
                for pave in the_ds.paves(an_e_index) {
                    let Some(a_v) = the_ds.shape(pave.index()).cloned() else {
                        continue;
                    };
                    if is_same(&a_v, &alone_vertices[0]) && !found1 {
                        param1 = pave.parameter();
                        found1 = true;
                    }
                    if is_same(&a_v, alone_vertices.last().unwrap()) && !found2 {
                        param2 = pave.parameter();
                        found2 = true;
                    }
                }
                if !(found1 && found2) {
                    return Ok(false); // `cxx:716-719`.
                }
                // `cxx:659-688`: re-cut the bound edge between the two paves.
                let src = if fit == 1 { &a_e1 } else { &a_e2 };
                let (v_first, v_last, lo, hi) = if param1 < param2 {
                    (
                        alone_vertices[0].clone(),
                        alone_vertices[alone_vertices.len() - 1].clone(),
                        param1,
                        param2,
                    )
                } else {
                    (
                        alone_vertices[alone_vertices.len() - 1].clone(),
                        alone_vertices[0].clone(),
                        param2,
                        param1,
                    )
                };
                let mut v1 = Vertex(v_first);
                let mut v2 = Vertex(v_last);
                v1.0.set_orientation(Orientation::Forward);
                v2.0.set_orientation(Orientation::Reversed);
                let mut new_bound = empty_copied_edge(src);
                new_bound.set_orientation(Orientation::Forward);
                builder.add(&mut new_bound, &v1.0);
                builder.add(&mut new_bound, &v2.0);
                GeometryRegistry::global().set_ranges_all(&new_bound, lo, hi);
                ordered.push(Edge(new_bound));
                // `cxx:692-714`: prepend the new edge when it starts at `aV2`.
                if let Some(ne) = new_edge.clone() {
                    let starts_at_v2 = explorer_shapes(
                        &ne.0,
                        ShapeType::Vertex,
                    )
                    .iter()
                    .any(|v| is_same(v, &v2.0));
                    if starts_at_v2 {
                        ordered.push(ne);
                    } else {
                        ordered.insert(0, ne);
                    }
                }
            }

            // `cxx:729-749`.
            if !ordered.is_empty() {
                let mut wire = Wire::new();
                for e in &ordered {
                    builder.add(&mut wire.0, &e.0);
                }
                let target = if fit == 1 { &a_e1 } else { &a_e2 };
                substitutor.replace(&forward(target), &wire.0);
            }

            // `cxx:749-752`.
            substitutor.apply(&a_face);
            let new_face = substitutor.value(&a_face).unwrap_or_else(|| a_face.clone());
            let mut new_face = new_face;
            new_face.set_orientation(a_face_ori);
            self.hist_map
                .insert(shape_map_key(&a_face), vec![new_face]);

            // `cxx:753-779`: every edge of the source face that the
            // substitution rewrote is re-recorded as its list of new edges.
            for e in explorer_shapes(&a_face, ShapeType::Edge) {
                let Some(a_new_value) = substitutor.value(&e) else {
                    continue;
                };
                if is_same(&a_new_value, &e) {
                    continue;
                }
                if self.hist_map.contains_key(&shape_map_key(&e)) {
                    continue;
                }
                let list =
                    explorer_shapes(&a_new_value, ShapeType::Edge);
                self.hist_map.insert(shape_map_key(&e), list);
            }
        }

        Ok(true) // `cxx:782`.
    }

    /// `MakeFacesSec` (`cxx:789-1019`).
    /// `MakeFacesSec` (`cxx:789-1018`).
    ///
    /// The two face columns share a real section curve, so the bound edges are
    /// replaced by the section-edges chain that runs from the common vertex of
    /// the left `u`-iso pair to the common vertex of the right pair
    /// (`FindFromUEdge` then `FindFromVEdge`).
    fn make_faces_sec(
        &mut self,
        the_index: usize,
        the_ds: &BopdsDS,
        the_face_index1: usize,
        the_face_index2: usize,
        the_ss_interf_index: usize,
    ) -> Result<bool, String> {
        // `cxx:792-797`.
        let a_b_curves = match the_ds.interf_ff().get(the_ss_interf_index) {
            Some(ff) => ff.curves().to_vec(),
            None => return Ok(false),
        };

        // `cxx:799-805`: `FilterSectionEdges(aBCurves, aSecPlane, theDS, aSecEdges)`.
        // `aSecPlane` is declared but never assigned, so it stays a null face
        // and the reference's `IntTools_BeanFaceIntersector` arm is dead.
        let mut a_sec_edges = filter_section_edges(&a_b_curves, the_ds);

        // `cxx:806-808`.
        let Some(u_edges) = self.u_edges.clone() else {
            return Ok(false);
        };
        let Some(bounds) = self.bounds.clone() else {
            return Ok(false);
        };

        // `cxx:808-813`.
        let (Some(left_e1), Some(left_e2), Some(right_e1), Some(right_e2)) = (
            u_edges.value(the_index, 1).cloned(),
            u_edges.value(the_index, 2).cloned(),
            u_edges.value(the_index + 1, 1).cloned(),
            u_edges.value(the_index + 1, 2).cloned(),
        ) else {
            return Ok(false);
        };

        // `cxx:815-818`.
        let (Some(index_left_e1), Some(index_left_e2), Some(index_right_e1), Some(index_right_e2)) = (
            the_ds.index(&left_e1),
            the_ds.index(&left_e2),
            the_ds.index(&right_e1),
            the_ds.index(&right_e2),
        ) else {
            return Ok(false);
        };

        // `cxx:819-838`: the two common vertices of the `u`-iso pairs.
        let mut first_vertex: Option<TopoShape> = None;
        let mut last_vertex: Option<TopoShape> = None;
        let mut _param_left_e1 = 0.0f64;
        let mut _param_left_e2 = 0.0f64;
        let mut _param_right_e1 = 0.0f64;
        let mut _param_right_e2 = 0.0f64;
        find_common_vertex(
            the_ds,
            index_left_e1,
            index_left_e2,
            &self.int_point_cross_dir,
            &mut first_vertex,
            &mut _param_left_e1,
            &mut _param_left_e2,
        );
        find_common_vertex(
            the_ds,
            index_right_e1,
            index_right_e2,
            &self.int_point_cross_dir,
            &mut last_vertex,
            &mut _param_right_e1,
            &mut _param_right_e2,
        );

        // `cxx:840-857`: `ChooseSection`; a found wire replaces the section
        // compound by its own edges and is recorded as the `v`-iso image.
        if let Some((sec_wire, _sec_plane, _is_singular)) = choose_section(
            self.transition,
            &self.axe_of_bis_plane,
            &a_sec_edges,
            first_vertex.as_ref(),
            last_vertex.as_ref(),
        ) {
            a_sec_edges = explorer_shapes(&sec_wire, ShapeType::Edge);
            store_vedge_in_hist_map(&self.v_edges, the_index, &sec_wire, &mut self.hist_map);
        }

        // `cxx:860`.
        let mut substitutor = MapReShape::new();

        for fit in 0..2usize {
            // `cxx:862-869`.
            let a_face_index = if fit == 0 {
                the_face_index1
            } else {
                the_face_index2
            };
            let Some(a_face) = the_ds.shape(a_face_index).cloned() else {
                return Ok(false);
            };
            let a_face_ori = a_face.orientation();
            let mut a_face_f = a_face.clone();
            a_face_f.set_orientation(Orientation::Forward);
            let Some(a_bound_edge) = bounds.value(the_index, 1 + fit).cloned() else {
                return Ok(false);
            };
            // `cxx:870`.
            let Some(a_bound_edge_index) = the_ds.index(&a_bound_edge) else {
                return Ok(false);
            };

            // `cxx:871-878`.
            let Some((a_ue1, a_ue2)) = get_u_edges(
                the_index,
                fit,
                &u_edges,
                &a_bound_edge,
                &a_face_f,
            ) else {
                return Ok(false);
            };
            let a_ue1_old = a_ue1.clone();
            let a_ue2_old = a_ue2.clone();

            // `cxx:882-911`: substitute each `u`-iso edge through the map,
            // matching the reference's `lst.First().Oriented(edge.Orientation())`.
            let mut a_ue1 = a_ue1;
            let mut a_ue2 = a_ue2;
            for (slot, old) in [(&mut a_ue1, &a_ue1_old), (&mut a_ue2, &a_ue2_old)] {
                if let Some(an_edge) = hist_first(&self.hist_map, &old.0) {
                    let an_edge = oriented_to(&an_edge, old.0.orientation());
                    if !is_same(&old.0, &an_edge) {
                        substitutor.replace(&forward(&old.0), &forward(&an_edge));
                    }
                    *slot = Edge(an_edge);
                }
            }

            // `cxx:912-919`.
            let mut comp_of_sec_edges = a_sec_edges.clone();
            let mut list_of_wire_edges: Vec<TopoShape> = Vec::new();

            // `cxx:920-958`.
            let Some(u_result) = find_from_u_edge(
                &a_ue1_old.0,
                &a_ue2_old.0,
                &a_ue1.0,
                &a_ue2.0,
                &a_face_f,
                &a_sec_edges,
                fit,
                &a_bound_edge,
                a_bound_edge_index,
                the_ds,
                &self.hist_map,
                comp_of_sec_edges,
                list_of_wire_edges,
            ) else {
                return Ok(false);
            };
            comp_of_sec_edges = u_result.comp_of_sec_edges;
            list_of_wire_edges = u_result.list_of_wire_edges;

            // `cxx:933-957`.
            let Some(v_result) = find_from_v_edge(
                &u_result.found_pave,
                u_result.is_on_u_edge,
                &a_ue1_old.0,
                &a_ue2_old.0,
                &a_face_f,
                &comp_of_sec_edges,
                fit,
                &a_bound_edge,
                a_bound_edge_index,
                the_ds,
                &self.hist_map,
            ) else {
                return Ok(false);
            };
            if !v_result.is_section_found && list_of_wire_edges.is_empty() {
                return Ok(false);
            }
            list_of_wire_edges.extend(v_result.list_of_wire_edges);

            // `cxx:965-981`: a raw `BRep_Builder` wire (`MakeWire` + `Add`), so
            // the edges keep the order and orientation of the list.
            if !list_of_wire_edges.is_empty() {
                let mut edges: Vec<Edge> = Vec::new();
                for s in &list_of_wire_edges {
                    if !is_same(&a_bound_edge, s) {
                        edges.push(Edge(s.clone()));
                    }
                }
                let a_w = TopoBuilder::new().make_wire(&edges);
                substitutor.replace(&forward(&a_bound_edge), &a_w.0);
            }

            // `cxx:984-988`.
            substitutor.apply(&a_face);
            let Some(mut new_face) = substitutor.value(&a_face) else {
                return Ok(false);
            };
            new_face.set_orientation(a_face_ori);
            self.hist_map
                .insert(shape_map_key(&a_face), vec![new_face]);

            // `cxx:990-1015`.
            for e in explorer_shapes(&a_face, ShapeType::Edge) {
                let Some(a_new_value) = substitutor.value(&e) else {
                    continue; // `cxx:996-999`.
                };
                if is_same(&a_new_value, &e) {
                    continue;
                }
                if self.hist_map.contains_key(&shape_map_key(&e)) {
                    continue; // `cxx:1001-1004`.
                }
                let list = explorer_shapes(&a_new_value, ShapeType::Edge);
                self.hist_map.insert(shape_map_key(&e), list);
            }
        }

        Ok(true) // `cxx:1018`.
    }
}

// ---------------------------------------------------------------------------
// Statics
// ---------------------------------------------------------------------------

/// `myHistMap.Find(shape)` restricted to `lst.IsEmpty()` being false
/// (`cxx:479-497`, `cxx:556-566`): the first recorded image of `shape`.
fn hist_first(
    the_hist_map: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    the_shape: &TopoShape,
) -> Option<TopoShape> {
    the_hist_map
        .get(&shape_map_key(the_shape))
        .and_then(|list| list.first().cloned())
}

/// `shape.Oriented(TopAbs_FORWARD)` (`TopoDS_Shape::Oriented`): a copy put
/// FORWARD, leaving the argument untouched.
fn forward(the_shape: &TopoShape) -> TopoShape {
    let mut s = the_shape.clone();
    s.set_orientation(Orientation::Forward);
    s
}

/// `BRepLib_MakeEdge(V1, V2)` (`BRepLib_MakeEdge.cxx:185-198`): a straight
/// segment between two existing vertices. `None` is the reference's
/// `BRepLib_LineThroughIdenticPoints` (`cxx:191-195`); `cxx:196` then runs
/// `Init(GL, V1, V2, 0, l)`.
fn make_edge_vv(v1: &Vertex, v2: &Vertex) -> Option<Edge> {
    let p1 = BRepTool::vertex_point(v1);
    let p2 = BRepTool::vertex_point(v2);
    let l = p1.distance(&p2);
    if l <= REAL_SMALL {
        return None;
    }
    let dir = GpDir::from_vec(&GpVec::from_pnts(&p1, &p2)).ok()?;
    let line: Arc<dyn Curve> = Arc::new(GeomLine::new(GpLin::from_pnt_dir(p1, dir)));
    let build = init_3d_vertices(&line, Some(v1), Some(v2), 0.0, l)?;
    let b = TopoBuilder::new();
    let mut e = b.make_edge(line, build.first, build.last);
    add_endpoint_vertices(&b, &mut e, &build);
    Some(e)
}

/// `TopoDS_Iterator(shape)`'s components: `TopoDS_Iterator.cxx:32-46` walks
/// the shape's own child list, so the port reads it directly. The components
/// keep the orientation `TopoDS_Builder::Add` stored.
fn shape_children(the_shape: &TopoShape) -> Vec<TopoShape> {
    let ts = the_shape.tshape.read().expect("poisoned TShape lock");
    ts.children.clone()
}

/// `BRep_Builder::Remove(theComp, theEdge)` on a compound held as an edge
/// list (`TopoDS_Builder::Remove`, `TopoDS_Builder.cxx:107-133`): every
/// component that `IsSame`s `the_edge` is dropped.
fn comp_remove(the_comp: &mut Vec<TopoShape>, the_edge: &TopoShape) {
    the_comp.retain(|e| !is_same(e, the_edge));
}

/// `ComputeAveragePlaneAndMaxDeviation` (`cxx:2620-2658`).
///
/// `IsSingular` reports a degenerate point set (all points on a line), in
/// which case the reference returns `-1` (`cxx:2645-2648`) and leaves
/// `thePlane` at whatever `GeomLib::AxeOfInertia` produced.
fn compute_average_plane_and_max_deviation(
    the_wire: &TopoShape,
    the_plane: &mut GpPln,
    is_singular: &mut bool,
) -> f64 {
    // `cxx:2623-2624`.
    const N: usize = 40;
    *is_singular = false;
    let reg = GeometryRegistry::global();
    let mut pnts: Vec<GpPnt> = Vec::new();
    // `cxx:2628-2639`: `GCPnts_UniformAbscissa Distribution(aCurve, N + 1)`,
    // then `for (i = 1; i <= N; i++)`.
    for e in shape_children(the_wire) {
        let Some(curve) = reg.edge_curve(&e) else {
            continue;
        };
        let Ok(params) = uniform_abscissa(curve.as_ref(), N + 1) else {
            continue;
        };
        for i in 1..=N {
            if let Some(&par) = params.get(i) {
                pnts.push(curve.d0(par));
            }
        }
    }

    // `cxx:2641-2648`: `GeomLib::AxeOfInertia(Pnts, Axe, IsSingular)` with the
    // declaration's default `Tol = 1.0e-7` (`GeomLib.hxx`).
    let (axe, singular) = geom_lib_axe_of_inertia(&pnts, 1.0e-7);
    *is_singular = singular;
    if singular {
        return -1.0;
    }

    // `cxx:2649`.
    *the_plane = GpPln::new(axe.to_ax3());
    // `cxx:2650-2657`.
    let mut max_deviation = 0.0f64;
    for p in &pnts {
        let dist = the_plane.signed_distance_pnt(p).abs();
        if dist > max_deviation {
            max_deviation = dist;
        }
    }
    max_deviation
}

/// `UpdateSectionEdge` (`cxx:2661-2705`): replace one endpoint of
/// `the_edge` by `the_vertex` and shorten the edge range to `the_param`.
///
/// `F_Edge = theEdge` shares the TShape, so every write lands on `the_edge`
/// exactly as in the reference; `the_vertex`'s orientation is set in place
/// (`cxx:2700`).
fn update_section_edge(
    the_edge: &TopoShape,
    the_const_vertex: &TopoShape,
    the_vertex: &mut TopoShape,
    the_param: f64,
) {
    let reg = GeometryRegistry::global();
    // `cxx:2666-2667`: `F_Edge = theEdge; F_Edge.Orientation(TopAbs_FORWARD)`.
    let f_edge = forward(the_edge);
    let (v1, v2) = edge_vertices(&Edge(f_edge.clone()));
    // `cxx:2671-2684`: the vertex opposite `theConstVertex`.
    let (or_of_vertex, another_vertex) = match (&v1, &v2) {
        (Some(v1), Some(v2)) if is_same(&v1.0, the_const_vertex) => {
            (Orientation::Reversed, Some(v2.0.clone()))
        }
        (Some(v1), Some(_)) => (Orientation::Forward, Some(v1.0.clone())),
        _ => (Orientation::Forward, None),
    };

    // `cxx:2687-2698`: `BRep_Tool::Range` + `B.Begin`/`Range`.
    let (mut fpar, mut lpar) = reg.edge_parameters(&f_edge);
    if or_of_vertex == Orientation::Forward {
        fpar = the_param;
    } else {
        lpar = the_param;
    }
    reg.set_ranges_all(&f_edge, fpar, lpar);

    // `cxx:2699-2703`: `F_Edge.Free(true)`, `B.Remove(F_Edge, AnotherVertex)`,
    // `theVertex.Orientation(OrOfVertex)`, `B.Add(F_Edge, theVertex)`.
    f_edge.set_free(true);
    let builder = TopoBuilder::new();
    let mut target = f_edge.clone();
    if let Some(another) = another_vertex {
        builder.remove(&mut target, &another);
    }
    the_vertex.set_orientation(or_of_vertex);
    builder.add(&mut target, the_vertex);
}

/// `FindEdgeCloseToBisectorPlane` (`cxx:2708-2777`): the edge at
/// `the_vertex` whose direction / plane is closest to the bisector axis, taken
/// out of `the_comp`.
///
/// `the_comp` holds the compound's edges (the reference mutates a
/// `TopoDS_Compound`); `None` is OCCT's null edge.
fn find_edge_close_to_bisector_plane(
    the_vertex: &TopoShape,
    the_comp: &mut Vec<TopoShape>,
    the_axis: &GpAx1,
) -> Option<Edge> {
    // `cxx:2711-2713`: `MapShapesAndAncestors(theComp, VERTEX, EDGE)`.
    let comp = TopoBuilder::new().make_compound_of(the_comp);
    let ve_map = map_shapes_and_ancestors(&comp.0, ShapeType::Vertex, ShapeType::Edge);
    // `cxx:2717-2720`.
    let edges = ve_map.get(&shape_map_key(the_vertex))?;

    // `cxx:2725-2729`.
    let min_edge = if edges.len() == 1 {
        Edge(edges[0].clone())
    } else {
        // `cxx:2731-2774`.
        let reg = GeometryRegistry::global();
        let mut min_angle = f64::MAX;
        let mut best: Option<Edge> = None;
        for an_edge in edges {
            let wire = Wire(an_edge.clone());
            let mut a_pln = GpPln::new(GpAx2::standard().to_ax3());
            let mut issing = false;
            compute_average_plane_and_max_deviation(&wire.0, &mut a_pln, &mut issing);
            let an_angle = if issing {
                // `cxx:2744-2758`: the edge is a line segment.
                let Some(curve) = reg.edge_curve(an_edge) else {
                    continue;
                };
                let (first_par, last_par) = reg.edge_parameters(an_edge);
                let first_pnt = curve.d0(first_par);
                let last_pnt = curve.d0(last_par);
                let edge_vec = GpVec::from_pnts(&first_pnt, &last_pnt);
                let Ok(edge_dir) = GpDir::from_vec(&edge_vec) else {
                    continue;
                };
                let edge_axis = GpAx1::new(first_pnt, edge_dir);
                let mut angle = edge_axis.direction().angle(the_axis.direction());
                if angle > std::f64::consts::FRAC_PI_2 {
                    angle = std::f64::consts::PI - angle;
                }
                std::f64::consts::FRAC_PI_2 - angle
            } else {
                // `cxx:2759-2766`.
                let mut angle = a_pln.axis().direction().angle(the_axis.direction());
                if angle > std::f64::consts::FRAC_PI_2 {
                    angle = std::f64::consts::PI - angle;
                }
                angle
            };
            if an_angle < min_angle {
                min_angle = an_angle;
                best = Some(Edge(an_edge.clone()));
            }
        }
        best?
    };

    // `cxx:2775-2776`: `B.Remove(theComp, MinEdge)`.
    comp_remove(the_comp, &min_edge.0);
    Some(min_edge)
}

/// `FindMiddleEdges` (`cxx:2779-2819`): walk the section compound from
/// `the_vertex1` to `the_vertex2`, always taking the edge closest to the
/// bisector plane.
fn find_middle_edges(
    the_vertex1: &TopoShape,
    the_vertex2: &TopoShape,
    the_axis: &GpAx1,
    the_comp: &mut Vec<TopoShape>,
    the_e_list: &mut Vec<TopoShape>,
) -> bool {
    // `cxx:2786-2794`.
    let comp = TopoBuilder::new().make_compound_of(the_comp);
    let ve_map = map_shapes_and_ancestors(&comp.0, ShapeType::Vertex, ShapeType::Edge);
    if ve_map.is_empty() {
        return false;
    }
    if !ve_map.contains_key(&shape_map_key(the_vertex1))
        || !ve_map.contains_key(&shape_map_key(the_vertex2))
    {
        return false;
    }

    // `cxx:2796-2818`.
    let mut cur_vertex = the_vertex1.clone();
    loop {
        let Some(cur_edge) = find_edge_close_to_bisector_plane(&cur_vertex, the_comp, the_axis)
        else {
            return false;
        };
        let (v1, v2) = edge_vertices(&cur_edge);
        let Some(v1) = v1 else {
            return false;
        };
        let Some(v2) = v2 else {
            return false;
        };
        cur_vertex = if is_same(&v1.0, &cur_vertex) {
            v2.0.clone()
        } else {
            v1.0.clone()
        };
        the_e_list.push(cur_edge.0.clone());
        if is_same(&cur_vertex, the_vertex2) {
            return true;
        }
    }
}

/// `FindCommonVertex(theFirstEdge, theLastEdge, ...)` (`cxx:2821-2855`): the
/// shared vertex of two section edges, going through
/// `TopExp::CommonVertex` unless both ends coincide.
fn find_common_vertex_of_edges(
    the_first_edge: &Edge,
    the_last_edge: &Edge,
    the_first_vertex: Option<&TopoShape>,
    the_last_vertex: Option<&TopoShape>,
) -> Option<TopoShape> {
    // `cxx:2827-2831`: `!theFirstVertex.IsSame(theLastVertex)`.
    let same_ends = match (the_first_vertex, the_last_vertex) {
        (Some(a), Some(b)) => is_same(a, b),
        (None, None) => true,
        _ => false,
    };
    if !same_ends {
        return common_vertex(the_first_edge, the_last_edge);
    }

    // `cxx:2833-2853`.
    let first = the_first_vertex?;
    let (Some(v1), Some(v2)) = edge_vertices(the_first_edge) else {
        return None;
    };
    let (Some(v3), Some(v4)) = edge_vertices(the_last_edge) else {
        return None;
    };
    if is_same(&v1.0, first) {
        if is_same(&v2.0, &v3.0) || is_same(&v2.0, &v4.0) {
            return Some(v2.0.clone());
        }
    } else if is_same(&v1.0, &v3.0) || is_same(&v1.0, &v4.0) {
        return Some(v1.0.clone());
    }
    None
}

/// `RemoveEdges` (`cxx:2520-2547`): the edges of `the_source_comp` that
/// `IsSame` none of `the_list_to_remove`.
fn remove_edges(the_source_comp: &[TopoShape], the_list_to_remove: &[TopoShape]) -> Vec<TopoShape> {
    let mut out: Vec<TopoShape> = Vec::new();
    for an_edge in explorer_shapes_in(the_source_comp, ShapeType::Edge) {
        // `cxx:2531-2538`.
        let found = the_list_to_remove
            .iter()
            .any(|x| is_same(&an_edge, x));
        if !found {
            out.push(an_edge);
        }
    }
    out
}

/// `TopExp_Explorer(compound, TopAbs_EDGE)` over a compound given as an edge
/// list, in list order. The list only holds edges, so this is the identity
/// restricted to `Edge` components.
fn explorer_shapes_in(the_comp: &[TopoShape], to_find: ShapeType) -> Vec<TopoShape> {
    let mut out: Vec<TopoShape> = Vec::new();
    for s in the_comp {
        if s.shape_type() == to_find {
            out.push(s.clone());
        }
    }
    out
}

/// `FilterSectionEdges` (`cxx:2549-2617`).
///
/// `MakeFacesSec` passes `aSecPlane` while it is still null (`cxx:800`,
/// `cxx:802`), so `theSecPlane.IsNull()` holds and the
/// `IntTools_BeanFaceIntersector` arm (`cxx:2576-2606`) never runs: every
/// pave block's section edge is added to the resulting compound
/// (`cxx:2608-2612`).
fn filter_section_edges(the_b_curves: &[BopdsCurve], the_ds: &BopdsDS) -> Vec<TopoShape> {
    let mut out: Vec<TopoShape> = Vec::new();
    for a_b_curve in the_b_curves {
        for a_pb in a_b_curve.pave_blocks() {
            // `cxx:2571-2573`: `nSect = aPB->Edge(); aS = theDS->Shape(nSect)`.
            let Some(a_s) = the_ds.shape(a_pb.edge()) else {
                continue;
            };
            out.push(a_s.clone());
        }
    }
    out
}

/// `GetUEdges` (`cxx:1511-1575`): the two `u`-iso edges of the face column
/// `the_rank`, ordered so that the first one ends where the bound edge starts.
fn get_u_edges(
    the_index: usize,
    the_rank: usize,
    the_u_edges: &ShapeTable,
    the_bound_edge: &TopoShape,
    the_face: &TopoShape,
) -> Option<(Edge, Edge)> {
    // `cxx:1519-1520`.
    let a_ue1 = the_u_edges.value(the_index, 1 + the_rank)?;
    let a_ue2 = the_u_edges.value(the_index + 1, 1 + the_rank)?;

    // `cxx:1522-1538`: the face's own occurrence of each.
    let mut a_face = the_face.clone();
    a_face.set_orientation(Orientation::Forward);
    let mut e1: Option<TopoShape> = None;
    let mut e2: Option<TopoShape> = None;
    for cur in explorer_shapes(&a_face, ShapeType::Edge) {
        if e1.is_none() && is_same(a_ue1, &cur) {
            e1 = Some(cur);
        } else if e2.is_none() && is_same(a_ue2, &cur) {
            e2 = Some(cur);
        }
    }
    // `cxx:1539-1542`.
    let (e1, e2) = (e1?, e2?);

    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(&a_face);
    // `cxx:1545-1551`: `BRep_Tool::CurveOnSurface(E1, aFace, f, l)`.
    let c1 = reg.edge_pcurve(&e1, key)?;
    let (c1_f, c1_l) = reg.pcurve_range(&e1, key)?;
    let pu1 = if the_rank == 0 {
        c1.d0(c1_l)
    } else {
        c1.d0(c1_f)
    };
    // `cxx:1552-1557`.
    let c2 = reg.edge_pcurve(the_bound_edge, key)?;
    // `cxx:1558-1559`: the 3D range, not the pcurve range.
    let (bf, _) = reg.edge_parameters(the_bound_edge);
    let pf = c2.d0(bf);

    // `cxx:1560-1561`.
    let a_v = if the_rank == 0 {
        last_vertex(&Edge(e1.clone()))
    } else {
        first_vertex(&Edge(e1.clone()))
    };
    let a_v = a_v?;
    let a_tolerance = BRepTool::vertex_tolerance(&a_v);

    // `cxx:1562-1568`: `BRepAdaptor_Surface::UResolution(aTolerance)`.
    let Some(surf) = reg.face_surface(&a_face) else {
        return Some((Edge(e1), Edge(e2)));
    };
    let u_tol = u_resolution(surf.as_ref(), a_tolerance);
    let pf_dist = ((pf.x() - pu1.x()).powi(2) + (pf.y() - pu1.y()).powi(2)).sqrt();
    if pf_dist > u_tol {
        Some((Edge(e2), Edge(e1)))
    } else {
        Some((Edge(e1), Edge(e2)))
    }
}

/// `BRepLib_MakeWire(edges)` (`BRepLib_MakeWire.cxx:123-453`) over an explicit
/// edge sequence. `None` is `!IsDone()`, i.e. the build failed and OCCT keeps
/// an empty `Shape()`.
fn wire_from_edges(the_edges: &[Edge]) -> Option<Wire> {    let mut w = crate::brep_lib_make_wire::MakeWire::new();
    for e in the_edges {
        w.add(e).ok()?;
    }
    if w.is_done() {
        Some(w.wire())
    } else {
        None
    }
}

/// `BRepLib_MakeWire(W, E)` (`BRepLib_MakeWire.cxx:214-260`): the edges of
/// `the_wire` followed by `the_edge`.
fn wire_with_edge(the_wire: &TopoShape, the_edge: &Edge) -> Option<Wire> {
    let mut edges = edges_of_wire(&Wire(the_wire.clone()));
    edges.push(the_edge.clone());
    wire_from_edges(&edges)
}

/// `BRepExtrema_ExtCC(E1, E2)` (`BRepExtrema_ExtCC.cxx`): the curve/curve
/// extrema of the two edges over their own parameter ranges, as
/// `(point on E1, point on E2, parameter on E1, parameter on E2)`.
fn extrema_of_edges(the_e1: &Edge, the_e2: &Edge) -> Vec<(GpPnt, GpPnt, f64, f64)> {
    let reg = GeometryRegistry::global();
    let (Some(c1), Some(c2)) = (reg.edge_curve(&the_e1.0), reg.edge_curve(&the_e2.0)) else {
        return Vec::new();
    };
    let (a1, b1) = reg.edge_parameters(&the_e1.0);
    let (a2, b2) = reg.edge_parameters(&the_e2.0);
    occt_geom::extrema_cc::curve_curve_extrema_all_range(
        c1.as_ref(),
        c2.as_ref(),
        a1,
        b1,
        a2,
        b2,
    )
    .into_iter()
    .map(|p| (p.p1, p.p2, p.u1, p.u2))
    .collect()
}

/// `ChooseSection` (`cxx:1020-1285`): the continuous section wire closest to
/// the bisector plane, with its average plane.
///
/// `None` is the reference's `false`; the third field is `IsSingular`.
fn choose_section(
    the_transition: TransitionStyle,
    the_axe_of_bis_plane: &GpAx2,
    the_comp: &[TopoShape],
    the_first_vertex: Option<&TopoShape>,
    the_last_vertex: Option<&TopoShape>,
) -> Option<(TopoShape, GpPln, bool)> {
    // `cxx:1028`.
    let mut is_singular = false;
    let axis = the_axe_of_bis_plane.axis();

    // `cxx:1032-1137`: with `myTransition == BRepFill_Right` and both end
    // vertices known, the section wire runs from the first to the last one.
    if the_transition == TransitionStyle::Right
        && the_first_vertex.is_some()
        && the_last_vertex.is_some()
    {
        let first_v = the_first_vertex?;
        let last_v = the_last_vertex?;

        // `cxx:1039-1045`.
        let mut old_comp: Vec<TopoShape> = the_comp.to_vec();

        // `cxx:1047-1052`.
        let first_edge = find_edge_close_to_bisector_plane(first_v, &mut old_comp, axis)?;

        // `cxx:1054-1059`: the compound must not run dry before the last edge.
        if old_comp.is_empty() {
            if let Some(v) = the_comp.first() {
                old_comp.push(v.clone());
            }
        }

        // `cxx:1060-1064`.
        let last_edge = find_edge_close_to_bisector_plane(last_v, &mut old_comp, axis)?;

        // `cxx:1072`.
        let mut new_wire: Vec<TopoShape> = vec![first_edge.0.clone()];

        // `cxx:1074-1132`.
        if !is_same(&first_edge.0, &last_edge.0) {
            let common =
                find_common_vertex_of_edges(&first_edge, &last_edge, Some(first_v), Some(last_v));
            if common.is_some() {
                // `cxx:1079-1082`.
                new_wire.push(last_edge.0.clone());
            } else {
                // `cxx:1084-1089`: the far end of each edge.
                let (fv1, fv2) = edge_vertices(&first_edge);
                let (fv1, fv2) = (fv1?, fv2?);
                let vertex1 = if is_same(&fv1.0, first_v) {
                    fv2.0.clone()
                } else {
                    fv1.0.clone()
                };
                let (lv1, lv2) = edge_vertices(&last_edge);
                let (lv1, lv2) = (lv1?, lv2?);
                let vertex2 = if is_same(&lv1.0, last_v) {
                    lv2.0.clone()
                } else {
                    lv1.0.clone()
                };

                // `cxx:1090-1101`.
                if find_middle_edges(&vertex1, &vertex2, axis, &mut old_comp, &mut new_wire) {
                    new_wire.push(last_edge.0.clone());
                } else {
                    // `cxx:1103-1131`: trim both edges at their closest points
                    // and join them with a vertex at the midpoint.
                    let sols = extrema_of_edges(&first_edge, &last_edge);
                    if !sols.is_empty() {
                        let mut imin = 0usize;
                        for i in 1..sols.len() {
                            let d_i = sols[i].0.square_distance(&sols[i].1);
                            let d_min = sols[imin].0.square_distance(&sols[imin].1);
                            if d_i < d_min {
                                imin = i;
                            }
                        }
                        let (p_on_first, p_on_last, par_first, par_last) = sols[imin].clone();
                        let a_min_dist = p_on_first.distance(&p_on_last);
                        let mid = GpPnt::new(
                            (p_on_first.x() + p_on_last.x()) * 0.5,
                            (p_on_first.y() + p_on_last.y()) * 0.5,
                            (p_on_first.z() + p_on_last.z()) * 0.5,
                        );
                        // `cxx:1121-1122`.
                        let mut common_vertex = TopoBuilder::new().make_vertex(mid, CONFUSION);
                        common_vertex.set_tolerance(1.001 * a_min_dist / 2.0);
                        // `cxx:1124-1125`.
                        update_section_edge(
                            &first_edge.0,
                            first_v,
                            &mut common_vertex.0,
                            par_first,
                        );
                        update_section_edge(
                            &last_edge.0,
                            last_v,
                            &mut common_vertex.0,
                            par_last,
                        );
                        // `cxx:1127-1129`.
                        new_wire.push(last_edge.0.clone());
                    }
                }
            }
        }

        // `cxx:1134-1136`.
        let edges: Vec<Edge> = new_wire.iter().map(|s| Edge(s.clone())).collect();
        let wire = wire_from_edges(&edges)?;
        return Some((
            wire.0,
            GpPln::new(the_axe_of_bis_plane.to_ax3()),
            is_singular,
        ));
    }

    // `cxx:1139-1148`: the general case walks the section compound edge by edge.
    let mut old_comp: Vec<TopoShape> = the_comp.to_vec();
    let mut an_error = false;
    let mut wseq: Vec<TopoShape> = Vec::new();

    loop {
        // `cxx:1156-1160`.
        let Some(first_shape) = old_comp.first().cloned() else {
            break;
        };
        // `cxx:1161-1163`.
        let first_edge = Edge(first_shape.clone());
        let Some(mut new_wire) = wire_from_edges(&[first_edge.clone()]) else {
            an_error = true;
            break;
        };
        comp_remove(&mut old_comp, &first_shape);
        // `cxx:1164-1168`.
        if new_wire.0.is_closed() {
            wseq.push(new_wire.0);
            continue;
        }

        loop {
            // `cxx:1171-1177`.
            let (ext0, ext1) = wire_vertices(&new_wire);
            let (Some(ext0), Some(ext1)) = (ext0, ext1) else {
                an_error = true;
                break;
            };
            // `cxx:1178-1190`.
            let comp = TopoBuilder::new().make_compound_of(&old_comp);
            let ve_map = map_shapes_and_ancestors(&comp.0, ShapeType::Vertex, ShapeType::Edge);
            let mut vedges: [Vec<TopoShape>; 2] = [Vec::new(), Vec::new()];
            for (j, ext) in [&ext0, &ext1].into_iter().enumerate() {
                if let Some(list) = ve_map.get(&shape_map_key(&ext.0)) {
                    vedges[j] = list.clone();
                }
            }
            // `cxx:1191-1196`.
            if vedges[0].is_empty() && vedges[1].is_empty() {
                break;
            }
            // `cxx:1197-1207`: a single continuation is appended directly.
            let mut modified = false;
            for j in 0..2 {
                if vedges[j].len() == 1 {
                    let an_edge = Edge(vedges[j][0].clone());
                    let Some(w) = wire_with_edge(&new_wire.0, &an_edge) else {
                        an_error = true;
                        break;
                    };
                    new_wire = w;
                    comp_remove(&mut old_comp, &an_edge.0);
                    modified = true;
                }
            }
            if an_error {
                break;
            }
            // `cxx:1208-1233`: several continuations - take the one whose
            // average plane deviates least.
            if !modified {
                let ind = if vedges[0].is_empty() { 1 } else { 0 };
                let edges = vedges[ind].clone();
                let mut the_ind = 0usize;
                let mut min_deviation = f64::MAX;
                for (j, cand) in edges.iter().enumerate() {
                    let Some(a_wire) = wire_with_edge(&new_wire.0, &Edge(cand.clone())) else {
                        break;
                    };
                    let mut a_pln = GpPln::new(GpAx2::standard().to_ax3());
                    let mut issing = false;
                    let deviation =
                        compute_average_plane_and_max_deviation(&a_wire.0, &mut a_pln, &mut issing);
                    if deviation < min_deviation {
                        min_deviation = deviation;
                        the_ind = j;
                    }
                }
                let Some(w) = wire_with_edge(&new_wire.0, &Edge(edges[the_ind].clone())) else {
                    an_error = true;
                    break;
                };
                new_wire = w;
                // `cxx:1232-1233`.
                comp_remove(&mut old_comp, &edges[the_ind]);
            }
            // `cxx:1234-1237`.
            if new_wire.0.is_closed() {
                break;
            }
        }
        // `cxx:1238`.
        wseq.push(new_wire.0);
        if an_error {
            break;
        }
    }

    // `cxx:1246-1248`.
    let mut min_angle = f64::MAX;
    if an_error || !old_comp.is_empty() {
        return None;
    }
    // `cxx:1249-1254`: one wire is the answer.
    if wseq.len() == 1 {
        let res_wire = wseq[0].clone();
        let mut res_plane = GpPln::new(GpAx2::standard().to_ax3());
        compute_average_plane_and_max_deviation(&res_wire, &mut res_plane, &mut is_singular);
        return Some((res_wire, res_plane, is_singular));
    }
    // `cxx:1255-1281`: keep the wire whose average plane is closest to the
    // bisector plane.
    let mut res_wire: Option<TopoShape> = None;
    let mut res_plane = GpPln::new(GpAx2::standard().to_ax3());
    for a_wire in &wseq {
        let mut a_pln = GpPln::new(GpAx2::standard().to_ax3());
        let mut issing = false;
        compute_average_plane_and_max_deviation(a_wire, &mut a_pln, &mut issing);
        if issing {
            continue;
        }
        // `cxx:1266-1273`.
        let mut angle = a_pln.axis().direction().angle(axis.direction());
        if angle > std::f64::consts::FRAC_PI_2 {
            angle = std::f64::consts::PI - angle;
        }
        if angle < min_angle {
            min_angle = angle;
            res_wire = Some(a_wire.clone());
            res_plane = a_pln;
        }
    }
    // `cxx:1249-1283`: with an empty `Wseq` the reference still answers `true`
    // and hands back its untouched (null) wire; the port reports the empty
    // wire so `MakeFacesSec` sees the same empty section compound.
    Some((
        res_wire.unwrap_or_else(|| Wire::new().0),
        res_plane,
        is_singular,
    ))
}

/// `StoreVedgeInHistMap` (`cxx:1393-1407`): the `v`-iso edge of the face
/// starts with a default (empty) entry in the map, then the intersection of
/// the two consecutive faces is appended to it.
pub(super) fn store_vedge_in_hist_map(
    the_v_edges: &[Option<TopoShape>],
    the_index: usize,
    the_new_vshape: &TopoShape,
    the_hist_map: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
) {
    let Some(Some(a_v_edge)) = the_v_edges.get(the_index) else {
        return;
    };
    the_hist_map
        .entry(shape_map_key(a_v_edge))
        .or_default()
        .push(the_new_vshape.clone());
}

/// `FindFreeVertices(theShape, theVerticesToAvoid, theListOfVertex)`
/// (`cxx:1409-1435`): the vertices of `theShape` outside
/// `the_vertices_to_avoid` whose ancestor count is below two.
///
/// `TopExp::MapShapesAndAncestors(theShape, TopAbs_VERTEX, TopAbs_EDGE, aMap)`
/// is the non-unique variant: every `(vertex, edge)` occurrence is one list
/// entry, so an edge whose two endpoints are the same vertex contributes two.
pub(super) fn find_free_vertices(
    the_shape: &TopoShape,
    the_vertices_to_avoid: &[TopoShape],
    the_list_of_vertex: &mut Vec<TopoShape>,
) {
    the_list_of_vertex.clear();
    // Insertion-ordered ancestor counts, matching the `IndexedDataMap` walk of
    // the reference.
    let mut counts: Vec<(TopoShape, usize)> = Vec::new();
    let bump = |v: &TopoShape, counts: &mut Vec<(TopoShape, usize)>| {
        let key = shape_map_key(v);
        if let Some(slot) = counts
            .iter_mut()
            .find(|(s, _)| shape_map_key(s) == key)
        {
            slot.1 += 1;
        } else {
            counts.push((v.clone(), 1));
        }
    };
    for e in explorer_shapes(the_shape, ShapeType::Edge) {
        for v in explorer_shapes(&e, ShapeType::Vertex) {
            bump(&v, &mut counts);
        }
    }
    for (v, n) in counts {
        if the_vertices_to_avoid.iter().any(|a| is_same(a, &v)) {
            continue;
        }
        if n < 2 {
            the_list_of_vertex.push(v);
        }
    }
}

/// `FindCommonVertex(theDS, theEIndex1, theEIndex2, theCrossDirection,
/// theCommonVertex, theParamOnE1, theParamOnE2)` (`cxx:1440-1508`).
///
/// Walks the `InterfEE` records of the pair, keeps the ones whose
/// `IntTools_CommonPrt::Type()` is `TopAbs_VERTEX`, takes the vertex parameters
/// through `IntTools_Tools::VertexParameters` and accepts the record only when
/// the cross product of the two edge tangents follows `theCrossDirection`
/// (`cxx:1494-1501`).
pub(super) fn find_common_vertex(
    the_ds: &BopdsDS,
    the_e_index1: usize,
    the_e_index2: usize,
    the_cross_direction: &GpVec,
    the_common_vertex: &mut Option<TopoShape>,
    the_param_on_e1: &mut f64,
    the_param_on_e2: &mut f64,
) -> bool {
    let Some(e1s) = the_ds.shape(the_e_index1).cloned() else {
        return false;
    };
    let Some(e2s) = the_ds.shape(the_e_index2).cloned() else {
        return false;
    };
    let (Some(c1), Some(c2)) = (
        BRepTool::edge_curve(&Edge(e1s)),
        BRepTool::edge_curve(&Edge(e2s)),
    ) else {
        return false;
    };

    for ee in the_ds.interf_ee() {
        if !((the_e_index1 == ee.index1() && the_e_index2 == ee.index2())
            || (the_e_index1 == ee.index2() && the_e_index2 == ee.index1()))
        {
            continue;
        }
        // `cxx:1469-1472`: records without a new vertex are skipped.
        let Some(n_new) = ee.new_index() else {
            continue;
        };
        let Some(cp) = ee.common_part() else {
            continue;
        };
        if cp.part_type() != CommonPartType::Vertex {
            continue;
        }
        let Some(vs) = the_ds.shape(n_new).cloned() else {
            continue;
        };

        // `cxx:1478-1484`: `IntTools_Tools::VertexParameters(aCP, t1, t2)`,
        // with the argument order following `Index1()`.
        let range1 = cp.range;
        let range2 = cp.ranges2().first().copied().unwrap_or_default();
        let vp1 = cp.vertex_parameter1().unwrap_or(f64::NAN);
        let vp2 = cp.vertex_parameter2().unwrap_or(f64::NAN);
        let (t1, t2) = if the_e_index1 == ee.index1() {
            crate::int_tools_curve_box::curve_box_vertex_parameters(range1, vp1, range2, vp2)
        } else {
            crate::int_tools_curve_box::curve_box_vertex_parameters(range2, vp2, range1, vp1)
        };

        // `cxx:1488-1501`: the intersection direction must follow the path.
        let dir1 = c1.d1(t1).1;
        let dir2 = c2.d1(t2).1;
        let cross = dir1.crossed(&dir2);
        if cross.dot(the_cross_direction) > SQUARE_CONFUSION {
            *the_common_vertex = Some(vs);
            *the_param_on_e1 = t1;
            *the_param_on_e2 = t2;
            return true;
        }
    }
    false
}

/// `SplitUEdges(theUEdges, theDS, theCrossDirection, theHistMap)`
/// (`cxx:1289-1392`).
///
/// For every row (the pair of consecutive `u`-iso edges of one face) that is
/// not already replaced in the map: look for a common vertex on the two edges;
/// when the `EE` search fails, fall back to the `VV` interference of the row's
/// `LastVertex(e1)` / `FirstVertex(e2)`. On success the two edges are re-cut at
/// the new vertex and bound in the map.
fn split_u_edges(
    the_u_edges: &ShapeTable,
    the_ds: &BopdsDS,
    the_cross_direction: &GpVec,
    the_hist_map: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
) -> bool {
    let low_col = 1usize;
    let up_col = the_u_edges.nb_cols();
    let up_row = the_u_edges.nb_rows();

    for ueit in 1..=up_row {
        let Some(a_e1) = the_u_edges.value(ueit, low_col).cloned() else {
            continue;
        };
        let Some(a_e2) = the_u_edges.value(ueit, up_col).cloned() else {
            continue;
        };
        if the_hist_map.contains_key(&shape_map_key(&a_e1))
            || the_hist_map.contains_key(&shape_map_key(&a_e2))
        {
            continue;
        }
        let (Some(an_e_index1), Some(an_e_index2)) = (the_ds.index(&a_e1), the_ds.index(&a_e2))
        else {
            continue;
        };

        let mut a_common_vertex: Option<TopoShape> = None;
        let mut apar1 = 0.0f64;
        let mut apar2 = 0.0f64;
        let mut b_vertex_found = find_common_vertex(
            the_ds,
            an_e_index1,
            an_e_index2,
            the_cross_direction,
            &mut a_common_vertex,
            &mut apar1,
            &mut apar2,
        );

        if !b_vertex_found {
            // `cxx:1327-1353`: `VV` fallback on the row's end/start vertices.
            let v1 = edge_vertices(&Edge(a_e1.clone())).1;
            let v2 = edge_vertices(&Edge(a_e2.clone())).0;
            if let (Some(v1), Some(v2)) = (v1, v2) {
                let vindex1 = the_ds.index(&v1.0).unwrap_or(usize::MAX);
                let vindex2 = the_ds.index(&v2.0).unwrap_or(usize::MAX);
                for vv in the_ds.interf_vv() {
                    let matches = (vindex1 == vv.index1() && vindex2 == vv.index2())
                        || (vindex1 == vv.index2() && vindex2 == vv.index1());
                    if !matches {
                        continue;
                    }
                    let Some(n_new) = vv.new_index() else {
                        continue;
                    };
                    a_common_vertex = the_ds.shape(n_new).cloned();
                    // `cxx:1349-1350`: `BRep_Tool::Parameter(V1, Edge(aE1))`
                    // and the same on `aE2`.
                    apar1 = BRepTool::parameter_of_vertex(&a_e1, &v1.0);
                    apar2 = BRepTool::parameter_of_vertex(&a_e2, &v2.0);
                    b_vertex_found = true;
                    break;
                }
            }
        }

        if !b_vertex_found {
            continue;
        }
        let Some(a_common_vertex) = a_common_vertex else {
            continue;
        };

        // `cxx:1357-1385`: `EmptyCopied` + `Add(E, V1)` + `Add(E, V2)` +
        // `Range(E, f, l)`. The copy keeps the geometry of the source edge, so
        // the port rebuilds it on the same curve with the new bounds.
        let (Some(curve1), Some(curve2)) = (
            BRepTool::edge_curve(&Edge(a_e1.clone())),
            BRepTool::edge_curve(&Edge(a_e2.clone())),
        ) else {
            continue;
        };
        let (f1, _) = BRepTool::edge_parameters(&Edge(a_e1.clone()));
        let (_, l2) = BRepTool::edge_parameters(&Edge(a_e2.clone()));
        let verts1 = edge_vertices(&Edge(a_e1.clone()));
        let verts2 = edge_vertices(&Edge(a_e2.clone()));
        let (Some(a_v1), Some(a_v2)) = (verts1.0, verts2.1) else {
            continue;
        };

        let Some(new_e1) = rebuild_bounded_edge(
            &curve1,
            &a_v1,
            &Vertex(a_common_vertex.clone()),
            f1,
            apar1,
        ) else {
            continue;
        };
        let Some(new_e2) = rebuild_bounded_edge(
            &curve2,
            &Vertex(a_common_vertex.clone()),
            &a_v2,
            apar2,
            l2,
        ) else {
            continue;
        };

        the_hist_map
            .insert(shape_map_key(&a_e1), vec![new_e1.0.clone()]);
        the_hist_map
            .insert(shape_map_key(&a_e2), vec![new_e2.0.clone()]);
    }
    true
}

/// The `B.Begin/Add(E, V)/Range(E, f, l)` sequence of
/// `BRepFill_TrimShellCorner.cxx:1360-1385` on the curve of an existing edge.
fn rebuild_bounded_edge(
    curve: &std::sync::Arc<dyn occt_geom::curve::Curve>,
    v_first: &Vertex,
    v_last: &Vertex,
    first: f64,
    last: f64,
) -> Option<Edge> {
    let b = TopoBuilder::new();
    let build = init_3d_vertices(curve, Some(v_first), Some(v_last), first, last)?;
    let mut e = b.make_edge(curve.clone(), build.first, build.last);
    add_endpoint_vertices(&b, &mut e, &build);
    Some(e)
}

// ---------------------------------------------------------------------------
// Section-walk statics (`FillGap` and its chain, `cxx:1577-2519`)
// ---------------------------------------------------------------------------

/// `shape.Oriented(theOrientation)` (`TopoDS_Shape::Oriented`): a copy of
/// `shape` carrying `the_orientation`.
fn oriented_to(shape: &TopoShape, the_orientation: Orientation) -> TopoShape {
    let mut s = shape.clone();
    s.set_orientation(the_orientation);
    s
}

/// `TopExp::Vertices(E, Vfirst, Vlast, /*CumOri=*/true)` (`TopExp.cxx:214-252`).
///
/// The component orientations are composed with the edge's own orientation
/// (`TopoDS_Iterator(aShape, CumOri=true)`, `TopoDS_Iterator.cxx:32-46`) before
/// the FORWARD/REVERSED split, so a REVERSED edge reports the *stored* REVERSED
/// vertex first. `CheckAndOrientEdges` relies on exactly that (`cxx:1785`).
fn edge_vertices_cum_ori(edge: &Edge) -> (Option<Vertex>, Option<Vertex>) {
    let ori = edge.0.orientation();
    let stored = edge
        .0
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .clone();
    let mut first: Option<Vertex> = None;
    let mut last: Option<Vertex> = None;
    for child in stored {
        if child.shape_type() != ShapeType::Vertex {
            continue;
        }
        let cum = Orientation::compose(ori, child.orientation());
        if cum == Orientation::Forward && first.is_none() {
            first = Some(Vertex(child));
        } else if cum == Orientation::Reversed && last.is_none() {
            last = Some(Vertex(child));
        }
    }
    (first, last)
}

/// `BRepAdaptor_Surface(theFace).UResolution(tol)` / `VResolution(tol)`
/// (`BRepAdaptor_Surface.cxx:1410-1430` delegate to
/// `GeomAdaptor_Surface::UResolution`, `:1818-1896`): the reference takes the
/// larger of the two as the working 2D tolerance (`cxx:1729-1735`).
fn surface_tolerance(the_face: &TopoShape, tol: f64) -> f64 {
    let reg = GeometryRegistry::global();
    let Some(surf) = reg.face_surface(the_face) else {
        return tol;
    };
    let u = occt_geom::approx_same_parameter::u_resolution(surf.as_ref(), tol);
    let v = occt_geom::approx_same_parameter::v_resolution(surf.as_ref(), tol);
    if u > v {
        u
    } else {
        v
    }
}

/// `CheckAndOrientEdges` (`cxx:1705-1827`): orient the section edges into one
/// chain and report whether their pcurve ends reach both gap endpoints.
///
/// `the_ordered_list` is the chain from `FindNextEdge`; `true` means both the
/// first and the last 2D point of the chain land on the gap's endpoints within
/// the surface resolution of the adjacent vertex tolerance.
fn check_and_orient_edges(
    the_ordered_list: &[TopoShape],
    the_first_point: &GpPnt2d,
    the_last_point: &GpPnt2d,
    the_face: &TopoShape,
    the_oriented_list: &mut Vec<TopoShape>,
) -> bool {
    // `cxx:1712-1715`.
    if the_ordered_list.is_empty() {
        return true;
    }
    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(the_face);

    // `cxx:1718-1721`.
    let mut a_e_prev = Edge(the_ordered_list[0].clone());

    // `cxx:1723-1735`.
    let Some(a_curve) = reg.edge_pcurve(&a_e_prev.0, key) else {
        return false;
    };
    let (pf, pl) = reg.pcurve_range(&a_e_prev.0, key).unwrap_or((0.0, 0.0));
    let (vf, vl) = edge_vertices(&a_e_prev);
    let a_tolerance1 = match &vf {
        Some(v) => surface_tolerance(the_face, BRepTool::vertex_tolerance(v)),
        None => surface_tolerance(the_face, CONFUSION),
    };
    let a_tolerance2 = match &vl {
        Some(v) => surface_tolerance(the_face, BRepTool::vertex_tolerance(v)),
        None => surface_tolerance(the_face, CONFUSION),
    };

    // `cxx:1737-1777`.
    let mut ap = a_curve.d0(pf);
    let mut first_found = false;
    let mut last_found = false;
    if pnt2d_distance(&ap, the_first_point) < a_tolerance1 {
        if the_oriented_list.is_empty() {
            the_oriented_list.push(oriented_to(&a_e_prev.0, Orientation::Forward));
        }
        first_found = true;
    } else if pnt2d_distance(&ap, the_last_point) < a_tolerance1 {
        if the_oriented_list.is_empty() {
            the_oriented_list.push(oriented_to(&a_e_prev.0, Orientation::Reversed));
        }
        last_found = true;
    }
    ap = a_curve.d0(pl);
    if pnt2d_distance(&ap, the_last_point) < a_tolerance2 {
        if the_oriented_list.is_empty() {
            the_oriented_list.push(oriented_to(&a_e_prev.0, Orientation::Forward));
        }
        last_found = true;
    } else if pnt2d_distance(&ap, the_first_point) < a_tolerance2 {
        if the_oriented_list.is_empty() {
            the_oriented_list.push(oriented_to(&a_e_prev.0, Orientation::Reversed));
        }
        first_found = true;
    }

    // `cxx:1779-1782`.
    if !the_oriented_list.is_empty() {
        a_e_prev = Edge(the_oriented_list[the_oriented_list.len() - 1].clone());
    }

    // `cxx:1784-1823`.
    for a_shape in the_ordered_list.iter().skip(1) {
        let a_e = Edge(a_shape.clone());
        let (a_v11, a_v12) = edge_vertices_cum_ori(&a_e_prev);
        let (a_v21, a_v22) = edge_vertices(&a_e);

        // `cxx:1788-1790`.
        let closed = match (&a_v12, &a_v21, &a_v11, &a_v22) {
            (Some(v12), Some(v21), _, _) if is_same(&v12.0, &v21.0) => true,
            (_, _, Some(v11), Some(v22)) if is_same(&v11.0, &v22.0) => true,
            _ => false,
        };
        let an_ori = if closed {
            Orientation::Forward
        } else {
            Orientation::Reversed
        };
        the_oriented_list.push(oriented_to(&a_e.0, an_ori));
        a_e_prev = Edge(the_oriented_list[the_oriented_list.len() - 1].clone());

        // `cxx:1794-1801`.
        let a_tolerance1 = match &a_v21 {
            Some(v) => surface_tolerance(the_face, BRepTool::vertex_tolerance(v)),
            None => surface_tolerance(the_face, CONFUSION),
        };
        let a_tolerance2 = match &a_v22 {
            Some(v) => surface_tolerance(the_face, BRepTool::vertex_tolerance(v)),
            None => surface_tolerance(the_face, CONFUSION),
        };
        // `cxx:1802-1822`.
        let Some(a_curve) = reg.edge_pcurve(&a_e.0, key) else {
            return false;
        };
        let (pf, pl) = reg.pcurve_range(&a_e.0, key).unwrap_or((0.0, 0.0));
        let ap = a_curve.d0(pf);
        if pnt2d_distance(&ap, the_first_point) < a_tolerance1 {
            first_found = true;
        } else if pnt2d_distance(&ap, the_last_point) < a_tolerance1 {
            last_found = true;
        }
        let ap = a_curve.d0(pl);
        if pnt2d_distance(&ap, the_first_point) < a_tolerance2 {
            first_found = true;
        } else if pnt2d_distance(&ap, the_last_point) < a_tolerance2 {
            last_found = true;
        }
    }

    first_found && last_found
}

/// `gp_Pnt2d::Distance` (`gp_Pnt2d.cxx`).
fn pnt2d_distance(a: &GpPnt2d, b: &GpPnt2d) -> f64 {
    let dx = a.x() - b.x();
    let dy = a.y() - b.y();
    (dx * dx + dy * dy).sqrt()
}

/// `FindVertex` (`cxx:1829-1887`): the vertex and pave of a `u`-iso edge at
/// rank `the_rank`, taken from its first `myHistMap` image.
fn find_vertex(
    the_edge: &TopoShape,
    the_rank: usize,
    the_ds: &BopdsDS,
    the_hist_map: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    the_vertex: &mut Option<TopoShape>,
    the_pave: &mut BopdsPave,
) -> bool {
    // `cxx:1839-1848`.
    let Some(lst) = the_hist_map.get(&shape_map_key(the_edge)) else {
        return false;
    };
    let Some(first) = lst.first() else {
        return false;
    };

    // `cxx:1851-1866`.
    let a_new_edge = Edge(first.clone());
    let (f, l) = GeometryRegistry::global().edge_parameters(&a_new_edge.0);
    let v = if the_rank == 0 {
        the_pave.set_parameter(l);
        last_vertex(&a_new_edge)
    } else {
        the_pave.set_parameter(f);
        first_vertex(&a_new_edge)
    };
    *the_vertex = v.map(|v| v.0);
    let vertex = match the_vertex.as_ref() {
        Some(v) => v,
        None => return true, // `cxx:2001-2004` guards the null vertex.
    };

    // `cxx:1867-1881`: the index in the DS, searching the new shapes when the
    // vertex is not one of the source shapes.
    let an_index = match the_ds.index(vertex) {
        Some(i) => i,
        None => {
            let mut found = usize::MAX;
            for i in the_ds.nb_source_shapes()..the_ds.nb_shapes() {
                if let Some(s) = the_ds.shape(i) {
                    if is_same(s, vertex) {
                        found = i;
                        break;
                    }
                }
            }
            found
        }
    };
    if an_index == usize::MAX {
        // The reference stores `-1`; no caller can then read a shape back from
        // `theDS->Shape(-1)`, so the walk reports failure where OCCT reaches a
        // null vertex. UNPORTED: `BOPDS_DS::Shape(-1)` has no Rust equivalent.
        the_pave.set_index(usize::MAX);
    } else {
        the_pave.set_index(an_index);
    }
    true
}

/// `FindNextVertex` (`cxx:1889-1931`): the first pave after `the_prev_pave` on
/// the bound edge whose shape is a *new* shape of the DS.
fn find_next_vertex(
    the_edge_index: usize,
    the_prev_pave: &BopdsPave,
    the_ds: &BopdsDS,
    the_next_vertex: &mut Option<TopoShape>,
    the_pave: &mut BopdsPave,
) -> bool {
    // `cxx:1900-1902`: `BOPDS_Pave anullpave;` is the default-constructed
    // (index 0, parameter 0) pave.
    let null_pave = BopdsPave::new(0, 0.0);
    let mut b_take_pave = the_prev_pave.is_equal(&null_pave);
    let mut b_found = false;

    // `cxx:1904-1922`.
    for a_tmp_pave in the_ds.paves(the_edge_index) {
        if b_take_pave {
            if the_ds.is_new_shape(a_tmp_pave.index()) {
                *the_next_vertex = the_ds.shape(a_tmp_pave.index()).cloned();
                *the_pave = a_tmp_pave;
                b_found = true;
                break;
            }
        } else if a_tmp_pave.is_equal(the_prev_pave) {
            b_take_pave = true;
        }
    }
    b_found
}

/// `GetPave` (`cxx:1933-1958`): the first or the last pave of the edge.
fn get_pave(
    the_edge_index: usize,
    is_first: bool,
    the_ds: &BopdsDS,
    the_pave: &mut BopdsPave,
) -> bool {
    let a_lp = the_ds.paves(the_edge_index);
    if a_lp.is_empty() {
        return false;
    }
    *the_pave = if is_first {
        a_lp[0]
    } else {
        a_lp[a_lp.len() - 1]
    };
    true
}

/// `FindNextEdge` (`cxx:1616-1703`): the edge chain that walks from
/// `the_first_vertex` to `the_last_vertex` inside `the_map_ve`.
///
/// The recursion carries its own copy of `the_map_to_avoid` (`cxx:1626-1627`),
/// so a dead end unwinds with the caller's list cleared and its cursor
/// restored (`cxx:1669-1672`).
#[allow(clippy::too_many_arguments)]
fn find_next_edge(
    the_first_vertex: &TopoShape,
    the_last_vertex: &TopoShape,
    the_map_ve: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    the_map_to_avoid: &[ShapeMapKey],
    the_ordered_list: &mut Vec<TopoShape>,
) -> bool {
    // `cxx:1626-1629`.
    let mut a_cur_vertex = the_first_vertex.clone();
    let mut a_map_to_avoid: Vec<ShapeMapKey> = the_map_to_avoid.to_vec();
    let mut a_list_of_edge: Vec<TopoShape> = Vec::new();
    let extent = the_map_ve.len();

    // `cxx:1630-1698`.
    for _ in 0..extent {
        let Some(lste) = the_map_ve.get(&shape_map_key(&a_cur_vertex)) else {
            break; // `cxx:1632-1635`.
        };
        let mut b_efound = false;
        for an_edge in lste.clone() {
            let a_save_cur_vertex = a_cur_vertex.clone();
            if a_map_to_avoid.contains(&shape_map_key(&an_edge)) {
                continue;
            }
            // `cxx:1647-1659`: move the cursor to the far end.
            let (v1, v2) = edge_vertices(&Edge(an_edge.clone()));
            if let Some(v1) = &v1 {
                if !is_same(&a_cur_vertex, &v1.0) {
                    a_cur_vertex = v1.0.clone();
                } else if let Some(v2) = &v2 {
                    if !is_same(&a_cur_vertex, &v2.0) {
                        a_cur_vertex = v2.0.clone();
                    }
                }
            } else if let Some(v2) = &v2 {
                if !is_same(&a_cur_vertex, &v2.0) {
                    a_cur_vertex = v2.0.clone();
                }
            }
            a_map_to_avoid.push(shape_map_key(&an_edge));
            b_efound = true;
            a_list_of_edge.push(an_edge);

            // `cxx:1663-1680`.
            if !is_same(&a_cur_vertex, the_last_vertex) {
                let mut a_list_tmp: Vec<TopoShape> = Vec::new();
                if !find_next_edge(
                    &a_cur_vertex,
                    the_last_vertex,
                    the_map_ve,
                    &a_map_to_avoid,
                    &mut a_list_tmp,
                ) {
                    a_list_of_edge.clear();
                    a_cur_vertex = a_save_cur_vertex;
                    continue;
                }
                a_list_of_edge.extend(a_list_tmp);
                the_ordered_list.extend(a_list_of_edge);
                return true;
            }
            break;
        }

        // `cxx:1683-1686`.
        if is_same(&a_cur_vertex, the_last_vertex) {
            break;
        }
        // `cxx:1688-1691`.
        if !b_efound {
            return false;
        }
    }

    // `cxx:1694-1697`.
    if is_same(&a_cur_vertex, the_last_vertex) {
        the_ordered_list.extend(a_list_of_edge);
        return true;
    }
    false
}

/// `FillGap` (`cxx:1577-1613`): the oriented section edges that bridge the 2D
/// points `the_first_point` / `the_last_point` between the two vertices.
fn fill_gap(
    the_first_vertex: &TopoShape,
    the_last_vertex: &TopoShape,
    the_first_point: &GpPnt2d,
    the_last_point: &GpPnt2d,
    the_face: &TopoShape,
    the_section_edges: &[TopoShape],
    the_ordered_list: &mut Vec<TopoShape>,
) -> bool {
    // `cxx:1587-1593`: the `(VERTEX, EDGE)` ancestors of the section compound.
    let comp = TopoBuilder::new().make_compound_of(the_section_edges);
    let a_map = map_shapes_and_ancestors(&comp.0, ShapeType::Vertex, ShapeType::Edge);
    if a_map.is_empty() {
        return false;
    }
    // `cxx:1595-1598`.
    if !a_map.contains_key(&shape_map_key(the_first_vertex))
        || !a_map.contains_key(&shape_map_key(the_last_vertex))
    {
        return false;
    }

    // `cxx:1599-1612`.
    let mut a_list_of_edge: Vec<TopoShape> = Vec::new();
    let a_map_to_avoid: Vec<ShapeMapKey> = Vec::new();
    if find_next_edge(
        the_first_vertex,
        the_last_vertex,
        &a_map,
        &a_map_to_avoid,
        &mut a_list_of_edge,
    ) && !a_list_of_edge.is_empty()
    {
        return check_and_orient_edges(
            &a_list_of_edge,
            the_first_point,
            the_last_point,
            the_face,
            the_ordered_list,
        );
    }
    false
}
/// `FindFromUEdge` outputs (`cxx:1962-2070`).
struct FromUEdge {
    /// `theSecEdgesNew`.
    comp_of_sec_edges: Vec<TopoShape>,
    /// `theListOfWireEdges`.
    list_of_wire_edges: Vec<TopoShape>,
    /// `theFoundPave`.
    found_pave: BopdsPave,
    /// `isOnUEdge`.
    is_on_u_edge: bool,
}

/// `FindFromUEdge` (`cxx:1962-2070`): the section edges that bridge the first
/// `u`-iso edge's end to the first reachable pave of the bound edge.
///
/// `None` is the reference's `false`. The `cxx:1998-2000` early `return true`
/// leaves the outputs at the caller's values, i.e. the incoming section
/// compound and an empty wire-edge list, with `theFoundPave` reset to `(0, 0)`.
#[allow(clippy::too_many_arguments)]
fn find_from_u_edge(
    the_ue1_old: &TopoShape,
    the_ue2_old: &TopoShape,
    the_ue1_new: &TopoShape,
    the_ue2_new: &TopoShape,
    the_face: &TopoShape,
    the_sec_edges: &[TopoShape],
    the_rank: usize,
    the_bound_edge: &TopoShape,
    the_bound_edge_index: usize,
    the_ds: &BopdsDS,
    the_hist_map: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    the_sec_edges_in: Vec<TopoShape>,
    the_list_in: Vec<TopoShape>,
) -> Option<FromUEdge> {
    // `cxx:1980-1982`.
    let mut found_pave = BopdsPave::new(0, 0.0);
    let mut is_on_u_edge = true;

    // `cxx:1984-1988`.
    let mut a_face_f = the_face.clone();
    a_face_f.set_orientation(Orientation::Forward);
    let mut comp_of_sec_edges = the_sec_edges_in;
    let mut list_of_wire_edges = the_list_in;

    // `cxx:1992-1996`.
    let mut a_pave1 = BopdsPave::new(0, 0.0);
    let mut a_pave2 = BopdsPave::new(0, 0.0);
    let mut a_prev_vertex: Option<TopoShape> = None;
    let mut a_next_vertex: Option<TopoShape> = None;
    let mut atmp_pave = BopdsPave::new(0, 0.0);

    // `cxx:1998-2000`: an unrecorded `u`-iso edge leaves every output alone.
    if !find_vertex(
        the_ue1_old,
        the_rank,
        the_ds,
        the_hist_map,
        &mut a_prev_vertex,
        &mut atmp_pave,
    ) {
        return Some(FromUEdge {
            comp_of_sec_edges,
            list_of_wire_edges,
            found_pave,
            is_on_u_edge,
        });
    }
    // `cxx:2002-2005`.
    let Some(a_prev) = a_prev_vertex.clone() else {
        return None;
    };

    // `cxx:2007-2010`.
    let a_first_v = a_prev.clone();
    let mut b_sec_found = false;
    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(&a_face_f);
    let Some(a_c1) = reg.edge_pcurve(the_ue1_new, key) else {
        return None;
    };
    let (c1_f, c1_l) = reg.pcurve_range(the_ue1_new, key).unwrap_or((0.0, 0.0));
    let mut p1 = if the_rank == 0 {
        a_c1.d0(c1_l)
    } else {
        a_c1.d0(c1_f)
    };

    // `cxx:2013-2017`.
    let nbpave = the_ds.paves(the_bound_edge_index).len();
    let mut pit = 0usize;

    // `cxx:2018-2039`: walk the bound edge's paves.
    let Some(a_c2) = reg.edge_pcurve(the_bound_edge, key) else {
        return None;
    };
    while find_next_vertex(
        the_bound_edge_index,
        &a_pave1,
        the_ds,
        &mut a_next_vertex,
        &mut a_pave2,
    ) && pit < nbpave
    {
        let Some(a_last_v) = a_next_vertex.clone() else {
            return None;
        };
        let p2 = a_c2.d0(a_pave2.parameter());
        let mut a_ordered_list: Vec<TopoShape> = Vec::new();

        if fill_gap(
            &a_first_v,
            &a_last_v,
            &p1,
            &p2,
            &a_face_f,
            &comp_of_sec_edges,
            &mut a_ordered_list,
        ) {
            comp_of_sec_edges = remove_edges(&comp_of_sec_edges, &a_ordered_list);
            list_of_wire_edges.extend(a_ordered_list);
            found_pave = a_pave2;
            is_on_u_edge = false;
            b_sec_found = true;
            break;
        }
        if let Some(v) = a_next_vertex.clone() {
            a_prev_vertex = Some(v);
        }
        a_pave1 = a_pave2;
        pit += 1;
    }

    // `cxx:2041-2060`: fall back to the second `u`-iso edge.
    if !b_sec_found {
        let mut a_next2: Option<TopoShape> = None;
        let mut a_pave2b = BopdsPave::new(0, 0.0);
        if find_vertex(
            the_ue2_old,
            the_rank,
            the_ds,
            the_hist_map,
            &mut a_next2,
            &mut a_pave2b,
        ) {
            if let Some(a_last_v) = a_next2.clone() {
                let Some(a_c3) = reg.edge_pcurve(the_ue2_new, key) else {
                    return None;
                };
                let p2 = a_c3.d0(a_pave2b.parameter());
                let mut a_ordered_list: Vec<TopoShape> = Vec::new();
                if fill_gap(
                    &a_first_v,
                    &a_last_v,
                    &p1,
                    &p2,
                    &a_face_f,
                    &comp_of_sec_edges,
                    &mut a_ordered_list,
                ) {
                    comp_of_sec_edges = remove_edges(&comp_of_sec_edges, &a_ordered_list);
                    list_of_wire_edges.extend(a_ordered_list);
                    found_pave = a_pave2b;
                    b_sec_found = true;
                    is_on_u_edge = true;
                }
            }
        }
    }

    // `cxx:2062-2068`: without a section the outputs stay as they came in.
    if !b_sec_found {
        comp_of_sec_edges = the_sec_edges.to_vec();
        list_of_wire_edges.clear();
    }
    Some(FromUEdge {
        comp_of_sec_edges,
        list_of_wire_edges,
        found_pave,
        is_on_u_edge,
    })
}

/// `FindFromVEdge` outputs (`cxx:2073-2513`).
struct FromVEdge {
    /// `theListOfWireEdges`.
    list_of_wire_edges: Vec<TopoShape>,
    /// `isSectionFound`.
    is_section_found: bool,
}

/// `FindFromVEdge` (`cxx:2073-2513`): the section edges that complete the wire
/// from the pave found by `FindFromUEdge` down to the bound edge's last pave.
///
/// The walk exists twice in the reference: a "by pairs non continuously" pass
/// (`cxx:2144-2226`) that collects candidate `(pave1, pave2, flag)` triples and
/// a "by pairs continuously" pass (`cxx:2228-2511`) that consumes them. Both
/// are reproduced in order; `anArrayOfListOfSec` stays 1-based as in OCCT.
#[allow(clippy::too_many_arguments)]
fn find_from_v_edge(
    the_prev_pave: &BopdsPave,
    is_on_u_edge: bool,
    the_ue1_old: &TopoShape,
    the_ue2_old: &TopoShape,
    the_face: &TopoShape,
    the_sec_edges: &[TopoShape],
    the_rank: usize,
    the_bound_edge: &TopoShape,
    the_bound_edge_index: usize,
    the_ds: &BopdsDS,
    the_hist_map: &HashMap<ShapeMapKey, Vec<TopoShape>>,
) -> Option<FromVEdge> {
    // `cxx:2090-2097`.
    let mut list_of_wire_edges: Vec<TopoShape> = Vec::new();
    let mut is_section_found = false;
    let mut a_face_f = the_face.clone();
    a_face_f.set_orientation(Orientation::Forward);
    let mut comp_of_sec_edges: Vec<TopoShape> = the_sec_edges.to_vec();

    // `cxx:2102-2113`: when the previous pave sits on the bound edge, a `u`-iso
    // end matching it means the walk is already closed.
    if is_on_u_edge {
        let mut atmp_vertex: Option<TopoShape> = None;
        let mut a_pave_of_e2 = BopdsPave::new(0, 0.0);
        if find_vertex(
            the_ue2_old,
            the_rank,
            the_ds,
            the_hist_map,
            &mut atmp_vertex,
            &mut a_pave_of_e2,
        ) && the_prev_pave.is_equal(&a_pave_of_e2)
        {
            return Some(FromVEdge {
                list_of_wire_edges,
                is_section_found,
            });
        }
    }

    let reg = GeometryRegistry::global();
    let key = GeometryRegistry::shape_key(&a_face_f);
    let Some(a_c2) = reg.edge_pcurve(the_bound_edge, key) else {
        return None;
    };

    // `cxx:2123-2136`.
    let mut a_pave1 = *the_prev_pave;
    if is_on_u_edge {
        let mut atmp_pave = BopdsPave::new(0, 0.0);
        if !get_pave(the_bound_edge_index, true, the_ds, &mut atmp_pave) {
            return None;
        }
        a_pave1 = atmp_pave;
    }
    let mut p1 = a_c2.d0(a_pave1.parameter());
    let mut a_prev_vertex: Option<TopoShape> = the_ds.shape(a_pave1.index()).cloned();

    // `cxx:2138-2142`.
    let nbpave = the_ds.paves(the_bound_edge_index).len();
    let mut pit = 0usize;
    // `NCollection_Array1<NCollection_List<TopoDS_Shape>>(1, nbpave)`.
    let mut an_array_of_list_of_sec: Vec<Vec<TopoShape>> = vec![Vec::new(); nbpave + 1];

    // `by pairs non continuously. begin` (`cxx:2144-2226`).
    let mut a_first_pave = a_pave1;
    let mut a_first_vertex = a_prev_vertex.clone();
    let mut apfirst = p1;
    let mut a_first_paves: Vec<BopdsPave> = Vec::new();
    let mut a_last_paves: Vec<BopdsPave> = Vec::new();
    let mut a_list_of_flags: Vec<i32> = Vec::new();
    let mut apaircounter = 1usize;

    for _k in 0..nbpave {
        a_pave1 = a_first_pave;
        p1 = apfirst;
        a_prev_vertex = a_first_vertex.clone();
        let mut bfound = false;
        pit = 0;
        let mut a_pave2 = BopdsPave::new(0, 0.0);
        let mut a_next_vertex: Option<TopoShape> = None;

        while find_next_vertex(
            the_bound_edge_index,
            &a_pave1,
            the_ds,
            &mut a_next_vertex,
            &mut a_pave2,
        ) && pit < nbpave
        {
            let Some(a_first_v) = a_prev_vertex.clone() else {
                return None;
            };
            let Some(a_last_v) = a_next_vertex.clone() else {
                return None;
            };
            let p2 = a_c2.d0(a_pave2.parameter());
            let mut a_ordered_list: Vec<TopoShape> = Vec::new();

            if fill_gap(
                &a_first_v,
                &a_last_v,
                &p1,
                &p2,
                &a_face_f,
                &comp_of_sec_edges,
                &mut a_ordered_list,
            ) {
                comp_of_sec_edges = remove_edges(&comp_of_sec_edges, &a_ordered_list);
                if apaircounter <= nbpave {
                    an_array_of_list_of_sec[apaircounter].extend(a_ordered_list);
                }
                apaircounter += 1;
                a_first_paves.push(a_first_pave);
                a_last_paves.push(a_pave2);
                a_list_of_flags.push(1);
                a_first_pave = a_pave2;
                a_first_vertex = a_next_vertex.clone();
                apfirst = p2;
                a_prev_vertex = a_next_vertex.clone();
                is_section_found = true;
                bfound = true;
            }
            a_pave1 = a_pave2;
            pit += 1;
        }

        // `cxx:2190-2210`: the closing pair against the second `u`-iso edge.
        let mut a_pave2v = BopdsPave::new(0, 0.0);
        let mut a_next2: Option<TopoShape> = None;
        if find_vertex(
            the_ue2_old,
            the_rank,
            the_ds,
            the_hist_map,
            &mut a_next2,
            &mut a_pave2v,
        ) {
            if let (Some(a_first_v), Some(a_last_v)) =
                (a_prev_vertex.clone(), a_next2.clone())
            {
                let Some(a_c3) = reg.edge_pcurve(the_ue2_old, key) else {
                    return None;
                };
                let p2 = a_c3.d0(a_pave2v.parameter());
                let mut a_ordered_list: Vec<TopoShape> = Vec::new();
                if fill_gap(
                    &a_first_v,
                    &a_last_v,
                    &p1,
                    &p2,
                    &a_face_f,
                    &comp_of_sec_edges,
                    &mut a_ordered_list,
                ) {
                    comp_of_sec_edges = remove_edges(&comp_of_sec_edges, &a_ordered_list);
                    if apaircounter <= nbpave {
                        an_array_of_list_of_sec[apaircounter].extend(a_ordered_list);
                    }
                    apaircounter += 1;
                    a_first_paves.push(a_first_pave);
                    a_last_paves.push(a_pave2v);
                    a_list_of_flags.push(0);
                    is_section_found = true;
                    break;
                }
            }
        }

        // `cxx:2212-2224`: advance the anchor pave.
        if !bfound {
            let mut a_next3: Option<TopoShape> = None;
            let mut a_pave2c = BopdsPave::new(0, 0.0);
            if !find_next_vertex(
                the_bound_edge_index,
                &a_first_pave,
                the_ds,
                &mut a_next3,
                &mut a_pave2c,
            ) {
                break;
            }
            a_first_pave = a_pave2c;
            apfirst = a_c2.d0(a_pave2c.parameter());
            a_first_vertex = a_next3;
        }
    }
    // `by pairs non continuously. end`

    // `by pairs continuously. begin` (`cxx:2228-2511`).
    a_pave1 = *the_prev_pave;
    if is_on_u_edge {
        let mut atmp_pave = BopdsPave::new(0, 0.0);
        if !get_pave(the_bound_edge_index, true, the_ds, &mut atmp_pave) {
            return None;
        }
        a_pave1 = atmp_pave;
    }
    p1 = a_c2.d0(a_pave1.parameter());
    a_prev_vertex = the_ds.shape(a_pave1.index()).cloned();

    pit = 0;
    let mut a_next_vertex: Option<TopoShape> = None;
    let mut a_pave2 = BopdsPave::new(0, 0.0);
    while find_next_vertex(
        the_bound_edge_index,
        &a_pave1,
        the_ds,
        &mut a_next_vertex,
        &mut a_pave2,
    ) && pit < nbpave
    {
        let Some(a_first_v) = a_prev_vertex.clone() else {
            return None;
        };
        let Some(a_last_v) = a_next_vertex.clone() else {
            return None;
        };
        let p2 = a_c2.d0(a_pave2.parameter());

        // `cxx:2250-2318`.
        let (bisinside, apbindex) = pair_inside_continuous(
            &a_pave1,
            &a_pave2,
            &p2,
            &a_first_paves,
            &a_last_paves,
            &a_list_of_flags,
            the_ue2_old,
            &a_face_f,
            key,
            the_ds,
        )?;

        if !bisinside {
            let mut a_ordered_list: Vec<TopoShape> = Vec::new();
            if fill_gap(
                &a_first_v,
                &a_last_v,
                &p1,
                &p2,
                &a_face_f,
                &comp_of_sec_edges,
                &mut a_ordered_list,
            ) {
                comp_of_sec_edges = remove_edges(&comp_of_sec_edges, &a_ordered_list);
                list_of_wire_edges.extend(a_ordered_list);
                is_section_found = true;
            } else {
                // `cxx:2333-2353`: fall back to the split edge of the pave pair.
                if let Some(a_e_split) = find_split_edge(
                    the_bound_edge_index,
                    &a_pave1,
                    &a_pave2,
                    the_ds,
                ) {
                    list_of_wire_edges.push(a_e_split);
                }
            }
        } else if apbindex > 0 {
            // `cxx:2357-2362`.
            if let Some(a_list_of_sec) = an_array_of_list_of_sec.get(apbindex) {
                list_of_wire_edges.extend(a_list_of_sec.clone());
            }
        }
        a_pave1 = a_pave2;
        a_prev_vertex = a_next_vertex.clone();
        p1 = p2;
        pit += 1;
    }

    // `cxx:2368-2478`: the closing step against the second `u`-iso edge.
    let mut a_next2: Option<TopoShape> = None;
    let mut a_pave2v = BopdsPave::new(0, 0.0);
    if find_vertex(
        the_ue2_old,
        the_rank,
        the_ds,
        the_hist_map,
        &mut a_next2,
        &mut a_pave2v,
    ) {
        if let (Some(a_first_v), Some(a_last_v)) = (a_prev_vertex.clone(), a_next2.clone()) {
            let Some(a_c3) = reg.edge_pcurve(the_ue2_old, key) else {
                return None;
            };
            let p2 = a_c3.d0(a_pave2v.parameter());

            // `cxx:2377-2426`: the `flag == 0` arm here is unconditional
            // (`cxx:2410-2414`), unlike the in-loop variant above.
            let (bisinside, apbindex) = pair_inside_closing(
                &a_pave1,
                &a_pave2v,
                &a_first_paves,
                &a_last_paves,
                &a_list_of_flags,
            );

            if !bisinside {
                let mut a_ordered_list: Vec<TopoShape> = Vec::new();
                if fill_gap(
                    &a_first_v,
                    &a_last_v,
                    &p1,
                    &p2,
                    &a_face_f,
                    &comp_of_sec_edges,
                    &mut a_ordered_list,
                ) {
                    comp_of_sec_edges = remove_edges(&comp_of_sec_edges, &a_ordered_list);
                    list_of_wire_edges.extend(a_ordered_list);
                    is_section_found = true;
                } else {
                    // `cxx:2442-2468`: the last pave replaces `aPave2`.
                    let mut a_last_pave = BopdsPave::new(0, 0.0);
                    if !get_pave(the_bound_edge_index, false, the_ds, &mut a_last_pave) {
                        return None;
                    }
                    if let Some(a_e_split) =
                        find_split_edge(the_bound_edge_index, &a_pave1, &a_last_pave, the_ds)
                    {
                        list_of_wire_edges.push(a_e_split);
                    }
                }
            } else if apbindex > 0 {
                if let Some(a_list_of_sec) = an_array_of_list_of_sec.get(apbindex) {
                    list_of_wire_edges.extend(a_list_of_sec.clone());
                }
            }
        }
    } else {
        // `cxx:2479-2509`.
        let mut a_last_pave = BopdsPave::new(0, 0.0);
        if !get_pave(the_bound_edge_index, false, the_ds, &mut a_last_pave) {
            return None;
        }
        if let Some(a_e_split) =
            find_split_edge(the_bound_edge_index, &a_pave1, &a_last_pave, the_ds)
        {
            list_of_wire_edges.push(a_e_split);
        }
    }

    Some(FromVEdge {
        list_of_wire_edges,
        is_section_found,
    })
}

/// `FindFromVEdge`'s in-loop "is this pair already covered" test
/// (`cxx:2254-2315`).
///
/// Returns `(bisinside, apbindex)`; `None` is an unreachable null pcurve, which
/// the reference would dereference.
#[allow(clippy::too_many_arguments)]
fn pair_inside_continuous(
    the_pave1: &BopdsPave,
    the_pave2: &BopdsPave,
    p2: &GpPnt2d,
    the_first_paves: &[BopdsPave],
    the_last_paves: &[BopdsPave],
    the_flags: &[i32],
    the_ue2_old: &TopoShape,
    the_face: &TopoShape,
    key: usize,
    the_ds: &BopdsDS,
) -> Option<(bool, usize)> {
    let mut bisinside = false;
    let mut apbindex = 0usize;
    let count = the_first_paves
        .len()
        .min(the_last_paves.len())
        .min(the_flags.len());

    // `cxx:2258-2262`: `apbcounter` starts at 1.
    for apbcounter in 1..=count {
        let a_pit1 = &the_first_paves[apbcounter - 1];
        let a_pit2 = &the_last_paves[apbcounter - 1];
        let a_flag = the_flags[apbcounter - 1];

        // `cxx:2267-2275`.
        let bfin = if the_pave1.is_equal(a_pit1) {
            true
        } else {
            the_pave1.parameter() > a_pit1.parameter()
        };

        // `cxx:2277-2305`.
        let blin = if a_flag != 0 {
            if the_pave2.is_equal(a_pit2) {
                true
            } else {
                the_pave2.parameter() < a_pit2.parameter()
            }
        } else if the_pave2.index() == a_pit2.index() && the_pave2.index() > 0 {
            let reg = GeometryRegistry::global();
            let pc = reg.edge_pcurve(the_ue2_old, key)?;
            let p3 = pc.d0(a_pit2.parameter());
            let a_v = the_ds.shape(the_pave2.index())?;
            let a_tolerance = surface_tolerance(the_face, BRepTool::vertex_tolerance(&Vertex(a_v.clone())));
            pnt2d_distance(p2, &p3) < a_tolerance
        } else {
            false
        };

        if bfin && blin {
            apbindex = apbcounter;
            bisinside = true;
            break;
        }
    }
    Some((bisinside, apbindex))
}

/// `FindFromVEdge`'s closing "is this pair already covered" test
/// (`cxx:2388-2425`): same shape as [`pair_inside_continuous`] except that the
/// `flag == 0` arm is unconditional (`cxx:2410-2414`).
fn pair_inside_closing(
    the_pave1: &BopdsPave,
    the_pave2: &BopdsPave,
    the_first_paves: &[BopdsPave],
    the_last_paves: &[BopdsPave],
    the_flags: &[i32],
) -> (bool, usize) {
    let mut bisinside = false;
    let mut apbindex = 0usize;
    let count = the_first_paves
        .len()
        .min(the_last_paves.len())
        .min(the_flags.len());
    for apbcounter in 1..=count {
        let a_pit1 = &the_first_paves[apbcounter - 1];
        let a_pit2 = &the_last_paves[apbcounter - 1];
        let a_flag = the_flags[apbcounter - 1];

        let bfin = if the_pave1.is_equal(a_pit1) {
            true
        } else {
            the_pave1.parameter() > a_pit1.parameter()
        };
        let blin = if a_flag != 0 {
            if the_pave2.is_equal(a_pit2) {
                true
            } else {
                the_pave2.parameter() < a_pit2.parameter()
            }
        } else {
            true
        };
        if bfin && blin {
            apbindex = apbcounter;
            bisinside = true;
            break;
        }
    }
    (bisinside, apbindex)
}

/// The `aESplit` lookup of `FindFromVEdge` (`cxx:2333-2353`, `cxx:2450-2463`):
/// the split edge of the pave block of the bound edge whose ends are exactly
/// `(the_pave1, the_pave2)`.
fn find_split_edge(
    the_bound_edge_index: usize,
    the_pave1: &BopdsPave,
    the_pave2: &BopdsPave,
    the_ds: &BopdsDS,
) -> Option<TopoShape> {
    for pb in the_ds.pave_blocks(the_bound_edge_index) {
        // `cxx:2338-2339`: `OriginalEdge() == theBoundEdgeIndex` plus both ends.
        if pb.original_edge() != the_bound_edge_index
            || !pb.pave1().is_equal(the_pave1)
            || !pb.pave2().is_equal(the_pave2)
        {
            continue;
        }
        // `cxx:2341-2346`: `aPB1->Edge() > 0`.
        if pb.has_edge() && pb.edge() > 0 {
            if let Some(s) = the_ds.shape(pb.edge()) {
                return Some(s.clone());
            }
        }
    }
    None
}