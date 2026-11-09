//! `Update` / `StoreShell` / `Normal`. Source: `HLRBRep_PolyAlgo::Update`
//! (`HLRBRep_PolyAlgo.cxx:149-201`), `StoreShell` (`cxx:289-623`) and `Normal`
//! (`cxx:626-715`).
//!
//! This is the batch that turns the per-face triangulations into
//! `HLRAlgo_PolyShellData`: it initialises `TMat`/`TLoc` and `TIMa`/`TILo`
//! (`cxx:302-341`), decides whether the shell is closed (`cxx:347-375`),
//! allocates one `HLRAlgo_PolyInternalData` per face (`cxx:447-466`), fills the
//! node records with `TTMultiply`-projected points (`cxx:485-497`), computes the
//! node normals through `Normal`, orients the triangles, assembles the bi-points
//! and finally hands everything to the outline passes and to
//! `UpdateEdgesBiPoints` / `UpdatePolyData`.
//!
//! Surface access seams:
//! * `BRep_Tool::Surface(aFace)` (`cxx:454`) is
//!   [`crate::tgeometry::GeometryRegistry::face_surface`], the *local* surface,
//!   which is also what `myGSurf` holds (`cxx:3009`);
//! * `myBSurf.Initialize(aFace, false)` (`cxx:500`) is the face surface plus the
//!   face location: `BRepAdaptor_Surface` derives from
//!   `GeomAdaptor_TransformedSurface` (`BRepAdaptor_Surface.hxx:50`) and its
//!   `EvalD1`/`EvalD2` transform the local derivatives by that location
//!   (`GeomAdaptor_TransformedSurface.cxx:248-297`), so [`HlrBRepPolyAlgo::normal`]
//!   applies [`GpVec::transform`] to every derivative vector;
//! * `GeomAdaptor_Surface::GetType() == GeomAbs_Plane` after unwrapping a
//!   `Geom_RectangularTrimmedSurface` (`cxx:455-466`) is
//!   [`crate::brep_surface::is_plane_type`].

use std::collections::HashSet;

use occt_core::cslib::normal::{normal_d1, normal_d2, CSLibDerivativeStatus};
use occt_core::gp::{GpTrsf, GpXyz};
use occt_core::precision::{ANGULAR, CONFUSION};
use occt_geom::Surface;

use super::{
    has_triangulation, nmask, node_adjacency::orient_triangle_in, HlrBRepPolyAlgo,
    PolyInternalDataArray, ShapeIndexedMap,
};
use crate::abs::{Orientation, ShapeType};
use crate::brep_surface::is_plane_type;
use crate::brep_tool::BRepTool;
use crate::hlr_algo::{HlrAlgoPolyData, HlrAlgoPolyInternalData, HlrAlgoPolyShellData};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{
    explorer_shapes, explorer_shapes_avoiding, map_shapes_and_ancestors, shape_map_key, ShapeMapKey,
};

impl HlrBRepPolyAlgo {
    /// `Update()` (`cxx:149-201`).
    pub fn update(&mut self) {
        self.algo.clear();
        self.e_map.clear();
        self.f_map.clear();

        // `cxx:156-160`.
        let Some(a_shape) = self.make_shape() else {
            return;
        };

        // `cxx:162-167`.
        let (a_nb_shell, is_isoled_f, is_isoled_e) = self.init_shape(&a_shape);
        if a_nb_shell <= 0 {
            return;
        }

        // `cxx:169-170`: `TopExp::MapShapes(aShape, TopAbs_EDGE/FACE, ...)`, the
        // global indices the bi-points and the poly data carry.
        for e in explorer_shapes(&a_shape, ShapeType::Edge) {
            self.e_map.add(&e);
        }
        for f in explorer_shapes(&a_shape, ShapeType::Face) {
            self.f_map.add(&f);
        }
        let a_nb_edge = self.e_map.extent();
        let a_nb_face = self.f_map.extent();

        // `cxx:171-176`: `anES(0, aNbEdge)`, `aPD(0, aNbFace)`,
        // `aPID(0, aNbFace)`, slot 0 unused.
        let mut an_es: Vec<i32> = vec![0; a_nb_edge + 1];
        let mut a_pd: Vec<HlrAlgoPolyData> = (0..=a_nb_face).map(|_| HlrAlgoPolyData::new()).collect();
        let mut a_pid: PolyInternalDataArray = vec![None; a_nb_face + 1];
        let mut a_shape_map1: HashSet<ShapeMapKey> = HashSet::new();
        let mut a_shape_map2: HashSet<ShapeMapKey> = HashSet::new();

        // `cxx:177`: `aShell` is `myAlgo->ChangePolyShell()`. Taking the vector
        // out keeps the `StoreShell(..., theShell, ...)` signature without
        // borrowing `myAlgo` for the whole loop; it is put back before
        // `myAlgo->Update()`.
        let mut a_shell = std::mem::take(self.algo.change_poly_shell());
        let mut i_shell = 0usize;

        // `cxx:178-190`.
        for a_shell_iter in explorer_shapes(&a_shape, ShapeType::Shell) {
            self.store_shell(
                &a_shell_iter,
                &mut i_shell,
                &mut a_shell,
                false,
                false,
                &mut an_es,
                &mut a_pd,
                &mut a_pid,
                &mut a_shape_map1,
                &mut a_shape_map2,
            );
        }
        // `cxx:191-194`.
        if is_isoled_f {
            self.store_shell(
                &a_shape,
                &mut i_shell,
                &mut a_shell,
                true,
                false,
                &mut an_es,
                &mut a_pd,
                &mut a_pid,
                &mut a_shape_map1,
                &mut a_shape_map2,
            );
        }
        // `cxx:195-198`.
        if is_isoled_e {
            self.store_shell(
                &a_shape,
                &mut i_shell,
                &mut a_shell,
                false,
                true,
                &mut an_es,
                &mut a_pd,
                &mut a_pid,
                &mut a_shape_map1,
                &mut a_shape_map2,
            );
        }

        *self.algo.change_poly_shell() = a_shell;

        // `cxx:199`.
        self.algo.update();
    }

    /// `StoreShell(...)` (`cxx:289-623`): fill one `HLRAlgo_PolyShellData`.
    ///
    /// `the_es` (`anES`) is the per-edge shell index OCCT fills at `cxx:579` and
    /// `cxx:615`; nothing in `Update` reads it back, exactly as in the source.
    #[allow(clippy::too_many_arguments)]
    pub fn store_shell(
        &mut self,
        the_shape: &TopoShape,
        the_i_shell: &mut usize,
        the_shell: &mut Vec<HlrAlgoPolyShellData>,
        the_isoled_f: bool,
        the_isoled_e: bool,
        the_es: &mut [i32],
        the_pd: &mut [HlrAlgoPolyData],
        the_pid: &mut PolyInternalDataArray,
        the_shape_map1: &mut HashSet<ShapeMapKey>,
        the_shape_map2: &mut HashSet<ShapeMapKey>,
    ) {
        let mut a_nb_face_shell = 0i32;
        let mut is_closed = false;

        // `cxx:302-341`: `TMat`/`TLoc` and `TIMa`/`TILo` from the projector.
        self.read_projection_matrices();

        // `cxx:342-394`: `aFaceExp` is only initialised and walked when the
        // shell is not an isolated-edge shell; OCCT leaves the explorer
        // untouched in that arm, so the port keeps it local to this block.
        if !the_isoled_e {
            let a_face_exp: Vec<TopoShape> = if !the_isoled_f {
                // `cxx:346`: `TopoDS_Shape::Closed()`.
                is_closed = the_shape.is_closed();
                if !is_closed {
                    // `cxx:348-375`: an edge used once forward and never
                    // reversed leaves the shell open.
                    let a_edges = explorer_shapes(the_shape, ShapeType::Edge);
                    let mut an_em = ShapeIndexedMap::new();
                    for e in &a_edges {
                        an_em.add(e);
                    }
                    let a_nb_edge = an_em.extent();
                    // `NCollection_Array1<int> aFlagArray(1, max(aNbEdge, 1))`.
                    let mut a_flag_array: Vec<i32> = vec![0; a_nb_edge + 1];
                    for e in &a_edges {
                        let an_edge = Edge(e.clone());
                        let an_edge_index = an_em.find_index(e);
                        let an_edge_orient = e.orientation();
                        if !BRepTool::is_degenerated(&an_edge) {
                            if an_edge_orient == Orientation::Forward {
                                a_flag_array[an_edge_index] += 1;
                            } else if an_edge_orient == Orientation::Reversed {
                                a_flag_array[an_edge_index] -= 1;
                            }
                        }
                    }
                    is_closed = true;
                    let mut an_edge_iter = 1usize;
                    while an_edge_iter <= a_nb_edge && is_closed {
                        is_closed = a_flag_array[an_edge_iter] == 0;
                        an_edge_iter += 1;
                    }
                }
                // `cxx:377`.
                explorer_shapes(the_shape, ShapeType::Face)
            } else {
                // `cxx:382`.
                explorer_shapes_avoiding(the_shape, ShapeType::Face, ShapeType::Shell)
            };

            // `cxx:385-394`.
            for a_face in &a_face_exp {
                if has_triangulation(a_face) && the_shape_map1.insert(shape_map_key(a_face)) {
                    a_nb_face_shell += 1;
                }
            }
        }

        // `cxx:396-401`.
        if a_nb_face_shell > 0 || the_isoled_e {
            *the_i_shell += 1;
            the_shell[*the_i_shell - 1] = HlrAlgoPolyShellData::new(a_nb_face_shell as usize);
        }

        if a_nb_face_shell > 0 {
            // `cxx:408-411`: `aFaceExp.Init(...)` for the same shape as above.
            let a_face_iter: Vec<TopoShape> = if !the_isoled_f {
                explorer_shapes(the_shape, ShapeType::Face)
            } else {
                explorer_shapes_avoiding(the_shape, ShapeType::Face, ShapeType::Shell)
            };
            // `aPsd->PolyData().SetValue(iFace, thePD.Value(aFaceIndex))`
            // (`cxx:448`) aliases the shell slot with the face-indexed array; the
            // port records the pair and copies the finished poly data below.
            let mut shell_faces: Vec<(usize, usize)> = Vec::with_capacity(a_nb_face_shell as usize);
            let mut i_face = 0usize;

            for a_face in &a_face_iter {
                let face = Face(a_face.clone());
                // `cxx:417-420`.
                let Some(a_tr) = BRepTool::face_triangulation(&face) else {
                    continue;
                };
                // `cxx:422`.
                if !the_shape_map2.insert(shape_map_key(a_face)) {
                    continue;
                }
                // `cxx:424-425`.
                let a_face_index = self.f_map.find_index(a_face);
                if a_face_index == 0 {
                    // `myFMap.FindIndex` cannot miss a face of the compound
                    // `Update` mapped; OCCT would write slot 0 of a 0-based
                    // array here, i.e. raise `Standard_OutOfRange`.
                    continue;
                }
                i_face += 1;
                let is_reversed = a_face.orientation().is_reversed();

                // `cxx:427-446`: `TTMa`/`TTLo` from `aLoc * aProjTrsf`.
                self.read_face_matrices(&a_face.location().transformation());

                let a_nb_nodes = a_tr.nb_nodes() as i32;
                let a_nb_tris = a_tr.nb_triangles() as i32;

                // `cxx:447-449`.
                the_pd[a_face_index] = HlrAlgoPolyData::new();
                shell_faces.push((i_face, a_face_index));

                // `cxx:450`.
                let mut a_pid = HlrAlgoPolyInternalData::new(a_nb_nodes, a_nb_tris);

                // `cxx:452-466`.
                let a_surf = BRepTool::face_surface(&face);
                let is_planar = a_surf
                    .as_deref()
                    .map(is_plane_type)
                    .unwrap_or(false);
                a_pid.set_planar(is_planar);

                // `cxx:470-483`: `HLRAlgo_TriangleData::Flags = 0` and the node
                // triple, reversed when the face is.
                for a_tri_iter in 1..=a_nb_tris {
                    let a_poly_tri = a_tr.triangles[(a_tri_iter - 1) as usize];
                    let a_tri_data = &mut a_pid.change_t_data()[a_tri_iter as usize];
                    a_tri_data.flags = 0;
                    if is_reversed {
                        a_tri_data.node3 = a_poly_tri.n0 as i32 + 1;
                        a_tri_data.node2 = a_poly_tri.n1 as i32 + 1;
                        a_tri_data.node1 = a_poly_tri.n2 as i32 + 1;
                    } else {
                        a_tri_data.node1 = a_poly_tri.n0 as i32 + 1;
                        a_tri_data.node2 = a_poly_tri.n1 as i32 + 1;
                        a_tri_data.node3 = a_poly_tri.n2 as i32 + 1;
                    }
                }

                // `cxx:485-497`: the point is the triangulation node mapped
                // through `TTMultiply`.
                for a_node_iter in 1..=a_nb_nodes {
                    let a_pnt = a_tr.nodes[(a_node_iter - 1) as usize];
                    let node = &mut a_pid.change_pi_nod()[a_node_iter as usize];
                    node.change_indices().nd_sg = 0;
                    node.change_indices().flag = 0;
                    let mut a_point = GpXyz::new(a_pnt.x(), a_pnt.y(), a_pnt.z());
                    self.ttmultiply_xyz(&mut a_point, false);
                    node.change_data().point = a_point;
                }

                // `cxx:498`.
                a_pid.update_links();

                // `cxx:499-525`.
                if a_tr.uv_nodes.is_some() {
                    let has_surf = a_surf.is_some();
                    let a_loc = a_face.location().transformation();
                    for a_node_iter in 1..=a_nb_nodes {
                        let idx = (a_node_iter - 1) as usize;
                        {
                            let node = &mut a_pid.change_pi_nod()[a_node_iter as usize];
                            if let Some(uv) = a_tr.uv_nodes.as_ref() {
                                node.change_data().uv =
                                    occt_core::gp::GpXY::new(uv[idx].0, uv[idx].1);
                            }
                            if let Some(normals) = a_tr.normals.as_ref() {
                                let n = normals[idx];
                                node.change_data().normal = GpXyz::new(n.x(), n.y(), n.z());
                            }
                        }
                        let has_normals = a_tr.normals.is_some();
                        let has_uv = a_tr.uv_nodes.is_some();
                        let a_ok = (has_normals || (has_surf && has_uv))
                            && self.normal(
                                &mut a_pid,
                                a_node_iter,
                                a_surf.as_deref(),
                                &a_loc,
                                false,
                            );
                        let node = &mut a_pid.change_pi_nod()[a_node_iter as usize];
                        if a_ok {
                            node.change_indices().flag |= nmask::NORM;
                        } else {
                            node.change_indices().flag &= !nmask::NORM;
                            node.change_data().scal = 0.0;
                        }
                    }
                }

                // `cxx:536-551`: `OrientTriangle` over every triangle of the
                // face (`cxx:3699-3886`).
                for a_tri_iter in 1..=a_nb_tris {
                    orient_triangle_in(&mut a_pid, a_tri_iter, &self.proj);
                }

                // `cxx:450`: `thePID.SetValue(aFaceIndex, ...)`.
                the_pid[a_face_index] = Some(a_pid);
            }

            // `cxx:563-568`.
            let a_nb_faces = self.f_map.extent();
            let ef = map_shapes_and_ancestors(the_shape, ShapeType::Edge, ShapeType::Face);

            // `cxx:569-593`: every edge of every face that carries internal
            // data contributes its bi-points once (`theShapeMap1`).
            let a_shell_idx = *the_i_shell - 1;
            for a_face_iter in 1..=a_nb_faces {
                if the_pid[a_face_iter].is_none() {
                    continue;
                }
                let Some(a_face) = self.f_map.find_key(a_face_iter).cloned() else {
                    continue;
                };
                for an_edge in explorer_shapes(&a_face, ShapeType::Edge) {
                    if !the_shape_map1.insert(shape_map_key(&an_edge)) {
                        continue;
                    }
                    let an_edge_index = self.e_map.find_index(&an_edge);
                    if an_edge_index == 0 || an_edge_index >= the_es.len() {
                        continue;
                    }
                    // `cxx:579`.
                    the_es[an_edge_index] = *the_i_shell as i32;
                    // `cxx:580-590`: `EF.FindIndex(anEdge) > 0` selects the
                    // ancestor list, else the empty list and no connexity.
                    let a_ls: &[TopoShape] = match ef.get(&shape_map_key(&an_edge)) {
                        Some(list) => list.as_slice(),
                        None => &[],
                    };
                    let a_is_connex = ef.contains_key(&shape_map_key(&an_edge));
                    let a_rg = if a_ls.len() == 2 {
                        // `cxx:1080`: `BRep_Tool::Continuity(theEdge, aF1, aF2)`.
                        crate::brep_lib_continuity::stored_continuity(
                            &Edge(an_edge.clone()),
                            &Face(a_ls[0].clone()),
                            &Face(a_ls[1].clone()),
                        )
                    } else {
                        0
                    };
                    let a_list = the_shell[a_shell_idx].change_edges();
                    self.init_bi_points_with_connexity(
                        an_edge_index,
                        &an_edge,
                        a_list,
                        the_pid,
                        a_ls,
                        a_is_connex,
                        a_rg,
                    );
                }
            }

            // `cxx:594-599`.
            {
                let a_list = the_shell[a_shell_idx].change_edges();
                self.insert_on_out_line(the_pid);
                self.check_fr_back_triangles(a_list, the_pid);
                self.update_out_lines(a_list, the_pid);
                self.update_edges_bi_points(a_list, the_pid, is_closed);
                self.update_poly_data(the_pd, the_pid, is_closed);
            }

            // `cxx:448`: publish the finished poly data into the shell slot the
            // face was stored at.
            for (i_face, a_face_index) in shell_faces {
                the_shell[a_shell_idx].change_poly_data()[i_face - 1] =
                    the_pd[a_face_index].clone();
            }

            // `cxx:601-604`: `thePID.ChangeValue(aFaceIter).Nullify()`.
            for a_face_iter in 1..=a_nb_faces {
                if a_face_iter < the_pid.len() {
                    the_pid[a_face_iter] = None;
                }
            }
        } else if the_isoled_e {
            // `cxx:606-621`: isolated edges, no face data at all.
            let a_shell_idx = *the_i_shell - 1;
            let a_list = the_shell[a_shell_idx].change_edges();
            for an_edge in explorer_shapes_avoiding(the_shape, ShapeType::Edge, ShapeType::Face) {
                if !the_shape_map1.insert(shape_map_key(&an_edge)) {
                    continue;
                }
                let an_edge_index = self.e_map.find_index(&an_edge);
                if an_edge_index == 0 || an_edge_index >= the_es.len() {
                    continue;
                }
                // `cxx:615`.
                the_es[an_edge_index] = *the_i_shell as i32;
                self.init_bi_points_with_connexity(
                    an_edge_index,
                    &an_edge,
                    a_list,
                    the_pid,
                    &[],
                    false,
                    0,
                );
            }
        }
    }

    /// `Normal(theNodeIndex, theNodIndices, theNod1RValues, theTriData,
    /// thePISeg, thePINod, theToOrient)` (`cxx:626-715`).
    ///
    /// The two node records arrive as the `the_data` owner plus the node index
    /// because `theTriData`/`thePISeg`/`thePINod` are three views of that same
    /// `HLRAlgo_PolyInternalData` in OCCT. `the_loc` is the face location folded
    /// into `myBSurf` (`BRepAdaptor_Surface.cxx:63-80`).
    ///
    /// Returns `false` exactly where OCCT returns `false` (`CSLib::Normal` could
    /// not produce a direction), and `true` otherwise; the caller owns the
    /// `NMsk_Norm` bit.
    pub fn normal(
        &self,
        the_data: &mut HlrAlgoPolyInternalData,
        the_node_index: i32,
        the_surf: Option<&dyn Surface>,
        the_loc: &GpTrsf,
        the_to_orient: bool,
    ) -> bool {
        let mut a_nod = *the_data.pi_nod()[the_node_index as usize].data();
        let mut a_flag = the_data.pi_nod()[the_node_index as usize].indices().flag;

        // `cxx:635-655`.
        if a_nod.normal.square_modulus() < CONFUSION {
            // UNPORTED: `BRepAdaptor_Surface::Initialize` returns before
            // `mySurf.Load` when the face has no surface
            // (`BRepAdaptor_Surface.cxx:66-69`), so OCCT would evaluate the
            // m_B Surf of the previous face here. The port has no such slot; this
            // arm answers `false`, the same value OCCT returns when
            // `CSLib::Normal` fails (`cxx:648-651`).
            let Some(a_surf) = the_surf else {
                return false;
            };
            let (a_u, a_v) = (a_nod.uv.x(), a_nod.uv.y());
            let (_, mut a_d1u, mut a_d1v) = a_surf.d1(a_u, a_v);
            a_d1u.transform(the_loc);
            a_d1v.transform(the_loc);
            let (a_status, a_dir) = normal_d1(&a_d1u, &a_d1v, ANGULAR);
            let a_norm = if a_status == CSLibDerivativeStatus::Done {
                match a_dir {
                    Some(d) => d,
                    None => return false,
                }
            } else {
                // `cxx:643-651`: the D2 fallback.
                let (_, mut a_d1u, mut a_d1v, mut a_d2u, mut a_d2v, mut a_d2uv) =
                    a_surf.d2(a_u, a_v);
                a_d1u.transform(the_loc);
                a_d1v.transform(the_loc);
                a_d2u.transform(the_loc);
                a_d2v.transform(the_loc);
                a_d2uv.transform(the_loc);
                let (is_ok, _, a_dir) =
                    normal_d2(&a_d1u, &a_d1v, &a_d2u, &a_d2v, &a_d2uv, ANGULAR);
                if !is_ok {
                    return false;
                }
                match a_dir {
                    Some(d) => d,
                    None => return false,
                }
            };
            a_nod.normal = *a_norm.xyz();
        }

        // `cxx:657`: `TMultiply(normal, myProj.Perspective())`.
        let a_persp = self.proj.perspective();
        let mut a_x = a_nod.normal.x();
        let mut a_y = a_nod.normal.y();
        let mut a_z = a_nod.normal.z();
        self.tmultiply(&mut a_x, &mut a_y, &mut a_z, a_persp);
        a_nod.normal = GpXyz::new(a_x, a_y, a_z);

        // `cxx:659-692`.
        if let Some(an_aver_norm) = Self::average_normal(the_data, the_node_index) {
            if a_nod.normal.dot(&an_aver_norm) < 0.0 {
                a_nod.normal.reverse();
            }
            let mut an_eye_dir = GpXyz::new(0.0, 0.0, -1.0);
            if a_persp {
                // `myProj.Focus()` raises when the projection is not
                // perspective; this arm is inside `if (Perspective())`.
                let a_focus = self.proj.focus().unwrap_or(0.0);
                an_eye_dir = GpXyz::new(a_nod.point.x(), a_nod.point.y(), a_nod.point.z() - a_focus);
                let an_eye_mod = an_eye_dir.modulus();
                if an_eye_mod > 0.0 {
                    an_eye_dir = an_eye_dir.multiplied(1.0 / an_eye_mod);
                }
            }
            a_nod.scal = a_nod.normal.dot(&an_eye_dir);
        } else {
            a_nod.scal = 0.0;
            a_nod.normal = GpXyz::new(1.0, 0.0, 0.0);
        }

        // `cxx:693-708`.
        if a_nod.scal > 0.0 {
            if a_nod.scal < self.tol_angular {
                a_nod.scal = 0.0;
                a_flag |= nmask::OUT_L;
            }
        } else if -a_nod.scal < self.tol_angular {
            a_nod.scal = 0.0;
            a_flag |= nmask::OUT_L;
        }

        {
            let node = &mut the_data.change_pi_nod()[the_node_index as usize];
            *node.change_data() = a_nod;
            node.change_indices().flag = a_flag;
        }

        // `cxx:710-713`.
        if the_to_orient {
            self.update_around_node(the_data, the_node_index);
        }
        true
    }
}
