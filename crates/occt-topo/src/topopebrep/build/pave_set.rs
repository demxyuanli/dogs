//! `TopOpeBRepBuild_PaveSet` — an exploration of a set of vertices to build
//! edges.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_PaveSet.hxx:30-72`,
//! `.cxx:28-489`.
//!
//! It is similar to `LoopSet` from `TopOpeBRepBuild` where `Loop` is `Pave`.
//!
//! Translation notes:
//! - `myVertices` is a list of `handle<Pave>`; the port keeps
//!   `Vec<LoopHandle>` because `Pave` is the
//!   [`LoopKind::Pave`](super::loop_shape::LoopKind::Pave) payload of a
//!   [`Loop`](super::loop_shape::Loop), exactly as OCCT stores one homogeneous
//!   list of base handles.
//! - `NCollection_List::Iterator` is a cursor field.

use occt_core::precision::PCONFUSION;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::shape::TopoShape;
use crate::topo_tools::{edge_vertices, explorer_shapes, is_equal, is_same};
use crate::Edge;

use std::sync::{Arc, RwLock};

use super::loop_set::LoopSetSource;
use super::loop_shape::{Loop, LoopHandle};

/// `class TopOpeBRepBuild_PaveSet : public TopOpeBRepBuild_LoopSet`
/// (`hxx:33`).
#[derive(Debug, Clone)]
pub struct PaveSet {
    /// The inherited `TopOpeBRepBuild_LoopSet::myListOfLoop` (`LoopSet.hxx:46`).
    /// Left unused: every iterator method of the base is overridden
    /// (`PaveSet.hxx:44-50`).
    my_list_of_loop: Vec<LoopHandle>,
    /// `myEdge` (`hxx:57`).
    my_edge: TopoShape,
    /// `myVertices` (`hxx:58`).
    my_vertices: Vec<LoopHandle>,
    /// `myVerticesIt` (`hxx:59`) — a cursor into `myVertices`.
    my_vertices_it: usize,
    /// `myHasEqualParameters` (`hxx:60`).
    my_has_equal_parameters: bool,
    /// `myEqualParameters` (`hxx:61`).
    my_equal_parameters: f64,
    /// `myClosed` (`hxx:62`).
    my_closed: bool,
    /// `myPrepareDone` (`hxx:63`).
    my_prepare_done: bool,
    /// `myRemovePV` (`hxx:64`).
    my_remove_pv: bool,
}

impl PaveSet {
    /// `TopOpeBRepBuild_PaveSet(E)` (`cxx:32-39`).
    pub fn new(e: &TopoShape) -> Self {
        Self {
            my_list_of_loop: Vec::new(),
            my_edge: e.clone(),
            my_vertices: Vec::new(),
            my_vertices_it: 0,
            my_has_equal_parameters: false,
            my_equal_parameters: 0.0,
            my_closed: false,
            my_prepare_done: false,
            my_remove_pv: true,
        }
    }

    /// `RemovePV(B)` (`cxx:43-46`).
    pub fn set_remove_pv(&mut self, b: bool) {
        self.my_remove_pv = b;
    }

    /// `Append(PV)` (`cxx:50-54`).
    pub fn append(&mut self, pv: LoopHandle) {
        self.my_vertices.push(pv);
        self.my_prepare_done = false;
    }

    /// `SortPave(List, SortedList)` (`cxx:62-127`).
    ///
    /// `cxx:64-65`: sorting a list of items by a given function is not
    /// available, so the original is an `n^2` selection of the smallest
    /// `Parameter()` remaining, followed by the `FORWARD`-head reordering of
    /// `cxx:93-124`.
    pub fn sort_pave(list: &[LoopHandle]) -> Vec<LoopHandle> {
        let n_pv = list.len();
        let mut t = vec![false; n_pv + 1];
        let mut sorted_list: Vec<LoopHandle> = Vec::with_capacity(n_pv);

        let mut pv1: Option<LoopHandle> = None;
        let mut i_pv = 0usize;
        for _i in 1..=n_pv {
            let mut parmin = f64::MAX;
            for (itest, pv2) in list.iter().enumerate() {
                let itest = itest + 1;
                if !t[itest] {
                    let par = match pv2.read() {
                        Ok(p) => p.as_pave().map(|q| q.parameter()).unwrap_or(0.0),
                        Err(_) => 0.0,
                    };
                    if par < parmin {
                        parmin = par;
                        pv1 = Some(pv2.clone());
                        i_pv = itest;
                    }
                }
            }
            if let Some(pv) = pv1.clone() {
                sorted_list.push(pv);
            }
            t[i_pv] = true;
        }

        // tete = FORWARD
        // modifier TopOpeBRepBuild_DataStructure::SortOnParameter
        let mut found = false;
        let mut l1: Vec<LoopHandle> = Vec::new();
        let mut l2: Vec<LoopHandle> = Vec::new();
        for pv in sorted_list.iter() {
            let o = match pv.read() {
                Ok(p) => p
                    .as_pave()
                    .map(|q| q.vertex.orientation())
                    .unwrap_or(Orientation::External),
                Err(_) => Orientation::External,
            };
            if !found {
                if o == Orientation::Forward {
                    found = true;
                    l1.push(pv.clone());
                } else {
                    l2.push(pv.clone());
                }
            } else {
                l1.push(pv.clone());
            }
        }

        let mut out = Vec::with_capacity(l1.len() + l2.len());
        out.extend(l1);
        out.extend(l2);
        out
    }

    /// `Prepare()` (`cxx:142-289`).
    ///
    /// Adds the edge vertices to the list of interferences; if an edge vertex
    /// `VE` is already in the list as interference `VI` (`cxx:143-148`):
    /// `VE` is not added, `VI` takes the orientation of `VE` when `VI` is
    /// `INTERNAL`, and `VI` is removed when it is `EXTERNAL` or when `VE` and
    /// `VI` have opposite orientations.
    pub fn prepare(&mut self) {
        if self.my_prepare_done {
            return;
        }

        let is_ed = BRepTool::is_degenerated(&Edge(self.my_edge.clone()));
        let mut edge_vertex_count = 0i32;

        if self.my_remove_pv {
            // jyl + 980217
            let edge_vertices_exp = explorer_shapes(&self.my_edge, ShapeType::Vertex);
            for ve in edge_vertices_exp {
                // VE = edge vertex
                let ve_ori = ve.orientation();
                let ve_bound = ve_ori == Orientation::Forward || ve_ori == Orientation::Reversed;

                let mut edge_vertex_index = 0i32;
                let mut add_ve = true;
                let mut add = false; // ofv
                let mut remove_at: Option<usize> = None;

                for (i, handle) in self.my_vertices.iter().enumerate() {
                    edge_vertex_index += 1; // skip edge vertices inserted at the head
                    if edge_vertex_index <= edge_vertex_count {
                        continue;
                    }

                    // PV = Parametrized vertex, VI = interference vertex
                    let mut pv = match handle.write() {
                        Ok(p) => p,
                        Err(p) => p.into_inner(),
                    };
                    let p = match pv.as_pave_mut() {
                        Some(p) => p,
                        None => continue,
                    };
                    let vi_ori = p.vertex.orientation();
                    let vi_same_ve = is_same(&p.vertex, &ve);
                    let vsd_same_ve = if p.has_same_domain {
                        match p.same_domain.as_ref() {
                            Some(vsd) => is_same(vsd, &ve),
                            None => false,
                        }
                    } else {
                        false
                    };
                    let same_vertex_processing = (vi_same_ve || vsd_same_ve) && !is_ed;

                    if same_vertex_processing {
                        // if (VEbound) { // xpu: 29-05-97
                        if ve_bound || vsd_same_ve {
                            match vi_ori {
                                Orientation::External => {
                                    remove_at = Some(i);
                                }
                                Orientation::Internal => {
                                    p.vertex.set_orientation(ve_ori);
                                }
                                Orientation::Forward | Orientation::Reversed => {
                                    // ofv:
                                    // if (VIori != VEori) myVertices.Remove(it);
                                    if vi_ori != ve_ori {
                                        remove_at = Some(i);
                                        let is_look = fun_is_look(&self.my_edge);
                                        if ve_bound && (vsd_same_ve || vi_same_ve) && is_look {
                                            add = true;
                                        }
                                    }
                                }
                            }
                        }
                        // ofv:
                        // addVE = false;
                        add_ve = add;
                        break;
                    }
                }

                if let Some(i) = remove_at {
                    self.my_vertices.remove(i);
                }
                // if VE not found in the list, add it
                if add_ve {
                    let par_ve = BRepTool::parameter_of_vertex(&self.my_edge, &ve);
                    let new_pv = Loop::new_pave(&ve, par_ve, true);
                    self.my_vertices
                        .insert(0, Arc::new(RwLock::new(new_pv)));
                    edge_vertex_count += 1;
                }
            }
        } // myRemovePV

        let ll = self.my_vertices.len() as i32;

        // if no more interferences vertices, clear the list
        if ll == edge_vertex_count {
            self.my_vertices.clear();
        } else if ll >= 2 {
            // sort the parametrized vertices on Parameter() value.
            let list = self.my_vertices.clone();
            self.my_vertices.clear();
            self.my_vertices = Self::sort_pave(&list);
        }

        self.my_prepare_done = true;
    }

    /// `Edge()` (`cxx:339-342`).
    pub fn edge(&self) -> &TopoShape {
        &self.my_edge
    }

    /// `HasEqualParameters()` (`cxx:346-419`).
    pub fn has_equal_parameters(&mut self) -> bool {
        self.my_has_equal_parameters = false;
        let n = self.my_vertices.len();

        for i in 0..n {
            if self.my_has_equal_parameters {
                break;
            }
            let (v1, p1) = match self.my_vertices[i].read() {
                Ok(p) => match p.as_pave() {
                    Some(q) => (q.vertex.clone(), q.parameter()),
                    None => continue,
                },
                Err(_) => continue,
            };

            for j in 0..n {
                if self.my_has_equal_parameters {
                    break;
                }
                let (v2, p2) = match self.my_vertices[j].read() {
                    Ok(p) => match p.as_pave() {
                        Some(q) => (q.vertex.clone(), q.parameter()),
                        None => continue,
                    },
                    Err(_) => continue,
                };
                if is_equal(&v2, &v1) {
                    continue;
                }
                let d = (p1 - p2).abs();
                if d < PCONFUSION {
                    self.my_has_equal_parameters = true;
                    self.my_equal_parameters = p1;
                }
            }
        }

        if !self.my_has_equal_parameters {
            let rd;
            let mut f = 0.0;
            {
                let e = Edge(self.my_edge.clone());
                let c_my_edge = BRepTool::edge_curve(&e);
                if c_my_edge.is_none() {
                    rd = false;
                } else {
                    f = BRepTool::edge_parameters(&e).0;
                    rd = true;
                }
            }
            if rd {
                for i in 0..n {
                    if self.my_has_equal_parameters {
                        break;
                    }
                    let p1 = match self.my_vertices[i].read() {
                        Ok(p) => match p.as_pave() {
                            Some(q) => q.parameter(),
                            None => continue,
                        },
                        Err(_) => continue,
                    };
                    let d = (p1 - f).abs();
                    if d < PCONFUSION {
                        self.my_has_equal_parameters = true;
                        self.my_equal_parameters = f;
                    }
                }
            }
        }

        self.my_has_equal_parameters
    }

    /// `EqualParameters()` (`cxx:423-433`).
    ///
    /// `cxx:428-430` raises `NYI` when no equal parameter was found; the port
    /// answers 0 like the `return 0.; // windowsNT` fallback of `cxx:432`.
    pub fn equal_parameters(&self) -> f64 {
        if self.my_has_equal_parameters {
            self.my_equal_parameters
        } else {
            0.0
        }
    }

    /// `ClosedVertices()` (`cxx:437-466`).
    pub fn closed_vertices(&mut self) -> bool {
        if self.my_vertices.is_empty() {
            return false;
        }

        let mut vmin: Option<TopoShape> = None;
        let mut vmax: Option<TopoShape> = None;
        let mut parmin = f64::MAX;
        let mut parmax = -f64::MAX;
        for handle in self.my_vertices.iter() {
            let (v, par) = match handle.read() {
                Ok(p) => match p.as_pave() {
                    Some(q) => (q.vertex.clone(), q.parameter()),
                    None => continue,
                },
                Err(_) => continue,
            };
            if par > parmax {
                vmax = Some(v.clone());
                parmax = par;
            }
            if par < parmin {
                vmin = Some(v.clone());
                parmin = par;
            }
        }

        self.my_closed = match (vmin.as_ref(), vmax.as_ref()) {
            (Some(a), Some(b)) => is_same(a, b),
            _ => false,
        };
        self.my_closed
    }
}

impl LoopSetSource for PaveSet {
    /// `ChangeListOfLoop()` of the inherited `LoopSet` — never consulted, since
    /// every iterator method is overridden below.
    fn change_list_of_loop(&mut self) -> &mut Vec<LoopHandle> {
        &mut self.my_list_of_loop
    }

    /// `InitLoop()` (`cxx:293-300`).
    fn init_loop(&mut self) {
        if !self.my_prepare_done {
            self.prepare();
        }
        self.my_vertices_it = 0;
    }

    /// `MoreLoop()` (`cxx:304-308`).
    fn more_loop(&self) -> bool {
        self.my_vertices_it < self.my_vertices.len()
    }

    /// `NextLoop()` (`cxx:312-315`).
    fn next_loop(&mut self) {
        self.my_vertices_it += 1;
    }

    /// `Loop()` (`cxx:319-322`).
    fn loop_(&self) -> Option<LoopHandle> {
        self.my_vertices.get(self.my_vertices_it).cloned()
    }
}

/// `static bool FUN_islook(const TopoDS_Edge& e)` (`cxx:129-138`): the edge is
/// "looking" when its two vertices are further apart than `1.e-8`.
///
/// `TopExp::Vertices` (`cxx:131`) yields optional handles in the port, so an
/// edge without both vertices answers false where OCCT would dereference a
/// null handle.
fn fun_is_look(e: &TopoShape) -> bool {
    let (v1, v2) = edge_vertices(&Edge(e.clone()));
    let (v1, v2) = match (v1, v2) {
        (Some(a), Some(b)) => (a, b),
        _ => return false,
    };
    let p1 = BRepTool::vertex_point(&v1);
    let p2 = BRepTool::vertex_point(&v2);
    let dp1p2 = p1.distance(&p2);
    dp1p2.abs() > 1.0e-8
}
