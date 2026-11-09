//! Set of triangles of one face, and the hiding test against them. Source:
//! `HLRAlgo_PolyData`
//! (`ModelingAlgorithms/TKHLR/HLRAlgo/HLRAlgo_PolyData.hxx:35-110`,
//! `HLRAlgo_PolyData.lxx:26-66`, `HLRAlgo_PolyData.cxx:30-1127`).
//!
//! Storage convention: the OCCT arrays are `NCollection_Array1` whose usable
//! range starts at `1`; the port stores plain `Vec`s indexed from `0`, so a
//! record that carries an OCCT 1-based node index (a `HLRAlgo_TriangleData`
//! node field) is read at `nodes[index - 1]`.

use occt_core::bnd::BndBox;
use occt_core::gp::{GpPnt, GpXY, GpXyz};

use super::bi_point::{BiPointIndices, BiPointPoints};
use super::edge_status::HlrAlgoEdgeStatus;
use super::poly_records::{poly_mask, HlrAlgoPolyHidingData, HlrAlgoTriangleData};

/// `HLRAlgo_PolyData::FaceIndices` (`hxx:40-51`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PolyFaceIndices {
    pub index: i32,
    pub min: i32,
    pub max: i32,
}

/// `HLRAlgo_PolyData::Triangle` (`hxx:53-57`) — the per-edge working triangle
/// handed down from `HLRAlgo_PolyAlgo`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct PolyTriangle {
    pub v1: GpXY,
    pub v2: GpXY,
    pub v3: GpXY,
    pub param: f64,
    pub tol_param: f64,
    pub tol_ang: f64,
    pub tolerance: f64,
}

/// The two 31-bit depth-field masks used by the bounding-box rejection tests
/// (`PolyData.cxx:85-86`, `:102-103`; `PolyAlgo.cxx:533-534`).
pub const POLY_BBOX_MASK_A: i32 = 0x8010_0200u32 as i32;
pub const POLY_BBOX_MASK_B: i32 = 0x8010_0000u32 as i32;

/// Data structure of a set of triangles. Source: `HLRAlgo_PolyData`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct HlrAlgoPolyData {
    face_indices: PolyFaceIndices,
    /// `myHNodes`; empty stands for the null handle.
    nodes: Vec<GpXyz>,
    /// `myHTData`; empty stands for the null handle.
    t_data: Vec<HlrAlgoTriangleData>,
    /// `myHPHDat`; `None` stands for the null handle, which is exactly what
    /// `Hiding()` tests (`lxx:60-64`).
    ph_dat: Option<Vec<HlrAlgoPolyHidingData>>,
}

impl HlrAlgoPolyData {
    /// `HLRAlgo_PolyData()` (`cxx:32`) — all handles null.
    pub fn new() -> Self {
        Self::default()
    }

    /// `HNodes(HNodes)` (`cxx:36-39`).
    pub fn set_h_nodes(&mut self, nodes: Vec<GpXyz>) {
        self.nodes = nodes;
    }

    /// `HTData(HTData)` (`cxx:41-45`).
    pub fn set_h_t_data(&mut self, t_data: Vec<HlrAlgoTriangleData>) {
        self.t_data = t_data;
    }

    /// `HPHDat(HPHDat)` (`cxx:47-53`).
    pub fn set_h_ph_dat(&mut self, ph_dat: Option<Vec<HlrAlgoPolyHidingData>>) {
        self.ph_dat = ph_dat;
    }

    /// `FaceIndex(I)` (`lxx:26-29`).
    pub fn set_face_index(&mut self, i: i32) {
        self.face_indices.index = i;
    }

    /// `FaceIndex()` (`lxx:32-35`).
    pub fn face_index(&self) -> i32 {
        self.face_indices.index
    }

    /// `Nodes()` (`lxx:38-42`); `nodes()[i]` is `Nodes(i)` (1-based).
    pub fn nodes(&self) -> &[GpXyz] {
        &self.nodes
    }

    /// `TData()` (`lxx:45-49`); `t_data()[i]` is `TData(i)` (1-based).
    pub fn t_data(&self) -> &[HlrAlgoTriangleData] {
        &self.t_data
    }

    /// `TData()` mutable form (the OCCT `ChangeArray1()` accessor).
    pub fn t_data_mut(&mut self) -> &mut [HlrAlgoTriangleData] {
        &mut self.t_data
    }

    /// `Nodes()` mutable form (the OCCT `ChangeArray1()` accessor).
    pub fn nodes_mut(&mut self) -> &mut [GpXyz] {
        &mut self.nodes
    }

    /// `PHDat()` mutable form (the OCCT `ChangeArray1()` accessor).
    pub fn ph_dat_mut(&mut self) -> Option<&mut Vec<HlrAlgoPolyHidingData>> {
        self.ph_dat.as_mut()
    }

    /// `PHDat()` (`lxx:52-56`); `None` when the handle is null.
    pub fn ph_dat(&self) -> Option<&[HlrAlgoPolyHidingData]> {
        self.ph_dat.as_deref()
    }

    /// `Hiding()` (`lxx:60-64`): `!myHPHDat.IsNull()`.
    pub fn hiding(&self) -> bool {
        self.ph_dat.is_some()
    }

    /// `Indices()` (`hxx:89`).
    pub fn indices(&self) -> &PolyFaceIndices {
        &self.face_indices
    }

    /// `Indices()` mutable form (the OCCT non-const accessor).
    pub fn change_indices(&mut self) -> &mut PolyFaceIndices {
        &mut self.face_indices
    }

    /// `UpdateGlobalMinMax(theBox)` (`cxx:56-77`): add the three nodes of
    /// every triangle flagged `FMskHiding`.
    pub fn update_global_min_max(&self, the_box: &mut BndBox) {
        let nbt = self.t_data.len();
        for i in 1..=nbt {
            let td = &self.t_data[i - 1];
            if td.flags & poly_mask::FMSK_HIDING != 0 {
                let p1 = &self.nodes[td.node1 as usize - 1];
                let p2 = &self.nodes[td.node2 as usize - 1];
                let p3 = &self.nodes[td.node3 as usize - 1];
                the_box.add_point(&GpPnt::new(p1.x(), p1.y(), p1.z()));
                the_box.add_point(&GpPnt::new(p2.x(), p2.y(), p2.z()));
                the_box.add_point(&GpPnt::new(p3.x(), p3.y(), p3.z()));
            }
        }
    }

    /// `HideByPolyData(thePoints, theTriangle, theIndices, HidingShell,
    /// status)` (`cxx:80-217`). `theIndices` is declared non-const in OCCT
    /// (`hxx:84`) but never written.
    pub fn hide_by_poly_data(
        &self,
        the_points: &BiPointPoints,
        the_triangle: &mut PolyTriangle,
        the_indices: &BiPointIndices,
        hiding_shell: bool,
        status: &mut HlrAlgoEdgeStatus,
    ) {
        if ((self.face_indices.max - the_indices.min_seg) & POLY_BBOX_MASK_A) != 0
            || ((the_indices.max_seg - self.face_indices.min) & POLY_BBOX_MASK_B) != 0
        {
            return;
        }
        let Some(ph_dat) = self.ph_dat.as_deref() else {
            return;
        };
        let h2 = ph_dat.len();
        for h in 1..=h2 {
            let ph = &ph_dat[h - 1];
            let a_triangle_indices = ph.indices();
            if ((a_triangle_indices.max - the_indices.min_seg) & POLY_BBOX_MASK_A) != 0
                || ((the_indices.max_seg - a_triangle_indices.min) & POLY_BBOX_MASK_B) != 0
            {
                continue;
            }
            let a_triangle = &self.t_data[a_triangle_indices.index as usize - 1];
            let mut not_connex = true;
            if hiding_shell {
                if self.face_indices.index == the_indices.face_conex1 {
                    if the_indices.face1_pt1 == a_triangle.node1 {
                        not_connex = the_indices.face1_pt2 != a_triangle.node2
                            && the_indices.face1_pt2 != a_triangle.node3;
                    } else if the_indices.face1_pt1 == a_triangle.node2 {
                        not_connex = the_indices.face1_pt2 != a_triangle.node3
                            && the_indices.face1_pt2 != a_triangle.node1;
                    } else if the_indices.face1_pt1 == a_triangle.node3 {
                        not_connex = the_indices.face1_pt2 != a_triangle.node1
                            && the_indices.face1_pt2 != a_triangle.node2;
                    }
                } else if self.face_indices.index == the_indices.face_conex2 {
                    if the_indices.face2_pt1 == a_triangle.node1 {
                        not_connex = the_indices.face2_pt2 != a_triangle.node2
                            && the_indices.face2_pt2 != a_triangle.node3;
                    } else if the_indices.face2_pt1 == a_triangle.node2 {
                        not_connex = the_indices.face2_pt2 != a_triangle.node3
                            && the_indices.face2_pt2 != a_triangle.node1;
                    } else if the_indices.face2_pt1 == a_triangle.node3 {
                        not_connex = the_indices.face2_pt2 != a_triangle.node1
                            && the_indices.face2_pt2 != a_triangle.node2;
                    }
                }
            }
            if !not_connex {
                continue;
            }
            let a_plane = ph.plane();
            let d1 = a_plane.normal.dot(&the_points.pntp1) - a_plane.d;
            let d2 = a_plane.normal.dot(&the_points.pntp2) - a_plane.d;
            let p1 = self.nodes[a_triangle.node1 as usize - 1];
            let p2 = self.nodes[a_triangle.node2 as usize - 1];
            let p3 = self.nodes[a_triangle.node3 as usize - 1];
            let (is_crossing, hide_before) = if d1 > the_triangle.tolerance {
                if d2 < -the_triangle.tolerance {
                    the_triangle.param = d1 / (d1 - d2);
                    (true, false)
                } else {
                    continue;
                }
            } else if d1 < -the_triangle.tolerance {
                if d2 > the_triangle.tolerance {
                    the_triangle.param = d1 / (d1 - d2);
                    (true, true)
                } else {
                    (false, false)
                }
            } else if d2 < -the_triangle.tolerance {
                (false, false)
            } else {
                continue;
            };
            // `cxx:159-167`, `:178-186`, `:188-199`, `:202-212`: the three
            // arms that reach `hideByOneTriangle` all refresh the projected
            // vertices first.
            the_triangle.v1 = GpXY::new(p1.x(), p1.y());
            the_triangle.v2 = GpXY::new(p2.x(), p2.y());
            the_triangle.v3 = GpXY::new(p3.x(), p3.y());
            hide_by_one_triangle(
                the_points,
                the_triangle,
                is_crossing,
                hide_before,
                a_triangle.flags,
                status,
            );
        }
    }
}

/// The mutable cursor state shared by the three edge blocks of
/// `hideByOneTriangle` (`cxx:229-238`).
struct HideState {
    cros_seg: bool,
    n1: i32,
    pd1: f64,
    pd2: f64,
    nn1: i32,
    nn2: i32,
    pend: f64,
    psta: f64,
    pp: f64,
    pdp: f64,
    npi: i32,
    o: [bool; 2],
    m: [bool; 2],
    p: [f64; 2],
    npi_rej: i32,
}

impl HideState {
    fn new() -> Self {
        Self {
            cros_seg: false,
            n1: 0,
            pd1: 0.0,
            pd2: 0.0,
            nn1: 0,
            nn2: 0,
            pend: 1.0,
            psta: 0.0,
            pp: 0.0,
            pdp: 0.0,
            npi: -1,
            o: [false, false],
            m: [false, false],
            p: [0.0, 0.0],
            npi_rej: 0,
        }
    }
}

/// `HLRAlgo_PolyData::hideByOneTriangle` (`cxx:222-1127`).
///
/// OCCT repeats the same 300-line block three times, once per triangle edge
/// (`V1->V2` at `:240-504`, `V2->V3` at `:506-771`, `V3->V1` at `:773-1037`).
/// The blocks differ only in the base vertex, the next vertex and the
/// `EMskGrALin{1,2,3}` flag tested, so the port factors them into
/// [`hide_by_one_edge`]; the order of the assignments and the two early
/// `return`s of each block are preserved.
fn hide_by_one_triangle(
    the_points: &BiPointPoints,
    the_triangle: &mut PolyTriangle,
    crossing: bool,
    hide_before: bool,
    tr_flags: i32,
    status: &mut HlrAlgoEdgeStatus,
) {
    let mut st = HideState::new();
    if !hide_by_one_edge(&mut st, the_triangle, the_points, the_triangle.v1, the_triangle.v2, tr_flags, poly_mask::EMSK_GR_A_LIN1) {
        return;
    }
    if !hide_by_one_edge(&mut st, the_triangle, the_points, the_triangle.v2, the_triangle.v3, tr_flags, poly_mask::EMSK_GR_A_LIN2) {
        return;
    }
    if !hide_by_one_edge(&mut st, the_triangle, the_points, the_triangle.v3, the_triangle.v1, tr_flags, poly_mask::EMSK_GR_A_LIN3) {
        return;
    }

    // `cxx:1039-1073`.
    if st.npi == -1 {
        if st.npi_rej >= 2 {
            return;
        }
    } else if st.npi == 0 {
        if st.o[0] {
            st.psta = st.p[0];
            st.pend = 1.0;
        } else {
            st.psta = 0.0;
            st.pend = st.p[0];
        }
    } else if st.npi == 1 {
        if st.p[0] > st.p[1] {
            st.psta = st.p[1];
            st.pend = st.p[0];
        } else {
            st.psta = st.p[0];
            st.pend = st.p[1];
        }
    }

    // `cxx:1075-1097`.
    if crossing {
        if hide_before {
            if the_triangle.param - st.psta < the_triangle.tol_param {
                return;
            } else if the_triangle.param < st.pend {
                st.pend = the_triangle.param;
            }
        } else if st.pend - the_triangle.param < the_triangle.tol_param {
            return;
        } else if st.psta < the_triangle.param {
            st.psta = the_triangle.param;
        }
    }

    // `cxx:1099-1121`.
    let mut total = if st.psta > 0.0 {
        st.psta < the_triangle.tol_param
    } else {
        st.psta > -the_triangle.tol_param
    };
    if total {
        let pfin = st.pend - 1.0;
        total = if pfin > 0.0 {
            pfin < the_triangle.tol_param
        } else {
            pfin > -the_triangle.tol_param
        };
    }
    if total {
        status.hide_all();
    } else {
        status.hide(
            st.psta,
            the_triangle.tol_param as f32,
            st.pend,
            the_triangle.tol_param as f32,
            false,
            false,
        );
    }
}

/// One edge block of `hideByOneTriangle`: clip the segment against the
/// infinite line `v_base -> v_next` and record the resulting parameter, exactly
/// as one of the three copies at `cxx:240-504` / `:506-771` / `:773-1037`.
/// Returns `false` where the OCCT block executes `return`.
fn hide_by_one_edge(
    st: &mut HideState,
    the_triangle: &PolyTriangle,
    the_points: &BiPointPoints,
    v_base: GpXY,
    v_next: GpXY,
    tr_flags: i32,
    gr_a_mask: i32,
) -> bool {
    let a_d = v_next.subtracted(&v_base);
    let a_a = GpXY::new(-a_d.y(), a_d.x()).multiplied_scalar(1.0 / a_d.modulus());
    let a_dot = a_a.dot(&v_base);
    let d1 = a_a.dot(&the_points.pntp12d()) - a_dot;
    let d2 = a_a.dot(&the_points.pntp22d()) - a_dot;
    if d1 > the_triangle.tolerance {
        if d2 < -the_triangle.tolerance {
            st.n1 = 2;
            st.cros_seg = true;
        } else {
            st.cros_seg = false;
        }
    } else if d1 < -the_triangle.tolerance {
        if d2 > the_triangle.tolerance {
            st.n1 = -1;
            st.cros_seg = true;
        } else {
            return false;
        }
    } else if d2 > the_triangle.tolerance {
        st.cros_seg = false;
    } else if d2 < -the_triangle.tolerance {
        return false;
    } else {
        st.cros_seg = false;
        if tr_flags & gr_a_mask != 0 {
            st.pd1 = (the_points.pntp1.x() - v_base.x()) / a_d.x();
            st.pd2 = (the_points.pntp2.x() - v_base.x()) / a_d.x();
        } else {
            st.pd1 = (the_points.pntp1.y() - v_base.y()) / a_d.y();
            st.pd2 = (the_points.pntp2.y() - v_base.y()) / a_d.y();
        }
        if st.pd1 < -the_triangle.tol_param {
            st.nn1 = 1;
        } else if st.pd1 < the_triangle.tol_param {
            st.nn1 = 2;
        } else if st.pd1 - 1.0 < -the_triangle.tol_param {
            st.nn1 = 3;
        } else if st.pd1 - 1.0 < the_triangle.tol_param {
            st.nn1 = 4;
        } else {
            st.nn1 = 5;
        }
        if st.pd2 < -the_triangle.tol_param {
            st.nn2 = 1;
        } else if st.pd2 < the_triangle.tol_param {
            st.nn2 = 2;
        } else if st.pd2 - 1.0 < -the_triangle.tol_param {
            st.nn2 = 3;
        } else if st.pd2 - 1.0 < the_triangle.tol_param {
            st.nn2 = 4;
        } else {
            st.nn2 = 5;
        }
        if st.nn1 == 3 {
            if st.nn2 == 1 {
                st.pend = st.pd1 / (st.pd1 - st.pd2);
            } else if st.nn2 == 5 {
                st.pend = (1.0 - st.pd1) / (st.pd2 - st.pd1);
            }
        } else if st.nn1 == 1 {
            if st.nn2 <= 2 {
                return false;
            } else {
                st.psta = -st.pd1 / (st.pd2 - st.pd1);
                if st.nn2 == 5 {
                    st.pend = (1.0 - st.pd1) / (st.pd2 - st.pd1);
                }
            }
        } else if st.nn1 == 5 {
            if st.nn2 >= 4 {
                return false;
            } else {
                st.psta = (st.pd1 - 1.0) / (st.pd1 - st.pd2);
                if st.nn2 == 1 {
                    st.pend = st.pd1 / (st.pd1 - st.pd2);
                }
            }
        } else if st.nn1 == 2 {
            if st.nn2 == 1 {
                return false;
            } else if st.nn2 == 5 {
                st.pend = (1.0 - st.pd1) / (st.pd2 - st.pd1);
            }
        } else if st.nn1 == 4 {
            if st.nn2 == 5 {
                return false;
            } else if st.nn2 == 1 {
                st.pend = st.pd1 / (st.pd1 - st.pd2);
            }
        }
    }
    if !st.cros_seg {
        return true;
    }
    let ad1 = if d1 < 0.0 { -d1 } else { d1 };
    let ad2 = if d2 < 0.0 { -d2 } else { d2 };
    st.pp = ad1 / (ad1 + ad2);
    if tr_flags & gr_a_mask != 0 {
        st.pdp = (the_points.pntp1.x()
            + (the_points.pntp2.x() - the_points.pntp1.x()) * st.pp
            - v_base.x())
            / a_d.x();
    } else {
        st.pdp = (the_points.pntp1.y()
            + (the_points.pntp2.y() - the_points.pntp1.y()) * st.pp
            - v_base.y())
            / a_d.y();
    }
    let mut out_side_p = false;
    let mut multiple = false;
    if st.pdp < -the_triangle.tol_param {
        out_side_p = true;
    } else if st.pdp < the_triangle.tol_param {
        multiple = true;
        out_side_p = st.multiple_cleanup();
    } else if st.pdp - 1.0 < -the_triangle.tol_param {
        // `cxx:452-454`: empty arm.
    } else if st.pdp - 1.0 < the_triangle.tol_param {
        multiple = true;
        out_side_p = st.multiple_cleanup();
    } else {
        out_side_p = true;
    }
    if out_side_p {
        st.npi_rej += 1;
    } else {
        st.npi += 1;
        if st.npi < 2 {
            st.p[st.npi as usize] = st.pp;
            st.o[st.npi as usize] = st.n1 == -1;
            st.m[st.npi as usize] = multiple;
        }
    }
    true
}

impl HideState {
    /// The `Multiple` cleanup loop duplicated at `cxx:443-462` and `:721-745`:
    /// drop every already registered hit that sits on the opposite side.
    /// Returns the `OutSideP` flag the OCCT loop computes.
    fn multiple_cleanup(&mut self) -> bool {
        let mut out_side_p = false;
        // `cxx:433` / `:460` / `:700` / `:727` / `:965` / `:992`: the loop
        // condition re-reads `npi`, which the body decrements, so the bound
        // shrinks while iterating.
        let mut l = 0;
        while l <= self.npi {
            if self.m[l as usize] {
                out_side_p = true;
                if self.o[l as usize] != (self.n1 == -1) {
                    if l == 0 && self.npi == 1 {
                        self.p[0] = self.p[1];
                        self.o[0] = self.o[1];
                        self.m[0] = self.m[1];
                    }
                    self.npi -= 1;
                    self.npi_rej += 1;
                }
            }
            l += 1;
        }
        out_side_p
    }
}
