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

/// Crate imports shared by `mod.rs` and every submodule.
pub(super) use std::collections::HashMap;
pub(super) use std::sync::Arc;
pub(super) use occt_core::gp::{GpAx1, GpAx2, GpDir, GpLin, GpPln, GpPnt, GpPnt2d, GpVec};
pub(super) use occt_core::precision::{CONFUSION, REAL_SMALL, SQUARE_CONFUSION};
pub(super) use occt_geom::approx_same_parameter::u_resolution;
pub(super) use occt_geom::curve::Curve;
pub(super) use occt_geom::gcpnts::uniform_abscissa;
pub(super) use occt_geom::geom_lib::geom_lib_axe_of_inertia;
pub(super) use occt_geom::line::GeomLine;
pub(super) use crate::abs::{Orientation, ShapeType};
pub(super) use crate::bopds::BopdsDS;
pub(super) use crate::bopds_pave::BopdsPave;
pub(super) use crate::bopds_ff::BopdsCurve;
pub(super) use crate::brep_lib_make_edge::{add_endpoint_vertices, init_3d_vertices};
pub(super) use crate::brep_tool::BRepTool;
pub(super) use crate::builder::TopoBuilder;
pub(super) use crate::inttools_data::CommonPartType;
pub(super) use crate::pave_filler::PaveFiller;
pub(super) use crate::shape::{Edge, TopoShape, Vertex, Wire};
pub(super) use crate::shape_fix_compose_shell::{empty_copied_edge, MapReShape};
pub(super) use crate::shhealing::{first_vertex, last_vertex, ShapeFixEdge};
pub(super) use crate::tgeometry::GeometryRegistry;
pub(super) use crate::topo_tools::{
    common_vertex, edge_vertices, edges_of_wire, explorer_shapes, is_same,
    map_shapes_and_ancestors, shape_map_key, vertices_of, wire_vertices, ShapeMapKey,
};


pub(super) use super::helpers::ShapeTable;
pub(super) use super::TransitionStyle;

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
}

// ---------------------------------------------------------------------------
// Submodules
// ---------------------------------------------------------------------------

mod choose_section;
mod from_edges;
mod helpers;
mod make_faces;
mod section_walk;
mod split_u_edges;

use self::split_u_edges::split_u_edges;
