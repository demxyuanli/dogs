//! `BRepMAT2d_BisectingLocus` (`BRepMAT2d_BisectingLocus.cxx:35-359`,
//! `.hxx:35-115`): runs `MAT2d_Mat2d` over the contours extracted by
//! [`crate::brep_mat2d_explorer::BRepMat2dExplorer`] and feeds the resulting
//! roots into `MAT_Graph::Perform`.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `Compute` | `cxx:57-166` | [`BRepMat2dBisectingLocus::compute`] |
//! | `RenumerationAndFusion` | `cxx:170-208` | `renumeration_and_fusion` |
//! | `IsDone` / `Graph` / `NumberOfContours` / `NumberOfElts` | `cxx:214-243` | same names |
//! | `NumberOfSections` / `BasicElt` / `GeomBis` / `GeomElt` | `cxx:248-316` | same names |
//! | `CutSketch` | `cxx:325-357` | [`cut_sketch`] |
//!
//! Deviations:
//!
//! * OCCT keeps `theGraph` (`MAT_Graph`) separate from the `MAT2d_Mat2d` that
//!   produced the bisectors, and `MAT_Graph::Perform` reads them through shared
//!   handles. The port's `MatGraph` is also the `MAT_Bisector` / `MAT_Edge`
//!   arena, so [`BRepMat2dBisectingLocus`] owns the `Mat2dMat2d` and exposes its
//!   graph through [`BRepMat2dBisectingLocus::graph`]; `theGraph = new MAT_Graph()`
//!   in `cxx:68` becomes a fresh `Mat2dMat2d`.
//! * `MAT_Side` / `GeomAbs_JoinType` are the `occt_topo::mat::MatSide` /
//!   `occt_core::kernel::geomabs::JoinType` enums.
//! * `CutSketch` replaces each cut curve in place with its pieces; the
//!   `InsertAfter` + `Remove` dance of `cxx:347-352` is one `splice`.

use std::collections::HashMap;
use std::rc::Rc;

use occt_core::kernel::geomabs::JoinType;
use occt_core::precision::INFINITE;
use occt_geom2d::bisector::BisectorBisec;
use occt_geom2d::curve::Curve2d;

use crate::brep_mat2d_explorer::BRepMat2dExplorer;
use crate::mat::{BasicEltId, BisectorId, MatArc, MatBasicElt, MatGraph, MatSide};
use crate::mat2d::{
    Mat2dBiInt, Mat2dCircuit, Mat2dCutCurve, Mat2dFigure, Mat2dItem, Mat2dMat2d, Mat2dTool2d,
};

/// `class BRepMAT2d_BisectingLocus` (`BRepMAT2d_BisectingLocus.hxx:35-115`).
#[derive(Default)]
pub struct BRepMat2dBisectingLocus {
    /// `isDone` (`.hxx`).
    is_done: bool,
    /// `nbContours` (`.hxx`): the number of contours of the last `Compute`.
    nb_contours: i32,
    /// `nbSect` (`.hxx`), an `NCollection_DataMap<MAT2d_BiInt, int>`.
    nb_sect: HashMap<Mat2dBiInt, i32>,
    /// `theTool` (`.hxx`).
    the_tool: Mat2dTool2d,
    /// `TheMAT` of the last `Compute`, owning the `theGraph` arena.
    the_mat: Mat2dMat2d,
}

impl BRepMat2dBisectingLocus {
    /// `BRepMAT2d_BisectingLocus()` (`cxx:44-47`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Compute(anExplo, IndexLine, aSide, aJoinType, IsOpenResult)`
    /// (`cxx:57-166`): cuts the sketch at its curvature extrema, builds the
    /// `MAT2d_Circuit`, runs `MAT2d_Mat2d`, and hands the roots to
    /// `MAT_Graph::Perform`.
    pub fn compute(
        &mut self,
        an_explo: &mut BRepMat2dExplorer,
        index_line: i32,
        a_side: MatSide,
        a_join_type: JoinType,
        is_open_result: bool,
    ) {
        // `cxx:67-68`: `nbSect.Clear(); theGraph = new MAT_Graph()`.
        self.nb_sect.clear();
        self.the_mat = Mat2dMat2d::new(is_open_result);

        // `cxx:69-72`.
        self.nb_contours = an_explo.number_of_contours() as i32;
        if self.nb_contours == 0 {
            return;
        }

        // `cxx:77-86`: the figure, one line per contour.
        let mut figure: Mat2dFigure = Vec::new();
        for i in 1..=an_explo.number_of_contours() {
            let mut line: Vec<Mat2dItem> = Vec::new();
            an_explo.init(i);
            while an_explo.more() {
                line.push(Mat2dItem::Curve(an_explo.value().clone()));
                an_explo.next();
            }
            figure.push(line);
        }

        // `cxx:90`: curve splitting.
        cut_sketch(&mut figure, &mut self.nb_sect);

        // `cxx:97-99`: the circuit the map is computed on.
        let mut a_circuit = Mat2dCircuit::new(a_join_type, is_open_result);
        a_circuit.perform(
            &mut figure,
            an_explo.get_is_closed(),
            index_line,
            a_side == MatSide::Left,
        );
        let a_circuit = Rc::new(a_circuit);

        // `cxx:104-106`.
        self.the_tool.sense(a_side);
        self.the_tool.set_join_type(a_join_type);
        self.the_tool.init_items(a_circuit);

        // `cxx:110-117`.
        if is_open_result {
            self.the_mat.create_mat_open(&mut self.the_tool);
        } else {
            self.the_mat.create_mat(&mut self.the_tool);
        }

        // `cxx:120-123`.
        self.is_done = self.the_mat.is_done();
        if !self.is_done {
            return;
        }

        // `cxx:129-133`: `TheRoots->BackAdd(TheMAT.Bisector())` in `Init` /
        // `More` / `Next` order.
        let mut the_roots: Vec<BisectorId> = Vec::new();
        self.the_mat.init();
        while self.the_mat.more() {
            the_roots.push(self.the_mat.bisector());
            self.the_mat.next();
        }

        // `cxx:135-138`.
        let semi_infinite = self.the_mat.semi_infinite();
        let nb_items = self.the_tool.number_of_items();
        let nb_bisectors = self.the_mat.number_of_bisectors();
        self.the_mat
            .graph_mut()
            .perform(semi_infinite, &the_roots, nb_items, nb_bisectors);

        // `cxx:143-165`: fusion of the basic elements duplicated when the
        // explorer holds several lines.
        if an_explo.number_of_contours() > 1 {
            let mut new_map: std::collections::BTreeMap<i32, BasicEltId> =
                std::collections::BTreeMap::new();
            let mut index_last = 1;
            for i in 1..=an_explo.number_of_contours() as i32 {
                let length_line = self.the_tool.circuit().line_length(i);
                self.renumeration_and_fusion(i, length_line, &mut index_last, &mut new_map);
            }
            self.the_mat.graph_mut().change_basic_elts(&new_map);
            self.the_mat.graph_mut().compact_arcs();
            self.the_mat.graph_mut().compact_nodes();
        }
    }

    /// `RenumerationAndFusion(ILine, LengthLine, IndexLast, NewMap)`
    /// (`cxx:170-208`): renumbers the basic elements in contour order (dropping
    /// duplicates) and fuses the arcs that became shared.
    fn renumeration_and_fusion(
        &mut self,
        i_line: i32,
        length_line: i32,
        index_last: &mut i32,
        new_map: &mut std::collections::BTreeMap<i32, BasicEltId>,
    ) {
        for i in 1..=length_line {
            // `cxx:186`: `theTool.Circuit()->RefToEqui(ILine, i)`.
            let s: Vec<i32> = self.the_tool.circuit().ref_to_equi(i_line, i).to_vec();

            // `cxx:188-190`.
            let ind_first = s[0];
            new_map.insert(*index_last, self.the_mat.graph().change_basic_elt(ind_first));
            *index_last += 1;

            // `cxx:192-208`.
            for &elt in &s[1..] {
                let mut merge_arc1 = false;
                let mut merge_arc2 = false;
                let mut geom_index_arc1 = 0;
                let mut geom_index_arc2 = 0;
                let mut geom_index_arc3 = 0;
                let mut geom_index_arc4 = 0;
                self.the_mat.graph_mut().fusion_of_basic_elts(
                    ind_first,
                    elt,
                    &mut merge_arc1,
                    &mut geom_index_arc1,
                    &mut geom_index_arc2,
                    &mut merge_arc2,
                    &mut geom_index_arc3,
                    &mut geom_index_arc4,
                );
                if merge_arc1 {
                    self.the_tool.bisec_fusion(geom_index_arc1, geom_index_arc2);
                }
                if merge_arc2 {
                    self.the_tool.bisec_fusion(geom_index_arc3, geom_index_arc4);
                }
            }
        }
    }

    /// `IsDone()` (`cxx:214-217`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `Graph()` (`cxx:222-225`).
    pub fn graph(&self) -> &MatGraph {
        self.the_mat.graph()
    }

    /// `NumberOfContours()` (`cxx:232-235`).
    pub fn number_of_contours(&self) -> i32 {
        self.nb_contours
    }

    /// `NumberOfElts(IndLine)` (`cxx:240-243`).
    pub fn number_of_elts(&self, ind_line: i32) -> i32 {
        self.the_tool.circuit().line_length(ind_line)
    }

    /// `NumberOfSections(IndLine, Index)` (`cxx:248-252`).
    pub fn number_of_sections(&self, ind_line: i32, index: i32) -> i32 {
        let b = Mat2dBiInt::new(ind_line, index);
        *self.nb_sect.get(&b).unwrap_or(&0)
    }

    /// `BasicElt(IndLine, Index)` (`cxx:257-270`): `Index` is the 1-based index
    /// inside `IndLine`, the answer the 1-based graph index.
    pub fn basic_elt(&self, ind_line: i32, index: i32) -> BasicEltId {
        let mut ind = index;
        for i in 1..ind_line {
            ind += self.the_tool.circuit().line_length(i);
        }
        self.the_mat.graph().basic_elt(ind)
    }

    /// `GeomBis(anArc, Reverse)` (`cxx:274-298`): the arc's bisector, with
    /// `reverse` telling whether its start point is the far end of the arc.
    ///
    /// `Reverse` is `true` when the curve starts at `-Precision::Infinite()`,
    /// and otherwise when the first node is closer to the last point of a
    /// bounded bisector than to its first point.
    pub fn geom_bis(&self, an_arc: &MatArc, reverse: &mut bool) -> BisectorBisec {
        *reverse = false;

        let bis = self.the_tool.geom_bis(an_arc.geom_index());
        let bis_curve = bis
            .value()
            .unwrap_or_else(|| panic!("BRepMAT2d_BisectingLocus::GeomBis: null bisector"))
            .clone();

        if bis_curve.first_parameter() <= -INFINITE {
            // `cxx:284-287`.
            *reverse = true;
        } else if bis_curve.last_parameter() < INFINITE {
            // `cxx:288-296`.
            let p_f = bis_curve.d0(bis_curve.first_parameter());
            let p_l = bis_curve.d0(bis_curve.last_parameter());
            let node = an_arc
                .first_node()
                .unwrap_or_else(|| panic!("BRepMAT2d_BisectingLocus::GeomBis: null node"));
            let p_node = self.geom_elt_of_node(node);
            if p_node.square_distance(&p_f) > p_node.square_distance(&p_l) {
                *reverse = true;
            }
        }

        // `cxx:298`: the bisector itself, not the temporary curve above.
        self.the_tool.geom_bis(an_arc.geom_index()).clone()
    }

    /// `GeomElt(aBasicElt)` (`cxx:305-308`).
    pub fn geom_elt(&self, a_basic_elt: &MatBasicElt) -> &Mat2dItem {
        self.the_tool.geom_elt(a_basic_elt.geom_index())
    }

    /// `GeomElt(aNode)` (`cxx:314-317`): `theTool.GeomPnt(aNode->GeomIndex())`.
    pub fn geom_elt_of_node(&self, a_node: crate::mat::NodeId) -> occt_core::gp::GpPnt2d {
        let geom_index = self.the_mat.graph().node_by_id(a_node).geom_index();
        self.the_tool.geom_pnt(geom_index)
    }

    /// The bisector geometry of `arc` (`theTool.GeomBis(anArc->GeomIndex())`).
    pub fn the_tool(&self) -> &Mat2dTool2d {
        &self.the_tool
    }

    /// The arena of the last `Compute`.
    pub fn the_mat(&self) -> &Mat2dMat2d {
        &self.the_mat
    }
}

/// `CutSketch(Figure, NbSect)` (`cxx:325-357`): splits every curve of the figure
/// at its curvature extrema and inflection points, recording how many sections
/// each source curve became.
fn cut_sketch(figure: &mut Mat2dFigure, nb_sect: &mut HashMap<Mat2dBiInt, i32>) {
    let mut cuter = Mat2dCutCurve::new();

    for i in 1..=figure.len() {
        let contour = &mut figure[i - 1];
        let mut i_curve_init = 0i32;
        let mut j = 1usize;

        while j <= contour.len() {
            i_curve_init += 1;

            // `cxx:340`: `Cuter.Perform(down_cast<Geom2d_Curve>(...))`. Every
            // item of a figure built by `BRepMAT2d_Explorer` is a curve.
            cuter.perform(contour[j - 1].expect_curve());

            let mut nb_section = 1;
            if !cuter.un_modified() {
                // `cxx:343-353`: the pieces replace the source curve in place.
                nb_section = cuter.nb_curves() as i32;
                let pieces: Vec<Mat2dItem> = (1..=nb_section)
                    .map(|k| Mat2dItem::Curve(cuter.value(k).clone()))
                    .collect();
                contour.splice(j - 1..j, pieces);
                // The C++ `j--` after `Remove(ico)` then the loop's `j++`: the
                // cursor lands on the last inserted piece, so `j` is
                // `j - 1 + nb_section` here.
                j += nb_section as usize - 1;
            }

            // `cxx:354-355`.
            let b = Mat2dBiInt::new(i as i32, i_curve_init);
            nb_sect.insert(b, nb_section);
            j += 1;
        }
    }
}
