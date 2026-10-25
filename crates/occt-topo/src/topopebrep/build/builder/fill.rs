//! `TopOpeBRepBuild_Builder` — the fill family that loads a shape set.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Builder.cxx:1876-1949`
//! (`FillShape`), `:1951-1977` (`FillFace`), `:1979-1989` (`FillSolid`),
//! `:1991-2020` (`FillVertexSet`), `:2023-2057` (`FillVertexSetOnValue`)
//! (OCCT 8.0.0).
//!
//! `FillShape` walks the sub-shapes of `S1` one level above the ones that end
//! up in the set: the wires of a face, the shells of a solid. A sub-shape the
//! DS does not know is classified against the same-domain partners of the
//! other operand and, when it survives, enters the set whole; a sub-shape with
//! geometry is split into its own sub-shapes by
//! [`Builder::split_shapes`](super::splits).
//!
//! A same-domain shape whose reference orientation differs
//! (`TopOpeBRepDS_DIFFORIENTED`) flips `RevOri`, so every shape taken from it
//! is reversed once more.

use std::sync::{Arc, RwLock};

use crate::abs::{Orientation, ShapeType, State};
use crate::shape::TopoShape;
use crate::topopebrep::ds::kind::Config;
use crate::topopebrep::ds::point_iterator::PointIterator;
use crate::topopebrep::tool::shape_explorer::ShapeExplorer;

use super::super::loop_shape::Loop;
use super::super::pave_set::PaveSet;
use super::super::shape_set::ShapeSetSource;
use super::super::wire_edge_set::WireEdgeSet;
use super::Builder;

impl Builder {
    /// `FillShape(S1, ToBuild1, LS2, ToBuild2, aSet, In_RevOri)`
    /// (`cxx:1876-1949`).
    pub fn fill_shape(
        &mut self,
        s1: &TopoShape,
        to_build1: State,
        ls2: &[TopoShape],
        to_build2: State,
        a_set: &mut dyn ShapeSetSource,
        in_rev_ori: bool,
    ) {
        let mut rev_ori = in_rev_ori;
        let t = s1.shape_type();

        // cxx:1886-1895: the explored type and the type of its sub-shapes.
        let (t1, t11) = match t {
            ShapeType::Face => (ShapeType::Wire, ShapeType::Edge),
            ShapeType::Solid | ShapeType::Shell => (ShapeType::Shell, ShapeType::Face),
            _ => (ShapeType::Compound, ShapeType::Compound),
        };

        // cxx:1900-1916: the orientation of a same-domain shape relative to
        // the reference one of its group.
        let hsd = self
            .data_structure()
            .map(|h| h.has_same_domain(s1, true))
            .unwrap_or(false);
        if hsd {
            if let Some(hds) = self.data_structure() {
                if hds.same_domain_orientation(s1) == Config::DiffOriented {
                    rev_ori = !rev_ori;
                }
            }
        }

        // cxx:1918-1919: work on a FORWARD shape <aShape>.
        let mut a_shape = s1.clone();
        self.my_build_tool
            .set_orientation(&mut a_shape, Orientation::Forward);

        // cxx:1924-1947: explore the sub-shapes of type <t1>.
        let mut ex1 = ShapeExplorer::with_shape(&a_shape, t1, ShapeType::Shape);
        while ex1.more() {
            let Some(a_sub_shape) = ex1.current().cloned() else {
                break;
            };

            let has = self
                .data_structure()
                .map(|h| h.has_shape(&a_sub_shape, true))
                .unwrap_or(false);

            if !has {
                // cxx:1929-1938: not in the DS, classify it against LS2.
                let keep = self.keep_shape(&a_sub_shape, ls2, to_build1);
                if keep {
                    let newori =
                        Self::orient(self.my_build_tool.orientation(&a_sub_shape), rev_ori);
                    let mut a_sub_shape = a_sub_shape;
                    self.my_build_tool.set_orientation(&mut a_sub_shape, newori);
                    a_set.add_shape(&a_sub_shape);
                }
            } else {
                // cxx:1939-1945: split the sub-shapes of type <t11>.
                let mut ex11 = ShapeExplorer::with_shape(&a_sub_shape, t11, ShapeType::Shape);
                self.split_shapes(&mut ex11, to_build1, to_build2, a_set, rev_ori);
            }

            ex1.next();
        }
    }

    /// `FillFace(F1, ToBuild1, LF2, ToBuild2, WES, RevOri)` (`cxx:1951-1977`).
    ///
    /// `myListOfFace` is the 2d-pure context `KeepShape` reads, so it lives
    /// only for the duration of the `FillShape` call.
    pub fn fill_face(
        &mut self,
        f1: &TopoShape,
        to_build1: State,
        lf2: &[TopoShape],
        to_build2: State,
        wes: &mut WireEdgeSet,
        rev_ori: bool,
    ) {
        self.my_list_of_face = lf2.to_vec();
        self.fill_shape(f1, to_build1, lf2, to_build2, wes, rev_ori);
        self.my_list_of_face.clear();
    }

    /// `FillSolid(S1, ToBuild1, LS2, ToBuild2, aSet, RevOri)` (`cxx:1979-1989`).
    pub fn fill_solid(
        &mut self,
        s1: &TopoShape,
        to_build1: State,
        ls2: &[TopoShape],
        to_build2: State,
        a_set: &mut dyn ShapeSetSource,
        rev_ori: bool,
    ) {
        self.fill_shape(s1, to_build1, ls2, to_build2, a_set, rev_ori);
    }

    /// `FillVertexSet(IT, ToBuild, PVS)` (`cxx:1991-2000`).
    pub fn fill_vertex_set(&self, it: &mut PointIterator, to_build: State, pvs: &mut PaveSet) {
        while it.more() {
            self.fill_vertex_set_on_value(it, to_build, pvs);
            it.next();
        }
    }

    /// `FillVertexSetOnValue(IT, ToBuild, PVS)` (`cxx:2023-2057`).
    ///
    /// `IT.Current()` is either an index in the new-point array of the builder
    /// or a shape index in the DS; `NewVertex` answers the first case and
    /// `DataStructure::Shape` the second (`cxx:2010-2021`).
    pub fn fill_vertex_set_on_value(
        &self,
        it: &PointIterator,
        to_build: State,
        pvs: &mut PaveSet,
    ) {
        let ind = it.geometry();
        let ispoint = it.is_point();

        let Some(hds) = self.data_structure() else {
            return;
        };

        let v = if ispoint && ind <= hds.nb_points() {
            self.new_vertex(ind).cloned()
        } else {
            hds.shape(ind, false)
        };
        // `cxx:2029`: `keep` is left permanently true in OCCT (the
        // `EXTERNAL` / `INTERNAL` test above it is commented out).
        let keep = true;

        // `cxx:2018`: `V` is a `TopoDS_Shape`; a null `NewVertex` answer stays
        // a null shape in OCCT, which the port spells as an empty vertex.
        let mut v = v.unwrap_or_else(|| TopoShape::new(ShapeType::Vertex));
        let par = it.parameter().unwrap_or(0.0);
        let ori = it.orientation(to_build);

        if keep {
            self.my_build_tool.set_orientation(&mut v, ori);
            let pv = Arc::new(RwLock::new(Loop::new_pave(&v, par, false)));
            pvs.append(pv);
        }
    }
}