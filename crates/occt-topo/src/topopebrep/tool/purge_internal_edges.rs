//! `TopOpeBRepTool_PurgeInternalEdges` — drop edges that no face carries.
//!
//! Source: `TKBool/TopOpeBRepTool/TopOpeBRepTool_PurgeInternalEdges.cxx`
//! (OCCT 8.0.0, 264 lines) and `.hxx:43-88`.
//!
//! | OCCT | lines |
//! |---|---|
//! | constructor | `cxx:36-48` |
//! | `Faces` | `cxx:52-63` |
//! | `NbEdges` | `cxx:67-91` |
//! | `Shape` | `cxx:95-100` |
//! | `BuildList` | `cxx:108-227` |
//! | `Perform` | `cxx:231-266` |
//!
//! `BuildList` fills `myMapFacLstEdg` (face -> internal edges) from
//! `myMapEdgLstFac` (edge -> ancestor faces, built by
//! `TopExp::MapShapesAndAncestors`). An edge qualifies when **every** face it
//! touches carries it as `TopAbs_INTERNAL` or `TopAbs_EXTERNAL` — a single
//! `FORWARD` / `REVERSED` occurrence keeps it (`cxx:156-177`).
//!
//! Four OCCT behaviours carry the semantics and are reproduced literally:
//!
//! 1. `cxx:146`: the `mapUniqEdg` guard is declared *inside* the per-edge
//!    block, so after the first face of `LstFacToTreat` takes the edge the
//!    guard suppresses every later face. Only `LstFacToTreat.first()` is ever
//!    written. This does not change `Perform`: the face key there is used only
//!    for `Find` (`cxx:255-256`), and `BRepTools_Substitution::Build` rebuilds
//!    every ancestor of the removed edge anyway.
//! 2. `cxx:189`: on an edge with **no** face ancestor the exhausted list
//!    iterator hands back a null `TopoDS_Shape`, which is then used as a
//!    `DataMap` key. A port has no null `TopoDS_Shape`, so those edges go to
//!    [`PurgeInternalEdges::my_orphan_edg`] and are kept out of [`Self::faces`]
//!    (a null-face entry would mean nothing to a caller) while still being
//!    substituted away by [`Self::perform`] — the outcome of `cxx:258`.
//! 3. `cxx:182-197`: the `else` arm covers `Extent() == 0` as well as `== 1`,
//!    so the null-face key above is reachable from two directions.
//! 4. `NbEdges` (`cxx:67-91`) is `const` and does **not** trigger
//!    `BuildList`; a `PerformNow == false` instance reports 0 until `Faces` or
//!    `Perform` runs.
//!
//! UNPORTED: none. The `Standard_NullObject_Raise_if` guards (`cxx:41`,
//! `cxx:71`, `cxx:98`) are unreachable here — a `TopoShape` is never null.

use std::collections::HashMap;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tools_substitution::BRepToolsSubstitution;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key, ShapeMapKey};

use super::fuse_edges::IndexedAncestors;

/// `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
/// TopTools_ShapeMapHasher>` (`hxx:86-87`) in insertion order.
type FaceEdges = (TopoShape, Vec<TopoShape>);

/// `class TopOpeBRepTool_PurgeInternalEdges` (`hxx:43-88`).
pub struct PurgeInternalEdges {
    /// `myShape` (`hxx:84`).
    my_shape: TopoShape,
    /// `myIsDone` (`hxx:85`).
    my_is_done: bool,
    /// `myMapEdgLstFac` (`hxx:74-76`): `EDGE -> FACE`, non-unique ancestors.
    my_map_edg_lst_fac: IndexedAncestors,
    /// `myMapFacLstEdg` (`hxx:86-87`): `FACE -> [EDGE]`.
    my_map_fac_lst_edg: Vec<FaceEdges>,
    /// `myMapFacLstEdg.IsBound` / `FindIndex` for [`Self::my_map_fac_lst_edg`].
    fac_index: HashMap<ShapeMapKey, usize>,
    /// `myMapFacLstEdg` entries keyed by the null `TopoDS_Shape` (`cxx:189`).
    my_orphan_edg: Vec<TopoShape>,
}

impl PurgeInternalEdges {
    /// Constructor (`cxx:36-48`). `PerformNow` selects the `Perform()` call of
    /// `cxx:45-47`.
    pub fn new(shape: &TopoShape, perform_now: bool) -> Self {
        let mut this = Self {
            my_shape: shape.clone(),
            my_is_done: false,
            my_map_edg_lst_fac: IndexedAncestors::new(),
            my_map_fac_lst_edg: Vec::new(),
            fac_index: HashMap::new(),
            my_orphan_edg: Vec::new(),
        };
        if perform_now {
            this.perform();
        }
        this
    }

    /// `Faces(theMapFacLstEdg)` (`cxx:52-63`): the internal edges per face.
    ///
    /// Runs `BuildList` when it has not run yet (`cxx:57-60`). The null-face
    /// entries OCCT would hand back are omitted (see module notes).
    pub fn faces(&mut self) -> Vec<FaceEdges> {
        if !self.my_is_done {
            self.build_list();
        }
        self.my_map_fac_lst_edg.clone()
    }

    /// `NbEdges()` (`cxx:67-91`): candidates to remove, empty list included.
    ///
    /// `const` in OCCT, so it does not build the list.
    pub fn nb_edges(&self) -> usize {
        self.my_map_fac_lst_edg.iter().map(|(_, l)| l.len()).sum::<usize>()
            + self.my_orphan_edg.len()
    }

    /// `Shape()` (`cxx:95-100`).
    pub fn shape(&self) -> &TopoShape {
        &self.my_shape
    }

    /// `IsDone()` (`hxx:68`).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// `BuildList()` (`cxx:108-227`).
    fn build_list(&mut self) {
        // `cxx:119-120`
        self.my_map_edg_lst_fac = IndexedAncestors::map_shapes_and_ancestors(
            &self.my_shape,
            ShapeType::Edge,
            ShapeType::Face,
        );
        self.my_map_fac_lst_edg.clear();
        self.fac_index.clear();
        self.my_orphan_edg.clear();

        // `cxx:132`: `LstFacToTreat`. `None` is the null `TopoDS_Shape` that
        // `cxx:189` appends when the ancestor list is empty.
        let mut lst_fac_to_treat: Vec<Option<TopoShape>> = Vec::new();

        // `cxx:134-226`
        for i_edg in 1..=self.my_map_edg_lst_fac.extent() {
            let edgecur = match self.my_map_edg_lst_fac.find_key(i_edg) {
                Some(e) => e.clone(),
                None => continue,
            };
            let lmap_fac: Vec<TopoShape> = self
                .my_map_edg_lst_fac
                .find_from_key(&edgecur)
                .map(|l| l.to_vec())
                .unwrap_or_default();

            // `cxx:143-144`
            lst_fac_to_treat.clear();

            if lmap_fac.len() > 1 {
                // `cxx:148-183`
                let mut to_keep = false;
                for facecur in &lmap_fac {
                    if to_keep {
                        break;
                    }
                    // `cxx:157`: `TopExp_Explorer(facecur, TopAbs_EDGE)` — the
                    // edge carries its orientation in this face.
                    for edge in explorer_shapes(facecur, ShapeType::Edge) {
                        let orien = edge.orientation();
                        if is_same(&edge, &edgecur) {
                            // `cxx:166-178`
                            if orien == Orientation::Internal || orien == Orientation::External {
                                lst_fac_to_treat.push(Some(facecur.clone()));
                            } else {
                                lst_fac_to_treat.clear();
                                to_keep = true;
                            }
                            break;
                        }
                    }
                }
            } else {
                // `cxx:185-192`
                let orien = edgecur.orientation();
                if orien == Orientation::Internal || orien == Orientation::External {
                    match lmap_fac.first() {
                        Some(f) => lst_fac_to_treat.push(Some(f.clone())),
                        None => lst_fac_to_treat.push(None),
                    }
                }
            }

            // `cxx:194-225`: OCCT writes only the first face (module note 1).
            match lst_fac_to_treat.first() {
                Some(Some(face)) => {
                    let face = face.clone();
                    self.bind_or_append(&face, &edgecur);
                }
                Some(None) => self.my_orphan_edg.push(edgecur.clone()),
                None => {}
            }
        }

        // `cxx:226`
        self.my_is_done = true;
    }

    /// `myMapFacLstEdg.Bind(face, [edge])` / `ChangeFind(face).Append(edge)`
    /// (`cxx:204-223`), collapsed to one call as derived in module note 1.
    fn bind_or_append(&mut self, face: &TopoShape, edge: &TopoShape) {
        let key = shape_map_key(face);
        if let Some(&i) = self.fac_index.get(&key) {
            // `cxx:216-221`
            self.my_map_fac_lst_edg[i].1.push(edge.clone());
        } else {
            // `cxx:206-213`
            self.fac_index.insert(key, self.my_map_fac_lst_edg.len());
            self.my_map_fac_lst_edg.push((face.clone(), vec![edge.clone()]));
        }
    }

    /// `Perform()` (`cxx:231-266`): substitute every internal edge by an empty
    /// list, then rebuild `myShape` from `Bsub.Copy(myShape).First()`.
    pub fn perform(&mut self) {
        // `cxx:235-238`
        if !self.my_is_done {
            self.build_list();
        }

        // `cxx:241`
        if self.my_map_fac_lst_edg.is_empty() && self.my_orphan_edg.is_empty() {
            return;
        }

        // `cxx:247-261`
        let mut bsub = BRepToolsSubstitution::new();
        let empty: [TopoShape; 0] = [];
        for (_, lst) in &self.my_map_fac_lst_edg {
            for e in lst {
                bsub.substitute(e, &empty);
            }
        }
        for e in &self.my_orphan_edg {
            bsub.substitute(e, &empty);
        }

        // `cxx:263-267`
        bsub.build(&self.my_shape);
        if bsub.is_copied(&self.my_shape) {
            if let Some(first) = bsub.copy(&self.my_shape).and_then(|l| l.first()) {
                self.my_shape = first.clone();
            }
        }
    }
}
