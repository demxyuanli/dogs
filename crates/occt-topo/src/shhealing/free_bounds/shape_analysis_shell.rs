//! `ShapeAnalysis_Shell` (`ShapeAnalysis_Shell.cxx`).

use super::*;

/// `ShapeAnalysis_Shell` (`ShapeAnalysis_Shell.hxx:26-131`).
pub struct ShapeAnalysisShell {
    /// `myShells`.
    shells: IndexedShapeMap,
    /// `myBad`.
    bad: IndexedShapeMap,
    /// `myFree`.
    free: IndexedShapeMap,
    /// `myConex`.
    conex: bool,
}

impl Default for ShapeAnalysisShell {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeAnalysisShell {
    /// `ShapeAnalysis_Shell()` (`cxx:28-31`).
    pub fn new() -> Self {
        Self {
            shells: IndexedShapeMap::new(),
            bad: IndexedShapeMap::new(),
            free: IndexedShapeMap::new(),
            conex: false,
        }
    }

    /// `Clear()` (`cxx:36-42`).
    pub fn clear(&mut self) {
        self.shells = IndexedShapeMap::new();
        self.bad = IndexedShapeMap::new();
        self.free = IndexedShapeMap::new();
        self.conex = false;
    }

    /// `LoadShells(shape)` (`cxx:46-65`).
    pub fn load_shells(&mut self, shape: &TopoShape) {
        if shape.shape_type() == ShapeType::Shell {
            self.shells.add(shape);
            return;
        }
        for sh in crate::topo_tools::explorer_shapes(shape, ShapeType::Shell) {
            self.shells.add(&sh);
        }
    }

    /// `CheckOrientedShells(shape, alsofree, checkinternaledges)`
    /// (`cxx:139-245`).
    pub fn check_oriented_shells(
        &mut self,
        shape: &TopoShape,
        alsofree: bool,
        checkinternaledges: bool,
    ) -> bool {
        // `cxx:143`.
        self.conex = false;

        let mut dirs = IndexedShapeMap::new();
        let mut revs = IndexedShapeMap::new();
        let mut ints = IndexedShapeMap::new();

        // `cxx:151-163`.
        let mut res = false;
        for sh in crate::topo_tools::explorer_shapes(shape, ShapeType::Shell) {
            if check_edges(&sh, &mut self.bad, &mut dirs, &mut revs, &mut ints) {
                if self.shells.add(&sh) {
                    res = true;
                }
            }
        }

        // `cxx:165-168`.
        if !alsofree {
            return res;
        }

        // `cxx:170-205`: every FORWARD-only edge is free.
        let nb = dirs.extent();
        for i in 1..=nb {
            let sh = dirs.find_key(i).cloned();
            let Some(sh) = sh else { continue };
            if !self.bad.contains(&sh) {
                if !revs.contains(&sh) {
                    if checkinternaledges {
                        if !ints.contains(&sh) {
                            self.free.add(&sh);
                        } else {
                            self.conex = true; // `cxx:228`.
                        }
                    } else {
                        self.free.add(&sh);
                    }
                } else {
                    self.conex = true; // `cxx:198`.
                }
            } else {
                self.conex = true; // `cxx:202`.
            }
        }

        // `cxx:207-241`: every REVERSED-only edge is free.
        let nb = revs.extent();
        for i in 1..=nb {
            let sh = revs.find_key(i).cloned();
            let Some(sh) = sh else { continue };
            if !self.bad.contains(&sh) {
                if !dirs.contains(&sh) {
                    if checkinternaledges {
                        if !ints.contains(&sh) {
                            self.free.add(&sh);
                        } else {
                            self.conex = true; // `cxx:222`.
                        }
                    } else {
                        self.free.add(&sh);
                    }
                } else {
                    self.conex = true; // `cxx:233`.
                }
            } else {
                self.conex = true; // `cxx:237`.
            }
        }

        res
    }

    /// `IsLoaded(shape)` (`cxx:248-255`).
    pub fn is_loaded(&self, shape: &TopoShape) -> bool {
        self.shells.contains(shape)
    }

    /// `NbLoaded()` (`cxx:258-261`).
    pub fn nb_loaded(&self) -> usize {
        self.shells.extent()
    }

    /// `Loaded(num)` (`cxx:264-267`).
    pub fn loaded(&self, num: usize) -> Option<&TopoShape> {
        self.shells.find_key(num)
    }

    /// `HasBadEdges()` (`cxx:272-275`).
    pub fn has_bad_edges(&self) -> bool {
        self.bad.extent() > 0
    }

    /// `BadEdges()` (`cxx:280-291`).
    pub fn bad_edges(&self) -> Compound {
        compound_from_seq(&self.bad.items)
    }

    /// `HasFreeEdges()` (`cxx:294-297`).
    pub fn has_free_edges(&self) -> bool {
        self.free.extent() > 0
    }

    /// `FreeEdges()` (`cxx:302-312`).
    pub fn free_edges(&self) -> Compound {
        compound_from_seq(&self.free.items)
    }

    /// `HasConnectedEdges()` (`cxx:317-320`).
    pub fn has_connected_edges(&self) -> bool {
        self.conex
    }
}

/// `CheckEdges(shape, bads, dirs, revs, ints)` (`ShapeAnalysis_Shell.cxx:72-135`).
fn check_edges(
    shape: &TopoShape,
    bads: &mut IndexedShapeMap,
    dirs: &mut IndexedShapeMap,
    revs: &mut IndexedShapeMap,
    ints: &mut IndexedShapeMap,
) -> bool {
    let mut res = false;
    if shape.shape_type() != ShapeType::Edge {
        for child in cumulated_children(shape) {
            if check_edges(&child, bads, dirs, revs, ints) {
                res = true;
            }
        }
        return res;
    }

    let edge = Edge(shape.clone());
    if BRepTool::is_degenerated(&edge) {
        return false; // `cxx:94-97`.
    }

    if shape.orientation() == Orientation::Forward {
        // `cxx:99-110`.
        if dirs.find_index(shape) == 0 {
            dirs.add(shape);
        } else {
            bads.add(shape);
            res = true;
        }
    }
    if shape.orientation() == Orientation::Reversed {
        // `cxx:111-123`.
        if revs.find_index(shape) == 0 {
            revs.add(shape);
        } else {
            bads.add(shape);
            res = true;
        }
    }
    if shape.orientation() == Orientation::Internal {
        // `cxx:124-131`: a repeated INTERNAL edge is deliberately not bad.
        if ints.find_index(shape) == 0 {
            ints.add(shape);
        }
    }
    res
}
