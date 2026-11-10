//! `TopOpeBRepBuild_BlockBuilder` — the connexity blocks of an element set.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_BlockBuilder.hxx:30-72`,
//! `.cxx:22-284`.
//!
//! Logic (`cxx:30-48`):
//! - a block is a set of connex elements of `SS`, and every element of `SS`
//!   appears in only one connexity block;
//! - all the elements of all the blocks are stored in `myOrientedShapeMap`;
//! - a connexity block is a segment `[f, l]` of element indices of that map;
//! - `myBlocks` holds the index of the first element of each block, so
//!   `f = myBlocks(i)` and `l = myBlocks(i+1) - 1`.
//!
//! `NCollection_Sequence<int>` is 1-based (`Value(i)`), hence the `- 1` shifts
//! below.

use std::collections::HashMap;

use crate::shape::TopoShape;
use crate::topopebrep::ds::data_structure::IndexedShapeSet;

use super::block_iterator::BlockIterator;
use super::shape_set::ShapeSetSource;

/// `class TopOpeBRepBuild_BlockBuilder` (`hxx:40`).
pub struct BlockBuilder {
    /// `myOrientedShapeMapIsValid` (`hxx:64`).
    my_oriented_shape_map_is_valid: HashMap<i32, i32>,
    /// `myOrientedShapeMap` (`hxx:65`).
    my_oriented_shape_map: IndexedShapeSet,
    /// `myBlocks` (`hxx:66`), 1-based.
    my_blocks: Vec<i32>,
    /// `myBlockIndex` (`hxx:67`).
    my_block_index: i32,
    /// `myIsDone` (`hxx:68`).
    my_is_done: bool,
    /// `myBlocksIsRegular` (`hxx:69`), 1-based.
    my_blocks_is_regular: Vec<i32>,
}

impl BlockBuilder {
    /// `TopOpeBRepBuild_BlockBuilder()` (`cxx:26-29`).
    pub fn new() -> Self {
        Self {
            my_oriented_shape_map_is_valid: HashMap::new(),
            my_oriented_shape_map: IndexedShapeSet::new(),
            my_blocks: Vec::new(),
            my_block_index: 0,
            my_is_done: false,
            my_blocks_is_regular: Vec::new(),
        }
    }

    /// `TopOpeBRepBuild_BlockBuilder(SS)` (`cxx:33-37`).
    pub fn with_shape_set(ss: &mut dyn ShapeSetSource) -> Self {
        let mut b = Self::new();
        b.make_block(ss);
        b
    }

    /// `MakeBlock(SS)` (`cxx:41-146`).
    ///
    /// `SS` is taken by reference and driven through its virtual members, as
    /// the C++ signature `MakeBlock(TopOpeBRepBuild_ShapeSet& SS)` does.
    pub fn make_block(&mut self, ss: &mut dyn ShapeSetSource) {
        self.my_oriented_shape_map.clear();
        self.my_oriented_shape_map_is_valid.clear();
        self.my_blocks.clear();
        self.my_blocks_is_regular.clear();

        let mut is_regular;
        let mut cur_nei;
        let mut m_extent;
        let mut e_index;

        ss.init_start_elements();
        while ss.more_start_elements() {
            let e = match ss.start_element() {
                Some(s) => s.clone(),
                None => break,
            };
            m_extent = self.my_oriented_shape_map.extent();
            e_index = self.add_element(&e);

            // E = current element of the element set SS
            // Mextent = index of last element stored in map M, before E is
            // added. Eindex = index of E added to M: Eindex > Mextent means E
            // is new in M.
            let e_new_in_m = e_index > m_extent;
            if e_new_in_m {
                // make a new block starting at element Eindex
                self.my_blocks.push(e_index);
                is_regular = true;
                // `_cxx:83` also sets `CurNei = 0`, which `_cxx:98` overwrites
                // before any read: an OCCT dead store.

                m_extent = self.my_oriented_shape_map.extent();
                let mut search_neighbours = e_index <= m_extent;
                while search_neighbours {
                    // E1 = element of M on which neighbours must be searched
                    let e1 = match self.my_oriented_shape_map.find_key(e_index) {
                        Some(s) => s.clone(),
                        None => break,
                    };
                    cur_nei = ss.max_number_sub_shape(&e1);
                    let mut cond_regu = true;
                    if cur_nei > 2 {
                        cond_regu = false;
                    }
                    is_regular = is_regular && cond_regu;
                    // compute neighbours of E1: add them to M to increase
                    // M.Extent().
                    ss.init_neighbours(&e1);

                    while ss.more_neighbours() {
                        let n = match ss.neighbour() {
                            Some(s) => s.clone(),
                            None => break,
                        };
                        self.add_element(&n);
                        ss.next_neighbour();
                    }

                    e_index += 1;
                    m_extent = self.my_oriented_shape_map.extent();
                    search_neighbours = e_index <= m_extent;
                }

                let iiregu = if is_regular { 1 } else { 0 };
                self.my_blocks_is_regular.push(iiregu);
            }

            ss.next_start_element();
        }

        // To value the l bound of the last connexity block created above, we
        // create an artificial block of value = myOrientedShapeMap.Extent() + 1.
        // The real number of connexity blocks is myBlocks.Length() - 1.
        m_extent = self.my_oriented_shape_map.extent();
        self.my_blocks.push(m_extent + 1);
        self.my_is_done = true;
    }

    /// `InitBlock()` (`cxx:150-153`).
    pub fn init_block(&mut self) {
        self.my_block_index = 1;
    }

    /// `MoreBlock()` (`cxx:157-163`): the length of `myBlocks` is one more
    /// than the number of connexity blocks.
    pub fn more_block(&self) -> bool {
        let l = self.my_blocks.len() as i32;
        self.my_block_index < l
    }

    /// `NextBlock()` (`cxx:167-170`).
    pub fn next_block(&mut self) {
        self.my_block_index += 1;
    }

    /// `BlockIterator()` (`cxx:174-180`).
    ///
    /// `myBlocks(myBlockIndex)` is 1-based, so `myBlockIndex + 1` is the
    /// `myBlockIndex`-th slot of the zero-based `Vec`.
    pub fn block_iterator(&self) -> BlockIterator {
        let lower = self
            .my_blocks
            .get((self.my_block_index - 1) as usize)
            .copied()
            .unwrap_or(0);
        let upper = self
            .my_blocks
            .get(self.my_block_index as usize)
            .copied()
            .unwrap_or(1)
            - 1;
        BlockIterator::with_bounds(lower, upper)
    }

    /// `Element(BI)` (`cxx:184-196`).
    pub fn element_of_iterator(&self, bi: &BlockIterator) -> Option<&TopoShape> {
        if !bi.more() {
            return None;
        }
        self.my_oriented_shape_map.find_key(bi.value())
    }

    /// `Element(int index)` (`cxx:198-208`).
    pub fn element_at(&self, index: i32) -> Option<&TopoShape> {
        if !self.my_oriented_shape_map_is_valid.contains_key(&index) {
            return None;
        }
        self.my_oriented_shape_map.find_key(index)
    }

    /// `Element(const TopoDS_Shape& E)` (`cxx:210-220`).
    pub fn element_index(&self, e: &TopoShape) -> i32 {
        if !self.my_oriented_shape_map.contains(e) {
            return 0;
        }
        self.my_oriented_shape_map.find_index(e)
    }

    /// `ElementIsValid(BI)` (`cxx:224-238`).
    pub fn element_is_valid_of_iterator(&self, bi: &BlockIterator) -> bool {
        if !bi.more() {
            return false;
        }
        let s_index = bi.value();
        let isb = *self
            .my_oriented_shape_map_is_valid
            .get(&s_index)
            .unwrap_or(&0);
        isb == 1
    }

    /// `ElementIsValid(int Sindex)` (`cxx:240-250`).
    pub fn element_is_valid(&self, s_index: i32) -> bool {
        if !self.my_oriented_shape_map_is_valid.contains_key(&s_index) {
            return false;
        }
        let isb = *self
            .my_oriented_shape_map_is_valid
            .get(&s_index)
            .unwrap_or(&0);
        isb == 1
    }

    /// `AddElement(S)` (`cxx:254-260`).
    pub fn add_element(&mut self, s: &TopoShape) -> i32 {
        let s_index = self.my_oriented_shape_map.add(s);
        self.my_oriented_shape_map_is_valid.insert(s_index, 1);
        s_index
    }

    /// `SetValid(BI, isvalid)` (`cxx:264-276`).
    pub fn set_valid_of_iterator(&mut self, bi: &BlockIterator, isvalid: bool) {
        if !bi.more() {
            return;
        }
        let s_index = bi.value();
        let i = if isvalid { 1 } else { 0 };
        self.my_oriented_shape_map_is_valid.insert(s_index, i);
    }

    /// `SetValid(int Sindex, isvalid)` (`cxx:278-288`).
    pub fn set_valid(&mut self, s_index: i32, isvalid: bool) {
        if !self.my_oriented_shape_map_is_valid.contains_key(&s_index) {
            return;
        }
        let i = if isvalid { 1 } else { 0 };
        self.my_oriented_shape_map_is_valid.insert(s_index, i);
    }

    /// `CurrentBlockIsRegular()` (`cxx:292-301`).
    pub fn current_block_is_regular(&self) -> bool {
        let i = *self
            .my_blocks_is_regular
            .get((self.my_block_index - 1).max(0) as usize)
            .unwrap_or(&0);
        i == 1
    }
}

impl Default for BlockBuilder {
    fn default() -> Self {
        Self::new()
    }
}
