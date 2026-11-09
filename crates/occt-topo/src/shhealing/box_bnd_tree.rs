//! `ShapeAnalysis_BoxBndTreeSelector` (`ShapeAnalysis_BoxBndTree.cxx:28-196`,
//! `ShapeAnalysis_BoxBndTree.hxx:32-111`) and the two container classes it
//! selects over: `NCollection_UBTree` (`NCollection_UBTree.hxx:132-517`) and
//! `NCollection_UBTreeFiller` (`NCollection_UBTreeFiller.hxx:76-173`).
//!
//! The selector is what `ShapeAnalysis_FreeBounds::connectWiresToWiresImpl`
//! (`ShapeAnalysis_FreeBounds.cxx:235-299`) uses to find the next wire to
//! chain at the current wire's tail or head. Its answer is *order dependent*:
//! `NCollection_UBTree::Select` walks the leaves in the order the tree was
//! built (`hxx:480-517`), keeps the last accepted object for a head
//! connection (`myStop` stays false, `cxx:104-112`), and stops at the first
//! tail connection. A port that searched the wires in sequence order would
//! chain different wires, so the tree is reproduced exactly:
//!
//! * `NCollection_UBTree::Add` (`hxx:421-475`) inserts a leaf by descending
//!   the branch whose union box grows least, gemmating a leaf into a branch
//!   when the branch box is out or the node is a leaf (`hxx:179-199`).
//! * `NCollection_UBTreeFiller::Fill` (`hxx:147-173`) inserts the items in
//!   Fisher-Yates order driven by `std::mt19937` default-seeded with 5489
//!   (`hxx:80`). [`Mt19937`] is that generator, so the tree shape (and hence
//!   the leaf order, and hence the chaining) is identical.
//!
//! Port notes:
//!
//! * `NCollection_Array1<int> myArrIndices(1, 2)` (`hxx:100`) is a two-slot
//!   `[usize; 2]`: slot `0` is `First`, slot `1` is `Last` (`cxx:50-56`).
//! * `NCollection_Map<int> myList` (`hxx:106`) only answers `Contains`
//!   (`cxx:47`, `:86`); it is a `HashSet`.
//! * `myStop` (`NCollection_UBTree::Selector::myStop`, `hxx:122`) is not
//!   initialized by the OCCT selector constructor; `SetStop` (`hxx:77`)
//!   resets it before every `Select`. The port initializes it to `false` and
//!   keeps `set_stop`.
//! * `throw Standard_NoSuchObject` (`cxx:41-45`) has no counterpart; an index
//!   outside `1 ..= seq.Length()` answers `false`, which is a value a
//!   selection over a well formed tree can never need.
//! * `ShapeExtend::EncodeStatus` / `DecodeStatus` (`ShapeExtend.cxx:57-115`)
//!   are the bit constants [`SE_OK`] ... [`SE_FAIL`] and
//!   [`ShapeAnalysisBoxBndTreeSelector::last_check_status`].
//! * A null `myFVertex` / `myLVertex` makes `IsSame` answer false in OCCT, so
//!   the port's [`same_vertex`] takes an `Option`.
//! * `NCollection_UBTree::Select` calls `Reject` on every visited node box,
//!   leaves included (`hxx:496-505`); the port keeps that order.

use std::collections::HashSet;

use occt_core::bnd::BndBox;
use occt_core::gp::GpPnt;
use occt_core::precision::REAL_SMALL;

use crate::brep_tool::BRepTool;
use crate::shape::{Vertex, Wire};
use crate::topo_tools::is_same;

use super::free_bounds::find_bounds;

/// `ShapeExtend_OK` (`ShapeExtend.cxx:59-60`).
pub(in crate::shhealing) const SE_OK: i32 = 0x0000;
/// `ShapeExtend_DONE1` (`ShapeExtend.cxx:61-62`).
pub(in crate::shhealing) const SE_DONE1: i32 = 0x0001;
/// `ShapeExtend_DONE2` (`ShapeExtend.cxx:63-64`).
pub(in crate::shhealing) const SE_DONE2: i32 = 0x0002;
/// `ShapeExtend_DONE3` (`ShapeExtend.cxx:65-66`).
pub(in crate::shhealing) const SE_DONE3: i32 = 0x0004;
/// `ShapeExtend_DONE4` (`ShapeExtend.cxx:67-68`).
pub(in crate::shhealing) const SE_DONE4: i32 = 0x0008;
/// `ShapeExtend_FAIL1` (`ShapeExtend.cxx:74-75`).
pub(in crate::shhealing) const SE_FAIL1: i32 = 0x0100;
/// `ShapeExtend_FAIL2` (`ShapeExtend.cxx:76-77`).
pub(in crate::shhealing) const SE_FAIL2: i32 = 0x0200;
/// `ShapeExtend_FAIL`, the mask of every `FAILi` (`ShapeExtend.cxx:94-95`).
pub(in crate::shhealing) const SE_FAIL: i32 = 0xff00;

/// `ShapeExtend::DecodeStatus(flag, status)` (`ShapeExtend.cxx:107-115`):
/// `ShapeExtend_OK` asks for "no flag at all", every other status asks for its
/// bit.
pub(in crate::shhealing) fn decode_status(flag: i32, status: i32) -> bool {
    if status == SE_OK {
        flag == 0
    } else {
        flag & status != 0
    }
}

/// `TopoDS_Shape::IsSame` on a possibly null occurrence: OCCT's null
/// `TopoDS_Vertex::IsSame` answers false because a null shape has no
/// `TShape`.
fn same_vertex(a: Option<&Vertex>, b: &Vertex) -> bool {
    a.map(|v| is_same(&v.0, &b.0)).unwrap_or(false)
}

/// `std::mt19937` seeded with 5489, the generator
/// `NCollection_UBTreeFiller` constructs at `hxx:80-81`.
pub(in crate::shhealing) struct Mt19937 {
    state: [u32; 624],
    index: usize,
}

impl Default for Mt19937 {
    fn default() -> Self {
        Self::new()
    }
}

impl Mt19937 {
    /// `std::mt19937(5489u)`, the canonical MT19937 seeding
    /// (`mersenne_twister_engine::seed`), which leaves `index == 624` so the
    /// first draw twists the state.
    pub(in crate::shhealing) fn new() -> Self {
        let mut state = [0u32; 624];
        state[0] = 5489;
        for i in 1..624 {
            let prev = state[i - 1];
            state[i] = 1_812_433_253u32
                .wrapping_mul(prev ^ (prev >> 30))
                .wrapping_add(i as u32);
        }
        Self { state, index: 624 }
    }

    /// `mersenne_twister_engine::_M_gen_rand`.
    fn twist(&mut self) {
        for i in 0..624 {
            let y = (self.state[i] & 0x8000_0000) | (self.state[(i + 1) % 624] & 0x7fff_ffff);
            let mut next = self.state[(i + 397) % 624] ^ (y >> 1);
            if y & 1 != 0 {
                next ^= 0x9908_b0df;
            }
            self.state[i] = next;
        }
        self.index = 0;
    }

    /// `mersenne_twister_engine::operator()`: temper the next state word.
    pub(in crate::shhealing) fn next(&mut self) -> u32 {
        if self.index >= 624 {
            self.twist();
        }
        let mut y = self.state[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }
}

/// `NCollection_UBTree<int, Bnd_Box>::TreeNode` (`hxx:132-200`).
struct UbNode {
    bnd: BndBox,
    /// `myObject`; `None` is the nullified object of a gemmated branch
    /// (`hxx:198`).
    object: Option<usize>,
    /// `myChildren`; `None` is `IsLeaf()`.
    children: Option<[usize; 2]>,
    parent: Option<usize>,
}

impl UbNode {
    fn leaf(bnd: BndBox, object: usize) -> Self {
        Self {
            bnd,
            object: Some(object),
            children: None,
            parent: None,
        }
    }
}

/// `NCollection_UBTree::Selector` (`hxx:78-123`), the three methods the tree
/// calls. `Accept` receives `mySeq` because the OCCT array is mutated (a wire
/// is reversed at `ShapeAnalysis_FreeBounds.cxx:308`) while the selector holds
/// it.
pub(in crate::shhealing) trait UbSelector {
    fn reject(&self, the_bnd: &BndBox) -> bool;
    fn accept(&mut self, the_obj: usize, seq: &[Wire]) -> bool;
    fn stop(&self) -> bool;
}

/// `NCollection_UBTree<int, Bnd_Box>` (`hxx:64-418`): the `Add` / `Select`
/// pair over an arena, which stands in for the `TreeNode` pointers.
pub(in crate::shhealing) struct UbTree {
    nodes: Vec<UbNode>,
    root: Option<usize>,
}

impl Default for UbTree {
    fn default() -> Self {
        Self::new()
    }
}

impl UbTree {
    pub(in crate::shhealing) fn new() -> Self {
        Self {
            nodes: Vec::new(),
            root: None,
        }
    }

    fn push(&mut self, node: UbNode) -> usize {
        self.nodes.push(node);
        self.nodes.len() - 1
    }

    fn node(&self, i: usize) -> &UbNode {
        &self.nodes[i]
    }

    /// `NCollection_UBTree::Add(theObj, theBnd)` (`hxx:421-475`).
    pub(in crate::shhealing) fn add(&mut self, the_obj: usize, the_bnd: BndBox) {
        let Some(mut branch) = self.root else {
            // `hxx:425-431`: accepting the first object.
            self.root = Some(self.push(UbNode::leaf(the_bnd, the_obj)));
            return;
        };
        let mut is_out_of_branch = self.node(branch).bnd.is_out_box(&the_bnd);
        loop {
            // `hxx:438-447`: the walk stops when the branch box is out or the
            // node is a leaf, and the leaf is gemmated aside on this level.
            if is_out_of_branch || self.node(branch).children.is_none() {
                let mut new_bnd = the_bnd;
                new_bnd.add_box(&self.node(branch).bnd);
                let old_bnd = self.node(branch).bnd;
                let old_object = self.node(branch).object;
                let old_children = self.node(branch).children;
                // `hxx:184-197`: `children[0] = *this`, `children[1]` = new.
                let first = self.push(UbNode {
                    bnd: old_bnd,
                    object: old_object,
                    children: old_children,
                    parent: Some(branch),
                });
                let second = self.push(UbNode {
                    bnd: the_bnd,
                    object: Some(the_obj),
                    children: None,
                    parent: Some(branch),
                });
                if let Some([g0, g1]) = old_children {
                    // `hxx:191-195`: the moved children point at the copy.
                    self.nodes[g0].parent = Some(first);
                    self.nodes[g1].parent = Some(first);
                }
                self.nodes[branch].bnd = new_bnd;
                self.nodes[branch].object = None;
                self.nodes[branch].children = Some([first, second]);
                return;
            }

            // `hxx:449-450`: the branch box grows.
            self.nodes[branch].bnd.add_box(&the_bnd);

            // `hxx:452-468`: pick the only child that is not out, or else the
            // child whose union with `theBnd` is smallest.
            let [c0, c1] = self.node(branch).children.unwrap();
            let out0 = self.node(c0).bnd.is_out_box(&the_bnd);
            let out1 = self.node(c1).bnd.is_out_box(&the_bnd);
            let best = if out0 != out1 {
                if out0 { 1 } else { 0 }
            } else {
                let mut union0 = the_bnd;
                union0.add_box(&self.node(c0).bnd);
                let mut union1 = the_bnd;
                union1.add_box(&self.node(c1).bnd);
                if union0.square_extent() > union1.square_extent() {
                    1
                } else {
                    0
                }
            };
            is_out_of_branch = if best == 0 { out0 } else { out1 };
            branch = if best == 0 { c0 } else { c1 };
        }
    }

    /// `NCollection_UBTree::Select(theBranch, theSelector)` (`hxx:479-517`): a
    /// depth first walk that visits `child(0)` before `child(1)`, calls
    /// `Reject` on every visited node and `Accept` on every visited leaf, and
    /// stops as soon as the selector answers `Stop`.
    pub(in crate::shhealing) fn select(
        &self,
        selector: &mut dyn UbSelector,
        seq: &[Wire],
    ) -> usize {
        let Some(root) = self.root else {
            return 0;
        };
        let mut stack = vec![root];
        let mut n_sel = 0usize;
        while let Some(index) = stack.pop() {
            let node = &self.nodes[index];
            if selector.reject(&node.bnd) {
                continue;
            }
            match node.children {
                Some([c0, c1]) => {
                    // `hxx:511-513`: push child(1) first so child(0) pops
                    // first.
                    stack.push(c1);
                    stack.push(c0);
                }
                None => {
                    if let Some(obj) = node.object {
                        if selector.accept(obj, seq) {
                            n_sel += 1;
                        }
                    }
                    if selector.stop() {
                        break;
                    }
                }
            }
        }
        n_sel
    }
}

/// `NCollection_UBTreeFiller<int, Bnd_Box>` (`hxx:56-142`).
pub(in crate::shhealing) struct UbTreeFiller {
    items: Vec<(usize, BndBox)>,
    random: Mt19937,
    full_random: bool,
}

impl Default for UbTreeFiller {
    fn default() -> Self {
        Self::new()
    }
}

impl UbTreeFiller {
    /// `NCollection_UBTreeFiller(theTree, theAlloc, isFullRandom = true)`
    /// (`hxx:76-83`).
    pub(in crate::shhealing) fn new() -> Self {
        Self {
            items: Vec::new(),
            random: Mt19937::new(),
            full_random: true,
        }
    }

    /// `Add(theObj, theBnd)` (`hxx:86-89`).
    pub(in crate::shhealing) fn add(&mut self, the_obj: usize, the_bnd: BndBox) {
        self.items.push((the_obj, the_bnd));
    }

    /// `Fill()` (`hxx:147-173`): insert in Fisher-Yates order, drawing the
    /// index from `std::mt19937 % i` and moving the last item into the hole.
    pub(in crate::shhealing) fn fill(&mut self, tree: &mut UbTree) -> usize {
        let nb_add = self.items.len();
        for i in (1..=nb_add).rev() {
            let draw = self.random.next() % i as u32;
            let ind = if self.full_random {
                draw as usize
            } else {
                i - draw as usize - 1
            };
            let (obj, bnd) = self.items[ind];
            tree.add(obj, bnd);
            self.items[ind] = self.items[i - 1];
        }
        self.items.clear();
        nb_add
    }
}

/// `ShapeAnalysis_BoxBndTreeSelector` (`hxx:32-111`).
pub(in crate::shhealing) struct ShapeAnalysisBoxBndTreeSelector {
    /// `myFBox` / `myLBox`.
    f_box: BndBox,
    l_box: BndBox,
    /// `myShared`.
    shared: bool,
    /// `myNb`.
    nb: usize,
    /// `myFVertex` / `myLVertex`.
    f_vertex: Option<Vertex>,
    l_vertex: Option<Vertex>,
    /// `myFPnt` / `myLPnt`.
    f_pnt: GpPnt,
    l_pnt: GpPnt,
    /// `myList`.
    list: HashSet<usize>,
    /// `myTol`.
    tol: f64,
    /// `myMin3d`.
    min3d: f64,
    /// `myArrIndices(1, 2)`; slot `0` is `First`, slot `1` is `Last`.
    arr_indices: [usize; 2],
    /// `myStatus`.
    status: i32,
    /// `myStop`.
    stop: bool,
}

impl ShapeAnalysisBoxBndTreeSelector {
    /// `ShapeAnalysis_BoxBndTreeSelector(theSeq, theShared)` (`hxx:36-47`).
    pub(in crate::shhealing) fn new(shared: bool) -> Self {
        Self {
            f_box: BndBox::new(),
            l_box: BndBox::new(),
            shared,
            nb: 0,
            f_vertex: None,
            l_vertex: None,
            f_pnt: GpPnt::new(0.0, 0.0, 0.0),
            l_pnt: GpPnt::new(0.0, 0.0, 0.0),
            list: HashSet::new(),
            tol: 1e-7,
            min3d: 1e-7,
            arr_indices: [0, 0],
            status: SE_OK,
            stop: false,
        }
    }

    /// `DefineBoxes(theFBox, theLBox)` (`hxx:49-55`).
    pub(in crate::shhealing) fn define_boxes(&mut self, the_f_box: BndBox, the_l_box: BndBox) {
        self.f_box = the_f_box;
        self.l_box = the_l_box;
        self.arr_indices = [0, 0];
    }

    /// `DefineVertexes(theVf, theVl)` (`hxx:57-63`).
    pub(in crate::shhealing) fn define_vertexes(&mut self, vf: Option<Vertex>, vl: Option<Vertex>) {
        self.f_vertex = vf;
        self.l_vertex = vl;
        self.status = SE_OK;
    }

    /// `DefinePnt(theFPnt, theLPnt)` (`hxx:65-69`).
    pub(in crate::shhealing) fn define_pnt(&mut self, the_f_pnt: GpPnt, the_l_pnt: GpPnt) {
        self.f_pnt = the_f_pnt;
        self.l_pnt = the_l_pnt;
        self.status = SE_OK;
    }

    /// `GetNb()` (`hxx:71`).
    pub(in crate::shhealing) fn get_nb(&self) -> usize {
        self.nb
    }

    /// `SetNb(theNb)` (`hxx:73`).
    pub(in crate::shhealing) fn set_nb(&mut self, the_nb: usize) {
        self.nb = the_nb;
    }

    /// `LoadList(elem)` (`hxx:75`): the map is only ever added to, never
    /// cleared (`ShapeAnalysis_FreeBounds.cxx:298`, `:407`).
    pub(in crate::shhealing) fn load_list(&mut self, elem: usize) {
        self.list.insert(elem);
    }

    /// `SetStop()` (`hxx:77`).
    pub(in crate::shhealing) fn set_stop(&mut self) {
        self.stop = false;
    }

    /// `SetTolerance(theTol)` (`hxx:79-84`).
    pub(in crate::shhealing) fn set_tolerance(&mut self, the_tol: f64) {
        self.tol = the_tol;
        self.min3d = the_tol;
        self.status = SE_OK;
    }

    /// `ContWire(nbWire)` (`hxx:86`).
    pub(in crate::shhealing) fn cont_wire(&self, nb_wire: usize) -> bool {
        self.list.contains(&nb_wire)
    }

    /// `LastCheckStatus(theStatus)` (`hxx:88-91`).
    pub(in crate::shhealing) fn last_check_status(&self, the_status: i32) -> bool {
        decode_status(self.status, the_status)
    }
}

impl UbSelector for ShapeAnalysisBoxBndTreeSelector {
    /// `Reject(theBnd)` (`cxx:28-34`): the object is skipped only when its box
    /// is out of *both* endpoint boxes.
    fn reject(&self, the_bnd: &BndBox) -> bool {
        let fch = self.f_box.is_out_box(the_bnd);
        let lch = self.l_box.is_out_box(the_bnd);
        fch && lch
    }

    /// `Accept(theObj)` (`cxx:37-195`). The `myShared` arm (`cxx:61-112`) is
    /// the one `ShapeAnalysis_FreeBounds` reaches: its shell constructor chains
    /// with `shared = true` (`ShapeAnalysis_FreeBounds.cxx:126`) and
    /// `SplitWire` inherits the same flag (`:586`). The point-distance arm
    /// (`cxx:114-192`) belongs to the sewing constructor's `shared = false`
    /// call (`:91`) and is ported for completeness.
    fn accept(&mut self, the_obj: usize, seq: &[Wire]) -> bool {
        if the_obj < 1 || the_obj > seq.len() {
            // `cxx:41-45`: `Standard_NoSuchObject`.
            return false;
        }
        if self.list.contains(&the_obj) {
            return false; // `cxx:47-50`.
        }

        let (v1, v2) = find_bounds(&seq[the_obj - 1].0);
        let (Some(v1), Some(v2)) = (v1, v2) else {
            // A wire without two endpoint vertices cannot connect; OCCT would
            // dereference a null shape in `BRep_Tool::Pnt` (`cxx:60`, `:117`).
            self.status = SE_FAIL2;
            self.stop = false;
            return false;
        };

        if self.shared {
            // `cxx:61-112`.
            let is_accept;
            if same_vertex(self.l_vertex.as_ref(), &v1) {
                self.status = SE_DONE1;
                self.arr_indices[1] = the_obj;
                is_accept = true;
            } else if same_vertex(self.l_vertex.as_ref(), &v2) {
                self.status = SE_DONE2;
                self.arr_indices[1] = the_obj;
                is_accept = true;
            } else if same_vertex(self.f_vertex.as_ref(), &v2) {
                self.status = SE_DONE3;
                self.arr_indices[0] = the_obj;
                is_accept = true;
            } else if same_vertex(self.f_vertex.as_ref(), &v1) {
                self.status = SE_DONE4;
                self.arr_indices[0] = the_obj;
                is_accept = true;
            } else {
                self.status = SE_FAIL2;
                is_accept = false;
            }

            if is_accept {
                // `cxx:101-108`.
                self.set_nb(the_obj);
                if self.arr_indices[1] != 0 {
                    self.stop = true;
                }
                return true;
            }
            self.stop = false; // `cxx:110-111`.
            return false; // `cxx:194`.
        }

        // `cxx:114-192`: the distance arm used by the `shared = false` call of
        // the sewing constructor.
        let p1 = BRepTool::vertex_point(&v1);
        let p2 = BRepTool::vertex_point(&v2);
        let tail_head = p1.distance(&self.l_pnt);
        let tail_tail = p2.distance(&self.l_pnt);
        let head_head = p1.distance(&self.f_pnt);
        let head_tail = p2.distance(&self.f_pnt);
        let mut dm1 = tail_head;
        let mut dm2 = head_tail;
        let mut res1 = 0;
        let mut res2 = 0;
        if tail_head > tail_tail {
            res1 = 1;
            dm1 = tail_tail;
        }
        if head_tail > head_head {
            res2 = 1;
            dm2 = head_head;
        }
        let mut result = res1;
        let min3d = dm1.min(dm2);
        if min3d > self.min3d {
            return false; // `cxx:140-143`.
        }

        let min_ind = if dm1 > dm2 { 0 } else { 1 };
        let max_ind = if dm1 > dm2 { 1 } else { 0 };
        self.arr_indices[min_ind] = the_obj;
        if (min3d - self.min3d) > REAL_SMALL {
            self.arr_indices[max_ind] = 0;
        }
        self.min3d = min3d; // `cxx:143-155`.
        if min3d > self.tol {
            self.status = SE_FAIL2;
            return false; // `cxx:156-159`.
        }

        let an_obj = if self.arr_indices[1] != 0 {
            self.arr_indices[1]
        } else {
            self.arr_indices[0]
        };
        self.set_nb(an_obj); // `cxx:160-161`.

        if min3d == 0.0 && min_ind == 1 {
            self.stop = true; // `cxx:163-166`.
        }

        if dm1 > dm2 {
            dm1 = dm2;
            result = res2 + 2;
        }
        let _ = dm1;
        if an_obj == the_obj {
            self.status = match result {
                0 => SE_DONE1,
                1 => SE_DONE2,
                2 => SE_DONE3,
                _ => SE_DONE4,
            };
        }
        true
    }

    fn stop(&self) -> bool {
        self.stop
    }
}
