//! Make loops from a set of connected links: closes an unordered set of edges
//! into closed loops, using graph cycle search and handling hanging chains.
//! Source: `Poly_MakeLoops.hxx` / `Poly_MakeLoops.cxx`.

use std::collections::{BTreeSet, HashMap};

use crate::gp::{GpDir, GpDir2d, GpVec};

/// Orientation flags that can be attached to a link.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LinkFlag(pub u8);

pub const LF_NONE: LinkFlag = LinkFlag(0);
pub const LF_FWD: LinkFlag = LinkFlag(1);  // forward orientation
pub const LF_REV: LinkFlag = LinkFlag(2);  // reversed orientation
pub const LF_BOTH: LinkFlag = LinkFlag(3); // both ways oriented
pub const LF_REVERSED: LinkFlag = LinkFlag(4); // the link is reversed

/// A link between two nodes, represented by a pair of node indices.
/// Node indices are 1-based (0 denotes a null link), matching OCCT.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Link {
    pub node1: usize,
    pub node2: usize,
    pub flags: u8,
}

impl Link {
    pub fn new(node1: usize, node2: usize) -> Self { Self { node1, node2, flags: LF_FWD.0 } }
    pub fn reverse(&mut self) { self.flags ^= LF_REVERSED.0; }
    pub fn is_reversed(&self) -> bool { self.flags & LF_REVERSED.0 != 0 }
    pub fn nullify(&mut self) { self.node1 = 0; self.node2 = 0; }
    pub fn is_null(&self) -> bool { self.node1 == 0 || self.node2 == 0 }
    /// Orientation-insensitive equality (undirected pair).
    pub fn same_pair(&self, other: &Link) -> bool {
        (other.node1 == self.node1 && other.node2 == self.node2)
            || (other.node1 == self.node2 && other.node2 == self.node1)
    }
}

impl Default for Link {
    fn default() -> Self { Self { node1: 0, node2: 0, flags: 0 } }
}

/// Result codes from `Perform`.
pub const RC_LOOPS_DONE: u8 = 1;
pub const RC_HANGING_LINKS: u8 = 2;
pub const RC_FAILURE: u8 = 4;

/// Helper providing adjacent links for each node, used by the loop search.
pub trait MakeLoopsHelper {
    /// Returns the links adjacent to the given node.
    fn get_adjacent_links(&self, node: usize) -> &[Link];
    /// Optional hook called from `add_link`.
    fn on_add_link(&self, _num: usize, _link: &Link) {}
    /// `Poly_MakeLoops::Helper::GetNormal` (`Poly_MakeLoops.hxx:104-118`).
    /// OCCT's base implementation returns `false`, which makes `chooseLeftWay`
    /// fall back to `theLstIndS.First()`; `None` is that `false`.
    fn get_normal(&self, _node: usize) -> Option<GpDir> {
        None
    }
    /// `Poly_MakeLoops3D::Helper::GetFirstTangent` (`Poly_MakeLoops.hxx:276`).
    fn get_first_tangent(&self, _link: &Link) -> Option<GpDir> {
        None
    }
    /// `Poly_MakeLoops3D::Helper::GetLastTangent` (`Poly_MakeLoops.hxx:279`).
    fn get_last_tangent(&self, _link: &Link) -> Option<GpDir> {
        None
    }
    /// `Poly_MakeLoops2D::Helper::GetFirstTangent` (`Poly_MakeLoops.hxx:317`).
    fn get_first_tangent_2d(&self, _link: &Link) -> Option<GpDir2d> {
        None
    }
    /// `Poly_MakeLoops2D::Helper::GetLastTangent` (`Poly_MakeLoops.hxx:320`).
    fn get_last_tangent_2d(&self, _link: &Link) -> Option<GpDir2d> {
        None
    }
}

/// Make loops from a set of connected links.
///
/// A link is represented by a pair of node indices (1-based). The algorithm
/// collects the links, then repeatedly walks from a starting link following
/// the adjacency provided by the helper until a closed loop is formed.
///
/// Node indices within `Link` are **1-based** (0 = null) to match OCCT
/// exactly; callers must pass 1-based node ids to `AddLink`.
pub struct PolyMakeLoops<'a> {
    helper: &'a dyn MakeLoopsHelper,
    /// Stored links; `link_index` maps an unordered pair to a 1-based index.
    map_links: Vec<Link>,
    link_index: HashMap<(usize, usize), usize>,
    loops: Vec<Vec<Link>>,
    start_indices: BTreeSet<isize>,
    hang_indices: BTreeSet<isize>,
    /// `Poly_MakeLoops2D` flag: `myRightWay = !theLeftWay`
    /// (`Poly_MakeLoops.cxx:682-687`). The 3D class has no such flag; the port
    /// collapses `Poly_MakeLoops3D` and `Poly_MakeLoops2D` into one struct, and
    /// this boolean selects which `chooseLeftWay` runs.
    two_d: bool,
    right_way: bool,
}

impl<'a> PolyMakeLoops<'a> {
    /// Constructor of the 3D flavour (`Poly_MakeLoops3D`). The helper supplies
    /// adjacency; a null helper would make the algorithm return a wrong result.
    pub fn new(helper: &'a dyn MakeLoopsHelper) -> Self {
        Self {
            helper,
            map_links: Vec::new(),
            link_index: HashMap::new(),
            loops: Vec::new(),
            start_indices: BTreeSet::new(),
            hang_indices: BTreeSet::new(),
            two_d: false,
            right_way: false,
        }
    }

    /// Constructor of the 2D flavour: `Poly_MakeLoops2D(theLeftWay, helper,
    /// alloc)` stores `myRightWay = !theLeftWay` (`Poly_MakeLoops.cxx:682-687`).
    pub fn new_2d(helper: &'a dyn MakeLoopsHelper, left_way: bool) -> Self {
        let mut ml = Self::new(helper);
        ml.two_d = true;
        ml.right_way = !left_way;
        ml
    }

    /// Resets the algorithm to its initial state.
    pub fn reset(&mut self, helper: &'a dyn MakeLoopsHelper) {
        self.helper = helper;
        self.map_links.clear();
        self.link_index.clear();
        self.loops.clear();
        self.start_indices.clear();
        self.hang_indices.clear();
    }

    /// Adds a link to the set. The orientation flags of a link are OR-ed when
    /// the same undirected pair is added twice.
    pub fn add_link(&mut self, link: &Link) {
        if link.node1 == link.node2 || link.node1 == 0 || link.node2 == 0 {
            return;
        }
        let key = if link.node1 < link.node2 { (link.node1, link.node2) } else { (link.node2, link.node1) };
        let index = *self.link_index.get(&key).unwrap_or(&0);
        if index == 0 {
            self.link_index.insert(key, self.map_links.len() + 1);
            self.map_links.push(*link);
        } else {
            let stored = &mut self.map_links[index - 1];
            stored.flags |= link.flags;
            stored.node1 = link.node1;
            stored.node2 = link.node2;
        }
        self.helper.on_add_link(self.map_links.len(), link);
    }

    /// Replaces one link with another (e.g. to change the order of nodes).
    pub fn replace_link(&mut self, link: &Link, new_link: &Link) {
        if new_link.node1 == new_link.node2 || new_link.node1 == 0 || new_link.node2 == 0 {
            return;
        }
        if let Some(index) = self.find_link_index(link) {
            self.map_links[index - 1] = *new_link;
            // Rebuild the index for the replaced pair.
            let key = if new_link.node1 < new_link.node2 { (new_link.node1, new_link.node2) } else { (new_link.node2, new_link.node1) };
            self.link_index.insert(key, index);
            self.helper.on_add_link(index, new_link);
        }
    }

    /// Sets a new orientation flag for a previously added link. Returns the
    /// old orientation (the `LF_BOTH` part).
    pub fn set_link_orientation(&mut self, link: &Link, orient: LinkFlag) -> LinkFlag {
        if let Some(index) = self.find_link_index(link) {
            let stored = &mut self.map_links[index - 1];
            let old = LinkFlag(stored.flags & LF_BOTH.0);
            stored.flags = orient.0;
            self.helper.on_add_link(index, &self.map_links[index - 1]);
            old
        } else {
            LF_NONE
        }
    }

    /// Finds the stored link by value (orientation-insensitive).
    pub fn find_link(&self, link: &Link) -> Link {
        match self.find_link_index(link) {
            Some(i) => self.map_links[i - 1],
            None => Link::default(),
        }
    }

    /// Returns the 1-based index of the stored link, or `None`.
    fn find_link_index(&self, link: &Link) -> Option<usize> {
        if link.node1 == 0 || link.node2 == 0 { return None; }
        let key = if link.node1 < link.node2 { (link.node1, link.node2) } else { (link.node2, link.node1) };
        self.link_index.get(&key).copied()
    }

    /// Runs the loop-building algorithm. Returns the bitmask of result codes,
    /// or an error on internal failure.
    pub fn perform(&mut self) -> Result<u8, String> {
        // Prepare the set of start indices.
        self.start_indices.clear();
        for (i, link) in self.map_links.iter().enumerate() {
            let idx = (i + 1) as isize;
            if link.flags & LF_FWD.0 != 0 { self.start_indices.insert(idx); }
            if link.flags & LF_REV.0 != 0 { self.start_indices.insert(-idx); }
        }

        let mut result: u8 = 0;
        // Two-pass loop: hanging links found on the first pass are retried on
        // the second pass.
        for pass in 0..2 {
            self.hang_indices.clear();
            while !self.start_indices.is_empty() {
                let index_s = *self.start_indices.iter().next().unwrap();
                let (contour, start_number) = self.find_contour(index_s);
                if start_number == 0 {
                    return Err("MakeLoops: internal failure".to_string());
                }
                if start_number <= contour.len() {
                    // There is a closed loop in the contour.
                    self.accept_contour(&contour, start_number);
                }
                if start_number > 1 {
                    // It is required to mark hanging edges.
                    let node = if start_number <= contour.len() {
                        self.get_first_node(index_s)
                    } else {
                        // Open contour: mark from the end back to a bifurcation.
                        let last_s = contour[start_number - 2];
                        self.get_last_node(last_s)
                    };
                    self.mark_hang_chain(node, index_s);
                }
            }
            if pass == 0 {
                // Move hanging links to start indices for the second pass.
                let hang: Vec<isize> = self.hang_indices.iter().cloned().collect();
                for h in hang { self.start_indices.insert(h); }
            }
        }

        if !self.loops.is_empty() { result |= RC_LOOPS_DONE; }
        if !self.hang_indices.is_empty() { result |= RC_HANGING_LINKS; }
        Ok(result)
    }

    /// Returns the number of loops in the result.
    pub fn get_nb_loops(&self) -> usize { self.loops.len() }

    /// Returns the loop of the given index.
    pub fn get_loop(&self, index: usize) -> &[Link] { &self.loops[index] }

    /// Returns the number of detected hanging chains.
    pub fn get_nb_hanging(&self) -> usize { self.hang_indices.len() }

    /// Fills in the list of hanging links.
    pub fn get_hanging_links(&self) -> Vec<Link> {
        let mut out = Vec::new();
        for &index_s in &self.hang_indices {
            let mut link = self.map_links[index_s.unsigned_abs() as usize - 1];
            if index_s < 0 { link.reverse(); }
            out.push(link);
        }
        out
    }

    /// Chooses the next link at a branching node.
    ///
    /// `Poly_MakeLoops3D::chooseLeftWay` (`Poly_MakeLoops.cxx:614-678`): the
    /// incoming link's tangent is projected into the plane orthogonal to the
    /// node normal to give the reference direction, and the candidate whose
    /// (projected) first tangent makes the smallest signed
    /// `gp_Dir::AngleWithRef` with it wins. Every step that OCCT guards with a
    /// failed helper accessor falls back to `theLstIndS.First()` (`:621-624`,
    /// `:628-631`, `:636-640`), as does an all-rejected candidate list
    /// (`:677`).
    ///
    /// The 2D flavour (`new_2d`, `Poly_MakeLoops2D`) runs
    /// [`Self::choose_left_way_2d`] instead (`cxx:692-738`).
    pub fn choose_left_way(&self, node: usize, seg_index: isize, lst_ind_s: &[isize]) -> isize {
        if self.two_d {
            return self.choose_left_way_2d(seg_index, lst_ind_s);
        }
        let mut angle_min = std::f64::consts::PI * 2.0;
        let Some(normal) = self.helper.get_normal(node) else {
            return lst_ind_s[0];
        };
        let link = &self.map_links[seg_index.unsigned_abs() as usize - 1];
        let Some(tgt_ref) = self.helper.get_last_tangent(link) else {
            return lst_ind_s[0];
        };

        // Project the tangent onto the plane orthogonal to the normal to get the
        // reference direction: `Normal.CrossCrossed(TgtRef, Normal)` (`cxx:635`).
        let n = GpVec::from_xyz(normal.xyz());
        let t = GpVec::from_xyz(tgt_ref.xyz());
        let tgt_ref_xyz = n.cross_crossed(&t);
        if tgt_ref_xyz.square_magnitude() < 1e-14 {
            // A problem with defining the reference direction: take the first way.
            return lst_ind_s[0];
        }
        let Ok(tgt_ref) = GpDir::from_xyz(tgt_ref_xyz.xyz()) else {
            return lst_ind_s[0];
        };

        // Find the way with the minimal angle to the reference direction (the
        // angle is in the range ]-PI; PI]).
        let mut res_index: isize = 0;
        for &ind_s in lst_ind_s {
            let link = &self.map_links[ind_s.unsigned_abs() as usize - 1];
            let Some(tgt) = self.helper.get_first_tangent(link) else {
                continue;
            };
            let t = GpVec::from_xyz(tgt.xyz());
            let tgt_xyz = n.cross_crossed(&t);
            if tgt_xyz.square_magnitude() < 1e-14 {
                // Skip a problem way.
                continue;
            }
            let Ok(tgt) = GpDir::from_xyz(tgt_xyz.xyz()) else {
                continue;
            };
            let mut angle = tgt.angle_with_ref(&tgt_ref, &normal);
            if angle < 1e-4 - std::f64::consts::PI {
                angle = std::f64::consts::PI;
            }
            if angle < angle_min {
                angle_min = angle;
                res_index = ind_s;
            }
        }
        if res_index == 0 {
            lst_ind_s[0]
        } else {
            res_index
        }
    }

    /// `Poly_MakeLoops2D::chooseLeftWay` (`Poly_MakeLoops.cxx:692-738`): the same
    /// minimal-signed-angle selection as the 3D flavour, but using the links'
    /// 2D first/last tangents directly (no normal projection) and negating every
    /// candidate angle when `myRightWay` is set (`theLeftWay == false`,
    /// `cxx:682-687`). A failed `GetLastTangent` on the incoming link falls back
    /// to `theLstIndS.First()` (`cxx:700-704`), failing candidates are skipped
    /// (`cxx:716-720`), and an all-rejected list also falls back (`cxx:737`).
    fn choose_left_way_2d(&self, seg_index: isize, lst_ind_s: &[isize]) -> isize {
        let mut angle_min = std::f64::consts::PI * 2.0;
        let link = &self.map_links[seg_index.unsigned_abs() as usize - 1];
        let Some(tgt_ref) = self.helper.get_last_tangent_2d(link) else {
            // A problem with defining the reference direction: take the first way.
            return lst_ind_s[0];
        };
        let mut res_index: isize = 0;
        for &ind_s in lst_ind_s {
            let link = &self.map_links[ind_s.unsigned_abs() as usize - 1];
            let Some(tgt) = self.helper.get_first_tangent_2d(link) else {
                // Skip a problem way.
                continue;
            };
            let mut angle = tgt.angle(&tgt_ref);
            if self.right_way {
                angle = -angle;
            }
            if angle < 1e-4 - std::f64::consts::PI {
                angle = std::f64::consts::PI;
            }
            if angle < angle_min {
                angle_min = angle;
                res_index = ind_s;
            }
        }
        if res_index == 0 {
            lst_ind_s[0]
        } else {
            res_index
        }
    }

    /// Collects edges in a chain until they form a closed contour. Returns the
    /// 1-based start index within `contour` where the loop begins; a value of
    /// `contour.len() + 1` means the contour is open.
    fn find_contour(&self, start_index_s: isize) -> (Vec<isize>, usize) {
        let mut contour: Vec<isize> = Vec::new();
        let mut node_link: HashMap<usize, isize> = HashMap::new();
        let mut index_s = start_index_s;
        let mut last_node = self.get_last_node(index_s);
        let start_number;
        loop {
            contour.push(index_s);
            node_link.insert(self.get_first_node(index_s), index_s);
            let index = index_s.unsigned_abs() as usize;

            // Collect the list of links from this node able to participate.
            let mut lst_ind_s: Vec<isize> = Vec::new();
            for link in self.helper.get_adjacent_links(last_node) {
                let Some(ind) = self.find_link_index(link) else { continue };
                if ind == index { continue; }
                let mut ind_s = ind as isize;
                if self.get_first_node(ind as isize) != last_node {
                    ind_s = -ind_s;
                }
                if self.can_link_be_taken(ind_s) {
                    lst_ind_s.push(ind_s);
                }
            }

            if lst_ind_s.is_empty() {
                // No more ways: open contour.
                start_number = contour.len() + 1;
                break;
            }

            let index_s_next = if lst_ind_s.len() == 1 {
                lst_ind_s[0]
            } else {
                self.choose_left_way(last_node, index_s, &lst_ind_s)
            };
            index_s = index_s_next;

            if index_s == 0 {
                // No more ways: open contour.
                start_number = contour.len() + 1;
                break;
            }
            if let Some(pos) = contour.iter().position(|&x| x == index_s) {
                // Entering the loop a second time: stop search.
                start_number = pos + 1;
                break;
            }
            if let Some(pos) = contour.iter().position(|&x| x == -index_s) {
                // Leaving the loop: stop search.
                start_number = pos + 2;
                break;
            }

            last_node = self.get_last_node(index_s);

            if let Some(&bound) = node_link.get(&last_node) {
                // Closing the loop: stop search.
                contour.push(index_s);
                start_number = contour.iter().position(|&x| x == bound).unwrap() + 1;
                break;
            }
        }
        (contour, start_number)
    }

    /// Builds a loop from the contour and appends it to the result list.
    fn accept_contour(&mut self, contour: &[isize], start_number: usize) {
        let mut new_loop: Vec<Link> = Vec::new();
        for i in start_number..=contour.len() {
            let index_s = contour[i - 1];
            let index = index_s.unsigned_abs() as usize;
            let mut oriented = self.map_links[index - 1];
            if index_s < 0 { oriented.reverse(); }
            new_loop.push(oriented);
            self.start_indices.remove(&index_s);
        }
        self.loops.push(new_loop);
    }

    /// Returns the first node of the link, taking orientation (the sign of the
    /// index) into account.
    fn get_first_node(&self, index_s: isize) -> usize {
        let link = &self.map_links[index_s.unsigned_abs() as usize - 1];
        if index_s > 0 { link.node1 } else { link.node2 }
    }

    /// Returns the last node of the link, taking orientation into account.
    fn get_last_node(&self, index_s: isize) -> usize {
        let link = &self.map_links[index_s.unsigned_abs() as usize - 1];
        if index_s > 0 { link.node2 } else { link.node1 }
    }

    /// Marks hanging links starting from the given node, removing them from
    /// the start indices.
    fn mark_hang_chain(&mut self, start_node: usize, start_index_s: isize) {
        let mut node1 = start_node;
        let mut index_s = start_index_s;
        let mut index = index_s.unsigned_abs() as usize;
        let is_out = node1 == self.get_first_node(index_s);
        loop {
            // Check whether the current link is hanging: count other ways.
            let links = self.helper.get_adjacent_links(node1).to_vec();
            let mut n_edges = 0;
            for link in &links {
                let Some(ind) = self.find_link_index(link) else { continue };
                if ind == index { continue; }
                let mut a_ind = ind as isize;
                if (is_out && node1 == link.node1) || (!is_out && node1 == link.node2) {
                    a_ind = -a_ind;
                }
                if self.can_link_be_taken(a_ind) { n_edges += 1; }
            }
            if n_edges > 0 {
                // Leave this chain.
                break;
            }

            // Mark the current link as hanging.
            self.start_indices.remove(&index_s);
            self.hang_indices.insert(index_s);

            // Get the other node of the link.
            if is_out { node1 = self.get_last_node(index_s); }
            else { node1 = self.get_first_node(index_s); }

            // Find the next link in the chain.
            let next_links = self.helper.get_adjacent_links(node1).to_vec();
            let mut next_index_s: isize = 0;
            for link in &next_links {
                let Some(ind) = self.find_link_index(link) else { continue };
                if ind == index { continue; }
                let mut a_ind = ind as isize;
                if (is_out && node1 == link.node2) || (!is_out && node1 == link.node1) {
                    a_ind = -a_ind;
                }
                if self.can_link_be_taken(a_ind) {
                    if next_index_s == 0 {
                        next_index_s = a_ind;
                    } else {
                        // More than one way: stop the chain.
                        next_index_s = 0;
                        break;
                    }
                }
            }
            if next_index_s == 0 { break; }
            index_s = next_index_s;
            index = index_s.unsigned_abs() as usize;
        }
    }

    /// Returns whether the link appointed by the signed index can participate
    /// in a loop in the given orientation. A boundary edge can be taken once.
    fn can_link_be_taken(&self, index_s: isize) -> bool {
        self.start_indices.contains(&index_s)
    }
}

/// Helper that provides adjacency from a fixed table (used in tests and simple
/// callers). Node ids are 1-based.
pub struct FixedLinkHelper {
    adjacent: HashMap<usize, Vec<Link>>,
}

impl FixedLinkHelper {
    pub fn new() -> Self { Self { adjacent: HashMap::new() } }
    pub fn add_link(&mut self, link: Link) {
        if link.node1 == 0 || link.node2 == 0 { return; }
        self.adjacent.entry(link.node1).or_default().push(link);
        self.adjacent.entry(link.node2).or_default().push(link);
    }
}

impl Default for FixedLinkHelper {
    fn default() -> Self { Self::new() }
}

impl MakeLoopsHelper for FixedLinkHelper {
    fn get_adjacent_links(&self, node: usize) -> &[Link] {
        self.adjacent.get(&node).map(|v| v.as_slice()).unwrap_or(&[])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn six_disjoint_quads_form_six_loops() {
        // Six disjoint quadrilaterals; each quad uses 4 unique nodes (1-based).
        let mut helper = FixedLinkHelper::new();
        let mut links: Vec<Link> = Vec::new();
        for q in 0..6 {
            let base = q * 4;
            let a = base + 1;
            let b = base + 2;
            let c = base + 3;
            let d = base + 4;
            for (x, y) in [(a, b), (b, c), (c, d), (d, a)] {
                let l = Link::new(x, y);
                helper.add_link(l);
                links.push(l);
            }
        }
        let mut ml = PolyMakeLoops::new(&helper);
        for l in &links {
            ml.add_link(l);
        }
        let rc = ml.perform().unwrap();
        assert_ne!(rc & RC_LOOPS_DONE, 0);
        assert_eq!(ml.get_nb_loops(), 6);
        for i in 0..6 {
            let l = ml.get_loop(i);
            assert_eq!(l.len(), 4, "loop {i} has 4 links");
            // Each loop is closed: last link's last node == first link's first node.
            let first_first = if l[0].is_reversed() { l[0].node2 } else { l[0].node1 };
            let last_last = if l[l.len() - 1].is_reversed() { l[l.len() - 1].node1 } else { l[l.len() - 1].node2 };
            assert_eq!(last_last, first_first);
        }
    }

    #[test]
    fn square_single_loop() {
        let mut helper = FixedLinkHelper::new();
        for (x, y) in [(1, 2), (2, 3), (3, 4), (4, 1)] {
            helper.add_link(Link::new(x, y));
        }
        let mut ml = PolyMakeLoops::new(&helper);
        for (x, y) in [(1, 2), (2, 3), (3, 4), (4, 1)] {
            ml.add_link(&Link::new(x, y));
        }
        let rc = ml.perform().unwrap();
        assert_ne!(rc & RC_LOOPS_DONE, 0);
        assert_eq!(ml.get_nb_loops(), 1);
        assert_eq!(ml.get_loop(0).len(), 4);
        assert_eq!(rc & RC_HANGING_LINKS, 0);
    }

    #[test]
    fn dangling_chain_is_hanging() {
        let mut helper = FixedLinkHelper::new();
        // Triangle (1,2,3) plus a dangling edge (3,5).
        for (x, y) in [(1, 2), (2, 3), (3, 1), (3, 5)] {
            helper.add_link(Link::new(x, y));
        }
        let mut ml = PolyMakeLoops::new(&helper);
        for (x, y) in [(1, 2), (2, 3), (3, 1), (3, 5)] {
            ml.add_link(&Link::new(x, y));
        }
        let rc = ml.perform().unwrap();
        assert_ne!(rc & RC_LOOPS_DONE, 0);
        assert_ne!(rc & RC_HANGING_LINKS, 0);
        assert_eq!(ml.get_nb_loops(), 1);
        assert_eq!(ml.get_nb_hanging(), 1);
        assert_eq!(ml.get_hanging_links().len(), 1);
    }

    #[test]
    fn link_equality_and_reversal() {
        let a = Link::new(1, 2);
        let b = Link::new(2, 1);
        assert!(a.same_pair(&b));
        let mut c = a;
        c.reverse();
        assert!(c.is_reversed());
        assert_eq!(c.flags, LF_FWD.0 | LF_REVERSED.0);
        assert_eq!(c.flags & LF_BOTH.0, LF_FWD.0); // orientation part unchanged
    }

    #[test]
    fn set_orientation_excludes_link() {
        let mut helper = FixedLinkHelper::new();
        for (x, y) in [(1, 2), (2, 3), (3, 1)] {
            helper.add_link(Link::new(x, y));
        }
        let mut ml = PolyMakeLoops::new(&helper);
        ml.add_link(&Link::new(1, 2));
        ml.add_link(&Link::new(2, 3));
        ml.add_link(&Link::new(3, 1));
        // Exclude one link: the loop becomes a hanging open chain.
        ml.set_link_orientation(&Link::new(1, 2), LF_NONE);
        let rc = ml.perform().unwrap();
        assert_ne!(rc & RC_HANGING_LINKS, 0);
    }
}
