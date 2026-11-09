//! 3D `BVH_RadixSorter` and the `BVH_Set` services it needs.
//!
//! Source: `FoundationClasses/TKMath/BVH/BVH_RadixSorter.hxx`:
//! * `BVH::EncodeMortonCode` (`hxx:64-83`), the 10-bit interleave the port
//!   computes bit-wise instead of through the 256-entry `THE_MORTON_LUT`
//!   (`hxx:33-62`) - the produced code is identical.
//! * `BVH::RadixSorter::Sort` / `perform` (`hxx:197-227`): MSD radix sort on
//!   the 30-bit Morton code, `std::partition` per digit.
//! * `BVH_RadixSorter::Perform` (`hxx:232-297`): Morton code assignment,
//!   sort, and the in-place rearrangement that calls `BVH_Set::Swap`.
//! * `BVH_Set::Box()` (`BVH_Set.hxx:31-44`) is [`BvhSet3::box_all`].
#![allow(dead_code)]

use super::box_3d::{BvhBox3d, BvhVec3d};

/// The subset of `BVH_Set<double, 3>` the sorter and the builders call.
pub trait BvhSet3 {
    /// `Size()` (`BVH_Set.hxx:24`).
    fn size(&self) -> usize;

    /// `Box(theIndex)` (`BVH_Set.hxx:27`), 0-based.
    fn box_at(&self, index: usize) -> BvhBox3d;

    /// `Swap(theIndex1, theIndex2)` (`BVH_Set.hxx:30`).
    fn swap_elements(&mut self, index1: usize, index2: usize);

    /// `Center(theIndex, theAxis)` (`BVH_Set.hxx:36-39`).
    fn center_at(&self, index: usize, axis: usize) -> f64 {
        self.box_at(index).center_axis(axis)
    }

    /// `BVH_Set::Box()` (`BVH_Set.hxx:31-44`): union of all element boxes.
    fn box_all(&self) -> BvhBox3d {
        let mut a_box = BvhBox3d::new();
        for i in 0..self.size() {
            a_box.combine(&self.box_at(i));
        }
        a_box
    }
}

/// `BVH::THE_NODE_MIN_SIZE` (`BVH_Constants.hxx:39`).
pub const THE_NODE_MIN_SIZE: f64 = 1.0e-5;

/// Pair of Morton code and primitive id (`BVH_EncodedLink`,
/// `BVH_RadixSorter.hxx:28`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BvhEncodedLink {
    /// `first`
    pub code: u32,
    /// `second` - 0-based primitive index.
    pub prim: usize,
}

/// Spreads the 10 low bits of `value` to every third bit.
fn spread_10_bits(value: u32) -> u32 {
    let mut out = 0u32;
    for i in 0..10 {
        out |= ((value >> i) & 1) << (3 * i);
    }
    out
}

/// `BVH::EncodeMortonCode(x, y, z)` (`BVH_RadixSorter.hxx:64-83`).
pub fn encode_morton_code(voxel_x: u32, voxel_y: u32, voxel_z: u32) -> u32 {
    spread_10_bits(voxel_x & 0x3FF)
        | (spread_10_bits(voxel_y & 0x3FF) << 1)
        | (spread_10_bits(voxel_z & 0x3FF) << 2)
}

/// `std::partition` with `BitPredicate(theDigit)`: the links whose bit is 0
/// move to the front. Returns the offset of the first 1-bit element.
fn partition_by_bit(links: &mut [BvhEncodedLink], digit: u32) -> usize {
    let bit = 1u32 << digit;
    let mut zeros: Vec<BvhEncodedLink> = Vec::new();
    let mut ones: Vec<BvhEncodedLink> = Vec::new();
    for link in links.iter() {
        if link.code & bit == 0 {
            zeros.push(*link);
        } else {
            ones.push(*link);
        }
    }
    let offset = zeros.len();
    zeros.extend(ones);
    links.copy_from_slice(&zeros);
    offset
}

/// `BVH::RadixSorter::perform` (`BVH_RadixSorter.hxx:217-227`).
fn radix_perform(mut links: &mut [BvhEncodedLink], mut digit: i32) {
    while !links.is_empty() && digit >= 0 {
        let offset = partition_by_bit(links, digit as u32);
        digit -= 1;
        let (left, right) = links.split_at_mut(offset);
        radix_perform(left, digit);
        links = right;
    }
}

/// `BVH::RadixSorter::Sort` (`BVH_RadixSorter.hxx:197-215`).
pub fn radix_sort(links: &mut [BvhEncodedLink], digit: i32) {
    if digit < 24 {
        radix_perform(links, digit);
    } else {
        let offset = partition_by_bit(links, digit as u32);
        let (left, right) = links.split_at_mut(offset);
        radix_sort(left, digit - 1);
        radix_sort(right, digit - 1);
    }
}

/// `BVH_RadixSorter<double, 3>` (`BVH_RadixSorter.hxx:86-112`).
#[derive(Debug, Clone, Default)]
pub struct BvhRadixSorter3 {
    /// `myEncodedLinks` (`hxx:113`).
    links: Vec<BvhEncodedLink>,
}

impl BvhRadixSorter3 {
    /// `EncodedLinks()` (`hxx:105-108`).
    pub fn encoded_links(&self) -> &[BvhEncodedLink] {
        &self.links
    }

    /// `BVH_RadixSorter::Perform(theSet)` (`hxx:100-103`).
    pub fn perform<S: BvhSet3>(&mut self, set: &mut S) {
        let size = set.size();
        if size == 0 {
            return;
        }
        let scene_box = set.box_all();
        self.perform_range(set, 0, size - 1, &scene_box);
    }

    /// `BVH_RadixSorter::Perform(theSet, theStart, theFinal)`
    /// (`hxx:232-297`).
    pub fn perform_range<S: BvhSet3>(
        &mut self,
        set: &mut S,
        the_start: usize,
        the_final: usize,
        scene_box: &BvhBox3d,
    ) {
        let dimension: i32 = 1024;
        if !scene_box.is_valid() {
            return;
        }

        let scene_min = scene_box.corner_min();
        let scene_max = scene_box.corner_max();

        // `BVH::BoxMinMax::CwiseMax(aNodeMinSizeVecT, aSceneMax - aSceneMin)`
        // (`hxx:244-247`).
        let node_min_size = [
            THE_NODE_MIN_SIZE.max(scene_max.x - scene_min.x),
            THE_NODE_MIN_SIZE.max(scene_max.y - scene_min.y),
            THE_NODE_MIN_SIZE.max(scene_max.z - scene_min.z),
        ];
        let reverse_size = [
            f64::from(dimension) / node_min_size[0],
            f64::from(dimension) / node_min_size[1],
            f64::from(dimension) / node_min_size[2],
        ];

        let nb_links = the_final - the_start + 1;
        self.links = vec![BvhEncodedLink { code: 0, prim: 0 }; nb_links];

        // Step 1 -- Assign Morton code to each primitive (hxx:251-273)
        for prim in the_start..=the_final {
            let center = set.box_at(prim).center();
            let voxel_f = BvhVec3d::new(
                (center.x - scene_min.x) * reverse_size[0],
                (center.y - scene_min.y) * reverse_size[1],
                (center.z - scene_min.z) * reverse_size[2],
            );
            let voxel_x = clamp_voxel(voxel_f.x, dimension);
            let voxel_y = clamp_voxel(voxel_f.y, dimension);
            let voxel_z = clamp_voxel(voxel_f.z, dimension);

            let code = encode_morton_code(voxel_x as u32, voxel_y as u32, voxel_z as u32);
            self.links[prim - the_start] = BvhEncodedLink { code, prim };
        }

        // Step 2 -- Sort primitives by their Morton codes (hxx:275)
        radix_sort(&mut self.links, 29);

        // `aLinkMap(encodedLink.second) = aLinkIdx` (hxx:277-283). The OCCT
        // array is bounded by `theStart..theFinal`; the port stores the whole
        // primitive range so the map is indexed by the primitive id.
        let mut link_map = vec![0usize; the_final + 1];
        for link_idx in 0..nb_links {
            link_map[self.links[link_idx].prim] = link_idx + the_start;
        }

        // Step 3 -- Rearranging primitive list in place (hxx:286-297)
        let mut prim_idx = the_start;
        while prim_idx <= the_final {
            let sort_idx = link_map[prim_idx];
            if prim_idx != sort_idx {
                set.swap_elements(prim_idx, sort_idx);
                link_map.swap(prim_idx, sort_idx);
            } else {
                prim_idx += 1;
            }
        }
    }
}

/// `std::clamp(IntFloor(comp), 0, aDimension - 1)` (`hxx:257-263`).
fn clamp_voxel(value: f64, dimension: i32) -> i32 {
    let floored = value.floor() as i64;
    let clamped = floored.clamp(0, i64::from(dimension) - 1);
    clamped as i32
}
