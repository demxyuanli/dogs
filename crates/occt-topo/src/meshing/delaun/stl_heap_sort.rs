//! Faithful port of the C++ standard library heap algorithms that OCCT calls
//! in `BRepMesh_Delaun` before feeding vertices to `compute`:
//!
//! ```text
//! ComparatorOfIndexedVertexOfDelaun aCmp(myMeshData);
//! std::make_heap(v.begin(), v.end(), aCmp);
//! std::sort_heap(v.begin(), v.end(), aCmp);
//! ```
//!
//! Sources:
//! - `BRepMesh_Delaun.cxx:332-334` (`perform`) and `:2415-2417` (`AddVertices`).
//! - `BRepMesh_Delaun.cxx:49-75` — `ComparatorOfVertexOfDelaun` compares
//!   `Coord().X() + Coord().Y()` with a plain `<`.
//! - MSVC STL `__msvc_heap_algorithms.hpp:52-131` — `_Pop_heap_hole_by_index`,
//!   `_Push_heap_by_index`, `_Make_heap_unchecked`; and `algorithm:7147-7152`
//!   `_Sort_heap_unchecked` + `:98-106` `_Pop_heap_unchecked`.
//!
//! The algorithms are reproduced op-for-op (not `sort`): equal keys end up in
//! the STL heap order, which decides the UV diagonal of cocircular triangles.
//! Do not replace these with a generic sort.

/// `_Push_heap_by_index` (`__msvc_heap_algorithms.hpp:52-63`).
fn push_heap_by_index<T: Copy, F: Fn(&T, &T) -> bool>(
    v: &mut [T],
    hole0: isize,
    top: isize,
    val: T,
    less: &F,
) {
    let mut hole = hole0;
    while top < hole {
        let idx = (hole - 1) >> 1;
        if !less(&v[idx as usize], &val) {
            break;
        }
        v[hole as usize] = v[idx as usize];
        hole = idx;
    }
    v[hole as usize] = val;
}

/// `_Pop_heap_hole_by_index` (`__msvc_heap_algorithms.hpp:57-84`).
fn pop_heap_hole_by_index<T: Copy, F: Fn(&T, &T) -> bool>(
    v: &mut [T],
    hole0: isize,
    bottom: isize,
    val: T,
    less: &F,
) {
    let top = hole0;
    let mut hole = hole0;
    let mut idx = hole0;
    let max_sequence_non_leaf = (bottom - 1) >> 1;
    while idx < max_sequence_non_leaf {
        idx = 2 * idx + 2;
        if less(&v[idx as usize], &v[(idx - 1) as usize]) {
            idx -= 1;
        }
        v[hole as usize] = v[idx as usize];
        hole = idx;
    }
    if idx == max_sequence_non_leaf && bottom % 2 == 0 {
        v[hole as usize] = v[(bottom - 1) as usize];
        hole = bottom - 1;
    }
    push_heap_by_index(v, hole, top, val, less);
}

/// `_Make_heap_unchecked` (`__msvc_heap_algorithms.hpp:121-131`).
pub fn make_heap<T: Copy, F: Fn(&T, &T) -> bool>(v: &mut [T], less: &F) {
    let bottom = v.len() as isize;
    let mut hole = bottom >> 1;
    while hole > 0 {
        hole -= 1;
        let val = v[hole as usize];
        pop_heap_hole_by_index(v, hole, bottom, val, less);
    }
}

/// `_Sort_heap_unchecked` (`algorithm:7147-7152`) driving `_Pop_heap_unchecked`
/// (`__msvc_heap_algorithms.hpp:98-106`) — pops the max to `last - 1` and
/// re-heaps the remaining hole with the popped-back value.
pub fn sort_heap<T: Copy, F: Fn(&T, &T) -> bool>(v: &mut [T], less: &F) {
    let mut last = v.len() as isize;
    while last >= 2 {
        last -= 1;
        let val = v[last as usize];
        v[last as usize] = v[0];
        pop_heap_hole_by_index(v, 0, last, val, less);
    }
}

/// `make_heap` + `sort_heap` with the `X + Y` key of
/// `ComparatorOfIndexedVertexOfDelaun` (`BRepMesh_Delaun.cxx:58-75`).
pub fn heap_sort_by_xy_sum<T: Copy, K: Fn(&T) -> (f64, f64)>(v: &mut [T], key: K) {
    let less = |a: &T, b: &T| {
        let (ax, ay) = key(a);
        let (bx, by) = key(b);
        ax + ay < bx + by
    };
    make_heap(v, &less);
    sort_heap(v, &less);
}
