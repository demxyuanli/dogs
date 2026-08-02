//! Port of OCCT BRepMesh value types — Wave 2 Delaunay core.
//!
//! Sources (`src/ModelingAlgorithms/TKMesh/BRepMesh/`):
//! - `BRepMesh_Vertex.hxx` + `BRepMesh_DegreeOfFreedom.hxx`
//! - `BRepMesh_Triangle.hxx`
//! - `BRepMesh_Circle.hxx`
//! - `BRepMesh_Edge.hxx` / `BRepMesh_OrientedEdge.hxx`
//! - `BRepMesh_PairOfIndex.hxx`
//!
//! These are lightweight `Copy` value structs exchanged between the Delaunay
//! core (`Delaun`, `DataStructureOfDelaun`) and the mesh tools. The type and
//! field names are contractual — Wave 2 peers and Wave 3 consume them by the
//! exact names below.

use occt_core::gp::{GpPnt, GpPnt2d};
use occt_core::precision::PCONFUSION;

/// State/movability of a mesh vertex. Source: `BRepMesh_DegreeOfFreedom`.
///
/// Ordered exactly as the OCCT C enum (`Free = 0`, ..., `Deleted = 6`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum VertexState {
    Free,
    InVolume,
    OnSurface,
    OnCurve,
    Fixed,
    Frontier,
    Deleted,
}

impl VertexState {
    /// Numeric value as in the OCCT C enum.
    pub fn index(self) -> usize {
        self as usize
    }

    /// Human-readable name.
    pub fn to_str(self) -> &'static str {
        match self {
            Self::Free => "Free",
            Self::InVolume => "InVolume",
            Self::OnSurface => "OnSurface",
            Self::OnCurve => "OnCurve",
            Self::Fixed => "Fixed",
            Self::Frontier => "Frontier",
            Self::Deleted => "Deleted",
        }
    }
}

/// Vertex of the mesh in parametric space. Source: `BRepMesh_Vertex`.
///
/// OCCT keeps the 3d association as an external map index (`Location3d`); the
/// contract inlines the 3d point together with that index for the Rust port.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DelaunVertex {
    /// Position in parametric space (UV).
    pub location: GpPnt2d,
    /// Associated 3d point.
    pub p3d: GpPnt,
    /// Index of the 3d point / external map key.
    pub index: i32,
    /// Movability of the vertex.
    pub state: VertexState,
}

impl DelaunVertex {
    /// Full constructor (parametric position, 3d point, 3d index, state).
    pub const fn new(location: GpPnt2d, p3d: GpPnt, index: i32, state: VertexState) -> Self {
        Self { location, p3d, index, state }
    }

    /// Creates a vertex from explicit UV coordinates, no 3d association.
    /// Source: `BRepMesh_Vertex(const double, const double, ...)`.
    pub const fn new_parametric(u: f64, v: f64, state: VertexState) -> Self {
        Self { location: GpPnt2d::new(u, v), p3d: GpPnt::zero(), index: 0, state }
    }

    /// Returns the parametric position.
    pub const fn location(&self) -> GpPnt2d {
        self.location
    }

    /// Returns the associated 3d point.
    pub const fn p3d(&self) -> GpPnt {
        self.p3d
    }

    /// Returns the index of the associated 3d point.
    pub const fn index(&self) -> i32 {
        self.index
    }

    /// Returns the movability of the vertex.
    pub const fn state(&self) -> VertexState {
        self.state
    }

    /// Alias for [`Self::state`] (OCCT `Movability`).
    pub const fn movability(&self) -> VertexState {
        self.state
    }

    /// Sets the parametric position.
    pub fn set_location(&mut self, location: GpPnt2d) {
        self.location = location;
    }

    /// Sets the associated 3d point.
    pub fn set_p3d(&mut self, p3d: GpPnt) {
        self.p3d = p3d;
    }

    /// Sets the index of the associated 3d point.
    pub fn set_index(&mut self, index: i32) {
        self.index = index;
    }

    /// Sets the movability of the vertex. Source: `SetMovability`.
    pub fn set_state(&mut self, state: VertexState) {
        self.state = state;
    }

    /// Alias for [`Self::set_state`] (OCCT `SetMovability`).
    pub fn set_movability(&mut self, state: VertexState) {
        self.state = state;
    }

    /// Checks for equality with another vertex. Source: `BRepMesh_Vertex::IsEqual`.
    ///
    /// Two vertices are equal when neither is `Deleted` and their parametric
    /// locations coincide within `Precision::PConfusion`.
    pub fn is_equal(&self, other: &Self) -> bool {
        if self.state == VertexState::Deleted || other.state == VertexState::Deleted {
            return false;
        }
        self.location.distance(&other.location) <= PCONFUSION
    }

    /// Alias for [`Self::is_equal`] (OCCT `operator==`).
    pub fn equals(&self, other: &Self) -> bool {
        self.is_equal(other)
    }
}

impl Default for DelaunVertex {
    fn default() -> Self {
        Self::new(GpPnt2d::zero(), GpPnt::zero(), 0, VertexState::Free)
    }
}

/// Triangle of the mesh consisting of oriented links. Source: `BRepMesh_Triangle`.
///
/// OCCT stores three link indices plus a per-link orientation bit; the contract
/// flattens orientations into the three corner vertex indices (link `k` joins
/// vertex `k` and vertex `(k+1) % 3`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DelaunTriangle {
    /// Indices of the three links.
    pub link_indices: [i32; 3],
    /// Indices of the three corner vertices.
    pub vertex_indices: [i32; 3],
}

impl DelaunTriangle {
    /// Constructor from link and vertex index arrays.
    pub const fn new(link_indices: [i32; 3], vertex_indices: [i32; 3]) -> Self {
        Self { link_indices, vertex_indices }
    }

    /// Returns the link indices.
    pub const fn link_indices(&self) -> &[i32; 3] {
        &self.link_indices
    }

    /// Returns the vertex indices.
    pub const fn vertex_indices(&self) -> &[i32; 3] {
        &self.vertex_indices
    }

    /// Returns the link index at position `k` (`0..=2`).
    pub const fn link_at(&self, k: usize) -> i32 {
        self.link_indices[k]
    }

    /// Returns the vertex index at position `k` (`0..=2`).
    pub const fn vertex_at(&self, k: usize) -> i32 {
        self.vertex_indices[k]
    }

    /// Checks for equality with another triangle. Source: `BRepMesh_Triangle::IsEqual`.
    ///
    /// A triangle is equal regardless of which vertex is "first", so the link
    /// indices are compared under cyclic rotation (as OCCT does for its edges).
    pub fn is_equal(&self, other: &Self) -> bool {
        for r in 0..3 {
            if (0..3).all(|k| self.link_indices[(k + r) % 3] == other.link_indices[k]) {
                return true;
            }
        }
        false
    }

    /// Alias for [`Self::is_equal`] (OCCT `operator==`).
    pub fn equals(&self, other: &Self) -> bool {
        self.is_equal(other)
    }
}

impl Default for DelaunTriangle {
    fn default() -> Self {
        Self { link_indices: [0; 3], vertex_indices: [0; 3] }
    }
}

/// Describes a 2d circle using only a center and a squared radius.
/// Source: `BRepMesh_Circle`.
///
/// OCCT stores the plain radius (a negative value marks an invalid circle); the
/// contract stores the squared radius plus an explicit `is_created` flag.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DelaunCircle {
    /// Center of the circle.
    pub center: GpPnt2d,
    /// Squared radius.
    pub radius_sq: f64,
    /// Whether the circle has been created/validated.
    pub is_created: bool,
}

impl DelaunCircle {
    /// Constructor from a center and a squared radius.
    pub const fn new(center: GpPnt2d, radius_sq: f64) -> Self {
        Self { center, radius_sq, is_created: true }
    }

    /// Constructor from a center and a plain radius (stored squared).
    pub const fn with_radius(center: GpPnt2d, radius: f64) -> Self {
        Self { center, radius_sq: radius * radius, is_created: radius >= 0.0 }
    }

    /// Returns the center of the circle.
    pub const fn center(&self) -> GpPnt2d {
        self.center
    }

    /// Returns the squared radius.
    pub const fn radius_sq(&self) -> f64 {
        self.radius_sq
    }

    /// Returns the plain radius (sqrt of the squared radius).
    pub fn radius(&self) -> f64 {
        self.radius_sq.sqrt()
    }

    /// Returns whether the circle is created.
    pub const fn is_created(&self) -> bool {
        self.is_created
    }

    /// Sets the center. Source: `SetLocation`.
    pub fn set_center(&mut self, center: GpPnt2d) {
        self.center = center;
    }

    /// Sets the squared radius.
    pub fn set_radius_sq(&mut self, radius_sq: f64) {
        self.radius_sq = radius_sq;
    }

    /// Sets the plain radius (stored squared). Source: `SetRadius`.
    pub fn set_radius(&mut self, radius: f64) {
        self.radius_sq = radius * radius;
        self.is_created = radius >= 0.0;
    }

    /// Marks the circle as created or not.
    pub fn set_is_created(&mut self, is_created: bool) {
        self.is_created = is_created;
    }
}

impl Default for DelaunCircle {
    fn default() -> Self {
        Self { center: GpPnt2d::zero(), radius_sq: 0.0, is_created: false }
    }
}

/// Link (edge) of the mesh connecting two vertices. Source: `BRepMesh_Edge`
/// + `BRepMesh_OrientedEdge`.
///
/// The contract additionally carries the parameter range and the link index
/// (OCCT stores those in the surrounding data structure).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DelaunLink {
    /// Index of the first vertex.
    pub v1: i32,
    /// Index of the second vertex.
    pub v2: i32,
    /// Parameter of the first vertex on the link's curve.
    pub first_param: f64,
    /// Parameter of the second vertex on the link's curve.
    pub last_param: f64,
    /// Index of the link in the mesh.
    pub index: i32,
}

impl DelaunLink {
    /// Full constructor.
    pub const fn new(v1: i32, v2: i32, first_param: f64, last_param: f64, index: i32) -> Self {
        Self { v1, v2, first_param, last_param, index }
    }

    /// Constructor without parameter range (parameters default to `0.0`).
    pub const fn new_unparameterized(v1: i32, v2: i32, index: i32) -> Self {
        Self { v1, v2, first_param: 0.0, last_param: 0.0, index }
    }

    /// Returns the index of the first vertex.
    pub const fn v1(&self) -> i32 {
        self.v1
    }

    /// Returns the index of the second vertex.
    pub const fn v2(&self) -> i32 {
        self.v2
    }

    /// Alias for [`Self::v1`] (OCCT `FirstNode`).
    pub const fn first_node(&self) -> i32 {
        self.v1
    }

    /// Alias for [`Self::v2`] (OCCT `LastNode`).
    pub const fn last_node(&self) -> i32 {
        self.v2
    }

    /// Returns the parameter of the first vertex.
    pub const fn first_param(&self) -> f64 {
        self.first_param
    }

    /// Returns the parameter of the second vertex.
    pub const fn last_param(&self) -> f64 {
        self.last_param
    }

    /// Returns the index of the link.
    pub const fn index(&self) -> i32 {
        self.index
    }

    /// Sets the index of the link.
    pub fn set_index(&mut self, index: i32) {
        self.index = index;
    }

    /// Returns a copy of this link with the orientation reversed (vertices and
    /// their parameters are swapped).
    pub const fn reversed(&self) -> Self {
        Self { v1: self.v2, v2: self.v1, first_param: self.last_param, last_param: self.first_param, index: self.index }
    }

    /// Checks if this link and `other` have the same orientation.
    /// Source: `BRepMesh_OrientedEdge::IsEqual`.
    pub const fn is_same_orientation(&self, other: &Self) -> bool {
        self.v1 == other.v1 && self.v2 == other.v2
    }

    /// Checks for equality with another link (either orientation).
    /// Source: `BRepMesh_Edge::IsEqual`.
    pub const fn is_equal(&self, other: &Self) -> bool {
        self.is_same_orientation(other)
            || (self.v1 == other.v2 && self.v2 == other.v1)
    }

    /// Alias for [`Self::is_equal`] (OCCT `operator==`).
    pub const fn equals(&self, other: &Self) -> bool {
        self.is_equal(other)
    }
}

impl Default for DelaunLink {
    fn default() -> Self {
        Self::new(0, 0, 0.0, 0.0, 0)
    }
}

/// Reference to a link together with its orientation. Source:
/// `BRepMesh_OrientedEdge` (adapted to an index + direction flag).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DelaunOrientedEdge {
    /// Index of the link in the mesh.
    pub link_index: i32,
    /// Whether the link is traversed in its forward direction.
    pub is_forward: bool,
}

impl DelaunOrientedEdge {
    /// Constructor.
    pub const fn new(link_index: i32, is_forward: bool) -> Self {
        Self { link_index, is_forward }
    }

    /// Returns the index of the link.
    pub const fn link_index(&self) -> i32 {
        self.link_index
    }

    /// Returns whether the link is traversed forward.
    pub const fn is_forward(&self) -> bool {
        self.is_forward
    }

    /// Returns a copy with the opposite orientation.
    pub const fn reversed(&self) -> Self {
        Self { link_index: self.link_index, is_forward: !self.is_forward }
    }

    /// Checks for equality (same link, same orientation).
    pub const fn is_equal(&self, other: &Self) -> bool {
        self.link_index == other.link_index && self.is_forward == other.is_forward
    }

    /// Alias for [`Self::is_equal`].
    pub const fn equals(&self, other: &Self) -> bool {
        self.is_equal(other)
    }
}

impl Default for DelaunOrientedEdge {
    fn default() -> Self {
        Self { link_index: -1, is_forward: false }
    }
}

/// A pair of integer indices, used to store element indices connected to a
/// link. Source: `BRepMesh_PairOfIndex`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DelaunPairOfIndex {
    /// First index (`-1` when unset).
    pub first: i32,
    /// Second index (`-1` when unset).
    pub second: i32,
}

impl DelaunPairOfIndex {
    /// Constructor.
    pub const fn new(first: i32, second: i32) -> Self {
        Self { first, second }
    }

    /// Clears both indices. Source: `Clear`.
    pub fn clear(&mut self) {
        self.first = -1;
        self.second = -1;
    }

    /// Appends an index to the pair. Source: `Append`.
    ///
    /// Fills the first free slot; panics when the pair already holds two
    /// indices (OCCT throws `Standard_OutOfRange`).
    pub fn append(&mut self, the_index: i32) {
        if self.first < 0 {
            self.first = the_index;
        } else if self.second >= 0 {
            panic!("DelaunPairOfIndex::append, more than two index to store");
        } else {
            self.second = the_index;
        }
    }

    /// Prepends an index to the pair. Source: `Prepend`.
    ///
    /// Shifts the first index into the second slot; panics when the pair is
    /// already full.
    pub fn prepend(&mut self, the_index: i32) {
        if self.second >= 0 {
            panic!("DelaunPairOfIndex::prepend, more than two index to store");
        }
        self.second = self.first;
        self.first = the_index;
    }

    /// Returns whether the pair is empty. Source: `IsEmpty`.
    pub const fn is_empty(&self) -> bool {
        self.first < 0
    }

    /// Returns the number of initialized indices (`0`, `1` or `2`).
    /// Source: `Extent`.
    pub const fn extent(&self) -> usize {
        if self.first < 0 { 0 } else if self.second < 0 { 1 } else { 2 }
    }

    /// Returns the first index of the pair. Source: `FirstIndex`.
    pub const fn first_index(&self) -> i32 {
        self.first
    }

    /// Returns the last (or only) index of the pair. Source: `LastIndex`.
    pub const fn last_index(&self) -> i32 {
        if self.second < 0 { self.first } else { self.second }
    }

    /// Returns the index at the 1-based position `pos` (`1` or `2`).
    /// Source: `Index`. Panics for any other position, as OCCT throws.
    pub const fn index(&self, pos: usize) -> i32 {
        match pos {
            1 => self.first,
            2 => self.second,
            _ => panic!("DelaunPairOfIndex::index, requested index is out of range"),
        }
    }

    /// Sets the index at the 1-based position `pos` (`1` or `2`).
    /// Source: `SetIndex`.
    pub fn set_index(&mut self, pos: usize, the_index: i32) {
        match pos {
            1 => self.first = the_index,
            2 => self.second = the_index,
            _ => panic!("DelaunPairOfIndex::set_index, requested index is out of range"),
        }
    }

    /// Removes the index at the 1-based position `pos` (`1` or `2`).
    /// Source: `RemoveIndex`. Removing the first slot shifts the second into it.
    pub fn remove_index(&mut self, pos: usize) {
        match pos {
            1 => self.first = self.second,
            2 => {}
            _ => panic!("DelaunPairOfIndex::remove_index, requested index is out of range"),
        }
        self.second = -1;
    }
}

impl Default for DelaunPairOfIndex {
    fn default() -> Self {
        Self { first: -1, second: -1 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vertex_state_enum_matches_occt_order() {
        let all = [
            VertexState::Free,
            VertexState::InVolume,
            VertexState::OnSurface,
            VertexState::OnCurve,
            VertexState::Fixed,
            VertexState::Frontier,
            VertexState::Deleted,
        ];
        assert_eq!(all.len(), 7);
        for (i, v) in all.iter().enumerate() {
            assert_eq!(v.index(), i);
        }
        assert_eq!(VertexState::Free.index(), 0);
        assert_eq!(VertexState::Deleted.index(), 6);
        assert_eq!(VertexState::OnSurface.to_str(), "OnSurface");
    }

    #[test]
    fn vertex_construction_getters_and_equality() {
        let v = DelaunVertex::new(
            GpPnt2d::new(1.0, 2.0),
            GpPnt::new(3.0, 4.0, 5.0),
            7,
            VertexState::Fixed,
        );
        assert_eq!(v.location(), GpPnt2d::new(1.0, 2.0));
        assert_eq!(v.p3d(), GpPnt::new(3.0, 4.0, 5.0));
        assert_eq!(v.index(), 7);
        assert_eq!(v.state(), VertexState::Fixed);
        assert_eq!(v.movability(), VertexState::Fixed);

        // Same UV within PConfusion => equal.
        let near = DelaunVertex::new(
            GpPnt2d::new(1.0 + 1e-10, 2.0),
            GpPnt::new(9.0, 9.0, 9.0),
            0,
            VertexState::Free,
        );
        assert!(v.is_equal(&near), "vertices within PConfusion must be equal");

        // Different UV => not equal.
        let far = DelaunVertex::new_parametric(10.0, 20.0, VertexState::Free);
        assert!(!v.is_equal(&far));

        // A Deleted vertex never equals anything.
        let mut del = v;
        del.set_state(VertexState::Deleted);
        assert!(!v.is_equal(&del));
        assert!(!del.is_equal(&near));

        // Defaults.
        let d = DelaunVertex::default();
        assert_eq!(d.location(), GpPnt2d::zero());
        assert_eq!(d.state(), VertexState::Free);
        assert_eq!(d.index(), 0);
    }

    #[test]
    fn triangle_cyclic_equality() {
        let t = DelaunTriangle::new([1, 2, 3], [10, 11, 12]);
        // Cyclic rotation of the same triangle is equal.
        let rotated = DelaunTriangle::new([3, 1, 2], [12, 10, 11]);
        assert!(t.is_equal(&rotated));
        assert_eq!(t, t);
        // Reversed (mirror) order is NOT a cyclic rotation.
        let mirrored = DelaunTriangle::new([1, 3, 2], [10, 12, 11]);
        assert!(!t.is_equal(&mirrored));
        // Different links are not equal.
        let other = DelaunTriangle::new([5, 6, 7], [20, 21, 22]);
        assert!(!t.is_equal(&other));
        assert_eq!(t.link_at(1), 2);
        assert_eq!(t.vertex_at(2), 12);
        assert_eq!(DelaunTriangle::default().link_indices(), &[0, 0, 0]);
    }

    #[test]
    fn circle_stores_squared_radius() {
        let c = DelaunCircle::with_radius(GpPnt2d::new(1.0, 2.0), 3.0);
        assert_eq!(c.center(), GpPnt2d::new(1.0, 2.0));
        assert!((c.radius_sq() - 9.0).abs() < 1e-12);
        assert!((c.radius() - 3.0).abs() < 1e-12);
        assert!(c.is_created());

        let mut c = c;
        c.set_radius(-1.0); // OCCT invalidates with a negative radius.
        assert!(!c.is_created());
        assert!((c.radius_sq() - 1.0).abs() < 1e-12);

        let d = DelaunCircle::default();
        assert!(!d.is_created());
        assert_eq!(d.radius_sq(), 0.0);

        let s = DelaunCircle::new(GpPnt2d::zero(), 4.0);
        assert_eq!(s.radius(), 2.0);
    }

    #[test]
    fn link_orientation_and_reversal() {
        let l = DelaunLink::new_unparameterized(1, 2, 5);
        assert_eq!(l.v1(), 1);
        assert_eq!(l.v2(), 2);
        assert_eq!(l.first_node(), 1);
        assert_eq!(l.last_node(), 2);
        assert_eq!(l.index(), 5);
        assert_eq!(l.first_param(), 0.0);

        let r = l.reversed();
        assert_eq!(r.v1(), 2);
        assert_eq!(r.v2(), 1);
        // Undirected equality holds for the reversed link.
        assert!(l.is_equal(&r));
        assert!(!l.is_same_orientation(&r));
        assert!(l.is_same_orientation(&l));
        // Different link is not equal.
        let other = DelaunLink::new(1, 3, 0.0, 1.0, 6);
        assert!(!l.is_equal(&other));

        // Parameters swap with the vertices.
        let p = DelaunLink::new(1, 2, 0.25, 0.75, 9);
        let pr = p.reversed();
        assert!((pr.first_param() - 0.75).abs() < 1e-12);
        assert!((pr.last_param() - 0.25).abs() < 1e-12);
    }

    #[test]
    fn oriented_edge_forward_reverse() {
        let e = DelaunOrientedEdge::new(3, true);
        assert_eq!(e.link_index(), 3);
        assert!(e.is_forward());
        assert!(e.is_equal(&DelaunOrientedEdge::new(3, true)));
        assert!(!e.is_equal(&DelaunOrientedEdge::new(3, false)));

        let r = e.reversed();
        assert_eq!(r.link_index(), 3);
        assert!(!r.is_forward());

        let d = DelaunOrientedEdge::default();
        assert_eq!(d.link_index(), -1);
        assert!(!d.is_forward());
    }

    #[test]
    fn pair_of_index_composition() {
        let mut p = DelaunPairOfIndex::default();
        assert!(p.is_empty());
        assert_eq!(p.extent(), 0);

        p.append(10);
        assert!(!p.is_empty());
        assert_eq!(p.extent(), 1);
        assert_eq!(p.first_index(), 10);
        assert_eq!(p.last_index(), 10);

        // Prepend shifts the first index into the second slot.
        p.prepend(30);
        assert_eq!(p.extent(), 2);
        assert_eq!(p.first, 30);
        assert_eq!(p.second, 10);

        // Append to a single-element pair fills the second slot.
        let mut q = DelaunPairOfIndex::new(5, -1);
        q.append(7);
        assert_eq!(q.extent(), 2);
        assert_eq!(q.index(1), 5);
        assert_eq!(q.index(2), 7);
        assert_eq!(q.last_index(), 7);

        // SetIndex targets 1-based positions.
        p.set_index(2, 40);
        assert_eq!(p.index(2), 40);

        // RemoveIndex(1) promotes the second index into the first slot.
        p.remove_index(1);
        assert_eq!(p.first, 40);
        assert_eq!(p.second, -1);
        assert_eq!(p.extent(), 1);

        // Clear resets to empty.
        p.clear();
        assert!(p.is_empty());
        assert_eq!(p.extent(), 0);
    }

    #[test]
    #[should_panic(expected = "more than two index")]
    fn pair_of_index_append_overflows() {
        let mut p = DelaunPairOfIndex::new(1, 2);
        p.append(3);
    }
}
