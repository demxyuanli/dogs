//! Topology location — nested datum transforms. Source: `TopLoc/`
pub mod datum;
use std::sync::{Arc, OnceLock};
use crate::gp::{GpPnt, GpTrsf};

/// Nested chain of transformations (like a singly-linked list of gp_Trsf).
///
/// The chain is immutable, so it is shared instead of copied: a
/// [`TopLocLocation`] is one `Arc` to its outermost link. Copying a
/// `TopoDS_Shape` copies its location on every `oriented` / `located` /
/// `clone`, and that copy is now a reference count rather than a deep walk of
/// the chain.
///
/// Each link memoizes the two values OCCT recomputes on every query:
/// [`TopLocLocation::transformation`] (the product of the chain) and
/// [`TopLocLocation::key`] (the comparison unit of `IsEqual`).
///
/// Identity links are never stored. `TopLoc_Location::Multiplied` returns
/// `*this` / `Other` unchanged when either side is identity, so
/// [`TopLocLocation::composed`] with an identity transform is exactly the
/// parent; a chain therefore has one link per real datum, which is also what
/// [`TopLocLocation::elementary_transforms`] relied on when it skipped identity
/// links.
#[derive(Clone)]
pub struct TopLocLocation {
    head: Option<Arc<LocNode>>,
}

/// One link of the chain: the datum's transform plus the shared state derived
/// from the rest of the chain.
#[derive(Debug)]
struct LocNode {
    /// This datum's `Transformation()`.
    transform: GpTrsf,
    /// The next, older link (a `TopLoc_Location`'s `myNext`).
    parent: Option<Arc<LocNode>>,
    /// Links in this chain; `1` for a single datum.
    depth: u32,
    /// `Transformation()` of this chain, cached when its value is bit-for-bit
    /// the product OCCT's left-to-right fold produces. The fold multiplies
    /// outermost first (`((t1 * t2) * t3) ...`), which a parent's cached
    /// product can only reproduce for chains of one or two links, so deeper
    /// chains leave this empty and fall back to the fold.
    cumulative: Option<GpTrsf>,
    /// [`TopLocLocation::elementary_transforms`] as hash keys, outermost
    /// first; built on the first [`TopLocLocation::key`] and shared from then
    /// on, so map keys no longer rebuild the sequence per lookup.
    key: OnceLock<LocationKey>,
}

impl std::fmt::Debug for TopLocLocation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("TopLocLocation")
            .field(&self.elementary_transforms())
            .finish()
    }
}

impl TopLocLocation {
    /// Identity location (no transform).
    pub fn identity() -> Self { Self { head: None } }

    /// Location holding exactly one datum (`TopLoc_Location(theDatum)`).
    fn from_transform(t: &GpTrsf) -> Self {
        if t.form() == crate::gp::TrsfForm::Identity {
            return Self::identity();
        }
        Self {
            head: Some(Arc::new(LocNode {
                transform: t.clone(),
                parent: None,
                depth: 1,
                // The fold's first step is `gp_Trsf().Multiplied(t)`, and
                // `GpTrsf::multiply` returns `t` itself when the left side is
                // identity (`gp_Trsf.cxx:392-395`), so this is the exact value.
                cumulative: Some(t.clone()),
                key: OnceLock::new(),
            })),
        }
    }

    /// Compose a transform with an existing location (`t * parent`).
    pub fn composed(t: &GpTrsf, parent: TopLocLocation) -> Self {
        // `TopLoc_Location::Multiplied`: an identity datum is not stored, so
        // composing with one hands back the other side.
        if t.form() == crate::gp::TrsfForm::Identity {
            return parent;
        }
        let Some(parent_node) = parent.head else {
            return Self::from_transform(t);
        };
        // The fold multiplies the outermost link first, so a two-link chain is
        // `t.Multiplied(parent product)`. From three links on that grouping
        // differs from `((t1 * t2) * t3)`, so no cached product is kept.
        let (depth, cumulative) = if parent_node.depth == 1 {
            (2, parent_node.cumulative.as_ref().map(|c| t.multiplied(c)))
        } else {
            (parent_node.depth + 1, None)
        };
        Self {
            head: Some(Arc::new(LocNode {
                transform: t.clone(),
                parent: Some(parent_node),
                depth,
                cumulative,
                key: OnceLock::new(),
            })),
        }
    }

    /// Compute the cumulative transformation (product of chain).
    pub fn transformation(&self) -> GpTrsf {
        let Some(node) = &self.head else {
            return GpTrsf::identity();
        };
        if let Some(cumulative) = &node.cumulative {
            return cumulative.clone();
        }
        // Fallback for chains long enough that the cached grouping would
        // differ: the same outer-to-inner fold as before.
        let mut result = GpTrsf::identity();
        let mut current = Some(&**node);
        while let Some(link) = current {
            result = result.multiplied(&link.transform);
            current = link.parent.as_deref();
        }
        result
    }

    /// Invert the location.
    pub fn inverted(&self) -> Self {
        let Some(node) = &self.head else {
            return Self::identity();
        };
        let inv = node.transform.inverted().unwrap_or(GpTrsf::identity());
        match &node.parent {
            Some(parent) => Self::composed(
                &inv,
                Self { head: Some(parent.clone()) }.inverted(),
            ),
            None => Self::from_transform(&inv),
        }
    }

    /// Is identity (no transform)?
    pub fn is_identity(&self) -> bool {
        self.head.is_none()
    }

    /// The chain's elementary transforms, outermost link first. Identity links
    /// carry no datum in OCCT and are never stored, so this is exactly the
    /// series `TopLoc_Location::IsEqual` (`TopLoc_Location.hxx:136`) compares;
    /// it is also the comparison unit of [`TopLocLocation::key`].
    pub fn elementary_transforms(&self) -> Vec<GpTrsf> {
        let mut out = Vec::with_capacity(self.head.as_ref().map_or(0, |n| n.depth as usize));
        let mut cur = self.head.as_deref();
        while let Some(node) = cur {
            out.push(node.transform.clone());
            cur = node.parent.as_deref();
        }
        out
    }

    /// `TopLoc_Location::IsEqual` (`TopLoc_Location.hxx:136`): the two locations
    /// hold the same series of elementary transforms. Datums compare by value,
    /// which [`LocationKey`] captures exactly.
    pub fn is_equal(&self, other: &Self) -> bool {
        self.key() == other.key()
    }

    /// Hashable value form used as the second half of a `TopTools_ShapeMapHasher`
    /// key (`TopTools_ShapeMapHasher.hxx:44-52` hashes `TShape` and `Location`
    /// together; `IsSame` is `TShape` plus `Location`, `TopoDS_Shape.hxx:262-267`).
    /// Two locations produce equal keys exactly when
    /// [`TopLocLocation::is_equal`] holds, so a map keyed by
    /// `(TShape address, LocationKey)` compares like `TopoDS_Shape::IsSame`:
    /// two occurrences of one `TShape` under different locations stay distinct.
    pub fn key(&self) -> LocationKey {
        let Some(node) = &self.head else {
            return LocationKey::identity();
        };
        node.key
            .get_or_init(|| {
                let mut keys = Vec::with_capacity(node.depth as usize);
                let mut cur = Some(&**node);
                while let Some(link) = cur {
                    keys.push(TrsfKey::of(&link.transform));
                    cur = link.parent.as_deref();
                }
                LocationKey(Some(Arc::from(keys)))
            })
            .clone()
    }

    /// Transform a point through the full chain.
    pub fn transforms_point(&self, p: &GpPnt) -> GpPnt { p.transformed(&self.transformation()) }
}

impl Default for TopLocLocation { fn default() -> Self { Self::identity() } }

/// Hashable value form of a [`TopLocLocation`] ([`TopLocLocation::key`]).
/// Ordered so shape-map keys (`(TShape, LocationKey)`) can be sorted
/// deterministically; the order itself carries no OCCT meaning.
///
/// Backed by a shared slice so that copying a key (every shape-map lookup
/// copies one) does not rebuild or reallocate the elementary-transform
/// sequence. `None` is the identity location, whose sequence is empty.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LocationKey(Option<Arc<[TrsfKey]>>);

impl LocationKey {
    /// The identity location's key: an empty elementary-transform sequence.
    fn identity() -> Self {
        Self(None)
    }
}

/// Hashable value form of one elementary `gp_Trsf`.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct TrsfKey {
    scale: u64,
    form: u8,
    matrix: [[u64; 3]; 3],
    loc: [u64; 3],
}

impl TrsfKey {
    fn of(t: &GpTrsf) -> Self {
        // `-0.0` and `0.0` compare equal field-wise, so they must hash alike.
        fn bits(v: f64) -> u64 {
            if v == 0.0 {
                0.0f64.to_bits()
            } else {
                v.to_bits()
            }
        }
        let mut matrix = [[0u64; 3]; 3];
        for (i, row) in t.matrix.m.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                matrix[i][j] = bits(*v);
            }
        }
        Self {
            scale: bits(t.scale),
            form: form_code(t.shape),
            matrix,
            loc: [bits(t.loc.x), bits(t.loc.y), bits(t.loc.z)],
        }
    }
}

/// `gp_TrsfForm` discriminant, for the hash key.
fn form_code(f: crate::gp::TrsfForm) -> u8 {
    use crate::gp::TrsfForm::*;
    match f {
        Identity => 0,
        Rotation => 1,
        Translation => 2,
        PntMirror => 3,
        Ax1Mirror => 4,
        Ax2Mirror => 5,
        Scale => 6,
        CompoundTrsf => 7,
        Other => 8,
    }
}
