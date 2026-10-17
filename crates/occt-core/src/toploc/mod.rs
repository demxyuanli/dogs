//! Topology location — nested datum transforms. Source: `TopLoc/`
pub mod datum;
use crate::gp::{GpTrsf, GpPnt};

/// Nested chain of transformations (like a singly-linked list of gp_Trsf).
/// Each node stores one transform and optionally points to a parent.
#[derive(Debug, Clone)]
pub struct TopLocLocation {
    transform: GpTrsf,
    parent: Option<Box<TopLocLocation>>,
}

impl TopLocLocation {
    /// Identity location (no transform).
    pub fn identity() -> Self { Self { transform: GpTrsf::identity(), parent: None } }

    /// Compose a transform with an existing location.
    pub fn composed(t: &GpTrsf, parent: TopLocLocation) -> Self {
        Self { transform: t.clone(), parent: Some(Box::new(parent)) }
    }

    /// Compute the cumulative transformation (product of chain).
    pub fn transformation(&self) -> GpTrsf {
        let mut result = GpTrsf::identity();
        let mut current = self;
        loop {
            result = (&result).multiplied(&current.transform);
            match &current.parent {
                Some(p) => current = p,
                None => break,
            }
        }
        result
    }

    /// Invert the location.
    pub fn inverted(&self) -> Self {
        let inv = self.transform.inverted().unwrap_or(GpTrsf::identity());
        match &self.parent {
            Some(p) => TopLocLocation::composed(&inv, p.inverted()),
            None => Self { transform: inv, parent: None },
        }
    }

    /// Is identity (no transform)?
    pub fn is_identity(&self) -> bool {
        self.transform.form() == crate::gp::TrsfForm::Identity && self.parent.is_none()
    }

    /// The chain's elementary transforms, outermost link first, with identity
    /// links skipped: identity links carry no datum in OCCT, so
    /// `TopLoc_Location::IsEqual` (`TopLoc_Location.hxx:136`) never sees them.
    /// This is the comparison unit of [`TopLocLocation::is_equal`] and of
    /// [`TopLocLocation::key`].
    pub fn elementary_transforms(&self) -> Vec<GpTrsf> {
        let mut out = Vec::new();
        let mut cur = Some(self);
        while let Some(node) = cur {
            if node.transform.form() != crate::gp::TrsfForm::Identity {
                out.push(node.transform.clone());
            }
            cur = node.parent.as_deref();
        }
        out
    }

    /// `TopLoc_Location::IsEqual` (`TopLoc_Location.hxx:136`): the two locations
    /// hold the same series of elementary transforms. Identity links carry no
    /// datum in OCCT, so they are skipped. Datums compare by value here, as in
    /// `datum.rs`.
    pub fn is_equal(&self, other: &Self) -> bool {
        fn trsf_eq(a: &GpTrsf, b: &GpTrsf) -> bool {
            a.scale == b.scale && a.shape == b.shape && a.matrix == b.matrix && a.loc == b.loc
        }
        let a = self.elementary_transforms();
        let b = other.elementary_transforms();
        a.len() == b.len() && a.iter().zip(b.iter()).all(|(x, y)| trsf_eq(x, y))
    }

    /// Hashable value form used as the second half of a `TopTools_ShapeMapHasher`
    /// key (`TopTools_ShapeMapHasher.hxx:44-52` hashes `TShape` and `Location`
    /// together; `IsSame` is `TShape` plus `Location`, `TopoDS_Shape.hxx:262-267`).
    /// Two locations produce equal keys exactly when
    /// [`TopLocLocation::is_equal`] holds, so a map keyed by
    /// `(TShape address, LocationKey)` compares like `TopoDS_Shape::IsSame`:
    /// two occurrences of one `TShape` under different locations stay distinct.
    pub fn key(&self) -> LocationKey {
        LocationKey(
            self.elementary_transforms()
                .iter()
                .map(TrsfKey::of)
                .collect(),
        )
    }

    /// Transform a point through the full chain.
    pub fn transforms_point(&self, p: &GpPnt) -> GpPnt { p.transformed(&self.transformation()) }
}

impl Default for TopLocLocation { fn default() -> Self { Self::identity() } }

/// Hashable value form of a [`TopLocLocation`] ([`TopLocLocation::key`]).
/// Ordered so shape-map keys (`(TShape, LocationKey)`) can be sorted
/// deterministically; the order itself carries no OCCT meaning.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LocationKey(Vec<TrsfKey>);

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
