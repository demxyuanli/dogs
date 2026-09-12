//! Topological shape types and orientations. Source: `TopAbs_ShapeEnum`, `TopAbs_Orientation`

/// Kind of topological shape. Source: `TopAbs_ShapeEnum`
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ShapeType {
    Compound,
    CompSolid,
    Solid,
    Shell,
    Face,
    Wire,
    Edge,
    Vertex,
    Shape,
}

impl ShapeType {
    /// Human-readable name (matches OCCT TopAbs::ShapeTypeToString).
    pub fn to_str(&self) -> &'static str {
        match self {
            Self::Compound => "Compound",
            Self::CompSolid => "CompSolid",
            Self::Solid => "Solid",
            Self::Shell => "Shell",
            Self::Face => "Face",
            Self::Wire => "Wire",
            Self::Edge => "Edge",
            Self::Vertex => "Vertex",
            Self::Shape => "Shape",
        }
    }

    /// Parse from string.
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "Compound" => Some(Self::Compound),
            "CompSolid" => Some(Self::CompSolid),
            "Solid" => Some(Self::Solid),
            "Shell" => Some(Self::Shell),
            "Face" => Some(Self::Face),
            "Wire" => Some(Self::Wire),
            "Edge" => Some(Self::Edge),
            "Vertex" => Some(Self::Vertex),
            "Shape" => Some(Self::Shape),
            _ => None,
        }
    }
}

/// Orientation of a shape within its parent. Source: `TopAbs_Orientation`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    Forward,
    Reversed,
    Internal,
    External,
}

impl Orientation {
    pub fn is_forward(&self) -> bool { *self == Self::Forward }
    pub fn is_reversed(&self) -> bool { *self == Self::Reversed }
    pub fn reversed(&self) -> Self {
        match self {
            Self::Forward => Self::Reversed,
            Self::Reversed => Self::Forward,
            Self::Internal | Self::External => *self,
        }
    }

    /// `TopAbs::Compose(parent, child)` — accumulated orientation of a child
    /// inside a parent. Table is indexed `[child][parent]` to match OCCT
    /// `TopAbs.hxx` (`aTable[Or2][Or1]`).
    pub fn compose(parent: Self, child: Self) -> Self {
        const T: [[Orientation; 4]; 4] = [
            [Orientation::Forward, Orientation::Reversed, Orientation::Internal, Orientation::External],
            [Orientation::Reversed, Orientation::Forward, Orientation::Internal, Orientation::External],
            [Orientation::Internal, Orientation::Internal, Orientation::Internal, Orientation::Internal],
            [Orientation::External, Orientation::External, Orientation::External, Orientation::External],
        ];
        T[child.as_index()][parent.as_index()]
    }

    fn as_index(self) -> usize {
        match self {
            Self::Forward => 0,
            Self::Reversed => 1,
            Self::Internal => 2,
            Self::External => 3,
        }
    }

    pub fn to_str(&self) -> &'static str {
        match self { Self::Forward=>"FORWARD", Self::Reversed=>"REVERSED", Self::Internal=>"INTERNAL", Self::External=>"EXTERNAL" }
    }
}

/// Topology flag: free, modified, check, oriented, closed, infinite. Source: `TopoDS_TShape` bits
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ShapeFlags {
    pub free: bool,
    pub modified: bool,
    pub check: bool,
    pub oriented: bool,
    pub closed: bool,
    pub infinite: bool,
    pub convex: bool,
}

impl Default for ShapeFlags {
    /// OCCT `TopoDS_TShape` ctor: Free | Modified | Orientable.
    fn default() -> Self {
        Self {
            free: true,
            modified: true,
            check: false,
            oriented: true,
            closed: false,
            infinite: false,
            convex: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shape_type_roundtrip() {
        for t in [ShapeType::Vertex, ShapeType::Edge, ShapeType::Face, ShapeType::Solid, ShapeType::Compound] {
            assert_eq!(ShapeType::from_str(t.to_str()), Some(t));
        }
    }

    #[test]
    fn orientation_reverse() {
        assert_eq!(Orientation::Forward.reversed(), Orientation::Reversed);
        assert_eq!(Orientation::Reversed.reversed(), Orientation::Forward);
        assert_eq!(Orientation::Internal.reversed(), Orientation::Internal);
    }

    #[test]
    fn orientation_compose_matches_topabs() {
        assert_eq!(Orientation::compose(Orientation::Forward, Orientation::Reversed), Orientation::Reversed);
        assert_eq!(Orientation::compose(Orientation::Reversed, Orientation::Reversed), Orientation::Forward);
        assert_eq!(Orientation::compose(Orientation::Internal, Orientation::Forward), Orientation::Internal);
        assert_eq!(Orientation::compose(Orientation::External, Orientation::Reversed), Orientation::External);
    }

    #[test]
    fn tshape_flag_defaults_match_occt() {
        let f = ShapeFlags::default();
        assert!(f.free && f.modified && f.oriented);
        assert!(!f.check && !f.closed && !f.infinite && !f.convex);
    }
}
