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
    pub fn to_str(&self) -> &'static str {
        match self { Self::Forward=>"FORWARD", Self::Reversed=>"REVERSED", Self::Internal=>"INTERNAL", Self::External=>"EXTERNAL" }
    }
}

/// Topology flag: free, modified, check, oriented, closed, infinite. Source: `TopAbs_ShapeEnum` flags
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ShapeFlags {
    pub free: bool,
    pub modified: bool,
    pub check: bool,
    pub oriented: bool,
    pub closed: bool,
    pub infinite: bool,
    pub convex: bool,
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
}
