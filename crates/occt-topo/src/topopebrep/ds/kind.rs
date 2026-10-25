//! `TopOpeBRepDS_Kind`, `TopOpeBRepDS_Config` and the `TopOpeBRepDS` static
//! helpers.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_Kind.hxx:20-33`,
//! `TopOpeBRepDS_Config.hxx:19-24`, `TopOpeBRepDS.cxx:25-318`.

use crate::abs::{Orientation, ShapeType, State};

/// Different types of objects in the data structure.
/// Source: `TopOpeBRepDS_Kind.hxx:20-33`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Kind {
    Point,
    Curve,
    Surface,
    Vertex,
    Edge,
    Wire,
    Face,
    Shell,
    Solid,
    CompSolid,
    Compound,
    Unknown,
}

/// Orientation of a same-domain pair.
/// Source: `TopOpeBRepDS_Config.hxx:19-24`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Config {
    /// `TopOpeBRepDS_UNSHGEOMETRY`
    #[default]
    UnshGeometry,
    /// `TopOpeBRepDS_SAMEORIENTED`
    SameOriented,
    /// `TopOpeBRepDS_DIFFORIENTED`
    DiffOriented,
}

/// `TopOpeBRepDS::SPrint(const TopAbs_State)` (`cxx:25-45`): `IN`, `OU`, `ON`,
/// `UN`.
pub fn state_s_print(s: State) -> String {
    match s {
        State::In => "IN".to_string(),
        State::Out => "OU".to_string(),
        State::On => "ON".to_string(),
        State::Unknown => "UN".to_string(),
    }
}

/// `TopOpeBRepDS::SPrint(const TopOpeBRepDS_Kind)` (`cxx:55-97`). `UNKNOWN` is
/// the `default:` arm and prints the empty string, as in OCCT.
pub fn kind_s_print(k: Kind) -> String {
    match k {
        Kind::Point => "PO".to_string(),
        Kind::Curve => "CU".to_string(),
        Kind::Surface => "SU".to_string(),
        Kind::Vertex => "VE".to_string(),
        Kind::Edge => "ED".to_string(),
        Kind::Wire => "WI".to_string(),
        Kind::Face => "FA".to_string(),
        Kind::Shell => "SH".to_string(),
        Kind::Solid => "SO".to_string(),
        Kind::CompSolid => "CS".to_string(),
        Kind::Compound => "CO".to_string(),
        Kind::Unknown => String::new(),
    }
}

/// `TopOpeBRepDS::SPrint(K, I, S1, S2)` (`cxx:99-115`): `S1(kind, i)S2` with a
/// blank in front of a one-digit index.
pub fn kind_index_s_print(k: Kind, i: i32, before: &str, after: &str) -> String {
    let si = if (0..=9).contains(&i) { " " } else { "" };
    format!("{before}({},{}{}){after}", kind_s_print(k), si, i)
}

/// `TopOpeBRepDS::SPrint(const TopAbs_Orientation)` (`cxx:154-157`), i.e.
/// `TopAbs::ShapeOrientationToString`.
pub fn orientation_s_print(o: Orientation) -> String {
    o.to_str().to_string()
}

/// `TopOpeBRepDS::SPrint(const TopOpeBRepDS_Config)` (`cxx:161-179`).
pub fn config_s_print(c: Config) -> String {
    match c {
        Config::UnshGeometry => "UNSH".to_string(),
        Config::SameOriented => "SAME".to_string(),
        Config::DiffOriented => "DIFF".to_string(),
    }
}

/// `TopOpeBRepDS::IsTopology` (`cxx:187-218`).
pub fn is_topology(k: Kind) -> bool {
    matches!(
        k,
        Kind::Compound
            | Kind::CompSolid
            | Kind::Solid
            | Kind::Shell
            | Kind::Face
            | Kind::Wire
            | Kind::Edge
            | Kind::Vertex
    )
}

/// `TopOpeBRepDS::IsGeometry` (`cxx:222-238`).
pub fn is_geometry(k: Kind) -> bool {
    matches!(k, Kind::Surface | Kind::Curve | Kind::Point)
}

/// `TopOpeBRepDS::ShapeToKind` (`cxx:244-277`).
///
/// OCCT starts from the dummy `TopOpeBRepDS_SOLID` and throws
/// `Standard_ProgramError` in the `default:` arm, so every shape type other
/// than `TopAbs_SHAPE` maps and `TopAbs_SHAPE` is the error.
pub fn shape_to_kind(s: ShapeType) -> Result<Kind, String> {
    match s {
        ShapeType::Vertex => Ok(Kind::Vertex),
        ShapeType::Edge => Ok(Kind::Edge),
        ShapeType::Wire => Ok(Kind::Wire),
        ShapeType::Face => Ok(Kind::Face),
        ShapeType::Shell => Ok(Kind::Shell),
        ShapeType::Solid => Ok(Kind::Solid),
        ShapeType::CompSolid => Ok(Kind::CompSolid),
        ShapeType::Compound => Ok(Kind::Compound),
        ShapeType::Shape => Err("Standard_ProgramError: TopOpeBRepDS::ShapeToKind".into()),
    }
}

/// `TopOpeBRepDS::KindToShape` (`cxx:282-318`): the dummy `TopAbs_SHAPE` is
/// returned for every non-topological kind, as in OCCT.
pub fn kind_to_shape(k: Kind) -> ShapeType {
    if !is_topology(k) {
        return ShapeType::Shape;
    }
    match k {
        Kind::Vertex => ShapeType::Vertex,
        Kind::Edge => ShapeType::Edge,
        Kind::Wire => ShapeType::Wire,
        Kind::Face => ShapeType::Face,
        Kind::Shell => ShapeType::Shell,
        Kind::Solid => ShapeType::Solid,
        Kind::CompSolid => ShapeType::CompSolid,
        Kind::Compound => ShapeType::Compound,
        // `IsTopology` rejects every other kind above.
        _ => ShapeType::Shape,
    }
}

/// `TopOpeBRepDS::SPrint(const TopAbs_ShapeEnum)` (`cxx:136-139`).
///
/// The OCCT signature cannot fail (`ShapeToKind` throws there for
/// `TopAbs_SHAPE`); `None` reports that case.
pub fn shape_s_print(t: ShapeType) -> Option<String> {
    shape_to_kind(t).ok().map(kind_s_print)
}

/// `TopOpeBRepDS::SPrint(const TopAbs_ShapeEnum, const int)` (`cxx:141-144`).
pub fn shape_index_s_print(t: ShapeType, i: i32) -> Option<String> {
    shape_to_kind(t)
        .ok()
        .map(|k| kind_index_s_print(k, i, "", ""))
}
