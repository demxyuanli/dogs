//! OCCT `IntPolyh` package (`TKBRep` / `TKGeomAlgo`): polygonal intersection
//! of two surfaces, ported from `src/ModelingAlgorithms/TKGeomAlgo/IntPolyh`.
//!
//! Modules follow the `.hxx`/`.cxx` split of the package:
//! `point` / `edge_couple` / `triangle` / `start_point` / `section_line` are
//! the data classes, `arrays` is `IntPolyh_Array<Type>`, `tools` holds the
//! static helpers and `maillage_affinage` is `IntPolyh_MaillageAffinage`.

pub mod arrays;
pub mod edge_couple;
pub mod intersection;
pub mod maillage_affinage;
pub mod point;
pub mod section_line;
pub mod start_point;
pub mod tools;
pub mod triangle;

pub use arrays::{
    ArrayOfEdges, ArrayOfPoints, ArrayOfPointNormal, ArrayOfSectionLines, ArrayOfTangentZones,
    ArrayOfTriangles, IntPolyhArray, IntPolyhPointNormal,
};
pub use edge_couple::{IntPolyhCouple, IntPolyhEdge};
pub use intersection::{IntPolyhIntersection, IntPolyhLinePoint, IntPolyhTangentPoint};
pub use maillage_affinage::MaillageAffinage;
pub use point::IntPolyhPoint;
pub use section_line::IntPolyhSectionLine;
pub use start_point::IntPolyhStartPoint;
pub use triangle::IntPolyhTriangle;
