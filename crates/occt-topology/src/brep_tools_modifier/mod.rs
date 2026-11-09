//! `BRepTools_Modification` and `BRepTools_Modifier` (TKBRep).
//!
//! Sources:
//! - `BRepTools_Modification.hxx:41-145` / `BRepTools_Modification.cxx:20-40`
//!   — the abstract description of a geometric modification.
//! - `BRepTools_Modifier.hxx:53-176` / `BRepTools_Modifier.cxx:1-852` — the
//!   engine that replays a modification over a whole shape tree.
//! - `BRepTools_TrsfModification.hxx:39-150` / `.cxx:1-486` — the `gp_Trsf`
//!   driven modification, kept as the reference implementation of the trait.
//!
//! The engine is one recursive walk (`Modifier::rebuild`) over the child tree
//! of the input shape: every shape is answered from `myMap`, which `Put`
//! prefills with a bound-but-null entry for the whole tree. A shape is rebuilt
//! when the modification asks for new geometry for it, or when any of its
//! children was rebuilt; a face whose surface changed also rewrites the
//! pcurves of the rebuilt edges through `NewCurve2d`
//! (`BRepTools_Modifier.cxx:397-587`).

mod modification;
mod modifier;
mod new_info;
mod rebuild;
mod trsf_modification;

pub use modification::*;
pub use modifier::*;
pub use trsf_modification::*;
