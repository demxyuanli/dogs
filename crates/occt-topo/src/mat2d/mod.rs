//! `MAT2d` package (`TKTopAlgo/MAT2d`) and its `MAT` graph companion
//! (`TKTopAlgo/MAT`): the medial axis transform of a set of 2D contours.
//!
//! | OCCT | lines | file |
//! |---|---|---|
//! | `MAT2d_BiInt` | `MAT2d_BiInt.cxx:19-49` | `bi_int.rs` |
//! | `MAT2d_Connexion` | `MAT2d_Connexion.cxx:23-232` | `connexion.rs` |
//! | `MAT2d_CutCurve` | `MAT2d_CutCurve.cxx:27-118` | `cut_curve.rs` |
//! | `MAT2d_MiniPath` | `MAT2d_MiniPath.cxx:44-396` | `mini_path.rs` |
//! | `MAT2d_Circuit` | `MAT2d_Circuit.cxx:57-825` | `circuit/` |
//! | `MAT2d_Tool2d` | `MAT2d_Tool2d.cxx:82-1481` | `tool2d/` |
//! | `MAT2d_Mat2d` | `MAT2d_Mat2d.cxx:17-1911` | `mat2d/` |
//!
//! The `MAT` graph classes live in `crate::mat`. `BRepMAT2d` (not ported) is the
//! only entry point; `MAT2d_Tool2d` is the consumer of `Bisector_Inter`
//! (`MAT2d_Tool2d.cxx`).

mod bi_int;
mod circuit;
mod connexion;
mod cut_curve;
mod item;
mod mat2d;
mod mini_path;
mod tool2d;

pub use bi_int::Mat2dBiInt;
pub use circuit::Mat2dCircuit;
pub use connexion::Mat2dConnexion;
pub use cut_curve::Mat2dCutCurve;
pub use item::{Mat2dFigure, Mat2dItem, Mat2dLine};
pub use mat2d::Mat2dMat2d;
pub use mini_path::{ConnexionHandle, Mat2dMiniPath};
pub use tool2d::Mat2dTool2d;
