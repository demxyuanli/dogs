//! OCCT geometry-base layer - pure Rust port of Open CASCADE Technology TKGeomBase.
//!
//! Depends on `occt-gp`. Nothing here may reference topology.
//!
//! # Architecture
//! ```text
//! geom/              (Geom_* curves and surfaces, GeomAbs, adaptors)
//! intcurve/          (IntCurve adaptors)
//! intf/ intimpargen/ (curve/curve and implicit/parametric intersection bases)
//! adv_approx/        (AdvApprox approximation)
//! bvh/               (BVH over triangulations)
//! io/                (Poly/STL/PLY triangulation readers and writers)
//! ```

pub mod adv_approx;
pub mod bvh;
pub mod geom;
pub mod intcurve;
pub mod intf;
pub mod intimpargen;
pub mod io;
