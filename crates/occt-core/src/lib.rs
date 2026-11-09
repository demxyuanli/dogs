//! OCCT foundation facade - compatibility surface over the physical layers.
//!
//! The old monolithic `occt-core` crate is now three physical crates:
//! ```text
//! occt-foundation   TKernel   kernel, precision, message, quantity
//! occt-gp           TKMath    gp, bnd, elib, cslib, bspl, poly, gprop, ...
//! occt-geombase     TKGeomBase geom, intcurve, intf, intimpargen, adv_approx, bvh, io
//! ```
//! Every path that used to live at `occt_core::<module>` keeps resolving here.

pub use occt_foundation::{kernel, message, precision, quantity};

pub use occt_gp::{
    approx, bnd, bspl, convert, cslib, elib, gcpnts, gp, gprop, hull,
    int, intana2d, intres2d, math_bracketed_root, math_bullard, math_direct_poly_roots, math_fn,
    math_function_all_roots, math_function_roots, math_function_sample, math_newton_function_root,
    math_trig_roots, numeric, poly, toploc,
};

pub use occt_geombase::{adv_approx, bvh, geom, intcurve, intf, intimpargen, io};

pub use gp::*;
pub use kernel::*;
pub use bnd::*;
pub use precision::{
    Precision, ANGULAR, APPROXIMATION, CONFUSION, INFINITE, INTERSECTION, RESOLUTION,
};
