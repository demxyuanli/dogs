//! OCCT facade - one import path for the whole kernel.
//!
//! Every layer is re-exported, so callers can depend on `occt` alone:
//! ```text
//! occt::core       occt-core      foundation, gp, geombase
//! occt::topo       occt-topo      TKBRep / TKTopAlgo / algorithms / meshing
//! occt::hlr        occt-hlr       hidden line removal
//! occt::exchange   occt-exchange  STEP / IGES / XCAF / RWMesh
//! occt::visual     occt-visual    scenes and SVG output
//! ```
//! `draw` is the Draw-style command surface built on all of the above.

pub mod draw;

pub use occt_core as core;
pub use occt_exchange as exchange;
pub use occt_hlr as hlr;
pub use occt_topo as topo;
pub use occt_visual as visual;
