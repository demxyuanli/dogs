//! Port of ShapeFix_ComposeShell (ShapeFix_ComposeShell.cxx, 3603 lines) and
//! its auxiliary classes.

mod break_wires;
mod collect_wires;
mod composite_surface;
mod dispatch_wires;
mod helpers;
mod load_wires;
mod make_faces_on_patch;
mod perform;
mod reshape;
mod shell;
mod split_by_grid;
mod split_by_line;
mod split_wire;
mod wire_data;
mod wire_segment;

pub use composite_surface::{CompositeSurface, Parametrisation};
pub use helpers::*;
pub use load_wires::load_wires;
pub use reshape::{apply_context, IdentityReShape, MapReShape, ReShape};
pub use shell::ComposeShell;
pub use wire_data::{reverse_wire_data, reverse_wire_data_on_face};
pub use wire_segment::WireSegment;
