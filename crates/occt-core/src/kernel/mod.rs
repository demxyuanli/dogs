//! TKernel replacements — error types, handle protocol, collection aliases.
//!
//! These replace the OCCT classes that are purely infrastructure (no geometry logic):
//! - `Standard_Failure` + subclass exceptions → `OCCError` enum
//! - `opencascade::handle<T>` → `Arc<T>` (stdlib)
//! - `NCollection_*` templates → Rust stdlib types + thin wrappers

pub mod error;
pub mod handle;
pub mod containers;
