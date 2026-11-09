//! Handle protocol — replaces OCCT's `Handle(T)` / `opencascade::handle<T>`.
//!
//! ## Rust equivalent
//!
//! OCCT uses an intrusive smart pointer (`opencascade::handle<T>`) that requires
//! `T` to inherit from `Standard_Transient`. The handle is reference-counted,
//! supports null, and provides RTTI via `DynamicType()` / `IsKind()`.
//!
//! Rust's `std::sync::Arc<T>` replaces this directly:
//! - Reference counting is automatic (no base class needed)
//! - Null is represented by `Option<Arc<T>>`
//! - RTTI is replaced by trait objects (`Arc<dyn Trait>`) or enum dispatch
//! - Thread-safety: `Arc` is `Send + Sync` when `T` is
//!
//! ## Migration rules
//!
//! | OCCT | Rust |
//! |------|------|
//! | `Handle(Geom_Curve)` | `Arc<dyn Curve>` or `Arc<GeomCurve>` |
//! | `Handle.IsNull()` | `Option<Arc<T>>::is_none()` |
//! | `Handle.DownCast<Geom_Line>()` | `Arc::downcast::<GeomLine>` (if using Any) or match on enum |
//! | `Handle = new Geom_Line(...)` | `Arc::new(GeomLine::new(...))` |
//! | `Standard_Transient` base class | Not needed — `Arc` handles refcounting |
//!
//! ## When porting
//!
//! - If OCCT type inherits `Standard_Transient` → wrap in `Arc<T>` at use site
//! - If OCCT uses `DynamicType()` / `IsKind()` → use Rust enum + match, or `dyn Trait` + `Any`
//! - If OCCT stores `Handle` as member → store `Arc<T>` directly
//!
//! Ponteil note: No custom Handle type is needed. `Arc` is standard, well-tested,
//! and already integrated with Rust's type system. Adding a wrapper adds
//! indirection without benefit.

/// Re-export the standard library types that replace OCCT Handles.
pub use std::sync::Arc;

/// Alias for OCCT's `opencascade::handle<T>` — direct mapping to `Arc<T>`.
pub type Handle<T> = Arc<T>;
