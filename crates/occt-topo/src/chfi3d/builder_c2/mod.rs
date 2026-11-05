//! The `ChFi3d_Builder_C2.cxx` file statics.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_C2.cxx:63-99`
//! (OCCT 8.0.0): the two `Reduce` overloads of `PerformTwoCornerbyInter`.
//! The member bodies of that file are `PerformTwoCornerbyInter`
//! (`:108-806`, ported as
//! [`crate::chfi3d::builder::ChFi3dBuilder::perform_two_corner_by_inter`]) and
//! `UpdateTolesp` (`:814-828`, ported as
//! [`crate::chfi3d::builder::state::update_tolesp`]).
//!
//! `Reduce` re-cuts the two surface adaptors (or the one curve adaptor) at the
//! middle of their range so that a `(p1, p2)` that both lie on the *same*
//! surface no longer share it: the two halves then belong to two disjoint
//! parameter rectangles, which is what the following intersection needs.
//!
//! Translation notes:
//! - `ChFi3d_Builder_C2.cxx` declares both as file statics of type
//!   `void Reduce(double, double, Handle&...)`; the port takes the adaptor by
//!   `&mut` because the OCCT body mutates the adaptor object the handle points
//!   at and the callers own their adaptors by value,
//! - `surf->Bounds(ud, uf, vd, vf)` is the *underlying* `Geom_Surface`'s own
//!   range (not the adaptor's trim), and `s.Load(surf, ud, uf, v1, v2)` is
//!   [`crate::topopebrep::tool::adaptor::SurfaceAdaptor::of_surface_with_range`]
//!   over the same basis surface.

mod reduce;

pub use reduce::{reduce_curve, reduce_surface};
