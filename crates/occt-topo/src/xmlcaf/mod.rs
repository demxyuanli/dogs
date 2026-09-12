//! XML XCAF container — an XML 1.0 assembly/attributes document.
//! Source: `XmlXCAF` (XML `XCAFDoc_*` document driver).
//!
//! A `XmlXcafDoc` is a tree of entries (`XmlEntry`), each carrying an
//! optional `TopoShape`, a list of string attributes (name / color / layer /
//! material …) and child entries — the text-format complement to the binary
//! `bincaf` container. The shape is stored topologically with the same scheme
//! `bincaf` uses, but as XML elements: type, then per-type geometry snapshots
//! (vertex points, edge curves, face surfaces with parameter ranges, and the
//! child sub-shape lists). Curves and surfaces cannot be downcast from
//! `Arc<dyn Curve>` / `Arc<dyn Surface>`, so they are classified by sampling
//! invariants (zero second derivative ⇒ line, periodic ⇒ circle, planar ⇒
//! plane, equidistant samples ⇒ sphere), exactly like `bincaf` / STEP / IGES.
//!
//! Document layout:
//! ```text
//! <?xml version="1.0" encoding="UTF-8"?>
//! <xcaf version="1.0">
//!   <entry name="...">
//!     <attribute kind="name" value="..."/>...
//!     <shape type="Solid">
//!       <shape type="Shell">
//!         <shape type="Face"><surface .../><shape type="Wire">...</shape></shape>
//!       </shape>
//!     </shape>
//!     <entry name="...">...</entry>
//!   </entry>
//! </xcaf>
//! ```
//! Attribute values are XML-escaped (`&` `<` `>` `"` `'`); the hand-rolled
//! parser unescapes them on read and rejects malformed input with `Err`.
mod prelude {

pub(crate) use std::f64::consts::PI;
pub(crate) use std::sync::Arc;

pub(crate) use occt_core::gp::{GpAx1, GpAx2, GpAx3, GpCirc, GpDir, GpPln, GpPnt, GpSphere, GpVec};
pub(crate) use occt_geom::{Curve, GeomCircle, GeomLine, GeomPlane, GeomSphere, Surface};

pub(crate) use crate::abs::ShapeType;
pub(crate) use crate::brep_surface::{classify_surface, face_plane, sphere_center, SurfaceKind};
pub(crate) use crate::builder::TopoBuilder;
pub(crate) use crate::shape::{Edge, Face, Shell, TopoShape, Wire};
pub(crate) use crate::tgeometry::GeometryRegistry;

}


mod p01;
mod p02;
pub use p01::*;
pub use p02::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
