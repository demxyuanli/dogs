//! OCCT data exchange - pure Rust port of the STEP / IGES / XCAF / RWMesh layers.
//!
//! Sits directly above the modelling core. The BRep construction objects
//! (`BRepTools_ShapeSet` and friends) stay in `occt-topo`; this crate only owns
//! file-format reading and writing plus the document model those formats share.
//!
//! # Architecture
//! ```text
//! step/ iges/ iges_read/   (STEP and IGES readers and writers)
//! xmlcaf/ bincaf/ xcaf*/   (XCAF document, OCAF-lite storage and styles)
//! model/                   (BRepModel document model shared by all formats)
//! brep_exchange, image_texture, tdata_std_named_data
//! rwmesh/ gltf vrml        (mesh read/write formats)
//! brep_assembly            (assembly helpers that depend on the STEP reader)
//! ```

pub mod bincaf;
pub mod brep_assembly;
pub mod brep_exchange;
pub mod gltf;
pub mod iges;
pub mod iges_read;
pub mod image_texture;
pub mod model;
pub mod rwmesh;
pub mod step;
pub mod tdata_std_named_data;
pub mod vrml;
pub mod xcaf;
pub mod xcaf_doc_vis_material;
pub mod xcaf_prs;
pub mod xcaf_prs_document_explorer;
pub mod xmlcaf;
