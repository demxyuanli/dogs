//! OCCT visualization - scene assembly and image output on top of HLR.
//!
//! # Architecture
//! ```text
//! viz_scene     (scene graph + camera, renders through render_svg and HLR)
//! render_svg    (SVG rasterizer for projected shapes)
//! brep_scene    (BRep shape to scene-item extraction)
//! ```

pub mod brep_scene;
pub mod render_svg;
pub mod viz_scene;
