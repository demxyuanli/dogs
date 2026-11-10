//! Port of `Intf_Tool` (`src/ModelingAlgorithms/TKGeomAlgo/Intf/Intf_Tool.hxx` +
//! `Intf_Tool.cxx`).
//!
//! `Intf_Tool` builds the part of a bounding box cut by an *infinite* line,
//! hyperbola or parabola: for each "context" (a `Bnd_Box` / `Bnd_Box2d` domain)
//! it returns the clipped box plus the parameter ranges of the curve segments
//! that run through that domain (`NbSegments`, `BeginParam`, `EndParam`).
//!
//! Ported members (OCCT 8.0.0, tag V8_0_0):
//! * `Intf_Tool()`                        `Intf_Tool.cxx:38-48`
//! * `Lin2dBox`                           `Intf_Tool.cxx:52-206`
//! * `Hypr2dBox`                          `Intf_Tool.cxx:210-369`
//! * `Inters2d(gp_Hypr2d, Bnd_Box2d)`     `Intf_Tool.cxx:373-473`
//! * `Parab2dBox`                         `Intf_Tool.cxx:477-638`
//! * `Inters2d(gp_Parab2d, Bnd_Box2d)`    `Intf_Tool.cxx:642-742`
//! * `LinBox`                             `Intf_Tool.cxx:746-963`
//! * `NbSegments`                         `Intf_Tool.cxx:1634-1637`
//! * `BeginParam`                         `Intf_Tool.cxx:1641-1645`
//! * `EndParam`                           `Intf_Tool.cxx:1649-1653`
//!
//! UNPORTED members:
//! * `HyprBox`              `Intf_Tool.cxx:967-1119`
//! * `Inters3d(gp_Hypr)`    `Intf_Tool.cxx:1123-1302`
//! * `Inters3d(gp_Parab)`   `Intf_Tool.cxx:1306-1485`
//! * `ParabBox`             `Intf_Tool.cxx:1489-1630`
//!
//! Reason: all four call `IntAna_IntConicQuad` (`Intf_Tool.cxx:1133, 1156,
//! 1179, 1202, 1225, 1248, 1316, 1339, 1362, 1385, 1408, 1431`) together with
//! the 3D `ElCLib::D1(gp_Hypr/gp_Parab, ...)` (`Intf_Tool.cxx:1015,
//! 1534`). Neither `IntAna_IntConicQuad` nor the 3D conic
//! `ElCLib::D1`/`Value` overloads exist in this crate yet (only
//! `IntAna2d` and the 2D `elib::clib2d` overloads are ported), so no
//! faithful body can be written. They are deliberately absent rather than
//! approximated.
//!
//! Note: the surrounding `intf::mod` prose mentions `LinTetra` /
//! `LinSphere`, but neither symbol exists in this class - the real member
//! list is the ten entries above (`Intf_Tool.hxx:42-66`).
mod helpers;
mod hyperbola2d;
mod lin2d;
mod lin3d;
mod parabola2d;
mod segments;
mod structure;

pub use structure::IntfTool;
