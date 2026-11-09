//! `BRepOffset_Offset::Init(Edge, Offset)` (`BRepOffset_Offset.cxx:1698-1724`):
//! the pipe swept along a whole edge, used for a free border on which the
//! rolling ball of `BiTgte_Blend` runs.
//!
//! Control flow, all from the `.cxx`:
//!
//! * `:1700-1702` `myShape = Edge`, `myOffset = std::abs(Offset)`.
//! * `:1704-1711` `CP = Geom_TrimmedCurve(BRep_Tool::Curve(Edge, Loc, f, l),
//!   f, l)` then `CP->Transform(Loc.Transformation())` - the location is baked
//!   into the curve, exactly like [`super::BRepOffsetOffset::init_pipe`].
//! * `:1713-1718` `GeomFill_Pipe Pipe(CP, myOffset)` + `Perform()`; a failed
//!   pipe raises `Standard_ConstructionError`.
//! * `:1720-1721` `BRepLib_MakeFace MF(Pipe.Surface(), Precision::Confusion())`
//!   and `myFace = MF.Face()`.
//! * `:1722-1724` `if (Offset < 0.) myFace.Reverse()`.
//!
//! Deviation: OCCT never sets `myStatus` here, so it keeps whatever value the
//! default constructor left (indeterminate); the port keeps `Good`, the same
//! choice [`super::BRepOffsetOffset::default`] documents.

use occt_core::precision::CONFUSION;
use occt_geom::geom_fill::pipe::Pipe;

use crate::abs::Orientation;
use crate::brep_lib_make_face::make_face_from_surface;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};

use super::init_pipe::trimmed_world;
use super::BRepOffsetOffset;

impl BRepOffsetOffset {
    /// `BRepOffset_Offset::Init(const TopoDS_Edge&, const double)`
    /// (`cxx:1698-1724`).
    pub fn init_edge(&mut self, edge: &Edge, offset: f64) -> Result<(), String> {
        // `:1700-1702`.
        self.my_shape = Some(edge.0.clone());
        let my_offset = offset.abs();

        // `:1704-1711`.
        let (f, l) = BRepTool::edge_parameters(edge);
        let cp = BRepTool::edge_curve(edge)
            .ok_or("Standard_NullObject: BRepOffset_Offset::Init(Edge): no 3D curve")?;
        let cp = trimmed_world(edge, &cp, f, l);
        // `GeomFill_Pipe(Path, Radius)` only handles a line or a circle path
        // when the radius is constant; the port answers `Err` there, which maps
        // to the `Standard_ConstructionError` of `:1716-1717`.
        let failure = "Standard_ConstructionError: GeomFill_Pipe : Cannot make a surface".to_string();
        let mut pipe = Pipe::with_path_radius(cp, my_offset).map_err(|_| failure.clone())?;
        // `Perform()` with the header defaults `WithParameters = true`,
        // `Polynomial = false` (`GeomFill_Pipe.hxx`).
        pipe.perform(true, false).map_err(|_| failure.clone())?;
        if !pipe.is_done() {
            return Err(failure);
        }

        // `:1720-1721`.
        let surf = pipe.surface().ok_or(failure)?;
        let mut face = make_face_from_surface(surf, CONFUSION);

        // `:1722-1724`.
        if offset < 0.0 {
            let ori = if face.0.orientation() == Orientation::Forward {
                Orientation::Reversed
            } else {
                Orientation::Forward
            };
            face = Face(face.0.oriented(ori));
        }

        self.my_face = Some(face);
        Ok(())
    }
}
