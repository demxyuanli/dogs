//! `IntCurveSurface_TheInterferenceOfHInter`.
//!
//! Source: `IntCurveSurface_TheInterferenceOfHInter.hxx` (TKGeomAlgo). The
//! class is the `Polygon3d = ThePolygonOfHInter` /
//! `Polyhedron = ThePolyhedronOfHInter` instantiation of
//! `Intf_InterferencePolygonPolyhedron` (`Intf_InterferencePolygonPolyhedron.gxx`),
//! which is already ported in `occt_core::intf`. This module binds the two
//! tools and exposes the `Perform` overloads `HInter` uses.

use occt_core::gp::GpLin;
use occt_core::intf::{
    IntfInterference, IntfInterferencePolygonPolyhedron, IntfPolyhGrid,
};

use super::polygon::{ThePolygonOfHInter, ThePolygonToolOfHInter};
use super::polyhedron::{ThePolyhedronOfHInter, ThePolyhedronToolOfHInter};

/// `IntCurveSurface_TheInterferenceOfHInter` (`...hxx:37-152`).
pub struct TheInterferenceOfHInter {
    base: IntfInterferencePolygonPolyhedron<ThePolygonToolOfHInter, ThePolyhedronToolOfHInter>,
}

impl Default for TheInterferenceOfHInter {
    fn default() -> Self {
        Self::new()
    }
}

impl TheInterferenceOfHInter {
    /// `IntCurveSurface_TheInterferenceOfHInter()` (`...hxx:43`).
    pub fn new() -> Self {
        Self {
            base: IntfInterferencePolygonPolyhedron::new(),
        }
    }

    /// `IntCurveSurface_TheInterferenceOfHInter(thePolyg, thePolyh)`
    /// (`...hxx:48-50`).
    pub fn of_polygon_polyhedron(
        the_polyg: &ThePolygonOfHInter,
        the_polyh: &ThePolyhedronOfHInter,
    ) -> Self {
        Self {
            base: IntfInterferencePolygonPolyhedron::with_polygon_polyhedron(the_polyg, the_polyh),
        }
    }

    /// `IntCurveSurface_TheInterferenceOfHInter(theLin, thePolyh)`
    /// (`...hxx:54-56`).
    pub fn of_lin_polyhedron(the_lin: &GpLin, the_polyh: &ThePolyhedronOfHInter) -> Self {
        Self {
            base: IntfInterferencePolygonPolyhedron::with_lin_polyhedron(the_lin, the_polyh),
        }
    }

    /// `IntCurveSurface_TheInterferenceOfHInter(theLins, thePolyh)`
    /// (`...hxx:60-62`).
    pub fn of_lins_polyhedron(the_lins: &[GpLin], the_polyh: &ThePolyhedronOfHInter) -> Self {
        Self {
            base: IntfInterferencePolygonPolyhedron::with_lins_polyhedron(the_lins, the_polyh),
        }
    }

    /// `IntCurveSurface_TheInterferenceOfHInter(thePolyg, thePolyh, theBoundSB)`
    /// (`...hxx:80-83`).
    pub fn of_polygon_polyhedron_grid(
        the_polyg: &ThePolygonOfHInter,
        the_polyh: &ThePolyhedronOfHInter,
        the_bsb: &IntfPolyhGrid,
    ) -> Self {
        Self {
            base: IntfInterferencePolygonPolyhedron::with_polygon_polyhedron_grid(
                the_polyg, the_polyh, the_bsb,
            ),
        }
    }

    /// `IntCurveSurface_TheInterferenceOfHInter(theLin, thePolyh, theBoundSB)`
    /// (`...hxx:87-90`).
    pub fn of_lin_polyhedron_grid(
        the_lin: &GpLin,
        the_polyh: &ThePolyhedronOfHInter,
        the_bsb: &IntfPolyhGrid,
    ) -> Self {
        Self {
            base: IntfInterferencePolygonPolyhedron::with_lin_polyhedron_grid(
                the_lin, the_polyh, the_bsb,
            ),
        }
    }

    /// `IntCurveSurface_TheInterferenceOfHInter(theLins, thePolyh, theBoundSB)`
    /// (`...hxx:94-97`).
    pub fn of_lins_polyhedron_grid(
        the_lins: &[GpLin],
        the_polyh: &ThePolyhedronOfHInter,
        the_bsb: &IntfPolyhGrid,
    ) -> Self {
        Self {
            base: IntfInterferencePolygonPolyhedron::with_lins_polyhedron_grid(
                the_lins, the_polyh, the_bsb,
            ),
        }
    }

    /// `Perform(thePolyg, thePolyh)` (`...hxx:66-68`).
    pub fn perform_polygon_polyhedron(
        &mut self,
        the_polyg: &ThePolygonOfHInter,
        the_polyh: &ThePolyhedronOfHInter,
    ) {
        self.base.perform_polygon_polyhedron(the_polyg, the_polyh);
    }

    /// `Perform(theLin, thePolyh)` (`...hxx:71-72`).
    pub fn perform_lin_polyhedron(&mut self, the_lin: &GpLin, the_polyh: &ThePolyhedronOfHInter) {
        self.base.perform_lin_polyhedron(the_lin, the_polyh);
    }

    /// `Perform(theLins, thePolyh)` (`...hxx:76-77`).
    pub fn perform_lins_polyhedron(&mut self, the_lins: &[GpLin], the_polyh: &ThePolyhedronOfHInter) {
        self.base.perform_lins_polyhedron(the_lins, the_polyh);
    }

    /// `Perform(thePolyg, thePolyh, theBoundSB)` (`...hxx:101-103`).
    pub fn perform_polygon_polyhedron_grid(
        &mut self,
        the_polyg: &ThePolygonOfHInter,
        the_polyh: &ThePolyhedronOfHInter,
        the_bsb: &IntfPolyhGrid,
    ) {
        self.base
            .perform_polygon_polyhedron_grid(the_polyg, the_polyh, the_bsb);
    }

    /// `Perform(theLin, thePolyh, theBoundSB)` (`...hxx:107-112`).
    pub fn perform_lin_polyhedron_grid(
        &mut self,
        the_lin: &GpLin,
        the_polyh: &ThePolyhedronOfHInter,
        the_bsb: &IntfPolyhGrid,
    ) {
        self.base.perform_lin_polyhedron_grid(the_lin, the_polyh, the_bsb);
    }

    /// `Perform(theLins, thePolyh, theBoundSB)` (`...hxx:114-116`).
    pub fn perform_lins_polyhedron_grid(
        &mut self,
        the_lins: &[GpLin],
        the_polyh: &ThePolyhedronOfHInter,
        the_bsb: &IntfPolyhGrid,
    ) {
        self.base
            .perform_lins_polyhedron_grid(the_lins, the_polyh, the_bsb);
    }

    /// The `Intf_Interference` base, for `NbSectionPoints` / `PntValue` /
    /// `NbTangentZones` / `ZoneValue`.
    pub fn base(&self) -> &IntfInterference {
        &self.base.base
    }
}
