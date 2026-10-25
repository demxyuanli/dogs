//! `TopOpeBRepDS_Surface` — a geometric surface and a tolerance.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_Surface.hxx:26-56`,
//! `TopOpeBRepDS_Surface.cxx:22-48`.

use std::sync::Arc;

use occt_geom::Surface as GeomSurface;

/// Source: `TopOpeBRepDS_Surface.hxx:49-53`
#[derive(Clone, Default)]
pub struct Surface {
    /// `mySurface` (`occ::handle<Geom_Surface>`; a null handle is `None`).
    surface: Option<Arc<dyn GeomSurface>>,
    tolerance: f64,
    keep: bool,
}

impl std::fmt::Debug for Surface {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TopOpeBRepDS_Surface")
            .field("is_null", &self.surface.is_none())
            .field("tolerance", &self.tolerance)
            .field("keep", &self.keep)
            .finish()
    }
}

impl Surface {
    /// `TopOpeBRepDS_Surface()` (`cxx:24-27`): a null handle, zero tolerance,
    /// `myKeep = false`.
    pub fn new() -> Self {
        Self {
            surface: None,
            tolerance: 0.0,
            keep: false,
        }
    }

    /// `TopOpeBRepDS_Surface(const handle<Geom_Surface>& P, const double T)`
    /// (`cxx:31-36`).
    pub fn with_tolerance(surface: Arc<dyn GeomSurface>, tolerance: f64) -> Self {
        Self {
            surface: Some(surface),
            tolerance,
            keep: false,
        }
    }

    /// `Assign(Other)` (`cxx:43-48`).
    pub fn assign(&mut self, other: &Surface) {
        self.surface = other.surface.clone();
        self.tolerance = other.tolerance;
        self.keep = other.keep;
    }

    /// `Surface()` (`hxx:44-47`).
    pub fn surface(&self) -> Option<&Arc<dyn GeomSurface>> {
        self.surface.as_ref()
    }

    /// `Tolerance()` (`hxx:49-52`).
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// `Tolerance(theTol)` (`hxx:55`).
    pub fn set_tolerance(&mut self, tol: f64) {
        self.tolerance = tol;
    }

    /// `Keep()` (`hxx:57`).
    pub fn keep(&self) -> bool {
        self.keep
    }

    /// `ChangeKeep(theToKeep)` (`hxx:59`).
    pub fn change_keep(&mut self, to_keep: bool) {
        self.keep = to_keep;
    }
}
