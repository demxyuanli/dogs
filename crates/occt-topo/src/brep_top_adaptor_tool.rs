//! `BRepTopAdaptor_Tool` (`ModelingAlgorithms/TKTopAlgo/BRepTopAdaptor/BRepTopAdaptor_Tool.hxx:27-57`,
//! `.cxx:23-105`): the `(surface, topol tool)` pair the HLR machinery keys per
//! face.
//!
//! `HLRTopoBRep_DSFiller` fills a
//! `NCollection_DataMap<TopoDS_Shape, BRepTopAdaptor_Tool, TopTools_ShapeMapHasher>`
//! and `HLRTopoBRep_OutLiner::ProcessFace` re-binds the entry of a re-built face
//! (`HLRTopoBRep_OutLiner.cxx:302`); the class itself only caches the
//! `BRepAdaptor_Surface` and the `BRepTopAdaptor_TopolTool` built on it.
//!
//! Deviations:
//! * `Tol2d` is ignored by every constructor (`cxx:31`, `:45`, `:52`, `:64`:
//!   the parameter is commented out), so it is dropped from the port.
//! * `myTopolTool` is a shared `occ::handle` in OCCT: copying the tool copies
//!   the handle, so later `Init` / `SetTopolTool` calls are seen through every
//!   copy. The port stores `Arc<TopolTool>` for the same sharing; `init_*` and
//!   `set_topol_tool` replace the `Arc` instead of mutating the object in
//!   place, which differs only for a caller still holding the previous handle.
//! * `GetTopolTool()` / `GetSurface()` print an `OCCT_DEBUG` message and return
//!   the (still null) members when `myloaded` is false (`cxx:70-102`); the port
//!   answers the same members without the print.
//! * `Destroy()` (`cxx:105`) is an empty body.

use std::sync::Arc;

use occt_geom::Surface;

use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::shape::Face;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `BRepTopAdaptor_Tool`.
#[derive(Clone)]
pub struct BrepTopAdaptorTool {
    /// `bool myloaded` (`hxx:52`).
    my_loaded: bool,
    /// `occ::handle<BRepTopAdaptor_TopolTool> myTopolTool` (`hxx:53`).
    my_topol_tool: Arc<TopolTool>,
    /// `occ::handle<Adaptor3d_Surface> myHSurface` (`hxx:54`).
    my_h_surface: Option<Arc<dyn Surface>>,
}

impl Default for BrepTopAdaptorTool {
    fn default() -> Self {
        Self::new()
    }
}

impl BrepTopAdaptorTool {
    /// `BRepTopAdaptor_Tool()` (`cxx:23-28`): a fresh topol tool, `myloaded`
    /// left false.
    pub fn new() -> Self {
        Self {
            my_loaded: false,
            my_topol_tool: Arc::new(TopolTool::new()),
            my_h_surface: None,
        }
    }

    /// `BRepTopAdaptor_Tool(F, Tol2d)` (`cxx:30-41`): `BRepAdaptor_Surface
    /// Initialize(F, true)` - i.e. the **restricted** adaptor - then
    /// `myTopolTool->Initialize(aSurf)`.
    pub fn of_face(f: &Face) -> Self {
        let adaptor = SurfaceAdaptor::of_face_restricted(f);
        let mut topol_tool = TopolTool::new();
        // `Initialize` throws `Standard_ConstructionError` for a surface that
        // carries no face (`BRepTopAdaptor_TopolTool.cxx:76-80`); the
        // restricted adaptor always carries one.
        let surface = match adaptor.surface().cloned() {
            Some(surface) => {
                let _ = topol_tool.initialize(surface.clone(), f);
                Some(surface)
            }
            None => None,
        };
        Self {
            my_loaded: true,
            my_topol_tool: Arc::new(topol_tool),
            my_h_surface: surface,
        }
    }

    /// `BRepTopAdaptor_Tool(Surface, Tol2d)` (`cxx:43-49`).
    pub fn of_surface(surface: Arc<dyn Surface>) -> Self {
        let mut topol_tool = TopolTool::new();
        topol_tool.initialize_surface_only(surface.clone());
        Self {
            my_loaded: true,
            my_topol_tool: Arc::new(topol_tool),
            my_h_surface: Some(surface),
        }
    }

    /// `Init(F, Tol2d)` (`cxx:51-59`): the **unrestricted** adaptor this time
    /// (`BRepAdaptor_Surface::Initialize(F)` without the restriction flag).
    pub fn init_face(&mut self, f: &Face) {
        let adaptor = SurfaceAdaptor::of_face(f);
        let mut topol_tool = TopolTool::new();
        self.my_h_surface = match adaptor.surface().cloned() {
            Some(surface) => {
                let _ = topol_tool.initialize(surface.clone(), f);
                Some(surface)
            }
            None => None,
        };
        self.my_topol_tool = Arc::new(topol_tool);
        self.my_loaded = true;
    }

    /// `Init(Surface, Tol2d)` (`cxx:61-66`).
    pub fn init_surface(&mut self, surface: Arc<dyn Surface>) {
        let mut topol_tool = TopolTool::new();
        topol_tool.initialize_surface_only(surface.clone());
        self.my_topol_tool = Arc::new(topol_tool);
        self.my_h_surface = Some(surface);
        self.my_loaded = true;
    }

    /// `GetTopolTool()` (`cxx:68-86`).
    pub fn topol_tool(&self) -> &TopolTool {
        &self.my_topol_tool
    }

    /// The shared handle behind `GetTopolTool()`: copying it shares the tool
    /// exactly as the OCCT `occ::handle` does.
    pub fn topol_tool_shared(&self) -> Arc<TopolTool> {
        Arc::clone(&self.my_topol_tool)
    }

    /// The non-const `GetTopolTool()` (`cxx:68-86`): a handle call in OCCT, so
    /// the tool itself is shared and mutated in place.
    pub fn topol_tool_mut(&mut self) -> &mut TopolTool {
        Arc::make_mut(&mut self.my_topol_tool)
    }

    /// `SetTopolTool(TT)` (`cxx:100-103`).
    pub fn set_topol_tool(&mut self, tt: Arc<TopolTool>) {
        self.my_topol_tool = tt;
    }

    /// `GetSurface()` (`cxx:88-102`).
    pub fn surface(&self) -> Option<&Arc<dyn Surface>> {
        self.my_h_surface.as_ref()
    }

    /// `myloaded` (`hxx:52`), the flag both getters test before answering.
    pub fn is_loaded(&self) -> bool {
        self.my_loaded
    }

    /// `Destroy()` (`cxx:105`): empty body.
    pub fn destroy(&mut self) {}
}
