//! `PerformBoolean` (`BRepFill_AdvancedEvolved.cxx:74-100`) and `IsLid`
//! (`cxx:375-404`).

use crate::bop_builder2::BopBuilder;
use crate::brep_tool::BRepTool;
use crate::shape::{Face, TopoShape};

use super::AdvancedEvolved;

impl AdvancedEvolved {
    /// `PerformBoolean(theArgsList, theResult)`
    /// (`BRepFill_AdvancedEvolved.cxx:74-100`): the general fuse of the
    /// arguments through a `BOPAlgo_PaveFiller` + `BOPAlgo_Builder` pair.
    ///
    /// Answers `None` where the reference answers `false` (`cxx:82-84`,
    /// `cxx:92-95`); `theResult` is then left untouched, exactly as the
    /// reference leaves its output argument alone.
    pub(super) fn perform_boolean(&self, args: &[TopoShape]) -> Option<TopoShape> {
        if args.is_empty() {
            // `BOPAlgo_PaveFiller::Init` reports "too few arguments".
            return None;
        }

        // `cxx:76-86`: `BOPAlgo_PaveFiller`, arguments + parallel + fuzzy.
        // The port keeps the filler of the builder it is going to feed, which
        // is the same object `PerformWithFiller` consumes.
        let mut builder = BopBuilder::new();
        for a in args {
            builder.add_argument(a);
        }
        builder.filler_mut().set_arguments(args);
        builder.filler_mut().set_fuzzy_value(self.fuzzy_value);
        if builder.filler_mut().perform().is_err() {
            return None;
        }
        if builder.filler().has_errors() {
            return None;
        }

        // `cxx:88-95`: `BOPAlgo_Builder` over the same arguments, then
        // `PerformWithFiller(aPF)`.
        match builder.perform_with_filler() {
            Ok(shape) if !builder.has_errors() => Some(shape),
            _ => None,
        }
    }

    /// `IsLid(theF, theMapOfLids)` (`cxx:375-404`): whether `theF` lies on the
    /// same plane as one of the lids.
    ///
    /// The reference compares the `Geom_Plane` handles of the two faces
    /// (`aPlane == aPlnF`, `cxx:398`); the port keeps one `Arc<dyn Surface>`
    /// per `Geom_Surface`, so handle equality is `Arc::ptr_eq`.
    pub(super) fn is_lid(f: &Face, lids: &[TopoShape]) -> bool {
        // `cxx:378-382`.
        if lids.is_empty() {
            return false;
        }

        // `cxx:384-388`: null when the face is not on a `Geom_Plane`.
        let Some(surf_f) = BRepTool::face_surface(f) else {
            return false;
        };
        if surf_f.gp_pln().is_none() {
            return false;
        }

        // `cxx:391-401`.
        for lid in lids {
            let Some(surf) = BRepTool::face_surface(&Face(lid.clone())) else {
                continue;
            };
            if surf.gp_pln().is_none() {
                continue;
            }
            if std::sync::Arc::ptr_eq(&surf, &surf_f) {
                return true;
            }
        }

        false
    }
}
