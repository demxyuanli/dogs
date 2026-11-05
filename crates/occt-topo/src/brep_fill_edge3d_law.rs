//! `BRepFill_Edge3DLaw` — a location law made of one 3D curve per path edge.
//!
//! Source: `BRepFill_Edge3DLaw.cxx:36-88`, `.hxx:21-41`.
//!
//! The constructor walks the path with `BRepTools_WireExplorer`, stores every
//! non-degenerated edge with a curve in `myEdges`, and gives each of them a
//! *copy* of the caller's prototype law restricted to that edge's 3D curve.
//!
//! Translation notes:
//! - `BRep_Tool::Curve(E, First, Last)` (the two-argument form, which folds the
//!   edge location into the curve) is [`BRepTool::edge_curve_world`];
//!   `First` / `Last` come from [`BRepTool::edge_parameters`].
//! - `GeomAdaptor_Curve(C, First, Last)` is the port's
//!   [`GeomTrimmedCurve`] window, which reparameterises the trim onto `[0, 1]`
//!   — the convention the whole `occt-geom` crate follows.
//! - `Geom_TrimmedCurve::Reverse()` (`cxx:63-70`) swaps the window ends
//!   (`trimmed.rs` `Curve::reverse`), so the following `FirstParameter()` /
//!   `LastParameter()` read `0` / `1` exactly like any other trimmed curve of
//!   the port.

use std::sync::Arc;

use occt_geom::curve::Curve;
use occt_geom::geom_fill::LocationLaw as GeomFillLocationLaw;
use occt_geom::trimmed::GeomTrimmedCurve;

use crate::abs::Orientation;
use crate::brep_fill_location_law::LocationLawBase;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::shape::Wire;

/// `BRepFill_Edge3DLaw` (`BRepFill_Edge3DLaw.hxx:30-38`).
pub struct Edge3DLaw {
    /// The inherited `BRepFill_LocationLaw` state.
    base: LocationLawBase,
}

impl Edge3DLaw {
    /// `BRepFill_Edge3DLaw(Path, Law)` (`cxx:36-83`).
    pub fn new(path: &Wire, law: &Arc<dyn GeomFillLocationLaw>) -> Self {
        let mut base = LocationLawBase::new();
        base.init(path);

        let mut ipath = 0usize;
        let mut wexp = WireExplorer::new();
        wexp.init(path, None);
        while wexp.more() {
            let Some(e) = wexp.current_edge() else {
                wexp.next();
                continue;
            };
            if !BRepTool::is_degenerated(&e) {
                ipath += 1;
                base.set_edge(ipath, &e);

                let Some(c) = BRepTool::edge_curve_world(&e) else {
                    wexp.next();
                    continue;
                };
                let (mut first, mut last) = BRepTool::edge_parameters(&e);
                let mut curve: Arc<dyn Curve> = c;
                if e.0.orientation() == Orientation::Reversed {
                    // `CBis = new Geom_TrimmedCurve(C, First, Last); CBis->Reverse();`
                    // (`cxx:63-66`), then the range of the reversed curve.
                    let mut cbis = GeomTrimmedCurve::new(Arc::clone(&curve), first, last);
                    cbis.reverse();
                    let cbis: Arc<dyn Curve> = Arc::new(cbis);
                    curve = cbis;
                    // `First = C->FirstParameter(); Last = C->LastParameter();`
                    first = curve.first_parameter();
                    last = curve.last_parameter();
                }

                // `AC = new GeomAdaptor_Curve(C, First, Last)` (`cxx:75`).
                let ac: Arc<dyn Curve> = Arc::new(GeomTrimmedCurve::new(
                    Arc::clone(&curve),
                    first,
                    last,
                ));
                // `myLaws->SetValue(ipath, Law->Copy());
                //  myLaws->ChangeValue(ipath)->SetCurve(AC);` (`cxx:76-77`).
                let copy = law.copy();
                let _ = copy.set_curve(ac);
                base.set_law(ipath, copy);
            }
            wexp.next();
        }

        Self { base }
    }

    /// The inherited `BRepFill_LocationLaw` state.
    pub fn base(&self) -> &LocationLawBase {
        &self.base
    }

    /// The inherited `BRepFill_LocationLaw` state, mutably.
    pub fn base_mut(&mut self) -> &mut LocationLawBase {
        &mut self.base
    }
}

/// `BRepFill_Edge3DLaw` seen through the `BRepFill_LocationLaw` base class:
/// `BRepFill_Sweep` and `BRepFill_PipeShell` hold the base handle.
impl crate::brep_fill_location_law::LocationLaw for Edge3DLaw {
    fn base(&self) -> &LocationLawBase {
        &self.base
    }

    fn base_mut(&mut self) -> &mut LocationLawBase {
        &mut self.base
    }
}
