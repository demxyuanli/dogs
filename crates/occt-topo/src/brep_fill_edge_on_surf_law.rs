//! `BRepFill_EdgeOnSurfLaw` — a location law whose frames are the Darboux
//! trihedra of the path's pcurves on a support face.
//!
//! Source: `BRepFill_EdgeOnSurfLaw.cxx:40-105`, `.hxx:24-35`.
//!
//! For every non-degenerated path edge the constructor looks for the first face
//! of `Surf` that carries a pcurve of the edge, and gives that edge a copy of
//! the caller's `GeomFill_CurveAndTrihedron(GeomFill_Darboux)` law restricted to
//! the pcurve. When no face carries a pcurve the law gives up: `hasresult`
//! turns false and `myLaws` is nullified.
//!
//! Translation notes:
//! - `BRep_Tool::CurveOnSurface(E, F, First, Last)` is
//!   [`GeometryRegistry::edge_pcurve`] keyed by the face's shape key; the
//!   `First` / `Last` out-parameters are the edge's range
//!   ([`BRepTool::edge_parameters`]), which is where OCCT reads them too.
//! - `Geom2d_TrimmedCurve(C, First, Last)` + `Reverse()` (the `REVERSED` edge
//!   arm, `cxx:89-94`) is [`GeomTrimmedCurve2d`].
//! - `BRepAdaptor_Surface(F)` (`cxx:76`) is only handed to
//!   `Adaptor3d_CurveOnSurface`, so the port passes [`BRepTool::face_surface`]
//!   directly; `Adaptor3d_CurveOnSurface(AC2d, AS)` is
//!   [`occt_geom::curve_on_surface::CurveOnSurface`].
//! - `GeomFill_Darboux` reads the pcurve and the surface out of the
//!   `Adaptor3d_CurveOnSurface` it is given through `SetCurve`
//!   (`GeomFill_Darboux.cxx:369-372`); the port's
//!   [`occt_geom::geom_fill::Darboux::new`] takes those two parts directly, so
//!   `new GeomFill_Darboux()` + `Law->SetCurve(AC)` is one `Darboux::new`.

use std::sync::Arc;

use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::geom_fill::{CurveAndTrihedron, Darboux, LocationLaw as GeomFillLocationLaw};

use crate::abs::Orientation;
use crate::abs::ShapeType;
use crate::brep_fill_location_law::LocationLawBase;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::shape::{TopoShape, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::explorer_shapes;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

/// `BRepFill_EdgeOnSurfLaw` (`BRepFill_EdgeOnSurfLaw.hxx:30-35`).
pub struct EdgeOnSurfLaw {
    /// The inherited `BRepFill_LocationLaw` state.
    base: LocationLawBase,
    /// `hasresult` (`hxx:33`).
    has_result: bool,
}

impl EdgeOnSurfLaw {
    /// `BRepFill_EdgeOnSurfLaw(Path, Surf)` (`cxx:40-101`).
    pub fn new(path: &Wire, surf: &TopoShape) -> Self {
        let mut base = LocationLawBase::new();
        base.init(path);
        let mut has_result = true;

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

                // `cxx:69-78`: the first face of `Surf` carrying a pcurve of
                // the edge.
                let mut found = false;
                let mut pcurve = None;
                let mut surface = None;
                for face_shape in explorer_shapes(surf, ShapeType::Face) {
                    let Some(face) = crate::shape::Face::wrap(face_shape) else {
                        continue;
                    };
                    let key = GeometryRegistry::shape_key(&face.0);
                    let c = GeometryRegistry::global().edge_pcurve(&e.0, key);
                    if let Some(c) = c {
                        found = true;
                        surface = BRepTool::face_surface(&face);
                        pcurve = Some(c);
                        break;
                    }
                }
                if !found {
                    // `cxx:79-83`: "Impossible to construct the law."
                    has_result = false;
                    base.nullify_laws();
                    return Self { base, has_result };
                }
                let Some(pcurve) = pcurve else {
                    has_result = false;
                    base.nullify_laws();
                    return Self { base, has_result };
                };
                let Some(surface) = surface else {
                    // `BRepAdaptor_Surface` over a face without a surface: OCCT
                    // would hand a null adaptor to the curve-on-surface.
                    has_result = false;
                    base.nullify_laws();
                    return Self { base, has_result };
                };

                // `cxx:85-94`: `First` / `Last` are the edge's own range.
                let (mut first, mut last) = BRepTool::edge_parameters(&e);
                let mut curve2d = pcurve;
                if e.0.orientation() == Orientation::Reversed {
                    // `CBis = new Geom2d_TrimmedCurve(C, First, Last);
                    //  CBis->Reverse();` then the reversed curve's range.
                    let mut cbis = Geom2dTrimmedCurve::new(Arc::clone(&curve2d), first, last);
                    cbis.reverse();
                    let cbis: Arc<dyn occt_geom2d::curve::Curve2d> = Arc::new(cbis);
                    curve2d = cbis;
                    first = curve2d.first_parameter();
                    last = curve2d.last_parameter();
                }

                // `AC2d = new Geom2dAdaptor_Curve(C, First, Last)` (`cxx:96`):
                // the adaptor evaluates `C` over `[First, Last]`.
                let ac2d: Arc<dyn occt_geom2d::curve::Curve2d> =
                    Arc::new(Geom2dTrimmedCurve::new(Arc::clone(&curve2d), first, last));

                // `AC = new Adaptor3d_CurveOnSurface(AC2d, AS)` (`cxx:97`).
                let ac: Arc<dyn Curve> =
                    Arc::new(CurveOnSurface::new(Arc::clone(&ac2d), Arc::clone(&surface)));

                // `TLaw = new GeomFill_Darboux(); Law = new
                // GeomFill_CurveAndTrihedron(TLaw)` (`cxx:58-59`), then
                // `myLaws->SetValue(ipath, Law->Copy())` +
                // `->SetCurve(AC)` (`cxx:98-99`).
                let darboux: Arc<dyn occt_geom::geom_fill::TrihedronLaw> =
                    Arc::new(Darboux::new(Arc::clone(&ac2d), surface));
                let law: Arc<dyn GeomFillLocationLaw> = Arc::new(CurveAndTrihedron::new(darboux));
                let _ = law.set_curve(ac);
                base.set_law(ipath, law);
            }
            wexp.next();
        }

        Self { base, has_result }
    }

    /// `HasResult()` (`cxx:103-106`).
    pub fn has_result(&self) -> bool {
        self.has_result
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

impl crate::brep_fill_location_law::LocationLaw for EdgeOnSurfLaw {
    fn base(&self) -> &LocationLawBase {
        &self.base
    }

    fn base_mut(&mut self) -> &mut LocationLawBase {
        &mut self.base
    }
}
