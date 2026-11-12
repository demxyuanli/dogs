//! `BRep_Builder::Transfert` — copy curve representations (and vertex
//! parameters) from one edge onto another.
//!
//! Source: `BRep_Builder.cxx:1159-1198` (edge / edge), `:1457-1465` (edge +
//! vertex / edge + vertex).

use std::sync::Arc;

use occt_core::toploc::TopLocLocation;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::brep_builder_update_vertex::update_vertex_on_edge;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::tshape_continuity::write_continuity_on_surfaces;

/// Snapshot of one `BRep_GCurve` CurveOnSurface / CurveOnClosedSurface slot.
struct PcurveSlot {
    surface: Arc<dyn Surface>,
    curves: Vec<Arc<dyn Curve2d>>,
    location: TopLocLocation,
    /// `BRep_CurveOnClosedSurface::Continuity` when `curves.len() >= 2`.
    closed_continuity: Option<u8>,
}

/// Snapshot of one `BRep_CurveOn2Surfaces` regularity.
struct RegularitySlot {
    surface1: Option<Arc<dyn Surface>>,
    surface2: Option<Arc<dyn Surface>>,
    location1: TopLocLocation,
    location2: TopLocLocation,
    continuity: u8,
}

fn snapshot_ein(ein: &TopoShape) -> (f64, Vec<PcurveSlot>, Vec<RegularitySlot>) {
    let reg = GeometryRegistry::global();
    let tol = {
        let ts = ein.tshape.read().expect("poisoned TShape lock");
        ts.edge_core()
            .map(|c| c.tolerance)
            .unwrap_or(occt_core::precision::REAL_EPSILON)
    };

    let (pcurves, regularities) = {
        let ts = ein.tshape.read().expect("poisoned TShape lock");
        let pcs = ts.edge_pcurves().cloned();
        let regs = ts
            .edge_core()
            .map(|c| c.regularities.clone())
            .unwrap_or_default();
        (pcs, regs)
    };

    let mut slots = Vec::new();
    if let Some(g) = pcurves {
        for (&key, curves) in g.curves.iter() {
            let Some(surface) = reg.surface_at_ptr(key) else {
                continue;
            };
            let closed_continuity = if curves.len() >= 2 {
                Some(g.continuities.get(&key).copied().unwrap_or(0))
            } else {
                None
            };
            slots.push(PcurveSlot {
                surface,
                curves: curves.clone(),
                location: g.location(key),
                closed_continuity,
            });
        }
    }

    let regs: Vec<RegularitySlot> = regularities
        .into_iter()
        .map(|r| RegularitySlot {
            surface1: r.surface1,
            surface2: r.surface2,
            location1: r.location1,
            location2: r.location2,
            continuity: r.continuity,
        })
        .collect();

    (tol, slots, regs)
}

/// `BRep_Builder::Transfert(Ein, Eout)` (`BRep_Builder.cxx:1159-1198`).
///
/// Walks Ein's curve representations: every CurveOnSurface (including a
/// CurveOnClosedSurface, which inherits `IsCurveOnSurface()` and therefore
/// takes the single-pcurve arm — the `else if IsCurveOnClosedSurface` branch
/// at `cxx:1177` is unreachable in OCCT) is written onto Eout through
/// `UpdateEdge`; every `IsRegularity()` representation is written through
/// `Continuity`.
pub fn transfert_edge(ein: &Edge, eout: &Edge) {
    let ein_loc = ein.0.location().clone();
    let (tol, slots, regs) = snapshot_ein(&ein.0);
    let reg = GeometryRegistry::global();

    for slot in slots {
        // `cxx:1172-1175`: `IsCurveOnSurface()` — first pcurve only.
        let Some(pc0) = slot.curves.first().cloned() else {
            continue;
        };
        let abs_loc = ein_loc.multiplied(&slot.location);
        reg.update_edge_pcurve_on_surface_at(
            &eout.0,
            &slot.surface,
            Some(pc0),
            abs_loc.clone(),
            tol,
        );
        // `cxx:1187-1195`: CurveOnClosedSurface is also `IsRegularity()`.
        if let Some(cont) = slot.closed_continuity {
            write_continuity_on_surfaces(
                eout,
                Some(slot.surface.clone()),
                Some(slot.surface),
                abs_loc.clone(),
                abs_loc,
                cont,
            );
        }
    }

    for r in regs {
        // Stored locations are already relative to Ein; restore absolute
        // `Ein.Location() * CR->Location()` (`cxx:1192-1193`).
        let abs1 = ein_loc.multiplied(&r.location1);
        let abs2 = ein_loc.multiplied(&r.location2);
        write_continuity_on_surfaces(
            eout,
            r.surface1,
            r.surface2,
            abs1,
            abs2,
            r.continuity,
        );
    }
}

/// `BRep_Builder::Transfert(Ein, Eout, Vin, Vout)` (`BRep_Builder.cxx:1457-1465`):
/// copy Vin's parameter on Ein onto Vout on Eout, raising Vout's tolerance.
pub fn transfert_edge_vertex(ein: &Edge, eout: &Edge, vin: &Vertex, vout: &Vertex) {
    let tol = BRepTool::vertex_tolerance(vin);
    let par = BRepTool::parameter_of_vertex(&ein.0, &vin.0);
    let _ = update_vertex_on_edge(&vout.0, par, &eout.0, tol);
}

/// Convenience over [`transfert_edge`] for bare `TopoShape` edges.
pub fn transfert_edge_shapes(ein: &TopoShape, eout: &TopoShape) {
    if let (Some(a), Some(b)) = (Edge::wrap(ein.clone()), Edge::wrap(eout.clone())) {
        transfert_edge(&a, &b);
    }
}
