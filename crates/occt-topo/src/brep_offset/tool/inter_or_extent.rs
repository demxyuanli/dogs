//! `BRepOffset_Tool::InterOrExtent` (`BRepOffset_Tool.cxx:2010-2068`):
//! intersection of two faces' surfaces with a `Geom_RectangularTrimmedSurface`
//! wrapping a plane unwrapped to its basis, then one edge per resulting line
//! carrying a p-curve built on each face.
//!
//! Difference to [`super::pipe_inter`] (which is the same walk): OCCT passes the
//! `GeomInt_IntSS` defaults here (`Approx = true`, `ApproxS1 = ApproxS2 =
//! false`, `GeomInt_IntSS.hxx:46-51`) and calls `BuildPCurves` on both faces
//! instead of writing the intersector's own p-curves.
//!
//! The method has no caller in OCCT 8.0.0 (`BRepOffset_Tool.hxx:117` declares
//! it, nothing in `TKOffset` calls it); it is ported so the class surface
//! matches, like the reference does.

use std::sync::Arc;

use occt_core::precision::CONFUSION;
use occt_geom::Surface;

use crate::abs::State;
use crate::brep_lib_make_edge::{add_endpoint_vertices, init_3d};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::geom_int::IntSS;
use crate::shape::{Face, TopoShape};

use super::curve_checks::to_small;
use super::pcurve_builder::build_pcurves;

impl crate::brep_offset::BRepOffsetTool {
    /// `InterOrExtent(F1, F2, L1, L2, Side)` (`cxx:2010-2068`). Both lists are
    /// cleared first, as in the reference (`cxx:2019-2020`).
    pub fn inter_or_extent(
        f1: &Face,
        f2: &Face,
        l1: &mut Vec<TopoShape>,
        l2: &mut Vec<TopoShape>,
        side: State,
    ) {
        l1.clear();
        l2.clear();

        // `cxx:2021-2022`: `BRep_Tool::Surface(F1)` / `Surface(F2)`. OCCT
        // dereferences both handles unguarded; the port returns.
        let (Some(s1), Some(s2)) = (BRepTool::face_surface(f1), BRepTool::face_surface(f2)) else {
            return;
        };

        // `cxx:2024-2041`.
        let s1 = unwrap_trimmed_plane(s1);
        let s2 = unwrap_trimmed_plane(s2);

        // `cxx:2043`: the default `(Approx, ApproxS1, ApproxS2)` of
        // `GeomInt_IntSS.hxx:46-51`.
        let mut inter = IntSS::new();
        inter.perform(
            s1,
            s2,
            CONFUSION,
            /* approx */ true,
            /* approx_s1 */ false,
            /* approx_s2 */ false,
        );

        // `cxx:2045`.
        if !inter.is_done() {
            return;
        }

        // `cxx:2047-2066`.
        for i in 1..=inter.nb_lines() {
            // `cxx:2049`.
            let Some(line) = inter.line(i) else {
                continue;
            };
            let ci = line.curve.clone();

            // `cxx:2051-2054`.
            if to_small(ci.as_ref()) {
                continue;
            }

            // `cxx:2055`: `TopoDS_Edge E = BRepLib_MakeEdge(CI)`. A failed
            // `BRepLib_MakeEdge` leaves a null edge in OCCT, which the two
            // `L.Append` calls at the end would store; the port drops the line.
            let (first, last) = (ci.first_parameter(), ci.last_parameter());
            let Some(build) = init_3d(&ci, first, last) else {
                continue;
            };
            let builder = TopoBuilder::new();
            let mut e = builder.make_edge(ci.clone(), build.first, build.last);
            add_endpoint_vertices(&builder, &mut e, &build);
            if build.closed {
                e.0.set_closed(true);
            }

            // `cxx:2056-2057`.
            let _ = build_pcurves(&e, f1);
            let _ = build_pcurves(&e, f2);

            // `cxx:2058-2063`.
            let Some((mut o1, mut o2)) = Self::orient_section(&e, f1, f2) else {
                continue;
            };
            if side == State::Out {
                o1 = o1.reversed();
                o2 = o2.reversed();
            }
            l1.push(e.0.oriented(o1));
            l2.push(e.0.oriented(o2));
        }
    }
}

/// The `S->DynamicType() == Geom_RectangularTrimmedSurface &&
/// RTS->BasisSurface()->DynamicType() == Geom_Plane` unwrap of
/// `cxx:2024-2041`: a plane-trimmed surface is replaced by its basis plane, any
/// other surface is kept as is.
fn unwrap_trimmed_plane(s: Arc<dyn Surface>) -> Arc<dyn Surface> {
    match s.rectangular_trimmed_basis() {
        Some(basis) if basis.gp_pln().is_some() => basis,
        _ => s,
    }
}
