//! `BRepOffset_Tool::PipeInter` (`BRepOffset_Tool.cxx:747-804`): intersection
//! of two tube faces, i.e. of the surfaces that carry them, with the resulting
//! 3D curves turned into edges carrying a p-curve on each face.
//!
//! Called from `BRepOffset_Inter3d::FaceInter` (`:214`, `:246`) and
//! `BRepOffset_Inter3d::ConnexIntByArc` (`:379`) whenever both faces come from
//! edges and the edges share a vertex without a sphere.

use std::sync::Arc;

use occt_core::precision::CONFUSION;
use occt_geom2d::curve::Curve2d;

use crate::abs::State;
use crate::brep_lib_make_edge::{add_endpoint_vertices, init_3d};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::geom_int::IntSS;
use crate::shape::{Face, TopoShape};
use crate::tgeometry::GeometryRegistry;

use super::curve_checks::to_small;
use super::pcurve_builder::{build_pcurves, put_in_bounds};

impl crate::brep_offset::BRepOffsetTool {
    /// `PipeInter(F1, F2, L1, L2, Side)` (`cxx:747-804`). Both lists are
    /// cleared first, as in the reference (`cxx:756-757`), and only filled when
    /// `GeomInt_IntSS` reports `IsDone()`.
    pub(crate) fn pipe_inter(
        f1: &Face,
        f2: &Face,
        l1: &mut Vec<TopoShape>,
        l2: &mut Vec<TopoShape>,
        side: State,
    ) {
        l1.clear();
        l2.clear();

        // `cxx:759-760`: `BRep_Tool::Surface(F1)` / `Surface(F2)`. OCCT
        // dereferences both handles unguarded; the port returns.
        let (Some(s1), Some(s2)) = (BRepTool::face_surface(f1), BRepTool::face_surface(f2)) else {
            return;
        };

        // `cxx:762`: `GeomInt_IntSS Inter(S1, S2, Precision::Confusion(), true,
        // true, true)`.
        let mut inter = IntSS::new();
        inter.perform(
            s1,
            s2,
            CONFUSION,
            /* approx */ true,
            /* approx_s1 */ true,
            /* approx_s2 */ true,
        );

        // `cxx:764`.
        if !inter.is_done() {
            return;
        }

        let reg = GeometryRegistry::global();
        let f1_key = GeometryRegistry::shape_key(&f1.0);
        let f2_key = GeometryRegistry::shape_key(&f2.0);

        // `cxx:766-802`.
        for i in 1..=inter.nb_lines() {
            // `cxx:768`: `CI = Inter.Line(i)`.
            let Some(line) = inter.line(i) else {
                continue;
            };
            let ci = line.curve.clone();

            // `cxx:769-772`.
            if to_small(ci.as_ref()) {
                continue;
            }

            // `cxx:773`: `TopoDS_Edge E = BRepLib_MakeEdge(CI)`. A failed
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

            // `cxx:774-783`: the p-curve on `F1` comes from the intersection
            // when `GeomInt_IntSS` produced one, otherwise it is built.
            let pc1: Option<Arc<dyn Curve2d>> =
                inter.line_on_s1(i).map(|c| put_in_bounds(f1, &e, c.clone()));
            match pc1 {
                Some(c2) => reg.set_edge_pcurve(&e.0, f1_key, c2),
                None => {
                    let _ = build_pcurves(&e, f1);
                }
            }

            // `cxx:784-793`: the same for `F2`.
            let pc2: Option<Arc<dyn Curve2d>> =
                inter.line_on_s2(i).map(|c| put_in_bounds(f2, &e, c.clone()));
            match pc2 {
                Some(c2) => reg.set_edge_pcurve(&e.0, f2_key, c2),
                None => {
                    let _ = build_pcurves(&e, f2);
                }
            }

            // `cxx:794-801`: `OrientSection(E, F1, F2, O1, O2)`, both reversed
            // on the `TopAbs_OUT` side. `OrientSection` answers `None` where
            // the reference dereferences a null handle; the port drops the
            // line then.
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
