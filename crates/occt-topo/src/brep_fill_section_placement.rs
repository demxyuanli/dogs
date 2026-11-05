//! `BRepFill_SectionPlacement` — where a sweep section sits on the spine.
//!
//! Source: `ModelingAlgorithms/TKBool/BRepFill/BRepFill_SectionPlacement.cxx:44-281`,
//! `.hxx:30-90` (OCCT 8.0.0).
//!
//! The class turns the caller's section shape into one `Geom_Geometry` (a
//! curve, or a `Geom_CartesianPoint` for a punctual section), asks
//! `GeomFill_SectionPlacement` for the best parameter on the *whole* path, maps
//! that global parameter onto the `BRepFill_LocationLaw` interval it falls in,
//! and then re-runs the placement on that single law — once for the interval's
//! law and, when the parameter sits on a knot, once more for the neighbouring
//! law; the answer with the larger `Angle()` wins (`cxx:248-268`).
//!
//! Translation notes:
//! - The OCCT object keeps `myLaw` / `mySection`; both are only read inside
//!   `Perform` (`AbscissaOnPath` reads the law through the cached `myIndex` /
//!   `myParam`), so the port stores the four results and computes
//!   `AbscissaOnPath` eagerly — the law is immutable after `Init`, so the value
//!   is the same.
//! - `BRepAdaptor_CompCurve(myLaw->Wire())` (`cxx:170`) is
//!   [`CompCurveAdaptor::of_wire`] with `theIsAC = false`, i.e. the global
//!   parameter is the edge-rank space `[i-1, i]`, the same space
//!   `BRepFill_LocationLaw::Init` gives its laws.
//! - `GeomFill_SectionPlacement` is
//!   [`occt_geom::geom_fill::section_placement::SectionPlacement`];
//!   `Perform(Adaptor3d_Curve, Tol)` is `perform_on`, `Perform(Tol)` is
//!   `perform`, `Perform(Param, Tol)` is `perform_at`.
//! - `GeomConvert_CompCurveToBSplineCurve::Add(C, Tol)` defaults
//!   (`hxx:56-60`) are passed explicitly: `After = false`, `WithRatio = true`,
//!   `MinM = 0`. `Add` answering `false` is retried with
//!   `max(Confusion, vertexTolerance)` (`cxx:154-157`).
//! - `BRep_Tool::Curve(E, f, l)` (`BRep_Tool.cxx:205-219`) **applies** the edge
//!   location, so the port uses [`BRepTool::edge_curve_world`]; the range is
//!   unchanged by the transform.
//! - `TopoDS::Vertex(Current())` / `TopoDS::Edge(Current())` down-casts and the
//!   `Standard_ConstructionError("Interval is not found")` (`cxx:206`) become
//!   skipped shapes and `Err`, respectively.

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::gp::GpTrsf;
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::convert_bspl::CompCurveToBSplineCurve;
use occt_geom::curve::Curve;
use occt_geom::geom_fill::section_placement::{
    SectionGeometry, SectionPlacement as GeomFillSectionPlacement,
};
use occt_geom::trimmed::GeomTrimmedCurve;

use crate::abs::{Orientation, ShapeType};
use crate::brep_adaptor_comp_curve::CompCurveAdaptor;
use crate::brep_fill_location_law::LocationLawBase;
use crate::brep_tool::BRepTool;
use crate::iterator::ShapeExplorer;
use crate::shape::{Edge, TopoShape, Vertex};
use crate::topo_tools::{edge_vertices, is_same};

const INTERVAL_NOT_FOUND: &str = "BRepFill_SectionPlacement: Interval is not found";
const NO_LOCATION_LAW: &str = "BRepFill_SectionPlacement: null location law";

/// `SearchParam(Law, Ind, TheV)` (`cxx:44-61`): the law parameter of the
/// vertex `TheV` on law `Ind`.
fn search_param(law: &LocationLawBase, ind: usize, the_v: &TopoShape) -> Result<f64, String> {
    let e = law
        .edge(ind)
        .ok_or_else(|| format!("BRepFill_SectionPlacement: no edge {ind}"))?;
    // `t = BRep_Tool::Parameter(TheV, E)`.
    let mut t = BRepTool::parameter_of_vertex(&e.0, the_v);
    if e.0.orientation() == Orientation::Reversed {
        // `cxx:49-59`: `C` is read for `f` / `l` only.
        let (f, l) = BRepTool::edge_parameters(&e);
        let law_curve = law
            .law(ind)
            .and_then(|law| law.get_curve())
            .ok_or_else(|| NO_LOCATION_LAW.to_string())?;
        let lf = law_curve.first_parameter();
        let ll = law_curve.last_parameter();
        t = ll - (t - f) * (ll - lf) / (l - f);
    }
    Ok(t)
}

/// `BRepFill_SectionPlacement` (`BRepFill_SectionPlacement.hxx:30-90`).
#[derive(Debug, Clone)]
pub struct SectionPlacement {
    /// `myTrsf` (`hxx:85`).
    trsf: GpTrsf,
    /// `myIndex` (`hxx:86`).
    index: usize,
    /// `myParam` (`hxx:87`).
    param: f64,
    /// `myLaw->Abscissa(myIndex, myParam)` (`AbscissaOnPath`, `cxx:277-280`).
    abscissa: f64,
}

impl SectionPlacement {
    /// `BRepFill_SectionPlacement(Law, Section, WithContact, WithCorrection)`
    /// (`cxx:63-74`): the "no vertex" constructor, which forwards a null
    /// vertex.
    pub fn new(
        law: &LocationLawBase,
        section: &TopoShape,
        with_contact: bool,
        with_correction: bool,
    ) -> Result<Self, String> {
        Self::perform(law, section, None, with_contact, with_correction)
    }

    /// `BRepFill_SectionPlacement(Law, Section, Vertex, WithContact,
    /// WithCorrection)` (`cxx:76-84`).
    pub fn with_vertex(
        law: &LocationLawBase,
        section: &TopoShape,
        vertex: &TopoShape,
        with_contact: bool,
        with_correction: bool,
    ) -> Result<Self, String> {
        Self::perform(law, section, Some(vertex), with_contact, with_correction)
    }

    /// `Perform(WithContact, WithCorrection, Vertex)` (`cxx:86-270`).
    fn perform(
        law: &LocationLawBase,
        section: &TopoShape,
        vertex: Option<&TopoShape>,
        with_contact: bool,
        with_correction: bool,
    ) -> Result<Self, String> {
        // `cxx:93-166`: the section as one geometric object.
        let a_section = section_geometry(section)?;

        // `cxx:169-170`: `GeomFill_SectionPlacement(myLaw->Law(1), aSection)`.
        let law1 = law.law(1).ok_or_else(|| NO_LOCATION_LAW.to_string())?;
        let mut placement = GeomFillSectionPlacement::new(law1, a_section);

        // `cxx:170-171`: the whole path, then the global parameter.
        let path = CompCurveAdaptor::of_wire(law.wire(), false);
        placement.perform_on(&path, CONFUSION);
        let section_param = placement.parameter_on_path();
        let param_confusion = PCONFUSION;

        // `cxx:175-207`: the interval of the location-law parameterisation.
        let mut law_index1 = 0usize;
        let mut law_index2 = 0usize;
        let mut is_interval_found = false;
        for law_index in 1..=law.nb_law() {
            let curr_knot = (law_index - 1) as f64;
            let next_knot = law_index as f64;
            is_interval_found = (curr_knot - param_confusion <= section_param)
                && (next_knot + param_confusion >= section_param);
            if !is_interval_found {
                continue;
            }
            law_index1 = law_index;
            if (section_param - curr_knot).abs() < param_confusion && law_index > 1 {
                law_index2 = law_index - 1;
            } else if (section_param - next_knot).abs() < param_confusion && law_index < law.nb_law() {
                law_index2 = law_index + 1;
            }
        }
        if !is_interval_found {
            return Err(INTERVAL_NOT_FOUND.to_string());
        }

        // `cxx:210-227`: is the vertex an end of a law edge?
        let mut is_vertex_on_law = false;
        if let Some(the_v) = vertex {
            for current_law_index in 1..=law.nb_law() {
                let Some(an_edge) = law.edge(current_law_index) else {
                    continue;
                };
                let (v1, v2) = edge_vertices(&an_edge);
                if vertex_is(&v1, the_v) || vertex_is(&v2, the_v) {
                    is_vertex_on_law = true;
                    law_index1 = current_law_index;
                    law_index2 = 0;
                    break;
                }
            }
        }

        // `cxx:232-246`: positioning on the localized edge (or two edges).
        placement.set_location(law.law(law_index1).ok_or_else(|| NO_LOCATION_LAW.to_string())?);
        if !is_vertex_on_law {
            placement.perform(CONFUSION);
        } else {
            let the_v = vertex.expect("vertex flag without a vertex");
            placement.perform_at(search_param(law, law_index1, the_v)?, CONFUSION);
        }

        let mut trsf = placement.transformation(with_contact, with_correction)?;
        let mut index = law_index1;
        let mut param = placement.parameter_on_path();
        let angle = placement.angle();

        // `cxx:248-268`: the neighbouring law of a knot.
        if law_index2 != 0 {
            placement.set_location(
                law.law(law_index2)
                    .ok_or_else(|| NO_LOCATION_LAW.to_string())?,
            );
            if !is_vertex_on_law {
                placement.perform(CONFUSION);
            } else {
                let mut a_vertex = vertex.expect("vertex flag without a vertex").clone();
                if law_index1 == law_index2 {
                    // `cxx:257-260`: the same law read from the other side.
                    a_vertex.reverse();
                }
                placement.perform_at(search_param(law, law_index2, &a_vertex)?, CONFUSION);
            }
            if placement.angle() > angle {
                trsf = placement.transformation(with_contact, with_correction)?;
                index = law_index2;
                param = placement.parameter_on_path();
            }
        }

        // `AbscissaOnPath` (`cxx:277-280`) on the cached index / parameter.
        let abscissa = law.abscissa(index, param)?;

        Ok(Self {
            trsf,
            index,
            param,
            abscissa,
        })
    }

    /// `Transformation()` (`cxx:272-275`).
    pub fn transformation(&self) -> &GpTrsf {
        &self.trsf
    }

    /// `AbscissaOnPath()` (`cxx:277-280`).
    pub fn abscissa_on_path(&self) -> f64 {
        self.abscissa
    }

    /// `myIndex` (`hxx:86`), 1-based law index of the winning placement.
    pub fn index(&self) -> usize {
        self.index
    }

    /// `myParam` (`hxx:87`), the parameter on that law.
    pub fn parameter_on_path(&self) -> f64 {
        self.param
    }
}

/// `V1.IsSame(aVertex)` on an optional `TopExp::Vertices` end.
fn vertex_is(v: &Option<Vertex>, the_v: &TopoShape) -> bool {
    matches!(v, Some(vertex) if is_same(&vertex.0, the_v))
}

/// `cxx:93-166`: the first non-degenerated edge with a 3D curve, trimmed to its
/// range and concatenated with the remaining edges; a punctual section becomes
/// a `Geom_CartesianPoint`.
fn section_geometry(section: &TopoShape) -> Result<SectionGeometry, String> {
    let edges: Vec<Edge> = ShapeExplorer::new(vec![section.clone()], ShapeType::Edge)
        .filter_map(|shape| Edge::wrap(shape))
        .collect();

    // `cxx:96-111`: the first valid curve.
    let mut first_index: Option<usize> = None;
    let mut found: Option<(Arc<dyn Curve>, f64, f64)> = None;
    for (index, edge) in edges.iter().enumerate() {
        if BRepTool::is_degenerated(edge) {
            continue;
        }
        if let Some(curve) = BRepTool::edge_curve_world(edge) {
            let (first, last) = BRepTool::edge_parameters(edge);
            found = Some((curve, first, last));
            first_index = Some(index);
            break;
        }
    }

    // `cxx:113-120`: no edge found, the section is a vertex.
    let Some((curve, first, last)) = found else {
        let Some(the_v) = ShapeExplorer::new(vec![section.clone()], ShapeType::Vertex).next() else {
            return Err("BRepFill_SectionPlacement: empty section".to_string());
        };
        // `cxx:117-119`: `BRep_Tool::Pnt` is the located point.
        return Ok(SectionGeometry::Point(BRepTool::vertex_point_world(&Vertex(the_v))));
    };

    // `cxx:121-160`: trim the first edge, concatenate the rest.
    let first_index = first_index.expect("curve without an edge index");
    let mut a_curve: Arc<dyn Curve> = Arc::new(GeomTrimmedCurve::new(curve, first, last));

    let remaining = &edges[first_index + 1..];
    let has_more = remaining
        .iter()
        .any(|edge| !BRepTool::is_degenerated(edge) && BRepTool::edge_curve_world(edge).is_some());
    if has_more {
        let a_precision_tolerance = CONFUSION;
        let mut converter =
            CompCurveToBSplineCurve::from_curve(a_curve.as_ref(), ParameterisationType::TgtThetaOver2)
                .map_err(|e| format!("GeomConvert_CompCurveToBSplineCurve: {e:?}"))?;
        for edge in remaining {
            // `cxx:134-139`: avoid null, degenerated edges.
            if BRepTool::is_degenerated(edge) {
                continue;
            }
            let Some(c) = BRepTool::edge_curve_world(edge) else {
                continue;
            };
            let (f, l) = BRepTool::edge_parameters(edge);

            // `cxx:147-151`.
            let (v_first, v_last) = edge_vertices(edge);
            let tol_first = v_first.as_ref().map(BRepTool::vertex_tolerance).unwrap_or(0.0);
            let tol_last = v_last.as_ref().map(BRepTool::vertex_tolerance).unwrap_or(0.0);
            let vertex_tolerance = tol_first.max(tol_last);

            // `cxx:153-158`.
            let trimmed = GeomTrimmedCurve::new(c, f, l);
            let added = converter
                .add(
                    &trimmed,
                    a_precision_tolerance.min(vertex_tolerance),
                    false,
                    true,
                    0,
                )
                .map_err(|e| format!("GeomConvert_CompCurveToBSplineCurve::Add: {e:?}"))?;
            if !added {
                converter
                    .add(
                        &trimmed,
                        a_precision_tolerance.max(vertex_tolerance),
                        false,
                        true,
                        0,
                    )
                    .map_err(|e| format!("GeomConvert_CompCurveToBSplineCurve::Add: {e:?}"))?;
            }
        }
        // `cxx:159`: `aCurve = aBSplineConverter.BSplineCurve()`.
        a_curve = match converter.into_curve() {
            Some(bspline) => Arc::new(bspline),
            None => a_curve,
        };
    }

    Ok(SectionGeometry::Curve(a_curve))
}
