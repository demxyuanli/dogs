//! `BRepFill_ShapeLaw` — the section law of a single wire (or vertex).
//!
//! Source: `BRepFill_ShapeLaw.cxx:48-486`, `.hxx:38-86`.
//!
//! One `GeomFill_SectionLaw` per non-degenerated edge of the wire: a
//! `GeomFill_UniformSection` when the law is constant, a
//! `GeomFill_EvolvedSection` when a `Law_Function` scales the section.
//!
//! Translation notes:
//! - The inherited `BRepFill_SectionLaw` data (including the edge -> law index
//!   map) lives in [`SectionLawBase`]; `myShape`, `myEdges`, `TheLaw` and the
//!   protected `vertex` flag are this class's own fields (`hxx:75-81`).
//! - `myEdges` is a `HArray1<TopoDS_Shape>` (`hxx:79`) — a punctual law stores
//!   the vertex itself in slot 1 (`cxx:60`) — so the port keeps
//!   `Vec<Option<TopoShape>>` and down-casts in [`ShapeLaw::edge`].
//! - `BRep_Tool::Curve(E, First, Last)` is [`BRepTool::edge_curve_world`];
//!   `BRep_Tool::Parameter(V, E)` is [`BRepTool::parameter_of_vertex`];
//!   `BRep_Tool::IsClosed(E)` is `TopoDS_Shape::Closed()`, i.e. the TShape
//!   flag (`e.0.is_closed()`), *not* a parameter comparison.
//! - `TopExp::FirstVertex` / `LastVertex` use OCCT's default `CumOri = false`
//!   (`TopExp.hxx:116`, `:121`), i.e. the vertices stored `FORWARD` /
//!   `REVERSED` inside the edge.
//! - `Geom_TrimmedCurve(C, First, Last)` is [`GeomTrimmedCurve`]: the port
//!   reparameterises the window onto `[0, 1]`, so a trimmed section curve is
//!   evaluated on `[0, 1]` and `GeomFill_UniformSection` receives that range.
//! - `GeomFill_UniformSection(C)` / `GeomFill_EvolvedSection(C, L)` use OCCT's
//!   default bounds `0.` / `1.` (`GeomFill_UniformSection.hxx:41-43`), which is
//!   what the port passes explicitly.
//! - `Geom_BSplineCurve(Poles, Weights, Knots, Mults, Degree, Periodic)` is
//!   assembled through [`GeomBSplineCurve::from_parts`] after the same
//!   `updateKnots()` expansion
//!   [`GeomBSplineCurve::from_poles_knots_mults_periodic`] applies.
//! - `BRepLProp::Continuity` throws `Standard_Failure` when the two points are
//!   farther apart than the linear tolerance (`BRepLProp.cxx:67-69`); the
//!   `SectionLaw` trait has no error channel, so the port answers
//!   `GeomAbs_C0` — the same answer OCCT gives for an open law (`cxx:426-428`).

use std::sync::Arc;

use occt_core::gp::{GpDir, GpPnt, GpTrsf};
use occt_core::convert::ParameterisationType;
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::convert_bspl::CompCurveToBSplineCurve;
use occt_geom::curve::Curve;
use occt_geom::geom_fill::{
    EvolvedSection, SectionLaw as GeomFillSectionLaw, UniformSection,
};
use occt_geom::law::LawFunction;
use occt_geom::line::GeomLine;
use occt_geom::trimmed::GeomTrimmedCurve;

use crate::abs::Orientation;
use crate::brep_fill_section_law::{SectionLaw, SectionLawBase};
use crate::brep_lprop;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::shape_ops::transformed_copy;
use crate::topo_tools::{common_vertex, edge_vertices, is_same};

/// `BRepFill_ShapeLaw` (`BRepFill_ShapeLaw.hxx:38-86`).
pub struct ShapeLaw {
    /// The inherited `BRepFill_SectionLaw` state.
    base: SectionLawBase,
    /// `vertex` (`hxx:75`): the input shape is a vertex.
    vertex: bool,
    /// `myShape` (`hxx:78`).
    shape: TopoShape,
    /// `myEdges` (`hxx:79`), 1-based.
    edges: Vec<Option<TopoShape>>,
    /// `TheLaw` (`hxx:81`): null for a constant law.
    the_law: Option<Arc<dyn LawFunction>>,
}

impl ShapeLaw {
    /// `BRepFill_ShapeLaw(const TopoDS_Vertex&, Build = true)` (`cxx:52-75`):
    /// a constant law of a line segment starting at the vertex, as long as
    /// twice the vertex tolerance plus `Precision::PConfusion()`.
    pub fn with_vertex(v: &Vertex, build: bool) -> Result<Self, String> {
        let mut law = Self {
            base: SectionLawBase::new(),
            vertex: true,
            shape: v.0.clone(),
            edges: vec![Some(v.0.clone())],
            the_law: None,
        };
        // `uclosed = false; vclosed = true; // constant law` (`cxx:57-58`).
        law.base.set_u_closed(false);
        law.base.set_v_closed(true);

        if build {
            law.base.reset_laws(1);
            // `gp_Dir D(gp_Dir::D::X)` (`cxx:65`): the section line runs along X.
            let d = GpDir::new(1.0, 0.0, 0.0).map_err(|e| format!("gp_Dir: {e}"))?;
            let line = GeomLine::from_pnt_dir(BRepTool::vertex_point(v), d);
            let last = 2.0 * BRepTool::vertex_tolerance(v) + PCONFUSION;
            // `Geom_TrimmedCurve TC(L, 0, Last)` (`cxx:69`).
            let tc: Arc<dyn Curve> = Arc::new(GeomTrimmedCurve::new(Arc::new(line), 0.0, last));
            let section = UniformSection::new(tc, 0.0, 1.0)?;
            law.base.set_law_at(1, Some(Arc::new(section)));
        }
        law.base.set_done(true);
        Ok(law)
    }

    /// `BRepFill_ShapeLaw(const TopoDS_Wire&, Build = true)` (`cxx:79-85`).
    pub fn with_wire(w: &Wire, build: bool) -> Result<Self, String> {
        let mut law = Self::new_for_wire(w, None);
        law.init(build)?;
        law.base.set_done(true);
        Ok(law)
    }

    /// `BRepFill_ShapeLaw(const TopoDS_Wire&, const Handle(Law_Function)&, Build)`
    /// (`cxx:90-102`): an evolutive law.
    pub fn with_wire_and_law(
        w: &Wire,
        l: Arc<dyn LawFunction>,
        build: bool,
    ) -> Result<Self, String> {
        let mut law = Self::new_for_wire(w, Some(l));
        law.init(build)?;
        law.base.set_done(true);
        Ok(law)
    }

    /// The two wire constructors share this prologue (`cxx:79-102`).
    fn new_for_wire(w: &Wire, l: Option<Arc<dyn LawFunction>>) -> Self {
        Self {
            base: SectionLawBase::new(),
            vertex: false,
            shape: w.0.clone(),
            edges: Vec::new(),
            the_law: l,
        }
    }

    /// `Init(Build)` (`cxx:108-228`): one elementary law per non-degenerated
    /// edge that carries a curve, then the U-closedness check.
    fn init(&mut self, build: bool) -> Result<(), String> {
        self.base.set_v_closed(true);

        let wire = Wire::wrap(self.shape.clone())
            .ok_or_else(|| "Standard_TypeMismatch: BRepFill_ShapeLaw".to_string())?;

        // First pass (`cxx:115-126`): count the edges that have a curve.
        let mut nb_edge = 0usize;
        let mut wexp = WireExplorer::new();
        wexp.init(&wire, None);
        while wexp.more() {
            if let Some(e) = wexp.current_edge() {
                if !BRepTool::is_degenerated(&e) && BRepTool::edge_curve(&e).is_some() {
                    nb_edge += 1;
                }
            }
            wexp.next();
        }

        self.base.reset_laws(nb_edge);
        self.edges = vec![None; nb_edge];

        // Second pass (`cxx:133-184`): fill the tables.
        let mut ii = 1usize;
        let mut wexp = WireExplorer::new();
        wexp.init(&wire, None);
        while wexp.more() {
            if let Some(e) = wexp.current_edge() {
                if !BRepTool::is_degenerated(&e) && BRepTool::edge_curve(&e).is_some() {
                    let (mut first, mut last) = BRepTool::edge_parameters(&e);
                    if let Some(slot) = ii.checked_sub(1).and_then(|i| self.edges.get_mut(i)) {
                        *slot = Some(e.0.clone());
                    }
                    self.base.bind_index(&e.0, ii as i32);
                    if build {
                        let mut curve = BRepTool::edge_curve(&e)
                            .ok_or_else(|| "Standard_NullObject: BRep_Tool::Curve".to_string())?;
                        if e.0.orientation() == Orientation::Reversed {
                            // `C = C->Reversed()` with the parameters run through
                            // `ReversedParameter` (`cxx:147-157`).
                            let aux = curve.reversed_parameter(first);
                            first = curve.reversed_parameter(last);
                            last = aux;
                            curve = Arc::from(curve.reversed());
                        }

                        let mut is_closed = e.0.is_closed();
                        if is_closed && (curve.first_parameter() - first).abs() > PCONFUSION {
                            is_closed = false; // trimmed curve differs
                        }

                        if ii > 1 || !is_closed {
                            // `TC = new Geom_TrimmedCurve(C, First, Last)`
                            // (`cxx:167-168`).
                            curve = Arc::new(GeomTrimmedCurve::new(
                                Arc::clone(&curve),
                                first,
                                last,
                            ));
                        }

                        let section: Arc<dyn GeomFillSectionLaw> = match &self.the_law {
                            None => Arc::new(UniformSection::new(curve, 0.0, 1.0)?),
                            Some(l) => {
                                Arc::new(EvolvedSection::new(Arc::clone(&curve), l.clone())?)
                            }
                        };
                        self.base.set_law_at(ii, Some(section));
                    }
                    ii += 1;
                }
            }
            wexp.next();
        }

        // Is the law closed by U ? (`cxx:186-227`).
        let mut u_closed = wire.is_closed();
        if !u_closed {
            let edge1 = self.edge(self.edges.len());
            let edge2 = self.edge(1);
            if let (Some(edge1), Some(edge2)) = (edge1, edge2) {
                let v1 = if edge1.0.orientation() == Orientation::Reversed {
                    edge_vertices(&edge1).0
                } else {
                    edge_vertices(&edge1).1
                };
                let v2 = if edge2.0.orientation() == Orientation::Reversed {
                    edge_vertices(&edge2).1
                } else {
                    edge_vertices(&edge2).0
                };
                if let (Some(v1), Some(v2)) = (v1, v2) {
                    if is_same(&v1.0, &v2.0) {
                        u_closed = true;
                    } else {
                        let u1 = BRepTool::parameter_of_vertex(&edge1.0, &v1.0);
                        let u2 = BRepTool::parameter_of_vertex(&edge2.0, &v2.0);
                        let eps = BRepTool::vertex_tolerance(&v2)
                            + BRepTool::vertex_tolerance(&v1);
                        if let (Some(c1), Some(c2)) = (
                            BRepTool::edge_curve_world(&edge1),
                            BRepTool::edge_curve_world(&edge2),
                        ) {
                            u_closed = c1.d0(u1).is_equal_with_tolerance(&c2.d0(u2), eps);
                        }
                    }
                }
            }
        }
        self.base.set_u_closed(u_closed);
        Ok(())
    }

    /// `Edge(Index)` (`BRepFill_ShapeLaw.lxx:15-18`).
    pub fn edge(&self, index: usize) -> Option<Edge> {
        index
            .checked_sub(1)
            .and_then(|i| self.edges.get(i))
            .and_then(|slot| slot.clone())
            .and_then(Edge::wrap)
    }

    /// `IsVertex()` (`cxx:232-235`).
    pub fn is_vertex_shape(&self) -> bool {
        self.vertex
    }

    /// `IsConstant()` (`cxx:239-242`): the law has no `Law_Function`.
    pub fn is_constant_law(&self) -> bool {
        self.the_law.is_none()
    }

    /// `Vertex(Index, Param)` (`cxx:247-285`): the section vertex, scaled by
    /// the law when the section is evolutive.
    pub fn vertex_at(&self, index: usize, param: f64) -> Option<Vertex> {
        let nb = self.edges.len();
        let e = if index <= nb {
            self.edge(index)?
        } else if index == nb + 1 {
            self.edge(index - 1)?
        } else {
            return None;
        };
        let (first, last) = edge_vertices(&e);
        let v = if index <= nb {
            if e.0.orientation() == Orientation::Reversed {
                last
            } else {
                first
            }
        } else if e.0.orientation() == Orientation::Reversed {
            first
        } else {
            last
        };

        let mut v = v?;
        if let Some(law) = &self.the_law {
            // `T.SetScale(gp_Pnt(0,0,0), TheLaw->Value(Param))` (`cxx:278-282`).
            let mut t = GpTrsf::identity();
            if t.set_scale(&GpPnt::zero(), law.value(param)).is_ok() {
                if let Ok(moved) = transformed_copy(&v.0, &t) {
                    if let Some(scaled) = Vertex::wrap(moved) {
                        v = scaled;
                    }
                }
            }
        }
        Some(v)
    }

    /// `VertexTol(Index, Param)` (`cxx:290-351`): the gap between the end of
    /// section `I1` and the start of section `I2`, plus `Precision::Confusion()`.
    pub fn vertex_tol_at(&self, index: usize, param: f64) -> f64 {
        let mut tol = CONFUSION;
        let nb = self.edges.len();
        let (i1, i2) = if index == 0 || index == nb {
            if !self.base.is_u_closed() {
                return tol; // The least possible error
            }
            (nb, 1usize)
        } else {
            (index, index + 1)
        };

        let (Some(law1), Some(law2)) = (self.base.law(i1), self.base.law(i2)) else {
            return tol;
        };
        // `PFirst = BS->Value(Knots->Value(Knots->Length()))` of the first law
        // (`cxx:332`), then the distance to `BS->Value(Knots->Value(1))` of the
        // second (`cxx:349`).
        let Ok((bs1, knots1)) = section_bspline(&law1, param) else {
            return tol;
        };
        let Ok((bs2, knots2)) = section_bspline(&law2, param) else {
            return tol;
        };
        let (Some(u_last), Some(u_first)) = (knots1.last().copied(), knots2.first().copied())
        else {
            return tol;
        };
        tol += bs1.d0(u_last).distance(&bs2.d0(u_first));
        tol
    }

    /// `ConcatenedLaw()` (`cxx:355-416`): one law over the concatenation of
    /// all section curves, or the single law when the section has one edge.
    pub fn concatened_law(&self) -> Result<Option<Arc<dyn GeomFillSectionLaw>>, String> {
        if self.base.nb_law() == 1 {
            return Ok(self.base.law(1));
        }
        if !self.shape.is_wire() {
            return Ok(None);
        }

        // Concatenation of edges (`cxx:367-406`).
        let Some(first) = self.edge(1) else {
            return Ok(None);
        };
        let Some(c) = BRepTool::edge_curve_world(&first) else {
            return Ok(None);
        };
        let (f, l) = BRepTool::edge_parameters(&first);
        let tc: Arc<dyn Curve> = Arc::new(GeomTrimmedCurve::new(c, f, l));
        let mut concat =
            CompCurveToBSplineCurve::from_curve(&*tc, ParameterisationType::TgtThetaOver2)
                .map_err(|e| format!("GeomConvert_CompCurveToBSplineCurve: {e:?}"))?;

        let mut bof = true;
        let mut ii = 2usize;
        while ii <= self.edges.len() && bof {
            let (Some(prev), Some(cur)) = (self.edge(ii - 1), self.edge(ii)) else {
                // OCCT dereferences `Edge(ii - 1)` / `Edge(ii)` here; a null
                // slot ends the loop instead.
                break;
            };
            let (Some(c), (f, l)) = (
                BRepTool::edge_curve_world(&cur),
                BRepTool::edge_parameters(&cur),
            ) else {
                // OCCT builds a `Geom_TrimmedCurve` over the null curve.
                break;
            };
            let tc: Arc<dyn Curve> = Arc::new(GeomTrimmedCurve::new(c, f, l));
            let eps = match common_vertex(&prev, &cur) {
                Some(v) => BRepTool::vertex_tolerance(&Vertex(v)),
                None => 10.0 * PCONFUSION,
            };
            bof = concat.add(&*tc, eps, true, false, 20).unwrap_or(false);
            if !bof {
                bof = concat.add(&*tc, 200.0 * eps, true, false, 20).unwrap_or(false);
            }
            ii += 1;
        }

        let Some(composite) = concat.bspline_curve().cloned() else {
            return Ok(None);
        };
        let composite: Arc<dyn Curve> = Arc::new(composite);
        let law: Arc<dyn GeomFillSectionLaw> = match &self.the_law {
            None => Arc::new(UniformSection::new(composite, 0.0, 1.0)?),
            Some(l) => Arc::new(EvolvedSection::new(composite, l.clone())?),
        };
        Ok(Some(law))
    }

    /// `Continuity(Index, TolAngular)` (`cxx:421-474`).
    pub fn continuity_at(&self, index: usize, tol_angular: f64) -> u8 {
        let nb = self.edges.len();
        let (edge1, edge2) = if index == 0 || index == nb {
            if !self.base.is_u_closed() {
                return brep_lprop::C0; // The least possible error
            }
            (self.edge(nb), self.edge(1))
        } else {
            (self.edge(index), self.edge(index + 1))
        };
        let (Some(edge1), Some(edge2)) = (edge1, edge2) else {
            return brep_lprop::C0;
        };

        // Common vertex (`cxx:439-462`): `TopExp::Vertices` is `CumOri = false`.
        let (e1_first, e1_last) = edge_vertices(&edge1);
        let (e2_first, e2_last) = edge_vertices(&edge2);
        let pick = |a: &Option<Vertex>, b: &Option<Vertex>| -> Option<(Vertex, Vertex)> {
            let (a, b) = (a.as_ref()?, b.as_ref()?);
            if is_same(&a.0, &b.0) {
                Some((a.clone(), b.clone()))
            } else {
                None
            }
        };
        let pair = pick(&e1_first, &e2_first)
            .or_else(|| pick(&e1_first, &e2_last))
            .or_else(|| pick(&e1_last, &e2_first))
            .or_else(|| pick(&e1_last, &e2_last));
        let Some((v1, v2)) = pair else {
            return brep_lprop::C0;
        };

        let u1 = BRepTool::parameter_of_vertex(&edge1.0, &v1.0);
        let u2 = BRepTool::parameter_of_vertex(&edge2.0, &v2.0);
        let (Some(c1), Some(c2)) = (
            BRepTool::edge_curve_world(&edge1),
            BRepTool::edge_curve_world(&edge2),
        ) else {
            return brep_lprop::C0;
        };
        let eps = BRepTool::vertex_tolerance(&v2) + BRepTool::vertex_tolerance(&v1);
        brep_lprop::continuity(
            &*c1,
            Orientation::Forward,
            u1,
            (c1.first_parameter(), c1.last_parameter()),
            &*c2,
            Orientation::Forward,
            u2,
            (c2.first_parameter(), c2.last_parameter()),
            false,
            eps,
            tol_angular,
        )
        .unwrap_or(brep_lprop::C0)
    }

    /// `D0(U, S)` (`cxx:478-486`): the section shape, scaled by the law.
    pub fn d0_shape(&self, u: f64) -> Option<TopoShape> {
        let mut s = self.shape.clone();
        if let Some(law) = &self.the_law {
            let mut t = GpTrsf::identity();
            if t.set_scale(&GpPnt::zero(), law.value(u)).is_ok() {
                if let Ok(moved) = transformed_copy(&s, &t) {
                    s = moved;
                }
            }
        }
        Some(s)
    }
}

impl SectionLaw for ShapeLaw {
    fn base(&self) -> &SectionLawBase {
        &self.base
    }

    /// `IsConstant()` (`cxx:239-242`).
    fn is_constant(&self) -> bool {
        self.is_constant_law()
    }

    /// `IsVertex()` (`cxx:232-235`).
    fn is_vertex(&self) -> bool {
        self.vertex
    }

    /// `ConcatenedLaw()` (`cxx:355-416`).
    fn concatened_law(&self) -> Option<Arc<dyn GeomFillSectionLaw>> {
        self.concatened_law().ok().flatten()
    }

    /// `Continuity(Index, TolAngular)` (`cxx:421-474`).
    fn continuity(&self, index: usize, tol_angular: f64) -> u8 {
        self.continuity_at(index, tol_angular)
    }

    /// `VertexTol(Index, Param)` (`cxx:290-351`).
    fn vertex_tol(&self, index: usize, param: f64) -> f64 {
        self.vertex_tol_at(index, param)
    }

    /// `Vertex(Index, Param)` (`cxx:247-285`).
    fn vertex(&self, index: usize, param: f64) -> Option<Vertex> {
        self.vertex_at(index, param)
    }

    /// `D0(U, S)` (`cxx:478-486`).
    fn d0(&self, u: f64) -> Option<TopoShape> {
        self.d0_shape(u)
    }
}

/// Build the `Geom_BSplineCurve` `VertexTol` evaluates (`cxx:317-349`):
/// `Geom_BSplineCurve(Poles, Weights, Knots, Mults, Degree, IsUPeriodic)`.
fn section_bspline(
    law: &Arc<dyn GeomFillSectionLaw>,
    param: f64,
) -> Result<(GeomBSplineCurve, Vec<f64>), String> {
    let (_nb_poles, _nb_knots, degree) = law.section_shape();
    let (poles, weights) = law.d0(param)?;
    let knots = law.knots();
    let mults = law.mults();
    let periodic = law.is_u_periodic();
    let flat = if periodic {
        occt_core::bspl::knots::knot_sequence_periodic(&knots, &mults, degree)
    } else {
        occt_core::bspl::banded_interp::knot_sequence(&knots, &mults, degree)
    };
    let curve = GeomBSplineCurve::from_parts(poles, Some(weights), flat, degree as usize, periodic);
    Ok((curve, knots))
}
