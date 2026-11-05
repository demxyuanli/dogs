//! `BRepFill_NSections` — the section law of N section wires.
//!
//! Source: `BRepFill_NSections.cxx:1-881`, `.hxx:39-86`.
//!
//! The class builds one `GeomFill_NSections` law per non-degenerated section
//! edge, all of them sharing the single `Geom_BSplineSurface` that
//! [`total_surf`] approximates from the whole stack of sections.
//!
//! Translation notes:
//! - The inherited `BRepFill_SectionLaw` data lives in [`SectionLawBase`];
//!   `VFirst`, `VLast`, `myShapes`, `myTrsfs`, `myParams`, `myEdges` and
//!   `mySurface` (`hxx:76-83`) are this class's own fields.
//! - `NCollection_HArray2<TopoDS_Shape> myEdges(1, NbEdge, 1, NbSects)` is
//!   `Vec<Vec<Option<TopoShape>>>` indexed `[i - 1][j - 1]`.
//! - `myIndices.Bind(E, -ii)` for a `REVERSED` section edge
//!   (`BRepFill_NSections.cxx:489-495`) is why
//!   [`SectionLawBase::bind_index`] takes a signed index.
//! - `BRep_Tool::Curve(E, First, Last)` is [`BRepTool::edge_curve`] with
//!   [`BRepTool::edge_parameters`], i.e. the **local** curve, exactly as
//!   `BRepFill_ShapeLaw` uses it: the laws of one surface must stay in the same
//!   space as the `gp_Trsf` each of them carries.
//! - `GeomConvert_ApproxCurve` / `GeomFill_AppSurf` failures raise
//!   `StdFail_NotDone` in OCCT; the port returns `Err`. `D0`, `Vertex`,
//!   `VertexTol`, `Continuity` and `ConcatenedLaw` have no error channel in the
//!   [`SectionLaw`] trait, so they answer `None` / the `GeomAbs_C0` floor where
//!   OCCT would raise.
//! - `OCCT_DEBUG` in this file only carries an empty `if (Affich)`
//!   (`cxx:392-398`), so nothing is skipped there.

mod edge_to_bspline;
mod totalsurf;

use std::sync::Arc;

use occt_core::gp::GpTrsf;
use occt_core::precision::CONFUSION;
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::bspline_surface::GeomBSplineSurface;
use occt_geom::curve::Curve;
use occt_geom::geom_fill::SectionLaw as GeomFillSectionLaw;
use occt_geom::geom_fill::NSections as GeomFillNSections;
use occt_geom::trimmed::GeomTrimmedCurve;
use occt_geom::Surface;

use crate::abs::Orientation;
use crate::brep_fill_section_law::{SectionLaw, SectionLawBase};
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_lprop;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::topo_tools::{edge_vertices, is_same};

pub use edge_to_bspline::{edge_to_bspline, point_curve_of_edge};
pub use totalsurf::total_surf;

/// `BRepFill_NSections` (`BRepFill_NSections.hxx:39-86`).
pub struct NSections {
    /// The inherited `BRepFill_SectionLaw` state.
    base: SectionLawBase,
    /// `VFirst` (`hxx:76`).
    v_first: f64,
    /// `VLast` (`hxx:77`).
    v_last: f64,
    /// `myShapes` (`hxx:78`).
    shapes: Vec<TopoShape>,
    /// `myTrsfs` (`hxx:79`).
    trsfs: Vec<GpTrsf>,
    /// `myParams` (`hxx:80`).
    params: Vec<f64>,
    /// `myEdges` (`hxx:81`): edges by 1-based (edge, section) index.
    edges: Vec<Vec<Option<TopoShape>>>,
    /// `mySurface` (`hxx:82`).
    surface: Option<Arc<GeomBSplineSurface>>,
}

impl NSections {
    /// `BRepFill_NSections(S, Build = true)` (`cxx:363-381`): the parameters
    /// are `0, 1, .., Length - 1`.
    pub fn with_sections(shapes: Vec<TopoShape>, build: bool) -> Result<Self, String> {
        let mut this = Self::new_shapes(shapes);
        let par: Vec<f64> = (0..this.shapes.len()).map(|i| i as f64).collect();
        this.params = par.clone();
        this.init(&par, build)?;
        this.base.set_done(true);
        Ok(this)
    }

    /// `BRepFill_NSections(S, Trsfs, P, VF, VL, Build = true)`
    /// (`cxx:387-419`): `myDone` stays false when `P` is not strictly
    /// increasing, and then nothing is built.
    #[allow(clippy::too_many_arguments)]
    pub fn with_transformations(
        shapes: Vec<TopoShape>,
        trsfs: Vec<GpTrsf>,
        params: Vec<f64>,
        v_first: f64,
        v_last: f64,
        build: bool,
    ) -> Result<Self, String> {
        let mut this = Self::new_shapes(shapes);
        // `cxx:399-404`: `P.Value(iseq) < P.Value(iseq + 1)` for every pair.
        if !params.windows(2).all(|w| w[0] < w[1]) {
            return Ok(this);
        }
        this.params = params.clone();
        this.trsfs = trsfs;
        this.v_first = v_first;
        this.v_last = v_last;
        this.init(&params, build)?;
        this.base.set_done(true);
        Ok(this)
    }

    /// The common prologue of the two constructors: `myShapes = S`.
    fn new_shapes(shapes: Vec<TopoShape>) -> Self {
        Self {
            base: SectionLawBase::new(),
            v_first: 0.0,
            v_last: 1.0,
            shapes,
            trsfs: Vec::new(),
            params: Vec::new(),
            edges: Vec::new(),
            surface: None,
        }
    }

    /// `IsVertex()` (`cxx:655-658`).
    pub fn is_vertex_law(&self) -> bool {
        false
    }

    /// `IsConstant()` (`cxx:662-665`).
    pub fn is_constant_law(&self) -> bool {
        false
    }

    /// `myEdges->ColLength()`, the number of law columns.
    pub fn nb_edges(&self) -> usize {
        self.edges.len()
    }

    /// `mySurface` (`hxx:82`).
    pub fn surface(&self) -> Option<Arc<GeomBSplineSurface>> {
        self.surface.clone()
    }

    /// `myShapes` (`hxx:78`).
    pub fn shapes(&self) -> &[TopoShape] {
        &self.shapes
    }

    /// `myParams` (`hxx:80`).
    pub fn params(&self) -> &[f64] {
        &self.params
    }

    fn wire_at(&self, index: usize) -> Result<Wire, String> {
        self.shapes
            .get(index.wrapping_sub(1))
            .cloned()
            .and_then(Wire::wrap)
            .ok_or_else(|| {
                "Standard_TypeMismatch: BRepFill_NSections section is not a wire".to_string()
            })
    }

    /// `myEdges->Value(i, j)`, 1-based on both axes.
    fn edge_at(&self, i: usize, j: usize) -> Result<Edge, String> {
        self.edges
            .get(i.wrapping_sub(1))
            .and_then(|row| row.get(j.wrapping_sub(1)))
            .and_then(|slot| slot.clone())
            .and_then(Edge::wrap)
            .ok_or_else(|| "Standard_NullObject: BRepFill_NSections::myEdges".to_string())
    }

    /// The non-degenerated edges of section `index`, in wire-explorer order.
    fn section_edges(&self, index: usize) -> Result<Vec<Edge>, String> {
        let wire = self.wire_at(index)?;
        let mut out = Vec::new();
        let mut wexp = WireExplorer::new();
        wexp.init(&wire, None);
        while wexp.more() {
            let edge = wexp
                .current_edge()
                .ok_or_else(|| "Standard_NullObject: BRepFill_NSections".to_string())?;
            out.push(edge);
            wexp.next();
        }
        Ok(out)
    }

    /// `Init(P, Build)` (`cxx:424-651`).
    fn init(&mut self, params: &[f64], build: bool) -> Result<(), String> {
        let nb_sects = params.len();
        if nb_sects == 0 || self.shapes.len() < nb_sects {
            return Err("Standard_OutOfRange: BRepFill_NSections::Init".to_string());
        }

        let mut ideb = 1usize;
        let mut ifin = nb_sects;

        // `cxx:436-446`: is the first section punctual?
        let mut w1_point = true;
        for edge in self.section_edges(1)? {
            w1_point = w1_point && BRepTool::is_degenerated(&edge);
        }
        if w1_point {
            ideb += 1;
        }

        // `cxx:447-457`: and the last one?
        let mut w2_point = true;
        for edge in self.section_edges(nb_sects)? {
            w2_point = w2_point && BRepTool::is_degenerated(&edge);
        }
        if w2_point {
            ifin -= 1;
        }

        // `cxx:459`: `vclosed = myShapes(1).IsSame(myShapes(NbSects))`.
        let v_closed = is_same(&self.shapes[0], &self.shapes[nb_sects - 1]);
        self.base.set_v_closed(v_closed);

        // `cxx:461-471`: the non-degenerated edge count of section `ideb`.
        let nb_edge = self
            .section_edges(ideb)?
            .iter()
            .filter(|e| !BRepTool::is_degenerated(e))
            .count();

        // `cxx:473`: `myEdges = new HArray2(1, NbEdge, 1, NbSects)`.
        self.edges = vec![vec![None; nb_sects]; nb_edge];
        self.base.reset_laws(nb_edge);

        // `cxx:476-546`: fill the table, the edge -> signed index map, and the
        // U-closedness flag.
        let mut u_closed = true;
        for jj in ideb..=ifin {
            let wire = self.wire_at(jj)?;
            let mut wexp = WireExplorer::new();
            wexp.init(&wire, None);
            let mut ii = 1usize;
            while ii <= nb_edge && wexp.more() {
                let edge = wexp
                    .current_edge()
                    .ok_or_else(|| "Standard_NullObject: BRepFill_NSections::Init".to_string())?;
                if !BRepTool::is_degenerated(&edge) {
                    self.edges[ii - 1][jj - 1] = Some(edge.0.clone());
                    let index = if edge.0.orientation() == Orientation::Forward {
                        ii as i32
                    } else {
                        -(ii as i32)
                    };
                    self.base.bind_index(&edge.0, index);
                }
                wexp.next();
                ii += 1;
            }

            // `cxx:501-543`.
            let mut w_closed = wire.is_closed();
            if !w_closed {
                let edge1 = self.edge_at(nb_edge, jj)?;
                let edge2 = self.edge_at(1, jj)?;
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
                let (Some(v1), Some(v2)) = (v1, v2) else {
                    return Err(
                        "Standard_NullObject: BRepFill_NSections::Init boundary vertex".to_string()
                    );
                };
                if is_same(&v1.0, &v2.0) {
                    w_closed = true;
                } else {
                    let u1 = BRepTool::parameter_of_vertex(&edge1.0, &v1.0);
                    let u2 = BRepTool::parameter_of_vertex(&edge2.0, &v2.0);
                    let eps = BRepTool::vertex_tolerance(&v2) + BRepTool::vertex_tolerance(&v1);
                    let (Some(c1), Some(c2)) = (
                        BRepTool::edge_curve_world(&edge1),
                        BRepTool::edge_curve_world(&edge2),
                    ) else {
                        return Err(
                            "Standard_NullObject: BRepFill_NSections::Init curve".to_string()
                        );
                    };
                    w_closed = c1.d0(u1).is_equal_with_tolerance(&c2.d0(u2), eps);
                }
            }
            if !w_closed {
                u_closed = false;
            }
        }
        self.base.set_u_closed(u_closed);

        // `cxx:548-558`: a punctual first section replicates its single edge
        // over the whole column.
        if w1_point {
            let first = self.section_edges(1)?.into_iter().next().ok_or_else(|| {
                "Standard_NullObject: BRepFill_NSections::Init punctual wire".to_string()
            })?;
            for ii in 1..=nb_edge {
                self.edges[ii - 1][0] = Some(first.0.clone());
            }
        }

        // `cxx:560-569`.
        if w2_point {
            let last = self
                .section_edges(nb_sects)?
                .into_iter()
                .next()
                .ok_or_else(|| {
                    "Standard_NullObject: BRepFill_NSections::Init punctual wire".to_string()
                })?;
            for ii in 1..=nb_edge {
                self.edges[ii - 1][nb_sects - 1] = Some(last.0.clone());
            }
        }

        // `cxx:571-591`: the shared surface, with the V degree raised to 2 so
        // that `GeomFill_NSections` can answer D2 correctly.
        let mut surface = total_surf(
            &self.edges,
            self.shapes.len(),
            nb_edge,
            &self.params,
            w1_point,
            w2_point,
            u_closed,
            v_closed,
            CONFUSION,
        )?;
        if surface.deg_v() < 2 {
            surface.increase_degree(surface.deg_u(), 2)?;
        }
        let surface = Arc::new(surface);
        self.surface = Some(Arc::clone(&surface));

        // `cxx:593-650`: one law per column.
        if build {
            for ii in 1..=nb_edge {
                let mut curves: Vec<Arc<dyn Curve>> = Vec::with_capacity(nb_sects);
                for jj in 1..=nb_sects {
                    let edge = self.edge_at(ii, jj)?;
                    curves.push(self.section_curve(&edge, ii)?);
                }
                // `cxx:644-648`.
                let law = GeomFillNSections::with_reference_surface(
                    curves,
                    self.trsfs.clone(),
                    self.params.clone(),
                    (ii - 1) as f64,
                    ii as f64,
                    self.v_first,
                    self.v_last,
                    Arc::clone(&surface),
                )?;
                self.base.set_law_at(ii, Some(Arc::new(law)));
            }
        }

        Ok(())
    }

    /// The per-edge curve of `cxx:600-643`: a point curve for a degenerated
    /// edge, otherwise the edge's own curve, reversed and trimmed.
    fn section_curve(&self, edge: &Edge, ii: usize) -> Result<Arc<dyn Curve>, String> {
        if BRepTool::is_degenerated(edge) {
            // `cxx:604-621`.
            return Ok(Arc::new(point_curve_of_edge(edge)?));
        }
        let mut curve = BRepTool::edge_curve(edge)
            .ok_or_else(|| "Standard_NullObject: BRep_Tool::Curve".to_string())?;
        let (mut first, mut last) = BRepTool::edge_parameters(edge);
        if edge.0.orientation() == Orientation::Reversed {
            // `cxx:627-635`: `C->Reversed()` with both parameters run through
            // `ReversedParameter`.
            let aux = curve.reversed_parameter(first);
            first = curve.reversed_parameter(last);
            last = aux;
            curve = Arc::from(curve.reversed());
        }
        // `cxx:636-639`: cut, unless this is the seam edge of a closed section.
        if ii > 1 || !edge.0.is_closed() {
            Ok(Arc::new(GeomTrimmedCurve::new(curve, first, last)))
        } else {
            Ok(curve)
        }
    }

    /// `Vertex(Index, Param)` (`cxx:669-695`).
    pub fn vertex_at(&self, index: usize, param: f64) -> Option<Vertex> {
        let nb = self.nb_edges();
        let (law_index, at_last) = if index <= nb {
            (index, false)
        } else if index == nb + 1 {
            (index - 1, true)
        } else {
            return None;
        };
        let law = self.base.law(law_index)?;
        let curve = law.bspline_surface()?.v_iso(param).ok()?;
        let u = if at_last {
            curve.last_parameter()
        } else {
            curve.first_parameter()
        };
        // `B.UpdateVertex(V, P, Precision::Confusion())`.
        Some(TopoBuilder::new().make_vertex(curve.d0(u), CONFUSION))
    }

    /// `VertexTol(Index, Param)` (`cxx:700-760`).
    pub fn vertex_tol_at(&self, index: usize, param: f64) -> f64 {
        let mut tol = CONFUSION;
        let nb = self.nb_edges();
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
        let (Some((bs1, knots1)), Some((bs2, knots2))) = (
            law_section_bspline(&law1, param),
            law_section_bspline(&law2, param),
        ) else {
            return tol;
        };
        // `PFirst = BS->Value(Knots->Value(Knots->Length()))` of the first law
        // (`cxx:741`), the distance to `BS->Value(Knots->Value(1))` of the
        // second (`cxx:759`).
        let (Some(u_last), Some(u_first)) = (knots1.last().copied(), knots2.first().copied())
        else {
            return tol;
        };
        tol += bs1.d0(u_last).distance(&bs2.d0(u_first));
        tol
    }

    /// `ConcatenedLaw()` (`cxx:764-786`).
    pub fn concatened_law(&self) -> Option<Arc<dyn GeomFillSectionLaw>> {
        if self.base.nb_law() == 1 {
            return self.base.law(1);
        }
        let surface = self.surface.clone()?;
        let (u_first, u_last) = surface.u_range();
        let (v_first, v_last) = surface.v_range();
        let mut curves: Vec<Arc<dyn Curve>> = Vec::with_capacity(self.shapes.len());
        for jj in 1..=self.shapes.len() {
            let param = *self.params.get(jj - 1)?;
            curves.push(Arc::new(surface.v_iso(param).ok()?));
        }
        let law = GeomFillNSections::with_reference_surface(
            curves,
            self.trsfs.clone(),
            self.params.clone(),
            u_first,
            u_last,
            v_first,
            v_last,
            surface,
        )
        .ok()?;
        Some(Arc::new(law))
    }

    /// `Continuity(Index, TolAngular)` (`cxx:790-859`): the worst continuity
    /// over all the sections.
    pub fn continuity_at(&self, index: usize, tol_angular: f64) -> u8 {
        let nb = self.nb_edges();
        let mut cont = brep_lprop::C0;
        for jj in 1..=self.shapes.len() {
            let (edge1, edge2) = if index == 0 || index == nb {
                if !self.base.is_u_closed() {
                    return brep_lprop::C0; // The least possible error
                }
                (self.edge_at(nb, jj), self.edge_at(1, jj))
            } else {
                (self.edge_at(index, jj), self.edge_at(index + 1, jj))
            };
            let (Ok(edge1), Ok(edge2)) = (edge1, edge2) else {
                return brep_lprop::C0;
            };

            // `cxx:818-833`: `TopExp::FirstVertex` / `LastVertex`, `CumOri = false`.
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

            let cont_jj = if BRepTool::is_degenerated(&edge1)
                || BRepTool::is_degenerated(&edge2)
            {
                // `cxx:835-838`.
                brep_lprop::CN
            } else {
                let (Some(v1), Some(v2)) = (v1, v2) else {
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
                // `cxx:846`: `BRepLProp::Continuity(Curve1, Curve2, U1, U2,
                // Eps, TolAngular)`.
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
            };

            if jj == 1 {
                cont = cont_jj;
            }
            if cont > cont_jj {
                cont = cont_jj;
            }
        }
        cont
    }

    /// `D0(V, S)` (`cxx:864-881`): the section wire at `V`.
    pub fn d0_shape(&self, v: f64) -> Option<TopoShape> {
        let mut mw = MakeWire::new();
        for ii in 1..=self.base.nb_law() {
            let law = self.base.law(ii)?;
            let iso = law.bspline_surface()?.v_iso(v).ok()?;
            let first = iso.first_parameter();
            let last = iso.last_parameter();
            // `BRepLib_MakeEdge(Curve, first, last)` (`BRepLib_MakeEdge.cxx:217-230`
            // takes the `Init(C, p1, p2)` arm) then `MW.Add(E)`.
            let arc: Arc<dyn Curve> = Arc::new(iso);
            let builder = TopoBuilder::new();
            let build = crate::brep_lib_make_edge::init_3d(&arc, first, last)?;
            let mut edge = builder.make_edge(arc, build.first, build.last);
            crate::brep_lib_make_edge::add_endpoint_vertices(&builder, &mut edge, &build);
            // OCCT's `MW.Add` only records `myError` on failure; the wire is
            // still the one `MW.Wire()` hands back (`BRepLib_MakeWire.cxx:110-453`).
            let _ = mw.add(&edge);
        }
        let wire = mw.wire();
        let oriented = wire.0.oriented(Orientation::Forward);
        Some(oriented)
    }
}

impl SectionLaw for NSections {
    fn base(&self) -> &SectionLawBase {
        &self.base
    }

    /// `IsConstant()` (`cxx:662-665`).
    fn is_constant(&self) -> bool {
        self.is_constant_law()
    }

    /// `IsVertex()` (`cxx:655-658`).
    fn is_vertex(&self) -> bool {
        self.is_vertex_law()
    }

    /// `ConcatenedLaw()` (`cxx:764-786`).
    fn concatened_law(&self) -> Option<Arc<dyn GeomFillSectionLaw>> {
        self.concatened_law()
    }

    /// `Continuity(Index, TolAngular)` (`cxx:790-859`).
    fn continuity(&self, index: usize, tol_angular: f64) -> u8 {
        self.continuity_at(index, tol_angular)
    }

    /// `VertexTol(Index, Param)` (`cxx:700-760`).
    fn vertex_tol(&self, index: usize, param: f64) -> f64 {
        self.vertex_tol_at(index, param)
    }

    /// `Vertex(Index, Param)` (`cxx:669-695`).
    fn vertex(&self, index: usize, param: f64) -> Option<Vertex> {
        self.vertex_at(index, param)
    }

    /// `D0(U, S)` (`cxx:864-881`).
    fn d0(&self, u: f64) -> Option<TopoShape> {
        self.d0_shape(u)
    }
}

/// The `Geom_BSplineCurve` `VertexTol` assembles from one law's section data
/// (`cxx:727-742` and `:745-759`):
/// `Geom_BSplineCurve(Poles, Weights, Knots, Mults, Degree, IsUPeriodic)`.
fn law_section_bspline(
    law: &Arc<dyn GeomFillSectionLaw>,
    param: f64,
) -> Option<(GeomBSplineCurve, Vec<f64>)> {
    let (_nb_poles, _nb_knots, degree) = law.section_shape();
    let (poles, weights) = law.d0(param).ok()?;
    let knots = law.knots();
    let mults = law.mults();
    let periodic = law.is_u_periodic();
    let flat = if periodic {
        occt_core::bspl::knots::knot_sequence_periodic(&knots, &mults, degree)
    } else {
        occt_core::bspl::banded_interp::knot_sequence(&knots, &mults, degree)
    };
    let curve = GeomBSplineCurve::from_parts(
        poles,
        Some(weights),
        flat,
        degree as usize,
        periodic,
    );
    Some((curve, knots))
}
