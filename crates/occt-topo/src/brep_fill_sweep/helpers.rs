//! Static helpers of `BRepFill_Sweep` (`BRepFill_Sweep.cxx:96-1470`).
//!
//! Every function here is a file-local `static` of the `.cxx`; the port keeps
//! them in one module because their only consumer is `BRepFill_Sweep`.
//!
//! Translation notes:
//! - `NCollection_HArray2<TopoDS_Shape>` is [`ShapeTable`], whose row/column
//!   indices stay 1-based like `NCollection_Array2::Value(ii, jj)`.
//! - `BRep_TEdge::Curves()` / `ChangeCurves()` are
//!   [`GeometryRegistry::edge_pcurve_reps`], which answers the same
//!   `(surface, pcurve, first, last)` tuples, so `Couture` and `Oriente` keep
//!   the surface-handle lookup of `BRep_Tool::CurveOnSurface`.
//! - A null `TopoDS_Shape` is `None`; `TopoDS_Shape::IsSame` is
//!   [`crate::topo_tools::is_same`] (same TShape).

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::GpVec2d;
use occt_core::kernel::geomabs::{CurveType, Shape};
use occt_core::precision::CONFUSION;
use occt_geom::curve::Curve;
use occt_geom::gcpnts::curve_length_range;
use occt_geom::projlib::curve_type;
use occt_geom::surface::Surface;
use occt_geom2d::Curve2d;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{is_same, shape_map_key, ShapeMapKey};

/// `NCollection_HArray2<TopoDS_Shape>` (`BRepFill_Sweep.cxx:2103-2106`):
/// `Value(row, col)` / `SetValue(row, col, shape)` keep the 1-based bounds of
/// the reference table. A null shape is `None`.
#[derive(Debug, Clone, Default)]
pub struct ShapeTable {
    nb_rows: usize,
    nb_cols: usize,
    values: Vec<Option<TopoShape>>,
}

impl ShapeTable {
    /// `new NCollection_HArray2<TopoDS_Shape>(1, NbRows, 1, NbCols)`.
    pub fn new(nb_rows: usize, nb_cols: usize) -> Self {
        Self {
            nb_rows,
            nb_cols,
            values: vec![None; nb_rows * nb_cols],
        }
    }

    /// `Array2::RowLength()`.
    pub fn nb_rows(&self) -> usize {
        self.nb_rows
    }

    /// `Array2::ColLength()`.
    pub fn nb_cols(&self) -> usize {
        self.nb_cols
    }

    /// `Array2::Value(row, col)` (`row`, `col` 1-based).
    pub fn value(&self, row: usize, col: usize) -> Option<&TopoShape> {
        if row == 0 || col == 0 || row > self.nb_rows || col > self.nb_cols {
            return None;
        }
        self.values[(row - 1) * self.nb_cols + (col - 1)].as_ref()
    }

    /// `Array2::SetValue(row, col, shape)` (`row`, `col` 1-based).
    pub fn set_value(&mut self, row: usize, col: usize, shape: Option<TopoShape>) {
        if row == 0 || col == 0 || row > self.nb_rows || col > self.nb_cols {
            return;
        }
        self.values[(row - 1) * self.nb_cols + (col - 1)] = shape;
    }

    /// `Array2::ChangeValue(row, col)` (`row`, `col` 1-based): the stored shape
    /// is mutated in place, which is how the reference flips an orientation
    /// (`VEdge(isec, j).Reverse()` at `cxx:2390`).
    pub fn value_mut(&mut self, row: usize, col: usize) -> Option<&mut TopoShape> {
        if row == 0 || col == 0 || row > self.nb_rows || col > self.nb_cols {
            return None;
        }
        self.values[(row - 1) * self.nb_cols + (col - 1)].as_mut()
    }

    /// `Array2::ChangeValue(row, col).Reverse()`.
    pub fn reverse_value(&mut self, row: usize, col: usize) {
        if let Some(shape) = self.value_mut(row, col) {
            shape.reverse();
        }
    }

    /// `Array2::Value(row, col).IsNull()`.
    pub fn is_null(&self, row: usize, col: usize) -> bool {
        self.value(row, col).is_none()
    }
}

/// `NumberOfPoles(W)` (`cxx:96-152`): the pole budget of a wire, used by
/// `BuildFace` to skip `BRepLib_FindSurface` for heavy wires.
pub fn number_of_poles(w: &Wire) -> i32 {
    let mut nb_points = 0;
    for e in crate::topo_tools::edges_of(w.shape()) {
        let Some(c) = BRepTool::edge_curve(&e) else {
            continue;
        };
        let (f, l) = BRepTool::edge_parameters(&e);
        if f == l {
            // `cxx:107-111`: degenerate edge.
            continue;
        }
        match curve_type(c.as_ref()) {
            CurveType::BezierCurve => {
                if let Some(poles) = c.bezier_poles() {
                    if poles.len() >= 2 {
                        nb_points += poles.len() as i32;
                    }
                }
            }
            CurveType::BSplineCurve => {
                if let Some(poles) = c.bspline_poles() {
                    if poles.len() >= 2 {
                        nb_points += poles.len() as i32;
                    }
                }
            }
            CurveType::Line => nb_points += 2,
            CurveType::Circle
            | CurveType::Ellipse
            | CurveType::Hyperbola
            | CurveType::Parabola => nb_points += 4,
            _ => nb_points += 15 + c.nb_intervals(Shape::C3 as u8),
        }
    }
    nb_points
}

/// `HasPCurves(E)` (`cxx:157-173`): true when the edge carries at least one
/// `CurveOnSurface` representation.
pub fn has_pcurves(e: &Edge) -> bool {
    !GeometryRegistry::global().edge_pcurve_reps(&e.0).is_empty()
}

/// `Translate(ArrayIn, In, ArrayOut, Out)` (`cxx:177-190`): copy a column.
pub fn translate(array_in: &ShapeTable, in_col: usize, array_out: &mut ShapeTable, out_col: usize) {
    let nb = array_out.nb_cols();
    for ii in 1..=nb {
        array_out.set_value(ii, out_col, array_in.value(ii, in_col).cloned());
    }
}

/// The `(pcurve, first, last)` of `E` on the surface `S`, using the
/// `BRep_CurveRepresentation` lookup of `BRep_Tool::CurveOnSurface`.
fn pcurve_on_surface(
    e: &Edge,
    s: &Arc<dyn Surface>,
) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    let reps = GeometryRegistry::global().edge_pcurve_reps(&e.0);
    let key = Arc::as_ptr(s) as *const () as usize;
    let on_surface: Vec<(Arc<dyn Curve2d>, f64, f64)> = reps
        .into_iter()
        .filter(|(surf, _, _, _)| Arc::as_ptr(surf) as *const () as usize == key)
        .map(|(_, pc, a, b)| (pc, a, b))
        .collect();
    if on_surface.len() >= 2 && e.0.orientation() == Orientation::Reversed {
        return on_surface.into_iter().nth(1);
    }
    on_surface.into_iter().next()
}

/// `Couture(E, S, L)` (`cxx:212-245`): the pcurve of `E` on the surface `S`.
///
/// OCCT asks for `BRep_GCurve::PCurve2()` when the representation is a
/// `CurveOnClosedSurface` and the edge is `REVERSED`, and `PCurve()` otherwise.
/// The port stores the two pcurves of a closed-surface representation in
/// forward-then-reversed order (`GeometryRegistry::edge_pcurves`), which is the
/// same `PCurve()`/`PCurve2()` pairing.
pub fn couture(e: &Edge, s: &Arc<dyn Surface>) -> Option<Arc<dyn Curve2d>> {
    pcurve_on_surface(e, s).map(|(pc, _, _)| pc)
}

/// `IsDegen(S, Tol)` (`cxx:1823-1875`): true when the surface is squeezed to a
/// line, tested first along the iso-U direction then along iso-V.
pub fn is_degen(s: &Arc<dyn Surface>, tol: f64) -> bool {
    let nb = 5;
    let (umin, umax) = s.u_range();
    let (vmin, vmax) = s.v_range();

    // Check the length of Iso-U (`cxx:1837-1855`).
    let t = (umin + umax) / 2.0;
    let p1 = s.d0(t, vmin);
    let p2 = s.d0(t, (vmin + vmax) / 2.0);
    let p3 = s.d0(t, vmax);
    let mut b = (p1.distance(&p2) + p2.distance(&p3)) < tol;
    if b {
        let dt = (umax - umin) / (nb as f64 + 1.0);
        for ii in 1..=nb {
            let tt = umin + ii as f64 * dt;
            let Some(iso) = s.u_iso_curve(tt) else {
                b = false;
                break;
            };
            let l = curve_length_range(iso.as_ref(), iso.first_parameter(), iso.last_parameter(), tol / 4.0);
            b = l <= tol;
            if !b {
                break;
            }
        }
    }
    if b {
        return true;
    }

    // Check the length of Iso-V (`cxx:1858-1873`).
    let t = (vmin + vmax) / 2.0;
    let p1 = s.d0(umin, t);
    let p2 = s.d0((umin + umax) / 2.0, t);
    let p3 = s.d0(umax, t);
    b = (p1.distance(&p2) + p2.distance(&p3)) < tol;
    if b {
        let dt = (vmax - vmin) / (nb as f64 + 1.0);
        for ii in 1..=nb {
            let tt = vmin + ii as f64 * dt;
            let Some(iso) = s.v_iso_curve(tt) else {
                b = false;
                break;
            };
            let l = curve_length_range(iso.as_ref(), iso.first_parameter(), iso.last_parameter(), tol / 4.0);
            b = l <= tol;
            if !b {
                break;
            }
        }
    }
    b
}

/// `KeepEdge(Face, Edge, List)` (`cxx:1413-1434`): the edges of `Face` whose
/// 3D curve is the very same `Geom_Curve` handle (and location) as `Edge`'s.
pub fn keep_edge(face: &TopoShape, edge: &Edge) -> Vec<TopoShape> {
    let mut list = Vec::new();
    let Some(cref) = BRepTool::edge_curve(edge) else {
        return list;
    };
    for e in crate::topo_tools::edges_of(face) {
        let Some(c) = BRepTool::edge_curve(&e) else {
            continue;
        };
        if Arc::ptr_eq(&cref, &c) && e.location().is_equal(edge.location()) {
            list.push(e.into());
        }
    }
    list
}

/// `ReverseEdgeInFirstOrLastWire(theWire, theEdge)` (`cxx:1877-1900`): remove
/// the edge from the wire, reverse it and append it back at the end.
pub fn reverse_edge_in_first_or_last_wire(the_wire: &mut TopoShape, the_edge: &TopoShape) {
    let mut edge_to_reverse = None;
    for an_edge in crate::topo_tools::edges_of(the_wire) {
        if is_same(&an_edge, the_edge) {
            edge_to_reverse = Some(an_edge);
            break;
        }
    }
    if let Some(mut e) = edge_to_reverse {
        let b = TopoBuilder::new();
        the_wire.set_free(true);
        b.remove(the_wire, &e);
        e.0.reverse();
        b.add(the_wire, &e);
    }
}

/// `ReverseModifiedEdges(theWire, theEmap)` (`cxx:1902-1936`): re-add every
/// wire edge, reversed when it is in the `myVEdgesModified` data map.
pub fn reverse_modified_edges(the_wire: &mut TopoShape, the_emap: &HashMap<ShapeMapKey, TopoShape>) {
    if the_emap.is_empty() {
        return;
    }
    reverse_modified_edges_with(the_wire, &|e| the_emap.contains_key(&shape_map_key(e)));
}

/// The `NCollection_Map<TopoDS_Shape>` overload of `ReverseModifiedEdges`
/// (`cxx:1902-1936`), used for the `ReversedEdges` set of `BuildShell`.
pub fn reverse_modified_edges_set(
    the_wire: &mut TopoShape,
    the_emap: &std::collections::HashSet<ShapeMapKey>,
) {
    if the_emap.is_empty() {
        return;
    }
    reverse_modified_edges_with(the_wire, &|e| the_emap.contains(&shape_map_key(e)));
}

/// The shared body of the two `ReverseModifiedEdges` overloads.
fn reverse_modified_edges_with(the_wire: &mut TopoShape, contains: &dyn Fn(&TopoShape) -> bool) {
    let ledges = crate::topo_tools::edges_of(the_wire);
    let b = TopoBuilder::new();
    the_wire.set_free(true);
    for e in &ledges {
        b.remove(the_wire, e);
    }
    for e in ledges {
        let mut an_edge = e;
        if contains(&an_edge) {
            an_edge.0.reverse();
        }
        b.add(the_wire, &an_edge);
    }
}

/// `Oriente(S, E)` (`cxx:455-491`): orient a natural-restriction edge so that
/// it runs along the increasing iso parameter on the boundary of `S`.
pub fn oriente(s: &Arc<dyn Surface>, e: &mut Edge) {
    let (ufirst, _ulast) = s.u_range();
    let (vfirst, _vlast) = s.v_range();
    let Some((c, f, l)) = pcurve_on_surface(e, s) else {
        return;
    };
    let (p, d) = c.d1((f + l) / 2.0);
    let v_ref = GpVec2d::new(0.0, 1.0);
    let u_ref = GpVec2d::new(1.0, 0.0);
    let isuiso = d.is_parallel(&v_ref, 0.1);
    let (isfirst, isopposite);
    if isuiso {
        isfirst = (p.x() - ufirst).abs() < CONFUSION;
        isopposite = d.is_opposite(&v_ref, 0.1);
        e.0.set_orientation(Orientation::Reversed);
    } else {
        isfirst = (p.y() - vfirst).abs() < CONFUSION;
        isopposite = d.is_opposite(&u_ref, 0.1);
        e.0.set_orientation(Orientation::Forward);
    }
    if !isfirst {
        e.0.reverse();
    }
    if isopposite {
        e.0.reverse();
    }
}

/// `BuildVertex(Iso, isfirst, First, Last, Vertex)` (`cxx:1440-1458`): move
/// `Vertex` to the iso curve end point.
pub fn build_vertex(iso: &Arc<dyn Curve>, isfirst: bool, first: f64, last: f64, vertex: &mut Vertex) {
    let val = if isfirst { first } else { last };
    let p = iso.value(val);
    let b = TopoBuilder::new();
    *vertex = b.make_vertex(p, CONFUSION);
}

/// `NullEdge(Vertex)` (`cxx:1462-1472`): a degenerate edge between `Vertex`
/// and its reversal.
pub fn null_edge(vertex: &mut Vertex) -> Edge {
    let b = TopoBuilder::new();
    let mut e = Edge::new();
    vertex.0.set_orientation(Orientation::Forward);
    b.add(&mut e.0, &vertex.0);
    let mut reversed = vertex.clone();
    reversed.0.reverse();
    b.add(&mut e.0, &reversed.0);
    GeometryRegistry::global().set_degenerated(&e.0, true);
    e
}
