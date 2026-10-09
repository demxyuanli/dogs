use super::*;

/// `double tolint = 1.0e-10;` (`cxx:1105`).
pub(super) const TOLINT: f64 = 1.0e-10;

/// `NCollection_DataMap<TopoDS_Shape, Bnd_Box2d, TopTools_ShapeMapHasher>`
/// (`cxx:1056`): keyed by `IsSame`, i.e. the `TShape` address, which the port
/// models with `GeometryRegistry::shape_key`.
pub(super) type Boxes = HashMap<usize, BndBox2d>;

/// `ShapeFix_IntersectionTool` constructor `(context, preci, maxtol)`
/// (`cxx:51-58`). `ShapeFix_Wire` builds it as
/// `ShapeFix_IntersectionTool ITool(Context(), Precision())`
/// (`ShapeFix_Wire.cxx:1204`), so `maxtol` keeps the C++
/// `ShapeFix_IntersectionTool.hxx:67` default `1.0`.
pub const SHAPE_FIX_INTERSECTION_MAX_TOL: f64 = 1.0;

pub(super) fn skey(e: &Edge) -> usize {
    GeometryRegistry::shape_key(&e.0)
}

/// `BRep_Builder::UpdateVertex(V, Tol)` (`BRep_Builder.cxx:1442-1452`):
/// `BRep_TVertex::UpdateTolerance` is a max (`BRep_TVertex.lxx:33-37`).
pub(super) fn update_vertex_tolerance(v: &Vertex, tol: f64) {
    v.set_tolerance(BRepTool::vertex_tolerance(v).max(tol));
}

pub(super) fn wire_len(wire: &Wire) -> i64 {
    edges_of_wire(wire).len() as i64
}

pub(super) fn wire_edge_at(wire: &Wire, num: i64) -> Option<Edge> {
    if num < 1 {
        return None;
    }
    edges_of_wire(wire).get((num - 1) as usize).cloned()
}

/// `sewd->Add(edge)` / `sewd->Add(edge, at)` (`ShapeExtend_WireData.cxx:400-424`)
/// for an edge taken through the composed accessor: `at == 0` appends.
pub(super) fn wire_add_edge_composed(wire: &mut Wire, at: usize, edge: &Edge) {
    let mut e = edge.clone();
    if wire.0.orientation() == Orientation::Reversed {
        e.0.reverse();
    }
    wire_insert_edge_before(wire, at, &e);
}

/// `BndLib_Add2dCurve::Add(gac, Precision::Confusion(), box)`
/// (`cxx:317-332`): the adaptor range is the pcurve's own `[First, Last]` for a
/// B-spline whose COS range escapes it, `[cf, cl]` otherwise.
pub(super) fn curve2d_box(c2d: &dyn Curve2d, cf: f64, cl: f64) -> BndBox2d {
    let mut box2d = BndBox2d::new();
    let a_first = c2d.first_parameter();
    let a_last = c2d.last_parameter();
    if c2d.is_bspline2d() && (cf < a_first || cl > a_last) {
        // pdn avoiding problems with segment in Bnd_Box (`cxx:328`).
        add_geom2d(c2d, CONFUSION, &mut box2d);
    } else {
        add_geom2d_range(c2d, cf, cl, CONFUSION, &mut box2d);
    }
    box2d
}


/// Out-parameters of `ShapeFix_IntersectionTool::FixSelfIntersectWire`
/// (`cxx:1028-1036`): `NbSplit`, `NbCut`, `NbRemoved` and the function result.
#[derive(Clone, Copy, Default, Debug)]
pub struct FixSelfIntersectWireStatus {
    /// `NbSplit` (`cxx:1061`).
    pub nb_split: i64,
    /// `NbCut` (`cxx:1062`).
    pub nb_cut: i64,
    /// `NbRemoved` (`cxx:1819`).
    pub nb_removed: i64,
    /// `isDone = (NbSplit || NbCut || nbReplaced || NbRemoved)` (`cxx:1828`).
    pub done: bool,
}

/// `myContext->Replace(V12, NewV1)` / `Replace(V12, NewV1.Reversed())`
/// (`cxx:1680-1690`): record the substitution and return the vertex the caller
/// must keep using, orientation matched to the original one.
pub(super) fn replace_vertex_keep_orientation(ctx: &mut dyn ReShape, old: &Vertex, keep: &Vertex) -> Vertex {
    if old.0.orientation() == keep.0.orientation() {
        ctx.replace(&old.0, &keep.0);
        keep.clone()
    } else {
        let rev = Vertex(keep.0.oriented(Orientation::Reversed));
        ctx.replace(&old.0, &rev.0);
        rev
    }
}
