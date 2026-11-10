use occt_core::gp::GpPnt2d;

use crate::abs::Orientation;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;

use super::BRepTool;

impl BRepTool {
    /// `BRep_Tool::Range(E, S, L, f, l)` (`BRep_Tool.cxx:1004-1013`): the
    /// stored p-curve's range, or the 3D edge range when no p-curve exists.
    pub(super) fn pcurve_range(e: &Edge, f: &Face) -> (f64, f64) {
        let face_key = GeometryRegistry::shape_key(&f.0);
        if let Some((t1, t2)) = GeometryRegistry::global().pcurve_range(&e.0, face_key) {
            if t1.is_finite() && t2.is_finite() {
                return (t1, t2);
            }
        }
        Self::edge_parameters(e)
    }

    /// `BRep_Tool::UVPoints(E, F, PFirst, PLast)` (`BRep_Tool.cxx:1090-1106`):    /// the UV coordinates of the edge's first and last vertices on the face's
    /// surface.
    ///
    /// `cxx:1096-1101` reverses a local copy of the edge when the face is
    /// `REVERSED`; the representation walk (`cxx:1017-1036`) then answers
    /// `UVPoints2` (the seam p-curve) when that local edge is `REVERSED` and
    /// the representation is a closed-surface one. Both p-curves are evaluated
    /// at the same `First`/`Last` (`BRep_CurveOnClosedSurface.cxx:123-134`).
    ///
    /// The `Geom_RectangularTrimmedSurface` unwrapping of `cxx:1042-1053` is
    /// folded into `gp_pln`, which already answers the basis plane of a
    /// trimmed plane.
    pub fn uv_points(edge: &Edge, face: &Face) -> (GpPnt2d, GpPnt2d) {
        use occt_core::precision::Precision;

        let face_reversed = face.0.orientation() == Orientation::Reversed;
        let edge_reversed = edge.0.orientation().is_reversed();
        // `BRep_Tool.cxx:1098-1101`: `aLocalEdge.Reverse()`.
        let local_reversed = edge_reversed != face_reversed;

        let face_key = GeometryRegistry::shape_key(&face.0);
        let pcs = GeometryRegistry::global().edge_pcurves(&edge.0, face_key);
        if !pcs.is_empty() {
            let pc = if pcs.len() >= 2 && local_reversed {
                &pcs[1]
            } else {
                &pcs[0]
            };
            let (first, last) = GeometryRegistry::global()
                .pcurve_range(&edge.0, face_key)
                .unwrap_or_else(|| Self::edge_parameters(edge));
            let mut pfirst = GpPnt2d::new(0.0, 0.0);
            let mut plast = GpPnt2d::new(0.0, 0.0);
            if !Precision::is_negative_infinite(first) {
                pfirst = pc.d0(first);
            }
            if !Precision::is_positive_infinite(last) {
                plast = pc.d0(last);
            }
            return (pfirst, plast);
        }

        // `cxx:1039-1083`: no representation -- for a plane, project the two
        // vertices through `ElSLib::Parameters`.
        if let Some(surf) = GeometryRegistry::global().face_surface(&face.0) {
            if let Some(pln) = surf.gp_pln() {
                let pos = pln.position();
                let loc = face.0.location();
                let invert = if loc.is_identity() {
                    None
                } else {
                    loc.transformation().inverted().ok()
                };
                // `cxx:1059-1061`: `Vf.Move(L.Inverted())`.
                let project = |v: &Vertex| -> GpPnt2d {
                    let p = match &invert {
                        Some(t) => Self::vertex_point_world(v).transformed(t),
                        None => Self::vertex_point(v),
                    };
                    let (u, vv) = occt_core::elib::slib::plane_parameters(&pos, &p);
                    GpPnt2d::new(u, vv)
                };
                let (vf, vl) = crate::topo_tools::edge_vertices(edge);
                let pfirst = match vf {
                    Some(v) => project(&v),
                    None => GpPnt2d::new(0.0, 0.0),
                };
                let plast = match vl {
                    Some(v) => project(&v),
                    None => GpPnt2d::new(0.0, 0.0),
                };
                return (pfirst, plast);
            }
        }
        (GpPnt2d::new(0.0, 0.0), GpPnt2d::new(0.0, 0.0))
    }

    /// `BRep_Tool::Parameters(V, F)` (`BRep_Tool.cxx:1661-1703`): the UV
    /// coordinates of the vertex `V` in the parametric space of the face `F`.
    ///
    /// `cxx:1663-1677` first looks for a `BRep_PointRepresentation` of `V` that
    /// lies on the face's surface and answers its `Parameter` / `Parameter2`
    /// ("case non Manifold"); only `BRep_PointOnSurface` answers
    /// `IsPointOnSurface` true (`BRep_PointRepresentation.cxx:53-56`,
    /// `BRep_PointOnSurface.cxx:50-53`), so the pair is the stored UV of that
    /// representation. `TopLoc_Location` is unported - the predicate is the
    /// surface handle identity (`brep_point_representation` module note).
    ///
    /// When no representation matches, the search falls to the face's edges
    /// (`cxx:1679-1701`, the branch OCCT itself documents as the usual one).
    /// The walk returns the first edge starting or ending at `V` and answers
    /// its [`Self::uv_points`] first / last point, `Last` for a degenerated
    /// edge ("Ambiguity (natural)", `cxx:1698`).
    ///
    /// `Standard_NoSuchObject` (`cxx:1702`) is the `Err` of the port.
    pub fn parameters_of_vertex_on_face(
        vertex: &TopoShape,
        face: &Face,
    ) -> Result<GpPnt2d, String> {
        // `cxx:1663-1677`: the vertex's point representations.
        if let Ok(ts) = vertex.tshape.read() {
            if let Some(core) = ts.vertex_core() {
                if let Some(surf) = Self::face_surface(face) {
                    for pr in &core.points {
                        if pr.is_point_on_surface_of(&surf) {
                            if let Some(p2) = pr.parameter2() {
                                return Ok(GpPnt2d::new(pr.parameter(), p2));
                            }
                        }
                    }
                }
            }
        }
        for value in crate::topo_tools::edges_of(&face.0) {
            let (vf, vl) = crate::topo_tools::edge_vertices(&value);
            let at_first = vf
                .as_ref()
                .is_some_and(|v| crate::topo_tools::is_same(&v.0, vertex));
            let at_last = !at_first
                && vl
                    .as_ref()
                    .is_some_and(|v| crate::topo_tools::is_same(&v.0, vertex));
            if at_first || at_last {
                let (pf, pl) = Self::uv_points(&value, face);
                return Ok(if at_first { pf } else { pl });
            }
        }
        Err("BRep_Tool:: no parameters on surface".into())
    }
}
