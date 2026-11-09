use super::prelude::*;


/// `ChFiDS_FreeBound` of `BRepOffset_Analyse::Type(E).First()` (`cxx:563`).
pub(super) fn is_free_bound(kind: TypeOfConcavity) -> bool {
    kind == TypeOfConcavity::FreeBound
}


/// `OrientationOfEdgeInFace(E, F)` (`cxx:326-343`): the orientation of the
/// `IsSame` occurrence of `E` in `F`, `TopAbs_EXTERNAL` when there is none.
pub(super) fn orientation_of_edge_in_face(edge: &Edge, face: &Face) -> Orientation {
    for e in explorer_shapes(&face.0, ShapeType::Edge) {
        if e.same_tshape(&edge.0) {
            return e.orientation();
        }
    }
    Orientation::External
}


/// `gce_MakeCone(P1, P2, R1, R2)` (`gce_MakeCone.cxx:234-283`): the cone whose
/// axis runs from `P1` to `P2` with section radii `R1` / `R2`.
pub(super) fn gce_make_cone(p1: &GpPnt, p2: &GpPnt, r1: f64, r2: f64) -> Option<GpCone> {
    let dist = p1.distance(p2);
    if dist < REAL_SMALL {
        return None; // `gce_NullAxis`
    }
    if r1 < 0.0 || r2 < 0.0 {
        return None; // `gce_NegativeRadius`
    }
    let mut angle = ((r1 - r2) / dist).atan().abs();
    if (std::f64::consts::FRAC_PI_2 - angle).abs() < REAL_SMALL || angle.abs() < REAL_SMALL {
        return None; // `gce_NullAngle`
    }
    let d1 = GpDir::from_vec(&GpVec::from_pnts(p1, p2)).ok()?;
    let (x, y, z) = (d1.x(), d1.y(), d1.z());
    let d2 = if x.abs() > REAL_SMALL {
        GpDir::new(-y, x, 0.0).ok()?
    } else if y.abs() > REAL_SMALL {
        GpDir::new(0.0, -z, y).ok()?
    } else if z.abs() > REAL_SMALL {
        GpDir::new(z, 0.0, -x).ok()?
    } else {
        return None;
    };
    if r1 > r2 {
        angle = -angle;
    }
    let ax2 = GpAx2::new(*p1, d1, d2).ok()?;
    GpCone::new(ax2.to_ax3(), r1, angle).ok()
}


/// `BRepLib_MakeEdge(V1, V2)` (`BRepLib_MakeEdge.cxx:185-198`).
pub(super) fn make_edge_between(vertex1: &Vertex, vertex2: &Vertex) -> Option<Edge> {
    let p1 = BRepTool::vertex_point(vertex1);
    let p2 = BRepTool::vertex_point(vertex2);
    if p1.distance(&p2) <= REAL_SMALL {
        return None; // `BRepLib_LineThroughIdenticPoints`
    }
    Some(TopoBuilder::new().make_edge_segment_with_vertices(&p1, &p2, vertex1, vertex2))
}


/// A straight 2D line anchored at `p`.
pub(super) fn line2d_at(p: GpPnt2d, dir: GpDir2d) -> Geom2dLine {
    Geom2dLine::from_pnt_dir(p, dir)
}


/// `TopoDS_Shape::Oriented(Or)` / `Reversed()`: `Orientation(Compose(Or,
/// myOrient))`. The port's [`TopoShape::oriented`] *replaces* the orientation,
/// so every reference `Oriented` / `Reversed` call goes through this helper.
pub(super) fn oriented_by(shape: &TopoShape, o: Orientation) -> TopoShape {
    let mut out = shape.clone();
    out.set_orientation(Orientation::compose(o, shape.orientation()));
    out
}
