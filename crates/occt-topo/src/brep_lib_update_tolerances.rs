//! `BRepLib::UpdateTolerances` (`BRepLib.cxx:1744-1979`).
//!
//! `InternalUpdateTolerances` (`cxx:1744-1962`) harmonizes the tolerances with
//! the rule `Tolerance(VERTEX) >= Tolerance(EDGE) >= Tolerance(FACE)`. It works
//! through one shape-to-tolerance map and writes it back with `UpdShTol`
//! (`cxx:837-909`). The three phases are the reference's:
//!
//! 1. faces, only when `theVerifyFaceTolerance` (`cxx:1751-1823`): each face
//!    gets a tolerance derived from its surface type times the largest extent
//!    of its bounding box, capped at `0.99`;
//! 2. edges (`cxx:1826-1856`): `max` of the face tolerances, bound only when it
//!    is greater than the edge's own tolerance;
//! 3. vertices (`cxx:1858-1959`): `max` of the edge tolerances raised to the
//!    distance between the vertex and the edge representations at the vertex
//!    parameter.
//!
//! The port covers the `IsMutableInput = true` path, which is the one the
//! reference uses from [`crate::brep_lib_same_parameter`] and from
//! `BRepFill_Evolved` (`BRepFill_Evolved.cxx:540`, `UpdateTolerances(myShape,
//! false)`). With an empty `BRepTools_ReShape` and `IsMutableInput = true`,
//! `UpdShTol` takes `UseOldSh` for every entry (`cxx:852-857`) and writes the
//! tolerance on the shape itself (`cxx:877-899`), so the `EmptyCopied` /
//! `theReshaper.Replace` arm (`cxx:858-875`, `:903-906`) is not reached.
//! UNPORTED: the public overload taking a `BRepTools_ReShape`
//! (`BRepLib.cxx:1974-1979`, `IsMutableInput = false`) - it has no caller in
//! this port yet, and porting it means reproducing that copy arm.

use std::collections::{HashMap, HashSet};

use occt_core::bnd::BndBox;
use occt_core::precision::{epsilon, CONFUSION};

use crate::abs::ShapeType;
use crate::brep_bnd_lib;
use crate::brep_tool::BRepTool;
use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{
    faces_of, map_shapes_and_ancestors, map_shapes_and_unique_ancestors, shapes_of, shape_map_key,
    ShapeMapKey,
};

/// `theBigTol` of the vertex phase (`BRepLib.cxx:1859`).
const BIG_TOL: f64 = 1.0e10;

/// `BRepLib::UpdateTolerances(S, verifyFaceTolerance)` (`BRepLib.cxx:1964-1971`).
///
/// The reference builds a local `BRepTools_ReShape` and passes
/// `IsMutableInput = true`; the reshaper is only consulted through
/// `theReshaper.Value(aSh)`, which answers `aSh` for an untracked shape.
pub fn update_tolerances(shape: &TopoShape, verify_face_tolerance: bool) -> Result<(), String> {
    internal_update_tolerances(shape, verify_face_tolerance)
}

/// `InternalUpdateTolerances` (`BRepLib.cxx:1744-1962`) with
/// `IsMutableInput = true`.
fn internal_update_tolerances(
    shape: &TopoShape,
    verify_face_tolerance: bool,
) -> Result<(), String> {
    let reg = GeometryRegistry::global();

    // The reference's `NCollection_DataMap<TopoDS_Shape, double>`: insertion
    // order is kept because `UpdShTol` walks it in that order.
    let mut sh_to_tol: Vec<(TopoShape, f64)> = Vec::new();
    let mut index: HashMap<ShapeMapKey, usize> = HashMap::new();

    // `cxx:1750-1824`: harmonize the faces with their minimum value.
    if verify_face_tolerance {
        for (face, tol) in face_tolerances(shape)? {
            bind(&mut sh_to_tol, &mut index, face, tol);
        }
    }

    // `cxx:1826-1856`: process the edges against their face ancestors.
    let edge_parents = map_shapes_and_ancestors(shape, ShapeType::Edge, ShapeType::Face);
    for edge in shapes_of(shape, ShapeType::Edge) {
        let key = shape_map_key(&edge);
        let mut tol = 0.0_f64;
        if let Some(parents) = edge_parents.get(&key) {
            for parent in parents {
                let face = Face(parent.clone());
                let ftol = if verify_face_tolerance {
                    match bound(&sh_to_tol, &index, &shape_map_key(parent)) {
                        Some(t) => t,
                        None => BRepTool::face_tolerance(&face),
                    }
                } else {
                    BRepTool::face_tolerance(&face)
                };
                tol = tol.max(ftol);
            }
        }
        // Update can only increase a tolerance, so an edge already above its
        // faces is left alone (`cxx:1848-1855`).
        if tol > BRepTool::edge_tolerance(&Edge(edge.clone())) {
            bind(&mut sh_to_tol, &mut index, edge, tol);
        }
    }

    // `cxx:1858-1959`: process the vertices against their edge ancestors.
    let vertex_parents =
        map_shapes_and_unique_ancestors(shape, ShapeType::Vertex, ShapeType::Edge, false);
    let mut initialized: HashSet<usize> = HashSet::new();
    for vertex in shapes_of(shape, ShapeType::Vertex) {
        let key = shape_map_key(&vertex);
        let point = BRepTool::vertex_point_world(&Vertex(vertex.clone()));
        let mut tol = 0.0_f64;
        let mut max_dist_sq = 0.0_f64;
        if let Some(parents) = vertex_parents.get(&key) {
            for parent in parents {
                let edge = Edge(parent.clone());
                let etol = match bound(&sh_to_tol, &index, &shape_map_key(parent)) {
                    Some(t) => t,
                    None => BRepTool::edge_tolerance(&edge),
                };
                tol = tol.max(etol);
                if tol > BIG_TOL {
                    continue;
                }
                if !reg.same_range(parent) {
                    continue;
                }
                let par = BRepTool::parameter_of_vertex(parent, &vertex);
                // `L = Eloc * cr->Location()` (`cxx:1888-1892`). The registry
                // stores no per-representation location (see
                // `brep_lib_inner_tolerances`), so the edge location is the
                // whole of `L` for every representation.
                let l = parent.location().transformation();
                if let Some(curve) = reg.edge_curve(parent) {
                    // `cxx:1896-1906`: the non-degenerated 3D curve.
                    let p3d = curve.value(par).transformed(&l);
                    max_dist_sq = max_dist_sq.max(p3d.square_distance(&point));
                }
                // `cxx:1908-1937`: every CurveOnSurface representation, both
                // pcurves of a closed surface included.
                for (surface, pcurve, _, _) in reg.edge_pcurve_reps(parent) {
                    let p2d = pcurve.value(par);
                    let p3d = surface.value(p2d.x(), p2d.y()).transformed(&l);
                    max_dist_sq = max_dist_sq.max(p3d.square_distance(&point));
                }
            }
        }
        tol = tol.max(max_dist_sq.sqrt());
        tol += 2.0 * epsilon(tol);

        let v_tol = BRepTool::vertex_tolerance(&Vertex(vertex.clone()));
        let an_upd_tol = tol > v_tol;
        let mut to_add = false;
        if verify_face_tolerance {
            // `cxx:1945-1956`: `Initialized.Add(aTV) && aVTol != tol`. The key
            // is the shared `BRep_TVertex`, not the located occurrence.
            to_add = initialized.insert(GeometryRegistry::shape_key(&vertex)) && v_tol != tol;
        }
        if an_upd_tol || to_add {
            bind(&mut sh_to_tol, &mut index, vertex, tol);
        }
    }

    // `cxx:1961`: `UpdShTol(aShToTol, IsMutableInput, theReshaper, true)`.
    apply_tolerances(&sh_to_tol);
    Ok(())
}

/// `cxx:1751-1823`: the per-face minimum tolerance.
fn face_tolerances(shape: &TopoShape) -> Result<Vec<(TopoShape, f64)>, String> {
    let mut out = Vec::new();
    for face in faces_of(shape) {
        // `S = BRep_Tool::Surface(curf, l); if (!S.IsNull())` (`cxx:1765-1767`).
        let Some(surface) = BRepTool::face_surface(&face) else {
            continue;
        };
        // `BRepBndLib::Add(curf, aB)` (`cxx:1768`): the `useTriangulation`
        // default is `true` (`BRepBndLib.hxx:48`); the port stores no
        // triangulation, so the geometry arms are the whole path.
        let mut b = BndBox::new();
        brep_bnd_lib::add(&face.0, &mut b, true);
        // `GeomAdaptor_Surface AS(S)` (`cxx:1772-1790`); `surface_kind` answers
        // the type of the basis surface for a trimmed surface, as
        // `GeomAdaptor_Surface::GetType` does.
        let mut tol = match surface_kind(surface.as_ref()) {
            GeomAbsSurfaceType::Plane
            | GeomAbsSurfaceType::Cylinder
            | GeomAbsSurfaceType::Cone => CONFUSION,
            GeomAbsSurfaceType::Sphere | GeomAbsSurfaceType::Torus => CONFUSION * 2.0,
            _ => CONFUSION * 4.0,
        };
        if !b.is_whole() {
            // `aB.Get(...)`: a void box throws `Standard_ConstructionError`
            // (`Bnd_Box.cxx:214-217`), so the reference never falls through to
            // the extent arithmetic with an empty box.
            let (xmin, xmax, ymin, ymax, zmin, zmax) =
                b.get().ok_or_else(|| "Bnd_Box is void".to_string())?;
            let mut d_max = 1.0;
            if !b.is_open_xmin() && !b.is_open_xmax() {
                d_max = xmax - xmin;
            }
            // `cxx:1800-1808` reuses `aYmin` / `aZmin` for the extents.
            let mut ext_y = ymin;
            if !b.is_open_ymin() && !b.is_open_ymax() {
                ext_y = ymax - ymin;
            }
            let mut ext_z = zmin;
            if !b.is_open_zmin() && !b.is_open_zmax() {
                ext_z = zmax - zmin;
            }
            if ext_y > d_max {
                d_max = ext_y;
            }
            if ext_z > d_max {
                d_max = ext_z;
            }
            tol *= d_max;
            // Do not process tolerances > 1 (`cxx:1817-1820`).
            if tol > 1.0 {
                tol = 0.99;
            }
        }
        out.push((face.0, tol));
    }
    Ok(out)
}

/// `NCollection_DataMap::Bind` on the shape-to-tolerance map.
fn bind(
    map: &mut Vec<(TopoShape, f64)>,
    index: &mut HashMap<ShapeMapKey, usize>,
    shape: TopoShape,
    tol: f64,
) {
    let key = shape_map_key(&shape);
    match index.get(&key) {
        Some(&i) => map[i].1 = tol,
        None => {
            index.insert(key, map.len());
            map.push((shape, tol));
        }
    }
}

/// `NCollection_DataMap::Seek` on the shape-to-tolerance map.
fn bound(map: &[(TopoShape, f64)], index: &HashMap<ShapeMapKey, usize>, key: &ShapeMapKey) -> Option<f64> {
    index.get(key).map(|&i| map[i].1)
}

/// `UpdShTol` (`BRepLib.cxx:877-899`) for the `UseOldSh` arm: the tolerance is
/// written on the shape itself.
fn apply_tolerances(sh_to_tol: &[(TopoShape, f64)]) {
    let reg = GeometryRegistry::global();
    for (shape, tol) in sh_to_tol {
        match shape.shape_type() {
            // `aB.UpdateFace(aNsh, aTol)` (`cxx:879-881`): `BRep_Builder::
            // UpdateFace` writes the value (`BRep_Builder.cxx:597-606`).
            ShapeType::Face => reg.set_face_tolerance(shape, *tol),
            // `aB.UpdateEdge(aNsh, aTol)` (`cxx:883-885`): `BRep_TEdge::
            // UpdateTolerance` only raises.
            ShapeType::Edge => reg.raise_edge_tolerance(shape, *tol),
            // `theVForceUpdate` is `true` (`cxx:1961`), so the vertex
            // tolerance is set (`aTV->Tolerance(aTol)`, `cxx:888-892`).
            ShapeType::Vertex => set_vertex_tolerance(shape, *tol),
            _ => {}
        }
    }
}

/// `BRep_TVertex::Tolerance(tol)` (`BRep_TVertex.lxx:24-27`): direct write.
fn set_vertex_tolerance(v: &TopoShape, tol: f64) {
    let mut ts = v.tshape.write().expect("poisoned TShape lock");
    ts.vertex_core_mut().tolerance = tol;
}
