use std::collections::{HashMap, VecDeque};

use crate::xcaf_prs::XcafPrsStyle;
use occt_core::gp::GpPnt;
use occt_core::poly::triangulation::Triangle;
use occt_topo::tgeometry::GeometryRegistry;

/// One face's mesh data as `RWMesh_FaceIterator::Next` obtains it
/// (`RWMesh_FaceIterator.cxx:83-101`): the face's own
/// `BRep_Tool::Triangulation`, else the triangulation the export driver meshed
/// into the face.
pub(super) struct FaceSource {
    vertices: Vec<GpPnt>,
    triangles: Vec<Triangle>,
    uv: Vec<occt_core::gp::GpPnt2d>,
    /// `Poly_Triangulation::Normal` (`cxx:59-67`): the stored per-node normals.
    normals: Option<Vec<GpPnt>>,
}

/// `BRep_Tool::Triangulation(theFace, theLocation)` first (`cxx:87`), then the
/// mesh-pass entry for that face. `None` means the face has neither, which
/// `RWMesh_FaceIterator::Next` skips together with a zero-triangle
/// triangulation (`cxx:88-93`).
pub(super) fn face_source(
    face: &occt_topo::shape::Face,
    by_key: &mut HashMap<usize, VecDeque<usize>>,
    tris: &[occt_topo::meshing::triangulator::FaceTriangulation],
) -> Option<FaceSource> {
    if let Some(tri) = occt_topo::brep_tool::BRepTool::face_triangulation(face) {
        if tri.triangles.is_empty() {
            return None;
        }
        return Some(FaceSource {
            vertices: tri.nodes.clone(),
            triangles: tri.triangles.clone(),
            uv: tri
                .uv_nodes
                .clone()
                .unwrap_or_default()
                .into_iter()
                .map(|(u, v)| occt_core::gp::GpPnt2d::new(u, v))
                .collect(),
            normals: tri.normals.clone(),
        });
    }
    let key = GeometryRegistry::shape_key(&face.0);
    let entry = by_key.get_mut(&key)?.pop_front()?;
    let ft = &tris[entry];
    Some(FaceSource {
        vertices: ft.vertices.clone(),
        triangles: ft.triangles.clone(),
        uv: ft.uv.clone(),
        normals: None,
    })
}

/// `RWMesh_FaceIterator` for one face (`hxx:85-152`): orientation-applied
/// triangles, `HasNormals` / `HasTexCoords` and the transformed normals.
///
/// `BRepLProp_SLProps(1, 1e-12)` is the iterator's surface-properties tool
/// (`RWMesh_FaceIterator.cxx:35,47`); its resolution is the CSLib tolerance
/// `normal()` hands to the normal estimate.
pub(super) const SL_PROPS_LIN_TOL: f64 = 1e-12;

pub(super) fn face_node(face: &occt_topo::shape::Face, source: FaceSource, style: XcafPrsStyle) -> crate::rwmesh::ObjFaceNode {
    // `RWMesh_FaceIterator::TriangleOriented` (`hxx:97-105`): a
    // `TopAbs_REVERSED` face occurrence swaps the last two nodes. The
    // `myIsMirrored` term needs a negative-determinant occurrence transform,
    // which this crate's shapes do not carry.
    let reversed = face.0.orientation().is_reversed();
    let triangles: Vec<Triangle> = source
        .triangles
        .iter()
        .map(|t| if reversed { Triangle::new(t.n0, t.n2, t.n1) } else { *t })
        .collect();

    // `RWMesh_FaceIterator::initFace` (`cxx:104-125`): `HasNormals()` is set
    // when the triangulation carries normals, or when it has UV nodes and the
    // face's surface (forward-oriented, unlocated, as in `cxx:115-121`) is
    // non-null. `HasTexCoords()` is `Poly_Triangulation::HasUVNodes`
    // (`hxx:112`).
    let has_tex_coords = !source.uv.is_empty() && source.uv.len() == source.vertices.len();
    let surface = if has_tex_coords {
        occt_topo::brep_tool::BRepTool::face_surface(face)
    } else {
        None
    };
    let has_normals = source.normals.is_some() || surface.is_some();

    // `RWMesh_FaceIterator::NormalTransformed(i)` (`hxx:113-126`): the stored
    // normal, transformed by the face location (identity here) and reversed for
    // a `TopAbs_REVERSED` face. Without stored normals the value comes from
    // `normal()` (`cxx:59-77`): `BRepLProp_SLProps::Normal` at the node's UV,
    // and `gp::DZ()` when `IsNormalDefined()` is false (a singular point such as
    // a cone apex). `BRepLProp_SLProps(1, 1e-12)` is the iterator's tool
    // (`cxx:35,47`), so the CSLib tolerance is `1e-12`.
    let normals: Vec<GpPnt> = if let Some(mut ns) = source.normals {
        if reversed {
            for n in &mut ns {
                *n = GpPnt::new(-n.x(), -n.y(), -n.z());
            }
        }
        ns
    } else if has_normals {
        let mut ns = Vec::with_capacity(source.vertices.len());
        for uv in &source.uv {
            // `gp_Dir aNormal(gp::DZ());`
            let mut dir = GpPnt::new(0.0, 0.0, 1.0);
            if let Some(surf) = surface.as_deref() {
                let (_, du, dv) = surf.d1(uv.x(), uv.y());
                if let (occt_core::cslib::CSLibDerivativeStatus::Done, Some(n)) =
                    occt_core::cslib::normal_d1(&du, &dv, SL_PROPS_LIN_TOL)
                {
                    dir = GpPnt::new(n.x(), n.y(), n.z());
                }
            }
            if reversed {
                dir = GpPnt::new(-dir.x(), -dir.y(), -dir.z());
            }
            ns.push(dir);
        }
        ns
    } else {
        Vec::new()
    };

    crate::rwmesh::ObjFaceNode {
        shape_key: GeometryRegistry::shape_key(&face.0),
        vertices: source.vertices,
        triangles,
        has_normals,
        has_tex_coords,
        normals,
        tex_coords: source.uv,
        style,
    }
}
