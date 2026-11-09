//! File-local helpers of `BRepOffset_MakeOffset::MakeOffsetShape` and
//! `Intersection3D` (`BRepOffset_MakeOffset.cxx`), ported with their reference
//! bodies:
//!
//! * `RemoveCorks(S, Faces)` (`:712-738`): rebuild `S` as a compound without
//!   the cap faces, and move every cap to the tail of `Faces` ("to reset it
//!   with proper orientation").
//! * `IsConnectedShell(S)` (`:742-750`): true when `BRepTools_Quilt` yields at
//!   most one shell.
//! * `MakeList(OffsetFaces, myInitOffsetFace, myFaces)` (`:754-775`): collect
//!   the offset-face images of the roots that are not caps.
//! * `EvalMax(S, Tol)` (`:779-791`): raise `Tol` to the largest vertex
//!   tolerance of `S`.
//! * `IsSolid(S)` (`:5637-5641`) and `AppendToList(theList, theShape)`
//!   (`:5647-5659`).
//!
//! Port guards, each replacing an OCCT exception or a null handle:
//!
//! * `IsConnectedShell` (`:746`) feeds `BRepTools_Quilt::Add` a shape without
//!   shells in some inputs; the port answers `true` for an empty shell list,
//!   which is what the reference's `Explo.Next(); !Explo.More()` evaluates to.
//! * `EvalMax` (`:786`) casts every explored shape to `TopoDS_Vertex`; the
//!   port skips a non-vertex the explorer cannot produce.

use crate::abs::ShapeType;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_sweep::Quilt;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{TopoShape, Vertex};
use crate::topo_tools::{explorer_shapes, is_same};

use crate::brep_offset::inter3d::ShapeIndexedMap;

/// `RemoveCorks(S, Faces)` (`cxx:712-738`).
///
/// `S` becomes a compound of the faces that are not caps; every cap stays in
/// `Faces`, moved to the end by `RemoveKey` + `Add` (`:732-733`).
pub fn remove_corks(s: &mut TopoShape, faces: &mut ShapeIndexedMap) {
    let builder = TopoBuilder::new();
    let mut children: Vec<TopoShape> = Vec::new();
    for cork in explorer_shapes(s, ShapeType::Face) {
        if !faces.contains(&cork) {
            // `B.Add(SS, Cork)` (`:729`).
            children.push(cork);
        } else {
            // `Faces.RemoveKey(Cork); Faces.Add(Cork)` (`:732-733`).
            faces.remove_key(&cork);
            faces.add(&cork);
        }
    }
    *s = builder.make_compound_of(&children).0;
}

/// `IsConnectedShell(S)` (`cxx:742-750`): `true` when the quilted shape has at
/// most one shell (`Explo.Next(); return !Explo.More()`).
pub fn is_connected_shell(s: &TopoShape) -> bool {
    let mut glue = Quilt::new();
    glue.add(s);
    let shells = glue.shells();
    explorer_shapes(&shells, ShapeType::Shell).len() <= 1
}

/// `MakeList(OffsetFaces, myInitOffsetFace, myFaces)` (`cxx:754-775`): the
/// images of every `myInitOffsetFace` root that is not a cap.
pub fn make_list(
    offset_faces: &mut Vec<TopoShape>,
    init_offset_face: &BRepAlgoImage,
    faces: &ShapeIndexedMap,
) {
    for root in init_offset_face.roots() {
        if faces.contains(root) {
            continue;
        }
        if init_offset_face.has_image(root) {
            offset_faces.extend(init_offset_face.image(root));
        }
    }
}

/// `EvalMax(S, Tol)` (`cxx:779-791`): `Tol = max(Tol, max vertex tolerance)`.
pub fn eval_max(s: &TopoShape, tol: &mut f64) {
    for shape in explorer_shapes(s, ShapeType::Vertex) {
        let Some(v) = Vertex::wrap(shape) else {
            continue;
        };
        let tol_v = BRepTool::vertex_tolerance(&v);
        if tol_v > *tol {
            *tol = tol_v;
        }
    }
}

/// `IsSolid(S)` (`cxx:5637-5641`): `true` when `S` contains a solid.
pub fn is_solid(s: &TopoShape) -> bool {
    !explorer_shapes(s, ShapeType::Solid).is_empty()
}

/// `AppendToList(theList, theShape)` (`cxx:5647-5659`): append unless an
/// `IsSame` entry is already there.
pub fn append_to_list(list: &mut Vec<TopoShape>, shape: &TopoShape) {
    if list.iter().any(|s| is_same(s, shape)) {
        return;
    }
    list.push(shape.clone());
}

/// `static bool GetSubShapes(const TopoDS_Shape& theShape, const
/// TopAbs_ShapeEnum theSSType, TopoDS_Shape& theResult)` (`cxx:5386-5406`).
///
/// Answers the compound of every `theSSType` sub-shape, or `None` where the
/// reference returns `false` without touching `theResult` (`:5391-5394`).
pub fn get_sub_shapes(the_shape: &TopoShape, the_ss_type: ShapeType) -> Option<TopoShape> {
    let children = explorer_shapes(the_shape, the_ss_type);
    if children.is_empty() {
        return None;
    }
    let builder = TopoBuilder::new();
    Some(builder.make_compound_of(&children).0)
}
