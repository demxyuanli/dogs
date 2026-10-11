//! `IGESSolid` own-parameter readers for the B-rep family:
//! 502 `VertexList`, 504 `EdgeList`, 508 `Loop`, 510 `Face`, 514 `Shell`, and
//! 186 `ManifoldSolidBRep`.
//!
//! Each function mirrors the matching `IGESSolid_Tool<Class>::ReadOwnParams`.
//! OCCT reads the referenced entities with a `STANDARD_TYPE(...)` argument,
//! which only reports `IGESData_TypeError` and leaves the slot null; the
//! directory is not available while a single record is parsed, so the port
//! keeps the pointer and defers that check to `brep_entity` (see
//! `resolve_brep_ref`).

use super::iges_data_iges_model::EntityData;
use super::iges_data_param_reader::ParamReader;

/// One entry of `IGESSolid_EdgeList` (`IGESSolid_ToolEdgeList.cxx:74-200`).
#[derive(Clone, Debug)]
pub(in crate::iges_read) struct SolidEdgeItem {
    pub curvature: Option<usize>,
    pub start_vertex_list: Option<usize>,
    pub start_vertex_index: Option<i32>,
    pub end_vertex_list: Option<usize>,
    pub end_vertex_index: Option<i32>,
}

/// One entry of `IGESSolid_Loop` (`IGESSolid_ToolLoop.cxx:90-209`).
#[derive(Clone, Debug)]
pub(in crate::iges_read) struct SolidLoopItem {
    pub edge_type: Option<i32>,
    pub edge: Option<usize>,
    pub list_index: Option<i32>,
    pub orientation: bool,
    pub param_curves: Vec<Option<usize>>,
    /// `IGESSolid_Loop::IsoparametricFlag(i, j)` is stored by `Init` and is
    /// not read by `IGESToBRep_BRepEntity::TransferLoop`.
    #[allow(dead_code)]
    pub isoparametric_flags: Vec<bool>,
}

/// `IGESToBRep::IsBRepEntity` (`IGESToBRep.cxx:318-350`).
pub(in crate::iges_read) fn is_brep_entity(type_number: i32) -> bool {
    matches!(type_number, 510 | 514 | 186 | 502 | 504 | 508)
}

/// 502 `IGESSolid_ToolVertexList::ReadOwnParams` (`cxx:46-91`).
pub(in crate::iges_read) fn parse_vertex_list(
    pr: &mut ParamReader<'_>,
) -> Result<EntityData, String> {
    let count = pr.read_integer()?;
    let mut vertices = Vec::new();
    if count > 0 {
        for _ in 0..count {
            // `ReadXYZ(CurrentList(1, 3))` (`cxx:76`). A malformed triple leaves
            // the slot at its default rather than aborting the list.
            match pr.read_xyz() {
                Ok(p) => vertices.push(p),
                Err(_) => vertices.push([0.0, 0.0, 0.0]),
            }
        }
    } else {
        return Err("XSTEP_182: IGESSolid_ToolVertexList.cxx:85-88".into());
    }
    Ok(EntityData::VertexList { vertices })
}

/// 504 `IGESSolid_ToolEdgeList::ReadOwnParams` (`cxx:47-297`).
pub(in crate::iges_read) fn parse_edge_list(pr: &mut ParamReader<'_>) -> Result<EntityData, String> {
    let length = pr.read_integer()?;
    if length <= 0 {
        return Err("XSTEP_184: IGESSolid_ToolEdgeList.cxx:64-68".into());
    }
    let mut edges = Vec::with_capacity(length as usize);
    for _ in 0..length {
        // Model space curve (`cxx:81-107`).
        let mut curvature = None;
        if pr.defined_else_skip() {
            curvature = pr.read_entity(true)?;
        }
        // Start vertex list (`cxx:113-147`), then the start vertex index.
        let mut start_vertex_list = None;
        if pr.defined_else_skip() {
            start_vertex_list = pr.read_entity(true)?;
        }
        let start_vertex_index = pr.read_integer().ok();
        // End vertex list (`cxx:165-200`), then the end vertex index.
        let mut end_vertex_list = None;
        if pr.defined_else_skip() {
            end_vertex_list = pr.read_entity(true)?;
        }
        let end_vertex_index = pr.read_integer().ok();
        edges.push(SolidEdgeItem {
            curvature,
            start_vertex_list,
            start_vertex_index,
            end_vertex_list,
            end_vertex_index,
        });
    }
    Ok(EntityData::EdgeList { edges })
}

/// 508 `IGESSolid_ToolLoop::ReadOwnParams` (`cxx:53-232`).
pub(in crate::iges_read) fn parse_loop(pr: &mut ParamReader<'_>) -> Result<EntityData, String> {
    let nbedges = pr.read_integer()?;
    if nbedges <= 0 {
        // `cxx:216-223`: an empty loop fails the tool and `Init` is not called.
        return Err("XSTEP_184: IGESSolid_ToolLoop.cxx:216-223".into());
    }
    let mut edges = Vec::with_capacity(nbedges as usize);
    for _ in 0..nbedges {
        let edge_type = pr.read_integer().ok();
        let mut edge = None;
        if pr.defined_else_skip() {
            edge = pr.read_entity(true)?;
        }
        // `cxx:127-135`: the entity must be a 502 or a 504. A mismatch fails the
        // tool and leaves the slot null; the type check needs the directory, so
        // it is deferred to `resolve_brep_ref`.
        let list_index = pr.read_integer().ok();
        let orientation = pr.read_boolean()?;
        let nb_param = pr.read_integer().unwrap_or(0);
        let mut param_curves = Vec::new();
        let mut isoparametric_flags = Vec::new();
        if nb_param > 0 {
            for _ in 0..nb_param {
                isoparametric_flags.push(pr.read_boolean()?);
                let mut curve = None;
                if pr.defined_else_skip() {
                    curve = pr.read_entity(true)?;
                }
                param_curves.push(curve);
            }
        }
        edges.push(SolidLoopItem {
            edge_type,
            edge,
            list_index,
            orientation,
            param_curves,
            isoparametric_flags,
        });
    }
    Ok(EntityData::SolidLoop { edges })
}

/// 510 `IGESSolid_ToolFace::ReadOwnParams` (`cxx:47-...`).
pub(in crate::iges_read) fn parse_face(pr: &mut ParamReader<'_>) -> Result<EntityData, String> {
    let mut surface = None;
    if pr.defined_else_skip() {
        surface = pr.read_entity(true)?;
    }
    let nbloops = pr.read_integer()?;
    if nbloops <= 0 {
        return Err("XSTEP_197: IGESSolid_ToolFace.cxx:97-101".into());
    }
    // `cxx:104`: the outer loop flag is read but only kept by `Init`; the
    // transfer never reads it back.
    let outer_loop = pr.read_boolean()?;
    let mut loops = Vec::with_capacity(nbloops as usize);
    for _ in 0..nbloops {
        let mut child = None;
        if pr.defined_else_skip() {
            child = pr.read_entity(true)?;
        }
        loops.push(child);
    }
    Ok(EntityData::SolidFace {
        surface,
        loops,
        outer_loop,
    })
}

/// 514 `IGESSolid_ToolShell::ReadOwnParams` (`cxx:47-124`).
pub(in crate::iges_read) fn parse_shell(pr: &mut ParamReader<'_>) -> Result<EntityData, String> {
    let nbfaces = pr.read_integer()?;
    if nbfaces <= 0 {
        return Err("XSTEP_200: IGESSolid_ToolShell.cxx:114-118".into());
    }
    let mut faces = Vec::with_capacity(nbfaces as usize);
    let mut orientations = Vec::with_capacity(nbfaces as usize);
    for _ in 0..nbfaces {
        let mut face = None;
        if pr.defined_else_skip() {
            face = pr.read_entity(true)?;
        }
        faces.push(face);
        // `cxx:108-112`: a failed read leaves the orientation slot at 0.
        orientations.push(pr.read_boolean().unwrap_or(false));
    }
    Ok(EntityData::SolidShell {
        faces,
        orientations,
    })
}

/// 186 `IGESSolid_ToolManifoldSolid::ReadOwnParams` (`cxx:43-149`).
pub(in crate::iges_read) fn parse_manifold_solid(
    pr: &mut ParamReader<'_>,
) -> Result<EntityData, String> {
    let mut shell = None;
    if pr.defined_else_skip() {
        shell = pr.read_entity(true)?;
    }
    let oriented = pr.read_boolean().unwrap_or(false);
    let nbshells = pr.read_integer()?;
    if nbshells < 0 {
        return Err("XSTEP_181: IGESSolid_ToolManifoldSolid.cxx:85-88".into());
    }
    let mut void_shells = Vec::new();
    let mut void_flags = Vec::new();
    if nbshells > 0 {
        for _ in 0..nbshells {
            let mut void = None;
            if pr.defined_else_skip() {
                void = pr.read_entity(true)?;
            }
            void_shells.push(void);
            // `cxx:140-144`: stored for `VoidOrientationFlag`, which
            // `TransferManifoldSolid` does not read (it reuses the single
            // `isoriented` flag for every shell).
            void_flags.push(pr.read_boolean().unwrap_or(false));
        }
    }
    Ok(EntityData::ManifoldSolid {
        shell,
        oriented,
        void_shells,
        void_flags,
    })
}
