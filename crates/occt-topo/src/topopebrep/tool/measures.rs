//! Shape-level measures used by the interference reducers.
//!
//! Source: `TKBool/TopOpeBRepTool/TopOpeBRepTool_TOPOLOGY.cxx:1148-1214`.

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topexp::Explorer;
use crate::topopebrep::tool::shape_tool::tolerance_of;

/// `FUN_tool_maxtol(S, typ, maxtol)` (`TopOpeBRepTool_TOPOLOGY.cxx:1148-1181`).
///
/// `maxtol` is an in/out value: the maximum of the incoming value and the
/// tolerances of the sub-shapes of type `typ`. The return value is whether any
/// such sub-shape exists (`ex.More()` before the loop).
pub fn fun_tool_maxtol_into(s: &TopoShape, typ: ShapeType, maxtol: &mut f64) -> bool {
    let mut ex = Explorer::new(s, typ);
    let has_such_type = ex.more();
    while ex.more() {
        let tol = tolerance_of(ex.current()).unwrap_or(0.0);
        if tol > *maxtol {
            *maxtol = tol;
        }
        ex.next();
    }
    has_such_type
}

/// `FUN_tool_maxtol(S)` (`TopOpeBRepTool_TOPOLOGY.cxx:1183-1191`): the maximum
/// over faces, edges and vertices, starting from `0.`.
pub fn fun_tool_maxtol(s: &TopoShape) -> f64 {
    let mut maxtol = 0.0;
    fun_tool_maxtol_into(s, ShapeType::Face, &mut maxtol);
    fun_tool_maxtol_into(s, ShapeType::Edge, &mut maxtol);
    fun_tool_maxtol_into(s, ShapeType::Vertex, &mut maxtol);
    maxtol
}

/// `FUN_tool_nbshapes(S, typ)` (`TopOpeBRepTool_TOPOLOGY.cxx:1194-1203`).
pub fn fun_tool_nbshapes(s: &TopoShape, typ: ShapeType) -> usize {
    let mut ex = Explorer::new(s, typ);
    let mut i = 0;
    while ex.more() {
        i += 1;
        ex.next();
    }
    i
}
