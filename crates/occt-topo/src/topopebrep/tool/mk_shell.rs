//! `TopOpeBRepTool_TOOL` list helpers used by the shell regularization.
//!
//! Source (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRepTool`):
//! - `TopOpeBRepTool_TOOL::MkShell(lF, She)` (`TopOpeBRepTool_TOOL.cxx:1819-1827`)
//! - `TopOpeBRepTool_TOOL::Remove(loS, toremove)` (`TOOL.cxx:1831-1848`)

use crate::builder::TopoBuilder;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::is_equal;

/// `TopOpeBRepTool_TOOL::MkShell(lF, She)` (`cxx:1819-1827`): `BB.MakeShell`
/// followed by one `BB.Add` per face, in list order.
pub fn mk_shell(lf: &[TopoShape]) -> TopoShape {
    let faces: Vec<Face> = lf.iter().map(|f| Face(f.clone())).collect();
    TopoBuilder::new().make_shell(&faces).0
}

/// `TopOpeBRepTool_TOOL::Remove(loS, toremove)` (`cxx:1831-1848`): drop every
/// occurrence equal to `to_remove`; answers whether one was found.
pub fn remove_from_list(lo_s: &mut Vec<TopoShape>, to_remove: &TopoShape) -> bool {
    let mut found = false;
    lo_s.retain(|s| {
        if is_equal(s, to_remove) {
            found = true;
            false
        } else {
            true
        }
    });
    found
}
