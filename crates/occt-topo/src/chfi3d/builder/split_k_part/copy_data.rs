//! The `SurfData` copy statics of `ChFi3d_Builder_SpKP.cxx` and the small
//! adapters `SplitKPart` needs around them.
//!
//! Sources: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_SpKP.cxx`
//! (OCCT 8.0.0): `CpInterf` (`:117-138`), `CpSD` (`:145-158`), the
//! `ChFi3d_cherche_element` / `ChFi3d_EdgeState` switcher setup
//! (`:1207-1225`, `:1344-1363`), and the `NCollection_Sequence` removals.

use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;
use std::sync::Arc;

use crate::abs::Orientation;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi_ds::{FaceInterference, SpineHandle, State, SurfData};
use crate::geom2d_hatch::Hatcher;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topopebrep::ds::{Curve, DataStructure, Surface};

/// `NCollection_Sequence::Remove(theIndex, theCount)`: `theIndex` 1-based, the
/// tail beyond the end clamped the way OCCT's `Remove(ilast + 1, Length())`
/// and `Remove(1, i - 1)` calls rely on.
pub(super) fn remove_count<T>(
    v: &mut Vec<T>,
    index: usize,
    count: usize,
) -> Result<(), String> {
    if count == 0 {
        return Ok(());
    }
    if index < 1 {
        return Err("ChFi3d_Builder::SplitKPart: sequence Remove out of range".into());
    }
    let last = (index + count - 1).min(v.len());
    if index <= v.len() {
        v.drain(index - 1..last);
    }
    Ok(())
}

/// `I->Init(); while (I->More()) { ... }` (`cxx:776-790`, `:811-825`):
/// `AddElement(Bc, Bc.Edge().Orientation())`, or `AddElement(*Gc,
/// TopAbs_FORWARD)` for an arc that is not edge-based, plus the `M` binding
/// `FillSD` reads back.
pub(super) fn add_elements(
    h: &mut Hatcher,
    m: &mut HashMap<i32, Option<Edge>>,
    tool: &Rc<RefCell<TopolTool>>,
) {
    let mut tool = tool.borrow_mut();
    tool.init();
    while tool.more() {
        let curve = tool.value().and_then(|v| v.curve()).cloned();
        if let Some(curve) = curve {
            let (ie, entry) = match tool.edge() {
                Some(e) => (h.add_element(curve, e.orientation()), Some(e.clone())),
                None => (h.add_element(curve, Orientation::Forward), None),
            };
            m.insert(ie, entry);
        }
        tool.next();
    }
}

/// `ChFi3d_cherche_element(bout, support, F1, threeE[0], boutemp)` for `F1` and
/// `F2`, with `threeE[2] = support` (`cxx:1209-1212`, `:1347-1350`).
pub(super) fn three_edges(
    bout: &Option<TopoShape>,
    support: &Edge,
    f1: &Option<Face>,
    f2: &Option<Face>,
) -> Result<[Edge; 3], String> {
    let (Some(bout), Some(f1), Some(f2)) = (bout.as_ref(), f1.as_ref(), f2.as_ref()) else {
        return Err("ChFi3d_Builder::SplitKPart: no corner vertex".into());
    };
    let v = Vertex::wrap(bout.clone())
        .ok_or_else(|| "ChFi3d_Builder::SplitKPart: not a vertex".to_string())?;
    let (e0, _) =
        crate::chfi3d::builder_helpers::vertex_edges::cherche_element(&v, support, f1)?;
    let (e1, _) =
        crate::chfi3d::builder_helpers::vertex_edges::cherche_element(&v, support, f2)?;
    Ok([e0, e1, support.clone()])
}

/// `Spine->FirstStatus()`.
pub(super) fn spine_first_status(spine: &SpineHandle) -> Result<State, String> {
    Ok(spine
        .read()
        .map_err(|_| "ChFi3d_Builder::SplitKPart: poisoned spine".to_string())?
        .first_status())
}

/// `Spine->LastStatus()`.
pub(super) fn spine_last_status(spine: &SpineHandle) -> Result<State, String> {
    Ok(spine
        .read()
        .map_err(|_| "ChFi3d_Builder::SplitKPart: poisoned spine".to_string())?
        .last_status())
}

/// `CpInterf(DStr, FI)` (`cxx:117-138`): shares the faces, copies the surface
/// and the curves. The port's `Geom_Curve` / `Geom2d_Curve` are immutable
/// behind `Arc`, so `Copy()` is an `Arc::clone`.
fn cp_interf(
    ds: &mut DataStructure,
    fi: &FaceInterference,
) -> Result<FaceInterference, String> {
    let mut new_f = fi.clone();
    let line_index = fi.line_index();
    if line_index == 0 {
        return Err("ChFi3d_Builder::SplitKPart: CpInterf on a null curve".into());
    }
    let tol = ds.curve(line_index).tolerance();
    let new_c = ds.curve(line_index).curve().cloned();
    let index = ds.add_curve(&mut Curve::with_curve(new_c, tol, false));
    new_f.set_line_index(index);
    if let Some(pc) = fi.p_curve_on_face() {
        new_f.set_p_curve_on_face(Some(Arc::clone(pc)));
    }
    if let Some(pc) = fi.p_curve_on_surf() {
        new_f.set_p_curve_on_surf(Some(Arc::clone(pc)));
    }
    Ok(new_f)
}

/// `CpSD(DStr, Data)` (`cxx:145-158`): a new `SurfData` sharing the faces and
/// the surface/curve indexes of the data structure, with fresh copies of the
/// surface and of the two interferences. The vertices and the 2d points are
/// *not* copied — `FillSD` fills them again for each domain.
pub(super) fn cp_sd(ds: &mut DataStructure, data: &SurfData) -> Result<SurfData, String> {
    let (new_s, tol) = {
        let tos = ds.surface(data.surf());
        (
            tos.surface()
                .cloned()
                .ok_or_else(|| "ChFi3d_Builder::SplitKPart: null surface".to_string())?,
            tos.tolerance(),
        )
    };
    let mut new_data = SurfData::new();
    new_data.set_surf(ds.add_surface(&Surface::with_tolerance(new_s, tol)));
    new_data.set_index_of_s1(data.index_of_s1());
    new_data.set_index_of_s2(data.index_of_s2());
    new_data.set_orientation(data.orientation());
    let i1 = cp_interf(ds, data.interference_on_s1())?;
    let i2 = cp_interf(ds, data.interference_on_s2())?;
    new_data.set_interference_on_s1(i1);
    new_data.set_interference_on_s2(i2);
    Ok(new_data)
}
