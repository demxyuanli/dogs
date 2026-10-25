//! `ChFi3d_SolidIndex`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:2300-2321`
//! (OCCT 8.0.0).

use crate::chfi_ds::ChFiDSMap;
use crate::chfi_ds::SpineHandle;
use crate::shape::TopoShape;
use crate::topopebrep::ds::HDataStructure;

/// `ChFi3d_SolidIndex(sp, DStr, MapESo, MapESh)` (`cxx:2300-2321`).
///
/// The index in `DStr` of the solid the first edge of `sp` belongs to: the
/// edge's first shell-or-solid from `MapESo`, or — when that map has no entry —
/// from `MapESh`.
///
/// `Err` carries the `Standard_Failure("SolidIndex : Spine incomplete")`
/// (`cxx:2305-2308`) and the `MapESo(edref).First()` / `MapESh(edref).First()`
/// reads of an empty list (`cxx:2311-2318`), which OCCT raises as
/// `Standard_OutOfRange`.
pub fn ch_fi3d_solid_index(
    sp: &SpineHandle,
    ds: &mut HDataStructure,
    map_eso: &ChFiDSMap,
    map_esh: &ChFiDSMap,
) -> Result<i32, String> {
    // `if (sp.IsNull() || sp->NbEdges() == 0) throw ...` (`cxx:2305-2308`).
    let edref: TopoShape = {
        let guard = sp
            .read()
            .map_err(|_| "ChFi3d_SolidIndex: poisoned spine".to_string())?;
        if guard.nb_edges() == 0 {
            return Err("SolidIndex : Spine incomplete".into());
        }
        // `TopoDS_Shape edref = sp->Edges(1);` (`cxx:2309`).
        guard
            .edges(1)
            .cloned()
            .ok_or_else(|| "SolidIndex : Spine incomplete".to_string())?
    };

    // `if (!MapESo(edref).IsEmpty()) shellousolid = MapESo(edref).First();
    // else shellousolid = MapESh(edref).First();` (`cxx:2311-2318`).
    let from_eso = map_eso.find_from_key(&edref);
    let shell_or_solid = if !from_eso.is_empty() {
        from_eso[0].clone()
    } else {
        let from_esh = map_esh.find_from_key(&edref);
        if from_esh.is_empty() {
            return Err("ChFi3d_SolidIndex: Standard_OutOfRange (empty map list)".into());
        }
        from_esh[0].clone()
    };

    // `const int solidindex = DStr.AddShape(shellousolid);` (`cxx:2319-2320`).
    Ok(ds.change_ds().add_shape(&shell_or_solid))
}
