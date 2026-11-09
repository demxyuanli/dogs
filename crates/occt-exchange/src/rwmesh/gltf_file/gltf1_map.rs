//! glTF 1.0 keyed collections: rewritten into the glTF 2.0 array form once per
//! read (`RWGltf_GltfJsonParser` resolves either form through
//! `RWGltf_GltfJsonValue::FindChild`).
//!
//! The keyed form is an object from string id to item; the port converts the
//! collections the geometry path indexes into (`buffers`, `bufferViews`,
//! `accessors`, `meshes`, `nodes`, `scenes`), keeps the id lists for the glTF
//! 1.0 name fallback (`RWGltf_GltfJsonParser.cxx:2410-2418`), and resolves every
//! reference string to the array index it now addresses.

use super::gltf::JVal;
use super::prelude::*;

/// Looks up `key` in an object value, for mutation.
pub(super) fn member_mut<'a>(v: &'a mut JVal, key: &str) -> Option<&'a mut JVal> {
    match v {
        JVal::Obj(entries) => entries.iter_mut().find(|(k, _)| k == key).map(|(_, val)| val),
        _ => None,
    }
}

/// Every item of a collection, for mutation (array order preserved).
pub(super) fn items_mut(v: &mut JVal) -> Vec<&mut JVal> {
    match v {
        JVal::Arr(items) => items.iter_mut().collect(),
        JVal::Obj(entries) => entries.iter_mut().map(|(_, val)| val).collect(),
        _ => Vec::new(),
    }
}

/// glTF 1.0 stores every collection as an object keyed by a string id, glTF 2.0
/// as an array. `RWGltf_GltfJsonParser` resolves both forms through
/// `RWGltf_GltfJsonValue::FindChild` (array index or object key), so the port
/// rewrites the keyed form into the array form once and resolves the reference
/// strings with the id maps.
///
/// The collections converted here are the ones the geometry path indexes into:
/// `buffers`, `bufferViews`, `accessors`, `meshes`, `nodes`, `scenes`. Images,
/// textures and materials keep their keyed form because their parse path
/// already reads the ids (`gltfParseTexture`, `gltfParseMaterial`).
pub(super) fn normalize_gltf1(root: &mut JVal) -> (Vec<String>, Vec<String>) {
    let (buffer_map, _) = keyed_to_array(root, "buffers");
    let (view_map, _) = keyed_to_array(root, "bufferViews");
    let (accessor_map, _) = keyed_to_array(root, "accessors");
    let (mesh_map, mesh_ids) = keyed_to_array(root, "meshes");
    let (node_map, node_ids) = keyed_to_array(root, "nodes");
    let (scene_map, _) = keyed_to_array(root, "scenes");

    if let Some(v) = member_mut(root, "bufferViews") {
        for bv in items_mut(v) {
            resolve_member(bv, "buffer", &buffer_map);
        }
    }
    if let Some(v) = member_mut(root, "accessors") {
        for ac in items_mut(v) {
            resolve_member(ac, "bufferView", &view_map);
        }
    }
    if let Some(v) = member_mut(root, "meshes") {
        for mesh in items_mut(v) {
            let Some(prims) = member_mut(mesh, "primitives") else { continue };
            for prim in items_mut(prims) {
                if let Some(attrs) = member_mut(prim, "attributes") {
                    if let JVal::Obj(entries) = attrs {
                        for (_, val) in entries.iter_mut() {
                            resolve_value(val, &accessor_map);
                        }
                    }
                }
                resolve_member(prim, "indices", &accessor_map);
            }
        }
    }
    if let Some(v) = member_mut(root, "nodes") {
        if let JVal::Arr(items) = v {
            for nd in items.iter_mut() {
                // `mesh`/`meshes` point into the mesh collection, `children`
                // into the node collection.
                resolve_member(nd, "mesh", &mesh_map);
                resolve_member(nd, "meshes", &mesh_map);
                resolve_member(nd, "children", &node_map);
            }
        }
    }
    if let Some(v) = member_mut(root, "scenes") {
        if let JVal::Arr(items) = v {
            for sc in items.iter_mut() {
                resolve_member(sc, "nodes", &node_map);
            }
        }
    }
    // `scene` is a key in 1.0 and an index in 2.0.
    resolve_member(root, "scene", &scene_map);
    (node_ids, mesh_ids)
}

/// Convert an object collection to an array and return its id → index map plus
/// the id of every item in array order.
///
/// `RWGltf_GltfJsonValue::FindChild` addresses array items by index and object
/// items by key, so the port keeps the addressing form (array) and the original
/// ids side by side; the ids feed the glTF 1.0 name fallback
/// (`RWGltf_GltfJsonParser.cxx:2410-2418`).
pub(super) fn keyed_to_array(root: &mut JVal, key: &str) -> (HashMap<String, usize>, Vec<String>) {
    let Some(v) = member_mut(root, key) else { return (HashMap::new(), Vec::new()) };
    match v {
        JVal::Obj(entries) => {
            let mut map = HashMap::with_capacity(entries.len());
            let mut ids = Vec::with_capacity(entries.len());
            let mut arr = Vec::with_capacity(entries.len());
            for (i, (id, val)) in entries.drain(..).enumerate() {
                map.insert(id.clone(), i);
                ids.push(id);
                arr.push(val);
            }
            *v = JVal::Arr(arr);
            (map, ids)
        }
        // Already an array: the addressing key of an item is its index
        // (`getKeyString` on an array element answers the index as a string).
        JVal::Arr(items) => (
            items.iter().enumerate().map(|(i, _)| (i.to_string(), i)).collect(),
            (0..items.len()).map(|i| i.to_string()).collect(),
        ),
        _ => (HashMap::new(), Vec::new()),
    }
}

/// Resolve one member of `v` if it is a reference into `map`.
pub(super) fn resolve_member(v: &mut JVal, key: &str, map: &HashMap<String, usize>) {
    if let Some(slot) = member_mut(v, key) {
        resolve_value(slot, map);
    }
}

/// Replace a string id (or an array of them) by the matching index.
pub(super) fn resolve_value(v: &mut JVal, map: &HashMap<String, usize>) {
    match v {
        JVal::Str(s) => {
            if let Some(idx) = map.get(s.as_str()) {
                *v = JVal::Int(*idx as i64);
            }
        }
        JVal::Arr(items) => {
            for item in items.iter_mut() {
                resolve_value(item, map);
            }
        }
        _ => {}
    }
}
