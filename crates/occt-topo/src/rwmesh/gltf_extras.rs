//! glTF `extras` → `TDataStd_NamedData`, the parser OCCT keeps in an unbounded
//! namespace inside `RWGltf_GltfJsonParser.cxx` (`RWGltf_ExtrasParser`,
//! `cxx:95-330`).
//!
//! A `extras` value must be a JSON object. Every member is flattened into a
//! dotted attribute name (`parent.child`, `cxx:194-197`) and stored by type:
//! integer, real, string, array of integers or array of reals. Nested objects
//! recurse; an empty array becomes an empty string, and an array of strings is
//! stored as one string with `;` separators (`cxx:263-312`). Anything else is
//! skipped.

use crate::tdata_std_named_data::TDataStdNamedData;

use super::gltf::JVal;

/// `RWGltf_ExtrasParser::ParseExtras` (`cxx:168-179`): `None` when there is no
/// `extras` member or when it holds nothing storable — the C++ side returns a
/// null handle unless `parseObject` initialized the result.
pub(super) fn parse_extras(value: Option<&JVal>) -> Option<TDataStdNamedData> {
    let value = value?;
    let mut out = TDataStdNamedData::default();
    if !parse_object(&mut out, value, "") {
        return None;
    }
    Some(out)
}

/// `RWGltf_ExtrasParser::parseObject` (`cxx:183-207`).
fn parse_object(out: &mut TDataStdNamedData, value: &JVal, prefix: &str) -> bool {
    let Some(items) = value.as_obj() else {
        return false;
    };
    let mut any = false;
    for (key, item) in items {
        let name = if prefix.is_empty() {
            key.clone()
        } else {
            format!("{prefix}.{key}")
        };
        // `cxx:199-201`: the first succeeded parser wins.
        let processed = parse_number(out, item, &name)
            || parse_string(out, item, &name)
            || parse_array(out, item, &name)
            || parse_object(out, item, &name);
        any |= processed;
    }
    any
}

/// `RWGltf_ExtrasParser::parseNumber` (`cxx:211-236`).
fn parse_number(out: &mut TDataStdNamedData, value: &JVal, name: &str) -> bool {
    if let Some(n) = value.as_int() {
        out.set_integer(name, n as i32);
        return true;
    }
    if let Some(n) = value.as_double() {
        out.set_real(name, n);
        return true;
    }
    false
}

/// `RWGltf_ExtrasParser::parseString` (`cxx:238-251`).
fn parse_string(out: &mut TDataStdNamedData, value: &JVal, name: &str) -> bool {
    if let Some(s) = value.as_str() {
        out.set_string(name, s);
        return true;
    }
    false
}

/// `RWGltf_ExtrasParser::parseArray` (`cxx:255-318`).
fn parse_array(out: &mut TDataStdNamedData, value: &JVal, name: &str) -> bool {
    let Some(items) = value.as_arr() else {
        return false;
    };
    if items.is_empty() {
        // `cxx:263-270`: an empty array is stored as an empty string.
        out.set_string(name, "");
        return true;
    }
    if items[0].as_int().is_some() {
        // `cxx:273-285`. Note `NCollection_HArray1<int>(0, theValue.Size())` has
        // one slot more than the JSON array; the last one keeps its default.
        let mut array = vec![0i32; items.len() + 1];
        for (i, item) in items.iter().enumerate() {
            array[i] = item.as_int().unwrap_or(0) as i32;
        }
        out.set_array_of_integers(name, array);
        return true;
    }
    if items[0].as_double().is_some() {
        // `cxx:286-297`, same one-slot-longer array as the integer case.
        let mut array = vec![0.0f64; items.len() + 1];
        for (i, item) in items.iter().enumerate() {
            array[i] = item.as_double().unwrap_or(0.0);
        }
        out.set_array_of_reals(name, array);
        return true;
    }
    if items[0].as_str().is_some() {
        // `cxx:298-311`: an array of strings becomes one ';'-separated string;
        // the separator also precedes the first element.
        let mut text = String::new();
        for item in items {
            text.push(';');
            text.push_str(item.as_str().unwrap_or(""));
        }
        out.set_string(name, &text);
        return true;
    }
    // `cxx:314-317`: unsupported array type, skipped with a warning (the port
    // has no message stream).
    true
}
