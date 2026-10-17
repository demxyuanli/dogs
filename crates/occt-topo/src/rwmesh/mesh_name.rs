//! `RWMesh_NameFormat` and `RWMesh::FormatName`.
//!
//! Source: `RWMesh_NameFormat.hxx:21-45`, `RWMesh.cxx:22-96`.
//!
//! `RWMesh::FormatName` picks the display name of a document node out of the
//! `TDataStd_Name` attributes of two labels — the *instance* label
//! (`theLabel`, a component of an assembly) and the *product* label
//! (`theRefLabel`, the shape the component refers to) — plus the label's OCAF
//! entry string for the `ProductAndInstanceAndOcaf` flavour.
//!
//! The port has no `TDF_Label`/`TDataStd_Name`; the names a label carries are
//! the `String` fields the document readers already resolved
//! (`ModelLabel::name`), and the entry is the label index, the same stand-in
//! `XCAFPrs_DocumentExplorer::DefineChildId` uses
//! (`xcaf_prs_document_explorer::define_child_id`). A name the label does not
//! carry is the empty string, exactly what `ReadNameAttribute` returns when
//! `FindAttribute` fails (`RWMesh.cxx:22-27`).

/// `RWMesh_NameFormat` (`RWMesh_NameFormat.hxx:21-45`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MeshNameFormat {
    /// `RWMesh_NameFormat_Empty`: no name.
    #[default]
    Empty,
    /// `RWMesh_NameFormat_Product`: the Product label's name.
    Product,
    /// `RWMesh_NameFormat_Instance`: the instance label's name.
    Instance,
    /// `RWMesh_NameFormat_InstanceOrProduct`: the instance name, else the
    /// Product name.
    InstanceOrProduct,
    /// `RWMesh_NameFormat_ProductOrInstance`: the Product name, else the
    /// instance name.
    ProductOrInstance,
    /// `RWMesh_NameFormat_ProductAndInstance`: `"Product [Instance]"`.
    ProductAndInstance,
    /// `RWMesh_NameFormat_ProductAndInstanceAndOcaf`:
    /// `"Product [Instance] [Entry]"`.
    ProductAndInstanceAndOcaf,
}

/// `RWMesh::ReadNameAttribute` (`RWMesh.cxx:22-27`): the label's
/// `TDataStd_Name` value, the empty string when the label carries none.
///
/// The port takes the already-resolved name string; see the module docs.
fn read_name_attribute(name: &str) -> &str {
    name
}

/// `RWMesh::FormatName` (`RWMesh.cxx:32-96`).
///
/// `instance_name` is `ReadNameAttribute(theLabel)`, `product_name` is
/// `ReadNameAttribute(theRefLabel)`, `entry` is `TDF_Tool::Entry(theLabel)`
/// (the label index in the port).
pub fn format_name(
    format: MeshNameFormat,
    instance_name: &str,
    product_name: &str,
    entry: &str,
) -> String {
    let inst = read_name_attribute(instance_name);
    let prod = read_name_attribute(product_name);
    match format {
        MeshNameFormat::Empty => String::new(),
        // `theRefLabel.FindAttribute(...)` (`cxx:39-46`).
        MeshNameFormat::Product => prod.to_string(),
        // `theLabel.FindAttribute(...)` (`cxx:47-53`).
        MeshNameFormat::Instance => inst.to_string(),
        // `cxx:54-63`: the instance name when non-empty, else the Product name.
        MeshNameFormat::InstanceOrProduct => {
            if !inst.is_empty() {
                inst.to_string()
            } else {
                prod.to_string()
            }
        }
        // `cxx:64-75`: the Product name when non-empty, else the instance name.
        MeshNameFormat::ProductOrInstance => {
            if !prod.is_empty() {
                prod.to_string()
            } else {
                inst.to_string()
            }
        }
        // `cxx:76-83`.
        MeshNameFormat::ProductAndInstance => {
            if !inst.is_empty() && prod != inst {
                format!("{prod} [{inst}]")
            } else if !prod.is_empty() {
                prod.to_string()
            } else {
                String::new()
            }
        }
        // `cxx:84-92`.
        MeshNameFormat::ProductAndInstanceAndOcaf => {
            if !inst.is_empty() && prod != inst {
                format!("{prod} [{inst}] [{entry}]")
            } else {
                format!("{prod} [{entry}]")
            }
        }
    }
}
