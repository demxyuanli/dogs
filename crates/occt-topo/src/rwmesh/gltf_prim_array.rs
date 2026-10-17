//! Deferred primitive arrays: the glTF reader's two-phase data loading.
//!
//! Source: `RWGltf_GltfArrayType.hxx`, `RWGltf_GltfAccessor.hxx`,
//! `RWGltf_GltfPrimArrayData.hxx`, `RWGltf_GltfLatePrimitiveArray.hxx/.cxx`,
//! and `RWMesh_TriangulationSource.hxx/.cxx`.
//!
//! The parser never decodes primitive bytes. `RWGltf_GltfJsonParser::gltfParseAccessor`
//! (`cxx:2039-2198`) only validates the accessor, then `gltfParseBufferView`
//! (`cxx:2206-2270`) and `gltfParseBuffer` (`cxx:2278-2367`) append one
//! `RWGltf_GltfPrimArrayData` per attribute: the accessor description, where the
//! bytes start, and where they live (an in-memory buffer decoded from a base64
//! `data:` URI, or a file the reader re-opens later).
//!
//! `RWGltf_GltfLatePrimitiveArray` collects those elements and is put on the
//! proxy face in place of a `Poly_Triangulation`
//! (`RWGltf_GltfJsonParser.cxx:1879-1885`, `:2000-2007`). `RWGltf_CafReader::readLateData`
//! (`RWGltf_CafReader.cxx:367-404`) then walks the face list and fills the
//! triangulation. Only triangulations defer: lines and points are filled at
//! parse time (`cxx:1960-1967`) because their shape needs the actual points.
//!
//! `RWMesh_TriangulationSource` is the `Poly_Triangulation` subclass that
//! carries the deferred counts and the edge array; this port has no
//! subclassing, so [`GltfLatePrimitiveArray`] holds the same state, and
//! [`GltfResolvedPrim`] is what a reader produces as the destination
//! triangulation (including the edge array a `Poly_Triangulation` itself does
//! not have).

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use occt_core::poly::Triangulation;

use super::gltf::RwGltfPrimitiveMode;

/// `RWGltf_GltfArrayType` (`RWGltf_GltfArrayType.hxx:23-35`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GltfArrayType {
    /// `RWGltf_GltfArrayType_UNKNOWN`.
    Unknown,
    /// `RWGltf_GltfArrayType_Indices`.
    Indices,
    /// `RWGltf_GltfArrayType_Position`.
    Position,
    /// `RWGltf_GltfArrayType_Normal`.
    Normal,
    /// `RWGltf_GltfArrayType_Color`.
    Color,
    /// `RWGltf_GltfArrayType_TCoord0`.
    TCoord0,
    /// `RWGltf_GltfArrayType_TCoord1`.
    TCoord1,
    /// `RWGltf_GltfArrayType_Joint`.
    Joint,
    /// `RWGltf_GltfArrayType_Weight`.
    Weight,
}

/// `RWGltf_GltfParseAttribType` (`RWGltf_GltfArrayType.hxx:38-70`): only the
/// names OCCT maps become an array type; everything else stays unknown and is
/// ignored by the parser (`RWGltf_GltfJsonParser.cxx:1900-1907`).
pub fn parse_attrib_type(name: &str) -> GltfArrayType {
    match name {
        "POSITION" => GltfArrayType::Position,
        "NORMAL" => GltfArrayType::Normal,
        "COLOR" => GltfArrayType::Color,
        "TEXCOORD_0" => GltfArrayType::TCoord0,
        "TEXCOORD_1" => GltfArrayType::TCoord1,
        "JOINT" => GltfArrayType::Joint,
        "WEIGHT" => GltfArrayType::Weight,
        _ => GltfArrayType::Unknown,
    }
}

/// `RWGltf_GltfAccessorLayout` (`RWGltf_GltfAccessorLayout.hxx:22-68`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GltfAccessorLayout {
    /// `RWGltf_GltfAccessorLayout_UNKNOWN`.
    #[default]
    Unknown,
    /// `RWGltf_GltfAccessorLayout_Scalar`.
    Scalar,
    /// `RWGltf_GltfAccessorLayout_Vec2`.
    Vec2,
    /// `RWGltf_GltfAccessorLayout_Vec3`.
    Vec3,
    /// `RWGltf_GltfAccessorLayout_Vec4`.
    Vec4,
    /// `RWGltf_GltfAccessorLayout_Mat2`.
    Mat2,
    /// `RWGltf_GltfAccessorLayout_Mat3`.
    Mat3,
    /// `RWGltf_GltfAccessorLayout_Mat4`.
    Mat4,
}

impl GltfAccessorLayout {
    /// `RWGltf_GltfParseAccessorType` (`RWGltf_GltfAccessorLayout.hxx:33-67`).
    pub fn from_type_str(s: &str) -> Self {
        match s {
            "SCALAR" => Self::Scalar,
            "VEC2" => Self::Vec2,
            "VEC3" => Self::Vec3,
            "VEC4" => Self::Vec4,
            "MAT2" => Self::Mat2,
            "MAT3" => Self::Mat3,
            "MAT4" => Self::Mat4,
            _ => Self::Unknown,
        }
    }

    /// Component count of the layout.
    pub fn nb_components(self) -> usize {
        match self {
            Self::Scalar => 1,
            Self::Vec2 => 2,
            Self::Vec3 => 3,
            Self::Vec4 | Self::Mat2 => 4,
            Self::Mat3 => 9,
            Self::Mat4 => 16,
            Self::Unknown => 0,
        }
    }
}

/// `RWGltf_GltfAccessorCompType` (`RWGltf_GltfAccessorCompType.hxx:19-27`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i32)]
pub enum GltfAccessorCompType {
    /// `RWGltf_GltfAccessorCompType_UNKNOWN`.
    #[default]
    Unknown = 0,
    /// `GL_BYTE`.
    Int8 = 5120,
    /// `GL_UNSIGNED_BYTE`.
    UInt8 = 5121,
    /// `GL_SHORT`.
    Int16 = 5122,
    /// `GL_UNSIGNED_SHORT`.
    UInt16 = 5123,
    /// `GL_UNSIGNED_INT`.
    UInt32 = 5125,
    /// `GL_FLOAT`.
    Float32 = 5126,
}

impl GltfAccessorCompType {
    /// The `componentType` literal; anything else is `UNKNOWN`.
    pub fn from_int(v: i64) -> Self {
        match v {
            5120 => Self::Int8,
            5121 => Self::UInt8,
            5122 => Self::Int16,
            5123 => Self::UInt16,
            5125 => Self::UInt32,
            5126 => Self::Float32,
            _ => Self::Unknown,
        }
    }

    /// Size of one component in bytes.
    pub fn size(self) -> usize {
        match self {
            Self::Int8 | Self::UInt8 => 1,
            Self::Int16 | Self::UInt16 => 2,
            Self::UInt32 | Self::Float32 => 4,
            Self::Unknown => 0,
        }
    }

    /// `RWGltf_GltfJsonParser.cxx:2081-2090`: the accepted component types.
    pub fn is_valid(self) -> bool {
        !matches!(self, Self::Unknown)
    }
}

/// `RWGltf_GltfAccessor` (`RWGltf_GltfAccessor.hxx:29-61`), without the cached
/// bounding box (this port computes a bounding box from the nodes instead of
/// carrying `Graphic3d_BndBox3d` on the accessor).
#[derive(Debug, Clone, Copy, Default)]
pub struct GltfAccessor {
    /// `RWGltf_GltfAccessor::ByteOffset` — added to the bufferView offset.
    pub byte_offset: usize,
    /// `RWGltf_GltfAccessor::Count`.
    pub count: usize,
    /// `RWGltf_GltfAccessor::ByteStride`; `gltfParseBuffer` overwrites it with
    /// the *bufferView* stride (`RWGltf_GltfJsonParser.cxx:2300`, `:2317`), so
    /// this field carries `gltfParseBufferView`'s `aBuffView.ByteStride`.
    pub byte_stride: usize,
    /// `RWGltf_GltfAccessor::Type`.
    pub layout: GltfAccessorLayout,
    /// `RWGltf_GltfAccessor::ComponentType`.
    pub comp_type: GltfAccessorCompType,
}

/// One element of a primitive array: a vertex attribute or the index buffer.
/// `RWGltf_GltfPrimArrayData` (`RWGltf_GltfPrimArrayData.hxx:17-38`).
#[derive(Debug, Clone)]
pub struct GltfPrimArrayData {
    /// `StreamData` — a buffer body decoded at parse time, i.e. one stored as a
    /// base64 `data:` URI (`RWGltf_GltfJsonParser.cxx:2315-2331`).
    pub stream_data: Option<Arc<Vec<u8>>>,
    /// `StreamUri` — the file `RWGltf_TriangulationReader::readFileData` opens
    /// for a deferred load: the GLB itself or an external `.bin`
    /// (`RWGltf_GltfJsonParser.cxx:2294-2305`, `:2357-2359`).
    pub stream_uri: String,
    /// `StreamOffset` — absolute byte offset into `stream_data` / `stream_uri`.
    pub stream_offset: usize,
    /// `StreamLength` — `0` for stream data, the bufferView length otherwise.
    pub stream_length: usize,
    /// The accessor describing the element layout.
    pub accessor: GltfAccessor,
    /// The array this element belongs to.
    pub array_type: GltfArrayType,
}

/// The triangulation a reader produced for a late primitive array — what
/// `Poly_Triangulation` plus `RWMesh_TriangulationSource::myEdges` hold after
/// loading.
#[derive(Debug, Clone)]
pub struct GltfResolvedPrim {
    /// Nodes, indices and UV/normal arrays.
    pub triangulation: Arc<Triangulation>,
    /// `RWMesh_TriangulationSource::myEdges` — the node index of every edge
    /// entry, used by `RWGltf_GltfJsonParser.cxx:1986-1998` to build the line's
    /// `Poly_Polygon3D`.
    pub edges: Vec<usize>,
}

/// `RWGltf_GltfLatePrimitiveArray` (`RWGltf_GltfLatePrimitiveArray.hxx:30-104`)
/// with the deferred-data state of `RWMesh_TriangulationSource`
/// (`RWMesh_TriangulationSource.hxx:68-88`).
#[derive(Debug)]
pub struct GltfLatePrimitiveArray {
    id: String,
    name: String,
    prim_mode: RwGltfPrimitiveMode,
    data: Vec<GltfPrimArrayData>,
    nb_deferred_nodes: usize,
    nb_deferred_triangles: usize,
    degenerated_tri_nb: AtomicUsize,
    /// The triangulation currently attached to the proxy shape. `None` while
    /// the array still holds only descriptors (`myToSkipLateDataLoading`) —
    /// `RWMesh_TriangulationSource` itself reports `NbNodes() == 0` then.
    resolved: Mutex<Option<Arc<GltfResolvedPrim>>>,
}

impl GltfLatePrimitiveArray {
    /// `RWGltf_GltfLatePrimitiveArray(theId, theName)` (`cxx:23-31`).
    pub fn new(
        id: String,
        name: String,
        prim_mode: RwGltfPrimitiveMode,
        data: Vec<GltfPrimArrayData>,
    ) -> Self {
        Self {
            id,
            name,
            prim_mode,
            data,
            nb_deferred_nodes: 0,
            nb_deferred_triangles: 0,
            degenerated_tri_nb: AtomicUsize::new(0),
            resolved: Mutex::new(None),
        }
    }

    /// `Id()` (`hxx:42`).
    pub fn id(&self) -> &str {
        &self.id
    }

    /// `Name()` (`hxx:45`).
    pub fn name(&self) -> &str {
        &self.name
    }

    /// `PrimitiveMode()` (`hxx:51`).
    pub fn prim_mode(&self) -> RwGltfPrimitiveMode {
        self.prim_mode
    }

    /// `Data()` (`hxx:81`).
    pub fn data(&self) -> &[GltfPrimArrayData] {
        &self.data
    }

    /// `NbDeferredNodes()` (`RWMesh_TriangulationSource.hxx:72-76`).
    pub fn nb_deferred_nodes(&self) -> usize {
        self.nb_deferred_nodes
    }

    /// `SetNbDeferredNodes()` (`RWMesh_TriangulationSource.hxx:78-81`).
    pub fn set_nb_deferred_nodes(&mut self, nb: usize) {
        self.nb_deferred_nodes = nb;
    }

    /// `NbDeferredTriangles()` (`RWMesh_TriangulationSource.hxx:83-86`).
    pub fn nb_deferred_triangles(&self) -> usize {
        self.nb_deferred_triangles
    }

    /// `SetNbDeferredTriangles()` (`RWMesh_TriangulationSource.hxx:88-91`).
    pub fn set_nb_deferred_triangles(&mut self, nb: usize) {
        self.nb_deferred_triangles = nb;
    }

    /// `ChangeDegeneratedTriNb()` (`RWMesh_TriangulationSource.hxx:46-48`).
    pub fn degenerated_tri_nb(&self) -> usize {
        self.degenerated_tri_nb.load(Ordering::Relaxed)
    }

    /// `ChangeDegeneratedTriNb() += theNb` (`RWGltf_TriangulationReader.cxx:663`).
    pub fn add_degenerated_tri_nb(&self, nb: usize) {
        self.degenerated_tri_nb.fetch_add(nb, Ordering::Relaxed);
    }

    /// `HasDeferredData()` (`RWGltf_GltfLatePrimitiveArray.hxx:88-91`).
    pub fn has_deferred_data(&self) -> bool {
        !self.data.is_empty()
            && (self.nb_deferred_triangles > 0 || self.nb_deferred_nodes > 0)
    }

    /// The triangulation currently attached to the proxy shape, if any.
    pub fn resolved(&self) -> Option<Arc<GltfResolvedPrim>> {
        self.resolved.lock().ok().and_then(|r| r.clone())
    }

    /// `BRep_Builder::UpdateFace(theFace, theTriangulation)`
    /// (`RWGltf_CafReader.cxx:62-66`): replace the proxy triangulation by the
    /// loaded one.
    pub fn set_resolved(&self, data: Arc<GltfResolvedPrim>) {
        if let Ok(mut r) = self.resolved.lock() {
            *r = Some(data);
        }
    }

    /// The entry count `RWGltf_TriangulationReader::load` walks
    /// (`RWGltf_TriangulationReader.cxx:473-479`).
    pub fn has_stream_data(&self) -> bool {
        self.data.iter().any(|d| d.stream_data.is_some())
    }
}
