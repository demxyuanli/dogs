//! `RWGltf_TriangulationReader` — decodes a deferred primitive array into a
//! triangulation.
//!
//! Source: `RWGltf_TriangulationReader.cxx`. `ReadStream` (`:570-1009`) walks
//! the array elements in the order
//! `RWGltf_GltfLatePrimitiveArray::AddPrimArrayData` stored them — POSITION,
//! INDICES, then the remaining attributes — and dispatches on the array type:
//! positions are converted by the coordinate-system converter, normals through
//! its linear part, `TEXCOORD_0` is flipped in V, indices become triangles or
//! edges depending on the primitive mode. `load` (`:462-513`) reads the
//! file-backed entries, `loadStreamData` (`:192-226`) reads only the entries
//! whose buffer was embedded in the JSON stream, and the two callers
//! `Load`/`LoadStreamData` add `finalizeLoading` (`:515-552`).
//!
//! **UNPORTED** in this reader:
//! * `KHR_draco_mesh_compression` (`:230-457`) — the port has no Draco
//!   decoder, and the parser never marks an accessor `IsCompressed`;
//! * `Poly_Triangulation::SetDoublePrecision` (`:126`) — this port always keeps
//!   `f64` nodes, i.e. the double-precision representation.

use std::io::{Read, Seek, SeekFrom};
use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::poly::triangulation::Triangle;
use occt_core::poly::Triangulation;

use super::coordinate_system::RwMeshCoordinateSystemConverter;
use super::gltf::RwGltfPrimitiveMode;
use super::gltf_prim_array::{
    GltfAccessor, GltfAccessorCompType, GltfAccessorLayout, GltfArrayType, GltfLatePrimitiveArray,
    GltfPrimArrayData, GltfResolvedPrim,
};

/// `THE_NORMAL_PREC2` (`RWGltf_TriangulationReader.cxx:36`): a normal whose
/// square modulus is below this is replaced by `+Z`.
const NORMAL_PREC2: f32 = 0.001;

/// The destination triangulation `RWGltf_TriangulationReader` fills — a
/// `Poly_Triangulation` plus the `RWMesh_TriangulationSource::myEdges` array.
#[derive(Default)]
struct PrimDest {
    nodes: Vec<GpPnt>,
    triangles: Vec<Triangle>,
    uvs: Vec<(f64, f64)>,
    normals: Vec<GpPnt>,
    edges: Vec<usize>,
}

impl PrimDest {
    fn into_resolved(self) -> GltfResolvedPrim {
        let mut tri = Triangulation::new(self.nodes, self.triangles);
        if !self.uvs.is_empty() {
            tri = tri.with_uv(self.uvs);
        }
        if !self.normals.is_empty() {
            tri = tri.with_normals(self.normals);
        }
        GltfResolvedPrim { triangulation: Arc::new(tri), edges: self.edges }
    }
}

/// `RWGltf_TriangulationReader` (`RWGltf_TriangulationReader.hxx`).
pub(super) struct GltfTriangulationReader<'a> {
    converter: &'a RwMeshCoordinateSystemConverter,
    /// `myToSkipDegenerateTris` (`RWMesh_TriangulationReader.hxx:90-92`).
    /// `RWGltf_CafReader::createMeshReaderContext` sets it FALSE
    /// (`RWGltf_CafReader.cxx:362`).
    skip_degenerates: bool,
}

impl<'a> GltfTriangulationReader<'a> {
    /// `createMeshReaderContext` (`RWGltf_CafReader.cxx:355-365`).
    pub(super) fn new(converter: &'a RwMeshCoordinateSystemConverter) -> Self {
        Self { converter, skip_degenerates: false }
    }

    /// `RWGltf_TriangulationReader::LoadStreamData` (`cxx:118-139`) — the
    /// `CafReader_GltfStreamDataLoadingFunctor`, i.e. what
    /// `myToSkipLateDataLoading` still loads: entries whose buffer body is part
    /// of the JSON stream. Returns `None` when the array has no such entry, or
    /// when the result holds no node (`finalizeLoading`).
    pub(super) fn load_stream_data(
        &self,
        late: &GltfLatePrimitiveArray,
    ) -> Result<Option<Arc<GltfResolvedPrim>>, String> {
        if late.prim_mode() == RwGltfPrimitiveMode::Unknown {
            return Ok(None);
        }
        let mut dest = PrimDest::default();
        let mut was_loaded = false;
        for element in late.data() {
            // `cxx:210-213`: stream entries only.
            if element.stream_data.is_none() {
                continue;
            }
            let bytes = stream_bytes(element)?;
            self.read_stream(&mut dest, &bytes, &element.accessor, element.array_type, late)?;
            was_loaded = true;
        }
        if !was_loaded {
            return Ok(None);
        }
        self.finalize(late, dest)
    }

    /// `RWGltf_TriangulationReader::Load` in the `KeepLateData` arm
    /// (`RWGltf_CafReader.cxx:130-135`): `LoadDeferredData` fills the source
    /// object itself, so the proxy face keeps its late array.
    pub(super) fn load_deferred_data(
        &self,
        late: &GltfLatePrimitiveArray,
    ) -> Result<Arc<GltfResolvedPrim>, String> {
        self.load_file_entries(late)
    }

    /// `DetachedLoadDeferredData` (`RWGltf_CafReader.cxx:137`): the loaded
    /// triangulation is a new object the caller puts on the face through
    /// `BRep_Builder::UpdateFace`.
    pub(super) fn detached_load_deferred_data(
        &self,
        late: &GltfLatePrimitiveArray,
    ) -> Result<Arc<GltfResolvedPrim>, String> {
        self.load_file_entries(late)
    }

    /// Parse-time load of a non-triangulation primitive array. `readLateData`
    /// walks `aDoc.FaceList()` only, so line and point shapes are not part of
    /// the deferred pass; their nodes must exist as soon as the shape is built
    /// (`RWGltf_GltfJsonParser.cxx:1960-1967` keeps only triangulations
    /// deferred). Both embedded (base64) and file-backed buffers are read here,
    /// which is what `LoadStreamData` + `LoadDeferredData` cover together.
    pub(super) fn load_all_data(
        &self,
        late: &GltfLatePrimitiveArray,
    ) -> Result<Arc<GltfResolvedPrim>, String> {
        if late.prim_mode() == RwGltfPrimitiveMode::Unknown {
            return Err(format!("gltf: buffer '{}' has an unknown primitive mode", late.id()));
        }
        let mut dest = PrimDest::default();
        for element in late.data() {
            let bytes = match &element.stream_data {
                Some(_) => stream_bytes(element)?,
                None => {
                    if element.stream_uri.is_empty() {
                        return Err(format!("gltf: buffer '{}' does not define uri", late.id()));
                    }
                    file_bytes(element)?
                }
            };
            self.read_stream(&mut dest, &bytes, &element.accessor, element.array_type, late)?;
        }
        self.finalize(late, dest)?
            .ok_or_else(|| format!("gltf: buffer '{}' defines an empty array", late.id()))
    }

    /// `RWGltf_TriangulationReader::load` (`cxx:462-513`).
    fn load_file_entries(
        &self,
        late: &GltfLatePrimitiveArray,
    ) -> Result<Arc<GltfResolvedPrim>, String> {
        if late.prim_mode() == RwGltfPrimitiveMode::Unknown {
            return Err(format!("gltf: buffer '{}' has an unknown primitive mode", late.id()));
        }
        let mut dest = PrimDest::default();
        for element in late.data() {
            if element.stream_data.is_some() {
                // `cxx:481-486`: stream data cannot be part of a deferred load.
                continue;
            }
            if element.stream_uri.is_empty() {
                return Err(format!("gltf: buffer '{}' does not define uri", late.id()));
            }
            let bytes = file_bytes(element)?;
            self.read_stream(&mut dest, &bytes, &element.accessor, element.array_type, late)?;
        }
        self.finalize(late, dest)?
            .ok_or_else(|| format!("gltf: buffer '{}' defines an empty array", late.id()))
    }

    /// `RWMesh_TriangulationReader::Load`'s `finalizeLoading` combined with
    /// `RWGltf_TriangulationReader::finalizeLoading` (`cxx:515-552`): an empty
    /// node array fails the load, and a triangle primitive without an index
    /// accessor reconstructs its triangle list from consecutive nodes.
    fn finalize(
        &self,
        late: &GltfLatePrimitiveArray,
        mut dest: PrimDest,
    ) -> Result<Option<Arc<GltfResolvedPrim>>, String> {
        if dest.nodes.is_empty() {
            return Ok(None);
        }
        if dest.triangles.is_empty() && late.prim_mode() == RwGltfPrimitiveMode::Triangles {
            // `cxx:533-550`.
            let nb_tris = dest.nodes.len() / 3;
            let mut tris = Vec::with_capacity(nb_tris);
            for i in 0..nb_tris {
                let (a, b, c) = (3 * i, 3 * i + 1, 3 * i + 2);
                match self.set_triangle(&dest, a, b, c) {
                    Some(t) => tris.push(t),
                    None => continue,
                }
            }
            dest.triangles = tris;
        }
        Ok(Some(Arc::new(dest.into_resolved())))
    }

    /// `RWMesh_TriangulationReader::setTriangle` (`RWMesh_TriangulationReader.hxx:259-278`).
    fn set_triangle(&self, dest: &PrimDest, a: usize, b: usize, c: usize) -> Option<Triangle> {
        let nb = dest.nodes.len();
        if a >= nb || b >= nb || c >= nb {
            return None;
        }
        if self.skip_degenerates && (a == b || a == c || b == c) {
            return None;
        }
        Some(Triangle::new(a, b, c))
    }

    /// `RWMesh_TriangulationReader::setEdge` (`RWMesh_TriangulationReader.cxx:179-197`).
    fn set_edge(dest: &mut PrimDest, index: usize, edge: usize) -> bool {
        if edge >= dest.nodes.len() {
            return false;
        }
        if index >= dest.edges.len() {
            dest.edges.resize(index + 1, usize::MAX);
        }
        dest.edges[index] = edge;
        true
    }

    /// `RWGltf_TriangulationReader::ReadStream` (`cxx:570-1009`).
    fn read_stream(
        &self,
        dest: &mut PrimDest,
        bytes: &[u8],
        accessor: &GltfAccessor,
        array_type: GltfArrayType,
        late: &GltfLatePrimitiveArray,
    ) -> Result<(), String> {
        let mode = late.prim_mode();
        if mode != RwGltfPrimitiveMode::Triangles
            && mode != RwGltfPrimitiveMode::Lines
            && mode != RwGltfPrimitiveMode::Points
        {
            // `cxx:580-586`: an unsupported primitive array is skipped, not an
            // error.
            return Ok(());
        }
        let stride = |elem: usize| if accessor.byte_stride != 0 { accessor.byte_stride } else { elem };
        match array_type {
            GltfArrayType::Indices => {
                if accessor.layout != GltfAccessorLayout::Scalar {
                    return Ok(());
                }
                match accessor.comp_type {
                    GltfAccessorCompType::UInt16 => {
                        // `cxx:596-677`.
                        let is_triangles = mode == RwGltfPrimitiveMode::Triangles;
                        let counter = if is_triangles {
                            accessor.count / 3
                        } else {
                            accessor.count
                        };
                        let s = stride(2);
                        let mut invalid = 0usize;
                        for i in 0..counter {
                            if is_triangles {
                                let mut idx = [0usize; 3];
                                for (k, slot) in idx.iter_mut().enumerate() {
                                    let off = i * 3 * s + k * s;
                                    *slot = u16_at(bytes, off)? as usize;
                                }
                                match self.set_triangle(dest, idx[0], idx[1], idx[2]) {
                                    Some(t) => dest.triangles.push(t),
                                    None => invalid += 1,
                                }
                            } else if mode == RwGltfPrimitiveMode::Lines {
                                let off = i * s;
                                let idx = u16_at(bytes, off)? as usize;
                                if !Self::set_edge(dest, i, idx) {
                                    invalid += 1;
                                }
                            }
                        }
                        self.report_invalid(late, invalid, counter)?;
                    }
                    GltfAccessorCompType::UInt32 | GltfAccessorCompType::UInt8 => {
                        // `cxx:679-819`: only the triangle arm reads indices.
                        if mode == RwGltfPrimitiveMode::Triangles {
                            let nb_tris = accessor.count / 3;
                            let elem = accessor.comp_type.size();
                            let s = stride(elem);
                            let mut invalid = 0usize;
                            for i in 0..nb_tris {
                                let mut idx = [0usize; 3];
                                for (k, slot) in idx.iter_mut().enumerate() {
                                    let off = i * 3 * s + k * s;
                                    *slot = match accessor.comp_type {
                                        GltfAccessorCompType::UInt8 => u8_at(bytes, off)? as usize,
                                        _ => u32_at(bytes, off)? as usize,
                                    };
                                }
                                match self.set_triangle(dest, idx[0], idx[1], idx[2]) {
                                    Some(t) => dest.triangles.push(t),
                                    None => invalid += 1,
                                }
                            }
                            self.report_invalid(late, invalid, nb_tris)?;
                        }
                    }
                    _ => {}
                }
                Ok(())
            }
            GltfArrayType::Position => {
                // `cxx:825-884`: Float32 + Vec3 only.
                if accessor.comp_type != GltfAccessorCompType::Float32
                    || accessor.layout != GltfAccessorLayout::Vec3
                {
                    return Ok(());
                }
                dest.nodes = vec![GpPnt::new(0.0, 0.0, 0.0); accessor.count];
                let s = stride(12);
                for i in 0..accessor.count {
                    let off = i * s;
                    let mut p = GpPnt::new(
                        f32_at(bytes, off)? as f64,
                        f32_at(bytes, off + 4)? as f64,
                        f32_at(bytes, off + 8)? as f64,
                    );
                    if !self.converter.is_empty() {
                        // `cxx:854-868`: positions go through the converter.
                        self.converter.transform_position(&mut p.coord);
                    }
                    dest.nodes[i] = p;
                }
                Ok(())
            }
            GltfArrayType::Normal => {
                // `cxx:886-957`.
                if accessor.comp_type != GltfAccessorCompType::Float32
                    || accessor.layout != GltfAccessorLayout::Vec3
                {
                    return Ok(());
                }
                if dest.nodes.len() != accessor.count {
                    return Err(format!(
                        "gltf: buffer '{}' has a normal array of a different size",
                        late.id()
                    ));
                }
                dest.normals = vec![GpPnt::new(0.0, 0.0, 0.0); accessor.count];
                let s = stride(12);
                for i in 0..accessor.count {
                    let off = i * s;
                    let mut n = [
                        f32_at(bytes, off)?,
                        f32_at(bytes, off + 4)?,
                        f32_at(bytes, off + 8)?,
                    ];
                    if n[0] * n[0] + n[1] * n[1] + n[2] * n[2] >= NORMAL_PREC2 {
                        if !self.converter.is_empty() {
                            let mut v = occt_core::gp::GpXyz::new(
                                n[0] as f64, n[1] as f64, n[2] as f64,
                            );
                            self.converter.transform_normal(&mut v);
                            n = [v.x as f32, v.y as f32, v.z as f32];
                        }
                    } else {
                        // `cxx:925-929`: a degenerate normal becomes +Z.
                        n = [0.0, 0.0, 1.0];
                    }
                    dest.normals[i] = GpPnt::new(n[0] as f64, n[1] as f64, n[2] as f64);
                }
                Ok(())
            }
            GltfArrayType::TCoord0 => {
                // `cxx:958-997`: Float32 + Vec2 only.
                if accessor.comp_type != GltfAccessorCompType::Float32
                    || accessor.layout != GltfAccessorLayout::Vec2
                {
                    return Ok(());
                }
                if dest.nodes.len() != accessor.count {
                    return Err(format!(
                        "gltf: buffer '{}' has a TEXCOORD_0 array of a different size",
                        late.id()
                    ));
                }
                dest.uvs = vec![(0.0, 0.0); accessor.count];
                let s = stride(8);
                for i in 0..accessor.count {
                    let off = i * s;
                    let u = f32_at(bytes, off)? as f64;
                    let v = f32_at(bytes, off + 4)? as f64;
                    // `cxx:992-993`: glTF images are top-down, OCCT's bottom-up.
                    dest.uvs[i] = (u, 1.0 - v);
                }
                Ok(())
            }
            // `cxx:998-1002`: attributes this reader does not store.
            GltfArrayType::Color
            | GltfArrayType::TCoord1
            | GltfArrayType::Joint
            | GltfArrayType::Weight => Ok(()),
            // `cxx:1003-1005`.
            GltfArrayType::Unknown => Err(format!(
                "gltf: buffer '{}' has an unknown array type",
                late.id()
            )),
        }
    }

    /// `cxx:656-676` / `:728-747` / `:798-818`: a fully degenerate index array
    /// drops the primitive; a partial one keeps the valid triangles.
    fn report_invalid(
        &self,
        late: &GltfLatePrimitiveArray,
        invalid: usize,
        total: usize,
    ) -> Result<(), String> {
        if invalid == 0 {
            return Ok(());
        }
        late.add_degenerated_tri_nb(invalid);
        if invalid == total {
            return Err(format!(
                "gltf: buffer '{}' has been skipped (all elements are degenerate)",
                late.id()
            ));
        }
        Ok(())
    }
}

/// `RWGltf_TriangulationReader::readStreamData` (`cxx:144-157`): the bytes
/// start at `StreamOffset` inside the decoded buffer body.
fn stream_bytes(element: &GltfPrimArrayData) -> Result<Vec<u8>, String> {
    let body = element
        .stream_data
        .as_ref()
        .ok_or_else(|| "gltf: missing stream data".to_string())?;
    if element.stream_offset > body.len() {
        return Err("gltf: stream offset is out of range".to_string());
    }
    Ok(body[element.stream_offset..].to_vec())
}

/// `RWGltf_TriangulationReader::readFileData` (`cxx:160-189`): seeks to
/// `StreamOffset` in `StreamUri` and reads `StreamLength` bytes.
fn file_bytes(element: &GltfPrimArrayData) -> Result<Vec<u8>, String> {
    let mut file = std::fs::File::open(&element.stream_uri).map_err(|e| {
        format!(
            "gltf: buffer refers to invalid file '{}': {e}",
            element.stream_uri
        )
    })?;
    file.seek(SeekFrom::Start(element.stream_offset as u64)).map_err(|e| {
        format!("gltf: buffer '{}': seek failed: {e}", element.stream_uri)
    })?;
    let len = element.stream_length;
    if len == 0 {
        let mut out = Vec::new();
        file.read_to_end(&mut out)
            .map_err(|e| format!("gltf: buffer '{}': read failed: {e}", element.stream_uri))?;
        return Ok(out);
    }
    let mut out = vec![0u8; len];
    file.read_exact(&mut out)
        .map_err(|e| format!("gltf: buffer '{}': read failed: {e}", element.stream_uri))?;
    Ok(out)
}

fn f32_at(b: &[u8], off: usize) -> Result<f32, String> {
    let s = b
        .get(off..off + 4)
        .ok_or_else(|| "gltf: buffer reading error".to_string())?;
    Ok(f32::from_le_bytes([s[0], s[1], s[2], s[3]]))
}

fn u8_at(b: &[u8], off: usize) -> Result<u8, String> {
    b.get(off).copied().ok_or_else(|| "gltf: buffer reading error".to_string())
}

fn u16_at(b: &[u8], off: usize) -> Result<u16, String> {
    let s = b
        .get(off..off + 2)
        .ok_or_else(|| "gltf: buffer reading error".to_string())?;
    Ok(u16::from_le_bytes([s[0], s[1]]))
}

fn u32_at(b: &[u8], off: usize) -> Result<u32, String> {
    let s = b
        .get(off..off + 4)
        .ok_or_else(|| "gltf: buffer reading error".to_string())?;
    Ok(u32::from_le_bytes([s[0], s[1], s[2], s[3]]))
}
