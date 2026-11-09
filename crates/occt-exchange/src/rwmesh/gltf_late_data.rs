//! `RWGltf_CafReader::readLateData` — the second phase of a glTF read.
//!
//! Source: `RWGltf_CafReader.cxx:367-404` and the two loading functors
//! (`cxx:38-164`). After the parser has built the proxy faces, this pass walks
//! the face list and gives every face its triangulation:
//!
//! * `myToSkipLateDataLoading` — only the primitive arrays whose buffer was
//!   embedded in the JSON stream are loaded (`RWGltf_TriangulationReader::LoadStreamData`);
//!   everything else stays a proxy to be loaded later.
//! * otherwise each face gets the stream data first and, when there is none,
//!   its file-backed data — in place when `myToKeepLateData` is set
//!   (`LoadDeferredData`), through a new triangulation otherwise
//!   (`DetachedLoadDeferredData` + `BRep_Builder::UpdateFace`).
//! * `myToParallel` spreads the faces over the default thread pool
//!   (`cxx:380-384`, `:391-398`); every face is independent, so this is a
//!   partition of the face list.

use std::sync::Arc;

use super::coordinate_system::RwMeshCoordinateSystemConverter;
use super::gltf_prim_array::GltfLatePrimitiveArray;
use super::gltf_triangulation_reader::GltfTriangulationReader;

/// The `RWGltf_CafReader` flags that steer this pass
/// (`RWGltf_CafReader.cxx:173-184`).
#[derive(Debug, Clone, Copy)]
pub(super) struct LateDataOptions {
    /// `myToSkipLateDataLoading` (`cxx:174`, default false).
    pub(super) to_skip_late_data_loading: bool,
    /// `myToKeepLateData` (`cxx:178`, default true).
    pub(super) to_keep_late_data: bool,
    /// `myToParallel` (`cxx:173`, default false).
    pub(super) to_parallel: bool,
}

impl Default for LateDataOptions {
    fn default() -> Self {
        Self { to_skip_late_data_loading: false, to_keep_late_data: true, to_parallel: false }
    }
}

/// `readLateData` (`RWGltf_CafReader.cxx:367-404`). `faces` is `aDoc.FaceList()`:
/// one proxy face per triangle primitive array, in creation order.
pub(super) fn read_late_data(
    faces: &[Arc<GltfLatePrimitiveArray>],
    converter: &RwMeshCoordinateSystemConverter,
    opts: LateDataOptions,
) -> Result<(), String> {
    if faces.is_empty() {
        return Ok(());
    }
    let reader = GltfTriangulationReader::new(converter);
    if !opts.to_parallel {
        // `cxx:391-398`: a pool of one thread is the sequential path.
        for face in faces {
            load_face(&reader, face, opts)?;
        }
        return Ok(());
    }
    // `cxx:391-398`: `min(theFaces.Length(), NbDefaultThreadsToLaunch())`.
    let nb_threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
        .min(faces.len());
    if nb_threads <= 1 {
        for face in faces {
            load_face(&reader, face, opts)?;
        }
        return Ok(());
    }
    let chunk = faces.len().div_ceil(nb_threads);
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(nb_threads);
        for part in faces.chunks(chunk) {
            let reader = &reader;
            handles.push(scope.spawn(move || -> Result<(), String> {
                for face in part {
                    load_face(reader, face, opts)?;
                }
                Ok(())
            }));
        }
        for handle in handles {
            match handle.join() {
                Ok(res) => res?,
                Err(_) => return Err("gltf: triangulation loading thread panicked".to_string()),
            }
        }
        Ok(())
    })
}

/// `CafReader_GltfBaseLoadingFunctor::operator()` (`cxx:52-77`) plus the
/// `loadData` override (`cxx:117-137`).
fn load_face(
    reader: &GltfTriangulationReader<'_>,
    face: &GltfLatePrimitiveArray,
    opts: LateDataOptions,
) -> Result<(), String> {
    // `cxx:59-66`: `LoadStreamData` then `UpdateFace` when it produced data.
    if let Some(data) = reader.load_stream_data(face)? {
        face.set_resolved(data);
        return Ok(());
    }
    if opts.to_skip_late_data_loading {
        // `CafReader_GltfStreamDataLoadingFunctor::loadData` (`cxx:160-166`):
        // only stream data is loaded; the proxy keeps its deferred descriptors.
        return Ok(());
    }
    // `cxx:130-138`.
    let data = if opts.to_keep_late_data {
        reader.load_deferred_data(face)?
    } else {
        reader.detached_load_deferred_data(face)?
    };
    face.set_resolved(data);
    Ok(())
}
