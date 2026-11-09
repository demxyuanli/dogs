//! glTF container split: JSON chunk plus BIN chunk (`RWGltf_CafReader::performMesh`,
//! `RWGltf_CafReader.cxx:194-286`).
//!
//! The whole file is read into memory here, so the BIN body is kept as its own
//! byte slice and its file offset is recorded for the deferred readers.

/// A glTF payload: the JSON document plus the binary body of a `.glb`.
///
/// `RWGltf_CafReader::performMesh` (`RWGltf_CafReader.cxx:194-286`) reads the
/// container header first and only then hands the JSON chunk to the parser; the
/// BIN chunk never passes through the JSON layer, it is addressed by stream
/// offset (`RWGltf_GltfJsonParser.cxx:2286-2306`). Here the whole file is in
/// memory, so the BIN chunk is kept as a byte slice of its own and stands in for
/// the buffer that carries no `uri`.
pub(super) struct GltfContainer {
    pub(super) json: String,
    pub(super) bin: Option<Vec<u8>>,
    /// Byte offset of the BIN body inside the file (`myBinBodyOffset` in
    /// `RWGltf_GltfJsonParser`), used for file-backed texture slices.
    pub(super) bin_offset: usize,
}

/// Split a glTF `.gltf` text or `.glb` binary container into JSON and BIN.
///
/// `RWGltf_CafReader::performMesh` (`cxx:194-286`):
/// - 12-byte header: `glTF`, uint32 version, uint32 total length;
/// - version 1 (`cxx:210-238`): `contentLength` and `contentFormat` follow, JSON
///   starts at byte 20, the binary body runs to the declared length;
/// - otherwise (`cxx:239-281`): walk at most two chunk headers, taking the JSON
///   chunk (`0x4E4F534A`) and the BIN chunk (`0x004E4942`);
/// - a non-`glTF` prefix (`cxx:283-286`) rewinds and parses the file as JSON.
pub(super) fn gltf_container(bytes: &[u8]) -> Result<GltfContainer, String> {
    let u32_at = |off: usize| -> Result<u32, String> {
        let s = bytes
            .get(off..off + 4)
            .ok_or_else(|| "gltf: truncated container header".to_string())?;
        Ok(u32::from_le_bytes([s[0], s[1], s[2], s[3]]))
    };

    if bytes.len() >= 4 && &bytes[0..4] == b"glTF" {
        let version = u32_at(4)?;
        let length = u32_at(8)? as usize;
        let (mut json_off, mut json_len, mut bin_off, mut bin_len);
        if version == 1 {
            // `cxx:214-236`: `contentLength` < 20 means a broken container; the
            // content format must be 0 (JSON).
            if length < 20 {
                return Err("gltf: broken glTF 1.0 container".into());
            }
            let scene_len = u32_at(12)? as usize;
            let scene_format = u32_at(16)?;
            if scene_format != 0 {
                return Err("gltf: unsupported glTF 1.0 scene format".into());
            }
            json_off = 20;
            json_len = scene_len;
            bin_off = json_off + json_len;
            bin_len = length.saturating_sub(bin_off);
        } else {
            json_off = 0;
            json_len = 0;
            bin_off = 0;
            bin_len = 0;
            // `cxx:246-280`: at most two chunk headers, each `length` +
            // `type`, followed by `length` bytes of payload.
            let mut chunk_off = 12usize;
            for _ in 0..2 {
                if chunk_off + 8 > length || chunk_off + 8 > bytes.len() {
                    break;
                }
                let chunk_len = u32_at(chunk_off)? as usize;
                let chunk_type = u32_at(chunk_off + 4)?;
                let data_off = chunk_off + 8;
                if chunk_type == 0x4E4F_534A {
                    json_off = data_off;
                    json_len = chunk_len;
                } else if chunk_type == 0x004E_4942 {
                    bin_off = data_off;
                    bin_len = chunk_len;
                }
                if chunk_len == 0 {
                    break;
                }
                chunk_off = data_off + chunk_len;
            }
            if json_len == 0 {
                return Err("gltf: .glb defines no JSON chunk".into());
            }
        }

        let json_bytes = bytes
            .get(json_off..json_off + json_len)
            .ok_or_else(|| "gltf: JSON chunk extends past end of file".to_string())?;
        // The spec pads the JSON chunk with spaces; leading bytes could be a BOM.
        let json = String::from_utf8(json_bytes.to_vec())
            .map_err(|_| "gltf: JSON chunk is not UTF-8".to_string())?;
        let bin = if bin_len > 0 {
            Some(
                bytes
                    .get(bin_off..bin_off + bin_len)
                    .ok_or_else(|| "gltf: BIN chunk extends past end of file".to_string())?
                    .to_vec(),
            )
        } else {
            None
        };
        return Ok(GltfContainer { json, bin, bin_offset: if bin_len > 0 { bin_off } else { 0 } });
    }

    let json = String::from_utf8(bytes.to_vec())
        .map_err(|e| format!("rwmesh: gltf is not UTF-8: {e}"))?;
    Ok(GltfContainer { json, bin: None, bin_offset: 0 })
}
