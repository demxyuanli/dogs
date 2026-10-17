//! `Image_Texture` — texture image definition. Source: `Image_Texture.hxx:29-131`,
//! `Image_Texture.cxx`.
//!
//! The port carries the subset `RWMesh_MaterialMap::CopyTexture` and
//! `RWObj_ObjMaterialMap::DefineMaterial` read:
//!
//! * [`ImageTexture::new_file`] — `Image_Texture(theFileName)` (`cxx:24-39`).
//! * [`ImageTexture::new_file_part`] — `Image_Texture(theFileName, theOffset,
//!   theLength)` (`cxx:43-54`).
//! * [`ImageTexture::new_buffer`] — `Image_Texture(theBuffer, theId)`
//!   (`cxx:58-70`).
//! * [`ImageTexture::probe_image_file_format`] — `ProbeImageFileFormat`
//!   (`cxx:226-302`).
//! * [`ImageTexture::write_image`] — `WriteImage(const AsciiString&)`
//!   (`cxx:306-330`) over `WriteImage(ostream&, ...)` (`cxx:334-400`).
//! * [`ImageTexture::mime_type`] — `MimeType` (`cxx:202-222`), the `mimeType`
//!   the GLB image buffer views carry.
//!
//! **UNPORTED**: `ReadImage` / `ReadCompressedImage` (they hand back an
//! `Image_PixMap`, which this crate has no counterpart for) and `DumpJson`.
//! `TextureId()` *is* ported because the
//! `equal_to<handle<Image_Texture>>` specialisation (`hxx:105-118`) compares it,
//! which is how `RWMesh_MaterialMap`'s `myImageMap` / `myImageFailMap` decide
//! whether a texture was already handled.

use std::path::Path;

/// `Image_Texture` (`Image_Texture.hxx:29`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ImageTexture {
    /// `myTextureId` (`hxx:120`).
    texture_id: String,
    /// `myImagePath` (`hxx:121`).
    image_path: String,
    /// `myBuffer` (`hxx:122`) — the encoded image bytes, `None` for a file-backed
    /// texture.
    buffer: Option<Vec<u8>>,
    /// `myOffset` (`hxx:123`), `-1` when the whole file is the image.
    offset: i64,
    /// `myLength` (`hxx:124`), `-1` when the file length is not bounded.
    length: i64,
}

impl ImageTexture {
    /// `Image_Texture(const TCollection_AsciiString& theFileName)`
    /// (`Image_Texture.cxx:24-39`): `myOffset = myLength = -1`, and the id is
    /// `"texture://" + theFileName` so that two textures on the same file are
    /// equal (the comment at `cxx:33`).
    pub fn new_file(file_name: &str) -> Self {
        Self {
            texture_id: if file_name.is_empty() {
                String::new()
            } else {
                format!("texture://{file_name}")
            },
            image_path: file_name.to_string(),
            buffer: None,
            offset: -1,
            length: -1,
        }
    }

    /// `Image_Texture(theFileName, theOffset, theLength)` (`cxx:43-54`): the id
    /// appends `";<offset>,<length>"` (`cxx:50-53`) because the same file can
    /// hold several images.
    pub fn new_file_part(file_name: &str, offset: i64, length: i64) -> Self {
        Self {
            texture_id: if file_name.is_empty() {
                String::new()
            } else {
                format!("texture://{file_name};{offset},{length}")
            },
            image_path: file_name.to_string(),
            buffer: None,
            offset,
            length,
        }
    }

    /// `Image_Texture(theBuffer, theId)` (`cxx:58-70`): an in-memory encoded
    /// image, id `"texturebuf://" + theId`.
    pub fn new_buffer(data: Vec<u8>, id: &str) -> Self {
        Self {
            texture_id: if id.is_empty() {
                String::new()
            } else {
                format!("texturebuf://{id}")
            },
            image_path: String::new(),
            buffer: Some(data),
            offset: -1,
            length: -1,
        }
    }

    /// `Image_Texture::TextureId` (`hxx:56`).
    pub fn texture_id(&self) -> &str {
        &self.texture_id
    }

    /// `Image_Texture::FilePath` (`hxx:59`).
    pub fn file_path(&self) -> &str {
        &self.image_path
    }

    /// `Image_Texture::FileOffset` (`hxx:62`).
    pub fn file_offset(&self) -> i64 {
        self.offset
    }

    /// `Image_Texture::FileLength` (`hxx:65`).
    pub fn file_length(&self) -> i64 {
        self.length
    }

    /// `Image_Texture::ProbeImageFileFormat` (`Image_Texture.cxx:226-302`): the
    /// image format sniffed from the first 20 bytes, `""` when unknown.
    ///
    /// Every failure path in the C++ (unreadable file, bad offset, short read)
    /// returns the empty string after `Message::SendFail` (`cxx:241-265`); the
    /// port has no message stream, so the empty string is the whole result.
    pub fn probe_image_file_format(&self) -> String {
        // `static const size_t THE_PROBE_SIZE = 20;` (`cxx:228`).
        const THE_PROBE_SIZE: usize = 20;
        let mut probe = [0u8; THE_PROBE_SIZE];

        if let Some(buffer) = &self.buffer {
            // `memcpy(aBuffer, myBuffer->Data(), min(Size(), 20))`
            // (`cxx:231-235`). The C++ copies fewer than 20 bytes when the
            // buffer is shorter and leaves the tail uninitialized; the port
            // zero-fills it, which only affects the bytes past the signature.
            let n = buffer.len().min(THE_PROBE_SIZE);
            probe[..n].copy_from_slice(&buffer[..n]);
        } else {
            let mut file = match std::fs::File::open(&self.image_path) {
                Ok(f) => f,
                Err(_) => return String::new(),
            };
            if self.offset >= 0 {
                // `aFileIn->seekg(myOffset, beg)` + `!good()` (`cxx:247-257`).
                use std::io::Seek;
                if file.seek(std::io::SeekFrom::Start(self.offset as u64)).is_err() {
                    return String::new();
                }
            }
            // `aFileIn->read(aBuffer, 20)` (`cxx:259-264`): a short read fails,
            // so a file shorter than 20 bytes has no probed format.
            use std::io::Read;
            let mut got = 0usize;
            while got < THE_PROBE_SIZE {
                match file.read(&mut probe[got..]) {
                    Ok(0) => return String::new(),
                    Ok(n) => got += n,
                    Err(_) => return String::new(),
                }
            }
        }

        // The signature table (`cxx:267-301`), compared in the same order.
        if probe[0..8] == [0x89, b'P', b'N', b'G', b'\r', b'\n', 0x1A, b'\n'] {
            return "png".to_string();
        } else if probe[0..3] == [0xFF, 0xD8, 0xFF] {
            return "jpg".to_string();
        } else if &probe[0..6] == b"GIF87a" || &probe[0..6] == b"GIF89a" {
            return "gif".to_string();
        } else if probe[0..4] == [b'I', b'I', 0x2A, 0x00] || probe[0..4] == [b'M', b'M', 0x00, 0x2A] {
            return "tiff".to_string();
        } else if &probe[0..2] == b"BM" {
            return "bmp".to_string();
        } else if &probe[0..4] == b"RIFF" && &probe[8..12] == b"WEBP" {
            return "webp".to_string();
        } else if &probe[0..4] == b"DDS " {
            return "dds".to_string();
        }
        String::new()
    }

    /// `Image_Texture::WriteImage(const TCollection_AsciiString& theFile)`
    /// (`Image_Texture.cxx:306-330`): create/truncate the file, then delegate to
    /// the stream overload. Returns `false` when the file cannot be created, when
    /// the copy fails, or when the final flush is not good (`cxx:322-327`).
    pub fn write_image(&self, the_file: &str) -> bool {
        let mut out = match std::fs::File::create(the_file) {
            Ok(f) => f,
            Err(_) => return false,
        };
        if !self.write_image_stream(&mut out) {
            return false;
        }
        // `aFileOut->flush(); if (!aFileOut->good())` (`cxx:322-327`).
        use std::io::Write;
        out.flush().is_ok()
    }

    /// `Image_Texture::MimeType` (`Image_Texture.cxx:202-222`): the MIME type
    /// derived from the probed format, `""` when the format is unknown.
    pub fn mime_type(&self) -> String {
        let a_type = self.probe_image_file_format();
        if a_type == "jpg" {
            "image/jpeg".to_string()
        } else if a_type == "png"
            || a_type == "bmp"
            || a_type == "webp"
            || a_type == "gif"
            || a_type == "tiff"
        {
            format!("image/{a_type}")
        } else if a_type == "dds" {
            "image/vnd-ms.dds".to_string()
        } else if !a_type.is_empty() {
            format!("image/x-{a_type}")
        } else {
            String::new()
        }
    }

    /// `Image_Texture::WriteImage(std::ostream&, const AsciiString&)`
    /// (`Image_Texture.cxx:334-400`).
    pub(crate) fn write_image_stream<W: std::io::Write>(&self, out: &mut W) -> bool {
        if let Some(buffer) = &self.buffer {
            // `theStream.write(myBuffer->Data(), myBuffer->Size())`
            // (`cxx:336-345`).
            return out.write_all(buffer).is_ok();
        }

        let mut file = match std::fs::File::open(&self.image_path) {
            Ok(f) => f,
            Err(_) => return false,
        };

        use std::io::{Read, Seek, SeekFrom};
        // `int64_t aLen = myLength;` (`cxx:356`).
        let mut len = self.length;
        if self.offset >= 0 {
            // `aFileIn->seekg(myOffset, beg)` + `!good()` (`cxx:357-365`).
            if file.seek(SeekFrom::Start(self.offset as u64)).is_err() {
                return false;
            }
        } else {
            // `seekg(0, end); aLen = tellg(); seekg(0, beg)` (`cxx:369-373`).
            let end = match file.seek(SeekFrom::End(0)) {
                Ok(p) => p,
                Err(_) => return false,
            };
            len = end as i64;
            if file.seek(SeekFrom::Start(0)).is_err() {
                return false;
            }
        }

        // `int aChunkSize = 4096; NCollection_Array1<uint8_t> aBuffer(0, 4095);`
        // (`cxx:375-376`). The chunk loop is expressed as a plain read-to-end
        // copy: every write below is a straight `theStream.write` of exactly the
        // bytes that were read (`cxx:377-391`), so the byte stream is identical.
        let mut remaining = len.max(0) as u64;
        let mut chunk = vec![0u8; 4096];
        while remaining > 0 {
            let want = remaining.min(chunk.len() as u64) as usize;
            match file.read(&mut chunk[..want]) {
                Ok(0) => return false,
                Ok(n) => {
                    if out.write_all(&chunk[..n]).is_err() {
                        return false;
                    }
                    remaining -= n as u64;
                }
                Err(_) => return false,
            }
        }
        true
    }

    /// `OSD_Path::FolderAndFileFromPath` (`OSD_Path.cxx`): the folder part (empty
    /// when the path carries none) and the trailing file name. Shared by
    /// `RWMesh_MaterialMap`'s constructor and `CopyTexture`.
    pub(crate) fn folder_and_file(path: &str) -> (String, String) {
        let normalized = path.replace('\\', "/");
        match normalized.rfind('/') {
            Some(i) => (normalized[..i].to_string(), normalized[i + 1..].to_string()),
            None => (String::new(), normalized),
        }
    }

    /// `OSD_Path::FileNameAndExtension` (`OSD_Path.cxx`): the file name split at
    /// its last dot, extension without the dot.
    pub(crate) fn file_name_and_extension(file: &str) -> (String, String) {
        match file.rfind('.') {
            Some(i) => (file[..i].to_string(), file[i + 1..].to_string()),
            None => (file.to_string(), String::new()),
        }
    }

    /// `OSD_Directory::Exists` for a bare path (no `OSD_Path` wrapper needed).
    pub(crate) fn is_dir(path: &str) -> bool {
        Path::new(path).is_dir()
    }
}
