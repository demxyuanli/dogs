//! `RWObj_ObjWriterContext` — the low-level OBJ record writer.
//!
//! Source: `RWObj_ObjWriterContext.cxx`. Split out of `obj_caf_writer/mod.rs`
//! unchanged (2026-10-17).

// ---------------------------------------------------------------------------
// RWObj_ObjWriterContext
// ---------------------------------------------------------------------------

/// `RWObj_ObjWriterContext` (`RWObj_ObjWriterContext.hxx:30-105`): the
/// low-level OBJ record writer. Writes into an in-memory buffer rather than a
/// `FILE*`, but keeps every OCCT control branch that this crate can reach.
pub struct ObjWriterContext {
    /// `NbFaces` (`RWObj_ObjWriterContext.hxx:92`).
    pub nb_faces: i32,
    /// `myFile` — the in-memory port always "opens" successfully.
    content: String,
    /// `myActiveMaterial` (`RWObj_ObjWriterContext.hxx:97`).
    active_material: String,
    /// `myElemPosFirst` / `myElemNormFirst` / `myElemUVFirst`
    /// (`RWObj_ObjWriterContext.hxx:99-101`) — `NCollection_Vec4<int>` with
    /// `xyz()` used as the index shift.
    elem_pos_first: [i32; 3],
    elem_norm_first: [i32; 3],
    elem_uv_first: [i32; 3],
    /// `myHasNormals` / `myHasTexCoords` (`RWObj_ObjWriterContext.hxx:103-104`).
    has_normals: bool,
    has_tex_coords: bool,
}

impl ObjWriterContext {
    /// `RWObj_ObjWriterContext::RWObj_ObjWriterContext`
    /// (`RWObj_ObjWriterContext.cxx:60-76`).
    pub fn new(_name: &str) -> Self {
        Self {
            nb_faces: 0,
            content: String::new(),
            active_material: String::new(),
            elem_pos_first: [1, 1, 1],
            elem_norm_first: [1, 1, 1],
            elem_uv_first: [1, 1, 1],
            has_normals: false,
            has_tex_coords: false,
        }
    }

    /// `RWObj_ObjWriterContext::Close` (`cxx:92-98`).
    pub fn close(&mut self) -> bool {
        true
    }

    /// `RWObj_ObjWriterContext::HasNormals` (`hxx:46`).
    pub fn has_normals(&self) -> bool {
        self.has_normals
    }

    /// `RWObj_ObjWriterContext::SetNormals` (`hxx:48`).
    pub fn set_normals(&mut self, has: bool) {
        self.has_normals = has;
    }

    /// `RWObj_ObjWriterContext::HasTexCoords` (`hxx:54`).
    pub fn has_tex_coords(&self) -> bool {
        self.has_tex_coords
    }

    /// `RWObj_ObjWriterContext::SetTexCoords` (`hxx:57`).
    pub fn set_tex_coords(&mut self, has: bool) {
        self.has_tex_coords = has;
    }

    /// `RWObj_ObjWriterContext::ActiveMaterial` (`hxx:66`).
    pub fn active_material(&self) -> &str {
        &self.active_material
    }

    /// `RWObj_ObjWriterContext::WriteHeader` (`cxx:100-141`). `file_info` is
    /// empty on the `RWObj_CafWriter` path (the DRAWEXE call passes none).
    pub fn write_header(
        &mut self,
        nb_nodes: i32,
        nb_elems: i32,
        mat_lib: &str,
        file_info: &[(String, String)],
    ) -> bool {
        self.content.push_str(
            "# Exported by Open CASCADE Technology [dev.opencascade.org]\n",
        );
        self.content.push_str(&format!("#  Vertices: {nb_nodes}\n"));
        self.content.push_str(&format!("#     Faces: {nb_elems}\n"));
        for (key, value) in file_info {
            // `splitLines` (`cxx:22-57`) + the key/value join (`cxx:112-136`).
            let key_lines: Vec<&str> = split_lines(key);
            let val_lines: Vec<&str> = split_lines(value);
            for (i, l) in key_lines.iter().enumerate() {
                self.content.push_str(if i > 0 { "\n# " } else { "# " });
                self.content.push_str(l);
            }
            self.content.push_str(if key_lines.is_empty() { "# " } else { ":" });
            for (i, l) in val_lines.iter().enumerate() {
                self.content.push_str(if i > 0 { "\n# " } else { " " });
                self.content.push_str(l);
            }
            self.content.push('\n');
        }
        if !mat_lib.is_empty() {
            self.content.push_str(&format!("mtllib {mat_lib}\n"));
        }
        true
    }

    /// `RWObj_ObjWriterContext::WriteActiveMaterial` (`cxx:145-150`).
    pub fn write_active_material(&mut self, material: &str) -> bool {
        self.active_material = material.to_string();
        if !material.is_empty() {
            self.content.push_str(&format!("usemtl {material}\n"));
        } else {
            self.content.push_str("usemtl\n");
        }
        true
    }

    /// `RWObj_ObjWriterContext::WriteTriangle` (`cxx:154-206`).
    pub fn write_triangle(&mut self, tri: [i32; 3]) -> bool {
        let pos = [tri[0] + self.elem_pos_first[0], tri[1] + self.elem_pos_first[1], tri[2] + self.elem_pos_first[2]];
        if self.has_normals {
            let norm = [
                tri[0] + self.elem_norm_first[0],
                tri[1] + self.elem_norm_first[1],
                tri[2] + self.elem_norm_first[2],
            ];
            if self.has_tex_coords {
                let uv = [
                    tri[0] + self.elem_uv_first[0],
                    tri[1] + self.elem_uv_first[1],
                    tri[2] + self.elem_uv_first[2],
                ];
                self.content.push_str(&format!(
                    "f {}/{}/{} {}/{}/{} {}/{}/{}\n",
                    pos[0], uv[0], norm[0], pos[1], uv[1], norm[1], pos[2], uv[2], norm[2]
                ));
            } else {
                self.content.push_str(&format!(
                    "f {}//{} {}//{} {}//{}\n",
                    pos[0], norm[0], pos[1], norm[1], pos[2], norm[2]
                ));
            }
            return true;
        }
        if self.has_tex_coords {
            let uv = [
                tri[0] + self.elem_uv_first[0],
                tri[1] + self.elem_uv_first[1],
                tri[2] + self.elem_uv_first[2],
            ];
            self.content.push_str(&format!(
                "f {}/{} {}/{} {}/{}\n",
                pos[0], uv[0], pos[1], uv[1], pos[2], uv[2]
            ));
            return true;
        }
        self.content.push_str(&format!("f {} {} {}\n", pos[0], pos[1], pos[2]));
        true
    }

    /// `RWObj_ObjWriterContext::WriteVertex` (`cxx:274-277`).
    pub fn write_vertex(&mut self, value: [f32; 3]) -> bool {
        self.content
            .push_str(&format!("v {:.*} {:.*} {:.*}\n", 6, value[0], 6, value[1], 6, value[2]));
        true
    }

    /// `RWObj_ObjWriterContext::WriteNormal` (`cxx:280-283`).
    pub fn write_normal(&mut self, value: [f32; 3]) -> bool {
        self.content
            .push_str(&format!("vn {:.*} {:.*} {:.*}\n", 6, value[0], 6, value[1], 6, value[2]));
        true
    }

    /// `RWObj_ObjWriterContext::WriteTexCoord` (`cxx:287-290`).
    pub fn write_tex_coord(&mut self, value: [f32; 2]) -> bool {
        self.content.push_str(&format!("vt {:.*} {:.*}\n", 6, value[0], 6, value[1]));
        true
    }

    /// `RWObj_ObjWriterContext::WriteGroup` (`cxx:296-298`).
    pub fn write_group(&mut self, value: &str) -> bool {
        if !value.is_empty() {
            self.content.push_str(&format!("g {value}\n"));
        } else {
            self.content.push_str("g\n");
        }
        true
    }

    /// `RWObj_ObjWriterContext::FlushFace` (`cxx:303-312`).
    pub fn flush_face(&mut self, nb_nodes: i32) {
        for i in 0..3 {
            self.elem_pos_first[i] += nb_nodes;
        }
        if self.has_normals {
            for i in 0..3 {
                self.elem_norm_first[i] += nb_nodes;
            }
        }
        if self.has_tex_coords {
            for i in 0..3 {
                self.elem_uv_first[i] += nb_nodes;
            }
        }
    }

    /// The accumulated OBJ text.
    pub fn content(&self) -> &str {
        &self.content
    }
}

/// `splitLines` (`RWObj_ObjWriterContext.cxx:22-57`): split on `\r` / `\n`,
/// dropping empty pieces, right-adjusting each line.
fn split_lines(text: &str) -> Vec<&str> {
    if text.is_empty() {
        return Vec::new();
    }
    text.split(['\r', '\n']).filter(|l| !l.is_empty()).collect()
}
