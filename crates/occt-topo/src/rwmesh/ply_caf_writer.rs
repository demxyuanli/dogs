//! ASCII PLY write side. Source: `RWPly_CafWriter` / `RWPly_PlyWriterContext`.
//!
//! Faithful port of the OCCT document exporter:
//!
//! * [`PlyWriterContext`] — `RWPly_PlyWriterContext` (`RWPly_PlyWriterContext.cxx`):
//!   the header (vertex / face element declarations, optional comments),
//!   `WriteVertex` / `WriteTriangle` / `WriteQuad` and the header-count checks
//!   done on `Close`.
//! * [`PlyCafWriterOptions`] — `RWPly_CafWriter`'s setters
//!   (`RWPly_CafWriter.hxx:39-135`), with the class defaults.
//! * [`perform`] — `RWPly_CafWriter::Perform` (`RWPly_CafWriter.cxx:70-183`)
//!   over an in-memory node list (the same node list `RWObj_CafWriter` walks,
//!   see [`super::ObjExportNode`]).
//!
//! Notes carried over from the C++:
//!
//! * The export is two-pass over the *same* node list (`cxx:143-172`): the first
//!   pass writes every vertex, the second writes the faces with a running vertex
//!   offset (`cxx:270-294`). The surface (part) id restarts at zero in each pass
//!   and advances once per node when `HasPartId()` is set (`cxx:146-161`).
//! * Vertex attributes are decided by the writer options, not by the data: when
//!   `SetNormals(true)` a face without normals still writes three `float` zeros,
//!   and `SetColors(true)` writes the default white `uchar` triple. Only
//!   `WriteIndices` consults the face (`cxx:229-267`).
//! * Colours are the style's colour in **linear** RGB scaled by 255 and
//!   truncated (`cxx:236-244`); the sRGB conversion in the C++ is commented out.
//!
//! Not ported (no counterpart in this crate, kept out rather than faked):
//!
//! * `myCSTrsf` (`cxx:126-133`) — the converter only becomes active when the
//!   document carries a length unit; `BRepModel` has no such attribute, so the
//!   port writes the coordinates it is given (the C++ default when
//!   `XCAFDoc_DocumentTool::GetLengthUnit` fails is the identity as well).
//! * `Message_LazyProgressScope` (`cxx:139-141`) — progress reporting only.
//! * The `TDocStd_Document` traversal itself: [`perform`] takes the node list
//!   `brep_exchange::obj_caf_nodes` already produced from
//!   `XCAFPrs_DocumentExplorer(OnlyLeafNodes)`, exactly like the OBJ port.
//!
//! Number formatting: the C++ writes through `operator<<` with the default
//! stream precision (6 significant digits, `%g`), i.e. "0.1" rather than
//! "0.100000001". [`cxx_g6`] reproduces that formatting, with
//! `Standard_CLocaleSentry`'s C locale decimal point.

use occt_core::gp::GpPnt;

use crate::rwmesh::ObjExportNode;
use crate::xcaf_doc_vis_material::srgb_to_linear;

// ---------------------------------------------------------------------------
// std::ostream numeric formatting
// ---------------------------------------------------------------------------

/// `std::ostream` default float formatting: `printf("%g")` with 6 significant
/// digits, C locale (`Standard_CLocaleSentry`, `RWPly_CafWriter.cxx:124`).
pub fn cxx_g6(value: f64) -> String {
    cxx_g(value, 6)
}

/// `printf("%.<prec>g", value)` for a finite `value`.
pub fn cxx_g(value: f64, prec: usize) -> String {
    if value.is_nan() {
        return "nan".to_string();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.to_string();
    }
    if value == 0.0 {
        return "0".to_string();
    }
    let prec = prec.max(1);
    let mut exp = value.abs().log10().floor() as i32;
    if exp < -4 || exp >= prec as i32 {
        // Scientific: one digit before the point, `prec - 1` after it.
        let mut mantissa = value / 10f64.powi(exp);
        let digits = prec - 1;
        let mut text = format!("{:.*}", digits, mantissa);
        let rounded: f64 = text.parse().unwrap_or(mantissa);
        if rounded.abs() >= 10.0 {
            // Rounding pushed the mantissa to 10.0: renormalize.
            exp += 1;
            mantissa = rounded / 10.0;
            text = format!("{:.*}", digits, mantissa);
        }
        let mantissa = strip_trailing_zeros(&text);
        format!("{}e{}{:02}", mantissa, if exp < 0 { '-' } else { '+' }, exp.abs())
    } else {
        let decimals = (prec as i32 - 1 - exp).max(0) as usize;
        strip_trailing_zeros(&format!("{:.*}", decimals, value))
    }
}

/// `%g` drops the trailing zeros of the fractional part (and a bare point).
fn strip_trailing_zeros(text: &str) -> String {
    if !text.contains('.') {
        return text.to_string();
    }
    let trimmed = text.trim_end_matches('0');
    trimmed.trim_end_matches('.').to_string()
}

// ---------------------------------------------------------------------------
// RWPly_PlyWriterContext
// ---------------------------------------------------------------------------

/// `splitLines` (`RWPly_PlyWriterContext.cxx:22-58`): split on `\r`, `\n` or
/// `\r\n`, drop empty lines and right-adjust (`TCollection_AsciiString::
/// RightAdjust`, i.e. strip trailing spaces).
fn split_lines(text: &str) -> Vec<String> {
    let mut lines = Vec::new();
    for line in text.split('\n') {
        let line = line.strip_suffix('\r').unwrap_or(line);
        let line = line.trim_end();
        if !line.is_empty() {
            lines.push(line.to_string());
        }
    }
    lines
}

/// `RWPly_PlyWriterContext` (`RWPly_PlyWriterContext.hxx:34-148`), writing into
/// an in-memory buffer instead of an `std::ostream`.
pub struct PlyWriterContext {
    content: String,
    name: String,
    nb_header_verts: i32,
    nb_header_elems: i32,
    nb_verts: i32,
    nb_elems: i32,
    surf_id: i32,
    vert_offset: i32,
    is_double_prec: bool,
    has_normals: bool,
    has_colors: bool,
    has_tex_coords: bool,
    has_surf_id: bool,
    is_open: bool,
}

impl Default for PlyWriterContext {
    /// `RWPly_PlyWriterContext::RWPly_PlyWriterContext()` (`cxx:62-76`).
    fn default() -> Self {
        Self {
            content: String::new(),
            name: String::new(),
            nb_header_verts: 0,
            nb_header_elems: 0,
            nb_verts: 0,
            nb_elems: 0,
            surf_id: 0,
            vert_offset: 0,
            is_double_prec: false,
            has_normals: false,
            has_colors: false,
            has_tex_coords: false,
            has_surf_id: false,
            is_open: false,
        }
    }
}

impl PlyWriterContext {
    /// Empty constructor (`cxx:62-76`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `IsDoublePrecision()` (`hxx:44`).
    pub fn is_double_precision(&self) -> bool {
        self.is_double_prec
    }

    /// `SetDoublePrecision()` (`hxx:47`).
    pub fn set_double_precision(&mut self, value: bool) {
        self.is_double_prec = value;
    }

    /// `HasNormals()` (`hxx:52`).
    pub fn has_normals(&self) -> bool {
        self.has_normals
    }

    /// `SetNormals()` (`hxx:55`).
    pub fn set_normals(&mut self, value: bool) {
        self.has_normals = value;
    }

    /// `HasTexCoords()` (`hxx:58`).
    pub fn has_tex_coords(&self) -> bool {
        self.has_tex_coords
    }

    /// `SetTexCoords()` (`hxx:62`).
    pub fn set_tex_coords(&mut self, value: bool) {
        self.has_tex_coords = value;
    }

    /// `HasColors()` (`hxx:65`).
    pub fn has_colors(&self) -> bool {
        self.has_colors
    }

    /// `SetColors()` (`hxx:68`).
    pub fn set_colors(&mut self, value: bool) {
        self.has_colors = value;
    }

    /// `HasSurfaceId()` (`hxx:72`).
    pub fn has_surface_id(&self) -> bool {
        self.has_surf_id
    }

    /// `SetSurfaceId(bool)` (`hxx:75`).
    pub fn set_surface_id_enabled(&mut self, value: bool) {
        self.has_surf_id = value;
    }

    /// `IsOpened()` (`hxx:78`).
    pub fn is_opened(&self) -> bool {
        self.is_open
    }

    /// `Open()` (`cxx:84-106`): the in-memory port always "opens"
    /// successfully, `theName` is kept only for diagnostics.
    pub fn open(&mut self, name: &str) -> bool {
        self.name = name.to_string();
        self.nb_header_verts = 0;
        self.nb_header_elems = 0;
        self.nb_verts = 0;
        self.nb_elems = 0;
        self.is_open = true;
        true
    }

    /// The accumulated text.
    pub fn content(&self) -> &str {
        &self.content
    }

    /// `Close()` (`cxx:110-138`). The C++ reports a header/written count
    /// mismatch as a message but still returns the stream state, so the result
    /// stays `true`; [`Self::count_mismatch`] exposes the condition.
    pub fn close(&mut self, is_aborted: bool) -> bool {
        if !self.is_open {
            return false;
        }
        self.is_open = false;
        let _ = is_aborted;
        true
    }

    /// `myNbVerts != myNbHeaderVerts` / `myNbElems != myNbHeaderElems`
    /// (`cxx:127-134`): `(header verts, written verts, header elems, written
    /// elems)` when they differ.
    pub fn count_mismatch(&self) -> Option<(i32, i32, i32, i32)> {
        if self.nb_verts != self.nb_header_verts || self.nb_elems != self.nb_header_elems {
            Some((self.nb_header_verts, self.nb_verts, self.nb_header_elems, self.nb_elems))
        } else {
            None
        }
    }

    /// `WriteHeader()` (`cxx:142-222`).
    pub fn write_header(
        &mut self,
        nb_nodes: i32,
        nb_elems: i32,
        file_info: &[(String, String)],
    ) -> bool {
        if !self.is_open {
            return false;
        }
        self.nb_header_verts = nb_nodes;
        self.nb_header_elems = nb_elems;
        self.content.push_str("ply\n");
        self.content.push_str("format ascii 1.0\n");
        self.content
            .push_str("comment Exported by Open CASCADE Technology [dev.opencascade.org]\n");
        for (key, value) in file_info {
            let key_lines = split_lines(key);
            let val_lines = split_lines(value);
            for (i, line) in key_lines.iter().enumerate() {
                if i > 0 {
                    self.content.push('\n');
                }
                self.content.push_str("comment ");
                self.content.push_str(line);
            }
            // `*myStream << (!aKeyLines.IsEmpty() ? ":" : "comment ");`
            self.content.push_str(if !key_lines.is_empty() { ":" } else { "comment " });
            for (i, line) in val_lines.iter().enumerate() {
                if i > 0 {
                    self.content.push('\n');
                }
                self.content.push_str("comment ");
                self.content.push_str(line);
            }
            self.content.push('\n');
        }
        self.content.push_str(&format!("element vertex {nb_nodes}\n"));
        if self.is_double_prec {
            self.content.push_str("property double x\n");
            self.content.push_str("property double y\n");
            self.content.push_str("property double z\n");
        } else {
            self.content.push_str("property float x\n");
            self.content.push_str("property float y\n");
            self.content.push_str("property float z\n");
        }
        if self.has_normals {
            self.content.push_str("property float nx\n");
            self.content.push_str("property float ny\n");
            self.content.push_str("property float nz\n");
        }
        if self.has_tex_coords {
            self.content.push_str("property float s\n");
            self.content.push_str("property float t\n");
        }
        if self.has_colors {
            self.content.push_str("property uchar red\n");
            self.content.push_str("property uchar green\n");
            self.content.push_str("property uchar blue\n");
        }
        if nb_elems > 0 {
            self.content.push_str(&format!("element face {nb_elems}\n"));
            self.content.push_str("property list uchar uint vertex_indices\n");
            if self.has_surf_id {
                self.content.push_str("property uint SurfaceID\n");
            }
        }
        self.content.push_str("end_header\n");
        true
    }

    /// `WriteVertex()` (`cxx:226-265`).
    pub fn write_vertex(
        &mut self,
        point: &GpPnt,
        norm: [f32; 3],
        uv: [f32; 2],
        color: [u8; 3],
    ) -> bool {
        if !self.is_open {
            return false;
        }
        if self.is_double_prec {
            self.content.push_str(&format!(
                "{} {} {}",
                cxx_g6(point.x()),
                cxx_g6(point.y()),
                cxx_g6(point.z())
            ));
        } else {
            // `(float)thePoint.X()` — the coordinates are truncated to `float`
            // before the stream sees them.
            self.content.push_str(&format!(
                "{} {} {}",
                cxx_g6(point.x() as f32 as f64),
                cxx_g6(point.y() as f32 as f64),
                cxx_g6(point.z() as f32 as f64)
            ));
        }
        if self.has_normals {
            self.content.push_str(&format!(
                " {} {} {}",
                cxx_g6(norm[0] as f64),
                cxx_g6(norm[1] as f64),
                cxx_g6(norm[2] as f64)
            ));
        }
        if self.has_tex_coords {
            self.content.push_str(&format!(
                " {} {}",
                cxx_g6(uv[0] as f64),
                cxx_g6(uv[1] as f64)
            ));
        }
        if self.has_colors {
            self.content.push_str(&format!(" {} {} {}", color[0], color[1], color[2]));
        }
        self.content.push('\n');
        self.nb_verts += 1;
        // `throw Standard_OutOfRange` (`cxx:258-264`).
        self.nb_verts <= self.nb_header_verts
    }

    /// `VertexOffset()` (`hxx:109`).
    pub fn vertex_offset(&self) -> i32 {
        self.vert_offset
    }

    /// `SetVertexOffset()` (`hxx:112`).
    pub fn set_vertex_offset(&mut self, offset: i32) {
        self.vert_offset = offset;
    }

    /// `SurfaceId()` (`hxx:115`).
    pub fn surface_id(&self) -> i32 {
        self.surf_id
    }

    /// `SetSurfaceId(int)` (`hxx:118`).
    pub fn set_surface_id(&mut self, id: i32) {
        self.surf_id = id;
    }

    /// `NbWrittenVertices()` (`hxx:105`).
    pub fn nb_written_vertices(&self) -> i32 {
        self.nb_verts
    }

    /// `NbWrittenElements()` (`hxx:127`).
    pub fn nb_written_elements(&self) -> i32 {
        self.nb_elems
    }

    /// `WriteTriangle()` (`cxx:268-289`).
    pub fn write_triangle(&mut self, tri: [i32; 3]) -> bool {
        if !self.is_open {
            return false;
        }
        let a = tri[0] + self.vert_offset;
        let b = tri[1] + self.vert_offset;
        let c = tri[2] + self.vert_offset;
        self.content.push_str(&format!("3 {a} {b} {c}"));
        if self.has_surf_id {
            self.content.push_str(&format!(" {}", self.surf_id));
        }
        self.content.push('\n');
        self.nb_elems += 1;
        self.nb_elems <= self.nb_header_elems
    }

    /// `WriteQuad()` (`cxx:293-313`).
    pub fn write_quad(&mut self, quad: [i32; 4]) -> bool {
        if !self.is_open {
            return false;
        }
        let a = quad[0] + self.vert_offset;
        let b = quad[1] + self.vert_offset;
        let c = quad[2] + self.vert_offset;
        let d = quad[3] + self.vert_offset;
        self.content.push_str(&format!("4 {a} {b} {c} {d}"));
        if self.has_surf_id {
            self.content.push_str(&format!(" {}", self.surf_id));
        }
        self.content.push('\n');
        self.nb_elems += 1;
        self.nb_elems <= self.nb_header_elems
    }

    /// `myName` (`hxx:135`), kept for diagnostics.
    pub fn name(&self) -> &str {
        &self.name
    }
}

// ---------------------------------------------------------------------------
// RWPly_CafWriter
// ---------------------------------------------------------------------------

/// `RWPly_CafWriter`'s option setters (`RWPly_CafWriter.hxx`), with the
/// constructor defaults (`cxx:33-41`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PlyCafWriterOptions {
    /// `SetDoublePrecision()`; FALSE by default.
    pub is_double_prec: bool,
    /// `SetNormals()`; TRUE by default.
    pub has_normals: bool,
    /// `SetColors()`; TRUE by default.
    pub has_colors: bool,
    /// `SetTexCoords()`; FALSE by default.
    pub has_tex_coords: bool,
    /// `SetPartId()`; TRUE by default.
    pub has_part_id: bool,
    /// `SetFaceId()`; FALSE by default.
    pub has_face_id: bool,
}

impl Default for PlyCafWriterOptions {
    fn default() -> Self {
        Self {
            is_double_prec: false,
            has_normals: true,
            has_colors: true,
            has_tex_coords: false,
            has_part_id: true,
            has_face_id: false,
        }
    }
}

/// `RWPly_CafWriter::toSkipFaceMesh` (`cxx:50-53`): `RWMesh_FaceIterator::
/// IsEmptyMesh`, i.e. no triangulation nodes or no triangles.
fn to_skip_face_mesh(face: &crate::rwmesh::ObjFaceNode) -> bool {
    face.vertices.is_empty() || face.triangles.is_empty()
}

/// `RWPly_CafWriter::addFaceInfo` (`cxx:186-192`).
fn add_face_info(face: &crate::rwmesh::ObjFaceNode, nb_nodes: &mut i32, nb_elems: &mut i32) {
    *nb_nodes += face.vertices.len() as i32;
    *nb_elems += face.triangles.len() as i32;
}

/// `Quantity_ColorRGBA`'s linear components to the `uchar` triple
/// (`cxx:236-244`): `(unsigned char)int(component * 255.0f)`.
fn color_to_uchar(style: &crate::xcaf_prs::XcafPrsStyle) -> [u8; 3] {
    if !style.has_color() {
        // `NCollection_Vec4<uint8_t> aColorVec(255)` (`cxx:232`).
        return [255, 255, 255];
    }
    let srgb = style.color();
    let linear = [
        srgb_to_linear(srgb.r as f64),
        srgb_to_linear(srgb.g as f64),
        srgb_to_linear(srgb.b as f64),
    ];
    let mut out = [255u8; 3];
    for (i, v) in linear.iter().enumerate() {
        out[i] = (v * 255.0f32) as u8;
    }
    out
}

/// `RWPly_CafWriter::writeNodes` (`cxx:226-267`).
fn write_nodes(writer: &mut PlyWriterContext, face: &crate::rwmesh::ObjFaceNode) -> bool {
    let color = color_to_uchar(&face.style);
    for (i, vertex) in face.vertices.iter().enumerate() {
        let mut norm = [0.0f32; 3];
        if face.has_normals {
            if let Some(n) = face.normals.get(i) {
                norm = [n.x() as f32, n.y() as f32, n.z() as f32];
            }
        }
        let mut uv = [0.0f32; 2];
        if face.has_tex_coords {
            if let Some(t) = face.tex_coords.get(i) {
                uv = [t.x() as f32, t.y() as f32];
            }
        }
        if !writer.write_vertex(vertex, norm, uv, color) {
            return false;
        }
    }
    true
}

/// `RWPly_CafWriter::writeIndices` (`cxx:271-295`).
fn write_indices(
    writer: &mut PlyWriterContext,
    face: &crate::rwmesh::ObjFaceNode,
    options: &PlyCafWriterOptions,
) -> bool {
    if options.has_face_id {
        let id = writer.surface_id() + 1;
        writer.set_surface_id(id);
    }
    for tri in &face.triangles {
        if !writer.write_triangle([tri.n0 as i32, tri.n1 as i32, tri.n2 as i32]) {
            return false;
        }
    }
    let offset = writer.vertex_offset() + face.vertices.len() as i32;
    writer.set_vertex_offset(offset);
    true
}

/// `RWPly_CafWriter::writeShape` (`cxx:195-224`): one node, one pass.
fn write_shape(
    writer: &mut PlyWriterContext,
    step: i32,
    node: &ObjExportNode,
    options: &PlyCafWriterOptions,
) -> bool {
    for face in &node.faces {
        if to_skip_face_mesh(face) {
            continue;
        }
        if step == 0 && !write_nodes(writer, face) {
            return false;
        }
        if step == 1 && !write_indices(writer, face, options) {
            return false;
        }
    }
    true
}

/// `RWPly_CafWriter::Perform` (`RWPly_CafWriter.cxx:70-183`).
///
/// `path` supplies the context name only; the text is returned instead of being
/// written, like `rwmesh::obj_caf_writer::perform`. The two failure exits of the
/// C++ (`theRootLabels.IsEmpty()` / `aNbNodesAll == 0`) become `Err`.
pub fn perform(
    path: &str,
    nodes: &[ObjExportNode],
    file_info: &[(String, String)],
    options: &PlyCafWriterOptions,
) -> Result<String, String> {
    if nodes.is_empty() {
        return Err("Nothing to export into PLY file".to_string());
    }

    let mut nb_nodes_all = 0i32;
    let mut nb_elems_all = 0i32;
    for node in nodes {
        for face in &node.faces {
            if to_skip_face_mesh(face) {
                continue;
            }
            add_face_info(face, &mut nb_nodes_all, &mut nb_elems_all);
        }
    }
    if nb_nodes_all == 0 {
        return Err("No mesh data to save".to_string());
    }

    let mut writer = PlyWriterContext::new();
    writer.set_double_precision(options.is_double_prec);
    writer.set_normals(options.has_normals);
    writer.set_colors(options.has_colors);
    writer.set_tex_coords(options.has_tex_coords);
    writer.set_surface_id_enabled(options.has_part_id || options.has_face_id);
    if !writer.open(path) || !writer.write_header(nb_nodes_all, nb_elems_all, file_info) {
        return Err(format!("Failed to open PLY file\n{path}"));
    }

    // `for (int aStepIter = 0; aStepIter < 2; ++aStepIter)`
    // (`cxx:143-172`).
    for step in 0..2 {
        writer.set_surface_id(0);
        for node in nodes {
            if options.has_part_id {
                let id = writer.surface_id() + 1;
                writer.set_surface_id(id);
            }
            if !write_shape(&mut writer, step, node, options) {
                return Err(format!("Failed to write PLY file\n{path}"));
            }
        }
    }

    let is_closed = writer.close(false);
    if !is_closed {
        return Err(format!("Failed to write PLY file\n{path}"));
    }
    Ok(writer.content().to_string())
}
