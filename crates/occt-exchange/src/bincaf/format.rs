use super::prelude::*;


// ---- surface / curve kind tags ------------------------------------------

pub(super) const SURF_PLANE: u8 = 0;
pub(super) const SURF_SPHERE: u8 = 1;
pub(super) const CURVE_LINE: u8 = 0;
pub(super) const CURVE_CIRCLE: u8 = 1;

/// 8-byte file magic: `BINXCAF` plus a NUL pad byte.
pub(super) const MAGIC: &[u8; 8] = b"BINXCAF\0";

/// A named attribute attached to an assembly entry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct XcafAttribute {
    /// e.g. `"name"`, `"color"`, `"layer"`, `"material"`.
    pub kind: String,
    /// Attribute payload (free-form string, e.g. `"1.0,0.0,0.0"`).
    pub value: String,
}

/// One node of the assembly tree: an optional shape plus attributes/children.
#[derive(Debug, Clone, Default)]
pub struct BinXcafEntry {
    pub shape: Option<TopoShape>,
    pub attributes: Vec<XcafAttribute>,
    pub children: Vec<BinXcafEntry>,
}

/// Binary XCAF document — a root entry wrapping the whole assembly tree.
#[derive(Debug, Clone, Default)]
pub struct BinXcaf {
    pub root: BinXcafEntry,
}

/// Serialize a `BinXcaf` document to its binary form.
pub fn serialize_bincaf(doc: &BinXcaf) -> Result<Vec<u8>, String> {
    let mut w = BinWriter::new();
    w.bytes(MAGIC);
    w.u32(1); // version
    write_entry(&mut w, &doc.root);
    Ok(w.buf)
}

/// Deserialize a `BinXcaf` document, validating the magic and version and
/// rejecting corrupt lengths with `Err`.
pub fn deserialize_bincaf(bytes: &[u8]) -> Result<BinXcaf, String> {
    if bytes.len() < 12 || &bytes[0..8] != MAGIC {
        return Err("bincaf: bad magic".into());
    }
    let version = u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]);
    if version != 1 {
        return Err(format!("bincaf: unsupported version {version}"));
    }
    let mut r = BinReader::new(&bytes[12..]);
    let root = read_entry(&mut r)?;
    Ok(BinXcaf { root })
}

/// Write a `BinXcaf` document to `path`.
pub fn write_bincaf_file(doc: &BinXcaf, path: &str) -> Result<(), String> {
    let bytes = serialize_bincaf(doc)?;
    std::fs::write(path, &bytes).map_err(|e| format!("bincaf: write {path}: {e}"))
}

/// Read a `BinXcaf` document from `path`.
pub fn read_bincaf_file(path: &str) -> Result<BinXcaf, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("bincaf: read {path}: {e}"))?;
    deserialize_bincaf(&bytes)
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

pub(super) struct BinWriter {
    pub(super) buf: Vec<u8>,
}

impl BinWriter {
    pub(super) fn new() -> Self {
        Self { buf: Vec::new() }
    }
    pub(super) fn u8(&mut self, v: u8) {
        self.buf.push(v);
    }
    pub(super) fn u32(&mut self, v: u32) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }
    pub(super) fn i32(&mut self, v: i32) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }
    pub(super) fn f64(&mut self, v: f64) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }
    pub(super) fn bytes(&mut self, b: &[u8]) {
        self.buf.extend_from_slice(b);
    }
    pub(super) fn str(&mut self, s: &str) {
        let b = s.as_bytes();
        self.u32(b.len() as u32);
        self.bytes(b);
    }
    pub(super) fn pnt(&mut self, p: &GpPnt) {
        self.f64(p.x());
        self.f64(p.y());
        self.f64(p.z());
    }
    pub(super) fn dir(&mut self, d: &GpDir) {
        self.f64(d.x());
        self.f64(d.y());
        self.f64(d.z());
    }
}

pub(super) fn write_entry(w: &mut BinWriter, e: &BinXcafEntry) {
    match &e.shape {
        Some(s) => {
            w.u8(1);
            write_shape(w, s);
        }
        None => w.u8(0),
    }
    w.i32(e.attributes.len() as i32);
    for a in &e.attributes {
        w.str(&a.kind);
        w.str(&a.value);
    }
    w.u32(e.children.len() as u32);
    for c in &e.children {
        write_entry(w, c);
    }
}

pub(super) fn write_shape(w: &mut BinWriter, s: &TopoShape) {
    w.u8(shape_type_id(s.shape_type()));
    match s.shape_type() {
        ShapeType::Vertex => {
            w.pnt(&GeometryRegistry::global().vertex_point(s));
            w.f64(GeometryRegistry::global().vertex_tolerance(s));
        }
        ShapeType::Edge => {
            write_curve(w, s);
            let verts = children_of_type(s, ShapeType::Vertex);
            w.u32(verts.len() as u32);
            for v in &verts {
                write_shape(w, v);
            }
        }
        ShapeType::Wire => {
            let edges = children_of_type(s, ShapeType::Edge);
            w.u32(edges.len() as u32);
            for e in &edges {
                write_shape(w, e);
            }
        }
        ShapeType::Face => {
            write_surface(w, s);
            let wires = children_of_type(s, ShapeType::Wire);
            w.u32(wires.len() as u32);
            for wi in &wires {
                write_shape(w, wi);
            }
        }
        ShapeType::Shell => {
            let faces = children_of_type(s, ShapeType::Face);
            w.u32(faces.len() as u32);
            for f in &faces {
                write_shape(w, f);
            }
        }
        ShapeType::Solid => {
            let shells = children_of_type(s, ShapeType::Shell);
            w.u32(shells.len() as u32);
            for sh in &shells {
                write_shape(w, sh);
            }
        }
        // Compound, CompSolid, Shape: arbitrary child list.
        _ => {
            let kids: Vec<TopoShape> = s
                .tshape
                .read()
                .unwrap()
                .children
                .clone();
            w.u32(kids.len() as u32);
            for k in &kids {
                write_shape(w, k);
            }
        }
    }
}

pub(super) fn write_curve(w: &mut BinWriter, s: &TopoShape) {
    let (lo, hi) = finite_range(GeometryRegistry::global().edge_parameters(s));
    let Some(curve) = GeometryRegistry::global().edge_curve(s) else {
        w.u8(CURVE_LINE);
        w.pnt(&GpPnt::zero());
        w.dir(&dir_x());
        w.f64(lo);
        w.f64(hi);
        return;
    };
    match classify_curve(curve.as_ref(), lo, hi) {
        CurveKind::Line => {
            w.u8(CURVE_LINE);
            let origin = curve.d0(0.0);
            w.pnt(&origin);
            w.dir(&dir_of(&curve.d1(0.0).1));
        }
        CurveKind::Circle => {
            w.u8(CURVE_CIRCLE);
            let p0 = curve.d0(lo);
            let p1 = curve.d0(lo + PI / 2.0);
            let p2 = curve.d0(lo + PI);
            let center = circle_center3(&p0, &p1, &p2).unwrap_or_else(GpPnt::zero);
            let r = center.distance(&p0);
            let n = GpVec::from_pnts(&p0, &p1)
                .xyz()
                .crossed(GpVec::from_pnts(&p0, &p2).xyz());
            let normal = GpDir::from_xyz(&n).unwrap_or(dir_z());
            let xdir = GpDir::from_vec(&GpVec::from_pnts(&center, &p0)).unwrap_or(dir_x());
            w.pnt(&center);
            w.dir(&normal);
            w.dir(&xdir);
            w.f64(r);
        }
        CurveKind::Other => {
            // ponytail: B-spline curves serialize as their tangent line.
            w.u8(CURVE_LINE);
            let origin = curve.d0(lo);
            w.pnt(&origin);
            w.dir(&dir_of(&curve.d1(lo).1));
        }
    }
    w.f64(lo);
    w.f64(hi);
}

pub(super) fn write_surface(w: &mut BinWriter, s: &TopoShape) {
    let Some(surf) = GeometryRegistry::global().face_surface(s) else {
        w.u8(SURF_PLANE);
        w.pnt(&GpPnt::zero());
        w.dir(&dir_z());
        let (u0, u1, v0, v1) = (0.0, 1.0, 0.0, 1.0);
        w.f64(u0);
        w.f64(u1);
        w.f64(v0);
        w.f64(v1);
        return;
    };
    let (u0, u1, v0, v1) = clamp_uv(surf.u_range().0, surf.u_range().1, surf.v_range().0, surf.v_range().1);
    match classify_surface(surf.as_ref()) {
        SurfaceKind::Plane => {
            let pln = face_plane(&Face(s.clone())).unwrap_or_default();
            w.u8(SURF_PLANE);
            w.pnt(&pln.location());
            w.dir(&pln.axis().direction());
        }
        SurfaceKind::Sphere => {
            let center = sphere_center(surf.as_ref()).unwrap_or_else(GpPnt::zero);
            let r = surf.d0(u0, 0.5 * (v0 + v1)).distance(&center);
            w.u8(SURF_SPHERE);
            w.pnt(&center);
            w.f64(r);
        }
        _ => {
            // ponytail: cylinder/cone/torus surfaces store a plane fallback;
            // add ring-sampled parameters when a curved-face export needs them.
            let pln = face_plane(&Face(s.clone())).unwrap_or_default();
            w.u8(SURF_PLANE);
            w.pnt(&pln.location());
            w.dir(&pln.axis().direction());
        }
    }
    w.f64(u0);
    w.f64(u1);
    w.f64(v0);
    w.f64(v1);
}

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

pub(super) struct BinReader<'a> {
    pub(super) buf: &'a [u8],
    pub(super) pos: usize,
}

impl<'a> BinReader<'a> {
    pub(super) fn new(buf: &'a [u8]) -> Self {
        Self { buf, pos: 0 }
    }
    pub(super) fn u8(&mut self) -> Result<u8, String> {
        let v = *self.buf.get(self.pos).ok_or("bincaf: truncated u8")?;
        self.pos += 1;
        Ok(v)
    }
    pub(super) fn u32(&mut self) -> Result<u32, String> {
        let b = self
            .buf
            .get(self.pos..self.pos + 4)
            .ok_or("bincaf: truncated u32")?;
        self.pos += 4;
        Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }
    pub(super) fn i32(&mut self) -> Result<i32, String> {
        let b = self
            .buf
            .get(self.pos..self.pos + 4)
            .ok_or("bincaf: truncated i32")?;
        self.pos += 4;
        Ok(i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }
    pub(super) fn f64(&mut self) -> Result<f64, String> {
        let b = self
            .buf
            .get(self.pos..self.pos + 8)
            .ok_or("bincaf: truncated f64")?;
        self.pos += 8;
        Ok(f64::from_le_bytes([
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7],
        ]))
    }
    pub(super) fn str(&mut self) -> Result<String, String> {
        let len = self.u32()? as usize;
        let end = self.pos.checked_add(len).ok_or("bincaf: string length overflow")?;
        let s = self
            .buf
            .get(self.pos..end)
            .ok_or("bincaf: string length out of bounds")?;
        let s = std::str::from_utf8(s).map_err(|e| format!("bincaf: invalid utf-8: {e}"))?;
        self.pos = end;
        Ok(s.to_string())
    }
    pub(super) fn pnt(&mut self) -> Result<GpPnt, String> {
        let x = self.f64()?;
        let y = self.f64()?;
        let z = self.f64()?;
        Ok(GpPnt::new(x, y, z))
    }
    pub(super) fn dir(&mut self) -> Result<GpDir, String> {
        let x = self.f64()?;
        let y = self.f64()?;
        let z = self.f64()?;
        GpDir::new(x, y, z).map_err(|e| format!("bincaf: bad direction: {e}"))
    }
}

pub(super) fn read_entry(r: &mut BinReader) -> Result<BinXcafEntry, String> {
    let has = r.u8()?;
    let shape = if has == 1 { Some(read_shape(r)?) } else { None };
    let nattrs = r.i32()?;
    if nattrs < 0 {
        return Err("bincaf: negative attribute count".into());
    }
    let mut attributes = Vec::with_capacity(nattrs as usize);
    for _ in 0..nattrs {
        let kind = r.str()?;
        let value = r.str()?;
        attributes.push(XcafAttribute { kind, value });
    }
    let nchildren = r.u32()? as usize;
    let mut children = Vec::with_capacity(nchildren.min(1 << 16));
    for _ in 0..nchildren {
        children.push(read_entry(r)?);
    }
    Ok(BinXcafEntry {
        shape,
        attributes,
        children,
    })
}

pub(super) fn read_shape(r: &mut BinReader) -> Result<TopoShape, String> {
    let id = r.u8()?;
    let ty = shape_type_from_id(id).ok_or("bincaf: unknown shape type")?;
    let b = TopoBuilder::new();
    match ty {
        ShapeType::Vertex => {
            let p = r.pnt()?;
            let tol = r.f64()?;
            Ok(b.make_vertex(p, tol).0)
        }
        ShapeType::Edge => {
            let (curve, lo, hi) = read_curve(r)?;
            let mut e = b.make_edge(curve, lo, hi);
            let vcount = r.u32()? as usize;
            for _ in 0..vcount {
                let v = read_shape(r)?;
                b.add(&mut e.0, &v);
            }
            Ok(e.0)
        }
        ShapeType::Wire => {
            let mut edges = Vec::new();
            let ecount = r.u32()? as usize;
            for _ in 0..ecount {
                edges.push(Edge(read_shape(r)?));
            }
            Ok(b.make_wire(&edges).0)
        }
        ShapeType::Face => {
            let (surface, _uv) = read_surface(r)?;
            let mut wires = Vec::new();
            let wcount = r.u32()? as usize;
            for _ in 0..wcount {
                wires.push(Wire(read_shape(r)?));
            }
            Ok(b.make_face(surface, &wires).0)
        }
        ShapeType::Shell => {
            let mut faces = Vec::new();
            let fcount = r.u32()? as usize;
            for _ in 0..fcount {
                faces.push(Face(read_shape(r)?));
            }
            Ok(b.make_shell(&faces).0)
        }
        ShapeType::Solid => {
            let mut shells = Vec::new();
            let scount = r.u32()? as usize;
            for _ in 0..scount {
                shells.push(Shell(read_shape(r)?));
            }
            Ok(b.make_solid(&shells).0)
        }
        // Compound, CompSolid, Shape: arbitrary child list.
        _ => {
            let mut comp = TopoShape::new(ty);
            let kcount = r.u32()? as usize;
            for _ in 0..kcount {
                let k = read_shape(r)?;
                b.add(&mut comp, &k);
            }
            Ok(comp)
        }
    }
}

pub(super) fn read_curve(r: &mut BinReader) -> Result<(Arc<dyn Curve>, f64, f64), String> {
    let kind = r.u8()?;
    let curve: Arc<dyn Curve> = match kind {
        CURVE_CIRCLE => {
            let center = r.pnt()?;
            let normal = r.dir()?;
            let xdir = r.dir()?;
            let radius = r.f64()?;
            let ax2 = GpAx2::new(center, normal, xdir)
                .map_err(|e| format!("bincaf: bad circle axis: {e}"))?;
            Arc::new(GeomCircle::new(GpCirc::new(ax2, radius)))
        }
        _ => {
            // CURVE_LINE and any unknown tag decode as a line.
            let p = r.pnt()?;
            let d = r.dir()?;
            Arc::new(GeomLine::from_pnt_dir(p, d))
        }
    };
    let lo = r.f64()?;
    let hi = r.f64()?;
    Ok((curve, lo, hi))
}

pub(super) fn read_surface(r: &mut BinReader) -> Result<(Arc<dyn Surface>, (f64, f64, f64, f64)), String> {
    let kind = r.u8()?;
    let surf: Arc<dyn Surface> = match kind {
        SURF_SPHERE => {
            let center = r.pnt()?;
            let radius = r.f64()?;
            let ax3 = GpAx3::from_ax1(&GpAx1::new(center, dir_z()));
            let sphere = GpSphere::new(ax3, radius).map_err(|e| format!("bincaf: {e}"))?;
            Arc::new(GeomSphere::new(sphere))
        }
        _ => {
            // SURF_PLANE and any unknown tag decode as a plane.
            let origin = r.pnt()?;
            let normal = r.dir()?;
            let ax3 = GpAx3::from_ax1(&GpAx1::new(origin, normal));
            Arc::new(GeomPlane::new(GpPln::new(ax3)))
        }
    };
    let u0 = r.f64()?;
    let u1 = r.f64()?;
    let v0 = r.f64()?;
    let v1 = r.f64()?;
    Ok((surf, (u0, u1, v0, v1)))
}

// ---------------------------------------------------------------------------
// Geometry classification / helpers
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum CurveKind {
    Line,
    Circle,
    Other,
}

pub(super) fn classify_curve(c: &dyn Curve, a: f64, b: f64) -> CurveKind {
    let (lo, hi) = if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (0.0, 1.0)
    };
    let span = hi - lo;
    if span < 1e-12 {
        return CurveKind::Other;
    }
    let mut d2s = Vec::with_capacity(6);
    for i in 0..6 {
        let u = lo + span * i as f64 / 5.0;
        d2s.push(c.d2(u).2.magnitude());
    }
    let max_d2 = d2s.iter().cloned().fold(0.0_f64, f64::max);
    if max_d2 < 1e-9 {
        return CurveKind::Line;
    }
    let min_d2 = d2s.iter().cloned().fold(f64::INFINITY, f64::min);
    if c.is_periodic() || (max_d2 - min_d2) / max_d2 < 0.02 {
        CurveKind::Circle
    } else {
        CurveKind::Other
    }
}

pub(super) fn shape_type_id(t: ShapeType) -> u8 {
    match t {
        ShapeType::Compound => 0,
        ShapeType::CompSolid => 1,
        ShapeType::Solid => 2,
        ShapeType::Shell => 3,
        ShapeType::Face => 4,
        ShapeType::Wire => 5,
        ShapeType::Edge => 6,
        ShapeType::Vertex => 7,
        ShapeType::Shape => 8,
    }
}

pub(super) fn shape_type_from_id(id: u8) -> Option<ShapeType> {
    match id {
        0 => Some(ShapeType::Compound),
        1 => Some(ShapeType::CompSolid),
        2 => Some(ShapeType::Solid),
        3 => Some(ShapeType::Shell),
        4 => Some(ShapeType::Face),
        5 => Some(ShapeType::Wire),
        6 => Some(ShapeType::Edge),
        7 => Some(ShapeType::Vertex),
        8 => Some(ShapeType::Shape),
        _ => None,
    }
}

pub(super) fn children_of_type(s: &TopoShape, t: ShapeType) -> Vec<TopoShape> {
    s.tshape
        .read()
        .unwrap()
        .children
        .iter()
        .filter(|h| h.shape_type() == t)
        .cloned()
        .collect()
}

pub(super) fn finite_range((a, b): (f64, f64)) -> (f64, f64) {
    if a.is_finite() && b.is_finite() && b > a {
        (a, b)
    } else {
        (0.0, 1.0)
    }
}

pub(super) fn clamp_uv(u0: f64, u1: f64, v0: f64, v1: f64) -> (f64, f64, f64, f64) {
    let c = |a: f64, b: f64| if a.is_finite() && b.is_finite() && b > a { (a, b) } else { (-1.0, 1.0) };
    let (u0, u1) = c(u0, u1);
    let (v0, v1) = c(v0, v1);
    (u0, u1, v0, v1)
}

pub(super) fn dir_x() -> GpDir {
    GpDir::new(1.0, 0.0, 0.0).expect("unit x")
}
pub(super) fn dir_z() -> GpDir {
    GpDir::new(0.0, 0.0, 1.0).expect("unit z")
}

pub(super) fn dir_of(v: &GpVec) -> GpDir {
    GpDir::from_vec(v).unwrap_or(dir_x())
}

/// 3×3 determinant.
pub(super) fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Circumcenter of three non-collinear 3D points (perpendicular-bisector
/// system solved by Cramer's rule).
pub(super) fn circle_center3(a: &GpPnt, b: &GpPnt, c: &GpPnt) -> Option<GpPnt> {
    let d1 = GpVec::from_pnts(a, b);
    let d2 = GpVec::from_pnts(a, c);
    let n = d1.xyz().crossed(d2.xyz());
    if n.square_modulus() < 1e-30 {
        return None;
    }
    let n2 = |p: &GpPnt| p.x() * p.x() + p.y() * p.y() + p.z() * p.z();
    let mat = [
        [d1.x(), d1.y(), d1.z()],
        [d2.x(), d2.y(), d2.z()],
        [n.x, n.y, n.z],
    ];
    let rhs = [
        0.5 * (n2(b) - n2(a)),
        0.5 * (n2(c) - n2(a)),
        a.x() * n.x + a.y() * n.y + a.z() * n.z,
    ];
    let d = det3(&mat);
    if d.abs() < 1e-30 {
        return None;
    }
    let mut o = [0.0; 3];
    for k in 0..3 {
        let mut m = mat;
        for i in 0..3 {
            m[i][k] = rhs[i];
        }
        o[k] = det3(&m) / d;
    }
    Some(GpPnt::new(o[0], o[1], o[2]))
}
