use super::prelude::*;
use super::*;

impl<'a> Cursor<'a> {
    pub(super) fn rest(&self) -> &'a str {
        &self.s[self.pos..]
    }
    pub(super) fn skip_ws(&mut self) {
        while let Some(b) = self.s.as_bytes().get(self.pos) {
            if b.is_ascii_whitespace() {
                self.pos += 1;
            } else {
                break;
            }
        }
    }
    pub(super) fn parse_name(&mut self) -> Result<String, String> {
        let b = self.s.as_bytes();
        let start = self.pos;
        while let Some(&c) = b.get(self.pos) {
            if c.is_ascii_alphanumeric() || c == b'_' || c == b'-' || c == b'.' {
                self.pos += 1;
            } else {
                break;
            }
        }
        if self.pos == start {
            return Err("xmlcaf: expected a name".into());
        }
        Ok(self.s[start..self.pos].to_string())
    }
    pub(super) fn parse_element(&mut self) -> Result<XmlEl, String> {
        self.skip_ws();
        if !self.rest().starts_with('<') {
            return Err("xmlcaf: expected '<'".into());
        }
        self.pos += 1;
        let name = self.parse_name()?;
        let mut attrs = Vec::new();
        let mut self_closing = false;
        loop {
            self.skip_ws();
            let rest = self.rest();
            if rest.starts_with("/>") {
                self.pos += 2;
                self_closing = true;
                break;
            }
            if rest.starts_with('>') {
                self.pos += 1;
                break;
            }
            let aname = self.parse_name()?;
            self.skip_ws();
            if !self.rest().starts_with('=') {
                return Err("xmlcaf: expected '=' after attribute name".into());
            }
            self.pos += 1;
            self.skip_ws();
            if !self.rest().starts_with('"') {
                return Err("xmlcaf: expected '\"' before attribute value".into());
            }
            self.pos += 1;
            let end = self.s[self.pos..]
                .find('"')
                .ok_or("xmlcaf: unterminated attribute value")?;
            let raw = &self.s[self.pos..self.pos + end];
            let value = unescape_xml(raw)?;
            self.pos += end + 1;
            attrs.push((aname, value));
        }
        if self_closing {
            return Ok(XmlEl { name, attrs, children: vec![] });
        }
        let mut children = Vec::new();
        loop {
            self.skip_ws();
            let rest = self.rest();
            if rest.starts_with("</") {
                self.pos += 2;
                let cname = self.parse_name()?;
                self.skip_ws();
                if !self.rest().starts_with('>') {
                    return Err("xmlcaf: expected '>' in closing tag".into());
                }
                self.pos += 1;
                if cname != name {
                    return Err(format!("xmlcaf: mismatched closing tag </{cname}> for <{name}>"));
                }
                return Ok(XmlEl { name, attrs, children });
            }
            if rest.is_empty() {
                return Err(format!("xmlcaf: unexpected end of input inside <{name}>"));
            }
            children.push(self.parse_element()?);
        }
    }
}

/// Parse an optional `<?xml …?>` declaration followed by the single root
/// element. Trailing content after the root is rejected.
pub(super) fn parse_xml_doc(s: &str) -> Result<XmlEl, String> {
    let mut c = Cursor { s, pos: 0 };
    c.skip_ws();
    if c.rest().starts_with("<?xml") {
        let end = c.rest().find("?>").ok_or("xmlcaf: unterminated XML declaration")?;
        c.pos += end + 2;
    }
    c.skip_ws();
    let root = c.parse_element()?;
    c.skip_ws();
    if !c.rest().is_empty() {
        return Err("xmlcaf: trailing content after root element".into());
    }
    Ok(root)
}

/// Parse a sequence of sibling elements (used for the inner content of a
/// `<shape>` element).
pub(super) fn parse_element_list(s: &str) -> Result<Vec<XmlEl>, String> {
    let mut c = Cursor { s, pos: 0 };
    let mut els = Vec::new();
    loop {
        c.skip_ws();
        if c.rest().is_empty() {
            return Ok(els);
        }
        els.push(c.parse_element()?);
    }
}

/// Parse the name and attributes out of a single opening tag (which may end
/// in `>` or `/>`). Child content is not consumed.
pub(super) fn parse_attrs(s: &str) -> Result<(String, Vec<(String, String)>), String> {
    let mut c = Cursor { s, pos: 0 };
    c.skip_ws();
    if !c.rest().starts_with('<') {
        return Err("xmlcaf: expected '<'".into());
    }
    c.pos += 1;
    let name = c.parse_name()?;
    let mut attrs = Vec::new();
    loop {
        c.skip_ws();
        let rest = c.rest();
        if rest.starts_with("/>") || rest.starts_with('>') {
            break;
        }
        let aname = c.parse_name()?;
        c.skip_ws();
        if !c.rest().starts_with('=') {
            return Err("xmlcaf: expected '=' after attribute name".into());
        }
        c.pos += 1;
        c.skip_ws();
        if !c.rest().starts_with('"') {
            return Err("xmlcaf: expected '\"' before attribute value".into());
        }
        c.pos += 1;
        let end = c.s[c.pos..]
            .find('"')
            .ok_or("xmlcaf: unterminated attribute value")?;
        let value = unescape_xml(&c.s[c.pos..c.pos + end])?;
        c.pos += end + 1;
        attrs.push((aname, value));
    }
    Ok((name, attrs))
}

/// Length in bytes of the UTF-8 character whose leading byte is `lead`.
pub(super) fn utf8_len(lead: u8) -> usize {
    if lead < 0x80 {
        1
    } else if lead < 0xE0 {
        2
    } else if lead < 0xF0 {
        3
    } else {
        4
    }
}

pub(super) fn entry_from_element(el: &XmlEl) -> Result<XmlEntry, String> {
    if el.name != "entry" {
        return Err(format!("xmlcaf: expected <entry>, found <{}>", el.name));
    }
    let name = el.attr("name").unwrap_or("").to_string();
    let mut attributes = Vec::new();
    let mut shape = None;
    let mut children = Vec::new();
    for c in &el.children {
        match c.name.as_str() {
            "attribute" => {
                let kind = c.attr("kind").unwrap_or("").to_string();
                let value = c.attr("value").unwrap_or("").to_string();
                attributes.push(XmlAttribute { kind, value });
            }
            "shape" => {
                shape = Some(shape_from_element(c)?);
            }
            "entry" => {
                children.push(entry_from_element(c)?);
            }
            other => return Err(format!("xmlcaf: unexpected <{other}> inside <entry>")),
        }
    }
    Ok(XmlEntry {
        name,
        shape,
        attributes,
        children,
    })
}

// ---------------------------------------------------------------------------
// Geometry classification / helpers (mirrors bincaf.rs)
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
