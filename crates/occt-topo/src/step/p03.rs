use super::prelude::*;
use super::*;

/// Shared assembly serialization; `splines` selects the geometry classifier.
pub(super) fn write_assembly_inner(a: &StepAssembly, splines: bool) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = splines;
    let mut defs: HashMap<String, usize> = HashMap::new();

    // The assembly's own product definition is the root that children hang off.
    let top = esc_str(&a.name);
    let top_product = ctx.w.emit(format!("PRODUCT('{top}','{top}','',({}))", ctx.prod_ctx));
    let top_form = ctx.w.emit(format!("PRODUCT_DEFINITION_FORMATION('','',#{top_product})"));
    let top_def = ctx.w.emit(format!("PRODUCT_DEFINITION('','','',#{top_form},#{})", ctx.def_ctx));
    defs.insert(a.name.clone(), top_def);

    for (name, shape) in &a.products {
        let n = esc_str(name);
        let product = ctx.w.emit(format!("PRODUCT('{n}','{n}','',({}))", ctx.prod_ctx));
        let formation = ctx.w.emit(format!("PRODUCT_DEFINITION_FORMATION('','',#{product})"));
        let def = ctx.w.emit(format!("PRODUCT_DEFINITION('','','',#{formation},#{})", ctx.def_ctx));
        let items = ctx.emit_top(shape);
        if !items.is_empty() {
            let rep = ctx.w.emit(format!(
                "ADVANCED_BREP_SHAPE_REPRESENTATION('{n}',({}),#{})",
                join_refs(&items),
                ctx.geom_ctx
            ));
            let pds = ctx.w.emit(format!("PRODUCT_DEFINITION_SHAPE('','',#{product})"));
            ctx.w.emit(format!("PRODUCT_DEFINITION_SHAPE_REPRESENTATION('',#{pds},#{rep})"));
        }
        defs.insert(name.clone(), def);
    }

    for (parent, kids) in &a.children {
        let pdef = defs
            .get(parent)
            .copied()
            .ok_or_else(|| format!("write_step_assembly: unknown parent '{parent}'"))?;
        for kid in kids {
            let kdef = defs
                .get(kid)
                .copied()
                .ok_or_else(|| format!("write_step_assembly: unknown child '{kid}'"))?;
            let kn = esc_str(kid);
            ctx.w.emit(format!(
                "NEXT_ASSEMBLY_USAGE_OCCURRENCE('{kn}','{kn}','',#{pdef},#{kdef},'')"
            ));
        }
    }
    Ok(ctx.finish())
}

/// Serialize a model to STEP, emitting B-spline geometry for non-analytic
/// curves and surfaces (the spline-capable counterpart of [`write_step`]).
pub fn write_step_model_with_splines(model: &BRepModel) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = true;
    for ms in &model.shapes {
        ctx.write_shape_named(&ms.name, &ms.shape);
    }
    Ok(ctx.finish())
}

/// Serialize a model to STEP, attaching a `COLOUR_RGB` style to each shape
/// whose `ModelShape.color` is set.
///
/// Shapes without a colour are written plain. This is the model-level
/// counterpart of [`write_step_with_color`], reading the colours already
/// stored on the [`BRepModel`] (e.g. via `add_with_color`).
pub fn write_step_model_with_colors(model: &BRepModel) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    for ms in &model.shapes {
        let rep = ctx.write_shape_named(&ms.name, &ms.shape);
        if let (Some(rep), Some(c)) = (rep, ms.color) {
            let colour = ctx.w.emit(format!(
                "COLOUR_RGB('',{},{},{})",
                step_real(c.r as f64),
                step_real(c.g as f64),
                step_real(c.b as f64)
            ));
            let fill = ctx.w.emit(format!("SURFACE_STYLE_FILL_AREA('',#{colour})"));
            let usage = ctx.w.emit(format!("SURFACE_STYLE_USAGE('',#{fill})"));
            ctx.w.emit(format!("STYLED_ITEM('',(#{usage}),#{rep})"));
        }
    }
    Ok(ctx.finish())
}

/// Serialize a model with per-shape [`StepWriteOptions`].
///
/// Each shape in the model is written with the shared options; spline output
/// applies to every shape, and the colour/units/name defaults apply as in
/// [`write_step_with_options`]. Shape names come from the model, so `opts.name`
/// is only used as a fallback for unnamed shapes.
pub fn write_step_model_with_options(model: &BRepModel, opts: &StepWriteOptions) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = opts.splines;
    for ms in &model.shapes {
        let name = if ms.name.is_empty() { opts.name.clone() } else { ms.name.clone() };
        let rep = ctx.write_shape_named(&name, &ms.shape);
        if let (Some(rep), Some((r, g, b))) = (rep, opts.color) {
            let colour = ctx.w.emit(format!(
                "COLOUR_RGB('',{},{},{})",
                step_real(r),
                step_real(g),
                step_real(b)
            ));
            let fill = ctx.w.emit(format!("SURFACE_STYLE_FILL_AREA('',#{colour})"));
            let usage = ctx.w.emit(format!("SURFACE_STYLE_USAGE('',#{fill})"));
            ctx.w.emit(format!("STYLED_ITEM('',(#{usage}),#{rep})"));
        }
        emit_si_units(&mut ctx.w, &opts.units);
    }
    Ok(ctx.finish())
}

/// Read a STEP physical file from disk, also returning collected warnings for
/// skipped or unsupported entities.
pub fn read_step_file_with_warnings(path: &str) -> Result<(BRepModel, Vec<String>), String> {
    let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    read_step_with_warnings(&content)
}

/// Serialize a list of named shapes into one STEP file.
///
/// Each `(name, shape)` pair becomes its own representation + product pair,
/// exactly as [`write_step`] does for a [`BRepModel`], without requiring the
/// caller to build a model object first. The reader returns each shape under
/// its `name`.
pub fn write_step_shapes(shapes: &[(String, TopoShape)]) -> String {
    let mut ctx = WriteCtx::new();
    for (name, shape) in shapes {
        ctx.write_shape_named(name, shape);
    }
    ctx.finish()
}

/// Spline-capable counterpart of [`write_step_shapes`]: non-analytic curves
/// and surfaces are emitted as B-splines.
pub fn write_step_shapes_with_splines(shapes: &[(String, TopoShape)]) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = true;
    for (name, shape) in shapes {
        ctx.write_shape_named(name, shape);
    }
    Ok(ctx.finish())
}

/// Serialize a compound as its named children in one STEP file.
///
/// A [`TopoShape`] of type `Compound` is flattened: each child shape is written
/// as a named representation (`Child1`, `Child2`, ...), which the reader
/// reconstructs as separate model shapes. Non-compound shapes are written as a
/// single `"Shape"` representation.
pub fn write_step_compound(compound: &TopoShape) -> Result<String, String> {
    if !compound.is_compound() {
        return write_step_with_splines(compound);
    }
    let kids: Vec<TopoShape> = compound
        .tshape
        .read()
        .unwrap()
        .children
        .clone();
    if kids.is_empty() {
        return Err("write_step_compound: empty compound".into());
    }
    let named: Vec<(String, TopoShape)> = kids
        .into_iter()
        .enumerate()
        .map(|(i, k)| (format!("Child{}", i + 1), k))
        .collect();
    Ok(write_step_shapes_with_splines(&named)?)
}

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

/// A parsed `#N=TYPE(...)` data record.
#[derive(Debug, Clone)]
pub(super) struct Record {
    pub(super) type_name: String,
    pub(super) args: Vec<String>,
}

/// Split a comma-separated argument list at the top nesting level, respecting
/// strings (with `''` escapes), parens and brackets.
pub(super) fn split_top(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut depth = 0usize;
    let mut cur = String::new();
    let mut chars = s.chars().peekable();
    let mut in_str = false;
    while let Some(c) = chars.next() {
        if in_str {
            cur.push(c);
            if c == '\'' {
                if chars.peek() == Some(&'\'') {
                    cur.push(chars.next().unwrap());
                } else {
                    in_str = false;
                }
            }
        } else {
            match c {
                '\'' => {
                    in_str = true;
                    cur.push(c);
                }
                '(' | '[' => {
                    depth += 1;
                    cur.push(c);
                }
                ')' | ']' => {
                    depth -= 1;
                    cur.push(c);
                }
                ',' if depth == 0 => {
                    out.push(cur.trim().to_string());
                    cur.clear();
                }
                _ => cur.push(c),
            }
        }
    }
    if !cur.trim().is_empty() {
        out.push(cur.trim().to_string());
    }
    out
}

/// Parse the entity body `TYPE(a,b,...)` into its type name and argument list.
///
/// A STEP *complex* entity instance groups several subtypes as space-separated
/// members inside one outer paren pair — `( A() B(...) C(...) )` — optionally
/// followed by the compound's own (empty) attribute list. The attributes of the
/// member subtypes form the attribute list of the most derived subtype; this
/// merge reconstructs that record (currently for the B-spline curve/surface
/// families, whose members split the degree/poles/form, the knots, and the
/// weights). Non-B-spline complexes (unit/context metadata) are skipped by the
/// resolver, so they keep an empty type name.
pub(super) fn parse_entity_body(body: &str) -> (String, Vec<String>) {
    let body = body.trim();
    if body.starts_with('(') {
        return merge_complex_body(body);
    }
    let open = body.find('(').unwrap_or(body.len());
    let type_name = body[..open].trim().to_string();
    let close = body.rfind(')').unwrap_or(body.len());
    let inner = if close > open + 1 {
        &body[open + 1..close]
    } else {
        ""
    };
    (type_name, split_top(inner))
}

/// Split `( MEMBER1(args) MEMBER2(args) ... )` into its member type/args pairs.
///
/// Members are separated by whitespace at the top level (not commas); each is a
/// `TYPE(...)` token. The outer parens and any trailing empty compound list are
/// stripped first.
pub(super) fn split_complex_members(body: &str) -> Vec<(String, Vec<String>)> {
    let b = body.trim();
    // Strip the outer pair of parens.
    let inner = if b.starts_with('(') && b.ends_with(')') {
        &b[1..b.len() - 1]
    } else {
        b
    };
    let bytes = inner.as_bytes();
    let mut members: Vec<(String, Vec<String>)> = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        while i < bytes.len() && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= bytes.len() {
            break;
        }
        // Member type name: identifier until `(`.
        let t0 = i;
        while i < bytes.len() && (bytes[i].is_ascii_alphanumeric() || bytes[i] == b'_') {
            i += 1;
        }
        if i >= bytes.len() || bytes[i] != b'(' {
            // Skip stray tokens (not a TYPE(...) member).
            while i < bytes.len() && !bytes[i].is_ascii_whitespace() {
                i += 1;
            }
            continue;
        }
        let type_name = inner[t0..i].to_string();
        // Scan to the matching `)`, respecting strings and nesting.
        let mut depth = 0usize;
        let mut in_str = false;
        let a0 = i;
        while i < bytes.len() {
            let c = bytes[i] as char;
            if in_str {
                if c == '\'' {
                    if i + 1 < bytes.len() && bytes[i + 1] == b'\'' {
                        i += 2;
                        continue;
                    }
                    in_str = false;
                }
            } else {
                match c {
                    '\'' => in_str = true,
                    '(' => depth += 1,
                    ')' => {
                        depth -= 1;
                        if depth == 0 {
                            i += 1;
                            break;
                        }
                    }
                    _ => {}
                }
            }
            i += 1;
        }
        let member_body = &inner[a0..i];
        let (_, args) = {
            let open = member_body.find('(').unwrap_or(0);
            let close = member_body.rfind(')').unwrap_or(member_body.len());
            if close > open + 1 {
                (String::new(), split_top(&member_body[open + 1..close]))
            } else {
                (String::new(), Vec::new())
            }
        };
        members.push((type_name, args));
    }
    members
}

/// Reconstruct a single record from a STEP complex entity body.
///
/// See [`parse_entity_body`] for the format. B-spline curves and surfaces merge
/// their members into a `B_SPLINE_CURVE_WITH_KNOTS` / `B_SPLINE_SURFACE_WITH_KNOTS`
/// record whose argument layout matches the resolver's expectations (degree /
/// poles / form / closed from the base member, knots + multiplicities from the
/// `_WITH_KNOTS` member, weights from the `RATIONAL_*` member). Other complex
/// entities (unit / context metadata) return an empty type name — the resolver
/// skips them.
pub(super) fn merge_complex_body(body: &str) -> (String, Vec<String>) {
    let members = split_complex_members(body);
    if members.is_empty() {
        return (String::new(), Vec::new());
    }
    // Helper: first member of a given type.
    let member = |t: &str| members.iter().find(|(ty, _)| ty == t);
    // Helper: arg of a member by index.
    let arg = |t: &str, idx: usize| -> Option<String> {
        member(t).and_then(|(_, a)| a.get(idx)).cloned()
    };

    let is_curve = member("B_SPLINE_CURVE").is_some() || member("B_SPLINE_CURVE_WITH_KNOTS").is_some();
    let is_surface = member("B_SPLINE_SURFACE").is_some()
        || member("B_SPLINE_SURFACE_WITH_KNOTS").is_some();

    if is_curve {
        // Complex members carry only the subtype's *own* attributes, no entity
        // name. B_SPLINE_CURVE args: (degree, control_points, curve_form,
        // closed, self_intersect). B_SPLINE_CURVE_WITH_KNOTS adds
        // (multiplicities, knots, knot_spec); RATIONAL_B_SPLINE_CURVE adds
        // (weights).
        let base = member("B_SPLINE_CURVE").or_else(|| member("B_SPLINE_CURVE_WITH_KNOTS"));
        let knots = member("B_SPLINE_CURVE_WITH_KNOTS");
        let rational = member("RATIONAL_B_SPLINE_CURVE");
        let (b_args, k_args, r_args) = match (base, knots, rational) {
            (Some((_, ba)), Some((_, ka)), Some((_, ra))) => (ba, ka, ra),
            (Some((_, ba)), Some((_, ka)), None) => (ba, ka, &Vec::new()),
            (Some((_, ba)), None, Some((_, ra))) => (ba, &Vec::new(), ra),
            (Some((_, ba)), None, None) => (ba, &Vec::new(), &Vec::new()),
            _ => return (String::new(), Vec::new()),
        };
        let degree = b_args.get(0).cloned().unwrap_or_default();
        let control_points = b_args.get(1).cloned().unwrap_or_default();
        let curve_form = b_args.get(2).cloned().unwrap_or_else(|| ".UNSPECIFIED.".to_string());
        let closed = b_args.get(3).cloned().unwrap_or_else(|| ".F.".to_string());
        let self_intersect = b_args.get(4).cloned().unwrap_or_else(|| ".F.".to_string());
        let weights = r_args.first().cloned().unwrap_or_else(|| "SELF".to_string());
        let mults = k_args.get(0).cloned().unwrap_or_else(|| "()".to_string());
        let knots_list = k_args.get(1).cloned().unwrap_or_else(|| "()".to_string());
        let knot_spec = k_args.get(2).cloned().unwrap_or_else(|| ".UNSPECIFIED.".to_string());
        let type_name = if k_args.is_empty() {
            "B_SPLINE_CURVE"
        } else {
            "B_SPLINE_CURVE_WITH_KNOTS"
        };
        // Resolver layout (B_SPLINE_CURVE_WITH_KNOTS) — the reader uses
        // args[1]=degree, [2]=control_points, [3]=weights, [6]=multiplicities,
        // [7]=knots; the remaining slots (curve_form/closed) are carried for
        // fidelity but not read.
        return (
            type_name.to_string(),
            vec![
                "''".to_string(),
                degree,
                control_points,
                weights,
                curve_form,
                closed,
                mults,
                knots_list,
                knot_spec,
            ],
        );
    }

    if is_surface {
        // B_SPLINE_SURFACE args: (u_degree, v_degree, control_points grid,
        // surface_form, closed_u, closed_v, self_intersect).
        // B_SPLINE_SURFACE_WITH_KNOTS adds (u_mults, v_mults, u_knots,
        // v_knots, knot_spec); RATIONAL_B_SPLINE_SURFACE adds (weights).
        let base = member("B_SPLINE_SURFACE").or_else(|| member("B_SPLINE_SURFACE_WITH_KNOTS"));
        let knots = member("B_SPLINE_SURFACE_WITH_KNOTS");
        let rational = member("RATIONAL_B_SPLINE_SURFACE");
        let (b_args, k_args, r_args) = match (base, knots, rational) {
            (Some((_, ba)), Some((_, ka)), Some((_, ra))) => (ba, ka, ra),
            (Some((_, ba)), Some((_, ka)), None) => (ba, ka, &Vec::new()),
            (Some((_, ba)), None, Some((_, ra))) => (ba, &Vec::new(), ra),
            (Some((_, ba)), None, None) => (ba, &Vec::new(), &Vec::new()),
            _ => return (String::new(), Vec::new()),
        };
        let deg_u = b_args.get(0).cloned().unwrap_or_default();
        let deg_v = b_args.get(1).cloned().unwrap_or_default();
        let control_points = b_args.get(2).cloned().unwrap_or_default();
        let surface_form = b_args.get(3).cloned().unwrap_or_else(|| ".UNSPECIFIED.".to_string());
        let closed_u = b_args.get(4).cloned().unwrap_or_else(|| ".F.".to_string());
        let closed_v = b_args.get(5).cloned().unwrap_or_else(|| ".F.".to_string());
        let self_intersect = b_args.get(6).cloned().unwrap_or_else(|| ".F.".to_string());
        let weights = r_args.first().cloned().unwrap_or_else(|| "SELF".to_string());
        // Knot member layout: (u_mults, v_mults, u_knots, v_knots, knot_spec).
        let u_mults = k_args.get(0).cloned().unwrap_or_else(|| "()".to_string());
        let v_mults = k_args.get(1).cloned().unwrap_or_else(|| "()".to_string());
        let u_knots = k_args.get(2).cloned().unwrap_or_else(|| "()".to_string());
        let v_knots = k_args.get(3).cloned().unwrap_or_else(|| "()".to_string());
        let knot_spec = k_args.get(4).cloned().unwrap_or_else(|| ".UNSPECIFIED.".to_string());
        let type_name = if k_args.is_empty() {
            "B_SPLINE_SURFACE"
        } else {
            "B_SPLINE_SURFACE_WITH_KNOTS"
        };
        // Standard ISO 10303-42 layout: name, u_degree, v_degree,
        // control_points grid, surface_form, closed_u, closed_v, self_intersect,
        // u_multiplicities, v_multiplicities, u_knots, v_knots, knot_spec, then
        // the weight grid as the optional 14th arg (rational surfaces).
        return (
            type_name.to_string(),
            vec![
                "''".to_string(),
                deg_u,
                deg_v,
                control_points,
                surface_form,
                closed_u,
                closed_v,
                self_intersect,
                u_mults,
                v_mults,
                u_knots,
                v_knots,
                knot_spec,
                weights,
            ],
        );
    }

    // Non-B-spline complex (unit/context metadata): the resolver skips it.
    (String::new(), Vec::new())
}

/// Split the DATA section into `#id=TYPE(...)` records.
pub(super) fn parse_records(data: &str) -> Result<HashMap<usize, Record>, String> {
    let bytes = data.as_bytes();
    let mut records = HashMap::new();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'#' {
            let mut j = i + 1;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            if j == i + 1 {
                i += 1;
                continue;
            }
            let id: usize = data[i + 1..j]
                .parse()
                .map_err(|_| format!("bad entity id at offset {i}"))?;
            let mut k = j;
            while k < bytes.len() && bytes[k].is_ascii_whitespace() {
                k += 1;
            }
            if k >= bytes.len() || bytes[k] != b'=' {
                return Err(format!("malformed record #{id}: expected '='"));
            }
            let (body, next) = parse_entity_body_text(data, k + 1)?;
            let (type_name, args) = parse_entity_body(&body);
            records.insert(id, Record { type_name, args });
            i = next;
        } else {
            i += 1;
        }
    }
    Ok(records)
}

/// Scan one entity body starting after the `=`, stopping at the terminating
/// `;`. Returns the body text and the index just past the `;`.
pub(super) fn parse_entity_body_text(data: &str, start: usize) -> Result<(String, usize), String> {
    let bytes = data.as_bytes();
    let mut i = start;
    let mut depth = 0usize;
    let mut in_str = false;
    while i < bytes.len() {
        let c = bytes[i] as char;
        if in_str {
            if c == '\'' {
                if i + 1 < bytes.len() && bytes[i + 1] == b'\'' {
                    i += 2;
                    continue;
                }
                in_str = false;
            }
        } else {
            match c {
                '\'' => in_str = true,
                '(' => depth += 1,
                ')' => {
                    if depth == 0 {
                        return Err("unbalanced parens in entity body".into());
                    }
                    depth -= 1;
                    if depth == 0 {
                        let body = data[start..=i].to_string();
                        let mut j = i + 1;
                        // STEP complex entities end `( A() B() C() )()` — the
                        // compound's own (usually empty) attribute list follows
                        // the grouped members. Consume it before the `;`.
                        let mut k = j;
                        while k < bytes.len() && bytes[k].is_ascii_whitespace() {
                            k += 1;
                        }
                        if k + 1 < bytes.len() && bytes[k] == b'(' && bytes[k + 1] == b')' {
                            j = k + 2;
                        }
                        while j < bytes.len() && bytes[j].is_ascii_whitespace() {
                            j += 1;
                        }
                        if j >= bytes.len() || bytes[j] != b';' {
                            return Err("missing ';' after entity body".into());
                        }
                        return Ok((body, j + 1));
                    }
                }
                _ => {}
            }
        }
        i += 1;
    }
    Err("unterminated entity body".into())
}

pub(super) fn parse_ref(s: &str) -> Option<usize> {
    s.trim().strip_prefix('#')?.trim().parse().ok()
}

/// STEP logical: `.T.` is true, `.F.` is false. Missing / other → `default`.
pub(super) fn parse_logical(s: Option<&str>, default: bool) -> bool {
    match s.map(str::trim) {
        Some(".T.") => true,
        Some(".F.") => false,
        _ => default,
    }
}

pub(super) fn parse_ref_list(s: &str) -> Vec<usize> {
    let s = s.trim();
    if !(s.starts_with('(') && s.ends_with(')')) {
        return Vec::new();
    }
    split_top(&s[1..s.len() - 1])
        .into_iter()
        .filter_map(|a| parse_ref(&a))
        .collect()
}

pub(super) fn parse_f64(s: &str) -> Result<f64, String> {
    s.trim()
        .parse()
        .map_err(|_| format!("bad real literal '{s}'"))
}

pub(super) fn parse_xyz(s: &str) -> Result<GpXyz, String> {
    let s = s.trim();
    let inner = s.trim_start_matches('(').trim_end_matches(')');
    let parts: Vec<String> = split_top(inner)
        .into_iter()
        .map(|p| p.trim().to_string())
        .collect();
    if parts.len() != 3 {
        return Err(format!("expected 3-component tuple, got '{s}'"));
    }
    Ok(GpXyz::new(
        parse_f64(&parts[0])?,
        parse_f64(&parts[1])?,
        parse_f64(&parts[2])?,
    ))
}

/// Parse a 2D `(x, y)` tuple (STEP 2D curves use two-component CARTESIAN_POINTs).
pub(super) fn parse_xy(s: &str) -> Result<(f64, f64), String> {
    let s = s.trim();
    let inner = s.trim_start_matches('(').trim_end_matches(')');
    let parts: Vec<String> = split_top(inner)
        .into_iter()
        .map(|p| p.trim().to_string())
        .collect();
    if parts.len() != 2 {
        return Err(format!("expected 2-component tuple, got '{s}'"));
    }
    Ok((parse_f64(&parts[0])?, parse_f64(&parts[1])?))
}

pub(super) fn parse_str(s: &str) -> String {
    let s = s.trim();
    if s.len() >= 2 && s.starts_with('\'') && s.ends_with('\'') {
        s[1..s.len() - 1].replace("''", "'")
    } else {
        s.to_string()
    }
}

/// Reference resolver with memoization; unknown/unsupported entities are
/// recorded as warnings and skipped.
pub(super) struct Resolver<'a> {
    pub(super) records: &'a HashMap<usize, Record>,
    pub(super) b: TopoBuilder,
    pub(super) shape_cache: RefCell<HashMap<usize, TopoShape>>,
    pub(super) point_cache: RefCell<HashMap<usize, GpPnt>>,
    pub(super) dir_cache: RefCell<HashMap<usize, GpDir>>,
    pub(super) axis_cache: RefCell<HashMap<usize, GpAx2>>,
    pub(super) curve_cache: RefCell<HashMap<usize, Arc<dyn Curve>>>,
    pub(super) surface_cache: RefCell<HashMap<usize, Arc<dyn Surface>>>,
    pub(super) curve2d_cache: RefCell<HashMap<usize, (Arc<dyn Curve2d>, (f64, f64))>>,
    /// SURFACE_CURVE / SEAM_CURVE id → its `associated_geometry` (pcurve_or_surface)
    /// reference list, carried for the face-level pcurve association.
    pub(super) surface_curve_pcurves: RefCell<HashMap<usize, Vec<usize>>>,
    /// Edge TShape pointer → the curve entity id of its geometry (an EDGE_CURVE's
    /// `curve` attribute, which may be a SURFACE_CURVE carrying pcurves).
    pub(super) edge_curve_ref: RefCell<HashMap<usize, usize>>,
    pub(super) resolving: RefCell<HashSet<usize>>,
    pub(super) warnings: RefCell<Vec<String>>,
}
