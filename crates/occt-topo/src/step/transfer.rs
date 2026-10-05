use super::prelude::*;
use super::*;

/// Shared assembly serialization; `splines` selects the geometry classifier.
pub(super) fn write_assembly_inner(a: &StepAssembly, splines: bool) -> Result<String, String> {
    let mut ctx = WriteCtx::new();
    ctx.splines = splines;
    let mut defs: HashMap<String, usize> = HashMap::new();

    // The assembly's own product definition is the root that children hang off.
    let top = esc_str(&a.name);
    let top_product = ctx.w.emit(format!("PRODUCT('{top}','{top}','',(#{}))", ctx.prod_ctx));
    let top_form = ctx.w.emit(format!("PRODUCT_DEFINITION_FORMATION('','',#{top_product})"));
    let top_def = ctx.w.emit(format!("PRODUCT_DEFINITION('','','',#{top_form},#{})", ctx.def_ctx));
    defs.insert(a.name.clone(), top_def);

    for (name, shape) in &a.products {
        let n = esc_str(name);
        let product = ctx.w.emit(format!("PRODUCT('{n}','{n}','',(#{}))", ctx.prod_ctx));
        let formation = ctx.w.emit(format!("PRODUCT_DEFINITION_FORMATION('','',#{product})"));
        let def = ctx.w.emit(format!("PRODUCT_DEFINITION('','','',#{formation},#{})", ctx.def_ctx));
        let items = ctx.emit_top(shape);
        if !items.is_empty() {
            let rep = ctx.w.emit(format!(
                "ADVANCED_BREP_SHAPE_REPRESENTATION('{n}',({}),#{})",
                join_refs(&items),
                ctx.geom_ctx
            ));
            let pds = ctx.w.emit(format!("PRODUCT_DEFINITION_SHAPE('','',#{def})"));
            ctx.w
                .emit(format!("SHAPE_DEFINITION_REPRESENTATION(#{pds},#{rep})"));
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
    /// Members `(TYPE, args)` of a complex entity `#N=( M1() M2() ... )` in file
    /// order; empty for a simple record. `parse_entity_body` merges the B-spline
    /// families into `type_name` / `args` and leaves every other complex body
    /// with an empty type name, so the unit / context metadata members are only
    /// reachable here (`step_precision`).
    pub(super) members: Vec<(String, Vec<String>)>,
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
        // Member type name: identifier until `(`. Whitespace (including a line
        // break) may separate the member name from its argument list, e.g.
        // `RATIONAL_B_SPLINE_SURFACE` newline `((...))` in ATU01038.step - ISO
        // 10303-21 treats the newline as a token separator, so the name still
        // belongs to the parenthesized argument list that follows it.
        let t0 = i;
        while i < bytes.len() && (bytes[i].is_ascii_alphanumeric() || bytes[i] == b'_') {
            i += 1;
        }
        let name_end = i;
        let mut open = i;
        while open < bytes.len() && bytes[open].is_ascii_whitespace() {
            open += 1;
        }
        if open >= bytes.len() || bytes[open] != b'(' {
            // Skip stray tokens (not a TYPE(...) member).
            while i < bytes.len() && !bytes[i].is_ascii_whitespace() {
                i += 1;
            }
            continue;
        }
        let type_name = inner[t0..name_end].to_string();
        i = open;
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
    let _arg = |t: &str, idx: usize| -> Option<String> {
        member(t).and_then(|(_, a)| a.get(idx)).cloned()
    };

    let is_curve = member("B_SPLINE_CURVE").is_some()
        || member("B_SPLINE_CURVE_WITH_KNOTS").is_some()
        || member("BEZIER_CURVE").is_some()
        || member("UNIFORM_CURVE").is_some()
        || member("QUASI_UNIFORM_CURVE").is_some();
    let is_surface = member("B_SPLINE_SURFACE").is_some()
        || member("B_SPLINE_SURFACE_WITH_KNOTS").is_some()
        || member("BEZIER_SURFACE").is_some()
        || member("UNIFORM_SURFACE").is_some()
        || member("QUASI_UNIFORM_SURFACE").is_some();

    if is_curve {
        // Complex members carry only the subtype's *own* attributes, no entity
        // name. B_SPLINE_CURVE args: (degree, control_points, curve_form,
        // closed, self_intersect). B_SPLINE_CURVE_WITH_KNOTS adds
        // (multiplicities, knots, knot_spec); RATIONAL_B_SPLINE_CURVE adds
        // (weights). The Bezier / uniform / quasi-uniform families have the same
        // base attribute list (`StepToGeom.cxx:295-458` maps them onto a
        // B-spline with synthesized knots, which the resolver's arms do).
        if let Some((ty, ba)) = ["BEZIER_CURVE", "UNIFORM_CURVE", "QUASI_UNIFORM_CURVE"]
            .iter()
            .find_map(|t| member(t))
        {
            let r_args = member("RATIONAL_B_SPLINE_CURVE")
                .map(|(_, a)| a.clone())
                .unwrap_or_default();
            let degree = ba.get(0).cloned().unwrap_or_default();
            let control_points = ba.get(1).cloned().unwrap_or_default();
            let curve_form = ba.get(2).cloned().unwrap_or_else(|| ".UNSPECIFIED.".to_string());
            let closed = ba.get(3).cloned().unwrap_or_else(|| ".F.".to_string());
            let self_intersect = ba.get(4).cloned().unwrap_or_else(|| ".F.".to_string());
            let weights = r_args.first().cloned().unwrap_or_else(|| "SELF".to_string());
            return (
                ty.clone(),
                vec![
                    "''".to_string(),
                    degree,
                    control_points,
                    weights,
                    curve_form,
                    closed,
                    self_intersect,
                ],
            );
        }
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
        let _self_intersect = b_args.get(4).cloned().unwrap_or_else(|| ".F.".to_string());
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
        // The Bezier / uniform / quasi-uniform families share the base attribute
        // list (`StepToGeom.cxx:522-743` converts them to a B-spline with
        // synthesized knots, which the resolver's arms do) and take the weight
        // grid from `RATIONAL_B_SPLINE_SURFACE`.
        if let Some((ty, ba)) = ["BEZIER_SURFACE", "UNIFORM_SURFACE", "QUASI_UNIFORM_SURFACE"]
            .iter()
            .find_map(|t| member(t))
        {
            let r_args = member("RATIONAL_B_SPLINE_SURFACE")
                .map(|(_, a)| a.clone())
                .unwrap_or_default();
            let deg_u = ba.get(0).cloned().unwrap_or_default();
            let deg_v = ba.get(1).cloned().unwrap_or_default();
            let control_points = ba.get(2).cloned().unwrap_or_default();
            let surface_form = ba.get(3).cloned().unwrap_or_else(|| ".UNSPECIFIED.".to_string());
            let closed_u = ba.get(4).cloned().unwrap_or_else(|| ".F.".to_string());
            let closed_v = ba.get(5).cloned().unwrap_or_else(|| ".F.".to_string());
            let self_intersect = ba.get(6).cloned().unwrap_or_else(|| ".F.".to_string());
            let weights = r_args.first().cloned().unwrap_or_else(|| "SELF".to_string());
            return (
                ty.clone(),
                vec![
                    "''".to_string(),
                    deg_u,
                    deg_v,
                    control_points,
                    weights,
                    surface_form,
                    closed_u,
                    closed_v,
                    self_intersect,
                ],
            );
        }
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
            // A complex entity whose members are not a B-spline family comes back
            // with an empty type name, and `parse_entity_body` drops the members
            // it could not merge. Keep them so the unit / context metadata
            // (`LENGTH_UNIT`, `SI_UNIT`, `GLOBAL_UNIT_ASSIGNED_CONTEXT`,
            // `UNCERTAINTY_MEASURE_WITH_UNIT`'s unit) stays readable.
            let members = if type_name.is_empty() && body.trim_start().starts_with('(') {
                split_complex_members(&body)
            } else {
                Vec::new()
            };
            records.insert(
                id,
                Record {
                    type_name,
                    args,
                    members,
                },
            );
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

// ---------------------------------------------------------------------------
// Read precision: `STEPControl_ActorRead::myPrecision`
// ---------------------------------------------------------------------------

/// `Interface_StaticStandards.cxx:33-37`: `read.precision.mode` is initialised
/// `'e' ""` then `SetIVal(..., 0)`, so the default mode is 0.
const READ_PRECISION_MODE: i32 = 0;
/// `Interface_StaticStandards.cxx:39`: `read.precision.val 'r' "1.e-03"`.
const READ_PRECISION_VAL: f64 = 1.0e-3;
/// `StepData_Factors.cxx:24` `myCascadeUnit(1.)` - "length unit for current
/// transfer process (mm by default)".
const CASCADE_UNIT: f64 = 1.0;

/// `STEPConstruct_UnitContext::ConvertSiPrefix` (`STEPConstruct_UnitContext.cxx:210-250`).
/// An unknown prefix falls through to `1.` (`cxx:246-249`).
fn convert_si_prefix(prefix: &str) -> f64 {
    match prefix.trim().trim_matches('.').to_ascii_uppercase().as_str() {
        "EXA" => 1.0e18,
        "PETA" => 1.0e15,
        "TERA" => 1.0e12,
        "GIGA" => 1.0e9,
        "MEGA" => 1.0e6,
        "KILO" => 1.0e3,
        "HECTO" => 1.0e2,
        "DECA" => 1.0e1,
        "DECI" => 1.0e-1,
        "CENTI" => 1.0e-2,
        "MILLI" => 1.0e-3,
        "MICRO" => 1.0e-6,
        "NANO" => 1.0e-9,
        "PICO" => 1.0e-12,
        "FEMTO" => 1.0e-15,
        "ATTO" => 1.0e-18,
        _ => 1.0,
    }
}

/// `STEPConstruct_UnitContext::ComputeFactors` (`cxx:306-450`) for a unit whose
/// complex body is a `ConversionBasedUnitAndLengthUnit` or a
/// `SiUnitAndLengthUnit`: `lengthFactor = parameter * 1000. / aCascadeUnit`
/// (`cxx:424-428`, the `METER` build option is not defined, so the scaled
/// branch is the live one). Returns `None` when the members are not a length
/// unit.
///
/// `cxx:322-381` (conversion based, e.g. INCH): `theFactor = theSIPFactor *
/// theMVAL`, the conversion measure times the SI prefix of the unit it converts
/// to (`INCH` -> `LENGTH_MEASURE(25.4)` on `.MILLI.METRE.` -> `25.4`).
/// `cxx:384-417` (SI): `theFactor = theSIPFactor * theSIUNF` (the SI name
/// factor is `1.`). `cxx:344-348` / `cxx:370`: a missing conversion factor or a
/// non-SI target yields no parameter (`ComputeFactors` returns `-1` / `3` and
/// the pre-set default stays in place).
fn length_unit_factor(
    members: &[(String, Vec<String>)],
    records: &HashMap<usize, Record>,
) -> Option<f64> {
    if !members.iter().any(|(t, _)| t == "LENGTH_UNIT") {
        return None;
    }
    let parameter =
        if let Some((_, args)) = members.iter().find(|(t, _)| t == "CONVERSION_BASED_UNIT") {
            let conv = args.get(1).and_then(|a| parse_ref(a))?;
            let conv = records.get(&conv)?;
            let value = conv.args.first().and_then(|a| measure_value(a))?;
            let target = conv.args.get(1).and_then(|a| parse_ref(a))?;
            si_unit_member_scale(&records.get(&target)?.members)? * value
        } else {
            si_unit_member_scale(members)?
        };
    Some(parameter * 1000.0 / CASCADE_UNIT)
}

/// `STEPConstruct_UnitContext::ComputeFactors` (`cxx:384-417`) for the
/// `SI_UNIT(prefix, name)` member of a unit's complex body: `theFactor =
/// theSIPFactor * theSIUNF` (`cxx:405`). `SiUnitNameFactor` (`cxx:254-268`)
/// assigns `theSIUNFactor = 1.` before its switch, so `theSIUNF` is `1.` for
/// every name; an unrecognised name only sets `status = 11` at `cxx:397-400`
/// (`cxx:361-365` for the conversion-based arm) and computation continues, so
/// the prefix factor is returned for any name. `None` when the unit has no
/// `SI_UNIT` member at all, which is the `return 3` branch of `cxx:368-371`.
fn si_unit_member_scale(members: &[(String, Vec<String>)]) -> Option<f64> {
    let (_, args) = members.iter().find(|(t, _)| t == "SI_UNIT")?;
    let prefix = args.first().map(String::as_str).unwrap_or("$");
    Some(if prefix.trim() == "$" {
        1.0
    } else {
        convert_si_prefix(prefix)
    })
}

/// `STEPConstruct_UnitContext::ComputeFactors` (`cxx:306-450`) for a unit whose
/// complex body is a `ConversionBasedUnitAndPlaneAngleUnit` or a
/// `SiUnitAndPlaneAngleUnit` (`cxx:439-444`): the returned parameter is the
/// `planeAngleFactor`. `None` when the unit is not a plane angle unit.
///
/// `cxx:322-381` (conversion based): `theFactor = theSIPFactor * theMVAL`, the
/// conversion measure times the SI prefix of the unit it converts to
/// (`DEGREE` -> `PLANE_ANGLE_MEASURE(0.01745329252)` on `.RADIAN.`).
/// `cxx:344-348` / `cxx:370`: a missing conversion factor or a non-SI target
/// yields no parameter (`ComputeFactors` returns `-1` / `3` and the pre-set
/// default stays in place).
fn plane_angle_unit_factor(
    members: &[(String, Vec<String>)],
    records: &HashMap<usize, Record>,
) -> Option<f64> {
    if !members.iter().any(|(t, _)| t == "PLANE_ANGLE_UNIT") {
        return None;
    }
    if let Some((_, args)) = members.iter().find(|(t, _)| t == "CONVERSION_BASED_UNIT") {
        let conv = args.get(1).and_then(|a| parse_ref(a))?;
        let conv = records.get(&conv)?;
        let value = conv.args.first().and_then(|a| measure_value(a))?;
        let target = conv.args.get(1).and_then(|a| parse_ref(a))?;
        return Some(si_unit_member_scale(&records.get(&target)?.members)? * value);
    }
    si_unit_member_scale(members)
}

/// `STEPControl_ActorRead::PrepareUnits` -> `STEPConstruct_UnitContext::
/// ComputeFactors(theGUAC, ...)` (`STEPControl_ActorRead.cxx:2312-2355`): the
/// plane angle factor of the `GLOBAL_UNIT_ASSIGNED_CONTEXT`. It is the value
/// `StepData_Factors::PlaneAngleFactor()` hands to
/// `StepToGeom::MakeConicalSurface` (`StepToGeom.cxx:1316`) and the value
/// `FactorDegreeRadian()` hands to `GeomConvert_Units::DegreeToRadian`
/// (`StepToTopoDS_TranslateEdge.cxx:579`).
///
/// `ComputeFactors` pre-sets `planeAngleFactor = PI/180.` (`cxx:281`) and leaves
/// it when the context lists no plane angle unit, so a context that parses but
/// has no angle unit keeps `PI/180.`. A context whose unit list does not resolve
/// is the "Bad RepresentationContext, default unit taken" branch
/// (`STEPControl_ActorRead.cxx:2300-2305`), where `ResetUnits` (`cxx:2400`)
/// leaves every factor at `1.`. The port resolves one flat record map instead of
/// per-representation factors, so the first unit context in the file is used.
pub(super) fn context_plane_angle_factor(records: &HashMap<usize, Record>) -> f64 {
    let mut usable_context = false;
    for rec in records.values() {
        let Some((_, units)) = rec
            .members
            .iter()
            .find(|(t, _)| t == "GLOBAL_UNIT_ASSIGNED_CONTEXT")
        else {
            continue;
        };
        for unit_id in units.iter().flat_map(|a| parse_ref_list(a)) {
            let Some(unit) = records.get(&unit_id) else {
                continue;
            };
            usable_context = true;
            if let Some(f) = plane_angle_unit_factor(&unit.members, records) {
                return f;
            }
        }
    }
    if usable_context {
        PI / 180.0
    } else {
        1.0
    }
}

/// `STEPConstruct_UnitContext::ComputeFactors(theGUAC, ...)` (`cxx:272-302`):
/// the length factor comes from the length unit listed by the representation's
/// `GLOBAL_UNIT_ASSIGNED_CONTEXT`. `lengthFactor` is pre-set to `1.` (`cxx:280`).
pub(super) fn context_length_factor(records: &HashMap<usize, Record>) -> f64 {
    for rec in records.values() {
        let Some((_, units)) = rec
            .members
            .iter()
            .find(|(t, _)| t == "GLOBAL_UNIT_ASSIGNED_CONTEXT")
        else {
            continue;
        };
        for unit_id in units.iter().flat_map(|a| parse_ref_list(a)) {
            if let Some(unit) = records.get(&unit_id) {
                if let Some(f) = length_unit_factor(&unit.members, records) {
                    return f;
                }
            }
        }
    }
    1.0
}

/// `LENGTH_MEASURE(2.E-005)` -> `2.E-005`.
fn measure_value(arg: &str) -> Option<f64> {
    let open = arg.find('(')?;
    let close = arg.rfind(')')?;
    arg[open + 1..close].trim().parse().ok()
}

/// `STEPControl_ActorRead::myPrecision` (`STEPControl_ActorRead.cxx:2370-2384`),
/// the value `StepToTopoDS_TranslateEdgeLoop::Precision()` returns (`cxx:236`)
/// and that `ShapeFix_EdgeProjAux::Compute` (`cxx:844`) and
/// `XSAlgo_ShapeProcessor::CheckPCurve` (`cxx:875` -> `cxx:175`) receive.
///
///   if (ReadPrecisionMode == 1)        myPrecision = ReadPrecisionVal;
///   else if (myUnit.HasUncertainty())  myPrecision = Uncertainty() * LengthFactor();
///   else                               myPrecision = ReadPrecisionVal;
///
/// `myUnit` is the local `STEPConstruct_UnitContext` filled by
/// `STEPControl_ActorRead::PrepareUnits` (`cxx:2312-2384`): `ComputeFactors`
/// (`cxx:2347`) sets the length factor and `ComputeTolerance` (`cxx:2363`) the
/// uncertainty.
pub(super) fn step_precision(records: &HashMap<usize, Record>) -> f64 {
    if READ_PRECISION_MODE == 1 {
        return READ_PRECISION_VAL;
    }
    // `STEPConstruct_UnitContext::ComputeTolerance` (`cxx:480-541`): only a
    // `UNCERTAINTY_MEASURE_WITH_UNIT` whose unit component is an SI (`cxx:504`)
    // or conversion-based (`cxx:521`) length unit counts, and `theUncertainty`
    // (initialised to `RealLast`) is only ever lowered.
    let mut the_uncertainty = f64::MAX;
    let mut has_uncertainty = false;
    for rec in records.values() {
        if rec.type_name != "UNCERTAINTY_MEASURE_WITH_UNIT" {
            continue;
        }
        let Some(value) = rec.args.first().and_then(|a| measure_value(a)) else {
            continue;
        };
        let Some(unit_id) = rec.args.get(1).and_then(|a| parse_ref(a)) else {
            continue;
        };
        let Some(unit) = records.get(&unit_id) else {
            continue;
        };
        if !unit.members.iter().any(|(t, _)| t == "LENGTH_UNIT") {
            continue;
        }
        if the_uncertainty > value {
            the_uncertainty = value;
        }
        has_uncertainty = true;
    }
    if !has_uncertainty {
        return READ_PRECISION_VAL;
    }
    the_uncertainty * context_length_factor(records)
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
    /// `STEPControl_ActorRead::myPrecision` for this file
    /// (`STEPControl_ActorRead.cxx:2370-2384`). Reaches the pcurve post-pass as
    /// `StepToTopoDS_TranslateEdgeLoop::Precision()` (`cxx:236`, `cxx:875`).
    pub(super) precision: f64,
    /// `StepData_Factors::PlaneAngleFactor()` (`StepControl_ActorRead.cxx:2347-2355`,
    /// `StepConstruct_UnitContext.cxx:439-444`): the file's plane angle unit
    /// expressed in radians. Reaches `StepToGeom::MakeConicalSurface`
    /// (`StepToGeom.cxx:1316`) and `GeomConvert_Units::DegreeToRadian`
    /// (`StepToTopoDS_TranslateEdge.cxx:579`).
    pub(super) plane_angle_factor: f64,
    /// `StepData_Factors::LengthFactor()` (`StepControl_ActorRead.cxx:2347-2355`,
    /// `StepConstruct_UnitContext.cxx:424-428`): the file's length unit expressed
    /// in the cascade unit (millimetres). `1.` for a millimetre file. Reaches
    /// every `StepToGeom::Make*` that turns a STEP coordinate, vector magnitude,
    /// radius or distance into geometry (`StepToGeom.cxx:1179`, `:1222`, `:1315`,
    /// `:1452`, `:1549`, `:1616`, `:1706`, `:1899`, `:1947`, `:2111`, `:2577`) and
    /// `GeomConvert_Units::DegreeToRadian` (`StepToTopoDS_TranslateEdge.cxx:573-580`).
    pub(super) length_factor: f64,
    pub(super) b: TopoBuilder,
    pub(super) shape_cache: RefCell<HashMap<usize, TopoShape>>,
    /// `StepToTopoDS_TranslateTool::Bind` for vertices: a `VERTEX_POINT` record
    /// that another record was bound onto resolves to the bound `TopoDS_Vertex`
    /// (`StepToTopoDS_TranslateEdgeLoop.cxx:384-396`, `:466-477`).
    pub(super) vertex_bind: RefCell<HashMap<usize, TopoShape>>,
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
    /// `ShapeFix_Shape::myContext` (`ShapeFix_Shape.cxx:75`): the
    /// `ShapeBuild_ReShape` handle shared by every fix tool of one `FixShape`
    /// pass. `FromSTEP.exec.op` is exactly `FixShape`
    /// (`STEPControl_Controller.cxx:201`) and
    /// `ShapeProcess_OperLibrary.cxx:830` runs it once per transferred root,
    /// so one context covers a whole root shape. `ShapeFix_Shape::Perform`
    /// materialises it with `myResult = Context()->Apply(S)`
    /// (`ShapeFix_Shape.cxx:257`, `ShapeFix_Shell.cxx:139`), which is what
    /// rewrites a face whose *neighbour* was split (T0M face 1693 vs 1695).
    /// The stack is pushed per root in `resolve_shell` / `resolve_solid` and
    /// popped together with that `Apply`.
    pub(super) heal_context: RefCell<Vec<crate::shape_fix_compose_shell::SharedReShape>>,
}
