use super::super::prelude::*;
use super::super::*;

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

/// A parsed `#N=TYPE(...)` data record.
#[derive(Debug, Clone)]
pub(in crate::step) struct Record {
    pub(in crate::step) type_name: String,
    pub(in crate::step) args: Vec<String>,
    /// Members `(TYPE, args)` of a complex entity `#N=( M1() M2() ... )` in file
    /// order; empty for a simple record. `parse_entity_body` merges the B-spline
    /// families into `type_name` / `args` and leaves every other complex body
    /// with an empty type name, so the unit / context metadata members are only
    /// reachable here (`step_precision`).
    pub(in crate::step) members: Vec<(String, Vec<String>)>,
}

/// Split a comma-separated argument list at the top nesting level, respecting
/// strings (with `''` escapes), parens and brackets.
pub(in crate::step) fn split_top(s: &str) -> Vec<String> {
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
pub(in crate::step) fn parse_entity_body(body: &str) -> (String, Vec<String>) {
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
pub(in crate::step) fn split_complex_members(body: &str) -> Vec<(String, Vec<String>)> {
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
pub(in crate::step) fn merge_complex_body(body: &str) -> (String, Vec<String>) {
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
pub(in crate::step) fn parse_records(data: &str) -> Result<HashMap<usize, Record>, String> {
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
pub(in crate::step) fn parse_entity_body_text(data: &str, start: usize) -> Result<(String, usize), String> {
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
