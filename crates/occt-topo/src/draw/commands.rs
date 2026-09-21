use super::prelude::*;
use super::*;


/// `unset <name>` — remove the session variable `name`. Unknown names are an
/// error, matching Tcl's `unset` behaviour; [`eval_number`] on a name that no
/// longer exists falls through to a parse error rather than reading a stale
/// value.
pub(super) fn cmd_unset(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 1, "unset")?;
    if session.vars.remove(&args[0]).is_none() {
        return Err(format!("draw: unset: variable '{}' not found", args[0]));
    }
    session.log.push(format!("unset: removed {}", args[0]));
    Ok(())
}
/// `incr <name> [step]` — add `step` (default 1) to the numeric variable
/// `name` and store the result back. The variable must already hold a number
/// (`set i 0` first). The natural counter for `for` bodies and the scripting
/// layer's integer arithmetic workhorse; port of Tcl's `incr`.
pub(super) fn cmd_incr(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    if args.len() != 1 && args.len() != 2 {
        return Err(format!(
            "draw: incr: expected 1 or 2 arguments, got {}",
            args.len()
        ));
    }
    let cur = session
        .vars
        .get(&args[0])
        .ok_or_else(|| format!("draw: incr: variable '{}' not found", args[0]))?;
    let cur = cur
        .parse::<f64>()
        .map_err(|_| format!("draw: incr: '{}' is not a number", cur))?;
    let step = if args.len() == 2 {
        eval_number(session, &args[1])?
    } else {
        1.0
    };
    let next = cur + step;
    let text = fmt_num(next);
    set_var(session, &args[0], &text);
    session.log.push(format!("incr {} = {text}", args[0]));
    Ok(())
}

/// `for <var> <start> <end> { <body...> }` — run the body once per integer in
/// `[start, end]` inclusive, stepping `−1` when `start > end`.
///
/// Before each iteration `var` is set to the current integer and the body
/// tokens are `$`-expanded again, so `box b$i 1 1 1` builds `b1`, `b2`, ...
/// The body is any Draw command line; braces are optional (`for i 1 3
/// box b$i 1 1 1` works). The loop variable remains set to its final value
/// after the loop. This is the interpreter's only looping construct — a
/// minimal `Draw_Interpretor::Eval` over a variable-binding driver.
pub(super) fn cmd_for(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    if args.len() < 4 {
        return Err(format!(
            "draw: for: expected at least 4 arguments, got {}",
            args.len()
        ));
    }
    let var = &args[0];
    let start = eval_number(session, &args[1])? as i64;
    let end = eval_number(session, &args[2])? as i64;
    let body = command_block(args, "for", 3)?;
    let step = if start <= end { 1 } else { -1 };
    let mut i = start;
    loop {
        set_var(session, var, &fmt_num(i as f64));
        let line = expand_vars(session, body).join(" ");
        execute_line(session, &line)?;
        if i == end || session.stop {
            break;
        }
        i += step;
    }
    Ok(())
}

/// `if <a> <b> <op> { <body...> }` — numeric comparison and conditional
/// execution.
///
/// When `a op b` holds (`op` one of `== != < > <= >=`, in the postfix slot
/// matching `expr`, e.g. `if $x 3 >`) the body is executed as one line,
/// exactly like a `for` body; otherwise it is skipped. Braces are optional.
/// `if` is the interpreter's branch, mirroring `Draw_Interpretor`'s `if` on
/// Tcl expr strings.
pub(super) fn cmd_if(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    if args.len() < 4 {
        return Err(format!(
            "draw: if: expected at least 4 arguments, got {}",
            args.len()
        ));
    }
    let a = eval_number(session, &args[0])?;
    let b = eval_number(session, &args[1])?;
    let holds = match args[2].as_str() {
        "==" => a == b,
        "!=" => a != b,
        "<" => a < b,
        ">" => a > b,
        "<=" => a <= b,
        ">=" => a >= b,
        op => return Err(format!("draw: if: unsupported comparison '{op}'")),
    };
    if !holds {
        return Ok(());
    }
    let body = command_block(args, "if", 3)?;
    let line = expand_vars(session, body).join(" ");
    execute_line(session, &line)
}

// ---------------------------------------------------------------------------
// Shape transforms
// ---------------------------------------------------------------------------

/// `translate <name> <dx> <dy> <dz> [out]` — translate `name`'s registered
/// geometry by `(dx, dy, dz)` and register the result under `out`, or back
/// under `name` when `out` is omitted. Wraps
/// [`translate_shape`](crate::shape_ops::translate_shape), the OCCT
/// `BRep_Tool::Transform` geometry rewrite — world coordinates move.
///
/// With an `out` name the moved geometry is a deep [`translated_copy`], so the
/// source shape is left untouched (a `TopoShape` clone shares its `TShape`, so
/// an in-place transform would move both). Without `out` the transform is
/// applied in place, matching OCCT's `Draw` `translate` with one argument.
pub(super) fn cmd_translate(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    if args.len() != 4 && args.len() != 5 {
        return Err(format!(
            "draw: translate: expected 4 or 5 arguments, got {}",
            args.len()
        ));
    }
    let v = GpVec::new(
        eval_number(session, &args[1])?,
        eval_number(session, &args[2])?,
        eval_number(session, &args[3])?,
    );
    if let Some(out) = args.get(4) {
        let shape = shape_owned(session, &args[0])?;
        let result = translated_copy(&shape, &v)?;
        register_shape(session, out, result);
    } else {
        let mut shape = shape_owned(session, &args[0])?;
        translate_shape(&mut shape, &v)?;
        register_shape(session, &args[0], shape);
    }
    Ok(())
}

/// `rotate <name> <ax> <ay> <az> <angle_deg>` — rotate `name` by `angle_deg`
/// degrees about the axis through the origin with direction `(ax, ay, az)`,
/// in place. Wraps [`rotate_shape`](crate::shape_ops::rotate_shape)
/// (`gp_Trsf::SetRotation` + `BRep_Tool::Transform`); a 90° Z rotation maps
/// the corner `(1,0,0)` to `(0,1,0)`.
pub(super) fn cmd_rotate(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 5, "rotate")?;
    let mut shape = shape_owned(session, &args[0])?;
    let dir = GpDir::new(
        eval_number(session, &args[1])?,
        eval_number(session, &args[2])?,
        eval_number(session, &args[3])?,
    )
    .map_err(|e| format!("draw: rotate: {e}"))?;
    let angle = eval_number(session, &args[4])?.to_radians();
    let axis = GpAx1::new(GpPnt::zero(), dir);
    rotate_shape(&mut shape, &axis, angle)?;
    register_shape(session, &args[0], shape);
    Ok(())
}

/// `scale <name> <factor> [out]` — uniformly scale `name` about the origin by
/// `factor` and register the result under `out`, or back under `name`. A unit
/// box scaled by 2 becomes `[0,2]³`. Wraps
/// [`scale_shape`](crate::shape_ops::scale_shape).
///
/// As with `translate`, the `out` form deep-copies so the source keeps its
/// geometry; the no-`out` form scales in place.
pub(super) fn cmd_scale(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    if args.len() != 2 && args.len() != 3 {
        return Err(format!(
            "draw: scale: expected 2 or 3 arguments, got {}",
            args.len()
        ));
    }
    let factor = eval_number(session, &args[1])?;
    if let Some(out) = args.get(2) {
        let shape = shape_owned(session, &args[0])?;
        let mut t = GpTrsf::identity();
        t.set_scale(&GpPnt::zero(), factor)
            .map_err(|e| format!("draw: scale: {e}"))?;
        let result = transformed_copy(&shape, &t)?;
        register_shape(session, out, result);
    } else {
        let mut shape = shape_owned(session, &args[0])?;
        scale_shape(&mut shape, &GpPnt::zero(), factor)?;
        register_shape(session, &args[0], shape);
    }
    Ok(())
}

/// `copy <name> <out>` — deep-copy `name` into `out`.
///
/// The copy is fully independent: it gets fresh `TShape` handles and its own
/// re-registered geometry via
/// [`translated_copy`](crate::shape_ops::translated_copy) with the identity
/// translation, so later edits to `out` (fillet, translate, ...) never touch
/// `name`. The analogue of OCCT's `copy` / `TopoDS::Transformed` with the
/// identity transform.
pub(super) fn cmd_copy(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 2, "copy")?;
    let shape = shape_owned(session, &args[0])?;
    let copy = translated_copy(&shape, &GpVec::zero())?;
    register_shape(session, &args[1], copy);
    Ok(())
}

/// `transform <name> <a11> .. <a33> <tx> <ty> <tz>` — apply an arbitrary
/// affine map to `name`, in place.
///
/// The nine `a` arguments are the 3×3 linear part read row-major (a11 a12 a13
/// / a21 a22 a23 / a31 a32 a33) and `(tx, ty, tz)` is the translation. The
/// `GpTrsf` is assembled directly with a unit scale and a compound form, then
/// applied with [`transform_shape`](crate::shape_ops::transform_shape) — the
/// OCCT `BRep_Tool::Transform` full affine path. An identity matrix with zero
/// translation is the no-op.
pub(super) fn cmd_transform(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 13, "transform")?;
    let mut m = GpMat::identity();
    for (i, row) in m.m.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = eval_number(session, &args[1 + i * 3 + j])?;
        }
    }
    let t = GpTrsf {
        scale: 1.0,
        shape: TrsfForm::CompoundTrsf,
        matrix: m,
        loc: GpXyz::new(
            eval_number(session, &args[10])?,
            eval_number(session, &args[11])?,
            eval_number(session, &args[12])?,
        ),
    };
    let mut shape = shape_owned(session, &args[0])?;
    transform_shape(&mut shape, &t)?;
    register_shape(session, &args[0], shape);
    Ok(())
}

/// `mirror <name> <nx> <ny> <nz> [out]` — mirror `name` across the plane
/// through the origin whose normal is `(nx, ny, nz)`.
///
/// The reflector matrix is the Householder `I − 2·n·nᵀ` built from the unit
/// normal and applied as a unit-scale `GpTrsf` — the OCCT `gp_Trsf::SetMirror`
/// for the `Ax2` (plane) case. A unit box mirrored across the Z = 0 plane
/// (normal `0 0 1`) maps `[0,1]³` to `[0,1]²×[−1,0]`. The result is registered
/// under `out`, or back under `name` when `out` is omitted.
pub(super) fn cmd_mirror(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    if args.len() != 4 && args.len() != 5 {
        return Err(format!(
            "draw: mirror: expected 4 or 5 arguments, got {}",
            args.len()
        ));
    }
    let (nx, ny, nz) = (
        eval_number(session, &args[1])?,
        eval_number(session, &args[2])?,
        eval_number(session, &args[3])?,
    );
    let norm = (nx * nx + ny * ny + nz * nz).sqrt();
    if norm < 1e-15 {
        return Err("draw: mirror: zero normal".into());
    }
    let (nx, ny, nz) = (nx / norm, ny / norm, nz / norm);
    let m = GpMat::new(
        1.0 - 2.0 * nx * nx,
        -2.0 * nx * ny,
        -2.0 * nx * nz,
        -2.0 * ny * nx,
        1.0 - 2.0 * ny * ny,
        -2.0 * ny * nz,
        -2.0 * nz * nx,
        -2.0 * nz * ny,
        1.0 - 2.0 * nz * nz,
    );
    let t = GpTrsf {
        scale: 1.0,
        shape: TrsfForm::CompoundTrsf,
        matrix: m,
        loc: GpXyz::zero(),
    };
    if let Some(out) = args.get(4) {
        let shape = shape_owned(session, &args[0])?;
        let result = transformed_copy(&shape, &t)?;
        register_shape(session, out, result);
    } else {
        let mut shape = shape_owned(session, &args[0])?;
        transform_shape(&mut shape, &t)?;
        register_shape(session, &args[0], shape);
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// View commands (soft render via viz_scene)
// ---------------------------------------------------------------------------

/// `view <name> [w h] [file]` — soft-render `name` through the session camera
/// and write a shaded binary PPM.
///
/// Defaults to 256×256 and `view.ppm`. The renderer is the [`crate::viz_scene`]
/// ray-cast shaded pipeline ([`render_scene_ppm_shaded`]) with the default
/// `RenderSettings` (one key light, Phong shading); the log records the file
/// and dimensions. The analogue of `V3d_View::Dump` writing `Write_PPM`.
pub(super) fn cmd_view(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    if args.len() < 1 || args.len() > 4 {
        return Err(format!(
            "draw: view: expected 1 to 4 arguments, got {}",
            args.len()
        ));
    }
    let shape = shape_owned(session, &args[0])?;
    let mut w = 256;
    let mut h = 256;
    let mut file = "view.ppm".to_string();
    if args.len() >= 3 {
        w = eval_number(session, &args[1])? as usize;
        h = eval_number(session, &args[2])? as usize;
    }
    if args.len() >= 4 {
        file = args[3].clone();
    }
    if w == 0 || h == 0 {
        return Err(format!("draw: view: invalid dimensions {w}x{h}"));
    }
    let mut scene = VizScene::new();
    scene.add(SceneShape::new(shape));
    let bytes = render_scene_ppm_shaded(
        &scene,
        &session.camera,
        w,
        h,
        EXPORT_DEFLECTION,
        &RenderSettings::default(),
    );
    std::fs::write(&file, &bytes).map_err(|e| format!("draw: view: {e}"))?;
    session.log.push(format!("view: wrote {file} ({w}x{h})"));
    Ok(())
}

/// `view_camera <ex> <ey> <ez> <tx> <ty> <tz>` — set the session camera to a
/// look-at view from `(ex, ey, ez)` toward `(tx, ty, tz)` with `+Y` up.
/// Subsequent `view` commands render through it. The OCCT analogue is
/// `V3d_View::SetViewOrientation` / `Camera::LookAt`.
pub(super) fn cmd_view_camera(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 6, "view_camera")?;
    let eye = GpPnt::new(
        eval_number(session, &args[0])?,
        eval_number(session, &args[1])?,
        eval_number(session, &args[2])?,
    );
    let target = GpPnt::new(
        eval_number(session, &args[3])?,
        eval_number(session, &args[4])?,
        eval_number(session, &args[5])?,
    );
    session.camera = Camera::look_at(eye, target, GpVec::new(0.0, 1.0, 0.0));
    session.log.push(format!(
        "view_camera: eye ({},{},{}) target ({},{},{})",
        eye.x(),
        eye.y(),
        eye.z(),
        target.x(),
        target.y(),
        target.z()
    ));
    Ok(())
}

/// `orbit <dx_deg> <dy_deg>` — orbit the session camera's eye around its
/// target: yaw by `dx_deg` about the world-up axis, then pitch by `dy_deg`
/// (clamped near the poles). Ports `V3d_View::Rotate` for the orbit case.
pub(super) fn cmd_orbit(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 2, "orbit")?;
    let dx = eval_number(session, &args[0])?;
    let dy = eval_number(session, &args[1])?;
    session.camera.camera_orbit(dx, dy);
    session.log.push(format!("orbit: dx {dx} dy {dy}"));
    Ok(())
}

/// `zoom <factor>` — zoom the session camera by moving the eye toward its
/// target; `factor > 1` zooms in. Ports `V3d_View::SetZoom`.
pub(super) fn cmd_zoom(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 1, "zoom")?;
    let factor = eval_number(session, &args[0])?;
    session.camera.camera_zoom(factor);
    session.log.push(format!("zoom: factor {factor}"));
    Ok(())
}

/// `pan <dx> <dy>` — pan the session camera's eye and target together along
/// the camera right/up plane, in screen pixels at a nominal 600-px viewport.
/// Ports `V3d_View::Pan`.
pub(super) fn cmd_pan(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 2, "pan")?;
    let dx = eval_number(session, &args[0])?;
    let dy = eval_number(session, &args[1])?;
    session.camera.camera_pan(dx, dy);
    session.log.push(format!("pan: dx {dx} dy {dy}"));
    Ok(())
}

// ---------------------------------------------------------------------------
// Procedures: `proc` / `def` definitions and `call`
// ---------------------------------------------------------------------------

/// Peel a token into its structural edges: `(opens, core, closes, sep)`.
///
/// The whitespace-only parser glues `-;`, `$x;` and `};` into single tokens, so
/// a body token may carry leading `{`s, a core, trailing `}`s and a trailing
/// `;` all at once. This splits one token into those parts — the `;` is
/// stripped *before* counting trailing `}`s so `};` closes a brace group and
/// then separates commands.
pub(super) fn peel_token(t: &str) -> (usize, &str, usize, bool) {
    let opens = t.chars().take_while(|&c| c == '{').count();
    let rest = &t[opens..];
    let (rest, sep) = match rest.strip_suffix(';') {
        Some(r) => (r, true),
        None => (rest, false),
    };
    let closes = rest.chars().rev().take_while(|&c| c == '}').count();
    let mid = &rest[..rest.len() - closes];
    (opens, mid, closes, sep)
}

/// Split a procedure body token list into command lines on `;`.
///
/// A `;` at brace depth zero separates commands, so `set a 1; expr $a 2 *` is
/// two commands while a `;` inside an `if`/`for` body stays inside its segment.
/// Brace characters and `;` are recognised on token edges (see [`peel_token`]).
/// Segments keep their (now clean) tokens so they round-trip through
/// [`execute_line`].
pub(super) fn split_body(body: &[String]) -> Vec<Vec<String>> {
    let mut segments: Vec<Vec<String>> = Vec::new();
    let mut cur: Vec<String> = Vec::new();
    let mut depth = 0i32;
    for t in body {
        let (opens, mid, closes, sep) = peel_token(t);
        for _ in 0..opens {
            depth += 1;
            cur.push("{".to_string());
        }
        if !mid.is_empty() {
            cur.push(mid.to_string());
        }
        for _ in 0..closes {
            depth -= 1;
            cur.push("}".to_string());
        }
        if sep && depth == 0 && !cur.is_empty() {
            segments.push(std::mem::take(&mut cur));
        }
    }
    if !cur.is_empty() {
        segments.push(cur);
    }
    segments
}

/// `proc <name> { p1 p2 ... } { body }` — define a procedure.
///
/// `def` is an alias. Parameters are bound positionally at call time and the
/// body — a `;`-separated command sequence, stored verbatim — is executed with
/// them in scope. The body may call other procedures (including itself) by
/// name or via `call`. `return [value]` stops the body and stores `value` (or
/// the empty string) in `result`. The `{`...`}` groups are optional: without
/// them the parameter list is empty and the body is the rest of the line.
pub(super) fn cmd_proc(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    if args.len() < 2 {
        return Err("draw: proc: expected <name> { params } { body }".into());
    }
    // The name is the one token expanded here (`proc $name ...`); everything
    // else — parameters and body — stays verbatim so `$param` references are
    // resolved at call time.
    let name = expand_token(session, &args[0]);
    // Optional `{ p1 p2 ... }` parameter group.
    let mut i = 1;
    let mut params = Vec::new();
    if args.get(i).map(|s| s == "{").unwrap_or(false) {
        i += 1;
        while i < args.len() && args[i] != "}" {
            params.push(args[i].clone());
            i += 1;
        }
        if i >= args.len() {
            return Err(format!("draw: proc '{name}': unbalanced parameter braces"));
        }
        i += 1; // skip the closing brace
    }
    // Body: the remaining tokens. A single `{ ... }` group is unwrapped so the
    // stored body is the bare command sequence; otherwise the tail is verbatim.
    // Brace characters are counted on token edges (`};` closes a group), so the
    // whitespace parser gluing `;` to a brace does not confuse the match.
    let mut body = Vec::new();
    if i < args.len() {
        if args[i].starts_with('{') {
            let mut depth = 0i32;
            let mut end = None;
            for (k, t) in args[i..].iter().enumerate() {
                let (opens, _, closes, _) = peel_token(t);
                depth += opens as i32 - closes as i32;
                if depth == 0 {
                    end = Some(i + k);
                    break;
                }
            }
            match end {
                Some(end) => {
                    body.extend_from_slice(&args[i + 1..end]);
                    if end + 1 < args.len() {
                        return Err(format!("draw: proc '{name}': unexpected tokens after body"));
                    }
                }
                None => return Err(format!("draw: proc '{name}': unbalanced body braces")),
            }
        } else {
            body.extend_from_slice(&args[i..]);
        }
    }
    if body.is_empty() {
        return Err(format!("draw: proc '{name}': empty body"));
    }
    session.procs.insert(
        name.clone(),
        ProcDef {
            name: name.clone(),
            params,
            body,
        },
    );
    session.log.push(format!("proc: defined {name}"));
    Ok(())
}

/// `call <name> <arg...>` — invoke a user-defined procedure by name.
///
/// Arguments are already `$`-expanded by the caller's [`execute_line`] and are
/// bound positionally to the procedure's parameters. The procedure's final
/// `result` value (from `expr` or `return`) is visible to the caller through
/// the shared `result` variable.
pub(super) fn cmd_call(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    if args.is_empty() {
        return Err("draw: call: expected a procedure name".into());
    }
    call_proc(session, &args[0], &args[1..])
}

/// `return [value]` — stop the current procedure.
///
/// Inside a procedure body this aborts the remaining commands; a value (or the
/// empty string) is stored in `result` first. At the top level — where there is
/// no procedure to stop — it behaves like `set result <value>`. The control
/// flow is signalled with the [`PROC_RETURN`] sentinel, which [`call_proc`]
/// catches and turns into a normal return.
pub(super) fn cmd_return(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    if args.len() > 1 {
        return Err(format!(
            "draw: return: expected 0 or 1 argument, got {}",
            args.len()
        ));
    }
    if let Some(v) = args.first() {
        set_var(session, "result", v);
        session.log.push(format!("result = {v}"));
    }
    if session.proc_stack.is_empty() {
        Ok(())
    } else {
        Err(PROC_RETURN.to_string())
    }
}

/// Execute a procedure body against the session with the parameters already
/// bound. Each `;`-separated segment is one command line through
/// [`execute_line`]; [`PROC_RETURN`] stops the loop. Errors are wrapped with
/// the procedure name and the current call stack.
pub(super) fn run_proc_body(session: &mut DrawSession, name: &str, body: &[String]) -> Result<(), String> {
    for segment in split_body(body) {
        if segment.is_empty() {
            continue;
        }
        let line = segment.join(" ");
        match execute_line(session, &line) {
            Ok(()) => {}
            Err(e) if e == PROC_RETURN => break,
            Err(e) => {
                let stack = session.proc_stack.join(" -> ");
                return Err(format!(
                    "draw: in proc '{name}' (call stack: {stack}): {e}"
                ));
            }
        }
    }
    Ok(())
}

/// Look up `name`, bind the call arguments to its parameters, run its body,
/// and restore the caller's bindings. Shared by `call <name>` and direct-name
/// invocation from [`dispatch`]. Recursive calls work — each level pushes onto
/// [`DrawSession::proc_stack`] and restores its parameters on the way out.
pub(super) fn call_proc(session: &mut DrawSession, name: &str, args: &[String]) -> Result<(), String> {
    let proc = session
        .procs
        .get(name)
        .cloned()
        .ok_or_else(|| format!("draw: procedure '{name}' not found"))?;
    if args.len() != proc.params.len() {
        return Err(format!(
            "draw: proc '{name}': expected {} argument(s), got {}",
            proc.params.len(),
            args.len()
        ));
    }
    if session.proc_stack.len() >= PROC_MAX_DEPTH {
        return Err(format!(
            "draw: proc '{name}': max recursion depth {PROC_MAX_DEPTH} exceeded"
        ));
    }
    // Bind parameters, saving the previous values so the caller's bindings are
    // restored when the procedure returns (parameters shadow outer variables).
    let mut saved: Vec<Option<String>> = Vec::with_capacity(proc.params.len());
    for (param, arg) in proc.params.iter().zip(args.iter()) {
        saved.push(session.vars.get(param).cloned());
        session.vars.insert(param.clone(), arg.clone());
    }
    session.proc_stack.push(name.to_string());
    let body = proc.body.clone();
    let result = run_proc_body(session, name, &body);
    session.proc_stack.pop();
    for (param, old) in proc.params.iter().zip(saved) {
        match old {
            Some(v) => {
                session.vars.insert(param.clone(), v);
            }
            None => {
                session.vars.remove(param);
            }
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Arrays, lists and strings
// ---------------------------------------------------------------------------

/// `lappend <name> <value...>` — append values to the list variable `name`,
/// creating it if needed.
///
/// A list is a space-separated string, matching Tcl's default representation.
/// Because a variable name may be an array element (`arr(0)`), `lappend
/// arr(0) 5` appends to a single element. `len` counts the elements.
pub(super) fn cmd_lappend(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    if args.len() < 2 {
        return Err(format!(
            "draw: lappend: expected at least 2 arguments, got {}",
            args.len()
        ));
    }
    let name = &args[0];
    let existing = session.vars.get(name).cloned().unwrap_or_default();
    let mut parts: Vec<String> = Vec::new();
    if !existing.is_empty() {
        parts.push(existing);
    }
    parts.extend(args[1..].iter().cloned());
    let joined = parts.join(" ");
    set_var(session, name, &joined);
    session.log.push(format!("lappend {name} = {joined}"));
    Ok(())
}
