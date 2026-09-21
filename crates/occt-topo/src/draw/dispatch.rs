use super::prelude::*;
use super::*;


/// `len <name>` — count the elements of the list variable `name` and store the
/// count in `result`. An undefined or empty variable is length 0.
pub(super) fn cmd_len(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 1, "len")?;
    let value = session.vars.get(&args[0]).cloned().unwrap_or_default();
    let n = if value.is_empty() {
        0
    } else {
        value.split_whitespace().count()
    };
    let text = fmt_num(n as f64);
    set_var(session, "result", &text);
    session.log.push(format!("len {} = {n}", args[0]));
    Ok(())
}
/// `concat <a> <b> <out>` — concatenate the list variables `a` and `b` into
/// `out` (Tcl `concat` semantics: space-joined, empties skipped). For plain
/// string concatenation use `$` expansion: `set x $a$b`.
pub(super) fn cmd_concat(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 3, "concat")?;
    let a = session.vars.get(&args[0]).cloned().unwrap_or_default();
    let b = session.vars.get(&args[1]).cloned().unwrap_or_default();
    let joined = [a, b]
        .iter()
        .filter(|s| !s.is_empty())
        .cloned()
        .collect::<Vec<_>>()
        .join(" ");
    set_var(session, &args[2], &joined);
    session.log.push(format!("concat {} = {joined}", args[2]));
    Ok(())
}

// ---------------------------------------------------------------------------
// Command registry
// ---------------------------------------------------------------------------

/// One entry in the centralized command table.
pub(super) struct CommandEntry {
    /// The command name (also the dispatch key).
    pub(super) name: &'static str,
    /// The uniform handler: session + argument tokens.
    pub(super) handler: fn(&mut DrawSession, &[String]) -> Result<(), String>,
    /// One-line description shown by `help`.
    pub(super) summary: &'static str,
    /// Inclusive accepted argument-count range; `usize::MAX` means "or more".
    pub(super) min_args: usize,
    pub(super) max_args: usize,
}

/// Human-readable arity for a table entry: `3`, `2..3`, or `4 or more`.
pub(super) fn arity_text(e: &CommandEntry) -> String {
    if e.max_args == usize::MAX {
        format!("{} or more", e.min_args)
    } else if e.min_args == e.max_args {
        format!("{}", e.min_args)
    } else {
        format!("{}..{}", e.min_args, e.max_args)
    }
}

// Adapt the boolean and writer handlers (which take a fixed `op`/`kind`) and
// the `exit`/`quit` pair to the uniform `(session, args)` table signature.
pub(super) fn cmd_fuse(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    cmd_boolean(session, args, BoolOp::Fuse)
}
pub(super) fn cmd_cut(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    cmd_boolean(session, args, BoolOp::Cut)
}
pub(super) fn cmd_common(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    cmd_boolean(session, args, BoolOp::Common)
}
pub(super) fn cmd_step(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    cmd_write(session, args, WriteKind::Step)
}
pub(super) fn cmd_obj(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    cmd_write(session, args, WriteKind::Obj)
}
pub(super) fn cmd_iges(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    cmd_write(session, args, WriteKind::Iges)
}
pub(super) fn cmd_stl(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    cmd_write(session, args, WriteKind::Stl)
}
pub(super) fn cmd_exit(session: &mut DrawSession, _args: &[String]) -> Result<(), String> {
    session.stop = true;
    Ok(())
}
pub(super) fn cmd_quit(session: &mut DrawSession, _args: &[String]) -> Result<(), String> {
    session.stop = true;
    Ok(())
}

/// `commands` — list every command name (built-ins plus user procedures),
/// sorted, on one log line.
pub(super) fn cmd_commands(session: &mut DrawSession, _args: &[String]) -> Result<(), String> {
    let mut names: Vec<&str> = COMMANDS.iter().map(|e| e.name).collect();
    for name in session.procs.keys() {
        names.push(name);
    }
    names.sort();
    session.log.push(format!("commands: {}", names.join(" ")));
    Ok(())
}

/// `arity <cmd>` — report the accepted argument count of a built-in command or
/// the parameter count of a user procedure.
pub(super) fn cmd_arity(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 1, "arity")?;
    let name = &args[0];
    if let Some(proc) = session.procs.get(name) {
        let n = proc.params.len();
        session.log.push(format!("arity {name}: {n}"));
        return Ok(());
    }
    match COMMANDS.iter().find(|e| e.name == name.as_str()) {
        Some(e) => {
            session.log.push(format!("arity {name}: {}", arity_text(e)));
            Ok(())
        }
        None => Err(format!("draw: arity: unknown command '{name}'")),
    }
}

/// The centralized command table: every dispatchable command, its handler, a
/// one-line description, and its accepted argument-count range. `help`,
/// `commands` and `arity` read from here, and [`dispatch`] looks up handlers
/// here — the fixed `match` table of earlier phases, now data.
pub(super) const COMMANDS: &[CommandEntry] = &[
    CommandEntry { name: "box", handler: cmd_box, summary: "build an axis-aligned box [0,dx]x[0,dy]x[0,dz]", min_args: 4, max_args: 4 },
    CommandEntry { name: "cylinder", handler: cmd_cylinder, summary: "build a Z-axis cylinder of radius r, height h", min_args: 3, max_args: 3 },
    CommandEntry { name: "sphere", handler: cmd_sphere, summary: "build a sphere of radius r centred at the origin", min_args: 2, max_args: 2 },
    CommandEntry { name: "cone", handler: cmd_cone, summary: "build a cone (base r1, top 0), height h", min_args: 4, max_args: 4 },
    CommandEntry { name: "torus", handler: cmd_torus, summary: "build a torus, major radius r1, minor r2", min_args: 3, max_args: 3 },
    CommandEntry { name: "fuse", handler: cmd_fuse, summary: "boolean union a u b -> out", min_args: 3, max_args: 3 },
    CommandEntry { name: "cut", handler: cmd_cut, summary: "boolean difference a - b -> out", min_args: 3, max_args: 3 },
    CommandEntry { name: "common", handler: cmd_common, summary: "boolean intersection a n b -> out", min_args: 3, max_args: 3 },
    CommandEntry { name: "fillet", handler: cmd_fillet, summary: "blend an edge of a solid by radius", min_args: 2, max_args: 3 },
    CommandEntry { name: "offset", handler: cmd_offset, summary: "offset a shell/solid by dist -> out", min_args: 3, max_args: 3 },
    CommandEntry { name: "mesh", handler: cmd_mesh, summary: "tessellate and report the triangle count", min_args: 2, max_args: 2 },
    CommandEntry { name: "step", handler: cmd_step, summary: "write a STEP file", min_args: 2, max_args: 2 },
    CommandEntry { name: "obj", handler: cmd_obj, summary: "write an OBJ mesh file", min_args: 2, max_args: 2 },
    CommandEntry { name: "iges", handler: cmd_iges, summary: "write an IGES file", min_args: 2, max_args: 2 },
    CommandEntry { name: "stl", handler: cmd_stl, summary: "write a binary STL file", min_args: 2, max_args: 2 },
    CommandEntry { name: "info", handler: cmd_info, summary: "print vertex/edge/face counts + volume", min_args: 1, max_args: 1 },
    CommandEntry { name: "bbox", handler: cmd_bbox, summary: "print the axis-aligned bounding-box corners", min_args: 1, max_args: 1 },
    CommandEntry { name: "verts", handler: cmd_verts, summary: "print the distinct vertex coordinates", min_args: 1, max_args: 1 },
    CommandEntry { name: "ls", handler: cmd_ls, summary: "list session shape names", min_args: 0, max_args: 0 },
    CommandEntry { name: "rm", handler: cmd_rm, summary: "remove a named shape", min_args: 1, max_args: 1 },
    CommandEntry { name: "clear", handler: cmd_clear, summary: "empty the session", min_args: 0, max_args: 0 },
    CommandEntry { name: "help", handler: cmd_help, summary: "list commands with arity and description", min_args: 0, max_args: 0 },
    CommandEntry { name: "commands", handler: cmd_commands, summary: "list all command names", min_args: 0, max_args: 0 },
    CommandEntry { name: "arity", handler: cmd_arity, summary: "report a command's accepted argument count", min_args: 1, max_args: 1 },
    CommandEntry { name: "exit", handler: cmd_exit, summary: "set the session stop flag", min_args: 0, max_args: 0 },
    CommandEntry { name: "quit", handler: cmd_quit, summary: "set the session stop flag", min_args: 0, max_args: 0 },
    CommandEntry { name: "set", handler: cmd_set, summary: "store a string/number variable", min_args: 2, max_args: 2 },
    CommandEntry { name: "expr", handler: cmd_expr, summary: "evaluate a op b (postfix) into result", min_args: 3, max_args: 3 },
    CommandEntry { name: "echo", handler: cmd_echo, summary: "append expanded text to the log", min_args: 0, max_args: usize::MAX },
    CommandEntry { name: "vars", handler: cmd_vars, summary: "list session variables", min_args: 0, max_args: 0 },
    CommandEntry { name: "unset", handler: cmd_unset, summary: "remove a variable", min_args: 1, max_args: 1 },
    CommandEntry { name: "incr", handler: cmd_incr, summary: "add step (default 1) to a number var", min_args: 1, max_args: 2 },
    CommandEntry { name: "for", handler: cmd_for, summary: "run a body per integer in [start, end]", min_args: 4, max_args: usize::MAX },
    CommandEntry { name: "if", handler: cmd_if, summary: "run a body when the comparison holds", min_args: 4, max_args: usize::MAX },
    CommandEntry { name: "proc", handler: cmd_proc, summary: "define a procedure with parameters and body", min_args: 2, max_args: usize::MAX },
    CommandEntry { name: "def", handler: cmd_proc, summary: "define a procedure (alias of proc)", min_args: 2, max_args: usize::MAX },
    CommandEntry { name: "call", handler: cmd_call, summary: "call a procedure by name with arguments", min_args: 1, max_args: usize::MAX },
    CommandEntry { name: "return", handler: cmd_return, summary: "stop the current procedure, optionally setting result", min_args: 0, max_args: 1 },
    CommandEntry { name: "lappend", handler: cmd_lappend, summary: "append values to a list/array variable", min_args: 2, max_args: usize::MAX },
    CommandEntry { name: "len", handler: cmd_len, summary: "count elements of a list variable into result", min_args: 1, max_args: 1 },
    CommandEntry { name: "concat", handler: cmd_concat, summary: "concatenate two list variables into a third", min_args: 3, max_args: 3 },
    CommandEntry { name: "translate", handler: cmd_translate, summary: "translate a shape by (dx,dy,dz) [out]", min_args: 4, max_args: 5 },
    CommandEntry { name: "rotate", handler: cmd_rotate, summary: "rotate a shape about an axis by degrees", min_args: 5, max_args: 5 },
    CommandEntry { name: "scale", handler: cmd_scale, summary: "scale a shape by a factor [out]", min_args: 2, max_args: 3 },
    CommandEntry { name: "mirror", handler: cmd_mirror, summary: "mirror a shape across a plane [out]", min_args: 4, max_args: 5 },
    CommandEntry { name: "copy", handler: cmd_copy, summary: "independent deep copy", min_args: 2, max_args: 2 },
    CommandEntry { name: "transform", handler: cmd_transform, summary: "arbitrary affine map (3x3 + translation)", min_args: 13, max_args: 13 },
    CommandEntry { name: "view", handler: cmd_view, summary: "soft-render a shape to a PPM file", min_args: 1, max_args: 4 },
    CommandEntry { name: "view_camera", handler: cmd_view_camera, summary: "set the look-at camera", min_args: 6, max_args: 6 },
    CommandEntry { name: "orbit", handler: cmd_orbit, summary: "orbit the camera about its target", min_args: 2, max_args: 2 },
    CommandEntry { name: "zoom", handler: cmd_zoom, summary: "zoom the camera toward its target", min_args: 1, max_args: 1 },
    CommandEntry { name: "pan", handler: cmd_pan, summary: "pan the camera in screen pixels", min_args: 2, max_args: 2 },
];

// ---------------------------------------------------------------------------
// Inspection
// ---------------------------------------------------------------------------

/// `bbox <name>` — print the axis-aligned bounding box of `name` as two corner
/// points, taken from the registered geometry (so transformed shapes report
/// their world bounds).
pub(super) fn cmd_bbox(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 1, "bbox")?;
    let shape = shape_owned(session, &args[0])?;
    let bb = crate::bbox_from_geometry::shape_bbox(&shape);
    if let Some((x0, x1, y0, y1, z0, z1)) = bb.get() {
        session.log.push(format!(
            "bbox {}: min ({x0}, {y0}, {z0}) max ({x1}, {y1}, {z1})",
            args[0]
        ));
    } else {
        session.log.push(format!("bbox {}: empty", args[0]));
    }
    Ok(())
}

/// `verts <name>` — print the distinct vertex coordinates of `name`, one per
/// log line, using the registered geometry (world points after a transform).
/// The OCCT analogue is `vertices` / `vstats` in the standard Draw plugin.
pub(super) fn cmd_verts(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 1, "verts")?;
    let shape = shape_owned(session, &args[0])?;
    let verts = crate::topo_tools_full::vertices_of(&shape);
    session.log.push(format!("verts {}: {} vertices", args[0], verts.len()));
    for v in verts {
        let p = crate::brep_tool::BRepTool::vertex_point(&v);
        session.log.push(format!("  ({:.6}, {:.6}, {:.6})", p.x(), p.y(), p.z()));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Dispatch
// ---------------------------------------------------------------------------

/// Parse `line` into a command and dispatch it, recording the error in
/// `session.last_error` when it fails.
///
/// This is the entry point the batch and REPL drivers call for every line.
/// Blank lines are ignored; `exit`/`quit` set the session stop flag rather
/// than failing.
pub fn execute_line(session: &mut DrawSession, line: &str) -> Result<(), String> {
    let tokens = parse_command(line);
    if tokens.is_empty() {
        return Ok(());
    }
    // Expand $variables before dispatch so `set x 3; box b $x 2 2` behaves
    // like `box b 3 2 2`. Unknown names survive unchanged (`expand_vars` is
    // deliberately non-destructive), and `for`/`if` bodies re-expand per
    // iteration once their loop variable is bound.
    //
    // `proc`/`def` definitions are the one exception: their body is stored
    // verbatim so `$param` references are expanded at *call* time, after the
    // parameters are bound. Pre-expanding here would break a body whose
    // parameter name collides with a global variable.
    let is_proc_def = matches!(tokens[0].as_str(), "proc" | "def");
    let tokens = if is_proc_def {
        tokens
    } else {
        expand_vars(session, &tokens)
    };
    let command = DrawCommand::Tokens(tokens);
    let result = dispatch(session, &command);
    // `PROC_RETURN` is a control-flow signal, not a real error: it stops a
    // procedure body without poisoning `last_error`.
    if let Err(ref e) = result {
        if *e != PROC_RETURN {
            session.last_error = Some(e.clone());
        }
    }
    result
}

/// Route a parsed command to its handler.
///
/// User-defined procedures shadow built-ins: a name registered by `proc`/`def`
/// is invoked directly with the command's arguments. Everything else is looked
/// up in the centralized [`COMMANDS`] table; unknown names fall through to the
/// error arm.
pub(super) fn dispatch(session: &mut DrawSession, command: &DrawCommand) -> Result<(), String> {
    let name = command.name();
    if session.procs.contains_key(name) {
        return call_proc(session, name, command.args());
    }
    match COMMANDS.iter().find(|e| e.name == name) {
        Some(entry) => (entry.handler)(session, command.args()),
        None => Err(format!("draw: unknown command '{name}'")),
    }
}

// ---------------------------------------------------------------------------
// Batch drivers
// ---------------------------------------------------------------------------

/// Execute a multi-line script, one command per line.
///
/// Blank lines and lines whose first non-whitespace character is `#` are
/// comments and are skipped. Execution stops at the first failing command
/// (returning its error) or when `exit`/`quit` sets the stop flag. On success
/// the accumulated log is returned.
pub fn run_script(session: &mut DrawSession, script: &str) -> Result<Vec<String>, String> {
    for line in script.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Err(e) = execute_line(session, line) {
            return Err(e);
        }
        if session.stop {
            break;
        }
    }
    Ok(session.log.clone())
}

/// Read a script file and [`run_script`] it. Reading errors surface as command
/// errors so the caller only handles `Result<_, String>`.
pub fn run_script_file(session: &mut DrawSession, path: &str) -> Result<Vec<String>, String> {
    let script = std::fs::read_to_string(path)
        .map_err(|e| format!("draw: cannot read script '{path}': {e}"))?;
    run_script(session, &script)
}

// ---------------------------------------------------------------------------
// Interactive drivers
// ---------------------------------------------------------------------------

/// Interactive REPL over injected streams — the testable core of [`draw_repl`].
///
/// Reads one line at a time from `input`, executes it against a fresh session,
/// and writes the new log entries (or the error) to `output` prefixed with
/// `#`. Stops on `exit`/`quit` or end of input. Errors do not terminate the
/// loop; they are reported and the next line is read.
pub fn draw_repl_with(input: &mut dyn BufRead, output: &mut dyn Write) -> Result<(), String> {
    let mut session = DrawSession::default();
    let mut line = String::new();
    loop {
        line.clear();
        let n = input
            .read_line(&mut line)
            .map_err(|e| format!("draw: read: {e}"))?;
        if n == 0 {
            break; // EOF
        }
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let log_start = session.log.len();
        match execute_line(&mut session, line) {
            Ok(()) => {
                for entry in &session.log[log_start..] {
                    writeln!(output, "# {entry}").map_err(|e| format!("draw: write: {e}"))?;
                }
            }
            Err(e) => {
                writeln!(output, "# draw: error: {e}").map_err(|e| format!("draw: write: {e}"))?;
            }
        }
        if session.stop {
            break;
        }
    }
    Ok(())
}

/// Interactive REPL reading stdin and writing stdout.
///
/// This is the top-level entry point for a live shell. It is intentionally
/// separated from [`draw_repl_with`] so the loop can be driven from any
/// `BufRead`/`Write` pair (tests inject a `Cursor` and a `Vec<u8>`), and so
/// tests never block on the real stdin.
pub fn draw_repl() -> Result<(), String> {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    draw_repl_with(&mut stdin.lock(), &mut stdout.lock())
}

// ---------------------------------------------------------------------------
// Queries
// ---------------------------------------------------------------------------

/// Estimated enclosed volume of a closed shape, in cubic units.
///
/// This wraps [`shape_volume`], which integrates the divergence theorem over
/// the tessellated surface. For curved shapes (spheres, cylinders, toruses)
/// the estimate is reliable and converges to the exact volume as the
/// deflection shrinks — `info` uses it for that reason. For planar boxes the
/// tessellation is an exact triangulation of the boundary, so the estimate is
/// exact in practice too, but it is still a mesh estimate: prefer a primitive's
/// analytic `volume()` method when one is available, and treat the number from
/// an open or inconsistently oriented shell as meaningless.
pub fn volume(shape: &TopoShape) -> f64 {
    shape_volume(shape, VOLUME_DEFLECTION)
}
