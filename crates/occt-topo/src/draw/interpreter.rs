use super::prelude::*;
use super::*;

// ---------------------------------------------------------------------------
// Tuning constants
// ---------------------------------------------------------------------------

/// Tolerance passed to the boolean kernel (matches the `bop_builder` tests).

pub(super) const BOOL_TOL: f64 = 1e-6;

/// Deflection used by exporters that do not take a mesh parameter.
pub(super) const EXPORT_DEFLECTION: f64 = 0.1;

/// Deflection for the `info` volume estimate — fine enough to be reliable on
/// curved shapes without being slow.
pub(super) const VOLUME_DEFLECTION: f64 = 0.05;

/// Maximum `proc` call depth. Recursion is supported, but a runaway recursive
/// procedure must fail with a clear error instead of exhausting the Rust stack.
pub(super) const PROC_MAX_DEPTH: usize = 64;

/// Sentinel error used by `return` to stop a procedure body. It is caught by
/// [`call_proc`] (which turns it into a normal `Ok(())`) and never recorded as
/// a session error; the control characters make a collision with a real error
/// message impossible.
pub(super) const PROC_RETURN: &str = "\u{1}draw:proc:return\u{1}";

// ---------------------------------------------------------------------------
// Session state
// ---------------------------------------------------------------------------

/// A user-defined procedure: a parameter list and a body token list.
///
/// `proc <name> { p1 p2 ... } { body }` stores the body *verbatim* (variable
/// references are expanded at call time, once parameters are bound), so a body
/// may reference its own parameters even when a global variable of the same
/// name exists. The body is a command sequence: `;` separates commands, and
/// each is executed through [`execute_line`] with the parameters bound as
/// session variables.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProcDef {
    /// The procedure name (also the key in [`DrawSession::procs`]).
    pub name: String,
    /// Parameter names, bound positionally from call arguments.
    pub params: Vec<String>,
    /// Body tokens: a `;`-separated command sequence (outer braces stripped).
    pub body: Vec<String>,
}

/// Mutable state shared by every command: the named-shape table, the output
/// log, the last error, the interactive stop flag, a TCL-style variable
/// namespace and a session camera.
///
/// Mirrors what a `Draw_Interpretor` keeps between commands (its variable
/// namespace and output buffer), flattened into one struct the script and REPL
/// drivers thread through.
pub struct DrawSession {
    /// Named shapes registered by primitive / boolean / fillet / offset /
    /// transform commands. The `info`, `bbox`, `verts`, `mesh`, exporter and
    /// binary commands read from here; `ls`, `rm`, `clear` manage it.
    pub shapes: HashMap<String, TopoShape>,
    /// Output lines accumulated while executing commands. The REPL prints new
    /// entries after each line; `run_script` returns the whole log.
    pub log: Vec<String>,
    /// The error message of the most recent failed command, if any. Batch
    /// drivers abort on the first error; the REPL reports it and continues.
    pub last_error: Option<String>,
    /// Set by `exit`/`quit`. Script and REPL loops check it after every line
    /// and stop — the interpreter equivalent of leaving `Draw_Interpretor`.
    pub stop: bool,
    /// TCL-style variable namespace: `set <name> <value>` writes here, and
    /// `$name` references in any command line expand from here before
    /// dispatch. `expr` stores its result in the `result` entry; `for` writes
    /// its loop variable here for each iteration.
    pub vars: HashMap<String, String>,
    /// The session view camera used by `view`, and driven by `view_camera`,
    /// `orbit`, `zoom` and `pan`.
    pub camera: Camera,
    /// User-defined procedures (`proc`/`def`), by name. A call to a name in
    /// here dispatches to the procedure rather than the built-in table.
    pub procs: HashMap<String, ProcDef>,
    /// The active procedure call stack (innermost last). Non-empty while a
    /// `proc` body is executing; used by `return` to know it must stop the
    /// body, and by error messages to render the call chain.
    pub proc_stack: Vec<String>,
}

impl Default for DrawSession {
    fn default() -> Self {
        Self {
            shapes: HashMap::new(),
            log: Vec::new(),
            last_error: None,
            stop: false,
            vars: HashMap::new(),
            camera: Camera::default(),
            procs: HashMap::new(),
            proc_stack: Vec::new(),
        }
    }
}

// ---------------------------------------------------------------------------
// Parsed command representation
// ---------------------------------------------------------------------------

/// A parsed command line: the raw token stream after whitespace splitting and
/// quote handling. The first token names the command; the rest are its
/// arguments.
///
/// This port keeps the token list in a single `Tokens` variant rather than a
/// per-command enum because the dispatch table is fixed — [`execute_line`]
/// reads the command name directly and matches it, which is the OCCT
/// `Draw_Interpretor` equivalent of looking up a registered command handler.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DrawCommand {
    /// The token list, command name first.
    Tokens(Vec<String>),
}

impl DrawCommand {
    /// The command name (the first token), or `""` for an empty command.
    pub fn name(&self) -> &str {
        match self {
            DrawCommand::Tokens(t) => t.first().map(String::as_str).unwrap_or(""),
        }
    }

    /// The argument tokens — everything after the command name.
    pub fn args(&self) -> &[String] {
        match self {
            DrawCommand::Tokens(t) => &t[1..],
        }
    }
}

// ---------------------------------------------------------------------------
// Line parser
// ---------------------------------------------------------------------------

/// Split an input line into command tokens.
///
/// Tokens are separated by whitespace. A double-quoted section is kept as a
/// single token (the quotes are stripped) so that paths and shape names
/// containing spaces — e.g. `step b "my file.step"` — survive the split. This
/// is a deliberately small subset of what OCCT's `Draw` does with Tcl quoting:
/// no escapes, no nested quotes, no variable substitution.
pub fn parse_command(line: &str) -> Vec<String> {
    let mut tokens: Vec<String> = Vec::new();
    let mut cur = String::new();
    let mut quoted = false;
    for c in line.chars() {
        match c {
            '"' => quoted = !quoted,
            c if c.is_whitespace() && !quoted => {
                if !cur.is_empty() {
                    tokens.push(std::mem::take(&mut cur));
                }
            }
            c => cur.push(c),
        }
    }
    if !cur.is_empty() {
        tokens.push(cur);
    }
    tokens
}

// ---------------------------------------------------------------------------
// Argument helpers
// ---------------------------------------------------------------------------

/// Reject a command whose argument count does not match `n` exactly.
pub(super) fn need_args(args: &[String], n: usize, cmd: &str) -> Result<(), String> {
    if args.len() != n {
        return Err(format!(
            "draw: {cmd}: expected {n} arguments, got {}",
            args.len()
        ));
    }
    Ok(())
}

/// Parse argument `i` as an `f64`.
pub(super) fn arg_f64(args: &[String], i: usize, cmd: &str) -> Result<f64, String> {
    args[i]
        .parse::<f64>()
        .map_err(|_| format!("draw: {cmd}: invalid number '{}'", args[i]))
}

/// Human-readable name of a `BoolOp` for error messages and dispatch.
pub(super) fn bool_op_name(op: BoolOp) -> &'static str {
    match op {
        BoolOp::Fuse => "fuse",
        BoolOp::Cut => "cut",
        BoolOp::Common => "common",
    }
}

// ---------------------------------------------------------------------------
// Variables, numbers and expressions
// ---------------------------------------------------------------------------

/// Render a real number the way the Tcl layer prints one: integral values
/// drop the `.0` (`3`, not `3.0`), everything else keeps the compact `f64`
/// rendering. This keeps `expr 3 4 *` storing `result = 12` and round-trips
/// through [`eval_number`].
pub(super) fn fmt_num(v: f64) -> String {
    if v.fract() == 0.0 && v.abs() < 9.0e15 {
        format!("{}", v as i64)
    } else {
        format!("{v}")
    }
}

/// Substitute every `$name` reference inside `tok` from the session's
/// variable namespace.
///
/// A variable name is the longest run of alphanumeric/underscore characters
/// after the `$`. Whole tokens (`$x`), embedded references (`b$i`) and
/// repeated references (`a$i$j`) all work. Unknown names are left exactly as
/// written — a literal `$` in a script (say a file name) is preserved instead
/// of erroring, which is friendlier than Tcl and keeps scripts non-brittle.
pub(super) fn expand_token(session: &DrawSession, tok: &str) -> String {
    let mut out = String::new();
    let mut rest = tok;
    while let Some(idx) = rest.find('$') {
        out.push_str(&rest[..idx]);
        let after = &rest[idx + 1..];
        let mut name = String::new();
        let mut name_bytes = 0;
        for c in after.chars() {
            if c.is_alphanumeric() || c == '_' {
                name.push(c);
                name_bytes += c.len_utf8();
            } else {
                break;
            }
        }
        // Array element: `$arr(0)` names the variable `arr(0)`. An index group
        // is any balanced `(...)` immediately after the base name; it is folded
        // into the lookup key so `set arr(0) 5` / `$arr(0)` round-trip.
        if name_bytes > 0 && after[name_bytes..].starts_with('(') {
            if let Some(close) = after[name_bytes..].find(')') {
                let end = name_bytes + close + 1;
                name.push_str(&after[name_bytes..end]);
                name_bytes = end;
            }
        }
        if name_bytes == 0 {
            out.push('$');
            rest = after;
            continue;
        }
        match session.vars.get(&name) {
            Some(v) => out.push_str(v),
            None => {
                out.push('$');
                out.push_str(&name);
            }
        }
        rest = &after[name_bytes..];
    }
    out.push_str(rest);
    out
}

/// Expand `$name` references in a token list — the substitution step applied
/// to every command line before dispatch. See [`expand_token`].
pub fn expand_vars(session: &DrawSession, tokens: &[String]) -> Vec<String> {
    tokens.iter().map(|t| expand_token(session, t)).collect()
}

/// Parse a number, or a `$variable` holding one, as an `f64`.
///
/// This is the numeric-argument resolver for `expr`, `if`, `for`, `translate`
/// and friends. A `$result` reference (an `expr` output) round-trips through
/// [`fmt_num`], so `eval_number(session, "$result")` returns the computed
/// value. Non-numeric text is an error.
pub fn eval_number(session: &DrawSession, s: &str) -> Result<f64, String> {
    let expanded = expand_token(session, s);
    expanded
        .trim()
        .parse::<f64>()
        .map_err(|_| format!("draw: invalid number '{s}'"))
}

/// Store `value` in the session variable `name` (creating or overwriting it).
/// Values are kept as strings — a number stored here is parsed lazily by
/// [`eval_number`] when a command consumes it.
pub fn set_var(session: &mut DrawSession, name: &str, value: &str) {
    session.vars.insert(name.to_string(), value.to_string());
}

/// Split the trailing block of a `for`/`if` command off its argument list.
///
/// The header is `args[0..header]`; the body is everything after it. When a
/// `{` token appears exactly at `header` the body runs to the end, dropping a
/// trailing `}`. Without braces the body is simply the tail of `args`. An
/// empty body is an error so a bare `for i 1 3 { }` cannot silently do
/// nothing.
pub(super) fn command_block<'a>(args: &'a [String], keyword: &str, header: usize) -> Result<&'a [String], String> {
    let body = if args.get(header).map(|s| s == "{").unwrap_or(false) {
        let mut end = args.len();
        if args.last().map(|s| s == "}").unwrap_or(false) {
            end -= 1;
        }
        &args[header + 1..end]
    } else {
        &args[header..]
    };
    if body.is_empty() {
        return Err(format!("draw: {keyword}: empty body"));
    }
    Ok(body)
}

// ---------------------------------------------------------------------------
// Session helpers
// ---------------------------------------------------------------------------

/// Look up a named shape, cloning it out of the session so callers can then
/// mutate the session (registering results) without a borrow conflict.
///
/// This is the interpreter equivalent of resolving a variable name; unknown
/// names are an error, matching OCCT's `Draw` behaviour of failing a command
/// that references an undefined shape.
pub fn shape_owned(session: &DrawSession, name: &str) -> Result<TopoShape, String> {
    session
        .shapes
        .get(name)
        .cloned()
        .ok_or_else(|| format!("draw: shape '{name}' not found"))
}

/// Insert `shape` into the session under `name`, echoing a `name = type` line
/// into the log the way OCCT's `Draw` echoes the assigned variable.
pub fn register_shape(session: &mut DrawSession, name: &str, shape: TopoShape) {
    session
        .log
        .push(format!("{name} = {}", shape.shape_type().to_str()));
    session.shapes.insert(name.to_string(), shape);
}

// ---------------------------------------------------------------------------
// Command handlers
// ---------------------------------------------------------------------------

/// `box <name> <dx> <dy> <dz>` — build an axis-aligned box `[0,dx]×[0,dy]×[0,dz]`
/// and register it under `name`.
pub(super) fn cmd_box(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 4, "box")?;
    let dx = arg_f64(args, 1, "box")?;
    let dy = arg_f64(args, 2, "box")?;
    let dz = arg_f64(args, 3, "box")?;
    let prim = BRepPrimBox::make_box(dx, dy, dz);
    register_shape(session, &args[0], prim.solid.into());
    Ok(())
}

/// `cylinder <name> <r> <h>` — build a Z-axis cylinder of radius `r` and
/// height `h`, registering it under `name`.
pub(super) fn cmd_cylinder(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 3, "cylinder")?;
    let r = arg_f64(args, 1, "cylinder")?;
    let h = arg_f64(args, 2, "cylinder")?;
    let prim = BRepPrimCylinder::make_cylinder(r, h);
    register_shape(session, &args[0], prim.solid.into());
    Ok(())
}

/// `sphere <name> <r>` — build a sphere of radius `r` centred at the origin.
pub(super) fn cmd_sphere(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 2, "sphere")?;
    let r = arg_f64(args, 1, "sphere")?;
    let prim = BRepPrimSphere::make_sphere(r);
    register_shape(session, &args[0], prim.solid.into());
    Ok(())
}

/// `cone <name> <r1> <r2> <h>` — build a cone with base radius `r1` and
/// height `h`. The primitive in this crate is a right cone (apex above the
/// base centre), so a truncated cone — `r2 != 0` — is rejected.
pub(super) fn cmd_cone(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 4, "cone")?;
    let r1 = arg_f64(args, 1, "cone")?;
    let r2 = arg_f64(args, 2, "cone")?;
    let h = arg_f64(args, 3, "cone")?;
    if r2.abs() > 1e-9 {
        return Err(
            "draw: cone: truncated cones (top radius != 0) are not supported; use r2 = 0".into(),
        );
    }
    let prim = BRepPrimCone::make_cone(r1, h);
    register_shape(session, &args[0], prim.solid.into());
    Ok(())
}

/// `torus <name> <r1> <r2>` — build a torus with major radius `r1` and minor
/// (tube) radius `r2`.
pub(super) fn cmd_torus(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 3, "torus")?;
    let r1 = arg_f64(args, 1, "torus")?;
    let r2 = arg_f64(args, 2, "torus")?;
    let prim = BRepPrimTorus::make_torus(r1, r2);
    register_shape(session, &args[0], prim.solid.into());
    Ok(())
}

/// `fuse a b out` / `cut a b out` / `common a b out` — boolean operation.
///
/// Both operands are looked up by name and combined with `op`, storing the
/// result under `out`. Routing goes through [`curved_boolean_ext`], which
/// dispatches to the exact planar [`boolean`] when both operands are
/// planar-faced and to the curved / voxel paths otherwise, so a single entry
/// point covers all inputs (OCCT's `BRepAlgoAPI_BooleanOperation`).
pub(super) fn cmd_boolean(session: &mut DrawSession, args: &[String], op: BoolOp) -> Result<(), String> {
    need_args(args, 3, bool_op_name(op))?;
    let a = shape_owned(session, &args[0])?;
    let b = shape_owned(session, &args[1])?;
    let result = curved_boolean_ext(&a, &b, op, BOOL_TOL)
        .map_err(|e| format!("draw: {}: {e}", bool_op_name(op)))?;
    register_shape(session, &args[2], result.shape);
    Ok(())
}

/// `fillet <name> <radius>` or `fillet <name> <edge_idx> <radius>`.
///
/// With two arguments, the first edge of `name` is filleted. With three, the
/// edge at the 0-based index into [`edges_of`] is filleted. The rebuilt solid
/// replaces `name` in the session, matching OCCT's in-place `fillet` command.
pub(super) fn cmd_fillet(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    let (name, radius) = match args.len() {
        2 => (args[0].clone(), arg_f64(args, 1, "fillet")?),
        3 => (args[0].clone(), arg_f64(args, 2, "fillet")?),
        n => return Err(format!("draw: fillet: expected 2 or 3 arguments, got {n}")),
    };
    let shape = shape_owned(session, &name)?;
    let edges = edges_of(&shape);
    let edge = match args.len() {
        2 => edges.first().cloned().ok_or("draw: fillet: shape has no edges")?,
        _ => {
            let idx = args[1]
                .parse::<usize>()
                .map_err(|_| format!("draw: fillet: invalid edge index '{}'", args[1]))?;
            edges
                .get(idx)
                .cloned()
                .ok_or_else(|| format!("draw: fillet: edge index {idx} out of range ({})", edges.len()))?
        }
    };
    let result = fillet_edge(&shape, &edge, radius).map_err(|e| format!("draw: fillet: {e}"))?;
    register_shape(session, &name, result);
    Ok(())
}

/// `offset <name> <dist> <out>` — offset a shell or solid by `dist`, storing
/// the result under `out`.
pub(super) fn cmd_offset(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 3, "offset")?;
    let shape = shape_owned(session, &args[0])?;
    let dist = arg_f64(args, 1, "offset")?;
    let result = offset_shell(&shape, dist).map_err(|e| format!("draw: offset: {e}"))?;
    register_shape(session, &args[2], result);
    Ok(())
}

/// `mesh <name> <deflection>` — tessellate `name` and report the triangle
/// count. The mesh itself is discarded; this is a quick way to confirm a shape
/// tessellates and to gauge its size.
pub(super) fn cmd_mesh(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 2, "mesh")?;
    let shape = shape_owned(session, &args[0])?;
    let deflection = arg_f64(args, 1, "mesh")?;
    let mesh = mesh_shape(&shape, deflection);
    session
        .log
        .push(format!("mesh {}: {} triangles", args[0], mesh.triangles.len()));
    Ok(())
}

/// Which on-disk format a `step/obj/iges/stl` command writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum WriteKind {
    Step,
    Obj,
    Iges,
    Stl,
}

pub(super) fn write_kind_name(kind: WriteKind) -> &'static str {
    match kind {
        WriteKind::Step => "step",
        WriteKind::Obj => "obj",
        WriteKind::Iges => "iges",
        WriteKind::Stl => "stl",
    }
}

/// `step <name> <file>` / `obj <name> <file>` / `iges <name> <file>` /
/// `stl <name> <file>` — serialize `name` to `file` in the given format.
///
/// STEP and IGES are written from the analytic B-Rep; OBJ and STL are mesh
/// exports tessellated at [`EXPORT_DEFLECTION`]. The file is written
/// atomically enough for scripting purposes, and the log records the path.
pub(super) fn cmd_write(
    session: &mut DrawSession,
    args: &[String],
    kind: WriteKind,
) -> Result<(), String> {
    need_args(args, 2, write_kind_name(kind))?;
    let shape = shape_owned(session, &args[0])?;
    let file = &args[1];
    let kind_name = write_kind_name(kind);
    match kind {
        WriteKind::Step => {
            let text = write_shape_step(&shape);
            std::fs::write(file, text).map_err(|e| format!("draw: {kind_name}: {e}"))?;
        }
        WriteKind::Obj => {
            let text = brep_to_obj(&shape, EXPORT_DEFLECTION);
            std::fs::write(file, text).map_err(|e| format!("draw: {kind_name}: {e}"))?;
        }
        WriteKind::Iges => {
            let text = write_shape_iges(&shape);
            std::fs::write(file, text).map_err(|e| format!("draw: {kind_name}: {e}"))?;
        }
        WriteKind::Stl => {
            let bytes = brep_to_stl_binary(&shape, EXPORT_DEFLECTION);
            std::fs::write(file, bytes).map_err(|e| format!("draw: {kind_name}: {e}"))?;
        }
    }
    session.log.push(format!("{kind_name}: wrote '{file}'"));
    Ok(())
}

/// `info <name>` — print the distinct vertex / edge / face counts and an
/// estimated enclosed volume for `name`.
pub(super) fn cmd_info(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 1, "info")?;
    let shape = shape_owned(session, &args[0])?;
    let counts = shape_counts(&shape);
    let nv = counts.get(&ShapeType::Vertex).copied().unwrap_or(0);
    let ne = counts.get(&ShapeType::Edge).copied().unwrap_or(0);
    let nf = counts.get(&ShapeType::Face).copied().unwrap_or(0);
    let vol = volume(&shape);
    session.log.push(format!(
        "info {}: vertices {nv}, edges {ne}, faces {nf}, volume {vol:.6}",
        args[0]
    ));
    Ok(())
}

/// `ls` — list the session shape names, sorted, on one log line.
pub(super) fn cmd_ls(session: &mut DrawSession, _args: &[String]) -> Result<(), String> {
    let mut names: Vec<&String> = session.shapes.keys().collect();
    names.sort();
    let joined = names
        .iter()
        .map(|s| s.as_str())
        .collect::<Vec<_>>()
        .join(" ");
    session.log.push(format!("ls: {joined}"));
    Ok(())
}

/// `rm <name>` — remove a named shape. Unknown names are an error.
pub(super) fn cmd_rm(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 1, "rm")?;
    if session.shapes.remove(&args[0]).is_none() {
        return Err(format!("draw: rm: shape '{}' not found", args[0]));
    }
    session.log.push(format!("rm: removed {}", args[0]));
    Ok(())
}

/// `clear` — empty the session: all shapes and the log.
pub(super) fn cmd_clear(session: &mut DrawSession, _args: &[String]) -> Result<(), String> {
    session.shapes.clear();
    session.log.clear();
    session.log.push("clear: session empty".into());
    Ok(())
}

/// `help` — list every registered command with its arity and description, one
/// per log line. User-defined procedures are listed too, marked `(proc)`.
pub(super) fn cmd_help(session: &mut DrawSession, _args: &[String]) -> Result<(), String> {
    let mut entries: Vec<&CommandEntry> = COMMANDS.iter().collect();
    entries.sort_by_key(|e| e.name);
    session.log.push("commands (arity — description):".into());
    for e in entries {
        session
            .log
            .push(format!("  {:<12} {:<10} {}", e.name, arity_text(e), e.summary));
    }
    let mut procs: Vec<&String> = session.procs.keys().collect();
    procs.sort();
    for name in procs {
        session.log.push(format!("  {:<12} (proc)", name));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Variables and expressions
// ---------------------------------------------------------------------------

/// `set <name> <value>` — store `value` (a string or a number) in the session
/// variable namespace. The value is stored verbatim after `$`-expansion, so
/// `set b $a` copies `a`. Later `$name` references expand anywhere in a
/// command line — whole tokens (`$x`) and embedded in longer tokens (`b$i`).
pub(super) fn cmd_set(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 2, "set")?;
    set_var(session, &args[0], &args[1]);
    session.log.push(format!("set {} = {}", args[0], args[1]));
    Ok(())
}

/// `expr <a> <b> <op>` — evaluate `a op b` on real numbers with `op` one of
/// `+ - * /` and store the result in the `result` variable. The operator is
/// the final token (postfix, matching the Draw REPL convention `expr 3 4 *`),
/// and the operands may be literals or `$variables`; the stored text
/// round-trips through [`eval_number`], so `box c $result 1 1` consumes the
/// value directly.
pub(super) fn cmd_expr(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 3, "expr")?;
    let a = eval_number(session, &args[0])?;
    let b = eval_number(session, &args[1])?;
    let result = match args[2].as_str() {
        "+" => a + b,
        "-" => a - b,
        "*" => a * b,
        "/" => {
            if b.abs() < 1e-15 {
                return Err("draw: expr: division by zero".into());
            }
            a / b
        }
        op => return Err(format!("draw: expr: unsupported operator '{op}'")),
    };
    let text = fmt_num(result);
    set_var(session, "result", &text);
    session.log.push(format!("result = {text}"));
    Ok(())
}

/// `echo <text...>` — append the (already `$`-expanded) arguments joined by
/// single spaces to the log. The scripting printf, useful interactively and
/// from `if`/`for` bodies.
pub(super) fn cmd_echo(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    session.log.push(args.join(" "));
    Ok(())
}

/// `vars` — list every session variable as `name = value`, sorted by name, on
/// one log line. The analogue of Tcl's `info vars` and handy before a `for`
/// loop to confirm what the loop variable will shadow.
pub(super) fn cmd_vars(session: &mut DrawSession, _args: &[String]) -> Result<(), String> {
    let mut names: Vec<&String> = session.vars.keys().collect();
    names.sort();
    let joined = names
        .iter()
        .map(|n| format!("{n} = {}", session.vars[n.as_str()]))
        .collect::<Vec<_>>()
        .join("  ");
    session.log.push(format!("vars: {joined}"));
    Ok(())
}
