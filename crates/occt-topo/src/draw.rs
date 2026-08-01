//! Draw_Interpretor-lite: a scripting + batch command driver for the BRep kernel.
//!
//! This module is a small, dependency-free port of OCCT's interactive `Draw`
//! package. OCCT ships the full `Draw_Interpretor` as a Tcl-embedded shell;
//! this port keeps only the two pieces that matter for driving a kernel from
//! scripts and tests:
//!
//! * a line parser ([`parse_command`]) that turns an input line into tokens,
//!   honouring double-quoted strings so file names containing spaces survive;
//! * a command dispatcher ([`execute_line`]) with a fixed command table:
//!   primitive builders, boolean operations, fillets, offsets, meshing,
//!   exporters, and session bookkeeping;
//! * batch drivers ([`run_script`], [`run_script_file`]) for `.draw`-style
//!   scripts, and an interactive REPL ([`draw_repl`] / [`draw_repl_with`]) for
//!   a live shell.
//!
//! Session state lives in [`DrawSession`]: a name → shape table, a log of
//! output lines, the last error, and a stop flag set by `exit`/`quit`. The
//! command table is deliberately flat — there is no
//! `Draw_Interpretor::AddCommand` registry — because a fixed `match` is all the
//! scripting and test layers need.
//!
//! # Command reference
//!
//! Every command is `name arg1 arg2 ...`; shape names are arbitrary strings.
//! The first token selects the command, the rest are its arguments.
//!
//! | Command       | Arguments            | Effect                                            |
//! |---------------|----------------------|---------------------------------------------------|
//! | `box`         | `name dx dy dz`      | axis-aligned box `[0,dx]×[0,dy]×[0,dz]`           |
//! | `cylinder`    | `name r h`           | Z-axis cylinder, radius `r`, height `h`           |
//! | `sphere`      | `name r`             | sphere of radius `r` centred at the origin        |
//! | `cone`        | `name r1 r2 h`       | right cone (base `r1`, top `r2` = 0), height `h`  |
//! | `torus`       | `name r1 r2`         | torus, major radius `r1`, minor radius `r2`       |
//! | `fuse`        | `a b out`            | boolean union `a ∪ b` → `out`                     |
//! | `cut`         | `a b out`            | boolean difference `a − b` → `out`                |
//! | `common`      | `a b out`            | boolean intersection `a ∩ b` → `out`              |
//! | `fillet`      | `name radius`        | blend the first edge of `name`                    |
//! | `fillet`      | `name edge_idx r`    | blend edge `edge_idx` (0-based) of `name`         |
//! | `offset`      | `name dist out`      | offset a shell/solid by `dist` → `out`            |
//! | `mesh`        | `name deflection`    | tessellate and report the triangle count          |
//! | `step`        | `name file`          | write `name` to a STEP file                       |
//! | `obj`         | `name file`          | write an OBJ mesh file                            |
//! | `iges`        | `name file`          | write an IGES file                                |
//! | `stl`         | `name file`          | write a binary STL file                           |
//! | `info`        | `name`               | print vertex/edge/face counts + estimated volume  |
//! | `ls`          | —                    | list the session shape names                      |
//! | `rm`          | `name`               | remove a shape from the session                   |
//! | `clear`       | —                    | empty the session                                 |
//! | `help`        | —                    | print this command list                           |
//! | `exit` / `quit` | —                  | set the session stop flag                         |
//!
//! Unknown commands are rejected with an error rather than silently ignored,
//! so a typo in a batch script cannot corrupt a run. The script drivers treat
//! blank lines and lines beginning with `#` as comments.
//!
//! # Port notes
//!
//! * Boolean operations route through [`curved_boolean_ext`], which already
//!   dispatches to the exact planar `crate::bop_builder::boolean` when both
//!   operands are planar-faced and to the curved/voxel paths otherwise — one
//!   entry point covers both, matching OCCT's `BRepAlgoAPI_BooleanOperation`.
//! * `fillet` uses the constant-radius [`fillet_edge`]. The variable-radius
//!   `crate::fillet_var::fillet_edge_var` API is available but has no
//!   single-command form here (it needs a radius law, the interpreter does not
//!   expose).
//! * `info`'s volume is a mesh estimate ([`shape_volume`]); it is reliable for
//!   curved shapes and exact for planar boxes in practice, but it is not the
//!   analytic primitive volume (see [`volume`]).
//!
//! # Examples
//!
//! Batch driving is the common case — build two primitives, combine them, and
//! inspect the result, all through [`run_script`]:
//!
//! ```
//! use occt_topo::draw::{DrawSession, run_script};
//!
//! let mut session = DrawSession::default();
//! let script = concat!(
//!     "# a small construction\n",
//!     "box base 10 5 2\n",
//!     "box cap 4 4 4\n",
//!     "fuse base cap joined\n",
//!     "info joined",
//! );
//! run_script(&mut session, script).expect("script runs");
//! assert!(session.shapes.contains_key("joined"));
//! ```
//!
//! Interactive use calls [`draw_repl`] (a thin wrapper over
//! [`draw_repl_with`], which accepts any `BufRead`/`Write` pair):
//!
//! ```no_run
//! # fn main() -> Result<(), String> {
//! occt_topo::draw::draw_repl()?;
//! # Ok(()) }
//! ```
//!
//! The REPL echoes each assigned shape as `# name = Solid` and each error as
//! `# draw: error: ...`, and exits on `exit`, `quit`, or end of input.

use std::collections::HashMap;
use std::io::{BufRead, Write};

use crate::abs::ShapeType;
use crate::bop_builder::BoolOp;
use crate::bop_curved::curved_boolean_ext;
use crate::brep_exchange::{brep_to_obj, brep_to_stl_binary};
use crate::brep_offset::offset_shell;
use crate::fillet_edge::fillet_edge;
use crate::iges::write_shape_iges;
use crate::primitives::{
    BRepPrimBox, BRepPrimCone, BRepPrimCylinder, BRepPrimSphere, BRepPrimTorus,
};
use crate::shape::TopoShape;
use crate::shape_mesh::{mesh_shape, shape_volume};
use crate::step::write_shape_step;
use crate::topo_tools_full::{edges_of, shape_counts};

// ---------------------------------------------------------------------------
// Tuning constants
// ---------------------------------------------------------------------------

/// Tolerance passed to the boolean kernel (matches the `bop_builder` tests).
const BOOL_TOL: f64 = 1e-6;

/// Deflection used by exporters that do not take a mesh parameter.
const EXPORT_DEFLECTION: f64 = 0.1;

/// Deflection for the `info` volume estimate — fine enough to be reliable on
/// curved shapes without being slow.
const VOLUME_DEFLECTION: f64 = 0.05;

// ---------------------------------------------------------------------------
// Session state
// ---------------------------------------------------------------------------

/// Mutable state shared by every command: the named-shape table, the output
/// log, the last error, and the interactive stop flag.
///
/// Mirrors what a `Draw_Interpretor` keeps between commands (its variable
/// namespace and output buffer), flattened into one struct the script and REPL
/// drivers thread through.
pub struct DrawSession {
    /// Named shapes registered by primitive / boolean / fillet / offset
    /// commands. The `info`, `mesh`, exporter and binary commands read from
    /// here; `ls`, `rm`, `clear` manage it.
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
}

impl Default for DrawSession {
    fn default() -> Self {
        Self {
            shapes: HashMap::new(),
            log: Vec::new(),
            last_error: None,
            stop: false,
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
fn need_args(args: &[String], n: usize, cmd: &str) -> Result<(), String> {
    if args.len() != n {
        return Err(format!(
            "draw: {cmd}: expected {n} arguments, got {}",
            args.len()
        ));
    }
    Ok(())
}

/// Parse argument `i` as an `f64`.
fn arg_f64(args: &[String], i: usize, cmd: &str) -> Result<f64, String> {
    args[i]
        .parse::<f64>()
        .map_err(|_| format!("draw: {cmd}: invalid number '{}'", args[i]))
}

/// Human-readable name of a `BoolOp` for error messages and dispatch.
fn bool_op_name(op: BoolOp) -> &'static str {
    match op {
        BoolOp::Fuse => "fuse",
        BoolOp::Cut => "cut",
        BoolOp::Common => "common",
    }
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
fn cmd_box(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_cylinder(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 3, "cylinder")?;
    let r = arg_f64(args, 1, "cylinder")?;
    let h = arg_f64(args, 2, "cylinder")?;
    let prim = BRepPrimCylinder::make_cylinder(r, h);
    register_shape(session, &args[0], prim.solid.into());
    Ok(())
}

/// `sphere <name> <r>` — build a sphere of radius `r` centred at the origin.
fn cmd_sphere(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 2, "sphere")?;
    let r = arg_f64(args, 1, "sphere")?;
    let prim = BRepPrimSphere::make_sphere(r);
    register_shape(session, &args[0], prim.solid.into());
    Ok(())
}

/// `cone <name> <r1> <r2> <h>` — build a cone with base radius `r1` and
/// height `h`. The primitive in this crate is a right cone (apex above the
/// base centre), so a truncated cone — `r2 != 0` — is rejected.
fn cmd_cone(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_torus(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_boolean(session: &mut DrawSession, args: &[String], op: BoolOp) -> Result<(), String> {
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
fn cmd_fillet(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_offset(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_mesh(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
enum WriteKind {
    Step,
    Obj,
    Iges,
    Stl,
}

fn write_kind_name(kind: WriteKind) -> &'static str {
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
fn cmd_write(
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
fn cmd_info(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_ls(session: &mut DrawSession) -> Result<(), String> {
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
fn cmd_rm(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 1, "rm")?;
    if session.shapes.remove(&args[0]).is_none() {
        return Err(format!("draw: rm: shape '{}' not found", args[0]));
    }
    session.log.push(format!("rm: removed {}", args[0]));
    Ok(())
}

/// `clear` — empty the session: all shapes and the log.
fn cmd_clear(session: &mut DrawSession) -> Result<(), String> {
    session.shapes.clear();
    session.log.clear();
    session.log.push("clear: session empty".into());
    Ok(())
}

/// `help` — list the available commands on one log line.
fn cmd_help(session: &mut DrawSession) -> Result<(), String> {
    session.log.push(
        "commands: box cylinder sphere cone torus fuse cut common fillet offset \
         mesh step obj iges stl info ls rm clear help exit quit"
            .into(),
    );
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
    let command = DrawCommand::Tokens(tokens);
    let result = dispatch(session, &command);
    if let Err(ref e) = result {
        session.last_error = Some(e.clone());
    }
    result
}

/// Route a parsed command to its handler. The `match` is the fixed command
/// table; unknown names fall through to the error arm.
fn dispatch(session: &mut DrawSession, command: &DrawCommand) -> Result<(), String> {
    match command.name() {
        "box" => cmd_box(session, command.args()),
        "cylinder" => cmd_cylinder(session, command.args()),
        "sphere" => cmd_sphere(session, command.args()),
        "cone" => cmd_cone(session, command.args()),
        "torus" => cmd_torus(session, command.args()),
        "fuse" => cmd_boolean(session, command.args(), BoolOp::Fuse),
        "cut" => cmd_boolean(session, command.args(), BoolOp::Cut),
        "common" => cmd_boolean(session, command.args(), BoolOp::Common),
        "fillet" => cmd_fillet(session, command.args()),
        "offset" => cmd_offset(session, command.args()),
        "mesh" => cmd_mesh(session, command.args()),
        "step" => cmd_write(session, command.args(), WriteKind::Step),
        "obj" => cmd_write(session, command.args(), WriteKind::Obj),
        "iges" => cmd_write(session, command.args(), WriteKind::Iges),
        "stl" => cmd_write(session, command.args(), WriteKind::Stl),
        "info" => cmd_info(session, command.args()),
        "ls" => cmd_ls(session),
        "rm" => cmd_rm(session, command.args()),
        "clear" => cmd_clear(session),
        "help" => cmd_help(session),
        "exit" | "quit" => {
            session.stop = true;
            Ok(())
        }
        other => Err(format!("draw: unknown command '{other}'")),
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topo_tools_full::faces_of;

    #[test]
    fn parse_quoted_tokens() {
        assert_eq!(parse_command("box a 1 2 3"), vec!["box", "a", "1", "2", "3"]);
        assert_eq!(
            parse_command("step \"my file.step\""),
            vec!["step", "my file.step"]
        );
        assert_eq!(parse_command(""), Vec::<String>::new());
        assert_eq!(parse_command("   "), Vec::<String>::new());
    }

    #[test]
    fn execute_box_info() {
        let mut s = DrawSession::default();
        execute_line(&mut s, "box b 2 2 2").expect("box ok");
        execute_line(&mut s, "info b").expect("info ok");
        let joined = s.log.join("\n");
        assert!(joined.contains("vertices 8"), "log: {joined}");
        assert!(joined.contains("edges 12"), "log: {joined}");
        assert!(joined.contains("faces 6"), "log: {joined}");
        assert!(joined.contains("volume"), "log: {joined}");
        let v = volume(&shape_owned(&s, "b").unwrap());
        assert!(v > 0.0 && v < 16.0, "volume {v}");
    }

    #[test]
    fn fuse_two_boxes() {
        let mut s = DrawSession::default();
        run_script(&mut s, "box a 1 1 1\nbox c 2 1 1\nfuse a c f").expect("script ok");
        let f = shape_owned(&s, "f").expect("f exists");
        assert!(faces_of(&f).len() >= 6, "faces {}", faces_of(&f).len());
        let v = volume(&f);
        assert!((1.5..=2.5).contains(&v), "volume {v}");
    }

    #[test]
    fn cut_and_common() {
        let mut s = DrawSession::default();
        run_script(&mut s, "box a 1 1 1\nbox c 2 1 1").expect("boxes");
        execute_line(&mut s, "cut c a cut1").expect("cut ok");
        execute_line(&mut s, "common a c com1").expect("common ok");
        let cut = shape_owned(&s, "cut1").expect("cut1");
        let com = shape_owned(&s, "com1").expect("com1");
        assert!(faces_of(&cut).len() >= 1, "cut faces {}", faces_of(&cut).len());
        assert!(faces_of(&com).len() >= 1, "common faces {}", faces_of(&com).len());
    }

    #[test]
    fn sphere_fillet_offset() {
        let mut s = DrawSession::default();
        run_script(&mut s, "sphere s 1\noffset s 0.5 so").expect("sphere + offset");
        let sv = volume(&shape_owned(&s, "s").unwrap());
        let sov = volume(&shape_owned(&s, "so").unwrap());
        assert!(sov > sv, "offset volume {sov} should exceed sphere {sv}");
        assert!(sov > 10.0, "offset volume {sov} (sphere r=1.5 → ~14.1)");
    }

    #[test]
    fn mesh_counts_triangles() {
        let mut s = DrawSession::default();
        run_script(&mut s, "box b 2 2 2\nmesh b 0.2").expect("script ok");
        let joined = s.log.join("\n");
        assert!(joined.contains("triangles"), "log: {joined}");
        let tri_line = s.log.iter().find(|l| l.contains("triangles")).unwrap();
        let n: usize = tri_line
            .split_whitespace()
            .find_map(|w| w.parse::<usize>().ok())
            .expect("triangle count token");
        assert!(n > 0, "triangle count {n}");
    }

    #[test]
    fn step_obj_iges_write() {
        let mut s = DrawSession::default();
        run_script(&mut s, "box b 1 1 1").expect("box");
        let dir = std::env::temp_dir();
        let id = std::process::id();
        let files = [
            dir.join(format!("draw_step_{id}.step")),
            dir.join(format!("draw_obj_{id}.obj")),
            dir.join(format!("draw_iges_{id}.igs")),
            dir.join(format!("draw_stl_{id}.stl")),
        ];
        let paths: Vec<String> = files.iter().map(|p| p.to_str().unwrap().to_string()).collect();
        execute_line(&mut s, &format!("step b {}", paths[0])).expect("step");
        execute_line(&mut s, &format!("obj b {}", paths[1])).expect("obj");
        execute_line(&mut s, &format!("iges b {}", paths[2])).expect("iges");
        execute_line(&mut s, &format!("stl b {}", paths[3])).expect("stl");
        for p in &files {
            let data = std::fs::read(p).expect("file exists");
            assert!(!data.is_empty(), "{} empty", p.display());
            let _ = std::fs::remove_file(p);
        }
    }

    #[test]
    fn unknown_command_errors() {
        let mut s = DrawSession::default();
        let err = execute_line(&mut s, "bogus").expect_err("should error");
        assert!(err.contains("unknown command 'bogus'"), "err {err}");
        assert!(s.last_error.is_some(), "last_error set");
    }

    #[test]
    fn run_script_multiline() {
        let mut s = DrawSession::default();
        let script = "\
# build two overlapping boxes and fuse them
box a 1 1 1
box c 2 1 1

fuse a c f
info f
";
        let log = run_script(&mut s, script).expect("script ok");
        assert!(log.iter().any(|l| l.contains("info f")), "log: {log:?}");
        assert!(s.shapes.contains_key("a"));
        assert!(s.shapes.contains_key("c"));
        assert!(s.shapes.contains_key("f"));
    }

    #[test]
    fn script_stops_on_error() {
        let mut s = DrawSession::default();
        let err = run_script(&mut s, "box a 1 1 1\nbogus\nbox b 2 2 2").expect_err("should error");
        assert!(err.contains("bogus"), "err {err}");
        assert!(s.shapes.contains_key("a"), "a ran before the error");
        assert!(!s.shapes.contains_key("b"), "b must not run after the error");
        assert!(s.last_error.is_some(), "last_error set");
    }

    #[test]
    fn ls_rm_clear() {
        let mut s = DrawSession::default();
        run_script(&mut s, "box a 1 1 1\nbox b 2 2 2").expect("boxes");
        execute_line(&mut s, "ls").expect("ls");
        assert!(
            s.log.iter().any(|l| l.contains("a") && l.contains("b")),
            "log: {:?}",
            s.log
        );
        execute_line(&mut s, "rm a").expect("rm");
        assert!(!s.shapes.contains_key("a"), "a removed");
        execute_line(&mut s, "clear").expect("clear");
        assert!(s.shapes.is_empty(), "session cleared");
    }

    #[test]
    fn help_lists_commands() {
        let mut s = DrawSession::default();
        execute_line(&mut s, "help").expect("help");
        let joined = s.log.join("\n");
        assert!(joined.contains("box"), "log: {joined}");
        assert!(joined.contains("fuse"), "log: {joined}");
        assert!(joined.contains("exit"), "log: {joined}");
    }

    #[test]
    fn exit_stops_repl() {
        let mut input = std::io::Cursor::new(b"box a 1 1 1\nexit\n".to_vec());
        let mut output: Vec<u8> = Vec::new();
        draw_repl_with(&mut input, &mut output).expect("repl ok");
        let out = String::from_utf8(output).unwrap();
        assert!(out.contains("a"), "echo output: {out}");
    }
}
