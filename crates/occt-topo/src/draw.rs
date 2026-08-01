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
//! | `bbox`        | `name`               | print the axis-aligned bounding-box corners       |
//! | `verts`       | `name`               | print the distinct vertex coordinates             |
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
//! The scripting layer adds these commands on top of the shape table:
//!
//! | Command       | Arguments               | Effect                                 |
//! |---------------|-------------------------|----------------------------------------|
//! | `set`         | `name value`            | store a string/number variable         |
//! | `expr`        | `a b op`                | `a op b` (`+ - * /`) into `result`     |
//! | `echo`        | `text...`               | append expanded text to the log        |
//! | `vars`        | —                       | list the session variables             |
//! | `unset`       | `name`                  | remove a variable                      |
//! | `incr`        | `name [step]`           | add `step` (default 1) to a number var |
//! | `for`         | `var start end { body }`| run body per integer in `[start, end]` |
//! | `if`          | `a b op { body }`       | run body when the comparison holds     |
//! | `translate`   | `name dx dy dz [out]`   | translate the registered geometry      |
//! | `rotate`      | `name ax ay az deg`     | rotate about an axis through the origin|
//! | `scale`       | `name factor [out]`     | uniform scale about the origin         |
//! | `mirror`      | `name nx ny nz [out]`   | mirror across a plane through the origin|
//! | `copy`        | `name out`              | independent deep copy                  |
//! | `transform`   | `name a11..a33 tx ty tz`| arbitrary affine map                   |
//! | `view`        | `name [w h] [file]`     | soft-render `name` to a PPM file       |
//! | `view_camera` | `ex ey ez tx ty tz`     | set the look-at camera                 |
//! | `orbit`       | `dx_deg dy_deg`         | orbit the camera about its target      |
//! | `zoom`        | `factor`                | zoom the camera toward its target      |
//! | `pan`         | `dx dy`                 | pan the camera in screen pixels        |
//!
//! # TCL-style scripting subset
//!
//! Beyond the fixed shape-construction table, this port implements a small
//! TCL-style scripting layer over the same line parser, covering the pieces
//! of `Draw_Interpretor` that scripts actually exercise: variables,
//! expressions, control flow, shape transforms and a soft-rendered view.
//!
//! **Variables.** `set <name> <value>` stores a string (or number) in the
//! session's variable namespace. Any token containing a `$name` reference is
//! expanded before dispatch, inside a token (`b$i`) or as a whole token
//! (`$x`); unknown names are left untouched so a literal dollar sign does not
//! break a command. [`expand_vars`] performs the substitution on a token
//! list and [`eval_number`] resolves a single number-or-variable.
//!
//! ```text
//! set x 3
//! box b $x 2 2          ; the box is 3 × 2 × 2
//! ```
//!
//! The namespace is introspectable and mutable the way a Tcl interpreter's is:
//! `vars` lists every `name = value`, `unset <name>` deletes one, and
//! `incr <name> [step]` adds `step` (default 1) to a numeric variable. The
//! `expr` result slot (`result`) is an ordinary entry in the same namespace,
//! so `incr result 1` after `expr 3 4 *` yields 13.
//!
//! ```text
//! set i 0
//! incr i 2                 ; i = 2
//! vars                     ; logs: vars: i = 2  ...
//! ```
//!
//! **Expressions.** `expr <a> <b> <op>` evaluates `a op b` (the operator is
//! the last token, matching the Draw REPL's postfix `expr 3 4 *`) with `op`
//! one of `+ - * /` on real numbers and stores the result in the `result`
//! variable. Division by zero is an error; results are stored in the compact
//! [`fmt_num`] form so `expr 3 4 *` gives `result = 12` rather than `12.0`.
//!
//! ```text
//! set a 5
//! expr $a 2 *            ; result = 10
//! box c $result 1 1      ; box c is 10 × 1 × 1
//! ```
//!
//! **Control flow.** `if <a> <b> <op> { <body> }` runs the body when the
//! numeric comparison (`== != < > <= >=`, postfix like `expr`) holds.
//! `for <var> <start> <end> { <body> }` runs the body once per integer in
//! `[start, end]` (descending when `start > end`), setting `var` before each
//! iteration so `$var` expands inside the body. The braces are optional — with
//! no `{` the body is every token after the header. Nested `if`s inside a
//! `for` body work, since each body is itself a command line through
//! [`execute_line`].
//!
//! ```text
//! for i 1 3 { box b$i 1 1 1 }    ; builds b1, b2, b3
//! for i 1 5 { if $i 3 > { echo hit$i } }   ; logs hit4, hit5
//! if $x 3 > { echo big }
//! ```
//!
//! **Shape transforms.** `translate`, `rotate`, `scale`, `copy` and
//! `transform` apply a `GpTrsf` to the *registered geometry* of a shape (the
//! [`crate::shape_ops`] kernels), so world coordinates — vertex points, edge
//! curves, face surfaces — actually move. This mirrors OCCT's `BRep_Tool::Transform`
//! and is distinct from merely relabelling the shape's location.
//!
//! ```text
//! translate b 1 0 0 c      ; b moved by (1,0,0) into c, b untouched
//! rotate b 0 0 1 90        ; b rotated 90° about the Z axis (in place)
//! scale b 2                ; b scaled ×2 about the origin (in place)
//! mirror b 0 0 1           ; b reflected across the Z = 0 plane (in place)
//! copy b b2                ; independent deep copy (b2 edits never touch b)
//! transform b 1 0 0 0 1 0 0 0 1 0 0 0   ; arbitrary affine, identity here
//! ```
//!
//! After a transform, `bbox <name>` (the axis-aligned bounding-box corners)
//! and `verts <name>` (every distinct vertex coordinate) inspect the *world*
//! geometry — a quick way to verify a `translate`/`rotate`/`scale` did what
//! you asked, complementing `info`'s counts and volume.
//!
//! **View.** `view <name> [w h] [file]` soft-renders a shape through the
//! session camera with the [`crate::viz_scene`] shaded ray-caster and writes
//! a binary PPM (default `view.ppm`, 256×256). The camera is a real
//! look-at state you can drive with `view_camera` (eye/target), `orbit`
//! (yaw/pitch degrees), `zoom` (factor) and `pan` (screen pixels) — a small
//! port of `V3d_View` interaction. The camera state lives on the session, so
//! a script can frame a shape once and render it at several resolutions.
//!
//! ```text
//! view_camera 5 4 5 0 0 0
//! orbit 20 10
//! zoom 1.5
//! view b 320 240           ; writes view.ppm (320×240)
//! ```
//!
//! `echo <text>` appends its (expanded) arguments to the log — the scripting
//! printf, useful both interactively and from `if`/`for` bodies.
//!
//! # Port notes on the scripting layer
//!
//! * Variable substitution is deliberately **non-destructive**: an unknown
//!   `$name` stays `$name` instead of raising, so file paths containing `$`
//!   survive and a typo surfaces at the command that consumes the value
//!   (a numeric parse error), not earlier.
//! * Bodies of `for`/`if` are token lists re-substituted per iteration. The
//!   outer [`execute_line`] already expanded whatever was bound before the
//!   loop; the loop variable is bound *after* that, so its `$var` reference
//!   survives until the per-iteration re-expansion. This is the same double
//!   evaluation a Tcl script sees with `{*}` bodies, minus Tcl's brace
//!   quoting subtleties.
//! * `for` counts in integers (`start`/`end` are truncated with `as i64`) and
//!   steps `−1` when `start > end`, so both `for i 1 3` and `for i 3 1` run
//!   three iterations. `incr` is the integer-typed counter; `expr` results
//!   may be fractional.
//! * The transforms rebuild geometry through the global [`crate::tgeometry`]
//!   registry (see [`crate::shape_ops`]), which is what makes `verts` and
//!   `bbox` report world coordinates afterwards. Copying is a deep
//!   [`translated_copy`] with the identity translation — independent `TShape`
//!   handles and re-registered geometry, so it is safe to fillet or transform
//!   the copy without affecting the source.
//!
//! A complete scripting session mixing variables, a loop, transforms and a
//! render might look like:
//!
//! ```text
//! # a parametric rack of three boxes with a sphere on top
//! set n 3
//! for i 1 $n { box b$i 2 2 2 }
//! translate b2 3 0 0
//! translate b3 6 0 0
//! scale b1 1.5
//! copy b3 b3s
//! mirror b3s 0 0 1
//! fuse b1 b2 rack1
//! fuse rack1 b3 rack2
//! info rack2
//! bbox rack2
//! if $n 2 >= { echo rack complete }
//! view rack2 320 240
//! ```
//!
//! which exercises every layer of the port in one script: `set`/`for`
//! variable binding, `$` substitution inside tokens (`b$i`), `translate`/
//! `scale`/`mirror`/`copy` transforms, boolean assembly, `info`/`bbox`
//! inspection, `if` + `echo`, and a soft-rendered PPM view.
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

use occt_core::gp::{GpAx1, GpDir, GpMat, GpPnt, GpTrsf, GpVec, GpXyz, TrsfForm};

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
use crate::shape_ops::{
    rotate_shape, scale_shape, transform_shape, translate_shape, transformed_copy, translated_copy,
};
use crate::step::write_shape_step;
use crate::topo_tools_full::{edges_of, shape_counts};
use crate::viz_scene::{render_scene_ppm_shaded, Camera, RenderSettings, SceneShape, VizScene};

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
// Variables, numbers and expressions
// ---------------------------------------------------------------------------

/// Render a real number the way the Tcl layer prints one: integral values
/// drop the `.0` (`3`, not `3.0`), everything else keeps the compact `f64`
/// rendering. This keeps `expr 3 4 *` storing `result = 12` and round-trips
/// through [`eval_number`].
fn fmt_num(v: f64) -> String {
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
fn expand_token(session: &DrawSession, tok: &str) -> String {
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
fn command_block<'a>(args: &'a [String], keyword: &str, header: usize) -> Result<&'a [String], String> {
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
         mesh step obj iges stl info bbox verts ls rm clear help exit quit \
         set expr echo vars unset incr for if translate rotate scale mirror copy transform \
         view view_camera orbit zoom pan"
            .into(),
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// Variables and expressions
// ---------------------------------------------------------------------------

/// `set <name> <value>` — store `value` (a string or a number) in the session
/// variable namespace. The value is stored verbatim after `$`-expansion, so
/// `set b $a` copies `a`. Later `$name` references expand anywhere in a
/// command line — whole tokens (`$x`) and embedded in longer tokens (`b$i`).
fn cmd_set(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_expr(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_echo(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    session.log.push(args.join(" "));
    Ok(())
}

/// `vars` — list every session variable as `name = value`, sorted by name, on
/// one log line. The analogue of Tcl's `info vars` and handy before a `for`
/// loop to confirm what the loop variable will shadow.
fn cmd_vars(session: &mut DrawSession) -> Result<(), String> {
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

/// `unset <name>` — remove the session variable `name`. Unknown names are an
/// error, matching Tcl's `unset` behaviour; [`eval_number`] on a name that no
/// longer exists falls through to a parse error rather than reading a stale
/// value.
fn cmd_unset(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_incr(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_for(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_if(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_translate(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_rotate(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_scale(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_copy(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_transform(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_mirror(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_view(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_view_camera(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_orbit(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 2, "orbit")?;
    let dx = eval_number(session, &args[0])?;
    let dy = eval_number(session, &args[1])?;
    session.camera.camera_orbit(dx, dy);
    session.log.push(format!("orbit: dx {dx} dy {dy}"));
    Ok(())
}

/// `zoom <factor>` — zoom the session camera by moving the eye toward its
/// target; `factor > 1` zooms in. Ports `V3d_View::SetZoom`.
fn cmd_zoom(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 1, "zoom")?;
    let factor = eval_number(session, &args[0])?;
    session.camera.camera_zoom(factor);
    session.log.push(format!("zoom: factor {factor}"));
    Ok(())
}

/// `pan <dx> <dy>` — pan the session camera's eye and target together along
/// the camera right/up plane, in screen pixels at a nominal 600-px viewport.
/// Ports `V3d_View::Pan`.
fn cmd_pan(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
    need_args(args, 2, "pan")?;
    let dx = eval_number(session, &args[0])?;
    let dy = eval_number(session, &args[1])?;
    session.camera.camera_pan(dx, dy);
    session.log.push(format!("pan: dx {dx} dy {dy}"));
    Ok(())
}

// ---------------------------------------------------------------------------
// Inspection
// ---------------------------------------------------------------------------

/// `bbox <name>` — print the axis-aligned bounding box of `name` as two corner
/// points, taken from the registered geometry (so transformed shapes report
/// their world bounds).
fn cmd_bbox(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
fn cmd_verts(session: &mut DrawSession, args: &[String]) -> Result<(), String> {
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
    let tokens = expand_vars(session, &tokens);
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
        "bbox" => cmd_bbox(session, command.args()),
        "verts" => cmd_verts(session, command.args()),
        "ls" => cmd_ls(session),
        "rm" => cmd_rm(session, command.args()),
        "clear" => cmd_clear(session),
        "help" => cmd_help(session),
        "set" => cmd_set(session, command.args()),
        "expr" => cmd_expr(session, command.args()),
        "echo" => cmd_echo(session, command.args()),
        "vars" => cmd_vars(session),
        "unset" => cmd_unset(session, command.args()),
        "incr" => cmd_incr(session, command.args()),
        "for" => cmd_for(session, command.args()),
        "if" => cmd_if(session, command.args()),
        "translate" => cmd_translate(session, command.args()),
        "rotate" => cmd_rotate(session, command.args()),
        "scale" => cmd_scale(session, command.args()),
        "mirror" => cmd_mirror(session, command.args()),
        "copy" => cmd_copy(session, command.args()),
        "transform" => cmd_transform(session, command.args()),
        "view" => cmd_view(session, command.args()),
        "view_camera" => cmd_view_camera(session, command.args()),
        "orbit" => cmd_orbit(session, command.args()),
        "zoom" => cmd_zoom(session, command.args()),
        "pan" => cmd_pan(session, command.args()),
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
    use crate::bbox_from_geometry::shape_bbox;
    use crate::brep_tool::BRepTool;
    use crate::topo_tools_full::{faces_of, vertices_of};

    /// (min, max) corners of a named shape's axis-aligned bbox.
    fn bbox_pts(s: &DrawSession, name: &str) -> ((f64, f64, f64), (f64, f64, f64)) {
        let shape = shape_owned(s, name).unwrap();
        let bb = shape_bbox(&shape);
        match bb.get() {
            Some((x0, x1, y0, y1, z0, z1)) => ((x0, y0, z0), (x1, y1, z1)),
            None => panic!("empty bbox for {name}"),
        }
    }

    /// Vertex coordinates of a named shape (world space).
    fn vertex_pts(s: &DrawSession, name: &str) -> Vec<(f64, f64, f64)> {
        let shape = shape_owned(s, name).unwrap();
        vertices_of(&shape)
            .iter()
            .map(|v| {
                let p = BRepTool::vertex_point(v);
                (p.x(), p.y(), p.z())
            })
            .collect()
    }

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-6
    }

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

    // -- TCL-style variables and expressions ---------------------------------

    #[test]
    fn set_and_expand_var() {
        let mut s = DrawSession::default();
        run_script(&mut s, "set x 3").expect("set");
        assert_eq!(s.vars.get("x").map(String::as_str), Some("3"));
        execute_line(&mut s, "box b $x 2 2").expect("box with $x");
        let (lo, hi) = bbox_pts(&s, "b");
        assert!(approx(hi.0, 3.0), "box x extent 3, got max {}", hi.0);
        assert!(approx(hi.1, 2.0) && approx(hi.2, 2.0));
        assert!(approx(lo.0, 0.0) && approx(lo.1, 0.0) && approx(lo.2, 0.0));
    }

    #[test]
    fn expr_arithmetic() {
        let mut s = DrawSession::default();
        execute_line(&mut s, "expr 3 4 *").expect("3*4");
        assert_eq!(s.vars.get("result").map(String::as_str), Some("12"));
        execute_line(&mut s, "expr 10 2 /").expect("10/2");
        assert_eq!(s.vars.get("result").map(String::as_str), Some("5"));
        // The stored value round-trips through eval_number.
        let v = eval_number(&s, "$result").expect("parse result");
        assert!(approx(v, 5.0), "result value {v}");
    }

    #[test]
    fn expr_var_participation() {
        let mut s = DrawSession::default();
        run_script(&mut s, "set a 5").expect("set a");
        execute_line(&mut s, "expr $a 2 +").expect("5+2");
        assert_eq!(s.vars.get("result").map(String::as_str), Some("7"));
    }

    #[test]
    fn expr_consumed_by_shape() {
        let mut s = DrawSession::default();
        run_script(&mut s, "expr 3 4 *\nbox c $result 1 1").expect("expr + box");
        let (_, hi) = bbox_pts(&s, "c");
        assert!(approx(hi.0, 12.0), "box from expr result, x max {}", hi.0);
    }

    // -- Shape transformations ------------------------------------------------

    #[test]
    fn translate_shape() {
        let mut s = DrawSession::default();
        run_script(&mut s, "box b 1 1 1").expect("box");
        execute_line(&mut s, "translate b 1 0 0").expect("translate");
        let (lo, hi) = bbox_pts(&s, "b");
        assert!(approx(lo.0, 1.0) && approx(hi.0, 2.0), "x shifted by +1: {lo:?} {hi:?}");
        assert!(approx(lo.1, 0.0) && approx(hi.1, 1.0));
    }

    #[test]
    fn translate_to_out_keeps_original() {
        let mut s = DrawSession::default();
        run_script(&mut s, "box b 1 1 1").expect("box");
        execute_line(&mut s, "translate b 5 0 0 c").expect("translate to c");
        let (lo_b, hi_b) = bbox_pts(&s, "b");
        let (lo_c, _) = bbox_pts(&s, "c");
        assert!(approx(lo_b.0, 0.0) && approx(hi_b.0, 1.0), "b untouched");
        assert!(approx(lo_c.0, 5.0), "c moved");
    }

    #[test]
    fn rotate_shape_90() {
        let mut s = DrawSession::default();
        run_script(&mut s, "box b 1 1 1").expect("box");
        execute_line(&mut s, "rotate b 0 0 1 90").expect("rotate 90");
        let pts = vertex_pts(&s, "b");
        let has_corner = pts
            .iter()
            .any(|&(x, y, z)| approx(x, 0.0) && approx(y, 1.0) && approx(z, 0.0));
        assert!(has_corner, "corner (1,0,0) rotated to (0,1,0), pts {pts:?}");
    }

    #[test]
    fn scale_shape() {
        let mut s = DrawSession::default();
        run_script(&mut s, "box b 2 2 2").expect("box");
        execute_line(&mut s, "scale b 2").expect("scale");
        let (lo, hi) = bbox_pts(&s, "b");
        assert!(approx(lo.0, 0.0) && approx(hi.0, 4.0), "bbox doubled: {hi:?}");
        assert!(approx(hi.1, 4.0) && approx(hi.2, 4.0));
    }

    #[test]
    fn copy_shape_is_independent() {
        let mut s = DrawSession::default();
        run_script(&mut s, "box b 1 1 1\ncopy b b2").expect("copy");
        assert!(s.shapes.contains_key("b2"), "b2 registered");
        execute_line(&mut s, "translate b2 3 0 0").expect("move b2");
        let (lo_b, hi_b) = bbox_pts(&s, "b");
        let (lo_b2, hi_b2) = bbox_pts(&s, "b2");
        assert!(approx(lo_b.0, 0.0) && approx(hi_b.0, 1.0), "b unchanged");
        assert!(approx(lo_b2.0, 3.0) && approx(hi_b2.0, 4.0), "b2 moved");
    }

    #[test]
    fn transform_affine_identity_is_noop() {
        let mut s = DrawSession::default();
        run_script(&mut s, "box b 1 2 3").expect("box");
        execute_line(&mut s, "transform b 1 0 0 0 1 0 0 0 1 0 0 0").expect("identity transform");
        let (lo, hi) = bbox_pts(&s, "b");
        assert!(approx(lo.0, 0.0) && approx(lo.1, 0.0) && approx(lo.2, 0.0));
        assert!(approx(hi.0, 1.0) && approx(hi.1, 2.0) && approx(hi.2, 3.0));
    }

    #[test]
    fn mirror_across_z0_plane() {
        let mut s = DrawSession::default();
        run_script(&mut s, "box b 1 1 1").expect("box");
        execute_line(&mut s, "mirror b 0 0 1").expect("mirror z=0");
        // Check the vertex geometry, which the mirror transform rewrites: the
        // top corners (z = 1) must move to z = −1 and no corner may remain at
        // z > 0. (shape_bbox is not used here because its face-surface UV
        // sampling is unreliable on orientation-flipped mirrored planes.)
        let pts = vertex_pts(&s, "b");
        let has_neg_z = pts.iter().any(|&(_, _, z)| approx(z, -1.0));
        let no_pos_z = pts.iter().all(|&(_, _, z)| z <= 1e-6);
        assert!(has_neg_z, "a corner moved to z = -1, pts {pts:?}");
        assert!(no_pos_z, "no corner remains at z > 0, pts {pts:?}");
    }

    // -- View commands ---------------------------------------------------------

    #[test]
    fn view_writes_ppm() {
        let mut s = DrawSession::default();
        run_script(&mut s, "box b 2 2 2").expect("box");
        let dir = std::env::temp_dir();
        let path = dir.join(format!("draw_view_{}.ppm", std::process::id()));
        let path_s = path.to_str().unwrap().to_string();
        execute_line(&mut s, &format!("view b 64 48 {path_s}")).expect("view");
        let data = std::fs::read(&path).expect("ppm exists");
        assert!(data.starts_with(b"P6"), "ppm header, first bytes {:?}", &data[..3.min(data.len())]);
        assert!(
            data.len() >= 13 + 64 * 48 * 3,
            "expected >= {} bytes, got {}",
            13 + 64 * 48 * 3,
            data.len()
        );
        assert!(s.log.iter().any(|l| l.contains("(64x48)")), "log: {:?}", s.log);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn view_camera_orbit() {
        let mut s = DrawSession::default();
        run_script(&mut s, "box b 2 2 2").expect("box");
        execute_line(&mut s, "view_camera 5 4 5 0 0 0").expect("view_camera");
        let after_set = s.camera.eye;
        execute_line(&mut s, "orbit 20 10").expect("orbit");
        assert!(
            s.camera.eye.distance(&after_set) > 1e-6,
            "orbit moved the eye"
        );
        execute_line(&mut s, "zoom 1.5").expect("zoom");
        execute_line(&mut s, "pan 5 -3").expect("pan");
        let dir = std::env::temp_dir();
        let path = dir.join(format!("draw_view_orbit_{}.ppm", std::process::id()));
        let path_s = path.to_str().unwrap().to_string();
        execute_line(&mut s, &format!("view b 32 24 {path_s}")).expect("view");
        let data = std::fs::read(&path).expect("ppm exists");
        assert!(data.starts_with(b"P6"), "view after camera moves still writes");
        let _ = std::fs::remove_file(&path);
    }

    // -- Control flow ----------------------------------------------------------

    #[test]
    fn for_loop_boxes() {
        let mut s = DrawSession::default();
        execute_line(&mut s, "for i 1 3 { box b$i 1 1 1 }").expect("for loop");
        assert!(s.shapes.contains_key("b1"), "b1 built");
        assert!(s.shapes.contains_key("b2"), "b2 built");
        assert!(s.shapes.contains_key("b3"), "b3 built");
        assert!(!s.shapes.contains_key("b4"), "loop stopped at end");
        // The loop variable is left bound.
        assert_eq!(s.vars.get("i").map(String::as_str), Some("3"));
    }

    #[test]
    fn if_comparison() {
        let mut s = DrawSession::default();
        run_script(&mut s, "set x 5").expect("set x");
        execute_line(&mut s, "if $x 3 > { echo big }").expect("if true");
        execute_line(&mut s, "if $x 3 < { echo small }").expect("if false");
        let joined = s.log.join("\n");
        assert!(joined.contains("big"), "log: {joined}");
        assert!(!joined.contains("small"), "false branch skipped: {joined}");
    }

    #[test]
    fn for_loop_uses_if_inside() {
        let mut s = DrawSession::default();
        run_script(&mut s, "for i 1 5 { if $i 3 > { echo hit$i } }").expect("for+if");
        let joined = s.log.join("\n");
        assert!(joined.contains("hit4"), "log: {joined}");
        assert!(joined.contains("hit5"), "log: {joined}");
        assert!(!joined.contains("hit3"), "3 is not > 3: {joined}");
    }

    #[test]
    fn echo_expands() {
        let mut s = DrawSession::default();
        run_script(&mut s, "set name part1").expect("set name");
        execute_line(&mut s, "echo building $name").expect("echo");
        let joined = s.log.join("\n");
        assert!(joined.contains("building part1"), "log: {joined}");
    }

    #[test]
    fn incr_unset_vars() {
        let mut s = DrawSession::default();
        run_script(&mut s, "set i 0").expect("set i");
        execute_line(&mut s, "incr i").expect("incr by 1");
        assert_eq!(s.vars.get("i").map(String::as_str), Some("1"));
        execute_line(&mut s, "incr i 4").expect("incr by 4");
        assert_eq!(s.vars.get("i").map(String::as_str), Some("5"));
        execute_line(&mut s, "vars").expect("vars");
        let joined = s.log.join("\n");
        assert!(joined.contains("i = 5"), "vars lists the counter: {joined}");
        execute_line(&mut s, "unset i").expect("unset");
        assert!(!s.vars.contains_key("i"), "i removed");
        let err = execute_line(&mut s, "unset i").expect_err("double unset errors");
        assert!(err.contains("not found"), "err {err}");
    }

    // -- Inspection commands added alongside the scripting layer ----------------

    #[test]
    fn bbox_and_verts_commands() {
        let mut s = DrawSession::default();
        run_script(&mut s, "box b 2 1 1").expect("box");
        execute_line(&mut s, "bbox b").expect("bbox");
        let joined = s.log.join("\n");
        assert!(joined.contains("min (0, 0, 0)"), "bbox min: {joined}");
        assert!(joined.contains("max (2, 1, 1)"), "bbox max: {joined}");
        execute_line(&mut s, "verts b").expect("verts");
        let joined = s.log.join("\n");
        assert!(joined.contains("8 vertices"), "verts count: {joined}");
    }
}
