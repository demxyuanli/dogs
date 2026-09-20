//! Draw_Interpretor-lite: a scripting + batch command driver for the BRep kernel.
//!
//! **UNPORTED (audit A14)**: a port-local, dependency-free command driver; OCCT's
//! Draw lives in TKDraw (`Draw_Interpretor`, `Draw_Commands`, `DBRep`) and is not
//! translated here.
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
//! | `help`        | —                    | list commands with arity and description          |
//! | `commands`    | —                    | list all command names                            |
//! | `arity`       | `cmd`                | report a command's accepted argument count        |
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
//! | `proc` / `def`| `name { params } { body }`| define a procedure                  |
//! | `call`        | `name args...`          | invoke a procedure by name             |
//! | `return`      | `[value]`               | stop the current procedure, set `result`|
//! | `lappend`     | `name value...`         | append values to a list/array variable |
//! | `len`         | `name`                  | count list elements into `result`      |
//! | `concat`      | `a b out`               | concatenate two list variables         |
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
//! **Procedures.** `proc <name> { p1 p2 ... } { body }` (alias `def`) stores a
//! `;`-separated command sequence under `name`. Calling it — `call <name>
//! args...` or directly by name — binds each parameter to the matching
//! argument and runs the body through [`execute_line`], so parameters shadow
//! global variables for the duration of the call and are restored afterwards.
//! `return [value]` stops the body (storing `value` in `result`), and `expr`
//! leaves its result in `result`, so a procedure's value is read back through
//! that variable. Bodies may call other procedures — including themselves, so
//! recursion works (bounded by [`PROC_MAX_DEPTH`]). An error inside a body is
//! reported with the procedure name and the call stack.
//!
//! ```text
//! proc fact { n } { if $n 2 < { return $n } ; expr $n 1 - ; set m $result ; call fact $m ; expr $n $result * }
//! call fact 5                          ; result = 120
//! ```
//!
//! **Arrays and lists.** A variable may be indexed as an array element —
//! `set arr(0) 5` stores the key `arr(0)`, and `$arr(0)` resolves it — or hold
//! a space-separated list built with `lappend`. `len <name>` counts list
//! elements into `result`; `concat <a> <b> <out>` joins two list variables.
//!
//! ```text
//! lappend pts 1.0
//! lappend pts 2.0 3.0            ; pts = "1.0 2.0 3.0"
//! len pts                       ; result = 3
//! set arr(0) 10 ; set arr(1) 20 ; $arr(0) + $arr(1) = 30
//! ```
//!
//! **Command registry.** Every command lives in a centralized table (name →
//! handler + description + arity). `help` prints the whole table, `commands`
//! lists the names, and `arity <cmd>` queries a command's accepted argument
//! count — introspection that makes a session self-documenting, the analogue
//! of OCCT's `Draw_Interpretor::PrintCommands` and `PrintHelp`.
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
//!   dispatches to `crate::bop_builder::boolean` (`BOPAlgo_BOP`) when both
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
mod prelude {

pub(crate) use std::collections::HashMap;
pub(crate) use std::io::{BufRead, Write};

pub(crate) use occt_core::gp::{GpAx1, GpDir, GpMat, GpPnt, GpTrsf, GpVec, GpXyz, TrsfForm};

pub(crate) use crate::abs::ShapeType;
pub(crate) use crate::bop_builder::BoolOp;
pub(crate) use crate::bop_curved::curved_boolean_ext;
pub(crate) use crate::brep_exchange::{brep_to_obj, brep_to_stl_binary};
pub(crate) use crate::brep_offset::offset_shell;
pub(crate) use crate::fillet_edge::fillet_edge;
pub(crate) use crate::iges::write_shape_iges;
pub(crate) use crate::primitives::{
    BRepPrimBox, BRepPrimCone, BRepPrimCylinder, BRepPrimSphere, BRepPrimTorus,
};
pub(crate) use crate::shape::TopoShape;
pub(crate) use crate::shape_mesh::{mesh_shape, shape_volume};
pub(crate) use crate::shape_ops::{
    rotate_shape, scale_shape, transform_shape, translate_shape, transformed_copy, translated_copy,
};
pub(crate) use crate::step::write_shape_step;
pub(crate) use crate::topo_tools_full::{edges_of, shape_counts};
pub(crate) use crate::viz_scene::{render_scene_ppm_shaded, Camera, RenderSettings, SceneShape, VizScene};

}


mod p01;
mod p02;
mod p03;
pub use p01::*;
pub use p02::*;
pub use p03::*;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
