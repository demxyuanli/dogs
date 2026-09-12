use super::prelude::*;
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

    // -- Phase 12: procedures, arrays and the command registry ------------------

    #[test]
    fn proc_define_and_call() {
        let mut s = DrawSession::default();
        run_script(&mut s, "proc add { a b } { expr $a $b + }").expect("define add");
        assert!(s.procs.contains_key("add"), "add registered");
        execute_line(&mut s, "call add 3 4").expect("call add");
        assert_eq!(s.vars.get("result").map(String::as_str), Some("7"));
        // Direct-name invocation is a shorthand for `call`.
        execute_line(&mut s, "add 10 2").expect("direct call");
        assert_eq!(s.vars.get("result").map(String::as_str), Some("12"));
    }

    #[test]
    fn def_alias_and_direct_call() {
        let mut s = DrawSession::default();
        run_script(&mut s, "def sq { x } { expr $x $x * }").expect("def");
        execute_line(&mut s, "sq 9").expect("direct call by name");
        assert_eq!(s.vars.get("result").map(String::as_str), Some("81"));
    }

    #[test]
    fn proc_body_verbatim_survives_global_collision() {
        let mut s = DrawSession::default();
        // A global `a` must not be substituted into the body at definition time;
        // the body is stored verbatim and `$a` resolves to the parameter.
        run_script(&mut s, "set a 999").expect("set a");
        run_script(&mut s, "proc add { a b } { expr $a $b + }").expect("define add");
        execute_line(&mut s, "call add 2 3").expect("call");
        assert_eq!(s.vars.get("result").map(String::as_str), Some("5"));
    }

    #[test]
    fn proc_param_restored_after_call() {
        let mut s = DrawSession::default();
        run_script(&mut s, "set n 100").expect("set n");
        run_script(&mut s, "proc f { n } { expr $n 1 + }").expect("define f");
        execute_line(&mut s, "call f 5").expect("call f");
        assert_eq!(s.vars.get("result").map(String::as_str), Some("6"));
        // The parameter shadowed the global `n` only for the duration of the call.
        assert_eq!(s.vars.get("n").map(String::as_str), Some("100"), "global n restored");
    }

    #[test]
    fn proc_recursion_factorial() {
        let mut s = DrawSession::default();
        let script = concat!(
            "proc fact { n } { ",
            "if $n 2 < { return $n }; ",
            "expr $n 1 -; set m $result; ",
            "call fact $m; expr $n $result * ",
            "}\n",
            "call fact 5\n",
        );
        run_script(&mut s, script).expect("factorial script");
        assert_eq!(s.vars.get("result").map(String::as_str), Some("120"));
    }

    #[test]
    fn proc_nested_calls_with_parameters() {
        let mut s = DrawSession::default();
        let script = concat!(
            "proc double { x } { expr $x 2 * }\n",
            "proc add { a b } { expr $a $b + }\n",
            "proc quad { x } { call double $x; set a $result; call double $a; expr $a $result + }\n",
            "call quad 3\n",
        );
        run_script(&mut s, script).expect("nested script");
        // quad(3) = double(3) + double(double(3)) = 6 + 12 = 18.
        assert_eq!(s.vars.get("result").map(String::as_str), Some("18"));
        // A helper with disjoint parameter names still works.
        run_script(&mut s, "call add 2 3").expect("call add");
        assert_eq!(s.vars.get("result").map(String::as_str), Some("5"));
    }

    #[test]
    fn proc_error_reports_call_stack() {
        let mut s = DrawSession::default();
        let script = concat!(
            "proc inner { } { bogus_cmd }\n",
            "proc outer { } { call inner }\n",
            "call outer\n",
        );
        let err = run_script(&mut s, script).expect_err("inner fails");
        assert!(err.contains("in proc 'inner'"), "proc name: {err}");
        assert!(err.contains("outer -> inner"), "call stack in {err}");
        assert!(err.contains("unknown command 'bogus_cmd'"), "root cause: {err}");
    }

    #[test]
    fn proc_recursion_depth_limit() {
        let mut s = DrawSession::default();
        run_script(&mut s, "proc loop { n } { call loop $n }").expect("define loop");
        let err = execute_line(&mut s, "call loop 0").expect_err("runaway recursion");
        assert!(err.contains("max recursion depth"), "err {err}");
    }

    #[test]
    fn return_at_top_level_sets_result() {
        let mut s = DrawSession::default();
        execute_line(&mut s, "return 42").expect("top-level return");
        assert_eq!(s.vars.get("result").map(String::as_str), Some("42"));
        // No procedure is active, so no error signal leaks to the caller.
        assert!(s.last_error.is_none(), "no error recorded");
    }

    #[test]
    fn array_elements_add_update_query_remove() {
        let mut s = DrawSession::default();
        run_script(&mut s, "set arr(0) 10\nset arr(1) 20").expect("set array elements");
        assert_eq!(s.vars.get("arr(0)").map(String::as_str), Some("10"));
        // `$arr(i)` resolves the whole element name.
        execute_line(&mut s, "expr $arr(0) $arr(1) +").expect("array expr");
        assert_eq!(s.vars.get("result").map(String::as_str), Some("30"));
        // Update an element.
        execute_line(&mut s, "set arr(0) 5").expect("update element");
        execute_line(&mut s, "expr $arr(0) 2 *").expect("expr after update");
        assert_eq!(s.vars.get("result").map(String::as_str), Some("10"));
        // Remove an element.
        execute_line(&mut s, "unset arr(1)").expect("unset element");
        assert!(!s.vars.contains_key("arr(1)"), "element removed");
        assert!(s.vars.contains_key("arr(0)"), "other element intact");
    }

    #[test]
    fn lappend_len_concat() {
        let mut s = DrawSession::default();
        execute_line(&mut s, "lappend list a").expect("lappend create");
        execute_line(&mut s, "lappend list b c").expect("lappend append");
        assert_eq!(s.vars.get("list").map(String::as_str), Some("a b c"));
        execute_line(&mut s, "len list").expect("len");
        assert_eq!(s.vars.get("result").map(String::as_str), Some("3"));
        execute_line(&mut s, "lappend list d").expect("lappend again");
        execute_line(&mut s, "len list").expect("len again");
        assert_eq!(s.vars.get("result").map(String::as_str), Some("4"));
        // len of an undefined variable is 0.
        execute_line(&mut s, "len nope").expect("len undefined");
        assert_eq!(s.vars.get("result").map(String::as_str), Some("0"));
        // concat joins two list variables into a third.
        run_script(&mut s, "set l1 \"1 2\"\nset l2 \"3 4\"").expect("set lists");
        execute_line(&mut s, "concat l1 l2 out").expect("concat");
        assert_eq!(s.vars.get("out").map(String::as_str), Some("1 2 3 4"));
        // lappend also works on an array element.
        execute_line(&mut s, "lappend arr(0) 7").expect("lappend element");
        assert_eq!(s.vars.get("arr(0)").map(String::as_str), Some("7"));
    }

    #[test]
    fn commands_lists_every_command() {
        let mut s = DrawSession::default();
        run_script(&mut s, "proc userproc { } { echo hi }").expect("define userproc");
        execute_line(&mut s, "commands").expect("commands");
        let joined = s.log.join("\n");
        for cmd in [
            "box", "fuse", "proc", "def", "call", "return", "lappend", "len",
            "concat", "commands", "arity", "help", "exit", "userproc",
        ] {
            assert!(joined.contains(cmd), "commands lists {cmd}: {joined}");
        }
    }

    #[test]
    fn arity_queries() {
        let mut s = DrawSession::default();
        execute_line(&mut s, "arity box").expect("arity box");
        assert!(s.log.iter().any(|l| l.contains("arity box: 4")), "log: {:?}", s.log);
        execute_line(&mut s, "arity translate").expect("arity translate");
        assert!(s.log.iter().any(|l| l.contains("arity translate: 4..5")), "log: {:?}", s.log);
        execute_line(&mut s, "arity echo").expect("arity echo");
        assert!(s.log.iter().any(|l| l.contains("arity echo: 0 or more")), "log: {:?}", s.log);
        let err = execute_line(&mut s, "arity nope").expect_err("unknown arity");
        assert!(err.contains("unknown command 'nope'"), "err {err}");
    }

    #[test]
    fn help_lists_descriptions_and_user_procs() {
        let mut s = DrawSession::default();
        run_script(&mut s, "proc helper { } { echo hi }").expect("define helper");
        execute_line(&mut s, "help").expect("help");
        let joined = s.log.join("\n");
        assert!(joined.contains("build an axis-aligned box"), "help description: {joined}");
        assert!(joined.contains("helper"), "help lists user proc: {joined}");
        assert!(joined.contains("(proc)"), "proc marker: {joined}");
    }
