//! Structural self-consistency check for the IGES writer (no gate covers the
//! IGES text, so this example encodes the invariants the writer must satisfy).
//!
//! Run: `cargo run --manifest-path crates/occt-topo/Cargo.toml --offline \
//!       --example iges_check -- Cube Sphere Shape Shape-2 HoledPlate`
//! (arguments are `data/<name>.step`; a path containing `/` is used as-is).
//!
//! Checked, per model:
//!  1. every card is 80 columns;
//!  2. the P section's sequence numbers are contiguous `1..N`, each P line's
//!     directory pointer is `2i - 1` of the entity it belongs to
//!     (`IGESData_IGESWriter.cxx:903`), and every DE's `pstart`/`pcount`
//!     (`cxx:834-835`) addresses that entity's own parameter cards;
//!  3. the leading **pointer fields** of each composite entity resolve to
//!     existing directories (102 curve list, 142 surface/curve3d, 144
//!     surface/outer+inner, 402 entity list, 192/194/196/198 point/axis/refdir,
//!     120 axis/generatrix, 122 directrix);
//!  4. the Terminate card's four counts equal the actual section sizes
//!     (`cxx:942-947`: `nbs`, `nbg`, `nbd * 2`, last P sequence).
use std::collections::{BTreeMap, HashSet};
use std::path::PathBuf;

use occt_topo::builder::TopoBuilder;
use occt_topo::iges::write_shape_iges;
use occt_topo::step::read_step_file;

fn col(l: &str, a: usize, b: usize) -> String {
    l.get(a..b).unwrap_or("").trim().to_string()
}

fn main() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..").join("..");
    let names: Vec<String> = std::env::args().skip(1).collect();
    let names = if names.is_empty() {
        vec!["Cube".to_string(), "Sphere".to_string()]
    } else {
        names
    };
    let mut bad = 0usize;
    for name in names {
        let path = if name.contains('/') {
            root.join(&name)
        } else {
            root.join("data").join(format!("{name}.step"))
        };
        let model = match read_step_file(&path.to_string_lossy()) {
            Ok(m) => m,
            Err(e) => {
                println!("ERR {name}: {e}");
                bad += 1;
                continue;
            }
        };
        let shape = if model.shapes.len() == 1 {
            model.shapes[0].shape.clone()
        } else {
            let parts: Vec<_> = model.shapes.iter().map(|s| s.shape.clone()).collect();
            TopoBuilder::new().make_compound_of(&parts).0
        };
        let iges = write_shape_iges(&shape);
        let mut problems: Vec<String> = Vec::new();
        let mut counts: BTreeMap<char, usize> = BTreeMap::new();
        let mut p_by_seq: BTreeMap<usize, String> = BTreeMap::new();
        let mut p_de_pointer: BTreeMap<usize, usize> = BTreeMap::new();
        // (type, pstart, pcount, trsf, de)
        let mut des: Vec<(i32, usize, usize, usize, usize)> = Vec::new();
        let mut pending: Option<(i32, usize, usize, usize)> = None;
        let mut t_card: Option<String> = None;
        for (lineno, l) in iges.lines().enumerate() {
            if l.len() != 80 {
                problems.push(format!("line {} is {} cols", lineno + 1, l.len()));
                break;
            }
            let sec = l.as_bytes()[72] as char;
            *counts.entry(sec).or_insert(0) += 1;
            let seq: usize = col(l, 73, 80).parse().unwrap_or(0);
            match sec {
                'P' => {
                    p_by_seq.insert(seq, l[..64].to_string());
                    p_de_pointer.insert(seq, col(l, 65, 72).parse().unwrap_or(0));
                }
                'D' => {
                    if seq % 2 == 1 {
                        pending = Some((
                            col(l, 0, 8).parse().unwrap_or(-1),
                            col(l, 8, 16).parse().unwrap_or(0),
                            col(l, 48, 56).parse().unwrap_or(0),
                            (seq + 1) / 2,
                        ));
                    } else if let Some((ty, pstart, trsf, de)) = pending.take() {
                        let pcount: usize = col(l, 24, 32).parse().unwrap_or(0);
                        des.push((ty, pstart, pcount, trsf, de));
                    }
                }
                'T' => t_card = Some(l.to_string()),
                _ => {}
            }
        }
        // (2) P bookkeeping
        if let Some(max) = p_by_seq.keys().max() {
            if *max != p_by_seq.len() {
                problems.push(format!("P sequence max {max} != count {}", p_by_seq.len()));
            }
        }
        for (ty, pstart, pcount, _trsf, de) in &des {
            if *pstart == 0 || *pcount == 0 {
                problems.push(format!("DE{de} ty={ty} pstart={pstart} pcount={pcount}"));
                continue;
            }
            if p_de_pointer.get(pstart) != Some(&(2 * de - 1)) {
                problems.push(format!(
                    "DE{de} ty={ty}: P{pstart} points to {:?}",
                    p_de_pointer.get(pstart)
                ));
            }
            for k in 0..*pcount {
                if !p_by_seq.contains_key(&(pstart + k)) {
                    problems.push(format!("DE{de} ty={ty} missing P seq {}", pstart + k));
                    break;
                }
            }
        }
        // (3) pointer fields
        let de_set: HashSet<usize> = des.iter().map(|d| d.4).collect();
        // Directories nothing points at (roots are legitimately unreferenced).
        let mut referenced: HashSet<usize> = HashSet::new();
        for (ty, pstart, pcount, _trsf, de) in &des {
            let mut body = String::new();
            for k in 0..*pcount {
                if let Some(c) = p_by_seq.get(&(pstart + k)) {
                    body.push_str(c.trim_end());
                }
            }
            let f: Vec<&str> = body.trim_end_matches(';').split(',').collect();
            // The pointer fields of this type, by record position - the same
            // layout the writer's `emit_refs` records (`Entries in OwnShared`
            // order). Collecting them lets the check both validate the numbers and
            // report which directories nothing references (T-85's motivation: OCCT
            // only writes entities reachable from the root, this writer writes
            // everything it created).
            let mut ptr_idx: Vec<usize> = Vec::new();
            match ty {
                102 => {
                    if let Some(n) = f.get(1).and_then(|s| s.trim().parse::<usize>().ok()) {
                        ptr_idx.extend((0..n).map(|k| 2 + k));
                    }
                }
                // 142 is "142,creation_mode,#surface,#curve_uv,#curve_3d,preference;"
                // (IGESGeom_ToolCurveOnSurface::WriteOwnParams), so the surface,
                // the 2-D (UV) curve and the 3-D curve are at fields 2, 3 and 4.
                // Reading only [2, 4] left every emitted UV curve looking
                // unreferenced (same class as the T-85 pointer-table fix above).
                142 => ptr_idx.extend([2usize, 3, 4]),
                144 => {
                    if let Some(n) = f.get(3).and_then(|s| s.trim().parse::<usize>().ok()) {
                        ptr_idx.push(4);
                        ptr_idx.extend((0..n).map(|k| 5 + k));
                    }
                }
                402 => {
                    if let Some(n) = f.get(1).and_then(|s| s.trim().parse::<usize>().ok()) {
                        ptr_idx.extend((0..n).map(|k| 2 + k));
                    }
                }
                // `192/194/198` are `"<ty>,#0,#1,{radius|params},#2;"`, so their
                // DE pointers sit at fields 1, 2 and the last one.
                192 | 194 | 198 => ptr_idx.extend([1usize, 2, f.len() - 1]),
                // `196` (sphere) is `"196,#0,{radius},#1,#2;"` — the radius comes
                // *between* the location and the two directions
                // (`IGESGeom_ToolSphere::WriteOwnParams`, mirrored by
                // `occt-topo/src/iges.rs::emit_refs(196, …)`), so the axis
                // direction is at field 3, not 2. Reading it as `[1, 2, last]`
                // left the axis looking unreferenced (T-85 remainder: Sphere
                // reported `{123: 1, 144: 1}` when only the root is an orphan).
                196 => ptr_idx.extend([1usize, 3, f.len() - 1]),
                120 => ptr_idx.extend([1usize, 2]),
                122 => ptr_idx.push(1),
                _ => {}
            }
            for idx in ptr_idx {
                if let Some(txt) = f.get(idx) {
                    if let Ok(p) = txt.trim().parse::<usize>() {
                        if p == 0 {
                            continue;
                        }
                        referenced.insert(p);
                        if !de_set.contains(&p) {
                            problems.push(format!("DE{de} ty={ty} field {idx} -> {p} (no such DE)"));
                        }
                    }
                }
            }
        }
        // (4) T card
        match &t_card {
            Some(t) => {
                let s: usize = col(t, 1, 8).parse().unwrap_or(0);
                let g: usize = col(t, 9, 16).parse().unwrap_or(0);
                let d: usize = col(t, 17, 24).parse().unwrap_or(0);
                let p: usize = col(t, 25, 32).parse().unwrap_or(0);
                let actual = |c: char| *counts.get(&c).unwrap_or(&0);
                if (s, g, d, p) != (actual('S'), actual('G'), actual('D'), actual('P')) {
                    problems.push(format!(
                        "T card S{s}G{g}D{d}P{p} vs actual S{}G{}D{}P{}",
                        actual('S'),
                        actual('G'),
                        actual('D'),
                        actual('P')
                    ));
                }
            }
            None => problems.push("no Terminate card".into()),
        }
        let orphans: Vec<(usize, i32)> = des
            .iter()
            .filter(|(_, _, _, _, de)| !referenced.contains(de))
            .map(|(ty, _, _, _, de)| (*de, *ty))
            .collect();
        if problems.is_empty() {
            println!(
                "ok  {name}: DE={} P={} sections={counts:?} unreferenced={}",
                des.len(),
                p_by_seq.len(),
                orphans.len()
            );
            if !orphans.is_empty() {
                let mut by_ty: BTreeMap<i32, usize> = BTreeMap::new();
                for (_, ty) in &orphans {
                    *by_ty.entry(*ty).or_insert(0) += 1;
                }
                println!("    unreferenced by type: {by_ty:?} (first DEs {:?})", &orphans.iter().take(6).map(|(d, t)| (*d, *t)).collect::<Vec<_>>());
            }
        } else {
            bad += 1;
            println!("BAD {name}: {} problem(s)", problems.len());
            for p in problems.iter().take(6) {
                println!("    {p}");
            }
        }
    }
    if bad > 0 {
        std::process::exit(1);
    }
}
