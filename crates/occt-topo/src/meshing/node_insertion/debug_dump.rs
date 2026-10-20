//! TEMPORARY diagnostics used to replay OCCT's Delaun dumps.

use std::collections::HashMap;

use occt_core::gp::{GpPnt, GpPnt2d, GpXY};
use occt_core::precision::SQUARE_CONFUSION;
use occt_geom::geom_api::project_point_on_surface;
use occt_geom::Surface;

use crate::abs::Orientation;

use super::super::data_model::*;
use super::super::delaun::Delaun;
use super::super::delaun_data::DelaunDataStructure;
use super::super::delaun_types::{DelaunLink, DelaunTriangle, DelaunVertex, VertexState};
use super::super::face_discret::Classifier;
#[allow(unused_imports)]
use super::super::mesh_algo::*;
use super::super::face_discret::PointState;
use super::super::deflection_control::{
    DelaunayDeflectionControlMeshAlgo, DeflectionMesh,
};
use super::super::parameters::MeshParameters;
use super::super::range_splitter::{
    mesh_surface_type, create_range_splitter, factory_uses_deflection_control,
    grab_params_of_internal_edges, RangeSplitter, MeshSurfaceType,
};
use super::*;

use super::*;

impl DelaunayNodeInsertionMeshAlgo {
    /// TEMPORARY diagnostic (DELAUNDUMP): structure dump in the text format
    /// consumed by `occt_probe --delaunstruct <dir>` (`delaun_in_f<face>.txt`):
    /// `TOL`, `CELLS`, `N i u v movability`, `L i first last movability`, `V ...`.
    pub(super) fn dump_delaun_structure(
        structure: &DelaunDataStructure,
        dir: &str,
        face_index: usize,
        cells_u: i32,
        cells_v: i32,
        indices: &[i32],
        free_xy: &[(f64, f64)],
    ) {
        use std::fmt::Write as _;
        let (tol_u, tol_v) = structure.tolerance();
        let mut out = String::new();
        let _ = writeln!(out, "TOL {tol_u:.17e} {tol_v:.17e}");
        let _ = writeln!(out, "CELLS {cells_u} {cells_v}");
        // 3D bbox of the registered nodes: pairs this dump with the probe's
        // `--boundary` / `--facestats` face list (whose face order differs from
        // the port's, so the index alone is not a pairing key).
        {
            let (mut lo, mut hi) = (GpPnt::zero(), GpPnt::zero());
            let mut first = true;
            for i in 1..=structure.nb_nodes() as i32 {
                let p = structure.node(i).p3d;
                if first {
                    lo = p;
                    hi = p;
                    first = false;
                } else {
                    lo = GpPnt::new(lo.x().min(p.x()), lo.y().min(p.y()), lo.z().min(p.z()));
                    hi = GpPnt::new(hi.x().max(p.x()), hi.y().max(p.y()), hi.z().max(p.z()));
                }
            }
            if !first {
                let _ = writeln!(
                    out,
                    "BOX {:.6} {:.6} {:.6} {:.6} {:.6} {:.6}",
                    lo.x(),
                    lo.y(),
                    lo.z(),
                    hi.x(),
                    hi.y(),
                    hi.z()
                );
            }
        }
        for i in 1..=structure.nb_nodes() as i32 {
            let n = structure.node(i);
            let _ = writeln!(
                out,
                "N {i} {:.17e} {:.17e} {}",
                n.location.x(),
                n.location.y(),
                n.state.index()
            );
        }
        for i in 1..=structure.nb_links() as i32 {
            let l = structure.link(i);
            let _ = writeln!(
                out,
                "L {i} {} {} {}",
                l.first_node(),
                l.last_node(),
                structure.link_movability(i).index()
            );
        }
        let mut vline = String::from("V");
        for &idx in indices {
            let _ = write!(vline, " {idx}");
        }
        let _ = writeln!(out, "{vline}");
        // TEMPORARY diagnostic (DELAUNDUMP): the surface nodes that
        // `AddVertices` will insert are not registered in the structure yet, so
        // they are appended as `N` lines with ids continuing after the
        // registered ones, together with the phase-1/phase-2 index lists.
        let base = structure.nb_nodes() as i32;
        for (k, (x, y)) in free_xy.iter().enumerate() {
            let _ = writeln!(out, "N {} {x:.17e} {y:.17e} {}", base + k as i32 + 1, VertexState::Free.index());
        }
        let mut v1 = String::from("V1");
        for &idx in indices {
            let _ = write!(v1, " {idx}");
        }
        let _ = writeln!(out, "{v1}");
        let mut v2 = String::from("V2");
        for k in 0..free_xy.len() {
            let _ = write!(v2, " {}", base + k as i32 + 1);
        }
        let _ = writeln!(out, "{v2}");
        let path = std::path::Path::new(dir).join(format!("delaun_in_f{face_index}.txt"));
        if std::fs::write(&path, out).is_err() {
            eprintln!("DELAUNDUMP: cannot write {}", path.display());
        }
    }

    /// TEMPORARY diagnostic (DELAUNDUMP_GRID): dumps the pre-`AddVertices`
    /// structure (nodes and links of the base mesh plus the registered surface
    /// nodes) together with the two insertion phases as `V1` / `V2`, so the
    /// probe can replay `BRepMesh_Delaun(structure, V1, cellsU, cellsV)` followed
    /// by `AddVertices(V2)` -- the exact control flow of
    /// `BRepMesh_DelaunayNodeInsertionMeshAlgo::postProcessMesh`.
    pub(super) fn dump_delaun_grid(
        structure: &DelaunDataStructure,
        dir: &str,
        face_index: usize,
        cells_u: i32,
        cells_v: i32,
        phase1: &[i32],
        phase2: &[i32],
    ) {
        use std::fmt::Write as _;
        let (tol_u, tol_v) = structure.tolerance();
        let mut out = String::new();
        let _ = writeln!(out, "TOL {tol_u:.17e} {tol_v:.17e}");
        let _ = writeln!(out, "CELLS {cells_u} {cells_v}");
        for i in 1..=structure.nb_nodes() as i32 {
            let n = structure.node(i);
            let _ = writeln!(
                out,
                "N {i} {:.17e} {:.17e} {}",
                n.location.x(),
                n.location.y(),
                n.state.index()
            );
        }
        for i in 1..=structure.nb_links() as i32 {
            let l = structure.link(i);
            let _ = writeln!(
                out,
                "L {i} {} {} {}",
                l.first_node(),
                l.last_node(),
                structure.link_movability(i).index()
            );
        }
        for (tag, ids) in [("V1", phase1), ("V2", phase2)] {
            let mut line = String::from(tag);
            for &idx in ids {
                let _ = write!(line, " {idx}");
            }
            let _ = writeln!(out, "{line}");
        }
        let path = std::path::Path::new(dir).join(format!("delaun_in_f{face_index}.txt"));
        if std::fs::write(&path, out).is_err() {
            eprintln!("DELAUNDUMP_GRID: cannot write {}", path.display());
        }
    }

    /// TEMPORARY diagnostic (DELAUNDUMP_GRID): per-pass state of the deflection
    /// control, node table (`ND <id> <x> <y>`) plus domain triangles
    /// (`T <tid> <a> <b> <c>`), for coordinate-level comparison with the probe's
    /// `OCCTDC P<n>` dump.
    pub(super) fn dump_domain_state(structure: &DelaunDataStructure, dir: &str, face_index: usize, tag: &str) {
        use std::fmt::Write as _;
        let mut out = String::new();
        for i in 1..=structure.nb_nodes() as i32 {
            let n = structure.node(i);
            if n.state == VertexState::Deleted {
                continue;
            }
            let _ = writeln!(out, "ND {i} {:.17e} {:.17e}", n.location.x(), n.location.y());
        }
        let mut tris: Vec<(i32, [i32; 3])> = Vec::new();
        for &tid in structure.elements_of_domain() {
            let e = structure.element(tid);
            let mut n = structure.element_nodes(&e);
            n.sort_unstable();
            tris.push((tid, n));
        }
        tris.sort_unstable();
        let _ = writeln!(out, "TRIS {}", tris.len());
        for (tid, n) in &tris {
            let _ = writeln!(out, "T {} {} {} {}", tid, n[0], n[1], n[2]);
        }
        let path = std::path::Path::new(dir).join(format!("delaun_{tag}_f{face_index}.txt"));
        if std::fs::write(&path, out).is_err() {
            eprintln!("DELAUNDUMP_GRID: cannot write {}", path.display());
        }
    }

    /// TEMPORARY diagnostic (DELONESTEP): appends one insertion step's domain
    /// triangle set (`STEP <k>` header, then `T <a> <b> <c>` with sorted node
    /// ids) to the step file, for step-by-step diffing against the probe's
    /// per-vertex `BRepMesh_Delaun::AddVertices` replay.
    pub(super) fn append_step_triangles(out: &mut String, step: usize, structure: &DelaunDataStructure) {
        use std::fmt::Write as _;
        let _ = writeln!(out, "STEP {step}");
        let mut tris: Vec<[i32; 3]> = Vec::new();
        for &tid in structure.elements_of_domain() {
            let e = structure.element(tid);
            let mut n = structure.element_nodes(&e);
            n.sort_unstable();
            tris.push(n);
        }
        tris.sort_unstable();
        for n in &tris {
            let _ = writeln!(out, "T {} {} {}", n[0], n[1], n[2]);
        }
    }

    /// TEMPORARY diagnostic (DELAUNDUMP_GRID): dumps the domain triangles of the
    /// structure as `T <a> <b> <c>` with sorted node ids, to be diffed against
    /// the probe's replay of the same structure.
    pub(super) fn dump_domain_triangles(structure: &DelaunDataStructure, dir: &str, face_index: usize, tag: &str) {
        use std::fmt::Write as _;
        let mut out = String::new();
        let mut tris: Vec<(i32, [i32; 3])> = Vec::new();
        for &tid in structure.elements_of_domain() {
            let e = structure.element(tid);
            let mut n = structure.element_nodes(&e);
            n.sort_unstable();
            tris.push((tid, n));
        }
        tris.sort_unstable();
        let _ = writeln!(out, "TRIS {}", tris.len());
        for (tid, n) in &tris {
            let _ = writeln!(out, "T {} {} {} {}", tid, n[0], n[1], n[2]);
        }
        let path = std::path::Path::new(dir).join(format!("delaun_{tag}_f{face_index}.txt"));
        if std::fs::write(&path, out).is_err() {
            eprintln!("DELAUNDUMP_GRID: cannot write {}", path.display());
        }
    }

    /// TEMPORARY diagnostic (DELDIAG): dumps the pre-triangulation structure of
    /// the face whose boundary contains the node at (13.0176, -24.7899, 290).
    pub(super) fn diag_dump_structure(&self, face_index: usize) {
        if std::env::var("DELDIAG").is_err() {
            return;
        }
        let nb = self.structure.nb_nodes() as i32;
        let is_target = (1..=nb).any(|i| {
            let p = self.structure.node(i).p3d;
            (p.x() - 13.017621944).abs() < 1e-3
                && (p.y() + 24.789947988).abs() < 1e-3
                && (p.z() - 290.0).abs() < 1e-3
        });
        if !is_target {
            return;
        }
        let nb_links = self.structure.nb_links() as i32;
        eprintln!("DELDIAG face={face_index} nodes={nb} links={nb_links}");
        for i in 1..=nb {
            let n = self.structure.node(i);
            eprintln!(
                "DELDIAG n i={i} st={} uv=({:.6},{:.6}) p=({:.6},{:.6},{:.6})",
                n.state.to_str(),
                n.location.x(),
                n.location.y(),
                n.p3d.x(),
                n.p3d.y(),
                n.p3d.z()
            );
        }
        for i in 1..=nb_links {
            let l = self.structure.link(i);
            let a = self.structure.node(l.first_node());
            let b = self.structure.node(l.last_node());
            eprintln!(
                "DELDIAG l i={i} st={} n=({},{}) p=[({:.6},{:.6},{:.6}),({:.6},{:.6},{:.6})]",
                self.structure.link_movability(i).to_str(),
                l.first_node(),
                l.last_node(),
                a.p3d.x(),
                a.p3d.y(),
                a.p3d.z(),
                b.p3d.x(),
                b.p3d.y(),
                b.p3d.z()
            );
        }
    }

}
