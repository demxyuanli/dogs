//! Single entry point for the diagnostic switches that live in the environment.
//!
//! The probes are compiled in unconditionally and some are consulted from hot
//! paths: every mesh node registration, every STEP entity translated.
//! `std::env::var` reaches the process environment on every call (on Windows
//! that is `GetEnvironmentVariableW`), so a probe check inside a per-edge or
//! per-node loop costs more than the code it guards. Boolean switches are
//! resolved once, on their first query, and reused afterwards.
//!
//! Diagnostics are enabled through the process environment before the run
//! starts; no library code sets these names at run time. Value-carrying probes
//! (a face index, an output path) use [`value`], which is not cached because
//! the value, not just its presence, is the input.

use std::sync::OnceLock;

/// Presence check for a value-carrying diagnostic (`None` when unset).
pub fn value(name: &str) -> Option<String> {
    std::env::var(name).ok()
}

/// Current face / stage tag for the `POSTDUMP` lines. Set by the mesh
/// pipeline (base pass, then each deflection-control insertion round) so the
/// interleaved per-vertex verdicts can be attributed to one stage.
pub static POST_STAGE: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

pub fn set_post_stage(v: usize) {
    POST_STAGE.store(v, std::sync::atomic::Ordering::Relaxed);
}

pub fn post_stage() -> usize {
    POST_STAGE.load(std::sync::atomic::Ordering::Relaxed)
}

/// Current face index for the `POSTDUMP` lines (see [`set_post_stage`]).
pub static POST_FACE: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

pub fn set_post_face(v: usize) {
    POST_FACE.store(v, std::sync::atomic::Ordering::Relaxed);
}

pub fn post_face() -> usize {
    POST_FACE.load(std::sync::atomic::Ordering::Relaxed)
}

macro_rules! flag {
    ($(#[$meta:meta])* $name:ident => $env:literal) => {
        $(#[$meta])*
        pub fn $name() -> bool {
            static ON: OnceLock<bool> = OnceLock::new();
            *ON.get_or_init(|| std::env::var($env).is_ok())
        }
    };
}

flag! {
    /// `DELREGSKEL` — mesh link registration skeleton dump
    /// (`NodeInsertionMeshAlgo::initDataStructure`).
    del_reg_skel => "DELREGSKEL"
}

flag! {
    /// `PREDCDUMP` — the triangulation the deflection control starts from
    /// (probe counterpart: `OCCTDC POST_INSERT`).
    pre_dc_dump => "PREDCDUMP"
}

flag! {
    /// `POSTDUMP` — one line per face with the final mesh size the export
    /// pipeline commits (`collectTriangles`: poly nodes / domain triangles).
    /// Counterpart of the probe's `FSTAT_OCC face= nodes= triangles=`.
    post_dump => "POSTDUMP"
}

flag! {
    /// `ZZSTEP_DBG` — STEP read node/record dump (`step/read_geometry.rs`).
    zzstep_dbg => "ZZSTEP_DBG"
}

flag! {
    /// `IGES_TRACE` — IGES read/write trace channel.
    iges_trace => "IGES_TRACE"
}

flag! {
    /// `IGES_TRACE_FACE` — IGES face-level trace channel.
    iges_trace_face => "IGES_TRACE_FACE"
}

flag! {
    /// `IGES_TRACE_WIRE` — IGES wire-level trace channel.
    iges_trace_wire => "IGES_TRACE_WIRE"
}

flag! {
    /// `IGES_TRACE_GROUP` — IGES group-level trace channel.
    iges_trace_group => "IGES_TRACE_GROUP"
}

flag! {
    /// `MBDIAG` — mesh builder diagnostics (`incremental_mesh/discret_root.rs`).
    mbdiag => "MBDIAG"
}

flag! {
    /// `DCDBG` — deflection-control diagnostics.
    dc_dbg => "DCDBG"
}

flag! {
    /// `DCDBGCTRL` — deflection-control post-process diagnostics.
    dc_dbg_ctrl => "DCDBGCTRL"
}

flag! {
    /// `SWEEP_DBG` — sweep diagnostics.
    sweep_dbg => "SWEEP_DBG"
}

flag! {
    /// `SHELLDBG` — shell check diagnostics.
    shell_dbg => "SHELLDBG"
}

flag! {
    /// `WDBGNEW` — wire building diagnostics.
    wire_dbg_new => "WDBGNEW"
}

flag! {
    /// `CONE_DBG` — cone fillet surface diagnostics.
    cone_dbg => "CONE_DBG"
}

flag! {
    /// `CLSWIRE` — 2D classifier wire registration dump.
    cls_wire => "CLSWIRE"
}

flag! {
    /// `CLSNODE` — 2D classifier node verdict dump.
    cls_node => "CLSNODE"
}

flag! {
    /// `BPUV` — boundary pcurve UV dump.
    bpuv => "BPUV"
}

flag! {
    /// `DELADDT` — Delaun vertex insertion dump.
    del_addt => "DELADDT"
}

flag! {
    /// `LENSDBG` — splitter edge-length dump.
    lens_dbg => "LENSDBG"
}

flag! {
    /// `SWEEP_FACE` — per-face sweep diagnostics.
    sweep_face => "SWEEP_FACE"
}

flag! {
    /// `PC2D` — pcurve discretization dump.
    pc2d => "PC2D"
}

flag! {
    /// `PCDEF` — pcurve deflection dump.
    pc_def => "PCDEF"
}

flag! {
    /// `PCW` — pcurve walk dump.
    pcw => "PCW"
}

flag! {
    /// `NURBSGRID` — NURBS splitter interval grid dump.
    nurbs_grid => "NURBSGRID"
}

flag! {
    /// `SPLITDBG` — splitter parameter interval dump.
    split_dbg => "SPLITDBG"
}

flag! {
    /// `PCSET` — pcurve registration dump.
    pc_set => "PCSET"
}

flag! {
    /// `PCS2DBG` — seam pcurve pair registration dump.
    pc_s2d => "PCS2DBG"
}
