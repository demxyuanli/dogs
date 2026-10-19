//! Cached switches for the temporary diagnostics that live in the environment.
//!
//! The probes are compiled in unconditionally and are consulted from hot paths:
//! every edge range write (`GeometryRegistry::set_edge*`), every mesh node
//! registration, every STEP entity translated. `std::env::var` reaches the
//! process environment on every call (on Windows that is
//! `GetEnvironmentVariableW`), so a probe check inside a per-edge or per-node
//! loop costs more than the code it guards.
//!
//! Each switch is resolved once, on its first query, and reused afterwards. The
//! diagnostics are enabled through the process environment before the run
//! starts; no library code sets these names at run time.

use std::sync::OnceLock;

/// `RANGEW` — edge-range write probe in `GeometryRegistry::probe_range_write`.
pub fn range_w() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("RANGEW").is_ok())
}

/// `DELREGSKEL` — mesh link registration skeleton dump
/// (`NodeInsertionMeshAlgo::initDataStructure`).
pub fn del_reg_skel() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("DELREGSKEL").is_ok())
}

/// `ZZSTEP_DBG` — STEP read node/record dump (`step/read_geometry.rs`).
pub fn zzstep_dbg() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("ZZSTEP_DBG").is_ok())
}

/// `IGES_TRACE` — IGES read/write trace channel.
pub fn iges_trace() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("IGES_TRACE").is_ok())
}

/// `IGES_TRACE_FACE` — IGES face-level trace channel.
pub fn iges_trace_face() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("IGES_TRACE_FACE").is_ok())
}

/// `IGES_TRACE_WIRE` — IGES wire-level trace channel.
pub fn iges_trace_wire() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("IGES_TRACE_WIRE").is_ok())
}

/// `MBDIAG` — mesh builder diagnostics (`incremental_mesh/discret_root.rs`).
pub fn mbdiag() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("MBDIAG").is_ok())
}

/// `DCDBG` — deflection-control diagnostics.
pub fn dc_dbg() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("DCDBG").is_ok())
}

/// `SWEEP_DBG` — sweep diagnostics.
pub fn sweep_dbg() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("SWEEP_DBG").is_ok())
}

/// `SHELLDBG` — shell check diagnostics.
pub fn shell_dbg() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var("SHELLDBG").is_ok())
}
