//! Options and alert (error/warning) reporting shared by the Boolean-component
//! algorithms.
//!
//! Port of `BOPAlgo_Options` (TKBO/BOPAlgo). Every algorithm in the Boolean
//! component derives its configuration from this class through
//! `BOPAlgo_Algo -> BOPAlgo_Options`. The ported surface keeps the four
//! per-run options:
//!
//! - **fuzzy tolerance** — an additional tolerance used to detect touching or
//!   coinciding cases during intersection. OCCT clamps it from below to
//!   `Precision::Confusion()` (1e-7);
//! - **parallel mode** — whether the algorithm may run on multiple threads.
//!   A process-global default (`BOPAlgo_Options::GetParallelMode`) is seeded
//!   at construction and can be overridden per instance;
//! - **OBB flag** — use Oriented Bounding Boxes to pre-filter the candidate
//!   intersection pairs;
//! - **alert report** — errors (fatal, the algorithm has failed) and warnings
//!   (non-fatal, degraded but usable result) collected while the algorithm
//!   runs. OCCT stores these in a `Message_Report` of typed alerts; here we
//!   keep two plain `Vec<String>` plus an `AlertKind` discriminator, which is
//!   all the ported callers need.

use std::sync::atomic::{AtomicBool, Ordering};

/// The process-global default for the parallel mode.
///
/// Mirrors OCCT's file-local `myGlobalRunParallel` in `BOPAlgo_Options.cxx`.
/// Relaxed ordering is enough: it is a hint flag, not a synchronisation
/// primitive.
static GLOBAL_PARALLEL_MODE: AtomicBool = AtomicBool::new(false);

/// Discriminator of an alert (report entry) produced by an algorithm.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlertKind {
    /// Fatal — the algorithm has failed and the result must not be used.
    Error,
    /// Non-fatal — the algorithm completed but with a degraded result.
    Warning,
}

impl AlertKind {
    /// Human-readable label, matching OCCT's `Message_Fail` / `Message_Warning`
    /// dump titles.
    pub fn as_str(self) -> &'static str {
        match self {
            AlertKind::Error => "Error",
            AlertKind::Warning => "Warning",
        }
    }
}

/// Per-run options shared by all Boolean-component algorithms.
///
/// Mirrors `BOPAlgo_Options`, with a plain-string alert report instead of
/// OCCT's typed `Message_Report`.
#[derive(Debug, Clone)]
pub struct BopAlgoOptions {
    /// Additional tolerance for touching/coinciding detection.
    fuzzy_value: f64,
    /// Run the computation on multiple threads when the algorithm supports it.
    run_parallel: bool,
    /// Use Oriented Bounding Boxes to pre-filter intersection pairs.
    use_obb: bool,
    /// Fatal alerts — a non-empty list means the algorithm has failed.
    errors: Vec<String>,
    /// Non-fatal alerts.
    warnings: Vec<String>,
}

impl Default for BopAlgoOptions {
    fn default() -> Self {
        Self::new()
    }
}

impl BopAlgoOptions {
    /// Empty constructor.
    ///
    /// `run_parallel` inherits the process-global parallel mode;
    /// `fuzzy_value` defaults to `Precision::Confusion()` (1e-7); OBB is off;
    /// the report is empty.
    pub fn new() -> Self {
        Self {
            fuzzy_value: 1e-7,
            run_parallel: get_parallel_mode(),
            use_obb: false,
            errors: Vec::new(),
            warnings: Vec::new(),
        }
    }

    /// Clears the collected errors and warnings.
    ///
    /// Exactly like `BOPAlgo_Options::Clear()` (which delegates to
    /// `myReport->Clear()`), the user-defined options (fuzzy/parallel/OBB) are
    /// *not* reset.
    pub fn clear(&mut self) {
        self.errors.clear();
        self.warnings.clear();
    }

    // -----------------------------------------------------------------------
    // Fuzzy tolerance
    // -----------------------------------------------------------------------

    /// Sets the additional tolerance.
    ///
    /// OCCT clamps the value from below with `Precision::Confusion()`:
    /// `myFuzzyValue = std::max(theFuzz, Precision::Confusion())`. A NaN is
    /// treated as the default confusion value.
    pub fn set_fuzzy_value(&mut self, v: f64) {
        self.fuzzy_value = if v.is_nan() { 1e-7 } else { v.max(1e-7) };
    }

    /// Returns the additional tolerance.
    pub fn fuzzy_value(&self) -> f64 {
        self.fuzzy_value
    }

    // -----------------------------------------------------------------------
    // Parallel processing mode
    // -----------------------------------------------------------------------

    /// Sets the per-instance parallel flag.
    pub fn set_run_parallel(&mut self, b: bool) {
        self.run_parallel = b;
    }

    /// Returns the per-instance parallel flag.
    pub fn run_parallel(&self) -> bool {
        self.run_parallel
    }

    // -----------------------------------------------------------------------
    // Oriented Bounding Boxes
    // -----------------------------------------------------------------------

    /// Enables/disables the use of Oriented Bounding Boxes for intersection
    /// pre-filtering.
    pub fn set_use_obb(&mut self, b: bool) {
        self.use_obb = b;
    }

    /// Returns the OBB flag.
    pub fn use_obb(&self) -> bool {
        self.use_obb
    }

    // -----------------------------------------------------------------------
    // Alert report
    // -----------------------------------------------------------------------

    /// Adds a fatal alert (the algorithm has failed).
    pub fn add_error(&mut self, msg: impl Into<String>) {
        self.errors.push(msg.into());
    }

    /// Adds a non-fatal alert.
    pub fn add_warning(&mut self, msg: impl Into<String>) {
        self.warnings.push(msg.into());
    }

    /// Adds an alert of the given kind.
    pub fn add_alert(&mut self, kind: AlertKind, msg: impl Into<String>) {
        match kind {
            AlertKind::Error => self.add_error(msg),
            AlertKind::Warning => self.add_warning(msg),
        }
    }

    /// Returns true if the algorithm has failed (at least one fatal alert).
    pub fn has_errors(&self) -> bool {
        !self.errors.is_empty()
    }

    /// Returns true if the algorithm has generated at least one warning.
    pub fn has_warnings(&self) -> bool {
        !self.warnings.is_empty()
    }

    /// Returns true if the algorithm has generated at least one alert of the
    /// given kind.
    pub fn has_alert(&self, kind: AlertKind) -> bool {
        match kind {
            AlertKind::Error => self.has_errors(),
            AlertKind::Warning => self.has_warnings(),
        }
    }

    /// Number of alerts of the given kind.
    pub fn alert_count(&self, kind: AlertKind) -> usize {
        match kind {
            AlertKind::Error => self.errors.len(),
            AlertKind::Warning => self.warnings.len(),
        }
    }

    /// Returns true when no fatal alert has been recorded.
    pub fn is_ok(&self) -> bool {
        !self.has_errors()
    }

    /// Returns the collected fatal alerts.
    pub fn errors(&self) -> &[String] {
        &self.errors
    }

    /// Returns the collected non-fatal alerts.
    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    /// Clears only the fatal alerts.
    pub fn clear_errors(&mut self) {
        self.errors.clear();
    }

    /// Clears only the non-fatal alerts.
    pub fn clear_warnings(&mut self) {
        self.warnings.clear();
    }

    /// Dumps the error status as a multi-line string (count + numbered list).
    ///
    /// Mirrors `BOPAlgo_Options::DumpErrors`, which forwards to
    /// `Message_Report::Dump(stream, Message_Fail)`.
    pub fn dump_errors(&self) -> String {
        dump_alerts("Error", &self.errors)
    }

    /// Dumps the warning statuses as a multi-line string.
    ///
    /// Mirrors `BOPAlgo_Options::DumpWarnings`.
    pub fn dump_warnings(&self) -> String {
        dump_alerts("Warning", &self.warnings)
    }
}

/// Formats a list of alerts like OCCT's `Message_Report::Dump`.
fn dump_alerts(kind: &str, alerts: &[String]) -> String {
    let mut out = String::new();
    out.push_str(&format!("{}s: {}\n", kind, alerts.len()));
    for (i, a) in alerts.iter().enumerate() {
        out.push_str(&format!("  {}. {}\n", i + 1, a));
    }
    out
}

// ---------------------------------------------------------------------------
// Process-global parallel mode
// ---------------------------------------------------------------------------

/// Returns the process-global parallel mode.
///
/// Mirrors `BOPAlgo_Options::GetParallelMode()`.
pub fn get_parallel_mode() -> bool {
    GLOBAL_PARALLEL_MODE.load(Ordering::Relaxed)
}

/// Sets the process-global parallel mode.
///
/// Mirrors `BOPAlgo_Options::SetParallelMode()`. Newly constructed
/// `BopAlgoOptions` instances inherit this value.
pub fn set_parallel_mode(v: bool) {
    GLOBAL_PARALLEL_MODE.store(v, Ordering::Relaxed);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults() {
        set_parallel_mode(false);
        let o = BopAlgoOptions::new();
        assert_eq!(o.fuzzy_value(), 1e-7, "fuzzy defaults to Precision::Confusion");
        assert_eq!(o.run_parallel(), get_parallel_mode(), "run_parallel inherits the global");
        assert!(!o.use_obb(), "OBB off by default");
        assert!(!o.has_errors());
        assert!(!o.has_warnings());
        assert!(o.is_ok());
    }

    #[test]
    fn fuzzy_value_roundtrip_and_clamp() {
        let mut o = BopAlgoOptions::new();
        o.set_fuzzy_value(0.5);
        assert_eq!(o.fuzzy_value(), 0.5);
        // OCCT clamps from below with Precision::Confusion().
        o.set_fuzzy_value(1e-9);
        assert_eq!(o.fuzzy_value(), 1e-7);
        o.set_fuzzy_value(0.0);
        assert_eq!(o.fuzzy_value(), 1e-7);
        // NaN is degraded to the default.
        o.set_fuzzy_value(f64::NAN);
        assert_eq!(o.fuzzy_value(), 1e-7);
    }

    #[test]
    fn parallel_and_obb_flags() {
        let mut o = BopAlgoOptions::new();
        assert!(!o.run_parallel());
        o.set_run_parallel(true);
        assert!(o.run_parallel());
        o.set_run_parallel(false);
        assert!(!o.run_parallel());

        assert!(!o.use_obb());
        o.set_use_obb(true);
        assert!(o.use_obb());
    }

    #[test]
    fn accumulate_errors_and_warnings() {
        let mut o = BopAlgoOptions::new();
        assert!(!o.has_errors());
        o.add_error("intersection failed");
        o.add_error("builder failed");
        o.add_warning("small edges ignored");
        assert!(o.has_errors());
        assert!(o.has_warnings());
        assert_eq!(o.errors().len(), 2);
        assert_eq!(o.warnings().len(), 1);
        assert_eq!(o.alert_count(AlertKind::Error), 2);
        assert_eq!(o.alert_count(AlertKind::Warning), 1);
        assert!(!o.is_ok());
        assert!(o.has_alert(AlertKind::Error));
        assert!(o.has_alert(AlertKind::Warning));
    }

    #[test]
    fn add_alert_by_kind() {
        let mut o = BopAlgoOptions::new();
        o.add_alert(AlertKind::Error, "E");
        o.add_alert(AlertKind::Warning, "W");
        assert_eq!(o.errors(), ["E"]);
        assert_eq!(o.warnings(), ["W"]);
        assert_eq!(AlertKind::Error.as_str(), "Error");
        assert_eq!(AlertKind::Warning.as_str(), "Warning");
    }

    #[test]
    fn clear_warnings_keeps_errors() {
        let mut o = BopAlgoOptions::new();
        o.add_error("E1");
        o.add_warning("W1");
        o.add_warning("W2");
        o.clear_warnings();
        assert!(!o.has_warnings());
        assert!(o.has_errors());
        assert_eq!(o.errors().len(), 1);
        // clear_errors leaves warnings alone
        o.clear_errors();
        assert!(!o.has_errors());
        assert!(o.is_ok());
    }

    #[test]
    fn clear_resets_report_not_options() {
        let mut o = BopAlgoOptions::new();
        o.set_fuzzy_value(0.25);
        o.set_run_parallel(true);
        o.set_use_obb(true);
        o.add_error("E");
        o.add_warning("W");
        o.clear();
        assert!(!o.has_errors());
        assert!(!o.has_warnings());
        // User-defined options survive Clear(), like BOPAlgo_Options::Clear().
        assert_eq!(o.fuzzy_value(), 0.25);
        assert!(o.run_parallel());
        assert!(o.use_obb());
    }

    #[test]
    fn dump_errors_and_warnings() {
        let mut o = BopAlgoOptions::new();
        o.add_error("E1");
        o.add_error("E2");
        let s = o.dump_errors();
        assert!(s.contains("Errors: 2"), "dump starts with count: {s}");
        assert!(s.contains("1. E1"));
        assert!(s.contains("2. E2"));

        o.add_warning("W1");
        let w = o.dump_warnings();
        assert!(w.contains("Warnings: 1"));
        assert!(w.contains("1. W1"));

        let empty = BopAlgoOptions::new();
        assert!(empty.dump_errors().contains("Errors: 0"));
    }

    #[test]
    fn global_parallel_mode_is_an_atomic() {
        set_parallel_mode(true);
        assert!(get_parallel_mode());
        // New instances inherit the global.
        let o = BopAlgoOptions::new();
        assert!(o.run_parallel());
        set_parallel_mode(false);
        assert!(!get_parallel_mode());
        let o2 = BopAlgoOptions::new();
        assert!(!o2.run_parallel());
    }
}
