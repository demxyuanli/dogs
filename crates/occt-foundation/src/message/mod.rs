//! Progress reporting and messaging. Source: `Message/`
//! OCCT's Message_Messenger → Rust logging/tracing equivalent.
//! Ponteil note: use the `log` or `tracing` crate in production. This is a minimal shim.

/// Severity level for messages.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Gravity { Trace, Info, Warning, Alarm, Fail }

/// Progress indicator for long-running operations.
#[derive(Debug, Clone)]
pub struct ProgressIndicator {
    name: String,
    total: f64,
    current: f64,
    start: std::time::Instant,
}

impl ProgressIndicator {
    pub fn new(name: &str, total: f64) -> Self { Self { name: name.to_string(), total: total.max(1.0), current: 0.0, start: std::time::Instant::now() } }
    pub fn name(&self) -> &str { &self.name }
    pub fn total(&self) -> f64 { self.total }
    pub fn current(&self) -> f64 { self.current }
    pub fn fraction(&self) -> f64 { (self.current / self.total).min(1.0) }
    pub fn elapsed_secs(&self) -> f64 { self.start.elapsed().as_secs_f64() }
    pub fn advance(&mut self, step: f64) { self.current = (self.current + step).min(self.total); }
    pub fn reset(&mut self) { self.current = 0.0; self.start = std::time::Instant::now(); }
    pub fn is_done(&self) -> bool { self.current >= self.total }
}

/// Sends a message at a given gravity level (stub — use `log` crate in production).
pub fn send_message(gravity: Gravity, msg: &str) {
    match gravity {
        Gravity::Trace | Gravity::Info => eprintln!("[info] {msg}"),
        Gravity::Warning => eprintln!("[warn] {msg}"),
        Gravity::Alarm | Gravity::Fail => eprintln!("[error] {msg}"),
    }
}

/// Macro-based progress scoping.
#[macro_export]
macro_rules! progress_scope {
    ($name:expr, $total:expr, $body:block) => {{
        let mut _pi = $crate::message::ProgressIndicator::new($name, $total);
        $body
    }};
}
