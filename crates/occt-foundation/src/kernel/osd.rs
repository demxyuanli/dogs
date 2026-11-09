//! OS layer abstraction over the Rust standard library.
//!
//! Maps the Open CASCADE `OSD` (Operating System Dependences) package onto
//! plain `std` functionality. No platform-specific code lives here; anything
//! that needs real platform hooks is a candidate for `cfg`-gated expansion.

/// Returns `true` if the file at `path` exists.
pub fn file_exists(path: &str) -> bool {
    std::path::Path::new(path).exists()
}

/// Reads the entire file at `path` into a string.
pub fn read_file(path: &str) -> std::io::Result<String> {
    std::fs::read_to_string(path)
}

/// Writes `contents` to the file at `path`, creating/truncating as needed.
pub fn write_file(path: &str, contents: &str) -> std::io::Result<()> {
    std::fs::write(path, contents)
}

/// Recursively creates every missing directory component of `path`.
pub fn create_dir_all(path: &str) -> std::io::Result<()> {
    std::fs::create_dir_all(path)
}

/// Deletes the file at `path`.
pub fn delete_file(path: &str) -> std::io::Result<()> {
    std::fs::remove_file(path)
}

/// Copies the file at `src` to `dst`.
pub fn copy_file(src: &str, dst: &str) -> std::io::Result<()> {
    std::fs::copy(src, dst).map(|_| ())
}

/// Returns the process's current working directory as a string.
pub fn current_dir() -> std::io::Result<String> {
    std::env::current_dir().map(|p| p.to_string_lossy().into_owned())
}

/// Returns the value of environment variable `name`, if set.
pub fn getenv(name: &str) -> Option<String> {
    std::env::var(name).ok()
}

/// Sets the environment variable `name` to `value`.
///
/// `std::env::set_var` became `unsafe` in edition 2024; the `unsafe` block
/// keeps this compiling on both edition 2021 and edition 2024, and the
/// `allow` silences the resulting `unused_unsafe` warning on 2021.
#[allow(unused_unsafe)]
pub fn setenv(name: &str, value: &str) {
    unsafe {
        std::env::set_var(name, value);
    }
}

/// Returns the machine host name (Windows `COMPUTERNAME`).
pub fn hostname() -> String {
    std::env::var("COMPUTERNAME").unwrap_or_else(|_| "unknown".into())
}

/// Returns the current user's login name (`USERNAME` or `USER`).
pub fn user_name() -> String {
    std::env::var("USERNAME")
        .or_else(|_| std::env::var("USER"))
        .unwrap_or_else(|_| "unknown".into())
}

/// Returns the operating system type, e.g. `"windows"` or `"linux"`.
pub fn os_type() -> String {
    std::env::consts::OS.into()
}

/// Returns the number of logical CPUs available to the process.
pub fn cpu_count() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

/// A simple wall-clock time wrapper, measured in seconds.
pub struct SystemTime {
    /// Seconds since the Unix epoch at construction time.
    pub seconds: f64,
}

impl SystemTime {
    /// Captures the current wall-clock time.
    pub fn now() -> Self {
        let d = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default();
        Self {
            seconds: d.as_secs_f64(),
        }
    }

    /// Seconds elapsed since this timestamp was created.
    pub fn elapsed(&self) -> f64 {
        let d = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default();
        d.as_secs_f64() - self.seconds
    }
}

/// Sleeps for `ms` milliseconds.
pub fn sleep_ms(ms: u64) {
    std::thread::sleep(std::time::Duration::from_millis(ms));
}

/// Computes a 64-bit FNV-1a hash of `data`, hex-encoded.
///
/// ponytail: placeholder named to stand in for `OSD`'s SHA-1 helper — this is
/// NOT cryptographic. Swap for a real SHA-1 implementation if any consumer
/// relies on actual digest values.
pub fn sha1_hex(data: &[u8]) -> String {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for b in data {
        hash ^= u64::from(*b);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    format!("{:016x}", hash)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn file_exists_is_false_for_missing() {
        assert!(!file_exists("definitely_not_a_real_file_xyz.tmp"));
    }

    #[test]
    fn cpu_count_is_positive() {
        assert!(cpu_count() > 0);
    }

    #[test]
    fn system_time_elapsed_is_nonnegative() {
        let t = SystemTime::now();
        assert!(t.elapsed() >= 0.0);
    }

    #[test]
    fn env_roundtrip() {
        setenv("OCCT_CORE_TEST_VAR", "hello");
        assert_eq!(getenv("OCCT_CORE_TEST_VAR").as_deref(), Some("hello"));
    }

    #[test]
    fn sha1_hex_is_stable_and_nonempty() {
        let a = sha1_hex(b"OpenCASCADE");
        let b = sha1_hex(b"OpenCASCADE");
        let c = sha1_hex(b"OpenCASCADe");
        assert_eq!(a, b);
        assert_ne!(a, c);
        assert_eq!(a.len(), 16);
    }
}
