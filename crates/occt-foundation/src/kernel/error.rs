//! OCCT exception hierarchy → Rust error types.
//! Source: `Standard_Failure.hxx`, `Standard_*Error.hxx`
//!
//! OCCT rule: constructors/operations throw subclasses of `Standard_Failure`.
//! Rust rule: return `Result<T, OCCError>` or use `panic!` for unrecoverable bugs.

use std::fmt;

/// Root error type — maps to OCCT's `Standard_Failure`.
/// All higher-level errors are variants of this enum.
#[derive(Debug, Clone, PartialEq)]
pub enum OCCError {
    /// Invalid construction arguments (zero-norm Dir, negative radius, etc.)
    /// Maps: Standard_ConstructionError
    Construction(String),

    /// Value outside function domain
    /// Maps: Standard_DomainError
    Domain(String),

    /// Index out of valid range
    /// Maps: Standard_OutOfRange
    OutOfRange(String),

    /// Null handle dereference
    /// Maps: Standard_NullObject
    NullObject(String),

    /// Key not found in collection
    /// Maps: Standard_NoSuchObject
    NoSuchObject(String),

    /// Division by zero in algorithm
    /// Maps: Standard_DivideByZero
    DivideByZero(String),

    /// Numeric computation failure (overflow, singular matrix, etc.)
    /// Maps: Standard_NumericError
    Numeric(String),

    /// Dimension mismatch (vectors of different sizes, non-square matrix, etc.)
    /// Maps: Standard_DimensionError + Standard_DimensionMismatch
    Dimension(String),

    /// Generic failure with message
    /// Maps: Standard_Failure (base class)
    Failure(String),
}

impl fmt::Display for OCCError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OCCError::Construction(s) => write!(f, "Construction error: {s}"),
            OCCError::Domain(s) => write!(f, "Domain error: {s}"),
            OCCError::OutOfRange(s) => write!(f, "Out of range: {s}"),
            OCCError::NullObject(s) => write!(f, "Null object: {s}"),
            OCCError::NoSuchObject(s) => write!(f, "No such object: {s}"),
            OCCError::DivideByZero(s) => write!(f, "Divide by zero: {s}"),
            OCCError::Numeric(s) => write!(f, "Numeric error: {s}"),
            OCCError::Dimension(s) => write!(f, "Dimension error: {s}"),
            OCCError::Failure(s) => write!(f, "{s}"),
        }
    }
}

impl std::error::Error for OCCError {}

/// Convenience: return Err on condition (equivalent to `*_Raise_if` macros).
/// Usage: `raise_if(radius < 0.0, "radius must be >= 0")` → `Err(OCCError::Construction(...))`
#[macro_export]
macro_rules! raise_construction {
    ($cond:expr, $msg:expr) => {
        if $cond { return Err($crate::kernel::error::OCCError::Construction($msg.into())); }
    };
}

#[macro_export]
macro_rules! raise_domain {
    ($cond:expr, $msg:expr) => {
        if $cond { return Err($crate::kernel::error::OCCError::Domain($msg.into())); }
    };
}

#[macro_export]
macro_rules! raise_out_of_range {
    ($cond:expr, $msg:expr) => {
        if $cond { return Err($crate::kernel::error::OCCError::OutOfRange($msg.into())); }
    };
}

#[macro_export]
macro_rules! raise_dimension {
    ($cond:expr, $msg:expr) => {
        if $cond { return Err($crate::kernel::error::OCCError::Dimension($msg.into())); }
    };
}

/// Result type alias used throughout the OCCT port.
pub type OCCResult<T> = Result<T, OCCError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn construction_error() {
        let err = OCCError::Construction("bad input".into());
        assert_eq!(err.to_string(), "Construction error: bad input");
    }

    #[test]
    fn raise_macro() {
        fn bad_func(x: f64) -> OCCResult<f64> {
            raise_construction!(x < 0.0, "x must be >= 0");
            Ok(x * 2.0)
        }
        assert!(bad_func(-1.0).is_err());
        assert_eq!(bad_func(5.0).unwrap(), 10.0);
    }

    #[test]
    fn dimension_error() {
        let err = OCCError::Dimension("size mismatch".into());
        assert_eq!(err.to_string(), "Dimension error: size mismatch");
    }
}
