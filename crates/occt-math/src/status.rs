//! Solver result codes. Source: `math_Status.hxx`

/// Return code for iterative solvers and optimizers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MathStatus {
    /// Solution found within tolerance
    Ok,
    /// Maximum iterations reached without convergence
    TooManyIterations,
    /// Function evaluation error (singular, out of domain, etc.)
    FunctionError,
    /// Direction search failed (line search cannot progress)
    DirectionSearchError,
    /// Root not bracketed (signs at endpoints are the same)
    NotBracketed,
}
