use super::super::*;

impl BopBuilder {
    // -----------------------------------------------------------------------
    // Alert report
    // -----------------------------------------------------------------------

    /// True if the algorithm has failed (at least one fatal alert).
    pub fn has_errors(&self) -> bool {
        !self.errors.is_empty()
    }

    /// Returns the collected fatal alerts.
    pub fn errors(&self) -> &[String] {
        &self.errors
    }

    /// Returns the collected non-fatal alerts.
    pub fn warnings(&self) -> &[String] {
        &self.warnings
    }

    /// Adds a fatal alert.
    pub fn add_error(&mut self, msg: String) {
        self.errors.push(msg);
    }

    /// Adds a non-fatal alert.
    pub fn add_warning(&mut self, msg: String) {
        self.warnings.push(msg);
    }
}
