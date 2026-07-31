//! Resource file parser for the Open CASCADE `.res` format.
//!
//! A `.res` file is a flat list of `KEY : value` or `KEY value` lines.
//! Blank lines and `//` comments are ignored. Values are stored as raw
//! strings; typed accessors parse on demand.

use std::collections::HashMap;

/// A flat key/value store parsed from a `.res` resource file.
pub struct ResourceManager {
    entries: HashMap<String, String>,
}

impl ResourceManager {
    /// Creates an empty resource manager.
    pub fn new() -> Self {
        Self {
            entries: HashMap::new(),
        }
    }

    /// Parses the `.res` file at `path`, merging its entries into this store.
    pub fn load_from_file(&mut self, path: &str) -> std::io::Result<()> {
        let content = std::fs::read_to_string(path)?;
        self.load_from_str(&content);
        Ok(())
    }

    /// Parses `.res`-format text, merging its entries into this store.
    pub fn load_from_str(&mut self, content: &str) {
        for line in content.lines() {
            let line = line.trim();
            // Skip blank lines and `//` comments.
            if line.is_empty() || line.starts_with("//") {
                continue;
            }
            let (key, value) = if let Some(pos) = line.find(':') {
                (line[..pos].trim(), line[pos + 1..].trim())
            } else {
                let trimmed = line.trim_start();
                match trimmed.find(char::is_whitespace) {
                    Some(pos) => (&trimmed[..pos], trimmed[pos..].trim()),
                    None => continue, // single token, no value
                }
            };
            if !key.is_empty() {
                self.entries.insert(key.to_string(), value.to_string());
            }
        }
    }

    /// Returns the raw value stored under `key`.
    pub fn get(&self, key: &str) -> Option<&str> {
        self.entries.get(key).map(|s| s.as_str())
    }

    /// Parses the value under `key` as a boolean (`true`/`false`/`1`/`0`).
    pub fn get_bool(&self, key: &str) -> Option<bool> {
        match self.get(key)?.trim().to_ascii_lowercase().as_str() {
            "true" | "1" => Some(true),
            "false" | "0" => Some(false),
            _ => None,
        }
    }

    /// Parses the value under `key` as a 32-bit integer.
    pub fn get_int(&self, key: &str) -> Option<i32> {
        self.get(key)?.trim().parse().ok()
    }

    /// Parses the value under `key` as a double-precision float.
    pub fn get_double(&self, key: &str) -> Option<f64> {
        self.get(key)?.trim().parse().ok()
    }

    /// Inserts or overwrites the entry `key` with `value`.
    pub fn set(&mut self, key: &str, value: &str) {
        self.entries.insert(key.to_string(), value.to_string());
    }

    /// Returns `true` if `key` has an entry.
    pub fn contains(&self, key: &str) -> bool {
        self.entries.contains_key(key)
    }

    /// Returns all stored keys.
    pub fn keys(&self) -> Vec<&str> {
        self.entries.keys().map(|k| k.as_str()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_basic_resource_text() {
        let mut rm = ResourceManager::new();
        rm.load_from_str(
            "NAME : OpenCASCADE\n\
             FLOAT_VAL 1.5\n\
             BOOL true\n\
             FLAG 1\n\
             // a comment\n\
             \n",
        );

        assert_eq!(rm.get("NAME"), Some("OpenCASCADE"));
        assert_eq!(rm.get_double("FLOAT_VAL"), Some(1.5));
        assert_eq!(rm.get_bool("BOOL"), Some(true));
        assert_eq!(rm.get_bool("FLAG"), Some(true));
        assert!(rm.contains("NAME"));
        assert!(!rm.contains("comment"));

        let keys = rm.keys();
        assert_eq!(keys.len(), 4);
    }

    #[test]
    fn set_overwrites_and_get_int() {
        let mut rm = ResourceManager::new();
        rm.set("COUNT", "3");
        assert_eq!(rm.get_int("COUNT"), Some(3));
        rm.set("COUNT", "42");
        assert_eq!(rm.get_int("COUNT"), Some(42));
        assert_eq!(rm.get_int("MISSING"), None);
    }
}
