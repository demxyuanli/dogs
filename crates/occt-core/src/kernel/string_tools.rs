//! String utilities matching TCollection_AsciiString helpers.
//! OCCT's TCollection_AsciiString is essentially a String; these helpers
//! cover the common operations used across data exchange.
use std::collections::HashMap;

/// Split a string on a delimiter, trimming whitespace.
pub fn split_trim(s: &str, delim: char) -> Vec<String> {
    s.split(delim).map(|p| p.trim().to_string()).filter(|p| !p.is_empty()).collect()
}

/// Parse a "key = value" line (also handles "key : value").
/// Returns (key, value) trimmed.
pub fn parse_key_value(line: &str) -> Option<(String, String)> {
    let line = line.trim();
    if line.is_empty() || line.starts_with("//") || line.starts_with('#') { return None; }
    for sep in ['=', ':'] {
        if let Some(idx) = line.find(sep) {
            let k = line[..idx].trim();
            let v = line[idx + 1..].trim();
            if !k.is_empty() { return Some((k.to_string(), v.to_string())); }
        }
    }
    // Bare value: key = whole token if single token, else skip
    if line.contains(' ') { return None; }
    Some((line.to_string(), String::new()))
}

/// Read a properties-style file into a map (skips comments/blanks).
pub fn read_properties(content: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for line in content.lines() {
        if let Some((k, v)) = parse_key_value(line) { map.insert(k, v); }
    }
    map
}

/// OCCT RealToInt — clamp double to i32 range (matches Standard::RealToInt).
pub fn real_to_int(r: f64) -> i32 {
    if r >= i32::MAX as f64 { i32::MAX }
    else if r <= i32::MIN as f64 { i32::MIN }
    else { r as i32 }
}

/// Format a double for serialization with enough digits to round-trip.
pub fn fmt_real(r: f64) -> String {
    if r == r.trunc() && r.abs() < 1e15 {
        format!("{r:.1}")
    } else {
        format!("{r:.17}")
    }
}

/// Is the string a valid floating-point literal?
pub fn is_real_literal(s: &str) -> bool {
    s.parse::<f64>().is_ok()
}

/// Pad or truncate a string to width.
pub fn fixed_width(s: &str, width: usize) -> String {
    if s.len() >= width { s[..width].to_string() }
    else { format!("{}{}", s, " ".repeat(width - s.len())) }
}

/// Case-insensitive contains.
pub fn contains_ignore_case(haystack: &str, needle: &str) -> bool {
    haystack.to_lowercase().contains(&needle.to_lowercase())
}

/// Strip trailing zeros from a decimal representation (e.g. "1.500" → "1.5").
pub fn strip_trailing_zeros(s: &str) -> String {
    if !s.contains('.') { return s.to_string(); }
    let t = s.trim_end_matches('0');
    let t = t.trim_end_matches('.');
    t.to_string()
}

/// Convert bytes to a hex string.
pub fn to_hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// Parse hex string to bytes.
pub fn from_hex(s: &str) -> Option<Vec<u8>> {
    if s.len() % 2 != 0 { return None; }
    (0..s.len()).step_by(2).map(|i| u8::from_str_radix(&s[i..i+2], 16).ok()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn key_value_parsing() {
        assert_eq!(parse_key_value("KEY = value"), Some(("KEY".into(), "value".into())));
        assert_eq!(parse_key_value("A : 1.5"), Some(("A".into(), "1.5".into())));
        assert_eq!(parse_key_value("// comment"), None);
        assert_eq!(parse_key_value("  "), None);
    }

    #[test]
    fn props_file() {
        let content = "name = box\n// comment\nsize : 5.0\n";
        let m = read_properties(content);
        assert_eq!(m.get("name"), Some(&"box".to_string()));
        assert_eq!(m.get("size"), Some(&"5.0".to_string()));
        assert_eq!(m.len(), 2);
    }

    #[test]
    fn real_to_int_clamp() {
        assert_eq!(real_to_int(3.7), 3);
        assert_eq!(real_to_int(1e30), i32::MAX);
        assert_eq!(real_to_int(-1e30), i32::MIN);
    }

    #[test]
    fn strip_zeros() {
        assert_eq!(strip_trailing_zeros("1.500"), "1.5");
        assert_eq!(strip_trailing_zeros("1.000"), "1");
        assert_eq!(strip_trailing_zeros("2"), "2");
    }
}
