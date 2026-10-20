//! STEP argument list parsers and the `DATA` section extractor.

use super::super::prelude::*;
use super::*;

pub(in crate::step) fn parse_usize_list(s: &str) -> Vec<usize> {
    let s = s.trim();
    if !(s.starts_with('(') && s.ends_with(')')) {
        return Vec::new();
    }
    split_top(&s[1..s.len() - 1])
        .into_iter()
        .filter_map(|a| a.trim().parse().ok())
        .collect()
}

pub(in crate::step) fn parse_real_list(s: &str) -> Vec<f64> {
    let s = s.trim();
    if !(s.starts_with('(') && s.ends_with(')')) {
        return Vec::new();
    }
    split_top(&s[1..s.len() - 1])
        .into_iter()
        .filter_map(|a| a.trim().parse().ok())
        .collect()
}

/// Split a nested `((...),(...),...)` argument into its inner lists.
pub(in crate::step) fn parse_nested_lists(s: &str) -> Vec<Vec<String>> {
    let s = s.trim();
    if !(s.starts_with('(') && s.ends_with(')')) {
        return Vec::new();
    }
    split_top(&s[1..s.len() - 1])
        .into_iter()
        .filter_map(|row| {
            let row = row.trim();
            if row.starts_with('(') && row.ends_with(')') {
                Some(split_top(&row[1..row.len() - 1]))
            } else {
                None
            }
        })
        .collect()
}

/// Parse a nested list of entity references (e.g. a B-spline surface pole grid).
pub(in crate::step) fn parse_nested_ref_list(s: &str) -> Vec<Vec<usize>> {
    parse_nested_lists(s)
        .into_iter()
        .map(|row| row.into_iter().filter_map(|a| parse_ref(&a)).collect())
        .collect()
}

/// Parse a nested list of real values (e.g. a B-spline surface weight grid).
pub(in crate::step) fn parse_nested_real_list(s: &str) -> Vec<Vec<f64>> {
    parse_nested_lists(s)
        .into_iter()
        .map(|row| row.into_iter().filter_map(|a| a.trim().parse().ok()).collect())
        .collect()
}

/// Parse a STEP physical file, returning the model plus collected warnings.
/// Extract the `DATA` section of a physical file (between `DATA;` and the
/// closing `ENDSEC;`), validating the file terminator.
pub(in crate::step) fn data_section(content: &str) -> Result<&str, String> {
    if !content.contains("END-ISO-10303-21") {
        return Err("STEP file: missing END-ISO-10303-21 terminator".into());
    }
    let data_start = content
        .find("DATA;")
        .ok_or("STEP file: missing DATA section")?;
    let after_data = &content[data_start + 5..];
    let data_end = after_data
        .find("ENDSEC;")
        .ok_or("STEP file: DATA section not closed by ENDSEC")?;
    Ok(&after_data[..data_end])
}
