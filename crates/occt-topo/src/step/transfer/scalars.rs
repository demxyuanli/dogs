use super::super::prelude::*;
use super::super::*;

pub(in crate::step) fn parse_ref(s: &str) -> Option<usize> {
    s.trim().strip_prefix('#')?.trim().parse().ok()
}

/// STEP logical: `.T.` is true, `.F.` is false. Missing / other → `default`.
pub(in crate::step) fn parse_logical(s: Option<&str>, default: bool) -> bool {
    match s.map(str::trim) {
        Some(".T.") => true,
        Some(".F.") => false,
        _ => default,
    }
}

pub(in crate::step) fn parse_ref_list(s: &str) -> Vec<usize> {
    let s = s.trim();
    if !(s.starts_with('(') && s.ends_with(')')) {
        return Vec::new();
    }
    split_top(&s[1..s.len() - 1])
        .into_iter()
        .filter_map(|a| parse_ref(&a))
        .collect()
}

pub(in crate::step) fn parse_f64(s: &str) -> Result<f64, String> {
    s.trim()
        .parse()
        .map_err(|_| format!("bad real literal '{s}'"))
}

pub(in crate::step) fn parse_xyz(s: &str) -> Result<GpXyz, String> {
    let s = s.trim();
    let inner = s.trim_start_matches('(').trim_end_matches(')');
    let parts: Vec<String> = split_top(inner)
        .into_iter()
        .map(|p| p.trim().to_string())
        .collect();
    if parts.len() != 3 {
        return Err(format!("expected 3-component tuple, got '{s}'"));
    }
    Ok(GpXyz::new(
        parse_f64(&parts[0])?,
        parse_f64(&parts[1])?,
        parse_f64(&parts[2])?,
    ))
}

/// Parse a 2D `(x, y)` tuple (STEP 2D curves use two-component CARTESIAN_POINTs).
pub(in crate::step) fn parse_xy(s: &str) -> Result<(f64, f64), String> {
    let s = s.trim();
    let inner = s.trim_start_matches('(').trim_end_matches(')');
    let parts: Vec<String> = split_top(inner)
        .into_iter()
        .map(|p| p.trim().to_string())
        .collect();
    if parts.len() != 2 {
        return Err(format!("expected 2-component tuple, got '{s}'"));
    }
    Ok((parse_f64(&parts[0])?, parse_f64(&parts[1])?))
}

pub(in crate::step) fn parse_str(s: &str) -> String {
    let s = s.trim();
    if s.len() >= 2 && s.starts_with('\'') && s.ends_with('\'') {
        // `StepData_StepReaderData::cleanText` (`StepData_StepReaderData.cxx:116-343`)
        // runs on every quoted string parameter before a reader sees it, in
        // place of the plain `''` collapse this function used to do.
        super::clean_text::clean_text(s)
    } else {
        s.to_string()
    }
}

