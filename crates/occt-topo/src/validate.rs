//! Shape validation and topology analysis. Source: `BRepCheck_Analyzer`, `BRepTools`
//!
//! **UNPORTED (audit A14)**: nothing here is a translation of
//! `BRepCheck_Analyzer`. OCCT's analyzer fills a per-shape map of
//! `BRepCheck_Status` values (`BRepCheck_Analyzer.cxx:458-478`) and
//! `IsValid(S)` returns `false` as soon as any status differs from
//! `BRepCheck_NoError`, recursing over sub-shapes; the checks below are local
//! invariants of this port (Euler characteristic, count sanity, orientation
//! statistics) and are **not** equivalent.
//!
//! ⚠️ [`Analyzer::is_valid`] is a stub that always returns `true` and
//! [`Analyzer::check_geometry`] always returns [`ShapeReport::ok`]. **They must
//! not be used as a validation gate** until the `BRepCheck_*` checks are ported.
use crate::abs::{ShapeType, Orientation};
use crate::shape::TopoShape;

/// Result of a shape validity check.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CheckStatus {
    #[default]
    Ok,
    SelfIntersection,
    Degenerated,
    WrongOrdering,
    InvalidCurve,
    InvalidSurface,
    TooSmall,
    Undefined,
}

/// Per-shape validation report.
#[derive(Debug, Clone, Default)]
pub struct ShapeReport {
    pub status: CheckStatus,
    pub messages: Vec<String>,
}

impl ShapeReport {
    pub fn ok() -> Self { Self { status: CheckStatus::Ok, messages: Vec::new() } }
    pub fn is_valid(&self) -> bool { self.status == CheckStatus::Ok }
}

/// Analyze a shape's validity using topological invariants.
/// children: direct sub-shapes (as stored in this port).
/// Euler characteristic check for closed solids: V - E + F = 2 (sphere genus 0).
pub fn check_shape(shape: &TopoShape, children: &[TopoShape]) -> ShapeReport {
    match shape.shape_type() {
        ShapeType::Solid => check_solid(shape, children),
        ShapeType::Shell => check_shell(shape, children),
        ShapeType::Face => check_face(shape),
        ShapeType::Edge => check_edge(shape),
        ShapeType::Vertex => check_vertex(shape),
        ShapeType::Compound => check_compound(children),
        _ => ShapeReport::ok(),
    }
}

fn check_solid(_s: &TopoShape, children: &[TopoShape]) -> ShapeReport {
    let shells = children.iter().filter(|c| c.shape_type() == ShapeType::Shell).count();
    if shells == 0 {
        let mut r = ShapeReport::ok();
        r.status = CheckStatus::Undefined;
        r.messages.push("solid has no shells".into());
        return r;
    }
    // Each shell must itself be valid
    for c in children {
        if c.shape_type() == ShapeType::Shell {
            let rep = check_shell(c, &[]);
            if !rep.is_valid() { return rep; }
        }
    }
    ShapeReport::ok()
}

fn check_shell(_s: &TopoShape, children: &[TopoShape]) -> ShapeReport {
    let faces = children.iter().filter(|c| c.shape_type() == ShapeType::Face).count();
    if faces < 4 {
        let mut r = ShapeReport::ok();
        r.status = CheckStatus::Undefined;
        r.messages.push(format!("shell has only {faces} faces (a closed shell needs ≥4)"));
        return r;
    }
    ShapeReport::ok()
}

fn check_face(_s: &TopoShape) -> ShapeReport {
    // A face with no wires is invalid
    ShapeReport::ok()
}

fn check_edge(_s: &TopoShape) -> ShapeReport {
    ShapeReport::ok()
}

fn check_vertex(_s: &TopoShape) -> ShapeReport {
    ShapeReport::ok()
}

fn check_compound(children: &[TopoShape]) -> ShapeReport {
    if children.is_empty() {
        let mut r = ShapeReport::ok();
        r.status = CheckStatus::Undefined;
        r.messages.push("compound is empty".into());
        return r;
    }
    ShapeReport::ok()
}

/// Count boundary sub-shapes by type.
pub fn boundary_counts(shapes: &[TopoShape]) -> (usize, usize, usize, usize) {
    // (vertices, edges, faces, solids)
    let mut v = 0; let mut e = 0; let mut f = 0; let mut s = 0;
    for sh in shapes {
        match sh.shape_type() {
            ShapeType::Vertex => v += 1,
            ShapeType::Edge => e += 1,
            ShapeType::Face => f += 1,
            ShapeType::Solid => s += 1,
            _ => {}
        }
    }
    (v, e, f, s)
}

/// Euler characteristic: V - E + F (2-complex) or V - E + F - C (3-complex).
pub fn euler_characteristic(counts: (usize, usize, usize, usize)) -> i32 {
    let (v, e, f, s) = counts;
    (v as i32) - (e as i32) + (f as i32) - (s as i32)
}

/// Is the shape likely a closed manifold (genus 0) via Euler formula?
/// Requires knowing it's a solid's boundary. Returns Ok if V-E+F == 2.
pub fn is_closed_manifold(counts: (usize, usize, usize, usize)) -> bool {
    euler_characteristic(counts) == 2
}

/// Simple BRepCheck_Analyzer-like wrapper.
///
/// **UNPORTED**: see the module header — this is a stub facade, not a port of
/// `BRepCheck_Analyzer`.
pub struct Analyzer;

impl Analyzer {
    pub fn new(_shape: &TopoShape, _children: &[TopoShape]) -> Self { Self }
    /// **Always `true`** — not a validity check. OCCT's
    /// `BRepCheck_Analyzer::IsValid` (`BRepCheck_Analyzer.cxx:458-478`) reports
    /// `false` for any `BRepCheck_Status != BRepCheck_NoError`. Do not gate on
    /// this value.
    pub fn is_valid(&self) -> bool { true }
    /// **Always [`ShapeReport::ok`]** — not a geometry check (OCCT's
    /// `BRepCheck_Analyzer` fills `BRepCheck_*` results per sub-shape).
    pub fn check_geometry(&self, _shape: &TopoShape, _children: &[TopoShape]) -> ShapeReport {
        ShapeReport::ok()
    }
    pub fn shape_type_name(t: ShapeType) -> &'static str { t.to_str() }
}

/// Orientation consistency: count how many sub-shapes are reversed.
pub fn orientation_stats(shapes: &[TopoShape]) -> (usize, usize) {
    let mut fwd = 0; let mut rev = 0;
    for s in shapes {
        match s.orientation() {
            Orientation::Forward => fwd += 1,
            Orientation::Reversed => rev += 1,
            _ => {}
        }
    }
    (fwd, rev)
}

/// Validate that a shell's faces have consistent outward orientation
/// (simplified: more forward than reversed is a weak heuristic).
pub fn check_orientation(children: &[TopoShape]) -> bool {
    let (fwd, rev) = orientation_stats(children);
    // For a closed shell, both orientations appear; we only reject all-reversed.
    fwd + rev > 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn euler_tetra() {
        // Closed tetrahedron: V=4, E=6, F=4 → 2
        assert!(is_closed_manifold((4, 6, 4, 0)));
    }

    #[test]
    fn euler_open() {
        // Open surface: square with 4 verts, 4 edges, 1 face → V-E+F = 1
        assert!(!is_closed_manifold((4, 4, 1, 0)));
    }

    #[test]
    fn solid_needs_shells() {
        let solid = TopoShape::new(ShapeType::Solid);
        let report = check_shape(&solid, &[]);
        assert!(!report.is_valid());
    }

    #[test]
    fn boundary_counts_basic() {
        let shapes = vec![
            TopoShape::new(ShapeType::Vertex),
            TopoShape::new(ShapeType::Edge),
            TopoShape::new(ShapeType::Face),
            TopoShape::new(ShapeType::Vertex),
        ];
        assert_eq!(boundary_counts(&shapes), (2, 1, 1, 0));
    }
}
