//! `TopOpeBRep_TypeLineCurve` — the kind of curve an intersection line is.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_TypeLineCurve.hxx:19-30`.

/// `enum TopOpeBRep_TypeLineCurve` (`hxx:19-30`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TypeLineCurve {
    /// `TopOpeBRep_ANALYTIC`
    Analytic,
    /// `TopOpeBRep_RESTRICTION`
    Restriction,
    /// `TopOpeBRep_WALKING`
    Walking,
    /// `TopOpeBRep_LINE`
    Line,
    /// `TopOpeBRep_CIRCLE`
    Circle,
    /// `TopOpeBRep_ELLIPSE`
    Ellipse,
    /// `TopOpeBRep_PARABOLA`
    Parabola,
    /// `TopOpeBRep_HYPERBOLA`
    Hyperbola,
    /// `TopOpeBRep_OTHERTYPE`
    #[default]
    OtherType,
}
