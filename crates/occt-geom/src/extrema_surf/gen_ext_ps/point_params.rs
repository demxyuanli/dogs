use super::super::prelude::*;

/// Precision::Infinite() (Precision.hxx:371) - the default of Extrema_POnSurf
/// (Extrema_POnSurf.hxx:33-36).
const PRECISION_INFINITE: f64 = occt_core::precision::INFINITE;

// ---------------------------------------------------------------------------
// Extrema_ElementType.hxx:20-26
// ---------------------------------------------------------------------------

/// Extrema_ElementType (Extrema_ElementType.hxx:20-26): which grid element a
/// [ExtremaPOnSurfParams] sits on.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum ExtremaElementType {
    #[default]
    Node,
    UIsoEdge,
    VIsoEdge,
    Face,
}

// ---------------------------------------------------------------------------
// Extrema_POnSurf.hxx:26-71
// ---------------------------------------------------------------------------

/// Extrema_POnSurf (Extrema_POnSurf.hxx:26-71): a 3d point and its (U, V)
/// parameters.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExtremaPOnSurf {
    u: f64,
    v: f64,
    p: GpPnt,
}

impl ExtremaPOnSurf {
    /// Extrema_POnSurf(theU, theV, theP) (Extrema_POnSurf.hxx:41-46).
    pub fn new(u: f64, v: f64, p: GpPnt) -> Self {
        ExtremaPOnSurf { u, v, p }
    }

    /// Value() (Extrema_POnSurf.hxx:49).
    pub fn value(&self) -> GpPnt {
        self.p
    }

    /// SetParameters (Extrema_POnSurf.hxx:53-58).
    pub fn set_parameters(&mut self, u: f64, v: f64, p: GpPnt) {
        self.u = u;
        self.v = v;
        self.p = p;
    }

    /// Parameter(U, V) (Extrema_POnSurf.hxx:61-65).
    pub fn parameter(&self) -> (f64, f64) {
        (self.u, self.v)
    }
}

impl Default for ExtremaPOnSurf {
    /// Extrema_POnSurf() (Extrema_POnSurf.hxx:32-37): (Infinite, Infinite) and
    /// P(Infinite, Infinite, Infinite).
    fn default() -> Self {
        ExtremaPOnSurf {
            u: PRECISION_INFINITE,
            v: PRECISION_INFINITE,
            p: GpPnt::new(PRECISION_INFINITE, PRECISION_INFINITE, PRECISION_INFINITE),
        }
    }
}

// ---------------------------------------------------------------------------
// Extrema_POnSurfParams.hxx:29-88
// ---------------------------------------------------------------------------

/// Extrema_POnSurfParams (Extrema_POnSurfParams.hxx:29-88): an [ExtremaPOnSurf]
/// plus the square distance to the projected point, the element type and the
/// grid indices it belongs to.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExtremaPOnSurfParams {
    point: ExtremaPOnSurf,
    sqr_distance: f64,
    element_type: ExtremaElementType,
    index_u: i32,
    index_v: i32,
}

impl ExtremaPOnSurfParams {
    /// Extrema_POnSurfParams() (Extrema_POnSurfParams.hxx:35-41).
    pub fn new_default() -> Self {
        ExtremaPOnSurfParams {
            point: ExtremaPOnSurf::default(),
            sqr_distance: 0.0,
            element_type: ExtremaElementType::Node,
            index_u: 0,
            index_v: 0,
        }
    }

    /// Extrema_POnSurfParams(theU, theV, thePnt)
    /// (Extrema_POnSurfParams.hxx:45-52).
    pub fn new(u: f64, v: f64, p: GpPnt) -> Self {
        ExtremaPOnSurfParams {
            point: ExtremaPOnSurf::new(u, v, p),
            sqr_distance: 0.0,
            element_type: ExtremaElementType::Node,
            index_u: 0,
            index_v: 0,
        }
    }

    /// SetSqrDistance (Extrema_POnSurfParams.hxx:56).
    pub fn set_sqr_distance(&mut self, d: f64) {
        self.sqr_distance = d;
    }

    /// GetSqrDistance (Extrema_POnSurfParams.hxx:59).
    pub fn sqr_distance(&self) -> f64 {
        self.sqr_distance
    }

    /// SetElementType (Extrema_POnSurfParams.hxx:62).
    pub fn set_element_type(&mut self, t: ExtremaElementType) {
        self.element_type = t;
    }

    /// GetElementType (Extrema_POnSurfParams.hxx:65).
    pub fn element_type(&self) -> ExtremaElementType {
        self.element_type
    }

    /// SetIndices (Extrema_POnSurfParams.hxx:69-73).
    pub fn set_indices(&mut self, iu: i32, iv: i32) {
        self.index_u = iu;
        self.index_v = iv;
    }

    /// GetIndices (Extrema_POnSurfParams.hxx:77-81).
    pub fn indices(&self) -> (i32, i32) {
        (self.index_u, self.index_v)
    }

    /// Extrema_POnSurf::Value() through the inheritance.
    pub fn value(&self) -> GpPnt {
        self.point.value()
    }

    /// Extrema_POnSurf::SetParameters() through the inheritance.
    pub fn set_parameters(&mut self, u: f64, v: f64, p: GpPnt) {
        self.point.set_parameters(u, v, p);
    }

    /// Extrema_POnSurf::Parameter() through the inheritance.
    pub fn parameter(&self) -> (f64, f64) {
        self.point.parameter()
    }
}

impl Default for ExtremaPOnSurfParams {
    fn default() -> Self {
        Self::new_default()
    }
}
