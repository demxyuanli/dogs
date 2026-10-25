//! `TopOpeBRepDS_Point` — a geometric point and a tolerance.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_Point.hxx:24-56`,
//! `TopOpeBRepDS_Point.cxx:20-82`.

use occt_core::gp::GpPnt;

/// Source: `TopOpeBRepDS_Point.hxx:49-53`
#[derive(Debug, Clone)]
pub struct Point {
    point: GpPnt,
    tolerance: f64,
    keep: bool,
}

impl Default for Point {
    /// `TopOpeBRepDS_Point()` (`cxx:22-26`).
    fn default() -> Self {
        Self::new()
    }
}

impl Point {
    /// `TopOpeBRepDS_Point()` (`cxx:22-26`): only `myKeep` is in the init list,
    /// so `gp_Pnt` is default-constructed (origin) and `myTolerance` is
    /// zero-initialized.
    pub fn new() -> Self {
        Self {
            point: GpPnt::new(0.0, 0.0, 0.0),
            tolerance: 0.0,
            keep: true,
        }
    }

    /// `TopOpeBRepDS_Point(const gp_Pnt& P, const double T)` (`cxx:30-35`).
    pub fn with_tolerance(p: GpPnt, t: f64) -> Self {
        Self {
            point: p,
            tolerance: t,
            keep: true,
        }
    }

    /// `TopOpeBRepDS_Point(const TopoDS_Shape& S)` (`cxx:39-43`):
    /// `Point(TopOpeBRepTool_ShapeTool::Pnt(S), ShapeTool::Tolerance(S))`.
    ///
    /// Both callees raise `Standard_ProgramError` when `S` is not a vertex,
    /// which is `Err` here.
    pub fn of_shape(s: &crate::shape::TopoShape) -> Result<Self, String> {
        let p = crate::topopebrep::tool::pnt_of(s)?;
        let t = crate::topopebrep::tool::tolerance_of(s)?;
        Ok(Self::with_tolerance(p, t))
    }

    /// `IsEqual(other)` (`cxx:47-52`): `myPoint.IsEqual(P.Point(),
    /// max(myTolerance, P.Tolerance()))`.
    pub fn is_equal(&self, other: &Point) -> bool {
        let t = self.tolerance.max(other.tolerance);
        self.point.distance(&other.point) <= t
    }

    /// `Point()` (`cxx:56-59`).
    pub fn point(&self) -> &GpPnt {
        &self.point
    }

    /// `ChangePoint()` (`cxx:63-66`).
    pub fn change_point(&mut self) -> &mut GpPnt {
        &mut self.point
    }

    /// `Tolerance()` (`cxx:70-73`).
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// `Tolerance(Tol)` (`cxx:77-80`).
    pub fn set_tolerance(&mut self, tol: f64) {
        self.tolerance = tol;
    }

    /// `Keep()` (`cxx:84-87`).
    pub fn keep(&self) -> bool {
        self.keep
    }

    /// `ChangeKeep(B)` (`cxx:90-93`).
    pub fn change_keep(&mut self, b: bool) {
        self.keep = b;
    }
}
