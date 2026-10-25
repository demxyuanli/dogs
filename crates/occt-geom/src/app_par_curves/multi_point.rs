//! `AppParCurves_MultiPoint` — the 3d / 2d pole bundle of one Bezier pole
//! index of a multi-curve approximation.
//!
//! Source: `ModelingData/TKGeomBase/AppParCurves/AppParCurves_MultiPoint.cxx:1-196`
//! and `AppParCurves_MultiPoint.hxx:37-136`, `.lxx:25-49` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `ttabPoint` / `ttabPoint2d` (`hxx:131-134`) are sparse OCCT handles that
//!   may be null when the matching count is zero (`cxx:36-48`); the port keeps
//!   two `Vec`s whose emptiness is the same condition, and `nbP` / `nbP2d`
//!   are the `Vec` lengths.
//! - `NCollection_Array1` with `Lower()` not 1 is copied into the 1-based
//!   storage by the `cxx:49-99` loops, which the port's slices already are.
//! - `Standard_OutOfRange` becomes a panic, OCCT's exception counterpart
//!   (same convention as `GeomFill_DraftTrihedron`).
//! - `Transform` (`cxx:105-122`) writes `tabPoint->SetValue(CuIndex, newP)`,
//!   i.e. it indexes the pole array with the **curve** index, not the pole
//!   index. That is the OCCT behaviour and it is kept verbatim, including the
//!   `OutOfRange` it raises when `CuIndex > nbP`.

use occt_core::gp::{GpPnt, GpPnt2d};

/// `AppParCurves_MultiPoint` (`AppParCurves_MultiPoint.hxx:38-136`).
#[derive(Debug, Clone)]
pub struct MultiPoint {
    /// `tabPoint` (`hxx:131`): the 3d poles, `nbP` of them.
    points: Vec<GpPnt>,
    /// `tabPoint2d` (`hxx:132`): the 2d poles, `nbP2d` of them.
    points2d: Vec<GpPnt2d>,
}

impl Default for MultiPoint {
    fn default() -> Self {
        Self::new()
    }
}

impl MultiPoint {
    /// `AppParCurves_MultiPoint()` (`cxx:28-30`): `nbP = nbP2d = 0`, both
    /// handles null.
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            points2d: Vec::new(),
        }
    }

    /// `AppParCurves_MultiPoint(NbPoles, NbPoles2d)` (`cxx:32-48`): the poles
    /// are default-constructed and filled later through `SetPoint` /
    /// `SetPoint2d`.
    pub fn with_size(nb_points: usize, nb_points2d: usize) -> Self {
        Self {
            points: vec![GpPnt::zero(); nb_points],
            points2d: vec![GpPnt2d::zero(); nb_points2d],
        }
    }

    /// `AppParCurves_MultiPoint(tabP)` (`cxx:49-60`): 3d poles only.
    pub fn from_3d(tab: &[GpPnt]) -> Self {
        Self {
            points: tab.to_vec(),
            points2d: Vec::new(),
        }
    }

    /// `AppParCurves_MultiPoint(tabP2d)` (`cxx:62-76`): 2d poles only.
    pub fn from_2d(tab: &[GpPnt2d]) -> Self {
        Self {
            points: Vec::new(),
            points2d: tab.to_vec(),
        }
    }

    /// `AppParCurves_MultiPoint(tabP, tabP2d)` (`cxx:78-100`).
    pub fn from_3d_2d(tab: &[GpPnt], tab2d: &[GpPnt2d]) -> Self {
        Self {
            points: tab.to_vec(),
            points2d: tab2d.to_vec(),
        }
    }

    /// `NbPoints` (`lxx:25-30`): `nbP`.
    pub fn nb_points(&self) -> usize {
        self.points.len()
    }

    /// `NbPoints2d` (`lxx:32-35`): `nbP2d`.
    pub fn nb_points2d(&self) -> usize {
        self.points2d.len()
    }

    /// `Dimension(Index)` (`lxx:37-48`): 3 for `Index <= nbP`, 2 above it.
    pub fn dimension(&self, index: i32) -> usize {
        let nb = (self.points.len() + self.points2d.len()) as i32;
        if index < 0 || index > nb {
            panic!("Standard_OutOfRange: AppParCurves_MultiPoint::Dimension");
        }
        if index <= self.points.len() as i32 {
            3
        } else {
            2
        }
    }

    /// `SetPoint(Index, Point)` (`cxx:140-146`).
    pub fn set_point(&mut self, index: i32, point: &GpPnt) {
        if index <= 0 || index > self.points.len() as i32 {
            panic!("Standard_OutOfRange: AppParCurves_MultiPoint::SetPoint");
        }
        self.points[(index - 1) as usize] = *point;
    }

    /// `Point(Index)` (`cxx:148-153`).
    pub fn point(&self, index: i32) -> &GpPnt {
        if index <= 0 || index > self.points.len() as i32 {
            panic!("Standard_OutOfRange: AppParCurves_MultiPoint::Point");
        }
        &self.points[(index - 1) as usize]
    }

    /// `SetPoint2d(Index, Point)` (`cxx:155-159`): `Index` counts 3d first.
    pub fn set_point2d(&mut self, index: i32, point: &GpPnt2d) {
        let nb_p = self.points.len() as i32;
        if index <= nb_p || index > nb_p + self.points2d.len() as i32 {
            panic!("Standard_OutOfRange: AppParCurves_MultiPoint::SetPoint2d");
        }
        self.points2d[(index - nb_p - 1) as usize] = *point;
    }

    /// `Point2d(Index)` (`cxx:161-166`).
    pub fn point2d(&self, index: i32) -> &GpPnt2d {
        let nb_p = self.points.len() as i32;
        if index <= nb_p || index > nb_p + self.points2d.len() as i32 {
            panic!("Standard_OutOfRange: AppParCurves_MultiPoint::Point2d");
        }
        &self.points2d[(index - nb_p - 1) as usize]
    }

    /// `Transform(CuIndex, x, dx, y, dy, z, dz)` (`cxx:105-122`):
    /// `newP = (x + P.X()*dx, y + P.Y()*dy, z + P.Z()*dz)`.
    ///
    /// `cxx:120` stores the result at the **curve** index, exactly as OCCT
    /// does; see the module notes.
    pub fn transform(&mut self, cu_index: i32, x: f64, dx: f64, y: f64, dy: f64, z: f64, dz: f64) {
        if self.dimension(cu_index) != 3 {
            panic!("Standard_OutOfRange: AppParCurves_MultiPoint::Transform");
        }
        let p = *self.point(cu_index);
        let new_p = GpPnt::new(x + p.x() * dx, y + p.y() * dy, z + p.z() * dz);
        if cu_index <= 0 || cu_index > self.points.len() as i32 {
            panic!("Standard_OutOfRange: AppParCurves_MultiPoint::Transform");
        }
        self.points[(cu_index - 1) as usize] = new_p;
    }

    /// `Transform2d(CuIndex, x, dx, y, dy)` (`cxx:124-138`).
    pub fn transform2d(&mut self, cu_index: i32, x: f64, dx: f64, y: f64, dy: f64) {
        if self.dimension(cu_index) != 2 {
            panic!("Standard_OutOfRange: AppParCurves_MultiPoint::Transform2d");
        }
        let p = *self.point2d(cu_index);
        let new_p = GpPnt2d::new(x + p.x() * dx, y + p.y() * dy);
        self.set_point2d(cu_index, &new_p);
    }

    /// `Dump(o)` (`cxx:168-196`): the same lines, one `String`.
    pub fn dump(&self) -> String {
        let mut out = String::new();
        out.push_str("AppParCurves_MultiPoint dump:\n");
        let nb_3d = self.nb_points();
        let nb_2d = self.nb_points2d();
        out.push_str(&format!(
            "It contains {} 3d points and {} 2d points.\n",
            nb_3d, nb_2d
        ));
        for (i, p) in self.points.iter().enumerate() {
            out.push_str(&format!("3D-Point #{}\n", i + 1));
            out.push_str(&format!(" Pole x = {}\n", p.x()));
            out.push_str(&format!(" Pole y = {}\n", p.y()));
            out.push_str(&format!(" Pole z = {}\n", p.z()));
        }
        for (i, p) in self.points2d.iter().enumerate() {
            out.push_str(&format!("2D-Point #{}\n", i + 1));
            out.push_str(&format!(" Pole x = {}\n", p.x()));
            out.push_str(&format!(" Pole y = {}\n", p.y()));
        }
        out
    }
}
