//! Data carried by the `Contap_TheIWalking` walking engine:
//! `IntSurf_PathPoint` (`IntSurf_PathPoint.hxx:31-75`, `.lxx`, `.cxx:19-40`),
//! `IntSurf_InteriorPoint` (`IntSurf_InteriorPoint.hxx`, `.lxx`) and
//! `IntWalk_WalkingData` (`IntWalk_WalkingData.hxx`), plus the small numeric
//! helpers the `IntWalk_IWalking.gxx` bodies keep as locals.

use occt_core::gp::{GpDir2d, GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_math::MathVector;

use occt_topo::int_tools_wline::PntOn2S;

/// `IntSurf_PathPoint` (`IntSurf_PathPoint.hxx:31-75`): a point given at input
/// on which a walking line starts or stops.
#[derive(Clone)]
pub struct PathPoint {
    /// `pt` (`hxx:67`).
    pt: GpPnt,
    /// `ispass` (`hxx:68`).
    is_pass: bool,
    /// `istgt` (`hxx:69`).
    is_tgt: bool,
    /// `vectg` (`hxx:70`).
    vec_tg: GpVec,
    /// `dirtg` (`hxx:71`).
    dir_tg: GpDir2d,
    /// `sequv` (`hxx:72`); the first entry is the parametric position
    /// (`Value2d`), the following ones the other positions of a multiple point
    /// (`Parameters`).
    se_uv: Vec<(f64, f64)>,
}

impl PathPoint {
    /// `IntSurf_PathPoint()` (`cxx:23-27`): `ispass = istgt = true`, no
    /// parameter recorded.
    pub fn new() -> Self {
        Self {
            pt: GpPnt::default(),
            is_pass: true,
            is_tgt: true,
            vec_tg: GpVec::zero(),
            dir_tg: GpDir2d::new(1.0, 0.0).expect("+X"),
            se_uv: Vec::new(),
        }
    }

    /// `IntSurf_PathPoint(P, U, V)` (`cxx:29-37`).
    pub fn with_uv(p: GpPnt, u: f64, v: f64) -> Self {
        let mut r = Self::new();
        r.se_uv.push((u, v));
        r.pt = p;
        r
    }

    /// `SetValue(P, U, V)` (`cxx:39-44`): re-binds the point and resets the
    /// parameter list to the single `(U, V)`.
    pub fn set_value(&mut self, p: GpPnt, u: f64, v: f64) {
        self.pt = p;
        self.se_uv.clear();
        self.se_uv.push((u, v));
    }

    /// `AddUV(U, V)` (`lxx:19-22`).
    pub fn add_uv(&mut self, u: f64, v: f64) {
        self.se_uv.push((u, v));
    }

    /// `SetDirections(V, D)` (`lxx:24-29`): also clears the tangency flag.
    pub fn set_directions(&mut self, v: GpVec, d: GpDir2d) {
        self.is_tgt = false;
        self.vec_tg = v;
        self.dir_tg = d;
    }

    /// `SetTangency(Tang)` (`lxx:31-35`).
    pub fn set_tangency(&mut self, tang: bool) {
        self.is_tgt = tang;
    }

    /// `SetPassing(Pass)` (`lxx:37-41`).
    pub fn set_passing(&mut self, pass: bool) {
        self.is_pass = pass;
    }

    /// `Value()` (`lxx:43-46`).
    pub fn value(&self) -> &GpPnt {
        &self.pt
    }

    /// `Value2d(U, V)` (`lxx:48-53`): the first recorded parameter pair. An
    /// empty list answers `(0, 0)` where OCCT dereferences a null handle.
    pub fn value_2d(&self) -> (f64, f64) {
        self.se_uv.first().copied().unwrap_or((0.0, 0.0))
    }

    /// `IsPassingPnt()` (`lxx:55-58`).
    pub fn is_passing_pnt(&self) -> bool {
        self.is_pass
    }

    /// `IsTangent()` (`lxx:60-63`).
    pub fn is_tangent(&self) -> bool {
        self.is_tgt
    }

    /// `Direction3d()` (`lxx:65-72`); `None` where OCCT raises
    /// `StdFail_UndefinedDerivative`.
    pub fn direction_3d(&self) -> Option<&GpVec> {
        if self.is_tgt {
            None
        } else {
            Some(&self.vec_tg)
        }
    }

    /// `Direction2d()` (`lxx:74-81`).
    pub fn direction_2d(&self) -> Option<GpDir2d> {
        if self.is_tgt {
            None
        } else {
            Some(self.dir_tg)
        }
    }

    /// `Multiplicity()` (`lxx:83-87`): `sequv->Length() - 1`.
    pub fn multiplicity(&self) -> i32 {
        self.se_uv.len() as i32 - 1
    }

    /// `Parameters(Index, U, V)` (`lxx:89-94`): `sequv->Value(Index + 1)`, so
    /// `Index` is 1-based over the entries after the first one.
    pub fn parameters(&self, index: i32) -> Option<(f64, f64)> {
        self.se_uv.get(index.max(0) as usize).copied()
    }
}

impl Default for PathPoint {
    fn default() -> Self {
        Self::new()
    }
}

/// `IntSurf_InteriorPoint` (`IntSurf_InteriorPoint.hxx:35-102`): a point in the
/// interior of the parametric domain where a closed line starts, or that a
/// walking line passes through.
#[derive(Clone, Copy)]
pub struct InteriorPoint {
    /// `point`.
    pub point: GpPnt,
    /// `paramu`.
    pub paramu: f64,
    /// `paramv`.
    pub paramv: f64,
    /// `direc`.
    pub direc: GpVec,
    /// `direc2d`.
    pub direc2d: GpVec2d,
}

impl InteriorPoint {
    /// `IntSurf_InteriorPoint(P, U, V, Direc, Direc2d)`.
    pub fn new(p: GpPnt, u: f64, v: f64, direc: GpVec, direc2d: GpVec2d) -> Self {
        Self {
            point: p,
            paramu: u,
            paramv: v,
            direc,
            direc2d,
        }
    }

    /// `Value()`.
    pub fn value(&self) -> &GpPnt {
        &self.point
    }

    /// `Parameters(U, V)`.
    pub fn parameters(&self) -> (f64, f64) {
        (self.paramu, self.paramv)
    }

    /// `UParameter()`.
    pub fn u_parameter(&self) -> f64 {
        self.paramu
    }

    /// `VParameter()`.
    pub fn v_parameter(&self) -> f64 {
        self.paramv
    }

    /// `Direction()`.
    pub fn direction(&self) -> &GpVec {
        &self.direc
    }

    /// `Direction2d()`.
    pub fn direction_2d(&self) -> &GpVec2d {
        &self.direc2d
    }
}

/// `IntWalk_WalkingData` (`IntWalk_WalkingData.hxx`): the state of one start
/// point. `etat` encodes "tangent / passes" as described in
/// `IntWalk_IWalking.gxx:124-139`.
#[derive(Clone, Copy, Default)]
pub struct WalkingData {
    pub etat: i32,
    pub ustart: f64,
    pub vstart: f64,
}

/// `IntWalk_StatusDeflection` (`IntWalk_StatusDeflection.hxx`).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum StatusDeflection {
    /// `IntWalk_OK`.
    Ok,
    /// `IntWalk_PasTropGrand`.
    PasTropGrand,
    /// `IntWalk_ArretSurPointPrecedent`.
    ArretSurPointPrecedent,
    /// `IntWalk_ArretSurPoint`.
    ArretSurPoint,
    /// `IntWalk_PointConfondu`.
    PointConfondu,
}

/// `IntSurf_PntOn2S::ValueOnSurface(OnFirst)` (`IntSurf_PntOn2S.cxx:...`): the
/// `(u, v)` pair of one of the two surfaces. `reversed` is the walking flag:
/// when it is set the parameters of the first surface are the current ones.
pub(crate) fn uv_on(p: &PntOn2S, reversed: bool) -> (f64, f64) {
    if reversed {
        p.parameters_on_s1()
    } else {
        p.parameters_on_s2()
    }
}

/// `IntSurf_PntOn2S::ValueOnSurface(reversed)` as a 2D point.
pub(crate) fn pnt2d_on(p: &PntOn2S, reversed: bool) -> GpPnt2d {
    let (u, v) = uv_on(p, reversed);
    GpPnt2d::new(u, v)
}

/// `IntSurf_PntOn2S::SetValue(P, OnFirst, U, V)` parametric half: the 3D
/// point is set by the caller through [`PntOn2S::set_point`].
pub(crate) fn set_uv(p: &mut PntOn2S, on_first: bool, u: f64, v: f64) {
    p.set_uv_on(on_first, u, v);
}

/// `math_Vector(1, 2)` holding `(a, b)`.
pub(crate) fn math_vec2(a: f64, b: f64) -> MathVector {
    let mut v = MathVector::new(1, 2);
    v.set_value(1, a);
    v.set_value(2, b);
    v
}

/// `CutVectorByTolerances` (`gxx:391-397`).
pub(crate) fn cut_vector_by_tolerances(v: &mut GpVec2d, tolerance: [f64; 2]) {
    if v.x().abs() < tolerance[0] {
        v.set_coord(0.0, v.y());
    }
    if v.y().abs() < tolerance[1] {
        v.set_coord(v.x(), 0.0);
    }
}

/// `Bnd_Range` (`Bnd_Range.hxx:27-324`): the parts the walking engine uses.
/// A range is void when `myLast < myFirst`.
#[derive(Clone, Copy)]
pub(crate) struct BndRange {
    /// `myFirst` (`hxx:325`).
    first: f64,
    /// `myLast` (`hxx:326`).
    last: f64,
}

impl BndRange {
    /// `Bnd_Range()` (`hxx:43-47`): void, `(0, -1)`.
    pub(crate) fn void() -> Self {
        Self {
            first: 0.0,
            last: -1.0,
        }
    }

    /// `Bnd_Range(theMin, theMax)` (`hxx:49-56`).
    pub(crate) fn new(min: f64, max: f64) -> Self {
        assert!(max >= min, "Bnd_Range: Last < First");
        Self {
            first: min,
            last: max,
        }
    }

    /// `IsVoid()` (`hxx:217-218`).
    pub(crate) fn is_void(&self) -> bool {
        self.last < self.first
    }

    /// `Delta()` (`hxx:213-214`): `Max - Min`, negative for a void range.
    pub(crate) fn delta(&self) -> f64 {
        self.last - self.first
    }

    /// `Add(theParameter)` (`hxx:100-111`).
    pub(crate) fn add(&mut self, parameter: f64) {
        if self.is_void() {
            self.first = parameter;
            self.last = parameter;
            return;
        }
        self.first = self.first.min(parameter);
        self.last = self.last.max(parameter);
    }

    /// `Enlarge(theDelta)` (`hxx:224-234`).
    pub(crate) fn enlarge(&mut self, delta: f64) {
        if self.is_void() {
            return;
        }
        self.first -= delta;
        self.last += delta;
    }

    /// `Common(theOther)` (`hxx:58`, `cxx:...`).
    pub(crate) fn common(&mut self, other: &Self) {
        if other.is_void() {
            self.first = 0.0;
            self.last = -1.0;
            return;
        }
        if self.is_void() {
            return;
        }
        self.first = self.first.max(other.first);
        self.last = self.last.min(other.last);
    }

    // `SetVoid()` (`hxx:220-224`) has no caller in the ported walking bodies;
    // a void range is produced by `Bnd_Range::void()`.
}
