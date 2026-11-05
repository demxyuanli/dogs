//! `BRepOffsetAPI_MakeFilling` — the modeling-API wrapper over the N-side
//! filling algorithm.
//!
//! Source: `ModelingAlgorithms/TKOffset/BRepOffsetAPI/BRepOffsetAPI_MakeFilling.cxx:22-212`
//! and `.hxx:66-240` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `(Degree, NbPtsOnCur, NbIter, Anisotropie, Tol2d, Tol3d, TolAng, TolCurv, MaxDeg, MaxSegments)` | `cxx:22-41` | [`BRepOffsetApiMakeFilling::with_params`] |
//! | the header's all-defaults constructor (`hxx:87-96`) | — | [`BRepOffsetApiMakeFilling::new`] |
//! | `SetConstrParam` / `SetResolParam` / `SetApproxParam` | `cxx:45-68` | the same names |
//! | `LoadInitSurface` | `cxx:72-75` | [`BRepOffsetApiMakeFilling::load_init_surface`] |
//! | the five `Add` overloads | `cxx:80-147` | `add*` |
//! | `Build(theRange)` | `cxx:152-156` | [`BRepOffsetApiMakeFilling::build`] |
//! | `IsDone` | `cxx:160-163` | [`BRepOffsetApiMakeFilling::is_done`] |
//! | `Generated(S)` | `cxx:170-173` | [`BRepOffsetApiMakeFilling::generated`] |
//! | `G0Error` / `G1Error` / `G2Error` (+ `Index`) | `cxx:179-212` | the same names |
//!
//! # Translation notes
//!
//! * Every member function of the reference is a one-line forward to
//!   `myFilling`, so the port keeps the same shape: [`Filling`] does all the
//!   work and this type only carries `myShape` — the `BRepBuilderAPI_MakeShape`
//!   member `Build` writes from `myFilling.Face()`.
//! * `IsDone()` answers `myFilling.IsDone()` (`cxx:162`), not the
//!   `BRepBuilderAPI_MakeShape` flag, so there is no separate done flag here.
//! * The header defaults `IsBound = true` on the two edge `Add` overloads
//!   (`hxx:170`, `hxx:188`). Rust has no default arguments, so each gets the
//!   two-argument form (`add`, `add_with_support`) plus the explicit-flag form
//!   (`add_with_flag`, `add_with_support_and_flag`).
//! * `Filling::Build` reports the reference's `Standard_Failure` throws as an
//!   `Err`, so `build` does too; `myShape` is written only when the filling
//!   completed (`cxx:155-156`).
//! * `Generated(S)` answers the filling's own `myGenerated` list
//!   (`cxx:172`), not a copy, so the port needs `&mut self` — the same borrow
//!   the [`Filling::generated`] method takes.

use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::Shape;

use crate::brep_fill_filling::Filling;
use crate::shape::{Edge, Face, TopoShape};

/// `BRepOffsetAPI_MakeFilling` (`BRepOffsetAPI_MakeFilling.hxx:66-240`).
pub struct BRepOffsetApiMakeFilling {
    /// `myFilling` (`hxx:238`).
    my_filling: Filling,
    /// `myShape` (the `BRepBuilderAPI_MakeShape` member `Build` writes).
    my_shape: Option<TopoShape>,
}

impl Default for BRepOffsetApiMakeFilling {
    /// The header's all-defaults constructor (`hxx:87-96`): `(3, 15, 2, false,
    /// 0.00001, 0.0001, 0.01, 0.1, 8, 9)`.
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetApiMakeFilling {
    /// `BRepOffsetAPI_MakeFilling()` (`hxx:87-96` defaults).
    pub fn new() -> Self {
        Self::with_params(3, 15, 2, false, 0.00001, 0.0001, 0.01, 0.1, 8, 9)
    }

    /// `BRepOffsetAPI_MakeFilling(Degree, NbPtsOnCur, NbIter, Anisotropie,
    /// Tol2d, Tol3d, TolAng, TolCurv, MaxDeg, MaxSegments)` (`cxx:22-41`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_params(
        degree: i32,
        nb_pts_on_cur: i32,
        nb_iter: i32,
        anisotropie: bool,
        tol2d: f64,
        tol3d: f64,
        tol_ang: f64,
        tol_curv: f64,
        max_deg: i32,
        max_segments: i32,
    ) -> Self {
        // `myFilling(...)` (`cxx:30-40`).
        Self {
            my_filling: Filling::new(
                degree,
                nb_pts_on_cur,
                nb_iter,
                anisotropie,
                tol2d,
                tol3d,
                tol_ang,
                tol_curv,
                max_deg,
                max_segments,
            ),
            my_shape: None,
        }
    }

    /// `myFilling` — the wrapped algorithm, for the queries this wrapper does
    /// not mirror (`Face()`, `IsDone()` on the raw filling, ...).
    pub fn filling(&self) -> &Filling {
        &self.my_filling
    }

    /// `myFilling`, mutable.
    pub fn filling_mut(&mut self) -> &mut Filling {
        &mut self.my_filling
    }

    /// `SetConstrParam(Tol2d, Tol3d, TolAng, TolCurv)` (`cxx:45-51`);
    /// the header's defaults are `hxx:107-109`.
    pub fn set_constr_param(&mut self, tol2d: f64, tol3d: f64, tol_ang: f64, tol_curv: f64) {
        self.my_filling
            .set_constr_param(tol2d, tol3d, tol_ang, tol_curv);
    }

    /// `SetResolParam(Degree, NbPtsOnCur, NbIter, Anisotropie)` (`cxx:56-62`);
    /// the header's defaults are `hxx:126-129`.
    pub fn set_resol_param(
        &mut self,
        degree: i32,
        nb_pts_on_cur: i32,
        nb_iter: i32,
        anisotropie: bool,
    ) {
        self.my_filling
            .set_resol_param(degree, nb_pts_on_cur, nb_iter, anisotropie);
    }

    /// `SetApproxParam(MaxDeg, MaxSegments)` (`cxx:67-70`); the header's
    /// defaults are `hxx:137`.
    pub fn set_approx_param(&mut self, max_deg: i32, max_segments: i32) {
        self.my_filling.set_approx_param(max_deg, max_segments);
    }

    /// `LoadInitSurface(Surf)` (`cxx:72-75`).
    pub fn load_init_surface(&mut self, surf: &Face) {
        self.my_filling.load_init_surface(surf);
    }

    /// `Add(Constr, Order, IsBound = true)` (`cxx:80-84`, `hxx:169-170`).
    pub fn add(&mut self, constr: &Edge, order: Shape) -> i32 {
        self.my_filling.add(constr, order, true)
    }

    /// `Add(Constr, Order, IsBound)` with the flag given explicitly.
    pub fn add_with_flag(&mut self, constr: &Edge, order: Shape, is_bound: bool) -> i32 {
        self.my_filling.add(constr, order, is_bound)
    }

    /// `Add(Constr, Support, Order, IsBound = true)` (`cxx:92-97`,
    /// `hxx:186-188`).
    pub fn add_with_support(&mut self, constr: &Edge, support: &Face, order: Shape) -> i32 {
        self.my_filling
            .add_with_support(constr, support, order, true)
    }

    /// `Add(Constr, Support, Order, IsBound)` with the flag given explicitly.
    pub fn add_with_support_and_flag(
        &mut self,
        constr: &Edge,
        support: &Face,
        order: Shape,
        is_bound: bool,
    ) -> i32 {
        self.my_filling
            .add_with_support(constr, support, order, is_bound)
    }

    /// `Add(Support, Order)` (`cxx:104-107`): the free constraint.
    pub fn add_free(&mut self, support: &Face, order: Shape) -> i32 {
        self.my_filling.add_free(support, order)
    }

    /// `Add(Point)` (`cxx:112-115`).
    pub fn add_point(&mut self, point: &GpPnt) -> Result<i32, String> {
        self.my_filling.add_point(point)
    }

    /// `Add(U, V, Support, Order)` (`cxx:123-132`).
    pub fn add_point_on_face(
        &mut self,
        u: f64,
        v: f64,
        support: &Face,
        order: Shape,
    ) -> Result<i32, String> {
        self.my_filling.add_point_on_face(u, v, support, order)
    }

    /// `Build(theRange)` (`cxx:152-156`): run the filling, then publish its
    /// face as `myShape`.
    pub fn build(&mut self) -> Result<(), String> {
        // `cxx:154`. The reference throws out of `Build` on a failure inside
        // `AddConstraints`, which leaves `myShape` untouched.
        self.my_filling.build()?;
        // `cxx:155`: `myShape = myFilling.Face()`.
        self.my_shape = self.my_filling.face().map(|f| f.0.clone());
        Ok(())
    }

    /// `IsDone()` (`cxx:160-163`): the filling's flag.
    pub fn is_done(&self) -> bool {
        self.my_filling.is_done()
    }

    /// `Shape()` — the `myShape` member. `None` before `Build` and when the
    /// filling produced no face.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `Generated(S)` (`cxx:170-173`): the filling's own list.
    pub fn generated(&mut self, s: &TopoShape) -> &[TopoShape] {
        self.my_filling.generated(s)
    }

    /// `G0Error()` (`cxx:179-182`).
    pub fn g0_error(&self) -> f64 {
        self.my_filling.g0_error()
    }

    /// `G1Error()` (`cxx:188-191`).
    pub fn g1_error(&self) -> f64 {
        self.my_filling.g1_error()
    }

    /// `G2Error()` (`cxx:197-200`).
    pub fn g2_error(&self) -> f64 {
        self.my_filling.g2_error()
    }

    /// `G0Error(Index)` (`cxx:210-213`): `Index` is the reference's 1-based
    /// constraint index.
    pub fn g0_error_of(&mut self, index: usize) -> Result<f64, String> {
        self.my_filling.g0_error_of(index)
    }

    /// `G1Error(Index)` (`cxx:219-222`), 1-based `Index`.
    pub fn g1_error_of(&mut self, index: usize) -> Result<f64, String> {
        self.my_filling.g1_error_of(index)
    }

    /// `G2Error(Index)` (`cxx:228-231`), 1-based `Index`.
    pub fn g2_error_of(&mut self, index: usize) -> Result<f64, String> {
        self.my_filling.g2_error_of(index)
    }
}
