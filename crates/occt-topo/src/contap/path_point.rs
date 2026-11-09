//! `Contap_ThePathPointOfTheSearch` - a solution point found on a restriction
//! arc by [`ContapTheSearch`](super::ContapTheSearch).
//!
//! Sources: `Contap_ThePathPointOfTheSearch.hxx:26-141`, `_0.cxx:23-57`.
//!
//! Deviation (handle -> value): `vtx` / `arc` are `occ::handle`s in OCCT; the
//! port stores the arc by value and the vertex in an `Option`. `Vertex()`
//! answers `None` where OCCT raises `Standard_DomainError` for a new point.

use std::sync::Arc;

use occt_core::gp::GpPnt;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_top_adaptor_hvertex::HVertex;
use crate::tgeometry::GeometryRegistry;

/// `Contap_ThePathPointOfTheSearch` (`hxx:26-76`).
#[derive(Clone)]
pub struct ContapThePathPointOfTheSearch {
    /// `point` (`hxx:70`).
    point: GpPnt,
    /// `tol` (`hxx:71`).
    tol: f64,
    /// `isnew` (`hxx:72`).
    is_new: bool,
    /// `vtx` (`hxx:73`): `None` iff `is_new`.
    vtx: Option<HVertex>,
    /// `arc` (`hxx:74`).
    arc: Arc<Curve2dAdaptor>,
    /// `param` (`hxx:75`).
    param: f64,
}

impl ContapThePathPointOfTheSearch {
    /// `Contap_ThePathPointOfTheSearch()` (`_0.cxx:23-29`): `tol = 0`,
    /// `isnew = true`, `param = 0`; `point` and `arc` are default-constructed.
    pub fn new() -> Self {
        Self {
            point: GpPnt::default(),
            tol: 0.0,
            is_new: true,
            vtx: None,
            arc: Arc::new(Curve2dAdaptor::new()),
            param: 0.0,
        }
    }

    /// `Contap_ThePathPointOfTheSearch(P, Tol, V, A, Parameter)`
    /// (`_0.cxx:31-44`): the point falls on the arc vertex `V`.
    pub fn on_vertex(
        p: GpPnt,
        tol: f64,
        v: HVertex,
        a: Arc<Curve2dAdaptor>,
        param: f64,
    ) -> Self {
        Self {
            point: p,
            tol,
            is_new: false,
            vtx: Some(v),
            arc: a,
            param,
        }
    }

    /// `Contap_ThePathPointOfTheSearch(P, Tol, A, Parameter)` (`_0.cxx:46-57`).
    pub fn at(p: GpPnt, tol: f64, a: Arc<Curve2dAdaptor>, param: f64) -> Self {
        Self {
            point: p,
            tol,
            is_new: true,
            vtx: None,
            arc: a,
            param,
        }
    }

    /// `SetValue(P, Tol, V, A, Parameter)` (`hxx:80-94`), inline in the class.
    pub fn set_value_on_vertex(
        &mut self,
        p: GpPnt,
        tol: f64,
        v: HVertex,
        a: Arc<Curve2dAdaptor>,
        param: f64,
    ) {
        self.is_new = false;
        self.point = p;
        self.tol = tol;
        self.vtx = Some(v);
        self.arc = a;
        self.param = param;
    }

    /// `SetValue(P, Tol, A, Parameter)` (`hxx:96-105`), inline in the class.
    pub fn set_value(&mut self, p: GpPnt, tol: f64, a: Arc<Curve2dAdaptor>, param: f64) {
        self.is_new = true;
        self.point = p;
        self.tol = tol;
        self.arc = a;
        self.param = param;
    }

    /// `Value()` (`hxx:107-110`).
    pub fn value(&self) -> &GpPnt {
        &self.point
    }

    /// `Tolerance()` (`hxx:112-115`).
    pub fn tolerance(&self) -> f64 {
        self.tol
    }

    /// `IsNew()` (`hxx:117-120`).
    pub fn is_new(&self) -> bool {
        self.is_new
    }

    /// `Vertex()` (`hxx:122-129`): `None` where OCCT raises
    /// `Standard_DomainError` for a new point.
    pub fn vertex(&self) -> Option<&HVertex> {
        self.vtx.as_ref()
    }

    /// `Arc()` (`hxx:131-134`).
    pub fn arc(&self) -> &Arc<Curve2dAdaptor> {
        &self.arc
    }

    /// `Parameter()` (`hxx:136-139`).
    pub fn parameter(&self) -> f64 {
        self.param
    }

    /// `ptsol.Arc() == A` (`IntStart_SearchOnBoundaries.gxx:937`, `:978`): OCCT
    /// compares the `Adaptor2d_Curve2d` handle addresses, and the domain hands
    /// out one arc per edge, so edge identity is the same test.
    pub fn is_on_arc(&self, a: &Curve2dAdaptor) -> bool {
        let mine = GeometryRegistry::shape_key(&self.arc.edge().0);
        let other = GeometryRegistry::shape_key(&a.edge().0);
        mine == other
    }
}

impl Default for ContapThePathPointOfTheSearch {
    fn default() -> Self {
        Self::new()
    }
}
