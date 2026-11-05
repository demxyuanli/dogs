//! `BRepFill_DraftLaw` — an `Edge3DLaw` whose frames are G0-stitched across the
//! path's vertices.
//!
//! Source: `BRepFill_DraftLaw.cxx:30-73`, `.hxx:21-33`.
//!
//! The class adds nothing but [`DraftLaw::clean_law`] on top of
//! `BRepFill_Edge3DLaw`: it walks the elementary laws and, wherever the two
//! normals at a path vertex are parallel, pre-multiplies the following law so
//! its frame starts exactly on the previous one.
//!
//! Translation notes:
//! - `occ::handle<GeomFill_LocationDraft>` (`cxx:39-43`) is stored as the
//!   generic `Arc<dyn GeomFillLocationLaw>` handle `BRepFill_Edge3DLaw`
//!   already consumes; the constructor body is only the inherited `Init` plus
//!   the per-edge `SetCurve`, so no `LocationDraft`-specific call is left.
//! - `ToG0` (`cxx:34-38`) is the same static helper as
//!   `BRepFill_LocationLaw.cxx:73-78`; the port shares
//!   [`crate::brep_fill_location_law::to_g0`] between the two.
//! - `T1` / `T2` (`cxx:62-63`) are filled and never read; the port keeps the
//!   two `Column(3)` reads in a comment only.

use std::sync::Arc;

use occt_core::gp::GpMat;
use occt_geom::geom_fill::LocationLaw as GeomFillLocationLaw;

use crate::brep_fill_edge3d_law::Edge3DLaw;
use crate::brep_fill_location_law::{to_g0, LocationLawBase};
use crate::shape::Wire;

/// `BRepFill_DraftLaw` (`BRepFill_DraftLaw.hxx:26-31`).
pub struct DraftLaw {
    /// The inherited `BRepFill_Edge3DLaw` state.
    inner: Edge3DLaw,
}

impl DraftLaw {
    /// `BRepFill_DraftLaw(Path, Law)` (`cxx:39-43`): the inherited
    /// `BRepFill_Edge3DLaw(Path, Law)` constructor and nothing else.
    pub fn new(path: &Wire, law: &Arc<dyn GeomFillLocationLaw>) -> Self {
        Self {
            inner: Edge3DLaw::new(path, law),
        }
    }

    /// `CleanLaw(TolAngular)` (`cxx:45-72`): give law `ipath` the rotation that
    /// takes its start frame onto the end frame of law `ipath - 1`, but only
    /// where the two normals are parallel within `tol_angular`.
    pub fn clean_law(&mut self, tol_angular: f64) -> Result<(), String> {
        let base = self.inner.base_mut();
        let nb = base.nb_law();
        if nb == 0 {
            // `myLaws->Value(1)` on a null array is out of range in OCCT.
            return Ok(());
        }
        // `myLaws->Value(1)->GetDomain(First, Last)` (`cxx:52`).
        let mut domain = base
            .law(1)
            .ok_or_else(|| NULL_LAW.to_string())?
            .get_domain();

        for ipath in 2..=nb {
            // `myLaws->Value(ipath - 1)->D0(Last, M1, V)` (`cxx:59`).
            let (m1, _v) = base
                .law(ipath - 1)
                .ok_or_else(|| NULL_LAW.to_string())?
                .d0(domain.1)?;
            // `myLaws->Value(ipath)->GetDomain(First, Last)` (`cxx:60`).
            let law = base.law(ipath).ok_or_else(|| NULL_LAW.to_string())?;
            domain = law.get_domain();
            // `myLaws->Value(ipath)->D0(First, M2, V)` (`cxx:61`).
            let (m2, _v) = law.d0(domain.0)?;

            // `T1.SetXYZ(M1.Column(3)); T2.SetXYZ(M2.Column(3));` (`cxx:62-63`):
            // both are written and never read, exactly as in the reference.
            let _t1 = m1.column(3);
            let _t2 = m2.column(3);
            let n1 = occt_core::gp::GpVec::from_xyz(&m1.column(1));
            let n2 = occt_core::gp::GpVec::from_xyz(&m2.column(1));

            // `if (N1.IsParallel(N2, TolAngular))` (`cxx:66`).
            if n1.is_parallel_ang(&n2, tol_angular) {
                // `ToG0(M1, M2, Trsf); myLaws->Value(ipath)->SetTrsf(Trsf)`
                // (`cxx:68-69`).
                let trsf: GpMat = to_g0(&m1, &m2)?;
                law.set_trsf(&trsf);
            }
        }
        Ok(())
    }

    /// The inherited `BRepFill_LocationLaw` state.
    pub fn base(&self) -> &LocationLawBase {
        self.inner.base()
    }

    /// The inherited `BRepFill_LocationLaw` state, mutably.
    pub fn base_mut(&mut self) -> &mut LocationLawBase {
        self.inner.base_mut()
    }
}

impl crate::brep_fill_location_law::LocationLaw for DraftLaw {
    fn base(&self) -> &LocationLawBase {
        self.inner.base()
    }

    fn base_mut(&mut self) -> &mut LocationLawBase {
        self.inner.base_mut()
    }
}

const NULL_LAW: &str = "Standard_NullObject: BRepFill_DraftLaw::myLaws";
