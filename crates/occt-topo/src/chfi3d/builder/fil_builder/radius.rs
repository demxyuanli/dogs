//! The radius / law API of `ChFi3d_FilBuilder`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder.cxx:234-398`
//! (OCCT 8.0.0). Every member starts with the same `if (IC <= NbElements())`
//! guard and forwards to the `ChFiDS_FilSpine` of the contour; the guard's
//! `else` arms answer `false` / `-1.`.
//!
//! Divergences that are not heuristics:
//! - `occ::handle<Law_Function>` is [`LawFunction`] and `Law_Composite` is
//!   [`LawComposite`]; `ChangeLaw(E)` returns a reference to the elementary law
//!   at the middle parameter of `E` (`ChFiDS_FilSpine.cxx:874-889`), which the
//!   port spells `change_fil_law` → `elementary_law[_mut](w)`,
//! - the two `Standard_DomainError`s `ChangeLaw` raises (`cxx:876-884`) become
//!   the `Err` of [`ChFi3dFilBuilder::get_bounds`] / [`get_law`] / [`set_law`],
//! - `double& F, double& L` becomes a returned pair, and the `bool` answer
//!   ("the law is not null") becomes `Option`.

use std::sync::{Arc, RwLock};

use occt_core::gp::GpXY;
use occt_geom::law::LawComposite;
use occt_geom::law::LawFunction;

use super::ChFi3dFilBuilder;
use crate::shape::{Edge, TopoShape, Vertex};

impl ChFi3dFilBuilder {
    /// `SetRadius(const handle(Law_Function)& C, const int IC, const int IinC)`
    /// (`cxx:234-241`).
    pub fn set_radius_law(
        &mut self,
        c: Arc<dyn LawFunction>,
        ic: usize,
        i_in_c: i32,
    ) -> Result<(), String> {
        if !self.ic_in_range(ic) {
            return Ok(());
        }
        let spine = self.fil_spine(ic)?;
        spine
            .write()
            .map_err(|_| "SetRadius: poisoned spine".to_string())?
            .set_radius_law(c, i_in_c);
        Ok(())
    }

    /// `IsConstant(const int IC)` (`cxx:245-253`).
    pub fn is_constant(&self, ic: usize) -> Result<bool, String> {
        if !self.ic_in_range(ic) {
            return Ok(false);
        }
        let spine = self.fil_spine(ic)?;
        let guard = spine
            .read()
            .map_err(|_| "IsConstant: poisoned spine".to_string())?;
        Ok(guard.is_fil_constant())
    }

    /// `Radius(const int IC)` (`cxx:257-265`): the constant radius of the whole
    /// contour, `-1.` when `IC` is out of range.
    ///
    /// `ChFiDS_FilSpine::Radius()` raises `Standard_DomainError` when the
    /// contour is not constant (`ChFiDS_FilSpine.cxx:358-364`); that is the
    /// port's `Err`.
    pub fn radius(&self, ic: usize) -> Result<f64, String> {
        if !self.ic_in_range(ic) {
            return Ok(-1.0);
        }
        let spine = self.fil_spine(ic)?;
        let guard = spine
            .read()
            .map_err(|_| "Radius: poisoned spine".to_string())?;
        guard
            .fil_radius()
            .ok_or_else(|| "Radius: ChFiDS_FilSpine::Radius raised Standard_DomainError".to_string())
    }

    /// `ResetContour(const int IC)` (`cxx:269-276`): `fsp->Reset(true)`.
    pub fn reset_contour(&mut self, ic: usize) -> Result<(), String> {
        if !self.ic_in_range(ic) {
            return Ok(());
        }
        let spine = self.fil_spine(ic)?;
        spine
            .write()
            .map_err(|_| "ResetContour: poisoned spine".to_string())?
            .reset_fil(true);
        Ok(())
    }

    /// `SetRadius(const double Radius, const int IC, const TopoDS_Edge& E)`
    /// (`cxx:280-287`).
    pub fn set_radius_on_edge(&mut self, radius: f64, ic: usize, e: &Edge) -> Result<(), String> {
        if !self.ic_in_range(ic) {
            return Ok(());
        }
        let spine = self.fil_spine(ic)?;
        spine
            .write()
            .map_err(|_| "SetRadius: poisoned spine".to_string())?
            .set_radius_on_edge(radius, &e.0);
        Ok(())
    }

    /// `UnSet(const int IC, const TopoDS_Edge& E)` (`cxx:291-298`):
    /// `fsp->UnSetRadius(E)`.
    pub fn unset_on_edge(&mut self, ic: usize, e: &Edge) -> Result<(), String> {
        if !self.ic_in_range(ic) {
            return Ok(());
        }
        let spine = self.fil_spine(ic)?;
        spine
            .write()
            .map_err(|_| "UnSet: poisoned spine".to_string())?
            .unset_radius_on_edge(&e.0);
        Ok(())
    }

    /// `SetRadius(const double Radius, const int IC, const TopoDS_Vertex& V)`
    /// (`cxx:302-309`).
    pub fn set_radius_on_vertex(
        &mut self,
        radius: f64,
        ic: usize,
        v: &Vertex,
    ) -> Result<(), String> {
        if !self.ic_in_range(ic) {
            return Ok(());
        }
        let spine = self.fil_spine(ic)?;
        spine
            .write()
            .map_err(|_| "SetRadius: poisoned spine".to_string())?
            .set_radius_on_vertex(radius, &v.0);
        Ok(())
    }

    /// `UnSet(const int IC, const TopoDS_Vertex& V)` (`cxx:313-320`):
    /// `fsp->UnSetRadius(V)`.
    pub fn unset_on_vertex(&mut self, ic: usize, v: &Vertex) -> Result<(), String> {
        if !self.ic_in_range(ic) {
            return Ok(());
        }
        let spine = self.fil_spine(ic)?;
        spine
            .write()
            .map_err(|_| "UnSet: poisoned spine".to_string())?
            .unset_radius_on_vertex(&v.0);
        Ok(())
    }

    /// `SetRadius(const gp_XY& UandR, const int IC, const int IinC)`
    /// (`cxx:324-331`).
    pub fn set_radius_at(&mut self, uand_r: GpXY, ic: usize, i_in_c: i32) -> Result<(), String> {
        if !self.ic_in_range(ic) {
            return Ok(());
        }
        let spine = self.fil_spine(ic)?;
        spine
            .write()
            .map_err(|_| "SetRadius: poisoned spine".to_string())?
            .set_radius_at(uand_r, i_in_c);
        Ok(())
    }

    /// `IsConstant(const int IC, const TopoDS_Edge& E)` (`cxx:335-343`):
    /// `fsp->IsConstant(fsp->Index(E))`.
    pub fn is_constant_of_edge(&self, ic: usize, e: &Edge) -> Result<bool, String> {
        if !self.ic_in_range(ic) {
            return Ok(false);
        }
        let spine = self.fil_spine(ic)?;
        let guard = spine
            .read()
            .map_err(|_| "IsConstant: poisoned spine".to_string())?;
        guard
            .is_fil_constant_on_edge(guard.index_of_edge(&e.0))
            .ok_or_else(|| "IsConstant: ChFiDS_FilSpine::IsConstant raised Standard_OutOfRange".to_string())
    }

    /// `Radius(const int IC, const TopoDS_Edge& E)` (`cxx:347-355`):
    /// `fsp->Radius(E)`, `-1.` when `IC` is out of range.
    pub fn radius_of_edge(&self, ic: usize, e: &Edge) -> Result<f64, String> {
        if !self.ic_in_range(ic) {
            return Ok(-1.0);
        }
        let spine = self.fil_spine(ic)?;
        let guard = spine
            .read()
            .map_err(|_| "Radius: poisoned spine".to_string())?;
        guard
            .fil_radius_of_edge(&e.0)
            .ok_or_else(|| "Radius: ChFiDS_FilSpine::Radius raised Standard_DomainError".to_string())
    }

    /// `GetBounds(const int IC, const TopoDS_Edge& E, double& F, double& L)`
    /// (`cxx:359-372`).
    ///
    /// `Ok(None)` is OCCT's `false` (the elementary law is null); the
    /// `Standard_DomainError` of `ChangeLaw` is the `Err`.
    pub fn get_bounds(&self, ic: usize, e: &Edge) -> Result<Option<(f64, f64)>, String> {
        let Some((law, w)) = self.change_law(ic, &e.0)? else {
            return Ok(None);
        };
        let guard = law.read().unwrap_or_else(|e| e.into_inner());
        Ok(guard.elementary_law(w).map(|l| l.bounds()))
    }

    /// `GetLaw(const int IC, const TopoDS_Edge& E)` (`cxx:376-384`): the
    /// elementary law serving `E`; OCCT's null handle is `None`.
    pub fn get_law(&self, ic: usize, e: &Edge) -> Result<Option<Arc<dyn LawFunction>>, String> {
        let Some((law, w)) = self.change_law(ic, &e.0)? else {
            return Ok(None);
        };
        let mut guard = law.write().unwrap_or_else(|e| e.into_inner());
        Ok(guard.elementary_law_mut(w).cloned())
    }

    /// `SetLaw(const int IC, const TopoDS_Edge& E, const handle(Law_Function)& L)`
    /// (`cxx:388-398`): `fsp->ChangeLaw(E) = L`.
    ///
    /// OCCT assigns through the reference returned by
    /// `Law_Composite::ChangeElementaryLaw(W)`; a composite without a law for
    /// `W` makes that `ChangeValue` raise `Standard_OutOfRange`, which the port
    /// answers as `Err`.
    pub fn set_law(
        &self,
        ic: usize,
        e: &Edge,
        l: Arc<dyn LawFunction>,
    ) -> Result<(), String> {
        let Some((law, w)) = self.change_law(ic, &e.0)? else {
            return Ok(());
        };
        let mut guard = law.write().unwrap_or_else(|e| e.into_inner());
        match guard.elementary_law_mut(w) {
            Some(slot) => {
                *slot = l;
                Ok(())
            }
            None => Err(
                "SetLaw: Law_Composite::ChangeElementaryLaw raised Standard_OutOfRange"
                    .to_string(),
            ),
        }
    }

    /// The `ChangeLaw(E)` prologue shared by [`Self::get_bounds`] /
    /// [`Self::get_law`] / [`Self::set_law`]:
    /// `ChFiDS_FilSpine::ChangeLaw` (`ChFiDS_FilSpine.cxx:874-889`).
    ///
    /// `Ok(None)` is the out-of-range `IC` guard of the three callers. Every
    /// `None` of `ChFiDS_FilSpine::ChangeLaw` itself is a
    /// `Standard_DomainError` in OCCT, so it is an `Err` here.
    fn change_law(
        &self,
        ic: usize,
        e: &TopoShape,
    ) -> Result<Option<(Arc<RwLock<LawComposite>>, f64)>, String> {
        if !self.ic_in_range(ic) {
            return Ok(None);
        }
        let spine = self.fil_spine(ic)?;
        let guard = spine
            .read()
            .map_err(|_| "ChangeLaw: poisoned spine".to_string())?;
        match guard.change_fil_law(e) {
            Some(v) => Ok(Some(v)),
            None => Err(
                "ChFiDS_FilSpine::ChangeLaw raised Standard_DomainError (limits not up to date, \
                 or no law on constant edges)"
                    .to_string(),
            ),
        }
    }
}
