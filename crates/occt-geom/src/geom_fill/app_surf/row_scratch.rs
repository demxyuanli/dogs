//! The row buffers `AppBlend_AppSurf::InternalPerform` and
//! `Perform(Lin, SecGen, NbMaxP)` keep alive across sections.
//!
//! Source: `TKGeomAlgo/AppBlend/AppBlend_AppSurf.gxx:220-230` (`tabAppP`,
//! `tabAppV`, `tabP2d`, `tabV2d`, `tabW`, `tabDW`, `tabAppP2d`, `tabAppV2d`),
//! `:626-647` and the three `F.Section(...)` write sites (`:236-237`,
//! `:298-299`, `:336`).

use occt_core::gp::{GpPnt, GpPnt2d, GpVec, GpVec2d};

use crate::app_def::multi_point_constraint::MultiPointConstraint;
use crate::app_par_curves::constraint::Constraint;
use crate::geom_fill::app_surf::SectionValues;

/// `static bool scal = 1` (`gxx:36`): the row scaling switch.
pub(super) const SCAL: bool = true;

/// The row buffers `InternalPerform` keeps alive across sections: `tabAppP`,
/// `tabAppV`, `tabP2d`, `tabV2d`, `tabW`, `tabDW`, `tabAppP2d`, `tabAppV2d`
/// (`gxx:220-230`, `:626-647`).
pub(super) struct RowScratch {
    /// `NbUPoles`.
    pub(super) nb_upoles: usize,
    /// `NbPoles2d`.
    pub(super) nb_poles2d: usize,
    /// `tabAppP`.
    pub(super) tab_app_p: Vec<GpPnt>,
    /// `tabAppV`.
    pub(super) tab_app_v: Vec<GpVec>,
    /// `tabP2d`, sized `max(1, NbPoles2d)`.
    pub(super) tab_p2d: Vec<GpPnt2d>,
    /// `tabV2d`, sized `max(1, NbPoles2d)`.
    pub(super) tab_v2d: Vec<GpVec2d>,
    /// `tabW`.
    pub(super) tab_w: Vec<f64>,
    /// `tabDW`.
    pub(super) tab_dw: Vec<f64>,
    /// `tabAppP2d`, sized `NbPoles2d + NbUPoles` (points then weights).
    pub(super) tab_app_p2d: Vec<GpPnt2d>,
    /// `tabAppV2d`.
    pub(super) tab_app_v2d: Vec<GpVec2d>,
}

impl RowScratch {
    /// The array declarations of `gxx:220-230`.
    pub(super) fn new(nb_upoles: usize, nb_poles2d: usize) -> Self {
        let n2 = nb_poles2d.max(1);
        Self {
            nb_upoles,
            nb_poles2d,
            tab_app_p: vec![GpPnt::new(0.0, 0.0, 0.0); nb_upoles],
            tab_app_v: vec![GpVec::zero(); nb_upoles],
            tab_p2d: vec![GpPnt2d::new(0.0, 0.0); n2],
            tab_v2d: vec![GpVec2d::zero(); n2],
            tab_w: vec![0.0; nb_upoles],
            tab_dw: vec![0.0; nb_upoles],
            tab_app_p2d: vec![GpPnt2d::new(0.0, 0.0); nb_poles2d + nb_upoles],
            tab_app_v2d: vec![GpVec2d::zero(); nb_poles2d + nb_upoles],
        }
    }

    /// `F.Section(...)` writing into the fixed-size arrays
    /// (`gxx:236-237`, `:298-299`, `:336`).
    ///
    /// OCCT's `NCollection_Array1` out-parameters are sized by `GetShape`, so a
    /// `Section` returning a different count is an out-of-range write there and
    /// an error here.
    pub(super) fn load(&mut self, values: &SectionValues, with_deriv: bool) -> Result<(), String> {
        if values.poles.len() != self.nb_upoles || values.weights.len() != self.nb_upoles {
            return Err("Standard_OutOfRange: AppBlend_AppSurf::Section".into());
        }
        self.tab_app_p.copy_from_slice(&values.poles);
        self.tab_w.copy_from_slice(&values.weights);
        if self.nb_poles2d != 0 {
            if values.poles2d.len() != self.nb_poles2d {
                return Err("Standard_OutOfRange: AppBlend_AppSurf::Section".into());
            }
            self.tab_p2d[..self.nb_poles2d].copy_from_slice(&values.poles2d);
        }
        if with_deriv {
            if values.d_poles.len() != self.nb_upoles || values.d_weights.len() != self.nb_upoles {
                return Err("Standard_OutOfRange: AppBlend_AppSurf::Section".into());
            }
            self.tab_app_v.copy_from_slice(&values.d_poles);
            self.tab_dw.copy_from_slice(&values.d_weights);
            if self.nb_poles2d != 0 {
                if values.d_poles2d.len() != self.nb_poles2d {
                    return Err("Standard_OutOfRange: AppBlend_AppSurf::Section".into());
                }
                self.tab_v2d[..self.nb_poles2d].copy_from_slice(&values.d_poles2d);
            }
        }
        Ok(())
    }

    /// The pole / weight block of `gxx:243-263` (also `:301-322`, `:341-359`):
    /// multiply each pole by its weight and append the weight as a 2d column.
    pub(super) fn weight_poles(&mut self, with_deriv: bool) {
        for j in 1..=self.nb_poles2d {
            self.tab_app_p2d[j - 1] = self.tab_p2d[j - 1];
            if with_deriv {
                self.tab_app_v2d[j - 1] = self.tab_v2d[j - 1];
            }
        }
        for j in 1..=self.nb_upoles {
            let w = self.tab_w[j - 1];
            if with_deriv {
                let dw = self.tab_dw[j - 1];
                self.tab_app_v2d[self.nb_poles2d + j - 1] = GpVec2d::new(dw, 0.0);
                let mut new_dv = GpPnt::new(0.0, 0.0, 0.0).xyz().clone();
                // `newDv.SetLinearForm(tabDW(j), tabAppP(j).XYZ(), tabW(j),
                // tabAppV(j).XYZ())` (`gxx:257`).
                new_dv.set_linear_form(
                    dw,
                    self.tab_app_p[j - 1].xyz(),
                    w,
                    &self.tab_app_v[j - 1].coord,
                );
                self.tab_app_v[j - 1] = GpVec::from_xyz(&new_dv);
            }
            self.tab_app_p[j - 1] = GpPnt::from_xyz(&self.tab_app_p[j - 1].xyz().multiply(w));
            self.tab_app_p2d[self.nb_poles2d + j - 1] = GpPnt2d::new(w, 0.0);
        }
    }

    /// The `AppDef_MultiPointConstraint` construction of `gxx:265-271` (also
    /// `:323-329`, `:361-370`) and the matching row constraint.
    pub(super) fn multi_point(&self, with_deriv: bool) -> (MultiPointConstraint, Constraint) {
        if with_deriv {
            (
                MultiPointConstraint::from_3d_2d_tang(
                    &self.tab_app_p,
                    &self.tab_app_p2d,
                    &self.tab_app_v,
                    &self.tab_app_v2d,
                ),
                Constraint::TangencyPoint,
            )
        } else {
            (
                MultiPointConstraint::from_3d_2d(&self.tab_app_p, &self.tab_app_p2d),
                Constraint::PassPoint,
            )
        }
    }
}
