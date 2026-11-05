//! `BRepFill_ACRLaw` — an `BRepFill_LocationLaw` whose elementary laws are
//! `GeomFill_LocationGuide`s re-based on the path's reduced curvilinear
//! abscissae.
//!
//! Source: `BRepFill_ACRLaw.cxx:35-103`, `.hxx:24-37`.
//!
//! The constructor computes `BRepFill::ComputeACR(Path)` once, then gives path
//! edge `ipath` a copy of the caller's guide law whose *origin* interval is
//! `[ACR(ipath - 1), ACR(ipath)]` — that is what makes the guide follow the
//! auxiliary spine at the same reduced abscissa as the spine.
//!
//! Translation notes:
//! - `NCollection_HArray1<double> OrigParam(0, NbEdge)` is the `Vec` returned by
//!   [`crate::brep_fill::compute_acr`], which keeps OCCT's `0`-based lower bound.
//! - `GeomFill_LocationGuide::SetOrigine(t1, t2)` is
//!   [`occt_geom::geom_fill::LocationGuide::set_origine`]; `Loc->Copy()` is the
//!   `LocationLaw::copy` of the same trait.
//! - `GeomAdaptor_Curve(C, First, Last)` is [`GeomTrimmedCurve`], the trimmed
//!   window the rest of `occt-geom` uses.

use std::sync::Arc;

use occt_geom::curve::Curve;
use occt_geom::geom_fill::{LocationGuide, LocationLaw as GeomFillLocationLaw};
use occt_geom::trimmed::GeomTrimmedCurve;

use crate::abs::Orientation;
use crate::brep_fill::compute_acr;
use crate::brep_fill_location_law::LocationLawBase;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::shape::Wire;

/// `BRepFill_ACRLaw` (`BRepFill_ACRLaw.hxx:30-37`).
pub struct AcrLaw {
    /// The inherited `BRepFill_LocationLaw` state.
    base: LocationLawBase,
    /// `OrigParam` (`hxx:35`): the reduced abscissae, `0`-based like OCCT's
    /// `NCollection_HArray1<double>`.
    orig_param: Vec<f64>,
}

impl AcrLaw {
    /// `BRepFill_ACRLaw(Path, theLaw)` (`cxx:35-102`).
    ///
    /// `theLaw` stays the concrete `GeomFill_LocationGuide` because
    /// `SetOrigine` (`GeomFill_LocationGuide.cxx:1432-1436`) is not part of the
    /// `GeomFill_LocationLaw` base class.
    pub fn new(path: &Wire, law: &LocationGuide) -> Result<Self, String> {
        let mut base = LocationLawBase::new();
        base.init(path);

        // `cxx:42-49`: the number of path edges.
        let mut nb_edge = 0usize;
        let mut wexp = WireExplorer::new();
        wexp.init(path, None);
        while wexp.more() {
            nb_edge += 1;
            wexp.next();
        }

        // `cxx:51-55`: `OrigParam(0) = 0`, `Orig = BRepFill::ComputeACR(Path)`.
        let orig = compute_acr(path)?;
        let mut orig_param = vec![0.0; nb_edge + 1];
        // `OrigParam->SetValue(0, 0)` (`cxx:66`) stays as built, then
        // `OrigParam->SetValue(ipath, Orig(ipath))` (`cxx:68-71`).
        for ipath in 1..=nb_edge {
            orig_param[ipath] = orig.get(ipath).copied().unwrap_or(0.0);
        }

        // `cxx:73-102`.
        let mut ipath = 0usize;
        wexp.init(path, None);
        while wexp.more() {
            let Some(e) = wexp.current_edge() else {
                wexp.next();
                continue;
            };
            if !BRepTool::is_degenerated(&e) {
                ipath += 1;
                base.set_edge(ipath, &e);

                let Some(c) = BRepTool::edge_curve_world(&e) else {
                    wexp.next();
                    continue;
                };
                let (mut first, mut last) = BRepTool::edge_parameters(&e);
                let mut curve: Arc<dyn Curve> = c;
                if e.0.orientation() == Orientation::Reversed {
                    // `CBis = new Geom_TrimmedCurve(C, First, Last);
                    //  CBis->Reverse();` (`cxx:83-88`).
                    let mut cbis = GeomTrimmedCurve::new(Arc::clone(&curve), first, last);
                    cbis.reverse();
                    let cbis: Arc<dyn Curve> = Arc::new(cbis);
                    curve = cbis;
                    first = curve.first_parameter();
                    last = curve.last_parameter();
                }

                // `AC = new GeomAdaptor_Curve(C, First, Last)` (`cxx:90`).
                let ac: Arc<dyn Curve> =
                    Arc::new(GeomTrimmedCurve::new(Arc::clone(&curve), first, last));

                // `cxx:92-100`: `Loc = theLaw; Loc->SetOrigine(t1, t2);
                //  myLaws->SetValue(ipath, Loc->Copy()); ...
                //  ->SetCurve(AC)`.
                let t1 = orig_param[ipath - 1];
                let t2 = orig_param[ipath];
                law.set_origine(t1, t2);
                let copy = law.copy();
                let _ = copy.set_curve(ac);
                base.set_law(ipath, copy);
            }
            wexp.next();
        }

        Ok(Self { base, orig_param })
    }

    /// `OrigParam` (`hxx:35`): the reduced abscissae, `0`-based, entry `ipath`
    /// being the end of path edge `ipath`.
    pub fn orig_param(&self) -> &[f64] {
        &self.orig_param
    }

    /// The inherited `BRepFill_LocationLaw` state.
    pub fn base(&self) -> &LocationLawBase {
        &self.base
    }

    /// The inherited `BRepFill_LocationLaw` state, mutably.
    pub fn base_mut(&mut self) -> &mut LocationLawBase {
        &mut self.base
    }
}

impl crate::brep_fill_location_law::LocationLaw for AcrLaw {
    fn base(&self) -> &LocationLawBase {
        &self.base
    }

    fn base_mut(&mut self) -> &mut LocationLawBase {
        &mut self.base
    }
}
