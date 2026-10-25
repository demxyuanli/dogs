//! `TopOpeBRepTool_REGUS::NearestF` and the `FUN_vectors` helper.
//!
//! Source (OCCT 8.0.0, `TopOpeBRepTool_REGUS.cxx`):
//! - `FUN_vectors(f, e, pare, nt, xx, tola, approx)` (`cxx:853-887`)
//! - `NearestF(e, lof, ffound)` (`cxx:891-1039`)
//!
//! `NearestF` decides which untouched face continues the current block: among
//! the faces of `lof` that *reduce* the 3d matter of the reference face `myf`,
//! it keeps the one with the smallest matter angle. When that angle falls
//! under `tola` the faces are tangent on `e` and the angle is re-sampled off
//! the edge by `MatterKPtg`.

use occt_core::gp::GpDir;

use crate::abs::Orientation;
use crate::shape::{Face, TopoShape};
use crate::topopebrep::tool::curvature::get_xx;
use crate::topopebrep::tool::matter::{matter_dirs, matter_k_ptg};
use crate::topopebrep::tool::ng_app::try_ng_app;
use crate::topopebrep::tool::project_curve::fun_tool_bounds;
use crate::topopebrep::tool::topology::{fun_tool_ng_geom_f, fun_tool_ng_geom_f_par, fun_tool_par_f};

use super::Regus;

/// `FUN_vectors(f, e, pare, nt, xx, tola, approx)` (`cxx:853-887`): the normal
/// `nt` and the in-face direction `xx` of `f` at `pare` on `e`.
fn fun_vectors(
    f: &TopoShape,
    e: &TopoShape,
    pare: f64,
    tola: f64,
    approx: bool,
) -> Option<(GpDir, GpDir)> {
    // `cxx:860-881`: `<nt>`.
    let mut nt = if approx {
        try_ng_app(pare, e, f, tola)?
    } else {
        fun_tool_ng_geom_f_par(pare, e, f)?
    };
    if f.orientation() == Orientation::Reversed {
        nt.reverse();
    }
    // `cxx:883-886`: `<xx>`.
    let xx = fun_tool_getxx(f, e, pare)?;
    Some((nt, xx))
}

/// `FUN_tool_getxx(Fi, Ei, parEi, XX)` (`TopOpeBRepTool_matter.cxx:132-147`):
/// project `parEi` on `Fi`, take the geometric normal there and orient the
/// result inside `Fi`.
pub(crate) fn fun_tool_getxx(f: &TopoShape, e: &TopoShape, par_e: f64) -> Option<GpDir> {
    let uv = fun_tool_par_f(e, par_e, f)?;
    let ng = fun_tool_ng_geom_f(&uv, &Face(f.clone()))?;
    get_xx(e, f, par_e, &ng)
}

impl Regus {
    /// `NearestF(e, lof, ffound)` (`cxx:891-1039`): the face of `lof` whose
    /// matter keeps the block continuous.
    ///
    /// `None` covers `FUN_vectors` failing on the reference or on a candidate
    /// (`cxx:926-931`, `:939-943`, `:997-1001`), `MatterKPtg` failing
    /// (`cxx:964-967`, `:1019-1022`) and "no face reduces the matter"
    /// (`cxx:976-980`).
    pub fn nearest_f(&self, e: &TopoShape, lof: &[TopoShape]) -> Option<TopoShape> {
        let fref = self.my_f.clone()?;

        // `cxx:912-918`: the reference parameter.
        let (f, l) = fun_tool_bounds(e);
        let eps = 0.45678;
        let pare = (1.0 - eps) * f + eps * l;

        // `cxx:919-922`: RONd (x, y, z) = (xxref, ntref, x^y).
        let tola = occt_core::precision::Precision::ANGULAR * 1.0e3;
        let (y, x) = fun_vectors(&fref, e, pare, tola, false)?;

        // `cxx:933-975`: the first candidate that reduces the matter decides
        // the starting angle.
        let mut ffound: Option<TopoShape> = None;
        let mut angfound = 0.0f64;
        let mut next = 0usize;
        for (i, fi) in lof.iter().enumerate() {
            let (nti, xxi) = fun_vectors(fi, e, pare, tola, false)?;
            let Some(mut angi) = matter_dirs(&x, &y, &xxi, &nti, tola) else {
                continue;
            };
            ffound = Some(fi.clone());

            if angi < tola {
                // `cxx:958-968`: the faces are tangent on `e`.
                angi = matter_k_ptg(&fref, fi, e)?;
            }
            angfound = angi;
            next = i + 1;
            break;
        }
        let ffound = ffound?;

        // `cxx:982-987`.
        if next >= lof.len() {
            return Some(ffound);
        }

        // `cxx:991-1037`: keep the smallest angle.
        let mut ffound = ffound;
        for fi in &lof[next..] {
            let (nti, xxi) = fun_vectors(fi, e, pare, tola, false)?;
            let Some(mut angi) = matter_dirs(&x, &y, &xxi, &nti, tola) else {
                continue;
            };
            if angi < tola {
                angi = matter_k_ptg(&fref, fi, e)?;
            }
            if angi > angfound {
                continue;
            }
            angfound = angi;
            ffound = fi.clone();
        }
        let _ = angfound;
        Some(ffound)
    }
}
