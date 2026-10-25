//! `GeomFill_Sweep::BuildKPart` — recognition of the elementary surface a
//! (location law, section law) pair sweeps out.
//!
//! Source: `TKGeomAlgo/GeomFill/GeomFill_Sweep.cxx` (OCCT 8.0.0):
//! `IsSweepParallelSpine` `cxx:486-570`, `BuildKPart` `cxx:577-1167`
//! (branch (1) translation `cxx:594-798`, branch (2) rotation `cxx:802-1112`,
//! trimming `cxx:1114-1167`).
//!
//! Layout:
//! - [`state`] — the locals both branches share (`S`, `Ok`, `UFirst`, ...),
//! - [`translation`] — branch (1), [`rotation`] — branch (2).

mod rotation;
mod state;
mod translation;

use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::gp::GpVec;

pub(crate) use state::KPartState;

use super::Sweep;
use crate::geom_fill::location_law::LocationLaw;
use crate::geom_fill::section_law::SectionLaw;
use crate::rectangular_trimmed::GeomRectangularTrimmedSurface;
use crate::surface::Surface;

impl Sweep {
    /// `GeomFill_Sweep::BuildKPart` (`cxx:577-1167`).
    pub(crate) fn build_k_part(&mut self) -> bool {
        // `cxx:579-591`: the initialisers.
        let mut st = KPartState::new(
            self.tol3d.min(self.bound_tol),
            self.first,
            self.last,
            self.s_error,
        );

        // `(1)` `cxx:594-798`: rectilinear trajectory.
        let mut error = 0.0;
        if self.loc.is_translation(&mut error) {
            st.error = error;
            if let Some(result) = self.kpart_translation(&mut st) {
                self.write_back(&st);
                return result;
            }
        }

        // `(2)` `cxx:802-1112`: circular trajectory. `IsRotation` overwrites
        // `error` the way `cxx:804` does.
        if self.loc.is_rotation(&mut error) {
            st.error = error;
            if let Some(result) = self.kpart_rotation(&mut st) {
                self.write_back(&st);
                return result;
            }
        }

        // `cxx:1114-1167`: trim the recognised surface.
        let ok = st.ok;
        if ok && !self.trim_k_part(&mut st) {
            self.write_back(&st);
            return false;
        }
        self.write_back(&st);
        ok
    }

    /// `cxx:1114-1167`: swap U/V for the exchanged surfaces, then apply the
    /// `Geom_RectangularTrimmedSurface` window that matches the periodicity
    /// flags.
    ///
    /// `false` is the `Standard_ConstructionError` of the
    /// `Geom_RectangularTrimmedSurface` constructor (`cxx:1131-1150`), which in
    /// OCCT unwinds out of `BuildKPart`; the port reports "no KPart" so
    /// `Build` falls back to the global approximation.
    fn trim_k_part(&mut self, st: &mut KPartState) -> bool {
        if st.exch_uv {
            std::mem::swap(&mut st.is_u_periodic, &mut st.is_v_periodic);
            std::mem::swap(&mut st.u_first, &mut st.v_first);
            std::mem::swap(&mut st.u_last, &mut st.v_last);
        }
        let s = match st.s.clone() {
            Some(s) => s,
            None => return false,
        };

        let trimmed: Arc<dyn Surface> = if !st.is_u_periodic && !st.is_v_periodic {
            match GeomRectangularTrimmedSurface::uv(
                s,
                st.u_first,
                st.u_last,
                st.v_first,
                st.v_last,
            ) {
                Ok(t) => Arc::new(t),
                Err(_) => return false,
            }
        } else if st.is_u_periodic {
            if st.is_v_periodic {
                s
            } else {
                match GeomRectangularTrimmedSurface::one_param(
                    s,
                    st.v_first,
                    st.v_last,
                    false,
                ) {
                    Ok(t) => Arc::new(t),
                    Err(_) => return false,
                }
            }
        } else {
            match GeomRectangularTrimmedSurface::one_param(s, st.u_first, st.u_last, true) {
                Ok(t) => Arc::new(t),
                Err(_) => return false,
            }
        };

        self.surface = Some(trimmed);
        true
    }

    /// Copy the `KPartState` locals back onto the sweep.
    fn write_back(&mut self, st: &KPartState) {
        self.s_error = st.s_error;
        self.exch_uv = st.exch_uv;
        self.u_reversed = st.u_reversed;
    }
}

/// `IsSweepParallelSpine` (`cxx:486-570`, OCC5073): the spine has to translate
/// the section along its own first-to-last chord, i.e. the section normal and
/// the spine direction are parallel.
fn is_sweep_parallel_spine(
    loc: &Arc<dyn LocationLaw>,
    sec: &dyn SectionLaw,
    tol: f64,
) -> bool {
    let (a_first, a_last) = loc.get_domain();

    // `cxx:499-519`: the first transformation.
    let (m_begin, v_begin) = match loc.d0(a_first) {
        Ok(v) => v,
        Err(_) => return false,
    };
    let tf_begin = match state::trsf_from_matrix(&m_begin, &v_begin) {
        Ok(t) => t,
        Err(_) => return false,
    };

    // `cxx:522-539`: the last transformation.
    let (m_end, v_end) = match loc.d0(a_last) {
        Ok(v) => v,
        Err(_) => return false,
    };
    let tf_end = match state::trsf_from_matrix(&m_end, &v_end) {
        Ok(t) => t,
        Err(_) => return false,
    };

    // `cxx:541-546`: the section surface.
    let a_surf = match sec.bspline_surface() {
        Some(s) => s,
        None => return false,
    };
    let v_max = a_surf.v_range().1;

    // `cxx:548-556`: the first section, transformed.
    let first_section = match sec.constant_section() {
        Ok(c) => c,
        Err(_) => return false,
    };
    let u_first = state::occt_infinite(first_section.first_parameter());
    let mut line = match first_section.gp_line() {
        Some(l) => l,
        None => return false,
    };
    line.transform(&tf_begin);

    // `cxx:558-564`: the first point of the last V iso, transformed.
    let a_last_section = match a_surf.v_iso_curve(v_max) {
        Some(c) => c,
        None => return false,
    };
    let a_first_parameter = a_last_section.first_parameter();
    let mut a_pnt_last_sec = a_last_section.value(a_first_parameter);
    a_pnt_last_sec.transform(&tf_end);

    // `cxx:566-570`.
    let a_pnt_first_sec = clib::line_value(&line, u_first);
    let a_vec_sec = GpVec::from_pnts(&a_pnt_first_sec, &a_pnt_last_sec);
    let a_vec_spine = v_end.subtracted(&v_begin);
    a_vec_sec.is_parallel_ang(&a_vec_spine, tol)
}
