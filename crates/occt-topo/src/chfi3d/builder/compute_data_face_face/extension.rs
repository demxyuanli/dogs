//! The extension half of the face / face `ComputeData`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_6.cxx:2336-2505`
//! (OCCT 8.0.0). This is the part [`super::ChFi3dBuilder::compute_data_face_face`]
//! runs once the walk has produced a presentable line: both ends are extended
//! through `SearchFace` + `IsObst`, the sides that ran into a face are
//! truncated with `SearchIndex` + `UpdateLine`, and the `Gd1` / `Gd2` / `Gf1` /
//! `Gf2` flags plus the two extension values of `SurfData` are written.
//!
//! Split out of `super` only for file length. The two `if reverse` / `if
//! forward` blocks are contiguous in the C++ and contiguous here; nothing of
//! the control flow moved across the split, and the state the block reads and
//! writes is the caller's locals bound by [`ExtensionState`].
//!
//! The one port-only local is `target`: OCCT reuses the `Target` variable of
//! the walk (`cxx:1756-1775`), which is dead after this point, so each branch
//! opens its own `let mut target`.

use crate::brep_blend::line::Line;
use crate::chfi_ds::SurfData;
use crate::chfi3d::builder_6_helpers::search_index;

use super::super::compute_data::point_parameter;
use super::super::update_line;
use super::super::ChFi3dBuilder;

/// The locals of [`super::ChFi3dBuilder::compute_data_face_face`] that the
/// extension half (`cxx:2336-2505`) reads and writes.
pub(super) struct ExtensionState<'a> {
    /// `Lin`.
    pub line: &'a mut Line,
    /// `Data`.
    pub data: &'a mut SurfData,
    /// `First`.
    pub first: &'a mut f64,
    /// `Last`.
    pub last: &'a mut f64,
    /// `SpFirst` as the walk left it.
    pub sp_first: f64,
    /// `SpLast` as the walk left it.
    pub sp_last: f64,
    /// `rab = 0.03 * (SpLast - SpFirst)` (`cxx:2330`).
    pub rab: f64,
    /// `Reverse`.
    pub reverse: bool,
    /// `Forward`.
    pub forward: bool,
    /// `intf`.
    pub intf: &'a mut i32,
    /// `intl`.
    pub intl: &'a mut i32,
    /// `Gd1`.
    pub gd1: &'a mut bool,
    /// `Gd2`.
    pub gd2: &'a mut bool,
    /// `Gf1`.
    pub gf1: &'a mut bool,
    /// `Gf2`.
    pub gf2: &'a mut bool,
    /// `debarc1`.
    pub debarc1: bool,
    /// `debarc2`.
    pub debarc2: bool,
    /// `debcas1`.
    pub debcas1: bool,
    /// `debcas2`.
    pub debcas2: bool,
    /// `finarc1`.
    pub finarc1: bool,
    /// `finarc2`.
    pub finarc2: bool,
    /// `fincas1`.
    pub fincas1: bool,
    /// `fincas2`.
    pub fincas2: bool,
    /// `backwContinueFailed`.
    pub backw_continue_failed: bool,
    /// `forwContinueFailed`.
    pub forw_continue_failed: bool,
    /// `unseulsuffitdeb`.
    pub unseulsuffitdeb: bool,
    /// `unseulsuffitfin`.
    pub unseulsuffitfin: bool,
}

impl ChFi3dBuilder {
    /// `cxx:2336-2505`: the extensions by tangent plane and the two `intf` /
    /// `intl` ladders. Answers nothing, like the C++ (the caller only reads the
    /// state it wrote).
    pub(super) fn apply_extensions(&self, st: &mut ExtensionState<'_>) {
        let mut line = &mut *st.line;
        let data = &mut *st.data;
        let first = &mut *st.first;
        let last = &mut *st.last;
        let intf = &mut *st.intf;
        let intl = &mut *st.intl;
        let gd1 = &mut *st.gd1;
        let gd2 = &mut *st.gd2;
        let gf1 = &mut *st.gf1;
        let gf2 = &mut *st.gf2;
        let sp_first = st.sp_first;
        let sp_last = st.sp_last;
        let rab = st.rab;
        let reverse = st.reverse;
        let forward = st.forward;
        let debarc1 = st.debarc1;
        let debarc2 = st.debarc2;
        let debcas1 = st.debcas1;
        let debcas2 = st.debcas2;
        let finarc1 = st.finarc1;
        let finarc2 = st.finarc2;
        let fincas1 = st.fincas1;
        let fincas2 = st.fincas2;
        let backw_continue_failed = st.backw_continue_failed;
        let forw_continue_failed = st.forw_continue_failed;
        let unseulsuffitdeb = st.unseulsuffitdeb;
        let unseulsuffitfin = st.unseulsuffitfin;

        // `Nbpnt` is the same shared local in the C++, so the truncation of the
        // beginning branch feeds the index bound of the end branch.
        let mut nbpnt = line.nb_points();

        // ============= POSSIBLE INVALIDATION =============
        // Preparation of the extensions by tangent plane (`cxx:2336-2420`).
        if reverse && *intf != 0 {
            *gd1 = debcas1; // skv(occ67)
            *gd2 = debcas2;
            if (debarc1 ^ debarc2) && !unseulsuffitdeb && (*first != sp_first) {
                // Case of incomplete path: the result is truncated instead of
                // the call exiting (`cxx:2344-2371`).
                let sortie = if debarc1 {
                    data.vertex_first_on_s1().parameter()
                } else {
                    data.vertex_first_on_s2().parameter()
                };
                if sortie - *first > self.tolesp {
                    let mut ind = search_index(sortie, &line) as isize;
                    if line
                        .point(ind as usize)
                        .is_some_and(|p| p.parameter() == sortie)
                    {
                        ind -= 1;
                    }
                    if ind >= 1 {
                        line.remove(1, ind as usize);
                        let _ = update_line(&mut line, true);
                    }
                    nbpnt = line.nb_points();
                    *first = point_parameter(&line, 1);
                }
            } else if (*intf >= 5) && !debarc1 && !debarc2 && (*first != sp_first) {
                let sortie = (2.0 * *first + *last) / 3.0;
                if sortie - *first > self.tolesp {
                    let mut ind = search_index(sortie, &line) as isize;
                    if line
                        .point(ind as usize)
                        .is_some_and(|p| p.parameter() == sortie)
                    {
                        ind -= 1;
                    }
                    if (nbpnt as isize) - ind < 3 {
                        ind = nbpnt as isize - 3;
                    }
                    if ind >= 1 {
                        line.remove(1, ind as usize);
                        let _ = update_line(&mut line, true);
                    }
                    nbpnt = line.nb_points();
                    *first = point_parameter(&line, 1);
                }
            }
            if *gd1 && *gd2 {
                let mut target = (point_parameter(&line, 1) - rab).min(*first);
                target = target.max(sp_first);
                data.set_first_extension_value((point_parameter(&line, 1) - target).abs());
            }
            if *intf != 0 && !unseulsuffitdeb {
                *intf = if (*gd1 && *gd2) || backw_continue_failed {
                    1
                } else {
                    0
                };
            } else if *intf != 0 && unseulsuffitdeb && (*intf < 5) {
                *intf = if *gd1 || *gd2 { 1 } else { 0 };
                // It is checked if there is no new face (`cxx:2410-2414`).
                if *intf != 0 && ((!debcas1 && debarc1) || (!debcas2 && debarc2)) {
                    *intf = 0;
                }
            } else if *intf < 5 {
                *intf = 0;
            }
        }

        if forward && *intl != 0 {
            *gf1 = fincas1; // skv(occ67)
            *gf2 = fincas2;
            if (finarc1 ^ finarc2) && !unseulsuffitfin && (*last != sp_last) {
                // Case of incomplete path: the result is truncated instead of
                // the call exiting (`cxx:2428-2455`).
                let sortie = if finarc1 {
                    data.vertex_last_on_s1().parameter()
                } else {
                    data.vertex_last_on_s2().parameter()
                };
                if *last - sortie > self.tolesp {
                    let mut ind = search_index(sortie, &line) as isize;
                    if line
                        .point(ind as usize)
                        .is_some_and(|p| p.parameter() == sortie)
                    {
                        ind += 1;
                    }
                    if ind <= nbpnt as isize {
                        line.remove(ind as usize, nbpnt);
                        let _ = update_line(&mut line, false);
                    }
                    nbpnt = line.nb_points();
                    *last = point_parameter(&line, nbpnt);
                }
            } else if (*intl >= 5) && !finarc1 && !finarc2 && (*last != sp_last) {
                // The same in the case when the entire `Lin` is an extension
                // (`cxx:2456-2480`).
                let sortie = (*first + 2.0 * *last) / 3.0;
                if *last - sortie > self.tolesp {
                    let mut ind = search_index(sortie, &line) as isize;
                    if line
                        .point(ind as usize)
                        .is_some_and(|p| p.parameter() == sortie)
                    {
                        ind += 1;
                    }
                    if ind < 3 {
                        ind = 3;
                    }
                    if ind <= nbpnt as isize {
                        line.remove(ind as usize, nbpnt);
                        let _ = update_line(&mut line, false);
                    }
                    nbpnt = line.nb_points();
                    *last = point_parameter(&line, nbpnt);
                }
            }
            if *gf1 && *gf2 {
                let mut target = (point_parameter(&line, nbpnt) + rab).max(*last);
                target = target.min(sp_last);
                data.set_last_extension_value((target - point_parameter(&line, nbpnt)).abs());
            }

            if *intl != 0 && !unseulsuffitfin {
                *intl = if (*gf1 && *gf2) || forw_continue_failed {
                    1
                } else {
                    0
                };
            } else if *intl != 0 && unseulsuffitfin && (*intl < 5) {
                *intl = if *gf1 || *gf2 { 1 } else { 0 };
                // It is checked if there is no new face (`cxx:2496-2499`).
                if *intl != 0 && ((!fincas1 && finarc1) || (!fincas2 && finarc2)) {
                    *intl = 0;
                }
            } else if *intl < 5 {
                *intl = 0;
            }
        }
    }
}
