//! The four `Localize*` scans shared by `LocOpe_CSIntersector`
//! (`LocOpe_CSIntersector.cxx:316-513`) and `LocOpe_CurveShapeIntersector`
//! (`LocOpe_CurveShapeIntersector.cxx:75-354`).
//!
//! The two classes carry the same four algorithms: `CSIntersector` forwards its
//! `Tol` argument, `CurveShapeIntersector` always passes
//! `Precision::Confusion()`. The difference is only that argument, so the scans
//! live here once and both classes call them.
//!
//! Each function walks the sorted point list from a parameter or from an index,
//! merges points whose parameters differ by at most `tol`, and reports the first
//! merged group whose orientation is not `TopAbs_EXTERNAL`. Indices are the
//! reference's, i.e. **1-based** (`IndFrom` / `IndTo` index the sequence the
//! callers hold).

use crate::abs::Orientation;

use super::pnt_face::PntFace;

/// One localization result: `(Or, IndFrom, IndTo)` with 1-based indices.
pub type Localized = (Orientation, usize, usize);

/// `LocAfter(Spt, From, Tol, Or, IndFrom, IndTo)`
/// (`CSIntersector.cxx:316-378`, `CurveShapeIntersector.cxx:75-133`).
pub(crate) fn loc_after_param(points: &[PntFace], from: f64, tol: f64) -> Option<Localized> {
    let nb = points.len();
    let from_eps = from - tol;
    let mut ifirst = 1usize;
    while ifirst <= nb && points[ifirst - 1].parameter() < from_eps {
        ifirst += 1;
    }
    if ifirst > nb {
        return None;
    }
    let mut i = ifirst;
    let mut ind_from = ifirst;
    loop {
        if i > nb {
            return None;
        }
        let mut or = points[i - 1].orientation();
        let param = points[i - 1].parameter();
        i += 1;
        while i <= nb {
            if points[i - 1].parameter() - param <= tol {
                if or != Orientation::External && or != points[i - 1].orientation() {
                    or = Orientation::External;
                }
                i += 1;
            } else {
                break;
            }
        }
        if or == Orientation::External {
            ind_from = i;
            if i > nb {
                return None;
            }
        } else {
            let ind_to = i - 1;
            return Some((or, ind_from, ind_to));
        }
    }
}

/// `LocBefore(Spt, From, Tol, Or, IndFrom, IndTo)`
/// (`CSIntersector.cxx:380-438`, `CurveShapeIntersector.cxx:136-197`).
pub(crate) fn loc_before_param(points: &[PntFace], from: f64, tol: f64) -> Option<Localized> {
    let nb = points.len();
    let from_eps = from + tol;
    let mut ifirst = nb;
    while ifirst >= 1 && points[ifirst - 1].parameter() > from_eps {
        ifirst -= 1;
    }
    if ifirst < 1 {
        return None;
    }
    let mut i = ifirst;
    let mut ind_to = ifirst;
    loop {
        if i < 1 || i > nb {
            return None;
        }
        let mut or = points[i - 1].orientation();
        let param = points[i - 1].parameter();
        i -= 1;
        while i >= 1 {
            if param - points[i - 1].parameter() <= tol {
                if or != Orientation::External && or != points[i - 1].orientation() {
                    or = Orientation::External;
                }
                i -= 1;
            } else {
                break;
            }
        }
        if or == Orientation::External {
            ind_to = i;
            if i < 1 {
                return None;
            }
        } else {
            let ind_from = i + 1;
            return Some((or, ind_from, ind_to));
        }
    }
}

/// `LocAfter(Spt, FromInd, Tol, Or, IndFrom, IndTo)`
/// (`CSIntersector.cxx:442-513`; `CurveShapeIntersector.cxx:201-275`).
pub(crate) fn loc_after_index(points: &[PntFace], from_ind: usize, tol: f64) -> Option<Localized> {
    let nb = points.len();
    if from_ind >= nb {
        return None;
    }
    let ifirst = if from_ind >= 1 {
        let from_eps = points[from_ind - 1].parameter() - tol;
        let mut i = from_ind + 1;
        while i <= nb && points[i - 1].parameter() < from_eps {
            i += 1;
        }
        i
    } else {
        1
    };
    if ifirst > nb {
        return None;
    }
    let mut i = ifirst;
    let mut ind_from = ifirst;
    loop {
        if i > nb {
            return None;
        }
        let mut or = points[i - 1].orientation();
        let param = points[i - 1].parameter();
        i += 1;
        while i <= nb {
            if points[i - 1].parameter() - param <= tol {
                if or != Orientation::External && or != points[i - 1].orientation() {
                    or = Orientation::External;
                }
                i += 1;
            } else {
                break;
            }
        }
        if or == Orientation::External {
            ind_from = i;
            if i > nb {
                return None;
            }
        } else {
            let ind_to = i - 1;
            return Some((or, ind_from, ind_to));
        }
    }
}

/// `LocBefore(Spt, FromInd, Tol, Or, IndFrom, IndTo)`
/// (`CSIntersector.cxx:441-513`; `CurveShapeIntersector.cxx:279-354`).
pub(crate) fn loc_before_index(points: &[PntFace], from_ind: usize, tol: f64) -> Option<Localized> {
    let nb = points.len();
    if from_ind <= 1 {
        return None;
    }
    let ifirst = if from_ind <= nb {
        let from_eps = points[from_ind - 1].parameter() + tol;
        let mut i = from_ind - 1;
        while i >= 1 && points[i - 1].parameter() > from_eps {
            i -= 1;
        }
        i
    } else {
        nb
    };
    if ifirst < 1 {
        return None;
    }
    let mut i = ifirst;
    let mut ind_to = ifirst;
    loop {
        if i < 1 || i > nb {
            return None;
        }
        let mut or = points[i - 1].orientation();
        let param = points[i - 1].parameter();
        i -= 1;
        while i >= 1 {
            if param - points[i - 1].parameter() <= tol {
                if or != Orientation::External && or != points[i - 1].orientation() {
                    or = Orientation::External;
                }
                i -= 1;
            } else {
                break;
            }
        }
        if or == Orientation::External {
            ind_to = i;
            if i < 1 {
                return None;
            }
        } else {
            let ind_from = i + 1;
            return Some((or, ind_from, ind_to));
        }
    }
}
