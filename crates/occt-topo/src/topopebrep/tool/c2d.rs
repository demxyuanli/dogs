//! `TopOpeBRepTool_2d.cxx` — the `FC2D_*` p-curve helpers of the classic
//! boolean engine.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_2d.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `FC2D_Prepare` | `96-135` | [`fc2d_prepare`] |
//! | `FC2D_HasC3D` | `137-145` | [`fc2d_has_c3d`] |
//! | `FC2D_HasCurveOnSurface` | `146-155` | [`fc2d_has_curve_on_surface`] |
//! | `FC2D_HasOldCurveOnSurface` | `156-178` | [`fc2d_has_old_curve_on_surface`] |
//! | `FC2D_PNewCurveOnSurface` | `181-206` | private |
//! | `FC2D_HasNewCurveOnSurface` | `207-231` | [`fc2d_has_new_curve_on_surface`] |
//! | `FC2D_AddNewCurveOnSurface` | `233-254` | [`fc2d_add_new_curve_on_surface`] |
//! | `FC2D_make2d` (no `EF`) | `258-336` | [`fc2d_make_curve_on_surface`] |
//! | `FC2D_MakeCurveOnSurface` | `342-352` | [`fc2d_make_curve_on_surface`] |
//! | `FC2D_CurveOnSurface` | `356-375` | [`fc2d_curve_on_surface`] |
//! | `FC2D_EditableCurveOnSurface` | `379-406` | [`fc2d_editable_curve_on_surface`] |
//! | `FC2D_translate` | `410-447` | [`fc2d_translate`] |
//! | `FC2D_CurveOnSurface` (`EF` overload) | `450-566` | [`fc2d_curve_on_surface_with_ef`] |
//!
//! OCCT keeps two p-curve stores: the one `BRep_Tool::CurveOnSurface` reads
//! (the edge's own geometry, "old"), and `GLOBAL_pmosloc2df`, a temporary
//! per-edge list of curves built during an operation ("new"). They are
//! different answers for the same `(E, F)` pair, so the port mirrors the
//! second one in [`NEW_CURVES`] instead of reusing the `GeometryRegistry`
//! p-curve table, which stands for the first.
//!
//! `FC2D_AddNewCurveOnSurface` `Bind`s the edge key, which *replaces* any list
//! already stored for that edge, and only then appends (`cxx:250-253`): an edge
//! therefore holds one "new" curve per call, not one per face. The port keeps
//! that shape of the store (a `Vec` of length one per edge).
//!
//! UNPORTED:
//! - the "edge without a 3d curve" arm of `FC2D_make2d` (`cxx:307-328`): it
//!   needs `FC2D_FancestorE` (`TopOpeBRepTool_2d.cxx:57-93`, unported), and
//!   `BRepAdaptor_Curve(E, FE)` (`BRepAdaptor_Curve.cxx:91-120`) with no 3d
//!   curve behind it. The arm answers an edge that has neither a 3d nor a 2d
//!   representation — a degenerated point edge. The port answers `None` there,
//!   which is what `make2d` returns when `FC2D_FancestorE` is null at
//!   `cxx:312-315`.
//! - `MakePCurve` (`cxx:58-94`) and the `OCCT_DEBUG` blocks.

use std::cell::RefCell;
use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::{GpTrsf2d, GpVec2d};
use occt_geom2d::Curve2d;

use crate::abs::Orientation;
use crate::shape::{Edge, Face, TopoShape};
use crate::boptools_2d;
use crate::brep_tool::BRepTool;
use crate::topo_tools::{shape_map_key, ShapeMapKey};
use crate::topopebrep::tool::iso::uv_iso_curve;

/// `TopOpeBRepTool_C2DF` (`TopOpeBRepTool_C2DF.hxx:53-63`): a p-curve with its
/// `2d` range, tolerance and the face it belongs to.
#[derive(Clone)]
pub struct C2dF {
    /// `myPC`
    pub curve: Option<Arc<dyn Curve2d>>,
    /// `myf2d`
    pub first: f64,
    /// `myl2d`
    pub last: f64,
    /// `mytol`
    pub tol: f64,
    /// `myFace`
    pub face: Option<TopoShape>,
}

impl C2dF {
    /// `TopOpeBRepTool_C2DF(PC, f2d, l2d, tol, F)` (`.cxx:34-45`).
    pub fn new(curve: Arc<dyn Curve2d>, first: f64, last: f64, tol: f64, face: &TopoShape) -> Self {
        Self {
            curve: Some(curve),
            first,
            last,
            tol,
            face: Some(face.clone()),
        }
    }

    /// `PC(f2d, l2d, tol)` (`.cxx:66-74`): the curve, with its stored range and
    /// tolerance.
    pub fn pc(&self) -> (Option<Arc<dyn Curve2d>>, f64, f64, f64) {
        (self.curve.clone(), self.first, self.last, self.tol)
    }

    /// `SetPC(PC, f2d, l2d, tol)` (`.cxx:47-57`).
    pub fn set_pc(&mut self, curve: Arc<dyn Curve2d>, first: f64, last: f64, tol: f64) {
        self.curve = Some(curve);
        self.first = first;
        self.last = last;
        self.tol = tol;
    }

    /// `IsFace(F)` (`.cxx:88-93`): `F.IsEqual(myFace)`.
    pub fn is_face(&self, f: &TopoShape) -> bool {
        match &self.face {
            Some(face) => crate::topo_tools::is_equal(face, f),
            None => false,
        }
    }
}

thread_local! {
    /// `GLOBAL_pmosloc2df` (`cxx:60-64`): per-edge "new" p-curves.
    static NEW_CURVES: RefCell<HashMap<ShapeMapKey, Vec<C2dF>>> = RefCell::new(HashMap::new());
}

/// `FC2D_Prepare(S1, S2)` (`cxx:96-135`): clears the "new" store. `S1`/`S2`
/// are only kept for the debug dump in OCCT (`GLOBAL_pS1c2df` /
/// `GLOBAL_pS2c2df`, read nowhere else).
pub fn fc2d_prepare() {
    NEW_CURVES.with(|m| m.borrow_mut().clear());
}

/// `FC2D_HasC3D(E)` (`cxx:137-145`).
pub fn fc2d_has_c3d(e: &TopoShape) -> bool {
    BRepTool::edge_curve_world(&Edge(e.clone())).is_some()
}

/// `FC2D_HasCurveOnSurface(E, F)` (`cxx:146-155`).
pub fn fc2d_has_curve_on_surface(e: &TopoShape, f: &TopoShape) -> bool {
    fc2d_has_old_curve_on_surface(e, f).is_some() || fc2d_has_new_curve_on_surface(e, f).is_some()
}

/// `FC2D_HasOldCurveOnSurface(E, F, C2D, f2d, l2d, tol)` (`cxx:156-169`):
/// `BRep_Tool::CurveOnSurface(E, F, f2d, l2d)` with `tol = Tolerance(E)`.
pub fn fc2d_has_old_curve_on_surface(e: &TopoShape, f: &TopoShape) -> Option<C2dF> {
    let tol = BRepTool::edge_tolerance(&Edge(e.clone()));
    let (curve, first, last) =
        boptools_2d::curve_on_surface_range(&Edge(e.clone()), &Face(f.clone()))?;
    Some(C2dF {
        curve: Some(curve),
        first,
        last,
        tol,
        face: Some(f.clone()),
    })
}

/// `FC2D_PNewCurveOnSurface(E, F)` (`cxx:181-206`).
fn fc2d_new_curve_on_surface(e: &TopoShape, f: &TopoShape) -> Option<C2dF> {
    NEW_CURVES.with(|m| {
        let map = m.borrow();
        let list = map.get(&shape_map_key(e))?;
        list.iter().find(|c2df| c2df.is_face(f)).cloned()
    })
}

/// `FC2D_HasNewCurveOnSurface(E, F, C2D, f2d, l2d, tol)` (`cxx:207-231`).
pub fn fc2d_has_new_curve_on_surface(e: &TopoShape, f: &TopoShape) -> Option<C2dF> {
    fc2d_new_curve_on_surface(e, f)
}

/// `FC2D_AddNewCurveOnSurface(C2D, E, F, f2d, l2d, tol)` (`cxx:233-254`).
///
/// Answers `1` for a null curve or an unprepared store, `0` otherwise.
pub fn fc2d_add_new_curve_on_surface(
    curve: &Arc<dyn Curve2d>,
    e: &TopoShape,
    f: &TopoShape,
    first: f64,
    last: f64,
    tol: f64,
) -> i32 {
    let c2df = C2dF::new(curve.clone(), first, last, tol, f);
    NEW_CURVES.with(|m| {
        let mut map = m.borrow_mut();
        map.insert(shape_map_key(e), vec![c2df]);
    });
    0
}

/// `FC2D_make2d(E, F, f2d, l2d, tol)` (`cxx:258-336`), without the
/// `FC2D_AddNewCurveOnSurface` that `FC2D_MakeCurveOnSurface` adds.
///
/// The curve is only built when neither store has one (`cxx:270-274`), and
/// only for an edge that carries a 3d curve (`cxx:282-305`). `trim3d` picks
/// the 3d range as the projection window (`cxx:294-300`); the port's
/// [`boptools_2d::make_2d`] projects over the whole edge either way, which is
/// what `MakePCurveOnFace` does for an untrimmed input.
fn fc2d_make2d(e: &TopoShape, f: &TopoShape, trim3d: bool) -> Option<C2dF> {
    let _ = trim3d;
    if let Some(old) = fc2d_has_old_curve_on_surface(e, f) {
        return Some(old);
    }
    let (f3d, l3d) = BRepTool::edge_parameters(&Edge(e.clone()));
    // `cxx:282`: no 3d curve, no projection.
    if !fc2d_has_c3d(e) {
        return None;
    }
    let curve = boptools_2d::make_2d(&Edge(e.clone()), &Face(f.clone())).ok()?;
    let tol = BRepTool::edge_tolerance(&Edge(e.clone()));
    Some(C2dF::new(curve, f3d, l3d, tol, f))
}

/// `FC2D_make2d(E, F, f2d, l2d, tol)` (`cxx:258-336`) followed by
/// `FC2D_AddNewCurveOnSurface` — `FC2D_MakeCurveOnSurface` (`cxx:342-352`).
pub fn fc2d_make_curve_on_surface(
    e: &TopoShape,
    f: &TopoShape,
    trim3d: bool,
) -> Option<C2dF> {
    let c2df = fc2d_make2d(e, f, trim3d)?;
    if let Some(curve) = &c2df.curve {
        fc2d_add_new_curve_on_surface(curve, e, f, c2df.first, c2df.last, c2df.tol);
    }
    Some(c2df)
}

/// `FC2D_CurveOnSurface(E, F, f, l, tol, trim3d)` (`cxx:356-375`): old, then
/// new, then build.
pub fn fc2d_curve_on_surface(e: &TopoShape, f: &TopoShape, trim3d: bool) -> Option<C2dF> {
    if let Some(old) = fc2d_has_old_curve_on_surface(e, f) {
        return Some(old);
    }
    if let Some(new) = fc2d_has_new_curve_on_surface(e, f) {
        return Some(new);
    }
    fc2d_make_curve_on_surface(e, f, trim3d)
}

/// `FC2D_EditableCurveOnSurface(E, F, f, l, tol, trim3d)` (`cxx:379-406`):
/// the old curve is handed out as a *copy*, because the caller edits it.
pub fn fc2d_editable_curve_on_surface(e: &TopoShape, f: &TopoShape, trim3d: bool) -> Option<C2dF> {
    // `cxx:386-395`
    if let Some(mut old) = fc2d_has_old_curve_on_surface(e, f) {
        if let Some(pc) = &old.curve {
            old.curve = Some(Arc::from(pc.clone_dyn()));
        }
        return Some(old);
    }
    // `cxx:396-403`
    if let Some(new) = fc2d_has_new_curve_on_surface(e, f) {
        return Some(new);
    }
    // `cxx:404`
    fc2d_make_curve_on_surface(e, f, trim3d)
}

/// `FC2D_translate(C2D, E, F, EF)` (`cxx:410-447`): on a periodic surface,
/// shift an iso-parametric p-curve so its origin matches the origin that `EF`'s
/// p-curve on `F` starts at.
///
/// Both early exits of OCCT are kept: `UVISO` (`cxx:421`) answers nothing for a
/// curve that is not a line parallel to a parametric axis, and a null `EF`
/// stops the block at `cxx:423`. The `|factor| > 1.e-6` gate (`cxx:435`) then
/// skips an already-aligned curve.
///
/// OCCT mutates the handle in place (`C2D->Translate`); the port's curves are
/// immutable, so the translated copy is returned and the caller re-stores it
/// where OCCT's `SetPC` would.
pub fn fc2d_translate(
    c2d: Arc<dyn Curve2d>,
    f: &TopoShape,
    ef: Option<&TopoShape>,
) -> Arc<dyn Curve2d> {
    let Some(surface) = BRepTool::face_surface_world(&Face(f.clone())) else {
        return c2d;
    };
    // `cxx:417-418`: `S1->IsUPeriodic() || S1->IsVPeriodic()`.
    let is_perio = surface.is_u_periodic() || surface.is_v_periodic();
    let Some(iso) = uv_iso_curve(c2d.as_ref()) else {
        return c2d;
    };
    if !is_perio {
        return c2d;
    }
    let Some(ef) = ef else {
        return c2d;
    };

    // `cxx:428-429`: `TopoDS_Face FFOR = F; FFOR.Orientation(TopAbs_FORWARD);`
    let mut ffor = f.clone();
    ffor.set_orientation(Orientation::Forward);
    let (p1, _p2) = BRepTool::uv_points(&Edge(ef.clone()), &Face(ffor));
    // `cxx:432-434`
    let p_ef = if iso.iso_u { p1.x() } else { p1.y() };
    let p_c2d = if iso.iso_u {
        iso.origin.x()
    } else {
        iso.origin.y()
    };
    let factor = p_ef - p_c2d;
    if factor.abs() <= 1e-6 {
        return c2d;
    }
    // `cxx:438-443`: `gp_Vec2d(1., 0.)`, or `(0., 1.)` for an iso-V curve.
    let transl = if iso.iso_v {
        GpVec2d::new(0.0, factor)
    } else {
        GpVec2d::new(factor, 0.0)
    };
    let mut a_tr = GpTrsf2d::identity();
    a_tr.set_translation_vec(&transl);
    Arc::from(c2d.transformed(&a_tr))
}

/// `pc2df->SetPC(C2D, f2d, l2d, tol)` (`TopOpeBRepTool_C2DF.cxx:47-57`) on the
/// entry `FC2D_PNewCurveOnSurface` returned.
fn fc2d_store_new_curve(c2df: &C2dF, e: &TopoShape) {
    NEW_CURVES.with(|m| {
        m.borrow_mut().insert(shape_map_key(e), vec![c2df.clone()]);
    });
}

/// `FC2D_CurveOnSurface(E, F, EF, f2d, l2d, tol, trim3d)` (`cxx:450-566`): the
/// seven-argument overload used by `TopOpeBRepBuild_Section.cxx:684`.
///
/// `FC2D_make2d` has a second body in OCCT whose only difference is the
/// `FC2D_translate` calls at `cxx:499` and `cxx:524`. Both sit after the
/// "already stored" early return (`cxx:467-470`), so one [`fc2d_make2d`] plus
/// the translation here reproduces the same control flow.
pub fn fc2d_curve_on_surface_with_ef(
    e: &TopoShape,
    f: &TopoShape,
    ef: Option<&TopoShape>,
    trim3d: bool,
) -> Option<C2dF> {
    // `cxx:548-551`
    if let Some(old) = fc2d_has_old_curve_on_surface(e, f) {
        return Some(old);
    }
    // `cxx:554-561`: the stored "new" curve is translated, then re-stored.
    if let Some(mut new) = fc2d_new_curve_on_surface(e, f) {
        if let Some(pc) = &new.curve {
            new.curve = Some(fc2d_translate(pc.clone(), f, ef));
        }
        fc2d_store_new_curve(&new, e);
        return Some(new);
    }
    // `cxx:563-565`: `FC2D_make2d(E, F, EF, ...)` then `Add`.
    let mut c2df = fc2d_make2d(e, f, trim3d)?;
    if let Some(pc) = &c2df.curve {
        let translated = fc2d_translate(pc.clone(), f, ef);
        c2df.curve = Some(translated.clone());
        fc2d_add_new_curve_on_surface(&translated, e, f, c2df.first, c2df.last, c2df.tol);
    }
    Some(c2df)
}
