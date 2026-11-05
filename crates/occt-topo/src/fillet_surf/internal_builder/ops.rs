//! The four `FilletSurf_InternalBuilder` overrides of the `ChFi3d_FilBuilder`
//! virtuals.
//!
//! Source: `ModelingAlgorithms/TKFillet/FilletSurf/FilletSurf_InternalBuilder.cxx:341-596`
//! (OCCT 8.0.0).
//!
//! OCCT overrides the four `PerformSurf` overloads of `ChFi3d_FilBuilder`
//! (`FilletSurf_InternalBuilder.hxx:131-235`):
//! - the face / face one runs a constant-radius walk (`cxx:341-503`),
//! - the three edge-contact ones only `throw Standard_DomainError(
//!   "BlendFunc_CSConstRad::Section : Not implemented")` (`cxx:505-596`).
//!
//! The port carries the overridden half in [`FilletSurfOps`] and forwards
//! everything else to [`FilSurfOps`], which is the `ChFi3d_FilBuilder` half
//! OCCT leaves untouched here.
//!
//! Divergences that are not heuristics:
//! - `throw Standard_ConstructionError("PerformSurf : this is not the spine of
//!   a fillet")` (`cxx:366-369`), `throw Standard_ConstructionError(
//!   "PerformSurf : no variable radiuses")` (`cxx:372-375`) and `throw
//!   Standard_Failure("PerformSurf : Failed approximation!")` (`cxx:456-459`)
//!   are the trait's `false`, which is what the caller's recovery path reacts
//!   to,
//! - `SeqData(1)` on an empty sequence raises `Standard_OutOfRange` in OCCT
//!   (`cxx:364`); the port answers `false` on the same arm,
//! - `S1->Face().Orientation()` raises `Standard_NullObject` for an adaptor
//!   without a face; the port answers `false` on the same arm,
//! - the three edge-contact overrides are the trait defaults, which already
//!   answer `Err` on every call (`ChFi3d_Builder_NotImp.cxx:125-220` is the
//!   body OCCT's base class has for them), so nothing is written here,
//! - `BRepBlend_ConstRad` / `BRepBlend_ConstRadInv` are [`ConstRad`] /
//!   [`ConstRadInv`] plus `Func.Set(myShape)`, which the port spells
//!   `set_section_shape`.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_geom::Curve;
use occt_math::MathVector;

use crate::blend::blend_func::BlendFuncSectionShape;
use crate::blend::{ConstRad, ConstRadInv};
use crate::brep_blend::Line;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi3d::builder::FilSurfOps;
use crate::chfi3d::builder::simul_data_face_face;
use crate::chfi3d::builder::{guide_as_curve, ChFi3dBuilder, ChFi3dSurfOps};
use crate::chfi3d::builder_helpers::fil_common_point::fil_common_point;
use crate::chfi_ds::{ElSpine, SpineClass, SpineHandle, SurfData};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `FilletSurf_InternalBuilder`'s share of `ChFi3d_FilBuilder`.
///
/// `myShape` (`ChFi3d_FilBuilder.hxx:377`) is the only member the subclass
/// body reads; it is carried here for the same reason [`FilSurfOps`] carries
/// it — see the module note of `chfi3d::builder::fil_builder`.
pub struct FilletSurfOps {
    /// `BlendFunc_SectionShape myShape;` (`ChFi3d_FilBuilder.hxx:377`).
    shape: BlendFuncSectionShape,
    /// The untouched `ChFi3d_FilBuilder` half of the virtuals.
    fil: FilSurfOps,
}

impl FilletSurfOps {
    /// A carrier for `FilSurfOps` built over the same `myShape`.
    pub fn new(shape: BlendFuncSectionShape) -> Self {
        Self {
            shape,
            fil: FilSurfOps::new(shape),
        }
    }
}

impl ChFi3dSurfOps for FilletSurfOps {
    /// `SimulKPart(SD)` — not overridden by `FilletSurf_InternalBuilder`.
    fn simul_k_part(
        &mut self,
        sd: &mut SurfData,
        ds: &crate::topopebrep::ds::HDataStructure,
    ) -> Result<(), String> {
        self.fil.simul_k_part(sd, ds)
    }

    /// `SimulSurf(...)` on two faces — not overridden here.
    #[allow(clippy::too_many_arguments)]
    fn simul_surf(
        &mut self,
        data: &mut SurfData,
        guide: &mut Arc<ElSpine>,
        spine: Option<&SpineHandle>,
        choix: i32,
        s1: &SurfaceAdaptor,
        i1: &Rc<RefCell<TopolTool>>,
        s2: &SurfaceAdaptor,
        i2: &Rc<RefCell<TopolTool>>,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        approx: bool,
        forward: bool,
        rec_on_s1: bool,
        rec_on_s2: bool,
        soldep: &MathVector,
        intf: &mut i32,
        intl: &mut i32,
        ctx: &crate::chfi3d::builder::ChFi3dCtx<'_>,
    ) -> bool {
        self.fil.simul_surf(
            data, guide, spine, choix, s1, i1, s2, i2, tol_guide, first, last, inside, approx,
            forward, rec_on_s1, rec_on_s2, soldep, intf, intl, ctx,
        )
    }

    /// `PerformFirstSection(...)` — not overridden here.
    #[allow(clippy::too_many_arguments)]
    fn perform_first_section(
        &mut self,
        spine: &SpineHandle,
        guide: &Arc<ElSpine>,
        choix: i32,
        s1: &SurfaceAdaptor,
        s2: &SurfaceAdaptor,
        i1: &Rc<RefCell<TopolTool>>,
        i2: &Rc<RefCell<TopolTool>>,
        par: f64,
        sol_dep: &mut MathVector,
        pos1: &mut crate::abs::State,
        pos2: &mut crate::abs::State,
        ctx: &crate::chfi3d::builder::ChFi3dCtx<'_>,
    ) -> bool {
        self.fil
            .perform_first_section(spine, guide, choix, s1, s2, i1, i2, par, sol_dep, pos1, pos2, ctx)
    }

    /// `PerformSurf(SeqData, Guide, Spine, Choix, S1, I1, S2, I2, MaxStep,
    /// Fleche, TolGuide, First, Last, Inside, Appro, Forward, RecOnS1, RecOnS2,
    /// Soldep, Intf, Intl)` (`cxx:341-503`).
    #[allow(clippy::too_many_arguments)]
    fn perform_surf(
        &mut self,
        base: &mut ChFi3dBuilder,
        data: &mut Vec<Arc<SurfData>>,
        guide: &mut Arc<ElSpine>,
        spine: Option<&SpineHandle>,
        choix: i32,
        s1: &SurfaceAdaptor,
        i1: &Rc<RefCell<TopolTool>>,
        s2: &SurfaceAdaptor,
        i2: &Rc<RefCell<TopolTool>>,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        approx: bool,
        forward: bool,
        rec_on_s1: bool,
        rec_on_s2: bool,
        soldep: &MathVector,
        intf: &mut i32,
        intl: &mut i32,
    ) -> bool {
        // `occ::handle<ChFiDS_SurfData> Data = SeqData(1);` (`cxx:365`).
        if data.is_empty() {
            return false;
        }
        // `occ::handle<ChFiDS_FilSpine> fsp =
        // occ::down_cast<ChFiDS_FilSpine>(Spine); if (fsp.IsNull()) throw;`
        // (`cxx:366-369`).
        let Some(spine_handle) = spine else {
            return false;
        };
        let (constant, radius) = match spine_handle.read() {
            Ok(guard) => {
                if guard.class() != SpineClass::Fil {
                    return false;
                }
                // `if (!fsp->IsConstant()) throw ...;` (`cxx:372-375`).
                (guard.is_fil_constant(), guard.fil_radius())
            }
            Err(_) => return false,
        };
        if !constant {
            return false;
        }
        let Some(radius) = radius else {
            return false;
        };

        // `TopAbs_Orientation Or = S1->Face().Orientation();` (`cxx:371`).
        let Some(or) = s1.face().map(|f| f.0.orientation()) else {
            return false;
        };

        // `BRepBlend_ConstRad Func(S1, S2, Guide); BRepBlend_ConstRadInv
        // FInv(S1, S2, Guide); Func.Set(fsp->Radius(), Choix);
        // FInv.Set(fsp->Radius(), Choix);` (`cxx:380-384`).
        let Some(s1_surf) = s1.surface().cloned() else {
            return false;
        };
        let Some(s2_surf) = s2.surface().cloned() else {
            return false;
        };
        let guide_curve: Arc<dyn Curve> = guide_as_curve(guide);
        let mut func = ConstRad::new(
            Arc::clone(&s1_surf),
            Arc::clone(&s2_surf),
            Arc::clone(&guide_curve),
        );
        let mut finv = ConstRadInv::new(Arc::clone(&s1_surf), Arc::clone(&s2_surf), guide_curve);
        func.set_radius(radius, choix);
        finv.set_radius(radius, choix);
        // `switch (GetFilletShape()) { ... Func.Set(...); }` (`cxx:385-394`).
        func.set_section_shape(self.shape);

        // `double PFirst = First;` (`cxx:395`).
        let p_first = *first;

        // `occ::handle<ChFiDS_ElSpine> EmptyGuide;` (`cxx:378`).
        let mut empty_guide: Option<Arc<ElSpine>> = None;
        let mut lin: Option<Rc<RefCell<Line>>> = None;

        // `done = SimulData(Data, Guide, EmptyGuide, lin, S1, I1, S2, I2, Func,
        // FInv, PFirst, MaxStep, Fleche, TolGuide, First, Last, Inside, Appro,
        // Forward, Soldep, 20, RecOnS1, RecOnS2);` (`cxx:396-419`).
        let done = simul_data_face_face(
            base.tolapp3d(),
            guide,
            empty_guide.as_mut(),
            &mut lin,
            s1_surf,
            Rc::clone(i1),
            s2_surf,
            Rc::clone(i2),
            &mut func,
            &mut finv,
            p_first,
            max_step,
            fleche,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            soldep,
            20,
            rec_on_s1,
            rec_on_s2,
        );
        // `if (!done) { return false; }` (`cxx:419-422`).
        if !done {
            return false;
        }
        let Some(lin) = lin else {
            return false;
        };

        // The four `ChFi3d_FilCommonPoint(...)` calls (`cxx:423-455`).
        let (trans1, trans2) = {
            let line = lin.borrow();
            (line.transition_on_s1(), line.transition_on_s2())
        };
        let (Ok(trans1), Ok(trans2)) = (trans1, trans2) else {
            return false;
        };
        {
            let line = lin.borrow();
            let data = Arc::make_mut(&mut data[0]);
            if line.start_point_on_first().nb_point_on_rst() != 0 {
                fil_common_point(
                    line.start_point_on_first(),
                    trans1,
                    true,
                    data.vertex_first_on_s1_mut(),
                    base.tolapp3d(),
                );
            }
            if line.end_point_on_first().nb_point_on_rst() != 0 {
                fil_common_point(
                    line.end_point_on_first(),
                    trans1,
                    false,
                    data.vertex_last_on_s1_mut(),
                    base.tolapp3d(),
                );
            }
            if line.start_point_on_second().nb_point_on_rst() != 0 {
                fil_common_point(
                    line.start_point_on_second(),
                    trans2,
                    true,
                    data.vertex_first_on_s2_mut(),
                    base.tolapp3d(),
                );
            }
            if line.end_point_on_second().nb_point_on_rst() != 0 {
                fil_common_point(
                    line.end_point_on_second(),
                    trans2,
                    false,
                    data.vertex_last_on_s2_mut(),
                    base.tolapp3d(),
                );
            }
        }

        // `done = CompleteData(Data, Func, lin, S1, S2, Or, false, false,
        // false, false);` (`cxx:456`). `Reversed` keeps its `false` default.
        let func = Rc::new(RefCell::new(func));
        let done = {
            let data = Arc::make_mut(&mut data[0]);
            base.complete_data_app_func(
                data,
                Rc::clone(&func),
                Rc::clone(&lin),
                s1,
                Some(s2),
                or,
                false,
                false,
                false,
                false,
                false,
            )
        };
        // `if (!done) { throw Standard_Failure("PerformSurf : Failed
        // approximation!"); }` (`cxx:456-459`) — see the module note.
        if !matches!(done, Ok(true)) {
            return false;
        }

        // `if (!Forward) { ... }` (`cxx:461-479`).
        let mut ok = false;
        if !forward {
            *intf = 0;
            let cp1 = data[0].vertex_first_on_s1().clone();
            if cp1.is_on_arc() {
                let f1 = s1.face().cloned();
                let mut bid = None;
                let found = base.search_face(Some(spine_handle), &cp1, f1.as_ref(), &mut bid);
                *intf = if found { 0 } else { 1 };
                ok = *intf != 0;
            }
            let cp2 = data[0].vertex_first_on_s2().clone();
            if cp2.is_on_arc() && !ok {
                let f2 = s2.face().cloned();
                let mut bid = None;
                let found = base.search_face(Some(spine_handle), &cp2, f2.as_ref(), &mut bid);
                *intf = if found { 0 } else { 1 };
            }
        }
        // `Intl = 0; ok = false; const ChFiDS_CommonPoint& cpl1 =
        // Data->VertexLastOnS1(); ...` (`cxx:480-497`).
        *intl = 0;
        ok = false;
        let cp1 = data[0].vertex_last_on_s1().clone();
        if cp1.is_on_arc() {
            let f1 = s1.face().cloned();
            let mut bid = None;
            let found = base.search_face(Some(spine_handle), &cp1, f1.as_ref(), &mut bid);
            *intl = if found { 0 } else { 1 };
            ok = *intl != 0;
        }
        let cp2 = data[0].vertex_last_on_s2().clone();
        if cp2.is_on_arc() && !ok {
            let f2 = s2.face().cloned();
            let mut bid = None;
            let found = base.search_face(Some(spine_handle), &cp2, f2.as_ref(), &mut bid);
            *intl = if found { 0 } else { 1 };
        }

        // `Data->FirstSpineParam(First); Data->LastSpineParam(Last);`
        // (`cxx:498-499`).
        {
            let data = Arc::make_mut(&mut data[0]);
            data.set_first_spine_param(*first);
            data.set_last_spine_param(*last);
        }
        true
    }
}
