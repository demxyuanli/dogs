//! `ChFi3d_FilBuilder::SimulKPart` and the `myShape` carrier.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder.cxx:475-554`
//! (OCCT 8.0.0).
//!
//! OCCT implements `SimulKPart` on the subclass because it needs `myDS` (the
//! base data structure) and nothing else; the four `SimulSurf` / four
//! `PerformSurf` bodies read only `myShape` (`cxx:558-2273`, verified by
//! grep). The port therefore splits the subclass member out of the per-call
//! carrier: [`FilSurfOps`] holds `myShape` (and whatever the later batches
//! need), while the base data structure is passed to `SimulKPart` explicitly.
//!
//! UNPORTED — nothing of `ChFi3dSurfOps` is left: all four `PerformSurf` and
//! the two edge / face `SimulSurf` are ported, and the two edge / edge
//! `SimulSurf` overloads keep the base class' `Err` because OCCT does not
//! override them for the fillet builder either.
//!
//! `SimulSurf` on two faces (`cxx:558-787`) lives in
//! [`super::simul_surf`], the two edge / face ones (`cxx:787-1248`) in
//! [`super::simul_surf_rst`], `PerformFirstSection` (`cxx:1500-1537`) in
//! [`super::perform_first_section`], `PerformSurf` on two faces
//! (`cxx:1538-1717`) in [`super::perform_surf`], its two edge / face
//! overloads (`cxx:1721-2063`) in [`super::perform_surf_rst`] and the edge /
//! edge one (`cxx:2065-2270`) in [`super::perform_surf_rst_rst`].

use std::sync::Arc;

use occt_core::elib::slib;
use occt_core::gp::GpPnt2d;
use crate::abs::State;
use crate::blend::blend_func::BlendFuncSectionShape;
use crate::brep_surface::{classify_surface, SurfaceKind};
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi_ds::{CircSectionArray, ElSpine, SpineHandle, StripeHandle, SurfData};
use crate::shape::TopoShape;
use crate::topopebrep::ds::HDataStructure;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;
use occt_math::MathVector;

use super::super::corner_ops::ChFi3dCornerOps;
use super::super::surf_ops::ChFi3dSurfOps;
use super::super::{ChFi3dBuilder, OneCornerExit};
use super::perform_first_section::perform_first_section;
use super::perform_surf::perform_surf_face_face;
use super::perform_surf_rst::{perform_surf_on_s1, perform_surf_on_s2};
use super::perform_surf_rst_rst::perform_surf_contact;
use super::simul_surf::simul_surf_face_face;
use super::simul_surf_rst::{simul_surf_on_s1, simul_surf_on_s2};
use crate::chfi3d::builder::ChFi3dCtx;

/// The subclass half of `ChFi3d_FilBuilder`: `myShape` (`hxx:377`) plus, once
/// the walkers are ported, the scratch the `PerformSurf` bodies keep.
///
/// It is a separate object rather than a field of
/// [`super::ChFi3dFilBuilder`] so the base can be mutably borrowed while it
/// holds this as `&mut dyn ChFi3dSurfOps` — see the module note.
pub struct FilSurfOps {
    /// `BlendFunc_SectionShape myShape;` (`hxx:377`).
    shape: BlendFuncSectionShape,
}

impl FilSurfOps {
    /// The carrier of `myShape`.
    pub fn new(shape: BlendFuncSectionShape) -> Self {
        Self { shape }
    }

    /// The current `myShape`.
    pub fn shape(&self) -> BlendFuncSectionShape {
        self.shape
    }

    /// `SetFilletShape` written through the carrier.
    pub fn set_shape(&mut self, shape: BlendFuncSectionShape) {
        self.shape = shape;
    }
}

/// `ChFi3d_FilBuilder::SimulKPart(const handle(ChFiDS_SurfData)& SD)`
/// (`cxx:475-554`): fills `SD`'s simplification sections from the DS surface.
///
/// `DStr.Surface(SD->Surf()).Surface()` — a null surface leaves `sec` null in
/// OCCT (the `switch` is skipped and `SetSimul(null)` is stored); the port
/// answers `Err` because a null surface has no `GetType()` to switch on.
///
/// `sec->ChangeValue(i)` with `1 <= i <= n` is bit-for-bit the port's
/// `change_value`.
fn simul_k_part_impl(sd: &mut SurfData, ds: &HDataStructure) -> Result<(), String> {
    let surface = ds
        .ds()
        .surface(sd.surf())
        .surface()
        .ok_or_else(|| "SimulKPart: DStr.Surface(SD->Surf()).Surface() is null".to_string())?
        .clone();

    // `gp_Pnt2d p1f = SD->InterferenceOnS1().PCurveOnSurf()->Value(...);`
    // (and the three siblings, `cxx:479-486`).
    let p1f = p_curve_point(sd, 1, true)?;
    let p1l = p_curve_point(sd, 1, false)?;
    let p2f = p_curve_point(sd, 2, true)?;
    let p2l = p_curve_point(sd, 2, false)?;

    let mut sec: Option<CircSectionArray> = None;
    match classify_surface(surface.as_ref()) {
        // `case GeomAbs_Cylinder:` (`cxx:493-505`).
        SurfaceKind::Cylinder => {
            let u1 = p1f.x();
            let u2 = p2f.x();
            let v1 = p1f.y().max(p2f.y());
            let v2 = p1l.y().min(p2l.y());
            let cy = surface
                .gp_cylinder()
                .ok_or_else(|| "SimulKPart: GeomAdaptor_Surface::Cylinder() on a cylinder".to_string())?;
            let mut arr = CircSectionArray::new(1, 2);
            arr.change_value(1).unwrap().set_circ(
                &slib::cylinder_v_iso(&cy.position(), cy.radius(), v1),
                u1,
                u2,
            );
            arr.change_value(2).unwrap().set_circ(
                &slib::cylinder_v_iso(&cy.position(), cy.radius(), v2),
                u1,
                u2,
            );
            sec = Some(arr);
        }
        // `case GeomAbs_Torus:` (`cxx:506-527`).
        SurfaceKind::Torus => {
            let v1 = p1f.y();
            let v2 = p2f.y();
            let u1 = p1f.x().max(p2f.x());
            let u2 = p1l.x().min(p2l.x());
            let ang = u2 - u1;
            let to = surface
                .gp_torus()
                .ok_or_else(|| "SimulKPart: GeomAdaptor_Surface::Torus() on a torus".to_string())?;
            let (majr, minr) = (to.major_radius, to.minor_radius);
            let n = section_count(ang);
            let mut arr = CircSectionArray::new(1, n);
            for i in 1..=n {
                let u = u1 + f64::from(i - 1) * (u2 - u1) / f64::from(n - 1);
                arr.change_value(i).unwrap().set_circ(
                    &slib::torus_u_iso(&to.pos, majr, minr, u),
                    v1,
                    v2,
                );
            }
            sec = Some(arr);
        }
        // `case GeomAbs_Sphere:` (`cxx:528-549`).
        SurfaceKind::Sphere => {
            let v1 = p1f.y();
            let v2 = p2f.y();
            let u1 = p1f.x().max(p2f.x());
            let u2 = p1l.x().min(p2l.x());
            let ang = u2 - u1;
            let sp = surface
                .gp_sphere()
                .ok_or_else(|| "SimulKPart: GeomAdaptor_Surface::Sphere() on a sphere".to_string())?;
            let rad = sp.radius;
            let n = section_count(ang);
            let mut arr = CircSectionArray::new(1, n);
            for i in 1..=n {
                let u = u1 + f64::from(i - 1) * (u2 - u1) / f64::from(n - 1);
                arr.change_value(i)
                    .unwrap()
                    .set_circ(&slib::sphere_u_iso(&sp.pos, rad, u), v1, v2);
            }
            sec = Some(arr);
        }
        // `default: break;` (`cxx:550-551`): `sec` stays null.
        _ => {}
    }

    // `SD->SetSimul(sec);` (`cxx:553`): `sec` is null on the `default` arm, and
    // `SetSimul` stores the handle unconditionally.
    match sec {
        Some(arr) => sd.set_simul(Arc::new(arr)),
        None => sd.reset_simul(),
    }
    Ok(())
}

/// `(int)(36. * ang / M_PI + 1)`, clamped to at least two (`cxx:514-518`,
/// `:536-540`).
fn section_count(ang: f64) -> i32 {
    let n = (36.0 * ang / std::f64::consts::PI + 1.0) as i32;
    if n < 2 {
        2
    } else {
        n
    }
}

/// The `Value()` of one of the four `PCurveOnSurf`s of `SD` (`cxx:479-486`).
fn p_curve_point(sd: &SurfData, on_s: i32, first: bool) -> Result<GpPnt2d, String> {
    let interf = sd.interference(on_s);
    let pc = interf.p_curve_on_surf().ok_or_else(|| {
        format!("SimulKPart: PCurveOnSurf of interference {on_s} is null")
    })?;
    let par = if first {
        interf.first_parameter()
    } else {
        interf.last_parameter()
    };
    Ok(pc.value(par))
}

impl ChFi3dSurfOps for FilSurfOps {
    /// `SimulKPart(SD)` (`cxx:475-554`).
    fn simul_k_part(&mut self, sd: &mut SurfData, ds: &HDataStructure) -> Result<(), String> {
        simul_k_part_impl(sd, ds)
    }

    /// `SimulSurf(...)` on two faces (`cxx:558-787`).
    #[allow(clippy::too_many_arguments)]
    fn simul_surf(
        &mut self,
        data: &mut SurfData,
        guide: &mut Arc<ElSpine>,
        spine: Option<&SpineHandle>,
        choix: i32,
        s1: &SurfaceAdaptor,
        i1: &std::rc::Rc<std::cell::RefCell<TopolTool>>,
        s2: &SurfaceAdaptor,
        i2: &std::rc::Rc<std::cell::RefCell<TopolTool>>,
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
        ctx: &ChFi3dCtx<'_>,
    ) -> bool {
        simul_surf_face_face(
            self.shape,
            data,
            guide,
            spine,
            choix,
            s1,
            i1,
            s2,
            i2,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            rec_on_s1,
            rec_on_s2,
            soldep,
            intf,
            intl,
            ctx,
        )
        // `throw Standard_ConstructionError("SimulSurf : this is not the spine
        // of the fillet")` (`cxx:578-582`) is the `Err`, which the caller sees
        // as the same "did not simulate".
        .unwrap_or(false)
    }

    /// `SimulSurf(...)` on `S1` (`cxx:788-1017`).
    #[allow(clippy::too_many_arguments)]
    fn simul_surf_on_s1(
        &mut self,
        data: &mut SurfData,
        guide: &Arc<ElSpine>,
        spine: Option<&SpineHandle>,
        choix: i32,
        s1: &SurfaceAdaptor,
        i1: &std::rc::Rc<std::cell::RefCell<TopolTool>>,
        pc1: Option<&crate::brep_adaptor_curve2d::Curve2dAdaptor>,
        sref1: &SurfaceAdaptor,
        pcref1: &crate::brep_adaptor_curve2d::Curve2dAdaptor,
        decroch1: &mut bool,
        s2: &SurfaceAdaptor,
        i2: &std::rc::Rc<std::cell::RefCell<TopolTool>>,
        or2: crate::abs::Orientation,
        _fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        approx: bool,
        forward: bool,
        rec_p: bool,
        rec_s: bool,
        rec_rst: bool,
        soldep: &MathVector,
        ctx: &ChFi3dCtx<'_>,
    ) -> Result<(), String> {
        // `Fleche` is commented out in the OCCT parameter list (`cxx:801`):
        // the body takes its step from `SimulParams` only.
        simul_surf_on_s1(
            self.shape,
            data,
            guide,
            spine,
            choix,
            s1,
            i1,
            pc1,
            sref1,
            pcref1,
            decroch1,
            s2,
            i2,
            or2,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            rec_p,
            rec_s,
            rec_rst,
            soldep,
            ctx,
        )
    }

    /// `SimulSurf(...)` on `S2` (`cxx:1020-1248`).
    #[allow(clippy::too_many_arguments)]
    fn simul_surf_on_s2(
        &mut self,
        data: &mut SurfData,
        guide: &Arc<ElSpine>,
        spine: Option<&SpineHandle>,
        choix: i32,
        s1: &SurfaceAdaptor,
        i1: &std::rc::Rc<std::cell::RefCell<TopolTool>>,
        or1: crate::abs::Orientation,
        s2: &SurfaceAdaptor,
        i2: &std::rc::Rc<std::cell::RefCell<TopolTool>>,
        pc2: Option<&crate::brep_adaptor_curve2d::Curve2dAdaptor>,
        sref2: &SurfaceAdaptor,
        pcref2: &crate::brep_adaptor_curve2d::Curve2dAdaptor,
        decroch2: &mut bool,
        _fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        approx: bool,
        forward: bool,
        rec_p: bool,
        rec_s: bool,
        rec_rst: bool,
        soldep: &MathVector,
        ctx: &ChFi3dCtx<'_>,
    ) -> Result<(), String> {
        simul_surf_on_s2(
            self.shape,
            data,
            guide,
            spine,
            choix,
            s1,
            i1,
            or1,
            s2,
            i2,
            pc2,
            sref2,
            pcref2,
            decroch2,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            rec_p,
            rec_s,
            rec_rst,
            soldep,
            ctx,
        )
    }

    /// `PerformFirstSection(...)` (`cxx:1500-1537`).
    #[allow(clippy::too_many_arguments)]
    fn perform_first_section(
        &mut self,
        spine: &SpineHandle,
        guide: &Arc<ElSpine>,
        choix: i32,
        s1: &SurfaceAdaptor,
        s2: &SurfaceAdaptor,
        i1: &std::rc::Rc<std::cell::RefCell<TopolTool>>,
        i2: &std::rc::Rc<std::cell::RefCell<TopolTool>>,
        par: f64,
        sol_dep: &mut MathVector,
        pos1: &mut State,
        pos2: &mut State,
        ctx: &ChFi3dCtx<'_>,
    ) -> bool {
        perform_first_section(
            self.shape, spine, guide, choix, s1, s2, i1, i2, par, sol_dep, pos1, pos2, ctx,
        )
        // The `Standard_ConstructionError` of a non-fillet spine
        // (`cxx:1512-1517`) reaches `StartSol` here as "section not found".
        .unwrap_or(false)
    }

    /// `PerformSurf(...)` on two faces (`cxx:1538-1717`).
    #[allow(clippy::too_many_arguments)]
    fn perform_surf(
        &mut self,
        base: &mut crate::chfi3d::builder::ChFi3dBuilder,
        data: &mut Vec<Arc<SurfData>>,
        guide: &mut Arc<ElSpine>,
        spine: Option<&SpineHandle>,
        choix: i32,
        s1: &SurfaceAdaptor,
        i1: &std::rc::Rc<std::cell::RefCell<TopolTool>>,
        s2: &SurfaceAdaptor,
        i2: &std::rc::Rc<std::cell::RefCell<TopolTool>>,
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
        // The two `throw` arms (`cxx:1565-1568`, `:1641-1644`) and the
        // `Standard_NullObject` of `S1->Face()` reach the caller as the same
        // "did not compute".
        perform_surf_face_face(
            base,
            self.shape,
            data,
            guide,
            spine,
            choix,
            s1,
            i1,
            s2,
            i2,
            max_step,
            fleche,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            rec_on_s1,
            rec_on_s2,
            soldep,
            intf,
            intl,
        )
    }

    /// `PerformSurf(...)` on `S1` (`cxx:1721-1890`).
    #[allow(clippy::too_many_arguments)]
    fn perform_surf_on_s1(
        &mut self,
        base: &mut crate::chfi3d::builder::ChFi3dBuilder,
        data: &mut Vec<Arc<SurfData>>,
        guide: &Arc<ElSpine>,
        spine: Option<&SpineHandle>,
        choix: i32,
        s1: &SurfaceAdaptor,
        i1: &std::rc::Rc<std::cell::RefCell<TopolTool>>,
        pc1: Option<&crate::brep_adaptor_curve2d::Curve2dAdaptor>,
        sref1: &SurfaceAdaptor,
        pcref1: &crate::brep_adaptor_curve2d::Curve2dAdaptor,
        decroch1: &mut bool,
        s2: &SurfaceAdaptor,
        i2: &std::rc::Rc<std::cell::RefCell<TopolTool>>,
        or2: crate::abs::Orientation,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        approx: bool,
        forward: bool,
        rec_p: bool,
        rec_s: bool,
        rec_rst: bool,
        soldep: &MathVector,
    ) -> Result<(), String> {
        // The two `throw`s of the body are the trait's `false`; `CallPerformSurf`
        // reacts to it exactly as it reacts to the C++ exception.
        let done = perform_surf_on_s1(
            base,
            self.shape,
            data,
            guide,
            spine,
            choix,
            s1,
            i1,
            pc1,
            sref1,
            pcref1,
            decroch1,
            s2,
            i2,
            or2,
            max_step,
            fleche,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            rec_p,
            rec_s,
            rec_rst,
            soldep,
        );
        if done {
            Ok(())
        } else {
            Err("PerformSurf : Failed processing!".to_string())
        }
    }

    /// `PerformSurf(...)` on `S2` (`cxx:1894-2063`).
    #[allow(clippy::too_many_arguments)]
    fn perform_surf_on_s2(
        &mut self,
        base: &mut crate::chfi3d::builder::ChFi3dBuilder,
        data: &mut Vec<Arc<SurfData>>,
        guide: &Arc<ElSpine>,
        spine: Option<&SpineHandle>,
        choix: i32,
        s1: &SurfaceAdaptor,
        i1: &std::rc::Rc<std::cell::RefCell<TopolTool>>,
        or1: crate::abs::Orientation,
        s2: &SurfaceAdaptor,
        i2: &std::rc::Rc<std::cell::RefCell<TopolTool>>,
        pc2: Option<&crate::brep_adaptor_curve2d::Curve2dAdaptor>,
        sref2: &SurfaceAdaptor,
        pcref2: &crate::brep_adaptor_curve2d::Curve2dAdaptor,
        decroch2: &mut bool,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        approx: bool,
        forward: bool,
        rec_p: bool,
        rec_s: bool,
        rec_rst: bool,
        soldep: &MathVector,
    ) -> Result<(), String> {
        let done = perform_surf_on_s2(
            base,
            self.shape,
            data,
            guide,
            spine,
            choix,
            s1,
            i1,
            or1,
            s2,
            i2,
            pc2,
            sref2,
            pcref2,
            decroch2,
            max_step,
            fleche,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            rec_p,
            rec_s,
            rec_rst,
            soldep,
        );
        if done {
            Ok(())
        } else {
            Err("PerformSurf : Failed processing!".to_string())
        }
    }

    /// `PerformSurf(...)` edge / edge (`cxx:2065-2270`).
    #[allow(clippy::too_many_arguments)]
    fn perform_surf_contact(
        &mut self,
        base: &mut crate::chfi3d::builder::ChFi3dBuilder,
        data: &mut Vec<Arc<SurfData>>,
        guide: &Arc<ElSpine>,
        spine: Option<&SpineHandle>,
        choix: i32,
        s1: &SurfaceAdaptor,
        i1: &std::rc::Rc<std::cell::RefCell<TopolTool>>,
        pc1: Option<&crate::brep_adaptor_curve2d::Curve2dAdaptor>,
        sref1: &SurfaceAdaptor,
        pcref1: &crate::brep_adaptor_curve2d::Curve2dAdaptor,
        decroch1: &mut bool,
        or1: crate::abs::Orientation,
        s2: &SurfaceAdaptor,
        i2: &std::rc::Rc<std::cell::RefCell<TopolTool>>,
        pc2: Option<&crate::brep_adaptor_curve2d::Curve2dAdaptor>,
        sref2: &SurfaceAdaptor,
        pcref2: &crate::brep_adaptor_curve2d::Curve2dAdaptor,
        decroch2: &mut bool,
        or2: crate::abs::Orientation,
        max_step: f64,
        fleche: f64,
        tol_guide: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        approx: bool,
        forward: bool,
        rec_p1: bool,
        rec_rst1: bool,
        rec_p2: bool,
        rec_rst2: bool,
        soldep: &MathVector,
    ) -> Result<(), String> {
        let done = perform_surf_contact(
            base,
            self.shape,
            data,
            guide,
            spine,
            choix,
            s1,
            i1,
            pc1,
            sref1,
            pcref1,
            decroch1,
            or1,
            s2,
            i2,
            pc2,
            sref2,
            pcref2,
            decroch2,
            or2,
            max_step,
            fleche,
            tol_guide,
            first,
            last,
            inside,
            approx,
            forward,
            rec_p1,
            rec_rst1,
            rec_p2,
            rec_rst2,
            soldep,
        );
        if done {
            Ok(())
        } else {
            Err("PerformSurf : Failed processing!".to_string())
        }
    }
}

/// The corner-side overrides of `ChFi3d_FilBuilder` (`ChFi3dCornerOps`).
///
/// `ExtentOneCorner` / `ExtentTwoCorner` / `ExtentThreeCorner`
/// (`ChFi3d_FilBuilder.cxx:2440-2591`), `SetRegul` (`cxx:2595-2634`),
/// `PerformTwoCorner` (`ChFi3d_FilBuilder_C2.cxx:143-1153`) and
/// `PerformThreeCorner` (`ChFi3d_FilBuilder_C3.cxx:241-1331`) are ported.
///
/// `PerformOneCorner` / `PerformMoreSurfdata` are non-virtual base members in
/// OCCT (`ChFi3d_Builder.hxx:780`, `:784`), so the overrides forward to the
/// inherent bodies; `PerformOneCorner` also runs the
/// `PerformIntersectionAtEnd` exit it takes at its two early returns
/// (`ChFi3d_Builder_C1.cxx:821-825`, `:910-914`).
impl ChFi3dCornerOps for FilSurfOps {
    /// `SetRegul()` (`cxx:2595-2634`) — see [`super::set_regul`].
    fn set_regul(&mut self, base: &ChFi3dBuilder) -> Result<(), String> {
        super::set_regul::set_regul(base)
    }

    /// `PerformTwoCorner(Index)` (`ChFi3d_FilBuilder_C2.cxx:143-1153`) —
    /// see [`super::perform_two_corner::perform_two_corner`].
    fn perform_two_corner(
        &mut self,
        base: &mut ChFi3dBuilder,
        index: usize,
    ) -> Result<(), String> {
        super::perform_two_corner::perform_two_corner(base, self.shape, index)
    }

    /// `PerformThreeCorner(Index)` (`ChFi3d_FilBuilder_C3.cxx:241-1331`) —
    /// see [`super::perform_three_corner::perform_three_corner`].
    fn perform_three_corner(
        &mut self,
        base: &mut ChFi3dBuilder,
        index: usize,
    ) -> Result<(), String> {
        super::perform_three_corner::perform_three_corner(base, self.shape, index)
    }

    /// `ExtentOneCorner(V, S)` (`cxx:2440-2469`) — see
    /// [`super::extent_corner::extent_one_corner`].
    fn extent_one_corner(
        &mut self,
        _base: &ChFi3dBuilder,
        v: &TopoShape,
        s: StripeHandle,
    ) -> Result<(), String> {
        super::extent_corner::extent_one_corner(v, &s)
    }

    /// `ExtentTwoCorner(V, LS)` (`cxx:2474-2543`).
    fn extent_two_corner(
        &mut self,
        _base: &ChFi3dBuilder,
        v: &TopoShape,
        ls: &[StripeHandle],
    ) -> Result<(), String> {
        super::extent_corner::extent_two_corner(v, ls)
    }

    /// `ExtentThreeCorner(V, LS)` (`cxx:2548-2594`).
    fn extent_three_corner(
        &mut self,
        _base: &ChFi3dBuilder,
        v: &TopoShape,
        ls: &[StripeHandle],
    ) -> Result<(), String> {
        super::extent_corner::extent_three_corner(v, ls)
    }

    /// `PerformOneCorner(Index, PrepareOnSame)` (`ChFi3d_Builder_C1.cxx:611-1647`).
    fn perform_one_corner(
        &mut self,
        base: &mut ChFi3dBuilder,
        index: usize,
        prepare_on_same: bool,
    ) -> Result<(), String> {
        match base.perform_one_corner(index, prepare_on_same)? {
            // `cxx:1647`.
            OneCornerExit::Done => Ok(()),
            // `PerformIntersectionAtEnd(Index); return;` (`cxx:821-825`,
            // `:910-914`).
            OneCornerExit::IntersectionAtEnd => base.dispatch_end_intersection_exit(index),
        }
    }

    /// `PerformMoreSurfdata(Index)` (`ChFi3d_Builder_C1.cxx:3771-4531`).
    fn perform_more_surfdata(
        &mut self,
        base: &mut ChFi3dBuilder,
        index: usize,
    ) -> Result<(), String> {
        base.perform_more_surfdata(index)
    }
}