//! `ChFi3d_Builder::CompleteData` — the four overloads that turn a filling
//! result into a `ChFiDS_SurfData`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_6.cxx:393-624`
//! (OCCT 8.0.0). The four C++ overloads become four methods because Rust has
//! no overloading; each keeps the C++ control flow:
//!
//! - [`ChFi3dBuilder::complete_data_geom`] (`cxx:393-505`), the surface /
//!   iso-curve overload,
//! - [`ChFi3dBuilder::complete_data_app_func`] (`cxx:514-550`), `Blend_Function`
//!   through `BRepBlend_AppFunc` + `BRepBlend_AppSurface`,
//! - [`ChFi3dBuilder::complete_data_rst`] (`cxx:557-586`), the
//!   `Blend_SurfRstFunction` overload,
//! - [`ChFi3dBuilder::complete_data_rst_rst`] (`cxx:594-624`), the
//!   `Blend_RstRstFunction` overload, whose `StoreData` call has no `Reversed`
//!   argument (it is `false`).
//!
//! Divergences that are not heuristics:
//! - `occ::handle<T>&` becomes `&mut T` (the `SurfData`) or an owned
//!   `Option<Arc<..>>`; `handle<Geom_Surface>` is the port's
//!   `Arc<dyn GeomSurface>`,
//! - the `Blend_Function&` argument of the second overload becomes
//!   `Rc<RefCell<F>>`, the ownership the port's `BRepBlend_AppFunc` needs; the
//!   third overload takes `BlendSurfRstFunction` the same way,
//! - the null `Adaptor3d_Surface` handle `StoreData` accepts is `Option`, as in
//!   [`store_data`],
//! - `ChFi3d_SameParameter` is called with the adaptor OCCT builds from the
//!   iso curve; the port's [`same_parameter`] reads the curve parameters
//!   directly, which are the adaptor parameters for an untrimmed
//!   `GeomAdaptor_Curve(Crv3d)`.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d};
use occt_geom::curve::Curve;
use occt_geom::surface::Surface as GeomSurface;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::Orientation;
use crate::blend::{BlendFunction, BlendRstRstFunction, BlendSurfRstFunction};
use crate::brep_blend::app_func::AppFuncRoot;
use crate::brep_blend::app_surface::AppSurface;
use crate::brep_blend::line::Line;
use crate::chfi3d::builder_helpers::same_parameter;
use crate::chfi_ds::surf_data::SurfData;
use crate::topopebrep::ds::{Curve as DsCurve, Surface as DsSurface};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::{store_data, ChFi3dBuilder};

/// `BRepBlend_Line::Point(1)` / `Point(NbPoints())` parameters (`cxx:531-532`,
/// `:568-569`); OCCT raises `Standard_OutOfRange` on an empty line.
fn line_bounds(lin: &Line) -> Result<(f64, f64), String> {
    let out_of_range = || "BRepBlend_Line::Point: Standard_OutOfRange".to_string();
    let first = lin.point(1).ok_or_else(out_of_range)?.parameter();
    let last = lin.point(lin.nb_points()).ok_or_else(out_of_range)?.parameter();
    Ok((first, last))
}

impl ChFi3dBuilder {
    /// `CompleteData(Data, Surfcoin, S1, PC1, S2, PC2, Or, On1, Gd1, Gd2, Gf1,
    /// Gf2)` (`ChFi3d_Builder_6.cxx:393-505`).
    ///
    /// Calculates the curves and the common points of a `ChFiDS_SurfData` from
    /// the surface `Surfcoin` and the p-curves `PC1` / `PC2` that the filling
    /// produced on `S1` / `S2`.
    #[allow(clippy::too_many_arguments)]
    pub fn complete_data_geom(
        &mut self,
        data: &mut SurfData,
        surfcoin: &Arc<dyn GeomSurface>,
        s1: &SurfaceAdaptor,
        pc1: Option<&Arc<dyn Curve2d>>,
        s2: &SurfaceAdaptor,
        pc2: Option<&Arc<dyn Curve2d>>,
        or: Orientation,
        on1: bool,
        gd1: bool,
        gd2: bool,
        gf1: bool,
        gf2: bool,
    ) -> Result<bool, String> {
        // `Data->ChangeSurf(DStr.AddSurface(...))` (`cxx:408-409`).
        let surf_index = self
            .my_ds
            .change_ds()
            .add_surface(&DsSurface::with_tolerance(Arc::clone(surfcoin), self.tolapp3d));
        data.set_surf(surf_index);

        let (u_first, u_last) = surfcoin.u_range();
        let (v_first, v_last) = surfcoin.v_range();
        if !gd1 {
            data.vertex_mut(true, 1)
                .set_point(&surfcoin.value(u_first, v_first));
        }
        if !gd2 {
            data.vertex_mut(true, 2)
                .set_point(&surfcoin.value(u_first, v_last));
        }
        if !gf1 {
            data.vertex_mut(false, 1)
                .set_point(&surfcoin.value(u_last, v_first));
        }
        if !gf2 {
            data.vertex_mut(false, 2)
                .set_point(&surfcoin.value(u_last, v_last));
        }

        // Calculate curves side S1 (`cxx:430-468`).
        let mut crv3d1: Option<Arc<dyn Curve>> = None;
        if pc1.is_some() {
            crv3d1 = surfcoin.v_iso_curve(v_first);
        }
        let pd1 = GpPnt2d::new(u_first, v_first);
        let pf1 = GpPnt2d::new(u_last, v_first);
        let dir1 = GpDir2d::new(pf1.x() - pd1.x(), pf1.y() - pd1.y())
            .map_err(|_| "gp_Dir2d: null vector".to_string())?;
        let mut p_curve_on_surf = Arc::new(Geom2dLine::from_pnt_dir(pd1, dir1)) as Arc<dyn Curve2d>;
        let mut tra1 = Orientation::Forward;
        let mut orsurf = or;
        let w = 0.5 * (u_first + u_last);
        let mut tolreached = 1.0e-5;
        // `c2dtrim` is declared before the S1 block and reused by the S2 one
        // (`cxx:447`, `:482-483`): when `PC2` is null the S1 trimmed curve stays
        // in the S2 interference, exactly as in OCCT.
        let mut c2dtrim: Option<Arc<dyn Curve2d>> = None;
        if let Some(pc1) = pc1 {
            let crv3d1 = crv3d1
                .as_ref()
                .ok_or_else(|| "Geom_Surface::VIso failed".to_string())?;
            c2dtrim = Some(Arc::new(Geom2dTrimmedCurve::new(
                Arc::clone(pc1),
                u_first,
                u_last,
            )));
            let s1_surface = s1
                .surface()
                .ok_or_else(|| "BRepAdaptor_Surface: no surface".to_string())?;
            let (_, tol) =
                same_parameter(crv3d1.as_ref(), &mut c2dtrim, s1_surface.as_ref(), self.tolapp3d);
            tolreached = tol;
            let puv = c2dtrim
                .as_ref()
                .ok_or_else(|| "Geom2d_TrimmedCurve: null".to_string())?
                .value(w);
            let (_, du, dv) = s1
                .d1(puv.x(), puv.y())
                .ok_or_else(|| "Adaptor3d_Surface::D1 failed".to_string())?;
            let nf = du.crossed(&dv);
            let (_, du_s, dv_s) = surfcoin.d1(w, v_first);
            let ns = du_s.crossed(&dv_s);
            if nf.dot(&ns) > 0.0 {
                tra1 = Orientation::Reversed;
            } else if on1 {
                orsurf = orsurf.reversed();
            }
        }
        let index1_of_curve = self
            .my_ds
            .change_ds()
            .add_curve(&mut DsCurve::with_curve(crv3d1.clone(), tolreached, false));
        {
            let fint1 = data.interference_mut(1);
            fint1.set_first_parameter(u_first);
            fint1.set_last_parameter(u_last);
            fint1.set_interference(index1_of_curve, tra1, c2dtrim.clone(), Some(p_curve_on_surf));
        }

        // Calculate curves side S2 (`cxx:470-504`).
        let mut crv3d2: Option<Arc<dyn Curve>> = None;
        if pc2.is_some() {
            crv3d2 = surfcoin.v_iso_curve(v_last);
        }
        let pd2 = GpPnt2d::new(u_first, v_last);
        let pf2 = GpPnt2d::new(u_last, v_last);
        let dir2 = GpDir2d::new(pf2.x() - pd2.x(), pf2.y() - pd2.y())
            .map_err(|_| "gp_Dir2d: null vector".to_string())?;
        p_curve_on_surf = Arc::new(Geom2dLine::from_pnt_dir(pd2, dir2)) as Arc<dyn Curve2d>;
        let mut tra2 = Orientation::Forward;
        if let Some(pc2) = pc2 {
            let crv3d2 = crv3d2
                .as_ref()
                .ok_or_else(|| "Geom_Surface::VIso failed".to_string())?;
            c2dtrim = Some(Arc::new(Geom2dTrimmedCurve::new(
                Arc::clone(pc2),
                u_first,
                u_last,
            )));
            let s2_surface = s2
                .surface()
                .ok_or_else(|| "BRepAdaptor_Surface: no surface".to_string())?;
            let (_, tol) =
                same_parameter(crv3d2.as_ref(), &mut c2dtrim, s2_surface.as_ref(), self.tolapp3d);
            tolreached = tol;
            let puv = c2dtrim
                .as_ref()
                .ok_or_else(|| "Geom2d_TrimmedCurve: null".to_string())?
                .value(w);
            let (_, du, dv) = s2
                .d1(puv.x(), puv.y())
                .ok_or_else(|| "Adaptor3d_Surface::D1 failed".to_string())?;
            let np = du.crossed(&dv);
            let (_, du_s, dv_s) = surfcoin.d1(w, v_last);
            let ns = du_s.crossed(&dv_s);
            if np.dot(&ns) < 0.0 {
                tra2 = Orientation::Reversed;
                if !on1 {
                    orsurf = orsurf.reversed();
                }
            }
        }
        let index2_of_curve = self
            .my_ds
            .change_ds()
            .add_curve(&mut DsCurve::with_curve(crv3d2, tolreached, false));
        {
            let fint2 = data.interference_mut(2);
            fint2.set_first_parameter(u_first);
            fint2.set_last_parameter(u_last);
            fint2.set_interference(index2_of_curve, tra2, c2dtrim, Some(p_curve_on_surf));
        }
        data.set_orientation(orsurf);
        Ok(true)
    }

    /// `CompleteData(Data, Func, lin, S1, S2, Or1, Gd1, Gd2, Gf1, Gf2,
    /// Reversed)` (`ChFi3d_Builder_6.cxx:514-550`): approximates the filling
    /// line with `BRepBlend_AppSurface` and hands the result to [`store_data`].
    #[allow(clippy::too_many_arguments)]
    pub fn complete_data_app_func<F: BlendFunction + 'static>(
        &mut self,
        data: &mut SurfData,
        func: Rc<RefCell<F>>,
        lin: Rc<RefCell<Line>>,
        s1: &SurfaceAdaptor,
        s2: Option<&SurfaceAdaptor>,
        or1: Orientation,
        gd1: bool,
        gd2: bool,
        gf1: bool,
        gf2: bool,
        reversed: bool,
    ) -> Result<bool, String> {
        // `BRepBlend_AppFunc(lin, Func, tolapp3d, 1.e-5)` (`cxx:528`).
        let the_func = AppFuncRoot::new_surf_surf(Rc::clone(&lin), func, self.tolapp3d, 1.0e-5);
        let (first, last) = {
            let line = lin.borrow();
            line_bounds(&line)?
        };
        // `BRepBlend_AppSurface approx(TheFunc, first, last, tolapp3d, 1.e-5,
        // tolappangle, myConti, Degmax = 20, Segmax = 5000)` (`cxx:530-539`).
        let approx = AppSurface::new(
            Arc::new(the_func),
            first,
            last,
            self.tolapp3d,
            1.0e-5,
            self.tolappangle,
            self.my_conti,
            20,
            5000,
        );
        if !approx.is_done() {
            return Ok(false);
        }
        let line = lin.borrow();
        store_data(
            data,
            &approx,
            &line,
            s1,
            s2,
            or1,
            gd1,
            gd2,
            gf1,
            gf2,
            reversed,
            self.my_ds.change_ds(),
        )
    }

    /// `CompleteData(Data, Func, lin, S1, S2, Or, Reversed)`
    /// (`ChFi3d_Builder_6.cxx:557-586`): the surf / rst overload. The
    /// approximation uses the `BRepBlend_AppSurface` defaults `(C0, 11, 50)`
    /// (`BRepBlend_AppSurface.hxx:36-40`), and the four `Gd` / `Gf` flags are
    /// `false`, as in the source.
    #[allow(clippy::too_many_arguments)]
    pub fn complete_data_rst<F: BlendSurfRstFunction + 'static>(
        &mut self,
        data: &mut SurfData,
        func: Rc<RefCell<F>>,
        lin: Rc<RefCell<Line>>,
        s1: &SurfaceAdaptor,
        s2: Option<&SurfaceAdaptor>,
        or: Orientation,
        reversed: bool,
    ) -> Result<bool, String> {
        // `BRepBlend_AppFuncRst(lin, Func, tolapp3d, 1.e-5)` (`cxx:566-567`).
        let the_func = AppFuncRoot::new_surf_rst(Rc::clone(&lin), func, self.tolapp3d, 1.0e-5);
        let (first, last) = {
            let line = lin.borrow();
            line_bounds(&line)?
        };
        let approx = AppSurface::new(
            Arc::new(the_func),
            first,
            last,
            self.tolapp3d,
            1.0e-5,
            self.tolappangle,
            occt_core::kernel::geomabs::Shape::C0,
            11,
            50,
        );
        if !approx.is_done() {
            return Ok(false);
        }
        let line = lin.borrow();
        store_data(
            data,
            &approx,
            &line,
            s1,
            s2,
            or,
            false,
            false,
            false,
            false,
            reversed,
            self.my_ds.change_ds(),
        )
    }

    /// `CompleteData(Data, Func, lin, S1, S2, Or)` (`ChFi3d_Builder_6.cxx:594-624`):
    /// the rst / rst overload. Like [`Self::complete_data_rst`] it uses the
    /// `BRepBlend_AppSurface` defaults and passes the four `Gd` / `Gf` flags as
    /// `false`; unlike it, the `StoreData` call passes no `Reversed`, so the
    /// argument stays `false`.
    #[allow(clippy::too_many_arguments)]
    pub fn complete_data_rst_rst<F: BlendRstRstFunction + 'static>(
        &mut self,
        data: &mut SurfData,
        func: Rc<RefCell<F>>,
        lin: Rc<RefCell<Line>>,
        s1: &SurfaceAdaptor,
        s2: Option<&SurfaceAdaptor>,
        or: Orientation,
    ) -> Result<bool, String> {
        // `BRepBlend_AppFuncRstRst(lin, Func, tolapp3d, 1.e-5)` (`cxx:603-604`).
        let the_func = AppFuncRoot::new_rst_rst(Rc::clone(&lin), func, self.tolapp3d, 1.0e-5);
        let (first, last) = {
            let line = lin.borrow();
            line_bounds(&line)?
        };
        let approx = AppSurface::new(
            Arc::new(the_func),
            first,
            last,
            self.tolapp3d,
            1.0e-5,
            self.tolappangle,
            occt_core::kernel::geomabs::Shape::C0,
            11,
            50,
        );
        if !approx.is_done() {
            return Ok(false);
        }
        let line = lin.borrow();
        store_data(
            data,
            &approx,
            &line,
            s1,
            s2,
            or,
            false,
            false,
            false,
            false,
            false,
            self.my_ds.change_ds(),
        )
    }
}
