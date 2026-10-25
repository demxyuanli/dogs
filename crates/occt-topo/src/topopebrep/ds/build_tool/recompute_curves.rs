//! `TopOpeBRepDS_BuildTool` — the curve / pcurve recomputation family.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepDS/TopOpeBRepDS_BuildTool.cxx`
//! (OCCT 8.0.0): `FUN_getUV` `:557-573`, `FUN_UisoLineOnSphe` `:62-83`
//! (ported in the parent module), `FUN_makeUisoLineOnSphe` `:602-651`,
//! `ComputePCurves` `:654-789`, `PutPCurves` `:793-813`, `RecomputeCurves`
//! `:816-890`.
//!
//! UNPORTED: nothing.

use std::sync::Arc;

use occt_core::elib::clib::adjust_periodic;
use occt_core::gp::{GpAx2d, GpDir2d, GpPnt2d, GpVec2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom::trimmed::GeomTrimmedCurveBasis;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::tool::curve_tool::CurveTool;
use crate::topopebrep::tool::shape_tool::resolution3d_face;

use super::super::curve::Curve as CurveData;
use super::BuildTool;

/// `gp::DY2d()` (`gp.hxx`).
const DY2D: GpDir2d = GpDir2d { x: 0.0, y: 1.0 };
/// `gp::NY2d()` (`gp.hxx`).
const NY2D: GpDir2d = GpDir2d { x: 0.0, y: -1.0 };

impl BuildTool {
    /// `ComputePCurves(C, E, newC, comppc1, comppc2, compc3d)`
    /// (`cxx:654-789`): project `C`'s 3D curve on the two faces of `newC`,
    /// re-range the edge vertices on the new curve and store the pcurves back
    /// into `newC` through its two `SurfaceCurveInterference`s.
    ///
    /// `newC` must already be the `HDS->ChangeCurve(inewC)` of
    /// [`Self::recompute_curves`]: the store goes through `newC`'s SCI
    /// handles, which are the ones `HDataStructure::make_curve` registered.
    pub fn compute_pcurves(
        &self,
        c: &CurveData,
        e: &TopoShape,
        new_c: &mut CurveData,
        comppc1: bool,
        comppc2: bool,
        compc3d: bool,
    ) {
        // cxx:661-662: `TopoDS::Face(newC.Shape1())` / `(Shape2())`.
        let Some(f1) = new_c.shape1().and_then(|s| Face::wrap(s.clone())) else {
            return;
        };
        let Some(f2) = new_c.shape2().and_then(|s| Face::wrap(s.clone())) else {
            return;
        };

        let c3d = c.curve().cloned();

        // cxx:665-670: bounding vertices `Vmin`, `Vmax` of `E` and their
        // parameters. `GetOrientedEdgeVertices` leaves `parmin` / `parmax` at
        // their `0.` initialisation when a vertex is null.
        let edge = Edge::wrap(e.clone()).expect("TopOpeBRepDS_BuildTool::ComputePCurves: not an edge");
        let (vmin, _vmax, params) = self.get_oriented_edge_vertices(&edge).unwrap_or((
            TopoShape::new(ShapeType::Vertex),
            TopoShape::new(ShapeType::Vertex),
            None,
        ));
        let (mut parmin, mut parmax) = params.unwrap_or((0.0, 0.0));

        // cxx:672-700: a null 3D curve ("sharp patches"). Both
        // `MakePCurveOnFace` calls in this arm hand the null handle to
        // `GeomAdaptor_Curve::Load`, which OCCT dereferences; the port cannot
        // project without a curve, so `PC1new` / `PC2new` stay null and only
        // the tolerance is rebuilt from `Resolution3d` of the two faces.
        let Some(c3d) = c3d else {
            let r1 = resolution3d_face(&f1.0, CONFUSION).unwrap_or(0.0);
            let r2 = resolution3d_face(&f2.0, CONFUSION).unwrap_or(0.0);
            let tol = CONFUSION.max(r1).max(r2);
            new_c.set_tolerance(tol);
            return;
        };

        let mut c3dnew: Arc<dyn Curve> = c3d.clone();

        // cxx:705-727: periodize `parmin` / `parmax` on a periodic curve and
        // trim the 3D curve to that window when `CompC3D` is on.
        if c3d.is_periodic() {
            let period = c3d.last_parameter() - c3d.first_parameter();
            let (f, l) = if vmin.orientation() == Orientation::Forward {
                (parmin, parmax)
            } else {
                (parmax, parmin)
            };
            parmin = f;
            parmax = l;
            adjust_periodic(f, f + period, PCONFUSION, &mut parmin, &mut parmax);
            if compc3d {
                c3dnew = Arc::new(GeomTrimmedCurveBasis::new(c3d.clone(), parmin, parmax));
            }
        }

        // cxx:730-741
        let mut tolreached2d1 = c.tolerance();
        let mut tolreached2d2 = c.tolerance();
        let mut pc1new = None;
        let mut pc2new = None;
        if comppc1 {
            if let Some(pc) = project_on_face(&f1, c3dnew.clone(), &mut tolreached2d1) {
                pc1new = Some(pc);
            }
        }
        if comppc2 {
            if let Some(pc) = project_on_face(&f2, c3dnew.clone(), &mut tolreached2d2) {
                pc2new = Some(pc);
            }
        }

        // cxx:743-755
        let (newtol, newparmin, newparmax) = self.update_edge_curve_tol(e, Some(c3dnew.clone()));

        // cxx:760-772: a `u`-isoline of the sphere is rebuilt from the
        // projected 3D samples so that it stays an exact 2D line.
        if let Some(pc) = &pc1new {
            if super::fun_u_iso_line_on_sphere(&f1, pc) {
                if let Some(newpc) = make_uiso_line_on_sphere(&f1, c3dnew.as_ref(), pc, newtol) {
                    pc1new = Some(newpc);
                }
            }
        }
        if let Some(pc) = &pc2new {
            if super::fun_u_iso_line_on_sphere(&f2, pc) {
                if let Some(newpc) = make_uiso_line_on_sphere(&f2, c3dnew.as_ref(), pc, newtol) {
                    pc2new = Some(newpc);
                }
            }
        }

        // cxx:774-788
        new_c.set_curve(Some(c3dnew), newtol);
        new_c.set_range(newparmin, newparmax);
        if let Some(pc) = pc1new {
            new_c.set_curve1(Some(pc));
        }
        if let Some(pc) = pc2new {
            new_c.set_curve2(Some(pc));
        }
    }

    /// `PutPCurves(newC, E, comppc1, comppc2)` (`cxx:793-813`): install
    /// `newC`'s two pcurves on `E` through [`BuildTool::p_curve`].
    pub fn put_pcurves(&self, new_c: &CurveData, e: &TopoShape, comppc1: bool, comppc2: bool) {
        if comppc1 {
            if let (Some(f1), Some(pc1)) = (new_c.shape1(), new_c.curve1()) {
                self.p_curve(f1, e, &pc1);
            }
        }
        if comppc2 {
            if let (Some(f2), Some(pc2)) = (new_c.shape2(), new_c.curve2()) {
                self.p_curve(f2, e, &pc2);
            }
        }
    }

    /// `RecomputeCurves(C, oldE, E, inewC, HDS)` (`cxx:816-890`).
    ///
    /// `C` is taken by reference while `HDS` is borrowed mutably, so the caller
    /// passes a clone of `HDS->Curve(oldE)` — OCCT reads `C` and writes `HDS`
    /// through the same object only in `MakeCurve`, which the port keeps
    /// non-aliasing.
    pub fn recompute_curves(
        &self,
        c: &CurveData,
        e: &TopoShape,
        inew_c: &mut i32,
        hds: &mut HDataStructure,
    ) {
        // cxx:825-831
        let gt = self.geom_tool;
        let compc3d = gt.comp_c3d();
        let comppc1 = gt.comp_pc1();
        let comppc2 = gt.comp_pc2();
        let comppc = comppc1 || comppc2;
        let iswalk = c.is_walk();
        let approx = self.approximation();
        let c3d = c.curve().cloned();

        // cxx:834-843: `Standard_ProgramError` guards; the port returns.
        if comppc1 && c.shape1().is_none() {
            return;
        }
        if comppc2 && c.shape2().is_none() {
            return;
        }
        // cxx:844-851: `TopExp::Vertices(E, Vmin, Vmax)`, both non-null.
        let verts = crate::topo_tools::explorer_shapes(e, ShapeType::Vertex);
        if verts.first().is_none() || verts.last().is_none() {
            return;
        }

        if iswalk && approx {
            // cxx:851-861
            if compc3d && c3d.is_none() {
                return;
            }
            self.approx_curves(c, e, inew_c, hds);
            let new_c = hds.change_curve(*inew_c);
            self.put_pcurves(new_c, e, comppc1, comppc2);
            return;
        }

        if !comppc {
            return;
        }

        // cxx:871-889
        let mut new_c1 = CurveData::new();
        *inew_c = hds
            .make_curve(c, &mut new_c1)
            .expect("TopOpeBRepDS_HDataStructure::MakeCurve");
        if iswalk && !approx {
            // cxx:876-882: keep the pcurves of the walking curve.
            if compc3d && c3d.is_none() {
                return;
            }
            let pc1 = c.curve1();
            let pc2 = c.curve2();
            let new_c = hds.change_curve(*inew_c);
            new_c.set_curve1(pc1);
            new_c.set_curve2(pc2);
        } else {
            // cxx:884-887
            self.compute_pcurves(
                c,
                e,
                hds.change_curve(*inew_c),
                comppc1,
                comppc2,
                compc3d,
            );
        }
        // cxx:888
        self.put_pcurves(hds.change_curve(*inew_c), e, comppc1, comppc2);
    }
}

/// `TopOpeBRepTool_CurveTool::MakePCurveOnFace(F, C3D, TolReached2d)`
/// (`CurveTool.cxx:1040-1156`) bound to a face: the surface, the face's
/// `BRepTools::UVBounds` and the whole curve (`first = last = 0`, so the
/// `trim` arm of `cxx:1046-1052` is skipped).
fn project_on_face(
    f: &Face,
    c3d: Arc<dyn Curve>,
    tol_reached_2d: &mut f64,
) -> Option<Arc<dyn Curve2d>> {
    let reg = GeometryRegistry::global();
    let surface = reg.face_surface(&f.0)?;
    let bounds = BRepTool::uv_bounds(f);
    let (pc, tol) = CurveTool::make_pcurve_on_face(surface, c3d, 0.0, 0.0, Some(bounds));
    *tol_reached_2d = tol;
    pc
}

/// `FUN_getUV(surf, C3D, par3d, u0, v0)` (`cxx:557-573`): the `(u, v)` of the
/// point of `C3D` at `par3d` on `surf`, through `GeomAPI_ProjectPointOnSurf`
/// and its `LowerDistanceParameters`. `None` where OCCT answers `false`
/// (`NbPoints() < 1`).
fn fun_get_uv(surf: &dyn Surface, c3d: &dyn Curve, par3d: f64) -> Option<(f64, f64)> {
    let p3d = c3d.d0(par3d);
    let p = occt_geom::geom_api::project_point_on_surface(surf, &p3d, CONFUSION)?;
    Some((p.u, p.v))
}

/// `FUN_makeUisoLineOnSphe(F, C3D, PCnew, tol3d)` (`cxx:602-651`): rebuild the
/// pcurve of a sphere's `u`-isoline as the exact 2D line through the projected
/// sample at `0.243234` of the 3D range, growing in `+V` or `-V` depending on
/// the sample at `0.543345`.
///
/// OCCT mutates `PCnew` in place through `Geom2d_Line::SetLin2d` and answers
/// `false` when the two samples disagree in `U`; the port answers the new line
/// and `None` for that `false`.
fn make_uiso_line_on_sphere(
    f: &Face,
    c3d: &dyn Curve,
    pcnew: &Arc<dyn Curve2d>,
    tol3d: f64,
) -> Option<Arc<dyn Curve2d>> {
    // cxx:608-614
    let p3df = c3d.first_parameter();
    let p3dl = c3d.last_parameter();
    const DELTAINF: f64 = 0.243234;
    const DELTASUP: f64 = 0.543345;
    let par3dinf = (1.0 - DELTAINF) * p3df + DELTAINF * p3dl;
    let par3dsup = (1.0 - DELTASUP) * p3df + DELTASUP * p3dl;

    let surf = GeometryRegistry::global().face_surface(&f.0)?;
    let (uinf, vinf) = fun_get_uv(surf.as_ref(), c3d, par3dinf)?;
    let (usup, vsup) = fun_get_uv(surf.as_ref(), c3d, par3dsup)?;

    // `Precision::Parametric(tol3d)` = `tol3d * 0.01` (`Precision.hxx:328`).
    let tol = tol3d * 0.01;
    if (uinf - usup).abs() > tol {
        return None;
    }

    // cxx:632-640
    let vdir = if vsup - vinf > -tol { DY2D } else { NY2D };
    let mut origin = GpPnt2d::new(uinf, vinf);
    origin.translate_vec(&GpVec2d::new(vdir.x, vdir.y).multiplied_scalar(p3df - par3dinf));

    // cxx:643-649: `BASISCURVE2D(PCnew)` then `down_cast<Geom2d_Line>`. A null
    // down-cast would be dereferenced by `SetLin2d`; the port refuses instead.
    let basis: &dyn Curve2d = pcnew.trimmed_basis().unwrap_or(pcnew.as_ref());
    basis.gp_lin2d()?;
    let line = Geom2dLine::new(GpAx2d::new(origin, vdir));
    Some(Arc::new(line))
}
