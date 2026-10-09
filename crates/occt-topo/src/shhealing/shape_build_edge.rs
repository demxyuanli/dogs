//! `ShapeBuild_Edge` subset needed by `ShapeFix_ComposeShell::DispatchWires`:
//! `Copy` (`ShapeBuild_Edge.cxx:417-426`), `RemovePCurve` (`:430-443`),
//! `SetRange3d` (`:338-356`), `ReassignPCurve` (`:530-592`) and the
//!
//! `TransformPCurve` (`:596-698`): the `uFact == 1` return, the `Line` arm
//! (`:624-641`), the `Geom2d_BezierCurve` arm (`:643-655`, returns the Bezier
//! with its poles rewritten), the `Geom2d_BSplineCurve` arm (`:687-697`) and
//! the `Geom2d_Conic` arm through `Geom2dConvert_ApproxCurve` (`:659-678`) are
//! ported.

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::gp::{GpAx2, GpAx2d, GpDir2d, GpGTrsf2d, GpLin2d, GpPnt2d, GpTrsf2d, GpVec2d, TrsfForm};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{APPROXIMATION, CONFUSION};
use occt_geom::geom_lib;
use occt_geom::{Curve, Surface};
use occt_geom2d::convert_approx_curve::Geom2dConvertApproxCurve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::Orientation;
use crate::boptools_2d;
use crate::brep_lib_same_range::{check_same_range, same_range, set_rep_ranges};
use crate::meshing::edge_discret::CurveOnSurface;
use crate::shape::{Edge, Face};
use crate::tgeometry::{EdgeGeom, GeometryRegistry};

/// `ShapeBuild_Edge` (`ShapeBuild_Edge.hxx`).
pub struct ShapeBuildEdge;

/// `evaluateMaxSegment` (`BRepLib.cxx:275-295`): the `MaxSegment` handed to
/// `GeomLib::BuildCurve3d` when the caller passes 0.
fn evaluate_max_segment(
    max_segment: i32,
    surface: &Arc<dyn Surface>,
    pcurve: &Arc<dyn Curve2d>,
) -> i32 {
    if max_segment != 0 {
        return max_segment;
    }
    let mut nb_s_knots = 0usize;
    if surface.is_bspline_surface() {
        let nu = surface.bspline_surface_uknots().map_or(0, |k| k.len());
        let nv = surface.bspline_surface_vknots().map_or(0, |k| k.len());
        nb_s_knots = nu.max(nv);
    }
    let nb_c2d_knots = pcurve.bspline_nb_knots().unwrap_or(0);
    (30 + nb_s_knots.max(nb_c2d_knots)) as i32
}

/// `edge.Reversed()` (`TopoDS_Shape::Reversed`): same TShape, flipped
/// orientation.
fn reversed(edge: &Edge) -> Edge {
    let mut e = edge.clone();
    let o = edge.0.orientation();
    e.0.set_orientation(match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        other => other,
    });
    e
}

/// `CountPCurves` (`ShapeBuild_Edge.cxx:510-528`): number of pcurves stored
/// on `edge` for `face` (2 for a seam).
fn count_pcurves(edge: &Edge, face: &Face) -> usize {
    let key = GeometryRegistry::shape_key(&face.0);
    GeometryRegistry::global().edge_pcurves(&edge.0, key).len()
}

impl ShapeBuildEdge {
    /// `Copy(edge, sharepcurves)` (`cxx:417-426`).
    pub fn copy(&self, edge: &Edge, share_pcurves: bool) -> Edge {
        let new_edge = crate::shhealing::copy_replace_vertices(edge);
        if !share_pcurves {
            crate::shhealing::copy_pcurves(&new_edge, edge);
        }
        new_edge
    }

    /// `RemovePCurve(edge, face)` (`cxx:430-443`).
    pub fn remove_pcurve(&self, edge: &Edge, face: &Face) {
        GeometryRegistry::global().remove_pcurves_on_surface(&edge.0, &face.0);
    }

    /// `SetRange3d(edge, first, last)` (`cxx:338-356`).
    pub fn set_range3d(&self, edge: &Edge, first: f64, last: f64) {
        GeometryRegistry::global().set_edge_range(&edge.0, first, last);
    }

    /// `BuildCurve3d(edge)` (`ShapeBuild_Edge.cxx:714-775`), i.e.
    /// `BRepLib::BuildCurve3d(edge, max(1.e-5, BRep_Tool::Tolerance(edge)))`
    /// (`BRepLib.cxx:301-455`) with its defaults `GeomAbs_C1`, `MaxDegree = 14`
    /// and `MaxSegment = 0` (`BRepLib.hxx:90-94`, `ShapeBuild_Edge.cxx:723`).
    ///
    /// The plane arm (`BRepLib.cxx:358-375`) converts the 2D curve analytically
    /// through `GeomLib::To3d`; the general arm (`cxx:377-447`) approximates the
    /// curve-on-surface through `GeomLib::BuildCurve3d`. Both replace the edge's
    /// 3D curve, widen the edge tolerance to at least `max(1.e-5, Tol)` and, for
    /// a single curve-on-surface representation, force `SameParameter`.
    ///
    /// UNPORTED: `BRep_Tool::CurveOnSurface` is indexed over the pcurves of the
    /// edge's CurveOnSurface representations (`BRep_Tool.cxx:488-533`); the
    /// port's `GeometryRegistry::edge_pcurve_reps` resolves the surface of each
    /// representation through `surface_by_ptr`, so a representation whose
    /// surface is not registered is invisible here (OCCT cannot build that
    /// representation at all).
    pub fn build_curve3d(&self, edge: &Edge) -> bool {
        let reg = GeometryRegistry::global();
        // `BRepLib.cxx:320-325`: an edge that already has a 3D curve returns
        // true without touching it.
        if reg.edge_curve(&edge.0).is_some() {
            return true;
        }
        // `ShapeBuild_Edge.cxx:723`.
        let tolerance = reg.edge_tolerance(&edge.0).max(1.0e-5);
        // `BRepLib.cxx:330-333`.
        if !check_same_range(&edge.0, CONFUSION) {
            same_range(&edge.0, tolerance);
        }

        // `BRepLib.cxx:335-357`: search a curve on a plane.
        let reps = reg.edge_pcurve_reps(&edge.0);
        let mut plane: Option<(Arc<dyn Curve2d>, f64, f64, GpAx2)> = None;
        for (surf, pc, f, l) in &reps {
            let basis = surf.rectangular_trimmed_basis().unwrap_or_else(|| surf.clone());
            if let Some(pln) = basis.gp_pln() {
                plane = Some((pc.clone(), *f, *l, pln.position().ax2()));
                break;
            }
        }

        if let Some((pc, f, l, axes)) = plane {
            // `cxx:358-368`: `GeomLib::To3d(axes, PC)`; a null handle is the
            // `Standard_NotImplemented` throw and returns false.
            let Some(c3d) = geom_lib::to_3d(&axes, pc.as_ref()) else {
                return false;
            };
            // `cxx:369-373`: `B.UpdateEdge(AnEdge, C3d, LocalLoc, 0.0e0)` then
            // `BRep_Tool::Range(AnEdge, S, LC, First, Last)` +
            // `B.Range(AnEdge, First, Last)`. `UpdateEdge` cannot lower the edge
            // tolerance, so it is carried over together with the edge flags.
            let mut g = EdgeGeom::new(c3d, f, l);
            g.tolerance = reg.edge_tolerance(&edge.0);
            g.same_parameter = reg.same_parameter(&edge.0);
            g.same_range = reg.same_range(&edge.0);
            g.degenerated = reg.is_degenerated_edge(&edge.0);
            reg.set_edge(&edge.0, g);
            set_rep_ranges(&edge.0, f, l);
        } else {
            // `cxx:377-452`.
            if reg.is_degenerated_edge(&edge.0) {
                return false; // `cxx:449-452`
            }
            // `cxx:385-408`: `BRep_Tool::CurveOnSurface` at index 1 and 2; `jj`
            // counts the representations found.
            let jj = reps.len().min(2);
            let Some((surf, pc, f, l)) = reps.first().cloned() else {
                return false;
            };
            // `cxx:410-418`: `Adaptor3d_CurveOnSurface` over the pcurve and the
            // surface; the port's `CurveOnSurface` is the same evaluator.
            let cos: Arc<dyn Curve> = Arc::new(CurveOnSurface::new(pc.clone(), surf.clone(), f, l));
            // `evaluateMaxSegment(MaxSegment, CurveOnSurface)` (`cxx:430`).
            let max_segment = evaluate_max_segment(0, &surf, &pc);
            let result = geom_lib::build_curve3d(
                tolerance,
                cos.as_ref(),
                &surf,
                &pc,
                f,
                l,
                Shape::C1,
                14,
                max_segment,
            );
            // `cxx:433-438`. The commented-out `max(tolerance, max_deviation)`
            // of this OCCT revision is a no-op: `Tolerance` is already
            // `max(1.e-5, BRep_Tool::Tolerance(edge))`.
            let max_deviation = reg.edge_tolerance(&edge.0).max(tolerance);
            let Some(c3d) = result.curve else {
                return false; // `cxx:437-440`
            };
            let mut g = EdgeGeom::new(c3d, f, l);
            g.tolerance = max_deviation;
            g.same_parameter = reg.same_parameter(&edge.0);
            g.same_range = reg.same_range(&edge.0);
            g.degenerated = reg.is_degenerated_edge(&edge.0);
            reg.set_edge(&edge.0, g);
            // `cxx:441-447`: with only one curve-on-surface representation the
            // edge "can be qualified sameparameter".
            if jj == 1 {
                reg.set_same_parameter(&edge.0, true);
            }
        }

        // `ShapeBuild_Edge.cxx:728-733`.
        if reg.same_range(&edge.0) {
            let (f, l) = reg.edge_parameters(&edge.0);
            set_rep_ranges(&edge.0, f, l);
        }
        // `cxx:734-738`.
        let Some(c3d) = reg.edge_curve(&edge.0) else {
            return false;
        };
        // `cxx:741-758` (OCC966): a non-periodic curve whose own range is
        // narrower than the edge range clamps the edge range and clears
        // SameRange.
        if !c3d.is_periodic() {
            let (mut f, mut l) = reg.edge_parameters(&edge.0);
            let mut is_less = false;
            if f < c3d.first_parameter() {
                is_less = true;
                f = c3d.first_parameter();
            }
            if l > c3d.last_parameter() {
                is_less = true;
                l = c3d.last_parameter();
            }
            if is_less {
                self.set_range3d(edge, f, l);
                reg.set_same_range(&edge.0, false);
            }
        }
        true
    }

    /// `ReassignPCurve(edge, old, sub)` (`cxx:530-592`).
    pub fn reassign_pcurve(&self, edge: &Edge, old: &Face, sub: &Face) -> bool {
        let mut npcurves = count_pcurves(edge, old);
        let (pc, f, l) = match boptools_2d::curve_on_surface_range(edge, old) {
            Some(v) => v,
            None => return false, // cxx:540-543
        };
        if npcurves == 0 {
            npcurves = 1; // cxx:544-547
        }

        let reg = GeometryRegistry::global();
        let old_key = GeometryRegistry::shape_key(&old.0);
        let sub_key = GeometryRegistry::shape_key(&sub.0);

        // cxx:552-564: if the pcurve was only one, remove; else leave second one.
        if npcurves > 1 {
            let erev = reversed(edge);
            if let Some((pc2, f2, l2)) = boptools_2d::curve_on_surface_range(&erev, old) {
                reg.set_edge_pcurve(&edge.0, old_key, pc2);
                reg.set_pcurve_range(&edge.0, old_key, f2, l2);
            }
        } else {
            self.remove_pcurve(edge, old);
        }

        // cxx:566-587: if edge does not have yet pcurves on sub, just add; else add as first.
        let npcs = count_pcurves(edge, sub);
        if npcs < 1 {
            reg.set_edge_pcurve(&edge.0, sub_key, pc.clone());
        } else {
            let erev = reversed(edge);
            if let Some((pcs, _cf, _cl)) = boptools_2d::curve_on_surface_range(&erev, sub) {
                if edge.0.orientation() == Orientation::Reversed {
                    reg.set_edge_pcurves(&edge.0, sub_key, vec![pcs, pc.clone()]);
                } else {
                    reg.set_edge_pcurves(&edge.0, sub_key, vec![pc.clone(), pcs]);
                }
            } else {
                reg.set_edge_pcurve(&edge.0, sub_key, pc.clone());
            }
        }
        reg.set_pcurve_range(&edge.0, sub_key, f, l); // cxx:589
        true
    }


    /// `TransformPCurve(pcurve, trans, uFact, aFirst, aLast)`
    /// (`ShapeBuild_Edge.cxx:596-698`): the `uFact == 1` return and the
    /// `uFact != 1` affinity branch (Line / Bezier / BSpline / Conic).
    pub fn transform_pcurve(
        &self,
        pcurve: &Arc<dyn Curve2d>,
        trans: &GpTrsf2d,
        u_fact: f64,
        a_first: &mut f64,
        a_last: &mut f64,
    ) -> Arc<dyn Curve2d> {
        let mut result = pcurve.clone_dyn();
        if trans.form() != TrsfForm::Identity {
            result.transform(trans); // cxx:605
            *a_first = result.transformed_parameter(*a_first, trans); // cxx:606
            *a_last = result.transformed_parameter(*a_last, trans); // cxx:607
        }
        if u_fact == 1.0 {
            return Arc::from(result); // cxx:609-612
        }
        // cxx:614-618: unwrap a trimmed curve to its basis.
        // cxx:620-623: `gp_GTrsf2d tMatu; tMatu.SetAffinity(gp::OY2d(), uFact);`
        let mut t_matu = GpGTrsf2d::identity();
        t_matu.set_affinity(
            &GpAx2d::new(GpPnt2d::new(0.0, 0.0), GpDir2d::new(0.0, 1.0).expect("dir")),
            u_fact,
        );
        // cxx:624-641: Line.
        let reference: &dyn Curve2d = result.trimmed_basis().unwrap_or(result.as_ref());
        if reference.gp_lin2d().is_some() {
            let pf = reference.d0(*a_first);
            let pl = reference.d0(*a_last);
            let pf2 = t_matu.transforms(&pf);
            let pl2 = t_matu.transforms(&pl);
            let v = GpVec2d::new(pl2.x() - pf2.x(), pl2.y() - pf2.y());
            if let Ok(dir) = GpDir2d::from_vec2d(&v) {
                let new_line = GpLin2d::from_pnt_dir(pf2, dir);
                *a_first = occt_core::elib::clib2d::line_parameter_ax2d(new_line.position(), &pf2);
                *a_last = occt_core::elib::clib2d::line_parameter_ax2d(new_line.position(), &pl2);
                return Arc::new(occt_geom2d::line::Geom2dLine::new(*new_line.position()));
            }
        }
        // cxx:643-655: a Bezier keeps its own type; `down_cast` returns the
        // same `result` (the unwrapped basis), whose poles the affinity is
        // applied to in place, and the Bezier itself is returned.
        if reference.is_bezier2d() {
            let mut bez = reference.clone_dyn();
            if let Some(poles) = bez.poles2d() {
                let scaled: Vec<GpPnt2d> = poles.iter().map(|p| t_matu.transforms(p)).collect();
                bez.set_poles2d(&scaled);
            }
            return Arc::from(bez);
        }
        // cxx:657-678: a conic is trimmed to `[aFirst, aLast]` and approximated
        // by `Geom2dConvert_ApproxCurve`; when that has no result,
        // `CurveToBSplineCurve(Convert_QuasiAngular)` is the fallback. The
        // resulting range is copied back into `aFirst`/`aLast`.
        let is_conic = reference.gp_circ2d().is_some()
            || reference.gp_elips2d().is_some()
            || reference.gp_hypr2d().is_some()
            || reference.gp_parab2d().is_some();
        if is_conic {
            let tcurve = Geom2dTrimmedCurve::new_sense(
                Arc::from(reference.clone_dyn()),
                *a_first,
                *a_last,
                true,
                false,
            );
            let approx = Geom2dConvertApproxCurve::new(&tcurve, APPROXIMATION, Shape::C1, 100, 6);
            let mut bs = match approx.curve() {
                Some(c) => c.clone(),
                None => match occt_geom2d::geom2d_convert::curve_to_bspline_curve_bspl(
                    &tcurve,
                    ParameterisationType::QuasiAngular,
                ) {
                    Some(b) => b,
                    None => return Arc::from(result),
                },
            };
            *a_first = bs.first_parameter();
            *a_last = bs.last_parameter();
            if let Some(poles) = bs.poles2d() {
                let scaled: Vec<GpPnt2d> = poles.iter().map(|p| t_matu.transforms(p)).collect();
                bs.set_poles2d(&scaled);
            }
            return Arc::new(bs);
        }
        // cxx:679-687: any other non-B-spline curve goes through
        // `CurveToBSplineCurve(Convert_QuasiAngular)`.
        if !reference.is_bspline2d() {
            if let Some(mut bs) = occt_geom2d::geom2d_convert::curve_to_bspline_curve_bspl(
                reference,
                ParameterisationType::QuasiAngular,
            ) {
                if let Some(poles) = bs.poles2d() {
                    let scaled: Vec<GpPnt2d> = poles.iter().map(|p| t_matu.transforms(p)).collect();
                    bs.set_poles2d(&scaled);
                }
                return Arc::new(bs);
            }
        }
        // cxx:687-697: a B-spline (`down_cast<Geom2d_BSplineCurve>(result)`)
        // transforms its own poles.
        let mut bs = reference.clone_dyn();
        if let Some(poles) = bs.poles2d() {
            let scaled: Vec<GpPnt2d> = poles.iter().map(|p| t_matu.transforms(p)).collect();
            bs.set_poles2d(&scaled);
        }
        Arc::from(bs)
    }
}