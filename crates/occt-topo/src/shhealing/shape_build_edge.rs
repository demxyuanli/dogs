//! `ShapeBuild_Edge` subset needed by `ShapeFix_ComposeShell::DispatchWires`:
//! `Copy` (`ShapeBuild_Edge.cxx:417-426`), `RemovePCurve` (`:430-443`),
//! `SetRange3d` (`:338-356`), `ReassignPCurve` (`:530-592`) and the
//!
//! `TransformPCurve` (`:596-698`): the `uFact == 1` return and the `Line` arm of
//! the `uFact != 1` affinity branch (`:624-641`) are ported; the
//! Bezier/BSpline/Conic pole rewriting (`:642-698`, needs
//! `Geom2dConvert_ApproxCurve` / `CurveToBSplineCurve` / pole setters) is UNPORTED.

use std::sync::Arc;

use occt_core::gp::{GpAx2d, GpDir2d, GpGTrsf2d, GpLin2d, GpPnt2d, GpTrsf2d, GpVec2d, TrsfForm};
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::boptools_2d;
use crate::shape::{Edge, Face};
use crate::tgeometry::GeometryRegistry;

/// `ShapeBuild_Edge` (`ShapeBuild_Edge.hxx`).
pub struct ShapeBuildEdge;

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
    /// (`ShapeBuild_Edge.cxx:596-698`), through the `uFact == 1` return.
    ///
    /// UNPORTED: the `uFact != 1` affine branch (`:614-698`).
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
        // cxx:642-698: a Bezier is converted to a BSpline through
        // `Geom2dConvert::CurveToBSplineCurve` (the exact poles/knots), then
        // the affinity is applied to the BSpline's poles.
        if reference.is_bezier2d() || reference.is_bspline2d() {
            let mut bs = if reference.is_bspline2d() {
                result.clone_dyn()
            } else {
                match occt_geom2d::geom2d_convert::curve_to_bspline_curve(reference) {
                    Some(b) => b,
                    None => return Arc::from(result),
                }
            };
            if let Some(poles) = bs.poles2d() {
                let scaled: Vec<GpPnt2d> = poles.iter().map(|p| t_matu.transforms(p)).collect();
                bs.set_poles2d(&scaled);
            }
            return Arc::from(bs);
        }
        // UNPORTED (cxx:657-686): the `Geom2d_Conic` path converts through
        // `Geom2dConvert_ApproxCurve` / `CurveToBSplineCurve(Convert_QuasiAngular)`.
        Arc::from(result)
    }
}