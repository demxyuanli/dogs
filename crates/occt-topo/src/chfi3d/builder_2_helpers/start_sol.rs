//! The `StartSol` file statics of `ChFi3d_Builder_2.cxx`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx` (OCCT
//! 8.0.0): `BonVoisin` (`:213-318`) and `Projection` (`:320-370`).
//!
//! | OCCT | here |
//! |---|---|
//! | `BonVoisin(Point, HS, F, plane, cured, XDep, YDep, EFMap, tol3d)` (`:213-318`) | [`bon_voisin`] |
//! | `Projection(PExt, P, C, W, Tol)` (`:320-370`) | [`projection`] |
//!
//! Divergences that are not heuristics:
//! - `occ::handle<BRepAdaptor_Surface>& HS` is `&mut SurfaceAdaptor`: OCCT
//!   re-points the adaptor at a neighbour face (`HS->Initialize(ff)`), which the
//!   port spells as a reassignment,
//! - `occ::handle<GeomAdaptor_Surface> plane` is `&dyn Surface` (the
//!   `Geom_Plane` the caller builds from the guideline point and tangent),
//! - `Extrema_ExtPC& PExt` is the curve / range pair the port's
//!   [`extrema_ext_pc_range_done`] takes: OCCT initialises the object once
//!   (`StartSol`, `cxx:935-940`) and calls `Perform(P)` inside `Projection`, so
//!   the port passes the same `(els, u_inf, u_sup)` every time,
//! - `Extrema_LocateExtPC ext(P, C, W, Tol/10)` is
//!   [`LocateExtPC::initialize`] + [`LocateExtPC::perform`] over the curve's own
//!   range (`Extrema_GLocateExtPC.hxx:110-122`),
//! - a null p-curve (`getCurveOnSurface` returning a null handle) or a failed
//!   adaptor evaluation answers `false` / leaves the output untouched, where
//!   OCCT would dereference the null handle.

use occt_core::gp::{GpPnt, GpVec};
use occt_core::kernel::geomabs::Shape;
use occt_geom::extrema_pc::{extrema_ext_pc_range_done, LocateExtPC};
use occt_geom::{Curve, Surface};

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_tool::BRepTool;
use crate::brep_tools;
use crate::chfi3d::builder_0_helpers::inter_plane_edge;
use crate::chfi3d::tangent_faces::is_tangent_faces;
use crate::chfi_ds::ChFiDSMap;
use crate::shape::{Edge, Face};
use crate::topexp::Explorer;
use crate::topo_tools::is_same;
use crate::topopebrep::tool::adaptor::{CurveAdaptor, SurfaceAdaptor};

/// `Projection(PExt, P, C, W, Tol)` (`cxx:320-370`).
///
/// Answers whether `w` holds a usable projection of `p` on `els` (the guideline
/// `ChFiDS_ElSpine`): `true` when the current parameter is already within `tol`,
/// otherwise after the local search (`Extrema_LocateExtPC`) or the global one
/// (`PExt.Perform`).
pub fn projection(
    els: &dyn Curve,
    u_inf: f64,
    u_sup: f64,
    p: &GpPnt,
    w: &mut f64,
    tol: f64,
) -> bool {
    // `Dist2 = C.Value(W).SquareDistance(P)` (`cxx:327`).
    let mut dist2 = els.d0(*w).square_distance(p);

    // It is checked if it is not already a solution (`cxx:329-332`).
    if dist2 < tol * tol {
        return true;
    }

    let mut ok = false;

    // `Extrema_LocateExtPC ext(P, C, W, Tol / 10)` (`cxx:338-353`). The
    // four-argument constructor calls `Initialize(C, TolF)`, which takes the
    // curve's own range.
    let mut ext = LocateExtPC::default();
    ext.initialize(els, els.first_parameter(), els.last_parameter(), tol / 10.0);
    ext.perform(p, *w);
    if ext.is_done() {
        let u = ext.point().parameter();
        let daux2 = els.d0(u).square_distance(p);
        if daux2 < dist2 {
            *w = u;
            dist2 = daux2;
            ok = true;
            if dist2 < tol * tol {
                return true;
            }
        }
    }

    // Global resolution (`cxx:356-368`): `PExt.Perform(P)`, then the scan over
    // `NbExt` keeping the strictly smaller square distance. `ii` starts at 1 in
    // OCCT; the port iterates the solutions in the same order.
    let (sols, done) = extrema_ext_pc_range_done(els, p, u_inf, u_sup);
    if done {
        for sol in &sols {
            if sol.sq_dist < dist2 {
                dist2 = sol.sq_dist;
                *w = sol.u;
                ok = true;
            }
        }
    }
    ok
}

/// `BonVoisin(Point, HS, F, plane, cured, XDep, YDep, EFMap, tol3d)`
/// (`cxx:213-318`): tries to move the contact point `(XDep, YDep)` of the
/// support `HS` onto a face tangent to `F` across an edge that the plane built
/// on `point` cuts.
///
/// `F` and `HS` are paired (`HS` is the adaptor of `F`, `HS->Face()`); only
/// `HS` is re-pointed, and `(XDep, YDep)` become the p-curve parameters on the
/// new face. `false` is OCCT's `bonvoisin = false` (no neighbour found).
pub fn bon_voisin(
    point: &GpPnt,
    hs: &mut SurfaceAdaptor,
    f: &Face,
    plane: &dyn Surface,
    cured: &Edge,
    x_dep: &mut f64,
    y_dep: &mut f64,
    ef_map: &ChFiDSMap,
    tol3d: f64,
) -> bool {
    // `gp_Pnt papp = HS->Value(XDep, YDep); double dist = RealLast();`
    // (`cxx:225-226`).
    let Some(surface) = hs.surface().cloned() else {
        return false;
    };
    let papp = surface.d0(*x_dep, *y_dep);
    let mut dist = f64::MAX;
    let mut found = false;

    // `for (Ex.Init(F, TopAbs_EDGE); Ex.More(); Ex.Next())` (`cxx:231-312`).
    let mut ex = Explorer::new(f.shape(), ShapeType::Edge);
    'edges: while ex.more() {
        let ecur = Edge(ex.current().clone());
        if !is_same(&ecur.0, &cured.0) {
            let hc = CurveAdaptor::of_edge(&ecur);
            let tolc = hc.resolution(tol3d);
            let mut winter = 0.0;
            if inter_plane_edge(plane, &hc, &mut winter, true, tolc) {
                let Some(np) = hc.d0(winter) else {
                    ex.next();
                    continue;
                };
                let ndist = np.square_distance(&papp);
                if ndist < dist {
                    // `bool isclosed = BRep_Tool::IsClosed(ecur, F);`
                    // `bool isreallyclosed = BRepTools::IsReallyClosed(ecur, F);`
                    let is_closed = BRepTool::is_closed_edge_face(&ecur, f);
                    let is_really_closed = brep_tools::is_really_closed(&ecur, f);
                    for ff_shape in ef_map.find_from_key(&ecur.0) {
                        let Some(ff) = Face::wrap(ff_shape.clone()) else {
                            continue;
                        };
                        let issame = is_same(&ff.0, &f.0);
                        // `bool istg = ChFi3d::IsTangentFaces(ecur, ff, F);`
                        // (`cxx:256-258`, default `GeomAbs_G1`).
                        let istg =
                            is_tangent_faces(&ecur, &ff, f, Shape::G1).unwrap_or(false);
                        if (!issame || (issame && is_really_closed)) && istg {
                            found = true;
                            let mut newe = Edge(ecur.0.oriented(Orientation::Forward));
                            dist = ndist;
                            *hs = SurfaceAdaptor::of_face(&ff);
                            let pc = if is_closed && !is_really_closed {
                                // `fff.Orientation(TopAbs_FORWARD)` then the edge of
                                // `fff` that `newe.IsSame` identifies (`cxx:269-279`).
                                let fff = Face(ff.0.oriented(Orientation::Forward));
                                let mut matched = None;
                                let mut ex2 = Explorer::new(fff.shape(), ShapeType::Edge);
                                while ex2.more() {
                                    let candidate = ex2.current().clone();
                                    if is_same(&newe.0, &candidate) {
                                        newe = Edge(candidate);
                                        matched = curve_on_surface_range(&newe, &fff);
                                        break;
                                    }
                                    ex2.next();
                                }
                                matched
                            } else {
                                curve_on_surface_range(&newe, &ff)
                            };
                            let Some((pc, _, _)) = pc else {
                                break 'edges;
                            };
                            let uv = pc.d0(winter);
                            *x_dep = uv.x();
                            *y_dep = uv.y();
                            if issame {
                                // `HS->D1(XDep, YDep, spt, sdu, sdv)` (`cxx:287-289`).
                                let Some((_, sdu, sdv)) = hs.d1(*x_dep, *y_dep) else {
                                    break 'edges;
                                };
                                let nors = sdu.crossed(&sdv);
                                let Some((cpt, cd)) = hc.d1(winter) else {
                                    break 'edges;
                                };
                                let vref = GpVec::from_pnts(point, &cpt);
                                let fff = Face(ff.0.oriented(Orientation::Forward));
                                if vref.dot(&nors.crossed(&cd)) < 0.0 {
                                    // `newe.Orientation(TopAbs_REVERSED)` (`cxx:299`).
                                    newe = Edge(newe.0.oriented(Orientation::Reversed));
                                }
                                let Some((pc, _, _)) = curve_on_surface_range(&newe, &fff)
                                else {
                                    break 'edges;
                                };
                                let uv = pc.d0(winter);
                                *x_dep = uv.x();
                                *y_dep = uv.y();
                            }
                            break;
                        }
                    }
                }
            }
        }
        ex.next();
    }
    let mut bonvoisin = true;
    if !found {
        // `if (!found) { bonvoisin = false; }` (`cxx:313-316`).
        bonvoisin = false;
    }
    bonvoisin
}
