//! `ChFi3d_Builder::StartSol`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:900-1108`
//! (OCCT 8.0.0).
//!
//! `StartSol` parses about ten parameters along the guideline, builds the two
//! contact points of the first section on the faces adjacent to the spine edge,
//! and, when that fails, retries on the neighbour faces found across a plane
//! normal to the guideline (`ChFi3d_Builder_2.cxx:1018-1105`).
//!
//! Divergences that are not heuristics:
//! - `occ::handle<ChFiDS_Spine>& Spine` is a [`SpineHandle`]: every
//!   `Spine->Member` call takes the write lock for the duration of that call,
//!   exactly where OCCT's reference keeps the object alive between calls,
//! - `occ::handle<BRepAdaptor_Surface>& HS1 / HS2` are `&mut SurfaceAdaptor`
//!   because `ConexFaces` re-points both handles (`cxx:864-888`),
//! - `PerformFirstSection` is pure virtual (`hxx:577-588`); it is dispatched
//!   through [`ChFi3dSurfOps`],
//! - the trailing `throw Standard_Failure("StartSol echec")` is the `Err` of
//!   the port after `Spine->SetErrorStatus(ChFiDS_StartsolFailure)`,
//! - a null p-curve (`getCurveOnSurface` returning a null handle) or a failing
//!   adaptor evaluation stops the attempt and answers the same failure as OCCT's
//!   `throw`, instead of dereferencing the null handle,
//! - `AS.UResolution(TolE)` is
//!   [`occt_geom::approx_same_parameter::u_resolution`] on the face's surface,
//!   the same `GeomAdaptor_Surface` the OCCT `BRepAdaptor_Surface` forwards to.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpPln, GpPnt2d, GpVec2d};
use occt_geom::approx_same_parameter::{u_resolution, v_resolution};
use occt_geom::{Curve, GeomPlane};
use occt_math::MathVector;

use crate::abs::{Orientation, State};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi3d::builder_2_helpers::{bon_voisin, projection};
use crate::chfi3d::side::next_side;
use crate::chfi_ds::{ElSpine, ErrorStatus, StripeHandle};
use crate::shape::{Edge, Face};
use crate::topo_tools::is_same;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::call_perform_surf::initialize_topol_tool;
use super::surf_ops::ChFi3dSurfOps;
use super::ChFi3dBuilder;

/// `gp_Vec2d::Rotate(M_PI / 2)` (`gp_Vec2d.cxx`): OCCT rotates in place by
/// `(x cos - y sin, x sin + y cos)`; with the exact quadrant angle the cosine is
/// `0` and the sine is `1`.
fn rotate_half_pi(v: GpVec2d) -> GpVec2d {
    GpVec2d::new(-v.y(), v.x())
}

impl ChFi3dBuilder {
    /// `StartSol(Stripe, HGuide, HS1, HS2, I1, I2, P1, P2, First)`
    /// (`cxx:900-1108`).
    ///
    /// On success `first` is the guideline parameter of the retained section and
    /// `(p1, p2)` its two support parameters; on failure the spine's error
    /// status is `ChFiDS_StartsolFailure` and `Err` carries the OCCT message.
    ///
    /// `Or1` / `Or2` (`cxx:932`) and the `iedge = 0` reset (`cxx:1020`) are
    /// written before their first read in OCCT as well; the port keeps the
    /// statements and allows the dead-assignment lint they raise.
    #[allow(clippy::too_many_arguments, unused_assignments)]
    pub fn start_sol_el_spine(
        &self,
        ops: &mut dyn ChFi3dSurfOps,
        stripe: &StripeHandle,
        hguide: &Arc<ElSpine>,
        hs1: &mut SurfaceAdaptor,
        hs2: &mut SurfaceAdaptor,
        i1: &Rc<RefCell<TopolTool>>,
        i2: &Rc<RefCell<TopolTool>>,
        p1: &mut GpPnt2d,
        p2: &mut GpPnt2d,
        first: &mut f64,
    ) -> Result<(), String> {
        // `occ::handle<ChFiDS_Spine>& Spine = Stripe->ChangeSpine();`
        // `ChFiDS_ElSpine& els = *HGuide;` (`cxx:910-911`).
        let spine_handle = {
            let guard = stripe
                .read()
                .map_err(|_| "StartSol: stripe lock poisoned".to_string())?;
            guard
                .spine()
                .cloned()
                .ok_or_else(|| "StartSol: stripe without spine".to_string())?
        };
        let els: &dyn Curve = hguide.as_ref();

        // `int nbed = Spine->NbEdges(); int nbessaimax = 3 * nbed;
        // if (nbessaimax < 10) nbessaimax = 10; double unsurnbessaimax = 1. /
        // nbessaimax;` (`cxx:912-918`).
        let (nbed, first_par_1, last_par_1, first_par_nbed, last_par_nbed) = {
            let guard = spine_handle
                .read()
                .map_err(|_| "StartSol: spine lock poisoned".to_string())?;
            (
                guard.nb_edges() as i32,
                guard.first_parameter_of(1),
                guard.last_parameter_of(1),
                guard.last_parameter_of(guard.nb_edges() as i32),
                guard.first_parameter_of(guard.nb_edges() as i32),
            )
        };
        let (first_par_1, last_par_1, first_par_nbed, last_par_nbed) = (
            first_par_1.ok_or_else(|| "StartSol: empty spine".to_string())?,
            last_par_1.ok_or_else(|| "StartSol: empty spine".to_string())?,
            first_par_nbed.ok_or_else(|| "StartSol: empty spine".to_string())?,
            last_par_nbed.ok_or_else(|| "StartSol: empty spine".to_string())?,
        );
        if nbed < 1 {
            return Err("StartSol: empty spine".into());
        }
        let nbessaimax = (3 * nbed).max(10);
        let unsurnbessaimax = 1.0 / f64::from(nbessaimax);
        let wf = 0.9981 * first_par_1 + 0.0019 * last_par_1;
        let wl = 0.9973 * last_par_nbed + 0.0027 * first_par_nbed;

        // `int iedge = 0; int RC = Stripe->Choix();` (`cxx:925-927`).
        let (rc, stripe_or1, stripe_or2) = {
            let guard = stripe
                .read()
                .map_err(|_| "StartSol: stripe lock poisoned".to_string())?;
            (
                guard.choix(),
                guard.orientation_on_face1(),
                guard.orientation_on_face2(),
            )
        };

        let mut iedge: i32 = 0;
        let mut or1 = Orientation::Forward;
        let mut or2 = Orientation::Forward;
        let mut choix: i32 = 0;
        // `math_Vector SolDep(1, 4)`.
        let mut sol_dep = MathVector::new(1, 4);
        let mut pos1 = State::Unknown;
        let mut pos2 = State::Unknown;

        // `PExt.Initialize(els, Spine->FirstParameter(1),
        // Spine->LastParameter(nbed), Precision::Confusion())` (`cxx:936-940`):
        // the port passes the same curve and range to `Projection`.
        let (u_inf, u_sup) = (first_par_1, last_par_nbed);

        // `for (nbessai = 0; nbessai <= nbessaimax; nbessai++)` (`cxx:942-1017`).
        for nbessai in 0..=nbessaimax {
            let t = f64::from(nbessai) * unsurnbessaimax;
            let w = wf * (1.0 - t) + wl * t;
            let ie = {
                let guard = spine_handle
                    .read()
                    .map_err(|_| "StartSol: spine lock poisoned".to_string())?;
                guard.index(w, true)
            };
            let Some(ie) = ie else {
                continue;
            };

            let mut cured: Option<Edge> = None;
            let mut faces: Option<(Face, Face)> = None;
            let mut tol_e = 1.0e-7;
            if iedge != ie {
                iedge = ie;
                // `cured = Spine->Edges(iedge);` (`cxx:950`).
                let Some(edge) = ({
                    let guard = spine_handle
                        .read()
                        .map_err(|_| "StartSol: spine lock poisoned".to_string())?;
                    guard.edges(iedge as usize).cloned()
                }) else {
                    continue;
                };
                let edge = Edge(edge);
                cured = Some(edge);
                // `TolE = BRep_Tool::Tolerance(cured);` (`cxx:952`).
                tol_e = crate::brep_tool::BRepTool::edge_tolerance(cured.as_ref().unwrap());
                // `ConexFaces(Spine, iedge, HS1, HS2)` (`cxx:953-958`).
                let Some((a1, a2)) = self.cone_x_faces(&spine_handle, iedge as usize) else {
                    continue;
                };
                *hs1 = a1;
                *hs2 = a2;
                let (a, b) = match (hs1.face().cloned(), hs2.face().cloned()) {
                    (Some(a), Some(b)) => (a, b),
                    _ => continue,
                };
                or1 = a.0.orientation();
                or2 = b.0.orientation();
                choix = next_side(&mut or1, &mut or2, stripe_or1, stripe_or2, rc);
                faces = Some((a, b));
            }
            let (Some(cured), Some((f1, f2))) = (cured.as_ref(), faces.as_ref()) else {
                continue;
            };

            // `Spine->Parameter(iedge, w, woned, true)` (`cxx:961`).
            let Some(woned) = ({
                let mut guard = spine_handle
                    .write()
                    .map_err(|_| "StartSol: spine lock poisoned".to_string())?;
                guard.parameter_of(iedge, w, true)
            }) else {
                continue;
            };

            // `cured.Orientation(TopAbs_FORWARD);
            // f1forward.Orientation(TopAbs_FORWARD);` (`cxx:962-965`).
            let cured_forward = Edge(cured.0.oriented(Orientation::Forward));
            let f1_forward = Face(f1.0.oriented(Orientation::Forward));
            let f2_forward = Face(f2.0.oriented(Orientation::Forward));
            // `PC = getCurveOnSurface(cured, f1forward)` (`cxx:966`).
            let Some((pc, _, _)) = curve_on_surface_range(&cured_forward, &f1_forward) else {
                continue;
            };
            // `I1->Initialize((const handle(Adaptor3d_Surface)&)HS1)` (`cxx:967`).
            if initialize_topol_tool(i1, hs1).is_err() {
                continue;
            }
            // `PC->D1(woned, P1, derive)` (`cxx:968`).
            let (mut point1, derive) = pc.d1(woned);
            *p1 = point1;

            // `if (derive.Magnitude() > Precision::PConfusion())` (`cxx:970-992`).
            if derive.magnitude() > occt_core::precision::PCONFUSION {
                let mut derive = rotate_half_pi(derive);
                // `AS.Initialize(f1); ResU = AS.UResolution(TolE); ResV =
                // AS.VResolution(TolE);` (`cxx:975-978`).
                let Some(face_surface) = crate::brep_tool::BRepTool::face_surface(&f1) else {
                    continue;
                };
                let res_u = u_resolution(face_surface.as_ref(), tol_e);
                let res_v = v_resolution(face_surface.as_ref(), tol_e);
                derive.multiply_scalar(2.0 * (derive.x().abs() * res_u + derive.y().abs() * res_v));
                let res_min = res_u.min(res_v);
                // `P2d = P1.Translated(derive); if (I1->Classify(P2d,
                // std::min(ResU, ResV), false) == TopAbs_IN) P1 = P2d;
                // else { P2d = P1.Translated(-derive); ... }` (`cxx:979-991`).
                let p2d = point1.translated_vec(&derive);
                if i1.borrow_mut().classify(p2d, res_min, false) == crate::fclass2d::FaceState::In {
                    point1 = p2d;
                } else {
                    let p2d = point1.translated_vec(&derive.reversed());
                    if i1.borrow_mut().classify(p2d, res_min, false)
                        == crate::fclass2d::FaceState::In
                    {
                        point1 = p2d;
                    }
                }
                *p1 = point1;
            }

            // `if (f1.IsSame(f2)) cured.Orientation(TopAbs_REVERSED);`
            // (`cxx:993-996`).
            let cured_for_f2 = if is_same(&f1.0, &f2.0) {
                Edge(cured.0.oriented(Orientation::Reversed))
            } else {
                cured.clone()
            };
            // `PC = getCurveOnSurface(cured, f2forward); P2 = PC->Value(woned);`
            // (`cxx:997-998`).
            let Some((pc2, _, _)) = curve_on_surface_range(&cured_for_f2, &f2_forward) else {
                continue;
            };
            *p2 = pc2.d0(woned);
            // `I2->Initialize(HSon2)` (`cxx:1000`).
            if initialize_topol_tool(i2, hs2).is_err() {
                continue;
            }

            // `SolDep(1..4) = P1.X(), P1.Y(), P2.X(), P2.Y()` (`cxx:1002-1005`).
            sol_dep.set_value(1, p1.x());
            sol_dep.set_value(2, p1.y());
            sol_dep.set_value(3, p2.x());
            sol_dep.set_value(4, p2.y());

            // `const BRepAdaptor_Curve& Ced = Spine->CurrentElementarySpine(iedge);
            // gp_Pnt pnt = Ced.Value(woned);` (`cxx:1006-1007`).
            let Some(pnt) = ({
                let mut guard = spine_handle
                    .write()
                    .map_err(|_| "StartSol: spine lock poisoned".to_string())?;
                guard.current_elementary_spine(iedge).and_then(|c| c.d0(woned))
            }) else {
                continue;
            };

            // `if (Projection(PExt, pnt, els, w, tolapp3d) &&
            // PerformFirstSection(Spine, HGuide, Choix, HS1, HS2, I1, I2, w,
            // SolDep, Pos1, Pos2))` (`cxx:1009-1016`).
            let mut w_proj = w;
            if projection(els, u_inf, u_sup, &pnt, &mut w_proj, self.tolapp3d)
                && ops.perform_first_section(
                    &spine_handle,
                    hguide,
                    choix,
                    hs1,
                    hs2,
                    i1,
                    i2,
                    w_proj,
                    &mut sol_dep,
                    &mut pos1,
                    &mut pos2,
                )
            {
                *p1 = GpPnt2d::new(sol_dep.value(1), sol_dep.value(2));
                *p2 = GpPnt2d::new(sol_dep.value(3), sol_dep.value(4));
                *first = w_proj;
                return Ok(());
            }
        }

        // `iedge = 0;` then the neighbour-face retry (`cxx:1020-1105`).
        iedge = 0;
        for nbessai in 0..=nbessaimax {
            let t = f64::from(nbessai) * unsurnbessaimax;
            let w = wf * (1.0 - t) + wl * t;
            let Some(iedge_new) = ({
                let guard = spine_handle
                    .read()
                    .map_err(|_| "StartSol: spine lock poisoned".to_string())?;
                guard.index(w, true)
            }) else {
                continue;
            };
            iedge = iedge_new;
            let Some(cured) = ({
                let guard = spine_handle
                    .read()
                    .map_err(|_| "StartSol: spine lock poisoned".to_string())?;
                guard.edges(iedge as usize).cloned()
            }) else {
                continue;
            };
            let cured = Edge(cured);
            // `ConexFaces(Spine, iedge, HS1, HS2)` (`cxx:1027`).
            let Some((a1, a2)) = self.cone_x_faces(&spine_handle, iedge as usize) else {
                continue;
            };
            *hs1 = a1;
            *hs2 = a2;
            let (mut f1, mut f2) = match (hs1.face().cloned(), hs2.face().cloned()) {
                (Some(a), Some(b)) => (a, b),
                _ => continue,
            };
            or1 = f1.0.orientation();
            or2 = f2.0.orientation();
            choix = next_side(&mut or1, &mut or2, stripe_or1, stripe_or2, rc);

            // `Spine->Parameter(iedge, w, woned, true)` (`cxx:1033`).
            let Some(woned) = ({
                let mut guard = spine_handle
                    .write()
                    .map_err(|_| "StartSol: spine lock poisoned".to_string())?;
                guard.parameter_of(iedge, w, true)
            }) else {
                continue;
            };

            let f1_forward = Face(f1.0.oriented(Orientation::Forward));
            let f2_forward = Face(f2.0.oriented(Orientation::Forward));
            // `PC = getCurveOnSurface(cured, f1forward); P1 = PC->Value(woned);`
            // (`cxx:1039-1041`).
            let Some((pc, _, _)) = curve_on_surface_range(&cured, &f1_forward) else {
                continue;
            };
            *p1 = pc.d0(woned);
            let Some((pc2, _, _)) = curve_on_surface_range(&cured, &f2_forward) else {
                continue;
            };
            *p2 = pc2.d0(woned);
            if initialize_topol_tool(i1, hs1).is_err() || initialize_topol_tool(i2, hs2).is_err() {
                continue;
            }

            sol_dep.set_value(1, p1.x());
            sol_dep.set_value(2, p1.y());
            sol_dep.set_value(3, p2.x());
            sol_dep.set_value(4, p2.y());

            let Some(pnt) = ({
                let mut guard = spine_handle
                    .write()
                    .map_err(|_| "StartSol: spine lock poisoned".to_string())?;
                guard.current_elementary_spine(iedge).and_then(|c| c.d0(woned))
            }) else {
                continue;
            };

            // `if (Projection(PExt, pnt, els, w, tolapp3d))` (`cxx:1057-1104`).
            let mut w_proj = w;
            if !projection(els, u_inf, u_sup, &pnt, &mut w_proj, self.tolapp3d) {
                continue;
            }
            ops.perform_first_section(
                &spine_handle,
                hguide,
                choix,
                hs1,
                hs2,
                i1,
                i2,
                w_proj,
                &mut sol_dep,
                &mut pos1,
                &mut pos2,
            );

            // `gp_Pnt P; gp_Vec V; HGuide->D1(w, P, V);` (`cxx:1059-1061`).
            let (p, v) = hguide.d1(w_proj);
            let Ok(dir) = GpDir::from_vec(&v) else {
                continue;
            };
            let Ok(ax3) = GpAx3::from_point_dir(p, dir) else {
                continue;
            };
            // `handle(Geom_Plane) pl = new Geom_Plane(P, V);` +
            // `handle(GeomAdaptor_Surface) plane = new
            // GeomAdaptor_Surface(pl);` (`cxx:1062-1063`).
            let plane = GeomPlane::new(GpPln::new(ax3));

            // `bool bonvoisin = true, found = false; for (NbChangement = 1;
            // bonvoisin && (!found) && (NbChangement < 5); NbChangement++)`
            // (`cxx:1065-1099`).
            let mut bonvoisin = true;
            let mut found = false;
            let mut nb_changement = 1;
            while bonvoisin && !found && nb_changement < 5 {
                if pos1 != State::In {
                    // `BonVoisin(P, HS1, f1, plane, cured, SolDep(1), SolDep(2),
                    // myEFMap, tolapp3d)`: the two coordinates live in `SolDep`
                    // in OCCT, so the port reads them out and writes them back.
                    let mut sd1 = sol_dep.value(1);
                    let mut sd2 = sol_dep.value(2);
                    bonvoisin = bon_voisin(
                        &p,
                        hs1,
                        &f1,
                        &plane,
                        &cured,
                        &mut sd1,
                        &mut sd2,
                        &self.my_ef_map,
                        self.tolapp3d,
                    );
                    sol_dep.set_value(1, sd1);
                    sol_dep.set_value(2, sd2);
                }
                if pos2 != State::In && bonvoisin {
                    let mut sd3 = sol_dep.value(3);
                    let mut sd4 = sol_dep.value(4);
                    bonvoisin = bon_voisin(
                        &p,
                        hs2,
                        &f2,
                        &plane,
                        &cured,
                        &mut sd3,
                        &mut sd4,
                        &self.my_ef_map,
                        self.tolapp3d,
                    );
                    sol_dep.set_value(3, sd3);
                    sol_dep.set_value(4, sd4);
                }
                if bonvoisin {
                    // `f1 = HS1->Face(); f2 = HS2->Face();`
                    // (`cxx:1079-1080`).
                    let (nf1, nf2) = match (hs1.face().cloned(), hs2.face().cloned()) {
                        (Some(a), Some(b)) => (a, b),
                        _ => break,
                    };
                    f1 = nf1;
                    f2 = nf2;
                    or1 = f1.0.orientation();
                    or2 = f2.0.orientation();
                    choix = next_side(&mut or1, &mut or2, stripe_or1, stripe_or2, rc);
                    if initialize_topol_tool(i1, hs1).is_err()
                        || initialize_topol_tool(i2, hs2).is_err()
                    {
                        break;
                    }
                    if ops.perform_first_section(
                        &spine_handle,
                        hguide,
                        choix,
                        hs1,
                        hs2,
                        i1,
                        i2,
                        w_proj,
                        &mut sol_dep,
                        &mut pos1,
                        &mut pos2,
                    ) {
                        *p1 = GpPnt2d::new(sol_dep.value(1), sol_dep.value(2));
                        *p2 = GpPnt2d::new(sol_dep.value(3), sol_dep.value(4));
                        *first = w_proj;
                        found = true;
                    }
                }
                nb_changement += 1;
            }
            if found {
                return Ok(());
            }
        }

        // `Spine->SetErrorStatus(ChFiDS_StartsolFailure); throw
        // Standard_Failure("StartSol echec");` (`cxx:1106-1108`).
        if let Ok(mut guard) = spine_handle.write() {
            guard.set_error_status(ErrorStatus::StartsolFailure);
        }
        Err("StartSol echec".into())
    }
}

