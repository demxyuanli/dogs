//! `Draft_Modification::Perform`, second block: "Calculate new edges"
//! (`Draft_Modification_1.cxx:826-1445`).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | edge loop head, `BRep_Tool::Curve`, tangent probe | `_1.cxx:826-842` | [`DraftModification::perform_edges`] |
//! | `detrompeur` projections | `_1.cxx:852-889` | `detrompeur` |
//! | `KPart` fast arm (circle translated along the extrusion) | `_1.cxx:906-960` | `recompute_edge_curve` |
//! | surface / surface intersection arm | `_1.cxx:961-1331` | [`super::edges_intersection`] |
//! | tangency arm | `_1.cxx:1333-1425` | `recompute_tangent_edge` |
//! | store the curve / unwrap `Geom_TrimmedCurve` | `_1.cxx:1426-1444` | [`DraftModification::perform_edges`] |
//!
//! # Translation notes
//!
//! * `BRep_Tool::Curve(theEdge, L, f, l)` + `Transformed(L.Transformation())`
//!   (`cxx:838-839`) is [`BRepTool::edge_curve_world`].
//! * `TopExp::FirstVertex` / `LastVertex` (`cxx:853-855`) is
//!   [`crate::topo_tools::edge_vertices`], and `BRep_Tool::Parameter(V, E)`
//!   (`cxx:856-857`) is [`BRepTool::parameter_of_vertex`].
//! * `GeomAPI_ProjectPointOnSurf(P, S, Tol)` (`cxx:868-871`) is
//!   `geom_api::project_point_on_surface`.
//! * `Geom_Plane::Position()` (`cxx:911`, `:919`) is [`Surface::gp_pln`] +
//!   `GpPln::position()`; `Geom_SurfaceOfLinearExtrusion::BasisCurve` /
//!   `Direction` (`cxx:912-913`, `:920-921`) are
//!   [`Surface::extrusion_basis_curve`] / [`Surface::extrusion_direction`].
//! * The reference leaves `gp_Dir dirextr` default-constructed when neither
//!   face is a cylinder or an extrusion (`cxx:1350-1406`); `gp_Dir`'s default is
//!   `(0, 0, 1)`, which is what the port initialises.

use std::sync::Arc;

use occt_core::elib::clib;
use occt_core::gp::{GpDir, GpDir2d, GpPnt, GpPnt2d, GpVec};
use occt_core::precision::{Precision, CONFUSION};
use occt_geom::curve::Curve;
use occt_geom::geom_api::project_point_on_surface;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use crate::brep_tool::BRepTool;
use crate::draft::error_status::DraftErrorStatus;
use crate::shape::{Edge, Face};
use crate::topo_tools::edge_vertices;

use super::DraftModification;

impl DraftModification {
    /// `Draft_Modification_1.cxx:826-1445`.
    pub(super) fn perform_edges(&mut self) -> Result<(), ()> {
        for ii in 1..=self.my_e_map.extent() {
            let key = self.my_e_map.find_key(ii).clone();
            let Some(the_edge) = Edge::wrap(key) else {
                continue;
            };
            // `cxx:838-839`
            let Some(curve) = BRepTool::edge_curve_world(&the_edge) else {
                continue;
            };

            // `cxx:841-842` and the `IsTangent(ptfixe)` probe of `cxx:844-846`.
            let (new_geometry, geometry_is_null, is_tangent, tangent_point, first_face, second_face) = {
                let einf = self.my_e_map.find_from_index(ii);
                let mut p = GpPnt::new(0.0, 0.0, 0.0);
                let tgt = einf.is_tangent(&mut p);
                (
                    einf.new_geometry(),
                    einf.geometry().is_none(),
                    tgt,
                    p,
                    einf.first_face().cloned(),
                    einf.second_face().cloned(),
                )
            };

            if new_geometry && geometry_is_null {
                let new_c = if is_tangent {
                    self.recompute_tangent_edge(
                        &the_edge,
                        &curve,
                        &tangent_point,
                        first_face,
                        second_face,
                        ii,
                    )?
                } else {
                    self.recompute_edge_curve(&the_edge, &curve, first_face, second_face, ii)?
                };
                // `cxx:1426-1433`: `Geom_TrimmedCurve` unwrapped to its basis.
                let stored = new_c.trim_basis().unwrap_or_else(|| new_c.clone());
                *self.my_e_map.change_from_index(ii).change_geometry() = Some(stored);
            } else if !new_geometry {
                // `cxx:1435-1443`: keep the existing 3D curve.
                let stored = curve.trim_basis().unwrap_or_else(|| curve.clone());
                *self.my_e_map.change_from_index(ii).change_geometry() = Some(stored);
            }
        }
        Ok(())
    }

    /// `cxx:837-1331`: the two faces' new surfaces cut each other, except for
    /// the circle-on-plane `KPart` case which is solved in closed form.
    fn recompute_edge_curve(
        &mut self,
        the_edge: &Edge,
        c: &Arc<dyn Curve>,
        first_face: Option<Face>,
        second_face: Option<Face>,
        ii: usize,
    ) -> Result<Arc<dyn Curve>, ()> {
        // `cxx:847-850`
        let (Some(first_face), Some(second_face)) = (first_face, second_face) else {
            return self.fail_edge(the_edge);
        };
        let s1 = self
            .my_f_map
            .find_from_key(&first_face.0)
            .and_then(|fi| fi.geometry())
            .cloned();
        let s2 = self
            .my_f_map
            .find_from_key(&second_face.0)
            .and_then(|fi| fi.geometry())
            .cloned();
        let (Some(s1), Some(s2)) = (s1, s2) else {
            return self.fail_edge(the_edge);
        };

        // `cxx:852-889`: which of the two faces the edge really sits on.
        let detrompeur = self.detrompeur(the_edge, c, &s1, &s2, &first_face, &second_face);

        // `cxx:897-904`
        let s1 = s1.rectangular_trimmed_basis().unwrap_or(s1);
        let s2 = s2.rectangular_trimmed_basis().unwrap_or(s2);

        // `cxx:906-937`: the `KPart` test.
        let mut k_part = false;
        let mut pc1 = true;
        let mut axis_pos = None;
        let mut the_new_curve = None;
        let mut the_dir_extr = None;
        if s1.is_surface_of_linear_extrusion() && s2.gp_pln().is_some() {
            k_part = true;
            axis_pos = s2.gp_pln().map(|pl| pl.position());
            the_new_curve = s1.extrusion_basis_curve();
            the_dir_extr = s1.extrusion_direction();
        } else if s2.is_surface_of_linear_extrusion() && s1.gp_pln().is_some() {
            k_part = true;
            pc1 = false;
            axis_pos = s1.gp_pln().map(|pl| pl.position());
            the_new_curve = s2.extrusion_basis_curve();
            the_dir_extr = s2.extrusion_direction();
        }

        // `cxx:923-936`
        let mut a_circ = None;
        if k_part {
            let basis = the_new_curve.clone();
            match basis.as_ref().and_then(|b| b.gp_circ()) {
                Some(circ) => {
                    let axis = axis_pos.expect("KPart axis");
                    let ax_of_circ = circ.position().direction();
                    k_part = ax_of_circ.is_parallel_tol(&axis.direction(), Precision::ANGULAR);
                    if k_part {
                        a_circ = Some(circ);
                    }
                }
                None => k_part = false,
            }
        }

        if k_part {
            // `cxx:939-960`: `newC = TheNewCurve->Translated(VV)` plus its pcurve.
            let circ = a_circ.expect("KPart circle");
            let axis = axis_pos.expect("KPart axis");
            let dir_extr = the_dir_extr.expect("KPart direction");
            let to_axis = GpVec::from_pnts(&circ.location(), &axis.location());
            let a_local_real = to_axis.dot(&GpVec::from_xyz(dir_extr.xyz()));
            let cos = dir_extr.dot(&axis.direction());
            let vv = GpVec::new(dir_extr.x(), dir_extr.y(), dir_extr.z())
                .multiplied_scalar(a_local_real / cos);
            let Some(basis) = the_new_curve else {
                return self.fail_edge(the_edge);
            };
            let new_c: Arc<dyn Curve> = Arc::from(basis.translated(&vv));
            // `cxx:949-950`: `Geom2d_Line(gp_Pnt2d(0., aLocalReal / Cos), gp::DX2d())`.
            let dx = GpDir2d::new(1.0, 0.0).expect("gp::DX2d");
            let line = Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, a_local_real / cos), dx);
            let line: Arc<dyn Curve2d> = Arc::new(line);
            // `cxx:952-959`
            let ei = self.my_e_map.change_from_index(ii);
            if pc1 {
                *ei.change_first_pc() = Some(line);
            } else {
                *ei.change_second_pc() = Some(line);
            }
            return Ok(new_c);
        }

        // `cxx:961-1331`
        super::edges_intersection::intersect_edge_faces(self, the_edge, c, detrompeur, ii)
    }

    /// `cxx:852-889`: project the edge endpoints on both faces; `1` means both
    /// ends already sit on `S1`, `2` on `S2`, `0` neither.
    fn detrompeur(
        &self,
        the_edge: &Edge,
        c: &Arc<dyn Curve>,
        s1: &Arc<dyn Surface>,
        s2: &Arc<dyn Surface>,
        first_face: &Face,
        second_face: &Face,
    ) -> i32 {
        // `cxx:853-861`
        let (Some(fv), Some(lv)) = edge_vertices(the_edge) else {
            return 0;
        };
        let prmfv = BRepTool::parameter_of_vertex(&the_edge.0, &fv.0);
        let prmlv = BRepTool::parameter_of_vertex(&the_edge.0, &lv.0);
        let (pfv, d1fv) = c.d1(prmfv);
        let (plv, d1lv) = c.d1(prmlv);
        let _ = (d1fv, d1lv);

        // `cxx:863-871`
        let tol_f1 = BRepTool::face_tolerance(first_face);
        let tol_f2 = BRepTool::face_tolerance(second_face);
        let proj1 = project_point_on_surface(s1.as_ref(), &pfv, tol_f1);
        let proj2 = project_point_on_surface(s1.as_ref(), &plv, tol_f1);
        let proj3 = project_point_on_surface(s2.as_ref(), &pfv, tol_f2);
        let proj4 = project_point_on_surface(s2.as_ref(), &plv, tol_f2);

        // `cxx:873-890`
        let mut detrompeur = 0;
        if let (Some(a), Some(b)) = (&proj1, &proj2) {
            if a.distance <= CONFUSION && b.distance <= CONFUSION {
                detrompeur = 1;
            }
        }
        if let (Some(a), Some(b)) = (&proj3, &proj4) {
            if a.distance <= CONFUSION && b.distance <= CONFUSION {
                detrompeur = 2;
            }
        }
        detrompeur
    }

    /// `cxx:1333-1425`: the edge is tangent at `ptfixe`, so the new curve is the
    /// line through that point along the extrusion direction of a cylindrical
    /// or extruded face.
    fn recompute_tangent_edge(
        &mut self,
        the_edge: &Edge,
        c: &Arc<dyn Curve>,
        ptfixe: &GpPnt,
        first_face: Option<Face>,
        second_face: Option<Face>,
        ii: usize,
    ) -> Result<Arc<dyn Curve>, ()> {
        // `cxx:1334-1352`
        let (Some(first_face), Some(second_face)) = (first_face, second_face) else {
            return self.fail_edge(the_edge);
        };
        let a_local_s1 = self
            .my_f_map
            .find_from_key(&first_face.0)
            .and_then(|fi| fi.geometry())
            .cloned();
        let a_local_s2 = self
            .my_f_map
            .find_from_key(&second_face.0)
            .and_then(|fi| fi.geometry())
            .cloned();
        let (Some(a_local_s1), Some(a_local_s2)) = (a_local_s1, a_local_s2) else {
            return self.fail_edge(the_edge);
        };
        let a_local_s1 = a_local_s1
            .rectangular_trimmed_basis()
            .unwrap_or(a_local_s1);
        let a_local_s2 = a_local_s2
            .rectangular_trimmed_basis()
            .unwrap_or(a_local_s2);

        // `cxx:1353-1407`
        let mut dir_extr = GpDir::new(0.0, 0.0, 1.0).expect("gp_Dir default");
        let mut pc1: Option<Arc<dyn Curve2d>> = None;
        let mut pc2: Option<Arc<dyn Curve2d>> = None;
        if let Some(cy) = a_local_s1.gp_cylinder() {
            dir_extr = *cy.axis().direction();
        } else if let Some(dir) = a_local_s1.extrusion_direction() {
            dir_extr = dir;
            // `cxx:1370-1379`
            if let Some(basis) = a_local_s1.extrusion_basis_curve() {
                if let Some(circ) = basis.gp_circ() {
                    let u = clib::circle_parameter(&circ.position(), ptfixe);
                    let dy = GpDir2d::new(0.0, 1.0).expect("gp::DY2d");
                    pc1 = Some(Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(u, 0.0), dy)));
                }
            }
        } else if let Some(cy) = a_local_s2.gp_cylinder() {
            dir_extr = *cy.axis().direction();
        } else if let Some(dir) = a_local_s2.extrusion_direction() {
            dir_extr = dir;
            // `cxx:1393-1402`
            if let Some(basis) = a_local_s2.extrusion_basis_curve() {
                if let Some(circ) = basis.gp_circ() {
                    let u = clib::circle_parameter(&circ.position(), ptfixe);
                    let dy = GpDir2d::new(0.0, 1.0).expect("gp::DY2d");
                    pc2 = Some(Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(u, 0.0), dy)));
                }
            }
        }
        {
            let ei = self.my_e_map.change_from_index(ii);
            if pc1.is_some() {
                *ei.change_first_pc() = pc1;
            }
            if pc2.is_some() {
                *ei.change_second_pc() = pc2;
            }
        }

        // `cxx:1407-1408`
        let mut new_c: Arc<dyn Curve> =
            Arc::new(occt_geom::line::GeomLine::from_pnt_dir(*ptfixe, dir_extr));

        // `cxx:1410-1425`
        let (_, d1fv) = c.d1(0.0);
        let (_, new_d1) = new_c.d1(0.0);
        if d1fv.dot(&new_d1) < 0.0 {
            new_c = Arc::from(new_c.reversed());
            let ei = self.my_e_map.change_from_index(ii);
            if let Some(p) = ei.first_pc().cloned() {
                let rev: Arc<dyn Curve2d> = Arc::from(p.reversed());
                *ei.change_first_pc() = Some(rev);
            }
            if let Some(p) = ei.second_pc().cloned() {
                let rev: Arc<dyn Curve2d> = Arc::from(p.reversed());
                *ei.change_second_pc() = Some(rev);
            }
        }
        Ok(new_c)
    }

    pub(super) fn fail_edge(&mut self, e: &Edge) -> Result<Arc<dyn Curve>, ()> {
        self.set_edge_failure(e);
        Err(())
    }

    /// The `errStat` / `badShape` writes of a `Draft_EdgeRecomputation` exit.
    pub(super) fn set_edge_failure(&mut self, e: &Edge) {
        self.err_stat = DraftErrorStatus::EdgeRecomputation;
        self.bad_shape = Some(e.0.clone());
    }
}
