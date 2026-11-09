//! `Draft_Modification::InternalAdd` (`Draft_Modification_1.cxx:116-420`): the
//! per-face step that records the new surface of `F` and the new 3D curves of
//! the edges whose second face is tangent or absent, recursing into that
//! neighbour.
//!
//! # Translation notes
//!
//! * `BRep_Tool::Surface(F, Lo)` followed by `Transformed(Lo.Transformation())`
//!   is [`BRepTool::face_surface_world`]; the `Geom_RectangularTrimmedSurface`
//!   unwrap of `cxx:134-137` is [`Surface::rectangular_trimmed_basis`], and like
//!   the reference `down_cast` it leaves the surface untouched when the dynamic
//!   type is something else.
//! * The postponed arm (`cxx:141-198`) needs `Flag == false`, which is only
//!   reachable from `Perform`; it records the neutral-plane circle instead of a
//!   surface, and `Propagate` later turns that circle into an extrusion.
//! * `EInf.Add(F)` with a null `OtherF` (`cxx:328`) is a no-op in
//!   `Draft_FaceInfo::Add` (`cxx:52-62` tests `IsNull()` first), so the port
//!   only calls [`DraftFaceInfo::add`] when the neighbour exists.
//! * `Standard_Type` comparisons against `Geom_CylindricalSurface` /
//!   `Geom_ConicalSurface` / `Geom_SurfaceOfLinearExtrusion`
//!   (`cxx:143-146`, `cxx:211-213`, `cxx:396-398`) become the
//!   [`Surface`] type accessors `gp_cylinder`, `gp_cone` and
//!   `is_surface_of_linear_extrusion`, whose defaults answer exactly those
//!   dynamic types.

use std::sync::Arc;

use occt_core::elib::slib;
use occt_core::gp::{GpDir, GpPln, GpVec};
use occt_core::precision::{Precision, ANGULAR, CONFUSION};
use occt_geom::circle::GeomCircle;
use occt_geom::intana::IntAnaIntConicQuad;
use occt_geom::rectangular_trimmed::GeomRectangularTrimmedSurface;
use occt_geom::Curve;
use occt_geom::Surface;

use crate::abs::ShapeType;
use crate::brep_lib_continuity::stored_continuity;
use crate::brep_tool::BRepTool;
use crate::brep_tools::{is_really_closed, uv_bounds};
use crate::draft::edge_info::DraftEdgeInfo;
use crate::draft::face_info::DraftFaceInfo;
use crate::draft::error_status::DraftErrorStatus;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same};

use super::new_curve::new_curve;
use super::new_surface::new_surface;
use super::{DraftModification, shape_face_orientation};

impl DraftModification {
    /// `InternalAdd(F, Direction, Angle, NeutralPlane, Flag)`
    /// (`Draft_Modification_1.cxx:116-420`).
    pub(super) fn internal_add(
        &mut self,
        f: &Face,
        direction: &GpDir,
        angle: f64,
        neutral_plane: &GpPln,
        flag: bool,
    ) -> bool {
        // `cxx:123-126`
        if self.my_f_map.contains(&f.0) {
            return self.bad_shape.is_none();
        }

        // `cxx:128`
        let oris = shape_face_orientation(&self.my_shape, f);

        // `cxx:131-137`
        let Some(mut s) = BRepTool::face_surface_world(f) else {
            // `BRep_Tool::Surface` answered a null handle, which the reference
            // dereferences at `cxx:133`; the port reports the face instead.
            self.bad_shape = Some(f.0.clone());
            self.err_stat = DraftErrorStatus::FaceRecomputation;
            return false;
        };
        if let Some(basis) = s.rectangular_trimmed_basis() {
            s = basis;
        }

        // `cxx:138-139`
        let mut new_s: Option<Arc<dyn Surface>> = None;
        let mut the_circle: Option<Arc<dyn Curve>> = None;

        // `cxx:141-198`
        let mut postponed = !flag;
        if postponed {
            if let Some(cyl) = s.gp_cylinder() {
                // `cxx:148-157`
                let axcyl = cyl.axis();
                let mut cir = slib::cylinder_v_iso(&cyl.position(), cyl.radius(), 0.0);
                let vv = GpVec::from_pnts(&cyl.location(), &neutral_plane.location());
                let dir_v = GpVec::from_xyz(axcyl.direction().xyz());
                cir.translate_vec(&dir_v.multiplied_scalar(vv.dot(&dir_v)));
                the_circle = Some(Arc::new(GeomCircle::new(cir)));
            } else if s.is_surface_of_linear_extrusion() {
                // `cxx:158-190`
                let Some(mut basis) = s.extrusion_basis_curve() else {
                    self.bad_shape = Some(f.0.clone());
                    self.err_stat = DraftErrorStatus::FaceRecomputation;
                    return false;
                };
                let Some(the_dirextr) = s.extrusion_direction() else {
                    self.bad_shape = Some(f.0.clone());
                    self.err_stat = DraftErrorStatus::FaceRecomputation;
                    return false;
                };
                if basis.is_geom_trimmed() {
                    if let Some(b) = basis.trim_basis() {
                        basis = b;
                    }
                }
                let Some(circ_basis) = basis.gp_circ() else {
                    self.bad_shape = Some(f.0.clone());
                    self.err_stat = DraftErrorStatus::FaceRecomputation;
                    return false;
                };
                let mut cir = circ_basis;
                let circ_axis = cir.axis();
                let dircir = circ_axis.direction();
                // `cxx:170-175`
                if !direction.is_parallel_tol(&dircir, ANGULAR) {
                    self.bad_shape = Some(f.0.clone());
                    self.err_stat = DraftErrorStatus::FaceRecomputation;
                    return false;
                }
                // `cxx:185-190`
                let axis = neutral_plane.position();
                let axis_loc = axis.location();
                let l = GpVec::from_pnts(&cir.location(), &axis_loc)
                    .dot(&GpVec::from_xyz(axis.direction().xyz()));
                let cos = the_dirextr.dot(&axis.direction());
                let vv = GpVec::from_xyz(the_dirextr.xyz()).multiplied_scalar(l / cos);
                cir.translate_vec(&vv);
                the_circle = Some(Arc::new(GeomCircle::new(cir)));
            } else {
                // `cxx:193-196`
                postponed = false;
            }
        }

        // `cxx:200-262`
        if !postponed {
            new_s = new_surface(&s, oris, direction, angle, neutral_plane);
            let Some(ns) = new_s.clone() else {
                self.bad_shape = Some(f.0.clone());
                self.err_stat = DraftErrorStatus::FaceRecomputation;
                return false;
            };
            // `cxx:210-260`: keep enough of the cylinder / cone domain to survive
            // the intersection with the neutral plane.
            let is_cyl = ns.gp_cylinder().is_some();
            let is_cone = ns.gp_cone().is_some();
            if is_cyl || is_cone {
                let (_, _, vmin, vmax) = uv_bounds(f);
                if !Precision::is_negative_infinite(vmin) && !Precision::is_positive_infinite(vmax) {
                    let deltav = 10.0 * (vmax - vmin);
                    let (mut vmin, mut vmax) = (vmin, vmax);
                    if is_cyl {
                        // `cxx:219-222`
                        vmin -= deltav;
                        vmax += deltav;
                    } else {
                        // `cxx:224-258`
                        let co = ns.gp_cone().unwrap();
                        let vapex = -co.radius() / co.semi_angle().sin();
                        if vmin < vapex {
                            if vmax + deltav > vapex {
                                vmax = vapex;
                                vmin -= 10.0 * (vmax - vmin);
                                vmax -= CONFUSION;
                            } else {
                                vmin -= deltav;
                                vmax += deltav;
                            }
                        } else if vmin - deltav < vapex {
                            vmin = vapex;
                            vmax += 10.0 * (vmax - vmin);
                            vmin += CONFUSION;
                        } else {
                            vmin -= deltav;
                            vmax += deltav;
                        }
                    }
                    // `cxx:259`
                    if let Ok(trimmed) = GeomRectangularTrimmedSurface::uv(
                        ns.clone(),
                        0.0,
                        2.0 * std::f64::consts::PI,
                        vmin,
                        vmax,
                    ) {
                        new_s = Some(Arc::new(trimmed));
                    }
                }
            }
        }

        // `cxx:264-274`
        let keep = new_s
            .as_ref()
            .map(|ns| !Arc::ptr_eq(&s, ns))
            .unwrap_or(true);
        if keep {
            let mut fi = DraftFaceInfo::with_surface(new_s.as_ref(), true);
            fi.root_face(&self.cur_face);
            self.my_f_map.add(&f.0, fi);
            if postponed {
                if let Some(fi) = self.my_f_map.change_from_key(&f.0) {
                    *fi.change_curve() = the_circle;
                }
            }
        }

        // `cxx:276-418`
        let mut map_of_e: Vec<TopoShape> = Vec::new();
        for edg in explorer_shapes(&f.0, ShapeType::Edge) {
            if self.bad_shape.is_some() {
                break;
            }
            let edg = Edge(edg);
            if self.my_e_map.contains(&edg.0) {
                continue;
            }
            let mut addedg = false;
            let mut addface = false;
            let mut other_f: Option<Face> = None;
            // `cxx:286-288`
            if is_really_closed(&edg, f) {
                addedg = true;
            } else {
                // `cxx:292-309`: find the other face(s) containing the edge.
                let mut nbother = 0;
                // `cxx:295`: `myEFMap.FindFromKey(edg)` throws
                // `Standard_NoSuchObject` for an edge outside the map, which
                // only happens when `Add` follows `Perform` (`cxx:98-101`
                // clears the map without rebuilding it).
                let ancestors = self
                    .my_ef_map
                    .find_from_key(&edg.0)
                    .expect("Standard_NoSuchObject: NCollection_IndexedDataMap::FindFromKey");
                for a in ancestors {
                    if !is_same(a, &f.0) {
                        if other_f.is_none() {
                            other_f = Some(Face(a.clone()));
                        }
                        nbother += 1;
                    }
                }
                if nbother >= 2 {
                    // `cxx:310-314`
                    self.bad_shape = Some(edg.0.clone());
                    self.err_stat = DraftErrorStatus::EdgeRecomputation;
                } else if let Some(of) = &other_f {
                    // `cxx:315-319`
                    if stored_continuity(&edg, f, of) >= 1 {
                        addface = true;
                        addedg = true;
                    }
                }
                // `cxx:319-322`: `nbother == 0` sets nothing (commented out).
            }

            if addedg {
                // `cxx:326-331`
                if postponed {
                    if let (Some(of), Some(fi)) = (&other_f, self.my_f_map.change_from_key(&f.0)) {
                        fi.add(of);
                    }
                }
                // `cxx:332-338`
                let Some(next_c) = BRepTool::edge_curve_world(&edg) else {
                    self.bad_shape = Some(edg.0.clone());
                    self.err_stat = DraftErrorStatus::EdgeRecomputation;
                    break;
                };
                let mut c = next_c;
                if c.is_geom_trimmed() {
                    if let Some(b) = c.trim_basis() {
                        c = b;
                    }
                }

                // `cxx:339-375`
                let mut einf = DraftEdgeInfo::with_new_geometry(true);
                let mut new_c: Option<Arc<dyn Curve>> = None;
                if postponed {
                    // `cxx:340-367`
                    einf.add(f);
                    if let Some(of) = &other_f {
                        einf.add(of);
                    }
                    let Some(local_geom) = c.gp_line() else {
                        self.bad_shape = Some(edg.0.clone());
                        self.err_stat = DraftErrorStatus::EdgeRecomputation;
                        break;
                    };
                    let ilipl =
                        IntAnaIntConicQuad::line_plane(&local_geom, neutral_plane, ANGULAR, 0.0, 0.0);
                    if ilipl.is_done() && ilipl.nb_points() != 0 {
                        einf.tangent(ilipl.point(1));
                    } else {
                        self.bad_shape = Some(edg.0.clone());
                        self.err_stat = DraftErrorStatus::EdgeRecomputation;
                        break;
                    }
                } else {
                    // `cxx:368-375`
                    match new_curve(&c, &s, oris, direction, angle, neutral_plane) {
                        Some(nc) => new_c = Some(nc),
                        None => {
                            self.bad_shape = Some(edg.0.clone());
                            self.err_stat = DraftErrorStatus::EdgeRecomputation;
                            break;
                        }
                    }
                }

                // `cxx:378-383`
                if let Some(nc) = &new_c {
                    if nc.is_geom_trimmed() {
                        if let Some(b) = nc.trim_basis() {
                            new_c = Some(b);
                        }
                    }
                }

                // `cxx:384-387`
                *einf.change_geometry() = new_c;
                einf.root_face(&self.cur_face);
                self.my_e_map.add(&edg.0, einf);
                map_of_e.push(edg.0.clone());

                // `cxx:388-413`
                if addface {
                    let fl = flag;
                    let Some(of) = other_f.clone() else {
                        break;
                    };
                    let alocal = BRepTool::face_surface(&of);
                    if let Some(alocal) = alocal {
                        let typ = alocal.rectangular_trimmed_basis().unwrap_or(alocal);
                        // `cxx:395-398`
                        let postponed_kind = typ.gp_cylinder().is_some()
                            || typ.is_surface_of_linear_extrusion();
                        if postponed_kind && self.my_f_map.contains(&f.0) && !flag && !postponed {
                            // `cxx:400-411`
                            self.my_f_map.remove_key(&f.0);
                            for e in &map_of_e {
                                self.my_e_map.remove_key(e);
                            }
                        }
                    }
                    // `cxx:412`
                    self.internal_add(&of, direction, angle, neutral_plane, fl);
                }
            }
        }

        // `cxx:419`
        self.bad_shape.is_none()
    }
}
