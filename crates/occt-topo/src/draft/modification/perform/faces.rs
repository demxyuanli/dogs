//! `Draft_Modification::Perform`, first block: "Calculate eventual faces"
//! (`Draft_Modification_1.cxx:724-824`).
//!
//! Every face record that `InternalAdd` left with `NewGeom == true` and a null
//! geometry is a postponed cylinder / extrusion; the two faces that share it
//! are planes, and their intersection line is the extrusion direction. The face
//! is the circle `InternalAdd` kept (`Finf.Curve()`) swept along that line,
//! trimmed to one and a bit turns.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | face loop, plane / plane direction | `_1.cxx:724-798` | [`DraftModification::perform_faces`] |
//! | extrusion surface + `UVBounds` widening | `_1.cxx:800-824` | [`DraftModification::perform_faces`] |
//!
//! # Translation notes
//!
//! * `IntAna_QuadQuadGeo(pp1, pp2, Angular, Confusion)` (`cxx:767`) is
//!   [`quadric_quadric_planes`]; every arm other than `IntAna_Line` takes the
//!   `Draft_FaceRecomputation` exit (`cxx:768-772`).
//! * `down_cast<Geom_Plane>` (`cxx:757-758`) is [`Surface::gp_pln`], and the
//!   `Geom_RectangularTrimmedSurface` unwrapping (`cxx:749-756`) is
//!   [`Surface::rectangular_trimmed_basis`].
//! * `Finf.Curve()` is down-cast to `Geom_Circle` in the reference
//!   (`cxx:803`); a record without a circle would build a null basis surface,
//!   so the port reports `Draft_FaceRecomputation` instead (documented
//!   deviation, no OCCT rule is added).

use std::sync::Arc;

use occt_core::gp::GpDir;
use occt_core::precision::{Precision, ANGULAR, CONFUSION};
use occt_geom::intana::{quadric_quadric_planes, QuadricIntersection};
use occt_geom::rectangular_trimmed::GeomRectangularTrimmedSurface;
use occt_geom::surface_of_linear_extrusion::GeomSurfaceOfLinearExtrusion;
use occt_geom::Surface;

use crate::brep_tool::BRepTool;
use crate::brep_tools::uv_bounds;
use crate::draft::error_status::DraftErrorStatus;
use crate::shape::{Face, TopoShape};

use super::DraftModification;

impl DraftModification {
    /// `Draft_Modification_1.cxx:724-824`.
    pub(super) fn perform_faces(&mut self) -> Result<(), ()> {
        for i in 1..=self.my_f_map.extent() {
            let fk = self.my_f_map.find_key(i).clone();
            // `cxx:728-730`
            let (new_geometry, geometry_null, first_face, second_face) = {
                let finf = self.my_f_map.change_from_index(i);
                (
                    finf.new_geometry(),
                    finf.geometry().is_none(),
                    finf.first_face().cloned(),
                    finf.second_face().cloned(),
                )
            };
            if !(new_geometry && geometry_null) {
                continue;
            }

            // `cxx:733-738`
            let (Some(f1), Some(f2)) = (first_face, second_face) else {
                return self.fail_face(&fk);
            };

            // `cxx:741-748`
            let s1 = self
                .my_f_map
                .find_from_key(&f1.0)
                .and_then(|fi| fi.geometry())
                .cloned();
            let s2 = self
                .my_f_map
                .find_from_key(&f2.0)
                .and_then(|fi| fi.geometry())
                .cloned();
            let (Some(s1), Some(s2)) = (s1, s2) else {
                return self.fail_face(&fk);
            };

            // `cxx:749-756`
            let s1 = s1.rectangular_trimmed_basis().unwrap_or(s1);
            let s2 = s2.rectangular_trimmed_basis().unwrap_or(s2);

            // `cxx:757-761`
            let (Some(p1), Some(p2)) = (s1.gp_pln(), s2.gp_pln()) else {
                return self.fail_face(&fk);
            };

            // `cxx:765-773`
            let QuadricIntersection::Line(line) =
                quadric_quadric_planes(&p1, &p2, ANGULAR, CONFUSION)
            else {
                return self.fail_face(&fk);
            };
            let mut extrdir = line.direction();

            // `cxx:777-797`: preserve the direction of the base face.
            let Some(k_face) = Face::wrap(fk.clone()) else {
                return self.fail_face(&fk);
            };
            let ref_surf = BRepTool::face_surface(&k_face)
                .map(|s| s.rectangular_trimmed_basis().unwrap_or(s));
            let mut dir_ref = GpDir::new(0.0, 0.0, 1.0).expect("GpDir::default");
            if let Some(ref_surf) = &ref_surf {
                if let Some(cy) = ref_surf.gp_cylinder() {
                    dir_ref = cy.position().direction();
                } else if let Some(d) = ref_surf.extrusion_direction() {
                    dir_ref = d;
                }
            }
            // `cxx:796-798`
            if extrdir.dot(&dir_ref) < 0.0 {
                extrdir.reverse();
            }

            // `cxx:802-804`: `new Geom_SurfaceOfLinearExtrusion(CCir, extrdir)`.
            let Some(circle) = self
                .my_f_map
                .find_from_index(i)
                .curve()
                .cloned()
            else {
                return self.fail_face(&fk);
            };
            if circle.gp_circ().is_none() {
                return self.fail_face(&fk);
            }
            let new_s: Arc<dyn Surface> =
                Arc::new(GeomSurfaceOfLinearExtrusion::new(circle, extrdir));

            // `cxx:806-821`
            let (_, _, vmin, vmax) = uv_bounds(&k_face);
            let (vmin, vmax) = if !Precision::is_negative_infinite(vmin)
                && !Precision::is_positive_infinite(vmax)
            {
                let deltav = 2.0 * (vmax - vmin);
                (vmin - deltav, vmax + deltav)
            } else {
                (-300.0, 300.0)
            };

            // `cxx:823`: `new Geom_RectangularTrimmedSurface(NewS, 0., 1.9*PI, vmin, vmax)`.
            let trimmed = GeomRectangularTrimmedSurface::uv(
                new_s,
                0.0,
                1.9 * std::f64::consts::PI,
                vmin,
                vmax,
            );
            let Ok(trimmed) = trimmed else {
                return self.fail_face(&fk);
            };
            let trimmed: Arc<dyn Surface> = Arc::new(trimmed);
            *self.my_f_map.change_from_index(i).change_geometry() = Some(trimmed);
        }
        Ok(())
    }

    /// `cxx:735-737`, `:743-747`, `:759-761`, `:768-771`.
    fn fail_face(&mut self, fk: &TopoShape) -> Result<(), ()> {
        self.err_stat = DraftErrorStatus::FaceRecomputation;
        self.bad_shape = Some(fk.clone());
        Err(())
    }
}
