//! The entry points of `BRepFill_Evolved`:
//! the two `Perform` overloads (`BRepFill_Evolved.cxx:309-331`),
//! `PrivatePerform` (`cxx:334-547`), `SetWork` (`cxx:614-618`),
//! `FindLocation` (`cxx:2404-2437`) and `TransformInitWork`
//! (`cxx:2441-2445`).
//!
//! `PrivatePerform` is the reference's shared driver: it moves the spine and
//! the profile into the workspace of `AxeProf`, projects the profile and cuts
//! the spine, then builds one `Vevo` per profile wire on the left side and one
//! per spine wire on the right side, glues them with a `BRepTools_Quilt`, and
//! finally transfers the result back into the initial space.

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpPnt, GpTrsf};
use occt_core::kernel::geomabs::JoinType;
use occt_core::toploc::TopLocLocation;

use crate::abs::Orientation;
use crate::brep_lib_find_surface::FindSurface;
use crate::brep_lib_make_face::make_face_from_wire;
use crate::brep_mat2d_bisecting_locus::BRepMat2dBisectingLocus;
use crate::brep_mat2d_explorer::BRepMat2dExplorer;
use crate::brep_mat2d_link_topo_bilo::BRepMat2dLinkTopoBilo;
use crate::brep_sweep::Quilt;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::mat::MatSide;
use crate::shape::{Face, TopoShape, Wire};
use crate::topo_tools::wires_of_face;

use super::evolved::BRepFillEvolved;
use super::prepare_profile::prepare_profile;
use super::prepare_spine::prepare_spine;
use super::shape_map::ShapeMap;
use super::side::side;

impl BRepFillEvolved {
    /// `Perform(const TopoDS_Wire& Spine, ...)` (`cxx:309-318`).
    pub fn perform_wire(
        &mut self,
        spine: &Wire,
        profile: &Wire,
        axe_prof: &GpAx3,
        join: JoinType,
        solid: bool,
    ) -> Result<(), String> {
        // `cxx:312-315`.
        self.my_spine_type = false;
        let face = make_face_from_wire(spine, true)
            .ok_or_else(|| "BRepLib_MakeFace: the wire is not planar".to_string())?;
        self.private_perform(&face, profile, axe_prof, join, solid)
    }

    /// `Perform(const TopoDS_Face& Spine, ...)` (`cxx:322-331`).
    pub fn perform_face(
        &mut self,
        spine: &Face,
        profile: &Wire,
        axe_prof: &GpAx3,
        join: JoinType,
        solid: bool,
    ) -> Result<(), String> {
        // `cxx:328`.
        self.my_spine_type = true;
        self.private_perform(spine, profile, axe_prof, join, solid)
    }

    /// `PrivatePerform(Spine, Profile, AxeProf, Join, Solid)`
    /// (`cxx:334-547`).
    pub fn private_perform(
        &mut self,
        spine: &Face,
        profile: &Wire,
        axe_prof: &GpAx3,
        join: JoinType,
        solid: bool,
    ) -> Result<(), String> {
        // `cxx:340-348`.
        self.my_spine = Some(Face(spine.0.oriented(Orientation::Forward)));
        self.my_profile = Some(Wire(profile.0.oriented(Orientation::Forward)));
        self.my_join_type = join;
        self.my_map.clear();

        // `cxx:349-352`.
        if join != JoinType::Arc {
            return Err("Standard_NotImplemented".into());
        }

        // `cxx:353-368`: positioning of the spine and the profile.
        let my_spine = self.my_spine.clone().expect("spine was just set");
        let my_profile = self.my_profile.clone().expect("profile was just set");

        let l_spine = self.find_location(&my_spine)?;
        let mut t = GpTrsf::identity();
        t.set_transformation(axe_prof);
        let l_profile = TopLocLocation::composed(&t, TopLocLocation::identity());
        let init_ls = my_spine.0.location().clone();
        let init_lp = my_profile.0.location().clone();
        self.transform_init_work(&l_spine, &l_profile);

        // `cxx:372-375`.
        let mut work_prof: Vec<TopoShape> = Vec::new();
        let mut map_prof: ShapeMap<TopoShape> = ShapeMap::new();
        let mut work_spine = Face(TopoShape::new(crate::abs::ShapeType::Face));
        let mut map_spine: ShapeMap<TopoShape> = ShapeMap::new();
        prepare_profile(
            self.my_profile.as_ref().expect("profile was just set"),
            &mut work_prof,
            &mut map_prof,
        );
        prepare_spine(
            self.my_spine.as_ref().expect("spine was just set"),
            &mut work_spine,
            &mut map_spine,
        );

        // `cxx:377-396`: which sides the profile parts sit on.
        let tol = super::side::FILL_CONFUSION;
        let mut ya_left = false;
        let mut ya_right = false;
        for wp in &work_prof {
            let sp = Wire(wp.clone());
            if side(&sp, tol) < 4 {
                ya_left = true;
            } else {
                ya_right = true;
            }
            if ya_left && ya_right {
                break;
            }
        }

        // `cxx:399-420`: the cut volume and its working profile wire.
        let mut cut_vevo = BRepFillEvolved::new();
        let builder = TopoBuilder::new();
        let mut wp = Wire(builder.make_shape(crate::abs::ShapeType::Wire));
        for sp_shape in &work_prof {
            let mut exp = WireExplorer::new();
            exp.init(&Wire(sp_shape.clone()), None);
            while exp.more() {
                if let Some(e) = exp.current() {
                    builder.add(&mut wp.0, &e.0);
                }
                exp.next();
            }
        }
        cut_vevo.set_work(&work_spine, &wp);

        // `cxx:422`.
        let mut glue = Quilt::new();

        // `cxx:424-463`: the vevos on the left.
        if ya_left {
            let mut exp = BRepMat2dExplorer::new();
            exp.perform(&work_spine);
            let mut locus = BRepMat2dBisectingLocus::new();
            locus.compute(&mut exp, 1, MatSide::Left, JoinType::Arc, false);
            let mut link = BRepMat2dLinkTopoBilo::new();
            link.perform(&exp, &locus);

            for sp_shape in &work_prof {
                let sp = Wire(sp_shape.clone());
                let c_side = side(&sp, tol);
                let mut vevo = BRepFillEvolved::new();
                if c_side == 1 {
                    vevo.elementary_perform(&work_spine, &sp, &locus, &mut link)?;
                } else if c_side == 2 {
                    vevo.planar_perform(&work_spine, &sp, &locus, &mut link, join);
                } else if c_side == 3 {
                    vevo.vertical_perform(&work_spine, &sp, &locus, &mut link, join)?;
                }
                cut_vevo.add_with_glue(&vevo, &sp, &mut glue);
            }
        }

        // `cxx:465-507`: the vevos on the right, one per spine wire.
        if ya_right {
            for spine_wire in wires_of_face(&work_spine) {
                // `cxx:478-481`: the face of the plane `z = 0` bounded by the
                // reversed spine wire.
                let mut face = builder.make_face_plane(&super::side::plane_at_z(0.0));
                builder.add(
                    &mut face.0,
                    &spine_wire
                        .0
                        .oriented(super::side::complement(spine_wire.0.orientation())),
                );

                let mut exp = BRepMat2dExplorer::new();
                exp.perform(&face);
                let mut locus = BRepMat2dBisectingLocus::new();
                locus.compute(&mut exp, 1, MatSide::Left, JoinType::Arc, false);
                let mut link = BRepMat2dLinkTopoBilo::new();
                link.perform(&exp, &locus);

                for sp_shape in &work_prof {
                    let sp = Wire(sp_shape.clone());
                    let c_side = side(&sp, tol);
                    let mut vevo = BRepFillEvolved::new();
                    if c_side == 4 {
                        vevo.elementary_perform(&face, &sp, &locus, &mut link)?;
                    } else if c_side == 5 {
                        vevo.planar_perform(&face, &sp, &locus, &mut link, join);
                    } else if c_side == 6 {
                        vevo.vertical_perform(&face, &sp, &locus, &mut link, join)?;
                    }
                    cut_vevo.add_with_glue(&vevo, &sp, &mut glue);
                }
            }
        }

        // `cxx:509-512`.
        if solid {
            cut_vevo.add_top_and_bottom(&mut glue);
        }

        // `cxx:514-519`.
        cut_vevo.continuity_on_offset_edge();

        // `cxx:521-527`: the shape via the quilt.
        cut_vevo.my_shape = Some(glue.shells());

        // `cxx:528-538`: transfer the generated map and the shape, then
        // reposition in the initial space.
        self.transfert(
            &cut_vevo,
            &map_prof,
            &map_spine,
            &l_spine.inverted(),
            &init_ls,
            &init_lp,
        );

        // `cxx:536-538`.
        if solid {
            self.make_solid();
        }

        // `cxx:540-545`: `BRepLib::UpdateTolerances(myShape, false)`.
        if let Some(shape) = self.my_shape.clone() {
            crate::brep_lib_update_tolerances::update_tolerances(&shape, false)?;
        }
        self.my_is_done = true;
        Ok(())
    }

    /// `SetWork(Sp, Pr)` (`cxx:614-618`).
    pub fn set_work(&mut self, sp: &Face, pr: &Wire) {
        self.my_spine = Some(sp.clone());
        self.my_profile = Some(pr.clone());
    }

    /// `FindLocation(Face)` (`cxx:2404-2437`): the transform that carries the
    /// plane of `Face` onto the reference frame
    /// `gp_Ax3(gp_Pnt(0,0,0), gp::DZ(), gp::DX())`.
    pub fn find_location(&self, face: &Face) -> Result<TopLocLocation, String> {
        // `cxx:2406-2408`: `BRep_Tool::Surface(Face, L)`.
        let mut surface = BRepTool::face_surface(face);
        let mut loc = face.0.location().clone();
        // `if (!S->IsKind(STANDARD_TYPE(Geom_Plane)))` (`cxx:2410`): the exact
        // type test, so a `Geom_RectangularTrimmedSurface` of a plane takes the
        // `BRepLib_FindSurface` branch just as it does in the reference.
        let is_plain_plane = surface
            .as_ref()
            .is_some_and(|s| s.rectangular_trimmed_basis().is_none() && s.gp_pln().is_some());
        if !is_plain_plane {
            // `cxx:2411-2419`.
            let fs = FindSurface::new(&face.0, -1.0, true, false);
            if !fs.is_found() {
                return Err("BRepFill_Evolved : The Face is not planar".into());
            }
            surface = fs.surface();
            // `L = FS.Location()` (`cxx:2417`). The port answers that location
            // as the identity because it does not store representation
            // locations; OCCT's `myLocation` (`BRepLib_FindSurface.cxx:289`) is
            // `E.Location() * cr->Location()`, and the edges the reference
            // walks come from `TopExp_Explorer(Face, EDGE)`, so it carries the
            // face's own location. Keep `L = fs_location * face_location` so
            // the single `S->Transformed(L.Transformation())` below
            // (`cxx:2425-2427`) yields the same world plane.
            loc = TopLocLocation::composed(
                &fs.location().transformation(),
                face.0.location().clone(),
            );
        }
        if !loc.is_identity() {
            // `cxx:2425-2427`.
            surface = surface.map(|s| Arc::from(s.transformed(&loc.transformation())));
        }

        // `cxx:2429-2430`: `P = down_cast<Geom_Plane>(S); Axis = P->Position()`.
        let Some(plane) = surface.as_ref().and_then(|s| s.gp_pln()) else {
            return Err("BRepFill_Evolved : The Face is not planar".into());
        };
        let axis = plane.position();

        // `cxx:2432-2436`.
        let axe_ref = GpAx3::new(
            GpPnt::new(0.0, 0.0, 0.0),
            GpDir::new(0.0, 0.0, 1.0).expect("gp::DZ"),
            &GpDir::new(1.0, 0.0, 0.0).expect("gp::DX"),
        )
        .expect("gp_Ax3: Z and X are perpendicular");
        let mut t = GpTrsf::identity();
        t.set_transformation_from_to(&axe_ref, &axis);
        Ok(TopLocLocation::composed(&t, TopLocLocation::identity()))
    }

    /// `TransformInitWork(LS, LP)` (`cxx:2441-2445`).
    pub fn transform_init_work(&mut self, ls: &TopLocLocation, lp: &TopLocLocation) {
        if let Some(spine) = self.my_spine.as_mut() {
            spine.0.move_location(ls);
        }
        if let Some(profile) = self.my_profile.as_mut() {
            profile.0.move_location(lp);
        }
    }
}
