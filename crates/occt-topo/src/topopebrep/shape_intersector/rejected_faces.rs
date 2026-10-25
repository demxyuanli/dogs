//! `TopOpeBRep_ShapeIntersector::RejectedFaces` and the two file-local statics
//! it relies on.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_ShapeIntersector.cxx:187-366`
//! (`RejectedFaces`), `:1085-1236` (`OneShapeIsHalfSpace`), `:1244-1341`
//! (`GetNewSolid`).
//!
//! Translation notes:
//! - `aListOfShape` is a `NCollection_List<TopoDS_Shape>` iterated with
//!   `Iterator`; the `.Append` order is reproduced, and the two
//!   `Remove(it)` + `break` sweeps (`cxx:276-309`) become
//!   `position` + `remove`.
//! - `BRepAdaptor_Surface FSurf(testFace)` (`cxx:1148`) and
//!   `BRepAdaptor_Surface ASurf(hsFace)` (`cxx:1269`) are the
//!   `Restriction == true` constructor, i.e. `BRepTools::UVBounds` on the
//!   face's pcurves. [`SurfaceAdaptor::of_face`] is the `Restriction == false`
//!   arm (see `tool/adaptor.rs:48-81`, where that arm is UNPORTED); for the
//!   untrimmed sphere / torus / revolution faces these two tests are written
//!   for, the surface's own range is the same interval.
//! - `cxx:1277` reads `LastUParameter()` for `MaxV`; that is kept verbatim.
//! - `BRep_Builder::MakeShell` + `Add` is [`TopoBuilder::make_shell`];
//!   `newShell.Closed(BRep_Tool::IsClosed(newShell))` is
//!   [`TopoShape::set_closed`] with [`BRepTool::is_closed`].
//! - `BRepLib_MakeWire(e1, e2, e3, e4)` is [`MakeWire`] (the class that
//!   reverses edges and unifies vertices), not `BRep_Builder`-level
//!   `TopoBuilder::make_wire`.

use std::f64::consts::PI;

use crate::abs::{Orientation, ShapeType};
use crate::brep_lib_make_face::make_face_from_wire;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::shape::{Face, Solid, TopoShape};
use crate::topo_tools::{
    explorer_shapes, is_same_occurrence, map_shapes_and_ancestors,
};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::ShapeIntersector;

impl ShapeIntersector {
    /// `RejectedFaces(anObj, aReference, aListOfShape)` (`cxx:187-366`).
    ///
    /// Fills `a_list_of_shape` with the faces / edges of `an_obj` whose box
    /// touches none of the solid boxes of `a_reference`.
    pub fn rejected_faces(
        &mut self,
        an_obj: &TopoShape,
        a_reference: &TopoShape,
        a_list_of_shape: &mut Vec<TopoShape>,
    ) -> Result<(), String> {
        let is_half_space = one_shape_is_half_space(an_obj, a_reference);

        if is_half_space != 0 {
            a_list_of_shape.clear();

            let new_reject_face = if is_half_space == 1 {
                let (solid, face) = get_new_solid(an_obj)?;
                self.init(&solid.0, a_reference);

                self.reject_unboxed(a_reference, &solid.0, ShapeType::Face, a_list_of_shape)?;
                self.reject_unboxed(a_reference, &solid.0, ShapeType::Edge, a_list_of_shape)?;

                face
            } else {
                let (solid, face) = get_new_solid(a_reference)?;
                self.init(an_obj, &solid.0);

                self.reject_unboxed(&solid.0, an_obj, ShapeType::Face, a_list_of_shape)?;
                self.reject_unboxed(&solid.0, an_obj, ShapeType::Edge, a_list_of_shape)?;

                face
            };

            // Remove all shapes of `newRejectFace` from the list (`cxx:271-290`).
            for edgef in explorer_shapes(&new_reject_face.0, ShapeType::Edge) {
                let mut found = None;
                for (i, shape) in a_list_of_shape.iter().enumerate() {
                    if shape.shape_type() != ShapeType::Edge {
                        continue;
                    }
                    if is_same_occurrence(shape, &edgef) {
                        found = Some(i);
                        break;
                    }
                }
                if let Some(i) = found {
                    a_list_of_shape.remove(i);
                }
            }

            // Remove `newRejectFace` itself (`cxx:293-309`).
            let mut found = None;
            for (i, shape) in a_list_of_shape.iter().enumerate() {
                if shape.shape_type() != ShapeType::Face {
                    continue;
                }
                if is_same_occurrence(shape, &new_reject_face.0) {
                    found = Some(i);
                    break;
                }
            }
            if let Some(i) = found {
                a_list_of_shape.remove(i);
            }

            self.init(an_obj, a_reference);
            return Ok(());
        }

        self.init(an_obj, a_reference);

        a_list_of_shape.clear();
        // Find faces to reject (`cxx:319-350`).
        self.reject_unboxed(a_reference, an_obj, ShapeType::Face, a_list_of_shape)?;
        self.reject_unboxed(a_reference, an_obj, ShapeType::Edge, a_list_of_shape)?;

        Ok(())
    }

    /// The `Compare(Current()).More()` loop repeated at `cxx:209-266` and
    /// `:326-350`: appends every `texplo` occurrence of `explore` whose box
    /// touches no scanned box of `scan`.
    ///
    /// `tscann` is `TopAbs_SOLID` at every call site (`cxx:202`, `:238`,
    /// `:320`).
    fn reject_unboxed(
        &mut self,
        scan: &TopoShape,
        explore: &TopoShape,
        texplo: ShapeType,
        out: &mut Vec<TopoShape>,
    ) -> Result<(), String> {
        self.face_scanner.clear();
        self.face_scanner
            .add_boxes_make_cob(scan, ShapeType::Solid, ShapeType::Shape)?;
        self.face_explorer.init(explore, texplo, ShapeType::Shape);

        while self.face_explorer.more() {
            let cur = match self.face_explorer.current() {
                Some(c) => c.clone(),
                None => break,
            };
            let touched = {
                let bs = self.face_scanner.change_box_sort();
                bs.compare(&cur)?;
                bs.more()
            };
            if !touched {
                out.push(cur);
            }
            self.face_explorer.next();
        }

        Ok(())
    }
}

/// `OneShapeIsHalfSpace(S1, S2)` (`cxx:1085-1236`).
///
/// Returns `0` (no half space), `1` (the half space is `S1`) or `2` (it is
/// `S2`).
fn one_shape_is_half_space(s1: &TopoShape, s2: &TopoShape) -> i32 {
    let mut result = 0;

    if s1.shape_type() != ShapeType::Solid || s2.shape_type() != ShapeType::Solid {
        return result;
    }

    let faces1 = explorer_shapes(s1, ShapeType::Face);
    let faces2 = explorer_shapes(s2, ShapeType::Face);
    let nb_faces_sol1 = faces1.len() as i32;
    let nb_faces_sol2 = faces2.len() as i32;

    if nb_faces_sol1 == 0 || nb_faces_sol2 == 0 {
        // Strange solids.
        return result;
    }
    if nb_faces_sol1 == 1 && nb_faces_sol2 == 1 {
        // Both shapes are half spaces.
        return result;
    }

    if !((nb_faces_sol1 == 1 && nb_faces_sol2 >= 2)
        || (nb_faces_sol2 == 1 && nb_faces_sol1 >= 2))
    {
        // `cxx:1218-1227`: both shells have more than two faces. The comment
        // there states the case is not classified; nothing is done.
        return result;
    }

    // If one solid has a shell made of a single face while the other has a
    // valid closed shell, the operation is a boolean with a half-space object.
    // Spheres, toruses and solids whose single face is built on a spherical
    // surface of revolution are not half spaces.
    let test_face = if nb_faces_sol1 == 1 {
        faces1.last().cloned()
    } else {
        faces2.last().cloned()
    };

    if let Some(test_face) = test_face.and_then(Face::wrap) {
        let f_surf = SurfaceAdaptor::of_face(&test_face);
        let mut solid_is_sphere_or_torus = false;

        if f_surf.get_type() == GeomAbsSurfaceType::Sphere
            || f_surf.get_type() == GeomAbsSurfaceType::Torus
        {
            let min_u = f_surf.first_u_parameter();
            let max_u = f_surf.last_u_parameter();
            let min_v = f_surf.first_v_parameter();
            let max_v = f_surf.last_v_parameter();
            let yes_u = (min_u - 0.0).abs() < 1.0e-9 && (max_u - 2.0 * PI).abs() < 1.0e-9;
            let yes_v = if f_surf.get_type() == GeomAbsSurfaceType::Sphere {
                (min_v - (-PI / 2.0)).abs() < 1.0e-9 && (max_v - PI / 2.0).abs() < 1.0e-9
            } else {
                (min_v - 0.0).abs() < 1.0e-9 && (max_v - 2.0 * PI).abs() < 1.0e-9
            };
            solid_is_sphere_or_torus = yes_u && yes_v;
        }

        if f_surf.get_type() == GeomAbsSurfaceType::SurfaceOfRevolution {
            if let Some(surf) = f_surf.surface() {
                let are_both_periodic = surf.is_u_periodic() && surf.is_v_periodic();
                if are_both_periodic {
                    let u_period = surf.u_period();
                    let v_period = surf.v_period();
                    let yes_u = (u_period - PI).abs() < 1.0e-9 || (u_period - 2.0 * PI).abs() < 1.0e-9;
                    let yes_v = (v_period - PI).abs() < 1.0e-9 || (v_period - 2.0 * PI).abs() < 1.0e-9;
                    solid_is_sphere_or_torus = yes_u && yes_v;
                }
            }
        }

        if solid_is_sphere_or_torus {
            return result;
        }

        // `SecondShellOk`: every edge of the valid shell must belong to exactly
        // two faces.
        let other = if nb_faces_sol1 == 1 { s2 } else { s1 };
        let a_map_ef = map_shapes_and_ancestors(other, ShapeType::Edge, ShapeType::Face);
        let mut second_shell_ok = true;
        for list_faces in a_map_ef.values() {
            let nb_faces = list_faces.len() as i32;
            if nb_faces != 2 {
                second_shell_ok = false;
                break;
            }
        }

        if second_shell_ok {
            result = if nb_faces_sol1 == 1 { 1 } else { 2 };
        }
    }

    result
}

/// `GetNewSolid(S, F)` (`cxx:1244-1341`).
///
/// Rebuilds the half-space solid with an extra "face on infinity", and returns
/// it together with that face — the caller removes the face (and its edges)
/// from the rejected list.
fn get_new_solid(s: &TopoShape) -> Result<(Solid, Face), String> {
    let hs_face = explorer_shapes(s, ShapeType::Face)
        .last()
        .cloned()
        .and_then(Face::wrap)
        .ok_or_else(|| "GetNewSolid : no face".to_string())?;

    let a_surf = SurfaceAdaptor::of_face(&hs_face);

    let min_u = a_surf.first_u_parameter();
    let max_u = a_surf.last_u_parameter();
    let min_v = a_surf.first_v_parameter();
    // `cxx:1277`: `double MaxV = ASurf.LastUParameter();`
    let max_v = a_surf.last_u_parameter();

    let mid_u = (max_u + min_u) * 0.5;
    let mid_v = (max_v + min_v) * 0.5;

    let surf = a_surf
        .surface()
        .ok_or_else(|| "GetNewSolid : no surface".to_string())?;
    let (_mid_p, surf_du, surf_dv) = surf.d1(mid_u, mid_v);

    let normal = surf_du.crossed(&surf_dv);
    let normal = if hs_face.0.orientation() == Orientation::Forward {
        normal.multiplied_scalar(-1.0e10)
    } else {
        normal.multiplied_scalar(1.0e10)
    };

    let pu1 = min_u + ((max_u - min_u) / 4.0).abs();
    let pu2 = min_u + ((max_u - min_u) / 4.0 * 3.0).abs();
    let pv1 = min_v + ((max_v - min_v) / 4.0).abs();
    let pv2 = min_v + ((max_v - min_v) / 4.0 * 3.0).abs();

    let p1 = surf.d0(pu1, pv1).translated_vec(&normal);
    let p2 = surf.d0(pu1, pv2).translated_vec(&normal);
    let p3 = surf.d0(pu2, pv1).translated_vec(&normal);
    let p4 = surf.d0(pu2, pv2).translated_vec(&normal);

    let b = TopoBuilder::new();
    let e1 = b.make_edge_segment(&p1, &p2);
    let e2 = b.make_edge_segment(&p2, &p4);
    let e3 = b.make_edge_segment(&p4, &p3);
    let e4 = b.make_edge_segment(&p3, &p1);

    let mut mkw = MakeWire::new();
    for e in [&e1, &e2, &e3, &e4] {
        mkw.add(e)
            .map_err(|_| "GetNewSolid : BRepLib_MakeWire".to_string())?;
    }
    let w = mkw.wire();

    let inf_face = make_face_from_wire(&w, false)
        .ok_or_else(|| "GetNewSolid : BRepLib_MakeFace".to_string())?;

    let new_shell = b.make_shell(&[hs_face, inf_face.clone()]);
    new_shell.0.set_closed(BRepTool::is_closed(&new_shell.0));

    let new_solid = b.make_solid(&[new_shell]);

    Ok((new_solid, inf_face))
}
