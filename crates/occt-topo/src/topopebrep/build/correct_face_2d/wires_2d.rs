//! `TopOpeBRepBuild_CorrectFace2d` — the 2d placement half.
//!
//! Sources (OCCT 8.0.0, `TopOpeBRepBuild_CorrectFace2d.cxx`):
//! `TranslateCurve2d` `:851-862`, `OuterWire` `:866-890`, `BndBoxWire`
//! `:894-907`, `MoveWire2d` `:911-977`, `MoveWires2d` `:982-1030` and the two
//! `UpdateEdge` overloads `:1033-1117`.
//!
//! `MoveWires2d` first shifts the outer wire into `[0, 2*PI]` and then shifts
//! every other wire's `u` range by whole periods until it falls inside the
//! outer wire's 2d box. The two `UpdateEdge` overloads are the "false builder"
//! of `cxx:539` / `:554`: they write the p-curve on the work copy *and* on the
//! original edge recorded in `myMap`, creating that original-side edge with
//! `EmptyCopied` + its vertices the first time it is needed.

use std::sync::Arc;

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpTrsf2d, GpVec2d};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::{Orientation, ShapeType};
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::{FClass2d, FaceState};
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools::{explorer_shapes, is_equal};
use crate::topopebrep::build::tools::write_edge_pcurve;
use crate::topopebrep::build::tools::write_edge_pcurves;
use crate::topopebrep::ds::data_structure::IndexedShapeSet;

use super::CorrectFace2d;

const TWO_PI: f64 = 2.0 * std::f64::consts::PI;

impl CorrectFace2d {
    /// `TranslateCurve2d(anE, aFace, aTrV, aC2DOut)` (`cxx:851-862`):
    /// `C2d->Copy()`, trimmed to the `CurveOnSurface` range, then translated.
    pub(super) fn translate_curve2d(
        &self,
        an_e: &TopoShape,
        a_face: &Face,
        a_tr_v: &GpVec2d,
    ) -> Option<Arc<dyn Curve2d>> {
        let e = Edge::wrap(an_e.clone())?;
        let (c2d, a_first, a_last) = crate::boptools_2d::curve_on_surface_range(&e, a_face)?;
        // `C2d->Copy()` has no port; `transformed(identity)` is the same copy.
        let a_tr_c: Arc<dyn Curve2d> = Arc::from(c2d.transformed(&GpTrsf2d::identity()));
        let trimmed = Geom2dTrimmedCurve::new(a_tr_c, a_first, a_last);
        let mut tr = GpTrsf2d::identity();
        tr.set_translation_vec(a_tr_v);
        let out: Box<dyn Curve2d> = trimmed.transformed(&tr);
        Some(Arc::from(out))
    }

    /// `OuterWire(anOuterWire)` (`cxx:866-890`): the first wire whose infinite
    /// point classifies `TopAbs_OUT` on a face built from that wire alone.
    pub(super) fn outer_wire(&self) -> Option<TopoShape> {
        let surf = BRepTool::face_surface_world(&self.my_corrected_face)?;
        let builder = TopoBuilder::new();

        for wire in explorer_shapes(&self.my_corrected_face.0, ShapeType::Wire) {
            let mut new_face =
                builder.make_face_with_tolerance(surf.clone(), &[], self.my_face_tolerance);
            let Some(a_wire) = Wire::wrap(wire.clone()) else {
                continue;
            };
            builder.add_wire(&mut new_face, &a_wire);

            let face = new_face;
            let Ok(a_class2d) = FClass2d::new(&face, self.my_face_tolerance) else {
                continue;
            };
            if a_class2d.perform_infinite_point() == FaceState::Out {
                return Some(wire);
            }
        }
        // Such a wire is not found.
        None
    }

    /// `BndBoxWire(aWire, B2d)` (`cxx:894-907`): the 2d box of every edge's
    /// p-curve on `myCorrectedFace`, enlarged by the edge tolerance.
    pub(super) fn bnd_box_wire(&self, a_wire: &TopoShape) -> BndBox2d {
        let mut a_b2d = BndBox2d::new();

        for an_edge in explorer_shapes(a_wire, ShapeType::Edge) {
            let Some(edge) = Edge::wrap(an_edge) else {
                continue;
            };
            let a_tol_e = BRepTool::edge_tolerance(&edge);
            // `BRepAdaptor_Curve2d aBAC2d(anEdge, myCorrectedFace)` +
            // `BndLib_Add2dCurve::Add(aBAC2d, TolE, aB2d)`: the adaptor's
            // range is the `CurveOnSurface` range.
            let a_bac2d = Curve2dAdaptor::of_edge_face(&edge, &self.my_corrected_face);
            if let Some(c2d) = a_bac2d.curve() {
                crate::bnd_lib_add2d::add_adaptor_range(
                    c2d.as_ref(),
                    a_bac2d.first_parameter(),
                    a_bac2d.last_parameter(),
                    a_tol_e,
                    &mut a_b2d,
                );
            }
        }
        a_b2d
    }

    /// `MoveWire2d(aWire, aTrV)` (`cxx:911-977`): translate the p-curves of one
    /// wire on `myCorrectedFace` (seam edges get both of their p-curves moved).
    pub(super) fn move_wire2d(&mut self, a_wire: &TopoShape, a_tr_v: &GpVec2d) {
        if a_tr_v.magnitude() < occt_core::precision::CONFUSION {
            return;
        }

        let face = self.my_corrected_face.clone();
        let tol = self.my_face_tolerance;
        let a_seq_edges: Vec<TopoShape> = explorer_shapes(a_wire, ShapeType::Edge);
        let mut an_ed_map = IndexedShapeSet::new();

        for an_edge_shape in a_seq_edges {
            let Some(an_edge) = Edge::wrap(an_edge_shape.clone()) else {
                continue;
            };

            if BRepTool::is_closed_edge_face(&an_edge, &face) {
                // a. Closed edge (`cxx:938-969`).
                if an_ed_map.contains(&an_edge_shape) {
                    continue;
                }
                an_ed_map.add(&an_edge_shape);

                let an_ori = an_edge_shape.orientation();
                let (an_ef, an_er) = if an_ori == Orientation::Forward {
                    let mut rev = an_edge_shape.clone();
                    rev.set_orientation(Orientation::Reversed);
                    (an_edge_shape.clone(), rev)
                } else {
                    let mut rev = an_edge_shape.clone();
                    rev.set_orientation(Orientation::Forward);
                    (rev, an_edge_shape.clone())
                };

                let Some(a_tr_cf) = self.translate_curve2d(&an_ef, &face, a_tr_v) else {
                    continue;
                };
                let Some(a_tr_cr) = self.translate_curve2d(&an_er, &face, a_tr_v) else {
                    continue;
                };
                self.update_edge_seam(&an_edge, a_tr_cf, a_tr_cr, &face, tol);
            } else {
                // b. Usual edge (`cxx:970-975`).
                let Some(a_tr_c) = self.translate_curve2d(&an_edge_shape, &face, a_tr_v) else {
                    continue;
                };
                self.update_edge(&an_edge, a_tr_c, &face, tol);
            }
        }
    }

    /// `MoveWires2d(anOuterWire)` (`cxx:982-1030`).
    pub(super) fn move_wires_2d(&mut self, an_outer_wire: &TopoShape) {
        // 1. Move the outer wire to `[0, 2*PI]` (`cxx:989-1000`).
        let b2d_outer = self.bnd_box_wire(an_outer_wire);
        let Some((outer_u1, _outer_v1, outer_u2, _outer_v2)) = b2d_outer.get() else {
            // `Bnd_Box2d::Get` on a void box throws in OCCT.
            return;
        };
        let mut a = 0.5 * (outer_u1 + outer_u2);
        let mut b = (-a / TWO_PI).trunc() as i32;
        let mut k = if a < 0.0 { 1 } else { 0 };
        let mut x = (b + k) as f64 * TWO_PI;
        self.move_wire2d(an_outer_wire, &GpVec2d::new(x, 0.0));

        let b2d_outer = self.bnd_box_wire(an_outer_wire);
        let Some((outer_u1, _v1, _u2, _v2)) = b2d_outer.get() else {
            return;
        };

        // 2. Move all other wires into the bounding box of the outer wire
        //    (`cxx:1004-1028`).
        let wires = explorer_shapes(&self.my_corrected_face.0, ShapeType::Wire);
        for wire in wires {
            self.my_current_wire = wire.clone();
            if is_equal(&self.my_current_wire, an_outer_wire) {
                continue;
            }
            let b2d = self.bnd_box_wire(&wire);
            let Some((u1, _v1, _u2, _v2)) = b2d.get() else {
                continue;
            };

            if b2d.is_out_box(&b2d_outer) {
                a = u1 - outer_u1;
                b = (-a / TWO_PI).trunc() as i32;
                k = if a < 0.0 { 1 } else { 0 };
                x = (b + k) as f64 * TWO_PI;
                self.move_wire2d(&wire, &GpVec2d::new(x, 0.0));
            }
        }
    }

    /// `UpdateEdge(ECopy, C, F, Tol)` (`cxx:1033-1073`).
    pub(super) fn update_edge(&mut self, e_copy: &Edge, c: Arc<dyn Curve2d>, f: &Face, tol: f64) {
        // `E` is the original edge from the original face.
        let original = self.my_ed_map_inversed.find_from_key(&e_copy.0).cloned();
        if let Some(e) = original {
            self.write_on_original(&e, c.clone(), f, tol);
        }
        write_edge_pcurve(e_copy, c, f, tol);
    }

    /// `UpdateEdge(ECopy, C1, C2, F, Tol)` (`cxx:1076-1117`).
    pub(super) fn update_edge_seam(
        &mut self,
        e_copy: &Edge,
        c1: Arc<dyn Curve2d>,
        c2: Arc<dyn Curve2d>,
        f: &Face,
        tol: f64,
    ) {
        let original = self.my_ed_map_inversed.find_from_key(&e_copy.0).cloned();
        if let Some(e) = original {
            self.write_on_original_seam(&e, c1.clone(), c2.clone(), f, tol);
        }
        write_edge_pcurves(e_copy, c1, c2, f, tol);
    }

    /// The `pMap` arm shared by both `UpdateEdge` overloads
    /// (`cxx:1041-1051` / `:1085-1096`): the original edge is created once with
    /// `EmptyCopied` + its vertices, then updated in place.
    fn write_on_original(&mut self, e: &TopoShape, c: Arc<dyn Curve2d>, f: &Face, tol: f64) {
        let has = self
            .my_map
            .as_ref()
            .map(|m| m.borrow().contains(e))
            .unwrap_or(false);

        if !has {
            let mut an_edge = e.empty_copied();
            let builder = TopoBuilder::new();
            for v in explorer_shapes(e, ShapeType::Vertex) {
                builder.add(&mut an_edge, &v);
            }
            if let Some(oe) = Edge::wrap(an_edge.clone()) {
                write_edge_pcurve(&oe, c, f, tol);
            }
            if let Some(map) = &self.my_map {
                map.borrow_mut().add(e, an_edge);
            }
        } else if let Some(map) = &self.my_map {
            let mut m = map.borrow_mut();
            if let Some(v) = m.change_from_key(e) {
                if let Some(ve) = Edge::wrap(v.clone()) {
                    write_edge_pcurve(&ve, c, f, tol);
                }
            }
        }
    }

    /// The `pMap` arm of the seam overload.
    fn write_on_original_seam(
        &mut self,
        e: &TopoShape,
        c1: Arc<dyn Curve2d>,
        c2: Arc<dyn Curve2d>,
        f: &Face,
        tol: f64,
    ) {
        let has = self
            .my_map
            .as_ref()
            .map(|m| m.borrow().contains(e))
            .unwrap_or(false);

        if !has {
            let mut an_edge = e.empty_copied();
            let builder = TopoBuilder::new();
            for v in explorer_shapes(e, ShapeType::Vertex) {
                builder.add(&mut an_edge, &v);
            }
            if let Some(oe) = Edge::wrap(an_edge.clone()) {
                write_edge_pcurves(&oe, c1, c2, f, tol);
            }
            if let Some(map) = &self.my_map {
                map.borrow_mut().add(e, an_edge);
            }
        } else if let Some(map) = &self.my_map {
            let mut m = map.borrow_mut();
            if let Some(v) = m.change_from_key(e) {
                if let Some(ve) = Edge::wrap(v.clone()) {
                    write_edge_pcurves(&ve, c1, c2, f, tol);
                }
            }
        }
    }
}
