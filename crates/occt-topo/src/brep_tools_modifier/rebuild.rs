//! `BRepTools_Modifier::Rebuild` (`BRepTools_Modifier.cxx:226-649`): the
//! recursive walk that turns the maps filled by `Perform` into the rebuilt
//! shape tree.
//!
//! The port keeps the reference's two "new geometry" notions apart, because
//! the pcurve pass needs both:
//! - `newgeom` (`cxx:340`) — this shape got new geometry in the `switch`.
//! - `the_new_geom` (`cxx:341`, the out parameter) — this shape or any of its
//!   children got new geometry.
//!
//! Geometry frames follow the reference's `myLoc.Predivided(S.Location())`
//! (`cxx:261`, `cxx:306`). This port keeps the geometry of a shape in the
//! shape's own frame and its location on the shape view, so the effective
//! transform the modification asks for is baked into the geometry and the
//! source location stays on the rebuilt shape view.

use std::sync::Arc;

use occt_core::gp::{GpDir2d, GpPnt2d};
use occt_core::toploc::TopLocLocation;
use occt_geom::{Curve, Surface};

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape};
use crate::shape_fix_compose_shell::{empty_copied_edge, empty_copied_face};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{explorer_shapes, is_same_occurrence, shape_map_key};

use super::modification::Modification;
use super::modifier::{children_keep_orientation, Modifier};

impl Modifier {
    /// `Rebuild(S, M, theNewGeom)` (`cxx:226-649`).
    pub(super) fn rebuild(
        &mut self,
        s: &TopoShape,
        m: &dyn Modification,
        the_new_geom: &mut bool,
    ) -> Result<bool, String> {
        let ts = s.shape_type();
        let key = shape_map_key(s);
        // `cxx:236-241`: an already rebuilt shape answers from the map.
        if let Some(Some(result)) = self.my_map.get(&key) {
            *the_new_geom = self.my_has_new_geom.contains(&key);
            return Ok(!is_same_occurrence(s, result));
        }

        let b = TopoBuilder::new();
        let reg = GeometryRegistry::global();
        let mut rebuild = false;
        let mut rev_wires = false;
        let mut res_or = Orientation::Forward;
        let mut no_3d_curve = false;
        let mut result: Option<TopoShape> = None;

        match ts {
            ShapeType::Face => {
                // `cxx:253-270`.
                let ns = self.my_ns_info.get(&key).cloned();
                rebuild = ns.is_some();
                if let Some(info) = &ns {
                    rev_wires = info.rev_wires;
                    let face = self.make_new_face(&b, s, info);
                    if info.rev_face {
                        res_or = Orientation::Reversed;
                    }
                    // `B.NaturalRestriction(result, BRep_Tool::NaturalRestriction(S))`
                    // (`cxx:269`).
                    reg.set_natural_restriction(&face.0, reg.natural_restriction(s));
                    result = Some(face.0);
                }
                // `cxx:272-286`: the triangulation pass.
                if let Some(tri) = m.new_triangulation(s) {
                    match result.as_ref().and_then(|r| Face::wrap(r.clone())) {
                        // The copied face already exists => update it.
                        Some(face) => b.update_face_triangulation(&face, Some(tri), true),
                        // Create a new face with a bare triangulation.
                        None => {
                            // `B.MakeFace(result, aTriangulation); result.Location(S.Location(), false);`
                            // (`cxx:282-283`).
                            let mut face = b.make_face_from_triangulation(tri);
                            face.0.set_location(s.location());
                            result = Some(face.0);
                        }
                    }
                    rebuild = true;
                }
            }
            ShapeType::Edge => {
                // `cxx:290-316`.
                let nc = self.my_nc_info.get(&key).cloned();
                rebuild = nc.is_some();
                if let Some(info) = &nc {
                    let src_edge = Edge(s.clone());
                    match &info.curve {
                        // `cxx:295-301`: `B.MakeEdge(E)` has no curve
                        // representation at all; `Degenerated` copies the flag
                        // and `UpdateEdge(E, Tol)` raises the tolerance
                        // (OCC217).
                        None => {
                            let edge = b.make_shape(ShapeType::Edge);
                            reg.set_degenerated(&edge, BRepTool::is_degenerated(&src_edge));
                            reg.raise_edge_tolerance(&edge, info.tolerance);
                            result = Some(edge);
                            no_3d_curve = true;
                        }
                        // `cxx:302-309`: `B.MakeEdge(E, C, L', Tol)`.
                        Some(curve) => {
                            let loc = info.location.predivided(s.location());
                            let curve = bake_curve(curve, &loc);
                            let (first, last) = (curve.first_parameter(), curve.last_parameter());
                            let edge = b.make_edge(curve, first, last);
                            reg.raise_edge_tolerance(&edge.0, info.tolerance);
                            result = Some(edge.0);
                            no_3d_curve = false;
                        }
                    }
                    let out = result.as_mut().expect("edge arm always builds a shape");
                    // `result.Location(S.Location(), false)` (`cxx:310`). The
                    // geometry was baked with `myLoc.Predivided(S.Location())`,
                    // so the shape keeps the source location.
                    out.set_location(s.location());
                    // `B.SameParameter` / `B.SameRange` (`cxx:314-315`).
                    reg.set_same_parameter(out, reg.same_parameter(s));
                    reg.set_same_range(out, reg.same_range(s));
                }
            }
            _ => {}
        }

        // `cxx:338-380`: rebuild the children, tracking both notions of new
        // geometry, then make the empty copy when only the children changed.
        let newgeom = rebuild;
        *the_new_geom = rebuild;
        let children = children_keep_orientation(s);
        for child in &children {
            let mut sub_new_geom = false;
            let subrebuilt = self.rebuild(child, m, &mut sub_new_geom)?;
            rebuild = subrebuilt || rebuild;
            *the_new_geom = *the_new_geom || sub_new_geom;
        }
        if *the_new_geom {
            self.my_has_new_geom.insert(key.clone());
        }
        if rebuild && !newgeom {
            // `cxx:376-380`: `S.EmptyCopied()` with the orientation forced
            // forward. `TopoDS_Face::EmptyCopy` / `TopoDS_Edge::EmptyCopy` carry
            // the geometry over (`BRep_TFace.cxx:36-43`, `BRep_TEdge.cxx:105-129`);
            // a VERTEX can never reach this branch, because `rebuild` for a
            // vertex is `newgeom` (a vertex has no children).
            let mut copy = match ts {
                ShapeType::Face => empty_copied_face(s),
                ShapeType::Edge => empty_copied_edge(s),
                _ => s.empty_copied(),
            };
            copy.set_orientation(Orientation::Forward);
            result = Some(copy);
        }

        // `cxx:382-637`: copy the sub-elements.
        if rebuild {
            let mut out = result
                .clone()
                .ok_or_else(|| "BRepTools_Modifier::Rebuild: rebuilt shape without result".to_string())?;
            for child in &children {
                let mut orient = child.orientation();
                let mapped = self.mapped(child).ok_or_else(|| {
                    "BRepTools_Modifier::Rebuild: unbound sub-shape".to_string()
                })?;
                if rev_wires || mapped.orientation() == Orientation::Reversed {
                    orient = orient.reversed();
                }
                b.add(&mut out, &mapped.oriented(orient));
            }

            if ts == ShapeType::Face {
                self.rewrite_pcurves(s, &out, rev_wires, newgeom, *the_new_geom, m)?;
            } else if ts == ShapeType::Edge && !no_3d_curve {
                // `cxx:590-630`.
                let edge = s;
                let mut edor = edge.orientation();
                if edor != Orientation::Reversed {
                    edor = Orientation::Forward;
                }
                for vertex in explorer_shapes(&edge.oriented(edor), ShapeType::Vertex) {
                    let (param, tol) = self.parameter_of(&vertex, edge, m);
                    let mut vtxrelat = vertex.orientation();
                    if edor == Orientation::Reversed {
                        vtxrelat = vtxrelat.reversed();
                    }
                    let local_v = self
                        .mapped(&vertex)
                        .ok_or_else(|| "BRepTools_Modifier::Rebuild: unbound vertex".to_string())?
                        .oriented(vtxrelat);
                    if self.my_mutable_input || !is_same_occurrence(&local_v, &vertex) {
                        update_vertex(&local_v, param, &out, tol);
                    }
                }
            }

            // `cxx:632-636`: update flags. `Orientable(S.Orientable())` has no
            // counterpart in `ShapeFlags`.
            out.set_closed(s.is_closed());
            out.set_infinite(s.infinite());
            result = Some(out);
        } else {
            // `cxx:638-641`.
            result = Some(s.clone());
        }

        // `cxx:643-646`.
        let mut out = result.expect("a result is always assigned");
        out.set_orientation(res_or);
        set_shape_flags(s, &out);
        self.my_map.insert(key, Some(out));
        Ok(rebuild)
    }

    /// `B.MakeFace(TopoDS::Face(result), mySurface, myLoc.Predivided(S.Location()), myToler)`
    /// (`cxx:259-263`) followed by `result.Location(S.Location(), false)`.
    fn make_new_face(&self, b: &TopoBuilder, s: &TopoShape, info: &super::modification::NewSurfaceInfo) -> Face {
        let loc = info.location.predivided(s.location());
        let surface = bake_surface(&info.surface, &loc);
        let mut face = b.make_face_with_tolerance(surface, &[], info.tolerance);
        face.0.set_location(s.location());
        face
    }

    /// `NewParameter(V, E, P, Tol)` with the reference's fallback to the
    /// stored parameter and the vertex tolerance
    /// (`BRepTools_Modifier.cxx:522-526`, `cxx:605-609`).
    fn parameter_of(&self, vertex: &TopoShape, edge: &TopoShape, m: &dyn Modification) -> (f64, f64) {
        match m.new_parameter(vertex, edge) {
            Some(info) => (info.parameter, info.tolerance),
            None => (
                BRepTool::parameter_of_vertex(edge, vertex),
                GeometryRegistry::global().vertex_tolerance(vertex),
            ),
        }
    }

    /// `cxx:397-587`: the pcurve pass of a rebuilt face.
    fn rewrite_pcurves(
        &mut self,
        s: &TopoShape,
        out: &TopoShape,
        rev_wires: bool,
        newgeom: bool,
        the_new_geom: bool,
        m: &dyn Modification,
    ) -> Result<(), String> {
        let reg = GeometryRegistry::global();
        let face = s;
        // `cxx:402-406`: `fcor` is the *source* face orientation with
        // INTERNAL/EXTERNAL folded to FORWARD.
        let mut fcor = face.orientation();
        if fcor != Orientation::Reversed {
            fcor = Orientation::Forward;
        }
        let out_face_key = GeometryRegistry::shape_key(out);
        let Some(out_surface) = reg.face_surface(out) else {
            return Ok(());
        };
        for edge in explorer_shapes(&face.oriented(fcor), ShapeType::Edge) {
            let mapped_edge = self
                .mapped(&edge)
                .ok_or_else(|| "BRepTools_Modifier::Rebuild: unbound edge".to_string())?;
            // `cxx:416-422`.
            if the_new_geom {
                if let Some(curve2d) = m.new_curve2d(&edge, face, &mapped_edge, out) {
                    // `cxx:433-473`: a true seam edge of this face, or an edge
                    // shared by two faces built on the same surface (OCC21772).
                    let mut is_closed = false;
                    let src_edge = Edge(edge.clone());
                    let src_face = Face(face.clone());
                    if BRepTool::is_closed_edge_face(&src_edge, &src_face) {
                        is_closed = !newgeom
                            || crate::brep_tools::is_really_closed(&src_edge, &src_face);
                        if !is_closed {
                            is_closed = self.another_face_shares(&edge, face);
                        }
                    }
                    if is_closed {
                        // `cxx:476-505`.
                        let cur_e = mapped_edge.clone();
                        let mut cur_f = out.clone();
                        cur_f.set_orientation(Orientation::Forward);
                        let forward = (!rev_wires && fcor != edge.orientation())
                            || (rev_wires && fcor == edge.orientation());
                        // `CurE.Orientation(TopAbs_FORWARD)` / `(REVERSED)`
                        // (`cxx:485`, `cxx:497`) decides which pcurve of the
                        // closed-surface representation
                        // `BRep_Tool::CurveOnSurface(CurE, CurF, f, l)`
                        // returns: `PCurve1` for a FORWARD edge, `PCurve2` for
                        // a REVERSED one (`BRep_Tool.cxx:334`, `:354-357`).
                        // `CurF` is oriented FORWARD (`cxx:477-478`), so
                        // `BRep_Tool.cxx:310-313` adds no reversal.
                        let slot = if forward { 0 } else { 1 };
                        let first_pc = reg
                            .edge_pcurves(&cur_e, out_face_key)
                            .get(slot)
                            .cloned()
                            .unwrap_or_else(origin_line2d);
                        let pair = if forward {
                            vec![first_pc, curve2d]
                        } else {
                            vec![curve2d, first_pc]
                        };
                        reg.set_edge_pcurves_on_surface(&cur_e, &out_surface, pair, 0.0);
                        // `currcurv = BRep_Tool::CurveOnSurface(CurE, CurF, f, l); B.Range(CurE, f, l);`
                        if let Some((f, l)) = reg.pcurve_range(&cur_e, out_face_key) {
                            reg.set_ranges_all(&cur_e, f, l);
                        }
                    } else {
                        // `cxx:508`: `B.UpdateEdge(myMap(ex.Current()), curve2d, result, 0.)`.
                        reg.update_edge_pcurve_on_surface(&mapped_edge, &out_surface, Some(curve2d), 0.0);
                    }

                    // `cxx:511-545`: a pcurve-only edge has no 3D curve, so
                    // the vertex parameters have to be written here.
                    if reg.edge_curve(&mapped_edge).is_none() {
                        for vertex in explorer_shapes(&edge, ShapeType::Vertex) {
                            let (param, tol) = self.parameter_of(&vertex, &edge, m);
                            let mut vtxrelat = vertex.orientation();
                            if edge.orientation() == Orientation::Reversed {
                                vtxrelat = vtxrelat.reversed();
                            }
                            let local_v = self
                                .mapped(&vertex)
                                .ok_or_else(|| "BRepTools_Modifier::Rebuild: unbound vertex".to_string())?
                                .oriented(vtxrelat);
                            update_vertex(&local_v, param, &mapped_edge, tol);
                        }
                    }
                }
            }
            // `cxx:548-585`: "Copy polygon on triangulation" — UNPORTED. The
            // port keeps no `Poly_PolygonOnTriangulation`
            // (`BRepTools_Modification.hxx:91-94`, see `modification.rs`), so
            // `NewPolygonOnTriangulation` and the two `B.UpdateEdge` seam
            // overloads have no counterpart.
        }
        Ok(())
    }

    /// `cxx:439-473`: another face of the modified shape that is built on the
    /// same surface and holds the very same edge makes this edge a true seam.
    ///
    /// The reference also compares the surfaces' locations (`aLoc.IsEqual(anOtherLoc)`);
    /// the port has no per-face surface location, so only the surface handle
    /// is compared (`aSurf == anOtherSurf`).
    fn another_face_shares(&self, edge: &TopoShape, face: &TopoShape) -> bool {
        let reg = GeometryRegistry::global();
        let resface = self.mapped_or_self(face);
        let Some(a_surf) = resface.as_ref().and_then(|f| reg.face_surface(f)) else {
            return false;
        };
        let Some(shape) = self.my_shape.clone() else { return false };
        for other in explorer_shapes(&shape, ShapeType::Face) {
            if is_same_occurrence(&other, face) {
                continue;
            }
            let Some(resface2) = self.mapped_or_self(&other) else { continue };
            let Some(other_surf) = reg.face_surface(&resface2) else { continue };
            if !Arc::ptr_eq(&a_surf, &other_surf) {
                continue;
            }
            for e in explorer_shapes(&other, ShapeType::Edge) {
                if is_same_occurrence(edge, &e) {
                    return true;
                }
            }
        }
        false
    }
}

/// The geometry frame `myLoc.Predivided(S.Location())` applied to a surface.
fn bake_surface(surface: &Arc<dyn Surface>, loc: &TopLocLocation) -> Arc<dyn Surface> {
    if loc.is_identity() {
        return surface.clone();
    }
    Arc::from(surface.transformed(&loc.transformation()))
}

/// The geometry frame `myLoc.Predivided(S.Location())` applied to a curve.
fn bake_curve(curve: &Arc<dyn Curve>, loc: &TopLocLocation) -> Arc<dyn Curve> {
    if loc.is_identity() {
        return curve.clone();
    }
    Arc::from(curve.transformed(&loc.transformation()))
}

/// `SetShapeFlags` (`cxx:834-842`).
///
/// UNPORTED: `Checked`, `Orientable` and `Convex` (`cxx:837`, `:838`, `:841`) —
/// `ShapeFlags` (`abs.rs`) carries only `Free`, `Closed`, `Infinite` and
/// `Modified`.
fn set_shape_flags(src: &TopoShape, dst: &TopoShape) {
    let (modified, closed, infinite) = {
        let ts = src.tshape.read().expect("poisoned TShape lock");
        (ts.modified(), ts.is_closed(), ts.infinite())
    };
    let mut ts = dst.tshape.write().expect("poisoned TShape lock");
    ts.set_modified(modified);
    ts.set_closed(closed);
    ts.set_infinite(infinite);
}

/// `new Geom2d_Line(gp::OX2d())` (`BRepTools_Modifier.cxx:489`, `cxx:499`).
///
/// UNPORTED: the port has no `gp::OX2d()` constant nor a `Geom2d_Line`
/// constructor from it in this crate; the pcurve is a placeholder OCCT only
/// installs for a seam edge that has no pcurve on the result face at all.
fn origin_line2d() -> Arc<dyn occt_geom2d::curve::Curve2d> {
    Arc::new(occt_geom2d::line::Geom2dLine::from_pnt_dir(
        GpPnt2d::new(0.0, 0.0),
        GpDir2d::new(1.0, 0.0).expect("unit X direction"),
    ))
}

/// `BRep_Builder::UpdateVertex(V, Par, E, Tol)` (`BRep_Builder.cxx:1220-1312`).
///
/// The edge's representation list is walked once (`cxx:1272-1304`) and every
/// `BRep_GCurve` gets `First(Par)` for a FORWARD vertex and `Last(Par)` for a
/// REVERSED one; the port's `set_ranges_all` is that whole-list write
/// (`tgeometry::registry_core`).
///
/// UNPORTED: the INTERNAL/EXTERNAL arm's `UpdatePoints` (`cxx:1287-1300`),
/// which appends a `BRep_PointRepresentation` to the `BRep_TVertex` — the port
/// keeps a single point plus one parameter per edge
/// (`set_vertex_parameter_on_edge`) and no representation list, the same
/// reduction the existing port of this builder makes
/// (`topopebrep::ds::build_tool` `update_vertex`).
fn update_vertex(v: &TopoShape, par: f64, e: &TopoShape, tol: f64) {
    let reg = GeometryRegistry::global();
    // `cxx:1241-1266`: the vertex's orientation inside the edge.
    // `TopoDS_Iterator itv(E.Oriented(TopAbs_FORWARD))` is the direct-children
    // walk, so the stored orientations are the ones compared.
    let mut ori = Orientation::Internal;
    let verts = children_keep_orientation(e);
    let degenerate = e.shape_type() == ShapeType::Edge && BRepTool::is_degenerated(&Edge(e.clone()));
    if verts.is_empty() && degenerate {
        ori = v.orientation();
    }
    for vcur in &verts {
        if is_same_occurrence(vcur, v) {
            ori = vcur.orientation();
            if ori == v.orientation() {
                break;
            }
        }
    }
    let (first, last) = reg.edge_parameters(e);
    match ori {
        Orientation::Forward => reg.set_ranges_all(e, par, last),
        Orientation::Reversed => reg.set_ranges_all(e, first, par),
        _ => {
            reg.set_vertex_parameter_on_edge(e, v, par);
            if let Some(c) = reg.edge_curve(e) {
                set_vertex_point(v, &c.d0(par));
            }
        }
    }
    if !matches!(ori, Orientation::Forward | Orientation::Reversed) {
        v.tshape.write().expect("poisoned TShape lock").set_modified(true);
    }
    // `TV->UpdateTolerance(Tol)` (`BRep_TVertex.cxx:29-35`): raise only.
    reg.raise_vertex_tolerance(v, tol);
    e.tshape.write().expect("poisoned TShape lock").set_modified(true);
}

/// The vertex point write of `UpdatePoints` (`BRep_Builder.cxx:1314-1390`):
/// the vertex keeps its tolerance.
fn set_vertex_point(v: &TopoShape, p: &occt_core::gp::GpPnt) {
    let reg = GeometryRegistry::global();
    let tol = reg.vertex_tolerance(v);
    reg.set_vertex(v, crate::tgeometry::VertexGeom { point: *p, tolerance: tol });
}
