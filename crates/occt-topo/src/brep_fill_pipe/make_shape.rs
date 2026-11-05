//! `BRepFill_Pipe::MakeShape`, `FindEdge`, `FindVertex`,
//! `RebuildTopOrBottomFace` and `UpdateTolFromTopOrBottomPCurve`
//! (`BRepFill_Pipe.cxx:96-146`, `:530-851`, `:860-990`, `:1058-1090`).

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::precision::CONFUSION;

use occt_geom::curve::Curve;
use occt_geom::curve_on_surface::CurveOnSurface;
use occt_geom::geom_fill::sweep::ApproxStyle;

use crate::abs::{Orientation, ShapeType};
use crate::brep_class3d::SolidClassifier;
use crate::brep_fill_pipe::Pipe;
use crate::brep_fill_section_law::SectionLaw;
use crate::brep_fill_shape_law::ShapeLaw;
use crate::brep_fill_sweep::helpers::ShapeTable;
use crate::brep_fill_sweep::{BrepFillSweep, TransitionStyle};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::iterator::ShapeIterator;
use crate::shape::{Edge, Face, Shell, Solid, TopoShape, Vertex, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, explorer_shapes, is_same, shape_map_key, ShapeMapKey};

/// The `UpdateMap` static of `BRepFill_Pipe.cxx:66-94`: append `the_value` to
/// the list bound to `the_key` unless it is already there, answering whether it
/// was appended.
pub(crate) fn update_map(
    the_map: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    the_key: &TopoShape,
    the_value: &TopoShape,
) -> bool {
    let list = the_map.entry(shape_map_key(the_key)).or_default();
    if list.iter().any(|value| is_same(value, the_value)) {
        return false;
    }
    list.push(the_value.clone());
    true
}

/// `NCollection_Array2` concatenation of `cxx:735-786`, one table at a time.
///
/// The reference copies `Old->Value(ii, jj)` for `ii` over `Old->ColLength()`
/// and `jj` over `Old->RowLength()`, i.e. it reads the old table with the two
/// bounds swapped relative to the destination (`Somme->SetValue(ii, jj)`).
/// Out-of-range reads are undefined in OCCT; [`ShapeTable`] answers `None` and
/// `set_value` ignores out-of-range writes.
pub(crate) fn concat_tables(old: &ShapeTable, aux: &ShapeTable) -> ShapeTable {
    let old_cols = old.nb_cols();
    let length = aux.nb_cols() + old_cols;
    let mut somme = ShapeTable::new(length, aux.nb_rows());
    for jj in 1..=old.nb_rows() {
        for ii in 1..=old_cols {
            let value = old.value(ii, jj).cloned();
            somme.set_value(ii, jj, value);
        }
        for kk in 1..=aux.nb_cols() {
            let ii = old_cols + kk;
            let value = aux.value(kk, jj).cloned();
            somme.set_value(ii, jj, value);
        }
    }
    somme
}

/// The `BRepFill_Sweep::Build` arguments of `cxx:684-687`, `:717-720` and
/// `:521`: `Modified` transition and `Location` approximation in every call
/// site, the class' own degree / segment / continuity limits, and the three
/// history maps the sweep fills.
pub(crate) fn run_sweep(pipe: &mut Pipe, mk_sw: &mut BrepFillSweep) -> Result<(), String> {
    let continuity = pipe.my_continuity;
    let degmax = pipe.my_degmax;
    let segmax = pipe.my_segmax;
    mk_sw.build(
        &mut pipe.my_reversed_edges,
        &mut pipe.my_tapes,
        &mut pipe.my_rails,
        TransitionStyle::Modified,
        continuity,
        ApproxStyle::Location,
        degmax,
        segmax,
    )
}

impl Pipe {
    ///
    /// `None` is the null shape the reference returns for the `default` arm of
    /// its type switch and for a sweep that produced nothing.
    pub(crate) fn make_shape(
        &mut self,
        s: &TopoShape,
        original_s: &TopoShape,
        first_shape: Option<&TopoShape>,
        last_shape: Option<&TopoShape>,
    ) -> Result<Option<TopoShape>, String> {
        let builder = TopoBuilder::new();
        let mut explode = false;
        let mut the_s = s.clone();
        let mut the_first = first_shape.cloned();
        let mut the_last = last_shape.cloned();
        // `cxx:540-545`.
        let initial_length = self.faces().map(|t| t.nb_cols()).unwrap_or(0);

        // `cxx:548`: `TopLoc_Location BackLoc(myTrsf.Inverted())` is built and
        // never read in 7.9.3.

        // `cxx:557-621`: the empty result of the right type.
        let mut result: Option<TopoShape> = match s.shape_type() {
            ShapeType::Vertex => Some(TopoShape::new(ShapeType::Wire)),
            ShapeType::Edge => {
                // `cxx:565-589`.
                let shell = TopoShape::new(ShapeType::Shell);
                let mut w = TopoShape::new(ShapeType::Wire);
                builder.add(&mut w, s);
                w.set_closed(BRepTool::is_closed(s));
                the_s = w;
                if let Some(first) = the_first.clone() {
                    let mut wf = TopoShape::new(ShapeType::Wire);
                    builder.add(&mut wf, &first);
                    wf.set_closed(BRepTool::is_closed(&first));
                    the_first = Some(wf);
                }
                if let Some(last) = the_last.clone() {
                    let mut wl = TopoShape::new(ShapeType::Wire);
                    builder.add(&mut wl, &last);
                    wl.set_closed(BRepTool::is_closed(&last));
                    the_last = Some(wl);
                }
                shell.set_closed(BRepTool::is_closed(&shell));
                Some(shell)
            }
            ShapeType::Wire => Some(TopoShape::new(ShapeType::Shell)),
            ShapeType::Face => {
                // `cxx:595-602`.
                let mut shell = TopoShape::new(ShapeType::Shell);
                explode = true;
                if !self.spine_closed() {
                    if let Some(first) = the_first.as_ref() {
                        let mut reversed = first.clone();
                        reversed.reverse();
                        builder.add(&mut shell, &reversed);
                    }
                }
                shell.set_closed(BRepTool::is_closed(&shell));
                Some(shell)
            }
            // `cxx:606-609`.
            ShapeType::Shell => {
                explode = true;
                Some(TopoShape::new(ShapeType::CompSolid))
            }
            // `cxx:611-614`.
            ShapeType::Solid | ShapeType::CompSolid => {
                return Err("Standard_DomainError: BRepFill_Pipe::profile contains solids".to_string())
            }
            // `cxx:617-620`.
            ShapeType::Compound => {
                explode = true;
                Some(TopoShape::new(ShapeType::Compound))
            }
            // `cxx:621-622`: `default: break;` leaves `result` null.
            _ => None,
        };

        if explode {
            // `cxx:629-668`.
            let first_children: Vec<TopoShape> = the_first
                .as_ref()
                .map(|shape| ShapeIterator::of_shape(shape).collect())
                .unwrap_or_default();
            let last_children: Vec<TopoShape> = the_last
                .as_ref()
                .map(|shape| ShapeIterator::of_shape(shape).collect())
                .unwrap_or_default();
            let children: Vec<TopoShape> = ShapeIterator::of_shape(&the_s).collect();
            let originals: Vec<TopoShape> = ShapeIterator::of_shape(original_s).collect();
            for (idx, child) in children.iter().enumerate() {
                // `itorig.Next()` runs in the `for` increment of `it`, so an
                // original tree with fewer children would read past its end;
                // `TopoDS_Iterator::Value()` on an exhausted iterator is
                // undefined in OCCT, so the port skips the pair.
                let Some(orig) = originals.get(idx) else {
                    continue;
                };
                let first = if the_first.is_some() {
                    first_children.get(idx)
                } else {
                    None
                };
                let last = if the_last.is_some() { last_children.get(idx) } else { None };
                if the_s.shape_type() == ShapeType::Face {
                    // `cxx:654-657`: the recursion's answer is dropped.
                    self.make_shape(child, orig, first, last)?;
                } else if let Some(sub) = self.make_shape(child, orig, first, last)? {
                    if let Some(res) = result.as_mut() {
                        builder.add(res, &sub);
                    }
                }
            }
        } else {
            match the_s.shape_type() {
                // `cxx:672-702`.
                ShapeType::Vertex => {
                    let Some(vertex) = Vertex::wrap(the_s.clone()) else {
                        return Err("Standard_TypeMismatch: BRepFill_Pipe::MakeShape".to_string());
                    };
                    let section: Arc<dyn SectionLaw> =
                        Arc::new(ShapeLaw::with_vertex(&vertex, true)?);
                    let location = self
                        .loc()
                        .ok_or_else(|| "BRepFill_Pipe::MakeShape: no location law".to_string())?;
                    let mut mk_sw = BrepFillSweep::new(section, location, true);
                    mk_sw.set_force_approx_c1(self.force_approx_c1());
                    let _ = run_sweep(self, &mut mk_sw);
                    let built = mk_sw.shape().cloned();
                    // `cxx:688`: a null result cannot be stored in the list.
                    if let Some(shape) = built.as_ref() {
                        update_map(self.gen_map_mut(), original_s, shape);
                    }
                    // `cxx:689`.
                    self.set_error_on_surf(mk_sw.error_on_surface());

                    // `cxx:692-699`.
                    let a_sections = mk_sw.sections();
                    let a_v_last = a_sections.nb_cols();
                    let a_first = a_sections.value(1, 1).cloned();
                    let a_last = a_sections.value(1, a_v_last).cloned();
                    self.set_first(a_first);
                    self.set_last(a_last);

                    // `cxx:701`.
                    self.build_history(&mk_sw, original_s)?;
                    return Ok(built);
                }

                // `cxx:704-690`: the wire section, including the labeling of
                // `myFaces` / `mySections` / `myEdges`.
                ShapeType::Wire => {
                    let Some(wire) = Wire::wrap(the_s.clone()) else {
                        return Err("Standard_TypeMismatch: BRepFill_Pipe::MakeShape".to_string());
                    };
                    let section: Arc<dyn SectionLaw> = Arc::new(ShapeLaw::with_wire(&wire, true)?);
                    let location = self
                        .loc()
                        .ok_or_else(|| "BRepFill_Pipe::MakeShape: no location law".to_string())?;
                    let mut mk_sw = BrepFillSweep::new(section, location, true);
                    mk_sw.set_bounds(
                        the_first.clone().and_then(Wire::wrap),
                        the_last.clone().and_then(Wire::wrap),
                    );
                    mk_sw.set_force_approx_c1(self.force_approx_c1());
                    let _ = run_sweep(self, &mut mk_sw);
                    result = mk_sw.shape().cloned();
                    self.set_error_on_surf(mk_sw.error_on_surface());

                    // `cxx:721-787`: "Labeling of elements".
                    let aux_faces = mk_sw.sub_shape().clone();
                    let aux_sections = mk_sw.sections().clone();
                    let aux_edges = mk_sw.inter_faces().clone();
                    // `cxx:722`: the guard is on `mySections`.
                    match self.sections().cloned() {
                        None => {
                            self.set_faces(Some(aux_faces));
                            self.set_sections(Some(aux_sections));
                            self.set_edges(Some(aux_edges));
                        }
                        Some(old_faces) => {
                            self.set_faces(Some(concat_tables(&old_faces, &aux_faces)));
                            let old_sections = self.sections().cloned().unwrap_or_default();
                            let old_sections_len = old_sections.nb_cols();
                            self.set_sections(Some(concat_tables(&old_sections, &aux_sections)));
                            // `cxx:761`: the first edge of the appended block.
                            self.set_cur_index_of_section_edge(old_sections_len + 1);
                            let old_edges = self.edges().cloned().unwrap_or_default();
                            self.set_edges(Some(concat_tables(&old_edges, &aux_edges)));
                        }
                    }

                    // `cxx:789`.
                    self.build_history(&mk_sw, original_s)?;
                }
                _ => {}
            }
        }

        // `cxx:793-848`: the face section is wrapped into a solid.
        if the_s.shape_type() == ShapeType::Face {
            let Some(result_shape) = result.clone() else {
                return Ok(None);
            };
            let Some(shell) = Shell::wrap(result_shape.clone()) else {
                return Ok(None);
            };
            // `cxx:796-804`: "jgv: top face was reversed".
            let faces_to_rebuild: Vec<TopoShape> =
                explorer_shapes(&result_shape, ShapeType::Face);
            for a_face in faces_to_rebuild {
                let mut reversed = a_face;
                reversed.reverse();
                self.rebuild_top_or_bottom_face(&reversed, true)?;
            }

            // `cxx:806-818`.
            let nb_cols = self.faces().map(|t| t.nb_cols()).unwrap_or(0);
            let nb_rows = self.faces().map(|t| t.nb_rows()).unwrap_or(0);
            let new_faces: Vec<TopoShape> = {
                let mut collected = Vec::new();
                if let Some(table) = self.faces() {
                    for ii in (initial_length + 1)..=nb_cols {
                        for jj in 1..=nb_rows {
                            let Some(value) = table.value(ii, jj) else {
                                continue;
                            };
                            if value.shape_type() == ShapeType::Face {
                                collected.push(value.clone());
                            }
                        }
                    }
                }
                collected
            };
            if let Some(res) = result.as_mut() {
                for a_face in &new_faces {
                    builder.add(res, a_face);
                }
            }

            // `cxx:820-829`: "if Spine is not closed add the last face".
            if !self.spine_closed() {
                let Some(the_last_shape) = the_last.clone() else {
                    return Ok(None);
                };
                // `cxx:826`: "jgv: bottom face".
                self.rebuild_top_or_bottom_face(&the_last_shape, false)?;
                let Some(last_face) = Face::wrap(the_last_shape) else {
                    return Err("Standard_TypeMismatch: BRepFill_Pipe::MakeShape".to_string());
                };
                if let Some(res) = result.as_mut() {
                    builder.add(res, &last_face.0);
                }
            }

            // `cxx:831-848`.
            let mut solid = Solid::new();
            shell.0.set_closed(true);
            builder.add_shell(&mut solid, &shell);
            let mut classifier = SolidClassifier::new();
            classifier.load(solid.0.clone());
            classifier.perform_infinite_point(CONFUSION);
            if classifier.state() == FaceState::In {
                let mut reversed_solid = Solid::new();
                let mut reversed = shell.0.clone();
                reversed.reverse();
                let Some(reversed_shell) = Shell::wrap(reversed) else {
                    return Err("Standard_TypeMismatch: BRepFill_Pipe::MakeShape".to_string());
                };
                builder.add_shell(&mut reversed_solid, &reversed_shell);
                solid = reversed_solid;
            }
            update_map(self.gen_map_mut(), original_s, &solid.0);
            return Ok(Some(solid.0));
        }

        Ok(result)
    }

    /// `FindEdge(S, E, InitialLength)` (`cxx:860-910`): the 1-based rank of
    /// `e` among the edges `S` contributes to the section law, or `0`.
    pub(crate) fn find_edge(
        &self,
        s: &TopoShape,
        e: &Edge,
        initial_length: &mut i32,
    ) -> Result<i32, String> {
        let mut result = 0i32;
        match s.shape_type() {
            // `cxx:868-876`.
            ShapeType::Edge => {
                *initial_length += 1;
                if is_same(s, &e.0) {
                    result = *initial_length;
                }
            }
            // `cxx:878-890`.
            ShapeType::Wire => {
                let wire = Wire::wrap(s.clone())
                    .ok_or_else(|| "Standard_TypeMismatch: BRepFill_Pipe::FindEdge".to_string())?;
                let section = ShapeLaw::with_wire(&wire, false)?;
                let nb_law = section.base().nb_law() as i32;
                let mut ii = 1;
                while ii <= nb_law && result == 0 {
                    if let Some(edge) = section.edge(ii as usize) {
                        if is_same(&edge.0, &e.0) {
                            result = *initial_length + ii;
                        }
                    }
                    ii += 1;
                }
                *initial_length += nb_law;
            }
            // `cxx:892-900`.
            ShapeType::Face | ShapeType::Shell | ShapeType::Compound => {
                for child in ShapeIterator::of_shape(s) {
                    if result != 0 {
                        break;
                    }
                    result = self.find_edge(&child, e, initial_length)?;
                }
            }
            // `cxx:902-905`.
            ShapeType::Solid | ShapeType::CompSolid => {
                return Err("Standard_DomainError: BRepFill_Pipe::SOLID or COMPSOLID".to_string())
            }
            _ => {}
        }
        Ok(result)
    }

    /// `FindVertex(S, V, InitialLength)` (`cxx:918-990`): the 1-based rank of
    /// `v` among the section vertices `S` contributes, or `0`.
    pub(crate) fn find_vertex(
        &self,
        s: &TopoShape,
        v: &Vertex,
        initial_length: &mut i32,
    ) -> Result<i32, String> {
        let mut result = 0i32;
        match s.shape_type() {
            // `cxx:927-934`.
            ShapeType::Vertex => {
                *initial_length += 1;
                if is_same(s, &v.0) {
                    result = *initial_length;
                }
            }
            // `cxx:936-956`.
            ShapeType::Edge => {
                let edge = Edge::wrap(s.clone()).ok_or_else(|| {
                    "Standard_TypeMismatch: BRepFill_Pipe::FindVertex".to_string()
                })?;
                let (mut first, mut last) = edge_vertices(&edge);
                if s.orientation() == Orientation::Reversed {
                    std::mem::swap(&mut first, &mut last);
                }
                let mut matched = false;
                if let Some(vertex) = first.as_ref() {
                    if is_same(&vertex.0, &v.0) {
                        result = *initial_length + 1;
                        matched = true;
                    }
                }
                if !matched {
                    if let Some(vertex) = last.as_ref() {
                        if is_same(&vertex.0, &v.0) {
                            result = *initial_length + 2;
                        }
                    }
                }
                *initial_length += 2;
            }
            // `cxx:958-970`.
            ShapeType::Wire => {
                let mut ii = *initial_length + 1;
                let wire = Wire::wrap(s.clone()).ok_or_else(|| {
                    "Standard_TypeMismatch: BRepFill_Pipe::FindVertex".to_string()
                })?;
                let section = ShapeLaw::with_wire(&wire, false)?;
                *initial_length += section.base().nb_law() as i32 + 1;
                while ii <= *initial_length && result == 0 {
                    if let Some(vertex) = section.vertex_at(ii as usize, 0.0) {
                        if is_same(&vertex.0, &v.0) {
                            result = ii;
                        }
                    }
                    ii += 1;
                }
            }
            // `cxx:972-980`.
            ShapeType::Face | ShapeType::Shell | ShapeType::Compound => {
                for child in ShapeIterator::of_shape(s) {
                    if result != 0 {
                        break;
                    }
                    result = self.find_vertex(&child, v, initial_length)?;
                }
            }
            // `cxx:982-985`.
            ShapeType::Solid | ShapeType::CompSolid => {
                return Err("Standard_DomainError: BRepFill_Pipe::SOLID or COMPSOLID".to_string())
            }
            _ => {}
        }
        Ok(result)
    }

    /// `RebuildTopOrBottomFace(aFace, IsTop)` (`cxx:1058-1090`): swap the
    /// wire's v-iso edges for the corresponding section edges.
    fn rebuild_top_or_bottom_face(
        &mut self,
        a_face: &TopoShape,
        is_top: bool,
    ) -> Result<(), String> {
        // `cxx:1060`: `IsTop ? 1 : mySections->RowLength()`. The second value is
        // used as a *column* index below; see the module note on the
        // reference's swapped bounds.
        let Some(sections) = self.sections().cloned() else {
            return Err("BRepFill_Pipe::RebuildTopOrBottomFace: no sections".to_string());
        };
        let index_of_section = if is_top { 1 } else { sections.nb_rows() };
        let builder = TopoBuilder::new();

        // `cxx:1063-1091`.
        for a_wire in ShapeIterator::of_shape(a_face) {
            let mut init_edges: Vec<TopoShape> = Vec::new();
            let mut res_edges: Vec<TopoShape> = Vec::new();
            for an_edge in ShapeIterator::of_shape(&a_wire) {
                for ii in self.cur_index_of_section_edge()..=sections.nb_cols() {
                    let Some(a_viso_edge) = sections.value(ii, index_of_section) else {
                        continue;
                    };
                    if is_same(&an_edge, a_viso_edge) {
                        init_edges.push(an_edge.clone());
                        res_edges.push(a_viso_edge.clone());
                        break;
                    }
                }
            }
            let mut wire = a_wire;
            wire.set_free(true);
            for ii in 0..init_edges.len() {
                builder.remove(&mut wire, &init_edges[ii]);
                update_tol_from_top_or_bottom_pcurve(a_face, &res_edges[ii])?;
                builder.add(&mut wire, &res_edges[ii]);
            }
        }
        Ok(())
    }
}

/// `UpdateTolFromTopOrBottomPCurve(aFace, anEdge)` (`cxx:96-146`): raise the
/// edge (and its vertices) to the largest gap between the 3D curve and the
/// pcurve lifted onto the face, sampled at `NCONTROL + 1` parameters.
///
/// `BRep_Tool::Curve(anEdge, fpar, lpar)` overwrites the range the pcurve
/// lookup returned (`cxx:108`), so both curves are sampled over the *3D*
/// range — which is what the reference passes to `Geom2dAdaptor_Curve`
/// (`cxx:112`).
fn update_tol_from_top_or_bottom_pcurve(
    a_face: &TopoShape,
    an_edge: &TopoShape,
) -> Result<(), String> {
    let geometry = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(a_face);
    // `cxx:98-104`: a null pcurve returns without touching the tolerance.
    let Some(p_curve) = geometry.edge_pcurve(an_edge, face_key) else {
        return Ok(());
    };
    // `cxx:106-109`.
    let Some(a_curve) = geometry.edge_curve(an_edge) else {
        return Ok(());
    };
    let (fpar, lpar) = BRepTool::edge_parameters(&crate::shape::Edge(an_edge.clone()));
    let Some(a_surf) = BRepTool::face_surface(&Face(a_face.clone())) else {
        return Err("BRepFill_Pipe::UpdateTolFromTopOrBottomPCurve: no surface".to_string());
    };

    // `cxx:111-117`.
    let con_s = CurveOnSurface::new(p_curve, a_surf);
    let mut tol = BRepTool::edge_tolerance(&crate::shape::Edge(an_edge.clone()));
    let init_tol = tol;
    let mut tol_tol = tol * tol;
    const NCONTROL: usize = 22;
    let delta = (lpar - fpar) / NCONTROL as f64;
    for i in 0..=NCONTROL {
        let par = fpar + i as f64 * delta;
        let pnt = a_curve.d0(par);
        let prj = con_s.d0(par);
        let sqdist = pnt.square_distance(&prj);
        if sqdist > tol_tol {
            tol_tol = sqdist;
        }
    }

    // `cxx:131-145`.
    tol = 1.00005 * tol_tol.sqrt();
    if tol >= init_tol {
        // `BB.UpdateEdge(anEdge, Tol)` and `BB.UpdateVertex(aVertex, Tol)` both
        // raise only (`BRep_TEdge.lxx:33-37`, `BRep_TVertex.lxx:33-37`).
        geometry.raise_edge_tolerance(an_edge, tol);
        for vertex in ShapeIterator::of_shape(an_edge) {
            if let Some(vertex) = Vertex::wrap(vertex) {
                geometry.raise_vertex_tolerance(&vertex.0, tol);
            }
        }
    }
    Ok(())
}
