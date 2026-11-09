//! `LocOpe_SplitShape::Add(const NCollection_List<TopoDS_Shape>& Lwires, const
//! TopoDS_Face& F)` (`LocOpe_SplitShape.cxx:256-598`) and
//! `Add(const TopoDS_Wire& W, const TopoDS_Face& F)` (`cxx:602-654`).
//!
//! `Add(Lwires, F)` walks the outer wire of the descendant face and chains the
//! section wires inserted at the vertices they meet, building the new boundary
//! wires of the pieces. The walk keeps two stacks: the vertices it broke at and
//! the wires those breaks belong to, so after a closed piece it resumes at the
//! next break instead of at the outer wire's first vertex.
//!
//! `Add(W, F)` is the single-wire entry point: it dispatches on `LocOpe::Closed`
//! and encodes `GeomAbs_CN` continuity on the section edges.

use crate::abs::{Orientation, ShapeType};
use crate::brep_lib_continuity::{has_continuity, set_stored_continuity};
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tools;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::locope::closed::closed_wire;
use crate::locope::shape_maps::{ShapeDataMap, ShapeValueMap};
use crate::locope::split_shape::direction::choose_direction;
use crate::locope::split_shape::is_inside::is_inside_wire;
use crate::locope::split_shape::{empty_copied_face, SplitShape};
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools::{explorer_shapes, is_same, wire_vertices};

impl SplitShape {
    /// `Add(W, F)` (`cxx:602-654`).
    ///
    /// The trailing `catch (Standard_Failure)` (`cxx:634-643`) has no
    /// counterpart: the two split routines answer `false` instead of throwing.
    pub fn add_wire(&mut self, w: &Wire, f: &Face) -> bool {
        // `cxx:608-610`.
        if self.done() {
            return false;
        }

        // `cxx:611-615`: an empty descendant list is rebuilt first.
        if self.map_of(&f.0).is_empty() {
            self.rebuild(&f.0);
        }

        // `cxx:616-632`.
        let ok = if !closed_wire(w, f) {
            self.add_open_wire(w, f)
        } else {
            self.add_closed_wire(w, f)
        };
        if !ok {
            return false;
        }

        // `cxx:645-652`: encode `GeomAbs_CN` regularity on the section edges.
        for sub in explorer_shapes(&w.0, ShapeType::Edge) {
            let Some(edg) = Edge::wrap(sub) else {
                continue;
            };
            if !has_continuity(&edg) {
                set_stored_continuity(&edg, f, f, 6);
            }
        }
        true
    }

    /// `Add(Lwires, F)` (`cxx:256-598`).
    ///
    /// The reference catches `Standard_Failure` one level up (`cxx:634-643`),
    /// so every place it would raise - `TopoDS::Face` / `TopoDS::Wire` on the
    /// wrong type, an unbound `myMap` / `VerSecMap` / `VerWireMap` lookup, an
    /// empty `BreakVertices` - answers `false` here.
    pub fn add_wires(&mut self, lwires: &[TopoShape], f: &Face) -> bool {
        // `cxx:259-262`.
        if self.done() {
            return false;
        }

        // `cxx:264-268`: `lf = myMap(F)`, rebuilt when empty. The reference
        // re-reads the list after `Rebuild` because the rebuild rewrites it.
        let mut lf = self.map_of(&f.0);
        if lf.is_empty() {
            self.rebuild(&f.0);
            lf = self.map_of(&f.0);
        }

        // `cxx:271-296`: the first descendant that contains any section wire.
        let mut inside = Vec::new();
        let mut found: Option<usize> = None;
        for (i, s) in lf.iter().enumerate() {
            let Some(fac) = Face::wrap(s.clone()) else {
                continue;
            };
            for w in lwires {
                let Some(wire) = Wire::wrap(w.clone()) else {
                    continue;
                };
                if is_inside_wire(&fac, &wire) {
                    inside.push(wire.0);
                }
            }
            if !inside.is_empty() {
                found = Some(i);
                break;
            }
        }
        let Some(found) = found else {
            return false;
        };

        // `cxx:299-302`: `FaceRef` is taken out of the descendant list.
        let Some(mut face_ref) = Face::wrap(lf[found].clone()) else {
            return false;
        };
        face_ref.0.set_orientation(Orientation::Forward);
        let mut new_lf: Vec<TopoShape> = lf[..found].to_vec();
        new_lf.extend(lf[found + 1..].iter().cloned());

        // `cxx:304-310`: each section is expected twice, once per side.
        let mut sections_times: ShapeValueMap<i32> = ShapeValueMap::new();
        for section in &inside {
            sections_times.bind(section, 2);
        }

        // `cxx:312-313`.
        let mut break_vertices: Vec<TopoShape> = Vec::new();
        let mut break_on_wires: Vec<TopoShape> = Vec::new();

        // `cxx:315-348`: the face wire the first vertex of every section lies
        // on. `aVer` in the reference starts as a null shape, so a wire without
        // vertices is never mistaken for a match.
        let mut ver_wire_map: ShapeDataMap = ShapeDataMap::new();
        for section in &inside {
            let (v1, v2) = wire_vertices(&Wire(section.clone()));
            for ver in [v1, v2] {
                let Some(ver) = ver else {
                    continue;
                };
                if ver_wire_map.is_bound(&ver.0) {
                    continue;
                }
                for wire in explorer_shapes(&face_ref.0, ShapeType::Wire) {
                    let mut hit = false;
                    for v in explorer_shapes(&wire, ShapeType::Vertex) {
                        if is_same(&v, &ver.0) {
                            hit = true;
                            break;
                        }
                    }
                    if hit {
                        ver_wire_map.bind(&ver.0, &wire);
                        break;
                    }
                }
            }
        }

        // `cxx:350-368`: every section, listed under both of its end vertices.
        let mut ver_sec_map: ShapeValueMap<Vec<TopoShape>> = ShapeValueMap::new();
        for section in &inside {
            let (v1, v2) = wire_vertices(&Wire(section.clone()));
            for ver in [v1, v2] {
                let Some(ver) = ver else {
                    continue;
                };
                if !ver_sec_map.is_bound(&ver.0) {
                    ver_sec_map.bind(&ver.0, Vec::new());
                }
                if let Some(list) = ver_sec_map.change(&ver.0) {
                    list.push(section.clone());
                }
            }
        }

        // `cxx:372-377`.
        let Some(outer) = brep_tools::outer_wire(&face_ref) else {
            return false;
        };
        let mut cur_wire = outer;
        let mut make_wire = MakeWire::new();
        let mut wexp = WireExplorer::of_wire_on_face(&cur_wire, &face_ref);
        let mut new_wires: Vec<TopoShape> = Vec::new();
        let builder = TopoBuilder::new();

        // `cxx:378-516`: the global walk.
        loop {
            // `cxx:380-384`: `theStartVertex` is only read after the closed
            // branch reassigns it, so the reference never uses this value.
            let Some(_start_vertex) = wexp.current_vertex().cloned() else {
                return false;
            };
            let Some(cur_edge) = wexp.current().cloned() else {
                return false;
            };
            let mut last_edge = cur_edge.clone();
            let _ = make_wire.add(&cur_edge);
            let mut section_wire: Option<Wire> = None;
            wexp.next();
            if !wexp.more() {
                wexp.init(&cur_wire, Some(&face_ref));
            }

            // `cxx:391-412`: extend the wire until it closes or a section
            // starts at the current vertex.
            let mut cur_vertex = None;
            loop {
                if make_wire.wire().is_closed() {
                    break;
                }
                let Some(vertex) = wexp.current_vertex().cloned() else {
                    return false;
                };
                cur_vertex = Some(vertex.clone());
                if ver_sec_map.is_bound(&vertex.0) {
                    // `cxx:399-403`: `ChooseDirection` picks the section that
                    // leaves `LastEdge` with the smallest left turn.
                    let sections = ver_sec_map.find(&vertex.0).cloned().unwrap_or_default();
                    let chosen = choose_direction(&last_edge.0, &vertex, &face_ref, &sections);
                    let Some(chosen) = chosen.and_then(Wire::wrap) else {
                        return false;
                    };
                    section_wire = Some(chosen);
                    break;
                }
                let Some(next_edge) = wexp.current().cloned() else {
                    return false;
                };
                let _ = make_wire.add(&next_edge);
                last_edge = next_edge;
                wexp.next();
                if !wexp.more() {
                    wexp.init(&cur_wire, Some(&face_ref));
                }
            }

            // `cxx:413-427`: a piece closed without using a section resumes at
            // the first pending break. `theStartVertex` is reassigned right
            // before the loop restarts, so the reference never reads it.
            if make_wire.wire().is_closed() {
                new_wires.push(make_wire.wire().0);
                let Some(next_break) = break_vertices.first().cloned() else {
                    return false;
                };
                break_vertices.remove(0);
                let Some(next_wire) = break_on_wires.first().cloned() else {
                    return false;
                };
                break_on_wires.remove(0);
                let Some(next_wire) = Wire::wrap(next_wire) else {
                    return false;
                };
                cur_wire = next_wire;
                wexp.init(&cur_wire, Some(&face_ref));
                // `cxx:422-425`: `wexp.Next()` runs until the break vertex is
                // reached; a wire without it would spin in OCCT, so the port
                // stops at the end of the round.
                loop {
                    let Some(vertex) = wexp.current_vertex() else {
                        return false;
                    };
                    if is_same(&vertex.0, &next_break) {
                        break;
                    }
                    wexp.next();
                    if !wexp.more() {
                        wexp.init(&cur_wire, Some(&face_ref));
                        break;
                    }
                }
                make_wire = MakeWire::new();
                continue;
            }

            // `cxx:429-432`.
            let Some(vertex) = cur_vertex else {
                return false;
            };
            break_vertices.push(vertex.0.clone());
            break_on_wires.push(cur_wire.0.clone());
            let mut break_vertex_now = vertex;
            let Some(mut section) = section_wire else {
                return false;
            };

            // `cxx:433-511`: follow sections until the wire closes.
            loop {
                let _ = make_wire.add_wire(&section);
                // `cxx:435-439`: `SectionsTimes(section)--`, unbound at zero.
                let Some(times) = sections_times.change(&section.0) else {
                    return false;
                };
                *times -= 1;
                if *times == 0 {
                    sections_times.unbind(&section.0);
                }

                if make_wire.wire().is_closed() {
                    // `cxx:440-459`.
                    new_wires.push(make_wire.wire().0);
                    if sections_times.is_empty() {
                        break;
                    }
                    let Some(next_break) = break_vertices.first().cloned() else {
                        return false;
                    };
                    break_vertices.remove(0);
                    let Some(next_wire) = break_on_wires.first().cloned() else {
                        return false;
                    };
                    break_on_wires.remove(0);
                    let Some(next_wire) = Wire::wrap(next_wire) else {
                        return false;
                    };
                    cur_wire = next_wire;
                    wexp.init(&cur_wire, Some(&face_ref));
                    loop {
                        let Some(v) = wexp.current_vertex() else {
                            return false;
                        };
                        if is_same(&v.0, &next_break) {
                            break;
                        }
                        wexp.next();
                        if !wexp.more() {
                            wexp.init(&cur_wire, Some(&face_ref));
                            break;
                        }
                    }
                    make_wire = MakeWire::new();
                    break;
                }

                // `cxx:461-511`.
                let (v1, v2) = wire_vertices(&section);
                let (Some(v1), Some(v2)) = (v1, v2) else {
                    return false;
                };
                let start_vertex = if is_same(&v1.0, &break_vertex_now.0) { v2 } else { v1 };
                let Some(wire) = ver_wire_map.find(&start_vertex.0).cloned() else {
                    return false;
                };
                let Some(wire) = Wire::wrap(wire) else {
                    return false;
                };
                cur_wire = wire;
                wexp.init(&cur_wire, Some(&face_ref));
                loop {
                    let Some(v) = wexp.current_vertex() else {
                        return false;
                    };
                    if is_same(&v.0, &start_vertex.0) {
                        break;
                    }
                    wexp.next();
                    if !wexp.more() {
                        wexp.init(&cur_wire, Some(&face_ref));
                        break;
                    }
                }

                // `cxx:472-477`.
                let Some(sections) = ver_sec_map.find(&start_vertex.0).cloned() else {
                    return false;
                };
                if sections.len() == 1 {
                    break;
                }

                // `cxx:479-480`: the other section on this vertex.
                let next_section = if is_same(&section.0, &sections[0]) {
                    sections.last().cloned()
                } else {
                    sections.first().cloned()
                };
                let Some(next_section) = next_section.and_then(Wire::wrap) else {
                    return false;
                };

                // `cxx:482-496`: how many sections leave vertices on this wire.
                // One means an inner touching wire, which the reference leaves
                // without a branch (`cxx:491-494`).
                let mut count = 0usize;
                for (_, value) in ver_wire_map.iter() {
                    if is_same(value, &cur_wire.0) {
                        count += 1;
                    }
                }
                if count != 1 {
                    // `cxx:497-508`.
                    let Some(start_edge) = wexp.current().cloned() else {
                        return false;
                    };
                    let ldirs = [start_edge.0.clone(), next_section.0.clone()];
                    let direction =
                        choose_direction(&section.0, &start_vertex, &face_ref, &ldirs);
                    if let Some(direction) = direction {
                        if is_same(&direction, &start_edge.0) {
                            break;
                        }
                    }
                }
                section = next_section;
                break_vertex_now = start_vertex;
            }

            // `cxx:512-515`.
            if sections_times.is_empty() {
                break;
            }
        }

        // `cxx:519-528`: one new face per rebuilt boundary wire.
        let mut new_faces: Vec<TopoShape> = Vec::new();
        for wire in &new_wires {
            let mut new_face = empty_copied_face(&face_ref);
            new_face.0.set_orientation(Orientation::Forward);
            builder.add(&mut new_face.0, wire);
            new_faces.push(new_face.0);
        }

        // `cxx:531-557`: the wires of `FaceRef` that no new wire uses are holes.
        let mut holes: Vec<TopoShape> = Vec::new();
        for wire in explorer_shapes(&face_ref.0, ShapeType::Wire) {
            let Some(hole_edge) = explorer_shapes(&wire, ShapeType::Edge).into_iter().next() else {
                continue;
            };
            let mut used = false;
            for new_wire in &new_wires {
                for edge in explorer_shapes(new_wire, ShapeType::Edge) {
                    if is_same(&edge, &hole_edge) {
                        used = true;
                        break;
                    }
                }
                if used {
                    break;
                }
            }
            if !used {
                holes.push(wire);
            }
        }

        // `cxx:559-572`: each hole goes into the new face that contains it.
        for hole in &holes {
            let Some(hole_wire) = Wire::wrap(hole.clone()) else {
                continue;
            };
            for i in 0..new_faces.len() {
                let Some(new_face) = Face::wrap(new_faces[i].clone()) else {
                    continue;
                };
                if is_inside_wire(&new_face, &hole_wire) {
                    builder.add(&mut new_faces[i], hole);
                    break;
                }
            }
        }

        // `cxx:575`: `lf.Append(NewFaces)`.
        new_lf.extend(new_faces);
        self.map_bind(&f.0, new_lf);

        // `cxx:578-582`: every wire of `F` loses its descendants.
        for wire in explorer_shapes(&f.0, ShapeType::Wire) {
            self.set_map_of(&wire, Vec::new());
        }

        // `cxx:585-595`: `GeomAbs_CN` regularity of the section edges.
        for section in &inside {
            for sub in explorer_shapes(section, ShapeType::Edge) {
                let Some(edg) = Edge::wrap(sub) else {
                    continue;
                };
                if !has_continuity(&edg) {
                    set_stored_continuity(&edg, f, f, 6);
                }
            }
        }
        true
    }
}
