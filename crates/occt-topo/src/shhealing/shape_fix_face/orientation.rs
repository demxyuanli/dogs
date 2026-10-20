//! `ShapeFix_Face::FixOrientation` (`ShapeFix_Face.cxx:1165-1646`).

use super::prelude::*;
use super::*;

impl ShapeFixFace {
    /// `ShapeFix_Face::FixOrientation` (`ShapeFix_Face.cxx:1165-1646`).
    /// One wire is `cxx:1255-1274`; several wires are `cxx:1276-1606`.
    ///
    /// The single reversed wire that `StepToTopoDS_TranslateFace` stores for a
    /// `FACE_BOUND(..., .F.)` makes `ShapeAnalysis::IsOuterBound` negative, so
    /// the wire is reversed (`ShapeExtend_WireData::Reverse(myFace)`,
    /// `cxx:1266-1268`) and the face is rebuilt over it with orientation
    /// FORWARD (`cxx:1611-1637`). Without it the mesher's
    /// `ShapeExtend_WireData(wire, chained=true, manifold=false)` list is the
    /// reverse of the wire's stored order and the Delaunay frontier starves
    /// (`Shap*Fix*.top` torus faces).
    ///
    /// Several wires take `orient_several_wires` (`cxx:1276-1606`).
    pub fn fix_orientation(&mut self) -> bool {
        // cxx:1114-1116: `FixOrientation()` clears `MapWires` and delegates to
        // the overload Perform calls (`cxx:689-694`).
        self.map_wires.clear();
        let mut done = false; // cxx:1169.

        // cxx:1170-1176: `myFace = TopoDS::Face(Context()->Apply(myFace))`.
        let mut face = match self.face.clone() {
            Some(f) => f,
            None => return false,
        };
        {
            let applied = self.context.apply(&face.0);
            if applied.is_face() {
                face = Face(applied);
                self.face = Some(face.clone());
            }
        }

        // cxx:1180-1235: `TopoDS_Iterator(myFace, false)` splits the direct
        // children into `ws` (wires to rebuild with), `allSubShapes` (vertices
        // and unoriented children, re-added unchanged) and `VerySmallWires`.
        let mut ws: Vec<Wire> = Vec::new();
        let mut all_sub: Vec<TopoShape> = Vec::new();
        let mut very_small = 0usize;
        let children: Vec<TopoShape> = {
            let ts = face.0.tshape.read().unwrap();
            ts.children.clone()
        };
        for child in children {
            let ori = child.orientation();
            // cxx:1182-1188.
            if child.shape_type() == ShapeType::Vertex
                || (ori != Orientation::Forward && ori != Orientation::Reversed)
            {
                all_sub.push(child);
                continue;
            }
            // cxx:1191-1226: a wire holding a single edge is measured by a
            // 10-point polyline of `ShapeAnalysis_Edge::Curve3d`; a multi-edge
            // wire keeps `RealLast()` and is always "long".
            let wire_children: Vec<TopoShape> = {
                let ts = child.tshape.read().unwrap();
                ts.children.clone()
            };
            let length = if wire_children.is_empty() {
                0.0 // cxx:1222-1225.
            } else if wire_children.len() == 1 {
                let edge = Edge(wire_children[0].clone());
                let mut len = 0.0;
                if !BRepTool::is_degenerated(&edge) {
                    if let Some(c3d) = BRepTool::edge_curve(&edge) {
                        let (mut first, mut last) = BRepTool::edge_parameters(&edge);
                        // `ShapeAnalysis_Edge::Curve3d(..., orient=true)`.
                        if wire_children[0].orientation().is_reversed() {
                            std::mem::swap(&mut first, &mut last);
                        }
                        const NB_CONTROL: usize = 10;
                        let mut prev = c3d.d0(first);
                        for j in 1..NB_CONTROL {
                            let prm = ((NB_CONTROL - 1 - j) as f64 * first
                                + j as f64 * last)
                                / (NB_CONTROL - 1) as f64;
                            let curr = c3d.d0(prm);
                            len += curr.distance(&prev);
                            prev = curr;
                        }
                    }
                }
                len
            } else {
                f64::MAX // cxx:1192 `RealLast()`.
            };
            // cxx:1227-1234.
            if length > CONFUSION {
                ws.push(Wire(child.clone()));
                all_sub.push(child);
            } else {
                very_small += 1;
            }
        }
        // cxx:1236-1238.
        if very_small > 0 {
            done = true;
        }

        let nb = ws.len();
        let nb_all = all_sub.len();
        // cxx:1243-1250: no usable wire, nothing to rebuild.
        if nb == 0 {
            return false;
        }

        let is_add_natural_bounds = self.is_need_add_natural_bound(&ws); // cxx:1252
        let mut nb_reversed = 0usize; // `aSeqReversed`

        // cxx:1255-1274.
        if nb == 1 {
            let builder = TopoBuilder::new();
            let mut af = match BRepTool::face_surface(&face) {
                Some(surf) => builder.make_face(surf, &[ws[0].clone()]), // cxx:1258-1261
                None => Face::new(),
            };
            af.0.set_orientation(Orientation::Forward); // cxx:1260
            if !is_add_natural_bounds && !is_outer_bound(&af) {
                // cxx:1266: `ShapeExtend_WireData(ws.Value(1))` -> the oriented
                // edge list (`chained=true` keeps `TopoDS_Iterator` order, a
                // REVERSED wire prepends), non-manifold edges kept aside.
                let (mut data, non_manifold): (Vec<Edge>, Vec<Edge>) =
                    edges_of_wire(&ws[0]).into_iter().partition(|e| {
                        let o = e.0.orientation();
                        o == Orientation::Forward || o == Orientation::Reversed
                    });
                if ws[0].0.orientation().is_reversed() {
                    data.reverse(); // `ShapeExtend_WireData.cxx:115-122`
                }
                // cxx:1267: `sbdw->Reverse(myFace)`.
                reverse_wire_data_on_face(&mut data, &face);
                // cxx:1268: `ws.SetValue(1, sbdw->Wire())`.
                ws[0] = wire_from_wire_data(&data, &non_manifold);
                nb_reversed += 1;
                done = true;
            }
        }
        // cxx:1276-1606: several wires, classify each against the others.
        if nb != 1 {
            if let Some(surf) = self.surf.clone() {
                let nrev = super::super::shape_fix_face_wires::orient_several_wires(
                    &face,
                    surf.as_ref(),
                    &mut ws,
                    &all_sub,
                    &mut self.map_wires,
                );
                if nrev > 0 {
                    nb_reversed += nrev;
                    done = true;
                }
            }
        }

        // cxx:1607-1609.
        if is_add_natural_bounds && nb == nb_reversed {
            done = false;
        }

        // cxx:1610-1646: rebuild the face from `ws` plus the untouched children.
        if done {
            let builder = TopoBuilder::new();
            let mut new_face = match BRepTool::face_surface(&face) {
                Some(surf) => builder.make_face(surf, &ws),
                None => Face::new(),
            };
            new_face.0.set_orientation(Orientation::Forward); // cxx:1613
            // `myFace.EmptyCopied()` (`cxx:1612`) carries the location and the
            // whole `BRep_TFace` (surface, tolerance, natural restriction).
            new_face.0.set_location(face.0.location());
            if let Some(geom) = GeometryRegistry::global().face_geom(&face.0) {
                GeometryRegistry::global().set_face(&new_face.0, geom);
            }
            // cxx:1620-1633.
            if nb < nb_all {
                for s2 in &all_sub {
                    let ori = s2.orientation();
                    if s2.shape_type() != ShapeType::Wire
                        || (ori != Orientation::Forward && ori != Orientation::Reversed)
                    {
                        builder.add(&mut new_face.0, s2);
                    }
                }
            }
            // cxx:1634-1636.
            if !self.fwd {
                new_face.0.set_orientation(Orientation::Reversed);
            }
            // cxx:1639-1641: `Context()->Replace(myFace, S)`.
            self.context.replace(&face.0, &new_face.0);
            self.face = Some(new_face);
        }
        done
    }
}
