//! The two wire-splitting steps of `ShapeFix_Face::Perform`'s second-part
//! loop: the file-local `SplitWire` (`ShapeFix_Face.cxx:242-341`), run from the
//! `NeedCheckSplitWire` block (`cxx:629-674`), and
//! `ShapeFix_Face::FixSplitFace` (`cxx:2905-3009`), run last (`cxx:711-716`).

use occt_core::precision::PCONFUSION;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::{FClass2d, FaceState};
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools_full::{edges_of_wire, is_same};

use super::face_geom_helpers::empty_copied_face;
use super::ShapeFixFace;
use crate::shape_fix_compose_shell::ReShape;

/// `SplitWire(face, wire, aResWires)` (`ShapeFix_Face.cxx:242-341`): re-chains
/// the wire's edges by shared vertices and returns every chain that closes
/// (`cxx:319-322`) or runs out of a connecting edge (`cxx:329-332`).
///
/// The C++ returns `true` on every path (`cxx:340`), so the chains are the
/// whole result.
pub(in crate::shhealing) fn split_wire(face: &Face, wire: &Wire) -> Vec<Wire> {
    // cxx:247-249.
    let edges: Vec<Edge> = edges_of_wire(wire);
    let nb = edges.len();
    let mut used: Vec<bool> = vec![false; nb];
    let mut res: Vec<Wire> = Vec::new();
    // cxx:250-337.
    for i in 0..nb {
        if used[i] {
            continue; // cxx:252-255
        }
        // cxx:257-261.
        let e1 = edges[i].clone();
        used[i] = true;
        let v0 = crate::shhealing::first_vertex(&e1);
        let mut v1 = crate::shhealing::last_vertex(&e1);
        let mut chain: Vec<Edge> = vec![e1.clone()];
        let mut is_connected_edge = true;
        let mut abandoned = false;
        // cxx:267-324.
        let mut j = 1usize;
        while j < nb && is_connected_edge {
            // cxx:269-285: rescan from index 2 for the first unused edge whose
            // first vertex is `V1`.
            let mut found: Option<usize> = None;
            for k in 1..nb {
                if used[k] {
                    continue;
                }
                let f2 = crate::shhealing::first_vertex(&edges[k]);
                if match (&f2, &v1) {
                    (Some(a), Some(b)) => is_same(&a.0, &b.0),
                    _ => false,
                } {
                    found = Some(k);
                    break;
                }
            }
            // cxx:287-291.
            let Some(k) = found else {
                is_connected_edge = false;
                break;
            };
            let e2 = edges[k].clone();
            chain.push(e2.clone());
            used[k] = true;
            v1 = crate::shhealing::last_vertex(&e2); // cxx:282
            // cxx:292-324: when the two chain endpoints coincide, the wire is
            // closed in 3D; the 2D pcurve ends then have to agree as well.
            let closed = match (&v1, &v0) {
                (Some(a), Some(b)) => is_same(&a.0, &b.0),
                _ => false,
            };
            if closed {
                // cxx:294-298: `BRep_Tool::CurveOnSurface(E1/E2, face, a/b)`.
                let c1 = crate::boptools_2d::curve_on_surface_oriented(&e1, face, false)
                    .map(|(pc, a, _b)| (pc, a));
                let c2 = crate::boptools_2d::curve_on_surface_oriented(&e2, face, false)
                    .map(|(pc, _a, b)| (pc, b));
                let (Some((curve1, mut a1)), Some((curve2, mut b2))) = (c1, c2) else {
                    abandoned = true; // cxx:299-302
                    break;
                };
                // cxx:305-312: a reversed edge is evaluated at its other end.
                if e1.0.orientation() == Orientation::Reversed {
                    a1 = curve1.last_parameter();
                }
                if e2.0.orientation() == Orientation::Reversed {
                    b2 = curve2.first_parameter();
                }
                let p0 = curve1.d0(a1); // cxx:313
                let p1 = curve2.d0(b2); // cxx:314
                // cxx:315-318: `GeomAdaptor_Surface::UResolution/VResolution`
                // times two, against the larger endpoint tolerance.
                let (Some(surf), Some(vv0), Some(vv1)) =
                    (BRepTool::face_surface(face), v0.as_ref(), v1.as_ref())
                else {
                    is_connected_edge = false;
                    break;
                };
                let tol = BRepTool::vertex_tolerance(vv0).max(BRepTool::vertex_tolerance(vv1));
                let res_u = occt_geom::approx_same_parameter::u_resolution(surf.as_ref(), tol);
                let res_v = occt_geom::approx_same_parameter::v_resolution(surf.as_ref(), tol);
                let max_resolution = 2.0 * res_u.max(res_v);
                if p0.square_distance(&p1) < max_resolution {
                    // cxx:319-322.
                    res.push(crate::shhealing::wire_from_wire_data(&chain, &[]));
                    break;
                }
            }
            j += 1;
        }
        // cxx:326-328.
        if abandoned {
            continue;
        }
        // cxx:329-332: a chain that ran out of a connecting edge is kept as an
        // open wire.
        if !is_connected_edge {
            res.push(crate::shhealing::wire_from_wire_data(&chain, &[]));
        }
        // cxx:333-336.
        if used.iter().all(|u| *u) {
            break;
        }
    }
    res
}

impl ShapeFixFace {
    /// `ShapeFix_Face::FixSplitFace(MapWires)` (`ShapeFix_Face.cxx:2905-3009`):
    /// turns each outer wire of `MapWires` into its own face, re-adding each
    /// nested wire inside the one that contains it. Returns true only when more
    /// than one face came out (`cxx:2986`), in which case `myResult` is a
    /// compound (`cxx:2988-2993`).
    ///
    /// `MapWires` is [`ShapeFixFace::map_wires`], filled by `FixOrientation`'s
    /// multi-wire nesting branch (`cxx:1276-1606`, `MapWires.Bind` at `cxx:1537`,
    /// `:1542`, `:1579`, `:1583`). A face with no nesting leaves the map empty,
    /// `NbWiresNew` stays below `NbWires` and the method returns false.
    pub fn fix_split_face(&mut self) -> bool {
        let builder = TopoBuilder::new();
        // cxx:2911-2915: `S = Context()->Apply(myFace)`.
        let face = match self.face.clone() {
            Some(f) => f,
            None => return false,
        };
        let applied = self.context.apply(&face.0);
        let s = if applied.is_face() { Face(applied) } else { face.clone() };

        // cxx:2916-2925.
        let (wires, _) = super::face_geom_helpers::split_face_children(&s);
        let mut nb_wires = 0usize;
        let mut nb_wires_new = 0usize;
        let mut faces: Vec<Face> = Vec::new();
        for wire in &wires {
            nb_wires += 1;
            let bound = self
                .map_wires
                .iter()
                .find(|(k, _)| is_same(&k.0, &wire.0))
                .map(|(_, v)| v.clone());
            let Some(int_wires) = bound else {
                continue; // cxx:2927
            };
            // cxx:2930-2934.
            let ents = edges_of_wire(wire);
            if ents.is_empty() {
                continue;
            }
            // cxx:2936-2947: only a wire closed in 3D can bound a face.
            let e1 = ents[0].clone();
            let e2 = ents[ents.len() - 1].clone();
            let v1 = crate::shhealing::first_vertex(&e1);
            let v2 = crate::shhealing::last_vertex(&e2);
            let same = match (&v1, &v2) {
                (Some(a), Some(b)) => is_same(&a.0, &b.0),
                _ => false,
            };
            if !same {
                return false; // cxx:2943-2946
            }
            // cxx:2948-2954.
            let mut tmp_face = empty_copied_face(&s);
            tmp_face.0.set_orientation(Orientation::Forward);
            builder.add(&mut tmp_face.0, &wire.0);
            nb_wires_new += 1;
            // cxx:2955-2972.
            for iw in &int_wires {
                let mut a_face = empty_copied_face(&tmp_face);
                a_face.0.set_orientation(Orientation::Forward);
                builder.add(&mut a_face.0, &iw.0);
                // `BRepTopAdaptor_FClass2d clas(aFace, PConfusion())` then
                // `PerformInfinitePoint()`. OCCT's constructor cannot fail; the
                // port's returns `Err` for a face whose UV definition is
                // degenerate, and that is folded into the `else` (Reversed)
                // arm rather than given an invented state.
                let inside = matches!(
                    FClass2d::new(&a_face, PCONFUSION),
                    Ok(clas) if clas.perform_infinite_point() == FaceState::In
                );
                if inside {
                    builder.add(&mut tmp_face.0, &iw.0); // cxx:2965-2967
                } else {
                    let mut rev = iw.clone();
                    rev.0.reverse();
                    builder.add(&mut tmp_face.0, &rev.0); // cxx:2968-2970
                }
                nb_wires_new += 1;
            }
            // cxx:2973-2976.
            if !self.my_fwd {
                tmp_face.0.set_orientation(Orientation::Reversed);
            }
            faces.push(tmp_face);
        }

        // cxx:2981-2984.
        if nb_wires != nb_wires_new {
            return false;
        }
        // cxx:2986.
        if faces.len() <= 1 {
            return false;
        }
        // cxx:2988-2993.
        let parts: Vec<TopoShape> = faces.iter().map(|f| f.0.clone()).collect();
        let comp: TopoShape = builder.make_compound_of(&parts).0;
        self.result = Some(comp.clone());
        self.context.replace(&face.0, &comp); // cxx:2996-2999
        // cxx:3000-3006: `myFace` ends as the last face of the compound.
        for f in crate::topo_tools_full::faces_of(&comp) {
            self.face = Some(f);
        }
        true
    }
}
