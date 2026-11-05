//! `PrepareSpine` (`BRepFill_OffsetWire.cxx:1195-1290`) with the statics it
//! uses: `CutEdge` (`cxx:1930-2070`), `CutCurve` (`cxx:2076-2125`) and
//! `EdgeVertices` (`cxx:127-136`).
//!
//! `PrepareSpine` rebuilds the spine as a face over the *same* surface and
//! location, whose wires are the input wires with each edge cut at the extrema
//! of curvature and at the inflexion points (`MAT2d_CutCurve`), so that the
//! bisecting locus sees monotone pieces. `myMapSpine` records which original
//! edge — and which original vertices — each piece came from, which
//! `GeneratedShapes` (`cxx:390-444`) later folds into `myMap`.
//!
//! `CutEdge` answers the number of pieces and leaves `Cuts` empty when the edge
//! is not modified. `ForceCut` overrides `MAT2d_CutCurve`:
//! * `0` (`cxx:1968-1984`): no override,
//! * `1` (`cxx:1986-1999`): two parts even when the curve is unmodified,
//! * `2` (`cxx:2001-2033`): three parts even when the curve is unmodified.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_fill_offset_wire::BRepFillOffsetWire;
use crate::brep_tool::BRepTool;
use occt_geom2d::curve::Curve2d;
use crate::builder::TopoBuilder;
use crate::mat2d::Mat2dCutCurve;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{
    edge_vertices, edges_of, is_same, shape_map_key, wire_vertices, wires_of_face,
};

impl BRepFillOffsetWire {
    /// `PrepareSpine()` (`cxx:1197-1290`): cut the spine edges, filling
    /// `myWorkSpine` and `myMapSpine`.
    pub fn prepare_spine(&mut self) {
        let builder = TopoBuilder::new();
        let spine = self
            .my_spine
            .clone()
            .expect("BRepFill_OffsetWire::PrepareSpine: null spine");

        self.my_work_spine = None;
        self.my_map_spine.clear();

        // `cxx:1206-1208`: the face keeps the spine's surface, location and
        // tolerance, with no wire.
        let surface = BRepTool::face_surface(&spine).expect("BRep_Tool::Surface(mySpine)");
        let tol_f = BRepTool::face_tolerance(&spine);
        let mut work = builder.make_face_with_tolerance(surface, &[], tol_f);
        work.0.set_location(spine.0.location());

        // `cxx:1210`: `TopoDS_Iterator IteF(mySpine)`.
        for wire in wires_of_face(&spine) {
            // `cxx:1214-1215`: `B.MakeWire(NW)` — an empty wire.
            let mut new_wire = builder.make_wire(&[]);

            // `cxx:1217-1228`.
            let mut forced_cut = 0i32;
            let mut nb_res_edges = -1i32;
            let edges = edges_of(&wire.0);
            let nb_edges = edges.len() as i32;

            if nb_edges == 1 && !self.my_is_open_result {
                forced_cut = 2;
            }

            // `cxx:1230`: `TopoDS_Iterator IteW(IteF.Value())`.
            for edge in &edges {
                // `cxx:1235-1237`.
                let (v1, v2) = edge_vertices(edge);
                if let Some(v) = &v1 {
                    self.bind_spine_vertex(&v.0);
                }
                if let Some(v) = &v2 {
                    self.bind_spine_vertex(&v.0);
                }

                // `cxx:1243-1248`.
                let num_curves_in_edge = num_curves_in_edge(&edge.0);
                if nb_edges == 2 && nb_res_edges == 0 && num_curves_in_edge > 1 {
                    forced_cut = 1;
                }

                // `cxx:1249`: the cut works on the forward copy.
                let forward = Edge(edge.0.oriented(Orientation::Forward));
                let mut cuts: Vec<TopoShape> = Vec::new();
                nb_res_edges = cut_edge(&forward, &spine, forced_cut, &mut cuts);

                if cuts.is_empty() {
                    // `cxx:1251-1255`.
                    builder.add_edge(&mut new_wire, edge);
                    self.my_map_spine
                        .insert(shape_map_key(&edge.0), (edge.0.clone(), edge.0.clone()));
                } else {
                    // `cxx:1257-1275`.
                    for cut in cuts {
                        let mut ne = Edge(cut);
                        ne.0.set_orientation(edge.0.orientation());
                        builder.add_edge(&mut new_wire, &ne);
                        self.my_map_spine
                            .insert(shape_map_key(&ne.0), (ne.0.clone(), edge.0.clone()));

                        let (cv1, cv2) = edge_vertices(&ne);
                        if let Some(v) = cv1 {
                            self.bind_spine_vertex_if_absent(&v.0, &edge.0);
                        }
                        if let Some(v) = cv2 {
                            self.bind_spine_vertex_if_absent(&v.0, &edge.0);
                        }
                    }
                }
            }

            // `cxx:1278-1284`: `NW.Closed(aV1.IsSame(aV2))`.
            let (a_v1, a_v2) = wire_vertices(&wire);
            let closed = match (&a_v1, &a_v2) {
                (Some(a), Some(b)) => is_same(&a.0, &b.0),
                _ => false,
            };
            new_wire.0.set_closed(closed);

            builder.add_wire(&mut work, &new_wire);
        }

        self.my_work_spine = Some(work);
    }

    /// `myMapSpine.Bind(V, V)` (`cxx:1236-1237`).
    fn bind_spine_vertex(&mut self, v: &TopoShape) {
        self.my_map_spine
            .insert(shape_map_key(v), (v.clone(), v.clone()));
    }

    /// `if (!myMapSpine.IsBound(V)) myMapSpine.Bind(V, E)` (`cxx:1264-1272`).
    fn bind_spine_vertex_if_absent(&mut self, v: &TopoShape, e: &TopoShape) {
        let key = shape_map_key(v);
        self.my_map_spine
            .entry(key)
            .or_insert_with(|| (v.clone(), e.clone()));
    }
}

/// `BRep_TEdge::Curves().Length()` — the number of curve representations of the
/// edge: the 3D curve (whose representation exists even when the curve handle
/// is null) plus one per pcurve.
fn num_curves_in_edge(e: &TopoShape) -> usize {
    let ts = e.tshape.read().expect("poisoned TShape lock");
    let mut n = 0usize;
    if let Some(core) = ts.edge_core() {
        // `BRep_Curve3D` is in the list when it was installed, either with a
        // curve or as the null representation (`tshape.rs` `curve3d_null`).
        if core.curve.is_some() || core.curve3d_null {
            n += 1;
        }
    }
    if let Some(pc) = ts.edge_pcurves() {
        n += pc.curves.iter().map(|(_, cs)| cs.len()).sum::<usize>();
    }
    n
}

/// `CutEdge(E, F, ForceCut, Cuts)` (`cxx:1930-2070`): the number of pieces,
/// appended to `cuts` (empty when the edge is not cut).
fn cut_edge(e: &Edge, f: &Face, force_cut: i32, cuts: &mut Vec<TopoShape>) -> i32 {
    cuts.clear();
    let mut cuter = Mat2dCutCurve::new();

    // `cxx:1950-1956`: the p-curve range is read first, then *overwritten* by
    // the 3D curve range, which is what the trimmed curve below is built over.
    let (c2d, _, _) = curve_on_surface_range(e, f)
        .unwrap_or_else(|| panic!("BRep_Tool::CurveOnSurface raised"));
    let a_c = BRepTool::edge_curve(e).unwrap_or_else(|| panic!("BRep_Tool::Curve raised"));
    let (f_par, l_par) = BRepTool::edge_parameters(e);
    let ct2d = Geom2dTrimmedCurve::new(c2d, f_par, l_par);

    // `cxx:1958-1961`: a circle spanning more than half a turn is not cut.
    if ct2d.basis_curve().is_geom2d_circle() && (f_par - l_par).abs() >= std::f64::consts::PI {
        return 0;
    }

    // `cxx:1966`: `Cuter.Perform(CT2d)`.
    cuter.perform(&ct2d);

    let mut the_curves: Vec<Geom2dTrimmedCurve> = Vec::new();
    match force_cut {
        0 => {
            // `cxx:1968-1984`.
            if cuter.un_modified() {
                return 0;
            }
            for k in 1..=cuter.nb_curves() {
                the_curves.push(cuter.value(k as i32).clone());
            }
        }
        1 => {
            // `cxx:1986-1999`.
            if cuter.un_modified() {
                cut_curve(&ct2d, 2, &mut the_curves);
            } else {
                for k in 1..=cuter.nb_curves() {
                    the_curves.push(cuter.value(k as i32).clone());
                }
            }
        }
        2 => {
            // `cxx:2001-2033`.
            if cuter.un_modified() {
                cut_curve(&ct2d, 3, &mut the_curves);
            } else if cuter.nb_curves() == 2 {
                let cc = cuter.value(1);
                if cc.last_parameter() > (l_par + f_par) / 2.0 {
                    cut_curve(cc, 2, &mut the_curves);
                    the_curves.push(cuter.value(2).clone());
                } else {
                    the_curves.push(cc.clone());
                    cut_curve(cuter.value(2), 2, &mut the_curves);
                }
            } else {
                for k in 1..=cuter.nb_curves() {
                    the_curves.push(cuter.value(k as i32).clone());
                }
            }
        }
        _ => panic!("CutEdge: unexpected ForceCut {force_cut}"),
    }

    // `cxx:2039-2068`: one edge per piece; only the last one keeps the original
    // end vertex, the intermediate vertices sit at the piece's last parameter.
    let a_tol = BRepTool::edge_tolerance(e);
    let builder = TopoBuilder::new();
    let reg = GeometryRegistry::global();
    let (v1, v2) = edge_vertices(e);
    let v1 = v1.expect("CutEdge: edge without a first vertex");
    let v2 = v2.expect("CutEdge: edge without a second vertex");
    let mut vf = v1;

    for (k, cc) in the_curves.iter().enumerate() {
        let vl = if k + 1 == the_curves.len() {
            v2.clone()
        } else {
            // `cxx:2051-2057`: `BRepLib_MakeVertex(P)` + `B.UpdateVertex(VL, aTol)`.
            let p = a_c.d0(cc.last_parameter());
            builder.make_vertex(GpPnt::new(p.x(), p.y(), p.z()), a_tol)
        };

        // `cxx:2059-2067`: `E.EmptyCopied()` keeps the curve representations,
        // which is what `B.Range` then re-trims to the piece.
        let ne = e.0.empty_copied();
        if let Some(g) = reg.edge_geom(e) {
            reg.set_edge(&ne, g);
        }
        reg.copy_edge_pcurve_slots(&ne, &e.0);
        let mut ne = ne;
        ne.set_orientation(Orientation::Forward);
        builder.add(&mut ne, &vf.0.oriented(Orientation::Forward));
        builder.add(&mut ne, &vl.0.oriented(Orientation::Reversed));
        reg.set_ranges_all(&ne, cc.first_parameter(), cc.last_parameter());
        cuts.push(ne.oriented(e.0.orientation()));

        vf = vl;
    }

    the_curves.len() as i32
}

/// `CutCurve(C, nbParts, theCurves)` (`cxx:2076-2125`): split `C` into
/// `nbParts` pieces at the parameters that are neither at the start nor
/// already at a coincident point.
fn cut_curve(c: &Geom2dTrimmedCurve, nb_parts: i32, the_curves: &mut Vec<Geom2dTrimmedCurve>) {
    let ptol = PCONFUSION * 10.0;
    let tol = CONFUSION * 10.0;
    let mut ya_cut = false;

    let mut uf = c.first_parameter();
    let ul = c.last_parameter();
    let mut pf = c.d0(uf);
    let pl = c.d0(ul);
    let step = (ul - uf) / nb_parts as f64;

    for i in 1..nb_parts {
        let uc = uf + i as f64 * step;
        let pc = c.d0(uc);

        if uc - uf > ptol && pc.distance(&pf) > tol {
            if ul - uc < ptol || pl.distance(&pc) < tol {
                continue;
            }
            the_curves.push(Geom2dTrimmedCurve::new(
                Arc::new(c.clone()),
                uf,
                uc,
            ));
            uf = uc;
            pf = pc;
            ya_cut = true;
        }
    }

    if ya_cut {
        the_curves.push(Geom2dTrimmedCurve::new(Arc::new(c.clone()), uf, ul));
    } else {
        the_curves.push(c.clone());
    }
}


