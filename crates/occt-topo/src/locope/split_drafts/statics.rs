//! The four file-local helpers of `LocOpe_SplitDrafts`
//! (`LocOpe_SplitDrafts.cxx:1461-1976`).
//!
//! * [`new_plane`] (`:1461-1538`): the inclined draft plane and its normal axis.
//! * [`make_face`] (`:1542-1713`): close a bag of edges into a face on `F`'s
//!   surface, inventing the missing iso p-curves.
//! * [`contains`] (`:1716-1725`): `IsSame` membership.
//! * [`new_edge`] (`:1731-1976`): the section edge of `F`'s surface with a draft
//!   plane, between two given vertices.

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx2d, GpDir, GpDir2d, GpLin, GpPln, GpPnt, GpVec2d};
use occt_core::precision::{Precision, CONFUSION, SQUARE_CONFUSION};
use occt_geom::intana::plane_plane_line;
use occt_geom::extrema_pc::extrema_ext_pc_range_done;
use occt_geom::{Curve, Surface};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::Geom2dLine;

use crate::abs::Orientation;
use crate::brep_gprop_full::surface_properties;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::geom_int::IntSS;
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::tgeometry::{FaceGeom, GeometryRegistry};
use crate::topo_tools::{edge_vertices, is_same};

/// `Contains` (`cxx:1716-1725`).
pub(crate) fn contains(list: &[TopoShape], s: &TopoShape) -> bool {
    list.iter().any(|x| is_same(x, s))
}

/// `TopoDS_Shape::EmptyCopied()` on an edge together with the geometry
/// `BRep_TEdge::EmptyCopy` (`BRep_TEdge.cxx:105-133`) carries over: curve,
/// range, tolerance, flags and every `CurveOnSurface` representation - only the
/// children are dropped (`cxx:798`, `cxx:1694`).
pub(crate) fn empty_copied_edge(e: &Edge) -> Edge {
    let registry = GeometryRegistry::global();
    let copy = Edge::new();
    if let Some(geom) = registry.edge_geom(&e.0) {
        registry.set_edge(&copy.0, geom);
    }
    copy
}

/// `TopoDS_Shape::EmptyCopied()` on a face with what `BRep_TFace::EmptyCopy`
/// (`BRep_TFace.cxx:36-43`) carries over: surface, location, tolerance.
pub(crate) fn empty_copied_face(f: &Face) -> TopoShape {
    let registry = GeometryRegistry::global();
    let copy = f.0.empty_copied();
    if let Some(geom) = registry.face_geom(&f.0) {
        let mut new_geom = FaceGeom::new(geom.surface);
        new_geom.tolerance = geom.tolerance;
        new_geom.natural_restriction = false;
        registry.set_face(&copy, new_geom);
    }
    copy
}

/// `NewPlane` (`cxx:1461-1538`): `Newpl` is the plane of `F` rotated about the
/// line where it meets the neutral plane so that the draft direction makes
/// `Ang` with it. `Modify == false` keeps the plane and only orients the normal.
pub(crate) fn new_plane(
    f: &TopoShape,
    extr: &GpDir,
    neutr: &GpPln,
    ang: f64,
    newpl: &mut GpPln,
    normal_f: &mut GpAx1,
    modify: bool,
) -> bool {
    let Some(surface) = BRepTool::face_surface(&Face(f.clone())) else {
        return false;
    };
    // `cxx:1476-1479`.
    let surface: Arc<dyn Surface> = match surface.rectangular_trimmed_basis() {
        Some(basis) => basis,
        None => surface,
    };
    let Some(plorig) = surface.gp_pln() else {
        return false;
    };

    // `cxx:1482-1493`.
    if !modify {
        *newpl = plorig;
        *normal_f = newpl.axis();
        let oris = f.orientation();
        if (newpl.is_direct() && oris == Orientation::Reversed)
            || (!newpl.is_direct() && oris == Orientation::Forward)
        {
            normal_f.reverse();
        }
        return true;
    }

    // `cxx:1496-1500`.
    let Some(line) = plane_plane_line(&plorig, neutr, Precision::ANGULAR, CONFUSION) else {
        return false;
    };

    // `cxx:1502-1515`.
    let nx = line.direction();
    *normal_f = plorig.axis();
    let ny = normal_f
        .direction()
        .crossed(&nx)
        .unwrap_or_else(|_| GpDir::default_dir());
    let a = extr.dot(&nx);

    if a.abs() <= 1.0 - Precision::ANGULAR {
        let mut b = extr.dot(&ny);
        let mut c = extr.dot(normal_f.direction());
        let oris = f.orientation();
        if (plorig.is_direct() && oris == Orientation::Reversed)
            || (!plorig.is_direct() && oris == Orientation::Forward)
        {
            b = -b;
            c = -c;
            normal_f.reverse();
        }
        // `cxx:1519-1532`.
        let denom = (1.0 - a * a).sqrt();
        let sina = ang.sin();
        if denom > sina.abs() {
            let phi = (b / denom).atan2(c / denom);
            let theta0 = (sina / denom).acos();
            let mut theta = theta0 - phi;
            if theta.cos() < 0.0 {
                theta = -theta0 - phi;
            }
            let axe = *line.position();
            *newpl = plorig.rotated(&axe, theta);
            return true;
        }
    }
    false
}

/// `MakeFace` (`cxx:1542-1713`): give every edge of `ledg` a p-curve on `f`
/// (missing ones are `U`-isos, reconstructed from a neighbour's p-curve), chain
/// the edges into closed wires and add the single resulting wire to `f`.
pub(crate) fn make_face(f: &mut TopoShape, ledg: &mut Vec<TopoShape>) {
    let b = TopoBuilder::new();
    let registry = GeometryRegistry::global();
    let face_key = GeometryRegistry::shape_key(f);

    // `cxx:1549-1611`: the p-curves that are missing.
    for i in 0..ledg.len() {
        let edg = Edge(ledg[i].clone());
        if registry.edge_pcurve(&edg.0, face_key).is_some() {
            continue;
        }
        let (f_rng, l_rng) = BRepTool::edge_parameters(&edg);
        let (v1, v2) = edge_vertices(&edg);
        let (Some(v1), Some(v2)) = (v1, v2) else {
            continue;
        };
        for j in 0..ledg.len() {
            if j == i || is_same(&ledg[j], &ledg[i]) {
                continue;
            }
            let edg2 = Edge(ledg[j].clone());
            let (vp1, vp2) = edge_vertices(&edg2);
            let (Some(vp1), Some(vp2)) = (vp1, vp2) else {
                continue;
            };
            if !(is_same(&vp1, &v1) || is_same(&vp2, &v1) || is_same(&vp1, &v2) || is_same(&vp2, &v2)) {
                continue;
            }
            let (f2, l2) = BRepTool::edge_parameters(&edg2);
            let Some(c22d) = registry.edge_pcurve(&edg2.0, face_key) else {
                continue;
            };
            // `cxx:1580-1603`: the 2D point of the shared vertex, shifted by the
            // other edge's parameter so the iso keeps the right constant.
            let mut pt2d = if is_same(&vp1, &v1) {
                c22d.d0(f2)
            } else if is_same(&vp2, &v1) {
                c22d.d0(l2)
            } else if is_same(&vp1, &v2) {
                c22d.d0(f2)
            } else {
                c22d.d0(l2)
            };
            let delta = if is_same(&vp1, &v1) || is_same(&vp2, &v1) {
                f_rng
            } else {
                l_rng
            };
            pt2d.set_y(pt2d.y() - delta);
            // `cxx:1601-1604`: `Geom2d_Line(pt2d, gp::DY2d())`.
            let c2d: Arc<dyn Curve2d> = Arc::new(Geom2dLine::new(GpAx2d::new(
                pt2d,
                GpDir2d::new(0.0, 1.0).expect("gp::DY2d"),
            )));
            registry.update_edge_pcurve_on_surface(
                &edg.0,
                &match BRepTool::face_surface(&Face(f.clone())) {
                    Some(s) => s,
                    None => return,
                },
                Some(c2d),
                BRepTool::edge_tolerance(&edg),
            );
            break;
        }
    }

    // `cxx:1614-1687`: chain the edges into closed wires.
    let mut lwires: Vec<Wire> = Vec::new();
    while !ledg.is_empty() {
        let mut wnew_edges: Vec<Edge> = Vec::new();
        let first = ledg.remove(0);
        let edg = Edge(first.clone());
        // `cxx:1626-1636`: `TopExp::Vertices` composes the edge orientation with
        // the children's, and the `REVERSED` branch swaps the two outputs, so
        // `VFirst` / `VLast` always come out as the stored `FORWARD` /
        // `REVERSED` children.
        let (Some(vf0), Some(vl0)) = edge_vertices(&edg) else {
            wnew_edges.push(edg);
            lwires.push(b.make_wire(&wnew_edges));
            continue;
        };
        wnew_edges.push(edg);
        let mut vf = vf0.0;
        let mut vl = vl0.0;
        let mut wdone = ledg.is_empty() || is_same(&vf, &vl);
        while !wdone {
            let mut hit: Option<(usize, Orientation)> = None;
            for (j, cand) in ledg.iter().enumerate() {
                let forward = cand.oriented(Orientation::Forward);
                let e2 = Edge(forward);
                let (Some(a), Some(c)) = edge_vertices(&e2) else {
                    continue;
                };
                let (a, c) = (a.0, c.0);
                // `cxx:1648-1669`, in the reference's order.
                if is_same(&a, &vl) {
                    vl = c;
                    hit = Some((j, Orientation::Forward));
                    break;
                } else if is_same(&c, &vf) {
                    vf = a;
                    hit = Some((j, Orientation::Forward));
                    break;
                } else if is_same(&a, &vf) {
                    vf = c;
                    hit = Some((j, Orientation::Reversed));
                    break;
                } else if is_same(&c, &vl) {
                    vl = a;
                    hit = Some((j, Orientation::Reversed));
                    break;
                }
            }
            match hit {
                None => wdone = true,
                Some((j, oredg)) => {
                    let taken = ledg.remove(j);
                    wnew_edges.push(Edge(taken.oriented(oredg)));
                    wdone = ledg.is_empty() || is_same(&vf, &vl);
                }
            }
        }
        lwires.push(b.make_wire(&wnew_edges));
    }

    // `cxx:1690-1712`.
    f.set_orientation(Orientation::Forward);
    for wire in &mut lwires {
        let mut new_face = empty_copied_face(&Face(f.clone()));
        b.add(&mut new_face, &wire.0);
        if let Ok((props, _area)) = surface_properties(&new_face) {
            if props.dim < 0.0 {
                wire.0.reverse();
            }
        }
    }
    if lwires.len() == 1 {
        let w = lwires[0].clone();
        b.add(f, &w.0);
    }
}

/// `NewEdge` (`cxx:1731-1976`): the intersection of `F`'s surface with the plane
/// `new_s`, trimmed between `v1` and `v2`. `None` is the reference's null edge.
pub(crate) fn new_edge(
    edg: &Edge,
    face: &TopoShape,
    new_s: &Arc<dyn Surface>,
    v1: &Vertex,
    v2: &Vertex,
) -> Option<Edge> {
    let s1 = BRepTool::face_surface(&Face(face.clone()))?;
    let app_s1 = s1.gp_pln().is_none();

    // `cxx:1745-1749`.
    let mut i2s = IntSS::new();
    i2s.perform(s1.clone(), new_s.clone(), CONFUSION, true, app_s1, false);
    if !i2s.is_done() || i2s.nb_lines() <= 0 {
        return None;
    }

    // `cxx:1753-1816`: locate `v1` and `v2` on the intersection lines.
    let pvf = BRepTool::vertex_point(v1);
    let pvl = BRepTool::vertex_point(v2);
    let mut prmf: f64 = 0.0;
    let mut prml: f64 = 0.0;
    let mut hit: Option<i32> = None;
    for i in 1..=i2s.nb_lines() {
        let Some(line) = i2s.line(i).map(|l| l.curve.clone()) else {
            continue;
        };
        let (first, last) = (line.first_parameter(), line.last_parameter());
        let (sols, done) = extrema_ext_pc_range_done(line.as_ref(), &pvf, first, last);
        if !done {
            continue;
        }
        let d_last = pvf.square_distance(&line.d0(last));
        let d_first = pvf.square_distance(&line.d0(first));
        let mut dist2min = d_first;
        let mut thepmin = first;
        if d_last < dist2min && !line.is_periodic() {
            dist2min = d_last;
            thepmin = last;
        }
        for sol in &sols {
            if sol.sq_dist < dist2min {
                dist2min = sol.sq_dist;
                thepmin = sol.u;
            }
        }
        if dist2min > SQUARE_CONFUSION {
            continue;
        }
        prmf = thepmin;

        // `cxx:1788-1814`: the same search for the far vertex.
        let (sols2, done2) = extrema_ext_pc_range_done(line.as_ref(), &pvl, first, last);
        if !done2 {
            continue;
        }
        let d_first2 = pvl.square_distance(&line.d0(first));
        let d_last2 = pvl.square_distance(&line.d0(last));
        // `cxx:1797-1801`: `TrimmedSquareDistances(Dist2, Dist2Min, ...)` puts
        // the *last* endpoint distance into `Dist2Min`, and the `!IsClosed`
        // switch is towards the *first* parameter.
        let mut dist2min2 = d_last2;
        let mut thepmin2 = last;
        if d_first2 < dist2min2 && !line.is_closed() {
            dist2min2 = d_first2;
            thepmin2 = first;
        }
        for sol in &sols2 {
            if sol.sq_dist < dist2min2 {
                dist2min2 = sol.sq_dist;
                thepmin2 = sol.u;
            }
        }
        if dist2min2 <= SQUARE_CONFUSION {
            prml = thepmin2;
            hit = Some(i);
            break;
        }
    }
    let i = hit?;

    // `cxx:1819-1946`.
    let mut rev = false;
    let mut vf = v1.clone();
    let mut vl = v2.clone();
    let cimg: Arc<dyn Curve> = i2s.line(i)?.curve.clone();
    let mut cimg2d: Option<Arc<dyn Curve2d>> = if app_s1 {
        i2s.line_on_s1(i).cloned()
    } else {
        None
    };

    if cimg.is_periodic() {
        let period = cimg.period();
        let imf = cimg.first_parameter();
        let iml = cimg.last_parameter();
        let (ef, el) = BRepTool::edge_parameters(edg);
        let delt = el - ef;
        let delt1 = (prml - prmf).abs();
        let delt2 = (period - delt1).abs();
        if delt1 == 0.0 || delt2 == 0.0 {
            prmf = imf;
            prml = iml;
        } else if (delt1 - delt).abs() > (delt2 - delt).abs() {
            if prml > prmf {
                if prml < iml {
                    prmf += period;
                } else {
                    prml -= period;
                }
            } else if prmf < iml {
                prml += period;
            } else {
                prmf -= period;
            }
        } else if (delt1 - delt).abs() < (delt2 - delt).abs() {
            if prmf >= iml && prml >= iml {
                prmf -= period;
                prml -= period;
            } else if prmf <= imf && prml <= imf {
                prmf += period;
                prml += period;
            }
        } else {
            if prmf > prml {
                prmf -= period;
            }
            if prmf >= iml && prml >= iml {
                prmf -= period;
                prml -= period;
            } else if prmf <= imf && prml <= imf {
                prmf += period;
                prml += period;
            }
        }
    }

    // `cxx:1919-1945`: shift the p-curve into the surface's U period.
    if s1.is_u_periodic() {
        if let Some(c2d) = &cimg2d {
            let speriod = s1.u_period();
            let pf = c2d.d0(prmf);
            let pl = c2d.d0(prml);
            let (uf, ul) = (pf.x(), pl.x());
            let mut ptra = 0.0;
            let mut ustart = uf.min(ul);
            while ustart < -Precision::PCONFUSION {
                ustart += speriod;
                ptra += speriod;
            }
            while ustart > speriod - Precision::PCONFUSION {
                ustart -= speriod;
                ptra -= speriod;
            }
            if ptra != 0.0 {
                let mut t = occt_core::gp::GpTrsf2d::identity();
                t.set_translation_vec(&GpVec2d::new(ptra, 0.0));
                cimg2d = Some(Arc::from(c2d.transformed(&t)));
            }
        }
    }

    // `cxx:1947-1958`.
    if prmf < prml {
        vf.0.set_orientation(Orientation::Forward);
        vl.0.set_orientation(Orientation::Reversed);
    } else {
        vf.0.set_orientation(Orientation::Reversed);
        vl.0.set_orientation(Orientation::Forward);
        rev = true;
    }

    // `cxx:1959-1968`.
    let b = TopoBuilder::new();
    let registry = GeometryRegistry::global();
    let first = cimg.first_parameter();
    let last = cimg.last_parameter();
    let mut new_edg = b.make_edge(cimg.clone(), first, last);
    registry.set_edge_tolerance(&new_edg.0, CONFUSION);
    b.add(&mut new_edg.0, &vf.0);
    b.add(&mut new_edg.0, &vl.0);
    registry.set_vertex_parameter_on_edge(&new_edg.0, &vf.0, prmf);
    registry.set_vertex_parameter_on_edge(&new_edg.0, &vl.0, prml);
    if app_s1 {
        if let Some(c2d) = &cimg2d {
            if let Some(s1) = BRepTool::face_surface(&Face(face.clone())) {
                registry.update_edge_pcurve_on_surface(&new_edg.0, &s1, Some(c2d.clone()), CONFUSION);
            }
        }
    }
    let _ = (GpLin::default, GpPnt::default);
    if rev {
        new_edg.0.set_orientation(Orientation::Reversed);
    }
    Some(new_edg)
}
