//! `HLRTopoBRep_DSFiller` (`ModelingAlgorithms/TKHLR/HLRTopoBRep/HLRTopoBRep_DSFiller.hxx:29-72`,
//! `.cxx:56-760`): turns the contour lines of every face into edges of
//! `HLRTopoBRep_Data`.
//!
//! `Insert` (`cxx:62-111`) walks the faces, runs `Contap_Contour::Perform`
//! (`cxx:96`) and stores the result through `InsertFace` (`cxx:120-592`); the
//! outline vertices it creates are then used by `ProcessEdges`
//! (`cxx:714-759`) to split the real edges at every outline vertex.
//!
//! Deviations:
//! * `TopoDS_Vertex` / `TopoDS_Edge` value handles are `Vertex` / `Edge`
//!   newtypes over `TopoShape`; an orientation write is
//!   [`TopoShape::set_orientation`] on the wrapped shape.
//! * `BRep_Builder::MakeEdge(E, C, Tol)` + `Range(E, First, Last)` is
//!   [`TopoBuilder::make_edge`] (which stores the range) plus
//!   `set_edge_tolerance`.
//! * `BRep_Builder::UpdateEdge(E, C2d, F, Tol)` is
//!   `update_edge_pcurve_on_surface` keyed by the face's stored surface, the
//!   same handle `BRep_Tool::Surface(F)` answers.
//! * `TopoDS_Edge::EmptyCopy()` (`BRep_TEdge.cxx:105-128`) is
//!   [`empty_copied_edge`]: a new `TShape` carrying the same curve
//!   representation and pcurve slots.
//! * `Standard_ProgramError` in the `Contap_Restriction` arm of the line
//!   switch (`cxx:510-512`) becomes `Err`.

use std::collections::HashSet;
use std::sync::Arc;

use occt_core::approx::ApproxParametrizationType;
use occt_core::bspl::banded_interp::knot_sequence;
use occt_core::precision::CONFUSION;
use occt_geom::bspline_curve::GeomBSplineCurve;
use occt_geom::geom_proj_lib;
use occt_geom::{Curve, GeomCircle, GeomLine};
use occt_geom2d::bspline_curve::Geom2dBSplineCurve;
use occt_geom2d::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::brep_approx::{ApproxLine, BRepApproxApprox};
use crate::brep_top_adaptor_tool::BrepTopAdaptorTool;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::contap::{ContapContour, ContapIType, ContapLine, ContapPoint};
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{edge_vertices, explorer_shapes, shape_map_key, ShapeMapKey};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::data::HlrTopoBRepData;
use super::face_iso_liner;
use super::out_liner::ShapeToTopAdaptorToolMap;

/// `RealLast()` (`Standard_Real.hxx`), the initial bound of the walking
/// approximation extents (`cxx:360-361`).
const REAL_LAST: f64 = f64::MAX;

/// `HLRTopoBRep_DSFiller::Insert(S, FO, DS, MST, nbIso)` (`cxx:62-111`).
///
/// `FO` must already carry its contour mode (`Contap_Contour::Init`); this
/// function only calls `Perform` per face, as OCCT does.
pub fn insert(
    s: &TopoShape,
    fo: &mut ContapContour,
    ds: &mut HlrTopoBRepData,
    mst: &mut ShapeToTopAdaptorToolMap,
    nb_iso: i32,
) -> Result<(), String> {
    // `NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> ShapeMap;`
    let mut shape_map: HashSet<ShapeMapKey> = HashSet::new();
    ds.clear(); // `DS.Clear()` (`cxx:71`)
    // `bool withPCurve = true; // instead of nbIso != 0;` (`cxx:72`).
    let with_pcurve = true;
    let mut f: i32 = 0;

    // `TopExp_Explorer ex(S, TopAbs_FACE);` (`cxx:70`).
    for face_shape in explorer_shapes(s, ShapeType::Face) {
        if !shape_map.insert(shape_map_key(&face_shape)) {
            continue; // `if (ShapeMap.Add(ex.Current()))` (`cxx:78`)
        }
        f += 1;

        // `TopoDS_Face S1 = TopoDS::Face(ex.Current()); S1.Orientation(FORWARD);`
        // (`cxx:81-82`).
        let mut s1 = face_shape;
        s1.set_orientation(Orientation::Forward);
        let face = Face::wrap(s1.clone())
            .ok_or_else(|| "HLRTopoBRep_DSFiller::Insert: not a face".to_string())?;

        // `if (MST.IsBound(S1)) { ... } else { BRepTopAdaptor_Tool BRT(S1,
        // Precision::PConfusion()); MST.Bind(S1, BRT); }` (`cxx:85-95`).
        let key = shape_map_key(&s1);
        if !mst.contains_key(&key) {
            mst.insert(key.clone(), BrepTopAdaptorTool::of_face(&face));
        }

        // `FO.Perform(Surface, Domain);` (`cxx:96`). `BRepTopAdaptor_Tool(F,
        // Tol2d)` builds the restricted `BRepAdaptor_Surface` (`cxx:30-41`).
        let surf = SurfaceAdaptor::of_face_restricted(&face);
        {
            let tool = mst
                .get_mut(&key)
                .ok_or_else(|| "HLRTopoBRep_DSFiller::Insert: unbound face".to_string())?;
            let domain = tool.topol_tool_mut();
            fo.perform(surf, domain)?;
        }

        // `if (FO.IsDone()) { if (!FO.IsEmpty()) InsertFace(...); }`
        // (`cxx:98-103`).
        if fo.is_done() && !fo.is_empty() {
            insert_face(&face, fo, ds, with_pcurve)?;
        }
        // `if (nbIso != 0) HLRTopoBRep_FaceIsoLiner::Perform(f, S1, DS, nbIso);`
        // (`cxx:105-108`).
        if nb_iso != 0 {
            face_iso_liner::perform(f, &face, ds, nb_iso);
        }
    }

    // `ProcessEdges(DS);` (`cxx:111`).
    process_edges(ds);
    Ok(())
}

/// `TopoDS_Edge::EmptyCopy()` (`BRep_TEdge.cxx:105-128`): a new `TShape` with
/// the same orientation, location, curve representation, tolerance and pcurves,
/// and no sub-shapes.
fn empty_copied_edge(e: &Edge) -> Edge {
    let copy = Edge::wrap(e.shape().empty_copied())
        .expect("HLRTopoBRep_DSFiller: EmptyCopy of an edge is an edge");
    let reg = GeometryRegistry::global();
    if let Some(g) = reg.edge_geom(e.shape()) {
        reg.set_edge(copy.shape(), g);
    }
    reg.copy_edge_pcurve_slots(copy.shape(), e.shape());
    copy
}

/// `HLRTopoBRep_DSFiller::InsertFace(FI, F, FO, DS, withPCurve)` (`cxx:120-592`).
fn insert_face(
    f: &Face,
    fo: &ContapContour,
    ds: &mut HlrTopoBRepData,
    with_pcurve: bool,
) -> Result<(), String> {
    // `const double tol = BRep_Tool::Tolerance(F);` (`cxx:128`).
    let tol = BRepTool::face_tolerance(f);
    // `DS.AddIntL(F); DS.AddOutL(F);` (`cxx:129-130`).
    ds.add_int_l(f);
    ds.add_out_l(f);

    // `NCollection_List<TopoDS_Shape>& IntL = DS.AddIntL(F);` (`cxx:129`):
    // the pieces appended below, mirrored locally for the correction pass.
    let mut int_l: Vec<TopoShape> = Vec::new();

    // `for (CurLine = 1; CurLine <= NbLines; CurLine++)` (`cxx:148-149`).
    let nb_lines = fo.nb_lines();
    for cur_line in 1..=nb_lines {
        let Some(line) = fo.line(cur_line) else {
            continue;
        };
        let nb_points = line.nb_vertex();
        if line.type_contour() == ContapIType::Restriction {
            // `TopoDS_Edge E = Line.Arc()->Edge(); OutL.Append(E);
            // TopExp::Vertices(E, VF, VL);` (`cxx:157-159`).
            let e = line.arc()?.edge().clone();
            ds.add_out_l(f).push(e.shape().clone());
            let (vf, vl) = edge_vertices(&e);

            // `for (CurPoint = 1; CurPoint <= NbPoints; CurPoint++)`
            // (`cxx:161`).
            for cur_point in 1..=nb_points {
                let Some(p) = line.vertex(cur_point) else {
                    continue;
                };
                if !p.is_internal() {
                    continue;
                }
                // `P.Value().IsEqual(Pnt(VF), Tol(VF)) &&
                // P.Value().IsEqual(Pnt(VL), Tol(VL))` (`cxx:167-175`).
                let (Some(vf), Some(vl)) = (&vf, &vl) else {
                    continue;
                };
                let at_vf = p
                    .value()
                    .is_equal_with_tolerance(&BRepTool::vertex_point(vf), BRepTool::vertex_tolerance(vf));
                let at_vl = p
                    .value()
                    .is_equal_with_tolerance(&BRepTool::vertex_point(vl), BRepTool::vertex_tolerance(vl));
                if at_vf && at_vl {
                    insert_vertex(p, tol, &e, ds);
                }
            }
            continue;
        }

        // `else { for (CurPoint = 1; CurPoint <= NbPoints; CurPoint++) }`
        // (`cxx:177-179`).
        let mut vf: Option<Vertex> = None;
        for cur_point in 1..=nb_points {
            let Some(pf) = line.vertex(cur_point) else {
                continue;
            };
            // `if (PF.IsInternal() && CurPoint != 1) VF = VL; else
            // VF = MakeVertex(PF, tol, DS);` (`cxx:182-190`).
            if pf.is_internal() && cur_point != 1 {
                // `VF = VL`: the `VL` of the previous iteration already points
                // at this vertex, so nothing to do.
            } else {
                vf = Some(make_vertex(pf, tol, ds)?);
            }
            let parf = pf.parameter_on_line();

            if cur_point >= nb_points {
                continue;
            }
            let Some(pl) = line.vertex(cur_point + 1) else {
                continue;
            };
            let vl = make_vertex(pl, tol, ds)?;
            let parl = pl.parameter_on_line();

            // `if ((parL - parF) > Precision::PConfusion())` (`cxx:199`).
            if !(parl - parf > CONFUSION) {
                continue;
            }

            let mut first = parf;
            let mut last = parl;
            let mut insuffisant_number_of_points = false;
            let mut c: Option<Arc<dyn Curve>> = None;
            let mut c2d: Option<Arc<dyn Curve2d>> = None;

            // `switch (Line.TypeContour())` (`cxx:208-517`).
            match line.type_contour() {
                ContapIType::Lin => {
                    // `C = new Geom_Line(Line.Line()); C2d =
                    // GeomProjLib::Curve2d(C, first, last, BRep_Tool::Surface(F), 1e-7);`
                    // (`cxx:211-219`).
                    c = Some(Arc::new(GeomLine::new(line.line()?)));
                    if with_pcurve {
                        let s = BRepTool::face_surface(f)
                            .ok_or_else(|| "HLRTopoBRep_DSFiller: face without surface".to_string())?;
                        let mut t = 1e-7;
                        c2d = geom_proj_lib::curve_2d_range(
                            c.as_ref().unwrap(),
                            first,
                            last,
                            &s,
                            &mut t,
                        );
                    }
                }
                ContapIType::Circle => {
                    // `C = new Geom_Circle(Line.Circle());` (`cxx:223`), then
                    // the located surface (`cxx:226-234`).
                    c = Some(Arc::new(GeomCircle::new(line.circle()?)));
                    if with_pcurve {
                        let s = located_face_surface(f)?;
                        let mut t = 1e-7;
                        c2d = geom_proj_lib::curve_2d_range(
                            c.as_ref().unwrap(),
                            first,
                            last,
                            &s,
                            &mut t,
                        );
                    }
                }
                ContapIType::Walking => {
                    // `int ipF = int(parF); int ipL = int(parL);` (`cxx:240-241`).
                    let ipf = parf as i32;
                    let ipl = parl as i32;
                    if ipl - ipf < 1 {
                        // `InsuffisantNumberOfPoints = true;` (`cxx:243-246`).
                        insuffisant_number_of_points = true;
                    } else if ipl - ipf < 5 {
                        // `cxx:321-350`: degree-1 interpolation of `nbp` points.
                        let nbp = (ipl - ipf + 1) as usize;
                        let mut knots = Vec::with_capacity(nbp);
                        let mut mults = vec![1i32; nbp];
                        let mut points = Vec::with_capacity(nbp);
                        for i in 1..=nbp {
                            knots.push(i as f64);
                            points.push(walking_point(line, i as i32 + ipf - 1)?.value());
                        }
                        mults[0] = 2;
                        mults[nbp - 1] = 2;
                        let curve = GeomBSplineCurve::from_poles_knots_mults(
                            points,
                            knots.clone(),
                            mults.clone(),
                            1,
                        )
                        .map_err(str::to_string)?;
                        c = Some(Arc::new(curve));

                        if with_pcurve {
                            let (xs, ys) = walking_uvs(line, nbp, ipf)?;
                            let flat = knot_sequence(&knots, &mults, 1);
                            let curve2d = Geom2dBSplineCurve::from_flat(xs, ys, None, flat, 1, false)
                                .map_err(str::to_string)?;
                            c2d = Some(Arc::new(curve2d));
                        }
                        first = 1.0;
                        last = nbp as f64;
                    } else {
                        // `cxx:351-508`: approximate `nbp` walked points.
                        let nbp = (ipl - ipf + 1) as usize;
                        let mut knots = Vec::with_capacity(nbp);
                        let mut mults = vec![1i32; nbp];
                        let mut points = Vec::with_capacity(nbp);
                        let (mut minx, mut miny, mut minz) = (REAL_LAST, REAL_LAST, REAL_LAST);
                        let (mut maxx, mut maxy, mut maxz) = (-REAL_LAST, -REAL_LAST, -REAL_LAST);
                        for i in 1..=nbp {
                            knots.push(i as f64);
                            let p = walking_point(line, i as i32 + ipf - 1)?.value();
                            minx = minx.min(p.x());
                            miny = miny.min(p.y());
                            minz = minz.min(p.z());
                            maxx = maxx.max(p.x());
                            maxy = maxy.max(p.y());
                            maxz = maxz.max(p.z());
                            points.push(p);
                        }
                        mults[0] = 2;
                        mults[nbp - 1] = 2;
                        let app_c = GeomBSplineCurve::from_poles_knots_mults(
                            points,
                            knots.clone(),
                            mults.clone(),
                            1,
                        )
                        .map_err(str::to_string)?;

                        let mut app_c2d: Option<Geom2dBSplineCurve> = None;
                        let (mut minu, mut minv) = (REAL_LAST, REAL_LAST);
                        let (mut maxu, mut maxv) = (-REAL_LAST, -REAL_LAST);
                        if with_pcurve {
                            let (xs, ys) = walking_uvs(line, nbp, ipf)?;
                            for i in 0..nbp {
                                minu = minu.min(xs[i]);
                                minv = minv.min(ys[i]);
                                maxu = maxu.max(xs[i]);
                                maxv = maxv.max(ys[i]);
                            }
                            let flat = knot_sequence(&knots, &mults, 1);
                            app_c2d = Some(
                                Geom2dBSplineCurve::from_flat(xs, ys, None, flat, 1, false)
                                    .map_err(str::to_string)?,
                            );
                        }
                        // `cxx:361-362`: `first = 1; last = nbp;` - OCCT sets
                        // them here and again inside the `!Approx.IsDone()`
                        // arm below.
                        #[allow(unused_assignments)]
                        {
                            first = 1.0;
                            last = nbp as f64;
                        }

                        // `BRepApprox_ApproxLine AppLine(AppC, AppC2d, CNull)`
                        // (`cxx:365-368`).
                        let app_line = ApproxLine::new(
                            Some(Arc::new(app_c.clone())),
                            app_c2d.clone().map(|c| Arc::new(c) as Arc<dyn Curve2d>),
                            None,
                        );

                        // `cxx:433-479`: the extents and the two tolerances.
                        let (mut maxx, maxy, maxz) =
                            (maxx - minx, maxy - miny, maxz - minz);
                        let (mut maxu, maxv) = (maxu - minu, maxv - minv);
                        if maxy > maxx {
                            maxx = maxy;
                        }
                        if maxz > maxx {
                            // Faithful to `cxx:449-451`, which assigns `Maxy`.
                            maxx = maxy;
                        }
                        if maxv > maxu {
                            maxu = maxv;
                        }
                        const TOL: f64 = 0.0001;
                        let tol3d = (TOL * maxx).clamp(1e-12, 0.1);
                        let tol2d = (TOL * maxu).clamp(1e-12, 0.1);

                        let mut approx = BRepApproxApprox::new();
                        // `Approx.SetParameters(TOL3d, TOL2d, 4, 8, 0, 30, false)`
                        // (`cxx:478`; the defaults are `NbPntMax = 30`,
                        // `ApproxWithTangency = true` and the call passes
                        // `tg = false`, `Parametrization = Approx_ChordLength`).
                        approx.set_parameters(
                            tol3d,
                            tol2d,
                            4,
                            8,
                            0,
                            30,
                            false,
                            ApproxParametrizationType::ChordLength,
                        );
                        // `Approx.Perform(AppLine, true, true, false, 1, nbp)`
                        // (`cxx:479`).
                        approx.perform_line(&app_line, true, true, false, 1, nbp as i32);

                        if !approx.is_done() {
                            // `C = AppC; C2d = AppC2d; first = 1; last = nbp;`
                            // (`cxx:481-486`).
                            c = Some(Arc::new(app_c));
                            c2d = app_c2d.map(|c| Arc::new(c) as Arc<dyn Curve2d>);
                            first = 1.0;
                            last = nbp as f64;
                        } else {
                            // `cxx:488-506`.
                            let app_val = approx.value();
                            let poles3d = app_val.base().curve3d(1);
                            let c3d = GeomBSplineCurve::from_poles_knots_mults(
                                poles3d,
                                app_val.knots().to_vec(),
                                app_val.multiplicities().to_vec(),
                                app_val.degree() as usize,
                            )
                            .map_err(str::to_string)?;
                            c = Some(Arc::new(c3d));

                            let poles2d = app_val.base().curve2d(1);
                            let xs: Vec<f64> = poles2d.iter().map(|p| p.x()).collect();
                            let ys: Vec<f64> = poles2d.iter().map(|p| p.y()).collect();
                            let flat = knot_sequence(
                                app_val.knots(),
                                app_val.multiplicities(),
                                app_val.degree(),
                            );
                            let c2d_curve = Geom2dBSplineCurve::from_flat(
                                xs,
                                ys,
                                None,
                                flat,
                                app_val.degree() as usize,
                                false,
                            )
                            .map_err(str::to_string)?;
                            first = c2d_curve.first_parameter();
                            last = c2d_curve.last_parameter();
                            c2d = Some(Arc::new(c2d_curve));
                        }
                    }
                }
                ContapIType::Restriction => {
                    // `throw Standard_ProgramError("... Restriction")`
                    // (`cxx:510-512`).
                    return Err(
                        "Standard_ProgramError: HLRTopoBRep_DSFiller::InsertFace : Restriction"
                            .to_string(),
                    );
                }
            }

            // `if (!InsuffisantNumberOfPoints) { ... }` (`cxx:519-541`).
            if insuffisant_number_of_points {
                continue;
            }
            let Some(curve) = c else {
                continue;
            };
            let Some(mut vf_v) = vf.clone() else {
                continue;
            };
            let mut vl_v = vl;

            // `B.MakeEdge(E, C, tol);` (`cxx:522-524`).
            let maker = TopoBuilder;
            let mut e = maker.make_edge(curve, first, last);
            GeometryRegistry::global().set_edge_tolerance(&e.0, tol);
            // `VF.Orientation(FORWARD); VL.Orientation(REVERSED);` (`cxx:525-526`).
            vf_v.0.set_orientation(Orientation::Forward);
            vl_v.0.set_orientation(Orientation::Reversed);
            // `B.Add(E, VF); B.Add(E, VL);` (`cxx:527-528`).
            maker.add(&mut e.0, &vf_v.0);
            maker.add(&mut e.0, &vl_v.0);

            // `if (!C2d.IsNull()) B.UpdateEdge(E, C2d, F, Tolerance(F));`
            // (`cxx:530-533`). `BRep_Builder::UpdateEdge(E, C, F, Tol)` binds
            // the pcurve to `BRep_Tool::Surface(F)` (`BRep_Builder.cxx:987-993`).
            if let Some(c2d_curve) = c2d {
                if let Some(surface) = BRepTool::face_surface(f) {
                    GeometryRegistry::global()
                        .update_edge_pcurve_on_surface(&e.0, &surface, Some(c2d_curve), tol);
                }
            }

            // `IntL.Append(E);` (`cxx:538`).
            int_l.push(e.0.clone());
        }
    }

    // `Correction of internal outlines: unite coinciding vertices`
    // (`cxx:547-591`).
    let sq_tol = tol * tol;
    let nb = int_l.len();
    for i1 in 0..nb {
        let (Some(v1_0), Some(v1_1)) = int_l_vertices(&int_l[i1]) else {
            continue;
        };
        for i2 in i1..nb {
            let line2 = int_l[i2].clone();
            if crate::topo_tools::is_same(&line2, &int_l[i1]) {
                continue; // (`cxx:561-564`)
            }
            let (Some(v2_0), Some(v2_1)) = int_l_vertices(&line2) else {
                continue;
            };
            let a_ver = [v1_0.clone(), v1_1.clone()];
            let a_ver2 = [v2_0, v2_1];
            for i in 0..2 {
                if i == 1 && crate::topo_tools::is_same(a_ver[0].shape(), a_ver[1].shape()) {
                    continue; // (`cxx:569-572`)
                }
                let pnt1 = BRepTool::vertex_point(&a_ver[i]);
                let mut a_ver_i = a_ver[i].clone();
                for j in 0..2 {
                    if crate::topo_tools::is_same(a_ver_i.shape(), a_ver2[j].shape()) {
                        continue; // (`cxx:578-581`)
                    }
                    let pnt2 = BRepTool::vertex_point(&a_ver2[j]);
                    if pnt1.square_distance(&pnt2) > sq_tol {
                        continue;
                    }
                    // `aBB.Remove(anIntLine2, aVer2[j]); aVer[i].Orientation(
                    // (j == 0) ? FORWARD : REVERSED); aBB.Add(anIntLine2, aVer[i]);`
                    // (`cxx:584-588`). The orientation lands on the local
                    // `aVer[i]`, so the next `j` starts from it.
                    let mut line2_mut = line2.clone();
                    TopoBuilder.remove(&mut line2_mut, a_ver2[j].shape());
                    a_ver_i.0.set_orientation(if j == 0 {
                        Orientation::Forward
                    } else {
                        Orientation::Reversed
                    });
                    TopoBuilder.add(&mut line2_mut, &a_ver_i.0);
                }
            }
        }
    }

    // Publish the pieces in `DS.AddIntL(F)` (the reference `IntL` of
    // `cxx:129`).
    ds.add_int_l(f).extend(int_l);
    Ok(())
}

/// `TopExp::Vertices(E, V1, V2)` for an edge whose children are the two
/// outline vertices added by `InsertFace`.
fn int_l_vertices(e: &TopoShape) -> (Option<Vertex>, Option<Vertex>) {
    match Edge::wrap(e.clone()) {
        Some(edge) => edge_vertices(&edge),
        None => (None, None),
    }
}

/// `BRep_Tool::Surface(F, Loc)` with `Loc` applied when it is not the identity
/// (`cxx:226-232`).
fn located_face_surface(f: &Face) -> Result<Arc<dyn occt_geom::Surface>, String> {
    if f.shape().location().is_identity() {
        BRepTool::face_surface(f)
            .ok_or_else(|| "HLRTopoBRep_DSFiller: face without surface".to_string())
    } else {
        BRepTool::face_surface_world(f)
            .ok_or_else(|| "HLRTopoBRep_DSFiller: face without surface".to_string())
    }
}

/// `Line.Point(Index)` of a `Contap_Walking` line (`cxx:331`, `:345`, `:368`).
fn walking_point<'a>(line: &'a ContapLine, index: i32) -> Result<&'a crate::int_tools_wline::PntOn2S, String> {
    if index < 1 {
        return Err("HLRTopoBRep_DSFiller: walking point index below 1".to_string());
    }
    line.point(index as usize).map_err(str::to_string)
}

/// `Line.Point(i + ipF - 1).ParametersOnS2(u, v)` for `i = 1..=nbp`
/// (`cxx:340-347`, `:401-424`).
fn walking_uvs(line: &ContapLine, nbp: usize, ipf: i32) -> Result<(Vec<f64>, Vec<f64>), String> {
    let mut xs = Vec::with_capacity(nbp);
    let mut ys = Vec::with_capacity(nbp);
    for i in 1..=nbp {
        let p = walking_point(line, i as i32 + ipf - 1)?;
        let (u, v) = p.parameters_on_s2();
        xs.push(u);
        ys.push(v);
    }
    Ok((xs, ys))
}

/// `HLRTopoBRep_DSFiller::MakeVertex(P, tol, DS)` (`cxx:598-659`).
fn make_vertex(p: &ContapPoint, tol: f64, ds: &mut HlrTopoBRepData) -> Result<Vertex, String> {
    let maker = TopoBuilder;
    if p.is_vertex() {
        // `V = P.Vertex()->Vertex(); DS.AddOutV(V);` (`cxx:603-608`).
        let v = p.vertex().map_err(str::to_string)?.vertex().clone();
        ds.add_out_v(&v);
        return Ok(v);
    }

    if p.is_on_arc() {
        // `cxx:612-643`: insert the vertex in the arc's vertex list.
        let e = p.arc().map_err(str::to_string)?.edge().clone();
        let par = p.parameter_on_arc().map_err(str::to_string)?;
        let p3d = *p.value();

        let mut found: Option<Vertex> = None;
        let mut inserted = false;
        ds.init_vertex(&e);
        while ds.more_vertex() {
            let cur_v = ds
                .vertex()
                .cloned()
                .ok_or_else(|| "HLRTopoBRep_DSFiller::MakeVertex: missing vertex".to_string())?;
            let cur_p = ds.parameter();
            let cur_v = Vertex::wrap(cur_v)
                .ok_or_else(|| "HLRTopoBRep_DSFiller::MakeVertex: not a vertex".to_string())?;
            let ppp = BRepTool::vertex_point(&cur_v);
            let ttt = BRepTool::vertex_tolerance(&cur_v);
            if p3d.is_equal_with_tolerance(&ppp, ttt) {
                found = Some(cur_v);
                break;
            } else if par < cur_p {
                let v = maker.make_vertex(p3d, tol);
                ds.insert_before(&v.0, par);
                found = Some(v);
                inserted = true;
                break;
            }
            ds.next_vertex();
        }
        if !inserted && !ds.more_vertex() {
            // `if (!DS.MoreVertex()) { B.MakeVertex(V, P.Value(), tol);
            // DS.Append(V, Par); }` (`cxx:637-641`).
            let v = maker.make_vertex(p3d, tol);
            ds.append(&v.0, par);
            found = Some(v);
        }
        let v = found
            .ok_or_else(|| "HLRTopoBRep_DSFiller::MakeVertex: vertex not resolved".to_string())?;
        // `DS.AddOutV(V);` (`cxx:642`).
        ds.add_out_v(&v);
        return Ok(v);
    }

    // `cxx:645-655`: an internal point gets a fresh vertex.
    let v = maker.make_vertex(*p.value(), tol);
    if p.is_internal() {
        ds.add_int_v(&v);
    } else {
        ds.add_out_v(&v);
    }
    Ok(v)
}

/// `HLRTopoBRep_DSFiller::InsertVertex(P, tol, E, DS)` (`cxx:668-706`).
fn insert_vertex(p: &ContapPoint, tol: f64, e: &Edge, ds: &mut HlrTopoBRepData) {
    let maker = TopoBuilder;
    let mut v: Option<Vertex> = None;
    if p.is_vertex() {
        // `V = P.Vertex()->Vertex();` (`cxx:676-678`).
        v = p.vertex().ok().map(|hv| hv.vertex().clone());
    } else {
        // `double Par = P.ParameterOnLine();` (`cxx:681`).
        let par = p.parameter_on_line();
        let mut inserted = false;
        ds.init_vertex(e);
        while ds.more_vertex() {
            let Some(cur_v) = ds.vertex().cloned() else {
                break;
            };
            let cur_p = ds.parameter();
            let Some(cur_v) = Vertex::wrap(cur_v) else {
                break;
            };
            if p.value()
                .is_equal_with_tolerance(&BRepTool::vertex_point(&cur_v), BRepTool::vertex_tolerance(&cur_v))
            {
                v = Some(cur_v);
                break;
            } else if par < cur_p {
                let nv = maker.make_vertex(*p.value(), tol);
                ds.insert_before(&nv.0, par);
                v = Some(nv);
                inserted = true;
                break;
            }
            ds.next_vertex();
        }
        if !inserted && !ds.more_vertex() {
            // `if (!DS.MoreVertex()) { B.MakeVertex(V, P.Value(), tol);
            // DS.Append(V, Par); }` (`cxx:698-703`).
            let nv = maker.make_vertex(*p.value(), tol);
            ds.append(&nv.0, par);
            v = Some(nv);
        }
    }
    // `DS.AddIntV(V);` (`cxx:705`).
    if let Some(v) = v {
        ds.add_int_v(&v);
    }
}

/// `HLRTopoBRep_DSFiller::ProcessEdges(DS)` (`cxx:714-759`): split every edge
/// that received outline vertices.
fn process_edges(ds: &mut HlrTopoBRepData) {
    let maker = TopoBuilder;

    // `for (DS.InitEdge(); DS.MoreEdge(); DS.NextEdge())` (`cxx:720`).
    ds.init_edge();
    while ds.more_edge() {
        let Some(e_shape) = ds.edge().cloned() else {
            ds.next_edge();
            continue;
        };
        ds.next_edge();
        let Some(e) = Edge::wrap(e_shape) else {
            continue;
        };

        // `NCollection_List<TopoDS_Shape>& SplE = DS.AddSplE(E);` (`cxx:723`).
        ds.add_spl_e(&e);

        // `VF = TopExp::FirstVertex(E); VL = TopExp::LastVertex(E);
        // BRep_Tool::Range(E, PF, PL);` (`cxx:724-726`).
        let (mut vf, vl) = edge_vertices(&e);
        let (mut pf, pl) = BRepTool::edge_parameters(&e);
        // `VF.Orientation(FORWARD); VL.Orientation(REVERSED);` (`cxx:728-729`).
        if let Some(v) = vf.as_mut() {
            v.0.set_orientation(Orientation::Forward);
        }

        // `for (DS.InitVertex(E); DS.MoreVertex(); DS.NextVertex())`
        // (`cxx:730`).
        let mut pieces: Vec<TopoShape> = Vec::new();
        let mut cursor = 0usize;
        ds.init_vertex(&e);
        while ds.more_vertex() {
            cursor += 1;
            let vi = ds.vertex().cloned();
            let pi = ds.parameter();
            ds.next_vertex();

            let Some(mut vi) = vi.and_then(Vertex::wrap) else {
                break;
            };
            let Some(vf_v) = vf.as_ref() else {
                break;
            };
            let Some(vl_v) = vl.as_ref() else {
                break;
            };
            // `VI.Orientation(REVERSED);` (`cxx:733`).
            vi.0.set_orientation(Orientation::Reversed);

            // `newE = E; newE.EmptyCopy(); newE.Orientation(FORWARD);
            // B.Add(newE, VF); ...` (`cxx:734-743`).
            let mut new_e = empty_copied_edge(&e);
            new_e.0.set_orientation(Orientation::Forward);
            add_vertex_with_parameter(&maker, &mut new_e, vf_v, pf, vf_v);
            add_vertex_with_parameter(&maker, &mut new_e, &vi, pi, &vi);

            // `newE.Orientation(E.Orientation()); SplE.Append(newE);`
            // (`cxx:744-746`).
            new_e.0.set_orientation(e.shape().orientation());
            pieces.push(new_e.0.clone());

            // `VF = VI; PF = PI; VF.Orientation(FORWARD);` (`cxx:747-748`).
            vf = Some(vi.clone());
            pf = pi;
            if let Some(v) = vf.as_mut() {
                v.0.set_orientation(Orientation::Forward);
            }
            let _ = vl_v;
        }
        let _ = cursor;

        // `newE = E; newE.EmptyCopy(); newE.Orientation(FORWARD);
        // B.Add(newE, VF); ... SplE.Append(newE);` (`cxx:749-757`).
        let (Some(vf_v), Some(vl_v)) = (vf.as_ref(), vl.as_ref()) else {
            continue;
        };
        let mut new_e = empty_copied_edge(&e);
        new_e.0.set_orientation(Orientation::Forward);
        add_vertex_with_parameter(&maker, &mut new_e, vf_v, pf, vf_v);
        add_vertex_with_parameter(&maker, &mut new_e, vl_v, pl, vl_v);
        new_e.0.set_orientation(e.shape().orientation());
        pieces.push(new_e.0.clone());

        ds.add_spl_e(&e).extend(pieces);
    }
}

/// `B.Add(newE, V); B.UpdateVertex(V, P, newE, BRep_Tool::Tolerance(V));`
/// (`cxx:738-742`, `:753-756`).
fn add_vertex_with_parameter(
    maker: &TopoBuilder,
    new_e: &mut Edge,
    vertex: &Vertex,
    parameter: f64,
    tolerance_source: &Vertex,
) {
    maker.add(&mut new_e.0, &vertex.0);
    // `BRep_Builder::UpdateVertex(V, P, E, Tol)` raises the edge tolerance
    // (`BRep_TEdge::UpdateTolerance`, `BRep_TEdge.lxx:33-37`) and stores the
    // vertex parameter.
    let tol = BRepTool::vertex_tolerance(tolerance_source);
    let reg = GeometryRegistry::global();
    if reg.edge_tolerance(new_e.shape()) < tol {
        reg.set_edge_tolerance(new_e.shape(), tol);
    }
    reg.set_vertex_parameter_on_edge(new_e.shape(), vertex.shape(), parameter);
}
