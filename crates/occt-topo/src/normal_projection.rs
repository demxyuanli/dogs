//! `BRepAlgo_NormalProjection` - normal projection of edges onto faces.
//!
//! Source: `ModelingAlgorithms/TKBool/BRepAlgo/BRepAlgo_NormalProjection.cxx`
//! (OCCT 8.0.0). `Build` is `cxx:163-585`, `IsElementary` `cxx:623-637`,
//! `BuildWire` `cxx:642-677`, the accessors `cxx:588-620`, the setters
//! `cxx:101-158`.
//!
//! Everything `Build` does is ported: the edge / face exploration
//! (`cxx:202-212`), the `ProjLib_HCompProjectedCurve` of each (edge, face) pair
//! (`cxx:238-239`), the `IsSinglePnt` arm (`cxx:257-269`), the two iso arms
//! (`cxx:275-311`), the `Only2d` test (`cxx:317-320`), the elementary shortcut
//! (`cxx:324-329`), the `Approx_CurveOnSurface` call with its `1.e3 * myTol3d`
//! reject (`cxx:335-341`), the degeneracy scan (`cxx:366-425`), the
//! degenerated-edge construction (`cxx:429-437`), the `UpdateEdge` /
//! vertex-tolerance writes (`cxx:443-451`) and the face-bounds block
//! (`cxx:454-535`): the non-degenerated arm runs `BRepAlgoAPI_Section(Face, prj)`
//! (`cxx:465-478`) through [`crate::brep_algo_api_section::SectionAlgo`], and
//! the two classification arms (`cxx:480-498` section-failure fallback,
//! `cxx:501-527` degenerated) go through [`NormalProjection::push_if_inside`].
//!
//! UNPORTED: `BOPAlgo_PaveFiller_6.cxx:506` - `IntTools_FaceFace::SetParameters`
//! is never called by the port's F/F stage, so the
//! `BOPAlgo_SectionAttribute(false, false, false)` that
//! `BRepAlgoAPI_Section::SetAttributes` (`BRepAlgoAPI_Section.cxx:209-213`)
//! installs has no effect: the section run always keeps the intersector's own
//! pcurves (the *builder* default). Same boundary as every other ported
//! `BOPAlgo_*` operation.
//!
//! `myFaceBounds` starts `true` (`cxx:75-78`), so a caller that does not want
//! the trimming passes `set_limit(false)` - which is what
//! `ChFi3d_Builder::ComputeCurve2d` does (`ChFi3d_Builder_C1.cxx:470`);
//! `CurveHermite::OrtProj` (`ChFi3d_Builder_CnCrn.cxx:406-410`) leaves it on.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::kernel::geomabs::Shape as GeomAbsShape;
use occt_core::precision::CONFUSION;
use occt_geom::approx::curve_on_surface::ApproxCurveOnSurface;
use occt_geom::approx_same_parameter::{u_resolution, v_resolution};
use occt_geom::projlib::{CompProjectedCurve, ProjectedCurve2d};
use occt_geom::trimmed::GeomTrimmedCurve;
use occt_geom::Curve;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;
use occt_geom2d::Geom2dBSplineCurve;

use crate::abs::ShapeType;
use crate::brep_algo_api_section::SectionAlgo;
use crate::brep_lib_make_edge;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::{FaceState, FClass2d};
use crate::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use crate::shape::{Edge, Face, TopoShape, Vertex, Wire};
use crate::topo_tools::{edges_of, explorer_shapes, faces_of, ShapeMapKey};
use crate::topopebrep::build::tools::write_edge_pcurve;

/// `BuildWire` (`cxx:642-677`) needs `myRes` to be a set of edges only.

/// `BRepLib_MakeEdge(C, U1, U2)` / `BRepLib_MakeEdge(C)`: an edge over the
/// parametric range of `curve` plus its two endpoint vertices
/// (`BRepLib_MakeEdge.cxx:606-796`, `:783-791`).
fn make_edge_over(builder: &TopoBuilder, curve: &Arc<dyn Curve>, first: f64, last: f64) -> Edge {
    let mut e = builder.make_edge(curve.clone(), first, last);
    if let Some(build) = brep_lib_make_edge::init_3d(curve, first, last) {
        brep_lib_make_edge::add_endpoint_vertices(builder, &mut e, &build);
    }
    e
}

/// `BRep_Builder::UpdateVertex(V, Tol)` (`BRep_Builder.cxx:1113-1125`): the
/// stored tolerance is the larger of `Tol` and the current one.
fn update_vertex_tolerance(v: &Vertex, tol: f64) {
    if tol > BRepTool::vertex_tolerance(v) {
        v.set_tolerance(tol);
    }
}

/// `BB.UpdateVertex(TopExp::FirstVertex(E), Tol)` + `LastVertex`
/// (`BRepAlgo_NormalProjection.cxx:327-328`, `:447-448`).
fn update_edge_vertex_tolerances(e: &Edge, tol: f64) {
    let (v1, v2) = crate::topo_tools::edge_vertices(e);
    if let Some(v) = v1 {
        update_vertex_tolerance(&v, tol);
    }
    if let Some(v) = v2 {
        update_vertex_tolerance(&v, tol);
    }
}

/// `BRepAlgo_NormalProjection` (`hxx:33-206`).
pub struct NormalProjection {
    /// `myIsDone` (`cxx:588-591`).
    is_done: bool,
    /// `myShape` - the shape projected *onto* (`Init`, `cxx:101-104`).
    shape: TopoShape,
    /// `myToProj` - the compound the projected shapes are added to.
    to_proj: TopoShape,
    /// `myRes` - the projection compound (`cxx:595-598`).
    res: TopoShape,
    /// `myAncestorMap` (`cxx:602-605`).
    ancestor_map: HashMap<ShapeMapKey, TopoShape>,
    /// `myCorresp` (`cxx:607-611`).
    corresp: HashMap<ShapeMapKey, TopoShape>,
    /// `myDescendants` (`cxx:615-619`).
    descendants: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `myTol3d` / `myTol2d` (`cxx:116-124`).
    tol3d: f64,
    tol2d: f64,
    /// `myContinuity` / `myMaxDegree` / `myMaxSeg` (`cxx:116-124`).
    continuity: u8,
    max_degree: i32,
    max_seg: i32,
    /// `myMaxDist` (`cxx:149-151`).
    max_dist: f64,
    /// `myWith3d` (`cxx:156-158`).
    with3d: bool,
    /// `myFaceBounds` (`cxx:142-145`).
    face_bounds: bool,
}

impl NormalProjection {
    /// `BRepAlgo_NormalProjection()` (`cxx:75-83`).
    pub fn new() -> Self {
        let mut this = Self {
            is_done: false,
            shape: TopoShape::new(ShapeType::Compound),
            to_proj: TopoShape::new(ShapeType::Compound),
            res: TopoShape::new(ShapeType::Compound),
            ancestor_map: HashMap::new(),
            corresp: HashMap::new(),
            descendants: HashMap::new(),
            tol3d: 0.0,
            tol2d: 0.0,
            continuity: GeomAbsShape::C1 as u8,
            max_degree: 0,
            max_seg: 0,
            max_dist: -1.0,
            with3d: true,
            face_bounds: true,
        };
        this.set_default_params();
        this
    }

    /// `BRepAlgo_NormalProjection(S)` (`cxx:86-97`).
    pub fn with_shape(shape: TopoShape) -> Self {
        let mut this = Self::new();
        this.init(shape);
        this
    }

    /// `Init(S)` (`cxx:101-104`).
    pub fn init(&mut self, shape: TopoShape) {
        self.shape = shape;
    }

    /// `Add(ToProj)` (`cxx:108-112`): `BB.Add(myToProj, ToProj)`.
    pub fn add(&mut self, to_proj: &TopoShape) {
        TopoBuilder::new().add(&mut self.to_proj, to_proj);
    }

    /// `SetParams` (`cxx:116-125`).
    pub fn set_params(
        &mut self,
        tol3d: f64,
        tol2d: f64,
        continuity: u8,
        max_degree: i32,
        max_seg: i32,
    ) {
        self.tol3d = tol3d;
        self.tol2d = tol2d;
        self.continuity = continuity;
        self.max_degree = max_degree;
        self.max_seg = max_seg;
    }

    /// `SetDefaultParams` (`cxx:131-138`).
    pub fn set_default_params(&mut self) {
        self.tol3d = 1e-4;
        self.tol2d = self.tol3d.powf(2.0 / 3.0);
        self.continuity = GeomAbsShape::C1 as u8;
        self.max_degree = 14;
        self.max_seg = 16;
    }

    /// `SetLimit(FaceBounds)` (`cxx:142-145`).
    pub fn set_limit(&mut self, face_bounds: bool) {
        self.face_bounds = face_bounds;
    }

    /// `SetMaxDistance(MaxDist)` (`cxx:149-151`).
    pub fn set_max_distance(&mut self, max_dist: f64) {
        self.max_dist = max_dist;
    }

    /// `Compute3d(With3d)` (`cxx:156-158`).
    pub fn compute3d(&mut self, with3d: bool) {
        self.with3d = with3d;
    }

    /// `IsDone()` (`cxx:588-591`).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// `Projection()` (`cxx:595-598`).
    pub fn projection(&self) -> &TopoShape {
        &self.res
    }

    /// `Ancestor(E)` (`cxx:602-605`).
    pub fn ancestor(&self, e: &Edge) -> Option<&TopoShape> {
        self.ancestor_map.get(&crate::bop_occt_util::shape_key(e.shape()))
    }

    /// `Couple(E)` (`cxx:607-611`).
    pub fn couple(&self, e: &Edge) -> Option<&TopoShape> {
        self.corresp.get(&crate::bop_occt_util::shape_key(e.shape()))
    }

    /// `Generated(S)` (`cxx:615-619`).
    pub fn generated(&self, s: &TopoShape) -> Option<&Vec<TopoShape>> {
        self.descendants.get(&crate::bop_occt_util::shape_key(s))
    }

    /// `IsElementary(C)` (`cxx:623-637`): a line, circle, ellipse, hyperbola
    /// or parabola.
    pub fn is_elementary(&self, c: &dyn Curve) -> bool {
        matches!(
            curve_kind(c),
            GeomAbsCurveType::Line
                | GeomAbsCurveType::Circle
                | GeomAbsCurveType::Ellipse
                | GeomAbsCurveType::Hyperbola
                | GeomAbsCurveType::Parabola
        )
    }

    /// `BRepTopAdaptor_FClass2d` midpoint classification shared by the two
    /// non-section arms of the face-bounds block: `cxx:480-498` (the section
    /// run failed) and `cxx:501-527` (the projection is degenerated). Appends
    /// `prj` to `myRes` / `DescenList` and binds the relations when the
    /// parameter-space midpoint of the p-curve is IN or ON the face.
    fn push_if_inside(
        &mut self,
        builder: &TopoBuilder,
        face: &Face,
        prj: &Edge,
        pcur2d: &Option<Arc<dyn Curve2d>>,
        descen_list: &mut Vec<TopoShape>,
        e: &Edge,
    ) -> Result<(), String> {
        // `cxx:483-490` / `:506-513` read `PCur2d` unconditionally.
        let pc = pcur2d
            .clone()
            .ok_or_else(|| "BRepAlgo_NormalProjection: null p-curve".to_string())?;
        let classifier = FClass2d::new(face, CONFUSION)?;
        let pmil = (pc.first_parameter() + pc.last_parameter()) / 2.0;
        let state = classifier.perform(pc.d0(pmil));
        if matches!(state, FaceState::In | FaceState::On) {
            builder.add(&mut self.res, &prj.0);
            descen_list.push(prj.0.clone());
            self.ancestor_map
                .insert(crate::bop_occt_util::shape_key(&prj.0), e.0.clone());
            self.corresp
                .insert(crate::bop_occt_util::shape_key(&prj.0), face.0.clone());
        }
        Ok(())
    }

    /// `Build()` (`cxx:163-585`).
    pub fn build(&mut self) -> Result<(), String> {
        self.is_done = false;
        self.res = TopoShape::new(ShapeType::Compound);
        let mut vertex_res = TopoShape::new(ShapeType::Compound);
        let mut ya_vertex_res = false;
        self.ancestor_map.clear();
        self.corresp.clear();
        self.descendants.clear();

        let builder = TopoBuilder::new();
        let edges = edges_of(&self.to_proj);
        let faces = faces_of(&self.shape);

        for e in &edges {
            let (u1, u2) = BRepTool::edge_parameters(e);
            let Some(raw) = BRepTool::edge_curve_world(e) else {
                continue;
            };
            // `BRepAdaptor_Curve(Edge)` (`cxx:220`): the edge curve restricted
            // to the edge range, which is the parameter space the projector
            // and `Approx_CurveOnSurface` report.
            let curve3d: Arc<dyn Curve> = Arc::new(GeomTrimmedCurve::new(raw, u1, u2));
            let elementary = self.is_elementary(curve3d.as_ref());
            let mut descen_list: Vec<TopoShape> = Vec::new();

            for f in &faces {
                let Some(surf) = BRepTool::face_surface_world(f) else {
                    continue;
                };
                // `cxx:230`: `TolU = UResolution(myTol3d) / 20`, idem `TolV`.
                let tol_u = u_resolution(surf.as_ref(), self.tol3d) / 20.0;
                let tol_v = v_resolution(surf.as_ref(), self.tol3d) / 20.0;
                // `cxx:238-239`: `ProjLib_HCompProjectedCurve(hsur, hcur, ...)`.
                let projector = CompProjectedCurve::with_tolerances_and_max_dist(
                    surf.clone(),
                    curve3d.clone(),
                    tol_u,
                    tol_v,
                    self.max_dist,
                );

                for k in 1..=projector.nb_curves() {
                    let Some((p2d, is_single)) = projector.is_single_pnt(k) else {
                        continue;
                    };
                    if is_single {
                        // `cxx:257-269`.
                        let p = surf.d0(p2d.x(), p2d.y());
                        let v = builder.make_vertex(p, CONFUSION);
                        let shape = v.0.clone();
                        descen_list.push(shape.clone());
                        builder.add(&mut vertex_res, &shape);
                        ya_vertex_res = true;
                        self.ancestor_map
                            .insert(crate::bop_occt_util::shape_key(&shape), e.0.clone());
                        continue;
                    }

                    // `cxx:273`: `Only2d = Only3d = false`.
                    let mut only2d = false;
                    let mut only3d = false;
                    let Some((udeb, ufin)) = projector.bounds(k) else {
                        continue;
                    };
                    let mut pcur2d: Option<Arc<dyn Curve2d>> = None;
                    let h_pcur: Arc<dyn Curve2d>;
                    if projector.is_u_iso(k).is_some_and(|(_, iso)| iso) {
                        // `cxx:275-292`.
                        let pdeb = projector.d0(udeb)?;
                        let pfin = projector.d0(ufin)?;
                        let bs2d = Geom2dBSplineCurve::new(
                            vec![pdeb.x(), pfin.x()],
                            vec![pdeb.y(), pfin.y()],
                            vec![udeb, udeb, ufin, ufin],
                            1,
                        )
                        .map_err(|e| e.to_string())?;
                        let trimmed = Arc::new(Geom2dTrimmedCurve::new(Arc::new(bs2d), udeb, ufin));
                        pcur2d = Some(trimmed.clone());
                        h_pcur = trimmed;
                        only3d = true;
                    } else if projector.is_v_iso(k).is_some_and(|(_, iso)| iso) {
                        // `cxx:298-311`.
                        let pdeb = projector.d0(udeb)?;
                        let pfin = projector.d0(ufin)?;
                        let bs2d = Geom2dBSplineCurve::new(
                            vec![pdeb.x(), pfin.x()],
                            vec![pdeb.y(), pfin.y()],
                            vec![udeb, udeb, ufin, ufin],
                            1,
                        )
                        .map_err(|e| e.to_string())?;
                        let trimmed = Arc::new(Geom2dTrimmedCurve::new(Arc::new(bs2d), udeb, ufin));
                        pcur2d = Some(trimmed.clone());
                        h_pcur = trimmed;
                        only3d = true;
                    } else {
                        // `cxx:315`: `HPCur = HProjector`.
                        h_pcur = Arc::new(ProjectedCurve2d::new(projector.clone()));
                    }

                    // `cxx:317-320`.
                    if (!self.with3d || elementary)
                        && projector.max_distance(k).unwrap_or(-1.0) <= self.tol3d
                    {
                        only2d = true;
                    }

                    let mut degenerated = false;
                    let prj: Edge = if only2d && only3d {
                        // `cxx:324-329`: `BRepLib_MakeEdge(GeomAdaptor::MakeCurve(*hcur),
                        // Udeb, Ufin)` then `UpdateEdge(prj, PCur2d, Face(j), myTol3d)`.
                        let e1 = make_edge_over(&builder, &curve3d, udeb, ufin);
                        let pc = pcur2d
                            .clone()
                            .ok_or_else(|| "BRepAlgo_NormalProjection: null p-curve".to_string())?;
                        write_edge_pcurve(&e1, pc, f, self.tol3d);
                        update_edge_vertex_tolerances(&e1, self.tol3d);
                        e1
                    } else {
                        // `cxx:335-341`: the approximation and its `1.e3` reject.
                        let mut appr =
                            ApproxCurveOnSurface::new(h_pcur, surf.clone(), udeb, ufin, self.tol3d);
                        appr.perform(
                            self.max_seg,
                            self.max_degree,
                            self.continuity,
                            only3d,
                            only2d,
                        )?;
                        if appr.max_error3d() > 1e3 * self.tol3d {
                            continue;
                        }
                        let err3d = appr.max_error3d();
                        if !only3d {
                            // `cxx:357-360`.
                            pcur2d = appr.curve2d().map(|c| {
                                Arc::new(c.clone()) as Arc<dyn Curve2d>
                            });
                        }
                        let e1 = if only2d {
                            // `cxx:361-365`.
                            let e1 = make_edge_over(&builder, &curve3d, udeb, ufin);
                            write_edge_pcurve(
                                &e1,
                                pcur2d
                                    .clone()
                                    .ok_or_else(|| {
                                        "BRepAlgo_NormalProjection: null p-curve".to_string()
                                    })?,
                                f,
                                err3d,
                            );
                            e1
                        } else {
                            // `cxx:366-437`.
                            let bs3d = appr
                                .curve3d()
                                .cloned()
                                .ok_or_else(|| {
                                    "BRepAlgo_NormalProjection: null 3D curve".to_string()
                                })?;
                            let bs3d_curve: Arc<dyn Curve> = Arc::new(bs3d.clone());
                            let fp = bs3d_curve.first_parameter();
                            let lp = bs3d_curve.last_parameter();
                            let (deg, point) = is_degenerated_3d(bs3d_curve.as_ref(), fp, lp, self.tol3d);
                            degenerated = deg;
                            let e1 = if deg {
                                // `cxx:429-437`: `BB.MakeVertex(VV)` +
                                // `UpdateVertex(VV, P1, myTol3d)`,
                                // `BB.MakeEdge(prj)` and the two vertex
                                // occurrences (FORWARD, then REVERSED), then
                                // `BB.Degenerated(prj, true)`.
                                let vv = builder.make_vertex(point, self.tol3d);
                                let mut e1 = builder.make_degenerated_edge(&point);
                                builder.add(
                                    &mut e1.0,
                                    &vv.0.oriented(crate::abs::Orientation::Forward),
                                );
                                builder.add(
                                    &mut e1.0,
                                    &vv.0.oriented(crate::abs::Orientation::Reversed),
                                );
                                e1
                            } else {
                                // `cxx:438-439`: `BRepLib_MakeEdge(BS3d)`.
                                make_edge_over(&builder, &bs3d_curve, fp, lp)
                            };
                            // `cxx:443-448`: `UpdateEdge(prj, PCur2d, Face(j),
                            // appr.MaxError3d())` and the two endpoint vertex
                            // tolerances.
                            //
                            // `cxx:449-451` (`BB.Range(prj, Face(j), Udeb,
                            // Ufin)`) has no counterpart and needs none: the
                            // p-curve of `cxx:357-360` is the approximation's
                            // own 2D curve, built over `[Udeb, Ufin]`
                            // (`cxx:335-336`), so the range OCCT would record
                            // is the curve's own range.
                            let pc = pcur2d.clone().ok_or_else(|| {
                                "BRepAlgo_NormalProjection: null p-curve".to_string()
                            })?;
                            write_edge_pcurve(&e1, pc, f, err3d);
                            update_edge_vertex_tolerances(&e1, err3d);
                            e1
                        };
                        e1
                    };

                    if self.face_bounds {
                        // `cxx:454-528`.
                        if !degenerated {
                            // `cxx:465-478`: `BRepAlgoAPI_Section(Faces->Value(j), prj)`.
                            let mut a_section =
                                SectionAlgo::with_shapes(&f.0, &prj.0, false);
                            if a_section.build().is_ok() {
                                // `cxx:469-478`: every EDGE of the section result
                                // is added. OCCT does *not* append to DescenList
                                // on this arm.
                                for a_e in
                                    explorer_shapes(a_section.shape(), ShapeType::Edge)
                                {
                                    builder.add(&mut self.res, &a_e);
                                    self.ancestor_map.insert(
                                        crate::bop_occt_util::shape_key(&a_e),
                                        e.0.clone(),
                                    );
                                    self.corresp.insert(
                                        crate::bop_occt_util::shape_key(&a_e),
                                        f.0.clone(),
                                    );
                                }
                            } else {
                                // `cxx:480-498`: classification fallback.
                                self.push_if_inside(
                                    &builder,
                                    f,
                                    &prj,
                                    &pcur2d,
                                    &mut descen_list,
                                    e,
                                )?;
                            }
                        } else {
                            // `cxx:501-527`: the degenerate arm classifies the
                            // midpoint of the p-curve.
                            self.push_if_inside(
                                &builder,
                                f,
                                &prj,
                                &pcur2d,
                                &mut descen_list,
                                e,
                            )?;
                        }
                    } else {
                        // `cxx:529-535`.
                        builder.add(&mut self.res, &prj.0);
                        descen_list.push(prj.0.clone());
                        self.ancestor_map
                            .insert(crate::bop_occt_util::shape_key(&prj.0), e.0.clone());
                        self.corresp
                            .insert(crate::bop_occt_util::shape_key(&prj.0), f.0.clone());
                    }
                }
            }
            // `cxx:537`: `myDescendants.Bind(Edges->Value(i), DescenList)`.
            self.descendants
                .insert(crate::bop_occt_util::shape_key(e.shape()), descen_list);
        }

        // `cxx:544-548`.
        if ya_vertex_res {
            builder.add(&mut self.res, &vertex_res);
        }
        self.is_done = true;
        Ok(())
    }

    /// `BuildWire(ListOfWire)` (`cxx:642-677`).
    pub fn build_wire(&self) -> Option<Wire> {
        let edges = edges_of(&self.res);
        if edges.is_empty() {
            return None;
        }
        let mut mw = MakeWire::new();
        for e in &edges {
            let _ = mw.add(e);
        }
        if !mw.is_done() {
            return None;
        }
        let wire = mw.wire();
        // `cxx:663-671`: the wire must hold exactly the edges handed in.
        if edges_of(wire.shape()).len() == edges.len() {
            Some(wire)
        } else {
            None
        }
    }
}

impl Default for NormalProjection {
    fn default() -> Self {
        Self::new()
    }
}

/// `cxx:366-425`: the degeneracy scan. Three points, then - if they all lie
/// inside `tol3d` of their own centroid - a second pass with ten points.
/// Answers `(degenerated, centroid)`.
fn is_degenerated_3d(c: &dyn Curve, first: f64, last: f64, tol3d: f64) -> (bool, GpPnt) {
    let pass = |nb_point: usize| -> (bool, GpPnt) {
        let dpar = (last - first) / (nb_point as f64 - 1.0);
        let mut p1 = GpPnt::new(0.0, 0.0, 0.0);
        for ii in 0..nb_point {
            let par = first + ii as f64 * dpar;
            let pp = c.value(par);
            p1 = GpPnt::new(
                p1.x() + pp.x() / nb_point as f64,
                p1.y() + pp.y() / nb_point as f64,
                p1.z() + pp.z() / nb_point as f64,
            );
        }
        let mut degenerated = true;
        for ii in 0..nb_point {
            let par = first + ii as f64 * dpar;
            let pp = c.value(par);
            if p1.distance(&pp) > tol3d {
                degenerated = false;
                break;
            }
        }
        (degenerated, p1)
    };
    let (deg, p1) = pass(3);
    if deg {
        pass(10)
    } else {
        (false, p1)
    }
}
