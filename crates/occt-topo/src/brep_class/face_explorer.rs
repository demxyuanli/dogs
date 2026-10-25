//! Port of `BRepClass_FaceExplorer` (`BRepClass_FaceExplorer.hxx:26-137`,
//! `BRepClass_FaceExplorer.cxx:26-235`): the face exploration used by
//! `BRepClass_FClassifier` for the UV classification and for the
//! `Segment` / `OtherSegment` probing lines.
//!
//! `Reject`, `RejectWire` and `RejectEdge` are the three constant `false`
//! bodies of `cxx:107-110`, `cxx:186-189` and `cxx:216-219`; they are kept as
//! methods so the caller mirrors `TopClass_FaceClassifier.pxx:47-48`,
//! `:75-77` and `:80-82`.

use occt_core::gp::{GpDir2d, GpLin2d, GpPnt2d, GpVec2d};
use occt_core::precision::{epsilon, Precision, CONFUSION};

use occt_geom2d::curve::Curve2d;
use occt_geom2d::extrema2d::ExtremaExtPC2d;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_surface::face_uv_bounds;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, Vertex};
use crate::topo_tools::{explorer_shapes, is_same, map_shapes_and_ancestors, shape_map_key};

use super::edge::BRepClassEdge;

/// `static const double Probing_Start` (`cxx:41`).
const PROBING_START: f64 = 0.123;
/// `static const double Probing_End` (`cxx:42`).
const PROBING_END: f64 = 0.7;
/// `static const double Probing_Step` (`cxx:43`).
const PROBING_STEP: f64 = 0.2111;
/// `const double SmallAngle = 0.001` (`cxx:167`).
const SMALL_ANGLE: f64 = 0.001;

/// `BRepClass_FaceExplorer`.
pub struct FaceExplorer {
    face: Face,
    /// `myCurEdgeInd`, 1-based.
    cur_edge_ind: usize,
    /// `myCurEdgePar`.
    cur_edge_par: f64,
    max_tolerance: f64,
    use_bnd_box: bool,
    /// `myWExplorer` / `myEExplorer` cursors. `None` is "not started".
    wire_index: Option<usize>,
    edge_index: Option<usize>,
    /// `myMapVE` built by `InitEdges`, as `(vertex, edges)` pairs.
    map_ve: Vec<(Vertex, Vec<Edge>)>,
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
}

impl FaceExplorer {
    /// `BRepClass_FaceExplorer(const TopoDS_Face&)` (`cxx:45-62`).
    pub fn new(face: &Face) -> Self {
        Self {
            face: Face(face.0.clone()),
            cur_edge_ind: 1,
            cur_edge_par: PROBING_START,
            max_tolerance: 0.1,
            use_bnd_box: false,
            wire_index: None,
            edge_index: None,
            map_ve: Vec::new(),
            u_min: f64::INFINITY,
            u_max: f64::NEG_INFINITY,
            v_min: f64::INFINITY,
            v_max: f64::NEG_INFINITY,
        }
    }

    /// `ComputeFaceBounds` (`cxx:66-77`): the surface's own range, and
    /// `BRepTools::UVBounds` when that range is infinite.
    fn compute_face_bounds(&mut self) {
        let (mut umin, mut umax, mut vmin, mut vmax) = match BRepTool::face_surface(&self.face) {
            Some(s) => {
                let (a, b) = s.u_range();
                let (c, d) = s.v_range();
                (a, b, c, d)
            }
            None => (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY),
        };
        if Precision::is_infinite(umin) || Precision::is_infinite(umax) || Precision::is_infinite(vmin) || Precision::is_infinite(vmax) {
            let (a, b, c, d) = face_uv_bounds(&self.face);
            umin = a;
            umax = b;
            vmin = c;
            vmax = d;
        }
        self.u_min = umin;
        self.u_max = umax;
        self.v_min = vmin;
        self.v_max = vmax;
    }

    /// `CheckPoint(gp_Pnt2d&)` (`cxx:81-116`). Returns `true` when the point
    /// was not changed.
    pub fn check_point(&mut self, the_point: &mut GpPnt2d) -> bool {
        if self.u_min > self.u_max {
            self.compute_face_bounds();
        }

        if Precision::is_infinite(self.u_min)
            || Precision::is_infinite(self.u_max)
            || Precision::is_infinite(self.v_min)
            || Precision::is_infinite(self.v_max)
        {
            return true;
        }

        let center = GpPnt2d::new(
            (self.u_min + self.u_max) / 2.0,
            (self.v_min + self.v_max) / 2.0,
        );
        let distance = center.distance(the_point);
        if Precision::is_infinite(distance) {
            the_point.set_coord(
                self.u_min - (self.u_max - self.u_min),
                self.v_min - (self.v_max - self.v_min),
            );
            return false;
        }

        let an_epsilon = epsilon(distance);
        if an_epsilon > (self.u_max - self.u_min).max(self.v_max - self.v_min) {
            // `gp_Vec2d aLinVec(aCenterPnt, thePoint); gp_Dir2d aLinDir(aLinVec);`
            // A null vector throws in `gp_Dir2d`; `anEpsilon > 0` means
            // `distance > 0`, so the vector is never null here.
            let vec = GpVec2d::new(the_point.x() - center.x(), the_point.y() - center.y());
            let dir = GpDir2d::from_vec2d(&vec).expect("CheckPoint: non-null center vector");
            the_point.set_coord(
                center.x() + dir.x() * (2.0 * an_epsilon),
                center.y() + dir.y() * (2.0 * an_epsilon),
            );
            return false;
        }

        true
    }

    /// `Reject(const gp_Pnt2d&)` (`cxx:120-123`).
    pub fn reject(&self, _p: &GpPnt2d) -> bool {
        false
    }

    /// `Segment` (`cxx:127-133`): restarts the scan at edge 1, probe 0.123.
    pub fn segment(&mut self, p: &GpPnt2d, l: &mut GpLin2d, par: &mut f64) -> bool {
        self.cur_edge_ind = 1;
        self.cur_edge_par = PROBING_START;
        self.other_segment(p, l, par)
    }

    /// `OtherSegment` (`cxx:137-260`).
    pub fn other_segment(&mut self, p: &GpPnt2d, l: &mut GpLin2d, par: &mut f64) -> bool {
        let tol_par_conf2 = CONFUSION * CONFUSION;
        let edges = explorer_shapes(&self.face.0, ShapeType::Edge);

        // `for (i = 1; anExpF.More(); anExpF.Next(), i++)`
        for (index0, shape) in edges.iter().enumerate() {
            let i = index0 + 1;
            if i != self.cur_edge_ind {
                continue;
            }

            let orientation = shape.orientation();
            if orientation == Orientation::Forward || orientation == Orientation::Reversed {
                let Some(edge) = Edge::wrap(shape.clone()) else {
                    continue;
                };
                let Some((c2d, mut f_par, mut l_par)) = curve_on_surface_range(&edge, &self.face)
                else {
                    // `if (!aC2d.IsNull()) { ... }` false: fall through to the
                    // edge advance below.
                    self.cur_edge_ind += 1;
                    self.cur_edge_par = PROBING_START;
                    continue;
                };

                // Treatment of infinite cases (`cxx:166-179`).
                if Precision::is_negative_infinite(f_par) {
                    if Precision::is_positive_infinite(l_par) {
                        f_par = -1.0;
                        l_par = 1.0;
                    } else {
                        f_par = l_par - 1.0;
                    }
                } else if Precision::is_positive_infinite(l_par) {
                    l_par = f_par + 1.0;
                }

                // `for (; myCurEdgePar < Probing_End; myCurEdgePar += Probing_Step)`:
                // the labelled block is the body, so `continue` in the C++ body
                // still reaches the loop increment.
                while self.cur_edge_par < PROBING_END {
                    'probe: {
                        let param_in =
                            self.cur_edge_par * f_par + (1.0 - self.cur_edge_par) * l_par;

                        let (p_on_c, tan_vec) = c2d.d1(param_in);
                        *par = p_on_c.square_distance(p);

                        if *par <= tol_par_conf2 {
                            break 'probe;
                        }
                        // `gp_Vec2d aLinVec(P, aPOnC); gp_Dir2d aLinDir(aLinVec);`
                        let lin_dir = GpDir2d::from_vec2d(&GpVec2d::new(
                            p_on_c.x() - p.x(),
                            p_on_c.y() - p.y(),
                        ))
                        .expect("OtherSegment: non-null probe vector");

                        let tan_mod = tan_vec.x() * tan_vec.x() + tan_vec.y() * tan_vec.y();
                        if tan_mod < tol_par_conf2 {
                            break 'probe;
                        }
                        let tan_norm = tan_mod.sqrt();
                        let tx = tan_vec.x() / tan_norm;
                        let ty = tan_vec.y() / tan_norm;
                        let sin_a = tx * lin_dir.y() - ty * lin_dir.x();
                        let mut is_small_angle = false;
                        if sin_a.abs() < SMALL_ANGLE {
                            is_small_angle = true;
                            if self.cur_edge_par + PROBING_STEP < PROBING_END {
                                break 'probe;
                            }
                        }

                        *l = GpLin2d::from_pnt_dir(*p, lin_dir);

                        // Check if ends of a curve lie on a line.
                        let fp_on_c = c2d.d0(f_par);
                        if square_distance_lin(l, &fp_on_c) <= tol_par_conf2 {
                            break 'probe;
                        }
                        let lp_on_c = c2d.d0(l_par);
                        if square_distance_lin(l, &lp_on_c) <= tol_par_conf2 {
                            break 'probe;
                        }

                        let mut par_cur = *par;
                        if is_small_angle {
                            // Try to find minimal distance between curve and line.
                            if let Some((min_dist, nearest)) =
                                project_on_curve_range(p, c2d.as_ref(), f_par, l_par)
                            {
                                let f_dist = p.square_distance(&fp_on_c);
                                let l_dist = p.square_distance(&lp_on_c);
                                let mut min_dist = min_dist;
                                let mut nearest = nearest;
                                if min_dist > f_dist {
                                    min_dist = f_dist;
                                    nearest = fp_on_c;
                                }
                                if min_dist > l_dist {
                                    min_dist = l_dist;
                                    nearest = lp_on_c;
                                }
                                if min_dist < par_cur {
                                    par_cur = min_dist;
                                    if par_cur < tol_par_conf2 {
                                        break 'probe;
                                    }
                                    let dir = GpDir2d::from_vec2d(&GpVec2d::new(
                                        nearest.x() - p.x(),
                                        nearest.y() - p.y(),
                                    ));
                                    if let Ok(dir) = dir {
                                        *l = GpLin2d::from_pnt_dir(*p, dir);
                                    }
                                }
                            }
                        }
                        self.cur_edge_par += PROBING_STEP;
                        if self.cur_edge_par >= PROBING_END {
                            self.cur_edge_ind += 1;
                            self.cur_edge_par = PROBING_START;
                        }

                        *par = par_cur.sqrt();
                        return true;
                    }
                    self.cur_edge_par += PROBING_STEP;
                }
            }

            // This curve is not valid for line construction. Go to another edge.
            self.cur_edge_ind += 1;
            self.cur_edge_par = PROBING_START;
        }

        // nothing found, return an horizontal line
        *par = f64::MAX;
        *l = GpLin2d::from_pnt_dir(*p, GpDir2d::new(1.0, 0.0).expect("unit X direction"));
        false
    }

    /// `InitWires` (`cxx:266-269`).
    pub fn init_wires(&mut self) {
        self.wire_index = Some(0);
    }

    /// `MoreWires`.
    pub fn more_wires(&self) -> bool {
        match self.wire_index {
            Some(i) => i < self.wire_shapes().len(),
            None => false,
        }
    }

    /// `NextWire`.
    pub fn next_wire(&mut self) {
        if let Some(i) = self.wire_index {
            self.wire_index = Some(i + 1);
        }
    }

    /// `RejectWire` (`cxx:274-277`).
    pub fn reject_wire(&self, _l: &GpLin2d, _par: f64) -> bool {
        false
    }

    /// `InitEdges` (`cxx:282-287`): start the edge scan of the current wire
    /// and rebuild `myMapVE` from it.
    pub fn init_edges(&mut self) {
        self.edge_index = Some(0);
        self.map_ve.clear();
        let Some(wire) = self.current_wire() else {
            return;
        };
        let map = map_shapes_and_ancestors(&wire.0, ShapeType::Vertex, ShapeType::Edge);
        for shape in explorer_shapes(&wire.0, ShapeType::Vertex) {
            let key = shape_map_key(&shape);
            let Some(ancestors) = map.get(&key) else {
                continue;
            };
            let Some(vertex) = Vertex::wrap(shape.clone()) else {
                continue;
            };
            if self.map_ve.iter().any(|(v, _)| is_same(v, &vertex)) {
                continue;
            }
            let edges: Vec<Edge> = ancestors.iter().filter_map(|a| Edge::wrap(a.clone())).collect();
            self.map_ve.push((vertex, edges));
        }
    }

    /// `MoreEdges`.
    pub fn more_edges(&self) -> bool {
        match self.edge_index {
            Some(i) => i < self.current_wire_edges().len(),
            None => false,
        }
    }

    /// `NextEdge`.
    pub fn next_edge(&mut self) {
        if let Some(i) = self.edge_index {
            self.edge_index = Some(i + 1);
        }
    }

    /// `RejectEdge` (`cxx:293-296`).
    pub fn reject_edge(&self, _l: &GpLin2d, _par: f64) -> bool {
        false
    }

    /// `CurrentEdge(BRepClass_Edge&, TopAbs_Orientation&)` (`cxx:299-307`).
    pub fn current_edge(&self) -> Option<(BRepClassEdge, Orientation)> {
        let edges = self.current_wire_edges();
        let index = self.edge_index?;
        let shape = edges.get(index)?;
        let edge = Edge::wrap(shape.clone())?;
        let orientation = edge.orientation();
        let mut out = BRepClassEdge::from_edge_face(edge, Face(self.face.0.clone()));
        out.set_next_edge(&self.map_ve);
        out.set_max_tolerance(self.max_tolerance);
        out.set_use_bnd_box(self.use_bnd_box);
        Some((out, orientation))
    }

    /// `MaxTolerance()`.
    pub fn max_tolerance(&self) -> f64 {
        self.max_tolerance
    }

    /// `SetMaxTolerance(theValue)`.
    pub fn set_max_tolerance(&mut self, value: f64) {
        self.max_tolerance = value;
    }

    /// `UseBndBox()`.
    pub fn use_bnd_box(&self) -> bool {
        self.use_bnd_box
    }

    /// `SetUseBndBox(theValue)`.
    pub fn set_use_bnd_box(&mut self, value: bool) {
        self.use_bnd_box = value;
    }

    fn wire_shapes(&self) -> Vec<crate::shape::TopoShape> {
        explorer_shapes(&self.face.0, ShapeType::Wire)
    }

    fn current_wire(&self) -> Option<crate::shape::Wire> {
        let index = self.wire_index?;
        let shapes = self.wire_shapes();
        crate::shape::Wire::wrap(shapes.get(index)?.clone())
    }

    fn current_wire_edges(&self) -> Vec<crate::shape::TopoShape> {
        match self.current_wire() {
            Some(wire) => explorer_shapes(&wire.0, ShapeType::Edge),
            None => Vec::new(),
        }
    }
}

/// `P.SquareDistance(aPOnC)` where `P` is on the line: the port evaluates the
/// squared distance from the line's normalised coefficients
/// (`gp_Lin2d::Coefficients`, `gp_Lin2d.hxx:91-96`).
fn square_distance_lin(l: &GpLin2d, p: &GpPnt2d) -> f64 {
    let d = l.distance(p);
    d * d
}

/// `Geom2dAPI_ProjectPointOnCurve::Init(P, C, Uinf, Usup)` followed by
/// `NbPoints() > 0`, `LowerDistance()` and `NearestPoint()`
/// (`Geom2dAPI_ProjectPointOnCurve.cxx:59-96`): `Extrema_ExtPC2d` over the
/// restricted range, the smallest `SquareDistance` with the first index
/// winning ties.
fn project_on_curve_range(
    p: &GpPnt2d,
    c: &dyn Curve2d,
    u_inf: f64,
    u_sup: f64,
) -> Option<(f64, GpPnt2d)> {
    let ext = ExtremaExtPC2d::new_range(c, p, u_inf, u_sup);
    if !ext.is_done() || ext.nb_ext() == 0 {
        return None;
    }
    let mut best = ext.square_distance(1);
    let mut best_index = 1usize;
    for i in 2..=ext.nb_ext() {
        let d2 = ext.square_distance(i);
        if d2 < best {
            best = d2;
            best_index = i;
        }
    }
    let (_u, point) = ext.point(best_index);
    Some((best, point))
}
