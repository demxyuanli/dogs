//! `BRepLib_FindSurface` - find a surface through the edges of a shape.
//!
//! Source: `BRepLib_FindSurface.cxx:52-621` (`Controle` at 52, `Is2DConnected`
//! at 75, `Is2DClosed` at 107, `fillParams` at 178, `fillPoints` at 217,
//! `Init` at 246, `Found`/`Surface`/`Tolerance`/`ToleranceReached`/`Existed`/
//! `Location` at 583-621) and `BRepLib_FindSurface.hxx:48-101`.
//!
//! Caller: `BRepLib_MakeFace(const TopoDS_Wire&, bool OnlyPlane)`
//! (`BRepLib_MakeFace.cxx:188-256`) reaches `Init(S, -1, OnlyPlane, true)`;
//! `BRepBuilderAPI_MakeFace(const TopoDS_Wire&, bool OnlyPlane)`
//! (`BRepBuilderAPI_MakeFace.cxx:215-224`) delegates to it, and
//! `BRepOffsetAPI_MakeOffset::BuildDomains` (`BRepOffsetAPI_MakeOffset.cxx:184-192`)
//! calls that with `OnlyPlane = true` for the spine wire of an offset.

use std::sync::Arc;

use occt_core::gp::{GpAx3, GpDir, GpPln, GpPnt};
use occt_core::precision::{epsilon, Precision, REAL_SMALL, RESOLUTION};
use occt_geom::{Curve, GeomPlane, Surface};
use occt_geom2d::curve::Curve2d;
use occt_math::{Jacobi, MathMatrix};

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::{FClass2d, FaceState};
use crate::int_tools_wline::{u_resolution, v_resolution};
use crate::shape::{Edge, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edge_vertices, edges_of, is_same, wires_of};

/// `GeomAbs_C3` (`GeomAbs_Shape.hxx:33`), the continuity
/// `BRepAdaptor_Curve::NbIntervals(GeomAbs_C3)` is asked for at
/// `BRepLib_FindSurface.cxx:410`.
const GEOM_ABS_C3: u8 = 5;

/// `IsEqual` (`Standard_Real.hxx:148-151`): `|a - b| < RealSmall()`.
fn is_equal(a: f64, b: f64) -> bool {
    (a - b).abs() < REAL_SMALL
}

/// `Controle` (`cxx:52-68`): the largest distance from `points` to `plane`.
fn controle(points: &[GpPnt], plane: &GpPln) -> f64 {
    let (a, b, c, d) = plane.coefficients();
    let mut max_dist = 0.0;
    for p in points {
        let dist = (a * p.x() + b * p.y() + c * p.z() + d).abs();
        if dist > max_dist {
            max_dist = dist;
        }
    }
    max_dist
}

/// `BRep_Tool::CurveOnSurface(edge, surface, location, f, l)`: the edge's
/// CurveOnSurface representation whose surface is `surface`, with its COS
/// range. `theLocation` is always identity in the port.
///
/// `null` stands for the null handle OCCT returns when the edge has no
/// representation on that surface.
fn curve_on_surface(
    edge: &Edge,
    surface: &Arc<dyn Surface>,
) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    for (s, pc, f, l) in GeometryRegistry::global().edge_pcurve_reps(&edge.0) {
        if Arc::ptr_eq(&s, surface) {
            return Some((pc, f, l));
        }
    }
    None
}

/// `Is2DConnected` (`cxx:75-99`): the last vertex of `e1` coincides with the
/// first vertex of `e2` in the parametric space of `surface`.
///
/// `TopExp::FirstVertex(theEdge2, true)` is the port's
/// [`edge_vertices`] first element.
fn is_2d_connected(e1: &Edge, e2: &Edge, surface: &Arc<dyn Surface>) -> bool {
    let p1 = match curve_on_surface(e1, surface) {
        Some((c, f, l)) => c.d0(if e1.0.orientation() == Orientation::Forward {
            l
        } else {
            f
        }),
        None => return false,
    };
    let p2 = match curve_on_surface(e2, surface) {
        Some((c, f, l)) => c.d0(if e2.0.orientation() == Orientation::Forward {
            f
        } else {
            l
        }),
        None => return false,
    };
    let (Some(v), _) = edge_vertices(e2) else {
        return false;
    };
    let tol3d = BRepTool::vertex_tolerance(&v);
    let tol2d = u_resolution(surface.as_ref(), tol3d) + v_resolution(surface.as_ref(), tol3d);
    p1.square_distance(&p2) < tol2d * tol2d
}

/// `Is2DClosed` (`cxx:107-154`, `catch` at 151-154): the edges of `shape` form
/// a closed wire in the parametric space of `surface`.
fn is_2d_closed(shape: &TopoShape, surface: &Arc<dyn Surface>) -> bool {
    // TopExp_Explorer(theShape, TopAbs_WIRE); the shape is the spine wire here.
    let Some(wire) = wires_of(shape).into_iter().next() else {
        return false;
    };
    // a tmp face (cxx:122)
    let tmp_face = crate::brep_lib_make_face::make_face_from_surface(
        surface.clone(),
        Precision::PCONFUSION,
    );
    // BRepTools_WireExplorer(aWire, aTmpFace)
    let edges = crate::meshing::model_builder::wire_builder::wire_edges_explorer(&wire, &tmp_face);
    let Some(first) = edges.first() else {
        return false;
    };
    let mut last = first.clone();
    let mut prev = first.clone();
    for e in &edges[1..] {
        if !is_2d_connected(&prev, e, surface) {
            return false;
        }
        last = e.clone();
        prev = e.clone();
    }
    let (Some(v1), _) = edge_vertices(first) else {
        return false;
    };
    let (_, Some(v2)) = edge_vertices(&last) else {
        return false;
    };
    is_same(&v1.0, &v2.0) && is_2d_connected(&last, first, surface)
}

/// `fillParams` (`cxx:178-215`): the sample parameters of one edge.
fn fill_params(knots: &[f64], degree: i32, par_min: f64, par_max: f64, out: &mut Vec<f64>) {
    let p_confusion = Precision::PCONFUSION;
    let mut prev = par_min;
    out.push(prev);
    let nb_p = degree.max(1);
    let mut i = 1usize;
    while i < knots.len() && knots[i - 1] < (par_max - p_confusion) {
        if knots[i] < par_min + p_confusion {
            i += 1;
            continue;
        }
        let step = (knots[i] - knots[i - 1]) / nb_p as f64;
        for k in 1..=nb_p {
            let par = knots[i - 1] + k as f64 * step;
            if par > par_max - p_confusion {
                break;
            }
            if par > prev + p_confusion {
                out.push(par);
                prev = par;
            }
        }
        i += 1;
    }
    out.push(par_max);
}

/// `fillPoints` (`cxx:217-243`): append the sampled points and their weights
/// (the sum of the distances to the neighbouring points *on the same edge*).
fn fill_points(
    curve: &dyn Curve,
    params: &[f64],
    points: &mut Vec<GpPnt>,
    weights: &mut Vec<f64>,
) {
    if params.is_empty() {
        return;
    }
    let mut dist_prev = 0.0;
    let mut p_prev = curve.d0(params[0]);
    let mut p_next = p_prev;
    for i in 0..params.len() {
        let dist_next = if i + 1 < params.len() {
            p_next = curve.d0(params[i + 1]);
            p_prev.distance(&p_next)
        } else {
            0.0
        };
        points.push(p_prev);
        weights.push(dist_prev + dist_next);
        dist_prev = dist_next;
        p_prev = p_next;
    }
}

/// `GeomAdaptor_Curve::GetType()` as `BRepLib_FindSurface.cxx:377-417` reads
/// it, with the data each arm needs.
enum AdaptorKind {
    /// `GeomAbs_BezierCurve`: `GC->FirstParameter()` / `LastParameter()` /
    /// `Degree()` of the *basis* Bezier (`cxx:379-387`).
    Bezier { degree: i32, first: f64, last: f64 },
    /// `GeomAbs_BSplineCurve`: `GC->Knots()` / `Degree()` of the basis
    /// (`cxx:388-391`).
    BSpline { knots: Vec<f64>, degree: i32 },
    /// `GeomAbs_Line`: exactly the two end parameters (`cxx:392-396`).
    Line,
    /// `GeomAbs_Circle` / `Ellipse` / `Hyperbola` / `Parabola`: four points
    /// (`cxx:398-404`).
    Conic,
    /// `default` (`cxx:405-418`): `15 + NbIntervals(GeomAbs_C3)` points.
    Other { nb_intervals: i32 },
}

/// `BRepAdaptor_Curve(edge)` as `GeomAdaptor_Curve::load` builds it
/// (`GeomAdaptor_Curve.cxx:252-254`): a `Geom_TrimmedCurve` is unwrapped to its
/// basis while the range stays in the basis parameter space, so `GetType()`,
/// `Bezier()`, `BSpline()` and `Value()` all resolve to the innermost basis.
struct CurveAdaptor {
    curve: Arc<dyn Curve>,
    first: f64,
    last: f64,
    kind: AdaptorKind,
}

impl CurveAdaptor {
    fn of_edge(edge: &Edge) -> Option<Self> {
        let curve = BRepTool::edge_curve(edge)?;
        let (curve, first, last) = match curve.untrimmed_basis() {
            Some((basis, f, l)) => (basis, f, l),
            None => (curve.clone(), curve.first_parameter(), curve.last_parameter()),
        };
        let kind = if let Some(poles) = curve.bezier_poles() {
            AdaptorKind::Bezier {
                degree: (poles.len() as i32 - 1).max(1),
                first: curve.first_parameter(),
                last: curve.last_parameter(),
            }
        } else if let Some(knots) = curve.bspline_knots() {
            AdaptorKind::BSpline {
                knots: knots.to_vec(),
                degree: curve.nurbs_degree().unwrap_or(1) as i32,
            }
        } else if curve.gp_line().is_some() {
            AdaptorKind::Line
        } else if curve.gp_circ().is_some()
            || curve.gp_ellipse().is_some()
            || curve.gp_hyperbola().is_some()
            || curve.gp_parabola().is_some()
        {
            AdaptorKind::Conic
        } else {
            AdaptorKind::Other {
                nb_intervals: curve.nb_intervals(GEOM_ABS_C3),
            }
        };
        Some(Self {
            curve,
            first,
            last,
            kind,
        })
    }
}

/// `BRepLib_FindSurface` (`BRepLib_FindSurface.hxx:48-101`).
pub struct FindSurface {
    my_tolerance: f64,
    my_tol_reached: f64,
    is_existed: bool,
    my_surface: Option<Arc<dyn Surface>>,
}

impl FindSurface {
    /// `BRepLib_FindSurface()` (`cxx:158-163`).
    pub fn empty() -> Self {
        Self {
            my_tolerance: 0.0,
            my_tol_reached: 0.0,
            is_existed: false,
            my_surface: None,
        }
    }

    /// `BRepLib_FindSurface(S, Tol, OnlyPlane, OnlyClosed)` (`cxx:166-172`).
    pub fn new(shape: &TopoShape, tol: f64, only_plane: bool, only_closed: bool) -> Self {
        let mut fs = Self::empty();
        fs.init(shape, tol, only_plane, only_closed);
        fs
    }

    /// `Found()` (`cxx:583-586`).
    pub fn found(&self) -> bool {
        self.my_surface.is_some()
    }

    /// `Surface()` (`cxx:590-593`).
    pub fn surface(&self) -> Option<Arc<dyn Surface>> {
        self.my_surface.clone()
    }

    /// `Tolerance()` (`cxx:595-598`).
    pub fn tolerance(&self) -> f64 {
        self.my_tolerance
    }

    /// `ToleranceReached()` (`cxx:602-605`).
    pub fn tolerance_reached(&self) -> f64 {
        self.my_tol_reached
    }

    /// `Existed()` (`cxx:611-614`).
    pub fn existed(&self) -> bool {
        self.is_existed
    }

    /// `Init(S, Tol, OnlyPlane, OnlyClosed)` (`cxx:246-578`).
    ///
    /// `myLocation.Identity()` (`cxx:255`) is the identity `TopLoc_Location`
    /// the port always carries.
    pub fn init(&mut self, shape: &TopoShape, tol: f64, only_plane: bool, only_closed: bool) {
        self.my_tolerance = tol;
        self.my_tol_reached = 0.0;
        self.is_existed = false;
        self.my_surface = None;

        // compute the tolerance (cxx:260-268)
        let edges = edges_of(shape);
        for e in &edges {
            let t = BRepTool::edge_tolerance(e);
            if t > self.my_tolerance {
                self.my_tolerance = t;
            }
        }

        // search an existing surface (cxx:270-348)
        let Some(first_edge) = edges.first().cloned() else {
            return; // no edges ....
        };
        let reg = GeometryRegistry::global();
        // iterate on the representations of the first edge
        for (surface, _, _, _) in reg.edge_pcurve_reps(&first_edge.0) {
            let mut candidate = Some(surface.clone());
            // check the other edges
            for other in &edges {
                if is_same(&other.0, &first_edge.0) {
                    continue;
                }
                let mut found = false;
                for (s, _, _, _) in reg.edge_pcurve_reps(&other.0) {
                    if Arc::ptr_eq(&s, &surface) {
                        found = true;
                        break;
                    }
                }
                if !found {
                    candidate = None;
                    break;
                }
            }
            // if OnlyPlane, eval if mySurface is a plane (cxx:324-332).
            // `occ::down_cast<Geom_Plane>` yields a null handle on a type
            // mismatch instead of throwing, so a non-plane candidate simply
            // rejects this representation.
            if only_plane {
                if let Some(s) = candidate.take() {
                    let basis = s.rectangular_trimmed_basis().unwrap_or(s);
                    candidate = if basis.gp_pln().is_some() {
                        Some(basis)
                    } else {
                        None
                    };
                }
            }
            if let Some(s) = candidate {
                // if S is e.g. the bottom face of a cylinder, mySurface can be
                // the lateral (cylindrical) face of the cylinder; reject an
                // improper mySurface (cxx:333-340)
                if !only_closed || is_2d_closed(shape, &s) {
                    self.my_surface = Some(s);
                    break;
                }
            }
        }

        if self.my_surface.is_some() {
            self.is_existed = true;
            return;
        }

        // no existing surface, search a plane (cxx:350-356).
        // 1. Collect the points along all edges of the shape. For each point
        //    calculate the WEIGHT = sum of distances from neighboring points
        //    (_only_ same edge).
        // 2. Minimizing the weighed sum of squared deviations compute
        //    coefficients of the sought plane.
        let mut points: Vec<GpPnt> = Vec::new();
        let mut weights: Vec<f64> = Vec::new();

        // ======================= Step #1 (cxx:361-424)
        for e in &edges {
            let Some(adaptor) = CurveAdaptor::of_edge(e) else {
                continue;
            };
            let (df_uf, df_ul) = (adaptor.first, adaptor.last);
            if is_equal(df_uf, df_ul) {
                continue; // Degenerate
            }
            let mut params: Vec<f64> = Vec::new();
            match &adaptor.kind {
                AdaptorKind::Bezier {
                    degree,
                    first,
                    last,
                } => {
                    fill_params(&[*first, *last], *degree, df_uf, df_ul, &mut params);
                }
                AdaptorKind::BSpline { knots, degree } => {
                    fill_params(knots, *degree, df_uf, df_ul, &mut params);
                }
                AdaptorKind::Line => {
                    params.push(df_uf);
                    params.push(df_ul);
                }
                AdaptorKind::Conic => {
                    fill_params(&[df_uf, df_ul], 3, df_uf, df_ul, &mut params);
                }
                AdaptorKind::Other { nb_intervals } => {
                    let nb_points = 15 + nb_intervals;
                    fill_params(&[df_uf, df_ul], nb_points - 1, df_uf, df_ul, &mut params);
                }
            }
            fill_points(adaptor.curve.as_ref(), &params, &mut points, &mut weights);
        }

        if points.len() < 3 {
            return;
        }

        // ======================= Step #2 (cxx:430-577)
        let mut mat = MathMatrix::new(1, 3, 1, 3);
        // Find the barycenter and normalize weights
        let mut max_weight = 0.0;
        let (mut bary_x, mut bary_y, mut bary_z) = (0.0, 0.0, 0.0);
        let mut sum_weight = 0.0;
        for (p, w) in points.iter().zip(&weights) {
            bary_x += w * p.x();
            bary_y += w * p.y();
            bary_z += w * p.z();
            sum_weight += w;
            if *w > max_weight {
                max_weight = *w;
            }
        }
        bary_x /= sum_weight;
        bary_y /= sum_weight;
        bary_z /= sum_weight;

        // Fill the matrix and the right vector
        for (p, w) in points.iter().zip(&weights) {
            let px = p.x() - bary_x;
            let py = p.y() - bary_y;
            let pz = p.z() - bary_z;
            let w = w / max_weight;
            mat.set_value(1, 1, mat.value(1, 1) + w * px * px);
            mat.set_value(1, 2, mat.value(1, 2) + w * px * py);
            mat.set_value(1, 3, mat.value(1, 3) + w * px * pz);
            mat.set_value(2, 2, mat.value(2, 2) + w * py * py);
            mat.set_value(2, 3, mat.value(2, 3) + w * py * pz);
            mat.set_value(3, 3, mat.value(3, 3) + w * pz * pz);
        }
        mat.set_value(2, 1, mat.value(1, 2));
        mat.set_value(3, 1, mat.value(1, 3));
        mat.set_value(3, 2, mat.value(2, 3));

        let eig = Jacobi::new(&mat);
        let mut is_solved = eig.is_done;
        let mut isol = 0usize;
        let mut vec3 = [0.0f64; 3];
        if is_solved {
            // We need vector with eigenvalue ~ 0.
            let mut e_min = f64::MAX;
            let mut e_max = -f64::MAX;
            for i in 1..=3 {
                let e = eig.values.value(i).abs();
                if e_min > e {
                    e_min = e;
                    isol = i as usize;
                }
                if e_max < e {
                    e_max = e;
                }
            }

            if isol == 0 {
                is_solved = false;
            } else {
                let eps = epsilon(e_max);
                if e_min <= eps {
                    for i in 1..=3 {
                        vec3[i - 1] = eig.vectors.value(i, isol);
                    }
                } else {
                    // try using vector product of other axes
                    let mut ind = [0usize; 2];
                    for i in 1..=3 {
                        if i == isol {
                            continue;
                        }
                        if ind[0] == 0 {
                            ind[0] = i;
                            continue;
                        }
                        if ind[1] == 0 {
                            ind[1] = i;
                            continue;
                        }
                    }
                    let v1 = [
                        eig.vectors.value(1, ind[0]),
                        eig.vectors.value(2, ind[0]),
                        eig.vectors.value(3, ind[0]),
                    ];
                    let v2 = [
                        eig.vectors.value(1, ind[1]),
                        eig.vectors.value(2, ind[1]),
                        eig.vectors.value(3, ind[1]),
                    ];
                    vec3 = [
                        v1[1] * v2[2] - v1[2] * v2[1],
                        v1[2] * v2[0] - v1[0] * v2[2],
                        v1[0] * v2[1] - v1[1] * v2[0],
                    ];
                }
                if vec3.iter().map(|v| v * v).sum::<f64>() < RESOLUTION {
                    is_solved = false;
                }
            }
        }

        if !is_solved {
            return;
        }

        // Removing very small values
        let max_v = vec3.iter().fold(0.0f64, |a, v| a.max(v.abs()));
        let eps = epsilon(max_v);
        for v in vec3.iter_mut() {
            if v.abs() <= eps {
                *v = 0.0;
            }
        }
        let Ok(norm) = GpDir::new(vec3[0], vec3[1], vec3[2]) else {
            return;
        };
        let Ok(ax3) = GpAx3::from_point_dir(GpPnt::new(bary_x, bary_y, bary_z), norm) else {
            return;
        };
        let plane = GpPln::new(ax3);
        let plane_surface: Arc<dyn Surface> = Arc::new(GeomPlane::new(plane.clone()));
        self.my_tol_reached = controle(&points, &plane);
        let weakness = 5.0;
        if self.my_tol_reached <= self.my_tolerance
            || (tol < 0.0 && self.my_tol_reached < self.my_tolerance * weakness)
        {
            self.my_surface = Some(plane_surface.clone());
            // If S is wire, try to orient surface according to orientation of
            // wire (cxx:564-576).
            if shape.shape_type() == ShapeType::Wire && shape.closed() {
                if let Some(w) = wires_of(shape).into_iter().next() {
                    let tmp_face = crate::brep_lib_make_face::make_face_from_surface(
                        plane_surface,
                        Precision::CONFUSION,
                    );
                    let mut tmp_face = tmp_face;
                    TopoBuilder::new().add_wire(&mut tmp_face, &w);
                    if let Ok(fclass) = FClass2d::new(&tmp_face, 0.0) {
                        if fclass.perform_infinite_point() == FaceState::In {
                            let loc = plane.location();
                            let dir = plane.axis().direction().reversed();
                            if let Ok(ax3) = GpAx3::from_point_dir(loc, dir) {
                                self.my_surface =
                                    Some(Arc::new(GeomPlane::new(GpPln::new(ax3))));
                            }
                        }
                    }
                }
            }
        }
    }
}
