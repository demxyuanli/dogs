//! `BRepClass_Intersector` — the 2D line/edge intersector the classifier hands
//! to `AreFacesCoincideInArea`.
//!
//! Source (OCCT 8.0.0):
//! - `ModelingAlgorithms/TKTopAlgo/BRepClass/BRepClass_Intersector.cxx:141-441`
//!   (`CheckOn`, `CheckSkip`, `RefineTolerance`, `Perform`),
//! - `BRepClass_Edge.cxx:27-30` (the defaults `MaxTolerance = Precision::Infinite()`,
//!   `UseBndBox = false`),
//! - `IntRes2d_Intersection`'s `Append` / `SetValues`, which the port gets from
//!   [`occt_core::intres2d`].
//!
//! The single caller in the port is
//! [`crate::topopebrep::build::builder::on_patches::are_faces_coincide_in_area`],
//! which builds a `BRepClass_Edge` with only `Edge()` and `Face()` bound, i.e.
//! `UseBndBox() == false` and `NextEdge()` null.
//!
//! UNPORTED, each guarded by exactly the OCCT condition that makes it
//! unreachable for that caller:
//! - the `anUseBndBox` arms (`cxx:351-358`, `:383-396`): they need
//!   `BndLib_Add2dCurve::Add` and `Bnd_Box2d::IsOut(gp_Lin2d)`, which the port
//!   does not carry. `UseBndBox()` is false for every `BRepClass_Edge`
//!   constructor and nothing in the port calls `SetUseBndBox(true)`.
//! - the body of `CheckSkip` after its `NextEdge().IsNull()` guard
//!   (`cxx:222-230`): the caller never calls `SetNextEdge`.

use occt_core::gp::lin2d::GpLin2d;
use occt_core::gp::vec2d::GpVec2d;
use occt_core::intres2d::{IntRes2dDomain, IntRes2dIntersectionPoint, IntRes2dIntersectionSegment};
use occt_core::precision::{CONFUSION, INFINITE, PCONFUSION, PINTERSECTION};
use occt_geom2d::geom2d_int::Geom2dIntGInter;
use occt_geom2d::line::Geom2dLine;

use crate::geom_bnd_lib_elclib2d::line_value_2d;
use crate::shape::Edge;
use crate::shape::Face;
use crate::topopebrep::build::tools::fc2d_curve_on_surface;

/// `BRepClass_Edge` (`hxx:39-95`) with the fields
/// [`BRepClassIntersector::perform`] reads.
pub struct BRepClassEdge {
    /// `myEdge` (`hxx:88`).
    pub edge: Option<Edge>,
    /// `myFace` (`hxx:89`).
    pub face: Option<Face>,
    /// `myNextEdge` (`hxx:90`); only `SetNextEdge` (`cxx:38-75`) fills it.
    pub next_edge: Option<Edge>,
    /// `myMaxTolerance` (`hxx:91`): `Precision::Infinite()` by default.
    pub max_tolerance: f64,
    /// `myUseBndBox` (`hxx:92`): `false` by default.
    pub use_bnd_box: bool,
}

impl Default for BRepClassEdge {
    /// `BRepClass_Edge()` (`BRepClass_Edge.cxx:27-30`).
    fn default() -> Self {
        Self {
            edge: None,
            face: None,
            next_edge: None,
            max_tolerance: INFINITE,
            use_bnd_box: false,
        }
    }
}

impl BRepClassEdge {
    /// `BRepClass_Edge(E, F)` (`hxx:47`), the two-shape constructor OCCT
    /// declares inline.
    pub fn new(e: Edge, f: Face) -> Self {
        Self {
            edge: Some(e),
            face: Some(f),
            ..Default::default()
        }
    }
}

/// `BRepClass_Intersector`'s `IntRes2d_Intersection` base (`done`, `lpnt`,
/// `lseg`).
pub struct BRepClassIntersector {
    done: bool,
    points: Vec<IntRes2dIntersectionPoint>,
    segments: Vec<IntRes2dIntersectionSegment>,
}

impl Default for BRepClassIntersector {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepClassIntersector {
    pub fn new() -> Self {
        Self {
            done: false,
            points: Vec::new(),
            segments: Vec::new(),
        }
    }

    /// `IsDone()` (`IntRes2d_Intersection.hxx`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `NbPoints()`.
    pub fn nb_points(&self) -> usize {
        self.points.len()
    }

    /// `Point(I)` (1-based in OCCT).
    pub fn point(&self, i: usize) -> Option<&IntRes2dIntersectionPoint> {
        self.points.get(i - 1)
    }

    /// `NbSegments()`.
    pub fn nb_segments(&self) -> usize {
        self.segments.len()
    }

    /// `Segment(I)` (1-based in OCCT).
    pub fn segment(&self, i: usize) -> Option<&IntRes2dIntersectionSegment> {
        self.segments.get(i - 1)
    }

    /// `Perform(L, P, Tol, E)` (`cxx:327-441`).
    pub fn perform(&mut self, l: &GpLin2d, p: f64, tol: f64, e: &BRepClassEdge) {
        self.done = false;
        self.points.clear();
        self.segments.clear();

        let (Some(ee), Some(f)) = (e.edge.as_ref(), e.face.as_ref()) else {
            return;
        };
        let mut a_tol_z = tol;

        // `cxx:342-347`: no p-curve, no answer.
        let Some((c2d, deb, fin)) = fc2d_curve_on_surface(ee, f) else {
            return;
        };

        // UNPORTED (`cxx:350-358`): the bound-box arm; see the module notes.
        debug_assert!(!e.use_bnd_box);

        // `cxx:361-381`: the "ON" case, checked directly against the edge.
        {
            let (a_deb_tol, a_fin_tol) = if a_tol_z > CONFUSION {
                (deb - a_tol_z, fin + a_tol_z)
            } else {
                (deb, fin)
            };
            if let Some(pnt) = check_on(f, l, c2d.as_ref(), &mut a_tol_z, fin, deb, a_deb_tol, a_fin_tol)
            {
                self.points.push(pnt);
                self.done = true;
                return;
            }
        }

        // UNPORTED (`cxx:383-396`): the `anUseBndBox` re-check of the tolerance
        // and of the box against the ray; same reason as above.

        let pdeb = c2d.d0(deb);
        let pfin = c2d.d0(fin);
        let toldeb = 1.0e-5;
        let tolfin = 1.0e-5;

        // `cxx:404-417`: the ray's own domain.
        let a_loc = l.location();
        let dl = if p != f64::MAX {
            IntRes2dDomain::bounded(&a_loc, 0.0, PCONFUSION, &line_value_2d(p, l), p, PCONFUSION)
        } else {
            IntRes2dDomain::semi_infinite(&a_loc, 0.0, PCONFUSION, true)
        };

        // `cxx:419-427`: the edge's domain, made periodic when the curve is.
        let mut de = IntRes2dDomain::bounded(&pdeb, deb, toldeb, &pfin, fin, tolfin);
        if c2d.is_periodic() {
            let first = c2d.first_parameter();
            let last = c2d.last_parameter();
            de.set_equivalent_parameters(first, first + last - first);
        }

        // `cxx:429-430`: the ray as a `Geom2d_Line`.
        let gl = Geom2dLine::new(*l.position());
        let mut inter = Geom2dIntGInter::new();
        inter.perform(&gl, &dl, c2d.as_ref(), &de, PCONFUSION, PINTERSECTION);

        // `cxx:435-438`: retry through a vertex with a high tolerance.
        if inter.nb_points() == 0 && inter.nb_segments() == 0 {
            check_skip(e);
        }

        // `cxx:440`: `SetValues(Inter)`.
        self.done = inter.is_done();
        for i in 1..=inter.nb_points() {
            self.points.push(inter.point(i));
        }
        for i in 1..=inter.nb_segments() {
            self.segments.push(inter.segment(i));
        }
    }
}

/// `CheckOn(thePntInter, theF, theL, theCur, theTolZ, theFin, theDeb)`
/// (`cxx:141-197`).
///
/// `theCur` is the edge's p-curve over the tolerance-expanded domain
/// [`a_deb_tol`, `a_fin_tol`] (`cxx:373`).
///
/// UNPORTED: `Extrema_ExtPC2d` runs over the adaptor's own range
/// (`cxx:148`); the port's [`occt_geom2d::extrema2d::point_curve_extrema2d_all`]
/// runs over the curve's stored range, so the `aTolZ` expansion of `cxx:373-375`
/// does not widen the search.
fn check_on(
    f: &Face,
    l: &GpLin2d,
    cur: &dyn occt_geom2d::curve::Curve2d,
    tol_z: &mut f64,
    fin: f64,
    deb: f64,
    _a_deb_tol: f64,
    _a_fin_tol: f64,
) -> Option<IntRes2dIntersectionPoint> {
    let ext = occt_geom2d::extrema2d::point_curve_extrema2d_all(cur, &l.location());
    let mut a_min_dist = f64::MAX;
    let mut a_min: Option<(f64, f64)> = None;
    for e in ext.iter() {
        let d = e.distance;
        if d < a_min_dist {
            a_min_dist = d;
            a_min = Some((e.u1, 0.0));
        }
    }
    let Some((a_par, _)) = a_min else {
        return None;
    };

    if a_min_dist > *tol_z {
        return None;
    }

    let a_pnt_exact = cur.d0(a_par);
    refine_tolerance(f, cur, a_par, tol_z);
    if a_min_dist > *tol_z {
        return None;
    }

    // `cxx:179-193`: the position of the point on the curve.
    let mut a_pos_on_curve = occt_core::intres2d::IntRes2dPosition::Middle;
    if (a_par - deb).abs() <= CONFUSION || a_par < deb {
        a_pos_on_curve = occt_core::intres2d::IntRes2dPosition::Head;
    } else if (a_par - fin).abs() <= CONFUSION || a_par > fin {
        a_pos_on_curve = occt_core::intres2d::IntRes2dPosition::End;
    }

    let a_tr_on_lin = occt_core::intres2d::IntRes2dTransition::undecided(
        occt_core::intres2d::IntRes2dPosition::Head,
    );
    let a_tr_on_curve =
        occt_core::intres2d::IntRes2dTransition::undecided(a_pos_on_curve);

    Some(IntRes2dIntersectionPoint::with_transitions(
        &a_pnt_exact,
        0.0,
        a_par,
        &a_tr_on_lin,
        &a_tr_on_curve,
        false,
    ))
}

/// `RefineTolerance(aF, aC, aT, aTolZ)` (`cxx:479-517`): only a cylinder has a
/// tolerance that follows the p-curve direction.
fn refine_tolerance(
    f: &Face,
    cur: &dyn occt_geom2d::curve::Curve2d,
    a_t: f64,
    tol_z: &mut f64,
) {
    use crate::topopebrep::tool::adaptor::SurfaceAdaptor;
    use occt_geom::approx_same_parameter::{u_resolution, v_resolution};

    let a_bas = SurfaceAdaptor::of_face(f);
    if a_bas.get_type() != crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType::Cylinder {
        return;
    }
    let Some(surf) = a_bas.surface() else {
        return;
    };
    let a_u_res = u_resolution(surf.as_ref(), *tol_z);
    let a_v_res = v_resolution(surf.as_ref(), *tol_z);
    let (_, a_v2d) = cur.d1(a_t);
    let Ok(a_d2d) = occt_core::gp::dir2d::GpDir2d::new(a_v2d.x(), a_v2d.y()) else {
        return;
    };
    let mut a_tol_x = a_u_res * a_d2d.y.abs() + a_v_res * a_d2d.x.abs();
    if a_tol_x < 0.0 {
        a_tol_x = -a_tol_x;
    }
    if a_tol_x < CONFUSION {
        a_tol_x = CONFUSION;
    }
    if a_tol_x < *tol_z {
        *tol_z = a_tol_x;
    }
}

/// `CheckSkip` (`cxx:200-323`) truncated at its first guard.
fn check_skip(e: &BRepClassEdge) {
    let (Some(edge), Some(face)) = (e.edge.as_ref(), e.face.as_ref()) else {
        return;
    };
    let _ = (edge, face);
    // `cxx:228-230`: without a next edge there is nothing to skip.
    if e.next_edge.is_none() {
        return;
    }
    // UNPORTED (`cxx:231-322`): the high-tolerance vertex skip needs
    // `TopExp::LastVertex`, `GC_MakeSegment2d` and a second `Geom2dInt_GInter`
    // with the segment. No caller in the port sets `NextEdge`.
    let _ = GpVec2d::zero();
}
