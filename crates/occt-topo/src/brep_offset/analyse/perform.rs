//! `BRepOffset_Analyse::Perform` and the two file-local helpers of
//! `BRepOffset_Analyse.cxx`: `EdgeAnalyse` (`cxx:81-144`) and
//! `CheckMixedContinuity` (`cxx:148-271`), plus `BuildAncestors`
//! (`cxx:276-284`).

use occt_core::kernel::geomabs::Shape;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_lib_continuity::stored_continuity;
use crate::brep_tool::BRepTool;
use crate::brep_tools::is_really_closed;
use crate::chfi3d::{define_connect_type, is_tangent_faces};
use crate::chfi_ds::TypeOfConcavity;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::local_analysis::SurfaceContinuity;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{
    explorer_shapes, explorer_shapes_avoiding, is_same, is_same_occurrence, shape_map_key,
};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::{Analyse, Interval};

/// `BRep_Tool::Range(E, F, f, l)` (`BRep_Tool.cxx:996-1001` ->
/// `Range(E, S, L, f, l)` at `cxx:964-992`): the range of the edge's
/// curve-on-surface representation of `face`, and the 3-D edge range when that
/// representation is absent. Orientation plays no part.
fn range_on_face(e: &Edge, f: &Face) -> (f64, f64) {
    let face_key = GeometryRegistry::shape_key(&f.0);
    match GeometryRegistry::global().pcurve_range(&e.0, face_key) {
        Some(range) => range,
        None => BRepTool::edge_parameters(e),
    }
}

/// `aSurfType == GeomAbs_BSplineSurface || aSurfType == GeomAbs_BezierSurface`
/// (`cxx:116-118`).
fn is_spline(t: GeomAbsSurfaceType) -> bool {
    matches!(
        t,
        GeomAbsSurfaceType::BSplineSurface | GeomAbsSurfaceType::BezierSurface
    )
}

impl Analyse {
    /// `Perform(S, Angle)` (`cxx:288-367`).
    ///
    /// `Err` carries the `Standard_Failure` OCCT lets escape from
    /// `ChFi3d::DefineConnectType` / `ChFi3d::IsTangentFaces` /
    /// `TreatTangentFaces`. The `Message_ProgressRange` early returns
    /// (`cxx:310-313`, `:361-365`) only fire on a user break; there is no
    /// progress range here.
    pub fn perform(&mut self, s: &TopoShape, angle: f64) -> Result<(), String> {
        self.shape = Some(s.clone());
        self.new_faces.clear();
        self.generated.clear();
        self.replacement.clear();
        self.descendants.borrow_mut().clear();

        self.angle = angle;
        let sin_tol = angle.sin().abs();

        // Build ancestors (`cxx:302`).
        self.build_ancestors(s);

        let mut a_letang: Vec<TopoShape> = Vec::new();
        for cur in explorer_shapes(&s.oriented(Orientation::Forward), ShapeType::Edge) {
            let e = Edge(cur);
            if self.map_edge_type.contains_key(&shape_map_key(&e.0)) {
                continue;
            }
            // `myMapEdgeType.Bind(E, LI)` with an empty list (`cxx:317-319`).
            self.bind_edge_types(&e, Vec::new());

            let l = self.ancestors(&e.0).to_vec();
            if l.is_empty() {
                continue;
            }

            if l.len() == 2 {
                let f1 = Face(l[0].clone());
                let f2 = Face(l[1].clone());
                self.edge_analyse(&e, &f1, &f2, sin_tol)?;

                // For tangent faces add artificial perpendicular face to close
                // the gap between them (`cxx:332-338`).
                let is_tangent = self
                    .edge_types(&e)
                    .and_then(|v| v.last())
                    .map(|i| i.kind() == TypeOfConcavity::Tangential)
                    .unwrap_or(false);
                if is_tangent {
                    a_letang.push(e.0.clone());
                }
            } else if l.len() == 1 {
                let f = Face(l[0].clone());
                let (u1, u2) = range_on_face(&e, &f);
                let mut inter = Interval::new(u1, u2, TypeOfConcavity::Other);
                if !is_really_closed(&e, &f) {
                    inter.set_type(TypeOfConcavity::FreeBound);
                }
                self.edge_types_mut(&e).push(inter);
            }
            // `cxx:352-357`: an edge shared by more than two faces is only
            // reported under `OCCT_DEBUG`.
        }

        self.treat_tangent_faces(&a_letang)?;
        self.done = true;
        Ok(())
    }

    /// `BuildAncestors(S, MA)` (`cxx:276-284`): `MapShapesAndUniqueAncestors`
    /// for `(VERTEX, EDGE)` then `(EDGE, FACE)`, on one map.
    fn build_ancestors(&mut self, s: &TopoShape) {
        self.ancestors.clear();
        self.map_unique_ancestors(s, ShapeType::Vertex, ShapeType::Edge);
        self.map_unique_ancestors(s, ShapeType::Edge, ShapeType::Face);
    }

    /// `TopExp::MapShapesAndUniqueAncestors(S, TS, TA, M)`
    /// (`TopExp.cxx:124-177`) with `useOrientation = false`: one ancestor per
    /// sub-shape occurrence, `IsSame` comparison, and an empty list for every
    /// `TS` shape reachable without passing through a `TA` shape.
    fn map_unique_ancestors(&mut self, s: &TopoShape, ts: ShapeType, ta: ShapeType) {
        for anc in explorer_shapes(s, ta) {
            for sub in explorer_shapes(&anc, ts) {
                let list = self.ancestors.change_or_insert(&sub);
                if !list.iter().any(|x| is_same_occurrence(x, &anc)) {
                    list.push(anc.clone());
                }
            }
        }
        for sub in explorer_shapes_avoiding(s, ts, ta) {
            self.ancestors.change_or_insert(&sub);
        }
    }

    /// `EdgeAnalyse(E, F1, F2, SinTol, LI)` (`cxx:81-144`).
    ///
    /// For two planes only the strong `BRep_Tool::Continuity` test is used
    /// (`cxx:103-114`); otherwise a pair of B-spline / Bezier surfaces first
    /// runs `CheckMixedContinuity` with the hard-coded angle tolerance `0.1`
    /// (`cxx:122-123`), then the weak `ChFi3d::IsTangentFaces` test
    /// (`cxx:128`).
    pub(super) fn edge_analyse(
        &mut self,
        e: &Edge,
        f1: &Face,
        f2: &Face,
        sin_tol: f64,
    ) -> Result<(), String> {
        let (f, l) = range_on_face(e, f1);
        let mut i = Interval::new(f, l, TypeOfConcavity::Other);

        let type1 = SurfaceAdaptor::of_face(f1).get_type();
        let type2 = SurfaceAdaptor::of_face(f2).get_type();
        let is_two_planes =
            type1 == GeomAbsSurfaceType::Plane && type2 == GeomAbsSurfaceType::Plane;

        let connect_type = if is_two_planes {
            if stored_continuity(e, f1, f2) > Shape::C0 as u8 {
                TypeOfConcavity::Tangential
            } else {
                define_connect_type(e, f1, f2, sin_tol, false)?
            }
        } else {
            let is_two_splines = is_spline(type1) && is_spline(type2);
            let mut is_mixed_concavity = false;
            if is_two_splines {
                let an_ang_tol = 0.1;
                is_mixed_concavity = check_mixed_continuity(e, f1, f2, an_ang_tol)?;
            }

            if !is_mixed_concavity {
                if is_tangent_faces(e, f1, f2, Shape::G1)? {
                    TypeOfConcavity::Tangential
                } else {
                    define_connect_type(e, f1, f2, sin_tol, false)?
                }
            } else {
                TypeOfConcavity::Mixed
            }
        };

        i.set_type(connect_type);
        self.edge_types_mut(e).push(i);
        Ok(())
    }
}

/// `CheckMixedContinuity(theEdge, theFace1, theFace2, theAngTol)`
/// (`cxx:148-271`): whether the two spline faces are convex in some sampled
/// regions along the edge and concave in others.
///
/// `LocalAnalysis_SurfaceContinuity` is built with `EpsNul = 0.001`,
/// `EpsC0 = TolC0`, `EpsC1 = EpsC2 = EpsG1 = theAngTol` (`cxx:237-248`) and the
/// header defaults `Percent = 0.01`, `Maxlen = 10000`
/// (`LocalAnalysis_SurfaceContinuity.hxx:146-147`).
fn check_mixed_continuity(
    edge: &Edge,
    face1: &Face,
    face2: &Face,
    ang_tol: f64,
) -> Result<bool, String> {
    let mut a_mixed_cont = false;
    if stored_continuity(edge, face1, face2) > Shape::C0 as u8 {
        // `cxx:156-162`: the faces are tangent along the whole edge.
        return Ok(a_mixed_cont);
    }
    // `cxx:163`: a stored C0 is the default value and cannot be trusted.
    let tol_c0 = 0.001_f64.max(1.5 * BRepTool::edge_tolerance(edge));

    let (c2d1, c2d2, a_first, a_last);
    if !is_same(&face1.0, &face2.0)
        && BRepTool::is_closed_edge_face(edge, face1)
        && BRepTool::is_closed_edge_face(edge, face2)
    {
        // Find the edge in the face 1: this edge will have correct orientation
        // (`cxx:170-189`).
        let mut a_face1 = face1.clone();
        a_face1.0.set_orientation(Orientation::Forward);
        let mut found: Option<Edge> = None;
        for cur in explorer_shapes(&a_face1.0, ShapeType::Edge) {
            if is_same(&cur, &edge.0) {
                found = Some(Edge(cur));
                break;
            }
        }
        let Some(mut an_edge_in_face1) = found else {
            return Ok(a_mixed_cont);
        };

        // `cxx:211-214`: a null p-curve answers "not mixed".
        let Some((c1, _, _)) = curve_on_surface_range(&an_edge_in_face1, &a_face1) else {
            return Ok(a_mixed_cont);
        };
        let mut a_face2 = face2.clone();
        a_face2.0.set_orientation(Orientation::Forward);
        an_edge_in_face1.0.reverse();
        let Some((c2, first2, last2)) = curve_on_surface_range(&an_edge_in_face1, &a_face2) else {
            return Ok(a_mixed_cont);
        };
        c2d1 = c1;
        c2d2 = c2;
        a_first = first2;
        a_last = last2;
    } else {
        // Obtaining of pcurves of edge on two faces (`cxx:198-209`).
        let Some((c1, _, _)) = curve_on_surface_range(edge, face1) else {
            return Ok(a_mixed_cont);
        };
        let mut ee = edge.clone();
        if is_same(&face1.0, &face2.0) {
            ee.0.reverse();
        }
        let Some((c2, first2, last2)) = curve_on_surface_range(&ee, face2) else {
            return Ok(a_mixed_cont);
        };
        c2d1 = c1;
        c2d2 = c2;
        a_first = first2;
        a_last = last2;
    }

    // `cxx:220-223`: a null surface answers "not mixed".
    let Some(a_surf1) = BRepTool::face_surface(face1) else {
        return Ok(a_mixed_cont);
    };
    let Some(a_surf2) = BRepTool::face_surface(face2) else {
        return Ok(a_mixed_cont);
    };

    let a_nb_samples = 23usize;
    let a_delta = (a_last - a_first) / (a_nb_samples - 1) as f64;
    let mut a_has_convex = false;
    let mut a_has_concave = false;
    let mut a_nb_valid = 0;

    // `cxx:234-260`: convex in some regions and concave in others.
    for i in 1..=a_nb_samples {
        let a_par = if i == a_nb_samples {
            a_last
        } else {
            a_first + (i - 1) as f64 * a_delta
        };

        let a_cont = SurfaceContinuity::from_curves(
            c2d1.as_ref(),
            c2d2.as_ref(),
            a_par,
            a_surf1.as_ref(),
            a_surf2.as_ref(),
            Shape::G1,
            0.001,
            tol_c0,
            ang_tol,
            ang_tol,
            ang_tol,
            0.01,
            10000.0,
        )?;
        if !a_cont.is_done() {
            continue;
        }

        a_nb_valid += 1;

        if !a_cont.is_g1()? && (!a_has_convex || !a_has_concave) {
            let an_angle = a_cont.c0_value()?;
            a_has_convex = a_has_convex || (an_angle > std::f64::consts::FRAC_PI_2 + ang_tol);
            a_has_concave = a_has_concave || (an_angle < std::f64::consts::FRAC_PI_2 - ang_tol);
        }
    }

    if a_nb_valid < a_nb_samples / 2 {
        return Ok(a_mixed_cont);
    }

    // Mixed connectivity: both convex and concave regions exist (`cxx:267-269`).
    a_mixed_cont = a_has_convex && a_has_concave;

    Ok(a_mixed_cont)
}
