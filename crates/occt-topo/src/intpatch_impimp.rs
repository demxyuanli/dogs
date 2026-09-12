//! `IntPatch_ImpImpIntersection` — analytic intersection of two implicit
//! quadrics (plane / cylinder / cone / sphere / torus).
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_ImpImpIntersection.cxx`
//! Perform / SetQuad / IntPP..IntSpSp plus TheSOnBounds, PutPointsOnLine,
//! ProcessSegments and ProcessRLine on UV-box restriction arcs.
//!
//! Pairs without a closed form in `occt_geom::intana` fall back to
//! `IntAna_IntQuadQuad` + `ProcessBounds` (ALine), then to walking.
//! Torus pairs (Cy/Co/Sp/To) use QuadQuadGeo circles; non-geometric
//! cylinder-cylinder uses `CyCyNoGeometric`.

use occt_core::gp::GpPnt2d;
use occt_geom::Surface;

use crate::fclass2d::FaceState;
use crate::geom_int::TopolTool;
use crate::geom_int::GeomIntLine;
use crate::int_tools_wline::PatchPoint;

#[path = "intpatch_impimp_quad.rs"]
mod quad;
#[path = "intpatch_impimp_glines.rs"]
mod glines;
#[path = "intpatch_impimp_pairs.rs"]
mod pairs;
#[path = "intpatch_impimp_sonb.rs"]
mod sonb;
#[path = "intpatch_impimp_put.rs"]
mod put;
#[path = "intpatch_impimp_rline.rs"]
mod rline;
#[path = "intpatch_impimp_cycy.rs"]
mod cycy;
#[path = "intpatch_impimp_bounds.rs"]
mod bounds;

use glines::PairOutcome;
use pairs::intersect_pair;
use put::{compute_vertex_parameters, put_points_on_line};
use quad::set_quad;
pub(crate) use sonb::{search_on_bounds, PathPoint};
pub(crate) use rline::{process_rline, process_segments};

pub(crate) use quad::{distance, gradient, quadric_tolerance, set_quad as try_set_quad, ImplicitQuad};

/// `IntPatch_ImpImpIntersection::IntStatus`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntStatus {
    Ok,
    InfiniteSectionCurve,
    Fail,
}

/// Analytic intersection of two implicit quadrics.
#[derive(Clone)]
pub struct ImpImpIntersection {
    status: IntStatus,
    empty: bool,
    tangent_faces: bool,
    opposite: bool,
    slin: Vec<GeomIntLine>,
    spnt: Vec<PatchPoint>,
}

impl ImpImpIntersection {
    pub fn new() -> Self {
        Self {
            status: IntStatus::Fail,
            empty: true,
            tangent_faces: false,
            opposite: false,
            slin: Vec::new(),
            spnt: Vec::new(),
        }
    }

    /// `Perform(S1, D1, S2, D2, TolArc, TolTang)`.
    pub fn perform(
        &mut self,
        s1: &dyn Surface,
        d1: &TopolTool,
        s2: &dyn Surface,
        d2: &TopolTool,
        tol_arc: f64,
        tol_tang: f64,
    ) {
        self.status = IntStatus::Fail;
        self.empty = true;
        self.tangent_faces = false;
        self.opposite = false;
        self.slin.clear();
        self.spnt.clear();

        let Some(q1) = set_quad(s1) else {
            return;
        };
        let Some(q2) = set_quad(s2) else {
            return;
        };

        let uv1 = cycy::CylUv {
            u0: d1.u_bounds().0,
            u1: d1.u_bounds().1,
            v0: d1.v_bounds().0,
            v1: d1.v_bounds().1,
        };
        let uv2 = cycy::CylUv {
            u0: d2.u_bounds().0,
            u1: d2.u_bounds().1,
            v0: d2.v_bounds().0,
            v1: d2.v_bounds().1,
        };
        match intersect_pair(&q1, &q2, tol_tang, uv1, uv2) {
            PairOutcome::Fail => {}
            PairOutcome::Infinite => {
                self.status = IntStatus::InfiniteSectionCurve;
                self.empty = true;
            }
            PairOutcome::Empty => {
                self.status = IntStatus::Ok;
                self.empty = true;
            }
            PairOutcome::Same => {
                self.status = IntStatus::Ok;
                self.empty = false;
                self.tangent_faces = true;
                if let Some(p) = sample_ref(&q1) {
                    let n1 = quad::normale(&q1, &p);
                    let n2 = quad::normale(&q2, &p);
                    self.opposite = n1.dot(&n2) < 0.0;
                }
            }
            PairOutcome::Result { lines, points } => {
                self.status = IntStatus::Ok;
                self.slin = lines;
                self.spnt = points;
                self.empty = self.slin.is_empty() && self.spnt.is_empty();
                if !self.empty {
                    self.postprocess(s1, d1, s2, d2, &q1, &q2, tol_arc);
                }
            }
        }
    }

    fn postprocess(
        &mut self,
        s1: &dyn Surface,
        d1: &TopolTool,
        s2: &dyn Surface,
        d2: &TopolTool,
        q1: &quad::ImplicitQuad,
        q2: &quad::ImplicitQuad,
        tol_arc: f64,
    ) {
        let sol1 = search_on_bounds(s1, q2, d1, tol_arc, tol_arc);
        let sol2 = search_on_bounds(s2, q1, d2, tol_arc, tol_arc);
        if !sol1.done || !sol2.done {
            self.status = IntStatus::Fail;
            return;
        }

        let mut all1 = sol1.all_arc_solution && q1.code() == q2.code();
        let mut all2 = sol2.all_arc_solution && q1.code() == q2.code();
        let nosolon_s1 = sol1.points.is_empty() && sol1.segments.is_empty();
        let nosolon_s2 = sol2.points.is_empty() && sol2.segments.is_empty();
        if nosolon_s1 && all1 {
            all1 = false;
        }
        if nosolon_s2 && all2 {
            all2 = false;
        }
        if all1 && all2 {
            self.empty = false;
            self.tangent_faces = true;
            self.slin.clear();
            self.spnt.clear();
            if let Some(p) = sample_ref(q1) {
                let n1 = quad::normale(q1, &p);
                let n2 = quad::normale(q2, &p);
                self.opposite = n1.dot(&n2) < 0.0;
            }
            self.status = IntStatus::Ok;
            return;
        }

        if !nosolon_s1 || !nosolon_s2 {
            self.empty = false;
            put_points_on_line(s1, s2, &sol1.points, &mut self.slin, true, q1, q2, tol_arc);
            put_points_on_line(s1, s2, &sol2.points, &mut self.slin, false, q2, q1, tol_arc);
            if !sol1.segments.is_empty() {
                process_segments(&sol1.segments, &mut self.slin, s1, s2, true, tol_arc);
            }
            if !sol2.segments.is_empty() {
                process_segments(&sol2.segments, &mut self.slin, s1, s2, false, tol_arc);
            }
            if !sol1.segments.is_empty() || !sol2.segments.is_empty() {
                process_rline(&mut self.slin, s1, s2, tol_arc, false);
            }
        } else {
            self.empty = self.slin.is_empty() && self.spnt.is_empty();
        }

        self.spnt.retain(|ip| {
            d1.classify(GpPnt2d::new(ip.u1, ip.v1), tol_arc) != FaceState::Out
                && d2.classify(GpPnt2d::new(ip.u2, ip.v2), tol_arc) != FaceState::Out
        });
        compute_vertex_parameters(&mut self.slin, tol_arc);
        self.empty = self.slin.is_empty() && self.spnt.is_empty();
    }

    /// `IsDone` — true unless the algorithm failed.
    pub fn is_done(&self) -> bool {
        self.status != IntStatus::Fail
    }

    pub fn status(&self) -> IntStatus {
        self.status
    }

    pub fn is_empty(&self) -> bool {
        self.empty
    }

    pub fn tangent_faces(&self) -> bool {
        self.tangent_faces
    }

    pub fn opposite_faces(&self) -> bool {
        self.opposite
    }

    pub fn nb_pnts(&self) -> i32 {
        self.spnt.len() as i32
    }

    pub fn point(&self, index: i32) -> Option<&PatchPoint> {
        self.spnt.get((index as usize).saturating_sub(1))
    }

    pub fn nb_lines(&self) -> i32 {
        self.slin.len() as i32
    }

    pub fn line(&self, index: i32) -> Option<&GeomIntLine> {
        self.slin.get((index as usize).saturating_sub(1))
    }

    pub fn lines(&self) -> &[GeomIntLine] {
        &self.slin
    }
}

impl Default for ImpImpIntersection {
    fn default() -> Self {
        Self::new()
    }
}

fn sample_ref(q: &quad::ImplicitQuad) -> Option<occt_core::gp::GpPnt> {
    Some(match q {
        quad::ImplicitQuad::Plane(p) => p.location(),
        quad::ImplicitQuad::Cylinder(c) => c.location(),
        quad::ImplicitQuad::Sphere(s) => s.location(),
        quad::ImplicitQuad::Cone(c) => c.apex(),
        quad::ImplicitQuad::Torus(t) => t.location(),
    })
}
