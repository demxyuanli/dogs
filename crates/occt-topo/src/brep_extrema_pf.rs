//! `BRepExtrema_ExtPF` — extrema between a vertex and a *bounded* face.
//!
//! Source: `TKTopAlgo/BRepExtrema/BRepExtrema_ExtPF.hxx:30-90`,
//! `BRepExtrema_ExtPF.cxx:20-100` (OCCT 8.0.0).
//!
//! `Initialize` builds the `Extrema_ExtPS` window over
//! `BRepTools::UVBounds(face)`; `Perform` keeps the surface extrema whose
//! `(u, v)` the face classifier reports `ON` or `IN`. The vertex of the C++
//! call is reduced to its point (`BRep_Tool::Pnt(TheVertex)`, `cxx:73`), which
//! is all the extrema engine reads.

use occt_core::gp::GpPnt;
use occt_geom::extrema_surf::{ext_ps_surface_type, ExtPs, ExtPsSurfaceType};

use crate::brep_tool::BRepTool;
use crate::brep_tools_uv;
use crate::fclass2d::{FaceState, FClass2d};
use crate::shape::Face;

/// `BRepExtrema_ExtPF` (`hxx:30-90`): the surviving solutions of a vertex/face
/// extrema search, in the order `Perform` appends them (`cxx:93-98`).
#[derive(Debug, Clone, Default)]
pub struct ExtPF {
    /// `mySqDist` (`cxx:70`).
    sq_dist: Vec<f64>,
    /// `myPoints` (`cxx:71`), as `(u, v, point)`.
    points: Vec<(f64, f64, GpPnt)>,
    /// `myExtPS.IsDone()`; `false` covers both the `OtherSurface` early return
    /// of `Initialize` (`cxx:50-53`) and a not-done `Perform`.
    done: bool,
}

impl ExtPF {
    /// `BRepExtrema_ExtPF(TheVertex, TheFace)` (`cxx:31-38`): `Initialize` then
    /// `Perform`.
    pub fn perform(p: &GpPnt, face: &Face) -> Self {
        let mut this = ExtPF {
            sq_dist: Vec::new(),
            points: Vec::new(),
            done: false,
        };

        // `Initialize(TheFace, ...)` (`cxx:42-65`).
        let Some(surf) = BRepTool::face_surface_world(face) else {
            return this;
        };
        if ext_ps_surface_type(surf.as_ref()) == ExtPsSurfaceType::Other {
            // `cxx:50-53`: protect against a non-geometric type.
            return this;
        }
        let tol = BRepTool::face_tolerance(face).min(occt_core::precision::CONFUSION);
        let (ur, vr) = surf
            .uv_resolution(tol)
            .unwrap_or((0.0, 0.0));
        let tol_u = ur.max(occt_core::precision::PCONFUSION);
        let tol_v = vr.max(occt_core::precision::PCONFUSION);
        let (u1, u2, v1, v2) = brep_tools_uv::uv_bounds_face(face);

        let mut ext_ps = ExtPs::new();
        ext_ps.initialize(surf.as_ref(), u1, u2, v1, v2, tol_u, tol_v);

        // `Perform(TheVertex, TheFace)` (`cxx:68-100`).
        ext_ps.perform(p);
        if ext_ps.is_done() {
            // `BRepClass_FaceClassifier classifier; classifier.Perform(F, uv, Tol)`
            // (`cxx:84-92`), with `Tol = BRep_Tool::Tolerance(TheFace)`.
            let tol_class = BRepTool::face_tolerance(face);
            let Ok(classifier) = FClass2d::new(face, tol_class) else {
                return this;
            };
            for i in 1..=ext_ps.nb_ext() {
                let (u, v, point) = ext_ps.point(i);
                let state = classifier.perform(occt_core::gp::GpPnt2d::new(u, v));
                if state == FaceState::On || state == FaceState::In {
                    this.sq_dist.push(ext_ps.square_distance(i));
                    this.points.push((u, v, point));
                }
            }
        }
        // `myExtPS.IsDone()` of `cxx:82` gates the append loop only; the answer
        // list is what `NbExt()`/`SquareDistance()` read.
        this.done = ext_ps.is_done();
        this
    }

    /// `IsDone()`: the search ran (the C++ `Perform` has no explicit failure
    /// flag beyond `myExtPS.IsDone()`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `NbExt()` (`hxx:60-63`).
    pub fn nb_ext(&self) -> usize {
        self.points.len()
    }

    /// `SquareDistance(i)` (`hxx:65-68`), 1-based like OCCT.
    pub fn square_distance(&self, i: usize) -> Option<f64> {
        self.sq_dist.get(i.checked_sub(1)?).copied()
    }

    /// `Parameter(i, U, V)` (`hxx:70-73`), 1-based like OCCT.
    pub fn parameter(&self, i: usize) -> Option<(f64, f64)> {
        self.points
            .get(i.checked_sub(1)?)
            .map(|(u, v, _)| (*u, *v))
    }
}
