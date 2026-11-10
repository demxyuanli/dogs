//! `HLRBRep_Intersector` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_Intersector.hxx:35-97`,
//! `.cxx:57-753`).
//!
//! "The Intersector computes 2D intersections of the projections of 3D curves.
//! It can also compute the intersection of a 3D line and a surface."
//! (`hxx:34-35`). Three results are possible and `myTypePerform` says which one
//! the accessors must answer (`cxx:653-745`):
//!
//! * `0` - a single simulated point stored in `mySinglePoint` (`cxx:482-507`),
//! * `1` - the 2D edge/edge result of `myIntersector` (`cxx:93-478`),
//! * `2` - the line/surface result of `myCSIntersector` (`cxx:514-648`).
//!
//! ## Ported
//!
//! * `Perform(theEdge1, theDa1, theDb1)` (`cxx:93-131`) - edge auto
//!   intersection,
//! * `Perform(theNA, ..., theNB, ..., theNoBound)` (`cxx:135-478`) - edge/edge
//!   with the `aPasBon` decalage loop that re-shrinks a domain whose endpoints
//!   coincide within `1000 * tol`,
//! * `SimulateOnePoint` (`cxx:482-507`),
//! * `Load` (`cxx:510-515`), `Perform(L, P)` (`cxx:518-648`), `CSPoint`
//!   (`cxx:706-710`) and `CSSegment` (`cxx:738-742`) - the
//!   `myCSIntersector` = `HLRBRep_InterCSurf` half (`hxx:93`), see
//!   [`crate::hlrbrep_inter_csurf`],
//! * `IsDone` / `NbPoints` / `Point` / `NbSegments` / `Segment` (`cxx:653-737`)
//!   for all three types.
//!
//! The 2D work is `myIntersector` = `HLRBRep_CInter` (`hxx:92`), whose
//! constructor sets `SetMinNbSamples(4)` (`cxx:84-87`); see
//! [`crate::hlrbrep_c_inter`] and [`crate::hlrbrep_curve_tool`].
//!
//! ## Deviations
//!
//! * `mySurface` (`hxx:80`) and `myPolyhedron` (`hxx:81`) are a raw pointer and
//!   a `new`/`delete`-managed pointer in OCCT. The port owns a
//!   [`HlrBRepSurface`] clone taken by [`HlrBRepIntersector::load`] and an
//!   `Option<ThePolyhedronOfHInter>` that `load` drops, so the caching
//!   (`cxx:534-546` builds the polyhedron once per loaded surface) is
//!   preserved. `Destroy()` (`cxx:748-752`), which only deletes `myPolyhedron`,
//!   therefore has no counterpart.
//! * `myTypePerform == 43` (`cxx:671-674`) is a value no assignment in the
//!   class produces; the port keeps the test verbatim.

use occt_core::gp::{GpLin, GpPnt2d, GpVec2d};
use occt_core::intimpargen::gen::determine_transition_simple;
use occt_core::intres2d::{
    IntRes2dDomain, IntRes2dIntersectionPoint, IntRes2dIntersectionSegment, IntRes2dPosition,
};
use occt_core::precision::CONFUSION;

use crate::hlrbrep_c_inter::HlrBRepCInter;
use crate::hlrbrep_edge_data::HlrBRepEdgeData;
use crate::hlrbrep_inter_csurf::HlrBRepInterCSurf;
use crate::hlrbrep_surface::HlrBRepSurface;
use occt_topo::intcurvesurface::{
    IntersectionPoint, IntersectionSegment, ThePolygonOfHInter, ThePolyhedronOfHInter,
};

/// `HLRBRep_Intersector` (`HLRBRep_Intersector.hxx:35-97`).
#[derive(Clone)]
pub struct HlrBRepIntersector {
    /// `IntRes2d_IntersectionPoint mySinglePoint` (`hxx:91`).
    my_single_point: IntRes2dIntersectionPoint,
    /// `int myTypePerform` (`hxx:92`): `0` simulated point, `1` curve/curve,
    /// `2` curve/surface.
    my_type_perform: i32,
    /// `HLRBRep_CInter myIntersector` (`hxx:93`).
    my_intersector: HlrBRepCInter,
    /// `HLRBRep_InterCSurf myCSIntersector` (`hxx:94`).
    my_cs_intersector: HlrBRepInterCSurf,
    /// `HLRBRep_Surface* mySurface` (`hxx:80`): owned clone, see the module
    /// notes.
    my_surface: Option<HlrBRepSurface>,
    /// `HLRBRep_ThePolyhedronOfInterCSurf* myPolyhedron` (`hxx:81`): built once
    /// per loaded surface (`cxx:534-546`), dropped by `Load`.
    my_polyhedron: Option<ThePolyhedronOfHInter>,
}

impl Default for HlrBRepIntersector {
    fn default() -> Self {
        Self::new()
    }
}

impl HlrBRepIntersector {
    /// `HLRBRep_Intersector()` (`cxx:57-88`): `myTypePerform(0)`,
    /// `mySurface(nullptr)`, `myPolyhedron(nullptr)`, then
    /// `myIntersector.SetMinNbSamples(4)` (`cxx:84-87`).
    pub fn new() -> Self {
        let mut my_intersector = HlrBRepCInter::new();
        my_intersector.set_min_nb_samples(4);
        Self {
            my_single_point: IntRes2dIntersectionPoint::new(),
            my_type_perform: 0,
            my_intersector,
            my_cs_intersector: HlrBRepInterCSurf::new(),
            my_surface: None,
            my_polyhedron: None,
        }
    }

    /// `Perform(theEdge1, theDa1, theDb1)` (`cxx:93-131`): auto intersection of
    /// one edge, its domain cut by `theDa1 * (b - a)` at the start and
    /// `theDb1 * (b - a)` at the end.
    pub fn perform_auto(&mut self, the_edge1: &HlrBRepEdgeData, the_da1: f64, the_db1: f64) {
        self.my_type_perform = 1;

        let (mut a, ta, mut b, tb) = the_edge1.status().bounds();
        let d = b - a;
        if the_da1 != 0.0 {
            a += d * the_da1;
        }
        if the_db1 != 0.0 {
            b -= d * the_db1;
        }
        let c1 = the_edge1.geometry();
        let pa = c1.d0(a);
        let pb = c1.d0(b);
        let a2 = c1.parameter2d(a);
        let b2 = c1.parameter2d(b);
        let d1 = IntRes2dDomain::bounded(&pa, a2, ta as f64, &pb, b2, tb as f64);

        // `tol = Precision::Confusion()` (`cxx:126-129`: the edge tolerance was
        // replaced by `Precision::Confusion()` for OCC27341).
        let tol = CONFUSION;
        self.my_intersector.perform_curve(c1, &d1, tol, tol);
    }

    /// `Perform(theNA, theEdge1, theDa1, theDb1, theNB, theEdge2, theDa2,
    /// theDb2, theNoBound)` (`cxx:135-478`): the two-edge intersection.
    ///
    /// `theNA` / `theNB` are only read by the `PERF` counters in OCCT
    /// (`cxx:136-137`, `:138-139`), which the port drops.
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        the_edge1: &HlrBRepEdgeData,
        the_da1: f64,
        the_db1: f64,
        the_edge2: &HlrBRepEdgeData,
        the_da2: f64,
        the_db2: f64,
        the_en_bout: bool,
    ) {
        let c1 = the_edge1.geometry();
        let c2 = the_edge2.geometry();

        self.my_type_perform = 1;

        // `tol1 = tol2 = Precision::Confusion(); tol = max(tol1, tol2)`
        // (`cxx:160-175`).
        let tol = CONFUSION;

        let mut decalage_a1 = 100.0;
        let mut decalage_a2 = 100.0;
        let mut decalage_b1 = 100.0;
        let mut decalage_b2 = 100.0;

        let mut pas_bon;
        let mut d1;
        let mut d2;
        loop {
            pas_bon = false;

            // ---- edge 1 ---------------------------------------------------
            let (mut a1, mut ta, mut b1, mut tb) = the_edge1.status().bounds();
            // `pdist = max(tol, 0.0000001)` (`cxx:194-199`).
            let pdist = if tol < 0.000_000_1 { 0.000_000_1 } else { tol };
            let d = b1 - a1;

            if the_da1 != 0.0 {
                let (_, va1) = c1.d1(a1);
                let qwe = va1.magnitude();
                if qwe > 1e-12 {
                    let dd = pdist * decalage_a1 / qwe;
                    if dd < d * 0.4 {
                        a1 += dd;
                    } else {
                        a1 += d * the_da1;
                        decalage_a1 = -1.0;
                    }
                } else {
                    a1 += d * the_da1;
                    decalage_a1 = -1.0;
                }
            }

            if the_db1 != 0.0 {
                let (_, vb1) = c1.d1(b1);
                let qwe = vb1.magnitude();
                if qwe > 1e-12 {
                    let dd = pdist * decalage_b1 / qwe;
                    if dd < d * 0.4 {
                        b1 -= dd;
                    } else {
                        b1 -= d * the_db1;
                        decalage_b1 = -1.0;
                    }
                } else {
                    b1 -= d * the_db1;
                    decalage_b1 = -1.0;
                }
            }

            let pa1 = c1.d0(a1);
            let pb1 = c1.d0(b1);
            let a1 = c1.parameter2d(a1);
            let b1 = c1.parameter2d(b1);

            if the_en_bout {
                ta = -1.0;
                tb = -1.0;
            }
            if ta as f64 > tol {
                ta = tol as f32;
            }
            if tb as f64 > tol {
                tb = tol as f32;
            }
            d1 = IntRes2dDomain::bounded(&pa1, a1, ta as f64, &pb1, b1, tb as f64);

            // ---- edge 2 ---------------------------------------------------
            let (mut a2, mut ta, mut b2, mut tb) = the_edge2.status().bounds();
            let d = b2 - a2;

            if the_da2 != 0.0 {
                let (_, va2) = c2.d1(a2);
                let qwe = va2.magnitude();
                if qwe > 1e-12 {
                    let dd = pdist * decalage_a2 / qwe;
                    if dd < d * 0.4 {
                        a2 += dd;
                    } else {
                        a2 += d * the_da2;
                        decalage_a2 = -1.0;
                    }
                } else {
                    a2 += d * the_da2;
                    decalage_a2 = -1.0;
                }
            }

            if the_db2 != 0.0 {
                let (_, vb2) = c2.d1(b2);
                let qwe = vb2.magnitude();
                if qwe > 1e-12 {
                    let dd = pdist * decalage_b2 / qwe;
                    if dd < d * 0.4 {
                        b2 -= dd;
                    } else {
                        b2 -= d * the_db2;
                        decalage_b2 = -1.0;
                    }
                } else {
                    b2 -= d * the_db2;
                    decalage_b2 = -1.0;
                }
            }

            let pa2 = c2.d0(a2);
            let pb2 = c2.d0(b2);
            let a2 = c2.parameter2d(a2);
            let b2 = c2.parameter2d(b2);

            if the_en_bout {
                ta = -1.0;
                tb = -1.0;
            }
            if ta as f64 > tol {
                ta = tol as f32;
            }
            if tb as f64 > tol {
                tb = tol as f32;
            }
            d2 = IntRes2dDomain::bounded(&pa2, a2, ta as f64, &pb2, b2, tb as f64);

            // `cxx:370-431`: when the caller forbade the extremities and the
            // nearest pair of domain corners is closer than `1000 * tol`, the
            // domains are shrunk further instead of being intersected.
            if the_en_bout {
                let mut a1a2 = f64::MAX;
                let mut a1b2 = f64::MAX;
                let mut b1a2 = f64::MAX;
                let mut b1b2 = f64::MAX;
                if the_da1 != 0.0 || the_da2 != 0.0 {
                    a1a2 = pa1.distance(&pa2);
                }
                if the_da1 != 0.0 || the_db2 != 0.0 {
                    a1b2 = pa1.distance(&pb2);
                }
                if the_db1 != 0.0 || the_da2 != 0.0 {
                    b1a2 = pb1.distance(&pa2);
                }
                if the_db1 != 0.0 || the_db2 != 0.0 {
                    b1b2 = pb1.distance(&pb2);
                }

                let mut cote = 1;
                let mut mindist = a1a2;
                if mindist > a1b2 {
                    mindist = a1b2;
                    cote = 2;
                }
                if mindist > b1a2 {
                    mindist = b1a2;
                    cote = 3;
                }
                if mindist > b1b2 {
                    mindist = b1b2;
                    cote = 4;
                }

                if mindist < tol * 1000.0 {
                    pas_bon = true;
                    match cote {
                        1 => {
                            decalage_a1 *= 2.0;
                            decalage_a2 *= 2.0;
                        }
                        2 => {
                            decalage_a1 *= 2.0;
                            decalage_b2 *= 2.0;
                        }
                        3 => {
                            decalage_b1 *= 2.0;
                            decalage_a2 *= 2.0;
                        }
                        _ => {
                            decalage_b1 *= 2.0;
                            decalage_b2 *= 2.0;
                        }
                    }
                    if decalage_a1 < 0.0
                        || decalage_a2 < 0.0
                        || decalage_b1 < 0.0
                        || decalage_b2 <= 0.0
                    {
                        pas_bon = false;
                    }
                }
            }

            if !pas_bon {
                self.my_intersector.perform(c1, &d1, c2, &d2, tol, tol);
            }
            if !pas_bon {
                break;
            }
        }
    }

    /// `SimulateOnePoint(theEdge1, theU, theEdge2, theV)` (`cxx:482-507`):
    /// synthesize the middle point of the two curves at `(theU, theV)`.
    pub fn simulate_one_point(
        &mut self,
        the_edge1: &HlrBRepEdgeData,
        the_u: f64,
        the_edge2: &HlrBRepEdgeData,
        the_v: f64,
    ) {
        let c1 = the_edge1.geometry();
        let c2 = the_edge2.geometry();
        let u3 = c1.parameter3d(the_u);
        let v3 = c2.parameter3d(the_v);
        let (p13, t13): (GpPnt2d, GpVec2d) = c1.d1(u3);
        let (_p23, t23): (GpPnt2d, GpVec2d) = c2.d1(v3);

        // `Pos1 = Pos2 = IntRes2d_Middle` (`cxx:500-501`) then
        // `IntImpParGen::DetermineTransition(Pos1, T13, Tr1, Pos2, T23, Tr2, 0.)`
        // (`cxx:503`). The return value is discarded in OCCT too.
        let pos1 = IntRes2dPosition::Middle;
        let pos2 = IntRes2dPosition::Middle;
        let mut tr1 = occt_core::intres2d::IntRes2dTransition::new();
        let mut tr2 = occt_core::intres2d::IntRes2dTransition::new();
        let _ = determine_transition_simple(pos1, &t13, &mut tr1, pos2, &t23, &mut tr2, 0.0);

        self.my_type_perform = 0;
        // `mySinglePoint.SetValues(P13, theU, theV, Tr1, Tr2, false)`
        // (`cxx:505`).
        self.my_single_point
            .set_values(&p13, the_u, the_v, &tr1, &tr2, false);
    }

    /// `Load(theSurface)` (`cxx:510-515`): remember the surface and drop the
    /// cached polyhedron, which `Perform(L, P)` rebuilds on demand.
    pub fn load(&mut self, the_surface: &HlrBRepSurface) {
        self.my_surface = Some(the_surface.clone());
        self.my_polyhedron = None;
    }

    /// `Perform(L, P)` (`cxx:518-648`): line/surface intersection.
    ///
    /// A plane / cylinder / cone / sphere / torus goes to
    /// `myCSIntersector.Perform(L, mySurface)` (`cxx:524-531`), i.e. the
    /// analytic arms of `HLRBRep_InterCSurf::PerformConicSurfLine`. Every other
    /// type samples the surface once per loaded surface
    /// (`cxx:534-546`), clips `[pmin, pmax]` - the parameter range of the
    /// polyhedron's eight bounding-box corners on `L` - against the caller's
    /// `P` (`cxx:547-631`), and runs the polygon/polyhedron overload
    /// (`cxx:632-637`).
    pub fn perform_line(&mut self, the_l: &GpLin, the_p: f64) {
        self.my_type_perform = 2;
        let Some(surface) = self.my_surface.as_ref() else {
            // `mySurface` is a raw pointer in OCCT; `HLRBRep_SurfaceTool` would
            // dereference null here.
            return;
        };

        match crate::hlrbrep_surface_tool::get_type(surface) {
            occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType::Plane
            | occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType::Cylinder
            | occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType::Cone
            | occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType::Sphere
            | occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType::Torus => {
                self.my_cs_intersector.perform(the_l, surface);
            }
            _ => {
                if self.my_polyhedron.is_none() {
                    let u1 = crate::hlrbrep_surface_tool::first_u_parameter(surface);
                    let v1 = crate::hlrbrep_surface_tool::first_v_parameter(surface);
                    let u2 = crate::hlrbrep_surface_tool::last_u_parameter(surface);
                    let v2 = crate::hlrbrep_surface_tool::last_v_parameter(surface);
                    let nbsu = crate::hlrbrep_surface_tool::nb_samples_u_range(surface, u1, u2);
                    let nbsv = crate::hlrbrep_surface_tool::nb_samples_v_range(surface, v1, v2);
                    let Some(view) = crate::hlrbrep_surface::HlrBRepSurfaceView::new(surface)
                    else {
                        return;
                    };
                    self.my_polyhedron = Some(ThePolyhedronOfHInter::new(
                        &view,
                        nbsu.max(1) as usize,
                        nbsv.max(1) as usize,
                        u1,
                        v1,
                        u2,
                        v2,
                    ));
                }
                let Some(polyhedron) = self.my_polyhedron.as_ref() else {
                    return;
                };
                let Some((x0, y0, z0, x1, y1, z1)) = polyhedron.bounding().get() else {
                    return;
                };

                // `cxx:547-618`: parameter of all eight box corners on the line.
                let mut p = occt_core::elib::clib::line_parameter(the_l, &occt_core::gp::GpPnt::new(x0, y0, z0));
                let (mut pmin, mut pmax) = (p, p);
                for corner in [
                    occt_core::gp::GpPnt::new(x0, y0, z1),
                    occt_core::gp::GpPnt::new(x1, y0, z0),
                    occt_core::gp::GpPnt::new(x1, y0, z1),
                    occt_core::gp::GpPnt::new(x0, y1, z0),
                    occt_core::gp::GpPnt::new(x0, y1, z1),
                    occt_core::gp::GpPnt::new(x1, y1, z0),
                    occt_core::gp::GpPnt::new(x1, y1, z1),
                ] {
                    p = occt_core::elib::clib::line_parameter(the_l, &corner);
                    if pmin > p {
                        pmin = p;
                    }
                    if pmax < p {
                        pmax = p;
                    }
                }
                pmin -= 0.000_001;
                pmax += 0.000_001;

                if pmin > the_p {
                    // `cxx:621-624`: push the box past `P` so the interference
                    // rejects every point.
                    pmin = pmax + 1.0;
                    pmax += 2.0;
                } else if pmax > the_p {
                    // `cxx:626-630`.
                    pmax = the_p + 0.000_000_1;
                }

                let polygon = ThePolygonOfHInter::with_range(
                    &crate::hlrbrep_inter_csurf::HlrBRepLineCurve::new(*the_l),
                    pmin,
                    pmax,
                    3,
                );
                self.my_cs_intersector
                    .perform_polygon_polyhedron(the_l, &polygon, surface, polyhedron);
            }
        }
    }

    /// `CSPoint(N)` (`cxx:706-710`): the curve/surface point of
    /// `myCSIntersector`.
    pub fn cs_point(&self, n: usize) -> &IntersectionPoint {
        self.my_cs_intersector.point(n)
    }

    /// `CSSegment(N)` (`cxx:738-742`).
    pub fn cs_segment(&self, n: usize) -> &IntersectionSegment {
        self.my_cs_intersector.segment(n)
    }

    /// `IsDone()` (`cxx:653-665`).
    pub fn is_done(&self) -> bool {
        match self.my_type_perform {
            1 => self.my_intersector.is_done(),
            // `cxx:659-662`: `myCSIntersector.IsDone()`.
            2 => self.my_cs_intersector.is_done(),
            // `cxx:661-664`: `return true;`.
            _ => true,
        }
    }

    /// `NbPoints()` (`cxx:669-688`).
    pub fn nb_points(&self) -> usize {
        if self.my_type_perform == 43 {
            // `cxx:671-674`.
            return 0;
        }
        match self.my_type_perform {
            1 => self.my_intersector.nb_points(),
            // `cxx:680-683`.
            2 => self.my_cs_intersector.nb_points(),
            // `cxx:684-687`: `return 1;`.
            _ => 1,
        }
    }

    /// `Point(N)` (`cxx:692-702`).
    pub fn point(&self, n: usize) -> IntRes2dIntersectionPoint {
        if self.my_type_perform == 0 {
            // `cxx:694-697`: `return mySinglePoint;`.
            self.my_single_point
        } else {
            // `cxx:698-701`: `myIntersector.Point(N)`.
            self.my_intersector.point(n)
        }
    }

    /// `NbSegments()` (`cxx:715-727`).
    pub fn nb_segments(&self) -> usize {
        match self.my_type_perform {
            1 => self.my_intersector.nb_segments(),
            // `cxx:720-723`.
            2 => self.my_cs_intersector.nb_segments(),
            // `cxx:724-727`: `return 0;`.
            _ => 0,
        }
    }

    /// `Segment(N)` (`cxx:731-735`).
    pub fn segment(&self, n: usize) -> IntRes2dIntersectionSegment {
        self.my_intersector.segment(n)
    }
}
