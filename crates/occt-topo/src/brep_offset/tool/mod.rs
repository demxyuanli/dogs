//! `BRepOffset_Tool` (`BRepOffset_Tool.hxx` / `BRepOffset_Tool.cxx`).
//!
//! Ported here so far (control flow taken from
//! `src/ModelingAlgorithms/TKOffset/BRepOffset/BRepOffset_Tool.cxx`):
//!
//! - `EdgeVertices` (`:141-151`)
//! - `Gabarit` (`:313-323`)
//! - `OrientSection` (`:486-568`)
//! - `FindCommonShapes` (`:572-618`) - `common_shapes.rs`
//! - `MapVertexEdges` (`:3685-3716`) - `common_shapes.rs`
//! - `BuildNeighbour` (`:3720-3789`) - `common_shapes.rs`
//! - `CheckPlanesNormals` (`:4509-4533`)
//! - `CheckBounds` (`:3309-3433`) together with its file-local
//!   `CompactUVBounds` (`:3267-3305`)
//! - `EnLargeFace` (`:3437-3634`) - `enlarge_face.rs`, with
//!   `EnlargeGeometry` (`:2923-3222`) - `enlarge_geometry.rs` - and
//!   `MakeFace` (`:2580-2919`) / `UpdatePCurves` (`:3230-3263`) -
//!   `extended_face_builder.rs`
//! - `ExtentEdge` (`:2072-2157`), `ProjectVertexOnEdge` (`:2161-2249`),
//!   `SelectEdge` (`:2536-2576`) and `TryParameter` (`:3638-3681`) -
//!   `extent_edge.rs`
//! - `BuildPCurves` (`:327-483`) with `FindPeriod` (`:155-190`) and
//!   `PutInBounds` (`:197-...`) - `pcurve_builder.rs`
//! - `IsOnSurface` (`:642-743`) and `ToSmall` (`:622-638`) -
//!   `curve_checks.rs`
//! - `TryProject` (`:1945-2006`) - `try_project.rs`
//! - `PipeInter` (`:747-804`) - `pipe_inter.rs`
//! - `Inter2d` (`:2253-2532`) - `inter2d.rs`
//! - `UpdateVertexTolerances` (`:4613-4658`), `IsAutonomVertex` (`:811-916`,
//!   both overloads), `AreConnex` (`:918-928`), `AreClosed` (`:930-939`),
//!   `BSplineEdges` (`:942-987`), `AngleWireEdge` (`:991-1032`),
//!   `ReconstructPCurves` (`:1035-1057`), `ConcatPCurves` (`:1062-1160`),
//!   `Glue` (`:1164-1259`), `AssembleEdge` (`:1345-1433`) and
//!   `GeomProjLib::Curve2d` (`GeomProjLib.cxx:54-132`) -
//!   `inter3d_statics.rs`
//! - `CheckIntersFF` (`:1264-1343`) together with the file-local
//!   `BOPTools_AlgoTools::MakeConnexityBlocks(compound, VERTEX, EDGE, ...)`
//!   walk (`BOPTools_AlgoTools.cxx:105-183`) - `check_inters_ff.rs`
//! - `PerformPlanes` (`:4537-4602`) and `IsInf` (`:4606-4611`) -
//!   `perform_planes.rs`
//! - `Inter3D` (`:1441-1941`) - `inter3d.rs`
//! - `ExtentFace` (`:3793-4333`) - `extent_face.rs`
//! - `InterOrExtent` (`:2010-2068`) - `inter_or_extent.rs`
//! - `Deboucle3D` (`:4337-4418`) - `deboucle3d.rs`
//! - `CorrectOrientation` (`:4442-4505`) with its file-local `IsInOut`
//!   (`:4422-4438`) - `correct_orientation.rs`
//!
//! Nothing of `BRepOffset_Tool.cxx` is left unported.

mod check_inters_ff;
mod common_shapes;
mod correct_orientation;
mod curve_checks;
mod deboucle3d;
mod enlarge_face;
mod enlarge_geometry;
mod extended_face_builder;
mod extent_edge;
mod extent_face;
mod inter2d;
mod inter3d;
mod inter3d_statics;
mod inter_or_extent;
mod pcurve_builder;
mod perform_planes;
mod pipe_inter;
mod try_project;

pub use enlarge_face::EnLargeFaceOptions;
pub use extent_face::ShapeMap;

pub(crate) use check_inters_ff::make_connexity_blocks;

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpDir, GpDir2d, GpLin2d};
use occt_core::precision::{ANGULAR, CONFUSION, REAL_SMALL};
use occt_geom::geom_bnd_lib::box_curve;
use occt_geom::{Curve, GeomTrimmedCurveBasis};
use occt_geom2d::{Curve2d, Geom2dLine};

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d;
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_tool::BRepTool;
use crate::brep_offset::analyse::Analyse;
use crate::chfi_ds::TypeOfConcavity;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::shape::{Edge, Face, Vertex};
use crate::shape_custom_curve2d::ShapeCustomCurve2d;
use crate::topo_tools::{edge_vertices, explorer_shapes};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// Static helpers of `BRepOffset_Tool`.
pub struct BRepOffsetTool;

impl BRepOffsetTool {
    /// `BRepOffset_Tool::EdgeVertices(E, V1, V2)` (`cxx:141-151`):
    /// `TopExp::Vertices(E, V1, V2)`, swapped when the edge is `REVERSED`, so
    /// `V1` is always the vertex the edge starts at.
    ///
    /// `(None, None)` is the pair of null vertices of an edge without endpoint
    /// children (`TopExp::Vertices` fills neither).
    pub fn edge_vertices(e: &Edge) -> (Option<Vertex>, Option<Vertex>) {
        let (v1, v2) = edge_vertices(e);
        if e.0.orientation() == Orientation::Reversed {
            (v2, v1)
        } else {
            (v1, v2)
        }
    }

    /// `BRepOffset_Tool::Gabarit(aCurve)` (`cxx:313-323`): the largest extent
    /// of the curve's bounding box.
    ///
    /// `BndLib_Add3dCurve::Add(GeomAdaptor_Curve(aCurve),
    /// Precision::Confusion(), aBox)` is `occt_geom::geom_bnd_lib::box_curve`
    /// over the curve's own parameter range. A void box makes `Bnd_Box::Get`
    /// raise `Standard_ConstructionError` in OCCT; the port answers `0.0`.
    pub fn gabarit(curve: &dyn Curve) -> f64 {
        let (u1, u2) = (curve.first_parameter(), curve.last_parameter());
        let b = box_curve(curve, u1, u2, CONFUSION);
        match b.get() {
            Some((xmin, xmax, ymin, ymax, zmin, zmax)) => {
                let mut dist = (xmax - xmin).max(ymax - ymin);
                dist = dist.max(zmax - zmin);
                dist
            }
            // Port guard: OCCT raises on a void box.
            None => 0.0,
        }
    }

    /// `BRepOffset_Tool::OrientSection(E, F1, F2, O1, O2)` (`cxx:486-568`):
    /// the orientation of the section edge `E` on `F1` and on `F2`.
    ///
    /// `None` groups OCCT's four unguarded null handles
    /// (`BRep_Tool::Surface(F1/F2)`, `CurveOnSurface(E, F1/F2)`,
    /// `Curve(E, L, f, l)`, each dereferenced at `:527`, `:534`, `:516`): the
    /// port answers `None` instead of faulting.
    pub fn orient_section(
        e: &Edge,
        f1: &Face,
        f2: &Face,
    ) -> Option<(Orientation, Orientation)> {
        let s1 = BRepTool::face_surface(f1)?;
        let s2 = BRepTool::face_surface(f2)?;
        let (c1, _f1, _l1) = boptools_2d::curve_on_surface_range(e, f1)?;
        let (c2, _f2, _l2) = boptools_2d::curve_on_surface_range(e, f2)?;
        // `BRep_Tool::Curve(E, L, f, l)` with `L` applied: `:516` transforms
        // the first derivative by `L.Transformation()`, i.e. world coordinates.
        let c = BRepTool::edge_curve_world(e)?;
        let (f, l) = BRepTool::edge_parameters(e);

        // `:501-514`: the abscissa midpoint `GCPnts_AbscissaPoint(BAcurve,
        // Length(BAcurve) / 2, f)`, else `BOPTools_AlgoTools2D::IntermediatePoint`.
        let par_on_c = Self::abscissa_midpoint(&c, f, l);

        // `:516-520`: `T1 = C->DN(ParOnC, 1).Transformed(L.Transformation())`.
        let mut t1 = c.eval_dn(par_on_c, 1);
        if t1.square_magnitude() > REAL_SMALL {
            t1.normalize();
        }

        // `:522-531`: the face normal of `F1` at the pcurve point.
        let p = c1.d0(par_on_c);
        let (_, d1u, d1v) = s1.d1(p.x(), p.y());
        let mut dn1 = d1u.crossed(&d1v);
        if f1.0.orientation() == Orientation::Reversed {
            dn1.reverse();
        }

        // `:532-539`.
        let p = c2.d0(par_on_c);
        let (_, d1u, d1v) = s2.d1(p.x(), p.y());
        let mut dn2 = d1u.crossed(&d1v);
        if f2.0.orientation() == Orientation::Reversed {
            dn2.reverse();
        }

        // `:540-549`.
        let pro_vec = dn2.crossed(&t1);
        let prod = dn1.dot(&pro_vec);
        let mut o1 = if prod < 0.0 {
            Orientation::Forward
        } else {
            Orientation::Reversed
        };

        // `:550-558`.
        let pro_vec = dn1.crossed(&t1);
        let prod = dn2.dot(&pro_vec);
        let mut o2 = if prod < 0.0 {
            Orientation::Forward
        } else {
            Orientation::Reversed
        };

        // `:560-567`.
        if f1.0.orientation() == Orientation::Reversed {
            o1 = o1.reversed();
        }
        if f2.0.orientation() == Orientation::Reversed {
            o2 = o2.reversed();
        }
        Some((o1, o2))
    }

    /// `GCPnts_AbscissaPoint(BAcurve, GCPnts_AbscissaPoint::Length(BAcurve) /
    /// 2.0, f)` with the `!IsDone()` fallback of `cxx:514`
    /// (`BOPTools_AlgoTools2D::IntermediatePoint(f, l)`).
    ///
    /// `BAcurve` is a `BRepAdaptor_Curve`, i.e. the edge's curve *bounded* by
    /// the edge range `[f, l]`; a stored `Geom_Line` / `Geom_Circle` has an
    /// infinite own range, so the bracket is rebuilt as a
    /// [`GeomTrimmedCurveBasis`] view (same parameterisation) exactly like
    /// `GeomAdaptor_Curve::Load(C, First, Last)`.
    fn abscissa_midpoint(c: &std::sync::Arc<dyn Curve>, f: f64, l: f64) -> f64 {
        let view: std::sync::Arc<dyn Curve> =
            if c.first_parameter() == f && c.last_parameter() == l {
                c.clone()
            } else {
                std::sync::Arc::new(GeomTrimmedCurveBasis::new(c.clone(), f, l))
            };
        let length = occt_geom::gcpnts::curve_length_between(view.as_ref(), f, l, None);
        occt_geom::gcpnts::abscissa_point(view.as_ref(), length / 2.0, f)
            .unwrap_or_else(|_| boptools_2d::intermediate_point(f, l))
    }

    /// `BRepOffset_Tool::CheckPlanesNormals(theFace1, theFace2, theTolAng)`
    /// (`cxx:4509-4533`): both faces planar and their oriented normals within
    /// `theTolAng`.
    pub fn check_planes_normals(f1: &Face, f2: &Face, tol_ang: f64) -> bool {
        let a1 = SurfaceAdaptor::of_face(f1);
        let a2 = SurfaceAdaptor::of_face(f2);
        if a1.get_type() != GeomAbsSurfaceType::Plane || a2.get_type() != GeomAbsSurfaceType::Plane
        {
            return false;
        }
        let (Some(p1), Some(p2)) = (a1.plane(), a2.plane()) else {
            // Port guard: `GeomAdaptor_Surface::Plane()` raises on a non-plane.
            return false;
        };
        let mut dn1: GpDir = p1.position().direction();
        if f1.0.orientation() == Orientation::Reversed {
            dn1.reverse();
        }
        let mut dn2: GpDir = p2.position().direction();
        if f2.0.orientation() == Orientation::Reversed {
            dn2.reverse();
        }
        dn1.angle(&dn2) < tol_ang
    }

    /// `BRepOffset_Tool::CheckBounds(F, Analyse, enlargeU, enlargeVfirst,
    /// enlargeVlast)` (`cxx:3309-3433`): which of the three enlargement
    /// directions are still free given the isolines already on the face.
    ///
    /// Returns `(enlargeU, enlargeVfirst, enlargeVlast)`; all three start at
    /// `true` (`cxx:3316-3318`) and only ever flip to `false`.
    ///
    /// `cxx:3343` reads `L.First().Type()` on a list that OCCT leaves empty
    /// only when the edge is degenerated; `NCollection_List::First` then raises
    /// `Standard_OutOfRange`. The port skips such an edge (documented above
    /// `type_of`: `Analyse::Perform` binds every edge of the shape).
    pub fn check_bounds(f: &Face, analyse: &Analyse) -> (bool, bool, bool) {
        let mut enlarge_u = true;
        let mut enlarge_vfirst = true;
        let mut enlarge_vlast = true;

        let mut u_bound = 0;
        let mut v_bound = 0;
        // `RealLast()` / `RealFirst()` (`cxx:3320-3321`).
        let mut u_first = f64::MAX;
        let mut u_last = -f64::MAX;
        let mut v_first = f64::MAX;
        let mut v_last = -f64::MAX;

        // `cxx:3321-3324`.
        let (uf1, uf2, vf1, vf2) = Self::compact_uv_bounds(f);

        // `cxx:3326-3330`: strip a `Geom_RectangularTrimmedSurface`.
        let mut the_surf = BRepTool::face_surface(f);
        if let Some(s) = the_surf.as_ref() {
            if let Some(basis) = s.rectangular_trimmed_basis() {
                the_surf = Some(basis);
            }
        }
        let Some(the_surf) = the_surf else {
            // Port guard: `theSurf->DynamicType()` at `cxx:3327` faults on a
            // null handle.
            return (enlarge_u, enlarge_vfirst, enlarge_vlast);
        };

        // `cxx:3332-3338`: `Geom_SurfaceOfLinearExtrusion` /
        // `Geom_SurfaceOfRevolution` / `Geom_BezierSurface` /
        // `Geom_BSplineSurface`.
        let is_checked_surface = the_surf.is_surface_of_linear_extrusion()
            || the_surf.is_surface_of_revolution()
            || the_surf.is_bezier_surface()
            || the_surf.is_bspline_surface();

        if is_checked_surface {
            for shape in explorer_shapes(&f.0, ShapeType::Edge) {
                let an_edge = Edge(shape);
                let l = analyse.type_of(&an_edge);
                let degenerated = BRepTool::is_degenerated(&an_edge);
                if l.is_empty() && !degenerated {
                    continue;
                }
                let Some(first) = l.first() else {
                    // Port guard: OCCT raises `Standard_OutOfRange` here.
                    continue;
                };
                let ot = first.kind();
                if ot != TypeOfConcavity::Tangential && !degenerated {
                    continue;
                }

                // `cxx:3347-3352`: the pcurve, unwrapped from
                // `Geom2d_TrimmedCurve`.
                let Some((a_curve, fpar, lpar)) = boptools_2d::curve_on_surface_range(&an_edge, f)
                else {
                    // Port guard: `aCurve->DynamicType()` faults on a null
                    // handle.
                    continue;
                };
                let a_curve: &dyn Curve2d = match a_curve.trimmed_basis() {
                    Some(basis) => basis,
                    None => a_curve.as_ref(),
                };

                // `cxx:3355-3372`: a `Geom2d_Line` pass-through, or the
                // BSpline/Bezier poles converted by
                // `ShapeCustom_Curve2d::ConvertToLine2d`.
                let the_line: Option<Geom2dLine> = if a_curve.is_geom2d_line() {
                    a_curve
                        .gp_lin2d()
                        .map(|lin| Self::line_from_lin2d(&lin))
                } else if a_curve.is_bezier2d() || a_curve.is_bspline2d() {
                    ShapeCustomCurve2d::convert_to_line2d(a_curve, fpar, lpar, CONFUSION)
                        .map(|(line, _, _, _)| line)
                } else {
                    None
                };
                let Some(the_line) = the_line else {
                    continue;
                };

                // `cxx:3374-3413`.
                let the_dir = *the_line.lin().direction();
                let the_loc = the_line.lin().location();
                if the_dir.is_parallel(&GpDir2d::default(), ANGULAR) {
                    // DX2d: an isoline of constant V.
                    v_bound += 1;
                    if degenerated {
                        if (the_loc.y() - vf1).abs() <= CONFUSION {
                            enlarge_vfirst = false;
                        } else {
                            // theLine->Location().Y() is near VF2
                            enlarge_vlast = false;
                        }
                    } else {
                        if the_loc.y() < v_first {
                            v_first = the_loc.y();
                        }
                        if the_loc.y() > v_last {
                            v_last = the_loc.y();
                        }
                    }
                } else if the_dir.is_parallel(&GpDir2d::new(0.0, 1.0).unwrap_or_default(), ANGULAR)
                {
                    // DY2d: an isoline of constant U.
                    u_bound += 1;
                    if the_loc.x() < u_first {
                        u_first = the_loc.x();
                    }
                    if the_loc.x() > u_last {
                        u_last = the_loc.x();
                    }
                }
            }
        }

        // `cxx:3419-3432`.
        if u_bound >= 2 || v_bound >= 2 {
            if u_bound >= 2
                && (uf1 - u_first).abs() <= CONFUSION
                && (uf2 - u_last).abs() <= CONFUSION
            {
                enlarge_u = false;
            }
            if v_bound >= 2
                && (vf1 - v_first).abs() <= CONFUSION
                && (vf2 - v_last).abs() <= CONFUSION
            {
                enlarge_vfirst = false;
                enlarge_vlast = false;
            }
        }

        (enlarge_u, enlarge_vfirst, enlarge_vlast)
    }

    /// `Geom2d_Line(alin)` from the `gp_Lin2d` the pcurve carries.
    fn line_from_lin2d(lin: &GpLin2d) -> Geom2dLine {
        Geom2dLine::from_pnt_dir(lin.location(), *lin.direction())
    }

    /// `CompactUVBounds(F, UMin, UMax, VMin, VMax)` (`cxx:3267-3305`): the
    /// tight UV box of the face's pcurves, sampled at `N = 33` parameters per
    /// edge, with the surface's own bounds as the fall-back.
    fn compact_uv_bounds(f: &Face) -> (f64, f64, f64, f64) {
        let n = 33.0;
        let mut b = BndBox2d::new();
        for shape in explorer_shapes(&f.0, ShapeType::Edge) {
            let e = Edge(shape);
            let c = Curve2dAdaptor::of_edge_face(&e, f);
            // `BRep_Tool::Range(E, U1, U2)`, the edge's 3D parameter range.
            let (u1, u2) = BRepTool::edge_parameters(&e);
            let mut u = u1;
            let du = (u2 - u1) / (n - 1.0);
            for _ in 1..33 {
                if let Some(p) = c.d0(u) {
                    b.add_point(&p);
                }
                u += du;
            }
            if let Some(p) = c.d0(u2) {
                b.add_point(&p);
            }
        }
        match b.get() {
            // `Bnd_Box2d::Get(UMin, VMin, UMax, VMax)`.
            Some((umin, vmin, umax, vmax)) => (umin, umax, vmin, vmax),
            None => match BRepTool::face_surface(f) {
                Some(s) => {
                    let (umin, umax) = s.u_range();
                    let (vmin, vmax) = s.v_range();
                    (umin, umax, vmin, vmax)
                }
                // Port guard: `BRep_Tool::Surface(F)->Bounds` faults on a null
                // handle.
                None => (0.0, 0.0, 0.0, 0.0),
            },
        }
    }
}
