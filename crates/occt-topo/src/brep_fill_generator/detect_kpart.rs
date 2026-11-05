//! `DetectKPart` and `CreateKPart` — the analytic particular cases of a ruled
//! patch (`BRepFill_Generator.cxx:62-355` and `:362-569`).
//!
//! `DetectKPart` classifies the pair (section edge, section edge) into the six
//! ruled quadrics `BRepFill_Generator` builds without `GeomFill_Generator`:
//! `1` cylinder, `2` cone, `3` torus, `4` plane, `5` extrusion, `-2` cone with
//! the apex on the first section. `3` and `5` are demoted to `0` at the end
//! (`cxx:349-353`), so only `1`, `2`, `4` and `-2` reach `CreateKPart`.

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx3, GpCirc, GpDir, GpLin, GpPln, GpPnt, GpVec};
use occt_core::precision::{ANGULAR, CONFUSION};
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom::{GeomCone, GeomCylinder, GeomPlane, GeomRectangularTrimmedSurface};

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use crate::shape::Edge;
use crate::topo_tools::edge_vertices;

/// `GeomAdaptor_Curve::Circle()` on a trimmed `Geom_Circle`: the adaptor loads
/// the basis curve (`GeomAdaptor_Curve.cxx:238-313`), so the query has to
/// recurse through `Geom_TrimmedCurve`.
fn circle_of(c: &dyn Curve) -> Option<GpCirc> {
    if let Some(circ) = c.gp_circ() {
        return Some(circ);
    }
    let (basis, _, _) = c.untrimmed_basis()?;
    circle_of(basis.as_ref())
}

/// `GeomAdaptor_Curve::Line()`, mirror of [`circle_of`].
fn line_of(c: &dyn Curve) -> Option<GpLin> {
    if let Some(line) = c.gp_line() {
        return Some(line);
    }
    let (basis, _, _) = c.untrimmed_basis()?;
    line_of(basis.as_ref())
}

/// `gp_Vec(A, B).Dot(gp_Vec(Dir))`.
fn dot_dir(v: &GpVec, d: &GpDir) -> f64 {
    v.dot(&GpVec::new(d.x(), d.y(), d.z()))
}

/// `DetectKPart(Edge1, Edge2)` (`cxx:62-355`). `-1` is the null-3D-curve error
/// the caller turns into `BRepFill_ThruSectionErrorStatus_Null3DCurve`.
pub fn detect_kpart(edge1: &Edge, edge2: &Edge) -> i32 {
    let degen1 = BRepTool::is_degenerated(edge1);
    let mut itype: i32 = 0;

    // `cxx:67-80`.
    let mut first1 = 0.0;
    let mut last1 = 0.0;
    let mut pos1 = GpPnt::zero();
    let mut dist1 = 0.0;
    let mut axe1: Option<GpAx1> = None;
    // The `GeomAdaptor_Curve AdC1` of the first section (`cxx:73`).
    let mut adc1: Option<Box<dyn Curve>> = None;

    if degen1 {
        // `cxx:82-88`.
        itype = -2;
        let (v1, _v2) = edge_vertices(edge1);
        pos1 = BRepTool::vertex_point_world(&v1.expect("Degenerated edge without vertices"));
    } else {
        // `cxx:90-127`.
        let Some(world) = BRepTool::edge_curve_world(edge1) else {
            return -1;
        };
        let (ff, ll) = BRepTool::edge_parameters(edge1);
        let mut owned = world.clone_dyn();
        first1 = ff;
        last1 = ll;
        if edge1.0.orientation() == Orientation::Reversed {
            owned.reverse();
            first1 = owned.reversed_parameter(ll);
            last1 = owned.reversed_parameter(ff);
        }
        match curve_kind(owned.as_ref()) {
            GeomAbsCurveType::Circle => {
                // `cxx:104-111`.
                let circ = circle_of(owned.as_ref()).expect("GeomAdaptor_Curve::Circle");
                itype = 1;
                pos1 = circ.location();
                dist1 = circ.radius();
                axe1 = Some(circ.axis());
            }
            GeomAbsCurveType::Line => {
                // `cxx:112-122`.
                let p_first = owned.d0(first1);
                let p_last = owned.d0(last1);
                itype = 4;
                pos1 = line_of(owned.as_ref())
                    .expect("GeomAdaptor_Curve::Line")
                    .location();
                dist1 = p_first.distance(&p_last);
                let vec = GpVec::from_pnts(&p_first, &p_last);
                let dir = GpDir::from_vec(&vec).expect("gp_Dir(first section chord)");
                axe1 = Some(GpAx1::new(p_first, dir));
            }
            _ => itype = 0,
        }
        adc1 = Some(owned);
    }

    if itype != 0 {
        let degen2 = BRepTool::is_degenerated(edge2);
        if degen2 {
            // `cxx:135-168`.
            let (v1, _v2) = edge_vertices(edge2);
            let pos = BRepTool::vertex_point_world(&v1.expect("Degenerated edge without vertices"));
            if itype == 1 {
                if pos.is_equal_with_tolerance(&pos1, CONFUSION) {
                    // The apex mixes with the circle centre (`cxx:141-145`).
                    itype = 0;
                } else {
                    let dir =
                        GpDir::from_vec(&GpVec::from_pnts(&pos1, &pos)).expect("gp_Dir(apex)");
                    let axe = GpAx1::new(pos1, dir);
                    let axis1 = axe1.expect("cxx:104-111 sets axe1 for IType 1");
                    itype = if axe.direction().is_parallel(axis1.direction()) {
                        2
                    } else {
                        0
                    };
                }
            } else if itype != 4 {
                // Not a plane (`cxx:163-167`).
                itype = 0;
            }
        } else {
            // `cxx:170-345`.
            let Some(world) = BRepTool::edge_curve_world(edge2) else {
                return -1;
            };
            let (f2, l2) = BRepTool::edge_parameters(edge2);
            let mut owned = world.clone_dyn();
            let mut first2 = f2;
            let mut last2 = l2;
            if edge2.0.orientation() == Orientation::Reversed {
                owned.reverse();
                first2 = owned.reversed_parameter(l2);
                last2 = owned.reversed_parameter(f2);
            }

            if itype > 0 && itype < 4 {
                if curve_kind(owned.as_ref()) != GeomAbsCurveType::Circle {
                    // `cxx:189-193`.
                    itype = 0;
                } else {
                    let circ = circle_of(owned.as_ref()).expect("GeomAdaptor_Curve::Circle");
                    let axis1 = axe1.expect("cxx:104-111 sets axe1 for IType 1");
                    if circ.axis().is_coaxial_tol(&axis1, ANGULAR, CONFUSION) {
                        if (circ.radius() - dist1).abs() < CONFUSION {
                            // `cxx:197-221`: cylinder or a piece of one.
                            let h1 = (last1 - first1).abs();
                            let h2 = (last2 - first2).abs();
                            let same_param_length = (h1 - h2).abs() < occt_core::precision::PCONFUSION;
                            let m1 = (first1 + last1) / 2.0;
                            let m2 = (first2 + last2) / 2.0;
                            let (p1, du) = adc1
                                .as_ref()
                                .expect("AdC1 loaded for a circular first section")
                                .d1(m1);
                            let p2 = owned.d0(m2);
                            let same = same_param_length
                                && GpVec::from_pnts(&p1, &p2).is_normal(&du);
                            itype = if same { 1 } else { 0 };
                        } else {
                            // `cxx:222-242`: truncated cone.
                            let h1 = (last1 - first1).abs();
                            let h2 = (last2 - first2).abs();
                            let same_param_length = (h1 - h2).abs() < occt_core::precision::PCONFUSION;
                            let m1 = (first1 + last1) / 2.0;
                            let m2 = (first2 + last2) / 2.0;
                            let (p1, du) = adc1
                                .as_ref()
                                .expect("AdC1 loaded for a circular first section")
                                .d1(m1);
                            let p2 = owned.d0(m2);
                            let same = same_param_length
                                && GpVec::from_pnts(&p1, &p2).is_normal(&du);
                            itype = if same { 2 } else { 0 };
                        }
                        if circ.location().is_equal_with_tolerance(&pos1, CONFUSION) {
                            // `cxx:243-248`: the centres are mixed.
                            itype = 0;
                        }
                    } else if circ.radius() == dist1 {
                        // `cxx:250-258`: torus (demoted to 0 below).
                        itype = 3;
                    } else {
                        itype = 0;
                    }
                }
            } else if itype >= 4 {
                if curve_kind(owned.as_ref()) != GeomAbsCurveType::Line {
                    // `cxx:266-271`.
                    itype = 0;
                } else {
                    let p_first = owned.d0(first2);
                    let p_last = owned.d0(last2);
                    let dist = p_first.distance(&p_last);
                    let vec = GpVec::from_pnts(&p_first, &p_last);
                    let dir_a = GpDir::from_vec(&vec).expect("gp_Dir(second section chord)");
                    let axe = GpAx1::new(p_first, dir_a);
                    let axis1 = axe1.expect("cxx:112-122 sets axe1 for IType 4");
                    if axe.direction().is_parallel(axis1.direction()) {
                        if (dist - dist1).abs() < CONFUSION {
                            let p_start = adc1
                                .as_ref()
                                .expect("AdC1 loaded for a linear first section")
                                .d0(first1);
                            let dir = GpDir::from_vec(&GpVec::from_pnts(&p_start, &p_first))
                                .expect("gp_Dir(first section to second section)");
                            itype = if dir.is_normal(&dir_a) { 4 } else { 5 };
                        } else {
                            itype = 0;
                        }
                    } else {
                        itype = 0;
                    }
                }
            } else if itype == -2 {
                // `cxx:308-344`: apex on the first section.
                if curve_kind(owned.as_ref()) == GeomAbsCurveType::Line {
                    itype = 4;
                } else if curve_kind(owned.as_ref()) == GeomAbsCurveType::Circle {
                    let circ = circle_of(owned.as_ref()).expect("GeomAdaptor_Curve::Circle");
                    let pos = circ.location();
                    let axe = circ.axis();
                    if pos1.is_equal_with_tolerance(&pos, CONFUSION) {
                        itype = 0;
                    } else {
                        let dir =
                            GpDir::from_vec(&GpVec::from_pnts(&pos1, &pos)).expect("gp_Dir(apex)");
                        let apex_axis = GpAx1::new(pos1, dir);
                        itype = if axe.direction().is_parallel(apex_axis.direction()) {
                            -2
                        } else {
                            0
                        };
                    }
                } else {
                    itype = 0;
                }
            }
        }
    }

    // `cxx:349-353`: torus and extrusion are not particular cases.
    if itype == 3 || itype == 5 {
        itype = 0;
    }
    itype
}

/// `CreateKPart(Edge1, Edge2, IType, Surf)` (`cxx:362-569`): the analytic
/// surface for a detected particular case. `false` is the null-3D-curve error.
pub fn create_kpart(
    edge1: &Edge,
    edge2: &Edge,
    itype: i32,
) -> (bool, Option<Arc<dyn Surface>>) {
    // `cxx:366-405`: first section.
    let degen1 = BRepTool::is_degenerated(edge1);
    let (mut v1f, mut v1l) = edge_vertices(edge1);
    let mut c1: Option<Box<dyn Curve>> = None;
    let mut aa = 0.0;
    let mut bb = 0.0;
    if !degen1 {
        let Some(world) = BRepTool::edge_curve_world(edge1) else {
            return (false, None);
        };
        let (a1, b1) = BRepTool::edge_parameters(edge1);
        let mut owned = world.clone_dyn();
        aa = a1;
        bb = b1;
        if edge1.0.orientation() == Orientation::Reversed {
            owned.reverse();
            aa = owned.reversed_parameter(b1);
            bb = owned.reversed_parameter(a1);
            std::mem::swap(&mut v1f, &mut v1l);
        }
        c1 = Some(owned);
    }

    // `cxx:407-441`: second section.
    let degen2 = BRepTool::is_degenerated(edge2);
    let (mut v2f, mut v2l) = edge_vertices(edge2);
    let mut c2: Option<Box<dyn Curve>> = None;
    if !degen2 {
        let Some(world) = BRepTool::edge_curve_world(edge2) else {
            return (false, None);
        };
        let (a2, b2) = BRepTool::edge_parameters(edge2);
        let mut owned = world.clone_dyn();
        if degen1 {
            // `cxx:426-439`: `aa`/`bb` come from the second edge when the first
            // one is degenerated.
            aa = a2;
            bb = b2;
        }
        if edge2.0.orientation() == Orientation::Reversed {
            owned.reverse();
            std::mem::swap(&mut v2f, &mut v2l);
        }
        c2 = Some(owned);
    }

    // `cxx:444-463`: the adaptors only exist over a present curve.
    let adc1 = c1.as_ref().map(|c| c.as_ref());
    let adc2 = c2.as_ref().map(|c| c.as_ref());

    // `cxx:466-564`.
    let surface: Option<Arc<dyn Surface>> = if itype == 1 {
        // `cxx:469-483`: cylinder.
        let circ1 = circle_of(adc1.expect("cylinder needs the first circle"))
            .expect("GeomAdaptor_Curve::Circle");
        let circ2 = circle_of(adc2.expect("cylinder needs the second circle"))
            .expect("GeomAdaptor_Curve::Circle");
        let mut ac1 = GpAx3::from_ax2(&circ1.position());
        let mut v = dot_dir(
            &GpVec::from_pnts(&circ1.location(), &circ2.location()),
            &ac1.direction(),
        );
        if v < 0.0 {
            ac1.z_reverse();
            v = -v;
        }
        let cyl: Arc<dyn Surface> = Arc::new(GeomCylinder::new(
            occt_core::gp::GpCylinder::new(ac1, circ1.radius()).expect("gp_Cylinder"),
        ));
        Some(Arc::new(
            GeomRectangularTrimmedSurface::uv(cyl, aa, bb, v.min(0.0), v.max(0.0))
                .expect("Geom_RectangularTrimmedSurface(cylinder)"),
        ))
    } else if itype == 2 {
        // `cxx:485-508`: cone.
        let k1 = circle_of(adc1.expect("cone needs the first circle"))
            .expect("GeomAdaptor_Curve::Circle");
        let mut ak1 = GpAx3::from_ax2(&k1.position());
        let (mut v, rad) = if degen2 {
            let apex = BRepTool::vertex_point_world(&v2f.expect("degenerated edge vertex"));
            (
                dot_dir(&GpVec::from_pnts(&k1.location(), &apex), &ak1.direction()),
                -k1.radius(),
            )
        } else {
            let k2 = circle_of(adc2.expect("cone needs the second circle"))
                .expect("GeomAdaptor_Curve::Circle");
            (
                dot_dir(&GpVec::from_pnts(&k1.location(), &k2.location()), &ak1.direction()),
                k2.radius() - k1.radius(),
            )
        };
        if v < 0.0 {
            ak1.z_reverse();
            v = -v;
        }
        let ang = (rad / v).atan();
        let cone = occt_core::gp::GpCone::new(ak1, k1.radius(), ang)
            .expect("gp_Cone(truncated cone)");
        let v_trim = v / ang.cos();
        let surf: Arc<dyn Surface> = Arc::new(GeomCone::new(cone));
        Some(Arc::new(
            GeomRectangularTrimmedSurface::uv(surf, aa, bb, v_trim.min(0.0), v_trim.max(0.0))
                .expect("Geom_RectangularTrimmedSurface(cone)"),
        ))
    } else if itype == -2 {
        // `cxx:511-528`: cone with the apex on the first section.
        let k2 = circle_of(adc2.expect("cone needs the second circle"))
            .expect("GeomAdaptor_Curve::Circle");
        let apex = BRepTool::vertex_point_world(&v1f.expect("degenerated edge vertex"));
        let mut ak2 = GpAx3::from_ax2(&k2.position());
        ak2.set_location(apex);
        let mut v = dot_dir(&GpVec::from_pnts(&apex, &k2.location()), &ak2.direction());
        let rad = k2.radius();
        if v < 0.0 {
            ak2.z_reverse();
            v = -v;
        }
        let ang = (rad / v).atan();
        let cone = occt_core::gp::GpCone::new(ak2, 0.0, ang).expect("gp_Cone(apex cone)");
        let v_trim = v / ang.cos();
        let surf: Arc<dyn Surface> = Arc::new(GeomCone::new(cone));
        Some(Arc::new(
            GeomRectangularTrimmedSurface::uv(surf, aa, bb, v_trim.min(0.0), v_trim.max(0.0))
                .expect("Geom_RectangularTrimmedSurface(apex cone)"),
        ))
    } else if itype == 4 {
        // `cxx:532-558`: plane.
        let l1 = (!degen1)
            .then(|| line_of(adc1.expect("plane needs the first line")))
            .flatten();
        let l2 = (!degen2)
            .then(|| line_of(adc2.expect("plane needs the second line")))
            .flatten();
        let a_line = l2.or(l1).expect("IType 4 has at least one straight section");
        let p1 = if degen1 {
            BRepTool::vertex_point_world(&v1f.expect("degenerated edge vertex"))
        } else {
            l1.expect("plane needs the first line").location()
        };
        let p2 = if degen2 {
            BRepTool::vertex_point_world(&v2f.expect("degenerated edge vertex"))
        } else {
            l2.expect("plane needs the second line").location()
        };
        let p1p2 = GpVec::from_pnts(&p1, &p2);
        let d1 = a_line.direction();
        let normal = GpDir::from_vec(
            &GpVec::new(d1.x(), d1.y(), d1.z()).crossed(&p1p2),
        )
        .expect("gp_Dir(normal of the plane)");
        let ax = GpAx3::new(a_line.location(), normal, &d1).expect("gp_Ax3(plane)");
        Some(Arc::new(GeomPlane::new(GpPln::new(ax))))
    } else {
        None
    };
    (true, surface)
}
