//! `LineConstructor` (`Contap_Contour.cxx:283-511`): cuts a contour line into
//! the pieces that lie inside the surface domain, one piece per pair of
//! consecutive vertices.

use occt_core::elib::clib;
use occt_core::elib::slib;
use occt_core::gp::GpPnt2d;
use occt_core::precision::{PCONFUSION, RESOLUTION};

use occt_topo::brep_top_adaptor_topol_tool::TopolTool;
use occt_topo::fclass2d::FaceState;
use occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use occt_topo::int_tools_wline::WLine;
use occt_topo::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::super::line::ContapLine;
use super::super::types::ContapIType;
use super::recadre::recadre;

/// `LineConstructor(slin, Domain, L, Surf)` (`cxx:283-511`).
///
/// The three type arms append to `slin` and never modify `l`: OCCT takes
/// `Contap_Line&` only because `Vertex(i)` is not `const`-qualified.
pub(crate) fn line_constructor(
    slin: &mut Vec<ContapLine>,
    domain: &mut TopolTool,
    l: &ContapLine,
    surf: &SurfaceAdaptor,
) {
    // `cxx:291`: `constexpr double Tol = Precision::PConfusion()`.
    let tol = PCONFUSION;
    let typ_l = l.type_contour();

    if typ_l == ContapIType::Walking {
        let nbvtx = l.nb_vertex();
        // `for (i = 1; i < nbvtx; i++)` (`cxx:299`).
        for i in 1..nbvtx {
            let firstp = l.vertex(i).map(|v| v.parameter_on_line()).unwrap_or(0.0) as i32;
            let lastp = l
                .vertex(i + 1)
                .map(|v| v.parameter_on_line())
                .unwrap_or(0.0) as i32;
            if firstp == lastp {
                continue;
            }
            // `cxx:305-307`: the integer midpoint parameter and its point.
            let pmid = (firstp + lastp) / 2;
            let Ok(pmid_pt) = l.point(pmid as usize) else {
                continue;
            };
            let (_u1, _v1, mut u2, mut v2) = pmid_pt.parameters();
            recadre(surf, &mut u2, &mut v2);
            let in2 = domain.classify(GpPnt2d::new(u2, v2), tol, true);
            if in2 == FaceState::Out {
                continue;
            }

            // `cxx:317-327`.
            let mut line_on_2s = WLine::new();
            if let Ok(p) = l.point(firstp as usize) {
                line_on_2s.add(*p);
            }
            for j in (firstp + 1)..=lastp {
                let Ok(pj) = l.point(j as usize) else { break };
                let Ok(pjm1) = l.point((j - 1) as usize) else {
                    break;
                };
                let a_sq_dist = pj.value().square_distance(&pjm1.value());
                if a_sq_dist > RESOLUTION {
                    line_on_2s.add(*pj);
                }
            }
            // `cxx:328-331`.
            if line_on_2s.nb_pnts() < 2 {
                continue;
            }
            let mut line = ContapLine::new();
            line.set_line_on_2s(line_on_2s);

            // `cxx:333-339`.
            if let Some(v) = l.vertex(i) {
                let mut pvtx = v.clone();
                pvtx.set_parameter(1.0);
                line.add_vertex(pvtx);
            }
            if let Some(v) = l.vertex(i + 1) {
                let mut pvtx = v.clone();
                let nbp = line.line_on_2s().map(|w| w.nb_pnts()).unwrap_or(0);
                pvtx.set_parameter(f64::from(nbp));
                line.add_vertex(pvtx);
            }
            line.set_transition_on_s(l.transition_on_s());
            slin.push(line);
        }
    } else if typ_l == ContapIType::Lin {
        let nbvtx = l.nb_vertex();
        let Ok(the_line) = l.line() else {
            return;
        };
        for i in 1..nbvtx {
            let firstp = l.vertex(i).map(|v| v.parameter_on_line()).unwrap_or(0.0);
            let lastp = l
                .vertex(i + 1)
                .map(|v| v.parameter_on_line())
                .unwrap_or(0.0);
            if firstp == lastp {
                continue;
            }
            let pmid = (firstp + lastp) * 0.5;
            let p_mid = clib::line_value(&the_line, pmid);
            // `cxx:359-370`: `ElSLib::Parameters` for a cylinder or a cone;
            // any other surface leaves `u2`/`v2` untouched (the source only
            // prints a diagnostic).
            let (mut u2, mut v2) = (0.0, 0.0);
            match surf.get_type() {
                GeomAbsSurfaceType::Cylinder => {
                    if let Some(cy) = surf.cylinder() {
                        let (u, v) = slib::cylinder_parameters(&cy.position(), &p_mid);
                        u2 = u;
                        v2 = v;
                    }
                }
                GeomAbsSurfaceType::Cone => {
                    if let Some(co) = surf.cone() {
                        let (u, v) =
                            slib::cone_parameters(&co.position(), co.radius(), co.semi_angle(), &p_mid);
                        u2 = u;
                        v2 = v;
                    }
                }
                _ => {}
            }
            recadre(surf, &mut u2, &mut v2);
            let in2 = domain.classify(GpPnt2d::new(u2, v2), tol, true);
            if in2 == FaceState::Out {
                continue;
            }

            // `cxx:381-389`.
            let mut line = ContapLine::new();
            line.set_value_line(&the_line);
            if let Some(v) = l.vertex(i) {
                line.add_vertex(v.clone());
            }
            if let Some(v) = l.vertex(i + 1) {
                line.add_vertex(v.clone());
            }
            line.set_transition_on_s(l.transition_on_s());
            slin.push(line);
        }
    } else if typ_l == ContapIType::Circle {
        let nbvtx = l.nb_vertex();
        let Ok(the_circle) = l.circle() else {
            return;
        };
        let mut novtx = true;
        if nbvtx != 0 {
            novtx = false;
        }
        // `for (i = 1; i < nbvtx || novtx; i++)` (`cxx:405`).
        let mut i = 1usize;
        while i < nbvtx || novtx {
            let mut firstp = 0.0f64;
            let mut lastp = std::f64::consts::TAU;
            if !novtx {
                firstp = l.vertex(i).map(|v| v.parameter_on_line()).unwrap_or(0.0);
                lastp = l
                    .vertex(i + 1)
                    .map(|v| v.parameter_on_line())
                    .unwrap_or(0.0);
            }
            if (firstp - lastp).abs() > 0.000000001 {
                let pmid = (firstp + lastp) * 0.5;
                let p_mid = clib::circle_value(&the_circle, pmid);
                let (mut u2, mut v2) = (0.0, 0.0);
                match surf.get_type() {
                    GeomAbsSurfaceType::Cylinder => {
                        if let Some(cy) = surf.cylinder() {
                            let (u, v) = slib::cylinder_parameters(&cy.position(), &p_mid);
                            u2 = u;
                            v2 = v;
                        }
                    }
                    GeomAbsSurfaceType::Cone => {
                        if let Some(co) = surf.cone() {
                            let (u, v) = slib::cone_parameters(
                                &co.position(),
                                co.radius(),
                                co.semi_angle(),
                                &p_mid,
                            );
                            u2 = u;
                            v2 = v;
                        }
                    }
                    GeomAbsSurfaceType::Sphere => {
                        if let Some(sp) = surf.sphere() {
                            let (u, v) = slib::sphere_parameters(sp.position(), &p_mid);
                            u2 = u;
                            v2 = v;
                        }
                    }
                    _ => {}
                }
                recadre(surf, &mut u2, &mut v2);
                let in2 = domain.classify(GpPnt2d::new(u2, v2), tol, true);
                if in2 != FaceState::Out {
                    // `cxx:441-452`.
                    let mut line = ContapLine::new();
                    line.set_value_circle(&the_circle);
                    if !novtx {
                        if let Some(v) = l.vertex(i) {
                            line.add_vertex(v.clone());
                        }
                        if let Some(v) = l.vertex(i + 1) {
                            line.add_vertex(v.clone());
                        }
                    }
                    line.set_transition_on_s(l.transition_on_s());
                    slin.push(line);
                }
            }
            novtx = false;
            i += 1;
        }

        // `cxx:457-503`: the closing piece from the last vertex to the first.
        if nbvtx != 0 {
            let firstp = l
                .vertex(nbvtx)
                .map(|v| v.parameter_on_line())
                .unwrap_or(0.0);
            let lastp = l.vertex(1).map(|v| v.parameter_on_line()).unwrap_or(0.0)
                + std::f64::consts::TAU;
            if (firstp - lastp).abs() > 0.0000000001 {
                let pmid = (firstp + lastp) * 0.5;
                let p_mid = clib::circle_value(&the_circle, pmid);
                let (mut u2, mut v2) = (0.0, 0.0);
                match surf.get_type() {
                    GeomAbsSurfaceType::Cylinder => {
                        if let Some(cy) = surf.cylinder() {
                            let (u, v) = slib::cylinder_parameters(&cy.position(), &p_mid);
                            u2 = u;
                            v2 = v;
                        }
                    }
                    GeomAbsSurfaceType::Cone => {
                        if let Some(co) = surf.cone() {
                            let (u, v) = slib::cone_parameters(
                                &co.position(),
                                co.radius(),
                                co.semi_angle(),
                                &p_mid,
                            );
                            u2 = u;
                            v2 = v;
                        }
                    }
                    GeomAbsSurfaceType::Sphere => {
                        if let Some(sp) = surf.sphere() {
                            let (u, v) = slib::sphere_parameters(sp.position(), &p_mid);
                            u2 = u;
                            v2 = v;
                        }
                    }
                    _ => {}
                }
                recadre(surf, &mut u2, &mut v2);
                let in2 = domain.classify(GpPnt2d::new(u2, v2), tol, true);
                if in2 != FaceState::Out {
                    // `cxx:492-500`.
                    let mut line = ContapLine::new();
                    line.set_value_circle(&the_circle);
                    if let Some(v) = l.vertex(nbvtx) {
                        line.add_vertex(v.clone());
                    }
                    if let Some(v) = l.vertex(1) {
                        let mut pvtx = v.clone();
                        pvtx.set_parameter(pvtx.parameter_on_line() + std::f64::consts::TAU);
                        line.add_vertex(pvtx);
                    }
                    line.set_transition_on_s(l.transition_on_s());
                    slin.push(line);
                }
            }
        }
    } else {
        // `cxx:504-509`: neither walked line, line nor circle.
        slin.push(l.clone());
    }
}
