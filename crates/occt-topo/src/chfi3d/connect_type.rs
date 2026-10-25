//! `ChFi3d::DefineConnectType` and its file-local `Correct2dPoint`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d.cxx:39-156` and
//! `:650-691`, declared in `ChFi3d.hxx:38-41`.

use occt_core::gp::GpPnt2d;
use occt_core::precision::{Precision, PCONFUSION, REAL_SMALL};

use crate::abs::Orientation;
use crate::boptools_2d::{curve_on_surface_range, intermediate_point};
use crate::brep_tool::BRepTool;
use crate::brep_tools::ori_edge_in_face;
use crate::chfi_ds::TypeOfConcavity;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::shape::{Edge, Face};
use crate::topopebrep::tool::adaptor::{CurveAdaptor, SurfaceAdaptor};

/// `ChFi3d::DefineConnectType(E, F1, F2, SinTol, CorrectPoint)`
/// (`cxx:39-156`): `Convex`, `Concave` or `Tangential` from the edge tangent
/// and the two surface normals at the middle of the edge.
///
/// `Err` carries the `Standard_ConstructionError` of
/// `BRepTools::OriEdgeInFace` (`cxx:79`) or the null handle that OCCT
/// dereferences (`cxx:96-97`, `:108-109`).
pub fn define_connect_type(
    e: &Edge,
    f1: &Face,
    f2: &Face,
    sin_tol: f64,
    correct_point: bool,
) -> Result<TypeOfConcavity, String> {
    let s1 = BRepTool::face_surface_world(f1)
        .ok_or_else(|| "ChFi3d::DefineConnectType: null surface F1".to_string())?;
    let s2 = BRepTool::face_surface_world(f2)
        .ok_or_else(|| "ChFi3d::DefineConnectType: null surface F2".to_string())?;
    // `BRep_Tool::CurveOnSurface(E, F1, f, l)`, then the seam arm
    // (`cxx:51-58`) re-reads the same edge reversed on F2 when F1 is F2.
    let c1 = curve_on_surface_range(e, f1).map(|(c, _, _)| c);
    let mut ee = e.clone();
    if crate::topo_tools::is_same(&f1.0, &f2.0) {
        ee.0.reverse();
    }
    let c2 = curve_on_surface_range(&ee, f2).map(|(c, _, _)| c);
    let (c1, c2) = match (c1, c2) {
        (Some(c1), Some(c2)) => (c1, c2),
        // `if (C1.IsNull() || C2.IsNull()) return ChFiDS_Other;` (`cxx:58-61`).
        _ => return Ok(TypeOfConcavity::Other),
    };

    let curve = CurveAdaptor::of_edge(e);
    let f = curve.first_parameter();
    let l = curve.last_parameter();
    let mut par_on_c = 0.5 * (f + l);
    let d1 = curve
        .d1(par_on_c)
        .ok_or_else(|| "ChFi3d::DefineConnectType: no edge curve".to_string())?;
    let mut t1 = d1.1;
    if t1.square_magnitude() <= REAL_SMALL {
        par_on_c = intermediate_point(f, l);
        let d1b = curve
            .d1(par_on_c)
            .ok_or_else(|| "ChFi3d::DefineConnectType: no edge curve".to_string())?;
        t1 = d1b.1;
    }
    if t1.square_magnitude() > REAL_SMALL {
        t1.normalize();
    }

    // `BRepTools::OriEdgeInFace(E, F1) == TopAbs_REVERSED` (`cxx:79-82`).
    if ori_edge_in_face(e, f1)? == Orientation::Reversed {
        t1.reverse();
    }
    if f1.0.orientation() == Orientation::Reversed {
        t1.reverse();
    }

    let mut p = c1.value(par_on_c);
    if correct_point {
        correct2d_point(f1, &mut p);
    }
    let (_, d1u, d1v) = s1.d1(p.x(), p.y());
    let mut dn1 = d1u.crossed(&d1v);
    if f1.0.orientation() == Orientation::Reversed {
        dn1.reverse();
    }

    let mut p2 = c2.value(par_on_c);
    if correct_point {
        correct2d_point(f2, &mut p2);
    }
    let (_, d1u, d1v) = s2.d1(p2.x(), p2.y());
    let mut dn2 = d1u.crossed(&d1v);
    if f2.0.orientation() == Orientation::Reversed {
        dn2.reverse();
    }

    dn1.normalize();
    dn2.normalize();

    let pro_vec = dn1.crossed(&dn2);
    let norm_pro_vec = pro_vec.magnitude();
    if norm_pro_vec < sin_tol {
        // plane (`cxx:123-138`).
        if dn1.dot(&dn2) > 0.0 {
            Ok(TypeOfConcavity::Tangential)
        } else {
            // "Mixed not finished!" (`cxx:132-136`): OCCT still answers Convex.
            Ok(TypeOfConcavity::Convex)
        }
    } else {
        let pro = if norm_pro_vec > REAL_SMALL {
            pro_vec.divided(norm_pro_vec)
        } else {
            pro_vec
        };
        if t1.dot(&pro) > 0.0 {
            Ok(TypeOfConcavity::Convex)
        } else {
            Ok(TypeOfConcavity::Concave)
        }
    }
}

/// `Correct2dPoint(theF, theP2d)` (`ChFi3d.cxx:650-691`): nudge a `(u, v)`
/// that sits within `1%` of a Bezier/B-spline parameter bound off that bound.
pub fn correct2d_point(face: &Face, p2d: &mut GpPnt2d) {
    let a_bas = SurfaceAdaptor::of_face(face);
    // `if (aBAS.GetType() < GeomAbs_BezierSurface) return;` (`cxx:653-656`).
    // The enum order of `GeomAbs_SurfaceType` is the one `GeomAbsSurfaceType`
    // keeps (`geom_bnd_lib_surface3d.rs:43-55`).
    if (a_bas.get_type() as u8) < (GeomAbsSurfaceType::BezierSurface as u8) {
        return;
    }
    const COEFF: f64 = 0.01;
    // `aBAS.Initialize(theF, true)` (`cxx:663`).
    let bas = SurfaceAdaptor::of_face_restricted(face);
    let u1 = bas.first_u_parameter();
    let u2 = bas.last_u_parameter();
    let v1 = bas.first_v_parameter();
    let v2 = bas.last_v_parameter();
    if !(Precision::is_infinite(u1) || Precision::is_infinite(u2)) {
        let eps = (COEFF * (u2 - u1)).max(PCONFUSION);
        if (p2d.x() - u1).abs() < eps {
            p2d.set_x(u1 + eps);
        }
        if (p2d.x() - u2).abs() < eps {
            p2d.set_x(u2 - eps);
        }
    }
    if !(Precision::is_infinite(v1) || Precision::is_infinite(v2)) {
        let eps = (COEFF * (v2 - v1)).max(PCONFUSION);
        if (p2d.y() - v1).abs() < eps {
            p2d.set_y(v1 + eps);
        }
        if (p2d.y() - v2).abs() < eps {
            p2d.set_y(v2 - eps);
        }
    }
}
