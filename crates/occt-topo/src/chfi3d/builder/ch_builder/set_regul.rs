//! `ChFi3d_ChBuilder::SetRegul` — the regularity encoding of the chamfer
//! result.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_ChBuilder.cxx:2218-2274`
//! (OCCT 8.0.0). This is **not** the fillet body
//! (`ChFi3d_FilBuilder.cxx:2595-2634`, ported in
//! [`crate::chfi3d::builder::fil_builder::set_regul`]): the fillet version
//! calls `ChFi3d_evalconti` on every regular curve, the chamfer version first
//! measures the two surface normals with `BRepAdaptor_Surface` /
//! `BRepAdaptor_Curve2d` and only stores the continuity when the two normals
//! are parallel to within `M_PI / 360` (`cxx:2262-2269`).
//!
//! Translation notes:
//! - `S.Initialize(F, false)` is the unrestricted adaptor
//!   ([`SurfaceAdaptor::of_face`], `BRepAdaptor_Surface.cxx:77-80`),
//! - `PC.Value(t).Coord(u, v)` is [`Curve2dAdaptor::d0`]; a `PC.Value` that
//!   fails is `Err`, OCCT would read an undefined `gp_Pnt2d`,
//! - `n1.Crossed(n2).SquareMagnitude()` is
//!   [`occt_core::gp::GpVec::cross_square_magnitude`],
//! - `reg.IsSurface1() && reg.IsSurface2()` is required by the body
//!   (`cxx:2240`); unlike the fillet version there is no `Merged` fallback, so
//!   a regular curve without two surface indices is skipped entirely,
//! - `myCoup->NewFaces(reg.S1())` takes the `S1()` index of the `ChFiDS_Regul`
//!   exactly as OCCT writes it (`NewFaces(int)` normally receives a *face*
//!   index; the chamfer path passes the surface index),
//! - `BRep_Builder::Continuity(E, F1, F2, C)` is
//!   [`crate::brep_lib_continuity::set_stored_continuity`].

use std::f64::consts::PI;

use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_lib_continuity::set_stored_continuity;
use crate::chfi3d::builder::ChFi3dBuilder;
use crate::chfi3d::builder_helpers::eval_conti;
use crate::shape::{Edge, Face};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

/// `SetRegul()` (`cxx:2218-2274`).
pub(super) fn set_regul(base: &ChFi3dBuilder) -> Result<(), String> {
    // `double Seuil = M_PI / 360.; double Seuil2 = Seuil * Seuil;`
    // (`cxx:2231-2232`).
    const SEUIL: f64 = PI / 360.0;
    let seuil2 = SEUIL * SEUIL;

    for reg in &base.my_regul {
        // `itc.Initialize(myCoup->NewEdges(reg.Curve())); if (itc.More())`
        // (`cxx:2238-2239`).
        let Some(new_edge) = base.my_coup.new_edges(reg.curve()).first().cloned() else {
            continue;
        };
        // `TopoDS_Edge E = TopoDS::Edge(itc.Value());` (`cxx:2240`).
        let e = Edge::wrap(new_edge).ok_or_else(|| {
            "ChSetRegul: the new edge of the regular curve is not an edge".to_string()
        })?;

        // `if (reg.IsSurface1() && reg.IsSurface2())` (`cxx:2241`).
        if !(reg.is_surface1() && reg.is_surface2()) {
            continue;
        }
        // `its1.Initialize(myCoup->NewFaces(reg.S1()));`
        // `its2.Initialize(myCoup->NewFaces(reg.S2())); if (its1.More() &&
        // its2.More())` (`cxx:2243-2246`).
        let (Some(faces1), Some(faces2)) =
            (base.my_coup.new_faces(reg.s1()), base.my_coup.new_faces(reg.s2()))
        else {
            continue;
        };
        let (Some(first1), Some(first2)) = (faces1.first().cloned(), faces2.first().cloned()) else {
            continue;
        };
        let f1 = Face::wrap(first1).ok_or_else(|| {
            "ChSetRegul: the first support of the regular curve is not a face".to_string()
        })?;
        let f2 = Face::wrap(first2).ok_or_else(|| {
            "ChSetRegul: the second support of the regular curve is not a face".to_string()
        })?;

        // `S.Initialize(F1, false); PC.Initialize(E, F1); t = 0.5 *
        // (PC.FirstParameter() + PC.LastParameter()); PC.Value(t).Coord(u, v);
        // S.D1(u, v, p, du, dv); n1 = du.Crossed(dv);` (`cxx:2248-2255`).
        let n1 = normal_at_half(&f1, &e)?;
        // `S.Initialize(F2, false); PC.Initialize(E, F2); ... n2 =
        // du.Crossed(dv);` (`cxx:2257-2259`).
        let n2 = normal_at_half(&f2, &e)?;

        // `if (n1.SquareMagnitude() > 1.e-14 && n2.SquareMagnitude() >
        // 1.e-14) { n1.Normalize(); n2.Normalize(); double sina2 =
        // n1.Crossed(n2).SquareMagnitude(); if (sina2 < Seuil2) { ... } }`
        // (`cxx:2261-2273`).
        if n1.square_magnitude() > 1.0e-14 && n2.square_magnitude() > 1.0e-14 {
            let n1 = n1.normalized();
            let n2 = n2.normalized();
            let sina2 = n1.cross_square_magnitude(&n2);
            if sina2 < seuil2 {
                // `GeomAbs_Shape cont = ChFi3d_evalconti(E, F1, F2);`
                // `B.Continuity(E, F1, F2, cont);` (`cxx:2268-2269`).
                let cont = eval_conti(&f1, &f2);
                set_stored_continuity(&e, &f1, &f2, cont as u8);
            }
        }
    }
    Ok(())
}

/// The unnormalised normal `du ^ dv` of `f` at the mid-parameter of the
/// p-curve of `e` on `f` (`cxx:2248-2259`).
fn normal_at_half(f: &Face, e: &Edge) -> Result<occt_core::gp::GpVec, String> {
    let s = SurfaceAdaptor::of_face(f);
    let pc = Curve2dAdaptor::of_edge_face(e, f);
    let t = 0.5 * (pc.first_parameter() + pc.last_parameter());
    let uv = pc
        .d0(t)
        .ok_or_else(|| "ChSetRegul: BRepAdaptor_Curve2d::Value".to_string())?;
    let (_, du, dv) = s
        .d1(uv.x(), uv.y())
        .ok_or_else(|| "ChSetRegul: BRepAdaptor_Surface::D1".to_string())?;
    Ok(du.crossed(&dv))
}
