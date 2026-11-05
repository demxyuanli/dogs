//! The file statics of `ChFi3d_Builder_CnCrn.cxx` that the `GeomPlate` section
//! of `PerformMoreThreeCorner` calls.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx`
//! (OCCT 8.0.0):
//!
//! - `Calcul_Orientation` (`cxx:205-218`),
//! - `CalculDroite` (`cxx:534-543`),
//! - `CalculBatten` (`cxx:548-677`),
//! - `OrientationIcNonVive` (`cxx:683-715`),
//! - `OrientationIcplusNonVive` (`cxx:721-754`),
//! - `OrientationAreteViveConsecutive` (`cxx:762-806`),
//! - `ParametrePlate` (`cxx:1084-1108`),
//! - `PlateOrientation` (`cxx:1180-1228`).
//!
//! Translation notes:
//! - the four `TopAbs_Orientation&` out-parameters become the return value; the
//!   `void` statics keep `Result` only for the reads that can raise in OCCT
//!   (`GeomLProp_CLProps2d::Tangent` on a non-C1 point, a null p-curve),
//! - `GeomAdaptor_Surface::GetType() == GeomAbs_Plane` is
//!   [`surface_kind`], the same reading `GeomAdaptor_Surface::Load` performs for
//!   a trimmed plane (`GeomAdaptor_Surface.cxx:...` / `surface_kind`'s note),
//! - `BndLib_Add2dCurve::Add(Geom2dAdaptor_Curve(pcurve), 0, bc)` is
//!   [`add_adaptor`], the Geom2d-handle overload,
//! - a null handle in `thePCArr` raises in OCCT the moment `Value(i)` is
//!   dereferenced; here that is the `Err` of [`plate_orientation`].

use std::sync::Arc;

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpDir2d, GpPnt, GpPnt2d, GpVec, GpVec2d};
use occt_geom::geom_plate::BuildPlateSurface;
use occt_geom::Surface;
use occt_geom2d::cl_props2d::GeomLPropCLProps2d;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::fair_curve::Batten;
use occt_geom2d::line::Geom2dLine;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::Orientation;
use crate::bnd_lib_add2d::add_adaptor;
use crate::brep_tool::BRepTool;
use crate::chfi_ds::Stripe;
use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};
use crate::iterator::cumulated_children;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{edge_vertices, is_same};

/// `Calcul_Orientation(stripe, jfposit, indice, orient)` (`cxx:205-218`): the
/// transition of the `S1` interference when `jfposit == 1`, of the `S2` one
/// otherwise. `indice` is the 1-based surf data index.
pub fn calcul_orientation(
    stripe: &Stripe,
    jfposit: i32,
    indice: i32,
) -> Result<Orientation, String> {
    if indice < 1 {
        return Err(format!("Calcul_Orientation: surf data index {indice} is out of range"));
    }
    let data = stripe
        .set_of_surf_data()
        .get(indice as usize - 1)
        .ok_or_else(|| format!("Calcul_Orientation: surf data index {indice} is out of range"))?;
    let fi = if jfposit == 1 {
        data.interference_on_s1()
    } else {
        data.interference_on_s2()
    };
    Ok(fi.transition())
}

/// `CalculDroite(p2d1, xdir, ydir, pcurve)` (`cxx:534-543`): the 2D line from
/// `p2d1` along `(xdir, ydir)` over `[0, |(xdir, ydir)|]`.
pub fn calcul_droite(
    p2d1: &GpPnt2d,
    xdir: f64,
    ydir: f64,
) -> Result<Arc<dyn Curve2d>, String> {
    let dir1 = GpDir2d::new(xdir, ydir)
        .map_err(|_| "CalculDroite: gp_Dir2d on a null vector".to_string())?;
    let line = Geom2dLine::from_pnt_dir(*p2d1, dir1);
    let l0 = (xdir * xdir + ydir * ydir).sqrt();
    Ok(Arc::new(Geom2dTrimmedCurve::new(Arc::new(line), 0.0, l0)))
}

/// `CalculBatten(ASurf, Face, xdir, ydir, p2d1, p2d2, contraint1, contraint2,
/// curv2d1, curv2d2, picicplus, picplusic, inverseic, inverseicplus, pcurve)`
/// (`cxx:548-677`).
///
/// The straight line is the fallback on every arm OCCT abandons the batten:
/// a sharp angle against a non-planar surface, a `Compute` failure, and a
/// resulting curve that leaves the face's UV bounds by more than `1.e-7`.
#[allow(clippy::too_many_arguments)]
pub fn calcul_batten(
    asurf: &Arc<dyn Surface>,
    face: &Face,
    xdir: f64,
    ydir: f64,
    p2d1: &GpPnt2d,
    p2d2: &GpPnt2d,
    contraint1: bool,
    contraint2: bool,
    curv2d1: &Arc<dyn Curve2d>,
    curv2d2: &Arc<dyn Curve2d>,
    picicplus: f64,
    picplusic: f64,
    inverseic: bool,
    inverseicplus: bool,
) -> Result<Arc<dyn Curve2d>, String> {
    // `bool anglebig = false; isplane = ASurf->GetType() == GeomAbs_Plane;`
    // (`cxx:565-568`).
    let mut isplane = surface_kind(asurf.as_ref()) == GeomAbsSurfaceType::Plane;
    let dir1 = GpDir2d::new(xdir, ydir)
        .map_err(|_| "CalculBatten: gp_Dir2d on a null vector".to_string())?;

    // `GeomLProp_CLProps2d CL1(curv2d1, picicplus, 1, 1.e-4); CL1.Tangent(dir3);`
    // (`cxx:570-574`).
    let mut cl1 = GeomLPropCLProps2d::new(curv2d1.as_ref(), picicplus, 1, 1.0e-4);
    if !cl1.is_tangent_defined() {
        return Err("CalculBatten: the 2D tangent is not defined on curv2d1".to_string());
    }
    let mut dir3 = cl1.tangent();
    // `GeomLProp_CLProps2d CL2(curv2d2, picplusic, 1, 1.e-4); CL2.Tangent(dir4);`
    // (`cxx:571-575`).
    let mut cl2 = GeomLPropCLProps2d::new(curv2d2.as_ref(), picplusic, 1, 1.0e-4);
    if !cl2.is_tangent_defined() {
        return Err("CalculBatten: the 2D tangent is not defined on curv2d2".to_string());
    }
    let mut dir4 = cl2.tangent();

    // `cxx:576-583`.
    if inverseic {
        dir3.reverse();
    }
    if inverseicplus {
        dir4.reverse();
    }

    // `cxx:584-596`: the two angles against `dir1`. `ang2` is the reflection of
    // `dir1.Angle(dir4)` into `]-pi, pi]` around `+/-pi`.
    let ang1 = dir1.angle(&dir3);
    let mut anglebig = false;
    let a4 = dir1.angle(&dir4);
    let ang2 = if a4 > 0.0 {
        std::f64::consts::PI - a4
    } else {
        -std::f64::consts::PI - a4
    };

    // `cxx:597-608`.
    if contraint1 && contraint2 {
        anglebig = ang1.abs() > 1.2 || ang2.abs() > 1.2;
    } else if contraint1 {
        anglebig = ang1.abs() > 1.2;
    } else if contraint2 {
        anglebig = ang2.abs() > 1.2;
    }

    // `cxx:609-612`.
    if isplane
        && (ang1.abs() > std::f64::consts::FRAC_PI_2
            || ang2.abs() > std::f64::consts::FRAC_PI_2)
    {
        isplane = false;
    }

    // `cxx:613-616`.
    if anglebig && !isplane {
        return calcul_droite(p2d1, xdir, ydir);
    }

    // `cxx:617-630`: the batten with the requested end constraints.
    let h = p2d2.distance(p2d1) / 20.0;
    let mut bat = Batten::new(*p2d1, *p2d2, h, 0.0)
        .map_err(|e| format!("CalculBatten: FairCurve_Batten: {e}"))?;
    // `Bat.SetFreeSliding(true)` (`cxx:586`).
    bat.set_free_sliding(true);
    if contraint1 {
        bat.set_angle1(ang1);
    } else {
        bat.set_constraint_order1(0);
    }
    if contraint2 {
        bat.set_angle2(ang2);
    } else {
        bat.set_constraint_order2(0);
    }

    // `Ok = Bat.Compute(Iana, 25, 1.e-2);` (`cxx:634`).
    let (mut ok, _iana) = bat.compute(25, 1.0e-2);

    // `cxx:643-671`: the result must stay inside the face's UV box.
    let mut pcurve: Option<Arc<dyn Curve2d>> = None;
    if ok {
        if let Some(c) = bat.curve() {
            let curve: Arc<dyn Curve2d> = Arc::new(c);
            let (umin, umax, vmin, vmax) = BRepTool::uv_bounds(face);
            let mut bc = BndBox2d::new();
            add_adaptor(curve.as_ref(), 0.0, &mut bc);
            let mut bf = BndBox2d::new();
            bf.update(umin, vmin, umax, vmax);
            match bc.get() {
                Some((uminc, vminc, umaxc, vmaxc)) => {
                    if uminc < umin - 1.0e-7 {
                        ok = false;
                    }
                    if umaxc > umax + 1.0e-7 {
                        ok = false;
                    }
                    if vminc < vmin - 1.0e-7 {
                        ok = false;
                    }
                    if vmaxc > vmax + 1.0e-7 {
                        ok = false;
                    }
                }
                // `bc` holds no point: `Bnd_Box2d::Get` raises on a void box in
                // OCCT (`Bnd_Box2d.lxx:...`); the port keeps the curve.
                None => {}
            }
            if ok {
                pcurve = Some(curve);
            }
        } else {
            ok = false;
        }
    }

    // `if (!Ok) CalculDroite(p2d1, xdir, ydir, pcurve);` (`cxx:671-674`).
    match (ok, pcurve) {
        (true, Some(c)) => Ok(c),
        _ => calcul_droite(p2d1, xdir, ydir),
    }
}

/// `OrientationIcNonVive(CDic, jfic, icicplus, sensic, orien)` (`cxx:683-715`):
/// the curve orientation between `ic` and `icplus` when `ic` is not a living
/// edge. `sensic == 1` reverses the interference transition.
pub fn orientation_ic_non_vive(
    stripe: &Stripe,
    jfic: i32,
    icicplus: i32,
    sensic: i32,
) -> Result<Orientation, String> {
    let orinterf = calcul_orientation(stripe, jfic, icicplus)?;
    Ok(if (sensic != 1) == orinterf.is_forward() {
        Orientation::Forward
    } else {
        Orientation::Reversed
    })
}

/// `OrientationIcplusNonVive(CDicplus, jficplus, icplusic, sensicplus, orien)`
/// (`cxx:721-754`): the `icplus` mirror, with the `jfp = 3 - jficplus` flip
/// and the opposing sense test.
pub fn orientation_icplus_non_vive(
    stripe: &Stripe,
    jficplus: i32,
    icplusic: i32,
    sensicplus: i32,
) -> Result<Orientation, String> {
    let jfp = 3 - jficplus;
    let orinterf = calcul_orientation(stripe, jfp, icplusic)?;
    Ok(if (sensicplus == 1) == orinterf.is_forward() {
        Orientation::Forward
    } else {
        Orientation::Reversed
    })
}

/// `OrientationAreteViveConsecutive(Fviveicicplus, Eviveic, V1, orien)`
/// (`cxx:762-806`): the orientation of the curve between two consecutive living
/// edges, read from the orientation of `Eviveic` in the forward `Fvive` face.
pub fn orientation_arete_vive_consecutive(
    fviveicicplus: &TopoShape,
    eviveic: &Edge,
    v1: &TopoShape,
) -> Result<Orientation, String> {
    // `TopoDS_Face F = TopoDS::Face(Fviveicicplus);` (`cxx:769`).
    let face = Face::wrap(fviveicicplus.clone())
        .ok_or_else(|| "OrientationAreteViveConsecutive: Fvive is not a face".to_string())?;
    // `for (ex.Init(F.Oriented(TopAbs_FORWARD), TopAbs_EDGE); ...)` (`cxx:772-779`).
    let mut orinterf = Orientation::Forward;
    for child in cumulated_children(&face.0.oriented(Orientation::Forward)) {
        if child.shape_type() == crate::abs::ShapeType::Edge && is_same(&child, &eviveic.0) {
            orinterf = child.orientation();
            break;
        }
    }
    // `vl = TopExp::LastVertex(E);` (`cxx:782-783`).
    let (_first, last) = edge_vertices(eviveic);
    let vl_is_v1 = last.as_ref().is_some_and(|v| is_same(&v.0, v1));
    Ok(if vl_is_v1 == orinterf.is_forward() {
        Orientation::Forward
    } else {
        Orientation::Reversed
    })
}

/// `ParametrePlate(n3d, PSurf, Surf, point, apperror, uv)` (`cxx:1084-1108`):
/// searches the ends of the `n3d` constraint p-curves for the UV point whose
/// 3D image is `point` within `apperror`.
///
/// OCCT leaves `uv` at the last value it computed when no end matches; the port
/// returns that same value (`Gp_Pnt2d()` — the origin — when `n3d == 0`).
pub fn parametre_plate(
    n3d: usize,
    psurf: &BuildPlateSurface,
    surf: &Arc<dyn Surface>,
    point: &GpPnt,
    apperror: f64,
) -> Result<GpPnt2d, String> {
    let curves = psurf.curves2d()?;
    let mut uv = GpPnt2d::new(0.0, 0.0);
    let mut trouve = false;
    for ip in 1..=n3d {
        if trouve {
            break;
        }
        let pc = curves
            .get(ip - 1)
            .and_then(|c| c.as_ref())
            .ok_or_else(|| format!("ParametrePlate: Curves2d({ip}) is null"))?;
        // `par = FirstParameter(); D0(par, uv); Surf->D0(uv, P1);
        // trouve = P1.IsEqual(point, apperror);` (`cxx:1096-1101`).
        uv = pc.d0(pc.first_parameter());
        let p1 = surf.d0(uv.x(), uv.y());
        trouve = p1.is_equal_with_tolerance(point, apperror);
        if !trouve {
            // `par = LastParameter(); D0(par, uv); ...` (`cxx:1102-1107`).
            uv = pc.d0(pc.last_parameter());
            let p1 = surf.d0(uv.x(), uv.y());
            trouve = p1.is_equal_with_tolerance(point, apperror);
        }
    }
    Ok(uv)
}

/// `PlateOrientation(thePlateSurf, thePCArr, theRefDir)` (`cxx:1180-1228`):
/// answers the plate orientation by comparing the plate's normal against the
/// sum of the face normals at the corner, sampled along every constraint curve.
///
/// OCCT computes `pp2` at `cxx:1198` and overwrites it with the `D1` point of
/// the first loop pass at `cxx:1205` without reading it; the port keeps the
/// same statement order, so the first assignment is dead by construction.
#[allow(unused_assignments)]
pub fn plate_orientation(
    plate_surf: &Arc<dyn Surface>,
    pc_arr: &[Option<Arc<dyn Curve2d>>],
    ref_dir: &GpVec,
) -> Result<Orientation, String> {
    // `int nb = thePCArr->Upper(); aPC = thePCArr->Value(nb);` (`cxx:1192-1193`).
    let nb = pc_arr.len();
    if nb == 0 {
        return Err("PlateOrientation: the p-curve array is empty".to_string());
    }
    // `fpar = aPC->FirstParameter(); lpar = aPC->LastParameter();
    //  aPC->D0((fpar + lpar) / 2., uv); thePlateSurf->D0(uv, pp1);
    //  aPC->D0(lpar, uv); thePlateSurf->D0(uv, pp2);` (`cxx:1194-1198`).
    let last_pc = pc_arr[nb - 1]
        .as_ref()
        .ok_or_else(|| format!("PlateOrientation: Curves2d({nb}) is null"))?;
    let (fpar, lpar) = (last_pc.first_parameter(), last_pc.last_parameter());
    let uv = last_pc.d0((fpar + lpar) / 2.0);
    let mut pp1 = plate_surf.d0(uv.x(), uv.y());
    let uv = last_pc.d0(lpar);
    let mut pp2 = plate_surf.d0(uv.x(), uv.y());

    let mut sum_scal1 = 0.0f64;
    let mut sum_scal2 = 0.0f64;

    // `for (i = 1; i <= nb; i++)` (`cxx:1200-1221`).
    for i in 1..=nb {
        let a_pc = pc_arr[i - 1]
            .as_ref()
            .ok_or_else(|| format!("PlateOrientation: Curves2d({i}) is null"))?;
        let fpar = a_pc.first_parameter();
        let lpar = a_pc.last_parameter();

        // `aPC->D0(fpar, uv); thePlateSurf->D1(uv, pp2, du, dv); n1 = du ^ dv;
        //  n1.Normalize();` (`cxx:1205-1208`).
        let uv = a_pc.d0(fpar);
        let (p, du, dv) = plate_surf.d1(uv.x(), uv.y());
        pp2 = p;
        let mut n1 = du.crossed(&dv);
        n1.normalize();

        // `aPC->D0((fpar + lpar) / 2., uv); thePlateSurf->D0(uv, pp3);
        //  vv1(pp2, pp1); vv2(pp2, pp3); n2 = vv2 ^ vv1; n2.Normalize();
        //  SumScal1 += n1 * n2; SumScal2 += n2 * theRefDir; pp1 = pp3;`
        // (`cxx:1210-1220`).
        let uv = a_pc.d0((fpar + lpar) / 2.0);
        let pp3 = plate_surf.d0(uv.x(), uv.y());
        let vv1 = GpVec::from_pnts(&pp2, &pp1);
        let vv2 = GpVec::from_pnts(&pp2, &pp3);
        let mut n2 = vv2.crossed(&vv1);
        n2.normalize();

        sum_scal1 += n1.dot(&n2);
        sum_scal2 += n2.dot(ref_dir);
        pp1 = pp3;
    }

    // `cxx:1222-1228`.
    Ok(if sum_scal2 * sum_scal1 > 0.0 {
        Orientation::Forward
    } else {
        Orientation::Reversed
    })
}
