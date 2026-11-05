//! The three free functions at the top of `BRepFill_Draft.cxx`:
//! `ComputeTrsf` (`cxx:68-104`), `Longueur` (`cxx:108-135`) and
//! `GoodOrientation` (`cxx:140-190`).

use occt_core::bnd::BndBox;
use occt_core::gp::{GpAx3, GpDir, GpPnt, GpTrsf, GpVec, GpXyz};

use crate::brep_fill_location_law::LocationLawBase;
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::shape::Wire;

/// `Bnd_Box::Get(Xmin, Ymin, Zmin, Xmax, Ymax, Zmax)` (`Bnd_Box.cxx`
/// `Get`): the six corners in the reference's argument order.
///
/// [`BndBox::get`] reports `(xmin, xmax, ymin, ymax, zmin, zmax)` and `None`
/// for a void box; the reference's `Get` would report `RealFirst()` /
/// `RealLast()` on an empty box, so a void box is an error here.
pub(super) fn bounds(b: &BndBox) -> Result<(f64, f64, f64, f64, f64, f64), String> {
    let (xmin, xmax, ymin, ymax, zmin, zmax) = b
        .get()
        .ok_or_else(|| "Bnd_Box::Get on an empty box".to_string())?;
    Ok((xmin, ymin, zmin, xmax, ymax, zmax))
}

/// `ComputeTrsf(W, D, Box, Tf)` (`cxx:68-104`): the frame whose `Z` axis is `D`
/// and whose origin is the wire's vertex barycenter, plus the wire's box in
/// that frame.
pub(super) fn compute_trsf(w: &Wire, d: &GpDir) -> Result<(BndBox, GpTrsf), String> {
    // `cxx:71-84`: the barycenter of the explorer's current vertices.
    let mut exp = WireExplorer::of_wire(w);
    let mut bary = GpXyz::new(0.0, 0.0, 0.0);
    let mut nb = 0usize;
    while exp.more() {
        if let Some(v) = exp.current_vertex() {
            bary = bary.added(BRepTool::vertex_point(v).xyz());
        }
        nb += 1;
        exp.next();
    }
    // `cxx:83`: `Bary /= nb`. A wire without an edge divides by zero, exactly
    // as the reference does.
    let bary = bary.divided(nb as f64);

    // `cxx:87-88`: `gp_Ax3 N(Bary, D); Tf.SetTransformation(N);`
    let n = GpAx3::from_point_dir(GpPnt::from_xyz(&bary), *d)
        .map_err(|e| format!("gp_Ax3: {e}"))?;
    let mut tf = GpTrsf::identity();
    tf.set_transformation(&n);

    // `cxx:92-95`: `TheW.Location(Loc)`.
    let the_w = crate::shape_ops::transformed_copy(&w.0, &tf)?;

    // `cxx:97-102`: `BndLib_Add3dCurve::Add(AC, 0.1, Box)` per edge.
    let mut box_ = BndBox::new();
    box_.set_void();
    for e in crate::topo_tools::edges_of(&the_w) {
        let Some(curve) = BRepTool::edge_curve_world(&e) else {
            continue;
        };
        let (u1, u2) = BRepTool::edge_parameters(&e);
        let b = crate::geom_bnd_lib_curve3d::box_curve(curve.as_ref(), u1, u2, 0.1);
        box_.add_box(&b);
    }
    Ok((box_, tf))
}

/// `Longueur(WBox, SBox, D, P)` (`cxx:108-135`): the rule length from the
/// wire's box to the stopping box, flipping `D` when the skin would run the
/// other way.
pub(super) fn longueur(w_box: &BndBox, s_box: &BndBox, d: &mut GpDir, p: &mut GpPnt) -> f64 {
    let (_, _, wzmin, _, _, wzmax) = bounds(w_box).expect("wire box is filled by ComputeTrsf");
    let (sxmin, symin, szmin, sxmax, symax, szmax) =
        bounds(s_box).expect("stop box is filled by the caller");

    // `cxx:120`.
    *p = GpPnt::new((sxmin + sxmax) / 2.0, (symin + symax) / 2.0, szmax);

    if szmax < wzmin {
        // `cxx:123-127`: the skin points the wrong way.
        d.reverse();
        p.set_z(szmin);
        wzmax - szmin
    } else {
        // `cxx:129-132`.
        szmax - wzmin
    }
}

/// `GoodOrientation(B, Law, D)` (`cxx:140-190`): `true` when the law's frames
/// turn counter-clockwise around `D` sampled along the path, i.e. when they
/// give an exterior skin.
pub(super) fn good_orientation(
    b: &BndBox,
    law: &LocationLawBase,
    d: &GpDir,
) -> Result<bool, String> {
    // `cxx:145-150`.
    let (a_xmin, a_ymin, a_zmin, a_xmax, a_ymax, a_zmax) = bounds(b)?;
    let p1 = GpPnt::new(a_xmin, a_ymin, a_zmin);
    let p2 = GpPnt::new(a_xmax, a_ymax, a_zmax);
    let v = GpVec::from_xyz(&p2.xyz().subtracted(p1.xyz()));

    // `cxx:151-153`: `Law->CurvilinearBounds(Law->NbLaw(), f, l); r = V.Magnitude() / l;`
    let (_f, l) = law.curvilinear_bounds(law.nb_law())?;
    let mut r = v.magnitude() / l;

    // `cxx:161-164`: `Nb = (int)(4 + 10 * r); r = l / Nb; Nb++;`
    let n0 = (4.0 + 10.0 * r) as i32;
    r = l / n0 as f64;
    let nb = (n0 + 1) as usize;

    // `cxx:166-178`: the sampled points and their barycenter.
    let mut pnts: Vec<GpPnt> = Vec::with_capacity(nb);
    let mut bary = GpXyz::new(0.0, 0.0, 0.0);
    for ii in 1..=nb {
        let (ind, t) = law.parameter((ii - 1) as f64 * r)?;
        let curve = law
            .law(ind)
            .ok_or_else(|| "Standard_NullObject: BRepFill_DraftLaw::Law".to_string())?
            .get_curve()
            .ok_or_else(|| "Standard_NullObject: BRepFill_DraftLaw::GetCurve".to_string())?;
        let pnt = curve.d0(t);
        bary = bary.added(pnt.xyz());
        pnts.push(pnt);
    }
    let center = GpPnt::from_xyz(&bary.divided(nb as f64));

    // `cxx:178-187`.
    let normal = GpVec::from_xyz(d.xyz());
    let mut angle = 0.0;
    let mut reference = GpVec::from_xyz(&pnts[0].xyz().subtracted(center.xyz()));
    for pnt in pnts.iter().skip(1) {
        let current = GpVec::from_xyz(&pnt.xyz().subtracted(center.xyz()));
        angle += reference.angle_with_ref(&current, &normal);
        reference = current;
    }

    // `cxx:189`: `return (Angle >= 0);`
    Ok(angle >= 0.0)
}
