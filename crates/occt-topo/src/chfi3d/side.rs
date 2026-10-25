//! `ChFi3d::ConcaveSide`, `ChFi3d::NextSide`, `ChFi3d::SameSide`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d.cxx:300-647`, declared in
//! `ChFi3d.hxx:60-89`.

use occt_core::gp::GpVec;

use crate::abs::Orientation;
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_tool::BRepTool;
use crate::chfi3d::coefficient::chfi3d_coefficient;
use crate::shape::Edge;
use crate::topo_tools::{explorer_shapes, is_same};
use crate::topopebrep::tool::adaptor::{CurveAdaptor, SurfaceAdaptor};

use std::f64::consts::PI;

/// `ChFi3d::ConcaveSide(S1, S2, E, Or1, Or2)` (`cxx:305-535`).
///
/// Returns the "choice of fillet" code, `0` when the edge occurrence is not
/// found in either face (`cxx:355`, `:371`) and `10` for the face that ends in
/// the extension of one of the two base faces (`cxx:464`).
///
/// `Err` carries the null-handle dereferences OCCT performs unconditionally
/// (`cxx:313`, `:379-380`, `:384-385`).
pub fn concave_side(
    s1: &SurfaceAdaptor,
    s2: &SurfaceAdaptor,
    e: &Edge,
    or1: &mut Orientation,
    or2: &mut Orientation,
) -> Result<i32, String> {
    *or1 = Orientation::Forward;
    *or2 = Orientation::Forward;
    let ce = CurveAdaptor::of_edge(e);
    let first = ce.first_parameter();
    let last = ce.last_parameter();
    let par = 0.691254 * first + 0.308746 * last;

    let f1 = s1
        .face()
        .ok_or_else(|| "ChFi3d::ConcaveSide: S1 without face".to_string())?
        .clone();
    let f2 = s2
        .face()
        .ok_or_else(|| "ChFi3d::ConcaveSide: S2 without face".to_string())?
        .clone();

    let (_, mut tge) = ce
        .d1(par)
        .ok_or_else(|| "ChFi3d::ConcaveSide: no edge curve".to_string())?;
    tge.normalize();
    let mut tge1 = tge;
    let mut tge2 = tge;
    if e.0.orientation() == Orientation::Reversed {
        tge.reverse();
    }

    let mut e1 = e.clone();
    let mut e2 = e.clone();
    e1.0.set_orientation(Orientation::Forward);
    e2.0.set_orientation(Orientation::Forward);

    if is_same(&f1.0, &f2.0) && BRepTool::is_closed_edge_face(e, &f1) {
        e2.0.set_orientation(Orientation::Reversed);
        tge2.reverse();
    } else {
        let mut found = false;
        for occ in explorer_shapes(&f1.0, crate::abs::ShapeType::Edge) {
            if is_same(&occ, &e.0) {
                if occ.orientation() == Orientation::Reversed {
                    tge1.reverse();
                }
                found = true;
                break;
            }
        }
        if !found {
            return Ok(0);
        }
        found = false;
        for occ in explorer_shapes(&f2.0, crate::abs::ShapeType::Edge) {
            if is_same(&occ, &e.0) {
                if occ.orientation() == Orientation::Reversed {
                    tge2.reverse();
                }
                found = true;
                break;
            }
        }
        if !found {
            return Ok(0);
        }
    }

    let pc1 = Curve2dAdaptor::of_edge_face(&e1, &f1);
    let pc2 = Curve2dAdaptor::of_edge_face(&e2, &f2);
    let mut p2d1 = pc1
        .d0(par)
        .ok_or_else(|| "ChFi3d::ConcaveSide: no pcurve on F1".to_string())?;
    let mut p2d2 = pc2
        .d0(par)
        .ok_or_else(|| "ChFi3d::ConcaveSide: no pcurve on F2".to_string())?;

    let (pt1, mut du1, mut dv1) = s1
        .d1(p2d1.x(), p2d1.y())
        .ok_or_else(|| "ChFi3d::ConcaveSide: S1->D1 failed".to_string())?;
    let mut ns1 = du1.crossed(&dv1);
    ns1.normalize();
    if f1.0.orientation() == Orientation::Reversed {
        ns1.reverse();
    }
    let (pt2, mut du2, mut dv2) = s2
        .d1(p2d2.x(), p2d2.y())
        .ok_or_else(|| "ChFi3d::ConcaveSide: S2->D1 failed".to_string())?;
    let mut ns2 = du2.crossed(&dv2);
    ns2.normalize();
    if f2.0.orientation() == Orientation::Reversed {
        ns2.reverse();
    }

    let mut dint1 = ns1.crossed(&tge1);
    let mut dint2 = ns2.crossed(&tge2);
    let mut ang = ns1.cross_magnitude(&ns2);
    if ang > 0.0001 * PI {
        let scal = ns2.dot(&dint1);
        if scal <= 0.0 {
            ns2.reverse();
            *or2 = Orientation::Reversed;
        }
        let scal = ns1.dot(&dint2);
        if scal <= 0.0 {
            ns1.reverse();
            *or1 = Orientation::Reversed;
        }
    } else if dint1.dot(&dint2) < 0.0 {
        // This is a forgotten regularity (`cxx:420-462`).
        let (_, du1b, dv1b, ddu, ddv, _) = s1
            .d2(p2d1.x(), p2d1.y())
            .ok_or_else(|| "ChFi3d::ConcaveSide: S1->D2 failed".to_string())?;
        du1 = du1b;
        dv1 = dv1b;
        let su = if du1.dot(&dint1) < 0.0 { ddu.reversed() } else { ddu };
        let sv = if dv1.dot(&dint1) < 0.0 { ddv.reversed() } else { ddv };
        du1 = du1.added(&su);
        dv1 = dv1.added(&sv);
        ns1 = du1.crossed(&dv1);
        ns1.normalize();
        if f1.0.orientation() == Orientation::Reversed {
            ns1.reverse();
        }

        let (_, du2b, dv2b, ddu, ddv, _) = s2
            .d2(p2d2.x(), p2d2.y())
            .ok_or_else(|| "ChFi3d::ConcaveSide: S2->D2 failed".to_string())?;
        du2 = du2b;
        dv2 = dv2b;
        let su = if du2.dot(&dint2) < 0.0 { ddu.reversed() } else { ddu };
        let sv = if dv2.dot(&dint2) < 0.0 { ddv.reversed() } else { ddv };
        du2 = du2.added(&su);
        dv2 = dv2.added(&sv);
        ns2 = du2.crossed(&dv2);
        ns2.normalize();
        if f2.0.orientation() == Orientation::Reversed {
            ns2.reverse();
        }

        dint1 = ns1.crossed(&tge1);
        dint2 = ns2.crossed(&tge2);
        ang = ns1.cross_magnitude(&ns2);
        if ang > 0.0001 * PI {
            let scal = ns2.dot(&dint1);
            if scal <= 0.0 {
                ns2.reverse();
                *or2 = Orientation::Reversed;
            }
            let scal = ns1.dot(&dint2);
            if scal <= 0.0 {
                ns1.reverse();
                *or1 = Orientation::Reversed;
            }
        } else {
            // This 10 shows that the face at end is in the extension of one of
            // the two base faces (`cxx:463-465`).
            return Ok(10);
        }
    } else {
        // Here it turns back, the points are taken in faces neither too close
        // nor too far as much as possible (`cxx:469-503`).
        let (u, v) = chfi3d_coefficient(&dint1, &du1, &dv1);
        p2d1.set_x(p2d1.x() + u);
        p2d1.set_y(p2d1.y() + v);
        let (u, v) = chfi3d_coefficient(&dint1, &du2, &dv2);
        p2d2.set_x(p2d2.x() + u);
        p2d2.set_y(p2d2.y() + v);
        let (pt1b, du1b, dv1b) = s1
            .d1(p2d1.x(), p2d1.y())
            .ok_or_else(|| "ChFi3d::ConcaveSide: S1->D1 failed".to_string())?;
        ns1 = du1b.crossed(&dv1b);
        if f1.0.orientation() == Orientation::Reversed {
            ns1.reverse();
        }
        let (pt2b, du2b, dv2b) = s2
            .d1(p2d2.x(), p2d2.y())
            .ok_or_else(|| "ChFi3d::ConcaveSide: S2->D1 failed".to_string())?;
        ns2 = du2b.crossed(&dv2b);
        if f2.0.orientation() == Orientation::Reversed {
            ns2.reverse();
        }
        let vref = GpVec::from_pnts(&pt1b, &pt2b);
        if ns1.dot(&vref) < 0.0 {
            *or1 = Orientation::Reversed;
        }
        if ns2.dot(&vref) > 0.0 {
            *or2 = Orientation::Reversed;
        }
        let _ = (pt1, pt2);
    }

    let mut choix_conge = match (*or1, *or2) {
        (Orientation::Forward, Orientation::Forward) => 1,
        (Orientation::Forward, _) => 7,
        (_, Orientation::Forward) => 3,
        _ => 5,
    };
    if ns1.crossed(&ns2).dot(&tge) >= 0.0 {
        choix_conge += 1;
    }
    Ok(choix_conge)
}

/// `ChFi3d::NextSide(Or1, Or2, OrSave1, OrSave2, ChoixSave)` (`cxx:537-601`).
pub fn next_side(
    or1: &mut Orientation,
    or2: &mut Orientation,
    or_save1: Orientation,
    or_save2: Orientation,
    choix_save: i32,
) -> i32 {
    *or1 = if *or1 == Orientation::Forward {
        or_save1
    } else {
        or_save1.reversed()
    };
    *or2 = if *or2 == Orientation::Forward {
        or_save2
    } else {
        or_save2.reversed()
    };

    let mut choix_conge = match (*or1, *or2) {
        (Orientation::Forward, Orientation::Forward) => 1,
        (Orientation::Forward, _) => {
            if choix_save < 0 {
                3
            } else {
                7
            }
        }
        (_, Orientation::Forward) => {
            if choix_save < 0 {
                7
            } else {
                3
            }
        }
        _ => 5,
    };
    if choix_save.abs() % 2 == 0 {
        choix_conge += 1;
    }
    choix_conge
}

/// `ChFi3d::NextSide(Or, OrSave, OrFace)` (`cxx:605-616`): the single-edge
/// form.
pub fn next_side_single(or: &mut Orientation, or_save: Orientation, or_face: Orientation) {
    *or = if *or == or_face {
        or_save
    } else {
        or_save.reversed()
    };
}

/// `ChFi3d::SameSide(Or, OrSave1, OrSave2, OrFace1, OrFace2)` (`cxx:622-647`).
pub fn same_side(
    or: Orientation,
    or_save1: Orientation,
    or_save2: Orientation,
    or_face1: Orientation,
    or_face2: Orientation,
) -> bool {
    let o1 = if or == or_face1 {
        or_save1
    } else {
        or_save1.reversed()
    };
    let o2 = if or == or_face2 {
        or_save2
    } else {
        or_save2.reversed()
    };
    o1 == o2
}
