//! `BRepFeat_RibSlot::ChoiceOfFaces` and `HeightMax`
//! (`ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_RibSlot.cxx:645-740`).
//!
//! `ChoiceOfFaces` picks the support face of a profile supported on an edge:
//! eight lines rotated by `i * PI / 9` around the profile tangent are cast
//! against every candidate face, and the face whose first hit has the smallest
//! non-negative parameter wins. `HeightMax` sizes the tool from the bounding
//! box of the basis and the limit.

use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::bnd::BndBox;
use occt_core::gp::{GpAx1, GpPln, GpPnt, GpTrsf};
use occt_geom::curve::Curve;
use occt_geom::line::GeomLine;

use crate::brep_bnd_lib;
use crate::locope::CSIntersector;
use crate::shape::{Face, TopoShape};

/// `ChoiceOfFaces(faces, cc, par, bnd, Pln)` (`cxx:645-704`).
///
/// `Ok(None)` is the null `TopoDS_Face` the reference returns when no candidate
/// is hit (`cxx:703`); a negative first parameter never wins (`cxx:695-697`),
/// so the answer is the hit closest to `cc(par)`. An `Err` is the
/// `Standard_ConstructionError` of `LocOpe_CSIntersector::Perform`.
pub fn choice_of_faces(
    faces: &[TopoShape],
    cc: &Arc<dyn Curve>,
    par: f64,
    _bound: f64,
    pln: &GpPln,
) -> Result<Option<Face>, String> {
    let (pp, tgt) = cc.d1(par);

    // `occ::handle<Geom_Line> l1 = new Geom_Line(pp, tgt);` (`cxx:663`).
    // `gp_Dir(tgt)` raises on a null tangent.
    let dir = occt_core::gp::GpDir::from_vec(&tgt)
        .map_err(|_| "Standard_ConstructionError: null profile tangent".to_string())?;
    let l1 = GeomLine::new(occt_core::gp::GpLin::from_pnt_dir(pp, dir));

    // `gp_Ax1 Axe(pp, Pln->Position().Direction())` and the eight rotated
    // `l1->Rotated(Axe, i * M_PI / 9.)` (`cxx:669-677`).
    let axe = GpAx1::new(pp, pln.position().direction());
    let mut scur: Vec<Option<Arc<dyn Curve>>> = Vec::with_capacity(8);
    for i in 1..=8 {
        let mut trsf = GpTrsf::identity();
        trsf.set_rotation_ax1(&axe, i as f64 * PI / 9.0)
            .map_err(|e| e.to_string())?;
        scur.push(Some(Arc::from(l1.transformed(&trsf))));
    }
    let counter = scur.len();

    // `double Par = RealLast();` (`cxx:679`).
    let mut par_best = f64::MAX;
    let mut chosen: Option<Face> = None;
    for f in faces {
        // `TopoDS::Face(it.Value())` (`cxx:684`).
        let Some(face) = Face::wrap(f.clone()) else {
            continue;
        };
        let mut asi = CSIntersector::with_shape(&face.0);
        asi.perform_curves(&scur)?;
        if !asi.is_done() {
            continue;
        }
        // `for (int jj = 1; jj <= Counter; jj++)` (`cxx:689-701`).
        for jj in 1..=counter {
            if asi.nb_points(jj).unwrap_or(0) >= 1 {
                let app = asi.point(jj, 1).map(|p| p.parameter()).unwrap_or(f64::MAX);
                if app >= 0.0 && app < par_best {
                    par_best = app;
                    chosen = Some(face.clone());
                }
            }
        }
    }
    Ok(chosen)
}

/// `HeightMax(theSbase, theSUntil, p1, p2)` (`cxx:712-740`): the two corners of
/// a box twice the size of the `theSbase` / `theSUntil` bounds, and the largest
/// absolute coordinate of that bounds.
///
/// The reference fills `c[0..5]` with
/// `Box.Get(c[0], c[2], c[4], c[1], c[3], c[5])` (`cxx:727`), i.e. the slots are
/// `xmin, xmax, ymin, ymax, zmin, zmax`, and then builds `p1` from `c[0..2]` and
/// `p2` from `c[3..5]` (`cxx:735-736`). The port keeps that slot layout: it is
/// what the reference computes, not a sensible corner pair.
pub fn height_max(sbase: &TopoShape, suntil: Option<&TopoShape>) -> (GpPnt, GpPnt, f64) {
    let mut b = BndBox::new();
    brep_bnd_lib::add(sbase, &mut b, true);
    if let Some(s) = suntil {
        brep_bnd_lib::add(s, &mut b, true);
    }
    let Some((xmin, ymin, zmin, xmax, ymax, zmax)) = b.get() else {
        // `Bnd_Box::Get` raises on a void box in the reference.
        return (GpPnt::default(), GpPnt::default(), 0.0);
    };
    let c = [xmin, xmax, ymin, ymax, zmin, zmax];
    let bnd = c.iter().fold(c[0], |acc, &v| if v > acc { v } else { acc });
    let p1 = GpPnt::new(c[0] - 2.0 * bnd, c[1] - 2.0 * bnd, c[2] - 2.0 * bnd);
    let p2 = GpPnt::new(c[3] + 2.0 * bnd, c[4] + 2.0 * bnd, c[5] + 2.0 * bnd);
    (p1, p2, bnd)
}
