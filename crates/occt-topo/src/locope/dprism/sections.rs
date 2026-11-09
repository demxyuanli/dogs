//! `LocOpe_DPrism::Curves` (`cxx:620-690`) and `BarycCurve` (`cxx:694-753`).
//!
//! Both read the sweep direction out of `myProfile2` (the `(dy, dz)` of its two
//! endpoints) and the normal of the planar spine, then:
//!
//! - `Curves`: one line per sampled point of every non-degenerated spine edge,
//!   through the point, along `dy * (Normal ^ d1) + dz * Normal`, trimmed to
//!   `[-2|Height|, 2|Height|]` (`cxx:657-689`).
//! - `BarycCurve`: a single line through the barycentre of
//!   `LocOpe::SampleEdges(myFirstShape)` (or of the spine when that cap is
//!   null), along `dz * Normal` (`cxx:734-752`).
//!
//! `NECHANT` is the reference's `#define NECHANT 7` (`cxx:57`).
//!
//! Deviations, all in the "`gp` raises" family already recorded in this tree:
//! a null cross product (`gp_Dir` construction, `cxx:681`), a null direction or
//! a null line vector (`gp_Vec::Normalize`, `gp_Lin` construction) skips that
//! sample instead of throwing, and an empty sample list answers `None` where
//! the reference divides by zero (`cxx:749`).

use std::sync::Arc;

use occt_core::gp::{GpDir, GpLin, GpPnt, GpVec};
use occt_geom::curve::Curve;
use occt_geom::line::GeomLine;
use occt_geom::trimmed::GeomTrimmedCurveBasis;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::locope::dprism::DPrism;
use crate::locope::sample_edges::sample_edges;
use crate::locope::shape_maps::ShapeSet;
use crate::shape::Edge;
use crate::topo_tools::{edge_vertices, explorer_shapes};

/// `#define NECHANT 7` (`cxx:57`).
const NECHANT: usize = 7;

/// `TopExp::Vertices(myProfile2, V1, V2)` followed by `BRep_Tool::Pnt`
/// (`cxx:623-627`, `cxx:697-699`): the two endpoint points of the middle
/// profile edge, i.e. the `(dy, dz)` the reference reads out of it.
fn profile2_points(profile2: &Edge) -> Option<(GpPnt, GpPnt)> {
    let (v1, v2) = edge_vertices(profile2);
    let (Some(v1), Some(v2)) = (v1, v2) else {
        return None;
    };
    Some((BRepTool::vertex_point(&v1), BRepTool::vertex_point(&v2)))
}

/// `Curves(Scurves)` (`cxx:620-690`).
pub(super) fn curves(m: &DPrism) -> Vec<Arc<dyn Curve>> {
    let mut scurves: Vec<Arc<dyn Curve>> = Vec::new();

    // `cxx:621-627`: the `(dy, dz)` of `myProfile2`.
    let Some(profile2) = m.my_profile2.clone() else {
        return scurves;
    };
    let Some((p1, p2)) = profile2_points(&profile2) else {
        return scurves;
    };
    let dy = p2.y() - p1.y();
    let dz = p2.z() - p1.z();

    // `cxx:630-648`: the normal of the spine's plane.
    let Some(surface) = BRepTool::face_surface(&m.my_spine) else {
        return scurves;
    };
    let surface = surface.rectangular_trimmed_basis().unwrap_or(surface);
    // `occ::down_cast<Geom_Plane>(S).IsNull()` -> `Standard_ConstructionError`.
    let Some(plane) = surface.gp_pln() else {
        return scurves;
    };
    let axis = plane.axis();
    let mut normale = axis.direction().clone();
    if !plane.is_direct() {
        normale.reverse();
    }

    // `cxx:649-657`.
    let mut the_map = ShapeSet::new();
    for e in explorer_shapes(&m.my_spine.0.oriented(Orientation::Forward), ShapeType::Edge) {
        let edg = Edge(e);
        // `if (!theMap.Add(edg)) continue;`
        if !the_map.add(&edg.0) {
            continue;
        }
        if BRepTool::is_degenerated(&edg) {
            continue;
        }

        // `cxx:665-666`: `BRep_Tool::Curve` followed by `Transformed(Loc)`.
        let Some(c) = BRepTool::edge_curve_world(&edg) else {
            continue;
        };
        let (f, l) = BRepTool::edge_parameters(&edg);

        // `cxx:668-669`.
        let u1 = -2.0 * m.my_height.abs();
        let u2 = 2.0 * m.my_height.abs();

        // `cxx:671-688`.
        for i in 0..=NECHANT {
            let prm = ((NECHANT - i) as f64 * f + i as f64 * l) / NECHANT as f64;
            let (pt, mut d1) = c.d1(prm);
            if edg.0.orientation() == Orientation::Reversed {
                d1.reverse();
            }
            // `d1.Normalize()`.
            let Ok(d1) = GpDir::from_vec(&d1) else {
                continue;
            };
            // `gp_Dir locy = Normale.Crossed(d1)`.
            let Ok(locy) = normale.crossed(&d1) else {
                continue;
            };
            // `gp_Vec ldir = dy * locy.XYZ() + dz * Normale.XYZ()`.
            let ldir = GpVec::from_xyz(locy.xyz())
                .multiplied_scalar(dy)
                .added(&GpVec::from_xyz(normale.xyz()).multiplied_scalar(dz));
            // `gp_Lin lin(pt, ldir)`.
            let Ok(dir) = GpDir::from_vec(&ldir) else {
                continue;
            };
            let the_lin: Arc<dyn Curve> = Arc::new(GeomLine::new(GpLin::from_pnt_dir(pt, dir)));
            // `Geom_TrimmedCurve(Lin, u1, u2, true)`.
            if let Ok(trlin) = GeomTrimmedCurveBasis::with_sense(the_lin, u1, u2, true) {
                scurves.push(Arc::new(trlin));
            }
        }
    }

    scurves
}

/// `BarycCurve()` (`cxx:694-753`).
pub(super) fn baryc_curve(m: &DPrism) -> Option<GeomLine> {
    // `cxx:695-700`: `(dz)` of `myProfile2`.
    let profile2 = m.my_profile2.clone()?;
    let (p1, p2) = profile2_points(&profile2)?;
    let dz = p2.z() - p1.z();

    // `cxx:702-718`.
    let surface = BRepTool::face_surface(&m.my_spine)?;
    let surface = surface.rectangular_trimmed_basis().unwrap_or(surface);
    let plane = surface.gp_pln()?;
    let axis = plane.axis();
    let mut normale = axis.direction().clone();
    if !plane.is_direct() {
        normale.reverse();
    }

    // `cxx:720-732`: the reversed-spine branch only traces (`cts20871`), the
    // `Normale.Reverse()` inside it is commented out in the reference.

    // `cxx:733`.
    let vec = GpVec::from_xyz(normale.xyz()).multiplied_scalar(dz);

    // `cxx:735-746`: the barycentre of the sampled cap edges.
    let spt = match &m.my_first_shape {
        Some(first) => sample_edges(first),
        None => sample_edges(&m.my_spine.0),
    };
    if spt.is_empty() {
        return None;
    }
    let mut bar = GpPnt::new(0.0, 0.0, 0.0);
    for p in &spt {
        let c = bar.xyz();
        bar = GpPnt::new(c.x + p.x(), c.y + p.y(), c.z + p.z());
    }
    let n = spt.len() as f64;
    let c = bar.xyz();
    let bar = GpPnt::new(c.x / n, c.y / n, c.z / n);

    // `cxx:750-752`: `gp_Ax1 newAx(bar, Vec)` then `Geom_Line`.
    let dir = GpDir::from_vec(&vec).ok()?;
    Some(GeomLine::new(GpLin::from_pnt_dir(bar, dir)))
}
