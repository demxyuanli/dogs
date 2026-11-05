//! `TrimEdge` (`BRepFill_OffsetWire.cxx:2379-2570`): cuts the parallel edge
//! `E` at the vertices collected from the bisectors, dropping the ones the
//! detromp list rejects.

use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_core::gp::GpPnt;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_fill_offset_wire::maps::ShapeToShapeMap;
use crate::brep_fill_offset_wire::vertices::{
    double_or_not_inside, is_inner_edge, store_in_map,
};

use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::edge_vertices;

/// `TrimEdge(E, ProE, AllSpine, Detromp, TheVer, ThePar, S, MapVV, IndOfE)`
/// (`cxx:2379-2570`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn trim_edge(
    e: &Edge,
    pro_e: &TopoShape,
    all_spine: &Face,
    detromp: &[TopoShape],
    the_ver: &mut Vec<Vertex>,
    the_par: &mut Vec<f64>,
    s: &mut Vec<TopoShape>,
    map_vv: &mut ShapeToShapeMap,
    ind_of_e: i32,
) {
    s.clear();

    // `cxx:2396-2408`: keep `ThePar` sorted, moving `TheVer` alongside.
    let mut changed = true;
    while changed {
        changed = false;
        for i in 0..the_par.len().saturating_sub(1) {
            if the_par[i] > the_par[i + 1] {
                the_par.swap(i, i + 1);
                the_ver.swap(i, i + 1);
                changed = true;
            }
        }
    }

    let degenerated = GeometryRegistry::global().is_degenerated_edge(&e.0);

    // `cxx:2413-2426`: drop the vertices the detromp list rejects.
    if !degenerated {
        let mut k = 0;
        while k < the_ver.len() {
            if double_or_not_inside(detromp, &the_ver[k]) {
                the_ver.remove(k);
                the_par.remove(k);
            } else {
                k += 1;
            }
        }
    }

    // `cxx:2431-2467`: one representation per coincident vertex pair.
    if !degenerated {
        let par_tol = 2.0 * PCONFUSION;
        let mut k = 0;
        while k + 1 < the_ver.len() {
            if crate::topo_tools::is_same(&the_ver[k].0, &the_ver[k + 1].0)
                || (the_par[k] - the_par[k + 1]).abs() <= par_tol
            {
                if k + 1 == the_ver.len() - 1 {
                    let (a, b) = (the_ver[k].0.clone(), the_ver[k + 1].0.clone());
                    store_in_map(&a, &b, map_vv);
                    the_ver.remove(k);
                    the_par.remove(k);
                } else {
                    let (a, b) = (the_ver[k + 1].0.clone(), the_ver[k].0.clone());
                    store_in_map(&a, &b, map_vv);
                    the_ver.remove(k + 1);
                    the_par.remove(k + 1);
                }
                k = k.saturating_sub(1);
            } else {
                k += 1;
            }
        }
    }

    // `cxx:2469-2547`: the extreme edge of an open result.
    if ind_of_e == 1 || ind_of_e == -1 {
        let mut new_edge = e.0.empty_copied();
        if let Some(geom) = GeometryRegistry::global().edge_geom(&e.0) {
            GeometryRegistry::global().set_edge(&new_edge, geom);
        }
        let (v1, v2) = edge_vertices(e);
        let mut v1 = v1.map(|v| v.0);
        let mut v2 = v2.map(|v| v.0);
        let Some((pcurve, mut fpar, mut lpar)) = curve_on_surface_range(e, all_spine) else {
            panic!("BRep_Tool::CurveOnSurface(E, aPlane, aLoc, fpar, lpar) raised");
        };
        let trimmed = is_inner_edge(pro_e, all_spine);
        let (tr_par1, tr_par2) = match trimmed {
            Some((a, b)) => (a, b),
            None => (fpar, lpar),
        };
        let to_trim_as_origin = trimmed.is_some();

        let builder = TopoBuilder::new();
        let forward = new_edge.orientation() != Orientation::Reversed;
        if ind_of_e == 1 {
            // `cxx:2487-2501`.
            if forward {
                if to_trim_as_origin {
                    fpar = tr_par1;
                    let tr = pcurve.d0(fpar);
                    v1 = Some(
                        builder
                            .make_vertex(GpPnt::new(tr.x(), tr.y(), 0.0), CONFUSION)
                            .0,
                    );
                }
                if let Some(v) = &v1 {
                    builder.add(&mut new_edge, &v.oriented(Orientation::Forward));
                }
                if !the_ver.is_empty() {
                    builder.add(&mut new_edge, &the_ver[0].0.oriented(Orientation::Reversed));
                }
                GeometryRegistry::global().set_ranges_all(
                    &new_edge,
                    fpar,
                    the_par.first().copied().unwrap_or(fpar),
                );
            } else {
                if to_trim_as_origin {
                    lpar = tr_par2;
                    let tr = pcurve.d0(lpar);
                    v2 = Some(
                        builder
                            .make_vertex(GpPnt::new(tr.x(), tr.y(), 0.0), CONFUSION)
                            .0,
                    );
                }
                if !the_ver.is_empty() {
                    builder.add(&mut new_edge, &the_ver[0].0.oriented(Orientation::Reversed));
                }
                if let Some(v) = &v2 {
                    builder.add(&mut new_edge, &v.oriented(Orientation::Forward));
                }
                GeometryRegistry::global().set_ranges_all(
                    &new_edge,
                    the_par.first().copied().unwrap_or(lpar),
                    lpar,
                );
            }
        } else if forward {
            // `cxx:2514-2528` (last edge of an open wire, FORWARD).
            if to_trim_as_origin {
                lpar = tr_par2;
                let tr = pcurve.d0(lpar);
                v2 = Some(
                    builder
                        .make_vertex(GpPnt::new(tr.x(), tr.y(), 0.0), CONFUSION)
                        .0,
                );
            }
            if !the_ver.is_empty() {
                builder.add(&mut new_edge, &the_ver[0].0.oriented(Orientation::Forward));
            }
            if let Some(v) = &v2 {
                builder.add(&mut new_edge, &v.oriented(Orientation::Reversed));
            }
            GeometryRegistry::global().set_ranges_all(
                &new_edge,
                the_par.first().copied().unwrap_or(fpar),
                lpar,
            );
        } else {
            // `cxx:2531-2543`.
            if to_trim_as_origin {
                fpar = tr_par1;
                let tr = pcurve.d0(fpar);
                v1 = Some(
                    builder
                        .make_vertex(GpPnt::new(tr.x(), tr.y(), 0.0), CONFUSION)
                        .0,
                );
            }
            if let Some(v) = &v1 {
                builder.add(&mut new_edge, &v.oriented(Orientation::Reversed));
            }
            if !the_ver.is_empty() {
                builder.add(&mut new_edge, &the_ver[0].0.oriented(Orientation::Forward));
            }
            GeometryRegistry::global().set_ranges_all(
                &new_edge,
                fpar,
                the_par.first().copied().unwrap_or(lpar),
            );
        }
        s.push(new_edge);
        return;
    }

    // `cxx:2549-2569`: consecutive pairs become the trimmed edges.
    let builder = TopoBuilder::new();
    let mut k = 0;
    while k + 1 < the_ver.len() {
        let mut new_edge = e.0.empty_copied();
        if let Some(geom) = GeometryRegistry::global().edge_geom(&e.0) {
            GeometryRegistry::global().set_edge(&new_edge, geom);
        }
        let reversed = new_edge.orientation() == Orientation::Reversed;
        if reversed {
            builder.add(&mut new_edge, &the_ver[k].0.oriented(Orientation::Reversed));
            builder.add(&mut new_edge, &the_ver[k + 1].0.oriented(Orientation::Forward));
        } else {
            builder.add(&mut new_edge, &the_ver[k].0.oriented(Orientation::Forward));
            builder.add(&mut new_edge, &the_ver[k + 1].0.oriented(Orientation::Reversed));
        }
        GeometryRegistry::global().set_ranges_all(&new_edge, the_par[k], the_par[k + 1]);
        s.push(new_edge);
        k += 2;
    }
}
