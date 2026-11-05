//! `ChFi3d_FilBuilder::SplitSurf` — split a fillet patch at the parameter
//! where the two sides of the band pinch together.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder.cxx:2274-2436`
//! (OCCT 8.0.0).
//!
//! OCCT declares the method on `ChFi3d_FilBuilder`, but the body reads only the
//! base class state (`myDS`, `tolapp3d`, `tol2d`) and the `SurfData` sequence,
//! so it is a method of [`ChFi3dBuilder`] here — the subclass' `PerformSurf`
//! bodies have no second handle on the base to call it through (see the
//! `FilSurfOps` note in `fil_builder/mod.rs`).
//!
//! The two halves of the body are kept statement by statement:
//!
//! 1. `(1)` scans the line for the parameters where the distance between
//!    `PointOnS1` and `PointOnS2` has a local minimum (`cxx:2299-2378`): a
//!    coarse scan over the line points, then a `math_FunctionRoot` on
//!    `ChFi3d_SearchSing` for each tripped minimum,
//! 2. `(2)` duplicates the first `SurfData` once per retained parameter and
//!    re-links the two halves at it (`cxx:2380-2435`).
//!
//! Divergences that are not heuristics:
//! - `ChFiDS_SurfData::Copy` is [`SurfData::copy_from`] (it leaves
//!   `indexOfC1` / `indexOfC2` / `isoncurv1` / `isoncurv2` alone, as OCCT
//!   does), so the duplicates start from `SurfData::new()` rather than from a
//!   Rust `clone`,
//! - `ref` is a handle in OCCT, i.e. the element of `SeqData` itself. The port
//!   edits a local copy and writes it back into the same `Arc` at the end
//!   (`Arc::make_mut`), which is observationally the same because nothing
//!   reads `SeqData` in between,
//! - `Line->Point(ii)` with `ii == NbPoints() + 1` reads past the end of the
//!   sequence (`cxx:2330`, reachable when the inner scan at `cxx:2324-2328`
//!   runs off the end). OCCT's `NCollection_Sequence::Value` does not check the
//!   index; the port stops the scan there instead of reading undefined memory.

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_core::precision::PCONFUSION;
use occt_geom::curve::Curve;
use occt_math::function_set_root::{FunctionWithDerivative, MathFunctionRoot};

use crate::brep_blend::Line;
use crate::chfi3d::search_sing::SearchSing;
use crate::chfi_ds::SurfData;

use super::ChFi3dBuilder;

/// `P.PointOnS1().Distance(P.PointOnS2())` (`cxx:2315`, `:2318`, `:2332`).
fn gap(point: &crate::blend::BlendPoint) -> f64 {
    point.point_on_s1().distance(point.point_on_s2())
}

impl ChFi3dBuilder {
    /// `SplitSurf(SeqData, Line)` (`ChFi3d_FilBuilder.cxx:2274-2436`).
    pub fn split_surf(
        &mut self,
        seq_data: &mut Vec<Arc<SurfData>>,
        line: &Line,
    ) -> Result<(), String> {
        let nbpnt = line.nb_points();
        if nbpnt < 3 {
            return Ok(());
        }
        let Some(first) = seq_data.first() else {
            // `SeqData(1)` on an empty sequence raises in OCCT.
            return Ok(());
        };
        let mut ref_data = SurfData::new();
        ref_data.copy_from(first);

        // `Surf = DStr.Surface(ref->Surf()).Surface(); Surf->Bounds(...);
        // Courbe1 = Surf->UIso(UFirst); Courbe2 = Surf->UIso(ULast);`
        // (`cxx:2288-2292`).
        let (u_first, u_last, surf) = {
            let dstr = self.my_ds.ds();
            let surf_index = ref_data.surf();
            if surf_index < 1 || surf_index > dstr.nb_surfaces() {
                // `TopOpeBRepDS_DataStructure::Surface` raises on a bad index.
                return Err("ChFi3d::SplitSurf: SurfData::Surf out of range".to_string());
            }
            let geometry = dstr
                .surface(surf_index)
                .surface()
                .cloned()
                .ok_or_else(|| "ChFi3d::SplitSurf: null surface in the DS".to_string())?;
            let (u_first, u_last) = geometry.u_range();
            (u_first, u_last, geometry)
        };
        let courbe1: Arc<dyn Curve> = surf
            .u_iso_curve(u_first)
            .ok_or_else(|| "Geom_Surface::UIso: Standard_NoSuchObject".to_string())?;
        let courbe2: Arc<dyn Curve> = surf
            .u_iso_curve(u_last)
            .ok_or_else(|| "Geom_Surface::UIso: Standard_NoSuchObject".to_string())?;
        let mut fonc = SearchSing::new(Arc::clone(&courbe1), Arc::clone(&courbe2));

        // `VFirst = min(Intf1.FirstParameter(), Intf2.FirstParameter());`
        // `VLast = max(Intf1.LastParameter(), Intf2.LastParameter());`
        // (`cxx:2300-2303`).
        let v_first = ref_data
            .interference_on_s1()
            .first_parameter()
            .min(ref_data.interference_on_s2().first_parameter());
        let v_last = ref_data
            .interference_on_s1()
            .last_parameter()
            .max(ref_data.interference_on_s2().last_parameter());

        // (1) Finds vi so that iso v=vi is punctual.
        let mut les_vi: Vec<f64> = Vec::new();

        // (1.1) Finds the first point inside (`cxx:2306-2319`).
        let parameter_at = |index: usize| line.point(index).map(|p| p.parameter());
        let mut ii = 1usize;
        while ii <= nbpnt {
            match parameter_at(ii) {
                Some(p) if p < v_first => ii += 1,
                _ => break,
            }
        }
        if ii == 1 {
            ii += 1;
        }
        let point_ii = line
            .point(ii)
            .ok_or_else(|| "ChFi3d::SplitSurf: BRepBlend_Line::Point out of range".to_string())?;
        let mut b = point_ii.parameter();
        let mut courant = gap(point_ii);
        let point_prev = line.point(ii - 1).ok_or_else(|| {
            "ChFi3d::SplitSurf: BRepBlend_Line::Point out of range".to_string()
        })?;
        let mut a = point_prev.parameter();
        let mut precedant = gap(point_prev);
        ii += 1;

        // (1.2) Find a minimum by "points" (`cxx:2322-2378`).
        while ii <= nbpnt {
            let Some(outer_parameter) = parameter_at(ii) else {
                break;
            };
            if outer_parameter > v_last {
                break;
            }
            loop {
                let Some(inner_parameter) = parameter_at(ii) else {
                    break;
                };
                if inner_parameter < v_last && inner_parameter - b < PCONFUSION {
                    ii += 1;
                } else {
                    break;
                }
            }
            // `const Blend_Point& pnt = Line->Point(ii);` (`cxx:2330`).
            let Some(pnt) = line.point(ii) else {
                break;
            };
            let c = pnt.parameter();
            let suivant = gap(pnt);
            if courant < precedant && courant < suivant {
                // (1.3) Find the exact minimum (`cxx:2336-2365`).
                let resol = MathFunctionRoot::new_with_bounds(
                    &mut fonc,
                    0.5 * (a + c),
                    self.tol2d,
                    a,
                    c,
                    50,
                );
                if resol.is_done() {
                    let racine = resol.root();
                    let mut val = 0.0;
                    fonc.value(racine, &mut val);
                    if val < self.tolapp3d {
                        // the solution (avoiding the risks of confusion)
                        let lower = match les_vi.last() {
                            Some(last) => *last + self.tol2d,
                            None => v_first + self.tol2d,
                        };
                        if racine > lower && racine < v_last - self.tol2d {
                            les_vi.push(racine);
                        }
                    }
                }
                // `else` is only the `#ifdef CHFI3D_DEB` trace (`cxx:2366-2371`).
            }
            // update if non duplication
            a = b;
            precedant = courant;
            b = c;
            courant = suivant;
            ii += 1;
        }

        // (2) Update of the sequence of SurfData (`cxx:2380-2435`).
        if les_vi.is_empty() {
            return Ok(());
        }
        let mut inserted: Vec<Arc<SurfData>> = Vec::with_capacity(les_vi.len());
        for t in &les_vi {
            let t = *t;
            // (2.0) copy and insertion.
            let mut sd = SurfData::new();
            sd.copy_from(&ref_data);
            let dstru = self.my_ds.change_ds();
            let surface = dstru.surface(sd.surf()).clone();
            let new_surf = dstru.add_surface(&surface);
            sd.set_surf(new_surf);
            let mut c1 = dstru
                .curve(sd.interference_on_s1().line_index())
                .clone();
            let new_c1 = dstru.add_curve(&mut c1);
            sd.interference_on_s1_mut().set_line_index(new_c1);
            let mut c2 = dstru
                .curve(sd.interference_on_s2().line_index())
                .clone();
            let new_c2 = dstru.add_curve(&mut c2);
            sd.interference_on_s2_mut().set_line_index(new_c2);

            // (2.1) Modification of common Point (`cxx:2405-2423`).
            sd.vertex_last_on_s1_mut().reset();
            sd.vertex_last_on_s2_mut().reset();
            ref_data.vertex_first_on_s1_mut().reset();
            ref_data.vertex_first_on_s2_mut().reset();
            let p1 = courbe1.value(t);
            let p2 = courbe2.value(t);
            let p3d = GpPnt::new(
                0.5 * (p1.x() + p2.x()),
                0.5 * (p1.y() + p2.y()),
                0.5 * (p1.z() + p2.z()),
            );
            let vertex_tol = p1.distance(&p2) + c1.tolerance().max(c2.tolerance());

            sd.vertex_last_on_s1_mut().set_point(&p3d);
            sd.vertex_last_on_s2_mut().set_point(&p3d);
            ref_data.vertex_first_on_s1_mut().set_point(&p3d);
            ref_data.vertex_first_on_s2_mut().set_point(&p3d);
            sd.vertex_last_on_s1_mut().set_tolerance(vertex_tol);
            sd.vertex_last_on_s2_mut().set_tolerance(vertex_tol);
            ref_data.vertex_first_on_s1_mut().set_tolerance(vertex_tol);
            ref_data.vertex_first_on_s2_mut().set_tolerance(vertex_tol);

            // (2.2) Modification of interferences (`cxx:2425-2433`).
            sd.interference_on_s1_mut().set_last_parameter(t);
            sd.interference_on_s2_mut().set_last_parameter(t);
            ref_data.interference_on_s1_mut().set_first_parameter(t);
            ref_data.interference_on_s2_mut().set_first_parameter(t);

            // Parameters on ElSpine.
            sd.set_last_spine_param(t);
            ref_data.set_first_spine_param(t);

            inserted.push(Arc::new(sd));
        }

        // `SeqData.InsertBefore(ii, SD)` in ascending `ii` puts the new
        // `SurfData` before the original one, in the order of `LesVi`.
        let count = inserted.len();
        seq_data.splice(0..0, inserted);
        if let Some(slot) = seq_data.get_mut(count) {
            *Arc::make_mut(slot) = ref_data;
        }
        Ok(())
    }
}
