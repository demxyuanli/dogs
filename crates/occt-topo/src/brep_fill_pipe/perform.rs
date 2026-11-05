//! `BRepFill_Pipe::Perform` and `DefineRealSegmax`
//! (`BRepFill_Pipe.cxx:193-311`, `:999-1049`).

use std::sync::Arc;

use occt_core::gp::{GpMat, GpTrsf, GpVec};
use occt_core::precision::PCONFUSION;
use occt_core::toploc::TopLocLocation;

use occt_geom::curve::Curve;
use occt_geom::geom_fill::{
    CorrectedFrenet, CurveAndTrihedron, DiscreteTrihedron, Frenet, GeomFillTrihedron,
    LocationLaw as GeomFillLocationLaw, TrihedronLaw,
};

use crate::abs::ShapeType;
use crate::brep_fill_edge3d_law::Edge3DLaw;
use crate::brep_fill_pipe::Pipe;
use crate::brep_fill_section_placement::SectionPlacement;
use crate::brep_tool::BRepTool;
use crate::iterator::ShapeIterator;
use crate::shape::{TopoShape, Wire};
use crate::shape_ops::transformed_copy;
use crate::shape_upgrade::RemoveLocations;

impl Pipe {
    /// `Perform(Spine, Profile, GeneratePartCase)` (`cxx:193-311`).
    ///
    /// `GeneratePartCase` is unused by the reference (`cxx:195` comments the
    /// parameter out).
    pub fn perform(
        &mut self,
        spine: &Wire,
        profile: &TopoShape,
        generate_part_case: bool,
    ) -> Result<(), String> {
        let _ = generate_part_case;

        // `cxx:198-205`.
        self.set_sections(None);
        self.set_faces(None);
        self.set_edges(None);
        self.set_spine(Some(spine.clone()));
        self.set_profile(Some(profile.clone()));

        self.define_real_segmax();

        // `cxx:210-226`: the elementary trihedron law, then the location law
        // over the spine.
        let trihedron: Arc<dyn TrihedronLaw> = match self.mode() {
            GeomFillTrihedron::Frenet => Arc::new(Frenet::new()),
            GeomFillTrihedron::CorrectedFrenet => Arc::new(CorrectedFrenet::new()),
            GeomFillTrihedron::DiscreteTrihedron => Arc::new(DiscreteTrihedron::new()),
            // `cxx:223-224`: the `default` arm leaves `TLaw` null and the
            // construction of `GeomFill_CurveAndTrihedron` then fails. The
            // constructor only stores the three modes above, so the arm is
            // unreachable; the port answers an error instead of a null law.
            _ => return Err("BRepFill_Pipe::Perform: null trihedron law".to_string()),
        };
        let geom_loc: Arc<dyn GeomFillLocationLaw> = Arc::new(CurveAndTrihedron::new(trihedron));
        let mut law = Edge3DLaw::new(spine, &geom_loc);
        if law.base().nb_law() == 0 {
            return Ok(()); // Degenerated case (`cxx:228-230`)
        }
        law.base_mut().transform_in_g0_law()?;

        // `cxx:232-234`: the section placement, whose transform is `myTrsf`.
        let place = SectionPlacement::new(law.base(), profile, false, false)?;
        self.set_trsf(place.transformation().clone());
        self.set_loc(Some(Arc::new(law)));

        // `cxx:236-240`: the profile moved into the local frame.
        let loc1 = profile.location().clone();
        let mut the_prof = profile.clone();
        the_prof.set_location(&TopLocLocation::composed(self.trsf(), loc1));

        // `cxx:243-273`: "Construct First && Last Shape".
        let location = self.location()?;
        let first_law = location
            .base()
            .law(1)
            .ok_or_else(|| "BRepFill_Pipe::Perform: no location law 1".to_string())?;
        let (first, _last) = first_law.get_domain();
        let (m, v) = first_law.d0(first)?;
        let mut fila = frame_trsf(&m, &v)?;
        fila.multiply(self.trsf());
        let loc_first = TopLocLocation::composed(&fila, TopLocLocation::identity());

        let mut my_first = profile.clone();
        if !loc_first.is_identity() {
            // `cxx:270`: `myFirst = BRepBuilderAPI_Transform(myProfile, fila, true)`.
            my_first = transformed_copy(profile, &fila)?;
        }
        // `cxx:274-278`: `ShapeUpgrade_RemoveLocations(COMPOUND)`.
        let mut rem_loc = RemoveLocations::new();
        rem_loc.set_remove_level(ShapeType::Compound);
        rem_loc.remove(&my_first);
        let my_first = rem_loc.result().cloned();

        // `cxx:280-292`: the frame at the end of the last law.
        let nb_law = location.base().nb_law();
        let last_law = location
            .base()
            .law(nb_law)
            .ok_or_else(|| "BRepFill_Pipe::Perform: no last location law".to_string())?;
        let (_first, last) = last_law.get_domain();
        let (m, v) = last_law.d0(last)?;
        let mut fila_last = frame_trsf(&m, &v)?;
        fila_last.multiply(self.trsf());
        let loc_last = TopLocLocation::composed(&fila_last, TopLocLocation::identity());

        // `cxx:298`: `if (!myLoc->IsClosed() || LocFirst != LocLast)`.
        // `TopLoc_Location::operator!=` (`TopLoc_Location.hxx:118-126`) compares
        // the two `TopLoc_Datum3D` handles, and the two datums above are
        // distinct objects, so the disjunction is constantly true and the
        // `else myLast = myFirst` arm (`cxx:305-307`) is unreachable. The port
        // therefore keeps the reference's first arm only; there is no datum
        // identity to test here (`TopLocLocation::is_equal` compares values,
        // which is `IsEqual`, not the operator OCCT applies).
        let mut my_last = profile.clone();
        if !loc_last.is_identity() {
            // `cxx:299-302`.
            my_last = transformed_copy(profile, &fila_last)?;
        }
        // `cxx:309-311`: the same `RemoveLocations` instance, as in the
        // reference (its `myMapNewShapes` is reused).
        rem_loc.remove(&my_last);
        let my_last = rem_loc.result().cloned();

        self.set_first(my_first.clone());
        self.set_last(my_last.clone());

        // `cxx:313`.
        let original = profile.clone();
        let shape = self.make_shape(&the_prof, &original, my_first.as_ref(), my_last.as_ref())?;
        self.set_shape(shape);
        Ok(())
    }

    /// `DefineRealSegmax()` (`cxx:999-1049`): raise `mySegmax` to the number of
    /// spine B-spline knot spans the trimmed edges actually use.
    fn define_real_segmax(&mut self) {
        let Some(spine) = self.spine().cloned() else {
            return;
        };
        let mut real_segmax = 0i32;
        for edge in ShapeIterator::of_shape(&spine.0) {
            let Some(the_edge) = crate::shape::Edge::wrap(edge) else {
                continue;
            };
            let (first, last) = BRepTool::edge_parameters(&the_edge);
            let Some(curve) = BRepTool::edge_curve(&the_edge) else {
                continue; // `cxx:1009-1012`: a null curve is skipped
            };

            // `cxx:1011-1023`: unwrap `Geom_TrimmedCurve` / `Geom_OffsetCurve`
            // down to the basis curve, checking the offset arm on the curve the
            // trim arm just produced (the two `if`s share one loop body).
            let mut curve = curve;
            loop {
                if curve.is_geom_trimmed() {
                    match curve.trim_basis() {
                        Some(basis) => curve = basis,
                        None => break,
                    }
                    if let Some((basis, _offset)) = curve.offset_curve() {
                        curve = basis;
                    }
                    continue;
                }
                match curve.offset_curve() {
                    Some((basis, _offset)) => curve = basis,
                    None => break,
                }
            }

            // `cxx:1025-1046`.
            let Some(bs) = curve.as_bspline_curve() else {
                continue;
            };
            let nb_knots = bs.nb_knots() as i32;
            let mut real_nb_knots = nb_knots;
            if first > bs.first_parameter() {
                let (i1, _i2) = locate_u(bs, first, PCONFUSION);
                real_nb_knots -= i1 - 1;
            }
            if last < bs.last_parameter() {
                let (i1, i2) = locate_u(bs, last, PCONFUSION);
                let _ = i1;
                real_nb_knots -= nb_knots - i2;
            }
            real_segmax += real_nb_knots - 1;
        }

        // `cxx:1048-1049`.
        if self.segmax() < real_segmax {
            self.set_segmax(real_segmax);
        }
    }
}

/// `gp_Trsf::SetValues(M(1,1)..M(3,3), V.XYZ())` (`cxx:250-262`, `:281-292`).
fn frame_trsf(m: &GpMat, v: &GpVec) -> Result<GpTrsf, String> {
    let mut trsf = GpTrsf::identity();
    trsf.set_values(
        m.value(1, 1),
        m.value(1, 2),
        m.value(1, 3),
        v.x(),
        m.value(2, 1),
        m.value(2, 2),
        m.value(2, 3),
        v.y(),
        m.value(3, 1),
        m.value(3, 2),
        m.value(3, 3),
        v.z(),
    )
    .map_err(|e| format!("gp_Trsf::SetValues: {e}"))?;
    Ok(trsf)
}

/// `Geom_BSplineCurve::LocateU` without knot repetition
/// (`Geom_BSplineCurve_1.cxx:702-760`, `WithKnotRepetition = false`): the
/// `(I1, I2)` bracket of `u` in 1-based distinct-knot indices.
///
/// Same body as the 2D port
/// (`occt-geom2d/src/bspline_curve/curve.rs::locate_u`), which is why
/// `BSplCLib` is not involved: `LocateU` is a pure knot-window search.
fn locate_u(bs: &occt_geom::GeomBSplineCurve, u: f64, parametric_tolerance: f64) -> (i32, i32) {
    let uknots = bs.distinct_knots();
    if uknots.is_empty() {
        return (0, 0);
    }
    let mut new_u = u;
    if bs.periodic() {
        // `PeriodicNormalization` (`cxx:711`).
        new_u = occt_core::bspl::locate::in_period(
            new_u,
            bs.first_parameter(),
            bs.last_parameter(),
        );
    }
    let u_first = uknots[0];
    let u_last = uknots[uknots.len() - 1];
    let tol = parametric_tolerance.abs();
    if (new_u - u_first).abs() <= tol {
        return (1, 1);
    }
    if (new_u - u_last).abs() <= tol {
        let n = uknots.len() as i32;
        return (n, n);
    }
    if new_u < u_first {
        return (0, 1);
    }
    if new_u > u_last {
        let n = uknots.len() as i32;
        return (n, n + 1);
    }
    let mut i1 = occt_core::bspl::locate::hunt_occt(uknots, new_u);
    i1 = i1.max(1).min(uknots.len() as i32);
    while (i1 as usize) < uknots.len() && (uknots[i1 as usize] - new_u).abs() <= tol {
        i1 += 1;
    }
    if (uknots[(i1 - 1) as usize] - new_u).abs() <= tol {
        (i1, i1)
    } else {
        (i1, i1 + 1)
    }
}
