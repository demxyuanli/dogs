//! `BRepMAT2d_Explorer` (`BRepMAT2d_Explorer.cxx:63-610`,
//! `.hxx:36-160`): extracts the 2D contours of a face's wires as a figure of
//! `Geom2d_TrimmedCurve`s, and repairs the gaps between consecutive edges.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `Perform` | `cxx:78-100` | [`BRepMat2dExplorer::perform`] |
//! | `Add(Wire, Face, Face)` | `cxx:105-342` | [`BRepMat2dExplorer::add_wire`] |
//! | `Clear` / `NewContour` / `Add(Curve)` | `cxx:372-400` | `clear` / `new_contour` / [`BRepMat2dExplorer::add_curve`] |
//! | `NumberOfContours` / `NumberOfCurves` / `Init` / `More` / `Next` / `Value` | `cxx:402-441` | same names |
//! | `Shape` / `Contour` / `IsModified` / `ModifiedShape` / `GetIsClosed` | `cxx:445-492` | same names |
//! | `MakeEdge` | `cxx:502-519` | [`make_edge`] |
//! | `GetCurveType` | `cxx:526-570` | [`curve_type_rank`] |
//! | `AdjustCurveEnd` | `cxx:576-609` | [`adjust_curve_end`] |
//!
//! Deviations, all forced by the port's data model:
//!
//! * `theCurves` holds `Geom2d_TrimmedCurve` values, not `handle(Geom2d_Curve)`:
//!   every curve this class builds is a trimmed curve, and the consumer
//!   `MAT2d_Circuit` down-casts each entry to `Geom2d_TrimmedCurve` anyway
//!   (`MAT2d_Circuit.cxx:103-105`).
//! * `GetCurveType` answers OCCT's `GeomAbs_CurveType` as a rank rather than
//!   `occt_core::kernel::geomabs::CurveType`: OCCT's enum starts at
//!   `GeomAbs_OtherCurve = 0`, while the port's `CurveType` puts `OtherCurve`
//!   last, and `cxx:194` / `:257` compare the values with `<=`.
//! * `NCollection_IndexedDataMap` (`anOldNewE`, `cxx:116`) becomes a `Vec` in
//!   insertion order, which is what its `FindKey(i)` / `FindFromIndex(i)` walk
//!   (`cxx:312-323`) needs.
//! * `BRep_Tool::Surface(F)` returning a null handle makes OCCT's
//!   `BRepBuilderAPI_MakeFace` dereference it; the port panics there instead.

use std::sync::Arc;

use occt_core::convert::ParameterisationType;
use occt_core::gp::GpPnt2d;
use occt_core::precision::CONFUSION;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::gce2d;
use occt_geom2d::geom2d_convert::curve_to_bspline_curve_bspl;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_lib_build_curves3d::build_curves3d;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{is_same, shape_map_key, wires_of_face, ShapeMapKey};
use crate::topopebrep::build::tools2d::{first_vertex, last_vertex};

/// The tolerance `Add(Wire, ...)` uses to detect a gap between two consecutive
/// edges (`cxx:139`: `aTolConf = 1.e-8`, deliberately *not*
/// `Precision::Confusion()`).
const GAP_TOL: f64 = 1.0e-8;

/// `GeomAbs_CurveType` ranks of [`curve_type_rank`], with OCCT's numeric order
/// (`GeomAbs_CurveType.hxx`): `OtherCurve` is `0`, so an unrecognised curve
/// ranks below every analytic one.
const RANK_OTHER: u8 = 0;
const RANK_LINE: u8 = 1;

/// `class BRepMAT2d_Explorer` (`BRepMAT2d_Explorer.hxx:36-160`).
#[derive(Clone, Default)]
pub struct BRepMat2dExplorer {
    /// `theCurves` (`.hxx`): one sequence of curves per contour.
    the_curves: Vec<Vec<Geom2dTrimmedCurve>>,
    /// `currentContour` (`.hxx`), 1-based like OCCT's `int`.
    current_contour: usize,
    /// `current` (`.hxx`), 1-based cursor inside `theCurves(currentContour)`.
    current: usize,
    /// `myIsClosed` (`.hxx`).
    my_is_closed: Vec<bool>,
    /// `myModifShapes` (`.hxx`), an `NCollection_DataMap<Shape, Shape>`:
    /// `(key, original, replacement)` in insertion order.
    my_modif_shapes: Vec<(ShapeMapKey, TopoShape, TopoShape)>,
    /// `myShape` (`.hxx`).
    my_shape: Option<TopoShape>,
}

impl BRepMat2dExplorer {
    /// `BRepMAT2d_Explorer()` (`cxx:65-68`).
    pub fn new() -> Self {
        let mut explorer = Self::default();
        explorer.clear();
        explorer
    }

    /// `BRepMAT2d_Explorer(aFace)` (`cxx:72-75`).
    pub fn with_face(a_face: &Face) -> Self {
        let mut explorer = Self::new();
        explorer.perform(a_face);
        explorer
    }

    /// `Perform(aFace)` (`cxx:78-100`).
    pub fn perform(&mut self, a_face: &Face) {
        self.clear();
        self.my_shape = Some(a_face.0.clone());

        // `cxx:82-84`: `TopoDS::Face(aFace); F.Orientation(TopAbs_FORWARD)`.
        let mut f = a_face.0.clone();
        f.set_orientation(Orientation::Forward);
        let f = Face(f);

        // `cxx:87-88`: `BRep_Tool::Surface(F)` then
        // `BRepBuilderAPI_MakeFace(aSurf, Precision::Confusion())`. OCCT
        // dereferences the surface handle without a null check.
        let surf = GeometryRegistry::global()
            .face_surface(&f.0)
            .unwrap_or_else(|| panic!("BRepMAT2d_Explorer::Perform: BRep_Tool::Surface is null"));
        let mut new_face = TopoBuilder::new().make_face_with_tolerance(surf, &[], CONFUSION);

        // `cxx:89-93`: `TopExp_Explorer(F, TopAbs_WIRE)`.
        let wires = wires_of_face(&f);
        for wire in wires {
            self.add_wire(&wire, &f, &mut new_face);
        }

        // `cxx:95`: `BRepLib::BuildCurves3d(aNewF)`.
        build_curves3d(&new_face.0);

        // `cxx:97`: `myModifShapes.Add(aFace, aNewF)`.
        self.modif_shapes_set(&a_face.0, &new_face.0);
    }

    /// `Add(Spine, aFace, aNewFace)` (`cxx:105-342`): walks the wire's edges,
    /// turns each into a 2D trimmed curve of `aFace`, and closes the gaps
    /// between them by either moving the current curve's start pole or the
    /// previous curve's end pole. When anything was repaired, a new wire of new
    /// edges is added to `a_new_face` in place of `spine`.
    pub fn add_wire(&mut self, spine: &Wire, a_face: &Face, a_new_face: &mut Face) {
        // `cxx:112-113`: `BRepTools_WireExplorer anExp(Spine, aFace)`. The spine
        // is *not* oriented FORWARD first (see the `skv` note at `cxx:108-114`).
        let mut an_exp = WireExplorer::of_wire_on_face(spine, a_face);
        // `cxx:116`: `anOldNewE`, an `IndexedDataMap<Shape, Shape>`.
        let mut an_old_new_e: Vec<(ShapeMapKey, TopoShape, TopoShape)> = Vec::new();
        if !an_exp.more() {
            return;
        }

        // `cxx:124-126`.
        self.new_contour();
        let contour = self.current_contour;
        self.my_is_closed[contour - 1] = spine.0.is_closed();

        // `cxx:128-135`.
        let a_first_edge = an_exp
            .current_edge()
            .expect("BRepMAT2d_Explorer::Add: empty wire explorer");
        let mut a_prev_edge = a_first_edge.clone();
        let a_tol_conf = GAP_TOL;
        let mut is_modif = false;

        // `cxx:143-145`: the first edge of the wire.
        indexed_add(&mut an_old_new_e, &a_first_edge.0, &a_first_edge.0);
        let mut ct2d = trimmed_on_face(&a_first_edge, a_face);
        if a_first_edge.0.orientation() == Orientation::Reversed {
            // `cxx:147-150`.
            ct2d.reverse();
        }

        // `cxx:152-153`.
        let a_p_first = ct2d.d0(ct2d.first_parameter());
        let mut a_p_last = ct2d.d0(ct2d.last_parameter());

        // `cxx:155-157`.
        self.add_curve(ct2d.clone());
        let mut a_first_curve = ct2d;
        an_exp.next();

        // `cxx:160-242`: the next edges.
        while an_exp.more() {
            let an_edge = an_exp
                .current_edge()
                .expect("BRepMAT2d_Explorer::Add: wire explorer lost its edge");

            indexed_add(&mut an_old_new_e, &an_edge.0, &an_edge.0);
            let mut ct2d = trimmed_on_face(&an_edge, a_face);
            if an_edge.0.orientation() == Orientation::Reversed {
                // `cxx:168-171`.
                ct2d.reverse();
            }

            // `cxx:173-179`.
            let a_p_cur_first = ct2d.d0(ct2d.first_parameter());
            let a_d = a_p_last.distance(&a_p_cur_first);
            if a_d > a_tol_conf {
                // `cxx:186-190`: the last curve of the contour and its type.
                let a_nb_c = self.the_curves[contour - 1].len();
                let mut c_prev = self.the_curves[contour - 1][a_nb_c - 1].clone();
                let t_c_prev = curve_type_rank(&c_prev);
                let t_c_curr = curve_type_rank(&ct2d);
                is_modif = true;
                if t_c_curr <= t_c_prev {
                    // `cxx:195-213`: move the *current* curve's first pole onto
                    // the previous curve's last point, and rebuild the current
                    // edge.
                    ct2d = adjust_curve_end(&ct2d, &a_p_last, true);
                    let a_vf = first_vertex(&an_edge.0);
                    let a_vl = last_vertex(&an_edge.0);
                    let new_edge = if an_edge.0.orientation() == Orientation::Forward {
                        make_edge(&ct2d, a_new_face, a_vf.as_ref(), a_vl.as_ref())
                    } else {
                        let mut reversed = ct2d.clone();
                        reversed.reverse();
                        make_edge(&reversed, a_new_face, a_vf.as_ref(), a_vl.as_ref())
                    };
                    let mut new_edge = new_edge;
                    new_edge
                        .0
                        .set_orientation(an_edge.0.orientation());
                    indexed_set(&mut an_old_new_e, &an_edge.0, &new_edge.0);
                } else {
                    // `cxx:214-236`: move the *previous* curve's last pole onto
                    // the current curve's first point, and rebuild the previous
                    // edge.
                    let a_p = ct2d.d0(ct2d.first_parameter());
                    c_prev = adjust_curve_end(&c_prev, &a_p, false);
                    self.the_curves[contour - 1][a_nb_c - 1] = c_prev.clone();
                    let a_vf = first_vertex(&a_prev_edge.0);
                    let a_vl = last_vertex(&a_prev_edge.0);
                    let new_edge = if a_prev_edge.0.orientation() == Orientation::Forward {
                        make_edge(&c_prev, a_new_face, a_vf.as_ref(), a_vl.as_ref())
                    } else {
                        let mut reversed = c_prev.clone();
                        reversed.reverse();
                        make_edge(&reversed, a_new_face, a_vf.as_ref(), a_vl.as_ref())
                    };
                    let mut new_edge = new_edge;
                    new_edge
                        .0
                        .set_orientation(a_prev_edge.0.orientation());
                    indexed_set(&mut an_old_new_e, &a_prev_edge.0, &new_edge.0);
                }
            }

            // `cxx:239-241`.
            a_p_last = ct2d.d0(ct2d.last_parameter());
            self.add_curve(ct2d);
            a_prev_edge = an_edge;
            // `cxx:160`: the `for` loop's `anExp.Next()`.
            an_exp.next();
        }

        // `cxx:246-304`: distance between the first and the last point of a
        // closed wire.
        if self.my_is_closed[contour - 1] && a_p_last.distance(&a_p_first) > a_tol_conf {
            let a_nb_c = self.the_curves[contour - 1].len();
            let mut c_prev = self.the_curves[contour - 1][a_nb_c - 1].clone();
            let t_c_prev = curve_type_rank(&c_prev);
            let t_c_curr = curve_type_rank(&a_first_curve);
            is_modif = true;
            if t_c_curr <= t_c_prev {
                // `cxx:260-278`.
                a_first_curve = adjust_curve_end(&a_first_curve, &a_p_last, true);
                self.the_curves[contour - 1][0] = a_first_curve.clone();
                let a_vf = first_vertex(&a_first_edge.0);
                let a_vl = last_vertex(&a_first_edge.0);
                let new_edge = if a_first_edge.0.orientation() == Orientation::Forward {
                    make_edge(&a_first_curve, a_new_face, a_vf.as_ref(), a_vl.as_ref())
                } else {
                    let mut reversed = a_first_curve.clone();
                    reversed.reverse();
                    make_edge(&reversed, a_new_face, a_vf.as_ref(), a_vl.as_ref())
                };
                let mut new_edge = new_edge;
                new_edge
                    .0
                    .set_orientation(a_first_edge.0.orientation());
                indexed_set(&mut an_old_new_e, &a_first_edge.0, &new_edge.0);
            } else {
                // `cxx:281-302`.
                let a_p = a_first_curve.d0(a_first_curve.first_parameter());
                c_prev = adjust_curve_end(&c_prev, &a_p, false);
                self.the_curves[contour - 1][a_nb_c - 1] = c_prev.clone();
                let a_vf = first_vertex(&a_prev_edge.0);
                let a_vl = last_vertex(&a_prev_edge.0);
                let new_edge = if a_prev_edge.0.orientation() == Orientation::Forward {
                    make_edge(&c_prev, a_new_face, a_vf.as_ref(), a_vl.as_ref())
                } else {
                    let mut reversed = c_prev.clone();
                    reversed.reverse();
                    make_edge(&reversed, a_new_face, a_vf.as_ref(), a_vl.as_ref())
                };
                let mut new_edge = new_edge;
                new_edge
                    .0
                    .set_orientation(a_prev_edge.0.orientation());
                indexed_set(&mut an_old_new_e, &a_prev_edge.0, &new_edge.0);
            }
        }

        // `cxx:305-339`. The new wire is the spine itself when nothing was
        // repaired.
        let a_new_wire = if is_modif {
            let mut new_wire = Wire::new();
            for (_, a_key, a_new_edge) in &an_old_new_e {
                // `cxx:320` / `cxx:322`.
                TopoBuilder::new().add_edge(&mut new_wire, &Edge(a_new_edge.clone()));
                self.modif_shapes_set(a_key, a_new_edge);
            }
            if self.my_is_closed[contour - 1] {
                new_wire.0.set_closed(true);
            }
            self.modif_shapes_set(&spine.0, &new_wire.0);
            new_wire
        } else {
            spine.clone()
        };

        // `cxx:341`: `aBuilder.Add(aNewFace, aNewWire)`.
        TopoBuilder::new().add_wire(a_new_face, &a_new_wire);
    }

    /// `Clear()` (`cxx:372-379`).
    pub fn clear(&mut self) {
        self.the_curves.clear();
        self.current_contour = 0;
        self.my_is_closed.clear();
        self.my_modif_shapes.clear();
    }

    /// `NewContour()` (`cxx:384-392`).
    fn new_contour(&mut self) {
        self.the_curves.push(Vec::new());
        self.my_is_closed.push(false);
        self.current_contour += 1;
    }

    /// `Add(aCurve)` (`cxx:395-398`).
    fn add_curve(&mut self, a_curve: Geom2dTrimmedCurve) {
        self.the_curves[self.current_contour - 1].push(a_curve);
    }

    /// `NumberOfContours()` (`cxx:403-406`).
    pub fn number_of_contours(&self) -> usize {
        self.the_curves.len()
    }

    /// `NumberOfCurves(IndexContour)` (`cxx:410-413`).
    pub fn number_of_curves(&self, index_contour: usize) -> usize {
        self.the_curves[index_contour - 1].len()
    }

    /// `Init(IndexContour)` (`cxx:417-420`).
    pub fn init(&mut self, index_contour: usize) {
        self.current_contour = index_contour;
        self.current = 1;
    }

    /// `More()` (`cxx:425-428`).
    pub fn more(&self) -> bool {
        self.current <= self.number_of_curves(self.current_contour)
    }

    /// `Next()` (`cxx:432-435`).
    pub fn next(&mut self) {
        self.current += 1;
    }

    /// `Value()` (`cxx:440-443`).
    pub fn value(&self) -> &Geom2dTrimmedCurve {
        &self.the_curves[self.current_contour - 1][self.current - 1]
    }

    /// `Shape()` (`cxx:447-450`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `Contour(IC)` (`cxx:454-458`).
    pub fn contour(&self, ic: usize) -> &[Geom2dTrimmedCurve] {
        &self.the_curves[ic - 1]
    }

    /// `IsModified(aShape)` (`cxx:462-473`).
    pub fn is_modified(&self, a_shape: &TopoShape) -> bool {
        match self.modif_shapes_get(a_shape) {
            Some(new_shape) => !is_same(&new_shape, a_shape),
            None => false,
        }
    }

    /// `ModifiedShape(aShape)` (`cxx:477-487`).
    pub fn modified_shape(&self, a_shape: &TopoShape) -> TopoShape {
        self.modif_shapes_get(a_shape)
            .unwrap_or_else(|| a_shape.clone())
    }

    /// `GetIsClosed()` (`cxx:491-494`).
    pub fn get_is_closed(&self) -> &[bool] {
        &self.my_is_closed
    }

    // ---------------------------------------------------------------
    // `myModifShapes`, an `NCollection_DataMap<Shape, Shape>`
    // ---------------------------------------------------------------

    /// `myModifShapes.Add(key, value)`: a bound key keeps its position and takes
    /// the new value.
    fn modif_shapes_set(&mut self, key: &TopoShape, value: &TopoShape) {
        let k = shape_map_key(key);
        self.modif_shapes_set_with_key(k, key, value);
    }

    fn modif_shapes_set_with_key(&mut self, k: ShapeMapKey, key: &TopoShape, value: &TopoShape) {
        match self.my_modif_shapes.iter_mut().find(|(kk, _, _)| *kk == k) {
            Some(entry) => entry.2 = value.clone(),
            None => self
                .my_modif_shapes
                .push((k, key.clone(), value.clone())),
        }
    }

    /// `myModifShapes.FindFromKey(aShape)`.
    fn modif_shapes_get(&self, a_shape: &TopoShape) -> Option<TopoShape> {
        let k = shape_map_key(a_shape);
        self.my_modif_shapes
            .iter()
            .find(|(kk, _, _)| *kk == k)
            .map(|(_, _, v)| v.clone())
    }
}

/// The 2D curve of `edge` on `face` as a `Geom2d_TrimmedCurve`
/// (`cxx:144-145`, `:165-166`): `BRep_Tool::CurveOnSurface(edge, face, UFirst,
/// ULast)` wrapped in a trim.
fn trimmed_on_face(edge: &Edge, face: &Face) -> Geom2dTrimmedCurve {
    let (curve, first, last) = curve_on_surface_range(edge, face).unwrap_or_else(|| {
        panic!("BRepMAT2d_Explorer::Add: BRep_Tool::CurveOnSurface returned a null curve")
    });
    Geom2dTrimmedCurve::new(curve, first, last)
}

/// `MakeEdge(theCurve, theFace, theVFirst, theVLast)` (`cxx:502-519`):
/// `BRep_Builder::MakeEdge` + `UpdateEdge(E, C2d, F, Confusion)` + the two
/// vertices + `Range(E, First, Last)`.
fn make_edge(
    the_curve: &Geom2dTrimmedCurve,
    the_face: &Face,
    the_v_first: Option<&TopoShape>,
    the_v_last: Option<&TopoShape>,
) -> Edge {
    let builder = TopoBuilder::new();
    let reg = GeometryRegistry::global();
    let mut edge = Edge::new();
    let face_key = GeometryRegistry::shape_key(&the_face.0);
    let a_f_par = the_curve.first_parameter();
    let a_l_par = the_curve.last_parameter();

    // `cxx:511-512`: `UpdateEdge(E, C2d, F, Precision::Confusion())`.
    reg.set_edge_pcurve(&edge.0, face_key, Arc::new(the_curve.clone()));
    reg.raise_edge_tolerance(&edge.0, CONFUSION);

    // `cxx:513-514`: the first vertex FORWARD, the last one REVERSED.
    if let Some(v) = the_v_first {
        builder.add(&mut edge.0, &v.oriented(Orientation::Forward));
    }
    if let Some(v) = the_v_last {
        builder.add(&mut edge.0, &v.oriented(Orientation::Reversed));
    }

    // `cxx:515`: `Range(E, First, Last)`.
    reg.set_ranges_all(&edge.0, a_f_par, a_l_par);
    edge
}

/// `GetCurveType(theC2d)` (`cxx:526-570`) as the numeric value of OCCT's
/// `GeomAbs_CurveType`.
///
/// OCCT unwraps exactly one `Geom2d_TrimmedCurve` level (`cxx:529-531`) and then
/// compares `DynamicType()`; the port replaces the `STANDARD_TYPE` chain with
/// the `Curve2d` exact-type predicates.
fn curve_type_rank(c: &Geom2dTrimmedCurve) -> u8 {
    let t: &dyn Curve2d = c.trimmed_basis().unwrap_or(c);
    if t.is_geom2d_circle() {
        2 // GeomAbs_Circle
    } else if t.is_geom2d_line() {
        RANK_LINE
    } else if t.gp_elips2d().is_some() {
        3 // GeomAbs_Ellipse
    } else if t.gp_parab2d().is_some() {
        5 // GeomAbs_Parabola
    } else if t.gp_hypr2d().is_some() {
        4 // GeomAbs_Hyperbola
    } else if t.is_bezier2d() {
        6 // GeomAbs_BezierCurve
    } else if t.is_bspline2d() {
        7 // GeomAbs_BSplineCurve
    } else if t.offset_basis().is_some() {
        8 // GeomAbs_OffsetCurve
    } else {
        RANK_OTHER
    }
}

/// `AdjustCurveEnd(theC2d, theP, isFirst)` (`cxx:576-609`).
///
/// A line is rebuilt as a segment between the moved point and the untouched end;
/// every other curve is converted to a `Convert_QuasiAngular` B-spline whose
/// first (or last) pole is moved. `CurveToBSplineCurve` answering null is a
/// `Standard_DomainError` in OCCT (`Geom2dConvert.cxx:442-444`).
fn adjust_curve_end(
    the_c2d: &Geom2dTrimmedCurve,
    the_p: &GpPnt2d,
    is_first: bool,
) -> Geom2dTrimmedCurve {
    if curve_type_rank(the_c2d) == RANK_LINE {
        // `cxx:580-592`.
        if is_first {
            let a_p = the_c2d.d0(the_c2d.last_parameter());
            gce2d::make_segment(the_p, &a_p)
        } else {
            let a_p = the_c2d.d0(the_c2d.first_parameter());
            gce2d::make_segment(&a_p, the_p)
        }
    } else {
        // `cxx:594-608`.
        let mut b_curve = curve_to_bspline_curve_bspl(the_c2d, ParameterisationType::QuasiAngular)
            .unwrap_or_else(|| {
                panic!("BRepMAT2d_Explorer::AdjustCurveEnd: CurveToBSplineCurve failed")
            });
        if is_first {
            b_curve.set_pole(1, the_p);
        } else {
            let nb = b_curve.nb_poles();
            b_curve.set_pole(nb, the_p);
        }
        let (first, last) = (b_curve.first_parameter(), b_curve.last_parameter());
        Geom2dTrimmedCurve::new(Arc::new(b_curve), first, last)
    }
}

/// `IndexedDataMap::Add(key, value)`.
fn indexed_add(
    map: &mut Vec<(ShapeMapKey, TopoShape, TopoShape)>,
    key: &TopoShape,
    value: &TopoShape,
) {
    let k = shape_map_key(key);
    match map.iter_mut().find(|(kk, _, _)| *kk == k) {
        Some(entry) => entry.2 = value.clone(),
        None => map.push((k, key.clone(), value.clone())),
    }
}

/// `IndexedDataMap::ChangeFromKey(key) = value`.
fn indexed_set(
    map: &mut [(ShapeMapKey, TopoShape, TopoShape)],
    key: &TopoShape,
    value: &TopoShape,
) {
    let k = shape_map_key(key);
    let entry = map
        .iter_mut()
        .find(|(kk, _, _)| *kk == k)
        .unwrap_or_else(|| {
            panic!("BRepMAT2d_Explorer::Add: ChangeFromKey on an unbound edge")
        });
    entry.2 = value.clone();
}
