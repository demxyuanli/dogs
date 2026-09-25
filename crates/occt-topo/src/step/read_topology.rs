use super::prelude::*;
use super::*;

/// Direct child edges of a wire **without** composing the wire orientation
/// (`TopoDS_Iterator` with `cumOri=false`). Used for `TranslateEdgeLoop`
/// SelectForwardSeam `EdgeO` (`cxx:700`).
fn edges_stored_on_wire(wire: &Wire) -> Vec<Edge> {
    wire.0
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .iter()
        .filter(|s| s.shape_type() == ShapeType::Edge)
        .cloned()
        .map(Edge)
        .collect()
}

/// The pole-type-independent part of `StepToGeom::MakeBSplineCurveCommon`
/// (`StepToGeom.cxx:776-894`), shared by the 3D and 2D arms (OCCT implements
/// both with one template, `:938-944` vs `:952-963`).
struct BsplineDescriptor {
    uknots: Vec<f64>,
    umults: Vec<i32>,
    /// Pole window `[start, end)` left after the multiplicity trim.
    pole_range: std::ops::Range<usize>,
    should_be_periodic: bool,
}

fn bspline_descriptor(
    degree: usize,
    n_poles: usize,
    mults: &[usize],
    knots: &[f64],
) -> Result<BsplineDescriptor, String> {
    // 1. unique knots + summed multiplicities (`cxx:784-821`).
    let mut uknots: Vec<f64> = Vec::new();
    let mut umults: Vec<i32> = Vec::new();
    for (i, &k) in knots.iter().enumerate() {
        let m = mults.get(i).copied().unwrap_or(1) as i32;
        if let Some(&last) = uknots.last() {
            if k - last <= occt_core::precision::epsilon(last.abs()) {
                *umults.last_mut().unwrap() += m;
                continue;
            }
        }
        uknots.push(k);
        umults.push(m);
    }
    if uknots.len() <= 1 {
        // `if (NbUniqueKnots <= 1) return nullptr;` (`cxx:795-798`).
        return Err("B_SPLINE_CURVE: at most one unique knot".into());
    }

    // 2. clamp multiplicities above `degree + 1` (`cxx:823-845`).
    let deg1 = degree as i32 + 1;
    let n_unique = umults.len();
    let mut first_diff = 0usize;
    let mut last_diff = 0usize;
    for i in 0..n_unique {
        if umults[i] > deg1 {
            if i == 0 {
                first_diff = (umults[i] - deg1) as usize;
            }
            if i == n_unique - 1 {
                last_diff = (umults[i] - deg1) as usize;
            }
            umults[i] = deg1;
        }
    }

    // 3. drop the extra poles (`cxx:847-871`).
    let diff = first_diff + last_diff;
    let n_unique_poles = n_poles
        .checked_sub(diff)
        .filter(|n| *n > 0)
        .ok_or_else(|| "B_SPLINE_CURVE: no poles left after multiplicity trim".to_string())?;

    // 4. periodic-looking descriptor (`cxx:873-894`).
    let summary: i32 = umults.iter().sum();
    let should_be_periodic = summary != (n_unique_poles as i32 + degree as i32 + 1)
        && umults[0] == umults[n_unique - 1]
        && (summary - umults[0]) == n_unique_poles as i32;

    Ok(BsplineDescriptor {
        uknots,
        umults,
        pole_range: first_diff..first_diff + n_unique_poles,
        should_be_periodic,
    })
}

/// `StepToGeom::MakeBSplineCurveCommon` (`StepToGeom.cxx:750-927`), the part
/// shared by the `B_SPLINE_CURVE_WITH_KNOTS` and Bezier / uniform /
/// quasi-uniform arms once the knot and multiplicity lists are known:
///
/// 1. duplicate knots are merged within `Epsilon(|lastKnot|)` and their
///    multiplicities summed (`cxx:784-821`);
/// 2. multiplicities above `degree + 1` are clamped to `degree + 1` and the
///    corresponding leading / trailing poles (and weights) are dropped
///    (`cxx:823-871`, `cxx:902-906`);
/// 3. a descriptor that "looks periodic" (`cxx:873-894`: the total multiplicity
///    is *not* `NbPoles + degree + 1`, the end multiplicities are equal and
///    `Σmults - mults(1) == NbPoles`) is built in the periodic representation,
///    i.e. `Geom_BSplineCurve(Poles, Knots, Mults, Degree, true)` whose flat
///    sequence is the periodic `BSplCLib::KnotSequence`;
/// 4. a closed curve (`ClosedCurve() && Degree() > 1 && IsClosed()`,
///    `cxx:920-926`) is forced periodic with `SetPeriodic()`.
///
/// `IsClosed()` is `StartPoint().SquareDistance(EndPoint()) <=
/// Precision::Computational()` (`Geom_BSplineCurve.cxx:146-149`).
fn make_bspline_curve_with_knots(
    degree: usize,
    poles: Vec<GpPnt>,
    weights: Option<Vec<f64>>,
    mults: &[usize],
    knots: &[f64],
    closed: bool,
) -> Result<GeomBSplineCurve, String> {
    let d = bspline_descriptor(degree, poles.len(), mults, knots)?;
    // 5. drop the extra poles / weights (`cxx:847-871`).
    let weights = match weights {
        Some(w) if w.len() >= d.pole_range.end => Some(w[d.pole_range.clone()].to_vec()),
        Some(_) => return Err("B_SPLINE_CURVE: weight count mismatch".into()),
        None => None,
    };
    let poles = poles[d.pole_range.clone()].to_vec();
    let (uknots, umults) = (&d.uknots, &d.umults);

    let mut curve = if d.should_be_periodic {
        // `new TBSplineCurve(Poles, [Weights,] UniqueKnots, UniqueMults, Degree,
        // true)`: `CheckCurveData` (`Geom_BSplineCurve.cxx:91-94`) requires
        // `NbPoles(Degree, true, Mults) == Poles.Length()`.
        if occt_core::bspl::knots::nb_poles(degree as i32, true, umults) as usize != poles.len() {
            return Err("B_SPLINE_CURVE: periodic pole/degree mismatch".into());
        }
        let flat = occt_core::bspl::knots::knot_sequence_periodic(uknots, umults, degree as i32);
        GeomBSplineCurve { poles, weights, knots: flat, degree, periodic: true }
    } else {
        let flat = occt_core::bspl::banded_interp::knot_sequence(uknots, umults, degree as i32);
        match weights {
            Some(w) => {
                GeomBSplineCurve::rational(poles, w, flat, degree).map_err(|e| e.to_string())?
            }
            None => GeomBSplineCurve::new(poles, flat, degree).map_err(|e| e.to_string())?,
        }
    };

    // 6. force periodicity on closed curves (`cxx:920-926`).
    if closed && degree > 1 {
        let sp = curve.d0(curve.first_parameter());
        let ep = curve.d0(curve.last_parameter());
        if sp.square_distance(&ep) <= occt_core::precision::COMPUTATIONAL {
            curve.set_periodic();
        }
    }
    Ok(curve)
}

/// The 2D arm of `StepToGeom::MakeBSplineCurveCommon`
/// (`StepToGeom::MakeBSplineCurve2d`, `StepToGeom.cxx:952-963`): identical
/// knot/multiplicity surgery, built into `Geom2d_BSplineCurve` (non-rational:
/// this port's 2D B-spline stores `xs`/`ys` without weights).
fn make_bspline_curve_2d_with_knots(
    degree: usize,
    xs: Vec<f64>,
    ys: Vec<f64>,
    mults: &[usize],
    knots: &[f64],
    closed: bool,
) -> Result<Geom2dBSplineCurve, String> {
    let d = bspline_descriptor(degree, xs.len(), mults, knots)?;
    let xs = xs[d.pole_range.clone()].to_vec();
    let ys = ys[d.pole_range.clone()].to_vec();
    let (uknots, umults) = (&d.uknots, &d.umults);

    let mut curve = if d.should_be_periodic {
        if occt_core::bspl::knots::nb_poles(degree as i32, true, umults) as usize != xs.len() {
            return Err("B_SPLINE_CURVE_2D: periodic pole/degree mismatch".into());
        }
        let flat = occt_core::bspl::knots::knot_sequence_periodic(uknots, umults, degree as i32);
        Geom2dBSplineCurve { xs, ys, knots: flat, degree, periodic: true }
    } else {
        let flat = occt_core::bspl::banded_interp::knot_sequence(uknots, umults, degree as i32);
        Geom2dBSplineCurve::new(xs, ys, flat, degree).map_err(|e| e.to_string())?
    };

    if closed && degree > 1 {
        let sp = curve.d0(curve.first_parameter());
        let ep = curve.d0(curve.last_parameter());
        if sp.square_distance(&ep) <= occt_core::precision::COMPUTATIONAL {
            curve.set_periodic();
        }
    }
    Ok(curve)
}

/// The STEP `closed_curve` / `closed` flag of a curve record.
///
/// The resolver sees two layouts: a plain instance `(name, degree, poles,
/// curve_form, closed, …)` with `closed` at index 4, and a merged rational
/// complex `(name, degree, poles, weights, curve_form, closed, …)` with
/// `closed` at index 5 (`transfer.rs::merge_complex_body`). The weights slot
/// distinguishes them: the merged body puts `SELF` or a `(…)` list there, a
/// plain body a `.FORM.` marker.
fn curve_record_closed(rec: &Record) -> bool {
    let merged = rec
        .args
        .get(3)
        .map(|s| {
            let t = s.trim();
            t == "SELF" || t.starts_with('(')
        })
        .unwrap_or(false);
    let idx = if merged { 5 } else { 4 };
    parse_logical(rec.args.get(idx).map(|s| s.as_str()), false)
}

impl<'a> Resolver<'a> {

    pub(super) fn new(records: &'a HashMap<usize, Record>) -> Self {
        Self {
            records,
            precision: step_precision(records),
            plane_angle_factor: context_plane_angle_factor(records),
            length_factor: context_length_factor(records),
            b: TopoBuilder::new(),
            shape_cache: RefCell::new(HashMap::new()),
            vertex_bind: RefCell::new(HashMap::new()),
            point_cache: RefCell::new(HashMap::new()),
            dir_cache: RefCell::new(HashMap::new()),
            axis_cache: RefCell::new(HashMap::new()),
            curve_cache: RefCell::new(HashMap::new()),
            surface_cache: RefCell::new(HashMap::new()),
            curve2d_cache: RefCell::new(HashMap::new()),
            surface_curve_pcurves: RefCell::new(HashMap::new()),
            edge_curve_ref: RefCell::new(HashMap::new()),
            resolving: RefCell::new(HashSet::new()),
            warnings: RefCell::new(Vec::new()),
        }
    }

    pub(super) fn record(&self, id: usize) -> Result<&'a Record, String> {
        self.records
            .get(&id)
            .ok_or_else(|| format!("reference to undefined entity #{id}"))
    }

    pub(super) fn warn(&self, msg: String) {
        self.warnings.borrow_mut().push(msg);
    }

    pub(super) fn resolve_shape(&self, id: usize) -> Result<TopoShape, String> {
        if let Some(s) = self.shape_cache.borrow().get(&id) {
            return Ok(s.clone());
        }
        let rec = self.record(id)?;
        // `StepToTopoDS_TranslateTool::Bind` — a `VERTEX_POINT` another record was
        // bound onto (`StepToTopoDS_TranslateEdgeLoop.cxx:384-396`, `:466-477`)
        // resolves to that vertex instead of a fresh one.
        if rec.type_name == "VERTEX_POINT" {
            if let Some(v) = self.vertex_bind.borrow().get(&id) {
                return Ok(v.clone());
            }
        }
        if !self.resolving.borrow_mut().insert(id) {
            return Err(format!("cyclic entity reference #{id}"));
        }
        let shape = match rec.type_name.as_str() {
            "VERTEX_POINT" => self.resolve_vertex(rec),
            "EDGE_CURVE" => self.resolve_edge(rec),
            "ORIENTED_EDGE" => self.resolve_oriented_edge(rec),
            "EDGE_LOOP" => self.resolve_loop(rec),
            "VERTEX_LOOP" => self.resolve_vertex_loop(rec),
            "FACE_OUTER_BOUND" => self.resolve_outer_bound(rec),
            "FACE_BOUND" => self.resolve_outer_bound(rec),
            "ADVANCED_FACE" => self.resolve_face(rec),
            "CLOSED_SHELL" => self.resolve_shell(rec),
            "MANIFOLD_SOLID_BREP" => self.resolve_solid(rec),
            "BREP_WITH_VOIDS" => self.resolve_brep_with_voids(rec),
            "MAPPED_ITEM" => self.resolve_mapped_item(rec),
            other => {
                self.warn(format!("unsupported topological entity {other} (#{id})"));
                Err(format!("unsupported entity {other} (#{id})"))
            }
        };
        self.resolving.borrow_mut().remove(&id);
        match shape {
            Ok(s) => {
                self.shape_cache.borrow_mut().insert(id, s.clone());
                Ok(s)
            }
            Err(e) => Err(e),
        }
    }

    pub(super) fn resolve_vertex(&self, rec: &'a Record) -> Result<TopoShape, String> {
        if rec.args.len() < 2 {
            return Err("VERTEX_POINT: bad args".into());
        }
        let pid = parse_ref(&rec.args[1]).ok_or("VERTEX_POINT: bad point ref")?;
        let p = self.resolve_point(pid)?;
        Ok(self.b.make_vertex(p, 0.0).0)
    }

    pub(super) fn resolve_edge(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let start_ref = parse_ref(&rec.args[1]).ok_or("EDGE_CURVE: bad start ref")?;
        let end_ref = parse_ref(&rec.args[2]).ok_or("EDGE_CURVE: bad end ref")?;
        let curve_ref = parse_ref(&rec.args[3]).ok_or("EDGE_CURVE: bad curve ref")?;
        // `StepToTopoDS_TranslateEdge.cxx:290-322`: `same_sense` picks which
        // VERTEX is the edge's geometric first / last. `.T.` keeps
        // `edge_start -> edge_end`; `.F.` swaps them, so the projected range
        // and the vertex children both follow the curve's own direction.
        let same_sense = parse_logical(rec.args.get(4).map(String::as_str), true);
        let (first_ref, last_ref) = if same_sense {
            (start_ref, end_ref)
        } else {
            (end_ref, start_ref)
        };
        let v1 = self.resolve_shape(first_ref)?;
        let v2 = self.resolve_shape(last_ref)?;
        if !v1.is_vertex() || !v2.is_vertex() {
            return Err("EDGE_CURVE: endpoints are not vertices".into());
        }
        let curve = self.resolve_curve(curve_ref)?;
        // `MakeFromCurve3D` (`TranslateEdge.cxx:433-441`): the projection points
        // are the two endpoint vertices; when the two translated vertices are the
        // same `TopoDS_Shape` (`V1.IsSame(V2)`), `GetCartesianPoints` re-reads
        // them from the EDGE_CURVE's own `edge_start` / `edge_end` vertex geometry
        // (STEP order, i.e. **not** the `same_sense`-swapped order).
        let mut p1 = GeometryRegistry::global().vertex_point(&v1);
        let mut p2 = GeometryRegistry::global().vertex_point(&v2);
        if crate::topo_tools_full::is_same(&v1, &v2) {
            let raw1 = self.resolve_shape(start_ref)?;
            let raw2 = self.resolve_shape(end_ref)?;
            p1 = GeometryRegistry::global().vertex_point(&raw1);
            p2 = GeometryRegistry::global().vertex_point(&raw2);
        }
        let (curve, first, last) = edge_from_curve3d(curve, &p1, &p2);
        // `MakeFromCurve3D` (`TranslateEdge.cxx:452-455`): distance at projected
        // params after `UpdateParam3d` (and any displaced-Line shift).
        let temp1 = curve.d0(first).distance(&p1);
        let temp2 = curve.d0(last).distance(&p2);
        let mut e = self.b.make_edge(curve, first, last);
        // `BRep_Builder` / `BRepLib_MakeEdge`: first vertex FORWARD, last REVERSED
        // (`TopoDS_Builder::Add`). Needed so `TopExp::LastVertex` exists on a
        // closed EDGE_CURVE (same TVertex stored twice) and
        // `ShapeAnalysis_Edge::FirstVertex` works on REVERSED seam uses.
        let v1 = Vertex(v1);
        let v2 = Vertex(v2);
        self.b.add_edge_vertices(&mut e, &v1, &v2);
        // `MakeFromCurve3D` (`TranslateEdge.cxx:478-479`): grow vertex tolerance
        // to cover Project residual (`BRep_Builder::UpdateVertex` only raises).
        // Enabled after ExtPC Project temps on Shape f7 seams are ~1e-7..1e-13
        // (well under preci=1e-3); invent sampler residuals previously densified.
        let t1 = crate::brep_tool::BRepTool::vertex_tolerance(&v1).max(1.000001 * temp1);
        let t2 = crate::brep_tool::BRepTool::vertex_tolerance(&v2).max(1.000001 * temp2);
        v1.set_tolerance(t1);
        v2.set_tolerance(t2);
        // Remember the edge's curve entity for the face-level pcurve association
        // (a SURFACE_CURVE's pcurve is matched to the face's surface by ref).
        self.edge_curve_ref
            .borrow_mut()
            .insert(Arc::as_ptr(&e.0.tshape) as usize, curve_ref);
        Ok(e.0)
    }

    pub(super) fn resolve_oriented_edge(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let edge_ref = parse_ref(&rec.args[3]).ok_or("ORIENTED_EDGE: bad edge ref")?;
        let mut s = self.resolve_shape(edge_ref)?;
        // `StepToTopoDS_TranslateEdgeLoop.cxx:546-556`: the occurrence
        // orientation is the ORIENTED_EDGE boolean combined with the
        // EDGE_CURVE's `same_sense`:
        //   FORWARD  when (Orientation && SameSense) || (!Orientation && !SameSense)
        //   REVERSED otherwise.
        // The edge's stored forward direction already follows `same_sense`
        // (`resolve_edge`), so this is the XNOR of the two booleans. Without it
        // an `EDGE_CURVE(...,.F.)` occurrence comes back backwards.
        let ori = parse_logical(rec.args.get(4).map(String::as_str), true);
        let ec_same_sense = self.edge_curve_same_sense(edge_ref);
        s.set_orientation(if ori == ec_same_sense {
            Orientation::Forward
        } else {
            Orientation::Reversed
        });
        Ok(s)
    }

    /// `same_sense` of the `EDGE_CURVE` an ORIENTED_EDGE references.
    /// `StepToTopoDS_TranslateEdgeLoop.cxx:306-308` follows a nested
    /// ORIENTED_EDGE down to its EDGE_CURVE (bug #29979). A chain that does not
    /// end at an EDGE_CURVE defaults to `.T.`.
    fn edge_curve_same_sense(&self, mut id: usize) -> bool {
        for _ in 0..8 {
            let Ok(rec) = self.record(id) else {
                return true;
            };
            match rec.type_name.as_str() {
                "EDGE_CURVE" => {
                    return parse_logical(rec.args.get(4).map(String::as_str), true);
                }
                "ORIENTED_EDGE" => match parse_ref(&rec.args[3]) {
                    Some(next) => id = next,
                    None => return true,
                },
                _ => return true,
            }
        }
        true
    }

    pub(super) fn resolve_loop(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let items = parse_ref_list(&rec.args[1]);
        // `StepToTopoDS_TranslateEdgeLoop.cxx:288-403` and `:405-491` bind
        // (confuse) the loop's vertices before the oriented edges are mapped, so a
        // loop whose edges reference distinct `VERTEX_POINT` entities at the same
        // point — or adjacent edges that share no vertex at all — still comes out
        // connected.
        self.bind_edge_loop_vertices(&items);
        let mut edges = Vec::with_capacity(items.len());
        for &it in &items {
            let s = self.resolve_shape(it)?;
            if !s.is_edge() {
                return Err(format!("#{it}: expected EDGE in EDGE_LOOP"));
            }
            edges.push(Edge(s));
        }
        // Keep EDGE_LOOP order for seam SelectForwardSeam (`TranslateEdgeLoop`
        // iterates OrientedEdges in list order; last UpdateEdge wins). Reorder
        // for connectivity after pcurve association in `resolve_face`.
        Ok(self.b.make_wire(&edges).0)
    }

    /// The `EDGE_CURVE` record an `ORIENTED_EDGE` reaches — following a nested
    /// `ORIENTED_EDGE` (bug #29979, `StepToTopoDS_TranslateEdgeLoop.cxx:301-307`)
    /// — together with the occurrence `Orientation` (`:429-433`).
    fn oriented_edge_parts(&self, id: usize) -> Option<(usize, bool)> {
        let mut cur = id;
        for _ in 0..8 {
            let rec = self.record(cur).ok()?;
            if rec.type_name != "ORIENTED_EDGE" {
                return None;
            }
            let ori = parse_logical(rec.args.get(4).map(String::as_str), true);
            let next = parse_ref(&rec.args[3])?;
            let nrec = self.record(next).ok()?;
            if nrec.type_name == "ORIENTED_EDGE" {
                cur = next;
                continue;
            }
            if nrec.type_name != "EDGE_CURVE" {
                return None;
            }
            return Some((next, ori));
        }
        None
    }

    /// The two vertex-binding passes of `StepToTopoDS_TranslateEdgeLoop`:
    ///
    /// 1. `cxx:288-403` (bug PRO7656): for every oriented edge whose
    ///    `EDGE_CURVE`'s start/end vertices translate to points within
    ///    `Precision::Confusion()`, bind one STEP vertex onto the other's
    ///    `TopoDS_Vertex` — `Vend → V1` when `Vend` was not yet bound, else
    ///    `Vstart → V2` when `Vstart` was not, else `Vend → V1` again
    ///    (`cxx:384-396`).
    /// 2. `cxx:405-491` (bug BUC50070 #3815): for every adjacent pair `(j, j+1)`
    ///    whose *meeting* vertices (`OrEdge1->Orientation() ? EdgeEnd : EdgeStart`
    ///    and the mirror for edge 2, `cxx:429-433`) are different STEP entities
    ///    and translate to points within `Precision()` (`cxx:464-478`), bind
    ///    `Vs1 → V2` when `EC1` is not yet translated, else `Vs2 → V1` when `EC2`
    ///    is not.
    ///
    /// The port's "translate tool" is `shape_cache` (an entity is bound once
    /// `resolve_shape` has cached it), and the binding itself is
    /// [`Resolver::vertex_bind`], consulted by `resolve_shape` for `VERTEX_POINT`.
    fn bind_edge_loop_vertices(&self, oriented: &[usize]) {
        let point_of = |v: &TopoShape| BRepTool::vertex_point(&Vertex(v.clone()));
        // Pass 1.
        for &oe in oriented {
            let Some((ec_id, _)) = self.oriented_edge_parts(oe) else {
                continue;
            };
            let Ok(ec) = self.record(ec_id) else { continue };
            let (Some(a), Some(b)) = (parse_ref(&ec.args[1]), parse_ref(&ec.args[2])) else {
                continue;
            };
            let same_sense = parse_logical(ec.args.get(4).map(String::as_str), true);
            // `cxx:355-365`.
            let (vstart, vend) = if same_sense { (a, b) } else { (b, a) };
            // `cxx:367-368`: the pre-call bound state.
            let ist_v = self.shape_cache.borrow().contains_key(&vstart);
            let ise_v = self.shape_cache.borrow().contains_key(&vend);
            let (Ok(v1), Ok(v2)) = (self.resolve_shape(vstart), self.resolve_shape(vend)) else {
                continue;
            };
            if point_of(&v1).distance(&point_of(&v2)) <= occt_core::precision::CONFUSION {
                let (from, to) = if !ise_v {
                    (vend, v1)
                } else if !ist_v {
                    (vstart, v2)
                } else {
                    (vend, v1)
                };
                self.vertex_bind.borrow_mut().insert(from, to);
            }
        }
        // Pass 2.
        let n = oriented.len();
        for j in 0..n {
            if n < 2 {
                break;
            }
            let Some((ec1_id, ori1)) = self.oriented_edge_parts(oriented[j]) else {
                continue;
            };
            let Some((ec2_id, ori2)) = self.oriented_edge_parts(oriented[(j + 1) % n]) else {
                continue;
            };
            let (Ok(ec1), Ok(ec2)) = (self.record(ec1_id), self.record(ec2_id)) else {
                continue;
            };
            let (Some(s1), Some(e1)) = (parse_ref(&ec1.args[1]), parse_ref(&ec1.args[2])) else {
                continue;
            };
            let (Some(s2), Some(e2)) = (parse_ref(&ec2.args[1]), parse_ref(&ec2.args[2])) else {
                continue;
            };
            // `cxx:429-433`.
            let vs1 = if ori1 { e1 } else { s1 };
            let vs2 = if ori2 { s2 } else { e2 };
            let vs11 = if ori1 { s1 } else { e1 };
            let vs22 = if ori2 { e2 } else { s2 };
            // `cxx:435-438`: already the same STEP vertex on one of the four ends.
            if vs1 == vs2 || vs1 == vs22 || vs2 == vs11 || vs22 == vs11 {
                continue;
            }
            let (Ok(v1), Ok(v2)) = (self.resolve_shape(vs1), self.resolve_shape(vs2)) else {
                continue;
            };
            if Arc::ptr_eq(&v1.tshape, &v2.tshape) {
                continue;
            }
            if point_of(&v1).distance(&point_of(&v2)) <= self.precision {
                // `cxx:466-477`.
                if !self.shape_cache.borrow().contains_key(&ec1_id) {
                    self.vertex_bind.borrow_mut().insert(vs1, v2);
                } else if !self.shape_cache.borrow().contains_key(&ec2_id) {
                    self.vertex_bind.borrow_mut().insert(vs2, v1);
                }
            }
        }
    }

    /// A `VERTEX_LOOP(name, vertex)` is the boundary of a degenerate face — a
    /// loop reduced to a single vertex (sphere pole, cone apex). It maps to a
    /// wire containing no edges (the face still references it via FACE_BOUND).
    pub(super) fn resolve_vertex_loop(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let _ = parse_ref(&rec.args[1]).ok_or("VERTEX_LOOP: bad vertex ref")?;
        Ok(self.b.make_wire(&[]).0)
    }

    pub(super) fn resolve_outer_bound(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let loop_ref = parse_ref(&rec.args[1]).ok_or("FACE_OUTER_BOUND: bad loop ref")?;
        let mut s = self.resolve_shape(loop_ref)?;
        if !s.is_wire() {
            return Err("FACE_OUTER_BOUND: loop is not a wire".into());
        }
        // `StepShape_FaceBound::Orientation`: `.F.` reverses the wire
        // (`StepToTopoDS_TranslateFace.cxx` PolyLoop / EdgeLoop bind).
        if !parse_logical(rec.args.get(2).map(String::as_str), true) {
            s.set_orientation(Orientation::Reversed);
        }
        Ok(s)
    }

    pub(super) fn resolve_face(&self, rec: &'a Record) -> Result<TopoShape, String> {
        // Two argument layouts occur in the wild:
        //  * STEP-214 (ISO standard, written by OCCT/FreeCAD and by this port's
        //    writer since the cylinder-validation fix): surface is the third
        //    argument — `ADVANCED_FACE(name, bounds, surface, same_sense)`.
        //  * legacy files from this port's earlier writer: surface is the
        //    second argument — `ADVANCED_FACE('', #surface, (bounds), .T.)`.
        // Detect by checking which argument holds a surface reference.
        let surf_ref = if let Some(r) = parse_ref(&rec.args[2]) {
            r
        } else {
            parse_ref(&rec.args[1]).ok_or("ADVANCED_FACE: bad surface ref")?
        };
        let surface = self.resolve_surface(surf_ref)?;
        // Bounds are the other of the two leading arguments.
        let bounds = if rec.args.get(2).map(|s| s.starts_with('(')).unwrap_or(false) {
            parse_ref_list(&rec.args[2])
        } else {
            parse_ref_list(&rec.args[1])
        };
        // `StepToTopoDS_TranslateFace.cxx:572-573` / `716-752`.
        let same_sense = parse_logical(rec.args.last().map(String::as_str), true);
        let a_same_sense = if self.step_surface_is_reversed(surf_ref) {
            !same_sense
        } else {
            same_sense
        };
        // `StepToTopoDS_TranslateFace.cxx:608-634` — a lone VertexLoop on a
        // sphere / BSpline / revolution is the whole closed surface; add
        // `BRepLib_MakeFace` natural bounds and skip the vertex loop itself.
        // `cxx:629` passes `Precision()`, the translator's `myPrecision`
        // (`STEPControl_ActorRead.cxx:2370-2384`) = file uncertainty x
        // `LengthFactor` (`step_precision`), not `Precision::Confusion()`.
        if bounds.len() == 1
            && self.bound_loop_is_vertex_loop(bounds[0])
            && (surface.gp_sphere().is_some()
                || surface.is_bspline_surface()
                || surface.is_surface_of_revolution())
        {
            let mut face = crate::brep_lib_make_face::make_face_from_surface(
                surface,
                self.precision,
            );
            face.0.set_orientation(if same_sense {
                Orientation::Forward
            } else {
                Orientation::Reversed
            });
            return Ok(face.0);
        }
        let mut wires = Vec::with_capacity(bounds.len());
        // `TranslateEdgeLoop.cxx:272`: `ForwardWire = FaceBound->Orientation()`
        // for SelectForwardSeam Compose — captured before sameSense reverse.
        let mut wire_bound_oris = Vec::with_capacity(bounds.len());
        for &b in &bounds {
            let mut s = self.resolve_shape(b)?;
            if !s.is_wire() {
                return Err(format!("#{b}: expected wire in face bounds"));
            }
            wire_bound_oris.push(s.orientation());
            // Bound orientation is already on the wire from `resolve_outer_bound`.
            // A reversed `Face_Surface` (or negative-major torus) flips it again
            // so the stored wire matches CAS.CADE (`cxx:714-723`).
            if !a_same_sense {
                s.reverse();
            }
            wires.push(Wire(s));
        }
        let mut face = self.b.make_face(surface, &wires);
        face.0.set_orientation(if same_sense {
            Orientation::Forward
        } else {
            Orientation::Reversed
        });
        let face_key = GeometryRegistry::shape_key(&face.0);
        // Associate each wire edge's SURFACE_CURVE pcurve with this face's
        // surface (`BRep_Builder::UpdateEdge(edge, pcurve, face, tol)`).
        let face_ori = face.0.orientation();
        for (w, &wire_o) in wires.iter_mut().zip(wire_bound_oris.iter()) {
            // EdgeO is ORIENTED_EDGE storage (`cxx:700`), not wire-composed.
            // Iterate EDGE_LOOP order so last seam UpdateEdge matches cxx.
            let loop_edges = edges_stored_on_wire(w);
            // `StepToTopoDS_GeometricTool.cxx:104-116`: the seam test counts how
            // many oriented edges of the loop reference the same step edge.
            let loop_refs: Vec<usize> = loop_edges
                .iter()
                .map(|x| {
                    self.edge_curve_ref
                        .borrow()
                        .get(&(Arc::as_ptr(&x.0.tshape) as usize))
                        .copied()
                        .unwrap_or(0)
                })
                .collect();
            for e in &loop_edges {
                let key = Arc::as_ptr(&e.0.tshape) as usize;
                let curve_ref = self.edge_curve_ref.borrow().get(&key).copied().unwrap_or(0);
                let nb_oe = loop_refs.iter().filter(|&&r| r == curve_ref).count();
                self.associate_edge_pcurve(e, surf_ref, face_key, face_ori, wire_o, nb_oe)?;
            }
            // `TranslateEdgeLoop.cxx:844-868` EdgeProjAux before CheckPCurves.
            // cxx:815 is only `B.Add(W,E)`; after CheckPCurves (cxx:875) the
            // function returns — no OrientEdgesOnWire / wire child reorder.
            // `TranslateEdgeLoop.cxx:236` `preci = Precision()`, set by
            // `STEPControl_ActorRead` from the file's
            // `GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT` uncertainty or from
            // `read.precision.val` (default 1.e-03); see `step_precision`.
            let preci = self.precision;
            // `StepToTopoDS_TranslateEdgeLoop.cxx:875` -> `CheckPCurves`
            // (`:105-177`) -> `XSAlgo_ShapeProcessor::CheckPCurve`
            // (`XSAlgo_ShapeProcessor.cxx:344-401`): the reader's "advanced check"
            // drops a pcurve whose U/V span wraps the surface or whose end
            // points disagree with the 3D curve beyond the read precision; the
            // dropped pcurve is re-projected by `ShapeFix_Edge::FixAddPCurve`
            // in the `ShapeFix` pass below (OCCT's `XSAlgo` pipeline).
            for e in &loop_edges {
                crate::shhealing::xsalgo_check_pcurve(e, &face, self.precision);
            }
            crate::shhealing::project_wire_pcurve_ranges(w, &face, preci);
            crate::shhealing::check_pcurves_and_shift(w, &face, preci);
        }
        // `BRepLib_MakeFace.cxx:860-866` forced SameParameter is already in
        // `make_face_uv` for natural-bound Offset faces. Wiring it on STEP
        // Offset faces (after EdgeProjAux) densifies Shape 6141 -> 6213 vs
        // occ 6150; TranslateFace.cxx has no equivalent call.
        Ok(face.0)
    }

    /// True when `bound_id` is a `FACE_BOUND` / `FACE_OUTER_BOUND` whose loop
    /// is a `VERTEX_LOOP` (`StepToTopoDS_TranslateFace.cxx:608`).
    fn bound_loop_is_vertex_loop(&self, bound_id: usize) -> bool {
        let Ok(rec) = self.record(bound_id) else {
            return false;
        };
        match rec.type_name.as_str() {
            "FACE_BOUND" | "FACE_OUTER_BOUND" => parse_ref(&rec.args[1])
                .and_then(|id| self.record(id).ok())
                .map(|r| r.type_name == "VERTEX_LOOP")
                .unwrap_or(false),
            "VERTEX_LOOP" => true,
            _ => false,
        }
    }

    /// `StepToTopoDS_TranslateFace.cxx:479-491` — SolidWorks torus with a
    /// negative major radius is treated as a reversed `Face_Surface`.
    pub(super) fn step_surface_is_reversed(&self, surf_id: usize) -> bool {
        let Ok(rec) = self.record(surf_id) else {
            return false;
        };
        match rec.type_name.as_str() {
            "RECTANGULAR_TRIMMED_SURFACE" => rec
                .args
                .get(1)
                .and_then(|s| parse_ref(s))
                .map(|basis| self.step_surface_is_reversed(basis))
                .unwrap_or(false),
            "TOROIDAL_SURFACE" => rec
                .args
                .get(2)
                .and_then(|s| parse_f64(s).ok())
                .map(|maj| maj < 0.0)
                .unwrap_or(false),
            _ => false,
        }
    }

    /// True when `assoc_ref` is a `PCURVE` whose basis surface is `surf_ref`
    /// (`StepToTopoDS_GeometricTool.cxx:57-63` / `:94-97` / `:143-146`).
    fn pcurve_lies_on(&self, assoc_ref: usize, surf_ref: usize) -> bool {
        matches!(self.resolve_pcurve(assoc_ref), Ok((basis, _)) if basis == surf_ref)
    }

    /// The raw STEP origin and direction ratios of a `PCURVE` whose underlying
    /// 2D curve is a `LINE` (`StepToTopoDS_GeometricTool.cxx:167-181`), as
    /// `(Pnt.X, Pnt.Y, Dir.X, Dir.Y)`. `None` for any other entity.
    fn step_2d_line(&self, assoc_ref: usize) -> Option<(f64, f64, f64, f64)> {
        let rec = self.record(assoc_ref).ok()?;
        if rec.type_name != "PCURVE" {
            return None;
        }
        let dri = parse_ref(&rec.args[2])?;
        let dri_rec = self.record(dri).ok()?;
        if dri_rec.type_name != "DEFINITIONAL_REPRESENTATION" {
            return None;
        }
        let line_ref = *parse_ref_list(&dri_rec.args[1]).first()?;
        let line = self.record(line_ref).ok()?;
        if line.type_name != "LINE" {
            return None;
        }
        let pnt = self.record(parse_ref(&line.args[1])?).ok()?;
        let dir = self.record(parse_ref(&line.args[2])?).ok()?;
        if pnt.type_name != "CARTESIAN_POINT" || dir.type_name != "DIRECTION" {
            return None;
        }
        let p = parse_xy(&pnt.args[1]).ok()?;
        let d = parse_xy(&dir.args[1]).ok()?;
        Some((p.0, p.1, d.0, d.1))
    }

    /// `StepToTopoDS_GeometricTool::IsSeamCurve` (`cxx:78-121`). `nb_oe` is how
    /// many oriented edges of the loop reference the same step edge.
    fn is_seam_curve(&self, curve_ref: usize, surf_ref: usize, nb_oe: usize) -> bool {
        if let Ok(rec) = self.record(curve_ref) {
            if rec.type_name == "SEAM_CURVE" {
                return true;
            }
        }
        let assoc = self
            .surface_curve_pcurves
            .borrow()
            .get(&curve_ref)
            .cloned()
            .unwrap_or_default();
        if assoc.len() != 2 {
            return false;
        }
        if !self.pcurve_lies_on(assoc[0], surf_ref) || !self.pcurve_lies_on(assoc[1], surf_ref) {
            return false;
        }
        // `cxx:104-116`: two oriented edges of the same wire share this edge.
        nb_oe == 2
    }

    /// `StepToTopoDS_GeometricTool::IsLikeSeam` (`cxx:133-218`): the two pcurves
    /// lie on the same surface but the edge is used once by the loop, and both
    /// are `LINE`s sharing an origin coordinate and a direction (CATIA BRep).
    fn is_like_seam(&self, curve_ref: usize, surf_ref: usize, nb_oe: usize) -> bool {
        let assoc = self
            .surface_curve_pcurves
            .borrow()
            .get(&curve_ref)
            .cloned()
            .unwrap_or_default();
        if assoc.len() != 2 {
            return false;
        }
        if !self.pcurve_lies_on(assoc[0], surf_ref) || !self.pcurve_lies_on(assoc[1], surf_ref) {
            return false;
        }
        // `cxx:160-165`: the two oriented edges are not in the same wire.
        if nb_oe != 1 {
            return false;
        }
        let (Some(l1), Some(l2)) = (self.step_2d_line(assoc[0]), self.step_2d_line(assoc[1]))
        else {
            return false;
        };
        // `cxx:183-196`, `preci2d = Precision::PConfusion()`.
        let preci2d = occt_core::precision::PCONFUSION;
        let delta_x = (l1.0 - l2.0).abs();
        let delta_y = (l1.1 - l2.1).abs();
        let delta_dir_x = (l1.2 - l2.2).abs();
        let delta_dir_y = (l1.3 - l2.3).abs();
        if delta_x < preci2d || delta_y < preci2d {
            delta_dir_x < preci2d && delta_dir_y < preci2d
        } else {
            false
        }
    }

    /// `TranslateEdgeLoop.cxx:699-734`: order the two seam pcurves so the first
    /// entry is the FORWARD one. `ShapeAnalysis_Curve::SelectForwardSeam`, then
    /// flip when the cumulative edge/wire/face orientation is reversed.
    fn order_seam_pcurves(
        &self,
        edge: &Edge,
        matched: &mut [Arc<dyn Curve2d>],
        face_ori: Orientation,
        wire_o: Orientation,
    ) {
        if matched.len() != 2 {
            return;
        }
        // WireO is FaceBound->Orientation (`cxx:272`), not the post-sameSense wire.
        let mut fwd =
            crate::pcurve_full::select_forward_seam(matched[0].as_ref(), matched[1].as_ref());
        if fwd != 0 {
            let edge_o = edge.0.orientation();
            let cumul = Orientation::compose(edge_o, wire_o);
            let cumul = Orientation::compose(cumul, face_ori);
            if cumul != Orientation::Forward {
                fwd = 3 - fwd;
            }
            if fwd == 2 {
                matched.swap(0, 1);
            }
        }
    }

    /// Match a wire edge's SURFACE_CURVE pcurve to `surf_ref` and attach its 2D
    /// curve to the edge for `face_key`. Source:
    /// `StepToTopoDS_GeometricTool::PCurve` + `StepToTopoDS_TranslateEdge::MakePCurve`.
    ///
    /// `nb_oe` is the number of oriented edges of the current loop that
    /// reference this edge; it selects the seam arm
    /// (`TranslateEdgeLoop.cxx:611-667`, `:689-767`, `:779-796`).
    pub(super) fn associate_edge_pcurve(
        &self,
        edge: &Edge,
        surf_ref: usize,
        face_key: usize,
        face_ori: Orientation,
        wire_o: Orientation,
        nb_oe: usize,
    ) -> Result<(), String> {
        let key = Arc::as_ptr(&edge.0.tshape) as usize;
        let Some(&curve_ref) = self.edge_curve_ref.borrow().get(&key) else {
            return Ok(());
        };
        let Some(pcurves) = self.surface_curve_pcurves.borrow().get(&curve_ref).cloned() else {
            return Ok(());
        };
        // Each associated_geometry entry is a pcurve_or_surface: a PCURVE (with
        // a 2D curve) or a SURFACE (intersection curve). Only PCURVEs carry the
        // pcurve; a SURFACE entry fails to resolve and is skipped. The raw 2D
        // curve is stored with its own parameterization; the meshing maps the
        // edge's 3D parameter onto it (`BRepMesh_EdgeParameterProvider`).
        let mut matched: Vec<Arc<dyn Curve2d>> = Vec::new();
        for pc in pcurves {
            if let Ok((basis_surf, c2d)) = self.resolve_pcurve(pc) {
                if basis_surf == surf_ref {
                    matched.push(c2d);
                }
            }
        }
        let is_seam = self.is_seam_curve(curve_ref, surf_ref, nb_oe);
        let is_like_seam = !is_seam && self.is_like_seam(curve_ref, surf_ref, nb_oe);
        let stored: Vec<Arc<dyn Curve2d>> = if is_seam {
            // `cxx:736-767`: a seam edge carries two pcurves on the same face
            // (one per side of the seam, e.g. the cone's `u = 0` and `u = 2pi`),
            // stored forward-then-reversed so the wire assembly can hand each
            // traversal its own side.
            self.order_seam_pcurves(edge, &mut matched, face_ori, wire_o);
            matched
        } else if is_like_seam {
            // `cxx:748-767`: `else UpdateEdge(E, C2d2, Face, 0.)` - only the
            // forward pcurve is stored for CATIA-like seams.
            self.order_seam_pcurves(edge, &mut matched, face_ori, wire_o);
            matched.into_iter().take(1).collect()
        } else {
            // `cxx:645-667` walks every matching pcurve and keeps the last, then
            // `cxx:786-793` does `B.UpdateEdge(E, C2d, Face, 0.)` with it.
            matched.pop().into_iter().collect()
        };
        if !stored.is_empty() {
            GeometryRegistry::global().set_edge_pcurves(&edge.0, face_key, stored);
        }
        Ok(())
    }

    pub(super) fn resolve_shell(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let faces = parse_ref_list(&rec.args[1]);
        let mut fs = Vec::with_capacity(faces.len());
        for &it in &faces {
            let s = self.resolve_shape(it)?;
            if !s.is_face() {
                return Err(format!("#{it}: expected FACE in CLOSED_SHELL"));
            }
            fs.push(Face(s));
        }
        Ok(self.b.make_shell(&fs).0)
    }

    pub(super) fn resolve_solid(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let outer = parse_ref(&rec.args[1]).ok_or("MANIFOLD_SOLID_BREP: bad shell ref")?;
        let s = self.resolve_shape(outer)?;
        if !s.is_shell() {
            return Err("MANIFOLD_SOLID_BREP: outer is not a shell".into());
        }
        Ok(self.b.make_solid(&[Shell(s)]).0)
    }

    /// `StepToTopoDS_Builder::Init(BrepWithVoids)` (`StepToTopoDS_Builder.cxx:
    /// 179-251`), reached from `STEPControl_ActorRead::TransferEntity`
    /// (`cxx:1827-1831`, checked before `MANIFOLD_SOLID_BREP` because
    /// `BREP_WITH_VOIDS` is a subtype). The outer `CLOSED_SHELL` becomes the
    /// solid's first shell; every `ORIENTED_CLOSED_SHELL` of the `voids` set
    /// becomes an inner shell, reversed when its `orientation` attribute is
    /// `.F.` (`cxx:237-241`).
    pub(super) fn resolve_brep_with_voids(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let outer_ref = parse_ref(&rec.args[1]).ok_or("BREP_WITH_VOIDS: bad outer ref")?;
        let outer = self.resolve_shape(outer_ref)?;
        if !outer.is_shell() {
            // `cxx:210-214`: "OuterShell from BrepWithVoids not mapped to TopoDS".
            return Err("BREP_WITH_VOIDS: outer is not a shell".into());
        }
        let mut shells = vec![Shell(outer)];
        for void in parse_ref_list(&rec.args[2]) {
            match self.resolve_oriented_closed_shell(void) {
                Ok(s) => shells.push(Shell(s)),
                // `cxx:244-247`: "A Void from BrepWithVoids not mapped to
                // TopoDS" is a warning; the solid keeps the shells that mapped.
                Err(e) => self.warn(format!("void #{void}: {e}")),
            }
        }
        Ok(self.b.make_solid(&shells).0)
    }

    /// One void of a `BREP_WITH_VOIDS`: `ORIENTED_CLOSED_SHELL(name, *,
    /// closed_shell, orientation)`. Its faces are the derived
    /// `SELF\connected_face_set.cfs_faces := closed_shell.cfs_faces`
    /// (`StepShape_OrientedClosedShell.cxx:81-87`), so the shell is translated
    /// through the referenced `CLOSED_SHELL`.
    fn resolve_oriented_closed_shell(&self, id: usize) -> Result<TopoShape, String> {
        let rec = self.record(id)?;
        if rec.type_name != "ORIENTED_CLOSED_SHELL" {
            return Err(format!("expected ORIENTED_CLOSED_SHELL at #{id}"));
        }
        let shell_ref = parse_ref(&rec.args[2]).ok_or("ORIENTED_CLOSED_SHELL: bad shell ref")?;
        let mut s = self.resolve_shape(shell_ref)?;
        if !s.is_shell() {
            return Err(format!("#{shell_ref}: void is not a shell"));
        }
        // `StepToTopoDS_Builder.cxx:239-241`:
        // `if (!anOCShell->Orientation()) aShape.Reverse();`
        if !parse_logical(rec.args.get(3).map(String::as_str), true) {
            s.set_orientation(Orientation::Reversed);
        }
        Ok(s)
    }

    /// `STEPControl_ActorRead::TransferEntity(StepRepr_MappedItem)`
    /// (`STEPControl_ActorRead.cxx:1971-2045`): an assembly instance. The shape
    /// of `mapping_source.mapped_representation` is placed by the transform
    /// from `mapping_source.mapping_origin` to `mapping_target`, or, when the
    /// mapping target is a `CARTESIAN_TRANSFORMATION_OPERATOR_3D`, by the
    /// transform carried by that operator (`cxx:2013-2032`). `MAPPED_ITEM(name,
    /// mapping_source, mapping_target)` (`RWStepRepr_RWMappedItem.cxx:37-63`);
    /// `REPRESENTATION_MAP(mapping_origin, mapped_representation)`
    /// (`RWStepRepr_RWRepresentationMap.cxx:34-57`).
    pub(super) fn resolve_mapped_item(&self, rec: &'a Record) -> Result<TopoShape, String> {
        let source_ref = parse_ref(&rec.args[1]).ok_or("MAPPED_ITEM: bad mapping_source ref")?;
        let target_ref = parse_ref(&rec.args[2]).ok_or("MAPPED_ITEM: bad mapping_target ref")?;
        let source = self.record(source_ref)?;
        if source.type_name != "REPRESENTATION_MAP" {
            return Err(format!(
                "MAPPED_ITEM: mapping_source #{source_ref} is {}",
                source.type_name
            ));
        }
        let origin_ref =
            parse_ref(&source.args[0]).ok_or("REPRESENTATION_MAP: bad mapping_origin ref")?;
        let maprep_ref = parse_ref(&source.args[1])
            .ok_or("REPRESENTATION_MAP: bad mapped_representation ref")?;
        // `cxx:1989-2001`: transfer the mapped representation, then warn when it
        // produced no shape.
        let (_, shapes) = self.resolve_representation(maprep_ref)?;
        if shapes.is_empty() {
            return Err(format!(
                "MAPPED_ITEM: mapped representation #{maprep_ref} produced no shape"
            ));
        }
        let mut mapped = if shapes.len() == 1 {
            shapes.into_iter().next().unwrap()
        } else {
            TopoBuilder::new().make_compound_of(&shapes).0
        };
        // `cxx:2013-2032`: two placement formulas.
        let target = self.record(target_ref)?;
        let trsf = if target.type_name == "CARTESIAN_TRANSFORMATION_OPERATOR_3D" {
            self.make_transformation3d(target_ref).ok()
        } else {
            match (self.resolve_axis2(origin_ref), self.resolve_axis2(target_ref)) {
                (Ok(ax_orig), Ok(ax_targ)) => Some(self.compute_axis_transform(&ax_orig, &ax_targ)),
                _ => None,
            }
        };
        match trsf {
            // `cxx:2038 ApplyTransformation`: `shape.Move(TopLoc_Location(Trsf))`.
            // This port's geometry registry is location-blind (the mesh path
            // reads the registered geometry, not `TopoShape::location`), so the
            // placement is baked into a copy, the same mechanism `Assembly`
            // export uses (`shape_ops::transformed_copy`).
            Some(t) if t.form() != occt_core::gp::TrsfForm::Identity => {
                match crate::shape_ops::transformed_copy(&mapped, &t) {
                    Ok(copy) => mapped = copy,
                    Err(e) => self.warn(format!("MAPPED_ITEM: transform failed: {e}")),
                }
            }
            Some(_) => {}
            // `cxx:2042-2044`: "Mapped Item, case not recognized, location
            // ignored" - keep the untransformed shape.
            None => self.warn(format!(
                "MAPPED_ITEM #{target_ref}: case not recognized, location ignored"
            )),
        }
        Ok(mapped)
    }

    /// `StepToGeom::MakeTransformation3d` (`StepToGeom.cxx:2153-2214`) for a
    /// `CARTESIAN_TRANSFORMATION_OPERATOR_3D`. The reader starts at parameter 3
    /// (two `functionally_defined_transformation` names are skipped,
    /// `RWStepGeom_RWCartesianTransformationOperator.cxx:36-49`), so axis1 is
    /// param 4, axis2 param 5, local_origin param 6, scale param 7 and axis3
    /// param 8 (`RWStepGeom_RWCartesianTransformationOperator3d.cxx:82-128`).
    pub(super) fn make_transformation3d(&self, id: usize) -> Result<occt_core::gp::GpTrsf, String> {
        let rec = self.record(id)?;
        let origin_ref = parse_ref(&rec.args[5]).ok_or("CTO3D: bad local_origin ref")?;
        // `cxx:2158`: `MakeCartesianPoint(LocalOrigin)` - scaled by LengthFactor.
        let p = self.resolve_point(origin_ref)?;
        // `cxx:2163-2198`: default X / Y, overridden by axis1 / axis2.
        let mut d1 = dir_x();
        if let Some(r) = parse_ref(&rec.args[3]) {
            if let Ok(d) = self.resolve_direction(r) {
                d1 = d;
            }
        }
        let mut d2 = dir_y();
        if let Some(r) = parse_ref(&rec.args[4]) {
            if let Ok(d) = self.resolve_direction(r) {
                d2 = d;
            }
        }
        // `cxx:2186-2202`: axis3, or `D1.Crossed(D2)` when absent.
        let d3 = match parse_ref(&rec.args[7]) {
            Some(r) => self
                .resolve_direction(r)
                .unwrap_or_else(|_| d1.crossed(&d2).unwrap_or_else(|_| dir_z())),
            None => d1.crossed(&d2).unwrap_or_else(|_| dir_z()),
        };
        // `gp_Ax3(P, D3, D1)`: Vy = D3 ^ D1, X re-orthogonalized as Vy ^ D3.
        let ydir = d3.crossed(&d1).map_err(|e| format!("CTO3D: {e}"))?;
        let xdir = ydir.crossed(&d3).map_err(|e| format!("CTO3D: {e}"))?;
        let ax3 = GpAx3 {
            axis: GpAx1::new(p, d3),
            vxdir: xdir,
            vydir: ydir,
        };
        let mut t = occt_core::gp::GpTrsf::identity();
        t.set_transformation(&ax3);
        // `cxx:2205-2208`: `if (HasScale) CT.SetScaleFactor(Scale())`.
        if let Some(s) = rec.args.get(6).and_then(|a| parse_f64(a).ok()) {
            t.scale = s;
        }
        // `cxx:2210`: `CT = CT.Inverted()`.
        t.invert().map_err(|e| format!("CTO3D: {e}"))?;
        Ok(t)
    }

    /// `STEPControl_ActorRead::ComputeTransformation` axis-pair arm
    /// (`STEPControl_ActorRead.cxx:2471-2489`):
    /// `Trsf.SetTransformation(ax3Targ, ax3Orig)`
    /// (`gp_Trsf.cxx:172-192`) maps `ax3Orig`-frame coordinates into the
    /// `ax3Targ` frame, i.e. `world->ax3Orig` composed with the inverse of
    /// `world->ax3Targ`.
    pub(super) fn compute_axis_transform(&self, orig: &GpAx2, targ: &GpAx2) -> occt_core::gp::GpTrsf {
        let ax3_orig = GpAx3 {
            axis: orig.axis,
            vxdir: orig.vxdir,
            vydir: orig.vydir,
        };
        let ax3_targ = GpAx3 {
            axis: targ.axis,
            vxdir: targ.vxdir,
            vydir: targ.vydir,
        };
        let mut t_orig = occt_core::gp::GpTrsf::identity();
        t_orig.set_transformation(&ax3_orig);
        let mut t_targ = occt_core::gp::GpTrsf::identity();
        t_targ.set_transformation(&ax3_targ);
        let inv_targ = t_targ
            .inverted()
            .unwrap_or_else(|_| occt_core::gp::GpTrsf::identity());
        t_orig.multiplied(&inv_targ)
    }

    pub(super) fn resolve_point(&self, id: usize) -> Result<GpPnt, String> {
        if let Some(p) = self.point_cache.borrow().get(&id) {
            return Ok(*p);
        }
        let rec = self.record(id)?;
        if rec.type_name != "CARTESIAN_POINT" {
            self.warn(format!("expected CARTESIAN_POINT, got {} (#{id})", rec.type_name));
            return Err(format!("expected CARTESIAN_POINT at #{id}"));
        }
        let v = parse_xyz(&rec.args[1]).map_err(|e| format!("CARTESIAN_POINT #{id}: {e}"))?;
        // `StepToGeom::MakeCartesianPoint` (`StepToGeom.cxx:1173-1186`): every
        // coordinate is multiplied by `LF = theLocalFactors.LengthFactor()`
        // (`cxx:1179`). The 2D reader (`MakeCartesianPoint2d`, `cxx:1191-1204`)
        // deliberately does NOT scale, so `resolve_point_2d` stays unchanged.
        let p = GpPnt::from_xyz(&v.multiplied(self.length_factor));
        self.point_cache.borrow_mut().insert(id, p);
        Ok(p)
    }

    pub(super) fn resolve_direction(&self, id: usize) -> Result<GpDir, String> {
        if let Some(d) = self.dir_cache.borrow().get(&id) {
            return Ok(*d);
        }
        let rec = self.record(id)?;
        if rec.type_name != "DIRECTION" {
            self.warn(format!("expected DIRECTION, got {} (#{id})", rec.type_name));
            return Err(format!("expected DIRECTION at #{id}"));
        }
        let v = parse_xyz(&rec.args[1]).map_err(|e| format!("DIRECTION #{id}: {e}"))?;
        let d = GpDir::new(v.x, v.y, v.z).map_err(|e| format!("DIRECTION #{id}: {e}"))?;
        self.dir_cache.borrow_mut().insert(id, d);
        Ok(d)
    }

    pub(super) fn resolve_vector(&self, id: usize) -> Result<GpVec, String> {
        let rec = self.record(id)?;
        if rec.type_name != "VECTOR" {
            self.warn(format!("expected VECTOR, got {} (#{id})", rec.type_name));
            return Err(format!("expected VECTOR at #{id}"));
        }
        let dir_ref = parse_ref(&rec.args[1]).ok_or("VECTOR: bad direction ref")?;
        let mag = parse_f64(&rec.args[2])?;
        let dir = self.resolve_direction(dir_ref)?;
        // `StepToGeom::MakeVectorWithMagnitude` (`StepToGeom.cxx:2569-2581`):
        // `V = D->Dir().XYZ() * SV->Magnitude() * LengthFactor` (`cxx:2577`).
        // `MakeVectorWithMagnitude2d` (`cxx:2586-2597`) does not scale, so
        // `resolve_vector_2d` stays unchanged.
        Ok(GpVec::from_xyz(
            &dir.xyz().multiplied(mag * self.length_factor),
        ))
    }

    /// `AXIS1_PLACEMENT(name, location, axis_direction)` → `GpAx1`.
    pub(super) fn resolve_axis1(&self, id: usize) -> Result<GpAx1, String> {
        let rec = self.record(id)?;
        if rec.type_name != "AXIS1_PLACEMENT" {
            self.warn(format!("expected AXIS1_PLACEMENT, got {} (#{id})", rec.type_name));
            return Err(format!("expected AXIS1_PLACEMENT at #{id}"));
        }
        let loc_ref = parse_ref(&rec.args[1]).ok_or("AXIS1: bad location ref")?;
        let dir_ref = parse_ref(&rec.args[2]).ok_or("AXIS1: bad direction ref")?;
        let loc = self.resolve_point(loc_ref)?;
        let dir = self.resolve_direction(dir_ref)?;
        Ok(GpAx1::new(loc, dir))
    }

    pub(super) fn resolve_axis2(&self, id: usize) -> Result<GpAx2, String> {
        if let Some(a) = self.axis_cache.borrow().get(&id) {
            return Ok(*a);
        }
        let rec = self.record(id)?;
        if rec.type_name != "AXIS2_PLACEMENT_3D" {
            self.warn(format!(
                "expected AXIS2_PLACEMENT_3D, got {} (#{id})",
                rec.type_name
            ));
            return Err(format!("expected AXIS2_PLACEMENT_3D at #{id}"));
        }
        let loc = self
            .resolve_point(parse_ref(&rec.args[1]).ok_or("AXIS2: bad location ref")?)?;
        let axis = match parse_ref(&rec.args[2]) {
            Some(r) => self.resolve_direction(r)?,
            None => dir_z(),
        };
        let refd = match parse_ref(&rec.args[3]) {
            Some(r) => self.resolve_direction(r)?,
            None => dir_x(),
        };
        let ax2 = GpAx2::new(loc, axis, refd).map_err(|e| format!("AXIS2_PLACEMENT_3D: {e}"))?;
        self.axis_cache.borrow_mut().insert(id, ax2);
        Ok(ax2)
    }

    /// One side of `StepToGeom::ExtractParameter` (`StepToGeom.cxx:2221-2317`):
    /// the master representation decides which select is authoritative; a
    /// parameter select becomes `Shift + Factor * value` (`cxx:2238`, `:2286`),
    /// a point select is projected onto the curve (`cxx:2249-2251`,
    /// `ShapeAnalysis_Curve::Project`; the port uses its Extrema-based
    /// projection, `int_tools_vertex_line::project_point_on_curve_param`).
    fn extract_parameter(
        &self,
        curve: &dyn Curve,
        params: &[f64],
        points: &[usize],
        master_rep: i32,
        fact: f64,
        shift: f64,
    ) -> Result<Option<f64>, String> {
        if master_rep == 2 {
            if let Some(p) = params.first() {
                return Ok(Some(shift + fact * p));
            }
        } else if master_rep == 1 {
            if let Some(&r) = points.first() {
                let p = self.resolve_point(r)?;
                return Ok(crate::int_tools_vertex_line::project_point_on_curve_param(curve, &p));
            }
        }
        // `cxx:2278-2314`: with an unspecified master representation the parameter
        // is preferred, and the point is projected only when no parameter exists.
        if let Some(p) = params.first() {
            return Ok(Some(shift + fact * p));
        }
        if let Some(&r) = points.first() {
            let p = self.resolve_point(r)?;
            return Ok(crate::int_tools_vertex_line::project_point_on_curve_param(curve, &p));
        }
        Ok(None)
    }

    /// `StepGeom_Axis2Placement3d::HasRefDirection` (`cxx:2404`): the STEP record
    /// carries `$` in the reference-direction slot when it is absent.
    pub(super) fn axis2_has_ref_direction(&self, id: usize) -> Result<bool, String> {
        let rec = self.record(id)?;
        if rec.type_name != "AXIS2_PLACEMENT_3D" {
            return Ok(true);
        }
        Ok(match rec.args.get(3).map(|s| s.trim()) {
            None | Some("") | Some("$") | Some(".F.") => false,
            _ => true,
        })
    }

    pub(super) fn resolve_curve(&self, id: usize) -> Result<Arc<dyn Curve>, String> {
        if let Some(c) = self.curve_cache.borrow().get(&id) {
            return Ok(c.clone());
        }
        let rec = self.record(id)?;
        let curve: Arc<dyn Curve> = match rec.type_name.as_str() {
            "LINE" => {
                let pnt = parse_ref(&rec.args[1]).ok_or("LINE: bad point ref")?;
                let vec = parse_ref(&rec.args[2]).ok_or("LINE: bad vector ref")?;
                let p = self.resolve_point(pnt)?;
                let v = self.resolve_vector(vec)?;
                let d = GpDir::from_vec(&v).map_err(|e| format!("LINE: {e}"))?;
                Arc::new(GeomLine::new(GpLin::from_pnt_dir(p, d)))
            }
            "CIRCLE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("CIRCLE: bad axis ref")?;
                // `StepToGeom::MakeCircle` (`StepToGeom.cxx:1212-1225`):
                // `SC->Radius() * LengthFactor()` (`cxx:1222`).
                let r = parse_f64(&rec.args[2])? * self.length_factor;
                Arc::new(GeomCircle::new(GpCirc::new(self.resolve_axis2(ax)?, r)))
            }
            "ELLIPSE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("ELLIPSE: bad axis ref")?;
                // `StepToGeom::MakeEllipse` (`StepToGeom.cxx:1536-1566`):
                // `majorR = SemiAxis1 * LF`, `minorR = SemiAxis2 * LF`
                // (`cxx:1549-1550`); when `majorR - minorR < 0` OCCT turns the X
                // direction by `A.XDirection() ^ A.Direction()` and swaps the radii
                // (`cxx:1557-1561`).
                let maj = parse_f64(&rec.args[2])? * self.length_factor;
                let min = parse_f64(&rec.args[3])? * self.length_factor;
                let mut ax2 = self.resolve_axis2(ax)?;
                let (maj, min) = if maj - min >= 0.0 {
                    (maj, min)
                } else {
                    if let Ok(xd) = ax2.x_direction().crossed(&ax2.direction()) {
                        ax2.set_x_direction(xd);
                    }
                    (min, maj)
                };
                Arc::new(GeomEllipse::new(GpElips::new(ax2, maj, min)))
            }
            "HYPERBOLA" => {
                let ax = parse_ref(&rec.args[1]).ok_or("HYPERBOLA: bad axis ref")?;
                // `StepToGeom::MakeHyperbola` (`StepToGeom.cxx:1604-1622`):
                // `SemiAxis() * LF`, `SemiImagAxis() * LF` (`cxx:1616-1617`).
                let maj = parse_f64(&rec.args[2])? * self.length_factor;
                let min = parse_f64(&rec.args[3])? * self.length_factor;
                Arc::new(GeomHyperbola::new(GpHypr::new(
                    self.resolve_axis2(ax)?,
                    maj,
                    min,
                )))
            }
            "PARABOLA" => {
                let ax = parse_ref(&rec.args[1]).ok_or("PARABOLA: bad axis ref")?;
                // `StepToGeom::MakeParabola` (`StepToGeom.cxx:1695-1710`):
                // `SC->FocalDist() * LengthFactor()` (`cxx:1706`).
                let f = parse_f64(&rec.args[2])? * self.length_factor;
                Arc::new(GeomParabola::new(GpParab::new(self.resolve_axis2(ax)?, f)))
            }
            "SURFACE_CURVE" | "SEAM_CURVE" => {
                // SURFACE_CURVE/SEAM_CURVE(name, curve_3d, pcurves, master_rep):
                // the 3D curve is the second argument; the pcurve list
                // (associated_geometry) is carried for the face-level pcurve
                // association. SEAM_CURVE is the seam of a closed surface.
                let c3d = parse_ref(&rec.args[1]).ok_or("SURFACE_CURVE: bad 3D curve ref")?;
                self.surface_curve_pcurves
                    .borrow_mut()
                    .insert(id, parse_ref_list(&rec.args[2]));
                self.resolve_curve(c3d)?
            }
            "BEZIER_CURVE" | "UNIFORM_CURVE" | "QUASI_UNIFORM_CURVE" => {
                // `StepToGeom::MakeBSplineCurve` (`StepToGeom.cxx:295-450`): the STEP
                // Bezier / uniform / quasi-uniform curve is converted into a
                // `BSplineCurveWithKnots` before being mapped:
                //   Bezier       -> knots {0,1}, multiplicities degree+1 (`:310-317`)
                //   Uniform      -> n_poles + degree + 1 knots i-1, mults 1 (`:338-347`)
                //   QuasiUniform -> n_poles - degree + 1 knots i-1, mults 1 except
                //                   both ends = degree+1 (`:362-384`)
                // The `_AND_RATIONAL_B_SPLINE_CURVE` complex forms add the weights
                // (`:385-421`, `:422-458`); a rational **Bezier** curve is the
                // exception — OCCT has no Bezier+rational arm, and its plain
                // `BezierCurve` arm (`:295-320`) never reads weights.
                let degree = parse_f64(&rec.args[1])? as usize;
                let poles: Vec<GpPnt> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point(r))
                    .collect::<Result<Vec<_>, _>>()?;
                // `MakeBSplineCurveCommon` synthesises the knot / multiplicity
                // lists this family implies (`StepToGeom.cxx:310-317`,
                // `:338-347`, `:362-384`) and then runs the shared periodic /
                // closed-curve logic.
                let n = poles.len();
                let (mults, knots): (Vec<usize>, Vec<f64>) = match rec.type_name.as_str() {
                    "BEZIER_CURVE" => {
                        (vec![degree + 1, degree + 1], vec![0.0, 1.0])
                    }
                    "UNIFORM_CURVE" => {
                        let nb = n + degree + 1;
                        (vec![1; nb], (0..nb).map(|i| i as f64).collect())
                    }
                    _ => {
                        let nb = n.saturating_sub(degree) + 1;
                        let mut m = vec![1usize; nb];
                        m[0] = degree + 1;
                        m[nb - 1] = degree + 1;
                        (m, (0..nb).map(|i| i as f64).collect())
                    }
                };
                // A merged rational complex carries the weights at index 3; a plain
                // entity has `curve_form` there (a marker like `.UNSPECIFIED.`).
                let weights_arg = rec.args.get(3).map(|s| s.trim().to_string());
                let weights = match weights_arg.as_deref() {
                    Some(w) if w.starts_with('(') && rec.type_name != "BEZIER_CURVE" => {
                        Some(parse_real_list(w))
                    }
                    _ => None,
                };
                if let Some(w) = weights.as_ref() {
                    if w.len() != poles.len() {
                        return Err(format!("{}: weight count mismatch", rec.type_name));
                    }
                }
                let curve = make_bspline_curve_with_knots(
                    degree,
                    poles,
                    weights,
                    &mults,
                    &knots,
                    curve_record_closed(rec),
                )?;
                Arc::new(curve)
            }
            "CURVE_REPLICA" => {
                // `StepToGeom::MakeCurve` CurveReplica arm (`StepToGeom.cxx:1351-1371`):
                // `C1 = MakeCurve(ParentCurve)`, then `C1->Transform(T1)` with
                // `T1 = MakeTransformation3d(Transformation)`. The guard
                // `!T.IsNull() && PC != SC` (`cxx:1358`) rejects a cyclic replica.
                let parent_ref = parse_ref(&rec.args[1]).ok_or("CURVE_REPLICA: bad parent ref")?;
                let trsf_ref =
                    parse_ref(&rec.args[2]).ok_or("CURVE_REPLICA: bad transformation ref")?;
                if parent_ref == id {
                    return Err(format!("CURVE_REPLICA: cyclic parent (#{id})"));
                }
                let parent = self.resolve_curve(parent_ref)?;
                let t = self.make_transformation3d(trsf_ref)?;
                Arc::from(parent.transformed(&t))
            }
            "B_SPLINE_CURVE_WITH_KNOTS" => {
                // Layout (10 args): name, degree, control_points, weights|SELF,
                // curve_form, closed, self_intersect, knots, multiplicities, knot_spec.
                let degree = parse_f64(&rec.args[1])? as usize;
                let poles: Vec<GpPnt> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point(r))
                    .collect::<Result<Vec<_>, _>>()?;
                let weights_arg = rec.args.get(3).map(|s| s.trim().to_string());
                // B_SPLINE_CURVE_WITH_KNOTS layout:
                // (name, degree, control_points, curve_form|weights, curve_form,
                //  closed, self_intersect, knot_multiplicities, knots, knot_spec).
                // `MakeBSplineCurveCommon` (`StepToGeom.cxx:776-927`) merges
                // duplicate knots, clamps the multiplicities, trims the poles and
                // then decides the periodic / closed-curve representation.
                let knot_mults = parse_usize_list(rec.args.get(6).map(|s| s.as_str()).unwrap_or("()"));
                let knot_values = parse_real_list(rec.args.get(7).map(|s| s.as_str()).unwrap_or("()"));
                let weights = match weights_arg.as_deref() {
                    // No weights: non-rational curve (curve_form is UNSPECIFIED
                    // or a non-rational flag like CIRCULAR/LINEAR).
                    None | Some("SELF") | Some(".UNSPECIFIED.") | Some(".CIRCULAR.")
                    | Some(".LINEAR.") => None,
                    Some(w) if w.starts_with('.') => None,
                    Some(w) => {
                        let weights = parse_real_list(w);
                        if weights.len() != poles.len() {
                            return Err("B_SPLINE_CURVE: weight count mismatch".into());
                        }
                        Some(weights)
                    }
                };
                let curve = make_bspline_curve_with_knots(
                    degree,
                    poles,
                    weights,
                    &knot_mults,
                    &knot_values,
                    curve_record_closed(rec),
                )?;
                Arc::new(curve)
            }
            "POLYLINE" => {
                // A connected sequence of CARTESIAN_POINTs; represent it as a
                // degree-1 B-spline that passes through every vertex.
                let pts: Vec<GpPnt> = parse_ref_list(&rec.args[1])
                    .into_iter()
                    .map(|r| self.resolve_point(r))
                    .collect::<Result<Vec<_>, _>>()?;
                if pts.len() < 2 {
                    return Err("POLYLINE: need at least 2 points".into());
                }
                let knots = uniform_knots_for(pts.len(), 1);
                Arc::new(
                    GeomBSplineCurve::new(pts, knots, 1)
                        .map_err(|e| format!("POLYLINE: {e}"))?,
                )
            }
            "B_SPLINE_CURVE" => {
                // Plain B-spline (no explicit knots): the knot vector is the
                // clamped uniform one implied by the pole count and degree.
                let degree = parse_f64(&rec.args[1])? as usize;
                let poles: Vec<GpPnt> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point(r))
                    .collect::<Result<Vec<_>, _>>()?;
                let knots = uniform_knots_for(poles.len(), degree);
                Arc::new(
                    GeomBSplineCurve::new(poles, knots, degree)
                        .map_err(|e| format!("B_SPLINE_CURVE: {e}"))?,
                )
            }
            "TRIMMED_CURVE" => {
                // `StepToGeom::MakeTrimmedCurve` (`StepToGeom.cxx:2323-2496`) with
                // `ExtractParameter` (`:2221-2317`).
                // Layout: name, basis_curve, trim_1, trim_2, sense_agreement,
                // master_representation.
                let basis_ref = parse_ref(&rec.args[1]).ok_or("TRIMMED_CURVE: bad basis ref")?;
                let basis_rec = self.record(basis_ref)?;
                let basis_kind = basis_rec.type_name.clone();
                let basis = self.resolve_curve(basis_ref)?;
                let (params1, points1) = parse_trimming_select(&rec.args[2]);
                let (params2, points2) = parse_trimming_select(&rec.args[3]);
                let sense = !matches!(
                    rec.args.get(4).map(|s| s.trim()),
                    Some(".F.") | Some("F") | Some("false")
                );
                // `cxx:2339-2350`: `.CARTESIAN.` = 1, `.PARAMETER.` = 2, else 0.
                let master_rep: i32 = match rec.args.get(5).map(|s| s.trim()) {
                    Some(".CARTESIAN.") => 1,
                    Some(".PARAMETER.") => 2,
                    _ => 0,
                };
                // `cxx:2352-2373`: with an unspecified master representation (or a
                // parameter one carrying two selects on each side), a Cartesian
                // point on both trims means the trims are points, not parameters.
                let n1 = params1.len() + points1.len();
                let n2 = params2.len() + points2.len();
                let is_point = (master_rep == 0 || (master_rep == 2 && n1 > 1 && n2 > 1))
                    && !points1.is_empty()
                    && !points2.is_empty();
                // `cxx:2375-2392`: parameter scaling per basis type.
                let (fact, shift) = match basis_kind.as_str() {
                    "LINE" => {
                        let vec_ref =
                            parse_ref(&basis_rec.args[2]).ok_or("LINE: bad vector ref")?;
                        // `cxx:2380`: `Dir()->Magnitude() * LengthFactor` — the port's
                        // `resolve_vector` already applies `LengthFactor`
                        // (`StepToGeom::MakeVectorWithMagnitude`, `cxx:2577`).
                        (self.resolve_vector(vec_ref)?.xyz().modulus(), 0.0)
                    }
                    "CIRCLE" | "ELLIPSE" => {
                        // `cxx:2386`: `PlaneAngleFactor`; `cxx:2387-2392`: a π/2 shift
                        // for an ellipse whose `SemiAxis1 - SemiAxis2 < 0`.
                        let mut shift = 0.0;
                        if basis_kind == "ELLIPSE" {
                            let s1 = parse_f64(&basis_rec.args[2])?;
                            let s2 = parse_f64(&basis_rec.args[3])?;
                            if s1 - s2 < 0.0 {
                                shift = 0.5 * PI;
                            }
                        }
                        (self.plane_angle_factor, shift)
                    }
                    _ => (1.0, 0.0),
                };
                // `cxx:2394-2425`: a conic whose placement has no reference direction
                // cannot be trimmed by parameters; OCCT returns the full period with
                // the sense of `SenseAgreement` (the port's `GeomTrimmedCurve` carries
                // no sense flag, so the full-period span is what is modelled).
                if matches!(basis_kind.as_str(), "CIRCLE" | "ELLIPSE") && !is_point && master_rep != 1 {
                    let ax_ref = parse_ref(&basis_rec.args[1]).ok_or("conic: bad axis ref")?;
                    if !self.axis2_has_ref_direction(ax_ref)? {
                        return Ok(Arc::new(GeomTrimmedCurve::new(basis, 0.0, 2.0 * PI)));
                    }
                }
                let t1 = self.extract_parameter(&*basis, &params1, &points1, master_rep, fact, shift)?;
                let t2 = self.extract_parameter(&*basis, &params2, &points2, master_rep, fact, shift)?;
                let (Some(mut trim1), Some(mut trim2)) = (t1, t2) else {
                    // `cxx:2495`: no parameter on either side → OCCT returns null.
                    return Err(format!("TRIMMED_CURVE: no trimming parameters (#{id})"));
                };
                let cf = basis.first_parameter();
                let cl = basis.last_parameter();
                // `cxx:2438-2459`: clamp into the basis range when it is not periodic.
                if !basis.is_periodic() {
                    trim1 = trim1.clamp(cf, cl);
                    trim2 = trim2.clamp(cf, cl);
                }
                if (trim1 - trim2).abs() < occt_core::precision::PCONFUSION {
                    if basis.is_periodic() {
                        // `cxx:2462-2465`: `ElCLib::AdjustPeriodic(cf, cl, PConfusion, …)`.
                        occt_core::elib::clib2d::adjust_periodic(
                            cf,
                            cl,
                            occt_core::precision::PCONFUSION,
                            &mut trim1,
                            &mut trim2,
                        );
                    } else {
                        // `cxx:2466-2480`: the closed (non-periodic) basis arm is
                        // UNPORTED — the `Curve` trait exposes no `IsClosed`, so the
                        // null-result branch OCCT takes otherwise is returned here.
                        return Err(format!("TRIMMED_CURVE: degenerate trim on closed basis (#{id})"));
                    }
                }
                // `cxx:2486-2493`: `SenseAgreement` selects the argument order. The
                // port's `GeomTrimmedCurve` stores only the span (no `Sense` flags,
                // `trimmed.rs`), so an agreed trim is written as-is and an opposite one
                // is written reversed — the reversal flags themselves are UNPORTED.
                if sense {
                    Arc::new(GeomTrimmedCurve::new(basis, trim1, trim2))
                } else {
                    Arc::new(GeomTrimmedCurve::new(basis, trim2, trim1))
                }
            }
            "OFFSET_CURVE_3D" => {
                // Layout: name, basis_curve, direction, distance, self_intersect,
                // curve_form.
                let basis_ref = parse_ref(&rec.args[1]).ok_or("OFFSET_CURVE_3D: bad curve ref")?;
                let dir_ref = parse_ref(&rec.args[2]).ok_or("OFFSET_CURVE_3D: bad dir ref")?;
                let offset = parse_f64(&rec.args[3])?;
                let basis = self.resolve_curve(basis_ref)?;
                let dir = self.resolve_direction(dir_ref)?;
                Arc::new(GeomOffsetCurve::new(basis, offset, dir))
            }
            other => {
                self.warn(format!("unsupported curve entity {other} (#{id})"));
                return Err(format!("unsupported curve entity {other} (#{id})"));
            }
        };
        self.curve_cache.borrow_mut().insert(id, curve.clone());
        Ok(curve)
    }

    /// 2D `CARTESIAN_POINT` → `GpPnt2d` (STEP pcurves use two-component points).
    pub(super) fn resolve_point_2d(&self, id: usize) -> Result<GpPnt2d, String> {
        let rec = self.record(id)?;
        if rec.type_name != "CARTESIAN_POINT" {
            return Err(format!("expected CARTESIAN_POINT at #{id}"));
        }
        let (x, y) = parse_xy(&rec.args[1]).map_err(|e| format!("CARTESIAN_POINT #{id}: {e}"))?;
        Ok(GpPnt2d::new(x, y))
    }

    /// 2D `DIRECTION` → `GpDir2d`.
    pub(super) fn resolve_direction_2d(&self, id: usize) -> Result<GpDir2d, String> {
        let rec = self.record(id)?;
        if rec.type_name != "DIRECTION" {
            return Err(format!("expected DIRECTION at #{id}"));
        }
        let (x, y) = parse_xy(&rec.args[1]).map_err(|e| format!("DIRECTION #{id}: {e}"))?;
        GpDir2d::new(x, y).map_err(|e| format!("DIRECTION #{id}: {e}"))
    }

    /// 2D `VECTOR` → `GpVec2d`.
    pub(super) fn resolve_vector_2d(&self, id: usize) -> Result<GpVec2d, String> {
        let rec = self.record(id)?;
        if rec.type_name != "VECTOR" {
            return Err(format!("expected VECTOR at #{id}"));
        }
        let dir_ref = parse_ref(&rec.args[1]).ok_or("VECTOR: bad direction ref")?;
        let mag = parse_f64(&rec.args[2])?;
        let dir = self.resolve_direction_2d(dir_ref)?;
        Ok(GpVec2d::new(dir.x * mag, dir.y * mag))
    }

    /// 2D `AXIS2_PLACEMENT_2D` → `GpAx22d` (the ref_direction is X; Y is its
    /// counter-clockwise normal).
    pub(super) fn resolve_axis22d(&self, id: usize) -> Result<GpAx22d, String> {
        let rec = self.record(id)?;
        if rec.type_name != "AXIS2_PLACEMENT_2D" {
            return Err(format!("expected AXIS2_PLACEMENT_2D at #{id}"));
        }
        let loc = self.resolve_point_2d(parse_ref(&rec.args[1]).ok_or("AXIS2_2D: bad loc ref")?)?;
        let vx = match parse_ref(&rec.args[2]) {
            Some(r) => self.resolve_direction_2d(r)?,
            None => GpDir2d::new(1.0, 0.0).unwrap(),
        };
        let vy = GpDir2d::new(-vx.y, vx.x).map_err(|e| format!("AXIS2_2D: {e}"))?;
        GpAx22d::new(loc, vx, vy).map_err(|e| format!("AXIS2_PLACEMENT_2D: {e}"))
    }

    /// 2D curve → `Arc<dyn Curve2d>` (STEP pcurve geometry). Source:
    /// `StepToGeom::MakeCurve2d`.
    pub(super) fn resolve_curve_2d(&self, id: usize) -> Result<(Arc<dyn Curve2d>, (f64, f64)), String> {
        if let Some(c) = self.curve2d_cache.borrow().get(&id) {
            return Ok(c.clone());
        }
        let rec = self.record(id)?;
        let (curve, range): (Arc<dyn Curve2d>, (f64, f64)) = match rec.type_name.as_str() {
            "LINE" => {
                let pnt = parse_ref(&rec.args[1]).ok_or("LINE: bad point ref")?;
                let vec = parse_ref(&rec.args[2]).ok_or("LINE: bad vector ref")?;
                let p = self.resolve_point_2d(pnt)?;
                let v = self.resolve_vector_2d(vec)?;
                let d = GpDir2d::new(v.x(), v.y()).map_err(|e| format!("LINE: {e}"))?;
                let mag = (v.x() * v.x() + v.y() * v.y()).sqrt();
                (Arc::new(Geom2dLine::new(GpAx2d::new(p, d))), (0.0, mag))
            }
            "CIRCLE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("CIRCLE: bad axis ref")?;
                let r = parse_f64(&rec.args[2])?;
                let ax22 = self.resolve_axis22d(ax)?;
                let c: Arc<dyn Curve2d> = Arc::new(Geom2dCircle::new(GpCirc2d::new(ax22, r)));
                let range = (c.first_parameter(), c.last_parameter());
                (c, range)
            }
            "ELLIPSE" => {
                let ax = parse_ref(&rec.args[1]).ok_or("ELLIPSE: bad axis ref")?;
                // `StepToGeom::MakeEllipse2d` (`StepToGeom.cxx:1571-1598`): no length
                // factor in the 2-D path (`cxx:1583-1584`), and when
                // `majorR - minorR < 0` the X direction is mirrored in Y —
                // `gp_Dir2d(X.X(), -X.Y())` (`cxx:1591-1593`) — with the radii swapped.
                let maj = parse_f64(&rec.args[2])?;
                let min = parse_f64(&rec.args[3])?;
                let mut ax22 = self.resolve_axis22d(ax)?;
                let (maj, min) = if maj - min >= 0.0 {
                    (maj, min)
                } else {
                    let x = ax22.x_direction();
                    if let Ok(mirror) = GpDir2d::new(x.x(), -x.y()) {
                        ax22.set_x_direction(mirror);
                    }
                    (min, maj)
                };
                let c: Arc<dyn Curve2d> =
                    Arc::new(Geom2dEllipse::new(GpElips2d::new(ax22, maj, min)));
                let range = (c.first_parameter(), c.last_parameter());
                (c, range)
            }
            "B_SPLINE_CURVE_WITH_KNOTS" => {
                let degree = parse_f64(&rec.args[1])? as usize;
                let pts: Vec<GpPnt2d> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point_2d(r))
                    .collect::<Result<Vec<_>, _>>()?;
                // `StepToGeom::MakeBSplineCurve2d` (`StepToGeom.cxx:952-963`) runs
                // the same `MakeBSplineCurveCommon`: duplicate knots merged,
                // multiplicities clamped, poles trimmed, periodic descriptors and
                // closed curves made periodic. UNPORTED here: the weights of a 2D
                // rational complex (`:902-906`) are ignored — this port's 2D
                // B-spline is non-rational (`xs`/`ys` only).
                let knot_mults =
                    parse_usize_list(rec.args.get(6).map(|s| s.as_str()).unwrap_or("()"));
                let knot_values =
                    parse_real_list(rec.args.get(7).map(|s| s.as_str()).unwrap_or("()"));
                let (xs, ys): (Vec<f64>, Vec<f64>) =
                    pts.iter().map(|p| (p.x(), p.y())).unzip();
                let c: Arc<dyn Curve2d> = Arc::new(make_bspline_curve_2d_with_knots(
                    degree,
                    xs,
                    ys,
                    &knot_mults,
                    &knot_values,
                    curve_record_closed(rec),
                )?);
                let range = (c.first_parameter(), c.last_parameter());
                (c, range)
            }
            "POLYLINE" => {
                let pts: Vec<GpPnt2d> = parse_ref_list(&rec.args[1])
                    .into_iter()
                    .map(|r| self.resolve_point_2d(r))
                    .collect::<Result<Vec<_>, _>>()?;
                if pts.len() < 2 {
                    return Err("POLYLINE: need at least 2 points".into());
                }
                let knots = uniform_knots_for(pts.len(), 1);
                let (xs, ys): (Vec<f64>, Vec<f64>) =
                    pts.iter().map(|p| (p.x(), p.y())).unzip();
                let c: Arc<dyn Curve2d> = Arc::new(
                    Geom2dBSplineCurve::new(xs, ys, knots, 1)
                        .map_err(|e| format!("POLYLINE: {e}"))?,
                );
                let range = (c.first_parameter(), c.last_parameter());
                (c, range)
            }
            "B_SPLINE_CURVE" => {
                let degree = parse_f64(&rec.args[1])? as usize;
                let pts: Vec<GpPnt2d> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point_2d(r))
                    .collect::<Result<Vec<_>, _>>()?;
                let knots = uniform_knots_for(pts.len(), degree);
                let (xs, ys): (Vec<f64>, Vec<f64>) =
                    pts.iter().map(|p| (p.x(), p.y())).unzip();
                let c: Arc<dyn Curve2d> = Arc::new(
                    Geom2dBSplineCurve::new(xs, ys, knots, degree)
                        .map_err(|e| format!("B_SPLINE_CURVE_2D: {e}"))?,
                );
                let range = (c.first_parameter(), c.last_parameter());
                (c, range)
            }
            "TRIMMED_CURVE" => {
                // `StepToGeom::MakeTrimmedCurve2d` (`StepToGeom.cxx:2505-2563`).
                let basis_ref = parse_ref(&rec.args[1]).ok_or("TRIMMED_CURVE: bad basis ref")?;
                let basis_rec = self.record(basis_ref)?;
                let basis_kind = basis_rec.type_name.clone();
                let (basis, basis_range) = self.resolve_curve_2d(basis_ref)?;
                // `cxx:2514-2517`: a basis that already is a 2-D B-spline curve is
                // returned **untrimmed**.
                if matches!(
                    basis_kind.as_str(),
                    "B_SPLINE_CURVE" | "B_SPLINE_CURVE_WITH_KNOTS"
                ) {
                    return Ok((basis, basis_range));
                }
                let (params1, _) = parse_trimming_select(&rec.args[2]);
                let (params2, _) = parse_trimming_select(&rec.args[3]);
                // `cxx:2523-2525`: both trims must be single parameter selects,
                // otherwise OCCT returns a null handle.
                if params1.len() != 1 || params2.len() != 1 {
                    return Err(format!("TRIMMED_CURVE(2d): trims are not single parameters (#{id})"));
                }
                let (u1, u2) = (params1[0], params2[0]);
                // `cxx:2528-2551`: Line → `Dir()->Magnitude()` (**no** LengthFactor in
                // the 2-D path); Circle/Ellipse → `PlaneAngleFactor`, plus a π/2 shift
                // for an ellipse with `SemiAxis1 - SemiAxis2 < 0`; parabola/hyperbola
                // are left as a TODO in OCCT itself (`cxx:2547-2551`).
                let (fact, shift) = match basis_kind.as_str() {
                    "LINE" => {
                        let vec_ref = parse_ref(&basis_rec.args[2]).ok_or("LINE: bad vector ref")?;
                        let vec_rec = self.record(vec_ref)?;
                        (parse_f64(&vec_rec.args[2])?, 0.0)
                    }
                    "CIRCLE" | "ELLIPSE" => {
                        let mut shift = 0.0;
                        if basis_kind == "ELLIPSE" {
                            let s1 = parse_f64(&basis_rec.args[2])?;
                            let s2 = parse_f64(&basis_rec.args[3])?;
                            if s1 - s2 < 0.0 {
                                shift = 0.5 * PI;
                            }
                        }
                        (self.plane_angle_factor, shift)
                    }
                    _ => (1.0, 0.0),
                };
                let new_u1 = shift + u1 * fact;
                let new_u2 = shift + u2 * fact;
                // UNPORTED: `Geom2dConvert::CurveToBSplineCurve(theTrimmed)`
                // (`cxx:2560`) is not ported (cf. A8/T-44 for the 3-D counterpart), so
                // the trimmed 2-D curve is kept as the trimmed curve itself; the
                // `SenseAgreement` argument of `Geom2d_TrimmedCurve` is likewise not
                // modelled by the port's `Geom2dTrimmedCurve` (span only).
                (
                    Arc::new(Geom2dTrimmedCurve::new(basis, new_u1, new_u2)),
                    (new_u1, new_u2),
                )
            }
            other => {
                self.warn(format!("unsupported 2D curve entity {other} (#{id})"));
                return Err(format!("unsupported 2D curve entity {other} (#{id})"));
            }
        };
        self.curve2d_cache.borrow_mut().insert(id, (curve.clone(), range));
        Ok((curve, range))
    }

    /// `PCURVE(name, basis_surface, reference_to_curve)` → the pcurve's basis
    /// surface reference and its 2D curve (the first item of the
    /// DEFINITIONAL_REPRESENTATION). Source: `StepToTopoDS_TranslateEdge::MakePCurve`.
    pub(super) fn resolve_pcurve(&self, id: usize) -> Result<(usize, Arc<dyn Curve2d>), String> {
        let rec = self.record(id)?;
        if rec.type_name != "PCURVE" {
            return Err(format!("expected PCURVE at #{id}"));
        }
        let basis_surf = parse_ref(&rec.args[1]).ok_or("PCURVE: bad basis surface ref")?;
        let dri = parse_ref(&rec.args[2]).ok_or("PCURVE: bad reference_to_curve ref")?;
        let dri_rec = self.record(dri)?;
        if dri_rec.type_name != "DEFINITIONAL_REPRESENTATION" {
            return Err(format!("expected DEFINITIONAL_REPRESENTATION at #{dri}"));
        }
        let items = parse_ref_list(&dri_rec.args[1]);
        let curve_ref = *items.first().ok_or("PCURVE: empty reference_to_curve")?;
        let (c2d, _range) = self.resolve_curve_2d(curve_ref)?;
        // `StepToTopoDS_TranslateEdge::MakePCurve` then
        // `GeomConvert_Units::DegreeToRadian` (`GeomConvert_Units.cxx:173-302`).
        let surf = self.resolve_surface(basis_surf)?;
        let c2d = self.degree_to_radian(c2d, curve_ref, surf.as_ref());
        Ok((basis_surf, c2d))
    }

    /// `GeomConvert_Units::DegreeToRadian`. STEP stores cylinder/cone/sphere/torus
    /// pcurve `U` in the file angle unit and cone `V` as axis length; OCCT's
    /// `ElSLib` cone uses `V` along the generatrix (`Length / cos(semiAngle)`).
    /// The length factor is the file's `StepData_Factors::LengthFactor()`
    /// (`StepToTopoDS_TranslateEdge.cxx:573-580` passes it to
    /// `GeomConvert_Units::DegreeToRadian`); the angle factor is
    /// `theLocalFactors.FactorDegreeRadian()`.
    ///
    /// Line location/direction are scaled then the direction is renormalized
    /// (`GeomConvert_Units.cxx:260-276`), so a V-iso generatrix keeps parameter
    /// speed 1 while a U-iso circle at written `V = height` moves to ElSLib `V`.
    /// BSpline poles are scaled in place (`cxx:278-294`). Plane Circle/Ellipse
    /// use `SetScale(origin, LengthFact)` (`cxx:215-222`). Circle/Ellipse on a
    /// revolution surface go through `Geom2dConvert::CurveToBSplineCurve` first
    /// (`cxx:229-236`) — that convert is not in this crate; those entities keep
    /// the raw 2D curve. Parabola/Hyperbola return unchanged (`cxx:237-252`).
    pub(super) fn degree_to_radian(
        &self,
        c2d: Arc<dyn Curve2d>,
        curve_ref: usize,
        surf: &dyn Surface,
    ) -> Arc<dyn Curve2d> {
        let length_fact = self.length_factor;
        let angle_fact = self.plane_angle_factor;
        let is_plane = classify_surface(surf) == SurfaceKind::Plane;
        // `GeomConvert_Units.cxx:191-227`. Offset (and other non-analytic
        // kinds) return the pcurve unchanged. Revolution is U-angle only.
        let (u_fact, v_fact) = if let Some((_, alpha)) = surf.cone_ref() {
            (angle_fact, length_fact / alpha.cos())
        } else if surf.is_surface_of_revolution() {
            (angle_fact, length_fact)
        } else {
            match classify_surface(surf) {
                SurfaceKind::Sphere | SurfaceKind::Torus => (angle_fact, angle_fact),
                SurfaceKind::Cylinder => (angle_fact, length_fact),
                SurfaceKind::Plane => (length_fact, length_fact),
                SurfaceKind::Cone | SurfaceKind::Other => return c2d,
            }
        };
        let Ok(rec) = self.record(curve_ref) else {
            return c2d;
        };
        match rec.type_name.as_str() {
            "CIRCLE" | "ELLIPSE" if is_plane => {
                let mut t = occt_core::gp::GpTrsf2d::default();
                if t.set_scale(&GpPnt2d::new(0.0, 0.0), length_fact).is_err() {
                    return c2d;
                }
                let mut scaled = c2d.clone_dyn();
                scaled.transform(&t);
                Arc::from(scaled)
            }
            "PARABOLA" | "HYPERBOLA" => c2d,
            _ if (u_fact - 1.0).abs() <= 1e-15 && (v_fact - 1.0).abs() <= 1e-15 => c2d,
            "LINE" => {
                let loc = c2d.d0(0.0);
                let new_loc = GpPnt2d::new(loc.x() * u_fact, loc.y() * v_fact);
                let (_, tan) = c2d.d1(0.0);
                let Ok(new_dir) = GpDir2d::new(tan.x() * u_fact, tan.y() * v_fact) else {
                    return c2d;
                };
                Arc::new(Geom2dLine::from_pnt_dir(new_loc, new_dir))
            }
            "B_SPLINE_CURVE_WITH_KNOTS" | "B_SPLINE_CURVE" | "POLYLINE" => self
                .scale_pcurve_bspline_poles(curve_ref, u_fact, v_fact)
                .unwrap_or(c2d),
            _ => c2d,
        }
    }

    /// `GeomConvert_Units.cxx:278-294` — affinity on BSpline poles.
    fn scale_pcurve_bspline_poles(
        &self,
        curve_ref: usize,
        u_fact: f64,
        v_fact: f64,
    ) -> Option<Arc<dyn Curve2d>> {
        let rec = self.record(curve_ref).ok()?;
        let (pts, knots, degree): (Vec<GpPnt2d>, Vec<f64>, usize) = match rec.type_name.as_str() {
            "B_SPLINE_CURVE_WITH_KNOTS" => {
                let degree = parse_f64(&rec.args[1]).ok()? as usize;
                let pts: Vec<GpPnt2d> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point_2d(r))
                    .collect::<Result<Vec<_>, _>>()
                    .ok()?;
                let knots = expand_knots(
                    &parse_usize_list(rec.args.get(6).map(|s| s.as_str()).unwrap_or("()")),
                    &parse_real_list(rec.args.get(7).map(|s| s.as_str()).unwrap_or("()")),
                );
                (pts, knots, degree)
            }
            "B_SPLINE_CURVE" => {
                let degree = parse_f64(&rec.args[1]).ok()? as usize;
                let pts: Vec<GpPnt2d> = parse_ref_list(&rec.args[2])
                    .into_iter()
                    .map(|r| self.resolve_point_2d(r))
                    .collect::<Result<Vec<_>, _>>()
                    .ok()?;
                let knots = uniform_knots_for(pts.len(), degree);
                (pts, knots, degree)
            }
            "POLYLINE" => {
                let pts: Vec<GpPnt2d> = parse_ref_list(&rec.args[1])
                    .into_iter()
                    .map(|r| self.resolve_point_2d(r))
                    .collect::<Result<Vec<_>, _>>()
                    .ok()?;
                if pts.len() < 2 {
                    return None;
                }
                let knots = uniform_knots_for(pts.len(), 1);
                (pts, knots, 1)
            }
            _ => return None,
        };
        let xs: Vec<f64> = pts.iter().map(|p| p.x() * u_fact).collect();
        let ys: Vec<f64> = pts.iter().map(|p| p.y() * v_fact).collect();
        Geom2dBSplineCurve::new(xs, ys, knots, degree)
            .ok()
            .map(|c| Arc::new(c) as Arc<dyn Curve2d>)
    }
}

