use super::prelude::*;


/// State of a 2D point relative to the face region. Source: `TopAbs_State`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]

pub enum FaceState {
    In,
    Out,
    On,
    Unknown,
}

/// A face boundary region decomposed into an outer ring and hole rings.
#[derive(Debug, Clone, Default)]
pub struct FaceRegion {
    pub outer: Vec<GpPnt2d>,
    pub holes: Vec<Vec<GpPnt2d>>,
}

impl FaceRegion {
    /// Classify a 2D point against the region.
    ///
    /// `On` when the point is within `on_tol` of a boundary ring segment; `In`
    /// when it is inside the outer ring and outside every hole; `Out` otherwise.
    /// An empty region (no boundary) classifies every point `In`.
    pub(super) fn point_in_region(&self, p: &GpPnt2d, on_tol: f64) -> FaceState {
        if self.outer.is_empty() {
            return FaceState::In;
        }
        if point_near_ring(&self.outer, p, on_tol) {
            return FaceState::On;
        }
        for hole in &self.holes {
            if point_near_ring(hole, p, on_tol) {
                return FaceState::On;
            }
        }
        if !point_in_polygon2d(&self.outer, p) {
            return FaceState::Out;
        }
        for hole in &self.holes {
            if point_in_polygon2d(hole, p) {
                return FaceState::Out;
            }
        }
        FaceState::In
    }
}

/// 2D point-in-face classifier. Source: `IntTools_FClass2d`.
#[derive(Debug, Clone)]
pub struct FClass2d {
    /// UV tolerance used for the `On` boundary test.
    pub(super) tol: f64,
    /// The face whose boundary was sampled.
    pub(super) face: Face,
    /// The sampled boundary region (None when the face has no boundary wires).
    pub(super) region: Option<FaceRegion>,
    /// Whether the face surface is periodic in `u` / `v`.
    pub(super) is_u_periodic: bool,
    pub(super) is_v_periodic: bool,
    /// Surface periods (0.0 for non-periodic directions).
    pub(super) u_period: f64,
    pub(super) v_period: f64,
    /// Bounding box of the sampled boundary rings (`IntTools_FClass2d::Umin`…).
    pub(super) umin: f64,
    pub(super) umax: f64,
    pub(super) vmin: f64,
    pub(super) vmax: f64,
    /// Whether the face is a "hole" (its largest boundary loop winds clockwise
    /// in the face's UV). Source: `IntTools_FClass2d::IsHole`.
    pub(super) my_is_hole: bool,
    /// Per-wire 2D classifiers (`IntTools_FClass2d::TabClass`).
    #[allow(dead_code)]
    pub(super) tab_class: Vec<Class2d>,
    /// Per-wire orientation (`IntTools_FClass2d::TabOrien`): 1 = area>0,
    /// 0 = area<0, -1 = bad wire.
    #[allow(dead_code)]
    pub(super) tab_orien: Vec<i32>,
    /// CSLib classifier of the winding-independent outer ring.
    pub(super) outer_class: Option<Class2d>,
    /// CSLib classifiers of hole rings.
    pub(super) hole_classes: Vec<Class2d>,
}

impl FClass2d {
    /// Build a classifier for `face` using the UV tolerance `tol`
    /// (`IntTools_FClass2d(F, TolUV)`).
    pub fn new(face: &Face, tol: f64) -> Result<Self, String> {
        let mut c = FClass2d {
            tol,
            face: Face::new(),
            region: None,
            is_u_periodic: false,
            is_v_periodic: false,
            u_period: 0.0,
            v_period: 0.0,
            umin: f64::INFINITY,
            umax: f64::NEG_INFINITY,
            vmin: f64::INFINITY,
            vmax: f64::NEG_INFINITY,
            my_is_hole: true,
            tab_class: Vec::new(),
            tab_orien: Vec::new(),
            outer_class: None,
            hole_classes: Vec::new(),
        };
        c.init(face, tol)?;
        Ok(c)
    }

    /// (Re)initialize the classifier from `face` and tolerance `tol`
    /// (`IntTools_FClass2d::Init`).
    ///
    /// Samples every boundary wire into UV rings. Bad wires (open, degenerate,
    /// or with a zero-area ring) are skipped rather than failing the call, so a
    /// face with several wires still initializes when one is unusable.
    pub fn init(&mut self, face: &Face, tol: f64) -> Result<(), String> {
        self.tol = tol;
        self.face = face.clone();
        let Some(surf) = GeometryRegistry::global().face_surface(&face.0) else {
            return Err("FClass2d::init: face has no registered surface".into());
        };
        self.is_u_periodic = surf.is_u_periodic();
        self.is_v_periodic = surf.is_v_periodic();

        // Surface period from the natural parametric range; 2π fallback for an
        // unbounded periodic direction (cylinder/cone/sphere/torus use 2π).
        let (su0, su1, sv0, sv1) = face_uv_bounds(face);
        self.u_period = if self.is_u_periodic {
            if su0.is_finite() && su1.is_finite() {
                su1 - su0
            } else {
                2.0 * PI
            }
        } else {
            0.0
        };
        self.v_period = if self.is_v_periodic {
            if sv0.is_finite() && sv1.is_finite() {
                sv1 - sv0
            } else {
                2.0 * PI
            }
        } else {
            0.0
        };

        let mut rings: Vec<(Vec<GpPnt2d>, f64)> = Vec::new();
        let mut umin = f64::INFINITY;
        let mut umax = f64::NEG_INFINITY;
        let mut vmin = f64::INFINITY;
        let mut vmax = f64::NEG_INFINITY;

        for wire in wires_of_face(face) {
            // BRepTopAdaptor_FClass2d::Init walks the wire with BRepTools_WireExplorer
            // (cxx:352), so edges are visited in connection order, not storage order.
            let edges = crate::meshing::model_builder::wire_builder::wire_edges_explorer(&wire, face);
            let mut polylines: Vec<Vec<GpPnt2d>> = Vec::new();
            for edge in &edges {
                let or = edge.orientation();
                if or != Orientation::Forward && or != Orientation::Reversed {
                    continue;
                }
                let mut pl = edge_points(edge, face);
                if pl.len() < 2 {
                    continue;
                }
                // BRepAdaptor_Curve2d(edge, face) walks the pcurve in the
                // edge's orientation; reverse the sample so REVERSED edges
                // contribute the opposite UV traversal (IntTools_FClass2d::Init).
                if or == Orientation::Reversed {
                    pl.reverse();
                }
                polylines.push(pl);
            }
            if polylines.is_empty() {
                continue;
            }
            // BRepTopAdaptor_FClass2d::Init (cxx:370-410): samples are appended in
            // wire order, the shared first sample of each later edge is skipped, and
            // the closing chord is implicit. No chaining or closure test applies.
            let ring = concat_wire_polylines(&polylines);
            if ring.len() < 3 {
                continue;
            }
            for p in &ring {
                umin = umin.min(p.x());
                umax = umax.max(p.x());
                vmin = vmin.min(p.y());
                vmax = vmax.max(p.y());
            }
            let area = polygon_area2d(&ring);
            if area.abs() < SQUARE_CONFUSION {
                continue;
            }
            rings.push((ring, area));
        }

        self.umin = umin;
        self.umax = umax;
        self.vmin = vmin;
        self.vmax = vmax;

        // Outer ring: the one with the largest |area|; every other ring is a
        // hole. This is winding-independent, so it stays correct for faces
        // whose UV loops happen to wind clockwise (e.g. inward-normal planes).
        let mut region = FaceRegion::default();
        if let Some((outer_idx, _)) = rings
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.1.abs().total_cmp(&b.1.abs()))
        {
            region.outer = rings[outer_idx].0.clone();
            for (i, (ring, _)) in rings.iter().enumerate() {
                if i != outer_idx {
                    region.holes.push(ring.clone());
                }
            }
        }
        self.region = if region.outer.is_empty() { None } else { Some(region) };

        // TabClass / TabOrien: one CSLib_Class2d per sampled wire, using the
        // face UV box and a Fleche of at least Toluv (OCCT floors FlecheU/V).
        self.tab_class.clear();
        self.tab_orien.clear();
        let fleche_u = self.tol;
        let fleche_v = self.tol;
        for (ring, area) in &rings {
            self.tab_class.push(Class2d::new(
                ring, fleche_u, fleche_v, umin, vmin, umax, vmax,
            ));
            self.tab_orien.push(if *area > 0.0 { 1 } else { 0 });
        }
        self.outer_class = self.region.as_ref().map(|r| {
            Class2d::new(&r.outer, fleche_u, fleche_v, umin, vmin, umax, vmax)
        });
        self.hole_classes = self
            .region
            .as_ref()
            .map(|r| {
                r.holes
                    .iter()
                    .map(|h| Class2d::new(h, fleche_u, fleche_v, umin, vmin, umax, vmax))
                    .collect()
            })
            .unwrap_or_default();

        // A face is a "hole" when its largest boundary loop winds clockwise
        // (material on the outside of the loop). An empty face counts as a hole.
        self.my_is_hole = match rings
            .iter()
            .max_by(|a, b| a.1.abs().total_cmp(&b.1.abs()))
        {
            Some((_, area)) => *area <= 0.0,
            None => true,
        };
        Ok(())
    }

    /// Classify the 2D point `puv` (`IntTools_FClass2d::Perform`).
    ///
    /// Periodic surface coordinates are folded by whole periods before
    /// classifying, and the periodic images are tried in order.
    pub fn perform(&self, puv: GpPnt2d) -> FaceState {
        self.perform_recadre(puv, true)
    }

    /// `Perform` with explicit `RecadreOnPeriodic` control.
    pub fn perform_recadre(&self, puv: GpPnt2d, recadre_on_periodic: bool) -> FaceState {
        self.perform_internal(puv, self.tol, recadre_on_periodic, false)
    }

    /// State of the infinite (far bottom-left) point: `In` when the face has
    /// no boundary (a closed periodic face such as a full sphere/torus),
    /// otherwise the classification of the UV-domain-outside corner point
    /// (`IntTools_FClass2d::PerformInfinitePoint`).
    pub fn perform_infinite_point(&self) -> FaceState {
        if !self.umin.is_finite()
            || !self.umax.is_finite()
            || !self.vmin.is_finite()
            || !self.vmax.is_finite()
        {
            return FaceState::In;
        }
        let p = GpPnt2d::new(
            self.umin - (self.umax - self.umin),
            self.vmin - (self.vmax - self.vmin),
        );
        self.perform_recadre(p, false)
    }

    /// `BRepTopAdaptor_FClass2d::Perform` with `RecadreOnPeriodic = false`
    /// (`BRepTopAdaptor_FClass2d.cxx:595-620`), the path `PerformInfinitePoint`
    /// and `ShapeFix_Face::FixOrientation` use.
    ///
    /// A point is IN only when it is inside every positive wire (`TabOrien == 1`,
    /// signed area above 0, the same case as OCCT `square < 0` at `cxx:445`)
    /// and outside every negative wire (`TabOrien == 0`). A clockwise
    /// loop therefore classifies its exterior as IN. `classify_tab` does not:
    /// it treats the geometric interior as IN regardless of winding, which
    /// would reverse the wrong wires in `FixOrientation`.
    ///
    /// `TabOrien(1) == -1` falls through to `BRepClass_FaceClassifier`
    /// (`cxx:640-642`). That classifier is not this type; an uncertain sample
    /// is reported `On`, which `FixOrientation` ignores (`cxx:1436`).
    pub fn perform_tab_orien(&self, puv: GpPnt2d) -> FaceState {
        if self.tab_class.is_empty() {
            return FaceState::In;
        }
        if self.tab_orien.first().copied().unwrap_or(-1) < 0 {
            return FaceState::On;
        }
        let mut dedans = 1i32;
        for (n, clas) in self.tab_class.iter().enumerate() {
            let orien = self.tab_orien.get(n).copied().unwrap_or(0);
            match clas.si_dans(&puv) {
                Class2dResult::Inside => {
                    if orien == 0 {
                        dedans = -1;
                        break;
                    }
                }
                Class2dResult::Outside => {
                    if orien == 1 {
                        dedans = -1;
                        break;
                    }
                }
                Class2dResult::Uncertain => {
                    dedans = 0;
                    break;
                }
            }
        }
        match dedans {
            1 => FaceState::In,
            0 => FaceState::On,
            _ => FaceState::Out,
        }
    }

    /// `PerformInfinitePoint` (`BRepTopAdaptor_FClass2d.cxx:515-522`) using
    /// [`Self::perform_tab_orien`] rather than the winding-independent classifier.
    pub fn perform_infinite_point_tab_orien(&self) -> FaceState {
        if !self.umin.is_finite()
            || !self.umax.is_finite()
            || !self.vmin.is_finite()
            || !self.vmax.is_finite()
        {
            return FaceState::In;
        }
        let p = GpPnt2d::new(
            self.umin - (self.umax - self.umin),
            self.vmin - (self.vmax - self.vmin),
        );
        self.perform_tab_orien(p)
    }

    /// Test whether `puv` lies on the face boundary restriction within `tol`
    /// (`IntTools_FClass2d::TestOnRestriction`). `On` when within `tol` of a
    /// boundary ring, `In`/`Out` otherwise.
    pub fn test_on_restriction(&self, puv: GpPnt2d, tol: f64) -> FaceState {
        self.perform_internal(puv, tol, true, true)
    }

    /// Whether the face is a "hole" (`IntTools_FClass2d::IsHole`).
    pub fn is_hole(&self) -> bool {
        self.my_is_hole
    }

    /// Classify `puv` against the sampled boundary region using the UV
    /// tolerance `on_tol`.
    pub fn point_in_region(&self, puv: GpPnt2d) -> FaceState {
        match &self.region {
            Some(r) => r.point_in_region(&puv, self.tol),
            None => FaceState::In,
        }
    }

    /// The sampled outer ring, if any (exposed for diagnostics/tests).
    pub fn outer_ring(&self) -> Option<&[GpPnt2d]> {
        self.region.as_ref().map(|r| r.outer.as_slice())
    }

    /// The sampled hole rings, if any.
    pub fn hole_rings(&self) -> &[Vec<GpPnt2d>] {
        self.region.as_ref().map(|r| r.holes.as_slice()).unwrap_or(&[])
    }

    /// `Perform` core: the periodic-image search loop from the C++ source.
    pub(super) fn perform_internal(&self, puv: GpPnt2d, on_tol: f64, recadre: bool, on_mode: bool) -> FaceState {
        if self.tab_class.is_empty() {
            return FaceState::In;
        }
        let mut u = puv.x();
        let mut v = puv.y();
        let (mut uu, mut vv) = (u, v);
        if recadre {
            if self.is_u_periodic && self.u_period > 0.0 {
                uu = adjust_periodic(u, self.umin, self.umax, self.u_period).0;
            }
            if self.is_v_periodic && self.v_period > 0.0 {
                vv = adjust_periodic(v, self.vmin, self.vmax, self.v_period).0;
            }
        }
        let mut urecadre = false;
        let mut vrecadre = false;
        let mut a_status = FaceState::Unknown;
        for _ in 0..64 {
            a_status = self.classify_tab(&GpPnt2d::new(u, v), on_tol, on_mode);

            if !recadre || (!self.is_u_periodic && !self.is_v_periodic) {
                return a_status;
            }
            if a_status == FaceState::In || a_status == FaceState::On {
                return a_status;
            }
            if !urecadre {
                u = uu;
                urecadre = true;
            } else if self.is_u_periodic {
                u += self.u_period;
            }
            if u > self.umax || !self.is_u_periodic {
                if !vrecadre {
                    v = vv;
                    vrecadre = true;
                } else if self.is_v_periodic {
                    v += self.v_period;
                }
                u = uu;
                if v > self.vmax || !self.is_v_periodic {
                    return a_status;
                }
            }
        }
        a_status
    }

    /// Classify with CSLib_Class2d on the winding-independent outer/hole rings.
    /// TabOrien is stored (Init) but not used here: combining by signed area
    /// would flip clockwise single-wire faces (unit-box bottom) to Out.
    pub(super) fn classify_tab(&self, puv: &GpPnt2d, on_tol: f64, on_mode: bool) -> FaceState {
        let Some(outer) = self.outer_class.as_ref() else {
            return FaceState::In;
        };
        let si = |c: &Class2d| {
            if on_mode {
                c.si_dans_on_mode(puv, on_tol)
            } else {
                c.si_dans(puv)
            }
        };
        match si(outer) {
            Class2dResult::Uncertain => FaceState::On,
            Class2dResult::Outside => FaceState::Out,
            Class2dResult::Inside => {
                for hole in &self.hole_classes {
                    match si(hole) {
                        Class2dResult::Uncertain => return FaceState::On,
                        Class2dResult::Inside => return FaceState::Out,
                        Class2dResult::Outside => {}
                    }
                }
                FaceState::In
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Edge sampling and ring chaining
// ---------------------------------------------------------------------------

/// Sample the UV pcurve of `edge` on `face` into a polyline (in the edge's
/// natural curve direction).
///
/// Preference: a pcurve stored on the edge for this face (matching
/// `BRep_Tool::CurveOnSurface`), then [`make_pcurve_full`], then a sampled
/// projection via [`edge_pcurve_on_face`]. Empty when the edge has no usable
/// curve or a degenerate parameter range.
pub(super) fn edge_points(edge: &Edge, face: &Face) -> Vec<GpPnt2d> {
    let (a, b) = GeometryRegistry::global().edge_parameters(&edge.0);
    if !a.is_finite() || !b.is_finite() || b - a < 1e-15 {
        return Vec::new();
    }
    // BRepAdaptor_Curve2d(edge, face) reads BRep_Tool::CurveOnSurface
    // (BRep_Tool.cxx:301-315), which returns PCurve2 for a REVERSED edge of a
    // representation on a closed surface (cxx:347-357). Reading only the first
    // stored pcurve made a reversed seam sample the u=2*pi line backwards
    // instead of the u=0 line, so the band ring could not be chained.
    if let Some(pc) = crate::boptools_2d::curve_on_surface(edge, face) {
        return sample_pcurve(pc.as_ref(), a, b);
    }
    if let Ok(pc) = make_pcurve_full(edge, face) {
        return sample_pcurve(pc.as_ref(), a, b);
    }
    edge_pcurve_on_face(edge, face, 32)
}

/// Sample a pcurve uniformly over the edge range `[a, b]`.
pub(super) fn sample_pcurve(pc: &dyn Curve2d, a: f64, b: f64) -> Vec<GpPnt2d> {
    let n = sample_count(pc, a, b);
    (0..n)
        .map(|i| {
            let t = a + (b - a) * i as f64 / (n.max(1) - 1) as f64;
            pc.d0(t)
        })
        .collect()
}

/// Number of boundary samples for a pcurve, mirroring
/// `BRepTopAdaptor_FClass2d::Perform` (`BRepTopAdaptor_FClass2d.cxx:179-185`):
///
/// ```text
/// Standard_Integer nbs = Geom2dInt_Geom2dCurveTool::NbSamples(aCurveAdaptor2D);
/// if (nbs > 2)
///   nbs *= 4;
/// ```
///
/// `first`/`last` are the adaptor's `FirstParameter()`/`LastParameter()` (the
/// edge pcurve range) and only matter for the circle arm of
/// [`nb_samples`]. `du` is then
/// `(plbid - pfbid) / (nbs - 1)` (`BRepTopAdaptor_FClass2d.cxx:186`), i.e.
/// uniform over `[a, b]`.
pub(super) fn sample_count(pc: &dyn Curve2d, first: f64, last: f64) -> usize {
    // `BRepTopAdaptor_FClass2d.cxx:180` -> `Geom2dInt_Geom2dCurveTool::NbSamples`
    // (`Geom2dInt_Geom2dCurveTool.cxx:73-91`) -> `Geom2dAdaptor_Curve::NbSamples`
    // (`Geom2dAdaptor_Curve.cxx:1391-1394`, body `cxx:1351-1389`).
    let mut nbs = nb_samples(pc, first, last);
    if nbs > 2 {
        // `BRepTopAdaptor_FClass2d.cxx:182-184`.
        nbs *= 4;
    }
    nbs.max(2)
}

/// Concatenate edge polylines in wire order (BRepTopAdaptor_FClass2d::Init,
/// cxx:370-410). Each later polyline drops its first sample, which repeats the
/// shared vertex of the previous edge; the ring is not required to close.
pub(super) fn concat_wire_polylines(polylines: &[Vec<GpPnt2d>]) -> Vec<GpPnt2d> {
    let mut ring: Vec<GpPnt2d> = Vec::new();
    for (k, pl) in polylines.iter().enumerate() {
        let skip = if k == 0 { 0 } else { 1 };
        ring.extend(pl.iter().skip(skip).copied());
    }
    ring
}
/// Fold `u` into `[umin, umax]` by adding/subtracting whole periods
/// (`GeomInt::AdjustPeriodic` semantics). Returns `(folded, offset)`.
pub(super) fn adjust_periodic(u: f64, umin: f64, umax: f64, period: f64) -> (f64, f64) {
    if period <= 0.0
        || !period.is_finite()
        || !umin.is_finite()
        || !umax.is_finite()
        || umax <= umin
    {
        return (u, 0.0);
    }
    if u >= umin && u <= umax {
        return (u, 0.0);
    }
    let du = if u < umin {
        ((umin - u) / period).ceil() * period
    } else {
        -((u - umax) / period).ceil() * period
    };
    (u + du, du)
}

/// Minimum distance from `p` to any segment of a closed ring.
pub(super) fn point_near_ring(ring: &[GpPnt2d], p: &GpPnt2d, tol: f64) -> bool {
    let n = ring.len();
    if n < 2 {
        return false;
    }
    for i in 0..n {
        let a = &ring[i];
        let b = &ring[(i + 1) % n];
        if point_segment_distance(p, a, b) <= tol {
            return true;
        }
    }
    false
}

/// Distance from a point to the segment `[a, b]`.
pub(super) fn point_segment_distance(p: &GpPnt2d, a: &GpPnt2d, b: &GpPnt2d) -> f64 {
    let dx = b.x() - a.x();
    let dy = b.y() - a.y();
    let len2 = dx * dx + dy * dy;
    if len2 < 1e-24 {
        return p.distance(a);
    }
    let t = (((p.x() - a.x()) * dx + (p.y() - a.y()) * dy) / len2).clamp(0.0, 1.0);
    let qx = a.x() + t * dx;
    let qy = a.y() + t * dy;
    ((p.x() - qx).powi(2) + (p.y() - qy).powi(2)).sqrt()
}
