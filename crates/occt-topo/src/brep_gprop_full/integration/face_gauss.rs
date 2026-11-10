use super::*;


impl FaceGauss {
    pub(in crate::brep_gprop_full) fn new(face: &Face) -> Result<Self, String> {
        let surface = BRepTool::face_surface_world(face).ok_or("brep_gprop_full: face has no surface")?;
        let kind = classify_surface_kind(surface.as_ref());
        let is_reversed = face.orientation().is_reversed();
        // `BRepGProp_Domain::Init(const TopoDS_Face& F)` (`BRepGProp_Domain.lxx:
        // 39-42`; ctor `:29-32`) explores `F.Oriented(TopAbs_FORWARD)`: the
        // domain's boundary edges carry only the *wire/edge* orientation, never
        // the face's. The face orientation enters the integral in exactly one
        // place - `BRepGProp_Face::mySReverse` (`BRepGProp_Face.cxx:196`, used
        // by `Normal` `:201-210`).
        //
        // Composing the face orientation into the edge views as well (what a
        // plain `TopExp_Explorer(face, TopAbs_EDGE)` does) makes `Dul =
        // Vuv.Y()` (`compute_domain`, `:524`/`:527`, from the p-curve reversal
        // of `BRepGProp_Face::Load(const TopoDS_Edge&)` `:164-185`) flip *with*
        // the normal, so the two cancel and the face orientation has no effect
        // on the integral at all. The fix therefore belongs here, at the
        // boundary source, and NOT in `compute_domain` or `normal`: every
        // single-factor change there just moves the sign between the two
        // cancelling terms.
        let fwd_face = Face(face.0.oriented(crate::abs::Orientation::Forward));
        let wires = wires_of_face(&fwd_face);
        let natural = wires.is_empty();

        // Build the UV inverse map. The frame is taken from the surface's own
        // parameterization so pcurve UV coordinates live in the same space the
        // `normal`/integrand evaluators use.
        let map = match kind {
            SurfKind::Plane => {
                let (o, xd, yd) =
                    plane_frame(surface.as_ref()).ok_or("brep_gprop_full: planar face frame")?;
                UVMap::Plane { o, xd, yd }
            }
            SurfKind::Cylinder => {
                let (o, xd, yd, z, r) =
                    cylinder_frame(surface.as_ref()).ok_or("brep_gprop_full: bad cylinder frame")?;
                UVMap::Cylinder { o, xd, yd, z, r }
            }
            SurfKind::Cone => {
                let (o, xd, yd, z, alpha) =
                    cone_frame(surface.as_ref()).ok_or("brep_gprop_full: bad cone frame")?;
                UVMap::Cone { o, xd, yd, z, alpha }
            }
            _ => UVMap::Generic,
        };

        // Boundary arcs from every edge of every wire.
        let mut arcs: Vec<BoundaryArc> = Vec::new();
        let mut seen: Vec<*const ()> = Vec::new();
        let mut has_repeated = false;
        for w in &wires {
            for e in edges_of_wire(w) {
                let key = Arc::as_ptr(&e.0.tshape) as *const ();
                if seen.contains(&key) {
                    has_repeated = true;
                }
                seen.push(key);
                // BRepGProp_Domain hands BRepGProp_Face::Load(const
                // TopoDS_Edge&) (BRepGProp_Face.cxx:164-185) an edge whose
                // orientation is the wire/edge composition only - the face was
                // oriented FORWARD by Init (BRepGProp_Domain.lxx:41), and that
                // is the orientation the p-curve reversal keys off
                // (BRep_Tool::CurveOnSurface, BRep_Tool.cxx:334,354-361).
                if let Some(arc) = build_arc(&e, face, &map) {
                    arcs.push(arc);
                }
            }
        }

        // UV bounds from the pcurves (BRepTools::UVBounds).
        let (u1, u2, v1, v2) = if natural {
            let (a, b) = surface.u_range();
            let (c, d) = surface.v_range();
            (a, b, c, d)
        } else {
            uv_bounds(surface.as_ref(), &arcs)
        };

        Ok(Self { surface, kind, is_reversed, natural, u1, u2, v1, v2, arcs, has_repeated_edges: has_repeated })
    }

    /// `BRepGProp_Face::Bounds` (`BRepGProp_Face.cxx:154-160`): the
    /// `BRepAdaptor_Surface` adaptor parameters. `BRepGProp_Face::Load(const
    /// TopoDS_Face&)` initialises that adaptor as `Initialize(F)` with its
    /// default `Restriction = true` (`BRepAdaptor_Surface.hxx:59`,
    /// `BRepGProp_Face.cxx:189-197`), so `Initialize` calls `Load(surface,
    /// umin, umax, vmin, vmax, ...)` with `BRepTools::UVBounds(F)`
    /// (`BRepAdaptor_Surface.cxx:56-81`). `GeomAdaptor_Surface::Load` stores
    /// those as `myUFirst`/`myULast`/`myVFirst`/`myVLast`
    /// (`GeomAdaptor_Surface.hxx:163-181`, `GeomAdaptor_Surface.cxx:404-416`),
    /// so `FirstUParameter()`/`LastUParameter()` report the **face** UV box for
    /// every surface type - there is no plane/other split.
    ///
    /// This matters: `BRepGProp_Gauss::Compute` clamps the boundary U into
    /// `[BU1, BU2]` (`BRepGProp_Gauss.cxx:752-759`). Using the natural periodic
    /// range (0, 2*pi) instead of the face box clips any face whose p-curves
    /// live in the negative branch of the periodic parameter (e.g. a cylinder
    /// face at u in [-0.794, 0.794], T0M face 987) and halves its area.
    pub(in crate::brep_gprop_full) fn surface_ranges(&self) -> (f64, f64, f64, f64) {
        self.bounds()
    }

    pub(in crate::brep_gprop_full) fn bounds(&self) -> (f64, f64, f64, f64) {
        (self.u1, self.u2, self.v1, self.v2)
    }

    /// Unnormalised surface normal (D1U × D1V), flipped for a REVERSED face.
    pub(in crate::brep_gprop_full) fn normal(&self, u: f64, v: f64) -> (GpPnt, GpVec) {
        let (p, du, dv) = surface_d1(self.surface.as_ref(), u, v);
        let mut n = du.crossed(&dv);
        if self.is_reversed {
            n.reverse();
        }
        (p, n)
    }

    /// Number of Gauss points for the U direction of the face.
    pub(in crate::brep_gprop_full) fn u_integration_order(&self) -> usize {
        match self.kind {
            SurfKind::Plane => 8,
            _ => 18,
        }
    }

    /// Number of Gauss points for the V direction of the face.
    pub(in crate::brep_gprop_full) fn v_integration_order(&self) -> usize {
        match self.kind {
            SurfKind::Plane => 8,
            _ => 18,
        }
    }

    /// Number of Gauss points along a boundary arc.
    pub(in crate::brep_gprop_full) fn arc_integration_order(&self, arc: &BoundaryArc) -> usize {
        match arc.kind {
            ArcKind::Line => 4,
            ArcKind::Circle | ArcKind::Other => 18,
        }
    }

    /// Surface integration order for the adaptive Gauss (SIntOrder).
    pub(in crate::brep_gprop_full) fn s_int_order(&self, eps: f64) -> usize {
        let (nu, nv) = match self.kind {
            SurfKind::Plane => (1, 1),
            SurfKind::Cylinder | SurfKind::Cone => (2, 1),
            SurfKind::Sphere | SurfKind::Torus | SurfKind::Other => (2, 2),
        };
        let sc = s_coeff(eps);
        let n = ((sc * (nu.max(nv) + 1) as f64).ceil() as usize).clamp(1, GPM);
        n
    }

    /// Number of U subintervals (SUIntSubs).
    #[allow(dead_code)]
    pub(in crate::brep_gprop_full) fn s_u_int_subs(&self) -> usize {
        match self.kind {
            SurfKind::Plane => 1,
            SurfKind::Cylinder | SurfKind::Cone | SurfKind::Sphere | SurfKind::Torus => 3,
            SurfKind::Other => 1,
        }
    }

    /// Number of V subintervals (SVIntSubs).
    #[allow(dead_code)]
    pub(in crate::brep_gprop_full) fn s_v_int_subs(&self) -> usize {
        match self.kind {
            SurfKind::Plane | SurfKind::Cylinder | SurfKind::Cone => 1,
            SurfKind::Sphere => 2,
            SurfKind::Torus => 3,
            SurfKind::Other => 1,
        }
    }

    /// U knot values (UKnots).
    pub(in crate::brep_gprop_full) fn u_knots(&self) -> Vec<f64> {
        match self.kind {
            SurfKind::Cylinder | SurfKind::Cone | SurfKind::Sphere | SurfKind::Torus => {
                vec![0.0, 2.0 * PI / 3.0, 4.0 * PI / 3.0, 2.0 * PI]
            }
            _ => vec![self.u1, self.u2],
        }
    }

    /// V knot values (VKnots).
    pub(in crate::brep_gprop_full) fn v_knots(&self) -> Vec<f64> {
        match self.kind {
            SurfKind::Sphere => vec![-PI / 2.0, 0.0, PI / 2.0],
            SurfKind::Torus => vec![0.0, 2.0 * PI / 3.0, 4.0 * PI / 3.0, 2.0 * PI],
            _ => vec![self.v1, self.v2],
        }
    }

    /// Boundary-arc L integration order (LIntOrder).
    #[allow(dead_code)]
    pub(in crate::brep_gprop_full) fn l_int_order(&self, _eps: f64, arc: &BoundaryArc) -> usize {
        let nl = match arc.kind {
            ArcKind::Line => 1,
            ArcKind::Circle => 6,
            ArcKind::Other => 9,
        };
        let ns = self.s_int_order(1.0) as f64;
        let nl = (nl as f64).max(ns);
        ((nl + 1.0).ceil() as usize).clamp(1, GPM)
    }

    /// Boundary-arc L subintervals (LIntSubs).
    #[allow(dead_code)]
    pub(in crate::brep_gprop_full) fn l_int_subs(&self, arc: &BoundaryArc) -> usize {
        match arc.kind {
            ArcKind::Line | ArcKind::Other => 1,
            ArcKind::Circle => 3,
        }
    }

    /// Boundary-arc L knots (LKnots).
    pub(in crate::brep_gprop_full) fn l_knots(&self, arc: &BoundaryArc) -> Vec<f64> {
        match arc.kind {
            ArcKind::Circle => vec![0.0, 2.0 * PI / 3.0, 4.0 * PI / 3.0, 2.0 * PI],
            _ => vec![arc.a, arc.b],
        }
    }

    /// Signed area of the UV boundary polygon (sampled from the pcurves).
    /// Positive ⇒ CCW, negative ⇒ CW.
    pub(in crate::brep_gprop_full) fn uv_polygon_signed_area(&self) -> f64 {
        let mut pts: Vec<GpPnt2d> = Vec::new();
        for arc in &self.arcs {
            for k in 0..=8 {
                let t = arc.a + (arc.b - arc.a) * k as f64 / 8.0;
                pts.push(arc.value(self.surface.as_ref(), t));
            }
        }
        let n = pts.len();
        if n < 3 {
            return 0.0;
        }
        let mut area = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            area += pts[i].x() * pts[j].y() - pts[j].x() * pts[i].y();
        }
        0.5 * area
    }

    /// +1 when the wire is CCW in UV, −1 when CW.
    pub(in crate::brep_gprop_full) fn wire_sign(&self) -> f64 {
        if self.uv_polygon_signed_area() >= 0.0 {
            1.0
        } else {
            -1.0
        }
    }
}
