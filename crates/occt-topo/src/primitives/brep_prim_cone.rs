use super::prelude::*;
use super::statics::*;


/// Cone primitive. Source: `BRepPrimAPI_MakeCone`.
#[derive(Debug, Clone)]
pub struct BRepPrimCone {
    pub solid: Solid,
    pub radius: f64,
    pub height: f64,
}


impl BRepPrimCone {
    /// Right cone of base radius `radius` at z = 0 and apex at z = `height`
    /// (`BRepPrimAPI_MakeCone(R1, R2, H)` with `R2 = 0`, `BRepPrim_Cone.hxx:71`).
    ///
    /// Boundary: 2 faces (base cap + lateral), 3 edges (base circle, seam,
    /// degenerate apex) and 2 vertices (base seam point + apex). That is
    /// `BRepPrim_OneAxis` with `HasBottom()` true and `HasTop()` false:
    /// `MeridianOnAxis(VMax)` holds at the apex, so the edge at VMax is the
    /// degenerate apex edge (`BRepPrim_OneAxis.cxx:1205-1213`), and `HasSides()`
    /// is false for a full revolution, so `StartEdge()` aliases `EndEdge()`
    /// (`cxx:976-979`).
    pub fn make_cone(radius: f64, height: f64) -> Self {
        assert!(
            radius > 0.0 && height > 0.0,
            "BRepPrimCone::make_cone: radius and height must be positive"
        );
        let b = TopoBuilder::new();
        let ax = GpAx3::new(GpPnt::zero(), dir(0.0, 0.0, 1.0), &dir(1.0, 0.0, 0.0))
            .expect("cone axis");

        // `BRepPrim_Cone::SetParameters` (`BRepPrim_Cone.cxx:170-173`):
        // `myHalfAngle = ATan((R2 - R1) / H)` (negative here: the cone closes
        // upward, `gp_Cone.hxx:32-34` puts the apex on the positive side),
        // `myRadius = R1`, `myVMin = 0`, `myVMax = Sqrt(H^2 + (R2 - R1)^2)` = the
        // slant length. `ElSLib::ConeDN` (`ElSLib.cxx:182-215`) parameterizes V
        // **along the generatrix**: `R(V) = Radius + V*sin(SAngle)`,
        // `Z(V) = V*cos(SAngle)`.
        let semi_angle = (-radius / height).atan();
        let v_max = (height * height + radius * radius).sqrt();

        let apex = GpPnt::new(0.0, 0.0, height);
        let base_p = GpPnt::new(radius, 0.0, 0.0);
        let v_apex = b.make_vertex(apex, occt_core::precision::CONFUSION);
        let v_base = b.make_vertex(base_p, occt_core::precision::CONFUSION);

        // BottomEdge: the base circle at z = 0 (`BRepPrim_OneAxis::BottomEdge`).
        let mut base_circle = b.make_edge(circle_curve(&ax, radius), 0.0, 2.0 * PI);
        b.add_edge_vertices(&mut base_circle, &v_base, &v_base);

        // StartEdge() == EndEdge(): the meridian at angle 0, the seam from the
        // base seam point up to the apex, ranged on the slant length.
        let mut seam = b.make_edge(line_curve(&base_p, &apex), 0.0, v_max);
        b.add_edge_vertices(&mut seam, &v_base, &v_apex);

        // TopEdge: the degenerate apex edge (`ETOP`, `BRepPrim_OneAxis.cxx:1205
        // -1213`); its 3D range is the angle range [0, 2*pi]
        // (`SetParameters`, `cxx:407`, `:418`).
        let mut top = b.make_degenerated_edge(&apex);
        b.add_edge_vertices(&mut top, &v_apex, &v_apex);
        GeometryRegistry::global().set_edge_range(&top.0, 0.0, 2.0 * PI);

        // Base cap — same construction as the cylinder's bottom cap
        // (`BRepPrim_OneAxis::BottomFace`, `BRepPrim_OneAxis.cxx:488-503`;
        // `BottomWire`, `cxx:751-770`): +Z plane, circle stored REVERSED, wire
        // REVERSED before the add and `ReverseFace` after it, which is the
        // FORWARD boundary view `BRepGProp_Domain` integrates. See
        // `make_cylinder` for why a -Z plane is not equivalent.
        let mut base_circle_reversed = base_circle.clone();
        base_circle_reversed.0.reverse();
        let mut base_wire = b.make_wire(&[base_circle_reversed]);
        base_wire.0.reverse();
        let mut base_face = b.make_face(
            Arc::new(GeomPlane::new(plane(GpPnt::zero(), dir(0.0, 0.0, 1.0)))),
            &[base_wire],
        );
        base_face.0.reverse();

        // Lateral: conical surface. `LateralWire` (`BRepPrim_OneAxis.cxx:660
        // -684`) is `AddWireEdge(TopEdge(), false)`, `AddWireEdge(EndEdge(),
        // true)`, `AddWireEdge(BottomEdge(), true)`, `AddWireEdge(StartEdge(),
        // false)`, and `BRepPrim_Builder::AddWireEdge` (`BRepPrim_Builder.cxx:
        // 184-192`) reverses the edge when `direct` is false: the wire is
        // [apex REVERSED, seam FORWARD, base circle FORWARD, seam REVERSED], the
        // seam being one TShape used with both orientations.
        let mut top_rev = top.clone();
        top_rev.0.reverse();
        let mut seam_back = seam.clone();
        seam_back.0.reverse();
        let lat_wire = b.make_wire(&[top_rev, seam.clone(), base_circle.clone(), seam_back]);

        let lateral_face = b.make_face(Arc::new(GeomCone::new(
            GpCone::new(ax, radius, semi_angle).expect("cone half angle"),
        )), &[lat_wire]);

        // pcurves (`BRepPrim_OneAxis::LateralFace`, `BRepPrim_OneAxis.cxx:388
        // -439`), with `myVMin = 0`, `myVMax = v_max`, `myAngle = 2*pi` and
        // `myMeridianOffset = 0`:
        //   * `ETOP` / `EBOTTOM` are `gp_Lin2d((0, VMax|VMin), +X)` — the u
        //     direction, one full period;
        //   * `ESTART` / `EEND` (the seam) use the *closed* form
        //     `SetPCurve(E, F, c1, c2)` with `c1 = gp_Lin2d((myAngle, -offset),
        //     +Y)` and `c2 = gp_Lin2d((0, -offset), +Y)` (`cxx:434-439`).
        {
            let reg = GeometryRegistry::global();
            let face_key = GeometryRegistry::shape_key(&lateral_face.0);
            let dir_u = GpDir2d::new(1.0, 0.0).expect("u direction");
            let dir_v = GpDir2d::new(0.0, 1.0).expect("v direction");
            let attach = |edge: &Edge, pc: Geom2dLine, range: (f64, f64)| {
                reg.set_edge_pcurve(&edge.0, face_key, Arc::new(pc));
                reg.set_pcurve_range(&edge.0, face_key, range.0, range.1);
            };
            attach(
                &top,
                Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, v_max), dir_u),
                (0.0, 2.0 * PI),
            );
            attach(
                &base_circle,
                Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, 0.0), dir_u),
                (0.0, 2.0 * PI),
            );
            reg.set_edge_pcurves(
                &seam.0,
                face_key,
                vec![
                    Arc::new(Geom2dLine::from_pnt_dir(
                        GpPnt2d::new(2.0 * PI, 0.0),
                        dir_v,
                    )),
                    Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, 0.0), dir_v)),
                ],
            );
            reg.set_pcurve_range(&seam.0, face_key, 0.0, v_max);

            // Cap pcurve: a circle of the cap radius in the cap plane's own UV
            // (`BRepPrim_OneAxis::BottomFace`, `cxx:465-468`). Without it the
            // closed cap wire degenerates to a straight UV line, whose boundary
            // term vanishes and the cap contributes nothing to `BRepGProp`.
            let cap_key = GeometryRegistry::shape_key(&base_face.0);
            reg.set_edge_pcurve(
                &base_circle.0,
                cap_key,
                Arc::new(Geom2dCircle::new(GpCirc2d::new(
                    GpAx22d::new(GpPnt2d::new(0.0, 0.0), dir_u, dir_v).expect("cap frame"),
                    radius,
                ))),
            );
            reg.set_pcurve_range(&base_circle.0, cap_key, 0.0, 2.0 * PI);
        }

        let shell = b.make_shell(&[base_face, lateral_face]);
        let solid = b.make_solid(&[shell]);
        Self { solid, radius, height }
    }

    pub fn volume(&self) -> f64 {
        PI * self.radius * self.radius * self.height / 3.0
    }
}
