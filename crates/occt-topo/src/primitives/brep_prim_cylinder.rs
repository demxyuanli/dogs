use super::prelude::*;
use super::statics::*;


/// Cylinder primitive. Source: `BRepPrimAPI_MakeCylinder`.
#[derive(Debug, Clone)]
pub struct BRepPrimCylinder {
    pub solid: Solid,
    pub radius: f64,
    pub height: f64,
}


impl BRepPrimCylinder {
    /// Cylinder of given radius and height, axis along Z.
    /// Boundary: 3 faces (bottom/top caps + lateral), 3 edges
    /// (two cap circles + a seam line), 2 vertices (seam endpoints).
    pub fn make_cylinder(radius: f64, height: f64) -> Self {
        assert!(
            radius > 0.0 && height > 0.0,
            "BRepPrimCylinder::make_cylinder: radius and height must be positive"
        );
        let b = TopoBuilder::new();
        let ax = GpAx3::new(GpPnt::zero(), dir(0.0, 0.0, 1.0), &dir(1.0, 0.0, 0.0))
            .expect("cylinder axis");

        let bottom = GpPnt::new(radius, 0.0, 0.0);
        let top = GpPnt::new(radius, 0.0, height);
        // BRepPrim_Builder::MakeVertex (BRepPrim_Builder.cxx:129-132) is
        // BRep_Builder::MakeVertex(V, P, Precision::Confusion()), so the seam
        // vertices carry the 1e-7 tolerance — it is what
        // BOPAlgo_WireSplitter::Tolerance2D (_1.cxx:859-881) turns into the UV
        // closing tolerance. A zero tolerance makes every anIsSameV2d test fail
        // (aD2 < 0), so the seam loops cannot close.
        let v_bottom = b.make_vertex(bottom, occt_core::precision::CONFUSION);
        let v_top = b.make_vertex(top, occt_core::precision::CONFUSION);

        // Ring circles: parameter [0, 2π], seam at θ = 0 → shared vertices.
        let mut bottom_circle = b.make_edge(circle_curve(&ax, radius), 0.0, 2.0 * PI);
        b.add_edge_vertices(&mut bottom_circle, &v_bottom, &v_bottom);

        let mut top_ax = ax;
        top_ax.set_location(GpPnt::new(0.0, 0.0, height));
        let mut top_circle = b.make_edge(circle_curve(&top_ax, radius), 0.0, 2.0 * PI);
        b.add_edge_vertices(&mut top_circle, &v_top, &v_top);

        // Seam line connecting the two seam vertices.
        let mut seam = b.make_edge(line_curve(&bottom, &top), 0.0, height);
        b.add_edge_vertices(&mut seam, &v_bottom, &v_top);

        // Bottom cap. `BRepPrim_OneAxis::BottomFace` (`BRepPrim_OneAxis.cxx:488
        // -503`) builds the plane from the *untranslated* axes — normal +Z —
        // then `ReverseFace`s it, and `BottomWire` (`cxx:751-770`) stores the
        // circle REVERSED (`AddWireEdge(..., false)`, `cxx:761`). Building the
        // cap with a -Z normal instead leaves the face Forward but mirrors the
        // plane's UV frame (`plane()` picks X = (0,0,1) for a -Z normal), so the
        // cap's circle p-curve runs *opposite* to the edge's 3-D
        // parameterization; `BRep_Tool::CurveOnSurface` requires the same sense
        // and a reversed p-curve makes every direction test downstream misfire
        // (`GetFaceDir`/`GetEdgeOff` pick the wrong face in
        // `BOPAlgo_ShellSplitter::SplitBlock`, T-82).
        let bottom_circle_reversed = Edge(bottom_circle.0.clone().oriented(Orientation::Reversed));
        let mut bottom_wire = b.make_wire(&[bottom_circle_reversed]);
        // `BRepPrim_OneAxis::BottomFace` (`BRepPrim_OneAxis.cxx:501-503`) is
        // `MakeFace` -> **`ReverseFace`** -> `AddFaceWire(BottomWire())`.
        // `TopoDS_Builder::Add` (`TopoDS_Builder.cxx:77-80`) composes the
        // parent's REVERSED orientation into the child it stores, so OCCT ends
        // up with **face REV / wire REV / edge REV**. Reversing the wire here
        // (before the add) reproduces that; reversing only the face afterwards
        // would leave the wire FORWARD, and the boundary view seen through
        // `BRepGProp_Domain::Init(F.Oriented(TopAbs_FORWARD))`
        // (`BRepGProp_Domain.lxx:39-42`) would then be REV instead of FORWARD
        // (a REVERSED cap integrates to -pi instead of +pi).
        bottom_wire.0.reverse();
        let mut bottom_face = b.make_face(
            Arc::new(GeomPlane::new(plane(GpPnt::zero(), dir(0.0, 0.0, 1.0)))),
            &[bottom_wire],
        );
        bottom_face.0.reverse();

        // Top cap: planar face bounded by the top circle.
        let top_wire = b.make_wire(&[top_circle.clone()]);
        let top_face = b.make_face(
            Arc::new(GeomPlane::new(plane(GpPnt::new(0.0, 0.0, height), dir(0.0, 0.0, 1.0)))),
            &[top_wire],
        );

        // Lateral face: the cylinder surface. Wire order and orientations follow
        // `BRepPrim_OneAxis::LateralWire` (`BRepPrim_OneAxis.cxx:660-684`):
        // `AddWireEdge(TopEdge(), false)`, `AddWireEdge(EndEdge(), true)`,
        // `AddWireEdge(BottomEdge(), true)`, `AddWireEdge(StartEdge(), false)`.
        // `BRepPrim_Builder::AddWireEdge` (`BRepPrim_Builder.cxx:184-192`)
        // applies `EE.Reverse()` when `direct == false`, so the wire is
        // **top circle REVERSED, seam FORWARD (u = 2*pi), bottom circle FORWARD,
        // seam REVERSED (u = 0)**. For a full revolution `HasSides()` is false
        // (`BRepPrim_OneAxis.cxx:331-334`), so `StartEdge()` aliases
        // `EndEdge()` (`cxx:976-979`, `:1028-1031`) and the same seam TShape
        // appears twice. The U winding of the loop is zero (bottom -2pi + top
        // +2pi, the seam contributing +h -h): the occurrence orientation written
        // by `TopoDSToStep_MakeStepWire.cxx:262`
        // (`OrientedEdge->Init(..., anEdge.Orientation() == TopAbs_FORWARD)`)
        // keeps the pcurve loop free of any periodic shift that
        // `ShapeFix_Wire::FixShifted` (`ShapeFix_Wire.cxx:1661-2125`) would have
        // to remove.
        //
        // The order also drives `BOPAlgo_WireSplitter`: the edges of the face
        // enter `aLE` in this order, which fixes the vertex slots of the
        // `mySmartMap` and therefore the entry edge of each `Path` walk.
        let mut lateral_top = top_circle.clone();
        lateral_top.0.reverse();
        let mut lateral_seam_back = seam.clone();
        lateral_seam_back.0.reverse();
        let lateral_wire = b.make_wire(&[
            lateral_top,
            seam.clone(),
            bottom_circle.clone(),
            lateral_seam_back,
        ]);
        let lateral_face = b.make_face(Arc::new(GeomCylinder::new(
            GpCylinder::new(ax, radius).expect("cylinder radius"),
        )), &[lateral_wire]);

        // pcurves on the lateral face (`BRepPrim_OneAxis::LateralFace`,
        // `BRepPrim_OneAxis.cxx:388-439`, reached through
        // `BRepPrim_Cylinder::MakeEmptyLateralFace`). `myVMin = 0`,
        // `myVMax = height`, `myAngle = 2*pi`, `myMeridianOffset = 0` and
        // `HasSides()` is false for a full revolution, so:
        //   * `ETOP`/`EBOTTOM` (the cap circles) are `gp_Lin2d((0, VMax|VMin),
        //     +X)` — u free;
        //   * `ESTART`/`EEND` (the seam) use the *closed* form
        //     `SetPCurve(E, F, c1, c2)` with `c1 = gp_Lin2d((myAngle,
        //     -offset), +Y)` and `c2 = gp_Lin2d((0, -offset), +Y)`
        //     (`cxx:434-439`). The two pcurves are what make
        //     `BRep_Tool::IsClosed(seam, lateral)` true
        //     (`BRep_Tool.cxx:820-840`: a curve representation on the closed
        //     surface) and what `BRep_Tool::CurveOnSurface` selects between for
        //     a REVERSED edge (`cxx:354-356`).
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
                &bottom_circle,
                Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, 0.0), dir_u),
                (0.0, 2.0 * PI),
            );
            attach(
                &top_circle,
                Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, height), dir_u),
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
            reg.set_pcurve_range(&seam.0, face_key, 0.0, height);

            // Cap pcurves. `BRepPrim_OneAxis::TopFace` / `BottomFace`
            // (`BRepPrim_OneAxis.cxx:465-468`) put a *circle* of the cap radius
            // in the cap plane's own UV (`gp_Circ2d(gp_Ax2d(gp_Pnt2d(0,0), X),
            // R)`). Without it the closed cap wire degenerates to a straight UV
            // line (`dv/dl == 0`), whose boundary term vanishes and the cap
            // contributes nothing to `BRepGProp`.
            for (cap, circle) in [(&bottom_face, &bottom_circle), (&top_face, &top_circle)] {
                let cap_key = GeometryRegistry::shape_key(&cap.0);
                reg.set_edge_pcurve(
                    &circle.0,
                    cap_key,
                    Arc::new(Geom2dCircle::new(GpCirc2d::new(GpAx22d::new(GpPnt2d::new(0.0, 0.0), dir_u, dir_v).expect("cap frame"), radius))),
                );
                reg.set_pcurve_range(&circle.0, cap_key, 0.0, 2.0 * PI);
            }
        }

        let shell = b.make_shell(&[bottom_face, top_face, lateral_face]);
        let solid = b.make_solid(&[shell]);
        Self { solid, radius, height }
    }

    pub fn volume(&self) -> f64 {
        PI * self.radius * self.radius * self.height
    }

    pub fn lateral_area(&self) -> f64 {
        2.0 * PI * self.radius * self.height
    }
}
