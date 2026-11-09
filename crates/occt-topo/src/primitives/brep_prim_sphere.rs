use super::prelude::*;
use super::statics::*;


/// Sphere primitive. Source: `BRepPrimAPI_MakeSphere`.
#[derive(Debug, Clone)]
pub struct BRepPrimSphere {
    pub solid: Solid,
    pub radius: f64,
}


impl BRepPrimSphere {
    /// Sphere of given radius centered at the origin.
    ///
    /// The single lateral face carries the `BRepPrim_OneAxis::LateralWire`
    /// structure of a revolution primitive (`BRepPrim_OneAxis.cxx:660-684`):
    /// `AddWireEdge(TopEdge(), false)`, `AddWireEdge(EndEdge(), true)`,
    /// `AddWireEdge(BottomEdge(), true)`, `AddWireEdge(StartEdge(), false)`,
    /// where the reversed flags leave [north pole REVERSED, seam FORWARD, south
    /// pole FORWARD, seam REVERSED]. The pole edges are degenerate and ranged
    /// `[0, myAngle]` (`SetParameters`, `cxx:407`, `:418`); `HasSides()` is false
    /// for a full revolution, so the meridian is one edge (`StartEdge()` ==
    /// `EndEdge()`, `cxx:976-979`). pcurves are the UV
    /// lines of `BRepPrim_OneAxis::LateralFace` (`cxx:389-439`):
    /// `gp_Lin2d((0, ±π/2), +X)` over `[0, 2π]` for the poles — that is what
    /// makes the UV loop enclose the whole sphere — and the closed form
    /// `gp_Lin2d((0 or 2π, -2π), +Y)` over `[3π/2, 5π/2]` for the meridian.
    /// Meridian geometry follows `BRepPrim_Sphere::SetMeridian`
    /// (`BRepPrim_Sphere.cxx:71-85`).
    pub fn make_sphere(radius: f64) -> Self {
        assert!(
            radius > 0.0,
            "BRepPrimSphere::make_sphere: radius must be positive"
        );
        let b = TopoBuilder::new();
        let r = radius;
        let north = GpPnt::new(0.0, 0.0, r);
        let south = GpPnt::new(0.0, 0.0, -r);
        let v_north = b.make_vertex(north, occt_core::precision::CONFUSION);
        let v_south = b.make_vertex(south, occt_core::precision::CONFUSION);

        let mer_ax = GpAx3::new(GpPnt::zero(), dir(0.0, -1.0, 0.0), &dir(1.0, 0.0, 0.0))
            .expect("meridian axis");
        let (t0, t1) = (1.5 * PI, 2.5 * PI);
        // StartEdge() == EndEdge(): for a full revolution `HasSides()` is false
        // (`BRepPrim_OneAxis.cxx:331-334`), so `StartEdge()` takes the sharing
        // branch `myEdges[ESTART] = myEdges[EEND]`
        // (`BRepPrim_OneAxis.cxx:976-979`) and the wire carries that one seam
        // TShape twice, with opposite orientations.
        let mut seam = b.make_edge(circle_curve(&mer_ax, r), t0, t1);
        b.add_edge_vertices(&mut seam, &v_south, &v_north);

        let mut top = b.make_degenerated_edge(&north);
        b.add_edge_vertices(&mut top, &v_north, &v_north);
        let mut bottom = b.make_degenerated_edge(&south);
        b.add_edge_vertices(&mut bottom, &v_south, &v_south);
        GeometryRegistry::global().set_edge_range(&top.0, 0.0, 2.0 * PI);
        GeometryRegistry::global().set_edge_range(&bottom.0, 0.0, 2.0 * PI);

        // `LateralWire` order and orientations (`cxx:660-684`):
        // `AddWireEdge(TopEdge(), false)` -> REVERSED, `AddWireEdge(EndEdge(),
        // true)` -> FORWARD, `AddWireEdge(BottomEdge(), true)` -> FORWARD,
        // `AddWireEdge(StartEdge(), false)` -> REVERSED (`BRepPrim_Builder::
        // AddWireEdge`, `BRepPrim_Builder.cxx:184-192`).
        let mut top_r = top.clone();
        top_r.0.reverse();
        let mut seam_back = seam.clone();
        seam_back.0.reverse();
        let wire = b.make_wire(&[top_r, seam.clone(), bottom.clone(), seam_back]);

        let surface: Arc<dyn Surface> = Arc::new(GeomSphere::new(
            GpSphere::new(GpAx3::standard(), radius).expect("sphere radius"),
        ));
        let face = b.make_face(surface, &[wire]);

        // pcurves (`BRepPrim_OneAxis::LateralFace`, `cxx:389-439`).
        let face_key = GeometryRegistry::shape_key(&face.0);
        let reg = GeometryRegistry::global();
        let dir_u = GpDir2d::new(1.0, 0.0).expect("u direction");
        let dir_v = GpDir2d::new(0.0, 1.0).expect("v direction");
        let attach = |edge: &Edge, pc: Geom2dLine, range: (f64, f64)| {
            reg.set_edge_pcurve(&edge.0, face_key, Arc::new(pc));
            reg.set_pcurve_range(&edge.0, face_key, range.0, range.1);
        };
        attach(
            &top,
            Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, 0.5 * PI), dir_u),
            (0.0, 2.0 * PI),
        );
        attach(
            &bottom,
            Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, -0.5 * PI), dir_u),
            (0.0, 2.0 * PI),
        );
        // The seam is one edge (`StartEdge()` == `EndEdge()`), so it carries the
        // *closed* form `SetPCurve(E, F, c1, c2)` with `c1 = gp_Lin2d((myAngle,
        // -offset), +Y)` and `c2 = gp_Lin2d((0, -offset), +Y)`
        // (`BRepPrim_OneAxis.cxx:434-439`): the surface V is the meridian
        // parameter shifted by `myMeridianOffset`, i.e. `t - 2*pi`.
        reg.set_edge_pcurves(
            &seam.0,
            face_key,
            vec![
                Arc::new(Geom2dLine::from_pnt_dir(
                    GpPnt2d::new(2.0 * PI, -2.0 * PI),
                    dir_v,
                )),
                Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, -2.0 * PI), dir_v)),
            ],
        );
        reg.set_pcurve_range(&seam.0, face_key, t0, t1);

        let shell = b.make_shell(&[face]);
        let solid = b.make_solid(&[shell]);
        Self { solid, radius }
    }

    pub fn volume(&self) -> f64 {
        4.0 / 3.0 * PI * self.radius.powi(3)
    }

    pub fn surface_area(&self) -> f64 {
        4.0 * PI * self.radius * self.radius
    }
}
