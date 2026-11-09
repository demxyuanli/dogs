use super::prelude::*;
use super::statics::*;


/// Torus primitive. Source: `BRepPrimAPI_MakeTorus`.
#[derive(Debug, Clone)]
pub struct BRepPrimTorus {
    pub solid: Solid,
    pub major_radius: f64,
    pub minor_radius: f64,
}


impl BRepPrimTorus {
    /// Torus with given major (centerline) and minor (tube) radii.
    ///
    /// The single lateral face carries the `BRepPrim_OneAxis::LateralWire`
    /// structure of a revolution primitive (`BRepPrim_OneAxis.cxx:660-684`).
    /// `BRepPrim_Torus` is a `BRepPrim_Revolution(Axes, 0, 2*pi)`
    /// (`BRepPrim_Torus.cxx:31-43`) whose meridian is a **closed** circle
    /// (`SetMeridian`, `cxx:63-75`: `Geom_Circle(gp_Ax2(Axes.Location +
    /// major*Axes.X, -Axes.Y, Axes.X), minor)`), so `MeridianClosed()` is true and
    /// `TopEdge()`/`BottomEdge()` are the *same* edge as are `StartEdge()`/
    /// `EndEdge()`. DRAW `dump` on `ptorus t 5 2` gives the wire
    /// `-6 +5 +6 -5` — [equator(E6) reversed, seam(E5) forward, equator forward,
    /// seam reversed] — with 2 edges, 1 vertex and 4 pcurves:
    ///
    /// * E6 the equator: `Circle(center 0, axis +Z, radius major+minor)`,
    ///   pcurves `gp_Lin2d((0, 0), +U)` and `gp_Lin2d((0, 2*pi), +U)`;
    /// * E5 the seam: `Circle(center (major, 0, 0), axis -Y, radius minor)`,
    ///   pcurves `gp_Lin2d((2*pi, 0), +V)` and `gp_Lin2d((0, 0), +V)`;
    ///
    /// both over one full period, sharing the vertex at
    /// `(major + minor, 0, 0)`.
    pub fn make_torus(major: f64, minor: f64) -> Self {
        assert!(
            major > 0.0 && minor > 0.0,
            "BRepPrimTorus::make_torus: radii must be positive"
        );
        let b = TopoBuilder::new();

        // `BRepPrim_Torus::SetMeridian` (`BRepPrim_Torus.cxx:63-75`).
        let mer_ax = GpAx3::new(
            GpPnt::new(major, 0.0, 0.0),
            dir(0.0, -1.0, 0.0),
            &dir(1.0, 0.0, 0.0),
        )
        .expect("torus meridian axis");
        let eq_ax = GpAx3::new(GpPnt::zero(), dir(0.0, 0.0, 1.0), &dir(1.0, 0.0, 0.0))
            .expect("torus axis");
        let v_seam = b.make_vertex(GpPnt::new(major + minor, 0.0, 0.0), occt_core::precision::CONFUSION);

        let mut seam = b.make_edge(circle_curve(&mer_ax, minor), 0.0, 2.0 * PI);
        b.add_edge_vertices(&mut seam, &v_seam, &v_seam);
        let mut equator = b.make_edge(circle_curve(&eq_ax, major + minor), 0.0, 2.0 * PI);
        b.add_edge_vertices(&mut equator, &v_seam, &v_seam);

        // `LateralWire` (`cxx:666-679`): TopEdge(false -> REVERSED),
        // EndEdge(true -> FORWARD), BottomEdge(true -> FORWARD),
        // StartEdge(false -> REVERSED) = the DRAW wire `-6 +5 +6 -5`.
        let mut equator_rev = equator.clone();
        equator_rev.0.reverse();
        let mut seam_rev = seam.clone();
        seam_rev.0.reverse();
        let wire = b.make_wire(&[equator_rev, seam.clone(), equator.clone(), seam_rev]);

        let face = b.make_face(
            Arc::new(GeomTorus::new(
                GpTorus::new(GpAx3::standard(), major, minor).expect("torus radii"),
            )),
            &[wire],
        );

        // pcurves (`BRepPrim_OneAxis::LateralFace`, `cxx:388-439`; the DRAW dump
        // lists Curve2d 1/2 on the equator and 3/4 on the seam).
        {
            let reg = GeometryRegistry::global();
            let fk = GeometryRegistry::shape_key(&face.0);
            let dir_u = GpDir2d::new(1.0, 0.0).expect("u direction");
            let dir_v = GpDir2d::new(0.0, 1.0).expect("v direction");
            reg.set_edge_pcurves(
                &equator.0,
                fk,
                vec![
                    Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, 0.0), dir_u)),
                    Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, 2.0 * PI), dir_u)),
                ],
            );
            reg.set_pcurve_range(&equator.0, fk, 0.0, 2.0 * PI);
            reg.set_edge_pcurves(
                &seam.0,
                fk,
                vec![
                    Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(2.0 * PI, 0.0), dir_v)),
                    Arc::new(Geom2dLine::from_pnt_dir(GpPnt2d::new(0.0, 0.0), dir_v)),
                ],
            );
            reg.set_pcurve_range(&seam.0, fk, 0.0, 2.0 * PI);
        }

        let shell = b.make_shell(&[face]);
        let solid = b.make_solid(&[shell]);
        Self { solid, major_radius: major, minor_radius: minor }
    }

    pub fn volume(&self) -> f64 {
        2.0 * PI * PI * self.major_radius * self.minor_radius.powi(2)
    }

    pub fn surface_area(&self) -> f64 {
        4.0 * PI * PI * self.major_radius * self.minor_radius
    }
}
