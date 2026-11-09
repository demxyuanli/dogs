//! `Contap_ContAna` - the analytic contour lines of the three quadrics that
//! have them in closed form: sphere, cylinder and cone.
//!
//! Sources: `Contap_ContAna.hxx:35-91`, `.cxx:29-483`, `.lxx:24-64`.
//!
//! Rust cannot overload `Perform`, so each OCCT overload is a separately named
//! method (`perform_sphere_dir`, `perform_sphere_dir_angle`,
//! `perform_sphere_eye`, and the same trio for the cylinder and the cone).
//! The `StdFail_NotDone` / `Standard_DomainError` / `Standard_OutOfRange`
//! guards of `NbContours` / `TypeContour` / `Circle` / `Line` answer `Err`.

use occt_core::gp::{GpAx2, GpCirc, GpCone, GpCylinder, GpDir, GpLin, GpPnt, GpSphere, GpXyz};

use occt_topo::geom_abs_curve_type::GeomAbsCurveType;

/// `static const double Tolpetit = 1.e-8;` (`cxx:29`).
const TOLPETIT: f64 = 1e-8;

/// `0.9999999999999` - the collinearity threshold of the sphere arms
/// (`cxx:45`, `:66`, `:102`).
const COLLINEAR: f64 = 0.9999999999999;

/// `Contap_ContAna` (`hxx:35-91`).
#[derive(Debug, Clone)]
pub struct ContapContAna {
    /// `done` (`hxx:74`).
    done: bool,
    /// `nbSol` (`hxx:75`).
    nb_sol: i32,
    /// `typL` (`hxx:76`).
    typ_l: GeomAbsCurveType,
    /// `pt1` .. `pt4` (`hxx:77-80`).
    pt1: GpPnt,
    pt2: GpPnt,
    pt3: GpPnt,
    pt4: GpPnt,
    /// `dir1` .. `dir4` (`hxx:81-84`).
    dir1: GpDir,
    dir2: GpDir,
    dir3: GpDir,
    dir4: GpDir,
    /// `prm` (`hxx:86`): the circle radius, or the `sqrt` factor of the
    /// silhouette generators.
    prm: f64,
}

impl Default for ContapContAna {
    /// `Contap_ContAna()` (`cxx:31-37`).
    fn default() -> Self {
        Self::new()
    }
}

impl ContapContAna {
    /// `Contap_ContAna()` (`cxx:31-37`): not done, no solution,
    /// `GeomAbs_OtherCurve`, `prm = 0`.
    pub fn new() -> Self {
        Self {
            done: false,
            nb_sol: 0,
            typ_l: GeomAbsCurveType::OtherCurve,
            pt1: GpPnt::default(),
            pt2: GpPnt::default(),
            pt3: GpPnt::default(),
            pt4: GpPnt::default(),
            dir1: GpDir::default(),
            dir2: GpDir::default(),
            dir3: GpDir::default(),
            dir4: GpDir::default(),
            prm: 0.0,
        }
    }

    /// `Perform(S, D)` (`cxx:39-56`): the sphere seen along `D`.
    pub fn perform_sphere_dir(&mut self, s: &GpSphere, d: &GpDir) {
        self.done = false;
        self.typ_l = GeomAbsCurveType::Circle;
        self.pt1 = s.location();
        self.dir1 = *d;
        let xaxis = *s.x_axis().direction();
        let yaxis = *s.y_axis().direction();
        self.dir2 = if d.dot(&xaxis).abs() < COLLINEAR {
            d.crossed(&xaxis).unwrap_or_default()
        } else {
            d.crossed(&yaxis).unwrap_or_default()
        };
        self.prm = s.radius;
        self.nb_sol = 1;
        self.done = true;
    }

    /// `Perform(S, D, Angle)` (`cxx:58-77`): the sphere seen along `D` at a
    /// draft `Angle`.
    pub fn perform_sphere_dir_angle(&mut self, s: &GpSphere, d: &GpDir, angle: f64) {
        self.done = false;
        self.typ_l = GeomAbsCurveType::Circle;
        self.dir1 = *d;
        let xaxis = *s.x_axis().direction();
        let yaxis = *s.y_axis().direction();
        self.dir2 = if d.dot(&xaxis).abs() < COLLINEAR {
            d.crossed(&xaxis).unwrap_or_default()
        } else {
            d.crossed(&yaxis).unwrap_or_default()
        };
        let alpha = if s.is_direct() { angle } else { -angle };
        // `pt1 = S.Location() - S.Radius() * sin(alpha) * D`.
        self.pt1 = GpPnt::from_xyz(
            &s.location()
                .xyz()
                .subtracted(&d.xyz().multiplied(s.radius * alpha.sin())),
        );
        self.prm = s.radius * alpha.cos();
        self.nb_sol = 1;
        self.done = true;
    }

    /// `Perform(S, Eye)` (`cxx:79-113`): the sphere seen from a point.
    pub fn perform_sphere_eye(&mut self, s: &GpSphere, eye: &GpPnt) {
        self.done = false;
        let radius = s.radius;
        let dist = eye.distance(&s.location());
        if dist <= radius {
            self.nb_sol = 0;
        } else {
            self.prm = radius * (1.0 - radius * radius / (dist * dist)).sqrt();
            if self.prm < TOLPETIT {
                self.nb_sol = 0;
            } else {
                let locxyz = *s.location().xyz();
                self.dir1 = GpDir::from_xyz(&eye.xyz().subtracted(&locxyz)).unwrap_or_default();
                self.pt1 =
                    GpPnt::from_xyz(&locxyz.added(&self.dir1.xyz().multiplied(radius * radius / dist)));
                let xaxis = *s.x_axis().direction();
                let yaxis = *s.y_axis().direction();
                self.dir2 = if self.dir1.dot(&xaxis).abs() < COLLINEAR {
                    self.dir1.crossed(&xaxis).unwrap_or_default()
                } else {
                    self.dir1.crossed(&yaxis).unwrap_or_default()
                };
                self.nb_sol = 1;
            }
            self.typ_l = GeomAbsCurveType::Circle;
        }
        self.done = true;
    }

    /// `Perform(C, D)` (`cxx:115-138`): the two cylinder silhouette
    /// generators along `D`.
    pub fn perform_cylinder_dir(&mut self, c: &GpCylinder, d: &GpDir) {
        self.done = false;
        let pos = c.position();
        let mut normale = pos.direction().xyz().crossed(d.xyz());
        if normale.modulus() <= 1e-15 {
            self.nb_sol = 0;
        } else {
            normale.normalize();
            self.typ_l = GeomAbsCurveType::Line;
            self.dir1 = pos.direction();
            self.dir2 = self.dir1;
            let loc = *c.location().xyz();
            self.pt1 = GpPnt::from_xyz(&loc.added(&normale.multiplied(c.radius)));
            self.pt2 = GpPnt::from_xyz(&loc.subtracted(&normale.multiplied(c.radius)));
            self.nb_sol = 2;
        }
        self.done = true;
    }

    /// `Perform(C, D, Angle)` (`cxx:140-198`): the two cylinder generators of a
    /// draft contour.
    pub fn perform_cylinder_dir_angle(&mut self, c: &GpCylinder, d: &GpDir, angle: f64) {
        self.done = false;
        let pos = c.position();
        let mut coefcos = d.dot(&pos.x_direction());
        let mut coefsin = d.dot(&pos.y_direction());
        let coefcst = (std::f64::consts::FRAC_PI_2 + angle).cos();

        let norm1 = coefcos * coefcos + coefsin * coefsin;

        if coefcst.abs() < norm1.sqrt() {
            self.typ_l = GeomAbsCurveType::Line;
            self.nb_sol = 2;
            self.dir1 = pos.direction();
            self.dir2 = self.dir1;

            if !c.is_direct() {
                // `cxx:157-161`: the normal is inverted.
                coefcos = -coefcos;
                coefsin = -coefsin;
            }

            self.prm = (norm1 - coefcst * coefcst).sqrt();
            let cost0 = (coefcos * coefcst - coefsin * self.prm) / norm1;
            let cost1 = (coefcos * coefcst + coefsin * self.prm) / norm1;
            let sint0 = (coefcos * self.prm + coefsin * coefcst) / norm1;
            let sint1 = (-coefcos * self.prm + coefsin * coefcst) / norm1;

            let xdir = *pos.x_direction().xyz();
            let ydir = *pos.y_direction().xyz();
            let loc = *c.location().xyz();
            let dirxyz = linear2(cost0, &xdir, sint0, &ydir).multiplied(c.radius);
            self.pt1 = GpPnt::from_xyz(&loc.added(&dirxyz));
            let dirxyz = linear2(cost1, &xdir, sint1, &ydir).multiplied(c.radius);
            self.pt2 = GpPnt::from_xyz(&loc.added(&dirxyz));
        } else {
            self.nb_sol = 0;
        }

        self.done = true;
    }

    /// `Perform(C, Eye)` (`cxx:200-226`): the two cylinder silhouette
    /// generators seen from a point.
    pub fn perform_cylinder_eye(&mut self, c: &GpCylinder, eye: &GpPnt) {
        self.done = false;
        let radius = c.radius;
        let theaxis = GpLin::new(c.axis());
        let dist = theaxis.distance(eye);
        if dist <= radius {
            self.nb_sol = 0;
        } else {
            self.typ_l = GeomAbsCurveType::Line;
            self.prm = radius * (1.0 - radius * radius / (dist * dist)).sqrt();
            self.dir1 = theaxis.direction();
            self.dir2 = self.dir1;
            // `gp_XYZ axeye(theaxis.Normal(Eye).Direction().XYZ())`.
            let axeye = *theaxis.normal(eye).direction().xyz();
            // `gp_XYZ normale((theaxis.Direction().Crossed(axeye)).XYZ())`.
            let normale = theaxis.direction().xyz().crossed(&axeye);
            let loc = *c.location().xyz();
            self.pt1 = GpPnt::from_xyz(&loc.added(&axeye.multiplied(radius * radius / dist)));
            self.pt2 = GpPnt::from_xyz(&self.pt1.xyz().subtracted(&normale.multiplied(self.prm)));
            self.pt1 = GpPnt::from_xyz(&self.pt1.xyz().added(&normale.multiplied(self.prm)));
            self.nb_sol = 2;
        }
        self.done = true;
    }

    /// `Perform(C, D)` (`cxx:228-287`): the two cone silhouette generators
    /// along `D`, both through the apex.
    pub fn perform_cone_dir(&mut self, c: &GpCone, d: &GpDir) {
        self.done = false;
        let tgtalpha = c.semi_angle().tan();
        let pos = c.position();

        let coefcos = d.dot(&pos.x_direction());
        let coefsin = d.dot(&pos.y_direction());
        let coefcst = d.dot(&c.axis().direction()) * tgtalpha;

        let norm1 = coefcos * coefcos + coefsin * coefsin;

        if coefcst.abs() < norm1.sqrt() {
            self.typ_l = GeomAbsCurveType::Line;
            self.nb_sol = 2;
            self.pt1 = c.apex();
            self.pt2 = self.pt1;

            self.prm = (norm1 - coefcst * coefcst).sqrt();
            let cost0 = (coefcos * coefcst - coefsin * self.prm) / norm1;
            let cost1 = (coefcos * coefcst + coefsin * self.prm) / norm1;
            let sint0 = (coefcos * self.prm + coefsin * coefcst) / norm1;
            let sint1 = (-coefcos * self.prm + coefsin * coefcst) / norm1;

            let xdir = *pos.x_direction().xyz();
            let ydir = *pos.y_direction().xyz();
            let zdir = *c.axis().direction().xyz();
            let dirxyz = linear3(cost0, &xdir, sint0, &ydir, 1.0 / tgtalpha, &zdir);
            self.dir1 = GpDir::from_xyz(&dirxyz).unwrap_or_default();
            self.pt1 = GpPnt::from_xyz(&self.pt1.xyz().added(&dirxyz));
            let dirxyz = linear3(cost1, &xdir, sint1, &ydir, 1.0 / tgtalpha, &zdir);
            self.dir2 = GpDir::from_xyz(&dirxyz).unwrap_or_default();
            self.pt2 = GpPnt::from_xyz(&self.pt2.xyz().added(&dirxyz));
        } else {
            self.nb_sol = 0;
        }
        self.done = true;
    }

    /// `Perform(C, D, Angle)` (`cxx:289-390`): up to four cone generators for
    /// a draft contour (two for each side of the apex).
    pub fn perform_cone_dir_angle(&mut self, c: &GpCone, d: &GpDir, angle: f64) {
        self.done = false;
        self.nb_sol = 0;

        let ang = c.semi_angle();
        let cosa = ang.cos();
        let sina = ang.sin();

        let pos = c.position();
        let coefcos = d.dot(&pos.x_direction());
        let coefsin = d.dot(&pos.y_direction());

        let coefcst1 = (std::f64::consts::FRAC_PI_2 + angle).cos();

        let norm1 = coefcos * coefcos + coefsin * coefsin;

        let coefnz = d.dot(&c.axis().direction()) * sina;
        let coefcst = (coefcst1 + coefnz) / cosa;

        if coefcst.abs() < norm1.sqrt() {
            self.typ_l = GeomAbsCurveType::Line;
            self.nb_sol += 2;
            self.pt1 = c.apex();
            self.pt2 = self.pt1;

            self.prm = (norm1 - coefcst * coefcst).sqrt();
            let cost0 = (coefcos * coefcst - coefsin * self.prm) / norm1;
            let cost1 = (coefcos * coefcst + coefsin * self.prm) / norm1;
            let sint0 = (coefcos * self.prm + coefsin * coefcst) / norm1;
            let sint1 = (-coefcos * self.prm + coefsin * coefcst) / norm1;

            let xdir = *pos.x_direction().xyz();
            let ydir = *pos.y_direction().xyz();
            let mut zdir = *c.axis().direction().xyz();
            if !c.is_direct() {
                zdir.reverse();
            }
            let dirxyz = linear3(cost0, &xdir, sint0, &ydir, cosa / sina, &zdir);
            self.dir1 = GpDir::from_xyz(&dirxyz).unwrap_or_default();
            self.pt1 = GpPnt::from_xyz(&self.pt1.xyz().added(&dirxyz));
            let dirxyz = linear3(cost1, &xdir, sint1, &ydir, cosa / sina, &zdir);
            self.dir2 = GpDir::from_xyz(&dirxyz).unwrap_or_default();
            self.pt2 = GpPnt::from_xyz(&self.pt2.xyz().added(&dirxyz));
        }

        let coefcst = (coefcst1 - coefnz) / cosa;

        if coefcst.abs() < norm1.sqrt() {
            self.typ_l = GeomAbsCurveType::Line;
            self.nb_sol += 2;
            self.pt3 = c.apex();
            self.pt4 = self.pt3;

            self.prm = (norm1 - coefcst * coefcst).sqrt();
            let cost0 = (coefcos * coefcst - coefsin * self.prm) / norm1;
            let cost1 = (coefcos * coefcst + coefsin * self.prm) / norm1;
            let sint0 = (coefcos * self.prm + coefsin * coefcst) / norm1;
            let sint1 = (-coefcos * self.prm + coefsin * coefcst) / norm1;

            let xdir = *pos.x_direction().xyz();
            let ydir = *pos.y_direction().xyz();
            let mut zdir = *c.axis().direction().xyz();
            if !c.is_direct() {
                zdir.reverse();
            }
            let dirxyz = linear3(cost0, &xdir, sint0, &ydir, -cosa / sina, &zdir);
            self.dir3 = GpDir::from_xyz(&dirxyz).unwrap_or_default();
            self.pt3 = GpPnt::from_xyz(&self.pt3.xyz().added(&dirxyz));
            let dirxyz = linear3(cost1, &xdir, sint1, &ydir, -cosa / sina, &zdir);
            self.dir4 = GpDir::from_xyz(&dirxyz).unwrap_or_default();
            self.pt4 = GpPnt::from_xyz(&self.pt4.xyz().added(&dirxyz));
            if self.nb_sol == 2 {
                self.pt1 = self.pt3;
                self.pt2 = self.pt4;
                self.dir1 = self.dir3;
                self.dir2 = self.dir4;
            }
        }

        self.done = true;
    }

    /// `Perform(C, Eye)` (`cxx:392-455`): the two cone silhouette generators
    /// seen from a point.
    pub fn perform_cone_eye(&mut self, c: &GpCone, eye: &GpPnt) {
        self.done = false;
        let tgtalpha = c.semi_angle().tan();
        let pos = c.position();

        let apexeye = eye.xyz().subtracted(c.apex().xyz());

        let coefcos = apexeye.dot(pos.x_direction().xyz());
        let coefsin = apexeye.dot(pos.y_direction().xyz());
        let coefcst = apexeye.dot(c.axis().direction().xyz()) * tgtalpha;

        let norm1 = coefcos * coefcos + coefsin * coefsin;
        let norm2 = (coefcos * coefcos + coefsin * coefsin).sqrt();

        if coefcst.abs() < norm2 {
            self.typ_l = GeomAbsCurveType::Line;
            self.nb_sol = 2;
            self.pt1 = c.apex();
            self.pt2 = self.pt1;

            self.prm = (norm1 - coefcst * coefcst).sqrt();
            let cost0 = (coefcos * coefcst - coefsin * self.prm) / norm1;
            let cost1 = (coefcos * coefcst + coefsin * self.prm) / norm1;
            let sint0 = (coefcos * self.prm + coefsin * coefcst) / norm1;
            let sint1 = (-coefcos * self.prm + coefsin * coefcst) / norm1;

            let xdir = *pos.x_direction().xyz();
            let ydir = *pos.y_direction().xyz();
            let zdir = *c.axis().direction().xyz();
            let dirxyz = linear3(cost0, &xdir, sint0, &ydir, 1.0 / tgtalpha, &zdir);
            self.dir1 = GpDir::from_xyz(&dirxyz).unwrap_or_default();
            self.pt1 = GpPnt::from_xyz(&self.pt1.xyz().added(&dirxyz));
            let dirxyz = linear3(cost1, &xdir, sint1, &ydir, 1.0 / tgtalpha, &zdir);
            self.dir2 = GpDir::from_xyz(&dirxyz).unwrap_or_default();
            self.pt2 = GpPnt::from_xyz(&self.pt2.xyz().added(&dirxyz));
        } else {
            self.nb_sol = 0;
        }
        self.done = true;
    }

    /// `IsDone()` (`lxx:26-29`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `NbContours()` (`lxx:31-37`): `StdFail_NotDone` before `Perform`.
    pub fn nb_contours(&self) -> Result<i32, &'static str> {
        if !self.done {
            return Err("Contap_ContAna::NbContours: StdFail_NotDone");
        }
        Ok(self.nb_sol)
    }

    /// `TypeContour()` (`lxx:39-48`): `StdFail_NotDone` before `Perform`,
    /// `Standard_DomainError` for zero solutions.
    pub fn type_contour(&self) -> Result<GeomAbsCurveType, &'static str> {
        if !self.done {
            return Err("Contap_ContAna::TypeContour: StdFail_NotDone");
        }
        if self.nb_sol == 0 {
            return Err("Contap_ContAna::TypeContour: Standard_DomainError");
        }
        Ok(self.typ_l)
    }

    /// `Circle()` (`lxx:52-64`).
    pub fn circle(&self) -> Result<GpCirc, &'static str> {
        if !self.done {
            return Err("Contap_ContAna::Circle: StdFail_NotDone");
        }
        if self.typ_l != GeomAbsCurveType::Circle || self.nb_sol == 0 {
            return Err("Contap_ContAna::Circle: Standard_DomainError");
        }
        let pos = GpAx2::new(self.pt1, self.dir1, self.dir2)
            .map_err(|_| "Contap_ContAna::Circle: degenerate gp_Ax2")?;
        Ok(GpCirc::new(pos, self.prm))
    }

    /// `Line(Index)` (`cxx:458-483`), 1-based.
    pub fn line(&self, index: i32) -> Result<GpLin, &'static str> {
        if !self.done {
            return Err("Contap_ContAna::Line: StdFail_NotDone");
        }
        if self.typ_l != GeomAbsCurveType::Line || self.nb_sol == 0 {
            return Err("Contap_ContAna::Line: Standard_DomainError");
        }
        if index <= 0 || index > self.nb_sol {
            return Err("Contap_ContAna::Line: Standard_OutOfRange");
        }
        let (p, d) = match index {
            1 => (self.pt1, self.dir1),
            2 => (self.pt2, self.dir2),
            3 => (self.pt3, self.dir3),
            4 => (self.pt4, self.dir4),
            _ => return Err("Contap_ContAna::Line: Standard_OutOfRange"),
        };
        Ok(GpLin::from_pnt_dir(p, d))
    }
}

/// `gp_XYZ::SetLinearForm` with two terms (`cxx:184-190`).
fn linear2(a1: f64, v1: &GpXyz, a2: f64, v2: &GpXyz) -> GpXyz {
    v1.multiplied(a1).added(&v2.multiplied(a2))
}

/// `gp_XYZ::SetLinearForm` with three terms (`cxx:274-280`).
fn linear3(a1: f64, v1: &GpXyz, a2: f64, v2: &GpXyz, a3: f64, v3: &GpXyz) -> GpXyz {
    v1.multiplied(a1)
        .added(&v2.multiplied(a2))
        .added(&v3.multiplied(a3))
}
