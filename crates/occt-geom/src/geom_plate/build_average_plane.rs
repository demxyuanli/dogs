//! `GeomPlate_BuildAveragePlane` — the inertial / max-flux average plane.
//!
//! Source: `TKGeomAlgo/GeomPlate/GeomPlate_BuildAveragePlane.hxx`,
//! `GeomPlate_BuildAveragePlane.cxx:47-750` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `NCollection_HArray1<gp_Pnt> myPts` is `Vec<GpPnt>` (index `i - 1`).
//! - `GeomLib::AxeOfInertia` is [`geom_lib_axe_of_inertia`]; the two-argument
//!   overload (`cxx:236`) uses OCCT's default tolerance `1.0e-7`
//!   (`GeomLib.hxx`, `AxeOfInertia(..., const double Tol = 1.0e-7)`).
//! - `ElSLib::Parameters(const gp_Pln&, const gp_Pnt&, U, V)` is
//!   [`plane_parameters`] on the plane's position.
//! - `math_Jacobi::Vector(k, V)` is column `k` of [`Jacobi::vectors`].
//! - `gp_Vec::Rotated` / `Rotate` are a `GpTrsf` rotation about the axis, which
//!   is the same linear map.
//! - `Precision::SquareConfusion()` is [`SQUARE_CONFUSION`].
//! - `RealFirst()` / `RealLast()` are `-f64::MAX` / `f64::MAX`.
//! - The `Standard_NoSuchObject` raises of `Plane()` (`cxx:285-288`) and
//!   `Line()` (`cxx:493-497`) become the `Err` of those accessors.

use occt_core::elib::slib::plane_parameters;
use occt_core::gp::{GpAx2, GpAx3, GpDir, GpPln, GpPnt, GpVec};
use occt_core::precision::SQUARE_CONFUSION;
use occt_math::{Jacobi, MathMatrix};

use crate::geom_lib::geom_lib_axe_of_inertia;
use crate::geom_plate::aij::Aij;
use crate::line::GeomLine;
use crate::plane::GeomPlane;

/// `GeomLib::AxeOfInertia`'s default `Tol` (`GeomLib.hxx`).
const AXE_OF_INERTIA_TOL: f64 = 1.0e-7;

/// `class GeomPlate_BuildAveragePlane` (`hxx:35-105`): "computes an average
/// inertial plane with an array of points".
pub struct BuildAveragePlane {
    /// `myPts` (`hxx:93`).
    pts: Vec<GpPnt>,
    /// `myUmax` (`hxx:94`).
    u_max: f64,
    /// `myVmax` (`hxx:95`).
    v_max: f64,
    /// `myVmin` (`hxx:96`).
    v_min: f64,
    /// `myUmin` (`hxx:97`).
    u_min: f64,
    /// `myPlane` (`hxx:98`).
    plane: Option<GeomPlane>,
    /// `myTol` (`hxx:99`).
    tol: f64,
    /// `myLine` (`hxx:100`).
    line: Option<GeomLine>,
    /// `myOX` (`hxx:101`).
    ox: GpVec,
    /// `myOY` (`hxx:102`).
    oy: GpVec,
    /// `myG` (`hxx:103`).
    g: GpPnt,
    /// `myNbBoundPoints` (`hxx:104`).
    nb_bound_points: i32,
}

impl BuildAveragePlane {
    /// `GeomPlate_BuildAveragePlane(Pts, NbBoundPoints, Tol, POption, NOption)`
    /// (`cxx:47-106`).
    ///
    /// `POption`: 1 = automatic parametrisation, 2 = parametrisation by eigen
    /// vectors. `NOption`: 1 = inertial plane, 2 = plane of max flux.
    pub fn from_points(
        pts: Vec<GpPnt>,
        nb_bound_points: i32,
        tol: f64,
        p_option: i32,
        n_option: i32,
    ) -> Result<Self, String> {
        let mut this = Self {
            pts,
            u_max: 0.0,
            v_max: 0.0,
            v_min: 0.0,
            u_min: 0.0,
            plane: None,
            tol,
            line: None,
            ox: GpVec::zero(),
            oy: GpVec::zero(),
            g: GpPnt::zero(),
            nb_bound_points,
        };
        let oz = this.def_plan(n_option)?;

        if oz.square_magnitude() > 0.0 {
            if p_option == 1 {
                let axe = GpAx2::from_axis(this.g, GpDir::from_vec(&oz)?);
                this.plane = Some(GeomPlane::new(GpPln::new(axe.to_ax3())));
                let pln = this.plane.as_ref().expect("just built").pln();
                this.ox = GpVec::from_xyz(pln.x_axis().direction().xyz());
                this.oy = GpVec::from_xyz(pln.y_axis().direction().xyz());
            } else {
                this.base_plan(&oz)?;
                let n_dir = GpDir::from_vec(&this.ox.crossed(&this.oy))?;
                let u_dir = GpDir::from_vec(&this.ox)?;
                let ax3 = GpAx3::new(this.g, n_dir, &u_dir)?;
                this.plane = Some(GeomPlane::new(GpPln::new(ax3)));
            }
            let pln = this.plane.as_ref().expect("just built").pln().clone();
            let (u, v) = plane_parameters(&pln.position(), &this.g);
            this.u_max = u;
            this.v_max = v;
            this.u_min = u;
            this.v_min = v;
            for i in 0..this.pts.len() {
                let (u, v) = plane_parameters(&pln.position(), &this.pts[i]);
                if this.u_max < u {
                    this.u_max = u;
                }
                if this.u_min > u {
                    this.u_min = u;
                }
                if this.v_max < v {
                    this.v_max = v;
                }
                if this.v_min > v {
                    this.v_min = v;
                }
            }
        }

        if this.is_line() {
            let dir = GpDir::from_vec(&this.ox)?;
            this.line = Some(GeomLine::from_pnt_dir(this.g, dir));
        }
        Ok(this)
    }

    /// `GeomPlate_BuildAveragePlane(Normals, Pts)` (`cxx:108-278`): "creates the
    /// plane from the best vector".
    pub fn from_normals(normals: &[GpVec], pts: Vec<GpPnt>) -> Result<Self, String> {
        let mut this = Self {
            pts,
            u_max: 0.0,
            v_max: 0.0,
            v_min: 0.0,
            u_min: 0.0,
            plane: None,
            tol: 0.0,
            line: None,
            ox: GpVec::zero(),
            oy: GpVec::zero(),
            g: GpPnt::zero(),
            nb_bound_points: 0,
        };
        let nn = normals.len();
        let best_vec = if nn == 1 {
            normals[0]
        } else if nn == 2 {
            let mut best = normals[0].added(&normals[1]);
            let sq_magn = best.square_magnitude();
            if sq_magn < SQUARE_CONFUSION {
                let sq1 = normals[0].square_magnitude();
                let sq2 = normals[1].square_magnitude();
                if sq1 > sq2 {
                    normals[0].normalized()
                } else {
                    normals[1].normalized()
                }
            } else {
                best = best.divide(sq_magn.sqrt());
                best
            }
        } else {
            // The common case (`cxx:146-227`).
            let mut max_angle = 0.0f64;
            for i in 0..nn - 1 {
                for j in i + 1..nn {
                    let angle = normals[i].angle(&normals[j]);
                    if angle > max_angle {
                        max_angle = angle;
                    }
                }
            }
            max_angle *= 1.2;
            max_angle /= 2.0;
            let n_int = 50;
            let count = nn * (nn - 1) / 2;
            let mut opt_vec = vec![GpVec::zero(); count];
            let mut opt_scal = vec![-f64::MAX; count];

            let mut k = 0usize;
            for i in 0..nn - 1 {
                for j in i + 1..nn {
                    let step = max_angle / n_int as f64;
                    let mut vec = normals[i].added(&normals[j]);
                    let sq_magn = vec.square_magnitude();
                    if sq_magn < SQUARE_CONFUSION {
                        k += 1;
                        continue;
                    }
                    vec = vec.divide(sq_magn.sqrt());

                    let cross1 = normals[i].crossed(&normals[j]);
                    let cross2 = vec.crossed(&cross1);
                    let axis = occt_core::gp::GpAx1::new(GpPnt::zero(), GpDir::from_vec(&cross2)?);

                    let mut vec1 = rotated(&vec, &axis, -max_angle)?;
                    for _n in 0..=2 * n_int {
                        vec1 = rotated(&vec1, &axis, step)?;
                        let mut min_scal = f64::MAX;
                        for m in 0..nn {
                            let scal = vec1.dot(&normals[m]);
                            if scal < min_scal {
                                min_scal = scal;
                            }
                        }
                        if min_scal > opt_scal[k] {
                            opt_scal[k] = min_scal;
                            opt_vec[k] = vec1;
                        }
                    }
                    k += 1;
                }
            }

            let mut best_scal = -f64::MAX;
            let mut index = 0usize;
            for k in 0..opt_scal.len() {
                if opt_scal[k] > best_scal {
                    best_scal = opt_scal[k];
                    index = k;
                }
            }
            opt_vec[index]
        };

        // Making the plane myPlane (`cxx:228-278`).
        let (axe, _is_singular) = geom_lib_axe_of_inertia(&this.pts, AXE_OF_INERTIA_TOL);
        let best_dir = GpDir::from_vec(&best_vec)?;
        let x_dir = best_dir.crossed(axe.x_direction())?;
        let x_dir = x_dir.crossed(&best_dir)?;

        let axe3 = GpAx3::new(axe.location(), best_dir, &x_dir)?;
        this.plane = Some(GeomPlane::new(GpPln::new(axe3)));

        // Initializing myUmin, myVmin, myUmax, myVmax (`cxx:245-275`).
        let pln = this.plane.as_ref().expect("just built").pln().clone();
        let (u, v) = plane_parameters(&pln.position(), &axe.location());
        this.u_max = u;
        this.v_max = v;
        this.u_min = u;
        this.v_min = v;
        let norm_vec = GpVec::from_xyz(pln.axis().direction().xyz());
        for i in 0..this.pts.len() {
            let a_vec = GpVec::from_pnts(&pln.location(), &this.pts[i]);
            let norm_vec = norm_vec.multiplied_scalar(a_vec.dot(&norm_vec));
            let projected = this.pts[i].translated_vec(&norm_vec.reversed());
            let (u, v) = plane_parameters(&pln.position(), &projected);
            if u > this.u_max {
                this.u_max = u;
            }
            if u < this.u_min {
                this.u_min = u;
            }
            if v > this.v_max {
                this.v_max = v;
            }
            if v < this.v_min {
                this.v_min = v;
            }
        }

        // Initializing myOX, myOY (`cxx:277-278`).
        this.ox = GpVec::from_xyz(pln.x_axis().direction().xyz());
        this.oy = GpVec::from_xyz(pln.y_axis().direction().xyz());
        Ok(this)
    }

    /// `Plane()` (`cxx:282-288`). Raises when the object is a line.
    pub fn plane(&self) -> Result<&GeomPlane, String> {
        if self.is_line() {
            return Err(
                "Cannot use the function 'GeomPlate_BuildAveragePlane::Plane()', the Object is a \
                 'Geom_Line'"
                    .to_string(),
            );
        }
        self.plane
            .as_ref()
            .ok_or_else(|| "GeomPlate_BuildAveragePlane::Plane is null".to_string())
    }

    /// `Line()` (`cxx:491-497`). Raises when the object is a plane.
    pub fn line(&self) -> Result<&GeomLine, String> {
        if self.is_plane() {
            return Err(
                "Cannot use the function 'GeomPlate_BuildAveragePlane::Line()', the Object is a \
                 'Geom_Plane'"
                    .to_string(),
            );
        }
        self.line
            .as_ref()
            .ok_or_else(|| "GeomPlate_BuildAveragePlane::Line is null".to_string())
    }

    /// `IsPlane()` (`cxx:502-506`).
    pub fn is_plane(&self) -> bool {
        self.ox.crossed(&self.oy).square_magnitude() != 0.0
    }

    /// `IsLine()` (`cxx:508-513`).
    pub fn is_line(&self) -> bool {
        self.ox.crossed(&self.oy).square_magnitude() == 0.0
    }

    /// `MinMaxBox(Umin, Umax, Vmin, Vmax)` (`cxx:292-302`).
    pub fn min_max_box(&self) -> (f64, f64, f64, f64) {
        (self.u_min, self.u_max, self.v_min, self.v_max)
    }

    /// `DefPlan(NOption)` (`cxx:304-382`): defines the average plane's normal
    /// and, for `NOption == 1`, the in-plane axes.
    fn def_plan(&mut self, n_option: i32) -> Result<GpVec, String> {
        let nb = self.pts.len();
        if nb == 0 {
            return Err("GeomPlate_BuildAveragePlane: no points".to_string());
        }
        let mut gb = GpPnt::zero();
        for p in &self.pts {
            gb = GpPnt::new(gb.x() + p.x(), gb.y() + p.y(), gb.z() + p.z());
        }
        self.g = GpPnt::new(gb.x() / nb as f64, gb.y() / nb as f64, gb.z() / nb as f64);

        let mut oz = GpVec::zero();
        if n_option == 1 {
            let (axe, _is_singular) = geom_lib_axe_of_inertia(&self.pts, self.tol);
            self.ox = GpVec::from_xyz(axe.x_direction().xyz());
            self.oy = GpVec::from_xyz(axe.y_direction().xyz());
            oz = GpVec::from_xyz(axe.direction().xyz());

            if self.nb_bound_points != 0 && self.pts.len() as i32 != self.nb_bound_points {
                let sum = self.bound_cross_sum()?;
                if oz.square_magnitude() == 0.0 || sum.square_magnitude() == 0.0 {
                    return Err("GeomPlate_BuildAveragePlane::DefPlan: null direction".to_string());
                }
                let mut the_angle = oz.angle(&sum);
                if the_angle > std::f64::consts::PI / 2.0 {
                    the_angle = std::f64::consts::PI - the_angle;
                }
                if the_angle > std::f64::consts::PI / 3.0 {
                    oz = sum;
                }
            }
        } else if n_option == 2 {
            oz = self.bound_cross_sum()?;
        }
        Ok(oz)
    }

    /// The `A = sum_{i=3..NbBoundPoints} (P(i-1) - P(1)) ^ (P(i) - P(1))` loop
    /// shared by `DefPlan`'s two options (`cxx:336-350`, `:365-379`).
    fn bound_cross_sum(&self) -> Result<GpVec, String> {
        let mut a = GpVec::zero();
        if self.nb_bound_points > self.pts.len() as i32 {
            return Err("Standard_OutOfRange: GeomPlate_BuildAveragePlane::DefPlan".to_string());
        }
        for i in 3..=self.nb_bound_points as usize {
            let b = GpVec::from_pnts(&self.pts[0], &self.pts[i - 2]);
            let c = GpVec::from_pnts(&self.pts[0], &self.pts[i - 1]);
            a = a.added(&b.crossed(&c));
        }
        Ok(a)
    }

    /// `BasePlan(OZ)` (`cxx:387-487`): the in-plane base from the eigenvectors
    /// of the projected inertia matrix.
    fn base_plan(&mut self, oz: &GpVec) -> Result<(), String> {
        let mut m = MathMatrix::new(1, 3, 1, 3);
        let oz_sq = oz.dot(oz);
        if oz_sq == 0.0 {
            return Err("GeomPlate_BuildAveragePlane::BasePlan: null normal".to_string());
        }
        for p in &self.pts {
            let mut proj = GpVec::from_pnts(&self.g, p);
            let scal = proj.dot(oz) / oz_sq;
            proj = proj.subtracted(&oz.multiplied_scalar(scal));
            let (x, y, z) = (proj.x(), proj.y(), proj.z());
            m.set_value(1, 1, m.value(1, 1) + x * x);
            m.set_value(2, 2, m.value(2, 2) + y * y);
            m.set_value(3, 3, m.value(3, 3) + z * z);
            m.set_value(1, 2, m.value(1, 2) + x * y);
            m.set_value(1, 3, m.value(1, 3) + x * z);
            m.set_value(2, 3, m.value(2, 3) + y * z);
        }
        m.set_value(2, 1, m.value(1, 2));
        m.set_value(3, 1, m.value(1, 3));
        m.set_value(3, 2, m.value(2, 3));

        let j = Jacobi::new(&m);
        let n1 = j.values.value(1);
        let n2 = j.values.value(2);
        let n3 = j.values.value(3);

        // Eigenvalue ordering (`cxx:421-471`).
        let r1 = n1.min(n2).min(n3);
        let (m1, m2, m3);
        if r1 == n1 {
            m1 = 1;
            let r2 = n2.min(n3);
            if r2 == n2 {
                m2 = 2;
                m3 = 3;
            } else {
                m2 = 3;
                m3 = 2;
            }
        } else if r1 == n2 {
            m1 = 2;
            let r2 = n1.min(n3);
            if r2 == n1 {
                m2 = 1;
                m3 = 3;
            } else {
                m2 = 3;
                m3 = 1;
            }
        } else {
            m1 = 3;
            let r2 = n1.min(n2);
            if r2 == n1 {
                m2 = 1;
                m3 = 2;
            } else {
                m2 = 2;
                m3 = 1;
            }
        }
        let _ = m1;
        let v2 = column(&j, m2);
        let v3 = column(&j, m3);

        if (n1.abs() <= self.tol && n2.abs() <= self.tol)
            || (n2.abs() <= self.tol && n3.abs() <= self.tol)
            || (n1.abs() <= self.tol && n3.abs() <= self.tol)
        {
            self.ox = v3;
            self.oy = oz.crossed(&self.ox);
        } else {
            self.ox = v3;
            self.oy = v2;
        }
        Ok(())
    }
}

/// `math_Jacobi::Vector(k, V)`: column `k` of the eigenvector matrix.
fn column(j: &Jacobi, k: usize) -> GpVec {
    GpVec::new(
        j.vectors.value(1, k),
        j.vectors.value(2, k),
        j.vectors.value(3, k),
    )
}

/// `gp_Vec::Rotated(const gp_Ax1&, const double)`: a rotation about the axis.
fn rotated(v: &GpVec, axis: &occt_core::gp::GpAx1, angle: f64) -> Result<GpVec, String> {
    let mut t = occt_core::gp::GpTrsf::identity();
    t.set_rotation_ax1(axis, angle).map_err(|e| e.to_string())?;
    Ok(v.transformed(&t))
}

/// `GeomPlate_BuildAveragePlane::HalfSpace` (`cxx:515-750`): narrows the normal
/// sequence to a half-space basis.
///
/// `Normals` and `Bset` are in-out: on `false` OCCT restores them from the
/// snapshot taken at `cxx:524-525`, so the caller sees the values it passed in.
pub fn half_space(
    new_normals: &[GpVec],
    normals: &mut Vec<GpVec>,
    bset: &mut Vec<Aij>,
    lin_tol: f64,
    ang_tol: f64,
) -> bool {
    let square_tol = lin_tol * lin_tol;

    // 1 (`cxx:523-525`).
    let save_normals = normals.clone();
    let save_bset = bset.clone();

    let null_vec = GpVec::zero();
    let mut b1set: Vec<Aij> = Vec::new();
    let mut b2set: Vec<Aij> = Vec::new();

    let mut i = 0usize;
    if normals.is_empty() {
        if new_normals.len() == 1 {
            normals.push(new_normals[0]);
            return true;
        }
        // 2 (`cxx:541-553`).
        let cross = new_normals[0].crossed(&new_normals[1]);
        if cross.square_magnitude() <= square_tol {
            return false;
        }
        let cross = cross.normalized();
        bset.push(Aij::from_parts(1, 2, cross));
        bset.push(Aij::from_parts(2, 1, cross.reversed()));
        normals.push(new_normals[0]);
        normals.push(new_normals[1]);
        i = 2;
    }

    while i < new_normals.len() {
        // 3 (`cxx:557-566`).
        for entry in bset.iter() {
            let scal = entry.vec.dot(&new_normals[i]);
            if scal >= -lin_tol {
                b2set.push(*entry);
            }
        }

        let ii = normals.len() + 1;
        let mut restore = false;
        for j in 0..ii - 1 {
            if normals[j].square_magnitude() == 0.0 {
                continue;
            }
            // 4 (`cxx:575-632`).
            let mut cross = new_normals[i].crossed(&normals[j]);
            if cross.square_magnitude() <= square_tol {
                restore = true;
                break;
            }
            cross = cross.normalized();

            match resolve(&cross, &mut b2set, normals, new_normals, i, ii, square_tol, ang_tol) {
                Verdict::Reject => {
                    restore = true;
                    break;
                }
                Verdict::Handled => {}
                Verdict::New => b1set.push(Aij::from_parts(ii as i32, (j + 1) as i32, cross)),
            }

            cross.reverse();
            match resolve(&cross, &mut b2set, normals, new_normals, i, ii, square_tol, ang_tol) {
                Verdict::Reject => {
                    restore = true;
                    break;
                }
                Verdict::Handled => {}
                Verdict::New => b1set.push(Aij::from_parts(ii as i32, (j + 1) as i32, cross)),
            }
        }
        if restore {
            *normals = save_normals;
            *bset = save_bset;
            return false;
        }

        // 5 (`cxx:690-711`).
        for entry in b1set.iter() {
            let mut is_ge_null = true;
            for n in normals.iter() {
                if n.square_magnitude() == 0.0 {
                    continue;
                }
                if entry.vec.dot(n) < -lin_tol {
                    is_ge_null = false;
                    break;
                }
            }
            if is_ge_null {
                b2set.push(*entry);
            }
        }

        // 6 (`cxx:713-718`).
        if b2set.is_empty() {
            *normals = save_normals;
            *bset = save_bset;
            return false;
        }

        // 7 (`cxx:720-725`).
        std::mem::swap(bset, &mut b2set);
        b2set.clear();
        b1set.clear();
        normals.push(new_normals[i]);

        // 8 (`cxx:727-745`).
        for j in 0..normals.len() {
            if normals[j].square_magnitude() == 0.0 {
                continue;
            }
            let found = bset
                .iter()
                .any(|e| e.ind1 == (j + 1) as i32 || e.ind2 == (j + 1) as i32);
            if !found {
                normals[j] = null_vec;
            }
        }

        i += 1;
    }

    true
}

/// The outcome of one `HalfSpace` "is the cross already in the basis" probe.
enum Verdict {
    /// The caller must restore the snapshot and return `false`.
    Reject,
    /// `isNew = false; break;` — nothing to append.
    Handled,
    /// `isNew` stayed true — append `Aij(ii, j, cross)`.
    New,
}

/// The `for (k = 1; k <= B2set.Length(); k++)` body of `cxx:585-631` (and its
/// `Cross.Reverse()` twin at `cxx:638-684`), which mutates `B2set` in place.
#[allow(clippy::too_many_arguments)]
fn resolve(
    cross: &GpVec,
    b2set: &mut [Aij],
    normals: &[GpVec],
    new_normals: &[GpVec],
    i: usize,
    ii: usize,
    square_tol: f64,
    ang_tol: f64,
) -> Verdict {
    for k in 0..b2set.len() {
        if !b2set[k].vec.is_opposite_ang(&cross.reversed(), ang_tol) {
            continue;
        }
        let ind1 = b2set[k].ind1 as usize;
        let ind2 = b2set[k].ind2 as usize;
        if ind1 == ii || ind2 == ii {
            return Verdict::Handled;
        }
        let cross1 = normals[ind1 - 1].crossed(&new_normals[i]);
        let cross2 = normals[ind2 - 1].crossed(&new_normals[i]);
        if cross1.square_magnitude() <= square_tol || cross2.square_magnitude() <= square_tol {
            return Verdict::Reject;
        }
        if cross1.is_opposite_ang(&cross2, ang_tol) {
            let cross2 = normals[ind1 - 1].crossed(&normals[ind2 - 1]);
            if cross1.is_opposite_ang(&cross2, ang_tol) {
                return Verdict::Reject;
            }
        } else if new_normals[i].angle(&normals[ind1 - 1]) > new_normals[i].angle(&normals[ind2 - 1])
        {
            b2set[k].ind2 = ind1 as i32;
            b2set[k].ind1 = ii as i32;
        } else {
            b2set[k].ind1 = ii as i32;
        }
        return Verdict::Handled;
    }
    Verdict::New
}
