//! `Plate_Plate` — the variational spline solver: a two-variable function
//! satisfying a set of constraints while minimising an energy criterion.
//!
//! Source: `Plate_Plate.hxx:48-167`, `Plate_Plate.cxx:34-2001`,
//! `Plate_Plate.lxx`.
//!
//! Layout:
//! - this file: the data (`cxx:143-165`), the `Load` overloads (`cxx:195-277`),
//!   `SolveTI` dispatch (`cxx:286-364`), `Init` / `IsDone` / `Continuity` /
//!   `CoefPol` / `UVBox` / `UVConstraints` / `SetPolynomialPartOnly`, and
//!   `Polm` (`lxx:20-47`),
//! - [`fundamental`]: `SolEm` (`cxx:1375-1883`),
//! - [`evaluate`]: `Evaluate` / `EvaluateDerivative` (`cxx:1257-1331`),
//! - [`solve_ti1`] / [`solve_ti2`] / [`solve_ti3`]: the three solve cases
//!   (`cxx:370-1056`), built on [`fill_xyz`]'s `fillXYZmatrix` (`cxx:1060-1216`).
//!
//! Translation notes:
//! - OCCT keeps `solution` / `points` / `deru` / `derv` as raw `void*` arrays
//!   sized inside `SolveTI1/2/3`; the port uses `Vec`s whose emptiness stands
//!   for the `nullptr` checks in `Evaluate` (`cxx:1260`, `:1299`) and whose
//!   length is `n_dim` / `n_el`.
//! - `SolEm` caches `Uold` / `Vold` / `U2` / `R` / `L` through a
//!   `const_cast` (`cxx:1380`). The port keeps the same cache in `Cell`s so
//!   `evaluate` stays `&self` while writing the cache, exactly like OCCT.
//! - `Message_ProgressScope` (`cxx:443`, `:561`, `:916`) has no counterpart;
//!   the `aScope.UserBreak()` arms (`cxx:446-450`, `:564-568`, `:919-923`)
//!   are therefore UNPORTED. They only trigger on a user break.
//! - `math_Gauss` maps to [`occt_math::MathGauss`]; the OCCT failure arms that
//!   add `1.0e-8` to the diagonal and retry with `pivot_max = 1.0e-18`
//!   (`cxx:451-467`, `:570-586`, `:925-944`) are reproduced verbatim.

mod evaluate;
mod fill_xyz;
mod fundamental;
mod solve_ti1;
mod solve_ti2;
mod solve_ti3;

use std::cell::Cell;

use occt_core::gp::{GpXY, GpXyz};

use crate::plate::{
    FreeGtoCConstraint, GlobalTranslationConstraint, GtoCConstraint, LinearScalarConstraint,
    LinearXyzConstraint, LineConstraint, PinpointConstraint, PlaneConstraint,
    SampledCurveConstraint,
};

/// `Plate_Plate` (`hxx:50-167`).
///
/// `Clone` is `Plate_Plate(const Plate_Plate&)` (`cxx:55-123`): the vectors are
/// copied and the constraint sequences are copied, exactly like the OCCT copy
/// constructor. The `Uold` / `Vold` memoisation fields are copied rather than
/// reset to `1.0e20`; that only changes how often `SolEm` recomputes, never the
/// values it returns.
#[derive(Debug, Clone)]
pub struct Plate {
    order: i32,
    n_el: usize,
    n_dim: usize,
    solution: Vec<GpXyz>,
    points: Vec<GpXY>,
    deru: Vec<i32>,
    derv: Vec<i32>,
    ok: bool,
    my_constraints: Vec<PinpointConstraint>,
    my_lxyz_constraints: Vec<LinearXyzConstraint>,
    my_lscalar_constraints: Vec<LinearScalarConstraint>,
    ddu: [f64; 10],
    ddv: [f64; 10],
    max_constraint_order: i32,
    polynomial_part_only: bool,
    uold: Cell<f64>,
    vold: Cell<f64>,
    u2: Cell<f64>,
    r_cache: Cell<f64>,
    l: Cell<f64>,
}

impl Default for Plate {
    /// `Plate_Plate()` (`cxx:34-53`): `OK = false`, `Uold = Vold = 1.0e20`,
    /// `U2 = R = L = 0`, zeroed `ddu` / `ddv`, `PolynomialPartOnly = false`.
    fn default() -> Self {
        Self {
            order: 0,
            n_el: 0,
            n_dim: 0,
            solution: Vec::new(),
            points: Vec::new(),
            deru: Vec::new(),
            derv: Vec::new(),
            ok: false,
            my_constraints: Vec::new(),
            my_lxyz_constraints: Vec::new(),
            my_lscalar_constraints: Vec::new(),
            ddu: [0.0; 10],
            ddv: [0.0; 10],
            max_constraint_order: 0,
            polynomial_part_only: false,
            uold: Cell::new(1.0e20),
            vold: Cell::new(1.0e20),
            u2: Cell::new(0.0),
            r_cache: Cell::new(0.0),
            l: Cell::new(0.0),
        }
    }
}

impl Plate {
    /// `Plate_Plate()` (`cxx:34-53`).
    pub fn new() -> Self {
        Self::default()
    }

    // ---- Load (`cxx:195-277`) ----

    /// `Load(const Plate_PinpointConstraint&)` (`cxx:195-204`).
    pub fn load_pinpoint(&mut self, p_const: PinpointConstraint) {
        self.ok = false;
        self.n_el += 1;
        self.my_constraints.push(p_const);
        let ordre_const = p_const.idu() + p_const.idv();
        if self.max_constraint_order < ordre_const {
            self.max_constraint_order = ordre_const;
        }
    }

    /// `Load(const Plate_LinearXYZConstraint&)` (`cxx:206-221`).
    pub fn load_linear_xyz(&mut self, l_xyz_const: LinearXyzConstraint) {
        self.ok = false;
        self.n_el += l_xyz_const.row_length();
        for ppc in l_xyz_const.ppc() {
            let ordre_const = ppc.idu() + ppc.idv();
            if self.max_constraint_order < ordre_const {
                self.max_constraint_order = ordre_const;
            }
        }
        self.my_lxyz_constraints.push(l_xyz_const);
    }

    /// `Load(const Plate_LinearScalarConstraint&)` (`cxx:223-237`).
    pub fn load_linear_scalar(&mut self, l_scalar_const: LinearScalarConstraint) {
        self.ok = false;
        self.n_el += l_scalar_const.row_length();
        for ppc in l_scalar_const.ppc() {
            let ordre_const = ppc.idu() + ppc.idv();
            if self.max_constraint_order < ordre_const {
                self.max_constraint_order = ordre_const;
            }
        }
        self.my_lscalar_constraints.push(l_scalar_const);
    }

    /// `Load(const Plate_LineConstraint&)` (`cxx:239-242`).
    pub fn load_line(&mut self, l_const: &LineConstraint) {
        self.load_linear_scalar(l_const.lsc().clone());
    }

    /// `Load(const Plate_PlaneConstraint&)` (`cxx:244-247`).
    pub fn load_plane(&mut self, p_const: &PlaneConstraint) {
        self.load_linear_scalar(p_const.lsc().clone());
    }

    /// `Load(const Plate_SampledCurveConstraint&)` (`cxx:249-252`).
    pub fn load_sampled_curve(&mut self, sc_const: &SampledCurveConstraint) {
        self.load_linear_xyz(sc_const.lxyzc().clone());
    }

    /// `Load(const Plate_GtoCConstraint&)` (`cxx:254-259`).
    pub fn load_gto_c(&mut self, gto_c_const: &GtoCConstraint) {
        for i in 0..gto_c_const.nb_ppc() {
            self.load_pinpoint(*gto_c_const.get_ppc(i));
        }
    }

    /// `Load(const Plate_FreeGtoCConstraint&)` (`cxx:261-272`).
    pub fn load_free_gto_c(&mut self, f_gto_c_const: &FreeGtoCConstraint) {
        for i in 0..f_gto_c_const.nb_ppc() {
            self.load_pinpoint(*f_gto_c_const.get_ppc(i));
        }
        for i in 0..f_gto_c_const.nb_lsc() {
            self.load_linear_scalar(f_gto_c_const.lsc(i).clone());
        }
    }

    /// `Load(const Plate_GlobalTranslationConstraint&)` (`cxx:274-277`).
    pub fn load_global_translation(&mut self, gt_const: &GlobalTranslationConstraint) {
        self.load_linear_xyz(gt_const.lxyzc().clone());
    }

    // ---- SolveTI (`cxx:286-364`) ----

    /// `SolveTI(ord = 4, anisotropie = 1.0)` (`cxx:286-364`).
    ///
    /// The OCCT default arguments are `ord = 4` / `anisotropie = 1.0`; pass
    /// them explicitly here.
    ///
    /// Deviation (documented, no control-flow change): the `UVBox` call
    /// (`cxx:314`) reads `myConstraints(1)` even when only `LinearXYZ` /
    /// `LinearScalar` constraints were loaded, which raises
    /// `Standard_OutOfRange` in OCCT. The port leaves `OK = false` and returns
    /// instead of panicking; every other guard (`ord <= 1`, `ord > 9`,
    /// `n_el < 1`, the `anisotropie` range, `du < 1.0e-10`, `dv < 1.0e-10`) is
    /// reproduced.
    pub fn solve_ti(&mut self, ord: i32, anisotropie: f64) {
        let iteration_number = 0;
        self.ok = false;
        self.order = ord;
        if ord <= 1 {
            return;
        }
        if ord > 9 {
            return;
        }
        if self.n_el < 1 {
            return;
        }
        if anisotropie < 1.0e-6 {
            return;
        }
        if anisotropie > 1.0e+6 {
            return;
        }

        if self.my_constraints.is_empty() {
            return;
        }
        // computation of the bounding box of the 2d PPconstraints
        let (xmin, xmax, ymin, ymax) = self.uv_box();

        let mut du = 0.5 * (xmax - xmin);
        if anisotropie > 1.0 {
            du *= anisotropie;
        }
        if du < 1.0e-10 {
            return;
        }
        self.ddu[0] = 1.0;
        for i in 1..=9 {
            self.ddu[i] = self.ddu[i - 1] / du;
        }

        let mut dv = 0.5 * (ymax - ymin);
        if anisotropie < 1.0 {
            dv /= anisotropie;
        }
        if dv < 1.0e-10 {
            return;
        }
        self.ddv[0] = 1.0;
        for i in 1..=9 {
            self.ddv[i] = self.ddv[i - 1] / dv;
        }

        if self.my_lscalar_constraints.is_empty() {
            if self.my_lxyz_constraints.is_empty() {
                self.solve_ti1(iteration_number);
            } else {
                self.solve_ti2(iteration_number);
            }
        } else {
            self.solve_ti3(iteration_number);
        }
    }

    // ---- accessors ----

    /// `IsDone()` (`cxx:1218-1221`).
    pub fn is_done(&self) -> bool {
        self.ok
    }

    /// `Order()` — the `order` member set by `SolveTI`.
    pub fn get_order(&self) -> i32 {
        self.order
    }

    /// `Continuity()` (`cxx:1367-1370`).
    pub fn continuity(&self) -> i32 {
        2 * self.order - 3 - self.max_constraint_order
    }

    /// `SetPolynomialPartOnly(PPOnly)` (`cxx:1997-2001`).
    pub fn set_polynomial_part_only(&mut self, pp_only: bool) {
        self.polynomial_part_only = pp_only;
    }

    /// `CoefPol()` (`cxx:1339-1357`): the `order` x `order` array of power
    /// basis coefficients of the polynomial part, `[iu][iv]`.
    pub fn coef_pol(&self) -> Vec<Vec<GpXyz>> {
        let o = self.order.max(0) as usize;
        let mut coefs = vec![vec![GpXyz::zero(); o]; o];
        let mut i = self.n_el;
        for iu in 0..o {
            let mut iv = 0usize;
            while iu + iv < o {
                coefs[iu][iv] =
                    self.solution[i].multiplied(self.ddu[iu] * self.ddv[iv]);
                i += 1;
                iv += 1;
            }
        }
        coefs
    }

    /// `UVBox(UMin, UMax, VMin, VMax)` (`cxx:1887-1981`): the bounding box of
    /// every constraint's 2d point, widened to `Bmin = 1.0e-3` per axis.
    pub fn uv_box(&self) -> (f64, f64, f64, f64) {
        let bmin = 1.0e-3;
        let first = self.my_constraints[0].pnt2d();
        let mut umin = first.x();
        let mut umax = umin;
        let mut vmin = first.y();
        let mut vmax = vmin;

        for c in self.my_constraints.iter().skip(1) {
            let x = c.pnt2d().x();
            if x < umin {
                umin = x;
            }
            if x > umax {
                umax = x;
            }
            let y = c.pnt2d().y();
            if y < vmin {
                vmin = y;
            }
            if y > vmax {
                vmax = y;
            }
        }

        for c in &self.my_lxyz_constraints {
            for ppc in c.ppc() {
                let x = ppc.pnt2d().x();
                if x < umin {
                    umin = x;
                }
                if x > umax {
                    umax = x;
                }
                let y = ppc.pnt2d().y();
                if y < vmin {
                    vmin = y;
                }
                if y > vmax {
                    vmax = y;
                }
            }
        }

        for c in &self.my_lscalar_constraints {
            for ppc in c.ppc() {
                let x = ppc.pnt2d().x();
                if x < umin {
                    umin = x;
                }
                if x > umax {
                    umax = x;
                }
                let y = ppc.pnt2d().y();
                if y < vmin {
                    vmin = y;
                }
                if y > vmax {
                    vmax = y;
                }
            }
        }

        if umax - umin < bmin {
            let um = 0.5 * (umin + umax);
            umin = um - 0.5 * bmin;
            umax = um + 0.5 * bmin;
        }
        if vmax - vmin < bmin {
            let vm = 0.5 * (vmin + vmax);
            vmin = vm - 0.5 * bmin;
            vmax = vm + 0.5 * bmin;
        }
        (umin, umax, vmin, vmax)
    }

    /// `UVConstraints(Seq)` (`cxx:1985-1994`): the `(0, 0)` order pinpoint
    /// constraints.
    pub fn uv_constraints(&self) -> Vec<GpXY> {
        let mut seq = Vec::new();
        for c in &self.my_constraints {
            if c.idu() == 0 && c.idv() == 0 {
                seq.push(*c.pnt2d());
            }
        }
        seq
    }

    /// `Init()` (`cxx:1232-1254`): back to the post-construction state,
    /// including `OK = true` (`cxx:1252`).
    pub fn init(&mut self) {
        self.my_constraints.clear();
        self.my_lxyz_constraints.clear();
        self.my_lscalar_constraints.clear();
        self.solution.clear();
        self.points.clear();
        self.deru.clear();
        self.derv.clear();
        self.order = 0;
        self.n_el = 0;
        self.n_dim = 0;
        self.ok = true;
        self.max_constraint_order = 0;
    }

    /// `Polm(point2d, iu, iv, idu, idv)` (`lxx:20-47`).
    fn polm(&self, point2d: GpXY, iu: i32, iv: i32, idu: i32, idv: i32) -> f64 {
        if idu > iu {
            return 0.0;
        }
        if idv > iv {
            return 0.0;
        }
        let u = point2d.x();
        let v = point2d.y();

        let mut value = 1.0;

        let degu = iu - idu;
        for _ in 0..degu {
            value *= u;
        }
        for i in (degu + 1..=iu).rev() {
            value *= i as f64;
        }

        let degv = iv - idv;
        for _ in 0..degv {
            value *= v;
        }
        for i in (degv + 1..=iv).rev() {
            value *= i as f64;
        }

        // The `ddu[iu] * ddv[iv]` product is kept: it does not change the
        // computed values for the polynomial part of `Sol`, but it improves
        // the conditioning of the matrix (`lxx:41-45`).
        value * self.ddu[iu as usize] * self.ddv[iv as usize]
    }

    // ---- raw array accessors (`Plate_Plate.lxx:49-81`) ----

    fn deru(&self, index: usize) -> i32 {
        self.deru[index]
    }

    fn derv(&self, index: usize) -> i32 {
        self.derv[index]
    }

    /// `Solution(index)`.
    fn solution_at(&self, index: usize) -> GpXyz {
        self.solution[index]
    }

    /// `Points(index)`.
    fn point_at(&self, index: usize) -> GpXY {
        self.points[index]
    }
}
