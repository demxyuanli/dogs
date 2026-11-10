//! `Contap_TheIWalking` - the `IntWalk_IWalking` instantiation of
//! `Contap_TheIWalking_0.cxx:22-56` (with `Contap_SurfFunction` as
//! `TheIWFunction`, `IntSurf_PathPoint` / `IntSurf_InteriorPoint` as the start
//! points and `Contap_TheIWLineOfTheIWalking` as `IntWalk_TheIWLine`).
//!
//! Sources: `Contap_TheIWalking.hxx:46-227` (state at `:202-227`),
//! `IntWalk_IWalking.lxx:19-59`, `IntWalk_IWalking.gxx:81-3140`.
//!
//! Control flow kept as in OCCT:
//! - `Perform` (`gxx:146-302`) loads `wd1` / `wd2`, sets the resolutions, the
//!   estimated parameter ranges, then runs `ComputeOpenLine` and
//!   `ComputeCloseLine`, and finally the `ToFillHoles` retry loop.
//! - the `wd1` / `wd2` vectors keep the dummy element OCCT appends in `Clear`
//!   (`gxx:110-123`) so that every `wd1[i]` / `wd2[i]` of the template maps to
//!   the same 1-based index here.
//!
//! UNPORTED: `Perform(Pnts1, Func, Caro, Reversed)` (`gxx:307-384`), the
//! overload without interior points. Its only callers are the `IntPatch`
//! walks (`IntPatch_TheIWalking`), which the port drives through
//! `occt_topo::intpatch_iwalking`; the contour walking of `Contap_Contour.cxx:1700`
//! always passes the interior points.

pub(crate) mod types;

pub mod cadrage;
mod close_line;
mod deflection;
mod line_ops;
mod open_line;
mod stop_tests;

use std::collections::HashMap;

use occt_core::gp::{GpDir2d, GpVec};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_math::MathVector;

use occt_topo::int_tools_wline::{u_resolution, v_resolution, PntOn2S};
use occt_topo::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::iw_line::ContapIWLineOfTheIWalking;
use super::surf_function::ContapSurfFunction;
use types::{BndRange, InteriorPoint, PathPoint, WalkingData};

pub use types::{InteriorPoint as ContapWalkingInteriorPoint, PathPoint as ContapWalkingPathPoint};

/// `Contap_TheIWalking` (`hxx:46-227`).
pub struct ContapTheIWalking {
    /// `done` (`hxx:203`).
    done: bool,
    /// `seqSingle` (`hxx:204`): the start points on which no line starts or
    /// ends.
    seq_single: Vec<PathPoint>,
    /// `fleche` (`hxx:205`).
    fleche: f64,
    /// `pas` (`hxx:206`).
    pas: f64,
    /// `tolerance` (`hxx:207`): `(UResolution, VResolution)`.
    vec_tolerance: [f64; 2],
    /// `epsilon` (`hxx:208`): the squared 3D confusion.
    epsilon: f64,
    /// `reversed` (`hxx:209`).
    reversed: bool,
    /// `wd1` (`hxx:210`), 1-based with the dummy of `Clear` at index 0.
    wd1: Vec<WalkingData>,
    /// `wd2` (`hxx:211`), 1-based with the dummy of `Clear` at index 0.
    wd2: Vec<WalkingData>,
    /// `nbMultiplicities` (`hxx:212`), 1-based with the dummy of `Clear`.
    nb_multiplicities: Vec<i32>,
    /// `mySRangeU` (`hxx:213`).
    s_range_u: BndRange,
    /// `mySRangeV` (`hxx:214`).
    s_range_v: BndRange,
    /// `Um` (`hxx:215`).
    um: f64,
    /// `UM` (`hxx:216`).
    umax: f64,
    /// `Vm` (`hxx:217`).
    vm: f64,
    /// `VM` (`hxx:218`).
    vmax: f64,
    /// `previousPoint` (`hxx:219`).
    previous_point: PntOn2S,
    /// `previousd3d` (`hxx:220`).
    previous_d3d: GpVec,
    /// `previousd2d` (`hxx:221`).
    previous_d2d: GpDir2d,
    /// `seqAjout` (`hxx:222`): the lines whose end point may be continued.
    seq_ajout: Vec<i32>,
    /// `seqAlone` (`hxx:223`): the lines with a free end.
    seq_alone: Vec<i32>,
    /// `PointLineLine` (`hxx:224`): `(start point, lines it belongs to)`.
    point_line_line: HashMap<i32, Vec<i32>>,
    /// `lines` (`hxx:225`).
    lines: Vec<ContapIWLineOfTheIWalking>,
    /// `ToFillHoles` (`hxx:226`).
    to_fill_holes: bool,
}

impl ContapTheIWalking {
    /// `Contap_TheIWalking(Epsilon, Deflection, Step, theToFillHoles)`
    /// (`hxx:61-66`, `gxx:82-101`): `tolerance` is allocated empty, `epsilon`
    /// is the squared `Epsilon`.
    pub fn new(epsilon: f64, deflection: f64, step: f64, to_fill_holes: bool) -> Self {
        Self {
            done: false,
            seq_single: Vec::new(),
            fleche: deflection,
            pas: step,
            vec_tolerance: [0.0, 0.0],
            epsilon: epsilon * epsilon,
            reversed: false,
            wd1: Vec::new(),
            wd2: Vec::new(),
            nb_multiplicities: Vec::new(),
            s_range_u: BndRange::void(),
            s_range_v: BndRange::void(),
            um: 0.0,
            umax: 0.0,
            vm: 0.0,
            vmax: 0.0,
            previous_point: PntOn2S {
                p: occt_core::gp::GpPnt::default(),
                u1: 0.0,
                v1: 0.0,
                u2: 0.0,
                v2: 0.0,
            },
            previous_d3d: GpVec::zero(),
            previous_d2d: GpDir2d::new(1.0, 0.0).expect("+X"),
            seq_ajout: Vec::new(),
            seq_alone: Vec::new(),
            point_line_line: HashMap::new(),
            lines: Vec::new(),
            to_fill_holes,
        }
    }

    /// `SetTolerance(Epsilon, Deflection, Step)` (`lxx:20-26`).
    pub fn set_tolerance(&mut self, epsilon: f64, deflection: f64, step: f64) {
        self.fleche = deflection;
        self.pas = step;
        self.epsilon = epsilon * epsilon;
    }

    /// `IsDone()` (`lxx:28-31`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `NbLines()` (`lxx:33-38`); `None` where OCCT raises `StdFail_NotDone`.
    pub fn nb_lines(&self) -> Option<i32> {
        if self.done {
            Some(self.lines.len() as i32)
        } else {
            None
        }
    }

    /// `Value(Index)` (`lxx:40-46`), 1-based.
    pub fn value(&self, index: i32) -> Option<&ContapIWLineOfTheIWalking> {
        if !self.done {
            return None;
        }
        self.lines.get((index.max(1) - 1) as usize)
    }

    /// `NbSinglePnts()` (`lxx:48-53`).
    pub fn nb_single_pnts(&self) -> Option<i32> {
        if self.done {
            Some(self.seq_single.len() as i32)
        } else {
            None
        }
    }

    /// `SinglePnt(Index)` (`lxx:55-60`), 1-based.
    pub fn single_pnt(&self, index: i32) -> Option<&PathPoint> {
        if !self.done {
            return None;
        }
        self.seq_single.get((index.max(1) - 1) as usize)
    }

    /// `Clear()` (`gxx:110-123`): reset the containers and re-add the dummy
    /// data that keeps the 1-based indices of the template.
    fn clear(&mut self) {
        self.wd1.clear();
        self.wd2.clear();
        let dummy = WalkingData {
            etat: -10,
            ustart: 0.0,
            vstart: 0.0,
        };
        self.wd1.push(dummy);
        self.wd2.push(dummy);
        self.nb_multiplicities.clear();
        self.nb_multiplicities.push(-1);

        self.done = false;
        self.seq_ajout.clear();
        self.lines.clear();
    }

    /// `Perform(Pnts1, Pnts2, Func, Caro, Reversed)` (`gxx:146-302`).
    pub fn perform(
        &mut self,
        pnts1: &[PathPoint],
        pnts2: &[InteriorPoint],
        func: &mut ContapSurfFunction,
        caro: &SurfaceAdaptor,
        reversed: bool,
    ) {
        let nb_pnts1 = pnts1.len();
        let nb_pnts2 = pnts2.len();

        self.clear();
        self.reversed = reversed;

        self.um = caro.first_u_parameter();
        self.vm = caro.first_v_parameter();
        self.umax = caro.last_u_parameter();
        self.vmax = caro.last_v_parameter();

        if self.umax < self.um {
            std::mem::swap(&mut self.umax, &mut self.um);
        }
        if self.vmax < self.vm {
            std::mem::swap(&mut self.vmax, &mut self.vm);
        }

        let step_u = self.pas * (self.umax - self.um);
        let step_v = self.pas * (self.vmax - self.vm);

        let mut umult: Vec<f64> = Vec::new();
        let mut vmult: Vec<f64> = Vec::new();

        // `gxx:186-218`: `etat1`, `ustart`, `vstart` and the multiplicities.
        for path_pnt in pnts1.iter() {
            let mut a_wd1 = WalkingData {
                etat: 1,
                ustart: 0.0,
                vstart: 0.0,
            };
            if !path_pnt.is_passing_pnt() {
                a_wd1.etat = 11;
            }
            if !path_pnt.is_tangent() {
                a_wd1.etat += 1;
            }
            if a_wd1.etat == 2 {
                a_wd1.etat = 11;
            }

            let (u, v) = path_pnt.value_2d();
            a_wd1.ustart = u;
            a_wd1.vstart = v;
            self.s_range_u.add(a_wd1.ustart);
            self.s_range_v.add(a_wd1.vstart);

            self.wd1.push(a_wd1);
            let a_nb_mult = path_pnt.multiplicity();
            self.nb_multiplicities.push(a_nb_mult);

            for j in 1..=a_nb_mult {
                if let Some((u, v)) = path_pnt.parameters(j) {
                    umult.push(u);
                    vmult.push(v);
                }
            }
        }

        // `gxx:220-234`: `etat2` of the interior points.
        for an_ip in pnts2.iter() {
            let (u, v) = an_ip.parameters();
            let mut a_wd2 = WalkingData {
                etat: 1,
                ustart: u,
                vstart: v,
            };
            self.s_range_u.add(a_wd2.ustart);
            self.s_range_v.add(a_wd2.vstart);

            if !is_tangent_ext_check(
                func,
                a_wd2.ustart,
                a_wd2.vstart,
                step_u,
                step_v,
                self.um,
                self.umax,
                self.vm,
                self.vmax,
            ) {
                a_wd2.etat = 13;
            }
            self.wd2.push(a_wd2);
        }

        let surf = caro.surface().map(|s| s.as_ref());
        self.vec_tolerance = [
            surf.map(|s| u_resolution(s, CONFUSION)).unwrap_or(0.0),
            surf.map(|s| v_resolution(s, CONFUSION)).unwrap_or(0.0),
        ];

        func.set_surface(caro.clone());

        // `gxx:240-258`.
        if self.s_range_u.delta() > self.vec_tolerance[0].max(PCONFUSION) {
            self.s_range_u.enlarge(self.s_range_u.delta());
            self.s_range_u
                .common(&BndRange::new(self.um, self.umax));
        } else {
            self.s_range_u = BndRange::new(self.um, self.umax);
        }

        if self.s_range_v.delta() > self.vec_tolerance[1].max(PCONFUSION) {
            self.s_range_v.enlarge(self.s_range_v.delta());
            self.s_range_v
                .common(&BndRange::new(self.vm, self.vmax));
        } else {
            self.s_range_v = BndRange::new(self.vm, self.vmax);
        }

        let mut rajout = false;

        // `gxx:260-267`.
        if nb_pnts1 != 0 {
            self.compute_open_line(&umult, &vmult, pnts1, func, &mut rajout);
        }
        if nb_pnts2 != 0 {
            self.compute_close_line(&umult, &vmult, pnts1, pnts2, func, &mut rajout);
        }

        // `gxx:269-292`.
        if self.to_fill_holes {
            let max_nb_iter = 10;
            let mut nb_iter = 0;
            while self.seq_alone.len() > 1 && nb_iter < max_nb_iter {
                nb_iter += 1;
                let mut pnts_in_holes: Vec<InteriorPoint> = Vec::new();
                let mut copy_seq_alone = self.seq_alone.clone();
                self.fill_pnts_in_holes(func, &mut copy_seq_alone, &mut pnts_in_holes);
                self.wd2.clear();
                self.wd2.push(WalkingData {
                    etat: -10,
                    ustart: 0.0,
                    vstart: 0.0,
                });
                for ip in pnts_in_holes.iter() {
                    let (u, v) = ip.parameters();
                    self.wd2.push(WalkingData {
                        etat: 13,
                        ustart: u,
                        vstart: v,
                    });
                }
                self.compute_close_line(&umult, &vmult, pnts1, &pnts_in_holes, func, &mut rajout);
            }
        }

        // `gxx:295-300`.
        for i in 1..=nb_pnts1 {
            if self.wd1[i].etat > 0 {
                self.seq_single.push(pnts1[i - 1].clone());
            }
        }
        self.done = true;
    }

    /// The tolerances as a `math_Vector(1, 2)` for `math_FunctionSetRoot`.
    fn tolerance_vector(&self) -> MathVector {
        let mut t = MathVector::new(1, 2);
        t.set_value(1, self.vec_tolerance[0]);
        t.set_value(2, self.vec_tolerance[1]);
        t
    }
}

/// `IsTangentExtCheck` (`gxx:40-79`): an interior point is dropped when `F`
/// stays below its tolerance in all four one-step directions.
fn is_tangent_ext_check(
    func: &mut ContapSurfFunction,
    u: f64,
    v: f64,
    step_u: f64,
    step_v: f64,
    uinf: f64,
    usup: f64,
    vinf: f64,
    vsup: f64,
) -> bool {
    use occt_math::function_set_root::MathFunctionSetWithDerivatives;

    let a_tol = func.tolerance();
    let par_u = [
        (u + step_u).min(usup),
        (u - step_u).max(uinf),
        u,
        u,
    ];
    let par_v = [v, v, (v + step_v).min(vsup), (v - step_v).max(vinf)];

    let mut a_x = MathVector::new(1, 2);
    let mut a_val = MathVector::new(1, 1);

    for i in 0..4 {
        a_x.set_value(1, par_u[i]);
        a_x.set_value(2, par_v[i]);

        if !func.value(&a_x, &mut a_val) {
            continue;
        }
        if func.root().abs() > a_tol {
            return false;
        }
    }
    true
}
