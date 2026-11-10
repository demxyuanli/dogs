//! `BRepTopAdaptor_FClass2d` — UV point classification against a face's sampled
//! boundary wires (`CSLib_Class2d` based).
//!
//! Source: `ModelingAlgorithms/TKTopAlgo/BRepTopAdaptor/BRepTopAdaptor_FClass2d.hxx/.cxx`.
//! Every boundary wire is sampled into a UV polygon and wrapped in a
//! `CSLib_Class2d`; the winding sign of each polygon becomes `TabOrien`, so a
//! point is inside the face when it is inside every positive-orientation
//! polygon and outside every negative-orientation one (`cxx:590-645`).
//!
//! Ported pieces:
//! - the constructor (`cxx:84-513`), including the per-wire sampling with
//!   `Geom2dInt_Geom2dCurveTool::NbSamples`, the unused-3D-point filter
//!   (`cxx:228-258`), the polygon-deviation measure `FlecheU`/`FlecheV`
//!   (`cxx:262-283`) and the re-discretization loop (`cxx:344-436`);
//! - `PerformInfinitePoint` (`cxx:515-523`), `Perform` (`cxx:525-683`).
//!
//! Not ported:
//! - `TestOnRestriction` (`cxx:685-841`). No caller in the ported tree uses
//!   it (only `BRepExtrema`/`BOPTools` callers of OCCT).
//! - `Copy` (`cxx:850-860`), which throws in OCCT.
//!
//! Translation notes:
//! - `BRepAdaptor_Surface::GetType()` is
//!   [`crate::geom_bnd_lib_surface3d::surface_kind`], which reports the same
//!   `GeomAbs_SurfaceType` (a trimmed surface reports its basis).
//! - `BRep_Tool::CurveOnSurface(E, F, first, last)` is
//!   [`crate::boptools_2d::curve_on_surface_range`]; `BRep_Tool::Range(E, F, ...)`
//!   in the re-discretization loop reads the same representation range.
//! - `GCPnts_QuasiUniformDeflection(C, defl)` on the 2D adaptor is
//!   [`QuasiUniformDeflection2d`] (the `Adaptor2d_Curve2d` instantiation).
//! - `BRepClass_FaceClassifier` is [`FClassifier`] over [`FaceExplorer`], the
//!   faithful `BRepClass_FClassifier` port (the same class
//!   `BRepClass_FaceClassifier::Perform` builds).

use std::sync::Arc;

use occt_core::cslib::{Class2d, Class2dResult};
use occt_core::elib::clib2d::parameter_lin2d;
use occt_core::elib::clib::line2d_value;
use occt_core::gp::{GpDir2d, GpLin2d, GpPnt, GpPnt2d, GpVec2d};
use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom2d::gcpnts_quasi_uniform_deflection::QuasiUniformDeflection2d;

use crate::abs::Orientation;
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_class::{FClassifier, FaceExplorer};
use crate::brep_tool::BRepTool;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::fclass2d::FaceState;
use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};
use crate::shape::{Edge, Face, Wire};
use crate::topo_tools::{edge_vertices, edges_of_wire, wires_of_face};

/// `RealLast()`.
const REAL_LAST: f64 = 1.7e308;

/// `constexpr double eps = 1.e-10` (`cxx:110`, `cxx:439`).
const EPS: f64 = 1.0e-10;

/// `constexpr double inv_eps = 1.e-9` (`cxx:146`, `cxx:363`): the parametric
/// range below which an edge is skipped.
const SKIP_RANGE: f64 = 1.0e-9;

/// `safeIncrement` (`cxx:49-60`): `theValue + theIncrement`, or the next
/// representable double towards `theDirection` when the sum does not move.
fn safe_increment(value: f64, direction: f64, increment: f64) -> f64 {
    let next = value + increment;
    if next == value {
        next_after(value, direction)
    } else {
        next
    }
}

/// `std::nextafter`.
fn next_after(from: f64, to: f64) -> f64 {
    if from == to {
        return to;
    }
    let bits = from.to_bits();
    let next = if (to > from) == (from > 0.0) {
        bits + 1
    } else {
        bits - 1
    };
    f64::from_bits(next)
}

/// `isDegenerated` (`cxx:58-82`): every sampled point lies within
/// `Precision::Confusion()` of the start point.
fn is_degenerated(curve: &Arc<dyn occt_geom::curve::Curve>, start: f64, end: f64) -> bool {
    let step = (end - start) * 0.1;
    let first = curve.d0(start);
    let mut param = start;
    while param < end {
        let p = curve.d0(param);
        if first.square_distance(&p) > CONFUSION {
            return false;
        }
        param = safe_increment(param, end, step);
    }
    true
}

/// Build `gp_Dir2d(gp_Vec2d(from, to))` (`cxx:266`).
fn dir_between(from: &GpPnt2d, to: &GpPnt2d) -> Option<GpDir2d> {
    let v = GpVec2d::new(to.x() - from.x(), to.y() - from.y());
    GpDir2d::from_vec2d(&v).ok()
}

/// The `FlecheU` / `FlecheV` update shared by the sampling loop
/// (`cxx:262-283`) and the re-discretization loop (`cxx:392-408`): the deviation
/// of the middle sampled point from the chord of the last two points.
fn update_fleche(
    seq: &[GpPnt2d],
    fleche_u: &mut f64,
    fleche_v: &mut f64,
) {
    let len = seq.len();
    if len < 3 {
        return;
    }
    let (p0, p1, p2) = (seq[len - 3], seq[len - 2], seq[len - 1]);
    let Some(dir) = dir_between(&p0, &p2) else {
        return;
    };
    let lin = GpLin2d::from_pnt_dir(p0, dir);
    let ul = parameter_lin2d(&lin, &p1);
    let pp = line2d_value(&lin, ul);
    let du = (pp.x() - p1.x()).abs();
    let dv = (pp.y() - p1.y()).abs();
    if du > *fleche_u {
        *fleche_u = du;
    }
    if dv > *fleche_v {
        *fleche_v = dv;
    }
}

/// `BRepClass_FaceClassifier::Perform(F, Puv, Tol)` → `State()`.
fn classify_point(face: &Face, puv: &GpPnt2d, tol: f64) -> FaceState {
    let mut fexp = FaceExplorer::new(face);
    FClassifier::with_args(&mut fexp, puv, tol).state()
}

/// Re-discretization of one wire at `defl` (`BRepTopAdaptor_FClass2d.cxx:349-411`).
///
/// Returns the sampled points plus the running `(FlecheU, FlecheV)`. A wire
/// whose discretization fails (`aDiscr.IsDone() == false`) stops the walk and
/// keeps the points gathered so far, which is the `break` at `cxx:373`.
fn rediscretize_wire(
    wire: &Wire,
    face: &Face,
    defl: f64,
) -> (Vec<GpPnt2d>, f64, f64) {
    let mut seq: Vec<GpPnt2d> = Vec::new();
    let mut fleche_u = 0.0;
    let mut fleche_v = 0.0;
    let mut firstpoint = 1;

    let walk: Vec<(Edge, Option<Orientation>)> = {
        let mut we = WireExplorer::of_wire_on_face(wire, face);
        let mut out = Vec::new();
        while we.more() {
            out.push((we.current_edge().unwrap_or_else(Edge::new), we.orientation().ok()));
            we.next();
        }
        out
    };

    for (edge, or) in walk {
        let Some(or) = or else { continue };
        if or == Orientation::Internal || or == Orientation::External {
            continue;
        }
        let Some((_pc, pfbid, plbid)) = curve_on_surface_range(&edge, face) else {
            continue;
        };
        if (plbid - pfbid).abs() < SKIP_RANGE {
            continue;
        }
        let adaptor = Curve2dAdaptor::of_edge_face(&edge, face);
        let Some(pc) = adaptor.curve() else { continue };

        let mut discr = QuasiUniformDeflection2d::new();
        // `GCPnts_QuasiUniformDeflection(C, aDiscrDefl)`: the adaptor's own
        // `[First, Last]` are the `CurveOnSurface` range (`pfbid`, `plbid`).
        discr.initialize(pc.as_ref(), defl, pfbid, plbid, 1);
        if !discr.is_done() {
            // `cxx:371-374`: `break` keeps the points gathered so far.
            break;
        }
        let nbp = discr.nb_points();
        // `cxx:375-383`.
        let (mut i, i_step, i_end) = if or == Orientation::Reversed {
            (nbp, -1, 0)
        } else {
            (1, 1, nbp + 1)
        };
        if firstpoint == 2 {
            i += i_step;
        }
        while i != i_end {
            if let Some(param) = discr.parameter(i) {
                // `C.Value(aDiscr.Parameter(i))`.
                if let Some(p) = adaptor.d0(param) {
                    seq.push(p);
                }
            }
            i += i_step;
        }
        if nbp > 2 {
            update_fleche(&seq, &mut fleche_u, &mut fleche_v);
        }
        firstpoint = 2;
    }

    (seq, fleche_u, fleche_v)
}

/// `BRepTopAdaptor_FClass2d`.
pub struct BrepTopAdaptorFClass2d {
    /// `myToluv`
    toluv: f64,
    /// `myFace` — the constructor's `Face` copy, oriented `TopAbs_FORWARD`
    /// (`cxx:97`).
    face: Face,
    /// `myTabClass`
    tab_class: Vec<Class2d>,
    /// `myTabOrien`: 1 = area > 0, 0 = area < 0 (a hole), -1 = bad wire.
    tab_orien: Vec<i32>,
    /// `myUmin` / `myVmin` / `myUmax` / `myVmax`
    umin: f64,
    vmin: f64,
    umax: f64,
    vmax: f64,
    /// `myU1` / `myV1` / `myU2` / `myV2`: the surface's UV period window.
    u1: f64,
    v1: f64,
    u2: f64,
    v2: f64,
}

impl BrepTopAdaptorFClass2d {
    /// `BRepTopAdaptor_FClass2d(aFace, TolUV)` (`cxx:84-513`).
    pub fn new(face: &Face, toluv: f64) -> Self {
        let wires = wires_of_face(face);
        Self::new_with_wires(face, &wires, toluv)
    }

    /// The constructor over an explicit wire list.
    ///
    /// OCCT builds a throwaway face with `BRep_Builder` for
    /// `BRepCheck_Face::ClassifyWires` (`BRepCheck_Face.cxx:355-360`):
    /// `myShape.EmptyCopied()` (which copies the surface, `BRep_TFace.cxx:36-43`)
    /// plus one wire. `BRep_Tool::CurveOnSurface` matches a representation by
    /// *surface*, so that throwaway face resolves to the same pcurves as the
    /// original; this port keys pcurves by face, so the same object is expressed
    /// as "this face, only these wires".
    pub fn new_with_wires(face: &Face, wires: &[Wire], toluv: f64) -> Self {
        let face = Face(face.0.oriented(Orientation::Forward));
        let mut this = Self {
            toluv,
            face,
            tab_class: Vec::new(),
            tab_orien: Vec::new(),
            umin: 0.0,
            vmin: 0.0,
            umax: 0.0,
            vmax: 0.0,
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
        };

        let mut a_nb_e;
        let mut is_bad_wire = false;

        for wire in wires {
            if is_bad_wire {
                break;
            }
            let mut nbpnts = 0i32;
            let mut seq: Vec<GpPnt2d> = Vec::new();
            let mut firstpoint = 1i32;
            let mut fleche_u = 0.0;
            let mut fleche_v = 0.0;
            let mut wire_is_not_empty = false;

            let nb_edges_total = edges_of_wire(wire).len() as i32;
            let mut nb_edges = nb_edges_total;
            a_nb_e = nb_edges_total;

            let mut ancien_pnt3d = GpPnt::new(0.0, 0.0, 0.0);
            let mut ancien_pnt3d_initialised = false;

            let walk: Vec<(Edge, Option<Orientation>)> = {
                let mut we = WireExplorer::of_wire_on_face(wire, &this.face);
                let mut out = Vec::new();
                while we.more() {
                    out.push((we.current_edge().unwrap_or_else(Edge::new), we.orientation().ok()));
                    we.next();
                }
                out
            };

            let mut abort = false;
            for (edge, or) in walk {
                nb_edges -= 1;
                let Some(or) = or else { continue };
                if or == Orientation::Internal || or == Orientation::External {
                    continue;
                }

                let Some((pc, pfbid, plbid)) = curve_on_surface_range(&edge, &this.face) else {
                    // `cxx:142-145`: `return;` leaves `TabClass` empty, so every
                    // later `Perform` answers IN (`cxx:530-533`).
                    abort = true;
                    break;
                };

                if (plbid - pfbid).abs() < SKIP_RANGE {
                    continue;
                }

                let mut degenerated = BRepTool::is_degenerated(&edge)
                    || BRepTool::is_closed_edge_face(&edge, &this.face);
                let (va, vb) = edge_vertices(&edge);
                if va.is_none() || vb.is_none() {
                    degenerated = true;
                }

                let curve3d = BRepTool::edge_curve(&edge);
                if !degenerated {
                    if let Some(c3) = curve3d.as_ref() {
                        // `cxx:171-176`: `isDegenerated(aCurveAdaptor3D, pfbid, plbid)`.
                        degenerated = is_degenerated(c3, pfbid, plbid);
                    }
                }

                // `cxx:180-196`: `NbSamples(aCurveAdaptor2D)`, quadrupled above 2.
                let mut nbs = crate::curve_sampling_2d::geom2d_nb_samples(pc.as_ref(), pfbid, plbid) as i32;
                if nbs > 2 {
                    nbs *= 4;
                }
                let mut du = (plbid - pfbid) / (nbs - 1) as f64;
                let mut u = if or == Orientation::Forward { pfbid } else { plbid };
                if or != Orientation::Forward {
                    du = -du;
                }

                if firstpoint == 2 {
                    u += du;
                }
                let avant = nbpnts;
                for _e in firstpoint..=nbs {
                    let p2d = pc.d0(u);
                    if p2d.x() < this.umin {
                        this.umin = p2d.x();
                    }
                    if p2d.x() > this.umax {
                        this.umax = p2d.x();
                    }
                    if p2d.y() < this.vmin {
                        this.vmin = p2d.y();
                    }
                    if p2d.y() > this.vmax {
                        this.vmax = p2d.y();
                    }

                    let mut dist = 1.0e20;
                    let mut p3d = GpPnt::new(0.0, 0.0, 0.0);
                    if !degenerated {
                        if let Some(c3) = curve3d.as_ref() {
                            p3d = c3.d0(u);
                            if nbpnts > 1 && ancien_pnt3d_initialised {
                                dist = p3d.distance(&ancien_pnt3d);
                            }
                        }
                    }

                    let mut is_real_curve3d = true;
                    if dist < CONFUSION {
                        if let Some(c3) = curve3d.as_ref() {
                            let mid_p3d = c3.d0(u - du / 2.0);
                            if p3d.distance(&mid_p3d) < CONFUSION {
                                is_real_curve3d = false;
                            }
                        }
                    }

                    if is_real_curve3d {
                        if !degenerated {
                            ancien_pnt3d = p3d;
                            ancien_pnt3d_initialised = true;
                        }
                        nbpnts += 1;
                        seq.push(p2d);
                    }

                    u += du;
                    let ii = nbpnts;
                    if ii > avant + 4 {
                        let a = seq[(ii - 2 - 1) as usize];
                        let b = seq[(ii - 1) as usize];
                        if a.square_distance(&b) != 0.0 {
                            update_fleche(&seq, &mut fleche_u, &mut fleche_v);
                        }
                    }
                }

                if firstpoint == 1 {
                    firstpoint = 2;
                }
                wire_is_not_empty = true;
            }

            if abort {
                // The constructor returns with whatever `TabClass` it had.
                return this;
            }

            if nb_edges != 0 {
                // `cxx:290-301`: leftover edges (a bad wire) get a degenerate
                // 2-point classifier and mark the whole face unclassifiable.
                let p = GpPnt2d::new(0.0, 0.0);
                this.tab_class.push(Class2d::new(
                    &[p, p],
                    fleche_u,
                    fleche_v,
                    this.umin,
                    this.vmin,
                    this.umax,
                    this.vmax,
                ));
                is_bad_wire = true;
                this.tab_orien.push(-1);
            } else if wire_is_not_empty {
                if nbpnts > 3 {
                    let mut square = 0.0;
                    let mut a_per = 0.0;
                    let mut im0 = 1i32;
                    let mut im1 = nbpnts - 1;
                    for _ii in 1..nbpnts {
                        if im1 >= nbpnts {
                            im1 = 1;
                        }
                        let p0 = seq[(im0 - 1) as usize];
                        let p1 = seq[(im1 - 1) as usize];
                        square += (p0.x() - p1.x()) * (p0.y() + p1.y()) * 0.5;
                        a_per += GpVec2d::new(p0.x() - p1.x(), p0.y() - p1.y()).magnitude();
                        im0 += 1;
                        im1 += 1;
                    }

                    let mut an_exp_thick = (2.0 * square.abs() / a_per).max(1e-7);
                    let mut a_defl = fleche_u.max(fleche_v);
                    let mut a_discr_defl = (a_defl * 0.1).min(an_exp_thick * 10.0);

                    while a_defl > an_exp_thick && a_discr_defl > 1e-7 {
                        // `cxx:349-436`: the polygon's deviation is too large for
                        // its area/perimeter ratio; re-discretize more tightly.
                        let (new_seq, new_fu, new_fv) =
                            rediscretize_wire(wire, &this.face, a_discr_defl);
                        seq = new_seq;
                        fleche_u = new_fu;
                        fleche_v = new_fv;
                        nbpnts = seq.len() as i32;
                        if nbpnts < 2 {
                            break;
                        }
                        square = 0.0;
                        a_per = 0.0;
                        im0 = 1;
                        im1 = nbpnts - 1;
                        for _ii in 1..nbpnts {
                            if im1 >= nbpnts {
                                im1 = 1;
                            }
                            let p0 = seq[(im0 - 1) as usize];
                            let p1 = seq[(im1 - 1) as usize];
                            square += (p0.x() - p1.x()) * (p0.y() + p1.y()) * 0.5;
                            a_per += GpVec2d::new(p0.x() - p1.x(), p0.y() - p1.y()).magnitude();
                            im0 += 1;
                            im1 += 1;
                        }
                        if a_per <= 0.0 {
                            break;
                        }
                        an_exp_thick = (2.0 * square.abs() / a_per).max(1e-7);
                        a_defl = fleche_u.max(fleche_v);
                        a_discr_defl = (a_discr_defl * 0.1).min(an_exp_thick * 10.0);
                    }

                    let orien = if a_nb_e == 1 && fleche_u < EPS && fleche_v < EPS && square.abs() < EPS
                    {
                        1
                    } else if square < 0.0 {
                        1
                    } else {
                        0
                    };
                    this.tab_orien.push(orien);

                    if fleche_u < this.toluv {
                        fleche_u = this.toluv;
                    }
                    if fleche_v < this.toluv {
                        fleche_v = this.toluv;
                    }
                    this.tab_class.push(Class2d::new(
                        &seq[..nbpnts as usize],
                        fleche_u,
                        fleche_v,
                        this.umin,
                        this.vmin,
                        this.umax,
                        this.vmax,
                    ));
                } else {
                    // `cxx:458-466`.
                    is_bad_wire = true;
                    this.tab_orien.push(-1);
                    let p1 = seq.first().copied().unwrap_or_else(|| GpPnt2d::new(0.0, 0.0));
                    let p2 = seq.get(1).copied().unwrap_or(p1);
                    this.tab_class.push(Class2d::new(
                        &[p1, p2],
                        fleche_u,
                        fleche_v,
                        this.umin,
                        this.vmin,
                        this.umax,
                        this.vmax,
                    ));
                }
            }
        }

        if !this.tab_class.is_empty() {
            if is_bad_wire {
                // `cxx:475-478`.
                this.tab_orien[0] = -1;
            }

            let surf = BRepTool::face_surface(&this.face);
            let kind = surf
                .as_ref()
                .map(|s| surface_kind(s.as_ref()))
                .unwrap_or(GeomAbsSurfaceType::OtherSurface);

            // `cxx:480-512`: surfaces periodic in `u` get the full `[0, 2pi]`
            // window; a torus is periodic in `v` as well.
            if matches!(
                kind,
                GeomAbsSurfaceType::Cone
                    | GeomAbsSurfaceType::Cylinder
                    | GeomAbsSurfaceType::Torus
                    | GeomAbsSurfaceType::Sphere
                    | GeomAbsSurfaceType::SurfaceOfRevolution
            ) {
                this.u1 = 0.0;
                this.u2 = 2.0 * std::f64::consts::PI;
            } else {
                this.u1 = 0.0;
                this.u2 = 0.0;
            }

            if kind == GeomAbsSurfaceType::Torus {
                this.v1 = 0.0;
                this.v2 = 2.0 * std::f64::consts::PI;
            } else {
                this.v1 = 0.0;
                this.v2 = 0.0;
            }
        }

        this
    }

    /// `PerformInfinitePoint()` (`cxx:515-523`).
    pub fn perform_infinite_point(&self) -> FaceState {
        if self.umax == -REAL_LAST
            || self.vmax == -REAL_LAST
            || self.umin == REAL_LAST
            || self.vmin == REAL_LAST
        {
            return FaceState::In;
        }
        let p = GpPnt2d::new(
            self.umin - (self.umax - self.umin),
            self.vmin - (self.vmax - self.vmin),
        );
        self.perform(p, false)
    }

    /// `Perform(Puv, RecadreOnPeriodic)` (`cxx:525-683`).
    pub fn perform(&self, puv: GpPnt2d, recadre_on_periodic: bool) -> FaceState {
        let nbtabclass = self.tab_class.len();
        if nbtabclass == 0 {
            return FaceState::In;
        }

        let mut u = puv.x();
        let mut v = puv.y();
        let mut uu = u;
        let mut vv = v;

        let surf = BRepTool::face_surface(&self.face);
        let (is_u_per, is_v_per, uperiod, vperiod) = match surf.as_ref() {
            Some(s) => (
                s.is_u_periodic(),
                s.is_v_periodic(),
                if s.is_u_periodic() { s.u_period() } else { 0.0 },
                if s.is_v_periodic() { s.v_period() } else { 0.0 },
            ),
            None => (false, false, 0.0, 0.0),
        };

        let mut a_status = FaceState::Unknown;
        let mut urecadre = false;
        let mut vrecadre = false;

        if recadre_on_periodic {
            if is_u_per {
                if uu < self.umin {
                    while uu < self.umin {
                        uu += uperiod;
                    }
                } else {
                    while uu >= self.umin {
                        uu -= uperiod;
                    }
                    uu += uperiod;
                }
            }
            if is_v_per {
                if vv < self.vmin {
                    while vv < self.vmin {
                        vv += vperiod;
                    }
                } else {
                    while vv >= self.vmin {
                        vv -= vperiod;
                    }
                    vv += vperiod;
                }
            }
        }

        loop {
            let mut dedans = 1;
            let p = GpPnt2d::new(u, v);

            if self.tab_orien[0] != -1 {
                for n in 0..nbtabclass {
                    let cur = self.tab_class[n].si_dans(&p);
                    match cur {
                        Class2dResult::Inside => {
                            if self.tab_orien[n] == 0 {
                                dedans = -1;
                                break;
                            }
                        }
                        Class2dResult::Outside => {
                            if self.tab_orien[n] == 1 {
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
                if dedans == 0 {
                    // `cxx:622-629`: `m_Toluv = min(Toluv, 4.0)`.
                    let m_toluv = if self.toluv > 4.0 { 4.0 } else { self.toluv };
                    a_status = classify_point(&self.face, &p, m_toluv);
                }
                if dedans == 1 {
                    a_status = FaceState::In;
                }
                if dedans == -1 {
                    a_status = FaceState::Out;
                }
            } else {
                // `TabOrien(1) == -1`: a bad wire, classify with the general
                // face classifier (`cxx:639-644`).
                a_status = classify_point(&self.face, &p, self.toluv);
            }

            if !recadre_on_periodic || (!is_u_per && !is_v_per) {
                return a_status;
            }
            if a_status == FaceState::In || a_status == FaceState::On {
                return a_status;
            }

            if !urecadre {
                u = uu;
                urecadre = true;
            } else if is_u_per {
                u += uperiod;
            }
            if u > self.umax || !is_u_per {
                if !vrecadre {
                    v = vv;
                    vrecadre = true;
                } else if is_v_per {
                    v += vperiod;
                }

                u = uu;

                if v > self.vmax || !is_v_per {
                    return a_status;
                }
            }
        }
    }

    /// `Toluv()`.
    pub fn toluv(&self) -> f64 {
        self.toluv
    }

    /// `U1()` / `U2()` / `V1()` / `V2()` — the periodic window (`cxx:490-511`).
    pub fn period_window(&self) -> (f64, f64, f64, f64) {
        (self.u1, self.v1, self.u2, self.v2)
    }

    /// `Umin()` / `Umax()` / `Vmin()` / `Vmax()` — the sampled box.
    pub fn uv_box(&self) -> (f64, f64, f64, f64) {
        (self.umin, self.umax, self.vmin, self.vmax)
    }
}

/// `Precision::PConfusion()` re-exported for [`BRepCheck_Face::IsInside`]'s
/// small-parametric-range skip.
///
/// [`BRepCheck_Face::IsInside`]: crate::brep_check
pub fn p_confusion() -> f64 {
    PCONFUSION
}
