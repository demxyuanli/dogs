//! Port of ShapeFix_ComposeShell (the class itself).
//! ShapeFix_ComposeShell.cxx:79-205 (ctor + Init), 647-832 (ComputeCode).
//! T-97: items below are faithful ports of the named OCCT source, but their
//! OCCT-side consumers are not all ported yet, so parts are not called from this
//! crate. The `dead_code` allowance is deliberate: **pending wiring**, not dead
//! code. Do not delete them to silence warnings (see
//! specs/_a3n00_gap_analysis.md §9.309/§9.310); wire the consumer instead.
#![allow(dead_code)]

use occt_core::gp::{GpLin2d, GpPnt2d};
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};

use super::composite_surface::CompositeSurface;
use super::helpers::*;
use super::reshape::{MapReShape, SharedReShape};

/// ShapeExtend_Status bits: ShapeExtend::EncodeStatus (ShapeExtend.cxx:54-98).
pub const SHAPEEXTEND_OK: i32 = 0x0000;
pub const SHAPEEXTEND_DONE1: i32 = 0x0001;
pub const SHAPEEXTEND_DONE2: i32 = 0x0002;
pub const SHAPEEXTEND_DONE3: i32 = 0x0004;
pub const SHAPEEXTEND_DONE4: i32 = 0x0008;
pub const SHAPEEXTEND_DONE5: i32 = 0x0010;
pub const SHAPEEXTEND_DONE6: i32 = 0x0020;
pub const SHAPEEXTEND_DONE7: i32 = 0x0040;
pub const SHAPEEXTEND_DONE8: i32 = 0x0080;
pub const SHAPEEXTEND_DONE: i32 = 0x00ff;
pub const SHAPEEXTEND_FAIL1: i32 = 0x0100;
pub const SHAPEEXTEND_FAIL2: i32 = 0x0200;
pub const SHAPEEXTEND_FAIL3: i32 = 0x0400;
pub const SHAPEEXTEND_FAIL4: i32 = 0x0800;
pub const SHAPEEXTEND_FAIL5: i32 = 0x1000;
pub const SHAPEEXTEND_FAIL6: i32 = 0x2000;
pub const SHAPEEXTEND_FAIL7: i32 = 0x4000;
pub const SHAPEEXTEND_FAIL8: i32 = 0x8000;
pub const SHAPEEXTEND_FAIL: i32 = 0xff00;

/// ShapeExtend::DecodeStatus (ShapeExtend.cxx:102-109).
pub fn decode_status(flag: i32, status: i32) -> bool {
    if status == SHAPEEXTEND_OK {
        return flag == 0;
    }
    (flag & status) != 0
}

/// ShapeFix_ComposeShell (ShapeFix_ComposeShell.hxx).
pub struct ComposeShell {
    pub(super) face: Option<Face>,
    pub(super) grid: CompositeSurface,
    pub(super) u_closed: bool,
    pub(super) v_closed: bool,
    pub(super) u_period: f64,
    pub(super) v_period: f64,
    pub(super) u_resolution: f64,
    pub(super) v_resolution: f64,
    pub(super) closed_mode: bool,
    pub(super) status: i32,
    pub(super) orient: Orientation,
    pub(super) precision: f64,
    pub(super) min_tolerance: f64,
    pub(super) max_tolerance: f64,
    /// myResult.
    pub(super) result: Option<TopoShape>,
    /// myContext: the ShapeBuild_ReShape shared by every SplitWire call
    /// (`ShapeFix_Root::Context()`).
    pub(super) context: SharedReShape,
    /// myInvertEdgeStatus (cxx:84 sets true in the ctor; Perform resets it at
    /// cxx:209 and CollectWires raises it at cxx:2731; MakeFacesOnPatch reads
    /// it at cxx:2997).
    pub(super) invert_edge_status: bool,
}

impl Default for ComposeShell {
    fn default() -> Self {
        Self::new()
    }
}

impl ComposeShell {
    /// ShapeFix_ComposeShell() (cxx:79-94); tolerance bounds come from
    /// ShapeFix_Root defaults (MinTolerance = Precision::Confusion(),
    /// MaxTolerance = FromSTEP.FixShape.MaxTolerance3d = 1.0, see
    /// shhealing/wire_fix.rs:26-35 for that constant provenance).
    pub fn new() -> Self {
        Self {
            face: None,
            grid: CompositeSurface::new(),
            u_closed: false,
            v_closed: false,
            u_period: 0.0,
            v_period: 0.0,
            u_resolution: 0.0,
            v_resolution: 0.0,
            closed_mode: false,
            status: SHAPEEXTEND_OK,
            orient: Orientation::Forward,
            precision: CONFUSION,
            min_tolerance: CONFUSION,
            max_tolerance: 1.0,
            result: None,
            context: SharedReShape::new(),
            invert_edge_status: true,
        }
    }

    /// ShapeFix_Root::LimitTolerance (ShapeFix_Root.hxx:80-81, "Returns
    /// tolerance limited by [myMinTol,myMaxTol]").
    pub fn limit_tolerance(&self, toler: f64) -> f64 {
        toler.max(self.min_tolerance).min(self.max_tolerance)
    }

    pub fn set_max_tolerance(&mut self, t: f64) {
        self.max_tolerance = t;
    }

    pub fn max_tolerance(&self) -> f64 {
        self.max_tolerance
    }

    /// CompShell.ClosedMode() (cxx:2258).
    pub fn set_closed_mode(&mut self, on: bool) {
        self.closed_mode = on;
    }

    pub fn closed_mode(&self) -> bool {
        self.closed_mode
    }

    pub fn status(&self) -> i32 {
        self.status
    }

    /// ShapeFix_ComposeShell::Result() (cxx:274-277).
    pub fn result(&self) -> Option<&TopoShape> {
        self.result.as_ref()
    }

    /// myResult assignment (cxx:245/249).
    pub fn set_result(&mut self, r: TopoShape) {
        self.result = Some(r);
    }

    /// ShapeFix_ComposeShell::Status(status) (cxx:281-284).
    pub fn status_of(&self, status: i32) -> bool {
        decode_status(self.status, status)
    }

    pub fn face(&self) -> Option<&Face> {
        self.face.as_ref()
    }

    /// `myInvertEdgeStatus` (read by `MakeFacesOnPatch` cxx:2997).
    pub fn invert_edge_status(&self) -> bool {
        self.invert_edge_status
    }

    /// `ShapeFix_Root::Context()`.
    pub fn context(&self) -> &SharedReShape {
        &self.context
    }

    pub fn context_mut(&mut self) -> &mut SharedReShape {
        &mut self.context
    }

    /// `ShapeFix_Root::SetContext` (`ShapeFix_Root.cxx:...`).
    pub fn set_context(&mut self, ctx: SharedReShape) {
        self.context = ctx;
    }

    /// Init(Grid, L, Face, Prec) (cxx:96-202). `L` (a TopLoc_Location) is not
    /// kept: the port's shapes carry no location.
    pub fn init(&mut self, grid: CompositeSurface, face: &Face, prec: f64) {
        self.u_closed = grid.is_u_closed();
        self.v_closed = grid.is_v_closed();
        self.u_period = grid.u_joint_value(grid.nb_u_patches() + 1) - grid.u_joint_value(1);
        self.v_period = grid.v_joint_value(grid.nb_v_patches() + 1) - grid.v_joint_value(1);
        if let Some(the_surface) = BRepTool::face_surface(face) {
            let elementary = the_surface.gp_pln().is_some()
                || the_surface.gp_cylinder().is_some()
                || the_surface.gp_cone().is_some()
                || the_surface.gp_sphere().is_some()
                || the_surface.gp_torus().is_some();
            if elementary {
                self.u_closed = self.u_closed && the_surface.is_u_closed();
                self.v_closed = self.v_closed && the_surface.is_v_closed();
            } else {
                let (mut u0, mut u1, mut v0, mut v1) = (
                    the_surface.u_range().0,
                    the_surface.u_range().1,
                    the_surface.v_range().0,
                    the_surface.v_range().1,
                );
                // cxx:120-124: BRepTools::UVBounds(Face) supplies the bounds
                // used below whenever a surface range is infinite.
                let (gu0, gu1, gv0, gv1) = if u0.is_infinite()
                    || u1.is_infinite()
                    || v0.is_infinite()
                    || v1.is_infinite()
                {
                    BRepTool::uv_bounds(face)
                } else {
                    (0.0, 0.0, 0.0, 0.0)
                };
                if self.u_closed {
                    if v0.is_infinite() {
                        v0 = gv0;
                    }
                    if v1.is_infinite() {
                        v1 = gv1;
                    }
                }
                if self.v_closed {
                    if u0.is_infinite() {
                        u0 = gu0;
                    }
                    if u1.is_infinite() {
                        u1 = gu1;
                    }
                }
                if self.u_closed {
                    let p0 = the_surface.value(u0, 0.5 * (v0 + v1));
                    let p1 = the_surface.value(u1, 0.5 * (v0 + v1));
                    if p0.distance(&p1) > CONFUSION * 10.0 {
                        self.u_closed = false;
                    }
                }
                if self.v_closed {
                    let p0 = the_surface.value(0.5 * (u0 + u1), v0);
                    let p1 = the_surface.value(0.5 * (u0 + u1), v1);
                    if p0.distance(&p1) > CONFUSION * 10.0 {
                        self.v_closed = false;
                    }
                }
            }
        }
        self.orient = face.0.orientation();
        self.precision = prec;
        self.status = SHAPEEXTEND_OK;
        self.grid = grid;
        self.face = Some(face.clone());

        // cxx:170-193: the UV resolution is the smallest per-patch value
        // GeomAdaptor_Surface resolution scaled by the joint/patch range ratio.
        let mut ures_min = f64::INFINITY;
        let mut vres_min = f64::INFINITY;
        for i in 1..=self.grid.nb_u_patches() {
            let u_range = self.grid.u_joint_value(i + 1) - self.grid.u_joint_value(i);
            for j in 1..=self.grid.nb_v_patches() {
                let v_range = self.grid.v_joint_value(j + 1) - self.grid.v_joint_value(j);
                let Some(patch) = self.grid.patch(i, j) else {
                    continue;
                };
                let (u1, u2, v1, v2) = (
                    patch.u_range().0,
                    patch.u_range().1,
                    patch.v_range().0,
                    patch.v_range().1,
                );
                let ures = occt_geom::approx_same_parameter::u_resolution(patch.as_ref(), 1.0)
                    * u_range
                    / (u2 - u1);
                let vres = occt_geom::approx_same_parameter::v_resolution(patch.as_ref(), 1.0)
                    * v_range
                    / (v2 - v1);
                if ures > 0.0 && ures_min > ures {
                    ures_min = ures;
                }
                if vres > 0.0 && vres_min > vres {
                    vres_min = vres;
                }
            }
        }
        self.u_resolution = if ures_min == f64::INFINITY {
            // Precision::Parametric(1.) = 1. * 0.01 (Precision.hxx:328).
            1.0 * 0.01
        } else {
            ures_min
        };
        self.v_resolution = if vres_min == f64::INFINITY {
            // Precision::Parametric(1.) = 1. * 0.01 (Precision.hxx:328).
            1.0 * 0.01
        } else {
            vres_min
        };
    }

    /// ComputeCode (cxx:647-832): IOR side code of the wire segment between
    /// two intersections, by deviation of its pcurve from `line` measured at
    /// NPOINTS points. `edges` is the wire data (1-based `beg_ind`/`end_ind`).
    pub fn compute_code(
        &mut self,
        edges: &[Edge],
        line: &GpLin2d,
        beg_ind: usize,
        end_ind: usize,
        beg_par: f64,
        end_par: f64,
        is_internal: bool,
    ) -> i32 {
        const NPOINTS: usize = 5;
        let mut code = IOR_UNDEF;
        let face = match &self.face {
            Some(f) => f.clone(),
            None => return code,
        };
        let mut special: i32 = if beg_ind == end_ind {
            let o = edges[beg_ind - 1].0.orientation();
            let fwd = o == Orientation::Forward || o == Orientation::Internal;
            if fwd == (beg_par > end_par) { 1 } else { 0 }
        } else {
            0
        };
        if special == 0 && beg_ind == end_ind && beg_par == end_par && (self.closed_mode || is_internal) {
            special = 1;
        }
        let mut begin = true;
        let mut shift = 0.0f64;
        let mut p2d0 = GpPnt2d::new(0.0, 0.0);
        let nb = edges.len();
        let mut i = beg_ind;
        loop {
            if i > nb {
                i = 1;
            }
            let edge = &edges[i - 1];
            let pcurve = crate::boptools_2d::curve_on_surface_oriented(edge, &face, false);
            let Some((c2d, f, l)) = pcurve else {
                self.status |= SHAPEEXTEND_FAIL3;
                i += 1;
                continue;
            };
            let tol = self.limit_tolerance(BRepTool::edge_tolerance(edge));
            let isreversed = edge.0.orientation() == Orientation::Reversed;
            let par1 = if i == beg_ind && special >= 0 {
                beg_par
            } else if isreversed {
                l
            } else {
                f
            };
            let par2 = if i == end_ind && special <= 0 {
                end_par
            } else if isreversed {
                f
            } else {
                l
            };
            let dpar = (par2 - par1) / (NPOINTS as f64 - 1.0);
            let np = if dpar.abs() < PCONFUSION { 1 } else { NPOINTS };
            let mut j = 0usize;
            while j < np {
                let par = par1 + dpar * j as f64;
                let mut p2d = c2d.d0(par);
                if self.closed_mode {
                    if self.u_closed && line.direction().x().abs() < PCONFUSION {
                        if begin {
                            shift = crate::shhealing::adjust_by_period(
                                p2d.x(),
                                line.location().x(),
                                self.u_period,
                            );
                        } else if j == 0 {
                            shift = crate::shhealing::adjust_by_period(
                                p2d.x() - p2d0.x(),
                                0.0,
                                self.u_period,
                            );
                        }
                        p2d = GpPnt2d::new(p2d.x() + shift, p2d.y());
                    }
                    if self.v_closed && line.direction().y().abs() < PCONFUSION {
                        if begin {
                            shift = crate::shhealing::adjust_by_period(
                                p2d.y(),
                                line.location().y(),
                                self.v_period,
                            );
                        } else if j == 0 {
                            shift = crate::shhealing::adjust_by_period(
                                p2d.y() - p2d0.y(),
                                0.0,
                                self.v_period,
                            );
                        }
                        p2d = GpPnt2d::new(p2d.x(), p2d.y() + shift);
                    }
                    begin = false;
                }
                p2d0 = p2d;
                let (pos, _dev) = point_line_position(&p2d, line);
                if pos != IOR_UNDEF {
                    let p2dl = project_point_on_line(&p2d, line);
                    if !is_coincided(&p2d, &p2dl, self.u_resolution, self.v_resolution, tol) {
                        if !self.closed_mode {
                            code = pos;
                            break;
                        } else {
                            code |= pos;
                        }
                    }
                }
                j += 1;
            }
            if j < np {
                // not tangency
                i = 0;
                break;
            }
            if i == end_ind {
                if special <= 0 {
                    break;
                } else {
                    special = -1;
                }
            }
            i += 1;
        }
        if self.closed_mode {
            if code != IOR_UNDEF && !begin {
                // in closed mode, if segment is of 2*pi length, it is BOTH
                let dev = point_line_deviation(&p2d0, line);
                if self.u_closed && line.direction().x().abs() < PCONFUSION {
                    if (dev.abs() - self.u_period).abs() < 0.1 * self.u_period {
                        code = IOR_BOTH;
                        if dev > 0.0 {
                            code |= IOR_POS;
                        }
                    } else if code == IOR_BOTH {
                        code = IOR_UNDEF;
                    }
                }
                if self.v_closed && line.direction().y().abs() < PCONFUSION {
                    if (dev.abs() - self.v_period).abs() < 0.1 * self.v_period {
                        code = IOR_BOTH;
                        if dev > 0.0 {
                            code |= IOR_POS;
                        }
                    } else if code == IOR_BOTH {
                        code = IOR_UNDEF;
                    }
                }
            }
            return code;
        }
        if i != 0 {
            code = IOR_UNDEF; // tangency
        } else if code == IOR_BOTH {
            // parity error in intersector
            code = IOR_LEFT;
            self.status |= SHAPEEXTEND_FAIL2;
        }
        code
    }
}
