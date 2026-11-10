//! `TopOpeBRepTool_mkTondgE` — the transition of a degenerated edge on a face.
//!
//! Source (OCCT 8.0.0, `ModelingAlgorithms/TKBool/TopOpeBRepTool`):
//! - `Initialize(dgE, F, uvi, Fi)` (`TopOpeBRepTool_mkTondgE.cxx:58-101`)
//! - `SetclE(clE)` (`cxx:105-108`), `IsT2d()` (`cxx:113-116`)
//! - `SetRest(pari, Ei)` (`cxx:121-132`)
//! - `MkTonE(mkT, par1, par2)` (`cxx:272-329`)
//! - `MkTonE(ei, mkT, par1, par2)` (`cxx:334-530`)
//! - `FUN_getEc` (`cxx:192-215`) and `FUN_MkTonE` (`cxx:217-267`), the two
//!   file-local statics both `MkTonE` arms walk through
//!
//! UNPORTED: `GetAllRest(lEi)` (`cxx:135-191`). Its only OCCT caller is the
//! commented-out block at `TopOpeBRep_vprdeg.cxx:615-621`, so nothing ported
//! can reach it.

use occt_core::gp::{GpDir, GpPnt2d, GpVec};
use occt_core::precision::ANGULAR;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::shape::{Face, TopoShape};
use crate::topexp::Explorer;
use crate::topo_tools::is_same;
use crate::topopebrep::tool::closed::{closed_e, is_closing_e, FORWARD};
use crate::topopebrep::tool::curvature::xx;
use crate::topopebrep::tool::iso::orientation_is_forward;
use crate::topopebrep::tool::matter::matter;
use crate::topopebrep::tool::project_curve::fun_tool_bounds;
use crate::topopebrep::tool::tangent::{tg_geom_e, tg_inside};
use crate::topopebrep::tool::topology::{fun_tool_ng_geom_f, fun_tool_par_on_e_f, tol_p};

/// `NOI`: no interference to build.
pub const NOI: i32 = 0;
/// `MKI1`: one interference, at `par1`.
pub const MKI1: i32 = 1;
/// `MKI2`: one interference, at `par2`.
pub const MKI2: i32 = 2;
/// `MKI12`: both interferences.
pub const MKI12: i32 = 3;

/// `TopOpeBRepTool_mkTondgE`.
#[derive(Clone)]
pub struct MktondgE {
    is_t2d: bool,
    has_rest: bool,
    cl_e: Option<TopoShape>,
    /// `myEpari`: `NCollection_DataMap<TopoDS_Shape, double>`. Kept as an
    /// association list; it holds only the face's non-closing restrictions
    /// bound through [`MktondgE::set_rest`].
    epari: Vec<(TopoShape, f64)>,
    dg_e: TopoShape,
    f: TopoShape,
    ngf: GpDir,
    uvi: GpPnt2d,
    fi: TopoShape,
    ngfi: GpDir,
}

impl Default for MktondgE {
    fn default() -> Self {
        Self::new()
    }
}

impl MktondgE {
    /// `TopOpeBRepTool_mkTondgE()` (`cxx:56`, `= default`). The C++ members are
    /// uninitialised; the port starts from neutral values.
    pub fn new() -> Self {
        Self {
            is_t2d: false,
            has_rest: false,
            cl_e: None,
            epari: Vec::new(),
            dg_e: TopoShape::new(ShapeType::Edge),
            f: TopoShape::new(ShapeType::Face),
            ngf: GpDir::default_dir(),
            uvi: GpPnt2d::new(0.0, 0.0),
            fi: TopoShape::new(ShapeType::Face),
            ngfi: GpDir::default_dir(),
        }
    }

    /// `Initialize(dgE, F, uvi, Fi)` (`cxx:58-101`): the normal of `F` at the
    /// first vertex of `dgE`, the normal of `Fi` at `uvi`, and `isT2d` from
    /// their parallelism.
    pub fn initialize(
        &mut self,
        dg_e: &TopoShape,
        f: &TopoShape,
        uvi: &GpPnt2d,
        fi: &TopoShape,
    ) -> bool {
        self.is_t2d = false;
        self.has_rest = false;
        self.cl_e = None;
        self.epari.clear();
        self.dg_e = dg_e.clone();
        self.f = f.clone();

        let Some(v) = first_vertex(dg_e) else {
            return false;
        };
        let par = BRepTool::parameter_of_vertex(dg_e, &v);
        let Some(uv) = fun_tool_par_on_e_f(dg_e, par, f) else {
            return false;
        };
        let Some(ngf) = fun_tool_ng_geom_f(&uv, &Face(f.clone())) else {
            return false;
        };
        self.ngf = ngf;

        self.uvi = *uvi;
        self.fi = fi.clone();
        let Some(ngfi) = fun_tool_ng_geom_f(uvi, &Face(fi.clone())) else {
            return false;
        };
        self.ngfi = ngfi;

        let dot = self.ngf.dot(&self.ngfi);
        self.is_t2d = (1.0 - dot.abs()).abs() < ANGULAR;
        true
    }

    /// `SetclE(clE)` (`cxx:105-108`): always `true` in OCCT.
    pub fn set_cl_e(&mut self, cl_e: &TopoShape) -> bool {
        self.cl_e = Some(cl_e.clone());
        true
    }

    /// `IsT2d()` (`cxx:113-116`).
    pub fn is_t2d(&self) -> bool {
        self.is_t2d
    }

    /// `SetRest(pari, Ei)` (`cxx:121-132`): `false` for a closing restriction,
    /// which also clears `hasRest`.
    pub fn set_rest(&mut self, pari: f64, ei: &TopoShape) -> bool {
        self.has_rest = true;
        if is_closing_e(ei, &self.fi) {
            self.has_rest = false;
            return false;
        }
        self.epari.push((ei.clone(), pari));
        true
    }

    /// `myclE`, resolved on first use by `FUN_getEc` (`cxx:285-291`).
    fn cl_e_or(&mut self) -> Option<TopoShape> {
        if let Some(c) = &self.cl_e {
            return Some(c.clone());
        }
        let v = first_vertex(&self.dg_e)?;
        let cle = fun_get_ec(&self.f, &v)?;
        self.cl_e = Some(cle.clone());
        Some(cle)
    }

    /// `MkTonE(mkT, par1, par2)` (`cxx:272-329`): the `fi`-level arm, only
    /// meaningful when `isT2d` is `false`.
    pub fn mk_ton_e(&mut self, par1: &mut f64, par2: &mut f64) -> i32 {
        if self.is_t2d {
            return NOI;
        }
        *par1 = 1.0e7;
        *par2 = 1.0e7;

        let Some(cl_e) = self.cl_e_or() else {
            return NOI;
        };
        let Some(v) = first_vertex(&self.dg_e) else {
            return NOI;
        };
        let Some((dir_in_cle, ovcle)) = tg_inside(&v, &cl_e) else {
            return NOI;
        };
        let mut faxis = as_vec(&self.ngf);
        if ovcle == FORWARD {
            faxis.reverse();
        }
        let Some(xx) = fun_tool_ng_geom_f(&self.uvi, &Face(self.fi.clone())) else {
            return NOI;
        };
        let mut xxi = as_vec(&xx);
        if orientation_is_forward(&self.fi) {
            xxi.reverse();
        }
        match fun_mk_ton_e(&faxis, &dir_in_cle, &xxi, par1, par2) {
            Some(_) => MKI12,
            None => NOI,
        }
    }

    /// `MkTonE(ei, mkT, par1, par2)` (`cxx:334-530`): the restriction-level arm.
    ///
    /// The returned `i32` is the OCCT `mkT` out-parameter; `false` mirrors the
    /// OCCT `return false` arms (`hasRest` unbound, no tangent, no normal).
    pub fn mk_ton_e_on_rest(
        &mut self,
        ei: &TopoShape,
        par1: &mut f64,
        par2: &mut f64,
    ) -> (i32, bool) {
        *par1 = 1.0e7;
        *par2 = 1.0e7;
        let bound = self.epari.iter().find(|(e, _)| is_same(e, ei));
        self.has_rest = bound.is_some();
        let Some((_, pari)) = bound else {
            return (NOI, false);
        };
        let pari = *pari;

        let (pfi, pli) = fun_tool_bounds(ei);
        let Some(tolpi) = tol_p(ei, &self.fi) else {
            return (NOI, false);
        };
        let onfi = (pari - pfi).abs() < tolpi;
        let onli = (pari - pli).abs() < tolpi;
        let Some(mut tgin1di) = tg_geom_e(pari, ei) else {
            return (NOI, false);
        };
        if onli {
            tgin1di.reverse();
        }

        let Some(cl_e) = self.cl_e_or() else {
            return (NOI, false);
        };
        let Some(v) = first_vertex(&self.dg_e) else {
            return (NOI, false);
        };
        let Some((dir_in_cle, ovcle)) = tg_inside(&v, &cl_e) else {
            return (NOI, false);
        };

        let mut faxis = as_vec(&self.ngf);
        if ovcle == FORWARD {
            faxis.reverse();
        }
        let closedi = closed_e(ei).is_some();

        if self.is_t2d {
            // `cxx:417-442`
            let Some(x) = xx(&self.uvi, &self.fi, pari, ei) else {
                return (NOI, false);
            };
            let xxi = as_vec(&x);
            let Some(outin) = fun_mk_ton_e(&faxis, &dir_in_cle, &xxi, par1, par2) else {
                return (NOI, false);
            };
            if !onfi && !onli {
                return (MKI12, true);
            }
            if closedi {
                return (MKI12, true);
            }
            let xxri = GpDir::from_vec(&tgin1di).unwrap_or_else(|_| GpDir::default_dir());
            return (mk_from_keep_ang(&xxi, &faxis, &xxri, outin), true);
        }

        // `cxx:444-497`
        let Some(x) = fun_tool_ng_geom_f(&self.uvi, &Face(self.fi.clone())) else {
            return (NOI, false);
        };
        let mut xxi = as_vec(&x);
        if orientation_is_forward(&self.fi) {
            xxi.reverse();
        }
        let Some(outin) = fun_mk_ton_e(&faxis, &dir_in_cle, &xxi, par1, par2) else {
            return (NOI, false);
        };
        let Some(xxri) = xx(&self.uvi, &self.fi, pari, ei) else {
            return (NOI, false);
        };

        let tgi = xxi.crossed(&faxis);
        let mut dot = tgi.dot(&as_vec(&xxri));
        if dot.abs() < ANGULAR {
            if (!onfi && !onli) || closedi {
                return (MKI12, true);
            }
            dot = tgi.dot(&tgin1di);
        }
        let keep_ang = dot > 0.0;
        let mkt = if keep_ang == outin { MKI1 } else { MKI2 };
        (mkt, true)
    }
}

/// The `mkT` tail shared by both `MkTonE(ei, ...)` arms (`cxx:517-528`):
/// `tgi = xxi ^ faxis`, `keepang = tgi . xxri`, then `mkT` from `keepang` /
/// `outin`.
fn mk_from_keep_ang(xxi: &GpVec, faxis: &GpVec, xxri: &GpDir, outin: bool) -> i32 {
    let tgi = xxi.crossed(faxis);
    let dot = tgi.dot(&as_vec(xxri));
    let keep_ang = dot > 0.0;
    if keep_ang == outin {
        MKI1
    } else {
        MKI2
    }
}

/// `FUN_MkTonE(faxis, dirINcle, xxi, ngf, par1, par2, outin)`
/// (`cxx:217-267`): the `ang` / `Cang` pair. `Some(outin)` mirrors OCCT's
/// always-`true` return.
fn fun_mk_ton_e(
    faxis: &GpVec,
    dir_in_cle: &GpVec,
    xxi: &GpVec,
    par1: &mut f64,
    par2: &mut f64,
) -> Option<bool> {
    let tgi = xxi.crossed(faxis);
    let tola = ANGULAR;
    let dot1 = dir_in_cle.dot(xxi);
    let is_oni = dot1.abs() < tola;

    let ang;
    let outin;
    if is_oni {
        let dot = dir_in_cle.dot(&tgi);
        ang = if dot > 0.0 { 0.0 } else { std::f64::consts::PI };
        outin = true;
    } else {
        ang = matter(dir_in_cle, &tgi.reversed(), faxis);
        let dot = dir_in_cle.crossed(&tgi).dot(faxis);
        outin = if dot1 < 0.0 { dot > 0.0 } else { dot < 0.0 };
    }

    let cang = if ang > std::f64::consts::PI {
        ang - std::f64::consts::PI
    } else {
        ang + std::f64::consts::PI
    };
    *par1 = if outin { ang } else { cang };
    *par2 = if outin { cang } else { ang };
    Some(outin)
}

/// `FUN_getEc(f, v, cle)` (`cxx:192-215`): the closing edge of `f` that carries
/// `v`.
fn fun_get_ec(f: &TopoShape, v: &TopoShape) -> Option<TopoShape> {
    let mut ex = Explorer::new(f, ShapeType::Edge);
    while ex.more() {
        let e = ex.current();
        if is_closing_e(&e, f) {
            let mut exv = Explorer::new(&e, ShapeType::Vertex);
            while exv.more() {
                if is_same(&exv.current(), v) {
                    return Some(e.clone());
                }
                exv.next();
            }
        }
        ex.next();
    }
    None
}

/// `TopExp_Explorer exv(mydgE, TopAbs_VERTEX); exv.Current()` (`cxx:72`,
/// `:282`, `:374`).
fn first_vertex(e: &TopoShape) -> Option<TopoShape> {
    let ex = Explorer::new(e, ShapeType::Vertex);
    if ex.more() {
        Some(ex.current().clone())
    } else {
        None
    }
}

/// `gp_Vec(dir)` from a `gp_Dir`.
fn as_vec(d: &GpDir) -> GpVec {
    GpVec::new(d.x(), d.y(), d.z())
}
