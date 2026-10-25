//! `TopOpeBRepTool_CORRISO` — repair the 2-D representation of a shape built on
//! a periodic surface so that it becomes UV-connected.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_CORRISO.cxx`
//! (OCCT 8.0.0, 1,399 lines) and `TopOpeBRepTool_CORRISO.hxx:34-129`.
//!
//! `Fref` is a face on an `x`-periodic surface (`x` is `u` or `v`). `S` is built
//! on `Fref`'s geometry and should be UV-closed. An edge of `S` whose pcurve
//! falls outside the `Fref` window is "faulty"; the tool translates that pcurve
//! by one period and writes the corrected curves back (`GetnewS`).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `TopOpeBRepTool_CORRISO(Fref)` | `84-93` | [`CorrIso::new`] |
//! | `Refclosed` | `112-126` | [`CorrIso::ref_closed`] |
//! | `Init` | `128-193` | [`CorrIso::init`] |
//! | `UVClosed` | `269-281` | [`CorrIso::uv_closed`] |
//! | `Tol` | `286-290` | [`CorrIso::tol`] |
//! | `FUN_getx` | `291-300` | `fun_getx` |
//! | `PurgeFyClosingE` | `304-712` | [`CorrIso::purge_fy_closing_e`] |
//! | `FUN_tool_recadre` | `714-749` | `fun_tool_recadre` |
//! | `EdgeOUTofBoundsUV` | `753-855` | [`CorrIso::edge_out_of_bounds_uv`] |
//! | `EdgesOUTofBoundsUV` | `859-885` | [`CorrIso::edges_out_of_bounds_uv`] |
//! | `EdgeWithFaultyUV(E, Ivfaulty)` | `889-1069` | [`CorrIso::edge_with_faulty_uv`] |
//! | `EdgesWithFaultyUV` | `1073-1138` | [`CorrIso::edges_with_faulty_uv`] |
//! | `EdgeWithFaultyUV(list, n, fyE, I)` | `1142-1157` | [`CorrIso::edge_with_faulty_uv_first`] |
//! | `TrslUV` | `1161-1215` | [`CorrIso::trsl_uv`] |
//! | `GetnewS` | `1236-1287` | [`CorrIso::get_new_s`] |
//! | `AddNewConnexity` | `1293-1340` | [`CorrIso::add_new_connexity`] |
//! | `RemoveOldConnexity` | `1344-1399` | [`CorrIso::remove_old_connexity`] |
//!
//! Notes kept from the OCCT source rather than "fixed":
//! - `EdgeOUTofBoundsUV` answers `SPLITEDGE` (`0`) for an edge that is not
//!   bound to a 2-D representation, because `cxx:762-765` returns `false` from
//!   a function whose result is an `int`.
//! - `GetnewS` only writes the two-pcurve seam form for a `FORWARD` edge
//!   (`cxx:1278-1281`); a `REVERSED` occurrence of the same edge writes
//!   nothing, because the map it reads is orientation-insensitive.
//! - the `NYIxpu271198` arm of `EdgeOUTofBoundsUV` (`cxx:845-851`) answers
//!   `SPLITEDGE` without computing the split parameter; there is no
//!   intersection of the pcurve with the `x`-iso seam yet.
//!
//! UNPORTED:
//! - the `OCCT_DEBUG` blocks (`FUN_RaiseError`, `FUN_Raise`, `FUN_tool_trace`,
//!   `STATIC_PURGE_mapv` / `STATIC_PURGE_mapeds`); they only print.
//! - `FUN_tool_correctdgE` (`cxx:1219-1232`) is commented out in OCCT.
//! - the `#ifdef OCCT_DEBUG` copy around `UVISO` at `cxx:1252-1257` is a
//!   debug-only duplicate of the same call.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpPnt2d, GpTrsf2d, GpVec2d};
use occt_geom::approx_same_parameter::{u_resolution, v_resolution};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::trimmed::Geom2dTrimmedCurve;

use crate::abs::{Orientation, ShapeType};
use crate::bnd_lib_add2d;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{
    explorer_shapes, is_same_occurrence, shape_map_key, ShapeMapKey,
};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;
use crate::topopebrep::tool::c2d::{
    fc2d_add_new_curve_on_surface, fc2d_editable_curve_on_surface, fc2d_has_old_curve_on_surface,
    C2dF,
};
use crate::topopebrep::tool::closed::{par_e, vertex, vertices};
use crate::topopebrep::tool::iso::{uv_iso_c2df, uv_iso_curve, uvf};

/// `#define SPLITEDGE (0)` (`cxx:714`).
pub const SPLITEDGE: i32 = 0;
/// `#define INCREASE (1)` (`cxx:715`).
pub const INCREASE: i32 = 1;
/// `#define DECREASE (-1)` (`cxx:716`).
pub const DECREASE: i32 = -1;
/// The `recadre = 10; // INIT` sentinel (`cxx:721`, `:758`).
pub const INIT: i32 = 10;

/// `FUN_getx(E, c2df, uiso, par)` (`cxx:291-300`): the `x` coordinate the
/// closure test reads — `Y` for a `u`-iso line, `X` otherwise.
fn fun_getx(c2df: &C2dF, uiso: bool, par: f64) -> Option<f64> {
    let uv = uvf(par, c2df)?;
    Some(if uiso { uv.y() } else { uv.x() })
}

/// `FUN_tool_recadre(minx, maxx, xfirst, xlast, tolx, maxsup)` (`cxx:718-749`).
///
/// Answers `(recadre, maxsup)`; `maxsup` is the out-parameter OCCT fills even
/// when `recadre` is not `SPLITEDGE` (`cxx:725`).
fn fun_tool_recadre(
    minx: f64,
    maxx: f64,
    xfirst: f64,
    xlast: f64,
    tolx: f64,
) -> (i32, bool) {
    let maxinf = maxx < xfirst + tolx;
    let mininf = minx < xfirst - tolx;
    let maxsup = maxx > xlast + tolx;
    let minsup = minx > xlast - tolx;
    let maxok = xfirst - tolx < maxx && maxx < xlast + tolx;
    let minok = xfirst - tolx < minx && minx < xlast + tolx;

    let recadre = if maxinf {
        INCREASE
    } else if minsup {
        DECREASE
    } else if mininf && maxok {
        SPLITEDGE
    } else if minok && maxsup {
        SPLITEDGE
    } else {
        INIT
    };
    (recadre, maxsup)
}

/// `TopOpeBRepTool_TOOL::TrslUV(t2d, C2DF)` (`TopOpeBRepTool_TOOL.cxx:1312-1317`):
/// `PC->Translate(t2d)` followed by `C2DF.SetPC(PC, f, l, tol)`.
///
/// The port answers a translated copy instead of mutating in place; the caller
/// writes it back through [`CorrIso::set_uv_rep`], which is the same state.
fn trsl_uv_c2df(t2d: &GpVec2d, c2df: &mut C2dF) {
    let (pc, first, last, tol) = c2df.pc();
    let Some(pc) = pc else {
        return;
    };
    let mut tr = GpTrsf2d::identity();
    tr.set_translation_vec(t2d);
    let moved: Arc<dyn Curve2d> = Arc::from(pc.transformed(&tr));
    c2df.set_pc(moved, first, last, tol);
}

/// `TopOpeBRepTool_CORRISO` (`TopOpeBRepTool_CORRISO.hxx:53-125`).
///
/// `myGAS` is the adaptor of the **face's own** surface
/// (`cxx:91-92`), not of its basis surface — [`CorrIso::tol`] must answer the
/// resolution of the geometry `S` is built on.
pub struct CorrIso {
    /// `myFref`
    f_ref: TopoShape,
    /// `myGAS`
    gas: SurfaceAdaptor,
    /// `myUclosed`
    u_closed: bool,
    /// `myVclosed`
    v_closed: bool,
    /// `myUper`
    u_per: f64,
    /// `myVper`
    v_per: f64,
    /// `myS`
    s: Option<TopoShape>,
    /// `myEds`
    eds: Vec<TopoShape>,
    /// `myERep2d` — a `DataMap` over `TopTools_ShapeMapHasher`, i.e. keyed on
    /// `(TShape, Location)` with the orientation ignored.
    e_rep_2d: HashMap<ShapeMapKey, C2dF>,
    /// `myVEds`
    v_eds: HashMap<ShapeMapKey, Vec<TopoShape>>,
}

impl CorrIso {
    /// `TopOpeBRepTool_CORRISO(const TopoDS_Face& Fref)` (`cxx:86-93`).
    ///
    /// OCCT's `TopOpeBRepTool_CORRISO()` (`cxx:82`) is `= default`, which
    /// leaves the four closure members indeterminate; the port has no such
    /// constructor, every instance starts from a reference face.
    pub fn new(f_ref: &TopoShape) -> Self {
        let face = Face(f_ref.clone());
        let (closed, u_closed, u_per, v_closed, v_per) =
            crate::topopebrep::tool::geometry::fun_tool_closed_s(&face);
        let _ = closed;
        Self {
            f_ref: f_ref.clone(),
            gas: SurfaceAdaptor::of_face(&face),
            u_closed,
            v_closed,
            u_per,
            v_per,
            s: None,
            eds: Vec::new(),
            e_rep_2d: HashMap::new(),
            v_eds: HashMap::new(),
        }
    }

    /// `Fref()` (`cxx:96-99`).
    pub fn f_ref(&self) -> &TopoShape {
        &self.f_ref
    }

    /// `GASref()` (`cxx:102-106`).
    pub fn gas_ref(&self) -> &SurfaceAdaptor {
        &self.gas
    }

    /// `Refclosed(int x, double& xperiod)` (`cxx:112-126`): `Some(period)` for
    /// a closed direction, `None` for `false` (including an `x` other than 1/2).
    pub fn ref_closed(&self, x: i32) -> Option<f64> {
        match x {
            1 => {
                if self.u_closed {
                    Some(self.u_per)
                } else {
                    None
                }
            }
            2 => {
                if self.v_closed {
                    Some(self.v_per)
                } else {
                    None
                }
            }
            _ => None,
        }
    }

    /// `Init(S)` (`cxx:128-193`): collect the edges of `S`, their editable
    /// pcurve on `Fref`, and the vertex-to-edges connexity.
    ///
    /// `FC2D_EditableCurveOnSurface` is called with its default
    /// `trim3d = false` (`TopOpeBRepTool_2d.hxx:77-82`).
    pub fn init(&mut self, s: &TopoShape) -> bool {
        self.e_rep_2d.clear();
        self.eds.clear();
        self.v_eds.clear();

        if s.is_null() {
            return false;
        }
        self.s = Some(s.clone());

        for e in explorer_shapes(s, ShapeType::Edge) {
            self.eds.push(e.clone());

            let has_old = fc2d_has_old_curve_on_surface(&e, &self.f_ref).is_some();
            let Some(c2df) = fc2d_editable_curve_on_surface(&e, &self.f_ref, false) else {
                return false;
            };
            if !has_old {
                if let Some(pc) = c2df.curve.clone() {
                    fc2d_add_new_curve_on_surface(
                        &pc, &e, &self.f_ref, c2df.first, c2df.last, c2df.tol,
                    );
                }
            }
            let key = shape_map_key(&e);
            self.e_rep_2d.insert(key, c2df);

            for v in explorer_shapes(&e, ShapeType::Vertex) {
                self.v_eds.entry(shape_map_key(&v)).or_default().push(e.clone());
            }
        }
        true
    }

    /// `S()` (`cxx:198-201`).
    pub fn s(&self) -> Option<&TopoShape> {
        self.s.as_ref()
    }

    /// `Eds()` (`cxx:204-207`).
    pub fn eds(&self) -> &[TopoShape] {
        &self.eds
    }

    /// `UVClosed()` (`cxx:269-281`).
    pub fn uv_closed(&self) -> bool {
        let nfybounds = 3;
        let stopatfirst = true;
        let found = self.edges_with_faulty_uv(&self.eds, nfybounds, stopatfirst);
        found.is_empty()
    }

    /// `Tol(I, tol3d)` (`cxx:286-290`): `I == 1` is `UResolution`, everything
    /// else `VResolution`. A face without a surface has no adaptor geometry to
    /// resolve against.
    pub fn tol(&self, i: i32, tol3d: f64) -> f64 {
        let Some(s) = self.gas.surface() else {
            return 0.0;
        };
        if i == 1 {
            u_resolution(s.as_ref(), tol3d)
        } else {
            v_resolution(s.as_ref(), tol3d)
        }
    }

    /// `PurgeFyClosingE(ClEds, fyClEds)` (`cxx:304-712`): which of the closing
    /// edges of a wire have a faulty UV representation. `None` is OCCT's
    /// `false`; `Some(list)` is `true` with `fyClEds` filled.
    pub fn purge_fy_closing_e(&self, cl_eds: &[TopoShape]) -> Option<Vec<TopoShape>> {
        let mut fy_cl_eds: Vec<TopoShape> = Vec::new();

        // `cxx:311-316`. `tttol` (1.e-8) is the running tolerance of the third
        // block below; `tttuvF` is the reference face's parametric tolerance.
        let mut tttol: f64 = 1.0e-8;
        let tttol_s = BRepTool::face_tolerance(&Face(self.f_ref.clone()));
        let tolu = self.tol(1, tttol_s);
        let tolv = self.tol(2, tttol_s);
        let tttuv_f = tolu.max(tolv);

        let mapcl: Vec<ShapeMapKey> = cl_eds.iter().map(shape_map_key).collect();
        let contains = |e: &TopoShape| mapcl.contains(&shape_map_key(e));

        // `cxx:328-332`
        let fyceds = self.edges_with_faulty_uv(cl_eds, 3, false);
        if fyceds.is_empty() {
            return None;
        }

        if fyceds.len() == 1 {
            // `cxx:335-397` (`ivf == 3 : cto016G*`)
            let (c_e, ivf) = (fyceds[0].0.clone(), fyceds[0].1);
            let fyeds = self.edges_with_faulty_uv(&self.eds, 3, false);
            let nfy = fyeds.len();

            let o_oc_e = c_e.orientation().complement();
            let mut ison_c_e = contains(&c_e.clone().oriented(o_oc_e));
            if ison_c_e {
                // `cxx:350-352`: `vcE` (`Vertices(cE, vcE)`) and `vOcE`
                // (`cxx:368-369`) are filled and never read.
                let oc_e = c_e.orientation();
                let tttol_c_e = BRepTool::edge_tolerance(&Edge(c_e.clone()));
                let tttuv_c_e = self.tol(1, tttol_c_e).max(self.tol(2, tttol_c_e));
                let c_e_2d = self.uv_rep(&c_e)?;

                let oo_oc_e = oc_e.complement();
                let oc_e_shape = c_e.clone().oriented(oo_oc_e);
                let tttol_oc_e = BRepTool::edge_tolerance(&Edge(oc_e_shape.clone()));
                let tttuv_oc_e = self.tol(1, tttol_oc_e).max(self.tol(2, tttol_oc_e));
                let oc_e_2d = self.uv_rep(&oc_e_shape)?;

                let parvce1 = par_e(1, &c_e)?;
                let uvvce1 = uvf(parvce1, &c_e_2d)?;

                let parv_oc_e2 = par_e(2, &oc_e_shape)?;
                let uvv_oc_e2 = uvf(parv_oc_e2, &oc_e_2d)?;

                let tol = tttuv_c_e.max(tttuv_oc_e);
                ison_c_e = uvvce1.distance(&uvv_oc_e2) < tol;
                if ison_c_e && (nfy != 1) {
                    // `cxx:387-390` (`cto009L2`)
                    return None;
                }
            }

            if ivf == 3 {
                fy_cl_eds.push(c_e);
                return Some(fy_cl_eds);
            }
        } else if fyceds.len() > 1 {
            // `cxx:399-491` (`ivf == 1,2 : cto016E*`): the faulty edges
            // describe a closing edge whose first and last UV bounds are not
            // connexed; such edges are removed.
            let mut hasinit = false;
            let mut isou = false;
            let mut isov = false;
            let mut o2d = GpPnt2d::new(0.0, 0.0);
            let mut xinf = 1.0e7;
            let mut xsup = -1.0e7;
            let mut infdef = false;
            let mut supdef = false;

            for (c_e, ivf) in fyceds.iter() {
                let c2df = self.uv_rep(c_e)?;
                // `cxx:424`: the `d2d` out-parameter is read nowhere below.
                let iso = uv_iso_c2df(&c2df)?;

                if hasinit {
                    let onsamline = (isou && iso.iso_u) || (isov && iso.iso_v);
                    if !onsamline {
                        return None;
                    }
                }
                if !hasinit {
                    isou = iso.iso_u;
                    isov = iso.iso_v;
                    o2d = iso.origin;
                    hasinit = true;
                } else {
                    let mut onsamline = false;
                    if isou && iso.iso_u {
                        let du = o2d.x() - iso.origin.x();
                        onsamline = du.abs() < tolu;
                    }
                    if isov && iso.iso_v {
                        let dv = o2d.y() - iso.origin.y();
                        onsamline = dv.abs() < tolv;
                    }
                    if !onsamline {
                        return None;
                    }
                }

                for i in 1..=2 {
                    let pari = par_e(i, c_e)?;
                    let xi = fun_getx(&c2df, isou, pari)?;
                    let vifaulty = *ivf == i || *ivf == 3;
                    if xi < xinf {
                        xinf = xi;
                        infdef = vifaulty;
                    }
                    if xi > xsup {
                        xsup = xi;
                        supdef = vifaulty;
                    }
                }
                fy_cl_eds.push(c_e.clone());
            }
            let toremove = infdef && supdef;
            if !toremove {
                fy_cl_eds.clear();
            }
        }

        if !fy_cl_eds.is_empty() {
            return Some(fy_cl_eds);
        }

        // `cxx:495-711`: do the two closing edges have their 2-D
        // representations "confounded"?
        for c_e in cl_eds.iter() {
            let oc_e = c_e.orientation();
            let tttol_c_e = BRepTool::edge_tolerance(&Edge(c_e.clone()));
            let tttuv_c_e = self.tol(1, tttol_c_e).max(self.tol(2, tttol_c_e));
            let c_e_2d = self.uv_rep(c_e)?;

            let ison_oc_e_2d = {
                let oo_oc_e = oc_e.complement();
                let oc_e_shape = c_e.clone().oriented(oo_oc_e);
                if !contains(&oc_e_shape) {
                    // `cxx:530-533`: the closing edge appears twice.
                    continue;
                }
                let tttol_oc_e = BRepTool::edge_tolerance(&Edge(oc_e_shape.clone()));
                let tttuv_oc_e = self.tol(1, tttol_oc_e).max(self.tol(2, tttol_oc_e));
                let oc_e_2d = self.uv_rep(&oc_e_shape)?;

                let parvce1 = par_e(1, c_e)?;
                let uvvce1 = uvf(parvce1, &c_e_2d)?;

                let parv_oc_e2 = par_e(2, &oc_e_shape)?;
                let uvv_oc_e2 = uvf(parv_oc_e2, &oc_e_2d)?;
                let tol = tttuv_c_e.max(tttuv_oc_e);
                uvvce1.distance(&uvv_oc_e2) < tol
            };
            if !ison_oc_e_2d {
                continue;
            }

            let mut nvc_e_ok = 0;
            for ivce in 1..=2 {
                let vce = vertex(ivce, c_e)?;
                let loe = self.connexity(&vce)?;

                let parvce = par_e(ivce, c_e)?;
                let uvvce = uvf(parvce, &c_e_2d)?;
                let tttol_vce = BRepTool::vertex_tolerance(&crate::shape::Vertex(vce.clone()));
                let tttuv_vce = self.tol(1, tttol_vce).max(self.tol(2, tttol_vce));

                let mut vceok = false;
                for e in loe.iter() {
                    if contains(e) {
                        // `cxx:601-605`: no connexity check on closing edges
                        // (`xpu090399 cto016E1`).
                        continue;
                    }
                    let e_2d = self.uv_rep(e)?;
                    let tttol_e = BRepTool::edge_tolerance(&Edge(e.clone()));
                    let tttuv_e = self.tol(1, tttol_e).max(self.tol(2, tttol_e));

                    let (ve1, ve2) = vertices(e);
                    for ive in 1..=2 {
                        let ve = if ive == 1 { ve1.clone() } else { ve2.clone() };
                        let Some(ve) = ve else {
                            continue;
                        };
                        if !is_same_occurrence(&ve, &vce) {
                            continue;
                        }
                        let parve = par_e(ive, e)?;
                        let uvve = uvf(parve, &e_2d)?;
                        if ive == ivce {
                            // `cxx:637-640`: a vertex FORWARD-connexed to a
                            // REVERSED one.
                            continue;
                        }
                        let tttol_ve =
                            BRepTool::vertex_tolerance(&crate::shape::Vertex(ve.clone()));
                        let tttuv_ve = self.tol(1, tttol_ve).max(self.tol(2, tttol_ve));

                        tttol = tttol.max(tttuv_f.max(
                            tttuv_e.max(tttuv_c_e.max(tttuv_ve.max(tttuv_vce))),
                        ));

                        let mut dd = uvve.distance(&uvvce);
                        let mut sameuv = dd < tttol;
                        if self.u_closed {
                            let xperiod = self.u_per;
                            dd = uvve.x() - uvvce.x();
                            sameuv = sameuv || ((dd.abs() - xperiod).abs() < tttol);
                        }
                        if self.v_closed {
                            let xperiod = self.v_per;
                            dd = uvve.y() - uvvce.y();
                            sameuv = sameuv || ((dd.abs() - xperiod).abs() < tttol);
                        }
                        if sameuv {
                            vceok = true;
                        }
                        break;
                    }
                    if vceok {
                        break;
                    }
                }
                if vceok {
                    nvc_e_ok += 1;
                }
            }

            if nvc_e_ok == 0 {
                // `cxx:694-694`: neither bound is connexed to a non-closed edge.
                fy_cl_eds.push(c_e.clone());
            }
        }

        if fy_cl_eds.is_empty() {
            None
        } else {
            Some(fy_cl_eds)
        }
    }

    /// `UVRep(E, C2DF)` (`cxx:212-222`).
    pub fn uv_rep(&self, e: &TopoShape) -> Option<C2dF> {
        self.e_rep_2d.get(&shape_map_key(e)).cloned()
    }

    /// `SetUVRep(E, C2DF)` (`cxx:226-235`): only an already-bound edge is
    /// replaced.
    pub fn set_uv_rep(&mut self, e: &TopoShape, c2df: &C2dF) -> bool {
        match self.e_rep_2d.get_mut(&shape_map_key(e)) {
            Some(slot) => {
                *slot = c2df.clone();
                true
            }
            None => false,
        }
    }

    /// `Connexity(V, Eds)` (`cxx:239-250`).
    pub fn connexity(&self, v: &TopoShape) -> Option<Vec<TopoShape>> {
        self.v_eds.get(&shape_map_key(v)).cloned()
    }


    /// `SetConnexity(V, Eds)` (`cxx:254-265`): only an already-bound vertex is
    /// replaced.
    pub fn set_connexity(&mut self, v: &TopoShape, eds: &[TopoShape]) -> bool {
        match self.v_eds.get_mut(&shape_map_key(v)) {
            Some(slot) => {
                *slot = eds.to_vec();
                true
            }
            None => false,
        }
    }
}

/// `FyEds.Bind(E, ivf)` (`cxx:875`, `:1126`): a `DataMap` bind, so an edge
/// already in the list keeps its position and takes the new value.
pub(crate) fn bind_fy(fy_eds: &mut Vec<(TopoShape, i32)>, e: &TopoShape, v: i32) {
    let key = shape_map_key(e);
    if let Some(slot) = fy_eds.iter_mut().find(|(s, _)| shape_map_key(s) == key) {
        slot.1 = v;
        return;
    }
    fy_eds.push((e.clone(), v));
}

mod faulty_uv;
mod translate;

