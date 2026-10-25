//! Edge/vertex orientation-in-parent and parameter-transfer helpers of
//! `TopOpeBRepTool_TOPOLOGY` that `TopOpeBRepDS_FaceInterferenceTool` and
//! `TopOpeBRepDS_FIR` read.
//!
//! Sources (OCCT 8.0.0, `TKBool/TopOpeBRepTool/TopOpeBRepTool_TOPOLOGY.cxx`):
//! - `FUN_tool_orientVinE` (`:378-397`)
//! - `FUN_tool_orientEinF` (`:400-418`)
//! - `FUN_tool_orientEinFFORWARD` (`:421-439`)
//! - `FUN_tool_EboundF` (`:441-457`)
//! - `FUN_tool_nggeomF` (`:460-466`, the `(p2d, F)` overload)
//! - `FUN_tool_parF` (`:969-1002`, both overloads)
//! - `FUN_tool_paronEF` (`:1005-1036`, both overloads)
//! - `FUN_tool_findparinBAC` (`:1121-1126`)
//! - `FUN_tool_findPinBAC` (`:1113-1118`)
//! - `FUN_tool_findparinE` (`:1131-1136`)
//! - `FUN_tool_findPinE` (`:1138-1143`)
//! - `FUN_tool_value(UV, F, P)` (`:328-333`)
//! - `FUN_tool_parVonE(V, E, par)` (`:887-929`)
//! - `FUN_tool_nggeomF(paronE, E, F, ...)` (`:468-537`, both overloads)
//! - `FUN_tool_Eshared(V, F1, F2, Eshared)` (`:843-883`)
//! - `FUN_tool_EtgF` (`:541-558`) and `FUN_tool_EtgOOE` (`:561-583`)

use occt_core::bnd::BndBox;
use occt_core::gp::{GpDir, GpPnt, GpPnt2d, GpVec};
use occt_core::precision::CONFUSION;

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topexp::Explorer;

use super::adaptor::CurveAdaptor;
use super::closed::{basis_surface, vertices};
use super::geometry::fun_tool_ng_s;
use super::project::fun_tool_proj_pon_f;
use super::project_curve::fun_tool_value;

/// `FIRST (1)` (`cxx:374`).
pub const FIRST: i32 = 1;
/// `LAST (2)` (`cxx:375`).
pub const LAST: i32 = 2;
/// `CLOSING (3)` (`cxx:376`).
pub const CLOSING: i32 = 3;

/// `FUN_tool_orientVinE(v, e)` (`cxx:378-397`): `FIRST` / `LAST` / `CLOSING`
/// (both vertices), `0` when `v` is not a vertex of `e`.
pub fn fun_tool_orient_vin_e(v: &TopoShape, e: &TopoShape) -> i32 {
    let (vf, vl) = vertices(e);
    let visf = vf.as_ref().is_some_and(|x| crate::topo_tools::is_same(x, v));
    let visl = vl.as_ref().is_some_and(|x| crate::topo_tools::is_same(x, v));
    let mut result = 0;
    if visf {
        result = FIRST;
    }
    if visl {
        result = LAST;
    }
    if visf && visl {
        result = CLOSING;
    }
    result
}

/// `FUN_tool_orientEinF(E, F, oriEinF)` (`cxx:400-418`): the orientation of
/// `E` in the face `F` as found by `TopExp_Explorer(F, TopAbs_EDGE)`.
///
/// The C++ returns the orientation of the *first* occurrence with
/// `EF.IsSame(E)` and reports `ok = e.More()` — i.e. whether any occurrence was
/// seen. The port returns the same pair; `None` is `ok == false`, and the
/// caller sees `oriEinF` defaulted to `FORWARD` (`cxx:406`).
pub fn fun_tool_orient_ein_f(e: &TopoShape, f: &TopoShape) -> Option<Orientation> {
    let mut ori = Orientation::Forward;
    let mut found = false;
    let mut ex = Explorer::new(f, ShapeType::Edge);
    while ex.more() {
        let ef = ex.current();
        if crate::topo_tools::is_same(ef, e) {
            ori = ef.orientation();
            found = true;
            break;
        }
        ex.next();
    }
    if found {
        Some(ori)
    } else {
        None
    }
}

/// `FUN_tool_findAncestor(lF, E, Fanc)` (`cxx:1548-1564`): the first face of
/// `lF` that has `E` as a bound.
pub fn fun_tool_find_ancestor(lf: &[TopoShape], e: &TopoShape) -> Option<TopoShape> {
    for f in lf {
        if fun_tool_orient_ein_f(e, f).is_some() {
            return Some(f.clone());
        }
    }
    None
}

/// `FUN_tool_orientEinFFORWARD(E, F, oriEinF)` (`cxx:421-439`): the same on/// `F.Oriented(TopAbs_FORWARD)`; `oriEinF` is set to `EXTERNAL` when the edge
/// is not a bound of `F` (`cxx:431`).
///
/// The port has no `TopoDS_Shape::Oriented`; the face orientation is instead
/// composed out of the occurrence orientation, which is what
/// `TopoDS_Face FF = TopoDS::Face(F.Oriented(TopAbs_FORWARD))` changes.
pub fn fun_tool_orient_ein_f_forward(e: &TopoShape, f: &TopoShape) -> Option<Orientation> {
    let ori = fun_tool_orient_ein_f(e, f)?;
    Some(Orientation::compose(f.orientation(), ori))
}

/// `FUN_tool_EboundF(E, F)` (`cxx:441-457`): `true` when `E` is a closing edge
/// of `F` or its orientation in `F` is `INTERNAL`/`EXTERNAL`.
pub fn fun_tool_ebound_f(e: &TopoShape, f: &TopoShape) -> bool {
    let Some(ori) = fun_tool_orient_ein_f_forward(e, f) else {
        return false;
    };
    if BRepTool::is_closed_edge_face(&Edge(e.clone()), &Face(f.clone())) {
        return true;
    }
    ori == Orientation::Internal || ori == Orientation::External
}

/// `FUN_tool_nggeomF(p2d, F)` (`cxx:460-466`): `gp_Vec(FUN_tool_ngS(p2d,
/// BRep_Tool::Surface(F)))`.
pub fn fun_tool_ng_geom_f(p2d: &GpPnt2d, f: &Face) -> Option<GpDir> {
    let s = BRepTool::face_surface_world(f)?;
    let s = basis_surface(&s);
    Some(fun_tool_ng_s(p2d, s.as_ref()))
}

/// `FUN_tool_parF(E, par, F, UV, tol3d)` (`cxx:969-991`): the point at `par` on
/// `E`, projected on `F`, accepted when the projection distance is below
/// `tol3d`.
pub fn fun_tool_par_f_tol(
    e: &TopoShape,
    par: f64,
    f: &TopoShape,
    tol3d: f64,
) -> Option<GpPnt2d> {
    let p = fun_tool_value(par, e)?;
    let (u, v, dist) = fun_tool_proj_pon_f(&p, &Face(f.clone()))?;
    if dist < tol3d {
        Some(GpPnt2d::new(u, v))
    } else {
        None
    }
}

/// `FUN_tool_parF(E, par, F, UV)` (`cxx:994-1002`): `tol3d =
/// BRep_Tool::Tolerance(F) * 1.e2` (the `// KK xpu` provisional factor of
/// `cxx:1000`, kept verbatim).
pub fn fun_tool_par_f(e: &TopoShape, par: f64, f: &TopoShape) -> Option<GpPnt2d> {
    let tol3d = BRepTool::face_tolerance(&Face(f.clone())) * 1.0e2;
    fun_tool_par_f_tol(e, par, f, tol3d)
}

/// `FUN_tool_paronEF(E, par, F, UV, tol3d)` (`cxx:1005-1025`).
///
/// `FC2D_CurveOnSurface` is the port's stored-then-built p-curve:
/// `boptools_2d::curve_on_surface_range` is `FC2D_HasOldCurveOnSurface`
/// (`TopOpeBRepTool_2d.cxx:356-360` reads the `BRep_GCurve` range) and
/// `boptools_2d::make_2d` is the `HasOld`/`HasNew`/`MakeCurveOnSurface`
/// cascade of `FC2D_CurveOnSurface` (`2d.cxx:356-376`).
pub fn fun_tool_par_on_e_f_tol(
    e: &TopoShape,
    par: f64,
    f: &TopoShape,
    tol3d: f64,
) -> Option<GpPnt2d> {
    let edge = Edge(e.clone());
    let face = Face(f.clone());
    match crate::boptools_2d::curve_on_surface_range(&edge, &face) {
        None => fun_tool_par_f_tol(e, par, f, tol3d),
        Some((pc, first, last)) => {
            if !(first <= par && par <= last) {
                return None;
            }
            let p = pc.value(par);
            Some(GpPnt2d::new(p.x(), p.y()))
        }
    }
}

/// `FUN_tool_paronEF(E, par, F, UV)` (`cxx:1028-1036`): `tol3d =
/// BRep_Tool::Tolerance(F) * 1.e2`.
pub fn fun_tool_par_on_e_f(e: &TopoShape, par: f64, f: &TopoShape) -> Option<GpPnt2d> {
    let tol3d = BRepTool::face_tolerance(&Face(f.clone())) * 1.0e2;
    fun_tool_par_on_e_f_tol(e, par, f, tol3d)
}

/// `FUN_tool_findparinBAC(BAC, par)` (`cxx:1121-1126`): the interior
/// parameter `(1 - 0.34567237) * f + 0.34567237 * l`.
pub fn fun_tool_find_par_in_bac(bc: &CurveAdaptor) -> f64 {
    let fe = bc.first_parameter();
    let le = bc.last_parameter();
    let t = 0.34567237;
    (1.0 - t) * fe + t * le
}

/// `FUN_tool_findPinBAC(BAC, P, par)` (`cxx:1113-1118`): the point at
/// [`fun_tool_find_par_in_bac`]. Always `true` in OCCT.
pub fn fun_tool_find_p_in_bac(bc: &CurveAdaptor) -> Option<(GpPnt, f64)> {
    let par = fun_tool_find_par_in_bac(bc);
    Some((bc.d0(par)?, par))
}

/// `FUN_tool_findparinE(E, par)` (`cxx:1131-1136`).
pub fn fun_tool_find_par_in_e(e: &TopoShape) -> f64 {
    let bc = CurveAdaptor::of_edge(&Edge(e.clone()));
    fun_tool_find_par_in_bac(&bc)
}

/// `FUN_tool_findPinE(E, P, par)` (`cxx:1138-1143`). The C++ always answers
/// `true` (`cxx:1117`, `:1142`); `None` here means the adaptor could not
/// evaluate the curve, which OCCT would have thrown on.
pub fn fun_tool_find_p_in_e(e: &TopoShape) -> Option<(GpPnt, f64)> {
    let bc = CurveAdaptor::of_edge(&Edge(e.clone()));
    fun_tool_find_p_in_bac(&bc)
}

/// `FUN_tool_value(UV, F, P)` (`cxx:328-333`):
/// `BRepAdaptor_Surface(F).Value(UV.X(), UV.Y())`.
pub fn fun_tool_value_uv(uv: &GpPnt2d, f: &Face) -> Option<GpPnt> {
    let s = BRepTool::face_surface_world(f)?;
    Some(s.d0(uv.x(), uv.y()))
}

/// `FUN_tool_parVonE(V, E, par)` (`cxx:887-929`): the parameter of `v` on `e`.
///
/// Three arms, verbatim: (1) `v` occurs in `E` as a child vertex
/// (`IsSame`), (2) `v` shares the domain of a vertex of `E`
/// (`gp_Pnt::IsEqual(pt, Precision::Confusion())`), (3) `Extrema_ExtPC(pt,
/// BRepAdaptor_Curve(E))` and [`fun_tool_get_index`] over its solutions.
pub fn fun_tool_par_von_e(v: &TopoShape, e: &TopoShape) -> Option<f64> {
    let conf = occt_core::precision::CONFUSION;
    let mut is_v_of_e = false;
    let mut par = 0.0;
    let mut ex = Explorer::new(e, ShapeType::Vertex);
    while ex.more() {
        let cur = ex.current();
        is_v_of_e = crate::topo_tools::is_same(cur, v);
        if is_v_of_e {
            par = BRepTool::parameter_of_vertex(e, cur);
            break;
        }
        ex.next();
    }
    if !is_v_of_e {
        let pt = BRepTool::vertex_point(&Vertex(v.clone()));
        // <v> can share same domain with a vertex of <E>
        let mut ex = Explorer::new(e, ShapeType::Vertex);
        while ex.more() {
            let vex = ex.current();
            let ptex = BRepTool::vertex_point(&Vertex(vex.clone()));
            if ptex.distance(&pt) <= conf {
                return Some(BRepTool::parameter_of_vertex(e, vex));
            }
            ex.next();
        }

        let bc = CurveAdaptor::of_edge(&Edge(e.clone()));
        let curve = bc.curve()?;
        let (sols, done) = occt_geom::extrema_pc::extrema_ext_pc_range_done(
            curve.as_ref(),
            &pt,
            bc.first_parameter(),
            bc.last_parameter(),
        );
        if !done || sols.is_empty() {
            return None;
        }
        let i = super::project_curve::fun_tool_getindex(&sols)?;
        par = sols[i].u;
    }
    Some(par)
}

/// `FUN_tool_nggeomF(paronE, E, F, nggeomF, tol)` (`cxx:468-525`): the normal
/// of `F` at the uv of `paronE` on `E`.
///
/// `p2d` comes from the p-curve when `E` is a bound of `F`
/// (`FUN_tool_orientEinFFORWARD`) and that p-curve exists, otherwise from
/// `FUN_tool_projPonF`; the projection is retried through
/// `FC2D_CurveOnSurface` (`cxx:499-509`, "modified by NIZHNY-MZV") and the
/// projected parameter is rejected when `d(p(uv), p3d) > tol`
/// (`cxx:514-522`, "xpu010698").
pub fn fun_tool_ng_geom_f_par_tol(
    paron_e: f64,
    e: &TopoShape,
    f: &TopoShape,
    tol: f64,
) -> Option<GpDir> {
    let edge = Edge(e.clone());
    let face = Face(f.clone());
    let mut p2d: Option<GpPnt2d> = None;
    let mut project = true;
    let edonfa = fun_tool_orient_ein_f_forward(e, f).is_some();
    if edonfa {
        if let Some((pc, _first, _last)) = crate::boptools_2d::curve_on_surface_range(&edge, &face)
        {
            project = false;
            let p = pc.value(paron_e);
            p2d = Some(GpPnt2d::new(p.x(), p.y()));
        }
    }
    if project {
        let p3d = fun_tool_value(paron_e, e)?;
        let mut ok = true;
        let mut uv = GpPnt2d::new(0.0, 0.0);
        match fun_tool_proj_pon_f(&p3d, &face) {
            Some((u, v, _d)) => uv = GpPnt2d::new(u, v),
            None => {
                // modified by NIZHNY-MZV  Wed Dec  1 13:55:08 1999
                // if !ok try to compute new pcurve
                ok = false;
                if let Some((pc, _first, _last)) =
                    crate::boptools_2d::curve_on_surface_range(&edge, &face)
                {
                    ok = true;
                    let p = pc.value(paron_e);
                    uv = GpPnt2d::new(p.x(), p.y());
                }
            }
        }
        if !ok {
            return None;
        }
        // modified by NIZHNY-MZV  Wed Dec  1 13:56:14 1999
        // xpu010698
        let p3duv = fun_tool_value_uv(&uv, &face)?;
        if p3duv.distance(&p3d) > tol {
            return None;
        }
        p2d = Some(uv);
    }
    fun_tool_ng_geom_f(&p2d?, &face)
}

/// `FUN_tool_nggeomF(paronE, E, F, nggeomF)` (`cxx:528-537`): `tol3d =
/// BRep_Tool::Tolerance(F) * 1.e2`.
pub fn fun_tool_ng_geom_f_par(paron_e: f64, e: &TopoShape, f: &TopoShape) -> Option<GpDir> {
    let tol3d = BRepTool::face_tolerance(&Face(f.clone())) * 1.0e2;
    fun_tool_ng_geom_f_par_tol(paron_e, e, f, tol3d)
}

/// `FUN_tool_Eshared(V, F1, F2, Eshared)` (`cxx:843-883`): an edge shared by
/// `F1` and `F2` that has `V` as a bound.
///
/// The C++ walks every edge of `F1`, keeps those carrying `V` (`IsSame`) in
/// `e1s`, then returns the first edge of `F2` that `IsSame` as one of them.
pub fn fun_tool_eshared(v: &TopoShape, f1: &TopoShape, f2: &TopoShape) -> Option<TopoShape> {
    let mut e1s: Vec<TopoShape> = Vec::new();
    let mut ex = Explorer::new(f1, ShapeType::Edge);
    while ex.more() {
        let e1 = ex.current();
        let mut e1ok = false;
        let mut exv = Explorer::new(e1, ShapeType::Vertex);
        while exv.more() {
            if crate::topo_tools::is_same(exv.current(), v) {
                e1ok = true;
                break;
            }
            exv.next();
        }
        if e1ok {
            e1s.push(e1.clone());
        }
        ex.next();
    }

    let mut ex = Explorer::new(f2, ShapeType::Edge);
    while ex.more() {
        let e2 = ex.current();
        for e1 in &e1s {
            if crate::topo_tools::is_same(e1, e2) {
                return Some(e2.clone());
            }
        }
        ex.next();
    }
    None
}

/// `FUN_tool_EtgF(paronE, E, p2d, F, tola)` (`cxx:541-558`): `true` when the
/// tangent of `E` at `paronE` is orthogonal to the normal of `F` at `p2d`.
pub fn fun_tool_etg_f(
    par_on_e: f64,
    e: &TopoShape,
    p2d: &GpPnt2d,
    f: &TopoShape,
    tola: f64,
) -> bool {
    let Some(tg_e) = super::tangent::tg_geom_e(par_on_e, e) else {
        return false; // NYIRAISE
    };
    let Some(ng_f) = fun_tool_ng_geom_f(p2d, &Face(f.clone())) else {
        return false;
    };
    let prod = tg_e.dot(&GpVec::new(ng_f.x(), ng_f.y(), ng_f.z()));
    prod.abs() < tola
}

/// `FUN_tool_EtgOOE(paronE, E, paronOOE, OOE, tola)` (`cxx:561-583`): `true`
/// when `E` and `OOE` are tangent at their respective parameters.
pub fn fun_tool_etg_oo_e(
    par_on_e: f64,
    e: &TopoShape,
    par_on_oo_e: f64,
    oo_e: &TopoShape,
    tola: f64,
) -> bool {
    let Some(tg_oo_e) = super::tangent::tg_geom_e(par_on_oo_e, oo_e) else {
        return false; // NYIRAISE
    };
    let Some(tg_e) = super::tangent::tg_geom_e(par_on_e, e) else {
        return false; // NYIRAISE
    };
    let prod = tg_oo_e.dot(&tg_e);
    (1.0 - prod.abs()).abs() < tola
}

/// `FUN_tool_PinC(P, BAC, pmin, pmax, tol)` (`cxx:279-299`): `true` when the
/// closest point of `P` on `BAC[ pmin, pmax ]` is within `tol`.
pub fn fun_tool_pin_c_range(
    p: &GpPnt,
    bc: &CurveAdaptor,
    pmin: f64,
    pmax: f64,
    tol: f64,
) -> bool {
    let (sols, done) = match bc.curve() {
        Some(c) => {
            occt_geom::extrema_pc::extrema_ext_pc_range_done(c.as_ref(), p, pmin, pmax)
        }
        None => (Vec::new(), false),
    };
    if !done || sols.is_empty() {
        return false;
    }
    match crate::topopebrep::tool::project_curve::fun_tool_getindex(&sols) {
        Some(i) => sols[i].sq_dist <= tol * tol,
        None => false,
    }
}

/// `FUN_tool_staPinE(P, E, tol)` (`TopOpeBRepTool_TOPOLOGY.cxx:336-357`):
/// `Extrema_ExtPC(P, BRepAdaptor_Curve(E))`, then the closest solution within
/// `tol` answers `IN`, beyond it `OUT`; a failed search stays `UNKNOWN`
/// (`cxx:345`).
pub fn fun_tool_sta_pin_e_tol(p: &GpPnt, e: &TopoShape, tol: f64) -> State {
    let bc = CurveAdaptor::of_edge(&Edge(e.clone()));
    let (sols, done) = match bc.curve() {
        Some(c) => occt_geom::extrema_pc::extrema_ext_pc_range_done(
            c.as_ref(),
            p,
            bc.first_parameter(),
            bc.last_parameter(),
        ),
        None => (Vec::new(), false),
    };
    if !done || sols.is_empty() {
        return State::Unknown;
    }
    match crate::topopebrep::tool::project_curve::fun_tool_getindex(&sols) {
        Some(i) => {
            let d = sols[i].sq_dist.sqrt();
            if d < tol {
                State::In
            } else {
                State::Out
            }
        }
        None => State::Unknown,
    }
}

/// `FUN_tool_staPinE(P, E)` (`TopOpeBRepTool_TOPOLOGY.cxx:359-365`): the
/// overload that derives the tolerance from the edge, `* 1.e2` (`cxx:362`).
pub fn fun_tool_sta_pin_e(p: &GpPnt, e: &TopoShape) -> State {
    let tol3d = BRepTool::edge_tolerance(&Edge(e.clone())) * 1.0e2;
    fun_tool_sta_pin_e_tol(p, e, tol3d)
}

/// `FUN_tool_PinC(P, BAC, tol)` (`cxx:302-310`): the bounds of `BAC`.
pub fn fun_tool_pin_c(p: &GpPnt, bc: &CurveAdaptor, tol: f64) -> bool {
    let pmin = bc.first_parameter();
    let pmax = bc.last_parameter();
    fun_tool_pin_c_range(p, bc, pmin, pmax, tol)
}

/// `FUN_tool_maxtol(S)` (`TopOpeBRepTool_TOPOLOGY.cxx:1183-1191`): the largest
/// tolerance among `S`'s faces, edges and vertices.
///
/// OCCT calls the shape-typed overload three times on the same accumulator
/// (`cxx:1148-1178`); the port walks the three types in one pass.
pub fn fun_tool_max_tol(s: &TopoShape) -> f64 {
    let mut max_tol = 0.0f64;
    for st in [ShapeType::Face, ShapeType::Edge, ShapeType::Vertex] {
        let mut ex = Explorer::new(s, st);
        while ex.more() {
            let ss = ex.current();
            let tolss = match st {
                ShapeType::Face => BRepTool::face_tolerance(&Face(ss.clone())),
                ShapeType::Edge => BRepTool::edge_tolerance(&Edge(ss.clone())),
                ShapeType::Vertex => BRepTool::vertex_tolerance(&Vertex(ss.clone())),
                _ => 0.0,
            };
            if tolss > max_tol {
                max_tol = tolss;
            }
            ex.next();
        }
    }
    max_tol
}

/// `TopOpeBRepTool_TOOL::TolP(E, F)` (`TOOL.cxx:1935-1939`):
/// `BRepAdaptor_Curve2d(E, F).Resolution(BRep_Tool::Tolerance(E))`.
///
/// `BRepAdaptor_Curve2d(E, F)` requires the p-curve of `E` on `F`; when the
/// port cannot build it, `None` stands for the OCCT `Geom2dAdaptor_Curve`
/// raise that the C++ never catches.
pub fn tol_p(e: &TopoShape, f: &TopoShape) -> Option<f64> {
    let tol_e = BRepTool::edge_tolerance(&Edge(e.clone()));
    let pc = crate::boptools_2d::make_2d(&Edge(e.clone()), &Face(f.clone())).ok()?;
    Some(pc.resolution(tol_e))
}

/// `FUN_tool_inS(subshape, shape)` (`TopOpeBRepTool_TOPOLOGY.cxx:832-840`):
/// `shape` mapped down to `subshape`'s own type contains it.
pub fn fun_tool_in_s(subshape: &TopoShape, shape: &TopoShape) -> bool {
    let st = subshape.shape_type();
    let mut ex = Explorer::new(shape, st);
    while ex.more() {
        if crate::topo_tools::is_same(&ex.current(), subshape) {
            return true;
        }
        ex.next();
    }
    false
}

/// `FUN_tool_comparebndkole(sh1, sh2)`
/// (`TopOpeBRepTool_TOPOLOGY.cxx:1218-1282`): compare the two bounding boxes,
/// returning 1 or 2 for the shape contained in the other, 0 otherwise.
///
/// `Bnd_Box::Get` raises on a void box; the port answers `0`, the same answer
/// as the `IsOut` test that precedes it would give for an empty box.
pub fn fun_tool_comparebndkole(sh1: &TopoShape, sh2: &TopoShape) -> i32 {
    let mut bnd1 = BndBox::new();
    crate::brep_bnd_lib::add(sh1, &mut bnd1, false);
    bnd1.set_gap(0.0);
    let mut bnd2 = BndBox::new();
    crate::brep_bnd_lib::add(sh2, &mut bnd2, false);
    bnd2.set_gap(0.0);

    // cxx:1230-1233
    if bnd1.is_out_box(&bnd2) {
        return 0;
    }
    let (Some(xyz1), Some(xyz2)) = (bnd1.get(), bnd2.get()) else {
        return 0;
    };
    let xyz1 = [xyz1.0, xyz1.1, xyz1.2, xyz1.3, xyz1.4, xyz1.5];
    let xyz2 = [xyz2.0, xyz2.1, xyz2.2, xyz2.3, xyz2.4, xyz2.5];
    let tol = CONFUSION;

    // cxx:1239-1274
    let mut neq = 0;
    let mut n2sup = 0;
    for i in 0..3 {
        let d = xyz2[i] - xyz1[i];
        if d.abs() < tol {
            neq += 1;
            continue;
        }
        if d < 0.0 {
            n2sup += 1;
        }
    }
    for i in 3..6 {
        let d = xyz2[i] - xyz1[i];
        if d.abs() < tol {
            neq += 1;
            continue;
        }
        if d > 0.0 {
            n2sup += 1;
        }
    }
    if n2sup + neq != 6 {
        return 0;
    }
    if neq == 6 {
        return 0;
    }

    // cxx:1280
    if n2sup > 0 {
        2
    } else {
        1
    }
}

/// `FUN_tool_pcurveonF(F, E)` (`TopOpeBRepTool_TOPOLOGY.cxx:1347-1367`):
/// project `E`'s 3-D curve on `F` and attach the resulting p-curve to `E`.
///
/// OCCT calls `TopOpeBRepTool_CurveTool::MakePCurveOnFace` (`CurveTool.cxx:1040-1160`),
/// which is `ProjLib_ProjectedCurve` plus the UV-window shifting of
/// `CurveTool.cxx:1071-1152`; the port uses
/// [`crate::pcurve::make_pcurve_on_face`], the same projection, and writes the
/// range and the p-curve the way `BRep_Builder::UpdateEdge(E, C2d, F, Tol)`
/// does. The UV-window shift is UNPORTED.
pub fn fun_tool_pcurveon_f(f: &TopoShape, e: &TopoShape) -> bool {
    let (Some(face), Some(edge)) = (Face::wrap(f.clone()), Edge::wrap(e.clone())) else {
        return false;
    };
    // cxx:1350-1354: `BRep_Tool::Curve(E, f, l)` null.
    if BRepTool::edge_curve(&edge).is_none() {
        return false;
    }
    let Ok(c2d) = crate::pcurve::make_pcurve_on_face(&edge, &face) else {
        return false;
    };
    let reg = crate::tgeometry::GeometryRegistry::global();
    let face_key = crate::tgeometry::GeometryRegistry::shape_key(&face.0);
    let (fpar, lpar) = BRepTool::edge_parameters(&edge);
    reg.set_edge_pcurve(&edge.0, face_key, c2d);
    reg.set_pcurve_range(&edge.0, face_key, fpar, lpar);
    true
}

/// `FUN_tool_MakeWire(loE, newW)` (`TopOpeBRepTool_TOPOLOGY.cxx:1626-1638`):
/// a fresh wire holding the edges in list order.
pub fn fun_tool_make_wire(lo_e: &[TopoShape], new_w: &mut Option<TopoShape>) -> bool {
    *new_w = None;
    let mut w = TopoShape::new(ShapeType::Wire);
    let builder = crate::builder::TopoBuilder::new();
    for e in lo_e {
        builder.add(&mut w, e);
    }
    *new_w = Some(w);
    true
}
