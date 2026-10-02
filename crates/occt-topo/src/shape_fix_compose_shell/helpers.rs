//! ShapeFix_ComposeShell file-static helpers
//! (ShapeFix_ComposeShell.cxx:290-303 constants, 308-495 geometry + patch
//! index, 876-938 check/patch/grid helpers).

use occt_core::gp::{GpLin2d, GpPnt2d, GpVec2d};
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face};
use super::wire_segment::WireSegment;

/// ShapeFix_ComposeShell.cxx:290-303.
pub const TOLINT: f64 = 1.0e-10;
pub const IOR_UNDEF: i32 = 0;
pub const IOR_LEFT: i32 = 1;
pub const IOR_RIGHT: i32 = 2;
pub const IOR_BOTH: i32 = 3;
pub const IOR_POS: i32 = 4;
pub const ITP_INTER: i32 = 8;
pub const ITP_BEGSEG: i32 = 16;
pub const ITP_ENDSEG: i32 = 32;
pub const ITP_TANG: i32 = 64;

/// PointLineDeviation (cxx:308-313): signed distance of p from the line.
pub fn point_line_deviation(p: &GpPnt2d, line: &GpLin2d) -> f64 {
    let dir = line.direction();
    let n = GpVec2d::new(-dir.y(), dir.x());
    let loc = line.location();
    let v = GpVec2d::new(p.x() - loc.x(), p.y() - loc.y());
    n.dot(&v)
}

/// PointLinePosition (cxx:318-322): IOR_LEFT / IOR_RIGHT / IOR_UNDEF.
pub fn point_line_position(p: &GpPnt2d, line: &GpLin2d) -> (i32, f64) {
    let dev = point_line_deviation(p, line);
    let pos = if dev > TOLINT {
        IOR_LEFT
    } else if dev < -TOLINT {
        IOR_RIGHT
    } else {
        IOR_UNDEF
    };
    (pos, dev)
}

/// ParamPointOnLine (cxx:336-339).
pub fn param_point_on_line(p: &GpPnt2d, line: &GpLin2d) -> f64 {
    let dir = line.direction();
    let loc = line.location();
    GpVec2d::new(dir.x(), dir.y()).dot(&GpVec2d::new(p.x() - loc.x(), p.y() - loc.y()))
}

/// ParamPointsOnLine (cxx:344-369): parameter of the intersection of segment
/// p1 p2 with the line (or the on-line endpoint).
pub fn param_points_on_line(p1: &GpPnt2d, p2: &GpPnt2d, line: &GpLin2d) -> f64 {
    let dist1 = point_line_deviation(p1, line);
    let dist2 = point_line_deviation(p2, line);
    let pconf = occt_core::precision::PCONFUSION;
    if dist1.abs() < pconf {
        if dist2.abs() < pconf {
            return 0.5 * (param_point_on_line(p1, line) + param_point_on_line(p2, line));
        }
        return param_point_on_line(p1, line);
    }
    if dist2.abs() < pconf {
        return param_point_on_line(p2, line);
    }
    if dist2 * dist1 > 0.0 {
        return 0.5 * (param_point_on_line(p1, line) + param_point_on_line(p2, line));
    }
    (param_point_on_line(p1, line) * dist2 - param_point_on_line(p2, line) * dist1)
        / (dist2 - dist1)
}

/// ProjectPointOnLine (cxx:374-377).
pub fn project_point_on_line(p: &GpPnt2d, line: &GpLin2d) -> GpPnt2d {
    let t = param_point_on_line(p, line);
    let loc = line.location();
    let dir = line.direction();
    GpPnt2d::new(loc.x() + dir.x() * t, loc.y() + dir.y() * t)
}

/// IsCoincided (cxx:451-462).
pub fn is_coincided(p1: &GpPnt2d, p2: &GpPnt2d, u_resolution: f64, v_resolution: f64, tol: f64) -> bool {
    let u_tol = u_resolution * tol;
    let v_tol = v_resolution * tol;
    (p1.x() - p2.x()).abs() <= TOLINT.max(u_tol) && (p1.y() - p2.y()).abs() <= TOLINT.max(v_tol)
}

/// ShapeAnalysis::AdjustToPeriod (ShapeAnalysis.cxx:66-69).
pub fn adjust_to_period(val: f64, val_min: f64, val_max: f64) -> f64 {
    crate::shhealing::adjust_by_period(val, 0.5 * (val_min + val_max), val_max - val_min)
}

/// GetPatchIndex (cxx:467-495): 1-based index of the patch holding Param.
pub fn get_patch_index(param: f64, params: &[f64], is_closed: bool) -> i32 {
    let np = params.len() as i32;
    if np < 2 {
        return 0;
    }
    let period = params[(np - 1) as usize] - params[0];
    let shift = if is_closed {
        adjust_to_period(param, params[0], params[(np - 1) as usize])
    } else {
        0.0
    };
    let p = param + shift;
    let mut i = 2i32;
    while i < np {
        if p < params[(i - 1) as usize] {
            break;
        }
        i += 1;
    }
    i -= 1;
    let ish = shift / period;
    let ishift = if ish < 0.0 {
        (ish - 0.5) as i32
    } else {
        (ish + 0.5) as i32
    };
    i - ishift * (np - 1)
}

/// DefinePatch (cxx:896-925).
pub fn define_patch(
    wire: &mut WireSegment,
    code: i32,
    is_cut_by_u: bool,
    cut_index: i32,
    number: i32,
) {
    let nb = if number > 0 {
        number as usize
    } else {
        wire.nb_edges()
    };
    if is_cut_by_u {
        if code & IOR_LEFT == 0 {
            wire.define_iu_min(nb, cut_index);
        }
        if code & IOR_RIGHT == 0 {
            wire.define_iu_max(nb, cut_index);
        }
    } else {
        if code & IOR_RIGHT == 0 {
            wire.define_iv_min(nb, cut_index);
        }
        if code & IOR_LEFT == 0 {
            wire.define_iv_max(nb, cut_index);
        }
    }
}

/// GetGridResolution (cxx:929-938).
pub fn get_grid_resolution(split_values: &[f64], cut_index: i32) -> f64 {
    let nb = split_values.len() as i32;
    let ci = cut_index as usize;
    let left_len = if cut_index > 1 {
        split_values[ci - 1] - split_values[ci - 2]
    } else {
        split_values[(nb - 1) as usize] - split_values[(nb - 2) as usize]
    };
    let right_len = if cut_index < nb {
        split_values[ci] - split_values[ci - 1]
    } else {
        split_values[1] - split_values[0]
    };
    left_len.min(right_len) / 3.0
}
/// CheckByCurve3d (cxx:876-892): the 3d point pos must lie within tol of c3d
/// at param; the transform t is applied when it is not the identity.
pub fn check_by_curve_3d(
    pos: &occt_core::gp::GpPnt,
    c3d: Option<&std::sync::Arc<dyn occt_geom::Curve>>,
    param: f64,
    t: &occt_core::gp::GpTrsf,
    tol: f64,
) -> bool {
    let Some(c) = c3d else {
        return true;
    };
    let mut pt = c.d0(param);
    if t.form() != occt_core::gp::TrsfForm::Identity {
        pt = pt.transformed(t);
    }
    // pos.SquareDistance(p) <= tol * tol (cxx:891).
    pos.distance(&pt).powi(2) <= tol * tol
}
/// DistributeSplitPoints (cxx:838-872): after the context replaced edge
/// `index` by `nsplit` edges, shift the split indices of that edge so they
/// keep pointing at the same sub-edges. `indexes` and `values` are 1-based in
/// OCCT; here they are 0-based slices in the same order.
pub fn distribute_split_points(
    edges: &[Edge],
    index: usize,
    nsplit: usize,
    indexes: &mut [i32],
    values: &[f64],
) {
    let isreversed = nsplit > 0
        && edges[index - 1].0.orientation() == crate::abs::Orientation::Reversed;
    let reg = crate::tgeometry::GeometryRegistry::global();
    let mut params = vec![0.0f64; nsplit];
    for i in 0..nsplit {
        let (f, l) = reg.edge_parameters(&edges[index + i - 1].0);
        params[i] = if isreversed { l } else { f };
    }
    let mut i = 0usize;
    while i < indexes.len() && (indexes[i] as usize) < index {
        i += 1;
    }
    let mut shift = 1usize;
    while i < indexes.len() && indexes[i] as usize == index {
        while shift < nsplit && isreversed != (values[i] > params[shift]) {
            shift += 1;
        }
        indexes[i] = (index + shift - 1) as i32;
        i += 1;
    }
    while i < indexes.len() {
        indexes[i] += (nsplit - 1) as i32;
        i += 1;
    }
}
/// ShapeAnalysis_Curve::SearchForExtremum (ShapeAnalysis_Curve.cxx:741-786):
/// Newton search for an extremum of the dir-projected coordinates of c2d in
/// [first, last]. `par` is the seed and receives the found parameter; `res`
/// receives the point. Returns false when the Newton step escapes the
/// interval more than twice.
pub fn search_for_extremum(
    c2d: &dyn occt_geom2d::curve::Curve2d,
    first: f64,
    last: f64,
    dir: &GpVec2d,
    par: &mut f64,
    res: &mut GpPnt2d,
) -> bool {
    let mut nb_out = 0i32;
    for _ in 0..10 {
        let prev_par = *par;
        let (p, d1, d2) = c2d.d2(*par);
        *res = p;
        let det = d2.dot(dir);
        if det.abs() < 1e-10 {
            return true;
        }
        *par -= d1.dot(dir) / det;
        if (*par - prev_par).abs() < PCONFUSION {
            return true;
        }
        if *par < first {
            let escaped = nb_out > 2 || prev_par == first; // cxx:770: nbOut++ > 2
            nb_out += 1;
            if escaped {
                return false;
            }
            *par = first;
        }
        if *par > last {
            let escaped = nb_out > 2 || prev_par == last; // cxx:778
            nb_out += 1;
            if escaped {
                return false;
            }
            *par = last;
        }
    }
    true
}

/// ShapeAnalysis_Curve::FillBndBox (ShapeAnalysis_Curve.cxx:788-848). Both
/// arms are ported: Exact=false samples NPoints (cxx:795-805), Exact=true
/// walks the C2 intervals and adds the X/Y extrema found by
/// SearchForExtremum (cxx:808-847).
pub fn fill_bnd_box(
    c2d: &dyn occt_geom2d::curve::Curve2d,
    first: f64,
    last: f64,
    npoints: i32,
    exact: bool,
    box2d: &mut occt_core::bnd::BndBox2d,
) {
    if !exact {
        let nseg = if npoints < 2 { 1 } else { npoints - 1 };
        let step = (last - first) / nseg as f64;
        for i in 0..=nseg {
            let par = first + i as f64 * step;
            box2d.add_point(&c2d.d0(par));
        }
        return;
    }
    // cxx:809-825: Geom2dAdaptor_Curve anAC(C2d, First, Last) and its
    // NbIntervals(GeomAbs_C2). GeomAbs_C2 == 4 (GeomAbs_Shape.hxx).
    const GEOMABS_C2: u8 = 4;
    let nb_int = c2d.nb_intervals(GEOMABS_C2);
    let nb_samples = if nb_int < 2 { npoints - 1 } else { nb_int };
    let params: Vec<f64> = if nb_samples == nb_int {
        c2d.parameter_intervals(GEOMABS_C2)
    } else {
        let step = (last - first) / nb_samples as f64;
        (0..=nb_samples).map(|i| first + f64::from(i) * step).collect()
    };
    for i in 0..=nb_samples as usize {
        let a_par1 = params[i];
        box2d.add_point(&c2d.d0(a_par1));
        if i < nb_samples as usize {
            let a_par2 = params[i + 1];
            let par = 0.5 * (a_par1 + a_par2);
            let mut pextr = GpPnt2d::new(0.0, 0.0);
            let mut parextr = par;
            if search_for_extremum(c2d, a_par1, a_par2, &GpVec2d::new(1.0, 0.0), &mut parextr, &mut pextr) {
                box2d.add_point(&pextr);
            }
            parextr = par;
            if search_for_extremum(c2d, a_par1, a_par2, &GpVec2d::new(0.0, 1.0), &mut parextr, &mut pextr) {
                box2d.add_point(&pextr);
            }
        }
    }
}

/// GetMiddlePoint (cxx:2940-2974).
pub fn get_middle_point(seg: &WireSegment, face: &crate::shape::Face) -> GpPnt2d {
    if seg.is_vertex() {
        // cxx:2942-2949: ShapeAnalysis_Surface::ValueOfUV of the vertex point.
        let v = seg.get_vertex().expect("vertex segment without vertex");
        let p3d = BRepTool::vertex_point(v);
        if let Some(surf) = BRepTool::face_surface(face) {
            return crate::pcurve_full::surface_value_of_uv(surf.as_ref(), &p3d, CONFUSION);
        }
        return GpPnt2d::new(0.0, 0.0);
    }
    let mut box2d = occt_core::bnd::BndBox2d::new();
    for e in seg.edges() {
        if let Some((c2d, cf, cl)) = crate::boptools_2d::curve_on_surface_oriented(e, face, false) {
            fill_bnd_box(c2d.as_ref(), cf, cl, 3, false, &mut box2d);
        }
    }
    match box2d.get() {
        Some((xmin, ymin, xmax, ymax)) => {
            GpPnt2d::new(0.5 * (xmax + xmin), 0.5 * (ymax + ymin))
        }
        None => GpPnt2d::new(0.0, 0.0),
    }
}
/// IsShortSegment (cxx:2394-2447): 1 for a closed segment whose interior
/// collapses to its vertex, -1 when only the 2d check fails, 0 otherwise.
pub fn is_short_segment(
    seg: &WireSegment,
    face: &crate::shape::Face,
    grid_surface: &super::composite_surface::CompositeSurface,
    u_resolution: f64,
    v_resolution: f64,
) -> i32 {
    let (Some(vf), Some(vl)) = (seg.first_vertex(), seg.last_vertex()) else {
        return 0;
    };
    if !crate::topo_tools_full::is_same(&vf.0, &vl.0) {
        return 0;
    }
    let pnt = BRepTool::vertex_point(&vf);
    let tol = BRepTool::vertex_tolerance(&vf);
    let tol2 = tol * tol;
    let mut code = 1i32;
    for edge in seg.edges() {
        // `cxx:2417`: `sae.LastVertex(edge)` is orientation aware, so a REVERSED
        // edge in the wire must yield the other topological vertex.
        let Some(last) = crate::shhealing::last_vertex(edge) else {
            return 0;
        };
        if !crate::topo_tools_full::is_same(&vf.0, &last.0) {
            return 0;
        }
        // `cxx:2423`: `sae.PCurve(edge, myFace, c2d, f, l)` keeps the default
        // `CumOri = true`.
        let Some((c2d, f, l)) =
            crate::boptools_2d::curve_on_surface_oriented(edge, face, true)
        else {
            continue;
        };
        let end_pnt = c2d.d0(l);
        let mid_pnt = c2d.d0(0.5 * (f + l));
        if !is_coincided(&end_pnt, &mid_pnt, u_resolution, v_resolution, tol) {
            code = -1;
        }
        let mid3d = grid_surface.value_uv(mid_pnt.x(), mid_pnt.y());
        // myLoc is the identity here: port shapes carry no location (cxx:2438).
        if mid3d.distance(&pnt).powi(2) > tol2 {
            return 0;
        }
    }
    code
}
/// IsSamePatch (cxx:2452-2507): true when the segment patch indices, shifted
/// onto the same period as the current ones, still span at most one patch in
/// each direction. With `extend` the caller indices are widened to the union.
pub fn is_same_patch(
    wire: &WireSegment,
    nu: i32,
    nv: i32,
    iumin: &mut i32,
    iumax: &mut i32,
    ivmin: &mut i32,
    ivmax: &mut i32,
    extend: bool,
) -> bool {
    let Some((mut jumin, mut jumax, mut jvmin, mut jvmax)) = wire.get_patch_index(1) else {
        return false;
    };
    let (mut du, mut dv) = (0i32, 0i32);
    if jumin - *iumin > nu {
        du = -(jumin - *iumin) / nu;
    } else if *iumin - jumin > nu {
        du = (*iumin - jumin) / nu;
    }
    if jvmin - *ivmin > nv {
        dv = -(jvmin - *ivmin) / nv;
    } else if *ivmin - jvmin > nv {
        dv = (*ivmin - jvmin) / nv;
    }
    if du != 0 {
        jumin += du * nu;
        jumax += du * nu;
    }
    if dv != 0 {
        jvmin += dv * nv;
        jvmax += dv * nv;
    }
    let iun = (*iumin).min(jumin);
    let iux = (*iumax).max(jumax);
    let ivn = (*ivmin).min(jvmin);
    let ivx = (*ivmax).max(jvmax);
    let ok = (iun == iux || iun + 1 == iux) && (ivn == ivx || ivn + 1 == ivx);
    if ok && extend {
        *iumin = iun;
        *iumax = iux;
        *ivmin = ivn;
        *ivmax = ivx;
    }
    ok
}

/// ShapeAnalysis_Edge::GetEndTangent2d (ShapeAnalysis_Edge.cxx:269-366): the
pub fn get_end_tangent_2d(
    edge: &Edge,
    face: &Face,
    atend: bool,
    pnt: &mut GpPnt2d,
    v: &mut GpVec2d,
    dparam: f64,
) -> bool {
    // cxx:293 PCurve(edge, S, L, c2d, cf, cl) with the default orient=true.
    let Some((c2d, cf, cl)) = crate::boptools_2d::curve_on_surface_oriented(edge, face, true)
    else {
        *v = GpVec2d::new(0.0, 0.0);
        return false;
    };
    let mut dpnew = dparam;
    if dpnew > CONFUSION {
        let delta = (cl - cf) * dpnew;
        if delta.abs() < PCONFUSION {
            dpnew = 0.0; // cxx:304-307
        } else {
            let (par1, par2) = if atend {
                (cl, cl - delta)
            } else {
                (cf, cf + delta)
            };
            *pnt = c2d.d0(par1);
            let ptmp = c2d.d0(par2);
            // cxx:316/324: atend -> pnt - ptmp, else ptmp - pnt.
            *v = if atend {
                GpVec2d::new(pnt.x() - ptmp.x(), pnt.y() - ptmp.y())
            } else {
                GpVec2d::new(ptmp.x() - pnt.x(), ptmp.y() - pnt.y())
            };
            if v.square_magnitude() < PCONFUSION * PCONFUSION {
                dpnew = 0.0; // cxx:326-329
            }
        }
    }
    if dpnew <= CONFUSION {
        // cxx:333-357: non-null tangent through the 3rd derivative, or the
        // straight chord between the endpoints.
        let par = if atend { cl } else { cf };
        let (p1, d1) = c2d.d1(par);
        *pnt = p1;
        *v = d1;
        if v.square_magnitude() < PCONFUSION * PCONFUSION {
            let (p2, _d1, d2) = c2d.d2(par);
            *pnt = p2;
            *v = d2;
            if v.square_magnitude() < PCONFUSION * PCONFUSION {
                let (p3, _d1, _d2, d3) = c2d.d3(par);
                *pnt = p3;
                *v = d3;
                if v.square_magnitude() < PCONFUSION * PCONFUSION {
                    let pe = c2d.d0(if atend { cf } else { cl });
                    *v = GpVec2d::new(pe.x() - pnt.x(), pe.y() - pnt.y());
                    if v.square_magnitude() < PCONFUSION * PCONFUSION {
                        return false;
                    }
                }
            }
        }
        if edge.0.orientation().is_reversed() {
            v.reverse(); // cxx:358-361
        }
    }
    true
}
