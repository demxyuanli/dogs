//! ShapeFix_ComposeShell::SplitByLine (ShapeFix_ComposeShell.cxx:1433-1914),
//! the `ShapeFix_WireSegment&` overload.

use std::sync::Arc;

use occt_core::bnd::BndBox2d;
use occt_core::gp::{GpDir2d, GpLin2d, GpPnt2d, GpTrsf2d, GpVec2d};
use occt_core::intres2d::{IntRes2dDomain, IntRes2dPosition};
use occt_core::precision::{CONFUSION, INFINITE, PCONFUSION};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;
use occt_geom2d::line::Geom2dLine;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::boptools_2d;
use crate::pcurve_full::surface_value_of_uv;
use crate::shape::{Edge, Vertex};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::is_same;

use super::helpers::*;
use super::reshape::{apply_context, ReShape};
use super::shell::{ComposeShell, SHAPEEXTEND_FAIL4};
use super::wire_segment::WireSegment;

impl ComposeShell {
    /// ShapeFix_ComposeShell::SplitByLine (cxx:1433-1914).
    #[allow(clippy::too_many_arguments)]
    pub fn split_by_line(
        &mut self,
        wire: &mut WireSegment,
        line: &GpLin2d,
        is_cut_by_u: bool,
        cut_index: i32,
        split_line_par: &mut Vec<f64>,
        split_line_code: &mut Vec<i32>,
        split_line_vertex: &mut Vec<Vertex>,
    ) -> bool {
        // cxx:1443-1444: Geom2dAdaptor_Curve jGAC(Geom2d_Line(line)).
        let j_c2d = Geom2dLine::new(*line.position());

        // cxx:1446-1448.
        let mut int_edge_ind: Vec<i32> = Vec::new();
        let mut int_edge_par: Vec<f64> = Vec::new();
        let mut int_line_par: Vec<f64> = Vec::new();

        // cxx:1450.
        let is_nonmanifold = wire.orientation() == Orientation::Internal;

        // cxx:1452-1477: gka correction for non-manifold vertices.
        if wire.is_vertex() {
            let face = match self.face() {
                Some(f) => f.clone(),
                None => return false,
            };
            let surface = match BRepTool::face_surface(&face) {
                Some(s) => s,
                None => return false,
            };
            let a_vert = wire.get_vertex().cloned();
            let a_vert = match a_vert {
                Some(v) => v,
                None => return false,
            };
            let a_p3d = BRepTool::vertex_point(&a_vert);
            let a_p2d = surface_value_of_uv(surface.as_ref(), &a_p3d, CONFUSION);
            let (code, _dev) = point_line_position(&a_p2d, line);
            if code != IOR_UNDEF {
                return false;
            }
            let par = param_point_on_line(&a_p2d, line);
            split_line_par.push(par);
            split_line_code.push(ITP_TANG);
            let builder = TopoBuilder::new();
            let mut a_vert_new = builder.make_vertex(a_p3d, BRepTool::vertex_tolerance(&a_vert));
            a_vert_new.0.set_orientation(Orientation::Forward);
            self.context_mut().replace(&a_vert.0, &a_vert_new.0);
            split_line_vertex.push(a_vert_new.clone());
            wire.set_vertex(Some(a_vert_new));
            return true;
        }

        let nbe = wire.nb_edges();
        let face = match self.face() {
            Some(f) => f.clone(),
            None => return false,
        };

        // cxx:1483-1495: closed mode direction and half period.
        let mut closed_dir = 0i32;
        if self.closed_mode {
            if self.u_closed && line.direction().x().abs() < PCONFUSION {
                closed_dir = -1;
            } else if self.v_closed && line.direction().y().abs() < PCONFUSION {
                closed_dir = 1;
            }
        }
        let half_period = 0.5
            * (if closed_dir != 0 {
                if closed_dir < 0 { self.u_period } else { self.v_period }
            } else {
                0.0
            });

        // cxx:1499-1501.
        let mut first_code = 0i32;
        let mut prev_code = 0i32;
        let mut first_pos = GpPnt2d::new(0.0, 0.0);
        let mut prev_pos = GpPnt2d::new(0.0, 0.0);
        let mut first_dev = 0.0f64;
        let mut prev_dev = 0.0f64;

        for iedge in 1..=nbe {
            let e = wire.edge(iedge).expect("edge").clone();
            let is_reversed = e.0.orientation() == Orientation::Reversed;

            // cxx:1507-1512: sae.PCurve(E, myFace, c2d, f, l, false).
            let (mut c2d, f, l): (Arc<dyn Curve2d>, f64, f64) =
                match boptools_2d::curve_on_surface_range(&e, &face) {
                    Some(v) => v,
                    None => continue,
                };

            // cxx:1516-1517.
            let mut posf = c2d.d0(f);
            let mut posl = c2d.d0(l);
            let mut pppf = posf;
            let mut pppl = posl;

            // cxx:1524-1575: ClosedMode period adjustment.
            let mut nb_iter = 1i32;
            let mut shift_next = GpVec2d::new(0.0, 0.0);
            if self.closed_mode {
                let mut box2d = BndBox2d::new();
                fill_bnd_box(c2d.as_ref(), f, l, 41, true, &mut box2d);
                let (umin, vmin, umax, vmax) = box2d.get().unwrap_or((0.0, 0.0, 0.0, 0.0));
                if closed_dir < 0 {
                    let x = line.location().x();
                    let mut shift = adjust_to_period(umin, x - self.u_period, x);
                    if shift != 0.0 {
                        let mut nc = c2d.clone_dyn();
                        let mut t = GpTrsf2d::identity();
                        t.set_translation_vec(&GpVec2d::new(shift, 0.0));
                        nc.transform(&t);
                        c2d = Arc::from(nc);
                        pppf = GpPnt2d::new(pppf.x() + shift, pppf.y());
                        pppl = GpPnt2d::new(pppl.x() + shift, pppl.y());
                    }
                    let d_umax = umax + shift - x;
                    shift_next = GpVec2d::new(if d_umax > 0.0 { -self.u_period } else { self.u_period }, 0.0);
                    nb_iter = (1.0 + d_umax.abs() / self.u_period) as i32;
                    shift = crate::shhealing::adjust_by_period(posf.x(), x, self.u_period);
                    posf = GpPnt2d::new(posf.x() + shift, posf.y());
                    shift = crate::shhealing::adjust_by_period(posl.x(), x, self.u_period);
                    posl = GpPnt2d::new(posl.x() + shift, posl.y());
                } else if closed_dir > 0 {
                    let y = line.location().y();
                    let mut shift = adjust_to_period(vmin, y - self.v_period, y);
                    if shift != 0.0 {
                        let mut nc = c2d.clone_dyn();
                        let mut t = GpTrsf2d::identity();
                        t.set_translation_vec(&GpVec2d::new(0.0, shift));
                        nc.transform(&t);
                        c2d = Arc::from(nc);
                        pppf = GpPnt2d::new(pppf.x(), pppf.y() + shift);
                        pppl = GpPnt2d::new(pppl.x(), pppl.y() + shift);
                    }
                    let d_vmax = vmax + shift - y;
                    shift_next = GpVec2d::new(0.0, if d_vmax > 0.0 { -self.v_period } else { self.v_period });
                    nb_iter = (1.0 + d_vmax.abs() / self.v_period) as i32;
                    shift = crate::shhealing::adjust_by_period(posf.y(), y, self.v_period);
                    posf = GpPnt2d::new(posf.x(), posf.y() + shift);
                    shift = crate::shhealing::adjust_by_period(posl.y(), y, self.v_period);
                    posl = GpPnt2d::new(posl.x(), posl.y() + shift);
                }
            }

            // cxx:1577-1602.
            let mut pos = if is_reversed { posl } else { posf };
            let (code, dev) = point_line_position(&pos, line);
            if iedge == 1 {
                first_code = code;
                first_pos = pos;
                first_dev = dev;
            } else if code == IOR_UNDEF || code != prev_code {
                if closed_dir == 0 || (dev - prev_dev).abs() < half_period {
                    int_line_par.push(param_points_on_line(&pos, &prev_pos, line));
                    int_edge_par.push(if is_reversed { l } else { f });
                    int_edge_ind.push(iedge as i32);
                }
            }
            pos = if is_reversed { posf } else { posl };
            let (pc, pd) = point_line_position(&pos, line);
            prev_code = pc;
            prev_dev = pd;
            prev_pos = pos;

            // cxx:1605-1656.
            for iter in 1..=nb_iter {
                let i_dom = IntRes2dDomain::bounded(&pppf, f, TOLINT, &pppl, l, TOLINT);
                let mut inter = Geom2dIntGInter::new();
                inter.perform_with_d2(&j_c2d, c2d.as_ref(), &i_dom, TOLINT, TOLINT);
                if inter.is_done() {
                    for i in 1..=inter.nb_points() {
                        let ip = inter.point(i);
                        if ip.transition_of_second().position_on_curve() == IntRes2dPosition::Middle
                            || (code != IOR_UNDEF && prev_code != IOR_UNDEF)
                        {
                            int_line_par.push(ip.param_on_first());
                            int_edge_par.push(ip.param_on_second());
                        }
                    }
                    for i in 1..=inter.nb_segments() {
                        let seg = inter.segment(i);
                        if seg.has_first_point() {
                            int_line_par.push(seg.first_point().param_on_first());
                            int_edge_par.push(seg.first_point().param_on_second());
                        }
                        if seg.has_last_point() {
                            int_line_par.push(seg.last_point().param_on_first());
                            int_edge_par.push(seg.last_point().param_on_second());
                        }
                    }
                }
                if iter < nb_iter {
                    if iter == 1 {
                        c2d = Arc::from(c2d.clone_dyn());
                    }
                    pppf = GpPnt2d::new(pppf.x() + shift_next.x(), pppf.y() + shift_next.y());
                    pppl = GpPnt2d::new(pppl.x() + shift_next.x(), pppl.y() + shift_next.y());
                    // c2d->Translate(shiftNext).
                    let mut nc = c2d.clone_dyn();
                    let mut t = GpTrsf2d::identity();
                    t.set_translation_vec(&shift_next);
                    nc.transform(&t);
                    c2d = Arc::from(nc);
                }
            }

            // cxx:1658-1672: clamp parameters into [f, l].
            let start = int_edge_ind.len();
            for i in start..int_edge_par.len() {
                if int_edge_par[i] < f {
                    int_edge_par[i] = f;
                } else if int_edge_par[i] > l {
                    int_edge_par[i] = l;
                }
            }

            // cxx:1674-1686: sort by parameter on edge.
            let mut i = int_edge_par.len();
            while i > start + 1 {
                let mut j = start + 1;
                while j < i {
                    if is_reversed == (int_edge_par[j] < int_edge_par[j - 1]) {
                        j += 1;
                        continue;
                    }
                    int_line_par.swap(j - 1, j);
                    int_edge_par.swap(j - 1, j);
                    j += 1;
                }
                i -= 1;
            }

            // cxx:1688-1692.
            for _i in start..int_edge_par.len() {
                int_edge_ind.push(iedge as i32);
            }

            // cxx:1694-1706.
            if iedge == nbe
                && wire.orientation() != Orientation::External
                && wire.orientation() != Orientation::Internal
                && (prev_code == IOR_UNDEF || prev_code != first_code)
            {
                if closed_dir == 0 || (first_dev - prev_dev).abs() < half_period {
                    int_line_par.push(param_points_on_line(&pos, &first_pos, line));
                    int_edge_par.push(if is_reversed { f } else { l });
                    int_edge_ind.push(iedge as i32);
                }
            }
        }

        if int_edge_par.is_empty() {
            return false; // cxx:1709-1714
        }

        // cxx:1718-1719.
        let mut int_code: Vec<i32> = Vec::new();
        let mut segment_codes: Vec<i32> = Vec::new();

        // cxx:1722-1767: remove duplicated points in closed mode.
        if self.closed_mode && int_edge_par.len() > 1 {
            // `ShapeFix_ComposeShell.cxx:1722`: `int j = IntEdgePar.Length();`
            // with **1-based** `NCollection_Sequence` indexing; the port's
            // `Vec` is 0-based, so the last index is `len - 1` (using
            // `len` made `int_edge_ind[j]` read past the end below).
            let mut j = int_edge_par.len() - 1;
            let mut i = 0usize;
            while i < int_edge_par.len() {
                if i == j {
                    break;
                }
                let remove = if int_edge_ind[i] == int_edge_ind[j]
                    && (int_edge_par[i] - int_edge_par[j]).abs() < PCONFUSION
                {
                    true
                } else if nbe == 1 || int_edge_ind[i] == (int_edge_ind[j] % nbe as i32) + 1 {
                    let e1 = wire.edge(int_edge_ind[j] as usize).expect("e1").clone();
                    let e2 = wire.edge(int_edge_ind[i] as usize).expect("e2").clone();
                    let (a1, b1) = match boptools_2d::curve_on_surface_range(&e1, &face) {
                        Some((_, x, y)) => (x, y),
                        None => (0.0, 0.0),
                    };
                    let (a2, b2) = match boptools_2d::curve_on_surface_range(&e2, &face) {
                        Some((_, x, y)) => (x, y),
                        None => (0.0, 0.0),
                    };
                    let e1_ref = if e1.0.orientation() == Orientation::Forward { b1 } else { a1 };
                    let e2_ref = if e2.0.orientation() == Orientation::Forward { a2 } else { b2 };
                    (int_edge_par[j] - e1_ref).abs() < PCONFUSION
                        && (int_edge_par[i] - e2_ref).abs() < PCONFUSION
                } else {
                    false
                };
                if remove {
                    int_line_par.remove(i);
                    int_edge_par.remove(i);
                    int_edge_ind.remove(i);
                    if j > i {
                        j -= 1;
                    }
                    continue;
                }
                j = i;
                i += 1;
            }
        }

        // cxx:1771-1782: segment codes.
        let n_par = int_edge_par.len();
        if n_par == 0 {
            return false;
        }
        for i in 1..=n_par {
            let j = if i < n_par { i + 1 } else { 1 };
            let code = self.compute_code(
                wire.edges(),
                line,
                int_edge_ind[i - 1] as usize,
                int_edge_ind[j - 1] as usize,
                int_edge_par[i - 1],
                int_edge_par[j - 1],
                is_nonmanifold,
            );
            segment_codes.push(code);
        }

        let mut a_new_seg_codes: Vec<i32> = Vec::new();

        if wire.orientation() == Orientation::External {
            for i in 1..=n_par {
                int_code.push(ITP_TANG | IOR_BOTH);
                a_new_seg_codes.push(segment_codes[i - 1]);
            }
        } else {
            if wire.orientation() != Orientation::Internal {
                // cxx:1799-1822.
                let mut i = 1usize;
                while i <= int_edge_par.len() {
                    let j = if i > 1 { i - 1 } else { int_edge_par.len() };
                    let k = if i < int_edge_par.len() { i + 1 } else { 1 };
                    if segment_codes[j - 1] == IOR_UNDEF && segment_codes[i - 1] == IOR_UNDEF {
                        if self.closed_mode
                            && (int_line_par[i - 1] - int_line_par[j - 1])
                                * (int_line_par[k - 1] - int_line_par[i - 1])
                                <= 0.0
                        {
                            i += 1;
                            continue;
                        }
                        int_edge_ind.remove(i - 1);
                        int_edge_par.remove(i - 1);
                        int_line_par.remove(i - 1);
                        segment_codes.remove(i - 1);
                        // `cxx:1820` is `i--` inside `for (i = 1; i <= IntEdgePar.Length(); i++)`:
                        // the implicit `i++` cancels it, so `i` stays put and the element that
                        // shifted into position `i` is re-checked. The port's `while` has no
                        // implicit increment, so `i` must simply not move here.
                        continue;
                    }
                    i += 1;
                }
            }
            if int_edge_par.is_empty() {
                return false; // cxx:1825-1828
            }

            // cxx:1834-1897.
            let mut j = int_edge_par.len();
            for i in 1..=int_edge_par.len() {
                let mut codej = segment_codes[j - 1];
                let mut codei = segment_codes[i - 1];
                if self.closed_mode {
                    if (codej & IOR_BOTH) == IOR_BOTH {
                        codej = if codej & IOR_POS != 0 { IOR_RIGHT } else { IOR_LEFT };
                    }
                    if (codei & IOR_BOTH) == IOR_BOTH {
                        codei = if codei & IOR_POS != 0 { IOR_LEFT } else { IOR_RIGHT };
                    }
                    a_new_seg_codes.push(codei);
                    if int_edge_ind[i - 1] == int_edge_ind[j - 1] {
                        a_new_seg_codes.push(codej);
                    }
                } else {
                    a_new_seg_codes.push(codei);
                }
                let mut ipcode = codej | codei;
                if codej == IOR_UNDEF {
                    if int_line_par[i - 1] > int_line_par[j - 1] {
                        ipcode |= ITP_ENDSEG;
                    } else {
                        ipcode |= ITP_BEGSEG;
                    }
                } else if codei == IOR_UNDEF {
                    let k = if i < int_line_par.len() { i + 1 } else { 1 };
                    if int_line_par[k - 1] > int_line_par[i - 1] {
                        ipcode |= ITP_BEGSEG;
                    } else {
                        ipcode |= ITP_ENDSEG;
                    }
                } else if i == j {
                    ipcode |= if (ipcode & IOR_BOTH) == IOR_BOTH && !is_nonmanifold {
                        ITP_INTER
                    } else {
                        ITP_TANG
                    };
                } else if codei == codej || is_nonmanifold {
                    ipcode |= ITP_TANG;
                } else {
                    ipcode |= ITP_INTER;
                }
                int_code.push(ipcode);
                j = i;
            }
        }

        // cxx:1900-1903.
        let mut int_vertices: Vec<Vertex> = Vec::new();
        let mut indexes = int_edge_ind.clone();
        *wire = self.split_wire(
            wire,
            &mut indexes,
            &int_edge_par,
            &mut int_vertices,
            &a_new_seg_codes,
            is_cut_by_u,
            cut_index,
        );

        // cxx:1906-1911.
        for i in 1..=int_line_par.len() {
            split_line_par.push(int_line_par[i - 1]);
            split_line_code.push(int_code[i - 1]);
            split_line_vertex.push(int_vertices[i - 1].clone());
        }


        let _ = (first_dev, first_pos, GpDir2d::new(0.0, 1.0), BndBox2d::new());
        true
    }
    /// ShapeFix_ComposeShell::SplitByLine (cxx:1918-2127), the sequence
    /// overload: split every wire, sort the split points along the line, merge
    /// null-length tangential segments, then walk the line and create one
    /// external wire segment per interior span.
    pub fn split_by_line_wires(
        &mut self,
        wires: &mut Vec<WireSegment>,
        line: &GpLin2d,
        is_cut_by_u: bool,
        cut_index: i32,
    ) {
        let mut split_line_par: Vec<f64> = Vec::new();
        let mut split_line_code: Vec<i32> = Vec::new();
        let mut split_line_vertex: Vec<Vertex> = Vec::new();

        // cxx:1929-1932.
        for i in 0..wires.len() {
            let mut w = wires[i].clone();
            self.split_by_line(
                &mut w,
                line,
                is_cut_by_u,
                cut_index,
                &mut split_line_par,
                &mut split_line_code,
                &mut split_line_vertex,
            );
            wires[i] = w;
        }

        // cxx:1935-1946: sort along the cutting line.
        let mut i = split_line_par.len();
        while i > 1 {
            for j in 1..i {
                if split_line_par[j - 1] > split_line_par[j] {
                    split_line_par.swap(j - 1, j);
                    split_line_code.swap(j - 1, j);
                    split_line_vertex.swap(j - 1, j);
                }
            }
            i -= 1;
        }

        // cxx:1949-1965: merge null-length tangential segments.
        let mut i = 1usize;
        while i < split_line_par.len() {
            if (split_line_par[i] - split_line_par[i - 1]).abs() > PCONFUSION
                && !is_same(&split_line_vertex[i - 1].0, &split_line_vertex[i].0)
            {
                i += 1;
                continue;
            }
            if (split_line_code[i - 1] & ITP_ENDSEG != 0 && split_line_code[i] & ITP_BEGSEG != 0)
                || (split_line_code[i - 1] & ITP_BEGSEG != 0 && split_line_code[i] & ITP_ENDSEG != 0)
            {
                let code = (split_line_code[i - 1] | split_line_code[i]) & IOR_BOTH;
                split_line_code[i - 1] =
                    code | (if code == IOR_BOTH { ITP_INTER } else { ITP_TANG });
                split_line_par.remove(i);
                split_line_code.remove(i);
                split_line_vertex.remove(i);
                // `ShapeFix_ComposeShell.cxx:1949` is a `for (i = 1; i <
                // Length(); i++)`: OCCT advances `i` after the removal too, so
                // the next comparison is between the elements that were at
                // `i+1`/`i+2` before the erase. The port's `while` must do the
                // same or it compares one extra (already merged) pair.
                i += 1;
            } else {
                i += 1;
            }
        }

        // cxx:1969-2110.
        let builder = TopoBuilder::new();
        let reg = GeometryRegistry::global();
        let face = match self.face() {
            Some(f) => f.clone(),
            None => return,
        };
        let face_key = GeometryRegistry::shape_key(&face.0);
        let mut parity = 0i32;
        let mut halfparity = 0i32;
        let mut tanglevel = 0i32;
        for i in 1..=split_line_par.len() {
            let code = split_line_code[i - 1];
            let interior = tanglevel == 0 && parity % 2 != 0;
            if code & ITP_INTER != 0 {
                parity += 1;
            } else if code & ITP_BEGSEG != 0 {
                tanglevel += 1;
                if halfparity == 0 {
                    halfparity = code & IOR_BOTH;
                } else if halfparity != (code & IOR_BOTH) {
                    parity += 1;
                }
            } else if code & ITP_ENDSEG != 0 {
                tanglevel -= 1;
                if halfparity == 0 {
                    halfparity = code & IOR_BOTH;
                } else if halfparity != (code & IOR_BOTH) {
                    parity += 1;
                }
            }
            if tanglevel < 0 {
                // cxx:2005-2011: OCCT only prints under OCCT_DEBUG.
            }
            if !interior {
                continue;
            }

            // cxx:2019-2057.
            let tmp_v1 = self.context().apply(&split_line_vertex[i - 2].0);
            let tmp_v2 = self.context().apply(&split_line_vertex[i - 1].0);
            let mut v1 = Vertex(tmp_v1);
            let mut v2 = Vertex(tmp_v2);
            let can_be_merged = i - 1 > 1 || i < split_line_par.len();
            let mut a_max_tol = self.max_tolerance();
            if a_max_tol <= 2.0 * CONFUSION {
                a_max_tol = INFINITE;
            }
            let a_tol1 = BRepTool::vertex_tolerance(&v1).min(a_max_tol);
            let a_tol2 = BRepTool::vertex_tolerance(&v2).min(a_max_tol);
            let a_d = BRepTool::vertex_point(&v1).square_distance(&BRepTool::vertex_point(&v2));
            if split_line_par[i - 1] - split_line_par[i - 2] < PCONFUSION
                || (can_be_merged && (a_d <= a_tol1 * a_tol1 || a_d <= a_tol2 * a_tol2))
            {
                if !is_same(&v1.0, &v2.0) {
                    let v = crate::shhealing::combine_vertex(&v1, &v2, 1.0001);
                    let mut r1 = v.clone();
                    r1.0.set_orientation(v1.0.orientation());
                    let mut r2 = v.clone();
                    r2.0.set_orientation(v2.0.orientation());
                    self.context_mut().replace(&v1.0, &r1.0);
                    self.context_mut().replace(&v2.0, &r2.0);
                    v1 = v.clone();
                    v2 = v;
                }
                continue;
            }

            // cxx:2061-2070: a new edge carrying only the two pcurves on myFace.
            let mut edge = builder.make_shape(ShapeType::Edge);
            v1.0.set_orientation(Orientation::Forward);
            v2.0.set_orientation(Orientation::Reversed);
            builder.add(&mut edge, &v1.0);
            builder.add(&mut edge, &v2.0);
            let lin1: Arc<dyn Curve2d> = Arc::new(Geom2dLine::new(*line.position()));
            let lin2: Arc<dyn Curve2d> = Arc::new(Geom2dLine::new(*line.position()));
            reg.set_edge_pcurves(&edge, face_key, vec![lin1, lin2]);
            reg.set_pcurve_range(&edge, face_key, split_line_par[i - 2], split_line_par[i - 1]);
            reg.set_edge_range(&edge, split_line_par[i - 2], split_line_par[i - 1]);


            // cxx:2072-2074.
            let mut seg = WireSegment::with_edges(vec![Edge(edge)], Orientation::External);
            define_patch(&mut seg, IOR_UNDEF, is_cut_by_u, cut_index, 0);
            // cxx:2078-2108.
            if !is_cut_by_u {
                let shift_u = if self.closed_mode && self.u_closed {
                    adjust_to_period(
                        split_line_par[i - 2] - TOLINT,
                        self.grid.u_joint_value(1),
                        self.grid.u_joint_value(2),
                    )
                } else {
                    0.0
                };
                let a_par = split_line_par[i - 2] + shift_u;
                seg.define_iu_min(
                    1,
                    get_patch_index(a_par + PCONFUSION, self.grid.u_joint_values(), self.u_closed),
                );
                seg.define_iu_max(
                    1,
                    get_patch_index(a_par - PCONFUSION, self.grid.u_joint_values(), self.u_closed) + 1,
                );
            } else {
                let shift_v = if self.closed_mode && self.v_closed {
                    adjust_to_period(
                        split_line_par[i - 2] - TOLINT,
                        self.grid.v_joint_value(1),
                        self.grid.v_joint_value(2),
                    )
                } else {
                    0.0
                };
                let a_par = split_line_par[i - 2] + shift_v;
                seg.define_iv_min(
                    1,
                    get_patch_index(a_par + PCONFUSION, self.grid.v_joint_values(), self.v_closed),
                );
                seg.define_iv_max(
                    1,
                    get_patch_index(a_par - PCONFUSION, self.grid.v_joint_values(), self.v_closed) + 1,
                );
            }
            wires.push(seg);
        }
        if parity % 2 != 0 {
            self.status |= SHAPEEXTEND_FAIL4; // cxx:2111-2113
        }

        // cxx:2120-2126.
        for w in wires.iter_mut() {
            let mut j = 1usize;
            while j <= w.nb_edges() {
                j += apply_context(w, j, self.context()) as usize;
            }
        }
    }
}