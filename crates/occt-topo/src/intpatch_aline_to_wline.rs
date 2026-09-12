//! `IntPatch_ALineToWLine` — sample an analytic line into a walking line.

use occt_core::precision::{CONFUSION, PCONFUSION};
use occt_geom::Surface;

use crate::geom_int::surface_parameters;
use crate::geom_int::ALine;
use crate::int_tools_wline::{PntOn2S, WLine, WLineWay};

use super::special_points::{add_cross_uv_iso_point, add_singular_pole};

/// Converter from `IntPatch_ALine` to `IntPatch_WLine`.
pub struct ALineToWLine<'a> {
    s1: &'a dyn Surface,
    s2: &'a dyn Surface,
    nb_points: i32,
    tol_open: f64,
    tol3d: f64,
}

impl<'a> ALineToWLine<'a> {
    pub fn new(s1: &'a dyn Surface, s2: &'a dyn Surface, nb_points: i32) -> Self {
        Self {
            s1,
            s2,
            nb_points: nb_points.max(3),
            tol_open: 1.0e-9,
            tol3d: CONFUSION,
        }
    }

    /// `MakeWLine` over the full ALine parameter range.
    pub fn make_wline(&self, aline: &ALine) -> Vec<WLine> {
        let mut f = aline.curve.first_parameter();
        let mut l = aline.curve.last_parameter();
        if !f.is_finite() || !l.is_finite() || (l - f).abs() < PCONFUSION {
            if aline.has_first_point && aline.has_last_point {
                f = aline.vertex(aline.first_index).parameter_on_line();
                l = aline.vertex(aline.last_index).parameter_on_line();
            } else if aline.nb_vertex() >= 2 {
                f = aline.vertex(1).parameter_on_line();
                l = aline.vertex(aline.nb_vertex()).parameter_on_line();
            } else {
                return Vec::new();
            }
        }
        if !aline.has_first_point {
            f += self.tol_open;
        }
        if !aline.has_last_point {
            l -= self.tol_open;
        }
        if l <= f {
            return Vec::new();
        }
        self.make_wline_range(aline, f, l)
    }

    fn make_wline_range(&self, aline: &ALine, f: f64, l: f64) -> Vec<WLine> {
        let n = self.nb_points as usize;
        let mut points = Vec::with_capacity(n);
        let mut prev: Option<PntOn2S> = None;
        for i in 0..n {
            let t = f + (l - f) * (i as f64) / ((n - 1) as f64);
            let p3d = aline.value(t);
            let Some((u1, v1)) = surface_parameters(self.s1, &p3d) else {
                continue;
            };
            let Some((u2, v2)) = surface_parameters(self.s2, &p3d) else {
                continue;
            };
            let mut p = PntOn2S {
                p: p3d,
                u1,
                v1,
                u2,
                v2,
            };
            if let Some(ref r) = prev {
                if let Some(pole) = add_singular_pole(self.s1, self.s2, r, &p3d, false)
                    .or_else(|| add_singular_pole(self.s2, self.s1, r, &p3d, true))
                {
                    p = pole;
                } else if let Some(seam) =
                    add_cross_uv_iso_point(self.s1, self.s2, r, self.tol3d, false).or_else(|| {
                        add_cross_uv_iso_point(self.s2, self.s1, r, self.tol3d, true)
                    })
                {
                    p = seam;
                }
            }
            prev = Some(p);
            points.push(p);
        }
        if points.len() < 2 {
            return Vec::new();
        }
        let mut wl = WLine::new();
        wl.set_creating_way(WLineWay::ImpImp);
        for p in points {
            wl.add(p);
        }
        for v in &aline.vertices {
            wl.vertices.push(*v);
        }
        wl.ensure_end_vertices();
        vec![wl]
    }
}
