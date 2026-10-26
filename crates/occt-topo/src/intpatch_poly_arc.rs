//! `IntPatch_PolyArc` - the 2D polygon of a restriction arc, used as the
//! "other object" of a restriction / line interference.
//!
//! Source: `ModelingAlgorithms/TKGeomAlgo/IntPatch/IntPatch_PolyArc.hxx` and
//! `.cxx`.
//!
//! Ported members (OCCT 8.0.0):
//! * ctor `(A, NbSample, Pfirst, Plast, BoxOtherPolygon)` (`cxx:36-186`),
//!   including the `Standard_ConstructionError` guard (`cxx:50-53`), the
//!   "significant zone" loop that widens `r` when the arc falls entirely
//!   outside the other polygon's box (`cxx:82-177`) and the `1.2` error
//!   inflation with its `1e-8` floor (`cxx:178-182`).
//! * `Closed` (`cxx:188-191`), `NbPoints` (`cxx:193-196`), `Point`
//!   (`cxx:198-207`), `Parameter` (`cxx:209-212`), `SetOffset` (`cxx:214-228`).
//!
//! Translation notes:
//! * `NCollection_Array1<gp_Pnt2d> brise` / `NCollection_Array1<double> param`
//!   (`hxx:56-57`) are `1 .. NbSample`. The port keeps them in `Vec`s with a
//!   dummy element 0 so the 1-based OCCT indices index the port directly.
//! * `MinMax` (`cxx:22-33`) is a file-static helper; it is inlined below.

use occt_core::bnd::box2d::BndBox2d;
use occt_core::gp::GpPnt2d;
use occt_core::intf::IntfPolygon2d;
use occt_geom2d::Curve2d;

use crate::intpatch_polygo::IntPatchPolygo;

/// `RealFirst()` / `RealLast()` (`Standard_Real.hxx`).
const REAL_FIRST: f64 = -f64::MAX;
const REAL_LAST: f64 = f64::MAX;

/// `MinMax(a1, a2, amin, amax)` (`cxx:22-33`).
fn min_max(a1: f64, a2: f64) -> (f64, f64) {
    if a1 < a2 {
        (a1, a2)
    } else {
        (a2, a1)
    }
}

/// `IntPatch_PolyArc` (`hxx:32-62`).
pub struct IntPatchPolyArc {
    /// `myBox` (`Intf_Polygon2d::myBox`).
    my_box: BndBox2d,
    /// `myError` (`hxx:55`).
    my_error: f64,
    /// `brise` (`hxx:56`), 1-based.
    brise: Vec<GpPnt2d>,
    /// `param` (`hxx:57`), 1-based.
    param: Vec<f64>,
    /// `offsetx` (`hxx:58`).
    offset_x: f64,
    /// `offsety` (`hxx:59`).
    offset_y: f64,
    /// `ferme` (`hxx:60`).
    ferme: bool,
}

impl IntPatchPolyArc {
    /// `IntPatch_PolyArc(A, NbSample, Pfirst, Plast, BoxOtherPolygon)`
    /// (`cxx:36-186`).
    ///
    /// Returns `Err` where OCCT throws `Standard_ConstructionError`
    /// (`cxx:50-53`).
    pub fn new(
        line: &dyn Curve2d,
        nb_sample: usize,
        a_pdeb: f64,
        a_pfin: f64,
        box_other_polygon: &BndBox2d,
    ) -> Result<Self, String> {
        let nb_sample = nb_sample.max(1);
        let mut this = Self {
            my_box: BndBox2d::new(),
            my_error: 0.0,
            brise: vec![GpPnt2d::zero(); nb_sample + 1],
            param: vec![0.0; nb_sample + 1],
            offset_x: 0.0,
            offset_y: 0.0,
            ferme: false,
        };

        let mut pdeb = a_pdeb;
        let mut pfin = a_pfin;

        if pdeb == REAL_FIRST || pfin == REAL_LAST || nb_sample <= 1 {
            return Err("IntPatch_PolyArc: Standard_ConstructionError".into());
        }

        // `BoxOtherPolygon.Get(bxmin, bymin, bxmax, bymax);` (`cxx:71`)
        let (bxmin, bymin, bxmax, bymax) = box_other_polygon
            .get()
            .ok_or_else(|| "IntPatch_PolyArc: void BoxOtherPolygon".to_string())?;
        let mut r = (bxmax - bxmin) + (bymax - bymin);
        let bx0 = (bxmax + bxmin) * 0.5;
        let by0 = (bymax + bymin) * 0.5;

        let mut r2;
        r *= 0.8;
        r2 = r * r * 49.0;
        let mut nbloop = 0;

        let mut index_inf: i64 = nb_sample as i64 + 1;
        let mut index_sup: i64 = 0;

        loop {
            nbloop += 1;
            let pas = (pfin - pdeb) / (nb_sample as f64 - 1.0);
            this.param[1] = pdeb;
            let p2d = line.d0(pdeb);
            let mut xs = p2d.x();
            let mut ys = p2d.y();
            this.brise[1] = GpPnt2d::new(xs, ys);

            this.my_box.set_void();
            this.my_box.add_point(&this.brise[1]);
            this.my_error = 0.0;

            for i in 2..=nb_sample {
                this.param[i] = pdeb + (i as f64 - 1.0) * pas;
                let p2d = line.d0(this.param[i]);
                let x = p2d.x();
                let y = p2d.y();
                this.brise[i] = GpPnt2d::new(x, y);
                let xxs = 0.5 * (xs + x);
                let yys = 0.5 * (ys + y);

                // Significant-zone tracking (`cxx:110-141`).
                let is_mid_pt_in_box = (bx0 - xxs).abs() + (by0 - yys).abs() < r;
                let mut is_seg_out = true;
                if !is_mid_pt_in_box {
                    let d = (x - xs) * (x - xs) + (y - ys) * (y - ys);
                    if d > r2 {
                        let (xmin, xmax) = min_max(xs, x);
                        let (ymin, ymax) = min_max(ys, y);
                        is_seg_out =
                            xmax < bxmin || xmin > bxmax || ymax < bymin || ymin > bymax;
                    }
                }
                if is_mid_pt_in_box || !is_seg_out {
                    if index_inf > i as i64 {
                        index_inf = (i as i64 - 2).max(1);
                    }
                    if index_sup < i as i64 {
                        index_sup = (i as i64 + 1).min(nb_sample as i64);
                    }
                }

                this.my_box.add_point(&this.brise[i]);
                let pm = line.d0(this.param[i] - pas * 0.5);
                let xm = pm.x() - xxs;
                let ym = pm.y() - yys;
                let xm = (xm * xm + ym * ym).sqrt();
                this.my_error = this.my_error.max(xm);
                xs = x;
                ys = y;
            }

            if index_inf > index_sup {
                r += r;
                r2 = r * r * 49.0;
                // OCCT loops again; `nbloop <= 10` bounds it (`cxx:177`).
                if nbloop > 10 {
                    break;
                }
            } else {
                // If the number of significant points is insufficient,
                // resample once more (`cxx:158-176`).
                if (index_sup - index_inf) < (nb_sample as i64 / 2) {
                    nbloop = 10;
                    pdeb = this.param[index_inf as usize];
                    pfin = this.param[index_sup as usize];
                    index_inf = nb_sample as i64 + 1;
                    index_sup = 0;
                } else {
                    break;
                }
            }
        }

        this.my_error *= 1.2;
        if this.my_error < 0.00000001 {
            this.my_error = 0.00000001;
        }
        this.my_box.enlarge(this.my_error);

        // `ferme = (Line->Value(aPdeb).Distance(Line->Value(aPfin)) <= 1e-7);`
        this.ferme = line.value(a_pdeb).distance(&line.value(a_pfin)) <= 1e-7;
        Ok(this)
    }

    /// `Closed()` (`cxx:188-191`).
    pub fn is_closed(&self) -> bool {
        self.ferme
    }

    /// `NbPoints()` (`cxx:193-196`): `brise.Length()`.
    pub fn nb_points(&self) -> usize {
        self.brise.len() - 1
    }

    /// `Point(Index)` (`cxx:198-207`).
    pub fn point(&self, index: usize) -> GpPnt2d {
        if self.offset_x == 0.0 && self.offset_y == 0.0 {
            return self.brise[index];
        }
        let p = self.brise[index];
        GpPnt2d::new(p.x() + self.offset_x, p.y() + self.offset_y)
    }

    /// `Parameter(Index)` (`cxx:209-212`).
    pub fn parameter(&self, index: usize) -> f64 {
        self.param[index]
    }

    /// `SetOffset(ox, oy)` (`cxx:214-228`).
    pub fn set_offset(&mut self, ox: f64, oy: f64) {
        let (xmin, ymin, xmax, ymax) = match self.my_box.get() {
            Some(v) => v,
            None => return,
        };
        let g = self.my_box.gap();

        self.my_box.set_void();
        self.my_box
            .update(xmin - self.offset_x, ymin - self.offset_y, xmax - self.offset_x, ymax - self.offset_y);
        self.offset_x = ox;
        self.offset_y = oy;
        self.my_box
            .update(xmin + self.offset_x, ymin + self.offset_y, xmax + self.offset_x, ymax + self.offset_y);
        self.my_box.set_gap(g);
    }
}

impl IntfPolygon2d for IntPatchPolyArc {
    /// `Intf_Polygon2d::Bounding()`.
    fn bounding(&self) -> &BndBox2d {
        &self.my_box
    }

    /// `Intf_Polygon2d::Closed()` (`cxx:188-191`).
    fn is_closed(&self) -> bool {
        self.ferme
    }

    /// `DeflectionOverEstimation()` (`lxx:24-27`): `myError`.
    fn deflection_over_estimation(&self) -> f64 {
        self.my_error
    }

    /// `NbSegments()` (`lxx:31-34`).
    fn nb_segments(&self) -> usize {
        self.nb_points().saturating_sub(1)
    }

    /// `Segment(Index, theBegin, theEnd)` (`lxx:38-43`).
    fn segment(&self, index: usize, begin: &mut GpPnt2d, end: &mut GpPnt2d) {
        *begin = self.point(index);
        *end = self.point(index + 1);
    }
}

impl IntPatchPolygo for IntPatchPolyArc {
    fn error(&self) -> f64 {
        self.my_error
    }

    fn nb_points(&self) -> usize {
        IntPatchPolyArc::nb_points(self)
    }

    fn point(&self, index: usize) -> GpPnt2d {
        IntPatchPolyArc::point(self, index)
    }
}
