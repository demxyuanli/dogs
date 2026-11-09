//! Hidden line removal on triangulations. Source: `HLRAlgo_PolyAlgo`
//! (`ModelingAlgorithms/TKHLR/HLRAlgo/HLRAlgo_PolyAlgo.hxx:27-94`,
//! `HLRAlgo_PolyAlgo.cxx:29-597`).
//!
//! The OCCT `mySegListIt` is an `NCollection_List<HLRAlgo_BiPoint>::Iterator`;
//! the port keeps the same cursor as `(seg_shell, seg_pos)`, both fields of the
//! algorithm, so that `NextHide` / `NextShow` reproduce the "advance the
//! current list, then walk the following shells" control flow.

use occt_core::bnd::BndBox;
use occt_core::precision::INFINITE;

use super::bi_point::{BiPointPoints, HlrAlgoBiPoint};
use super::edge_status::HlrAlgoEdgeStatus;
use super::poly_data::{HlrAlgoPolyData, PolyTriangle, POLY_BBOX_MASK_A, POLY_BBOX_MASK_B};
use super::poly_records::poly_mask;
use super::poly_shell_data::HlrAlgoPolyShellData;

/// What `HLRAlgo_PolyAlgo::Hide` / `Show` report about the current segment.
/// OCCT returns the `PointsT&` and fills the five out-parameters
/// (`cxx:502-521`, `:583-596`); the port returns them together. `PointsT` is a
/// value type, so copying it is equivalent to handing out the reference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PolyBiPointInfo {
    pub points: BiPointPoints,
    pub index: i32,
    pub reg1: bool,
    pub regn: bool,
    pub outl: bool,
    pub intl: bool,
}

/// Removes hidden lines on triangulations. Source: `HLRAlgo_PolyAlgo`.
#[derive(Debug, Clone, PartialEq)]
pub struct HlrAlgoPolyAlgo {
    poly_shell: Vec<HlrAlgoPolyShellData>,
    triangle: PolyTriangle,
    nbr_shell: usize,
    /// `myCurShell`: 1-based, `0` before the first `NextHide` / `NextShow`.
    cur_shell: usize,
    found: bool,
    /// Shell that owns the segment list the cursor walks.
    seg_shell: usize,
    /// Index of the current segment inside `seg_shell`'s edge list.
    seg_pos: usize,
}

impl Default for HlrAlgoPolyAlgo {
    /// `HLRAlgo_PolyAlgo()` (`cxx:30-37`): only the two tolerances are set.
    fn default() -> Self {
        Self {
            poly_shell: Vec::new(),
            triangle: PolyTriangle {
                tol_param: 0.000_000_01,
                tol_ang: 0.0001,
                ..PolyTriangle::default()
            },
            nbr_shell: 0,
            cur_shell: 0,
            found: false,
            seg_shell: 0,
            seg_pos: 0,
        }
    }
}

impl HlrAlgoPolyAlgo {
    /// `HLRAlgo_PolyAlgo()` (`cxx:30-37`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Init(theNbShells)` (`cxx:42-45`).
    pub fn init(&mut self, nb_shells: usize) {
        self.poly_shell
            .resize(nb_shells, HlrAlgoPolyShellData::default());
    }

    /// `PolyShell()` (`hxx:35-38`).
    pub fn poly_shell(&self) -> &[HlrAlgoPolyShellData] {
        &self.poly_shell
    }

    /// `ChangePolyShell()` (`hxx:40`).
    pub fn change_poly_shell(&mut self) -> &mut Vec<HlrAlgoPolyShellData> {
        &mut self.poly_shell
    }

    /// `Clear()` (`cxx:47-52`).
    pub fn clear(&mut self) {
        self.poly_shell.clear();
        self.nbr_shell = 0;
    }

    /// `Update()` (`cxx:55-463`): build the integer bounding boxes of every
    /// segment, face and shell, and fill the hiding-triangle planes.
    pub fn update(&mut self) {
        let big = INFINITE;
        let mut a_box = BndBox::new();
        self.nbr_shell = self.poly_shell.len();
        for psd in &self.poly_shell {
            psd.update_global_min_max(&mut a_box);
        }

        // `cxx:80-100`.
        let (a_x_min, a_x_max, a_y_min, a_y_max, a_z_min, a_z_max) =
            a_box.get().unwrap_or((-big, big, -big, big, -big, big));
        let dx = a_x_max - a_x_min;
        let dy = a_y_max - a_y_min;
        let dz = a_z_max - a_z_min;
        let mut precad = dx;
        if precad < dy {
            precad = dy;
        }
        if precad < dz {
            precad = dz;
        }
        self.triangle.tolerance = precad * self.triangle.tol_param;
        precad *= 0.01;
        let sur_dx = 1020.0 / (dx + precad);
        let sur_dy = 1020.0 / (dy + precad);
        let sur_dz = 508.0 / (dz + precad);
        precad *= 0.5;
        let deca_x = -a_x_min + precad;
        let deca_y = -a_y_min + precad;
        let deca_z = -a_z_min + precad;

        for s in 0..self.poly_shell.len() {
            let mut x_shell_min = big;
            let mut y_shell_min = big;
            let mut z_shell_min = big;
            let mut x_shell_max = -big;
            let mut y_shell_max = -big;
            let mut z_shell_max = -big;

            // `cxx:112-186`: segment boxes.
            {
                let psd = &mut self.poly_shell[s];
                for bp in psd.change_edges().iter_mut() {
                    let p = bp.points();
                    let (x_seg_min, x_seg_max) = if p.pntp1.x() < p.pntp2.x() {
                        (p.pntp1.x(), p.pntp2.x())
                    } else {
                        (p.pntp2.x(), p.pntp1.x())
                    };
                    let (y_seg_min, y_seg_max) = if p.pntp1.y() < p.pntp2.y() {
                        (p.pntp1.y(), p.pntp2.y())
                    } else {
                        (p.pntp2.y(), p.pntp1.y())
                    };
                    let (z_seg_min, z_seg_max) = if p.pntp1.z() < p.pntp2.z() {
                        (p.pntp1.z(), p.pntp2.z())
                    } else {
                        (p.pntp2.z(), p.pntp1.z())
                    };
                    let nx_min = ((deca_x + x_seg_min) * sur_dx) as i32;
                    let ny_min = ((deca_y + y_seg_min) * sur_dy) as i32;
                    let nz_min = ((deca_z + z_seg_min) * sur_dz) as i32;
                    let nx_max = ((deca_x + x_seg_max) * sur_dx) as i32;
                    let ny_max = ((deca_y + y_seg_max) * sur_dy) as i32;
                    let nz_max = ((deca_z + z_seg_max) * sur_dz) as i32;
                    let ind = bp.change_indices();
                    ind.min_seg = (ny_min + (nx_min << 11)) << 10;
                    ind.min_seg += nz_min;
                    ind.max_seg = (ny_max + (nx_max << 11)) << 10;
                    ind.max_seg += nz_max + 0x0000_0200;
                    if x_shell_min > x_seg_min {
                        x_shell_min = x_seg_min;
                    }
                    if x_shell_max < x_seg_max {
                        x_shell_max = x_seg_max;
                    }
                    if y_shell_min > y_seg_min {
                        y_shell_min = y_seg_min;
                    }
                    if y_shell_max < y_seg_max {
                        y_shell_max = y_seg_max;
                    }
                    if z_shell_min > z_seg_min {
                        z_shell_min = z_seg_min;
                    }
                    if z_shell_max < z_seg_max {
                        z_shell_max = z_seg_max;
                    }
                }
            }

            // `cxx:187-429`: hiding faces.
            let nb_face = self.poly_shell[s].poly_data().len();
            let mut nb_fa_hi = 0usize;
            for j in 1..=nb_face {
                if !self.poly_shell[s].poly_data()[j - 1].hiding() {
                    continue;
                }
                nb_fa_hi += 1;
                let mut x_polyt_min = big;
                let mut y_polyt_min = big;
                let mut z_polyt_min = big;
                let mut x_polyt_max = -big;
                let mut y_polyt_max = -big;
                let mut z_polyt_max = -big;
                let mut nb_hide = 0usize;
                // `cxx:205`: the plane coefficients live across the triangle
                // loop of one face, so a degenerate triangle (`dn == 0`) keeps
                // the previous triangle's normal.
                let mut a = 0.0f64;
                let mut b = 0.0f64;
                let mut c = 0.0f64;
                // `d` is recomputed from `a, b, c` for every triangle, so it
                // needs no initial value (`cxx:205` zeroes it anyway).
                let mut d;

                let pd: &mut HlrAlgoPolyData = &mut self.poly_shell[s].change_poly_data()[j - 1];
                let nb_t = pd.t_data().len();
                for otheri in 1..=nb_t {
                    let td = pd.t_data()[otheri - 1];
                    if td.flags & poly_mask::FMSK_HIDING == 0 {
                        continue;
                    }
                    let p1 = pd.nodes()[td.node1 as usize - 1];
                    let p2 = pd.nodes()[td.node2 as usize - 1];
                    let p3 = pd.nodes()[td.node3 as usize - 1];
                    let (x1, y1, z1) = (p1.x(), p1.y(), p1.z());
                    let (x2, y2, z2) = (p2.x(), p2.y(), p2.z());
                    let (x3, y3, z3) = (p3.x(), p3.y(), p3.z());
                    let mut x_trian_min = x1;
                    let mut x_trian_max = x1;
                    let mut y_trian_min = y1;
                    let mut y_trian_max = y1;
                    let mut z_trian_min = z1;
                    let mut z_trian_max = z1;
                    // `cxx:231-277`: the else-if branches, not two ifs.
                    if x_trian_min > x2 {
                        x_trian_min = x2;
                    } else if x_trian_max < x2 {
                        x_trian_max = x2;
                    }
                    if y_trian_min > y2 {
                        y_trian_min = y2;
                    } else if y_trian_max < y2 {
                        y_trian_max = y2;
                    }
                    if z_trian_min > z2 {
                        z_trian_min = z2;
                    } else if z_trian_max < z2 {
                        z_trian_max = z2;
                    }
                    if x_trian_min > x3 {
                        x_trian_min = x3;
                    } else if x_trian_max < x3 {
                        x_trian_max = x3;
                    }
                    if y_trian_min > y3 {
                        y_trian_min = y3;
                    } else if y_trian_max < y3 {
                        y_trian_max = y3;
                    }
                    if z_trian_min > z3 {
                        z_trian_min = z3;
                    } else if z_trian_max < z3 {
                        z_trian_max = z3;
                    }
                    let nx_min = ((deca_x + x_trian_min) * sur_dx) as i32;
                    let ny_min = ((deca_y + y_trian_min) * sur_dy) as i32;
                    let nz_min = ((deca_z + z_trian_min) * sur_dz) as i32;
                    let nx_max = ((deca_x + x_trian_max) * sur_dx) as i32;
                    let ny_max = ((deca_y + y_trian_max) * sur_dy) as i32;
                    let nz_max = ((deca_z + z_trian_max) * sur_dz) as i32;
                    let min_trian = ((ny_min + (nx_min << 11)) << 10) + nz_min - 0x0000_0200;
                    let max_trian = ((ny_max + (nx_max << 11)) << 10) + nz_max;
                    let dx1 = x2 - x1;
                    let dy1 = y2 - y1;
                    let dz1 = z2 - z1;
                    let dx2 = x3 - x2;
                    let dy2 = y3 - y2;
                    let dz2 = z3 - z2;
                    let dx3 = x1 - x3;
                    let dy3 = y1 - y3;
                    let dnx = dy1 * dz2 - dy2 * dz1;
                    let dny = dz1 * dx2 - dz2 * dx1;
                    let dnz = dx1 * dy2 - dx2 * dy1;
                    let dn = (dnx * dnx + dny * dny + dnz * dnz).sqrt();
                    // `cxx:306-310`: only normalised when `dn > 0`.
                    if dn > 0.0 {
                        a = dnx / dn;
                        b = dny / dn;
                        c = dnz / dn;
                    }
                    d = a * x1 + b * y1 + c * z1;
                    nb_hide += 1;
                    if let Some(phdat) = pd.ph_dat_mut() {
                        phdat[nb_hide - 1].set(otheri as i32, min_trian, max_trian, a, b, c, d);
                    }
                    let adx1 = if dx1 < 0.0 { -dx1 } else { dx1 };
                    let ady1 = if dy1 < 0.0 { -dy1 } else { dy1 };
                    let adx2 = if dx2 < 0.0 { -dx2 } else { dx2 };
                    let ady2 = if dy2 < 0.0 { -dy2 } else { dy2 };
                    let adx3 = if dx3 < 0.0 { -dx3 } else { dx3 };
                    let ady3 = if dy3 < 0.0 { -dy3 } else { dy3 };
                    let mut flags = td.flags;
                    set_mask(&mut flags, poly_mask::EMSK_GR_A_LIN1, adx1 > ady1);
                    set_mask(&mut flags, poly_mask::EMSK_GR_A_LIN2, adx2 > ady2);
                    set_mask(&mut flags, poly_mask::EMSK_GR_A_LIN3, adx3 > ady3);
                    pd.t_data_mut()[otheri - 1].flags = flags;
                    if x_polyt_min > x_trian_min {
                        x_polyt_min = x_trian_min;
                    }
                    if x_polyt_max < x_trian_max {
                        x_polyt_max = x_trian_max;
                    }
                    if y_polyt_min > y_trian_min {
                        y_polyt_min = y_trian_min;
                    }
                    if y_polyt_max < y_trian_max {
                        y_polyt_max = y_trian_max;
                    }
                    if z_polyt_min > z_trian_min {
                        z_polyt_min = z_trian_min;
                    }
                    if z_polyt_max < z_trian_max {
                        z_polyt_max = z_trian_max;
                    }
                }

                // `cxx:394-405`.
                let nx_min = ((deca_x + x_polyt_min) * sur_dx) as i32;
                let ny_min = ((deca_y + y_polyt_min) * sur_dy) as i32;
                let nz_min = ((deca_z + z_polyt_min) * sur_dz) as i32;
                let nx_max = ((deca_x + x_polyt_max) * sur_dx) as i32;
                let ny_max = ((deca_y + y_polyt_max) * sur_dy) as i32;
                let nz_max = ((deca_z + z_polyt_max) * sur_dz) as i32;
                let fi = pd.change_indices();
                fi.min = ((ny_min + (nx_min << 11)) << 10) + nz_min - 0x0000_0200;
                fi.max = ((ny_max + (nx_max << 11)) << 10) + nz_max;
                if x_shell_min > x_polyt_min {
                    x_shell_min = x_polyt_min;
                }
                if x_shell_max < x_polyt_max {
                    x_shell_max = x_polyt_max;
                }
                if y_shell_min > y_polyt_min {
                    y_shell_min = y_polyt_min;
                }
                if y_shell_max < y_polyt_max {
                    y_shell_max = y_polyt_max;
                }
                if z_shell_min > z_polyt_min {
                    z_shell_min = z_polyt_min;
                }
                if z_shell_max < z_polyt_max {
                    z_shell_max = z_polyt_max;
                }
            }

            // `cxx:431-462`.
            let psd = &mut self.poly_shell[s];
            if nb_fa_hi > 0 {
                let nx_min = ((deca_x + x_shell_min) * sur_dx) as i32;
                let ny_min = ((deca_y + y_shell_min) * sur_dy) as i32;
                let nz_min = ((deca_z + z_shell_min) * sur_dz) as i32;
                let nx_max = ((deca_x + x_shell_max) * sur_dx) as i32;
                let ny_max = ((deca_y + y_shell_max) * sur_dy) as i32;
                let nz_max = ((deca_z + z_shell_max) * sur_dz) as i32;
                let si = psd.change_indices();
                si.min = ((ny_min + (nx_min << 11)) << 10) + nz_min - 0x0000_0200;
                si.max = ((ny_max + (nx_max << 11)) << 10) + nz_max;
                psd.update_hiding(nb_fa_hi);
                for j in 1..=nb_face {
                    if psd.poly_data()[j - 1].hiding() {
                        psd.push_hiding_poly_data(j - 1);
                    }
                }
            } else {
                psd.update_hiding(0);
                let si = psd.change_indices();
                si.min = 0;
                si.max = 0;
            }
        }
    }

    /// `InitHide()` (`hxx:48-52`).
    pub fn init_hide(&mut self) {
        self.cur_shell = 0;
        self.next_hide();
    }

    /// `MoreHide()` (`hxx:54`).
    pub fn more_hide(&self) -> bool {
        self.found
    }

    /// `NextHide()` (`cxx:467-496`).
    pub fn next_hide(&mut self) {
        self.found = false;
        if self.cur_shell != 0 {
            self.seg_pos += 1;
            if self.seg_more() {
                self.found = true;
            }
        }
        if !self.found {
            self.cur_shell += 1;
            while self.cur_shell <= self.nbr_shell && !self.found {
                self.seg_shell = self.cur_shell;
                self.seg_pos = 0;
                if self.seg_more() {
                    self.found = true;
                } else {
                    self.cur_shell += 1;
                }
            }
        }
    }

    /// `Hide(status, Index, reg1, regn, outl, intl)` (`cxx:502-547`).
    pub fn hide(&mut self, status: &mut HlrAlgoEdgeStatus) -> PolyBiPointInfo {
        let bp = &self.poly_shell[self.seg_shell - 1].edges()[self.seg_pos];
        let points = *bp.points();
        let indices = *bp.indices();
        *status = HlrAlgoEdgeStatus::from_bounds(
            0.0,
            self.triangle.tol_param as f32,
            1.0,
            self.triangle.tol_param as f32,
        );
        let info = PolyBiPointInfo {
            points,
            index: indices.shape_index,
            reg1: bp.rg1_line(),
            regn: bp.rgn_line(),
            outl: bp.out_line(),
            intl: bp.int_line(),
        };
        if bp.hidden() {
            status.hide_all();
            return info;
        }
        let cur_shell = self.cur_shell;
        let nbr_shell = self.nbr_shell;
        let triangle = &mut self.triangle;
        for s in 1..=nbr_shell {
            let psd = &self.poly_shell[s - 1];
            if !psd.hiding() {
                continue;
            }
            let si = psd.indices();
            if ((si.max - indices.min_seg) & POLY_BBOX_MASK_A) != 0
                || ((indices.max_seg - si.min) & POLY_BBOX_MASK_B) != 0
            {
                continue;
            }
            let is_hiding_shell = s == cur_shell;
            for pd in psd.hiding_poly_data() {
                pd.hide_by_poly_data(&points, triangle, &indices, is_hiding_shell, status);
            }
        }
        info
    }

    /// `InitShow()` (`hxx:67-71`).
    pub fn init_show(&mut self) {
        self.cur_shell = 0;
        self.next_show();
    }

    /// `MoreShow()` (`hxx:73`).
    pub fn more_show(&self) -> bool {
        self.found
    }

    /// `NextShow()` (`cxx:551-578`) — same walk as [`Self::next_hide`].
    pub fn next_show(&mut self) {
        self.next_hide();
    }

    /// `Show(Index, reg1, regn, outl, intl)` (`cxx:583-597`).
    pub fn show(&self) -> PolyBiPointInfo {
        let bp: &HlrAlgoBiPoint = &self.poly_shell[self.seg_shell - 1].edges()[self.seg_pos];
        let indices = *bp.indices();
        PolyBiPointInfo {
            points: *bp.points(),
            index: indices.shape_index,
            reg1: bp.rg1_line(),
            regn: bp.rgn_line(),
            outl: bp.out_line(),
            intl: bp.int_line(),
        }
    }

    /// `mySegListIt.More()`: the cursor of [`Self::seg_shell`].
    fn seg_more(&self) -> bool {
        self.seg_shell >= 1
            && self.seg_shell <= self.nbr_shell
            && self.seg_pos < self.poly_shell[self.seg_shell - 1].edges().len()
    }
}

fn set_mask(flags: &mut i32, mask: i32, on: bool) {
    if on {
        *flags |= mask;
    } else {
        *flags &= !mask;
    }
}
