//! A point / segment record of the polygon hidden line removal. Source:
//! `HLRAlgo_BiPoint`
//! (`ModelingAlgorithms/TKHLR/HLRAlgo/HLRAlgo_BiPoint.hxx:29-260`,
//! `HLRAlgo_BiPoint.cxx:24-240`).

use occt_core::gp::{GpXY, GpXyz};

/// `EMskFlags` (`hxx:248-255`).
const EMSK_RG1_LINE: i32 = 1;
const EMSK_RGN_LINE: i32 = 2;
const EMSK_OUT_LINE: i32 = 4;
const EMSK_INT_LINE: i32 = 8;
const EMSK_HIDDEN: i32 = 16;

/// `HLRAlgo_BiPoint::IndicesT` (`hxx:31-58`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BiPointIndices {
    pub shape_index: i32,
    pub face_conex1: i32,
    pub face1_pt1: i32,
    pub face1_pt2: i32,
    pub face_conex2: i32,
    pub face2_pt1: i32,
    pub face2_pt2: i32,
    pub min_seg: i32,
    pub max_seg: i32,
    pub seg_flags: i32,
}

impl Default for BiPointIndices {
    /// `IndicesT()` (`hxx:33-46`): `ShapeIndex = -1`, the rest `0`.
    fn default() -> Self {
        Self {
            shape_index: -1,
            face_conex1: 0,
            face1_pt1: 0,
            face1_pt2: 0,
            face_conex2: 0,
            face2_pt1: 0,
            face2_pt2: 0,
            min_seg: 0,
            max_seg: 0,
            seg_flags: 0,
        }
    }
}

/// `HLRAlgo_BiPoint::PointsT` (`hxx:60-72`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BiPointPoints {
    pub pnt1: GpXyz,
    pub pnt2: GpXyz,
    pub pntp1: GpXyz,
    pub pntp2: GpXyz,
}

impl Default for BiPointPoints {
    fn default() -> Self {
        Self {
            pnt1: GpXyz::new(0.0, 0.0, 0.0),
            pnt2: GpXyz::new(0.0, 0.0, 0.0),
            pntp1: GpXyz::new(0.0, 0.0, 0.0),
            pntp2: GpXyz::new(0.0, 0.0, 0.0),
        }
    }
}

impl BiPointPoints {
    /// `PntP12D()` (`hxx:66`).
    pub fn pntp12d(&self) -> GpXY {
        GpXY::new(self.pntp1.x, self.pntp1.y)
    }

    /// `PntP22D()` (`hxx:68`).
    pub fn pntp22d(&self) -> GpXY {
        GpXY::new(self.pntp2.x, self.pntp2.y)
    }
}

/// A segment / point of the polygon HLR, carrying the 3D and projected
/// endpoints plus the face connexity indices. Source: `HLRAlgo_BiPoint`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct HlrAlgoBiPoint {
    indices: BiPointIndices,
    points: BiPointPoints,
}

impl HlrAlgoBiPoint {
    /// `HLRAlgo_BiPoint() = default` (`hxx:75`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `HLRAlgo_BiPoint(X1..ZT2, Index, reg1, regn, outl, intl)`
    /// (`cxx:24-56`).
    #[allow(clippy::too_many_arguments)]
    pub fn from_flags(
        x1: f64,
        y1: f64,
        z1: f64,
        x2: f64,
        y2: f64,
        z2: f64,
        xt1: f64,
        yt1: f64,
        zt1: f64,
        xt2: f64,
        yt2: f64,
        zt2: f64,
        index: i32,
        reg1: bool,
        regn: bool,
        outl: bool,
        intl: bool,
    ) -> Self {
        let mut s = Self::default();
        s.points.pnt1 = GpXyz::new(x1, y1, z1);
        s.points.pnt2 = GpXyz::new(x2, y2, z2);
        s.points.pntp1 = GpXyz::new(xt1, yt1, zt1);
        s.points.pntp2 = GpXyz::new(xt2, yt2, zt2);
        s.indices.shape_index = index;
        s.indices.seg_flags = 0;
        s.set_rg1_line(reg1);
        s.set_rgn_line(regn);
        s.set_out_line(outl);
        s.set_int_line(intl);
        s.set_hidden(false);
        s
    }

    /// `HLRAlgo_BiPoint(X1..ZT2, Index, flag)` (`cxx:60-86`).
    #[allow(clippy::too_many_arguments)]
    pub fn with_flag(
        x1: f64,
        y1: f64,
        z1: f64,
        x2: f64,
        y2: f64,
        z2: f64,
        xt1: f64,
        yt1: f64,
        zt1: f64,
        xt2: f64,
        yt2: f64,
        zt2: f64,
        index: i32,
        flag: i32,
    ) -> Self {
        let mut s = Self::default();
        s.points.pnt1 = GpXyz::new(x1, y1, z1);
        s.points.pnt2 = GpXyz::new(x2, y2, z2);
        s.points.pntp1 = GpXyz::new(xt1, yt1, zt1);
        s.points.pntp2 = GpXyz::new(xt2, yt2, zt2);
        s.indices.shape_index = index;
        s.indices.seg_flags = flag;
        s.set_hidden(false);
        s
    }

    /// `HLRAlgo_BiPoint(X1..ZT2, Index, i1, i1p1, i1p2, reg1, regn, outl,
    /// intl)` (`cxx:88-125`).
    #[allow(clippy::too_many_arguments)]
    pub fn from_face1(
        x1: f64,
        y1: f64,
        z1: f64,
        x2: f64,
        y2: f64,
        z2: f64,
        xt1: f64,
        yt1: f64,
        zt1: f64,
        xt2: f64,
        yt2: f64,
        zt2: f64,
        index: i32,
        i1: i32,
        i1p1: i32,
        i1p2: i32,
        reg1: bool,
        regn: bool,
        outl: bool,
        intl: bool,
    ) -> Self {
        let mut s = Self::default();
        s.points.pnt1 = GpXyz::new(x1, y1, z1);
        s.points.pnt2 = GpXyz::new(x2, y2, z2);
        s.points.pntp1 = GpXyz::new(xt1, yt1, zt1);
        s.points.pntp2 = GpXyz::new(xt2, yt2, zt2);
        s.indices.shape_index = index;
        s.indices.face_conex1 = i1;
        s.indices.face1_pt1 = i1p1;
        s.indices.face1_pt2 = i1p2;
        s.indices.face_conex2 = 0;
        s.indices.face2_pt1 = 0;
        s.indices.face2_pt2 = 0;
        s.indices.seg_flags = 0;
        s.set_rg1_line(reg1);
        s.set_rgn_line(regn);
        s.set_out_line(outl);
        s.set_int_line(intl);
        s.set_hidden(false);
        s
    }

    /// `HLRAlgo_BiPoint(X1..ZT2, Index, i1, i1p1, i1p2, flag)`
    /// (`cxx:127-158`).
    #[allow(clippy::too_many_arguments)]
    pub fn from_face1_with_flag(
        x1: f64,
        y1: f64,
        z1: f64,
        x2: f64,
        y2: f64,
        z2: f64,
        xt1: f64,
        yt1: f64,
        zt1: f64,
        xt2: f64,
        yt2: f64,
        zt2: f64,
        index: i32,
        i1: i32,
        i1p1: i32,
        i1p2: i32,
        flag: i32,
    ) -> Self {
        let mut s = Self::default();
        s.points.pnt1 = GpXyz::new(x1, y1, z1);
        s.points.pnt2 = GpXyz::new(x2, y2, z2);
        s.points.pntp1 = GpXyz::new(xt1, yt1, zt1);
        s.points.pntp2 = GpXyz::new(xt2, yt2, zt2);
        s.indices.shape_index = index;
        s.indices.face_conex1 = i1;
        s.indices.face1_pt1 = i1p1;
        s.indices.face1_pt2 = i1p2;
        s.indices.face_conex2 = 0;
        s.indices.face2_pt1 = 0;
        s.indices.face2_pt2 = 0;
        s.indices.seg_flags = flag;
        s.set_hidden(false);
        s
    }

    /// `HLRAlgo_BiPoint(X1..ZT2, Index, i1, i1p1, i1p2, i2, i2p1, i2p2,
    /// reg1, regn, outl, intl)` (`cxx:160-202`).
    #[allow(clippy::too_many_arguments)]
    pub fn from_two_faces(
        x1: f64,
        y1: f64,
        z1: f64,
        x2: f64,
        y2: f64,
        z2: f64,
        xt1: f64,
        yt1: f64,
        zt1: f64,
        xt2: f64,
        yt2: f64,
        zt2: f64,
        index: i32,
        i1: i32,
        i1p1: i32,
        i1p2: i32,
        i2: i32,
        i2p1: i32,
        i2p2: i32,
        reg1: bool,
        regn: bool,
        outl: bool,
        intl: bool,
    ) -> Self {
        let mut s = Self::default();
        s.points.pnt1 = GpXyz::new(x1, y1, z1);
        s.points.pnt2 = GpXyz::new(x2, y2, z2);
        s.points.pntp1 = GpXyz::new(xt1, yt1, zt1);
        s.points.pntp2 = GpXyz::new(xt2, yt2, zt2);
        s.indices.shape_index = index;
        s.indices.face_conex1 = i1;
        s.indices.face1_pt1 = i1p1;
        s.indices.face1_pt2 = i1p2;
        s.indices.face_conex2 = i2;
        s.indices.face2_pt1 = i2p1;
        s.indices.face2_pt2 = i2p2;
        s.indices.seg_flags = 0;
        s.set_rg1_line(reg1);
        s.set_rgn_line(regn);
        s.set_out_line(outl);
        s.set_int_line(intl);
        s.set_hidden(false);
        s
    }

    /// `HLRAlgo_BiPoint(X1..ZT2, Index, i1, i1p1, i1p2, i2, i2p1, i2p2,
    /// flag)` (`cxx:204-240`).
    #[allow(clippy::too_many_arguments)]
    pub fn from_two_faces_with_flag(
        x1: f64,
        y1: f64,
        z1: f64,
        x2: f64,
        y2: f64,
        z2: f64,
        xt1: f64,
        yt1: f64,
        zt1: f64,
        xt2: f64,
        yt2: f64,
        zt2: f64,
        index: i32,
        i1: i32,
        i1p1: i32,
        i1p2: i32,
        i2: i32,
        i2p1: i32,
        i2p2: i32,
        flag: i32,
    ) -> Self {
        let mut s = Self::default();
        s.points.pnt1 = GpXyz::new(x1, y1, z1);
        s.points.pnt2 = GpXyz::new(x2, y2, z2);
        s.points.pntp1 = GpXyz::new(xt1, yt1, zt1);
        s.points.pntp2 = GpXyz::new(xt2, yt2, zt2);
        s.indices.shape_index = index;
        s.indices.face_conex1 = i1;
        s.indices.face1_pt1 = i1p1;
        s.indices.face1_pt2 = i1p2;
        s.indices.face_conex2 = i2;
        s.indices.face2_pt1 = i2p1;
        s.indices.face2_pt2 = i2p2;
        s.indices.seg_flags = flag;
        s.set_hidden(false);
        s
    }

    /// `Rg1Line()` (`hxx:192`).
    pub fn rg1_line(&self) -> bool {
        self.indices.seg_flags & EMSK_RG1_LINE != 0
    }

    /// `Rg1Line(B)` (`hxx:194-202`).
    pub fn set_rg1_line(&mut self, b: bool) {
        set_flag(&mut self.indices.seg_flags, EMSK_RG1_LINE, b);
    }

    /// `RgNLine()` (`hxx:204`).
    pub fn rgn_line(&self) -> bool {
        self.indices.seg_flags & EMSK_RGN_LINE != 0
    }

    /// `RgNLine(B)` (`hxx:206-214`).
    pub fn set_rgn_line(&mut self, b: bool) {
        set_flag(&mut self.indices.seg_flags, EMSK_RGN_LINE, b);
    }

    /// `OutLine()` (`hxx:216`).
    pub fn out_line(&self) -> bool {
        self.indices.seg_flags & EMSK_OUT_LINE != 0
    }

    /// `OutLine(B)` (`hxx:218-226`).
    pub fn set_out_line(&mut self, b: bool) {
        set_flag(&mut self.indices.seg_flags, EMSK_OUT_LINE, b);
    }

    /// `IntLine()` (`hxx:228`).
    pub fn int_line(&self) -> bool {
        self.indices.seg_flags & EMSK_INT_LINE != 0
    }

    /// `IntLine(B)` (`hxx:230-238`).
    pub fn set_int_line(&mut self, b: bool) {
        set_flag(&mut self.indices.seg_flags, EMSK_INT_LINE, b);
    }

    /// `Hidden()` (`hxx:240`).
    pub fn hidden(&self) -> bool {
        self.indices.seg_flags & EMSK_HIDDEN != 0
    }

    /// `Hidden(B)` (`hxx:242-250`).
    pub fn set_hidden(&mut self, b: bool) {
        set_flag(&mut self.indices.seg_flags, EMSK_HIDDEN, b);
    }

    /// `Indices()` (`hxx:252`).
    pub fn indices(&self) -> &BiPointIndices {
        &self.indices
    }

    /// `Indices()` mutable form (the OCCT non-const accessor).
    pub fn change_indices(&mut self) -> &mut BiPointIndices {
        &mut self.indices
    }

    /// `Points()` (`hxx:254`).
    pub fn points(&self) -> &BiPointPoints {
        &self.points
    }

    /// `Points()` mutable form (the OCCT non-const accessor).
    pub fn change_points(&mut self) -> &mut BiPointPoints {
        &mut self.points
    }
}

fn set_flag(flag: &mut i32, mask: i32, b: bool) {
    if b {
        *flag |= mask;
    } else {
        *flag &= !mask;
    }
}
