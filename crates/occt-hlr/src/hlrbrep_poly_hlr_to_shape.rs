//! Filtering of the polygonal HLR results into 2D shapes.
//!
//! Source: `HLRBRep_PolyHLRToShape`
//! (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_PolyHLRToShape.hxx:58-123`,
//! `HLRBRep_PolyHLRToShape.cxx:37-212`, `HLRBRep_PolyHLRToShape.lxx:21-143`).
//!
//! * [`HlrBRepPolyHlrToShape::update`] is `Update` (`cxx:37-79`): every bi-point
//!   of `HLRBRep_PolyAlgo` is projected by the projector transform, split by
//!   `HLRAlgo_EdgeIterator` into its visible and hidden parts, and stored as
//!   [`BiPnt2D`] segments.
//! * [`HlrBRepPolyHlrToShape::internal_compound`] is `InternalCompound`
//!   (`cxx:83-212`): the four `typ` filters of `cxx:125-140` / `:168-183`
//!   (`1` = outlines, `2` = G1 lines, `3` = G2 lines, anything else = the
//!   default sharp/silhouette set) select the segments, and each surviving
//!   segment becomes a 2D edge built by `BRepLib_MakeEdge2d(P1, P2)`
//!   (`cxx:154` / `:201`).
//! * The twelve `VCompound` / `Rg1LineVCompound` / `RgNLineVCompound` /
//!   `OutLineVCompound` / `HCompound` / ... entry points of `lxx:21-143` are
//!   ported with the same `(typ, visible)` pairs.
//!
//! `BRepLib_MakeEdge2d(P1, P2)` (`BRepLib_MakeEdge2d.cxx:129-140`, i.e.
//! `gp_Lin2d L(P1, gp_Vec2d(P1, P2))` + `Init(GL, P1, P2, 0, l)`,
//! `BRepLib_MakeEdge2d.cxx:387-403` then `:465-633`) becomes [`make_edge_2d`]:
//! a `Geom2d_Line` on `BRepLib::Plane()` (`BRepLib.cxx:138-145`, the `gp::XOY()`
//! plane) carrying a **p-curve only** - the C++ edge has no 3D curve either,
//! because `B.UpdateEdge(E, C, BRepLib::Plane(), TopLoc_Location(), preci)`
//! (`BRepLib_MakeEdge2d.cxx:621`) writes a `CurveOnSurface` representation and
//! nothing else.
//!
//! The vertex/range body is [`occt_topo::brep_lib_make_edge::init_2d`], which ports
//! `BRepLib_MakeEdge::Init(C2d, S, VV1, VV2, p1, p2)`
//! (`BRepLib_MakeEdge.cxx:905-1094`). Two differences against
//! `BRepLib_MakeEdge2d::Init` (`BRepLib_MakeEdge2d.cxx:465-633`) hold and are
//! intentional:
//!
//! * the generic `Init` writes `E.Closed(closed)` (`BRepLib_MakeEdge.cxx:1092`),
//!   `BRepLib_MakeEdge2d::Init` does not - [`make_edge_2d`] never writes it
//!   either, matching the 2D factory;
//! * `BRepLib_MakeEdge2d::Point(const gp_Pnt2d&)` lifts a 2D point to
//!   `(x, y, 0)` directly, while the generic `Init` evaluates
//!   `P = S->Value(C->Value(p))` (`BRepLib_MakeEdge.cxx:968-981`). On the
//!   `gp::XOY()` plane the two agree, so the endpoint vertices land on exactly
//!   the same points.
//!
//! Other deviations:
//!
//! * `gp::XOY()` (`BRepLib.cxx:140-144`) is one process-wide plane in OCCT; the
//!   port keeps one [`OnceLock`] handle so every produced edge shares the same
//!   `Geom_Surface` pointer and therefore the same p-curve registry key.
//! * `myAlgo` is an `occ::handle` in OCCT; the port keeps a clone of the
//!   algorithm (`HLRBRep_PolyAlgo` is `Clone`) so `Show` mode can re-walk it
//!   after `Update` has returned. The cursors `InitHide`/`InitShow` drive are
//!   per-clone, so the caller's algorithm is left at whatever position it had -
//!   OCCT leaves `myAlgo` on the last hidden segment instead
//!   (`cxx:51-78`). Every reader of the cursor calls `InitHide` / `InitShow`
//!   first (`cxx:51`, `:165`), so the walk itself is unaffected.
//! * `BRep_Builder::Add(Result, <null edge>)` cannot happen: `cxx:152` / `:199`
//!   only build an edge when the squared distance of the two 2D points is above
//!   `1e-20`, far above `gp::Resolution()`, so `BRepLib_LineThroughIdenticPoints`
//!   (`BRepLib_MakeEdge2d.cxx:131-136`) is out of reach.

use std::collections::HashSet;
use std::sync::{Arc, OnceLock};

use occt_core::gp::{GpDir2d, GpPln, GpPnt2d, GpTrsf, GpVec2d};
use occt_core::precision::{CONFUSION, REAL_SMALL};
use occt_geom::plane::GeomPlane;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;
use occt_geom2d::line::Geom2dLine;

use occt_topo::abs::{Orientation, ShapeType};
use occt_topo::brep_lib_make_edge;
use occt_topo::builder::TopoBuilder;
use crate::hlr_algo::{HlrAlgoEdgeIterator, HlrAlgoEdgeStatus};
use crate::hlrbrep_poly_algo::{HlrBRepPolyAlgo, PolyAlgoSegment};
use occt_topo::shape::{Edge, TopoShape};
use occt_topo::tgeometry::GeometryRegistry;
use occt_topo::topo_tools::{explorer_shapes, shape_map_key, ShapeMapKey};

/// `HLRBRep_BiPnt2D` (`HLRBRep_BiPnt2D.hxx:27-110`): a projected segment with
/// the shape it came from and the four line-class flags of the filter.
#[derive(Debug, Clone)]
pub struct BiPnt2D {
    p1: GpPnt2d,
    p2: GpPnt2d,
    shape: TopoShape,
    rg1_line: bool,
    rgn_line: bool,
    out_line: bool,
    int_line: bool,
}

impl BiPnt2D {
    /// `HLRBRep_BiPnt2D(x1, y1, x2, y2, S, reg1, regn, outl, intl)`
    /// (`HLRBRep_BiPnt2D.hxx:40-56`).
    pub fn new(
        p1: GpPnt2d,
        p2: GpPnt2d,
        shape: TopoShape,
        reg1: bool,
        regn: bool,
        outl: bool,
        intl: bool,
    ) -> Self {
        Self {
            p1,
            p2,
            shape,
            rg1_line: reg1,
            rgn_line: regn,
            out_line: outl,
            int_line: intl,
        }
    }

    /// `P1()` (`hxx:76`).
    pub fn p1(&self) -> &GpPnt2d {
        &self.p1
    }

    /// `P2()` (`hxx:78`).
    pub fn p2(&self) -> &GpPnt2d {
        &self.p2
    }

    /// `Shape()` (`hxx:80`).
    pub fn shape(&self) -> &TopoShape {
        &self.shape
    }

    /// `Shape(S)` (`hxx:82`).
    pub fn set_shape(&mut self, s: TopoShape) {
        self.shape = s;
    }

    /// `Rg1Line()` (`hxx:84`).
    pub fn rg1_line(&self) -> bool {
        self.rg1_line
    }

    /// `Rg1Line(B)` (`hxx:86`).
    pub fn set_rg1_line(&mut self, b: bool) {
        self.rg1_line = b;
    }

    /// `RgNLine()` (`hxx:88`).
    pub fn rgn_line(&self) -> bool {
        self.rgn_line
    }

    /// `RgNLine(B)` (`hxx:90`).
    pub fn set_rgn_line(&mut self, b: bool) {
        self.rgn_line = b;
    }

    /// `OutLine()` (`hxx:92`).
    pub fn out_line(&self) -> bool {
        self.out_line
    }

    /// `OutLine(B)` (`hxx:94`).
    pub fn set_out_line(&mut self, b: bool) {
        self.out_line = b;
    }

    /// `IntLine()` (`hxx:96`).
    pub fn int_line(&self) -> bool {
        self.int_line
    }

    /// `IntLine(B)` (`hxx:98`).
    pub fn set_int_line(&mut self, b: bool) {
        self.int_line = b;
    }
}

/// `BRepLib::Plane()` (`BRepLib.cxx:138-145`): the shared `gp::XOY()` plane.
///
/// OCCT builds it once into a static handle (`BRepLib.cxx:140-144`); the port
/// does the same so the p-curve registry key (the `Geom_Surface` pointer,
/// [`GeometryRegistry::update_edge_pcurve_on_surface`]) is stable across calls.
pub(crate) fn brep_lib_plane() -> &'static Arc<dyn Surface> {
    static PLANE: OnceLock<Arc<dyn Surface>> = OnceLock::new();
    PLANE.get_or_init(|| Arc::new(GeomPlane::new(GpPln::default())) as Arc<dyn Surface>)
}

/// `BRepLib_MakeEdge2d(const gp_Pnt2d& P1, const gp_Pnt2d& P2)`
/// (`BRepLib_MakeEdge2d.cxx:129-140`) followed by the `Init` chain of
/// `:387-403` and `:465-633`.
///
/// `None` is `BRepLib_LineThroughIdenticPoints` (`cxx:131-136`). The edge
/// carries a p-curve on [`brep_lib_plane`] and no 3D curve, with the range
/// `[0, l]` that `B.Range(E, p1, p2)` stores (`cxx:630`).
pub fn make_edge_2d(p1: &GpPnt2d, p2: &GpPnt2d) -> Option<Edge> {
    // `cxx:131`: `l = P1.Distance(P2)` against `gp::Resolution()`.
    let l = p1.distance(p2);
    if l <= REAL_SMALL {
        return None;
    }
    // `cxx:137`: `gp_Lin2d L(P1, gp_Vec2d(P1, P2))`.
    let dir = GpDir2d::new(p2.x() - p1.x(), p2.y() - p1.y()).ok()?;
    let line = Arc::new(Geom2dLine::from_pnt_dir(*p1, dir)) as Arc<dyn Curve2d>;
    // `cxx:139` / `:387-403` / `:465-633`: the two vertices at
    // `BRepLib::Precision()`, then the range/vertex checks and the swap when
    // `p1 >= p2` (`:499-519`).
    let plane = brep_lib_plane();
    let build = brep_lib_make_edge::init_2d(&line, plane, 0.0, l)?;

    // `cxx:619-630`: `B.MakeEdge(E)` + `B.UpdateEdge(E, C, Plane, Loc, preci)`
    // + the two endpoint vertices + `B.Range(E, p1, p2)`.
    let mut edge = Edge::new();
    let reg = GeometryRegistry::global();
    reg.set_ranges_all(&edge.0, build.first, build.last);
    reg.update_edge_pcurve_on_surface(&edge.0, plane, Some(line), CONFUSION);
    reg.set_pcurve_range_on_surface(&edge.0, plane, build.first, build.last);
    let b = TopoBuilder::new();
    brep_lib_make_edge::add_endpoint_vertices(&b, &mut edge, &build);
    if build.reverse {
        // `BRepLib_MakeEdge.cxx:1089-1092`; `BRepLib_MakeEdge2d.cxx:619-633`
        // has no such write, but `p1 = 0 < l = p2` keeps the flag false here.
        edge.0.set_orientation(Orientation::Reversed);
    }
    Some(edge)
}

/// A framework for filtering the computation results of an
/// `HLRBRep_PolyAlgo` algorithm by extraction. Source:
/// `HLRBRep_PolyHLRToShape`.
#[derive(Debug, Clone, Default)]
pub struct HlrBRepPolyHlrToShape {
    algo: Option<HlrBRepPolyAlgo>,
    bi_pnt_vis: Vec<BiPnt2D>,
    bi_pnt_hid: Vec<BiPnt2D>,
    hide_mode: bool,
}

impl HlrBRepPolyHlrToShape {
    /// `HLRBRep_PolyHLRToShape()` (`cxx:33`): OCCT's defaulted constructor
    /// leaves `myHideMode` (`hxx:122`) to its default `false`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Update(A)` (`cxx:37-79`).
    pub fn update(&mut self, algo: &HlrBRepPolyAlgo) {
        // `cxx:39-40`: `myAlgo = A; myHideMode = true;`.
        self.algo = Some(algo.clone());
        self.hide_mode = true;
        // `cxx:44-45`: both lists are cleared.
        self.bi_pnt_vis.clear();
        self.bi_pnt_hid.clear();

        // `cxx:48`: `const gp_Trsf& T = myAlgo->Projector().Transformation()`.
        let t = algo.projector().transformation().clone();
        let mut a = algo.clone();

        // `cxx:51-78`: `for (InitHide(); MoreHide(); NextHide())`.
        a.init_hide();
        while a.more_hide() {
            // `cxx:53`: `myAlgo->Hide(status, S, reg1, regn, outl, intl)`. OCCT
            // declares `status` outside the loop (`cxx:49`) but only reads it
            // before the next `Hide`; the port scopes it to the iteration so
            // the iterator borrow ends with it.
            let mut status = HlrAlgoEdgeStatus::new();
            if let Some(seg) = a.hide(&mut status) {
                self.append_parts(seg, &status, &t);
            }
            a.next_hide();
        }
    }

    /// `cxx:54-75`: the body of the `Hide` loop.
    fn append_parts(&mut self, seg: PolyAlgoSegment, status: &HlrAlgoEdgeStatus, t: &GpTrsf) {
        // `cxx:54-57`: `aSta = aPoints.Pnt1; aEnd = aPoints.Pnt2; T.Transforms`.
        let mut a_sta = seg.points.pnt1;
        let mut a_end = seg.points.pnt2;
        t.transforms_xyz(&mut a_sta);
        t.transforms_xyz(&mut a_end);
        // `cxx:58-60`.
        let sta_2d = GpPnt2d::new(a_sta.x(), a_sta.y());
        let end_2d = GpPnt2d::new(a_end.x(), a_end.y());
        let a_d = GpVec2d::new(end_2d.x() - sta_2d.x(), end_2d.y() - sta_2d.y());
        // `cxx:61`: `if (aD.Modulus() > 1.e-10)`.
        if a_d.coord.modulus() > 1.0e-10 {
            // `cxx:64-68`.
            let mut it = HlrAlgoEdgeIterator::new();
            it.init_visible(status);
            while it.more_visible() {
                let (sta, _tol_sta, end, _tol_end) = it.visible();
                self.bi_pnt_vis.push(BiPnt2D::new(
                    offset(&sta_2d, &a_d, sta),
                    offset(&sta_2d, &a_d, end),
                    seg.shape.clone(),
                    seg.reg1,
                    seg.regn,
                    seg.outl,
                    seg.intl,
                ));
                it.next_visible();
            }
            // `cxx:71-75`.
            let mut it = HlrAlgoEdgeIterator::new();
            it.init_hidden(status);
            while it.more_hidden() {
                let (sta, _tol_sta, end, _tol_end) = it.hidden();
                self.bi_pnt_hid.push(BiPnt2D::new(
                    offset(&sta_2d, &a_d, sta),
                    offset(&sta_2d, &a_d, end),
                    seg.shape.clone(),
                    seg.reg1,
                    seg.regn,
                    seg.outl,
                    seg.intl,
                ));
                it.next_hidden();
            }
        }
    }

    /// `Show()` (`lxx:21-24`).
    pub fn show_mode(&mut self) {
        self.hide_mode = false;
    }

    /// `Hide()` (`lxx:28-31`).
    pub fn hide_mode(&mut self) {
        self.hide_mode = true;
    }

    /// `InternalCompound(typ, visible, S)` (`cxx:83-212`).
    ///
    /// `None` is the null `TopoDS_Shape` OCCT returns when nothing was added
    /// (`cxx:207-210`).
    pub fn internal_compound(
        &self,
        typ: i32,
        visible: bool,
        s: Option<&TopoShape>,
    ) -> Option<TopoShape> {
        // `cxx:87-99`: the `NCollection_Map` of the S edges and faces.
        let mut map: HashSet<ShapeMapKey> = HashSet::new();
        if let Some(s) = s {
            for e in explorer_shapes(s, ShapeType::Edge) {
                map.insert(shape_map_key(&e));
            }
            for f in explorer_shapes(s, ShapeType::Face) {
                map.insert(shape_map_key(&f));
            }
        }

        // `cxx:100-104`.
        let mut result = TopoShape::new(ShapeType::Compound);
        let builder = TopoBuilder::new();
        let mut added = false;

        if self.hide_mode {
            // `cxx:106-158`.
            let list = if visible {
                &self.bi_pnt_vis
            } else {
                &self.bi_pnt_hid
            };
            for bp in list {
                // `cxx:121-140`.
                let mut todraw =
                    draw_filter(typ, bp.rg1_line, bp.rgn_line, bp.out_line, bp.int_line);
                // `cxx:141-147`.
                if s.is_some() {
                    todraw = map.contains(&shape_map_key(bp.shape()));
                }
                if todraw {
                    // `cxx:148-156`: `FirstP2d.SquareDistance(LastP2d) > 1.e-20`.
                    if bp.p1().square_distance(bp.p2()) > 1.0e-20 {
                        if let Some(e) = make_edge_2d(bp.p1(), bp.p2()) {
                            builder.add(&mut result, &e.0);
                        }
                        // `cxx:155`: `added = true` regardless of the builder's
                        // verdict, exactly as in the source.
                        added = true;
                    }
                }
            }
        } else {
            // `cxx:160-206`.
            let Some(algo) = self.algo.as_ref() else {
                return None;
            };
            // `cxx:162`.
            let t = algo.projector().transformation().clone();
            let mut a = algo.clone();
            // `cxx:165-205`.
            a.init_show();
            while a.more_show() {
                let Some(seg) = a.show() else {
                    a.next_show();
                    continue;
                };
                // `cxx:168-183`.
                let mut todraw = draw_filter(typ, seg.reg1, seg.regn, seg.outl, seg.intl);
                // `cxx:184-190`.
                if s.is_some() {
                    todraw = map.contains(&shape_map_key(&seg.shape));
                }
                if todraw {
                    // `cxx:193-199`: project `Pnt1`/`Pnt2` and test the
                    // **squared** modulus against `1.e-20`.
                    let mut a_sta = seg.points.pnt1;
                    let mut a_end = seg.points.pnt2;
                    t.transforms_xyz(&mut a_sta);
                    t.transforms_xyz(&mut a_end);
                    let sta_2d = GpPnt2d::new(a_sta.x(), a_sta.y());
                    let end_2d = GpPnt2d::new(a_end.x(), a_end.y());
                    let a_d = GpVec2d::new(end_2d.x() - sta_2d.x(), end_2d.y() - sta_2d.y());
                    if a_d.coord.square_modulus() > 1.0e-20 {
                        if let Some(e) = make_edge_2d(&sta_2d, &end_2d) {
                            builder.add(&mut result, &e.0);
                        }
                        added = true;
                    }
                }
                a.next_show();
            }
        }

        // `cxx:207-210`.
        if !added {
            return None;
        }
        Some(result)
    }

    /// `VCompound()` (`lxx:35-38`).
    pub fn v_compound(&self) -> Option<TopoShape> {
        self.internal_compound(4, true, None)
    }

    /// `VCompound(S)` (`lxx:42-45`).
    pub fn v_compound_of(&self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(4, true, Some(s))
    }

    /// `Rg1LineVCompound()` (`lxx:49-52`).
    pub fn rg1_line_v_compound(&self) -> Option<TopoShape> {
        self.internal_compound(2, true, None)
    }

    /// `Rg1LineVCompound(S)` (`lxx:56-59`).
    pub fn rg1_line_v_compound_of(&self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(2, true, Some(s))
    }

    /// `RgNLineVCompound()` (`lxx:63-66`).
    pub fn rgn_line_v_compound(&self) -> Option<TopoShape> {
        self.internal_compound(3, true, None)
    }

    /// `RgNLineVCompound(S)` (`lxx:70-73`).
    pub fn rgn_line_v_compound_of(&self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(3, true, Some(s))
    }

    /// `OutLineVCompound()` (`lxx:77-80`).
    pub fn out_line_v_compound(&self) -> Option<TopoShape> {
        self.internal_compound(1, true, None)
    }

    /// `OutLineVCompound(S)` (`lxx:84-87`).
    pub fn out_line_v_compound_of(&self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(1, true, Some(s))
    }

    /// `HCompound()` (`lxx:91-94`).
    pub fn h_compound(&self) -> Option<TopoShape> {
        self.internal_compound(4, false, None)
    }

    /// `HCompound(S)` (`lxx:98-101`).
    pub fn h_compound_of(&self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(4, false, Some(s))
    }

    /// `Rg1LineHCompound()` (`lxx:105-108`).
    pub fn rg1_line_h_compound(&self) -> Option<TopoShape> {
        self.internal_compound(2, false, None)
    }

    /// `Rg1LineHCompound(S)` (`lxx:112-115`).
    pub fn rg1_line_h_compound_of(&self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(2, false, Some(s))
    }

    /// `RgNLineHCompound()` (`lxx:119-122`).
    pub fn rgn_line_h_compound(&self) -> Option<TopoShape> {
        self.internal_compound(3, false, None)
    }

    /// `RgNLineHCompound(S)` (`lxx:126-129`).
    pub fn rgn_line_h_compound_of(&self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(3, false, Some(s))
    }

    /// `OutLineHCompound()` (`lxx:133-136`).
    pub fn out_line_h_compound(&self) -> Option<TopoShape> {
        self.internal_compound(1, false, None)
    }

    /// `OutLineHCompound(S)` (`lxx:140-143`).
    pub fn out_line_h_compound_of(&self, s: &TopoShape) -> Option<TopoShape> {
        self.internal_compound(1, false, Some(s))
    }

    /// The visible segments `Update` recorded, in `Append` order.
    pub fn bi_pnt_vis(&self) -> &[BiPnt2D] {
        &self.bi_pnt_vis
    }

    /// The hidden segments `Update` recorded, in `Append` order.
    pub fn bi_pnt_hid(&self) -> &[BiPnt2D] {
        &self.bi_pnt_hid
    }

    /// `myHideMode` (`hxx:122`).
    pub fn is_hide_mode(&self) -> bool {
        self.hide_mode
    }
}

/// `aSta2D + k * aD` (`cxx:68`, `:75`).
fn offset(sta: &GpPnt2d, dir: &GpVec2d, k: f64) -> GpPnt2d {
    GpPnt2d::new(sta.x() + k * dir.x(), sta.y() + k * dir.y())
}

/// The four `todraw` arms of `cxx:125-140`, repeated verbatim at
/// `cxx:168-183`.
fn draw_filter(typ: i32, reg1: bool, regn: bool, outl: bool, intl: bool) -> bool {
    match typ {
        // `cxx:125-128` / `:168-171`.
        1 => intl,
        // `cxx:129-132` / `:172-175`.
        2 => reg1 && !regn && !outl,
        // `cxx:133-136` / `:176-179`.
        3 => regn && !outl,
        // `cxx:137-140` / `:180-183`.
        _ => !intl && (!reg1 || outl),
    }
}
