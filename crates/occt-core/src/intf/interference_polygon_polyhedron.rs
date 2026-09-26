//! Port of `Intf_InterferencePolygonPolyhedron`
//! (`src/ModelingAlgorithms/TKGeomAlgo/Intf/Intf_InterferencePolygonPolyhedron.gxx`),
//! the generic 3D interference engine behind
//! `IntCurveSurface_TheInterferenceOfHInter` and
//! `HLRBRep_TheInterferenceOfInterCSurf`.
//!
//! The OCCT file is a `.gxx` template: every `ToolPolygon3d::` /
//! `ToolPolyh::` static call is resolved by the including `_0.cxx`.
//! The port keeps that split as two traits, following the style already used by
//! `intf::polygon2d::IntfPolygon2d`:
//! * [`IntfPolygon3dTool`] mirrors
//!   `IntCurveSurface_ThePolygonToolOfHInter`
//!   (`IntCurveSurface_ThePolygonToolOfHInter.hxx:33-65`).
//! * [`IntfPolyhedronTool`] mirrors
//!   `IntCurveSurface_ThePolyhedronToolOfHInter`
//!   (`IntCurveSurface_ThePolyhedronToolOfHInter.hxx:36-109`).
//! The two `Intersect` overloads have identical bodies in OCCT (the second
//! is the first with the plane data passed in, `gxx:1071-1075`); the port
//! keeps one body and lets the first overload compute the plane and delegate.
//!
//! Ported members (OCCT 8.0.0, tag V8_0_0):
//! * default ctor                       `gxx:60-65`
//! * polygon/polyhedron ctors           `gxx:73-106`
//! * line/polyhedron ctor               `gxx:114-147`
//! * lines/polyhedron ctor              `gxx:155-193`
//! * `Perform` x3                    `gxx:197-288`
//! * `Interference` (polyg, polyh)   `gxx:296-350`
//! * plane-receiving ctor overloads     `gxx:358-430`
//! * `Perform` x3 (with grid)        `gxx:434-523`
//! * `Interference` (polyg, polyh, grid) `gxx:531-623`
//! * `Intersect` (plane recomputed)  `gxx:741-1053`
//! * `Intersect` (plane supplied)    `gxx:1057-1368`
//!
//! UNPORTED in this module:
//! * the `#if 0` dead overload (`gxx:630-739`) - dead code in OCCT
//!   itself; it is not compiled there either.
//! * the trailing proximity block of both `Intersect` overloads,
//!   `gxx:981-1052` and `gxx:1296-1367`: it builds an
//!   `Extrema_ExtElC(LinPol, LinTri, 1e-8)` and derives extra section points
//!   from `SquareDistance` / `Points` / `IsInSegment`. Neither
//!   `Extrema_ExtElC` nor `Extrema_POnCurv` is ported in this crate
//!   (grep for `Extrema_ExtElC` returns nothing), so the block is skipped
//!   and marked here instead of being approximated. `IsInSegment`
//!   (`gxx:36-56`) has no other caller and is likewise not ported.
//!
//! `Bnd_BoundSortBox` (`gxx:127, 171, 224, 266, 303`) is
//! [`IntfPolyhGrid`] below, which now delegates to
//! [`crate::bnd::BoundSortBox`] -- the faithful port of `Bnd_BoundSortBox`
//! (`Bnd_BoundSortBox.cxx`) added by task T-89. Earlier this file carried a
//! brute-force surrogate because the crate's `bnd::sortbox` was a self-made
//! grid that returned an empty candidate set for a component box larger than
//! one cell; that implementation is gone. Components-bounding indices are
//! 1-based, as in the OCCT `NCollection_HArray1` handed to
//! `Bnd_BoundSortBox::Initialize`.

use std::marker::PhantomData;

use crate::bnd::BndBox;
use crate::gp::{GpDir, GpLin, GpPnt, GpVec, GpXyz};
use crate::precision::{epsilon, REAL_SMALL};

use super::interference::IntfInterference;
use super::plane::plane_equation;
use super::section_point::{IntfPIType, IntfSectionPoint};
use super::tool::IntfTool;

/// `Pourcent3[4]` (`gxx:34`).
const POURCENT3: [usize; 4] = [0, 1, 2, 0];

/// `ToolPolygon3d` - the static polygon services the gxx calls.
/// Compare `IntCurveSurface_ThePolygonToolOfHInter.hxx:33-65`.
///
/// All indices are 1-based, as in OCCT (`gxx:310, 175`).
pub trait IntfPolygon3dTool {
    /// The polygon type (`Polygon3d` in the gxx).
    type Polygon3d;

    /// `Bounding` (`IntCurveSurface_ThePolygonToolOfHInter.hxx:33-36`).
    fn bounding(the_polyg: &Self::Polygon3d) -> &BndBox;

    /// `DeflectionOverEstimation` (`...hxx:38-41`).
    fn deflection_over_estimation(the_polyg: &Self::Polygon3d) -> f64;

    /// `Closed` (`...hxx:43-46`).
    fn closed(the_polyg: &Self::Polygon3d) -> bool;

    /// `NbSegments` (`...hxx:48-51`).
    fn nb_segments(the_polyg: &Self::Polygon3d) -> usize;

    /// `BeginOfSeg(index)` (`...hxx:54-58`), 1-based.
    fn begin_of_seg(the_polyg: &Self::Polygon3d, index: usize) -> GpPnt;

    /// `EndOfSeg(index)` (`...hxx:61-65`), 1-based.
    fn end_of_seg(the_polyg: &Self::Polygon3d, index: usize) -> GpPnt;
}

/// `ToolPolyh` - the static polyhedron services the gxx calls.
/// Compare `IntCurveSurface_ThePolyhedronToolOfHInter.hxx:36-109`.
pub trait IntfPolyhedronTool {
    /// The polyhedron type (`Polyhedron` in the gxx).
    type Polyhedron;

    /// `Bounding` (`...hxx:36-39`).
    fn bounding(the_polyh: &Self::Polyhedron) -> &BndBox;

    /// `ComponentsBounding` (`...hxx:43-47`): one box per triangle.
    fn components_bounding(the_polyh: &Self::Polyhedron) -> &[BndBox];

    /// `DeflectionOverEstimation` (`...hxx:50-53`).
    fn deflection_over_estimation(the_polyh: &Self::Polyhedron) -> f64;

    /// `NbTriangles` (`...hxx:56-59`).
    fn nb_triangles(the_polyh: &Self::Polyhedron) -> usize;

    /// `Triangle(Index, P1, P2, P3)` (`...hxx:63-70`).
    fn triangle(the_polyh: &Self::Polyhedron, index: i32) -> (i32, i32, i32);

    /// `Point(Index)` (`...hxx:73-76`).
    fn point(the_polyh: &Self::Polyhedron, index: i32) -> GpPnt;

    /// `TriConnex` (`...hxx:84-92`). Returns
    /// `(return_value, TriCon, OtherP)`; the gxx discards all three
    /// (`gxx:915, 1230`).
    fn tri_connex(
        the_polyh: &Self::Polyhedron,
        triang: i32,
        pivot: i32,
        pedge: i32,
    ) -> (i32, i32, i32);

    /// `IsOnBound(Index1, Index2)` (`...hxx:98-103`).
    fn is_on_bound(the_polyh: &Self::Polyhedron, index1: i32, index2: i32) -> bool;

    /// `GetBorderDeflection` (`...hxx:106-109`).
    fn get_border_deflection(the_polyh: &Self::Polyhedron) -> f64;
}

/// The `Bnd_BoundSortBox` the gxx builds
/// (`gxx:126-127, 170-171, 223-224, 265-266, 302-303`), holding the
/// `ComponentsBounding` array and answering with the boxes touched by the
/// query -- [`crate::bnd::BoundSortBox`], the faithful port of
/// `Bnd_BoundSortBox` (`Bnd_BoundSortBox.cxx`), so `Compare` is the OCCT
/// voxel-grid voxel candidate set (`cxx:410-504`) rather than a re-scan.
///
/// Indices returned by [`IntfPolyhGrid::compare`] are 1-based, so they can
/// be passed straight to `ToolPolyh::Triangle` / `Point`, as in OCCT.
#[derive(Clone, Debug, Default)]
pub struct IntfPolyhGrid {
    /// `myEnclosingBox` + `myBoxes` (`Bnd_BoundSortBox.hxx:139-147`).
    /// `BoundSortBox::Compare` is `&mut self` in the port (it caches the last
    /// result, `myLastResult`), and the gxx passes the grid as a shared
    /// reference, so the cell is what bridges the two.
    grid: std::cell::RefCell<crate::bnd::BoundSortBox>,
}

impl IntfPolyhGrid {
    /// `Bnd_BoundSortBox::Initialize(enclosing, setOfBoxes)`
    /// (`Bnd_BoundSortBox.cxx:352-368`).
    pub fn initialize<TH: IntfPolyhedronTool>(the_polyh: &TH::Polyhedron) -> Self {
        Self {
            grid: {
                let mut g = crate::bnd::BoundSortBox::new();
                g.initialize_with_enclosing(
                    TH::bounding(the_polyh),
                    &TH::components_bounding(the_polyh).to_vec(),
                );
                std::cell::RefCell::new(g)
            },
        }
    }

    /// `Bnd_BoundSortBox::Compare(const Bnd_Box&)`
    /// (`Bnd_BoundSortBox.cxx:410-504`): the boxes touched by `the_box`.
    /// OCCT's indices are 1-based; `BoundSortBox` already returns them that
    /// way, so only the width is narrowed here.
    pub fn compare(&self, the_box: &BndBox) -> Vec<i32> {
        self.grid.borrow_mut().compare(the_box).iter().map(|i| *i as i32).collect()
    }
}

/// `Intf_InterferencePolygonPolyhedron`
/// (`IntCurveSurface_TheInterferenceOfHInter.hxx:37-152`).
///
/// The OCCT class derives from `Intf_Interference`; the port composes the
/// base as [`IntfInterference`] (the same choice as the 2D port).
#[derive(Clone, Debug)]
pub struct IntfInterferencePolygonPolyhedron<TP: IntfPolygon3dTool, TH: IntfPolyhedronTool> {
    /// The `Intf_Interference` base sub-object.
    pub base: IntfInterference,
    /// `BeginOfClosedPolygon` (`...hxx:150`).
    begin_of_closed_polygon: bool,
    /// `iLin` (`...hxx:151`), 1-based polygon segment (0 for lines).
    i_lin: i32,
    _phantom: PhantomData<(TP, TH)>,
}

impl<TP: IntfPolygon3dTool, TH: IntfPolyhedronTool> Default
    for IntfInterferencePolygonPolyhedron<TP, TH>
{
    fn default() -> Self {
        Self::new()
    }
}

impl<TP: IntfPolygon3dTool, TH: IntfPolyhedronTool>
    IntfInterferencePolygonPolyhedron<TP, TH>
{
    /// `Intf_InterferencePolygonPolyhedron()` (`gxx:60-65`).
    pub fn new() -> Self {
        Self {
            base: IntfInterference::new(false),
            begin_of_closed_polygon: false,
            i_lin: 0,
            _phantom: PhantomData,
        }
    }

    // ------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------

    /// `Intf_InterferencePolygonPolyhedron(thePolyg, thePolyh)`
    /// (`gxx:73-88`).
    pub fn with_polygon_polyhedron(
        the_polyg: &TP::Polygon3d,
        the_polyh: &TH::Polyhedron,
    ) -> Self {
        let mut me = Self::new();
        me.perform_polygon_polyhedron(the_polyg, the_polyh);
        me
    }

    /// `Intf_InterferencePolygonPolyhedron(thePolyg, thePolyh, PolyhGrid)`
    /// (`gxx:90-106`).
    pub fn with_polygon_polyhedron_grid(
        the_polyg: &TP::Polygon3d,
        the_polyh: &TH::Polyhedron,
        polyh_grid: &IntfPolyhGrid,
    ) -> Self {
        let mut me = Self::new();
        me.perform_polygon_polyhedron_grid(the_polyg, the_polyh, polyh_grid);
        me
    }

    /// `Intf_InterferencePolygonPolyhedron(theLin, thePolyh)`
    /// (`gxx:114-147`).
    pub fn with_lin_polyhedron(the_lin: &GpLin, the_polyh: &TH::Polyhedron) -> Self {
        let mut me = Self::new();
        me.perform_lin_polyhedron(the_lin, the_polyh);
        me
    }

    /// `Intf_InterferencePolygonPolyhedron(theLin, thePolyh, PolyhGrid)`
    /// (`gxx:358-388`).
    pub fn with_lin_polyhedron_grid(
        the_lin: &GpLin,
        the_polyh: &TH::Polyhedron,
        polyh_grid: &IntfPolyhGrid,
    ) -> Self {
        let mut me = Self::new();
        me.perform_lin_polyhedron_grid(the_lin, the_polyh, polyh_grid);
        me
    }

    /// `Intf_InterferencePolygonPolyhedron(theLins, thePolyh)`
    /// (`gxx:155-193`).
    pub fn with_lins_polyhedron(the_lins: &[GpLin], the_polyh: &TH::Polyhedron) -> Self {
        let mut me = Self::new();
        me.perform_lins_polyhedron(the_lins, the_polyh);
        me
    }

    /// `Intf_InterferencePolygonPolyhedron(theLins, thePolyh, PolyhGrid)`
    /// (`gxx:396-430`).
    pub fn with_lins_polyhedron_grid(
        the_lins: &[GpLin],
        the_polyh: &TH::Polyhedron,
        polyh_grid: &IntfPolyhGrid,
    ) -> Self {
        let mut me = Self::new();
        me.perform_lins_polyhedron_grid(the_lins, the_polyh, polyh_grid);
        me
    }

    // ------------------------------------------------------------------
    // Perform
    // ------------------------------------------------------------------

    /// `Perform(thePolyg, thePolyh)` (`gxx:197-210`).
    pub fn perform_polygon_polyhedron(
        &mut self,
        the_polyg: &TP::Polygon3d,
        the_polyh: &TH::Polyhedron,
    ) {
        self.base.self_interference(false);
        self.base.tolerance = TP::deflection_over_estimation(the_polyg)
            + TH::deflection_over_estimation(the_polyh);
        if self.base.tolerance == 0.0 {
            self.base.tolerance = epsilon(1000.0);
        }

        if !TP::bounding(the_polyg).is_out_box(TH::bounding(the_polyh)) {
            self.interference_polygon_polyhedron(the_polyg, the_polyh);
        }
    }

    /// `Perform(theLin, thePolyh)` (`gxx:214-245`).
    pub fn perform_lin_polyhedron(&mut self, the_lin: &GpLin, the_polyh: &TH::Polyhedron) {
        self.base.self_interference(false);
        self.base.tolerance = TH::deflection_over_estimation(the_polyh);
        if self.base.tolerance == 0.0 {
            self.base.tolerance = epsilon(1000.0);
        }

        self.begin_of_closed_polygon = false;

        let polyh_grid = IntfPolyhGrid::initialize::<TH>(the_polyh);

        self.i_lin = 0;

        let mut bof_lin = BndBox::new();
        let mut btoo = IntfTool::new();
        btoo.lin_box(the_lin, TH::bounding(the_polyh), &mut bof_lin);

        for ind_tri in polyh_grid.compare(&bof_lin) {
            let p0 = the_lin.location();
            let dir = the_lin.direction();
            let p1 = p0.translated_vec(&GpVec::new(dir.x(), dir.y(), dir.z()));
            self.intersect(&p0, &p1, true, ind_tri as i32, the_polyh);
        }
    }

    /// `Perform(theLins, thePolyh)` (`gxx:253-288`).
    pub fn perform_lins_polyhedron(&mut self, the_lins: &[GpLin], the_polyh: &TH::Polyhedron) {
        self.base.self_interference(false);
        self.base.tolerance = TH::deflection_over_estimation(the_polyh);
        if self.base.tolerance == 0.0 {
            self.base.tolerance = epsilon(1000.0);
        }

        let mut bof_lin = BndBox::new();
        let mut btoo = IntfTool::new();
        self.begin_of_closed_polygon = false;

        let polyh_grid = IntfPolyhGrid::initialize::<TH>(the_polyh);

        for i_lin in 1..=the_lins.len() {
            self.i_lin = i_lin as i32;

            btoo.lin_box(&the_lins[i_lin - 1], TH::bounding(the_polyh), &mut bof_lin);

            for ind_tri in polyh_grid.compare(&bof_lin) {
                let p0 = the_lins[i_lin - 1].location();
                let dir = the_lins[i_lin - 1].direction();
                let p1 = p0.translated_vec(&GpVec::new(dir.x(), dir.y(), dir.z()));
                self.intersect(&p0, &p1, true, ind_tri as i32, the_polyh);
            }
        }
    }

    /// `Perform(thePolyg, thePolyh, PolyhGrid)` (`gxx:434-448`).
    pub fn perform_polygon_polyhedron_grid(
        &mut self,
        the_polyg: &TP::Polygon3d,
        the_polyh: &TH::Polyhedron,
        polyh_grid: &IntfPolyhGrid,
    ) {
        self.base.self_interference(false);
        self.base.tolerance = TP::deflection_over_estimation(the_polyg)
            + TH::deflection_over_estimation(the_polyh);
        if self.base.tolerance == 0.0 {
            self.base.tolerance = epsilon(1000.0);
        }

        if !TP::bounding(the_polyg).is_out_box(TH::bounding(the_polyh)) {
            self.interference_polygon_polyhedron_grid(the_polyg, the_polyh, polyh_grid);
        }
    }

    /// `Perform(theLin, thePolyh, PolyhGrid)` (`gxx:452-482`).
    pub fn perform_lin_polyhedron_grid(
        &mut self,
        the_lin: &GpLin,
        the_polyh: &TH::Polyhedron,
        polyh_grid: &IntfPolyhGrid,
    ) {
        self.base.self_interference(false);
        self.base.tolerance = TH::deflection_over_estimation(the_polyh);
        if self.base.tolerance == 0.0 {
            self.base.tolerance = epsilon(1000.0);
        }

        self.begin_of_closed_polygon = false;

        self.i_lin = 0;

        let mut bof_lin = BndBox::new();
        let mut btoo = IntfTool::new();
        btoo.lin_box(the_lin, TH::bounding(the_polyh), &mut bof_lin);

        for ind_tri in polyh_grid.compare(&bof_lin) {
            let p0 = the_lin.location();
            let dir = the_lin.direction();
            let p1 = p0.translated_vec(&GpVec::new(dir.x(), dir.y(), dir.z()));
            self.intersect(&p0, &p1, true, ind_tri as i32, the_polyh);
        }
    }

    /// `Perform(theLins, thePolyh, PolyhGrid)` (`gxx:490-523`).
    pub fn perform_lins_polyhedron_grid(
        &mut self,
        the_lins: &[GpLin],
        the_polyh: &TH::Polyhedron,
        polyh_grid: &IntfPolyhGrid,
    ) {
        self.base.self_interference(false);
        self.base.tolerance = TH::deflection_over_estimation(the_polyh);
        if self.base.tolerance == 0.0 {
            self.base.tolerance = epsilon(1000.0);
        }

        let mut bof_lin = BndBox::new();
        let mut btoo = IntfTool::new();
        self.begin_of_closed_polygon = false;

        for i_lin in 1..=the_lins.len() {
            self.i_lin = i_lin as i32;

            btoo.lin_box(&the_lins[i_lin - 1], TH::bounding(the_polyh), &mut bof_lin);

            for ind_tri in polyh_grid.compare(&bof_lin) {
                let p0 = the_lins[i_lin - 1].location();
                let dir = the_lins[i_lin - 1].direction();
                let p1 = p0.translated_vec(&GpVec::new(dir.x(), dir.y(), dir.z()));
                self.intersect(&p0, &p1, true, ind_tri as i32, the_polyh);
            }
        }
    }

    // ------------------------------------------------------------------
    // Interference
    // ------------------------------------------------------------------

    /// `Interference(thePolyg, thePolyh)` (`gxx:296-350`).
    fn interference_polygon_polyhedron(
        &mut self,
        the_polyg: &TP::Polygon3d,
        the_polyh: &TH::Polyhedron,
    ) {
        let polyh_grid = IntfPolyhGrid::initialize::<TH>(the_polyh);

        self.begin_of_closed_polygon = TP::closed(the_polyg);

        let def_ph = TH::deflection_over_estimation(the_polyh);

        for i_lin in 1..=TP::nb_segments(the_polyg) {
            self.i_lin = i_lin as i32;

            let mut bof_seg = BndBox::new();
            bof_seg.set_void();
            bof_seg.add_point(&TP::begin_of_seg(the_polyg, i_lin));
            bof_seg.add_point(&TP::end_of_seg(the_polyg, i_lin));
            bof_seg.enlarge(TP::deflection_over_estimation(the_polyg));

            let maliste = polyh_grid.compare(&bof_seg);
            for ind_tri in maliste {
                let p1 = TP::begin_of_seg(the_polyg, i_lin);
                let p2 = TP::end_of_seg(the_polyg, i_lin);
                let (t0, t1, t2) = TH::triangle(the_polyh, ind_tri as i32);
                let pa = TH::point(the_polyh, t0);
                let pb = TH::point(the_polyh, t1);
                let pc = TH::point(the_polyh, t2);
                let pa_pb = GpVec::from_pnts(&pa, &pb);
                let pa_pc = GpVec::from_pnts(&pa, &pc);
                let mut normale = pa_pb.crossed(&pa_pc);
                let norm_normale = normale.magnitude();
                if norm_normale < 1e-14 {
                    continue;
                }
                normale = normale.multiplied_scalar(def_ph / norm_normale);
                let p1m = p1.translated_vec(&normale.reversed());
                let p1p = p1.translated_vec(&normale);
                let p2m = p2.translated_vec(&normale.reversed());
                let p2p = p2.translated_vec(&normale);
                self.intersect(&p1m, &p2p, false, ind_tri as i32, the_polyh);
                self.intersect(&p1p, &p2m, false, ind_tri as i32, the_polyh);
            }
            self.begin_of_closed_polygon = false;
        }
    }

    /// `Interference(thePolyg, thePolyh, PolyhGrid)` (`gxx:531-623`,
    /// the "Modified by MKK - Thu Oct 25 12:40:11 2007" branch).
    fn interference_polygon_polyhedron_grid(
        &mut self,
        the_polyg: &TP::Polygon3d,
        the_polyh: &TH::Polyhedron,
        polyh_grid: &IntfPolyhGrid,
    ) {
        self.begin_of_closed_polygon = TP::closed(the_polyg);

        for i_lin in 1..=TP::nb_segments(the_polyg) {
            self.i_lin = i_lin as i32;

            let mut bof_seg = BndBox::new();
            bof_seg.set_void();
            bof_seg.add_point(&TP::begin_of_seg(the_polyg, i_lin));
            bof_seg.add_point(&TP::end_of_seg(the_polyg, i_lin));
            bof_seg.enlarge(TP::deflection_over_estimation(the_polyg));

            let def_ph = TH::deflection_over_estimation(the_polyh);
            let maliste = polyh_grid.compare(&bof_seg);

            let mut p1 = GpPnt::new(0.0, 0.0, 0.0);
            let mut p2 = GpPnt::new(0.0, 0.0, 0.0);
            let mut beg0 = GpPnt::new(0.0, 0.0, 0.0);
            let mut end0 = GpPnt::new(0.0, 0.0, 0.0);
            if !maliste.is_empty() {
                p1 = TP::begin_of_seg(the_polyg, i_lin);
                p2 = TP::end_of_seg(the_polyg, i_lin);
                beg0 = p1;
                end0 = p2;
            }

            for ind_tri in maliste {
                let (t0, t1, t2) = TH::triangle(the_polyh, ind_tri as i32);
                let p_tri = [t0, t1, t2];
                let pa = TH::point(the_polyh, p_tri[0]);
                let pb = TH::point(the_polyh, p_tri[1]);
                let pc = TH::point(the_polyh, p_tri[2]);
                let (tri_nor, tri_dp) = plane_equation(&pa, &pb, &pc);

                // enlarge boundary segment
                if i_lin == 1 {
                    let dif = p1.xyz().subtracted(p2.xyz());
                    let dist = dif.modulus();
                    if dist > REAL_SMALL {
                        let dif = dif.divided(dist);
                        let a_cos = dif.dot(&tri_nor).abs();
                        if a_cos > REAL_SMALL {
                            let shift = def_ph / a_cos;
                            beg0 = GpPnt::from_xyz(&p1.xyz().added(&dif.multiplied(shift)));
                        }
                    }
                } else if i_lin == TP::nb_segments(the_polyg) {
                    let dif = p2.xyz().subtracted(p1.xyz());
                    let dist = dif.modulus();
                    if dist > REAL_SMALL {
                        let dif = dif.divided(dist);
                        let a_cos = dif.dot(&tri_nor).abs();
                        if a_cos > REAL_SMALL {
                            let shift = def_ph / a_cos;
                            end0 = GpPnt::from_xyz(&p2.xyz().added(&dif.multiplied(shift)));
                        }
                    }
                }
                let d_beg_tri = tri_nor.dot(beg0.xyz()) - tri_dp;
                let d_end_tri = tri_nor.dot(end0.xyz()) - tri_dp;

                self.intersect_with_plane(
                    &beg0, &end0, false, ind_tri as i32, the_polyh, &tri_nor, tri_dp, d_beg_tri,
                    d_end_tri,
                );
            }
            self.begin_of_closed_polygon = false;
        }
    }

    // ------------------------------------------------------------------
    // Intersect
    // ------------------------------------------------------------------

    /// `Intersect(BegO, EndO, Infinite, TTri, thePolyh)`
    /// (`gxx:741-1053`).
    ///
    /// OCCT duplicates the whole body in this overload and in the
    /// plane-supplying one (`gxx:1057-1368`); the only difference is that
    /// the latter receives `TriNormal` / `TriDp` / `dBegTri` /
    /// `dEndTri` from its caller (`gxx:1071-1075`). The port folds
    /// them into [`Self::intersect_with_plane`].
    fn intersect(
        &mut self,
        beg_o: &GpPnt,
        end_o: &GpPnt,
        infinite: bool,
        t_tri: i32,
        the_polyh: &TH::Polyhedron,
    ) {
        let (t0, t1, t2) = TH::triangle(the_polyh, t_tri);
        let pa = TH::point(the_polyh, t0);
        let pb = TH::point(the_polyh, t1);
        let pc = TH::point(the_polyh, t2);
        let (tri_nor, tri_dp) = plane_equation(&pa, &pb, &pc);
        let d_beg_tri = tri_nor.dot(beg_o.xyz()) - tri_dp;
        let d_end_tri = tri_nor.dot(end_o.xyz()) - tri_dp;
        self.intersect_with_plane(
            beg_o, end_o, infinite, t_tri, the_polyh, &tri_nor, tri_dp, d_beg_tri, d_end_tri,
        );
    }

    /// `Intersect(BegO, EndO, Infinite, TTri, thePolyh, TriNormal, TriDp,
    /// dBegTri, dEndTri)` (`gxx:1057-1368`).
    ///
    /// UNPORTED: the trailing proximity block `gxx:1296-1367` (see the
    /// module header).
    #[allow(clippy::too_many_arguments)]
    fn intersect_with_plane(
        &mut self,
        beg_o: &GpPnt,
        end_o: &GpPnt,
        infinite: bool,
        t_tri: i32,
        the_polyh: &TH::Polyhedron,
        tri_normal: &GpXyz,
        _tri_dp: f64,
        d_beg_tri: f64,
        d_end_tri: f64,
    ) {
        let mut typ_on_g = IntfPIType::Edge;
        let (t0, t1, t2) = TH::triangle(the_polyh, t_tri);
        let p_tri = [t0, t1, t2];
        let tri_nor = *tri_normal;
        let mut no_intersection_with_triangle = false;

        let t = d_beg_tri - d_end_tri;
        let mut param = if t >= 1.0e-16 || t <= -1.0e-16 {
            d_beg_tri / t
        } else {
            d_beg_tri
        };
        let floatgap = epsilon(1000.0);

        if !infinite {
            if d_beg_tri <= floatgap && d_beg_tri >= -floatgap {
                param = 0.0;
                typ_on_g = IntfPIType::Vertex;
                if self.begin_of_closed_polygon {
                    no_intersection_with_triangle = false;
                }
            } else if d_end_tri <= floatgap && d_end_tri >= -floatgap {
                param = 1.0;
                typ_on_g = IntfPIType::Vertex;
                no_intersection_with_triangle = false;
            }
            if param < 0.0 || param > 1.0 {
                no_intersection_with_triangle = true;
            }
        }

        if !no_intersection_with_triangle {
            let sp_lieu = beg_o
                .xyz()
                .added(&end_o.xyz().subtracted(beg_o.xyz()).multiplied(param));
            let mut d_pi_e = [0.0f64; 3];
            let mut d_pt_pi = [0.0f64; 3];
            let mut is = 0i32;
            let mut s_edge = -1i32;
            let mut s_vertex = -1i32;
            let mut tbreak = 0i32;

            { // is = 0
                let pt0 = TH::point(the_polyh, p_tri[0]);
                let pt1 = TH::point(the_polyh, p_tri[1]);
                let seg_t = pt1.xyz().subtracted(pt0.xyz());
                let vec_p = sp_lieu.subtracted(pt0.xyz());
                d_pt_pi[0] = vec_p.modulus();
                if d_pt_pi[0] <= floatgap {
                    s_vertex = 0;
                    is = 0;
                    tbreak = 1;
                } else {
                    let seg_t_x_vec_p = seg_t.crossed(&vec_p);
                    let modulus_seg_t_x_vec_p = seg_t_x_vec_p.modulus();
                    let mut sigd = seg_t_x_vec_p.dot(&tri_nor);
                    if sigd > floatgap {
                        sigd = 1.0;
                    } else if sigd < -floatgap {
                        sigd = -1.0;
                    } else {
                        sigd = 0.0;
                    }
                    d_pi_e[0] = sigd * (modulus_seg_t_x_vec_p / seg_t.modulus());
                    if d_pi_e[0] <= floatgap && d_pi_e[0] >= -floatgap {
                        s_edge = 0;
                        is = 0;
                        tbreak = 1;
                    }
                }
            }

            if tbreak == 0 { // is = 1
                let pt1 = TH::point(the_polyh, p_tri[1]);
                let pt2 = TH::point(the_polyh, p_tri[2]);
                let seg_t = pt2.xyz().subtracted(pt1.xyz());
                let vec_p = sp_lieu.subtracted(pt1.xyz());
                d_pt_pi[1] = vec_p.modulus();
                if d_pt_pi[1] <= floatgap {
                    s_vertex = 1;
                    is = 1;
                    tbreak = 1;
                } else {
                    let seg_t_x_vec_p = seg_t.crossed(&vec_p);
                    let modulus_seg_t_x_vec_p = seg_t_x_vec_p.modulus();
                    let mut sigd = seg_t_x_vec_p.dot(&tri_nor);
                    if sigd > floatgap {
                        sigd = 1.0;
                    } else if sigd < -floatgap {
                        sigd = -1.0;
                    } else {
                        sigd = 0.0;
                    }
                    d_pi_e[1] = sigd * (modulus_seg_t_x_vec_p / seg_t.modulus());
                    if d_pi_e[1] <= floatgap && d_pi_e[1] >= -floatgap {
                        s_edge = 1;
                        is = 1;
                        tbreak = 1;
                    }
                }
            }
            if tbreak == 0 { // is = 2
                let pt2 = TH::point(the_polyh, p_tri[2]);
                let pt0 = TH::point(the_polyh, p_tri[0]);
                let seg_t = pt0.xyz().subtracted(pt2.xyz());
                let vec_p = sp_lieu.subtracted(pt2.xyz());
                d_pt_pi[2] = vec_p.modulus();
                if d_pt_pi[2] <= floatgap {
                    s_vertex = 2;
                    is = 2;
                }
                let seg_t_x_vec_p = seg_t.crossed(&vec_p);
                let modulus_seg_t_x_vec_p = seg_t_x_vec_p.modulus();
                let mut sigd = seg_t_x_vec_p.dot(&tri_nor);
                if sigd > floatgap {
                    sigd = 1.0;
                } else if sigd < -floatgap {
                    sigd = -1.0;
                } else {
                    sigd = 0.0;
                }
                d_pi_e[2] = sigd * (modulus_seg_t_x_vec_p / seg_t.modulus());
                if d_pi_e[2] <= floatgap && d_pi_e[2] >= -floatgap {
                    s_edge = 2;
                    is = 2;
                }
            }
            let _ = d_pt_pi;
            // fin for i=0 to 2

            if s_vertex > -1 {
                // OCCT: triCon = TTri; pedg = pTri[Pourcent3[sVertex + 1]];
                // (the TriConnex call is commented out at gxx:905-909 and both
                // values are unused).
                let _ = p_tri[POURCENT3[(s_vertex + 1) as usize]];
                let sp = IntfSectionPoint::with_3d(
                    &GpPnt::from_xyz(&sp_lieu),
                    typ_on_g,
                    0,
                    self.i_lin,
                    param,
                    IntfPIType::Vertex,
                    p_tri[is as usize],
                    0,
                    0.0,
                    1.0,
                );
                self.base.my_s_poins.push(sp);
            } else if s_edge > -1 {
                let pivot = p_tri[s_edge as usize];
                let pedge = p_tri[POURCENT3[(s_edge + 1) as usize]];
                // gxx:915 / 1230 discard the return value and both out-params
                // (triCon, pedg), so none of them is bound here.
                let _ = TH::tri_connex(the_polyh, t_tri, pivot, pedge);
                let sp = IntfSectionPoint::with_3d(
                    &GpPnt::from_xyz(&sp_lieu),
                    typ_on_g,
                    0,
                    self.i_lin,
                    param,
                    IntfPIType::Edge,
                    pivot.min(pedge),
                    pivot.max(pedge),
                    0.0,
                    1.0,
                );
                self.base.my_s_poins.push(sp);
            } else if d_pi_e[0] > 0.0 && d_pi_e[1] > 0.0 && d_pi_e[2] > 0.0 {
                let sp = IntfSectionPoint::with_3d(
                    &GpPnt::from_xyz(&sp_lieu),
                    typ_on_g,
                    0,
                    self.i_lin,
                    param,
                    IntfPIType::Face,
                    t_tri,
                    0,
                    0.0,
                    1.0,
                );
                self.base.my_s_poins.push(sp);
            } else {
                // Modified by Sergey KHROMOV - Fri Dec 7 14:40:11 2001
                for i in 1..=3usize {
                    let ind_p1 = if i == 3 { p_tri[0] } else { p_tri[i] };
                    let ind_p2 = p_tri[i - 1];

                    if TH::is_on_bound(the_polyh, ind_p1, ind_p2) {
                        let deflection = TH::get_border_deflection(the_polyh);
                        let beg_p = TH::point(the_polyh, ind_p1);
                        let end_p = TH::point(the_polyh, ind_p2);
                        let vec_tri = GpVec::from_pnts(&beg_p, &end_p);
                        // gp_Dir DirTri(VecTri) raises on a null vector in OCCT;
                        // the port skips the degenerate boundary edge.
                        if let Ok(dir_tri) = GpDir::from_vec(&vec_tri) {
                            let lin_tri = GpLin::from_pnt_dir(beg_p, dir_tri);
                            let a_p_on_e = GpPnt::from_xyz(&sp_lieu);
                            let a_dist = lin_tri.distance(&a_p_on_e);

                            if a_dist <= deflection {
                                let a_v_loc_p_on_e = GpVec::from_pnts(&beg_p, &a_p_on_e);
                                let a_vec_dir_tri =
                                    GpVec::new(dir_tri.x(), dir_tri.y(), dir_tri.z());
                                let a_par = a_v_loc_p_on_e.dot(&a_vec_dir_tri);
                                let a_max_par = vec_tri.magnitude();

                                if a_par >= 0.0 && a_par <= a_max_par {
                                    let sp = IntfSectionPoint::with_3d(
                                        &GpPnt::from_xyz(&sp_lieu),
                                        typ_on_g,
                                        0,
                                        self.i_lin,
                                        param,
                                        IntfPIType::Face,
                                        t_tri,
                                        0,
                                        0.0,
                                        1.0,
                                    );
                                    self.base.my_s_poins.push(sp);
                                }
                            }
                        }
                    }
                }
            }
        } // if NoIntersectionWithTriangle == false

        // UNPORTED: `gxx:1296-1367` - Extrema_ExtElC proximity block
        // (see the module header).
    }
}
