//! `TopOpeBRep_EdgesIntersector` — the intersection of two edges that lie on
//! the same surface.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_EdgesIntersector.hxx:45-204`,
//! `TopOpeBRep_EdgesIntersector.cxx:70-1199`,
//! `TopOpeBRep_EdgesIntersector_1.cxx:32-463`.
//!
//! Translation notes:
//! - `Geom2dAdaptor_Curve` has no port class (see
//!   [`crate::topopebrep::interf::hctxee2d`]): `myCurve1/2` store the
//!   `Arc<dyn Curve2d>` together with its own range, so `Curve(Index)` returns
//!   that curve and `GetType` / `Circle` / `Resolution` read it directly.
//! - `FC2D_CurveOnSurface` and `FC2D_MakeCurveOnSurface`
//!   (`TopOpeBRepTool_2d.cxx:342-376`) both end in
//!   `FC2D_make2d` (`cxx:263-338`), whose first step is
//!   `BRep_Tool::CurveOnSurface` and whose fallback is
//!   `TopOpeBRepTool_CurveTool::MakePCurveOnFace`; [`super::interf::hctxee2d`]'s
//!   `project_edge_on_face` is the translation of the same control flow
//!   (`Hctxee2d.cxx:47-197` is a copy of `cxx:316-495`), so both OCCT entry
//!   points map onto `curve_on_surface_range` + that helper.
//! - `mySurface1/2` are `BRepAdaptor_Surface` handles; [`SurfaceAdaptor`] is
//!   the port's adaptor, and `EvalD0(u, v)` is its surface's `d0`.
//! - `mysp2d` / `mylpnt` / `mylseg` are `NCollection_Sequence` (1-based);
//!   they are `Vec`s here and every index is translated explicitly.
//! - `TopOpeBRep_EdgesIntersector::Dump` (`cxx:1164-1199`) is `OCCT_DEBUG`
//!   only, so it is not ported.

mod perform;
mod point_ops;

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_core::intres2d::{IntRes2dDomain, IntRes2dIntersectionPoint, IntRes2dIntersectionSegment};
use occt_geom2d::curve::Curve2d;
use occt_geom2d::geom2d_int::Geom2dIntGInter;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::is_same;
use crate::topopebrep::interf::point2d::Point2d;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;
use crate::topopebrep::tool::shape_tool;
use crate::topopebrep::tool::tol::ftol_face_tolerances2d;

/// `class TopOpeBRep_EdgesIntersector` (`hxx:45-204`).
pub struct EdgesIntersector {
    /// `myFace1` (`hxx:167`)
    face1: Face,
    /// `myFace2` (`hxx:168`)
    face2: Face,
    /// `mySurface1` (`hxx:169`)
    surface1: SurfaceAdaptor,
    /// `mySurface2` (`hxx:170`)
    surface2: SurfaceAdaptor,
    /// `mySurfaceType1` (`hxx:171`)
    surface_type1: GeomAbsSurfaceType,
    /// `mySurfaceType2` (`hxx:172`)
    surface_type2: GeomAbsSurfaceType,
    /// `mySurfacesSameOriented` (`hxx:173`)
    surfaces_same_oriented: bool,
    /// `myFacesSameOriented` (`hxx:174`)
    faces_same_oriented: bool,
    /// `myDomain1` (`hxx:175`)
    domain1: IntRes2dDomain,
    /// `myDomain2` (`hxx:176`)
    domain2: IntRes2dDomain,
    /// `myEdge1` (`hxx:177`)
    edge1: Edge,
    /// `myEdge2` (`hxx:178`)
    edge2: Edge,
    /// `myCurve1` (`hxx:179`)
    curve1: Option<Arc<dyn Curve2d>>,
    /// `myCurve2` (`hxx:180`)
    curve2: Option<Arc<dyn Curve2d>>,
    /// `myTol1` (`hxx:181`)
    tol1: f64,
    /// `myTol2` (`hxx:182`)
    tol2: f64,
    /// `myTolForced` (`hxx:183`)
    tol_forced: bool,
    /// `myIntersector` (`hxx:185`)
    intersector: Geom2dIntGInter,
    /// `mylpnt` (`hxx:186`)
    lpnt: Vec<IntRes2dIntersectionPoint>,
    /// `mylseg` (`hxx:187`)
    lseg: Vec<IntRes2dIntersectionSegment>,
    /// `myNbPoints` (`hxx:188`)
    nb_points: i32,
    /// `myNbSegments` (`hxx:189`)
    nb_segments: i32,
    /// `myTrueNbPoints` (`hxx:190`)
    true_nb_points: i32,
    /// `myPointIndex` (`hxx:191`), 1-based as in OCCT.
    point_index: i32,
    /// `myIsVertexPointIndex` (`hxx:192`)
    is_vertex_point_index: i32,
    /// `myIsVertexIndex` (`hxx:193`)
    is_vertex_index: i32,
    /// `myIsVertexValue` (`hxx:194`)
    is_vertex_value: bool,
    /// `myIsVertexVertex` (`hxx:195`)
    is_vertex_vertex: Option<crate::shape::Vertex>,
    /// `myDimension` (`hxx:196`)
    dimension: i32,
    /// `myHasSegment` (`hxx:197`)
    has_segment: bool,
    /// `mySameDomain` (`hxx:198`)
    same_domain: bool,
    /// `myf1surf1F_sameoriented` (`hxx:199`)
    f1surf1f_sameoriented: bool,
    /// `myf2surf1F_sameoriented` (`hxx:200`)
    f2surf1f_sameoriented: bool,
    /// `mysp2d` (`hxx:201`)
    sp2d: Vec<Point2d>,
    /// `myip2d` (`hxx:202`), 1-based as in OCCT.
    ip2d: i32,
    /// `mynp2d` (`hxx:203`)
    np2d: i32,
    /// `myselectkeep` (`hxx:204`)
    select_keep: bool,
}

impl EdgesIntersector {
    /// `TopOpeBRep_EdgesIntersector()` (`cxx:70-94`).
    pub fn new() -> Self {
        let mut me = Self {
            // `mySurface1 = new BRepAdaptor_Surface()` (`cxx:72-73`).
            face1: Face::new(),
            face2: Face::new(),
            surface1: SurfaceAdaptor::of_face(&Face::new()),
            surface2: SurfaceAdaptor::of_face(&Face::new()),
            surface_type1: GeomAbsSurfaceType::OtherSurface,
            surface_type2: GeomAbsSurfaceType::OtherSurface,
            surfaces_same_oriented: false,
            faces_same_oriented: false,
            domain1: IntRes2dDomain::new(),
            domain2: IntRes2dDomain::new(),
            edge1: Edge::new(),
            edge2: Edge::new(),
            curve1: None,
            curve2: None,
            // `myTol1 = 0.; myTol2 = 0.` (`cxx:77-78`).
            tol1: 0.0,
            tol2: 0.0,
            tol_forced: false,
            intersector: Geom2dIntGInter::new(),
            lpnt: Vec::new(),
            lseg: Vec::new(),
            nb_points: 0,
            // `myNbSegments = 0; myHasSegment = false` (`cxx:83-84`).
            nb_segments: 0,
            true_nb_points: 0,
            point_index: 0,
            is_vertex_point_index: 0,
            is_vertex_index: 0,
            is_vertex_value: false,
            is_vertex_vertex: None,
            // `myDimension = 2` (`cxx:79`).
            dimension: 2,
            has_segment: false,
            same_domain: false,
            // `myf1surf1F_sameoriented = true; myf2surf1F_sameoriented = true`
            // (`cxx:80-81`).
            f1surf1f_sameoriented: true,
            f2surf1f_sameoriented: true,
            sp2d: Vec::new(),
            // `myip2d = mynp2d = 0` (`cxx:90`).
            ip2d: 0,
            np2d: 0,
            // `myselectkeep = true` (`cxx:91`).
            select_keep: true,
        };
        // `SetSameDomain(false)` (`cxx:85`).
        me.set_same_domain(false);
        me
    }

    /// `SetFaces(F1, F2)` (`cxx:98-103`).
    pub fn set_faces(&mut self, f1: &TopoShape, f2: &TopoShape) {
        let (b1, b2) = (occt_core::bnd::BndBox::new(), occt_core::bnd::BndBox::new());
        self.set_faces_with_boxes(f1, f2, &b1, &b2);
    }

    /// `SetFaces(F1, F2, B1, B2)` (`cxx:106-173`).
    pub fn set_faces_with_boxes(
        &mut self,
        f1: &TopoShape,
        f2: &TopoShape,
        b1: &occt_core::bnd::BndBox,
        b2: &occt_core::bnd::BndBox,
    ) {
        // `bool computerestriction = false` (`cxx:110`).
        self.f1surf1f_sameoriented = true;
        self.f2surf1f_sameoriented = true;
        self.surfaces_same_oriented = true;
        self.faces_same_oriented = true;

        // `myFace1 = TopoDS::Face(F1); S1.Initialize(myFace1, false)` (`cxx:120-122`).
        self.face1 = Face::wrap(f1.clone()).unwrap_or_else(Face::new);
        self.surface1 = SurfaceAdaptor::of_face(&self.face1);
        self.surface_type1 = self.surface1.get_type();

        self.face2 = Face::wrap(f2.clone()).unwrap_or_else(Face::new);
        self.surface2 = SurfaceAdaptor::of_face(&self.face2);
        self.surface_type2 = self.surface2.get_type();

        // `TopoDS_Face face1forward = myFace1; face1forward.Orientation(FORWARD)`
        // (`cxx:130-132`).
        let mut face1forward = self.face1.clone();
        face1forward.0.set_orientation(Orientation::Forward);

        // `so11 = FacesSameOriented(face1forward, myFace1)` (`cxx:134-135`).
        self.f1surf1f_sameoriented = shape_tool::faces_same_oriented(
            &face1forward.0,
            &self.face1.0,
        );
        // `so21 = FacesSameOriented(face1forward, myFace2)` (`cxx:137-138`).
        self.f2surf1f_sameoriented =
            shape_tool::faces_same_oriented(&face1forward.0, &self.face2.0);

        // `mySurfacesSameOriented = SurfacesSameOriented(S1, S2)` (`cxx:139`).
        self.surfaces_same_oriented =
            shape_tool::surfaces_same_oriented(&self.surface1, &self.surface2);
        // `myFacesSameOriented = FacesSameOriented(myFace1, myFace2)` (`cxx:140`).
        self.faces_same_oriented =
            shape_tool::faces_same_oriented(&self.face1.0, &self.face2.0);

        if !self.tol_forced {
            // `FTOL_FaceTolerances2d(B1, B2, myFace1, myFace2, S1, S2, myTol1,
            // myTol2)` (`cxx:143-144`).
            let (mut t1, mut t2) = (self.tol1, self.tol2);
            ftol_face_tolerances2d(
                b1,
                b2,
                &self.face1,
                &self.face2,
                &self.surface1,
                &self.surface2,
                &mut t1,
                &mut t2,
            );
            // `myTol1 = (myTol1 > 1.e-4) ? 1.e-4 : myTol1` (`cxx:145-146`).
            self.tol1 = if t1 > 1.0e-4 { 1.0e-4 } else { t1 };
            self.tol2 = if t2 > 1.0e-4 { 1.0e-4 } else { t2 };
        }
    }

    /// `ForceTolerances(Tol1, Tol2)` (`cxx:177-182`).
    pub fn force_tolerances(&mut self, tol1: f64, tol2: f64) {
        self.tol1 = tol1;
        self.tol2 = tol2;
        self.tol_forced = true;
    }

    /// `Dimension(D)` (`cxx:691-697`).
    pub fn set_dimension(&mut self, dim: i32) {
        if dim == 1 || dim == 2 {
            self.dimension = dim;
        }
    }

    /// `Dimension()` (`cxx:701-704`).
    pub fn dimension(&self) -> i32 {
        self.dimension
    }

    /// `IsEmpty()` (`cxx:951-955`).
    pub fn is_empty(&self) -> bool {
        self.np2d == 0
    }

    /// `HasSegment()` (`cxx:958-961`).
    pub fn has_segment(&self) -> bool {
        self.has_segment
    }

    /// `SameDomain()` (`cxx:965-968`).
    pub fn same_domain(&self) -> bool {
        self.same_domain
    }

    /// `Edge(Index)` (`cxx:972-986`).
    pub fn edge(&self, index: i32) -> Result<&TopoShape, String> {
        match index {
            1 => Ok(&self.edge1.0),
            2 => Ok(&self.edge2.0),
            _ => Err("TopOpeBRep_EdgesIntersector::Edge".to_string()),
        }
    }

    /// `Curve(Index)` (`cxx:990-1004`).
    pub fn curve(&self, index: i32) -> Result<&Arc<dyn Curve2d>, String> {
        let c = match index {
            1 => self.curve1.as_ref(),
            2 => self.curve2.as_ref(),
            _ => return Err("TopOpeBRep_EdgesIntersector::Curve".to_string()),
        };
        c.ok_or_else(|| "TopOpeBRep_EdgesIntersector::Curve".to_string())
    }

    /// `Face(Index)` (`cxx:1008-1022`).
    pub fn face(&self, index: i32) -> Result<&Face, String> {
        match index {
            1 => Ok(&self.face1),
            2 => Ok(&self.face2),
            _ => Err("TopOpeBRep_EdgesIntersector::Face".to_string()),
        }
    }

    /// `Surface(Index)` (`cxx:1026-1040`).
    pub fn surface(&self, index: i32) -> Result<&SurfaceAdaptor, String> {
        match index {
            1 => Ok(&self.surface1),
            2 => Ok(&self.surface2),
            _ => Err("TopOpeBRep_EdgesIntersector::Surface".to_string()),
        }
    }

    /// `SurfacesSameOriented()` (`cxx:1044-1047`).
    pub fn surfaces_same_oriented(&self) -> bool {
        self.surfaces_same_oriented
    }

    /// `FacesSameOriented()` (`cxx:1051-1054`).
    pub fn faces_same_oriented(&self) -> bool {
        self.faces_same_oriented
    }

    /// `ToleranceMax()` (`cxx:1134-1138`).
    pub fn tolerance_max(&self) -> f64 {
        self.tol1.max(self.tol2)
    }

    /// `Tolerances(tol1, tol2)` (`cxx:1142-1146`).
    pub fn tolerances(&self) -> (f64, f64) {
        (self.tol1, self.tol2)
    }

    /// `NbPoints()` (`cxx:1150-1153`).
    pub fn nb_points(&self) -> i32 {
        self.nb_points
    }

    /// `NbSegments()` (`cxx:1157-1160`).
    pub fn nb_segments(&self) -> i32 {
        self.nb_segments
    }

    /// `InitPoint(selectkeep)` (`cxx:1058-1063`).
    pub fn init_point(&mut self, select_keep: bool) {
        self.select_keep = select_keep;
        self.ip2d = 1;
        self.np2d = self.sp2d.len() as i32;
        self.find();
    }

    /// `MorePoint()` (`cxx:1068-1071`).
    pub fn more_point(&self) -> bool {
        self.ip2d <= self.np2d
    }

    /// `NextPoint()` (`cxx:1076-1079`).
    pub fn next_point(&mut self) {
        self.ip2d += 1;
        self.find();
    }

    /// `Find()` (`cxx:1084-1105`).
    fn find(&mut self) {
        while self.ip2d <= self.np2d {
            if self.select_keep {
                let kf = self
                    .sp2d
                    .get((self.ip2d - 1) as usize)
                    .map(|p| p.keep())
                    .unwrap_or(false);
                if kf {
                    break;
                } else {
                    self.ip2d += 1;
                }
            } else {
                break;
            }
        }
    }

    /// `Points()` (`cxx:1109-1112`).
    pub fn points(&self) -> &[Point2d] {
        &self.sp2d
    }

    /// `Point()` (`cxx:1116-1119`).
    pub fn point(&self) -> Result<&Point2d, String> {
        self.sp2d
            .get((self.ip2d - 1) as usize)
            .ok_or_else(|| "TopOpeBRep_EdgesIntersector::Point".to_string())
    }

    /// `Point(I)` (`cxx:1123-1130`): raises when `I` is out of `1..=Length`.
    pub fn point_at(&self, i: i32) -> Result<&Point2d, String> {
        if i < 1 || i > self.sp2d.len() as i32 {
            return Err("TopOpeBRep_EdgesIntersector::Point(I)".to_string());
        }
        Ok(&self.sp2d[(i - 1) as usize])
    }

    /// The mutable form of `Point()` used by the `Perform` tail
    /// (`cxx:694-691` region: `mysp2d(myip2d)`).
    fn point_mut(&mut self) -> &mut Point2d {
        &mut self.sp2d[(self.ip2d - 1) as usize]
    }
}

impl Default for EdgesIntersector {
    fn default() -> Self {
        Self::new()
    }
}

/// Whether two faces share one surface handle and one location
/// (`S1 == S2 && L1 == L2`, `cxx:369-375`).
pub(super) fn same_support(f1: &Face, f2: &Face) -> bool {
    let (Some(s1), Some(s2)) = (BRepTool::face_surface(f1), BRepTool::face_surface(f2)) else {
        return false;
    };
    Arc::ptr_eq(&s1, &s2) && f1.0.location().key() == f2.0.location().key()
}

/// `myFace1.IsSame(myFace2)` (`cxx:367`).
pub(super) fn same_face(f1: &Face, f2: &Face) -> bool {
    is_same(&f1.0, &f2.0)
}

/// `BRep_Tool::UVPoints(E, F, pfirst, plast)` (`BRep_Tool.cxx`), typed as the
/// two `gp_Pnt2d` the domains are built from.
pub(super) fn uv_points(e: &Edge, f: &Face) -> (GpPnt2d, GpPnt2d) {
    BRepTool::uv_points(e, f)
}
