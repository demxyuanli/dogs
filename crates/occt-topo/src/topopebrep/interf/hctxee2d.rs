//! `TopOpeBRep_Hctxee2d` — the edge/edge intersection context: the 2D images of
//! two edges on one of the two faces, with their `IntRes2d` domains.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_Hctxee2d.hxx:27-58`,
//! `.cxx:44-251`.
//!
//! Translation notes:
//! - `Geom2dAdaptor_Curve` has no port class: the port's 2D curves carry their
//!   own parameter range (`Arc<dyn Curve2d>`), so `myCurve1/2` store the curve
//!   directly.
//! - `SetEdges`'s `FC2D_CurveOnSurface(E, F, first, last, tolpc)` is
//!   [`crate::boptools_2d::curve_on_surface_range`], the port's
//!   `BRep_Tool::CurveOnSurface(E, F, first, last)` plus the stored
//!   representation's range (`TopOpeBRepTool_2d.cxx:356-376`).
//! - `TopoDS_Shape& Edge(Index)` returns a `TopoDS_Edge` member upcast in OCCT;
//!   here the accessor returns `&Edge`.
//! - `FUN_tool_projPonF` / `FUN_tool_onapex` are
//!   [`crate::topopebrep::tool::fun_tool_proj_pon_f`] /
//!   [`crate::topopebrep::tool::fun_tool_onapex`].
//! - The `BRep_Tool::Surface(F, L)` identity test `S1 == S2 && L1 == L2`
//!   (`cxx:86-92`) is the registered surface's `Arc` identity plus the two
//!   faces' `LocationKey`s.
//! - `TopOpeBRepTool_CurveTool::MakePCurveOnFace(F1, NC, tolreached2d)`
//!   (`CurveTool.cxx:1040-1136`) is **UNPORTED**: its body is
//!   `ProjLib_ProjectedCurve` + `::MakePCurve` plus a sphere quasi-periodic shift
//!   (`cxx:1084-1117`) and a UV clip (`cxx:1119-1136`), and neither
//!   `ProjLib_ProjectedCurve` nor `BRepApprox_Approx` exists in this port yet.
//!   `project_edge_on_face` therefore delegates the projection to
//!   [`crate::pcurve::make_pcurve_on_face`] (the port's
//!   `BOPTools_AlgoTools2D::MakePCurveOnFace` projection of the edge's 3D curve
//!   onto the face) and carries neither `TolReached2d` nor the sphere shift.

use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_core::intres2d::IntRes2dDomain;
use occt_geom2d::curve::Curve2d;

use crate::abs::ShapeType;
use crate::brep_tool::BRepTool;
use crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::shape::{Edge, Face};
use crate::topo_tools::{edge_vertices, is_same, map_shapes_and_ancestors, shape_map_key};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;
use crate::topopebrep::tool::fun_tool_onapex;
use crate::topopebrep::tool::fun_tool_proj_pon_f;

/// `class TopOpeBRep_Hctxee2d` (`hxx:31-58`).
#[derive(Clone)]
pub struct Hctxee2d {
    /// `myEdge1` (`hxx:51`)
    edge1: Option<Edge>,
    /// `myCurve1` (`hxx:52`)
    curve1: Option<Arc<dyn Curve2d>>,
    /// `myDomain1` (`hxx:53`)
    domain1: IntRes2dDomain,
    /// `myEdge2` (`hxx:54`)
    edge2: Option<Edge>,
    /// `myCurve2` (`hxx:55`)
    curve2: Option<Arc<dyn Curve2d>>,
    /// `myDomain2` (`hxx:56`)
    domain2: IntRes2dDomain,
}

impl Default for Hctxee2d {
    /// `TopOpeBRep_Hctxee2d::TopOpeBRep_Hctxee2d() = default` (`cxx:44`).
    fn default() -> Self {
        Self::new()
    }
}

impl Hctxee2d {
    /// `TopOpeBRep_Hctxee2d()` (`cxx:44`).
    pub fn new() -> Self {
        Self {
            edge1: None,
            curve1: None,
            domain1: IntRes2dDomain::new(),
            edge2: None,
            curve2: None,
            domain2: IntRes2dDomain::new(),
        }
    }

    /// `SetEdges(E1, E2, BAS1, BAS2)` (`cxx:47-197`).
    pub fn set_edges(
        &mut self,
        e1: &Edge,
        e2: &Edge,
        bas1: &SurfaceAdaptor,
        bas2: &SurfaceAdaptor,
    ) -> Result<(), String> {
        let f1 = bas1
            .face()
            .cloned()
            .ok_or_else(|| "TopOpeBRep_Hctxee2d::SetEdges : no face 1".to_string())?;
        let st1 = bas1.get_type();
        let f2 = bas2
            .face()
            .cloned()
            .ok_or_else(|| "TopOpeBRep_Hctxee2d::SetEdges : no face 2".to_string())?;

        self.edge1 = Some(e1.clone());
        self.edge2 = Some(e2.clone());

        // `cxx:62-78`: PC1 = FC2D_CurveOnSurface(myEdge1, F1, first, last,
        // tolpc); a null curve raises (`cxx:63-66`).
        let (pc1, first, last) = crate::boptools_2d::curve_on_surface_range(e1, &f1)
            .ok_or_else(|| "TopOpeBRep_Hctxee2d::SetEdges : no 2d curve".to_string())?;
        self.curve1 = Some(pc1);
        let (pfirst, plast) = BRepTool::uv_points(e1, &f1);
        let tole = BRepTool::edge_tolerance(e1);
        self.domain1 = IntRes2dDomain::bounded(&pfirst, first, tole, &plast, last, tole);

        let memesfaces = is_same(&f1.0, &f2.0);
        let memesupport = same_support(&f1, &f2);

        // `first` / `last` are by-reference in OCCT and the arms below overwrite
        // them; the `myDomain2` construction at `cxx:176-179` reads the values
        // the chosen arm left behind.
        let mut first = first;
        let mut last = last;

        if st1 == GeomAbsSurfaceType::Plane || memesfaces || memesupport {
            // `cxx:93-111`: edge 2 read on F1 directly.
            let Some((pc2, first2, last2)) = crate::boptools_2d::curve_on_surface_range(e2, &f1)
            else {
                return Err("TopOpeBRep_Hctxee2d::SetEdges : no 2d curve".to_string());
            };
            self.curve2 = Some(pc2);
            let (pfirst, plast) = BRepTool::uv_points(e2, &f1);
            let tole = BRepTool::edge_tolerance(e2);
            self.domain2 = IntRes2dDomain::bounded(&pfirst, first2, tole, &plast, last2, tole);
            return Ok(());
        }

        // `cxx:112-175`.
        let pc2on1 = project_edge_on_face(e1, e2, &f1, &mut first, &mut last)?;

        // `cxx:174-179`: `if (!PC2on1.IsNull())`.
        let Some(pc2on1) = pc2on1 else {
            return Ok(());
        };
        self.curve2 = Some(pc2on1.clone());
        let tole = BRepTool::edge_tolerance(e2);
        let pfirst = pc2on1.d0(first);
        let plast = pc2on1.d0(last);
        self.domain2 = IntRes2dDomain::bounded(&pfirst, first, tole, &plast, last, tole);
        Ok(())
    }

    /// `Edge(Index)` (`cxx:202-217`). OCCT returns the `TopoDS_Shape` base of
    /// `myEdge1/2`; a not-yet-set edge is an error here.
    pub fn edge(&self, index: i32) -> Result<&Edge, String> {
        match index {
            1 => self
                .edge1
                .as_ref()
                .ok_or_else(|| "TopOpeBRep_Hctxee2d::Edge".to_string()),
            2 => self
                .edge2
                .as_ref()
                .ok_or_else(|| "TopOpeBRep_Hctxee2d::Edge".to_string()),
            _ => Err("TopOpeBRep_Hctxee2d::Edge".to_string()),
        }
    }

    /// `Curve(Index)` (`cxx:219-236`).
    pub fn curve(&self, index: i32) -> Result<&Arc<dyn Curve2d>, String> {
        match index {
            1 => self
                .curve1
                .as_ref()
                .ok_or_else(|| "TopOpeBRep_Hctxee2d::Curve".to_string()),
            2 => self
                .curve2
                .as_ref()
                .ok_or_else(|| "TopOpeBRep_Hctxee2d::Curve".to_string()),
            _ => Err("TopOpeBRep_Hctxee2d::Curve".to_string()),
        }
    }

    /// `Domain(Index)` (`cxx:238-253`).
    pub fn domain(&self, index: i32) -> Result<&IntRes2dDomain, String> {
        match index {
            1 => Ok(&self.domain1),
            2 => Ok(&self.domain2),
            _ => Err("TopOpeBRep_Hctxee2d::Domain".to_string()),
        }
    }
}

/// `S1 == S2 && L1 == L2` (`cxx:86-92`): the two faces share one surface handle
/// and one location.
fn same_support(f1: &Face, f2: &Face) -> bool {
    let (Some(s1), Some(s2)) = (BRepTool::face_surface(f1), BRepTool::face_surface(f2)) else {
        return false;
    };
    Arc::ptr_eq(&s1, &s2) && f1.0.location().key() == f2.0.location().key()
}

/// The `else` arm of `SetEdges` (`cxx:112-175`): the p-curve of `e2` on `f1`
/// (OCCT's `PC2on1`), or `None` for the `if (!ok) return;` of `cxx:128-131`
/// and for the `NYIxpu210998` empty arm (`cxx:159-162`).
///
/// `first` / `last` are updated exactly where OCCT's by-reference parameters
/// are: `BRep_Tool::CurveOnSurface(ee, F1, first, last)` (`cxx:155`) and
/// `BRep_Tool::Curve(myEdge2, loc, first, last)` (`cxx:168`).
pub(crate) fn project_edge_on_face(
    e1: &Edge,
    e2: &Edge,
    f1: &Face,
    first: &mut f64,
    last: &mut f64,
) -> Result<Option<Arc<dyn Curve2d>>, String> {
    let Some(surf) = BRepTool::face_surface(f1) else {
        return Err("TopOpeBRep_Hctxee2d::SetEdges : no surface".to_string());
    };

    if BRepTool::is_degenerated(e2) {
        // xpu210998 : cto900Q3 (`cxx:119-162`).
        let (Some(v2), _) = edge_vertices(e2) else {
            return Ok(None);
        };
        let pt2 = BRepTool::vertex_point(&v2);
        let Some((u, v, _d)) = fun_tool_proj_pon_f(&pt2, f1) else {
            // `if (!ok) return;` (`cxx:128-131`, `nyiRaise`).
            return Ok(None);
        };
        let uv2 = GpPnt2d::new(u, v);
        if fun_tool_onapex(&uv2, surf.as_ref()) {
            // `cxx:135-157`: the apex arm walks the ancestors of one vertex of
            // edge 1 (`TopExp::Vertices(myEdge1, vf, vl)`).
            let (vf, vl) = edge_vertices(e1);
            let (Some(vf), Some(vl)) = (vf, vl) else {
                return Ok(None);
            };
            let ptf = BRepTool::vertex_point(&vf);
            let df = pt2.distance(&ptf);
            let tolf = BRepTool::vertex_tolerance(&vf);
            let onf = df < tolf;
            let v1 = if onf { vf } else { vl };
            let map = map_shapes_and_ancestors(&f1.0, ShapeType::Vertex, ShapeType::Edge);
            if let Some(ancestors) = map.get(&shape_map_key(&v1.0)) {
                for anc in ancestors {
                    let Some(ee) = Edge::wrap(anc.clone()) else {
                        continue;
                    };
                    if !BRepTool::is_degenerated(&ee) {
                        continue;
                    }
                    // `PC2on1 = BRep_Tool::CurveOnSurface(ee, F1, first, last)`
                    // (`cxx:155`): also rewrites the two parameters.
                    if let Some((pc, t1, t2)) =
                        crate::boptools_2d::curve_on_surface_range(&ee, f1)
                    {
                        *first = t1;
                        *last = t2;
                        return Ok(Some(pc));
                    }
                }
            }
        }
        // `NYIxpu210998` (`cxx:159-162`): nothing is built.
        return Ok(None);
    }

    // `cxx:165-172`: project the 3D curve of edge 2 onto the surface of face 1.
    // `BRep_Tool::Curve(myEdge2, loc, first, last)` rewrites the parameters.
    let (f, l) = BRepTool::edge_parameters(e2);
    *first = f;
    *last = l;
    match crate::pcurve::make_pcurve_on_face(e2, f1) {
        Ok(pc) => Ok(Some(pc)),
        Err(_) => Ok(None),
    }
}
