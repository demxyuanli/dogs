//! `BRepFill_Sweep::BuildFace` (`BRepFill_Sweep.cxx:526-762`): one natural
//! border face of the swept surface from its four boundary edges.
//!
//! Translation notes:
//! - The two candidate planar surfaces are `GeomLib_IsPlanarSurface` over `S`
//!   (`cxx:564-569`) and `BRepLib_FindSurface(theWire, -1, true)`
//!   (`cxx:582-591`), both held as `Arc<dyn Surface>` over a `Geom_Plane`.
//! - `BRepBuilderAPI_MakeWire::Add` is [`crate::brep_lib_make_wire::MakeWire`];
//!   the "last edge added to WW" of `cxx:616-625` is the wire's last child, and
//!   `IsSame` (TShape identity) decides whether it is a new edge.
//! - `BRepLib_MakeFace(thePlane, WW)`'s plane+wire constructor
//!   (`BRepLib_MakeFace.cxx:262-271`) is [`make_plane_face`]:
//!   `Init(Pl, Bound=false, Confusion)` (`cxx:366-385`) then `Add(W)`
//!   (`cxx:873-880`), then `CheckInside()` when `Inside && IsClosed(W)`.
//! - `BRep_TEdge::Tolerance(Precision::Confusion())` of `cxx:573-578` is
//!   [`GeometryRegistry::set_edge_tolerance`]; the restore of `cxx:590-593` is
//!   the same setter with the saved values.
//! - `BRep_Builder::UpdateEdge(E, null C2d, S, L, Tol)` (`cxx:721-724`) is
//!   [`GeometryRegistry::update_edge_pcurve_on_surface`] with `None`.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::GpVec;
use occt_core::precision::CONFUSION;
use occt_geom::geom_lib_is_planar_surface::IsPlanarSurface;
use occt_geom::plane::GeomPlane;
use occt_geom::surface::Surface;

use crate::brep_lib_find_surface::FindSurface;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tool::BRepTool;
use crate::brep_top_adaptor_fclass2d::BrepTopAdaptorFClass2d;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{is_same, shape_map_key, ShapeMapKey};

use super::helpers;
use super::same_parameter::update_edge_on_plane;

/// `BuildFace` (`cxx:526-762`).
///
/// `None` is the reference's null `F`: the caller passed an empty `TopoDS_Face`
/// and `BRepLib_MakeFace` errored (`cxx:710-718` leaves it untouched).
pub fn build_face(
    s: &Arc<dyn Surface>,
    e1: &Edge,
    e2: &Edge,
    e3: &Edge,
    e4: &Edge,
    ee_map: &mut HashMap<ShapeMapKey, TopoShape>,
    exch_uv: bool,
    u_reverse: bool,
) -> Option<Face> {
    let reg = GeometryRegistry::global();

    // `cxx:550-557`.
    let tol1 = BRepTool::edge_tolerance(e1);
    let tol2 = BRepTool::edge_tolerance(e2);
    let tol3 = BRepTool::edge_tolerance(e3);
    let tol4 = BRepTool::edge_tolerance(e4);
    let tol = tol1.min(tol2).min(tol3).min(tol4);

    let mut is_plan = false;
    // `thePlane`: the `Geom_Plane` the face is built on.
    let mut the_plane: Option<Arc<dyn Surface>> = None;

    // `cxx:562-595`.
    if !is_same(&e1.0, &e3.0) && !is_same(&e2.0, &e4.0) {
        let is_p = IsPlanarSurface::new(s, tol);
        if is_p.is_planar() {
            is_plan = true;
            if let Ok(pln) = is_p.plan() {
                the_plane = Some(Arc::new(GeomPlane::new(pln)));
            }
        } else {
            for e in [e1, e2, e3, e4] {
                reg.set_edge_tolerance(&e.0, CONFUSION);
            }

            // `BRepLib_MakeWire(E1, E2, E3, E4)` (`cxx:581`).
            let mut mw = MakeWire::new();
            for e in [e1, e2, e3, e4] {
                let _ = mw.add(e);
            }
            let the_wire = mw.wire();
            if helpers::number_of_poles(&the_wire) <= 100 {
                // `BRepLib_FindSurface FS(theWire, -1, true)` (`cxx:584`).
                let fs = FindSurface::new(&the_wire.0, -1.0, true, false);
                if fs.is_found() {
                    if let Some(surf) = fs.surface() {
                        // `down_cast<Geom_Plane>` (`cxx:588`).
                        if surf.gp_pln().is_some() {
                            is_plan = true;
                            the_plane = Some(surf);
                        }
                    }
                }
            }
            reg.set_edge_tolerance(&e1.0, tol1);
            reg.set_edge_tolerance(&e2.0, tol2);
            reg.set_edge_tolerance(&e3.0, tol3);
            reg.set_edge_tolerance(&e4.0, tol4);
        }
    }

    // `cxx:598-688`: the wire, edge by edge.
    let mut mw = MakeWire::new();

    let mut le1 = e1.clone();
    helpers::oriente(s, &mut le1);
    if !is_plan || !BRepTool::is_degenerated(&le1) {
        let _ = mw.add(&le1);
    }

    let mut le2 = e2.clone();
    helpers::oriente(s, &mut le2);
    if !is_plan || !BRepTool::is_degenerated(&le2) {
        add_and_record(&mut mw, &le2, ee_map);
    }

    let mut e = if is_same(&e3.0, &e1.0) {
        // `cxx:630-634`: `E = e1; E.Reverse();`
        let mut x = le1.clone();
        x.0.reverse();
        x
    } else {
        let mut x = e3.clone();
        helpers::oriente(s, &mut x);
        x
    };
    if !is_plan || !BRepTool::is_degenerated(&e) {
        add_and_record(&mut mw, &e, ee_map);
    }

    e = if is_same(&e4.0, &e2.0) {
        // `cxx:660-664`: `E = e2; E.Reverse();`
        let mut x = le2.clone();
        x.0.reverse();
        x
    } else {
        let mut x = e4.clone();
        helpers::oriente(s, &mut x);
        x
    };
    if !is_plan || !BRepTool::is_degenerated(&e) {
        add_and_record(&mut mw, &e, ee_map);
    }

    // `cxx:690`: `WW = B.Wire();`
    let ww = mw.wire();

    // `cxx:693-750`: the face.
    let mut f: Face;
    if is_plan {
        // `cxx:698-708`: orient the plane so that its normal follows
        // `S->D1` at the parametric centre.
        let mut plane = the_plane.clone();
        let (u_first, u_last) = s.u_range();
        let (v_first, v_last) = s.v_range();
        let (_, du, dv) = s.d1((u_first + u_last) / 2.0, (v_first + v_last) / 2.0);
        let ns = du.crossed(&dv);
        let reversed = match &plane {
            Some(p) => match p.gp_pln() {
                Some(pln) => {
                    let ax = pln.axis();
                    let np = ax.direction();
                    let np = GpVec::new(np.x(), np.y(), np.z());
                    ns.dot(&np) < 0.0
                }
                None => false,
            },
            None => false,
        };
        if reversed {
            // `thePlane->UReverse()` (`cxx:706`).
            plane = plane.as_ref().and_then(|p| p.u_reversed());
        }
        // `BRepLib_MakeFace MkF(thePlane, WW)` (`cxx:709`).
        let built = plane.and_then(|pl| make_plane_face(&pl, &ww));
        let Some(face) = built else {
            // `MkF.Error() != BRepLib_FaceDone` (`cxx:710-718`).
            return None;
        };
        // `cxx:720-724`: suspend the 2D representation on the swept surface.
        for (edge, t) in [(e1, tol1), (e2, tol2), (e3, tol3), (e4, tol4)] {
            reg.update_edge_pcurve_on_surface(&edge.0, s, None, t);
        }
        // `cxx:727-730`.
        update_edge_on_plane(&face, e1);
        update_edge_on_plane(&face, e2);
        update_edge_on_plane(&face, e3);
        update_edge_on_plane(&face, e4);
        f = face;
    } else {
        // `cxx:747-751`.
        let b = TopoBuilder::new();
        f = b.make_face_with_tolerance(s.clone(), &[], CONFUSION);
        b.add_wire(&mut f, &ww);
    }

    // `cxx:753-761`.
    if exch_uv {
        f.0.reverse();
    }
    if u_reverse {
        f.0.reverse();
    }
    Some(f)
}

/// `cxx:608-627` / `:641-657` / `:670-687`: add the edge and, when
/// `BRepBuilderAPI_MakeWire` replaced it with another one, record the new edge
/// in `EEmap`.
fn add_and_record(mw: &mut MakeWire, e: &Edge, ee_map: &mut HashMap<ShapeMapKey, TopoShape>) {
    let _ = mw.add(e);
    if BRepTool::is_degenerated(e) {
        return;
    }
    // "take the last edge added to WW".
    let ww = mw.wire();
    let new_edge = {
        let ts = ww.0.tshape.read().expect("poisoned TShape lock");
        ts.children.last().cloned()
    };
    let Some(new_edge) = new_edge else {
        return;
    };
    if !is_same(&e.0, &new_edge) {
        ee_map.insert(shape_map_key(&e.0), new_edge);
    }
}

/// `BRepLib_MakeFace(const gp_Pln& P, const TopoDS_Wire& W, const bool Inside = true)`
/// (`BRepLib_MakeFace.cxx:262-271`, `:366-385`, `:873-880`, `:907-924`).
fn make_plane_face(plane: &Arc<dyn Surface>, w: &Wire) -> Option<Face> {
    let reg = GeometryRegistry::global();
    let b = TopoBuilder::new();
    // `Init(Pl, false, Precision::Confusion())` (`cxx:264`): `MakeFace` plus
    // `NaturalRestriction(F, true)`.
    let mut face = b.make_face_with_tolerance(plane.clone(), &[], CONFUSION);
    reg.set_natural_restriction(&face.0, true);
    // `Add(W)` (`cxx:265`): add the wire and clear `NaturalRestriction`.
    b.add_wire(&mut face, w);
    reg.set_natural_restriction(&face.0, false);
    // `if (Inside && BRep_Tool::IsClosed(W)) CheckInside();` (`cxx:266-270`).
    if w.is_closed() {
        return Some(check_inside(&face));
    }
    Some(face)
}

/// `BRepLib_MakeFace::CheckInside` (`BRepLib_MakeFace.cxx:907-924`): rebuild the
/// face with every child reversed when the infinite point classifies `IN`.
///
/// `EmptyCopied` copies the surface, location and tolerance
/// (`BRep_TFace::EmptyCopy`, `BRep_TFace.cxx:36-43`) — but not
/// `NaturalRestriction`, which is why the copied [`crate::tgeometry::FaceGeom`]
/// clears it.
fn check_inside(face: &Face) -> Face {
    let fclass = BrepTopAdaptorFClass2d::new(face, 0.0);
    if fclass.perform_infinite_point() != FaceState::In {
        return face.clone();
    }
    let reg = GeometryRegistry::global();
    let mut s = face.0.empty_copied();
    if let Some(mut g) = reg.face_geom(&face.0) {
        g.natural_restriction = false;
        reg.set_face(&s, g);
    }
    let b = TopoBuilder::new();
    let children = {
        let ts = face.0.tshape.read().expect("poisoned TShape lock");
        ts.children.clone()
    };
    for mut c in children {
        c.reverse();
        b.add(&mut s, &c);
    }
    Face(s)
}
