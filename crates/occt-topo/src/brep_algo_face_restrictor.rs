//! `BRepAlgo_FaceRestrictor` — rebuild the faces bounded by a set of wires on
//! a reference face.
//!
//! Source: `ModelingAlgorithms/TKBool/BRepAlgo/BRepAlgo_FaceRestrictor.cxx:39-454`
//! and `.hxx` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `BRepAlgo_FaceRestrictor()` | `cxx:39` | [`BRepAlgoFaceRestrictor::new`] |
//! | `Init(F, Proj, CorrectionOrientation)` | `cxx:44-51` | [`BRepAlgoFaceRestrictor::init`] |
//! | `Add(W)` | `cxx:55-58` | [`BRepAlgoFaceRestrictor::add`] |
//! | `Clear` | `cxx:62-66` | [`BRepAlgoFaceRestrictor::clear`] |
//! | `ChangePCurve` | `cxx:69-83` | [`change_pcurve`] |
//! | `ProjCurve3d` | `cxx:87-105` | [`proj_curve3d`] |
//! | `Perform` | `cxx:109-180` | [`BRepAlgoFaceRestrictor::perform`] |
//! | `IsDone` / `More` / `Next` / `Current` | `cxx:184-207` | same names |
//! | `IsClosed` | `cxx:210-222` | [`is_closed_wire`] |
//! | `IsInside` | `cxx:226-265` | [`is_inside`] |
//! | `Store` | `cxx:269-289` | [`store`] |
//! | `BuildFaceIn` | `cxx:292-366` | [`build_face_in`] |
//! | `PerformWithCorrection` | `cxx:370-454` | [`BRepAlgoFaceRestrictor::perform_with_correction`] |
//!
//! `Perform` is the path `BRepOffsetAPI_MakeOffset::BuildDomains` takes
//! (`BRepOffsetAPI_MakeOffset.cxx:214`, `FR.Init(myFace, true)`): every wire
//! must carry a p-curve on the reference surface, `WireToFace` then rebuilds
//! the region faces. `PerformWithCorrection` is the other mode and classifies
//! the wires against each other (`cxx:370-454`).
//!
//! Translation notes:
//!
//! * `myFace` is a `TopoDS_Face` that can be null; here it is
//!   `Option<Face>` and a null face makes `Perform` answer "not done" at the
//!   `BRep_Tool::Surface` step (`cxx:125`), which is what a null handle does
//!   to `ChangePCurve` / `ProjCurve3d` in OCCT.
//! * `BRep_Tool::CurveOnSurface(E, S, L, f, l)` is keyed by the `Geom_Surface`
//!   handle; [`pcurve_on_surface`] reads the registry's CurveOnSurface
//!   representations and matches the surface by pointer.
//! * `BRep_Builder::UpdateEdge(E, C2, S, L, Tol)` is
//!   [`GeometryRegistry::update_edge_pcurve_on_surface`].
//! * `L.Inverted().Multiplied(LE)` of `ProjCurve3d` (`cxx:99-100`) is the
//!   identity in the port: the registry stores 3D curves and surfaces in one
//!   frame, so `BRep_Tool::Curve` and `BRep_Tool::Surface` report the same
//!   location and the composition vanishes. The `Geom_TrimmedCurve::Transform`
//!   is therefore not applied.
//! * `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>` over
//!   `TopTools_ShapeMapHasher` is `HashMap<ShapeMapKey, Vec<TopoShape>>`
//!   (`IsSame` identity), and `NCollection_List` is `Vec` in insertion order.
//! * `BRepTopAdaptor_FClass2d` is [`crate::brep_top_adaptor_fclass2d`], whose
//!   `perform` answers [`FaceState`].
//! * `TopOpeBRepBuild_WireToFace` is [`crate::topopebrep::build::WireToFace`];
//!   it needs a `Builder`, which the caller of `MakeFaces` owns in OCCT
//!   (`BRepAlgo_FaceRestrictor.cxx:136-177` creates a local `WTF` only), so
//!   [`BRepAlgoFaceRestrictor::perform`] builds a throwaway one on
//!   `BuildTool::new()` exactly where OCCT builds its throwaway `Builder`.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::GpPnt2d;
use occt_core::precision::Precision;
use occt_geom::Surface;
use occt_geom2d::curve::Curve2d;

use crate::abs::{Orientation, ShapeType};
use crate::boptools_2d::curve_on_surface_range;
use crate::brep_top_adaptor_fclass2d::BrepTopAdaptorFClass2d;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key, wire_vertices, ShapeMapKey};
use crate::topopebrep::build::{Builder, WireToFace};
use crate::topopebrep::ds::BuildTool;

/// `class BRepAlgo_FaceRestrictor` (`hxx`).
#[derive(Clone, Default)]
pub struct BRepAlgoFaceRestrictor {
    /// `myFace` (`hxx:70`).
    my_face: Option<Face>,
    /// `modeProj` (`hxx:72`).
    mode_proj: bool,
    /// `myCorrection` (`hxx:74`).
    my_correction: bool,
    /// `myDone` (`hxx:76`).
    my_done: bool,
    /// `wires` (`hxx:78`).
    wires: Vec<TopoShape>,
    /// `faces` (`hxx:80`).
    faces: Vec<TopoShape>,
    /// `keyIsIn` (`hxx:82`): each key maps to the wires it contains.
    key_is_in: HashMap<ShapeMapKey, Vec<TopoShape>>,
    /// `keyContains` (`hxx:84`): each key maps to the wires containing it.
    key_contains: HashMap<ShapeMapKey, Vec<TopoShape>>,
}

impl BRepAlgoFaceRestrictor {
    /// `BRepAlgo_FaceRestrictor()` (`cxx:39`): `= default`, i.e. a null face,
    /// `modeProj` / `myCorrection` / `myDone` false and empty lists.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Init(F, Proj, CorrectionOrientation)` (`cxx:44-51`).
    pub fn init(&mut self, f: Option<Face>, proj: bool, correction_orientation: bool) {
        // `cxx:48-50`
        self.my_face = f;
        self.mode_proj = proj;
        self.my_correction = correction_orientation;
    }

    /// `Add(W)` (`cxx:55-58`): `wires.Append(W)`.
    ///
    /// OCCT takes the wire by reference and never modifies it here, so the
    /// port clones.
    pub fn add(&mut self, w: &TopoShape) {
        self.wires.push(w.clone());
    }

    /// `Clear()` (`cxx:62-66`).
    pub fn clear(&mut self) {
        self.wires.clear();
        self.faces.clear();
    }

    /// `IsDone()` (`cxx:184-187`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `More()` (`cxx:191-194`): `!faces.IsEmpty()`.
    ///
    /// OCCT removes the head on every `Next`, so `More` / `Current` / `Next`
    /// walk a destructive list; the port keeps `faces` and walks it with a
    /// cursor-free `first`.
    pub fn more(&self) -> bool {
        !self.faces.is_empty()
    }

    /// `Next()` (`cxx:197-200`): `faces.RemoveFirst()`.
    pub fn next(&mut self) {
        if !self.faces.is_empty() {
            self.faces.remove(0);
        }
    }

    /// `Current()` (`cxx:204-207`): `TopoDS::Face(faces.First())`.
    pub fn current(&self) -> Option<Face> {
        self.faces.first().map(|f| Face(f.clone()))
    }

    /// The whole `faces` list, for callers that iterate it themselves.
    pub fn faces(&self) -> &[TopoShape] {
        &self.faces
    }

    /// `Perform()` (`cxx:109-180`).
    pub fn perform(&mut self) {
        // `cxx:111-115`
        if self.my_correction {
            self.perform_with_correction();
            return;
        }

        self.my_done = false;

        // `cxx:125`: `BRep_Tool::Surface(myFace, L)`. A null face or a face
        // without a surface leaves nothing to project onto; OCCT keeps going
        // with a null handle and `ChangePCurve` / `ProjCurve3d` then abort the
        // loop, so the port stops before the loop.
        let Some(face) = self.my_face.clone() else {
            return;
        };
        let Some(surface) = GeometryRegistry::global().face_surface(&face.0) else {
            return;
        };

        // `cxx:134-136`
        let mut wtf = WireToFace::new();

        // `cxx:138-176`
        for w in &self.wires {
            // `cxx:142`: `for (Exp.Init(W, TopAbs_EDGE); ...)`.
            for e in explorer_shapes(w, ShapeType::Edge) {
                let edge = Edge(e);
                if pcurve_on_surface(&edge, &surface).is_none() {
                    if self.mode_proj {
                        // `cxx:151-158`: project the 3D curve onto the surface.
                        if !proj_curve3d(&edge, &surface) {
                            return;
                        }
                    } else {
                        // `cxx:159-171`: re-glue the first p-curve on `S`,
                        // else project.
                        let no_pcurve = change_pcurve(&edge, &surface);
                        if no_pcurve && !proj_curve3d(&edge, &surface) {
                            return;
                        }
                    }
                }
            }
            // `cxx:174`
            wtf.add_wire(w);
        }

        // `cxx:175`
        // `cxx:176`: `WTF.MakeFaces(myFace, faces)`.
        let mut builder = Builder::new(BuildTool::new());
        wtf.make_faces(&face.0, &mut builder, &mut self.faces);

        // `cxx:178`
        self.my_done = true;
    }

    /// `PerformWithCorrection()` (`cxx:370-454`).
    pub fn perform_with_correction(&mut self) {
        // `cxx:372-374`
        self.my_done = false;
        let builder = TopoBuilder::new();
        let reg = GeometryRegistry::global();
        let Some(face) = self.my_face.clone() else {
            return;
        };
        let face_geom = reg.face_geom(&face.0);

        // `cxx:377-396`: reorient every closed wire to the left.
        for i in 0..self.wires.len() {
            let mut w = self.wires[i].clone();
            // `cxx:382-386`: `NF = myFace.EmptyCopied()` with
            // `Orientation(FORWARD)`, then `B.Add(NF, W)`.
            let mut nf = face.0.empty_copied();
            if let Some(g) = &face_geom {
                reg.set_face(&nf, g.clone());
            }
            nf.set_orientation(Orientation::Forward);
            builder.add(&mut nf, &w);

            // `cxx:388-395`
            if is_closed_wire(&w) {
                let fclass = BrepTopAdaptorFClass2d::new(&Face(nf), Precision::PCONFUSION);
                if fclass.perform_infinite_point() != FaceState::Out {
                    // `cxx:393`: `W.Reverse()`.
                    w = w.oriented(w.orientation().reversed());
                }
            }
            self.wires[i] = w;
        }

        // `cxx:399-423`: classify the wires against the closed ones.
        for i in 0..self.wires.len() {
            let w1 = self.wires[i].clone();
            if !is_closed_wire(&w1) {
                continue;
            }
            // `cxx:407-412`
            let mut nf = face.0.empty_copied();
            if let Some(g) = &face_geom {
                reg.set_face(&nf, g.clone());
            }
            nf.set_orientation(Orientation::Forward);
            let mut nf_shape = nf;
            builder.add(&mut nf_shape, &w1);
            let nf_face = Face(nf_shape);

            let fclass = BrepTopAdaptorFClass2d::new(&nf_face, Precision::PCONFUSION);
            for j in 0..self.wires.len() {
                let w2 = self.wires[j].clone();
                // `cxx:418`: `!W1.IsSame(W2) && IsInside(W2, NF, FClass2d)`.
                if !is_same(&w1, &w2) && is_inside(&w2, &nf_face, &fclass) {
                    store(&mut self.key_is_in, &mut self.key_contains, &w2, &w1);
                }
            }
        }

        // `cxx:425-433`: the exterior wires.
        let mut wire_ext: Vec<TopoShape> = Vec::new();
        for w in &self.wires {
            let k = shape_map_key(w);
            let empty = self.key_is_in.get(&k).map(|v| v.is_empty()).unwrap_or(true);
            if empty {
                wire_ext.push(w.clone());
            }
        }

        // `cxx:435-453`
        for w in &wire_ext {
            let k = shape_map_key(w);
            let empty = self.key_is_in.get(&k).map(|v| v.is_empty()).unwrap_or(true);
            if !empty {
                continue;
            }
            // `cxx:440-446`: `NewFace = myFace.EmptyCopied()` FORWARD, add W.
            let mut new_face = face.0.empty_copied();
            if let Some(g) = &face_geom {
                reg.set_face(&new_face, g.clone());
            }
            new_face.set_orientation(Orientation::Forward);
            builder.add(&mut new_face, w);
            // `cxx:447`
            self.faces.push(new_face.clone());
            // `cxx:451`
            build_face_in(
                Face(new_face),
                w,
                &mut self.key_contains,
                &mut self.key_is_in,
                Orientation::Forward,
                &mut self.faces,
                &face_geom,
            );
        }

        // `cxx:454`
        self.my_done = true;
    }
}

/// The curve-on-surface representation of `edge` on `surface`, i.e.
/// `BRep_Tool::CurveOnSurface(E, S, L, f, l, 1)` (`BRep_Tool.cxx:347-357`):
/// the first representation whose `Geom_Surface` handle is `surface`.
fn pcurve_on_surface(edge: &Edge, surface: &Arc<dyn Surface>) -> Option<(Arc<dyn Curve2d>, f64, f64)> {
    let reg = GeometryRegistry::global();
    reg.edge_pcurve_reps(&edge.0)
        .into_iter()
        .find(|(s, _, _, _)| Arc::ptr_eq(s, surface))
        .map(|(_, c, f, l)| (c, f, l))
}

/// `ChangePCurve(E, S, L)` (`BRepAlgo_FaceRestrictor.cxx:69-83`).
///
/// Reads the edge's first CurveOnSurface representation (on any surface) and
/// re-attaches it to `S`; answers `true` when there was none.
fn change_pcurve(edge: &Edge, surface: &Arc<dyn Surface>) -> bool {
    let reg = GeometryRegistry::global();
    // `cxx:78`: `BRep_Tool::CurveOnSurface(E, C2, SE, LE, f, l, 1)`.
    let c2 = reg.edge_pcurve_reps(&edge.0).into_iter().next().map(|(_, c, _, _)| c);
    let is_null = c2.is_none();
    if let Some(c) = c2 {
        // `cxx:81`: `BB.UpdateEdge(E, C2, S, L, Precision::Confusion())`.
        reg.update_edge_pcurve_on_surface(&edge.0, surface, Some(c), Precision::CONFUSION);
    }
    // `cxx:82`: `return (C2.IsNull())`.
    is_null
}

/// `ProjCurve3d(E, S, L)` (`BRepAlgo_FaceRestrictor.cxx:87-105`).
fn proj_curve3d(edge: &Edge, surface: &Arc<dyn Surface>) -> bool {
    let reg = GeometryRegistry::global();
    // `cxx:93`: `C = BRep_Tool::Curve(E, LE, f, l)`.
    let Some(curve) = reg.edge_curve(&edge.0) else {
        // `cxx:94-97`
        return false;
    };
    let (f, l) = reg.edge_parameters(&edge.0);

    // `cxx:98-100`: `CT = Geom_TrimmedCurve(C, f, l)` then
    // `CT->Transform(LL.Transformation())`. `LL = L.Inverted().Multiplied(LE)`
    // is the identity in the port (see the module header).
    let trimmed: Arc<dyn occt_geom::Curve> =
        Arc::new(occt_geom::GeomTrimmedCurve::new(curve, f, l));

    // `cxx:102`: `C2 = GeomProjLib::Curve2d(CT, S)`.
    let c2 = occt_geom::geom_proj_lib::curve_2d(&trimmed, surface);

    // `cxx:103`: `BB.UpdateEdge(E, C2, S, L, Precision::Confusion())`. A null
    // `C2` removes the representation on `S` (`BRep_Builder.cxx:149`); OCCT
    // returns `true` either way.
    reg.update_edge_pcurve_on_surface(&edge.0, surface, c2, Precision::CONFUSION);

    // `cxx:104`
    true
}

/// `IsClosed(W)` (`BRepAlgo_FaceRestrictor.cxx:210-222`).
fn is_closed_wire(w: &TopoShape) -> bool {
    if w.is_closed() {
        // `cxx:214-217`
        return true;
    }
    // `cxx:218-221`
    let (v1, v2) = wire_vertices(&crate::shape::Wire(w.clone()));
    match (v1, v2) {
        (Some(a), Some(b)) => is_same(&a.0, &b.0),
        _ => false,
    }
}

/// `IsInside(wir, F, FClass2d)` (`BRepAlgo_FaceRestrictor.cxx:226-265`).
///
/// The classifier OCCT passes in is ignored by the signature
/// (`BRepTopAdaptor_FClass2d& /*FClass2d*/`) and a fresh one is built on the
/// empty face `F` (`cxx:259`), so the port ignores it too.
fn is_inside(wir: &TopoShape, f: &Face, _fclass: &BrepTopAdaptorFClass2d) -> bool {
    // `cxx:229-231`
    let edges = explorer_shapes(wir, ShapeType::Edge);
    let Some(e) = edges.first() else {
        // `cxx:264`
        return false;
    };
    let edge = Edge(e.clone());

    // `cxx:235`: `BRep_Tool::CurveOnSurface(edg, F, f, l)`.
    let Some((c2d, f_par, l_par)) = curve_on_surface_range(&edge, f) else {
        return false;
    };

    // `cxx:238-255`: the evaluation parameter.
    let prm = if !Precision::is_negative_infinite(f_par) && !Precision::is_positive_infinite(l_par)
    {
        (f_par + l_par) / 2.0
    } else if Precision::is_negative_infinite(f_par) && Precision::is_positive_infinite(l_par) {
        0.0
    } else if Precision::is_negative_infinite(f_par) {
        l_par - 1.0
    } else {
        f_par + 1.0
    };

    // `cxx:258`
    let pt2d: GpPnt2d = c2d.d0(prm);

    // `cxx:259-261`
    let fclass = BrepTopAdaptorFClass2d::new(f, Precision::PCONFUSION);
    let st2 = fclass.perform(pt2d, false);
    st2 == FaceState::In
}

/// `Store(W2, W1, keyIsIn, keyContains)` (`BRepAlgo_FaceRestrictor.cxx:269-289`).
fn store(
    key_is_in: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    key_contains: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    w2: &TopoShape,
    w1: &TopoShape,
) {
    // `cxx:278-281`
    key_is_in.entry(shape_map_key(w2)).or_default().push(w1.clone());
    // `cxx:282-288`
    key_contains.entry(shape_map_key(w1)).or_default().push(w2.clone());
}

/// `BuildFaceIn(F, W, KeyContains, KeyIsIn, Orientation, Faces)`
/// (`BRepAlgo_FaceRestrictor.cxx:292-366`).
#[allow(clippy::too_many_arguments)]
fn build_face_in(
    f: Face,
    w: &TopoShape,
    key_contains: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    key_is_in: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    orientation: Orientation,
    faces: &mut Vec<TopoShape>,
    face_geom: &Option<crate::tgeometry::FaceGeom>,
) {
    // `cxx:302-305`
    let Some(contained) = key_contains.get(&shape_map_key(w)).cloned() else {
        return;
    };
    if contained.is_empty() {
        return;
    }

    // `cxx:312-327`: remove `W` from every `KeyIsIn(WI)`.
    for wi in &contained {
        let k = shape_map_key(wi);
        if let Some(l2) = key_is_in.get_mut(&k) {
            if let Some(pos) = l2.iter().position(|x| is_same(x, w)) {
                l2.remove(pos);
            }
        }
    }

    // `cxx:329-339`: the contained wires no one else contains.
    let mut wire_ext: Vec<TopoShape> = Vec::new();
    for wi in &contained {
        let k = shape_map_key(wi);
        let empty = key_is_in.get(&k).map(|v| v.is_empty()).unwrap_or(true);
        if empty {
            wire_ext.push(wi.clone());
        }
    }

    // `cxx:341-365`
    for wi in &wire_ext {
        let k = shape_map_key(wi);
        let empty = key_is_in.get(&k).map(|v| v.is_empty()).unwrap_or(true);
        if !empty {
            continue;
        }
        let reg = GeometryRegistry::global();
        let builder = TopoBuilder::new();
        if orientation == Orientation::Forward {
            // `cxx:348-353`: add WI reversed to F, then recurse REVERSED.
            let mut nwi = wi.clone();
            nwi = nwi.oriented(nwi.orientation().reversed());
            let mut ff = f.0.clone();
            builder.add(&mut ff, &nwi);
            build_face_in(
                Face(ff),
                wi,
                key_contains,
                key_is_in,
                Orientation::Reversed,
                faces,
                face_geom,
            );
        } else {
            // `cxx:355-364`: a new empty face over the same surface.
            let Some(first) = faces.first().cloned() else {
                return;
            };
            let mut nf = first.empty_copied();
            if let Some(g) = face_geom {
                reg.set_face(&nf, g.clone());
            }
            builder.add(&mut nf, wi);
            faces.push(nf.clone());
            build_face_in(
                Face(nf),
                wi,
                key_contains,
                key_is_in,
                Orientation::Forward,
                faces,
                face_geom,
            );
        }
    }
}
