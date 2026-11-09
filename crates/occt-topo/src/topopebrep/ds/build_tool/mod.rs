//! `TopOpeBRepDS_BuildTool` — partial translation.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepDS/TopOpeBRepDS_BuildTool.hxx:31-190`,
//! `.cxx:120-1980` (OCCT 8.0.0).
//!
//! Ported members:
//! - `Parameter(E, V, P)` (`cxx:950-993`) together with the
//!   `BRep_Builder::UpdateVertex(V, P, E, Tol)` call it ends on
//!   (`BRep_Builder.cxx:1220-1311`).
//! - `MakeShell` / `MakeSolid` / `CopyEdge` / `CopyFace` / `AddEdgeVertex` /
//!   `AddWireEdge` / `AddFaceWire` / `AddShellFace` / `AddSolidShell` /
//!   `Range` / `Curve3D` / `Orientation` / `Closed` / `UpdateSurface` (both
//!   overloads) / `TranslateOnPeriodic` / `PCurve(F, E, PC)`.
//! - the `Make*` family and the edge-update family live in the sibling files
//!   [`make_shapes`] (`MakeEdge` ×4, `MakeFace`, `MakeVertex`) and
//!   [`edge_update`] (`GetOrientedEdgeVertices`, `UpdateEdgeCurveTol`,
//!   `UpdateEdge`, `Project`, `Parameter(C, E, V)`, `PCurve(F, E, CDS, PC)`).
//!
//! UNPORTED, with their OCCT location:
//! - `MakeEdge(E, C, BDS)` (`cxx:230-239`) keeps the surface-keyed pcurve of
//!   its null-curve arm unported: the port keys pcurves by face and has no
//!   `BRep_Builder::UpdateEdge(E, PC, Surface, L, Tol)` representation.
//!
//! `RecomputeCurve` (`cxx:1344-1376`), `RecomputeBSpline1Curve`
//! (`cxx:1377-1453`) and `RecomputeCurveOnCone` (`cxx:1454-1580`) are inside a
//! `/* ... */` block in OCCT 8.0.0 (`cxx:1342`, closed at `cxx:1580` by the
//! `- merge 04-07-97` comment), so they are dead code and not ported.
//!
//! `ApproxCurves` (`cxx:471-653`) lives in the sibling file [`approx_curves`];
//! `ComputePCurves` (`cxx:654-789`), `PutPCurves` (`cxx:793-813`) and
//! `RecomputeCurves` (`cxx:816-890`) live in [`recompute_curves`].

mod approx_curves;
mod edge_update;
mod make_shapes;
mod recompute_curves;

use occt_core::elib::clib::in_period;
use occt_core::gp::{GpPnt, GpTrsf2d, GpVec2d};
use occt_core::precision::CONFUSION;
use occt_geom::curve::Curve;
use occt_geom::surface::Surface;
use occt_geom2d::curve::Curve2d;
use std::sync::Arc;

use crate::abs::{Orientation, ShapeType};

use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::tgeometry::{FaceGeom, GeometryRegistry};
use crate::topo_tools::{children_keep_orientation, explorer_shapes, is_same};use crate::topopebrep::tool::geom_tool::GeomTool;
use crate::topopebrep::tool::out_curve_type::OutCurveType;
use crate::topopebrep::tool::shape_tool::adjust_on_periodic;

/// `class TopOpeBRepDS_BuildTool` (`hxx:31`).
///
/// `myBuilder` (`hxx:188`) is a `BRep_Builder`, which has no state; the port
/// keeps `myCurveTool`'s `TopOpeBRepTool_GeomTool` (`hxx:210`, read through
/// `GetGeomTool` / `ChangeGeomTool`) and the two flags `myOverWrite` /
/// `myTranslate` (`hxx:186-187`).
#[derive(Clone, Copy)]
pub struct BuildTool {
    /// `myCurveTool.myGeomTool` (`hxx:210`), `TopOpeBRepTool_GeomTool(APPROX)`
    /// at construction (`cxx:87-92`): the `OutCurveType` overload of
    /// `TopOpeBRepTool_CurveTool` keeps `CompC3D` / `CompPC1` / `CompPC2` at
    /// their `true` defaults (`.cxx:75-79`, `GeomTool.cxx:26-35`).
    geom_tool: GeomTool,
    /// `myOverWrite` (`hxx:186`), `true` at construction (`cxx:89`).
    my_over_write: bool,
    /// `myTranslate` (`hxx:187`), `true` at construction (`cxx:90`).
    my_translate: bool,
}

impl Default for BuildTool {
    fn default() -> Self {
        Self::new()
    }
}

impl BuildTool {
    /// `TopOpeBRepDS_BuildTool()` (`cxx:87-92`).
    pub fn new() -> Self {
        Self {
            geom_tool: GeomTool::new(OutCurveType::Approx, true, true, true),
            my_over_write: true,
            my_translate: true,
        }
    }

    /// `GetGeomTool()` (`cxx:134-138`).
    pub fn geom_tool(&self) -> &GeomTool {
        &self.geom_tool
    }

    /// `ChangeGeomTool()` (`cxx:142-147`).
    pub fn change_geom_tool(&mut self) -> &mut GeomTool {
        &mut self.geom_tool
    }

    /// `Approximation()` (`cxx:1315-1318`): `TypeC3D() != BSPLINE1`.
    pub fn approximation(&self) -> bool {
        self.geom_tool.type_c3d() != OutCurveType::BSpline1
    }

    /// `OverWrite()` (`cxx:112-115`).
    pub fn over_write(&self) -> bool {
        self.my_over_write
    }

    /// `OverWrite(O)` (`cxx:117-120`).
    pub fn set_over_write(&mut self, o: bool) {
        self.my_over_write = o;
    }

    /// `Translate()` (`cxx:122-125`).
    pub fn translate(&self) -> bool {
        self.my_translate
    }

    /// `Translate(T)` (`cxx:127-130`).
    pub fn set_translate(&mut self, t: bool) {
        self.my_translate = t;
    }


    /// `Parameter(E, V, P)` (`cxx:950-993`).
    ///
    /// Maps `P` into the period of a periodic edge curve when `V` sits
    /// `REVERSED` in `E` and `P` lies before the curve's first parameter, then
    /// hands the parameter to `BRep_Builder::UpdateVertex`.
    pub fn parameter(&self, e: &TopoShape, v: &TopoShape, p: f64) {
        let edge = Edge::wrap(e.clone()).expect("TopOpeBRepDS_BuildTool::Parameter: not an edge");
        let vert = Vertex::wrap(v.clone()).expect("TopOpeBRepDS_BuildTool::Parameter: not a vertex");
        let mut p = p;

        // cxx:959-989
        let (f, _l) = BRepTool::edge_parameters(&edge);
        if let Some(c) = BRepTool::edge_curve(&edge) {
            if c.is_periodic() {
                let per = c.period();

                let mut o_v = Orientation::Forward;
                let mut ex_v = explorer_shapes(e, ShapeType::Vertex).into_iter();
                let mut found = false;
                while let Some(vofe) = ex_v.next() {
                    if is_same(&vofe, v) {
                        o_v = vofe.orientation();
                        found = true;
                        break;
                    }
                }
                if found && o_v == Orientation::Reversed && p < f {
                    p = in_period(p, f, f + per);
                }
            }
        }

        // cxx:991-992: `myBuilder.UpdateVertex(v, p, e, 0)`.
        self.update_vertex(&vert, p, &edge);
    }

    /// `MakeVertex(V, P)` (`cxx:150-155`).
    pub fn make_vertex(&self, p: GpPnt, tol: f64) -> TopoShape {
        TopoBuilder::new().make_vertex(p, tol).0
    }

    /// `MakeWire(W)` (`cxx:265-268`) over `TopoDS_Builder::MakeWire`
    /// (`TopoDS_Builder.hxx:77`): a new empty `TopoDS_TWire`, `Location`
    /// identity and `Orientation` FORWARD (`TopoDS_Builder::MakeShape`,
    /// `TopoDS_Builder.cxx:28-32`).
    pub fn make_wire(&self, w: &mut TopoShape) {
        *w = TopoShape::new(ShapeType::Wire);
    }

    /// `MakeShell(Sh)` (`cxx:279-282`).
    pub fn make_shell(&self, sh: &mut TopoShape) {
        *sh = TopoShape::new(ShapeType::Shell);
    }

    /// `MakeSolid(S)` (`cxx:286-289`).
    pub fn make_solid(&self, s: &mut TopoShape) {
        *s = TopoShape::new(ShapeType::Solid);
    }

    /// `CopyEdge(Ein, Eou)` (`cxx:293-317`).
    ///
    /// `Eou = Ein.EmptyCopied()` copies `BRep_TEdge`'s tolerance, every
    /// `BRep_GCurve` / `BRep_CurveOn2Surfaces` representation and the
    /// Degenerated / SameParameter / SameRange flags, but no polygons
    /// (`BRep_TEdge.cxx:105-133`); the port re-registers that side-table
    /// state because the new `TShape` has a fresh key. `cxx:313-316` then
    /// redefines the range on `Eou`.
    pub fn copy_edge(&self, ein: &TopoShape) -> TopoShape {
        let e1 = Edge::wrap(ein.clone()).expect("TopOpeBRepDS_BuildTool::CopyEdge: not an edge");
        let (f, l) = BRepTool::edge_parameters(&e1);
        let eou = ein.empty_copied();
        let reg = GeometryRegistry::global();
        if let Some(g) = reg.edge_geom(&e1.0) {
            reg.set_edge(&eou, g);
        }
        reg.set_ranges_all(&eou, f, l);
        eou
    }

    /// `CopyFace(Fin, Fou)` (`cxx:895-899`) over `BRep_TFace::EmptyCopy`
    /// (`BRep_TFace.cxx:36-43`): the surface, the location and the tolerance.
    pub fn copy_face(&self, fin: &TopoShape) -> TopoShape {
        let fou = fin.empty_copied();
        let reg = GeometryRegistry::global();
        if let Some(g) = reg.face_geom(fin) {
            reg.set_face(&fou, g);
        }
        fou
    }

    /// `AddEdgeVertex(E, V)` (`cxx:915-919`).
    ///
    /// The three-argument overload `AddEdgeVertex(Ein, Eou, V)`
    /// (`cxx:902-911`, which also runs `BRep_Builder::Transfert`) has no
    /// caller in the ported `TopOpeBRepBuild` sources.
    pub fn add_edge_vertex(&self, e: &mut TopoShape, v: &TopoShape) {
        TopoBuilder::new().add(e, v);
    }

    /// `AddWireEdge(W, E)` (`cxx:922-925`).
    pub fn add_wire_edge(&self, w: &mut TopoShape, e: &TopoShape) {
        TopoBuilder::new().add(w, e);
    }

    /// `AddFaceWire(F, W)` (`cxx:929-932`).
    pub fn add_face_wire(&self, f: &mut TopoShape, w: &TopoShape) {
        TopoBuilder::new().add(f, w);
    }

    /// `AddShellFace(Sh, F)` (`cxx:936-939`).
    pub fn add_shell_face(&self, sh: &mut TopoShape, f: &TopoShape) {
        TopoBuilder::new().add(sh, f);
    }

    /// `AddSolidShell(S, Sh)` (`cxx:943-946`).
    pub fn add_solid_shell(&self, s: &mut TopoShape, sh: &TopoShape) {
        TopoBuilder::new().add(s, sh);
    }

    /// `Range(E, first, last)` (`cxx:997-1002`) over
    /// `BRep_Builder::Range` (`BRep_Builder.cxx:1091-1117`): every
    /// `BRep_GCurve` representation gets the range.
    pub fn range(&self, e: &TopoShape, first: f64, last: f64) {
        GeometryRegistry::global().set_ranges_all(e, first, last);
    }

    /// `Curve3D(E, C, Tol)` (`cxx:1069-1073`) over
    /// `BRep_Builder::UpdateEdge(E, C, Tol)` (`BRep_Builder.cxx:1180-1191`).
    pub fn curve3d(&self, e: &TopoShape, c: Arc<dyn Curve>, tol: f64) {
        GeometryRegistry::global().set_edge_curve3d(e, c, tol);
    }

    /// `Orientation(S)` (`cxx:1301-1304`).
    pub fn orientation(&self, s: &TopoShape) -> Orientation {
        s.orientation()
    }

    /// `Orientation(S, O)` (`cxx:1294-1297`).
    pub fn set_orientation(&self, s: &mut TopoShape, o: Orientation) {
        s.set_orientation(o);
    }

    /// `Closed(S, B)` (`cxx:1308-1311`).
    pub fn closed(&self, s: &TopoShape, b: bool) {
        s.set_closed(b);
    }

    /// `UpdateSurface(F, SU)` (`cxx:1320-1327`): `BRep_Builder::UpdateFace`
    /// replaces the surface and keeps the face tolerance.
    pub fn update_surface(&self, f: &TopoShape, su: Arc<dyn Surface>) {
        let reg = GeometryRegistry::global();
        let tolerance = reg.face_tolerance(f);
        let natural_restriction = reg.natural_restriction(f);
        reg.set_face(
            f,
            FaceGeom {
                surface: su,
                tolerance,
                natural_restriction,
            },
        );
    }

    /// `UpdateSurface(E, oldF, newF)` (`cxx:1329-1339`): re-attach the pcurve
    /// of `E` on `oldF` to `newF`, with `oldF`'s tolerance as the new edge
    /// tolerance (`BRep_Builder::UpdateEdge(E, PC, F, Tol)`).
    pub fn update_surface_edge(&self, e: &TopoShape, old_f: &TopoShape, new_f: &TopoShape) {
        let reg = GeometryRegistry::global();
        let old_key = GeometryRegistry::shape_key(old_f);
        let new_key = GeometryRegistry::shape_key(new_f);
        let Some(pc) = reg.edge_pcurve(e, old_key) else {
            return;
        };
        let tol = reg.face_tolerance(old_f);
        reg.set_edge_pcurve(e, new_key, pc);
        reg.raise_edge_tolerance(e, tol);
    }

    /// `TranslateOnPeriodic(F, E, PC)` (`cxx:1077-1144`).
    ///
    /// The pcurve is shifted in `(u, v)` so that the point at the middle of
    /// the edge's 3D range lands on the face's periodic window. `newv` is a
    /// constant `true` in OCCT (`cxx:1103`), so the `!newv` arm (`cxx:1132-1135`)
    /// is unreachable and is kept only as the `else` of the same flag.
    pub fn translate_on_periodic(
        &self,
        f: &Face,
        e: &Edge,
        pc: &Arc<dyn Curve2d>,
    ) -> Arc<dyn Curve2d> {
        let (first, mut last) = BRepTool::edge_parameters(e);
        if let Some(c3d) = BRepTool::edge_curve(e) {
            if c3d.is_periodic() && last < first {
                // cxx:1090-1095
                last += (first - last).abs();
            }
        }

        let uiso_line_on_sphere = fun_u_iso_line_on_sphere(f, pc);
        let newv = true;

        // cxx:1107-1109
        let t = (first + last) * 0.5;
        let ptest = pc.d0(t);
        let u1 = ptest.x();
        let v1 = ptest.y();
        let mut u2 = u1;
        let mut v2 = v1;

        if newv {
            if uiso_line_on_sphere {
                // cxx:1116-1125
                if let (Some(c3d), Some(surf)) = (
                    BRepTool::edge_curve(e),
                    GeometryRegistry::global().face_surface(&f.0),
                ) {
                    let p3dtest = c3d.d0(t);
                    if let Some(ps) = occt_geom::geom_api::project_point_on_surface(
                        surf.as_ref(),
                        &p3dtest,
                        CONFUSION,
                    ) {
                        u2 = ps.u;
                        v2 = ps.v;
                    }
                }
            } else {
                adjust_on_periodic(f, &mut u2, &mut v2);
            }
        }
        let du = u2 - u1;
        let dv = v2 - v1;
        if du != 0.0 || dv != 0.0 {
            let mut tr = GpTrsf2d::identity();
            tr.set_translation_vec(&GpVec2d::new(du, dv));
            return match pc.transformed(&tr) {
                b => Arc::from(b),
            };
        }
        pc.clone()
    }

    /// `PCurve(F, E, PC)` (`cxx:1188-1239`).
    ///
    /// The trailing INTERNAL-vertex reparametrisation (`cxx:1222-1237`) walks
    /// the edge's vertices, keeps the `TopAbs_INTERNAL` ones, reads
    /// `BRep_Tool::Parameter(vi, EE)` and calls
    /// `BRep_Builder::UpdateVertex(vi, newpar, EE, FF, tolvi)`
    /// ([`crate::brep_builder_update_vertex::update_vertex_on_edge_on_surface`]).
    /// `TopOpeBRepTool_ShapeTool::Tolerance(V)` (`TOPOLOGY.cxx:377-383`) is
    /// `BRep_Tool::Tolerance(V)` when the vertex is a `TopoDS_Vertex`.
    ///
    /// The `UPDATEVERTEX_1` / `UPDATEVERTEX_2` path counters and the
    /// `OCCT_DEBUG` print are unported (`cxx:1226-1234`).
    pub fn p_curve(&self, f: &TopoShape, e: &TopoShape, pc: &Arc<dyn Curve2d>) {
        let ff = Face::wrap(f.clone()).expect("TopOpeBRepDS_BuildTool::PCurve: not a face");
        let ee = Edge::wrap(e.clone()).expect("TopOpeBRepDS_BuildTool::PCurve: not an edge");

        // cxx:1196-1211
        let tran = self.my_translate;
        let uiso_line_on_sphere = fun_u_iso_line_on_sphere(&ff, pc);
        let overwrite = if uiso_line_on_sphere {
            false
        } else {
            self.my_over_write
        };
        let mut pct = pc.clone();
        if tran {
            pct = self.translate_on_periodic(&ff, &ee, &pct);
        }

        if overwrite {
            // cxx:1213-1216: `myBuilder.UpdateEdge(EE, PCT, FF, 0)`.
            let face_key = GeometryRegistry::shape_key(f);
            GeometryRegistry::global().set_edge_pcurve(e, face_key, pct);
        } else {
            // cxx:1217-1220
            self.set_the_pcurve(&ee, &ff, ee.0.orientation(), &pct);
        }

        // cxx:1222-1237: reparametrise every INTERNAL vertex on the new pcurve.
        let Some(face_surface) = BRepTool::face_surface(&ff) else {
            return;
        };
        for vi in children_keep_orientation(e) {
            if vi.shape_type() != ShapeType::Vertex || vi.orientation() != Orientation::Internal {
                continue;
            }
            let tolvi = GeometryRegistry::global().vertex_tolerance(&vi);
            let newpar = BRepTool::parameter_of_vertex(e, &vi);
            let _ = crate::brep_builder_update_vertex::update_vertex_on_edge_on_surface(
                &vi,
                newpar,
                e,
                &face_surface,
                tolvi,
            );
        }
    }

    /// `TopOpeBRepDS_SetThePCurve(B, E, F, O, C)` (`cxx:1147-1184`).
    ///
    /// A face with a plane surface always takes `UpdateEdge(E, C, F, tol)`;
    /// otherwise an existing pcurve turns the call into the two-curve
    /// `UpdateEdge(E, C1, C2, F, tol)` (`BRep_Builder.cxx:1205-1240`), the
    /// order depending on `O`.
    pub(super) fn set_the_pcurve(&self, e: &Edge, f: &Face, o: Orientation, c: &Arc<dyn Curve2d>) {
        let reg = GeometryRegistry::global();
        let face_key = GeometryRegistry::shape_key(&f.0);
        let is_plane = reg
            .face_surface(&f.0)
            .map(|s| s.gp_pln().is_some())
            .unwrap_or(false);
        let oc = if is_plane {
            None
        } else {
            reg.edge_pcurve(&e.0, face_key)
        };
        match oc {
            None => reg.set_edge_pcurve(&e.0, face_key, c.clone()),
            Some(oc) => {
                if BRepTool::is_degenerated(e) {
                    return;
                }
                let pair = if o == Orientation::Reversed {
                    vec![oc, c.clone()]
                } else {
                    vec![c.clone(), oc]
                };
                reg.set_edge_pcurves(&e.0, face_key, pair);
            }
        }
    }

    /// `BRep_Builder::UpdateVertex(V, Par, E, Tol)` (`BRep_Builder.cxx:1220-1311`).
    ///
    /// `Tol` is 0 at the only call site (`BuildTool.cxx:992`), and
    /// `BRep_TVertex::UpdateTolerance(T)` (`BRep_TVertex.cxx:29-35`) only raises
    /// the tolerance, so the tolerance call has no effect; the shared port
    /// ([`crate::brep_builder_update_vertex::update_vertex_on_edge`]) runs it
    /// anyway, because a raise to 0 is a no-op.
    fn update_vertex(&self, v: &Vertex, par: f64, e: &Edge) {
        let _ = crate::brep_builder_update_vertex::update_vertex_on_edge(&v.0, par, &e.0, 0.0);
    }
}

/// `FUN_UisoLineOnSphe(F, PC)` (`cxx:62-83`): true when `PC` is a 2D line
/// parallel to the `u` axis and the face's basis surface is a sphere, i.e. the
/// pcurve is a `u`-isoline of a sphere.
///
/// `TopOpeBRepTool_ShapeTool::BASISSURFACE` unwraps one
/// `Geom_RectangularTrimmedSurface` level, and `BASISCURVE2D` one
/// `Geom2d_TrimmedCurve` level.
fn fun_u_iso_line_on_sphere(f: &Face, pc: &Arc<dyn Curve2d>) -> bool {
    let reg = GeometryRegistry::global();
    let Some(surface) = reg.face_surface(&f.0) else {
        return false;
    };
    let basis = surface
        .rectangular_trimmed_basis()
        .unwrap_or_else(|| surface.clone());
    let curve_basis: &dyn Curve2d = pc.trimmed_basis().unwrap_or(pc.as_ref());

    let is_sphere = basis.gp_sphere().is_some();
    let is_line_2d = curve_basis.is_geom2d_line();
    if !(is_sphere && is_line_2d) {
        return false;
    }
    let Some(line) = curve_basis.gp_lin2d() else {
        return false;
    };
    // `Precision::Parametric(Precision::Confusion())` = `Confusion * 0.01`
    // (`Precision.hxx:328`).
    let d = *line.direction();
    d.x().abs() < CONFUSION * 0.01
}
