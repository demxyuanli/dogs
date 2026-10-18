//! `BRepLib::EncodeRegularity` and `BRepLib::ContinuityOfFaces`.
//!
//! Source: `BRepLib.cxx:2105-2590`. `EncodeRegul` (`IGESToBRep_Actor.cxx:44-66`)
//! calls `EncodeRegularity(S, read.encoderegularity.angle)` with the default
//! angle `0.01` (`Interface_StaticStandards.cxx:51`).
//!
//! A failed normal or curvature is the `Standard_Failure` catch at
//! `BRepLib.cxx:2582-2588`: that edge is left unmarked. A surface whose
//! `Geom_Surface::D2` is only the port finite-difference default
//! (`occt-geom/src/surface.rs`) is treated the same way when the walk needs
//! second derivatives (`GeomLProp_SurfaceUtils.hxx:376` and the
//! `Extrema_GFuncExtPC` derivative at `Extrema_GFuncExtPC.hxx:344`).

use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

use occt_core::cslib::normal::{normal_d1_mag, CSLibNormalStatus};
use occt_core::gp::{GpDir, GpPnt, GpPnt2d, GpTrsf, GpVec, GpVec2d};
use occt_core::math_direct_poly_roots::DirectPolynomialRoots;
use occt_core::precision::{CONFUSION, PCONFUSION, SQUARE_CONFUSION};
use occt_geom::Surface;
use occt_geom::extrema_pc::LocateExtPC;
use occt_geom2d::Curve2d;

use crate::abs::ShapeType;
use crate::boptools_2d::curve_on_surface_oriented;
use crate::brep_tool::BRepTool;
use crate::meshing::edge_discret::CurveOnSurface;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::{edges_of_wire, faces_of, is_equal, is_same, wires_of_face};

/// `GeomAbs_Shape` (`GeomAbs_Shape.hxx`). Order matches the comparisons in
/// `ContinuityOfFaces` (`BRepLib.cxx:2300-2408`).
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum AbsShape {
    C0 = 0,
    G1 = 1,
    C1 = 2,
    G2 = 3,
    C2 = 4,
    /// Kept so `Cn` stays `GeomAbs_CN` (6). `ContinuityOfFaces` does not
    /// return `C3`.
    #[allow(dead_code)]
    C3 = 5,
    Cn = 6,
}

/// `(edge, face_a, face_b)` with `face_a <= face_b`. `BRep_Tool::Continuity`
/// matches either face order (`BRep_Tool.cxx` continuity query).
fn pair_key(edge: &Edge, f1: &Face, f2: &Face) -> (usize, usize, usize) {
    let e = GeometryRegistry::shape_key(&edge.0);
    let mut a = GeometryRegistry::shape_key(&f1.0);
    let mut b = GeometryRegistry::shape_key(&f2.0);
    if a > b {
        std::mem::swap(&mut a, &mut b);
    }
    (e, a, b)
}

fn continuity_store() -> &'static Mutex<HashMap<(usize, usize, usize), u8>> {
    static STORE: OnceLock<Mutex<HashMap<(usize, usize, usize), u8>>> = OnceLock::new();
    STORE.get_or_init(|| Mutex::new(HashMap::new()))
}

/// Stored face-pair continuity, or `GeomAbs_C0` when the edge was never encoded.
pub fn stored_continuity(edge: &Edge, f1: &Face, f2: &Face) -> u8 {
    let guard = continuity_store()
        .lock()
        .unwrap_or_else(|err| err.into_inner());
    guard.get(&pair_key(edge, f1, f2)).copied().unwrap_or(0)
}

/// `BRep_Builder::Continuity(E, F1, F2, C)` (`BRep_Builder.cxx:1012-1021`).
/// `TShape::Locked` (`cxx:1034-1036`) is not stored.
pub fn set_stored_continuity(edge: &Edge, f1: &Face, f2: &Face, cont: u8) {
    let shape = match cont {
        1 => AbsShape::G1,
        2 => AbsShape::C1,
        3 => AbsShape::G2,
        4 => AbsShape::C2,
        5 => AbsShape::C3,
        6 => AbsShape::Cn,
        _ => AbsShape::C0,
    };
    write_continuity(edge, f1, f2, shape);
}

fn write_continuity(edge: &Edge, f1: &Face, f2: &Face, cont: AbsShape) {
    let mut guard = continuity_store()
        .lock()
        .unwrap_or_else(|err| err.into_inner());
    guard.insert(pair_key(edge, f1, f2), cont as u8);
}

fn unwrap_trimmed(surf: std::sync::Arc<dyn Surface>) -> std::sync::Arc<dyn Surface> {
    surf.rectangular_trimmed_basis().unwrap_or(surf)
}

/// Elementary surfaces are plane, cylinder, cone, sphere, torus
/// (`Geom_ElementarySurface`). A surface of revolution is not elementary.
fn is_elementary(surf: &dyn Surface) -> bool {
    surf.gp_pln().is_some()
        || surf.gp_cylinder().is_some()
        || surf.gp_cone().is_some()
        || surf.gp_sphere().is_some()
        || surf.gp_torus().is_some()
}

/// True when `Surface::d2` is an analytic override, not the finite-difference
/// default in `occt-geom/src/surface.rs`.
fn has_analytic_d2(surf: &dyn Surface) -> bool {
    if is_elementary(surf)
        || surf.is_bspline_surface()
        || surf.is_surface_of_revolution()
        || surf.is_surface_of_linear_extrusion()
    {
        return true;
    }
    if let Some(basis) = surf.rectangular_trimmed_basis() {
        return has_analytic_d2(basis.as_ref());
    }
    false
}

fn dir_transformed(d: &GpDir, t: &GpTrsf) -> Option<GpDir> {
    let v = GpVec::new(d.x(), d.y(), d.z()).transformed(t);
    GpDir::from_vec(&v).ok()
}

struct FaceProps<'a> {
    surf: &'a dyn Surface,
    c2d: &'a dyn Curve2d,
    trsf: &'a GpTrsf,
    reversed: bool,
    analytic_d2: bool,
    uv: GpPnt2d,
    tangent: GpVec2d,
    d1u: GpVec,
    d1v: GpVec,
}

impl<'a> FaceProps<'a> {
    fn new(
        surf: &'a dyn Surface,
        c2d: &'a dyn Curve2d,
        trsf: &'a GpTrsf,
        reversed: bool,
    ) -> Self {
        Self {
            surf,
            c2d,
            trsf,
            reversed,
            analytic_d2: has_analytic_d2(surf),
            uv: GpPnt2d::new(0.0, 0.0),
            tangent: GpVec2d::zero(),
            d1u: GpVec::new(0.0, 0.0, 0.0),
            d1v: GpVec::new(0.0, 0.0, 0.0),
        }
    }

    /// `SurfaceProperties::Calculate` (`BRepLib.cxx:2120-2125`).
    fn calculate(&mut self, t: f64) {
        let (uv, tan) = self.c2d.d1(t);
        let (_p, d1u, d1v) = self.surf.d1(uv.x(), uv.y());
        self.uv = uv;
        self.tangent = tan;
        self.d1u = d1u;
        self.d1v = d1v;
    }

    /// `SurfaceProperties::Value` (`BRepLib.cxx:2128`).
    fn value(&self) -> GpPnt {
        self.surf
            .d0(self.uv.x(), self.uv.y())
            .transformed(self.trsf)
    }

    /// `SurfaceProperties::Derivative` (`BRepLib.cxx:2131-2148`).
    fn derivative(&self) -> GpVec {
        let mut ortho = GpVec2d::new(-self.tangent.y(), self.tangent.x());
        let len = ortho.magnitude();
        if len < CONFUSION {
            return GpVec::new(0.0, 0.0, 0.0);
        }
        ortho = ortho.divided(len);
        if self.reversed {
            ortho.reverse();
        }
        let deriv = self
            .d1u
            .multiplied_scalar(ortho.x())
            .added(&self.d1v.multiplied_scalar(ortho.y()));
        deriv.transformed(self.trsf)
    }

    /// `SurfaceProperties::Normal` (`BRepLib.cxx:2150-2155`).
    /// `GeomLProp_SLProps` is built with `Precision::Confusion()`
    /// (`BRepLib.cxx:2112`), passed to `CSLib::Normal` as `MagTol`
    /// (`GeomLProp_SurfaceUtils.hxx:350-357`). Undefined raises
    /// `LProp_NotDefined`.
    fn normal(&self) -> Option<GpDir> {
        let (status, dir) = normal_d1_mag(&self.d1u, &self.d1v, CONFUSION);
        if status != CSLibNormalStatus::Defined {
            return None;
        }
        dir_transformed(&dir?, self.trsf)
    }

    /// `SurfaceProperties::Curvature` (`BRepLib.cxx:2159-2179`).
    /// Returns `(max_dir, max_curv, min_dir, min_curv)`.
    fn curvature(&self) -> Option<(GpDir, f64, GpDir, f64)> {
        if !self.analytic_d2 {
            return None;
        }
        let (_p, _d1u, _d1v, d2u, d2v, d2uv) = self.surf.d2(self.uv.x(), self.uv.y());
        let (status, normal) = normal_d1_mag(&self.d1u, &self.d1v, CONFUSION);
        if status != CSLibNormalStatus::Defined {
            return None;
        }
        let normal = normal?;
        if self.d1u.magnitude() <= CONFUSION || self.d1v.magnitude() <= CONFUSION {
            return None;
        }
        let (min_c, max_c, dir_min, dir_max) =
            compute_surf_curvatures(&self.d1u, &self.d1v, &d2u, &d2v, &d2uv, &normal)?;
        let mut c1 = max_c;
        let mut c2 = min_c;
        if self.reversed {
            c1 = -c1;
            c2 = -c2;
        }
        if self.trsf.is_negative() {
            c1 = -c1;
            c2 = -c2;
        }
        let d1 = dir_transformed(&dir_max, self.trsf)?;
        let d2 = dir_transformed(&dir_min, self.trsf)?;
        Some((d1, c1, d2, c2))
    }
}

/// `GeomLProp_SurfaceUtils::ComputeSurfCurvatures` (`hxx:376-479`).
/// `None` is the `return false` path (curvature undefined).
fn compute_surf_curvatures(
    d1u: &GpVec,
    d1v: &GpVec,
    d2u: &GpVec,
    d2v: &GpVec,
    duv: &GpVec,
    normal: &GpDir,
) -> Option<(f64, f64, GpDir, GpDir)> {
    let norm = GpVec::new(normal.x(), normal.y(), normal.z());
    let e = d1u.square_magnitude();
    let f = d1u.dot(d1v);
    let g = d1v.square_magnitude();
    let l = norm.dot(d2u);
    let m = norm.dot(duv);
    let n = norm.dot(d2v);
    let a0 = e * m - f * l;
    let b0 = e * n - g * l;
    let c0 = f * n - g * m;
    let max_abc = a0.abs().max(b0.abs()).max(c0.abs());
    let eps = f64::EPSILON;
    if max_abc < eps {
        if g < eps {
            return None;
        }
        let curv = n / g;
        let dir_min = GpDir::from_vec(d1u).ok()?;
        let crossed = d1u.crossed(&norm);
        let dir_max = GpDir::from_vec(&crossed).ok()?;
        return Some((curv, curv, dir_min, dir_max));
    }
    let a = a0 / max_abc;
    let b = b0 / max_abc;
    let c = c0 / max_abc;
    let (c1, c2, v1, v2) = if a.abs() > eps {
        let roots = DirectPolynomialRoots::new3(a, b, c);
        if !roots.is_done() || roots.infinite_roots() || roots.nb_solutions() != 2 {
            return None;
        }
        let r1 = roots.value(1);
        let r2 = roots.value(2);
        let curv1 = ((l * r1 + 2.0 * m) * r1 + n) / ((e * r1 + 2.0 * f) * r1 + g);
        let curv2 = ((l * r2 + 2.0 * m) * r2 + n) / ((e * r2 + 2.0 * f) * r2 + g);
        let vect1 = d1u.multiplied_scalar(r1).added(d1v);
        let vect2 = d1u.multiplied_scalar(r2).added(d1v);
        (curv1, curv2, vect1, vect2)
    } else if c.abs() > eps {
        let roots = DirectPolynomialRoots::new3(c, b, a);
        if !roots.is_done() || roots.infinite_roots() || roots.nb_solutions() != 2 {
            return None;
        }
        let r1 = roots.value(1);
        let r2 = roots.value(2);
        let curv1 = ((n * r1 + 2.0 * m) * r1 + l) / ((g * r1 + 2.0 * f) * r1 + e);
        let curv2 = ((n * r2 + 2.0 * m) * r2 + l) / ((g * r2 + 2.0 * f) * r2 + e);
        let vect1 = d1u.added(&d1v.multiplied_scalar(r1));
        let vect2 = d1u.added(&d1v.multiplied_scalar(r2));
        (curv1, curv2, vect1, vect2)
    } else {
        (l / e, n / g, *d1u, *d1v)
    };
    let (min_c, max_c, min_v, max_v) = if c1 < c2 {
        (c1, c2, v1, v2)
    } else {
        (c2, c1, v2, v1)
    };
    let dir_min = GpDir::from_vec(&min_v).ok()?;
    let dir_max = GpDir::from_vec(&max_v).ok()?;
    Some((min_c, max_c, dir_min, dir_max))
}

/// Edges of `face` with the face forced `FORWARD` (`BRepLib.cxx:2213-2215`).
fn edges_on_forward_face(face: &Face) -> Vec<Edge> {
    let stored = face
        .0
        .tshape
        .read()
        .expect("poisoned TShape lock")
        .children
        .clone();
    let mut out = Vec::new();
    for child in stored {
        if child.shape_type() != ShapeType::Wire {
            continue;
        }
        out.extend(edges_of_wire(&Wire(child)));
    }
    out
}

/// Pcurves for `ContinuityOfFaces` (`BRepLib.cxx:2201-2248`).
fn face_pcurves(
    edge: &Edge,
    face1: &Face,
    face2: &Face,
) -> Option<(std::sync::Arc<dyn Curve2d>, std::sync::Arc<dyn Curve2d>)> {
    let seam_closed = !is_same(&face1.0, &face2.0)
        && BRepTool::is_closed_edge_face(edge, face1)
        && BRepTool::is_closed_edge_face(edge, face2);
    if seam_closed {
        let edge_in_face1 = edges_on_forward_face(face1)
            .into_iter()
            .find(|e| is_same(&e.0, &edge.0))?;
        let (c1, _, _) = curve_on_surface_oriented(&edge_in_face1, face1, false)?;
        let mut edge_in_face2 = edge_in_face1.clone();
        edge_in_face2.0.reverse();
        let (c2, _, _) = curve_on_surface_oriented(&edge_in_face2, face2, false)?;
        Some((c1, c2))
    } else {
        let (c1, _, _) = curve_on_surface_oriented(edge, face1, false)?;
        let mut edge2 = edge.clone();
        if is_same(&face1.0, &face2.0) {
            edge2.0.reverse();
        }
        let (c2, _, _) = curve_on_surface_oriented(&edge2, face2, false)?;
        Some((c1, c2))
    }
}

/// `BRepLib::ContinuityOfFaces` (`BRepLib.cxx:2197-2418`).
/// `None` is the caught `Standard_Failure`: do not write continuity.
fn continuity_of_faces(
    edge: &Edge,
    face1: &Face,
    face2: &Face,
    angle_tol: f64,
) -> Option<AbsShape> {
    // `BRepLib.cxx:2202`: `isSeam = theFace1.IsEqual(theFace2)`.
    let is_seam = is_equal(&face1.0, &face2.0);
    let (c1, c2) = match face_pcurves(edge, face1, face2) {
        Some(v) => v,
        None => return Some(AbsShape::C0),
    };
    let s1 = unwrap_trimmed(BRepTool::face_surface(face1)?);
    let s2 = unwrap_trimmed(BRepTool::face_surface(face2)?);
    let elementary = is_elementary(s1.as_ref()) && is_elementary(s2.as_ref());
    if is_seam && elementary {
        return Some(AbsShape::Cn);
    }
    let trsf1 = face1.0.location().transformation();
    let trsf2 = face2.0.location().transformation();
    let mut sp1 = FaceProps::new(
        s1.as_ref(),
        c1.as_ref(),
        &trsf1,
        face1.0.orientation().is_reversed(),
    );
    let mut sp2 = FaceProps::new(
        s2.as_ref(),
        c2.as_ref(),
        &trsf2,
        face2.0.orientation().is_reversed(),
    );
    let (mut f, mut l) = BRepTool::edge_parameters(edge);
    let eps = (l - f) / 100.0;
    f += eps;
    l -= eps;
    // `aHC2 = new BRepAdaptor_Curve(anEdgeInFace2, theFace2)`
    // (`cxx:2344-2346`), built lazily once and reused by every sample.
    // `BRepAdaptor_Curve` values live in the face's world frame while
    // `BRep_Tool::Surface` is unlocated (`FaceProps` applies `trsf2` the same
    // way), so the basis surface is transformed here. The adaptor window is
    // the shrunk `[f, l]` of `cxx:2287-2288`; `Initialize` is stateless, so
    // re-initializing per sample replays `cxx:2347` exactly.
    let s2_world: std::sync::Arc<dyn Surface> = std::sync::Arc::from(s2.transformed(&trsf2));
    let cos2 = CurveOnSurface::new(c2.clone(), s2_world, f, l);
    let angle_tol2 = angle_tol * angle_tol;
    let mut cont = if elementary {
        AbsShape::Cn
    } else {
        AbsShape::C2
    };
    for i in 0..=20 {
        if cont <= AbsShape::C0 {
            break;
        }
        let u = f + (l - f) * (i as f64) / 20.0;
        sp1.calculate(u);
        sp2.calculate(u);
        let der1 = sp1.derivative();
        let sq1 = der1.square_magnitude();
        let mut der2 = sp2.derivative();
        let mut sq2 = der2.square_magnitude();
        let mut smooth = der1.cross_square_magnitude(&der2) <= angle_tol2 * sq1 * sq2;
        if smooth {
            let mut n1 = sp1.normal()?;
            if face1.0.orientation().is_reversed() {
                n1.reverse();
            }
            let mut n2 = sp2.normal()?;
            if face2.0.orientation().is_reversed() {
                n2.reverse();
            }
            if n1.dot(&n2) < 0.0 {
                return Some(AbsShape::C0);
            }
        }
        if !smooth {
            // `cxx:2348-2355`: `ext.Perform(aSP1.Value(), u)`; on
            // `IsDone() && IsMin()` re-evaluate the second surface derivative
            // at the projected parameter. `aHC2` is an
            // `Adaptor3d_CurveOnSurface`, whose `GetType()` is the `myType`
            // set by `EvalKPart` (`Adaptor3d_CurveOnSurface.cxx:1374-1377`
            // and `:1552-1720`): `GeomAbs_OtherCurve` unless the pcurve is a
            // plane conic or an iso-arc. The port's `CurveOnSurface` answers
            // `gp_circ()` for the `GeomAbs_Circle` arms but exposes no
            // `gp_line()`, so an `EvalKPart` `GeomAbs_Line` cos (`cxx:1565`,
            // `:1687`, `:1701`) types as `OtherCurve` here; either way the
            // non-iso case runs `Extrema_GLocateExtPC`'s free-form branch
            // (`Extrema_GLocateExtPC.hxx:125-272` -> `Extrema_GenLocateExtPC`,
            // `Extrema_GenLocateExtPC.hxx:108-126`), where `|F| >= 1.e-07`
            // clears `IsDone`.
            let ext = LocateExtPC::new_range(&sp1.value(), &cos2, u, f, l, PCONFUSION);
            if ext.is_done() && ext.is_min() {
                sp2.calculate(ext.point().parameter());
                der2 = sp2.derivative();
                sq2 = der2.square_magnitude();
            }
            smooth = der1.cross_square_magnitude(&der2) <= angle_tol2 * sq1 * sq2;
        }
        let mut cur = if smooth {
            let mut g = AbsShape::G1;
            if (sq1.sqrt() - sq2.sqrt()).abs() < CONFUSION && der1.dot(&der2) > SQUARE_CONFUSION {
                g = AbsShape::C1;
            }
            g
        } else {
            return Some(AbsShape::C0);
        };
        if cont < AbsShape::G2 {
            if cur < cont {
                cont = cur;
            }
            continue;
        }
        let (d1a, c1a, d1b, c1b) = sp1.curvature()?;
        let (d2a, c2a, d2b, c2b) = sp2.curvature()?;
        let dirs1 = [d1a, d1b];
        let lens1 = [c1a, c1b];
        let dirs2 = [d2a, d2b];
        let lens2 = [c2a, c2b];
        for step in 0..2 {
            let other = 1 - step;
            if dirs1[0].xyz().cross_square_magnitude(dirs2[step].xyz()) <= SQUARE_CONFUSION
                && (lens1[0] - lens2[step]).abs() < CONFUSION
                && dirs1[1].xyz().cross_square_magnitude(dirs2[other].xyz()) <= SQUARE_CONFUSION
                && (lens1[1] - lens2[other]).abs() < CONFUSION
            {
                if cur == AbsShape::C1
                    && dirs1[0].dot(&dirs2[step]) > CONFUSION
                    && dirs1[1].dot(&dirs2[other]) > CONFUSION
                {
                    cur = AbsShape::C2;
                } else {
                    cur = AbsShape::G2;
                }
                break;
            }
        }
        if cur < cont {
            cont = cur;
        }
    }
    if elementary && cont == AbsShape::C2 {
        cont = AbsShape::Cn;
    }
    Some(cont)
}

fn sewing_edge(edge: &Edge, face: &Face) -> bool {
    let ori = edge.0.orientation();
    for wire in wires_of_face(face) {
        for cur in edges_of_wire(&wire) {
            if is_same(&edge.0, &cur.0) && ori != cur.0.orientation() {
                return true;
            }
        }
    }
    false
}

/// `BRepLib::EncodeRegularity(E, F1, F2, TolAng)` (`BRepLib.cxx:2570-2589`).
fn encode_one(edge: &Edge, f1: &Face, f2: &Face, tol_ang: f64) {
    if stored_continuity(edge, f1, f2) > AbsShape::C0 as u8 {
        return;
    }
    let Some(cont) = continuity_of_faces(edge, f1, f2, tol_ang) else {
        return;
    };
    write_continuity(edge, f1, f2, cont);
}

/// `BRepLib::EncodeRegularity(S, TolAng)` (`BRepLib.cxx:2428-2536`).
/// The list overload (`cxx:2545`) is not used by `EncodeRegul`.
pub fn encode_regularity(shape: &TopoShape, tol_ang: f64) {
    if shape.shape_type() == ShapeType::Compound || shape.shape_type() == ShapeType::CompSolid {
        let children = shape.tshape.read().expect("poisoned TShape lock").children.clone();
        for child in &children {
            encode_regularity(child, tol_ang);
        }
        return;
    }
    let mut uses: HashMap<usize, (Edge, Vec<Face>)> = HashMap::new();
    for face in faces_of(shape) {
        for wire in wires_of_face(&face) {
            for edge in edges_of_wire(&wire) {
                let key = GeometryRegistry::shape_key(&edge.0);
                let entry = uses.entry(key).or_insert_with(|| (edge.clone(), Vec::new()));
                if !entry
                    .1
                    .iter()
                    .any(|f| is_same(&f.0, &face.0))
                {
                    entry.1.push(face.clone());
                }
            }
        }
    }
    for (_key, (edge, faces)) in &uses {
        if faces.len() >= 2 {
            encode_one(edge, &faces[0], &faces[1], tol_ang);
        } else if faces.len() == 1 && sewing_edge(edge, &faces[0]) {
            encode_one(edge, &faces[0], &faces[0], tol_ang);
        }
    }
}
