//! `IsInput` / `IsG1` / `SearchFaceOnV` — the three file-static helpers that
//! `ChFi3d_Builder::SearchFace` walks a vertex of the guideline with.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_2.cxx:399-589`
//! (OCCT 8.0.0).
//!
//! `IsInput(Vec, Ve, Fa)` (`cxx:399-471`) decides whether the direction `Vec`
//! leaves the face `Fa` through the vertex `Ve`: it collects the two edges of
//! `Fa` that start (or end) at `Ve`, builds the plane of their tangents and
//! checks that the angle of `Vec`, projected into that plane, falls between
//! them.
//!
//! `IsG1` (`cxx:475-521`) answers the face of `TheMap(E)` other than `FRef`
//! that is tangent to it through `E`, with the cutting-edge (`FVoi == FRef`)
//! fallback of `cxx:496-520`.
//!
//! `SearchFaceOnV` (`cxx:532-589`) is the two-face form used when the common
//! point carries its own vector.
//!
//! Divergences that are not heuristics:
//! - `TopExp_Explorer(Fa, TopAbs_WIRE)` / `BRepTools_WireExplorer` are
//!   [`topo_tools::wires_of_face`] / [`topo_tools::edges_of_wire`]; the port
//!   visits the wires in the same order and resets the counter per wire, which
//!   is what `cxx:414` (`for (Trouve = 0, WireExp.Init(W); ...)`) does,
//! - `BRepAdaptor_Curve Cb(E); Cb.D1(Parameter(Ve, E), P, V)` is
//!   [`CurveAdaptor::of_edge`] + [`CurveAdaptor::d1`], `BRep_Tool::Parameter`
//!   is [`BRepTool::parameter_of_vertex`]; an edge without a 3D curve raises
//!   `Standard_NullObject` in OCCT and contributes no direction here,
//! - `gp_Trsf::SetTransformation(gp_Ax3)` / `Transforms(gp_XYZ)` / `Invert()` are
//!   the [`GpTrsf`] methods of the same names, so the projection of the
//!   direction into the tangent plane is OCCT's own sequence and not a rewritten
//!   projection,
//! - the `gp_Dir(gp_Vec)` conversions OCCT makes implicitly inside
//!   `gp_Ax3(Point, Normal, Vec3d[0])` are [`GpDir::new`], whose
//!   `Standard_ConstructionError` on a null vector answers `false` here — the
//!   raise OCCT would take at that point,
//! - `occ::down_cast<TopoDS::Edge>` / `TopoDS::Face` are [`Edge::wrap`] /
//!   [`Face::wrap`]; the `IsSame` identity is [`topo_tools::is_same`].

use occt_core::gp::{GpAx3, GpDir, GpPnt, GpTrsf, GpVec};
use occt_core::kernel::geomabs::Shape;
use occt_core::precision::CONFUSION;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::chfi_ds::{ChFiDSMap, CommonPoint};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools;
use crate::topopebrep::tool::adaptor::CurveAdaptor;

/// `IsInput(Vec, Ve, Fa)` (`cxx:399-471`): `Vec` points out of `Fa` at `Ve`.
///
/// `fa` is an `Option` because a null `TopoDS_Face` has no counterpart in the
/// port; it answers `false` exactly like an empty explorer.
pub fn is_input(vec: &GpVec, ve: &TopoShape, fa: Option<&Face>) -> bool {
    // `gp_Vec Vec3d[2]; gp_Pnt Point;` (`cxx:407-408`).
    let mut vec3d = [GpVec::new(0.0, 0.0, 0.0), GpVec::new(0.0, 0.0, 0.0)];
    let mut point: Option<GpPnt> = None;
    let mut trouve = 0usize;

    // `for (; (FaceExp.More() && (Trouve < 2)); FaceExp.Next())` (`cxx:411`).
    if let Some(fa) = fa {
        for wire in topo_tools::wires_of_face(fa) {
            // `for (Trouve = 0, WireExp.Init(W); ...)` (`cxx:414`): the counter
            // restarts at every wire.
            trouve = 0;
            for e in topo_tools::edges_of_wire(&wire) {
                if trouve >= 2 {
                    break;
                }
                // `TopExp::Vertices(E, Vf, Vl)` (`cxx:417`).
                let (vf, vl) = topo_tools::edge_vertices(&e);
                let at_first = vf
                    .as_ref()
                    .is_some_and(|v| topo_tools::is_same(&v.0, ve));
                // `else if (Vl.IsSame(Ve))` (`cxx:424`).
                let at_last = !at_first
                    && vl
                        .as_ref()
                        .is_some_and(|v| topo_tools::is_same(&v.0, ve));
                if !at_first && !at_last {
                    continue;
                }
                // `BRepAdaptor_Curve Cb(E); Cb.D1(BRep_Tool::Parameter(Ve, E),
                // Point, Vec3d[Trouve])` (`cxx:420-428`).
                let cb = CurveAdaptor::of_edge(&e);
                let par = BRepTool::parameter_of_vertex(&e.0, ve);
                let Some((p, d1)) = cb.d1(par) else {
                    // `Standard_NullObject` in OCCT; no direction to add.
                    continue;
                };
                point = Some(p);
                vec3d[trouve] = if at_last {
                    // `Vec3d[Trouve].Reverse()` (`cxx:428`).
                    d1.reversed()
                } else {
                    d1
                };
                trouve += 1;
            }
            if trouve >= 2 {
                break;
            }
        }
    }

    // `if (Trouve < 2) return false;` (`cxx:433-436`).
    if trouve < 2 {
        return false;
    }
    let Some(point) = point else {
        return false;
    };

    // `Normal = Vec3d[0] ^ Vec3d[1]` (`cxx:438-439`).
    let normal = vec3d[0].crossed(&vec3d[1]);
    // `if (Normal.SquareMagnitude() < Precision::Confusion())` (`cxx:440`).
    if normal.square_magnitude() < CONFUSION {
        // `return Vec.IsParallel(Vec3d[0], Precision::Confusion())` (`cxx:442`).
        return vec.is_parallel_ang(&vec3d[0], CONFUSION);
    }

    // `amax = Vec3d[1].AngleWithRef(Vec3d[0], Normal)` (`cxx:446-455`).
    let mut amax = vec3d[1].angle_with_ref(&vec3d[0], &normal);
    let mut amin = 0.0;
    if amax < 0.0 {
        amin = amax;
        amax = 0.0;
    }

    // `gp_Ax3 Axe(Point, Normal, Vec3d[0])` (`cxx:458`).
    let (Ok(z_dir), Ok(x_dir)) = (
        GpDir::new(normal.x(), normal.y(), normal.z()),
        GpDir::new(vec3d[0].x(), vec3d[0].y(), vec3d[0].z()),
    ) else {
        // `Standard_ConstructionError` of the implicit `gp_Dir(gp_Vec)`.
        return false;
    };
    let Ok(axe) = GpAx3::new(point, z_dir, &x_dir) else {
        return false;
    };

    // `gp_Trsf Transf; Transf.SetTransformation(Axe); gp_XYZ coord = Vec.XYZ();
    // Transf.Transforms(coord); coord.SetZ(0); Transf.Invert();
    // Transf.Transforms(coord); gp_Vec theProj(coord);` (`cxx:459-466`).
    let mut transf = GpTrsf::identity();
    transf.set_transformation(&axe);
    let mut coord = *vec.xyz();
    transf.transforms_xyz(&mut coord);
    coord.set_z(0.0);
    if transf.invert().is_err() {
        return false;
    }
    transf.transforms_xyz(&mut coord);
    let the_proj = GpVec::new(coord.x(), coord.y(), coord.z());

    // `Angle = theProj.AngleWithRef(Vec3d[0], Normal); return ((Angle >= amin)
    // && (Angle <= amax));` (`cxx:469-470`).
    let angle = the_proj.angle_with_ref(&vec3d[0], &normal);
    angle >= amin && angle <= amax
}

/// `IsG1(TheMap, E, FRef, FVoi)` (`cxx:475-521`): the face of `TheMap(E)`
/// other than `FRef` that is `G1` to it through `E`, or `FRef` itself when `E`
/// is a cutting edge.
///
/// `f_voi` is the out-parameter; it is left untouched when the answer is
/// `false`, as in OCCT.
pub fn is_g1(
    ef_map: &ChFiDSMap,
    e: &Edge,
    f_ref: Option<&Face>,
    f_voi: &mut Option<Face>,
) -> bool {
    let Some(f_ref) = f_ref else {
        // `TopoDS::Face(It.Value()).IsSame(FRef)` is always false against a null
        // `FRef`, which leaves `ChFi3d::IsTangentFaces` without a face to test.
        return false;
    };

    // `for (It.Initialize(TheMap(E)); It.More(); It.Next())` (`cxx:482-495`).
    for value in ef_map.get_from_key(&e.0) {
        if topo_tools::is_same(value, &f_ref.0) {
            continue;
        }
        let Some(candidate) = Face::wrap(value.clone()) else {
            continue;
        };
        *f_voi = Some(candidate.clone());
        // `if (ChFi3d::IsTangentFaces(E, FRef, FVoi))` (`cxx:489`); the failure
        // OCCT raises is the `Err` of the port, answered `false`.
        if crate::chfi3d::is_tangent_faces(e, f_ref, &candidate, Shape::G1)
            .unwrap_or(false)
        {
            return true;
        }
    }

    // Cutting edge: `E` appears twice in `FRef` with opposite orientations
    // (`cxx:498-519`).
    let mut or_set = false;
    let mut orient = Orientation::Forward;
    for ed in topo_tools::edges_of(&f_ref.0) {
        if !topo_tools::is_same(&ed.0, &e.0) {
            continue;
        }
        if !or_set {
            orient = ed.0.orientation();
            or_set = true;
        } else if ed.0.orientation() == orient.reversed() {
            *f_voi = Some(f_ref.clone());
            // `return ChFi3d::IsTangentFaces(E, FRef, FRef)` (`cxx:517`).
            return crate::chfi3d::is_tangent_faces(e, f_ref, f_ref, Shape::G1)
                .unwrap_or(false);
        }
    }
    false
}

/// `SearchFaceOnV(Pc, FRef, VEMap, EFMap, F1, F2)` (`cxx:532-589`): up to two
/// faces around the vertex of `Pc` that the vector of `Pc` leaves.
pub fn search_face_on_v(
    pc: &CommonPoint,
    f_ref: Option<&Face>,
    ve_map: &ChFiDSMap,
    ef_map: &ChFiDSMap,
    f1: &mut Option<Face>,
    f2: &mut Option<Face>,
) -> usize {
    // `bool FindFace = IsInput(Pc.Vector(), Pc.Vertex(), FRef); if (FindFace)
    // FindFace = IsInput(Pc.Vector().Reversed(), Pc.Vertex(), FRef);`
    // (`cxx:540-544`).
    let (Some(vector), Some(vertex)) = (pc.vector(), pc.vertex()) else {
        // OCCT dereferences the two handles of the common point; the port sees
        // the absent value instead.
        return 0;
    };
    let mut find_face = is_input(vector, vertex, f_ref);
    if find_face {
        find_face = is_input(&vector.reversed(), vertex, f_ref);
    }
    // `if (FindFace) { F1 = FRef; return 1; }` (`cxx:546-550`).
    if find_face {
        *f1 = f_ref.cloned();
        return 1;
    }

    // `for (ItE.Initialize(VEMap(Pc.Vertex())); ItE.More() && (Num < 2);
    // ItE.Next())` (`cxx:557-587`).
    let mut num = 0usize;
    for value in ve_map.get_from_key(vertex) {
        if num >= 2 {
            break;
        }
        let Some(e) = Edge::wrap(value.clone()) else {
            continue;
        };
        // `for (ItF.Initialize(EFMap(E)), Trouve = false; ...)`.
        let mut trouve = f_ref.is_some_and(|f_ref| {
            ef_map
                .get_from_key(&e.0)
                .iter()
                .any(|f| topo_tools::is_same(f, &f_ref.0))
        });
        if trouve {
            let mut f_voi = None;
            trouve = is_g1(ef_map, &e, f_ref, &mut f_voi);
            if trouve {
                // `Trouve = IsInput(Pc.Vector(), Pc.Vertex(), FVoi)`.
                trouve = match f_voi {
                    Some(ref f_voi) => is_input(vector, vertex, Some(f_voi)),
                    None => false,
                };
                if trouve {
                    if num == 0 {
                        *f1 = f_voi;
                    } else {
                        *f2 = f_voi;
                    }
                    num += 1;
                }
            }
        }
    }
    num
}
