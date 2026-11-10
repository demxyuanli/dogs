use crate::shape::TopoShape;
use occt_core::gp::GpDir;

use super::super::data_structure::DataStructure;
use super::super::h_data_structure::HDataStructure;
use super::super::kind::Config;


/// `FDS_aresamdom(BDS, ES, F1, F2)` (`cxx:282-321`): true when the faces
/// `F1`, `F2`, issued from the same origin shape and connected by the edge
/// `ES`, are tangent.
///
/// 1/ the SameDomain coding of the DS; 2/ when the DS carries none — the
/// SameDomain information is not looked for between the faces of one origin
/// shape — the geometric test on the two plane normals at the UV points of
/// `ES` (`cxx:295-318`).
pub fn fds_aresamdom(
    bds: &DataStructure,
    es: &TopoShape,
    f1: &TopoShape,
    f2: &TopoShape,
) -> bool {
    for ss in bds.shape_same_domain(f1) {
        if crate::topo_tools::is_same(ss, f2) {
            return true;
        }
    }

    let face1 = crate::shape::Face(f1.clone());
    let face2 = crate::shape::Face(f2.clone());
    let (Some(su1), Some(su2)) = (
        crate::brep_tool::BRepTool::face_surface(&face1),
        crate::brep_tool::BRepTool::face_surface(&face2),
    ) else {
        return false;
    };
    let su1 = crate::topopebrep::tool::basis_surface(&su1);
    let su2 = crate::topopebrep::tool::basis_surface(&su2);
    let st1 = crate::geom_bnd_lib_surface3d::surface_kind(&*su1);
    let st2 = crate::geom_bnd_lib_surface3d::surface_kind(&*su2);
    let plpl = st1 == crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType::Plane
        && st2 == crate::geom_bnd_lib_surface3d::GeomAbsSurfaceType::Plane;
    if !plpl {
        return false;
    }

    let edge = crate::shape::Edge(es.clone());
    let (p2d1, p2d2) = crate::brep_tool::BRepTool::uv_points(&edge, &face1);
    let d1 = crate::topopebrep::tool::fun_tool_ng_s(&p2d1, &*su1);
    let d2 = crate::topopebrep::tool::fun_tool_ng_s(&p2d2, &*su2);
    let tola = occt_core::precision::ANGULAR;
    let dot = d1.dot(&d2);
    (1.0 - dot.abs()).abs() < tola
}

/// `FDS_aresamdom(BDS, SI, isb1, isb2)` (`cxx:330-357`): the shape-index
/// overload, which filters the three shapes by type first.
pub fn fds_aresamdom_from_indices(bds: &DataStructure, si: i32, isb1: i32, isb2: i32) -> bool {
    let Some(es) = bds.shape(si, true) else {
        return false;
    };
    if es.shape_type() != crate::abs::ShapeType::Edge {
        return false;
    }
    let Some(f1) = bds.shape(isb1, true) else {
        return false;
    };
    if f1.shape_type() != crate::abs::ShapeType::Face {
        return false;
    }
    let Some(f2) = bds.shape(isb2, true) else {
        return false;
    };
    if f2.shape_type() != crate::abs::ShapeType::Face {
        return false;
    }
    fds_aresamdom(bds, &es, &f1, &f2)
}

/// `FDS_EdgeIsConnexToSameDomainFaces(E, HDS)` (`cxx:357-391`): at least two
/// faces connected through `E`, one pair of which `FDS_aresamdom` reports as
/// same-domain.
///
/// The inner loop starts at the outer cursor (`cxx:376`), so the first pair
/// tested is `(f1, f1)`; `samdom` is overwritten on every inner iteration and
/// the loops stop on the first `true` (`cxx:379-388`). `HDS` is unused by
/// OCCT (`cxx:359`) beyond reaching the DS.
pub fn fds_edge_is_connex_to_same_domain_faces(
    connex: &super::super::connex::Connexity,
    e: &TopoShape,
    hds: &HDataStructure,
) -> bool {
    let bds = hds.ds();
    let lf = connex.edge_connexity_same_shape(e, hds);
    if lf.len() < 2 {
        return false;
    }
    for i1 in 0..lf.len() {
        let f1 = &lf[i1];
        for j in i1..lf.len() {
            let f2 = &lf[j];
            if fds_aresamdom(bds, e, f1, f2) {
                return true;
            }
        }
    }
    false
}

/// `FDS_HasSameDomain3d(BDS, E, PLSD)` (`cxx:498-537`): `E` has a
/// same-domain partner whose `SameDomainOri` is `UNSHGEOMETRY` (3-D geometry
/// not shared).
///
/// The C++ `PLSD` out-parameter collects every such partner; `None` is a null
/// `PLSD`, in which case OCCT `break`s on the first hit (`cxx:529-532`). The
/// port always collects, and the `break`/`continue` split is reproduced by the
/// `collect_all` flag so a caller emitting only the boolean reads the same
/// early exit.
pub fn fds_has_same_domain3d(
    bds: &DataStructure,
    e: &TopoShape,
    collect_all: bool,
) -> (bool, Vec<TopoShape>) {
    let lssd = bds.shape_same_domain(e);
    let mut out = Vec::new();
    if lssd.is_empty() {
        return (false, out);
    }
    let mut hsd3d = false;
    for esd in lssd.iter() {
        let c = bds.same_domain_ori_of(esd);
        if c == Config::UnshGeometry {
            hsd3d = true;
            out.push(esd.clone());
            if !collect_all {
                break;
            }
        }
    }
    (hsd3d, out)
}

/// `FDS_Config3d(E1, E2, c)` (`cxx:539-575`): `SAMEORIENTED` /
/// `DIFFORIENTED` from the sign of the dot product of the two 3-D tangents at
/// the midpoint of `E1` and its projection on `E2`. `None` is `cxx:563-566`
/// (`!ok1 || !ok2`).
pub fn fds_config3d(e1: &TopoShape, e2: &TopoShape) -> Option<Config> {
    let (pe1, par_e1) = crate::topopebrep::tool::fun_tool_find_p_in_e(e1)?;
    let ve1 = crate::topopebrep::tool::tg_geom_e(par_e1, e1)?;
    let (par_e2, _d_e2) = crate::topopebrep::tool::fun_tool_proj_pon_e(&pe1, e2)?;
    let ve2 = crate::topopebrep::tool::tg_geom_e(par_e2, e2)?;
    let de1 = GpDir::from_vec(&ve1).ok()?;
    let de2 = GpDir::from_vec(&ve2).ok()?;
    let dot = de1.dot(&de2);
    Some(if dot > 0.0 {
        Config::SameOriented
    } else {
        Config::DiffOriented
    })
}

/// `FDS_HasSameDomain2d(BDS, E, PLSD)` (`cxx:577-615`): the 2-D counterpart,
/// keeping partners whose `SameDomainOri` is `SAMEORIENTED` or `DIFFORIENTED`.
pub fn fds_has_same_domain2d(
    bds: &DataStructure,
    e: &TopoShape,
    collect_all: bool,
) -> (bool, Vec<TopoShape>) {
    let lssd = bds.shape_same_domain(e);
    let mut out = Vec::new();
    if lssd.is_empty() {
        return (false, out);
    }
    let mut hsd2d = false;
    for esd in lssd.iter() {
        let c = bds.same_domain_ori_of(esd);
        if c == Config::SameOriented || c == Config::DiffOriented {
            hsd2d = true;
            out.push(esd.clone());
            if !collect_all {
                break;
            }
        }
    }
    (hsd2d, out)
}

/// `FUN_ds_ONesd(BDS, IE, EspON, IEsd)` (`cxx:3308-3345`): the same-domain
/// partner of `E = BDS.Shape(IE)` that carries the ON part `EspON`.
///
/// The probe point is taken at `x = 0.456789` of `EspON`'s range
/// (`cxx:3316-3319`); a partner is kept when the probe projects on it within
/// `1.e3 * tolerance` (`cxx:3332-3339`). The answer is the DS index of that
/// partner, written through `iesd`.
pub fn fun_ds_on_esd(
    bds: &DataStructure,
    ie: i32,
    esp_on: &TopoShape,
    iesd: &mut i32,
) -> bool {
    // `cxx:3313`: `BDS.Shape(IE)` — the default `findWithoutKeep = true`.
    let Some(e) = bds.shape(ie, true) else {
        return false;
    };
    let (f, l) = crate::topopebrep::tool::fun_tool_bounds(esp_on);
    let x = 0.456789;
    let par = (1.0 - x) * f + x * l;
    let Some(p3d) = crate::topopebrep::tool::fun_tool_value(par, esp_on) else {
        return false;
    };

    let lesd = bds.shape_same_domain(&e).to_vec();
    for esd in lesd {
        let Some((_parp, d)) = crate::topopebrep::tool::fun_tool_proj_pon_e(&p3d, &esd) else {
            continue;
        };
        let tolesd = crate::brep_tool::BRepTool::edge_tolerance(&crate::shape::Edge(esd.clone()));
        if !(d < tolesd * 1.0e3) {
            continue;
        }
        // `cxx:3341`: `IEsd = BDS.Shape(esd)` — again `findWithoutKeep = true`.
        *iesd = bds.shape_index(&esd, true);
        return true;
    }
    false
}

