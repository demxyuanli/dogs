//! `TopOpeBRepDS_FaceInterferenceTool` — the complex face/curve transition of
//! a face interference.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_FaceInterferenceTool.hxx:34-80`
//! and `TopOpeBRepDS_FaceInterferenceTool.cxx:1-453` (OCCT 8.0.0).
//!
//! Translation notes:
//! - `STATIC_TOREVERSE` (`cxx:35`, `xpu150498`) is a `thread_local` `bool` in
//!   OCCT. It is set by `Init` (`cxx:192`, `:207`) and read by `Add`
//!   (`cxx:310`) and `Transition` (`cxx:417`). The port carries it as the
//!   `to_reverse` field: OCCT's own use is single-instance (the tool is a
//!   member of `TopOpeBRepDS_Reducer`, driven one interfered edge at a time),
//!   so the field reproduces the same lifetime without shared mutable state.
//! - `myPBDS` is `TopOpeBRepDS_PDataStructure`, a raw `DataStructure*`
//!   (`TopOpeBRepDS_PDataStructure.hxx:21`); here it is a shared borrow.
//! - `FUN_Parameters` (`cxx:46-68`) and the two `FUN_ComputeGeomData`
//!   overloads (`cxx:72-76`, `:86-161`) are `Standard_EXPORT` free functions of
//!   this translation unit and are re-declared by `TopOpeBRepDS_FIR.cxx`; they
//!   are ported here as [`fun_parameters`] / [`fun_compute_geom_data`] /
//!   [`fun_compute_geom_data_full`].
//! - `TopTrans_SurfaceTransition` is `crate::top_trans_surface_transition`.
//! - The `Standard_ProgramError` throws become `Err(String)`.

use occt_core::gp::{GpDir, GpPnt, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::extrema_surf::{ExtPs, ExtremaExtAlgo, ExtremaExtFlag};
use occt_geom::geom_lprop::SlProps;

use crate::abs::{Orientation, State};
use crate::brep_tool::BRepTool;
use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};
use crate::shape::{Edge, Face, TopoShape};
use crate::top_trans_surface_transition::SurfaceTransition;
use crate::topopebrep::ds::data_structure::DataStructure;
use crate::topopebrep::ds::export::{fds_config3d, fds_data, fds_has_same_domain3d};
use crate::topopebrep::ds::interference::Interference;
use crate::topopebrep::ds::kind::{Config, Kind};
use crate::topopebrep::tool::{
    basis_surface, fun_tool_find_p_in_e, fun_tool_ng_geom_f, fun_tool_orient_ein_f_forward,
    fun_tool_par_on_e_f, fun_tool_plane, fun_tool_proj_pon_f,
};

/// `FUN_Parameters(Pnt, F, u, v)` (`cxx:46-68`): the 2-D parameters of the
/// projection of `Pnt` on `F`'s surface, accepted when the squared distance is
/// below `tolF^2 * 1.e6` (`cxx:66-68`, the `// NYINYI` factor of `cxx:67`).
///
/// `Extrema_ExtPS(Pnt, Surf, fu, lu, fv, lv, uvtol, uvtol, Extrema_ExtFlag_MIN)`
/// is [`ExtPs`] with the flag set to `Min` and the adaptor's own bounds.
/// `!IsDone()` and `NbExt() == 0` are the two `false` arms (`cxx:55-63`).
pub fn fun_parameters(p: &GpPnt, f: &Face) -> Option<(f64, f64)> {
    let s = BRepTool::face_surface_world(f)?;
    let (fu, lu) = s.u_range();
    let (fv, lv) = s.v_range();
    let uvtol = BRepTool::face_tolerance(f);
    let mut ex = ExtPs::new();
    ex.set_flag(ExtremaExtFlag::Min);
    ex.set_algo(ExtremaExtAlgo::Grad);
    ex.initialize(s.as_ref(), fu, lu, fv, lv, uvtol, uvtol);
    ex.perform(p);
    if !ex.is_done() || ex.nb_ext() == 0 {
        return None;
    }
    let (u, v, _) = ex.point(1);

    // xpu281098 : CTS21216 (FIR, f4,e7on)
    let d2 = ex.square_distance(1);
    let tol_f = BRepTool::face_tolerance(f);
    if d2 < tol_f * tol_f * 1.0e6 {
        Some((u, v))
    } else {
        None
    }
}

/// `FUN_ComputeGeomData(F, uv, Norm)` (`cxx:72-76`).
pub fn fun_compute_geom_data(f: &Face, uv: (f64, f64)) -> Option<GpDir> {
    fun_tool_ng_geom_f(&occt_core::gp::GpPnt2d::new(uv.0, uv.1), f)
}

/// `FUN_sphere(F)` (`cxx:78-84`): the basis surface's adaptor type is
/// `GeomAbs_Sphere`.
fn fun_sphere(base: &dyn occt_geom::surface::Surface) -> bool {
    surface_kind(base) == GeomAbsSurfaceType::Sphere
}

/// `FUN_ComputeGeomData(F, uv, Norm, D1, D2, Cur1, Cur2)` (`cxx:86-161`).
///
/// The `!curdef` throw (`cxx:104-107`) and the `Standard_Failure` of the
/// umbilic non-plane non-sphere arm (`cxx:129-132`) are the two `Err` cases.
pub fn fun_compute_geom_data_full(
    f: &Face,
    uv: (f64, f64),
) -> Result<(GpDir, GpDir, GpDir, f64, f64), String> {
    let uu = uv.0;
    let vv = uv.1;
    let Some(s) = BRepTool::face_surface_world(f) else {
        return Err("TopOpeBRepDS_FaceInterferenceTool::Init".to_string());
    };
    let base = basis_surface(&s);

    let sphere = fun_sphere(base.as_ref());
    let mut plane = fun_tool_plane(f);

    // Getting the principal directions, the normal and the curvatures.
    let mut props = SlProps::new(s.as_ref(), uu, vv, 2, CONFUSION)?;
    let curv_def = props.is_curvature_defined(s.as_ref());
    if !curv_def {
        return Err("TopOpeBRepDS_FaceInterferenceTool::Init".to_string());
    }
    let umbilic = props.is_umbilic(s.as_ref())?;

    let (cur1, cur2, d1, d2, norm);
    if umbilic {
        let m = props.mean_curvature(s.as_ref())?;
        cur1 = m;
        cur2 = m;

        // xpu030998 : cto901A3
        let toll = 1.0e-8;
        let ooplane = cur1.abs() < toll && cur2.abs() < toll;
        plane = plane || ooplane;

        let n = if plane {
            fun_tool_ng_geom_f(&occt_core::gp::GpPnt2d::new(uu, vv), f)
                .ok_or_else(|| "FUN_ComputeGeomData".to_string())?
        } else if sphere {
            let sph = base
                .gp_sphere()
                .ok_or_else(|| "FUN_ComputeGeomData".to_string())?;
            let center = sph.location();
            let value = s.d0(uu, vv);
            GpDir::from_vec(&GpVec::from_pnts(&center, &value))
                .map_err(|e| format!("FUN_ComputeGeomData: {e}"))?
        } else {
            return Err("FUN_ComputeGeomData".to_string());
        };

        let d1v = n;
        let (x, y, z) = (d1v.x(), d1v.y(), d1v.z());
        let tol = CONFUSION;
        let nullx = x.abs() < tol;
        let nully = y.abs() < tol;
        let nullz = z.abs() < tol;
        let d2v = if nullx && nully {
            GpDir::new(1.0, 0.0, 0.0).map_err(|e| format!("gp_Dir: {e}"))?
        } else if nullx && nullz {
            GpDir::new(1.0, 0.0, 0.0).map_err(|e| format!("gp_Dir: {e}"))?
        } else if nully && nullz {
            GpDir::new(0.0, 1.0, 0.0).map_err(|e| format!("gp_Dir: {e}"))?
        } else {
            GpDir::new(y * z, x * z, -2.0 * x * y).map_err(|e| format!("gp_Dir: {e}"))?
        };
        norm = n;
        d1 = d1v;
        d2 = d2v;
    } else {
        cur1 = props.max_curvature(s.as_ref())?;
        cur2 = props.min_curvature(s.as_ref())?;
        let (max_d, min_d) = props.curvature_directions(s.as_ref())?;
        d1 = max_d;
        d2 = min_d;
        norm = fun_tool_ng_geom_f(&occt_core::gp::GpPnt2d::new(uu, vv), f)
            .ok_or_else(|| "FUN_ComputeGeomData".to_string())?;
    }

    Ok((norm, d1, d2, cur1, cur2))
}

/// `TopOpeBRepDS_FaceInterferenceTool` (`hxx:34-80`).
pub struct FaceInterferenceTool<'a> {
    /// `myPBDS`.
    p_bds: &'a DataStructure,
    /// `myrefdef`.
    ref_def: bool,
    /// `myFaceOrientation`.
    face_orientation: Orientation,
    /// `myFaceOriented`.
    face_oriented: i32,
    /// `myTool`.
    tool: SurfaceTransition,
    /// `myEdge`.
    edge: TopoShape,
    /// `isLine`.
    is_line: bool,
    /// `myPntOnEd`.
    pnt_on_ed: GpPnt,
    /// `myParOnEd`.
    par_on_ed: f64,
    /// `myOnEdDef`.
    on_ed_def: bool,
    /// `myTole`.
    tole: f64,
    /// `STATIC_TOREVERSE` (`cxx:35`, `xpu150498`).
    to_reverse: bool,
}

impl<'a> FaceInterferenceTool<'a> {
    /// `TopOpeBRepDS_FaceInterferenceTool(PBDS)` (`cxx:164-171`): `myrefdef =
    /// myOnEdDef = false`.
    pub fn new(p_bds: &'a DataStructure) -> Self {
        Self {
            p_bds: p_bds,
            ref_def: false,
            face_orientation: Orientation::Forward,
            face_oriented: 0,
            tool: SurfaceTransition::new(),
            edge: TopoShape::new(crate::abs::ShapeType::Shape),
            is_line: false,
            pnt_on_ed: GpPnt::zero(),
            par_on_ed: 0.0,
            on_ed_def: false,
            tole: 0.0,
            to_reverse: false,
        }
    }

    /// `Init(FFI, EE, EEisnew, Iin)` (`cxx:175-273`).
    ///
    /// `I.IsNull()` is `!i.is_shape_shape()` (`cxx:181-186`), since
    /// `TopOpeBRepDS_ShapeShapeInterference` is the type
    /// `occ::down_cast` keeps.
    pub fn init(
        &mut self,
        ffi: &TopoShape,
        ee: &TopoShape,
        ee_is_new: bool,
        i: &Interference,
    ) -> Result<(), String> {
        if !i.is_shape_shape() {
            return Ok(());
        }
        let fi = Face(ffi.clone());
        let e = Edge(ee.clone());

        // xpu150498
        self.to_reverse = false;
        if ee_is_new {
            let g = i.geometry();
            let eg = match self.p_bds.shape(g, true) {
                Some(s) => Edge(s),
                None => return Err("TopOpeBRepDS_FaceInterferenceTool".to_string()),
            };
            let cf = fds_config3d(&e.0, &eg.0)
                .ok_or_else(|| "TopOpeBRepDS_FaceInterferenceTool".to_string())?;
            if cf == Config::DiffOriented {
                self.to_reverse = true;
            }
        } // xpu150498

        self.face_orientation = fi.0.orientation();
        self.face_oriented = i.support();

        self.edge = e.0.clone();
        // Get a middle point on `E`; the geometric data is described locally
        // around this point.
        let _o_e_in_fi = fun_tool_orient_ein_f_forward(&e.0, &fi.0);
        // isLine = FUN_tool_line(E);
        self.is_line = false;

        if !self.on_ed_def {
            let (p, par) = fun_tool_find_p_in_e(&e.0)
                .ok_or_else(|| "TopOpeBRepDS_FaceInterferenceTool".to_string())?;
            self.pnt_on_ed = p;
            self.par_on_ed = par;
        }

        self.tole = occt_core::precision::ANGULAR;
        let uv = match fun_tool_par_on_e_f(&e.0, self.par_on_ed, &fi.0) {
            Some(uv) => uv,
            None => {
                let (u, v, _d) = fun_tool_proj_pon_f(&self.pnt_on_ed, &fi)
                    .ok_or_else(|| "TopOpeBRepDS_FaceInterferenceTool".to_string())?;
                occt_core::gp::GpPnt2d::new(u, v)
            }
        };

        let tmp = crate::topopebrep::tool::tg_geom_e(self.par_on_ed, &e.0)
            .ok_or_else(|| "TopOpeBRepDS_FaceInterferenceTool".to_string())?;
        let tgt = GpDir::from_vec(&tmp).map_err(|e| format!("gp_Dir: {e}"))?;
        let uv = (uv.x(), uv.y());
        if self.is_line {
            let norm = fun_compute_geom_data(&fi, uv)
                .ok_or_else(|| "FUN_ComputeGeomData".to_string())?;
            self.tool.reset_line(&tgt, &norm);
        } else {
            let (norm, d1, d2, cur1, cur2) = fun_compute_geom_data_full(&fi, uv)?;
            self.tool.reset(&tgt, &norm, &d1, &d2, cur1, cur2);
        }
        self.ref_def = true;
        Ok(())
    }

    /// `Add(FFI, FFT, EE, EEisnew, Iin)` (`cxx:277-383`).
    pub fn add(
        &mut self,
        ffi: &TopoShape,
        fft: &TopoShape,
        ee: &TopoShape,
        ee_is_new: bool,
        i: &Interference,
    ) -> Result<(), String> {
        if !i.is_shape_shape() {
            return Ok(());
        }
        let fi = Face(ffi.clone());
        let ft = Face(fft.clone());
        let e = Edge(ee.clone());
        if !crate::topo_tools::is_same(&e.0, &self.edge) {
            return Err("TopOpeBRepDS_FaceInterferenceTool".to_string());
        }

        if !self.ref_def {
            // First interference on the oriented face: Init.
            return self.init(ffi, ee, ee_is_new, i);
        }

        let (_gt, g, _st, _s) = fds_data(i);
        let eg = match self.p_bds.shape(g, true) {
            Some(s) => Edge(s),
            None => return Err("TopOpeBRepDS_FaceInterferenceTool".to_string()),
        };
        let _hsd = fds_has_same_domain3d(self.p_bds, &eg.0, true);
        let same = !self.to_reverse; // xpu150498

        let mut ori_loc = i.transition().orientation(State::In, crate::abs::ShapeType::Face);
        // xpu150498 : CTS20205 : sp(e5) = sp(e4 of rank=1) and c3d(e5) c3d(e4)
        // are diff oriented. As transitions on face<iFI> are given relative to
        // the geometry of e5, we have to complement them. cto 016 E1
        let rev = !same && (ori_loc == Orientation::Forward || ori_loc == Orientation::Reversed);
        if rev {
            ori_loc = ori_loc.complement(); // xpu150498
        }

        let ori_eft = fun_tool_orient_ein_f_forward(&eg.0, &ft.0);
        let ori_efi = fun_tool_orient_ein_f_forward(&eg.0, &fi.0);
        let (ori_tan, egofft) = match ori_eft {
            Some(o) => {
                let mut oritan = o;
                if ee_is_new && !same {
                    oritan = o.complement();
                }
                (oritan, true)
            }
            None => match ori_efi {
                Some(o) => {
                    let mut oritan = o;
                    if ee_is_new && !same {
                        oritan = o.complement();
                    }
                    (oritan, false)
                }
                None => return Err("TopOpeBRepDS_FaceInterferenceTool".to_string()),
            },
        };

        let uv = if egofft {
            match fun_tool_par_on_e_f(&e.0, self.par_on_ed, &ft.0) {
                Some(uv) => uv,
                None => {
                    let (u, v, _d) = fun_tool_proj_pon_f(&self.pnt_on_ed, &ft)
                        .ok_or_else(|| "TopOpeBRepDS_FaceInterferenceTool".to_string())?;
                    occt_core::gp::GpPnt2d::new(u, v)
                }
            }
        } else {
            let (u, v, _d) = fun_tool_proj_pon_f(&self.pnt_on_ed, &ft)
                .ok_or_else(|| "TopOpeBRepDS_FaceInterferenceTool".to_string())?;
            occt_core::gp::GpPnt2d::new(u, v)
        };

        let uv = (uv.x(), uv.y());
        if self.is_line {
            let norm =
                fun_compute_geom_data(&ft, uv).ok_or_else(|| "FUN_ComputeGeomData".to_string())?;
            self.tool.compare_plane(self.tole, &norm, ori_loc, ori_tan);
        } else {
            let (norm, d1, d2, cur1, cur2) = fun_compute_geom_data_full(&ft, uv)?;
            self.tool
                .compare_curved(self.tole, &norm, &d1, &d2, cur1, cur2, ori_loc, ori_tan);
        }
        Ok(())
    }

    /// `Add(E, C, I)` (`cxx:386-392`): `// NYI`.
    pub fn add_curve(&mut self, _e: &TopoShape, _i: &Interference) {
        // UNPORTED: `TopOpeBRepDS_FaceInterferenceTool.cxx:386-392` is an empty
        // body with the sole comment `// NYI`.
    }

    /// `Transition(I)` (`cxx:396-424`).
    pub fn transition(&self, i: &mut Interference) {
        if self.face_orientation == Orientation::Internal {
            let t = i.change_transition();
            t.set_state_before(State::In);
            t.set_state_after(State::In);
        } else if self.face_orientation == Orientation::External {
            let t = i.change_transition();
            t.set_state_before(State::Out);
            t.set_state_after(State::Out);
        } else {
            i.set_support(self.face_oriented);
            let stb = self.tool.state_before();
            let sta = self.tool.state_after();
            let t = i.change_transition();
            // `T.Set(stb, sta)`: the shapes default to `TopAbs_FACE`
            // (`TopOpeBRepDS_Transition.hxx:42-45`).
            t.set(
                stb,
                sta,
                crate::abs::ShapeType::Face,
                crate::abs::ShapeType::Face,
            );
            // xpu150498
            let mut o = t.orientation(State::In, crate::abs::ShapeType::Face);
            let rev = self.to_reverse && (o == Orientation::Forward || o == Orientation::Reversed);
            if rev {
                o = o.complement();
            }
            t.set_orientation(o);
            // xpu150498
        }
    }

    /// `SetEdgePntPar(P, p)` (`cxx:428-433`).
    pub fn set_edge_pnt_par(&mut self, p: GpPnt, par: f64) {
        self.pnt_on_ed = p;
        self.par_on_ed = par;
        self.on_ed_def = true;
    }

    /// `GetEdgePntPar(P, p)` (`cxx:436-446`): throws `Standard_ProgramError`
    /// when the point was never set.
    pub fn get_edge_pnt_par(&self) -> Result<(GpPnt, f64), String> {
        if !self.on_ed_def {
            return Err("GetEdgePntPar".to_string());
        }
        Ok((self.pnt_on_ed, self.par_on_ed))
    }

    /// `IsEdgePntParDef()` (`cxx:450-453`).
    pub fn is_edge_pnt_par_def(&self) -> bool {
        self.on_ed_def
    }
}

/// `TopOpeBRepDS_Kind` values this tool reads as support/geometry kinds
/// (`MDSkf`, `MDSke`) — kept named so the `Add` guards read like `cxx:449-451`.
#[allow(dead_code)]
const _KINDS: (Kind, Kind) = (Kind::Face, Kind::Edge);
