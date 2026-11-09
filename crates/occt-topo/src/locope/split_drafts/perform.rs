//! `LocOpe_SplitDrafts::Perform` (`LocOpe_SplitDrafts.cxx:94-1436`).
//!
//! The stage order is the reference's:
//!
//! 1. seed `myMap` with every face of `myShape` and record `F`'s orientation
//!    (`cxx:132-153`);
//! 2. build the two draft planes (`NewPlane`, `cxx:155-165`);
//! 3. split `W` with the line where those planes meet (`cxx:177-351`);
//! 4. split `F` along `theW` with `LocOpe_Spliter` (`cxx:353-403`);
//! 5. map every edge of the split face to the face drafted from it, or to the
//!    other face that shares it (`cxx:405-470`);
//! 6. walk the vertices of those edges, building the new section edges
//!    (`cxx:471-866`);
//! 7. build the new boundary edges of `theW` (`cxx:868-1100`);
//! 8. rebuild every piece of `F` on its draft plane (`cxx:1102-1358`);
//! 9. rewrite `myMap` with the descendants and assemble with
//!    `LocOpe_BuildShape` (`cxx:1360-1435`).
//!
//! `Err` is the `Standard_NullObject` of `cxx:122`, the
//! `Standard_ConstructionError` of `cxx:129` / `cxx:152`, and the
//! `Standard_ConstructionError("GeomFill_Pipe : Cannot make a surface")` of
//! `cxx:192` / `cxx:438`. Every other early `return` in the reference leaves
//! `myResult` null, so it becomes `Ok(())` with no result here.

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpDir, GpPln};
use occt_core::precision::{Precision, CONFUSION};
use occt_geom::geom_fill::pipe::Pipe;
use occt_geom::{GeomLine, GeomPlane, Surface};

use crate::abs::{Orientation, ShapeType};
use crate::geom_int::IntSS;
use crate::locope::shape_maps::{ShapeIndexedMap, ShapeSet, ShapeValueMap};
use crate::locope::split_drafts::state::PerformState;
use crate::locope::split_drafts::statics::new_plane;
use crate::locope::split_drafts::SplitDrafts;
use crate::shape::{Face, TopoShape, Wire};
use crate::topo_tools::{self, AncestorMap};

/// `Standard_NullObject` (`cxx:122`).
pub const NULL_OBJECT: &str = "Standard_NullObject";
/// `Standard_ConstructionError` (`cxx:129`, `cxx:152`).
pub const CONSTRUCTION_ERROR: &str = "Standard_ConstructionError";
/// `Standard_ConstructionError("GeomFill_Pipe : Cannot make a surface")`
/// (`cxx:192`, `:438`).
pub const PIPE_ERROR: &str = "GeomFill_Pipe : Cannot make a surface";

/// The eight argument slots of the ten-parameter `Perform` (`hxx:70-79`).
#[derive(Debug, Clone)]
pub struct PerformArgs {
    pub extract_g: GpDir,
    pub n_pl_g: GpPln,
    pub angle_g: f64,
    pub extract_d: GpDir,
    pub n_pl_d: GpPln,
    pub angle_d: f64,
    pub modify_left: bool,
    pub modify_right: bool,
}

impl SplitDrafts {
    /// `Perform(F, W, Extr, NPl, Angle)` (`cxx:94-101`): the one-sided overload
    /// forwards with `Extrd == Extrg`, `NPld == NPlg`, `Angled == Angleg`,
    /// `ModifyLeft = true`, `ModifyRight = false`.
    pub fn perform_one_side(
        &mut self,
        f: &Face,
        w: &Wire,
        extr: &GpDir,
        n_pl: &GpPln,
        angle: f64,
    ) -> Result<(), String> {
        self.perform(
            f,
            w,
            &PerformArgs {
                extract_g: *extr,
                n_pl_g: n_pl.clone(),
                angle_g: angle,
                extract_d: *extr,
                n_pl_d: n_pl.clone(),
                angle_d: angle,
                modify_left: true,
                modify_right: false,
            },
        )
    }

    /// `Perform(F, W, Extrg, NPlg, Angleg, Extrd, NPld, Angled, ModLeft,
    /// ModRight)` (`cxx:105-1436`).
    pub fn perform(&mut self, f: &Face, w: &Wire, args: &PerformArgs) -> Result<(), String> {
        // `cxx:120-121`.
        self.my_result = None;
        self.my_map.clear();

        // `cxx:122-125`.
        let Some(my_shape) = self.my_shape.clone() else {
            return Err(NULL_OBJECT.to_string());
        };
        if f.0.is_null() || w.0.is_null() {
            return Err(NULL_OBJECT.to_string());
        }
        // `cxx:127-130`.
        if !args.modify_left && !args.modify_right {
            return Err(CONSTRUCTION_ERROR.to_string());
        }

        // `cxx:132-148`: seed `myMap` with every face, remember `F`'s
        // orientation.
        let mut ori_f = Orientation::Forward;
        let mut fin_s = false;
        for fac in topo_tools::faces_of(&my_shape) {
            self.my_map.bind(&fac.0, Vec::new());
            if topo_tools::is_same(&fac.0, &f.0) {
                ori_f = fac.0.orientation();
                fin_s = true;
            }
        }
        // `cxx:150-153`.
        if !fin_s {
            return Err(CONSTRUCTION_ERROR.to_string());
        }

        // `cxx:155-165`: the local face, oriented as `myShape` stores it, then
        // the two draft planes. `aLocalFace` is what both `NewPlane` calls see.
        let a_local_face = f.0.oriented(ori_f);

        let mut new_plg = GpPln::default();
        let mut normal_fg = GpAx1::default();
        let mut new_pld = GpPln::default();
        let mut normal_fd = GpAx1::default();
        let ok_g = new_plane(
            &a_local_face,
            &args.extract_g,
            &args.n_pl_g,
            args.angle_g,
            &mut new_plg,
            &mut normal_fg,
            args.modify_left,
        );
        let ok_d = new_plane(
            &a_local_face,
            &args.extract_d,
            &args.n_pl_d,
            args.angle_d,
            &mut new_pld,
            &mut normal_fd,
            args.modify_right,
        );
        if !ok_g || !ok_d {
            return Ok(());
        }

        // `cxx:172-175`.
        let new_sg: Arc<dyn Surface> = Arc::new(GeomPlane::new(new_plg));
        let new_sd: Arc<dyn Surface> = Arc::new(GeomPlane::new(new_pld));
        let line_pipe: Arc<dyn occt_geom::Curve> = Arc::new(GeomLine::new(
            occt_core::gp::GpLin::new(normal_fg),
        ));
        let mut i2s = IntSS::new();
        i2s.perform(new_sg.clone(), new_sd.clone(), CONFUSION, true, false, false);

        let mut st = PerformState {
            new_sg,
            new_sd,
            normal_fg,
            line_pipe,
            i2s,
            the_w: w.clone(),
            res: None,
            the_left: Vec::new(),
            map_w: ShapeValueMap::new(),
            map_ev: ShapeIndexedMap::new(),
            the_map_ef: AncestorMap::default(),
            map_sg: ShapeValueMap::new(),
            map_sd: ShapeValueMap::new(),
            map_son_s: ShapeValueMap::new(),
            the_map: ShapeSet::new(),
            the_faces: Vec::new(),
            subs: crate::brep_tools_substitution::BRepToolsSubstitution::new(),
        };

        // `cxx:177-351`.
        if !self.split_wire(w, f, &mut st)? {
            return Ok(());
        }
        // `cxx:353-403`.
        if !self.split_face(f, &mut st)? {
            return Ok(());
        }
        // `cxx:405-470`.
        if !self.build_generated(f, &mut st)? {
            return Ok(());
        }
        // `cxx:471-866`.
        if !self.build_vertices(f, &mut st)? {
            return Ok(());
        }
        // `cxx:868-1100`.
        if !self.build_wire_edges(&mut st)? {
            return Ok(());
        }
        // `cxx:1102-1358`.
        if !self.build_draft_faces(f, &mut st)? {
            return Ok(());
        }
        // `cxx:1360-1435`.
        self.build_descendants(f, &mut st)
    }
}

/// The `GeomFill_Pipe` sequence shared by `cxx:184-195` and `cxx:432-441`:
/// `GenerateParticularCase(true)`, `Init(theLinePipe, C)` (the
/// `GeomFill_IsCorrectedFrenet` default of `GeomFill_Pipe.hxx:191-193`),
/// `Perform(true)`.
pub(crate) fn pipe_surface(
    line_pipe: &Arc<dyn occt_geom::Curve>,
    section: Arc<dyn occt_geom::Curve>,
) -> Result<Arc<dyn Surface>, String> {
    let mut the_pipe = Pipe::new();
    the_pipe.set_generate_particular_case(true);
    the_pipe
        .init_path_section_option(
            line_pipe.clone(),
            section,
            occt_geom::geom_fill::GeomFillTrihedron::CorrectedFrenet,
        )
        .map_err(|_| PIPE_ERROR.to_string())?;
    the_pipe
        .perform(true, false)
        .map_err(|_| PIPE_ERROR.to_string())?;
    if !the_pipe.is_done() {
        return Err(PIPE_ERROR.to_string());
    }
    the_pipe.surface().ok_or_else(|| PIPE_ERROR.to_string())
}

/// `TopAbs::Reverse` (`TopAbs.hxx:74-96`).
pub(crate) fn reverse_orientation(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        other => other,
    }
}

/// `TopExp_Explorer(shape, TopAbs_EDGE)` (`cxx:200`, `:335`, `:388`, ...).
pub(crate) fn edges_of_shape(shape: &TopoShape) -> Vec<TopoShape> {
    topo_tools::explorer_shapes(shape, ShapeType::Edge)
}

/// `IntCurveSurface_HInter::Perform(HAC, HAS)` over the adaptors'
/// `GeomAdaptor_Curve::Load(C)` / `GeomAdaptor_Surface::Load(S)` ranges
/// (`cxx:208-213`, `cxx:729-733`, `cxx:1745-1762`).
///
/// The surface window follows
/// `IntCurveSurface_InterUtils::EstLimForInfSurf` (`InterUtils.pxx:137-143`):
/// only when a bound is infinite are all four clamped to `+-1.e10`.
pub(crate) fn intcs_perform(
    curve: &dyn occt_geom::Curve,
    surface: &dyn Surface,
) -> Option<crate::intcurvesurface::HInterResult> {
    let cu_range = (curve.first_parameter(), curve.last_parameter());
    let (mut u1, mut u2) = surface.u_range();
    let (mut v1, mut v2) = surface.v_range();
    if Precision::is_infinite(u1)
        || Precision::is_infinite(u2)
        || Precision::is_infinite(v1)
        || Precision::is_infinite(v2)
    {
        u1 = u1.max(-1.0e10);
        u2 = u2.min(1.0e10);
        v1 = v1.max(-1.0e10);
        v2 = v2.min(1.0e10);
    }
    crate::intcurvesurface::perform_curve_surface(curve, surface, cu_range, (u1, u2, v1, v2)).ok()
}
