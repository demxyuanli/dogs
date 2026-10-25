//! `GeomFill_Pipe::Init` (`GeomFill_Pipe.cxx:364-728`, OCCT 8.0.0).
//!
//! Every `Init` stores the path adaptor / first / last section and builds the
//! `LocationLaw` + `SectionLaw` pair the sweep needs. Which of those arms is
//! reachable depends on classes that this tree has not ported yet; those arms
//! are marked UNPORTED with the OCCT source they need, rather than replaced by
//! an approximation.
//!
//! UNPORTED:
//! - `Init(const Geom2d_Curve&, const Geom_Surface&, const Geom_Curve&)`
//!   (`cxx:564-581`): needs `Adaptor3d_CurveOnSurface`, which does not exist in
//!   `occt-geom` (`occt-topo` has `BRepMesh_CurveOnSurface`, a downstream
//!   class). The `GeomFill_Darboux` arm itself is ported.

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx2, GpCirc, GpDir, GpPnt, GpTrsf, GpVec};
use occt_core::precision::{CONFUSION, PCONFUSION};

use crate::circle::GeomCircle;
use crate::curve::Curve;
use crate::geom_fill::constant_bi_normal::ConstantBiNormal;
use crate::geom_fill::corrected_frenet::CorrectedFrenet;
use crate::geom_fill::curve_and_trihedron::CurveAndTrihedron;
use crate::geom_fill::fixed::Fixed;
use crate::geom_fill::frenet::Frenet;
use crate::geom_fill::guide_trihedron_ac::GuideTrihedronAC;
use crate::geom_fill::guide_trihedron_plan::GuideTrihedronPlan;
use crate::geom_fill::location_guide::LocationGuide;
use crate::geom_fill::location_law::LocationLaw;
use crate::geom_fill::n_sections::NSections;
use crate::geom_fill::pipe::check_sense::check_sense;
use crate::geom_fill::pipe::Pipe;
use crate::geom_fill::pipe_error::GeomFillPipeError;
use crate::geom_fill::section_law::SectionLaw;
use crate::geom_fill::section_placement::{SectionGeometry, SectionPlacement};
use crate::geom_fill::trihedron::GeomFillTrihedron;
use crate::geom_fill::trihedron_law::TrihedronLaw;
use crate::geom_fill::trihedron_with_guide::TrihedronWithGuide;
use crate::geom_fill::uniform_section::UniformSection;
use crate::geom_lprop::cl_props::ClProps;

/// `gp::XOY()` (`gp.hxx`): the standard 3D coordinate system.
fn xoy() -> GpAx2 {
    GpAx2::standard()
}

/// `gp::OZ()` (`gp.hxx`).
fn oz() -> GpAx1 {
    GpAx1::new(GpPnt::new(0.0, 0.0, 0.0), GpDir::from_axis(occt_core::gp::dir::DirAxis::Z))
}

impl Pipe {
    /// `GeomFill_Pipe::Init(Path, Radius)` (`cxx:431-448`).
    pub fn init_path_radius(
        &mut self,
        path: Arc<dyn Curve>,
        radius: f64,
    ) -> Result<(), String> {
        self.pipe_type = 1;
        self.error = 0.0;
        self.radius = radius;

        self.adp_path = Some(path.clone());

        // `new Geom_Circle(gp::XOY(), Radius); C->Rotate(gp::OZ(), M_PI/2)`
        // (`cxx:439-440`).
        let mut circ = GpCirc::new(xoy(), radius);
        let mut rot = GpTrsf::identity();
        rot.set_rotation_ax1(&oz(), std::f64::consts::FRAC_PI_2)
            .map_err(|e| format!("gp_Trsf::SetRotation: {e}"))?;
        circ.transform(&rot);

        let section: Arc<dyn Curve> = Arc::new(GeomCircle::new(circ));
        self.sec = Some(Arc::new(UniformSection::new(
            section,
            path.first_parameter(),
            path.last_parameter(),
        )?));

        let tlaw: Arc<dyn TrihedronLaw> = Arc::new(CorrectedFrenet::new());
        let loc = CurveAndTrihedron::new(tlaw);
        loc.set_curve(path);
        self.loc = Some(Arc::new(loc));
        Ok(())
    }

    /// `GeomFill_Pipe::Init(Path, FirstSect, Option)` (`cxx:452-558`).
    pub fn init_path_section_option(
        &mut self,
        path: Arc<dyn Curve>,
        first_sect: Arc<dyn Curve>,
        option: GeomFillTrihedron,
    ) -> Result<(), String> {
        self.adp_path = Some(path.clone());
        let mut param = path.first_parameter();

        let tlaw: Option<Arc<dyn TrihedronLaw>>;

        match option {
            GeomFillTrihedron::CorrectedFrenet => {
                tlaw = Some(Arc::new(CorrectedFrenet::new()));
            }
            // OCCT falls through `case GeomFill_IsDarboux:` into
            // `case GeomFill_IsFrenet:` (`cxx:466-477`); the debug print in
            // between is `#ifdef OCCT_DEBUG` only.
            GeomFillTrihedron::Darboux | GeomFillTrihedron::Frenet => {
                tlaw = Some(Arc::new(Frenet::new()));
            }
            GeomFillTrihedron::Fixed => {
                let eps = 1.0e-9;
                let (v1, v2) = tangent_normal_vectors(path.as_ref(), param, eps)?;
                tlaw = Some(Arc::new(Fixed::new(&v1, &v2)?));
            }
            GeomFillTrihedron::ConstantNormal => {
                let frenet: Arc<dyn TrihedronLaw> = Arc::new(Frenet::new());
                let loc = CurveAndTrihedron::new(frenet);
                loc.set_curve(path.clone());
                let loc: Arc<dyn LocationLaw> = Arc::new(loc);
                let mut place =
                    SectionPlacement::new(loc.clone(), SectionGeometry::Curve(first_sect.clone()));
                place.perform(CONFUSION);
                let ponsec = place.parameter_on_section();

                let eps = 1.0e-9;
                let v = constant_normal_vector(first_sect.as_ref(), ponsec, eps)?;
                let tlaw_cb: Arc<dyn TrihedronLaw> = Arc::new(ConstantBiNormal::new(
                    &GpDir::from_vec(&v)
                        .map_err(|_| "gp_Dir: null magnitude".to_string())?,
                ));
                tlaw = Some(tlaw_cb);
                // `myLoc` is replaced by the block below, but OCCT leaves the
                // Frenet-based location in place until then (`cxx:509-511`).
                self.loc = Some(loc);
            }
            _ => {
                return Err("GeomFill::Init : Unknown Option".into());
            }
        }

        if let Some(tlaw) = tlaw {
            let loc = CurveAndTrihedron::new(tlaw);
            loc.set_curve(path.clone());
            let loc: Arc<dyn LocationLaw> = Arc::new(loc);
            let mut place = SectionPlacement::new(loc.clone(), SectionGeometry::Curve(first_sect));
            place.perform(CONFUSION);
            param = place.parameter_on_path();
            let sect = place
                .section(false)?
                .ok_or_else(|| "GeomFill_SectionPlacement: null section".to_string())?;

            self.sec = Some(Arc::new(UniformSection::new(
                sect,
                path.first_parameter(),
                path.last_parameter(),
            )?));
            self.loc = Some(loc);
        }
        // `param` is only read by the `Place.ParameterOnPath()` above; OCCT
        // never stores it (`cxx:552`), so the local stays unused there too.
        let _ = param;
        Ok(())
    }

    /// `GeomFill_Pipe::Init(Path, FirstSect, Direction)` (`cxx:586-604`).
    pub fn init_path_first_section_dir(
        &mut self,
        path: Arc<dyn Curve>,
        first_sect: Arc<dyn Curve>,
        direction: GpDir,
    ) -> Result<(), String> {
        self.init();

        self.adp_path = Some(path.clone());

        let tlaw: Arc<dyn TrihedronLaw> = Arc::new(ConstantBiNormal::new(&direction));
        let loc = CurveAndTrihedron::new(tlaw);
        loc.set_curve(path.clone());
        let loc: Arc<dyn LocationLaw> = Arc::new(loc);

        let mut place = SectionPlacement::new(loc.clone(), SectionGeometry::Curve(first_sect));
        place.perform(CONFUSION);
        let sect = place
            .section(false)?
            .ok_or_else(|| "GeomFill_SectionPlacement: null section".to_string())?;

        self.sec = Some(Arc::new(UniformSection::new(
            sect,
            path.first_parameter(),
            path.last_parameter(),
        )?));
        self.loc = Some(loc);
        Ok(())
    }

    /// `GeomFill_Pipe::Init(Path, FirstSect, LastSect)` (`cxx:679-728`).
    pub fn init_path_first_last(
        &mut self,
        path: Arc<dyn Curve>,
        first_sect: Arc<dyn Curve>,
        last_sect: Arc<dyn Curve>,
    ) -> Result<(), String> {
        // `cxx:683-685`. OCCT does not call `Init()` here, so `myKPart` /
        // `myStatus` / `mySurface` keep their previous values.
        self.pipe_type = 3;
        self.error = 0.0;
        self.radius = 0.0;
        let first = path.first_parameter();
        let last = path.last_parameter();
        self.adp_path = Some(path.clone());

        // `cxx:687-690`.
        let tlaw: Arc<dyn TrihedronLaw> = Arc::new(CorrectedFrenet::new());

        // `cxx:692-697`: a failed `SetCurve` only records
        // `GeomFill_ImpossibleContact`; it is not an exception.
        let loc = CurveAndTrihedron::new(tlaw);
        let loc: Arc<dyn LocationLaw> = Arc::new(loc);
        self.loc = Some(loc.clone());
        if !loc.set_curve(path) {
            self.status = GeomFillPipeError::ImpossibleContact;
            return Ok(());
        }

        // `cxx:702-712`: both placements use `Perform(first, Confusion)`, so
        // even the last section is placed at the *first* path parameter.
        let mut seq_c: Vec<Arc<dyn Curve>> = Vec::with_capacity(2);
        for section in [first_sect, last_sect] {
            let mut place = SectionPlacement::new(loc.clone(), SectionGeometry::Curve(section));
            place.perform_at(first, CONFUSION);
            seq_c.push(
                place
                    .section(false)?
                    .ok_or_else(|| "GeomFill_SectionPlacement: null section".to_string())?,
            );
        }
        let seq_p = vec![first, last];

        // `cxx:717-721`.
        let (sense_changed, new_seq) = check_sense(&seq_c);
        if sense_changed {
            seq_c = new_seq;
        }

        // `cxx:723-727`.
        if seq_c.is_empty() {
            return Err("Standard_ConstructionError: GeomFill_Pipe::Init: empty sections".into());
        }
        let deb = seq_c[0].first_parameter();
        let fin = seq_c[0].last_parameter();
        self.sec = Some(Arc::new(NSections::with_bounds(
            seq_c, seq_p, deb, fin, first, last,
        )?));
        Ok(())
    }

    /// `GeomFill_Pipe::Init(Path, NSections)` (`cxx:606-674`).
    pub fn init_path_sections(
        &mut self,
        path: Arc<dyn Curve>,
        nsections: Vec<Arc<dyn Curve>>,
    ) -> Result<(), String> {
        // `cxx:611-613`. As in `Init(Path, FirstSect, LastSect)`, OCCT does not
        // call `Init()` first.
        self.pipe_type = 3;
        self.error = 0.0;
        self.radius = 0.0;
        self.adp_path = Some(path.clone());

        // `cxx:616-622`.
        let tlaw: Arc<dyn TrihedronLaw> = Arc::new(CorrectedFrenet::new());
        let loc = CurveAndTrihedron::new(tlaw);
        let loc: Arc<dyn LocationLaw> = Arc::new(loc);
        self.loc = Some(loc.clone());
        loc.set_curve(path.clone());

        // `cxx:623-633`: place every section and keep its path parameter.
        let mut seq_c: Vec<Arc<dyn Curve>> = Vec::with_capacity(nsections.len());
        let mut seq_p: Vec<f64> = Vec::with_capacity(nsections.len());
        for section in nsections {
            let mut place = SectionPlacement::new(loc.clone(), SectionGeometry::Curve(section));
            place.perform(CONFUSION);
            seq_p.push(place.parameter_on_path());
            seq_c.push(
                place
                    .section(false)?
                    .ok_or_else(|| "GeomFill_SectionPlacement: null section".to_string())?,
            );
        }

        // `cxx:635-640`.
        let (sense_changed, new_seq) = check_sense(&seq_c);
        if sense_changed {
            seq_c = new_seq;
        }

        // `cxx:642-659`: order the sections by path parameter. OCCT's bubble
        // loop restarts after every exchange, so replicate the flag.
        let mut play_again = true;
        while play_again {
            play_again = false;
            let n = seq_c.len();
            for i in 0..n {
                for j in i..n {
                    if seq_p[i] > seq_p[j] {
                        seq_p.swap(i, j);
                        seq_c.swap(i, j);
                        play_again = true;
                    }
                }
            }
        }

        // `cxx:660-667`.
        for i in 0..seq_c.len().saturating_sub(1) {
            if (seq_p[i + 1] - seq_p[i]).abs() < PCONFUSION {
                return Err("Standard_ConstructionError: GeomFill_Pipe::Init with NSections : \
                            invalid parameters"
                    .into());
            }
        }

        // `cxx:669-672`.
        if seq_c.is_empty() {
            return Err("Standard_ConstructionError: GeomFill_Pipe::Init: empty sections".into());
        }
        let first = path.first_parameter();
        let last = path.last_parameter();
        let deb = seq_c[0].first_parameter();
        let fin = seq_c[0].last_parameter();
        self.sec = Some(Arc::new(NSections::with_bounds(
            seq_c, seq_p, deb, fin, first, last,
        )?));
        Ok(())
    }

    /// `GeomFill_Pipe::Init(const Geom2d_Curve&, const Geom_Surface&, const
    /// Geom_Curve&)` (`cxx:564-581`).
    ///
    /// UNPORTED: `myAdpPath` is an `Adaptor3d_CurveOnSurface`
    /// (`cxx:570-571`), which `occt-geom` does not have; `GeomFill_Darboux` and
    /// `GeomFill_SectionPlacement::Perform(Path, Tol)` are ported.
    pub fn init_path2d_support_section(
        &mut self,
        _path: Arc<dyn Curve>,
        _support: Arc<dyn crate::surface::Surface>,
        _first_sect: Arc<dyn Curve>,
    ) -> Result<(), String> {
        Err(
            "UNPORTED: GeomFill_Pipe::Init(Geom2d_Curve, Geom_Surface, Geom_Curve) needs \
             Adaptor3d_CurveOnSurface (GeomFill_Pipe.cxx:564-581)"
                .into(),
        )
    }

    /// `GeomFill_Pipe::Init(Path, Guide, FirstSect, ByACR, rotat)`
    /// (`cxx:364-411`).
    pub fn init_with_guide(
        &mut self,
        path: Arc<dyn Curve>,
        guide: Arc<dyn Curve>,
        first_sect: Arc<dyn Curve>,
        by_acr: bool,
        rotat: bool,
    ) -> Result<(), String> {
        self.init();
        self.adp_path = Some(path.clone());

        // `cxx:378-391`: the trihedron law driven by the guide.
        let tlaw: Arc<dyn TrihedronWithGuide> = if by_acr {
            Arc::new(GuideTrihedronAC::new(guide))
        } else {
            Arc::new(GuideTrihedronPlan::new(guide))
        };
        tlaw.set_curve(path.clone());

        // `cxx:394-395`.
        let the_loc = Arc::new(LocationGuide::new(tlaw)?);
        the_loc.set_curve(path.clone());

        // `cxx:397-398`.
        let mut place = SectionPlacement::new(
            the_loc.clone() as Arc<dyn LocationLaw>,
            SectionGeometry::Curve(first_sect),
        );
        place.perform(CONFUSION);
        let sect = place
            .section(false)?
            .ok_or_else(|| "GeomFill_SectionPlacement: null section".to_string())?;

        // `cxx:401-403`.
        let sec = Arc::new(UniformSection::new(
            sect,
            path.first_parameter(),
            path.last_parameter(),
        )?);
        self.sec = Some(sec.clone() as Arc<dyn SectionLaw>);

        if rotat {
            // `cxx:405-408`: `Set(mySec, rotat, first, last, 0., angle)`.
            let _angle = the_loc.set(
                Some(sec as Arc<dyn SectionLaw>),
                true,
                path.first_parameter(),
                path.last_parameter(),
                0.0,
            )?;
        }
        self.loc = Some(the_loc as Arc<dyn LocationLaw>);
        Ok(())
    }

    /// `GeomFill_Pipe::Init(const Adaptor3d_Curve&, const Adaptor3d_Curve&,
    /// const Adaptor3d_Curve&, Radius)` (`cxx:732-742`): the type 4 "pipe
    /// between two guide lines" setup, no law construction.
    pub fn init_adaptors(
        &mut self,
        path: Arc<dyn Curve>,
        curve1: Arc<dyn Curve>,
        curve2: Arc<dyn Curve>,
        radius: f64,
    ) {
        self.pipe_type = 4;
        self.error = 0.0;
        self.radius = radius;
        self.adp_path = Some(path);
        self.adp_first_sect = Some(curve1);
        self.adp_last_sect = Some(curve2);
    }

    /// `myAdpPath` (`hxx:304`).
    pub(crate) fn adp_path(&self) -> Option<&Arc<dyn Curve>> {
        self.adp_path.as_ref()
    }

    /// `myAdpFirstSect` (`hxx:305`).
    pub(crate) fn adp_first_sect(&self) -> Option<&Arc<dyn Curve>> {
        self.adp_first_sect.as_ref()
    }

    /// `myAdpLastSect` (`hxx:306`).
    pub(crate) fn adp_last_sect(&self) -> Option<&Arc<dyn Curve>> {
        self.adp_last_sect.as_ref()
    }

    /// `myType` (`hxx:313`).
    pub(crate) fn pipe_type(&self) -> i32 {
        self.pipe_type
    }

    /// `myRadius` (`hxx:302`).
    pub(crate) fn radius(&self) -> f64 {
        self.radius
    }

    /// `mySec` (`hxx:312`).
    pub(crate) fn sec(&self) -> Option<&Arc<dyn SectionLaw>> {
        self.sec.as_ref()
    }

    /// `myLoc` (`hxx:311`).
    pub(crate) fn loc(&self) -> Option<&Arc<dyn LocationLaw>> {
        self.loc.as_ref()
    }

    /// `myKPart` (`hxx:315`).
    pub(crate) fn k_part(&self) -> bool {
        self.k_part
    }
}

/// The `GeomLProp_CLProps` + `Curvature` block of the `GeomFill_IsFixed` arm
/// (`cxx:479-504`): returns the `(tangent, normal)` pair handed to
/// `GeomFill_Fixed`.
fn tangent_normal_vectors(
    curve: &dyn Curve,
    param: f64,
    eps: f64,
) -> Result<(GpVec, GpVec), String> {
    let mut v1 = GpVec::new(0.0, 0.0, 1.0);
    let mut v2 = GpVec::new(0.0, 1.0, 0.0);
    let mut cp = ClProps::new(curve, param, 2, eps)?;
    if cp.is_tangent_defined() {
        let d = cp.tangent()?;
        v1 = GpVec::from_xyz(d.xyz());
        v1.normalize();
        if cp.curvature()? > eps {
            let d = cp.normal()?;
            v2 = GpVec::from_xyz(d.xyz());
            v2.normalize();
        } else {
            let axe = GpAx2::from_axis(GpPnt::new(0.0, 0.0, 0.0), d);
            let dx = *axe.x_direction();
            v2 = GpVec::from_xyz(dx.xyz());
            v2.normalize();
        }
    }
    Ok((v1, v2))
}

/// The `GeomLProp_CLProps` block of the `GeomFill_IsConstantNormal` arm
/// (`cxx:513-538`): the constant bi-normal `V` starts as `(0, 1, 0)` and is
/// replaced by the section normal when that normal is defined.
fn constant_normal_vector(curve: &dyn Curve, param: f64, eps: f64) -> Result<GpVec, String> {
    let mut v = GpVec::new(0.0, 1.0, 0.0);
    let mut cp = ClProps::new(curve, param, 2, eps)?;
    if cp.is_tangent_defined() {
        let d = cp.tangent()?;
        if cp.curvature()? > eps {
            let n = cp.normal()?;
            v = GpVec::from_xyz(n.xyz());
            v.normalize();
        } else {
            let axe = GpAx2::from_axis(GpPnt::new(0.0, 0.0, 0.0), d);
            let dx = *axe.x_direction();
            v = GpVec::from_xyz(dx.xyz());
            v.normalize();
        }
    }
    Ok(v)
}
