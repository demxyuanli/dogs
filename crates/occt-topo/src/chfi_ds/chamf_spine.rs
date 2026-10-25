//! `ChFiDS_ChamfSpine` — the guideline of a chamfer, i.e. `ChFiDS_Spine` plus
//! the section distances.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_ChamfSpine.hxx:20-66`,
//! `ChFiDS_ChamfSpine.cxx:27-126` (OCCT 8.0.0).
//!
//! `ChFiDS_ChamfSpine` derives from `ChFiDS_Spine`; a `SpineHandle` is
//! `Arc<RwLock<Spine>>` and cannot hold a derived class, so the four members
//! OCCT declares at `hxx:53-56` live in [`super::spine::Spine`] as
//! [`ChamfSpineData`] and the methods below are written on `Spine` (see
//! [`super::spine::SpineClass`]).
//!
//! Divergence that is not a heuristic: the accessors that raise
//! `Standard_Failure` when the chamfer is of another kind (`cxx:47-52`,
//! `:66-72`, `:88-95`) answer `Option` instead.

use super::enums::{ChamfMethod, ChamfMode};
use super::spine::Spine;

/// `d1`, `d2`, `angle`, `mChamf` (`ChFiDS_ChamfSpine.hxx:53-56`).
#[derive(Debug, Clone, Copy)]
pub struct ChamfSpineData {
    /// `d1`
    d1: f64,
    /// `d2`
    d2: f64,
    /// `angle`
    angle: f64,
    /// `mChamf`
    mchamf: ChamfMethod,
}

impl ChamfSpineData {
    /// The member values `ChFiDS_ChamfSpine::ChFiDS_ChamfSpine`
    /// (`cxx:27-42`) writes: three zeros and `ChFiDS_Sym`.
    pub fn new() -> Self {
        Self {
            d1: 0.0,
            d2: 0.0,
            angle: 0.0,
            mchamf: ChamfMethod::Sym,
        }
    }

    /// `d1`.
    pub fn d1(&self) -> f64 {
        self.d1
    }

    /// `d2`.
    pub fn d2(&self) -> f64 {
        self.d2
    }

    /// `angle`.
    pub fn angle(&self) -> f64 {
        self.angle
    }

    /// `mChamf` — what [`Spine::is_chamfer`] answers.
    pub fn method(&self) -> ChamfMethod {
        self.mchamf
    }
}

impl Default for ChamfSpineData {
    fn default() -> Self {
        Self::new()
    }
}

impl Spine {
    /// `GetDist(Dis)` (`cxx:47-53`): the symmetric distance.
    pub fn get_dist(&self) -> Option<f64> {
        if self.chamf.mchamf != ChamfMethod::Sym {
            return None;
        }
        Some(self.chamf.d1)
    }

    /// `SetDist(Dis)` (`cxx:57-61`): selects the symmetric kind.
    pub fn set_dist(&mut self, dis: f64) {
        self.chamf.mchamf = ChamfMethod::Sym;
        self.chamf.d1 = dis;
    }

    /// `Dists(Dis1, Dis2)` (`cxx:65-73`): the two distances.
    pub fn dists(&self) -> Option<(f64, f64)> {
        if self.chamf.mchamf != ChamfMethod::TwoDist {
            return None;
        }
        Some((self.chamf.d1, self.chamf.d2))
    }

    /// `SetDists(Dis1, Dis2)` (`cxx:77-83`): selects the two-distance kind.
    pub fn set_dists(&mut self, dis1: f64, dis2: f64) {
        self.chamf.mchamf = ChamfMethod::TwoDist;
        self.chamf.d1 = dis1;
        self.chamf.d2 = dis2;
    }

    /// `GetDistAngle(Dis, Angle)` (`cxx:87-97`): the distance and the angle.
    pub fn get_dist_angle(&self) -> Option<(f64, f64)> {
        if self.chamf.mchamf != ChamfMethod::DistAngle {
            return None;
        }
        Some((self.chamf.d1, self.chamf.angle))
    }

    /// `SetDistAngle(Dis, Angle)` (`cxx:102-112`): selects the distance-angle
    /// kind.
    pub fn set_dist_angle(&mut self, dis: f64, angle: f64) {
        self.chamf.mchamf = ChamfMethod::DistAngle;
        self.chamf.d1 = dis;
        self.chamf.angle = angle;
    }

    /// `SetMode(theMode)` (`cxx:115-118`): `myMode`.
    pub fn set_chamf_mode(&mut self, the_mode: ChamfMode) {
        self.set_mode(the_mode);
    }

    /// `IsChamfer()` (`cxx:122-126`): `mChamf`, despite the name — OCCT returns
    /// the kind of section, not a boolean (`hxx:62-64`).
    pub fn is_chamfer(&self) -> ChamfMethod {
        self.chamf.mchamf
    }
}
