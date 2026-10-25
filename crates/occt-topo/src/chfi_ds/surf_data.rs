//! `ChFiDS_SurfData` — the data of one fillet patch, shared by the two faces
//! it sits between.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_SurfData.hxx:30-181`,
//! `ChFiDS_SurfData.lxx:18-250`, `ChFiDS_SurfData.cxx:24-281` (OCCT 8.0.0).

use std::sync::Arc;

use occt_core::gp::GpPnt2d;

use crate::abs::Orientation;

use super::circ_section::CircSectionArray;
use super::common_point::CommonPoint;
use super::face_interference::FaceInterference;

/// Source: `ChFiDS_SurfData.hxx:32-181`.
///
/// In OCCT this derives from `Standard_Transient` because
/// `ChFiDS_Stripe::SetOfSurfData` holds handles; the port uses `Arc<SurfData>`
/// at those sites instead.
#[derive(Debug, Clone)]
pub struct SurfData {
    /// `pfirstOnS1`
    pfirst_on_s1: CommonPoint,
    /// `plastOnS1`
    plast_on_s1: CommonPoint,
    /// `pfirstOnS2`
    pfirst_on_s2: CommonPoint,
    /// `plastOnS2`
    plast_on_s2: CommonPoint,
    /// `intf1`
    intf1: FaceInterference,
    /// `intf2`
    intf2: FaceInterference,
    /// `p2df1`
    p2df1: GpPnt2d,
    /// `p2dl1`
    p2dl1: GpPnt2d,
    /// `p2df2`
    p2df2: GpPnt2d,
    /// `p2dl2`
    p2dl2: GpPnt2d,
    /// `ufspine`
    ufspine: f64,
    /// `ulspine`
    ulspine: f64,
    /// `myfirstextend`
    myfirstextend: f64,
    /// `mylastextend`
    mylastextend: f64,
    /// `simul` — `occ::handle<Standard_Transient>`; the only thing OCCT ever
    /// stores in it is an `NCollection_HArray1<ChFiDS_CircSection>`
    /// (`ChFi3d_FilBuilder.cxx:711`, read back by `down_cast` at `:466`).
    simul: Option<Arc<CircSectionArray>>,
    /// `indexOfS1`
    index_of_s1: i32,
    /// `indexOfC1`
    index_of_c1: i32,
    /// `indexOfS2`
    index_of_s2: i32,
    /// `indexOfC2`
    index_of_c2: i32,
    /// `indexOfConge`
    index_of_conge: i32,
    /// `isoncurv1`
    isoncurv1: bool,
    /// `isoncurv2`
    isoncurv2: bool,
    /// `twistons1`
    twistons1: bool,
    /// `twistons2`
    twistons2: bool,
    /// `orientation`
    orientation: Orientation,
}

impl Default for SurfData {
    /// `ChFiDS_SurfData()` (`cxx:24-40`).
    fn default() -> Self {
        Self {
            pfirst_on_s1: CommonPoint::new(),
            plast_on_s1: CommonPoint::new(),
            pfirst_on_s2: CommonPoint::new(),
            plast_on_s2: CommonPoint::new(),
            intf1: FaceInterference::new(),
            intf2: FaceInterference::new(),
            p2df1: GpPnt2d::default(),
            p2dl1: GpPnt2d::default(),
            p2df2: GpPnt2d::default(),
            p2dl2: GpPnt2d::default(),
            ufspine: 0.0,
            ulspine: 0.0,
            myfirstextend: 0.0,
            mylastextend: 0.0,
            simul: None,
            index_of_s1: 0,
            index_of_c1: 0,
            index_of_s2: 0,
            index_of_c2: 0,
            index_of_conge: 0,
            isoncurv1: false,
            isoncurv2: false,
            twistons1: false,
            twistons2: false,
            orientation: Orientation::Forward,
        }
    }
}

impl SurfData {
    /// `ChFiDS_SurfData()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `Copy(Other)` (`cxx:44-73`).
    ///
    /// OCCT copies `indexOfS1` / `indexOfS2` / `indexOfConge` but **not**
    /// `indexOfC1` / `indexOfC2` / `isoncurv1` / `isoncurv2` — the two
    /// curve-support flags are left at their current value. That asymmetry is
    /// in OCCT and is kept here.
    pub fn copy_from(&mut self, other: &SurfData) {
        self.index_of_s1 = other.index_of_s1;
        self.index_of_s2 = other.index_of_s2;
        self.index_of_conge = other.index_of_conge;
        self.orientation = other.orientation;
        self.intf1 = other.intf1.clone();
        self.intf2 = other.intf2.clone();

        self.pfirst_on_s1 = other.pfirst_on_s1.clone();
        self.plast_on_s1 = other.plast_on_s1.clone();
        self.pfirst_on_s2 = other.pfirst_on_s2.clone();
        self.plast_on_s2 = other.plast_on_s2.clone();

        self.ufspine = other.ufspine;
        self.ulspine = other.ulspine;

        self.simul = other.simul.clone();
        self.p2df1 = other.p2df1;
        self.p2dl1 = other.p2dl1;
        self.p2df2 = other.p2df2;
        self.p2dl2 = other.p2dl2;

        self.myfirstextend = other.myfirstextend;
        self.mylastextend = other.mylastextend;

        self.twistons1 = other.twistons1;
        self.twistons2 = other.twistons2;
    }

    /// `IndexOfS1()` (`lxx:19-22`).
    pub fn index_of_s1(&self) -> i32 {
        self.index_of_s1
    }

    /// `IndexOfS2()` (`lxx:51-54`).
    pub fn index_of_s2(&self) -> i32 {
        self.index_of_s2
    }

    /// `IsOnCurve1()` (`lxx:27-30`).
    pub fn is_on_curve1(&self) -> bool {
        self.isoncurv1
    }

    /// `IsOnCurve2()` (`lxx:57-60`).
    pub fn is_on_curve2(&self) -> bool {
        self.isoncurv2
    }

    /// `IndexOfC1()` (`lxx:34-38`): OCCT raises `Standard_Failure` when the
    /// interference is not on a curve.
    pub fn index_of_c1(&self) -> Option<i32> {
        if !self.isoncurv1 {
            return None;
        }
        Some(self.index_of_c1)
    }

    /// `IndexOfC2()` (`lxx:63-67`).
    pub fn index_of_c2(&self) -> Option<i32> {
        if !self.isoncurv2 {
            return None;
        }
        Some(self.index_of_c2)
    }

    /// `Surf()` (`lxx:81-84`).
    pub fn surf(&self) -> i32 {
        self.index_of_conge
    }

    /// `Orientation()` (`lxx:88-91`).
    pub fn orientation(&self) -> Orientation {
        self.orientation
    }

    /// `InterferenceOnS1()` (`lxx:95-98`).
    pub fn interference_on_s1(&self) -> &FaceInterference {
        &self.intf1
    }

    /// `InterferenceOnS2()` (`lxx:101-104`).
    pub fn interference_on_s2(&self) -> &FaceInterference {
        &self.intf2
    }

    /// `ChangeInterferenceOnS1()` (`lxx:164-167`) read as the reference OCCT
    /// returns, so that `Data->ChangeInterferenceOnS1().SetInterference(...)`
    /// (`ChFiKPart_ComputeData_FilPlnPln.cxx:140-143`) writes in place.
    pub fn interference_on_s1_mut(&mut self) -> &mut FaceInterference {
        &mut self.intf1
    }

    /// `ChangeInterferenceOnS2()` (`lxx:170-173`), the `intf2` mirror of
    /// [`Self::interference_on_s1_mut`].
    pub fn interference_on_s2_mut(&mut self) -> &mut FaceInterference {
        &mut self.intf2
    }

    /// `VertexFirstOnS1()` (`lxx:108-111`).
    pub fn vertex_first_on_s1(&self) -> &CommonPoint {
        &self.pfirst_on_s1
    }

    /// `VertexFirstOnS2()` (`lxx:126-129`).
    pub fn vertex_first_on_s2(&self) -> &CommonPoint {
        &self.pfirst_on_s2
    }

    /// `VertexLastOnS1()` (`lxx:115-118`).
    pub fn vertex_last_on_s1(&self) -> &CommonPoint {
        &self.plast_on_s1
    }

    /// `VertexLastOnS2()` (`lxx:133-136`).
    pub fn vertex_last_on_s2(&self) -> &CommonPoint {
        &self.plast_on_s2
    }

    /// `ChangeVertexFirstOnS1()` (`lxx:108-111`) read as the reference it
    /// returns, so that `Data->ChangeVertexFirstOnS1().SetPoint(...)`
    /// (`ChFiKPart_ComputeData.cxx:701`) writes in place.
    pub fn vertex_first_on_s1_mut(&mut self) -> &mut CommonPoint {
        &mut self.pfirst_on_s1
    }

    /// `ChangeVertexLastOnS1()` (`lxx:115-118`).
    pub fn vertex_last_on_s1_mut(&mut self) -> &mut CommonPoint {
        &mut self.plast_on_s1
    }

    /// `ChangeVertexFirstOnS2()` (`lxx:126-129`).
    pub fn vertex_first_on_s2_mut(&mut self) -> &mut CommonPoint {
        &mut self.pfirst_on_s2
    }

    /// `ChangeVertexLastOnS2()` (`lxx:133-136`).
    pub fn vertex_last_on_s2_mut(&mut self) -> &mut CommonPoint {
        &mut self.plast_on_s2
    }

    /// `ChangeIndexOfS1(Index)` (`lxx:138-141`).
    pub fn set_index_of_s1(&mut self, index: i32) {
        self.index_of_s1 = index;
    }

    /// `ChangeIndexOfS2(Index)` (`lxx:144-147`).
    pub fn set_index_of_s2(&mut self, index: i32) {
        self.index_of_s2 = index;
    }

    /// `ChangeSurf(Index)` (`lxx:151-154`).
    pub fn set_surf(&mut self, index: i32) {
        self.index_of_conge = index;
    }

    /// `SetIndexOfC1(Index)` (`lxx:42-46`): a non-zero index also flips
    /// `isoncurv1` on.
    pub fn set_index_of_c1(&mut self, index: i32) {
        self.index_of_c1 = index;
        self.isoncurv1 = index != 0;
    }

    /// `SetIndexOfC2(Index)` (`lxx:72-76`).
    pub fn set_index_of_c2(&mut self, index: i32) {
        self.index_of_c2 = index;
        self.isoncurv2 = index != 0;
    }

    /// `ChangeOrientation()` (`lxx:158-161`).
    pub fn set_orientation(&mut self, orientation: Orientation) {
        self.orientation = orientation;
    }

    /// `ChangeInterferenceOnS1()` (`lxx:164-167`).
    pub fn set_interference_on_s1(&mut self, intf: FaceInterference) {
        self.intf1 = intf;
    }

    /// `ChangeInterferenceOnS2()` (`lxx:170-173`).
    pub fn set_interference_on_s2(&mut self, intf: FaceInterference) {
        self.intf2 = intf;
    }

    /// `ChangeVertexFirstOnS1()` (`lxx:177-180`).
    pub fn set_vertex_first_on_s1(&mut self, p: CommonPoint) {
        self.pfirst_on_s1 = p;
    }

    /// `ChangeVertexLastOnS1()` (`lxx:184-187`).
    pub fn set_vertex_last_on_s1(&mut self, p: CommonPoint) {
        self.plast_on_s1 = p;
    }

    /// `ChangeVertexFirstOnS2()` (`lxx:192-195`).
    pub fn set_vertex_first_on_s2(&mut self, p: CommonPoint) {
        self.pfirst_on_s2 = p;
    }

    /// `ChangeVertexLastOnS2()` (`lxx:199-202`).
    pub fn set_vertex_last_on_s2(&mut self, p: CommonPoint) {
        self.plast_on_s2 = p;
    }

    /// `Interference(OnS)` (`cxx:90-102`): `OnS == 1` selects `intf1`,
    /// anything else `intf2`.
    pub fn interference(&self, on_s: i32) -> &FaceInterference {
        if on_s == 1 {
            &self.intf1
        } else {
            &self.intf2
        }
    }

    /// `ChangeInterference(OnS)` (`cxx:106-115`).
    pub fn interference_mut(&mut self, on_s: i32) -> &mut FaceInterference {
        if on_s == 1 {
            &mut self.intf1
        } else {
            &mut self.intf2
        }
    }

    /// `Index(OfS)` (`cxx:76-86`).
    pub fn index(&self, of_s: i32) -> i32 {
        if of_s == 1 {
            self.index_of_s1
        } else {
            self.index_of_s2
        }
    }

    /// `Vertex(First, OnS)` (`cxx:117-136`).
    pub fn vertex(&self, first: bool, on_s: i32) -> &CommonPoint {
        if first && on_s == 1 {
            &self.pfirst_on_s1
        } else if first && on_s == 2 {
            &self.pfirst_on_s2
        } else if !first && on_s == 1 {
            &self.plast_on_s1
        } else {
            &self.plast_on_s2
        }
    }

    /// `ChangeVertex(First, OnS)` (`cxx:140-158`).
    pub fn vertex_mut(&mut self, first: bool, on_s: i32) -> &mut CommonPoint {
        if first && on_s == 1 {
            &mut self.pfirst_on_s1
        } else if first && on_s == 2 {
            &mut self.pfirst_on_s2
        } else if !first && on_s == 1 {
            &mut self.plast_on_s1
        } else {
            &mut self.plast_on_s2
        }
    }

    /// `IsOnCurve(OnS)` (`lxx:205-210`).
    pub fn is_on_curve(&self, on_s: i32) -> bool {
        if on_s == 1 {
            return self.isoncurv1;
        }
        self.isoncurv2
    }

    /// `IndexOfC(OnS)` (`lxx:213-228`): OCCT raises `Standard_Failure` when the
    /// interference is not on a curve.
    pub fn index_of_c(&self, on_s: i32) -> Option<i32> {
        if on_s == 1 {
            if !self.isoncurv1 {
                return None;
            }
            return Some(self.index_of_c1);
        }
        if !self.isoncurv2 {
            return None;
        }
        Some(self.index_of_c2)
    }

    /// `FirstSpineParam()` (`cxx:163-166`).
    pub fn first_spine_param(&self) -> f64 {
        self.ufspine
    }

    /// `LastSpineParam()` (`cxx:168-171`).
    pub fn last_spine_param(&self) -> f64 {
        self.ulspine
    }

    /// `FirstSpineParam(Par)` (`cxx:175-178`).
    pub fn set_first_spine_param(&mut self, par: f64) {
        self.ufspine = par;
    }

    /// `LastSpineParam(Par)` (`cxx:181-184`).
    pub fn set_last_spine_param(&mut self, par: f64) {
        self.ulspine = par;
    }

    /// `FirstExtensionValue()` (`cxx:189-192`).
    pub fn first_extension_value(&self) -> f64 {
        self.myfirstextend
    }

    /// `LastExtensionValue()` (`cxx:195-198`).
    pub fn last_extension_value(&self) -> f64 {
        self.mylastextend
    }

    /// `FirstExtensionValue(Extend)` (`cxx:201-204`).
    pub fn set_first_extension_value(&mut self, extend: f64) {
        self.myfirstextend = extend;
    }

    /// `LastExtensionValue(Extend)` (`cxx:207-210`).
    pub fn set_last_extension_value(&mut self, extend: f64) {
        self.mylastextend = extend;
    }

    /// `Simul()` (`cxx:216-219`).
    pub fn simul(&self) -> Option<&Arc<CircSectionArray>> {
        self.simul.as_ref()
    }

    /// `SetSimul(S)` (`cxx:223-226`).
    pub fn set_simul(&mut self, s: Arc<CircSectionArray>) {
        self.simul = Some(s);
    }

    /// `ResetSimul()` (`cxx:230-233`).
    pub fn reset_simul(&mut self) {
        self.simul = None;
    }

    /// `Get2dPoints(First, OnS)` (`cxx:254-268`).
    pub fn get_2d_point(&self, first: bool, on_s: i32) -> GpPnt2d {
        if first && on_s == 1 {
            return self.p2df1;
        } else if !first && on_s == 1 {
            return self.p2dl1;
        } else if first && on_s == 2 {
            return self.p2df2;
        }
        self.p2dl2
    }

    /// `Get2dPoints(P2df1, P2dl1, P2df2, P2dl2)` (`cxx:239-249`).
    pub fn get_2d_points(&self) -> (GpPnt2d, GpPnt2d, GpPnt2d, GpPnt2d) {
        (self.p2df1, self.p2dl1, self.p2df2, self.p2dl2)
    }

    /// `Set2dPoints(P2df1, P2dl1, P2df2, P2dl2)` (`cxx:271-281`).
    pub fn set_2d_points(
        &mut self,
        p2df1: GpPnt2d,
        p2dl1: GpPnt2d,
        p2df2: GpPnt2d,
        p2dl2: GpPnt2d,
    ) {
        self.p2df1 = p2df1;
        self.p2dl1 = p2dl1;
        self.p2df2 = p2df2;
        self.p2dl2 = p2dl2;
    }

    /// `TwistOnS1()` (`lxx:232-235`).
    pub fn twist_on_s1(&self) -> bool {
        self.twistons1
    }

    /// `TwistOnS2()` (`lxx:237-240`).
    pub fn twist_on_s2(&self) -> bool {
        self.twistons2
    }

    /// `TwistOnS1(T)` (`lxx:242-245`).
    pub fn set_twist_on_s1(&mut self, t: bool) {
        self.twistons1 = t;
    }

    /// `TwistOnS2(T)` (`lxx:247-250`).
    pub fn set_twist_on_s2(&mut self, t: bool) {
        self.twistons2 = t;
    }
}
