//! `Blend_Point` — one point of a blending line.
//!
//! Source: `ModelingAlgorithms/TKFillet/Blend/Blend_Point.hxx:31-308`,
//! `Blend_Point.lxx:19-200`, `Blend_Point.cxx:21-483` (OCCT 8.0.0).
//!
//! The nine OCCT constructors and the ten `SetValue` overloads all land on
//! [`BlendPoint::set_value_*`] methods: the class keeps one flat set of members
//! and C++ overloading only chooses which of them to write, so one named method
//! per overload keeps the control flow identical without an enum discriminant.
//! The `hass1` / `hass2` / `hasc1` / `hasc2` / `istgt` flags each overload
//! leaves are preserved exactly as written in `.cxx:21-483` (they are *not*
//! uniformly "set everything").
//!
//! Divergences that are not heuristics:
//! - `Standard_DomainError` from the accessors that guard a missing tangent or
//!   parameter (`lxx:39-197`) becomes `Err`,
//! - `gp_Vec2d` becomes [`GpVec2d`] and `gp_Pnt` / `gp_Vec` become [`GpPnt`] /
//!   [`GpVec`]; `Tangent2dOnS1` / `Tangent2dOnS2` / `Tangent2d` rebuild the 2d
//!   vector from the stored components, like `lxx:82-95`,
//! - `Blend_Point()` (`cxx:21-24`) leaves every member but `istgt` uninitialized
//!   in OCCT; the port starts the four flags at `false` and the numbers at `0.`
//!   (the `gp_Pnt` / `gp_Vec` members are `(0,0,0)` in both).

use occt_core::gp::{GpPnt, GpVec, GpVec2d};

/// `class Blend_Point` (`hxx:31-308`).
#[derive(Debug, Clone, Copy)]
pub struct BlendPoint {
    /// `pt1` (`hxx:283`)
    pt1: GpPnt,
    /// `pt2` (`hxx:284`)
    pt2: GpPnt,
    /// `tg1` (`hxx:285`)
    tg1: GpVec,
    /// `tg2` (`hxx:286`)
    tg2: GpVec,
    /// `prm` (`hxx:287`)
    prm: f64,
    /// `u1` (`hxx:288`)
    u1: f64,
    /// `v1` (`hxx:289`)
    v1: f64,
    /// `u2` (`hxx:290`)
    u2: f64,
    /// `v2` (`hxx:291`)
    v2: f64,
    /// `pc1` (`hxx:292`)
    pc1: f64,
    /// `pc2` (`hxx:293`)
    pc2: f64,
    /// `utg12d` (`hxx:294`)
    utg12d: f64,
    /// `vtg12d` (`hxx:295`)
    vtg12d: f64,
    /// `utg22d` (`hxx:296`)
    utg22d: f64,
    /// `vtg22d` (`hxx:297`)
    vtg22d: f64,
    /// `hass1` (`hxx:298`)
    hass1: bool,
    /// `hass2` (`hxx:299`)
    hass2: bool,
    /// `hasc1` (`hxx:300`)
    hasc1: bool,
    /// `hasc2` (`hxx:301`)
    hasc2: bool,
    /// `istgt` (`hxx:302`)
    istgt: bool,
}

impl Default for BlendPoint {
    fn default() -> Self {
        Self::new()
    }
}

impl BlendPoint {
    /// `Blend_Point()` (`cxx:21-24`): only `istgt` is assigned in OCCT.
    pub fn new() -> Self {
        Self {
            pt1: GpPnt::zero(),
            pt2: GpPnt::zero(),
            tg1: GpVec::zero(),
            tg2: GpVec::zero(),
            prm: 0.0,
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
            pc1: 0.0,
            pc2: 0.0,
            utg12d: 0.0,
            vtg12d: 0.0,
            utg22d: 0.0,
            vtg22d: 0.0,
            hass1: false,
            hass2: false,
            hasc1: false,
            hasc2: false,
            istgt: true,
        }
    }

    /// `Blend_Point()` + `SetValue(P1, P2, Param, U1, V1, U2, V2, Tg1, Tg2,
    /// Tg12d, Tg22d)` (`cxx:26-54`, `:79-109`): a point on two surfaces, with
    /// tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn on_two_surfaces_with_tangents(
        p1: GpPnt,
        p2: GpPnt,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tg1: GpVec,
        tg2: GpVec,
        tg12d: GpVec2d,
        tg22d: GpVec2d,
    ) -> Self {
        let mut point = Self::new();
        point.set_on_two_surfaces_with_tangents(
            p1, p2, param, u1, v1, u2, v2, tg1, tg2, tg12d, tg22d,
        );
        point
    }

    /// `Blend_Point()` + `SetValue(P1, P2, Param, U1, V1, U2, V2)`
    /// (`cxx:56-77`, `:111-131`): a point on two surfaces, without tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn on_two_surfaces(
        p1: GpPnt,
        p2: GpPnt,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) -> Self {
        let mut point = Self::new();
        point.set_on_two_surfaces(p1, p2, param, u1, v1, u2, v2);
        point
    }

    /// `Blend_Point()` + `SetValue(Ps, Pc, Param, U, V, W, Tgs, Tgc, Tg2d)`
    /// (`cxx:133-158`, `:180-205`): a point on a surface and a curve, with
    /// tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn on_surface_and_curve_with_tangents(
        ps: GpPnt,
        pc: GpPnt,
        param: f64,
        u: f64,
        v: f64,
        w: f64,
        tgs: GpVec,
        tgc: GpVec,
        tg2d: GpVec2d,
    ) -> Self {
        let mut point = Self::new();
        point.set_on_surface_and_curve_with_tangents(ps, pc, param, u, v, w, tgs, tgc, tg2d);
        point
    }

    /// `Blend_Point()` + `SetValue(Ps, Pc, Param, U, V, W)` (`cxx:160-178`,
    /// `:207-225`): a point on a surface and a curve, without tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn on_surface_and_curve(
        ps: GpPnt,
        pc: GpPnt,
        param: f64,
        u: f64,
        v: f64,
        w: f64,
    ) -> Self {
        let mut point = Self::new();
        point.set_on_surface_and_curve(ps, pc, param, u, v, w);
        point
    }

    /// `Blend_Point()` + `SetValue(P1, P2, Param, U1, V1, U2, V2, PC, Tg1, Tg2,
    /// Tg12d, Tg22d)` (`cxx:227-259`, `:285-317`): a point on a surface and a
    /// curve on surface, with tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn on_surface_and_pcurve_with_tangents(
        p1: GpPnt,
        p2: GpPnt,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pc: f64,
        tg1: GpVec,
        tg2: GpVec,
        tg12d: GpVec2d,
        tg22d: GpVec2d,
    ) -> Self {
        let mut point = Self::new();
        point.set_on_surface_and_pcurve_with_tangents(
            p1, p2, param, u1, v1, u2, v2, pc, tg1, tg2, tg12d, tg22d,
        );
        point
    }

    /// `Blend_Point()` + `SetValue(P1, P2, Param, U1, V1, U2, V2, PC)`
    /// (`cxx:261-283`, `:319-341`): a point on a surface and a curve on
    /// surface, without tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn on_surface_and_pcurve(
        p1: GpPnt,
        p2: GpPnt,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pc: f64,
    ) -> Self {
        let mut point = Self::new();
        point.set_on_surface_and_pcurve(p1, p2, param, u1, v1, u2, v2, pc);
        point
    }

    /// `Blend_Point()` + `SetValue(P1, P2, Param, U1, V1, U2, V2, PC1, PC2,
    /// Tg1, Tg2, Tg12d, Tg22d)` (`cxx:343-376`, `:403-439`): a point on two
    /// curves on surfaces, with tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn on_two_pcurves_with_tangents(
        p1: GpPnt,
        p2: GpPnt,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pc1: f64,
        pc2: f64,
        tg1: GpVec,
        tg2: GpVec,
        tg12d: GpVec2d,
        tg22d: GpVec2d,
    ) -> Self {
        let mut point = Self::new();
        point.set_on_two_pcurves_with_tangents(
            p1, p2, param, u1, v1, u2, v2, pc1, pc2, tg1, tg2, tg12d, tg22d,
        );
        point
    }

    /// `Blend_Point()` + `SetValue(P1, P2, Param, U1, V1, U2, V2, PC1, PC2)`
    /// (`cxx:378-401`, `:441-465`): a point on two curves on surfaces, without
    /// tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn on_two_pcurves(
        p1: GpPnt,
        p2: GpPnt,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pc1: f64,
        pc2: f64,
    ) -> Self {
        let mut point = Self::new();
        point.set_on_two_pcurves(p1, p2, param, u1, v1, u2, v2, pc1, pc2);
        point
    }

    /// `SetValue(P1, P2, Param, U1, V1, U2, V2, Tg1, Tg2, Tg12d, Tg22d)`
    /// (`cxx:79-109`).
    #[allow(clippy::too_many_arguments)]
    pub fn set_on_two_surfaces_with_tangents(
        &mut self,
        p1: GpPnt,
        p2: GpPnt,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tg1: GpVec,
        tg2: GpVec,
        tg12d: GpVec2d,
        tg22d: GpVec2d,
    ) {
        self.pt1 = p1;
        self.pt2 = p2;
        self.prm = param;
        self.u1 = u1;
        self.v1 = v1;
        self.hass1 = true;
        self.u2 = u2;
        self.v2 = v2;
        self.hass2 = true;
        self.hasc1 = false;
        self.hasc2 = false;
        self.istgt = false;
        self.tg1 = tg1;
        self.tg2 = tg2;
        self.utg12d = tg12d.coord.x;
        self.vtg12d = tg12d.coord.y;
        self.utg22d = tg22d.coord.x;
        self.vtg22d = tg22d.coord.y;
    }

    /// `SetValue(P1, P2, Param, U1, V1, U2, V2)` (`cxx:111-131`).
    #[allow(clippy::too_many_arguments)]
    pub fn set_on_two_surfaces(
        &mut self,
        p1: GpPnt,
        p2: GpPnt,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        self.pt1 = p1;
        self.pt2 = p2;
        self.prm = param;
        self.u1 = u1;
        self.v1 = v1;
        self.hass1 = true;
        self.u2 = u2;
        self.v2 = v2;
        self.hass2 = true;
        self.hasc1 = false;
        self.hasc2 = false;
        self.istgt = true;
    }

    /// `SetValue(Ps, Pc, Param, U, V, W, Tgs, Tgc, Tg2d)` (`cxx:180-205`).
    #[allow(clippy::too_many_arguments)]
    pub fn set_on_surface_and_curve_with_tangents(
        &mut self,
        ps: GpPnt,
        pc: GpPnt,
        param: f64,
        u: f64,
        v: f64,
        w: f64,
        tgs: GpVec,
        tgc: GpVec,
        tg2d: GpVec2d,
    ) {
        self.pt1 = ps;
        self.pt2 = pc;
        self.prm = param;
        self.u1 = u;
        self.v1 = v;
        self.hass1 = true;
        self.hass2 = false;
        self.hasc1 = false;
        self.pc2 = w;
        self.hasc2 = true;
        self.istgt = false;
        self.tg1 = tgs;
        self.tg2 = tgc;
        self.utg12d = tg2d.coord.x;
        self.vtg12d = tg2d.coord.y;
    }

    /// `SetValue(Ps, Pc, Param, U, V, W)` (`cxx:207-225`).
    #[allow(clippy::too_many_arguments)]
    pub fn set_on_surface_and_curve(
        &mut self,
        ps: GpPnt,
        pc: GpPnt,
        param: f64,
        u: f64,
        v: f64,
        w: f64,
    ) {
        self.pt1 = ps;
        self.pt2 = pc;
        self.prm = param;
        self.u1 = u;
        self.v1 = v;
        self.hass1 = true;
        self.hass2 = false;
        self.hasc1 = false;
        self.pc2 = w;
        self.hasc2 = true;
        self.istgt = true;
    }

    /// `SetValue(P1, P2, Param, U1, V1, U2, V2, PC, Tg1, Tg2, Tg12d, Tg22d)`
    /// (`cxx:285-317`).
    #[allow(clippy::too_many_arguments)]
    pub fn set_on_surface_and_pcurve_with_tangents(
        &mut self,
        p1: GpPnt,
        p2: GpPnt,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pc: f64,
        tg1: GpVec,
        tg2: GpVec,
        tg12d: GpVec2d,
        tg22d: GpVec2d,
    ) {
        self.pt1 = p1;
        self.pt2 = p2;
        self.prm = param;
        self.u1 = u1;
        self.v1 = v1;
        self.hass1 = true;
        self.u2 = u2;
        self.v2 = v2;
        self.hass2 = true;
        self.hasc1 = false;
        self.pc2 = pc;
        self.hasc2 = true;
        self.istgt = false;
        self.tg1 = tg1;
        self.tg2 = tg2;
        self.utg12d = tg12d.coord.x;
        self.vtg12d = tg12d.coord.y;
        self.utg22d = tg22d.coord.x;
        self.vtg22d = tg22d.coord.y;
    }

    /// `SetValue(P1, P2, Param, U1, V1, U2, V2, PC)` (`cxx:319-341`).
    #[allow(clippy::too_many_arguments)]
    pub fn set_on_surface_and_pcurve(
        &mut self,
        p1: GpPnt,
        p2: GpPnt,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pc: f64,
    ) {
        self.pt1 = p1;
        self.pt2 = p2;
        self.prm = param;
        self.u1 = u1;
        self.v1 = v1;
        self.hass1 = true;
        self.u2 = u2;
        self.v2 = v2;
        self.hass2 = true;
        self.hasc1 = false;
        self.pc2 = pc;
        self.hasc2 = true;
        self.istgt = true;
    }

    /// `SetValue(P1, P2, Param, U1, V1, U2, V2, PC1, PC2, Tg1, Tg2, Tg12d,
    /// Tg22d)` (`cxx:403-439`).
    #[allow(clippy::too_many_arguments)]
    pub fn set_on_two_pcurves_with_tangents(
        &mut self,
        p1: GpPnt,
        p2: GpPnt,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pc1: f64,
        pc2: f64,
        tg1: GpVec,
        tg2: GpVec,
        tg12d: GpVec2d,
        tg22d: GpVec2d,
    ) {
        self.pt1 = p1;
        self.pt2 = p2;
        self.prm = param;
        self.u1 = u1;
        self.v1 = v1;
        self.hass1 = true;
        self.u2 = u2;
        self.v2 = v2;
        self.hass2 = true;
        self.pc1 = pc1;
        self.hasc1 = true;
        self.pc2 = pc2;
        self.hasc2 = true;
        self.istgt = false;
        self.tg1 = tg1;
        self.tg2 = tg2;
        self.utg12d = tg12d.coord.x;
        self.vtg12d = tg12d.coord.y;
        self.utg22d = tg22d.coord.x;
        self.vtg22d = tg22d.coord.y;
    }

    /// `SetValue(P1, P2, Param, U1, V1, U2, V2, PC1, PC2)` (`cxx:441-465`).
    #[allow(clippy::too_many_arguments)]
    pub fn set_on_two_pcurves(
        &mut self,
        p1: GpPnt,
        p2: GpPnt,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pc1: f64,
        pc2: f64,
    ) {
        self.pt1 = p1;
        self.pt2 = p2;
        self.prm = param;
        self.u1 = u1;
        self.v1 = v1;
        self.hass1 = true;
        self.u2 = u2;
        self.v2 = v2;
        self.hass2 = true;
        self.pc1 = pc1;
        self.hasc1 = true;
        self.pc2 = pc2;
        self.hasc2 = true;
        self.istgt = true;
    }

    /// `SetValue(P1, P2, Param, PC1, PC2)` (`cxx:467-483`): the only overload
    /// that leaves `hass1` / `hass2` false and writes neither `u`/`v` pair.
    pub fn set_on_two_curves(&mut self, p1: GpPnt, p2: GpPnt, param: f64, pc1: f64, pc2: f64) {
        self.pt1 = p1;
        self.pt2 = p2;
        self.prm = param;
        self.hass1 = false;
        self.hass2 = false;
        self.pc1 = pc1;
        self.hasc1 = true;
        self.pc2 = pc2;
        self.hasc2 = true;
        self.istgt = true;
    }

    /// `SetParameter(Param)` (`lxx:19-21`).
    pub fn set_parameter(&mut self, param: f64) {
        self.prm = param;
    }

    /// `Parameter()` (`lxx:32-35`).
    pub fn parameter(&self) -> f64 {
        self.prm
    }

    /// `PointOnS1()` (`lxx:23-26`).
    pub fn point_on_s1(&self) -> &GpPnt {
        &self.pt1
    }

    /// `PointOnS2()` (`lxx:28-30`).
    pub fn point_on_s2(&self) -> &GpPnt {
        &self.pt2
    }

    /// `ParametersOnS1(U1, V1)` (`lxx:37-45`): `Standard_DomainError` when the
    /// point was not built on a surface.
    pub fn parameters_on_s1(&self) -> Result<(f64, f64), String> {
        if !self.hass1 {
            return Err("Blend_Point::ParametersOnS1: Standard_DomainError".to_string());
        }
        Ok((self.u1, self.v1))
    }

    /// `ParametersOnS2(U2, V2)` (`lxx:47-55`).
    pub fn parameters_on_s2(&self) -> Result<(f64, f64), String> {
        if !self.hass2 {
            return Err("Blend_Point::ParametersOnS2: Standard_DomainError".to_string());
        }
        Ok((self.u2, self.v2))
    }

    /// `IsTangencyPoint()` (`lxx:57-60`).
    pub fn is_tangency_point(&self) -> bool {
        self.istgt
    }

    /// `TangentOnS1()` (`lxx:62-71`).
    pub fn tangent_on_s1(&self) -> Result<&GpVec, String> {
        if self.istgt {
            return Err("Blend_Point::TangentOnS1: Standard_DomainError".to_string());
        }
        Ok(&self.tg1)
    }

    /// `TangentOnS2()` (`lxx:73-79`).
    pub fn tangent_on_s2(&self) -> Result<&GpVec, String> {
        if self.istgt {
            return Err("Blend_Point::TangentOnS2: Standard_DomainError".to_string());
        }
        Ok(&self.tg2)
    }

    /// `Tangent2dOnS1()` (`lxx:82-89`).
    pub fn tangent2d_on_s1(&self) -> Result<GpVec2d, String> {
        if self.istgt || !self.hass1 {
            return Err("Blend_Point::Tangent2dOnS1: Standard_DomainError".to_string());
        }
        Ok(GpVec2d::new(self.utg12d, self.vtg12d))
    }

    /// `Tangent2dOnS2()` (`lxx:91-96`).
    pub fn tangent2d_on_s2(&self) -> Result<GpVec2d, String> {
        if self.istgt || !self.hass2 {
            return Err("Blend_Point::Tangent2dOnS2: Standard_DomainError".to_string());
        }
        Ok(GpVec2d::new(self.utg22d, self.vtg22d))
    }

    /// `PointOnS()` (`lxx:98-101`).
    pub fn point_on_s(&self) -> &GpPnt {
        &self.pt1
    }

    /// `PointOnC()` (`lxx:103-106`).
    pub fn point_on_c(&self) -> &GpPnt {
        &self.pt2
    }

    /// `ParametersOnS(U1, V1)` (`lxx:108-116`): the same guard as
    /// `ParametersOnS1`.
    pub fn parameters_on_s(&self) -> Result<(f64, f64), String> {
        if !self.hass1 {
            return Err("Blend_Point::ParametersOnS: Standard_DomainError".to_string());
        }
        Ok((self.u1, self.v1))
    }

    /// `ParameterOnC()` (`lxx:118-123`): guards `hasc2`, not `hass2`.
    pub fn parameter_on_c(&self) -> Result<f64, String> {
        if !self.hasc2 {
            return Err("Blend_Point::ParameterOnC: Standard_DomainError".to_string());
        }
        Ok(self.pc2)
    }

    /// `TangentOnS()` (`lxx:125-134`): `istgt || !hass1`.
    pub fn tangent_on_s(&self) -> Result<&GpVec, String> {
        if self.istgt || !self.hass1 {
            return Err("Blend_Point::TangentOnS: Standard_DomainError".to_string());
        }
        Ok(&self.tg1)
    }

    /// `TangentOnC()` (`lxx:136-143`): `istgt` only.
    pub fn tangent_on_c(&self) -> Result<&GpVec, String> {
        if self.istgt {
            return Err("Blend_Point::TangentOnC: Standard_DomainError".to_string());
        }
        Ok(&self.tg2)
    }

    /// `Tangent2d()` (`lxx:145-152`): `istgt || !hass1`.
    pub fn tangent2d(&self) -> Result<GpVec2d, String> {
        if self.istgt || !self.hass1 {
            return Err("Blend_Point::Tangent2d: Standard_DomainError".to_string());
        }
        Ok(GpVec2d::new(self.utg12d, self.vtg12d))
    }

    /// `PointOnC1()` (`lxx:154-157`).
    pub fn point_on_c1(&self) -> &GpPnt {
        &self.pt1
    }

    /// `PointOnC2()` (`lxx:159-162`).
    pub fn point_on_c2(&self) -> &GpPnt {
        &self.pt2
    }

    /// `ParameterOnC1()` (`lxx:164-171`).
    pub fn parameter_on_c1(&self) -> Result<f64, String> {
        if !self.hasc1 {
            return Err("Blend_Point::ParameterOnC1: Standard_DomainError".to_string());
        }
        Ok(self.pc1)
    }

    /// `ParameterOnC2()` (`lxx:173-182`).
    pub fn parameter_on_c2(&self) -> Result<f64, String> {
        if !self.hasc2 {
            return Err("Blend_Point::ParameterOnC2: Standard_DomainError".to_string());
        }
        Ok(self.pc2)
    }

    /// `TangentOnC1()` (`lxx:184-190`): `istgt || !hass1`.
    pub fn tangent_on_c1(&self) -> Result<&GpVec, String> {
        if self.istgt || !self.hass1 {
            return Err("Blend_Point::TangentOnC1: Standard_DomainError".to_string());
        }
        Ok(&self.tg1)
    }

    /// `TangentOnC2()` (`lxx:192-198`): `istgt` only.
    pub fn tangent_on_c2(&self) -> Result<&GpVec, String> {
        if self.istgt {
            return Err("Blend_Point::TangentOnC2: Standard_DomainError".to_string());
        }
        Ok(&self.tg2)
    }
}
