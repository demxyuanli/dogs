//! `Bisector_BisecCC::ChangeGuide` (`Bisector_BisecCC.cxx:399-449`).
//!
//! Returns a copy of the bisector re-parameterised along the *second* support
//! curve: the two curves, signs and convexity flags are swapped, and the
//! polygon is rebuilt with `ParamOnC1` / `ParamOnC2` exchanged and
//! `ParamOnBis` taken from the old `ParamOnC2`. The vertex order is reversed
//! when `sign1 == sign2`, because in that case the two parameterisations run in
//! opposite directions (`cxx:415-443`).
//!
//! `Bisector_Inter::NeighbourPerform` uses it to put two bisectors that share
//! one support curve onto the same parameter (`Bisector_Inter.cxx:352`).

use crate::bisector::{BisectorPointOnBis, BisectorPolyBis};

use super::BisectorBisecCC;

impl BisectorBisecCC {
    /// `ChangeGuide()` (`cxx:399-449`).
    pub fn change_guide(&self) -> BisectorBisecCC {
        let mut c = BisectorBisecCC::new();
        // `cxx:403-408`: `Curve(1, curve2)`, `Curve(2, curve1)`, `Sign(I, ...)`,
        // `IsConvex(I, ...)`.
        c.curve1 = self.curve2.clone();
        c.curve2 = self.curve1.clone();
        c.sign1 = self.sign2;
        c.sign2 = self.sign1;
        c.is_convex1 = self.is_convex2;
        c.is_convex2 = self.is_convex1;

        // `cxx:414-443`: rebuild the polygon.
        let mut poly = BisectorPolyBis::new();
        let length = self.my_polygon.length();
        if self.sign1 == self.sign2 {
            // `cxx:420-429`: same-sign bisectors reverse the vertex order.
            for i in (1..=length).rev() {
                poly.append(swapped_point(&self.my_polygon.value(i)));
            }
        } else {
            // `cxx:433-442`: opposite signs keep the order.
            for i in 1..=length {
                poly.append(swapped_point(&self.my_polygon.value(i)));
            }
        }
        c.my_polygon = poly;

        // `cxx:445-446`: `FirstParameter(Poly.First().ParamOnBis())` /
        // `LastParameter(Poly.Last().ParamOnBis())`. OCCT reads `First()` /
        // `Last()` of an empty polygon out of bounds; the port leaves the
        // interval sequences empty in that case.
        if !c.my_polygon.is_empty() {
            let first = c.my_polygon.first().param_on_bis();
            c.append_first_parameter(first);
            let last = c.my_polygon.last().param_on_bis();
            c.append_last_parameter(last);
        }
        c
    }
}

/// `Bisector_PointOnBis(P.ParamOnC2(), P.ParamOnC1(), P.ParamOnC2(),
/// P.Distance(), P.Point())` (`cxx:423-427`, `:436-440`).
fn swapped_point(p: &BisectorPointOnBis) -> BisectorPointOnBis {
    BisectorPointOnBis::with_params(
        p.param_on_c2(),
        p.param_on_c1(),
        p.param_on_c2(),
        p.distance(),
        p.point(),
    )
}
