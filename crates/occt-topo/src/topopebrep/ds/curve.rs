//! `TopOpeBRepDS_Curve` — a geometric curve, its tolerance, its range, the two
//! shapes it separates and the two surface/curve interferences.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_Curve.hxx:24-145`,
//! `TopOpeBRepDS_Curve.cxx:22-269`.

use std::sync::Arc;

use occt_core::precision::PCONFUSION;
use occt_geom::Curve as GeomCurve;
use occt_geom2d::curve::Curve2d;

use crate::shape::TopoShape;

use super::interference::InterferenceHandle;

/// Source: `TopOpeBRepDS_Curve.hxx:122-135`
#[derive(Clone, Default)]
pub struct Curve {
    curve: Option<Arc<dyn GeomCurve>>,
    first: f64,
    last: f64,
    range_defined: bool,
    tolerance: f64,
    is_walk: bool,
    s1: Option<TopoShape>,
    s2: Option<TopoShape>,
    sci1: Option<InterferenceHandle>,
    sci2: Option<InterferenceHandle>,
    keep: bool,
    mother: i32,
    ds_index: i32,
}

impl std::fmt::Debug for Curve {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TopOpeBRepDS_Curve")
            .field("is_null", &self.curve.is_none())
            .field("first", &self.first)
            .field("last", &self.last)
            .field("range_defined", &self.range_defined)
            .field("tolerance", &self.tolerance)
            .field("is_walk", &self.is_walk)
            .field("keep", &self.keep)
            .field("mother", &self.mother)
            .field("ds_index", &self.ds_index)
            .finish()
    }
}

impl Curve {
    /// `TopOpeBRepDS_Curve()` (`cxx:26-38`): a null curve, an undefined range,
    /// `myTolerance = Precision::Confusion()`, `myKeep = true`.
    pub fn new() -> Self {
        Self {
            curve: None,
            first: 0.0,
            last: 0.0,
            range_defined: false,
            tolerance: PCONFUSION,
            is_walk: false,
            s1: None,
            s2: None,
            sci1: None,
            sci2: None,
            keep: true,
            mother: 0,
            ds_index: 0,
        }
    }

    /// `TopOpeBRepDS_Curve(C, T, IsWalk)` (`cxx:42-52`): the tolerance and the
    /// walk flag come from `DefineCurve`, everything else from the member-init
    /// list (note that `myTolerance` is *not* initialized there).
    pub fn with_curve(curve: Option<Arc<dyn GeomCurve>>, t: f64, is_walk: bool) -> Self {
        let mut c = Self {
            curve: None,
            first: 0.0,
            last: 0.0,
            range_defined: false,
            tolerance: 0.0,
            is_walk: false,
            s1: None,
            s2: None,
            sci1: None,
            sci2: None,
            keep: true,
            mother: 0,
            ds_index: 0,
        };
        c.define_curve(curve, t, is_walk);
        c
    }

    /// `DefineCurve(C, T, IsWalk)` (`cxx:56-62`).
    pub fn define_curve(&mut self, curve: Option<Arc<dyn GeomCurve>>, t: f64, is_walk: bool) {
        self.curve = curve;
        self.tolerance = t;
        self.is_walk = is_walk;
    }

    /// `Tolerance(T)` (`cxx:66-69`).
    pub fn set_tolerance(&mut self, t: f64) {
        self.tolerance = t;
    }

    /// `SetSCI(SCI1, SCI2)` (`cxx:73-78`).
    pub fn set_sci(&mut self, sci1: Option<InterferenceHandle>, sci2: Option<InterferenceHandle>) {
        self.sci1 = sci1;
        self.sci2 = sci2;
    }

    /// `GetSCI1()` (`cxx:90-93`).
    pub fn get_sci1(&self) -> Option<&InterferenceHandle> {
        self.sci1.as_ref()
    }

    /// `GetSCI2()` (`cxx:97-100`).
    pub fn get_sci2(&self) -> Option<&InterferenceHandle> {
        self.sci2.as_ref()
    }

    /// `SetShapes(S1, S2)` (`cxx:110-114`).
    pub fn set_shapes(&mut self, s1: Option<TopoShape>, s2: Option<TopoShape>) {
        self.s1 = s1;
        self.s2 = s2;
    }

    /// `GetShapes(S1, S2)` (`cxx:118-122`).
    pub fn get_shapes(&self) -> (Option<&TopoShape>, Option<&TopoShape>) {
        (self.s1.as_ref(), self.s2.as_ref())
    }

    /// `Shape1()` (`cxx:126-129`).
    pub fn shape1(&self) -> Option<&TopoShape> {
        self.s1.as_ref()
    }

    /// `ChangeShape1()` (`cxx:133-136`).
    pub fn change_shape1(&mut self) -> &mut Option<TopoShape> {
        &mut self.s1
    }

    /// `Shape2()` (`cxx:140-143`).
    pub fn shape2(&self) -> Option<&TopoShape> {
        self.s2.as_ref()
    }

    /// `ChangeShape2()` (`cxx:147-150`).
    pub fn change_shape2(&mut self) -> &mut Option<TopoShape> {
        &mut self.s2
    }

    /// `Curve()` (`cxx:166-169`).
    pub fn curve(&self) -> Option<&Arc<dyn GeomCurve>> {
        self.curve.as_ref()
    }

    /// `ChangeCurve()` (`cxx:154-157`).
    pub fn change_curve(&mut self) -> &mut Option<Arc<dyn GeomCurve>> {
        &mut self.curve
    }

    /// `SetRange(First, Last)` (`cxx:173-178`).
    pub fn set_range(&mut self, first: f64, last: f64) {
        self.first = first;
        self.last = last;
        self.range_defined = true;
    }

    /// `Range(First, Last)` (`cxx:182-190`): returns `myRangeDefined` and
    /// writes the range only when it is defined.
    pub fn range(&self) -> Option<(f64, f64)> {
        if self.range_defined {
            Some((self.first, self.last))
        } else {
            None
        }
    }

    /// `Tolerance()` (`cxx:194-197`).
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// `Curve(C3D, Tol)` (`cxx:201-205`).
    pub fn set_curve(&mut self, c3d: Option<Arc<dyn GeomCurve>>, tol: f64) {
        self.curve = c3d;
        self.tolerance = tol;
    }

    /// `Curve1()` (`cxx:209-220`): the pcurve of `mySCI1` when it is a
    /// `SurfaceCurveInterference`, the null handle otherwise.
    pub fn curve1(&self) -> Option<Arc<dyn Curve2d>> {
        let sci1 = self.sci1.as_ref()?;
        let guard = match sci1.read() {
            Ok(g) => g,
            Err(p) => p.into_inner(),
        };
        guard.p_curve().cloned()
    }

    /// `Curve1(PC1)` (`cxx:224-230`): ignored when `mySCI1` is null or is not a
    /// `SurfaceCurveInterference` (the OCCT `down_cast` would be null).
    pub fn set_curve1(&mut self, pc1: Option<Arc<dyn Curve2d>>) {
        if let Some(sci1) = self.sci1.as_ref() {
            let mut guard = match sci1.write() {
                Ok(g) => g,
                Err(p) => p.into_inner(),
            };
            guard.set_p_curve(pc1);
        }
    }

    /// `Curve2()` (`cxx:234-245`).
    pub fn curve2(&self) -> Option<Arc<dyn Curve2d>> {
        let sci2 = self.sci2.as_ref()?;
        let guard = match sci2.read() {
            Ok(g) => g,
            Err(p) => p.into_inner(),
        };
        guard.p_curve().cloned()
    }

    /// `Curve2(PC2)` (`cxx:249-255`).
    pub fn set_curve2(&mut self, pc2: Option<Arc<dyn Curve2d>>) {
        if let Some(sci2) = self.sci2.as_ref() {
            let mut guard = match sci2.write() {
                Ok(g) => g,
                Err(p) => p.into_inner(),
            };
            guard.set_p_curve(pc2);
        }
    }

    /// `IsWalk()` (`cxx:259-262`).
    pub fn is_walk(&self) -> bool {
        self.is_walk
    }

    /// `ChangeIsWalk(B)` (`cxx:266-269`).
    pub fn change_is_walk(&mut self, b: bool) {
        self.is_walk = b;
    }

    /// `Keep()` (`cxx:273-276`).
    pub fn keep(&self) -> bool {
        self.keep
    }

    /// `ChangeKeep(B)` (`cxx:280-283`).
    pub fn change_keep(&mut self, b: bool) {
        self.keep = b;
    }

    /// `Mother()` (`cxx:287-290`).
    pub fn mother(&self) -> i32 {
        self.mother
    }

    /// `ChangeMother(I)` (`cxx:294-297`).
    pub fn change_mother(&mut self, i: i32) {
        self.mother = i;
    }

    /// `DSIndex()` (`cxx:301-304`).
    pub fn ds_index(&self) -> i32 {
        self.ds_index
    }

    /// `ChangeDSIndex(I)` (`cxx:308-311`).
    pub fn change_ds_index(&mut self, i: i32) {
        self.ds_index = i;
    }
}
