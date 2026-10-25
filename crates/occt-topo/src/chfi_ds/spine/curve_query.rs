//! The adaptor queries of `ChFiDS_Spine`: abscissa, parameter, `D0`-`D2`,
//! and the current analytic curve.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_Spine.cxx:537-932` and
//! the `const_cast` accesses to `myCurve` (`:548-552`, `:584-590`, `:727-734`,
//! `:776-781`, `:826-831`), OCCT 8.0.0.

use super::*;

impl Spine {
    /// `Absc(U)` (`cxx:537-540`): the abscissa of parameter `U` on the current
    /// elementary spine.
    pub fn absc(&mut self, u: f64) -> Option<f64> {
        let index = self.index_of_curve;
        self.absc_on(u, index)
    }

    /// `Absc(U, I)` (`cxx:544-563`), `I` 1-based.
    pub fn absc_on(&mut self, u: f64, i: i32) -> Option<f64> {
        self.ensure_current(i)?;
        let curve: &dyn Curve = self.my_curve.as_ref()?;
        let mut l = self.first_parameter_of(i)?;
        let reversed = self
            .spine
            .get(usize::try_from(i - 1).ok()?)?
            .orientation()
            .is_reversed();
        if reversed {
            l += curve_length_between(curve, u, curve.last_parameter(), None);
        } else {
            l += curve_length_between(curve, curve.first_parameter(), u, None);
        }
        Some(l)
    }

    /// `Parameter(AbsC, U, Oriented = true)` (`cxx:567-578`).
    pub fn parameter(&mut self, abs_c: f64, oriented: bool) -> Option<f64> {
        let len = self.abscissa_len()?;
        let mut index = len;
        for i in 1..len {
            if abs_c < self.abscissa_value(i)? {
                index = i;
                break;
            }
        }
        self.parameter_of(index, abs_c, oriented)
    }

    /// `Parameter(Index, AbsC, U, Oriented = true)` (`cxx:582-612`).
    ///
    /// `GCPnts_AbscissaPoint(myCurve, L, myCurve.FirstParameter(), uapp)`
    /// (`cxx:609`) is the four-argument constructor
    /// (`GCPnts_AbscissaPoint.cxx:511-517`), so the caller's `uapp` is used as
    /// the initial guess and the epsilon is the adaptor's
    /// `Resolution(Precision::Confusion())`.
    pub fn parameter_of(&mut self, index: i32, abs_c: f64, oriented: bool) -> Option<f64> {
        self.ensure_current(index)?;
        let curve: &dyn Curve = self.my_curve.as_ref()?;
        let reversed = self
            .spine
            .get(usize::try_from(index - 1).ok()?)?
            .orientation()
            .is_reversed();
        let l = if reversed {
            self.abscissa_value(self.index_of_curve)? - abs_c
        } else if self.index_of_curve == 1 {
            abs_c
        } else {
            abs_c - self.abscissa_value(self.index_of_curve - 1)?
        };
        let t = l / self.length_of(index)?;
        let uapp = (1.0 - t) * curve.first_parameter() + t * curve.last_parameter();
        let mut u = abscissa_point_with_guess(curve, l, curve.first_parameter(), uapp)?;
        if reversed && oriented {
            u = (curve.last_parameter() + curve.first_parameter()) - u;
        }
        Some(u)
    }

    /// `Value(AbsC)` (`cxx:702-739`).
    pub fn value(&mut self, abs_c: f64) -> Option<GpPnt> {
        let mut l = abs_c;
        let index = self.prepare(&mut l)?;

        if index == -1 {
            // `cxx:711-717`.
            let mut p = self.firstori;
            let vp = self.firsttgt.multiplied_scalar(l);
            p.translate_vec(&vp);
            return Some(p);
        }
        if index == self.abscissa_len()? + 1 {
            // `cxx:719-726`.
            let mut p = self.lastori;
            let vp = self.lasttgt.multiplied_scalar(l);
            p.translate_vec(&vp);
            return Some(p);
        }
        self.ensure_current(index)?;
        let curve: &dyn Curve = self.my_curve.as_ref()?;
        let t = l / self.length_of(index)?;
        let uapp = (1.0 - t) * curve.first_parameter() + t * curve.last_parameter();
        let u = abscissa_point_with_guess(curve, l, curve.first_parameter(), uapp)?;
        Some(curve.d0(u))
    }

    /// `D0(AbsC, P)` (`cxx:744-747`).
    pub fn d0(&mut self, abs_c: f64) -> Option<GpPnt> {
        self.value(abs_c)
    }

    /// `D1(AbsC, P, V1)` (`cxx:751-793`).
    pub fn d1(&mut self, abs_c: f64) -> Option<(GpPnt, GpVec)> {
        let mut l = abs_c;
        let index = self.prepare(&mut l)?;
        let len = self.abscissa_len()?;

        if index == -1 {
            // `cxx:758-764`: the tangent itself is *not* normalized here.
            let mut p = self.firstori;
            let v1 = self.firsttgt;
            let vp = self.firsttgt.multiplied_scalar(l);
            p.translate_vec(&vp);
            return Some((p, v1));
        }
        if index == len + 1 {
            // `cxx:766-773`.
            let mut p = self.lastori;
            let v1 = self.lasttgt;
            let vp = self.lasttgt.multiplied_scalar(l);
            p.translate_vec(&vp);
            return Some((p, v1));
        }
        self.ensure_current(index)?;
        let curve: &dyn Curve = self.my_curve.as_ref()?;
        let t = l / self.length_of(index)?;
        let uapp = (1.0 - t) * curve.first_parameter() + t * curve.last_parameter();
        let u = abscissa_point_with_guess(curve, l, curve.first_parameter(), uapp)?;
        let (p, mut v1) = curve.d1(u);
        let mut d1 = 1.0 / v1.magnitude();
        if self
            .spine
            .get(usize::try_from(index - 1).ok()?)?
            .orientation()
            .is_reversed()
        {
            d1 = -d1;
        }
        v1 = v1.multiplied_scalar(d1);
        Some((p, v1))
    }

    /// `D2(AbsC, P, V1, V2)` (`cxx:797-852`).
    pub fn d2(&mut self, abs_c: f64) -> Option<(GpPnt, GpVec, GpVec)> {
        let mut l = abs_c;
        let index = self.prepare(&mut l)?;
        let len = self.abscissa_len()?;

        if index == -1 {
            // `cxx:807-814`.
            let mut p = self.firstori;
            let v1 = self.firsttgt;
            let v2 = GpVec::zero();
            let vp = self.firsttgt.multiplied_scalar(l);
            p.translate_vec(&vp);
            return Some((p, v1, v2));
        }
        if index == len + 1 {
            // `cxx:816-824`.
            let mut p = self.lastori;
            let v1 = self.lasttgt;
            let v2 = GpVec::zero();
            let vp = self.lasttgt.multiplied_scalar(l);
            p.translate_vec(&vp);
            return Some((p, v1, v2));
        }
        self.ensure_current(index)?;
        let curve: &dyn Curve = self.my_curve.as_ref()?;
        let t = l / self.length_of(index)?;
        let uapp = (1.0 - t) * curve.first_parameter() + t * curve.last_parameter();
        let u = abscissa_point_with_guess(curve, l, curve.first_parameter(), uapp)?;
        let (p, mut v1, mut v2) = curve.d2(u);
        // `cxx:838-849`, verbatim: the arc-length reparametrisation of D2.
        let mut n1 = v1.square_magnitude();
        let d2 = -(v1.dot(&v2)) * (1.0 / n1) * (1.0 / n1);
        v2 = v2.multiplied_scalar(1.0 / n1);
        n1 = n1.sqrt();
        let va = v1.multiplied_scalar(d2);
        v2 = v2.added(&va);
        let mut d1 = 1.0 / n1;
        if self
            .spine
            .get(usize::try_from(index - 1).ok()?)?
            .orientation()
            .is_reversed()
        {
            d1 = -d1;
        }
        v1 = v1.multiplied_scalar(d1);
        Some((p, v1, v2))
    }

    /// `SetCurrent(Index)` (`cxx:856-862`).
    pub fn set_current(&mut self, index: i32) -> Option<()> {
        self.ensure_current(index)
    }

    /// `CurrentElementarySpine(Index)` (`cxx:866-876`) — returns `myCurve`
    /// after switching it to `Index`. OCCT hands back a const reference to its
    /// single adaptor; the port hands back a copy of it.
    pub fn current_elementary_spine(&mut self, index: i32) -> Option<CurveAdaptor> {
        self.ensure_current(index)?;
        self.my_curve.clone()
    }

    /// `CurrentIndexOfElementarySpine()` (`lxx:152-155`).
    pub fn current_index_of_elementary_spine(&self) -> i32 {
        self.index_of_curve
    }

    /// `GetType()` (`cxx:879-882`): the current `BRepAdaptor_Curve` type.
    pub fn get_type(&self) -> Option<GeomAbsCurveType> {
        self.my_curve.as_ref().map(|c| c.get_type())
    }

    /// `Line()` (`cxx:885-897`): the current edge's line, with the origin moved
    /// to the end the orientation starts from, and the direction reversed for a
    /// `REVERSED` edge.
    pub fn line(&self) -> Option<GpLin> {
        let curve: &dyn Curve = self.my_curve.as_ref()?;
        let mut ll = curve.gp_line()?;
        let (first, last) = (curve.first_parameter(), curve.last_parameter());
        let reversed = self
            .spine
            .get(usize::try_from(self.index_of_curve - 1).ok()?)?
            .orientation()
            .is_reversed();
        if reversed {
            // `gp_Lin::Reverse()` flips the axis direction, keeping its origin.
            let mut ax = *ll.position();
            ax.reverse();
            ll.set_position(ax);
            ll.set_location(curve.d0(last));
        } else {
            ll.set_location(curve.d0(first));
        }
        Some(ll)
    }

    /// `Circle()` (`cxx:901-915`): the current edge's circle re-axed so that
    /// its X direction points at the first point of the edge (at the last point
    /// and about a reversed Z for a `REVERSED` edge).
    pub fn circle(&self) -> Option<GpCirc> {
        let curve: &dyn Curve = self.my_curve.as_ref()?;
        let circ = curve.gp_circ()?;
        let ac = circ.position();
        let reversed = self
            .spine
            .get(usize::try_from(self.index_of_curve - 1).ok()?)?
            .orientation()
            .is_reversed();
        let p = if reversed {
            curve.d0(curve.last_parameter())
        } else {
            curve.d0(curve.first_parameter())
        };
        let dc = GpDir::from_vec(&GpVec::from_pnts(&ac.location(), &p)).ok()?;
        let mut zz = ac.direction();
        if reversed {
            zz.reverse();
        }
        let a = GpAx2::new(ac.location(), zz, dc).ok()?;
        Some(GpCirc::new(a, circ.radius()))
    }
}
