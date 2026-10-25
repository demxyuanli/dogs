//! `IntSurf_Transition`, `IntSurf_TypeTrans` and `IntSurf_Situation` — the
//! transition of an intersection line relative to a restriction curve.
//!
//! Source: `src/ModelingAlgorithms/TKGeomAlgo/IntSurf/IntSurf_Transition.hxx:31-92`,
//! `.lxx:16-74`, `.cxx:17-44`; `IntSurf_TypeTrans.hxx:19-26`;
//! `IntSurf_Situation.hxx:19-24`; `IntSurf.cxx:28-74`.
//!
//! The three accessors raise `Standard_DomainError` outside their valid
//! transition type (`lxx:47-53`, `:57-64`, `:67-73`); the port returns `Err`.

use occt_core::gp::{GpDir, GpVec};
use occt_core::precision::ANGULAR;
use occt_core::precision::CONFUSION;

/// `enum IntSurf_TypeTrans` (`TypeTrans.hxx:19-26`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum IntSurfTypeTrans {
    /// `IntSurf_In`
    In,
    /// `IntSurf_Out`
    Out,
    /// `IntSurf_Touch`
    Touch,
    /// `IntSurf_Undecided`
    #[default]
    Undecided,
}

/// `enum IntSurf_Situation` (`Situation.hxx:19-24`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum IntSurfSituation {
    /// `IntSurf_Inside`
    Inside,
    /// `IntSurf_Outside`
    Outside,
    /// `IntSurf_Unknown`
    #[default]
    Unknown,
}

/// `class IntSurf_Transition` (`hxx:31-92`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntSurfTransition {
    /// `tangent` (`hxx:88`)
    tangent: bool,
    /// `typetra` (`hxx:89`)
    typetra: IntSurfTypeTrans,
    /// `situat` (`hxx:90`)
    situat: IntSurfSituation,
    /// `oppos` (`hxx:91`)
    oppos: bool,
}

impl Default for IntSurfTransition {
    /// `IntSurf_Transition()` (`cxx:36-44`): an UNDECIDED transition.
    fn default() -> Self {
        Self::new()
    }
}

impl IntSurfTransition {
    /// `IntSurf_Transition()` (`cxx:36-44`).
    pub fn new() -> Self {
        Self {
            tangent: false,
            typetra: IntSurfTypeTrans::Undecided,
            situat: IntSurfSituation::Unknown,
            oppos: false,
        }
    }

    /// `IntSurf_Transition(Tangent, Type)` (`cxx:17-25`): an IN or OUT
    /// transition.
    pub fn with_type(tangent: bool, type_trans: IntSurfTypeTrans) -> Self {
        Self {
            tangent,
            typetra: type_trans,
            situat: IntSurfSituation::Unknown,
            oppos: false,
        }
    }

    /// `IntSurf_Transition(Tangent, Situ, Oppos)` (`cxx:27-34`): a TOUCH
    /// transition.
    pub fn with_situation(tangent: bool, situ: IntSurfSituation, oppos: bool) -> Self {
        Self {
            tangent,
            typetra: IntSurfTypeTrans::Touch,
            situat: situ,
            oppos,
        }
    }

    /// `SetValue(Tangent, Type)` (`lxx:16-21`).
    pub fn set_value_with_type(&mut self, tangent: bool, type_trans: IntSurfTypeTrans) {
        self.tangent = tangent;
        self.typetra = type_trans;
    }

    /// `SetValue(Tangent, Situ, Oppos)` (`lxx:23-33`).
    pub fn set_value_with_situation(&mut self, tangent: bool, situ: IntSurfSituation, oppos: bool) {
        self.tangent = tangent;
        self.typetra = IntSurfTypeTrans::Touch;
        self.situat = situ;
        self.oppos = oppos;
    }

    /// `SetValue()` (`lxx:35-40`): an UNDECIDED transition. Only `typetra` is
    /// written, as in OCCT.
    pub fn set_undecided(&mut self) {
        self.typetra = IntSurfTypeTrans::Undecided;
    }

    /// `TransitionType()` (`lxx:42-47`).
    pub fn transition_type(&self) -> IntSurfTypeTrans {
        self.typetra
    }

    /// `IsTangent()` (`lxx:49-55`).
    pub fn is_tangent(&self) -> Result<bool, String> {
        if self.typetra == IntSurfTypeTrans::Undecided {
            return Err("IntSurf_Transition::IsTangent: Standard_DomainError".to_string());
        }
        Ok(self.tangent)
    }

    /// `Situation()` (`lxx:57-64`).
    pub fn situation(&self) -> Result<IntSurfSituation, String> {
        if self.typetra != IntSurfTypeTrans::Touch {
            return Err("IntSurf_Transition::Situation: Standard_DomainError".to_string());
        }
        Ok(self.situat)
    }

    /// `IsOpposite()` (`lxx:66-73`).
    pub fn is_opposite(&self) -> Result<bool, String> {
        if self.typetra != IntSurfTypeTrans::Touch {
            return Err("IntSurf_Transition::IsOpposite: Standard_DomainError".to_string());
        }
        Ok(self.oppos)
    }
}

/// `IntSurf::MakeTransition(TgFirst, TgSecond, Normale, TFirst, TSecond)`
/// (`IntSurf.cxx:28-74`): the transition of the line and of the restriction
/// from the mixed product of the normal and the two tangents.
///
/// The `0.0001` threshold and the `Precision::Angular()` guard are the ones in
/// `cxx:45-71`; they are not parameters of the algorithm.
pub fn make_transition(
    tg_first: &GpVec,
    tg_second: &GpVec,
    normale: &GpDir,
    t_first: &mut IntSurfTransition,
    t_second: &mut IntSurfTransition,
) {
    // `gp_Vec pvect(TgSecond.Crossed(TgFirst))` (`cxx:39`).
    let pvect = tg_second.crossed(tg_first);
    let n_tg_second = tg_second.magnitude();
    let n_tg_first = tg_first.magnitude();
    let n_tg_second_n_tg_first_angular = n_tg_second * n_tg_first * ANGULAR;

    if n_tg_first <= CONFUSION {
        t_first.set_value_with_type(true, IntSurfTypeTrans::Undecided);
        t_second.set_value_with_type(true, IntSurfTypeTrans::Undecided);
    } else if n_tg_second <= CONFUSION || pvect.magnitude() <= n_tg_second_n_tg_first_angular {
        let oppos = tg_first.dot(tg_second) < 0.0;
        t_first.set_value_with_situation(true, IntSurfSituation::Unknown, oppos);
        t_second.set_value_with_situation(true, IntSurfSituation::Unknown, oppos);
    } else {
        let mut yu = pvect.dot(&GpVec::from_xyz(normale.xyz()));
        yu /= n_tg_second * n_tg_first;
        if yu > 0.0001 {
            t_first.set_value_with_type(false, IntSurfTypeTrans::In);
            t_second.set_value_with_type(false, IntSurfTypeTrans::Out);
        } else if yu < -0.0001 {
            t_first.set_value_with_type(false, IntSurfTypeTrans::Out);
            t_second.set_value_with_type(false, IntSurfTypeTrans::In);
        } else {
            t_first.set_value_with_type(true, IntSurfTypeTrans::Undecided);
            t_second.set_value_with_type(true, IntSurfTypeTrans::Undecided);
        }
    }
}
