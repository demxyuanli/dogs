//! `HLRBRep_Algo` (`ModelingAlgorithms/TKHLR/HLRBRep/HLRBRep_Algo.hxx:27-92`,
//! `.cxx:29-98`): the `TopoDS_Shape`-flavoured facade over
//! [`HlrBRepInternalAlgo`]. Every `Add` wraps the shape in an
//! `HLRTopoBRep_OutLiner` and forwards to `Load`.
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `HLRBRep_Algo()` | `cxx:33` | [`HlrBRepAlgo::new`] |
//! | `HLRBRep_Algo(A)` | `cxx:37-41` | [`HlrBRepAlgo::from`] |
//! | `Add(S, SData, nbIso)` | `cxx:45-51` | [`HlrBRepAlgo::add_with_data`] |
//! | `Add(S, nbIso)` | `cxx:55-58` | [`HlrBRepAlgo::add`] |
//! | `Index(S)` | `cxx:62-78` | [`HlrBRepAlgo::index`] |
//! | `OutLinedShapeNullify()` | `cxx:82-92` | [`HlrBRepAlgo::out_lined_shape_nullify`] |
//!
//! # Translation notes
//!
//! * `ShapeBounds(i).Shape()->OriginalShape() == S` (`cxx:68`) and
//!   `OutLinedShape() == S` (`cxx:72`) compare `TopoDS_Shape` values; the port
//!   compares the ported `TopoShape` the same way
//!   ([`occt_topo::shape::TopoShape`] `PartialEq`).
//! * `ShapeBounds(i).Shape()->OutLinedShape(TopoDS_Shape())` (`cxx:90`) is
//!   `set_out_lined_shape(None)`; `DataStructure().Clear()` (`cxx:91`) is
//!   [`crate::hlrtopobrep::data::HlrTopoBRepData::clear`].
//! * `DEFINE_STANDARD_RTTIEXT` / `IMPLEMENT_STANDARD_RTTIEXT` have no port
//!   counterpart.

use std::any::Any;
use std::sync::Arc;

use crate::hlrbrep_internal_algo::HlrBRepInternalAlgo;
use crate::hlrtopobrep::out_liner::HlrTopoBRepOutLiner;
use occt_topo::shape::TopoShape;
use occt_topo::topo_tools::is_same;

/// `HLRBRep_Algo`.
pub struct HlrBRepAlgo {
    /// The `HLRBRep_InternalAlgo` base sub-object.
    base: HlrBRepInternalAlgo,
}

impl HlrBRepAlgo {
    /// `HLRBRep_Algo()` (`cxx:33`): `= default`, i.e. the base default
    /// constructor.
    pub fn new() -> Self {
        Self {
            base: HlrBRepInternalAlgo::new(),
        }
    }

    /// `HLRBRep_Algo(const handle(HLRBRep_Algo)& A)` (`cxx:37-41`): the base
    /// copy constructor.
    pub fn from(a: &HlrBRepAlgo) -> Self {
        Self {
            base: HlrBRepInternalAlgo::from(&a.base),
        }
    }

    /// `Add(S, SData, nbIso)` (`cxx:45-51`):
    /// `Load(new HLRTopoBRep_OutLiner(S), SData, nbIso)`.
    pub fn add_with_data(
        &mut self,
        s: &TopoShape,
        s_data: Option<Arc<dyn Any + Send + Sync>>,
        nb_iso: i32,
    ) {
        self.base.load_with_data(
            Arc::new(HlrTopoBRepOutLiner::with_original(s.clone())),
            s_data,
            nb_iso,
        );
    }

    /// `Add(S, nbIso = 0)` (`cxx:55-58`): `Load(new HLRTopoBRep_OutLiner(S), nbIso)`.
    pub fn add(&mut self, s: &TopoShape, nb_iso: i32) {
        self.base.load(
            Arc::new(HlrTopoBRepOutLiner::with_original(s.clone())),
            nb_iso,
        );
    }

    /// `Index(S)` (`cxx:62-78`): the 1-based index of the shape or of its
    /// outlined copy, `0` when neither is found.
    pub fn index(&self, s: &TopoShape) -> i32 {
        let n = self.base.nb_shapes();
        for i in 1..=n {
            let Some(sb) = self.base.shape_bounds(i) else {
                continue;
            };
            let Some(shape) = sb.shape() else {
                continue;
            };
            // `ShapeBounds(i).Shape()->OriginalShape() == S` (`cxx:68`):
            // `TopoDS_Shape::operator==` is `IsSame` (`TopoDS_Shape.hxx:254-258`).
            if shape
                .original_shape()
                .map(|o| is_same(o, s))
                .unwrap_or(false)
            {
                return i;
            }
            // `ShapeBounds(i).Shape()->OutLinedShape() == S` (`cxx:72`).
            if shape
                .out_lined_shape()
                .map(|o| is_same(o, s))
                .unwrap_or(false)
            {
                return i;
            }
        }
        0 // `cxx:77`
    }

    /// `OutLinedShapeNullify()` (`cxx:82-92`): clear the outlined shape and the
    /// outline data structure of every loaded shape.
    ///
    /// OCCT writes through the shared `handle<HLRTopoBRep_OutLiner>`
    /// (`cxx:90-91`), so every holder of that handle observes the clear. The
    /// port stores the liner behind `Arc` inside
    /// [`HlrBRepShapeBounds`]; the two writes are applied to a copy of the
    /// liner whose slot then takes it, which is the same state for the bounds
    /// entry (`HlrTopoBRepOutLiner` and `HlrTopoBRep_Data` are both `Clone`).
    pub fn out_lined_shape_nullify(&mut self) {
        let n = self.base.nb_shapes();
        for i in 1..=n {
            let Some(sb) = self.base.shape_bounds_mut(i) else {
                continue;
            };
            let Some(shape) = sb.shape() else {
                continue;
            };
            let mut liner = (**shape).clone();
            // `Shape()->OutLinedShape(TopoDS_Shape())` (`cxx:90`).
            liner.set_out_lined_shape(None);
            // `Shape()->DataStructure().Clear()` (`cxx:91`).
            liner.data_structure_mut().clear();
            sb.set_shape(Arc::new(liner));
        }
    }

    /// The `HLRBRep_InternalAlgo` base: `Projector`, `Update`, `Hide`,
    /// `NbShapes`, `ShapeBounds`, ...
    pub fn base(&self) -> &HlrBRepInternalAlgo {
        &self.base
    }

    /// The mutable `HLRBRep_InternalAlgo` base.
    pub fn base_mut(&mut self) -> &mut HlrBRepInternalAlgo {
        &mut self.base
    }
}

impl Default for HlrBRepAlgo {
    fn default() -> Self {
        Self::new()
    }
}
