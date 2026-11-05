//! `BRepFill_CompatibleWires` — put a sequence of wires into a state where
//! every wire has the same number of edges and a coherent origin/orientation,
//! so the surface ruled through them is not twisted.
//!
//! Source: `BRepFill_CompatibleWires.cxx:772-2629`, `.hxx:33-102`.
//!
//! Port status, keeping the OCCT control flow:
//!
//! * Ported: the data layer, both constructors (`cxx:772-786`), `Init`
//!   (`cxx:788-796`), `SetPercent` (`cxx:800-808`), `IsDone` (`cxx:810-813`),
//!   `Shape` (`cxx:817-820`), `GeneratedShapes` (`cxx:824-836`), `Generated`
//!   (`cxx:1001-1004`), `IsDegeneratedFirstSection` / `IsDegeneratedLastSection`
//!   (`cxx:841-851`), the analysis half of `Perform` (`cxx:855-958`: the
//!   FORWARD normalisation, `WireContinuity` per wire, the edge counts,
//!   `report`, the `myMap` initialisation, `myDegen1` / `myDegen2` and the
//!   `allClosed` / `allOpen` split), `WireContinuity` (`cxx:300-364`) and the
//!   nine file statics (`cxx:131-770`).
//! * Ported, each in its own module:
//!   - [`polar_method::same_number_by_polar_method`] (`cxx:1008-1503`)
//!   - [`acr_method::same_number_by_acr`] (`cxx:1505-1785`, with
//!     `BRepFill::InsertACR`)
//!   - [`compute_origin::compute_origin`] (`cxx:1787-2126`; the
//!     `OCCT_DEBUG_EFV` block `cxx:2128-2449` is not part of a default build)
//!   - [`search_origin::search_origin`] (`cxx:2453-2629`)
//!
//! Translation notes:
//! - `NCollection_Sequence<TopoDS_Shape> myInit / myWork` are `Vec<TopoShape>`.
//! - `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>` is a
//!   `HashMap<ShapeMapKey, Vec<TopoShape>>`; `Bind` overwrites, which is what
//!   `HashMap::insert` does.
//! - `myWork(i).Oriented(TopAbs_FORWARD)` (`cxx:872-874`) copies the wire with
//!   the FORWARD orientation and stores it back.
//! - `BRepFill_ThruSectionErrorStatus` is the enum already carried by
//!   [`crate::brep_fill_generator`], the other user of the same header.

pub mod acr_method;
pub mod compute_origin;
pub mod continuity;
pub mod perform;
pub mod polar_method;
pub mod search_origin;
pub mod statics;

use std::collections::HashMap;

use crate::brep_fill_generator::ThruSectionErrorStatus;
use crate::shape::TopoShape;
use crate::topo_tools::{shape_map_key, ShapeMapKey};

/// `BRepFill_CompatibleWires` (`BRepFill_CompatibleWires.hxx:33-102`).
pub struct CompatibleWires {
    /// `myInit` (`hxx:94`): the sections as given.
    pub(crate) my_init: Vec<TopoShape>,
    /// `myWork` (`hxx:95`): the sections being made compatible.
    pub(crate) my_work: Vec<TopoShape>,
    /// `myPercent` (`hxx:96`): the `SetPercent` cut ratio, `0.1` after `Init`.
    pub(crate) my_percent: f64,
    /// `myDegen1` (`hxx:97`): the first section is punctual.
    pub(crate) my_degen1: bool,
    /// `myDegen2` (`hxx:98`): the last section is punctual.
    pub(crate) my_degen2: bool,
    /// `myStatus` (`hxx:99`).
    pub(crate) my_status: ThruSectionErrorStatus,
    /// `myMap` (`hxx:100`): every section edge to the edges it was split into.
    pub(crate) my_map: HashMap<ShapeMapKey, Vec<TopoShape>>,
}

impl Default for CompatibleWires {
    fn default() -> Self {
        Self::new()
    }
}

impl CompatibleWires {
    /// `BRepFill_CompatibleWires()` (`cxx:772-776`): only `myStatus` is set,
    /// so `myPercent` and the two degeneracy flags keep their C++ defaults.
    pub fn new() -> Self {
        Self {
            my_init: Vec::new(),
            my_work: Vec::new(),
            my_percent: 0.1,
            my_degen1: false,
            my_degen2: false,
            my_status: ThruSectionErrorStatus::NotDone,
            my_map: HashMap::new(),
        }
    }

    /// `BRepFill_CompatibleWires(Sections)` (`cxx:780-786`).
    pub fn with_sections(sections: &[TopoShape]) -> Self {
        let mut me = Self::new();
        me.init(sections);
        me
    }

    /// `Init(Sections)` (`cxx:788-796`): keep the input, reset the work copy,
    /// the cut ratio, the status and the map.
    pub fn init(&mut self, sections: &[TopoShape]) {
        self.my_init = sections.to_vec();
        self.my_work = sections.to_vec();
        self.my_percent = 0.1;
        self.my_status = ThruSectionErrorStatus::NotDone;
        self.my_map.clear();
    }

    /// `SetPercent(Percent)` (`cxx:800-808`): ignored outside `(0, 1)`.
    pub fn set_percent(&mut self, percent: f64) {
        if 0.0 < percent && percent < 1.0 {
            self.my_percent = percent;
        }
    }

    /// `IsDone()` (`cxx:810-813`).
    pub fn is_done(&self) -> bool {
        self.my_status == ThruSectionErrorStatus::Done
    }

    /// `GetStatus()` (`hxx:54`).
    pub fn get_status(&self) -> ThruSectionErrorStatus {
        self.my_status
    }

    /// `Shape()` (`cxx:817-820`): the compatible wires.
    pub fn shape(&self) -> &[TopoShape] {
        &self.my_work
    }

    /// `GeneratedShapes(SubSection)` (`cxx:824-836`): an unknown edge answers
    /// the empty list, like OCCT's `static Empty`.
    pub fn generated_shapes(&self, sub_section: &TopoShape) -> &[TopoShape] {
        let key = shape_map_key(sub_section);
        self.my_map.get(&key).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// `Generated()` (`cxx:1001-1004`).
    pub fn generated(&self) -> &HashMap<ShapeMapKey, Vec<TopoShape>> {
        &self.my_map
    }

    /// `IsDegeneratedFirstSection()` (`cxx:841-845`).
    pub fn is_degenerated_first_section(&self) -> bool {
        self.my_degen1
    }

    /// `IsDegeneratedLastSection()` (`cxx:848-851`).
    pub fn is_degenerated_last_section(&self) -> bool {
        self.my_degen2
    }

    /// `Perform(WithRotation)` (`cxx:855-995`).
    pub fn perform(&mut self, with_rotation: bool) -> Result<(), String> {
        perform::perform(self, with_rotation)
    }
}
