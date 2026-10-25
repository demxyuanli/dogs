//! `TopOpeBRepDS_HDataStructure`'s `StoreInterference` / `StoreInterferences` /
//! `ClearStoreInterferences` family.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_HDataStructure.cxx:669-813`.
//!
//! `StoreInterference(I, LI)` appends `I` to `LI` and then, by
//! `I->GeometryType()`, to the DS list attached to the *geometry*
//! (`cxx:680-705`). The `SOLID` / `FACE` / `EDGE` / `VERTEX` and `POINT` arms
//! are commented out in OCCT, so they only append to `LI`.
//!
//! For the overloads whose `LI` *is* `myDS.ChangeShapeInterferences(...)`
//! (`cxx:720-813`), the single `LI.Append(I)` is the shape-list append; the
//! port therefore pushes to the shape list once and then runs only the
//! geometry switch.

use crate::shape::TopoShape;
use crate::topopebrep::ds::interference::InterferenceHandle;
use crate::topopebrep::ds::interference_iterator::InterferenceList;
use crate::topopebrep::ds::kind::Kind;

use super::HDataStructure;

impl HDataStructure {
    /// `StoreInterference(I, LI)` (`cxx:669-706`).
    ///
    /// `LI` is a plain list output here, not the DS list, so the geometry
    /// append is the only side effect on the DS.
    pub fn store_interference_in_list(&mut self, i: &InterferenceHandle, li: &mut InterferenceList) {
        // append I to list LI
        li.push(i.clone());

        // append I to list of interference connected to G = I->Geometry()
        self.store_geometry_interference(i);
    }

    /// `StoreInterference(I, S)` (`cxx:710-721`): `Standard_ProgramError` when
    /// `S` is not in the DS.
    pub fn store_interference_on_shape(
        &mut self,
        i: &InterferenceHandle,
        s: &TopoShape,
    ) -> Result<(), String> {
        let h = self.ds().has_shape(s, true);
        if !h {
            return Err("StoreInterference on shape out of DS".into());
        }
        // StoreInterference(I, myDS.ChangeShapeInterferences(S))
        self.change_ds().change_shape_interferences(s).push(i.clone());
        self.store_geometry_interference(i);
        Ok(())
    }

    /// `StoreInterference(I, IS)` (`cxx:725-737`): `Standard_ProgramError` when
    /// `IS` is outside `1..=NbShapes()`.
    pub fn store_interference_on_index(
        &mut self,
        i: &InterferenceHandle,
        is: i32,
    ) -> Result<(), String> {
        let n = self.ds().nb_shapes();
        if is < 1 || is > n {
            return Err("StoreInterference on index out of DS".into());
        }
        self.change_ds().change_shape_interferences_at(is).push(i.clone());
        self.store_geometry_interference(i);
        Ok(())
    }

    /// `StoreInterferences(LI, IS)` (`cxx:742-756`). No index check in OCCT.
    pub fn store_interferences_on_index(&mut self, li: &InterferenceList, is: i32) {
        for i in li {
            self.change_ds().change_shape_interferences_at(is).push(i.clone());
            self.store_geometry_interference(i);
        }
    }

    /// `StoreInterferences(LI, S)` (`cxx:760-774`).
    pub fn store_interferences_on_shape(&mut self, li: &InterferenceList, s: &TopoShape) {
        for i in li {
            self.change_ds().change_shape_interferences(s).push(i.clone());
            self.store_geometry_interference(i);
        }
    }

    /// `ClearStoreInterferences(LI, IS)` (`cxx:778-793`): the DS list is
    /// cleared first, then rebuilt from `LI`.
    pub fn clear_store_interferences_on_index(&mut self, li: &InterferenceList, is: i32) {
        self.change_ds().change_shape_interferences_at(is).clear();
        for i in li {
            self.change_ds().change_shape_interferences_at(is).push(i.clone());
            self.store_geometry_interference(i);
        }
    }

    /// `ClearStoreInterferences(LI, S)` (`cxx:798-813`).
    pub fn clear_store_interferences_on_shape(&mut self, li: &InterferenceList, s: &TopoShape) {
        self.change_ds().change_shape_interferences(s).clear();
        for i in li {
            self.change_ds().change_shape_interferences(s).push(i.clone());
            self.store_geometry_interference(i);
        }
    }

    /// `switch (I->GeometryType())` of `StoreInterference` (`cxx:680-705`):
    /// `SURFACE` and `CURVE` append to the geometry list, everything else is
    /// commented out.
    fn store_geometry_interference(&mut self, i: &InterferenceHandle) {
        let (g, gt) = {
            let h = i.read().expect("interference lock");
            (h.geometry(), h.geometry_type())
        };
        match gt {
            Kind::Solid | Kind::Face | Kind::Edge | Kind::Vertex => {}
            Kind::Surface => {
                self.change_ds().change_surface_interferences(g).push(i.clone());
            }
            Kind::Curve => {
                self.change_ds().change_curve_interferences(g).push(i.clone());
            }
            Kind::Point => {}
            _ => {}
        }
    }
}
