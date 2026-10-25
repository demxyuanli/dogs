//! `TopOpeBRepTool_connexity` — the edges a vertex is bound to, by orientation.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_connexity.cxx`
//! (OCCT 8.0.0, 244 lines): `Item` (`:87-93`), `AllItems` (`:96-108`),
//! `AddItem` (`:110-123`), `RemoveItem` (`:127-159`), `ChangeItem` (`:163-167`),
//! `IsMultiple` (`:171-181`), `IsFaulty` (`:185-191`) and `IsInternal`
//! (`:193-243`).
//!
//! `theItems` is an `NCollection_Array1` over 1..5 (`cxx:49`, `:57`) indexed by
//! the orientation codes of [`super::closed`]: FORWARD, REVERSED, INTERNAL,
//! EXTERNAL, CLOSING.

use crate::shape::TopoShape;
use crate::topo_tools::is_equal;
use crate::topopebrep::tool::closed::{
    ori_in_sor, CLOSING, EXTERNAL, FORWARD, INTERNAL, REVERSED,
};

/// `class TopOpeBRepTool_connexity` (`hxx:29-64`).
#[derive(Debug, Clone, Default)]
pub struct Connexity {
    /// `theKey` (`hxx:62`).
    the_key: Option<TopoShape>,
    /// `theItems` (`hxx:63`), indexed `1..=5` as in OCCT.
    the_items: [Vec<TopoShape>; 5],
}

impl Connexity {
    /// `TopOpeBRepTool_connexity()` (`cxx:47-50`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `TopOpeBRepTool_connexity(Key)` (`cxx:54-58`).
    pub fn with_key(key: &TopoShape) -> Self {
        Self {
            the_key: Some(key.clone()),
            the_items: Default::default(),
        }
    }

    /// `SetKey(Key)` (`cxx:63-66`).
    pub fn set_key(&mut self, key: &TopoShape) {
        self.the_key = Some(key.clone());
    }

    /// `Key()` (`cxx:71-74`).
    pub fn key(&self) -> Option<&TopoShape> {
        self.the_key.as_ref()
    }

    fn slot(ori_key: i32) -> Option<usize> {
        // `theItems(OriKey)` is 1-based (`cxx:49`).
        if (1..=5).contains(&ori_key) {
            Some((ori_key - 1) as usize)
        } else {
            None
        }
    }

    /// `Item(OriKey, Item)` (`cxx:87-93`).
    pub fn item(&self, ori_key: i32, item: &mut Vec<TopoShape>) -> usize {
        item.clear();
        if let Some(i) = Self::slot(ori_key) {
            item.extend(self.the_items[i].iter().cloned());
        }
        item.len()
    }

    /// `AllItems(Item)` (`cxx:96-108`): slots `1..4` (CLOSING excluded).
    pub fn all_items(&self, item: &mut Vec<TopoShape>) -> usize {
        item.clear();
        for i in 0..4 {
            item.extend(self.the_items[i].iter().cloned());
        }
        item.len()
    }

    /// `AddItem(OriKey, Item)` (`cxx:110-116`).
    pub fn add_items(&mut self, ori_key: i32, item: &[TopoShape]) {
        if let Some(i) = Self::slot(ori_key) {
            self.the_items[i].extend(item.iter().cloned());
        }
    }

    /// `AddItem(OriKey, Item)` (`cxx:118-123`).
    pub fn add_item(&mut self, ori_key: i32, item: &TopoShape) {
        if let Some(i) = Self::slot(ori_key) {
            self.the_items[i].push(item.clone());
        }
    }

    /// `RemoveItem(OriKey, Item)` (`cxx:127-142`): `TopoDS_Shape::IsEqual`.
    pub fn remove_item_ori(&mut self, ori_key: i32, item: &TopoShape) -> bool {
        let Some(i) = Self::slot(ori_key) else {
            return false;
        };
        match self.the_items[i].iter().position(|x| is_equal(x, item)) {
            Some(pos) => {
                self.the_items[i].remove(pos);
                true
            }
            None => false,
        }
    }

    /// `RemoveItem(Item)` (`cxx:146-159`): every slot `1..5`.
    pub fn remove_item(&mut self, item: &TopoShape) -> bool {
        let mut removed = false;
        for i in 1..=5 {
            if self.remove_item_ori(i, item) {
                removed = true;
            }
        }
        removed
    }

    /// `ChangeItem(OriKey)` (`cxx:163-167`); `None` for an out-of-range key,
    /// which OCCT would index out of bounds.
    pub fn change_item(&mut self, ori_key: i32) -> Option<&mut Vec<TopoShape>> {
        Self::slot(ori_key).map(|i| &mut self.the_items[i])
    }

    /// `IsMultiple()` (`cxx:171-181`): more than one of FORWARD + CLOSING.
    ///
    /// `INTERNAL` is deliberately not counted (`cxx:175`: `NOT VALID`).
    pub fn is_multiple(&self) -> bool {
        let mut lfound = Vec::new();
        let mut nkeyitem = self.item(FORWARD, &mut lfound);
        nkeyitem += self.item(CLOSING, &mut lfound);
        nkeyitem > 1
    }

    /// `IsFaulty()` (`cxx:185-191`): the FORWARD and REVERSED counts differ.
    pub fn is_faulty(&self) -> bool {
        let mut lfound = Vec::new();
        let nkey_rintem = self.item(FORWARD, &mut lfound);
        let nkey_fitem = self.item(REVERSED, &mut lfound);
        nkey_rintem != nkey_fitem
    }

    /// `IsInternal(Item)` (`cxx:193-243`).
    ///
    /// Keeps the INTERNAL items whose forward form is INTERNAL in the key, then
    /// appends the EXTERNAL items whose forward form is INTERNAL. Items that
    /// fail the orientation tests are dropped from the working copies.
    pub fn is_internal(&self, item: &mut Vec<TopoShape>) -> usize {
        item.clear();
        let Some(key) = self.the_key.as_ref() else {
            return 0;
        };

        let mut l_int: Vec<TopoShape> = self.the_items[INTERNAL as usize - 1].clone();
        let mut kept: Vec<TopoShape> = Vec::new();
        for item1 in &l_int {
            if item1.orientation() != crate::abs::Orientation::Internal {
                kept.push(item1.clone());
                continue;
            }
            let forward = item1.oriented(crate::abs::Orientation::Forward);
            let o_key1 = ori_in_sor(&forward, key, true);
            if o_key1 != INTERNAL {
                // `cxx:214`: `lINT.Remove(it1)`
                continue;
            }
            kept.push(item1.clone());
        }
        l_int = kept;

        let l_ext: Vec<TopoShape> = self.the_items[EXTERNAL as usize - 1].clone();
        for item2 in &l_ext {
            if item2.orientation() != crate::abs::Orientation::External {
                continue;
            }
            let forward = item2.oriented(crate::abs::Orientation::Forward);
            let o_key2 = ori_in_sor(&forward, key, true);
            if o_key2 == INTERNAL {
                l_int.push(item2.clone());
            }
        }

        item.extend(l_int);
        item.len()
    }
}
