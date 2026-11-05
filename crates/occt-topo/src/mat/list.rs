//! `MAT_TList` (`MAT_TList.gxx:19-345`), the doubly linked list that
//! `MAT_ListOfEdge` (`MAT_ListOfEdge_0.cxx`) and `MAT_ListOfBisector`
//! (`MAT_ListOfBisector_0.cxx`) instantiate.
//!
//! OCCT's list is a chain of `MAT_TListNode`s with a cursor
//! (`thecurrentnode` / `thecurrentindex`). `Loop()` (`gxx:317-321`) closes the
//! chain into a ring, and every cursor move then wraps; a list that was never
//! looped walks off the end into a null node, so `More()` answers false.
//!
//! The port keeps the same state machine over a `Vec`: `cursor` is the 0-based
//! index of `thecurrentnode` (`None` is the null node), and `tombstone` holds
//! the item of the node the cursor still points at after `Unlink`
//! (`gxx:153-179`), where OCCT leaves the cursor on the detached node.
//! `Loop()` only flips `looped`, because wrap-around is a cursor rule here and
//! the storage order never changes.

/// `MAT_TList` (`MAT_TList.gxx:19-345`) over the items it stores.
#[derive(Debug, Clone, PartialEq)]
pub struct MatTList<T> {
    /// The live nodes, in list order.
    items: Vec<T>,
    /// `thecurrentindex` (`gxx:21`): the 1-based position of the current node,
    /// `0` while the cursor is null. `Unlink` decrements it and `LinkBefore` /
    /// `FrontAdd` / `Permute` increment it, exactly as OCCT does.
    current_index: usize,
    /// `thecurrentnode`, as a 0-based index; `None` is the null node.
    cursor: Option<usize>,
    /// The item of the node `Unlink` detached, which the cursor still points
    /// at until the next `Next()`.
    tombstone: Option<T>,
    /// `Loop()` was called (`gxx:317-321`): `Next` / `Previous` wrap.
    looped: bool,
}

impl<T> Default for MatTList<T> {
    /// `MAT_TList()` (`gxx:19-23`): null first / last node, index and count 0.
    fn default() -> Self {
        Self::new()
    }
}

impl<T> MatTList<T> {
    /// `MAT_TList()` (`gxx:19-23`).
    pub fn new() -> Self {
        Self {
            items: Vec::new(),
            current_index: 0,
            cursor: None,
            tombstone: None,
            looped: false,
        }
    }

    /// `Index()` (`MAT_TList.hxx`): the 1-based position of the current node.
    pub fn index(&self) -> usize {
        self.current_index
    }

    /// `First()` (`gxx:27-29`): the cursor moves to the first node.
    pub fn first(&mut self) {
        self.tombstone = None;
        self.cursor = if self.items.is_empty() { None } else { Some(0) };
        self.current_index = 1;
    }

    /// `Last()` (`gxx:33-35`).
    pub fn last(&mut self) {
        self.tombstone = None;
        self.cursor = self.items.len().checked_sub(1);
        self.current_index = self.items.len();
    }

    /// `Init(anitem)` (`gxx:41-53`): `First()` then walk until the node whose
    /// item is `anitem` is current. OCCT compares the handles themselves, so
    /// the comparison is the caller's `PartialEq` (identity for slot ids).
    pub fn init(&mut self, an_item: &T)
    where
        T: PartialEq,
    {
        self.first();
        while self.more() {
            if self.current() == an_item {
                break;
            }
            self.next();
        }
    }

    /// `Next()` (`gxx:57-63`): off the end of a linear list the cursor becomes
    /// null; a looped list wraps (`thecurrentindex % thenumberofitems + 1`).
    pub fn next(&mut self) {
        let len = self.items.len();
        if self.tombstone.is_some() {
            // The cursor already sits on the node that followed the detached
            // one (`gxx:65`: `thecurrentnode = thecurrentnode->Next()`); only
            // `thecurrentindex` still has to be advanced. `Next()` does nothing
            // once the list is empty (`gxx:59`), so the cursor stays on the
            // detached node and `Current()` keeps answering its item.
            if len == 0 {
                return;
            }
            self.tombstone = None;
            if self.cursor.is_some() {
                self.current_index = self.current_index % len + 1;
            }
            return;
        }
        let Some(index) = self.cursor else {
            return;
        };
        if len == 0 {
            self.cursor = None;
        } else if self.looped {
            self.cursor = Some((index + 1) % len);
        } else if index + 1 < len {
            self.cursor = Some(index + 1);
        } else {
            self.cursor = None;
        }
        if self.cursor.is_some() {
            self.current_index = self.current_index % len + 1;
        }
    }

    /// `Previous()` (`gxx:67-71`).
    pub fn previous(&mut self) {
        let Some(index) = self.cursor else {
            return;
        };
        let len = self.items.len();
        if len == 0 {
            self.cursor = None;
        } else if self.looped {
            self.cursor = Some((index + len - 1) % len);
        } else if index > 0 {
            self.cursor = Some(index - 1);
        } else {
            self.cursor = None;
        }
        if self.cursor.is_some() {
            self.current_index = (self.current_index + len - 2) % len + 1;
        }
    }

    /// `More()` (`gxx:79-81`): `!thecurrentnode.IsNull()`.
    pub fn more(&self) -> bool {
        self.cursor.is_some() || self.tombstone.is_some()
    }

    /// `Current()` (`gxx:86-89`).
    pub fn current(&self) -> &T {
        match &self.tombstone {
            Some(item) => item,
            None => &self.items[self.cursor.expect("MAT_TList::Current: null node")],
        }
    }

    /// `Current(anitem)` (`gxx:93-96`).
    pub fn set_current(&mut self, an_item: T) {
        match &mut self.tombstone {
            Some(item) => *item = an_item,
            None => {
                let index = self.cursor.expect("MAT_TList::Current: null node");
                self.items[index] = an_item;
            }
        }
    }

    /// `FirstItem()` (`gxx:99-101`).
    pub fn first_item(&self) -> &T {
        &self.items[0]
    }

    /// `LastItem()` (`gxx:105-107`).
    pub fn last_item(&self) -> &T {
        &self.items[self.items.len() - 1]
    }

    /// `PreviousItem()` (`gxx:113-117`): `thecurrentnode->Previous()`; a looped
    /// list always has one, a linear one is undefined at the first node.
    pub fn previous_item(&self) -> &T {
        if let Some(item) = &self.tombstone {
            if self.items.is_empty() {
                // The detached node was the last of the ring, so its own
                // `Previous()` is the node itself (`gxx:118`).
                return item;
            }
        }
        let index = self
            .cursor
            .expect("MAT_TList::PreviousItem: null node");
        let len = self.items.len();
        if self.looped {
            &self.items[(index + len - 1) % len]
        } else {
            &self.items[index - 1]
        }
    }

    /// `NextItem()` (`gxx:121-125`).
    pub fn next_item(&self) -> &T {
        let len = self.items.len();
        if self.tombstone.is_some() {
            return &self.items[self.cursor.expect("MAT_TList::NextItem: null node")];
        }
        let index = self.cursor.expect("MAT_TList::NextItem: null node");
        if self.looped {
            &self.items[(index + 1) % len]
        } else {
            &self.items[index + 1]
        }
    }

    /// `Brackets(anindex)` (`gxx:129-149`): walks the cursor to the 1-based
    /// `anindex` and returns that item.
    pub fn brackets(&mut self, an_index: usize) -> &T {
        assert!(an_index >= 1 && an_index <= self.items.len());
        self.tombstone = None;
        self.cursor = Some(an_index - 1);
        self.current_index = an_index;
        &self.items[an_index - 1]
    }

    /// `Unlink()` (`gxx:153-179`): the current node leaves the list; the cursor
    /// stays on it so that the next `Next()` steps to its successor.
    pub fn unlink(&mut self) {
        let index = self.cursor.expect("MAT_TList::Unlink: null node");
        assert!(self.tombstone.is_none());
        let item = self.items.remove(index);
        self.tombstone = Some(item);
        self.cursor = if index < self.items.len() {
            Some(index)
        } else if self.looped && !self.items.is_empty() {
            Some(0)
        } else {
            None
        };
        // `gxx:176`: `thecurrentindex--`.
        self.current_index = self.current_index.saturating_sub(1);
    }

    /// `LinkBefore(anitem)` (`gxx:183-207`): the new node is inserted before the
    /// current one, which stays current.
    pub fn link_before(&mut self, an_item: T) {
        let index = self.cursor.expect("MAT_TList::LinkBefore: null node");
        self.items.insert(index, an_item);
        self.cursor = Some(index + 1);
        if self.current_index > 0 {
            self.current_index += 1;
        }
    }

    /// `LinkAfter(anitem)` (`gxx:211-233`): the new node follows the current one,
    /// which stays current, so the next `Next()` lands on the new node.
    pub fn link_after(&mut self, an_item: T) {
        let index = self.cursor.expect("MAT_TList::LinkAfter: null node");
        self.items.insert(index + 1, an_item);
    }

    /// `FrontAdd(anitem)` (`gxx:237-254`).
    pub fn front_add(&mut self, an_item: T) {
        self.items.insert(0, an_item);
        if self.cursor.is_some() {
            self.cursor = self.cursor.map(|index| index + 1);
        }
        if self.current_index > 0 {
            self.current_index += 1;
        }
    }

    /// `BackAdd(anitem)` (`gxx:258-274`).
    pub fn back_add(&mut self, an_item: T) {
        self.items.push(an_item);
    }

    /// `Permute()` (`gxx:278-313`): exchanges the current node with its
    /// successor.
    pub fn permute(&mut self) {
        let index = self.cursor.expect("MAT_TList::Permute: null node");
        let len = self.items.len();
        let next = if index + 1 < len { index + 1 } else { 0 };
        self.items.swap(index, next);
        self.cursor = Some(next);
        // `gxx:309`: `thecurrentindex++`.
        self.current_index += 1;
    }

    /// `Loop()` (`gxx:317-321`): closes the chain into a ring. A `Vec` cursor
    /// wraps by rule, so only the flag has to change.
    pub fn loop_list(&mut self) {
        self.looped = true;
    }

    /// `Number()` (`MAT_TList.hxx`).
    pub fn number(&self) -> usize {
        self.items.len()
    }

    /// `IsEmpty()` (`MAT_TList.hxx`).
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// The items in list order, for callers that only read the list.
    pub fn items(&self) -> &[T] {
        &self.items
    }
}
