//! `BRepAlgo_Image` — the link table between a shape and the shapes it became.
//!
//! Source: `ModelingAlgorithms/TKBool/BRepAlgo/BRepAlgo_Image.cxx:1-353` and
//! `.hxx:36-101` (OCCT 8.0.0).
//!
//! | OCCT | lines | port |
//! |---|---|---|
//! | `BRepAlgo_Image()` | `cxx:27` | [`BRepAlgoImage::new`] |
//! | `SetRoot` | `cxx:31-34` | [`BRepAlgoImage::set_root`] |
//! | `Bind(OldS, NewS)` | `cxx:38-50` | [`BRepAlgoImage::bind`] |
//! | `Bind(OldS, L)` | `cxx:54-73` | [`BRepAlgoImage::bind_list`] |
//! | `Add(OldS, NewS)` | `cxx:87-94` | [`BRepAlgoImage::add`] |
//! | `Add(OldS, L)` | `cxx:98-106` | [`BRepAlgoImage::add_list`] |
//! | `Clear` | `cxx:76-83` | [`BRepAlgoImage::clear`] |
//! | `Remove` | `cxx:109-132` | [`BRepAlgoImage::remove`] |
//! | `RemoveRoot` | `cxx:299-330` | [`BRepAlgoImage::remove_root`] |
//! | `ReplaceRoot` | `cxx:334-353` | [`BRepAlgoImage::replace_root`] |
//! | `Roots` | `cxx:137-141` | [`BRepAlgoImage::roots`] |
//! | `IsImage` | `cxx:145-149` | [`BRepAlgoImage::is_image`] |
//! | `ImageFrom` | `cxx:153-160` | [`BRepAlgoImage::image_from`] |
//! | `Root` | `cxx:164-189` | [`BRepAlgoImage::root`] |
//! | `HasImage` | `cxx:193-196` | [`BRepAlgoImage::has_image`] |
//! | `Image` | `cxx:204-211` | [`BRepAlgoImage::image`] |
//! | `LastImage` | `cxx:218-238` | [`BRepAlgoImage::last_image`] |
//! | `Compact` | `cxx:242-266` | [`BRepAlgoImage::compact`] |
//! | `Filter` | `cxx:270-295` | [`BRepAlgoImage::filter`] |
//!
//! Translation notes:
//!
//! * `down` is `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`
//!   over `TopTools_ShapeMapHasher`, i.e. `HashMap<ShapeMapKey, Vec<TopoShape>>`
//!   keyed by [`shape_map_key`] (`IsSame` identity) with `Vec` in insertion
//!   order.
//! * `up` is the same map without a list. `NCollection_DataMap` has no
//!   iteration order in OCCT; [`BRepAlgoImage::filter`] is the only walk over
//!   it and it removes one violating key per pass until none is left, so the
//!   port keeps `up_order` (insertion order) instead of a hash order and
//!   reaches the same fixed point.
//! * `Image(S)` answers a `static NCollection_List` in the reference
//!   (`cxx:208-210`), which accumulates across calls — a reference defect.
//!   The port answers a fresh one-element `Vec`.
//! * `Bind` / `Add` / `Remove` / `ImageFrom` / `Root` raise
//!   `Standard_ConstructionError`; the port answers `Err` with the same text.

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key, ShapeMapKey};

/// `class BRepAlgo_Image` (`hxx:35-101`).
#[derive(Clone, Default)]
pub struct BRepAlgoImage {
    /// `roots` (`hxx:95`).
    roots: Vec<TopoShape>,
    /// `up` (`hxx:96`): image -> its generator.
    up: HashMap<ShapeMapKey, TopoShape>,
    /// Insertion order of `up` keys (OCCT iterates it unordered).
    up_order: Vec<ShapeMapKey>,
    /// `down` (`hxx:97`): generator -> its images.
    down: HashMap<ShapeMapKey, Vec<TopoShape>>,
}

impl BRepAlgoImage {
    pub fn new() -> Self {
        Self::default()
    }

    /// `SetRoot` (`cxx:31-34`).
    pub fn set_root(&mut self, s: &TopoShape) {
        self.roots.push(s.clone());
    }

    /// `Bind(OldS, NewS)` (`cxx:38-50`).
    pub fn bind(&mut self, old_s: &TopoShape, new_s: &TopoShape) -> Result<(), String> {
        let key = shape_map_key(old_s);
        if self.down.contains_key(&key) {
            return Err(" BRepAlgo_Image::Bind".into());
        }
        self.down.insert(key, vec![new_s.clone()]);
        self.set_up(new_s, old_s);
        Ok(())
    }

    /// `Bind(OldS, L)` (`cxx:54-73`).
    pub fn bind_list(&mut self, old_s: &TopoShape, list: &[TopoShape]) -> Result<(), String> {
        if self.has_image(old_s) {
            return Err(" BRepAlgo_Image::Bind".into());
        }
        for v in list {
            if !self.has_image(old_s) {
                self.bind(old_s, v)?;
            } else {
                self.add(old_s, v)?;
            }
        }
        Ok(())
    }

    /// `Add(OldS, NewS)` (`cxx:87-94`).
    pub fn add(&mut self, old_s: &TopoShape, new_s: &TopoShape) -> Result<(), String> {
        if !self.has_image(old_s) {
            return Err(" BRepAlgo_Image::Add".into());
        }
        self.down
            .entry(shape_map_key(old_s))
            .or_default()
            .push(new_s.clone());
        self.set_up(new_s, old_s);
        Ok(())
    }

    /// `Add(OldS, L)` (`cxx:98-106`).
    pub fn add_list(&mut self, old_s: &TopoShape, list: &[TopoShape]) -> Result<(), String> {
        for v in list {
            self.add(old_s, v)?;
        }
        Ok(())
    }

    /// `Clear` (`cxx:76-83`).
    pub fn clear(&mut self) {
        self.roots.clear();
        self.up.clear();
        self.up_order.clear();
        self.down.clear();
    }

    /// `Remove(S)` (`cxx:109-132`).
    pub fn remove(&mut self, s: &TopoShape) -> Result<(), String> {
        let key = shape_map_key(s);
        let Some(old_s) = self.up.get(&key).cloned() else {
            return Err(" BRepAlgo_Image::Remove".into());
        };
        let old_key = shape_map_key(&old_s);
        let mut emptied = false;
        if let Some(list) = self.down.get_mut(&old_key) {
            if let Some(i) = list.iter().position(|x| is_same(x, s)) {
                list.remove(i);
            }
            emptied = list.is_empty();
        }
        if emptied {
            self.down.remove(&old_key);
        }
        self.up.remove(&key);
        self.up_order.retain(|k| *k != key);
        Ok(())
    }

    /// `RemoveRoot(Root)` (`cxx:299-330`).
    pub fn remove_root(&mut self, root: &TopoShape) {
        let Some(pos) = self.roots.iter().position(|x| is_same(x, root)) else {
            return;
        };
        self.roots.remove(pos);

        let root_key = shape_map_key(root);
        if let Some(new_s) = self.down.get(&root_key).cloned() {
            for img in &new_s {
                let img_key = shape_map_key(img);
                let owned = self.up.get(&img_key).map(|o| is_same(o, root)).unwrap_or(false);
                if owned {
                    self.up.remove(&img_key);
                    self.up_order.retain(|k| *k != img_key);
                }
            }
            self.down.remove(&root_key);
        }
    }

    /// `ReplaceRoot(OldRoot, NewRoot)` (`cxx:334-353`).
    pub fn replace_root(
        &mut self,
        old_root: &TopoShape,
        new_root: &TopoShape,
    ) -> Result<(), String> {
        if !self.has_image(old_root) {
            return Ok(());
        }
        let list = self.image(old_root);
        if self.has_image(new_root) {
            self.add_list(new_root, &list)?;
        } else {
            self.bind_list(new_root, &list)?;
        }
        self.set_root(new_root);
        self.remove_root(old_root);
        Ok(())
    }

    /// `Roots` (`cxx:137-141`).
    pub fn roots(&self) -> &[TopoShape] {
        &self.roots
    }

    /// `IsImage(S)` (`cxx:145-149`).
    pub fn is_image(&self, s: &TopoShape) -> bool {
        self.up.contains_key(&shape_map_key(s))
    }

    /// `ImageFrom(S)` (`cxx:153-160`).
    pub fn image_from(&self, s: &TopoShape) -> Result<TopoShape, String> {
        match self.up.get(&shape_map_key(s)) {
            Some(o) => Ok(o.clone()),
            None => Err(" BRepAlgo_Image::ImageFrom".into()),
        }
    }

    /// `Root(S)` (`cxx:164-189`).
    pub fn root(&self, s: &TopoShape) -> Result<TopoShape, String> {
        let first = match self.up.get(&shape_map_key(s)) {
            Some(o) => o.clone(),
            None => return Err(" BRepAlgo_Image::FirstImageFrom".into()),
        };
        let mut s1 = first;
        let mut s2 = s.clone();
        if is_same(&s1, &s2) {
            return Ok(self.up.get(&shape_map_key(s)).cloned().unwrap_or_else(|| s.clone()));
        }
        while let Some(next) = self.up.get(&shape_map_key(&s1)).cloned() {
            s2 = s1;
            s1 = next;
            if is_same(&s1, &s2) {
                break;
            }
        }
        Ok(self.up.get(&shape_map_key(&s2)).cloned().unwrap_or(s2))
    }

    /// `HasImage(S)` (`cxx:193-196`).
    pub fn has_image(&self, s: &TopoShape) -> bool {
        self.down.contains_key(&shape_map_key(s))
    }

    /// `Image(S)` (`cxx:204-211`).
    pub fn image(&self, s: &TopoShape) -> Vec<TopoShape> {
        match self.down.get(&shape_map_key(s)) {
            Some(l) => l.clone(),
            None => vec![s.clone()],
        }
    }

    /// `LastImage(S, L)` (`cxx:218-238`): append to `l` the deepest images of
    /// `s`, stopping at a self-image.
    pub fn last_image(&self, s: &TopoShape, l: &mut Vec<TopoShape>) {
        match self.down.get(&shape_map_key(s)) {
            None => l.push(s.clone()),
            Some(list) => {
                for v in list.clone() {
                    if is_same(&v, s) {
                        l.push(s.clone());
                    } else {
                        self.last_image(&v, l);
                    }
                }
            }
        }
    }

    /// `Compact` (`cxx:242-266`): keep only the root -> last-image links.
    pub fn compact(&mut self) -> Result<(), String> {
        let mut m: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();
        for s in self.roots.clone() {
            let mut li = Vec::new();
            if self.has_image(&s) {
                self.last_image(&s, &mut li);
            }
            m.insert(shape_map_key(&s), li);
        }
        self.up.clear();
        self.up_order.clear();
        self.down.clear();
        for s in self.roots.clone() {
            if let Some(li) = m.get(&shape_map_key(&s)).cloned() {
                self.bind_list(&s, &li)?;
            }
        }
        Ok(())
    }

    /// `Filter(S, T)` (`cxx:270-295`): drop every image of type `T` that is not
    /// under `S`, one per pass, until none is left.
    pub fn filter(&mut self, s: &TopoShape, t: ShapeType) -> Result<(), String> {
        let keep: Vec<TopoShape> = explorer_shapes(s, t);
        let mut changed = true;
        while changed {
            changed = false;
            for key in self.up_order.clone() {
                let Some(a_s) = self.up.get(&key).cloned() else {
                    continue;
                };
                if a_s.shape_type() == t && !keep.iter().any(|k| is_same(k, &a_s)) {
                    self.remove(&a_s)?;
                    changed = true;
                    break;
                }
            }
        }
        Ok(())
    }

    fn set_up(&mut self, new_s: &TopoShape, old_s: &TopoShape) {
        let key = shape_map_key(new_s);
        if self.up.insert(key.clone(), old_s.clone()).is_none() {
            self.up_order.push(key);
        }
    }
}
