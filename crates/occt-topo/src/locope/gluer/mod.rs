//! `LocOpe_Gluer` (`LocOpe_Gluer.hxx:33-82`, `.cxx:52-555`, `.lxx:18-91`):
//! glues the faces of a feature shape (`Snew`) onto the faces of a basis shape
//! (`Sbase`) and rebuilds the basis around them.
//!
//! Module table (`ModelingAlgorithms/TKFeat/LocOpe`):
//!
//! | OCCT | Rust |
//! |---|---|
//! | ctors / `Init` / `IsDone` / `ResultingShape` / `BasisShape` / `GluedShape` / `OpeType` / `Edges` / `TgtEdges` (`.lxx:18-91`) | this file |
//! | `Bind(Face, Face)` / `Bind(Edge, Edge)` / `DescendantFaces` (`cxx:52-152`, `cxx:335-351`) | this file |
//! | `Perform` (`cxx:154-333`) | [`perform`] |
//! | `GetOrientation` / `Contains` / `AddEdges` (file-local statics, `cxx:354-555`) | [`get_orientation`], [`statics`] |
//!
//! `Perform` is the package's assembly point: it feeds the bound face/edge pairs
//! into `LocOpe_WiresOnShape` (`cxx:171-194`), splits the basis with
//! `LocOpe_Spliter` (`cxx:203-207`), then rebuilds it with `LocOpe_Generator`
//! driven by a `LocOpe_GluedShape` (`cxx:227-231`). Once the generator is done
//! it re-encodes regularity across the glued boundaries (`cxx:260-298`) and
//! records the section edges that meet two target faces (`cxx:300-332`).
//!
//! Deviations, all for shapes the reference leaves null or raises on:
//! - `mySb` / `mySn` / `myRes` are `Option`s; a bind or `Perform` on them
//!   answers the `Standard_ConstructionError` of `cxx:72`, `cxx:82`, `cxx:95`,
//!   `cxx:149` and `cxx:165` as `Err`.
//! - `myMapEF` values are `Option<TopoShape>` because `cxx:135` nullifies an
//!   edge that borders two differently-bound faces; `Perform` skips those
//!   (`cxx:178-184`).
//! - `ResultingShape` (`cxx:44-51` in `.lxx`) and `DescendantFaces`
//!   (`cxx:339-342`) raise `StdFail_NotDone` before `Perform` succeeded; the
//!   ports answer `Err` / `None`.
//! - `AddEdges` (`cxx:471-555`) computes a mask of vertices near the new faces
//!   and then discards it: its only consumer is the empty `if (flag == 1) { }`
//!   at `cxx:549-551`, and it never touches `myRes`, `myEdges` or `myTgtEdges`.
//!   The port keeps the call site and documents the omission in [`statics`].

mod get_orientation;
mod perform;
mod statics;

use crate::abs::{Orientation, ShapeType};
use crate::locope::operation::LocOpeOperation;
use crate::locope::shape_maps::{ShapeDataMap, ShapeIndexedMap, ShapeValueMap};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{explorer_shapes, is_same};

/// `Standard_ConstructionError` raised by `LocOpe_Gluer::Bind` / `Perform`.
pub(super) const CONSTRUCTION_ERROR: &str = "LocOpe_Gluer: Standard_ConstructionError";

/// `LocOpe_Gluer` (`hxx:33-82`).
#[derive(Debug, Clone)]
pub struct Gluer {
    /// `bool myDone` (`hxx:69`).
    my_done: bool,
    /// `TopoDS_Shape mySb` (`hxx:70`): the basis shape.
    my_sb: Option<TopoShape>,
    /// `TopoDS_Shape mySn` (`hxx:71`): the glued (feature) shape.
    my_sn: Option<TopoShape>,
    /// `TopoDS_Shape myRes` (`hxx:72`).
    my_res: Option<TopoShape>,
    /// `TopAbs_Orientation myOri` (`hxx:73`): the decided orientation, starts
    /// at `TopAbs_INTERNAL` (`.lxx:26`).
    my_ori: Orientation,
    /// `LocOpe_Operation myOpe` (`hxx:74`).
    my_ope: LocOpeOperation,
    /// `myMapEF` (`hxx:76`): the glued faces keyed with their basis face, and
    /// the basis face of every edge of a glued face (`cxx:122-141`).
    my_map_ef: ShapeIndexedMap<Option<TopoShape>>,
    /// `myMapEE` (`hxx:77`): `Bind(Edge, Edge)` pairs (`cxx:145-152`).
    my_map_ee: ShapeDataMap,
    /// `myDescF` (`hxx:78`): the descendant list of every face of `mySb` and
    /// `mySn` (`cxx:209-225`, `cxx:239-258`).
    my_desc_f: ShapeValueMap<Vec<TopoShape>>,
    /// `myEdges` (`hxx:79`).
    my_edges: Vec<TopoShape>,
    /// `myTgtEdges` (`hxx:80`).
    my_tgt_edges: Vec<TopoShape>,
}

impl Default for Gluer {
    fn default() -> Self {
        Self::new()
    }
}

impl Gluer {
    /// `LocOpe_Gluer()` (`.lxx:20-26`).
    pub fn new() -> Self {
        Self {
            my_done: false,
            my_sb: None,
            my_sn: None,
            my_res: None,
            my_ori: Orientation::Internal,
            my_ope: LocOpeOperation::Invalid,
            my_map_ef: ShapeIndexedMap::new(),
            my_map_ee: ShapeDataMap::new(),
            my_desc_f: ShapeValueMap::new(),
            my_edges: Vec::new(),
            my_tgt_edges: Vec::new(),
        }
    }

    /// `LocOpe_Gluer(Sbase, Snew)` (`.lxx:30-37`).
    pub fn with_shapes(sbase: &TopoShape, snew: &TopoShape) -> Self {
        let mut this = Self::new();
        this.init(sbase, snew);
        this
    }

    /// `Init(Sbase, Snew)` (`cxx:52-62`).
    pub fn init(&mut self, sbase: &TopoShape, snew: &TopoShape) {
        self.my_sb = Some(sbase.clone());
        self.my_sn = Some(snew.clone());
        self.my_res = None;
        self.my_map_ef.clear();
        self.my_map_ee.clear();
        self.my_desc_f.clear();
        self.my_edges.clear();
        self.my_tgt_edges.clear();
        self.my_done = false;
        self.my_ori = Orientation::Internal;
        self.my_ope = LocOpeOperation::Invalid;
    }

    /// `IsDone()` (`.lxx:41-44`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `ResultingShape()` (`.lxx:48-54`); `Err` is `StdFail_NotDone`.
    pub fn resulting_shape(&self) -> Result<&TopoShape, &'static str> {
        if !self.my_done {
            return Err("StdFail_NotDone");
        }
        self.my_res.as_ref().ok_or("StdFail_NotDone")
    }

    /// `BasisShape()` (`.lxx:58-61`).
    pub fn basis_shape(&self) -> Option<&TopoShape> {
        self.my_sb.as_ref()
    }

    /// `GluedShape()` (`.lxx:65-68`).
    pub fn glued_shape(&self) -> Option<&TopoShape> {
        self.my_sn.as_ref()
    }

    /// `OpeType()` (`.lxx:73-76`).
    pub fn ope_type(&self) -> LocOpeOperation {
        self.my_ope
    }

    /// `Edges()` (`.lxx:80-83`).
    pub fn edges(&self) -> &[TopoShape] {
        &self.my_edges
    }

    /// `TgtEdges()` (`.lxx:87-90`).
    pub fn tgt_edges(&self) -> &[TopoShape] {
        &self.my_tgt_edges
    }

    /// `Bind(const TopoDS_Face& Fnew, const TopoDS_Face& Fbase)` (`cxx:66-142`).
    ///
    /// `Err` is the `Standard_ConstructionError` of `cxx:80` / `cxx:95`: both
    /// faces must already belong to the glued / basis shape.
    pub fn bind_face(&mut self, f_new: &Face, f_base: &Face) -> Result<(), String> {
        // `cxx:68-81`: find Fnew inside mySn.
        let f_new_or = self.oriented_member(f_new, true)?;
        // `cxx:87-98`: find Fbase inside mySb.
        let f_base_or = self.oriented_member(f_base, false)?;

        // `cxx:101`.
        let ori = get_orientation::get_orientation(&f_new_or, &f_base_or);

        // `cxx:103-121`.
        if self.my_ori == Orientation::Internal {
            self.my_ori = ori;
            if self.my_ori == Orientation::Reversed {
                if let Some(sn) = self.my_sn.as_mut() {
                    sn.reverse();
                }
                self.my_ope = LocOpeOperation::Cut;
            } else {
                self.my_ope = LocOpeOperation::Fuse;
            }
        } else if ori != Orientation::Forward {
            self.my_ope = LocOpeOperation::Invalid;
        }

        // `cxx:122-138`: every edge of the glued face records its basis face; an
        // edge shared by two glued faces bound to different basis faces is
        // nullified instead.
        for edg in explorer_shapes(&f_new_or.0, ShapeType::Edge) {
            if !self.my_map_ef.contains(&edg) {
                self.my_map_ef.add(&edg, Some(f_base_or.0.clone()));
                continue;
            }
            let same = self
                .my_map_ef
                .find_from_key(&edg)
                .and_then(|v| v.as_ref())
                .map(|v| is_same(v, &f_base_or.0))
                .unwrap_or(false);
            if !same {
                if let Some(slot) = self.my_map_ef.change_from_key(&edg) {
                    *slot = None;
                }
            }
        }
        // `cxx:141`.
        self.my_map_ef.add(&f_new_or.0, Some(f_base_or.0.clone()));
        Ok(())
    }

    /// `Bind(const TopoDS_Edge& Enew, const TopoDS_Edge& Ebase)` (`cxx:144-152`).
    ///
    /// `Err` is the `Standard_ConstructionError` of `cxx:149`: one edge of the
    /// glued shape cannot be bound to two different basis edges.
    pub fn bind_edge(&mut self, e_new: &Edge, e_base: &Edge) -> Result<(), String> {
        if self.my_map_ee.is_bound(&e_new.0) {
            let same = self
                .my_map_ee
                .find(&e_new.0)
                .map(|v| is_same(v, &e_base.0))
                .unwrap_or(false);
            if !same {
                return Err(CONSTRUCTION_ERROR.to_string());
            }
        }
        self.my_map_ee.bind(&e_new.0, &e_base.0);
        Ok(())
    }

    /// `DescendantFaces(const TopoDS_Face& F)` (`cxx:335-351`).
    ///
    /// `Err` is the `StdFail_NotDone` raise of `cxx:339-342`; an unbound face
    /// answers the reference's `static nullList` as `Ok(&[])`.
    pub fn descendant_faces(&self, f: &Face) -> Result<&[TopoShape], &'static str> {
        if !self.my_done {
            return Err("StdFail_NotDone");
        }
        Ok(self.my_desc_f.find(&f.0).map(|v| v.as_slice()).unwrap_or(&[]))
    }

    /// `TopExp_Explorer(mySn | mySb, TopAbs_FACE)` up to `F`, answered with the
    /// orientation the shape gives it (`cxx:68-101`, `cxx:87-100`).
    fn oriented_member(&self, f: &Face, new_side: bool) -> Result<Face, String> {
        let source = if new_side {
            self.my_sn.as_ref()
        } else {
            self.my_sb.as_ref()
        };
        let found = match source {
            Some(shape) => explorer_shapes(shape, ShapeType::Face)
                .into_iter()
                .find(|cur| is_same(cur, &f.0)),
            None => None,
        };
        let Some(cur) = found else {
            return Err(CONSTRUCTION_ERROR.to_string());
        };
        Ok(Face(f.0.oriented(cur.orientation())))
    }
}
