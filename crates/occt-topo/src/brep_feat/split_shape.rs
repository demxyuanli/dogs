//! `BRepFeat_SplitShape` — splitting a shape with wires or edges of a face.
//!
//! Source: `ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_SplitShape.hxx:1-121`,
//! `BRepFeat_SplitShape.cxx:26-90`, `BRepFeat_SplitShape.lxx:22-82`.
//!
//! The class is the `BRepBuilderAPI_MakeShape` facade of `LocOpe_Spliter`
//! driven by a `LocOpe_WiresOnShape`; every `Add` overload forwards to one of
//! the `WiresOnShape::Bind` methods (`lxx:59-82`).

use crate::abs::ShapeType;
use crate::locope::{Spliter, WiresOnShape};
use crate::shape::{Compound, Edge, Face, TopoShape, Wire};
use crate::topo_tools;

/// `BRepFeat_SplitShape` (`hxx:50-121`).
#[derive(Debug, Clone)]
pub struct BRepFeatSplitShape {
    my_s_shape: Spliter,
    /// `Handle(LocOpe_WiresOnShape) myWOnShape` — a null handle until `Init`
    /// or the `S` constructor (`lxx:19`, `:28-31`, `:41-51`).
    my_w_on_shape: Option<WiresOnShape>,
    my_right: Vec<TopoShape>,
    my_shape: Option<TopoShape>,
    my_done: bool,
}

impl Default for BRepFeatSplitShape {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepFeatSplitShape {
    /// `BRepFeat_SplitShape()` (`lxx:19`): both members default, the
    /// `WiresOnShape` handle is null.
    pub fn new() -> Self {
        Self {
            my_s_shape: Spliter::new(),
            my_w_on_shape: None,
            my_right: Vec::new(),
            my_shape: None,
            my_done: false,
        }
    }

    /// `BRepFeat_SplitShape(S)` (`lxx:23-28`).
    pub fn with_shape(s: &TopoShape) -> Self {
        Self {
            my_s_shape: Spliter::with_shape(s),
            my_w_on_shape: Some(WiresOnShape::new(s.clone())),
            my_right: Vec::new(),
            my_shape: None,
            my_done: false,
        }
    }

    /// `Add(const NCollection_Sequence<TopoDS_Shape>&)` (`lxx:33-36`): the
    /// edges / wires are put on the faces of the initial shape they lie on.
    ///
    /// A null `myWOnShape` dereferences in the reference; the port answers the
    /// `Standard_NullObject` failure instead.
    pub fn add_edges(&mut self, the_edges: &[TopoShape]) -> Result<bool, String> {
        match self.my_w_on_shape.as_mut() {
            Some(w) => w.add(the_edges),
            None => Err("Standard_NullObject: WiresOnShape".to_string()),
        }
    }

    /// `Init(S)` (`lxx:40-51`): re-initializes the spliter and the
    /// `WiresOnShape` (creating the handle when it is null).
    pub fn init(&mut self, s: &TopoShape) {
        self.my_s_shape.init(s);
        match self.my_w_on_shape.as_mut() {
            Some(w) => w.init(s.clone()),
            None => self.my_w_on_shape = Some(WiresOnShape::new(s.clone())),
        }
        self.my_right.clear();
        self.my_shape = None;
        self.my_done = false;
    }

    /// `SetCheckInterior(ToCheckInterior)` (`lxx:55-58`). A null handle is a
    /// dereference in the reference; the port keeps it a no-op.
    pub fn set_check_interior(&mut self, to_check_interior: bool) {
        if let Some(w) = self.my_w_on_shape.as_mut() {
            w.set_check_interior(to_check_interior);
        }
    }

    /// `Add(const TopoDS_Wire& W, const TopoDS_Face& F)` (`lxx:62-65`).
    pub fn add_wire(&mut self, w: &Wire, f: &Face) -> Result<(), String> {
        match self.my_w_on_shape.as_mut() {
            Some(ww) => ww.bind_wire(w, f),
            None => Err("Standard_NullObject: WiresOnShape".to_string()),
        }
    }

    /// `Add(const TopoDS_Edge& E, const TopoDS_Face& F)` (`lxx:69-72`).
    pub fn add_edge(&mut self, e: &Edge, f: &Face) -> Result<(), String> {
        match self.my_w_on_shape.as_mut() {
            Some(ww) => ww.bind_edge(e, f),
            None => Err("Standard_NullObject: WiresOnShape".to_string()),
        }
    }

    /// `Add(const TopoDS_Compound& Comp, const TopoDS_Face& F)` (`lxx:76-79`).
    pub fn add_compound(&mut self, comp: &Compound, f: &Face) -> Result<(), String> {
        match self.my_w_on_shape.as_mut() {
            Some(ww) => ww.bind_compound(comp, f),
            None => Err("Standard_NullObject: WiresOnShape".to_string()),
        }
    }

    /// `Add(const TopoDS_Edge& E, const TopoDS_Edge& EOn)` (`lxx:82-85`).
    pub fn add_edge_on_edge(&mut self, e: &Edge, e_on: &Edge) -> Result<(), String> {
        match self.my_w_on_shape.as_mut() {
            Some(ww) => {
                ww.bind_edge_to_edge(e, e_on);
                Ok(())
            }
            None => Err("Standard_NullObject: WiresOnShape".to_string()),
        }
    }

    /// `DirectLeft()` (`cxx:37-40`).
    pub fn direct_left(&self) -> Option<&[TopoShape]> {
        self.my_s_shape.direct_left()
    }

    /// `Left()` (`cxx:45-48`). The reference raises `StdFail_NotDone` when the
    /// spliter is not done; the port answers `None`.
    pub fn left(&self) -> Option<&[TopoShape]> {
        self.my_s_shape.left()
    }

    /// `Right()` (`cxx:51-73`): the faces of the result that are not in
    /// `Left()`, computed once and cached in `myRight`.
    pub fn right(&mut self) -> &[TopoShape] {
        if self.my_right.is_empty() {
            let mut map_of_left = crate::locope::ShapeSet::new();
            if let Some(left) = self.my_s_shape.left() {
                for s in left {
                    map_of_left.add(s);
                }
            }
            if let Some(shape) = self.my_shape.as_ref() {
                for face in topo_tools::explorer_shapes(shape, ShapeType::Face) {
                    if !map_of_left.contains(&face) {
                        self.my_right.push(face);
                    }
                }
            }
        }
        &self.my_right
    }

    /// `Build()` (`cxx:26-35`): `mySShape.Perform(myWOnShape)`, then `Done()`
    /// + `myShape` and `myRight.Clear()` when the spliter is done.
    pub fn build(&mut self) -> Result<(), String> {
        let mut w = self.my_w_on_shape.take();
        let res = match w.as_mut() {
            Some(ws) => self.my_s_shape.perform(ws),
            None => Err("Standard_NullObject: WiresOnShape".to_string()),
        };
        self.my_w_on_shape = w;
        res?;
        if self.my_s_shape.is_done() {
            self.my_done = true;
            self.my_shape = self.my_s_shape.resulting_shape().cloned();
            self.my_right.clear();
        }
        Ok(())
    }

    /// The inherited `IsDone()` of `BRepBuilderAPI_MakeShape`.
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// The inherited `Shape()`.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_shape.as_ref()
    }

    /// `IsDeleted(F)` (`cxx:76-82`): the descendant list is empty.
    pub fn is_deleted(&self, f: &TopoShape) -> bool {
        self.my_s_shape.descendant_shapes(f).is_empty()
    }

    /// `Modified(F)` (`cxx:85-88`) — `mySShape.DescendantShapes(F)`.
    pub fn modified(&self, f: &TopoShape) -> Vec<TopoShape> {
        self.my_s_shape.descendant_shapes(f)
    }
}
