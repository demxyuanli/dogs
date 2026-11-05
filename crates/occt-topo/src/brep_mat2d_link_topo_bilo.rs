//! `BRepMAT2d_LinkTopoBilo` (`BRepMAT2d_LinkTopoBilo.cxx:39-322`,
//! `.hxx:20-90`): the link between the topology of the spine (edges and
//! vertices of a face) and the `MAT_BasicElt`s of its bisecting locus.
//!
//! `Perform` walks the wires of the explored face and, for each contour,
//! `LinkToContour` decides which curve of the contour each basic element comes
//! from (`KC`, positive when the element follows the contour direction); the
//! basic element is then bound to the edge `TopoSeq[|KC|]` or, when the element
//! is a `Geom2d_CartesianPoint`, to one of that edge's vertices.
//!
//! `myMap` is an `NCollection_DataMap<TopoDS_Shape, List<handle(MAT_BasicElt)>>`
//! with the default `TopTools_ShapeMapHasher`, i.e. keyed by `IsSame`
//! (`TShape` + location, orientation ignored): [`shape_map_key`].
//!
//! Deviation: OCCT iterates `NCollection_DataMap<int, int> LinkBECont` in hash
//! bucket order (`cxx:169-170`); the port keeps the same bindings in a
//! `BTreeMap` and iterates keys ascending. Each basic element is bound to one
//! shape at most, so the order only matters if two elements of the same contour
//! share a shape, which the contour walk below cannot produce.

use std::collections::{BTreeMap, HashMap};

use crate::abs::{Orientation, ShapeType};
use crate::brep_mat2d_bisecting_locus::BRepMat2dBisectingLocus;
use crate::brep_mat2d_explorer::BRepMat2dExplorer;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::mat::BasicEltId;
use crate::mat2d::Mat2dItem;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools::{edge_vertices, shape_map_key, wires_of_face, ShapeMapKey};

/// `class BRepMAT2d_LinkTopoBilo` (`BRepMAT2d_LinkTopoBilo.hxx:30-90`).
#[derive(Clone, Default)]
pub struct BRepMat2dLinkTopoBilo {
    /// `myMap` (`.hxx`): `NCollection_DataMap<TopoDS_Shape,
    /// NCollection_List<handle(MAT_BasicElt)>>`.
    my_map: HashMap<ShapeMapKey, Vec<BasicEltId>>,
    /// `myBEShape` (`.hxx`): `NCollection_DataMap<handle(MAT_BasicElt),
    /// TopoDS_Shape>`.
    my_be_shape: HashMap<BasicEltId, TopoShape>,
    /// `myKey` (`.hxx`).
    my_key: Option<ShapeMapKey>,
    /// `current` (`.hxx`): 1-based cursor inside `myMap(myKey)`.
    current: usize,
    /// `isEmpty` (`.hxx`).
    is_empty: bool,
}

impl BRepMat2dLinkTopoBilo {
    /// `BRepMAT2d_LinkTopoBilo()` (`cxx:42-46`).
    pub fn new() -> Self {
        Self {
            my_map: HashMap::new(),
            my_be_shape: HashMap::new(),
            my_key: None,
            current: 0,
            is_empty: true,
        }
    }

    /// `BRepMAT2d_LinkTopoBilo(Explo, BiLo)` (`cxx:49-53`).
    pub fn with_locus(explo: &BRepMat2dExplorer, bilo: &BRepMat2dBisectingLocus) -> Self {
        let mut link = Self::new();
        link.perform(explo, bilo);
        link
    }

    /// `Perform(Explo, BiLo)` (`cxx:58-81`). A spine that is not a face raises
    /// `Standard_ConstructionError` in OCCT (`cxx:79`); the port panics.
    pub fn perform(&mut self, explo: &BRepMat2dExplorer, bilo: &BRepMat2dBisectingLocus) {
        self.my_map.clear();
        self.my_be_shape.clear();

        let shape = explo
            .shape()
            .expect("BRepMAT2d_LinkTopoBilo::Perform: null shape")
            .clone();
        if shape.shape_type() != ShapeType::Face {
            panic!("BRepMAT2d_LinkTopoBilo::Perform: Standard_ConstructionError");
        }

        // `cxx:68-77`: `TopExp_Explorer Exp(S, TopAbs_WIRE)`, `IndContour`
        // starts at 1.
        for (i, w) in wires_of_face(&Face(shape)).iter().enumerate() {
            self.link_to_wire(w, explo, i as i32 + 1, bilo);
        }
    }

    /// `Init(S)` (`cxx:85-97`).
    pub fn init(&mut self, s: &TopoShape) {
        self.is_empty = false;
        self.current = 1;
        let key = shape_map_key(s);
        if self.my_map.contains_key(&key) {
            self.my_key = Some(key);
        } else {
            self.is_empty = true;
            self.my_key = None;
        }
    }

    /// `More()` (`cxx:102-107`).
    pub fn more(&self) -> bool {
        if self.is_empty {
            return false;
        }
        match self.my_key.as_ref().and_then(|k| self.my_map.get(k)) {
            // `cxx:107`: `current <= myMap(myKey).Length()`.
            Some(list) => self.current <= list.len(),
            // `myMap(myKey)` on an unbound key raises in OCCT; `Init` only sets
            // `myKey` for a bound key, so this arm is unreachable.
            None => false,
        }
    }

    /// `Next()` (`cxx:113-116`).
    pub fn next(&mut self) {
        self.current += 1;
    }

    /// `Value()` (`cxx:119-123`): `myMap(myKey).Value(current)`.
    pub fn value(&self) -> BasicEltId {
        let key = self
            .my_key
            .as_ref()
            .expect("BRepMAT2d_LinkTopoBilo::Value: no current key");
        let list = self
            .my_map
            .get(key)
            .expect("BRepMAT2d_LinkTopoBilo::Value: key is not bound");
        list[self.current - 1]
    }

    /// `GeneratingShape(BE)` (`cxx:127-130`): `myBEShape(BE)`.
    pub fn generating_shape(&self, be: BasicEltId) -> TopoShape {
        self.my_be_shape
            .get(&be)
            .unwrap_or_else(|| panic!("BRepMAT2d_LinkTopoBilo::GeneratingShape: not bound"))
            .clone()
    }

    /// `LinkToWire(W, Explo, IndC, BiLo)` (`cxx:137-210`).
    fn link_to_wire(
        &mut self,
        w: &Wire,
        explo: &BRepMat2dExplorer,
        ind_c: i32,
        bilo: &BRepMat2dBisectingLocus,
    ) {
        // `cxx:143`: `BRepTools_WireExplorer TheExp(W)` — no face.
        let mut exp = WireExplorer::of_wire(w);
        // `cxx:152-157`: `TopoSeq` collects `TheExp.Current()` in wire order.
        let mut topo_seq: Vec<TopoShape> = Vec::new();
        while exp.more() {
            let e = exp
                .current_edge()
                .expect("BRepMAT2d_LinkTopoBilo::LinkToWire: null current edge");
            topo_seq.push(e.0.clone());
            exp.next();
        }

        // `cxx:163`.
        let link_be_cont = link_to_contour(explo, ind_c, bilo);

        // `cxx:169-209`.
        for (be_index, kc) in link_be_cont {
            let be = bilo.graph().basic_elt(be_index);
            let basic_elt = bilo.graph().basic_elt_by_id(be);
            let is_point = matches!(bilo.geom_elt(basic_elt), Mat2dItem::Point(_));
            let kc_idx = kc.unsigned_abs() as usize;
            // `cxx:175`: `S = TopoSeq.Value(std::abs(KC))` — 1-based.
            let mut s = topo_seq
                .get(kc_idx - 1)
                .unwrap_or_else(|| panic!("BRepMAT2d_LinkTopoBilo::LinkToWire: TopoSeq out of range"))
                .clone();

            if is_point {
                // `cxx:178-192`: the point element belongs to one of the edge's
                // vertices. `TopExp::Vertices` runs with its default
                // `CumOri = false` (`TopExp.cxx:225`), i.e. it answers the
                // stored FORWARD / REVERSED children, and `cxx:178-187` swaps
                // the two output arguments when the edge itself is REVERSED -
                // so `VF` / `VL` below are the edge's *oriented* start / end
                // vertices, whichever way the wire traverses it.
                let (vf, vl) = edge_vertices(&Edge(s.clone()));
                let (walk_first, walk_last) = if s.orientation() == Orientation::Reversed {
                    (vl, vf)
                } else {
                    (vf, vl)
                };
                let vertex = if kc > 0 {
                    walk_last
                } else {
                    walk_first
                }
                .unwrap_or_else(|| {
                    panic!("BRepMAT2d_LinkTopoBilo::LinkToWire: edge without vertices")
                });
                s = vertex.0;
            }

            // `cxx:194-199`.
            self.my_map
                .entry(shape_map_key(&s))
                .or_default()
                .push(be);

            // `cxx:201-208`.
            let oriented = if kc < 0 {
                let reversed = if s.orientation() == Orientation::Reversed {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
                s.oriented(reversed)
            } else {
                s
            };
            self.my_be_shape.insert(be, oriented);
        }
    }
}

/// `LinkToContour(Explo, IndC, BiLo, Link)` (`cxx:217-322`): the index of the
/// contour curve each basic element comes from, negated when the element runs
/// against the contour direction.
fn link_to_contour(
    explo: &BRepMat2dExplorer,
    ind_c: i32,
    bilo: &BRepMat2dBisectingLocus,
) -> BTreeMap<i32, i32> {
    let mut link: BTreeMap<i32, i32> = BTreeMap::new();

    // `cxx:236`: `const NCollection_Sequence<handle(Geom2d_Curve)>& Cont =
    // Explo.Contour(IndC)`.
    let cont_len = explo.contour(ind_c as usize).len();

    // `cxx:240-243`.
    let mut ind_on_cont = 1i32;
    let mut prec_ind_on_cont = -1i32;
    let mut nb_sect = bilo.number_of_sections(ind_c, 1);
    let mut isect = 0i32;
    let mut direct_sense = true;
    let mut last_point = false;

    // `cxx:249-321`.
    for i in 1..=bilo.number_of_elts(ind_c) {
        let be = bilo.basic_elt(ind_c, i);
        let basic_elt = bilo.graph().basic_elt_by_id(be);
        let is_point = matches!(bilo.geom_elt(basic_elt), Mat2dItem::Point(_));

        if !is_point {
            isect += 1;
            // `cxx:261-269`.
            if direct_sense {
                link.insert(basic_elt.index(), ind_on_cont);
            } else {
                link.insert(basic_elt.index(), -ind_on_cont);
            }
        } else {
            // `cxx:274-282`.
            if direct_sense || last_point {
                link.insert(basic_elt.index(), prec_ind_on_cont);
            } else {
                link.insert(basic_elt.index(), -prec_ind_on_cont);
            }
        }

        prec_ind_on_cont = ind_on_cont;

        // `cxx:291-319`.
        if !is_point && isect == nb_sect {
            if ind_on_cont < cont_len as i32 && direct_sense {
                ind_on_cont += 1;
                nb_sect = bilo.number_of_sections(ind_c, ind_on_cont);
                isect = 0;
            } else {
                if !direct_sense {
                    ind_on_cont -= 1;
                    if ind_on_cont != 0 {
                        nb_sect = bilo.number_of_sections(ind_c, ind_on_cont);
                    }
                    last_point = false;
                } else {
                    last_point = true;
                }
                isect = 0;
                direct_sense = false;
            }
        }
    }
    link
}
