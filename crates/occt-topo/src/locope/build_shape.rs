//! `LocOpe_BuildShape` (`LocOpe_BuildShape.hxx:25-46`, `.cxx:54-463`).
//!
//! Rebuilds solids and shells out of a loose list of faces: the faces are
//! connected through their shared edges into shells (`Add` / `Propagate`), each
//! shell is oriented against the material side, accepted as a solid when every
//! one of its edges is manifold (exactly two faces) or degenerated, and the
//! accepted solids are finally nested one level deep (`IsInside`).
//!
//! Module table (`ModelingAlgorithms/TKFeat/LocOpe`):
//!
//! | OCCT | Rust |
//! |---|---|
//! | ctor / `Perform` (`cxx:54-320`) / `Shape` (`lxx`) | this file |
//! | `static Add` (`cxx:323-362`) | [`add`] |
//! | `static Propagate` (`cxx:367-436`) | [`propagate`] |
//! | `static IsInside` (`cxx:441-463`) | [`is_inside`] |
//!
//! Container mapping:
//! - `mapF` (`NCollection_IndexedMap<TopoDS_Shape>`, 1-based) is
//!   [`ShapeIndexedMap`] holding the face itself, so `mapF(j)` is
//!   [`ShapeIndexedMap::find_key`] and `mapF.FindIndex(F)` is
//!   [`ShapeIndexedMap::find_index`].
//! - `mapI` / `mapIf` (`NCollection_Map<int>`) are `BTreeSet<usize>`. The
//!   reference walks hash buckets over the integer value; every insertion is
//!   ascending and `Propagate` never inserts, so ascending order is the only
//!   order those sets ever have.
//! - `theMapEF` / `theMapEFbis` (`TopExp::MapShapesAndAncestors`) are
//!   [`AncestorMap`].
//! - `mapSh` (`cxx:88`) is declared and never read in the reference, so it is
//!   not built here either.
//! - `imbSh` (`NCollection_DataMap<Shape, List<Shape>>`) is
//!   [`ShapeValueMap`].
//!
//! Deviations, all forced by the language or by an OCCT raise:
//! - `myRes` has no null state, so it is `Option<TopoShape>`; `Shape()` is
//!   `None` until `Perform` produces a shape.
//! - The two `Standard_ConstructionError` raises of `Add` (`cxx:333`,
//!   `cxx:353`) become the `Err` of [`Perform`](BuildShape::perform).

use std::collections::BTreeSet;

use crate::abs::{Orientation, ShapeType};
use crate::brep_class3d::SolidClassifier;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::fclass2d::FaceState;
use crate::locope::shape_maps::{ShapeIndexedMap, ShapeValueMap};
use crate::shape::{Shell, Solid, TopoShape, Vertex};
use crate::topo_tools::{explorer_shapes, is_same, AncestorMap};

/// `LocOpe_BuildShape` (`hxx:28-46`).
#[derive(Debug, Clone, Default)]
pub struct BuildShape {
    /// `TopoDS_Shape myRes`
    my_res: Option<TopoShape>,
}

impl BuildShape {
    /// `LocOpe_BuildShape()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// `LocOpe_BuildShape(L)`: the workload constructor.
    pub fn with_faces(&mut self, l: &[TopoShape]) -> Result<(), String> {
        self.perform(l)
    }

    /// `Perform(L)` (`cxx:54-320`): build shape(s) from the faces of `L`.
    pub fn perform(&mut self, l: &[TopoShape]) -> Result<(), String> {
        // `cxx:60`: `myRes.Nullify()`.
        self.my_res = None;

        let builder = TopoBuilder::new();

        // `cxx:66-77`: the compound of the distinct faces, added as they are.
        let mut map_f: ShapeIndexedMap<()> = ShapeIndexedMap::new();
        let mut faces: Vec<TopoShape> = Vec::new();
        for s in l {
            if s.shape_type() != ShapeType::Face || map_f.contains(s) {
                continue;
            }
            map_f.add(s, ());
            faces.push(s.clone());
        }

        // `cxx:79-82`: no face, no shape.
        if map_f.extent() == 0 {
            return Ok(());
        }

        // `cxx:86`: `TopExp::MapShapesAndAncestors(C, EDGE, FACE, theMapEF)`.
        let compound = builder.make_compound_of(&faces);
        let map_ef = AncestorMap::map_shapes_and_ancestors(&compound.0, ShapeType::Edge, ShapeType::Face);
        let nbedges = map_ef.length();

        // `cxx:90-94`.
        let mut map_i: BTreeSet<usize> = BTreeSet::new();
        let mut map_if: BTreeSet<usize> = BTreeSet::new();
        let mut lshell: Vec<TopoShape> = Vec::new();
        let mut lresult: Vec<TopoShape> = Vec::new();

        // `cxx:96-210`: one shell per connected group of faces.
        loop {
            // `cxx:98-102`: the first edge index not processed yet.
            let mut i = 1usize;
            while i <= nbedges {
                if !map_i.contains(&i) {
                    break;
                }
                i += 1;
            }
            if i <= nbedges {
                // `cxx:106-108`.
                map_f.clear();
                map_if.clear();
                add(i, &mut map_i, &mut map_f, &map_ef)?;

                // `cxx:110-126`: manifoldness and the reference face.
                let mut manifold = true;
                let mut face_ref: Option<TopoShape> = None;
                for j in 1..=map_f.extent() {
                    let orient = map_f.find_key(j).orientation();
                    if orient == Orientation::Internal || orient == Orientation::External {
                        manifold = false;
                    } else if face_ref.is_none() {
                        face_ref = Some(map_f.find_key(j).clone());
                    }
                    map_if.insert(j);
                }

                let mut new_sh = Shell::new();
                match (&face_ref, manifold) {
                    // `cxx:130-138`: a bundle of unoriented faces.
                    (None, _) => {
                        for j in 1..=map_f.extent() {
                            builder.add(&mut new_sh.0, map_f.find_key(j));
                        }
                    }
                    (Some(face_ref), _) => {
                        // `cxx:141-152`: the unorientable faces go in as they are.
                        if !manifold {
                            for j in 1..=map_f.extent() {
                                let o = map_f.find_key(j).orientation();
                                if o == Orientation::Internal || o == Orientation::External {
                                    builder.add(&mut new_sh.0, map_f.find_key(j));
                                    map_if.remove(&j);
                                }
                            }
                        }
                        // `cxx:154-155`.
                        builder.add(&mut new_sh.0, face_ref);
                        propagate(face_ref, &mut new_sh.0, &map_f, &mut map_if);
                    }
                }

                // `cxx:157`: `newSh.Closed(BRep_Tool::IsClosed(newSh))`.
                let closed = BRepTool::is_closed(&new_sh.0);
                new_sh.0.set_closed(closed);

                if !manifold {
                    // `cxx:158-161`.
                    lshell.push(new_sh.0.oriented(Orientation::Internal));
                } else {
                    // `cxx:162-190`: every edge must be manifold or degenerated.
                    let map_ef_bis = AncestorMap::map_shapes_and_ancestors(
                        &new_sh.0,
                        ShapeType::Edge,
                        ShapeType::Face,
                    );
                    let mut ok = true;
                    for k in 1..=map_ef_bis.length() {
                        let Some(ed) = crate::shape::Edge::wrap(map_ef_bis.find_key(k).clone()) else {
                            continue;
                        };
                        let ori_ed = ed.0.orientation();
                        if ori_ed == Orientation::Internal || ori_ed == Orientation::External {
                            continue;
                        }
                        let nbfac = map_ef_bis.ancestors(k).len();
                        if nbfac > 2 {
                            // `cxx:176-179`: "peu probable".
                            ok = false;
                            break;
                        }
                        if nbfac == 1 && !BRepTool::is_degenerated(&ed) {
                            ok = false;
                            break;
                        }
                    }

                    if ok {
                        // `cxx:190-204`: a closed, manifold shell is a solid; the
                        // classifier decides whether it must be reversed.
                        let new_so = builder.make_solid(&[new_sh.clone()]);
                        let mut class = SolidClassifier::new();
                        class.load(new_so.0.clone());
                        class.perform_infinite_point(occt_core::precision::CONFUSION);
                        if class.state() == FaceState::In {
                            lresult.push(new_sh.0.oriented(Orientation::Reversed));
                        } else {
                            lresult.push(new_sh.0.clone());
                        }
                    } else {
                        // `cxx:205-208`.
                        lshell.push(new_sh.0.oriented(Orientation::Internal));
                    }
                }
            }

            // `cxx:210`: `} while (mapI.Extent() < Nbedges);`
            if map_i.len() >= nbedges {
                break;
            }
        }

        // `cxx:212-253`: one level of nesting between the built solids.
        let mut imb_sh: ShapeValueMap<Vec<TopoShape>> = ShapeValueMap::new();
        let mut l_intern: Vec<TopoShape> = Vec::new();
        for sh in &lresult {
            // `cxx:223-224`: a throwaway solid carrying just this shell.
            let tempo = builder.make_solid(&[Shell(sh.clone())]);
            imb_sh.bind(sh, Vec::new());
            for sh2 in &lresult {
                if is_same(sh2, sh) {
                    continue;
                }
                if is_inside(sh2, &tempo.0) {
                    l_intern.push(sh2.clone());
                    if let Some(v) = imb_sh.change(sh) {
                        v.push(sh2.clone());
                    }
                }
            }
        }

        // `cxx:247-253`: the nested shells are not solids of their own.
        for sh in &l_intern {
            if imb_sh.is_bound(sh) {
                imb_sh.unbind(sh);
            }
        }

        // `cxx:257-295`: solids out of the nesting table.
        let mut lsolid: Vec<TopoShape> = Vec::new();
        while imb_sh.extent() != 0 {
            // `cxx:263-271`: the first entry with a non-empty nested list.
            let first_with_inner = imb_sh
                .iter()
                .find(|(_, v)| !v.is_empty())
                .map(|(k, v)| (k.clone(), v.clone()));
            match first_with_inner {
                Some((key, inner)) => {
                    // `cxx:272-282`.
                    let mut new_so = Solid::new();
                    builder.add(&mut new_so.0, &key);
                    for it in &inner {
                        let mut r = it.clone();
                        r.reverse();
                        builder.add(&mut new_so.0, &r);
                    }
                    lsolid.push(new_so.0);
                    imb_sh.unbind(&key);
                }
                None => {
                    // `cxx:283-293`.
                    let keys: Vec<TopoShape> = imb_sh.iter().map(|(k, _)| k.clone()).collect();
                    for k in keys {
                        let mut new_so = Solid::new();
                        builder.add(&mut new_so.0, &k);
                        lsolid.push(new_so.0);
                    }
                    imb_sh.clear();
                }
            }
        }

        // `cxx:297-319`.
        let nbsol = lsolid.len();
        let nbshl = lshell.len();
        if nbsol == 1 && nbshl == 0 {
            self.my_res = Some(lsolid[0].clone());
        } else if nbsol == 0 && nbshl == 1 {
            self.my_res = Some(lshell[0].clone());
        } else {
            let mut items: Vec<TopoShape> = Vec::new();
            items.extend(lsolid.iter().cloned());
            items.extend(lshell.iter().cloned());
            self.my_res = Some(builder.make_compound_of(&items).0);
        }
        Ok(())
    }

    /// `Shape()` (`lxx`): `None` is the null `myRes` before a shape is built.
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_res.as_ref()
    }
}

/// `static void Add(ind, mapI, mapF, mapEF)` (`cxx:323-362`): flood over the
/// faces connected to edge `ind` (1-based), collecting them in `map_f`.
fn add(
    ind: usize,
    map_i: &mut BTreeSet<usize>,
    map_f: &mut ShapeIndexedMap<()>,
    map_ef: &AncestorMap,
) -> Result<(), String> {
    // `cxx:332-335`: a repeated index is `Standard_ConstructionError`.
    if !map_i.insert(ind) {
        return Err("LocOpe_BuildShape: edge index visited twice".to_string());
    }

    // `cxx:337-360`.
    for face in map_ef.ancestors(ind - 1) {
        if map_f.contains(face) {
            continue;
        }
        map_f.add(face, ());
        for edg in explorer_shapes(face, ShapeType::Edge) {
            // `cxx:351-355`: the edge must be in the map, else
            // `Standard_ConstructionError`.
            let Some(indedg) = map_ef.find_index(&edg) else {
                return Err("LocOpe_BuildShape: edge not in the ancestor map".to_string());
            };
            let indedg = indedg + 1;
            if !map_i.contains(&indedg) {
                add(indedg, map_i, map_f, map_ef)?;
            }
        }
    }
    Ok(())
}

/// `static void Propagate(F, Sh, mapF, mapIf)` (`cxx:367-436`): orient `F`'s
/// neighbours against the shared edge and add them to `Sh`, recursively.
fn propagate(
    f: &TopoShape,
    sh: &mut TopoShape,
    map_f: &ShapeIndexedMap<()>,
    map_if: &mut BTreeSet<usize>,
) {
    let builder = TopoBuilder::new();
    // `cxx:373-378`.
    let Some(indf) = map_f.find_index(f) else {
        return;
    };
    if !map_if.contains(&indf) {
        return;
    }
    map_if.remove(&indf);
    // `cxx:379-382`.
    if map_if.is_empty() {
        return;
    }

    // `cxx:384-435`.
    for edg in explorer_shapes(f, ShapeType::Edge) {
        let ored1 = edg.orientation();
        // `cxx:390-393`.
        if ored1 == Orientation::Internal || ored1 == Orientation::External {
            continue;
        }
        let mut ored2 = Orientation::Forward;
        let mut found: Option<TopoShape> = None;
        // `cxx:396-412`: the first pending face carrying the same edge.
        for &key in map_if.iter() {
            let new_f = map_f.find_key(key);
            let mut hit = false;
            for e2 in explorer_shapes(new_f, ShapeType::Edge) {
                if is_same(&e2, &edg) {
                    ored2 = e2.orientation();
                    hit = true;
                    break;
                }
            }
            if hit {
                found = Some(new_f.clone());
                break;
            }
        }

        // `cxx:413-435`.
        if let Some(mut ftoadd) = found {
            let mut added = false;
            if ored2 == ored1 {
                // `cxx:418-423`: same orientation on both sides - reverse.
                ftoadd.reverse();
                builder.add(sh, &ftoadd);
                added = true;
            } else if ored2 == ored1.reversed() {
                // `cxx:425-429`.
                builder.add(sh, &ftoadd);
                added = true;
            }
            if added {
                propagate(&ftoadd, sh, map_f, map_if);
            }
        }
    }
}

/// `static bool IsInside(S1, S2)` (`cxx:441-463`): true when the first vertex
/// of `S1` that the classifier does not put ON lies IN `S2`. The trailing
/// `return true` is the reference's answer when every vertex is ON.
fn is_inside(s1: &TopoShape, s2: &TopoShape) -> bool {
    let mut class = SolidClassifier::new();
    class.load(s2.clone());
    for s in explorer_shapes(s1, ShapeType::Vertex) {
        let vtx = Vertex(s);
        // `cxx:448-450`: `BRep_Tool::Pnt` / `Tolerance` then `Class.Perform`.
        let pt = BRepTool::vertex_point_world(&vtx);
        let tol = BRepTool::vertex_tolerance(&vtx);
        class.perform(&pt, tol);
        match class.state() {
            FaceState::In => return true,
            FaceState::Out => return false,
            _ => {}
        }
    }
    true
}
