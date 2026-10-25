//! Edge/face connexity tables of the classic boolean engine — a port of
//! `TKBool/TopOpeBRepDS/TopOpeBRepDS_connex.cxx` (433 lines).
//!
//! OCCT keeps the three tables in `thread_local` globals (`GLOBAL_elf1`,
//! `GLOBAL_elf2`, `GLOBAL_fle`, plus the `GLOBAL_los` empty-list answer and the
//! `GLOBAL_FDSCNX_prepared` flag) so that `FDSCNX_Prepare` can be called once
//! and the `FDSCNX_*` free functions queried many times without threading a
//! handle through every caller (`cxx:25-35`, the `NYI to CDLize as a tool of
//! DS` comment). The port keeps the same "prepare once, query many" shape but
//! makes the state explicit in [`Connexity`]; the free functions become
//! methods so no process-wide mutable state is introduced.
//!
//! Table contents (`FDSCNX_Prepare`, `cxx:127-220`): for every `FACE` of the
//! DS that carries an `AncestorRank` (`1` or `2`, `cxx:180-185`), each of its
//! `EDGE` sub-shapes that the DS knows (`BDS.HasShape`, `cxx:198-202`) is
//! appended to
//! - `fle[face]` — the edges of one face (`cxx:205-209`), and
//! - `elf{rank}[edge]` — the faces sharing one edge (`cxx:211-217`).
//!
//! Both appends happen once per *occurrence* yielded by `TopExp_Explorer`, so
//! a seam edge of a closed face appears twice in `fle[face]`.
//!
//! UNPORTED: `FDSCNX_DumpIndex` (`cxx:298-355`), `FDSCNX_Dump(HDS, I)`
//! (`cxx:356-419`) and `FDSCNX_Dump(HDS)` (`cxx:420-433`) — they print
//! `TopoDS_Shape` through `std::cout <<`, and the port has no stream formatter
//! for `TopoDS_Shape`; adding one would invent an output format OCCT does not
//! have.

use std::collections::HashMap;

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, is_same, shape_map_key, ShapeMapKey};

use super::h_data_structure::HDataStructure;

/// `TopTools_ShapeMapHasher` ↔ list of shapes, the `NCollection_DataMap`
/// used by the three tables.
type ConnexMap = HashMap<ShapeMapKey, Vec<TopoShape>>;

/// The `GLOBAL_elf1` / `GLOBAL_elf2` / `GLOBAL_fle` / `GLOBAL_los` /
/// `GLOBAL_FDSCNX_prepared` state of `TopOpeBRepDS_connex.cxx` (`cxx:25-35`).
#[derive(Default)]
pub struct Connexity {
    /// `GLOBAL_elf1` (`cxx:26`): edge -> faces of rank 1.
    elf1: ConnexMap,
    /// `GLOBAL_elf2` (`cxx:29`): edge -> faces of rank 2.
    elf2: ConnexMap,
    /// `GLOBAL_fle` (`cxx:32`): face -> its DS edges.
    fle: ConnexMap,
    /// `GLOBAL_FDSCNX_prepared` (`cxx:35`). OCCT also tests the three
    /// pointers for `nullptr` on the paths that read a table without going
    /// through `FDSCNX_Prepare` (`cxx:70-72`, `:246-249`, `:118-122`); since
    /// the maps always exist here, this flag is what those tests become.
    prepared: bool,
}

impl Connexity {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// `FDSCNX_Close()` (`cxx:41-58`): drop every table and the prepared
    /// flag.
    pub fn close(&mut self) {
        self.elf1.clear();
        self.elf2.clear();
        self.fle.clear();
        self.prepared = false;
    }

    /// `FDSCNX_Prepare(S1, S2, HDS)` (`cxx:127-220`). `S1`/`S2` are unused in
    /// OCCT (`cxx:127-128`).
    ///
    /// `HDS.IsNull()` (`cxx:131-135`) clears the flag: the port takes a
    /// borrowed `HDS`, so that arm is [`Connexity::close`] instead.
    pub fn prepare(&mut self, hds: &HDataStructure) {
        self.elf1.clear();
        self.elf2.clear();
        self.fle.clear();
        self.prepared = false;

        let bds = hds.ds();
        let n = bds.nb_shapes();
        for i in 1..=n {
            let Some(f) = bds.shape(i, true) else {
                continue;
            };
            if f.shape_type() != ShapeType::Face {
                continue;
            }
            let rf = bds.ancestor_rank_of(&f);
            if rf == 0 {
                continue;
            }
            // `TopExp_Explorer exe; for (exe.Init(f, TopAbs_EDGE); ...)`
            // (`cxx:191-192`): occurrences, not a set.
            let edges: Vec<TopoShape> = explorer_shapes(&f, ShapeType::Edge)
                .into_iter()
                .filter(|e| bds.has_shape(e, true))
                .collect();
            if edges.is_empty() {
                continue;
            }

            let f_key = shape_map_key(&f);
            self.fle
                .entry(f_key)
                .or_default()
                .extend(edges.iter().cloned());

            let elf = if rf == 1 {
                &mut self.elf1
            } else {
                &mut self.elf2
            };
            for e in &edges {
                elf.entry(shape_map_key(e)).or_default().push(f.clone());
            }
        }
        self.prepared = true;
    }

    /// `FDSCNX_EdgeConnexityShapeIndex(E, HDS, SI)` (`cxx:60-99`): the faces of
    /// rank `SI` (1 or 2) sharing the edge `E`.
    ///
    /// Every early return answers `GLOBAL_los`, the shared empty list
    /// (`cxx:67-91`), which is `&[]` here.
    pub fn edge_connexity_shape_index(
        &self,
        e: &TopoShape,
        hds: &HDataStructure,
        si: i32,
    ) -> &[TopoShape] {
        if !self.prepared {
            return &[];
        }
        if si != 1 && si != 2 {
            return &[];
        }
        if e.shape_type() != ShapeType::Edge {
            return &[];
        }
        if !self.has_connex_face(e, hds) {
            return &[];
        }
        if hds.ds().ancestor_rank_of(e) == 0 {
            return &[];
        }
        let elf = if si == 1 { &self.elf1 } else { &self.elf2 };
        match elf.get(&shape_map_key(e)) {
            Some(l) => l.as_slice(),
            None => &[],
        }
    }

    /// `FDSCNX_EdgeConnexitySameShape(S, HDS)` (`cxx:103-125`).
    ///
    /// `S = edge` answers the faces connected through `S`; `S = face` answers
    /// the edges of `S` that have at least one other connected face. Other
    /// types answer the empty list.
    pub fn edge_connexity_same_shape(
        &self,
        s: &TopoShape,
        hds: &HDataStructure,
    ) -> &[TopoShape] {
        match s.shape_type() {
            ShapeType::Edge => {
                let si = hds.ds().ancestor_rank_of(s);
                self.edge_connexity_shape_index(s, hds, si)
            }
            ShapeType::Face => {
                if !self.prepared {
                    return &[];
                }
                match self.fle.get(&shape_map_key(s)) {
                    Some(l) => l.as_slice(),
                    None => &[],
                }
            }
            _ => &[],
        }
    }

    /// `FDSCNX_HasConnexFace(S, HDS)` (`cxx:222-251`).
    pub fn has_connex_face(&self, s: &TopoShape, hds: &HDataStructure) -> bool {
        let t = s.shape_type();
        if t != ShapeType::Face && t != ShapeType::Edge {
            return false;
        }
        let rs = hds.ds().ancestor_rank_of(s);
        if rs == 0 {
            return false;
        }
        if !self.prepared {
            return false;
        }
        let elf = if rs == 1 { &self.elf1 } else { &self.elf2 };
        let key = shape_map_key(s);
        if t == ShapeType::Edge {
            elf.contains_key(&key)
        } else {
            self.fle.contains_key(&key)
        }
    }

    /// `FDSCNX_FaceEdgeConnexFaces(F, E, HDS, LF)` (`cxx:253-296`): the faces
    /// sharing `E` with `F`, `F` itself excluded.
    ///
    /// The first loop (`cxx:262-280`) only checks that `E` is one of `F`'s
    /// connectivity edges; the answer is built from the second loop
    /// (`cxx:282-295`).
    pub fn face_edge_connex_faces(
        &self,
        f: &TopoShape,
        e: &TopoShape,
        hds: &HDataStructure,
        lf: &mut Vec<TopoShape>,
    ) {
        lf.clear();
        let loe = self.edge_connexity_same_shape(f, hds);
        if loe.is_empty() {
            return;
        }
        if !loe.iter().any(|x| is_same(x, e)) {
            return;
        }
        let lof = self.edge_connexity_same_shape(e, hds);
        if lof.is_empty() {
            return;
        }
        for other in lof {
            if !is_same(other, f) {
                lf.push(other.clone());
            }
        }
    }
}
