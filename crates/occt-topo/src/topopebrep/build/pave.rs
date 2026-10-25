//! `TopOpeBRepBuild_Pave` — a parametrized vertex on an edge.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Pave.hxx:30-78`,
//! `.cxx:22-112`.
//!
//! `TopOpeBRepBuild_Pave` is a `TopOpeBRepBuild_Loop` subclass; the port keeps
//! the subclass fields in [`PaveData`] and stores them in the
//! [`LoopKind::Pave`](super::loop_shape::LoopKind::Pave) payload of a
//! [`Loop`] (see `loop_shape.rs`).

use crate::shape::TopoShape;
use crate::topopebrep::ds::kind::Kind;

use super::loop_shape::{Loop, LoopKind};

/// The `TopOpeBRepBuild_Pave` fields (`hxx:72-77`).
#[derive(Debug, Clone)]
pub struct PaveData {
    /// `myVertex` (`hxx:72`).
    pub vertex: TopoShape,
    /// `myParam` (`hxx:73`) — the parameter of `myVertex` on the edge.
    pub param: f64,
    /// `myIsShape` (`hxx:74`) — `bound` of the constructor: true when the
    /// vertex is an old (edge) vertex, false when it is a new one.
    pub is_shape: bool,
    /// `myHasSameDomain` (`hxx:75`).
    pub has_same_domain: bool,
    /// `mySameDomain` (`hxx:76`), null until `SameDomain(VSD)` is called.
    pub same_domain: Option<TopoShape>,
    /// `myIntType` (`hxx:77`), `TopOpeBRepDS_FACE` at construction.
    pub int_type: Kind,
}

impl Loop {
    /// `TopOpeBRepBuild_Pave(V, P, bound)` (`cxx:26-34`).
    ///
    /// The constructor first runs the `TopOpeBRepBuild_Loop(V)` base
    /// constructor (`cxx:27`), which sets `myIsShape = true`,
    /// `myShape = V` and `myBlockIterator(0, 0)`.
    pub fn new_pave(v: &TopoShape, p: f64, bound: bool) -> Self {
        let mut l = Loop::new(v);
        l.kind = LoopKind::Pave(PaveData {
            vertex: v.clone(),
            param: p,
            is_shape: bound,
            has_same_domain: false,
            same_domain: None,
            int_type: Kind::Face,
        });
        l
    }
}

impl PaveData {
    /// `HasSameDomain(B)` (`cxx:38-41`).
    pub fn set_has_same_domain(&mut self, b: bool) {
        self.has_same_domain = b;
    }

    /// `SameDomain(VSD)` (`cxx:45-48`).
    pub fn set_same_domain(&mut self, vsd: &TopoShape) {
        self.same_domain = Some(vsd.clone());
    }

    /// `HasSameDomain()` (`cxx:52-55`).
    pub fn has_same_domain(&self) -> bool {
        self.has_same_domain
    }

    /// `SameDomain()` (`cxx:59-62`).
    pub fn same_domain(&self) -> Option<&TopoShape> {
        self.same_domain.as_ref()
    }

    /// `Vertex()` (`cxx:66-69`).
    pub fn vertex(&self) -> &TopoShape {
        &self.vertex
    }

    /// `ChangeVertex()` (`cxx:73-76`).
    pub fn change_vertex(&mut self) -> &mut TopoShape {
        &mut self.vertex
    }

    /// `Parameter()` (`cxx:80-83`).
    pub fn parameter(&self) -> f64 {
        self.param
    }

    /// `Parameter(Par)` (`cxx:88-91`).
    pub fn set_parameter(&mut self, par: f64) {
        self.param = par;
    }

    /// `InterferenceType()` (`cxx:107-110`) — returns `TopOpeBRepDS_Kind&`.
    pub fn interference_type_mut(&mut self) -> &mut Kind {
        &mut self.int_type
    }

    /// `InterferenceType()`, read-only view.
    pub fn interference_type(&self) -> Kind {
        self.int_type
    }
}
