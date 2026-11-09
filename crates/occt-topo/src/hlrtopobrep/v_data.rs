//! `HLRTopoBRep_VData` (`ModelingAlgorithms/TKHLR/HLRTopoBRep/HLRTopoBRep_VData.hxx:26-43`,
//! `HLRTopoBRep_VData.cxx:22-28`, `HLRTopoBRep_VData.lxx:19-36`): a vertex of an
//! edge plus the parameter at which it was inserted while splitting the edge
//! for the hidden-line removal.
//!
//! OCCT stores a `TopoDS_Shape` that defaults to the null shape; the port
//! stores `Option<TopoShape>` so "no vertex yet" stays expressible.

use crate::shape::TopoShape;

/// `HLRTopoBRep_VData`.
#[derive(Debug, Clone, Default)]
pub struct HlrTopoBRepVData {
    /// `double myParameter` (`hxx:40`).
    my_parameter: f64,
    /// `TopoDS_Shape myVertex` (`hxx:41`).
    my_vertex: Option<TopoShape>,
}

impl HlrTopoBRepVData {
    /// `HLRTopoBRep_VData()` (`lxx:19-22`): `myParameter = 0.0`, null vertex.
    pub fn new() -> Self {
        Self::default()
    }

    /// `HLRTopoBRep_VData(P, V)` (`cxx:22-28`).
    pub fn with_vertex(parameter: f64, vertex: &TopoShape) -> Self {
        Self {
            my_parameter: parameter,
            my_vertex: Some(vertex.clone()),
        }
    }

    /// `Parameter()` (`lxx:26-29`).
    pub fn parameter(&self) -> f64 {
        self.my_parameter
    }

    /// `Vertex()` (`lxx:33-36`).
    pub fn vertex(&self) -> Option<&TopoShape> {
        self.my_vertex.as_ref()
    }
}
