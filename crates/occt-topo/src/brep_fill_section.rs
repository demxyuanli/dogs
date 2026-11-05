//! `BRepFill_Section` — one cross-section of a swept shell.
//!
//! Source: `BRepFill_Section.cxx:25-137`, `.hxx:25-66`, `.lxx:15-50`.
//!
//! A section stores the caller's profile, the wire the sweep actually uses
//! (locations removed, `BRepFill_Section.cxx:46-71`), the vertex that law
//! functions take as the section's origin, and the `WithContact` /
//! `WithCorrection` flags of `BRepFill_Generator`'s section placement.
//!
//! A punctual profile (a single vertex) becomes a closed wire holding one
//! degenerated edge whose two vertices are that vertex in both orientations
//! (`cxx:55-70`) — the shape a sweep of a point needs to stay a legal wire.

use crate::abs::{Orientation, ShapeType};
use crate::builder::TopoBuilder;
use crate::iterator::cumulated_children;
use crate::iterator::ShapeExplorer;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::shape_upgrade::RemoveLocations;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::is_same;

/// `BRepFill_Section` (`BRepFill_Section.hxx:25-66`).
#[derive(Debug, Clone, Default)]
pub struct Section {
    /// `myOriginalShape` — the profile as the caller passed it, location and
    /// all. Null only for a default-constructed section.
    original_shape: Option<TopoShape>,
    /// `wire` — the profile with its locations removed.
    wire: Option<Wire>,
    /// `vertex` — the section vertex a law is evaluated at.
    vertex: Option<Vertex>,
    is_law: bool,
    is_punctual: bool,
    contact: bool,
    correction: bool,
}

impl Section {
    /// `BRepFill_Section()` (`cxx:25-33`): everything null / false.
    pub fn new() -> Self {
        Self::default()
    }

    /// `BRepFill_Section(Profile, V, WithContact, WithCorrection)`
    /// (`cxx:35-73`). A profile that is neither a wire nor a vertex is the
    /// `Standard_Failure` of `cxx:74`.
    pub fn with_profile(
        profile: &TopoShape,
        v: &Vertex,
        with_contact: bool,
        with_correction: bool,
    ) -> Result<Self, String> {
        let mut section = Self {
            original_shape: Some(profile.clone()),
            wire: None,
            vertex: Some(v.clone()),
            is_law: false,
            is_punctual: false,
            contact: with_contact,
            correction: with_correction,
        };

        // `cxx:48-50`: `SetRemoveLevel(TopAbs_COMPOUND)` then `Remove`.
        let mut rem_loc = RemoveLocations::new();
        rem_loc.set_remove_level(ShapeType::Compound);
        rem_loc.remove(profile);
        let a_profile = rem_loc.result().cloned().unwrap_or_else(|| profile.clone());

        // `cxx:52-71`.
        match a_profile.shape_type() {
            ShapeType::Wire => {
                section.wire = Wire::wrap(a_profile);
            }
            ShapeType::Vertex => {
                section.is_punctual = true;
                let Some(a_vertex) = Vertex::wrap(a_profile) else {
                    return Err(BAD_SECTION.to_string());
                };
                // `cxx:60-68`: one degenerated edge over the vertex in both
                // orientations.
                let builder = TopoBuilder::new();
                let mut deg_edge = Edge::new();
                let mut forward = a_vertex.0.clone();
                forward.set_orientation(Orientation::Forward);
                let mut reversed = a_vertex.0.clone();
                reversed.set_orientation(Orientation::Reversed);
                builder.add(&mut deg_edge.0, &forward);
                builder.add(&mut deg_edge.0, &reversed);
                // `BB.Degenerated(DegEdge, true)` (`BRep_Builder.cxx:1073-1087`):
                // the flag plus a `BRep_Curve3D` representation holding a null
                // curve, which the port records as `curve3d_null`.
                GeometryRegistry::global().set_degenerated(&deg_edge.0, true);
                GeometryRegistry::global().set_curve3d_null(&deg_edge.0, true);

                // `cxx:69-70`.
                let mut wire = Wire::new();
                builder.add_edge(&mut wire, &deg_edge);
                wire.0.set_closed(true);
                section.wire = Some(wire);
            }
            _ => return Err(BAD_SECTION.to_string()),
        }
        Ok(section)
    }

    /// `Set(IsLaw)` (`cxx:78-81`).
    pub fn set(&mut self, is_law: bool) {
        self.is_law = is_law;
    }

    /// `OriginalShape` (`.lxx:17-20`).
    pub fn original_shape(&self) -> Option<&TopoShape> {
        self.original_shape.as_ref()
    }

    /// `Wire` (`.lxx:22-25`).
    pub fn wire(&self) -> Option<&Wire> {
        self.wire.as_ref()
    }

    /// `Vertex` (`.lxx:27-30`).
    pub fn vertex(&self) -> Option<&Vertex> {
        self.vertex.as_ref()
    }

    /// `IsLaw` (`.lxx:32-35`).
    pub fn is_law(&self) -> bool {
        self.is_law
    }

    /// `IsPunctual` (`.lxx:37-40`).
    pub fn is_punctual(&self) -> bool {
        self.is_punctual
    }

    /// `WithContact` (`.lxx:42-45`).
    pub fn with_contact(&self) -> bool {
        self.contact
    }

    /// `WithCorrection` (`.lxx:47-50`).
    pub fn with_correction(&self) -> bool {
        self.correction
    }

    /// `ModifiedShape(theShape)` (`cxx:83-137`): the counterpart of an original
    /// sub-shape inside the location-free wire. `None` where OCCT answers a null
    /// shape.
    pub fn modified_shape(&self, the_shape: &TopoShape) -> Option<TopoShape> {
        let original = self.original_shape.as_ref()?;
        let wire = self.wire.as_ref()?;
        match the_shape.shape_type() {
            // `cxx:89-93`.
            ShapeType::Wire => {
                if is_same(the_shape, original) {
                    Some(wire.0.clone())
                } else {
                    None
                }
            }
            // `cxx:94-109`: walk both wires' children in parallel.
            ShapeType::Edge => {
                let originals = cumulated_children(original);
                let modified = cumulated_children(&wire.0);
                for (an_original_edge, an_edge) in originals.iter().zip(modified.iter()) {
                    if is_same(an_original_edge, the_shape) {
                        return Some(an_edge.clone());
                    }
                }
                None
            }
            // `cxx:110-131`.
            ShapeType::Vertex => {
                if is_same(the_shape, original) {
                    // `TopExp_Explorer(wire, TopAbs_VERTEX).Current()`.
                    return ShapeExplorer::new(cumulated_children(&wire.0), ShapeType::Vertex).next();
                }
                let exp_orig = ShapeExplorer::new(cumulated_children(original), ShapeType::Vertex);
                let exp_wire =
                    ShapeExplorer::new(cumulated_children(&wire.0), ShapeType::Vertex);
                for (an_original_vertex, a_vertex) in exp_orig.zip(exp_wire) {
                    if is_same(&an_original_vertex, the_shape) {
                        return Some(a_vertex);
                    }
                }
                None
            }
            // `cxx:132-133`.
            _ => None,
        }
    }
}

const BAD_SECTION: &str = "BRepFill_Section: bad shape type of section";
