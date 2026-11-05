//! `ChFi3d_Builder` result / fault / contour accessors.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_1.cxx:384-623`,
//! `:1181-1356` (OCCT 8.0.0). Split out of `builder/mod.rs` so the
//! orchestration layer (`Compute` and friends) can grow without pushing the
//! module over the file-length limit.

use std::sync::Arc;

use occt_core::kernel::geomabs::Shape;
use occt_geom::Surface as GeomSurface;

use super::{edge_index_in_spine, is_tangent_faces, spine_of, ChFi3dBuilder};
use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::chfi_ds::{ChFiDSMap, ErrorStatus, SpineHandle, StripeHandle};
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::is_same;
use crate::topopebrep::build::HBuilder;

impl ChFi3dBuilder {
    /// `IsDone()` (`ChFi3d_Builder_1.cxx:384-387`).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// `fleche` (`ChFi3d_Builder.hxx:376`), the deflection `SetParams` stores
    /// (`ChFi3d_Builder_1.cxx:358-372`). OCCT reads the protected member
    /// directly from `PerformTwoCorner` / `PerformThreeCorner`; the port exposes
    /// it to those sibling modules.
    pub fn fleche(&self) -> f64 {
        self.fleche
    }

    /// `Shape()` (`ChFi3d_Builder_1.cxx:391-395`): `Standard_NoSuchObject` when
    /// the computation did not run, which the port answers as `Err`; the
    /// returned `Option` is the OCCT null-shape case.
    pub fn shape(&self) -> Result<Option<TopoShape>, String> {
        if !self.done {
            return Err("ChFi3d_Builder::Shape() - no result".to_string());
        }
        Ok(self.my_shape_result.clone())
    }

    /// `NbFaultyContours()` (`ChFi3d_Builder_1.cxx:399-402`).
    pub fn nb_faulty_contours(&self) -> usize {
        self.bad_stripes.len()
    }

    /// `FaultyContour(I)` (`ChFi3d_Builder_1.cxx:405-434`): the `I`-th faulty
    /// contour is the position of its stripe in `myListStripe`, `0` when the
    /// stripe is not there (the OCCT sentinel).
    pub fn faulty_contour(&self, i: usize) -> usize {
        let Some(stripe) = i.checked_sub(1).and_then(|k| self.bad_stripes.get(k)) else {
            return 0;
        };
        for (k, candidate) in self.my_list_stripe.iter().enumerate() {
            if Arc::ptr_eq(candidate, stripe) {
                return k + 1;
            }
        }
        0
    }

    /// `NbComputedSurfaces(IC)` (`ChFi3d_Builder_1.cxx:438-467`).
    pub fn nb_computed_surfaces(&self, ic: usize) -> usize {
        let Some(stripe) = self.stripe_at(ic) else {
            return 0;
        };
        let Ok(guard) = stripe.read() else {
            return 0;
        };
        if guard.spine().is_none() {
            return 0;
        }
        guard.set_of_surf_data().len()
    }

    /// `ComputedSurface(IC, IS)` (`ChFi3d_Builder_1.cxx:470-488`): the surface
    /// of the `IS`-th `SurfData` of the `IC`-th contour, read through the DS.
    pub fn computed_surface(&self, ic: usize, is: usize) -> Option<Arc<dyn GeomSurface>> {
        let stripe = self.stripe_at(ic)?;
        let guard = stripe.read().ok()?;
        let data = guard.set_of_surf_data().get(is.checked_sub(1)?)?;
        self.my_ds.surface(data.surf()).surface().cloned()
    }

    /// `NbFaultyVertices()` (`ChFi3d_Builder_1.cxx:491-494`).
    pub fn nb_faulty_vertices(&self) -> usize {
        self.bad_vertices.len()
    }

    /// `FaultyVertex(IV)` (`ChFi3d_Builder_1.cxx:497-513`), `IV` `1`-based.
    pub fn faulty_vertex(&self, iv: usize) -> Option<Vertex> {
        let shape = iv.checked_sub(1).and_then(|k| self.bad_vertices.get(k))?;
        Vertex::wrap(shape.clone())
    }

    /// `HasResult()` (`ChFi3d_Builder_1.cxx:517-520`).
    pub fn has_result(&self) -> bool {
        self.has_result
    }

    /// `BadShape()` (`ChFi3d_Builder_1.cxx:524-528`): `Standard_NoSuchObject`
    /// when no partial result was stored, which the port answers as `Err`.
    pub fn bad_shape(&self) -> Result<Option<TopoShape>, String> {
        if !self.has_result {
            return Err("ChFi3d_Builder::BadShape() - no result".to_string());
        }
        Ok(self.bad_shape.clone())
    }

    /// `StripeStatus(IC)` (`ChFi3d_Builder_1.cxx:532-548`).
    pub fn stripe_status(&self, ic: usize) -> Option<ErrorStatus> {
        let stripe = self.stripe_at(ic)?;
        let spine = spine_of(&stripe)?;
        let guard = spine.read().ok()?;
        Some(guard.error_status())
    }

    /// `Remove(E)` (`ChFi3d_Builder_1.cxx:1181-1197`): drops the first stripe
    /// one of whose edges is `E`.
    pub fn remove(&mut self, e: &Edge) {
        let mut found = None;
        for (k, stripe) in self.my_list_stripe.iter().enumerate() {
            let Some(spine) = spine_of(stripe) else {
                continue;
            };
            if edge_index_in_spine(&spine, &e.0).is_some() {
                found = Some(k);
                break;
            }
        }
        if let Some(k) = found {
            self.my_list_stripe.remove(k);
        }
    }

    /// `Value(I)` (`ChFi3d_Builder_1.cxx:1201-1210`).
    pub fn value(&self, i: usize) -> Option<SpineHandle> {
        let stripe = self.stripe_at(i)?;
        spine_of(&stripe)
    }

    /// `NbElements()` (`ChFi3d_Builder_1.cxx:1213-1226`): the stripes up to the
    /// first one whose `Spine()` is null.
    pub fn nb_elements(&self) -> usize {
        let mut n = 0;
        for stripe in &self.my_list_stripe {
            if spine_of(stripe).is_none() {
                break;
            }
            n += 1;
        }
        n
    }

    /// `Contains(E)` (`ChFi3d_Builder_1.cxx:1231-1252`): the `1`-based index of
    /// the contour holding `E`, `0` when no contour does.
    pub fn contains(&self, e: &Edge) -> usize {
        self.contains_with_index(e).0
    }

    /// `Contains(E, IndexInSpine)` (`ChFi3d_Builder_1.cxx:1255-1277`): the
    /// contour index and the `1`-based index of `E` inside that contour.
    pub fn contains_with_index(&self, e: &Edge) -> (usize, usize) {
        for (k, stripe) in self.my_list_stripe.iter().enumerate() {
            let Some(spine) = spine_of(stripe) else {
                break;
            };
            if let Some(j) = edge_index_in_spine(&spine, &e.0) {
                return (k + 1, j);
            }
        }
        (0, 0)
    }

    /// `Length(IC)` (`ChFi3d_Builder_1.cxx:1281-1290`):
    /// `Value(IC)->LastParameter(NbEdges())`. The `-1` of an out-of-range `IC`
    /// and the `Standard_OutOfRange` of a not-yet-`Load`ed guideline both
    /// answer `None`.
    pub fn length(&self, ic: usize) -> Option<f64> {
        if ic == 0 || ic > self.nb_elements() {
            return None;
        }
        let spine = self.value(ic)?;
        let guard = spine.read().ok()?;
        guard.last_parameter_of(guard.nb_edges() as i32)
    }

    /// `FirstVertex(IC)` (`ChFi3d_Builder_1.cxx:1293-1300`).
    pub fn first_vertex(&self, ic: usize) -> Option<TopoShape> {
        if ic == 0 || ic > self.nb_elements() {
            return None;
        }
        let spine = self.value(ic)?;
        let guard = spine.read().ok()?;
        guard.first_vertex()
    }

    /// `LastVertex(IC)` (`ChFi3d_Builder_1.cxx:1303-1311`).
    pub fn last_vertex(&self, ic: usize) -> Option<TopoShape> {
        if ic == 0 || ic > self.nb_elements() {
            return None;
        }
        let spine = self.value(ic)?;
        let guard = spine.read().ok()?;
        guard.last_vertex()
    }

    /// `Abscissa(IC, V)` (`ChFi3d_Builder_1.cxx:1315-1322`): `-1.` both when
    /// `IC` is out of range and when `V` is not on the contour, as in OCCT.
    pub fn abscissa(&self, ic: usize, v: &Vertex) -> f64 {
        if ic == 0 || ic > self.nb_elements() {
            return -1.0;
        }
        let Some(spine) = self.value(ic) else {
            return -1.0;
        };
        let Ok(guard) = spine.read() else {
            return -1.0;
        };
        guard.absc_of_vertex(&v.0)
    }

    /// `RelativeAbscissa(IC, V)` (`ChFi3d_Builder_1.cxx:1325-1334`).
    pub fn relative_abscissa(&self, ic: usize, v: &Vertex) -> f64 {
        if ic == 0 || ic > self.nb_elements() {
            return -1.0;
        }
        self.abscissa(ic, v) / self.length(ic).unwrap_or(-1.0)
    }

    /// `Closed(IC)` (`ChFi3d_Builder_1.cxx:1337-1345`).
    pub fn closed(&self, ic: usize) -> bool {
        if ic == 0 || ic > self.nb_elements() {
            return false;
        }
        let Some(spine) = self.value(ic) else {
            return false;
        };
        let Ok(guard) = spine.read() else {
            return false;
        };
        // `IsClosed()` raises on a guideline without edges; the port answers
        // `None`, which is "not closed" here.
        guard.is_closed().unwrap_or(false)
    }

    /// `ClosedAndTangent(IC)` (`ChFi3d_Builder_1.cxx:1348-1355`): OCCT asks
    /// `IsPeriodic()`, i.e. `firstState == ChFiDS_Closed`.
    pub fn closed_and_tangent(&self, ic: usize) -> bool {
        if ic == 0 || ic > self.nb_elements() {
            return false;
        }
        let Some(spine) = self.value(ic) else {
            return false;
        };
        let Ok(guard) = spine.read() else {
            return false;
        };
        guard.is_periodic()
    }

    /// `Builder()` (`ChFi3d_Builder_1.cxx:552-555`).
    pub fn builder(&self) -> &HBuilder {
        &self.my_coup
    }

    /// The protected `myEFMap` (`hxx:750`), read by the subclasses
    /// (`ChFi3d_ChBuilder.cxx:205`, `ChFi3d_FilBuilder_C2.cxx:434`).
    pub fn ef_map(&self) -> &ChFiDSMap {
        &self.my_ef_map
    }

    /// The protected `myESoMap` (`hxx:751`).
    pub fn eso_map(&self) -> &ChFiDSMap {
        &self.my_eso_map
    }

    /// The protected `myEShMap` (`hxx:752`).
    pub fn esh_map(&self) -> &ChFiDSMap {
        &self.my_esh_map
    }

    /// The protected `myVFMap` (`hxx:753`).
    pub fn vf_map(&self) -> &ChFiDSMap {
        &self.my_vf_map
    }

    /// The protected `myShape` (`hxx:842`), the input shape of the builder.
    pub fn input_shape(&self) -> &TopoShape {
        &self.my_shape
    }

    /// The protected `myListStripe` (`hxx:755`), in creation order.
    ///
    /// OCCT reads it directly from the API classes that wrap a builder
    /// (`BRepFilletAPI_MakeFillet.cxx`, `FilletSurf_InternalBuilder.cxx`); the
    /// port exposes the slice to those siblings.
    pub fn stripes(&self) -> &[StripeHandle] {
        &self.my_list_stripe
    }

    /// The protected `tolesp` (`hxx:744`), the tolerance on the spine.
    pub fn tolesp(&self) -> f64 {
        self.tolesp
    }

    /// The protected `tolapp3d` (`hxx:747`), the 3d approximation tolerance.
    pub fn tolapp3d(&self) -> f64 {
        self.tolapp3d
    }

    /// The protected `myConti` (`hxx:748`), the internal continuity.
    pub fn continuity(&self) -> Shape {
        self.my_conti
    }

    /// `FaceTangency(E0, E1, V)` (`ChFi3d_Builder_1.cxx:564-623`): the faces
    /// opposite to `E0` and `E1` are tangent through the faces starting at `V`.
    ///
    /// `Err` carries the `Standard_ConstructionError("ChFi3d_Builder:only 2
    /// faces")` OCCT raises when an edge has more than two faces, and the
    /// failure `ChFi3d::IsTangentFaces` may raise.
    pub fn face_tangency(&self, e0: &Edge, e1: &Edge, v: &Vertex) -> Result<bool, String> {
        // It is checked if the connection is not on a regular edge.
        let mut faces: Vec<Face> = Vec::with_capacity(2);
        for value in self.my_ef_map.get_from_key(&e1.0) {
            if faces.len() > 1 {
                return Err("ChFi3d_Builder:only 2 faces".to_string());
            }
            let face = Face::wrap(value.clone())
                .ok_or_else(|| "ChFi3d_Builder:only 2 faces".to_string())?;
            faces.push(face);
        }
        if faces.len() < 2 {
            return Ok(false);
        }
        // `ChFi3d::IsTangentFaces(E1, F[0], F[1])` with its `GeomAbs_G1`
        // default (`ChFi3d.hxx:44-48`).
        if is_tangent_faces(e1, &faces[0], &faces[1], Shape::G1)? {
            return Ok(false);
        }

        for value in self.my_ve_map.get_from_key(&v.0) {
            let Some(ec) = Edge::wrap(value.clone()) else {
                continue;
            };
            if is_same(&ec.0, &e0.0) || is_same(&ec.0, &e1.0) {
                continue;
            }
            let orientation = ec.0.orientation();
            if orientation == Orientation::Internal || orientation == Orientation::External {
                continue;
            }
            if BRepTool::is_degenerated(&ec) {
                continue;
            }
            let mut faces: Vec<Face> = Vec::with_capacity(2);
            for candidate in self.my_ef_map.get_from_key(&ec.0) {
                if faces.len() > 1 {
                    return Err("ChFi3d_Builder:only 2 faces".to_string());
                }
                let face = Face::wrap(candidate.clone())
                    .ok_or_else(|| "ChFi3d_Builder:only 2 faces".to_string())?;
                faces.push(face);
            }
            if faces.len() < 2 {
                return Ok(false);
            }
            if !is_tangent_faces(&ec, &faces[0], &faces[1], Shape::G1)? {
                return Ok(false);
            }
        }
        Ok(true)
    }}
