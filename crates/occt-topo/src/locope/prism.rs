//! `LocOpe_Prism` (`LocOpe_Prism.hxx:25-90`, `.cxx:48-295`).
//!
//! `BRepSweep_Prism` plus the feature-side bookkeeping `LocOpe_Generator` needs:
//! `FirstShape` / `LastShape`, and the descendant list of every edge of the
//! base (`Shapes`).
//!
//! Module table (`ModelingAlgorithms/TKFeat/LocOpe`):
//!
//! | OCCT | Rust |
//! |---|---|
//! | ctors (`cxx:48-71`) / `Perform` (`cxx:75-103`) | this file |
//! | `IntPerf` (`cxx:107-226`) | [`Prism::int_perf`] |
//! | `Shape` / `FirstShape` / `LastShape` / `Shapes` (`cxx:232-259`) | accessors |
//! | `Curves` (`cxx:264-279`) / `BarycCurve` (`cxx:284-295`) | [`Prism::curves`] / [`Prism::baryc_curve`] |
//!
//! Container mapping: `myMap` (`NCollection_DataMap<Shape,
//! NCollection_List<Shape>>`) is [`ShapeValueMap`]; the reference's `myMap(S)`
//! raises `Standard_NoSuchObject` on an unbound key, so [`Prism::shapes`]
//! answers `None` there.
//!
//! Deviations, all forced by the language or already recorded for this tree:
//! - `myRes` has no null state, so it is `Option<TopoShape>` and `Shape()` is
//!   `None` where OCCT raises `StdFail_NotDone`.
//! - `Curves` builds `Geom_TrimmedCurve` through [`GeomTrimmedCurve::new`],
//!   which reparameterises to `[0, 1]` (the tree-wide convention recorded in
//!   T-493); OCCT keeps the line's own `[-2h, 2h]` range.

use std::sync::Arc;

use occt_core::gp::{GpDir, GpLin, GpPnt, GpTrsf, GpVec};
use occt_geom::curve::Curve;
use occt_geom::line::GeomLine;
use occt_geom::trimmed::GeomTrimmedCurve;

use crate::abs::ShapeType;
use crate::brep_sweep::Prism as SweepPrism;
use crate::brep_tools_modifier::{Modifier, TrsfModification};
use crate::locope::build_shape::BuildShape;
use crate::locope::sample_edges::sample_edges;
use crate::locope::shape_maps::ShapeValueMap;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{explorer_shapes, AncestorMap};

/// `LocOpe_Prism` (`hxx:38-90`).
#[derive(Debug, Clone, Default)]
pub struct Prism {
    /// `TopoDS_Shape myBase`
    my_base: Option<TopoShape>,
    /// `gp_Vec myVec`
    my_vec: GpVec,
    /// `gp_Vec myTra`
    my_tra: GpVec,
    /// `bool myIsTrans`
    my_is_trans: bool,
    /// `bool myDone`
    my_done: bool,
    /// `TopoDS_Shape myRes`
    my_res: Option<TopoShape>,
    /// `TopoDS_Shape myFirstShape`
    my_first_shape: Option<TopoShape>,
    /// `TopoDS_Shape myLastShape`
    my_last_shape: Option<TopoShape>,
    /// `NCollection_DataMap<Shape, List<Shape>> myMap`
    my_map: ShapeValueMap<Vec<TopoShape>>,
}

impl Prism {
    /// `LocOpe_Prism()` (`cxx:48-53`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `LocOpe_Prism(Base, V)` (`cxx:57-60`).
    pub fn with_base(base: &TopoShape, v: &GpVec) -> Self {
        let mut this = Self::new();
        this.my_base = Some(base.clone());
        this.my_vec = *v;
        this.my_is_trans = false;
        this.int_perf();
        this
    }

    /// `LocOpe_Prism(Base, V, Vtra)` (`cxx:65-71`).
    pub fn with_base_and_translation(base: &TopoShape, v: &GpVec, vtra: &GpVec) -> Self {
        let mut this = Self::new();
        this.my_base = Some(base.clone());
        this.my_vec = *v;
        this.my_tra = *vtra;
        this.my_is_trans = true;
        this.int_perf();
        this
    }

    /// `Perform(Base, V)` (`cxx:75-87`).
    pub fn perform(&mut self, base: &TopoShape, v: &GpVec) {
        self.reset();
        self.my_base = Some(base.clone());
        self.my_vec = *v;
        self.my_is_trans = false;
        self.int_perf();
    }

    /// `Perform(Base, V, Vtra)` (`cxx:91-103`).
    pub fn perform_with_translation(&mut self, base: &TopoShape, v: &GpVec, vtra: &GpVec) {
        self.reset();
        self.my_base = Some(base.clone());
        self.my_vec = *v;
        self.my_tra = *vtra;
        self.my_is_trans = true;
        self.int_perf();
    }

    /// The `Clear` / `Nullify` prelude shared by both `Perform` overloads
    /// (`cxx:77-83`, `cxx:93-99`).
    fn reset(&mut self) {
        self.my_map.clear();
        self.my_first_shape = None;
        self.my_last_shape = None;
        self.my_base = None;
        self.my_res = None;
    }

    /// `IntPerf()` (`cxx:107-226`).
    fn int_perf(&mut self) {
        let Some(mut the_base) = self.my_base.clone() else {
            return;
        };

        // `cxx:111-121`: the translated copy of the base, when the feature is
        // given a translation.
        let mut modif: Option<Modifier> = None;
        if self.my_is_trans {
            let mut t = GpTrsf::identity();
            t.set_translation_vec(&self.my_tra);
            let modbase = TrsfModification::new(t);
            let mut m = Modifier::new(false);
            m.init(&the_base);
            if m.perform(&modbase).is_ok() {
                if let Ok(modified) = m.modified_shape(&the_base) {
                    the_base = modified;
                }
            }
            modif = Some(m);
        }

        // `cxx:123`: `BRepSweep_Prism thePrism(theBase, myVec)` with the
        // reference's default `Copy = true` / `Canonize = true`
        // (`BRepSweep_Prism.hxx:59-68`).
        let Ok(the_prism) = SweepPrism::from_vec(the_base.clone(), &self.my_vec, true, true) else {
            // The reference raises `Standard_ConstructionError` from the sweep
            // constructor on a null vector; `myDone` stays false.
            return;
        };

        // `cxx:125-126`.
        self.my_first_shape = the_prism.first_shape();
        self.my_last_shape = the_prism.last_shape();

        if the_base.shape_type() == ShapeType::Face {
            // `cxx:128-146`.
            for e in explorer_shapes(&the_base, ShapeType::Edge) {
                let edg = Edge(e);
                if self.my_map.is_bound(&edg.0) {
                    continue;
                }
                self.my_map.bind(&edg.0, Vec::new());
                // `cxx:136-141`: `thePrism.Shape(edg)` unless null.
                if let Some(desc) = the_prism.shape_of(&edg.0) {
                    if let Some(v) = self.my_map.change(&edg.0) {
                        v.push(desc);
                    }
                }
            }
            // `cxx:145`.
            self.my_res = the_prism.shape();
        } else {
            // `cxx:149-208`: the base is not a single face.
            let the_ef_map =
                AncestorMap::map_shapes_and_ancestors(&the_base, ShapeType::Edge, ShapeType::Face);
            let mut lfaces: Vec<TopoShape> = Vec::new();
            let mut toremove = false;
            // `cxx:157-177`.
            for i in 0..the_ef_map.length() {
                let edg = the_ef_map.find_key(i).clone();
                self.my_map.bind(&edg, Vec::new());
                if let Some(desc) = the_prism.shape_of(&edg) {
                    if the_ef_map.ancestors(i).len() >= 2 {
                        toremove = true;
                    } else {
                        if let Some(v) = self.my_map.change(&edg) {
                            v.push(desc.clone());
                        }
                        lfaces.push(desc);
                    }
                }
            }

            if toremove {
                // `cxx:178-189`: the faces of `FirstShape` and `LastShape` go
                // back into the face list, then `LocOpe_BuildShape` rebuilds.
                if let Some(first) = self.my_first_shape.clone() {
                    for f in explorer_shapes(&first, ShapeType::Face) {
                        lfaces.push(f);
                    }
                }
                if let Some(last) = self.my_last_shape.clone() {
                    for f in explorer_shapes(&last, ShapeType::Face) {
                        lfaces.push(f);
                    }
                }
                let mut bs = BuildShape::new();
                if bs.perform(&lfaces).is_ok() {
                    self.my_res = bs.shape().cloned();
                }
            } else {
                // `cxx:191-207`.
                for e in explorer_shapes(&the_base, ShapeType::Edge) {
                    let edg = Edge(e);
                    if self.my_map.is_bound(&edg.0) {
                        continue;
                    }
                    self.my_map.bind(&edg.0, Vec::new());
                    if let Some(desc) = the_prism.shape_of(&edg.0) {
                        if let Some(v) = self.my_map.change(&edg.0) {
                            v.push(desc);
                        }
                    }
                }
                self.my_res = the_prism.shape();
            }
        }

        // `cxx:211-225`: when the base was translated, move the descendant
        // lists onto the *original* edges and drop the translated keys.
        if let Some(modif) = modif {
            let base = self.my_base.clone();
            if let Some(base) = base {
                for e in explorer_shapes(&base, ShapeType::Edge) {
                    let edg = Edge(e);
                    let Some(edgbis) = modif.modified_shape(&edg.0).ok().and_then(Edge::wrap) else {
                        continue;
                    };
                    if !crate::topo_tools::is_same(&edgbis.0, &edg.0) && self.my_map.is_bound(&edgbis.0) {
                        let moved = self.my_map.find(&edgbis.0).cloned().unwrap_or_default();
                        self.my_map.bind(&edg.0, moved);
                        self.my_map.unbind(&edgbis.0);
                    }
                }
            }
        }

        // `cxx:225`.
        self.my_done = true;
    }

    /// `Shape()` (`cxx:232-237`); `None` is the `StdFail_NotDone` raise.
    pub fn shape(&self) -> Option<&TopoShape> {
        if self.my_done {
            self.my_res.as_ref()
        } else {
            None
        }
    }

    /// `IsDone()`: the reference has no such method; it is the `myDone` flag
    /// behind the `StdFail_NotDone` raise, exposed so a delegating wrapper
    /// (`LocOpe_LinearForm`) can mirror the same raise.
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `FirstShape()` (`cxx:243-247`).
    pub fn first_shape(&self) -> Option<&TopoShape> {
        self.my_first_shape.as_ref()
    }

    /// `LastShape()` (`cxx:252-256`).
    pub fn last_shape(&self) -> Option<&TopoShape> {
        self.my_last_shape.as_ref()
    }

    /// `Shapes(S)` (`cxx:256-259`).
    pub fn shapes(&self, s: &TopoShape) -> Option<&[TopoShape]> {
        self.my_map.find(s).map(|v| v.as_slice())
    }

    /// `Curves(Scurves)` (`cxx:264-279`): one trimmed line per sampled point of
    /// `FirstShape`, along the prism vector, over `[-2h, 2h]` with
    /// `h = |V|`.
    pub fn curves(&self) -> Vec<Arc<dyn Curve>> {
        let mut scurves: Vec<Arc<dyn Curve>> = Vec::new();
        let Some(first) = self.my_first_shape.clone() else {
            return scurves;
        };
        let spt = sample_edges(&first);
        let v = self.my_vec;
        let height = (v.x() * v.x() + v.y() * v.y() + v.z() * v.z()).sqrt();
        let u1 = -2.0 * height;
        let u2 = 2.0 * height;
        for p in &spt {
            // `gp_Ax1 theAx(spt(jj), myVec)` then `Geom_Line(theAx)`.
            let ax = axis_line(p, &v);
            let Some(ax) = ax else {
                continue;
            };
            let the_lin: Arc<dyn Curve> = Arc::new(ax);
            // `Geom_TrimmedCurve(theLin, u1, u2, true)`.
            scurves.push(Arc::new(GeomTrimmedCurve::new(the_lin, u1, u2)));
        }
        scurves
    }

    /// `BarycCurve()` (`cxx:284-295`): the line through the barycentre of the
    /// sampled points of `FirstShape`, along the prism vector.
    pub fn baryc_curve(&self) -> Option<GeomLine> {
        let first = self.my_first_shape.clone()?;
        let spt = sample_edges(&first);
        if spt.is_empty() {
            // `bar.ChangeCoord().Divide(0)` is a division by zero in the
            // reference; answer `None` instead of producing NaN coordinates.
            return None;
        }
        let mut bar = GpPnt::new(0.0, 0.0, 0.0);
        for p in &spt {
            let c = bar.xyz();
            bar = GpPnt::new(c.x + p.x(), c.y + p.y(), c.z + p.z());
        }
        let n = spt.len() as f64;
        let c = bar.xyz();
        let bar = GpPnt::new(c.x / n, c.y / n, c.z / n);
        // `gp_Ax1 newAx(bar, myVec)` then `Geom_Line(newAx)`.
        axis_line(&bar, &self.my_vec)
    }
}

/// `Geom_Line(gp_Ax1(P, V))`: the `gp_Ax1` direction is the normalised `V`, so
/// a null vector is the reference's `Standard_ConstructionError`.
fn axis_line(p: &GpPnt, v: &GpVec) -> Option<GeomLine> {
    let d = GpDir::from_vec(v).ok()?;
    Some(GeomLine::new(GpLin::from_pnt_dir(*p, d)))
}
