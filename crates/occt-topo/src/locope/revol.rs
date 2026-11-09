//! `LocOpe_Revol` (`LocOpe_Revol.hxx:25-95`, `.cxx:46-312`).
//!
//! `BRepSweep_Revol` plus the same feature-side bookkeeping as
//! [`super::prism::Prism`]: `FirstShape` / `LastShape` and the descendant list
//! of every edge of the base.
//!
//! Module table (`ModelingAlgorithms/TKFeat/LocOpe`):
//!
//! | OCCT | Rust |
//! |---|---|
//! | `Perform` (`cxx:57-90`) | this file |
//! | `IntPerf` (`cxx:94-212`) | [`Revol::int_perf`] |
//! | `Shape` / `FirstShape` / `LastShape` / `Shapes` (`cxx:216-245`) | accessors |
//! | `Curves` (`cxx:249-267`) / `BarycCurve` (`cxx:271-290`) | [`Revol::curves`] / [`Revol::baryc_curve`] |
//! | `static FindCircle` (`cxx:294-312`) | [`find_circle`] |
//!
//! Container mapping: `myMap` is [`ShapeValueMap`]; `myMap(S)` raises in the
//! reference on an unbound key, so [`Revol::shapes`] answers `None`.
//!
//! Deviations, all forced by the language:
//! - `myRes` has no null state, so it is `Option<TopoShape>` and `Shape()` is
//!   `None` where OCCT raises `StdFail_NotDone`.
//! - `Curves` / `BarycCurve` answer `GeomCircle` values instead of a
//!   `Geom_Curve` handle, because a circle is the only curve `FindCircle` can
//!   produce.

use std::sync::Arc;

use occt_core::gp::{GpAx1, GpAx2, GpCirc, GpDir, GpPnt, GpTrsf, GpVec};
use occt_core::precision::CONFUSION;
use occt_geom::circle::GeomCircle;
use occt_geom::curve::Curve;

use crate::abs::ShapeType;
use crate::brep_sweep::Revol as SweepRevol;
use crate::brep_tools_modifier::{Modifier, TrsfModification};
use crate::locope::build_shape::BuildShape;
use crate::locope::sample_edges::sample_edges;
use crate::locope::shape_maps::ShapeValueMap;
use crate::shape::{Edge, TopoShape};
use crate::topo_tools::{explorer_shapes, AncestorMap};

/// `LocOpe_Revol` (`hxx:38-95`).
#[derive(Debug, Clone, Default)]
pub struct Revol {
    /// `TopoDS_Shape myBase`
    my_base: Option<TopoShape>,
    /// `gp_Ax1 myAxis`
    my_axis: GpAx1,
    /// `double myAngle`
    my_angle: f64,
    /// `double myAngTra`
    my_ang_tra: f64,
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

impl Revol {
    /// `LocOpe_Revol()` (`cxx:46-54`).
    pub fn new() -> Self {
        Self {
            my_axis: GpAx1::default(),
            ..Default::default()
        }
    }

    /// `LocOpe_Revol(Base, Axis, Angle)` (`cxx` ctor delegating to `Perform`).
    pub fn with_base(base: &TopoShape, axis: &GpAx1, angle: f64) -> Self {
        let mut this = Self::new();
        this.perform(base, axis, angle);
        this
    }

    /// `LocOpe_Revol(Base, Axis, Angle, angledec)`.
    pub fn with_base_and_dec(base: &TopoShape, axis: &GpAx1, angle: f64, angledec: f64) -> Self {
        let mut this = Self::new();
        this.perform_with_dec(base, axis, angle, angledec);
        this
    }

    /// `Perform(Base, Axis, Angle)` (`cxx:57-70`).
    pub fn perform(&mut self, base: &TopoShape, axis: &GpAx1, angle: f64) {
        self.my_map.clear();
        self.my_first_shape = None;
        self.my_last_shape = None;
        self.my_base = Some(base.clone());
        self.my_angle = angle;
        self.my_axis = *axis;
        self.my_ang_tra = 0.0;
        self.my_is_trans = false;
        self.my_res = None;
        self.int_perf();
    }

    /// `Perform(Base, Axis, Angle, angledec)` (`cxx:74-90`).
    pub fn perform_with_dec(&mut self, base: &TopoShape, axis: &GpAx1, angle: f64, angledec: f64) {
        self.my_map.clear();
        self.my_first_shape = None;
        self.my_last_shape = None;
        self.my_base = Some(base.clone());
        self.my_angle = angle;
        self.my_axis = *axis;
        self.my_ang_tra = angledec;
        self.my_is_trans = true;
        self.my_res = None;
        self.int_perf();
    }

    /// `IntPerf()` (`cxx:94-212`).
    fn int_perf(&mut self) {
        let Some(mut the_base) = self.my_base.clone() else {
            return;
        };

        // `cxx:96-106`: the rotated copy of the base, when the feature is given
        // a start angle.
        let mut modif: Option<Modifier> = None;
        if self.my_is_trans {
            let mut t = GpTrsf::identity();
            if t.set_rotation_ax1(&self.my_axis, self.my_ang_tra).is_ok() {
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
        }

        // `cxx:108`: `BRepSweep_Revol theRevol(theBase, myAxis, myAngle)` with
        // the reference's default `Copy = true` (`BRepSweep_Revol.hxx:46-49`).
        let Ok(the_revol) = SweepRevol::new_with_angle(&the_base, &self.my_axis, self.my_angle, true)
        else {
            // The reference raises `Standard_ConstructionError` when
            // `|Angle| <= Precision::Angular()`; `myDone` stays false.
            return;
        };

        // `cxx:110-111`.
        self.my_first_shape = the_revol.first_shape();
        self.my_last_shape = the_revol.last_shape();

        if the_base.shape_type() == ShapeType::Face {
            // `cxx:114-131`.
            for e in explorer_shapes(&the_base, ShapeType::Edge) {
                let edg = Edge(e);
                if self.my_map.is_bound(&edg.0) {
                    continue;
                }
                self.my_map.bind(&edg.0, Vec::new());
                if let Some(desc) = the_revol.shape_of(&edg.0) {
                    if let Some(v) = self.my_map.change(&edg.0) {
                        v.push(desc);
                    }
                }
            }
            self.my_res = the_revol.shape();
        } else {
            // `cxx:134-195`: the base is not a single face.
            let the_ef_map =
                AncestorMap::map_shapes_and_ancestors(&the_base, ShapeType::Edge, ShapeType::Face);
            let mut lfaces: Vec<TopoShape> = Vec::new();
            let mut toremove = false;
            // `cxx:143-161`.
            for i in 0..the_ef_map.length() {
                let edg = the_ef_map.find_key(i).clone();
                self.my_map.bind(&edg, Vec::new());
                if let Some(desc) = the_revol.shape_of(&edg) {
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
                // `cxx:163-175`.
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
                // `cxx:176-193`.
                for e in explorer_shapes(&the_base, ShapeType::Edge) {
                    let edg = Edge(e);
                    if self.my_map.is_bound(&edg.0) {
                        continue;
                    }
                    self.my_map.bind(&edg.0, Vec::new());
                    if let Some(desc) = the_revol.shape_of(&edg.0) {
                        if let Some(v) = self.my_map.change(&edg.0) {
                            v.push(desc);
                        }
                    }
                }
                self.my_res = the_revol.shape();
            }
        }

        // `cxx:197-211`: when the base was rotated, move the descendant lists
        // onto the original edges and drop the rotated keys.
        if let Some(modif) = modif {
            if let Some(base) = self.my_base.clone() {
                for e in explorer_shapes(&base, ShapeType::Edge) {
                    let edg = Edge(e);
                    let Some(edgbis) = modif.modified_shape(&edg.0).ok().and_then(Edge::wrap) else {
                        continue;
                    };
                    if !crate::topo_tools::is_same(&edgbis.0, &edg.0)
                        && self.my_map.is_bound(&edgbis.0)
                    {
                        let moved = self.my_map.find(&edgbis.0).cloned().unwrap_or_default();
                        self.my_map.bind(&edg.0, moved);
                        self.my_map.unbind(&edgbis.0);
                    }
                }
            }
        }

        // `cxx:212`.
        self.my_done = true;
    }

    /// `Shape()` (`cxx:216-224`); `None` is the `StdFail_NotDone` raise.
    pub fn shape(&self) -> Option<&TopoShape> {
        if self.my_done {
            self.my_res.as_ref()
        } else {
            None
        }
    }

    /// The `myDone` flag behind the `StdFail_NotDone` raise, exposed for
    /// delegating wrappers (`LocOpe_RevolutionForm`).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// `FirstShape()` (`cxx:229-233`).
    pub fn first_shape(&self) -> Option<&TopoShape> {
        self.my_first_shape.as_ref()
    }

    /// `LastShape()` (`cxx:238-242`).
    pub fn last_shape(&self) -> Option<&TopoShape> {
        self.my_last_shape.as_ref()
    }

    /// `Shapes(S)` (`cxx:243-245`).
    pub fn shapes(&self, s: &TopoShape) -> Option<&[TopoShape]> {
        self.my_map.find(s).map(|v| v.as_slice())
    }

    /// `Curves(Scurves)` (`cxx:249-267`): one circle per sampled point of
    /// `FirstShape` that `FindCircle` accepts.
    pub fn curves(&self) -> Vec<Arc<dyn Curve>> {
        let mut scurves: Vec<Arc<dyn Curve>> = Vec::new();
        let Some(first) = self.my_first_shape.clone() else {
            return scurves;
        };
        let spt = sample_edges(&first);
        for p in &spt {
            if let Some(circ) = find_circle(&self.my_axis, p) {
                // `Geom_Circle(CAX.Position(), CAX.Radius())`.
                let ci: Arc<dyn Curve> = Arc::new(GeomCircle::new(circ));
                scurves.push(ci);
            }
        }
        scurves
    }

    /// `BarycCurve()` (`cxx:271-290`): the circle through the barycentre of
    /// the sampled points of `FirstShape`, when `FindCircle` accepts it.
    pub fn baryc_curve(&self) -> Option<GeomCircle> {
        let first = self.my_first_shape.clone()?;
        let spt = sample_edges(&first);
        if spt.is_empty() {
            // `bar.ChangeCoord().Divide(0)` divides by zero in the reference;
            // `FindCircle` on that NaN point answers `false`, so `None` is the
            // reference's answer.
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
        find_circle(&self.my_axis, &bar).map(GeomCircle::new)
    }
}

/// `static bool FindCircle(Ax, Pt, Ci)` (`cxx:294-312`): the circle of axis
/// `Ax` through `Pt`, or `None` when `Pt` sits on the axis
/// (`Radius < Precision::Confusion()`).
fn find_circle(ax: &GpAx1, pt: &GpPnt) -> Option<GpCirc> {
    let dax = ax.direction();
    let loc = ax.location();
    // `gp_Vec OP(Ax.Location(), Pt)` then `prm = OP.Dot(Dax)`.
    let op = GpVec::new(
        pt.x() - loc.x(),
        pt.y() - loc.y(),
        pt.z() - loc.z(),
    );
    let dax_v = GpVec::new(dax.x(), dax.y(), dax.z());
    let prm = op.dot(&dax_v);
    // `gp_Pnt prj(Ax.Location().XYZ().Added(prm * Dax.XYZ()))`.
    let prj = GpPnt::new(
        loc.x() + prm * dax.x(),
        loc.y() + prm * dax.y(),
        loc.z() + prm * dax.z(),
    );
    // `gp_Vec axx(prj, Pt)` then `Radius = axx.Magnitude()`.
    let axx = GpVec::new(
        pt.x() - prj.x(),
        pt.y() - prj.y(),
        pt.z() - prj.z(),
    );
    let radius = axx.magnitude();
    if radius < CONFUSION {
        return None;
    }
    // `Ci.SetRadius(Radius); Ci.SetPosition(gp_Ax2(prj, Dax, axx))`.
    let x_dir = GpDir::from_vec(&axx).ok()?;
    let pos = GpAx2::new(prj, *dax, x_dir).ok()?;
    Some(GpCirc::new(pos, radius))
}
