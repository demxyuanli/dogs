//! `LocOpe_Pipe` (`LocOpe_Pipe.hxx:1-70`, `.cxx:1-550`, `.lxx:1-45`): the
//! `BRepFill_Pipe` sweep with the modifications the Pipe feature needs.
//!
//! Module table:
//!
//! | OCCT member | Rust location | Reference |
//! |---|---|---|
//! | constructor, i.e. the face-fusion pass | `constructor.rs` | `cxx:51-273` |
//! | `Shape` / `Shapes` / `Spine` / `Profile` / `FirstShape` / `LastShape` | `mod.rs` | `cxx:279-324`, `.lxx:19-45` |
//! | `Curves` / `BarycCurve` (+ the file-local static `Orientation`) | `curves.rs` | `cxx:329-548` |
//!
//! The class declares `myFirstShape` / `myLastShape` (`hxx:67-68`) but never
//! assigns them: `FirstShape()` / `LastShape()` answer `myPipe.FirstShape()` /
//! `myPipe.LastShape()` (`.lxx:34-45`). The port carries no fields for them.
//!
//! Reference behaviours with no directly translatable counterpart, each
//! documented at its call site:
//!
//! * `BarycCurve` (`cxx:494-499`) appends the null handle and then reads
//!   `seq(1)` of the empty sequence, which raises `Standard_OutOfRange`; the
//!   port answers that error.
//! * `Curves` (`cxx:357`) assigns `P1 = C->Value(p2)` and never reads it again;
//!   the port keeps the evaluation and marks it dead.
//! * `Geom_BSplineCurve::NbKnots()` is the *distinct* knot count
//!   (`Geom_BSplineCurve_1.cxx:380-383`), so the merge reads the distinct
//!   knots, not the flat knot vector.

mod constructor;
mod curves;

use std::sync::Arc;

use occt_core::gp::GpPnt;
use occt_geom::geom_fill::GeomFillTrihedron;
use occt_geom::Curve;

use crate::abs::ShapeType;
use crate::brep_fill_pipe::Pipe as BRepFillPipe;
use crate::locope::shape_maps::ShapeValueMap;
use crate::locope::sample_edges;
use crate::shape::{Edge, TopoShape, Vertex, Wire};
use crate::topo_tools::{explorer_shapes, is_same};

/// `LocOpe_Pipe` (`LocOpe_Pipe.hxx:36-69`): the class data.
pub struct Pipe {
    /// `myPipe` (`hxx:62`).
    my_pipe: BRepFillPipe,
    /// `myMap` (`hxx:63`): profile edge -> the faces swept from it.
    my_map: ShapeValueMap<Vec<TopoShape>>,
    /// `myRes` (`hxx:64`).
    my_res: Option<TopoShape>,
    /// `myGShap` (`hxx:65`).
    my_gshap: Vec<TopoShape>,
    /// `myCrvs` (`hxx:66`).
    my_crvs: Vec<Option<Arc<dyn Curve>>>,
}

impl Pipe {
    /// `LocOpe_Pipe(Spine, Profile)` (`cxx:51-273`).
    ///
    /// `BRepFill_Pipe(Spine, Profile)` keeps the defaults of
    /// `BRepFill_Pipe.hxx:76`, i.e. `GeomFill_IsCorrectedFrenet`, no forced C1
    /// and `GeneratePartCase = false` (`BRepFill_Pipe.cxx:161-189`).
    pub fn new(spine: &Wire, profile: &TopoShape) -> Result<Self, String> {
        let my_pipe = BRepFillPipe::with_spine(
            spine,
            profile,
            GeomFillTrihedron::CorrectedFrenet,
            false,
            false,
        )?;
        let (my_map, my_res) = constructor::fuse_faces(&my_pipe, spine, profile)?;
        Ok(Self {
            my_pipe,
            my_map,
            my_res,
            my_gshap: Vec::new(),
            my_crvs: Vec::new(),
        })
    }

    /// `Spine()` (`.lxx:19-22`).
    pub fn spine(&self) -> Option<&TopoShape> {
        self.my_pipe.spine().map(|w| &w.0)
    }

    /// `Profile()` (`.lxx:26-29`).
    pub fn profile(&self) -> Option<&TopoShape> {
        self.my_pipe.profile()
    }

    /// `FirstShape()` (`.lxx:33-36`).
    pub fn first_shape(&self) -> Option<&TopoShape> {
        self.my_pipe.first_shape()
    }

    /// `LastShape()` (`.lxx:40-43`).
    pub fn last_shape(&self) -> Option<&TopoShape> {
        self.my_pipe.last_shape()
    }

    /// `Shape()` (`cxx:279-282`).
    pub fn shape(&self) -> Option<&TopoShape> {
        self.my_res.as_ref()
    }

    /// `Shapes(S)` (`cxx:287-324`): the shapes swept from one profile edge
    /// (the faces stored by the constructor) or from one profile vertex (the
    /// rails of every spine edge).
    pub fn shapes(&mut self, s: &TopoShape) -> Result<&[TopoShape], String> {
        // `cxx:288-292`.
        let typ_s = s.shape_type();
        if typ_s != ShapeType::Edge && typ_s != ShapeType::Vertex {
            return Err("Standard_DomainError: LocOpe_Pipe::Shapes".to_string());
        }
        // `cxx:293-303`: "the shape has to belong to the profile".
        let profile = self
            .my_pipe
            .profile()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone: LocOpe_Pipe::Shapes".to_string())?;
        if !explorer_shapes(&profile, typ_s)
            .iter()
            .any(|cur| is_same(cur, s))
        {
            return Err("Standard_NoSuchObject: LocOpe_Pipe::Shapes".to_string());
        }

        // `cxx:306`.
        self.my_gshap.clear();

        // `cxx:307-320`: "a vertex gives the rails".
        if typ_s == ShapeType::Vertex {
            let v_profile = Vertex::wrap(s.clone())
                .ok_or_else(|| "Standard_TypeMismatch: LocOpe_Pipe::Shapes".to_string())?;
            let spine = self
                .my_pipe
                .spine()
                .cloned()
                .ok_or_else(|| "StdFail_NotDone: LocOpe_Pipe::Shapes".to_string())?;
            for edsp in explorer_shapes(&spine.0, ShapeType::Edge) {
                let edsp = Edge::wrap(edsp)
                    .ok_or_else(|| "Standard_TypeMismatch: LocOpe_Pipe::Shapes".to_string())?;
                // `cxx:313-317`. `BRepFill_Pipe::Edge` raises instead of
                // answering a null edge when the rail slot is missing
                // (`BRepFill_Pipe.cxx:450-451`), so the reference's
                // `!resed.IsNull()` guard never filters here.
                let resed = self.my_pipe.edge(&edsp, &v_profile)?;
                self.my_gshap.push(resed.0);
            }
            return Ok(self.my_gshap.as_slice());
        }

        // `cxx:323-324`: "an edge gives the faces stored by the constructor".
        let e_profile = Edge::wrap(s.clone())
            .ok_or_else(|| "Standard_TypeMismatch: LocOpe_Pipe::Shapes".to_string())?;
        match self.my_map.find(&e_profile.0) {
            Some(list) => Ok(list.as_slice()),
            None => Err("Standard_NoSuchObject: LocOpe_Pipe::Shapes".to_string()),
        }
    }

    /// `Curves(Spt)` (`cxx:329-426`): one curve per point of `Spt`, the sweep
    /// of that point along the spine.
    ///
    /// `myCrvs.Clear()` first (`cxx:333`), as the reference does; a point whose
    /// swept wire carries no edge stores the null handle (`cxx:370-374`).
    pub fn curves(&mut self, spt: &[GpPnt]) -> Result<&[Option<Arc<dyn Curve>>], String> {
        self.my_crvs.clear();
        for point in spt {
            let (mut seq, max_deg) = curves::collect_bsplines(self, point)?;
            if seq.is_empty() {
                self.my_crvs.push(None);
                continue;
            }
            let merged = curves::merge_bsplines(&mut seq, max_deg)?;
            self.my_crvs.push(Some(Arc::new(merged)));
        }
        Ok(self.my_crvs.as_slice())
    }

    /// `BarycCurve()` (`cxx:449-549`): the sweep of the barycentre of
    /// `FirstShape`'s samples.
    pub fn baryc_curve(&mut self) -> Result<Option<Arc<dyn Curve>>, String> {
        // `cxx:452-462`: the barycentre of `LocOpe::SampleEdges(FirstShape)`.
        let base = self
            .first_shape()
            .cloned()
            .ok_or_else(|| "StdFail_NotDone: LocOpe_Pipe::BarycCurve".to_string())?;
        let spt = sample_edges(&base);
        let mut sx = 0.0;
        let mut sy = 0.0;
        let mut sz = 0.0;
        for pvt in &spt {
            sx += pvt.x();
            sy += pvt.y();
            sz += pvt.z();
        }
        let n = spt.len() as f64;
        // `bar.ChangeCoord().Divide(spt.Length())`; an empty sample set divides
        // by zero exactly as the reference does.
        let p1 = GpPnt::new(sx / n, sy / n, sz / n);

        // `cxx:467-491`.
        let (mut seq, max_deg) = curves::collect_bsplines(self, &p1)?;
        // `cxx:494-499`: the reference appends the (null) handle to `myCrvs`
        // and then reads `seq(1)`, which raises `Standard_OutOfRange` when the
        // swept wire carried no edge.
        if seq.is_empty() {
            self.my_crvs.push(None);
            return Err(
                "Standard_OutOfRange: LocOpe_Pipe::BarycCurve reads seq(1) of an empty sequence"
                    .to_string(),
            );
        }
        // `cxx:500-548`.
        let merged = curves::merge_bsplines(&mut seq, max_deg)?;
        Ok(Some(Arc::new(merged)))
    }

    /// `myPipe` (`hxx:62`), mutable: `PipeLine` updates the sweep state
    /// (`BRepFill_Pipe.cxx:502-526`).
    pub(super) fn pipe_mut(&mut self) -> &mut BRepFillPipe {
        &mut self.my_pipe
    }
}
