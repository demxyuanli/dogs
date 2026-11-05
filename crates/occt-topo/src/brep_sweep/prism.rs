//! `BRepSweep_Prism` (`BRepSweep_Prism.cxx:32-145`): the translated prism,
//! a thin front end over [`Translation`] (`BRepSweep_Translation`).
//!
//! The three private helpers of the reference are reproduced here because they
//! decide the directing indexation and the location of the sweep:
//! - `NumShape()` (`cxx:106-109`) = `Sweep_NumShape(2, TopAbs_EDGE)`, i.e. a
//!   bounded directing edge with two vertices; `closed`, `BegInf`, `EndInf`
//!   all take the `Sweep_NumShape` defaults (`Sweep_NumShape.hxx:49-54`).
//! - `NumShape(Inf)` (`cxx:111-122`): infinite → `Init(0, EDGE, false, true,
//!   true)`; semi-infinite along `D` → `Init(1, EDGE, false, false, true)`.
//! - `Location(V)` (`cxx:126-132`): `gp_Trsf::SetTranslation(V)`.
//!
//! The `gp_Dir` constructor (`cxx:43-49`) passes `gp_Vec(D)` in both the
//! `Location` and the `V` argument, i.e. a **unit** vector
//! (`gp_Dir::operator gp_Vec`, `gp_Dir.hxx`), which is why `Vec()` answers a
//! unitar vector for an infinite prism (`BRepSweep_Prism.hxx:81-83`).

use occt_core::gp::{GpDir, GpTrsf, GpVec};
use occt_core::precision::CONFUSION;
use occt_core::toploc::TopLocLocation;

use crate::abs::ShapeType;
use crate::shape::TopoShape;

use super::grid::SweepOps;
use super::num_shape::{NumShape, NumShapeTool};
use super::translation::Translation;

/// `BRepSweep_Prism` (`BRepSweep_Prism.hxx:36-105`).
pub struct Prism {
    /// `myTranslation`.
    pub translation: Translation,
}

impl Prism {
    /// `BRepSweep_Prism(S, V, Copy, Canonize)` (`cxx:32-39`).
    ///
    /// The reference constructs `myTranslation` first and only then raises on
    /// `V.Magnitude() <= Precision::Confusion()`; the base constructor raises
    /// on the weaker `V.Magnitude() < Precision::Confusion()`
    /// (`BRepSweep_Translation.cxx:130-131`), so the prism's stricter test is
    /// the one that decides and the port answers it up front.
    pub fn from_vec(
        shape: TopoShape,
        v: &GpVec,
        copy: bool,
        canonize: bool,
    ) -> Result<Self, String> {
        if v.magnitude() <= CONFUSION {
            return Err("BRepSweep_Prism::Constructor".into());
        }
        // `Location(V)` (`cxx:126-132`).
        let mut gpt = GpTrsf::identity();
        gpt.set_translation_vec(v);
        let location = TopLocLocation::composed(&gpt, TopLocLocation::identity());
        let translation = Translation::new(shape, Self::num_shape(), location, *v, copy, canonize)?;
        Ok(Self { translation })
    }

    /// `BRepSweep_Prism(S, D, Inf, Copy, Canonize)` (`cxx:43-49`).
    ///
    /// `D` converts to the unit vector `gp_Vec(D)`, used both as the location
    /// translation and as the sweep vector. The base constructor's
    /// `V.Magnitude() < Precision::Confusion()` can never fire on a unit
    /// vector, and the prism's own constructor has no check, so this path
    /// always builds.
    pub fn from_dir(
        shape: TopoShape,
        d: &GpDir,
        inf: bool,
        copy: bool,
        canonize: bool,
    ) -> Result<Self, String> {
        // `gp_Vec(D)` (`gp_Vec.hxx:64`): the coordinate triple of the
        // direction, i.e. a unit vector.
        let v = GpVec::new(d.x(), d.y(), d.z());
        let mut gpt = GpTrsf::identity();
        gpt.set_translation_vec(&v);
        let location = TopLocLocation::composed(&gpt, TopLocLocation::identity());
        let translation =
            Translation::new(shape, Self::num_shape_inf(inf), location, v, copy, canonize)?;
        Ok(Self { translation })
    }

    /// `NumShape()` (`cxx:106-109`): `Sweep_NumShape(2, TopAbs_EDGE)` with the
    /// `Sweep_NumShape` defaults for `Closed` / `BegInf` / `EndInf`.
    fn num_shape() -> NumShape {
        let mut n = NumShape::dummy();
        n.init(2, ShapeType::Edge, false, false, false);
        n
    }

    /// `NumShape(const bool Inf)` (`cxx:111-122`).
    fn num_shape_inf(inf: bool) -> NumShape {
        let mut n = NumShape::dummy();
        if inf {
            n.init(0, ShapeType::Edge, false, true, true);
        } else {
            n.init(1, ShapeType::Edge, false, false, true);
        }
        n
    }

    /// `Shape()` (`cxx:54-57`).
    pub fn shape(&self) -> Option<TopoShape> {
        self.translation.shape()
    }

    /// `Shape(aGenS)` (`cxx:61-64`): the reference answers a null shape unless
    /// `myGenShapeTool.Index(aGenS) != 0 && HasShape(aGenS, myDirWire)`
    /// (`BRepSweep_NumLinearRegularSweep.cxx`), which is `None` here.
    pub fn shape_of(&self, gen: &TopoShape) -> Option<TopoShape> {
        let tables = &self.translation.tables;
        if tables.gen_tool.index(gen) == 0 || !self.translation.has_shape(gen, tables.dir) {
            return None;
        }
        super::grid::shape_of(&self.translation, gen, tables.dir)
    }

    /// `FirstShape()` (`cxx:66-69`).
    pub fn first_shape(&self) -> Option<TopoShape> {
        let gen = self.translation.tables.gen.clone();
        self.first_shape_of(&gen)
    }

    /// `FirstShape(aGenS)` (`cxx:71-74`).
    pub fn first_shape_of(&self, gen: &TopoShape) -> Option<TopoShape> {
        self.directed_shape(gen, true)
    }

    /// `LastShape()` (`cxx:81-84`).
    pub fn last_shape(&self) -> Option<TopoShape> {
        let gen = self.translation.tables.gen.clone();
        self.last_shape_of(&gen)
    }

    /// `LastShape(aGenS)` (`cxx:86-89`).
    pub fn last_shape_of(&self, gen: &TopoShape) -> Option<TopoShape> {
        self.directed_shape(gen, false)
    }

    /// The shared body of `FirstShape` / `LastShape`
    /// (`BRepSweep_NumLinearRegularSweep.cxx`): `HasFirstVertex()` /
    /// `HasLastVertex()` gate, then `HasShape` + `Shape(gen, vertex)`.
    fn directed_shape(&self, gen: &TopoShape, first: bool) -> Option<TopoShape> {
        let tables = &self.translation.tables;
        let tool = NumShapeTool::new(tables.dir);
        let (has, dir) = if first {
            (tool.has_first_vertex(), tool.first_vertex())
        } else {
            (tool.has_last_vertex(), tool.last_vertex())
        };
        if !has || !self.translation.has_shape(gen, dir) {
            return None;
        }
        super::grid::shape_of(&self.translation, gen, dir)
    }

    /// `Vec()` (`cxx:93-96`).
    pub fn vec(&self) -> GpVec {
        self.translation.vec
    }

    /// `IsUsed(aGenS)` (`cxx:135-138` →
    /// `BRepSweep_NumLinearRegularSweep::IsUsed`).
    pub fn is_used(&self, gen: &TopoShape) -> bool {
        let tables = &self.translation.tables;
        let i_gen = tables.gen_tool.index(gen);
        if i_gen == 0 {
            return false;
        }
        let upper_col = tables.dir_tool.nb_shapes() as usize;
        let mut is_built = false;
        let mut is_used = false;
        for j in 2..=upper_col {
            is_built = is_built || tables.built.borrow()[i_gen][j];
            is_used = is_used || tables.used.borrow()[i_gen][j];
        }
        let closed = self.translation.tables.dir.closed;
        if is_used {
            if gen.shape_type() == ShapeType::Vertex && self.translation.is_invariant(gen) {
                return tables.used.borrow()[i_gen][1] || !closed;
            }
            return is_used;
        }
        if is_built {
            if !self.translation.has_shape(gen, tables.dir) && !closed {
                return true;
            }
            if gen.shape_type() == ShapeType::Vertex && !closed && !tables.built.borrow()[i_gen][1] {
                return true;
            }
        }
        false
    }

    /// `GenIsUsed(aGenS)` (`cxx:141-144` →
    /// `BRepSweep_NumLinearRegularSweep::GenIsUsed`).
    pub fn gen_is_used(&self, gen: &TopoShape) -> bool {
        let tables = &self.translation.tables;
        let i_gen = tables.gen_tool.index(gen);
        if i_gen == 0 {
            return false;
        }
        if i_gen == 1 {
            tables.built.borrow()[i_gen][1]
        } else {
            tables.built.borrow()[i_gen][1] && tables.used.borrow()[i_gen][1]
        }
    }
}
