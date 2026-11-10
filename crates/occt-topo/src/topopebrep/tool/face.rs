//! `TopOpeBRepTool_face` — the face a wire bounds on a reference face.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_face.cxx`
//! (OCCT 8.0.0, 129 lines): `FUN_reverse` (`:30-42`), `Init` (`:46-77`),
//! `IsDone` (`:82-85`), `Finite` (`:89-95`), `Ffinite` (`:99-105`), `W`
//! (`:110-113`) and `RealF` (`:122-129`).

use crate::abs::ShapeType;
use crate::builder::TopoBuilder;
use crate::fclass2d::{FClass2d, FaceState};
use crate::shape::{Face, TopoShape};
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::explorer_shapes;

/// `class TopOpeBRepTool_face` (`hxx:26-56`).
#[derive(Debug, Clone, Default)]
pub struct ToolFace {
    /// `myW` (`hxx:53`).
    my_w: Option<TopoShape>,
    /// `myfinite` (`hxx:54`).
    my_finite: bool,
    /// `myFfinite` (`hxx:55`); `None` is a null `TopoDS_Face`.
    my_f_finite: Option<TopoShape>,
}

/// `FUN_reverse(f, frev)` (`cxx:30-42`): empty copy of `f` with every
/// sub-shape added reversed.
///
/// `f.EmptyCopied()` on a face copies the surface, the tolerance and the
/// natural-restriction flag (`BRep_TFace.cxx:36-43`), so the copy keeps the
/// geometry of the source face.
fn fun_reverse(f: &TopoShape) -> TopoShape {
    let frev = f.empty_copied();
    let reg = GeometryRegistry::global();
    if let Some(geom) = reg.face_geom(f) {
        reg.set_face(&frev, geom);
    }
    let builder = TopoBuilder::new();
    let mut frev = frev;
    for child in explorer_shapes(f, ShapeType::Wire) {
        let mut reversed = child.clone();
        reversed.reverse();
        builder.add(&mut frev, &reversed);
    }
    frev
}

impl ToolFace {
    /// `TopOpeBRepTool_face()` (`hxx:42`, defaulted in `cxx:28`).
    pub fn new() -> Self {
        Self {
            my_w: None,
            my_finite: false,
            my_f_finite: None,
        }
    }

    /// `Init(W, Fref)` (`cxx:46-77`).
    ///
    /// `Fref.EmptyCopied()` carries the reference surface and tolerance; the
    /// wire is added to that copy and the copy is marked as a natural
    /// restriction before `BRepTopAdaptor_FClass2d` looks for the infinite
    /// point.
    pub fn init(&mut self, w: &TopoShape, f_ref: &TopoShape) -> bool {
        self.my_f_finite = None;
        self.my_w = Some(w.clone());

        // `cxx:55-58`: `TopoDS_Face fres = TopoDS::Face(Fref.EmptyCopied())`.
        let mut fres = f_ref.empty_copied();
        let reg = GeometryRegistry::global();
        if let Some(geom) = reg.face_geom(f_ref) {
            reg.set_face(&fres, geom);
        }

        // `cxx:58-59`: `B.Add(fres, W)`.
        TopoBuilder::new().add(&mut fres, w);

        // `cxx:60`: `B.NaturalRestriction(fres, true)`.
        reg.set_natural_restriction(&fres, true);

        // `cxx:63-65`: `BRepTopAdaptor_FClass2d FClass(fres, 0.)`;
        // `infinite = (FClass.PerformInfinitePoint() == TopAbs_IN)`.
        let infinite = match Face::wrap(fres.clone()) {
            Some(face) => match FClass2d::new(&face, 0.0) {
                Ok(class) => class.perform_infinite_point() == FaceState::In,
                // `FClass2d::new` is `BRepTopAdaptor_FClass2d`'s constructor,
                // which has no failing branch in OCCT; a rejected face here
                // leaves `PerformInfinitePoint` without a region to test.
                Err(_) => false,
            },
            None => false,
        };
        self.my_finite = !infinite;

        // `cxx:68-76`: `myFfinite = myfinite ? fres : FUN_reverse(fres, myFfinite)`.
        self.my_f_finite = Some(if self.my_finite {
            fres
        } else {
            fun_reverse(&fres)
        });
        true
    }

    /// `IsDone()` (`cxx:82-85`).
    pub fn is_done(&self) -> bool {
        self.my_f_finite.is_some()
    }

    /// `Finite()` (`cxx:89-95`).
    ///
    /// OCCT raises `Standard_Failure("TopOpeBRepTool_face NOT DONE")`; the
    /// port's callers check [`Self::is_done`] first, so the flag is returned
    /// as stored here.
    pub fn finite(&self) -> bool {
        self.my_finite
    }

    /// `Ffinite()` (`cxx:99-105`).
    pub fn f_finite(&self) -> Option<&TopoShape> {
        self.my_f_finite.as_ref()
    }

    /// `W()` (`cxx:110-113`).
    pub fn w(&self) -> Option<&TopoShape> {
        self.my_w.as_ref()
    }

    /// `RealF()` (`cxx:122-129`).
    pub fn real_f(&self) -> Option<TopoShape> {
        let f_finite = self.my_f_finite.as_ref()?;
        if self.my_finite {
            return Some(f_finite.clone());
        }
        Some(fun_reverse(f_finite))
    }
}
