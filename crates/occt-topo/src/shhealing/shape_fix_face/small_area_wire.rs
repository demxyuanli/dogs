//! `ShapeFix_Face::FixSmallAreaWire` (`ShapeFix_Face.cxx:2331-2394`), the tail
//! step of `ShapeFix_Face::Perform` (`cxx:730-736`).
//!
//! The step runs under `NeedFix(myFixSmallAreaWireMode, false)`, i.e. only when
//! the mode is explicitly `1`; both reader paths leave it at `ClearModes`' `-1`
//! (`ShapeProcess_OperLibrary.cxx:831` copies `FixSmallAreaWireMode`, whose
//! default is `FixMode::FixOrNot`), so it is inert on `FromSTEP.FixShape` /
//! `FromIGES.FixShape` and only reachable when a caller turns it on.
//!
//! UNPORTED: `Context()->Remove(myFace)` for the all-wires-removed case with
//! `theIsRemoveSmallFace` (`cxx:2381-2384`). `ShapeBuild_ReShape::Remove` has no
//! counterpart in the `ReShape` trait (which only carries `apply` / `replace`),
//! and `theIsRemoveSmallFace` is `NeedFix(myRemoveSmallAreaFaceMode, false)` =
//! false on both reader paths, so that arm is unreachable there.

use super::prelude::*;
use super::ShapeFixFace;

impl ShapeFixFace {
    /// `ShapeFix_Face::FixSmallAreaWire(theIsRemoveSmallFace)`
    /// (`ShapeFix_Face.cxx:2331-2394`).
    ///
    /// Walks `myFace`'s stored wires, drops every one whose
    /// `ShapeAnalysis_Wire::CheckSmallArea` reports a null area
    /// (`cxx:2360-2366`) and rebuilds the face from the survivors. Returns
    /// false when nothing was dropped or when every wire was.
    pub fn fix_small_area_wire(&mut self, the_is_remove_small_face: bool) -> bool {
        // cxx:2333-2337: run on `Context()->Apply(myFace)`.
        let mut face = match self.face.clone() {
            Some(f) => f,
            None => return false,
        };
        let mut ctx = self.context.clone();
        let applied = ctx.apply(&face.0);
        if applied.is_face() {
            face = Face(applied);
        }

        // `TopoDS_Iterator aWIt(myFace, false)` (`cxx:2345`): the stored wires,
        // not the ones composed through `myFace`'s orientation. OCCT's skip
        // test (`cxx:2349-2354`) only rejects a child that is neither a wire nor
        // FORWARD/REVERSED, so a well-formed face's extra children are the only
        // things it can drop; the port keeps the wire filter alone.
        let stored: Vec<Wire> = {
            let ts = face.0.tshape.read().expect("poisoned TShape lock");
            ts.children
                .iter()
                .filter(|s| s.shape_type() == ShapeType::Wire)
                .cloned()
                .map(Wire)
                .collect()
        };

        // cxx:2359-2369: `ShapeAnalysis_Wire(aWire, myFace, Precision())` then
        // `CheckSmallArea(aWire)`.
        let mut nb_removed = 0usize;
        let mut nb_wires = 0usize;
        let mut kept: Vec<Wire> = Vec::with_capacity(stored.len());
        for wire in stored {
            if crate::shhealing::check_small_area(&wire, &face) {
                nb_removed += 1; // cxx:2362-2365: null area, wire skipped
            } else {
                kept.push(wire);
                nb_wires += 1;
            }
        }

        if nb_removed == 0 {
            return false; // cxx:2371-2374
        }
        if nb_wires == 0 {
            if the_is_remove_small_face {
                // cxx:2381-2384 `Context()->Remove(myFace)` — UNPORTED, see the
                // module header.
            }
            return false;
        }

        // cxx:2386-2392: `aFace.Orientation(myFace.Orientation())`, then
        // `Context()->Replace(myFace, aFace)`, then `myFace = aFace`.
        //
        // The port rewrites the face's child list in place instead of
        // `EmptyCopied` + `Add`, which is the same substitution seen through
        // `Context()->Apply`: the face's own `BRep_TFace` (surface, tolerance)
        // and the surviving edges' pcurves stay bound to the shape, and keeping
        // `myFace`'s orientation leaves the `Replace` record an identity.
        {
            let mut ts = face.0.tshape.write().expect("poisoned TShape lock");
            ts.children = kept.into_iter().map(|w| w.0).collect();
        }
        self.face = Some(face);
        true
    }
}
