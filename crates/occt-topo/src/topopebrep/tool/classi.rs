//! `TopOpeBRepTool_CLASSI` — classify wires of a face against each other.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepTool/TopOpeBRepTool_CLASSI.cxx`
//! (OCCT 8.0.0, 616 lines): `Init2d` (`:52-55`), `HasInit2d` (`:59-62`),
//! `Add2d` (`:66-96`), `GetBox2d` (`:100-113`), `ClassiBnd2d` (`:118-217`),
//! `FUN_thegreatest` (`:219-250`), `Classip2d` (`:254-467`), `Getface`
//! (`:470-479`), `FUN_addOwlw` (`:484-500`) and `Classilist` (`:502-615`).
//!
//! The class answers five codes: `SAME`, `DIFF`, `UNKNOWN`, `ONE_IN_TWO` and
//! `TWO_IN_ONE` (`cxx:36-40`). They are the values `FUN_tool_ClassifW` of
//! `TopOpeBRepTool_RegularizeW.cxx` switches on.

use std::collections::HashMap;

use occt_core::bnd::BndBox2d;
use occt_core::gp::GpPnt2d;

use crate::abs::ShapeType;
use crate::bnd_lib_add2d;
use crate::brep_adaptor_curve2d::Curve2dAdaptor;
use crate::brep_class::FaceClassifier;
use crate::brep_class3d::find_a_point_in_the_face_uv;
use crate::brep_tool::BRepTool;
use crate::fclass2d::FaceState;
use crate::shape::{Edge, Face, TopoShape, Wire};
use crate::topo_tools::{explorer_shapes, shape_map_key, ShapeMapKey};
use crate::topopebrep::tool::c2d::fc2d_has_curve_on_surface;
use crate::topopebrep::tool::face::ToolFace;
use crate::topopebrep::tool::geometry::tol_uv;
use crate::topopebrep::tool::topology::fun_tool_par_on_e_f_tol;

/// `#define SAME (-1)` (`cxx:36`).
pub const SAME: i32 = -1;
/// `#define DIFF (-2)` (`cxx:37`).
pub const DIFF: i32 = -2;
/// `#define UNKNOWN (0)` (`cxx:38`).
pub const UNKNOWN: i32 = 0;
/// `#define oneINtwo (1)` (`cxx:39`).
pub const ONE_IN_TWO: i32 = 1;
/// `#define twoINone (2)` (`cxx:40`).
pub const TWO_IN_ONE: i32 = 2;

/// `FUN_addOwlw(Ow, lw, lresu)` (`cxx:484-500`): `Ow` when it has no splits,
/// its splits otherwise.
pub fn fun_add_owlw(ow: &TopoShape, lw: &[TopoShape], lresu: &mut Vec<TopoShape>) {
    if lw.is_empty() {
        lresu.push(ow.clone());
    } else {
        lresu.extend(lw.iter().cloned());
    }
}

/// `FUN_thegreatest(F1, class2)` (`cxx:219-250`).
///
/// `class2` already classified the point of `ffi2` against `ffi1`; this scans
/// the edges of `F1` for a mid-edge point that is IN (answers `oneINtwo`) or
/// OUT (answers `twoINone`) of `class2`.
fn fun_thegreatest(f1: &TopoShape, class2: &FaceClassifier) -> i32 {
    let tol_f1 = BRepTool::face_tolerance(&Face(f1.clone()));
    for e1 in explorer_shapes(f1, ShapeType::Edge) {
        let (f1_par, l1_par) = BRepTool::edge_parameters(&Edge(e1.clone()));
        let x = 0.45678;
        let p1 = (1.0 - x) * f1_par + x + l1_par;
        let Some(_uv1) = fun_tool_par_on_e_f_tol(&e1, p1, f1, tol_f1) else {
            continue;
        };
        match class2.state() {
            FaceState::In => return ONE_IN_TWO,
            FaceState::Out => return TWO_IN_ONE,
            _ => {}
        }
    }
    UNKNOWN
}

/// `class TopOpeBRepTool_CLASSI` (`hxx:35-71`).
#[derive(Default)]
pub struct Classi {
    /// `myFref` (`hxx:68`).
    my_f_ref: Option<TopoShape>,
    /// `mymapsbox2d` (`hxx:69`): `Add2d` order is not observable, only
    /// `Contains` / `FindFromKey` are used.
    my_maps_box2d: HashMap<ShapeMapKey, BndBox2d>,
    /// `mymapsface` (`hxx:70`).
    my_maps_face: HashMap<ShapeMapKey, ToolFace>,
}

impl Classi {
    /// `TopOpeBRepTool_CLASSI()` (`cxx:47`).
    pub fn new() -> Self {
        Self::default()
    }

    /// `Init2d(Fref)` (`cxx:52-55`).
    pub fn init_2d(&mut self, f_ref: &TopoShape) {
        self.my_f_ref = Some(f_ref.clone());
    }

    /// `HasInit2d()` (`cxx:59-62`).
    pub fn has_init_2d(&self) -> bool {
        self.my_f_ref.is_some()
    }

    /// `Add2d(S)` (`cxx:66-96`): the UV bounding box of `S` on `myFref`,
    /// cached. Answers `false` when the face has not been initialised or an
    /// edge of `S` has no p-curve on it (`cxx:69`, `:86-89`).
    pub fn add_2d(&mut self, s: &TopoShape) -> bool {
        let Some(f_ref) = self.my_f_ref.clone() else {
            return false;
        };
        let key = shape_map_key(s);
        if self.my_maps_box2d.contains_key(&key) {
            return true;
        }

        let mut b2d = BndBox2d::new();
        for e in explorer_shapes(s, ShapeType::Edge) {
            let tol_e = BRepTool::edge_tolerance(&Edge(e.clone()));
            // `cxx:87-89`
            if !fc2d_has_curve_on_surface(&e, &f_ref) {
                return false;
            }
            let adaptor = Curve2dAdaptor::of_edge_face(&Edge(e.clone()), &Face(f_ref.clone()));
            // `cxx:90-92`: `BC2d.Resolution(tolE)` then
            // `BndLib_Add2dCurve::Add(BC2d, tol2d, B2d)`. The `haspc` test
            // above already guarantees a p-curve; an adaptor that still has
            // none reports no resolution, which OCCT's `IsNull` check turns
            // into the same `false` (`cxx:88`).
            let Some(tol2d) = adaptor.resolution(tol_e) else {
                return false;
            };
            match adaptor.curve() {
                Some(c) => bnd_lib_add2d::add_adaptor(c.as_ref(), tol2d, &mut b2d),
                None => return false,
            }
        }
        self.my_maps_box2d.insert(key, b2d);
        true
    }

    /// `GetBox2d(S, B2d)` (`cxx:100-113`).
    pub fn get_box_2d(&mut self, s: &TopoShape, b2d: &mut BndBox2d) -> bool {
        let key = shape_map_key(s);
        if !self.my_maps_box2d.contains_key(&key) && !self.add_2d(s) {
            return false;
        }
        match self.my_maps_box2d.get(&key) {
            Some(b) => {
                *b2d = b.clone();
                true
            }
            None => false,
        }
    }

    /// `ClassiBnd2d(S1, S2, tol, chklarge)` (`cxx:118-217`).
    pub fn classi_bnd_2d(
        &mut self,
        s1: &TopoShape,
        s2: &TopoShape,
        tol: f64,
        chklarge: bool,
    ) -> i32 {
        let mut b: [BndBox2d; 2] = [BndBox2d::new(), BndBox2d::new()];
        for (i, s) in [s1, s2].iter().enumerate() {
            let key = shape_map_key(s);
            if !self.my_maps_box2d.contains_key(&key) && !self.add_2d(s) {
                // `cxx:129-131`, `:139-141`: `return false` from an `int`
                // function, i.e. `UNKNOWN` (0).
                return UNKNOWN;
            }
            match self.my_maps_box2d.get(&key) {
                Some(box2d) => b[i] = box2d.clone(),
                None => return UNKNOWN,
            }
        }

        // `cxx:145-151`: `(Umin(i), Vmin(i), Umax(i), Vmax(i))`.
        let mut uv = [[0.0_f64; 4]; 2];
        for i in 0..2 {
            match b[i].get() {
                // `Bnd_Box2d::Get` (no gap): `xmin, ymin, xmax, ymax`.
                Some((xmin, ymin, xmax, ymax)) => uv[i] = [xmin, ymin, xmax, ymax],
                // A void box has no bounds. OCCT would read garbage; there is
                // no OCCT branch for it, so the comparison cannot be made.
                None => return UNKNOWN,
            }
        }

        // `cxx:163-179`: disjoint boxes.
        for k in [0usize, 2] {
            for i in 0..2 {
                let ii = i;
                let jj = if i == 0 { 1 } else { 0 };
                //  diff = Umin<ii> - Umax<jj> : k = 0
                //  diff = Vmin<ii> - Vmax<jj> : k = 2
                let diff = uv[ii][k] - uv[jj][k + 1];
                let disjoint = if chklarge { diff >= -tol } else { diff > 0.0 };
                if disjoint {
                    return DIFF;
                }
            }
        }

        // `cxx:181-216`: which box is IN the other.
        for i in 0..2 {
            let ii = i;
            let jj = if i == 0 { 1 } else { 0 };
            let mut smaller = true;
            let mut same = true;

            for k in [0usize, 2] {
                let diff = uv[ii][k] - uv[jj][k];
                smaller = if chklarge { smaller && (diff > -tol) } else { smaller && (diff > 0.0) };
                same = same && (diff.abs() <= tol);
            }
            for k in [1usize, 3] {
                let diff = uv[ii][k] - uv[jj][k];
                smaller = if chklarge { smaller && (diff < tol) } else { smaller && (diff < 0.0) };
                same = same && (diff.abs() <= tol);
            }

            if same {
                return SAME;
            }
            if smaller {
                return if ii == 0 { ONE_IN_TWO } else { TWO_IN_ONE };
            }
        }
        UNKNOWN
    }

    /// `Classip2d(S1, S2, stabnd2d12)` (`cxx:254-467`).
    pub fn classip_2d(&mut self, s1: &TopoShape, s2: &TopoShape, stabnd2d12: i32) -> i32 {
        let Some(f_ref) = self.my_f_ref.clone() else {
            return UNKNOWN;
        };
        let bnd2d_unk = stabnd2d12 == UNKNOWN || stabnd2d12 == SAME;

        // `cxx:266-279`: `fa1.Init(Wire(S1), myFref)`.
        let key1 = shape_map_key(s1);
        let fa1 = match self.my_maps_face.get(&key1) {
            Some(fa) => fa.clone(),
            None => {
                let mut fa = ToolFace::new();
                let wire = match Wire::wrap(s1.clone()) {
                    Some(w) => w,
                    None => return UNKNOWN,
                };
                if !fa.init(&wire.0, &f_ref) {
                    return UNKNOWN;
                }
                self.my_maps_face.insert(key1, fa.clone());
                fa
            }
        };
        let Some(ffi1) = fa1.f_finite().cloned() else {
            return UNKNOWN;
        };
        let finite1 = fa1.finite();

        // `cxx:284-299`: `fa2.Init(Wire(S2), myFref)`.
        let key2 = shape_map_key(s2);
        let fa2 = match self.my_maps_face.get(&key2) {
            Some(fa) => fa.clone(),
            None => {
                let mut fa = ToolFace::new();
                let wire = match Wire::wrap(s2.clone()) {
                    Some(w) => w,
                    None => return UNKNOWN,
                };
                if !fa.init(&wire.0, &f_ref) {
                    return UNKNOWN;
                }
                self.my_maps_face.insert(key2, fa.clone());
                fa
            }
        };
        let Some(ffi2) = fa2.f_finite().cloned() else {
            return UNKNOWN;
        };
        let finite2 = fa2.finite();

        // `cxx:303-310`: `p2d1` in `ffi1`.
        let mut u1 = 0.0_f64;
        let mut v1 = 0.0_f64;
        if !find_a_point_in_the_face_uv(&Face(ffi1.clone()), &mut u1, &mut v1) {
            return UNKNOWN;
        }
        let p2d1 = GpPnt2d::new(u1, v1);

        // `cxx:315-318`: `class2` classifies `p2d1` against `ffi2`.
        let tol2d2 = tol_uv(&Face(ffi2.clone()), BRepTool::face_tolerance(&Face(ffi2.clone())));
        let mut class2 = FaceClassifier::new();
        class2.perform(&Face(ffi2.clone()), p2d1, tol2d2);
        let sta12 = class2.state();

        // `cxx:320-368`.
        let mut staffi12 = UNKNOWN;
        let mut try21 = false;
        if bnd2d_unk {
            try21 = true;
        } else if stabnd2d12 == DIFF {
            return DIFF;
        } else if sta12 == FaceState::In {
            if stabnd2d12 == ONE_IN_TWO {
                staffi12 = ONE_IN_TWO;
            } else if stabnd2d12 == TWO_IN_ONE {
                staffi12 = TWO_IN_ONE;
            } else {
                try21 = true;
            }
        } else if sta12 == FaceState::Out {
            if stabnd2d12 == TWO_IN_ONE {
                staffi12 = TWO_IN_ONE;
            } else {
                try21 = true;
            }
        } else {
            // `cxx:364-367`: `NYIxpu140199`
            return UNKNOWN;
        }

        if try21 {
            // `cxx:372-379`: `p2d2` in `ffi2`.
            let mut u2 = 0.0_f64;
            let mut v2 = 0.0_f64;
            if !find_a_point_in_the_face_uv(&Face(ffi2.clone()), &mut u2, &mut v2) {
                return UNKNOWN;
            }
            let p2d2 = GpPnt2d::new(u2, v2);

            // `cxx:382-385`
            let tol2d1 =
                tol_uv(&Face(ffi1.clone()), BRepTool::face_tolerance(&Face(ffi1.clone())));
            let mut class1 = FaceClassifier::new();
            class1.perform(&Face(ffi1.clone()), p2d2, tol2d1);
            let sta21 = class1.state();

            if bnd2d_unk {
                // `cxx:387-404`
                if sta12 == FaceState::Out && sta21 == FaceState::Out {
                    staffi12 = DIFF;
                } else if sta12 == FaceState::In && sta21 == FaceState::Out {
                    staffi12 = ONE_IN_TWO;
                } else if sta12 == FaceState::Out && sta21 == FaceState::In {
                    staffi12 = TWO_IN_ONE;
                } else if sta12 == FaceState::In && sta21 == FaceState::In {
                    staffi12 = fun_thegreatest(&ffi1, &class2);
                }
            } else if sta12 == FaceState::In {
                // `cxx:408-412`
                if sta21 == FaceState::Out {
                    staffi12 = ONE_IN_TWO;
                }
            } else if sta12 == FaceState::Out {
                // `cxx:414-427`
                if sta21 == FaceState::Out {
                    staffi12 = DIFF;
                } else if sta21 == FaceState::In {
                    staffi12 = TWO_IN_ONE;
                }
            } else {
                // `cxx:428`: `NYIxpu140199`
                return UNKNOWN;
            }
        }

        // `cxx:433-466`
        if staffi12 == DIFF {
            if finite1 && finite2 {
                return DIFF;
            }
            if finite1 && !finite2 {
                return TWO_IN_ONE;
            }
            if !finite1 && finite2 {
                return ONE_IN_TWO;
            }
            // `cxx:447`: `!!!!!!! or twoINone`
            return ONE_IN_TWO;
        } else if staffi12 == ONE_IN_TWO {
            if !finite1 && !finite2 {
                return TWO_IN_ONE;
            }
            return ONE_IN_TWO;
        } else if staffi12 == TWO_IN_ONE {
            if !finite1 && !finite2 {
                return ONE_IN_TWO;
            }
            return TWO_IN_ONE;
        }
        UNKNOWN
    }

    /// `Getface(S, fa)` (`cxx:470-479`).
    pub fn get_face(&self, s: &TopoShape) -> Option<ToolFace> {
        self.my_maps_face.get(&shape_map_key(s)).cloned()
    }

    /// `Classilist(lS, mapgreasma)` (`cxx:502-615`): group the wires of `lS`
    /// into `(greater, its smaller wires)`.
    ///
    /// `nite` is compared against `nitemax` but never incremented in OCCT
    /// (`cxx:508`, `:525`); the loop ends when `lw` shrinks to one wire, which
    /// the port keeps as-is.
    pub fn classilist(
        &mut self,
        ls: &[TopoShape],
        map_greater_smaller: &mut HashMap<ShapeMapKey, Vec<TopoShape>>,
    ) -> bool {
        let Some(f_ref) = self.my_f_ref.clone() else {
            return false;
        };
        let tol_ref = BRepTool::face_tolerance(&Face(f_ref.clone()));
        // `cxx:508`: `nyixpu : cdliser??`
        let tol_uv_ref = tol_uv(&Face(f_ref), tol_ref);

        map_greater_smaller.clear();
        for w in ls {
            map_greater_smaller.insert(shape_map_key(w), Vec::new());
        }

        let mut lw: Vec<TopoShape> = ls.to_vec();
        let nw_start = lw.len();
        if nw_start <= 1 {
            return true;
        }
        let nite_max = (nw_start * (nw_start - 1)) / 2;
        let nite = 0usize;

        while nite <= nite_max {
            if lw.len() <= 1 {
                break;
            }

            // `cxx:537-545`: `wi1` = the first of `lw` still bound in `map`.
            let mut k = 0usize;
            let mut wi1 = lw[0].clone();
            while k < lw.len() {
                wi1 = lw[k].clone();
                if map_greater_smaller.contains_key(&shape_map_key(&wi1)) {
                    break;
                }
                k += 1;
            }

            // `cxx:549-610`: compare `wi1` with the wires after it.
            loop {
                if k >= lw.len() {
                    break;
                }
                if !map_greater_smaller.contains_key(&shape_map_key(&wi1)) {
                    break;
                }
                // `itw.Next()` (`cxx:557`)
                k += 1;
                if k >= lw.len() {
                    break;
                }

                let mut sta12 = UNKNOWN;
                let mut out_all = false;
                let mut wi2 = lw[k].clone();
                while k < lw.len() {
                    wi2 = lw[k].clone();
                    if !map_greater_smaller.contains_key(&shape_map_key(&wi2)) {
                        k += 1;
                        continue;
                    }
                    let stabnd2d12 = self.classi_bnd_2d(&wi1, &wi2, tol_uv_ref, true);
                    sta12 = self.classip_2d(&wi1, &wi2, stabnd2d12);
                    if sta12 == DIFF {
                        out_all = true;
                        k += 1;
                        continue;
                    }
                    break;
                }

                if sta12 == ONE_IN_TWO {
                    // `cxx:588-595`: greater = shape two.
                    let mut lsma = Vec::new();
                    if let Some(item) = map_greater_smaller.get(&shape_map_key(&wi1)) {
                        fun_add_owlw(&wi1, item, &mut lsma);
                    }
                    if let Some(lwgre) = map_greater_smaller.get_mut(&shape_map_key(&wi2)) {
                        lwgre.extend(lsma);
                    }
                    map_greater_smaller.remove(&shape_map_key(&wi1));
                } else if sta12 == TWO_IN_ONE {
                    // `cxx:596-601`: greater = shape one.
                    let mut lsma = Vec::new();
                    if let Some(item) = map_greater_smaller.get(&shape_map_key(&wi2)) {
                        fun_add_owlw(&wi2, item, &mut lsma);
                    }
                    if let Some(lwgre) = map_greater_smaller.get_mut(&shape_map_key(&wi1)) {
                        lwgre.extend(lsma);
                    }
                    map_greater_smaller.remove(&shape_map_key(&wi2));
                } else if out_all {
                    // `cxx:603-606`: nothing's done.
                } else {
                    return false;
                }
            }

            lw.remove(0);
        }
        true
    }
}
