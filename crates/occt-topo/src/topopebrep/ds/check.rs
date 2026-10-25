//! `TopOpeBRepDS_Check` — a tool verifying the integrity and the structure of
//! the data structure.
//!
//! Source: `TKBool/TopOpeBRepDS/TopOpeBRepDS_Check.hxx:34-108`,
//! `TopOpeBRepDS_Check.cxx:20-656`, `TopOpeBRepDS_CheckStatus.hxx:19-23`.
//!
//! `ChkIntg` walks every interference list of the DS and, through
//! [`Check::check_ds`], records for each support / geometry whether it exists
//! and whether a topological index carries the shape type its `Kind` names.
//! `ChkIntgSamDom` verifies the same-domain tables and `OneVertexOnPnt` the
//! vertex / point duplication rules.
//!
//! Translation notes:
//! - `HDataStructure` is borrowed (`Check<'a>`), not held as a mutable handle:
//!   every method of the OCCT class only reads it, and `HDS()` returns the
//!   handle by reference. OCCT's default constructor leaves a null handle that
//!   every method dereferences, so the port only offers [`Check::new`].
//! - `NCollection_DataMap<int, TopOpeBRepDS_CheckStatus>` is a `BTreeMap`, so
//!   `PrintMap` walks the indices in ascending order; the OCCT hash order is
//!   unspecified and not reproducible (same choice as `ds/tki.rs`).
//! - `Standard_OStream&` becomes `&mut String`.

use std::collections::BTreeMap;

use crate::abs::ShapeType;
use crate::shape::TopoShape;
use crate::topopebrep::tool::shape_tool;

use super::h_data_structure::HDataStructure;
use super::interference::{InterferenceData, InterferenceHandle};
use super::interference_tool;
use super::kind::Kind;
use super::point_explorer::PointExplorer;

/// `TopOpeBRepDS_CheckStatus` (`CheckStatus.hxx:19-23`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckStatus {
    /// `TopOpeBRepDS_OK`
    Ok,
    /// `TopOpeBRepDS_NOK`
    Nok,
}

/// `class TopOpeBRepDS_Check` (`hxx:34-108`).
pub struct Check<'a> {
    /// `myHDS`
    hds: &'a HDataStructure,
    /// `myMapSurfaceStatus`
    map_surface_status: BTreeMap<i32, CheckStatus>,
    /// `myMapCurveStatus`
    map_curve_status: BTreeMap<i32, CheckStatus>,
    /// `myMapPointStatus`
    map_point_status: BTreeMap<i32, CheckStatus>,
    /// `myMapShapeStatus`
    map_shape_status: BTreeMap<i32, CheckStatus>,
}

impl<'a> Check<'a> {
    /// `TopOpeBRepDS_Check(const handle<TopOpeBRepDS_HDataStructure>&)`
    /// (`cxx:42-49`).
    pub fn new(hds: &'a HDataStructure) -> Self {
        Self {
            hds,
            map_surface_status: BTreeMap::new(),
            map_curve_status: BTreeMap::new(),
            map_point_status: BTreeMap::new(),
            map_shape_status: BTreeMap::new(),
        }
    }

    /// `ChkIntg()` (`cxx:64-103`).
    ///
    /// The first loop writes `bI` (plain assignment, `cxx:75`), every later
    /// loop *ands* into it (`cxx:82`, `:91`, `:98`); the result of
    /// `CheckEdgeParameter(myHDS)` (`cxx:100`) is discarded. Both are kept.
    pub fn chk_intg(&mut self) -> bool {
        let hds = self.hds;
        let ds = hds.ds();
        let mut b_i = false;

        let nshape = ds.nb_shapes();
        for i in 1..=nshape {
            b_i = self.chk_intg_interf(ds.shape_interferences_at(i, false));
        }
        let nsurface = ds.nb_surfaces();
        for i in 1..=nsurface {
            b_i = b_i && self.chk_intg_interf(ds.surface_interferences(i));
        }
        let ncurve = ds.nb_curves();
        for i in 1..=ncurve {
            b_i = b_i && self.chk_intg_interf(ds.curve_interferences(i));
        }
        let npoint = ds.nb_points();
        for i in 1..=npoint {
            b_i = b_i && self.chk_intg_interf(ds.point_interferences(i));
        }

        // `CheckEdgeParameter(myHDS);` (`cxx:100`): the result is dropped.
        let _ = check_edge_parameter(hds);

        b_i
    }

    /// `ChkIntgInterf(LI)` (`cxx:109-123`): both the support and the geometry
    /// of every interference are checked with [`Check::check_ds`].
    pub fn chk_intg_interf(&mut self, li: &[InterferenceHandle]) -> bool {
        let mut is_ok = true;
        for ih in li {
            let Ok(i1) = ih.read() else { continue };
            is_ok = is_ok && self.check_ds(i1.support(), i1.support_type());
            is_ok = is_ok && self.check_ds(i1.geometry(), i1.geometry_type());
        }
        is_ok
    }

    /// `CheckDS(I, K)` (`cxx:127-238`).
    ///
    /// `SURFACE` / `CURVE` / `POINT` are checked against the geometry counts and
    /// return; every other kind falls through to the topology test, where the
    /// `switch` (`cxx:199-221`) leaves `se = TopAbs_COMPOUND` for a kind it
    /// does not name.
    pub fn check_ds(&mut self, i: i32, k: Kind) -> bool {
        let nb = match k {
            Kind::Surface => Some(self.hds.nb_surfaces()),
            Kind::Curve => Some(self.hds.nb_curves()),
            Kind::Point => Some(self.hds.nb_points()),
            _ => None,
        };
        if let Some(nb) = nb {
            let map = match k {
                Kind::Surface => &mut self.map_surface_status,
                Kind::Curve => &mut self.map_curve_status,
                _ => &mut self.map_point_status,
            };
            if nb < i {
                map.remove(&i);
                map.insert(i, CheckStatus::Nok);
                return false;
            }
            if !map.contains_key(&i) {
                map.insert(i, CheckStatus::Ok);
            }
            return true;
        }

        if self.hds.nb_shapes() < i {
            self.map_shape_status.remove(&i);
            self.map_shape_status.insert(i, CheckStatus::Nok);
            return false;
        }
        let Some(s) = self.hds.shape(i, false) else {
            self.map_shape_status.remove(&i);
            self.map_shape_status.insert(i, CheckStatus::Nok);
            return false;
        };

        let se = match k {
            Kind::Solid => ShapeType::Solid,
            Kind::Shell => ShapeType::Shell,
            Kind::Face => ShapeType::Face,
            Kind::Wire => ShapeType::Wire,
            Kind::Edge => ShapeType::Edge,
            Kind::Vertex => ShapeType::Vertex,
            _ => ShapeType::Compound,
        };
        if s.shape_type() != se {
            self.map_shape_status.remove(&i);
            self.map_shape_status.insert(i, CheckStatus::Nok);
            return false;
        }
        if !self.map_shape_status.contains_key(&i) {
            self.map_shape_status.insert(i, CheckStatus::Ok);
        }
        true
    }

    /// `ChkIntgSamDom()` (`cxx:240-297`).
    ///
    /// `bb` is declared once for the whole function (`cxx:242`) and is never
    /// reset, so after the first `Sind` found inside `ShapeSameDomain(Sref)`
    /// the later iterations skip the `!bb` failure (`cxx:290-293`). Kept.
    pub fn chk_intg_sam_dom(&self) -> bool {
        let hds = self.hds;
        let ds = hds.ds();
        let mut b = true;
        let mut bb = false;
        let nb_sh = hds.nb_shapes();
        for i in 1..=nb_sh {
            let Some(sind) = hds.shape(i, false) else {
                continue;
            };
            if !self.check_shapes(ds.shape_same_domain_at(i)) {
                b = false;
            }

            // Verification de SameDomaineRef (`cxx:256-261`).
            let curr = ds.same_domain_ref(i);
            let loc = ds.same_domain_ref(curr);
            if curr != 0 && curr != loc {
                b = false;
            }

            if curr != 0 {
                let Some(sref) = hds.shape(curr, false) else {
                    b = false;
                    continue;
                };
                if sind.shape_type() != sref.shape_type() {
                    b = false;
                }
                if i != curr {
                    let losr = ds.shape_same_domain_at(curr);
                    for sh in losr {
                        let loc = hds.shape_index(sh, false);
                        if loc == i {
                            bb = true;
                            break;
                        }
                    }
                    if !bb {
                        b = false;
                    }
                }
            }
        }
        b
    }

    /// `CheckShapes(LS)` (`cxx:303-316`): every shape must be present in the
    /// DS; an empty list is `true`.
    pub fn check_shapes(&self, ls: &[TopoShape]) -> bool {
        for it_s in ls {
            if self.hds.shape_index(it_s, false) == 0 {
                return false;
            }
        }
        true
    }

    /// `OneVertexOnPnt()` (`cxx:322-400`).
    ///
    /// `j < NbPo` in the point / point loop (`cxx:389`) is strictly less, so the
    /// last point is never paired; the `IsEqual` body is empty in OCCT
    /// (`cxx:391-394`) and stays empty.
    ///
    /// OCCT raises `Standard_ProgramError` from
    /// `TopOpeBRepTool_ShapeTool::Tolerance` / `Pnt` (`ShapeTool.cxx:46-80`) for
    /// a shape that is neither face, edge nor vertex; the loop only ever holds
    /// vertices, and the port reports such a case as `false`.
    pub fn one_vertex_on_pnt(&self) -> bool {
        let hds = self.hds;
        let ds = hds.ds();
        let mut b = true;
        let nb_po = hds.nb_points();

        let mut vert: Vec<i32> = Vec::new();
        for i in 1..=hds.nb_shapes() {
            let Some(s) = hds.shape(i, false) else {
                continue;
            };
            if s.shape_type() == ShapeType::Vertex && hds.has_shape(&s, true) {
                vert.push(i);
            }
        }

        let nb_ve = vert.len();
        for ii in 0..nb_ve {
            let curr1 = vert[ii];
            let Some(s1) = hds.shape(curr1, false) else {
                continue;
            };
            let sdr1 = hds.same_domain_reference(&s1);
            for jj in (ii + 1)..nb_ve {
                let curr2 = vert[jj];
                let Some(s2) = hds.shape(curr2, false) else {
                    continue;
                };
                let sdr2 = hds.same_domain_reference(&s2);
                let Ok(tol1) = shape_tool::tolerance_of(&s1) else {
                    return false;
                };
                let Ok(tol2) = shape_tool::tolerance_of(&s2) else {
                    return false;
                };
                let (Ok(p1), Ok(p2)) = (shape_tool::pnt_of(&s1), shape_tool::pnt_of(&s2)) else {
                    return false;
                };
                let dist = p1.distance(&p2);
                if dist <= tol1 + tol2 {
                    if sdr1 != sdr2 {
                        b = false;
                    }
                } else if sdr1 == sdr2 {
                    b = false;
                }
            }

            let mut pe = PointExplorer::with_ds(ds, true);
            while pe.more() {
                if let Some(ds_pnt) = pe.point() {
                    let pnt1 = *ds_pnt.point();
                    let tol1 = ds_pnt.tolerance();
                    let Ok(tol2) = shape_tool::tolerance_of(&s1) else {
                        return false;
                    };
                    let Ok(pnt2) = shape_tool::pnt_of(&s1) else {
                        return false;
                    };
                    let dist = pnt1.distance(&pnt2);
                    if dist <= tol1 + tol2 {
                        b = false;
                    }
                }
                pe.next();
            }
        }

        for i in 1..=nb_po {
            let pe = PointExplorer::with_ds(ds, true);
            if pe.is_point(i) {
                let Some(ds_pnt1) = hds.point(i) else {
                    continue;
                };
                for j in (i + 1)..nb_po {
                    let Some(ds_pnt2) = hds.point(j) else {
                        continue;
                    };
                    if ds_pnt1.is_equal(ds_pnt2) {
                        // `cxx:391-394`: the OCCT body is empty.
                    }
                }
            }
        }

        b
    }

    /// `HDS()` (`cxx:648-651`).
    pub fn hds(&self) -> &'a HDataStructure {
        self.hds
    }

    /// `ChangeHDS()` (`cxx:655-658`) re-binds the handle in OCCT; the borrowed
    /// form re-binds the reference instead.
    pub fn set_hds(&mut self, hds: &'a HDataStructure) {
        self.hds = hds;
    }

    /// `PrintIntg(S)` (`cxx:462-510`).
    pub fn print_intg(&self, s: &mut String) {
        s.push('\n');
        s.push('\n');
        s.push_str("************************************************");
        s.push('\n');
        s.push_str("state of the DS : (only the tested elements)");
        s.push('\n');
        s.push('\n');

        self.print_map(&self.map_surface_status, "Surface", s);
        self.print_map(&self.map_curve_status, "Curve", s);
        self.print_map(&self.map_point_status, "Point", s);

        let mut map_vertex: BTreeMap<i32, CheckStatus> = BTreeMap::new();
        let mut map_edge: BTreeMap<i32, CheckStatus> = BTreeMap::new();
        let mut map_wire: BTreeMap<i32, CheckStatus> = BTreeMap::new();
        let mut map_face: BTreeMap<i32, CheckStatus> = BTreeMap::new();
        let mut map_solid: BTreeMap<i32, CheckStatus> = BTreeMap::new();
        for (i, st) in &self.map_shape_status {
            let Some(sh) = self.hds.shape(*i, false) else {
                continue;
            };
            match sh.shape_type() {
                ShapeType::Vertex => map_vertex.insert(*i, *st),
                ShapeType::Edge => map_edge.insert(*i, *st),
                ShapeType::Wire => map_wire.insert(*i, *st),
                ShapeType::Face => map_face.insert(*i, *st),
                ShapeType::Solid => map_solid.insert(*i, *st),
                _ => None,
            };
        }

        self.print_map(&map_vertex, "Vertex", s);
        self.print_map(&map_edge, "Edge", s);
        self.print_map(&map_wire, "Wire", s);
        self.print_map(&map_face, "Face", s);
        self.print_map(&map_solid, "Solid", s);

        s.push('\n');
        s.push('\n');
    }

    /// `PrintMap(MapStat, eltstr, S)` (`cxx:515-538`).
    pub fn print_map(&self, map_stat: &BTreeMap<i32, CheckStatus>, eltstr: &str, s: &mut String) {
        if map_stat.is_empty() {
            return;
        }
        let mut b = true;
        s.push(' ');
        s.push_str(eltstr);
        s.push_str("\t(/");
        s.push_str(&map_stat.len().to_string());
        s.push_str(")\tnumber ");

        self.print_elts(map_stat, CheckStatus::Ok, &mut b, s);
        if !b {
            s.push_str(" = OK");
            s.push('\n');
        }
        self.print_elts(map_stat, CheckStatus::Nok, &mut b, s);
        if !b {
            s.push_str(" = NOK");
            s.push('\n');
        }
    }

    /// `PrintElts(MapStat, Stat, b, S)` (`cxx:547-575`).
    ///
    /// `b` is in / out: `bb` starts as `!b`, `b` is cleared, the first matching
    /// element clears `b` and the "\t\tnumber " prefix is only written when the
    /// incoming `b` was `false`.
    pub fn print_elts(
        &self,
        map_stat: &BTreeMap<i32, CheckStatus>,
        stat: CheckStatus,
        b: &mut bool,
        s: &mut String,
    ) {
        let mut bb = !*b;
        *b = true;
        for (i, st) in map_stat {
            if *st == stat {
                if *b {
                    *b = false;
                }
                if bb {
                    s.push_str("\t\tnumber ");
                    bb = false;
                }
                s.push_str(&i.to_string());
                s.push(' ');
            }
        }
    }

    /// `Print(stat, S)` (`cxx:579-594`).
    pub fn print(&self, stat: CheckStatus, s: &mut String) {
        match stat {
            CheckStatus::Ok => s.push_str("OK"),
            CheckStatus::Nok => s.push_str("NOK"),
        }
    }

    /// `PrintShape(SE, S)` (`cxx:598-625`).
    pub fn print_shape(&self, se: ShapeType, s: &mut String) {
        match se {
            ShapeType::Solid => s.push_str("Solid "),
            ShapeType::Shell => s.push_str("Shell "),
            ShapeType::Face => s.push_str("Face  "),
            ShapeType::Wire => s.push_str("Wire  "),
            ShapeType::Edge => s.push_str("Edge  "),
            ShapeType::Vertex => s.push_str("Vertex"),
            _ => {}
        }
    }

    /// `PrintShape(index, S)` (`cxx:629-644`): the final `myHDS->Shape(index)`
    /// of `cxx:640` prints nothing, so only the two diagnostics reach the
    /// stream.
    pub fn print_shape_index(&self, index: i32, s: &mut String) {
        if self.hds.nb_shapes() < index {
            s.push_str("**PB**IN**TopOpeBRepDS_Check::PrintShape** ");
            return;
        }
        if !self.map_shape_status.contains_key(&index) {
            s.push_str("NO CHECK HAS PROCESSING");
            s.push('\n');
            return;
        }
        s.push(' ');
        let _ = self.hds.shape(index, false);
    }
}

/// `CheckEdgeParameter(const handle<TopOpeBRepDS_HDataStructure>&)`
/// (`cxx:406-456`).
///
/// OCCT casts the parameter to `int` before comparing it with `1.e50`
/// (`cxx:427` and `cxx:448`): a finite parameter can never survive that cast
/// and stay above `1.e50`, so the `IsOK = false` arms are unreachable. The port
/// keeps the cast and the comparison as written; Rust's `as i32` saturates where
/// the C++ cast is undefined.
pub fn check_edge_parameter(hds: &HDataStructure) -> bool {
    let ds = hds.ds();
    let mut is_ok = true;

    let nshape = ds.nb_shapes();
    for i in 1..=nshape {
        for ih in ds.shape_interferences_at(i, false) {
            let Ok(i1) = ih.read() else { continue };
            if let Some(evi) = i1.as_edge_vertex() {
                let param = evi.param as i32;
                if f64::from(param) > 1.0e50 {
                    is_ok = false;
                }
            }
        }
    }

    let ncurve = ds.nb_curves();
    for i in 1..=ncurve {
        for ih in ds.curve_interferences(i) {
            let Ok(i1) = ih.read() else { continue };
            if matches!(i1.data(), InterferenceData::CurvePoint { .. }) {
                let param = interference_tool::parameter(&i1) as i32;
                if f64::from(param) > 1.0e50 {
                    is_ok = false;
                }
            }
        }
    }
    is_ok
}
