//! `TopOpeBRepBuild_Builder` — the five KPart predicates.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_KPart.cxx`
//! (OCCT 8.0.0):
//! - `KPisdisj` (`:1389-1421`) / `KPisdisjsh` (`:1786-1848`),
//! - `KPisfafa` (`:1428-1433`) / `KPisfafash` (`:1940-2005`),
//! - `KPissoso` (`:1442-1447`) / `KPissososh` (`:1856-1931`),
//! - `KPiskoletge` (`:1253-1385`) / `KPiskoletgesh` (`:1620-1707`) with
//!   `sectionedgesON` (`:1195-1220`) and `allIonsectionedges` (`:1223-1246`),
//! - `KPiskoletgeanalyse` (`:2010-2096`).
//!
//! `FUNKP_KPiskolesh` (`kpkole.cxx:1033-1075`) is translated here because
//! `KPiskoletgesh` calls it; `KPiskolesh` / `KPiskole` / `MergeKPartiskole`
//! (`kpkole.cxx:1085-1157`, `:576-732`, `:56-572`) stay UNPORTED.

use std::collections::HashSet;

use crate::abs::{ShapeType, State};
use crate::brep_tools;
use crate::shape::{Face, TopoShape};
use crate::topo_tools::{
    explorer_shapes, explorer_shapes_avoiding, shape_map_key, ShapeMapKey,
};
use crate::topopebrep::ds::kind::Config;
use crate::topopebrep::tool::geometry::{fun_tool_cylinder, fun_tool_plane};
use crate::topopebrep::tool::topology::fun_tool_in_s;

use super::super::k_part::{RESNEWSOL, RESNULL, RESSHAPE1, RESSHAPE2, RESUNDEF};
use super::Builder;

/// `FUN_sortplcy` (`KPart.cxx:210-229`): split `lof` into its planes and its
/// cylinders.
pub(super) fn fun_sortplcy(lof: &[TopoShape]) -> (Vec<TopoShape>, Vec<TopoShape>) {
    let mut lopl = Vec::new();
    let mut locy = Vec::new();
    for ff in lof {
        let Some(face) = Face::wrap(ff.clone()) else {
            continue;
        };
        if fun_tool_plane(&face) {
            lopl.push(ff.clone());
        }
        if fun_tool_cylinder(&face) {
            locy.push(ff.clone());
        }
    }
    (lopl, locy)
}

/// `TopExp::MapShapes(w, TopAbs_EDGE, mape)` (`KPart.cxx:1333-1337`).
fn edge_set(s: &TopoShape) -> HashSet<ShapeMapKey> {
    explorer_shapes(s, ShapeType::Edge)
        .iter()
        .map(shape_map_key)
        .collect()
}

/// `BRepTools::OuterWire(face)`.
fn outer_wire_of(f: &TopoShape) -> Option<TopoShape> {
    let face = Face::wrap(f.clone())?;
    brep_tools::outer_wire(&face).map(|w| w.0)
}

impl Builder {
    /// `KPisdisj()` (`KPart.cxx:1389-1421`): no interference at all on either
    /// operand.
    pub fn kp_isdisj(&self) -> i32 {
        let Some(hds) = self.my_data_structure.as_ref() else {
            return 0;
        };
        // cxx:1400-1407
        let has1 = self
            .my_shape1
            .as_ref()
            .map(|s| hds.has_shape(s, true))
            .unwrap_or(false);
        if !has1 {
            return 0;
        }
        let has2 = self
            .my_shape2
            .as_ref()
            .map(|s| hds.has_shape(s, true))
            .unwrap_or(false);
        if !has2 {
            return 0;
        }

        // cxx:1409-1410
        let isdisj1 = self
            .my_shape1
            .as_ref()
            .map(|s| self.kp_isdisjsh(s))
            .unwrap_or(0);
        let isdisj2 = self
            .my_shape2
            .as_ref()
            .map(|s| self.kp_isdisjsh(s))
            .unwrap_or(0);

        // cxx:1419
        i32::from(isdisj1 != 0 && isdisj2 != 0)
    }

    /// `KPisfafa()` (`KPart.cxx:1428-1433`).
    pub fn kp_isfafa(&self) -> i32 {
        let is1 = self
            .my_shape1
            .as_ref()
            .map(|s| self.kp_isfafash(s))
            .unwrap_or(0);
        let is2 = self
            .my_shape2
            .as_ref()
            .map(|s| self.kp_isfafash(s))
            .unwrap_or(0);
        i32::from(is1 != 0 && is2 != 0)
    }

    /// `KPissoso()` (`KPart.cxx:1442-1447`).
    pub fn kp_issoso(&self) -> i32 {
        let is1 = self
            .my_shape1
            .as_ref()
            .map(|s| self.kp_issososh(s))
            .unwrap_or(0);
        let is2 = self
            .my_shape2
            .as_ref()
            .map(|s| self.kp_issososh(s))
            .unwrap_or(0);
        i32::from(is1 != 0 && is2 != 0)
    }

    /// `KPisdisjsh(Sarg)` (`KPart.cxx:1786-1848`).
    pub fn kp_isdisjsh(&self, sarg: &TopoShape) -> i32 {
        if sarg.is_null() {
            return 0;
        }

        // cxx:1796-1813: no solid, face or edge carries a geometry.
        for t in [ShapeType::Solid, ShapeType::Face, ShapeType::Edge] {
            if self.kp_lhg(sarg, t) != 0 {
                return 0;
            }
        }

        // cxx:1815-1832: one HasSameDomain level only — no face under a
        // same-domain solid carries a same-domain face.
        let mut lshsd = Vec::new();
        if self.kp_lhsd_into(sarg, ShapeType::Solid, &mut lshsd) != 0 {
            let solids = lshsd.clone();
            for s in solids {
                if self.kp_lhsd(&s, ShapeType::Face) != 0 {
                    return 0;
                }
            }
        }

        // cxx:1834-1845: no edge under a same-domain face carries a
        // same-domain edge.
        if self.kp_lhsd_into(sarg, ShapeType::Face, &mut lshsd) != 0 {
            let faces = lshsd.clone();
            for s in faces {
                if self.kp_lhsd(&s, ShapeType::Edge) != 0 {
                    return 0;
                }
            }
        }

        // cxx:1847
        1
    }

    /// `KPissososh(Sarg)` (`KPart.cxx:1856-1931`).
    pub fn kp_issososh(&self, sarg: &TopoShape) -> i32 {
        // cxx:1859-1874: every solid is a free solid.
        let nsol1 = explorer_shapes(sarg, ShapeType::Solid).len();
        let nsol2 = explorer_shapes_avoiding(sarg, ShapeType::Solid, ShapeType::CompSolid).len();
        if nsol1 != 0 && nsol1 != nsol2 {
            return 0;
        }

        // cxx:1878-1883: all of them are HasSameDomain.
        if self.kp_lhsd(sarg, ShapeType::Solid) as usize != nsol1 {
            return 0;
        }

        // cxx:1886-1930: no free shell / face / edge / vertex.
        for (find, avoid) in [
            (ShapeType::Shell, ShapeType::Solid),
            (ShapeType::Face, ShapeType::Shell),
            (ShapeType::Edge, ShapeType::Wire),
            (ShapeType::Vertex, ShapeType::Edge),
        ] {
            if !explorer_shapes_avoiding(sarg, find, avoid).is_empty() {
                return 0;
            }
        }

        // cxx:1932
        1
    }

    /// `KPisfafash(Sarg)` (`KPart.cxx:1940-2005`).
    pub fn kp_isfafash(&self, sarg: &TopoShape) -> i32 {
        // cxx:1942-1962: every face is a free face.
        let nfac1 = explorer_shapes(sarg, ShapeType::Face).len();
        let nfac2 = explorer_shapes_avoiding(sarg, ShapeType::Face, ShapeType::Shell).len();
        if nfac1 != 0 && nfac1 != nfac2 {
            return 0;
        }

        // cxx:1965-1969: all of them are HasSameDomain.
        if self.kp_lhsd(sarg, ShapeType::Face) as usize != nfac1 {
            return 0;
        }

        // cxx:1973-2004: no free wire / edge / vertex.
        for (find, avoid) in [
            (ShapeType::Wire, ShapeType::Face),
            (ShapeType::Edge, ShapeType::Wire),
            (ShapeType::Vertex, ShapeType::Edge),
        ] {
            if !explorer_shapes_avoiding(sarg, find, avoid).is_empty() {
                return 0;
            }
        }

        // cxx:2004
        1
    }

    /// `KPiskoletge()` (`KPart.cxx:1253-1385`): solids glued on tangent faces
    /// with tangent edges.
    pub fn kp_iskoletge(&mut self) -> i32 {
        // cxx:1259-1278
        let Some(shape1) = self.my_shape1.clone() else {
            return 0;
        };
        let mut lshsd1 = Vec::new();
        let mut lfhsd1 = Vec::new();
        if !self.kp_iskoletgesh(&shape1, &mut lshsd1, &mut lfhsd1) {
            return 0;
        }
        let (lplhsd1, lcyhsd1) = fun_sortplcy(&lfhsd1);
        if lplhsd1.len() != 1 || lcyhsd1.len() > 1 {
            return 0;
        }

        // cxx:1280-1292
        let Some(shape2) = self.my_shape2.clone() else {
            return 0;
        };
        let mut lshsd2 = Vec::new();
        let mut lfhsd2 = Vec::new();
        if !self.kp_iskoletgesh(&shape2, &mut lshsd2, &mut lfhsd2) {
            return 0;
        }
        let (lplhsd2, lcyhsd2) = fun_sortplcy(&lfhsd2);
        if lplhsd2.len() != 1 {
            return 0;
        }

        // cxx:1296-1302: several solids on either side fall back to the
        // general case.
        if lshsd1.len() > 1 || lshsd2.len() > 1 {
            return 0;
        }

        // cxx:1305-1313
        let f1 = lplhsd1[0].clone();
        let f2 = lplhsd2[0].clone();

        // cxx:1329-1337
        let Some(outerw1) = outer_wire_of(&f1) else {
            return 0;
        };
        let Some(outerw2) = outer_wire_of(&f2) else {
            return 0;
        };
        let mape1 = edge_set(&outerw1);
        let mape2 = edge_set(&outerw2);

        // cxx:1339-1347
        if !self.section_edges_on(&outerw1, &mape2) {
            return 0;
        }
        if !self.section_edges_on(&outerw2, &mape1) {
            return 0;
        }

        // cxx:1352-1359
        if !self.all_i_on_section_edges(&f1, &mape1, &mape2) {
            return 0;
        }
        if !self.all_i_on_section_edges(&f2, &mape2, &mape1) {
            return 0;
        }

        // cxx:1363-1383: KPcycy.
        if !lcyhsd1.is_empty() {
            if !(lcyhsd1.len() == 1 && lcyhsd2.len() == 1) {
                return 0;
            }
            if !fun_tool_in_s(&outerw1, &f1) {
                return 0;
            }
            if !fun_tool_in_s(&outerw2, &f2) {
                return 0;
            }
        }

        // cxx:1385
        1
    }

    /// `KPiskoletgesh(Sarg, lShsd, lfhsd)` (`KPart.cxx:1620-1707`).
    pub fn kp_iskoletgesh(
        &self,
        sarg: &TopoShape,
        lshsd: &mut Vec<TopoShape>,
        lfhsd: &mut Vec<TopoShape>,
    ) -> bool {
        // cxx:1628-1632
        if !self.fun_kp_iskolesh(sarg, lshsd) {
            return false;
        }

        // cxx:1637-1639
        self.kp_lhsd_into(sarg, ShapeType::Face, lfhsd);

        // cxx:1640-1677
        for fac in lfhsd.clone() {
            let Some(face) = Face::wrap(fac.clone()) else {
                return false;
            };
            let isplan = fun_tool_plane(&face);
            let iscylinder = fun_tool_cylinder(&face);
            if iscylinder {
                continue;
            }
            if !isplan {
                return false;
            }

            let Some(outerw) = brep_tools::outer_wire(&face) else {
                return false;
            };

            for ed in explorer_shapes(&outerw.0, ShapeType::Edge) {
                let isse = self
                    .my_data_structure
                    .as_ref()
                    .map(|hds| hds.ds().is_section_edge(&ed, true))
                    .unwrap_or(false);
                if self.splits(&ed, State::On).is_empty() {
                    return false;
                }
                if !isse {
                    return false;
                }
            }
        }

        // cxx:1707
        true
    }

    /// `FUNKP_KPiskolesh(BU, BDS, Sarg, lShsd, lfhsd)` (`kpkole.cxx:1033-1075`).
    ///
    /// `lfhsd` is unused in OCCT (`kpkole.cxx:1037`); `KPiskoletgesh` refills
    /// it through `KPlhsd` right after.
    pub(super) fn fun_kp_iskolesh(&self, sarg: &TopoShape, lshsd: &mut Vec<TopoShape>) -> bool {
        // kpkole.cxx:1044-1047
        if sarg.is_null() {
            return false;
        }

        // kpkole.cxx:1049-1053
        if self.kp_lhsd_into(sarg, ShapeType::Solid, lshsd) == 0 {
            return false;
        }
        let sol = lshsd[0].clone();

        // kpkole.cxx:1055-1072
        let mut lfhg = Vec::new();
        if self.kp_lhg_into(&sol, ShapeType::Face, &mut lfhg) != 0 {
            let Some(hds) = self.my_data_structure.as_ref() else {
                return false;
            };
            for f in lfhg {
                for i in hds.ds().shape_interferences(&f, true) {
                    if i.read()
                        .expect("interference lock")
                        .shape_shape_data()
                        .is_none()
                    {
                        return false;
                    }
                }
            }
        }

        // kpkole.cxx:1074
        true
    }

    /// `sectionedgesON(HDS, outerw1, mape2)` (`KPart.cxx:1195-1220`): every
    /// same-domain edge of every `outerw1` edge is in `mape2`.
    fn section_edges_on(&self, outerw1: &TopoShape, mape2: &HashSet<ShapeMapKey>) -> bool {
        let Some(hds) = self.my_data_structure.as_ref() else {
            return false;
        };
        for e1 in explorer_shapes(outerw1, ShapeType::Edge) {
            let sd = hds.same_domain(&e1);
            // cxx:1207-1209
            if sd.is_empty() {
                return false;
            }
            for e2 in sd {
                if !mape2.contains(&shape_map_key(e2)) {
                    return false;
                }
            }
        }
        true
    }

    /// `allIonsectionedges(HDS, f1, mape1, mape2)` (`KPart.cxx:1223-1246`).
    ///
    /// OCCT down-casts every interference of `f1` to
    /// `TopOpeBRepDS_ShapeShapeInterference` and dereferences it
    /// (`cxx:1234-1239`); a non-shape-shape interference would raise there. The
    /// port answers `false`, which is OCCT's answer for an interference whose
    /// geometry is outside both maps.
    fn all_i_on_section_edges(
        &self,
        f1: &TopoShape,
        mape1: &HashSet<ShapeMapKey>,
        mape2: &HashSet<ShapeMapKey>,
    ) -> bool {
        let Some(hds) = self.my_data_structure.as_ref() else {
            return false;
        };
        for i in hds.ds().shape_interferences(f1, true) {
            let i = i.read().expect("interference lock");
            let Some(ssd) = i.shape_shape_data() else {
                return false;
            };
            let g1 = i.geometry();
            let isgbound = ssd.g_bound;
            let Some(e1) = hds.shape(g1, true) else {
                return false;
            };
            let m = if isgbound { mape1 } else { mape2 };
            if !m.contains(&shape_map_key(&e1)) {
                return false;
            }
        }
        true
    }

    /// `KPiskoletgeanalyse(config2, Stsol1, Stsol2, ires)` (`KPart.cxx:2010-2096`).
    pub fn kp_iskoletge_analyse(&self, config2: Config, stsol1: State, stsol2: State) -> i32 {
        let mut ires = RESUNDEF;

        // cxx:2032-2058
        if config2 == Config::DiffOriented {
            if stsol1 == State::In && stsol2 == State::In {
                ires = RESNULL;
            }
            if stsol1 == State::Out && stsol2 == State::In {
                ires = RESSHAPE1;
            }
            if stsol2 == State::Out && stsol1 == State::In {
                ires = RESSHAPE2;
            }
            if stsol1 == State::Out && stsol2 == State::Out {
                ires = RESNEWSOL;
            }
        }

        // cxx:2060-2091
        if config2 == Config::SameOriented {
            if stsol1 == State::In && stsol2 == State::In {
                ires = RESSHAPE1;
            }
            if stsol1 == State::Out && stsol2 == State::In {
                ires = RESNULL;
            }
            if stsol2 == State::Out && stsol1 == State::In {
                ires = RESNEWSOL;
            }
            if stsol1 == State::Out && stsol2 == State::Out {
                ires = RESSHAPE2;
            }
        }

        ires
    }
}
