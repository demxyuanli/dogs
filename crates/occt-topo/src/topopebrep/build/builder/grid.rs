//! `TopOpeBRepBuild_Builder` — the `G*` helpers of `TopOpeBRepBuild_Grid.cxx`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Grid.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `GToMerge` | `95-120` |
//! | `GTakeCommonOfSame` | `125-147` |
//! | `GTakeCommonOfDiff` | `152-174` |
//! | `GFindSamDom(S, L1, L2)` | `182-190` |
//! | `GFindSamDomSODO(S, LSO, LDO)` | `260-268` |
//! | `GCopyList` (both overloads) | `441-464` |
//! | `GParamOnReference` | `1055-1116` |
//!
//! plus `TopOpeBRepBuild_FUN_aresamegeom` (`TopOpeBRepBuild_GridFF.cxx:132-141`),
//! which [`super::grid_ee`] needs.

use occt_core::elib::clib;
use occt_core::elib::slib;
use occt_core::gp::GpPnt2d;

use crate::abs::{Orientation, ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape, Vertex};
use crate::topo_tools::is_same;
use crate::topopebrep::ds::kind::Config;

use super::super::g_topo::GTopo;
use super::Builder;

impl Builder {
    /// `GToSplit(S, TB)` (`Grid.cxx:30-91`): `!IsSplit(S, TB) &&
    /// (S is a face of the DS || S has a geometry || S has same-domain)`.
    ///
    /// This is *not* `Builder::ToSplit` (`.cxx:308-326`, ported as
    /// [`Builder::to_split`]): the grid version adds the `facshap` term and
    /// asks `IsSplit` with the state `TB`.
    ///
    /// The `OCCT_DEBUG` block (`cxx:38-88`) only prints, and the
    /// `GetcontextSPEON` arm (`cxx:54-88`) is debug-only too.
    pub fn g_to_split(&self, s: &TopoShape, tb: State) -> bool {
        let issplit = self.is_split(s, tb);
        let facshap = self
            .data_structure()
            .map(|h| h.has_shape(s, true) && s.shape_type() == ShapeType::Face)
            .unwrap_or(false);
        let hasgeom = self
            .data_structure()
            .map(|h| h.has_geometry(s))
            .unwrap_or(false);
        let hassame = self
            .data_structure()
            .map(|h| h.has_same_domain(s, true))
            .unwrap_or(false);
        (!issplit) && (facshap || hasgeom || hassame)
    }

    /// `GToMerge(S)` (`Grid.cxx:95-120`): `!IsMerged(S, IN) && HasSameDomain(S)`.
    pub fn to_merge(&self, s: &TopoShape) -> bool {
        let ismerged = self.is_merged(s, State::In);
        let hassame = self
            .data_structure()
            .map(|h| h.has_same_domain(s, true))
            .unwrap_or(false);
        !ismerged && hassame
    }

    /// `GTakeCommonOfSame(G)` (`Grid.cxx:125-147`): both states are `ON` and
    /// the two are equal (`OUT`/`OUT` or `IN`/`IN`).
    pub fn take_common_of_same(g: &GTopo) -> bool {
        let (t1, t2) = g.states_on();
        matches!(
            (t1, t2),
            (State::Out, State::Out) | (State::In, State::In)
        )
    }

    /// `GTakeCommonOfDiff(G)` (`Grid.cxx:152-174`): the two states are `ON` and
    /// differ (`OUT`/`IN` or `IN`/`OUT`).
    pub fn take_common_of_diff(g: &GTopo) -> bool {
        let (t1, t2) = g.states_on();
        matches!(
            (t1, t2),
            (State::Out, State::In) | (State::In, State::Out)
        )
    }

    /// `GFindSamDom(S, L1, L2)` (`Grid.cxx:182-190`): clear both lists, seed
    /// `L1` with `S`, then [`Builder::find_same_domain`].
    pub fn find_same_domain_from(
        &self,
        s: &TopoShape,
        l1: &mut Vec<TopoShape>,
        l2: &mut Vec<TopoShape>,
    ) {
        l1.clear();
        l2.clear();
        l1.push(s.clone());
        self.find_same_domain(l1, l2);
    }

    /// `GFindSamDomSODO(S, LSO, LDO)` (`Grid.cxx:260-268`): clear both lists,
    /// seed `LSO` with `S`, then
    /// [`Builder::find_same_domain_same_orientation`].
    ///
    /// `GFindSamDomSODO` returns early (`cxx:277-280`) when `LSO` is empty after
    /// the seed, which cannot happen here; the guard is kept for symmetry.
    pub fn find_same_domain_same_orientation_from(
        &self,
        s: &TopoShape,
        lso: &mut Vec<TopoShape>,
        ldo: &mut Vec<TopoShape>,
    ) {
        lso.clear();
        ldo.clear();
        lso.push(s.clone());
        if lso.is_empty() {
            return;
        }
        let sref = lso[0].clone();
        let Some(hds) = self.data_structure() else {
            return;
        };
        let oref = hds.same_domain_orientation(&sref);
        self.find_same_domain(lso, ldo);

        let mut llso: Vec<TopoShape> = Vec::new();
        let mut lldo: Vec<TopoShape> = Vec::new();
        for s in lso.iter().chain(ldo.iter()) {
            let o = hds.same_domain_orientation(s);
            if o == oref && !Self::contains(s, &llso) {
                llso.push(s.clone());
            } else if o != oref && !Self::contains(s, &lldo) {
                lldo.push(s.clone());
            }
        }

        *lso = llso;
        *ldo = lldo;
    }

    /// `GCopyList(Lin, Lou)` (`Grid.cxx:459-464`): append the whole of `Lin`.
    pub fn copy_list(l_in: &[TopoShape], l_out: &mut Vec<TopoShape>) {
        Self::copy_list_range(l_in, 1, l_in.len() as i32, l_out);
    }

    /// `GCopyList(Lin, I1, I2, Lou)` (`Grid.cxx:441-456`): append the elements
    /// `[I1..I2]` of `Lin`, the first element of `Lin` being index 1.
    pub fn copy_list_range(
        l_in: &[TopoShape],
        i1: i32,
        i2: i32,
        l_out: &mut Vec<TopoShape>,
    ) {
        for (i, el) in l_in.iter().enumerate() {
            let i = i as i32 + 1;
            if i >= i1 && i <= i2 {
                l_out.push(el.clone());
            }
        }
    }

    /// `GParamOnReference(V, E, P)` (`Grid.cxx:1055-1116`): the parameter of `V`
    /// on `myEdgeReference` through its p-curve on `myFaceReference`.
    ///
    /// `cxx:1065-1076` rejects a non-planar reference face, `cxx:1082-1086`
    /// raises when the p-curve is missing, and `cxx:1108-1112` answers `false`
    /// for anything but the five conics. `cxx:1092-1108` switches on
    /// `Geom2dAdaptor_Curve::GetType()`:
    ///
    /// - `GeomAbs_Line` / `GeomAbs_Circle` use `ElCLib::Parameter` on the 2D
    ///   curve (`occt_core::elib::clib`).
    /// - UNPORTED: `GeomAbs_Ellipse` (`cxx:1098-1100`), `GeomAbs_Hyperbola`
    ///   (`cxx:1101-1103`) and `GeomAbs_Parabola` (`cxx:1104-1106`) need the 2D
    ///   `ElCLib::Parameter` overloads the port does not carry; they answer
    ///   `false` like OCCT's own `default` arm, leaving `p` untouched.
    pub fn param_on_reference(&self, v: &TopoShape, p: &mut f64) -> bool {
        let Some(face_ref) = self.my_face_reference.clone() else {
            return false;
        };
        let Some(edge_ref) = self.my_edge_reference.clone() else {
            return false;
        };

        let Some(su) = BRepTool::face_surface(&Face(face_ref.clone())) else {
            return false;
        };
        let Some(pln) = su.gp_pln() else {
            return false;
        };

        let p3 = BRepTool::vertex_point(&Vertex(v.clone()));
        let (u, vv) = slib::plane_parameters(&pln.position(), &p3);
        let p2 = GpPnt2d::new(u, vv);

        let Some((c2d, _f, _l)) = crate::boptools_2d::curve_on_surface_range(
            &Edge(edge_ref.clone()),
            &Face(face_ref.clone()),
        ) else {
            return false;
        };

        if let Some(l) = c2d.gp_lin2d() {
            *p = clib::line2d_parameter(&l, &p2);
            return true;
        }
        if let Some(c) = c2d.gp_circ2d() {
            *p = clib::circle2d_parameter(&c.position(), &p2);
            return true;
        }

        false
    }
}

/// `TopOpeBRepBuild_FUN_aresamegeom(S1, S2)` (`GridFF.cxx:132-141`): both
/// shapes FORWARD, then `TopOpeBRepTool_ShapeTool::ShapesSameOriented`.
pub fn fun_aresamegeom(s1: &TopoShape, s2: &TopoShape) -> bool {
    let mut sf1 = s1.clone();
    sf1.set_orientation(Orientation::Forward);
    let mut sf2 = s2.clone();
    sf2.set_orientation(Orientation::Forward);
    crate::topopebrep::tool::shape_tool::shapes_same_oriented(&sf1, &sf2)
}

/// `myDataStructure->DS().ShapeSameDomain(V).First()`, `None` when the list is
/// empty (OCCT dereferences the iterator, so it would raise there).
pub(super) fn first_same_domain(hds: &crate::topopebrep::ds::h_data_structure::HDataStructure,
                                v: &TopoShape) -> Option<TopoShape> {
    hds.same_domain(v).first().cloned()
}

/// `G1.Config1() == TopOpeBRepDS_DIFFORIENTED`.
pub(super) fn is_diff_oriented(g: &GTopo) -> bool {
    g.config1() == Config::DiffOriented
}
