//! `TopOpeBRep_FacesFiller` — fills the DS with the result of a face/face
//! intersection.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRep/TopOpeBRep_FacesFiller.hxx:42-311`,
//! `TopOpeBRep_FacesFiller.cxx:19-750`, `TopOpeBRep_FacesFiller_1.cxx:19-1508`,
//! `TopOpeBRep_ProcessGR.cxx:57-443`, `TopOpeBRep_ProcessSectionEdges.cxx:19-468`,
//! `TopOpeBRep_vpr.cxx`, `TopOpeBRep_vprclo.cxx`, `TopOpeBRep_vprdeg.cxx`.
//!
//! Translation notes:
//! - `myFacesIntersector`, `myHDS`, `myDS` and `myLine` (`hxx:290-295`) are raw
//!   pointers into the caller's intersector / data structure. OCCT keeps them
//!   across calls; the port passes `&mut FacesIntersector` and `&mut
//!   HDataStructure` to each method and keeps only the *line index*
//!   ([`FacesFiller::line_index`]). The line itself is reached through
//!   `FacesIntersector::change_line`, which is what `*myLine` denotes.
//! - `myPShapeClassifier` (`hxx:308`) is a `TopOpeBRepTool_PShapeClassifier`
//!   pointer; `nullptr` (the constructor, `cxx:90-97`) is `None` here.
//!   `Insert` allocates a `TopOpeBRepTool_ShapeClassifier` when it is unset
//!   (`cxx:128-131`).
//! - `myPointClassifier` (`hxx:307`) is a value member, stored by value here.
//! - `myDataforDegenEd` (`hxx:301-302`) is an `NCollection_DataMap<TopoDS_Shape,
//!   List<TopoDS_Shape>>`; the port stores an association list
//!   ([`DegenEdgeData`]) because `TopoShape` has no hash.
//! - Debug-only members and calls (`myHFFD`, `TopOpeBRep_FFDumper`, `debrest`,
//!   `debvprmess`, ...) are not translated.
//!
//! Ported in this file: the constructor, the classifier accessors, `Insert`,
//! `LoadLine`, `CheckLine`, the three `VP_Position` overloads, `PequalVPonR`,
//! `Lminmax`, `IsVPtransLok`, `TransvpOK` and the accessors.
//!
//! Ported in the sibling modules:
//! - `section.rs`: `ProcessSectionEdges` / `GetESL` / `KeepRLine` /
//!   `FUN_BREP_SDM_RE` (`TopOpeBRep_ProcessSectionEdges.cxx:19-468`),
//! - `line.rs`: `ProcessLine` / `FillLine` / `FillLineVPonR` / `AddShapesLine` /
//!   `GetGeometry` / `MakeGeometry` / `GetFFGeometry`
//!   (`TopOpeBRep_FacesFiller_1.cxx:460-1507`),
//! - `process_gr.rs`: `StBipVPonF` / `StateVPonFace` / `LSameDomainERL` /
//!   `VPParamOnER` / `EqualpPonR` (`TopOpeBRep_ProcessGR.cxx:59-444`),
//! - `vpr.rs`: `ProcessVPR` / `ProcessVPnotonR` / `FUN_VPIndex` /
//!   `FUN_VPgeometryfound` / `FUN_transForWL` (`TopOpeBRep_vpr.cxx:120-470`,
//!   `_1.cxx:507-645`),
//! - `vpclo.rs`: `GetEdgeTrans` / `ProcessVPonclosingR`
//!   (`TopOpeBRep_vprclo.cxx:66-296`),
//! - `vponr.rs`: the `ProcessVPonR` entry point plus `FUN_onedge` /
//!   `FUN_brep_ONfirstP`,
//! - `degen_data.rs`: `FUN_FillVof12` / `FUN_addmapve` / `FUN_GetdgData`
//!   (`TopOpeBRep_vprdeg.cxx:105-379`),
//! - `vprdeg_interf.rs`: `FUN_putInterfonDegenEd` / `ProcessVPondgE` and the
//!   `local_*` helpers around them (`vprdeg.cxx:381-862`),
//! - `vp_index.rs`: `FUN_VPIndex` (`TopOpeBRep_vpr.cxx:350-469`),
//! - `vponr.rs` / `vponr_ds.rs`: the two halves of `ProcessVPonR`
//!   (`TopOpeBRep_vpr.cxx:877-1732`),
//! - `vpr_trans_edge.rs`: `FUN_LineRestF` / `FUN_newtransEdge` and the
//!   `FUN_select*interference` helpers (`vpr.cxx:471-682`,
//!   `TopOpeBRepDS_ProcessInterferencesTool.cxx:111-226`),
//! - `vpr_cpi.rs`: `FUN_processCPI` (`vpr.cxx:782-853`),
//! - `rline.rs`: `ProcessRLine` (`_1.cxx:647-1125`).
//!
//! UNPORTED (explicitly, with the OCCT body to translate):
//! - `FC2D_HasC3D` / `FUN_tool_projPonC2D` (`TopOpeBRepTool_PROJECT.cxx:146-224`)
//!   in `fun_proj_pon_l`.
//! - The `processUNK` block of `ProcessRLine` (`_1.cxx:687-745`): compiled out
//!   in OCCT (`processUNK` is `false` at `_1.cxx:682`).
//!
//! `LineInter::has_v_pon_r()` is true for most real face/face sections, so
//! [`FacesFiller::insert`] runs `FillLineVPonR` → `ProcessRLine` /
//! `ProcessVPonR` on such a line.

use occt_core::gp::GpPnt;
use occt_core::precision::PCONFUSION;

use crate::abs::Orientation;
use crate::shape::{Face, TopoShape};
use crate::topopebrep::ds::h_data_structure::HDataStructure;
use crate::topopebrep::ds::transition::Transition;
use crate::topopebrep::faces_intersector::FacesIntersector;
use crate::topopebrep::ff_transition_tool::FFTransitionTool;
use crate::topopebrep::interf::v_point_inter::VPointInter;
use crate::topopebrep::interf::v_point_inter_classifier::VPointInterClassifier;
use crate::topopebrep::interf::PointClassifier;
use crate::topopebrep::line_inter::LineInter;
use crate::topopebrep::tool::shape_classifier::ShapeClassifier;
use crate::topopebrep::type_line_curve::TypeLineCurve;

pub mod degen_data;
pub mod line;
pub mod process_gr;
pub mod rline;
pub mod section;
pub mod vp_index;
pub mod vp_position;
pub mod vpclo;
pub mod vponr;
pub mod vponr_ds;
pub mod vpr;
pub mod vpr_cpi;
pub mod vpr_trans_edge;
pub mod vprdeg_interf;

/// `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
/// TopTools_ShapeMapHasher> myDataforDegenEd` (`hxx:301-302`).
///
/// The key is the face the degenerated data was collected for; the value is
/// the list `FUN_GetdgData` / `FUN_FillVof12` build
/// (`TopOpeBRep_ProcessSectionEdges.cxx`).
#[derive(Default)]
pub struct DegenEdgeData {
    pairs: Vec<(TopoShape, Vec<TopoShape>)>,
}

impl DegenEdgeData {
    /// `Clear()` (`TopOpeBRep_FacesFiller_1.cxx:1165`, commented out in OCCT).
    pub fn clear(&mut self) {
        self.pairs.clear();
    }

    /// `Bind(K, I)`.
    pub fn bind(&mut self, k: &TopoShape, v: &[TopoShape]) {
        match self.pairs.iter_mut().find(|(k0, _)| crate::topo_tools::is_same(k0, k)) {
            Some((_, dst)) => {
                dst.clear();
                dst.extend(v.iter().cloned());
            }
            None => self.pairs.push((k.clone(), v.to_vec())),
        }
    }

    /// `IsBound(K)` / `Find(K)`.
    pub fn find(&self, k: &TopoShape) -> Option<&[TopoShape]> {
        self.pairs
            .iter()
            .find(|(k0, _)| crate::topo_tools::is_same(k0, k))
            .map(|(_, v)| v.as_slice())
    }

    /// `IsBound(K)`.
    pub fn is_bound(&self, k: &TopoShape) -> bool {
        self.find(k).is_some()
    }

    /// `ChangeFind(K).Append(V)` (`FUN_addmapve`, `vprdeg.cxx:148`, `:154`,
    /// `:172-173`): appends to the existing list, binding `K` first when absent.
    pub fn change_append(&mut self, k: &TopoShape, v: &TopoShape) {
        match self.pairs.iter_mut().find(|(k0, _)| crate::topo_tools::is_same(k0, k)) {
            Some((_, dst)) => dst.push(v.clone()),
            None => self.pairs.push((k.clone(), vec![v.clone()])),
        }
    }

    /// `IsEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// `Iterator` over the `(key, list)` bindings, in insertion order (the
    /// `NCollection_DataMap::Iterator` order `FUN_GetdgData` relies on).
    pub fn iter(&self) -> impl Iterator<Item = (&TopoShape, &Vec<TopoShape>)> {
        self.pairs.iter().map(|(k, v)| (k, v))
    }
}

/// `class TopOpeBRep_FacesFiller` (`hxx:42-311`).
pub struct FacesFiller {
    /// `myF1` (`hxx:286`).
    f1: Option<Face>,
    /// `myF2` (`hxx:287`).
    f2: Option<Face>,
    /// `myF1ori` (`hxx:288`).
    f1_ori: Orientation,
    /// `myF2ori` (`hxx:289`).
    f2_ori: Orientation,
    /// `myFFfirstDSP` (`hxx:294`).
    ff_first_dsp: i32,
    /// `myLine` (`hxx:295`), as an index into the `FacesIntersector`.
    line_index: i32,
    /// `myLineOK` (`hxx:296`).
    line_ok: bool,
    /// `myLineINL` (`hxx:297`).
    line_inl: bool,
    /// `myLineIsonEdge` (`hxx:299`).
    line_is_on_edge: bool,
    /// `myERL` (`hxx:300`).
    erl: Vec<TopoShape>,
    /// `myDataforDegenEd` (`hxx:301-302`).
    data_for_degen_ed: DegenEdgeData,
    /// `myLastVPison0` (`hxx:303`).
    last_vp_is_on0: bool,
    /// `mykeptVPnbr` (`hxx:304`).
    kept_vp_nbr: i32,
    /// `myDSCIndex` (`hxx:305`).
    dsc_index: i32,
    /// `myDSCIL` (`hxx:306`).
    dscil: crate::topopebrep::ds::interference_iterator::InterferenceList,
    /// `myPointClassifier` (`hxx:307`).
    point_classifier: PointClassifier,
    /// `myPShapeClassifier` (`hxx:308`). OCCT holds a raw pointer shared with
    /// `TopOpeBRep_DSFiller`; the port shares it through
    /// `Rc<RefCell<..>>`.
    pshape_classifier: Option<std::rc::Rc<std::cell::RefCell<ShapeClassifier>>>,
    /// `myexF1` (`hxx:309`).
    ex_f1: i32,
    /// `myexF2` (`hxx:310`).
    ex_f2: i32,
}

impl Default for FacesFiller {
    fn default() -> Self {
        Self::new()
    }
}

impl FacesFiller {
    /// `TopOpeBRep_FacesFiller()` (`cxx:90-97`): `myPShapeClassifier` is null
    /// and the two trace indices are zero.
    pub fn new() -> Self {
        Self {
            f1: None,
            f2: None,
            f1_ori: Orientation::Forward,
            f2_ori: Orientation::Forward,
            ff_first_dsp: 0,
            line_index: 0,
            line_ok: false,
            line_inl: false,
            line_is_on_edge: false,
            erl: Vec::new(),
            data_for_degen_ed: DegenEdgeData::default(),
            last_vp_is_on0: false,
            kept_vp_nbr: 0,
            dsc_index: 0,
            dscil: Vec::new(),
            point_classifier: PointClassifier::new(),
            pshape_classifier: None,
            ex_f1: 0,
            ex_f2: 0,
        }
    }

    /// `PShapeClassifier()` (`cxx:102-105`).
    pub fn pshape_classifier(&self) -> Option<std::rc::Rc<std::cell::RefCell<ShapeClassifier>>> {
        self.pshape_classifier.clone()
    }

    /// `SetPShapeClassifier(PSC)` (`cxx:109-112`).
    pub fn set_pshape_classifier(
        &mut self,
        psc: Option<std::rc::Rc<std::cell::RefCell<ShapeClassifier>>>,
    ) {
        self.pshape_classifier = psc;
    }

    /// `ChangePointClassifier()` (`cxx:184-187`).
    pub fn change_point_classifier(&mut self) -> &mut PointClassifier {
        &mut self.point_classifier
    }

    /// `Insert(S1, S2, FACINT, HDS)` (`cxx:116-179`).
    ///
    /// The same-domain arm (`cxx:144-148`) records the two faces and returns.
    /// The rest walks the intersector's lines three times, exactly as OCCT
    /// does, then calls `ProcessSectionEdges`, `FUN_MakeERL` and the
    /// `LoadLine` / `ProcessLine` pass.
    pub fn insert(
        &mut self,
        s1: &TopoShape,
        s2: &TopoShape,
        facint: &mut FacesIntersector,
        hds: &mut HDataStructure,
    ) -> Result<(), String> {
        self.f1 = Face::wrap(s1.clone());
        self.f1_ori = s1.orientation();
        self.f2 = Face::wrap(s2.clone());
        self.f2_ori = s2.orientation();
        if self.pshape_classifier.is_none() {
            self.pshape_classifier = Some(std::rc::Rc::new(std::cell::RefCell::new(
                ShapeClassifier::new(),
            )));
        }

        if facint.same_domain()? {
            hds.change_ds().fill_shapes_same_domain(s1, s2, false)?;
            return Ok(());
        }

        // `cxx:149-154`: set the faces on every line.
        facint.init_line();
        while facint.more_line() {
            let i = facint.current_line_index();
            let (Some(f1), Some(f2)) = (self.f1.clone(), self.f2.clone()) else {
                return Err("TopOpeBRep_FacesFiller::Insert : null face".to_string());
            };
            facint.change_line(i).set_faces(&f1, &f2);
            facint.next_line();
        }

        // `cxx:156`
        self.vp_position_all(facint)?;

        // `cxx:158-163`
        facint.init_line();
        while facint.more_line() {
            let i = facint.current_line_index();
            let l = facint.change_line(i);
            l.set_has_v_pon_r();
            l.set_inl();
            l.set_is_v_closed();
            facint.next_line();
        }

        // `cxx:166` ProcessSectionEdges();
        self.process_section_edges(facint, hds)?;

        // `cxx:167` myFFfirstDSP = myDS->NbPoints() + 1;
        self.ff_first_dsp = hds.change_ds().nb_points() + 1;

        // `cxx:169` FUN_MakeERL((*myFacesIntersector), myERL);
        self.erl.clear();
        self.erl.extend(facint.restrictions().iter().cloned());

        // `cxx:171-178`
        facint.init_line();
        while facint.more_line() {
            let i = facint.current_line_index();
            self.load_line(facint, i)?;
            self.process_line(facint, hds)?;
            facint.next_line();
        }

        Ok(())
    }

    /// `LoadLine(L)` (`cxx:190-211`).
    pub fn load_line(
        &mut self,
        facint: &mut FacesIntersector,
        index: i32,
    ) -> Result<(), String> {
        self.line_index = index;
        let mut bchk = self.check_line(facint, index)?;
        let binl = facint.change_line(index).inl();
        self.line_inl = binl;
        let t = facint.change_line(index).type_line_curve();
        if !bchk && binl && t == TypeLineCurve::Line {
            bchk = true;
        }
        facint.change_line(index).set_ok(bchk);
        self.line_ok = bchk;
        if !self.line_ok {
            return Ok(());
        }
        facint.change_line(index).compute_face_face_transition();
        Ok(())
    }

    /// `CheckLine(L)` (`cxx:219-371`): `false` for a walking line with fewer
    /// than two VPoints, a two-point line whose two ends are the same vertex,
    /// and the degenerate conic cases.
    ///
    /// OCCT reads `myLine->VPBounds` in the `CIRCLE` / `ELLIPSE` arms
    /// (`cxx:315`, `:340`); its `myLine` is `L` by construction of `LoadLine`,
    /// so the port reads `l` there.
    pub fn check_line(
        &self,
        facint: &mut FacesIntersector,
        index: i32,
    ) -> Result<bool, String> {
        let (tol1, _tol2) = facint.get_tolerances();
        let l = facint.change_line(index);
        let mut check = true;
        let t = l.type_line_curve();
        let nbvp = l.nb_v_point();

        if t == TypeLineCurve::Walking {
            if nbvp < 2 {
                check = false;
            }
        } else if t == TypeLineCurve::Line {
            let mut np = 0;
            for k in 1..=nbvp {
                if l.v_point(k).map(|vp| vp.keep()).unwrap_or(false) {
                    np += 1;
                }
            }
            if np != 2 {
                return Ok(true);
            }

            let mut a: Option<&VPointInter> = None;
            let mut b: Option<&VPointInter> = None;
            for k in 1..=nbvp {
                let Some(vp) = l.v_point(k) else { continue };
                if !vp.keep() {
                    continue;
                }
                if a.is_none() {
                    a = Some(vp);
                } else {
                    b = Some(vp);
                }
            }
            let (Some(a), Some(b)) = (a, b) else {
                return Ok(check);
            };

            let mut v1: Option<TopoShape> = None;
            if a.is_vertex_on_s1().unwrap_or(false) {
                v1 = a.vertex_on_s1().ok().cloned();
            }
            if a.is_vertex_on_s2().unwrap_or(false) {
                v1 = a.vertex_on_s2().ok().cloned();
            }
            let mut v2: Option<TopoShape> = None;
            if b.is_vertex_on_s1().unwrap_or(false) {
                v2 = b.vertex_on_s1().ok().cloned();
            }
            if b.is_vertex_on_s2().unwrap_or(false) {
                v2 = b.vertex_on_s2().ok().cloned();
            }

            if let (Some(v1), Some(v2)) = (v1, v2) {
                if crate::topo_tools::is_same(&v1, &v2) {
                    return Ok(false);
                }
            }
        } else {
            let notrnotw =
                t != TypeLineCurve::Restriction && t != TypeLineCurve::Walking;
            if notrnotw {
                if t == TypeLineCurve::Circle {
                    // cto 012 D2, faces 6 et 1, line 3 incorrecte.
                    if let Some((i_inon1, i_inonn, n_inon)) = l.vp_bounds() {
                        if n_inon >= 2 {
                            let a_par = l
                                .v_point(i_inon1)
                                .and_then(|vp| vp.parameter_on_line().ok());
                            let b_par = l
                                .v_point(i_inonn)
                                .and_then(|vp| vp.parameter_on_line().ok());
                            if let (Some(pa), Some(pb)) = (a_par, b_par) {
                                if (pa - pb).abs() < tol1 {
                                    check = false;
                                }
                            }
                        }
                    }
                } else if t == TypeLineCurve::Hyperbola {
                    if let Some((_f, _ll, n_inon)) = l.vp_bounds() {
                        if n_inon < 2 {
                            check = false;
                        }
                    }
                } else if t == TypeLineCurve::Ellipse {
                    if let Some((i_inon1, i_inonn, n_inon)) = l.vp_bounds() {
                        if n_inon < 2 {
                            check = false;
                        } else {
                            let a_par = l
                                .v_point(i_inon1)
                                .and_then(|vp| vp.parameter_on_line().ok());
                            let b_par = l
                                .v_point(i_inonn)
                                .and_then(|vp| vp.parameter_on_line().ok());
                            if let (Some(pa), Some(pb)) = (a_par, b_par) {
                                if (pa - pb).abs() < tol1 {
                                    check = false;
                                }
                            }
                        }
                    }
                }
            }
        }

        Ok(check)
    }

    /// `Face(I)` (`cxx:690-702`): `Standard_ProgramError` for any other index.
    pub fn face(&self, i: i32) -> Result<&Face, String> {
        match i {
            1 => self
                .f1
                .as_ref()
                .ok_or_else(|| "TopOpeBRep_FacesFiller::Face".to_string()),
            2 => self
                .f2
                .as_ref()
                .ok_or_else(|| "TopOpeBRep_FacesFiller::Face".to_string()),
            _ => Err("TopOpeBRep_FacesFiller::Face".to_string()),
        }
    }

    /// `FaceFaceTransition(L, I)` (`cxx:706-711`).
    pub fn face_face_transition_of<'a>(
        &self,
        l: &'a LineInter,
        i: i32,
    ) -> Result<&'a Transition, String> {
        l.face_face_transition(i)
    }

    /// `FaceFaceTransition(I)` (`cxx:715-719`): the transition of the current
    /// line.
    pub fn face_face_transition<'a>(
        &self,
        facint: &'a FacesIntersector,
        i: i32,
    ) -> Result<&'a Transition, String> {
        let l = facint
            .lines()
            .and_then(|h| h.value(self.line_index))
            .ok_or_else(|| "TopOpeBRep_FacesFiller::FaceFaceTransition".to_string())?;
        l.face_face_transition(i)
    }

    /// `myLine->VPBounds(iINON1, iINONn, nINON)`: the two in/on bounds of the
    /// current line and how many in/on VPoints it has.
    ///
    /// `LineInter::vp_bounds` needs `&mut` because OCCT caches the two indices
    /// on first use; the port clones the line and caches on that clone, which is
    /// the same value the next call would compute.
    pub fn line_vp_bounds(&self, facint: &mut FacesIntersector) -> (i32, i32, i32) {
        let mut l = facint.change_line(self.line_index).clone();
        l.vp_bounds().unwrap_or((0, 0, 0))
    }

    /// `SetTraceIndex(exF1, exF2)` (`cxx:737-741`).
    pub fn set_trace_index(&mut self, ex_f1: i32, ex_f2: i32) {
        self.ex_f1 = ex_f1;
        self.ex_f2 = ex_f2;
    }

    /// `GetTraceIndex(exF1, exF2)` (`cxx:745-750`).
    pub fn get_trace_index(&self) -> (i32, i32) {
        (self.ex_f1, self.ex_f2)
    }
}