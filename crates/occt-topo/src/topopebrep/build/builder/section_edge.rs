//! `TopOpeBRepBuild_Builder` — `SplitSectionEdges`, the section-edge list
//! splitter of `TopOpeBRepBuild_Section.cxx`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Section.cxx:304-914`
//! (OCCT 8.0.0), plus the file-scope helpers it calls:
//!
//! | OCCT | lines |
//! |---|---|
//! | `FUN_periodicS` | `85-114` |
//! | `FUN_periodic` | `116-123` |
//! | `FUN_onboundsper(uv, F)` | `125-156` |
//! | `FUN_onboundsper(E, F, uv)` | `158-176` |
//! | `FUN_PinC` | `178-197` |
//! | `FUN_PinE` | `200-209` |
//! | `FUN_EstaEE` | `211-244` |
//! | `FUN_selLEE` | `257-277` |
//! | `FUN_FindEinSLOS` | `282-300` |
//!
//! Six passes, in order: split every section edge
//! ([`Builder::split_section_edge`]); collect the shared `Split(ON)` pieces of
//! same-domain section edges (the `MEIN` / `MEOUT` maps); rewrite the `Split(ON)`
//! lists of the edges that took part; give every `Split(ON)` piece of a
//! reference section edge the p-curve of every connected face; then two passes
//! over periodic faces — the first materialises missing p-curves, the second
//! translates the iso-parametric pieces back onto the face bounds.
//!
//! Deviations, all forced by the port's data model:
//! - `FDSCNX_EdgeConnexitySameShape` (`cxx:619`) reads OCCT's process-wide
//!   `GLOBAL_elf1` / `GLOBAL_elf2` / `GLOBAL_fle` tables, prepared by the DS
//!   filler. The port has no process-wide tables, so pass 4 rebuilds a
//!   [`Connexity`] from the same DS. `FDSCNX_Prepare`
//!   (`TopOpeBRepDS_connex.cxx:127-220`) derives both tables from the DS alone,
//!   so the contents match.
//! - `FUN_onboundsper(E, F, uv)` answers an out-parameter even when it returns
//!   `false`; pass 6 reads it only after a `false`. OCCT reads uninitialised
//!   stack memory there when the 3-D projection itself failed, so the port
//!   answers `None` and skips the translation (equivalent to a zero offset).
//! - UNPORTED: `FUN_onboundsper(E, F, uv)`'s 3-D projection out-parameter when
//!   the projection itself failed — see the previous bullet; every p-curve step
//!   of pass 4 runs the seven-argument `FC2D_CurveOnSurface`
//!   (`TopOpeBRepTool_2d.cxx:450-566`) with the `EF` argument OCCT looks up at
//!   `cxx:636-662`.

use std::collections::HashMap;
use std::sync::Arc;

use occt_core::gp::{GpPnt, GpPnt2d, GpVec2d};
use occt_geom2d::curve::Curve2d;

use crate::abs::{ShapeType, State};
use crate::brep_tool::BRepTool;
use crate::geom_bnd_lib_surface3d::{surface_kind, GeomAbsSurfaceType};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{explorer_shapes, is_equal, is_same, shape_map_key, ShapeMapKey};
use crate::topopebrep::build::ShapeListOfShape;
use crate::topopebrep::ds::connex::Connexity;
use crate::topopebrep::ds::export::fds_data;
use crate::topopebrep::ds::kind::Kind;
use crate::topopebrep::tool::c2d::fc2d_curve_on_surface_with_ef;
use crate::topopebrep::tool::iso::uv_iso_curve;
use crate::topopebrep::tool::project_curve::fun_tool_bounds;
use crate::topopebrep::tool::split_e::{split_e, trsl_uv_modif_e};
use crate::topopebrep::tool::topology::fun_tool_par_f;
use crate::topopebrep::build::tools::write_edge_pcurve;

use super::Builder;

/// `FUN_periodicS(F)` (`cxx:85-114`): a cylinder, a cone or a sphere.
///
/// OCCT goes through `TopOpeBRepTool_ShapeTool::BASISSURFACE` then
/// `GeomAdaptor_Surface::GetType`, and its `NYI : for Torus,
/// SurfaceOfRevolution..` comment keeps the torus out.
fn fun_periodic_s(f: &TopoShape) -> bool {
    if f.shape_type() != ShapeType::Face {
        return false;
    }
    let Some(su) = BRepTool::face_surface(&Face(f.clone())) else {
        return false;
    };
    matches!(
        surface_kind(su.as_ref()),
        GeomAbsSurfaceType::Cylinder | GeomAbsSurfaceType::Cone | GeomAbsSurfaceType::Sphere
    )
}

/// `FUN_periodic(F, uper, vper)` (`cxx:116-123`).
///
/// `per = uper || vper` is left to the caller, as in OCCT.
fn fun_periodic(f: &Face) -> Option<(bool, bool)> {
    let su = BRepTool::face_surface(f)?;
    Some((su.is_u_periodic(), su.is_v_periodic()))
}

/// `FUN_onboundsper(uv, F)` (`cxx:125-156`, the 2-D arm): `uv` sits on a bound
/// of a periodic direction, within `1.e-8 * 1.e-2`.
fn fun_onboundsper_uv(uv: &GpPnt2d, f: &TopoShape) -> bool {
    let Some(su) = BRepTool::face_surface(&Face(f.clone())) else {
        return false;
    };
    let uclo = su.is_u_periodic();
    let vclo = su.is_v_periodic();
    if !uclo && !vclo {
        return false;
    }
    let (u1, u2) = su.u_range();
    let (v1, v2) = su.v_range();
    let toluv = 1.0e-8 * 1.0e-2;
    if uclo {
        return (u1 - uv.x()).abs() < toluv || (u2 - uv.x()).abs() < toluv;
    }
    if vclo {
        return (v1 - uv.y()).abs() < toluv || (v2 - uv.y()).abs() < toluv;
    }
    false
}

/// `FUN_onboundsper(E, F, uv)` (`cxx:158-176`, the 3-D arm): projects the point
/// at `0.45678 * f + (1 - 0.45678) * l` of `e` on `f`, then asks the 2-D arm.
///
/// Answers `(on_bounds, uv)`, `uv` being the projection when it could be
/// computed — OCCT keeps it even when the answer is `false`.
fn fun_onboundsper_edge(e: &TopoShape, f: &TopoShape) -> (bool, Option<GpPnt2d>) {
    let (ff, ll) = fun_tool_bounds(e);
    let x = 0.45678;
    let par = ff * x + (1.0 - x) * ll;
    let Some(uv) = fun_tool_par_f(e, par, f) else {
        return (false, None);
    };
    let on = fun_onboundsper_uv(&uv, f);
    (on, Some(uv))
}

/// `FUN_PinC(P, C, pmin, pmax, tol)` (`cxx:178-197`): `C` projects `P` to
/// within `tol`.
fn fun_pin_c(p: &GpPnt, c: &dyn occt_geom::curve::Curve, tol: f64) -> bool {
    match occt_geom::geom_api::project_point_on_curve(c, p, tol) {
        Some(r) => r.distance <= tol,
        None => false,
    }
}

/// `FUN_PinE(P, E)` (`cxx:200-209`): [`fun_pin_c`] with the edge tolerance.
fn fun_pin_e(p: &GpPnt, e: &TopoShape) -> bool {
    let edge = Edge(e.clone());
    let Some(c) = BRepTool::edge_curve(&edge) else {
        return false;
    };
    let tol = BRepTool::edge_tolerance(&edge);
    fun_pin_c(p, c.as_ref(), tol)
}

/// `FUN_EstaEE(E, sta, EE)` (`cxx:211-244`): the point of `e` at
/// `0.417789 * f + (1 - 0.417789) * l` is in `ee` (or is not, for `OUT`).
///
/// A degenerated edge answers `sta != IN` without evaluating anything. OCCT
/// raises `Standard_ProgramError` for any other state; the port answers
/// `false`, which is what the `IN`/`OUT` callers can observe.
fn fun_esta_ee(e: &TopoShape, sta: State, ee: &TopoShape) -> bool {
    let edge = Edge(e.clone());
    if BRepTool::is_degenerated(&edge) {
        return sta != State::In;
    }
    let Some(c) = BRepTool::edge_curve(&edge) else {
        return false;
    };
    let (f, l) = BRepTool::edge_parameters(&edge);
    let t = 0.417789;
    let par = (1.0 - t) * f + t * l;
    let p = c.d0(par);
    let ein_ee = fun_pin_e(&p, ee);
    match sta {
        State::In => ein_ee,
        State::Out => !ein_ee,
        _ => false,
    }
}

/// `FUN_selLEE(LE, E, sta, LEsta)` (`cxx:257-277`): move the edges of `le`
/// located `sta` / `e` into `lesta` (which is appended to, never cleared).
fn fun_sel_lee(le: &mut Vec<TopoShape>, e: &TopoShape, sta: State, lesta: &mut Vec<TopoShape>) {
    let mut i = 0;
    while i < le.len() {
        let cur = le[i].clone();
        if fun_esta_ee(&cur, sta, e) {
            lesta.push(cur);
            le.remove(i);
        } else {
            i += 1;
        }
    }
}

/// `FUN_FindEinSLOS(E, loslos)` (`cxx:282-300`): one cell of `loslos` carries
/// `e` as its shape.
fn fun_find_e_in_slos(e: &TopoShape, loslos: &[ShapeListOfShape]) -> bool {
    loslos
        .iter()
        .any(|slos| slos.shape().is_some_and(|s| is_same(s, e)))
}

/// `FC2D_CurveOnSurface(E, F, f, l, tolpc, trim3d)` (`2d.cxx:356-376`): the
/// stored p-curve (`FC2D_HasOldCurveOnSurface` / `FC2D_HasNewCurveOnSurface`),
/// otherwise `FC2D_MakeCurveOnSurface` (`2d.cxx:342-353`), which builds the
/// curve and registers it on the edge.
pub(super) fn fc2d_curve_on_surface(e: &TopoShape, f: &TopoShape) -> Option<Arc<dyn Curve2d>> {
    crate::boptools_2d::build_pcurve_for_edge_on_face(&Edge(e.clone()), &Face(f.clone())).ok()
}

impl Builder {
    /// `SplitSectionEdges()` (`Section.cxx:304-915`).
    pub fn split_section_edges(&mut self) {
        if self.my_split_section_edges_done {
            return;
        }

        // `cxx:311-321`: the non-null section edges, in DS order.
        let nes = match self.data_structure() {
            Some(hds) => hds.ds().nb_section_edges(),
            None => 0,
        };
        let mut sedges: Vec<TopoShape> = Vec::with_capacity(nes.max(0) as usize);
        for i in 1..=nes {
            let e = self
                .data_structure()
                .and_then(|hds| hds.ds().section_edge(i, true));
            if let Some(e) = e {
                sedges.push(e);
            }
        }

        for e in &sedges {
            self.split_section_edge(e);
        }

        let (mein, meout) = self.collect_shared_on_pieces(&sedges);
        self.rewrite_on_lists(&sedges, &mein, &meout);
        self.attach_pcurves_on_connected_faces(&sedges);
        self.materialise_periodic_pcurves();
        self.translate_on_bounds_pieces();

        self.my_split_section_edges_done = true;
    }

    /// `cxx:326-514`, pass 2: walk the same-domain partners of every split
    /// section edge, remove from `LEoutLESD` the pieces that lie `IN` a partner
    /// and record, per pair `(E1, E2)`, the pieces they share.
    ///
    /// `MEIN` maps a section edge to `(other edge, shared pieces)` cells;
    /// `MEOUT` maps it to the pieces that stay `OUT` of every partner.
    #[allow(clippy::type_complexity)]
    fn collect_shared_on_pieces(
        &self,
        sedges: &[TopoShape],
    ) -> (
        HashMap<ShapeMapKey, Vec<ShapeListOfShape>>,
        HashMap<ShapeMapKey, Vec<TopoShape>>,
    ) {
        let mut mein: HashMap<ShapeMapKey, Vec<ShapeListOfShape>> = HashMap::new();
        let mut meout: HashMap<ShapeMapKey, Vec<TopoShape>> = HashMap::new();

        for e in sedges {
            let Some(hds) = self.data_structure() else {
                return (mein, meout);
            };
            let bds = hds.ds();
            let ie = hds.shape_index(e, true);
            let re = bds.ancestor_rank_of(e);

            // `cxx:348-352`
            if !self.is_split(e, State::On) {
                continue;
            }
            let lesd: Vec<TopoShape> = bds.shape_same_domain(e).to_vec();
            if lesd.is_empty() {
                continue;
            }

            let lespon: Vec<TopoShape> = self.splits(e, State::On).to_vec();
            let mut leoutlesd = lespon.clone();
            let iref = bds.same_domain_ref_of(e);

            for esd in &lesd {
                let iesd = hds.shape_index(esd, true);
                let lesdspon: Vec<TopoShape> = self.splits(esd, State::On).to_vec();

                // `cxx:390`: reduce `LEoutLESD` to the pieces OUT of `ESD`.
                let mut dummy_los: Vec<TopoShape> = Vec::new();
                fun_sel_lee(&mut leoutlesd, esd, State::In, &mut dummy_los);

                // `cxx:394-428`: the edge that owns the pieces.
                let n_lespon = lespon.len();
                let n_lesdspon = lesdspon.len();
                let (e1, e2, le1): (&TopoShape, &TopoShape, &Vec<TopoShape>);
                if n_lespon != 0 && n_lesdspon != 0 {
                    let take_e = (re == 1 && iesd != iref) || ie == iref;
                    if take_e {
                        e1 = e;
                        e2 = esd;
                        le1 = &lespon;
                    } else {
                        e1 = esd;
                        e2 = e;
                        le1 = &lesdspon;
                    }
                } else if n_lespon != 0 {
                    e1 = e;
                    e2 = esd;
                    le1 = &lespon;
                } else if n_lesdspon != 0 {
                    e1 = esd;
                    e2 = e;
                    le1 = &lesdspon;
                } else {
                    continue;
                }

                // `cxx:432-440`: the pieces of E1 that lie IN E2.
                let mut le1loc = le1.clone();
                let mut le1ine2: Vec<TopoShape> = Vec::new();
                fun_sel_lee(&mut le1loc, e2, State::In, &mut le1ine2);
                if le1ine2.is_empty() {
                    continue;
                }

                // `cxx:447-465`: the two cells of the pair are created even
                // when the pair is not recorded, because pass 3 tests for them.
                mein.entry(shape_map_key(e1)).or_default();
                mein.entry(shape_map_key(e2)).or_default();

                // `cxx:468-470`
                let e2_elem_of_e1 = fun_find_e_in_slos(
                    e2,
                    mein.get(&shape_map_key(e1)).map_or(&[][..], |v| v.as_slice()),
                );
                let e1_elem_of_e2 = fun_find_e_in_slos(
                    e1,
                    mein.get(&shape_map_key(e2)).map_or(&[][..], |v| v.as_slice()),
                );
                if e2_elem_of_e1 || e1_elem_of_e2 {
                    continue;
                }

                // `cxx:474-488`: regularise the shared pieces.
                let mut lr: Vec<TopoShape> = Vec::new();
                for piece in &le1ine2 {
                    let mut newle: Vec<TopoShape> = Vec::new();
                    if split_e(piece, &mut newle) {
                        lr.extend(newle);
                    } else {
                        lr.push(piece.clone());
                    }
                }
                // `cxx:490-503`
                mein
                    .entry(shape_map_key(e1))
                    .or_default()
                    .push(ShapeListOfShape::with_list(e2, lr.clone()));
                mein
                    .entry(shape_map_key(e2))
                    .or_default()
                    .push(ShapeListOfShape::with_list(e1, lr));
            }

            // `cxx:506-512`
            meout
                .entry(shape_map_key(e))
                .or_default()
                .extend(leoutlesd.iter().cloned());
        }

        (mein, meout)
    }

    /// `cxx:516-577`, pass 3: replace the `Split(ON)` list of every section
    /// edge that took part in pass 2 by its `MEOUT` pieces followed by every
    /// `MEIN` cell's pieces.
    fn rewrite_on_lists(
        &mut self,
        sedges: &[TopoShape],
        mein: &HashMap<ShapeMapKey, Vec<ShapeListOfShape>>,
        meout: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    ) {
        for e in sedges {
            let Some(hds) = self.data_structure() else {
                return;
            };
            if !self.is_split(e, State::On) {
                continue;
            }
            if hds.ds().shape_same_domain(e).is_empty() {
                continue;
            }
            let key = shape_map_key(e);
            let isb_meout = meout.contains_key(&key);
            let isb_mein = mein.contains_key(&key);
            if !isb_meout && !isb_mein {
                continue;
            }

            let mut acc: Vec<TopoShape> = Vec::new();
            if let Some(l) = meout.get(&key) {
                acc.extend(l.iter().cloned());
            }
            if let Some(loslos) = mein.get(&key) {
                for slos in loslos {
                    acc.extend(slos.list().iter().cloned());
                }
            }

            let lespon = self.change_split(e, State::On);
            lespon.clear();
            lespon.extend(acc);
        }
    }

    /// `cxx:579-695`, pass 4: for every reference section edge, give each of
    /// its `Split(ON)` pieces the p-curve of every face connected to one of its
    /// same-domain partners.
    fn attach_pcurves_on_connected_faces(&self, sedges: &[TopoShape]) {
        let Some(hds) = self.data_structure() else {
            return;
        };
        // See the module docs: the port rebuilds the connexity tables the DS
        // filler had prepared in OCCT (`cxx:619`).
        let mut connex = Connexity::new();
        connex.prepare(hds);

        for e in sedges {
            let Some(hds) = self.data_structure() else {
                return;
            };
            let bds = hds.ds();
            let lesd: Vec<TopoShape> = bds.shape_same_domain(e).to_vec();
            if lesd.is_empty() {
                continue;
            }
            // `cxx:600-608`
            let ie = bds.shape_index(e, true);
            let re = bds.same_domain_ref_of(e);
            if ie != re {
                continue;
            }

            // `cxx:610-622`: every face connected through a partner edge.
            let mut lf: Vec<TopoShape> = Vec::new();
            for esd in &lesd {
                lf.extend(connex.edge_connexity_same_shape(esd, hds).iter().cloned());
            }

            let lon: Vec<TopoShape> = self.splits(e, State::On).to_vec();
            for f in &lf {
                let if_ = bds.shape_index(f, true);
                let rf = bds.ancestor_rank(if_);

                // `cxx:636-662`: the partner edge of the same rank inside `F`,
                // i.e. the `EF` argument of the seven-argument
                // `FC2D_CurveOnSurface`.
                let mut esd_f: Option<TopoShape> = None;
                for esd in &lesd {
                    let iesd = bds.shape_index(esd, true);
                    if bds.ancestor_rank(iesd) != rf {
                        continue;
                    }
                    for ee in explorer_shapes(f, ShapeType::Edge) {
                        if is_equal(&ee, esd) {
                            esd_f = Some(esd.clone());
                            break;
                        }
                    }
                    if esd_f.is_some() {
                        break;
                    }
                }

                // `cxx:624`: `tolE = BRep_Tool::Tolerance(E)`, with `E` the
                // reference section edge.
                let tol_e = BRepTool::edge_tolerance(&Edge(e.clone()));

                for eon in &lon {
                    let edge = Edge(eon.clone());
                    let face = Face(f.clone());
                    // `cxx:668-672`
                    if crate::boptools_2d::has_curve_on_surface(&edge, &face) {
                        continue;
                    }
                    // `cxx:684-691`: a null p-curve raises
                    // `Standard_ProgramError("...SSE null PC on F")`.
                    let c2df = fc2d_curve_on_surface_with_ef(eon, f, esd_f.as_ref(), true)
                        .expect("TopOpeBRepBuild_Builder::SSE null PC on F");
                    // `cxx:692-693`: `tol = max(tolE, tolpc)` then
                    // `BB.UpdateEdge(eon, PC, F, tol)`.
                    let tol = if tol_e > c2df.tol { tol_e } else { c2df.tol };
                    if let Some(pc) = c2df.curve {
                        write_edge_pcurve(&Edge(eon.clone()), pc, &face, tol);
                    }
                }
            }
        }
    }

    /// `cxx:697-780`, pass 5: for every section edge meeting a periodic
    /// cylinder / cone / sphere through a double face interference, build the
    /// p-curve on that face. The result is dropped — the call only fills the
    /// p-curve table (`FC2D_MakeCurveOnSurface` registers what it builds).
    fn materialise_periodic_pcurves(&self) {
        let Some(hds) = self.data_structure() else {
            return;
        };
        let bds = hds.ds();
        let nsha = bds.nb_shapes();
        for i in 1..=nsha {
            let Some(ff) = hds.shape(i, true) else {
                continue;
            };
            if ff.shape_type() != ShapeType::Face {
                continue;
            }
            let li = bds.shape_interferences(&ff, true);
            if li.is_empty() {
                continue;
            }
            let mut jobs: Vec<(TopoShape, TopoShape)> = Vec::new();
            for handle in li {
                let it = handle.read().expect("poisoned Interference lock");
                // `SSI` is null for every non shape/shape interference.
                if it.g_bound().is_none() {
                    continue;
                }
                let (_, gi, st, si) = fds_data(&it);
                if st != Kind::Face {
                    continue;
                }
                let tfe = it.transition();
                if tfe.shape_after() != ShapeType::Face || tfe.shape_before() != ShapeType::Face
                {
                    continue;
                }
                let Some(fs) = hds.shape(si, true) else {
                    continue;
                };
                // `cxx:742-746`
                if !fun_periodic_s(&fs) {
                    continue;
                }
                let Some(eg) = hds.shape(gi, true) else {
                    continue;
                };
                // `cxx:752-756`
                if !bds.is_section_edge(&eg, true) {
                    continue;
                }
                jobs.push((eg, fs));
            }

            for (eg, fs) in jobs {
                // `cxx:760-764`
                if crate::boptools_2d::has_curve_on_surface(
                    &Edge(eg.clone()),
                    &Face(fs.clone()),
                ) {
                    continue;
                }
                // `cxx:766-769`
                if BRepTool::edge_curve(&Edge(eg.clone())).is_none() {
                    panic!("TopOpeBRepBuild_Builder::SSE EG without C3D");
                }
                // `cxx:770-778`
                if fc2d_curve_on_surface(&eg, &fs).is_none() {
                    panic!("TopOpeBRepBuild_Builder::SSE EG without PC on FS");
                }
            }
        }
    }

    /// `cxx:782-913`, pass 6: on a periodic face, translate the iso-parametric
    /// `Split(ON)` pieces of a section edge back onto the face bounds.
    fn translate_on_bounds_pieces(&mut self) {
        let Some(hds) = self.data_structure() else {
            return;
        };
        let bds = hds.ds();
        let nsha = bds.nb_shapes();
        // Collected first: the rewrite below needs `&mut self`.
        let mut pairs: Vec<(TopoShape, TopoShape, bool, bool)> = Vec::new();
        for i in 1..=nsha {
            let Some(ff) = hds.shape(i, true) else {
                continue;
            };
            if ff.shape_type() != ShapeType::Face {
                continue;
            }
            // `cxx:794-799`
            let Some((uper, vper)) = fun_periodic(&Face(ff.clone())) else {
                continue;
            };
            if !(uper || vper) {
                continue;
            }
            let li = bds.shape_interferences(&ff, true);
            if li.is_empty() {
                continue;
            }
            for handle in li {
                let it = handle.read().expect("poisoned Interference lock");
                if it.g_bound().is_none() {
                    continue;
                }
                let (_, gi, st, _) = fds_data(&it);
                if st != Kind::Face {
                    continue;
                }
                // `cxx:822-826`: `GBound() == 1` is skipped.
                if it.g_bound() == Some(true) {
                    continue;
                }
                let Some(eg) = hds.shape(gi, true) else {
                    continue;
                };
                if !bds.is_section_edge(&eg, true) {
                    continue;
                }
                pairs.push((eg, ff.clone(), uper, vper));
            }
        }

        for (eg, ff, uper, vper) in pairs {
            // `cxx:840-848`
            let Some(pc) = fc2d_curve_on_surface(&eg, &ff) else {
                // OCCT hands the null `PC` to `UVISO`, which dereferences it.
                continue;
            };
            let Some(iso) = uv_iso_curve(pc.as_ref()) else {
                continue;
            };

            // `cxx:850-910`: translate every piece, then replace the list.
            let new_on: Vec<TopoShape> = self.splits(&eg, State::On).to_vec();
            let mut acc: Vec<TopoShape> = Vec::new();
            for eon in new_on {
                // `cxx:857-861`
                let Some(pcon) = fc2d_curve_on_surface(&eon, &ff) else {
                    panic!("TopOpeBRepBuild_Builder::splitON");
                };
                let Some(iso_on) = uv_iso_curve(pcon.as_ref()) else {
                    panic!("TopOpeBRepBuild_Builder::splitON");
                };
                // `cxx:863`
                let pc_ko = (iso.iso_u && !iso_on.iso_u) || (iso.iso_v && !iso_on.iso_v);
                if pc_ko {
                    panic!("TopOpeBRepBuild_Builder::splitON");
                }

                // `cxx:870-874`
                let test = (uper && iso.iso_u) || (vper && iso.iso_v);
                if !test {
                    acc.push(eon);
                    continue;
                }

                // `cxx:876-881`
                let (is_on_clo, uv_ok) = fun_onboundsper_edge(&eon, &ff);
                if is_on_clo {
                    acc.push(eon);
                    continue;
                }

                // `cxx:884-905`
                if fun_onboundsper_uv(&iso_on.origin, &ff) {
                    if let Some(uv_ok) = uv_ok {
                        let dtr = if iso.iso_u {
                            uv_ok.x() - iso_on.origin.x()
                        } else {
                            uv_ok.y() - iso_on.origin.y()
                        };
                        let tr = if iso.iso_u {
                            GpVec2d::new(dtr, 0.0)
                        } else {
                            GpVec2d::new(0.0, dtr)
                        };
                        let mag = tr.magnitude();
                        let toluv = 1.0e-8 * 1.0e2; // `// NYINYI`
                        if mag > toluv {
                            trsl_uv_modif_e(&tr, &ff, &eon);
                        }
                    }
                }
                acc.push(eon);
            }

            let l = self.change_split(&eg, State::On);
            l.clear();
            l.extend(acc);
        }
    }
}
