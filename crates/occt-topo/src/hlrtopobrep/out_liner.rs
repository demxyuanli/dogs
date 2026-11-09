//! `HLRTopoBRep_OutLiner` (`ModelingAlgorithms/TKHLR/HLRTopoBRep/HLRTopoBRep_OutLiner.hxx:34-75`,
//! `.cxx:44-346`, `.lxx:19-55`): the shape whose outlines are extracted for the
//! exact HLR pass.
//!
//! `Fill` (`cxx:65-93`) runs the two algorithms this package delegates to:
//! * `Contap_Contour` (`cxx:74-86`) - `Init(Eye)` / `Init(DirZ)` and the
//!   contour extraction itself ([`ContapContour`], `contap/contour/`);
//! * `HLRTopoBRep_DSFiller::Insert(S, FO, myDS, MST, nbIso)` (`cxx:87`,
//!   [`super::ds_filler::insert`]).
//!
//! `ProcessFace` (`cxx:99-307`) and `BuildShape` (`cxx:309-346`) - the two
//! consumers of `HLRTopoBRep_Data` - are ported in full.
//!
//! Deviations:
//! * `TopoDS_Shape myOriginalShape` / `myOutLinedShape` (`hxx:72-73`) are
//!   nullable in OCCT (`IsNull()` at `cxx:69-70`); the port stores
//!   `Option<TopoShape>`.
//! * `BRep_Builder::MakeWire` / `Add` are [`TopoBuilder::make_shape`] /
//!   [`TopoBuilder::add`]; `MakeWire` + `Add` is the incremental builder form
//!   OCCT uses here (`cxx:128-129`, `:151`).
//! * `TopExp_Explorer(shape, ToFind, ToAvoid)` is
//!   [`explorer_shapes_avoiding`] (`cxx:331`, `:340`).
//! * `MST.ChangeFind(F)` (`cxx:302`) throws `Standard_OutOfRange` for an unbound
//!   face; the port only re-binds when the entry exists (the map is filled by
//!   the unported `DSFiller`).
//! * `OutL` (`cxx:155`) is a reference into `myDS`; the port re-enters
//!   [`HlrTopoBRepData::add_out_l`] for every append, which yields the same list
//!   content and order.

use std::collections::{HashMap, HashSet};

use occt_core::gp::{GpPnt, GpVec};
use occt_geom::extrema_pc::extrema_ext_pc_range_done;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::brep_top_adaptor_tool::BrepTopAdaptorTool;
use crate::builder::TopoBuilder;
use crate::contap::ContapContour;
use crate::geom_abs_curve_type::{curve_kind, GeomAbsCurveType};
use crate::hlr_algo::HlrAlgoProjector;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{
    edge_vertices, explorer_shapes, explorer_shapes_avoiding, is_same_occurrence,
    map_shapes_and_ancestors, shape_map_key, ShapeMapKey,
};

use super::data::HlrTopoBRepData;

/// `gp_Pnt::RealLast()` (`Standard_Real.hxx`), the initial `dist` of
/// `HLRBTopoBRep_OutLiner.cxx:203`.
const REAL_LAST: f64 = f64::MAX;

/// The `MST` argument type of `Fill` / `ProcessFace` / `BuildShape`:
/// `NCollection_DataMap<TopoDS_Shape, BRepTopAdaptor_Tool, TopTools_ShapeMapHasher>`.
pub type ShapeToTopAdaptorToolMap = HashMap<ShapeMapKey, BrepTopAdaptorTool>;

/// `HLRTopoBRep_OutLiner`.
#[derive(Debug, Clone, Default)]
pub struct HlrTopoBRepOutLiner {
    /// `TopoDS_Shape myOriginalShape` (`hxx:72`).
    my_original_shape: Option<TopoShape>,
    /// `TopoDS_Shape myOutLinedShape` (`hxx:73`).
    my_out_lined_shape: Option<TopoShape>,
    /// `HLRTopoBRep_Data myDS` (`hxx:74`).
    my_ds: HlrTopoBRepData,
}

impl HlrTopoBRepOutLiner {
    /// `HLRTopoBRep_OutLiner()` (`cxx:44`): `= default`, both shapes null.
    pub fn new() -> Self {
        Self::default()
    }

    /// `HLRTopoBRep_OutLiner(OriS)` (`cxx:46-51`).
    pub fn with_original(ori_s: TopoShape) -> Self {
        Self {
            my_original_shape: Some(ori_s),
            my_out_lined_shape: None,
            my_ds: HlrTopoBRepData::new(),
        }
    }

    /// `HLRTopoBRep_OutLiner(OriS, OutS)` (`cxx:53-59`).
    pub fn with_shapes(ori_s: TopoShape, out_s: TopoShape) -> Self {
        Self {
            my_original_shape: Some(ori_s),
            my_out_lined_shape: Some(out_s),
            my_ds: HlrTopoBRepData::new(),
        }
    }

    /// `OriginalShape(OriS)` (`lxx:19-22`).
    pub fn set_original_shape(&mut self, ori_s: Option<TopoShape>) {
        self.my_original_shape = ori_s;
    }

    /// `OriginalShape()` (`lxx:24-27`).
    pub fn original_shape(&self) -> Option<&TopoShape> {
        self.my_original_shape.as_ref()
    }

    /// The mutable `TopoDS_Shape&` of `lxx:24-27`.
    pub fn original_shape_mut(&mut self) -> Option<&mut TopoShape> {
        self.my_original_shape.as_mut()
    }

    /// `OutLinedShape(OutS)` (`lxx:29-32`).
    pub fn set_out_lined_shape(&mut self, out_s: Option<TopoShape>) {
        self.my_out_lined_shape = out_s;
    }

    /// `OutLinedShape()` (`lxx:34-37`).
    pub fn out_lined_shape(&self) -> Option<&TopoShape> {
        self.my_out_lined_shape.as_ref()
    }

    /// The mutable `TopoDS_Shape&` of `lxx:34-37`.
    pub fn out_lined_shape_mut(&mut self) -> Option<&mut TopoShape> {
        self.my_out_lined_shape.as_mut()
    }

    /// `DataStructure()` (`lxx:39-42`).
    pub fn data_structure(&self) -> &HlrTopoBRepData {
        &self.my_ds
    }

    /// The non-const `DataStructure()` of `lxx:39-42`.
    pub fn data_structure_mut(&mut self) -> &mut HlrTopoBRepData {
        &mut self.my_ds
    }

    /// `Fill(P, MST, nbIso)` (`cxx:65-93`).
    ///
    /// `cxx:69-91` in order:
    /// 1. `myOriginalShape.IsNull()` guard;
    /// 2. `myOutLinedShape.IsNull()` guard;
    /// 3. `Vecz = (0,0,1)` transformed by `P.Transformation().Inverted()`;
    /// 4. `Contap_Contour FO; FO.Init(Eye)` with `Eye = P.Focus() * Vecz.XYZ()`
    ///    under a perspective projector, `FO.Init(DirZ)` with `DirZ(Vecz)`
    ///    otherwise;
    /// 5. `HLRTopoBRep_DSFiller::Insert(myOriginalShape, FO, myDS, MST, nbIso)`;
    /// 6. `BuildShape(MST)`.
    pub fn fill(
        &mut self,
        p: &HlrAlgoProjector,
        mst: &mut ShapeToTopAdaptorToolMap,
        nb_iso: i32,
    ) -> Result<(), String> {
        let Some(original) = self.my_original_shape.clone() else {
            return Ok(()); // `cxx:69`
        };
        if self.my_out_lined_shape.is_some() {
            return Ok(()); // `cxx:70`
        }

        // `cxx:72-74`: `gp_Vec Vecz(0,0,1); gp_Trsf Tr(P.Transformation());
        // Tr.Invert(); Vecz.Transform(Tr);`
        let mut vecz = GpVec::new(0.0, 0.0, 1.0);
        let tr = p.transformation().inverted().map_err(|_| {
            "HLRTopoBRep_OutLiner::Fill: the projector transformation is singular".to_string()
        })?;
        vecz.transform(&tr);

        // `cxx:75-86`: `Contap_Contour FO;` and the mode `Init` it receives.
        let mut fo = ContapContour::new();
        if p.perspective() {
            let focus = p.focus().map_err(|_| {
                "HLRTopoBRep_OutLiner::Fill: projector is not perspective".to_string()
            })?;
            // `gp_Pnt Eye; Eye.SetXYZ(P.Focus() * Vecz.XYZ()); FO.Init(Eye);`
            fo.init_eye(&GpPnt::from_xyz(&vecz.coord.multiplied(focus)));
        } else {
            // `gp_Dir DirZ(Vecz); FO.Init(DirZ);` - `Contap_Contour::Init(const
            // gp_Dir&)` is the port's `init_direction` (`contour/mod.rs:154-161`),
            // which normalises through `GpDir::from_vec`, the same
            // `gp_Dir(const gp_Vec&)` construction.
            fo.init_direction(&vecz).map_err(|_| {
                "HLRTopoBRep_OutLiner::Fill: degenerate view direction".to_string()
            })?;
        }

        // `HLRTopoBRep_DSFiller::Insert(myOriginalShape, FO, myDS, MST, nbIso);`
        // (`cxx:87`).
        super::ds_filler::insert(&original, &mut fo, &mut self.my_ds, mst, nb_iso)?;

        // `BuildShape(MST);` (`cxx:88`).
        self.build_shape(mst);
        Ok(())
    }

    /// `ProcessFace(F, S, MST)` (`cxx:99-307`): build the outlined copy `NF` of
    /// face `F` from `myDS` and add it to `S`.
    pub fn process_face(
        &mut self,
        f: &Face,
        s: &mut TopoShape,
        mst: &mut ShapeToTopAdaptorToolMap,
    ) {
        let b = TopoBuilder;
        // `TopExp::MapShapesAndAncestors(F, TopAbs_VERTEX, TopAbs_EDGE, aVEMap)`
        // (`cxx:109-110`).
        let a_ve_map = map_shapes_and_ancestors(f.shape(), ShapeType::Vertex, ShapeType::Edge);

        // `NF = F.EmptyCopied()` (`cxx:125`).
        let mut nf = f.shape().empty_copied();

        // `for (exW.Init(F, TopAbs_WIRE); ...)` (`cxx:127-152`).
        for w_shape in explorer_shapes(f.shape(), ShapeType::Wire) {
            let mut w = TopoShape::new(ShapeType::Wire); // `B.MakeWire(W)` (`cxx:129`)
            for e_shape in explorer_shapes(&w_shape, ShapeType::Edge) {
                let Some(e) = Edge::wrap(e_shape) else {
                    continue;
                };
                if self.my_ds.edge_has_spl_e(&e) {
                    let splits: Vec<TopoShape> = self.my_ds.edge_spl_e(&e).to_vec();
                    for piece in splits {
                        let Some(mut new_e) = Edge::wrap(piece) else {
                            continue;
                        };
                        new_e.0.set_orientation(e.shape().orientation());
                        self.my_ds.add_old_s(new_e.shape(), e.shape());
                        b.add(&mut w, new_e.shape());
                    }
                } else {
                    b.add(&mut w, e.shape());
                }
            }
            b.add(&mut nf, &w); // add the new wire in the new face (`cxx:151`)
        }

        // `myDS.AddIntL(F); NCollection_List& OutL = myDS.AddOutL(F);`
        // (`cxx:154-155`). `AddOutL` creates the entry the `OutL.Append` below
        // targets.
        self.my_ds.add_int_l(f);
        self.my_ds.add_out_l(f);

        if self.my_ds.face_has_int_l(f) {
            self.append_outlines_on_restriction(f, &mut nf, &a_ve_map);
        }

        if self.my_ds.face_has_iso_l(f) {
            self.append_iso_lines(f, &mut nf);
        }

        // `myDS.AddOldS(NF, F); MST.Bind(NF, MST.ChangeFind(F));`
        // (`cxx:300-302`).
        self.my_ds.add_old_s(&nf, f.shape());
        if let Some(tool) = mst.get(&shape_map_key(f.shape())).cloned() {
            mst.insert(shape_map_key(&nf), tool);
        }

        b.add(s, &nf); // add the face in the shell (`cxx:304`)
    }

    /// The `myDS.FaceHasIntL(F)` arm of `ProcessFace` (`cxx:157-258`): the
    /// internal outlines, with the "outline coincides with a real edge" test at
    /// `cxx:168-222` and the internal-outline wire of `cxx:224-257`.
    fn append_outlines_on_restriction(
        &mut self,
        f: &Face,
        nf: &mut TopoShape,
        a_ve_map: &HashMap<ShapeMapKey, Vec<TopoShape>>,
    ) {
        let b = TopoBuilder;
        let mut w: Option<TopoShape> = None;
        let int_l: Vec<TopoShape> = self.my_ds.face_int_l(f).to_vec();

        for e_shape in int_l {
            let Some(mut e) = Edge::wrap(e_shape) else {
                continue;
            };
            e.0.set_orientation(Orientation::Internal);

            // `BRepAdaptor_Curve C(E); par = 0.34*First + 0.66*Last;
            // gp_Pnt P = C.Value(par);` (`cxx:168-170`).
            let (par, p) = match (BRepTool::edge_curve(&e), BRepTool::edge_parameters(&e)) {
                (Some(c), (u1, u2)) => {
                    let par = 0.34 * u1 + 0.66 * u2;
                    let p = c.value(par);
                    (Some(par), p)
                }
                _ => (None, GpPnt::new(0.0, 0.0, 0.0)),
            };
            let _ = par;

            let (v1, v2) = edge_vertices(&e); // `TopExp::Vertices(E, V1, V2)` (`cxx:172`)
            let mut same_edge = false;

            if let (Some(v1), Some(v2)) = (v1, v2) {
                // `if (!V1.IsNull() && aVEMap.Contains(V1))` (`cxx:175`)
                if let Some(candidates) = a_ve_map.get(&shape_map_key(v1.shape())) {
                    for a_e_shape in candidates {
                        let Some(a_e) = Edge::wrap(a_e_shape.clone()) else {
                            continue;
                        };
                        let (a_v1, a_v2) = edge_vertices(&a_e);
                        let (Some(a_v1), Some(a_v2)) = (a_v1, a_v2) else {
                            continue;
                        };
                        // `cxx:183`: the two endpoints must match in either
                        // order.
                        let same_ends = (is_same_occurrence(v1.shape(), a_v1.shape())
                            && is_same_occurrence(v2.shape(), a_v2.shape()))
                            || (is_same_occurrence(v1.shape(), a_v2.shape())
                                && is_same_occurrence(v2.shape(), a_v1.shape()));
                        if !same_ends {
                            continue;
                        }

                        let c_kind = BRepTool::edge_curve(&e)
                            .map(|c| curve_kind(c.as_ref()))
                            .unwrap_or(GeomAbsCurveType::OtherCurve);
                        let a_kind = BRepTool::edge_curve(&a_e)
                            .map(|c| curve_kind(c.as_ref()))
                            .unwrap_or(GeomAbsCurveType::OtherCurve);
                        if c_kind == GeomAbsCurveType::Line && a_kind == GeomAbsCurveType::Line {
                            // `cxx:186-190`: two lines sharing both endpoints
                            // are the same edge.
                            same_edge = true;
                            break;
                        }

                        // `Extrema_ExtPC anExt(P, aC)` (`cxx:194`), i.e.
                        // `Initialize(aC, First, Last)` then `Perform(P)`.
                        let (Some(a_c), (a_u1, a_u2)) =
                            (BRepTool::edge_curve(&a_e), BRepTool::edge_parameters(&a_e))
                        else {
                            continue;
                        };
                        let (sols, done) =
                            extrema_ext_pc_range_done(a_c.as_ref(), &p, a_u1, a_u2);
                        if done && !sols.is_empty() {
                            // `dist = std::min(dist, anExt.SquareDistance(ec))`
                            // over `ec = 1..NbExt` (`cxx:201-208`).
                            let mut dist = REAL_LAST;
                            for sol in &sols {
                                dist = dist.min(sol.sq_dist);
                            }
                            if dist <= 1.0e-14 {
                                same_edge = true; // `cxx:210-214`
                                break;
                            }
                        }
                    }
                }
            }

            if same_edge {
                // `OutL.Append(E)` (`cxx:225`).
                self.my_ds.add_out_l(f).push(e.shape().clone());
                continue;
            }

            if self.my_ds.edge_has_spl_e(&e) {
                let splits: Vec<TopoShape> = self.my_ds.edge_spl_e(&e).to_vec();
                for piece in splits {
                    let mut new_e = piece;
                    new_e.set_orientation(Orientation::Internal);
                    if w.is_none() {
                        w = Some(TopoShape::new(ShapeType::Wire));
                    }
                    self.my_ds.add_old_s(&new_e, f.shape());
                    b.add(w.as_mut().unwrap(), &new_e);
                }
            } else {
                if w.is_none() {
                    w = Some(TopoShape::new(ShapeType::Wire));
                }
                self.my_ds.add_old_s(e.shape(), f.shape());
                b.add(w.as_mut().unwrap(), e.shape());
            }
        }

        if let Some(w) = w {
            b.add(nf, &w); // add the new wire in the new face (`cxx:255`)
        }
    }

    /// The `myDS.FaceHasIsoL(F)` arm of `ProcessFace` (`cxx:260-299`): the iso
    /// lines of the face. OCCT notes they are never split (`cxx:270`) but
    /// still handles the split list.
    fn append_iso_lines(&mut self, f: &Face, nf: &mut TopoShape) {
        let b = TopoBuilder;
        let mut w: Option<TopoShape> = None;
        let iso_l: Vec<TopoShape> = self.my_ds.face_iso_l(f).to_vec();

        for e_shape in iso_l {
            let Some(mut e) = Edge::wrap(e_shape) else {
                continue;
            };
            e.0.set_orientation(Orientation::Internal);

            if self.my_ds.edge_has_spl_e(&e) {
                let splits: Vec<TopoShape> = self.my_ds.edge_spl_e(&e).to_vec();
                for piece in splits {
                    let mut new_e = piece;
                    new_e.set_orientation(Orientation::Internal);
                    if w.is_none() {
                        w = Some(TopoShape::new(ShapeType::Wire));
                    }
                    self.my_ds.add_old_s(&new_e, f.shape());
                    b.add(w.as_mut().unwrap(), &new_e);
                }
            } else {
                if w.is_none() {
                    w = Some(TopoShape::new(ShapeType::Wire));
                }
                self.my_ds.add_old_s(e.shape(), f.shape());
                b.add(w.as_mut().unwrap(), e.shape());
            }
        }

        if let Some(w) = w {
            b.add(nf, &w);
        }
    }

    /// `BuildShape(MST)` (`cxx:309-346`): make `myOutLinedShape` a compound of
    /// the re-built faces, shell by shell.
    pub fn build_shape(&mut self, mst: &mut ShapeToTopAdaptorToolMap) {
        let b = TopoBuilder;
        let Some(original) = self.my_original_shape.clone() else {
            return;
        };

        // `B.MakeCompound(TopoDS::Compound(myOutLinedShape));`
        let mut compound = TopoShape::new(ShapeType::Compound);
        // `NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> ShapeMap;`
        let mut shape_map: HashSet<ShapeMapKey> = HashSet::new();

        // `for (exshell.Init(myOriginalShape, TopAbs_SHELL); ...)` (`cxx:315-330`).
        for shell_shape in explorer_shapes(&original, ShapeType::Shell) {
            let mut the_shell = TopoShape::new(ShapeType::Shell);
            the_shell.set_closed(shell_shape.is_closed()); // `cxx:320`

            for face_shape in explorer_shapes(&shell_shape, ShapeType::Face) {
                if !shape_map.insert(shape_map_key(&face_shape)) {
                    continue; // `if (ShapeMap.Add(...))`
                }
                let Some(face) = Face::wrap(face_shape) else {
                    continue;
                };
                self.process_face(&face, &mut the_shell, mst);
            }
            b.add(&mut compound, &the_shell); // `B.Add(myOutLinedShape, theShell)`
        }

        // `for (exface.Init(myOriginalShape, TopAbs_FACE, TopAbs_SHELL); ...)`
        // (`cxx:332-338`): faces not in a shell.
        for face_shape in explorer_shapes_avoiding(&original, ShapeType::Face, ShapeType::Shell) {
            if !shape_map.insert(shape_map_key(&face_shape)) {
                continue;
            }
            let Some(face) = Face::wrap(face_shape) else {
                continue;
            };
            self.process_face(&face, &mut compound, mst);
        }

        // `for (exedge.Init(myOriginalShape, TopAbs_EDGE, TopAbs_FACE); ...)`
        // (`cxx:340-344`): edges not in a face.
        for edge_shape in explorer_shapes_avoiding(&original, ShapeType::Edge, ShapeType::Face) {
            b.add(&mut compound, &edge_shape);
        }

        self.my_out_lined_shape = Some(compound);
    }
}

