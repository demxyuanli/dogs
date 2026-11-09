//! Face-split side of `BRepOffset_BuildOffsetFaces`: the two split builders
//! and the methods that drive them
//! (`BuildSplitsOfTrimmedFace` `_1.cxx:411-433`,
//! `BuildSplitsOfFace` `_1.cxx:438-486`,
//! `BuildSplitsOfTrimmedFaces` `_1.cxx:1036-1086`,
//! `FillHistory` `_1.cxx:9413-9496`).

use crate::abs::{Orientation, ShapeType};
use crate::bop_builder2::BopBuilder;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_offset::inter3d::ShapeIndexedMap;
use crate::builder_area::AreaBuilder;
use crate::builder_face::FaceBuilder;
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::{explorer_shapes, shape_map_key};

use super::statics::has_data;
use super::BRepOffsetBuildOffsetFaces;

/// `static void BuildSplitsOfTrimmedFace(const TopoDS_Face&,
/// const TopoDS_Shape&, NCollection_List<TopoDS_Shape>&,
/// const Message_ProgressRange&)` (`_1.cxx:411-433`).
///
/// `BOPAlgo_Splitter aSplitter; aSplitter.AddArgument(theFace);
/// aSplitter.AddArgument(theEdges); aSplitter.SetToFillHistory(false);
/// aSplitter.Perform()` (`_1.cxx:415-421`) is carried by
/// [`BopBuilder`]: `BOPAlgo_Splitter::Perform`
/// (`BOPAlgo_Splitter.cxx:57-95`) builds the pave filler over
/// `arguments + tools` and then runs `PerformInternal` - the general-fuse
/// build, which is exactly [`BopBuilder::perform`] after
/// [`BopBuilder::set_arguments`]. Its `BuildResult`
/// (`BOPAlgo_Splitter.cxx:97-119`) only unwraps a single-argument result,
/// which cannot fire here because the call always passes two arguments.
/// `SetToFillHistory(false)` only skips the history fill, and this port never
/// reads the splitter's history (`[`BopBuilder::fill_history`] stays `false`
/// on a fresh builder).
///
/// The result loop (`_1.cxx:428-432`) appends every `TopAbs_FACE` of the
/// splitter's shape, i.e. the pieces the face was cut into; the edge arguments
/// cannot contribute faces.
pub(super) fn build_splits_of_trimmed_face(
    the_face: &Face,
    the_edges: &TopoShape,
    the_lf_images: &mut Vec<TopoShape>,
) {
    // `_1.cxx:415-421`.
    let mut a_splitter = BopBuilder::new();
    a_splitter.add_argument(&the_face.0);
    a_splitter.add_argument(the_edges);
    a_splitter.set_to_fill_history(false);
    if a_splitter.perform().is_err() || a_splitter.has_errors() {
        // `_1.cxx:422-425`: a failed splitter contributes nothing.
        return;
    }

    // `_1.cxx:428-432`.
    let a_shape = a_splitter.result();
    for a_f in explorer_shapes(a_shape, ShapeType::Face) {
        the_lf_images.push(a_f);
    }
}

/// `static void BuildSplitsOfFace(const TopoDS_Face&, const TopoDS_Shape&,
/// DataMap<Shape, Shape>& theFacesOrigins, List& theLFImages)`
/// (`_1.cxx:438-486`).
///
/// Two OCCT dependencies map onto existing ports:
///
/// * `BRepLib::BuildPCurveForEdgesOnPlane(aLE, aFF)` (`_1.cxx:463`) - the
///   inline loop of `BRepLib.hxx:122-131` over
///   `BRepLib::BuildPCurveForEdgeOnPlane` (`BRepLib_1.cxx:313-341`), which
///   reads `BRep_Tool::CurveOnSurface(aE, aF, ..., &isStored)` and stores the
///   computed p-curve when one had to be computed. That is
///   [`crate::boptools_2d::build_pcurve_for_edge_on_face`].
/// * `BOPAlgo_BuilderFace aBF; aBF.SetFace(aFF); aBF.SetShapes(aLE);
///   aBF.Perform()` (`_1.cxx:466-469`) - [`FaceBuilder`] with
///   [`AreaBuilder::set_face`] / [`AreaBuilder::set_shapes`] /
///   [`AreaBuilder::perform`]; `aBF.Areas()` (`_1.cxx:476`) is
///   [`AreaBuilder::areas`].
pub(super) fn build_splits_of_face(
    the_face: &Face,
    the_edges: &TopoShape,
    the_faces_origins: &mut std::collections::HashMap<crate::topo_tools::ShapeMapKey, TopoShape>,
    the_lf_images: &mut Vec<TopoShape>,
) {
    // `_1.cxx:444`.
    the_lf_images.clear();

    // `_1.cxx:446-455`: every edge enters twice, FORWARD then REVERSED.
    let mut a_le: Vec<TopoShape> = Vec::new();
    for a_e in explorer_shapes(the_edges, ShapeType::Edge) {
        let mut a_forward = a_e.clone();
        a_forward.set_orientation(Orientation::Forward);
        a_le.push(a_forward);
        let mut a_reversed = a_e;
        a_reversed.set_orientation(Orientation::Reversed);
        a_le.push(a_reversed);
    }

    // `_1.cxx:458-460`: the face is worked on in FORWARD orientation, the
    // original orientation is restored on every piece afterwards.
    let mut a_ff = the_face.clone();
    let an_or = the_face.0.orientation();
    a_ff.0.set_orientation(Orientation::Forward);

    // `_1.cxx:463`: build the p-curves on the plane.
    for a_e in &a_le {
        if let Some(a_edge) = Edge::wrap(a_e.clone()) {
            // A face without a registered surface cannot carry a p-curve; the
            // reference's `CurveOnSurface` would raise there, so the port skips
            // the edge instead of aborting the whole split.
            let _ = crate::boptools_2d::build_pcurve_for_edge_on_face(&a_edge, &a_ff);
        }
    }

    // `_1.cxx:466-474`.
    let mut a_bf = FaceBuilder::new();
    a_bf.set_face(&a_ff);
    a_bf.set_shapes(&a_le);
    if a_bf.perform().is_err() {
        return;
    }

    // `_1.cxx:476-485`.
    for a_f_sp in a_bf.areas() {
        let mut a_f_sp = a_f_sp.clone();
        a_f_sp.set_orientation(an_or);
        the_lf_images.push(a_f_sp.clone());
        the_faces_origins.insert(shape_map_key(&a_f_sp), the_face.0.clone());
    }
}

impl BRepOffsetBuildOffsetFaces {
    /// `void BRepOffset_BuildOffsetFaces::BuildSplitsOfTrimmedFaces(const
    /// Message_ProgressRange&)` (`_1.cxx:1036-1086`).
    ///
    /// `myImage` is passed in as `the_image` (see the module note of
    /// [`super::BRepOffsetBuildOffsetFaces`]). The reference's local
    /// `anEdgesOrigins` (`_1.cxx:1044-1045`, bound as a fallback when
    /// `myEdgesOrigins` is null, `_1.cxx:1046-1049`) has no counterpart: the
    /// port owns `my_edges_origins`, so it is never null.
    pub fn build_splits_of_trimmed_faces(&mut self, the_image: &mut BRepAlgoImage) {
        // `_1.cxx:1038-1041`.
        if !has_data(&self.my_faces) {
            return;
        }

        // `_1.cxx:1054`.
        self.intersect_trimmed_edges();

        // `_1.cxx:1056-1081`.
        for a_f_shape in self.my_faces.clone() {
            let Some(a_f) = Face::wrap(a_f_shape) else {
                continue;
            };
            // `_1.cxx:1064-1066`.
            let (a_ce, b_found) = self.get_edges(&a_f, None);

            // `_1.cxx:1069-1076`: nothing to split - keep the face as its own
            // image when the history does not know it yet.
            if !b_found {
                if !the_image.has_image(&a_f.0) {
                    // `myOFImages(myOFImages.Add(aF, List())).Append(aF)`:
                    // `Add` replaces the value under an existing key.
                    self.my_of_images.bind_list(&a_f.0, &[a_f.0.clone()]);
                }
                continue;
            }

            // `_1.cxx:1078-1080`.
            let mut a_lf_images: Vec<TopoShape> = Vec::new();
            build_splits_of_trimmed_face(&a_f, &a_ce, &mut a_lf_images);
            self.my_of_images.bind_list(&a_f.0, &a_lf_images);
        }

        // `_1.cxx:1083`.
        self.fill_history(the_image);
    }

    /// `void BRepOffset_BuildOffsetFaces::FillHistory()` (`_1.cxx:9413-9496`).
    ///
    /// `myImage->HasImage / Add / Bind` (`_1.cxx:9446-9454`, `:9477-9488`) map
    /// onto [`BRepAlgoImage::has_image`] / [`BRepAlgoImage::add`] /
    /// [`BRepAlgoImage::bind_list`]. `Add` on a root that is already bound
    /// throws in OCCT; the port ignores the `Err`, like every other
    /// `BRepOffset` module.
    ///
    /// Only the edges that appear in some face split are recorded
    /// (`_1.cxx:9476-9488`), which is what `anEdgesMap` collects at
    /// `_1.cxx:9450-9460`.
    pub fn fill_history(&self, the_image: &mut BRepAlgoImage) {
        // `_1.cxx:9415-9418`.
        let a_nb_f = self.my_of_images.extent();
        if a_nb_f == 0 {
            return;
        }

        // `_1.cxx:9428`: the edges kept in the splits.
        let mut an_edges_map = ShapeIndexedMap::new();

        // `_1.cxx:9431-9461`: history of the faces.
        for i in 1..=a_nb_f {
            let Some(a_f) = self.my_of_images.keys().get(i - 1).cloned() else {
                continue;
            };
            let a_lf_images: Vec<TopoShape> = match self.my_of_images.seek(&a_f) {
                Some(a_list) => a_list.to_vec(),
                None => continue,
            };
            // `_1.cxx:9434-9437`.
            if a_lf_images.is_empty() {
                continue;
            }

            // `_1.cxx:9440-9448`.
            if the_image.has_image(&a_f) {
                for a_item in &a_lf_images {
                    let _ = the_image.add(&a_f, a_item);
                }
            } else {
                let _ = the_image.bind_list(&a_f, &a_lf_images);
            }

            // `_1.cxx:9450-9460`.
            for a_f_im in &a_lf_images {
                for a_e in explorer_shapes(a_f_im, ShapeType::Edge) {
                    an_edges_map.add(&a_e);
                }
            }
        }

        // `_1.cxx:9464-9490`: history of the edges whose images survived.
        let a_oe: Vec<TopoShape> = self.my_oe_images.keys().to_vec();
        for a_e in &a_oe {
            let a_le_im: Vec<TopoShape> = match self.my_oe_images.seek(a_e) {
                Some(a_list) => a_list.to_vec(),
                None => continue,
            };
            let mut b_has_image = the_image.has_image(a_e);
            for a_e_im in &a_le_im {
                if an_edges_map.contains(a_e_im) {
                    if b_has_image {
                        let _ = the_image.add(a_e, a_e_im);
                    } else {
                        let _ = the_image.bind(a_e, a_e_im);
                        b_has_image = true;
                    }
                }
            }
        }
    }
}
