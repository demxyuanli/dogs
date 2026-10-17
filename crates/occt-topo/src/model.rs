//! BRep model — manage collections of shapes with names.
//! Provides a lightweight scene/document for CAD exchange.
use std::collections::HashMap;
use crate::shape::TopoShape;
use crate::abs::ShapeType;
use crate::xcaf_doc_vis_material::XcafDocVisMaterial;
use crate::xcaf_prs::XcafPrsStyle;
use occt_core::bnd::BndBox;

/// A named shape in the model.
#[derive(Debug, Clone)]
pub struct ModelShape {
    pub name: String,
    pub shape: TopoShape,
    /// The surface colour the `XCAFDoc_ColorTool` holds for this shape's own
    /// label (`XCAFDoc_ColorSurf`, `XCAFPrs_DocumentExplorer::mergedStyle`
    /// `XCAFPrs_DocumentExplorer.cxx:28-81`). This is the document-node style
    /// the OBJ writer receives as `aDocNode.Style`.
    pub color: Option<occt_core::quantity::Color>,
    pub layer: String,
    /// Index into [`BRepModel::materials`], i.e. the `VisMaterialRefGUID` tree
    /// node link `XCAFDoc_VisMaterialTool::SetShapeMaterial` creates
    /// (`XCAFDoc_VisMaterialTool.cxx:131-147`). `None` means no link.
    pub material: Option<usize>,
    /// `XCAFPrs::CollectStyleSettings` output for this node
    /// (`XCAFPrs.cxx:118-260`): the `(subshape, style)` settings carried by the
    /// node label and its `AddSubShape` sub-shape labels, in collection order.
    /// The OBJ writer turns it into the per-face style map with
    /// `RWMesh_ShapeIterator::dispatchStyles`
    /// (`rwmesh::obj_caf_writer`). Empty when the STEP file carries no style.
    pub style_settings: Vec<(TopoShape, XcafPrsStyle)>,
    /// `TDataStd_NamedData` attached to this node's label by
    /// `setShapeNamedData` (`RWGltf_CafReader.cxx:658`). The glTF reader fills
    /// it from the node's or mesh's `extras` object.
    pub named_data: Option<crate::tdata_std_named_data::TDataStdNamedData>,
}

/// The role a `TDF_Label` plays in the `XCAFDoc_ShapeTool` shape tree. These
/// are the four kinds of label `RWMesh_CafReader::addShapeIntoDoc`
/// (`RWMesh_CafReader.cxx:283-511`) creates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LabelRole {
    /// `XCAFDoc_ShapeTool::IsSimpleShape`: the label holds a plain shape
    /// (`AddShape(theShape, false)`).
    Shape,
    /// `XCAFDoc_ShapeTool::IsAssembly`: the label holds a compound whose
    /// children are reference labels (`AddShape(compound, true)`).
    Assembly,
    /// `XCAFDoc_ShapeTool::IsReference`: an instance label inside an assembly
    /// (`AddComponent`); [`ModelLabel::referred`] is its Product label.
    Reference,
    /// `XCAFDoc_ShapeTool::IsSubShape`: a face/edge/vertex label under a
    /// Product (`AddSubShape`), created only for sub-shapes carrying their own
    /// attributes (`RWMesh_CafReader.cxx:458-511`).
    SubShape,
}

/// One `TDF_Label` of the `XCAFDoc_ShapeTool` shape tree.
///
/// OCCT keys every document attribute by label; this crate has no `TDocStd`
/// document, so a label is addressed by its index in [`BRepModel::labels`].
/// The attributes below are exactly the ones the OBJ writer reads:
/// `TDataStd_Name` (`setShapeName`), `XCAFDoc_ColorTool` (`setShapeStyle`),
/// the `VisMaterialRefGUID` link and `TDataStd_NamedData`
/// (`setShapeNamedData`).
#[derive(Debug, Clone)]
pub struct ModelLabel {
    pub name: String,
    pub role: LabelRole,
    /// The label's shape (`XCAFDoc_ShapeTool::GetShape`). `None` for an
    /// assembly label whose compound has not been rebuilt and for a reference
    /// label that stores no shape of its own.
    pub shape: Option<TopoShape>,
    /// `XCAFDoc_ShapeTool::GetReferredShape`: the Product label of a reference
    /// label; `None` for every other role.
    pub referred: Option<usize>,
    pub parent: Option<usize>,
    pub children: Vec<usize>,
    /// `XCAFDoc_ColorTool::SetColor(..., XCAFDoc_ColorSurf)`.
    pub color: Option<occt_core::quantity::Color>,
    /// `XCAFDoc_ColorTool::SetColor(..., XCAFDoc_ColorCurv)`: the curve colour
    /// `setShapeStyle` stores for a line style
    /// (`RWMesh_CafReader.cxx:229-232`). Only edge labels carry it on the glTF
    /// path (`RWGltf_GltfJsonParser.cxx:2020-2023`); the OBJ writer reads the
    /// surface colour.
    pub color_curv: Option<occt_core::quantity::Color>,
    /// `XCAFDoc_VisMaterialTool::SetShapeMaterial` link (index into
    /// [`BRepModel::materials`]).
    pub material: Option<usize>,
    /// `TDataStd_NamedData` attached to this label.
    pub named_data: Option<crate::tdata_std_named_data::TDataStdNamedData>,
}

/// BRepModel: a set of named shapes with bounding box.
#[derive(Debug, Clone, Default)]
pub struct BRepModel {
    pub shapes: Vec<ModelShape>,
    /// Material table: the `XCAFDoc_VisMaterialTool` child labels
    /// (`XCAFDoc_VisMaterialTool.cxx:77-127`). An entry becomes `None` when
    /// `RemoveMaterial` forgets the label's attributes; the slot stays so the
    /// indices held by `ModelShape::material` remain label identities.
    pub materials: Vec<Option<XcafDocVisMaterial>>,
    /// `XCAFDoc_ShapeTool` label tree. Empty for models that were not filled
    /// from a document-style reader (`RWGltf_CafReader::fillDocument`,
    /// `RWGltf_CafReader.cxx:435-468`); [`BRepModel::shapes`] stays the flat
    /// geometry list either way.
    pub labels: Vec<ModelLabel>,
    /// `XCAFDoc_ShapeTool::GetFreeShapes` (`cxx:733-746`): the top-level labels
    /// in creation order.
    pub roots: Vec<usize>,
    /// `XCAFDoc_LengthUnit` on the document's main label root
    /// (`XCAFDoc_DocumentTool.cxx:374-408`): metres per model unit, i.e. the
    /// third `GetLengthUnit` argument's default (meter) scale. `STEPCAFControl_Reader::prepareUnits`
    /// (`cxx:508-522`) stores the file's cascade unit as
    /// `SetLengthUnit(doc, cascadeUnitMM, Millimeter)`, so a millimetre file
    /// records `0.001`.
    ///
    /// No reader of the port's geometry applies this scale — model coordinates
    /// stay in file units, exactly like an OCCT document — it only feeds
    /// consumers that need a physical unit, `RWGltf_CafWriter::Perform`'s
    /// `myCSTrsf.SetInputLengthUnit` (`RWGltf_CafWriter.cxx:645-649`) being the
    /// first. `None` is a document without the attribute, where
    /// `GetLengthUnit` returns FALSE and the converter stays unscaled.
    pub length_unit: Option<f64>,
    by_name: HashMap<String, usize>,
}

impl BRepModel {
    pub fn new() -> Self { Self::default() }

    pub fn add(&mut self, name: &str, shape: TopoShape) -> usize {
        let idx = self.shapes.len();
        self.shapes.push(ModelShape {
            name: name.into(),
            shape,
            color: None,
            layer: String::new(),
            material: None,
            style_settings: Vec::new(),
            named_data: None,
        });
        self.by_name.insert(name.into(), idx);
        idx
    }

    pub fn add_with_color(&mut self, name: &str, shape: TopoShape, color: occt_core::quantity::Color) -> usize {
        let idx = self.add(name, shape);
        self.shapes[idx].color = Some(color);
        idx
    }

    pub fn find(&self, name: &str) -> Option<&ModelShape> {
        self.by_name.get(name).map(|&i| &self.shapes[i])
    }

    pub fn find_mut(&mut self, name: &str) -> Option<&mut ModelShape> {
        let i = *self.by_name.get(name)?;
        self.shapes.get_mut(i)
    }

    pub fn get(&self, i: usize) -> Option<&ModelShape> { self.shapes.get(i) }

    pub fn remove(&mut self, name: &str) -> Option<TopoShape> {
        let idx = self.by_name.remove(name)?;
        // Compact: remove at idx, fix map indices
        let removed = self.shapes.remove(idx);
        for (_n, i) in self.by_name.iter_mut() {
            if *i > idx { *i -= 1; }
        }
        Some(removed.shape)
    }

    pub fn len(&self) -> usize { self.shapes.len() }
    pub fn is_empty(&self) -> bool { self.shapes.is_empty() }

    // -----------------------------------------------------------------
    // XCAFDoc_VisMaterialTool (`XCAFDoc_VisMaterialTool.cxx:77-236`)
    // -----------------------------------------------------------------
    //
    // OCCT keys the material table and every shape -> material link by
    // `TDF_Label`. This crate has no document, so a material is addressed by its
    // index in `materials` (the label's position in the tool's children) and a
    // shape by its name (the node label identity `by_name` already tracks).

    /// `XCAFDoc_VisMaterialTool::AddMaterial(theMat, theName)`
    /// (`cxx:77-88`): append the material and, when `theName` is non-empty, set
    /// the `TDataStd_Name` that `RWMesh_MaterialMap::AddMaterial` later reads as
    /// the `newmtl` key (`RWMesh_MaterialMap.cxx:67-77`).
    pub fn add_material(&mut self, mut material: XcafDocVisMaterial, name: &str) -> usize {
        // The attribute is attached to a *new* child label, so the name of that
        // label is whatever `theName` says — an empty name means the new label
        // carries no `TDataStd_Name` at all.
        material.set_label_name(if name.is_empty() { None } else { Some(name.to_string()) });
        self.materials.push(Some(material));
        self.materials.len() - 1
    }

    /// `XCAFDoc_VisMaterialTool::AddMaterial(theName)` (`cxx:92-103`): a fresh
    /// empty `XCAFDoc_VisMaterial` under a new child label.
    pub fn add_new_material(&mut self, name: &str) -> usize {
        self.add_material(XcafDocVisMaterial::new(), name)
    }

    /// `XCAFDoc_VisMaterialTool::RemoveMaterial` (`cxx:107-110`):
    /// `theLabel.ForgetAllAttributes(true)` — the label stays, its attributes go,
    /// so `GetMaterial` and `GetMaterials` stop reporting it while existing tree
    /// links keep pointing at the (now empty) label.
    pub fn remove_material(&mut self, index: usize) {
        if let Some(slot) = self.materials.get_mut(index) {
            *slot = None;
        }
    }

    /// `XCAFDoc_VisMaterialTool::GetMaterial` (`cxx:55-56`): null unless the
    /// label is in the material table.
    pub fn get_material(&self, index: usize) -> Option<&XcafDocVisMaterial> {
        self.materials.get(index).and_then(|m| m.as_ref())
    }

    /// `XCAFDoc_VisMaterialTool::IsMaterial` (`hxx:58-61`).
    pub fn is_material(&self, index: usize) -> bool {
        self.get_material(index).is_some()
    }

    /// `XCAFDoc_VisMaterialTool::GetMaterials` (`cxx:114-127`): the labels still
    /// carrying an `XCAFDoc_VisMaterial`, in table order.
    pub fn material_count(&self) -> usize {
        self.materials.iter().filter(|m| m.is_some()).count()
    }

    /// `XCAFDoc_VisMaterialTool::SetShapeMaterial(theShapeLabel, theMaterialLabel)`
    /// (`cxx:131-147`): a null material label forgets the link
    /// (`cxx:134-138`); otherwise the link is (re)pointed at `material`.
    /// Returns FALSE when `shape_name` names no shape, like the `TopoDS_Shape`
    /// overload (`cxx:191-202`).
    pub fn set_shape_material(&mut self, shape_name: &str, material: Option<usize>) -> bool {
        let Some(&index) = self.by_name.get(shape_name) else {
            return false;
        };
        self.shapes[index].material = material;
        true
    }

    /// `XCAFDoc_VisMaterialTool::UnSetShapeMaterial` (`cxx:151-154`).
    pub fn unset_shape_material(&mut self, shape_name: &str) -> bool {
        let Some(&index) = self.by_name.get(shape_name) else {
            return false;
        };
        self.shapes[index].material = None;
        true
    }

    /// `XCAFDoc_VisMaterialTool::IsSetShapeMaterial` (`cxx:158-162`): TRUE when
    /// the link's tree node has a father, i.e. when a material label was set —
    /// even if that label has since been emptied by `RemoveMaterial`.
    pub fn is_set_shape_material(&self, shape_name: &str) -> bool {
        self.find(shape_name).map(|s| s.material.is_some()).unwrap_or(false)
    }

    /// `XCAFDoc_VisMaterialTool::GetShapeMaterial` (`cxx:181-187`): the material
    /// of the linked label, or NULL when no link is set or the linked label no
    /// longer carries a material.
    pub fn get_shape_material(&self, shape_name: &str) -> Option<&XcafDocVisMaterial> {
        let index = self.find(shape_name)?.material?;
        self.get_material(index)
    }

    // -----------------------------------------------------------------
    // XCAFDoc_ShapeTool shape tree
    // (`XCAFDoc_ShapeTool.cxx:633-1130`, `:2137-2270`)
    // -----------------------------------------------------------------
    //
    // The tree is the `TDF_Label` hierarchy the document-style readers fill and
    // `XCAFPrs_DocumentExplorer` walks. A label owns either a shape
    // (`Shape`/`SubShape`), a compound rebuilt from its components
    // (`Assembly`), or nothing but a link to the Product it instances
    // (`Reference`).

    /// `XCAFDoc_ShapeTool::AddShape(theShape, theMakeAssembly)`
    /// (`cxx:633-652` + `addShape`): a new free label holding `shape`. An
    /// assembly label starts with the empty compound it is given here; its
    /// content is rebuilt by [`BRepModel::update_assemblies`].
    pub fn add_shape_label(&mut self, shape: TopoShape, make_assembly: bool) -> usize {
        let idx = self.new_shape_label(shape, make_assembly);
        self.roots.push(idx);
        idx
    }

    /// `XCAFDoc_ShapeTool::addShape(theLabel, ...)` (`cxx:496-540`): a shape
    /// label that is not a free shape, i.e. one a component refers to. `AddShape`
    /// publishes a label through `GetFreeShapes` only while no reference points at
    /// it, so a product added here is reachable through its reference only.
    fn new_shape_label(&mut self, shape: TopoShape, make_assembly: bool) -> usize {
        let idx = self.labels.len();
        self.labels.push(ModelLabel {
            name: String::new(),
            role: if make_assembly { LabelRole::Assembly } else { LabelRole::Shape },
            shape: Some(shape),
            referred: None,
            parent: None,
            children: Vec::new(),
            color: None,
            color_curv: None,
            material: None,
            named_data: None,
        });
        idx
    }

    /// `XCAFDoc_ShapeTool::AddComponent(theAssembly, theCompLabel, theLoc)`
    /// (`cxx:954-990`): a reference label under `assembly` pointing at the
    /// Product label `referred`. A simple shape promoted to an assembly keeps
    /// its shape and becomes an [`LabelRole::Assembly`].
    pub fn add_component(&mut self, assembly: usize, referred: usize) -> Option<usize> {
        // `cxx:959-972`: only an assembly (or a simple shape promoted to one)
        // can take components.
        match self.labels.get(assembly).map(|l| l.role) {
            Some(LabelRole::Assembly) => {}
            Some(LabelRole::Shape) => {
                if let Some(l) = self.labels.get_mut(assembly) {
                    l.role = LabelRole::Assembly;
                }
            }
            _ => return None,
        }
        if self.labels.get(referred).is_none() {
            return None;
        }
        Some(self.new_reference(Some(assembly), referred, None))
    }

    /// `XCAFDoc_ShapeTool::AddComponent(theAssembly, theShape, theMakeAssembly)`
    /// (`cxx:994-1007`): `AddShape` for the shape *without* its location, then a
    /// reference that carries the location. In this crate a shape's placement is
    /// baked into its geometry instead of a `TopLoc_Location`, so `unlocated` is
    /// the Product's definition and `placed` the instance's geometry; the
    /// reference label keeps the latter (a reference normally stores no shape of
    /// its own). Returns `(reference_label, product_label)`.
    pub fn add_component_shape(
        &mut self,
        assembly: usize,
        unlocated: TopoShape,
        placed: TopoShape,
        make_assembly: bool,
    ) -> Option<(usize, usize)> {
        let product = self.new_shape_label(unlocated, make_assembly);
        let label = match self.add_component(assembly, product) {
            Some(l) => l,
            None => {
                self.labels.remove(product);
                return None;
            }
        };
        self.labels[label].shape = Some(placed);
        Some((label, product))
    }

    /// `XCAFDoc_ShapeTool::AddShape(theShape, theMakeAssembly)`
    /// (`cxx:633-652` + `addShape` `cxx:496-540`): a free label for `shape`. A
    /// shape carrying a location is turned into a reference to a label holding
    /// the same shape without the location (`cxx:510-519`), so the free label is
    /// the reference and the Product label comes second. Without a location both
    /// are the same label (when `!split`).
    ///
    /// `split` is the port's stand-in for `!theShape.Location().IsIdentity()`.
    pub fn add_free_shape(
        &mut self,
        unlocated: TopoShape,
        placed: TopoShape,
        make_assembly: bool,
        split: bool,
    ) -> (usize, usize) {
        let product = self.new_shape_label(unlocated, make_assembly);
        if !split {
            self.roots.push(product);
            return (product, product);
        }
        let label = self.new_reference(None, product, Some(placed));
        self.roots.push(label);
        (label, product)
    }

    /// `XCAFDoc_ShapeTool::MakeReference` (`cxx:442-460`): a reference label
    /// under `parent` (a free label when `None`) pointing at `referred`.
    fn new_reference(&mut self, parent: Option<usize>, referred: usize, placed: Option<TopoShape>) -> usize {
        let idx = self.labels.len();
        self.labels.push(ModelLabel {
            name: String::new(),
            role: LabelRole::Reference,
            shape: placed,
            referred: Some(referred),
            parent,
            children: Vec::new(),
            color: None,
            color_curv: None,
            material: None,
            named_data: None,
        });
        if let Some(p) = parent {
            self.labels[p].children.push(idx);
        }
        idx
    }

    /// `XCAFDoc_ShapeTool::AddSubShape(theParentLabel, theShape)`
    /// (`cxx:1129-1150`): a sub-shape label under `parent`.
    pub fn add_sub_shape_label(&mut self, parent: usize, shape: TopoShape) -> Option<usize> {
        if !matches!(self.labels.get(parent).map(|l| l.role), Some(LabelRole::Shape) | Some(LabelRole::Assembly)) {
            return None;
        }
        let idx = self.labels.len();
        self.labels.push(ModelLabel {
            name: String::new(),
            role: LabelRole::SubShape,
            shape: Some(shape),
            referred: None,
            parent: Some(parent),
            children: Vec::new(),
            color: None,
            color_curv: None,
            material: None,
            named_data: None,
        });
        self.labels[parent].children.push(idx);
        Some(idx)
    }

    /// `XCAFDoc_ShapeTool::IsAssembly` (`hxx`): the label carries the
    /// `AssemblyGUID` attribute.
    pub fn is_assembly(&self, label: usize) -> bool {
        matches!(self.labels.get(label).map(|l| l.role), Some(LabelRole::Assembly))
    }

    /// `XCAFDoc_ShapeTool::IsReference` (`hxx`).
    pub fn is_reference(&self, label: usize) -> bool {
        matches!(self.labels.get(label).map(|l| l.role), Some(LabelRole::Reference))
    }

    /// `XCAFDoc_ShapeTool::GetReferredShape(theLabel, theRefLabel)`
    /// (`cxx:936-950`): the Product label a reference label points at. A
    /// non-reference label refers to itself.
    pub fn get_referred_shape(&self, label: usize) -> Option<usize> {
        let l = self.labels.get(label)?;
        Some(l.referred.unwrap_or(label))
    }

    /// `XCAFDoc_ShapeTool::GetReferredShape` exactly (`cxx:939-950`): `Some`
    /// only when `label` is a reference label, `None` otherwise.
    pub fn referred_label(&self, label: usize) -> Option<usize> {
        match self.labels.get(label) {
            Some(l) if l.role == LabelRole::Reference => l.referred,
            _ => None,
        }
    }

    /// The label's shape (`XCAFDoc_ShapeTool::GetShape`): its own shape, else
    /// the shape of the Product it refers to, else the compound rebuilt from
    /// its components.
    pub fn label_shape(&self, label: usize) -> Option<TopoShape> {
        let l = self.labels.get(label)?;
        if let Some(s) = &l.shape {
            if l.role != LabelRole::Assembly || !l.children.is_empty() {
                return Some(s.clone());
            }
        }
        if let Some(r) = l.referred {
            return self.label_shape(r);
        }
        None
    }

    pub fn label(&self, i: usize) -> Option<&ModelLabel> { self.labels.get(i) }
    pub fn label_mut(&mut self, i: usize) -> Option<&mut ModelLabel> { self.labels.get_mut(i) }

    /// `readNameAttribute` (`RWObj_CafWriter.cxx:60-71`): the label's
    /// `TDataStd_Name`, empty when unnamed.
    pub fn label_name(&self, label: usize) -> &str {
        self.labels.get(label).map(|l| l.name.as_str()).unwrap_or("")
    }

    /// `TDataStd_Name::Set(theLabel, theName)` followed by the reader's global
    /// "replace empty names" pass, which only fills labels that stayed unnamed.
    pub fn set_label_name(&mut self, label: usize, name: &str) {
        if let Some(l) = self.labels.get_mut(label) {
            l.name = name.to_string();
        }
    }

    /// `XCAFDoc_ShapeTool::UpdateAssemblies` (`cxx:1021-1040`): rebuild every
    /// root assembly compound from its components (`updateComponent`,
    /// `cxx:2137-2270`).
    pub fn update_assemblies(&mut self) {
        let roots = self.roots.clone();
        let mut updated = std::collections::HashSet::new();
        for root in roots {
            let Some(ref_label) = self.get_referred_shape(root) else { continue };
            self.update_component(ref_label, &mut updated);
        }
    }

    /// `XCAFDoc_ShapeTool::updateComponent` (`cxx:2137-2270`), restricted to
    /// the case this crate creates: an assembly label starts with the empty
    /// compound from `AddShape` (`cxx:2166-2174`, "a component may happen to be
    /// removed" — the topological child count never matches), so the compound
    /// is always rebuilt from the referred component shapes. The `isModified`
    /// bookkeeping only avoids rewriting an already up-to-date compound.
    fn update_component(
        &mut self,
        label: usize,
        updated: &mut std::collections::HashSet<usize>,
    ) -> Option<TopoShape> {
        if !self.is_assembly(label) {
            return None;
        }
        if !updated.insert(label) {
            return self.labels.get(label).and_then(|l| l.shape.clone());
        }
        let components = self.labels.get(label).map(|l| l.children.clone()).unwrap_or_default();
        let mut parts = Vec::with_capacity(components.len());
        for comp in components {
            // `cxx:2191-2204`: the component shape comes from the instance
            // label (`GetShape(aComponentLab)`), located by the component
            // placement. This crate's shapes carry no `TopLoc_Location`, so the
            // placement is already baked into the shape.
            let mut shape = self.labels.get(comp).and_then(|l| l.shape.clone());
            if shape.is_none() {
                if let Some(referred) = self.labels.get(comp).and_then(|l| l.referred) {
                    self.update_component(referred, updated);
                    shape = self.labels.get(referred).and_then(|l| l.shape.clone());
                }
            }
            if let Some(s) = shape {
                parts.push(s);
            }
        }
        let compound = crate::builder::TopoBuilder::new().make_compound_of(&parts).0;
        if let Some(l) = self.labels.get_mut(label) {
            l.shape = Some(compound.clone());
        }
        Some(compound)
    }

    pub fn names(&self) -> Vec<&str> { self.shapes.iter().map(|s| s.name.as_str()).collect() }

    /// Compute the model bounding box (union of per-shape transforms applied).
    /// children_fn: provides sub-shapes for geometry access (unused for bbox here).
    pub fn bounding_box(&self) -> BndBox {
        let b = BndBox::new();
        // Shapes have no embedded geometry in this port; bbox is empty unless
        // sub-shapes carry coordinates. Keep as a stub returning void box.
        for s in &self.shapes {
            if s.shape.shape_type() == ShapeType::Vertex {
                // No point data stored — skip.
            }
        }
        b
    }

    /// Filter shapes by type.
    pub fn by_type(&self, t: ShapeType) -> Vec<&ModelShape> {
        self.shapes.iter().filter(|s| s.shape.shape_type() == t).collect()
    }

    /// Count shapes by type.
    pub fn type_counts(&self) -> HashMap<ShapeType, usize> {
        let mut m = HashMap::new();
        for s in &self.shapes { *m.entry(s.shape.shape_type()).or_insert(0) += 1; }
        m
    }

    /// Apply a transform to all shapes.
    pub fn move_all(&mut self, t: &occt_core::gp::GpTrsf) {
        for s in &mut self.shapes {
            s.shape = crate::transform::move_shape(&s.shape, t);
        }
    }

    /// Rename a shape.
    pub fn rename(&mut self, old: &str, new: &str) -> bool {
        let idx = match self.by_name.remove(old) { Some(i) => i, None => return false };
        self.shapes[idx].name = new.into();
        self.by_name.insert(new.into(), idx);
        true
    }

    /// Collect all shapes into a single compound (top-level).
    pub fn as_compound(&self) -> TopoShape {
        let comp = TopoShape::new(ShapeType::Compound);
        // In this port compound children are not stored on TShape; return bare compound.
        comp
    }

    /// Merge another model into this one (name collision → keep both, suffix _2).
    pub fn merge(&mut self, other: &BRepModel) {
        for s in &other.shapes {
            let mut name = s.name.clone();
            let mut n = 2;
            while self.by_name.contains_key(&name) {
                name = format!("{}_{n}", s.name);
                n += 1;
            }
            self.add(&name, s.shape.clone());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_find_remove() {
        let mut m = BRepModel::new();
        let v = TopoShape::new(ShapeType::Vertex);
        m.add("v1", v.clone());
        assert!(m.find("v1").is_some());
        assert_eq!(m.len(), 1);
        let got = m.remove("v1").unwrap();
        assert!(got.is_vertex());
        assert!(m.find("v1").is_none());
    }

    #[test]
    fn rename_works() {
        let mut m = BRepModel::new();
        m.add("a", TopoShape::new(ShapeType::Face));
        assert!(m.rename("a", "b"));
        assert!(m.find("b").is_some());
        assert!(m.find("a").is_none());
    }

    #[test]
    fn merge_suffixes() {
        let mut m1 = BRepModel::new();
        m1.add("s", TopoShape::new(ShapeType::Vertex));
        let mut m2 = BRepModel::new();
        m2.add("s", TopoShape::new(ShapeType::Edge));
        m1.merge(&m2);
        assert_eq!(m1.len(), 2);
        assert!(m1.find("s").is_some());
        assert!(m1.find("s_2").is_some());
    }

    #[test]
    fn type_counts() {
        let mut m = BRepModel::new();
        m.add("e1", TopoShape::new(ShapeType::Edge));
        m.add("e2", TopoShape::new(ShapeType::Edge));
        m.add("f1", TopoShape::new(ShapeType::Face));
        let c = m.type_counts();
        assert_eq!(c[&ShapeType::Edge], 2);
        assert_eq!(c[&ShapeType::Face], 1);
    }
}
