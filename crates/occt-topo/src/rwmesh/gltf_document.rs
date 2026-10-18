//! glTF document fill: the `XCAFDoc_ShapeTool` label tree
//! `RWMesh_CafReader::addShapeIntoDoc` builds.
//!
//! Source:
//! * `RWGltf_GltfJsonParser::bindNamedShape` (`RWGltf_GltfJsonParser.cxx:2372-2506`)
//!   — the per-shape attributes `myAttribMap` collects while the graph is parsed;
//!   the node/mesh/primitive-array bind order is what separates a Product name
//!   from an Instance name.
//! * `RWGltf_CafReader::fillDocument` / `addShapeIntoDoc`
//!   (`RWGltf_CafReader.cxx:435-467`, `:470-688`) — the replay of those
//!   attributes into labels.
//! * `RWMesh_CafReader::setShapeName` / `setShapeStyle` / `setShapeNamedData` /
//!   `addSubShapeIntoDoc` (`RWMesh_CafReader.cxx:199-277`, `:456-511`).
//!
//! Two label kinds appear per placed mesh: the Product (the mesh definition,
//! shared by every node that instances it) and, under the assembly, a reference
//! label holding the placement and the Instance attributes. OCCT keys both the
//! attribute map and `CafDocumentTools::ComponentMap` by `TopoDS_Shape`, which
//! includes the `TopLoc_Location`; this crate bakes placements into the geometry
//! (the mesh pipeline does not apply shape locations), so:
//!
//! * the Product label keeps the *unplaced* shape (the primitive array's own
//!   geometry, `gltfParsePrimArray`'s shape at the identity location);
//! * the reference label keeps the placed copy, which stands in for the
//!   `XCAFDoc_Location` attribute the C++ stores on it. `updateComponent`
//!   (`XCAFDoc_ShapeTool.cxx:2137-2270`) therefore rebuilds an assembly from the
//!   placed components, the same compound the C++ builds from located ones.
//!
//! Not ported here, and why:
//! * `XCAFDoc_ShapeTool::prepareAssembly` (`cxx:590-626`): it gives every
//!   component of an assembly compound a non-identity location so `IsSame` keys
//!   stay apart. Locations are baked here, so every component already has its own
//!   shape.
//! * `AddShape`'s `FindShape` dedup (`cxx:501-505`) and `CafDocumentTools::
//!   OriginalShapeMap` (`RWMesh_CafReader.cxx:366-373`, write-only in the C++
//!   reader). The reader's own `ComponentMap` (ported below) already reuses the
//!   Product of every repeated mesh placement.
//! * `myShapeScaleMap` (`RWGltf_GltfJsonParser.cxx:2384-2400`, `:489-519`): a
//!   workaround for one shared face receiving different non-uniform scales. Each
//!   placement builds its own shape here, and `myToApplyScale` already scales the
//!   triangulation at parse time.

use std::collections::HashMap;

use occt_core::gp::{GpTrsf, TrsfForm};
use occt_core::quantity::Color;

use crate::abs::ShapeType;
use crate::builder::TopoBuilder;
use crate::model::BRepModel;
use crate::shape::TopoShape;
use crate::tdata_std_named_data::TDataStdNamedData;
use crate::tgeometry::GeometryRegistry;
use crate::xcaf_prs::XcafPrsStyle;

use super::coordinate_system::RwMeshCoordinateSystemConverter;
use super::gltf::{inst_to_local_node, inst_to_node, GltfNodeInst, RwGltfPrimitiveMode};
use super::gltf_shape_reader::{node_proxy_shape, prim_proxy_shape};
use super::gltf_types::{GlMesh, GlNode, GltfCtx};

/// `RWMesh_NodeAttributes` (`RWMesh_Reader.hxx:40-62`): the attributes
/// `bindNamedShape` binds to one shape occurrence.
#[derive(Clone, Default)]
struct NodeAttribs {
    /// `Name` — the `TDataStd_Name` the label ends up with.
    name: String,
    /// `RawName` — the primitive array's mesh name, the value an unnamed node
    /// copies (`RWGltf_GltfJsonParser.cxx:2444-2447`).
    raw_name: String,
    style: XcafPrsStyle,
    named_data: Option<TDataStdNamedData>,
}

/// One shape `addShapeIntoDoc` receives, plus the attributes bound to it.
struct DocShape {
    /// `theShape`: the placed shape. The node location is baked into the
    /// geometry here instead of living in a `TopLoc_Location`.
    shape: TopoShape,
    /// `theShape.Located(TopLoc_Location())` — OCCT's `aShapeNoLoc`, the Product
    /// key of `myAttribMap` and of `CafDocumentTools::ComponentMap`.
    unlocated: TopoShape,
    /// `!theShape.Location().IsIdentity()`: the Instance attributes are keyed
    /// apart from the Product ones. `bindNamedShape` also forces this apart with
    /// a fake identity location when a named node sits at the identity location
    /// (`RWGltf_GltfJsonParser.cxx:2449-2454`).
    split: bool,
    /// `myAttribMap[theShape]`.
    instance: Option<NodeAttribs>,
    /// `myAttribMap[aShapeNoLoc]`.
    product: Option<NodeAttribs>,
    children: Vec<DocShape>,
    /// `RWGltf_GltfLatePrimitiveArray::Name()` of a primitive-array shape: what
    /// an unnamed node copies instead of a mesh name
    /// (`RWGltf_GltfJsonParser.cxx:2461-2487`).
    late_name: String,
    /// This shape is a primitive array, so `late_name` takes part in that scan.
    has_late: bool,
}

impl DocShape {
    fn leaf(shape: TopoShape, unlocated: TopoShape) -> Self {
        Self {
            shape,
            unlocated,
            split: false,
            instance: None,
            product: None,
            children: Vec::new(),
            late_name: String::new(),
            has_late: false,
        }
    }
}

/// `RWMesh_CafReader::CafDocumentTools` (`RWMesh_CafReader.hxx:57-64`), reduced
/// to the state that is not a label attribute.
struct DocTools {
    /// `ComponentMap`: an unlocated shape already in the document, and the
    /// Product label it was added as.
    component_map: HashMap<usize, usize>,
    /// `XCAFDoc_VisMaterialTool` identity: `XCAFDoc_VisMaterial` handle -> entry
    /// in [`BRepModel::materials`]. Shares the table with the flat shape list.
    materials: HashMap<u64, usize>,
}

// ---------------------------------------------------------------------------
// Document fill
// ---------------------------------------------------------------------------

/// Fill `model` from a parsed glTF scene: the flat shape list (one entry per
/// placed primitive array, the `MeshScene` nodes) followed by the label tree of
/// `RWGltf_CafReader::fillDocument` (`cxx:435-467`).
pub(super) fn fill_document(
    model: &mut BRepModel,
    nodes: &[GlNode],
    meshes: &[GlMesh],
    scene_roots: &[Vec<usize>],
    insts: &[GltfNodeInst],
    ctx: &GltfCtx<'_>,
    skip_empty_nodes: bool,
) {
    let mut tools = DocTools { component_map: HashMap::new(), materials: HashMap::new() };
    let mut prim_inst: HashMap<(usize, usize, usize), usize> = HashMap::new();
    let mut prim_count: HashMap<(usize, usize), usize> = HashMap::new();

    for (i, inst) in insts.iter().enumerate() {
        // One flat shape per placed array, named and materialized like the
        // `MeshScene` node it becomes in `read_gltf_scene`.
        let node = inst_to_node(inst);
        let mut name = inst.name.clone();
        let mut n = 2;
        while model.find(&name).is_some() {
            name = format!("{}_{n}", inst.name);
            n += 1;
        }
        let idx = model.add(&name, node_proxy_shape(&node));
        if let Some(nd) = &inst.named_data {
            model.shapes[idx].named_data = Some(nd.clone());
        }
        if let Some(vis) = &inst.vis_material {
            let mi = *tools
                .materials
                .entry(vis.id())
                .or_insert_with(|| model.add_material(vis.clone(), vis.raw_name().unwrap_or("")));
            model.shapes[idx].material = Some(mi);
        }

        let prim = *prim_count.entry((inst.node, inst.mesh)).or_insert(0);
        prim_inst.insert((inst.node, inst.mesh, prim), i);
        prim_count.insert((inst.node, inst.mesh), prim + 1);
    }

    let mut builder = Builder {
        nodes,
        meshes,
        insts,
        prim_inst,
        unplaced_prim: HashMap::new(),
        visited: vec![false; nodes.len()],
        is_gltf1: ctx.is_gltf1,
        use_mesh_name_as_fallback: ctx.use_mesh_name_as_fallback,
        skip_empty_nodes,
        converter: ctx.converter,
    };
    for roots in scene_roots {
        for &r in roots {
            // `gltfParseSceneNodes` (`cxx:1536-1542`): a reference to a missing
            // node stops the remaining roots of that scene.
            if r >= nodes.len() {
                break;
            }
            let Some(doc) = builder.build_node(r, &GpTrsf::identity()) else { continue };
            add_shape_into_doc(model, &mut tools, &doc, None, "");
        }
    }
    // `fillDocument` (`cxx:464`).
    model.update_assemblies();
}

// ---------------------------------------------------------------------------
// Shape tree building (gltfParseMesh / gltfParseSceneNode)
// ---------------------------------------------------------------------------

struct Builder<'a> {
    nodes: &'a [GlNode],
    meshes: &'a [GlMesh],
    insts: &'a [GltfNodeInst],
    /// `(node, mesh, primitive)` -> index in `insts`.
    prim_inst: HashMap<(usize, usize, usize), usize>,
    /// The unplaced shape of a primitive array — `myShapeMap[ShapeMapGroup_
    /// PrimArray]` (`RWGltf_GltfJsonParser.cxx:1847-1877`), shared by every node
    /// that places it.
    unplaced_prim: HashMap<(usize, usize), TopoShape>,
    /// `findNodeShape` (`cxx:1573-1576`): a node is turned into a shape once.
    visited: Vec<bool>,
    is_gltf1: bool,
    use_mesh_name_as_fallback: bool,
    skip_empty_nodes: bool,
    converter: &'a RwMeshCoordinateSystemConverter,
}

impl Builder<'_> {
    /// `gltfParseSceneNode` (`RWGltf_GltfJsonParser.cxx:1566-1717`).
    fn build_node(&mut self, idx: usize, world: &GpTrsf) -> Option<DocShape> {
        if self.visited[idx] {
            return None;
        }
        self.visited[idx] = true;
        let n = &self.nodes[idx];
        // `cxx:493`, `:637`: the local location is conjugated by the converter.
        let mut local = n.trsf.clone();
        if !self.converter.is_empty() {
            self.converter.transform_transformation(&mut local);
        }
        let world = world.multiplied(&local);
        let local_identity = local.form() == TrsfForm::Identity;

        let mut sub: Vec<DocShape> = Vec::new();
        for &c in &n.children {
            if let Some(cs) = self.build_node(c, &world) {
                sub.push(cs);
            }
        }
        let had_children = !sub.is_empty();
        for m in node_meshes(n) {
            if let Some(ms) = self.build_mesh(idx, m) {
                sub.push(ms);
            }
        }

        let mut doc = if sub.is_empty() {
            if self.skip_empty_nodes {
                // `gltfParseSceneNodes` (`cxx:1554-1559`): a node shape without
                // face, edge or vertex is dropped.
                return None;
            }
            DocShape::leaf(TopoShape::new(ShapeType::Compound), TopoShape::new(ShapeType::Compound))
        } else if !had_children && sub.len() == 1 {
            // `cxx:1705-1708`: a single mesh and no child nodes — the node shape
            // is that sub-shape.
            sub.pop().unwrap()
        } else {
            let parts: Vec<TopoShape> = sub.iter().map(|d| d.shape.clone()).collect();
            let unplaced: Vec<TopoShape> = sub.iter().map(|d| d.unlocated.clone()).collect();
            let mut compound = DocShape::leaf(
                TopoBuilder::new().make_compound_of(&parts).0,
                TopoBuilder::new().make_compound_of(&unplaced).0,
            );
            compound.children = sub;
            compound
        };
        bind_node(&mut doc, n, local_identity, self.is_gltf1, self.use_mesh_name_as_fallback);
        Some(doc)
    }

    /// `gltfParseMesh` (`RWGltf_GltfJsonParser.cxx:1720-1778`): the mesh shape is
    /// the single primitive-array shape when there is one, else a compound.
    fn build_mesh(&mut self, node: usize, mesh_idx: usize) -> Option<DocShape> {
        let mesh = self.meshes.get(mesh_idx)?;
        let mut prims: Vec<DocShape> = Vec::with_capacity(mesh.prims.len());
        for prim_idx in 0..mesh.prims.len() {
            if let Some(doc) = self.build_prim(node, mesh_idx, prim_idx) {
                prims.push(doc);
            }
        }
        if prims.is_empty() {
            // `cxx:1754-1771`: no primitive array produced a shape, so the mesh
            // shape stays null and is not added to the node.
            return None;
        }
        let mut doc = if prims.len() == 1 {
            prims.pop().unwrap()
        } else {
            let parts: Vec<TopoShape> = prims.iter().map(|d| d.shape.clone()).collect();
            let unplaced: Vec<TopoShape> = prims.iter().map(|d| d.unlocated.clone()).collect();
            let mut compound = DocShape::leaf(
                TopoBuilder::new().make_compound_of(&parts).0,
                TopoBuilder::new().make_compound_of(&unplaced).0,
            );
            compound.children = prims;
            compound
        };
        bind_mesh(&mut doc, mesh, self.is_gltf1, self.use_mesh_name_as_fallback);
        Some(doc)
    }

    /// `gltfParsePrimArray`'s shape (`RWGltf_GltfJsonParser.cxx:1969-2013`).
    fn build_prim(&mut self, node: usize, mesh_idx: usize, prim_idx: usize) -> Option<DocShape> {
        let inst_idx = *self.prim_inst.get(&(node, mesh_idx, prim_idx))?;
        let inst = self.insts.get(inst_idx)?;
        // `cxx:1960`: only an array that carries data produces a shape.
        if inst.late.data().is_empty() {
            return None;
        }
        let placed = prim_proxy_shape(&inst_to_node(inst));
        let key = (mesh_idx, prim_idx);
        let unplaced = match self.unplaced_prim.get(&key) {
            Some(s) => s.clone(),
            None => {
                let s = prim_proxy_shape(&inst_to_local_node(inst));
                self.unplaced_prim.insert(key, s.clone());
                s
            }
        };
        let mesh = self.meshes.get(mesh_idx)?;
        let mut doc = DocShape::leaf(placed, unplaced);
        doc.late_name = inst.late.name().to_string();
        doc.has_late = true;
        // `cxx:2012-2028`: the primitive array's attributes — a curve color for
        // lines, the material otherwise — and the mesh name as `RawName`.
        doc.product = prim_style(inst).map(|style| NodeAttribs {
            name: String::new(),
            raw_name: mesh.name.clone(),
            style,
            named_data: None,
        });
        Some(doc)
    }
}

/// The mesh ids a node references: glTF 2.0 `mesh`, then the glTF 1.0 `meshes`
/// array (`RWGltf_GltfJsonParser.cxx:1644-1707`).
fn node_meshes(n: &GlNode) -> Vec<usize> {
    let mut ids = Vec::with_capacity(n.meshes.len() + 1);
    if let Some(m) = n.mesh {
        ids.push(m);
    }
    ids.extend(n.meshes.iter().copied());
    ids
}

/// `gltfParsePrimArray`'s style: lines carry their base color as the curve
/// color, everything else binds the `XCAFDoc_VisMaterial`
/// (`RWGltf_GltfJsonParser.cxx:2014-2028`). `None` when the array has no style
/// (`RWGltf_GltfLatePrimitiveArray::HasStyle`, `hxx:57`).
fn prim_style(inst: &GltfNodeInst) -> Option<XcafPrsStyle> {
    let mut style = XcafPrsStyle::default();
    if inst.late.prim_mode() == RwGltfPrimitiveMode::Lines {
        let mat = inst.material.as_ref()?;
        style.set_color_curv(Color::new(
            mat.diffuse.0 as f32,
            mat.diffuse.1 as f32,
            mat.diffuse.2 as f32,
        ));
    } else {
        style.set_material(Some(inst.vis_material.clone()?));
    }
    Some(style)
}

// ---------------------------------------------------------------------------
// bindNamedShape
// ---------------------------------------------------------------------------

/// `bindMeshShape` -> `bindNamedShape(theShape, ShapeMapGroup_Meshes, identity,
/// ...)` (`RWGltf_GltfJsonParser.hxx:333-347`): the mesh's own name, extras and
/// the primitive-array style. The identity location means the Product and
/// Instance attributes are the same map entry.
fn bind_mesh(doc: &mut DocShape, mesh: &GlMesh, is_gltf1: bool, fallback: bool) {
    let mut att = NodeAttribs {
        // `cxx:2410-2418`: the user name, else the glTF 1.0 id.
        name: if !mesh.name.is_empty() {
            mesh.name.clone()
        } else if is_gltf1 {
            mesh.id.clone()
        } else {
            String::new()
        },
        raw_name: if is_gltf1 { mesh.id.clone() } else { String::new() },
        style: XcafPrsStyle::default(),
        named_data: mesh.extras.clone(),
    };
    match doc.shape.shape_type() {
        ShapeType::Face | ShapeType::Edge | ShapeType::Vertex => {
            if let Some(prim) = doc.product.clone() {
                // `cxx:2435-2442`: the shape keeps the style of the primitive
                // array bound to the same face.
                if let Some(m) = prim.style.material() {
                    att.style.set_material(Some(m.clone()));
                }
                if prim.style.is_set_color_curv() {
                    att.style.set_color_curv(prim.style.color_curv());
                }
                if prim.style.is_set_color_surf() {
                    att.style.set_color_surf(prim.style.color_surf_rgba());
                }
                // `cxx:2444-2447`: fall back to the mesh name the array carries.
                if att.name.is_empty() && fallback {
                    att.name = prim.raw_name.clone();
                }
            }
        }
        _ => {
            if att.name.is_empty() && fallback {
                // `cxx:2459-2489`: a compound takes the name its faces' arrays
                // agree on.
                if let Some(name) = common_late_name(doc) {
                    att.name = name;
                }
            }
        }
    }
    doc.product = Some(att.clone());
    doc.instance = Some(att);
    doc.split = false;
}

/// `bindNodeShape` -> `bindNamedShape(theShape, ShapeMapGroup_Nodes, ...)`
/// (`RWGltf_GltfJsonParser.hxx:314-330`): the node's name and extras, the style
/// inherited from the shape, and the identity-location trick that keeps the
/// Product name and the Instance name in separate map entries.
fn bind_node(doc: &mut DocShape, node: &GlNode, local_identity: bool, is_gltf1: bool, fallback: bool) {
    let mut att = NodeAttribs {
        name: if !node.name.is_empty() {
            node.name.clone()
        } else if is_gltf1 {
            node.id.clone()
        } else {
            String::new()
        },
        raw_name: if is_gltf1 { node.id.clone() } else { String::new() },
        style: XcafPrsStyle::default(),
        named_data: node.extras.clone(),
    };
    let mut split_keys = false;
    match doc.shape.shape_type() {
        ShapeType::Face | ShapeType::Edge | ShapeType::Vertex => {
            if let Some(face) = doc.product.clone() {
                // `cxx:2435-2442`: the node keeps the style of the shape it
                // replaces in the attribute map.
                if let Some(m) = face.style.material() {
                    att.style.set_material(Some(m.clone()));
                }
                if face.style.is_set_color_curv() {
                    att.style.set_color_curv(face.style.color_curv());
                }
                if face.style.is_set_color_surf() {
                    att.style.set_color_surf(face.style.color_surf_rgba());
                }
                if att.name.is_empty() && fallback {
                    // `cxx:2444-2447`
                    att.name = face.raw_name.clone();
                } else if !face.name.is_empty() && local_identity {
                    // `cxx:2449-2454`: keep the Product name (from the mesh) and
                    // the Instance name (from the node) in separate entries.
                    split_keys = true;
                }
            }
        }
        _ => {
            if att.name.is_empty() && fallback {
                if let Some(name) = common_late_name(doc) {
                    att.name = name;
                }
            }
            // `cxx:2489-2497` splits the keys for a non-face node too, but only
            // when the unlocated shape already carries a name in the attribute
            // map. A node compound is built fresh here, so that lookup is always
            // empty and the branch never fires.
        }
    }
    doc.split = !local_identity || split_keys;
    if !doc.split {
        // One map entry: the node attributes replace the mesh ones.
        doc.product = Some(att.clone());
    }
    doc.instance = Some(att);
}

/// The single `RWGltf_GltfLatePrimitiveArray::Name()` every face of `doc`
/// carries, else `None` (`RWGltf_GltfJsonParser.cxx:2461-2487`).
fn common_late_name(doc: &DocShape) -> Option<String> {
    fn collect(doc: &DocShape, out: &mut Vec<String>) {
        if doc.has_late {
            out.push(doc.late_name.clone());
        }
        for c in &doc.children {
            collect(c, out);
        }
    }
    let mut names = Vec::new();
    collect(doc, &mut names);
    let first = names.first()?;
    if first.is_empty() || names.iter().any(|n| n != first) {
        return None;
    }
    Some(first.clone())
}

// ---------------------------------------------------------------------------
// addShapeIntoDoc
// ---------------------------------------------------------------------------

/// `RWGltf_CafReader::addShapeIntoDoc` (`RWGltf_CafReader.cxx:470-688`): add
/// `doc` to the document under `parent` (a free root when `None`) and recurse
/// into its components / sub-shapes.
fn add_shape_into_doc(
    model: &mut BRepModel,
    tools: &mut DocTools,
    doc: &DocShape,
    parent: Option<usize>,
    parent_name: &str,
) -> Option<usize> {
    let shape_type = doc.shape.shape_type();
    // `cxx:523-550`: a compound becomes an assembly when a direct component is
    // not a face, or is a face that already carries a name (the mesh name).
    let to_make_assembly = shape_type == ShapeType::Compound
        && doc.children.iter().any(|c| {
            c.shape.shape_type() != ShapeType::Face
                || c.product.as_ref().map(|a| !a.name.is_empty()).unwrap_or(false)
        });
    let empty_compound = || TopoShape::new(ShapeType::Compound);

    let mut reused = false;
    let (new_label, product_label) = match parent {
        None => {
            // `cxx:553-557`: `AddShape`. A located shape comes back as a
            // reference to an unlocated Product (`XCAFDoc_ShapeTool.cxx:510-519`).
            let product_shape = if to_make_assembly { empty_compound() } else { doc.unlocated.clone() };
            model.add_free_shape(product_shape, doc.shape.clone(), to_make_assembly, doc.split)
        }
        Some(parent_label) if model.is_assembly(parent_label) => {
            // `cxx:558-578`: a component. An unlocated shape already in the
            // document is one more instance of its Product.
            let key = GeometryRegistry::shape_key(&doc.unlocated);
            match tools.component_map.get(&key).copied() {
                Some(product) => {
                    reused = true;
                    let label = model.add_component(parent_label, product)?;
                    if let Some(l) = model.label_mut(label) {
                        l.shape = Some(doc.shape.clone());
                    }
                    (label, product)
                }
                None => {
                    let product_shape =
                        if to_make_assembly { empty_compound() } else { doc.unlocated.clone() };
                    let (label, product) = model.add_component_shape(
                        parent_label,
                        product_shape,
                        doc.shape.clone(),
                        to_make_assembly,
                    )?;
                    tools.component_map.insert(key, product);
                    (label, product)
                }
            }
        }
        Some(parent_label) => {
            // `cxx:579-588`: a sub-shape of a simple-shape Product.
            let label = model.add_sub_shape_label(parent_label, doc.shape.clone())?;
            (label, label)
        }
    };

    let product_attribs = doc.product.clone().unwrap_or_default();
    let parent_present = parent.is_some();
    // `cxx:604-646`: the Instance attributes override the Product ones.
    let mut has_product_name = false;
    if new_label != product_label {
        let instance = if doc.split { doc.instance.as_ref() } else { None };
        match instance {
            Some(inst) => {
                if !inst.style.is_equal(&product_attribs.style) {
                    set_shape_style(model, tools, new_label, &inst.style);
                }
                if inst.named_data != product_attribs.named_data {
                    set_named_data(model, new_label, inst.named_data.as_ref());
                }
                set_shape_name(model, new_label, shape_type, &inst.name, parent_present, parent_name);
                if product_attribs.name.is_empty() && !inst.name.is_empty() {
                    // `cxx:626-632`: an unnamed Product copies its first
                    // instance's name.
                    has_product_name = true;
                    set_shape_name(model, product_label, shape_type, &inst.name, parent_present, parent_name);
                } else if inst.name.is_empty() && !product_attribs.name.is_empty() {
                    set_shape_name(
                        model,
                        new_label,
                        shape_type,
                        &product_attribs.name,
                        parent_present,
                        parent_name,
                    );
                }
            }
            None => set_shape_name(
                model,
                new_label,
                shape_type,
                &product_attribs.name,
                parent_present,
                parent_name,
            ),
        }
    }

    if reused {
        // `cxx:418-422`: the Product is already in the document.
        return Some(new_label);
    }

    // `cxx:652-658`: the Product attributes, shared across its instances.
    if !has_product_name {
        set_shape_name(model, product_label, shape_type, &product_attribs.name, parent_present, parent_name);
    }
    set_shape_style(model, tools, product_label, &product_attribs.style);
    set_named_data(model, product_label, product_attribs.named_data.as_ref());

    if model.is_assembly(product_label) {
        // `cxx:660-677`: the components of an assembly are shapes of their own.
        for child in &doc.children {
            add_shape_into_doc(model, tools, child, Some(product_label), "");
        }
    } else {
        // `cxx:678-687`: a simple shape keeps the children that carry their own
        // attributes as sub-shape labels.
        for child in &doc.children {
            add_sub_shape_into_doc(model, tools, child, product_label);
        }
    }
    Some(new_label)
}

/// `RWMesh_CafReader::addSubShapeIntoDoc` (`RWMesh_CafReader.cxx:456-511`).
fn add_sub_shape_into_doc(model: &mut BRepModel, tools: &mut DocTools, doc: &DocShape, parent: usize) {
    let shape_type = doc.shape.shape_type();
    let has_attribs = doc.product.is_some();
    // `cxx:471-477`: a shell or wire without attributes may still hold shaped
    // children, so it is walked anyway.
    if !has_attribs && shape_type != ShapeType::Shell && shape_type != ShapeType::Wire {
        return;
    }
    for child in &doc.children {
        add_sub_shape_into_doc(model, tools, child, parent);
    }
    let Some(attribs) = doc.product.clone() else { return };
    let Some(label) = model.add_sub_shape_label(parent, doc.shape.clone()) else { return };
    set_shape_name(model, label, shape_type, &attribs.name, false, "");
    set_shape_style(model, tools, label, &attribs.style);
    set_named_data(model, label, attribs.named_data.as_ref());
}

// ---------------------------------------------------------------------------
// Label attribute helpers
// ---------------------------------------------------------------------------

/// `RWMesh_CafReader::setShapeName` (`RWMesh_CafReader.cxx:199-217`). The glTF
/// reader passes an empty `theParentName` (`RWGltf_CafReader.cxx:454`, `:674`),
/// so that arm never names a label here.
fn set_shape_name(
    model: &mut BRepModel,
    label: usize,
    shape_type: ShapeType,
    name: &str,
    parent_present: bool,
    parent_name: &str,
) {
    if !name.is_empty() {
        model.set_label_name(label, name);
    } else if parent_present {
        model.set_label_name(label, shape_type_to_string(shape_type));
    } else if !parent_name.is_empty() {
        model.set_label_name(label, parent_name);
    }
}

/// `RWMesh_CafReader::setShapeStyle` (`RWMesh_CafReader.cxx:221-246`): the
/// `XCAFDoc_ColorTool` roles and the `VisMaterialRefGUID` link.
fn set_shape_style(model: &mut BRepModel, tools: &mut DocTools, label: usize, style: &XcafPrsStyle) {
    if style.is_set_color_surf() {
        if let Some(l) = model.label_mut(label) {
            l.color = Some(style.color_surf());
        }
    }
    if style.is_set_color_curv() {
        if let Some(l) = model.label_mut(label) {
            l.color_curv = Some(style.color_curv());
        }
    }
    if let Some(mat) = style.material() {
        let index = match tools.materials.get(&mat.id()).copied() {
            Some(i) => i,
            None => {
                let i = model.add_material(mat.clone(), mat.raw_name().unwrap_or(""));
                tools.materials.insert(mat.id(), i);
                i
            }
        };
        if let Some(l) = model.label_mut(label) {
            l.material = Some(index);
        }
    }
}

/// `RWMesh_CafReader::setShapeNamedData` (`RWMesh_CafReader.cxx:248-277`).
fn set_named_data(model: &mut BRepModel, label: usize, named_data: Option<&TDataStdNamedData>) {
    let Some(nd) = named_data else { return };
    if let Some(l) = model.label_mut(label) {
        l.named_data = Some(nd.clone());
    }
}

/// `RWMesh_CafReader::shapeTypeToString` (`RWMesh_CafReader.hxx:317-322`):
/// `TopAbs::ShapeTypeToString` with the first letter capitalized.
fn shape_type_to_string(t: ShapeType) -> &'static str {
    match t {
        ShapeType::Compound => "Compound",
        ShapeType::CompSolid => "Compsolid",
        ShapeType::Solid => "Solid",
        ShapeType::Shell => "Shell",
        ShapeType::Face => "Face",
        ShapeType::Wire => "Wire",
        ShapeType::Edge => "Edge",
        ShapeType::Vertex => "Vertex",
        ShapeType::Shape => "Shape",
    }
}
