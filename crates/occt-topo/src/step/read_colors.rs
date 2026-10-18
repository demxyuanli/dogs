//! STEP colours onto a `BRepModel`. Source: `STEPCAFControl_Reader`
//! (`ReadColors` `cxx:1485-1547`, `SetStyle` `cxx:1223-1447`, `IsOverriden`
//! `cxx:1456-1480`) plus `XCAFPrs::CollectStyleSettings` (`XCAFPrs.cxx:118-330`).
//!
//! OCCT links a styled item to a *TDF label*; this crate has no document tree,
//! so the `BRepModel`'s shapes play the label roles:
//!
//! * one document shape per transferred representation, i.e. one model shape
//!   per `rep_index` entry (the reader-side result of `STEPCAFControl_Reader::
//!   AddShape(reader.OneShape(), ...)`, `cxx:750`), and
//! * `XCAFDoc_ShapeTool::SearchUsingMap(aS, aL, false, true)` (`cxx:134-195`):
//!   the model shape that *is* `aS` (`myShapeLabels`), else
//!   `FindMainShapeUsingMap` (`cxx:1193-1205`, the `mySubShapes` registry that
//!   `XCADoc_ShapeTool::addShape` fills with every strict sub-shape of a
//!   non-assembly shape, `cxx:565-583`) through `AddSubShape`
//!   (`cxx:1139-1187`), which needs `IsSimpleShape(shapeL) && IsTopLevel(shapeL)`.
//!
//! The two colour sinks are [`crate::model::ModelShape`]'s:
//!
//! * `color` — the `XCAFDoc_ColorSurf` of the shape's own label, i.e. the
//!   `XCAFPrs_DocumentNode::Style` the OBJ writer receives as `aDocNode.Style`
//!   (`XCAFPrs_DocumentExplorer::mergedStyle`, `cxx:28-81`); and
//! * `style_settings` — `XCAFPrs_Style` per stylised label,
//!   `XCAFPrs::CollectStyleSettings`' shape-keyed `theSettings` map
//!   (`XCAFPrs.cxx:307-330`): the labels `SetColor` / `SetVisibility` wrote plus
//!   the shape's own label, empty styles and empty compounds dropped.
//!
//! **UNPORTED**, with `SetStyle`'s own numbering:
//!
//! * the SHUO arms — the `anIsComponent` / `findStyledSR` `while` loop
//!   (`cxx:1301-1364`), `SetAssemblyComponentStyle` (`cxx:1250-1258`) and
//!   `CollectStyleSettings`' SHUO block (`XCAFPrs.cxx:221-303`): this crate has
//!   no NAUO/SHUO graph (`StepAssembly` carries the product tree only);
//! * the `XCAFDoc_VisMaterialTool` arms — every `Material()` stays null, both
//!   here and via `ReadMaterials` (`cxx:5596`), as documented by
//!   [`crate::xcaf_prs`];
//! * the `XCAFDoc_LayerTool` layer arms (`XCAFPrs.cxx:128-139`, `:176-216`):
//!   `ModelShape.layer` is a name, not a layer tree with its own colours;
//! * `propagateColorToParts` (`cxx:1037-1077`) after `ReadColors`: it only
//!   propagates a colour carried by an *assembly* label, and the port's model
//!   shapes are the composed roots, never assembly labels.

use super::construct_styles::{ConstructStyles, StyleColors};
use super::prelude::*;
use super::*;
use crate::iterator::cumulated_children;
use crate::model::ModelShape;
use crate::topo_tools_full::is_same;
use crate::xcaf_prs::XcafPrsStyle;
use occt_core::quantity::ColorRGBA;

/// `STEPCAFControl_Reader::ReadColors` (`cxx:1485-1547`).
///
/// `rep_index` maps a shape-representation entity (`ADVANCED_BREP_SHAPE_
/// REPRESENTATION` / `SHAPE_REPRESENTATION`, or the bare `MANIFOLD_SOLID_BREP`
/// root of the representation-less fallback) to the model shape the reader
/// built from it; it is the port's `AddShape` label registry.
pub(super) fn read_colors(
    model: &mut BRepModel,
    resolver: &Resolver,
    rep_index: &HashMap<usize, usize>,
) {
    let records = resolver.records;
    let mut styles = ConstructStyles::new(records);
    if !styles.load_styles() {
        // `if (!Styles.LoadStyles()) return false;` (`cxx:1490-1495`).
        return;
    }
    // `Styles.LoadInvisStyles(aHSeqOfInvisStyle)` (`cxx:1499-1500`).
    let invis = styles.load_invis_styles();
    // `const Interface_Graph& aGraph = Styles.Graph()` (`cxx:1514`).
    let referrers = reference_index(records);

    // "apply root styles earlier, as they can be overridden" (`cxx:1517-1528`).
    for i in 1..=styles.nb_root_styles() {
        let style = styles.root_style(i);
        if !is_overriden(&referrers, records, style, true) {
            set_style(model, resolver, rep_index, &styles, &invis, style);
        }
    }
    // "apply leaf styles, they can override root styles" (`cxx:1530-1541`).
    for i in 1..=styles.nb_styles() {
        let style = styles.style(i);
        if !is_overriden(&referrers, records, style, false) {
            set_style(model, resolver, rep_index, &styles, &invis, style);
        }
    }

    // `CollectStyleSettings` for every node (`XCAFPrs.cxx:118-330`), run after
    // the whole `SetColor` / `SetVisibility` pass left the labels final.
    for ms in &mut model.shapes {
        collect_style_settings(ms);
    }
}

/// `SetStyle` (`STEPCAFControl_Reader.cxx:1223-1447`).
fn set_style(
    model: &mut BRepModel,
    resolver: &Resolver,
    rep_index: &HashMap<usize, usize>,
    styles: &ConstructStyles<'_>,
    invis: &[usize],
    style_id: usize,
) {
    let Some(rec) = styles.record_opt(style_id) else { return };
    let is_over_riding = is_over_riding_styled_item_type(&rec.type_name);

    // `if (occ::handle<StepVisual_OverRidingStyledItem> anOverridingStyle = ...)`
    // (`cxx:1239-1253`): the over-ridden style is written first, so the
    // over-riding item can overwrite its colour slots. The
    // `ContextDependentOverRidingStyledItem` arm (`cxx:1250-1258`) is UNPORTED
    // (module header).
    if is_over_riding {
        if let Some(over_ridden) = rec.args.get(3).and_then(|a| parse_ref(a)) {
            set_style(model, resolver, rep_index, styles, &invis, over_ridden);
        }
    }

    // `bool anIsVisible = true; for (int si ...)` (`cxx:1262-1273`).
    let is_visible = !invis.contains(&style_id);

    // `if (!theStyles.GetColors(...) && anIsVisible) return;` (`cxx:1276-1283`).
    let colors = styles.colors(style_id);
    if !colors.is_defined() && is_visible {
        return;
    }

    // `int anIndex = aTP->MapIndex(anItems.Value(itemIt).Value()); ... aS =
    // TransferBRep::ShapeResult(aBinder);` (`cxx:1293-1300`): the item's
    // *transferred* shape. An item the transfer never bound keeps an index 0,
    // i.e. a null `aS`, and is skipped (`cxx:1368-1372`).
    let Some(item) = rec.args.get(2).and_then(|a| parse_ref(a)) else { return };
    let Some(s) = transferred_shape(model, resolver, rep_index, item) else {
        return;
    };

    // `TDF_Label aL; bool isFound = theSTool->SearchUsingMap(aS, aL, false, true);`
    // (`cxx:1377-1378`).
    let found = search_using_map(model, &s);

    // `if (!aSurfCol.IsNull() || !aBoundCol.IsNull() || !aCurveCol.IsNull() ||
    // aRenderProps.IsDefined()) { ... }` (`cxx:1379-1437`).
    if colors.is_defined() {
        if let Some(label) = &found {
            // `if (isFound) { theCTool->SetColor(aL, ...); }` (`cxx:1403-1415`).
            let style = style_from_colors(resolver.records, &colors, is_visible);
            add_setting(model, label, style);
        } else {
            // `else { for (TopoDS_Iterator it(aS); it.More(); it.Next()) {
            // if (theSTool->SearchUsingMap(it.Value(), aL1, false, true)) {
            // theCTool->SetColor(aL1, ...); } } }` (`cxx:1417-1437`). The child
            // label keeps the default visibility: only the colour is written.
            for child in cumulated_children(&s) {
                let Some(child_label) = search_using_map(model, &child) else { continue };
                let style = style_from_colors(resolver.records, &colors, true);
                add_setting(model, &child_label, style);
            }
        }
    }

    // `if (!anIsVisible) { if (isFound) theCTool->SetVisibility(aL, false); }`
    // (`cxx:1439-1443`) — the visibility attribute of the found label, applied
    // whether or not the style carried a colour.
    if !is_visible {
        if let Some(label) = &found {
            let mut style = XcafPrsStyle::default();
            style.set_visibility(false);
            add_setting(model, label, style);
        }
    }
}

/// The label `SearchUsingMap(aS, aL, false, true)` finds: a model shape's own
/// label, or the `AddSubShape` sub-shape label under it.
enum Found {
    Label { index: usize, key_shape: TopoShape },
}

/// `XCAFDoc_ShapeTool::SearchUsingMap` (`cxx:134-195`) restricted to the
/// branches this crate's label set can reach. The `mySimpleShapes` /
/// `GetUsers` arms (`cxx:145-176`) belong to located shapes and assemblies,
/// which have no counterpart here: a model shape is its own top-level label,
/// never a component reference (`cxx:849-856`), and no model shape is an
/// `XCAFDoc::AssemblyGUID` holder (`cxx:539-546`), because the assembly path
/// (`cxx:929-945`) builds a `StepAssembly`, not a `BRepModel`.
fn search_using_map(model: &BRepModel, s: &TopoShape) -> Option<Found> {
    // `if (myShapeLabels.IsBound(S)) { L = myShapeLabels.Find(S); return true; }`
    // (`cxx:139-143`). `NCollection_DataMap::Bind` keeps the last binding, so
    // the last model shape holding that exact shape wins.
    let mut exact: Option<usize> = None;
    for (i, ms) in model.shapes.iter().enumerate() {
        if is_same(&ms.shape, s) {
            exact = Some(i);
        }
    }
    if let Some(i) = exact {
        return Some(Found::Label { index: i, key_shape: model.shapes[i].shape.clone() });
    }

    // `TDF_Label mainL = FindMainShapeUsingMap(S); if (mainL.IsNull()) return
    // false; L = AddSubShape(mainL, S);` (`cxx:188-193`).
    let mut main: Option<usize> = None;
    for (i, ms) in model.shapes.iter().enumerate() {
        if contains_subshape(&ms.shape, s) {
            main = Some(i);
        }
    }
    let main = main?;
    // `AddSubShape`: `if (!IsSimpleShape(shapeL) || !IsTopLevel(shapeL)) return
    // false;` (`cxx:1143-1147`) — see the function doc for why both hold.
    // `if (!IsSubShape(shapeL, sub))` (`cxx:1152`): `FindMainShapeUsingMap`
    // already answered that.
    Some(Found::Label { index: main, key_shape: s.clone() })
}

/// `XCAFDoc_ShapeMapTool::SetShape` + `IsSubShape`
/// (`XCAFDoc_ShapeMapTool.cxx:89-108`): the registry holds the shape's
/// *children* and everything below them, so a shape is never its own sub-shape.
fn contains_subshape(main: &TopoShape, sub: &TopoShape) -> bool {
    cumulated_children(main).iter().any(|c| subtree_contains(c, sub))
}

/// `AddSubShape`'s recursive `myMap.Add` walk (`XCAFDoc_ShapeMapTool.cxx:89-97`).
fn subtree_contains(root: &TopoShape, sub: &TopoShape) -> bool {
    is_same(root, sub) || cumulated_children(root).iter().any(|c| subtree_contains(c, sub))
}

/// `XCAFPrs_Style` for one label: `fillStyleColors` (`XCAFPrs.cxx:34-52`) plus
/// `SetVisibility` (`cxx:172`, `XCAFPrs.cxx`). `SetStyle` only ever writes
/// `XCAFDoc_ColorSurf` / `XCAFDoc_ColorCurv` (`STEPCAFControl_Reader.cxx:
/// 1405-1414`), and `GetColor(L, type, color)` reads back the attribute of
/// exactly that type (`XCAFDoc_ColorTool.cxx:278-302`), so `Gen` — which
/// `fillStyleColors` would copy into *both* slots — is never set here.
fn style_from_colors(
    records: &HashMap<usize, Record>,
    colors: &StyleColors,
    is_visible: bool,
) -> XcafPrsStyle {
    let mut style = XcafPrsStyle::default();
    style.set_visibility(is_visible);
    // `Quantity_ColorRGBA aFullSCol; if (!aSurfCol.IsNull()) { DecodeColor(...);
    // aFullSCol = Quantity_ColorRGBA(aSCol); } ... if (aRenderProps.IsDefined())
    // aFullSCol = aRenderProps.GetRGBAColor();` (`cxx:1384-1399`): the rendering
    // properties' colour replaces the decoded one when both are present.
    let surface = if colors.rendering_props.is_defined() {
        Some(colors.rendering_props.rgba_color())
    } else {
        colors
            .surface_color
            .and_then(|c| decode_color(records, c))
            .map(|c| ColorRGBA { r: c.r, g: c.g, b: c.b, a: 1.0 })
    };
    if let Some(surface) = surface {
        style.set_color_surf(surface);
    }
    // `if (!aBoundCol.IsNull()) theCTool->SetColor(aL, aBCol, XCAFDoc_ColorCurv);`
    // then `if (!aCurveCol.IsNull()) SetColor(aL, aCCol, XCAFDoc_ColorCurv);`
    // (`cxx:1406-1414`) — both write the same attribute, so the curve colour
    // wins.
    let curve = colors
        .curve_color
        .and_then(|c| decode_color(records, c))
        .or_else(|| colors.boundary_color.and_then(|c| decode_color(records, c)));
    if let Some(curve) = curve {
        style.set_color_curv(curve);
    }
    style
}

/// `theCTool->SetColor(aL, ...)` / `SetVisibility(aL, false)`: write one label's
/// attribute. `SetColor` replaces the attribute of its own type, so a label
/// accumulates colour slots across styled items (`XCAFDoc_ColorTool.cxx:
/// 304-360`); `SetVisibility(false)` only ever sets the `InvisibleGUID`
/// attribute, so visibility merges with `&&`
/// (`XCAFDoc_ColorTool.cxx:442-462`).
fn add_setting(model: &mut BRepModel, label: &Found, style: XcafPrsStyle) {
    let Found::Label { index, key_shape } = label;
    let ms = &mut model.shapes[*index];
    if let Some(slot) = ms.style_settings.iter_mut().find(|(s, _)| is_same(s, key_shape)) {
        if style.is_set_color_surf() {
            slot.1.set_color_surf(style.color_surf_rgba());
        }
        if style.is_set_color_curv() {
            slot.1.set_color_curv(style.color_curv());
        }
        slot.1.set_visibility(slot.1.is_visible() && style.is_visible());
        return;
    }
    ms.style_settings.push((key_shape.clone(), style));
}

/// `XCAFPrs::CollectStyleSettings` (`XCAFPrs.cxx:118-330`) for one node label:
/// the shape-keyed settings map, restricted to the port's label set.
fn collect_style_settings(ms: &mut ModelShape) {
    // `if (aStyle.IsEmpty()) continue;` (`cxx:307-310`) — an invisible style is
    // not empty, a visible style with no colour is.
    ms.style_settings.retain(|(_, style)| !style.is_empty());
    // `if (aSubshape.ShapeType() == TopAbs_COMPOUND) { if (aSubshape.NbChildren()
    // == 0) continue; }` (`cxx:312-318`).
    ms.style_settings
        .retain(|(key, _)| !key.is_compound() || !cumulated_children(key).is_empty());
    // `XCAFPrs_DocumentExplorer::mergedStyle` (`cxx:44-57`): the label's
    // `XCAFDoc_ColorSurf` is the node style the OBJ writer reads.
    if let Some((_, style)) = ms
        .style_settings
        .iter()
        .find(|(key, _)| is_same(key, &ms.shape))
    {
        if style.is_set_color_surf() {
            ms.color = Some(style.color_surf());
        }
    }
}

/// The shape the transfer bound to a styled item's target
/// (`aTP->MapIndex` / `TransferBRep::ShapeResult`): a document shape for a
/// shape-representation target, else the resolver's transferred shape.
fn transferred_shape(
    model: &BRepModel,
    resolver: &Resolver,
    rep_index: &HashMap<usize, usize>,
    item: usize,
) -> Option<TopoShape> {
    if let Some(&i) = rep_index.get(&item) {
        return Some(model.shapes[i].shape.clone());
    }
    resolver.shape_cache.borrow().get(&item).cloned()
}

/// `IsOverriden` (`STEPCAFControl_Reader.cxx:1456-1480`).
fn is_overriden(
    referrers: &HashMap<usize, Vec<usize>>,
    records: &HashMap<usize, Record>,
    style: usize,
    is_root: bool,
) -> bool {
    let Some(sharers) = referrers.get(&style) else { return false };
    for &sharer in sharers {
        let Some(rec) = records.get(&sharer) else { continue };
        // `occ::down_cast<StepVisual_OverRidingStyledItem>(aSubs.Value())`
        // (`cxx:1464-1466`): the subtype family, whose only reference to a
        // styled item is `over_ridden_style` (attribute 4,
        // `RWStepVisual_RWOverRidingStyledItem.cxx:78-87`), i.e. exactly the
        // sharings that make `theStyle` an over-ridden style.
        if !is_over_riding_styled_item_type(&rec.type_name) {
            continue;
        }
        if !is_root {
            return true;
        }
        // "for root style returns true only if it is overridden by other root
        // style" (`cxx:1470-1477`): the over-riding item's AP242 target must be
        // a `StepShape_ShapeRepresentation`.
        if rec
            .args
            .get(2)
            .and_then(|a| parse_ref(a))
            .and_then(|t| records.get(&t))
            .map(|t| is_shape_representation_type(&t.type_name))
            .unwrap_or(false)
        {
            return true;
        }
    }
    false
}

/// `StepVisual_OverRidingStyledItem` and its subtype
/// (`StepVisual_ContextDependentOverRidingStyledItem`).
fn is_over_riding_styled_item_type(type_name: &str) -> bool {
    matches!(
        type_name,
        "OVER_RIDING_STYLED_ITEM" | "CONTEXT_DEPENDENT_OVER_RIDING_STYLED_ITEM"
    )
}

/// `Interface_Graph`'s reference relation restricted to what [`is_overriden`]
/// reads: referenced entity -> the records referencing it. Built once because
/// the crate's model is a plain record table instead of an `Interface_Graph`.
fn reference_index(records: &HashMap<usize, Record>) -> HashMap<usize, Vec<usize>> {
    let mut index: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut ids: Vec<usize> = records.keys().copied().collect();
    ids.sort_unstable();
    for id in ids {
        let Some(rec) = records.get(&id) else { continue };
        let mut targets: Vec<usize> = Vec::new();
        for arg in &rec.args {
            refs_in(arg, &mut targets);
        }
        for (_, member_args) in &rec.members {
            for arg in member_args {
                refs_in(arg, &mut targets);
            }
        }
        for target in targets {
            let sharers = index.entry(target).or_default();
            if !sharers.contains(&id) {
                sharers.push(id);
            }
        }
    }
    index
}

/// Every `#N` one argument mentions, nested `(a,b)` lists included.
fn refs_in(arg: &str, out: &mut Vec<usize>) {
    let trimmed = arg.trim();
    if trimmed.starts_with('(') && trimmed.ends_with(')') {
        for element in split_top(&trimmed[1..trimmed.len() - 1]) {
            refs_in(&element, out);
        }
    } else if let Some(r) = parse_ref(trimmed) {
        out.push(r);
    }
}
