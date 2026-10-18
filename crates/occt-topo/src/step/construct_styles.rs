//! `STEPConstruct_Styles` — reading and writing shape styles (color) to and
//! from a STEP file. Source: `STEPConstruct_Styles.cxx` / `.hxx`.
//!
//! Reader side, which is what the OBJ export needs:
//!
//! * [`ConstructStyles::load_styles`] — `LoadStyles` (`cxx:498-554`): collects
//!   the `STYLED_ITEM`s reachable from
//!   `MECHANICAL_DESIGN_GEOMETRIC_PRESENTATION_REPRESENTATION` / `DRAUGHTING_MODEL`
//!   (plus every free-standing styled item), splitting them into the
//!   representation-level *root* styles (`myRootStyles`) and the remaining
//!   *leaf* styles (`myStyles`).
//! * [`ConstructStyles::load_invis_styles`] — `LoadInvisStyles`
//!   (`cxx:559-592`).
//! * [`ConstructStyles::get_colors`] — `GetColors` (`cxx:735-770`) with its
//!   `ProcessAs*` helpers (`cxx:76-220`).
//! * [`ConstructStyles::decode_color`] — `DecodeColor` (`cxx:918-991`).
//!
//! **UNPORTED**: the writing half — `MakeColorPSA` (`cxx:597-713`),
//! `GetColorPSA` (`cxx:715-733`), `EncodeColor` (`cxx:774-916`),
//! `CreateMDGPR` (`cxx:341-379`), `CreateNAUOSRD` (`cxx:380-444`),
//! `FindContext` (`cxx:445-495`) and `AddStyle` (`cxx:283-340`). The crate's
//! STEP writer emits the plain `COLOUR_RGB` `STYLED_ITEM` form directly
//! (`step/transfer.rs::write_step_model_with_colors`), so none of them has a
//! caller; they are documented here rather than faked.
//!
//! Entity order: OCCT walks `Interface_InterfaceModel::Value(1..NbEntities())`,
//! i.e. the DATA-section order. This crate stores records in a `HashMap` keyed
//! by entity number, so the walk uses ascending entity numbers — the order
//! `read_step_impl` already relies on for its own representation scan. A file
//! whose entities are not numbered in file order would reorder `myStyles`.

use std::collections::HashMap;

use occt_core::quantity::Color;

use super::*;
use crate::step::rendering_properties::RenderingProperties;

/// `GetColors`' out-parameters (`cxx:735-741`). Each color is an entity id
/// (`occ::handle<StepVisual_Colour>` in OCCT); `None` is a null handle.
pub(in crate::step) struct StyleColors {
    /// `theSurfaceColour`.
    pub surface_color: Option<usize>,
    /// `theBoundaryColour`.
    pub boundary_color: Option<usize>,
    /// `theCurveColour`.
    pub curve_color: Option<usize>,
    /// `theRenderingProps`.
    pub rendering_props: RenderingProperties,
    /// `theIsComponent`.
    pub is_component: bool,
}

impl StyleColors {
    /// `GetColors`' return value (`cxx:768-769`).
    pub(in crate::step) fn is_defined(&self) -> bool {
        self.surface_color.is_some()
            || self.boundary_color.is_some()
            || self.curve_color.is_some()
            || self.rendering_props.is_defined()
    }
}

/// `STEPConstruct_Styles` (`STEPConstruct_Styles.hxx:54-192`), restricted to
/// the read path (`myStyles` / `myRootStyles`).
pub(in crate::step) struct ConstructStyles<'a> {
    records: &'a HashMap<usize, Record>,
    /// `myStyles` (`hxx:188`) — an `NCollection_IndexedMap<Standard_Transient>`,
    /// i.e. an ordered set of handles.
    styles: Vec<usize>,
    /// `myRootStyles` (`hxx:189`).
    root_styles: Vec<usize>,
}

impl<'a> ConstructStyles<'a> {
    /// `STEPConstruct_Styles::STEPConstruct_Styles`
    /// (`STEPConstruct_Styles.cxx:33-38`): an empty tool. `Init(WS)`
    /// (`cxx:236-244`) only stores the work session, so the port's constructor
    /// takes the record table directly.
    pub(in crate::step) fn new(records: &'a HashMap<usize, Record>) -> Self {
        Self { records, styles: Vec::new(), root_styles: Vec::new() }
    }

    /// The STEP model, as `StepRepr_Representation` lookups and `DecodeColor`
    /// need it.
    pub(in crate::step) fn records(&self) -> &'a HashMap<usize, Record> {
        self.records
    }

    pub(in crate::step) fn record_opt(&self, id: usize) -> Option<&'a Record> {
        self.records.get(&id)
    }

    /// `STEPConstruct_Styles::NbStyles` (`cxx:246-249`).
    pub(in crate::step) fn nb_styles(&self) -> usize {
        self.styles.len()
    }

    /// `STEPConstruct_Styles::Style` (`cxx:253-258`) — `i` is 1-based.
    pub(in crate::step) fn style(&self, i: usize) -> usize {
        self.styles[i - 1]
    }

    /// `STEPConstruct_Styles::NbRootStyles` (`cxx:260-263`).
    pub(in crate::step) fn nb_root_styles(&self) -> usize {
        self.root_styles.len()
    }

    /// `STEPConstruct_Styles::RootStyle` (`cxx:267-272`) — `i` is 1-based.
    pub(in crate::step) fn root_style(&self, i: usize) -> usize {
        self.root_styles[i - 1]
    }

    /// `STEPConstruct_Styles::LoadStyles` (`cxx:498-554`).
    ///
    /// `myPSA` (`cxx:500`) is only used by the writing half, so it is not
    /// carried. `tMDGPR` / `tDM` are exact `DynamicType()` matches
    /// (`cxx:514`), i.e. record type names; `tSI` is an `IsKind` match
    /// (`cxx:540`), i.e. the `StyledItem` family.
    pub(in crate::step) fn load_styles(&mut self) -> bool {
        self.styles.clear();
        self.root_styles.clear();

        let mut ids: Vec<usize> = self.records.keys().copied().collect();
        ids.sort_unstable();
        for id in ids {
            let Some(rec) = self.records.get(&id) else { continue };
            if rec.type_name == "MECHANICAL_DESIGN_GEOMETRIC_PRESENTATION_REPRESENTATION"
                || rec.type_name == "DRAUGHTING_MODEL"
            {
                // `container = down_cast<StepRepr_Representation>(enti)`
                // (`cxx:517-518`); `container->ItemsValue(j)` walks the `items`
                // attribute (`cxx:520-529`).
                for item in rec.args.get(1).map(|a| parse_ref_list(a)).unwrap_or_default() {
                    let Some(irec) = self.records.get(&item) else { continue };
                    if !is_styled_item_type(&irec.type_name) {
                        continue;
                    }
                    self.add_style(item);
                }
            } else if is_styled_item_type(&rec.type_name) {
                // `else if (enti->IsKind(tSI))` (`cxx:540-552`).
                self.add_style(id);
            }
        }
        !self.styles.is_empty()
    }

    /// The `cxx:526-552` classification and the `Contains` guards: a styled
    /// item whose AP242 target is a shape representation is a root style,
    /// everything else a leaf style.
    fn add_style(&mut self, id: usize) {
        let is_root = self
            .records
            .get(&id)
            .and_then(|r| r.args.get(2))
            .and_then(|a| parse_ref(a))
            .and_then(|t| self.records.get(&t))
            .map(|t| is_shape_representation_type(&t.type_name))
            .unwrap_or(false);
        if is_root {
            if !self.root_styles.contains(&id) {
                self.root_styles.push(id);
            }
        } else if !self.styles.contains(&id) {
            self.styles.push(id);
        }
    }

    /// `STEPConstruct_Styles::LoadInvisStyles` (`cxx:559-592`): every
    /// `STYLED_ITEM` referenced by an `INVISIBILITY` entity.
    pub(in crate::step) fn load_invis_styles(&self) -> Vec<usize> {
        let mut out: Vec<usize> = Vec::new();
        let mut ids: Vec<usize> = self.records.keys().copied().collect();
        ids.sort_unstable();
        for id in ids {
            // `if (enti->DynamicType() != tInvisibility) continue;`
            // (`cxx:569-572`).
            let Some(rec) = self.records.get(&id) else { continue };
            if rec.type_name != "INVISIBILITY" {
                continue;
            }
            // `container->InvisibleItemsValue(j)` -> `StepVisual_InvisibleItem`
            // select, case 1 when the reference is a `StepVisual_StyledItem`
            // (`StepVisual_InvisibleItem.cxx`, `CaseNum` -> 1 -> `StyledItem()`).
            let items = rec.args.first().map(|a| parse_ref_list(a)).unwrap_or_default();
            for item in items {
                let Some(irec) = self.records.get(&item) else { continue };
                if !is_styled_item_type(&irec.type_name) {
                    continue;
                }
                if !out.contains(&item) {
                    out.push(item);
                }
            }
        }
        out
    }

    /// `STEPConstruct_Styles::GetColors` (`cxx:735-770`).
    pub(in crate::step) fn colors(&self, style_id: usize) -> StyleColors {
        let mut res = StyleColors {
            surface_color: None,
            boundary_color: None,
            curve_color: None,
            rendering_props: RenderingProperties::default(),
            is_component: false,
        };
        let Some(rec) = self.records.get(&style_id) else { return res };

        // `for (int aPSAIndex = 1; aPSAIndex <= theStyle->NbStyles(); ++aPSAIndex)`
        // (`cxx:749-755`): the `styles` set of the styled item (attribute 2),
        // each entry a `StepVisual_PresentationStyleAssignment`.
        let psas = rec.args.get(1).map(|a| parse_ref_list(a)).unwrap_or_default();
        for psa in psas {
            let Some(prec) = self.records.get(&psa) else { continue };
            if prec.type_name != "PRESENTATION_STYLE_ASSIGNMENT" {
                continue;
            }
            let Some(pss_list) = prec.args.first() else { continue };
            let pss_ids = parse_ref_list(pss_list);
            if pss_ids.is_empty() {
                // `if (aPSA.IsNull() || aPSA->Styles().IsNull()) continue;`
                // (`cxx:751-754`).
                continue;
            }
            res.is_component = true;

            for pss in pss_ids {
                // `ProcessAsSurfaceStyleUsage(...) || ProcessAsCurveStyle(...)`
                // (`cxx:764-765`) — `||` short-circuits.
                if !process_as_surface_style_usage(
                    self.records,
                    pss,
                    &mut res.surface_color,
                    &mut res.boundary_color,
                    &mut res.rendering_props,
                ) {
                    process_as_curve_style(self.records, pss, &mut res.curve_color);
                }
            }
        }
        res
    }

    /// The `GetColors` return value, as `SetStyle` reads it (`cxx:1279-1282`).
    pub(in crate::step) fn has_colors(&self, style_id: usize) -> Option<StyleColors> {
        let colors = self.colors(style_id);
        if colors.is_defined() { Some(colors) } else { None }
    }
}

/// `STEPConstruct_Styles::DecodeColor` (`cxx:918-991`). `Color` carries
/// sRGB in this port, so `SetValues(..., Quantity_TOC_sRGB)` is a plain
/// component store.
pub(in crate::step) fn decode_color(
    records: &HashMap<usize, Record>,
    color: usize,
) -> Option<Color> {
    let rec = records.get(&color)?;
    // `Colour->IsKind(STANDARD_TYPE(StepVisual_ColourRgb))` (`cxx:921-942`).
    if rec.type_name == "COLOUR_RGB" {
        let r = parse_f64(rec.args.get(1).map(String::as_str).unwrap_or("$")).ok()?;
        let g = parse_f64(rec.args.get(2).map(String::as_str).unwrap_or("$")).ok()?;
        let b = parse_f64(rec.args.get(3).map(String::as_str).unwrap_or("$")).ok()?;
        // `if (rgb->Red() > 1. || rgb->Green() > 1. || rgb->Blue() > 1.)`
        // normalises by the largest component (`cxx:922-939`).
        let (r, g, b) = if r > 1.0 || g > 1.0 || b > 1.0 {
            let norm = r.max(g).max(b);
            (r / norm, g / norm, b / norm)
        } else {
            (r, g, b)
        };
        return Some(Color::new(r as f32, g as f32, b as f32));
    }
    // `Colour->IsKind(STANDARD_TYPE(StepVisual_PreDefinedColour))`
    // (`cxx:943-988`): a physical file names the concrete type
    // `PRE_DEFINED_COLOUR` or its subtype `DRAUGHTING_PRE_DEFINED_COLOUR`
    // (`StepVisual_DraughtingPreDefinedColour`).
    //
    // Both readers take the item's `name` as a single *string* parameter
    // (`RWStepVisual_RWPreDefinedColour.cxx:40-45`,
    // `RWStepVisual_RWDraughtingPreDefinedColour.cxx:40-45`:
    // `data->ReadString(num, 1, "name", ach, aName)`) and store it in the
    // entity's own `StepVisual_PreDefinedItem`
    // (`StepVisual_PreDefinedColour.cxx:23-26` constructs it inline in the
    // constructor). So `GetPreDefinedItem()->Name()` (`cxx:947-948`) is that
    // literal string: the file writes `DRAUGHTING_PRE_DEFINED_COLOUR('white')`,
    // one string argument, never an entity reference.
    if rec.type_name == "DRAUGHTING_PRE_DEFINED_COLOUR" || rec.type_name == "PRE_DEFINED_COLOUR" {
        let name = parse_str(rec.args.first()?);
        // `name.IsEqual(...)` against the eight predefined colors
        // (`cxx:949-987`); an unrecognised name returns false (`cxx:982-987`).
        return predefined_color(&name);
    }
    None
}

/// `Quantity_NOC_*` values for the eight names `DecodeColor` recognises
/// (`cxx:949-980`), as sRGB triples (`Quantity_Color`'s predefined table).
fn predefined_color(name: &str) -> Option<Color> {
    match name {
        "red" => Some(Color::RED),
        "green" => Some(Color::GREEN),
        "blue" => Some(Color::BLUE),
        "yellow" => Some(Color::YELLOW),
        "magenta" => Some(Color::MAGENTA),
        "cyan" => Some(Color::CYAN),
        "black" => Some(Color::BLACK),
        "white" => Some(Color::WHITE),
        _ => None,
    }
}

/// `StepVisual_StyledItem` and its subtypes (`MyISKIND(tSI)`, `cxx:540`).
/// `CONTEXT_DEPENDENT_OVER_RIDING_STYLED_ITEM` derives from
/// `OVER_RIDING_STYLED_ITEM` (`RWStepVisual_RWContextDependentOverRidingStyledItem.hxx`).
pub(in crate::step) fn is_styled_item_type(type_name: &str) -> bool {
    matches!(
        type_name,
        "STYLED_ITEM"
            | "OVER_RIDING_STYLED_ITEM"
            | "CONTEXT_DEPENDENT_OVER_RIDING_STYLED_ITEM"
    )
}

/// `anItem->IsKind(STANDARD_TYPE(StepShape_ShapeRepresentation))`
/// (`STEPConstruct_Styles.cxx:530`): `StepShape_ShapeRepresentation` itself and
/// every entity derived from it (`DataExchange/TKDESTEP/StepShape`, plus
/// `StepVisual_TessellatedShapeRepresentation`).
pub(in crate::step) fn is_shape_representation_type(type_name: &str) -> bool {
    matches!(
        type_name,
        "SHAPE_REPRESENTATION"
            | "ADVANCED_BREP_SHAPE_REPRESENTATION"
            | "COMPOUND_SHAPE_REPRESENTATION"
            | "CSG_SHAPE_REPRESENTATION"
            | "EDGE_BASED_WIREFRAME_SHAPE_REPRESENTATION"
            | "FACETED_BREP_SHAPE_REPRESENTATION"
            | "GEOMETRICALLY_BOUNDED_SURFACE_SHAPE_REPRESENTATION"
            | "GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION"
            | "MANIFOLD_SURFACE_SHAPE_REPRESENTATION"
            | "NON_MANIFOLD_SURFACE_SHAPE_REPRESENTATION"
            | "POINT_REPRESENTATION"
            | "SHAPE_DIMENSION_REPRESENTATION"
            | "SHAPE_REPRESENTATION_WITH_PARAMETERS"
            | "TRANSITIONAL_SHAPE_REPRESENTATION"
            | "TESSELLATED_SHAPE_REPRESENTATION"
    )
}

/// `ProcessAsSurfaceStyleRendering` (`STEPConstruct_Styles.cxx:85-100`).
fn process_as_surface_style_rendering(
    records: &HashMap<usize, Record>,
    sses: usize,
    rendering_props: &mut RenderingProperties,
) -> bool {
    let Some(rec) = records.get(&sses) else { return false };
    // `theSSES.SurfaceStyleRendering()` — select case 7, matched by kind
    // (`StepVisual_SurfaceStyleElementSelect.hxx:47`), so both
    // `SURFACE_STYLE_RENDERING` and its subtype
    // `SURFACE_STYLE_RENDERING_WITH_PROPERTIES` land here.
    if rec.type_name != "SURFACE_STYLE_RENDERING"
        && rec.type_name != "SURFACE_STYLE_RENDERING_WITH_PROPERTIES"
    {
        return false;
    }
    // `occ::down_cast<StepVisual_SurfaceStyleRenderingWithProperties>(aSSR)`
    // (`cxx:92-96`).
    if rec.type_name != "SURFACE_STYLE_RENDERING_WITH_PROPERTIES" {
        return true;
    }
    rendering_props.init(records, sses);    // `return theRenderingProps.IsDefined()` (`cxx:99`).
    rendering_props.is_defined()
}

/// `ProcessAsSurfaceStyleBoundary` (`STEPConstruct_Styles.cxx:110-126`).
fn process_as_surface_style_boundary(
    records: &HashMap<usize, Record>,
    sses: usize,
    boundary_color: &mut Option<usize>,
) -> bool {
    let Some(rec) = records.get(&sses) else { return false };
    if rec.type_name != "SURFACE_STYLE_BOUNDARY" {
        return false;
    }
    // `aSSB->StyleOfBoundary()` (attribute 1,
    // `RWStepVisual_RWSurfaceStyleBoundary.cxx:41-50`).
    let Some(cs) = rec.args.first().and_then(|a| parse_ref(a)) else { return true };
    let Some(crec) = records.get(&cs) else { return true };
    if crec.type_name != "CURVE_STYLE" {
        return true;
    }
    // `aCS->CurveColour()` (attribute 4,
    // `RWStepVisual_RWCurveStyle.cxx:68-71`).
    *boundary_color = crec.args.get(3).and_then(|a| parse_ref(a));
    true
}

/// `ProcessAsSurfaceStyleFillArea` (`STEPConstruct_Styles.cxx:137-162`).
fn process_as_surface_style_fill_area(
    records: &HashMap<usize, Record>,
    sses: usize,
    side: SurfaceSide,
    surface_color: &mut Option<usize>,
) -> bool {
    let Some(rec) = records.get(&sses) else { return false };
    if rec.type_name != "SURFACE_STYLE_FILL_AREA" {
        return false;
    }
    let Some(fas) = rec.args.first().and_then(|a| parse_ref(a)) else { return true };
    let Some(fas_rec) = records.get(&fas) else { return true };
    if fas_rec.type_name != "FILL_AREA_STYLE" {
        return true;
    }
    // `for (int aFSSIndex = 1; aFSSIndex <= aFAS->NbFillStyles(); aFSSIndex++)`
    // (`cxx:150-161`).
    let fill_styles = fas_rec.args.get(1).map(|a| parse_ref_list(a)).unwrap_or_default();
    for fss in fill_styles {
        let Some(frec) = records.get(&fss) else { continue };
        if frec.type_name != "FILL_AREA_STYLE_COLOUR" {
            continue;
        }
        let Some(fc) = frec.args.get(1).and_then(|a| parse_ref(a)) else { continue };
        // `if (!aFASC.IsNull() && (theSurfaceColour.IsNull() ||
        // theSide != StepVisual_ssNegative))` (`cxx:157-161`): a negative
        // side color only fills an unset surface color (`trj3_s1-pe.stp`).
        if surface_color.is_none() || side != SurfaceSide::Negative {
            *surface_color = Some(fc);
        }
    }
    true
}

/// `StepVisual_SurfaceSide` enum (`StepVisual_SurfaceSide.hxx`), only the
/// member the fill-area rule tests.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SurfaceSide {
    Positive,
    Negative,
    Both,
}

/// `ProcessAsSurfaceStyleUsage` (`STEPConstruct_Styles.cxx:176-200`).
fn process_as_surface_style_usage(
    records: &HashMap<usize, Record>,
    pss: usize,
    surface_color: &mut Option<usize>,
    boundary_color: &mut Option<usize>,
    rendering_props: &mut RenderingProperties,
) -> bool {
    let Some(rec) = records.get(&pss) else { return false };
    // `thePSS.SurfaceStyleUsage()` — select case 3, matched by kind.
    if rec.type_name != "SURFACE_STYLE_USAGE" {
        return false;
    }
    // `aSSU->Side()` / `aSSU->Style()` — attributes 1 and 2
    // (`RWStepVisual_RWSurfaceStyleUsage.cxx:60-64`).
    let side = match rec.args.first().map(|s| s.trim()) {
        Some(".NEGATIVE.") => SurfaceSide::Negative,
        Some(".BOTH.") => SurfaceSide::Both,
        _ => SurfaceSide::Positive,
    };
    let Some(sss) = rec.args.get(1).and_then(|a| parse_ref(a)) else { return true };
    let Some(sss_rec) = records.get(&sss) else { return true };
    if sss_rec.type_name != "SURFACE_SIDE_STYLE" {
        return true;
    }
    // `for (int aSSESIndex = 1; aSSESIndex <= aSSS->NbStyles(); ++aSSESIndex)`
    // (`cxx:188-196`).
    let elements = sss_rec.args.get(1).map(|a| parse_ref_list(a)).unwrap_or_default();
    for sses in elements {
        if !process_as_surface_style_fill_area(records, sses, side, surface_color)
            && !process_as_surface_style_boundary(records, sses, boundary_color)
        {
            process_as_surface_style_rendering(records, sses, rendering_props);
        }
    }
    true
}

/// `ProcessAsCurveStyle` (`STEPConstruct_Styles.cxx:209-220`).
fn process_as_curve_style(
    records: &HashMap<usize, Record>,
    pss: usize,
    curve_color: &mut Option<usize>,
) -> bool {
    let Some(rec) = records.get(&pss) else { return false };
    // `thePSS.CurveStyle()` — select case 2, matched by kind.
    if rec.type_name != "CURVE_STYLE" {
        return false;
    }
    *curve_color = rec.args.get(3).and_then(|a| parse_ref(a));
    true
}
