//! `IGESToBRep_CurveAndSurface::TransferGeometry` (`cxx:214-571`).
//!
//! 402 forms 1, 7, 14, and 15 are `IsKind(IGESBasic_Group)` (`cxx:369`).
//! `GroupWithoutBackP`, `OrderedGroup`, and `OrderedGroupWithoutBackP` inherit
//! `IGESBasic_Group`, so the later arm at `:446` is not reached.
//! Entities 408 and 308 are `:264-368`. The shared location tail is `:533-569`.
//! `HasTransf` uses `CompoundLocation`; this applies the entity's own 124, the
//! same gap as `finish_plane`. 402 form 1 does not run that tail.
//! A failed group member is tolerated: `cxx:411-427` (402) and `cxx:343-355`
//! (308) wrap each item transfer in `try/catch (Standard_Failure)`, nullify the
//! item and keep building. The port maps that to a skipped `Err`.

use std::collections::{HashMap, HashSet};

use occt_core::gp::{GpTrsf, GpVec};
use occt_core::toploc::TopLocLocation;

use crate::builder::TopoBuilder;
use crate::shape::TopoShape;

use super::iges_data_iges_model::{EntityData, IgesModel};
use super::iges_data_param_reader::ParamReader;
use super::iges_to_brep_basic_curve::{convert_location, entity_location};
use super::iges_to_brep_topo_curve::transfer_topo_curve;
use super::iges_to_brep_topo_surface::transfer_topo_surface;

/// The two OCCT binders `IGESToBRep_CurveAndSurface` keeps in the transfer
/// process (`cxx:575-658`).
///
/// `cache` / `result_set` are the `TransferBRep_ShapeBinder`: one result per
/// entity, present as soon as `SetShapeResult` runs even when the shape is null
/// (`SetResultPresent`, `TransferBRep_BinderOfShape.cxx:39-42`).
/// `list_cache` / `list_count` are the `TransferBRep_ShapeListBinder`: a
/// 1-based list of results for the 502 / 504 list entities, where a null item
/// still occupies its slot (`AddResult`, `TransferBRep_ShapeListBinder.cxx:56-59`).
pub(in crate::iges_read) struct Session<'a> {
    pub model: &'a IgesModel,
    pub cache: HashMap<usize, TopoShape>,
    /// Entities whose single result was bound, including a null one
    /// (`IGESToBRep_CurveAndSurface.cxx:605-612`, `:575-583`).
    pub result_set: HashSet<usize>,
    /// 1-based items added by `AddShapeResult` (`cxx:651-658`).
    pub list_cache: HashMap<(usize, usize), Option<TopoShape>>,
    /// `NbShapes()` of the list binder (`cxx:616-628`).
    pub list_count: HashMap<usize, usize>,
    /// `IGESToBRep_TopoSurface::TheULength` (`cxx:191`, `:732`).
    pub u_length: f64,
}

impl<'a> Session<'a> {
    pub(in crate::iges_read) fn new(model: &'a IgesModel) -> Self {
        Self {
            model,
            cache: HashMap::new(),
            result_set: HashSet::new(),
            list_cache: HashMap::new(),
            list_count: HashMap::new(),
            u_length: 1.0,
        }
    }

    /// `HasShapeResult` (`cxx:575-583`). `TransferBRep_ShapeBinder::SetResult`
    /// marks the result present even when the shape is null, which is why a
    /// failed `TransferFace` is not retried.
    pub(in crate::iges_read) fn has_shape_result(&self, index: usize) -> bool {
        self.result_set.contains(&index)
    }

    /// `SetShapeResult` (`cxx:605-612`).
    pub(in crate::iges_read) fn set_shape_result(&mut self, index: usize, shape: Option<TopoShape>) {
        self.result_set.insert(index);
        match shape {
            Some(s) => {
                self.cache.insert(index, s);
            }
            None => {
                self.cache.remove(&index);
            }
        }
    }

    /// `GetShapeResult` (`cxx:587-604`).
    pub(in crate::iges_read) fn get_shape_result(&self, index: usize) -> Option<TopoShape> {
        self.cache.get(&index).cloned()
    }

    /// `NbShapeResult` (`cxx:616-628`).
    pub(in crate::iges_read) fn nb_shape_result(&self, index: usize) -> usize {
        self.list_count.get(&index).copied().unwrap_or(0)
    }

    /// `AddShapeResult` (`cxx:651-658`). A `None` item is the OCCT
    /// "add null shape to avoid shift of indexing" slot.
    pub(in crate::iges_read) fn add_shape_result(&mut self, index: usize, shape: Option<TopoShape>) {
        let count = self.list_count.entry(index).or_insert(0);
        *count += 1;
        self.list_cache.insert((index, *count), shape);
    }

    /// `GetShapeResult(start, num)` (`cxx:630-648`). OCCT indexes 1-based and
    /// reads nothing above `NbShapes()`; `num == 0` is out of range there too.
    pub(in crate::iges_read) fn list_shape_result(&self, index: usize, num: usize) -> Option<TopoShape> {
        if num == 0 || num > self.nb_shape_result(index) {
            return None;
        }
        self.list_cache.get(&(index, num)).cloned().flatten()
    }
}

/// `IGESToBRep::IsBasicCurve` (`IGESToBRep.cxx:128-159`).
pub(in crate::iges_read) fn is_basic_curve(type_number: i32) -> bool {
    matches!(type_number, 126 | 110 | 100 | 104 | 106 | 112)
}

/// `IGESToBRep::IsBasicSurface` (`IGESToBRep.cxx:167-204`).
pub(in crate::iges_read) fn is_basic_surface(type_number: i32) -> bool {
    matches!(type_number, 128 | 114 | 190 | 192 | 194 | 196 | 198)
}

/// `IGESToBRep::IsTopoCurve` (`IGESToBRep.cxx:212-243`).
pub(in crate::iges_read) fn is_topo_curve(type_number: i32) -> bool {
    is_basic_curve(type_number) || matches!(type_number, 102 | 142 | 141 | 116 | 130)
}

/// `IGESToBRep::IsTopoSurface` (`IGESToBRep.cxx:251-310`).
///
/// The 402 form 9 arm (`cxx:291-309`) needs the entity itself: an
/// `IGESBasic_SingleParent` whose single parent and every child is an
/// `IGESGeom_Plane` (IGES type 402 form 9, a perforated face) is a topo surface.
pub(in crate::iges_read) fn is_topo_surface(model: &IgesModel, index: usize) -> bool {
    // `cxx:254-257`: a null handle is false.
    let Ok(ent) = model.entity(index) else {
        return false;
    };
    if is_basic_surface(ent.type_number) {
        return true;
    }
    if matches!(ent.type_number, 144 | 120 | 122 | 118 | 108 | 143 | 140) {
        return true;
    }
    // `cxx:291-309`, SingleParent special case (Perforated Face: contains only
    // PLANEs).
    let EntityData::SingleParent { parent, children } = &ent.data else {
        return false;
    };
    let is_plane = |i: usize| matches!(model.entity(i), Ok(e) if e.type_number == 108);
    // `cxx:297` dereferences `SingleParent()` without testing it; a null parent
    // handle is a null `IGESGeom_Plane`, which no consumer accepts
    // (`IGESToBRep_TopoSurface.cxx:1543-1551`), so the class test is false here.
    let Some(parent) = *parent else {
        return false;
    };
    if !is_plane(parent) {
        return false;
    }
    children.iter().all(|c| matches!(c, Some(i) if is_plane(*i)))
}

/// `IGESToBRep::IsCurveAndSurface` (`IGESToBRep.cxx:101-120`).
pub(in crate::iges_read) fn is_curve_and_surface(model: &IgesModel, index: usize) -> bool {
    // `cxx:105-108`: a null handle is false.
    let Ok(ent) = model.entity(index) else {
        return false;
    };
    is_topo_curve(ent.type_number)
        || is_topo_surface(model, index)
        || super::igessolid::is_brep_entity(ent.type_number)
}

pub(in crate::iges_read) fn transfer_geometry(
    session: &mut Session<'_>,
    index: usize,
) -> Result<Option<TopoShape>, String> {
    if let Some(hit) = session.cache.get(&index) {
        return Ok(Some(hit.clone()));
    }
    let ent = session.model.entity(index)?;
    let type_number = ent.type_number;
    let form_number = ent.form_number;
    // `read.iges.onlyvisible` is 0, so the blank-status returns at
    // `IGESToBRep_CurveAndSurface.cxx:246-248` and `:371-373` are not taken.
    if is_curve_and_surface(session.model, index) {
        let shape = transfer_curve_and_surface(session, index)?;
        return match shape {
            Some(s) => {
                session.cache.insert(index, s.clone());
                Ok(Some(s))
            }
            None => Ok(None),
        };
    }
    if type_number == 408 {
        return transfer_singular_subfigure(session, index);
    }
    if type_number == 308 {
        return transfer_subfigure_def(session, index);
    }
    if type_number == 402 && matches!(form_number, 1 | 7 | 14 | 15) {
        return transfer_group(session, index);
    }
    Ok(None)
}

fn transfer_curve_and_surface(
    session: &mut Session<'_>,
    index: usize,
) -> Result<Option<TopoShape>, String> {
    // `IGESToBRep_CurveAndSurface::TransferCurveAndSurface` (`cxx:160-209`).
    let type_number = session.model.entity(index)?.type_number;
    if is_topo_curve(type_number) {
        return Ok(Some(transfer_topo_curve(session, index)?));
    }
    if is_topo_surface(session.model, index) {
        return Ok(Some(transfer_topo_surface(session, index)?));
    }
    if super::igessolid::is_brep_entity(type_number) {
        // `TransferBRepEntity` (`IGESToBRep_BRepEntity.cxx:110-136`) fails the
        // 502 / 504 / 508 list entities with `IGES_1005` and returns a null
        // shape, so those entities transfer to nothing.
        return super::brep_entity::transfer_brep_entity(session, index);
    }
    Err(format!(
        "UNPORTED entity {type_number}: IGESToBRep_CurveAndSurface.cxx:185-189"
    ))
}

fn transfer_group(session: &mut Session<'_>, index: usize) -> Result<Option<TopoShape>, String> {
    // `IGESToBRep_CurveAndSurface.cxx:369-445`.
    let members = match &session.model.entity(index)?.data {
        EntityData::Group { entities } => entities.clone(),
        _ => return Err("IGES 402 has no entity list".into()),
    };
    if members.is_empty() {
        return Ok(None);
    }
    let mut shapes = Vec::new();
    // `ProblemInGroup` (`cxx:388`, set at `:430`). OCCT only uses it to send
    // `IGES_1030` as a warning (`cxx:439-443`); the port has no message channel,
    // so the flag is reported through `IGES_TRACE_GROUP` instead.
    let mut problem_in_group = false;
    let mut skipped: Vec<String> = Vec::new();
    for child in members {
        // `cxx:406-427`: a cached result is reused, otherwise the item is
        // transferred inside `try { OCC_CATCH_SIGNALS; TransferGeometry } catch
        // (Standard_Failure) { item.Nullify(); SendFail(IGES_1015); }`. The port
        // represents that failure as `Err`, and `transfer_geometry` already
        // serves the cache, so the child transfer is wrapped here instead of
        // propagated (`cxx:429-431`): an item that yields nothing sets
        // `ProblemInGroup` and is skipped, and the compound is still returned.
        let item = match transfer_geometry(session, child) {
            Ok(Some(item)) => item,
            Ok(None) => {
                problem_in_group = true;
                continue;
            }
            Err(reason) => {
                problem_in_group = true;
                skipped.push(reason);
                continue;
            }
        };
        shapes.push(item);
    }
    if problem_in_group && std::env::var("IGES_TRACE_GROUP").is_ok() {
        eprintln!(
            "TRACEGROUP index={index} problem_in_group kept={} skipped={skipped:?}",
            shapes.len()
        );
    }
    let compound = TopoBuilder::new().make_compound_of(&shapes).0;
    session.cache.insert(index, compound.clone());
    Ok(Some(compound))
}

pub(in crate::iges_read) fn parse_subfigure_def(pr: &mut ParamReader<'_>) -> Result<EntityData, String> {
    // `IGESBasic_ToolSubfigureDef.cxx:61-79`.
    let depth = pr.read_integer()?;
    let name = pr.read_text()?;
    let nb = pr.read_integer()?;
    if nb < 0 {
        return Err("XSTEP_210: IGESBasic_ToolSubfigureDef.cxx:75".into());
    }
    let mut entities = Vec::with_capacity(nb as usize);
    for _ in 0..nb {
        entities.push(pr.read_entity(true)?);
    }
    Ok(EntityData::SubfigureDef {
        depth,
        name,
        entities,
    })
}

pub(in crate::iges_read) fn parse_singular_subfigure(
    pr: &mut ParamReader<'_>,
) -> Result<EntityData, String> {
    // `IGESBasic_ToolSingularSubfigure.cxx:58-112`.
    let subfigure = pr
        .read_entity(false)?
        .ok_or_else(|| "XSTEP_212: IGESBasic_ToolSingularSubfigure.cxx:63".to_string())?;
    let translation = pr.read_xyz()?;
    let (has_scale, scale) = if pr.defined_else_skip() {
        (true, pr.read_real()?)
    } else {
        (false, 1.0)
    };
    Ok(EntityData::SingularSubfigure {
        subfigure,
        translation,
        has_scale,
        scale,
    })
}

fn transfer_singular_subfigure(
    session: &mut Session<'_>,
    index: usize,
) -> Result<Option<TopoShape>, String> {
    // `cxx:264-305`, then the location tail `cxx:533-569`.
    // `read.iges.onlyvisible` is 0, so the blank return at `:266-268` is not taken.
    let (subfigure, translation, has_scale, scale) = match &session.model.entity(index)?.data {
        EntityData::SingularSubfigure {
            subfigure,
            translation,
            has_scale,
            scale,
        } => (*subfigure, *translation, *has_scale, *scale),
        _ => return Err("IGES 408 has no singular-subfigure parameters".into()),
    };
    let Some(mut shape) = transfer_geometry(session, subfigure)? else {
        return Ok(None);
    };
    let unit = session.model.unit_factor;
    let mut placed = GpTrsf::identity();
    placed.set_translation_part(&GpVec::new(
        translation[0] * unit,
        translation[1] * unit,
        translation[2] * unit,
    ));
    if has_scale {
        placed
            .set_scale_factor(scale)
            .map_err(|e| format!("IGES_1015: IGESToBRep_CurveAndSurface.cxx:250 {e}"))?;
    }
    shape = place_result(session, index, shape, Some(placed))?;
    session.cache.insert(index, shape.clone());
    Ok(Some(shape))
}

fn transfer_subfigure_def(
    session: &mut Session<'_>,
    index: usize,
) -> Result<Option<TopoShape>, String> {
    // `cxx:308-368`. `NbEntities < 1` returns a null shape (`cxx:314-318`).
    let entities = match &session.model.entity(index)?.data {
        EntityData::SubfigureDef { entities, .. } => entities.clone(),
        _ => return Err("IGES 308 has no subfigure-definition parameters".into()),
    };
    if entities.is_empty() {
        return Ok(None);
    }
    let mut shapes = Vec::new();
    for child in entities {
        let Some(child) = child else {
            continue;
        };
        // `cxx:338-355`: a cached result is reused, otherwise the item is
        // transferred inside `try { OCC_CATCH_SIGNALS; TransferGeometry } catch
        // (Standard_Failure) { item.Nullify(); SendFail(IGES_1015); }`. As in
        // `transfer_group`, the port turns that failure into a skipped item
        // (`cxx:356-364`) instead of propagating it.
        match transfer_geometry(session, child) {
            Ok(Some(item)) => shapes.push(item),
            Ok(None) | Err(_) => continue,
        }
    }
    let mut shape = TopoBuilder::new().make_compound_of(&shapes).0;
    shape = place_result(session, index, shape, None)?;
    session.cache.insert(index, shape.clone());
    Ok(Some(shape))
}

/// Location tail (`cxx:533-569`). `t408` is set only for entity 408.
fn place_result(
    session: &Session<'_>,
    index: usize,
    mut shape: TopoShape,
    t408: Option<GpTrsf>,
) -> Result<TopoShape, String> {
    if let Some(g) = entity_location(session.model, index)? {
        if let Some(mut tr) = convert_location(1.0e-4, &g, session.model.unit_factor) {
            if let Some(extra) = t408.as_ref() {
                let tra = tr.translation_part();
                let trans = extra.translation_part();
                tr.set_translation_part(&GpVec::new(
                    tra.x() + trans.x(),
                    tra.y() + trans.y(),
                    tra.z() + trans.z(),
                ));
                let sc = tr.scale_factor() * extra.scale_factor();
                tr.set_scale_factor(sc)
                    .map_err(|e| format!("IGES_1015: IGESToBRep_CurveAndSurface.cxx:250 {e}"))?;
            }
            shape.move_location(&TopLocLocation::composed(&tr, TopLocLocation::identity()));
        }
    } else if let Some(extra) = t408.as_ref() {
        shape.move_location(&TopLocLocation::composed(extra, TopLocLocation::identity()));
    }
    Ok(shape)
}
