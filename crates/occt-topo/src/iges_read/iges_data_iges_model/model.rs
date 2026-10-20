use crate::iges_read::iges_data_param_reader::IgesParam;

use super::entities::EntityData;
use super::parse::parse_own;
use super::pointers::pointer_targets;
use super::units::unit_flag_value;

#[derive(Clone, Debug)]
pub(in crate::iges_read) struct IgesEntity {
    pub index: usize,
    pub type_number: i32,
    pub form_number: i32,
    /// `IGESData_IGESEntity::BlankStatus` (`cxx:256-259`), low 4 bits of the
    /// two-digit blank field. `00` visible, `01` blanked.
    pub blank_status: i32,
    /// Directory field 7. Zero means no 124 matrix.
    pub transform: i32,
    pub data: EntityData,
    pub shared: bool,
}

#[derive(Clone, Debug)]
pub(in crate::iges_read) struct IgesModel {
    pub entities: Vec<IgesEntity>,
    /// `IGESData_GlobalSection::UnitValue`.
    pub unit_factor: f64,
    /// Global parameter 19, minimum resolution.
    pub resolution: f64,
}

impl IgesModel {
    pub(in crate::iges_read) fn entity(&self, index: usize) -> Result<&IgesEntity, String> {
        self.entities
            .get(index)
            .filter(|e| e.index == index)
            .ok_or_else(|| format!("IGES entity {index} is out of range"))
    }
}

pub(in crate::iges_read) struct DirLine {
    pub type_number: i32,
    pub form_number: i32,
    pub blank_status: i32,
    pub transform: i32,
    pub de_seq: i32,
}

pub(in crate::iges_read) fn build_model(
    global: &[IgesParam],
    directory: &[DirLine],
    params_by_de: &std::collections::HashMap<i32, Vec<IgesParam>>,
) -> Result<IgesModel, String> {
    // Global parameters are 1-based and include the two delimiter fields.
    // Parameter 14 is the unit flag, 19 is the minimum resolution.
    let unit_flag = match global.get(13) {
        Some(IgesParam::Int(n)) => *n,
        Some(IgesParam::Real(r)) => *r as i32,
        _ => 2,
    };
    let resolution = match global.get(18) {
        Some(IgesParam::Real(r)) => *r,
        Some(IgesParam::Int(n)) => *n as f64,
        _ => 1.0e-7,
    };
    let unit_factor = unit_flag_value(unit_flag);
    if unit_factor == 0.0 {
        return Err(format!("IGES unit flag {unit_flag} has no UnitFlagValue"));
    }

    let mut entities = Vec::with_capacity(directory.len() + 1);
    entities.push(IgesEntity {
        index: 0,
        type_number: 0,
        form_number: 0,
        blank_status: 0,
        transform: 0,
        data: EntityData::Raw { type_number: 0 },
        shared: false,
    });
    for (i, dir) in directory.iter().enumerate() {
        let index = i + 1;
        let tokens = params_by_de.get(&dir.de_seq).ok_or_else(|| {
            format!(
                "IGES entity {index} (DE {}) has no parameter data",
                dir.de_seq
            )
        })?;
        let data = parse_own(dir.type_number, dir.form_number, tokens)?;
        entities.push(IgesEntity {
            index,
            type_number: dir.type_number,
            form_number: dir.form_number,
            blank_status: dir.blank_status,
            transform: dir.transform,
            data,
            shared: false,
        });
    }

    let mut shared = vec![false; entities.len()];
    for ent in entities.iter().skip(1) {
        for p in pointer_targets(&ent.data) {
            if p < shared.len() {
                shared[p] = true;
            } else {
                return Err(format!("IGES pointer to entity {p} is out of range"));
            }
        }
    }
    for (i, flag) in shared.iter().enumerate().skip(1) {
        entities[i].shared = *flag;
    }

    Ok(IgesModel {
        entities,
        unit_factor,
        resolution,
    })
}
