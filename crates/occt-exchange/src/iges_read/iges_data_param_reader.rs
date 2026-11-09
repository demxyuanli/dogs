//! `IGESData_ParamReader` field cursor.
//!
//! Own-parameter indexes are 1-based and skip the entity type, which sits at
//! index 0 (`IGESData_ParamReader.cxx:116-133`: `ParamValue(num)` is
//! `list[num + base]`, and `base` is the type parameter).

use super::iges_data_iges_model::de_pointer_to_entity;

#[derive(Clone, Debug)]
pub(in crate::iges_read) enum IgesParam {
    Void,
    Int(i32),
    Real(f64),
    Text(String),
}

pub(in crate::iges_read) struct ParamReader<'a> {
    params: &'a [IgesParam],
    curr: usize,
}

impl<'a> ParamReader<'a> {
    /// `params[0]` is the entity type. The cursor starts on the first own field.
    pub(in crate::iges_read) fn new(params: &'a [IgesParam]) -> Self {
        Self { params, curr: 1 }
    }

    pub(in crate::iges_read) fn remaining(&self) -> usize {
        self.params.len().saturating_sub(self.curr)
    }

    /// `IGESData_ParamReader::SetCurrentNumber` used to skip spline-surface
    /// arbitrary values (`IGESGeom_ToolSplineSurface.cxx:204-216`).
    pub(in crate::iges_read) fn skip(&mut self, n: usize) {
        self.curr = self.curr.saturating_add(n);
    }

    pub(in crate::iges_read) fn is_defined(&self) -> bool {
        // `IGESData_ParamReader::IsParamDefined` (`cxx:137-144`).
        self.curr < self.params.len() && !matches!(self.params[self.curr], IgesParam::Void)
    }

    fn take(&mut self) -> Result<&'a IgesParam, String> {
        if self.curr >= self.params.len() {
            return Err("IGES parameter number out of range".into());
        }
        let p = &self.params[self.curr];
        self.curr += 1;
        Ok(p)
    }

    pub(in crate::iges_read) fn read_integer(&mut self) -> Result<i32, String> {
        // `IGESData_ParamReader::ReadInteger`.
        match self.take()? {
            IgesParam::Int(n) => Ok(*n),
            IgesParam::Real(r) if r.fract() == 0.0 && r.abs() <= i32::MAX as f64 => Ok(*r as i32),
            _ => Err("expected IGES integer parameter".into()),
        }
    }

    /// `IGESData_ParamReader::ReadBoolean` (`cxx:395-435`): a void field reads
    /// as false and succeeds, an integer reads as `!= 0`, anything else fails.
    pub(in crate::iges_read) fn read_boolean(&mut self) -> Result<bool, String> {
        match self.take()? {
            IgesParam::Void => Ok(false),
            IgesParam::Int(n) => Ok(*n != 0),
            _ => Err("not an Integer (for Boolean)".into()),
        }
    }

    pub(in crate::iges_read) fn read_real(&mut self) -> Result<f64, String> {
        // `IGESData_ParamReader::ReadReal`.
        match self.take()? {
            IgesParam::Real(r) => Ok(*r),
            IgesParam::Int(n) => Ok(*n as f64),
            _ => Err("expected IGES real parameter".into()),
        }
    }

    pub(in crate::iges_read) fn read_xy(&mut self) -> Result<[f64; 2], String> {
        // `IGESData_ParamReader::ReadXY` of `CurrentList(1, 2)`.
        Ok([self.read_real()?, self.read_real()?])
    }

    /// `IGESData_ParamReader::ReadText` (`IGESBasic_ToolSubfigureDef.cxx:66`).
    pub(in crate::iges_read) fn read_text(&mut self) -> Result<String, String> {
        match self.take()? {
            IgesParam::Text(s) => Ok(s.clone()),
            _ => Err("XSTEP_209: IGESBasic_ToolSubfigureDef.cxx:66".into()),
        }
    }

    pub(in crate::iges_read) fn read_xyz(&mut self) -> Result<[f64; 3], String> {
        // `IGESData_ParamReader::ReadXYZ` of `CurrentList(1, 3)`.
        Ok([self.read_real()?, self.read_real()?, self.read_real()?])
    }

    /// `IGESData_ParamReader::DefinedElseSkip`: a void field is skipped and
    /// reported absent. A defined field stays under the cursor.
    pub(in crate::iges_read) fn defined_else_skip(&mut self) -> bool {
        if self.curr >= self.params.len() {
            return false;
        }
        if matches!(self.params[self.curr], IgesParam::Void) {
            self.curr += 1;
            return false;
        }
        true
    }

    /// Positive value is a directory sequence. Zero is a null pointer.
    pub(in crate::iges_read) fn read_entity(&mut self, allow_null: bool) -> Result<Option<usize>, String> {
        let de = self.read_integer()?;
        let ent = de_pointer_to_entity(de)?;
        if ent.is_none() && !allow_null {
            return Err(format!("null IGES entity pointer (DE {de})"));
        }
        Ok(ent)
    }
}
