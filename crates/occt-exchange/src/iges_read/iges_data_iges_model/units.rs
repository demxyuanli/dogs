/// `IGESData_BasicEditor::UnitFlagValue` (`IGESData_BasicEditor.cxx:492-520`)
/// divided by the default cascade unit 1 (`IGESData_GlobalSection::UnitValue`).
pub(in crate::iges_read) fn unit_flag_value(flag: i32) -> f64 {
    match flag {
        1 => 25.4,
        2 => 1.0,
        3 => 1.0,
        4 => 304.8,
        5 => 1_609_344.0,
        6 => 1000.0,
        7 => 1_000_000.0,
        8 => 0.0254,
        9 => 0.001,
        10 => 10.0,
        11 => 0.0000254,
        _ => 0.0,
    }
}

/// Directory sequence (1-based, odd for the first line of a pair) to the
/// 1-based entity rank. `IGESData_IGESReaderTool` uses `(1 + pointer) / 2`.
pub(in crate::iges_read) fn de_pointer_to_entity(de: i32) -> Result<Option<usize>, String> {
    if de == 0 {
        return Ok(None);
    }
    if de < 0 {
        return Err(format!(
            "negative IGES DE pointer {de} is not ported (IGESData_IGESReaderTool.cxx)"
        ));
    }
    Ok(Some(((de as usize) + 1) / 2))
}
