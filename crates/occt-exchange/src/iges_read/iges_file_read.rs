//! `IGESFile_Read` (`IGESFile_Read.cxx:76-94`) loads through `igesread`.

use std::fs;
use std::path::Path;

use super::iges_data_iges_model::{build_model, IgesModel};
use super::igesread::igesread;

pub(in crate::iges_read) fn iges_file_read(path: &Path) -> Result<IgesModel, String> {
    let data = fs::read(path).map_err(|_| format!("File not found : {}", path.display()))?;
    let read = igesread(&data)?;
    build_model(&read.global, &read.directory, &read.params_by_de)
}
