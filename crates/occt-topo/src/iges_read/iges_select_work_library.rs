//! `IGESSelect_WorkLibrary::ReadFile` (`IGESSelect_WorkLibrary.cxx:73-99`).

use std::path::Path;

use super::iges_data_iges_model::IgesModel;
use super::iges_file_read::iges_file_read;

pub(in crate::iges_read) fn work_library_read_file(path: &Path) -> Result<IgesModel, String> {
    // status < 0 is "file not found"; status > 0 is a read error. Both fail
    // the load and null the model (`cxx:84-98`). `iges_file_read` returns
    // those as `Err`.
    iges_file_read(path)
}
