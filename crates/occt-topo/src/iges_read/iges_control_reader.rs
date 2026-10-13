//! `IGESControl_Reader` (`IGESControl_Reader.cxx:52-121`).

use std::path::Path;

use crate::shape::TopoShape;

use super::iges_to_brep_actor::actor_recognize;
use super::xscontrol_reader::XsControlReader;

/// `read.iges.onlyvisible` default is 0 (`IGESData.cxx:172`), so the ctor
/// stores `theReadOnlyVisible = false` (`IGESControl_Reader.cxx:57-58`).
const READ_ONLY_VISIBLE: bool = false;

pub(in crate::iges_read) struct IgesControlReader {
    session: XsControlReader,
    roots_ready: bool,
    roots: Vec<usize>,
}

impl IgesControlReader {
    pub(in crate::iges_read) fn new() -> Self {
        Self {
            session: XsControlReader::new(),
            roots_ready: false,
            roots: Vec::new(),
        }
    }

    pub(in crate::iges_read) fn read_file(&mut self, path: &Path) -> Result<(), String> {
        self.roots_ready = false;
        self.roots.clear();
        self.session.read_file(path)
    }

    /// `IGESControl_Reader::NbRootsForTransfer` (`cxx:82-121`).
    pub(in crate::iges_read) fn nb_roots_for_transfer(&mut self) -> Result<usize, String> {
        if self.roots_ready {
            return Ok(self.roots.len());
        }
        self.roots_ready = true;
        self.roots.clear();
        let model = self
            .session
            .model
            .as_ref()
            .ok_or_else(|| "IGES model is null".to_string())?;
        for ent in model.entities.iter().skip(1) {
            if ent.shared || !actor_recognize(model, ent.index) {
                continue;
            }
            // Visible-only filter. Default `theReadOnlyVisible` is false, so
            // every recognized unshared entity is a root.
            if !READ_ONLY_VISIBLE || ent.blank_status == 0 {
                self.roots.push(ent.index);
            }
        }
        Ok(self.roots.len())
    }

    pub(in crate::iges_read) fn transfer_one_shape(&mut self) -> Result<TopoShape, String> {
        self.nb_roots_for_transfer()?;
        let roots = self.roots.clone();
        self.session.transfer_roots(&roots)?;
        self.session.one_shape()
    }
}

/// Read `path` and transfer its roots to one shape.
pub fn read_iges_shape(path: &str) -> Result<TopoShape, String> {
    let mut reader = IgesControlReader::new();
    reader.read_file(Path::new(path))?;
    reader.transfer_one_shape()
}
