//! `XSControl_Reader::ReadFile`, `TransferRoots`, and `OneShape`.

use std::path::Path;

use crate::builder::TopoBuilder;
use crate::shape::TopoShape;

use super::iges_data_iges_model::IgesModel;
use super::iges_select_work_library::work_library_read_file;
use super::iges_to_brep_actor::actor_transfer;

pub(in crate::iges_read) struct XsControlReader {
    pub model: Option<IgesModel>,
    pub shapes: Vec<TopoShape>,
}

impl XsControlReader {
    pub(in crate::iges_read) fn new() -> Self {
        Self {
            model: None,
            shapes: Vec::new(),
        }
    }

    /// `XSControl_Reader::ReadFile` (`XSControl_Reader.cxx:106-110`).
    pub(in crate::iges_read) fn read_file(&mut self, path: &Path) -> Result<(), String> {
        let model = work_library_read_file(path)?;
        self.model = Some(model);
        // `InitTransferReader(4)` clears the previous transfer.
        self.shapes.clear();
        Ok(())
    }

    /// `XSControl_Reader::TransferRoots` (`XSControl_Reader.cxx:273-297`).
    /// `roots` is the result of `NbRootsForTransfer`.
    pub(in crate::iges_read) fn transfer_roots(&mut self, roots: &[usize]) -> Result<usize, String> {
        self.shapes.clear();
        let produced = {
            let model = self
                .model
                .as_ref()
                .ok_or_else(|| "IGES transfer has no model".to_string())?;
            let mut produced = Vec::new();
            for &root in roots {
                if let Some(shape) = actor_transfer(model, root)? {
                    produced.push(shape);
                }
            }
            produced
        };
        let transferred = produced.len();
        self.shapes = produced;
        Ok(transferred)
    }

    /// `XSControl_Reader::OneShape` (`XSControl_Reader.cxx:330-349`).
    pub(in crate::iges_read) fn one_shape(&self) -> Result<TopoShape, String> {
        if self.shapes.is_empty() {
            return Err("IGES transfer produced no shape".into());
        }
        if self.shapes.len() == 1 {
            return Ok(self.shapes[0].clone());
        }
        let b = TopoBuilder::new();
        let compound = b.make_compound_of(&self.shapes);
        if compound.0.tshape.read().expect("lock").children.is_empty() {
            return Err("IGES transfer produced only null shapes".into());
        }
        Ok(compound.0)
    }
}
