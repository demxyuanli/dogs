//! `CreateSewedShape` (`cxx:5029-5244`) and `CreateOutputInformations`
//! (`cxx:5256-5363`).

use crate::abs::ShapeType;
use crate::brep_sweep::Quilt;
use crate::shape::TopoShape;

use super::shape_maps::{ShapeKeyedMap, ShapeSet};
use super::Sewing;

impl Sewing {
    /// `CreateSewedShape` (`cxx:5029-5244`): quilt the (reshaped) input shapes
    /// into shells, then assemble the result compound.
    pub(crate) fn create_sewed_shape(&mut self) {
        let mut quilt = Quilt::new();
        let is_local = self.shape.is_some();
        if is_local {
            let ns = self
                .shape
                .as_ref()
                .and_then(|s| self.reshape.apply_option(s))
                .unwrap_or_else(|| self.shape.clone().unwrap());
            quilt.add(&ns);
        }
        let keys: Vec<TopoShape> = self.old_shapes.keys();
        for key in keys {
            let Some(sh) = self.old_shapes.get(&key) else {
                continue;
            };
            if sh.is_null() {
                continue;
            }
            let applied = self.reshape.apply_option(&sh);
            let Some(applied) = applied else {
                continue;
            };
            self.old_shapes.bind(&key, applied.clone());
            if !is_local {
                quilt.add(&applied);
            }
        }
        let new_shape = quilt.shells();
        let mut numsh = 0usize;
        let mut old_shells = ShapeSet::new();
        let builder = crate::builder::TopoBuilder::new();
        let mut comp = TopoShape::new(ShapeType::Compound);
        for sh in sh_shape_children(&new_shape) {
            let mut has_edges = false;
            let mut sh = sh;
            if sh.shape_type() == ShapeType::Shell {
                if self.nonmanifold {
                    has_edges = !old_shells.contains(&sh);
                } else {
                    let mut face;
                    let mut numf = 0usize;
                    for f in crate::topo_tools::explorer_shapes(&sh, ShapeType::Face) {
                        face = f;
                        numf += 1;
                        if numf >= 2 {
                            break;
                        }
                        let _ = &face;
                    }
                    let mut faces = crate::topo_tools::explorer_shapes(&sh, ShapeType::Face);
                    let numf = faces.len().min(2);
                    if numf == 1 {
                        let f = faces.remove(0);
                        builder.add(&mut comp, &f);
                    } else if numf > 1 {
                        builder.add(&mut comp, &sh);
                    }
                    if numf > 0 {
                        numsh += 1;
                    }
                }
            } else if sh.shape_type() == ShapeType::Face {
                if self.nonmanifold {
                    let mut shell = TopoShape::new(ShapeType::Shell);
                    builder.add(&mut shell, &sh);
                    sh = shell;
                    has_edges = true;
                } else {
                    builder.add(&mut comp, &sh);
                    numsh += 1;
                }
            } else {
                builder.add(&mut comp, &sh);
                numsh += 1;
            }
            if has_edges {
                old_shells.add(&sh);
            }
        }
        // Process the collected shells (`cxx:4585-4676`).
        if self.nonmanifold {
            let nb_old_shells = old_shells.extent();
            if nb_old_shells == 1 {
                // Single shell - check for a single face.
                let sh = old_shells.key_at(1).unwrap().clone();
                let faces = crate::topo_tools::explorer_shapes(&sh, ShapeType::Face);
                let numf = faces.len().min(2);
                if numf == 1 {
                    builder.add(&mut comp, &faces[0]);
                } else if numf > 1 {
                    builder.add(&mut comp, &sh);
                }
                if numf > 0 {
                    numsh += 1;
                }
            } else if nb_old_shells > 0 {
                // Several shells should be merged (`cxx:4608-4676`).
                let mut index_merged = vec![false; nb_old_shells];
                let mut merged_count = 0usize;
                while merged_count < nb_old_shells {
                    let mut new_shell = TopoShape::new(ShapeType::Shell);
                    let mut shell_started = false;
                    let mut new_edges = ShapeSet::new();
                    for i in 1..=nb_old_shells {
                        if index_merged[i - 1] {
                            continue;
                        }
                        let shell = old_shells.key_at(i).unwrap().clone();
                        if !shell_started {
                            shell_started = true;
                            for child in sh_shape_children(&shell) {
                                builder.add(&mut new_shell, &child);
                            }
                            for edge in crate::topo_tools::explorer_shapes(&shell, ShapeType::Edge) {
                                new_edges.add(&edge);
                            }
                            index_merged[i - 1] = true;
                            merged_count += 1;
                        } else {
                            let mut has_shared_edge = false;
                            for edge in crate::topo_tools::explorer_shapes(&shell, ShapeType::Edge) {
                                if new_edges.contains(&edge) {
                                    has_shared_edge = true;
                                    break;
                                }
                            }
                            if has_shared_edge {
                                for edge in
                                    crate::topo_tools::explorer_shapes(&shell, ShapeType::Edge)
                                {
                                    new_edges.add(&edge);
                                }
                                for face in
                                    crate::topo_tools::explorer_shapes(&shell, ShapeType::Face)
                                {
                                    builder.add(&mut new_shell, &face);
                                }
                                index_merged[i - 1] = true;
                                merged_count += 1;
                            }
                        }
                    }
                    // Process the new shell (`cxx:4656-4675`).
                    let faces = crate::topo_tools::explorer_shapes(&new_shell, ShapeType::Face);
                    let numf = faces.len().min(2);
                    if numf == 1 {
                        builder.add(&mut comp, &faces[0]);
                    } else if numf > 1 {
                        builder.add(&mut comp, &new_shell);
                    }
                    if numf > 0 {
                        numsh += 1;
                    }
                }
            }
        }
        if numsh == 1 {
            // Extract the single component (`cxx:4677-4680`).
            let children = sh_shape_children(&comp);
            self.sewed_shape = children.into_iter().next();
        } else {
            self.sewed_shape = Some(comp);
        }
    }

    /// `CreateOutputInformations` (`cxx:5256-5363`).
    pub(crate) fn create_output_informations(&mut self) {
        // Construct edgeSections (`cxx:4702-4734`).
        let mut edge_sections: ShapeKeyedMap<Vec<TopoShape>> = ShapeKeyedMap::new();
        let bounds = self.bound_faces.keys();
        for bound in bounds {
            let mut lsect = self.bound_sections.get(&bound).unwrap_or_default();
            let applied = match self.reshape.apply_option(&bound) {
                Some(a) => a,
                None => continue,
            };
            for edge in crate::topo_tools::explorer_shapes(&applied, ShapeType::Edge) {
                let mut sec = bound.clone();
                for section in &lsect {
                    match self.reshape.apply_option(section) {
                        Some(applied_section)
                            if crate::topo_tools::is_same(&edge, &applied_section) =>
                        {
                            sec = section.clone();
                            break;
                        }
                        _ => {}
                    }
                }
                edge_sections.change(&edge).push(sec);
            }
            lsect.clear();
        }
        // Fill the maps of free, contiguous and multiple edges
        // (`cxx:4737-4757`).
        let entries: Vec<(TopoShape, Vec<TopoShape>)> = edge_sections
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        for (edge, list_section) in entries {
            if list_section.len() == 1 {
                if crate::brep_tool::BRepTool::is_degenerated(&crate::shape::Edge(edge.clone())) {
                    self.degenerated.add(&edge);
                } else {
                    self.free_edges.add(&edge);
                }
            } else if list_section.len() == 2 {
                self.contigous_edges.add(&edge, list_section);
            } else {
                self.multiple_edges.add(&edge);
            }
        }
        // Construct myContigSectBound (`cxx:4759-4780`).
        let bounds = self.bound_faces.keys();
        for bound in bounds {
            let Some(sections) = self.bound_sections.get(&bound) else {
                continue;
            };
            for section in &sections {
                if !self.merged_edges.contains(section) {
                    continue;
                }
                let Some(nedge) = self.reshape.apply_option(section) else {
                    continue;
                };
                if !crate::topo_tools::is_same(&bound, section) && self.contigous_edges.contains(&nedge)
                {
                    self.contig_sec_bound.bind(section, bound.clone());
                }
            }
        }
    }
}

/// `TopoDS_Iterator(shape, /*cumOri=*/false)` children.
fn sh_shape_children(shape: &TopoShape) -> Vec<TopoShape> {
    shape.tshape.read().expect("poisoned TShape lock").children.clone()
}

