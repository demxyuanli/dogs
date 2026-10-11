//! `BRepSweep_Tool` and `BRepSweep_NumLinearRegularSweep::Shape`.
//!
//! Sources: `BRepSweep_Tool.cxx`, `BRepSweep_NumLinearRegularSweep.cxx`.

use std::cell::RefCell;
use std::collections::HashSet;

use crate::abs::{Orientation, ShapeType};
use crate::builder::TopoBuilder;
use crate::iterator::cumulated_children;
use crate::shape::TopoShape;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools_full::is_same;

use super::num_shape::{NumShape, NumShapeIterator, NumShapeTool};

/// `BRepSweep_Tool` (`BRepSweep_Tool.cxx`): `TopExp::MapShapes` order.
#[derive(Debug, Clone)]
pub struct GenTool {
    /// 1-based. Slot 0 is unused.
    map: Vec<TopoShape>,
}

impl GenTool {
    pub fn new(shape: &TopoShape) -> Self {
        let mut map = vec![TopoShape::new(ShapeType::Shape)];
        let mut seen = HashSet::new();
        fn walk(s: &TopoShape, map: &mut Vec<TopoShape>, seen: &mut HashSet<usize>) {
            let id = std::sync::Arc::as_ptr(&s.tshape) as usize;
            if !seen.insert(id) {
                return;
            }
            map.push(s.clone());
            for c in cumulated_children(s) {
                walk(&c, map, seen);
            }
        }
        walk(shape, &mut map, &mut seen);
        Self { map }
    }

    pub fn nb_shapes(&self) -> usize {
        self.map.len() - 1
    }

    pub fn index(&self, shape: &TopoShape) -> usize {
        self.map
            .iter()
            .position(|s| is_same(s, shape))
            .unwrap_or(0)
    }

    pub fn shape(&self, index: usize) -> TopoShape {
        self.map[index].clone()
    }
}

/// Grid of generated shapes (`BRepSweep_NumLinearRegularSweep` fields).
pub struct SweepTables {
    pub gen: TopoShape,
    pub dir: NumShape,
    pub gen_tool: GenTool,
    pub dir_tool: NumShapeTool,
    pub shapes: RefCell<Vec<Vec<Option<TopoShape>>>>,
    pub built: RefCell<Vec<Vec<bool>>>,
    pub used: RefCell<Vec<Vec<bool>>>,
    pub error: RefCell<Option<String>>,
}

impl SweepTables {
    pub fn new(gen: TopoShape, dir: NumShape) -> Self {
        let gen_tool = GenTool::new(&gen);
        let dir_tool = NumShapeTool::new(dir);
        let ng = gen_tool.nb_shapes();
        let nd = dir_tool.nb_shapes() as usize;
        let row = vec![None; nd + 1];
        let brow = vec![false; nd + 1];
        Self {
            gen,
            dir,
            gen_tool,
            dir_tool,
            shapes: RefCell::new(vec![row.clone(); ng + 1]),
            built: RefCell::new(vec![brow.clone(); ng + 1]),
            used: RefCell::new(vec![brow; ng + 1]),
            error: RefCell::new(None),
        }
    }

    pub(crate) fn fail(&self, msg: &str) {
        *self.error.borrow_mut() = Some(msg.to_string());
    }
}

/// Geometric hooks `BRepSweep_Rotation` supplies to the grid.
pub trait SweepOps {
    fn tables(&self) -> &SweepTables;
    fn make_empty_vertex(&self, gen: &TopoShape, dir: NumShape) -> TopoShape;
    fn make_empty_directing_edge(&self, gen: &TopoShape, dir: NumShape) -> TopoShape;
    fn make_empty_generating_edge(&self, gen: &TopoShape, dir: NumShape) -> TopoShape;
    fn make_empty_face(&self, gen: &TopoShape, dir: NumShape) -> TopoShape;
    fn set_parameters(
        &self,
        new_face: &TopoShape,
        new_vertex: &TopoShape,
        gen_f: &TopoShape,
        gen_v: &TopoShape,
        dir: NumShape,
    );
    fn set_directing_parameter(
        &self,
        new_edge: &TopoShape,
        new_vertex: &TopoShape,
        gen_v: &TopoShape,
        dir_e: NumShape,
        dir_v: NumShape,
    );
    fn set_generating_parameter(
        &self,
        new_edge: &TopoShape,
        new_vertex: &TopoShape,
        gen_e: &TopoShape,
        gen_v: &TopoShape,
        dir: NumShape,
    );
    fn set_pcurve(
        &self,
        new_face: &TopoShape,
        new_edge: &TopoShape,
        gen_f: &TopoShape,
        gen_e: &TopoShape,
        dir: NumShape,
        orien: Orientation,
    );
    fn set_generating_pcurve(
        &self,
        new_face: &TopoShape,
        new_edge: &TopoShape,
        gen_e: &TopoShape,
        dir_e: NumShape,
        dir_v: NumShape,
        orien: Orientation,
    );
    fn set_directing_pcurve(
        &self,
        new_face: &TopoShape,
        new_edge: &TopoShape,
        gen_e: &TopoShape,
        gen_v: &TopoShape,
        dir_e: NumShape,
        orien: Orientation,
    );
    fn direct_solid(&self, gen: &TopoShape, dir: NumShape) -> Orientation;
    fn ggd_is_to_add(
        &self,
        new_shape: &TopoShape,
        new_sub: &TopoShape,
        gen: &TopoShape,
        sub_gen: &TopoShape,
        dir: NumShape,
    ) -> bool;
    fn gdd_is_to_add(
        &self,
        new_shape: &TopoShape,
        new_sub: &TopoShape,
        gen: &TopoShape,
        dir: NumShape,
        sub_dir: NumShape,
    ) -> bool;
    fn separated_wires(
        &self,
        new_shape: &TopoShape,
        new_sub: &TopoShape,
        gen: &TopoShape,
        sub_gen: &TopoShape,
        dir: NumShape,
    ) -> bool;
    fn split_shell(&self, shell: &TopoShape) -> TopoShape;
    fn set_continuity(&self, gen: &TopoShape, dir: NumShape);
    fn has_shape(&self, gen: &TopoShape, dir: NumShape) -> bool;
    fn is_invariant(&self, gen: &TopoShape) -> bool;
}

pub fn shape_of(host: &dyn SweepOps, gen: &TopoShape, dir: NumShape) -> Option<TopoShape> {
    let tables = host.tables();
    if tables.error.borrow().is_some() {
        return None;
    }
    let i_gen = tables.gen_tool.index(gen);
    let i_dir = tables.dir_tool.index(dir) as usize;
    if i_gen == 0 || i_dir == 0 {
        return None;
    }
    if tables.built.borrow()[i_gen][i_dir] {
        return tables.shapes.borrow()[i_gen][i_dir].clone();
    }
    let gen_ty = tables.gen_tool.shape(i_gen).shape_type();
    let dir_ty = dir.ty;
    let builder = TopoBuilder::new();
    if dir_ty == ShapeType::Vertex {
        let made = match gen_ty {
            ShapeType::Vertex => host.make_empty_vertex(gen, dir),
            ShapeType::Edge => host.make_empty_generating_edge(gen, dir),
            ShapeType::Wire => TopoShape::new(ShapeType::Wire),
            ShapeType::Face => host.make_empty_face(gen, dir),
            ShapeType::Shell => TopoShape::new(ShapeType::Shell),
            ShapeType::Compound => TopoShape::new(ShapeType::Compound),
            ShapeType::Solid | ShapeType::CompSolid => {
                tables.fail("Solids are not Processed");
                return None;
            }
            ShapeType::Shape => {
                tables.fail("Unknown Shape");
                return None;
            }
        };
        store(tables, i_gen, i_dir, made);
        let mut b_gen = gen.clone();
        b_gen.set_orientation(Orientation::Forward);
        for sub in cumulated_children(&b_gen) {
            let or = sub.orientation();
            if !host.has_shape(&sub, dir) {
                continue;
            }
            let Some(new_shape) = shape_of(host, &sub, dir) else {
                continue;
            };
            let i_new_gen = tables.gen_tool.index(&sub);
            if !host.ggd_is_to_add(&cell(tables, i_gen, i_dir), &new_shape, gen, &sub, dir) {
                continue;
            }
            add_into_cell(&builder, tables, i_gen, i_dir, &new_shape, or);
            mark_used(tables, i_new_gen, i_dir);
            let sub_ty = sub.shape_type();
            if gen_ty == ShapeType::Face {
                if sub_ty == ShapeType::Vertex {
                    host.set_parameters(&cell(tables, i_gen, i_dir), &new_shape, gen, &sub, dir);
                } else if sub_ty == ShapeType::Edge {
                    host.set_pcurve(&cell(tables, i_gen, i_dir), &new_shape, gen, &sub, dir, or);
                } else if sub_ty == ShapeType::Wire {
                    let mut c_gen = sub.clone();
                    c_gen.set_orientation(Orientation::Forward);
                    for subsub in cumulated_children(&c_gen) {
                        let pr = subsub.orientation();
                        if host.has_shape(&subsub, dir) {
                            if let Some(edge) = shape_of(host, &subsub, dir) {
                                host.set_pcurve(
                                    &cell(tables, i_gen, i_dir),
                                    &edge,
                                    gen,
                                    &subsub,
                                    dir,
                                    pr,
                                );
                            }
                        }
                    }
                }
            } else if gen_ty == ShapeType::Edge {
                host.set_generating_parameter(
                    &cell(tables, i_gen, i_dir),
                    &new_shape,
                    &b_gen,
                    &sub,
                    dir,
                );
            }
        }
    } else if dir_ty == ShapeType::Edge {
        let mut new_wire = TopoShape::new(ShapeType::Wire);
        let mut new_shell = TopoShape::new(ShapeType::Shell);
        let mut wire_seq: Vec<TopoShape> = Vec::new();
        let mut sepwires = false;
        let made = match gen_ty {
            ShapeType::Vertex => host.make_empty_directing_edge(gen, dir),
            ShapeType::Edge => host.make_empty_face(gen, dir),
            ShapeType::Wire => TopoShape::new(ShapeType::Shell),
            ShapeType::Face => {
                let solid = TopoShape::new(ShapeType::Solid);
                solid
            }
            ShapeType::Shell => TopoShape::new(ShapeType::CompSolid),
            ShapeType::Compound => TopoShape::new(ShapeType::Compound),
            ShapeType::Solid | ShapeType::CompSolid => {
                tables.fail("Solids are not Processed");
                return None;
            }
            ShapeType::Shape => {
                tables.fail("Unknown Shape");
                return None;
            }
        };
        store(tables, i_gen, i_dir, made);
        let mut b_gen = gen.clone();
        b_gen.set_orientation(Orientation::Forward);
        for sub in cumulated_children(&b_gen) {
            if !host.has_shape(&sub, dir) {
                continue;
            }
            let Some(new_shape) = shape_of(host, &sub, dir) else {
                continue;
            };
            let i_new_gen = tables.gen_tool.index(&sub);
            if !host.ggd_is_to_add(&cell(tables, i_gen, i_dir), &new_shape, gen, &sub, dir) {
                continue;
            }
            let sub_ty = sub.shape_type();
            if gen_ty == ShapeType::Edge {
                let or = sub.orientation();
                if host.separated_wires(&cell(tables, i_gen, i_dir), &new_shape, gen, &sub, dir) {
                    sepwires = true;
                    let mut wi = TopoShape::new(ShapeType::Wire);
                    add_oriented(&builder, &mut wi, &new_shape, or);
                    mark_used(tables, i_new_gen, i_dir);
                    wi.set_closed(shape_is_closed(&wi));
                    wire_seq.push(wi);
                } else {
                    add_oriented(&builder, &mut new_wire, &new_shape, or);
                    mark_used(tables, i_new_gen, i_dir);
                }
                host.set_directing_pcurve(
                    &cell(tables, i_gen, i_dir),
                    &new_shape,
                    &b_gen,
                    &sub,
                    dir,
                    or,
                );
            } else if gen_ty == ShapeType::Wire {
                let or = sub.orientation();
                add_into_cell(&builder, tables, i_gen, i_dir, &new_shape, or);
                mark_used(tables, i_new_gen, i_dir);
            } else if gen_ty == ShapeType::Face {
                let or = sub.orientation();
                if sub_ty == ShapeType::Wire {
                    for lt in cumulated_children(&new_shape) {
                        let composed = Orientation::compose(lt.orientation(), or);
                        add_oriented(&builder, &mut new_shell, &lt, composed);
                    }
                } else if sub_ty == ShapeType::Edge {
                    add_oriented(&builder, &mut new_shell, &new_shape, or);
                    mark_used(tables, i_new_gen, i_dir);
                }
            } else if gen_ty == ShapeType::Shell || gen_ty == ShapeType::Compound {
                add_into_cell(&builder, tables, i_gen, i_dir, &new_shape, Orientation::Forward);
                mark_used(tables, i_new_gen, i_dir);
            } else {
                let or = sub.orientation();
                add_into_cell(&builder, tables, i_gen, i_dir, &new_shape, or);
                mark_used(tables, i_new_gen, i_dir);
            }
        }
        let mut kt = NumShapeIterator::new(dir);
        while kt.more() {
            let sub_dir = kt.value();
            if host.has_shape(gen, sub_dir) {
                if let Some(new_shape) = shape_of(host, gen, sub_dir) {
                    let i_new_dir = tables.dir_tool.index(sub_dir) as usize;
                    if host.gdd_is_to_add(
                        &cell(tables, i_gen, i_dir),
                        &new_shape,
                        gen,
                        dir,
                        sub_dir,
                    ) {
                        if gen_ty == ShapeType::Edge {
                            let or = kt.orientation().reversed();
                            add_oriented(&builder, &mut new_wire, &new_shape, or);
                            mark_used(tables, i_gen, i_new_dir);
                            host.set_generating_pcurve(
                                &cell(tables, i_gen, i_dir),
                                &new_shape,
                                gen,
                                dir,
                                sub_dir,
                                or,
                            );
                        } else if gen_ty == ShapeType::Vertex {
                            let or = kt.orientation();
                            add_into_cell(&builder, tables, i_gen, i_dir, &new_shape, or);
                            mark_used(tables, i_gen, i_new_dir);
                            host.set_directing_parameter(
                                &cell(tables, i_gen, i_dir),
                                &new_shape,
                                gen,
                                dir,
                                sub_dir,
                            );
                        } else if gen_ty == ShapeType::Face {
                            let or = kt.orientation();
                            add_oriented(&builder, &mut new_shell, &new_shape, or);
                            mark_used(tables, i_gen, i_new_dir);
                        }
                    }
                }
            }
            kt.next();
        }
        if gen_ty == ShapeType::Edge {
            if sepwires {
                for wi in &wire_seq {
                    add_into_cell(&builder, tables, i_gen, i_dir, wi, wi.orientation());
                }
            } else {
                new_wire.set_closed(shape_is_closed(&new_wire));
                add_into_cell(&builder, tables, i_gen, i_dir, &new_wire, new_wire.orientation());
            }
            tables.built.borrow_mut()[i_gen][i_dir] = true;
            host.set_continuity(gen, dir);
        }
        if gen_ty == ShapeType::Wire {
            host.set_continuity(gen, dir);
        }
        if gen_ty == ShapeType::Face {
            new_shell.set_closed(shape_is_closed(&new_shell));
            let temp = host.split_shell(&new_shell);
            let shell_ori = host.direct_solid(gen, dir);
            let kids = cumulated_children(&temp);
            if kids.len() > 1 {
                for lt in &kids {
                    add_into_cell(&builder, tables, i_gen, i_dir, lt, shell_ori);
                }
            } else {
                add_into_cell(&builder, tables, i_gen, i_dir, &new_shell, shell_ori);
            }
        }
    } else if dir_ty == ShapeType::Wire {
        let made = match gen_ty {
            ShapeType::Vertex => TopoShape::new(ShapeType::Wire),
            ShapeType::Edge | ShapeType::Wire => TopoShape::new(ShapeType::Shell),
            ShapeType::Face | ShapeType::Shell => TopoShape::new(ShapeType::CompSolid),
            ShapeType::Compound => TopoShape::new(ShapeType::Compound),
            ShapeType::Solid | ShapeType::CompSolid => {
                tables.fail("Solids are not Processed");
                return None;
            }
            ShapeType::Shape => {
                tables.fail("Unknown Shape");
                return None;
            }
        };
        store(tables, i_gen, i_dir, made);
        let mut kt = NumShapeIterator::new(dir);
        while kt.more() {
            let sub_dir = kt.value();
            if host.has_shape(gen, sub_dir) {
                let i_new_dir = tables.dir_tool.index(sub_dir) as usize;
                let or = kt.orientation();
                if let Some(new_shape) = shape_of(host, gen, sub_dir) {
                    add_into_cell(&builder, tables, i_gen, i_dir, &new_shape, or);
                    mark_used(tables, i_gen, i_new_dir);
                }
            }
            kt.next();
        }
    }
    tables.built.borrow_mut()[i_gen][i_dir] = true;
    if let Some(s) = tables.shapes.borrow_mut()[i_gen][i_dir].as_mut() {
        if s.shape_type() == ShapeType::Wire || s.shape_type() == ShapeType::Shell {
            let closed = shape_is_closed(s);
            s.set_closed(closed);
        }
    }
    tables.shapes.borrow()[i_gen][i_dir].clone()
}

pub fn shape_result(host: &dyn SweepOps) -> Option<TopoShape> {
    let tables = host.tables();
    if host.has_shape(&tables.gen, tables.dir) {
        shape_of(host, &tables.gen.clone(), tables.dir)
    } else {
        None
    }
}

fn store(tables: &SweepTables, i: usize, j: usize, shape: TopoShape) {
    tables.shapes.borrow_mut()[i][j] = Some(shape);
}

fn cell(tables: &SweepTables, i: usize, j: usize) -> TopoShape {
    tables.shapes.borrow()[i][j]
        .clone()
        .unwrap_or_else(|| TopoShape::new(ShapeType::Shape))
}

fn mark_used(tables: &SweepTables, i: usize, j: usize) {
    if i > 0 && j > 0 {
        tables.used.borrow_mut()[i][j] = true;
    }
}

fn add_into_cell(
    builder: &TopoBuilder,
    tables: &SweepTables,
    i: usize,
    j: usize,
    child: &TopoShape,
    orien: Orientation,
) {
    let mut parent = cell(tables, i, j);
    add_oriented(builder, &mut parent, child, orien);
}

fn add_oriented(builder: &TopoBuilder, parent: &mut TopoShape, child: &TopoShape, orien: Orientation) {
    let mut c = child.clone();
    c.set_orientation(orien);
    builder.add(parent, &c);
}

/// `BRep_Tool::IsClosed` for a wire or a shell (`BRep_Tool.cxx:1707-1755`).
pub fn shape_is_closed(s: &TopoShape) -> bool {
    let mut has_bound = false;
    let mut odd: Vec<usize> = Vec::new();
    match s.shape_type() {
        ShapeType::Wire => {
            for v in crate::topo_tools_full::vertices_of(s) {
                if matches!(v.0.orientation(), Orientation::Internal | Orientation::External) {
                    continue;
                }
                has_bound = true;
                toggle(&mut odd, GeometryRegistry::shape_key(&v.0));
            }
        }
        ShapeType::Shell => {
            for e in crate::topo_tools_full::edges_of(s) {
                if crate::brep_tool::BRepTool::is_degenerated(&e)
                    || matches!(e.0.orientation(), Orientation::Internal | Orientation::External)
                {
                    continue;
                }
                has_bound = true;
                toggle(&mut odd, GeometryRegistry::shape_key(&e.0));
            }
        }
        _ => return s.closed(),
    }
    has_bound && odd.is_empty()
}

fn toggle(set: &mut Vec<usize>, k: usize) {
    match set.iter().position(|x| *x == k) {
        Some(i) => {
            set.remove(i);
        }
        None => set.push(k),
    }
}
