

use occt_core::gp::{
    GpAx1, GpDir, GpTrsf,
};

use occt_topo::abs::ShapeType;
use occt_topo::shape::TopoShape;


pub(super) fn collect_edges(shape: &TopoShape, out: &mut Vec<TopoShape>) {
    if shape.shape_type() == ShapeType::Edge {
        out.push(shape.clone());
        return;
    }
    let children = shape.tshape.read().expect("lock").children.clone();
    for child in children {
        collect_edges(&child, out);
    }
}

pub(super) fn moved_revolution_axis(axis: &GpAx1, tr: &GpTrsf) -> Result<GpAx1, String> {
    let loc = axis.location().transformed(tr);
    let mut d = *axis.direction().xyz();
    tr.transforms_xyz_dir(&mut d);
    let mut dir = GpDir::from_xyz(&d)
        .map_err(|_| "IGES 120 axis is null after location: gp_Dir::Transform".to_string())?;
    let det = tr.scale_factor() * tr.h_vectorial_part().determinant();
    if det < 0.0 {
        dir.reverse();
    }
    Ok(GpAx1::new(loc, dir))
}
