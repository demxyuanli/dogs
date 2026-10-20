//! TEMPORARY diagnostic probe (delete after use).
use occt_topo::meshing::range_splitter::MeshSurfaceType;
use occt_topo::model::BRepModel;
use occt_topo::shape::TopoShape;
use occt_topo::step::read_step_file;

#[path = "model_probe/diagnostics_census.rs"]
mod diagnostics_census;
#[path = "model_probe/diagnostics_faces.rs"]
mod diagnostics_faces;
#[path = "model_probe/diagnostics_mesh.rs"]
mod diagnostics_mesh;
#[path = "model_probe/diagnostics_topology.rs"]
mod diagnostics_topology;
#[path = "model_probe/report.rs"]
mod report;

pub(crate) struct Ctx {
    pub(crate) model: BRepModel,
    pub(crate) shape: TopoShape,
    pub(crate) override_lin: Option<f64>,
    pub(crate) mode_edges: bool,
    pub(crate) mode_dev: bool,
    pub(crate) mode_low: bool,
    pub(crate) mode_fixms: bool,
    pub(crate) mode_fstats: bool,
}

fn tname(t: MeshSurfaceType) -> &'static str {
    match t {
        MeshSurfaceType::Plane => "Plane", MeshSurfaceType::Cylinder => "Cylinder",
        MeshSurfaceType::Cone => "Cone", MeshSurfaceType::Sphere => "Sphere",
        MeshSurfaceType::Torus => "Torus", MeshSurfaceType::SurfaceOfRevolution => "Rev",
        MeshSurfaceType::SurfaceOfExtrusion => "Extr", MeshSurfaceType::BezierSurface => "Bezier",
        MeshSurfaceType::BSplineSurface => "BSpline", MeshSurfaceType::OffsetSurface => "Offset",
        MeshSurfaceType::OtherSurface => "Other",
    }
}

fn main() {
    let path = std::env::args().nth(1).expect("usage: model_probe <step> [lin_override|--edges]");
    let override_lin: Option<f64> = std::env::args().nth(2).and_then(|s| s.parse().ok());
    let model = read_step_file(&path).expect("read step");
    // Match the OCCT probe's `aReader.OneShape()`: the FIRST transfer root.
    // `model.shapes` can hold both a wrapper `Compound` and the `Solid` inside
    // it, so compounding every root would traverse the solid 2-3 times.
    let shape = model.shapes[0].shape.clone();
    let ctx = Ctx {
        model,
        shape,
        override_lin,
        mode_edges: std::env::args().any(|a| a == "--edges"),
        mode_dev: std::env::args().any(|a| a == "--dev"),
        mode_low: std::env::args().any(|a| a == "--low"),
        mode_fixms: std::env::args().any(|a| a == "--fixms"),
        mode_fstats: std::env::args().any(|a| a == "--fstats"),
    };
    if diagnostics_faces::run(&ctx) {
        return;
    }
    if diagnostics_topology::run(&ctx) {
        return;
    }
    if diagnostics_census::run(&ctx) {
        return;
    }
    if diagnostics_mesh::run(&ctx) {
        return;
    }
    report::run(&ctx);
}
