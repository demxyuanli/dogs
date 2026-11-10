use super::super::prelude::*;
use super::super::*;

impl BopBuilder {
    /// Empty constructor.
    pub fn new() -> Self {
        Self {
            filler: PaveFiller::new(),
            history: BopHistory::new(),
            origins: HashMap::new(),
            shapes_sd: HashMap::new(),
            sd_by_key: HashMap::new(),
            arguments: Vec::new(),
            objects: Vec::new(),
            tools: Vec::new(),
            obj_state: FaceState::Out,
            tools_state: FaceState::Out,
            errors: Vec::new(),
            warnings: Vec::new(),
            result_shape: TopoShape::new(ShapeType::Compound),
            in_parts: HashMap::new(),
            fill_history: true,
            check_open_solid: true,
        }
    }
}
