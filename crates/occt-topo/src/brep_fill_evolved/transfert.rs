//! `BRepFill_Evolved::Transfert` (`BRepFill_Evolved.cxx:1966-2053`): copies
//! the shape, the generated map and the top/bottom of the accumulated cut
//! volume into `self`, and repositions everything into the initial space.

use occt_core::toploc::TopLocLocation;

use crate::brep_lib_same_parameter::same_parameter_edge_inplace;
use crate::brep_tool::BRepTool;
use crate::shape::TopoShape;
use crate::tgeometry::GeometryRegistry;
use crate::topo_tools::edges_of;

use super::evolved::BRepFillEvolved;
use super::shape_map::ShapeMap;

impl BRepFillEvolved {
    /// `Transfert(Vevo, MapProf, MapSpine, LS, InitLS, InitLP)`
    /// (`cxx:1966-2053`).
    pub fn transfert(
        &mut self,
        vevo: &BRepFillEvolved,
        map_prof: &ShapeMap<TopoShape>,
        map_spine: &ShapeMap<TopoShape>,
        ls: &TopLocLocation,
        init_ls: &TopLocLocation,
        init_lp: &TopLocLocation,
    ) {
        // `cxx:1976-1980`.
        self.my_shape = vevo.shape().cloned();
        if let Some(spine) = self.my_spine.as_mut() {
            spine.0.set_location(init_ls);
        }
        if let Some(profile) = self.my_profile.as_mut() {
            profile.0.set_location(init_lp);
        }
        if let Some(shape) = self.my_shape.as_mut() {
            shape.move_location(ls);
        }

        // `cxx:1987-1996`: the SameRange / SameParameter flags are cleared and
        // `BRepLib::SameParameter(Edge)` recomputes them.
        if let Some(shape) = self.my_shape.clone() {
            let registry = GeometryRegistry::global();
            for e in edges_of(&shape) {
                registry.set_same_range(&e.0, false);
                registry.set_same_parameter(&e.0, false);
                let tol = BRepTool::edge_tolerance(&e);
                same_parameter_edge_inplace(&e, tol);
            }
        }

        // `cxx:2012-2045`: transfer `Vevo`'s map, repositioning both keys.
        let spine_keys: Vec<TopoShape> = vevo.my_map.keys().cloned().collect();
        for key in spine_keys {
            // OCCT reads `MapSpine(iteS.Key())`, which raises when the work
            // spine shape has no image; `PrepareSpine` binds every one.
            let Some(initial_spine) = map_spine.find(&key) else {
                continue;
            };
            let mut initial_spine = initial_spine.clone();
            initial_spine.move_location(ls);

            let Some(inner) = vevo.my_map.find(&key) else {
                continue;
            };
            let prof_keys: Vec<TopoShape> = inner.keys().cloned().collect();
            for prof_key in prof_keys {
                let Some(initial_prof) = map_prof.find(&prof_key) else {
                    continue;
                };
                let mut initial_prof = initial_prof.clone();
                initial_prof.set_location(init_lp);

                let Some(gen_shapes) = vevo.my_map.list(&key, &prof_key) else {
                    continue;
                };
                let target = self.my_map.list_mut(&initial_spine, &initial_prof);
                for s in gen_shapes.clone() {
                    let mut moved = s;
                    moved.move_location(ls);
                    target.push(moved);
                }
            }
        }

        // `cxx:2048-2052`.
        self.my_top = vevo.top().cloned().map(|mut s| {
            s.move_location(ls);
            s
        });
        self.my_bottom = vevo.bottom().cloned().map(|mut s| {
            s.move_location(ls);
            s
        });
    }
}
