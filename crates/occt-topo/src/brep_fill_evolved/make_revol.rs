//! `BRepFill_Evolved::MakeRevol` (`BRepFill_Evolved.cxx:2329-2400`): revolve
//! the profile placed at the end of one spine edge around the vertical axis
//! through `VLast`, and record the generated shapes in `myMap`.

use std::f64::consts::PI;

use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPnt, GpTrsf};
use occt_core::toploc::TopLocLocation;

use crate::brep_sweep::Revol;
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::shape::{Edge, Vertex, Wire};

use super::evolved::BRepFillEvolved;
use super::put_profil_at::put_profil_at;
use super::shape_map::ShapeMap;
use super::side::{edge_vertices_oriented, side, FILL_CONFUSION};

impl BRepFillEvolved {
    /// `MakeRevol(SE, VLast, AxeRef)` (`BRepFill_Evolved.cxx:2329-2400`).
    pub fn make_revol(
        &mut self,
        se: &Edge,
        v_last: &Vertex,
        axe_ref: &GpAx3,
    ) -> Result<(), String> {
        let profile = self
            .my_profile
            .clone()
            .ok_or_else(|| "BRepFill_Evolved::MakeRevol: null profile".to_string())?;
        let spine = self
            .my_spine
            .clone()
            .ok_or_else(|| "BRepFill_Evolved::MakeRevol: null spine".to_string())?;

        // `cxx:2337-2341`: the same right-side turn as `MakePipe`.
        let mut trsf = GpTrsf::identity();
        if side(&profile, FILL_CONFUSION) > 3 {
            let oz = GpAx1::new(
                GpPnt::new(0.0, 0.0, 0.0),
                GpDir::new(0.0, 0.0, 1.0).map_err(|e| e.to_string())?,
            );
            trsf.set_rotation_ax1(&oz, PI).map_err(|e| e.to_string())?;
        }
        let dum_loc = TopLocLocation::composed(&trsf, profile.0.location().clone());
        let moved = profile.0.located(&dum_loc);

        // `cxx:2342-2343`: `PutProfilAt(..., false)`.
        let mut gen_prof = put_profil_at(&moved, axe_ref, &se.0, &spine.0, false)?;

        // `cxx:2349`: `gp_Ax1 AxeRev(BRep_Tool::Pnt(VLast), -gp::DZ())`.
        let axis = GpAx1::new(
            crate::brep_tool::BRepTool::vertex_point(v_last),
            GpDir::new(0.0, 0.0, -1.0).map_err(|e| e.to_string())?,
        );

        // `cxx:2353-2357`: turn the profile so the bissectrices do not cross
        // the sewings; `GenProf.Move(DummyLoc)`.
        let mut dummy = GpTrsf::identity();
        dummy
            .set_rotation_ax1(&axis, 1.5 * PI)
            .map_err(|e| e.to_string())?;
        let dummy_loc = TopLocLocation::composed(&dummy, TopLocLocation::identity());
        gen_prof.move_location(&dummy_loc);

        // `cxx:2359`: `BRepSweep_Revol Rev(GenProf, AxeRev, true)`.
        let rev = Revol::new(&gen_prof, &axis, true)?;

        // `cxx:2368-2369`.
        self.my_map.bind(&v_last.0, ShapeMap::new());

        // `cxx:2371-2399`: walk the reference and the placed profile in
        // parallel and transfer each generated edge / face.
        let mut prof_exp = WireExplorer::new();
        let mut gen_exp = WireExplorer::new();
        prof_exp.init(&profile, None);
        gen_exp.init(&Wire::wrap(gen_prof.clone()).ok_or_else(|| {
            "BRepFill_Evolved::MakeRevol: profile is not a wire".to_string()
        })?, None);
        let mut first_vertex = true;
        while prof_exp.more() && gen_exp.more() {
            let Some(prof_edge) = prof_exp.current().cloned() else {
                break;
            };
            let Some(gen_edge) = gen_exp.current().cloned() else {
                break;
            };

            let (vf, vl) = edge_vertices_oriented(&prof_edge);
            let (vfg, vlg) = edge_vertices_oriented(&gen_edge);

            // `cxx:2379-2388`.
            if first_vertex {
                if let Some(vf) = vf.as_ref() {
                    self.my_map
                        .entry_or_bind(&v_last.0, ShapeMap::new)
                        .bind(&vf.0, Vec::new());
                    if let Some(vfg) = vfg.as_ref() {
                        let rv = rev
                            .shape_of(&vfg.0)
                            .ok_or_else(|| "BRepSweep_Revol::Shape".to_string())?;
                        let oo = rv.orientation();
                        self.my_map
                            .list_mut(&v_last.0, &vf.0)
                            .push(rv.oriented(oo));
                    }
                }
                first_vertex = false;
            }

            // `cxx:2389-2392`.
            self.my_map
                .entry_or_bind(&v_last.0, ShapeMap::new)
                .bind(&prof_edge.0, Vec::new());
            let rf = rev
                .shape_of(&gen_edge.0)
                .ok_or_else(|| "BRepSweep_Revol::Shape".to_string())?;
            let oo = crate::abs::Orientation::compose(rf.orientation(), gen_edge.orientation());
            self.my_map
                .list_mut(&v_last.0, &prof_edge.0)
                .push(rf.oriented(oo));

            // `cxx:2393-2398`.
            if let (Some(vl), Some(vlg)) = (vl.as_ref(), vlg.as_ref()) {
                self.my_map
                    .entry_or_bind(&v_last.0, ShapeMap::new)
                    .bind(&vl.0, Vec::new());
                let rv = rev
                    .shape_of(&vlg.0)
                    .ok_or_else(|| "BRepSweep_Revol::Shape".to_string())?;
                let oo = rv.orientation();
                self.my_map.list_mut(&v_last.0, &vl.0).push(rv.oriented(oo));
            }

            prof_exp.next();
            gen_exp.next();
        }

        Ok(())
    }
}
