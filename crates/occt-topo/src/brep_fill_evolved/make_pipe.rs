//! `BRepFill_Evolved::MakePipe` (`BRepFill_Evolved.cxx:2266-2327`): sweep the
//! profile along one edge of the spine with `BRepFill_Pipe`, and record the
//! generated faces / edges in `myMap`.

use std::f64::consts::PI;

use occt_core::gp::{GpAx1, GpAx3, GpDir, GpPnt, GpTrsf};
use occt_core::toploc::TopLocLocation;
use occt_geom::geom_fill::GeomFillTrihedron;

use crate::brep_fill_pipe::Pipe;
use crate::brep_lib_make_wire::MakeWire;
use crate::brep_tools_modifier::{Modifier, TrsfModification};
use crate::brep_tools_wire_explorer::WireExplorer;
use crate::shape::{Edge, Wire};

use super::evolved::BRepFillEvolved;
use super::put_profil_at::put_profil_at;
use super::shape_map::ShapeMap;
use super::side::{edge_vertices_oriented, side, FILL_CONFUSION};

impl BRepFillEvolved {
    /// `MakePipe(SE, AxeRef)` (`BRepFill_Evolved.cxx:2266-2327`).
    pub fn make_pipe(&mut self, se: &Edge, axe_ref: &GpAx3) -> Result<(), String> {
        let profile = self
            .my_profile
            .clone()
            .ok_or_else(|| "BRepFill_Evolved::MakePipe: null profile".to_string())?;
        let spine = self
            .my_spine
            .clone()
            .ok_or_else(|| "BRepFill_Evolved::MakePipe: null spine".to_string())?;

        // `cxx:2272-2276`: a right-side profile is turned by `Pi` around OZ
        // before it is placed.
        let mut trsf = GpTrsf::identity();
        if side(&profile, FILL_CONFUSION) > 3 {
            let oz = GpAx1::new(
                GpPnt::new(0.0, 0.0, 0.0),
                GpDir::new(0.0, 0.0, 1.0).map_err(|e| e.to_string())?,
            );
            trsf
                .set_rotation_ax1(&oz, PI)
                .map_err(|e| e.to_string())?;        }
        let dum_loc = TopLocLocation::composed(&trsf, profile.0.location().clone());
        let moved = profile.0.located(&dum_loc);

        // `cxx:2277-2278`.
        let dummy_prof = put_profil_at(&moved, axe_ref, &se.0, &spine.0, true)?;

        // `cxx:2285-2290`: rebuild a copy of the profile so the locations do
        // not accumulate on the edges of `myProfile`.
        let trsf_mod = TrsfModification::new(GpTrsf::identity());
        let modif = Modifier::on_shape_with(&dummy_prof, &trsf_mod)?;
        let gen_shape = modif.modified_shape(&dummy_prof)?;
        let gen_prof = Wire::wrap(gen_shape)
            .ok_or_else(|| "BRepFill_Evolved::MakePipe: profile is not a wire".to_string())?;

        // `cxx:2292`: `BRepFill_Pipe Pipe(BRepLib_MakeWire(SE), GenProf)`.
        let mut make_wire = MakeWire::new();
        let _ = make_wire.add(se);
        let spine_wire = make_wire.wire();
        let pipe = Pipe::with_spine(
            &spine_wire,
            &gen_prof.0,
            GeomFillTrihedron::CorrectedFrenet,
            false,
            false,
        )?;

        // `cxx:2305-2306`: `myMap.Bind(SE, P)`.
        self.my_map.bind(&se.0, ShapeMap::new());

        // `cxx:2308-2326`: walk the reference and the generated profile in
        // parallel.
        let mut prof_exp = WireExplorer::new();
        let mut gen_exp = WireExplorer::new();
        prof_exp.init(&profile, None);
        gen_exp.init(&gen_prof, None);
        let mut first_vertex = true;
        while prof_exp.more() && gen_exp.more() {
            let Some(prof_edge) = prof_exp.current().cloned() else {
                break;
            };
            let Some(gen_edge) = gen_exp.current().cloned() else {
                break;
            };

            // `cxx:2311-2312`.
            let (vf, vl) = edge_vertices_oriented(&prof_edge);
            let (vfg, vlg) = edge_vertices_oriented(&gen_edge);

            // `cxx:2314-2320`.
            if first_vertex {
                if let (Some(vf), Some(vfg)) = (vf.as_ref(), vfg.as_ref()) {
                    self.my_map
                        .entry_or_bind(&se.0, ShapeMap::new)
                        .bind(&vf.0, Vec::new());
                    let edge = pipe.edge(se, vfg)?;
                    self.my_map
                        .list_mut(&se.0, &vf.0)
                        .push(edge.0.clone());
                }
                first_vertex = false;
            }

            // `cxx:2321-2325`.
            if let (Some(vl), Some(vlg)) = (vl.as_ref(), vlg.as_ref()) {
                self.my_map
                    .entry_or_bind(&se.0, ShapeMap::new)
                    .bind(&vl.0, Vec::new());
                let edge = pipe.edge(se, vlg)?;
                self.my_map.list_mut(&se.0, &vl.0).push(edge.0.clone());
            }
            self.my_map
                .entry_or_bind(&se.0, ShapeMap::new)
                .bind(&prof_edge.0, Vec::new());
            let face = pipe.face(se, &gen_edge)?;
            self.my_map
                .list_mut(&se.0, &prof_edge.0)
                .push(face.0.clone());

            prof_exp.next();
            gen_exp.next();
        }

        Ok(())
    }
}
