//! `BRepOffset_MakeOffset::MakeSolid` (`BRepOffset_MakeOffset.cxx:3794-3894`)
//! together with its two file-local helpers `UpdateInitOffset`
//! (`:3110-3147`) and `CorrectSolid` (`:4304-4382`).
//!
//! `MakeSolid` first rebinds the history images to the shells of
//! `myOffsetShape` (`:3805-3807`), then walks every shell (`:3818-3845`):
//! closed shells go into the solid when the input was a solid or the mode is
//! thickening (`:3816`, `:3832`), the rest go into a compound, and `S1` keeps
//! the first of them (`:3838-3842`). A solid of more than one shell is checked
//! with `BRepCheck_Analyzer`; an invalid one is repaired by `CorrectSolid`
//! (`:3847-3867`). The final choice between the solid, the compound and `S1`
//! is `:3868-3893`.
//!
//! UNPORTED / port notes:
//!
//! * `BRepCheck_Analyzer` is [`crate::brep_check::Analyzer`]: `new` records a
//!   result per sub-shape (`brep_check/mod.rs:274-307`) and `is_valid` reads
//!   them (`:311-313`). What the recorded flags cover is the per-type checks of
//!   `brep_check` (`vertex_check`, `edge_check`, ...), not the full
//!   `BRepCheck_*` set; a repair is attempted whenever the recorded results
//!   report an error. `CorrectSolid` is ported so the `nbs > 1` arm
//!   (`:3847-3867`) matches the reference line for line.
//! * `BRepGProp::VolumeProperties(aSh, aVProps, true)` (`:4315`) filters to
//!   closed shells; the port's `brep_gprop_full::properties::volume_properties`
//!   has no `OnlyClosed` switch and integrates every face of the shell.
//! * `Message_ProgressScope` / `!More()` -> `BRepOffset_UserBreak`
//!   (`:3821-3825`).
//! * `BRepAlgo_Image::Filter` is `brep_algo_image::BRepAlgoImage::filter`
//!   (`cxx:3141-3143`).

use occt_core::precision::CONFUSION;

use crate::abs::ShapeType;
use crate::brep_algo_image::BRepAlgoImage;
use crate::brep_check;
use crate::brep_class3d::SolidClassifier;
use crate::brep_gprop_full::volume_properties;
use crate::brep_tool::BRepTool;
use crate::builder::TopoBuilder;
use crate::shape::{Shell, Solid, TopoShape, Vertex};
use crate::topo_tools::{explorer_shapes, shape_map_key};

use super::BRepOffsetMakeOffset;

/// `UpdateInitOffset(myInitOffset, myImageOffset, myOffsetShape,
/// theShapeType)` (`cxx:3110-3147`).
fn update_init_offset(
    init_offset: &mut BRepAlgoImage,
    image_offset: &BRepAlgoImage,
    offset_shape: &TopoShape,
    shape_type: ShapeType,
) {
    let mut niof = BRepAlgoImage::new();
    let roots = init_offset.roots().to_vec();
    // `cxx:3119-3121`.
    for root in &roots {
        niof.set_root(root);
    }
    // `cxx:3122-3138`.
    for si in &roots {
        let mut l1: Vec<TopoShape> = Vec::new();
        init_offset.last_image(si, &mut l1);
        let mut li: Vec<TopoShape> = Vec::new();
        for o1 in &l1 {
            let mut l2: Vec<TopoShape> = Vec::new();
            image_offset.last_image(o1, &mut l2);
            li.extend(l2);
        }
        let _ = niof.bind_list(si, &li);
    }
    // `cxx:3140-3143`.
    let _ = niof.filter(offset_shape, shape_type);
    *init_offset = niof;
}

/// `CorrectSolid(theSol, theSolList)` (`cxx:4304-4382`).
///
/// Answers the corrected solid and the extra solids appended by the
/// `TopAbs_OUT` arm (`:4366-4374`).
fn correct_solid(sol: &Solid) -> (Solid, Vec<TopoShape>) {
    let builder = TopoBuilder::new();
    let shells: Vec<TopoShape> = explorer_shapes(&sol.0, ShapeType::Shell);

    let mut outer_shell: Option<TopoShape> = None;
    let mut vols: Vec<f64> = Vec::new();
    let mut vol_max = 0.0f64;
    let mut outer_vol = 0.0f64;
    // `cxx:4312-4327`.
    for a_sh in &shells {
        let mass = volume_properties(a_sh)
            .map(|p| p.mass())
            .unwrap_or(0.0);
        if mass.abs() > vol_max {
            outer_vol = mass;
            vol_max = outer_vol.abs();
            outer_shell = Some(a_sh.clone());
        }
        vols.push(mass);
    }
    let mut out: Vec<TopoShape> = Vec::new();
    // `cxx:4328-4331`.
    if outer_vol.abs() < CONFUSION {
        return (sol.clone(), out);
    }
    // `cxx:4332-4335`: a negative outer volume means the outer shell faces
    // inwards.
    if outer_vol < 0.0 {
        if let Some(s) = outer_shell.as_mut() {
            s.reverse();
        }
    }
    let Some(outer) = outer_shell.clone() else {
        return (sol.clone(), out);
    };

    // `cxx:4336-4338`.
    let mut new_sol = builder.make_solid(&[Shell(outer.clone())]);
    new_sol.0.set_closed(true);
    // `cxx:4339`.
    let mut classifier = SolidClassifier::new();
    classifier.load(new_sol.0.clone());

    // `cxx:4340-4377`.
    for (a_sh, vol) in shells.iter().zip(vols.iter()) {
        let mut a_sh = a_sh.clone();
        if shape_map_key(&a_sh) == shape_map_key(&outer) {
            continue;
        }
        let vertices = explorer_shapes(&a_sh, ShapeType::Vertex);
        let Some(v_shape) = vertices.first().cloned() else {
            continue;
        };
        let Some(v) = Vertex::wrap(v_shape) else {
            continue;
        };
        let p = BRepTool::vertex_point(&v);
        classifier.perform(&p, BRepTool::vertex_tolerance(&v));
        if classifier.state() == crate::fclass2d::FaceState::In {
            if *vol > 0.0 {
                a_sh.reverse();
            }
            let mut shells_now: Vec<Shell> = explorer_shapes(&new_sol.0, ShapeType::Shell)
                .into_iter()
                .filter_map(Shell::wrap)
                .collect();
            shells_now.push(Shell(a_sh));
            new_sol = builder.make_solid(&shells_now);
            new_sol.0.set_closed(true);
        } else {
            if *vol < 0.0 {
                a_sh.reverse();
            }
            let a_sol = builder.make_solid(&[Shell(a_sh)]);
            a_sol.0.set_closed(true);
            out.push(a_sol.0);
        }
    }
    (new_sol, out)
}

impl BRepOffsetMakeOffset {
    /// `MakeSolid(theRange)` (`cxx:3794-3894`).
    pub fn make_solid(&mut self) {
        let Some(offset_shape) = self.my_offset_shape.clone() else {
            return;
        };

        // `cxx:3805-3807`: supporting history.
        update_init_offset(
            &mut self.my_init_offset_face,
            &self.my_image_offset,
            &offset_shape,
            ShapeType::Face,
        );
        update_init_offset(
            &mut self.my_init_offset_edge,
            &self.my_image_offset,
            &offset_shape,
            ShapeType::Edge,
        );

        let builder = TopoBuilder::new();
        let mut nb_shell = 0;
        let mut nc: Vec<TopoShape> = Vec::new();
        let mut s1: Option<TopoShape> = None;
        let mut sol_shells: Vec<Shell> = Vec::new();
        // `cxx:3816`: `aMakeSolid = (myShape.ShapeType() == TopAbs_SOLID) ||
        // myThickening`.
        let a_make_solid = self
            .my_shape
            .as_ref()
            .is_some_and(|s| s.shape_type() == ShapeType::Solid)
            || self.my_thickening;

        // `cxx:3818-3845`.
        for sh_shape in explorer_shapes(&offset_shape, ShapeType::Shell) {
            let Some(mut sh) = Shell::wrap(sh_shape) else {
                continue;
            };
            if self.my_thickening && self.my_offset > 0.0 {
                sh.0.reverse();
            }
            nb_shell += 1;
            if sh.0.is_closed() && a_make_solid {
                sol_shells.push(sh);
            } else {
                nc.push(sh.0.clone());
                if nb_shell == 1 {
                    s1 = Some(sh.0);
                }
            }
        }

        // `cxx:3846-3848`.
        let nbs = sol_shells.len();
        let mut sol = builder.make_solid(&sol_shells);
        sol.0.set_closed(true);
        let mut sol_is_null = nbs == 0;

        // `cxx:3849-3867`: checking solid.
        if nbs > 1 && !brep_check::is_valid(&sol.0) {
            let (new_sol, a_sol_list) = correct_solid(&sol);
            sol = new_sol;
            if !a_sol_list.is_empty() {
                nc.push(sol.0.clone());
                nc.extend(a_sol_list);
                sol_is_null = true;
            }
        }

        let nc_is_null = nc.is_empty();
        // `cxx:3868-3893`.
        if !sol_is_null && !nc_is_null {
            nc.push(sol.0);
            let comp = builder.make_compound_of(&nc);
            self.my_offset_shape = Some(comp.0);
        } else if sol_is_null && !nc_is_null {
            self.my_offset_shape = if nb_shell == 1 {
                s1
            } else {
                let comp = builder.make_compound_of(&nc);
                Some(comp.0)
            };
        } else if !sol_is_null && nc_is_null {
            self.my_offset_shape = Some(sol.0);
        } else {
            let comp = builder.make_compound_of(&nc);
            self.my_offset_shape = Some(comp.0);
        }
    }
}
