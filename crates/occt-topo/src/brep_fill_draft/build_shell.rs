//! `BRepFill_Draft::BuildShell` (`BRepFill_Draft.cxx:467-531`): the sweep plus
//! its orientation control and the optional stopping face.

use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;
use std::sync::Arc;

use occt_core::gp::GpDir;
use occt_core::precision::Precision;
use occt_geom::geom_fill::sweep::ApproxStyle;
use occt_geom::Surface;

use crate::abs::{Orientation, ShapeType};
use crate::brep_fill_location_law::LocationLaw;
use crate::brep_fill_sweep::helpers::ShapeTable;
use crate::brep_fill_sweep::BrepFillSweep;
use crate::brep_tool::BRepTool;
use crate::shape::Shell;
use crate::topo_tools::{explorer_shapes, ShapeMapKey};

use super::BrepFillDraft;

impl BrepFillDraft {
    /// `BuildShell(Surf, KeepOutSide)` (`cxx:467-531`).
    ///
    /// When the sweep fails the reference clears `myDone` and returns before
    /// the stopping face is added, which is what this does.
    pub(super) fn build_shell(
        &mut self,
        surf: Option<&Arc<dyn Surface>>,
        keep_out_side: bool,
    ) -> Result<(), String> {
        // `cxx:469-471`.
        let sec = self
            .my_sec
            .clone()
            .ok_or_else(|| "Standard_NullObject: BRepFill_Draft::mySec".to_string())?;
        let law = self
            .my_loc
            .clone()
            .ok_or_else(|| "Standard_NullObject: BRepFill_Draft::myLoc".to_string())?;
        let loc: Arc<dyn LocationLaw> = law;
        let mut sweep = BrepFillSweep::new(sec, loc, true);
        // `cxx:470`: `Sweep.SetTolerance(myTol)` replaces `Tol3d` only.
        sweep.set_tolerance_3d(self.my_tol);
        sweep.set_angular_control(self.angmin, self.angmax);

        // `cxx:473-481`: the three dummy containers of `Build`.
        let mut dummy: HashSet<ShapeMapKey> = HashSet::new();
        let mut dummy2: HashMap<ShapeMapKey, ShapeTable> = HashMap::new();
        let mut dummy3: HashMap<ShapeMapKey, ShapeTable> = HashMap::new();
        let _ = sweep.build(
            &mut dummy,
            &mut dummy2,
            &mut dummy3,
            self.my_style,
            self.my_cont,
            ApproxStyle::Location,
            11,
            30,
        );

        // `cxx:482-521`.
        if !sweep.is_done() {
            self.my_done = false;
            return Ok(());
        }
        self.my_shape = sweep.shape().cloned();
        self.my_shell = self.my_shape.clone().and_then(Shell::wrap);
        self.my_faces = sweep.sub_shape().clone();
        self.my_sections = sweep.sections().clone();
        self.my_done = true;

        // Control of the orientation (`cxx:490-514`).
        let mut out = true;
        let shell = self.my_shell.clone();
        if let Some(shell) = shell {
            if let Some(face) = explorer_shapes(&shell.0, ShapeType::Face).first().cloned() {
                let face = crate::shape::Face::wrap(face);
                if let Some(face) = face {
                    if let Some(surface) = BRepTool::face_surface_world(&face) {
                        // `BRepAdaptor_Surface SF(F)` is built with
                        // `Restriction = true` (`BRepAdaptor_Surface.cxx:56-81`):
                        // `Load(aSurface, umin, umax, vmin, vmax,
                        // L.Transformation())` with `BRepTools::UVBounds(F)`,
                        // so `FirstUParameter` / `FirstVParameter` answer the
                        // face's UV bounds, not the surface's own.
                        let (u, _, v, _) = crate::brep_tools::uv_bounds(&face);
                        let (_p, v1, v2) = surface.d1(u, v);
                        let mut normal = v1.crossed(&v2);
                        if face.0.orientation() == Orientation::Reversed {
                            normal.reverse();
                        }
                        if normal.magnitude() > 1.0e-10 {
                            // `cxx:509`: `out = myDir.Angle(V) > M_PI / 2;`
                            let dir = GpDir::from_vec(&normal)
                                .map_err(|e| format!("gp_Dir: {e}"))?;
                            out = self.my_dir.angle(&dir) > PI / 2.0;
                        }
                    }
                }
            }
        }
        if out == self.is_internal {
            if let Some(shell) = self.my_shell.as_mut() {
                shell.0.reverse();
            }
            if let Some(shape) = self.my_shape.as_mut() {
                shape.reverse();
            }
        }

        // `cxx:523-530`: "Add the face at end".
        if let Some(surf) = surf {
            let face = crate::brep_lib_make_face::make_face_from_surface(
                surf.clone(),
                Precision::CONFUSION,
            );
            self.fuse(&face.0, keep_out_side)?;
        }
        Ok(())
    }
}
