//! `TopOpeBRepBuild_GridSS.cxx` — the solid/shell half of the
//! solid/solid grid.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_GridSS.cxx`
//! (OCCT 8.0.0):
//!
//! | OCCT | lines |
//! |---|---|
//! | `GMergeSolids(LSO1, LSO2, G1)` | `282-350` |
//! | `GFillSolidsSFS(LS1, LS2, G1, SFS)` | `353-428` |
//! | `GFillSolidSFS(SO1, LSO2, G1, SFS)` | `433-494` |
//! | `GFillSurfaceTopologySFS(SO1, G1, SFS)` | `500-528` |
//! | `GFillSurfaceTopologySFS(SSit, G1, SFS)` | `532-590` |
//! | `GFillShellSFS(SH, LSO2, G1, SFS)` | `595-645` |
//!
//! `GFillSurfaceTopologySFS(SO1, G1, SFS)` is `NYI` in OCCT (`cxx:525`): the
//! body only prints. The port keeps the `StatesON`/`Type` reads that can throw
//! and nothing else.
//!
//! UNPORTED: the `OCCT_DEBUG` blocks (`cxx:295-304`, `:382-395`, `:406-419`,
//! `:444-451`, `:474-483`, `:511-527`, `:551-564`, `:574-586`, `:604-612`,
//! `:627-633`) — trace and dump hooks only.

use crate::abs::{Orientation, ShapeType, State};
use crate::shape::TopoShape;
use crate::topo_tools::explorer_shapes;
use crate::topopebrep::build::g_topo::GTopo;
use crate::topopebrep::build::shell_face_set::ShellFaceSet;
use crate::topopebrep::build::solid_builder::SolidBuilder;
use crate::topopebrep::ds::surface_iterator::SurfaceIterator;

use super::Builder;

impl Builder {
    /// `GMergeSolids(LSO1, LSO2, G1)` (`cxx:282-350`).
    pub fn g_merge_solids(&mut self, lso1: &[TopoShape], lso2: &[TopoShape], g1: &GTopo) {
        let Some(so1) = lso1.first().cloned() else {
            return;
        };
        let (tb1, tb2) = g1.states_on();

        self.my_solid_reference = Some(so1.clone());

        let mut sfs = ShellFaceSet::with_shape(&so1);
        self.g_fill_solids_sfs(lso1, lso2, g1, &mut sfs);

        // `cxx:311-319`: build the new solids on a `FORWARD` copy of `SO1`.
        let mut so1f = lso1[0].clone();
        so1f.set_orientation(Orientation::Forward);
        let mut sobu = SolidBuilder::with_shell_face_set(&mut sfs, true);

        let mut lsom: Vec<TopoShape> = Vec::new();
        self.gsobu_make_solids(&so1f, &mut sobu, &mut lsom);

        // `cxx:322-332`: record the new solids on the shape-1 solids.
        for a_so1 in lso1 {
            if self.is_merged(a_so1, tb1) {
                continue;
            }
            let sol = self.change_merged(a_so1, tb1);
            sol.clear();
            sol.extend(lsom.iter().cloned());
        }

        // `cxx:335-348`: and on the shape-2 solids.
        for so2 in lso2 {
            if self.is_merged(so2, tb2) {
                continue;
            }
            let sol = self.change_merged(so2, tb2);
            sol.clear();
            sol.extend(lsom.iter().cloned());
        }
    }

    /// `GFillSolidsSFS(LS1, LS2, G1, SFS)` (`cxx:353-428`).
    pub fn g_fill_solids_sfs(
        &mut self,
        ls1: &[TopoShape],
        ls2: &[TopoShape],
        g1: &GTopo,
        sfs: &mut ShellFaceSet,
    ) {
        if ls1.is_empty() {
            return;
        }
        let (tb1, tb2) = g1.states_on();
        self.my_process_on = self.opecom() || self.opefus();
        if self.my_process_on {
            self.my_on_faces_map.clear();
        }

        self.my_solid_reference = Some(ls1[0].clone());

        // `cxx:378-398`: shape 1 solids against `LS2`.
        let mut g = g1.clone();
        for s in ls1 {
            if !self.is_merged(s, tb1) {
                self.g_fill_solid_sfs(s, ls2, &g, sfs);
            }
        }

        // `cxx:400-421`: shape 2 solids against `LS1`, permuted grid.
        g = g1.copy_permuted();
        for s in ls2 {
            if !self.is_merged(s, tb2) {
                self.g_fill_solid_sfs(s, ls1, &g, sfs);
            }
        }

        if self.my_process_on {
            self.add_on_patches_sfs(g1, sfs);
            self.my_process_on = false;
        }
    }

    /// `GFillSolidSFS(SO1, LSO2, G1, SFS)` (`cxx:433-494`).
    ///
    /// Virtual in OCCT (`Builder.hxx`), overridden by
    /// `Builder1::GFillSolidSFS` (`Builder1.cxx:186-241`); the branch on
    /// [`Builder::is_builder1`] is where the vtable would land.
    pub fn g_fill_solid_sfs(
        &mut self,
        so1: &TopoShape,
        lso2: &[TopoShape],
        g1: &GTopo,
        sfs: &mut ShellFaceSet,
    ) {
        if self.is_builder1 {
            self.b1_g_fill_solid_sfs(so1, lso2, g1, sfs);
            return;
        }
        self.base_g_fill_solid_sfs(so1, lso2, g1, sfs);
    }

    /// `TopOpeBRepBuild_Builder::GFillSolidSFS` (`cxx:433-494`): the classic
    /// builder body, also what `Builder1::GFillSolidSFS` falls back to when
    /// `GLOBAL_USE_NEW_BUILDER` is false (`Builder1.cxx:192-196`).
    pub(crate) fn base_g_fill_solid_sfs(
        &mut self,
        so1: &TopoShape,
        lso2: &[TopoShape],
        g1: &GTopo,
        sfs: &mut ShellFaceSet,
    ) {
        let (tb1, _tb2) = g1.states_on();
        let rev_ori1 = g1.is_to_reverse1();

        // `cxx:453-457`: work on a `FORWARD` solid.
        let mut sof = so1.clone();
        sof.set_orientation(Orientation::Forward);
        self.my_solid_to_fill = Some(sof.clone());

        for mut sh in explorer_shapes(&sof, ShapeType::Shell) {
            let hasshape = self
                .data_structure()
                .map(|hds| hds.has_shape(&sh, true))
                .unwrap_or(false);

            if !hasshape {
                // `cxx:466-488`: the shell is not in the DS, classify it
                // against the `LSO2` solids.
                let keep = self.g_keep_shape(&sh, lso2, tb1);
                if keep {
                    let ori_sh = sh.orientation();
                    let newori_sh = Builder::orient(ori_sh, rev_ori1);
                    sh.set_orientation(newori_sh);
                    sfs.add_shape(&sh);
                }
            } else {
                // `cxx:489-492`: the shell has faces with geometry.
                self.g_fill_shell_sfs(&sh, lso2, g1, sfs);
            }
        }
    }

    /// `GFillSurfaceTopologySFS(SO1, G1, SFS)` (`cxx:500-528`): `NYI` in OCCT,
    /// the body only reads the grid and prints.
    pub fn g_fill_surface_topology_sfs(&self, g1: &GTopo) {
        let _ = g1.states_on();
        let _ = g1.type_();
    }

    /// `GFillSurfaceTopologySFS(SSit, G1, SFS)` (`cxx:532-590`).
    ///
    /// `SSit` is OCCT's `TopOpeBRepDS_SurfaceIterator`; the port passes the
    /// iterator by reference because it is borrowed from the data structure.
    pub fn g_fill_surface_topology_sfs_it(
        &self,
        ss_it: &SurfaceIterator,
        g1: &GTopo,
        sfs: &mut ShellFaceSet,
    ) {
        let (tb1, _tb2) = g1.states_on();
        let conf = g1.config1();
        let mut tb = tb1;
        if conf == crate::topopebrep::ds::kind::Config::DiffOriented {
            // `cxx:540-550`.
            if tb1 == State::Out {
                tb = State::In;
            } else if tb1 == State::In {
                tb = State::Out;
            }
        }

        // `cxx:567-568`: `iG` is the index of the new surface.
        let i_g = ss_it.geometry();
        let lnew_f = self.new_faces(i_g).map(|v| v.to_vec()).unwrap_or_default();
        for mut f in lnew_f {
            let ori = ss_it.orientation(tb);
            f.set_orientation(ori);
            sfs.add_element(&f);
        }
    }

    /// `GFillShellSFS(SH, LSO2, G1, SFS)` (`cxx:595-645`).
    #[allow(clippy::too_many_arguments)]
    pub fn g_fill_shell_sfs(
        &mut self,
        sh: &TopoShape,
        lso2: &[TopoShape],
        g1: &GTopo,
        sfs: &mut ShellFaceSet,
    ) {
        // Virtual in OCCT, overridden by `Builder1::GFillShellSFS`
        // (`Builder1.cxx:243-301`).
        if self.is_builder1 {
            self.b1_g_fill_shell_sfs(sh, lso2, g1, sfs);
            return;
        }
        self.base_g_fill_shell_sfs(sh, lso2, g1, sfs);
    }

    /// `TopOpeBRepBuild_Builder::GFillShellSFS` (`cxx:596-645`).
    pub(crate) fn base_g_fill_shell_sfs(
        &mut self,
        sh: &TopoShape,
        lso2: &[TopoShape],
        g1: &GTopo,
        sfs: &mut ShellFaceSet,
    ) {
        // `cxx:615-628`: every same-domain face first.
        for for_ in explorer_shapes(sh, ShapeType::Face) {
            let hsd = self
                .data_structure()
                .map(|hds| hds.has_same_domain(&for_, true))
                .unwrap_or(false);
            if hsd {
                self.g_fill_face_sfs(&for_, lso2, g1, sfs);
            }
        }

        // `cxx:635-645`: then the faces without same domain.
        for for_ in explorer_shapes(sh, ShapeType::Face) {
            let hsd = self
                .data_structure()
                .map(|hds| hds.has_same_domain(&for_, true))
                .unwrap_or(false);
            if !hsd {
                self.g_fill_face_sfs(&for_, lso2, g1, sfs);
            }
        }
    }
}
