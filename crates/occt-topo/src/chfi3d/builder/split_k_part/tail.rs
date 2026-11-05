//! The second half of `ChFi3d_Builder::SplitKPart` (OCCT
//! `ChFi3d_Builder_SpKP.cxx:1412-1509`): the spine parameters of every
//! `SurfData`, then the extension / truncation of the two spine extremities.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_SpKP.cxx`
//! (OCCT 8.0.0).

use std::sync::Arc;

use crate::chfi_ds::{SpineHandle, SurfData};
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::super::ChFi3dBuilder;
use super::copy_data::remove_count;
use super::spine_params::{param_on_spine, reset_prol};

impl ChFi3dBuilder {
    /// `// Add parameters of the spine on SurfDatas.` through
    /// `// Truncation at the end.` (`cxx:1412-1509`).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn split_k_part_tail(
        &mut self,
        set_data: &mut Vec<Arc<SurfData>>,
        spine: &SpineHandle,
        iedge: i32,
        s1: &SurfaceAdaptor,
        s2: &SurfaceAdaptor,
        on_s: i32,
        cntl_fi_on_s: i32,
        intf: bool,
        intl: bool,
    ) -> Result<(), String> {
        let tolesp = self.tolesp;
        // Add parameters of the spine on SurfDatas.
        let (periodic, nb_edges, period) = {
            let guard = spine
                .read()
                .map_err(|_| "ChFi3d_Builder::SplitKPart: poisoned spine".to_string())?;
            (
                guard.is_periodic(),
                guard.nb_edges() as i32,
                guard.period().unwrap_or(0.),
            )
        };
        for i in 1..=set_data.len() {
            let cd7 = Arc::make_mut(&mut set_data[i - 1]);
            let ftg = cd7.interference(on_s).first_parameter();
            let ltg = cd7.interference(on_s).last_parameter();
            let (mut fsp, pokdeb) = param_on_spine(
                self.my_ds.change_ds(),
                ftg,
                cd7,
                spine,
                iedge,
                intf,
                intl,
                tolesp,
            )?;
            if !pokdeb {
                fsp = reset_prol(self.my_ds.change_ds(), cd7, spine, iedge, true)?;
            }
            let (mut lsp, pokfin) = param_on_spine(
                self.my_ds.change_ds(),
                ltg,
                cd7,
                spine,
                iedge,
                intf,
                intl,
                tolesp,
            )?;
            if !pokfin {
                lsp = reset_prol(self.my_ds.change_ds(), cd7, spine, iedge, false)?;
            }
            if periodic && iedge == nb_edges && lsp < fsp {
                lsp += period;
            } else if periodic && iedge == 1 && lsp < fsp {
                fsp -= period;
            }
            cd7.set_first_spine_param(fsp);
            cd7.set_last_spine_param(lsp);
        }

        if intf && !set_data.is_empty() {
            // Extension of the spine.
            let fsp = set_data.first().expect("SplitKPart").first_spine_param();
            spine
                .write()
                .map_err(|_| "ChFi3d_Builder::SplitKPart: poisoned spine".to_string())?
                .set_first_parameter(fsp);
        } else {
            // Truncation at the beginning.
            let mut i = 1usize;
            while i <= set_data.len() {
                let f_spine = spine
                    .read()
                    .map_err(|_| {
                        "ChFi3d_Builder::SplitKPart: poisoned spine".to_string()
                    })?
                    .first_parameter_of(iedge)
                    .ok_or_else(|| {
                        "ChFi3d_Builder::SplitKPart: FirstParameter".to_string()
                    })?;
                let (fsp, lsp) = {
                    let cd8 = &set_data[i - 1];
                    (cd8.first_spine_param(), cd8.last_spine_param())
                };
                if lsp > f_spine {
                    if fsp > f_spine {
                        break;
                    }
                    let cd8 = Arc::make_mut(&mut set_data[i - 1]);
                    self.trunc(cd8, spine, s1, s2, iedge, true, cntl_fi_on_s)?;
                    break;
                }
                i += 1;
            }
            if i > 1 {
                remove_count(set_data, 1, i - 1)?;
            }
        }

        if intl && !set_data.is_empty() {
            // Extension of the spine.
            let lsp = set_data.last().expect("SplitKPart").last_spine_param();
            spine
                .write()
                .map_err(|_| "ChFi3d_Builder::SplitKPart: poisoned spine".to_string())?
                .set_last_parameter(lsp);
        } else {
            // Truncation at the end.
            let mut i = set_data.len();
            loop {
                if i == 0 {
                    break;
                }
                let l_spine = spine
                    .read()
                    .map_err(|_| {
                        "ChFi3d_Builder::SplitKPart: poisoned spine".to_string()
                    })?
                    .last_parameter_of(iedge)
                    .ok_or_else(|| {
                        "ChFi3d_Builder::SplitKPart: LastParameter".to_string()
                    })?;
                let (fsp, lsp) = {
                    let cd9 = &set_data[i - 1];
                    (cd9.first_spine_param(), cd9.last_spine_param())
                };
                if fsp < l_spine {
                    if lsp < l_spine {
                        break;
                    }
                    let cd9 = Arc::make_mut(&mut set_data[i - 1]);
                    self.trunc(cd9, spine, s1, s2, iedge, false, cntl_fi_on_s)?;
                    break;
                }
                i -= 1;
            }
            if i < set_data.len() {
                remove_count(set_data, i + 1, set_data.len())?;
            }
        }
        Ok(())
    }
}
