//! The "Nullify degenerated ChFi/Faces interferences" prelude of
//! `ChFi3d_FilDS`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_0.cxx:2611-2653`
//! (OCCT 8.0.0), with the two `eap occ293` / `eap occ354` comments of the
//! original kept at their sites.
//!
//! The pass has two jobs on a stripe with more than one `SurfData`:
//! 1. every face interference whose parameter range is degenerate
//!    (`|FirstParameter - LastParameter| <= Precision::PConfusion()`) loses its
//!    3D curve in the DS, because the curve is represented by the p-curve bound
//!    through the `SetSCI` of the surface interference,
//! 2. at the two ends of the stripe, the two common points of the degenerated
//!    pair swap their payloads: the arc is dropped and both take the second
//!    point's coordinates (`occ354`).

use std::sync::Arc;

use occt_core::precision::PCONFUSION;

use crate::chfi_ds::SurfData;
use crate::chfi3d::builder_0_helpers::fil_interf::reset_both_points;

use super::{write_vertex, ChFi3dBuilder};

impl ChFi3dBuilder {
    /// One `(j, onS)` cell of the nullify prelude (`cxx:2615-2652`).
    pub(super) fn nullify_degenerated_interference(
        &mut self,
        seq: &mut [Arc<SurfData>],
        j: usize,
        on_s: i32,
    ) -> Result<(), String> {
        // `const ChFiDS_FaceInterference& Fi = Fd->Interference(onS); IcFil1 =
        // Fi.LineIndex();` (`cxx:2621-2622`).
        let fd = super::seq_at(seq, j)?;
        let fi = fd.interference(on_s).clone();
        let ic_fil = fi.line_index();
        // `if (!IcFil1) continue;` (`cxx:2623-2626`).
        if ic_fil == 0 {
            return Ok(());
        }
        // `double FiLen = std::abs(Fi.FirstParameter() -
        // Fi.LastParameter()); if (FiLen > Precision::PConfusion()) continue;`
        // (`cxx:2627-2631`).
        let fi_len = (fi.first_parameter() - fi.last_parameter()).abs();
        if fi_len > PCONFUSION {
            return Ok(());
        }
        // `TopOpeBRepDS_Curve& cc = DStr.ChangeCurve(IcFil1);
        // cc.ChangeCurve().Nullify();` (`cxx:2632-2633`).
        *self.my_ds.change_ds().change_curve(ic_fil).change_curve() = None;

        // `if (j != 1 && j != SeqFil.Length()) continue;` (`cxx:2636-2639`).
        let len = seq.len();
        if j != 1 && j != len {
            return Ok(());
        }

        // `bool isfirst = (j == 1); int i = isfirst ? j + 1 : j - 1;
        // ChFiDS_CommonPoint& CP1 = SeqFil(i)->ChangeVertex(isfirst, onS);`
        // (`cxx:2640-2642`).
        let is_first = j == 1;
        let i = if is_first { j + 1 } else { j - 1 };

        // `if (Fd->Vertex(isfirst, onS).IsOnArc() && CP1.IsOnArc())`
        // (`cxx:2643`).
        let fd_on_arc = fd.vertex(is_first, on_s).is_on_arc();
        let cp1_on_arc = seq
            .get(i - 1)
            .ok_or_else(|| format!("ChFi3d_FilDS: SeqFil has no element {i}"))?
            .vertex(is_first, on_s)
            .is_on_arc();
        if !(fd_on_arc && cp1_on_arc) {
            return Ok(());
        }

        // `ChFiDS_CommonPoint& CP2 = Fd->ChangeVertex(!isfirst, onS);
        // CP1.Reset(); CP1.SetPoint(CP2.Point());
        // CP2.Reset(); CP2.SetPoint(CP1.Point());` (`cxx:2645-2649`).
        let mut cp1 = seq[i - 1].vertex(is_first, on_s).clone();
        let mut cp2 = fd.vertex(!is_first, on_s).clone();
        reset_both_points(&mut cp1, &mut cp2);
        write_vertex(seq, i, is_first, on_s, cp1);
        write_vertex(seq, j, !is_first, on_s, cp2);
        Ok(())
    }
}
