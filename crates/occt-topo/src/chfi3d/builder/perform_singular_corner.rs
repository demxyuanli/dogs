//! `ChFi3d_Builder::PerformSingularCorner` - the degenerated-vertex corner.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder.cxx:682-757`
//! (OCCT 8.0.0).
//!
//! At a vertex whose stripes all carry a degenerated `SurfData` (both common
//! points of the end coincide), the port stores the vertex once in the DS and
//! gives every such stripe a single arc: `ChFi3d_ComputeArete` builds the 3D
//! curve and the p-curve, `TopOpeBRepDS_Curve(C3d, tolreached)` is added to the
//! DS, and the stripe records the curve index, the parameters, the p-curve and
//! the point index on both sides.
//!
//! The body is the same shape as `ChFi3d_SingularExtremity`
//! (`ChFi3d_Builder_2.cxx:1635-1743`, ported in
//! [`crate::chfi3d::builder::make_extremities`]) with the stripe list of the
//! vertex instead of a single stripe.

use std::sync::Arc;

use occt_core::gp::GpPnt2d;

use crate::chfi3d::builder_0_helpers::{
    ch_fi3d_compute_arete, ch_fi3d_index_of_surf_data, ch_fi3d_index_point_in_ds, ComputeArete,
};
use crate::topopebrep::ds::Curve as DsCurve;

use super::ChFi3dBuilder;

impl ChFi3dBuilder {
    /// `PerformSingularCorner(Index)` (`cxx:682-757`).
    pub fn perform_singular_corner(&mut self, index: usize) -> Result<(), String> {
        // `const TopoDS_Vertex& Vtx = myVDataMap.FindKey(Index);`
        let vtx = self
            .my_v_data_map
            .find_key(index)
            .cloned()
            .ok_or_else(|| "PerformSingularCorner: myVDataMap index out of range".to_string())?;
        let stripes = self
            .my_v_data_map
            .find_from_index(index)
            .map(|l| l.to_vec())
            .unwrap_or_default();

        let mut ivtx = 0;
        for (i, stripe) in stripes.iter().enumerate() {
            // `int sens = 0; num = ChFi3d_IndexOfSurfData(Vtx, stripe, sens);
            // isfirst = (sens == 1); Fd = ...Sequence().Value(num);`
            let mut sens = 0;
            let num = ch_fi3d_index_of_surf_data(&vtx, stripe, &mut sens)?;
            let isfirst = sens == 1;
            let fd: Arc<crate::chfi_ds::SurfData> = {
                let guard = stripe
                    .read()
                    .map_err(|_| "PerformSingularCorner: poisoned stripe".to_string())?;
                let idx = usize::try_from(num)
                    .ok()
                    .and_then(|n| n.checked_sub(1))
                    .ok_or_else(|| "PerformSingularCorner: bad SurfData index".to_string())?;
                guard.set_of_surf_data().get(idx).cloned().ok_or_else(|| {
                    "PerformSingularCorner: SurfData index out of range".to_string()
                })?
            };

            // `const ChFiDS_CommonPoint& CV1 = Fd->Vertex(isfirst, 1);`
            let cv1 = fd.vertex(isfirst, 1);
            let cv2 = fd.vertex(isfirst, 2);

            // `if (CV1.Point().IsEqual(CV2.Point(), 0))` (`cxx:703`).
            if !cv1
                .point()
                .is_equal_with_tolerance(cv2.point(), 0.0)
            {
                continue;
            }

            // `if (i == 0) Ivtx = ChFi3d_IndexPointInDS(CV1, DStr);`
            if i == 0 {
                ivtx = ch_fi3d_index_point_in_ds(cv1, self.my_ds.change_ds());
            }

            // `VOnS1 = InterferenceOnS1().PCurveOnSurf()->Value(First/Last)`
            // (`cxx:718-731`).
            let fi1 = fd.interference_on_s1();
            let fi2 = fd.interference_on_s2();
            let pc1 = fi1.p_curve_on_surf().ok_or_else(|| {
                "PerformSingularCorner: InterferenceOnS1 has no p-curve".to_string()
            })?;
            let pc2 = fi2.p_curve_on_surf().ok_or_else(|| {
                "PerformSingularCorner: InterferenceOnS2 has no p-curve".to_string()
            })?;
            let (par1, par2) = if isfirst {
                (fi1.first_parameter(), fi2.first_parameter())
            } else {
                (fi1.last_parameter(), fi2.last_parameter())
            };
            let von_s1: GpPnt2d = pc1.value(par1);
            let von_s2: GpPnt2d = pc2.value(par2);

            // `ChFi3d_ComputeArete(CV1, VOnS1, CV2, VOnS2,
            // DStr.Surface(Fd->Surf()).Surface(), ...)` (`cxx:734-746`).
            let surface = self
                .my_ds
                .ds()
                .surface(fd.surf())
                .surface()
                .cloned()
                .ok_or_else(|| {
                    "PerformSingularCorner: DStr.Surface(Surf()).Surface() is null".to_string()
                })?;
            let mut arete = ComputeArete::default();
            let (tol3d, tol2d) = (self.tolapp3d, self.tolapp2d);
            ch_fi3d_compute_arete(
                &cv1, &von_s1, &cv2, &von_s2, &surface, &mut arete, tol3d, tol2d, 0,
            )?;

            // `Crv = TopOpeBRepDS_Curve(C3d, tolreached); Icurv =
            // DStr.AddCurve(Crv);` (`cxx:746-748`).
            let mut crv = DsCurve::with_curve(arete.c3d.clone(), arete.tolreached, false);
            let icurv = self.my_ds.change_ds().add_curve(&mut crv);

            // `stripe->SetCurve(Icurv, isfirst); SetParameters(isfirst,
            // Pardeb, Parfin); ChangePCurve(isfirst) = PCurv;
            // SetIndexPoint(Ivtx, isfirst, 1); SetIndexPoint(Ivtx, isfirst, 2);`
            // (`cxx:749-753`).
            let mut guard = stripe
                .write()
                .map_err(|_| "PerformSingularCorner: poisoned stripe".to_string())?;
            guard.set_curve(icurv, isfirst);
            guard.set_parameters(isfirst, arete.pardeb, arete.parfin);
            *guard.change_p_curve(isfirst) = arete.p_curve;
            guard.set_index_point(ivtx, isfirst, 1);
            guard.set_index_point(ivtx, isfirst, 2);
        }
        Ok(())
    }
}
