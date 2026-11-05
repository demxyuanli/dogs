//! `PSurf.Perform()` and `GeomPlate_MakeApprox` of the `GeomPlate` section, and
//! the DS storage of the plate surface.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_CnCrn.cxx`
//! (OCCT 8.0.0), `cxx:3494-3626`:
//!
//! - `PSurf.Perform()` and the `IsDone` guard that opens the whole tail,
//! - the `Disc2dContour(4, S2d)` / `Disc3dContour(4, 0, S3d)` sampling, the
//!   `PlateG0Criterion` at `max(tolapp, 10 * G0Error)` and
//!   `GeomPlate_MakeApprox(plate, crit, tolapp, 9, 8)`,
//! - `DStr.AddSurface(TopOpeBRepDS_Surface(Surf, Mapp.ApproxError()))` with the
//!   `apperror = CriterionError * 1.1` the tails carry, the `myEVIMap` bind of
//!   `V1`, and the `Isurf` registration,
//! - the `deuxconges` write-back of the `V1` vertex on every non-sharp spine
//!   arc, with the reversed orientation stored so that a top-level re-read
//!   removes the living edges,
//! - `PlateOrientation` and the `SolidSurfaceInterference` of the stripe's
//!   solid on the plate surface.
//!
//! The `else` of the `IsDone` guard (`cxx:3891-3926`, the partial-result
//! branch) belongs to the following stage; this one answers `false` so the
//! caller can run it.
//!
//! Translation notes:
//! - `gbPlate = PSurf.Surface()` is [`BuildPlateSurface::surface`]; a missing
//!   handle is the `Standard_NullObject` the OCCT constructor would raise on,
//! - `nbcarreau = 9` / `degmax = 8` (`cxx:3507-3508`), and the two defaults of
//!   the OCCT constructor (`hxx:45-46`): `GeomAbs_C1` and `1.1`,
//! - `TopOpeBRepDS_Surface Tsurf(Surf, Mapp.ApproxError())`
//!   (`cxx:3541`): the DS tolerance is the *approximation* error while
//!   `apperror` (kept for the tails, `cxx:3525-3526`) is the *criterion*
//!   error times `1.1`,
//! - `indice = Fd->Surf()` (`cxx:3601`) feeds only the commented-out block of
//!   `cxx:3602-3616` and is not read afterwards, so the port does not compute
//!   it,
//! - `myEVIMap.IsBound(V1) / Bind / ChangeFind(V1).Append(Isurf)`
//!   (`cxx:3544-3549`) is the state's `(shape, index list)` table, the same
//!   reading as [`super::super::ChFi3dBuilder::change_evi`],
//! - `DStr.ChangeShapeInterferences(SolInd)` is
//!   [`crate::topopebrep::ds::DataStructure::change_shape_interferences_at`].

use std::sync::{Arc, RwLock};

use occt_core::kernel::geomabs::Shape;
use occt_geom::geom_plate::{BuildPlateSurface, MakeApprox, PlateG0Criterion};
use occt_geom::Surface;

use crate::abs::{Orientation, ShapeType};
use crate::brep_tool::BRepTool;
use crate::chfi3d::builder_0_helpers::ch_fi3d_fil_vertex_in_ds;
use crate::shape::TopoShape;
use crate::topo_tools::{explorer_shapes, is_same};
use crate::topopebrep::ds::{DataStructure, Interference, Kind, Surface as DsSurface, Transition};

use super::plate_helpers::plate_orientation;
use super::surf_index::indices;
use super::{spine_edge_at, MoreThreeCorner};

/// `double tolapp = 1.e-4` (`cxx:1312`) — the `Tol3d` of the approximation.
const TOLAPP: f64 = 1.0e-4;
/// `int nbcarreau = 9` (`cxx:3507`).
const NBCARREAU: i32 = 9;
/// `int degmax = 8` (`cxx:3508`).
const DEGMAX: i32 = 8;
/// `double coef = 1.1` (`cxx:3524`).
const COEF: f64 = 1.1;

impl MoreThreeCorner {
    /// `cxx:3494-3626`.
    ///
    /// `evi` is the builder's `myEVIMap`. Answers `true` when `PSurf.IsDone()`,
    /// `false` when the plate failed and the partial-result tail of
    /// `cxx:3891-3926` has to run instead.
    pub(super) fn plate_solve_approx(
        &mut self,
        ds: &mut DataStructure,
        evi: &mut Vec<(TopoShape, Vec<i32>)>,
        psurf: &mut BuildPlateSurface,
    ) -> Result<bool, String> {
        // `PSurf.Perform();` (`cxx:3494`).
        psurf.perform()?;

        // `if (PSurf.IsDone())` (`cxx:3506`).
        if !psurf.is_done() {
            return Ok(false);
        }

        // `occ::handle<GeomPlate_Surface> gpPlate = PSurf.Surface();`
        // (`cxx:3510`).
        let gp_plate = psurf
            .surface()
            .cloned()
            .ok_or_else(|| "PerformMoreThreeCorner: PSurf.Surface() is null".to_string())?;

        // `PSurf.Disc2dContour(4, S2d); PSurf.Disc3dContour(4, 0, S3d);`
        // (`cxx:3516-3517`).
        let mut s2d = Vec::new();
        psurf.disc2d_contour(&mut s2d)?;
        let mut s3d = Vec::new();
        psurf.disc3d_contour(0, &mut s3d)?;

        // `seuil = std::max(tolapp, 10 * PSurf.G0Error());
        //  GeomPlate_PlateG0Criterion critere(S2d, S3d, seuil);`
        // (`cxx:3518-3519`).
        let seuil = TOLAPP.max(10.0 * psurf.g0_error());
        let critere = PlateG0Criterion::with_defaults(
            s2d.iter().map(|p| (p.x(), p.y())).collect(),
            s3d,
            seuil,
        );

        // `GeomPlate_MakeApprox Mapp(gpPlate, critere, tolapp, nbcarreau,
        //  degmax);` (`cxx:3520-3521`), with the two defaults of
        // `hxx:45-46`.
        let mapp = MakeApprox::from_criterion(
            Arc::clone(&gp_plate),
            &critere,
            TOLAPP,
            NBCARREAU,
            DEGMAX,
            Shape::C1,
            COEF,
        )?;

        // `double coef = 1.1, apperror; apperror = Mapp.CriterionError() *
        //  coef;` (`cxx:3524-3526`).
        let apperror = mapp.criterion_error() * COEF;
        // `TopOpeBRepDS_Surface Tsurf(Surf, Mapp.ApproxError());`
        // (`cxx:3541`).
        let approx_error = mapp.approx_error();
        let surf: Arc<dyn Surface> = Arc::new(
            mapp.into_surface()
                .ok_or_else(|| "PerformMoreThreeCorner: Mapp.Surface() is null".to_string())?,
        );

        // `int Isurf = DStr.AddSurface(Tsurf);` (`cxx:3542`).
        let isurf = ds.add_surface(&DsSurface::with_tolerance(Arc::clone(&surf), approx_error));

        // `if (!myEVIMap.IsBound(V1)) { NCollection_List<int> li;
        //  myEVIMap.Bind(V1, li); } myEVIMap.ChangeFind(V1).Append(Isurf);`
        // (`cxx:3544-3549`).
        match evi.iter_mut().position(|(key, _)| is_same(key, &self.vtx)) {
            Some(k) => evi[k].1.push(isurf),
            None => evi.push((self.vtx.clone(), vec![isurf])),
        }

        // `int SolInd = CD.Value(0)->SolidIndex();` (`cxx:3551`).
        let sol_ind = {
            let stripe = self.stripe_at(0)?.clone();
            let guard = stripe
                .read()
                .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
            guard.solid_index()
        };

        // `cxx:3558-3593`: with two fillets on the vertex, the living edges that
        // arrive at the top are stored with their reversed orientation so that
        // a top-level re-read removes them from the DS.
        if self.deuxconges {
            let nedge = self.nedge as usize;
            for ic in 0..nedge {
                if self.sharp[ic] {
                    continue;
                }
                // `nbedge = CD.Value(ic)->Spine()->NbEdges();
                //  if (sens.Value(ic) == 1) Arcspine = Edges(1);
                //  else Arcspine = Edges(nbedge);` (`cxx:3565-3573`).
                let stripe = self.stripe_at(ic)?.clone();
                let arcspine = spine_edge_at(&stripe, self.sens[ic] == 1)?;

                // `int IArcspine = DStr.AddShape(Arcspine);` (`cxx:3576`).
                let i_arcspine = ds.add_shape(&arcspine.0);

                // `OVtx = FORWARD; for (ex.Init(Arcspine.Oriented(FORWARD),
                //  VERTEX); ...) if (V1.IsSame(ex.Current())) { OVtx =
                //  ex.Current().Orientation(); break; } OVtx =
                //  TopAbs::Reverse(OVtx);` (`cxx:3577-3586`).
                let mut ovtx = Orientation::Forward;
                for v in explorer_shapes(&arcspine.0.oriented(Orientation::Forward), ShapeType::Vertex)
                {
                    if is_same(&v, &self.vtx) {
                        ovtx = v.orientation();
                        break;
                    }
                }
                ovtx = ovtx.reversed();

                // `parVtx = BRep_Tool::Parameter(V1, Arcspine);` (`cxx:3587`).
                let par_vtx = BRepTool::parameter_of_vertex(&arcspine, &self.vtx);

                // `interfv = ChFi3d_FilVertexInDS(OVtx, IArcspine, IVtx,
                //  parVtx); DStr.ChangeShapeInterferences(IArcspine).Append(
                //  interfv);` (`cxx:3588-3590`).
                let interfv = ch_fi3d_fil_vertex_in_ds(ovtx, i_arcspine, self.ivtx, par_vtx);
                ds.change_shape_interferences_at(i_arcspine).push(interfv);
            }
        }

        // `Indices(nedge, 0, icplus, icmoins); isfirst = (sens.Value(0) == 1);
        //  const ChFiDS_SurfData& Fd = CD.Value(0)->SetOfSurfData()->Value(
        //  i.Value(0, icmoins)); indice = Fd->Surf(); orsurfdata =
        //  Fd->Orientation();` (`cxx:3597-3617`).
        let (_icplus, icmoins) = indices(self.nedge as usize, 0);
        let isfirst = self.sens[0] == 1;
        let fd = {
            let stripe = self.stripe_at(0)?.clone();
            let guard = stripe
                .read()
                .map_err(|_| "PerformMoreThreeCorner: poisoned stripe".to_string())?;
            let idx = self.i[0][icmoins];
            let pos = usize::try_from(idx)
                .ok()
                .and_then(|i| i.checked_sub(1))
                .ok_or_else(|| {
                    format!("PerformMoreThreeCorner: surf data index {idx} is out of range")
                })?;
            Arc::clone(guard.set_of_surf_data().get(pos).ok_or_else(|| {
                format!("PerformMoreThreeCorner: surf data index {idx} is out of range")
            })?)
        };
        let orsurfdata = fd.orientation();

        // `orplate = PlateOrientation(Surf, PSurf.Curves2d(), SumFaceNormalAtV1);`
        // (`cxx:3618`).
        let pc_arr = psurf.curves2d()?;
        let orplate = plate_orientation(&surf, &pc_arr, &self.sum_face_normal_at_v1)?;

        // `SSI = new TopOpeBRepDS_SolidSurfaceInterference(
        //  TopOpeBRepDS_Transition(orplate), SOLID, SolInd, SURFACE, Isurf);
        //  SolidInterfs.Append(SSI);` (`cxx:3620-3626`).
        let ssi = Arc::new(RwLock::new(Interference::solid_surface(
            Transition::from_orientation(orplate),
            Kind::Solid,
            sol_ind,
            Kind::Surface,
            isurf,
        )));
        ds.change_shape_interferences_at(sol_ind).push(ssi);

        // The locals the following stages read (`cxx:3526`, `:3542`, `:3618`).
        self.plate_surface = Some(surf);
        self.isurf = isurf;
        self.apperror = apperror;
        self.or_surfdata = orsurfdata;
        self.or_plate = orplate;
        self.is_first_plate = isfirst;

        Ok(true)
    }
}
