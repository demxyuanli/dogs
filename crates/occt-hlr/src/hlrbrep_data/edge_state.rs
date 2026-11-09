//! Local 3D state and level queries of `HLRBRep_Data` - `EdgeState`
//! (`HLRBRep_Data.cxx:1580-1656`), `HidingStartLevel` (`:1658-1735`) and
//! `Compare` (`:1737-1744`).
//!
//! `EdgeState(p1, p2)` answers the two states of the current edge around the
//! crossing at parameter `p1` on `myLE` / `p2` on `myFE`: it lifts the surface
//! normal at `(p2)` - through `HLRBRep_EdgeFaceTool::UVPoint` and `mySLProps` -
//! and takes the sign of that normal against the edge tangent (`D1` at `p1`),
//! with the viewing vector flipped to point away from the eye. The sign is
//! thresholded by `myToler * 10`, so a tangent almost lying in the face answers
//! `ON`/`ON` rather than one of the two crossing states.
//!
//! `HidingStartLevel` finds the parameter of the first point of `ED` that can
//! still be hidden - the mid-point between the two "outermost" interference
//! parameters that do not exceed the edge end - classifies that point with
//! `Classify(..., LevelFlag = true, ...)`, and then walks the interference
//! list again, summing the levels of everything strictly before that
//! parameter (within `1%` of the edge range): a `FORWARD` interference
//! subtracts, a `REVERSED` one adds. The result is the number of hiding levels
//! at the start of the visible part.
//!
//! `Compare` is `Classify(..., LevelFlag = false, Level = 0, param = 0.)`: it
//! only answers whether the edge start is `IN`, `ON` or `OUT`.

use occt_core::gp::{GpDir, GpPnt2d};

use occt_topo::abs::{Orientation, State};
use crate::hlr_algo::HlrAlgoInterference;

use super::HlrBRepData;

impl HlrBRepData {
    /// `EdgeState(p1, p2, stbef, staft)` (`cxx:1580-1656`).
    ///
    /// Returns `(stbef, staft)` in place of the two OCCT out-parameters. Both
    /// the "no UV point" and the "no normal" arms answer `OUT`/`OUT`
    /// (`cxx:1637-1652`).
    pub fn edge_state(&mut self, p1: f64, p2: f64) -> (State, State) {
        let le = self.my_le as usize;
        let fe = self.my_fe as usize;

        // `if (HLRBRep_EdgeFaceTool::UVPoint(p2, myFEGeom, iFaceGeom, pu, pv))`
        // (`cxx:1590-1591`).
        let uv = match self.i_face_geom {
            Some(idx) => crate::hlrbrep_edge_face_tool::uv_point(
                p2,
                self.my_e_data[fe].geometry(),
                self.my_f_data[idx as usize].geometry(),
            ),
            None => None,
        };
        let Some((pu, pv)) = uv else {
            return (State::Out, State::Out);
        };

        self.my_sl_props.set_parameters(pu, pv);
        if !self.my_sl_props.is_normal_defined() {
            return (State::Out, State::Out);
        }
        let Some(&nrm_face) = self.my_sl_props.normal() else {
            return (State::Out, State::Out);
        };

        // `((HLRBRep_Curve*)myLEGeom)->D1(p1, Pbid, TngEdge)` (`cxx:1598-1599`).
        let (pbid, tng_edge) = self.my_e_data[le].geometry().d1_3d(p1);

        // `cxx:1601-1617`: the viewing vector, `(0, 0, NZ)` in the axonometric
        // case and the eye ray in the perspective one, then moved into the
        // model frame by the inverse transformation.
        let ti = self.my_proj.inverted_transformation();
        let mut v = if self.my_proj.perspective() {
            let mut p2d = GpPnt2d::new(0.0, 0.0);
            self.my_proj.project(&pbid, &mut p2d);
            let focus = match self.my_proj.focus() {
                Ok(f) => f,
                Err(_) => return (State::Out, State::Out),
            };
            match GpDir::new(p2d.x(), p2d.y(), -focus) {
                Ok(d) => d,
                Err(_) => return (State::Out, State::Out),
            }
        } else {
            // `gp_Dir(gp_Dir::D::NZ)` = `(0, 0, -1)`.
            match GpDir::new(0.0, 0.0, -1.0) {
                Ok(d) => d,
                Err(_) => return (State::Out, State::Out),
            }
        };
        v.transform(ti);

        let mut nrm = nrm_face;
        if nrm.dot(&v) > 0.0 {
            nrm.reverse();
        }

        // `scal = (TngEdge.SquareMagnitude() > 1.e-10) ? NrmFace.Dot(gp_Dir(TngEdge)) : 0.`
        // (`cxx:1619`).
        let scal = if tng_edge.square_magnitude() > 1.0e-10 {
            match GpDir::new(tng_edge.x(), tng_edge.y(), tng_edge.z()) {
                Ok(td) => nrm.dot(&td),
                Err(_) => 0.0,
            }
        } else {
            0.0
        };

        // `myToler * 10` is a `float` product in OCCT, promoted for the compare.
        let tol10 = (self.my_toler * 10.0f32) as f64;
        if scal > tol10 {
            (State::In, State::Out)
        } else if scal < -tol10 {
            (State::Out, State::In)
        } else {
            (State::On, State::On)
        }
    }

    /// `HidingStartLevel(E, ED, IL)` (`cxx:1658-1735`).
    ///
    /// `IL` is the interference list the hider accumulated for edge `E`
    /// (`NCollection_List<HLRAlgo_Interference>` in OCCT).
    pub fn hiding_start_level(&mut self, e: i32, il: &[HlrAlgoInterference]) -> i32 {
        let ed = e as usize;

        // `sta` / `end` are the **3D** parameters of the edge's 2D bounds
        // (`cxx:1665-1667`).
        let (mut sta, mut end) = {
            let ec = self.my_e_data[ed].geometry();
            (
                ec.parameter3d(ec.first_parameter()),
                ec.parameter3d(ec.last_parameter()),
            )
        };
        let tolpar = (end - sta) * 0.01;

        // `cxx:1671-1690`: keep shrinking the interval from the nearer side;
        // an interference past the edge end ends the scan.
        let mut more = true;
        for int in il {
            if !more {
                break;
            }
            let param = int.intersection().parameter();
            if param > end {
                more = false;
            } else if (param - sta).abs() > (param - end).abs() {
                end = param;
            } else {
                sta = param;
            }
        }

        let param = 0.5 * (sta + end);
        let mut level = 0;
        self.classify(e, true, &mut level, param);

        // `cxx:1696-1733`: replay the list and accumulate the levels strictly
        // before `param`, tolerating `tolpar`.
        more = true;
        for int in il {
            if !more {
                break;
            }
            let p = int.intersection().parameter();
            if p < param - tolpar {
                match int.transition() {
                    Orientation::Forward => level -= int.intersection().level(),
                    Orientation::Reversed => level += int.intersection().level(),
                    // `TopAbs_EXTERNAL` / `TopAbs_INTERNAL` / `default: break;`.
                    _ => {}
                }
            } else if p > param + tolpar {
                more = false;
            }
            // `|p - param| <= tolpar` is the "Bad Parameter" debug arm.
        }
        level
    }

    /// `Compare(E, ED)` (`cxx:1737-1744`): `Classify(E, ED, false, 0, 0.)`.
    pub fn compare(&mut self, e: i32) -> State {
        let mut level = 0;
        self.classify(e, false, &mut level, 0.0)
    }
}
