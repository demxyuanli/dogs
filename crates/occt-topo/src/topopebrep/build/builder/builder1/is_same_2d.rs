//! `TopOpeBRepBuild_Builder1::IsSame2d`.
//!
//! Source: `ModelingAlgorithms/TKBool/TopOpeBRepBuild/TopOpeBRepBuild_Builder1.cxx:2115-2223`
//! (OCCT 8.0.0).
//!
//! Only fully closed edges on a periodic face can be "the same in 2d": the
//! piece of the tool edge is translated in the `(u, v)` plane onto the piece of
//! the object edge, and when the two p-curves run the same way the translated
//! copy is stored as the first of the two seam representations.
//!
//! UNPORTED: nothing; the two `FC2D_CurveOnSurface` queries that OCCT
//! dereferences without a null test answer `1` here instead
//! (`cxx:2156`, `:2162`, `:2178`, `:2195`), see [`Builder1::is_same_2d`].

use std::sync::Arc;

use occt_core::gp::{GpTrsf2d, GpVec2d};
use occt_geom2d::curve::Curve2d;

use crate::abs::Orientation;
use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};
use crate::topopebrep::build::tools::{fc2d_curve_on_surface, write_edge_pcurves};

use super::super::Builder;

/// `gp_Resolution` (`cxx:2151`): the smallest translation worth applying.
const GP_RESOLUTION: f64 = 1.0e-10;

impl Builder {
    /// `IsSame2d(aSeq, aListOfPiecesOut2d)` (`cxx:2115-2223`).
    ///
    /// Returns `0` when the tool piece was moved (and appended to
    /// `a_list_of_pieces_out_2d`), `1` when nothing had to be done.
    pub(crate) fn is_same_2d(
        &self,
        a_seq: &[TopoShape],
        a_list_of_pieces_out_2d: &mut Vec<TopoShape>,
    ) -> i32 {
        if a_seq.len() < 6 {
            return 0;
        }

        let a_face_obj = a_seq[0].clone();
        let an_edge_obj = a_seq[1].clone();
        let a_piece_obj = a_seq[2].clone();
        let a_face_tool = a_seq[3].clone();
        let an_edge_tool = a_seq[4].clone();
        let a_piece_tool = a_seq[5].clone();

        let (Some(a_f_obj), Some(_a_f_tool)) = (Face::wrap(a_face_obj.clone()), Face::wrap(a_face_tool.clone()))
        else {
            return 1;
        };
        let (Some(an_e_obj), Some(an_e_tool)) = (
            Edge::wrap(an_edge_obj.clone()),
            Edge::wrap(an_edge_tool.clone()),
        ) else {
            return 1;
        };
        let (Some(a_p_obj), Some(a_p_tool)) = (
            Edge::wrap(a_piece_obj.clone()),
            Edge::wrap(a_piece_tool.clone()),
        ) else {
            return 1;
        };

        // `cxx:2138-2141`: the whole trick only applies to a periodic face.
        let Some(surf) = BRepTool::face_surface(&a_f_obj) else {
            return 1;
        };
        if !(surf.is_u_periodic() || surf.is_v_periodic()) {
            return 1;
        }

        // `cxx:2144-2147`: and only to fully closed edges (`Vf == Vl`).
        if !BRepTool::is_closed(&an_edge_obj) || !BRepTool::is_closed(&an_edge_tool) {
            return 1;
        }

        // `cxx:2156-2159`: the tool piece, sampled at `f1 * PAR_T + ...`.
        let Some((c2d_piece_tool, f1, l1)) = fc2d_curve_on_surface(&a_p_tool, &a_f_obj) else {
            return 1;
        };
        let par_p = f1 * super::PAR_T + (1.0 - super::PAR_T) * l1;
        let a_p_piece = c2d_piece_tool.d0(par_p);

        // `cxx:2162-2174`: where that point projects onto the tool edge.
        let Some((c2d, f, l)) = fc2d_curve_on_surface(&an_e_tool, &a_f_obj) else {
            return 1;
        };
        let Some(proj) = occt_geom2d::geom2d_api::project_point_on_curve_2d(
            c2d.as_ref(),
            &a_p_piece,
            occt_core::precision::CONFUSION,
        ) else {
            return 1;
        };
        let mut par_ori = proj.parameter;
        if par_ori < f {
            par_ori = 2.0 * std::f64::consts::PI + par_ori;
        }
        if par_ori > l {
            par_ori -= 2.0 * std::f64::consts::PI;
        }

        let a_uv2 = c2d.d0(par_ori);

        // `cxx:2179-2183`: the object piece, sampled the same way.
        let Some((c2d_piece_obj, f, l)) = fc2d_curve_on_surface(&a_p_obj, &a_f_obj) else {
            return 1;
        };
        let par = f * super::PAR_T + (1.0 - super::PAR_T) * l;
        let a_uv1 = c2d_piece_obj.d0(par);

        // `cxx:2184-2186`: `gp_Vec2d(aUV1, aUV2)`.
        let a_translate_v = GpVec2d::new(a_uv2.x() - a_uv1.x(), a_uv2.y() - a_uv1.y());

        if a_translate_v.magnitude() >= GP_RESOLUTION {
            // `cxx:2188-2191`: `C2DPieceTool->Copy()` then `Translate`.
            let mut tr = GpTrsf2d::identity();
            tr.set_translation_vec(&a_translate_v);
            let a_tr_c2d: Arc<dyn Curve2d> = Arc::from(c2d_piece_tool.transformed(&tr));

            let a_t_fuv = a_tr_c2d.d0(f1);
            let a_t_luv = a_tr_c2d.d0(l1);
            let a_tr_vec = GpVec2d::new(a_t_luv.x() - a_t_fuv.x(), a_t_luv.y() - a_t_fuv.y());

            // `cxx:2193-2203`: the direction of the original object edge.
            let Some((c2d_edge_obj, fo, lo)) = fc2d_curve_on_surface(&an_e_obj, &a_f_obj) else {
                return 1;
            };
            let a_o_fuv = c2d_edge_obj.d0(fo);
            let a_o_luv = c2d_edge_obj.d0(lo);
            let mut a_o_vec = GpVec2d::new(a_o_luv.x() - a_o_fuv.x(), a_o_luv.y() - a_o_fuv.y());
            if an_e_obj.0.orientation() == Orientation::Reversed {
                a_o_vec.reverse();
            }
            let is_tr_first = a_tr_vec.dot(&a_o_vec) <= 0.0;

            // `cxx:2205-2217`: two seam representations, order by direction.
            let tol_e = BRepTool::edge_tolerance(&a_p_tool);
            if is_tr_first {
                write_edge_pcurves(&a_p_tool, a_tr_c2d, c2d_piece_tool, &a_f_obj, tol_e);
            } else {
                write_edge_pcurves(&a_p_tool, c2d_piece_tool, a_tr_c2d, &a_f_obj, tol_e);
            }

            a_list_of_pieces_out_2d.push(a_p_tool.0.clone());
            return 0;
        }

        1
    }
}
