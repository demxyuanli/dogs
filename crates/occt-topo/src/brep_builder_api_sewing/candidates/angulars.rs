//! `EvaluateAngulars` (BRepBuilderAPI_Sewing.cxx:1189-1280).

use occt_core::gp::GpVec;
use occt_core::precision::RESOLUTION;

use crate::brep_tool::BRepTool;
use crate::shape::{Edge, Face, TopoShape};

use super::abscissa::uniform_abscissa_range;
use super::super::Sewing;


impl Sewing {
    /// `EvaluateAngulars(sequenceSec, secForward, tabAng, indRef)`
    /// (`cxx:1189-1280`): the mean angle between the surface normal of the
    /// reference section and each following one, sampled at `4` uniform
    /// abscissa points.
    #[allow(dead_code)] // cxx:1769 is commented out: the only OCCT caller is dead source.
    pub(crate) fn evaluate_angulars(
        &self,
        sequence_sec: &[TopoShape],
        sec_forward: &[bool],
        tab_ang: &mut [f64],
        ind_ref: usize,
    ) {
        for a in tab_ang.iter_mut() {
            *a = -1.0;
        }
        let npt = 4usize;
        let leng_sec = sequence_sec.len();
        let mut norm_ref: Vec<GpVec> = vec![GpVec::zero(); npt];
        for i in ind_ref..=leng_sec {
            let edge = sequence_sec[i - 1].clone();
            let mut bnd = edge.clone();
            if let Some(b) = self.section_bound.get(&edge) {
                bnd = b;
            }
            let (face, c2d_opt) = if let Some(faces) = self.bound_faces.get(&bnd) {
                let face = Face(faces[0].clone());
                let c2d = crate::boptools_2d::curve_on_surface_range(&Edge(edge.clone()), &face)
                    .map(|(c, _, _)| c);
                (Some(face), c2d)
            } else if i == ind_ref {
                return;
            } else {
                continue;
            };
            let Some(face) = face else { continue };
            let Some(c2d) = c2d_opt else { continue };
            let Some(surf) = self.face_surface_located(&face) else {
                continue;
            };
            let Some(c3d) = BRepTool::edge_curve_adaptor_world(&Edge(edge.clone())) else {
                continue;
            };
            let (first, last) = BRepTool::edge_parameters(&Edge(edge));
            let Some(params) = uniform_abscissa_range(c3d.as_ref(), npt, first, last) else {
                continue;
            };
            let mut cumulate_angular = 0.0;
            let mut nb_computed = 0usize;
            for j in 1..=npt {
                let jj = if sec_forward.get(i - 1).copied().unwrap_or(false) || i == ind_ref {
                    j
                } else {
                    npt - j + 1
                };
                let p = c2d.d0(params[jj - 1]);
                let (_pnt, w1, w2) = surf.d1(p.x(), p.y());
                let n = w1.crossed(&w2);
                if i == ind_ref {
                    norm_ref[j - 1] = n;
                } else if n.magnitude() > RESOLUTION && norm_ref[j - 1].magnitude() > RESOLUTION {
                    nb_computed += 1;
                    let cos = n.dot(&norm_ref[j - 1]) / (n.magnitude() * norm_ref[j - 1].magnitude());
                    let mut angular = cos.clamp(-1.0, 1.0).acos();
                    if angular > std::f64::consts::FRAC_PI_2 {
                        angular = std::f64::consts::PI - angular;
                    }
                    cumulate_angular += angular;
                }
            }
            if nb_computed > 0 {
                tab_ang[i - 1] = cumulate_angular / nb_computed as f64;
            }
        }
    }
}
