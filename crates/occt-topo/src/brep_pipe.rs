//! Phase 4 module: brep_pipe — sweep a profile along a path (MakePipe).
//!
//! Port of `BRepOffsetAPI_MakePipe`: extrude a planar profile polygon along a
//! polyline path, building ruled quad lateral faces between consecutive
//! sections and planar end caps.

use occt_core::gp::{GpPnt, GpVec, GpXyz};

use crate::abs::ShapeType;
use crate::builder::TopoBuilder;
use crate::loft;
use crate::shape::{Face, Shell, Solid, Wire};
use crate::topo_tools_full;

/// Result of a polyline sweep: the solid plus its boundary faces and the path.
#[derive(Debug, Clone)]
pub struct Pipe {
    pub solid: Option<Solid>,
    pub shell: Shell,
    pub faces: Vec<Face>,
    pub path: Vec<GpPnt>,
}

/// Shift every profile point by `offset`.
pub fn translate_profile(profile: &[GpPnt], offset: &GpVec) -> Vec<GpPnt> {
    profile
        .iter()
        .map(|p| GpPnt::new(p.x() + offset.x(), p.y() + offset.y(), p.z() + offset.z()))
        .collect()
}

/// Sweep a closed planar profile along a polyline path.
///
/// The profile is translated to each path point; consecutive sections are
/// connected with ruled quad faces and the two ends are capped. The profile
/// may be given with a duplicated closing point (first == last); the duplicate
/// is dropped.
pub fn pipe_along_polyline(profile: &[GpPnt], path: &[GpPnt]) -> Result<Pipe, String> {
    if path.len() < 2 {
        return Err("pipe_along_polyline: path needs >= 2 points".into());
    }
    let mut prof = profile.to_vec();
    if prof.len() >= 2 && prof[0].distance(prof.last().unwrap()) < 1e-9 {
        prof.pop();
    }
    if prof.len() < 3 {
        return Err("pipe_along_polyline: profile needs >= 3 points".into());
    }
    // One closed section per path point: the profile translated by that offset.
    let sections: Vec<Wire> = path
        .iter()
        .map(|p| {
            let sec = translate_profile(&prof, &GpVec::new(p.x(), p.y(), p.z()));
            loft::section_wire(&sec)
        })
        .collect::<Result<_, _>>()?;

    let lofted = loft::loft_sections(&sections)?;
    let faces = topo_tools_full::faces_of(&lofted.solid.0);
    let shell = topo_tools_full::shapes_of(&lofted.solid.0, ShapeType::Shell)
        .into_iter()
        .next()
        .map(Shell)
        .unwrap_or_else(|| TopoBuilder::new().make_shell(&faces));
    Ok(Pipe {
        solid: Some(lofted.solid),
        shell,
        faces,
        path: path.to_vec(),
    })
}

/// Volume of a swept pipe: the tessellated solid volume, with a bounding-box
/// fallback when the tessellation is unusable.
pub fn pipe_volume(pipe: &Pipe) -> f64 {
    if let Some(solid) = &pipe.solid {
        let v = crate::brep_gprop::volume(solid, 0.1);
        if v.is_finite() && v > 0.0 {
            return v;
        }
    }
    // ponytail: bbox fallback — only hit when the tessellated volume is unusable.
    let mut min = GpXyz::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut max = GpXyz::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for v in topo_tools_full::vertices_of(&pipe.shell.0) {
        let p = topo_tools_full::vertex_position(&v).coord;
        min.x = min.x.min(p.x);
        min.y = min.y.min(p.y);
        min.z = min.z.min(p.z);
        max.x = max.x.max(p.x);
        max.y = max.y.max(p.y);
        max.z = max.z.max(p.z);
    }
    (max.x - min.x).max(0.0) * (max.y - min.y).max(0.0) * (max.z - min.z).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topo_tools_full::faces_of;

    fn square_profile() -> Vec<GpPnt> {
        vec![
            GpPnt::new(0., 0., 0.),
            GpPnt::new(1., 0., 0.),
            GpPnt::new(1., 1., 0.),
            GpPnt::new(0., 1., 0.),
        ]
    }

    #[test]
    fn translate_profile_shifts_points() {
        let pts = [GpPnt::new(1., 2., 3.)];
        let shifted = translate_profile(&pts, &GpVec::new(1., 1., 1.));
        assert!(shifted[0].is_equal(&GpPnt::new(2., 3., 4.)));
    }

    #[test]
    fn pipe_square_along_two_point_path() {
        let profile = square_profile();
        let path = [GpPnt::new(0., 0., 0.), GpPnt::new(0., 0., 5.)];
        let pipe = pipe_along_polyline(&profile, &path).expect("pipe");
        assert_eq!(pipe.faces.len(), 6, "2 caps + 4 lateral");
        assert_eq!(faces_of(&pipe.shell.0).len(), 6);
        let vol = pipe_volume(&pipe);
        assert!((vol - 5.0).abs() < 0.2, "volume {vol} (expect 1×5)");
    }

    #[test]
    fn pipe_along_three_point_bent_path() {
        let profile = square_profile();
        let path = [
            GpPnt::new(0., 0., 0.),
            GpPnt::new(0., 0., 3.),
            GpPnt::new(0., 3., 6.),
        ];
        let pipe = pipe_along_polyline(&profile, &path).expect("bent pipe");
        assert_eq!(pipe.faces.len(), 10, "2 caps + 8 lateral");
        assert!(pipe_volume(&pipe) > 0.0);
    }

    #[test]
    fn pipe_rejects_short_inputs() {
        assert!(pipe_along_polyline(&square_profile(), &[GpPnt::zero()]).is_err());
        let two_pts = [GpPnt::zero(), GpPnt::new(1., 0., 0.)];
        assert!(pipe_along_polyline(&two_pts, &[GpPnt::zero(), GpPnt::new(0., 0., 1.)]).is_err());
    }
}
