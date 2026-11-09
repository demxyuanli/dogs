//! `HLRBRep_Surface`'s side / above predicates:
//! `SideRowsOfPoles` (`HLRBRep_Surface.cxx:72-146`), `IsSide`
//! (`HLRBRep_Surface.cxx:148-240`) and `IsAbove`
//! (`HLRBRep_Surface.cxx:244-306`).
//!
//! Notes on the port:
//! * The C++ `SideRowsOfPoles` dereferences `myProj` unconditionally. Every
//!   caller (`IsSide`) runs under a projector, so a missing projector answers
//!   `false` here instead of being undefined behaviour.
//! * The pole grids are 1-based `NCollection_Array2` in OCCT; the port keeps
//!   the same U-major indexing with 0-based slices (`poles[iu][iv]`).

use occt_core::gp::GpPnt;
use occt_geom::gprop_pequation::GPropPEquation;

use occt_topo::geom_abs_curve_type::GeomAbsCurveType;
use occt_topo::geom_bnd_lib_surface3d::GeomAbsSurfaceType;
use crate::hlrbrep_curve::HlrBRepCurve;

use super::HlrBRepSurface;

impl HlrBRepSurface {
    /// `SideRowsOfPoles(tol, nbuPoles, nbvPoles, Pnt)`
    /// (`HLRBRep_Surface.cxx:72-146`). `poles` is the projector-space pole grid
    /// in U-major order; it is transformed in place exactly as OCCT does.
    pub fn side_rows_of_poles(
        &self,
        tol: f64,
        nbu_poles: i32,
        nbv_poles: i32,
        poles: &mut [Vec<GpPnt>],
    ) -> bool {
        let Some(proj) = self.projector() else {
            return false;
        };
        if nbu_poles <= 0 || nbv_poles <= 0 {
            return false;
        }
        let nu = nbu_poles as usize;
        let nv = nbv_poles as usize;
        if poles.len() < nu || poles.iter().take(nu).any(|row| row.len() < nv) {
            return false;
        }
        let trsf = proj.transformation();
        let tole = tol;

        for iu in 0..nu {
            for iv in 0..nv {
                poles[iu][iv].transform(trsf);
            }
        }

        // Side iso u?
        let mut result = true;
        for iu in 0..nu {
            if !result {
                break;
            }
            let (x0, y0) = (poles[iu][0].x(), poles[iu][0].y());
            for iv in 1..nv {
                if !result {
                    break;
                }
                result = (poles[iu][iv].x() - x0).abs() < tole
                    && (poles[iu][iv].y() - y0).abs() < tole;
            }
        }
        if result {
            return result;
        }

        // Side iso v?
        result = true;
        for iv in 0..nv {
            if !result {
                break;
            }
            let (x0, y0) = (poles[0][iv].x(), poles[0][iv].y());
            for iu in 1..nu {
                if !result {
                    break;
                }
                result = (poles[iu][iv].x() - x0).abs() < tole
                    && (poles[iu][iv].y() - y0).abs() < tole;
            }
        }
        if result {
            return result;
        }

        // Are the poles in a side plane?
        let mut flat = Vec::with_capacity(nu * nv);
        for iu in 0..nu {
            for iv in 0..nv {
                flat.push(poles[iu][iv]);
            }
        }
        let pl = GPropPEquation::new(&flat, tol);
        if pl.is_planar() {
            if let Some(plane) = pl.plane() {
                result = plane.axis().direction().z().abs() < 0.0001;
            }
        }

        result
    }

    /// `IsSide(tolF, toler)` (`HLRBRep_Surface.cxx:148-240`).
    pub fn is_side(&self, tol_f: f64, toler: f64) -> bool {
        let Some(proj) = self.projector() else {
            return false;
        };
        let trsf = proj.transformation();

        match self.get_type() {
            GeomAbsSurfaceType::Plane => {
                // `gp_Pln Pl = Plane(); gp_Ax1 A = Pl.Axis();`
                let Some(pl) = self.plane() else {
                    return false;
                };
                let axis = pl.axis();
                let mut pt = *axis.location();
                let mut d = *axis.direction();
                pt.transform(trsf);
                d.transform(trsf);
                let r = if proj.perspective() {
                    let focus = proj.focus().unwrap_or(0.0);
                    d.z() * focus - (d.x() * pt.x() + d.y() * pt.y() + d.z() * pt.z())
                } else {
                    d.z()
                };
                r.abs() < toler
            }
            GeomAbsSurfaceType::Cylinder => {
                if proj.perspective() {
                    return false;
                }
                let Some(cyl) = self.cylinder() else {
                    return false;
                };
                let axis = cyl.axis();
                let mut d = *axis.direction();
                d.transform(trsf);
                let r = (d.x() * d.x() + d.y() * d.y()).sqrt();
                r < toler
            }
            GeomAbsSurfaceType::Cone => {
                if !proj.perspective() {
                    return false;
                }
                let Some(cone) = self.cone() else {
                    return false;
                };
                let mut pt = cone.apex();
                pt.transform(trsf);
                let focus = proj.focus().unwrap_or(0.0);
                const CONE_TOL: f64 = 0.001;
                // `gp_Pnt::IsEqual(P, tol)` is `Distance(P) <= tol`
                // (`gp_Pnt.cxx:75-78`), not a strict comparison.
                pt.distance(&GpPnt::new(0.0, 0.0, focus)) <= CONE_TOL
            }
            GeomAbsSurfaceType::BezierSurface => {
                if proj.perspective() {
                    return false;
                }
                let Some(surf) = self.surface() else {
                    return false;
                };
                let Some((flat, nu, nv)) = surf.bezier_control_grid() else {
                    return false;
                };
                let Some(mut poles) = grid_from_flat(&flat, nu, nv) else {
                    return false;
                };
                let (nu, nv) = (nu as i32, nv as i32);
                self.side_rows_of_poles(tol_f, nu, nv, &mut poles)
            }
            GeomAbsSurfaceType::BSplineSurface => {
                if proj.perspective() {
                    return false;
                }
                let Some(surf) = self.surface() else {
                    return false;
                };
                let Some(src) = surf.bspline_surface_poles() else {
                    return false;
                };
                let mut poles: Vec<Vec<GpPnt>> = src.to_vec();
                let nu = poles.len() as i32;
                let nv = poles.first().map(|r| r.len()).unwrap_or(0) as i32;
                self.side_rows_of_poles(tol_f, nu, nv, &mut poles)
            }
            _ => false,
        }
    }

    /// `IsAbove(back, A, tol)` (`HLRBRep_Surface.cxx:244-306`).
    pub fn is_above(&self, back: bool, a: &HlrBRepCurve, tol: f64) -> bool {
        if self.get_type() != GeomAbsSurfaceType::Plane {
            return false;
        }
        let Some(pl) = self.plane() else {
            return false;
        };
        let (ca, cb, cc, cd) = pl.coefficients();
        let signed = |p: &GpPnt, back: bool| {
            let dd = ca * p.x() + cb * p.y() + cc * p.z() + cd;
            if back {
                -dd
            } else {
                dd
            }
        };

        let u1 = a.parameter3d(a.first_parameter());
        let u2 = a.parameter3d(a.last_parameter());
        let mut u = u1;

        if signed(&a.d0_3d(u), back) < -tol {
            return false;
        }

        if a.get_type() != GeomAbsCurveType::Line {
            let nb_pnt = 30;
            let step = (u2 - u1) / (nb_pnt + 1) as f64;
            for _ in 1..=nb_pnt {
                u += step;
                if signed(&a.d0_3d(u), back) < -tol {
                    return false;
                }
            }
        }

        u = u2;
        if signed(&a.d0_3d(u), back) < -tol {
            return false;
        }

        true
    }
}

/// Row-major `poles[iu * nv + iv]` (`Geom_BezierSurface.cxx:444-469`) as the
/// U-major grid `SideRowsOfPoles` walks.
fn grid_from_flat(flat: &[GpPnt], nu: usize, nv: usize) -> Option<Vec<Vec<GpPnt>>> {
    if nu == 0 || nv == 0 || flat.len() < nu * nv {
        return None;
    }
    let mut grid = Vec::with_capacity(nu);
    for iu in 0..nu {
        grid.push(flat[iu * nv..iu * nv + nv].to_vec());
    }
    Some(grid)
}
