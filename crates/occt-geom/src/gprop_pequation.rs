//! `GProp_PEquation` (`ModelingData/TKGeomBase/GProp/GProp_PEquation.hxx:40-119`,
//! `.cxx:28-226`).
//!
//! Principal-axis analysis of a point cloud: eigendecomposition of the central
//! inertia matrix, extents of the cloud along the three principal axes, then a
//! dimension count against the tolerance.
//!
//! `MathLin::Jacobi(aMat, true)` (`MathLin_Jacobi.hxx:52-141`) delegates to
//! `::Jacobi` (`math_Recipes.cxx:91-200`); [`occt_math::jacobi::Jacobi`] is that
//! same routine - same threshold schedule (`0.2*sm/(n*n)` for the first three
//! sweeps, `0.0` afterwards), same `100*|a(iq,ip)|` skip test, same selection
//! sort on the eigenvalues - so it stands in for `MathLin::Jacobi` here.
//!
//! The class lives in `occt-geom` rather than beside `occt_core::gprop::gprops`:
//! that crate is the one that may use `occt-math`, and the `PEquation` exported
//! by `occt_core::gprop::gprops` is an explicit plane-equation helper, not this
//! class.

use occt_core::gprop::gprops::PGProps;
use occt_core::gp::{GpAx1, GpAx3, GpDir, GpLin, GpPln, GpPnt, GpVec};
use occt_math::jacobi::Jacobi;
use occt_math::MathMatrix;

/// `GProp_PEquation::Type` (`hxx:46-54`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PEquationType {
    /// `Type::None` - the eigendecomposition did not converge.
    None,
    /// All points coincident within tolerance.
    Point,
    /// Points collinear within tolerance.
    Line,
    /// Points coplanar within tolerance.
    Plane,
    /// Points span 3D space.
    Space,
}

/// `GProp_PEquation`.
#[derive(Debug, Clone)]
pub struct GPropPEquation {
    /// `Type myType` (`hxx:109`).
    type_: PEquationType,
    /// `gp_Pnt myG` (`hxx:110`): mean point, plane origin or box corner.
    g: GpPnt,
    /// `gp_Vec myV1` (`hxx:111`): plane normal / line direction / box edge 1.
    v1: GpVec,
    /// `gp_Vec myV2` (`hxx:112`): box edge 2 (Space only).
    v2: GpVec,
    /// `gp_Vec myV3` (`hxx:113`): box edge 3 (Space only).
    v3: GpVec,
    /// `gp_Pnt myBarycentre` (`hxx:115`), always valid after construction.
    barycentre: GpPnt,
    /// `gp_Vec myAxes[3]` (`hxx:116`): unit principal axes.
    axes: [GpVec; 3],
    /// `double myExtents[3]` (`hxx:117`): extent along each principal axis.
    extents: [f64; 3],
}

impl GPropPEquation {
    /// `GProp_PEquation(const NCollection_Array1<gp_Pnt>&, double theTol)`
    /// (`cxx:28-177`).
    pub fn new(points: &[GpPnt], tol: f64) -> Self {
        // `GProp_PGProps aProps(thePnts); myBarycentre = aProps.CentreOfMass();`
        // (`cxx:29-30`). `CenterOfMass` is the absolute centre; `MatrixOfInertia`
        // is the central matrix, as `GProp_GProps::MatrixOfInertia` documents.
        let props = PGProps::from_points(points);
        let barycentre = props.center_of_mass();
        let g = barycentre;
        let inertia = props.matrix_of_inertia();

        // `math_Matrix aMat(1, 3, 1, 3)` filled from `anInertia` (`cxx:34-42`).
        let mut mat = MathMatrix::new(1, 3, 1, 3);
        for i in 0..3 {
            for j in 0..3 {
                mat.set_value(i + 1, j + 1, inertia.m[i][j]);
            }
        }
        let eigen = Jacobi::new(&mat);

        // `if (!anEigen.IsDone())` (`cxx:45-55`): identity axes, zero extents.
        if !eigen.is_done {
            return Self {
                type_: PEquationType::None,
                g,
                v1: GpVec::new(1.0, 0.0, 0.0),
                v2: GpVec::new(0.0, 1.0, 0.0),
                v3: GpVec::new(0.0, 0.0, 1.0),
                barycentre,
                axes: [
                    GpVec::new(1.0, 0.0, 0.0),
                    GpVec::new(0.0, 1.0, 0.0),
                    GpVec::new(0.0, 0.0, 1.0),
                ],
                extents: [0.0; 3],
            };
        }

        // `myAxes[k]` is column `k+1` of the eigenvector matrix (`cxx:56-60`).
        let mut axes = [GpVec::new(0.0, 0.0, 0.0); 3];
        for (k, ax) in axes.iter_mut().enumerate() {
            *ax = GpVec::new(
                eigen.vectors.value(1, k + 1),
                eigen.vectors.value(2, k + 1),
                eigen.vectors.value(3, k + 1),
            );
        }

        // Projections of every point on the three axes, relative to the
        // barycentre (`cxx:62-108`).
        let (xg, yg, zg) = (g.x(), g.y(), g.z());
        let mut max = [f64::MIN, f64::MIN, f64::MIN];
        let mut min = [f64::MAX, f64::MAX, f64::MAX];
        for p in points {
            let dx = p.x() - xg;
            let dy = p.y() - yg;
            let dz = p.z() - zg;
            for k in 0..3 {
                let d = dx * axes[k].x() + dy * axes[k].y() + dz * axes[k].z();
                if d > max[k] {
                    max[k] = d;
                }
                if d < min[k] {
                    min[k] = d;
                }
            }
        }

        let extents = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];

        // `cxx:111-128`: dimension count and the `aDimCode` bitmap the type
        // switch below reads back.
        let mut dimension = 3;
        let mut dim_code = 0;
        if extents[0].abs() <= tol {
            dimension -= 1;
            dim_code = 1;
        }
        if extents[1].abs() <= tol {
            dimension -= 1;
            dim_code = 2 * (dim_code + 1);
        }
        if extents[2].abs() <= tol {
            dimension -= 1;
            dim_code = 3 * (dim_code + 1);
        }

        // `cxx:129-176`.
        let mut type_ = PEquationType::None;
        let mut v1 = GpVec::new(0.0, 0.0, 0.0);
        let mut v2 = GpVec::new(0.0, 0.0, 0.0);
        let mut v3 = GpVec::new(0.0, 0.0, 0.0);
        let mut g = g;
        match dimension {
            0 => type_ = PEquationType::Point,
            1 => {
                type_ = PEquationType::Line;
                v1 = if dim_code == 4 {
                    axes[2]
                } else if dim_code == 6 {
                    axes[1]
                } else {
                    axes[0]
                };
            }
            2 => {
                type_ = PEquationType::Plane;
                v1 = if dim_code == 1 {
                    axes[0]
                } else if dim_code == 2 {
                    axes[1]
                } else {
                    axes[2]
                };
            }
            3 => {
                type_ = PEquationType::Space;
                // `myG.SetXYZ(myG.XYZ() + aMin1*Axes[0] + ...)` (`cxx:169-170`).
                g = GpPnt::from_xyz(
                    &g.coord
                        .added(&axes[0].coord.multiplied(min[0]))
                        .added(&axes[1].coord.multiplied(min[1]))
                        .added(&axes[2].coord.multiplied(min[2])),
                );
                v1 = GpVec::new(
                    extents[0] * axes[0].x(),
                    extents[0] * axes[0].y(),
                    extents[0] * axes[0].z(),
                );
                v2 = GpVec::new(
                    extents[1] * axes[1].x(),
                    extents[1] * axes[1].y(),
                    extents[1] * axes[1].z(),
                );
                v3 = GpVec::new(
                    extents[2] * axes[2].x(),
                    extents[2] * axes[2].y(),
                    extents[2] * axes[2].z(),
                );
            }
            _ => {}
        }

        Self {
            type_,
            g,
            v1,
            v2,
            v3,
            barycentre,
            axes,
            extents,
        }
    }

    /// `GetType()` (`hxx:62`).
    pub fn get_type(&self) -> PEquationType {
        self.type_
    }

    /// `IsPlanar()` (`hxx:66`).
    pub fn is_planar(&self) -> bool {
        self.type_ == PEquationType::Plane
    }

    /// `IsLinear()` (`hxx:69`).
    pub fn is_linear(&self) -> bool {
        self.type_ == PEquationType::Line
    }

    /// `IsPoint()` (`hxx:72`).
    pub fn is_point(&self) -> bool {
        self.type_ == PEquationType::Point
    }

    /// `IsSpace()` (`hxx:75`).
    pub fn is_space(&self) -> bool {
        self.type_ == PEquationType::Space
    }

    /// `Plane()` (`cxx:182-188`): `gp_Pln(myG, myV1)`. `None` where OCCT throws
    /// `Standard_NoSuchObject` because the cloud is not planar.
    pub fn plane(&self) -> Option<GpPln> {
        if !self.is_planar() {
            return None;
        }
        let n = GpDir::from_xyz(&self.v1.coord).ok()?;
        Some(GpPln::new(GpAx3::from_ax1(&GpAx1::new(self.g, n))))
    }

    /// `Line()` (`cxx:192-198`): `gp_Lin(myG, gp_Dir(myV1))`. `None` where OCCT
    /// throws because the cloud is not linear.
    pub fn line(&self) -> Option<GpLin> {
        if !self.is_linear() {
            return None;
        }
        let d = GpDir::from_xyz(&self.v1.coord).ok()?;
        Some(GpLin::new(GpAx1::new(self.g, d)))
    }

    /// `Point()` (`cxx:202-210`). `None` where OCCT throws because the cloud is
    /// not coincident.
    pub fn point(&self) -> Option<GpPnt> {
        if !self.is_point() {
            return None;
        }
        Some(self.g)
    }

    /// `Box(P, V1, V2, V3)` (`cxx:214-225`). `None` where OCCT throws because
    /// the cloud spans less than 3D space.
    pub fn box_edges(&self) -> Option<(GpPnt, GpVec, GpVec, GpVec)> {
        if !self.is_space() {
            return None;
        }
        Some((self.g, self.v1, self.v2, self.v3))
    }

    /// `Barycentre()` (`hxx:97`).
    pub fn barycentre(&self) -> GpPnt {
        self.barycentre
    }

    /// `PrincipalAxis(theIndex)` (`hxx:100`): `theIndex` is 1, 2 or 3.
    pub fn principal_axis(&self, the_index: usize) -> Option<GpVec> {
        self.axes.get(the_index.checked_sub(1)?).copied()
    }

    /// `Extent(theIndex)` (`hxx:104`): `theIndex` is 1, 2 or 3.
    pub fn extent(&self, the_index: usize) -> Option<f64> {
        self.extents.get(the_index.checked_sub(1)?).copied()
    }
}
