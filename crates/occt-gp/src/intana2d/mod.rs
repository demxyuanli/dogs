//! `IntAna2d` — analytical intersection of 2D conics.
//! Source: `src/ModelingData/TKGeomBase/IntAna2d`.
pub mod ana_intersection;
pub mod conic;
pub mod int_point;
pub mod outils;

pub use ana_intersection::IntAna2dAnaIntersection;
pub use conic::IntAna2dConic;
pub use int_point::IntAna2dIntPoint;
pub use outils::{
    coord_ancien_repere, points_confondus, traitement_points_confondus, MyDirectPolynomialRoots,
};
