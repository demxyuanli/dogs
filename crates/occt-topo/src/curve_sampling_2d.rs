//! Re-export of the OCCT 2D-curve sample counts, now hosted in `occt-geom2d`
//! (`geom2d_int::curve_sampling`) so the `Geom2dInt` port there can use them.
//! See that module for the OCCT line references.
pub(crate) use occt_geom2d::geom2d_int::nb_samples;
