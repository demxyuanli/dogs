//! The three `Geom2dGcc` solution-class enums, which only `Geom2dGcc`
//! (`Geom2dGcc.cxx:22-37`) and the `Geom2dGcc_*` facades use to name the
//! branch they are in:
//!
//! | OCCT | lines | here |
//! |---|---|---|
//! | `Geom2dGcc_Type1` (`Geom2dGcc_Type1.hxx:20-28`) | 30 | [`Geom2dGccType1`] |
//! | `Geom2dGcc_Type2` (`Geom2dGcc_Type2.hxx:20-34`) | 36 | [`Geom2dGccType2`] |
//! | `Geom2dGcc_Type3` (`Geom2dGcc_Type3.hxx:20-26`) | 26 | [`Geom2dGccType3`] |
//!
//! The enumerators keep OCCT's names minus the `Geom2dGcc_` prefix, so
//! `Geom2dGcc_CuCuCu` is [`Geom2dGccType1::CuCuCu`]. `Type1` / `Type2` / `Type3`
//! are three separate C++ `enum`s whose values overlap (`Geom2dGcc_CuCu` /
//! `Geom2dGcc_CuCuOnCu` / `Geom2dGcc_CuCuCu` all start at 0); three Rust enums
//! reproduce that, where a single flat enum would have to rename them.
//!
//! UNPORTED: none.

/// `Geom2dGcc_Type1` (`Geom2dGcc_Type1.hxx:22-28`): the argument type triple
/// of `Geom2dGcc_Circ2d3Tan`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Geom2dGccType1 {
    /// `Geom2dGcc_CuCuCu`.
    CuCuCu,
    /// `Geom2dGcc_CiCuCu`.
    CiCuCu,
    /// `Geom2dGcc_CiCiCu`.
    CiCiCu,
    /// `Geom2dGcc_CiLiCu`.
    CiLiCu,
    /// `Geom2dGcc_LiLiCu`.
    LiLiCu,
    /// `Geom2dGcc_LiCuCu`.
    LiCuCu,
}

/// `Geom2dGcc_Type2` (`Geom2dGcc_Type2.hxx:22-34`): the argument type pair and
/// the carrier of `Geom2dGcc_Circ2d2TanOn`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Geom2dGccType2 {
    /// `Geom2dGcc_CuCuOnCu`.
    CuCuOnCu,
    /// `Geom2dGcc_CiCuOnCu`.
    CiCuOnCu,
    /// `Geom2dGcc_LiCuOnCu`.
    LiCuOnCu,
    /// `Geom2dGcc_CuPtOnCu`.
    CuPtOnCu,
    /// `Geom2dGcc_CuCuOnLi`.
    CuCuOnLi,
    /// `Geom2dGcc_CiCuOnLi`.
    CiCuOnLi,
    /// `Geom2dGcc_LiCuOnLi`.
    LiCuOnLi,
    /// `Geom2dGcc_CuPtOnLi`.
    CuPtOnLi,
    /// `Geom2dGcc_CuCuOnCi`.
    CuCuOnCi,
    /// `Geom2dGcc_CiCuOnCi`.
    CiCuOnCi,
    /// `Geom2dGcc_LiCuOnCi`.
    LiCuOnCi,
    /// `Geom2dGcc_CuPtOnCi`.
    CuPtOnCi,
}

/// `Geom2dGcc_Type3` (`Geom2dGcc_Type3.hxx:22-26`): the second-argument type
/// of `Geom2dGcc_Circ2d2TanRad`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Geom2dGccType3 {
    /// `Geom2dGcc_CuCu`.
    CuCu,
    /// `Geom2dGcc_CiCu`.
    CiCu,
}
