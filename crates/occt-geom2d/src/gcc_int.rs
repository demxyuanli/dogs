//! `GccInt` package (`TKGeomAlgo/GccInt`): elementary bisecting loci between
//! two simple 2D geometric objects.
//!
//! | OCCT | lines |
//! |---|---|
//! | `GccInt_IType` | `GccInt_IType.hxx:16-24` |
//! | `GccInt_Bisec` | `GccInt_Bisec.hxx:44-80`, `GccInt_Bisec.cxx:1-57` |
//! | `GccInt_BLine` | `GccInt_BLine.cxx:1-36` |
//! | `GccInt_BCirc` | `GccInt_BCirc.cxx:1-36` |
//! | `GccInt_BParab` | `GccInt_BParab.cxx:1-36` |
//! | `GccInt_BHyper` | `GccInt_BHyper.cxx:1-36` |
//! | `GccInt_BElips` | `GccInt_BElips.cxx:1-36` |
//! | `GccInt_BPoint` | `GccInt_BPoint.cxx:1-36` |
//!
//! OCCT models the family as a deferred `GccInt_Bisec` root plus one
//! `Standard_Transient` subclass per arc type. The four accessors of the root
//! (`Point`, `Line`, `Circle`, `Hyperbola`, `Parabola`, `Ellipse`) throw
//! `Standard_DomainError` when the arc type does not match; the Rust enum
//! answers `Option` instead, so a wrong query is a value rather than a panic.
//! This mirrors the `Geom2d_Curve::gp_*2d` accessors already used across the
//! crate (`curve.rs:149-169`).

use occt_core::gp::{GpCirc2d, GpElips2d, GpHypr2d, GpLin2d, GpParab2d, GpPnt2d};

/// `enum GccInt_IType` (`GccInt_IType.hxx:16-24`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GccIntIType {
    /// `GccInt_Lin`
    Lin,
    /// `GccInt_Cir`
    Cir,
    /// `GccInt_Ell`
    Ell,
    /// `GccInt_Par`
    Par,
    /// `GccInt_Hpr`
    Hpr,
    /// `GccInt_Pnt`
    Pnt,
}

/// `GccInt_Bisec` and its six concrete arcs (`GccInt_B{Line,Circ,Parab,Hyper,Elips,Point}`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GccIntBisec {
    /// `GccInt_BLine`
    Line(GpLin2d),
    /// `GccInt_BCirc`
    Circle(GpCirc2d),
    /// `GccInt_BElips`
    Ellipse(GpElips2d),
    /// `GccInt_BParab`
    Parabola(GpParab2d),
    /// `GccInt_BHyper`
    Hyperbola(GpHypr2d),
    /// `GccInt_BPoint`
    Point(GpPnt2d),
}

impl GccIntBisec {
    /// `GccInt_Bisec::ArcType` (`GccInt_BLine.cxx:24-27` and siblings).
    pub fn arc_type(&self) -> GccIntIType {
        match self {
            Self::Line(_) => GccIntIType::Lin,
            Self::Circle(_) => GccIntIType::Cir,
            Self::Ellipse(_) => GccIntIType::Ell,
            Self::Parabola(_) => GccIntIType::Par,
            Self::Hyperbola(_) => GccIntIType::Hpr,
            Self::Point(_) => GccIntIType::Pnt,
        }
    }

    /// `GccInt_Bisec::Point` (`GccInt_Bisec.cxx:29-39`); `None` where OCCT
    /// raises `Standard_DomainError`.
    pub fn point(&self) -> Option<GpPnt2d> {
        match self {
            Self::Point(p) => Some(*p),
            _ => None,
        }
    }

    /// `GccInt_Bisec::Line` (`GccInt_Bisec.cxx:41-51`).
    pub fn line(&self) -> Option<GpLin2d> {
        match self {
            Self::Line(l) => Some(*l),
            _ => None,
        }
    }

    /// `GccInt_Bisec::Circle` (`GccInt_Bisec.cxx:53-57`).
    pub fn circle(&self) -> Option<GpCirc2d> {
        match self {
            Self::Circle(c) => Some(*c),
            _ => None,
        }
    }

    /// `GccInt_Bisec::Hyperbola` (`GccInt_Bisec.cxx`).
    pub fn hyperbola(&self) -> Option<GpHypr2d> {
        match self {
            Self::Hyperbola(h) => Some(*h),
            _ => None,
        }
    }

    /// `GccInt_Bisec::Parabola` (`GccInt_Bisec.cxx`).
    pub fn parabola(&self) -> Option<GpParab2d> {
        match self {
            Self::Parabola(p) => Some(*p),
            _ => None,
        }
    }

    /// `GccInt_Bisec::Ellipse` (`GccInt_Bisec.cxx`).
    pub fn ellipse(&self) -> Option<GpElips2d> {
        match self {
            Self::Ellipse(e) => Some(*e),
            _ => None,
        }
    }
}
