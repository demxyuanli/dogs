//! Transformation form enum. Source: `gp_TrsfForm.hxx`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrsfForm { Identity, Rotation, Translation, PntMirror, Ax1Mirror, Ax2Mirror, Scale, CompoundTrsf, Other }
