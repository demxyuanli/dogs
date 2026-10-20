//! `IGESData_IGESModel` directory and own parameters for the planar family.
//!
//! Field order follows each tool's `ReadOwnParams`.

mod entities;
mod model;
mod parse;
mod pointers;
mod units;

pub(in crate::iges_read) use self::entities::*;
pub(in crate::iges_read) use self::model::*;
pub(in crate::iges_read) use self::parse::*;
pub(in crate::iges_read) use self::pointers::*;
pub(in crate::iges_read) use self::units::*;
