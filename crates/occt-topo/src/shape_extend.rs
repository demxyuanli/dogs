//! Port of the OCCT `ShapeExtend` status codec (`ShapeExtend.hxx:52-56`,
//! `ShapeExtend.cxx:54-109`). The healing algorithms accumulate a bit flag and
//! query it through `ShapeExtend::DecodeStatus(flag, status)`.

/// `ShapeExtend::EncodeStatus` (`ShapeExtend.cxx:54-98`).
pub const SHAPEEXTEND_OK: i32 = 0x0000;
pub const SHAPEEXTEND_DONE1: i32 = 0x0001;
pub const SHAPEEXTEND_DONE2: i32 = 0x0002;
pub const SHAPEEXTEND_DONE3: i32 = 0x0004;
pub const SHAPEEXTEND_DONE4: i32 = 0x0008;
pub const SHAPEEXTEND_DONE5: i32 = 0x0010;
pub const SHAPEEXTEND_DONE6: i32 = 0x0020;
pub const SHAPEEXTEND_DONE7: i32 = 0x0040;
pub const SHAPEEXTEND_DONE8: i32 = 0x0080;
/// `ShapeExtend_DONE` = any of DONE1..DONE8.
pub const SHAPEEXTEND_DONE: i32 = 0x00ff;
pub const SHAPEEXTEND_FAIL1: i32 = 0x0100;
pub const SHAPEEXTEND_FAIL2: i32 = 0x0200;
pub const SHAPEEXTEND_FAIL3: i32 = 0x0400;
pub const SHAPEEXTEND_FAIL4: i32 = 0x0800;
pub const SHAPEEXTEND_FAIL5: i32 = 0x1000;
pub const SHAPEEXTEND_FAIL6: i32 = 0x2000;
pub const SHAPEEXTEND_FAIL7: i32 = 0x4000;
pub const SHAPEEXTEND_FAIL8: i32 = 0x8000;
/// `ShapeExtend_FAIL` = any of FAIL1..FAIL8.
pub const SHAPEEXTEND_FAIL: i32 = 0xff00;

/// `ShapeExtend::DecodeStatus` (`ShapeExtend.cxx:102-109`): `OK` queries the
/// empty flag, every other status queries its encoded bit.
pub fn decode_status(flag: i32, status: i32) -> bool {
    if status == SHAPEEXTEND_OK {
        return flag == 0;
    }
    (flag & status) != 0
}
