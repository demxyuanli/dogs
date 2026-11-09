//! TColStd — standard collection aliases. Source: `TColStd/`
//! Maps OCCT's TColStd_Array1Of* → Rust Vec<T>.
//!
//! OCCT naming convention: TColStd_Array1OfReal → Vec<f64>
//!                          TColStd_Array1OfInteger → Vec<i32>
//!                          TColStd_SequenceOfReal → Vec<f64>
//!                          TColStd_HArray1OfReal → Arc<Vec<f64>>
//!
//! All 1-indexed in OCCT → 0-indexed in Rust. Subtract 1 from indices.

use std::sync::Arc;

// Sequence aliases (dynamically sized, ordered)
pub type Array1OfReal = Vec<f64>;
pub type Array1OfInteger = Vec<i32>;
pub type Array1OfBoolean = Vec<bool>;
pub type Array1OfByte = Vec<u8>;
pub type Array1OfCharacter = Vec<char>;
pub type Array1OfTransient = Vec<Arc<()>>;

pub type SequenceOfReal = Vec<f64>;
pub type SequenceOfInteger = Vec<i32>;
pub type SequenceOfAddress = Vec<usize>;

// Handle-based arrays (shared ownership)
pub type HArray1OfReal = Arc<Vec<f64>>;
pub type HArray1OfInteger = Arc<Vec<i32>>;
pub type HArray1OfBoolean = Arc<Vec<bool>>;

// Packed maps (bitset)
pub type PackedMapOfInteger = crate::kernel::containers::PackedMap;

// Data maps
pub type DataMapOfIntegerReal = std::collections::HashMap<i32, f64>;
pub type DataMapOfIntegerInteger = std::collections::HashMap<i32, i32>;
pub type DataMapOfStringInteger = std::collections::HashMap<String, i32>;
pub type IndexedDataMapOfTransientTransient = crate::kernel::containers::IndexedDataMap<usize, usize>;

// List
pub type ListOfInteger = Vec<i32>;
pub type ListOfReal = Vec<f64>;
pub type ListOfTransient = Vec<Arc<()>>;
