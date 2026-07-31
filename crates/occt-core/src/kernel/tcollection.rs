//! TCollection — typed collection aliases. Source: `TCollection/`
//! Maps OCCT's TCollection_* to Rust stdlib.
//!
//! OCCT naming: TCollection_AsciiString → Rust String
//!              TCollection_ExtendedString → Rust String (UTF-8)
//!              TCollection_HExtendedString → Box<[char16_t]>
//!
//! Ponteil note: Rust String is UTF-8, covers both OCCT string types. No HExtendedString needed.

pub type AsciiString = String;
pub type ExtendedString = String;
pub type HAsciiString = String;   // Handle(TCollection_HAsciiString) → String
pub type HExtendedString = String;
