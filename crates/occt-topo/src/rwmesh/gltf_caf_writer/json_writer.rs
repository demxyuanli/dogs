//! `RWGltf_GltfOStreamWriter` — the `rapidjson::Writer` the glTF exporter
//! serializes into.
//!
//! Source: `RWGltf_GltfOStreamWriter.hxx:25-33`, `RapidJSON`'s
//! `Writer`/`internal::dtoa` (the writer is a plain `rapidjson::Writer`, i.e.
//! the non-pretty one: no whitespace anywhere).
//!
//! Two behaviours matter for byte parity with OCCT's output:
//!
//! * **Separators.** `Writer::Prefix` emits `,` between the members of one
//!   array/object and never after a key, which is what [`GltfJsonWriter`]
//!   reproduces with a `(count, has_key)` frame per nesting level.
//! * **Numbers.** `Writer::Double` goes through `internal::dtoa`, whose
//!   `Prettify` step uses ECMAScript `Number::toString` layout: plain decimal
//!   for `1e-7 <= |v| < 1e21`, exponent notation outside, and a `.0` appended
//!   when the mantissa has no fractional part (`0.0`, `1.0`). [`format_double`]
//!   takes the shortest round-trip digits (`format!("{:e}")`) and lays them
//!   out the same way. `Writer::Int`/`Int64` print as plain integers and
//!   `Writer::Bool` as `true`/`false`.
//!
//! `Writer::String` writes UTF-8 through unchanged (the writer has
//! `kWriteDefaultFlags`, so no `\u` escaping of non-ASCII) and escapes only
//! `"`, `\` and the control characters, using the short forms where JSON has
//! them.

/// The nesting level state `rapidjson::Writer` keeps in `level_stack_`.
struct Level {
    /// `Level::valueCount` — members written so far, keys included.
    count: usize,
    /// `Writer::hasKey_`: the next value belongs to the key just written.
    has_key: bool,
}

/// Output buffer state of the glTF JSON document.
pub(crate) struct GltfJsonWriter {
    out: Vec<u8>,
    stack: Vec<Level>,
}

impl GltfJsonWriter {
    /// Fresh writer at the root.
    pub(crate) fn new() -> Self {
        Self { out: Vec::new(), stack: Vec::new() }
    }

    /// The serialized bytes.
    pub(crate) fn into_bytes(self) -> Vec<u8> {
        self.out
    }

    /// `Writer::EndValue` bookkeeping for the next member of the current level.
    fn prefix(&mut self) {
        let Some(level) = self.stack.last_mut() else {
            return;
        };
        if level.has_key {
            level.has_key = false;
        } else {
            if level.count > 0 {
                self.out.push(b',');
            }
            level.count += 1;
        }
    }

    /// `Writer::StartObject`.
    pub(crate) fn start_object(&mut self) {
        self.prefix();
        self.out.push(b'{');
        self.stack.push(Level { count: 0, has_key: false });
    }

    /// `Writer::EndObject`.
    pub(crate) fn end_object(&mut self) {
        self.out.push(b'}');
        self.stack.pop();
    }

    /// `Writer::StartArray`.
    pub(crate) fn start_array(&mut self) {
        self.prefix();
        self.out.push(b'[');
        self.stack.push(Level { count: 0, has_key: false });
    }

    /// `Writer::EndArray`.
    pub(crate) fn end_array(&mut self) {
        self.out.push(b']');
        self.stack.pop();
    }

    /// `Writer::Key`.
    pub(crate) fn key(&mut self, name: &str) {
        self.prefix();
        self.write_string(name);
        // `Writer::Key` writes the `Separator(':')` right after the name and
        // before marking `hasKey_`.
        self.out.push(b':');
        if let Some(level) = self.stack.last_mut() {
            level.has_key = true;
        }
    }

    /// `Writer::String`.
    pub(crate) fn string(&mut self, value: &str) {
        self.prefix();
        self.write_string(value);
    }

    /// Write `"<value>"` with the JSON escapes `Writer::WriteString` applies.
    fn write_string(&mut self, value: &str) {
        self.out.push(b'"');
        for ch in value.chars() {
            match ch {
                '"' => self.out.extend_from_slice(b"\\\""),
                '\\' => self.out.extend_from_slice(b"\\\\"),
                '\u{08}' => self.out.extend_from_slice(b"\\b"),
                '\u{0C}' => self.out.extend_from_slice(b"\\f"),
                '\n' => self.out.extend_from_slice(b"\\n"),
                '\r' => self.out.extend_from_slice(b"\\r"),
                '\t' => self.out.extend_from_slice(b"\\t"),
                c if (c as u32) < 0x20 => {
                    // `Writer::WriteString` -> `\\u00XX` in upper-case hex.
                    self.out.extend_from_slice(format!("\\u{:04X}", c as u32).as_bytes());
                }
                c => {
                    let mut buf = [0u8; 4];
                    self.out.extend_from_slice(c.encode_utf8(&mut buf).as_bytes());
                }
            }
        }
        self.out.push(b'"');
    }

    /// `Writer::Int`.
    pub(crate) fn int(&mut self, value: i32) {
        self.prefix();
        self.out.extend_from_slice(itoa(value as i64).as_bytes());
    }

    /// `Writer::Int64`.
    pub(crate) fn int64(&mut self, value: i64) {
        self.prefix();
        self.out.extend_from_slice(itoa(value).as_bytes());
    }

    /// `Writer::Double`.
    pub(crate) fn double(&mut self, value: f64) {
        self.prefix();
        self.out.extend_from_slice(format_double(value).as_bytes());
    }

    /// `Writer::Bool`.
    pub(crate) fn boolean(&mut self, value: bool) {
        self.prefix();
        self.out.extend_from_slice(if value { b"true" } else { b"false" });
    }
}

/// Plain decimal integer formatting.
fn itoa(value: i64) -> String {
    value.to_string()
}

/// `internal::dtoa` + `Prettify` from RapidJSON (`rapidjson/internal/dtoa.h`),
/// i.e. the `Number::toString` layout the glTF writer emits.
///
/// A non-finite value cannot reach this function: every number the exporter
/// writes comes from geometry or a material, and `Writer::Double` itself would
/// print `null` for it — a branch the exporter never exercises.
pub(crate) fn format_double(value: f64) -> String {
    if value == 0.0 {
        return if value.is_sign_negative() { "-0.0" } else { "0.0" }.to_string();
    }
    let negative = value < 0.0;
    let magnitude = value.abs();
    // Shortest round-trip digits: `{:e}` yields `D.DDDDe[+-]X`.
    let sci = format!("{magnitude:e}");
    let (mantissa, exponent) = sci
        .split_once('e')
        .map(|(m, e)| (m, e.parse::<i32>().unwrap_or(0)))
        .unwrap_or((sci.as_str(), 0));
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    let length = digits.len() as i32;
    // `value = digits * 10^k`, `kk = length + k` (`Prettify`'s exponent).
    let k = exponent - (length - 1);
    let kk = length + k;
    let digits = digits.as_str();

    let body = if (0..=21).contains(&k) && kk <= 21 {
        // 1234e7 -> 12340000000.0
        let mut text = String::with_capacity(kk as usize + 2);
        text.push_str(digits);
        for _ in length..kk {
            text.push('0');
        }
        text.push_str(".0");
        text
    } else if 0 < kk && kk <= 21 {
        // 1234e-2 -> 12.34
        let split = kk as usize;
        format!("{}.{}", &digits[..split], &digits[split..])
    } else if -6 < kk && kk <= 0 {
        // 1234e-6 -> 0.001234
        let mut text = String::with_capacity(digits.len() + 2 + (-kk) as usize);
        text.push_str("0.");
        for _ in 0..(-kk) {
            text.push('0');
        }
        text.push_str(digits);
        text
    } else if length == 1 {
        // 1e30
        format!("{digits}e{}", kk - 1)
    } else {
        // 1.234e30
        format!("{}.{}e{}", &digits[..1], &digits[1..], kk - 1)
    };

    if negative {
        format!("-{body}")
    } else {
        body
    }
}
