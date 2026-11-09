//! Fixed-width IGES cards: the 80-column record layout and the number/Hollerith formatters. Source: `IGESData_IGESWriter` (`IGESData_IGESWriter.cxx:763-880`) with `Interface_FloatWriter` and `Interface_LineBuffer`.
use super::prelude::*;

// ---------------------------------------------------------------------------
// Formatting helpers
// ---------------------------------------------------------------------------

/// Pad or truncate `s` to exactly 80 columns.
pub fn rec80(s: &str) -> String {
    let mut line = String::with_capacity(80);
    line.push_str(s);
    if line.len() > 80 {
        line.truncate(80);
    } else {
        line.push_str(&" ".repeat(80 - line.len()));
    }
    line
}

/// Data columns of a Start/Global/Directory card
/// (`IGESData_IGESWriter.cxx:38-39`: `MaxcarsG = 72`).
pub(super) const MAXCARS_G: usize = 72;
/// Data columns of a Parameter card (`cxx:39`: `MaxcarsP = 64`).
pub(super) const MAXCARS_P: usize = 64;

/// One IGES card: the data field left-justified in `width` columns, then
/// `section` and the card's sequence number in the last seven columns
/// (`IGESData_IGESWriter.cxx:836-880` for D, `:769-794` for G/S).
/// The sequence is `%7.7d` (`Sprintf(finlin, "G%7.7d", i)` and `D%7.7d`),
/// zero-filled, not blank-filled.
pub(super) fn sec_line(section: char, seq: usize, content: &str, width: usize) -> String {
    let mut s = String::with_capacity(80);
    for ch in content.chars().take(width) {
        s.push(ch);
    }
    s.push_str(&" ".repeat(width - s.len()));
    s.push(section);
    s.push_str(&format!("{seq:07}"));
    s
}

/// One Parameter card (`IGESData_IGESWriter.cxx:902-925`):
/// `data` (64 columns) + blank + the owning entity's directory-entry pointer
/// (`2*i - 1`) + `P` + the card's sequence number.
/// `Sprintf(finlin, " %7.7dP%7.7d", 2*i-1, j)` (`cxx:903`).
pub(super) fn param_line(content: &str, de_pointer: usize, seq: usize) -> String {
    let mut s = String::with_capacity(80);
    for ch in content.chars().take(MAXCARS_P) {
        s.push(ch);
    }
    s.push_str(&" ".repeat(MAXCARS_P - s.len()));
    s.push_str(&format!(" {de_pointer:07}P{seq:07}"));
    s
}

/// The `Interface_LineBuffer` a section writer fills
/// (`IGESData_IGESWriter.cxx:49-50`: `thecurr(MaxcarsG + 1)`, `:248`:
/// `thecurr.SetMax(MaxcarsP)` at the start of the parameter section): text is
/// appended while it fits on the current card, and a card that overflows is
/// dumped **before** the text that does not fit is written
/// (`Interface_LineBuffer.cxx:68-80` `CanGet`, `:150-157` `Keep`).
///
/// This is what decides where a Parameter card ends: not a fixed 64-column cut,
/// but the first field boundary that does not fit. `Send(v)`
/// (`IGESData_IGESWriter.cxx:560-565`, `:586-589`) writes the separator and the
/// value as two separate additions, so a card can also end on a separator whose
/// value went to the next card.
pub(super) struct CardBuffer {
    /// Cards already dumped, i.e. `thepars` of `IGESData_IGESWriter`.
    pub(super) cards: Vec<String>,
    /// The card being filled, i.e. `thecurr`.
    pub(super) line: String,
    /// `thecurr`'s maximum length: `MaxcarsP` (64) in the parameter section,
    /// `MaxcarsG` (72) in the global one (`IGESData_IGESWriter.cxx:39-40`).
    pub(super) max: usize,
}

impl CardBuffer {
    pub(super) fn new(max: usize) -> Self {
        Self {
            cards: Vec::new(),
            line: String::new(),
            max,
        }
    }

    /// `Interface_LineBuffer::Moved` through `IGESData_IGESWriter`'s flush
    /// (`IGESData_IGESWriter.cxx:508-512`, `:541-550`): the current card is
    /// appended to the section **as is** — trailing blanks are added at printing
    /// time (`cxx:913`), and a flush of a short card keeps it short.
    pub(super) fn flush(&mut self) {
        self.cards.push(std::mem::take(&mut self.line));
    }

    /// `IGESData_IGESWriter::AddString` (`cxx:497-533`). The `+ 1` of the
    /// `CanGet` call is OCCT's own reservation (comment at `cxx:504`): one
    /// character must remain free, so that the *next* separator is never the
    /// first character of a line.
    pub(super) fn add_string(&mut self, text: &str) {
        let bytes = text.as_bytes();
        if self.line.len() + bytes.len() + 1 > self.max {
            self.flush();
        }
        // Text longer than a card fills whole cards (`cxx:516-531`); the buffer
        // is empty here, so each iteration dumps exactly `max` characters.
        let mut rest = bytes;
        while rest.len() > self.max {
            let n = self.max - self.line.len();
            self.line
                .push_str(&String::from_utf8_lossy(&rest[..n.min(rest.len())]));
            self.flush();
            rest = &rest[self.max..];
        }
        let n = (self.max - self.line.len()).min(rest.len());
        self.line.push_str(&String::from_utf8_lossy(&rest[..n]));
    }

    /// `IGESData_IGESWriter::AddChar` (`cxx:535-553`), used for the parameter
    /// separators and the record delimiter.
    pub(super) fn add_char(&mut self, c: u8) {
        if self.line.len() + 1 > self.max {
            self.flush();
        }
        self.line.push(c as char);
    }
}

/// Right-justify a value in an 8-column IGES field.
pub(super) fn field8(v: &str) -> String {
    format!("{v:>8}")
}

/// Every real parameter of every entity goes through
/// `IGESData_IGESWriter::Send(const double)` (`IGESData_IGESWriter.cxx:582-589`),
/// which is `thefloatw.Write(val, lval)`, and the writer constructs `thefloatw`
/// as `Interface_FloatWriter(9)` (`IGESData_IGESWriter.cxx:49-50`, copy `:67`).
/// `Interface_FloatWriter(9)` runs `SetDefaults(9)`
/// (`Interface_FloatWriter.cxx:23-27` -> `:50-70`), which fixes the two `Sprintf`
/// formats and the range they are chosen by:
///
/// - `"%11.9f"` while `|val|` is in `[0.1, 1000)`, `"%11.9E"` outside it
///   (`cxx:139-147`);
/// - `thezerosup = true`, so `Convert` then strips trailing zeros from the
///   mantissa and drops an `E+00` exponent (`cxx:148-174`).
///
/// `%11.9E`, unlike Rust's `{:.9E}`, always writes the exponent sign and at
/// least two exponent digits (`1.5E+3` against `1.500000000E+03`), and the zero
/// suppression then leaves the shortest form (`150.`, `-0.`, `2.842170943E-14`).
pub(super) const FLOAT_RANGE_MIN: f64 = 0.1;
pub(super) const FLOAT_RANGE_MAX: f64 = 1000.0;

/// `Sprintf(text, "%11.9f", val)` (`cxx:143`). The field width is a minimum and
/// `"0.000000000"` already fills 11 columns, so no padding can ever appear.
pub(super) fn float_fixed(x: f64) -> String {
    format!("{x:.9}")
}

/// `Sprintf(text, "%11.9E", val)` (`cxx:145`), re-normalised to C's `%E`
/// layout: sign always present, exponent padded to two digits.
pub(super) fn float_sci(x: f64) -> String {
    let s = format!("{x:.9E}");
    let (mantissa, exp) = match s.split_once('E') {
        Some((m, e)) => (m, e),
        None => (s.as_str(), "0"),
    };
    let (sign, digits) = match exp.strip_prefix('-') {
        Some(d) => ('-', d),
        None => ('+', exp),
    };
    let pad = if digits.len() < 2 { "0" } else { "" };
    format!("{mantissa}E{sign}{pad}{digits}")
}

/// `Interface_FloatWriter::Convert`'s zero suppression (`cxx:148-174`): scan the
/// first 16 characters for the exponent marker, remember the exponent, cut the
/// mantissa there, drop its trailing zeros, and put the exponent back unless it
/// was exactly `+00`. Without a marker the cut point is the terminator, so a
/// fixed-format value only loses its trailing zeros (`150.000000000` -> `150.`).
pub(super) fn float_zero_suppress(text: &str) -> String {
    let chars: Vec<char> = text.chars().collect();
    let mut cut = chars.len();
    let mut exp = String::new();
    let mut i = 0;
    while i < 16 && i < chars.len() {
        if chars[i] == 'e' || chars[i] == 'E' {
            let mut e = String::from("E");
            e.extend(chars[i + 1..(i + 5).min(chars.len())].iter());
            // `cxx:159-162`: only the exponent 0 loses its marker
            // (`lxp[1] == '+' && lxp[2] == '0' && lxp[3] == '0' && lxp[4] == '\0'`).
            if e != "E+00" {
                exp = e;
            }
            cut = i;
            break;
        }
        i += 1;
    }
    let mut mantissa: String = chars[..cut].iter().collect();
    while mantissa.ends_with('0') {
        mantissa.pop();
    }
    mantissa + &exp
}

/// Render an f64 in IGES parameter syntax (always carries a decimal point or
/// exponent marker).
pub(super) fn num(x: f64) -> String {
    // `Sprintf` of a non-finite value yields `inf`/`nan`, which is not an IGES
    // number; no OCCT transfer produces one, and the port keeps the file
    // well-formed by writing the zero both formats would zero-suppress to.
    if !x.is_finite() {
        return "0.".into();
    }
    let in_range = (x >= FLOAT_RANGE_MIN && x < FLOAT_RANGE_MAX)
        || (x <= -FLOAT_RANGE_MIN && x > -FLOAT_RANGE_MAX);
    let text = if in_range { float_fixed(x) } else { float_sci(x) };
    float_zero_suppress(&text)
}

/// `nH` text parameter (`IGESData_GlobalSection::MakeHollerith`,
/// `IGESData_GlobalSection.cxx:47-69`). An empty string is an empty field.
pub(super) fn hollerith(text: &str) -> String {
    if text.is_empty() {
        return String::new();
    }
    format!("{}H{text}", text.len())
}

/// `IGESData_GlobalSection::NewDateString(..., year < 0)` (`IGESData.cxx:247`):
/// `YYYYMMDD.HHMMSS`, always 15 characters.
pub(super) fn iges_date_now() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0) as i64;
    let tod = secs.rem_euclid(86400);
    let days = secs.div_euclid(86400);
    let z = days + 719468;
    let era = (if z >= 0 { z } else { z - 146096 }) / 146097;
    let doe = (z - era * 146097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let mut y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = ((5 * doy + 2) / 153) as i64;
    let d = doy as i64 - (153 * mp + 2) / 5 + 1;
    let m = mp + if mp < 10 { 3 } else { -9 };
    if m <= 2 {
        y += 1;
    }
    let hh = tod / 3600;
    let mm = (tod % 3600) / 60;
    let ss = tod % 60;
    format!("{y:04}{m:02}{d:02}.{hh:02}{mm:02}{ss:02}")
}

/// Pointer placeholder for the `k`-th entry of an entity's reference list.
///
/// `#` never appears in IGES parameter data (strings use the `nH…` form), so a
/// `#k` token unambiguously marks a DE pointer position until
/// [`resolve_placeholders`] turns it into the final number.
pub(super) fn ph(k: usize) -> String {
    format!("#{k}")
}

/// `n` comma-separated placeholders: `#0,#1,…`.
pub(super) fn ph_list(n: usize) -> String {
    (0..n).map(ph).collect::<Vec<_>>().join(",")
}

/// True when every `#k` in `params` addresses an entry of the reference list.
pub(super) fn pointers_all_known(params: &str, n_refs: usize) -> bool {
    placeholder_indices(params).all(|k| k < n_refs)
}

/// The `k` of every `#k` token, in text order.
pub(super) fn placeholder_indices(params: &str) -> impl Iterator<Item = usize> + '_ {
    let bytes = params.as_bytes();
    let mut i = 0usize;
    let mut found: Vec<usize> = Vec::new();
    while i < bytes.len() {
        if bytes[i] == b'#' {
            let start = i + 1;
            let mut j = start;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            if j > start {
                if let Ok(k) = params[start..j].parse::<usize>() {
                    found.push(k);
                }
                i = j;
                continue;
            }
        }
        i += 1;
    }
    found.into_iter()
}

/// Replace every `#k` with the final **pointer** of the `k`-th referenced entity,
/// i.e. `Interface_InterfaceModel::DNum` = `2 * Number(entity) - 1`
/// (`Interface_InterfaceModel.cxx`), which is the sequence number of the
/// entity's **first Directory card**. `IGESData_IGESWriter::Send(handle, …)`
/// writes through `DNum` (`IGESData_IGESWriter.cxx:609-621`), so the D section
/// lays out two cards per entity and every in-model reference names an odd
/// D-section line. `refs` holds the referenced entities in parameter order and
/// `new_of` the reachability renumbering; a `#k` without a target keeps its
/// literal text (only reachable when an emitter records fewer refs than it
/// writes placeholders, which `emit_refs` already checks in debug builds).
pub(super) fn resolve_placeholders(params: &str, refs: &[usize], new_of: &HashMap<usize, usize>) -> String {
    let bytes = params.as_bytes();
    let mut out = String::with_capacity(params.len());
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes[i] == b'#' {
            let start = i + 1;
            let mut j = start;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            if j > start {
                if let Ok(k) = params[start..j].parse::<usize>() {
                    if let Some(n) = refs.get(k).and_then(|r| new_of.get(r)) {
                        out.push_str(&(2 * n - 1).to_string());
                        i = j;
                        continue;
                    }
                }
            }
        }
        out.push(bytes[i] as char);
        i += 1;
    }
    out
}

