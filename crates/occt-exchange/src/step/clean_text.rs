//! STEP string decoding. Source: `StepData_StepReaderData::cleanText`
//! (`StepData_StepReaderData.cxx:116-343`) and its helper
//! `convertCharacterTo16bit` (`cxx:67-112`).
//!
//! Every quoted string parameter of a STEP file goes through `cleanText` before
//! a reader sees it (`StepData_StepReaderData::ReadString` `cxx:1859`,
//! `ReadNext` `cxx:904` / `:1020` / `:1112` / `:1305`): the ISO 10303-21 control
//! directives `\X\`, `\X2\…\X0\`, `\X4\…\X0\`, `\S\`, `\P*\`, `\N\`, `\T\` and
//! the doubled `''` / `\\` escapes are decoded, line breaks inside the literal
//! are dropped, and the surrounding apostrophes are removed.
//!
//! OCCT runs the loop on UTF-16 code units of a string converted from the
//! file's *source code page* (`cxx:124-125`). The source code page defaults to
//! UTF-8 (`StepData_StepReaderData.hxx:62`), so this port works on the
//! already-decoded Rust `str`; `toConversion` (`cxx:124`) is always true here
//! and the splice branch (`cxx:276-290`) always runs.
//!
//! **UNPORTED**: the `\P{N}\` / `\S\` code page tables. `\S\` converts a single
//! ISO 8859 character through `Resource_Unicode::ConvertFormatToUnicode`
//! (`Resource_Unicode.cxx:655-706`, `THE_CODEPAGES_ANSI` for ISO 8859-1…9 and
//! CP1250…) which this crate does not carry; only the default ISO 8859-1 page
//! (`cxx:130`) is decodable, as a byte-to-code-point identity.

/// `convertCharacterTo16bit` (`cxx:67-112`): a hex digit, 0 when the character
/// is not a hexadecimal digit.
fn hex_digit(c: char) -> u32 {
    match c {
        '0'..='9' => c as u32 - '0' as u32,
        'A'..='F' => c as u32 - 'A' as u32 + 10,
        'a'..='f' => c as u32 - 'a' as u32 + 10,
        _ => 0,
    }
}

/// A UTF-16 code unit as a `char`: the port keeps OCCT's `?` for a value that
/// is not a scalar (an unpaired surrogate, which
/// `TCollection_ExtendedString::AssignCat(char16_t)` stores as-is but this
/// crate's `char` cannot).
fn from_code_unit(unit: u32) -> char {
    char::from_u32(unit).unwrap_or('?')
}

/// `StepData_StepReaderData::cleanText` (`cxx:116-343`). `quoted` includes both
/// surrounding apostrophes, as the file spells them.
pub(super) fn clean_text(quoted: &str) -> String {
    // `if (theVal->Length() == 2) { theVal->Clear(); return; }` (`cxx:118-122`):
    // the empty literal.
    if quoted.chars().count() == 2 {
        return String::new();
    }

    let mut res: Vec<char> = quoted.chars().skip(1).collect();
    // Indices stay 1-based as in the cxx; `size` is `aResStringSize`, the index
    // of the closing apostrophe (`cxx:126`). `res[0]` is the first character
    // after the opening apostrophe, which `cleanText` skips by passing
    // `theVal->ToCString() + 1` to the conversion (`cxx:125`).
    // Indices stay 1-based as in the cxx; `size` is `aResStringSize`, the index
    // of the closing apostrophe (`cxx:126`).
    let mut size = res.len() as i64 - 1;
    // `aTempExtString` (`cxx:127`): the characters a control directive decodes
    // to.
    let mut temp: Vec<char> = Vec::new();
    // `aLocalFormatType` (`cxx:130`): the code page a `\S\` directive uses,
    // ISO 8859-1 until a `\P*\` directive changes it.
    let mut local_page = 1i64;
    let mut i = 1i64;
    // `for (int aStringInd = 1; aStringInd <= aResStringSize; ++aStringInd)`
    // (`cxx:132`).
    loop {
        if i > size {
            break;
        }
        let cur = res[(i - 1) as usize];
        // `aSetCharInd = aStringInd` (`cxx:135`).
        let mut set_ind = i;

        // `if (aChar == '\\' && aStringInd <= aResStringSize - 3)`
        // (`cxx:136`).
        if cur == '\\' && i <= size - 3 {
            let dir = res[i as usize];
            let sec_slash = res[(i + 1) as usize] == '\\';
            let third_slash = res[(i + 2) as usize] == '\\';
            let mut converted = false;

            // `("\P{N}\")` (`cxx:148-160`).
            if dir == 'P' && third_slash {
                let page = res[(i + 1) as usize].to_ascii_uppercase();
                if ('A'..='I').contains(&page) {
                    // `aLocalFormatType = iso8859_1 + (aPageId - 'A')`
                    // (`cxx:152`).
                    local_page = page as i64 - 'A' as i64 + 1;
                }
                // UNPORTED: the "unsupported symbol in place of *" warning
                // (`cxx:155-159`); the port's `Resolver` collects no STEP
                // check messages.
                converted = true;
                i += 3;
            }
            // `("\S\")` (`cxx:165-172`).
            else if dir == 'S' && sec_slash {
                // `Resource_Unicode::ConvertFormatToUnicode(aLocalFormatType,
                // {aResChar | 0x80})` (`cxx:168-169`): only the ISO 8859-1 page
                // is ported (module header), where the conversion is the
                // byte-to-code-point identity. Any other page reaches OCCT's
                // `aRes == 0 -> '?'` fallback (`Resource_Unicode.cxx:701-703`)
                // because `THE_CODEPAGES_ANSI` is unported.
                let byte = res[(i + 2) as usize] as u32 | 0x80;
                temp = vec![if local_page == 1 { from_code_unit(byte) } else { '?' }];
                converted = true;
                i += 3;
            }
            // `("\X\")` (`cxx:178-188`): two hexadecimal digits, U+0000..U+00FF.
            else if dir == 'X' && i <= size - 4 && sec_slash {
                let hi = hex_digit(res[(i + 2) as usize]);
                let lo = hex_digit(res[(i + 3) as usize]);
                temp = vec![from_code_unit((hi << 4) | lo)];
                converted = true;
                i += 4;
            }
            // `("\X2\" or "\X4\")` (`cxx:193-273`).
            else if dir == 'X' && third_slash {
                let first = i + 3;
                let mut last = i;
                let mut closed = false;
                // `for (; aStringInd <= aResStringSize && !isClosed;
                // ++aStringInd)` (`cxx:199-209`): find the closing `\X0\`.
                let mut j = i;
                while j <= size {
                    if size - j > 2
                        && res[(j - 1) as usize] == '\\'
                        && res[j as usize] == 'X'
                        && res[(j + 1) as usize] == '0'
                        && res[(j + 2) as usize] == '\\'
                    {
                        last = j - 1;
                        // `aStringInd = aStringInd + 2` (`cxx:206`) plus the
                        // loop's `++aStringInd`: the index of the closing
                        // backslash of `\X0\`.
                        i = j + 3;
                        closed = true;
                        break;
                    }
                    j += 1;
                }
                if !closed {
                    // `if (!isClosed) aLastInd = aStringInd = aResStringSize`
                    // (`cxx:210-213`).
                    last = size;
                    i = size;
                }
                let str_len = last - first;

                // `"\X2\"` (`cxx:218-241`): UTF-16, groups of four or three
                // hexadecimal digits.
                if res[(first - 2) as usize] == '2' && size - first > 3 {
                    let step = if str_len % 4 == 0 { 4 } else { 3 };
                    if str_len % step != 0 {
                        // `aTempExtString.AssignCat('?')` (`cxx:223`).
                        temp = vec!['?'];
                    } else {
                        let mut unit: u32 = 0;
                        for k in 1..=str_len {
                            unit |= hex_digit(res[(k + first - 1) as usize]);
                            if k % step == 0 {
                                temp.push(from_code_unit(unit));
                                unit = 0;
                            }
                            unit <<= 4;
                        }
                    }
                }
                // `"\X4\"` (`cxx:245-272`): UTF-32, groups of eight hexadecimal
                // digits.
                else if res[(first - 2) as usize] == '4' && size - first > 7 {
                    if str_len % 8 != 0 {
                        temp = vec!['?'];
                    } else {
                        let mut unit: u32 = 0;
                        for k in 1..=str_len {
                            unit |= hex_digit(res[(k + first - 1) as usize]);
                            if k % 8 == 0 {
                                temp.push(from_code_unit(unit));
                                unit = 0;
                            }
                            unit <<= 4;
                        }
                    }
                }
                converted = true;
            }

            // `if (isConverted)` (`cxx:276-291`): replace the directive — the
            // characters `aSetCharInd..=aStringInd` — with `aTempExtString`.
            if converted {
                // `aResStringSize -= aStringInd - aSetCharInd -
                // aTempExtString.Length() + 1` (`cxx:280`).
                size -= i - set_ind - temp.len() as i64 + 1;
                let head = res[..(set_ind - 1) as usize].to_vec();
                let tail: Vec<char> = res[i as usize..].to_vec();
                res = head;
                res.extend(temp.iter().copied());
                res.extend(tail);
                set_ind += temp.len() as i64;
                res.truncate(size.max(0) as usize);
                // `aStringInd = aSetCharInd - 1` (`cxx:285`), then the loop's
                // `++aStringInd`.
                i = set_ind;
                temp.clear();
                continue;
            }
        }

        // `if (aStringInd <= aResStringSize - 1)` (`cxx:293`).
        if i <= size - 1 {
            let next = res[i as usize];
            // `if (aCharNext == aChar && (aChar == '\'' || aChar == '\\'))`
            // (`cxx:297-302`): a doubled apostrophe or backslash collapses into
            // one.
            if next == cur && (cur == '\'' || cur == '\\') {
                res.remove((i - 1) as usize);
                size -= 1;
                res.truncate(size.max(0) as usize);
            } else if cur == '\\' {
                // `const bool isDirective = ...` (`cxx:305-308`).
                let is_directive = i <= size - 2 && res[(i + 1) as usize] == '\\';
                if is_directive && next == 'N' {
                    // `aResString.SetValue(aSetCharInd++, '\n')` (`cxx:310-319`).
                    res[(i - 1) as usize] = '\n';
                    res.drain(i as usize..(i + 2) as usize);
                    size -= 2;
                    res.truncate(size.max(0) as usize);
                    i += 1;
                    continue;
                } else if is_directive && next == 'T' {
                    res[(i - 1) as usize] = '\t';
                    res.drain(i as usize..(i + 2) as usize);
                    size -= 2;
                    res.truncate(size.max(0) as usize);
                    i += 1;
                    continue;
                }
            }
        }

        // `if (aChar == '\n' || aChar == '\r')` (`cxx:331-337`): a line break
        // inside the literal is dropped.
        if cur == '\n' || cur == '\r' {
            res.remove((i - 1) as usize);
            size -= 1;
            res.truncate(size.max(0) as usize);
            // `aStringInd--` (`cxx:336`), then the loop's `++aStringInd`.
        }
        i += 1;
    }

    // `aResString.Trunc(aResStringSize)` (`cxx:340`) drops the closing
    // apostrophe; `TCollection_AsciiString(aResString, 0)` (`cxx:341-342`) emits
    // it as UTF-8, which is this crate's `String`.
    res.truncate(size.max(0) as usize);
    res.into_iter().collect()
}
