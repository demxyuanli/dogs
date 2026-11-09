//! Directory and parameter scan from `analiges.c`.

use super::liriges::c_strlen;
use super::structiges::IgesTables;

pub(in crate::iges_read) const ARG_VIDE: i32 = 0;
pub(in crate::iges_read) const ARG_QUID: i32 = 1;
pub(in crate::iges_read) const ARG_CHAR: i32 = 2;
pub(in crate::iges_read) const ARG_INT: i32 = 3;
pub(in crate::iges_read) const ARG_SIGN: i32 = 4;
pub(in crate::iges_read) const ARG_REAL: i32 = 5;
pub(in crate::iges_read) const ARG_EXP: i32 = 6;
pub(in crate::iges_read) const ARG_REXP: i32 = 7;
pub(in crate::iges_read) const ARG_MEXP: i32 = 8;

const BASES: [i32; 10] = [
    1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000,
];

/// `IGES_decode` (`analiges.c:52-64`).
pub(in crate::iges_read) fn iges_decode(line: &[u8], depuis: usize, tant: usize) -> i32 {
    let mut val = 0i32;
    let depart = depuis + tant - 1;
    for i in 0..tant {
        if depart < i || depart - i >= line.len() {
            break;
        }
        let uncar = line[depart - i];
        if uncar == b' ' {
            break;
        } else if uncar == b'+' {
            continue;
        } else if uncar == b'-' {
            val = -val;
        } else if uncar != b'0' {
            val += (uncar as i32 - 48) * BASES[i];
        }
    }
    val
}

fn copstr(line: &[u8], depuis: usize, tant: usize) -> String {
    let end = (depuis + tant).min(line.len());
    let start = depuis.min(end);
    let bytes = &line[start..end];
    let n = bytes.iter().position(|c| *c == 0).unwrap_or(bytes.len());
    String::from_utf8_lossy(&bytes[..n]).into_owned()
}

/// `iges_Dsect` (`analiges.c:75-113`).
pub(in crate::iges_read) fn iges_dsect(tables: &mut IgesTables, dstat: &mut i32, numsec: i32, line: &[u8]) {
    if *dstat == 0 {
        tables.new_part(numsec);
        if let Some(cur) = tables.current_part_mut() {
            cur.typ = iges_decode(line, 0, 8);
            cur.tra = iges_decode(line, 48, 8);
            cur.blk = iges_decode(line, 64, 2);
            let _poi = iges_decode(line, 8, 8);
            let _pdef = iges_decode(line, 16, 8);
            let _niv = iges_decode(line, 32, 8);
            let _vue = iges_decode(line, 40, 8);
            let _aff = iges_decode(line, 56, 8);
            let _sub = iges_decode(line, 66, 2);
            let _use = iges_decode(line, 68, 2);
            let _her = iges_decode(line, 70, 2);
        }
        *dstat = 1;
    } else if *dstat == 1 {
        if let Some(cur) = tables.current_part_mut() {
            cur.form = iges_decode(line, 32, 8);
            let _typ2 = iges_decode(line, 0, 8);
            let _epa = iges_decode(line, 8, 8);
            let _col = iges_decode(line, 16, 8);
            let _nbl = iges_decode(line, 24, 8);
            let _res1 = copstr(line, 40, 8);
            let _res2 = copstr(line, 48, 8);
            let _nom = copstr(line, 56, 8);
            let _num = copstr(line, 64, 8);
        }
        *dstat = 0;
    }
}

/// `iges_Psect` (`analiges.c:118-129`). `atoi` of the DE field at column 66.
pub(in crate::iges_read) fn iges_psect(tables: &mut IgesTables, line: &mut [u8; 100]) {
    let dnum = c_atoi(&line[65..]);
    line[64] = 0;
    tables.select_part(dnum);
}

fn c_atoi(bytes: &[u8]) -> i32 {
    let mut i = 0usize;
    while i < bytes.len() && (bytes[i] == b' ' || bytes[i] == b'\t') {
        i += 1;
    }
    let mut sign = 1i32;
    if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
        if bytes[i] == b'-' {
            sign = -1;
        }
        i += 1;
    }
    let mut val = 0i32;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        val = val.saturating_mul(10).saturating_add((bytes[i] - b'0') as i32);
        i += 1;
    }
    val * sign
}

pub(in crate::iges_read) struct ParamScan {
    nbcar_h: i32,
    numcar: i32,
    reste: i32,
    typarg: i32,
}

impl ParamScan {
    pub(in crate::iges_read) fn new() -> Self {
        Self {
            nbcar_h: 0,
            numcar: 0,
            reste: 0,
            typarg: ARG_VIDE,
        }
    }
}

/// `iges_param` (`analiges.c:146-261`).
pub(in crate::iges_read) fn iges_param(
    tables: &mut IgesTables,
    scan: &mut ParamScan,
    pstat: &mut i32,
    line: &[u8],
    c_separ: u8,
    c_fin: u8,
    lonlin: i32,
) {
    if *pstat == 0 {
        scan.reste = 0;
    }
    if *pstat != 2 {
        scan.numcar = 0;
    }
    if *pstat < 3 {
        scan.nbcar_h = 0;
    } else {
        scan.numcar = scan.nbcar_h;
        if scan.numcar > lonlin {
            let n = lonlin.max(0) as usize;
            let end = n.min(c_strlen(line));
            tables.add_param(&String::from_utf8_lossy(&line[..end]));
            scan.nbcar_h -= lonlin;
            return;
        } else {
            let n = scan.nbcar_h.max(0) as usize;
            let end = n.min(c_strlen(line));
            tables.add_param(&String::from_utf8_lossy(&line[..end]));
            scan.nbcar_h = 0;
        }
    }
    let mut i0 = 0usize;
    scan.typarg = ARG_VIDE;
    let mut param = [0u8; 80];
    let mut i = 0usize;
    loop {
        let at = scan.numcar as usize + i;
        if at >= line.len() || line[at] == 0 {
            break;
        }
        let unpar = line[at];
        if unpar == c_separ {
            *pstat = 2;
            param[i.min(79)] = 0;
            emit(tables, scan, &param, i0, i);
            scan.reste = 0;
            let mut j = i + 1;
            loop {
                let jat = scan.numcar as usize + j;
                if jat >= line.len() || line[jat] == 0 {
                    *pstat = 1;
                    return;
                }
                if line[jat] != b' ' {
                    scan.numcar += (i + 1) as i32;
                    return;
                }
                j += 1;
            }
        }
        if unpar == c_fin {
            *pstat = 1;
            param[i.min(79)] = 0;
            emit(tables, scan, &param, i0, i);
            scan.reste = 0;
            return;
        }
        if i < 80 {
            param[i] = unpar;
        }
        if (48..58).contains(&unpar) {
            if scan.typarg == ARG_INT {
            } else if scan.typarg == ARG_VIDE {
                scan.typarg = ARG_INT;
            } else if scan.typarg == ARG_EXP {
                scan.typarg = ARG_REXP;
            }
        } else if unpar == b'+' || unpar == b'-' {
            if scan.typarg == ARG_VIDE {
                scan.typarg = ARG_SIGN;
            } else if scan.typarg != ARG_EXP && scan.typarg != ARG_MEXP {
                scan.typarg = ARG_QUID;
            }
        } else if unpar == b'.' {
            if scan.typarg == ARG_VIDE || scan.typarg == ARG_INT || scan.typarg == ARG_SIGN {
                scan.typarg = ARG_REAL;
            } else {
                scan.typarg = ARG_QUID;
            }
        } else if unpar == b'E' || unpar == b'e' || unpar == b'D' || unpar == b'd' {
            if scan.typarg == ARG_REAL {
                scan.typarg = ARG_EXP;
            } else if scan.typarg == ARG_INT || scan.typarg == ARG_SIGN {
                scan.typarg = ARG_MEXP;
            } else {
                scan.typarg = ARG_QUID;
            }
        } else if unpar == b'H' {
            if scan.typarg != ARG_INT {
                scan.typarg = ARG_QUID;
            } else {
                scan.typarg = ARG_CHAR;
                scan.nbcar_h = 0;
                let mut j = i0;
                while j < i {
                    if (48..58).contains(&param[j]) {
                        scan.nbcar_h = scan.nbcar_h * 10 + (param[j] as i32 - 48);
                    } else {
                        scan.nbcar_h = 0;
                        break;
                    }
                    j += 1;
                }
                if scan.numcar + i as i32 + scan.nbcar_h >= lonlin {
                    let mut j = 1i32;
                    while j < lonlin - scan.numcar - i as i32 {
                        let src = (scan.numcar as usize) + i + j as usize;
                        let dst = i + j as usize;
                        if src < line.len() && dst < 80 {
                            param[dst] = line[src];
                        }
                        j += 1;
                    }
                    let end = (lonlin - scan.numcar).max(0) as usize;
                    if end < 80 {
                        param[end] = 0;
                    }
                    scan.nbcar_h = scan.numcar + i as i32 + scan.nbcar_h + 1 - lonlin;
                    *pstat = 3;
                    let text_end = end.min(80);
                    let text = String::from_utf8_lossy(&param[i0..text_end]);
                    let text = text.trim_end_matches('\0').to_string();
                    // `analiges.c:232` always starts a parameter; `reste` stays 1
                    // so the next physical line appends (`iges_addparam`).
                    tables.new_param(scan.typarg, &text);
                    scan.reste = 1;
                    return;
                } else {
                    let mut j = 1i32;
                    while j <= scan.nbcar_h {
                        let src = (scan.numcar as usize) + i + j as usize;
                        let dst = i + j as usize;
                        if src < line.len() && dst < 80 {
                            param[dst] = line[src];
                        }
                        j += 1;
                    }
                    i += scan.nbcar_h.max(0) as usize;
                }
            }
        } else if unpar == b' ' {
            if scan.typarg == ARG_VIDE {
                i0 = i + 1;
            } else {
                let mut j = i + 1;
                loop {
                    let jat = scan.numcar as usize + j;
                    if jat >= line.len() || line[jat] == 0 {
                        break;
                    }
                    let nxt = line[jat];
                    if nxt == c_separ || nxt == c_fin {
                        break;
                    }
                    if nxt != b' ' {
                        scan.typarg = ARG_QUID;
                        break;
                    }
                    j += 1;
                }
            }
        } else {
            scan.typarg = ARG_QUID;
        }
        i += 1;
        if i >= 80 {
            break;
        }
    }
    *pstat = 1;
    if i < 80 {
        param[i] = 0;
    }
    scan.reste = -1;
    if i > i0 {
        let text = String::from_utf8_lossy(&param[i0..i.min(80)]);
        tables.new_param(scan.typarg, text.trim_end_matches('\0'));
    }
}

fn emit(tables: &mut IgesTables, scan: &ParamScan, param: &[u8], i0: usize, i: usize) {
    let end = i.min(param.len());
    let start = i0.min(end);
    let text = String::from_utf8_lossy(&param[start..end]);
    let text = text.trim_end_matches('\0');
    if scan.reste == 0 {
        tables.new_param(scan.typarg, text);
    } else if scan.reste > 0 {
        tables.add_param(text);
    }
}
