//! `igesread` (`igesread.c:51-130`).
//!
//! `IGESFile_Read` calls this with `modefnes` false (`IGESFile_Read.cxx:94`).

use std::collections::HashMap;

use super::analiges::{
    iges_dsect, iges_param, iges_psect, ParamScan, ARG_CHAR, ARG_EXP, ARG_INT, ARG_MEXP, ARG_QUID,
    ARG_REAL, ARG_REXP, ARG_SIGN, ARG_VIDE,
};
use super::iges_data_iges_model::DirLine;
use super::iges_data_param_reader::IgesParam;
use super::liriges::{c_strlen, LineReader};
use super::structiges::{IgesTables, RawParam};

pub(in crate::iges_read) struct IgesRead {
    pub global: Vec<IgesParam>,
    pub directory: Vec<DirLine>,
    pub params_by_de: HashMap<i32, Vec<IgesParam>>,
}

pub(in crate::iges_read) fn igesread(data: &[u8]) -> Result<IgesRead, String> {
    let mut reader = LineReader::new(data.to_vec());
    let mut tables = IgesTables::new();
    let mut scan = ParamScan::new();
    let mut line = [0u8; 100];
    let mut numsec = 0i32;
    let mut lesect = [0i32; 6];
    let mut i0 = 0i32;
    let mut dstat = 0i32;
    let mut pstat = 0i32;
    let mut c_separ = b',';
    let mut c_fin = b';';
    let modefnes = false;

    loop {
        let i = reader.lire(&mut numsec, &mut line, modefnes);
        if i <= 0 || i < i0 {
            if i == 0 {
                break;
            }
            // `XSTEP_18` (`igesread.c:75`). Before any section the file is rejected.
            if i0 == 0 {
                return Err("IGES syntax error before the first section: igesread.c:78-81".into());
            }
            if (i0 as usize) < lesect.len() {
                lesect[i0 as usize] += 1;
            }
            continue;
        }
        lesect[i as usize] += 1;
        i0 = i;
        if numsec != lesect[i as usize] {
            // `XSTEP_19` (`igesread.c:87-91`). Reading continues.
        }
        if i == 1 {
            line[72] = 0;
            let n = c_strlen(&line).min(72);
            let text = String::from_utf8_lossy(&line[..n]);
            tables.new_param(ARG_VIDE, &text);
        }
        if i == 2 {
            tables.set_global();
            loop {
                if lesect[i as usize] == 1 {
                    let mut n0 = 0usize;
                    if line[0] != b',' {
                        c_separ = line[2];
                        n0 = 3;
                    }
                    if line[n0 + 1] != c_separ {
                        c_fin = line[n0 + 3];
                    }
                }
                iges_param(&mut tables, &mut scan, &mut pstat, &line, c_separ, c_fin, 72);
                if pstat != 2 {
                    break;
                }
            }
        }
        if i == 3 {
            iges_dsect(&mut tables, &mut dstat, numsec, &line);
        }
        if i == 4 {
            iges_psect(&mut tables, &mut line);
            loop {
                iges_param(&mut tables, &mut scan, &mut pstat, &line, c_separ, c_fin, 64);
                if pstat != 2 {
                    break;
                }
            }
        }
    }

    let global = tables
        .header
        .iter()
        .map(raw_to_param)
        .collect::<Result<Vec<_>, _>>()?;
    let mut directory = Vec::new();
    let mut params_by_de = HashMap::new();
    for part in tables.parts() {
        directory.push(DirLine {
            type_number: part.typ,
            form_number: part.form,
            blank_status: part.blk & 15,
            transform: part.tra,
            de_seq: part.numpart,
        });
        let tokens = part
            .params
            .iter()
            .map(raw_to_param)
            .collect::<Result<Vec<_>, _>>()?;
        params_by_de.insert(part.numpart, tokens);
    }
    Ok(IgesRead {
        global,
        directory,
        params_by_de,
    })
}

fn raw_to_param(raw: &RawParam) -> Result<IgesParam, String> {
    let text = raw.parval.trim();
    match raw.typarg {
        ARG_VIDE => Ok(IgesParam::Void),
        ARG_INT | ARG_SIGN => text
            .parse::<i32>()
            .map(IgesParam::Int)
            .map_err(|_| format!("bad IGES integer '{text}'")),
        ARG_REAL | ARG_REXP | ARG_EXP | ARG_MEXP => {
            let norm = text.replace(['D', 'd'], "E");
            norm.parse::<f64>()
                .map(IgesParam::Real)
                .map_err(|_| format!("bad IGES real '{text}'"))
        }
        ARG_CHAR => Ok(IgesParam::Text(hollerith_body(&raw.parval))),
        ARG_QUID => Ok(IgesParam::Text(text.to_string())),
        other => Err(format!(
            "UNPORTED IGES parameter type {other}: igesread.h:67-75"
        )),
    }
}

fn hollerith_body(parval: &str) -> String {
    let b = parval.as_bytes();
    let mut i = 0usize;
    let mut n = 0usize;
    while i < b.len() && b[i].is_ascii_digit() {
        n = n * 10 + (b[i] - b'0') as usize;
        i += 1;
    }
    if i < b.len() && b[i] == b'H' {
        i += 1;
        let end = (i + n).min(b.len());
        return String::from_utf8_lossy(&b[i..end]).into_owned();
    }
    parval.to_string()
}
