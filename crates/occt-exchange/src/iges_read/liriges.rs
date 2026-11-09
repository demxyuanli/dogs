//! `iges_lire` (`liriges.c:34-165`).
//!
//! One character is taken first so a file that uses only CR still advances
//! (`liriges.c:53-60`). `modefnes` is false for `IGESFile_Read`.

pub(in crate::iges_read) struct LineReader {
    data: Vec<u8>,
    pos: usize,
    eof: bool,
    relire: bool,
}

impl LineReader {
    pub(in crate::iges_read) fn new(data: Vec<u8>) -> Self {
        Self {
            data,
            pos: 0,
            eof: false,
            relire: false,
        }
    }

    /// Section number 1..5, 0 at EOF, -1 when the line is skipped.
    pub(in crate::iges_read) fn lire(
        &mut self,
        numsec: &mut i32,
        line: &mut [u8; 100],
        modefnes: bool,
    ) -> i32 {
        if !self.relire {
            if *numsec == 0 {
                line[72] = b' ';
                line[79] = b' ';
            }
            line[0] = 0;
            if !self.fill_line(line, modefnes) {
                return 0;
            }
            if *numsec == 0 && line[72] != b'S' && line[79] == b' ' {
                line[0] = 0;
                if !self.fill_line(line, modefnes) {
                    return 0;
                }
            }
            if (line[0] & 128) != 0 && modefnes {
                for i in 0..80 {
                    line[i] ^= 150u8.wrapping_add((i as u8) & 3);
                }
            }
        }
        if self.eof {
            return 0;
        }
        if let Some(cut) = line.iter().position(|c| *c == 0x1A) {
            line[cut] = 0;
            return 0;
        }
        self.relire = false;
        if line[0] == 0 || line[0] == b'\n' || line[0] == b'\r' {
            return self.lire(numsec, line, modefnes);
        }
        if let Some(n) = scan_int(line, 73) {
            *numsec = n;
            if let Some(sec) = section_of(line[72]) {
                line[72] = 0;
                return sec;
            }
            if c_strlen(line) == 80
                && (line[79] == b'\n' || line[79] == b'\r')
                && line[0].is_ascii_digit()
            {
                let mut index = 1usize;
                while index < 100 && line[index].is_ascii_digit() {
                    index += 1;
                }
                if index < 100 && (line[index] == b'D' || line[index] == b'd') {
                    for j in (1..80).rev() {
                        line[j] = line[j - 1];
                    }
                    line[0] = b'.';
                }
                if let Some(sec) = section_of(line[72]) {
                    line[72] = 0;
                    return sec;
                }
            }
        }
        let mut i = c_strlen(line);
        while i > 0 && (line[i] == 0 || line[i] == b'\n' || line[i] == b'\r' || line[i] == b' ') {
            i -= 1;
        }
        if i != c_strlen(line) && i + 1 < 100 {
            line[i + 1] = 0;
        }
        while i > 0 && line[i].is_ascii_digit() {
            i -= 1;
        }
        let Some(n) = scan_int(line, i + 1) else {
            return -1;
        };
        *numsec = n;
        while i > 0 && line[i] == b' ' {
            i -= 1;
        }
        let typesec = line[i];
        if let Some(sec) = section_of(typesec) {
            line[i] = 0;
            return sec;
        }
        -1
    }

    fn fill_line(&mut self, line: &mut [u8; 100], modefnes: bool) -> bool {
        if modefnes {
            return self.fgets_at(line, 0, 99);
        }
        loop {
            let mut one = [0u8; 4];
            if !self.fgets_at(&mut one, 0, 2) {
                return false;
            }
            if one[0] != b'\r' && one[0] != b'\n' {
                line[0] = one[0];
                break;
            }
        }
        self.fgets_at(line, 1, 80)
    }

    /// `fgets` into `dest[at..]`, reading at most `max_with_nul - 1` bytes.
    fn fgets_at(&mut self, dest: &mut [u8], at: usize, max_with_nul: usize) -> bool {
        if max_with_nul <= 1 || at >= dest.len() {
            return false;
        }
        if self.pos >= self.data.len() {
            self.eof = true;
            return false;
        }
        let cap = (max_with_nul - 1).min(dest.len() - at);
        let mut n = 0usize;
        while n < cap && self.pos < self.data.len() {
            let c = self.data[self.pos];
            self.pos += 1;
            dest[at + n] = c;
            n += 1;
            if c == b'\n' {
                break;
            }
        }
        if self.pos >= self.data.len() && (n == 0 || dest[at + n - 1] != b'\n') {
            self.eof = true;
        }
        if at + n < dest.len() {
            dest[at + n] = 0;
        }
        true
    }
}

fn section_of(typesec: u8) -> Option<i32> {
    match typesec {
        b'S' => Some(1),
        b'G' => Some(2),
        b'D' => Some(3),
        b'P' => Some(4),
        b'T' => Some(5),
        _ => None,
    }
}

pub(in crate::iges_read) fn c_strlen(line: &[u8]) -> usize {
    line.iter().position(|c| *c == 0).unwrap_or(line.len())
}

fn scan_int(line: &[u8], from: usize) -> Option<i32> {
    if from >= line.len() {
        return None;
    }
    let mut i = from;
    while i < line.len() && line[i] == b' ' {
        i += 1;
    }
    if i >= line.len() || line[i] == 0 {
        return None;
    }
    let mut sign = 1i32;
    if line[i] == b'+' || line[i] == b'-' {
        if line[i] == b'-' {
            sign = -1;
        }
        i += 1;
    }
    if i >= line.len() || !line[i].is_ascii_digit() {
        return None;
    }
    let mut val = 0i32;
    while i < line.len() && line[i].is_ascii_digit() {
        val = val.saturating_mul(10).saturating_add((line[i] - b'0') as i32);
        i += 1;
    }
    Some(val * sign)
}
