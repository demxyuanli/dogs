//! Directory and parameter lists from `structiges.c`.
//!
//! Page allocators (`Maxparts`, `Maxpar`, `Maxcar`) are `Vec` storage. The
//! current-part cursor follows `iges_curpart` (`structiges.c:152-179`).

#[derive(Clone, Debug)]
pub(in crate::iges_read) struct RawParam {
    pub typarg: i32,
    pub parval: String,
}

#[derive(Clone, Debug)]
pub(in crate::iges_read) struct DirPart {
    pub numpart: i32,
    pub typ: i32,
    pub tra: i32,
    pub blk: i32,
    pub form: i32,
    pub params: Vec<RawParam>,
}

enum Recording {
    Start,
    Header,
    Part,
}

pub(in crate::iges_read) struct IgesTables {
    pub starts: Vec<RawParam>,
    pub header: Vec<RawParam>,
    parts: Vec<DirPart>,
    cur_index: Option<usize>,
    recording: Recording,
}

impl IgesTables {
    pub(in crate::iges_read) fn new() -> Self {
        Self {
            starts: Vec::new(),
            header: Vec::new(),
            parts: Vec::new(),
            cur_index: None,
            recording: Recording::Start,
        }
    }

    /// `iges_setglobal` (`structiges.c:129-130`).
    pub(in crate::iges_read) fn set_global(&mut self) {
        self.recording = Recording::Header;
    }

    /// `iges_newpart` (`structiges.c:134-148`).
    pub(in crate::iges_read) fn new_part(&mut self, numsec: i32) {
        self.parts.push(DirPart {
            numpart: numsec,
            typ: 0,
            tra: 0,
            blk: 0,
            form: 0,
            params: Vec::new(),
        });
        self.cur_index = Some(self.parts.len() - 1);
        self.recording = Recording::Part;
    }

    pub(in crate::iges_read) fn current_part_mut(&mut self) -> Option<&mut DirPart> {
        let ix = self.cur_index?;
        self.parts.get_mut(ix)
    }

    /// `iges_curpart` (`structiges.c:152-179`).
    pub(in crate::iges_read) fn select_part(&mut self, dnum: i32) {
        let Some(cur) = self.cur_index else {
            return;
        };
        if self.parts[cur].numpart == dnum {
            self.recording = Recording::Part;
            return;
        }
        let next = if cur + 1 < self.parts.len() { cur + 1 } else { 0 };
        if self.parts.get(next).is_some_and(|p| p.numpart == dnum) {
            self.cur_index = Some(next);
            self.recording = Recording::Part;
            return;
        }
        if let Some(ix) = self.parts.iter().position(|p| p.numpart == dnum) {
            self.cur_index = Some(ix);
            self.recording = Recording::Part;
        } else {
            self.cur_index = None;
        }
    }

    /// `iges_newparam` (`structiges.c:185-212`).
    pub(in crate::iges_read) fn new_param(&mut self, typarg: i32, parval: &str) {
        let list = match self.recording {
            Recording::Start => &mut self.starts,
            Recording::Header => &mut self.header,
            Recording::Part => match self.cur_index {
                Some(ix) => &mut self.parts[ix].params,
                None => return,
            },
        };
        list.push(RawParam {
            typarg,
            parval: parval.to_string(),
        });
    }

    /// `iges_addparam` (`structiges.c:215-228`).
    pub(in crate::iges_read) fn add_param(&mut self, parval: &str) {
        if parval.is_empty() {
            return;
        }
        let list = match self.recording {
            Recording::Start => &mut self.starts,
            Recording::Header => &mut self.header,
            Recording::Part => match self.cur_index {
                Some(ix) => &mut self.parts[ix].params,
                None => return,
            },
        };
        if let Some(last) = list.last_mut() {
            last.parval.push_str(parval);
        }
    }

    pub(in crate::iges_read) fn parts(&self) -> &[DirPart] {
        &self.parts
    }
}
