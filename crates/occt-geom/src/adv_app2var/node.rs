//! `AdvApp2Var_Node` (`AdvApp2Var_Node.hxx:27-88`, `AdvApp2Var_Node.cxx:47-54`).

/// Constraints at one end of an iso. Points are addressed as `SetPoint(iu, iv)`.
#[derive(Clone)]
pub(super) struct App2VarNode {
    ord_u: i32,
    ord_v: i32,
    coord: [f64; 2],
    points: Vec<[f64; 3]>,
    /// `myErrors`, initialized to 0 (`AdvApp2Var_Node.cxx:30`).
    errors: Vec<f64>,
}

impl App2VarNode {
    pub(super) fn new(iu: i32, iv: i32) -> Self {
        Self::at(iu, iv, 0.0, 0.0)
    }

    pub(super) fn at(iu: i32, iv: i32, u: f64, v: f64) -> Self {
        let nu = iu.max(0) + 1;
        let nv = iv.max(0) + 1;
        Self {
            ord_u: iu,
            ord_v: iv,
            coord: [u, v],
            points: vec![[0.0; 3]; (nu * nv) as usize],
            errors: vec![0.0; (nu * nv) as usize],
        }
    }

    /// `SetPoint` (`AdvApp2Var_Node.hxx:54`).
    pub(super) fn set_point(&mut self, iu: i32, iv: i32, pt: [f64; 3]) -> Result<(), String> {
        let at = self.index(iu, iv)?;
        self.points[at] = pt;
        Ok(())
    }

    /// `Point` (`AdvApp2Var_Node.hxx:57`).
    pub(super) fn point(&self, iu: i32, iv: i32) -> Result<[f64; 3], String> {
        Ok(self.points[self.index(iu, iv)?])
    }

    /// `Error` (`AdvApp2Var_Node.hxx:66`). The array stays 0; `SetError` is not called.
    pub(super) fn error(&self, iu: i32, iv: i32) -> Result<f64, String> {
        Ok(self.errors[self.index(iu, iv)?])
    }

    fn index(&self, iu: i32, iv: i32) -> Result<usize, String> {
        let nu = self.ord_u.max(0) + 1;
        if iu < 0 || iv < 0 || iu > self.ord_u.max(0) || iv > self.ord_v.max(0) {
            return Err("AdvApp2Var_Node::Point: AdvApp2Var_Node.hxx:57".into());
        }
        Ok((iu + iv * nu) as usize)
    }

    pub(super) fn coord(&self) -> [f64; 2] {
        self.coord
    }

    pub(super) fn order_u(&self) -> i32 {
        self.ord_u
    }

    pub(super) fn order_v(&self) -> i32 {
        self.ord_v
    }
}

fn xyz(contr: &[f64], order: i32) -> Result<[f64; 3], String> {
    let at = (order * 3) as usize;
    if order < 0 || at + 2 >= contr.len() {
        return Err("AdvApp2Var_Node::SetPoint: AdvApp2Var_Iso.cxx:309".into());
    }
    Ok([contr[at], contr[at + 1], contr[at + 2]])
}

/// Extremity constraints of one derivative (`AdvApp2Var_Iso.cxx:307-332`).
pub(super) fn store_ends(
    begin: &mut App2VarNode,
    end: &mut App2VarNode,
    iapp: i32,
    extrem_order: i32,
    isofav: i32,
    contr1: &[f64],
    contr2: &[f64],
) -> Result<(), String> {
    if extrem_order < 0 {
        return Ok(());
    }
    for ider in 0..=extrem_order {
        let pb = xyz(contr1, ider)?;
        let pe = xyz(contr2, ider)?;
        if isofav == 2 {
            begin.set_point(ider, iapp, pb)?;
            end.set_point(ider, iapp, pe)?;
        } else {
            begin.set_point(iapp, ider, pb)?;
            end.set_point(iapp, ider, pe)?;
        }
    }
    Ok(())
}
