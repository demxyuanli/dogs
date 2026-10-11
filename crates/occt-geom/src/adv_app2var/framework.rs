//! `AdvApp2Var_Framework` for `InitGrid(1)` (`AdvApp2Var_ApproxAFunc2Var.cxx:239-331`).
//!
//! `NbInt` is 1, so the regular-cutting loop (`:323-329`) does not run.

use crate::adv_approx::Cutting;

use super::node::App2VarNode;

/// One frontier iso (`AdvApp2Var_Iso` constructor, `AdvApp2Var_Iso.cxx:58-88`).
#[derive(Clone)]
pub(super) struct GridIso {
    pub isofav: i32,
    pub constant: f64,
    pub u0: f64,
    pub u1: f64,
    pub v0: f64,
    pub v1: f64,
    pub position: i32,
    pub approximated: bool,
    pub errors: Vec<(f64, f64)>,
    pub somtab: Vec<f64>,
    pub diftab: Vec<f64>,
    /// Canonical iso polynomial (`Iso::Polynom`, `AdvApp2Var_Iso.cxx:366`).
    pub polynom: Vec<f64>,
    /// `Iso::NbCoeff` (`AdvApp2Var_Iso.cxx:543`). Stays 0 until `MakeApprox`.
    pub nb_coeff: i32,
    pub ord_u: i32,
    pub ord_v: i32,
}

pub(super) struct IsoStatus {
    pub approximated: bool,
    pub has_result: bool,
    pub errors: Vec<(f64, f64)>,
    pub somtab: Vec<f64>,
    pub diftab: Vec<f64>,
    pub polynom: Vec<f64>,
    pub nb_coeff: i32,
}

pub(super) struct Framework {
    nodes: Vec<App2VarNode>,
    /// IsoV frontiers. `FirstNotApprox` scans these first.
    u_strips: Vec<Vec<GridIso>>,
    /// IsoU frontiers.
    v_strips: Vec<Vec<GridIso>>,
}

impl Framework {
    /// Four corners and four frontier isos (`AdvApp2Var_ApproxAFunc2Var.cxx:256-318`).
    pub(super) fn init_grid(u0: f64, u1: f64, v0: f64, v1: f64, iu: i32, iv: i32) -> Self {
        let iso = |isofav, constant, position| GridIso {
            isofav,
            constant,
            u0,
            u1,
            v0,
            v1,
            position,
            approximated: false,
            errors: Vec::new(),
            somtab: Vec::new(),
            diftab: Vec::new(),
            polynom: Vec::new(),
            nb_coeff: 0,
            ord_u: iu,
            ord_v: iv,
        };
        Self {
            nodes: vec![
                App2VarNode::at(iu, iv, u0, v0),
                App2VarNode::at(iu, iv, u1, v0),
                App2VarNode::at(iu, iv, u0, v1),
                App2VarNode::at(iu, iv, u1, v1),
            ],
            u_strips: vec![vec![iso(2, v0, 1), iso(2, v1, 2)]],
            v_strips: vec![vec![iso(1, u0, 3), iso(1, u1, 4)]],
        }
    }

    /// `FirstNotApprox` (`AdvApp2Var_Framework.cxx:77-106`).
    fn first_not_approx(&self) -> Option<(bool, i32, i32)> {
        for (u_side, strips) in [(true, &self.u_strips), (false, &self.v_strips)] {
            for (strip_index, strip) in strips.iter().enumerate() {
                for (iso_index, iso) in strip.iter().enumerate() {
                    if !iso.approximated {
                        return Some((u_side, (iso_index + 1) as i32, (strip_index + 1) as i32));
                    }
                }
            }
        }
        None
    }

    fn nb_iso_in_v(&self) -> i32 {
        self.u_strips.len() as i32 + 1
    }

    /// `FirstNode` (`AdvApp2Var_Framework.cxx:111-121`).
    fn first_node(&self, isofav: i32, index_iso: i32, index_strip: i32) -> i32 {
        let count = self.nb_iso_in_v();
        if isofav == 1 {
            count * (index_strip - 1) + index_iso
        } else {
            count * (index_iso - 1) + index_strip
        }
    }

    /// `LastNode` (`AdvApp2Var_Framework.cxx:125-135`).
    fn last_node(&self, isofav: i32, index_iso: i32, index_strip: i32) -> i32 {
        let count = self.nb_iso_in_v();
        if isofav == 1 {
            count * index_strip + index_iso
        } else {
            count * (index_iso - 1) + index_strip + 1
        }
    }

    fn iso(&self, u_side: bool, index_iso: i32, index_strip: i32) -> &GridIso {
        let strips = if u_side { &self.u_strips } else { &self.v_strips };
        &strips[(index_strip - 1) as usize][(index_iso - 1) as usize]
    }

    /// `ChangeIso` (`AdvApp2Var_Framework.cxx:139-147`).
    fn change_iso(&mut self, u_side: bool, index_iso: i32, index_strip: i32, iso: GridIso) {
        let strips = if u_side {
            &mut self.u_strips
        } else {
            &mut self.v_strips
        };
        strips[(index_strip - 1) as usize][(index_iso - 1) as usize] = iso;
    }

    /// `ComputeConstraints` (`AdvApp2Var_ApproxAFunc2Var.cxx:644-718`).
    ///
    /// A failed iso asks `PrefAndRec::Value`. A cut updates the network and
    /// the framework. Otherwise `OverwriteApprox` keeps a curve that exists.
    pub(super) fn compute_constraints<F>(
        &mut self,
        net: &mut super::network::Network,
        max_patches: i32,
        u_choice: &crate::adv_approx::PrefAndRec,
        v_choice: &crate::adv_approx::PrefAndRec,
        mut approx: F,
    ) -> Result<(), String>
    where
        F: FnMut(&GridIso, &mut App2VarNode, &mut App2VarNode) -> Result<IsoStatus, String>,
    {
        loop {
            let Some((u_side, index_iso, index_strip)) = self.first_not_approx() else {
                return Ok(());
            };
            let iso = self.iso(u_side, index_iso, index_strip).clone();
            let n1 = self.first_node(iso.isofav, index_iso, index_strip);
            let n2 = self.last_node(iso.isofav, index_iso, index_strip);
            if n1 < 1 || n2 < 1 || n1 as usize > self.nodes.len() || n2 as usize > self.nodes.len()
            {
                return Err(
                    "AdvApp2Var_Framework::FirstNode: AdvApp2Var_Framework.cxx:111".into(),
                );
            }
            let mut begin = self.nodes[(n1 - 1) as usize].clone();
            let mut end = self.nodes[(n2 - 1) as usize].clone();
            let status = approx(&iso, &mut begin, &mut end)?;
            if status.approximated {
                self.store_iso(u_side, index_iso, index_strip, iso, status, begin, end);
            } else {
                // `AdvApp2Var_ApproxAFunc2Var.cxx:671-716`.
                let t0 = Self::t0(&iso);
                let t1 = Self::t1(&iso);
                let iso_v = iso.isofav == 2;
                let mut dec = 0.0;
                let more = if iso_v {
                    u_choice.value(t0, t1, &mut dec)
                } else {
                    v_choice.value(t0, t1, &mut dec)
                };
                let nb_patch = if iso_v {
                    (net.nb_patch_in_u() + 1).checked_mul(net.nb_patch_in_v())
                } else {
                    (net.nb_patch_in_v() + 1).checked_mul(net.nb_patch_in_u())
                }
                .ok_or_else(|| {
                    "AdvApp2Var_ApproxAFunc2Var::ComputeConstraints: AdvApp2Var_ApproxAFunc2Var.cxx:677"
                        .to_string()
                })?;
                if nb_patch <= max_patches && more {
                    if iso_v {
                        net.update_in_u(dec)?;
                        self.update_in_u(dec)?;
                    } else {
                        net.update_in_v(dec)?;
                        self.update_in_v(dec)?;
                    }
                } else if status.has_result {
                    // `OverwriteApprox` (`AdvApp2Var_Iso.cxx:428-433`).
                    self.store_iso(u_side, index_iso, index_strip, iso, status, begin, end);
                } else {
                    return Err(
                        "AdvApp2Var_ApproxAFunc2Var : Curve Approximation Error".into(),
                    );
                }
            }
        }
    }

    fn store_iso(
        &mut self,
        u_side: bool,
        index_iso: i32,
        index_strip: i32,
        mut iso: GridIso,
        status: IsoStatus,
        begin: App2VarNode,
        end: App2VarNode,
    ) {
        let n1 = self.first_node(iso.isofav, index_iso, index_strip);
        let n2 = self.last_node(iso.isofav, index_iso, index_strip);
        iso.approximated = true;
        iso.errors = status.errors;
        iso.somtab = status.somtab;
        iso.diftab = status.diftab;
        iso.polynom = status.polynom;
        iso.nb_coeff = status.nb_coeff;
        self.change_iso(u_side, index_iso, index_strip, iso);
        self.nodes[(n1 - 1) as usize] = begin;
        self.nodes[(n2 - 1) as usize] = end;
    }

    /// `Node(U, V)` (`AdvApp2Var_Framework.cxx:151-164`).
    pub(super) fn node_at(&self, u: f64, v: f64) -> Result<&App2VarNode, String> {
        for node in &self.nodes {
            let coord = node.coord();
            if coord[0] == u && coord[1] == v {
                return Ok(node);
            }
        }
        self.nodes
            .last()
            .ok_or_else(|| "AdvApp2Var_Framework::Node: AdvApp2Var_Framework.cxx:163".to_string())
    }

    fn t0(iso: &GridIso) -> f64 {
        if iso.isofav == 1 { iso.v0 } else { iso.u0 }
    }

    fn t1(iso: &GridIso) -> f64 {
        if iso.isofav == 1 { iso.v1 } else { iso.u1 }
    }

    /// `findStripByBounds` (`AdvApp2Var_Framework.cxx:31-44`).
    fn find_strip(strips: &[Vec<GridIso>], first: f64, last: f64) -> usize {
        let mut index = 0;
        while index < strips.len()
            && (Self::t0(&strips[index][0]) != first || Self::t1(&strips[index][0]) != last)
        {
            index += 1;
        }
        index
    }

    /// `findIsoByConstant` (`AdvApp2Var_Framework.cxx:46-56`).
    fn find_iso(strip: &[GridIso], nb_iso: usize, constant: f64) -> usize {
        let mut index = 0;
        while index < nb_iso && index < strip.len() && strip[index].constant != constant {
            index += 1;
        }
        index
    }

    /// `IsoU` (`AdvApp2Var_Framework.cxx:168-177`).
    pub(super) fn iso_u(&self, u: f64, v0: f64, v1: f64) -> Result<&GridIso, String> {
        let strip_index = Self::find_strip(&self.v_strips, v0, v1);
        let strip = self
            .v_strips
            .get(strip_index)
            .ok_or_else(|| "AdvApp2Var_Framework::IsoU: AdvApp2Var_Framework.cxx:172".to_string())?;
        let iso_index = Self::find_iso(strip, self.u_strips.len(), u);
        strip.get(iso_index).ok_or_else(|| {
            "AdvApp2Var_Framework::IsoU: AdvApp2Var_Framework.cxx:176".to_string()
        })
    }

    /// `IsoV` (`AdvApp2Var_Framework.cxx:181-189`).
    pub(super) fn iso_v(&self, u0: f64, u1: f64, v: f64) -> Result<&GridIso, String> {
        let strip_index = Self::find_strip(&self.u_strips, u0, u1);
        let strip = self
            .u_strips
            .get(strip_index)
            .ok_or_else(|| "AdvApp2Var_Framework::IsoV: AdvApp2Var_Framework.cxx:185".to_string())?;
        let iso_index = Self::find_iso(strip, self.v_strips.len(), v);
        strip.get(iso_index).ok_or_else(|| {
            "AdvApp2Var_Framework::IsoV: AdvApp2Var_Framework.cxx:188".to_string()
        })
    }

    /// `UpdateInU` (`AdvApp2Var_Framework.cxx:194-288`).
    pub(super) fn update_in_u(&mut self, cutting: f64) -> Result<(), String> {
        let mut strip_index = 0usize;
        while strip_index < self.u_strips.len() {
            let first = &self.u_strips[strip_index][0];
            if first.u0 <= cutting && first.u1 >= cutting {
                break;
            }
            strip_index += 1;
        }
        if strip_index >= self.u_strips.len() || self.u_strips[strip_index].is_empty() {
            return Err("AdvApp2Var_Framework::UpdateInU: AdvApp2Var_Framework.cxx:202".into());
        }
        let u_first = self.u_strips[strip_index][0].u0;
        let u_last = self.u_strips[strip_index][0].u1;
        for iso in &mut self.u_strips[strip_index] {
            iso.change_span(u_first, cutting);
            iso.reset_approx();
        }
        let new_strip: Vec<GridIso> = self.u_strips[strip_index]
            .iter()
            .map(|iso| iso.spawn_span(cutting, u_last))
            .collect();
        self.u_strips.insert(strip_index + 1, new_strip);
        let u_at = strip_index;
        for strip in &mut self.v_strips {
            if u_at >= strip.len() {
                return Err("AdvApp2Var_Framework::UpdateInU: AdvApp2Var_Framework.cxx:253".into());
            }
            let u0 = strip[u_at].u0;
            let v0 = strip[u_at].v0;
            let v1 = strip[u_at].v1;
            strip[u_at].change_domain(u0, cutting, v0, v1);
            let new_iso = strip[u_at].spawn_constant(cutting, u0, cutting, v0, v1);
            strip.insert(u_at + 1, new_iso);
            let next = u_at + 2;
            if next >= strip.len() {
                return Err("AdvApp2Var_Framework::UpdateInU: AdvApp2Var_Framework.cxx:268".into());
            }
            let u1 = strip[next].u1;
            let v0 = strip[next].v0;
            let v1 = strip[next].v1;
            strip[next].change_domain(cutting, u1, v0, v1);
        }
        self.insert_u_nodes(cutting)
    }

    /// `UpdateInV` (`AdvApp2Var_Framework.cxx:292-381`).
    pub(super) fn update_in_v(&mut self, cutting: f64) -> Result<(), String> {
        let mut strip_index = 0usize;
        while strip_index < self.v_strips.len() {
            let first = &self.v_strips[strip_index][0];
            if first.v0 <= cutting && first.v1 >= cutting {
                break;
            }
            strip_index += 1;
        }
        if strip_index >= self.v_strips.len() || self.v_strips[strip_index].is_empty() {
            return Err("AdvApp2Var_Framework::UpdateInV: AdvApp2Var_Framework.cxx:295".into());
        }
        let v_first = self.v_strips[strip_index][0].v0;
        let v_last = self.v_strips[strip_index][0].v1;
        for iso in &mut self.v_strips[strip_index] {
            iso.change_span(v_first, cutting);
            iso.reset_approx();
        }
        let new_strip: Vec<GridIso> = self.v_strips[strip_index]
            .iter()
            .map(|iso| iso.spawn_v_span(cutting, v_last))
            .collect();
        self.v_strips.insert(strip_index + 1, new_strip);
        let v_at = strip_index;
        for strip in &mut self.u_strips {
            if v_at >= strip.len() {
                return Err("AdvApp2Var_Framework::UpdateInV: AdvApp2Var_Framework.cxx:347".into());
            }
            let u0 = strip[v_at].u0;
            let u1 = strip[v_at].u1;
            let v0 = strip[v_at].v0;
            strip[v_at].change_domain(u0, u1, v0, cutting);
            let new_iso = strip[v_at].spawn_constant(cutting, u0, u1, v0, cutting);
            strip.insert(v_at + 1, new_iso);
            let next = v_at + 2;
            if next >= strip.len() {
                return Err("AdvApp2Var_Framework::UpdateInV: AdvApp2Var_Framework.cxx:362".into());
            }
            let u0 = strip[next].u0;
            let u1 = strip[next].u1;
            let v1 = strip[next].v1;
            strip[next].change_domain(u0, u1, cutting, v1);
        }
        self.insert_v_nodes(cutting)
    }

    fn insert_u_nodes(&mut self, cutting: f64) -> Result<(), String> {
        if self.nodes.is_empty() {
            return Err("AdvApp2Var_Framework::UpdateInU: AdvApp2Var_Framework.cxx:274".into());
        }
        let mut prev = self.nodes[0].clone();
        let mut index = 1i32;
        while index < self.nodes.len() as i32 {
            let next_at = index as usize;
            if next_at >= self.nodes.len() {
                return Err("AdvApp2Var_Framework::UpdateInU: AdvApp2Var_Framework.cxx:277".into());
            }
            let next = self.nodes[next_at].clone();
            if prev.coord()[0] < cutting
                && next.coord()[0] > cutting
                && prev.coord()[1] == next.coord()[1]
            {
                self.nodes.insert(
                    index as usize,
                    App2VarNode::at(prev.order_u(), prev.order_v(), cutting, prev.coord()[1]),
                );
            }
            prev = next;
            index += 1;
        }
        Ok(())
    }

    fn insert_v_nodes(&mut self, cutting: f64) -> Result<(), String> {
        let stride = self.u_strips.len() as i32 + 1;
        let mut start = 1i32;
        while start <= self.nodes.len() as i32
            && self.nodes[(start - 1) as usize].coord()[1] < cutting
        {
            start += stride;
        }
        for iso_index in 1..=stride {
            let at = (iso_index - 1) as usize;
            if at >= self.nodes.len() {
                return Err("AdvApp2Var_Framework::UpdateInV: AdvApp2Var_Framework.cxx:375".into());
            }
            let source = self.nodes[at].clone();
            let insert_at = (start + iso_index - 2) as usize;
            if insert_at > self.nodes.len() {
                return Err("AdvApp2Var_Framework::UpdateInV: AdvApp2Var_Framework.cxx:379".into());
            }
            self.nodes.insert(
                insert_at,
                App2VarNode::at(source.order_u(), source.order_v(), source.coord()[0], cutting),
            );
        }
        Ok(())
    }
}

impl GridIso {
    /// Two-argument `ChangeDomain` (`AdvApp2Var_Iso.cxx:380-392`).
    fn change_span(&mut self, first: f64, last: f64) {
        if self.isofav == 1 {
            self.v0 = first;
            self.v1 = last;
        } else {
            self.u0 = first;
            self.u1 = last;
        }
    }

    /// Four-argument `ChangeDomain` (`AdvApp2Var_Iso.cxx:396-402`).
    fn change_domain(&mut self, u0: f64, u1: f64, v0: f64, v1: f64) {
        self.u0 = u0;
        self.u1 = u1;
        self.v0 = v0;
        self.v1 = v1;
    }

    /// `ResetApprox` (`AdvApp2Var_Iso.cxx:420-423`).
    fn reset_approx(&mut self) {
        self.approximated = false;
    }

    fn spawn_span(&self, first: f64, last: f64) -> Self {
        let mut iso = self.clone();
        iso.u0 = first;
        iso.u1 = last;
        iso.position = 0;
        iso.approximated = false;
        iso.errors.clear();
        iso.somtab.clear();
        iso.diftab.clear();
        iso.polynom.clear();
        iso.nb_coeff = 0;
        iso
    }

    fn spawn_v_span(&self, first: f64, last: f64) -> Self {
        let mut iso = self.clone();
        iso.v0 = first;
        iso.v1 = last;
        iso.position = 0;
        iso.approximated = false;
        iso.errors.clear();
        iso.somtab.clear();
        iso.diftab.clear();
        iso.polynom.clear();
        iso.nb_coeff = 0;
        iso
    }

    fn spawn_constant(&self, constant: f64, u0: f64, u1: f64, v0: f64, v1: f64) -> Self {
        let mut iso = self.clone();
        iso.constant = constant;
        iso.u0 = u0;
        iso.u1 = u1;
        iso.v0 = v0;
        iso.v1 = v1;
        iso.position = 0;
        iso.approximated = false;
        iso.errors.clear();
        iso.somtab.clear();
        iso.diftab.clear();
        iso.polynom.clear();
        iso.nb_coeff = 0;
        iso
    }
}
