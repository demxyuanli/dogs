//! Port of `Intf_Interference`
//! (`src/ModelingAlgorithms/TKGeomAlgo/Intf/Intf_Interference.hxx/.cxx/.lxx`).

use super::section_line::IntfSectionLine;
use super::section_point::IntfSectionPoint;
use super::tangent_zone::IntfTangentZone;

/// `Intf_Interference` (`Intf_Interference.hxx:37-100`).
///
/// OCCT subclasses (`Intf_InterferencePolygon2d`) reach `mySPoins` / `mySLines`
/// / `myTZones` / `SelfIntf` / `Tolerance` directly; the port keeps them public
/// and composes rather than inherits.
#[derive(Clone, Debug)]
pub struct IntfInterference {
    /// `mySPoins` (`hxx:93`)
    pub my_s_poins: Vec<IntfSectionPoint>,
    /// `mySLines` (`hxx:94`)
    pub my_s_lines: Vec<IntfSectionLine>,
    /// `myTZones` (`hxx:95`)
    pub my_t_zones: Vec<IntfTangentZone>,
    /// `SelfIntf` (`hxx:96`)
    pub self_intf: bool,
    /// `Tolerance` (`hxx:97`)
    pub tolerance: f64,
}

impl Default for IntfInterference {
    fn default() -> Self {
        Self::new(false)
    }
}

impl IntfInterference {
    /// `Intf_Interference(Self)` (`Intf_Interference.cxx:26-30`).
    pub fn new(self_intf: bool) -> Self {
        Self {
            my_s_poins: Vec::new(),
            my_s_lines: Vec::new(),
            my_t_zones: Vec::new(),
            self_intf,
            tolerance: 0.0,
        }
    }

    /// `NbSectionPoints()` (`Interference.lxx:20-23`).
    pub fn nb_section_points(&self) -> usize {
        self.my_s_poins.len()
    }

    /// `PntValue(Index)` (`Interference.lxx:32-35`).
    pub fn pnt_value(&self, index: usize) -> &IntfSectionPoint {
        &self.my_s_poins[index]
    }

    /// `NbSectionLines()` (`Interference.lxx:41-44`).
    pub fn nb_section_lines(&self) -> usize {
        self.my_s_lines.len()
    }

    /// `LineValue(Index)` (`Interference.lxx:48-51`).
    pub fn line_value(&self, index: usize) -> &IntfSectionLine {
        &self.my_s_lines[index]
    }

    /// `NbTangentZones()` (`Interference.lxx:59-62`).
    pub fn nb_tangent_zones(&self) -> usize {
        self.my_t_zones.len()
    }

    /// `ZoneValue(Index)` (`Interference.lxx:66-69`).
    pub fn zone_value(&self, index: usize) -> &IntfTangentZone {
        &self.my_t_zones[index]
    }

    /// `GetTolerance()` (`Interference.lxx:73-76`).
    pub fn get_tolerance(&self) -> f64 {
        self.tolerance
    }

    /// `SelfInterference(Self)` (`Intf_Interference.cxx:37-43`).
    pub fn self_interference(&mut self, self_intf: bool) {
        self.self_intf = self_intf;
        self.my_s_poins.clear();
        self.my_s_lines.clear();
        self.my_t_zones.clear();
    }

    /// `Contains(LePnt)` (`Intf_Interference.cxx:287-304`).
    pub fn contains(&self, le_pnt: &IntfSectionPoint) -> bool {
        if self.my_s_lines.iter().any(|l| l.contains(le_pnt)) {
            return true;
        }
        self.my_t_zones.iter().any(|z| z.contains(le_pnt))
    }

    /// `Insert(const Intf_TangentZone&)` (`Intf_Interference.cxx:50-174`).
    pub fn insert_zone(&mut self, la_zone: &IntfTangentZone) -> bool {
        if self.my_t_zones.is_empty() {
            return false;
        }
        let mut lzin = 0usize;
        let mut lunp = 0usize;
        let mut lotp = 0usize;
        let mut lunl = 0usize;
        let mut lotl = 0usize;
        let mut same = false;
        let mut inserted = true;
        let nplz = la_zone.number_of_points();

        for iz in 1..=self.my_t_zones.len() {
            let npcz = self.my_t_zones[iz - 1].number_of_points();
            for ipz1 in 1..=npcz {
                let ipz0 = if ipz1 <= 1 { npcz } else { ipz1 - 1 };
                let ipz2 = (ipz1 % npcz) + 1;

                for ilz1 in 1..=nplz {
                    let ilz2 = (ilz1 % nplz) + 1;

                    let p_ipz1 = self.my_t_zones[iz - 1].get_point(ipz1);
                    if p_ipz1.is_equal(la_zone.get_point(ilz1)) {
                        let p_ipz0 = self.my_t_zones[iz - 1].get_point(ipz0);
                        if p_ipz0.is_equal(la_zone.get_point(ilz2)) {
                            lzin = iz;
                            lunp = ipz0;
                            lotp = ipz1;
                            lunl = ilz1;
                            lotl = ilz2;
                            same = false;
                            break;
                        } else {
                            let p_ipz2 = self.my_t_zones[iz - 1].get_point(ipz2);
                            if p_ipz2.is_equal(la_zone.get_point(ilz2)) {
                                lzin = iz;
                                lunp = ipz1;
                                lotp = ipz2;
                                lunl = ilz1;
                                lotl = ilz2;
                                same = true;
                                break;
                            } else {
                                lzin = iz;
                                lunp = ipz1;
                                lunl = ilz1;
                            }
                        }
                    }
                }
                if lotp != 0 {
                    break;
                }
            }
            if lotp != 0 {
                break;
            }
        }

        if lotp != 0 {
            let mut ilc = lotl + 1;
            while ((ilc - 1) % nplz) + 1 != lunl {
                let p = *la_zone.get_point(((ilc - 1) % nplz) + 1);
                self.my_t_zones[lzin - 1].insert_before(lotp, &p);
                if !same {
                    lotp += 1;
                }
                ilc += 1;
            }
        } else if lunp > 0 {
            let mut loop_ = false;
            let mut ilc = lunl;
            loop {
                let p = *la_zone.get_point(((ilc - 1) % nplz) + 1);
                self.my_t_zones[lzin - 1].insert_before(lunp, &p);
                lunp += 1;
                if loop_ && ((ilc - 1) % nplz) + 1 == lunl {
                    break;
                }
                loop_ = true;
                ilc += 1;
            }
        } else {
            inserted = false;
        }

        if inserted {
            let the_new = self.my_t_zones[lzin - 1].clone();
            self.my_t_zones.remove(lzin - 1);
            if !self.insert_zone(&the_new) {
                self.my_t_zones.push(the_new);
            }
        }
        inserted
    }

    /// `Insert(pdeb, pfin)` (`Intf_Interference.cxx:178-281`).
    pub fn insert_segment(&mut self, pdeb: &IntfSectionPoint, pfin: &IntfSectionPoint) {
        let mut inserted = false;
        let mut the_ls = 0usize;
        let mut begin = false;
        let mut the_bout = *pfin;

        for ils in 1..=self.my_s_lines.len() {
            let nd = self.my_s_lines[ils - 1].is_end(pdeb);
            let nf = self.my_s_lines[ils - 1].is_end(pfin);
            if nd == 1 {
                if nf > 1 {
                    self.my_s_lines[ils - 1].close();
                }
                inserted = true;
                the_ls = ils;
                begin = true;
                break;
            } else if nd > 1 {
                if nf == 1 {
                    self.my_s_lines[ils - 1].close();
                }
                inserted = true;
                the_ls = ils;
                begin = false;
                break;
            } else if nf == 1 {
                inserted = true;
                the_ls = ils;
                begin = true;
                the_bout = *pdeb;
                break;
            } else if nf > 1 {
                inserted = true;
                the_ls = ils;
                begin = false;
                the_bout = *pdeb;
                break;
            }
        }

        if !inserted {
            let mut la_ls = IntfSectionLine::new();
            la_ls.append_point(pdeb);
            la_ls.append_point(pfin);
            self.my_s_lines.push(la_ls);
        } else {
            let mut nd = 0usize;
            let mut merged = false;
            for ils in 1..=self.my_s_lines.len() {
                if ils != the_ls {
                    nd = self.my_s_lines[ils - 1].is_end(&the_bout);
                    if nd == 1 {
                        if begin {
                            self.my_s_lines[the_ls - 1].reverse();
                        }
                        let src = self.my_s_lines[the_ls - 1].clone();
                        self.my_s_lines[ils - 1].prepend_line(&src);
                        merged = true;
                        break;
                    } else if nd > 1 {
                        if !begin {
                            self.my_s_lines[the_ls - 1].reverse();
                        }
                        let src = self.my_s_lines[the_ls - 1].clone();
                        self.my_s_lines[ils - 1].append_line(&src);
                        merged = true;
                        break;
                    }
                }
            }
            if nd > 0 {
                if merged {
                    self.my_s_lines.remove(the_ls - 1);
                }
            } else if begin {
                self.my_s_lines[the_ls - 1].prepend_point(&the_bout);
            } else {
                self.my_s_lines[the_ls - 1].append_point(&the_bout);
            }
        }
    }

    /// `Dump()` (`Intf_Interference.cxx:310-327`). OCCT writes to `std::cout`;
    /// this port keeps the method as a no-op hook.
    pub fn dump(&self) {}
}
