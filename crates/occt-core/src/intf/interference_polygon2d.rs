//! Port of `Intf_InterferencePolygon2d`
//! (`src/ModelingAlgorithms/TKGeomAlgo/Intf/Intf_InterferencePolygon2d.hxx/.cxx`).

use crate::bnd::box2d::BndBox2d;
use crate::gp::GpPnt2d;
use crate::precision::{epsilon, Precision};

use super::interference::IntfInterference;
use super::polygon2d::IntfPolygon2d;
use super::section_point::{IntfPIType, IntfSectionPoint};
use super::tangent_zone::IntfTangentZone;

/// `PRCANG = Precision::Angular()` (`Intf_InterferencePolygon2d.cxx:30-32`).
const PRCANG: f64 = Precision::ANGULAR;

/// `Intf_InterferencePolygon2d` (`Intf_InterferencePolygon2d.hxx:33-76`).
#[derive(Clone, Debug)]
pub struct IntfInterferencePolygon2d {
    /// base `Intf_Interference`
    pub base: IntfInterference,
    /// `oClos` (`hxx:73`)
    o_clos: bool,
    /// `tClos` (`hxx:74`)
    t_clos: bool,
    /// `nbso` (`hxx:75`)
    nbso: usize,
}

impl Default for IntfInterferencePolygon2d {
    /// `Intf_InterferencePolygon2d()` (`cxx:35-42`).
    fn default() -> Self {
        Self::new()
    }
}

impl IntfInterferencePolygon2d {
    /// `Intf_InterferencePolygon2d()` (`cxx:35-42`).
    pub fn new() -> Self {
        Self {
            base: IntfInterference::new(false),
            o_clos: false,
            t_clos: false,
            nbso: 0,
        }
    }

    /// `Intf_InterferencePolygon2d(const Intf_Polygon2d& Obje1, const
    /// Intf_Polygon2d& Obje2)` (`cxx:48-68`).
    pub fn between(obje1: &dyn IntfPolygon2d, obje2: &dyn IntfPolygon2d) -> Self {
        let mut r = Self {
            base: IntfInterference::new(false),
            o_clos: false,
            t_clos: false,
            nbso: 0,
        };
        if !obje1.bounding().is_out_box(obje2.bounding()) {
            r.base.tolerance =
                obje1.deflection_over_estimation() + obje2.deflection_over_estimation();
            if r.base.tolerance == 0.0 {
                r.base.tolerance = epsilon(1000.0);
            }
            r.nbso = obje1.nb_segments();
            r.o_clos = obje1.closed();
            r.t_clos = obje2.closed();
            r.interference_two(obje1, obje2);
            r.clean();
        }
        r
    }

    /// `Intf_InterferencePolygon2d(const Intf_Polygon2d& Obje)` (`cxx:75-90`):
    /// self interference.
    pub fn self_intersection(obje: &dyn IntfPolygon2d) -> Self {
        let mut r = Self {
            base: IntfInterference::new(true),
            o_clos: false,
            t_clos: false,
            nbso: 0,
        };
        r.base.tolerance = obje.deflection_over_estimation() * 2.0;
        if r.base.tolerance == 0.0 {
            r.base.tolerance = epsilon(1000.0);
        }
        r.o_clos = obje.closed();
        r.t_clos = r.o_clos;
        r.interference_one(obje);
        r.clean();
        r
    }

    /// `Perform(Obje1, Obje2)` (`cxx:94-111`).
    pub fn perform_two(&mut self, obje1: &dyn IntfPolygon2d, obje2: &dyn IntfPolygon2d) {
        self.base.self_interference(false);
        if !obje1.bounding().is_out_box(obje2.bounding()) {
            self.base.tolerance =
                obje1.deflection_over_estimation() + obje2.deflection_over_estimation();
            if self.base.tolerance == 0.0 {
                self.base.tolerance = epsilon(1000.0);
            }
            self.nbso = obje1.nb_segments();
            self.o_clos = obje1.closed();
            self.t_clos = obje2.closed();
            self.interference_two(obje1, obje2);
            self.clean();
        }
    }

    /// `Perform(Obje)` (`cxx:115-126`).
    pub fn perform_one(&mut self, obje: &dyn IntfPolygon2d) {
        self.base.self_interference(true);
        self.base.tolerance = obje.deflection_over_estimation() * 2.0;
        if self.base.tolerance == 0.0 {
            self.base.tolerance = epsilon(1000.0);
        }
        self.o_clos = obje.closed();
        self.t_clos = self.o_clos;
        self.interference_one(obje);
        self.clean();
    }

    /// `Pnt2dValue(Index)` (`cxx:134-137`).
    pub fn pnt2d_value(&self, index: usize) -> GpPnt2d {
        let p = self.base.my_s_poins[index - 1].pnt();
        GpPnt2d::new(p.x(), p.y())
    }

    /// `Interference(Obje1, Obje2)` (`cxx:141-174`).
    fn interference_two(&mut self, obje1: &dyn IntfPolygon2d, obje2: &dyn IntfPolygon2d) {
        let n1 = self.nbso;
        let n2 = obje2.nb_segments();
        let d1 = obje1.deflection_over_estimation();
        let d2 = obje2.deflection_over_estimation();

        for i_obje1 in 1..=n1 {
            let mut b_so = BndBox2d::new();
            b_so.set_void();
            let (mut p1b, mut p1e) = (GpPnt2d::zero(), GpPnt2d::zero());
            obje1.segment(i_obje1, &mut p1b, &mut p1e);
            b_so.add_point(&p1b);
            b_so.add_point(&p1e);
            b_so.enlarge(d1);
            if !obje2.bounding().is_out_box(&b_so) {
                for i_obje2 in 1..=n2 {
                    let mut b_st = BndBox2d::new();
                    b_st.set_void();
                    let (mut p2b, mut p2e) = (GpPnt2d::zero(), GpPnt2d::zero());
                    obje2.segment(i_obje2, &mut p2b, &mut p2e);
                    b_st.add_point(&p2b);
                    b_st.add_point(&p2e);
                    b_st.enlarge(d2);
                    if !b_so.is_out_box(&b_st) {
                        self.intersect(i_obje1, i_obje2, &p1b, &p1e, &p2b, &p2e);
                    }
                }
            }
        }
    }

    /// `Interference(Obje)` (`cxx:178-210`).
    fn interference_one(&mut self, obje: &dyn IntfPolygon2d) {
        let n = obje.nb_segments();
        let d = obje.deflection_over_estimation();

        for i_obje1 in 1..=n {
            let mut b_so = BndBox2d::new();
            b_so.set_void();
            let (mut p1b, mut p1e) = (GpPnt2d::zero(), GpPnt2d::zero());
            obje.segment(i_obje1, &mut p1b, &mut p1e);
            b_so.add_point(&p1b);
            b_so.add_point(&p1e);
            b_so.enlarge(d);
            if !obje.bounding().is_out_box(&b_so) {
                for i_obje2 in (i_obje1 + 1)..=n {
                    let mut b_st = BndBox2d::new();
                    b_st.set_void();
                    let (mut p2b, mut p2e) = (GpPnt2d::zero(), GpPnt2d::zero());
                    obje.segment(i_obje2, &mut p2b, &mut p2e);
                    b_st.add_point(&p2b);
                    b_st.add_point(&p2e);
                    b_st.enlarge(d);
                    if !b_so.is_out_box(&b_st) {
                        self.intersect(i_obje1, i_obje2, &p1b, &p1e, &p2b, &p2e);
                    }
                }
            }
        }
    }

    /// `Clean()` (`cxx:214-305`).
    fn clean(&mut self) {
        let mut nb_it = self.base.my_t_zones.len();
        let mut decal = 0usize;
        let mut only1_seg = false;

        for ltz in 1..=nb_it {
            let mut tsp = 0usize;
            let mut tsps = 0usize;
            let (pr1mi, pr1ma) = self.base.my_t_zones[ltz - decal - 1].param_on_first();
            let delta1 = pr1ma - pr1mi;
            let (pr2mi, pr2ma) = self.base.my_t_zones[ltz - decal - 1].param_on_second();
            let delta2 = pr2ma - pr2mi;
            if delta1 < 1.0 && delta2 < 1.0 {
                only1_seg = true;
            }
            if delta1 == 0.0 || delta2 == 0.0 {
                only1_seg = true;
            }

            for lpi in 1..=self.base.my_t_zones[ltz - decal - 1].number_of_points() {
                let pi1 = *self.base.my_t_zones[ltz - decal - 1].get_point(lpi);
                if pi1.incidence() <= PRCANG {
                    tsp = 0;
                    tsps = 0;
                    break;
                }
                let (dim1, _addr1, _par1) = pi1.info_first();
                let (dim2, _addr2, _par2) = pi1.info_second();
                if dim1 == IntfPIType::Edge && dim2 == IntfPIType::Edge {
                    tsps = 0;
                    if tsp > 0 {
                        tsp = 0;
                        only1_seg = false;
                        break;
                    }
                    tsp = lpi;
                } else if dim1 != IntfPIType::External && dim2 != IntfPIType::External {
                    tsps = lpi;
                }
            }
            if tsp > 0 {
                let p = *self.base.my_t_zones[ltz - decal - 1].get_point(tsp);
                self.base.my_s_poins.push(p);
                self.base.my_t_zones.remove(ltz - decal - 1);
                decal += 1;
            } else if only1_seg && tsps != 0 {
                let p = *self.base.my_t_zones[ltz - decal - 1].get_point(tsps);
                self.base.my_s_poins.push(p);
                self.base.my_t_zones.remove(ltz - decal - 1);
                decal += 1;
            }
        }

        nb_it = self.base.my_s_poins.len();
        decal = 0;
        for lpi in 1..=nb_it {
            for ltz in 1..=self.base.my_t_zones.len() {
                let p = self.base.my_s_poins[lpi - decal - 1];
                if self.base.my_t_zones[ltz - 1].range_contains(&p) {
                    self.base.my_s_poins.remove(lpi - decal - 1);
                    decal += 1;
                    break;
                }
            }
        }
    }

    /// `Intersect(iObje1, iObje2, BegO, EndO, BegT, EndT)` (`cxx:308-819`).
    #[allow(clippy::too_many_arguments)]
    fn intersect(
        &mut self,
        i_obje1: usize,
        i_obje2: usize,
        beg_o: &GpPnt2d,
        end_o: &GpPnt2d,
        beg_t: &GpPnt2d,
        end_t: &GpPnt2d,
    ) {
        if self.base.self_intf && i_obje1.abs_diff(i_obje2) <= 1 {
            return;
        }

        let mut nbpi = 0usize;
        let mut par_o = [0.0f64; 9];
        let mut par_t = [0.0f64; 9];
        let mut the_pi: Vec<IntfSectionPoint> = Vec::new();
        let seg_t = end_t.xy().subtracted(beg_t.xy());
        let seg_o = end_o.xy().subtracted(beg_o.xy());

        let lg_t = seg_t.dot(&seg_t).sqrt();
        if lg_t <= 0.0 {
            return;
        }
        let lg_o = seg_o.dot(&seg_o).sqrt();
        if lg_o <= 0.0 {
            return;
        }

        let sig_ps = if seg_o.dot(&seg_t) > 0.0 { 1.0 } else { -1.0 };
        let floatgap = epsilon(lg_o + lg_t);
        let sin_teta = (seg_o.crossed(&seg_t) / lg_o) / lg_t;
        let ray_intf = if sin_teta > 0.0 {
            self.base.tolerance / sin_teta
        } else {
            0.0
        };

        // Interference <begO> <segT>
        let db_ot = beg_o.xy().subtracted(beg_t.xy()).crossed(&seg_t) / lg_t;
        let db_obt = beg_o.distance(beg_t);
        let db_oet = beg_o.distance(end_t);
        if db_ot.abs() <= self.base.tolerance {
            if db_obt <= self.base.tolerance {
                nbpi += 1;
                par_o[nbpi] = 0.0;
                par_t[nbpi] = 0.0;
                the_pi.push(IntfSectionPoint::with_2d(
                    beg_o,
                    IntfPIType::Vertex,
                    i_obje1 as i32,
                    0.0,
                    IntfPIType::Vertex,
                    i_obje2 as i32,
                    0.0,
                    sin_teta,
                ));
            }
            if db_oet <= self.base.tolerance {
                nbpi += 1;
                par_o[nbpi] = 0.0;
                par_t[nbpi] = 1.0;
                the_pi.push(IntfSectionPoint::with_2d(
                    beg_o,
                    IntfPIType::Vertex,
                    i_obje1 as i32,
                    0.0,
                    IntfPIType::Vertex,
                    (i_obje2 + 1) as i32,
                    0.0,
                    sin_teta,
                ));
            }
            if db_obt > self.base.tolerance
                && db_oet > self.base.tolerance
                && db_obt + db_oet <= lg_t + self.base.tolerance
            {
                nbpi += 1;
                par_o[nbpi] = 0.0;
                par_t[nbpi] = db_obt / lg_t;
                the_pi.push(IntfSectionPoint::with_2d(
                    beg_o,
                    IntfPIType::Vertex,
                    i_obje1 as i32,
                    0.0,
                    IntfPIType::Edge,
                    i_obje2 as i32,
                    par_t[nbpi],
                    sin_teta,
                ));
            }
        }

        // Interference <endO> <segT>
        let de_ot = end_o.xy().subtracted(beg_t.xy()).crossed(&seg_t) / lg_t;
        let de_obt = end_o.distance(beg_t);
        let de_oet = end_o.distance(end_t);
        if de_ot.abs() <= self.base.tolerance {
            if de_obt <= self.base.tolerance {
                nbpi += 1;
                par_o[nbpi] = 1.0;
                par_t[nbpi] = 0.0;
                the_pi.push(IntfSectionPoint::with_2d(
                    end_o,
                    IntfPIType::Vertex,
                    (i_obje1 + 1) as i32,
                    0.0,
                    IntfPIType::Vertex,
                    i_obje2 as i32,
                    0.0,
                    sin_teta,
                ));
            }
            if de_oet <= self.base.tolerance {
                nbpi += 1;
                par_o[nbpi] = 1.0;
                par_t[nbpi] = 1.0;
                the_pi.push(IntfSectionPoint::with_2d(
                    end_o,
                    IntfPIType::Vertex,
                    (i_obje1 + 1) as i32,
                    0.0,
                    IntfPIType::Vertex,
                    (i_obje2 + 1) as i32,
                    0.0,
                    sin_teta,
                ));
            }
            if de_obt > self.base.tolerance
                && de_oet > self.base.tolerance
                && de_obt + de_oet <= lg_t + self.base.tolerance
            {
                nbpi += 1;
                par_o[nbpi] = 1.0;
                par_t[nbpi] = de_obt / lg_t;
                the_pi.push(IntfSectionPoint::with_2d(
                    end_o,
                    IntfPIType::Vertex,
                    (i_obje1 + 1) as i32,
                    0.0,
                    IntfPIType::Edge,
                    i_obje2 as i32,
                    par_t[nbpi],
                    sin_teta,
                ));
            }
        }

        // Interference <begT> <segO>
        let db_to = beg_t.xy().subtracted(beg_o.xy()).crossed(&seg_o) / lg_o;
        if db_to.abs() <= self.base.tolerance
            && db_obt > self.base.tolerance
            && de_obt > self.base.tolerance
            && db_obt + de_obt <= lg_o + self.base.tolerance
        {
            nbpi += 1;
            par_o[nbpi] = db_obt / lg_o;
            par_t[nbpi] = 0.0;
            the_pi.push(IntfSectionPoint::with_2d(
                beg_t,
                IntfPIType::Edge,
                i_obje1 as i32,
                par_o[nbpi],
                IntfPIType::Vertex,
                i_obje2 as i32,
                0.0,
                sin_teta,
            ));
        }

        // Interference <endT> <segO>
        let de_to = end_t.xy().subtracted(beg_o.xy()).crossed(&seg_o) / lg_o;
        if de_to.abs() <= self.base.tolerance
            && db_oet > self.base.tolerance
            && de_oet > self.base.tolerance
            && db_oet + de_oet <= lg_o + self.base.tolerance
        {
            nbpi += 1;
            par_o[nbpi] = db_oet / lg_o;
            par_t[nbpi] = 1.0;
            the_pi.push(IntfSectionPoint::with_2d(
                end_t,
                IntfPIType::Edge,
                i_obje1 as i32,
                par_o[nbpi],
                IntfPIType::Vertex,
                (i_obje2 + 1) as i32,
                0.0,
                sin_teta,
            ));
        }

        let mut edge_sp = false;
        let mut par_osp = 0.0;
        let mut par_tsp = 0.0;

        if (db_ot - de_ot).abs() > floatgap && (db_to - de_to).abs() > floatgap {
            par_osp = db_ot / (db_ot - de_ot);
            par_tsp = db_to / (db_to - de_to);
            if db_ot * de_ot <= 0.0 && db_to * de_to <= 0.0 {
                edge_sp = true;
            } else if nbpi == 0 {
                return;
            }

            // If there is no interference it is necessary to take the points
            // segment by segment
            if nbpi == 0 && sin_teta > PRCANG {
                nbpi += 1;
                par_o[nbpi] = par_osp;
                par_t[nbpi] = par_tsp;
                the_pi.push(IntfSectionPoint::with_2d(
                    &GpPnt2d::new(beg_o.x() + seg_o.x() * par_osp, beg_o.y() + seg_o.y() * par_osp),
                    IntfPIType::Edge,
                    i_obje1 as i32,
                    par_osp,
                    IntfPIType::Edge,
                    i_obje2 as i32,
                    par_tsp,
                    sin_teta,
                ));
            } else if ray_intf >= self.base.tolerance {
                let delta_o = ray_intf / lg_o;
                let delta_t = ray_intf / lg_t;
                let mut par_odeb = par_osp - delta_o;
                let mut par_ofin = par_osp + delta_o;
                let mut par_tdeb = par_tsp - sig_ps * delta_t;
                let mut par_tfin = par_tsp + sig_ps * delta_t;
                if nbpi == 0 {
                    par_o[1] = par_odeb;
                    par_o[2] = par_ofin;
                    par_t[1] = par_tdeb;
                    par_t[2] = par_tfin;
                    while nbpi < 2 {
                        nbpi += 1;
                        let x = beg_o.x() + seg_o.x() * par_o[nbpi];
                        let y = beg_o.y() + seg_o.y() * par_o[nbpi];
                        the_pi.push(IntfSectionPoint::with_2d(
                            &GpPnt2d::new(x, y),
                            IntfPIType::External,
                            i_obje1 as i32,
                            par_o[nbpi],
                            IntfPIType::External,
                            i_obje2 as i32,
                            par_t[nbpi],
                            sin_teta,
                        ));
                    }
                } else if nbpi == 1 {
                    let mut ok = true;
                    if 0.0 < par_odeb && par_odeb < 1.0 && 0.0 < par_tdeb && par_tdeb < 1.0 {
                        par_o[nbpi + 1] = par_odeb;
                        par_t[nbpi + 1] = par_tdeb;
                    } else if 0.0 < par_ofin
                        && par_ofin < 1.0
                        && 0.0 < par_tfin
                        && par_tfin < 1.0
                    {
                        par_o[nbpi + 1] = par_ofin;
                        par_t[nbpi + 1] = par_tfin;
                    } else {
                        ok = false;
                    }

                    if ok {
                        let x = beg_o.x() + seg_o.x() * par_o[nbpi + 1];
                        let y = beg_o.y() + seg_o.y() * par_o[nbpi + 1];
                        if the_pi[0].pnt().distance(&crate::gp::GpPnt::new(x, y, 0.0))
                            >= self.base.tolerance / 4.0
                        {
                            nbpi += 1;
                            the_pi.push(IntfSectionPoint::with_2d(
                                &GpPnt2d::new(x, y),
                                IntfPIType::External,
                                i_obje1 as i32,
                                par_o[nbpi],
                                IntfPIType::External,
                                i_obje2 as i32,
                                par_t[nbpi],
                                sin_teta,
                            ));
                        }
                    }
                } else {
                    // more than one singularity
                    let mut par_omin = par_o[1];
                    let mut par_omax = par_o[1];
                    let mut par_tmin = par_t[1];
                    let mut par_tmax = par_t[1];
                    for i in 2..=nbpi {
                        par_omin = par_omin.min(par_o[i]);
                        par_omax = par_omax.max(par_o[i]);
                        par_tmin = par_tmin.min(par_t[i]);
                        par_tmax = par_tmax.max(par_t[i]);
                    }

                    let mut delta;
                    if par_odeb < 0.0 {
                        delta = -par_odeb;
                        par_odeb = 0.0;
                        par_tdeb += sig_ps * (delta * (delta_t / delta_o));
                    }
                    if par_ofin > 1.0 {
                        delta = par_ofin - 1.0;
                        par_ofin = 1.0;
                        par_tfin -= sig_ps * (delta * (delta_t / delta_o));
                    }
                    if sig_ps > 0.0 {
                        if par_tdeb < 0.0 {
                            delta = -par_tdeb;
                            par_tdeb = 0.0;
                            par_odeb += delta * (delta_o / delta_t);
                        }
                        if par_tfin > 1.0 {
                            delta = par_tfin - 1.0;
                            par_tfin = 1.0;
                            par_ofin -= delta * (delta_o / delta_t);
                        }
                    } else {
                        if par_tdeb > 1.0 {
                            delta = par_tdeb - 1.0;
                            par_tdeb = 1.0;
                            par_odeb += delta * (delta_o / delta_t);
                        }
                        if par_tfin < 0.0 {
                            delta = -par_tfin;
                            par_tfin = 0.0;
                            par_ofin -= delta * (delta_o / delta_t);
                        }
                    }

                    if (par_odeb < par_omin && par_omin > 0.0)
                        || (sig_ps > 0.0 && par_tdeb < par_tmin && par_tmin > 0.0)
                        || (sig_ps < 0.0 && par_tdeb > par_tmax && par_tmax < 1.0)
                    {
                        nbpi += 1;
                        par_o[nbpi] = par_odeb.max(0.0).min(1.0);
                        par_t[nbpi] = par_tdeb.max(0.0).min(1.0);
                        let x = beg_o.x() + seg_o.x() * par_o[nbpi];
                        let y = beg_o.y() + seg_o.y() * par_o[nbpi];
                        the_pi.push(IntfSectionPoint::with_2d(
                            &GpPnt2d::new(x, y),
                            IntfPIType::External,
                            i_obje1 as i32,
                            par_o[nbpi],
                            IntfPIType::External,
                            i_obje2 as i32,
                            par_t[nbpi],
                            sin_teta,
                        ));
                    }

                    if (par_ofin > par_omax && par_omax < 1.0)
                        || (sig_ps < 0.0 && par_tfin < par_tmin && par_tmin > 0.0)
                        || (sig_ps > 0.0 && par_tfin > par_tmax && par_tmax < 1.0)
                    {
                        nbpi += 1;
                        par_o[nbpi] = par_ofin.min(1.0).max(0.0);
                        par_t[nbpi] = par_tfin.min(1.0).max(0.0);
                        let x = beg_o.x() + seg_o.x() * par_o[nbpi];
                        let y = beg_o.y() + seg_o.y() * par_o[nbpi];
                        the_pi.push(IntfSectionPoint::with_2d(
                            &GpPnt2d::new(x, y),
                            IntfPIType::External,
                            i_obje1 as i32,
                            par_o[nbpi],
                            IntfPIType::External,
                            i_obje2 as i32,
                            par_t[nbpi],
                            sin_teta,
                        ));
                    }
                }
            }
        }

        // The points too close to each other are suspended
        let mut suppr;
        loop {
            suppr = false;
            let mut i = 2;
            while !suppr && i <= nbpi {
                let pim1 = *the_pi[i - 2].pnt();
                let pi = *the_pi[i - 1].pnt();
                let mut d = pi.distance(&pim1);
                d *= 50.0;
                if d < lg_t && d < lg_o {
                    for j in i..nbpi {
                        the_pi[j - 1] = the_pi[j];
                    }
                    nbpi -= 1;
                    suppr = true;
                }
                i += 1;
            }
            if !suppr {
                break;
            }
        }

        if nbpi == 1 {
            if edge_sp {
                the_pi[0] = IntfSectionPoint::with_2d(
                    &GpPnt2d::new(
                        beg_o.x() + seg_o.x() * par_osp,
                        beg_o.y() + seg_o.y() * par_osp,
                    ),
                    IntfPIType::Edge,
                    i_obje1 as i32,
                    par_osp,
                    IntfPIType::Edge,
                    i_obje2 as i32,
                    par_tsp,
                    sin_teta,
                );
                par_o[1] = par_osp;
                par_t[1] = par_tsp;
            }
            if !self.base.self_intf {
                let contains = self
                    .base
                    .my_s_poins
                    .iter()
                    .any(|p| the_pi[0].is_equal(p));
                if !contains {
                    self.base.my_s_poins.push(the_pi[0]);
                }
            } else if i_obje2 - i_obje1 != 1
                && (!self.o_clos || (i_obje1 != 1 && i_obje2 != self.nbso))
            {
                self.base.my_s_poins.push(the_pi[0]);
            }
        } else if nbpi >= 2 {
            let mut the_tz = IntfTangentZone::new();
            if nbpi == 2 {
                let p1 = the_pi[0];
                let p2 = the_pi[1];
                the_tz.polygon_insert(&p1);
                the_tz.polygon_insert(&p2);
            } else {
                let mut lmin = 1usize;
                let mut lmax = 1usize;
                for lpj in 2..=nbpi {
                    if par_o[lpj] < par_o[lmin] {
                        lmin = lpj;
                    } else if par_o[lpj] > par_o[lmax] {
                        lmax = lpj;
                    }
                }
                let p = the_pi[lmin - 1];
                the_tz.polygon_insert(&p);
                let p = the_pi[lmax - 1];
                the_tz.polygon_insert(&p);

                let mut ltmin = 1usize;
                let mut ltmax = 1usize;
                for lpj in 2..=nbpi {
                    if par_t[lpj] < par_t[ltmin] {
                        ltmin = lpj;
                    } else if par_t[lpj] > par_t[ltmax] {
                        ltmax = lpj;
                    }
                }
                if ltmin != lmin && ltmin != lmax {
                    let p = the_pi[ltmin - 1];
                    the_tz.polygon_insert(&p);
                }
                if ltmax != lmin && ltmax != lmax {
                    let p = the_pi[ltmax - 1];
                    the_tz.polygon_insert(&p);
                }
            }

            if edge_sp {
                let p = IntfSectionPoint::with_2d(
                    &GpPnt2d::new(
                        beg_o.x() + seg_o.x() * par_osp,
                        beg_o.y() + seg_o.y() * par_osp,
                    ),
                    IntfPIType::Edge,
                    i_obje1 as i32,
                    par_osp,
                    IntfPIType::Edge,
                    i_obje2 as i32,
                    par_tsp,
                    sin_teta,
                );
                the_tz.polygon_insert(&p);
            }

            let nbtz = self.base.my_t_zones.len();
            let mut l_index: Vec<usize> = Vec::new();
            for ltz in 1..=nbtz {
                if the_tz.has_common_range(&self.base.my_t_zones[ltz - 1]) {
                    l_index.push(ltz);
                }
            }
            if l_index.is_empty() {
                self.base.my_t_zones.push(the_tz);
            } else {
                let mut decal2 = 0usize;
                let indexfirst = l_index.remove(0);
                self.base.my_t_zones[indexfirst - 1].append_zone(&the_tz);
                while !l_index.is_empty() {
                    let index = l_index.remove(0);
                    let src = self.base.my_t_zones[index - decal2 - 1].clone();
                    self.base.my_t_zones[indexfirst - 1].append_zone(&src);
                    self.base.my_t_zones.remove(index - decal2 - 1);
                    decal2 += 1;
                }
            }
        }
    }
}
