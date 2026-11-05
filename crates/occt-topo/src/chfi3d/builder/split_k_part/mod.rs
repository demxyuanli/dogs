//! `ChFi3d_Builder::SplitKPart` and the `SurfData` copy statics of
//! `ChFi3d_Builder_SpKP.cxx`.
//!
//! Sources: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_SpKP.cxx`
//! (OCCT 8.0.0):
//! - `CpInterf` (`:117-138`), `CpSD` (`:145-158`),
//! - `SplitKPart` (`:749-1516`),
//! - `Trunc` in [`trunc`] (`:372-506`),
//! - `Tri` / `FillSD` / `AdjustParam` in [`domains`] (`:161-198`, `:547-744`),
//! - `ParamOnSpine` / `YaUnVoisin` / `ResetProl` / `ComputeAbscissa` in
//!   [`spine_params`] (`:202-368`, `:510-545`).
//!
//! The OCCT member reads `myDS`, `myVEMap`, `myEFMap`, `tolesp`, `tol2d`,
//! `tolapp3d` and calls `SearchFace`, so it is a method here, not a free
//! function. `SetData` is the port's `Vec<Arc<SurfData>>`: OCCT shares one
//! `ChFiDS_SurfData` between the caller and the sequence and copies with
//! `CpSD` between the entries; the port clones the incoming `SurfData` once and
//! then keeps the same append-then-`cp_sd` sequence, which leaves the entries
//! with exactly the parameter values OCCT ends up with.

use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;
use std::sync::Arc;

use occt_core::precision::Precision;

use crate::abs::ShapeType;
use crate::brep_top_adaptor_topol_tool::TopolTool;
use crate::chfi_ds::{SpineHandle, State, SurfData};
use crate::fclass2d::FaceState;
use crate::geom2d_hatch::{Hatcher, Intersector};
use crate::shape::{Edge, Face, TopoShape};
use crate::topo_tools::edge_vertices;
use crate::topopebrep::tool::adaptor::SurfaceAdaptor;

use super::ChFi3dBuilder;
use copy_data::{
    add_elements, cp_sd, remove_count, spine_first_status, spine_last_status, three_edges,
};
use domains::{adjust_param, fill_sd, tri};
use spine_params::{compute_abscissa, elementary_spine, spine_edge};

mod copy_data;
mod domains;
mod spine_params;
mod tail;
mod trunc;

impl ChFi3dBuilder {
    /// `SplitKPart(Data, SetData, Spine, Iedge, S1, I1, S2, I2, intf, intl)`
    /// (`cxx:749-1516`).
    ///
    /// Hatches the two `PCurveOnFace` against the tangency lines of `I1` / `I2`
    /// and splits `Data` per hatch domain. `intf` / `intl` are the "extension
    /// at the beginning / end of the spine" flags: they come in from
    /// `PerformSetOfKPart` and are re-derived by `SearchFace`.
    #[allow(clippy::too_many_arguments)]
    pub fn split_k_part(
        &mut self,
        data: &Arc<SurfData>,
        set_data: &mut Vec<Arc<SurfData>>,
        spine: &SpineHandle,
        iedge: i32,
        s1: &SurfaceAdaptor,
        i1: &Rc<RefCell<TopolTool>>,
        s2: &SurfaceAdaptor,
        i2: &Rc<RefCell<TopolTool>>,
        intf: &mut bool,
        intl: &mut bool,
    ) -> Result<bool, String> {
        // The hatching of each face is started by tangency lines.
        let pitol = Precision::PINTERSECTION;
        let tol2d = self.tol2d;
        let tolesp = self.tolesp;

        let mut m1: HashMap<i32, Option<Edge>> = HashMap::new();
        let mut m2: HashMap<i32, Option<Edge>> = HashMap::new();
        let mut ih1 = 0;
        let mut ih2 = 0;
        let mut nb1 = 1usize;
        let mut nb2 = 1usize;

        // Cutting of tangency lines (hatching).
        let inter = Intersector::new(pitol, pitol);
        let mut h1 = Hatcher::new(inter.clone(), tol2d, self.tolapp3d, false, false);
        let mut h2 = Hatcher::new(inter, tol2d, self.tolapp3d, false, false);

        // `Geom2dAdaptor_Curve ll1; ll1.Load(C1);` — the port keeps the p-curve
        // itself, whose own parameters are what `ll1.IsPeriodic()` /
        // `ll1.Period()` answer.
        let c1 = data.interference_on_s1().p_curve_on_face().cloned();
        if let Some(c1) = c1.as_ref() {
            add_elements(&mut h1, &mut m1, i1);
            ih1 = h1.trim_curve(Arc::clone(c1));
            h1.compute_domains_ind(ih1);
            if !h1.is_done_ind(ih1) {
                return Ok(false);
            }
            nb1 = h1.nb_domains(ih1);
            if nb1 == 0 {
                return Ok(false);
            }
        }

        let c2 = data.interference_on_s2().p_curve_on_face().cloned();
        if let Some(c2) = c2.as_ref() {
            add_elements(&mut h2, &mut m2, i2);
            ih2 = h2.trim_curve(Arc::clone(c2));
            h2.compute_domains_ind(ih2);
            if !h2.is_done_ind(ih2) {
                return Ok(false);
            }
            nb2 = h2.nb_domains(ih2);
            if nb2 == 0 {
                return Ok(false);
            }
        }

        // Return start and end vertexes of the Spine.
        let support = spine_edge(spine, iedge)?;
        let (v1, v2) = edge_vertices(&support);
        let mut bout1 = v1.map(|v| v.0);
        let mut bout2 = v2.map(|v| v.0);
        if support.orientation().is_reversed() {
            std::mem::swap(&mut bout1, &mut bout2);
        }
        let null_vertex = TopoShape::new(ShapeType::Vertex);
        let bout1_shape = bout1.as_ref().unwrap_or(&null_vertex);
        let bout2_shape = bout2.as_ref().unwrap_or(&null_vertex);

        // Return faces.
        let f1: Option<Face> = s1.face().cloned();
        let f2: Option<Face> = s2.face().cloned();
        let mut fbid: Option<Face> = None;

        // Restriction of SurfDatas by cut lines.
        let mut cd = (**data).clone();
        let (f1_shape, f2_shape) = (
            f1.as_ref()
                .map(|f| f.shape().clone())
                .unwrap_or_else(|| TopoShape::new(ShapeType::Face)),
            f2.as_ref()
                .map(|f| f.shape().clone())
                .unwrap_or_else(|| TopoShape::new(ShapeType::Face)),
        );
        cd.set_index_of_s1(self.my_ds.change_ds().add_shape(&f1_shape));
        cd.set_index_of_s2(self.my_ds.change_ds().add_shape(&f2_shape));

        let mut ind1 = vec![0usize; nb1];
        let mut ind2 = vec![0usize; nb2];
        let wref = 0.;

        // Switcher of access to surfs of SurfData (eap occ293).
        let mut on_s = 1;
        // FaceInterference to control length in OnSame situation (eap occ354).
        let mut cntl_fi_on_s = 0;

        if c1.is_none() && c2.is_none() {
            return Ok(false);
        } else if c1.is_none()
            || (nb1 == 1 && !h1.domain(ih1, 1).has_first_point())
        {
            // It is checked if the point 2d of the degenerated edge is in the
            // face.
            if c1.is_none() {
                let p2d1 = cd.get_2d_point(false, 1);
                let situ = i1.borrow_mut().classify(p2d1, 1.0e-8, false);
                if situ == FaceState::Out {
                    return Ok(false);
                }
            }

            // Parsing of domains by increasing parameters.
            if !tri(&mut h2, ih2, &mut ind2, &mut nb2, wref, 0., pitol) {
                return Ok(false);
            }
            // Filling of SurfData.
            for i in 1..=nb2 {
                let dom2 = h2.domain(ih2, ind2[i - 1] as i32).clone();
                fill_sd(
                    self.my_ds.change_ds(),
                    &mut cd,
                    &m2,
                    &dom2,
                    dom2.first_point().parameter(),
                    true,
                    2,
                    pitol,
                    bout1_shape,
                )?;
                fill_sd(
                    self.my_ds.change_ds(),
                    &mut cd,
                    &m2,
                    &dom2,
                    dom2.second_point().parameter(),
                    false,
                    2,
                    pitol,
                    bout2_shape,
                )?;
                set_data.push(Arc::new(cd.clone()));
                cd = cp_sd(self.my_ds.change_ds(), &cd)?;
            }
            if *intf {
                let sd = Arc::make_mut(&mut set_data[0]);
                let cp2 = sd.vertex_first_on_s2();
                let on_same = spine_first_status(spine)? == State::OnSame;
                if cp2.is_on_arc() && on_same {
                    *intf = !self.search_face(Some(spine), cp2, f2.as_ref(), &mut fbid);
                } else {
                    *intf = false;
                }
            }
            if *intl {
                let sd = Arc::make_mut(set_data.last_mut().expect("SplitKPart"));
                let cp2 = sd.vertex_last_on_s2();
                let on_same = spine_last_status(spine)? == State::OnSame;
                if cp2.is_on_arc() && on_same {
                    *intl = !self.search_face(Some(spine), cp2, f2.as_ref(), &mut fbid);
                } else {
                    *intl = false;
                }
            }
        } else if c2.is_none()
            || (nb2 == 1 && !h2.domain(ih2, 1).has_first_point())
        {
            // It is checked if the point 2d of the degenerated edge is in the
            // face.
            if c2.is_none() {
                let p2d2 = cd.get_2d_point(false, 2);
                let situ = i2.borrow_mut().classify(p2d2, 1.0e-8, false);
                if situ == FaceState::Out {
                    return Ok(false);
                }
            }

            // Parsing of domains by increasing parameters.
            if !tri(&mut h1, ih1, &mut ind1, &mut nb1, wref, 0., pitol) {
                return Ok(false);
            }
            // Filling of SurfData.
            for i in 1..=nb1 {
                let dom1 = h1.domain(ih1, ind1[i - 1] as i32).clone();
                fill_sd(
                    self.my_ds.change_ds(),
                    &mut cd,
                    &m1,
                    &dom1,
                    dom1.first_point().parameter(),
                    true,
                    1,
                    pitol,
                    bout1_shape,
                )?;
                fill_sd(
                    self.my_ds.change_ds(),
                    &mut cd,
                    &m1,
                    &dom1,
                    dom1.second_point().parameter(),
                    false,
                    1,
                    pitol,
                    bout2_shape,
                )?;
                set_data.push(Arc::new(cd.clone()));
                cd = cp_sd(self.my_ds.change_ds(), &cd)?;
            }
            if *intf {
                let sd = Arc::make_mut(&mut set_data[0]);
                let cp1 = sd.vertex_first_on_s1();
                let on_same = spine_first_status(spine)? == State::OnSame;
                if cp1.is_on_arc() && on_same {
                    *intf = !self.search_face(Some(spine), cp1, f1.as_ref(), &mut fbid);
                } else {
                    *intf = false;
                }
            }
            if *intl {
                let sd = Arc::make_mut(set_data.last_mut().expect("SplitKPart"));
                let cp1 = sd.vertex_last_on_s1();
                let on_same = spine_last_status(spine)? == State::OnSame;
                if cp1.is_on_arc() && on_same {
                    *intl = !self.search_face(Some(spine), cp1, f1.as_ref(), &mut fbid);
                } else {
                    *intl = false;
                }
            }
        } else {
            // Parsing of domains by increasing parameters; if there is a 2d
            // circle on a plane, one goes on 2D line of opposite face.
            let mut period1 = 0.;
            let mut period2 = 0.;
            let c1_periodic = c1.as_ref().map(|c| c.is_periodic()).unwrap_or(false);
            let c2_periodic = c2.as_ref().map(|c| c.is_periodic()).unwrap_or(false);
            if c1_periodic {
                if !tri(&mut h2, ih2, &mut ind2, &mut nb2, wref, 0., pitol) {
                    return Ok(false);
                }
                period1 = c1.as_ref().map(|c| c.period()).unwrap_or(0.);
                if !tri(&mut h1, ih1, &mut ind1, &mut nb1, wref, period1, pitol) {
                    return Ok(false);
                }
            } else {
                if !tri(&mut h1, ih1, &mut ind1, &mut nb1, wref, 0., pitol) {
                    return Ok(false);
                }
                if c2_periodic {
                    period2 = c2.as_ref().map(|c| c.period()).unwrap_or(0.);
                }
                if !tri(&mut h2, ih2, &mut ind2, &mut nb2, wref, period2, pitol) {
                    return Ok(false);
                }
            }

            // Filling of SurfData.
            let mut ion1: Vec<i32> = Vec::new();
            let mut ion2: Vec<i32> = Vec::new();
            for i in 1..=nb1 {
                let dom1 = h1.domain(ih1, ind1[i - 1] as i32).clone();
                let (mut f1p, mut l1p, acheval1) =
                    adjust_param(&dom1, wref, period1, pitol);
                let nbcoup1 = if acheval1 { 2 } else { 1 };
                for _icoup1 in 1..=nbcoup1 {
                    for j in 1..=nb2 {
                        let dom2 = h2.domain(ih2, j as i32).clone();
                        let (mut f2p, mut l2p, acheval2) =
                            adjust_param(&dom2, wref, period2, pitol);
                        let nbcoup2 = if acheval2 { 2 } else { 1 };
                        for _icoup2 in 1..=nbcoup2 {
                            if f2p <= l1p && f1p <= l2p {
                                if f1p >= f2p - tol2d {
                                    fill_sd(
                                        self.my_ds.change_ds(),
                                        &mut cd,
                                        &m1,
                                        &dom1,
                                        f1p,
                                        true,
                                        1,
                                        pitol,
                                        bout1_shape,
                                    )?;
                                }
                                if f2p >= f1p - tol2d {
                                    fill_sd(
                                        self.my_ds.change_ds(),
                                        &mut cd,
                                        &m2,
                                        &dom2,
                                        f2p,
                                        true,
                                        2,
                                        pitol,
                                        bout1_shape,
                                    )?;
                                }
                                if l1p >= l2p - tol2d {
                                    fill_sd(
                                        self.my_ds.change_ds(),
                                        &mut cd,
                                        &m2,
                                        &dom2,
                                        l2p,
                                        false,
                                        2,
                                        pitol,
                                        bout2_shape,
                                    )?;
                                }
                                if l2p >= l1p - tol2d {
                                    fill_sd(
                                        self.my_ds.change_ds(),
                                        &mut cd,
                                        &m1,
                                        &dom1,
                                        l1p,
                                        false,
                                        1,
                                        pitol,
                                        bout2_shape,
                                    )?;
                                }
                                set_data.push(Arc::new(cd.clone()));
                                cd = cp_sd(self.my_ds.change_ds(), &cd)?;
                                ion1.push(i as i32);
                                ion2.push(j as i32);
                            }
                            f2p += period2;
                            l2p += period2;
                        }
                    }
                    f1p += period1;
                    l1p += period1;
                }
            }

            // Processing of extensions. Do not truncate, otherwise, problems
            // of intersection for PerformCorner: after the call of SplitKPart
            // in PerformSetOfKPart, spines have been extended to the
            // extremities by Extend, and the extensions of SurfData are
            // preserved.
            if *intf {
                // We are at the beginning of the spine.
                let mut ifirst = 0usize;
                let mut dist = f64::MAX;
                let ed = elementary_spine(spine, iedge)?;
                for i1 in 1..=set_data.len() {
                    let cd1 = &set_data[i1 - 1];
                    let cp1 = cd1.vertex_first_on_s1();
                    let cp2 = cd1.vertex_first_on_s2();
                    if cp1.is_on_arc()
                        && !self.search_face(Some(spine), cp1, f1.as_ref(), &mut fbid)
                    {
                        let ptg = cd1.interference_on_s1().first_parameter();
                        let dsp = compute_abscissa(&ed, ptg);
                        if dsp.abs() < dist {
                            ifirst = i1;
                            dist = dsp.abs();
                        }
                    } else if cp2.is_on_arc()
                        && !self.search_face(Some(spine), cp2, f2.as_ref(), &mut fbid)
                    {
                        let ptg = cd1.interference_on_s2().first_parameter();
                        let dsp = compute_abscissa(&ed, ptg);
                        if dsp.abs() < dist {
                            ifirst = i1;
                            dist = dsp.abs();
                        }
                    }
                }
                if ifirst > 1 {
                    remove_count(set_data, 1, ifirst - 1)?;
                    remove_count(&mut ion1, 1, ifirst - 1)?;
                    remove_count(&mut ion2, 1, ifirst - 1)?;
                }
                if set_data.is_empty() {
                    return Ok(false);
                }

                let (cp1_arc, cp2_arc) = {
                    let cd2 = &set_data[0];
                    (
                        cd2.vertex_first_on_s1().is_on_arc(),
                        cd2.vertex_first_on_s2().is_on_arc(),
                    )
                };
                if cp1_arc && cp2_arc {
                    let cp1 = set_data[0].vertex_first_on_s1();
                    let found1 =
                        self.search_face(Some(spine), cp1, f1.as_ref(), &mut fbid);
                    let cp2 = set_data[0].vertex_first_on_s2();
                    let found2 =
                        self.search_face(Some(spine), cp2, f2.as_ref(), &mut fbid);
                    *intf = !found1 && !found2;
                } else if cp1_arc {
                    let sov = set_data[0].vertex_first_on_s2().clone();
                    let cp1 = set_data[0].vertex_first_on_s1();
                    if !self.search_face(Some(spine), cp1, f1.as_ref(), &mut fbid) {
                        let dom_index = ind2[(ion2[0] - 1) as usize];
                        let dom = h2.domain(ih2, dom_index as i32).clone();
                        let cd2 = Arc::make_mut(&mut set_data[0]);
                        fill_sd(
                            self.my_ds.change_ds(),
                            cd2,
                            &m2,
                            &dom,
                            dom.first_point().parameter(),
                            true,
                            2,
                            pitol,
                            bout1_shape,
                        )?;
                        let cp2 = cd2.vertex_first_on_s2();
                        let cp2_ok = !cp2.is_on_arc()
                            || (cp2.is_on_arc()
                                && self.search_face(
                                    Some(spine),
                                    cp2,
                                    f2.as_ref(),
                                    &mut fbid,
                                ));
                        if cp2_ok {
                            *cd2.vertex_first_on_s2_mut() = sov;
                            if spine_first_status(spine)? != State::OnSame {
                                let par = cd2.interference(1).parameter(true);
                                cd2.interference_mut(2).set_parameter(par, true);
                                *intf = false;
                            }
                        }
                    } else {
                        *intf = false;
                    }
                } else if cp2_arc {
                    let sov = set_data[0].vertex_first_on_s1().clone();
                    let cp2 = set_data[0].vertex_first_on_s2();
                    if !self.search_face(Some(spine), cp2, f2.as_ref(), &mut fbid) {
                        let dom_index = ind1[(ion1[0] - 1) as usize];
                        let dom = h1.domain(ih1, dom_index as i32).clone();
                        let cd2 = Arc::make_mut(&mut set_data[0]);
                        fill_sd(
                            self.my_ds.change_ds(),
                            cd2,
                            &m1,
                            &dom,
                            dom.first_point().parameter(),
                            true,
                            1,
                            pitol,
                            bout1_shape,
                        )?;
                        let cp1 = cd2.vertex_first_on_s1();
                        let cp1_ok = !cp1.is_on_arc()
                            || (cp1.is_on_arc()
                                && self.search_face(
                                    Some(spine),
                                    cp1,
                                    f1.as_ref(),
                                    &mut fbid,
                                ));
                        if cp1_ok {
                            *cd2.vertex_first_on_s1_mut() = sov;
                            if spine_first_status(spine)? != State::OnSame {
                                let par = cd2.interference(2).parameter(true);
                                cd2.interference_mut(1).set_parameter(par, true);
                                *intf = false;
                            }
                        }
                    } else {
                        *intf = false;
                    }
                }

                // Select `onS` switcher so that to get on spine params from
                // Interference with a face where both edges at corner are
                // OnSame (eap occ293).
                if *intf && spine_first_status(spine)? == State::OnSame {
                    let three_e = three_edges(&bout1, &support, &f1, &f2)?;
                    if crate::chfi3d::builder_helpers::faces::edge_state(
                        [&three_e[0], &three_e[1], &three_e[2]],
                        &self.my_ef_map,
                    )? == State::OnSame
                    {
                        on_s = 1;
                    } else {
                        on_s = 2;
                    }
                    cntl_fi_on_s = 3 - on_s;
                }
            }
            if *intl {
                // We are at the end of the spine.
                let mut ilast = 0usize;
                let mut dist = f64::MAX;
                let (f_spine, l_spine) = {
                    let guard = spine.read().map_err(|_| {
                        "ChFi3d_Builder::SplitKPart: poisoned spine".to_string()
                    })?;
                    (
                        guard.first_parameter_of(iedge).ok_or_else(|| {
                            "ChFi3d_Builder::SplitKPart: FirstParameter".to_string()
                        })?,
                        guard.last_parameter_of(iedge).ok_or_else(|| {
                            "ChFi3d_Builder::SplitKPart: LastParameter".to_string()
                        })?,
                    )
                };
                let ed = elementary_spine(spine, iedge)?;
                for i2 in 1..=set_data.len() {
                    let cd3 = &set_data[i2 - 1];
                    let cp1 = cd3.vertex_last_on_s1();
                    let cp2 = cd3.vertex_last_on_s2();
                    if cp1.is_on_arc()
                        && !self.search_face(Some(spine), cp1, f1.as_ref(), &mut fbid)
                    {
                        let ptg = cd3.interference_on_s1().last_parameter();
                        let dsp = -compute_abscissa(&ed, ptg) - f_spine + l_spine;
                        if dsp.abs() < dist {
                            ilast = i2;
                            dist = dsp.abs();
                        }
                    } else if cp2.is_on_arc()
                        && !self.search_face(Some(spine), cp2, f2.as_ref(), &mut fbid)
                    {
                        let ptg = cd3.interference_on_s2().last_parameter();
                        let dsp = -compute_abscissa(&ed, ptg) - f_spine + l_spine;
                        if dsp.abs() < dist {
                            ilast = i2;
                            dist = dsp.abs();
                        }
                    }
                }
                let lll = set_data.len();
                if ilast < lll {
                    remove_count(set_data, ilast + 1, lll)?;
                    remove_count(&mut ion1, ilast + 1, lll)?;
                    remove_count(&mut ion2, ilast + 1, lll)?;
                }
                if set_data.is_empty() {
                    return Ok(false);
                }

                let last = set_data.len();
                let (cp1_arc, cp2_arc) = {
                    let cd4 = &set_data[last - 1];
                    (
                        cd4.vertex_last_on_s1().is_on_arc(),
                        cd4.vertex_last_on_s2().is_on_arc(),
                    )
                };
                if cp1_arc && cp2_arc {
                    let cp1 = set_data[last - 1].vertex_last_on_s1();
                    let found1 =
                        self.search_face(Some(spine), cp1, f1.as_ref(), &mut fbid);
                    let cp2 = set_data[last - 1].vertex_last_on_s2();
                    let found2 =
                        self.search_face(Some(spine), cp2, f2.as_ref(), &mut fbid);
                    *intl = !found1 && !found2;
                } else if cp1_arc {
                    let sov = set_data[last - 1].vertex_last_on_s2().clone();
                    let cp1 = set_data[last - 1].vertex_last_on_s1();
                    if !self.search_face(Some(spine), cp1, f1.as_ref(), &mut fbid) {
                        let dom_index = ind2[(ion2[ion2.len() - 1] - 1) as usize];
                        let dom = h2.domain(ih2, dom_index as i32).clone();
                        let cd4 = Arc::make_mut(&mut set_data[last - 1]);
                        fill_sd(
                            self.my_ds.change_ds(),
                            cd4,
                            &m2,
                            &dom,
                            dom.second_point().parameter(),
                            false,
                            2,
                            pitol,
                            bout2_shape,
                        )?;
                        let cp2 = cd4.vertex_last_on_s2();
                        let cp2_ok = !cp2.is_on_arc()
                            || (cp2.is_on_arc()
                                && self.search_face(
                                    Some(spine),
                                    cp2,
                                    f2.as_ref(),
                                    &mut fbid,
                                ));
                        if cp2_ok {
                            *cd4.vertex_last_on_s2_mut() = sov;
                            if spine_last_status(spine)? != State::OnSame {
                                let par = cd4.interference(1).parameter(false);
                                cd4.interference_mut(2).set_parameter(par, false);
                                *intl = false;
                            }
                        }
                    } else {
                        *intl = false;
                    }
                } else if cp2_arc {
                    let sov = set_data[last - 1].vertex_last_on_s1().clone();
                    let cp2 = set_data[last - 1].vertex_last_on_s2();
                    if !self.search_face(Some(spine), cp2, f2.as_ref(), &mut fbid) {
                        let dom_index = ind1[(ion1[ion1.len() - 1] - 1) as usize];
                        let dom = h1.domain(ih1, dom_index as i32).clone();
                        let cd4 = Arc::make_mut(&mut set_data[last - 1]);
                        fill_sd(
                            self.my_ds.change_ds(),
                            cd4,
                            &m1,
                            &dom,
                            dom.second_point().parameter(),
                            false,
                            1,
                            pitol,
                            bout2_shape,
                        )?;
                        let cp1 = cd4.vertex_last_on_s1();
                        let cp1_ok = !cp1.is_on_arc()
                            || (cp1.is_on_arc()
                                && self.search_face(
                                    Some(spine),
                                    cp1,
                                    f1.as_ref(),
                                    &mut fbid,
                                ));
                        if cp1_ok {
                            *cd4.vertex_last_on_s1_mut() = sov;
                            if spine_last_status(spine)? != State::OnSame {
                                let par = cd4.interference(2).parameter(false);
                                cd4.interference_mut(1).set_parameter(par, false);
                                *intl = false;
                            }
                        }
                    } else {
                        *intl = false;
                    }
                }

                // Select `onS` switcher (eap occ293), end of the spine.
                if *intl && spine_last_status(spine)? == State::OnSame {
                    let three_e = three_edges(&bout2, &support, &f1, &f2)?;
                    if crate::chfi3d::builder_helpers::faces::edge_state(
                        [&three_e[0], &three_e[1], &three_e[2]],
                        &self.my_ef_map,
                    )? == State::OnSame
                    {
                        on_s = 1;
                    } else {
                        on_s = 2;
                    }
                    cntl_fi_on_s = 3 - on_s;
                }
            }
        }

        if !*intf {
            // SurfData are entirely suspended before the beginning of the edge.
            let mut i = 1usize;
            let mut okdoc = set_data.is_empty();
            while !okdoc {
                let ltg = set_data[i - 1].interference(on_s).last_parameter();
                let ed = elementary_spine(spine, iedge)?;
                let nl = compute_abscissa(&ed, ltg);
                if nl < -tolesp {
                    set_data.remove(i - 1);
                } else {
                    i += 1;
                }
                okdoc = set_data.is_empty() || i > set_data.len();
            }
        }
        if !*intl {
            // SurfData are entirely suspended after the end of the edge.
            let mut i = 1usize;
            let mut okdoc = set_data.is_empty();
            while !okdoc {
                let ftg = set_data[i - 1].interference(on_s).first_parameter();
                let (f_spine, l_spine) = {
                    let guard = spine.read().map_err(|_| {
                        "ChFi3d_Builder::SplitKPart: poisoned spine".to_string()
                    })?;
                    (
                        guard.first_parameter_of(iedge).ok_or_else(|| {
                            "ChFi3d_Builder::SplitKPart: FirstParameter".to_string()
                        })?,
                        guard.last_parameter_of(iedge).ok_or_else(|| {
                            "ChFi3d_Builder::SplitKPart: LastParameter".to_string()
                        })?,
                    )
                };
                let ed = elementary_spine(spine, iedge)?;
                let nl = compute_abscissa(&ed, ftg);
                if nl > (l_spine - f_spine + tolesp) {
                    set_data.remove(i - 1);
                } else {
                    i += 1;
                }
                okdoc = set_data.is_empty() || i > set_data.len();
            }
        }

        self.split_k_part_tail(
            set_data,
            spine,
            iedge,
            s1,
            s2,
            on_s,
            cntl_fi_on_s,
            *intf,
            *intl,
        )?;
        Ok(true)
    }
}
