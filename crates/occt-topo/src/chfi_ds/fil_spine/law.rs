//! `ChFiDS_FilSpine::ComputeLaw` — the radius law of one elementary spine.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_FilSpine.cxx:531-856`
//! (OCCT 8.0.0). The `mklaw` helper it calls is in [`super`].

use occt_core::elib::clib;
use occt_core::precision::CONFUSION;
use occt_geom::law::{LawComposite, LawConstant, LawFunction, LawInterpol};

use crate::chfi_ds::spine::Spine;
use super::*;

impl Spine {
    /// `ComputeLaw(Els)` (`cxx:531-856`): the composite of the radius
    /// evolution over the elementary spine `Els`.
    ///
    /// The spine is walked from the edge that starts `Els` (`Index(deb, true)`)
    /// to the one that ends it (`Index(fin, false)`), the run is wrapped when
    /// the guideline is periodic (`cxx:553-571`), the radii at both extremities
    /// are read from the neighbouring constants (`cxx:573-760`), and the inner
    /// edges are turned into laws by `mklaw` (`cxx:762-847`).
    pub fn compute_law(&self, els: &ElSpine) -> Option<LawComposite> {
        let tol3d = CONFUSION;
        let mut deb = els.first_parameter();
        let mut fin = els.last_parameter();
        let mut curdeb = deb;
        let mut curfin = fin;
        let ideb = self.index(deb, true)?;
        let ifin = self.index(fin, false)?;
        let len = self.nb_edges() as i32;
        // `cxx:539-541`: if the spine is periodic, attention to the index and
        // parameters.
        let mut spinedeb = self.first_parameter();
        let mut spinefin = self.last_parameter()?;

        let mut nbed = ifin - ideb + 1;
        let mut bidfin = ifin;

        let mut loi = LawComposite::new();
        let mut list: Vec<Arc<dyn LawFunction>> = Vec::new();
        let mut rdeb = 0.0f64;
        let mut rfin = 0.0f64;
        let mut icur = 1i32;
        let mut lawencours = false;

        if self.is_periodic() {
            // `cxx:553-564`.
            if deb < 0.0 && ideb > ifin {
                bidfin += len;
            } else if fin > self.last_parameter_of(len)? && ideb > ifin {
                bidfin += len;
            }
            nbed = bidfin - ideb + 1;
        }
        let mut ind = vec![0i32; usize::try_from(nbed).ok()?];
        // `cxx:566-571`.
        {
            let mut j = 1usize;
            for i in ideb..=bidfin {
                ind[j - 1] = ((i - 1) % len) + 1;
                j += 1;
            }
        }

        if els.is_periodic() {
            // `cxx:573-639`: a periodic composite is created at range, which is
            // eventually offset relatively to the elspine, to avoid a single
            // point at origin.
            loi.set_periodic();
            // `cxx:581-608`: is there a constant edge?
            let mut k = 1i32;
            while k <= len {
                if self.is_fil_constant_on_edge(k)? {
                    spinedeb = self.first_parameter_of(k)?;
                    deb = spinedeb;
                    curdeb = spinedeb;
                    spinefin = deb + self.period()?;
                    fin = spinefin;
                    for l in 1..=len {
                        ind[(l - 1) as usize] = ((k + l - 2) % len) + 1;
                    }
                    rdeb = self.fil_radius_of_edge_index(k)?;
                    rfin = rdeb;
                    icur += 1;
                    if len == 1 {
                        // `cxx:594-597`: because InPeriod will make 0.
                        curfin = self.last_parameter_of(k)?;
                    } else {
                        curfin = clib::in_period(self.last_parameter_of(k)?, spinedeb, spinefin);
                    }
                    let mut curloi = LawConstant::new();
                    curloi.set(rdeb, curdeb, curfin);
                    list.push(Arc::new(curloi));
                    curdeb = curfin;
                    break;
                }
                k += 1;
            }
            if k > len {
                // `cxx:609-639`: no constant edge.
                if self.fil.parandrad.is_empty() {
                    return None; // throw Standard_DomainError("Radius not defined")
                }
                let mut nbp = self.fil.parandrad.len() as i32;
                if nbp > 1 {
                    deb = self.fil.parandrad[0].x;
                    fin = deb + self.period()?;
                    if self.fil.parandrad[self.fil.parandrad.len() - 1].x - fin < -tol3d {
                        nbp += 1;
                    }
                } else {
                    nbp += 1;
                }
                let mut pr: Vec<(f64, f64)> = Vec::with_capacity(nbp as usize);
                for l in 1..nbp {
                    pr.push(self.fil.parandrad[(l - 1) as usize].coords());
                }
                // `cxx:634-635`: the table is closed at `fin` with the first
                // radius.
                let first_y = pr.first()?.1;
                pr.push((fin, first_y));
                // `cxx:636-638`: `Law_Interpol::Set(pr, true)`.
                let mut curloi = LawInterpol::new();
                curloi.set(&pr, true).ok()?;
                list.push(Arc::new(curloi));
                for law in list.drain(..) {
                    loi.add_law(law);
                }
                return Some(loi);
            }
        } else if self.is_periodic() {
            // `cxx:640-697`: the spine is periodic but the elspine is not.
            // start radius.
            if self.is_fil_constant_on_edge(ind[0])? {
                rdeb = self.fil_radius_of_edge_index(ind[0])?;
                curfin = self.last_parameter_of(ind[0])?;
                curfin = clib::in_period(curfin, spinedeb + tol3d, spinefin + tol3d);
                curfin = fin.min(curfin);
                let mut curloi = LawConstant::new();
                curloi.set(rdeb, curdeb, curfin);
                list.push(Arc::new(curloi));
                curdeb = curfin;
                icur += 1;
            } else {
                // `cxx:657-671`: there is inevitably kpart right before.
                let mut iprec = ind[0] - 1;
                if iprec == 0 {
                    iprec = len;
                }
                if self.is_fil_constant_on_edge(iprec)? {
                    rdeb = self.fil_radius_of_edge_index(iprec)?;
                } else {
                    return None; // throw "AppendLaw : previous constant is missing!"
                }
                lawencours = true;
            }
            // the radius at end.
            let ind_nbed = ind[(nbed - 1) as usize];
            if self.is_fil_constant_on_edge(ind_nbed)? {
                rfin = self.fil_radius_of_edge_index(ind_nbed)?;
            } else {
                // `cxx:680-695`: there is inevitably kpart right after.
                let mut isuiv = ind_nbed + 1;
                if isuiv == len + 1 {
                    isuiv = 1;
                }
                if self.is_fil_constant_on_edge(isuiv)? {
                    rfin = self.fil_radius_of_edge_index(isuiv)?;
                } else {
                    return None; // throw "AppendLaw : next constant is missing!"
                }
            }
        } else {
            // `cxx:698-760`: the spine is not periodic either.
            // the radius at start.
            if self.is_fil_constant_on_edge(ind[0])? {
                rdeb = self.fil_radius_of_edge_index(ind[0])?;
                curfin = fin.min(self.last_parameter_of(ind[0])?);
                let mut curloi = LawConstant::new();
                curloi.set(rdeb, curdeb, curfin);
                list.push(Arc::new(curloi));
                curdeb = curfin;
                icur += 1;
            } else if ind[0] > 1 {
                if self.is_fil_constant_on_edge(ind[0] - 1)? {
                    rdeb = self.fil_radius_of_edge_index(ind[0] - 1)?;
                } else {
                    return None; // throw "AppendLaw : previous constant is missing"
                }
                lawencours = true;
            } else if self.fil.parandrad.is_empty() {
                return None; // throw "AppendLaw : no radius on vertex"
            } else {
                rdeb = -1.0;
                lawencours = true;
            }
            // the radius at end.
            let ind_nbed = ind[(nbed - 1) as usize];
            if self.is_fil_constant_on_edge(ind_nbed)? {
                rfin = self.fil_radius_of_edge_index(ind_nbed)?;
            } else if ind_nbed < len {
                if self.is_fil_constant_on_edge(ind_nbed + 1)? {
                    rfin = self.fil_radius_of_edge_index(ind_nbed + 1)?;
                } else {
                    return None; // throw "AppendLaw : next constant is missing"
                }
            } else if self.fil.parandrad.is_empty() {
                return None; // throw "AppendLaw : no radius on vertex"
            } else {
                rfin = -1.0;
            }
        }

        // `cxx:762-847`: there are infos on the extremities of the elspine, all
        // edges are parsed.
        while icur <= nbed {
            let ind_icur = ind[(icur - 1) as usize];
            let ind_nbed = ind[(nbed - 1) as usize];
            if self.is_fil_constant_on_edge(ind_icur)? {
                let rcur = self.fil_radius_of_edge_index(ind_icur)?;
                if lawencours {
                    // `cxx:769-776`.
                    let mut temp: Vec<Arc<dyn LawFunction>> = Vec::new();
                    mklaw(
                        &mut temp,
                        &self.fil.parandrad,
                        curdeb,
                        curfin,
                        rdeb,
                        rcur,
                        self.is_periodic(),
                        spinedeb,
                        spinefin,
                        tol3d,
                    )?;
                    list.extend(temp);
                    lawencours = false;
                    curdeb = curfin;
                }
                // `cxx:777-800`.
                curfin = self.last_parameter_of(ind_icur)?;
                if self.is_periodic() {
                    curfin = clib::in_period(curfin, spinedeb + tol3d, spinefin + tol3d);
                    if ind_icur == ind_nbed {
                        // `cxx:781-793`: the curfin can be wrong if the last
                        // edge passes above the origin periodic spline.
                        let mut biddeb = self.first_parameter_of(ind_icur)?;
                        biddeb = clib::in_period(biddeb, spinedeb + tol3d, spinefin + tol3d);
                        if biddeb >= curfin {
                            curfin = fin;
                        } else {
                            curfin = fin.min(curfin);
                        }
                    } else {
                        curfin = fin.min(curfin);
                    }
                }
                if (curfin - curdeb) > tol3d {
                    // `cxx:801-808`.
                    rdeb = rcur;
                    let mut curloi = LawConstant::new();
                    curloi.set(rdeb, curdeb, curfin);
                    list.push(Arc::new(curloi));
                    curdeb = curfin;
                }
            } else {
                // `cxx:809-846`.
                curfin = self.last_parameter_of(ind_icur)?;
                if self.is_periodic() {
                    curfin = clib::in_period(curfin, spinedeb + tol3d, spinefin + tol3d);
                }
                curfin = fin.min(curfin);
                lawencours = true;
                if ind_icur == ind_nbed {
                    if self.is_periodic() {
                        // `cxx:822-836`: the curfin can be wrong if the last
                        // edge passes above the origin periodic spline.
                        let mut biddeb = self.first_parameter_of(ind_icur)?;
                        curfin = self.last_parameter_of(ind_icur)?;
                        biddeb = clib::in_period(biddeb, spinedeb + tol3d, spinefin + tol3d);
                        curfin = clib::in_period(curfin, spinedeb + tol3d, spinefin + tol3d);
                        if biddeb >= curfin {
                            curfin = fin;
                        } else {
                            curfin = fin.min(curfin);
                        }
                    } else if ind_icur == len {
                        // `cxx:838-842`: or if it is the end of spine with
                        // extension.
                        curfin = fin;
                    }
                    let mut temp: Vec<Arc<dyn LawFunction>> = Vec::new();
                    mklaw(
                        &mut temp,
                        &self.fil.parandrad,
                        curdeb,
                        curfin,
                        rdeb,
                        rfin,
                        self.is_periodic(),
                        spinedeb,
                        spinefin,
                        tol3d,
                    )?;
                    list.extend(temp);
                }
            }
            icur += 1;
        }

        // `cxx:550` / `:849-852`: OCCT declares `lastloi` and only tests it for
        // null, never assigning it, so the trailing append never runs.
        for law in list.drain(..) {
            loi.add_law(law);
        }
        Some(loi)
    }
}
