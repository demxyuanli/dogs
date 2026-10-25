//! `ChFiDS_FilSpine` — the guideline of a fillet: `ChFiDS_Spine` plus the
//! radius evolution law.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFiDS/ChFiDS_FilSpine.hxx:34-102`,
//! `ChFiDS_FilSpine.cxx:41-929` (OCCT 8.0.0).
//!
//! The radius is described twice, exactly as OCCT does: by the
//! `(parameter, radius)` table `parandrad` (`NCollection_Sequence<gp_XY>`,
//! `X` = abscissa on the guideline, `Y` = radius), and by one
//! `Law_Composite` per elementary spine in `laws`, built by [`Spine::compute_law`]
//! from that table.
//!
//! `ChFiDS_FilSpine` derives from `ChFiDS_Spine`; a `SpineHandle` is
//! `Arc<RwLock<Spine>>` and cannot hold a derived class, so the two members
//! OCCT declares at `hxx:98-99` live in [`super::spine::Spine`] as
//! [`FilSpineData`] and the methods below are written on `Spine` (see
//! [`super::spine::SpineClass`]). The names carry the `fil` qualifier where
//! OCCT would have resolved the override through the vtable.
//!
//! Divergences that are not heuristics:
//! - OCCT raises `Standard_DomainError` from `Radius(IE)` (`cxx:342`),
//!   `Radius()` (`cxx:360-363`), `ChangeLaw` (`cxx:876-885`), `mklaw`
//!   (`cxx:443-446`) and `ComputeLaw` (`cxx:611-614`, `:669-671`, `:691-694`,
//!   `:720-723`, `:724-726`, `:747-749`, `:751-753`); the port answers `None`,
//! - `occ::handle<Law_Function>` becomes `Arc<dyn LawFunction>`, and
//!   `occ::handle<Law_Composite>` becomes `Arc<RwLock<LawComposite>>` because
//!   `ChangeLaw` hands the composite's elementary law out for writing
//!   (`ChFi3d_FilBuilder.cxx:364`, `:396`),
//! - `ChangeLaw` returns the `(composite, W)` pair instead of
//!   `occ::handle<Law_Function>&`: Rust cannot return a reference into a
//!   `RwLock`, so the caller takes the write guard and calls
//!   `LawComposite::elementary_law_mut(W)` (`Law_Composite.cxx:222-226`).
//!
//! Module layout:
//! ```text
//! mod  FilSpineData, the radius table writers, Radius / IsConstant,
//!      AppendElSpine / AppendLaw / Law / ChangeLaw, MaxRadFromSeqAndLaws
//! law  ComputeLaw
//! ```

use std::sync::{Arc, RwLock};

use occt_core::elib::clib;
use occt_core::gp::GpXY;
use occt_core::precision::{CONFUSION, REAL_SMALL};
use occt_geom::law::{LawComposite, LawConstant, LawFunction, LawInterpol, LawS};

use crate::shape::TopoShape;

use super::el_spine::ElSpine;
use super::spine::Spine;

mod law;

/// `ChFiDS_FilSpine::parandrad` and `::laws` (`ChFiDS_FilSpine.hxx:98-99`).
#[derive(Clone)]
pub struct FilSpineData {
    /// `parandrad`: `X` is the abscissa on the spine, `Y` the radius
    /// (`NCollection_Sequence<gp_XY>`; a `Vec` keeps OCCT's 1-based reads
    /// translated at the accessor).
    parandrad: Vec<GpXY>,
    /// `laws`: one composite per elementary spine, appended by
    /// [`Spine::append_fil_el_spine`] and read by [`Spine::fil_law`].
    laws: Vec<Arc<RwLock<LawComposite>>>,
}

impl FilSpineData {
    /// The member values `ChFiDS_FilSpine::ChFiDS_FilSpine` (`cxx:41-46`)
    /// leaves behind: both sequences empty.
    pub fn new() -> Self {
        Self {
            parandrad: Vec::new(),
            laws: Vec::new(),
        }
    }

    /// `parandrad`.
    pub fn parandrad(&self) -> &[GpXY] {
        &self.parandrad
    }

    /// `laws`.
    pub fn laws(&self) -> &[Arc<RwLock<LawComposite>>] {
        &self.laws
    }
}

impl Default for FilSpineData {
    fn default() -> Self {
        Self::new()
    }
}

/// `Law_Composite` carries no `Debug` (it caches through a `Mutex`); the port
/// prints the shape of the members, like `Spine`'s other handle members.
impl std::fmt::Debug for FilSpineData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FilSpineData")
            .field("parandrad", &self.parandrad)
            .field("nb_laws", &self.laws.len())
            .finish()
    }
}

/// `mklaw` (`ChFiDS_FilSpine.cxx:383-527`), file-static in OCCT.
///
/// Appends to `res` the laws that carry the radius between `curdeb` and
/// `curfin`: a `Law_Constant`, a `Law_S`, or the interpolation of the
/// `parandrad` points that fall inside the range. `Rdeb` / `Rfin` `< 0` mean
/// "unknown, read it from the first / last point" — that is how OCCT marks an
/// extremity whose constant has not been found.
fn mklaw(
    res: &mut Vec<Arc<dyn LawFunction>>,
    pr: &[GpXY],
    curdeb: f64,
    curfin: f64,
    rdeb: f64,
    rfin: f64,
    recadre: bool,
    deb: f64,
    fin: f64,
    tol3d: f64,
) -> Option<()> {
    let mut npr: Vec<GpXY> = Vec::new();
    let mut rad = rdeb;
    let mut raf = rfin;
    let mut yaunpointsurledeb = false;
    let mut yaunpointsurlefin = false;
    for cur in pr {
        // `cxx:404-407`.
        let mut wcur = cur.x;
        if recadre {
            wcur = clib::in_period(wcur, deb, fin);
        }
        if curdeb - tol3d <= wcur && wcur <= curfin + tol3d {
            if wcur - curdeb < tol3d {
                // `cxx:410-421`.
                yaunpointsurledeb = true;
                let mut ncur = *cur;
                if rdeb < 0.0 {
                    rad = cur.y;
                }
                ncur.set_coord(curdeb, rad);
                npr.push(ncur);
            } else if curfin - wcur < tol3d {
                // `cxx:422-431`.
                yaunpointsurlefin = true;
                let mut ncur = *cur;
                if rfin < 0.0 {
                    raf = cur.y;
                }
                ncur.set_coord(curfin, raf);
                npr.push(ncur);
            } else {
                npr.push(GpXY::new(wcur, cur.y));
            }
        }
    }

    if npr.is_empty() {
        // `cxx:440-459`.
        if rdeb < 0.0 && rfin < 0.0 {
            return None; // throw Standard_DomainError("Impossible to create the law")
        } else if rdeb < 0.0 || rfin < 0.0 {
            let r = if rfin < 0.0 { rdeb } else { rfin };
            let mut loi = LawConstant::new();
            loi.set(r, curdeb, curfin);
            res.push(Arc::new(loi));
        } else {
            let mut loi = LawS::new();
            loi.set(curdeb, rdeb, curfin, rfin).ok()?;
            res.push(Arc::new(loi));
        }
        return Some(());
    }

    // `cxx:460-470`.
    if !yaunpointsurledeb && rdeb >= 0.0 {
        npr.push(GpXY::new(curdeb, rdeb));
    }
    if !yaunpointsurlefin && rfin >= 0.0 {
        npr.push(GpXY::new(curfin, rfin));
    }
    let mut nbp = npr.len() as i32;
    // `cxx:472-485`: the same double loop (a bubble sort by abscissa).
    let mut i = 1i32;
    while i < nbp {
        let mut j = i + 1;
        while j <= nbp {
            if npr[(i - 1) as usize].x > npr[(j - 1) as usize].x {
                npr.swap((i - 1) as usize, (j - 1) as usize);
            }
            j += 1;
        }
        i += 1;
    }
    // `cxx:486-500`: duplicates are removed.
    let mut fini = nbp <= 1;
    let mut i = 1i32;
    while !fini {
        if (npr[(i - 1) as usize].x - npr[i as usize].x).abs() < tol3d {
            npr.remove((i - 1) as usize);
            nbp -= 1;
        } else {
            i += 1;
        }
        fini = i >= nbp;
    }

    // `cxx:502-507`.
    if rad < 0.0 {
        let first = *npr.first()?;
        let mut loi = LawConstant::new();
        loi.set(first.y, curdeb, first.x);
        res.push(Arc::new(loi));
    }
    // `cxx:508-518`: `Law_Interpol::Set(tpr, 0., 0., false)`, i.e. the
    // `Set(ParAndRad, Dd, Df, Periodic)` overload with both tangents at zero
    // (`Law_Interpol.cxx:108-142`).
    if nbp > 1 {
        let mut tpr: Vec<(f64, f64)> = Vec::with_capacity(nbp as usize);
        for l in 1..=nbp {
            tpr.push(npr[(l - 1) as usize].coords());
        }
        let mut curloi = LawInterpol::new();
        curloi.set_with_derivatives(&tpr, 0.0, 0.0, false).ok()?;
        res.push(Arc::new(curloi));
    }
    // `cxx:519-524`.
    if raf < 0.0 {
        let last = *npr.last()?;
        let mut loi = LawConstant::new();
        loi.set(last.y, last.x, curfin);
        res.push(Arc::new(loi));
    }
    Some(())
}

impl Spine {
    /// `ChFiDS_FilSpine::Reset(AllData = false)` (`cxx:51-81`), which overrides
    /// the virtual `ChFiDS_Spine::Reset`.
    pub fn reset_fil(&mut self, all_data: bool) -> Option<()> {
        self.reset(all_data)?;
        self.fil.laws.clear();
        if all_data {
            self.fil.parandrad.clear();
        } else {
            // `cxx:59-77`: complete `parandrad` so that it spans the whole
            // spine, and close it for a periodic spine.
            let spinedeb = self.first_parameter();
            let spinefin = self.last_parameter()?;
            // `cxx:64-65`: `parandrad.First()` / `Last()` raise on an empty
            // sequence.
            let mut first_uand_r = *self.fil.parandrad.first()?;
            let mut last_uand_r = *self.fil.parandrad.last()?;
            if (spinedeb - first_uand_r.x).abs() > REAL_SMALL {
                first_uand_r.set_x(spinedeb);
                self.fil.parandrad.insert(0, first_uand_r);
            }
            if (spinefin - last_uand_r.x).abs() > REAL_SMALL {
                last_uand_r.set_x(spinefin);
                self.fil.parandrad.push(last_uand_r);
            }
            if self.is_periodic() {
                // `cxx:78-80`.
                let n = self.fil.parandrad.len();
                let first_y = self.fil.parandrad[0].y;
                self.fil.parandrad[n - 1].set_y(first_y);
            }
        }
        Some(())
    }

    /// `SetRadius(Radius, E)` (`cxx:85-92`): the constant radius on edge `E`,
    /// stored at the two ends (`X = 0` and `X = 1` of the edge's own span).
    pub fn set_radius_on_edge(&mut self, radius: f64, e: &TopoShape) -> Option<()> {
        self.splitdone = false;
        let ie = self.index_of_edge(e);
        let first_uand_r = GpXY::new(0.0, radius);
        let last_uand_r = GpXY::new(1.0, radius);
        self.set_radius_at(first_uand_r, ie)?;
        self.set_radius_at(last_uand_r, ie)
    }

    /// `UnSetRadius(const TopoDS_Edge&)` (`cxx:96-119`): drops every point of
    /// `parandrad` between the two ends of `E`.
    pub fn unset_radius_on_edge(&mut self, e: &TopoShape) -> Option<()> {
        self.splitdone = false;
        let ie = self.index_of_edge(e);

        let uf = self.first_parameter_of(ie)?;
        let ul = self.last_parameter_of(ie)?;
        let mut ifirst = 0usize;
        let mut ilast = 0usize;
        for i in 1..=self.fil.parandrad.len() {
            if (self.fil.parandrad[i - 1].x - uf).abs() <= REAL_SMALL {
                ifirst = i;
            }
            if (self.fil.parandrad[i - 1].x - ul).abs() <= REAL_SMALL {
                ilast = i;
            }
        }
        if ifirst != 0 && ilast != 0 {
            // `NCollection_Sequence::Remove(From, To)` removes the range.
            self.fil.parandrad.drain(ifirst - 1..ilast);
        }
        Some(())
    }

    /// `SetRadius(Radius, V)` (`cxx:124-128`).
    pub fn set_radius_on_vertex(&mut self, radius: f64, v: &TopoShape) -> Option<()> {
        let npar = self.absc_of_vertex(v);
        self.set_radius_at(GpXY::new(npar, radius), 0)
    }

    /// `SetRadius(Radius)` (`cxx:132-138`): the constant radius over the whole
    /// spine.
    pub fn set_radius(&mut self, radius: f64) -> Option<()> {
        self.fil.parandrad.clear();
        let first_uand_r = GpXY::new(self.first_parameter(), radius);
        let last_uand_r = GpXY::new(self.last_parameter()?, radius);
        self.set_radius_at(first_uand_r, 0)?;
        self.set_radius_at(last_uand_r, 0)
    }

    /// `SetRadius(const gp_XY& UandR, const int IinC)` (`cxx:143-217`): the
    /// core writer. `IinC == 0` reads `UandR.X()` as an abscissa on the whole
    /// spine, otherwise as a fraction of edge `IinC`'s span.
    ///
    /// When the spine is already split, the law of every elementary spine that
    /// serves `W` is rebuilt (`cxx:191-215`).
    pub fn set_radius_at(&mut self, uand_r: GpXY, i_in_c: i32) -> Option<()> {
        let w = if i_in_c == 0 {
            uand_r.x
        } else {
            let uf = self.first_parameter_of(i_in_c)?;
            let ul = self.last_parameter_of(i_in_c)?;
            uf + uand_r.x * (ul - uf)
        };

        let pr = GpXY::new(w, uand_r.y);
        let mut i = 1usize;
        while i <= self.fil.parandrad.len() {
            let cur = self.fil.parandrad[i - 1];
            if cur.x == w {
                // `cxx:161-174`.
                self.fil.parandrad[i - 1].set_y(uand_r.y);
                if !self.splitdone {
                    return Some(());
                }
                break;
            } else if cur.x > w {
                // `cxx:175-186`.
                self.fil.parandrad.insert(i - 1, pr);
                if !self.splitdone {
                    return Some(());
                }
                break;
            }
            i += 1;
        }
        if i == self.fil.parandrad.len() + 1 {
            // `cxx:187-189`.
            self.fil.parandrad.push(pr);
        }

        // `cxx:191-215`: after a split the law serving `W` has to be replayed.
        if self.splitdone {
            let elspines: Vec<Arc<ElSpine>> = self.el_spines().to_vec();
            if elspines.is_empty() {
                return Some(());
            }
            if elspines[0].is_periodic() {
                // `cxx:198-201`: only the first law is rebuilt.
                let law = self.compute_law(&elspines[0])?;
                self.fil.laws[0] = Arc::new(RwLock::new(law));
            } else {
                // `cxx:202-214`: every elementary spine whose span holds `W`.
                let n = elspines.len().min(self.fil.laws.len());
                for idx in 0..n {
                    let els = &elspines[idx];
                    let uf = els.first_parameter();
                    let ul = els.last_parameter();
                    if uf <= w && w <= ul {
                        let law = self.compute_law(els)?;
                        self.fil.laws[idx] = Arc::new(RwLock::new(law));
                    }
                }
            }
        }
        Some(())
    }

    /// `UnSetRadius(const TopoDS_Vertex&)` (`cxx:221-231`): drops the point at
    /// the vertex's abscissa.
    pub fn unset_radius_on_vertex(&mut self, v: &TopoShape) {
        let npar = self.absc_of_vertex(v);
        for i in 1..=self.fil.parandrad.len() {
            if self.fil.parandrad[i - 1].x == npar {
                self.fil.parandrad.remove(i - 1);
                break;
            }
        }
    }

    /// `SetRadius(const occ::handle<Law_Function>& C, const int IinC)` —
    /// `cxx:235-241`. OCCT builds the composite into the local `prout` and lets
    /// it go out of scope without storing it; only `splitdone` and the cleared
    /// `parandrad` survive (`IinC` is unused there too).
    pub fn set_radius_law(&mut self, c: Arc<dyn LawFunction>, _i_in_c: i32) {
        self.splitdone = false;
        let mut prout = LawComposite::new();
        prout.add_law(c);
        self.fil.parandrad.clear();
    }

    /// `IsConstant()` (`cxx:246-264`): `true` when every radius of `parandrad`
    /// is the same (and there is at least one).
    pub fn is_fil_constant(&self) -> bool {
        if self.fil.parandrad.is_empty() {
            return false;
        }
        let radius = self.fil.parandrad[0].y;
        for p in self.fil.parandrad.iter().skip(1) {
            if (radius - p.y).abs() > CONFUSION {
                return false;
            }
        }
        true
    }

    /// `IsConstant(const int IE)` (`cxx:269-304`): `true` when the radius is
    /// constant over edge `IE`.
    pub fn is_fil_constant_on_edge(&self, ie: i32) -> Option<bool> {
        let uf = self.first_parameter_of(ie)?;
        let ul = self.last_parameter_of(ie)?;

        let len = self.fil.parandrad.len();
        let mut start_rad = 0.0f64;
        let mut i = 1usize;
        while i < len {
            let par = self.fil.parandrad[i - 1].x;
            let rad = self.fil.parandrad[i - 1].y;
            let nextpar = self.fil.parandrad[i].x;
            if (uf - par).abs() <= REAL_SMALL
                || (par < uf && uf < nextpar && nextpar - uf > REAL_SMALL)
            {
                start_rad = rad;
                break;
            }
            i += 1;
        }
        i += 1;
        while i <= len {
            let par = self.fil.parandrad[i - 1].x;
            let rad = self.fil.parandrad[i - 1].y;
            if (rad - start_rad).abs() > CONFUSION {
                return Some(false);
            }
            if (ul - par).abs() <= REAL_SMALL {
                return Some(true);
            }
            if par > ul {
                return Some(true);
            }
            i += 1;
        }
        Some(true)
    }

    /// `Radius(const TopoDS_Edge&)` (`cxx:309-312`).
    pub fn fil_radius_of_edge(&self, e: &TopoShape) -> Option<f64> {
        self.fil_radius_of_edge_index(self.index_of_edge(e))
    }

    /// `Radius(const int IE)` (`cxx:316-354`): the constant radius of edge
    /// `IE`; OCCT raises `Standard_DomainError("Edge is not constant")` when the
    /// scan finds a different radius (`cxx:340-343`).
    pub fn fil_radius_of_edge_index(&self, ie: i32) -> Option<f64> {
        let uf = self.first_parameter_of(ie)?;
        let ul = self.last_parameter_of(ie)?;

        let len = self.fil.parandrad.len();
        let mut start_rad = 0.0f64;
        let mut i = 1usize;
        while i < len {
            let par = self.fil.parandrad[i - 1].x;
            let rad = self.fil.parandrad[i - 1].y;
            let nextpar = self.fil.parandrad[i].x;
            if (uf - par).abs() <= REAL_SMALL
                || (par < uf && uf < nextpar && nextpar - uf > REAL_SMALL)
            {
                start_rad = rad;
                break;
            }
            i += 1;
        }
        i += 1;
        while i <= len {
            let par = self.fil.parandrad[i - 1].x;
            let rad = self.fil.parandrad[i - 1].y;
            if (rad - start_rad).abs() > CONFUSION {
                return None;
            }
            if (ul - par).abs() <= REAL_SMALL {
                return Some(start_rad);
            }
            if par > ul {
                return Some(start_rad);
            }
            i += 1;
        }
        Some(start_rad)
    }

    /// `Radius()` (`cxx:358-364`): the radius of a spine that is constant all
    /// along; OCCT raises `Standard_DomainError("Spine is not constant")`.
    pub fn fil_radius(&self) -> Option<f64> {
        if !self.is_fil_constant() {
            return None;
        }
        Some(self.fil.parandrad.first()?.y)
    }

    /// `AppendElSpine(Els)` (`cxx:368-372`), the `ChFiDS_FilSpine` override:
    /// the base append plus the law of the new elementary spine.
    pub fn append_fil_el_spine(&mut self, els: Arc<ElSpine>) -> Option<()> {
        self.append_el_spine(els.clone());
        self.append_law(&els)
    }

    /// `AppendLaw(Els)` (`cxx:376-381`).
    pub fn append_law(&mut self, els: &ElSpine) -> Option<()> {
        let l = self.compute_law(els)?;
        self.fil.laws.push(Arc::new(RwLock::new(l)));
        Some(())
    }

    /// `Law(Els)` (`cxx:858-870`): the composite attached to elementary spine
    /// `Els`, found by walking `elspines` and `laws` in step (OCCT compares the
    /// handles themselves).
    pub fn fil_law(&self, els: &Arc<ElSpine>) -> Option<Arc<RwLock<LawComposite>>> {
        for (i, s) in self.el_spines().iter().enumerate() {
            if Arc::ptr_eq(s, els) {
                return self.fil.laws.get(i).cloned();
            }
        }
        None
    }

    /// `ChangeLaw(E)` (`cxx:874-890`): the composite that serves edge `E` and
    /// the parameter `W` of its elementary law, so that the caller can replace
    /// that component — OCCT returns `Law_Composite::ChangeElementaryLaw(W)`
    /// (`Law_Composite.cxx:222-226`) by reference.
    ///
    /// OCCT raises `Standard_DomainError` when the limits are not up to date
    /// (`cxx:876-880`) or when `E` is constant (`cxx:882-885`).
    pub fn change_fil_law(&self, e: &TopoShape) -> Option<(Arc<RwLock<LawComposite>>, f64)> {
        if !self.split_done() {
            return None;
        }
        let ie = self.index_of_edge(e);
        if self.is_fil_constant_on_edge(ie)? {
            return None;
        }
        let hsp = self.el_spine_of_edge_index(ie)?;
        let w = 0.5 * (self.first_parameter_of(ie)? + self.last_parameter_of(ie)?);
        let lc = self.fil_law(&hsp)?;
        Some((lc, w))
    }

    /// `MaxRadFromSeqAndLaws()` (`cxx:893-928`): the largest radius found in
    /// `parandrad` and, for every law, at its five equally spaced samples plus
    /// its upper bound (`cxx:909-925`).
    pub fn max_rad_from_seq_and_laws(&self) -> f64 {
        let mut max_rad = 0.0f64;

        for p in &self.fil.parandrad {
            if p.y > max_rad {
                max_rad = p.y;
            }
        }

        for law in &self.fil.laws {
            let guard = law.read().unwrap_or_else(|e| e.into_inner());
            let (fpar, lpar) = guard.bounds();
            let delta = (lpar - fpar) * 0.2;
            for i in 0..=4 {
                let par = fpar + f64::from(i) * delta;
                let rad = guard.value(par);
                if rad > max_rad {
                    max_rad = rad;
                }
            }
            let rad = guard.value(lpar);
            if rad > max_rad {
                max_rad = rad;
            }
        }

        max_rad
    }
}
