//! `Law_Composite` — a piecewise composition of `Law_Function` components.
//!
//! Source: `Law_Composite.hxx`/`Law_Composite.cxx` (OCCT 8.0.0).
//!
//! Translation notes:
//! - OCCT keeps the `occ::handle<Law_Function>` components in a
//!   `NCollection_List` and caches the component that served the last query in
//!   `curfunc`; the port keeps the components in a `Vec` and the cache as an
//!   index. `Law_Composite::Prepare` is `const`-free in OCCT but is called from
//!   the `const` `Value` / `D1` / `D2`, so the cache lives behind a `Mutex`
//!   (`Law_Composite.cxx:28-87`).
//! - `Law_Composite::Prepare` takes `double& W` and may rewrite it (periodic
//!   normalization); the port takes `&mut f64` for the same reason.
//! - `ElCLib::InPeriod` (`ElCLib.cxx:95-111`) is private here: `Law_Composite`
//!   is its only caller in this crate.

use std::sync::{Arc, Mutex};

use occt_core::kernel::geomabs::Shape;
use occt_core::precision::{epsilon, Precision};

use crate::law::function::LawFunction;

/// `ElCLib::InPeriod` (`ElCLib.cxx:95-111`) with `UFirst` / `ULast` as the
/// period bounds.
fn in_period(u: f64, u_first: f64, u_last: f64) -> f64 {
    // `cxx:97-101`: infinite bounds are left untouched to avoid overflowing
    // the period.
    if Precision::is_infinite(u) || Precision::is_infinite(u_first) || Precision::is_infinite(u_last)
    {
        return u;
    }
    let a_period = u_last - u_first;
    // `cxx:105-108`.
    if a_period < epsilon(u_last) {
        return u;
    }
    // `cxx:110`.
    let shifted = u + a_period * ((u_first - u) / a_period).ceil();
    if shifted > u_first {
        shifted
    } else {
        u_first
    }
}

/// The mutable part of `Law_Composite`: the cached component (`curfunc`) and
/// the domain read from the first / last component on the first `Prepare`.
struct CompositeCache {
    cur: Option<usize>,
    first: f64,
    last: f64,
}

/// `Law_Composite` (`Law_Composite.hxx:33-110`).
pub struct LawComposite {
    laws: Vec<Arc<dyn LawFunction>>,
    periodic: bool,
    t_first: f64,
    t_last: f64,
    ptol: f64,
    cache: Mutex<CompositeCache>,
}

impl Default for LawComposite {
    fn default() -> Self {
        Self::new()
    }
}

impl LawComposite {
    /// `Law_Composite::Law_Composite` (`Law_Composite.cxx:32-40`):
    /// `first = -1e100`, `last = 1e100`, `TFirst = -1e100`, `TLast = 1e100`,
    /// `PTol = 0.`.
    pub fn new() -> Self {
        Self {
            laws: Vec::new(),
            periodic: false,
            t_first: -1.0e100,
            t_last: 1.0e100,
            ptol: 0.0,
            cache: Mutex::new(CompositeCache {
                cur: None,
                first: -1.0e100,
                last: 1.0e100,
            }),
        }
    }

    /// `Law_Composite::Law_Composite(First, Last, Tol)`
    /// (`Law_Composite.cxx:44-54`): `first` / `last` keep their defaults until
    /// `Prepare` reads them from the components.
    pub fn new_with_bounds(t_first: f64, t_last: f64, ptol: f64) -> Self {
        let mut composite = Self::new();
        composite.t_first = t_first;
        composite.t_last = t_last;
        composite.ptol = ptol;
        composite
    }

    /// `Law_Composite::ChangeElementaryLaw(W)` (`Law_Composite.cxx:222-226`):
    /// the component that serves `W`.
    ///
    /// OCCT takes `W` by value, so the periodic normalization `Prepare` applies
    /// is not visible to the caller.
    pub fn elementary_law(&self, w: f64) -> Option<&dyn LawFunction> {
        let mut ww = w;
        let index = self.prepare(&mut ww)?;
        self.laws.get(index).map(|law| law.as_ref())
    }

    /// `Law_Composite::ChangeLaws` (`Law_Composite.cxx:230-233`).
    pub fn laws_mut(&mut self) -> &mut Vec<Arc<dyn LawFunction>> {
        &mut self.laws
    }

    /// `Law_Composite::ChangeElementaryLaw(W)` (`Law_Composite.cxx:222-226`),
    /// by mutable reference: the same component as [`Self::elementary_law`],
    /// handed out so the caller can replace it. OCCT returns
    /// `occ::handle<Law_Function>&` (a reference into `funclist`), which
    /// `ChFiDS_FilSpine::ChangeLaw` (`ChFiDS_FilSpine.cxx:874-890`) and
    /// `ChFi3d_FilBuilder` (`ChFi3d_FilBuilder.cxx:364`, `:396`) write
    /// through.
    pub fn elementary_law_mut(&mut self, w: f64) -> Option<&mut Arc<dyn LawFunction>> {
        let mut ww = w;
        let index = self.prepare(&mut ww)?;
        self.laws.get_mut(index)
    }

    /// The components, in insertion order.
    pub fn laws(&self) -> &[Arc<dyn LawFunction>] {
        &self.laws
    }

    /// Append a component to `funclist` (the port's replacement for reaching
    /// `Law_Composite::ChangeLaws` and pushing a handle).
    pub fn add_law(&mut self, law: Arc<dyn LawFunction>) {
        self.laws.push(law);
    }

    /// `Law_Composite::IsPeriodic` (`Law_Composite.cxx:236-239`).
    pub fn is_periodic(&self) -> bool {
        self.periodic
    }

    /// `Law_Composite::SetPeriodic` (`Law_Composite.cxx:244-247`).
    pub fn set_periodic(&mut self) {
        self.periodic = true;
    }

    fn lock(&self) -> std::sync::MutexGuard<'_, CompositeCache> {
        self.cache.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// `Law_Composite::Prepare` (`Law_Composite.cxx:164-216`): select and cache
    /// the component that serves `w`, normalizing `w` into the domain first
    /// when the composite is periodic. Answers the index of `curfunc`.
    fn prepare(&self, w: &mut f64) -> Option<usize> {
        if self.laws.is_empty() {
            return None;
        }
        // `cxx:167-175`: probe slightly inside the domain so that a query
        // exactly on a shared boundary picks the component that owns it.
        let eps = if *w - self.t_first < self.t_last - *w {
            self.ptol
        } else {
            -self.ptol
        };
        let mut cache = self.lock();
        if cache.cur.is_none() {
            // `cxx:177-184`: the domain is read from the last and first
            // components.
            let last_bounds = self.laws[self.laws.len() - 1].bounds();
            cache.last = last_bounds.1;
            let first_bounds = self.laws[0].bounds();
            cache.first = first_bounds.0;
            cache.cur = Some(0);
        }
        let first = cache.first;
        let last = cache.last;
        let mut w_test = *w + eps;
        if self.periodic {
            // `cxx:186-190`.
            w_test = in_period(w_test, first, last);
            *w = w_test - eps;
        }
        let cur = cache.cur.unwrap_or(0);
        // `cxx:192-196`.
        let (f, l) = self.laws[cur].bounds();
        if f <= w_test && w_test <= l {
            return Some(cur);
        }
        // `cxx:197-201`.
        let chosen = if *w <= first {
            Some(0)
        } else if *w >= last {
            Some(self.laws.len() - 1)
        } else {
            // `cxx:205-215`: the loop leaves `curfunc` on the last component
            // when no component covers `Wtest`.
            let mut cur_index = 0usize;
            for i in 0..self.laws.len() {
                cur_index = i;
                let (f, l) = self.laws[i].bounds();
                if f <= w_test && w_test <= l {
                    break;
                }
            }
            Some(cur_index)
        };
        cache.cur = chosen;
        chosen
    }
}

impl LawFunction for LawComposite {
    /// `Law_Composite::Continuity` (`Law_Composite.cxx:58-61`): not
    /// implemented in OCCT (`Standard_NotImplemented`).
    fn continuity(&self) -> Shape {
        panic!("Law_Composite::Continuity: Standard_NotImplemented")
    }

    /// `Law_Composite::NbIntervals` (`Law_Composite.cxx:66-79`).
    fn nb_intervals(&self, s: Shape) -> Result<i32, String> {
        let mut nbr_interval = 0;
        for law in &self.laws {
            nbr_interval += law.nb_intervals(s)?;
        }
        Ok(nbr_interval)
    }

    /// `Law_Composite::Intervals` (`Law_Composite.cxx:84-103`).
    ///
    /// The first component writes both of its bounds into `T(1)` / `T(2)`;
    /// every component then writes its own interior interval bounds from
    /// `IGlob = 2` on.
    fn intervals(&self, s: Shape) -> Result<Vec<f64>, String> {
        if self.laws.is_empty() {
            return Ok(Vec::new());
        }
        let (first, second) = self.laws[0].bounds();
        let mut result = vec![first, second];
        let mut iglob = 2usize;
        for law in &self.laws {
            let nb_index = law.nb_intervals(s)? as usize + 1;
            let local = law.intervals(s)?;
            // `cxx:99-102`: `Iloc` runs `2..=nb_index`; `IGlob` is the 1-based
            // index into `T`, so the vector slot is `IGlob - 1`.
            for iloc in 2..=nb_index {
                let slot = iglob - 1;
                if result.len() <= slot {
                    result.resize(slot + 1, 0.0);
                }
                result[slot] = local[iloc - 1];
                iglob += 1;
            }
        }
        Ok(result)
    }

    /// `Law_Composite::Value` (`Law_Composite.cxx:107-112`).
    fn value(&self, x: f64) -> f64 {
        let mut w = x;
        match self.prepare(&mut w) {
            Some(index) => self.laws[index].value(w),
            None => f64::NAN,
        }
    }

    /// `Law_Composite::D1` (`Law_Composite.cxx:116-121`).
    fn d1(&self, x: f64) -> (f64, f64) {
        let mut w = x;
        match self.prepare(&mut w) {
            Some(index) => self.laws[index].d1(w),
            None => (f64::NAN, f64::NAN),
        }
    }

    /// `Law_Composite::D2` (`Law_Composite.cxx:125-131`).
    fn d2(&self, x: f64) -> (f64, f64, f64) {
        let mut w = x;
        match self.prepare(&mut w) {
            Some(index) => self.laws[index].d2(w),
            None => (f64::NAN, f64::NAN, f64::NAN),
        }
    }

    /// `Law_Composite::Trim` (`Law_Composite.cxx:137-144`): only the component
    /// list is carried over; `first` / `last` / `periodic` stay at their
    /// constructor values and `Prepare` re-reads the domain on the first query.
    fn trim(&self, p_first: f64, p_last: f64, tol: f64) -> Result<Arc<dyn LawFunction>, String> {
        let mut trimmed = LawComposite::new_with_bounds(p_first, p_last, tol);
        trimmed.laws = self.laws.clone();
        Ok(Arc::new(trimmed))
    }

    /// `Law_Composite::Bounds` (`Law_Composite.cxx:148-152`).
    fn bounds(&self) -> (f64, f64) {
        let cache = self.lock();
        (cache.first, cache.last)
    }
}
