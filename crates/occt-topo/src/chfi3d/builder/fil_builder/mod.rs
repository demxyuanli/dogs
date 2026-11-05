//! `ChFi3d_FilBuilder` — the fillet subclass of `ChFi3d_Builder`.
//!
//! Source: `ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_FilBuilder.cxx` (OCCT
//! 8.0.0). This module is the **data layer**: the two file statics `MaxRad` /
//! `SimulParams`, the constructor, `SetFilletShape` / `GetFilletShape`, the two
//! `Add` overloads, the radius / law API, `Simulate`, `NbSurf`, `Sect` and
//! `SimulKPart`.
//!
//! Module layout:
//! ```text
//! statics                MaxRad (`cxx:77-110`), SimulParams (`cxx:112-143`)
//! radius                 SetRadius / IsConstant / Radius / ResetContour / UnSet /
//!                        GetBounds / GetLaw / SetLaw (`cxx:234-398`)
//! simul_k_part           SimulKPart (`cxx:475-554`), the section-shape carrier
//!                        and the `ChFi3dSurfOps` override
//! simul_surf             `SimulSurf` on two faces (`cxx:558-787`)
//! simul_surf_rst         `SimulSurf` on `S1` / on `S2` (`cxx:787-1248`)
//! perform_first_section  `PerformFirstSection` (`cxx:1500-1537`)
//! perform_surf           `PerformSurf` on two faces (`cxx:1538-1717`)
//! perform_surf_rst       `PerformSurf` on `S1` / on `S2` (`cxx:1721-2063`)
//! perform_surf_rst_rst   `PerformSurf` edge / edge (`cxx:2065-2270`)
//! perform_two_corner     `PerformTwoCorner` (`ChFi3d_FilBuilder_C2.cxx:143-1153`)
//! perform_three_corner   `PerformThreeCorner` (`ChFi3d_FilBuilder_C3.cxx:241-1331`)
//! extent_corner          ExtentOne/Two/ThreeCorner (`cxx:2440-2591`)
//! set_regul              SetRegul (`cxx:2595-2634`)
//! ```
//!
//! [`ChFi3dFilBuilder::simulate`] reaches [`FilSurfOps`], which answers all of
//! `SimulKPart`, two of the four `SimulSurf`, all four `PerformSurf` and all
//! eight of the `ChFi3dCornerOps` entry points; the two edge / edge `SimulSurf`
//! overloads keep the `UNPORTED` `Err` of the base class, exactly like OCCT.
//!
//! Divergences that are not heuristics:
//! - OCCT derives `ChFi3d_FilBuilder` from `ChFi3d_Builder`, and the protected
//!   virtuals of the base are member functions of the *same* object. The port
//!   keeps the base as a field (`base`) and puts the only subclass member
//!   (`myShape`) in a separate [`FilSurfOps`] carrier, so the base can take
//!   `&mut dyn ChFi3dSurfOps` while it is itself mutably borrowed — no OCCT
//!   member other than `myShape` is read by the four `SimulSurf` / four
//!   `PerformSurf` bodies (`cxx:558-2273`, verified by grep),
//! - `occ::handle<TopoDS_Edge>` / `const TopoDS_Edge&` become [`Edge`];
//!   `ChFiDS_FilSpine` is the [`Spine`] with `SpineClass::Fil`, so every
//!   `occ::down_cast<ChFiDS_FilSpine>(Value(IC))` becomes
//!   [`ChFi3dBuilder::value`] plus that class check,
//! - `occ::handle<NCollection_HArray1<ChFiDS_CircSection>>` is
//!   [`CircSectionArray`], and a null handle is `None`,
//! - OCCT's `SetRadius` / `UnSet` / ... are `void` and throw from inside
//!   `ChFiDS_FilSpine`; the port answers `Result` on the same arms.

mod extent_corner;
mod perform_first_section;
mod perform_surf;mod perform_surf_rst;
mod perform_surf_rst_rst;
mod perform_three_corner;
mod perform_two_corner;
mod radius;
mod set_regul;
mod simul_k_part;
mod simul_surf;
mod simul_surf_rst;
// `MaxRad` (`ChFi3d_FilBuilder.cxx:77-110`) and `SimulParams` (cxx:112-143):
// `SimulParams` is consumed by the four `SimulSurf` (cxx:558-1248); `MaxRad`
// is read by `ExtentTwoCorner` (`cxx:2499-2502`) through `extent_corner`.
mod statics;

pub use simul_k_part::FilSurfOps;

// `HGuide` as a `Curve` (`BRepBlend_ConstRad` and friends take the guide as an
// `Adaptor3d_Curve`): shared with `FilletSurf_InternalBuilder`, which builds
// the same function objects from its own `PerformSurf` override.
pub(crate) use simul_surf::guide_as_curve;

// The `TopAbs_State` / `FaceState` rename shared by the two subclass
// `PerformFirstSection` bodies (`ChFi3d_FilBuilder.cxx:1500-1537`,
// `ChFi3d_ChBuilder.cxx:1319-1614`).
pub(crate) use perform_first_section::{state_face, state_of};

use std::sync::{Arc, RwLock};

use crate::blend::blend_func::BlendFuncSectionShape;
use crate::chfi3d::builder::ChFi3dBuilder;
use crate::chfi3d::fillet_shape::FilletShape;
use crate::chfi_ds::{CircSectionArray, Spine, SpineClass, SpineHandle, Stripe, StripeHandle};
use crate::shape::{Edge, TopoShape};

/// `ChFi3d_FilBuilder` (`ChFi3d_FilBuilder.hxx:48-378`).
///
/// The base class half lives in `base`; the only member OCCT adds is `myShape`
/// (`hxx:377`), carried by [`FilSurfOps`] — see the module note.
pub struct ChFi3dFilBuilder {
    /// `ChFi3d_Builder` — the base class sub-object.
    base: ChFi3dBuilder,
    /// `myShape` (`hxx:377`).
    my_shape: BlendFuncSectionShape,
}

impl ChFi3dFilBuilder {
    /// `ChFi3d_FilBuilder(S, FShape, Ta)` (`cxx:147-153`):
    /// `ChFi3d_Builder(S, Ta)` then `SetFilletShape(FShape)`.
    pub fn new(shape: &TopoShape, fshape: FilletShape, ta: f64) -> Self {
        let mut this = Self {
            base: ChFi3dBuilder::new(shape, ta),
            my_shape: BlendFuncSectionShape::Rational,
        };
        this.set_fillet_shape(fshape);
        this
    }

    /// The `ChFi3d_Builder` sub-object, for the accessors the base owns
    /// (`IsDone`, `Shape`, `NbElements`, `Value`, `Compute`, ...).
    pub fn base(&self) -> &ChFi3dBuilder {
        &self.base
    }

    /// The mutable `ChFi3d_Builder` sub-object.
    pub fn base_mut(&mut self) -> &mut ChFi3dBuilder {
        &mut self.base
    }

    /// `SetFilletShape(FShape)` (`cxx:157-171`).
    ///
    /// OCCT's `switch` has no `default` arm, so an out-of-range enum value
    /// leaves `myShape` untouched; the port's enum cannot be out of range.
    pub fn set_fillet_shape(&mut self, fshape: FilletShape) {
        self.my_shape = match fshape {
            FilletShape::Rational => BlendFuncSectionShape::Rational,
            FilletShape::QuasiAngular => BlendFuncSectionShape::QuasiAngular,
            FilletShape::Polynomial => BlendFuncSectionShape::Polynomial,
        };
    }

    /// The subclass member `myShape` (`ChFi3d_FilBuilder.hxx:377`) in the form
    /// the blend function objects take: the `PerformSurf` bodies read it
    /// directly (`ChFi3d_FilBuilder.cxx:385-394`), and so does
    /// `FilletSurf_InternalBuilder::PerformSurf` (`FilletSurf_InternalBuilder.cxx:385-394`).
    pub fn fillet_shape(&self) -> BlendFuncSectionShape {
        self.my_shape
    }

    /// `GetFilletShape()` (`cxx:175-193`).
    pub fn get_fillet_shape(&self) -> FilletShape {
        match self.my_shape {
            BlendFuncSectionShape::Rational => FilletShape::Rational,
            BlendFuncSectionShape::QuasiAngular => FilletShape::QuasiAngular,
            BlendFuncSectionShape::Polynomial => FilletShape::Polynomial,
            // `default: break;` with `filshape = ChFi3d_Rational` (`cxx:189-192`).
            BlendFuncSectionShape::Linear => FilletShape::Rational,
        }
    }

    /// `Add(E)` (`cxx:197-218`): starts a new contour on `E` when `E` is not
    /// already in a contour and is an input edge.
    pub fn add(&mut self, e: &Edge) -> Result<(), String> {
        // `if (!Contains(E) && myEFMap.Contains(E))` (`cxx:201`).
        if self.base.contains(e) != 0 || !self.base.ef_map().contains(&e.0) {
            return Ok(());
        }

        // `occ::handle<ChFiDS_Stripe> Stripe = new ChFiDS_Stripe();`
        // `occ::handle<ChFiDS_Spine>& Sp = Stripe->ChangeSpine();`
        // `Sp = new ChFiDS_FilSpine(tolesp);` (`cxx:203-206`).
        let tolesp = self.base.tolesp;
        let spine: SpineHandle = Arc::new(RwLock::new(Spine::with_tolerance_fil(tolesp)));
        let stripe: StripeHandle = Arc::new(RwLock::new(Stripe::new()));
        stripe
            .write()
            .map_err(|_| "Add: poisoned stripe".to_string())?
            .set_spine(spine.clone());

        // `TopoDS_Edge E_wnt = E; E_wnt.Orientation(TopAbs_FORWARD);`
        // `Spine->SetEdges(E_wnt);` (`cxx:208-210`).
        let mut e_wnt = e.0.clone();
        e_wnt.set_orientation(crate::abs::Orientation::Forward);
        spine
            .write()
            .map_err(|_| "Add: poisoned spine".to_string())?
            .set_edges(&e_wnt);

        // `if (PerformElement(Spine, -1, dummy)) { PerformExtremity(Spine);
        // Spine->Load(); myListStripe.Append(Stripe); }` (`cxx:211-217`).
        if self.base.perform_element(&spine, -1.0, None)? {
            self.base.perform_extremity(&spine)?;
            spine
                .write()
                .map_err(|_| "Add: poisoned spine".to_string())?
                .load();
            self.base.my_list_stripe.push(stripe);
        }
        Ok(())
    }

    /// `Add(Radius, E)` (`cxx:222-230`): `Add(E)`, then a constant radius on
    /// the new contour.
    pub fn add_with_radius(&mut self, radius: f64, e: &Edge) -> Result<(), String> {
        self.add(e)?;
        // `int IC = Contains(E); if (IC) { SetRadius(Radius, IC, E); }`.
        let ic = self.base.contains(e);
        if ic != 0 {
            self.set_radius_on_edge(radius, ic, e)?;
        }
        Ok(())
    }

    /// `Simulate(IC)` (`cxx:402-435`): `PerformSetOfSurf(Stripe, true)` on the
    /// `IC`-th stripe.
    ///
    /// The `#ifdef OCCT_DEBUG` chronometer block (`cxx:404-412`, `:423-434`) is
    /// not translated.
    pub fn simulate(&mut self, ic: usize) -> Result<(), String> {
        // `for (itel.Initialize(myListStripe); itel.More(); itel.Next(), i++)`
        // `{ if (i == IC) { PerformSetOfSurf(itel.ChangeValue(), true); break; } }`.
        let some = ic
            .checked_sub(1)
            .and_then(|k| self.base.my_list_stripe.get(k).cloned());
        if let Some(stripe) = some {
            let mut ops = FilSurfOps::new(self.my_shape);
            self.base.perform_set_of_surf(&mut ops, &stripe, true)?;
        }
        Ok(())
    }

    /// `NbSurf(IC)` (`cxx:439-451`): the number of `SurfData` of the `IC`-th
    /// stripe, `0` when there is no such stripe.
    pub fn nb_surf(&self, ic: usize) -> usize {
        let some = ic
            .checked_sub(1)
            .and_then(|k| self.base.my_list_stripe.get(k))
            .and_then(|s| s.read().ok());
        match some {
            Some(guard) => guard.set_of_surf_data().len(),
            None => 0,
        }
    }

    /// `Sect(IC, IS)` (`cxx:455-471`): the simulated section of the `IS`-th
    /// `SurfData` of the `IC`-th stripe; `None` is OCCT's null handle.
    pub fn sect(&self, ic: usize, is: usize) -> Option<Arc<CircSectionArray>> {
        let stripe = ic
            .checked_sub(1)
            .and_then(|k| self.base.my_list_stripe.get(k))?;
        let guard = stripe.read().ok()?;
        let sd = guard.set_of_surf_data().get(is.checked_sub(1)?)?;
        sd.simul().cloned()
    }

    /// `ChFi3d_Builder::Compute()` reached on this subclass
    /// (`BRepFilletAPI_MakeFillet.cxx:374`: `myBuilder.Compute()`).
    ///
    /// Translation boundary: the port's base entry point
    /// [`ChFi3dBuilder::compute`] takes the two groups of `ChFi3d_Builder`
    /// virtuals as `&mut dyn` parameters, because Rust cannot alias one object
    /// into two `&mut` references. The two [`FilSurfOps`] carriers built here
    /// are therefore distinct objects holding the same `myShape`; `Compute`
    /// never writes that field (`set_shape` is only reached from
    /// `set_fillet_shape`), so the two carriers stay identical to the single
    /// `this` OCCT passes implicitly.
    pub fn compute(&mut self) -> Result<(), String> {
        let shape = self.my_shape;
        let mut corner_ops = FilSurfOps::new(shape);
        let mut surf_ops = FilSurfOps::new(shape);
        self.base.compute(&mut corner_ops, &mut surf_ops)
    }

    /// The spine of the `IC`-th stripe as a `ChFiDS_FilSpine`, the
    /// `occ::down_cast<ChFiDS_FilSpine>(Value(IC))` of `cxx:238` and friends.
    ///
    /// OCCT's `down_cast` on a non-fillet spine is undefined behaviour; the
    /// port answers `None` and the callers that would dereference it report
    /// `Err`.
    fn fil_spine(&self, ic: usize) -> Result<SpineHandle, String> {
        let spine = self
            .base
            .value(ic)
            .ok_or_else(|| format!("ChFi3d_FilBuilder: Value({ic}) has no spine"))?;
        let is_fil = spine
            .read()
            .map_err(|_| "ChFi3d_FilBuilder: poisoned spine".to_string())?
            .class()
            == SpineClass::Fil;
        if !is_fil {
            return Err(format!(
                "ChFi3d_FilBuilder: Value({ic}) is not a ChFiDS_FilSpine"
            ));
        }
        Ok(spine)
    }

    /// `IC <= NbElements()` — the guard every radius accessor starts with.
    fn ic_in_range(&self, ic: usize) -> bool {
        ic != 0 && ic <= self.base.nb_elements()
    }
}
